// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/INSStaggeredVelocityBcCoef.h>
#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/Database.h>

#include <CellData.h>
#include <CellVariable.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

// Solve a Stokes problem on the finest level of a hierarchy with the Stokes PETSc level solver and the velocity and
// boundary condition objects of the Navier-Stokes integrator, which read the solution vector as the target velocity
// whenever they are evaluated. The lower boundary of the first coordinate direction has a traction condition with zero
// data and the other boundaries have Dirichlet conditions with a constant value for each velocity component. The finest
// level touches the traction boundary, and a coarse-fine boundary of the level ends at it. A velocity that has the
// constant value of the Dirichlet data in each component, with zero pressure, solves the problem. The solution vector
// holds that velocity in all of its values and a pressure of one, and the solution must have zero pressure and the
// velocity.
int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app->getInputDatabase();
    const auto hierarchy_data = setup_hierarchy<NDIM>(app);
    Pointer<PatchHierarchy<NDIM>> hierarchy = std::get<0>(hierarchy_data);
    const int ln = hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("stokes_petsc_level_solver_bc_objects");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f");
    Pointer<CellVariable<NDIM, double>> h_var = new CellVariable<NDIM, double>("h");
    const int u_idx = var_db->registerVariableAndContext(u_var, context, IntVector<NDIM>(1));
    const int p_idx = var_db->registerVariableAndContext(p_var, context, IntVector<NDIM>(1));
    const int f_idx = var_db->registerVariableAndContext(f_var, context, IntVector<NDIM>(1));
    const int h_idx = var_db->registerVariableAndContext(h_var, context, IntVector<NDIM>(1));
    for (const int idx : { u_idx, p_idx, f_idx, h_idx })
    {
        level->allocatePatchData(idx);
    }

    // The boundary condition objects of the integrator. A Neumann condition on the velocity is a traction condition.
    Pointer<INSStaggeredHierarchyIntegrator> integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator", input_db->getDatabase("INSStaggeredHierarchyIntegrator"), false);
    std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> physical_storage;
    std::vector<RobinBcCoefStrategy<NDIM>*> physical_bcs(NDIM, nullptr);
    for (int axis = 0; axis < NDIM; ++axis)
    {
        physical_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("physical_bc", nullptr));
        for (int face = 0; face < 2 * NDIM; ++face)
        {
            if (face == 0)
            {
                physical_storage.back()->setBoundarySlope(face, 0.0);
            }
            else
            {
                physical_storage.back()->setBoundaryValue(face, axis + 1.0);
            }
        }
        physical_bcs[axis] = physical_storage.back().get();
    }
    std::vector<std::unique_ptr<INSStaggeredVelocityBcCoef>> velocity_storage;
    std::vector<RobinBcCoefStrategy<NDIM>*> velocity_bcs(NDIM, nullptr);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        velocity_storage.push_back(
            std::make_unique<INSStaggeredVelocityBcCoef>(axis, integrator, physical_bcs, TRACTION));
        velocity_bcs[axis] = velocity_storage.back().get();
    }
    INSStaggeredPressureBcCoef pressure_bc(integrator, physical_bcs, TRACTION);
    Pointer<StaggeredStokesPhysicalBoundaryHelper> helper = new StaggeredStokesPhysicalBoundaryHelper();
    helper->cacheBcCoefData(physical_bcs, 0.0, hierarchy);

    PoissonSpecifications coefficients("coefficients");
    coefficients.setCConstant(1.0);
    coefficients.setDConstant(-1.0);
    StaggeredStokesPETScLevelSolver solver("bc_objects_solver", input_db->getDatabase("solver_db"), "bc_objects_");
    solver.setVelocityPoissonSpecifications(coefficients);
    solver.setComponentsHaveNullSpace(false, false);
    solver.setPhysicalBcCoefs(velocity_bcs, &pressure_bc);
    solver.setPhysicalBoundaryHelper(helper);
    solver.setHomogeneousBc(false);

    // The velocity component axis has the constant value axis + 1, which solves -Laplace(u) + u + grad(p) = u with
    // p = 0, so the force equals the velocity. The solution vector holds the velocity in all of its values, which the
    // boundary condition objects read at the traction boundary, and a pressure of one.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
        Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
        Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
        p_data->fillAll(1.0);
        h_data->fillAll(0.0);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(u_data->getArrayData(axis).getBox()); it; it++)
            {
                const SideIndex<NDIM> i(it(), axis, SideIndex<NDIM>::Lower);
                (*u_data)(i) = axis + 1.0;
                (*f_data)(i) = axis + 1.0;
            }
        }
    }

    SAMRAIVectorReal<NDIM, double> x("x", hierarchy, ln, ln), b("b", hierarchy, ln, ln);
    x.addComponent(u_var, u_idx);
    x.addComponent(p_var, p_idx);
    b.addComponent(f_var, f_idx);
    b.addComponent(h_var, h_idx);
    solver.initializeSolverState(x, b);
    if (!solver.solveSystem(x, b))
    {
        TBOX_ERROR("The solver did not converge.\n");
    }
    solver.deallocateSolverState();

    double velocity_error = 0.0, pressure_error = 0.0;
    int number_of_sides = 0;
    bool solution_is_finite = true;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
            {
                const double velocity_difference =
                    std::abs((*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) - (axis + 1.0));
                solution_is_finite = solution_is_finite && std::isfinite(velocity_difference);
                velocity_error = std::max(velocity_error, velocity_difference);
                ++number_of_sides;
            }
        }
        for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
        {
            const double pressure_value = (*p_data)(CellIndex<NDIM>(it()));
            solution_is_finite = solution_is_finite && std::isfinite(pressure_value);
            pressure_error = std::max(pressure_error, std::abs(pressure_value));
        }
    }
    velocity_error = IBTK_MPI::maxReduction(velocity_error);
    pressure_error = IBTK_MPI::maxReduction(pressure_error);
    number_of_sides = IBTK_MPI::sumReduction(number_of_sides);
    if (!solution_is_finite)
    {
        TBOX_ERROR("The solution is not finite.\n");
    }
    if (!(velocity_error < 1.0e-9))
    {
        TBOX_ERROR("The error in the velocity is " << velocity_error << ".\n");
    }
    if (!(pressure_error < 1.0e-9))
    {
        TBOX_ERROR("The error in the pressure is " << pressure_error << ".\n");
    }

    plog << "number_of_patches = " << level->getNumberOfPatches() << '\n';
    plog << "number_of_sides = " << number_of_sides << '\n';
    return 0;
}
