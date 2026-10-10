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
#include <ibtk/SideSynchCopyFillPattern.h>

#include <tbox/Database.h>

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <RefineAlgorithm.h>
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

// Solve a Stokes problem twice with the Stokes PETSc level solver on a level of several patches that lie side by side
// along a traction boundary, with the real velocity and pressure boundary condition objects of the integrator. At a
// traction boundary the velocity objects read the normal velocity of the solution vector at two adjacent positions
// along the boundary, which at the end of a patch includes a side that belongs to the next patch. The level is periodic
// in the direction along the boundary and has a coarse-fine boundary opposite to it. The solution vector holds the same
// interior values at the start of both solves and different values in its ghost values that lie in other patches of
// the level. This checks that the two solutions agree, that every copy of a side that patches share holds the solution
// on return, that the ghost values at the coarse-fine and physical boundaries and all pressure ghost values are as they
// were given, and that the ghost values that lie in other patches of the level hold the interior values that the
// neighbors had at the start of the solve, not the solution.
namespace
{
constexpr double same_level_ghost[2] = { 7.0, -5.0 };
constexpr double coarse_fine_ghost = 2.0;
constexpr double physical_ghost = 3.5;
constexpr double pressure_ghost = 1.25;

// The velocity that the caller supplies in the interior, a function of the location of a side.
double
caller_velocity(const double* x, const int axis)
{
    return 0.5 * (axis + 1.0) * std::cos(2.0 * M_PI * x[0]) + 0.3 * (axis + 2.0) * x[1];
}

void
side_location(Pointer<Patch<NDIM>> patch, const hier::Index<NDIM>& i, const int axis, double* x)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    for (int d = 0; d < NDIM; ++d)
    {
        x[d] = pgeom->getXLower()[d] + (i(d) - patch->getBox().lower(d)) * pgeom->getDx()[d] +
               (axis == d ? 0.0 : 0.5 * pgeom->getDx()[d]);
    }
}

// 0: interior side of the patch, 1: ghost side that is a side of the level, 2: ghost side at a coarse-fine boundary,
// 3: ghost side outside a physical boundary.
int
side_category(Pointer<PatchLevel<NDIM>> level, Pointer<Patch<NDIM>> patch, const hier::Index<NDIM>& g, const int axis)
{
    if (SideGeometry<NDIM>::toSideBox(patch->getBox(), axis).contains(g))
    {
        return 0;
    }
    const IntVector<NDIM> shift = level->getGridGeometry()->getPeriodicShift(level->getRatio());
    const Box<NDIM> domain = level->getPhysicalDomain()[0];
    for (int d = 0; d < NDIM; ++d)
    {
        if (shift(d) == 0 && (g(d) < domain.lower(d) || g(d) > domain.upper(d) + (axis == d ? 1 : 0)))
        {
            return 3;
        }
    }
    const BoxArray<NDIM>& boxes = level->getBoxes();
    int number_of_images = 1;
    for (int d = 0; d < NDIM; ++d)
    {
        number_of_images *= 3;
    }
    for (int m = 0; m < number_of_images; ++m)
    {
        hier::Index<NDIM> h = g;
        int rest = m;
        bool valid = true;
        for (int d = 0; d < NDIM; ++d)
        {
            const int k = (rest % 3) - 1;
            rest /= 3;
            if (k != 0 && shift(d) == 0)
            {
                valid = false;
            }
            h(d) += k * shift(d);
        }
        if (!valid)
        {
            continue;
        }
        for (int b = 0; b < boxes.getNumberOfBoxes(); ++b)
        {
            if (SideGeometry<NDIM>::toSideBox(boxes[b], axis).contains(h))
            {
                return 1;
            }
        }
    }
    return 2;
}
} // namespace

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
    Pointer<VariableContext> context = var_db->getContext("stokes_level_solver_solution_contract");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f");
    Pointer<CellVariable<NDIM, double>> h_var = new CellVariable<NDIM, double>("h");
    Pointer<SideVariable<NDIM, double>> check_var = new SideVariable<NDIM, double>("check");
    const int u_idx = var_db->registerVariableAndContext(u_var, context, IntVector<NDIM>(1));
    const int p_idx = var_db->registerVariableAndContext(p_var, context, IntVector<NDIM>(1));
    const int f_idx = var_db->registerVariableAndContext(f_var, context, IntVector<NDIM>(1));
    const int h_idx = var_db->registerVariableAndContext(h_var, context, IntVector<NDIM>(1));
    const int check_idx = var_db->registerVariableAndContext(check_var, context, IntVector<NDIM>(1));
    for (const int idx : { u_idx, p_idx, f_idx, h_idx, check_idx })
    {
        level->allocatePatchData(idx);
    }

    // The boundary condition objects of the integrator: Neumann conditions on the velocity at the lower boundary are
    // traction conditions.
    Pointer<INSStaggeredHierarchyIntegrator> integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator", input_db->getDatabase("INSStaggeredHierarchyIntegrator"), false);
    std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> physical_storage;
    std::vector<RobinBcCoefStrategy<NDIM>*> physical_bcs(NDIM, nullptr);
    for (int axis = 0; axis < NDIM; ++axis)
    {
        physical_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("physical_bc", nullptr));
        for (int face = 0; face < 2 * NDIM; ++face)
        {
            physical_storage.back()->setBoundarySlope(face, 0.5);
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
    StaggeredStokesPETScLevelSolver solver("contract_solver", input_db->getDatabase("solver_db"), "contract_");
    solver.setVelocityPoissonSpecifications(coefficients);
    solver.setComponentsHaveNullSpace(false, false);
    solver.setPhysicalBcCoefs(velocity_bcs, &pressure_bc);
    solver.setPhysicalBoundaryHelper(helper);
    solver.setHomogeneousBc(false);

    SAMRAIVectorReal<NDIM, double> x("x", hierarchy, ln, ln), b("b", hierarchy, ln, ln);
    x.addComponent(u_var, u_idx);
    x.addComponent(p_var, p_idx);
    b.addComponent(f_var, f_idx);
    b.addComponent(h_var, h_idx);
    b.setToScalar(0.0, /*interior_only*/ false);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
            {
                double location[NDIM];
                side_location(patch, it(), axis, location);
                (*f_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) =
                    std::sin(2.0 * M_PI * location[0]) + 0.25 * (axis + 1.0) * location[1];
            }
        }
    }

    // The solution vector at the start of a solve.
    auto set_x = [&](const int state)
    {
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(u_data->getArrayData(axis).getBox()); it; it++)
                {
                    double location[NDIM];
                    side_location(patch, it(), axis, location);
                    double value = 0.0;
                    switch (side_category(level, patch, it(), axis))
                    {
                    case 0:
                        value = caller_velocity(location, axis);
                        break;
                    case 1:
                        value = same_level_ghost[state];
                        break;
                    case 2:
                        value = coarse_fine_ghost;
                        break;
                    default:
                        value = physical_ghost;
                    }
                    (*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) = value;
                }
            }
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
            for (Box<NDIM>::Iterator it(p_data->getGhostBox()); it; it++)
            {
                (*p_data)(CellIndex<NDIM>(it())) = patch->getBox().contains(it()) ? 0.0 : pressure_ghost;
            }
        }
    };

    std::vector<double> solution[2];
    solver.initializeSolverState(x, b);
    for (int state = 0; state < 2; ++state)
    {
        set_x(state);
        if (!solver.solveSystem(x, b))
        {
            TBOX_ERROR("The solver did not converge.\n");
        }

        // Copy the solution to the check data and synchronize the sides that patches share. This changes nothing if
        // every copy of a shared side holds the solution already.
        double synchronization_change = 0.0, distance_from_caller = 0.0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<SideData<NDIM, double>> check_data = patch->getPatchData(check_idx);
            check_data->copy(*u_data);
        }
        RefineAlgorithm<NDIM> synchronization;
        synchronization.registerRefine(check_idx, check_idx, check_idx, nullptr, new SideSynchCopyFillPattern());
        synchronization.createSchedule(level)->fillData(0.0);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<SideData<NDIM, double>> check_data = patch->getPatchData(check_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(u_data->getArrayData(axis).getBox()); it; it++)
                {
                    const SideIndex<NDIM> i(it(), axis, SideIndex<NDIM>::Lower);
                    double location[NDIM];
                    side_location(patch, it(), axis, location);
                    const double value = (*u_data)(i);
                    switch (side_category(level, patch, it(), axis))
                    {
                    case 0:
                        if (!std::isfinite(value))
                        {
                            TBOX_ERROR("A velocity solution value is not finite.\n");
                        }
                        synchronization_change = std::max(synchronization_change, std::abs(value - (*check_data)(i)));
                        distance_from_caller =
                            std::max(distance_from_caller, std::abs(value - caller_velocity(location, axis)));
                        solution[state].push_back(value);
                        break;
                    case 1:
                        if (std::abs(value - caller_velocity(location, axis)) > 1.0e-13)
                        {
                            TBOX_ERROR("A velocity ghost value that lies in another patch of the level is "
                                       << value << " instead of the interior value " << caller_velocity(location, axis)
                                       << " that the neighbor had at the start of the solve.\n");
                        }
                        break;
                    case 2:
                        if (value != coarse_fine_ghost)
                        {
                            TBOX_ERROR("A coarse-fine ghost value was changed to " << value << ".\n");
                        }
                        break;
                    default:
                        if (value != physical_ghost)
                        {
                            TBOX_ERROR("A physical boundary ghost value was changed to " << value << ".\n");
                        }
                    }
                }
            }
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
            for (Box<NDIM>::Iterator it(p_data->getGhostBox()); it; it++)
            {
                if (!patch->getBox().contains(it()) && (*p_data)(CellIndex<NDIM>(it())) != pressure_ghost)
                {
                    TBOX_ERROR("A pressure ghost value was changed to " << (*p_data)(CellIndex<NDIM>(it())) << ".\n");
                }
            }
        }
        synchronization_change = IBTK_MPI::maxReduction(synchronization_change);
        distance_from_caller = IBTK_MPI::maxReduction(distance_from_caller);
        if (synchronization_change != 0.0)
        {
            TBOX_ERROR("The copies of a side that patches share differ by " << synchronization_change << ".\n");
        }
        if (!(distance_from_caller > 1.0e-3))
        {
            TBOX_ERROR("The solution is too close to the interior values that the caller supplied.\n");
        }
    }
    solver.deallocateSolverState();

    double l2 = 0.0, difference = 0.0;
    for (size_t k = 0; k < solution[0].size(); ++k)
    {
        l2 += solution[0][k] * solution[0][k];
        difference = std::max(difference, std::abs(solution[0][k] - solution[1][k]));
    }
    l2 = std::sqrt(IBTK_MPI::sumReduction(l2));
    difference = IBTK_MPI::maxReduction(difference);
    if (!(difference < 1.0e-8))
    {
        TBOX_ERROR("The two solutions differ by " << difference << ".\n");
    }
    const int number_of_sides = IBTK_MPI::sumReduction(static_cast<int>(solution[0].size()));

    plog << "number_of_patches = " << level->getNumberOfPatches() << '\n';
    plog << "number_of_sides = " << number_of_sides << '\n';
    plog << "velocity_l2_norm = " << l2 << '\n';
    return 0;
}
