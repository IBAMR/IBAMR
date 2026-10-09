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

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PoissonUtilities.h>
#include <ibtk/SCPoissonPETScLevelSolver.h>

#include <tbox/Database.h>

#include <CoarseFineBoundary.h>
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
#include <string>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Solve a side-centered Poisson problem on the finest level of a hierarchy whose coarse-fine boundaries have convex and
// concave corners, with the side-centered PETSc level solver. The solution is a constant velocity, which the
// discretization reproduces exactly, with the matching constant Dirichlet or Neumann data at the physical boundary and
// exact values in the ghost cells at the coarse-fine boundaries. The solver reproduces the field only if the
// right-hand side is corrected for exactly the couplings that the matrix of the level does not contain.
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
    const bool dirichlet = input_db->getString("bc_type") == "DIRICHLET";

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("sc_poisson_cf_boundary");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f");
    const int u_idx = var_db->registerVariableAndContext(u_var, context, IntVector<NDIM>(1));
    const int f_idx = var_db->registerVariableAndContext(f_var, context, IntVector<NDIM>(0));
    level->allocatePatchData(u_idx);
    level->allocatePatchData(f_idx);

    PoissonSpecifications poisson_spec("poisson_spec");
    poisson_spec.setCZero();
    poisson_spec.setDConstant(-1.0);
    std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> bc_storage;
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM, nullptr);
    for (int axis = 0; axis < NDIM; ++axis)
    {
        bc_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("bc_coef", nullptr));
        for (int face = 0; face < 2 * NDIM; ++face)
        {
            if (dirichlet)
            {
                bc_storage.back()->setBoundaryValue(face, axis + 1.0);
            }
            else
            {
                bc_storage.back()->setBoundarySlope(face, 0.0);
            }
        }
        bc_coefs[axis] = bc_storage.back().get();
    }

    // The exact solution is the constant axis + 1 for each component. The ghost values are exact and the initial guess
    // is zero. The right-hand side vanishes except at Dirichlet sides on the physical boundary, which hold the data.
    const Box<NDIM>& domain_box = level->getPhysicalDomain()[0];
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
        f_data->fillAll(0.0);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            u_data->getArrayData(axis).fillAll(axis + 1.0);
            const Box<NDIM> domain_side_box = SideGeometry<NDIM>::toSideBox(domain_box, axis);
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); b; b++)
            {
                const SideIndex<NDIM> i(b(), axis, SideIndex<NDIM>::Lower);
                (*u_data)(i) = 0.0;
                const bool at_boundary =
                    b()(axis) == domain_side_box.lower(axis) || b()(axis) == domain_side_box.upper(axis);
                if (dirichlet && at_boundary)
                {
                    (*f_data)(i) = axis + 1.0;
                }
            }
        }
    }

    // The correction of the right-hand side at the coarse-fine boundary, for a right-hand side without ghost cells.
    CoarseFineBoundary<NDIM> cf_boundary(*hierarchy, ln, IntVector<NDIM>(1));
    double correction_l1 = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        const SAMRAI::tbox::Array<BoundaryBox<NDIM>>& type_1_cf_bdry =
            cf_boundary.getBoundaries(patch->getPatchNumber(), 1);
        const SAMRAI::tbox::Array<BoundaryBox<NDIM>>& type_2_cf_bdry =
            cf_boundary.getBoundaries(patch->getPatchNumber(), 2);
        SideData<NDIM, double> correction(patch->getBox(), 1, IntVector<NDIM>(0));
        correction.fillAll(0.0);
        PoissonUtilities::adjustRHSAtCoarseFineBoundary(
            correction, *u_data, patch, poisson_spec, type_1_cf_bdry, type_2_cf_bdry, bc_coefs, 0.0, false);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); b; b++)
            {
                correction_l1 += std::abs(correction(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower)));
            }
        }
    }
    correction_l1 = IBTK_MPI::sumReduction(correction_l1);

    SCPoissonPETScLevelSolver solver("solver", input_db->getDatabase("solver_db"), "solver_");
    solver.setPoissonSpecifications(poisson_spec);
    solver.setPhysicalBcCoefs(bc_coefs);
    solver.setHomogeneousBc(false);
    SAMRAIVectorReal<NDIM, double> u("u", hierarchy, ln, ln), f("f", hierarchy, ln, ln);
    u.addComponent(u_var, u_idx);
    f.addComponent(f_var, f_idx);
    solver.initializeSolverState(u, f);
    const bool converged = solver.solveSystem(u, f);
    solver.deallocateSolverState();
    if (!converged)
    {
        TBOX_ERROR("The solver did not converge.\n");
    }

    double error = 0.0;
    bool solution_is_finite = true;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); b; b++)
            {
                const double difference =
                    std::abs((*u_data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower)) - (axis + 1.0));
                solution_is_finite = solution_is_finite && std::isfinite(difference);
                error = std::max(error, difference);
            }
        }
    }
    error = IBTK_MPI::maxReduction(error);
    if (!solution_is_finite)
    {
        TBOX_ERROR("The solution is not finite.\n");
    }
    if (!(error < 1.0e-9))
    {
        TBOX_ERROR("The error of the solution is " << error << ".\n");
    }

    plog << "number_of_patches = " << level->getNumberOfPatches() << '\n';
    plog << "correction_l1 = " << correction_l1 << '\n';
    return 0;
}
