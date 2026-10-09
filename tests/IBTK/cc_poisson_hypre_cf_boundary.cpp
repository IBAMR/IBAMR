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
#include <ibtk/CCPoissonHypreLevelSolver.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/Database.h>

#include <CellData.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Solve a cell-centered Poisson problem on the finest level of a two-level hierarchy with
// CCPoissonHypreLevelSolver. The fine level has coarse-fine boundaries and does not touch the physical boundary. The
// input selects a multigrid configuration of the solver, which converges only if the matrix has no entries that couple
// to cells outside the level. The exact solution is a constant, which the discretization reproduces. The initial guess
// is zero in the interior and the exact value in the ghost cells, which hold the boundary data at the coarse-fine
// boundary. The input sets initial_guess_nonzero = FALSE.
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
    const double exact = input_db->getDouble("exact_solution");
    const double tolerance = input_db->getDouble("error_tolerance");

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("cc_poisson_hypre_cf_boundary");
    Pointer<CellVariable<NDIM, double>> u_var = new CellVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> f_var = new CellVariable<NDIM, double>("f");
    const int u_idx = var_db->registerVariableAndContext(u_var, context, IntVector<NDIM>(1));
    const int f_idx = var_db->registerVariableAndContext(f_var, context, IntVector<NDIM>(0));
    level->allocatePatchData(u_idx);
    level->allocatePatchData(f_idx);

    PoissonSpecifications poisson_spec("poisson_spec");
    poisson_spec.setCZero();
    poisson_spec.setDConstant(-1.0);
    LocationIndexRobinBcCoefs<NDIM> bc_coef("bc_coef", nullptr);
    for (int face = 0; face < 2 * NDIM; ++face)
    {
        bc_coef.setBoundaryValue(face, exact);
    }

    CCPoissonHypreLevelSolver solver("solver", input_db->getDatabase("solver_db"), "solver_");
    solver.setPoissonSpecifications(poisson_spec);
    solver.setPhysicalBcCoef(&bc_coef);
    solver.setHomogeneousBc(false);
    SAMRAIVectorReal<NDIM, double> u("u", hierarchy, ln, ln), f("f", hierarchy, ln, ln);
    u.addComponent(u_var, u_idx);
    f.addComponent(f_var, f_idx);
    solver.initializeSolverState(u, f);

    // Start from zero in the interior and the exact value in the ghost cells.
    f.setToScalar(0.0);
    u.setToScalar(exact, /*interior_only*/ false);
    u.setToScalar(0.0, /*interior_only*/ true);
    const bool converged = solver.solveSystem(u, f);
    solver.deallocateSolverState();
    if (!converged)
    {
        TBOX_ERROR("The solver did not converge.\n");
    }

    double u_min = std::numeric_limits<double>::max(), u_max = std::numeric_limits<double>::lowest();
    int number_of_cells = 0;
    bool solution_is_finite = true;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
        {
            const double value = (*u_data)(ci());
            solution_is_finite = solution_is_finite && std::isfinite(value);
            u_min = std::min(u_min, value);
            u_max = std::max(u_max, value);
            ++number_of_cells;
        }
    }
    u_min = IBTK_MPI::minReduction(u_min);
    u_max = IBTK_MPI::maxReduction(u_max);
    number_of_cells = IBTK_MPI::sumReduction(number_of_cells);
    if (!solution_is_finite)
    {
        TBOX_ERROR("The solution is not finite.\n");
    }
    if (!(std::abs(u_min - exact) < tolerance && std::abs(u_max - exact) < tolerance))
    {
        TBOX_ERROR("The solution ranges from " << u_min << " to " << u_max << " instead of being " << exact << ".\n");
    }

    plog << "number_of_patches = " << level->getNumberOfPatches() << '\n';
    plog << "number_of_cells = " << number_of_cells << '\n';
    plog << "solution_min = " << u_min << '\n';
    plog << "solution_max = " << u_max << '\n';
    return 0;
}
