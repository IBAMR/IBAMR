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
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/Database.h>

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <HierarchyCellDataOpsReal.h>
#include <PoissonSpecifications.h>
#include <RefineAlgorithm.h>
#include <SAMRAIVectorReal.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Solve the singular problem -Laplacian(u) = f on a periodic level with depth-1, depth-2 and depth-1 data in turn,
// using one CCPoissonHypreLevelSolver whose null space contains the constant vector. Each use must give, up to a
// constant, the solution of a fresh solver, and must satisfy the discrete equation with the mean of the right-hand side
// removed; the right-hand side of the caller must be unchanged.
namespace
{
// The part of component k of the right-hand side that has zero mean over the level.
double
oscillating_rhs(const double* x, const int k)
{
    const double m = k + 1.0;
    double value = m * std::sin(2.0 * M_PI * (x[0] + 0.1 * m)) * std::cos(2.0 * M_PI * m * x[1]);
    for (int d = 2; d < NDIM; ++d)
    {
        value *= 1.0 + 0.5 * std::sin(2.0 * M_PI * x[d]);
    }
    return value;
}

// The constant part of component k of the right-hand side.
double
constant_rhs(const int k)
{
    return 0.5 + 1.5 * k;
}

void
cell_location(Pointer<Patch<NDIM>> patch, const hier::Index<NDIM>& i, double* x)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    for (int d = 0; d < NDIM; ++d)
    {
        x[d] = pgeom->getXLower()[d] + (i(d) - patch->getBox().lower(d) + 0.5) * pgeom->getDx()[d];
    }
}

void
check_finite(const double value, const std::string& what)
{
    if (!std::isfinite(value))
    {
        TBOX_ERROR(what << " is not finite.\n");
    }
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
    const double tolerance = input_db->getDouble("error_tolerance");
    Pointer<Database> solver_db = input_db->getDatabase("solver_db");

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("cc_poisson_hypre_null_space");

    PoissonSpecifications poisson_spec("poisson_spec");
    poisson_spec.setCZero();
    poisson_spec.setDConstant(-1.0);

    // The control-volume weights of the level, and a scratch variable with a single depth component, because the data
    // operations act on all the depth components of data together.
    HierarchyMathOps math_ops("math_ops", hierarchy, ln, ln);
    const int wgt_idx = math_ops.getCellWeightPatchDescriptorIndex();
    HierarchyCellDataOpsReal<NDIM, double> level_ops(hierarchy, ln, ln);
    const double total_volume = level_ops.sumControlVolumes(wgt_idx, wgt_idx);
    Pointer<CellVariable<NDIM, double>> component_var = new CellVariable<NDIM, double>("component");
    const int component_idx = var_db->registerVariableAndContext(component_var, context, IntVector<NDIM>(0));
    level->allocatePatchData(component_idx);

    // Copy component k of the data with index idx to the scratch variable and return its control-volume-weighted mean.
    auto copy_component = [&](const int idx, const unsigned int k)
    {
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> data = patch->getPatchData(idx);
            Pointer<CellData<NDIM, double>> component_data = patch->getPatchData(component_idx);
            component_data->copyDepth(0, *data, k);
        }
        return level_ops.integral(component_idx, wgt_idx) / total_volume;
    };

    // Subtract the control-volume-weighted mean from component k of the data with index idx.
    auto remove_mean = [&](const int idx, const unsigned int k)
    {
        const double mean = copy_component(idx, k);
        level_ops.addScalar(component_idx, component_idx, -mean);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> data = patch->getPatchData(idx);
            Pointer<CellData<NDIM, double>> component_data = patch->getPatchData(component_idx);
            data->copyDepth(k, *component_data, 0);
        }
    };

    // Solve with the given solver and data of the given depth, and return the solution without its mean in the cells
    // of the patches of this process, the components of a cell next to each other. The control-volume-weighted mean of
    // component k of the right-hand side that the caller supplies is returned in rhs_mean[k].
    auto solve = [&](CCPoissonHypreLevelSolver& solver,
                     const unsigned int depth,
                     const std::string& tag,
                     std::vector<double>& solution,
                     std::vector<double>& rhs_mean)
    {
        Pointer<CellVariable<NDIM, double>> u_var = new CellVariable<NDIM, double>("u_" + tag, depth);
        Pointer<CellVariable<NDIM, double>> f_var = new CellVariable<NDIM, double>("f_" + tag, depth);
        const int u_idx = var_db->registerVariableAndContext(u_var, context, IntVector<NDIM>(1));
        const int f_idx = var_db->registerVariableAndContext(f_var, context, IntVector<NDIM>(0));
        level->allocatePatchData(u_idx);
        level->allocatePatchData(f_idx);

        // The right-hand side, and a copy of it in the order of the cells and components.
        std::vector<double> rhs;
        rhs_mean.assign(depth, 0.0);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<CellData<NDIM, double>> f_data = patch->getPatchData(f_idx);
            u_data->fillAll(0.0);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                double x[NDIM];
                cell_location(patch, ci(), x);
                for (unsigned int k = 0; k < depth; ++k)
                {
                    const double value = oscillating_rhs(x, k) + constant_rhs(k);
                    check_finite(value, "A right-hand side value in use " + tag);
                    (*f_data)(ci(), k) = value;
                    rhs.push_back(value);
                }
            }
        }
        for (unsigned int k = 0; k < depth; ++k)
        {
            rhs_mean[k] = copy_component(f_idx, k);
        }

        solver.setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(depth, nullptr));
        SAMRAIVectorReal<NDIM, double> u("u_" + tag, hierarchy, ln, ln), f("f_" + tag, hierarchy, ln, ln);
        u.addComponent(u_var, u_idx);
        f.addComponent(f_var, f_idx);
        solver.initializeSolverState(u, f);
        const bool converged = solver.solveSystem(u, f);
        solver.deallocateSolverState();
        if (!converged)
        {
            TBOX_ERROR("The solver did not converge in use " << tag << ".\n");
        }

        // The right-hand side of the caller is unchanged.
        size_t n = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> f_data = patch->getPatchData(f_idx);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                for (unsigned int k = 0; k < depth; ++k)
                {
                    if ((*f_data)(ci(), k) != rhs[n++])
                    {
                        TBOX_ERROR("The right-hand side of component " << k << " was changed by the solve in use "
                                                                       << tag << ".\n");
                    }
                }
            }
        }

        // Remove the mean of every component of the solution, which the solver does not specify.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                for (unsigned int k = 0; k < depth; ++k)
                {
                    check_finite((*u_data)(ci(), k), "A solution value in use " + tag);
                }
            }
        }
        for (unsigned int k = 0; k < depth; ++k)
        {
            remove_mean(u_idx, k);
        }
        solution.clear();
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                for (unsigned int k = 0; k < depth; ++k)
                {
                    const double value = (*u_data)(ci(), k);
                    check_finite(value, "A solution value without its mean in use " + tag);
                    solution.push_back(value);
                }
            }
        }

        // The solution satisfies the discrete equation with the mean of each component of the right-hand side removed.
        // A constant added to the solution does not change this.
        RefineAlgorithm<NDIM> ghost_fill;
        ghost_fill.registerRefine(u_idx, u_idx, u_idx, nullptr);
        ghost_fill.createSchedule(level)->fillData(0.0);
        double residual = 0.0;
        n = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                for (unsigned int k = 0; k < depth; ++k)
                {
                    double laplacian = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        CellIndex<NDIM> lower = ci(), upper = ci();
                        --lower(d);
                        ++upper(d);
                        laplacian += ((*u_data)(upper, k) - 2.0 * (*u_data)(ci(), k) + (*u_data)(lower, k)) /
                                     (pgeom->getDx()[d] * pgeom->getDx()[d]);
                    }
                    const double r = (rhs[n++] - rhs_mean[k]) + laplacian;
                    check_finite(r, "A residual value in use " + tag);
                    residual = std::max(residual, std::abs(r));
                }
            }
        }
        residual = IBTK_MPI::maxReduction(residual);
        if (!(residual < tolerance))
        {
            TBOX_ERROR("The residual in use " << tag << " is " << residual << ".\n");
        }

        level->deallocatePatchData(u_idx);
        level->deallocatePatchData(f_idx);
    };

    plog << "number_of_patches = " << level->getNumberOfPatches() << '\n';
    int number_of_cells = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        number_of_cells += level->getPatch(p())->getBox().size();
    }
    plog << "number_of_cells = " << IBTK_MPI::sumReduction(number_of_cells) << '\n';

    // One solver object for all three uses.
    CCPoissonHypreLevelSolver solver("reused_solver", solver_db, "reused_solver_");
    solver.setPoissonSpecifications(poisson_spec);
    solver.setNullSpace(true);
    const unsigned int depths[3] = { 1, 2, 1 };
    for (int use = 0; use < 3; ++use)
    {
        const unsigned int depth = depths[use];
        const std::string number = std::to_string(use + 1);
        std::vector<double> solution, fresh_solution, rhs_mean, fresh_rhs_mean;
        solve(solver, depth, "reused_" + number, solution, rhs_mean);

        CCPoissonHypreLevelSolver fresh_solver("fresh_solver_" + number, solver_db, "fresh_solver_" + number + "_");
        fresh_solver.setPoissonSpecifications(poisson_spec);
        fresh_solver.setNullSpace(true);
        solve(fresh_solver, depth, "fresh_" + number, fresh_solution, fresh_rhs_mean);

        if (solution.size() != fresh_solution.size())
        {
            TBOX_ERROR("The solutions in use " << number << " have different sizes.\n");
        }
        double difference = 0.0;
        std::vector<double> max_norm(depth, 0.0);
        for (size_t n = 0; n < solution.size(); ++n)
        {
            difference = std::max(difference, std::abs(solution[n] - fresh_solution[n]));
            max_norm[n % depth] = std::max(max_norm[n % depth], std::abs(solution[n]));
        }
        difference = IBTK_MPI::maxReduction(difference);
        if (!(difference < tolerance))
        {
            TBOX_ERROR("The solution in use " << number << " differs from the solution of a fresh solver by "
                                              << difference << ".\n");
        }

        plog << "use_" << number << "_depth = " << depth << '\n';
        for (unsigned int k = 0; k < depth; ++k)
        {
            plog << "use_" << number << "_rhs_mean_" << k << " = " << rhs_mean[k] << '\n';
            plog << "use_" << number << "_solution_max_norm_" << k << " = " << IBTK_MPI::maxReduction(max_norm[k])
                 << '\n';
        }
    }
    return 0;
}
