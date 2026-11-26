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
#include <ibtk/CCPoissonPETScLevelSolver.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/SCPoissonHypreLevelSolver.h>
#include <ibtk/SCPoissonPETScLevelSolver.h>

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
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Solve a Poisson problem on the finest level of a two-level hierarchy with a level solver. Then change the
// coefficients and the boundary values, initialize the solver again, and solve again. A solver of the same class that
// is initialized once with the changed coefficients and boundary values must give the same solution, bit for bit.
namespace
{
struct CellCase
{
    using Data = CellData<NDIM, double>;
    static Pointer<hier::Variable<NDIM>> make_variable(const std::string& name)
    {
        return new CellVariable<NDIM, double>(name);
    }
    template <class F>
    static void visit(Pointer<Patch<NDIM>> patch, const int idx, F f)
    {
        Pointer<Data> data = patch->getPatchData(idx);
        for (Box<NDIM>::Iterator b(data->getGhostBox()); b; b++)
        {
            f((*data)(CellIndex<NDIM>(b())), patch->getBox().contains(b()));
        }
    }
    static void set_bc(PoissonSolver& solver, LocationIndexRobinBcCoefs<NDIM>& bc_coef)
    {
        solver.setPhysicalBcCoef(&bc_coef);
    }
};

struct SideCase
{
    using Data = SideData<NDIM, double>;
    static Pointer<hier::Variable<NDIM>> make_variable(const std::string& name)
    {
        return new SideVariable<NDIM, double>(name);
    }
    template <class F>
    static void visit(Pointer<Patch<NDIM>> patch, const int idx, F f)
    {
        Pointer<Data> data = patch->getPatchData(idx);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
            for (Box<NDIM>::Iterator b(data->getArrayData(axis).getBox()); b; b++)
            {
                f((*data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower)), side_box.contains(b()));
            }
        }
    }
    static void set_bc(PoissonSolver& solver, LocationIndexRobinBcCoefs<NDIM>& bc_coef)
    {
        solver.setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, &bc_coef));
    }
};

// Counts the initializations of a PETSc level solver that compute all the data and the initializations that compute
// only the data that depend on the coefficients. The hypre solvers have no hook that shows which path ran, so the test
// does not check that their second initialization reuses the grid and the stencil.
template <class Solver>
class CheckedSolver : public Solver
{
public:
    using Solver::Solver;
    int d_number_of_full_initializations = 0;
    int d_number_of_coefficient_initializations = 0;

protected:
    void initializeSolverStateSpecialized(const SAMRAIVectorReal<NDIM, double>& x,
                                          const SAMRAIVectorReal<NDIM, double>& b) override
    {
        ++d_number_of_full_initializations;
        Solver::initializeSolverStateSpecialized(x, b);
    }
    void initializeCoefficientsSpecialized(const SAMRAIVectorReal<NDIM, double>& x,
                                           const SAMRAIVectorReal<NDIM, double>& b) override
    {
        ++d_number_of_coefficient_initializations;
        Solver::initializeCoefficientsSpecialized(x, b);
    }
};

template <class Solver>
using TestSolver = std::conditional_t<std::is_base_of_v<PETScLevelSolver, Solver>, CheckedSolver<Solver>, Solver>;

template <class Solver>
void
check_initializations(Solver&, int, int)
{
}

template <class Solver>
void
check_initializations(CheckedSolver<Solver>& solver, const int full, const int coefficient)
{
    if (solver.d_number_of_full_initializations != full ||
        solver.d_number_of_coefficient_initializations != coefficient)
    {
        TBOX_ERROR("The solver ran " << solver.d_number_of_full_initializations << " full initializations and "
                                     << solver.d_number_of_coefficient_initializations
                                     << " coefficient initializations instead of " << full << " and " << coefficient
                                     << ".\n");
    }
}

template <class Case, class Solver>
void
run_case(Pointer<PatchHierarchy<NDIM>> hierarchy, Pointer<Database> input_db)
{
    const int ln = hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("level_solver_reinitialize");
    Pointer<hier::Variable<NDIM>> u_var = Case::make_variable("u");
    Pointer<hier::Variable<NDIM>> f_var = Case::make_variable("f");
    const int u_idx = var_db->registerVariableAndContext(u_var, context, IntVector<NDIM>(1));
    const int f_idx = var_db->registerVariableAndContext(f_var, context, IntVector<NDIM>(0));
    level->allocatePatchData(u_idx);
    level->allocatePatchData(f_idx);
    SAMRAIVectorReal<NDIM, double> u("u", hierarchy, ln, ln), f("f", hierarchy, ln, ln);
    u.addComponent(u_var, u_idx);
    f.addComponent(f_var, f_idx);

    // The coefficients and the boundary values of the first and of the second problem.
    const double c[2] = { 1.0, 2.5 };
    const double d[2] = { -1.0, -1.5 };
    const double boundary_value[2] = { 2.0, 3.0 };
    PoissonSpecifications poisson_spec("poisson_spec");
    LocationIndexRobinBcCoefs<NDIM> bc_coef("bc_coef", nullptr);
    auto configure = [&](PoissonSolver& solver, const int problem)
    {
        poisson_spec.setCConstant(c[problem]);
        poisson_spec.setDConstant(d[problem]);
        for (int face = 0; face < 2 * NDIM; ++face)
        {
            bc_coef.setBoundaryValue(face, boundary_value[problem]);
        }
        solver.setPoissonSpecifications(poisson_spec);
        Case::set_bc(solver, bc_coef);
        solver.setHomogeneousBc(false);
    };

    // Solve with the initialized solver from the same vectors, and return the solution values in the interior.
    auto solve = [&](PoissonSolver& solver, int& iterations)
    {
        int counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Case::visit(level->getPatch(p()),
                        u_idx,
                        [&](double& value, const bool) { value = 1.0 + 0.125 * (counter++ % 13); });
            Case::visit(
                level->getPatch(p()), f_idx, [&](double& value, const bool) { value = 0.5 + 0.25 * (counter++ % 7); });
        }
        solver.solveSystem(u, f);
        iterations = solver.getNumIterations();
        std::vector<double> solution;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Case::visit(level->getPatch(p()),
                        u_idx,
                        [&](double& value, const bool interior)
                        {
                            if (interior)
                            {
                                solution.push_back(value);
                            }
                        });
        }
        return solution;
    };

    // The solver that is initialized twice.
    TestSolver<Solver> reinitialized_solver(
        "reinitialized_solver", input_db->getDatabase("solver_db"), "reinitialized_solver_");
    int iterations_first, iterations_reinitialized, iterations_fresh;
    configure(reinitialized_solver, 0);
    reinitialized_solver.initializeSolverState(u, f);
    check_initializations(reinitialized_solver, 1, 1);
    const std::vector<double> first = solve(reinitialized_solver, iterations_first);
    configure(reinitialized_solver, 1);
    reinitialized_solver.initializeSolverState(u, f);
    check_initializations(reinitialized_solver, 1, 2);
    const std::vector<double> reinitialized = solve(reinitialized_solver, iterations_reinitialized);
    reinitialized_solver.deallocateSolverState();

    // The solver that is initialized once.
    Solver fresh_solver("fresh_solver", input_db->getDatabase("solver_db"), "fresh_solver_");
    configure(fresh_solver, 1);
    fresh_solver.initializeSolverState(u, f);
    const std::vector<double> fresh = solve(fresh_solver, iterations_fresh);
    fresh_solver.deallocateSolverState();

    double max_difference = 0.0, max_change = 0.0, solution_min = std::numeric_limits<double>::max(),
           solution_max = std::numeric_limits<double>::lowest();
    for (std::size_t k = 0; k < fresh.size(); ++k)
    {
        if (!std::isfinite(fresh[k]) || !std::isfinite(reinitialized[k]))
        {
            TBOX_ERROR("A solution value is not finite.\n");
        }
        max_difference = std::max(max_difference, std::abs(reinitialized[k] - fresh[k]));
        max_change = std::max(max_change, std::abs(fresh[k] - first[k]));
        solution_min = std::min(solution_min, fresh[k]);
        solution_max = std::max(solution_max, fresh[k]);
    }
    max_difference = IBTK_MPI::maxReduction(max_difference);
    max_change = IBTK_MPI::maxReduction(max_change);
    solution_min = IBTK_MPI::minReduction(solution_min);
    solution_max = IBTK_MPI::maxReduction(solution_max);
    const int number_of_values = IBTK_MPI::sumReduction(static_cast<int>(fresh.size()));
    if (reinitialized != fresh || iterations_reinitialized != iterations_fresh)
    {
        TBOX_ERROR(
            "The solver that was initialized twice differs from a solver that was initialized once: the "
            "solutions differ by "
            << max_difference << ", and the iteration counts are " << iterations_reinitialized << " and "
            << iterations_fresh << ".\n");
    }
    if (!(max_change > 1.0e-3))
    {
        TBOX_ERROR("The solution did not change with the coefficients and boundary values: " << max_change << ".\n");
    }

    plog << "number_of_values = " << number_of_values << '\n';
    plog << "iterations_first = " << iterations_first << '\n';
    plog << "iterations_reinitialized = " << iterations_reinitialized << '\n';
    plog << "iterations_fresh = " << iterations_fresh << '\n';
    plog << "max_difference_reinitialized_fresh = " << max_difference << '\n';
    plog << "max_change_from_first_problem = " << max_change << '\n';
    plog << "solution_min = " << solution_min << '\n';
    plog << "solution_max = " << solution_max << '\n';
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
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(hierarchy->getFinestLevelNumber());
    plog << "number_of_patches = " << level->getNumberOfPatches() << '\n';

    const std::string solver_class = input_db->getDatabase("test")->getString("solver_class");
    if (solver_class == "CCPoissonPETScLevelSolver")
    {
        run_case<CellCase, CCPoissonPETScLevelSolver>(hierarchy, input_db);
    }
    else if (solver_class == "SCPoissonPETScLevelSolver")
    {
        run_case<SideCase, SCPoissonPETScLevelSolver>(hierarchy, input_db);
    }
    else if (solver_class == "CCPoissonHypreLevelSolver")
    {
        run_case<CellCase, CCPoissonHypreLevelSolver>(hierarchy, input_db);
    }
    else if (solver_class == "SCPoissonHypreLevelSolver")
    {
        run_case<SideCase, SCPoissonHypreLevelSolver>(hierarchy, input_db);
    }
    else
    {
        TBOX_ERROR("Unknown solver_class \"" << solver_class << "\"\n"
                                             << "  valid values are the four Poisson PETSc and hypre level solvers\n");
    }
    return 0;
}
