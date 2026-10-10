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
#include <ibtk/CCPoissonPETScLevelSolver.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
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
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Solve a Helmholtz problem with a constant solution on a single level that is periodic in every direction and has
// several patches, with the cell-centered and the side-centered PETSc level solvers. The solution vector starts with a
// garbage value in its interior and a sentinel value in all ghost cells. This checks that on return the interior values
// hold the solution in every patch's copy of a side that lies on a patch boundary, including the two copies of a side
// across a periodic boundary, and that every ghost value is as it was given. A second solver takes a single Richardson
// iteration, which shows that its default is to ignore the interior values of the solution vector: starting from
// garbage gives the same result as starting from zero.

namespace
{
constexpr double garbage = 7.0;
constexpr double sentinel = -99.0;

struct CellCase
{
    using Data = CellData<NDIM, double>;
    static constexpr const char* name = "cell";
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
    static void setup_solver(Pointer<PoissonSolver>& solver,
                             Pointer<Database> db,
                             const PoissonSpecifications& spec,
                             LocationIndexRobinBcCoefs<NDIM>& bc_coef)
    {
        solver = new CCPoissonPETScLevelSolver("cell_solver", db, "cell_solver_");
        solver->setPoissonSpecifications(spec);
        solver->setPhysicalBcCoef(&bc_coef);
    }
};

struct SideCase
{
    using Data = SideData<NDIM, double>;
    static constexpr const char* name = "side";
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
    static void setup_solver(Pointer<PoissonSolver>& solver,
                             Pointer<Database> db,
                             const PoissonSpecifications& spec,
                             LocationIndexRobinBcCoefs<NDIM>& bc_coef)
    {
        solver = new SCPoissonPETScLevelSolver("side_solver", db, "side_solver_");
        solver->setPoissonSpecifications(spec);
        solver->setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, &bc_coef));
    }
};

template <class Case>
void
run_case(Pointer<PatchHierarchy<NDIM>> hierarchy, Pointer<Database> input_db)
{
    const int ln = hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
    const double exact = input_db->getDouble("exact_solution");
    const double tolerance = input_db->getDouble("error_tolerance");

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(std::string("level_solver_solution_contract_") + Case::name);
    Pointer<hier::Variable<NDIM>> u_var = Case::make_variable(std::string("u_") + Case::name);
    Pointer<hier::Variable<NDIM>> f_var = Case::make_variable(std::string("f_") + Case::name);
    const int u_idx = var_db->registerVariableAndContext(u_var, context, IntVector<NDIM>(1));
    const int f_idx = var_db->registerVariableAndContext(f_var, context, IntVector<NDIM>(0));
    level->allocatePatchData(u_idx);
    level->allocatePatchData(f_idx);

    PoissonSpecifications poisson_spec("poisson_spec");
    poisson_spec.setCConstant(1.0);
    poisson_spec.setDConstant(-1.0);
    LocationIndexRobinBcCoefs<NDIM> bc_coef("bc_coef", nullptr);

    SAMRAIVectorReal<NDIM, double> u("u", hierarchy, ln, ln), f("f", hierarchy, ln, ln);
    u.addComponent(u_var, u_idx);
    f.addComponent(f_var, f_idx);
    f.setToScalar(exact);

    // The initial value of u: garbage in the interior (or zero) and the sentinel in all ghost cells. The garbage is not
    // constant, because a Richardson iteration for this operator maps every constant to the same result.
    auto set_u = [&](const bool use_garbage)
    {
        int counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Case::visit(level->getPatch(p()),
                        u_idx,
                        [&](double& value, const bool interior)
                        { value = interior ? (use_garbage ? garbage + 0.25 * (counter++ % 7) : 0.0) : sentinel; });
        }
    };

    // Solve to convergence from garbage values.
    {
        Pointer<PoissonSolver> solver;
        Case::setup_solver(solver, input_db->getDatabase("solver_db"), poisson_spec, bc_coef);
        solver->setHomogeneousBc(false);
        set_u(true);
        if (!solver->solveSystem(u, f))
        {
            TBOX_ERROR("The " << Case::name << " solver did not converge.\n");
        }
    }
    double u_min = std::numeric_limits<double>::max(), u_max = std::numeric_limits<double>::lowest();
    int number_of_values = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Case::visit(level->getPatch(p()),
                    u_idx,
                    [&](double& value, const bool interior)
                    {
                        if (interior)
                        {
                            if (!std::isfinite(value))
                            {
                                TBOX_ERROR("A " << Case::name << " solution value is not finite.\n");
                            }
                            u_min = std::min(u_min, value);
                            u_max = std::max(u_max, value);
                            ++number_of_values;
                        }
                        else if (value != sentinel)
                        {
                            TBOX_ERROR("A " << Case::name << " ghost value was changed to " << value << ".\n");
                        }
                    });
    }
    u_min = IBTK_MPI::minReduction(u_min);
    u_max = IBTK_MPI::maxReduction(u_max);
    number_of_values = IBTK_MPI::sumReduction(number_of_values);
    if (!(std::abs(u_min - exact) < tolerance && std::abs(u_max - exact) < tolerance))
    {
        TBOX_ERROR("The " << Case::name << " solution ranges from " << u_min << " to " << u_max << " instead of being "
                          << exact << ".\n");
    }

    // Take one iteration from garbage values and from zero.
    std::vector<double> results[2];
    for (int start = 0; start < 2; ++start)
    {
        Pointer<PoissonSolver> solver;
        Case::setup_solver(solver, input_db->getDatabase("one_iteration_solver_db"), poisson_spec, bc_coef);
        solver->setHomogeneousBc(false);
        set_u(start == 0);
        solver->solveSystem(u, f);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Case::visit(level->getPatch(p()),
                        u_idx,
                        [&](double& value, const bool interior)
                        {
                            if (interior)
                            {
                                results[start].push_back(value);
                            }
                        });
        }
    }
    if (results[0] != results[1])
    {
        TBOX_ERROR("The " << Case::name << " solver used the interior values of the initial guess.\n");
    }

    plog << Case::name << "_values = " << number_of_values << '\n';
    plog << Case::name << "_solution_min = " << u_min << '\n';
    plog << Case::name << "_solution_max = " << u_max << '\n';
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
    run_case<CellCase>(hierarchy, input_db);
    run_case<SideCase>(hierarchy, input_db);
    return 0;
}
