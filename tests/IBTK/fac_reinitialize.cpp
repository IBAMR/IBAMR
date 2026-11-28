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
#include <ibtk/CCPoissonBoxRelaxationFACOperator.h>
#include <ibtk/CCPoissonLevelRelaxationFACOperator.h>
#include <ibtk/CCPoissonPointRelaxationFACOperator.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PoissonFACPreconditioner.h>
#include <ibtk/SCPoissonPointRelaxationFACOperator.h>

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
#include <memory>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Apply a FAC preconditioner to a right-hand side on a two-level hierarchy. Then change the coefficients and the
// boundary conditions, initialize the preconditioner again, and apply it again. For the cell-centered box and level
// operators, a third problem with a variable diffusion coefficient takes the full path. A preconditioner of the same
// class that is initialized once with the coefficients and boundary conditions of a problem must give the same result
// as the preconditioner that was initialized again for that problem, bit for bit.
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
    static constexpr int s_number_of_bc_objects = 1;
    static void set_bc(PoissonSolver& solver, const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    {
        solver.setPhysicalBcCoef(bc_coefs[0]);
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
    static constexpr int s_number_of_bc_objects = NDIM;
    static void set_bc(PoissonSolver& solver, const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    {
        solver.setPhysicalBcCoefs(bc_coefs);
    }
};

// Counts the initializations of the operator that compute all the data and the initializations that rebuild only the
// data that depend on the coefficients and the boundary conditions.
template <class Operator>
class CheckedOperator : public Operator
{
public:
    using Operator::Operator;
    int d_number_of_full_initializations = 0;
    int d_number_of_coefficient_initializations = 0;

    void check(const int full, const int coefficient) const
    {
        if (d_number_of_full_initializations != full || d_number_of_coefficient_initializations != coefficient)
        {
            TBOX_ERROR("The operator ran " << d_number_of_full_initializations << " full initializations and "
                                           << d_number_of_coefficient_initializations
                                           << " coefficient initializations instead of " << full << " and "
                                           << coefficient << ".\n");
        }
    }

protected:
    void initializeOperatorStateSpecialized(const SAMRAIVectorReal<NDIM, double>& solution,
                                            const SAMRAIVectorReal<NDIM, double>& rhs,
                                            const int coarsest_reset_ln,
                                            const int finest_reset_ln) override
    {
        ++d_number_of_full_initializations;
        Operator::initializeOperatorStateSpecialized(solution, rhs, coarsest_reset_ln, finest_reset_ln);
    }
    void initializeOperatorCoefficientsSpecialized(const SAMRAIVectorReal<NDIM, double>& solution,
                                                   const SAMRAIVectorReal<NDIM, double>& rhs) override
    {
        ++d_number_of_coefficient_initializations;
        Operator::initializeOperatorCoefficientsSpecialized(solution, rhs);
    }
};

template <class Case, class Operator>
void
run_case(Pointer<PatchHierarchy<NDIM>> hierarchy, Pointer<Database> input_db)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();
    Pointer<Database> test_db = input_db->getDatabase("test");
    const int ghost_width = test_db->getIntegerWithDefault("ghost_width", 1);
    const int neumann_face = test_db->getIntegerWithDefault("neumann_face", -1);
    const int neumann_object = test_db->getIntegerWithDefault("neumann_object", 0);
    const bool variable_d = test_db->getBoolWithDefault("variable_d", false);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("fac_reinitialize");
    Pointer<hier::Variable<NDIM>> u_var = Case::make_variable("u");
    Pointer<hier::Variable<NDIM>> f_var = Case::make_variable("f");
    Pointer<SideVariable<NDIM, double>> d_var = new SideVariable<NDIM, double>("d");
    const int u_idx = var_db->registerVariableAndContext(u_var, context, IntVector<NDIM>(ghost_width));
    const int f_idx = var_db->registerVariableAndContext(f_var, context, IntVector<NDIM>(ghost_width));
    const int d_idx = var_db->registerVariableAndContext(d_var, context, IntVector<NDIM>(1));
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(u_idx);
        hierarchy->getPatchLevel(ln)->allocatePatchData(f_idx);
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_idx);
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        // The value depends only on the index of the face, so a face that two patches share has one value.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<SideData<NDIM, double>> d_data = level->getPatch(p())->getPatchData(d_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(d_data->getArrayData(axis).getBox()); b; b++)
                {
                    int face_hash = axis;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        face_hash += (2 * d + 3) * b()(d);
                    }
                    (*d_data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower)) =
                        -1.0 - 0.0625 * (((face_hash % 5) + 5) % 5);
                }
            }
        }
    }
    SAMRAIVectorReal<NDIM, double> u("u", hierarchy, 0, finest_ln), f("f", hierarchy, 0, finest_ln);
    u.addComponent(u_var, u_idx);
    f.addComponent(f_var, f_idx);

    // The problems. The first two have constant diffusion coefficients and differ in the coefficients and in the
    // boundary condition objects. The second problem has a Neumann condition at one face of one object, if the input
    // gives it. A third problem, if the input sets variable_d, has a variable diffusion coefficient, which changes the
    // kind of the fill pattern of the operator and so changes the layout key. The values of Dirichlet conditions do
    // not enter the error equation of the preconditioner.
    const int number_of_problems = variable_d ? 3 : 2;
    const double c[3] = { 1.0, 2.5, 1.7 };
    const double d[3] = { -1.0, -1.5, -1.0 };
    std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> bc_storage[2];
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs[2];
    for (int problem = 0; problem < 2; ++problem)
    {
        for (int object = 0; object < Case::s_number_of_bc_objects; ++object)
        {
            bc_storage[problem].push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("bc_coef", nullptr));
            for (int face = 0; face < 2 * NDIM; ++face)
            {
                bc_storage[problem].back()->setBoundaryValue(face, 0.0);
            }
            if (problem == 1 && neumann_face >= 0 && object == neumann_object)
            {
                bc_storage[problem].back()->setBoundarySlope(neumann_face, 0.0);
            }
            bc_coefs[problem].push_back(bc_storage[problem].back().get());
        }
    }
    PoissonSpecifications poisson_spec("poisson_spec");
    auto configure = [&](PoissonSolver& solver, const int problem)
    {
        poisson_spec.setCConstant(c[problem]);
        if (problem == 2)
        {
            poisson_spec.setDPatchDataId(d_idx);
        }
        else
        {
            poisson_spec.setDConstant(d[problem]);
        }
        solver.setPoissonSpecifications(poisson_spec);
        Case::set_bc(solver, bc_coefs[problem == 2 ? 0 : problem]);
        solver.setSolutionTime(0.0);
        solver.setTimeInterval(0.0, 1.0);
    };

    // Apply the preconditioner to the same right-hand side, and return the values in the interiors of the patches.
    auto apply = [&](PoissonFACPreconditioner& preconditioner)
    {
        int counter = 0;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Case::visit(level->getPatch(p()), u_idx, [&](double& value, const bool) { value = 0.0; });
                Case::visit(level->getPatch(p()),
                            f_idx,
                            [&](double& value, const bool) { value = 0.5 + 0.25 * (counter++ % 7); });
            }
        }
        preconditioner.solveSystem(u, f);
        std::vector<double> solution;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
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
        }
        return solution;
    };

    // The preconditioner that is initialized again for each problem. The first two problems have the same layout key,
    // so the operator reuses its layout for the second. The third problem takes the full path.
    Pointer<Database> operator_db = input_db->getDatabase("operator_db");
    Pointer<Database> preconditioner_db = input_db->getDatabase("preconditioner_db");
    Pointer<CheckedOperator<Operator>> reinitialized_operator =
        new CheckedOperator<Operator>("reinitialized_operator", operator_db, "reinitialized_operator_");
    PoissonFACPreconditioner reinitialized_preconditioner(
        "reinitialized_preconditioner", reinitialized_operator, preconditioner_db, "reinitialized_preconditioner_");
    std::vector<double> reinitialized[3];
    for (int problem = 0; problem < number_of_problems; ++problem)
    {
        configure(reinitialized_preconditioner, problem);
        reinitialized_preconditioner.initializeSolverState(u, f);
        reinitialized_operator->check(problem == 2 ? 2 : 1, problem == 0 ? 0 : 1);
        reinitialized[problem] = apply(reinitialized_preconditioner);
    }
    reinitialized_preconditioner.deallocateSolverState();

    // The preconditioners that are initialized once.
    std::vector<double> fresh[3];
    for (int problem = 0; problem < number_of_problems; ++problem)
    {
        Pointer<Operator> fresh_operator = new Operator("fresh_operator", operator_db, "fresh_operator_");
        PoissonFACPreconditioner fresh_preconditioner(
            "fresh_preconditioner", fresh_operator, preconditioner_db, "fresh_preconditioner_");
        configure(fresh_preconditioner, problem);
        fresh_preconditioner.initializeSolverState(u, f);
        fresh[problem] = apply(fresh_preconditioner);
        fresh_preconditioner.deallocateSolverState();
    }

    for (int problem = 0; problem < number_of_problems; ++problem)
    {
        double max_difference = 0.0, max_change = 0.0, solution_min = std::numeric_limits<double>::max(),
               solution_max = std::numeric_limits<double>::lowest();
        for (std::size_t k = 0; k < fresh[problem].size(); ++k)
        {
            if (!std::isfinite(fresh[problem][k]) || !std::isfinite(reinitialized[problem][k]))
            {
                TBOX_ERROR("A solution value is not finite.\n");
            }
            max_difference = std::max(max_difference, std::abs(reinitialized[problem][k] - fresh[problem][k]));
            if (problem > 0)
            {
                max_change = std::max(max_change, std::abs(fresh[problem][k] - fresh[problem - 1][k]));
            }
            solution_min = std::min(solution_min, fresh[problem][k]);
            solution_max = std::max(solution_max, fresh[problem][k]);
        }
        max_difference = IBTK_MPI::maxReduction(max_difference);
        max_change = IBTK_MPI::maxReduction(max_change);
        solution_min = IBTK_MPI::minReduction(solution_min);
        solution_max = IBTK_MPI::maxReduction(solution_max);
        const int number_of_values = IBTK_MPI::sumReduction(static_cast<int>(fresh[problem].size()));
        if (reinitialized[problem] != fresh[problem])
        {
            TBOX_ERROR("The preconditioner that was initialized again for problem "
                       << problem << " differs from a preconditioner that was initialized once: the results differ by "
                       << max_difference << ".\n");
        }
        if (problem > 0 && !(max_change > 1.0e-3))
        {
            TBOX_ERROR("The result did not change from problem " << problem - 1 << " to problem " << problem << ": "
                                                                 << max_change << ".\n");
        }

        plog << "problem = " << problem << '\n';
        plog << "number_of_values = " << number_of_values << '\n';
        plog << "max_difference_reinitialized_fresh = " << max_difference << '\n';
        plog << "max_change_from_previous_problem = " << max_change << '\n';
        plog << "solution_min = " << solution_min << '\n';
        plog << "solution_max = " << solution_max << '\n';
    }
    plog << "full_initializations = " << reinitialized_operator->d_number_of_full_initializations << '\n';
    plog << "coefficient_initializations = " << reinitialized_operator->d_number_of_coefficient_initializations << '\n';
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
    plog << "number_of_levels = " << hierarchy->getNumberOfLevels() << '\n';
    plog << "number_of_patches_on_finest_level = " << level->getNumberOfPatches() << '\n';

    const std::string operator_class = input_db->getDatabase("test")->getString("operator_class");
    if (operator_class == "CCPoissonPointRelaxationFACOperator")
    {
        run_case<CellCase, CCPoissonPointRelaxationFACOperator>(hierarchy, input_db);
    }
    else if (operator_class == "CCPoissonBoxRelaxationFACOperator")
    {
        run_case<CellCase, CCPoissonBoxRelaxationFACOperator>(hierarchy, input_db);
    }
    else if (operator_class == "CCPoissonLevelRelaxationFACOperator")
    {
        run_case<CellCase, CCPoissonLevelRelaxationFACOperator>(hierarchy, input_db);
    }
    else if (operator_class == "SCPoissonPointRelaxationFACOperator")
    {
        run_case<SideCase, SCPoissonPointRelaxationFACOperator>(hierarchy, input_db);
    }
    else
    {
        TBOX_ERROR("Unknown operator_class \"" << operator_class << "\"\n"
                                               << "  valid values are the CC point, box, and level relaxation and "
                                                  "the SC point relaxation FAC operators\n");
    }
    return 0;
}
