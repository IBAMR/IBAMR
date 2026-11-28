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
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/SCLaplaceOperator.h>
#include <ibtk/VCSCViscousOperator.h>

#include <tbox/Database.h>

#include <CellData.h>
#include <CellVariable.h>
#include <LocationIndexRobinBcCoefs.h>
#include <NodeData.h>
#include <NodeGeometry.h>
#include <NodeVariable.h>
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
#include <utility>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Apply a Laplace operator to a function on a two-level hierarchy. Then change the coefficients and the boundary
// condition objects, initialize the operator again, and apply it again. After the problems, deallocate the operator,
// initialize it for the first problem, and apply it. An operator of the same class that is initialized once for a
// problem must give the same result as the operator that was initialized again, bit for bit. The test stops with an
// error unless the operator reused its HierarchyMathOps object exactly when it should. The last step initializes for
// new vectors of the same structure after the old vectors are freed.
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
};

// Counts the initializations of the operator that reuse its HierarchyMathOps object and those that build one. The
// operator has no public access to the ghost-fill operator, so the test does not count those.
template <class Operator>
class CheckedOperator : public Operator
{
public:
    using Operator::Operator;
    int d_math_ops_reused = 0, d_math_ops_built = 0;

    void initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAIVectorReal<NDIM, double>& out) override
    {
        // Own the previous object during the call, so that a new object cannot have its address.
        Pointer<HierarchyMathOps> previous_math_ops = this->getHierarchyMathOps();
        Operator::initializeOperatorState(in, out);
        (this->getHierarchyMathOps() == previous_math_ops ? d_math_ops_reused : d_math_ops_built) += 1;
    }

    void check(const int math_ops_reused, const int math_ops_built) const
    {
        if (d_math_ops_reused != math_ops_reused || d_math_ops_built != math_ops_built)
        {
            TBOX_ERROR("The operator reused " << d_math_ops_reused << " HierarchyMathOps objects and built "
                                              << d_math_ops_built << ", instead of " << math_ops_reused << " and "
                                              << math_ops_built << ".\n");
        }
    }
};

// The three kinds of operator differ in the vectors, the diffusion coefficient, and the boundary condition objects.
enum class Kind
{
    CC,
    SC,
    VCSC
};

void
run_case(const Kind kind, Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();
    const bool cell = kind == Kind::CC;
    const int number_of_bc_objects = cell ? 1 : NDIM;

    // The problems. The first two have constant diffusion coefficients (a viscosity field, for VCSC) and differ in the
    // coefficients and in the boundary condition objects. The second problem has a Neumann condition at one face of
    // the first object. The third problem of the cell-centered operator has a variable diffusion coefficient, which
    // changes the kind of the fill pattern. The operator then initializes for the first problem again after it has
    // been deallocated, and for the second problem with new vectors after the first vectors are freed.
    struct Step
    {
        int d_problem, d_vectors;
        bool d_deallocate, d_free_old_vectors;
    };
    std::vector<Step> steps = { { 0, 0, false, false }, { 1, 0, false, false } };
    if (cell) steps.push_back({ 2, 0, false, false });
    steps.push_back({ 0, 0, true, false });
    steps.push_back({ 1, 1, false, true });
    const double c[3] = { 1.0, 2.5, 1.7 };
    const double d[3] = { -1.0, -1.5, -1.0 };
    const double mu_scale[2] = { 1.0, 1.6 };
    const int neumann_face = 0; // The lower x face.

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("operator_reinitialize");
    Pointer<SideVariable<NDIM, double>> d_var = new SideVariable<NDIM, double>("d");
    Pointer<NodeVariable<NDIM, double>> mu_var = new NodeVariable<NDIM, double>("mu");
    const int d_idx = var_db->registerVariableAndContext(d_var, context, IntVector<NDIM>(1));
    const int mu_idx = var_db->registerVariableAndContext(mu_var, context, IntVector<NDIM>(1));
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (const int idx : { d_idx, mu_idx })
        {
            level->allocatePatchData(idx);
        }
        // The values depend only on the index of the face or the node, so a face or node that two patches share has
        // one value.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> d_data = patch->getPatchData(d_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(d_data->getArrayData(axis).getBox()); b; b++)
                {
                    int face_hash = axis;
                    for (int dim = 0; dim < NDIM; ++dim)
                    {
                        face_hash += (2 * dim + 3) * b()(dim);
                    }
                    (*d_data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower)) =
                        -1.0 - 0.0625 * (((face_hash % 5) + 5) % 5);
                }
            }
        }
    }
    // Two sets of vectors of the same structure.
    struct Vectors
    {
        Pointer<hier::Variable<NDIM>> d_x_var, d_y_var;
        int d_x_idx = -1, d_y_idx = -1;
        std::unique_ptr<SAMRAIVectorReal<NDIM, double>> d_x, d_y;
    };
    Vectors vectors[2];
    for (int k = 0; k < 2; ++k)
    {
        const std::string suffix = std::to_string(k);
        Vectors& v = vectors[k];
        v.d_x_var = cell ? CellCase::make_variable("x" + suffix) : SideCase::make_variable("x" + suffix);
        v.d_y_var = cell ? CellCase::make_variable("y" + suffix) : SideCase::make_variable("y" + suffix);
        v.d_x_idx = var_db->registerVariableAndContext(v.d_x_var, context, IntVector<NDIM>(1));
        v.d_y_idx = var_db->registerVariableAndContext(v.d_y_var, context, IntVector<NDIM>(1));
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            hierarchy->getPatchLevel(ln)->allocatePatchData(v.d_x_idx);
            hierarchy->getPatchLevel(ln)->allocatePatchData(v.d_y_idx);
        }
        v.d_x = std::make_unique<SAMRAIVectorReal<NDIM, double>>("x" + suffix, hierarchy, 0, finest_ln);
        v.d_y = std::make_unique<SAMRAIVectorReal<NDIM, double>>("y" + suffix, hierarchy, 0, finest_ln);
        v.d_x->addComponent(v.d_x_var, v.d_x_idx);
        v.d_y->addComponent(v.d_y_var, v.d_y_idx);
    }

    // The viscosity is a positive function of the node index, scaled by the problem.
    auto set_mu = [&](const double scale)
    {
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<NodeData<NDIM, double>> mu_data = patch->getPatchData(mu_idx);
                const Box<NDIM> ghost_box = Box<NDIM>::grow(patch->getBox(), IntVector<NDIM>(1));
                for (Box<NDIM>::Iterator b(NodeGeometry<NDIM>::toNodeBox(ghost_box)); b; b++)
                {
                    int node_hash = 0;
                    for (int dim = 0; dim < NDIM; ++dim)
                    {
                        node_hash += (2 * dim + 3) * b()(dim);
                    }
                    (*mu_data)(NodeIndex<NDIM>(b(), NodeIndex<NDIM>::LowerLeft)) =
                        scale * (1.0 + 0.125 * (((node_hash % 5) + 5) % 5));
                }
            }
        }
    };

    std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> bc_storage[2];
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs[2];
    for (int problem = 0; problem < 2; ++problem)
    {
        for (int object = 0; object < number_of_bc_objects; ++object)
        {
            bc_storage[problem].push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("bc_coef", nullptr));
            for (int face = 0; face < 2 * NDIM; ++face)
            {
                bc_storage[problem].back()->setBoundaryValue(face, 0.5 * (problem + 1) * (object + 1));
            }
            if (problem == 1 && object == 0)
            {
                bc_storage[problem].back()->setBoundarySlope(neumann_face, 0.0);
            }
            bc_coefs[problem].push_back(bc_storage[problem].back().get());
        }
    }

    PoissonSpecifications poisson_spec("poisson_spec");
    auto configure = [&](LaplaceOperator& op, const int problem)
    {
        if (kind == Kind::VCSC)
        {
            if (problem == 0)
            {
                poisson_spec.setCZero();
            }
            else
            {
                poisson_spec.setCConstant(c[problem]);
            }
            poisson_spec.setDPatchDataId(mu_idx);
            set_mu(mu_scale[problem]);
        }
        else
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
        }
        op.setPoissonSpecifications(poisson_spec);
        const std::vector<RobinBcCoefStrategy<NDIM>*>& bcs = bc_coefs[problem == 2 ? 0 : problem];
        if (cell)
        {
            op.setPhysicalBcCoef(bcs[0]);
        }
        else
        {
            op.setPhysicalBcCoefs(bcs);
        }
        op.setHomogeneousBc(false);
        op.setSolutionTime(0.0);
        op.setTimeInterval(0.0, 1.0);
    };

    // Apply the operator to the same function, and return the values in the interiors of the patches.
    auto apply = [&](LaplaceOperator& op, Vectors& v)
    {
        int counter = 0;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                auto set = [&](double& value, const bool) { value = 0.5 + 0.25 * (counter++ % 7); };
                auto zero = [&](double& value, const bool) { value = 0.0; };
                if (cell)
                {
                    CellCase::visit(level->getPatch(p()), v.d_x_idx, set);
                    CellCase::visit(level->getPatch(p()), v.d_y_idx, zero);
                }
                else
                {
                    SideCase::visit(level->getPatch(p()), v.d_x_idx, set);
                    SideCase::visit(level->getPatch(p()), v.d_y_idx, zero);
                }
            }
        }
        op.apply(*v.d_x, *v.d_y);
        std::vector<double> result;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                auto collect = [&](double& value, const bool interior)
                {
                    if (interior) result.push_back(value);
                };
                if (cell)
                {
                    CellCase::visit(level->getPatch(p()), v.d_y_idx, collect);
                }
                else
                {
                    SideCase::visit(level->getPatch(p()), v.d_y_idx, collect);
                }
            }
        }
        return result;
    };

    auto make_operator = [&](const std::string& name) -> Pointer<LaplaceOperator>
    {
        switch (kind)
        {
        case Kind::CC:
            return new CheckedOperator<CCLaplaceOperator>(name, nullptr, false);
        case Kind::SC:
            return new CheckedOperator<SCLaplaceOperator>(name, nullptr, false);
        default:
            return new CheckedOperator<VCSCViscousOperator>(name, false);
        }
    };
    auto checked = [&](LaplaceOperator& op, const int math_ops_reused, const int math_ops_built)
    {
        switch (kind)
        {
        case Kind::CC:
            dynamic_cast<CheckedOperator<CCLaplaceOperator>&>(op).check(math_ops_reused, math_ops_built);
            break;
        case Kind::SC:
            dynamic_cast<CheckedOperator<SCLaplaceOperator>&>(op).check(math_ops_reused, math_ops_built);
            break;
        default:
            dynamic_cast<CheckedOperator<VCSCViscousOperator>&>(op).check(math_ops_reused, math_ops_built);
        }
    };

    // The counts of the HierarchyMathOps objects that an operator reused and built.
    auto counts = [&](LaplaceOperator& op)
    {
        switch (kind)
        {
        case Kind::CC:
        {
            auto& checked_op = dynamic_cast<CheckedOperator<CCLaplaceOperator>&>(op);
            return std::make_pair(checked_op.d_math_ops_reused, checked_op.d_math_ops_built);
        }
        case Kind::SC:
        {
            auto& checked_op = dynamic_cast<CheckedOperator<SCLaplaceOperator>&>(op);
            return std::make_pair(checked_op.d_math_ops_reused, checked_op.d_math_ops_built);
        }
        default:
        {
            auto& checked_op = dynamic_cast<CheckedOperator<VCSCViscousOperator>&>(op);
            return std::make_pair(checked_op.d_math_ops_reused, checked_op.d_math_ops_built);
        }
        }
    };

    // The operator that is initialized again at each step, and an operator that is initialized once for the same step.
    // The HierarchyMathOps object is built for the first step and after the deallocation, and reused otherwise.
    Pointer<LaplaceOperator> reinitialized_operator = make_operator("reinitialized_operator");
    std::vector<std::vector<double>> reinitialized, fresh;
    int expected_ops_reused = 0, expected_ops_built = 0;
    for (std::size_t step = 0; step < steps.size(); ++step)
    {
        const Step& current = steps[step];
        Vectors& v = vectors[current.d_vectors];
        if (current.d_deallocate) reinitialized_operator->deallocateOperatorState();
        if (current.d_free_old_vectors)
        {
            // Free the old vectors and remove their indices, as the owners of the vectors do.
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                hierarchy->getPatchLevel(ln)->deallocatePatchData(vectors[0].d_x_idx);
                hierarchy->getPatchLevel(ln)->deallocatePatchData(vectors[0].d_y_idx);
            }
            var_db->removePatchDataIndex(vectors[0].d_x_idx);
            var_db->removePatchDataIndex(vectors[0].d_y_idx);
        }
        configure(*reinitialized_operator, current.d_problem);
        reinitialized_operator->initializeOperatorState(*v.d_x, *v.d_y);
        ((step == 0 || current.d_deallocate) ? expected_ops_built : expected_ops_reused) += 1;
        checked(*reinitialized_operator, expected_ops_reused, expected_ops_built);
        reinitialized.push_back(apply(*reinitialized_operator, v));

        Pointer<LaplaceOperator> fresh_operator = make_operator("fresh_operator");
        configure(*fresh_operator, current.d_problem);
        fresh_operator->initializeOperatorState(*v.d_x, *v.d_y);
        checked(*fresh_operator, 0, 1);
        fresh.push_back(apply(*fresh_operator, v));
        fresh_operator->deallocateOperatorState();
    }
    const auto measured = counts(*reinitialized_operator);
    reinitialized_operator->deallocateOperatorState();

    for (std::size_t step = 0; step < steps.size(); ++step)
    {
        double max_difference = 0.0, max_change = 0.0, result_min = std::numeric_limits<double>::max(),
               result_max = std::numeric_limits<double>::lowest();
        for (std::size_t k = 0; k < fresh[step].size(); ++k)
        {
            if (!std::isfinite(fresh[step][k]) || !std::isfinite(reinitialized[step][k]))
            {
                TBOX_ERROR("A result value is not finite.\n");
            }
            max_difference = std::max(max_difference, std::abs(reinitialized[step][k] - fresh[step][k]));
            if (step > 0)
            {
                max_change = std::max(max_change, std::abs(fresh[step][k] - fresh[step - 1][k]));
            }
            result_min = std::min(result_min, fresh[step][k]);
            result_max = std::max(result_max, fresh[step][k]);
        }
        max_difference = IBTK_MPI::maxReduction(max_difference);
        max_change = IBTK_MPI::maxReduction(max_change);
        result_min = IBTK_MPI::minReduction(result_min);
        result_max = IBTK_MPI::maxReduction(result_max);
        const int number_of_values = IBTK_MPI::sumReduction(static_cast<int>(fresh[step].size()));
        if (reinitialized[step] != fresh[step])
        {
            TBOX_ERROR("The operator that was initialized again for step "
                       << step << " differs from an operator that was initialized once: the results differ by "
                       << max_difference << ".\n");
        }
        if (step > 0 && !(max_change > 1.0e-3))
        {
            TBOX_ERROR("The result did not change from step " << step - 1 << " to step " << step << ": " << max_change
                                                              << ".\n");
        }

        plog << "step = " << step << '\n';
        plog << "problem = " << steps[step].d_problem << '\n';
        plog << "vectors = " << steps[step].d_vectors << '\n';
        plog << "number_of_values = " << number_of_values << '\n';
        plog << "max_difference_reinitialized_fresh = " << max_difference << '\n';
        plog << "max_change_from_previous_step = " << max_change << '\n';
        plog << "result_min = " << result_min << '\n';
        plog << "result_max = " << result_max << '\n';
    }
    plog << "math_ops_reused = " << measured.first << '\n';
    plog << "math_ops_built = " << measured.second << '\n';
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
    if (operator_class == "CCLaplaceOperator")
    {
        run_case(Kind::CC, hierarchy);
    }
    else if (operator_class == "SCLaplaceOperator")
    {
        run_case(Kind::SC, hierarchy);
    }
    else if (operator_class == "VCSCViscousOperator")
    {
        run_case(Kind::VCSC, hierarchy);
    }
    else
    {
        TBOX_ERROR("Unknown operator_class \"" << operator_class << "\"\n"
                                               << "  valid values are CCLaplaceOperator, SCLaplaceOperator, and "
                                                  "VCSCViscousOperator\n");
    }
    return 0;
}
