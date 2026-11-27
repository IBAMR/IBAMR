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

#include <ibamr/StaggeredStokesBlockFactorizationPreconditioner.h>
#include <ibamr/StaggeredStokesBlockPreconditioner.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/StaggeredStokesProjectionPreconditioner.h>
#include <ibamr/VCStaggeredStokesOperator.h>
#include <ibamr/VCStaggeredStokesProjectionPreconditioner.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LinearSolver.h>
#include <ibtk/PoissonSolver.h>

#include <tbox/Database.h>
#include <tbox/Logger.h>

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
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

// Apply a Stokes operator or preconditioner on two levels, change the coefficients and boundary conditions, initialize
// it again, and apply it again; also with new vectors and after a deallocation. Each result must equal, bit for bit,
// that of an object initialized once, and the HierarchyMathOps and data operations must be reused exactly when
// expected.
namespace
{
enum class Kind
{
    STOKES,
    VC_STOKES,
    BLOCK,
    BLOCK_NULL_SPACE,
    BLOCK_FACTORIZATION,
    PROJECTION,
    VC_PROJECTION
};

// The counts of the initializations that reused the HierarchyMathOps object and the data operations, and of those that
// built them. Only the objects with a public or protected accessor are counted.
struct Counts
{
    int d_math_ops_reused = 0, d_math_ops_built = 0, d_data_ops_reused = 0, d_data_ops_built = 0;

    void check(const Counts& expected) const
    {
        if (d_math_ops_reused != expected.d_math_ops_reused || d_math_ops_built != expected.d_math_ops_built ||
            d_data_ops_reused != expected.d_data_ops_reused || d_data_ops_built != expected.d_data_ops_built)
        {
            TBOX_ERROR("The initializations reused "
                       << d_math_ops_reused << " HierarchyMathOps objects and built " << d_math_ops_built
                       << ", and reused " << d_data_ops_reused << " data operations and built " << d_data_ops_built
                       << ", instead of " << expected.d_math_ops_reused << ", " << expected.d_math_ops_built << ", "
                       << expected.d_data_ops_reused << ", and " << expected.d_data_ops_built << ".\n");
        }
    }
};

// Counts the initializations of an operator by whether they reuse its HierarchyMathOps object.
template <class Operator>
class CheckedOperator : public Operator
{
public:
    using Operator::Operator;
    Counts d_counts;

    void initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAIVectorReal<NDIM, double>& out) override
    {
        // Own the previous object during the call, so that a new object cannot have its address.
        Pointer<HierarchyMathOps> previous_math_ops = this->getHierarchyMathOps();
        Operator::initializeOperatorState(in, out);
        (this->getHierarchyMathOps() == previous_math_ops ? d_counts.d_math_ops_reused : d_counts.d_math_ops_built) +=
            1;
    }
};

// Counts the initializations of a preconditioner by whether they reuse its HierarchyMathOps object and its data
// operations. The velocity and pressure data operations count as one pair.
template <class Preconditioner>
class CheckedPreconditioner : public Preconditioner
{
public:
    using Preconditioner::Preconditioner;
    Counts d_counts;

    void initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                               const SAMRAIVectorReal<NDIM, double>& b) override
    {
        Pointer<HierarchyMathOps> previous_math_ops = this->d_hier_math_ops;
        Pointer<HierarchyDataOpsReal<NDIM, double>> previous_velocity_ops = this->d_velocity_data_ops;
        Pointer<HierarchyDataOpsReal<NDIM, double>> previous_pressure_ops = this->d_pressure_data_ops;
        Preconditioner::initializeSolverState(x, b);
        (this->d_hier_math_ops == previous_math_ops ? d_counts.d_math_ops_reused : d_counts.d_math_ops_built) += 1;
        const bool data_ops_reused =
            this->d_velocity_data_ops == previous_velocity_ops && this->d_pressure_data_ops == previous_pressure_ops;
        (data_ops_reused ? d_counts.d_data_ops_reused : d_counts.d_data_ops_built) += 1;
    }
};

// A subdomain solver that copies the right-hand side. It makes the results of the preconditioners depend on the
// preconditioner, and not on a solver.
class CopySolver : public PoissonSolver, public LinearSolver
{
public:
    bool solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b) override
    {
        x.copyVector(Pointer<SAMRAIVectorReal<NDIM, double>>(&b, false));
        return true;
    }
};

// A block preconditioner that copies the right-hand side, removes the null space vectors, and records the volume of the
// physical domain and the integrals of the result.
class BlockPreconditioner : public StaggeredStokesBlockPreconditioner
{
public:
    BlockPreconditioner(const bool needs_pressure_solver)
        : StaggeredStokesBlockPreconditioner(/*needs_velocity_solver*/ false, needs_pressure_solver)
    {
        GeneralSolver::init("BlockPreconditioner", /*homogeneous_bc*/ true);
    }

    double d_volume = 0.0, d_U_integral = 0.0, d_P_integral_before = 0.0, d_P_integral_after = 0.0;

    bool solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b) override
    {
        x.copyVector(Pointer<SAMRAIVectorReal<NDIM, double>>(&b, false));
        Pointer<SAMRAIVectorReal<NDIM, double>> U_vec =
            new SAMRAIVectorReal<NDIM, double>("U", d_hierarchy, d_coarsest_ln, d_finest_ln);
        U_vec->addComponent(
            x.getComponentVariable(0), x.getComponentDescriptorIndex(0), d_velocity_wgt_idx, d_velocity_data_ops);
        Pointer<SAMRAIVectorReal<NDIM, double>> P_vec =
            new SAMRAIVectorReal<NDIM, double>("P", d_hierarchy, d_coarsest_ln, d_finest_ln);
        P_vec->addComponent(
            x.getComponentVariable(1), x.getComponentDescriptorIndex(1), d_pressure_wgt_idx, d_pressure_data_ops);
        d_P_integral_before = d_pressure_data_ops->integral(x.getComponentDescriptorIndex(1), d_pressure_wgt_idx);
        correctNullSpace(U_vec, P_vec);
        d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();
        d_U_integral = d_velocity_data_ops->integral(x.getComponentDescriptorIndex(0), d_velocity_wgt_idx);
        d_P_integral_after = d_pressure_data_ops->integral(x.getComponentDescriptorIndex(1), d_pressure_wgt_idx);
        return true;
    }

    void initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                               const SAMRAIVectorReal<NDIM, double>& b) override
    {
        StaggeredStokesBlockPreconditioner::initializeSolverState(x, b);
        d_is_initialized = true;
    }

    void deallocateSolverState() override
    {
        StaggeredStokesBlockPreconditioner::deallocateSolverState();
        d_is_initialized = false;
    }
};

// The vectors of a problem: variables, data, and vectors with control volumes.
struct Vectors
{
    Pointer<SideVariable<NDIM, double>> d_u_var, d_f_var;
    Pointer<CellVariable<NDIM, double>> d_p_var, d_h_var;
    int d_u_idx = -1, d_p_idx = -1, d_f_idx = -1, d_h_idx = -1;
    std::unique_ptr<SAMRAIVectorReal<NDIM, double>> d_x, d_b;
};

std::unique_ptr<Vectors>
make_vectors(Pointer<PatchHierarchy<NDIM>> hierarchy, HierarchyMathOps& math_ops, const std::string& suffix)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("stokes_operator_reinitialize");
    auto vectors = std::make_unique<Vectors>();
    vectors->d_u_var = new SideVariable<NDIM, double>("u" + suffix);
    vectors->d_p_var = new CellVariable<NDIM, double>("p" + suffix);
    vectors->d_f_var = new SideVariable<NDIM, double>("f" + suffix);
    vectors->d_h_var = new CellVariable<NDIM, double>("h" + suffix);
    vectors->d_u_idx = var_db->registerVariableAndContext(vectors->d_u_var, context, IntVector<NDIM>(1));
    vectors->d_p_idx = var_db->registerVariableAndContext(vectors->d_p_var, context, IntVector<NDIM>(1));
    vectors->d_f_idx = var_db->registerVariableAndContext(vectors->d_f_var, context, IntVector<NDIM>(1));
    vectors->d_h_idx = var_db->registerVariableAndContext(vectors->d_h_var, context, IntVector<NDIM>(1));
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        for (const int idx : { vectors->d_u_idx, vectors->d_p_idx, vectors->d_f_idx, vectors->d_h_idx })
        {
            hierarchy->getPatchLevel(ln)->allocatePatchData(idx);
        }
    }
    const int side_wgt_idx = math_ops.getSideWeightPatchDescriptorIndex();
    const int cell_wgt_idx = math_ops.getCellWeightPatchDescriptorIndex();
    vectors->d_x = std::make_unique<SAMRAIVectorReal<NDIM, double>>("x" + suffix, hierarchy, 0, finest_ln);
    vectors->d_x->addComponent(vectors->d_u_var, vectors->d_u_idx, side_wgt_idx);
    vectors->d_x->addComponent(vectors->d_p_var, vectors->d_p_idx, cell_wgt_idx);
    vectors->d_b = std::make_unique<SAMRAIVectorReal<NDIM, double>>("b" + suffix, hierarchy, 0, finest_ln);
    vectors->d_b->addComponent(vectors->d_f_var, vectors->d_f_idx, side_wgt_idx);
    vectors->d_b->addComponent(vectors->d_h_var, vectors->d_h_idx, cell_wgt_idx);
    return vectors;
}

void
run_case(const Kind kind, Pointer<PatchHierarchy<NDIM>> hierarchy, Pointer<Database> pc_db_source)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();
    const bool is_operator = kind == Kind::STOKES || kind == Kind::VC_STOKES;
    const bool is_block = kind == Kind::BLOCK || kind == Kind::BLOCK_NULL_SPACE;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("stokes_operator_reinitialize");
    HierarchyMathOps math_ops("test_math_ops", hierarchy, 0, finest_ln);
    std::unique_ptr<Vectors> vectors[2] = { make_vectors(hierarchy, math_ops, "0"),
                                            make_vectors(hierarchy, math_ops, "1") };

    // The viscosity and the cell-centered coefficient for the variable-coefficient classes.
    Pointer<NodeVariable<NDIM, double>> mu_var = new NodeVariable<NDIM, double>("mu");
    Pointer<CellVariable<NDIM, double>> d_cc_var = new CellVariable<NDIM, double>("d_cc");
    const int mu_idx = var_db->registerVariableAndContext(mu_var, context, IntVector<NDIM>(1));
    const int d_cc_idx = var_db->registerVariableAndContext(d_cc_var, context, IntVector<NDIM>(1));
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(mu_idx);
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_cc_idx);
    }
    auto set_coefficient_data = [&](const double scale)
    {
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM> ghost_box = Box<NDIM>::grow(patch->getBox(), IntVector<NDIM>(1));
                Pointer<NodeData<NDIM, double>> mu_data = patch->getPatchData(mu_idx);
                for (Box<NDIM>::Iterator b(NodeGeometry<NDIM>::toNodeBox(ghost_box)); b; b++)
                {
                    int hash = 0;
                    for (int dim = 0; dim < NDIM; ++dim) hash += (2 * dim + 3) * b()(dim);
                    (*mu_data)(NodeIndex<NDIM>(b(), NodeIndex<NDIM>::LowerLeft)) =
                        scale * (1.0 + 0.125 * (((hash % 5) + 5) % 5));
                }
                Pointer<CellData<NDIM, double>> d_data = patch->getPatchData(d_cc_idx);
                for (Box<NDIM>::Iterator b(d_data->getGhostBox()); b; b++)
                {
                    int hash = 0;
                    for (int dim = 0; dim < NDIM; ++dim) hash += (2 * dim + 3) * b()(dim);
                    (*d_data)(CellIndex<NDIM>(b())) = scale * (1.0 + 0.0625 * (((hash % 7) + 7) % 7));
                }
            }
        }
    };

    // The two problems. The first has slip walls and Neumann pressure conditions. The second has other boundary
    // condition objects, with Dirichlet data and a Dirichlet pressure condition at one face.
    const double c[2] = { 1.0, 2.5 };
    const double d[2] = { -1.0, -1.5 };
    const double mu_scale[2] = { 1.0, 1.6 };
    std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> bc_storage;
    std::vector<RobinBcCoefStrategy<NDIM>*> velocity_bcs[2];
    RobinBcCoefStrategy<NDIM>* pressure_bcs[2] = { nullptr, nullptr };
    Pointer<StaggeredStokesPhysicalBoundaryHelper> helpers[2];
    for (int problem = 0; problem < 2; ++problem)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            bc_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("velocity_bc", nullptr));
            for (int face = 0; face < 2 * NDIM; ++face)
            {
                if (face / 2 == axis)
                {
                    bc_storage.back()->setBoundaryValue(face, problem == 0 ? 0.0 : 0.25 * (axis + 1));
                }
                else if (problem == 0)
                {
                    bc_storage.back()->setBoundarySlope(face, 0.0);
                }
                else
                {
                    bc_storage.back()->setBoundaryValue(face, 0.125);
                }
            }
            velocity_bcs[problem].push_back(bc_storage.back().get());
        }
        bc_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("pressure_bc", nullptr));
        for (int face = 0; face < 2 * NDIM; ++face)
        {
            if (problem == 1 && face == 0)
            {
                bc_storage.back()->setBoundaryValue(face, 0.3);
            }
            else
            {
                bc_storage.back()->setBoundarySlope(face, 0.0);
            }
        }
        pressure_bcs[problem] = bc_storage.back().get();
        helpers[problem] = new StaggeredStokesPhysicalBoundaryHelper();
        helpers[problem]->cacheBcCoefData(velocity_bcs[problem], 0.0, hierarchy);
    }

    PoissonSpecifications u_spec("u_spec"), p_spec("p_spec");
    auto configure = [&](StaggeredStokesSolver* solver, StaggeredStokesOperator* op, const int problem)
    {
        u_spec.setCConstant(c[problem]);
        p_spec.setCZero();
        p_spec.setDConstant(-1.0 / (1.0 + problem));
        if (kind == Kind::VC_STOKES)
        {
            u_spec.setDPatchDataId(mu_idx);
            set_coefficient_data(mu_scale[problem]);
        }
        else
        {
            u_spec.setDConstant(d[problem]);
        }
        if (kind == Kind::VC_PROJECTION) set_coefficient_data(mu_scale[problem]);
        if (op)
        {
            op->setVelocityPoissonSpecifications(u_spec);
            op->setPhysicalBcCoefs(velocity_bcs[problem], pressure_bcs[problem]);
            op->setPhysicalBoundaryHelper(helpers[problem]);
            op->setHomogeneousBc(false);
            op->setSolutionTime(1.0);
            op->setTimeInterval(0.0, 1.0);
        }
        else
        {
            solver->setVelocityPoissonSpecifications(u_spec);
            solver->setPhysicalBcCoefs(velocity_bcs[problem], pressure_bcs[problem]);
            solver->setPhysicalBoundaryHelper(helpers[problem]);
            solver->setSolutionTime(1.0);
            solver->setTimeInterval(0.0, 1.0);
        }
    };

    // Fill the data of one set of vectors with the same values every time.
    auto fill = [&](Vectors& v)
    {
        int counter = 0;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(v.d_u_idx);
                Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(v.d_p_idx);
                Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(v.d_f_idx);
                Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(v.d_h_idx);
                u_data->fillAll(0.0);
                p_data->fillAll(0.0);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(f_data->getArrayData(axis).getBox()); it; it++)
                    {
                        (*f_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) = 0.5 + 0.25 * (counter++ % 7);
                    }
                }
                for (Box<NDIM>::Iterator it(h_data->getGhostBox()); it; it++)
                {
                    (*h_data)(CellIndex<NDIM>(it())) = 0.125 * (counter++ % 5);
                }
                // Operators act on the left-hand side, so give it the same values.
                if (is_operator)
                {
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (Box<NDIM>::Iterator it(u_data->getArrayData(axis).getBox()); it; it++)
                        {
                            (*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) =
                                0.25 + 0.125 * (counter++ % 5);
                        }
                    }
                    for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
                    {
                        (*p_data)(CellIndex<NDIM>(it())) = 0.5 + 0.0625 * (counter++ % 3);
                    }
                }
            }
        }
    };
    // The kind of each collected value: 0 velocity, 1 pressure, 2 neither.
    std::vector<char> value_kinds;
    auto collect = [&](const std::vector<int>& u_indices, const std::vector<int>& p_indices)
    {
        std::vector<double> values;
        value_kinds.clear();
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                for (const int u_idx : u_indices)
                {
                    Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
                        {
                            values.push_back((*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)));
                            value_kinds.push_back(0);
                        }
                    }
                }
                for (const int p_idx : p_indices)
                {
                    Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
                    for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
                    {
                        values.push_back((*p_data)(CellIndex<NDIM>(it())));
                        value_kinds.push_back(1);
                    }
                }
            }
        }
        return values;
    };

    // The steps. Each step has a problem, a set of vectors, whether the object is deallocated first, and whether the
    // first vectors are freed first. Preconditioners initialize for new variables in the third step, and every object
    // initializes for the first problem after a deallocation. An operator ends with new vectors after the first
    // vectors are freed.
    struct Step
    {
        int d_problem, d_vectors;
        bool d_deallocate, d_free_old_vectors;
    };
    std::vector<Step> steps = { { 0, 0, false, false }, { 1, 0, false, false } };
    if (!is_operator) steps.push_back({ 1, 1, false, false });
    steps.push_back({ 0, 0, true, false });
    if (is_operator) steps.push_back({ 1, 1, false, true });

    // The expected counts. The HierarchyMathOps object is built for the first step and after the deallocation. The
    // data operations are built then, and for the new variables.
    Counts expected;
    for (std::size_t step = 0; step < steps.size(); ++step)
    {
        const bool built = step == 0 || steps[step].d_deallocate;
        (built ? expected.d_math_ops_built : expected.d_math_ops_reused) += 1;
        const bool new_variables = !is_operator && step > 0 && steps[step].d_vectors != steps[step - 1].d_vectors;
        (built || new_variables ? expected.d_data_ops_built : expected.d_data_ops_reused) += 1;
    }
    if (is_operator)
    {
        expected.d_data_ops_reused = 0;
        expected.d_data_ops_built = 0;
    }

    // The object under test. A subject initializes for a problem and a set of vectors, evaluates, and deallocates.
    struct Subject
    {
        std::function<void(int, Vectors&)> d_initialize;
        std::function<void()> d_deallocate;
        std::function<std::vector<double>(Vectors&)> d_evaluate;
        std::function<Counts()> d_counts;
    };
    auto make_subject = [&](const std::string& name) -> Subject
    {
        Subject subject;
        if (is_operator)
        {
            auto make = [&](auto op)
            {
                subject.d_initialize = [&, op](const int problem, Vectors& v)
                {
                    configure(nullptr, op.getPointer(), problem);
                    op->initializeOperatorState(*v.d_x, *v.d_b);
                };
                subject.d_deallocate = [&, op]() { op->deallocateOperatorState(); };
                subject.d_evaluate = [&, op](Vectors& v)
                {
                    fill(v);
                    op->apply(*v.d_x, *v.d_b);
                    return collect({ v.d_f_idx }, { v.d_h_idx });
                };
                subject.d_counts = [&, op]() { return op->d_counts; };
            };
            if (kind == Kind::STOKES)
            {
                make(Pointer<CheckedOperator<StaggeredStokesOperator>>(
                    new CheckedOperator<StaggeredStokesOperator>(name, false)));
            }
            else
            {
                make(Pointer<CheckedOperator<VCStaggeredStokesOperator>>(
                    new CheckedOperator<VCStaggeredStokesOperator>(name, false)));
            }
            return subject;
        }
        auto make = [&](auto pc, const bool null_space)
        {
            Pointer<CopySolver> velocity_solver = new CopySolver(), pressure_solver = new CopySolver();
            if (null_space) pressure_solver->setNullSpace(true);
            if (!is_block) pc->setVelocitySubdomainSolver(velocity_solver);
            if (kind != Kind::BLOCK) pc->setPressureSubdomainSolver(pressure_solver);
            subject.d_initialize = [&, pc](const int problem, Vectors& v)
            {
                configure(pc.getPointer(), nullptr, problem);
                pc->setPressurePoissonSpecifications(p_spec);
                pc->initializeSolverState(*v.d_x, *v.d_b);
            };
            subject.d_deallocate = [&, pc]() { pc->deallocateSolverState(); };
            subject.d_counts = [&, pc]() { return pc->d_counts; };
        };
        switch (kind)
        {
        case Kind::BLOCK:
        case Kind::BLOCK_NULL_SPACE:
        {
            Pointer<CheckedPreconditioner<BlockPreconditioner>> pc =
                new CheckedPreconditioner<BlockPreconditioner>(kind == Kind::BLOCK_NULL_SPACE);
            make(pc, kind == Kind::BLOCK_NULL_SPACE);
            subject.d_evaluate = [&, pc](Vectors& v)
            {
                fill(v);
                pc->solveSystem(*v.d_x, *v.d_b);
                std::vector<double> values = collect({ v.d_u_idx }, { v.d_p_idx });
                for (const double value :
                     { pc->d_volume, pc->d_U_integral, pc->d_P_integral_before, pc->d_P_integral_after })
                {
                    values.push_back(value);
                    value_kinds.push_back(2);
                }
                return values;
            };
            break;
        }
        case Kind::BLOCK_FACTORIZATION:
        {
            Pointer<CheckedPreconditioner<StaggeredStokesBlockFactorizationPreconditioner>> pc =
                new CheckedPreconditioner<StaggeredStokesBlockFactorizationPreconditioner>(name, pc_db_source, "pc_");
            make(pc, false);
            subject.d_evaluate = [&, pc](Vectors& v)
            {
                fill(v);
                pc->solveSystem(*v.d_x, *v.d_b);
                return collect({ v.d_u_idx }, { v.d_p_idx });
            };
            break;
        }
        case Kind::PROJECTION:
        {
            Pointer<CheckedPreconditioner<StaggeredStokesProjectionPreconditioner>> pc =
                new CheckedPreconditioner<StaggeredStokesProjectionPreconditioner>(name, pc_db_source, "pc_");
            make(pc, false);
            subject.d_evaluate = [&, pc](Vectors& v)
            {
                fill(v);
                pc->solveSystem(*v.d_x, *v.d_b);
                return collect({ v.d_u_idx }, { v.d_p_idx });
            };
            break;
        }
        default:
        {
            Pointer<CheckedPreconditioner<VCStaggeredStokesProjectionPreconditioner>> pc =
                new CheckedPreconditioner<VCStaggeredStokesProjectionPreconditioner>(name, pc_db_source, "pc_");
            pc->setVelocityCellCenteredDCoefficient(d_cc_idx);
            make(pc, false);
            subject.d_evaluate = [&, pc](Vectors& v)
            {
                fill(v);
                pc->solveSystem(*v.d_x, *v.d_b);
                return collect({ v.d_u_idx }, { v.d_p_idx });
            };
        }
        }
        return subject;
    };

    // The object that is initialized again at each step, and an object that is initialized once for the same step.
    Subject reinitialized_subject = make_subject("reinitialized_subject");
    std::vector<std::vector<double>> reinitialized, fresh;
    std::vector<std::vector<char>> fresh_value_kinds;
    Counts expected_fresh;
    expected_fresh.d_math_ops_built = 1;
    expected_fresh.d_data_ops_built = is_operator ? 0 : 1;
    for (std::size_t k = 0; k < steps.size(); ++k)
    {
        const Step& step = steps[k];
        Vectors& v = *vectors[step.d_vectors];
        if (step.d_deallocate) reinitialized_subject.d_deallocate();
        if (step.d_free_old_vectors)
        {
            // Free the first vectors and remove their indices, as the owners of the vectors do.
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                for (const int idx :
                     { vectors[0]->d_u_idx, vectors[0]->d_p_idx, vectors[0]->d_f_idx, vectors[0]->d_h_idx })
                {
                    hierarchy->getPatchLevel(ln)->deallocatePatchData(idx);
                }
            }
            for (const int idx : { vectors[0]->d_u_idx, vectors[0]->d_p_idx, vectors[0]->d_f_idx, vectors[0]->d_h_idx })
            {
                var_db->removePatchDataIndex(idx);
            }
        }
        reinitialized_subject.d_initialize(step.d_problem, v);
        reinitialized.push_back(reinitialized_subject.d_evaluate(v));

        // A new name for each step: a second StaggeredStokesBlockFactorizationPreconditioner with the same name stops
        // with an error.
        Subject fresh_subject = make_subject("fresh_subject_" + std::to_string(k));
        fresh_subject.d_initialize(step.d_problem, v);
        fresh.push_back(fresh_subject.d_evaluate(v));
        fresh_value_kinds.push_back(value_kinds);
        fresh_subject.d_counts().check(expected_fresh);
        fresh_subject.d_deallocate();
    }
    const Counts measured = reinitialized_subject.d_counts();
    measured.check(expected);
    reinitialized_subject.d_deallocate();

    for (std::size_t step = 0; step < steps.size(); ++step)
    {
        double max_difference = 0.0, max_change = 0.0;
        double range_min[2] = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
        double range_max[2] = { std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest() };
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
            const int kind_of_value = fresh_value_kinds[step][k];
            if (kind_of_value < 2)
            {
                range_min[kind_of_value] = std::min(range_min[kind_of_value], fresh[step][k]);
                range_max[kind_of_value] = std::max(range_max[kind_of_value], fresh[step][k]);
            }
        }
        max_difference = IBTK_MPI::maxReduction(max_difference);
        max_change = IBTK_MPI::maxReduction(max_change);
        for (int kind_of_value = 0; kind_of_value < 2; ++kind_of_value)
        {
            range_min[kind_of_value] = IBTK_MPI::minReduction(range_min[kind_of_value]);
            range_max[kind_of_value] = IBTK_MPI::maxReduction(range_max[kind_of_value]);
        }
        const int number_of_values = IBTK_MPI::sumReduction(static_cast<int>(fresh[step].size()));
        if (reinitialized[step] != fresh[step])
        {
            TBOX_ERROR("The object that was initialized again for step "
                       << step << " differs from an object that was initialized once: the results differ by "
                       << max_difference << ".\n");
        }
        // The result of a block preconditioner does not depend on the problem.
        if (!is_block && step > 0 && steps[step].d_vectors == steps[step - 1].d_vectors && !(max_change > 1.0e-3))
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
        plog << "velocity_min = " << range_min[0] << '\n';
        plog << "velocity_max = " << range_max[0] << '\n';
        plog << "pressure_min = " << range_min[1] << '\n';
        plog << "pressure_max = " << range_max[1] << '\n';
    }
    if (is_block)
    {
        // The last four values of a block result: the volume and the integrals.
        const std::vector<double>& last = fresh.back();
        const double volume = last[last.size() - 4], U_integral = last[last.size() - 3],
                     P_integral_before = last[last.size() - 2], P_integral_after = last[last.size() - 1];
        if (!(volume > 0.0)) TBOX_ERROR("The volume of the physical domain is " << volume << ".\n");
        if (kind == Kind::BLOCK_NULL_SPACE ? !(std::abs(P_integral_after) < 1.0e-10 * std::abs(P_integral_before)) :
                                             P_integral_after != P_integral_before)
        {
            TBOX_ERROR("The pressure integral is " << P_integral_after << " after the removal of the null space, from "
                                                   << P_integral_before << ".\n");
        }
        plog << "volume = " << volume << '\n';
        plog << "U_integral = " << U_integral << '\n';
        plog << "P_integral_before = " << P_integral_before << '\n';
        plog << "P_integral_after = " << P_integral_after << '\n';
    }
    plog << "math_ops_reused = " << measured.d_math_ops_reused << '\n';
    plog << "math_ops_built = " << measured.d_math_ops_built << '\n';
    plog << "data_ops_reused = " << measured.d_data_ops_reused << '\n';
    plog << "data_ops_built = " << measured.d_data_ops_built << '\n';
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app->getInputDatabase();
    const auto hierarchy_data = setup_hierarchy<NDIM>(app);
    Pointer<PatchHierarchy<NDIM>> hierarchy = std::get<0>(hierarchy_data);
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(hierarchy->getFinestLevelNumber());
    plog << "number_of_levels = " << hierarchy->getNumberOfLevels() << '\n';
    plog << "number_of_patches_on_finest_level = " << level->getNumberOfPatches() << '\n';

    const std::string case_name = input_db->getDatabase("test")->getString("case");
    Kind kind = Kind::STOKES;
    if (case_name == "stokes")
        kind = Kind::STOKES;
    else if (case_name == "vc_stokes")
        kind = Kind::VC_STOKES;
    else if (case_name == "block")
        kind = Kind::BLOCK;
    else if (case_name == "block_null_space")
        kind = Kind::BLOCK_NULL_SPACE;
    else if (case_name == "block_factorization")
        kind = Kind::BLOCK_FACTORIZATION;
    else if (case_name == "projection")
        kind = Kind::PROJECTION;
    else if (case_name == "vc_projection")
        kind = Kind::VC_PROJECTION;
    else
    {
        TBOX_ERROR("Unknown case \"" << case_name << "\"\n"
                                     << "  valid values are stokes, vc_stokes, block, block_null_space, "
                                        "block_factorization, projection, and vc_projection\n");
    }
    run_case(kind,
             hierarchy,
             input_db->keyExists("preconditioner_db") ? input_db->getDatabase("preconditioner_db") : nullptr);
    return 0;
}
