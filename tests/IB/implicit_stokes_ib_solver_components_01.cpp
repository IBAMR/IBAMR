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

// Nonlinear and Jacobian actions with a live IBMethod.
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/StaggeredStokesIBJacobianOperator.h>
#include <ibamr/StaggeredStokesIBOperator.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScMFFDJacobianOperator.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/ib_kernel_evaluators.h>

#include <tbox/Logger.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PatchHierarchy.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <tuple>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

namespace
{
using HierarchyVector = SAMRAIVectorReal<NDIM, double>;

class BoundaryCheckedStokesOperator : public StaggeredStokesOperator
{
public:
    BoundaryCheckedStokesOperator() : StaggeredStokesOperator("operator_test::stokes")
    {
    }
    void apply(HierarchyVector& x, HierarchyVector& y) override
    {
        ++evaluations;
        // Stokes fills the whole hierarchy even for a single-level vector.
        // Its coarsening transactions populate these coarse input values.
        std::vector<std::pair<int, int>> support;
        for (int ln = 0; ln < x.getCoarsestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = x.getPatchHierarchy()->getPatchLevel(ln);
            for (int component = 0; component < 2; ++component)
            {
                const int idx = x.getComponentDescriptorIndex(component);
                if (!level->checkAllocated(idx))
                {
                    level->allocatePatchData(idx);
                    support.emplace_back(ln, idx);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        if (component == 0)
                        {
                            Pointer<SideData<NDIM, double>> data = level->getPatch(p())->getPatchData(idx);
                            data->fillAll(0.0);
                        }
                        else
                        {
                            Pointer<CellData<NDIM, double>> data = level->getPatch(p())->getPatchData(idx);
                            data->fillAll(0.0);
                        }
                    }
                }
            }
        }
        StaggeredStokesOperator::apply(x, y);
        for (const std::pair<int, int>& allocation : support)
        {
            x.getPatchHierarchy()->getPatchLevel(allocation.first)->deallocatePatchData(allocation.second);
        }
    }
    void modifyRhsForBcs(HierarchyVector& y) override
    {
        ++rhs_calls;
        StaggeredStokesOperator::modifyRhsForBcs(y);
    }
    void imposeSolBcs(HierarchyVector& x) override
    {
        ++sol_calls;
        StaggeredStokesOperator::imposeSolBcs(x);
    }
    int evaluations = 0, rhs_calls = 0, sol_calls = 0;
};

void
generate_operator_structure(const unsigned int&, const int&, int& n, std::vector<IBTK::Point>& X, void* ctx)
{
    n = 16;
    X.resize(n);
    for (int k = 0; k < n; ++k)
    {
        const bool physical_boundary = ctx && *static_cast<bool*>(ctx);
        X[k](0) = (physical_boundary ? 0.1 : 0.5) + (physical_boundary ? 0.06 : 0.16) * std::cos(2.0 * M_PI * k / n);
        X[k](1) = 0.5 + 0.12 * std::sin(2.0 * M_PI * k / n);
    }
}

void
generate_operator_springs(
    const unsigned int&,
    const int&,
    std::multimap<int, IBRedundantInitializer::Edge>& edges,
    std::map<IBRedundantInitializer::Edge, IBRedundantInitializer::SpringSpec, IBRedundantInitializer::EdgeComp>& specs,
    void*)
{
    for (int k = 0; k < 16; ++k)
    {
        IBRedundantInitializer::Edge edge = { k, (k + 1) % 16 };
        if (edge.first > edge.second)
        {
            std::swap(edge.first, edge.second);
        }
        edges.emplace(edge.first, edge);
        IBRedundantInitializer::SpringSpec spring;
        spring.force_fcn_idx = 0;
        // Nonzero rest length makes the force Jacobian depend on the base state.
        spring.parameters = { 8.0, 0.02 };
        specs.emplace(edge, spring);
    }
}

void
set_operator_velocity(int idx, Pointer<PatchLevel<NDIM>> level, double amplitude, double phase)
{
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> data = patch->getPatchData(idx);
        Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(data->getGhostBox(), axis)); b; b++)
            {
                const int other = 1 - axis;
                const double q = geometry->getXLower()[other] +
                                 (b()(other) - patch->getBox().lower()(other) + 0.5) * geometry->getDx()[other];
                (*data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower)) =
                    amplitude * (axis == 0 ? 1.0 : -0.7) * std::sin(2.0 * M_PI * q + phase);
            }
        }
    }
}

int
run_operators(Pointer<AppInitializer> app)
{
    PetscErrorCode ierr;
    const int level_num = app->getInputDatabase()->getIntegerWithDefault("operator_level", 0);
    bool physical_boundary = app->getInputDatabase()->getBoolWithDefault("physical_boundary", false);
    const TimeSteppingType type =
        IBAMR::string_to_enum<TimeSteppingType>(app->getInputDatabase()->getString("time_stepping"));
    const double current = 0.25, dt = 0.125, next = current + dt;
    const bool midpoint = type == MIDPOINT_RULE;
    const double force_scale = type == TRAPEZOIDAL_RULE ? 0.5 : 1.0;
    const double position_scale = type == BACKWARD_EULER ? 1.0 : (midpoint ? 0.25 : 0.5);
    const double force_time = midpoint ? current + dt / 2 : next;
    // Check the time/state choices through actual positions and independently
    // composed residual/Jacobian actions below, not the implementation's table.
    bool time_valid = true;
    Pointer<IBMethod> method = new IBMethod("IBMethod", app->getComponentDatabase("IBMethod"));
    method->setUseFixedLEOperators(true);
    Pointer<IBStandardForceGen> force = new IBStandardForceGen();
    method->registerIBLagrangianForceFunction(force);
    Pointer<CartesianGridGeometry<NDIM>> geometry =
        new CartesianGridGeometry<NDIM>("CartesianGeometry", app->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", geometry);
    Pointer<StandardTagAndInitialize<NDIM>> tagger = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", method, app->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
        "GriddingAlgorithm", app->getComponentDatabase("GriddingAlgorithm"), tagger, boxes, balancer);
    Pointer<IBRedundantInitializer> initializer =
        new IBRedundantInitializer("IBRedundantInitializer", app->getComponentDatabase("IBRedundantInitializer"));
    initializer->setStructureNamesOnLevel(level_num, { "curve" });
    initializer->registerInitStructureFunction(generate_operator_structure, &physical_boundary);
    initializer->registerInitSpringDataFunction(generate_operator_springs);
    method->registerLInitStrategy(initializer);
    gridding->makeCoarsestLevel(hierarchy, current);
    for (int ln = 0; ln < level_num; ++ln)
    {
        gridding->makeFinerLevel(hierarchy, current, true, 1);
    }
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_num);
    TBOX_ASSERT(level->getLevelNumber() == level_num);
    if (IBTK_MPI::getNodes() != 1 || level->getNumberOfPatches() != 1)
    {
        TBOX_ERROR("Operator fixture requires one patch on one rank\n");
    }

    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = variables->getContext("operators");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    Pointer<SideVariable<NDIM, int>> u_dof_var = new SideVariable<NDIM, int>("u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_var = new CellVariable<NDIM, int>("p_dof");
    const IntVector<NDIM> ghosts = method->getMinimumGhostCellWidth();
    std::vector<int> allocated;
    auto register_data =
        [&](Pointer<SAMRAI::hier::Variable<NDIM>> variable, const std::string& name, IntVector<NDIM> width)
    {
        const int idx = variables->registerVariableAndContext(variable, variables->getContext(name), width);
        level->allocatePatchData(idx, current);
        allocated.push_back(idx);
        return idx;
    };
    const int u_current = register_data(u_var, "current", ghosts);
    const int scratch = register_data(u_var, "scratch", ghosts);
    const int f_scratch = register_data(u_var, "force", ghosts);
    const int u_dof = register_data(u_dof_var, "dofs", ghosts);
    const int p_dof = register_data(p_dof_var, "dofs", IntVector<NDIM>(0));
    set_operator_velocity(u_current, level, 0.03, 0.2);
    std::vector<Pointer<CoarsenSchedule<NDIM>>> synch(level_num + 1);
    std::vector<Pointer<RefineSchedule<NDIM>>> fill(level_num + 1), prolong(level_num + 1);
    std::vector<Pointer<LocationIndexRobinBcCoefs<NDIM>>> physical_coefs(NDIM);
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM, nullptr);
    Pointer<CartSideRobinPhysBdryOp> physical_bc;
    if (physical_boundary)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            physical_coefs[axis] = new LocationIndexRobinBcCoefs<NDIM>("velocity_bc", nullptr);
            for (int face = 0; face < 2 * NDIM; ++face)
            {
                physical_coefs[axis]->setBoundaryValue(face, axis == 1 && face < 2 ? 0.3 : 0.0);
            }
            bc_coefs[axis] = physical_coefs[axis];
        }
        physical_bc = new CartSideRobinPhysBdryOp(u_current, bc_coefs, false);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            physical_bc->setPhysicalBoundaryConditions(*level->getPatch(p()), current, ghosts);
        }
        method->beginDataRedistribution(hierarchy, gridding);
        method->endDataRedistribution(hierarchy, gridding);
    }
    method->initializePatchHierarchy(hierarchy, gridding, u_current, synch, fill, 0, current, true);
    method->freeLInitStrategy();
    initializer.setNull();
    method->preprocessIntegrateData(current, next, 1);
    method->updateFixedLEOperators();
    method->interpolateVelocity(u_current, synch, fill, current);
    if (physical_boundary)
    {
        physical_bc->setPatchDataIndex(scratch);
        RefineAlgorithm<NDIM> ghost_fill;
        ghost_fill.registerRefine(scratch, scratch, scratch, nullptr);
        fill[level_num] = ghost_fill.createSchedule(level, physical_bc.getPointer());
    }

    Pointer<HierarchySideDataOpsReal<NDIM, double>> side_ops =
        new HierarchySideDataOpsReal<NDIM, double>(hierarchy, level_num, level_num);
    Pointer<HierarchyCellDataOpsReal<NDIM, double>> cell_ops =
        new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, level_num, level_num);
    HierarchyMathOps math_ops("operator_test::math", hierarchy);
    const int u = variables->registerVariableAndContext(u_var, context, ghosts);
    const int p = variables->registerVariableAndContext(p_var, context, IntVector<NDIM>(1));
    Pointer<HierarchyVector> base = new HierarchyVector("base", hierarchy, level_num, level_num);
    base->addComponent(u_var, u, math_ops.getSideWeightPatchDescriptorIndex(), side_ops);
    base->addComponent(p_var, p, math_ops.getCellWeightPatchDescriptorIndex(), cell_ops);
    base->allocateVectorData();
    std::vector<Pointer<HierarchyVector>> vectors = { base };
    auto clone = [&](const std::string& name)
    {
        Pointer<HierarchyVector> vector = base->cloneVector(name);
        vector->allocateVectorData();
        vector->setToScalar(0.0);
        vectors.push_back(vector);
        return vector;
    };
    Pointer<HierarchyVector> direction = clone("direction"), residual = clone("residual"), expected = clone("expected"),
                             action = clone("action"), finite_difference = clone("finite_difference"),
                             plus = clone("plus"), minus = clone("minus"), work = clone("work"),
                             difference = clone("difference"), first_action = clone("first_action"),
                             mffd_base = clone("mffd_base"), mffd_direction = clone("mffd_direction");
    set_operator_velocity(direction->getComponentDescriptorIndex(0), level, 0.2, 0.8);
    cell_ops->setToScalar(direction->getComponentDescriptorIndex(1), -0.25);

    if (app->getInputDatabase()->getBoolWithDefault("supplied_action_only", false))
    {
        TBOX_ASSERT(level_num == 1 && direction->getCoarsestLevelNumber() == 1 &&
                    direction->getFinestLevelNumber() == 1);
        Pointer<PatchLevel<NDIM>> coarse_level = hierarchy->getPatchLevel(0);
        TBOX_ASSERT(!coarse_level->checkAllocated(u_dof));
        std::vector<int> counts;
        StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(counts, u_dof, p_dof, level);
        Mat coupling = nullptr;
        ierr = MatCreateAIJ(
            PETSC_COMM_WORLD, counts[0], counts[0], counts[0], counts[0], 1, nullptr, 0, nullptr, &coupling);
        IBTK_CHKERRQ(ierr);
        Pointer<SideData<NDIM, int>> velocity_dofs = level->getPatch(0)->getPatchData(u_dof);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(level->getPatch(0)->getBox(), axis)); b; b++)
            {
                const PetscInt row = (*velocity_dofs)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower));
                ierr = MatSetValue(coupling, row, row, 2.0, INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }
        ierr = MatAssemblyBegin(coupling, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(coupling, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        Pointer<BoundaryCheckedStokesOperator> stokes = new BoundaryCheckedStokesOperator();
        PoissonSpecifications coefs("supplied_level_coefs");
        coefs.setCConstant(2.0);
        coefs.setDConstant(-0.01);
        stokes->setVelocityPoissonSpecifications(coefs);
        StaggeredStokesIBOperator::Context ctx;
        ctx.stokes_op = stokes;
        ctx.u_dof_index_idx = u_dof;
        ctx.p_dof_index_idx = p_dof;
        StaggeredStokesIBJacobianOperator jacobian("supplied_level_jacobian");
        jacobian.setOperatorContext(ctx);
        jacobian.setTimeInterval(current, next);
        jacobian.setSolutionTime(force_time);
        jacobian.setIBCouplingJacobian(coupling);
        jacobian.initializeOperatorState(*direction, *action);
        jacobian.formJacobian(*direction);
        TBOX_ASSERT(jacobian.getBaseVector());
        jacobian.apply(*direction, *action);
        stokes->apply(*direction, *expected);
        side_ops->axpy(expected->getComponentDescriptorIndex(0),
                       2.0,
                       direction->getComponentDescriptorIndex(0),
                       expected->getComponentDescriptorIndex(0));
        difference->subtract(action, expected);
        const double error = difference->maxNorm();
        const double scale = std::max(1.0, expected->maxNorm());
        TBOX_ASSERT(std::isfinite(error) && error <= 1.0e-12 * scale);
        pout << "supplied_action_level = " << level_num << '\n' << "supplied_action_error = " << error << '\n';
        jacobian.deallocateOperatorState();
        ierr = MatDestroy(&coupling);
        IBTK_CHKERRQ(ierr);
        method->postprocessIntegrateData(current, next, 1);
        for (auto& vector : vectors)
        {
            free_vector_components(*vector);
        }
        for (int idx : allocated)
        {
            level->deallocatePatchData(idx);
            variables->removePatchDataIndex(idx);
        }
        return 0;
    }

    if (physical_boundary)
    {
        // Exercise the live strategy near the wall, without assembling an IB matrix.
        Pointer<StaggeredStokesOperator> stokes = new StaggeredStokesOperator("boundary_stokes");
        PoissonSpecifications coefs("boundary_coefs");
        coefs.setCConstant(2.0);
        coefs.setDConstant(-0.01);
        stokes->setVelocityPoissonSpecifications(coefs);
        stokes->setPhysicalBcCoefs(bc_coefs, nullptr);
        StaggeredStokesIBOperator::Context ctx;
        ctx.ib_implicit_ops = method;
        ctx.stokes_op = stokes;
        ctx.hier_velocity_data_ops = side_ops;
        ctx.u_phys_bdry_op = physical_bc;
        ctx.u_synch_scheds = synch;
        ctx.u_ghost_fill_scheds = fill;
        ctx.f_prolongation_scheds = prolong;
        ctx.u_idx = scratch;
        ctx.f_idx = f_scratch;
        ctx.u_current_idx = u_current;
        ctx.time_stepping_type = BACKWARD_EULER;
        StaggeredStokesIBOperator nonlinear("nonlinear");
        StaggeredStokesIBJacobianOperator jacobian("boundary_jacobian");
        nonlinear.setOperatorContext(ctx);
        jacobian.setOperatorContext(ctx);
        for (GeneralOperator* op :
             { static_cast<GeneralOperator*>(&nonlinear), static_cast<GeneralOperator*>(&jacobian) })
        {
            op->setTimeInterval(current, next);
            op->setSolutionTime(next);
        }
        Vec physical_position = nullptr;
        Vec X0 = method->getLDataManager()->getLData("X", level_num)->getVec();
        ierr = VecDuplicate(X0, &physical_position);
        IBTK_CHKERRQ(ierr);
        double fd_error = 0.0, boundary_position_change = 0.0, coupling_norm = 0.0;
        for (int cycle = 0; cycle < 2; ++cycle)
        {
            method->setUpdatedPosition(X0);
            nonlinear.initializeOperatorState(*base, *residual);
            jacobian.initializeOperatorState(*base, *residual);
            for (int state = 0; state < 2; ++state)
            {
                base->setToScalar(0.0);
                set_operator_velocity(u, level, state == 0 ? 0.12 : 0.24, state == 0 ? 0.1 : 0.6);
                nonlinear.apply(*base, *residual);
                std::vector<Pointer<LData>>* positions = nullptr;
                bool* needs_fill = nullptr;
                method->getPositionData(&positions, &needs_fill, TimePoint::NEW_TIME);
                ierr = VecCopy((*positions)[level_num]->getVec(), physical_position);
                IBTK_CHKERRQ(ierr);
                // A zero-boundary control proves the marker interpolation sees the wall data.
                for (int face = 0; face < 2; ++face)
                {
                    physical_coefs[1]->setBoundaryValue(face, 0.0);
                }
                nonlinear.apply(*base, *expected);
                ierr = VecAXPY(physical_position, -1.0, (*positions)[level_num]->getVec());
                IBTK_CHKERRQ(ierr);
                PetscReal position_change = 0.0;
                ierr = VecNorm(physical_position, NORM_INFINITY, &position_change);
                IBTK_CHKERRQ(ierr);
                boundary_position_change = std::max(boundary_position_change, position_change);
                for (int face = 0; face < 2; ++face)
                {
                    physical_coefs[1]->setBoundaryValue(face, 0.3);
                }
                nonlinear.apply(*base, *residual);
                jacobian.formJacobian(*base);
                jacobian.apply(*direction, *action);
                stokes->apply(*direction, *expected);
                difference->subtract(action, expected);
                const double ib_norm = difference->maxNorm();
                coupling_norm = std::max(coupling_norm, ib_norm);
                const double h = 1.0e-5;
                work->linearSum(1.0, base, h, direction);
                nonlinear.apply(*work, *plus);
                work->linearSum(1.0, base, -h, direction);
                nonlinear.apply(*work, *minus);
                finite_difference->linearSum(0.5 / h, plus, -0.5 / h, minus);
                difference->subtract(action, finite_difference);
                const double error = difference->maxNorm() / std::max(1.0, finite_difference->maxNorm());
                fd_error = std::max(fd_error, error);
                if (!std::isfinite(error) || error > 1.0e-6 || !std::isfinite(position_change) ||
                    position_change <= 1.0e-4 || !std::isfinite(ib_norm) || ib_norm <= 1.0e-4)
                {
                    TBOX_ERROR("The physical base state failed the finite-difference or nontrivial-action check.\n");
                }
            }
            jacobian.deallocateOperatorState();
            nonlinear.deallocateOperatorState();
        }
        pout << "physical_base_fd_error = " << std::scientific << std::setprecision(0) << fd_error << std::defaultfloat
             << std::setprecision(6) << '\n'
             << "boundary_position_change = " << boundary_position_change << '\n'
             << "nonzero_ib_action = " << coupling_norm << '\n';
        ierr = VecDestroy(&physical_position);
        IBTK_CHKERRQ(ierr);
        method->postprocessIntegrateData(current, next, 1);
        for (auto& vector : vectors)
        {
            free_vector_components(*vector);
        }
        for (int idx : allocated)
        {
            level->deallocatePatchData(idx);
            variables->removePatchDataIndex(idx);
        }
        return 0;
    }

    std::vector<int> counts;
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(counts, u_dof, p_dof, level);
    Mat J = nullptr, A = nullptr;
    Vec X0 = method->getLDataManager()->getLData("X", level_num)->getVec();
    const auto evaluator = IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::IB4{} };
    PETScMatUtilities::constructPatchLevelSCInterpOp(J, evaluator, X0, counts, u_dof, level);
    method->constructLagrangianForceJacobian(A, MATAIJ, force_time);
    Vec eulerian = nullptr, spread = nullptr, interpolated = nullptr, position = nullptr, expected_position = nullptr;
    ierr = MatCreateVecs(J, &eulerian, &interpolated);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(eulerian, &spread);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(interpolated, &position);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(position, &expected_position);
    IBTK_CHKERRQ(ierr);
    Pointer<LData> expected_force = new LData("expected_force", 16, NDIM);
    Pointer<CartesianPatchGeometry<NDIM>> patch_geometry = level->getPatch(0)->getPatchGeometry();
    const double cell_volume = patch_geometry->getDx()[0] * patch_geometry->getDx()[1];
    Pointer<BoundaryCheckedStokesOperator> stokes = new BoundaryCheckedStokesOperator();
    PoissonSpecifications coefs("coefs");
    coefs.setCConstant(2.0);
    coefs.setDConstant(-0.01);
    stokes->setVelocityPoissonSpecifications(coefs);
    stokes->setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr), nullptr);
    StaggeredStokesIBOperator::Context ctx;
    ctx.ib_implicit_ops = method;
    ctx.stokes_op = stokes;
    ctx.hier_velocity_data_ops = side_ops;
    ctx.u_synch_scheds = synch;
    ctx.u_ghost_fill_scheds = fill;
    ctx.f_prolongation_scheds = prolong;
    ctx.u_idx = scratch;
    ctx.f_idx = f_scratch;
    ctx.u_current_idx = u_current;
    ctx.u_dof_index_idx = u_dof;
    ctx.p_dof_index_idx = p_dof;
    ctx.time_stepping_type = type;
    StaggeredStokesIBOperator nonlinear("nonlinear");
    StaggeredStokesIBJacobianOperator jacobian("jacobian");
    PETScMFFDJacobianOperator mffd("mffd");
    const bool test_mffd = type == BACKWARD_EULER;
    nonlinear.setOperatorContext(ctx);
    jacobian.setOperatorContext(ctx);
    mffd.setOperator(stokes);
    const std::string context_check = app->getInputDatabase()->getStringWithDefault("context_check", "");
    if (!context_check.empty())
    {
        Pointer<Logger::Appender> abort_appender = new TestAppender();
        Logger::getInstance()->setAbortAppender(abort_appender);
    }
    if (context_check == "incomplete")
    {
        ctx.hier_velocity_data_ops.setNull();
        nonlinear.setOperatorContext(ctx);
        nonlinear.initializeOperatorState(*base, *residual);
        return 0;
    }
    if (context_check == "replace")
    {
        jacobian.setTimeInterval(current, next);
        jacobian.setSolutionTime(force_time);
        jacobian.initializeOperatorState(*base, *residual);
        jacobian.setOperatorContext(ctx);
        return 0;
    }
    for (GeneralOperator* op : { static_cast<GeneralOperator*>(&nonlinear),
                                 static_cast<GeneralOperator*>(&jacobian),
                                 static_cast<GeneralOperator*>(&mffd) })
    {
        op->setTimeInterval(current, next);
        op->setSolutionTime(force_time);
    }
    auto copy_to_petsc = [&](Vec v, Pointer<HierarchyVector> x)
    {
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            v, x->getComponentDescriptorIndex(0), u_dof, x->getComponentDescriptorIndex(1), p_dof, level);
    };
    auto close =
        [&](Pointer<HierarchyVector> lhs, Pointer<HierarchyVector> rhs, double tol, const char* label = nullptr)
    {
        difference->subtract(lhs, rhs);
        const double error = difference->maxNorm(), norm = rhs->maxNorm();
        const double bound = tol * std::max(1.0, norm);
        const bool valid = std::isfinite(error) && std::isfinite(norm) && error <= bound;
        if (!valid && label)
        {
            pout << label << ": error = " << std::setprecision(17) << error << ", bound = " << bound
                 << std::setprecision(6) << '\n';
        }
        return valid;
    };
    constexpr double RESIDUAL_TOL = 1.0e-11, JACOBIAN_TOL = 1.0e-9, FD_TOL = 1.0e-6, INITIALIZATION_TOL = 1.0e-12;
    bool residual_valid = true, derivative_valid = true, assembled_valid = true, base_valid = true,
         lifecycle_valid = true, boundary_valid = true, nontrivial = true;
    bool nonlinear_add_valid = true, jacobian_add_valid = true, mffd_valid = true;
    // Reinitialize all operator storage, and change the base twice per lifetime.
    for (int cycle = 0; cycle < 2; ++cycle)
    {
        ierr = VecCopy(X0, position);
        IBTK_CHKERRQ(ierr);
        // The second lifetime uses different fixed current/endpoint coupling
        // positions. Trapezoidal stepping must retain stored Lagrangian Ucurrent.
        if (cycle == 1)
        {
            ierr = VecShift(position, 0.01);
            IBTK_CHKERRQ(ierr);
        }
        method->setUpdatedPosition(position);
        nonlinear.initializeOperatorState(*base, *residual);
        jacobian.initializeOperatorState(*base, *residual);
        if (test_mffd)
        {
            mffd.initializeOperatorState(*mffd_base, *residual);
        }
        // Fixed endpoint coordinates were supplied above; midpoint coupling
        // uses their average with the unchanged current coordinates.
        ierr = VecCopy(position, expected_position);
        IBTK_CHKERRQ(ierr);
        if (midpoint)
        {
            ierr = VecAXPBY(expected_position, 0.5, 0.5, X0);
            IBTK_CHKERRQ(ierr);
        }
        PETScMatUtilities::constructPatchLevelSCInterpOp(J, evaluator, expected_position, counts, u_dof, level);
        for (int state = 0; state < 2; ++state)
        {
            base->setToScalar(0.0);
            set_operator_velocity(u, level, state == 0 ? 0.12 : 0.24, state == 0 ? 0.1 : 0.6);
            nonlinear.apply(*base, *residual);
            // Independent position/force composition using the live force law
            // and a directly assembled interpolation matrix, with no FAC helper.
            side_ops->linearSum(scratch, position_scale, u, midpoint ? position_scale : 0.0, u_current);
            side_ops->copyData(work->getComponentDescriptorIndex(0), scratch);
            cell_ops->setToScalar(work->getComponentDescriptorIndex(1), 0.0);
            copy_to_petsc(eulerian, work);
            ierr = MatMult(J, eulerian, interpolated);
            IBTK_CHKERRQ(ierr);
            ierr = VecWAXPY(position, dt, interpolated, X0);
            IBTK_CHKERRQ(ierr);
            if (type == TRAPEZOIDAL_RULE)
            {
                std::vector<Pointer<LData>>* U_current_data;
                method->getVelocityData(&U_current_data, TimePoint::CURRENT_TIME);
                ierr = VecAXPY(position, dt / 2, (*U_current_data)[level_num]->getVec());
                IBTK_CHKERRQ(ierr);
            }
            ierr = VecCopy(position, expected_position);
            IBTK_CHKERRQ(ierr);
            std::vector<Pointer<LData>>*X_data, *U_data;
            bool* X_ghost;
            const TimePoint force_point = midpoint ? TimePoint::HALF_TIME : TimePoint::NEW_TIME;
            method->getPositionData(&X_data, &X_ghost, force_point);
            method->getVelocityData(&U_data, force_point);
            ierr = VecAXPY(position, -1.0, (*X_data)[level_num]->getVec());
            IBTK_CHKERRQ(ierr);
            PetscReal position_error;
            ierr = VecNorm(position, NORM_INFINITY, &position_error);
            IBTK_CHKERRQ(ierr);
            time_valid = time_valid && std::isfinite(position_error) && position_error < 1.0e-12;
            stokes->apply(*base, *expected);
            ierr = VecSet(expected_force->getVec(), 0.0);
            IBTK_CHKERRQ(ierr);
            force->computeLagrangianForce(expected_force,
                                          (*X_data)[level_num],
                                          (*U_data)[level_num],
                                          hierarchy,
                                          level_num,
                                          force_time,
                                          method->getLDataManager());
            ierr = MatMultTranspose(J, expected_force->getVec(), spread);
            IBTK_CHKERRQ(ierr);
            copy_to_petsc(eulerian, expected);
            ierr = VecAXPY(eulerian, -force_scale / cell_volume, spread);
            IBTK_CHKERRQ(ierr);
            StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(eulerian,
                                                                    expected->getComponentDescriptorIndex(0),
                                                                    u_dof,
                                                                    expected->getComponentDescriptorIndex(1),
                                                                    p_dof,
                                                                    level,
                                                                    nullptr,
                                                                    nullptr);
            residual_valid = close(residual, expected, RESIDUAL_TOL, "nonlinear residual") && residual_valid;
            nontrivial = nontrivial && residual->maxNorm() > 1.0e-6;
            // Assemble -gamma*dt*alpha J^T A J at the actual force position.
            ierr = MatZeroEntries(A);
            IBTK_CHKERRQ(ierr);
            force->computeLagrangianForceJacobian(A,
                                                  MAT_FINAL_ASSEMBLY,
                                                  1.0,
                                                  (*X_data)[level_num],
                                                  0.0,
                                                  nullptr,
                                                  hierarchy,
                                                  level_num,
                                                  force_time,
                                                  method->getLDataManager());
            Mat coupling = nullptr;
            ierr = MatPtAP(A, J, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &coupling);
            IBTK_CHKERRQ(ierr);
            ierr = MatScale(coupling, -force_scale * dt * position_scale / cell_volume);
            IBTK_CHKERRQ(ierr);
            // Leave stale nonlinear state from a different input deliberately:
            // formJacobian must reconstruct its own supplied base.
            nonlinear.apply(*direction, *work);
            jacobian.formJacobian(*base);
            if (type != BACKWARD_EULER)
            {
                method->getPositionData(&X_data, &X_ghost, force_point);
                ierr = VecWAXPY(position, -1.0, expected_position, (*X_data)[level_num]->getVec());
                IBTK_CHKERRQ(ierr);
                ierr = VecNorm(position, NORM_INFINITY, &position_error);
                IBTK_CHKERRQ(ierr);
                time_valid = time_valid && std::isfinite(position_error) && position_error < 1.0e-12;
            }
            base_valid = close(jacobian.getBaseVector(), base, 0.0) && base_valid;
            jacobian.apply(*direction, *action);
            stokes->apply(*direction, *expected);
            nontrivial = nontrivial && !close(action, expected, 1.0e-5);
            if (state == 0)
            {
                first_action->copyVector(action);
            }
            else
            {
                nontrivial = nontrivial && !close(action, first_action, 1.0e-7);
            }
            jacobian.setIBCouplingJacobian(coupling);
            // Release the creator's reference, then reinstall the borrowed handle.
            // The operator must retain its own reference throughout replacement.
            Mat coupling_alias = coupling;
            ierr = MatDestroy(&coupling);
            IBTK_CHKERRQ(ierr);
            jacobian.setIBCouplingJacobian(coupling_alias);
            jacobian.apply(*direction, *expected);
            assembled_valid = close(action, expected, JACOBIAN_TOL, "assembled Jacobian") && assembled_valid;
            jacobian.setIBCouplingJacobian(nullptr);

            const double h = 1.0e-5;
            plus->linearSum(1.0, base, h, direction);
            minus->linearSum(1.0, base, -h, direction);
            nonlinear.apply(*plus, *expected);
            nonlinear.apply(*minus, *work);
            finite_difference->linearSum(0.5 / h, expected, -0.5 / h, work);
            derivative_valid =
                close(action, finite_difference, FD_TOL, "centered finite difference") && derivative_valid;
            if (test_mffd)
            {
                // A real Stokes action exercises SAMRAI-backed MFFD function storage.
                // Its derivative is linear: changed-base coverage comes from the
                // stored-vector and evaluation checks, not a changing derivative.
                side_ops->setToScalar(mffd_base->getComponentDescriptorIndex(0), state == 0 ? 0.5 : 1.0);
                cell_ops->setToScalar(mffd_base->getComponentDescriptorIndex(1), state == 0 ? -0.25 : -0.5);
                side_ops->setToScalar(mffd_direction->getComponentDescriptorIndex(0), 0.25);
                cell_ops->setToScalar(mffd_direction->getComponentDescriptorIndex(1), 0.5);
                stokes->setTimeInterval(current, next);
                stokes->setSolutionTime(force_time);
                stokes->setHomogeneousBc(true);
                const int evaluations_before_form = stokes->evaluations;
                mffd.formJacobian(*mffd_base);
                base_valid = base_valid && stokes->evaluations == evaluations_before_form + 1;
                base_valid = close(mffd.getBaseVector(), mffd_base, 0.0) && base_valid;
                mffd.apply(*mffd_direction, *expected);
                stokes->apply(*mffd_direction, *work);
                mffd_valid =
                    close(expected, work, FD_TOL, "MFFD Stokes action") && work->maxNorm() > 1.0e-12 && mffd_valid;
            }

            nonlinear.applyAdd(*base, *direction, *expected);
            work->add(residual, direction);
            nonlinear_add_valid = close(expected, work, RESIDUAL_TOL, "nonlinear applyAdd") && nonlinear_add_valid;
            expected->copyVector(direction);
            nonlinear.applyAdd(*base, *expected, *expected);
            nonlinear_add_valid = close(expected, work, RESIDUAL_TOL, "nonlinear applyAdd y=z") && nonlinear_add_valid;
            expected->copyVector(base);
            nonlinear.applyAdd(*expected, *direction, *expected);
            nonlinear_add_valid = close(expected, work, RESIDUAL_TOL, "nonlinear applyAdd x=z") && nonlinear_add_valid;
            jacobian.applyAdd(*direction, *base, *expected);
            work->add(action, base);
            jacobian_add_valid = close(expected, work, JACOBIAN_TOL, "Jacobian applyAdd") && jacobian_add_valid;
            expected->copyVector(base);
            jacobian.applyAdd(*direction, *expected, *expected);
            jacobian_add_valid = close(expected, work, JACOBIAN_TOL, "Jacobian applyAdd y=z") && jacobian_add_valid;
            expected->copyVector(direction);
            jacobian.applyAdd(*expected, *base, *expected);
            jacobian_add_valid = close(expected, work, JACOBIAN_TOL, "Jacobian applyAdd x=z") && jacobian_add_valid;
        }
        for (GeneralOperator* op :
             { static_cast<GeneralOperator*>(&nonlinear), static_cast<GeneralOperator*>(&jacobian) })
        {
            for (bool homogeneous : { false, true })
            {
                op->setHomogeneousBc(homogeneous);
                op->setSolutionTime(current + 0.03125);
                stokes->setTimeInterval(-2.0, -1.0);
                stokes->setSolutionTime(-1.0);
                stokes->setHomogeneousBc(!homogeneous);
                op->modifyRhsForBcs(*residual);
                boundary_valid = boundary_valid && stokes->getTimeInterval() == std::make_pair(current, next) &&
                                 stokes->getSolutionTime() == current + 0.03125 &&
                                 stokes->getHomogeneousBc() == homogeneous;
                stokes->setTimeInterval(-2.0, -1.0);
                stokes->setSolutionTime(-1.0);
                stokes->setHomogeneousBc(!homogeneous);
                op->imposeSolBcs(*base);
                boundary_valid = boundary_valid && stokes->getTimeInterval() == std::make_pair(current, next) &&
                                 stokes->getSolutionTime() == current + 0.03125 &&
                                 stokes->getHomogeneousBc() == homogeneous;
                op->setSolutionTime(force_time);
            }
        }
        std::vector<int> base_indices = { jacobian.getBaseVector()->getComponentDescriptorIndex(0),
                                          jacobian.getBaseVector()->getComponentDescriptorIndex(1) };
        if (test_mffd)
        {
            base_indices.push_back(mffd.getBaseVector()->getComponentDescriptorIndex(0));
            base_indices.push_back(mffd.getBaseVector()->getComponentDescriptorIndex(1));
            mffd.deallocateOperatorState();
        }
        jacobian.deallocateOperatorState();
        nonlinear.deallocateOperatorState();
        lifecycle_valid = lifecycle_valid && !mffd.getIsInitialized() && !jacobian.getIsInitialized() &&
                          !nonlinear.getIsInitialized() && !jacobian.getBaseVector() && !mffd.getBaseVector();
        for (const int idx : base_indices)
        {
            Pointer<SAMRAI::hier::Variable<NDIM>> variable;
            lifecycle_valid = !variables->mapIndexToVariable(idx, variable) && lifecycle_valid;
        }
        for (const int idx : { u, p, u_current, scratch, f_scratch })
        {
            Pointer<SAMRAI::hier::Variable<NDIM>> variable;
            lifecycle_valid = variables->mapIndexToVariable(idx, variable) && level->getPatch(0)->getPatchData(idx) &&
                              lifecycle_valid;
        }
    }
    // A zero supplied contribution distinguishes its action from the nonzero
    // strategy Jacobian without inspecting which implementation is selected.
    Mat zero_coupling = nullptr;
    ierr = MatPtAP(A, J, MAT_INITIAL_MATRIX, 1.0, &zero_coupling);
    IBTK_CHKERRQ(ierr);
    ierr = MatZeroEntries(zero_coupling);
    IBTK_CHKERRQ(ierr);
    jacobian.initializeOperatorState(*base, *residual);
    jacobian.setIBCouplingJacobian(zero_coupling);
    jacobian.apply(*direction, *first_action);
    stokes->apply(*direction, *work);
    const bool initial_supplied_valid = close(first_action, work, INITIALIZATION_TOL, "initial supplied action");

    PETScKrylovLinearSolver outer("initialization_sequence", nullptr, "initialization_sequence_");
    outer.setOperator(Pointer<LinearOperator>(&jacobian, false));
    outer.setTimeInterval(current, next);
    outer.setSolutionTime(force_time);
    outer.initializeSolverState(*base, *residual);
    lifecycle_valid = !jacobian.getBaseVector() && lifecycle_valid;
    // The installed matrix survives the outer solver's reinitialization, and
    // forming the Jacobian with it leaves the strategy untouched.
    jacobian.apply(*direction, *expected);
    const bool retained_valid = close(expected, first_action, INITIALIZATION_TOL, "retained supplied action");
    jacobian.formJacobian(*base);
    lifecycle_valid = jacobian.getBaseVector() && lifecycle_valid;
    jacobian.apply(*direction, *expected);
    const bool supplied_form_valid =
        close(expected, first_action, INITIALIZATION_TOL, "supplied action after formJacobian");

    // Clearing the matrix selects the strategy action.
    jacobian.setIBCouplingJacobian(nullptr);
    jacobian.formJacobian(*base);
    jacobian.apply(*direction, *action);
    difference->subtract(action, first_action);
    const double selection_difference = difference->maxNorm();
    nontrivial = !close(action, first_action, 1.0e-5) && nontrivial;

    // Reinstalling the matrix restores the zero-coupling (Stokes-only) action.
    jacobian.setIBCouplingJacobian(zero_coupling);
    jacobian.apply(*direction, *expected);
    const bool supplied_valid = close(expected, first_action, INITIALIZATION_TOL, "reinstalled supplied action");
    outer.deallocateSolverState();
    lifecycle_valid = !jacobian.getIsInitialized() && !jacobian.getBaseVector() && lifecycle_valid;
    ierr = MatDestroy(&zero_coupling);
    IBTK_CHKERRQ(ierr);
    pout << "strategy_supplied_action_difference = " << selection_difference << '\n';
    boundary_valid = boundary_valid && stokes->rhs_calls == 8 && stokes->sol_calls == 8;
    // Check operator accuracy, not equality of cancellation-sensitive errors.
    pout << "Accuracy checks use error_inf <= tolerance * max(1, reference_inf); attained roundoff may vary.\n";
    for (const auto& check : { std::make_tuple("nonlinear_residual", residual_valid, RESIDUAL_TOL),
                               std::make_tuple("assembled_jacobian", assembled_valid, JACOBIAN_TOL),
                               std::make_tuple("centered_fd", derivative_valid, FD_TOL),
                               std::make_tuple("nonlinear_apply_add", nonlinear_add_valid, RESIDUAL_TOL),
                               std::make_tuple("jacobian_apply_add", jacobian_add_valid, JACOBIAN_TOL),
                               std::make_tuple("initial_supplied_action", initial_supplied_valid, INITIALIZATION_TOL),
                               std::make_tuple("retained_supplied_action", retained_valid, INITIALIZATION_TOL),
                               std::make_tuple("supplied_form_jacobian", supplied_form_valid, INITIALIZATION_TOL),
                               std::make_tuple("reinstalled_supplied_action", supplied_valid, INITIALIZATION_TOL) })
    {
        pout << std::get<0>(check) << " = " << (std::get<1>(check) ? "true" : "false")
             << ", tolerance = " << std::get<2>(check) << '\n';
        if (!std::get<1>(check))
        {
            TBOX_ERROR("Failed check: " << std::get<0>(check) << ".\n");
        }
    }
    if (test_mffd)
    {
        pout << "mffd_stokes_action = " << (mffd_valid ? "true" : "false") << ", tolerance = " << FD_TOL << '\n';
        if (!mffd_valid)
        {
            TBOX_ERROR("Failed check: mffd_stokes_action.\n");
        }
    }
    for (const auto& check :
         std::vector<std::pair<std::string, bool>>{ { "time_state_scaling_valid", time_valid },
                                                    { "nontrivial_coupling_valid", nontrivial },
                                                    { "updated_base_state_valid", base_valid },
                                                    { "boundary_forwarding_valid", boundary_valid },
                                                    { "operator_lifecycle_valid", lifecycle_valid } })
    {
        pout << check.first << " = " << (check.second ? "true" : "false") << '\n';
        if (!check.second)
        {
            TBOX_ERROR("Failed check: " << check.first << ".\n");
        }
    }
    method->postprocessIntegrateData(current, next, 1);
    for (Vec* v : { &eulerian, &spread, &interpolated, &position, &expected_position })
    {
        ierr = VecDestroy(v);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatDestroy(&A);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&J);
    IBTK_CHKERRQ(ierr);
    for (auto& vector : vectors)
    {
        free_vector_components(*vector);
    }
    for (int idx : allocated)
    {
        level->deallocatePatchData(idx);
        variables->removePatchDataIndex(idx);
    }
    return 0;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Logger::getInstance()->setWarning(false);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "components.log");
    return run_operators(app);
}
