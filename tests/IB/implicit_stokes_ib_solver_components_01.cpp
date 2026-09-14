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
#include <ibamr/StaggeredStokesIBJacobianFACPreconditioner.h>
#include <ibamr/StaggeredStokesIBJacobianOperator.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesIBOperator.h>
#include <ibamr/StaggeredStokesLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CCPoissonPETScLevelSolver.h>
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
#include <ibtk/SCPoissonPETScLevelSolver.h>
#include <ibtk/ib_kernels.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/Logger.h>
#include <tbox/MemoryDatabase.h>

#if defined(PETSC_USE_LOG)
#include <petsclog.h>
#endif

#include <BergerRigoutsos.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PatchHierarchy.h>
#include <ProcessorMapping.h>
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
#include <limits>
#include <map>
#include <set>
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
        int failures = 0;
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
                failures += !std::isfinite(error) || error > 1.0e-6 || !std::isfinite(position_change) ||
                            position_change <= 1.0e-4 || !std::isfinite(ib_norm) || ib_norm <= 1.0e-4;
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
        pout << "test_failures = " << failures << std::endl;
        return failures;
    }

    std::vector<int> counts;
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(counts, u_dof, p_dof, level);
    Mat J = nullptr, A = nullptr;
    Vec X0 = method->getLDataManager()->getLData("X", level_num)->getVec();
    const auto evaluator = IBKernelEvaluatorTensorProduct{ IBKernels::IB4{} };
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
    jacobian.formJacobian(*base);
    jacobian.apply(*direction, *action);
    jacobian.setIBCouplingJacobian(nullptr);
    jacobian.apply(*direction, *expected);
    difference->subtract(action, first_action);
    const double selection_difference = difference->maxNorm();
    const bool strategy_valid = close(action, expected, INITIALIZATION_TOL, "outer initialization strategy");
    nontrivial = !close(action, first_action, 1.0e-5) && nontrivial;

    // Install the supplied matrix after the outer solver has initialized its
    // operator. This must restore the zero-coupling (Stokes-only) action.
    jacobian.setIBCouplingJacobian(zero_coupling);
    jacobian.apply(*direction, *expected);
    const bool supplied_valid =
        close(expected, first_action, INITIALIZATION_TOL, "post-initialization supplied action");
    outer.deallocateSolverState();
    lifecycle_valid = !jacobian.getIsInitialized() && !jacobian.getBaseVector() && lifecycle_valid;
    ierr = MatDestroy(&zero_coupling);
    IBTK_CHKERRQ(ierr);
    pout << "outer_initialization_action_change = " << selection_difference << '\n';
    boundary_valid = boundary_valid && stokes->rhs_calls == 8 && stokes->sol_calls == 8;
    // Check operator accuracy, not equality of cancellation-sensitive errors.
    pout << "Accuracy checks use error_inf <= tolerance * max(1, reference_inf); attained roundoff may vary.\n";
    int failures = 0;
    for (const auto& check : { std::make_tuple("nonlinear_residual", residual_valid, RESIDUAL_TOL),
                               std::make_tuple("assembled_jacobian", assembled_valid, JACOBIAN_TOL),
                               std::make_tuple("centered_fd", derivative_valid, FD_TOL),
                               std::make_tuple("nonlinear_apply_add", nonlinear_add_valid, RESIDUAL_TOL),
                               std::make_tuple("jacobian_apply_add", jacobian_add_valid, JACOBIAN_TOL),
                               std::make_tuple("initial_supplied_action", initial_supplied_valid, INITIALIZATION_TOL),
                               std::make_tuple("outer_initialization_strategy", strategy_valid, INITIALIZATION_TOL),
                               std::make_tuple("post_initialization_supplied", supplied_valid, INITIALIZATION_TOL) })
    {
        pout << std::get<0>(check) << " = " << (std::get<1>(check) ? "true" : "false")
             << ", tolerance = " << std::get<2>(check) << '\n';
        failures += !std::get<1>(check);
    }
    if (test_mffd)
    {
        pout << "mffd_stokes_action = " << (mffd_valid ? "true" : "false") << ", tolerance = " << FD_TOL << '\n';
        failures += !mffd_valid;
    }
    for (const auto& check :
         std::vector<std::pair<std::string, bool>>{ { "time_state_scaling_valid", time_valid },
                                                    { "nontrivial_coupling_valid", nontrivial },
                                                    { "updated_base_state_valid", base_valid },
                                                    { "boundary_forwarding_valid", boundary_valid },
                                                    { "operator_lifecycle_valid", lifecycle_valid } })
    {
        pout << check.first << " = " << (check.second ? "true" : "false") << '\n';
        failures += !check.second;
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
    pout << "test_failures = " << failures << std::endl;
    return failures;
}
} // namespace

namespace
{
// Expose existing protected state for lifecycle assertions, without changing
// the production API or replacing any solver operation.
template <class Solver>
class LevelSolverProbe : public Solver
{
public:
    LevelSolverProbe(const std::string& name, Pointer<Database> db) : Solver(name, db, "")
    {
    }
    Mat matrixBeforeKSP() const
    {
        return d_matrix_before_ksp;
    }
    PetscInt referencesBeforeKSP() const
    {
        return d_references_before_ksp;
    }
    bool shellStorageEmpty() const
    {
        return this->d_sub_x.empty() && this->d_sub_y.empty() && this->d_sub_ksp.empty() &&
               this->d_restriction.empty() && this->d_prolongation.empty();
    }
    std::vector<Vec> retainShellVectors()
    {
        std::vector<Vec> result = this->d_sub_x;
        result.insert(result.end(), this->d_sub_y.begin(), this->d_sub_y.end());
        for (Vec v : result)
        {
            PetscErrorCode ierr = PetscObjectReference(reinterpret_cast<PetscObject>(v));
            IBTK_CHKERRQ(ierr);
        }
        return result;
    }

protected:
    void initializeSolverStateSpecialized(const HierarchyVector& x, const HierarchyVector& b) override
    {
        Solver::initializeSolverStateSpecialized(x, b);
        d_matrix_before_ksp = this->d_petsc_mat;
        PetscErrorCode ierr =
            PetscObjectGetReference(reinterpret_cast<PetscObject>(d_matrix_before_ksp), &d_references_before_ksp);
        IBTK_CHKERRQ(ierr);
    }

private:
    Mat d_matrix_before_ksp = nullptr;
    PetscInt d_references_before_ksp = 0;
};

struct LevelFixture
{
    Pointer<PatchHierarchy<NDIM>> hierarchy;
    Pointer<PatchLevel<NDIM>> level;
    Pointer<HierarchyVector> x, b;
    std::vector<int> indices;
    std::vector<PetscInt> velocity_ids;
    std::vector<int> dof_counts, velocity_counts;
    int full_size = 0;

    explicit LevelFixture(Pointer<Database> geometry_db, int ln = 0, bool full = true, bool distributed = false)
    {
        if (IBTK_MPI::getNodes() != (distributed ? 2 : 1))
        {
            TBOX_ERROR("Unexpected rank count for level fixture\n");
        }
        Pointer<CartesianGridGeometry<NDIM>> geometry = new CartesianGridGeometry<NDIM>("level_geometry", geometry_db);
        hierarchy = new PatchHierarchy<NDIM>("level_hierarchy", geometry);
        const int ranks = IBTK_MPI::getNodes();
        BoxArray<NDIM> boxes(ranks);
        ProcessorMapping mapping(ranks);
        for (int rank = 0; rank < ranks; ++rank)
        {
            SAMRAI::hier::Index<NDIM> lower(0), upper(15);
            lower(0) = rank * 16 / ranks;
            upper(0) = (rank + 1) * 16 / ranks - 1;
            boxes[rank] = Box<NDIM>(lower, upper);
            mapping.setProcessorAssignment(rank, rank);
        }
        hierarchy->makeNewPatchLevel(0, IntVector<NDIM>(1), boxes, mapping);
        if (ln == 1)
        {
            boxes[0] = full ? Box<NDIM>(SAMRAI::hier::Index<NDIM>(0), SAMRAI::hier::Index<NDIM>(31)) :
                              Box<NDIM>(SAMRAI::hier::Index<NDIM>(8), SAMRAI::hier::Index<NDIM>(23));
            hierarchy->makeNewPatchLevel(1, IntVector<NDIM>(2), boxes, mapping);
        }
        level = hierarchy->getPatchLevel(ln);
        VariableDatabase<NDIM>* db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = db->getContext("level_fixture");
        Pointer<SideVariable<NDIM, double>> u = new SideVariable<NDIM, double>("level_u");
        Pointer<CellVariable<NDIM, double>> p = new CellVariable<NDIM, double>("level_p");
        if (db->checkVariableExists("level_u"))
        {
            u = db->getVariable("level_u");
        }
        if (db->checkVariableExists("level_p"))
        {
            p = db->getVariable("level_p");
        }
        const int ui = db->registerVariableAndContext(u, context, IntVector<NDIM>(1));
        const int pi = db->registerVariableAndContext(p, context, IntVector<NDIM>(1));
        x = new HierarchyVector("level_x", hierarchy, ln, ln);
        x->addComponent(u, ui, -1, new HierarchySideDataOpsReal<NDIM, double>(hierarchy, ln, ln));
        x->addComponent(p, pi, -1, new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, ln, ln));
        x->allocateVectorData();
        x->setToScalar(0.0);
        b = x->cloneVector("level_b");
        b->allocateVectorData();
        b->setToScalar(0.0);
        Pointer<SideVariable<NDIM, int>> ud = new SideVariable<NDIM, int>("level_ud");
        Pointer<CellVariable<NDIM, int>> pd = new CellVariable<NDIM, int>("level_pd");
        if (db->checkVariableExists("level_ud"))
        {
            ud = db->getVariable("level_ud");
        }
        if (db->checkVariableExists("level_pd"))
        {
            pd = db->getVariable("level_pd");
        }
        const int udi = db->registerVariableAndContext(ud, context, IntVector<NDIM>(1));
        const int pdi = db->registerVariableAndContext(pd, context, IntVector<NDIM>(1));
        indices = { udi, pdi };
        for (int idx : indices)
        {
            level->allocatePatchData(idx);
        }
        StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(dof_counts, udi, pdi, level);
        const int rank = IBTK_MPI::getRank();
        int first_dof = 0;
        for (int r = 0; r < ranks; ++r)
        {
            full_size += dof_counts[r];
            if (r < rank)
            {
                first_dof += dof_counts[r];
            }
        }
        std::set<int> velocity;
        Pointer<SideData<NDIM, int>> data = level->getPatch(rank)->getPatchData(udi);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator i(SideGeometry<NDIM>::toSideBox(level->getPatch(rank)->getBox(), axis)); i; i++)
            {
                const int id = (*data)(SideIndex<NDIM>(i(), axis, SideIndex<NDIM>::Lower));
                if (id >= first_dof && id < first_dof + dof_counts[rank])
                {
                    velocity.insert(id);
                }
            }
        }
        velocity_ids.assign(velocity.begin(), velocity.end());
        const int local_velocity = velocity_ids.size();
        velocity_counts.resize(ranks);
        int ierr = MPI_Allgather(&local_velocity, 1, MPI_INT, velocity_counts.data(), 1, MPI_INT, PETSC_COMM_WORLD);
        IBTK_CHKERRQ(ierr);
        std::vector<int> offsets(ranks, 0);
        for (int r = 1; r < ranks; ++r)
        {
            offsets[r] = offsets[r - 1] + velocity_counts[r - 1];
        }
        std::vector<PetscInt> all_velocity(offsets.back() + velocity_counts.back());
        ierr = MPI_Allgatherv(velocity_ids.data(),
                              local_velocity,
                              MPIU_INT,
                              all_velocity.data(),
                              velocity_counts.data(),
                              offsets.data(),
                              MPIU_INT,
                              PETSC_COMM_WORLD);
        IBTK_CHKERRQ(ierr);
        velocity_ids = std::move(all_velocity);
    }
    ~LevelFixture()
    {
        free_vector_components(*b);
        free_vector_components(*x);
        for (int idx : indices)
        {
            level->deallocatePatchData(idx);
            VariableDatabase<NDIM>::getDatabase()->removePatchDataIndex(idx);
        }
    }
};

Pointer<Database>
level_solver_database(const std::string& pc = "none", int overlap = 0)
{
    Pointer<Database> db = new MemoryDatabase("level_solver");
    db->putString("ksp_type", "gmres");
    db->putString("options_prefix", "level_");
    db->putString("pc_type", pc);
    db->putString("shell_pc_type", "additive");
    db->putBool("initial_guess_nonzero", false);
    db->putDouble("rel_residual_tol", 1.0e-12);
    db->putInteger("max_iterations", 100);
    const int size[NDIM] = { 4, 8 }, width[NDIM] = { overlap, overlap };
    db->putIntegerArray("subdomain_box_size", size, NDIM);
    db->putIntegerArray("subdomain_overlap_size", width, NDIM);
    return db;
}

Mat
level_test_matrix(PetscInt n, double shift)
{
    Mat matrix;
    PetscErrorCode ierr = MatCreateAIJ(PETSC_COMM_WORLD, n, n, n, n, 2, nullptr, 0, nullptr, &matrix);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = 0; i < n; ++i)
    {
        const PetscInt cols[2] = { i, (i + 3) % n };
        const PetscScalar vals[2] = { shift + 0.001 * i, 0.125 };
        ierr = MatSetValues(matrix, 1, &i, 2, cols, vals, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return matrix;
}

// Maximum errors across the matrix replacements and solver lifetimes in one case.
double max_matrix_error = 0.0, max_level_solve_error = 0.0, max_mapping_error = 0.0;

bool
matrices_equal(Mat a, Mat b)
{
    Mat diff;
    PetscErrorCode ierr = MatDuplicate(a, MAT_COPY_VALUES, &diff);
    IBTK_CHKERRQ(ierr);
    ierr = MatAXPY(diff, -1.0, b, DIFFERENT_NONZERO_PATTERN);
    IBTK_CHKERRQ(ierr);
    PetscReal norm;
    ierr = MatNorm(diff, NORM_INFINITY, &norm);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&diff);
    IBTK_CHKERRQ(ierr);
    max_matrix_error = std::max(max_matrix_error, norm);
    return std::isfinite(norm) && norm < 1.0e-12;
}

// Exercise the installed solver KSP and its actual configured matrix/PC.
bool
check_level_solve(PETScLevelSolver& solver)
{
    Mat matrix;
    PetscErrorCode ierr = KSPGetOperators(solver.getPETScKSP(), &matrix, nullptr);
    IBTK_CHKERRQ(ierr);
    Vec exact, rhs, solution;
    ierr = MatCreateVecs(matrix, &exact, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(exact, &solution);
    IBTK_CHKERRQ(ierr);
    PetscInt first, last;
    ierr = VecGetOwnershipRange(exact, &first, &last);
    IBTK_CHKERRQ(ierr);
    PetscScalar* values;
    ierr = VecGetArray(exact, &values);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = first; i < last; ++i)
    {
        values[i - first] = std::sin(0.13 * i) + 0.5;
    }
    ierr = VecRestoreArray(exact, &values);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(matrix, exact, rhs);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSolve(solver.getPETScKSP(), rhs, solution);
    IBTK_CHKERRQ(ierr);
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(solver.getPETScKSP(), &reason);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(solution, -1.0, exact);
    IBTK_CHKERRQ(ierr);
    PetscReal error;
    ierr = VecNorm(solution, NORM_INFINITY, &error);
    IBTK_CHKERRQ(ierr);
    for (Vec* v : { &exact, &rhs, &solution })
    {
        ierr = VecDestroy(v);
        IBTK_CHKERRQ(ierr);
    }
    max_level_solve_error = std::max(max_level_solve_error, error);
    return reason > 0 && std::isfinite(error) && error < 1.0e-9;
}

bool
check_stokes_vector_mapping(StaggeredStokesPETScLevelSolver& solver, LevelFixture& fixture)
{
    const int u = fixture.x->getComponentDescriptorIndex(0), p = fixture.x->getComponentDescriptorIndex(1);
    set_operator_velocity(u, fixture.level, 0.3, 0.4);
    Pointer<CellData<NDIM, double>> pressure = fixture.level->getPatch(0)->getPatchData(p);
    pressure->fill(0.25);
    Mat matrix;
    PetscErrorCode ierr = KSPGetOperators(solver.getPETScKSP(), &matrix, nullptr);
    IBTK_CHKERRQ(ierr);
    Vec exact, rhs, result;
    ierr = MatCreateVecs(matrix, &exact, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(exact, &result);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        exact, u, fixture.indices[0], p, fixture.indices[1], fixture.level);
    ierr = MatMult(matrix, exact, rhs);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(rhs,
                                                            fixture.b->getComponentDescriptorIndex(0),
                                                            fixture.indices[0],
                                                            fixture.b->getComponentDescriptorIndex(1),
                                                            fixture.indices[1],
                                                            fixture.level,
                                                            nullptr,
                                                            nullptr);
    fixture.x->setToScalar(0.0);
    const bool converged = solver.solveSystem(*fixture.x, *fixture.b);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        result, u, fixture.indices[0], p, fixture.indices[1], fixture.level);
    ierr = VecAXPY(result, -1.0, exact);
    IBTK_CHKERRQ(ierr);
    PetscReal error;
    ierr = VecNorm(result, NORM_INFINITY, &error);
    IBTK_CHKERRQ(ierr);
    for (Vec* v : { &exact, &rhs, &result })
    {
        ierr = VecDestroy(v);
        IBTK_CHKERRQ(ierr);
    }
    max_mapping_error = std::max(max_mapping_error, error);
    return converged && std::isfinite(error) && error < 1.0e-9;
}

PetscInt
matrix_references(Mat matrix)
{
    PetscInt references;
    PetscErrorCode ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(matrix), &references);
    IBTK_CHKERRQ(ierr);
    return references;
}

bool
check_matrix_reference_lifetime(LevelFixture& fixture)
{
    Mat a = level_test_matrix(fixture.full_size, 4.0), b = level_test_matrix(fixture.full_size, 5.0);
    bool valid = true;
    {
        StaggeredStokesPETScLevelSolver solver("uninitialized_lifetime", level_solver_database(), "");
        // Each setter owns a reference, even when both inputs are the same Mat.
        solver.setOperatorMat(a);
        solver.setAugmentedOperatorMat(a);
        valid = matrix_references(a) == 3 && valid;
        solver.setOperatorMat(a);
        solver.setAugmentedOperatorMat(a);
        valid = matrix_references(a) == 3 && valid;
        solver.setOperatorMat(b);
        valid = matrix_references(a) == 2 && matrix_references(b) == 2 && valid;
        solver.setAugmentedOperatorMat(b);
        valid = matrix_references(a) == 1 && matrix_references(b) == 3 && valid;
        solver.setOperatorMat(nullptr);
        solver.setOperatorMat(nullptr);
        valid = matrix_references(b) == 2 && valid;
        solver.setOperatorMat(a);
        valid = matrix_references(a) == 2 && valid;
        // Destruction must release inputs even without initialization.
    }
    valid = matrix_references(a) == 1 && matrix_references(b) == 1 && valid;
    {
        StaggeredStokesPETScLevelSolver solver("initialized_lifetime", level_solver_database(), "");
        solver.setTimeInterval(0.0, 1.0);
        solver.setSolutionTime(1.0);
        solver.setOperatorMat(a);
        solver.setAugmentedOperatorMat(b);
        solver.initializeSolverState(*fixture.x, *fixture.b);
        valid = check_level_solve(solver) && valid;
        // Destruction also tears down KSP state before releasing both inputs.
    }
    valid = matrix_references(a) == 1 && matrix_references(b) == 1 && valid;
    PetscErrorCode ierr = MatDestroy(&a);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&b);
    IBTK_CHKERRQ(ierr);
    return valid;
}

int
run_level_operator(Pointer<AppInitializer> app, bool augmentation)
{
    LevelFixture fixture(app->getComponentDatabase("CartesianGeometry"));
    LevelSolverProbe<StaggeredStokesPETScLevelSolver> solver("supplied_stokes", level_solver_database());
    solver.setTimeInterval(0.0, 1.0);
    solver.setSolutionTime(1.0);
    bool identity = true, creator_valid = true, values_valid = true, solves = true;
    bool references_valid = augmentation || check_matrix_reference_lifetime(fixture);
    for (int cycle = 0; cycle < (augmentation ? 4 : 2); ++cycle)
    {
        const bool full_augmentation = cycle % 2 == 0;
        Mat creator = level_test_matrix(fixture.full_size, 4.0 + cycle), original;
        PetscErrorCode ierr = MatDuplicate(creator, MAT_COPY_VALUES, &original);
        IBTK_CHKERRQ(ierr);
        solver.setOperatorMat(creator);
        references_valid = matrix_references(creator) == 2 && references_valid;
        Mat augmented = nullptr, expected = nullptr, augmented_original = nullptr;
        if (augmentation)
        {
            const PetscInt n = full_augmentation ? fixture.full_size : fixture.velocity_ids.size();
            augmented = level_test_matrix(n, 1.0);
            ierr = MatDuplicate(augmented, MAT_COPY_VALUES, &augmented_original);
            IBTK_CHKERRQ(ierr);
            ierr = MatDuplicate(creator, MAT_COPY_VALUES, &expected);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetOption(expected, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
            IBTK_CHKERRQ(ierr);
            // Independently embed known entries, using actual coupled DOF data.
            for (PetscInt row = 0; row < n; ++row)
            {
                const PetscInt next = (row + 3) % n;
                const PetscInt full_row = full_augmentation ? row : fixture.velocity_ids[row];
                const PetscInt cols[2] = { full_row, full_augmentation ? next : fixture.velocity_ids[next] };
                const PetscScalar vals[2] = { 1.0 + 0.001 * row, 0.125 };
                ierr = MatSetValues(expected, 1, &full_row, 2, cols, vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            ierr = MatAssemblyBegin(expected, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(expected, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
        }
        solver.setAugmentedOperatorMat(augmented);
        if (augmentation)
        {
            references_valid = matrix_references(augmented) == 2 && references_valid;
        }
        // Reject missing retention before testing caller release, so that a
        // regression reports failure instead of dereferencing a dangling Mat.
        if (!references_valid)
        {
            TBOX_ERROR("Installed matrix references were not retained correctly.\n");
        }
        const Mat operator_alias = creator, augmentation_alias = augmented;
        ierr = MatDestroy(&creator);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&augmented);
        IBTK_CHKERRQ(ierr);
        for (int lifetime = 0; lifetime < 2; ++lifetime)
        {
            // No caller-owned reference remains, even before initialization.
            references_valid = matrix_references(operator_alias) == 1 && references_valid;
            if (augmentation)
            {
                references_valid = matrix_references(augmentation_alias) == 1 && references_valid;
            }
            solver.initializeSolverState(*fixture.x, *fixture.b);
            Mat installed;
            ierr = KSPGetOperators(solver.getPETScKSP(), &installed, nullptr);
            IBTK_CHKERRQ(ierr);
            identity = identity && installed == solver.matrixBeforeKSP() &&
                       (augmentation ? installed != operator_alias : installed == operator_alias) &&
                       solver.referencesBeforeKSP() == 1;
            values_valid = matrices_equal(installed, augmentation ? expected : original) && values_valid;
            creator_valid = matrices_equal(operator_alias, original) && creator_valid;
            if (augmentation)
            {
                creator_valid = matrices_equal(augmentation_alias, augmented_original) && creator_valid;
            }
            solves = check_level_solve(solver) && solves;
            solves = check_stokes_vector_mapping(solver, fixture) && solves;
            solver.deallocateSolverState();
            // Same-handle replacement must also work with only the solver's
            // reference remaining. Aliases are used read-only throughout.
            solver.setOperatorMat(operator_alias);
            solver.setAugmentedOperatorMat(augmentation_alias);
            references_valid = matrix_references(operator_alias) == 1 && references_valid;
            if (augmentation)
            {
                references_valid = matrix_references(augmentation_alias) == 1 && references_valid;
            }
        }
        // Retain observer references to check that clearing releases exactly
        // the solver's references, without destroying another owner's data.
        creator = operator_alias;
        augmented = augmentation_alias;
        ierr = PetscObjectReference(reinterpret_cast<PetscObject>(creator));
        IBTK_CHKERRQ(ierr);
        ierr = PetscObjectReference(reinterpret_cast<PetscObject>(augmented));
        IBTK_CHKERRQ(ierr);
        solver.setOperatorMat(nullptr);
        solver.setAugmentedOperatorMat(nullptr);
        references_valid = matrix_references(creator) == 1 && references_valid;
        creator_valid = matrices_equal(creator, original) && creator_valid;
        if (augmentation)
        {
            references_valid = matrix_references(augmented) == 1 && references_valid;
            creator_valid = matrices_equal(augmented, augmented_original) && creator_valid;
        }
        for (Mat* m : { &creator, &original, &augmented, &expected, &augmented_original })
        {
            ierr = MatDestroy(m);
            IBTK_CHKERRQ(ierr);
        }
    }
    int failures = !identity + !creator_valid + !values_valid + !solves + !references_valid;
    pout << "matrix_identity_valid = " << (identity ? "true" : "false") << '\n'
         << "creator_lifetime_valid = " << (creator_valid ? "true" : "false") << '\n'
         << "retained_references_valid = " << (references_valid ? "true" : "false") << '\n'
         << "matrix_error = " << max_matrix_error << '\n'
         << "level_solve_error = " << max_level_solve_error << '\n'
         << "velocity_pressure_mapping_error = " << max_mapping_error << '\n'
         << "test_failures = " << failures << std::endl;
    return failures;
}

int
run_distributed_augmentation(Pointer<AppInitializer> app)
{
    LevelFixture fixture(app->getComponentDatabase("CartesianGeometry"), 0, true, true);
    const int rank = IBTK_MPI::getRank();
    Mat creator = nullptr, augmentation = nullptr;
    PetscErrorCode ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                                       fixture.dof_counts[rank],
                                       fixture.dof_counts[rank],
                                       PETSC_DECIDE,
                                       PETSC_DECIDE,
                                       1,
                                       nullptr,
                                       0,
                                       nullptr,
                                       &creator);
    IBTK_CHKERRQ(ierr);
    PetscInt first, last;
    ierr = MatGetOwnershipRange(creator, &first, &last);
    IBTK_CHKERRQ(ierr);
    for (PetscInt row = first; row < last; ++row)
    {
        ierr = MatSetValue(creator, row, row, 4.0, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(creator, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(creator, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    Mat expected = nullptr, original = nullptr;
    ierr = MatDuplicate(creator, MAT_COPY_VALUES, &expected);
    IBTK_CHKERRQ(ierr);
    ierr = MatDuplicate(creator, MAT_COPY_VALUES, &original);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(expected, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    IBTK_CHKERRQ(ierr);
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        fixture.velocity_counts[rank],
                        fixture.velocity_counts[rank],
                        PETSC_DECIDE,
                        PETSC_DECIDE,
                        1,
                        nullptr,
                        1,
                        nullptr,
                        &augmentation);
    IBTK_CHKERRQ(ierr);
    PetscInt vfirst, vlast;
    ierr = MatGetOwnershipRange(augmentation, &vfirst, &vlast);
    IBTK_CHKERRQ(ierr);
    int off_process = 0;
    for (PetscInt row = vfirst; row < vlast; ++row)
    {
        const PetscInt other = (row + fixture.velocity_ids.size() / 2) % fixture.velocity_ids.size();
        const PetscInt columns[2] = { row, other };
        const PetscScalar values[2] = { 1.0 + 0.001 * row, 0.125 };
        ierr = MatSetValues(augmentation, 1, &row, 2, columns, values, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
        const PetscInt full_row = fixture.velocity_ids[row];
        const PetscInt full_columns[2] = { full_row, fixture.velocity_ids[other] };
        ierr = MatSetValues(expected, 1, &full_row, 2, full_columns, values, ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        off_process += full_columns[1] < first || full_columns[1] >= last;
    }
    for (Mat mat : { augmentation, expected })
    {
        ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
    }
    Mat augmentation_original = nullptr;
    ierr = MatDuplicate(augmentation, MAT_COPY_VALUES, &augmentation_original);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScLevelSolver solver("distributed_augmentation", level_solver_database(), "");
    solver.setTimeInterval(0.0, 1.0);
    solver.setSolutionTime(1.0);
    solver.setOperatorMat(creator);
    solver.setAugmentedOperatorMat(augmentation);
    Mat retained_creator = creator, retained_augmentation = augmentation;
    ierr = MatDestroy(&creator);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&augmentation);
    IBTK_CHKERRQ(ierr);
    bool valid = IBTK_MPI::sumReduction(off_process) > 0;
    double action_error = 0.0;
    for (int cycle = 0; cycle < 2; ++cycle)
    {
        valid = matrix_references(retained_creator) == 1 && matrix_references(retained_augmentation) == 1 && valid;
        solver.initializeSolverState(*fixture.x, *fixture.b);
        Mat installed = nullptr;
        ierr = KSPGetOperators(solver.getPETScKSP(), &installed, nullptr);
        IBTK_CHKERRQ(ierr);
        valid = installed != retained_creator && matrices_equal(installed, expected) &&
                matrices_equal(retained_creator, original) &&
                matrices_equal(retained_augmentation, augmentation_original) && valid;
        Vec x = nullptr, result = nullptr, reference = nullptr;
        ierr = MatCreateVecs(expected, &x, &result);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(result, &reference);
        IBTK_CHKERRQ(ierr);
        PetscScalar* values;
        ierr = VecGetArray(x, &values);
        IBTK_CHKERRQ(ierr);
        for (PetscInt row = first; row < last; ++row)
        {
            values[row - first] = 1.0 + 0.01 * row;
        }
        ierr = VecRestoreArray(x, &values);
        IBTK_CHKERRQ(ierr);
        ierr = MatMult(installed, x, result);
        IBTK_CHKERRQ(ierr);
        ierr = MatMult(expected, x, reference);
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(result, -1.0, reference);
        IBTK_CHKERRQ(ierr);
        PetscReal error;
        ierr = VecNorm(result, NORM_INFINITY, &error);
        IBTK_CHKERRQ(ierr);
        valid = std::isfinite(error) && error <= 1.0e-12 && valid;
        action_error = std::max(action_error, error);
        for (Vec* vec : { &x, &result, &reference })
        {
            ierr = VecDestroy(vec);
            IBTK_CHKERRQ(ierr);
        }
        valid = check_level_solve(solver) && valid;
        solver.deallocateSolverState();
    }
    valid = matrix_references(retained_creator) == 1 && matrix_references(retained_augmentation) == 1 && valid;
    solver.setOperatorMat(nullptr);
    solver.setAugmentedOperatorMat(nullptr);
    for (Mat* mat : { &original, &expected, &augmentation_original })
    {
        ierr = MatDestroy(mat);
        IBTK_CHKERRQ(ierr);
    }
    pout << "distributed matrix error = " << max_matrix_error << '\n'
         << "distributed action error = " << action_error << '\n';
    if (!valid)
    {
        TBOX_ERROR("Distributed compact augmentation mapping, action, or creator lifetime failed.\n");
    }
    return 0;
}

int
run_initialized_matrix_setter(Pointer<AppInitializer> app, bool augmentation)
{
    // Use the repository's path-independent expected-error output.
    Pointer<Logger::Appender> abort_appender = new TestAppender();
    Logger::getInstance()->setAbortAppender(abort_appender);
    LevelFixture fixture(app->getComponentDatabase("CartesianGeometry"));
    StaggeredStokesPETScLevelSolver solver("initialized_setter", level_solver_database(), "");
    solver.setTimeInterval(0.0, 1.0);
    solver.setSolutionTime(1.0);
    Mat matrix = level_test_matrix(fixture.full_size, 4.0);
    solver.setOperatorMat(matrix);
    if (augmentation)
    {
        solver.setAugmentedOperatorMat(matrix);
    }
    solver.initializeSolverState(*fixture.x, *fixture.b);
    if (augmentation)
    {
        solver.setAugmentedOperatorMat(matrix);
    }
    else
    {
        solver.setOperatorMat(matrix);
    }
    // Even same-handle setters must reject initialized state. Returning zero
    // on unexpected continuation makes this expect_error=true case fail attest.
    solver.deallocateSolverState();
    PetscErrorCode ierr = MatDestroy(&matrix);
    IBTK_CHKERRQ(ierr);
    pout << "ERROR: initialized matrix setter unexpectedly returned.\n";
    return 0;
}

template <class Solver>
bool
check_shell_state(LevelSolverProbe<Solver>& solver,
                  Pointer<HierarchyVector> x,
                  Pointer<HierarchyVector> b,
                  PetscInt& overlap_total)
{
    bool valid = true;
    for (int cycle = 0; cycle < 2; ++cycle)
    {
        solver.initializeSolverState(*x, *b);
        std::vector<IS>*nonoverlap, *overlap;
        solver.getASMSubdomains(&nonoverlap, &overlap);
        valid = valid && overlap->size() == 8 && nonoverlap->size() == 8;
        overlap_total = 0;
        for (IS is : *overlap)
        {
            PetscInt n;
            PetscErrorCode ierr = ISGetSize(is, &n);
            IBTK_CHKERRQ(ierr);
            overlap_total += n;
        }
        valid = check_level_solve(solver) && valid;
        std::vector<Vec> retained = solver.retainShellVectors();
        valid = valid && retained.size() == 16;
        solver.deallocateSolverState();
        valid = solver.shellStorageEmpty() && valid;
        PetscInt max_references = 0;
        for (Vec& v : retained)
        {
            PetscInt references;
            PetscErrorCode ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(v), &references);
            IBTK_CHKERRQ(ierr);
            // Only this test reference may remain after solver teardown.
            valid = references == 1 && valid;
            max_references = std::max(max_references, references);
            ierr = VecDestroy(&v);
            IBTK_CHKERRQ(ierr);
        }
        if (max_references != 1)
        {
            pout << "unreleased_shell_vector_references = " << max_references << std::endl;
        }
    }
    return valid;
}

int
run_level_state(Pointer<AppInitializer> app)
{
    bool cc_valid = true, sc_valid = true, stokes_valid = true, domain_valid = true;
    {
        LevelFixture fixture(app->getComponentDatabase("CartesianGeometry"));
        Pointer<HierarchyVector> cc_x = new HierarchyVector("cc_x", fixture.hierarchy, 0, 0);
        Pointer<HierarchyVector> sc_x = new HierarchyVector("sc_x", fixture.hierarchy, 0, 0);
        cc_x->addComponent(fixture.x->getComponentVariable(1), fixture.x->getComponentDescriptorIndex(1));
        sc_x->addComponent(fixture.x->getComponentVariable(0), fixture.x->getComponentDescriptorIndex(0));
        Pointer<HierarchyVector> cc_b = cc_x->cloneVector("cc_b"), sc_b = sc_x->cloneVector("sc_b");
        cc_b->allocateVectorData();
        sc_b->allocateVectorData();
        PoissonSpecifications coefs("state_coefs");
        coefs.setCConstant(2.0);
        coefs.setDConstant(0.0);
        PetscInt cc_previous = 0, sc_previous = 0, stokes_previous = 0;
        for (int width : { 0, 2 })
        {
            Pointer<Database> db = level_solver_database("shell", width);
            // Exercise both existing shell compositions across these lifetimes.
            db->putString("shell_pc_type", width == 0 ? "multiplicative" : "additive");
            LevelSolverProbe<CCPoissonPETScLevelSolver> cc("state_cc", db);
            LevelSolverProbe<SCPoissonPETScLevelSolver> sc("state_sc", db);
            cc.setPoissonSpecifications(coefs);
            sc.setPoissonSpecifications(coefs);
            cc.setPhysicalBcCoef(nullptr);
            sc.setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr));
            cc.setTimeInterval(0.0, 1.0);
            sc.setTimeInterval(0.0, 1.0);
            PetscInt cc_total, sc_total, stokes_total;
            cc_valid = check_shell_state(cc, cc_x, cc_b, cc_total) && cc_valid;
            sc_valid = check_shell_state(sc, sc_x, sc_b, sc_total) && sc_valid;
            LevelSolverProbe<StaggeredStokesPETScLevelSolver> stokes("state_stokes", db);
            Mat creator = level_test_matrix(fixture.full_size, 4.0);
            stokes.setOperatorMat(creator);
            stokes.setTimeInterval(0.0, 1.0);
            stokes_valid = check_shell_state(stokes, fixture.x, fixture.b, stokes_total) && stokes_valid;
            stokes.setOperatorMat(nullptr);
            PetscErrorCode ierr = MatDestroy(&creator);
            IBTK_CHKERRQ(ierr);
            if (width == 2)
            {
                cc_valid = cc_valid && cc_total > cc_previous;
                sc_valid = sc_valid && sc_total > sc_previous;
                stokes_valid = stokes_valid && stokes_total > stokes_previous;
            }
            cc_previous = cc_total;
            sc_previous = sc_total;
            stokes_previous = stokes_total;
        }
        free_vector_components(*cc_b);
        free_vector_components(*sc_b);
    }
    for (int variant = 0; variant < 3; ++variant)
    {
        const bool full = variant != 2;
        LevelFixture fixture(app->getComponentDatabase("CartesianGeometry"), variant == 0 ? 0 : 1, full);
        LevelSolverProbe<StaggeredStokesPETScLevelSolver> solver("domain_stokes", level_solver_database());
        PoissonSpecifications coefs("domain_coefs");
        coefs.setCConstant(2.0);
        coefs.setDConstant(-0.01);
        solver.setVelocityPoissonSpecifications(coefs);
        solver.setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr), nullptr);
        solver.setComponentsHaveNullSpace(false, true);
        solver.setTimeInterval(0.0, 1.0);
        solver.setSolutionTime(1.0);
        solver.initializeSolverState(*fixture.x, *fixture.b);
        Mat matrix;
        MatNullSpace nullspace;
        PetscErrorCode ierr = KSPGetOperators(solver.getPETScKSP(), &matrix, nullptr);
        IBTK_CHKERRQ(ierr);
        ierr = MatGetNullSpace(matrix, &nullspace);
        IBTK_CHKERRQ(ierr);
        domain_valid = domain_valid && (full ? nullspace != nullptr : nullspace == nullptr);
        if (full && nullspace)
        {
            PetscBool is_nullspace;
            ierr = MatNullSpaceTest(nullspace, matrix, &is_nullspace);
            IBTK_CHKERRQ(ierr);
            domain_valid = domain_valid && is_nullspace;
        }
        solver.deallocateSolverState();
    }
    int failures = !cc_valid + !sc_valid + !stokes_valid + !domain_valid;
    pout << "shell_solve_error = " << max_level_solve_error << '\n'
         << "cc_state_valid = " << (cc_valid ? "true" : "false") << '\n'
         << "sc_state_valid = " << (sc_valid ? "true" : "false") << '\n'
         << "stokes_state_valid = " << (stokes_valid ? "true" : "false") << '\n'
         << "domain_nullspace_valid = " << (domain_valid ? "true" : "false") << '\n'
         << "test_failures = " << failures << std::endl;
    return failures;
}
} // namespace

namespace
{
#if defined(PETSC_USE_LOG)
PetscInt petsc_vec_creation_count = 0;
PetscInt petsc_vec_destruction_count = 0;

#if PETSC_VERSION_GE(3, 20, 0)
PetscErrorCode
ignore_petsc_log_event(PetscLogEvent, int, PetscObject, PetscObject, PetscObject, PetscObject)
{
    return 0;
}
#endif

PetscErrorCode
count_petsc_vec_creation(PetscObject object)
{
    PetscClassId class_id;
    const PetscErrorCode ierr = PetscObjectGetClassId(object, &class_id);
    if (ierr)
    {
        return ierr;
    }
    if (class_id == VEC_CLASSID)
    {
        ++petsc_vec_creation_count;
    }
    return 0;
}

PetscErrorCode
count_petsc_vec_destruction(PetscObject object)
{
    PetscClassId class_id;
    const PetscErrorCode ierr = PetscObjectGetClassId(object, &class_id);
    if (ierr)
    {
        return ierr;
    }
    if (class_id == VEC_CLASSID)
    {
        ++petsc_vec_destruction_count;
    }
    return 0;
}
#endif

struct StructureSpec
{
    int num_curve_points = 64;
    double ds = 1.0 / 64.0;
    double x_center = 0.5;
    double y_center = 0.5;
    double x_radius = 0.2;
    double y_radius = 0.2;
    double spring_stiffness = 2.0e2;
    int finest_ln = 0;
};

void
generate_structure(const unsigned int& strct_num,
                   const int& ln,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn,
                   void* ctx)
{
    auto* spec = static_cast<StructureSpec*>(ctx);
    if (!spec)
    {
        TBOX_ERROR("generate_structure(): missing structure specification context\n");
    }

    if (ln != spec->finest_ln || strct_num != 0)
    {
        num_vertices = 0;
        vertex_posn.resize(0);
        return;
    }

    num_vertices = spec->num_curve_points;
    vertex_posn.resize(num_vertices);
    for (int k = 0; k < num_vertices; ++k)
    {
        const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_vertices);
        vertex_posn[k](0) = spec->x_center + spec->x_radius * std::cos(theta);
        vertex_posn[k](1) = spec->y_center + spec->y_radius * std::sin(theta);
    }
}

void
generate_springs(
    const unsigned int& strct_num,
    const int& ln,
    std::multimap<int, IBRedundantInitializer::Edge>& spring_map,
    std::map<IBRedundantInitializer::Edge, IBRedundantInitializer::SpringSpec, IBRedundantInitializer::EdgeComp>&
        spring_spec,
    void* ctx)
{
    auto* spec = static_cast<StructureSpec*>(ctx);
    if (!spec)
    {
        TBOX_ERROR("generate_springs(): missing structure specification context\n");
    }
    if (ln != spec->finest_ln || strct_num != 0)
    {
        return;
    }

    for (int k = 0; k < spec->num_curve_points; ++k)
    {
        IBRedundantInitializer::Edge edge = { k, (k + 1) % spec->num_curve_points };
        if (edge.first > edge.second)
        {
            std::swap(edge.first, edge.second);
        }
        spring_map.insert(std::make_pair(edge.first, edge));

        IBRedundantInitializer::SpringSpec spec_data;
        spec_data.force_fcn_idx = 0;
        spec_data.parameters.resize(2);
        spec_data.parameters[0] = spec->spring_stiffness;
        spec_data.parameters[1] = 0.0;
        spring_spec.insert(std::make_pair(edge, spec_data));
    }
}

bool
side_l2_norm_is_finite(const Pointer<HierarchySideDataOpsReal<NDIM, double>>& side_data_ops,
                       const int data_idx,
                       const int weight_idx,
                       double& l2_norm)
{
    l2_norm = side_data_ops->L2Norm(data_idx, weight_idx);
    return std::isfinite(l2_norm);
}

bool
cell_l2_norm_is_finite(const Pointer<HierarchyCellDataOpsReal<NDIM, double>>& cell_data_ops,
                       const int data_idx,
                       const int weight_idx,
                       double& l2_norm)
{
    l2_norm = cell_data_ops->L2Norm(data_idx, weight_idx);
    return std::isfinite(l2_norm);
}

void
set_divergence_free_probe_velocity(const int u_idx, Pointer<PatchHierarchy<NDIM>> patch_hierarchy)
{
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();
            const double* const x_lower = patch_geom->getXLower();
            const SAMRAI::hier::Index<NDIM>& lower = patch_box.lower();

            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                for (Box<NDIM>::Iterator b(side_box); b; b++)
                {
                    const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                    double x[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        const double offset = (d == axis) ? 0.0 : 0.5;
                        x[d] = x_lower[d] + (static_cast<double>(i_s(d) - lower(d)) + offset) * dx[d];
                    }
                    const double val = (axis == 0) ? std::sin(2.0 * M_PI * x[1]) : -std::sin(2.0 * M_PI * x[0]);
                    (*u_data)(i_s) = val;
                }
            }
        }
    }
}

int
run_foundation(Pointer<AppInitializer> app_initializer)
{
    int test_failures = 0;
    {
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        const double current_time = 0.0;
        const double dt = input_db->getDoubleWithDefault("DT", 0.005);
        const double rho = input_db->getDoubleWithDefault("RHO", 1.0);
        const double mu = input_db->getDoubleWithDefault("MU", 1.0);
        const double new_time = current_time + dt;
        const bool use_fixed_le_operators = input_db->getBoolWithDefault("USE_FIXED_LE_OPERATORS", true);
        StructureSpec structure_spec;
        structure_spec.ds = input_db->getDoubleWithDefault("DS", 1.0 / 64.0);
        structure_spec.x_center = input_db->getDoubleWithDefault("X_CENTER", 0.5);
        structure_spec.y_center = input_db->getDoubleWithDefault("Y_CENTER", 0.5);
        structure_spec.x_radius = input_db->getDoubleWithDefault("X_RADIUS", 0.2);
        structure_spec.y_radius = input_db->getDoubleWithDefault("Y_RADIUS", 0.2);
        structure_spec.spring_stiffness = input_db->getDoubleWithDefault("SPRING_STIFFNESS", 2.0e2);
        if (!(structure_spec.ds > 0.0))
        {
            TBOX_ERROR("DS must be positive\n");
        }
        if (!(structure_spec.x_radius > 0.0) || !(structure_spec.y_radius > 0.0))
        {
            TBOX_ERROR("X_RADIUS and Y_RADIUS must be positive\n");
        }
        const double a = structure_spec.x_radius;
        const double b = structure_spec.y_radius;
        const double h = std::pow((a - b) / (a + b), 2.0);
        const double circumference =
            M_PI * (a + b) * (1.0 + 3.0 * h / (10.0 + std::sqrt(std::max(0.0, 4.0 - 3.0 * h))));
        structure_spec.num_curve_points = std::max(3, static_cast<int>(circumference / structure_spec.ds));
        if (structure_spec.num_curve_points < 3)
        {
            TBOX_ERROR("computed num_curve_points must be >= 3\n");
        }

        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        ib_method_ops->setUseFixedLEOperators(true);

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               ib_method_ops,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        structure_spec.finest_ln = input_db->getIntegerWithDefault("MAX_LEVELS", 1) - 1;
        ib_initializer->setStructureNamesOnLevel(structure_spec.finest_ln, { "curve2d" });
        ib_initializer->registerInitStructureFunction(generate_structure, &structure_spec);
        ib_initializer->registerInitSpringDataFunction(generate_springs, &structure_spec);
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, current_time);
        int tag_buffer = input_db->getIntegerWithDefault("TAG_BUFFER", 1);
        int level_number = 0;
        bool done = false;
        while (!done && gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, current_time, true, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> current_ctx = var_db->getContext("current_ctx");
        Pointer<VariableContext> scratch_ctx = var_db->getContext("scratch_ctx");
        Pointer<VariableContext> solver_ctx = var_db->getContext("solver_ctx");

        Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f_var");
        Pointer<CellVariable<NDIM, double>> g_var = new CellVariable<NDIM, double>("g_var");
        Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p_var");
        Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("p_dof_index");
        Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u_var");
        Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("u_dof_index");

        const IntVector<NDIM> ib_ghosts = ib_method_ops->getMinimumGhostCellWidth();
        const IntVector<NDIM> one_ghost = IntVector<NDIM>(1);

        const int u_current_idx = var_db->registerVariableAndContext(u_var, current_ctx, ib_ghosts);
        const int u_sol_idx = var_db->registerVariableAndContext(u_var, solver_ctx, one_ghost);
        const int f_rhs_idx = var_db->registerVariableAndContext(f_var, solver_ctx, one_ghost);
        const int p_sol_idx = var_db->registerVariableAndContext(p_var, solver_ctx, one_ghost);
        const int g_rhs_idx = var_db->registerVariableAndContext(g_var, solver_ctx, one_ghost);
        const int u_scratch_idx = var_db->registerVariableAndContext(u_var, scratch_ctx, ib_ghosts);
        const int f_scratch_idx = var_db->registerVariableAndContext(f_var, scratch_ctx, ib_ghosts);
        const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, scratch_ctx, ib_ghosts);
        // The independent MAC matrix assembly accesses neighboring pressure indices.
        const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, scratch_ctx, one_ghost);

        const std::vector<int> allocated_patch_data_indices = {
            u_current_idx, u_scratch_idx, f_scratch_idx, u_dof_index_idx, p_dof_index_idx
        };
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int data_idx : allocated_patch_data_indices)
            {
                level->allocatePatchData(data_idx, current_time);
            }
        }

        Pointer<HierarchySideDataOpsReal<NDIM, double>> hier_velocity_data_ops =
            new HierarchySideDataOpsReal<NDIM, double>(patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        Pointer<HierarchyCellDataOpsReal<NDIM, double>> hier_pressure_data_ops =
            new HierarchyCellDataOpsReal<NDIM, double>(patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            u_init->setDataOnPatchHierarchy(u_current_idx, u_var, patch_hierarchy, current_time);
        }
        else
        {
            hier_velocity_data_ops->setToScalar(u_current_idx, 0.0, false);
        }
        hier_velocity_data_ops->setToScalar(u_scratch_idx, 0.0, false);
        hier_velocity_data_ops->setToScalar(f_scratch_idx, 0.0, false);

        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        std::vector<Pointer<CoarsenSchedule<NDIM>>> u_synch_scheds(finest_ln + 1);
        std::vector<Pointer<RefineSchedule<NDIM>>> u_ghost_fill_scheds(finest_ln + 1);
        std::vector<Pointer<RefineSchedule<NDIM>>> f_prolongation_scheds(finest_ln + 1);
        RefineAlgorithm<NDIM> velocity_ghost_fill;
        velocity_ghost_fill.registerRefine(u_scratch_idx, u_scratch_idx, u_scratch_idx, nullptr);
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            u_ghost_fill_scheds[ln] = velocity_ghost_fill.createSchedule(patch_hierarchy->getPatchLevel(ln));
        }

        // Populate the Lagrangian ghost-node/periodic-image distribution before spreading.
        ib_method_ops->beginDataRedistribution(patch_hierarchy, gridding_algorithm);
        ib_method_ops->endDataRedistribution(patch_hierarchy, gridding_algorithm);
        ib_method_ops->initializePatchHierarchy(patch_hierarchy,
                                                gridding_algorithm,
                                                u_current_idx,
                                                u_synch_scheds,
                                                u_ghost_fill_scheds,
                                                0,
                                                current_time,
                                                true);
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();

        ib_method_ops->preprocessIntegrateData(current_time, new_time, /*num_cycles*/ 1);
        ib_method_ops->updateFixedLEOperators();

        std::vector<std::vector<int>> num_dofs_per_proc(finest_ln + 1);
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
                num_dofs_per_proc[ln], u_dof_index_idx, p_dof_index_idx, level);
        }

        Mat A = nullptr;
        ib_method_ops->constructLagrangianForceJacobian(A, MATAIJ, new_time);
        Mat J = nullptr;
        // Preprocessing starts with X_new = X_current, before any position update.
        const Vec initial_position = ib_method_ops->getLDataManager()->getLData("X", finest_ln)->getVec();
        PETScMatUtilities::constructPatchLevelSCInterpOp(J,
                                                         IBKernelEvaluatorTensorProduct{ IBKernels::IB4{} },
                                                         initial_position,
                                                         num_dofs_per_proc[finest_ln],
                                                         u_dof_index_idx,
                                                         patch_hierarchy->getPatchLevel(finest_ln));

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        Pointer<SAMRAIVectorReal<NDIM, double>> eul_sol_vec =
            new SAMRAIVectorReal<NDIM, double>("eul_sol_vec", patch_hierarchy, 0, finest_ln);
        eul_sol_vec->addComponent(u_var, u_sol_idx, wgt_sc_idx, hier_velocity_data_ops);
        eul_sol_vec->addComponent(p_var, p_sol_idx, wgt_cc_idx, hier_pressure_data_ops);
        eul_sol_vec->allocateVectorData();

        Pointer<SAMRAIVectorReal<NDIM, double>> eul_rhs_vec =
            new SAMRAIVectorReal<NDIM, double>("eul_rhs_vec", patch_hierarchy, 0, finest_ln);
        eul_rhs_vec->addComponent(f_var, f_rhs_idx, wgt_sc_idx, hier_velocity_data_ops);
        eul_rhs_vec->addComponent(g_var, g_rhs_idx, wgt_cc_idx, hier_pressure_data_ops);
        eul_rhs_vec->allocateVectorData();

        hier_velocity_data_ops->copyData(u_sol_idx, u_current_idx);
        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            p_init->setDataOnPatchHierarchy(p_sol_idx, p_var, patch_hierarchy, current_time);
        }
        else
        {
            hier_pressure_data_ops->setToScalar(p_sol_idx, 0.0, false);
        }
        hier_velocity_data_ops->setToScalar(f_rhs_idx, 0.0, false);
        hier_pressure_data_ops->setToScalar(g_rhs_idx, 0.0, false);

        const double lambda = 0.0;
        PoissonSpecifications U_problem_coefs("stokes_ib_solver_components::U_problem_coefs");
        U_problem_coefs.setCConstant(rho / dt + lambda);
        U_problem_coefs.setDConstant(-mu);

        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        if (periodic_shift.min() <= 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);
                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);
                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
        }

        Pointer<StaggeredStokesOperator> stokes_op =
            new StaggeredStokesOperator("stokes_ib_solver_components::stokes_op", false);
        stokes_op->setVelocityPoissonSpecifications(U_problem_coefs);
        stokes_op->setPhysicalBcCoefs(u_bc_coefs, nullptr);
        stokes_op->setTimeInterval(current_time, new_time);
        stokes_op->setSolutionTime(new_time);

        StaggeredStokesIBOperator::Context ctx;
        ctx.ib_implicit_ops = ib_method_ops;
        ctx.stokes_op = stokes_op;
        ctx.u_phys_bdry_op = nullptr;
        ctx.hier_velocity_data_ops = hier_velocity_data_ops;
        ctx.u_synch_scheds = u_synch_scheds;
        ctx.u_ghost_fill_scheds = u_ghost_fill_scheds;
        ctx.f_prolongation_scheds = f_prolongation_scheds;
        ctx.u_idx = u_scratch_idx;
        ctx.f_idx = f_scratch_idx;
        ctx.u_current_idx = u_current_idx;
        ctx.u_dof_index_idx = u_dof_index_idx;
        ctx.p_dof_index_idx = p_dof_index_idx;
        ctx.use_fixed_le_operators = use_fixed_le_operators;
        ctx.time_stepping_type = IBAMR::string_to_enum<TimeSteppingType>(
            input_db->getStringWithDefault("IB_TIME_STEPPING", "MIDPOINT_RULE"));

        StaggeredStokesIBOperator nonlinear_op("stokes_ib_solver_components::nonlinear_op", false);
        nonlinear_op.setOperatorContext(ctx);
        nonlinear_op.setTimeInterval(current_time, new_time);
        nonlinear_op.setSolutionTime(new_time);
        nonlinear_op.initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);

        Pointer<SAMRAIVectorReal<NDIM, double>> nonlinear_probe = eul_sol_vec->cloneVector("nonlinear_probe");
        nonlinear_probe->allocateVectorData();
        nonlinear_probe->setToScalar(0.0);
        set_divergence_free_probe_velocity(nonlinear_probe->getComponentDescriptorIndex(0), patch_hierarchy);

        Pointer<SAMRAIVectorReal<NDIM, double>> f_probe = eul_rhs_vec->cloneVector("f_probe");
        f_probe->allocateVectorData();
        f_probe->setToScalar(0.0);
        nonlinear_op.apply(*nonlinear_probe, *f_probe);

        double nonlinear_side_norm = std::numeric_limits<double>::quiet_NaN();
        double nonlinear_cell_norm = std::numeric_limits<double>::quiet_NaN();
        const bool expect_trivial_nonlinear = (std::abs(structure_spec.spring_stiffness) <= 1.0e-14) &&
                                              (std::abs(rho) <= 1.0e-14) && (std::abs(mu) <= 1.0e-14);
        if (!side_l2_norm_is_finite(
                hier_velocity_data_ops, f_probe->getComponentDescriptorIndex(0), wgt_sc_idx, nonlinear_side_norm) ||
            !cell_l2_norm_is_finite(
                hier_pressure_data_ops, f_probe->getComponentDescriptorIndex(1), wgt_cc_idx, nonlinear_cell_norm))
        {
            ++test_failures;
            pout << "nonlinear operator produced non-finite norm" << std::endl;
        }
        else if (nonlinear_side_norm <= 1.0e-14 && nonlinear_cell_norm <= 1.0e-14)
        {
            if (!expect_trivial_nonlinear)
            {
                ++test_failures;
                pout << "nonlinear operator action is trivial" << std::endl;
            }
        }
        else if (expect_trivial_nonlinear)
        {
            ++test_failures;
            pout << "nonlinear operator action is nontrivial when SPRING_STIFFNESS, RHO, and MU are zero" << std::endl;
        }

        Pointer<StaggeredStokesIBJacobianOperator> jac_op =
            new StaggeredStokesIBJacobianOperator("stokes_ib_solver_components::jacobian_op");
        jac_op->setOperatorContext(ctx);
        jac_op->setTimeInterval(current_time, new_time);
        jac_op->setSolutionTime(new_time);
        // Restore the physical base before initialization freezes the coupling
        // positions; the preceding nonlinear probe left a different endpoint.
        nonlinear_op.apply(*eul_sol_vec, *f_probe);
        jac_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);
        jac_op->formJacobian(*eul_sol_vec);

        Pointer<SAMRAIVectorReal<NDIM, double>> v = eul_sol_vec->cloneVector("v");
        v->allocateVectorData();
        v->setToScalar(0.0);
        set_divergence_free_probe_velocity(v->getComponentDescriptorIndex(0), patch_hierarchy);
        hier_pressure_data_ops->setToScalar(v->getComponentDescriptorIndex(1), -0.25, false);

        Pointer<SAMRAIVectorReal<NDIM, double>> jv = eul_rhs_vec->cloneVector("jv");
        jv->allocateVectorData();
        jv->setToScalar(0.0);
        jac_op->apply(*v, *jv);

        Pointer<SAMRAIVectorReal<NDIM, double>> diff = eul_rhs_vec->cloneVector("diff");
        diff->allocateVectorData();
        diff->setToScalar(0.0);

        double jv_side_norm = std::numeric_limits<double>::quiet_NaN();
        double jv_cell_norm = std::numeric_limits<double>::quiet_NaN();
        const bool jv_finite =
            side_l2_norm_is_finite(
                hier_velocity_data_ops, jv->getComponentDescriptorIndex(0), wgt_sc_idx, jv_side_norm) &&
            cell_l2_norm_is_finite(
                hier_pressure_data_ops, jv->getComponentDescriptorIndex(1), wgt_cc_idx, jv_cell_norm);
        if (!jv_finite)
        {
            ++test_failures;
            pout << "jacobian norms are non-finite" << std::endl;
        }
        else if (jv_side_norm <= 1.0e-14 && jv_cell_norm <= 1.0e-14)
        {
            ++test_failures;
            pout << "jacobian action is trivial" << std::endl;
        }
        nonlinear_op.apply(*eul_sol_vec, *f_probe);

        Pointer<Database> stokes_ib_precond_db =
            input_db->isDatabase("stokes_ib_precond_db") ? input_db->getDatabase("stokes_ib_precond_db") : nullptr;

        Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_op = new StaggeredStokesIBLevelRelaxationFACOperator(
            "stokes_ib_solver_components::fac_op", stokes_ib_precond_db, "stokes_ib_pc_");
        Pointer<StaggeredStokesIBJacobianFACPreconditioner> fac_pc = new StaggeredStokesIBJacobianFACPreconditioner(
            "stokes_ib_solver_components::fac_pc", fac_op, stokes_ib_precond_db, "stokes_ib_pc_");
        Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();

        fac_pc->setVelocityPoissonSpecifications(U_problem_coefs);
        fac_pc->setPhysicalBcCoefs(u_bc_coefs, nullptr);
        fac_pc->setPhysicalBoundaryHelper(bc_helper);
        fac_pc->setTimeInterval(current_time, new_time);
        fac_pc->setSolutionTime(new_time);
        fac_pc->setHomogeneousBc(true);
        fac_pc->setComponentsHaveNullSpace(false, true);
        fac_pc->setIBTimeSteppingType(ctx.time_stepping_type);
        fac_pc->setIBForceJacobian(A);
        fac_pc->setIBInterpOp(J);
        fac_pc->setIBImplicitStrategy(ib_method_ops);
        fac_pc->initializeSolverState(*eul_sol_vec, *eul_rhs_vec);

        const bool verify_galerkin_operator_borrowing =
            input_db->getBoolWithDefault("VERIFY_GALERKIN_OPERATOR_BORROWING", false);
        const bool rediscretize_stokes = stokes_ib_precond_db->getBoolWithDefault("rediscretize_stokes", true);
        const bool rediscretize_residual = stokes_ib_precond_db->getBoolWithDefault("res_rediscretized_stokes", true);
        std::vector<Mat> coupling_reference(finest_ln + 1, nullptr), rediscretized_reference(finest_ln + 1, nullptr),
            operator_reference(finest_ln + 1, nullptr);
        for (int ln = finest_ln; ln >= 0; --ln)
        {
            PetscErrorCode ierr;
            if (ln == finest_ln)
            {
                ierr = MatPtAP(A, J, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &coupling_reference[ln]);
                IBTK_CHKERRQ(ierr);
                const double beta = ctx.time_stepping_type == BACKWARD_EULER ? 1.0 : 0.5;
                double scale = -dt * beta * beta;
                const IntVector<NDIM> ratio = patch_hierarchy->getPatchLevel(ln)->getRatio();
                for (int d = 0; d < NDIM; ++d)
                {
                    scale *= ratio(d) / grid_geometry->getDx()[d];
                }
                ierr = MatScale(coupling_reference[ln], scale);
                IBTK_CHKERRQ(ierr);
            }
            else
            {
                ierr = MatPtAP(coupling_reference[ln + 1],
                               fac_op->getProlongationOp(ln),
                               MAT_INITIAL_MATRIX,
                               PETSC_DEFAULT,
                               &coupling_reference[ln]);
                IBTK_CHKERRQ(ierr);
                ierr = MatDiagonalScale(coupling_reference[ln], fac_op->getRestrictionScalingOp(ln), nullptr);
                IBTK_CHKERRQ(ierr);
            }
            StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(rediscretized_reference[ln],
                                                                             U_problem_coefs,
                                                                             u_bc_coefs,
                                                                             new_time,
                                                                             num_dofs_per_proc[ln],
                                                                             u_dof_index_idx,
                                                                             p_dof_index_idx,
                                                                             patch_hierarchy->getPatchLevel(ln));
            ierr = MatAXPY(rediscretized_reference[ln], 1.0, coupling_reference[ln], DIFFERENT_NONZERO_PATTERN);
            IBTK_CHKERRQ(ierr);
            if (rediscretize_stokes || ln == finest_ln)
            {
                ierr = MatDuplicate(rediscretized_reference[ln], MAT_COPY_VALUES, &operator_reference[ln]);
                IBTK_CHKERRQ(ierr);
            }
            else
            {
                AO ordering = nullptr;
                int u_offset = 0, p_offset = 0;
                StaggeredStokesPETScVecUtilities::constructPatchLevelAO(ordering,
                                                                        num_dofs_per_proc[ln],
                                                                        u_dof_index_idx,
                                                                        p_dof_index_idx,
                                                                        patch_hierarchy->getPatchLevel(ln),
                                                                        u_offset,
                                                                        p_offset);
                Mat prolongation = nullptr;
                Vec scaling = nullptr;
                // These fixtures use the standard RT0 velocity and conservative pressure transfers.
                StaggeredStokesPETScMatUtilities::constructProlongationOp(prolongation,
                                                                          "RT0",
                                                                          "CONSERVATIVE",
                                                                          u_dof_index_idx,
                                                                          p_dof_index_idx,
                                                                          num_dofs_per_proc[ln + 1],
                                                                          num_dofs_per_proc[ln],
                                                                          patch_hierarchy->getPatchLevel(ln + 1),
                                                                          patch_hierarchy->getPatchLevel(ln),
                                                                          ordering,
                                                                          u_offset,
                                                                          p_offset);
                PETScMatUtilities::constructRestrictionScalingOp(prolongation, scaling);
                ierr = MatPtAP(operator_reference[ln + 1],
                               prolongation,
                               MAT_INITIAL_MATRIX,
                               PETSC_DEFAULT,
                               &operator_reference[ln]);
                IBTK_CHKERRQ(ierr);
                ierr = MatDiagonalScale(operator_reference[ln], scaling, nullptr);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&scaling);
                IBTK_CHKERRQ(ierr);
                ierr = MatDestroy(&prolongation);
                IBTK_CHKERRQ(ierr);
                ierr = AODestroy(&ordering);
                IBTK_CHKERRQ(ierr);
            }
        }
        constexpr double FAC_ACTION_TOL = 1.0e-9;
        pout << "FAC checks use stated accuracy bounds; attained roundoff may vary.\n";
        auto check_fac_residual_work_vector_cache = [&](double& reuse_error)
        {
            Pointer<SAMRAIVectorReal<NDIM, double>> first_residual =
                eul_rhs_vec->cloneVector("fac_residual_workspace_first");
            Pointer<SAMRAIVectorReal<NDIM, double>> second_residual =
                eul_rhs_vec->cloneVector("fac_residual_workspace_second");
            first_residual->allocateVectorData();
            second_residual->allocateVectorData();

#if defined(PETSC_USE_LOG)
            petsc_vec_creation_count = 0;
            petsc_vec_destruction_count = 0;
#if PETSC_VERSION_GE(3, 20, 0)
            PetscLogHandler log_handler = nullptr;
            PetscErrorCode log_ierr = PetscLogHandlerCreateLegacy(PETSC_COMM_WORLD,
                                                                  ignore_petsc_log_event,
                                                                  ignore_petsc_log_event,
                                                                  count_petsc_vec_creation,
                                                                  count_petsc_vec_destruction,
                                                                  &log_handler);
            IBTK_CHKERRQ(log_ierr);
            log_ierr = PetscLogHandlerStart(log_handler);
            IBTK_CHKERRQ(log_ierr);
#else
            // Older PETSc versions expose the object callbacks directly.
            const auto saved_create = PetscLogPHC;
            const auto saved_destroy = PetscLogPHD;
            PetscLogPHC = count_petsc_vec_creation;
            PetscLogPHD = count_petsc_vec_destruction;
#endif
            // Verify that logging is active before trusting zero event counts.
            Vec logging_probe = nullptr;
            PetscErrorCode ierr = VecCreate(PETSC_COMM_WORLD, &logging_probe);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&logging_probe);
            IBTK_CHKERRQ(ierr);
            const bool logging_active = petsc_vec_creation_count == 1 && petsc_vec_destruction_count == 1;
            if (!logging_active)
            {
                pout << "PETSc Vec logging calibration failed" << std::endl;
            }
            petsc_vec_creation_count = 0;
            petsc_vec_destruction_count = 0;
#endif
            fac_op->computeResidual(*first_residual, *nonlinear_probe, *eul_rhs_vec, 0, finest_ln);
            fac_op->computeResidual(*second_residual, *nonlinear_probe, *eul_rhs_vec, 0, finest_ln);
#if defined(PETSC_USE_LOG)
#if PETSC_VERSION_GE(3, 20, 0)
            ierr = PetscLogHandlerStop(log_handler);
            IBTK_CHKERRQ(ierr);
            ierr = PetscLogHandlerDestroy(&log_handler);
            IBTK_CHKERRQ(ierr);
#else
            PetscLogPHC = saved_create;
            PetscLogPHD = saved_destroy;
#endif
            if (!logging_active || petsc_vec_creation_count != 0 || petsc_vec_destruction_count != 0)
            {
                ++test_failures;
                pout << "FAC residual Vec allocation check failed: creations = " << petsc_vec_creation_count
                     << ", destructions = " << petsc_vec_destruction_count << std::endl;
            }
#endif
            Pointer<SAMRAIVectorReal<NDIM, double>> stokes_reference =
                eul_rhs_vec->cloneVector("fac_residual_stokes_reference");
            stokes_reference->allocateVectorData();
            if (rediscretize_residual)
            {
                stokes_op->apply(*nonlinear_probe, *stokes_reference);
            }
            double operator_error = 0.0, composition_error = 0.0, elastic_action_norm = 0.0;
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                Mat installed = nullptr, matrix_difference = nullptr;
                PetscErrorCode check_ierr =
                    KSPGetOperators(fac_op->getStaggeredStokesPETScLevelSolver(ln)->getPETScKSP(), &installed, nullptr);
                IBTK_CHKERRQ(check_ierr);
                check_ierr = MatDuplicate(installed, MAT_COPY_VALUES, &matrix_difference);
                IBTK_CHKERRQ(check_ierr);
                check_ierr = MatAXPY(matrix_difference, -1.0, operator_reference[ln], DIFFERENT_NONZERO_PATTERN);
                IBTK_CHKERRQ(check_ierr);
                PetscReal error = 0.0;
                check_ierr = MatNorm(matrix_difference, NORM_INFINITY, &error);
                IBTK_CHKERRQ(check_ierr);
                operator_error = std::max(operator_error, error);
                check_ierr = MatDestroy(&matrix_difference);
                IBTK_CHKERRQ(check_ierr);
                Vec solution = nullptr, expected = nullptr, actual = nullptr;
                check_ierr = MatCreateVecs(operator_reference[ln], &solution, &expected);
                IBTK_CHKERRQ(check_ierr);
                check_ierr = VecDuplicate(expected, &actual);
                IBTK_CHKERRQ(check_ierr);
                const Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
                StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(solution,
                                                                      nonlinear_probe->getComponentDescriptorIndex(0),
                                                                      u_dof_index_idx,
                                                                      nonlinear_probe->getComponentDescriptorIndex(1),
                                                                      p_dof_index_idx,
                                                                      level);
                check_ierr = MatMult(coupling_reference[ln], solution, expected);
                IBTK_CHKERRQ(check_ierr);
                check_ierr = VecNorm(expected, NORM_2, &error);
                IBTK_CHKERRQ(check_ierr);
                elastic_action_norm = std::max(elastic_action_norm, error);
                if (rediscretize_residual)
                {
                    // Use the hierarchy Stokes action to retain coarse-fine synchronization.
                    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
                        actual,
                        stokes_reference->getComponentDescriptorIndex(0),
                        u_dof_index_idx,
                        stokes_reference->getComponentDescriptorIndex(1),
                        p_dof_index_idx,
                        level);
                    check_ierr = VecAXPY(expected, 1.0, actual);
                    IBTK_CHKERRQ(check_ierr);
                }
                else
                {
                    check_ierr = MatMult(operator_reference[ln], solution, expected);
                    IBTK_CHKERRQ(check_ierr);
                }
                StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(actual,
                                                                      eul_rhs_vec->getComponentDescriptorIndex(0),
                                                                      u_dof_index_idx,
                                                                      eul_rhs_vec->getComponentDescriptorIndex(1),
                                                                      p_dof_index_idx,
                                                                      level);
                check_ierr = VecAYPX(expected, -1.0, actual);
                IBTK_CHKERRQ(check_ierr);
                StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(actual,
                                                                      first_residual->getComponentDescriptorIndex(0),
                                                                      u_dof_index_idx,
                                                                      first_residual->getComponentDescriptorIndex(1),
                                                                      p_dof_index_idx,
                                                                      level);
                check_ierr = VecAXPY(actual, -1.0, expected);
                IBTK_CHKERRQ(check_ierr);
                check_ierr = VecNorm(actual, NORM_INFINITY, &error);
                IBTK_CHKERRQ(check_ierr);
                composition_error = std::max(composition_error, error);
                check_ierr = VecDestroy(&solution);
                IBTK_CHKERRQ(check_ierr);
                check_ierr = VecDestroy(&expected);
                IBTK_CHKERRQ(check_ierr);
                check_ierr = VecDestroy(&actual);
                IBTK_CHKERRQ(check_ierr);
            }
            const bool operator_valid = std::isfinite(operator_error) && operator_error <= FAC_ACTION_TOL;
            const bool composition_valid = std::isfinite(composition_error) && composition_error <= FAC_ACTION_TOL;
            pout << "fac_operator_valid = " << (operator_valid ? "true" : "false")
                 << ", residual_composition_valid = " << (composition_valid ? "true" : "false")
                 << ", absolute_inf_bound = " << FAC_ACTION_TOL << ", elastic_action_norm = " << elastic_action_norm
                 << std::endl;
            if (!operator_valid || !composition_valid || !std::isfinite(elastic_action_norm) ||
                elastic_action_norm <= 1.0e-12)
            {
                pout << "FAC comparison failed: operator_error = " << std::setprecision(17) << operator_error
                     << ", composition_error = " << composition_error
                     << ", elastic_action_norm = " << elastic_action_norm << std::setprecision(6) << std::endl;
                ++test_failures;
            }
            const double residual_norm = first_residual->maxNorm();
            first_residual->subtract(first_residual, second_residual);
            reuse_error = std::abs(first_residual->maxNorm());
            free_vector_components(*first_residual);
            free_vector_components(*second_residual);
            free_vector_components(*stokes_reference);
            return std::isfinite(residual_norm) && residual_norm > 1.0e-12 && std::isfinite(reuse_error) &&
                   reuse_error == 0.0;
        };

        double fac_residual_work_vector_reuse_error = 0.0;
        const bool fac_residual_repeat_valid =
            check_fac_residual_work_vector_cache(fac_residual_work_vector_reuse_error);
        if (!fac_residual_repeat_valid)
        {
            ++test_failures;
        }

        bool galerkin_operator_available_valid = !verify_galerkin_operator_borrowing;
        bool galerkin_operator_creator_lifetime_valid = !verify_galerkin_operator_borrowing;
        if (verify_galerkin_operator_borrowing)
        {
            Pointer<StaggeredStokesPETScLevelSolver> coarse_level_solver =
                fac_op->getStaggeredStokesPETScLevelSolver(0);
            Mat supplied_operator = nullptr;
            PetscErrorCode ierr = KSPGetOperators(coarse_level_solver->getPETScKSP(), &supplied_operator, nullptr);
            IBTK_CHKERRQ(ierr);
            PetscInt reference_count_before = 0;
            ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(supplied_operator), &reference_count_before);
            IBTK_CHKERRQ(ierr);
            galerkin_operator_available_valid = supplied_operator != nullptr && reference_count_before > 1;
            coarse_level_solver->deallocateSolverState();
            PetscInt n_rows = 0;
            PetscInt n_columns = 0;
            ierr = MatGetSize(supplied_operator, &n_rows, &n_columns);
            IBTK_CHKERRQ(ierr);
            PetscInt reference_count_after = 0;
            ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(supplied_operator), &reference_count_after);
            IBTK_CHKERRQ(ierr);
            // Solver-state teardown retains the installed input in addition to
            // the FAC creator's reference. Explicit clearing releases only the former.
            galerkin_operator_creator_lifetime_valid = n_rows > 0 && n_rows == n_columns && reference_count_after == 2;
            coarse_level_solver->setOperatorMat(nullptr);
            galerkin_operator_creator_lifetime_valid =
                matrix_references(supplied_operator) == 1 && galerkin_operator_creator_lifetime_valid;
            if (!galerkin_operator_available_valid || !galerkin_operator_creator_lifetime_valid)
            {
                ++test_failures;
            }
        }

        bool fac_reinitialization_valid = true;
        double fac_residual_work_vector_reinitialize_error = 0.0;
        fac_pc->deallocateSolverState();
        fac_pc->initializeSolverState(*eul_sol_vec, *eul_rhs_vec);
        const bool fac_residual_repeat_reinitialize_valid =
            check_fac_residual_work_vector_cache(fac_residual_work_vector_reinitialize_error);
        if (verify_galerkin_operator_borrowing)
        {
            Pointer<StaggeredStokesPETScLevelSolver> coarse_level_solver =
                fac_op->getStaggeredStokesPETScLevelSolver(0);
            Mat supplied_operator = nullptr;
            PetscErrorCode ierr = KSPGetOperators(coarse_level_solver->getPETScKSP(), &supplied_operator, nullptr);
            IBTK_CHKERRQ(ierr);
            PetscInt reference_count = 0;
            ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(supplied_operator), &reference_count);
            IBTK_CHKERRQ(ierr);
            fac_reinitialization_valid = supplied_operator != nullptr && reference_count > 1;
        }
        if (!fac_residual_repeat_reinitialize_valid || !fac_reinitialization_valid)
        {
            ++test_failures;
        }

        Pointer<PETScKrylovLinearSolver> linear_solver =
            new PETScKrylovLinearSolver("stokes_ib_solver_components::linear_solver", nullptr, "ib_");
        linear_solver->setOperator(jac_op);
        linear_solver->setPreconditioner(fac_pc);
        linear_solver->setTimeInterval(current_time, new_time);
        linear_solver->setSolutionTime(new_time);
        linear_solver->setInitialGuessNonzero(false);
        // Richardson self-scaling makes the FAC action nonlinear. Flexible GMRES
        // uses a right preconditioner and measures the unpreconditioned residual.
        linear_solver->setKSPType("fgmres");

        Pointer<SAMRAIVectorReal<NDIM, double>> linear_sol = eul_sol_vec->cloneVector("linear_sol");
        linear_sol->allocateVectorData();
        linear_sol->setToScalar(0.0);
        linear_solver->initializeSolverState(*linear_sol, *jv);
        // Outer initialization clears the Jacobian base and rebuilds FAC state.
        jac_op->formJacobian(*eul_sol_vec);
        if (finest_ln == 0)
        {
            jac_op->setIBCouplingJacobian(fac_op->getEulerianElasticityLevelOp(finest_ln));
        }
        // Multilevel application deliberately uses the strategy action.
        jac_op->apply(*v, *diff);
        diff->subtract(diff, jv);
        const double initialized_action_error = diff->L2Norm();
        constexpr double INITIALIZED_ACTION_TOL = 1.0e-9;
        const double initialized_action_bound = INITIALIZED_ACTION_TOL * std::max(1.0, jv->L2Norm());
        const bool initialized_action_valid =
            std::isfinite(initialized_action_error) && initialized_action_error <= initialized_action_bound;
        pout << "initialized_jacobian_action_valid = " << (initialized_action_valid ? "true" : "false")
             << ", L2_bound = " << initialized_action_bound << std::endl;
        if (!initialized_action_valid)
        {
            pout << "Initialized Jacobian action failed: error = " << std::setprecision(17) << initialized_action_error
                 << ", bound = " << initialized_action_bound << std::setprecision(6) << std::endl;
            ++test_failures;
        }
        PetscErrorCode linear_ierr = KSPSetPCSide(linear_solver->getPETScKSP(), PC_RIGHT);
        IBTK_CHKERRQ(linear_ierr);
        linear_ierr = KSPSetNormType(linear_solver->getPETScKSP(), KSP_NORM_UNPRECONDITIONED);
        IBTK_CHKERRQ(linear_ierr);
        const bool linear_success = linear_solver->solveSystem(*linear_sol, *jv);
        if (!linear_success)
        {
            ++test_failures;
            pout << "krylov linear solve failed" << std::endl;
        }
        pout << "krylov_linear_iterations = " << linear_solver->getNumIterations() << std::endl;
        jac_op->apply(*linear_sol, *diff);
        diff->subtract(diff, jv);
        const double actual_residual = diff->L2Norm();
        const double residual_limit =
            std::max(linear_solver->getAbsoluteTolerance(), linear_solver->getRelativeTolerance() * jv->L2Norm());
        KSPConvergedReason reason;
        PCSide side;
        KSPNormType norm_type;
        KSPType ksp_type;
        linear_ierr = KSPGetConvergedReason(linear_solver->getPETScKSP(), &reason);
        IBTK_CHKERRQ(linear_ierr);
        linear_ierr = KSPGetPCSide(linear_solver->getPETScKSP(), &side);
        IBTK_CHKERRQ(linear_ierr);
        linear_ierr = KSPGetNormType(linear_solver->getPETScKSP(), &norm_type);
        IBTK_CHKERRQ(linear_ierr);
        linear_ierr = KSPGetType(linear_solver->getPETScKSP(), &ksp_type);
        IBTK_CHKERRQ(linear_ierr);
        const double residual_bound = std::max(1.0e-12, 2.0 * residual_limit);
        const bool physical_residual_valid = std::isfinite(actual_residual) && actual_residual <= residual_bound;
        const bool krylov_linear_residual_valid = std::string(ksp_type) == KSPFGMRES && side == PC_RIGHT &&
                                                  norm_type == KSP_NORM_UNPRECONDITIONED && reason > 0 &&
                                                  std::isfinite(linear_solver->getResidualNorm()) &&
                                                  linear_solver->getResidualNorm() >= 0.0 && physical_residual_valid;
        // The integration check requires a bounded physical residual, not a
        // particular over-converged residual from the preconditioned solve.
        pout << "physical_residual_valid = " << (physical_residual_valid ? "true" : "false")
             << ", acceptance_bound = " << residual_bound << ", pc_side = " << side << ", norm_type = " << norm_type
             << ", ksp_type = " << ksp_type << ", reason = " << reason << std::endl;
        if (!krylov_linear_residual_valid)
        {
            pout << "Krylov residual check failed: actual_residual = " << std::setprecision(17) << actual_residual
                 << ", reported_residual = " << linear_solver->getResidualNorm() << ", bound = " << residual_bound
                 << std::setprecision(6) << std::endl;
            ++test_failures;
        }

        double linear_side_norm = std::numeric_limits<double>::quiet_NaN();
        double linear_cell_norm = std::numeric_limits<double>::quiet_NaN();
        if (!side_l2_norm_is_finite(
                hier_velocity_data_ops, linear_sol->getComponentDescriptorIndex(0), wgt_sc_idx, linear_side_norm) ||
            !cell_l2_norm_is_finite(
                hier_pressure_data_ops, linear_sol->getComponentDescriptorIndex(1), wgt_cc_idx, linear_cell_norm))
        {
            ++test_failures;
            pout << "krylov linear solve produced non-finite norm" << std::endl;
        }
        else if (linear_side_norm <= 1.0e-14 && linear_cell_norm <= 1.0e-14)
        {
            ++test_failures;
            pout << "krylov linear solve action is trivial" << std::endl;
        }

        linear_solver->deallocateSolverState();
        fac_pc->deallocateSolverState();
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            PetscErrorCode ierr = MatDestroy(&coupling_reference[ln]);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&rediscretized_reference[ln]);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&operator_reference[ln]);
            IBTK_CHKERRQ(ierr);
        }

        if (verify_galerkin_operator_borrowing)
        {
            pout << "galerkin_operator_available_valid = " << (galerkin_operator_available_valid ? "true" : "false")
                 << std::endl;
            pout << "galerkin_operator_creator_lifetime_valid = "
                 << (galerkin_operator_creator_lifetime_valid ? "true" : "false") << std::endl;
            pout << "galerkin_operator_reinitialization_valid = " << (fac_reinitialization_valid ? "true" : "false")
                 << std::endl;
        }
        pout << "FAC Vec allocation measurement: calibrated and enforced with PETSC_USE_LOG; unavailable without it."
             << std::endl;
        pout << "fac_residual_repeat_error = " << fac_residual_work_vector_reuse_error
             << ", reinitialized_repeat_error = " << fac_residual_work_vector_reinitialize_error << std::endl;

        jac_op->deallocateOperatorState();
        nonlinear_op.deallocateOperatorState();

        ib_method_ops->postprocessIntegrateData(current_time, new_time, /*num_cycles*/ 1);

        for (auto vec : { nonlinear_probe, f_probe, v, jv, diff, linear_sol })
        {
            free_vector_components(*vec);
        }

        deallocate_vector_data(*eul_sol_vec);
        deallocate_vector_data(*eul_rhs_vec);
        free_vector_components(*eul_sol_vec);
        free_vector_components(*eul_rhs_vec);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int data_idx : allocated_patch_data_indices)
            {
                if (level->checkAllocated(data_idx))
                {
                    level->deallocatePatchData(data_idx);
                }
            }
        }

        PetscErrorCode ierr = MatDestroy(&A);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&J);
        IBTK_CHKERRQ(ierr);

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }

        pout << "test_failures = " << test_failures << std::endl;
    }

    return test_failures;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Logger::getInstance()->setWarning(false);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "components.log");
    const std::string test_case = app->getInputDatabase()->getString("test_case");
    if (test_case == "foundation")
    {
        return run_foundation(app);
    }
    if (test_case == "foundation_coarse_solver")
    {
        Pointer<Logger::Appender> abort_appender = new TestAppender();
        Logger::getInstance()->setAbortAppender(abort_appender);
        return run_foundation(app);
    }
    if (test_case == "foundation_wrong_strategy")
    {
        Pointer<Logger::Appender> abort_appender = new TestAppender();
        Logger::getInstance()->setAbortAppender(abort_appender);
        Pointer<FACPreconditionerStrategy> strategy =
            new StaggeredStokesLevelRelaxationFACOperator("stokes_fac", nullptr, "");
        StaggeredStokesIBJacobianFACPreconditioner solver("wrong_strategy", strategy, nullptr, "");
        return 0;
    }
    if (test_case == "operators")
    {
        return run_operators(app);
    }
    if (test_case == "level_borrowing")
    {
        return run_level_operator(app, false);
    }
    if (test_case == "level_augmentation")
    {
        return run_level_operator(app, true);
    }
    if (test_case == "distributed_augmentation")
    {
        return run_distributed_augmentation(app);
    }
    if (test_case == "set_operator_initialized")
    {
        return run_initialized_matrix_setter(app, false);
    }
    if (test_case == "set_augmentation_initialized")
    {
        return run_initialized_matrix_setter(app, true);
    }
    if (test_case == "level_state")
    {
        return run_level_state(app);
    }
    TBOX_ERROR("Unknown component test case: " << test_case << '\n');
    return 1;
}
