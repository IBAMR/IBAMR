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

/*
 * This test targets the implicit Stokes-IB solver components directly rather
 * than the full hierarchy integrator. It builds the nonlinear residual,
 * compares the Jacobian against a matrix-free finite-difference action, and
 * exercises the FAC-preconditioned linear solve on controlled single-level and
 * multilevel configurations.
 */

#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/StaggeredStokesIBJacobianFACPreconditioner.h>
#include <ibamr/StaggeredStokesIBJacobianOperator.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesIBOperator.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScMFFDJacobianOperator.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <IntVector.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SAMRAI_config.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include "../navier_stokes/cav_raw_operator_comparator.h"

#include <ibamr/app_namespaces.h>

namespace
{
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

// This opt-in diagnostic follows the rank-one pressure-CAV fixture in its production patch order. Reconstructing each
// residual from the live level operator independently checks the incremental update used by the smoother.
struct LocalSolveDiagnostics
{
    Mat level_operator = nullptr;
    const std::vector<std::vector<int>>* patch_dofs = nullptr;
    Vec sweep_rhs = nullptr;
    Vec accumulated_correction = nullptr;
    Vec fresh_residual = nullptr;
    Vec residual_difference = nullptr;
    std::size_t expected_ordinal = 0;
    int sweep_count = 0;
    int solve_count = 0;
    double max_backward_error = 0.0;
    double max_incremental_fresh_error = 0.0;
    double max_local_rhs_error = 0.0;
    bool order_valid = true;
    bool values_finite = true;
    bool trace_enabled = false;
    int trace_patch_ordinal = 0;
    bool trace_artifacts_round_trip = true;
    bool trace_index_round_trip = true;
    bool trace_selection_valid = true;
    bool trace_mutation_detected = true;
    std::vector<IBAMR::TestSupport::CAVLocalSolveTraceRecord> trace_records;

    void initialize(Mat operator_mat, const std::vector<std::vector<int>>& patches)
    {
        level_operator = operator_mat;
        patch_dofs = &patches;
        PetscErrorCode ierr = MatCreateVecs(level_operator, &accumulated_correction, &sweep_rhs);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(sweep_rhs, &fresh_residual);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(sweep_rhs, &residual_difference);
        IBTK_CHKERRQ(ierr);
    }

    void enableTrace(const int patch_ordinal)
    {
        if (patch_ordinal < -1) TBOX_ERROR("CAV local trace patch ordinal must be -1 or nonnegative\n");
        trace_enabled = true;
        trace_patch_ordinal = patch_ordinal;
    }

    void observe(const int ordinal, Mat local_matrix, Vec local_rhs, Vec local_solution, Vec current_global_source)
    {
        TBOX_ASSERT(level_operator);
        TBOX_ASSERT(patch_dofs);
        TBOX_ASSERT(!patch_dofs->empty());
        TBOX_ASSERT(ordinal >= 0 && static_cast<std::size_t>(ordinal) < patch_dofs->size());

        PetscErrorCode ierr = 0;
        if (ordinal == 0)
        {
            order_valid = order_valid && expected_ordinal == 0;
            expected_ordinal = 0;
            ++sweep_count;
            ierr = VecCopy(current_global_source, sweep_rhs);
            IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(accumulated_correction);
            IBTK_CHKERRQ(ierr);
        }
        order_valid = order_valid && static_cast<std::size_t>(ordinal) == expected_ordinal;

        if (trace_enabled && (trace_patch_ordinal == -1 || ordinal == trace_patch_ordinal))
        {
            const int sweep = sweep_count - 1;
            const std::string artifact_stem =
                "cav_fgmres_local_sweep" + std::to_string(sweep) + "_patch" + std::to_string(ordinal);
            const IBAMR::TestSupport::CAVLocalSolveTraceRecord record = { sweep, ordinal, artifact_stem };
            trace_artifacts_round_trip = trace_artifacts_round_trip &&
                                         IBAMR::TestSupport::writeCAVLocalSolveTraceArtifacts(
                                             record, local_matrix, local_rhs, local_solution, current_global_source);
            trace_records.push_back(record);
        }

        ierr = MatMult(level_operator, accumulated_correction, fresh_residual);
        IBTK_CHKERRQ(ierr);
        ierr = VecAYPX(fresh_residual, -1.0, sweep_rhs);
        IBTK_CHKERRQ(ierr);
        ierr = VecCopy(current_global_source, residual_difference);
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(residual_difference, -1.0, fresh_residual);
        IBTK_CHKERRQ(ierr);
        double difference_inf = 0.0;
        double incremental_inf = 0.0;
        double fresh_inf = 0.0;
        ierr = VecNorm(residual_difference, NORM_INFINITY, &difference_inf);
        IBTK_CHKERRQ(ierr);
        ierr = VecNorm(current_global_source, NORM_INFINITY, &incremental_inf);
        IBTK_CHKERRQ(ierr);
        ierr = VecNorm(fresh_residual, NORM_INFINITY, &fresh_inf);
        IBTK_CHKERRQ(ierr);
        const double incremental_fresh_error = difference_inf / std::max({ 1.0, incremental_inf, fresh_inf });
        max_incremental_fresh_error = std::max(max_incremental_fresh_error, incremental_fresh_error);

        const std::vector<int>& dofs = (*patch_dofs)[static_cast<std::size_t>(ordinal)];
        std::vector<PetscInt> petsc_dofs(dofs.begin(), dofs.end());
        std::vector<PetscScalar> restricted_source(dofs.size());
        ierr = VecGetValues(current_global_source,
                            static_cast<PetscInt>(petsc_dofs.size()),
                            petsc_dofs.data(),
                            restricted_source.data());
        IBTK_CHKERRQ(ierr);
        const PetscScalar* local_rhs_values = nullptr;
        ierr = VecGetArrayRead(local_rhs, &local_rhs_values);
        IBTK_CHKERRQ(ierr);
        for (std::size_t k = 0; k < dofs.size(); ++k)
        {
            max_local_rhs_error = std::max(
                max_local_rhs_error, static_cast<double>(PetscAbsScalar(restricted_source[k] - local_rhs_values[k])));
        }
        ierr = VecRestoreArrayRead(local_rhs, &local_rhs_values);
        IBTK_CHKERRQ(ierr);

        Vec local_defect = nullptr;
        ierr = VecDuplicate(local_rhs, &local_defect);
        IBTK_CHKERRQ(ierr);
        ierr = MatMult(local_matrix, local_solution, local_defect);
        IBTK_CHKERRQ(ierr);
        ierr = VecAYPX(local_defect, -1.0, local_rhs);
        IBTK_CHKERRQ(ierr);
        double matrix_inf = 0.0;
        double solution_inf = 0.0;
        double rhs_inf = 0.0;
        double defect_inf = 0.0;
        ierr = MatNorm(local_matrix, NORM_INFINITY, &matrix_inf);
        IBTK_CHKERRQ(ierr);
        ierr = VecNorm(local_solution, NORM_INFINITY, &solution_inf);
        IBTK_CHKERRQ(ierr);
        ierr = VecNorm(local_rhs, NORM_INFINITY, &rhs_inf);
        IBTK_CHKERRQ(ierr);
        ierr = VecNorm(local_defect, NORM_INFINITY, &defect_inf);
        IBTK_CHKERRQ(ierr);
        const double backward_error = defect_inf / (matrix_inf * solution_inf + rhs_inf + 1.0e-30);
        max_backward_error = std::max(max_backward_error, backward_error);
        values_finite = values_finite && std::isfinite(incremental_fresh_error) && std::isfinite(backward_error) &&
                        std::isfinite(max_local_rhs_error);
        ierr = VecDestroy(&local_defect);
        IBTK_CHKERRQ(ierr);

        const PetscScalar* local_solution_values = nullptr;
        ierr = VecGetArrayRead(local_solution, &local_solution_values);
        IBTK_CHKERRQ(ierr);
        ierr = VecSetValues(accumulated_correction,
                            static_cast<PetscInt>(petsc_dofs.size()),
                            petsc_dofs.data(),
                            local_solution_values,
                            ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreArrayRead(local_solution, &local_solution_values);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyBegin(accumulated_correction);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(accumulated_correction);
        IBTK_CHKERRQ(ierr);

        ++solve_count;
        expected_ordinal = (static_cast<std::size_t>(ordinal) + 1) % patch_dofs->size();
    }

    bool complete() const
    {
        return patch_dofs && sweep_count > 0 && expected_ordinal == 0 &&
               solve_count == sweep_count * static_cast<int>(patch_dofs->size());
    }

    void finalizeTrace()
    {
        if (!trace_enabled) return;
        IBAMR::TestSupport::writeCAVLocalSolveTraceIndex(trace_records, "cav_fgmres_local_trace.txt");
        const auto reread_records = IBAMR::TestSupport::readCAVLocalSolveTraceIndex("cav_fgmres_local_trace.txt");
        trace_index_round_trip = IBAMR::TestSupport::sameCAVLocalSolveTraceIndex(trace_records, reread_records);
        const std::size_t expected_count = trace_patch_ordinal == -1 ?
                                               static_cast<std::size_t>(sweep_count * patch_dofs->size()) :
                                               static_cast<std::size_t>(sweep_count);
        trace_selection_valid = trace_records.size() == expected_count;
        auto mutated_records = reread_records;
        if (!mutated_records.empty()) ++mutated_records.front().patch_ordinal;
        trace_mutation_detected = !IBAMR::TestSupport::sameCAVLocalSolveTraceIndex(trace_records, mutated_records);
    }

    void deallocate()
    {
        PetscErrorCode ierr = VecDestroy(&residual_difference);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&fresh_residual);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&accumulated_correction);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&sweep_rhs);
        IBTK_CHKERRQ(ierr);
        level_operator = nullptr;
        patch_dofs = nullptr;
    }
};

std::string
fac_stage_name(const IBTK::FACPreconditioner::CycleStage stage)
{
    using CycleStage = IBTK::FACPreconditioner::CycleStage;
    switch (stage)
    {
    case CycleStage::PRE_SMOOTH_INPUT:
        return "pre-smooth-input";
    case CycleStage::PRE_SMOOTH_OUTPUT:
        return "pre-smooth-output";
    case CycleStage::COARSE_RHS:
        return "coarse-rhs";
    case CycleStage::COARSE_CORRECTION:
        return "coarse-correction";
    case CycleStage::POST_SMOOTH_INPUT:
        return "post-smooth-input";
    case CycleStage::POST_SMOOTH_OUTPUT:
        return "post-smooth-output";
    }
    TBOX_ERROR("Unsupported FAC cycle stage\n");
    return "";
}

struct FACTraceExporter
{
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_operator;
    Pointer<PatchHierarchy<NDIM>> hierarchy;
    int u_dof_index_idx = -1;
    int p_dof_index_idx = -1;
    bool active = false;
    bool artifacts_round_trip = true;
    bool index_round_trip = true;
    bool order_valid = true;
    bool mutation_detected = true;
    std::vector<IBAMR::TestSupport::CAVFACStageTraceRecord> records;

    void observe(const IBTK::FACPreconditioner::CycleStage stage,
                 const int level_num,
                 const SAMRAIVectorReal<NDIM, double>& solution,
                 const SAMRAIVectorReal<NDIM, double>& rhs)
    {
        if (!active) return;
        Pointer<StaggeredStokesPETScLevelSolver> level_solver =
            fac_operator->getStaggeredStokesPETScLevelSolver(level_num);
        const auto live_view = level_solver->getLiveOperatorStateView();
        if (!live_view.initialized || !live_view.operator_mat)
            TBOX_ERROR("FAC trace requires an initialized live level operator\n");

        // FAC stage vectors are hierarchy-native. Use the existing level DOF map for the transient serialized view;
        // do not introduce a second hierarchy-wide numbering or retain an algebraic copy.
        Vec solution_vec = nullptr;
        Vec rhs_vec = nullptr;
        PetscErrorCode ierr = MatCreateVecs(live_view.operator_mat, &solution_vec, &rhs_vec);
        IBTK_CHKERRQ(ierr);
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(solution_vec,
                                                              solution.getComponentDescriptorIndex(0),
                                                              u_dof_index_idx,
                                                              solution.getComponentDescriptorIndex(1),
                                                              p_dof_index_idx,
                                                              hierarchy->getPatchLevel(level_num));
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(rhs_vec,
                                                              rhs.getComponentDescriptorIndex(0),
                                                              u_dof_index_idx,
                                                              rhs.getComponentDescriptorIndex(1),
                                                              p_dof_index_idx,
                                                              hierarchy->getPatchLevel(level_num));

        const std::string stage_name = fac_stage_name(stage);
        const std::string artifact_stem =
            "cav_fac_stage" + std::to_string(records.size()) + "_" + stage_name + "_level" + std::to_string(level_num);
        IBAMR::TestSupport::writeCAVRawVectorMarket(solution_vec, artifact_stem + "_solution.mtx");
        IBAMR::TestSupport::writeCAVRawVectorMarket(rhs_vec, artifact_stem + "_rhs.mtx");
        artifacts_round_trip =
            artifacts_round_trip &&
            IBAMR::TestSupport::sameCAVRawVectorMarket(
                solution_vec, IBAMR::TestSupport::readCAVRawVectorMarket(artifact_stem + "_solution.mtx")) &&
            IBAMR::TestSupport::sameCAVRawVectorMarket(
                rhs_vec, IBAMR::TestSupport::readCAVRawVectorMarket(artifact_stem + "_rhs.mtx"));
        records.push_back({ stage_name, level_num, artifact_stem });

        ierr = VecDestroy(&rhs_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&solution_vec);
        IBTK_CHKERRQ(ierr);
    }

    void finalize(const std::vector<IBTK::FACPreconditioner::CycleStage>& expected_stages,
                  const std::vector<int>& expected_levels)
    {
        IBAMR::TestSupport::writeCAVFACStageTraceIndex(records, "cav_fac_trace.txt");
        const auto reread_records = IBAMR::TestSupport::readCAVFACStageTraceIndex("cav_fac_trace.txt");
        index_round_trip = IBAMR::TestSupport::sameCAVFACStageTraceIndex(records, reread_records);
        order_valid = records.size() == expected_stages.size() && expected_stages.size() == expected_levels.size();
        for (std::size_t k = 0; order_valid && k < records.size(); ++k)
        {
            order_valid =
                records[k].stage == fac_stage_name(expected_stages[k]) && records[k].level == expected_levels[k];
        }
        auto mutated_records = reread_records;
        if (!mutated_records.empty()) ++mutated_records.front().level;
        mutation_detected = !IBAMR::TestSupport::sameCAVFACStageTraceIndex(records, mutated_records);
    }
};

struct KrylovTraceContext
{
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_operator;
    Pointer<PatchHierarchy<NDIM>> hierarchy;
    int u_dof_index_idx = -1;
    int p_dof_index_idx = -1;
    bool artifacts_round_trip = true;
    std::vector<IBAMR::TestSupport::CAVKrylovTraceRecord> records;
};

PetscErrorCode
record_krylov_trace(KSP ksp, const PetscInt iteration, const PetscReal residual_norm, void* context)
{
    PetscFunctionBeginUser;
    auto* trace = static_cast<KrylovTraceContext*>(context);
    Vec solution = nullptr;
    PetscCall(KSPBuildSolution(ksp, nullptr, &solution));
    PetscBool is_samrai_vector = PETSC_FALSE;
    PetscCall(PetscObjectTypeCompare(reinterpret_cast<PetscObject>(solution), "Vec_SAMRAI", &is_samrai_vector));
    PetscCheck(is_samrai_vector,
               PetscObjectComm(reinterpret_cast<PetscObject>(ksp)),
               PETSC_ERR_ARG_WRONG,
               "CAV Krylov trace requires the live Vec_SAMRAI iterate");
    const std::string artifact_stem = "cav_fgmres_iteration" + std::to_string(iteration);
    // The outer iterate has no single contiguous hierarchy numbering. Unwrap the live hierarchy vector and export each
    // level through the same velocity-pressure global numbering used by its live level operator.
    Pointer<SAMRAIVectorReal<NDIM, PetscScalar>> hierarchy_solution;
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVectorRead(solution, &hierarchy_solution);
    try
    {
        for (int ln = hierarchy_solution->getCoarsestLevelNumber(); ln <= hierarchy_solution->getFinestLevelNumber();
             ++ln)
        {
            Pointer<StaggeredStokesPETScLevelSolver> level_solver =
                trace->fac_operator->getStaggeredStokesPETScLevelSolver(ln);
            const auto live_view = level_solver->getLiveOperatorStateView();
            Vec level_solution = nullptr;
            PetscCall(MatCreateVecs(live_view.operator_mat, &level_solution, nullptr));
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(level_solution,
                                                                  hierarchy_solution->getComponentDescriptorIndex(0),
                                                                  trace->u_dof_index_idx,
                                                                  hierarchy_solution->getComponentDescriptorIndex(1),
                                                                  trace->p_dof_index_idx,
                                                                  trace->hierarchy->getPatchLevel(ln));
            const std::string filename = artifact_stem + "_level" + std::to_string(ln) + "_solution.mtx";
            IBAMR::TestSupport::writeCAVRawVectorMarket(level_solution, filename);
            trace->artifacts_round_trip = trace->artifacts_round_trip &&
                                          IBAMR::TestSupport::sameCAVRawVectorMarket(
                                              level_solution, IBAMR::TestSupport::readCAVRawVectorMarket(filename));
            PetscCall(VecDestroy(&level_solution));
        }
        trace->records.push_back({ static_cast<int>(iteration), static_cast<double>(residual_norm), artifact_stem });
    }
    catch (const std::exception& exception)
    {
        IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(solution, &hierarchy_solution);
        SETERRQ(PetscObjectComm(reinterpret_cast<PetscObject>(ksp)), PETSC_ERR_FILE_WRITE, "%s", exception.what());
    }
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(solution, &hierarchy_solution);
    PetscFunctionReturn(PETSC_SUCCESS);
}

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
    if (ln != spec->finest_ln || strct_num != 0) return;

    for (int k = 0; k < spec->num_curve_points; ++k)
    {
        IBRedundantInitializer::Edge edge = { k, (k + 1) % spec->num_curve_points };
        if (edge.first > edge.second) std::swap(edge.first, edge.second);
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
#if (NDIM == 2)
                    const double val = (axis == 0) ? std::sin(2.0 * M_PI * x[1]) : -std::sin(2.0 * M_PI * x[0]);
#else
                    const double val = (axis == 0) ?
                                           std::sin(2.0 * M_PI * x[1]) :
                                           ((axis == 1) ? std::sin(2.0 * M_PI * x[2]) : std::sin(2.0 * M_PI * x[0]));
#endif
                    (*u_data)(i_s) = val;
                }
            }
        }
    }
}

} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Suppress warnings so the test output does not depend on build-path
    // details embedded in warning diagnostics.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);

    int test_failures = 0;

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const double current_time = 0.0;
        const double dt = input_db->getDoubleWithDefault("DT", 0.005);
        const double rho = input_db->getDoubleWithDefault("RHO", 1.0);
        const double mu = input_db->getDoubleWithDefault("MU", 1.0);
        const double new_time = current_time + dt;
        const bool use_fixed_le_operators = input_db->getBoolWithDefault("USE_FIXED_LE_OPERATORS", true);
        const bool verify_fgmres_physical_residuals = input_db->keyExists("VERIFY_FGMRES_PHYSICAL_RESIDUALS") ?
                                                          input_db->getBool("VERIFY_FGMRES_PHYSICAL_RESIDUALS") :
                                                          false;
        const bool verify_fgmres_local_diagnostics = input_db->keyExists("VERIFY_FGMRES_LOCAL_DIAGNOSTICS") ?
                                                         input_db->getBool("VERIFY_FGMRES_LOCAL_DIAGNOSTICS") :
                                                         false;
        const bool verify_cav_live_dynamic_trace_schema =
            input_db->keyExists("VERIFY_CAV_LIVE_DYNAMIC_TRACE_SCHEMA") ?
                input_db->getBool("VERIFY_CAV_LIVE_DYNAMIC_TRACE_SCHEMA") :
                false;
        const int cav_local_trace_patch_ordinal =
            verify_cav_live_dynamic_trace_schema ? input_db->getInteger("CAV_LOCAL_TRACE_PATCH_ORDINAL") : 0;
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
        const IntVector<NDIM> no_ghosts = IntVector<NDIM>(0);

        const int u_current_idx = var_db->registerVariableAndContext(u_var, current_ctx, ib_ghosts);
        const int u_sol_idx = var_db->registerVariableAndContext(u_var, solver_ctx, one_ghost);
        const int f_rhs_idx = var_db->registerVariableAndContext(f_var, solver_ctx, one_ghost);
        const int p_sol_idx = var_db->registerVariableAndContext(p_var, solver_ctx, one_ghost);
        const int g_rhs_idx = var_db->registerVariableAndContext(g_var, solver_ctx, one_ghost);
        const int u_scratch_idx = var_db->registerVariableAndContext(u_var, scratch_ctx, ib_ghosts);
        const int f_scratch_idx = var_db->registerVariableAndContext(f_var, scratch_ctx, ib_ghosts);
        const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, scratch_ctx, ib_ghosts);
        const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, scratch_ctx, no_ghosts);

        const std::vector<int> allocated_patch_data_indices = {
            u_current_idx, u_scratch_idx, f_scratch_idx, u_dof_index_idx, p_dof_index_idx
        };
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int data_idx : allocated_patch_data_indices) level->allocatePatchData(data_idx, current_time);
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
        ib_method_ops->constructInterpOp(J,
                                         PETScMatUtilities::ib_4_delta_fcn,
                                         PETScMatUtilities::ib_4_delta_stencil,
                                         num_dofs_per_proc[finest_ln],
                                         u_dof_index_idx,
                                         new_time);

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
        ctx.patch_level = patch_hierarchy->getPatchLevel(finest_ln);
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
        jac_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);
        jac_op->formJacobian(*eul_sol_vec);

        Pointer<SAMRAIVectorReal<NDIM, double>> v = eul_sol_vec->cloneVector("v");
        v->allocateVectorData();
        v->setToScalar(0.0);
        if (verify_fgmres_physical_residuals)
        {
            // A rigid translation lies in the elastic coupling nullspace and
            // cannot provide a well-scaled end-to-end IB residual check.
            set_divergence_free_probe_velocity(v->getComponentDescriptorIndex(0), patch_hierarchy);
        }
        else
        {
            hier_velocity_data_ops->setToScalar(v->getComponentDescriptorIndex(0), 1.0, false);
        }
        hier_pressure_data_ops->setToScalar(v->getComponentDescriptorIndex(1), -0.25, false);

        Pointer<SAMRAIVectorReal<NDIM, double>> jv = eul_rhs_vec->cloneVector("jv");
        jv->allocateVectorData();
        jv->setToScalar(0.0);
        jac_op->apply(*v, *jv);

        const double fd_rel_tol = input_db->getDoubleWithDefault("FD_REL_TOL", 5.0e-2);
        Pointer<PETScMFFDJacobianOperator> mffd_jac_op =
            new PETScMFFDJacobianOperator("stokes_ib_solver_components::mffd_jacobian_op", "ib_jac_mffd_");
        mffd_jac_op->setOperator(Pointer<GeneralOperator>(&nonlinear_op, false));
        mffd_jac_op->setTimeInterval(current_time, new_time);
        mffd_jac_op->setSolutionTime(new_time);
        mffd_jac_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);
        mffd_jac_op->formJacobian(*eul_sol_vec);

        Pointer<SAMRAIVectorReal<NDIM, double>> fd_jv = eul_rhs_vec->cloneVector("fd_jv");
        fd_jv->allocateVectorData();
        fd_jv->setToScalar(0.0);
        mffd_jac_op->apply(*v, *fd_jv);

        Pointer<SAMRAIVectorReal<NDIM, double>> diff = eul_rhs_vec->cloneVector("diff");
        diff->allocateVectorData();
        diff->subtract(fd_jv, jv);

        double jv_side_norm = std::numeric_limits<double>::quiet_NaN();
        double jv_cell_norm = std::numeric_limits<double>::quiet_NaN();
        double diff_side_norm = std::numeric_limits<double>::quiet_NaN();
        double diff_cell_norm = std::numeric_limits<double>::quiet_NaN();
        const bool jv_finite =
            side_l2_norm_is_finite(
                hier_velocity_data_ops, jv->getComponentDescriptorIndex(0), wgt_sc_idx, jv_side_norm) &&
            cell_l2_norm_is_finite(
                hier_pressure_data_ops, jv->getComponentDescriptorIndex(1), wgt_cc_idx, jv_cell_norm);
        const bool diff_finite =
            side_l2_norm_is_finite(
                hier_velocity_data_ops, diff->getComponentDescriptorIndex(0), wgt_sc_idx, diff_side_norm) &&
            cell_l2_norm_is_finite(
                hier_pressure_data_ops, diff->getComponentDescriptorIndex(1), wgt_cc_idx, diff_cell_norm);
        if (!jv_finite || !diff_finite)
        {
            ++test_failures;
            pout << "jacobian norms are non-finite" << std::endl;
        }
        else if (jv_side_norm <= 1.0e-14 && jv_cell_norm <= 1.0e-14)
        {
            ++test_failures;
            pout << "jacobian action is trivial" << std::endl;
        }
        else
        {
            const double rel_error =
                std::sqrt(diff_side_norm * diff_side_norm + diff_cell_norm * diff_cell_norm) /
                std::max(std::sqrt(jv_side_norm * jv_side_norm + jv_cell_norm * jv_cell_norm), 1.0e-14);
            pout << "fd_relative_error = " << rel_error << std::endl;
            if (!(rel_error <= fd_rel_tol))
            {
                ++test_failures;
                pout << "fd_relative_error exceeds tolerance: " << fd_rel_tol << std::endl;
            }
        }
        mffd_jac_op->deallocateOperatorState();

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
        const bool verify_pressure_cav =
            input_db->keyExists("VERIFY_PRESSURE_CAV") ? input_db->getBool("VERIFY_PRESSURE_CAV") : false;
        const bool verify_fac_cycle_observer =
            input_db->keyExists("VERIFY_FAC_CYCLE_OBSERVER") ? input_db->getBool("VERIFY_FAC_CYCLE_OBSERVER") : false;
        const bool verify_galerkin_operator_borrowing = input_db->keyExists("VERIFY_GALERKIN_OPERATOR_BORROWING") ?
                                                            input_db->getBool("VERIFY_GALERKIN_OPERATOR_BORROWING") :
                                                            false;
        const bool verify_cav_live_export_schema = input_db->keyExists("VERIFY_CAV_LIVE_EXPORT_SCHEMA") ?
                                                       input_db->getBool("VERIFY_CAV_LIVE_EXPORT_SCHEMA") :
                                                       false;
        if (verify_cav_live_export_schema && !verify_pressure_cav)
            TBOX_ERROR("The CAV live-export schema check requires VERIFY_PRESSURE_CAV = TRUE\n");
        if (verify_cav_live_dynamic_trace_schema &&
            (!verify_pressure_cav || !verify_fac_cycle_observer || !verify_fgmres_physical_residuals))
        {
            TBOX_ERROR(
                "The CAV dynamic-trace schema check requires pressure CAV, the FAC observer, FGMRES physical "
                "residuals\n");
        }
        if (verify_galerkin_operator_borrowing && !verify_fac_cycle_observer)
            TBOX_ERROR("The Galerkin borrowing reinitialization check requires VERIFY_FAC_CYCLE_OBSERVER = TRUE\n");
        Pointer<StaggeredStokesPETScLevelSolver> pressure_cav_level_solver;
        bool dynamic_level_operator_round_trip = !verify_cav_live_dynamic_trace_schema;
        bool dynamic_level_dof_map_round_trip = !verify_cav_live_dynamic_trace_schema;
        bool dynamic_level_cav_data_round_trip = !verify_cav_live_dynamic_trace_schema;
        auto galerkin_operator_borrowing_is_valid = [&]()
        {
            bool valid = finest_ln > 0;
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                const auto live_view = fac_op->getStaggeredStokesPETScLevelSolver(ln)->getLiveOperatorStateView();
                const bool expect_galerkin_operator = ln != finest_ln;
                valid = valid && live_view.initialized && live_view.operator_mat &&
                        live_view.operator_was_provided == expect_galerkin_operator &&
                        live_view.includes_augmented_operator != expect_galerkin_operator;
            }
            return valid;
        };
        fac_pc->initializeSolverState(*eul_sol_vec, *eul_rhs_vec);
        bool galerkin_operator_borrowing_valid = !verify_galerkin_operator_borrowing;
        if (verify_galerkin_operator_borrowing)
        {
            galerkin_operator_borrowing_valid = galerkin_operator_borrowing_is_valid();
            if (!galerkin_operator_borrowing_valid) ++test_failures;
        }
        if (verify_pressure_cav)
        {
            pressure_cav_level_solver = fac_op->getStaggeredStokesPETScLevelSolver(finest_ln);
            const auto& overlap_dofs = pressure_cav_level_solver->getASMSubdomains();
            const auto& nonoverlap_dofs = pressure_cav_level_solver->getASMNonoverlapSubdomains();
            const auto& pressure_seeds = pressure_cav_level_solver->getCouplingAwareASMPressureSeedDOFs();
            const int expected_patch_count = input_db->getInteger("EXPECTED_PRESSURE_CAV_PATCH_COUNT");
            if (verify_cav_live_dynamic_trace_schema &&
                cav_local_trace_patch_ordinal >= static_cast<int>(overlap_dofs.size()))
                TBOX_ERROR("CAV local trace patch ordinal is outside the live patch range\n");

            bool seed_order_valid = static_cast<int>(pressure_seeds.size()) == expected_patch_count &&
                                    pressure_seeds.size() == overlap_dofs.size();
            int enlarged_patch_count = 0;
            for (std::size_t k = 0; k < std::min(pressure_seeds.size(), overlap_dofs.size()); ++k)
            {
                seed_order_valid =
                    seed_order_valid &&
                    std::binary_search(overlap_dofs[k].begin(), overlap_dofs[k].end(), pressure_seeds[k]);
                if (overlap_dofs[k].size() > 2 * NDIM + 1) ++enlarged_patch_count;
            }
            const bool partition_absent = nonoverlap_dofs.empty();
            if (!seed_order_valid || enlarged_patch_count == 0 || !partition_absent) ++test_failures;

            pout << "pressure_cav_patch_count = " << pressure_seeds.size() << std::endl;
            pout << "pressure_cav_enlarged_patch_count = " << enlarged_patch_count << std::endl;
            pout << "pressure_cav_seed_order_valid = " << (seed_order_valid ? "true" : "false") << std::endl;
            pout << "pressure_cav_partition_absent = " << (partition_absent ? "true" : "false") << std::endl;
        }
        if (verify_cav_live_dynamic_trace_schema)
        {
            dynamic_level_operator_round_trip = true;
            dynamic_level_dof_map_round_trip = true;
            dynamic_level_cav_data_round_trip = true;
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                Pointer<StaggeredStokesPETScLevelSolver> level_solver = fac_op->getStaggeredStokesPETScLevelSolver(ln);
                const auto live_view = level_solver->getLiveOperatorStateView();
                const std::string prefix = "cav_dynamic_level" + std::to_string(ln);
                IBAMR::TestSupport::writeCAVRawMatrixMarket(live_view.operator_mat, prefix + "_A.mtx");
                dynamic_level_operator_round_trip =
                    dynamic_level_operator_round_trip &&
                    IBAMR::TestSupport::sameCAVRawMatrixMarket(
                        live_view.operator_mat, IBAMR::TestSupport::readCAVRawMatrixMarket(prefix + "_A.mtx"));
                const auto dof_map = IBAMR::TestSupport::captureCAVRawDofMap(
                    patch_hierarchy->getPatchLevel(ln), u_dof_index_idx, p_dof_index_idx);
                IBAMR::TestSupport::writeCAVRawDofMap(dof_map, prefix + "_dof_map.txt");
                dynamic_level_dof_map_round_trip =
                    dynamic_level_dof_map_round_trip &&
                    IBAMR::TestSupport::sameCAVRawDofMap(dof_map,
                                                         IBAMR::TestSupport::readCAVRawDofMap(prefix + "_dof_map.txt"));
                if (ln > 0)
                {
                    const auto& pressure_seeds = level_solver->getCouplingAwareASMPressureSeedDOFs();
                    const auto& overlap_dofs = level_solver->getASMSubdomains();
                    IBAMR::TestSupport::writeCAVRawIndexList(pressure_seeds, prefix + "_pressure_seeds.txt");
                    IBAMR::TestSupport::writeCAVRawIndexSets(overlap_dofs, prefix + "_patches.txt");
                    dynamic_level_cav_data_round_trip =
                        dynamic_level_cav_data_round_trip && !pressure_seeds.empty() &&
                        pressure_seeds.size() == overlap_dofs.size() &&
                        pressure_seeds ==
                            IBAMR::TestSupport::readCAVRawIndexList(prefix + "_pressure_seeds.txt") &&
                        overlap_dofs == IBAMR::TestSupport::readCAVRawIndexSets(prefix + "_patches.txt");
                }
            }
            if (!dynamic_level_operator_round_trip || !dynamic_level_dof_map_round_trip ||
                !dynamic_level_cav_data_round_trip)
                ++test_failures;
        }
        Mat SAJ = fac_op->getEulerianElasticityLevelOp(finest_ln);
        if (verify_cav_live_export_schema)
        {
            using namespace IBAMR::TestSupport;
            const auto live_view = pressure_cav_level_solver->getLiveOperatorStateView();
            const auto& pressure_seeds = pressure_cav_level_solver->getCouplingAwareASMPressureSeedDOFs();
            const auto& overlap_dofs = pressure_cav_level_solver->getASMSubdomains();
            const auto& nonoverlap_dofs = pressure_cav_level_solver->getASMNonoverlapSubdomains();
            const std::string prefix = "cav_live_export_contract";

            writeCAVRawMatrixMarket(live_view.operator_mat, prefix + "_A.mtx");
            writeCAVRawMatrixMarket(SAJ, prefix + "_E_h.mtx");
            const CAVRawDofMapData dof_map =
                captureCAVRawDofMap(patch_hierarchy->getPatchLevel(finest_ln), u_dof_index_idx, p_dof_index_idx);
            writeCAVRawDofMap(dof_map, prefix + "_dof_map.txt");
            writeCAVRawIndexList(pressure_seeds, prefix + "_pressure_seeds.txt");
            writeCAVRawIndexSets(overlap_dofs, prefix + "_patches.txt");
            writeCAVRawIndexSets(nonoverlap_dofs, prefix + "_nonoverlap.txt");

            Pointer<Database> level_solver_db = stokes_ib_precond_db->getDatabase("level_solver_db");
            if (level_solver_db->getString("coupling_aware_asm_patch_seed_type") != "PRESSURE_CELL")
            {
                TBOX_ERROR("The native live-export contract requires pressure-cell CAV construction\n");
            }
            const std::string shell_pc_type = level_solver_db->getString("shell_pc_type");
            std::string local_solver_backend;
            if (shell_pc_type == "multiplicative-blas-lapack")
            {
                const std::string solver_type = level_solver_db->getString("blas_lapack_subdomain_solver_type");
                if (solver_type != "lu" && solver_type != "svd")
                    TBOX_ERROR("The native CAV backend export supports only BLAS/LAPACK LU and SVD\n");
                local_solver_backend = "blas-lapack-" + solver_type;
            }
            else if (shell_pc_type == "multiplicative-eigen-reference")
            {
                const std::string solver_type =
                    level_solver_db->getStringWithDefault("eigen_subdomain_solver_type", "partial-piv-lu");
                local_solver_backend = "eigen-reference-" + solver_type;
            }
            else if (shell_pc_type == "multiplicative-eigen-schur-complement")
            {
                local_solver_backend = "eigen-schur-complement";
            }
            else
            {
                TBOX_ERROR("Unsupported native multiplicative CAV export backend: " << shell_pc_type << "\n");
            }
            CAVLiveExportManifest manifest;
            manifest.candidate_sha = "0123456789abcdef0123456789abcdef01234567";
            manifest.oracle_sha = "5b77344db6746269f8c77695c99e9043907ba74b";
            manifest.case_id = "native-live-construction-export";
            manifest.mpi_ranks = IBTK_MPI::getNodes();
            manifest.pressure_equation = "minus-div";
            manifest.pressure_equation_row_multiplier_to_oracle = -1.0;
            manifest.pressure_gauge = "zero-mean-correction";
            manifest.patch_seed_type = "PRESSURE_CELL";
            manifest.closure_policy = level_solver_db->getString("coupling_aware_asm_closure_policy");
            manifest.seed_stride = level_solver_db->getInteger("coupling_aware_asm_seed_stride");
            manifest.traversal_order = level_solver_db->getString("coupling_aware_asm_seed_traversal_order");
            manifest.composition = "multiplicative";
            manifest.local_solver_backend = local_solver_backend;
            writeCAVLiveExportManifest(manifest, prefix + "_manifest.txt");

            const bool operator_round_trip =
                sameCAVRawMatrixMarket(live_view.operator_mat, readCAVRawMatrixMarket(prefix + "_A.mtx"));
            const bool elasticity_round_trip = sameCAVRawMatrixMarket(SAJ, readCAVRawMatrixMarket(prefix + "_E_h.mtx"));
            const bool dof_map_round_trip = sameCAVRawDofMap(dof_map, readCAVRawDofMap(prefix + "_dof_map.txt"));
            const bool pressure_seed_round_trip = pressure_seeds == readCAVRawIndexList(prefix + "_pressure_seeds.txt");
            const bool patch_round_trip = overlap_dofs == readCAVRawIndexSets(prefix + "_patches.txt");
            const bool nonoverlap_round_trip = nonoverlap_dofs == readCAVRawIndexSets(prefix + "_nonoverlap.txt");
            const bool manifest_round_trip =
                sameCAVLiveExportManifest(manifest, readCAVLiveExportManifest(prefix + "_manifest.txt"));

            bool mutations_detected = !dof_map.dofs.empty() && pressure_seeds.size() >= 2 && !overlap_dofs.empty() &&
                                      !overlap_dofs.front().empty();
            if (mutations_detected)
            {
                CAVRawMatrixMarketData omitted_operator_entry = readCAVRawMatrixMarket(prefix + "_A.mtx");
                mutations_detected = !omitted_operator_entry.entries.empty();
                if (mutations_detected)
                {
                    omitted_operator_entry.entries.pop_back();
                    CAVRawDofMapData altered_dof_map = dof_map;
                    altered_dof_map.dofs.front().axis = altered_dof_map.dofs.front().axis == 0 ? 1 : 0;
                    std::vector<int> reordered_seeds = pressure_seeds;
                    std::swap(reordered_seeds[0], reordered_seeds[1]);
                    std::vector<std::vector<int>> omitted_patch_dof = overlap_dofs;
                    omitted_patch_dof.front().pop_back();
                    std::vector<std::vector<int>> altered_nonoverlap = nonoverlap_dofs;
                    altered_nonoverlap.push_back({ pressure_seeds.front() });
                    CAVLiveExportManifest altered_manifest = manifest;
                    altered_manifest.pressure_equation_row_multiplier_to_oracle = 1.0;
                    mutations_detected = !sameCAVRawMatrixMarket(live_view.operator_mat, omitted_operator_entry) &&
                                         !sameCAVRawDofMap(dof_map, altered_dof_map) &&
                                         reordered_seeds != pressure_seeds && omitted_patch_dof != overlap_dofs &&
                                         altered_nonoverlap != nonoverlap_dofs &&
                                         !sameCAVLiveExportManifest(manifest, altered_manifest);
                }
            }

            if (!operator_round_trip || !elasticity_round_trip || !dof_map_round_trip || !pressure_seed_round_trip ||
                !patch_round_trip || !nonoverlap_round_trip || !manifest_round_trip || !mutations_detected)
            {
                ++test_failures;
            }
            pout << "cav_live_export_operator_round_trip = " << (operator_round_trip ? "true" : "false") << std::endl;
            pout << "cav_live_export_elasticity_round_trip = " << (elasticity_round_trip ? "true" : "false")
                 << std::endl;
            pout << "cav_live_export_dof_map_round_trip = " << (dof_map_round_trip ? "true" : "false") << std::endl;
            pout << "cav_live_export_pressure_seed_round_trip = " << (pressure_seed_round_trip ? "true" : "false")
                 << std::endl;
            pout << "cav_live_export_patch_round_trip = " << (patch_round_trip ? "true" : "false") << std::endl;
            pout << "cav_live_export_nonoverlap_round_trip = " << (nonoverlap_round_trip ? "true" : "false")
                 << std::endl;
            pout << "cav_live_export_manifest_round_trip = " << (manifest_round_trip ? "true" : "false") << std::endl;
            pout << "cav_live_export_mutations_detected = " << (mutations_detected ? "true" : "false") << std::endl;
        }
        jac_op->setIBCouplingJacobian(SAJ);

        using CycleStage = IBTK::FACPreconditioner::CycleStage;
        std::vector<CycleStage> observed_fac_stages;
        std::vector<int> observed_fac_levels;
        bool fac_observer_views_valid = true;
        FACTraceExporter fac_trace_exporter;
        if (verify_cav_live_dynamic_trace_schema)
        {
            fac_trace_exporter.fac_operator = fac_op;
            fac_trace_exporter.hierarchy = patch_hierarchy;
            fac_trace_exporter.u_dof_index_idx = u_dof_index_idx;
            fac_trace_exporter.p_dof_index_idx = p_dof_index_idx;
        }
        IBTK::FACPreconditioner::CycleObserver fac_cycle_observer = [&](const CycleStage stage,
                                                                        const int level_num,
                                                                        const SAMRAIVectorReal<NDIM, double>& solution,
                                                                        const SAMRAIVectorReal<NDIM, double>& rhs)
        {
            observed_fac_stages.push_back(stage);
            observed_fac_levels.push_back(level_num);
            const double solution_norm = solution.L2Norm();
            const double rhs_norm = rhs.L2Norm();
            fac_observer_views_valid = fac_observer_views_valid && solution.getPatchHierarchy() == patch_hierarchy &&
                                       rhs.getPatchHierarchy() == patch_hierarchy &&
                                       solution.getCoarsestLevelNumber() == 0 && rhs.getCoarsestLevelNumber() == 0 &&
                                       solution.getFinestLevelNumber() == finest_ln &&
                                       rhs.getFinestLevelNumber() == finest_ln && std::isfinite(solution_norm) &&
                                       std::isfinite(rhs_norm);
            if (!std::isfinite(solution_norm) || !std::isfinite(rhs_norm))
            {
                pout << "fac_cycle_observer_nonfinite_stage = " << static_cast<int>(stage) << std::endl;
                pout << "fac_cycle_observer_nonfinite_level = " << level_num << std::endl;
                pout << "fac_cycle_observer_solution_norm = " << solution_norm << std::endl;
                pout << "fac_cycle_observer_rhs_norm = " << rhs_norm << std::endl;
            }
            fac_trace_exporter.observe(stage, level_num, solution, rhs);
        };

        std::vector<CycleStage> expected_fac_stages;
        std::vector<int> expected_fac_levels;
        for (int ln = finest_ln; ln > 0; --ln)
        {
            expected_fac_stages.push_back(CycleStage::PRE_SMOOTH_INPUT);
            expected_fac_stages.push_back(CycleStage::PRE_SMOOTH_OUTPUT);
            expected_fac_levels.push_back(ln);
            expected_fac_levels.push_back(ln);
        }
        expected_fac_stages.push_back(CycleStage::COARSE_RHS);
        expected_fac_stages.push_back(CycleStage::COARSE_CORRECTION);
        expected_fac_levels.push_back(0);
        expected_fac_levels.push_back(0);
        for (int ln = 1; ln <= finest_ln; ++ln)
        {
            expected_fac_stages.push_back(CycleStage::POST_SMOOTH_INPUT);
            expected_fac_stages.push_back(CycleStage::POST_SMOOTH_OUTPUT);
            expected_fac_levels.push_back(ln);
            expected_fac_levels.push_back(ln);
        }
        std::vector<CycleStage> expected_fac_no_pre_stages = { CycleStage::COARSE_RHS,
                                                                CycleStage::COARSE_CORRECTION };
        std::vector<int> expected_fac_no_pre_levels = { 0, 0 };
        for (int ln = 1; ln <= finest_ln; ++ln)
        {
            expected_fac_no_pre_stages.push_back(CycleStage::POST_SMOOTH_INPUT);
            expected_fac_no_pre_stages.push_back(CycleStage::POST_SMOOTH_OUTPUT);
            expected_fac_no_pre_levels.push_back(ln);
            expected_fac_no_pre_levels.push_back(ln);
        }
        Pointer<SAMRAIVectorReal<NDIM, double>> fac_observer_sol;
        Pointer<SAMRAIVectorReal<NDIM, double>> fac_observer_rhs;
        bool fac_observer_first_cycle_valid = true;
        bool fac_observer_disabled_silent = true;
        bool fac_observer_no_pre_cycle_valid = true;
        if (verify_fac_cycle_observer)
        {
            const int configured_num_pre_sweeps = fac_pc->getNumPreSmoothingSweeps();
            fac_observer_sol = eul_sol_vec->cloneVector("fac_observer_sol");
            fac_observer_sol->allocateVectorData();
            fac_observer_sol->setToScalar(0.0);
            fac_observer_rhs = jv->cloneVector("fac_observer_rhs");
            fac_observer_rhs->allocateVectorData();
            fac_observer_rhs->copyVector(jv);
            fac_pc->setCycleObserver(fac_cycle_observer);
            fac_trace_exporter.active = verify_cav_live_dynamic_trace_schema;
            const bool observer_solve_success = fac_pc->solveSystem(*fac_observer_sol, *fac_observer_rhs);
            fac_trace_exporter.active = false;
            fac_observer_first_cycle_valid = observer_solve_success && observed_fac_stages == expected_fac_stages &&
                                             observed_fac_levels == expected_fac_levels && fac_observer_views_valid &&
                                             std::isfinite(fac_observer_sol->L2Norm()) &&
                                             fac_observer_sol->L2Norm() > 1.0e-14;
            if (!fac_observer_first_cycle_valid) ++test_failures;
            if (verify_cav_live_dynamic_trace_schema)
            {
                fac_trace_exporter.finalize(expected_fac_stages, expected_fac_levels);
                if (!fac_trace_exporter.artifacts_round_trip || !fac_trace_exporter.index_round_trip ||
                    !fac_trace_exporter.order_valid || !fac_trace_exporter.mutation_detected)
                    ++test_failures;
            }

            const std::size_t observed_stage_count = observed_fac_stages.size();
            fac_pc->setCycleObserver({});
            fac_observer_sol->setToScalar(0.0);
            fac_observer_rhs->copyVector(jv);
            const bool disabled_solve_success = fac_pc->solveSystem(*fac_observer_sol, *fac_observer_rhs);
            fac_observer_disabled_silent = disabled_solve_success && observed_fac_stages.size() == observed_stage_count;
            if (!fac_observer_disabled_silent) ++test_failures;

            observed_fac_stages.clear();
            observed_fac_levels.clear();
            fac_observer_views_valid = true;
            fac_pc->setNumPreSmoothingSweeps(0);
            fac_pc->setCycleObserver(fac_cycle_observer);
            fac_observer_sol->setToScalar(0.0);
            // The no-presmoothing FAC path restricts its hierarchy RHS in
            // place. Keep this observer-only check from changing the physical
            // RHS used by the subsequent FGMRES validation.
            fac_observer_rhs->copyVector(jv);
            const bool no_pre_solve_success = fac_pc->solveSystem(*fac_observer_sol, *fac_observer_rhs);
            fac_observer_no_pre_cycle_valid =
                no_pre_solve_success && observed_fac_stages == expected_fac_no_pre_stages &&
                observed_fac_levels == expected_fac_no_pre_levels && fac_observer_views_valid &&
                std::isfinite(fac_observer_sol->L2Norm()) && fac_observer_sol->L2Norm() > 1.0e-14;
            if (!fac_observer_no_pre_cycle_valid)
            {
                ++test_failures;
                pout << "fac_cycle_observer_no_pre_solve_success = "
                     << (no_pre_solve_success ? "true" : "false") << std::endl;
                pout << "fac_cycle_observer_no_pre_stage_order_valid = "
                     << (observed_fac_stages == expected_fac_no_pre_stages ? "true" : "false") << std::endl;
                pout << "fac_cycle_observer_no_pre_level_order_valid = "
                     << (observed_fac_levels == expected_fac_no_pre_levels ? "true" : "false") << std::endl;
                pout << "fac_cycle_observer_no_pre_views_valid = "
                     << (fac_observer_views_valid ? "true" : "false") << std::endl;
                pout << "fac_cycle_observer_no_pre_solution_norm = " << fac_observer_sol->L2Norm() << std::endl;
            }
            fac_pc->setCycleObserver({});
            fac_pc->setNumPreSmoothingSweeps(configured_num_pre_sweeps);
        }

        {
            PetscErrorCode ierr = 0;
            Mat saj_unscaled = nullptr;
            ierr = MatPtAP(A, J, MAT_INITIAL_MATRIX, 1.0, &saj_unscaled);
            IBTK_CHKERRQ(ierr);
            Pointer<PatchLevel<NDIM>> finest_level = patch_hierarchy->getPatchLevel(finest_ln);
            Pointer<CartesianGridGeometry<NDIM>> finest_grid_geom = patch_hierarchy->getGridGeometry();
            const double* const dx0 = finest_grid_geom->getDx();
            const IntVector<NDIM>& ratio = finest_level->getRatio();
            double cell_volume = 1.0;
            for (unsigned d = 0; d < NDIM; ++d)
            {
                cell_volume *= dx0[d] / static_cast<double>(ratio(d));
            }
            const double theta_ds = 2.0 * M_PI / static_cast<double>(structure_spec.num_curve_points);

            Mat saj_cell_scaled = nullptr;
            ierr = MatDuplicate(saj_unscaled, MAT_COPY_VALUES, &saj_cell_scaled);
            IBTK_CHKERRQ(ierr);
            ierr = MatScale(saj_cell_scaled, -dt / cell_volume);
            IBTK_CHKERRQ(ierr);

            Mat saj_theta_scaled = nullptr;
            ierr = MatDuplicate(saj_unscaled, MAT_COPY_VALUES, &saj_theta_scaled);
            IBTK_CHKERRQ(ierr);
            ierr = MatScale(saj_theta_scaled, -dt * theta_ds / cell_volume);
            IBTK_CHKERRQ(ierr);

            auto compute_rel_mat_error = [](Mat lhs, Mat rhs) -> double
            {
                Mat diff = nullptr;
                PetscErrorCode ierr_local = MatDuplicate(lhs, MAT_COPY_VALUES, &diff);
                IBTK_CHKERRQ(ierr_local);
                ierr_local = MatAXPY(diff, -1.0, rhs, DIFFERENT_NONZERO_PATTERN);
                IBTK_CHKERRQ(ierr_local);
                PetscReal diff_norm = 0.0;
                PetscReal rhs_norm = 0.0;
                ierr_local = MatNorm(diff, NORM_FROBENIUS, &diff_norm);
                IBTK_CHKERRQ(ierr_local);
                ierr_local = MatNorm(lhs, NORM_FROBENIUS, &rhs_norm);
                IBTK_CHKERRQ(ierr_local);
                ierr_local = MatDestroy(&diff);
                IBTK_CHKERRQ(ierr_local);
                return static_cast<double>(diff_norm) / std::max(static_cast<double>(rhs_norm), 1.0e-14);
            };

            const double saj_cell_scaled_rel_error = compute_rel_mat_error(SAJ, saj_cell_scaled);
            const double saj_theta_scaled_rel_error = compute_rel_mat_error(SAJ, saj_theta_scaled);
            pout << "saj_cell_scaled_relative_error = " << saj_cell_scaled_rel_error << std::endl;
            pout << "saj_theta_scaled_relative_error = " << saj_theta_scaled_rel_error << std::endl;

            ierr = MatDestroy(&saj_unscaled);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&saj_cell_scaled);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&saj_theta_scaled);
            IBTK_CHKERRQ(ierr);
        }

        const bool run_saj_vector_compare = input_db->getBoolWithDefault("RUN_SAJ_VECTOR_COMPARE", false);
        if (run_saj_vector_compare)
        {
            Pointer<SAMRAIVectorReal<NDIM, double>> saj_jv = eul_rhs_vec->cloneVector("saj_jv");
            saj_jv->allocateVectorData();
            saj_jv->setToScalar(0.0);
            jac_op->apply(*v, *saj_jv);

            Pointer<SAMRAIVectorReal<NDIM, double>> saj_diff = eul_rhs_vec->cloneVector("saj_diff");
            saj_diff->allocateVectorData();
            saj_diff->subtract(saj_jv, jv);

            double saj_jv_side_norm = std::numeric_limits<double>::quiet_NaN();
            double saj_jv_cell_norm = std::numeric_limits<double>::quiet_NaN();
            double saj_diff_side_norm = std::numeric_limits<double>::quiet_NaN();
            double saj_diff_cell_norm = std::numeric_limits<double>::quiet_NaN();
            const bool saj_jv_finite =
                side_l2_norm_is_finite(
                    hier_velocity_data_ops, saj_jv->getComponentDescriptorIndex(0), wgt_sc_idx, saj_jv_side_norm) &&
                cell_l2_norm_is_finite(
                    hier_pressure_data_ops, saj_jv->getComponentDescriptorIndex(1), wgt_cc_idx, saj_jv_cell_norm);
            const bool saj_diff_finite =
                side_l2_norm_is_finite(
                    hier_velocity_data_ops, saj_diff->getComponentDescriptorIndex(0), wgt_sc_idx, saj_diff_side_norm) &&
                cell_l2_norm_is_finite(
                    hier_pressure_data_ops, saj_diff->getComponentDescriptorIndex(1), wgt_cc_idx, saj_diff_cell_norm);
            if (!saj_jv_finite || !saj_diff_finite)
            {
                ++test_failures;
                pout << "saj jacobian comparison norms are non-finite" << std::endl;
            }
            else
            {
                const double saj_rel_error =
                    std::sqrt(saj_diff_side_norm * saj_diff_side_norm + saj_diff_cell_norm * saj_diff_cell_norm) /
                    std::max(std::sqrt(saj_jv_side_norm * saj_jv_side_norm + saj_jv_cell_norm * saj_jv_cell_norm),
                             1.0e-14);
                pout << "saj_relative_error = " << saj_rel_error << std::endl;
                const double saj_rel_tol = input_db->getDoubleWithDefault("SAJ_REL_TOL", 5.0e-12);
                if (!(saj_rel_error <= saj_rel_tol))
                {
                    ++test_failures;
                    pout << "saj_relative_error exceeds tolerance: " << saj_rel_tol << std::endl;
                }
            }
        }

        Pointer<PETScKrylovLinearSolver> linear_solver =
            new PETScKrylovLinearSolver("stokes_ib_solver_components::linear_solver", nullptr, "ib_");
        linear_solver->setOperator(jac_op);
        linear_solver->setPreconditioner(fac_pc);
        linear_solver->setTimeInterval(current_time, new_time);
        linear_solver->setSolutionTime(new_time);
        linear_solver->setInitialGuessNonzero(false);

        Pointer<SAMRAIVectorReal<NDIM, double>> linear_sol = eul_sol_vec->cloneVector("linear_sol");
        linear_sol->allocateVectorData();
        linear_sol->setToScalar(0.0);
        LocalSolveDiagnostics local_solve_diagnostics;
        KrylovTraceContext krylov_trace;
        bool krylov_trace_index_round_trip = !verify_cav_live_dynamic_trace_schema;
        bool krylov_trace_order_valid = !verify_cav_live_dynamic_trace_schema;
        bool krylov_trace_history_valid = !verify_cav_live_dynamic_trace_schema;
        bool krylov_trace_mutation_detected = !verify_cav_live_dynamic_trace_schema;
        bool physical_summary_round_trip = !verify_cav_live_dynamic_trace_schema;
        bool physical_summary_mutation_detected = !verify_cav_live_dynamic_trace_schema;
        if (verify_fgmres_physical_residuals)
        {
            linear_solver->initializeSolverState(*linear_sol, *jv);
            SAJ = fac_op->getEulerianElasticityLevelOp(finest_ln);
            PetscErrorCode ierr = KSPSetResidualHistory(linear_solver->getPETScKSP(), nullptr, 201, PETSC_TRUE);
            IBTK_CHKERRQ(ierr);
            if (verify_cav_live_dynamic_trace_schema)
            {
                krylov_trace.fac_operator = fac_op;
                krylov_trace.hierarchy = patch_hierarchy;
                krylov_trace.u_dof_index_idx = u_dof_index_idx;
                krylov_trace.p_dof_index_idx = p_dof_index_idx;
                ierr = KSPMonitorSet(linear_solver->getPETScKSP(), record_krylov_trace, &krylov_trace, nullptr);
                IBTK_CHKERRQ(ierr);
            }
        }
        if (verify_fgmres_local_diagnostics)
        {
            if (!verify_fgmres_physical_residuals || !pressure_cav_level_solver)
            {
                TBOX_ERROR("FGMRES local diagnostics require the pressure-CAV physical-residual case\n");
            }
            Mat level_operator = nullptr;
            PetscErrorCode ierr = KSPGetOperators(pressure_cav_level_solver->getPETScKSP(), &level_operator, nullptr);
            IBTK_CHKERRQ(ierr);
            local_solve_diagnostics.initialize(level_operator, pressure_cav_level_solver->getASMSubdomains());
            if (verify_cav_live_dynamic_trace_schema)
                local_solve_diagnostics.enableTrace(cav_local_trace_patch_ordinal);
            pressure_cav_level_solver->setShellSubdomainSolveObserver(
                [&](const int ordinal, Mat local_matrix, Vec local_rhs, Vec local_solution, Vec current_global_source) {
                    local_solve_diagnostics.observe(
                        ordinal, local_matrix, local_rhs, local_solution, current_global_source);
                });
        }
        const bool linear_success = linear_solver->solveSystem(*linear_sol, *jv);
        if (verify_cav_live_dynamic_trace_schema)
        {
            PetscErrorCode ierr = KSPMonitorCancel(linear_solver->getPETScKSP());
            IBTK_CHKERRQ(ierr);
            IBAMR::TestSupport::writeCAVKrylovTraceIndex(krylov_trace.records, "cav_fgmres_trace.txt");
            const auto reread_records = IBAMR::TestSupport::readCAVKrylovTraceIndex("cav_fgmres_trace.txt");
            krylov_trace_index_round_trip =
                IBAMR::TestSupport::sameCAVKrylovTraceIndex(krylov_trace.records, reread_records);
            krylov_trace_order_valid = !krylov_trace.records.empty();
            for (std::size_t k = 0; krylov_trace_order_valid && k < krylov_trace.records.size(); ++k)
                krylov_trace_order_valid = krylov_trace.records[k].iteration == static_cast<int>(k);
            auto mutated_records = reread_records;
            if (!mutated_records.empty()) ++mutated_records.front().iteration;
            krylov_trace_mutation_detected =
                !IBAMR::TestSupport::sameCAVKrylovTraceIndex(krylov_trace.records, mutated_records);
            if (!krylov_trace.artifacts_round_trip || !krylov_trace_index_round_trip || !krylov_trace_order_valid ||
                !krylov_trace_mutation_detected)
                ++test_failures;
        }
        if (verify_fgmres_local_diagnostics)
        {
            pressure_cav_level_solver->setShellSubdomainSolveObserver({});
            local_solve_diagnostics.finalizeTrace();
            const double local_diagnostic_tol = input_db->getDouble("FGMRES_LOCAL_DIAGNOSTIC_TOL");
            const bool local_diagnostics_valid =
                local_solve_diagnostics.complete() && local_solve_diagnostics.order_valid &&
                local_solve_diagnostics.values_finite &&
                local_solve_diagnostics.max_backward_error <= local_diagnostic_tol &&
                local_solve_diagnostics.max_incremental_fresh_error <= local_diagnostic_tol &&
                local_solve_diagnostics.max_local_rhs_error <= local_diagnostic_tol;
            const bool local_trace_valid =
                local_solve_diagnostics.trace_artifacts_round_trip && local_solve_diagnostics.trace_index_round_trip &&
                local_solve_diagnostics.trace_selection_valid && local_solve_diagnostics.trace_mutation_detected;
            if (!local_diagnostics_valid || (verify_cav_live_dynamic_trace_schema && !local_trace_valid))
                ++test_failures;
            pout << std::setprecision(16);
            pout << "fgmres_local_sweep_count = " << local_solve_diagnostics.sweep_count << std::endl;
            pout << "fgmres_local_solve_count = " << local_solve_diagnostics.solve_count << std::endl;
            pout << "fgmres_local_backward_error_max = " << local_solve_diagnostics.max_backward_error << std::endl;
            pout << "fgmres_incremental_fresh_error_max = " << local_solve_diagnostics.max_incremental_fresh_error
                 << std::endl;
            pout << "fgmres_local_rhs_error_max = " << local_solve_diagnostics.max_local_rhs_error << std::endl;
            pout << "fgmres_local_diagnostics_valid = " << (local_diagnostics_valid ? "true" : "false") << std::endl;
            if (verify_cav_live_dynamic_trace_schema)
            {
                pout << "cav_local_trace_record_count = " << local_solve_diagnostics.trace_records.size() << std::endl;
                pout << "cav_local_trace_artifacts_round_trip = "
                     << (local_solve_diagnostics.trace_artifacts_round_trip ? "true" : "false") << std::endl;
                pout << "cav_local_trace_index_round_trip = "
                     << (local_solve_diagnostics.trace_index_round_trip ? "true" : "false") << std::endl;
                pout << "cav_local_trace_selection_valid = "
                     << (local_solve_diagnostics.trace_selection_valid ? "true" : "false") << std::endl;
                pout << "cav_local_trace_mutation_detected = "
                     << (local_solve_diagnostics.trace_mutation_detected ? "true" : "false") << std::endl;
            }
            pout << std::setprecision(6);
            local_solve_diagnostics.deallocate();
        }
        if (!linear_success)
        {
            ++test_failures;
            pout << "krylov linear solve failed" << std::endl;
        }
        pout << "krylov_linear_iterations = " << linear_solver->getNumIterations() << std::endl;
        const bool report_krylov_linear_residual_norm = !input_db->keyExists("REPORT_KRYLOV_LINEAR_RESIDUAL_NORM") ||
                                                        input_db->getBool("REPORT_KRYLOV_LINEAR_RESIDUAL_NORM");
        if (report_krylov_linear_residual_norm)
        {
            pout << "krylov_linear_residual_norm = " << linear_solver->getResidualNorm() << std::endl;
        }
        else
        {
            const bool residual_norm_valid =
                std::isfinite(linear_solver->getResidualNorm()) && linear_solver->getResidualNorm() >= 0.0;
            if (!residual_norm_valid) ++test_failures;
            pout << "krylov_linear_residual_valid = " << (residual_norm_valid ? "true" : "false") << std::endl;
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

        if (verify_fgmres_physical_residuals)
        {
            const double original_relative_residual_tol = input_db->getDouble("FGMRES_ORIGINAL_RELATIVE_RESIDUAL_TOL");
            const double ib_coupling_relative_residual_tol = input_db->getDouble("IB_COUPLING_RELATIVE_RESIDUAL_TOL");
            const char* ksp_type = nullptr;
            PetscErrorCode ierr = KSPGetType(linear_solver->getPETScKSP(), &ksp_type);
            IBTK_CHKERRQ(ierr);
            const bool fgmres_type_valid = ksp_type && std::string(ksp_type) == KSPFGMRES;
            KSPConvergedReason converged_reason = KSP_CONVERGED_ITERATING;
            ierr = KSPGetConvergedReason(linear_solver->getPETScKSP(), &converged_reason);
            IBTK_CHKERRQ(ierr);
            const bool fgmres_converged = converged_reason > 0;

            const PetscReal* residual_history = nullptr;
            PetscInt residual_history_size = 0;
            ierr = KSPGetResidualHistory(linear_solver->getPETScKSP(), &residual_history, &residual_history_size);
            IBTK_CHKERRQ(ierr);
            bool fgmres_history_valid = residual_history_size >= 2;
            for (PetscInt k = 0; k < residual_history_size; ++k)
            {
                fgmres_history_valid = fgmres_history_valid && std::isfinite(residual_history[k]);
            }
            if (verify_cav_live_dynamic_trace_schema)
            {
                krylov_trace_history_valid =
                    krylov_trace.records.size() == static_cast<std::size_t>(residual_history_size);
                for (PetscInt k = 0; krylov_trace_history_valid && k < residual_history_size; ++k)
                {
                    krylov_trace_history_valid =
                        krylov_trace.records[static_cast<std::size_t>(k)].residual_norm == residual_history[k];
                }
            }
            const double fgmres_history_final_relative =
                residual_history_size >= 2 ?
                    residual_history[residual_history_size - 1] / std::max(residual_history[0], 1.0e-30) :
                    std::numeric_limits<double>::quiet_NaN();
            fgmres_history_valid = fgmres_history_valid && std::isfinite(fgmres_history_final_relative) &&
                                   fgmres_history_final_relative <= original_relative_residual_tol;

            Pointer<SAMRAIVectorReal<NDIM, double>> linear_action = eul_rhs_vec->cloneVector("linear_action");
            Pointer<SAMRAIVectorReal<NDIM, double>> original_residual = eul_rhs_vec->cloneVector("original_residual");
            Pointer<SAMRAIVectorReal<NDIM, double>> ib_action_matrix_free =
                eul_rhs_vec->cloneVector("ib_action_matrix_free");
            for (const Pointer<SAMRAIVectorReal<NDIM, double>>& diagnostic_vec :
                 { linear_action, original_residual, ib_action_matrix_free })
            {
                diagnostic_vec->allocateVectorData();
                diagnostic_vec->setToScalar(0.0);
            }

            // Re-evaluate the original matrix-free Jacobian after the solve;
            // the KSP-reported residual is not used as a physical substitute.
            jac_op->apply(*linear_sol, *linear_action);
            original_residual->subtract(jv, linear_action);

            // Evaluate the matrix-free IB action directly. With zero velocity
            // coefficients and zero pressure, the auxiliary Jacobian's
            // momentum component contains only the live IB coupling action;
            // its divergence component is intentionally discarded.
            PoissonSpecifications ib_only_coefs("stokes_ib_solver_components::ib_only_coefs");
            ib_only_coefs.setCConstant(0.0);
            ib_only_coefs.setDConstant(0.0);
            Pointer<StaggeredStokesOperator> ib_only_stokes_op =
                new StaggeredStokesOperator("stokes_ib_solver_components::ib_only_stokes_op", false);
            ib_only_stokes_op->setVelocityPoissonSpecifications(ib_only_coefs);
            ib_only_stokes_op->setPhysicalBcCoefs(u_bc_coefs, nullptr);
            ib_only_stokes_op->setTimeInterval(current_time, new_time);
            ib_only_stokes_op->setSolutionTime(new_time);
            StaggeredStokesIBOperator::Context ib_only_ctx = ctx;
            ib_only_ctx.stokes_op = ib_only_stokes_op;
            Pointer<StaggeredStokesIBJacobianOperator> ib_only_jac_op =
                new StaggeredStokesIBJacobianOperator("stokes_ib_solver_components::ib_only_jacobian_op");
            ib_only_jac_op->setOperatorContext(ib_only_ctx);
            ib_only_jac_op->setTimeInterval(current_time, new_time);
            ib_only_jac_op->setSolutionTime(new_time);
            ib_only_jac_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);
            ib_only_jac_op->formJacobian(*eul_sol_vec);
            Pointer<SAMRAIVectorReal<NDIM, double>> ib_input = linear_sol->cloneVector("ib_input");
            ib_input->allocateVectorData();
            ib_input->copyVector(linear_sol);
            hier_pressure_data_ops->setToScalar(ib_input->getComponentDescriptorIndex(1), 0.0, false);
            ib_only_jac_op->apply(*ib_input, *ib_action_matrix_free);
            hier_pressure_data_ops->setToScalar(ib_action_matrix_free->getComponentDescriptorIndex(1), 0.0, false);

            const double rhs_momentum_norm =
                hier_velocity_data_ops->L2Norm(jv->getComponentDescriptorIndex(0), wgt_sc_idx);
            const double rhs_divergence_norm =
                hier_pressure_data_ops->L2Norm(jv->getComponentDescriptorIndex(1), wgt_cc_idx);
            const double momentum_residual_norm =
                hier_velocity_data_ops->L2Norm(original_residual->getComponentDescriptorIndex(0), wgt_sc_idx);
            const double divergence_residual_norm =
                hier_pressure_data_ops->L2Norm(original_residual->getComponentDescriptorIndex(1), wgt_cc_idx);
            const double rhs_norm = std::hypot(rhs_momentum_norm, rhs_divergence_norm);
            const double original_residual_norm = std::hypot(momentum_residual_norm, divergence_residual_norm);
            const double denominator_floor = std::sqrt(std::numeric_limits<double>::epsilon()) * rhs_norm;
            const double original_relative_residual = original_residual_norm / std::max(rhs_norm, 1.0e-30);
            const double momentum_relative_residual =
                momentum_residual_norm / std::max({ rhs_momentum_norm, denominator_floor, 1.0e-30 });
            const double divergence_relative_residual =
                divergence_residual_norm / std::max({ rhs_divergence_norm, denominator_floor, 1.0e-30 });

            Vec saj_input = nullptr;
            Vec saj_output = nullptr;
            Vec rhs_output = nullptr;
            Vec ib_matrix_free_output = nullptr;
            ierr = MatCreateVecs(SAJ, &saj_input, &saj_output);
            IBTK_CHKERRQ(ierr);
            ierr = VecDuplicate(saj_output, &rhs_output);
            IBTK_CHKERRQ(ierr);
            ierr = VecDuplicate(saj_output, &ib_matrix_free_output);
            IBTK_CHKERRQ(ierr);
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(saj_input,
                                                                  linear_sol->getComponentDescriptorIndex(0),
                                                                  u_dof_index_idx,
                                                                  linear_sol->getComponentDescriptorIndex(1),
                                                                  p_dof_index_idx,
                                                                  patch_hierarchy->getPatchLevel(finest_ln));
            ierr = MatMult(SAJ, saj_input, saj_output);
            IBTK_CHKERRQ(ierr);
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(rhs_output,
                                                                  jv->getComponentDescriptorIndex(0),
                                                                  u_dof_index_idx,
                                                                  jv->getComponentDescriptorIndex(1),
                                                                  p_dof_index_idx,
                                                                  patch_hierarchy->getPatchLevel(finest_ln));
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(ib_matrix_free_output,
                                                                  ib_action_matrix_free->getComponentDescriptorIndex(0),
                                                                  u_dof_index_idx,
                                                                  ib_action_matrix_free->getComponentDescriptorIndex(1),
                                                                  p_dof_index_idx,
                                                                  patch_hierarchy->getPatchLevel(finest_ln));
            double ib_matrix_free_action_norm = std::numeric_limits<double>::quiet_NaN();
            double ib_saj_action_norm = std::numeric_limits<double>::quiet_NaN();
            double rhs_algebraic_norm = std::numeric_limits<double>::quiet_NaN();
            ierr = VecNorm(rhs_output, NORM_2, &rhs_algebraic_norm);
            IBTK_CHKERRQ(ierr);
            ierr = VecNorm(ib_matrix_free_output, NORM_2, &ib_matrix_free_action_norm);
            IBTK_CHKERRQ(ierr);
            ierr = VecNorm(saj_output, NORM_2, &ib_saj_action_norm);
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPY(ib_matrix_free_output, -1.0, saj_output);
            IBTK_CHKERRQ(ierr);
            double ib_coupling_residual_norm = std::numeric_limits<double>::quiet_NaN();
            ierr = VecNorm(ib_matrix_free_output, NORM_2, &ib_coupling_residual_norm);
            IBTK_CHKERRQ(ierr);
            const double ib_coupling_relative_residual =
                ib_coupling_residual_norm /
                std::max({ ib_matrix_free_action_norm,
                           ib_saj_action_norm,
                           std::sqrt(std::numeric_limits<double>::epsilon()) * rhs_algebraic_norm,
                           1.0e-30 });

            if (verify_cav_live_dynamic_trace_schema)
            {
                IBAMR::TestSupport::CAVPhysicalResidualSummary summary;
                summary.rhs_l2 = rhs_norm;
                summary.rhs_momentum_l2 = rhs_momentum_norm;
                summary.rhs_divergence_l2 = rhs_divergence_norm;
                summary.residual_l2 = original_residual_norm;
                summary.residual_momentum_l2 = momentum_residual_norm;
                summary.residual_divergence_l2 = divergence_residual_norm;
                summary.ib_matrix_free_action_l2 = ib_matrix_free_action_norm;
                summary.ib_matrix_action_l2 = ib_saj_action_norm;
                summary.ib_residual_l2 = ib_coupling_residual_norm;
                summary.denominator_floor = denominator_floor;
                summary.relative_residual = original_relative_residual;
                summary.momentum_relative_residual = momentum_relative_residual;
                summary.divergence_relative_residual = divergence_relative_residual;
                summary.ib_relative_residual = ib_coupling_relative_residual;
                IBAMR::TestSupport::writeCAVPhysicalResidualSummary(summary, "cav_physical_residual.txt");
                const auto reread_summary =
                    IBAMR::TestSupport::readCAVPhysicalResidualSummary("cav_physical_residual.txt");
                physical_summary_round_trip =
                    IBAMR::TestSupport::sameCAVPhysicalResidualSummary(summary, reread_summary);
                auto mutated_summary = reread_summary;
                mutated_summary.residual_momentum_l2 += 1.0;
                physical_summary_mutation_detected =
                    !IBAMR::TestSupport::sameCAVPhysicalResidualSummary(summary, mutated_summary);
            }

            const bool original_residual_valid = std::isfinite(original_relative_residual) &&
                                                 original_relative_residual <= original_relative_residual_tol;
            const bool ib_coupling_residual_valid = std::isfinite(ib_coupling_relative_residual) &&
                                                    ib_coupling_relative_residual <= ib_coupling_relative_residual_tol;
            if (!fgmres_type_valid || !fgmres_converged || !fgmres_history_valid || !original_residual_valid ||
                !ib_coupling_residual_valid ||
                (verify_cav_live_dynamic_trace_schema &&
                 (!krylov_trace_history_valid || !physical_summary_round_trip || !physical_summary_mutation_detected)))
            {
                ++test_failures;
            }

            pout << "fgmres_type_valid = " << (fgmres_type_valid ? "true" : "false") << std::endl;
            pout << "fgmres_converged_reason = " << static_cast<int>(converged_reason) << std::endl;
            pout << "fgmres_converged = " << (fgmres_converged ? "true" : "false") << std::endl;
            pout << "fgmres_residual_history_size = " << residual_history_size << std::endl;
            pout << std::setprecision(16);
            for (PetscInt k = 0; k < residual_history_size; ++k)
            {
                pout << "fgmres_residual_history_" << k << " = " << residual_history[k] << std::endl;
            }
            pout << "physical_stiffness = " << structure_spec.spring_stiffness * structure_spec.ds << std::endl;
            pout << "original_rhs_l2 = " << rhs_norm << std::endl;
            pout << "original_rhs_momentum_l2 = " << rhs_momentum_norm << std::endl;
            pout << "original_rhs_divergence_l2 = " << rhs_divergence_norm << std::endl;
            pout << "original_residual_l2 = " << original_residual_norm << std::endl;
            pout << "ib_matrix_free_action_l2 = " << ib_matrix_free_action_norm << std::endl;
            pout << "ib_saj_action_l2 = " << ib_saj_action_norm << std::endl;
            pout << std::setprecision(6);
            pout << "fgmres_history_final_relative = " << fgmres_history_final_relative << std::endl;
            pout << "fgmres_history_valid = " << (fgmres_history_valid ? "true" : "false") << std::endl;
            pout << "original_momentum_residual_l2 = " << momentum_residual_norm << std::endl;
            pout << "original_momentum_relative_residual = " << momentum_relative_residual << std::endl;
            pout << "original_divergence_residual_l2 = " << divergence_residual_norm << std::endl;
            pout << "original_divergence_relative_residual = " << divergence_relative_residual << std::endl;
            pout << "original_relative_residual = " << original_relative_residual << std::endl;
            pout << "ib_coupling_residual_l2 = " << ib_coupling_residual_norm << std::endl;
            pout << "ib_coupling_relative_residual = " << ib_coupling_relative_residual << std::endl;
            pout << "fgmres_original_residual_valid = " << (original_residual_valid ? "true" : "false") << std::endl;
            pout << "ib_coupling_residual_valid = " << (ib_coupling_residual_valid ? "true" : "false") << std::endl;
            if (verify_cav_live_dynamic_trace_schema)
            {
                pout << "cav_dynamic_level_operator_round_trip = "
                     << (dynamic_level_operator_round_trip ? "true" : "false") << std::endl;
                pout << "cav_dynamic_level_dof_map_round_trip = "
                     << (dynamic_level_dof_map_round_trip ? "true" : "false") << std::endl;
                pout << "cav_dynamic_level_cav_data_round_trip = "
                     << (dynamic_level_cav_data_round_trip ? "true" : "false") << std::endl;
                pout << "cav_fac_trace_record_count = " << fac_trace_exporter.records.size() << std::endl;
                pout << "cav_fac_trace_artifacts_round_trip = "
                     << (fac_trace_exporter.artifacts_round_trip ? "true" : "false") << std::endl;
                pout << "cav_fac_trace_index_round_trip = " << (fac_trace_exporter.index_round_trip ? "true" : "false")
                     << std::endl;
                pout << "cav_fac_trace_order_valid = " << (fac_trace_exporter.order_valid ? "true" : "false")
                     << std::endl;
                pout << "cav_fac_trace_mutation_detected = "
                     << (fac_trace_exporter.mutation_detected ? "true" : "false") << std::endl;
                pout << "cav_fgmres_trace_record_count = " << krylov_trace.records.size() << std::endl;
                pout << "cav_fgmres_trace_artifacts_round_trip = "
                     << (krylov_trace.artifacts_round_trip ? "true" : "false") << std::endl;
                pout << "cav_fgmres_trace_index_round_trip = " << (krylov_trace_index_round_trip ? "true" : "false")
                     << std::endl;
                pout << "cav_fgmres_trace_order_valid = " << (krylov_trace_order_valid ? "true" : "false") << std::endl;
                pout << "cav_fgmres_trace_history_valid = " << (krylov_trace_history_valid ? "true" : "false")
                     << std::endl;
                pout << "cav_fgmres_trace_mutation_detected = " << (krylov_trace_mutation_detected ? "true" : "false")
                     << std::endl;
                pout << "cav_physical_summary_round_trip = " << (physical_summary_round_trip ? "true" : "false")
                     << std::endl;
                pout << "cav_physical_summary_mutation_detected = "
                     << (physical_summary_mutation_detected ? "true" : "false") << std::endl;
            }

            ierr = VecDestroy(&ib_matrix_free_output);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&rhs_output);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&saj_output);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&saj_input);
            IBTK_CHKERRQ(ierr);
            ib_input->deallocateVectorData();
            ib_only_jac_op->deallocateOperatorState();
            for (const Pointer<SAMRAIVectorReal<NDIM, double>>& diagnostic_vec :
                 { linear_action, original_residual, ib_action_matrix_free })
            {
                diagnostic_vec->deallocateVectorData();
            }
            linear_solver->deallocateSolverState();
        }

        bool fac_observer_reinitialize_valid = true;
        bool galerkin_operator_borrowing_reinitialize_valid = !verify_galerkin_operator_borrowing;
        if (verify_fac_cycle_observer)
        {
            observed_fac_stages.clear();
            observed_fac_levels.clear();
            fac_observer_views_valid = true;
            fac_pc->setCycleObserver(fac_cycle_observer);
        }
        fac_pc->deallocateSolverState();
        if (verify_fac_cycle_observer)
        {
            fac_pc->initializeSolverState(*eul_sol_vec, *eul_rhs_vec);
            fac_observer_sol->setToScalar(0.0);
            fac_observer_rhs->copyVector(jv);
            const bool reinitialized_solve_success = fac_pc->solveSystem(*fac_observer_sol, *fac_observer_rhs);
            fac_observer_reinitialize_valid =
                reinitialized_solve_success && observed_fac_stages == expected_fac_stages &&
                observed_fac_levels == expected_fac_levels && fac_observer_views_valid &&
                std::isfinite(fac_observer_sol->L2Norm()) && fac_observer_sol->L2Norm() > 1.0e-14;
            if (!fac_observer_reinitialize_valid) ++test_failures;
            if (verify_galerkin_operator_borrowing)
            {
                galerkin_operator_borrowing_reinitialize_valid = galerkin_operator_borrowing_is_valid();
                if (!galerkin_operator_borrowing_reinitialize_valid) ++test_failures;
            }
            fac_pc->setCycleObserver({});
            fac_pc->deallocateSolverState();

            pout << "fac_cycle_observer_first_cycle_valid = " << (fac_observer_first_cycle_valid ? "true" : "false")
                 << std::endl;
            pout << "fac_cycle_observer_disabled_silent = " << (fac_observer_disabled_silent ? "true" : "false")
                 << std::endl;
            pout << "fac_cycle_observer_no_pre_cycle_valid = " << (fac_observer_no_pre_cycle_valid ? "true" : "false")
                 << std::endl;
            pout << "fac_cycle_observer_reinitialize_valid = " << (fac_observer_reinitialize_valid ? "true" : "false")
                 << std::endl;
        }
        if (verify_galerkin_operator_borrowing)
        {
            pout << "galerkin_operator_borrowing_valid = " << (galerkin_operator_borrowing_valid ? "true" : "false")
                 << std::endl;
            pout << "galerkin_operator_borrowing_reinitialize_valid = "
                 << (galerkin_operator_borrowing_reinitialize_valid ? "true" : "false") << std::endl;
        }
        if (verify_pressure_cav)
        {
            const bool state_cleared = pressure_cav_level_solver->getCouplingAwareASMPressureSeedDOFs().empty();
            if (!state_cleared) ++test_failures;
            pout << "pressure_cav_state_cleared = " << (state_cleared ? "true" : "false") << std::endl;
        }

        jac_op->deallocateOperatorState();
        nonlinear_op.deallocateOperatorState();

        ib_method_ops->postprocessIntegrateData(current_time, new_time, /*num_cycles*/ 1);

        deallocate_vector_data(*eul_sol_vec);
        deallocate_vector_data(*eul_rhs_vec);
        free_vector_components(*eul_sol_vec);
        free_vector_components(*eul_rhs_vec);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int data_idx : allocated_patch_data_indices)
            {
                if (level->checkAllocated(data_idx)) level->deallocatePatchData(data_idx);
            }
        }

        PetscErrorCode ierr = MatDestroy(&A);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&J);
        IBTK_CHKERRQ(ierr);

        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

        plog << "Input database:\n";
        input_db->printClassData(plog);
        pout << "test_failures = " << test_failures << std::endl;
    }

    return test_failures;
}
