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

#include <ibamr/INSAveragingTurbulenceStatistics.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/snapshot_utilities.h>

#include <tbox/MemoryDatabase.h>

#include <CellData.h>
#include <CellIndex.h>
#include <CellIterator.h>
#include <HierarchyDataOpsManager.h>
#include <NodeData.h>
#include <NodeIndex.h>
#include <NodeIterator.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RobinBcCoefStrategy.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cctype>
#include <utility>

#include <ibamr/app_namespaces.h> // IWYU pragma: keep

namespace
{
void
allocate_patch_data(const int idx,
                    const double time,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(idx)) level->allocatePatchData(idx, time);
    }
}

void
deallocate_patch_data(const int idx, const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(idx)) level->deallocatePatchData(idx);
    }
}

SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>
build_velocity_product_db(const std::string& object_name, const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db)
{
    auto db = new SAMRAI::tbox::MemoryDatabase(object_name + "::VelocityProductDB");
    db->putDouble("period_start", input_db->getDouble("period_start"));
    db->putDouble("period_end", input_db->getDouble("period_end"));
    db->putDouble("threshold", input_db->getDouble("threshold"));
    db->putInteger("num_snapshots", input_db->getInteger("num_snapshots"));
    db->putBool("enable_logging", input_db->getBool("enable_logging"));
    db->putBool("output_data", false);
    db->putString("dir_dump_name", input_db->getStringWithDefault("dir_dump_name", "./"));
    db->putString("refine_type", input_db->getStringWithDefault("refine_type", "CONSERVATIVE_LINEAR_REFINE"));
    return db;
}

template <class DataType, class IndexType>
void
store_symmetric_tensor(DataType& data, const IndexType& idx, const double* entries)
{
#if (NDIM == 2)
    data(idx, 0) = entries[0];
    data(idx, 1) = entries[1];
    data(idx, 2) = entries[2];
#endif
#if (NDIM == 3)
    data(idx, 0) = entries[0];
    data(idx, 1) = entries[1];
    data(idx, 2) = entries[2];
    data(idx, 3) = entries[3];
    data(idx, 4) = entries[4];
    data(idx, 5) = entries[5];
#endif
}

std::string
normalize_centering(std::string centering)
{
    std::transform(centering.begin(),
                   centering.end(),
                   centering.begin(),
                   [](const unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return centering;
}
} // namespace

namespace IBAMR
{
INSAveragingTurbulenceStatistics::INSAveragingTurbulenceStatistics(std::string object_name,
                                                                   Pointer<SideVariable<NDIM, double>> U_var,
                                                                   Pointer<Database> input_db,
                                                                   Pointer<GridGeometry<NDIM>> grid_geom,
                                                                   const bool register_for_restart)
    : INSTurbulenceStatistics(std::move(object_name), input_db->getDoubleWithDefault("statistics_start_time", 0.0))
{
    d_refine_type = input_db->getStringWithDefault("refine_type", d_refine_type);
    d_analysis_centering =
        normalize_centering(input_db->getStringWithDefault("analysis_centering", d_analysis_centering));
    if (d_analysis_centering != "CELL" && d_analysis_centering != "NODE")
    {
        TBOX_ERROR(d_object_name << "::INSAveragingTurbulenceStatistics(): unsupported analysis centering "
                                 << d_analysis_centering << "\n");
    }

    d_U_sc_scratch_var = new SideVariable<NDIM, double>(d_object_name + "::U_sc_scratch");
    d_U_cc_var = new CellVariable<NDIM, double>(d_object_name + "::U_cc", NDIM);
    d_U_nc_var = new NodeVariable<NDIM, double>(d_object_name + "::U_nc", NDIM, false);
    d_UU_cc_var = new CellVariable<NDIM, double>(d_object_name + "::UU_cc", NDIM * (NDIM + 1) / 2);
    d_UU_nc_var = new NodeVariable<NDIM, double>(d_object_name + "::UU_nc", NDIM * (NDIM + 1) / 2, false);
    d_U_mean_cc_var = new CellVariable<NDIM, double>(d_object_name + "::U_mean_cc", NDIM);
    d_U_mean_nc_var = new NodeVariable<NDIM, double>(d_object_name + "::U_mean_nc", NDIM, false);
    d_UU_mean_cc_var = new CellVariable<NDIM, double>(d_object_name + "::UU_mean_cc", NDIM * (NDIM + 1) / 2);
    d_UU_mean_nc_var = new NodeVariable<NDIM, double>(d_object_name + "::UU_mean_nc", NDIM * (NDIM + 1) / 2, false);

    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext(d_object_name + "::CONTEXT");
    d_U_sc_scratch_idx = var_db->registerVariableAndContext(d_U_sc_scratch_var, ctx, IntVector<NDIM>(1));
    d_U_cc_idx = var_db->registerVariableAndContext(d_U_cc_var, ctx, IntVector<NDIM>(0));
    d_U_nc_idx = var_db->registerVariableAndContext(d_U_nc_var, ctx, IntVector<NDIM>(0));
    d_UU_cc_idx = var_db->registerVariableAndContext(d_UU_cc_var, ctx, IntVector<NDIM>(0));
    d_UU_nc_idx = var_db->registerVariableAndContext(d_UU_nc_var, ctx, IntVector<NDIM>(0));
    d_U_mean_cc_idx = var_db->registerVariableAndContext(d_U_mean_cc_var, ctx, IntVector<NDIM>(0));
    d_U_mean_nc_idx = var_db->registerVariableAndContext(d_U_mean_nc_var, ctx, IntVector<NDIM>(0));
    d_UU_mean_cc_idx = var_db->registerVariableAndContext(d_UU_mean_cc_var, ctx, IntVector<NDIM>(0));
    d_UU_mean_nc_idx = var_db->registerVariableAndContext(d_UU_mean_nc_var, ctx, IntVector<NDIM>(0));

    if (d_analysis_centering == "CELL")
    {
        d_U_avg_manager = new IBTK::HierarchyAveragedDataManager(
            d_object_name + "::VelocityAveraging", d_U_cc_var, input_db, grid_geom, register_for_restart);
        d_UU_avg_manager = new IBTK::HierarchyAveragedDataManager(d_object_name + "::VelocityProductAveraging",
                                                                  d_UU_cc_var,
                                                                  build_velocity_product_db(d_object_name, input_db),
                                                                  grid_geom,
                                                                  register_for_restart);
    }
    else
    {
        d_U_avg_manager = new IBTK::HierarchyAveragedDataManager(
            d_object_name + "::VelocityAveraging", d_U_nc_var, input_db, grid_geom, register_for_restart);
        d_UU_avg_manager = new IBTK::HierarchyAveragedDataManager(d_object_name + "::VelocityProductAveraging",
                                                                  d_UU_nc_var,
                                                                  build_velocity_product_db(d_object_name, input_db),
                                                                  grid_geom,
                                                                  register_for_restart);
    }
}

bool
INSAveragingTurbulenceStatistics::updateStatistics(const int U_idx,
                                                   const Pointer<SideVariable<NDIM, double>> U_var,
                                                   const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                                                   const double data_time,
                                                   const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                   const Pointer<IBTK::HierarchyMathOps> hier_math_ops)
{
    // Delay statistics accumulation until the configured start time.
    if (!shouldUpdateStatistics(data_time)) return false;

    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
    const int U_avg_idx = d_analysis_centering == "CELL" ? d_U_cc_idx : d_U_nc_idx;
    const int UU_avg_idx = d_analysis_centering == "CELL" ? d_UU_cc_idx : d_UU_nc_idx;
    const int wgt_idx = d_analysis_centering == "CELL" ? wgt_cc_idx : IBTK::invalid_index;

    allocate_patch_data(d_U_sc_scratch_idx, data_time, hierarchy);
    allocate_patch_data(U_avg_idx, data_time, hierarchy);
    allocate_patch_data(UU_avg_idx, data_time, hierarchy);

    Pointer<PatchHierarchy<NDIM>> hierarchy_nc = hierarchy;
    auto* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double>> hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(d_U_sc_scratch_var, hierarchy_nc, true);
    hier_sc_data_ops->copyData(d_U_sc_scratch_idx, U_idx);

    using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    IBTK::HierarchyGhostCellInterpolation ghost_fill;
    ghost_fill.initializeOperatorState(ITC(d_U_sc_scratch_idx,
                                           "CONSERVATIVE_LINEAR_REFINE",
                                           true,
                                           "CONSERVATIVE_COARSEN",
                                           "LINEAR",
                                           false,
                                           velocity_bc_coefs),
                                       hierarchy);
    ghost_fill.fillData(data_time);

    if (d_analysis_centering == "CELL")
    {
        hier_math_ops->interp(d_U_cc_idx, d_U_cc_var, d_U_sc_scratch_idx, d_U_sc_scratch_var, nullptr, data_time, true);
    }
    else
    {
        hier_math_ops->interp(
            d_U_nc_idx, d_U_nc_var, true, d_U_sc_scratch_idx, d_U_sc_scratch_var, nullptr, data_time, true);
    }

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            if (d_analysis_centering == "CELL")
            {
                Pointer<CellData<NDIM, double>> U_cc_data = patch->getPatchData(d_U_cc_idx);
                Pointer<CellData<NDIM, double>> UU_data = patch->getPatchData(d_UU_cc_idx);

                for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    const CellIndex<NDIM> idx = *ci;
#if (NDIM == 2)
                    const double uu[3] = { (*U_cc_data)(idx, 0) * (*U_cc_data)(idx, 0),
                                           (*U_cc_data)(idx, 1) * (*U_cc_data)(idx, 1),
                                           (*U_cc_data)(idx, 0) * (*U_cc_data)(idx, 1) };
#endif
#if (NDIM == 3)
                    const double uu[6] = {
                        (*U_cc_data)(idx, 0) * (*U_cc_data)(idx, 0), (*U_cc_data)(idx, 1) * (*U_cc_data)(idx, 1),
                        (*U_cc_data)(idx, 2) * (*U_cc_data)(idx, 2), (*U_cc_data)(idx, 1) * (*U_cc_data)(idx, 2),
                        (*U_cc_data)(idx, 0) * (*U_cc_data)(idx, 2), (*U_cc_data)(idx, 0) * (*U_cc_data)(idx, 1)
                    };
#endif
                    store_symmetric_tensor(*UU_data, idx, uu);
                }
            }
            else
            {
                Pointer<NodeData<NDIM, double>> U_nc_data = patch->getPatchData(d_U_nc_idx);
                Pointer<NodeData<NDIM, double>> UU_data = patch->getPatchData(d_UU_nc_idx);

                for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
                {
                    const NodeIndex<NDIM> idx = *ni;
#if (NDIM == 2)
                    const double uu[3] = { (*U_nc_data)(idx, 0) * (*U_nc_data)(idx, 0),
                                           (*U_nc_data)(idx, 1) * (*U_nc_data)(idx, 1),
                                           (*U_nc_data)(idx, 0) * (*U_nc_data)(idx, 1) };
#endif
#if (NDIM == 3)
                    const double uu[6] = {
                        (*U_nc_data)(idx, 0) * (*U_nc_data)(idx, 0), (*U_nc_data)(idx, 1) * (*U_nc_data)(idx, 1),
                        (*U_nc_data)(idx, 2) * (*U_nc_data)(idx, 2), (*U_nc_data)(idx, 1) * (*U_nc_data)(idx, 2),
                        (*U_nc_data)(idx, 0) * (*U_nc_data)(idx, 2), (*U_nc_data)(idx, 0) * (*U_nc_data)(idx, 1)
                    };
#endif
                    store_symmetric_tensor(*UU_data, idx, uu);
                }
            }
        }
    }

    const bool U_steady = d_U_avg_manager->updateTimeAveragedSnapshot(U_avg_idx, data_time, hierarchy, wgt_idx);
    const bool UU_steady = d_UU_avg_manager->updateTimeAveragedSnapshot(UU_avg_idx, data_time, hierarchy, wgt_idx);

    deallocate_patch_data(UU_avg_idx, hierarchy);
    deallocate_patch_data(U_avg_idx, hierarchy);
    deallocate_patch_data(d_U_sc_scratch_idx, hierarchy);
    return U_steady && UU_steady;
}

bool
INSAveragingTurbulenceStatistics::isAtSteadyState() const
{
    return d_U_avg_manager->isAtPeriodicSteadyState() && d_UU_avg_manager->isAtPeriodicSteadyState();
}

const std::string&
INSAveragingTurbulenceStatistics::getAnalysisCentering() const
{
    return d_analysis_centering;
}

IBTK::HierarchyAveragedDataManager&
INSAveragingTurbulenceStatistics::getAveragedVelocityManager()
{
    return *d_U_avg_manager;
}

const IBTK::HierarchyAveragedDataManager&
INSAveragingTurbulenceStatistics::getAveragedVelocityManager() const
{
    return *d_U_avg_manager;
}

IBTK::HierarchyAveragedDataManager&
INSAveragingTurbulenceStatistics::getAveragedVelocityProductManager()
{
    return *d_UU_avg_manager;
}

const IBTK::HierarchyAveragedDataManager&
INSAveragingTurbulenceStatistics::getAveragedVelocityProductManager() const
{
    return *d_UU_avg_manager;
}

void
INSAveragingTurbulenceStatistics::fillReynoldsStressSnapshot(const int R_idx,
                                                             const Pointer<CellVariable<NDIM, double>> /*R_var*/,
                                                             const double time,
                                                             const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                             const Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                             const double tol) const
{
    if (d_analysis_centering != "CELL")
    {
        TBOX_ERROR(d_object_name << "::fillReynoldsStressSnapshot(): cell-centered Reynolds stresses requested, but "
                                 << "analysis_centering = " << d_analysis_centering << "\n");
    }

    const double snapshot_time = d_U_avg_manager->getTimePoint(time, tol);
    allocate_patch_data(d_U_mean_cc_idx, snapshot_time, hierarchy);
    allocate_patch_data(d_UU_mean_cc_idx, snapshot_time, hierarchy);

    IBTK::fill_snapshot_on_hierarchy(
        d_U_avg_manager->getSnapshotCache(), d_U_mean_cc_idx, snapshot_time, hierarchy, d_refine_type, tol);
    IBTK::fill_snapshot_on_hierarchy(
        d_UU_avg_manager->getSnapshotCache(), d_UU_mean_cc_idx, snapshot_time, hierarchy, d_refine_type, tol);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<CellData<NDIM, double>> U_mean_data = patch->getPatchData(d_U_mean_cc_idx);
            Pointer<CellData<NDIM, double>> UU_mean_data = patch->getPatchData(d_UU_mean_cc_idx);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM> idx = *ci;
                const int depth = R_data->getDepth();
#if (NDIM == 2)
                const double r_sym[3] = { (*UU_mean_data)(idx, 0) - (*U_mean_data)(idx, 0) * (*U_mean_data)(idx, 0),
                                          (*UU_mean_data)(idx, 1) - (*U_mean_data)(idx, 1) * (*U_mean_data)(idx, 1),
                                          (*UU_mean_data)(idx, 2) - (*U_mean_data)(idx, 0) * (*U_mean_data)(idx, 1) };
                if (depth == 3)
                {
                    store_symmetric_tensor(*R_data, idx, r_sym);
                }
                else if (depth == 4)
                {
                    (*R_data)(idx, 0) = r_sym[0];
                    (*R_data)(idx, 1) = r_sym[2];
                    (*R_data)(idx, 2) = r_sym[2];
                    (*R_data)(idx, 3) = r_sym[1];
                }
                else
                {
                    TBOX_ERROR(d_object_name << "::fillReynoldsStressSnapshot(): unsupported 2D tensor depth " << depth
                                             << "\n");
                }
#endif
#if (NDIM == 3)
                const double r_sym[6] = { (*UU_mean_data)(idx, 0) - (*U_mean_data)(idx, 0) * (*U_mean_data)(idx, 0),
                                          (*UU_mean_data)(idx, 1) - (*U_mean_data)(idx, 1) * (*U_mean_data)(idx, 1),
                                          (*UU_mean_data)(idx, 2) - (*U_mean_data)(idx, 2) * (*U_mean_data)(idx, 2),
                                          (*UU_mean_data)(idx, 3) - (*U_mean_data)(idx, 1) * (*U_mean_data)(idx, 2),
                                          (*UU_mean_data)(idx, 4) - (*U_mean_data)(idx, 0) * (*U_mean_data)(idx, 2),
                                          (*UU_mean_data)(idx, 5) - (*U_mean_data)(idx, 0) * (*U_mean_data)(idx, 1) };
                if (depth == 6)
                {
                    store_symmetric_tensor(*R_data, idx, r_sym);
                }
                else if (depth == 9)
                {
                    (*R_data)(idx, 0) = r_sym[0];
                    (*R_data)(idx, 1) = r_sym[5];
                    (*R_data)(idx, 2) = r_sym[4];
                    (*R_data)(idx, 3) = r_sym[5];
                    (*R_data)(idx, 4) = r_sym[1];
                    (*R_data)(idx, 5) = r_sym[3];
                    (*R_data)(idx, 6) = r_sym[4];
                    (*R_data)(idx, 7) = r_sym[3];
                    (*R_data)(idx, 8) = r_sym[2];
                }
                else
                {
                    TBOX_ERROR(d_object_name << "::fillReynoldsStressSnapshot(): unsupported 3D tensor depth " << depth
                                             << "\n");
                }
#endif
            }
        }
    }

    deallocate_patch_data(d_UU_mean_cc_idx, hierarchy);
    deallocate_patch_data(d_U_mean_cc_idx, hierarchy);
}

void
INSAveragingTurbulenceStatistics::fillReynoldsStressSnapshot(const int R_idx,
                                                             const Pointer<NodeVariable<NDIM, double>> /*R_var*/,
                                                             const double time,
                                                             const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                             const Pointer<IBTK::HierarchyMathOps> /*hier_math_ops*/,
                                                             const double tol) const
{
    if (d_analysis_centering != "NODE")
    {
        TBOX_ERROR(d_object_name << "::fillReynoldsStressSnapshot(): node-centered Reynolds stresses requested, but "
                                 << "analysis_centering = " << d_analysis_centering << "\n");
    }

    const double snapshot_time = d_U_avg_manager->getTimePoint(time, tol);
    allocate_patch_data(d_U_mean_nc_idx, snapshot_time, hierarchy);
    allocate_patch_data(d_UU_mean_nc_idx, snapshot_time, hierarchy);

    IBTK::fill_snapshot_on_hierarchy(
        d_U_avg_manager->getSnapshotCache(), d_U_mean_nc_idx, snapshot_time, hierarchy, d_refine_type, tol);
    IBTK::fill_snapshot_on_hierarchy(
        d_UU_avg_manager->getSnapshotCache(), d_UU_mean_nc_idx, snapshot_time, hierarchy, d_refine_type, tol);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<NodeData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<NodeData<NDIM, double>> U_mean_data = patch->getPatchData(d_U_mean_nc_idx);
            Pointer<NodeData<NDIM, double>> UU_mean_data = patch->getPatchData(d_UU_mean_nc_idx);

            for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
            {
                const NodeIndex<NDIM> idx = *ni;
                const int depth = R_data->getDepth();
#if (NDIM == 2)
                const double r_sym[3] = { (*UU_mean_data)(idx, 0) - (*U_mean_data)(idx, 0) * (*U_mean_data)(idx, 0),
                                          (*UU_mean_data)(idx, 1) - (*U_mean_data)(idx, 1) * (*U_mean_data)(idx, 1),
                                          (*UU_mean_data)(idx, 2) - (*U_mean_data)(idx, 0) * (*U_mean_data)(idx, 1) };
                if (depth == 3)
                {
                    store_symmetric_tensor(*R_data, idx, r_sym);
                }
                else if (depth == 4)
                {
                    (*R_data)(idx, 0) = r_sym[0];
                    (*R_data)(idx, 1) = r_sym[2];
                    (*R_data)(idx, 2) = r_sym[2];
                    (*R_data)(idx, 3) = r_sym[1];
                }
                else
                {
                    TBOX_ERROR(d_object_name << "::fillReynoldsStressSnapshot(): unsupported 2D tensor depth " << depth
                                             << "\n");
                }
#endif
#if (NDIM == 3)
                const double r_sym[6] = { (*UU_mean_data)(idx, 0) - (*U_mean_data)(idx, 0) * (*U_mean_data)(idx, 0),
                                          (*UU_mean_data)(idx, 1) - (*U_mean_data)(idx, 1) * (*U_mean_data)(idx, 1),
                                          (*UU_mean_data)(idx, 2) - (*U_mean_data)(idx, 2) * (*U_mean_data)(idx, 2),
                                          (*UU_mean_data)(idx, 3) - (*U_mean_data)(idx, 1) * (*U_mean_data)(idx, 2),
                                          (*UU_mean_data)(idx, 4) - (*U_mean_data)(idx, 0) * (*U_mean_data)(idx, 2),
                                          (*UU_mean_data)(idx, 5) - (*U_mean_data)(idx, 0) * (*U_mean_data)(idx, 1) };
                if (depth == 6)
                {
                    store_symmetric_tensor(*R_data, idx, r_sym);
                }
                else if (depth == 9)
                {
                    (*R_data)(idx, 0) = r_sym[0];
                    (*R_data)(idx, 1) = r_sym[5];
                    (*R_data)(idx, 2) = r_sym[4];
                    (*R_data)(idx, 3) = r_sym[5];
                    (*R_data)(idx, 4) = r_sym[1];
                    (*R_data)(idx, 5) = r_sym[3];
                    (*R_data)(idx, 6) = r_sym[4];
                    (*R_data)(idx, 7) = r_sym[3];
                    (*R_data)(idx, 8) = r_sym[2];
                }
                else
                {
                    TBOX_ERROR(d_object_name << "::fillReynoldsStressSnapshot(): unsupported 3D tensor depth " << depth
                                             << "\n");
                }
#endif
            }
        }
    }

    deallocate_patch_data(d_UU_mean_nc_idx, hierarchy);
    deallocate_patch_data(d_U_mean_nc_idx, hierarchy);
}
} // namespace IBAMR
