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
#include <VisItDataWriter.h>

#include <array>
#include <cmath>
#include <utility>

#include <ibamr/app_namespaces.h> // IWYU pragma: keep

namespace
{
constexpr int SYM_TENSOR_DEPTH = NDIM * (NDIM + 1) / 2;

#if (NDIM == 2)
constexpr std::array<std::array<int, 2>, SYM_TENSOR_DEPTH> SYM_COMPONENTS = { { { 0, 0 }, { 1, 1 }, { 0, 1 } } };
#endif

#if (NDIM == 3)
constexpr std::array<std::array<int, 2>, SYM_TENSOR_DEPTH> SYM_COMPONENTS = {
    { { 0, 0 }, { 1, 1 }, { 2, 2 }, { 1, 2 }, { 0, 2 }, { 0, 1 } }
};
#endif

int
sym_component_index(const int comp0, const int comp1)
{
    const int i = comp0 <= comp1 ? comp0 : comp1;
    const int j = comp0 <= comp1 ? comp1 : comp0;
    for (int comp = 0; comp < SYM_TENSOR_DEPTH; ++comp)
    {
        if (SYM_COMPONENTS[comp][0] == i && SYM_COMPONENTS[comp][1] == j) return comp;
    }
    TBOX_ERROR("sym_component_index(): unsupported component pair (" << i << ", " << j << ")\n");
    return IBTK::invalid_index;
}

template <class DataType, class IndexType>
void
store_symmetric_tensor(DataType& data, const IndexType& idx, const std::array<double, SYM_TENSOR_DEPTH>& entries)
{
    for (int comp = 0; comp < SYM_TENSOR_DEPTH; ++comp) data(idx, comp) = entries[comp];
}

template <class DataType, class IndexType>
std::array<double, NDIM>
read_velocity(const DataType& data, const IndexType& idx)
{
    std::array<double, NDIM> velocity = {};
    for (int comp = 0; comp < NDIM; ++comp) velocity[comp] = data(idx, comp);
    return velocity;
}

std::array<double, SYM_TENSOR_DEPTH>
compute_symmetric_second_moment(const std::array<double, NDIM>& velocity)
{
    std::array<double, SYM_TENSOR_DEPTH> second_moment = {};
    for (int comp = 0; comp < SYM_TENSOR_DEPTH; ++comp)
    {
        const int i = SYM_COMPONENTS[comp][0];
        const int j = SYM_COMPONENTS[comp][1];
        second_moment[comp] = velocity[i] * velocity[j];
    }
    return second_moment;
}

template <class DataType, class IndexType>
void
fill_symmetric_second_moment(DataType& velocity_product_data, const DataType& velocity_data, const IndexType& idx)
{
    const auto velocity = read_velocity(velocity_data, idx);
    const auto second_moment = compute_symmetric_second_moment(velocity);
    store_symmetric_tensor(velocity_product_data, idx, second_moment);
}

template <class ReynoldsDataType, class MeanDataType, class IndexType>
void
fill_reynolds_tensor(ReynoldsDataType& reynolds_data,
                     const MeanDataType& mean_velocity_data,
                     const MeanDataType& mean_velocity_product_data,
                     const IndexType& idx,
                     const std::string& object_name)
{
    const auto mean_velocity = read_velocity(mean_velocity_data, idx);
    std::array<double, SYM_TENSOR_DEPTH> reynolds_sym = {};
    for (int comp = 0; comp < SYM_TENSOR_DEPTH; ++comp)
    {
        reynolds_sym[comp] = mean_velocity_product_data(idx, comp);
    }
    for (int comp = 0; comp < SYM_TENSOR_DEPTH; ++comp)
    {
        const int i = SYM_COMPONENTS[comp][0];
        const int j = SYM_COMPONENTS[comp][1];
        reynolds_sym[comp] -= mean_velocity[i] * mean_velocity[j];
    }

    const int depth = reynolds_data.getDepth();
    if (depth == SYM_TENSOR_DEPTH)
    {
        store_symmetric_tensor(reynolds_data, idx, reynolds_sym);
    }
    else if (depth == NDIM * NDIM)
    {
        for (int comp_i = 0; comp_i < NDIM; ++comp_i)
        {
            for (int comp_j = 0; comp_j < NDIM; ++comp_j)
            {
                const int component = sym_component_index(comp_i, comp_j);
                reynolds_data(idx, comp_i * NDIM + comp_j) = reynolds_sym[component];
            }
        }
    }
    else
    {
        TBOX_ERROR(object_name << ": unsupported Reynolds-tensor depth " << depth << "\n");
    }
}

void
copy_averaging_manager_db_entries(const Pointer<Database>& src_db, const Pointer<Database>& dst_db)
{
    if (!src_db) return;

    if (src_db->keyExists("period_start")) dst_db->putDouble("period_start", src_db->getDouble("period_start"));
    if (src_db->keyExists("period_end")) dst_db->putDouble("period_end", src_db->getDouble("period_end"));
    if (src_db->keyExists("threshold")) dst_db->putDouble("threshold", src_db->getDouble("threshold"));
    if (src_db->keyExists("num_snapshots")) dst_db->putInteger("num_snapshots", src_db->getInteger("num_snapshots"));
    if (src_db->keyExists("enable_logging")) dst_db->putBool("enable_logging", src_db->getBool("enable_logging"));
    if (src_db->keyExists("output_data")) dst_db->putBool("output_data", src_db->getBool("output_data"));
    if (src_db->keyExists("dir_dump_name")) dst_db->putString("dir_dump_name", src_db->getString("dir_dump_name"));
    if (src_db->keyExists("refine_type")) dst_db->putString("refine_type", src_db->getString("refine_type"));
    if (src_db->keyExists("gcw")) dst_db->putIntegerArray("gcw", src_db->getIntegerArray("gcw"));
}

Pointer<Database>
build_averaging_manager_db(const std::string& object_name,
                           const Pointer<Database>& input_db,
                           const std::string& sub_db_name)
{
    Pointer<Database> manager_db = new MemoryDatabase(object_name);
    copy_averaging_manager_db_entries(input_db, manager_db);
    if (input_db->isDatabase(sub_db_name))
    {
        copy_averaging_manager_db_entries(input_db->getDatabase(sub_db_name), manager_db);
    }
    return manager_db;
}

double
map_to_period(const double t_start, const double t_end, double time)
{
    const double period = t_end - t_start;
    if (IBTK::abs_equal_eps(period, 0.0)) return t_start;
#ifndef NDEBUG
    TBOX_ASSERT(period > 0.0);
#endif
    time -= t_start;
    time = std::abs(std::fmod(time, period));
    time += t_start;
    return time;
}

template <class VariableType>
void
fill_snapshot_at_requested_time(IBTK::SnapshotCache& snapshot_cache,
                                const int dst_idx,
                                const int scratch_idx,
                                const Pointer<VariableType>& var,
                                double time,
                                const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                const std::string& refine_type,
                                const double period_length,
                                const double tol)
{
    const auto snapshot = snapshot_cache.getSnapshot(time, tol);
    if (snapshot.second)
    {
        IBTK::fill_snapshot_on_hierarchy(snapshot_cache, dst_idx, snapshot.first, hierarchy, refine_type, tol);
        return;
    }

    Pointer<PatchHierarchy<NDIM>> hierarchy_nc = hierarchy;
    Pointer<HierarchyDataOpsReal<NDIM, double>> hier_data_ops =
        HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(var, hierarchy_nc, true);
    hier_data_ops->resetLevels(0, hierarchy_nc->getFinestLevelNumber());
    IBTK::fill_snapshot_at_time(
        snapshot_cache, dst_idx, time, scratch_idx, hierarchy_nc, refine_type, hier_data_ops, period_length);
}
} // namespace

namespace IBAMR
{
INSAveragingTurbulenceStatistics::INSAveragingTurbulenceStatistics(std::string object_name,
                                                                   Pointer<SideVariable<NDIM, double>> /*U_var*/,
                                                                   Pointer<Database> input_db,
                                                                   Pointer<GridGeometry<NDIM>> grid_geom,
                                                                   const bool register_for_restart)
    : INSTurbulenceStatistics(std::move(object_name), input_db->getDoubleWithDefault("statistics_start_time", 0.0))
{
    d_refine_type = input_db->getStringWithDefault("refine_type", d_refine_type);
    d_data_centering = string_to_enum<DataCentering>(input_db->getStringWithDefault(
        "analysis_centering", IBAMR::enum_to_string<IBAMR::DataCentering>(d_data_centering)));
    d_output_reynolds_stress = input_db->getBoolWithDefault("output_reynolds_stress", d_output_reynolds_stress);
    d_output_tke = input_db->getBoolWithDefault("output_tke", d_output_tke);
    d_reynolds_stress_output_name =
        input_db->getStringWithDefault("reynolds_stress_output_name", d_reynolds_stress_output_name);
    d_tke_output_name = input_db->getStringWithDefault("tke_output_name", d_tke_output_name);

    Pointer<Database> velocity_averaging_db =
        build_averaging_manager_db(d_object_name + "::VelocityAveragingDB", input_db, "VelocityAveraging");
    Pointer<Database> velocity_product_averaging_db = build_averaging_manager_db(
        d_object_name + "::VelocityProductAveragingDB", input_db, "VelocityProductAveraging");
    if (!velocity_product_averaging_db->keyExists("output_data") ||
        !(input_db->isDatabase("VelocityProductAveraging") &&
          input_db->getDatabase("VelocityProductAveraging")->keyExists("output_data")))
    {
        velocity_product_averaging_db->putBool("output_data", false);
    }

    d_period_start = velocity_averaging_db->getDouble("period_start");
    d_period_end = velocity_averaging_db->getDouble("period_end");
    d_period_length = d_period_end - d_period_start;

    d_velocity_side_scratch_var = new SideVariable<NDIM, double>(d_object_name + "::U_sc_scratch");
    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext(d_object_name + "::CONTEXT");
    Pointer<VariableContext> plot_ctx = var_db->getContext(d_object_name + "::PLOT_CONTEXT");
    d_velocity_side_scratch_idx =
        var_db->registerVariableAndContext(d_velocity_side_scratch_var, ctx, IntVector<NDIM>(1));

    switch (d_data_centering)
    {
    case DataCentering::CELL:
        d_velocity_cell_var = new CellVariable<NDIM, double>(d_object_name + "::U_cc", NDIM);
        d_velocity_product_cell_var = new CellVariable<NDIM, double>(d_object_name + "::UU_cc", SYM_TENSOR_DEPTH);
        d_velocity_mean_cell_var = new CellVariable<NDIM, double>(d_object_name + "::U_mean_cc", NDIM);
        d_velocity_product_mean_cell_var =
            new CellVariable<NDIM, double>(d_object_name + "::UU_mean_cc", SYM_TENSOR_DEPTH);
        d_velocity_cell_idx = var_db->registerVariableAndContext(d_velocity_cell_var, ctx, IntVector<NDIM>(0));
        d_velocity_product_cell_idx =
            var_db->registerVariableAndContext(d_velocity_product_cell_var, ctx, IntVector<NDIM>(0));
        d_velocity_mean_cell_idx =
            var_db->registerVariableAndContext(d_velocity_mean_cell_var, ctx, IntVector<NDIM>(0));
        d_velocity_product_mean_cell_idx =
            var_db->registerVariableAndContext(d_velocity_product_mean_cell_var, ctx, IntVector<NDIM>(0));
        if (d_output_reynolds_stress)
        {
            d_reynolds_stress_cell_var =
                new CellVariable<NDIM, double>(d_object_name + "::ReynoldsStress_cc", NDIM * NDIM);
            d_reynolds_stress_plot_idx =
                var_db->registerVariableAndContext(d_reynolds_stress_cell_var, plot_ctx, IntVector<NDIM>(0));
        }
        if (d_output_tke)
        {
            d_tke_cell_var = new CellVariable<NDIM, double>(d_object_name + "::TKE_cc");
            d_tke_plot_idx = var_db->registerVariableAndContext(d_tke_cell_var, plot_ctx, IntVector<NDIM>(0));
        }
        d_velocity_average_manager = new IBTK::HierarchyAveragedDataManager(d_object_name + "::VelocityAveraging",
                                                                            d_velocity_cell_var,
                                                                            velocity_averaging_db,
                                                                            grid_geom,
                                                                            register_for_restart);
        d_velocity_product_average_manager =
            new IBTK::HierarchyAveragedDataManager(d_object_name + "::VelocityProductAveraging",
                                                   d_velocity_product_cell_var,
                                                   velocity_product_averaging_db,
                                                   grid_geom,
                                                   register_for_restart);
        break;
    case DataCentering::NODE:
        d_velocity_node_var = new NodeVariable<NDIM, double>(d_object_name + "::U_nc", NDIM, false);
        d_velocity_product_node_var =
            new NodeVariable<NDIM, double>(d_object_name + "::UU_nc", SYM_TENSOR_DEPTH, false);
        d_velocity_mean_node_var = new NodeVariable<NDIM, double>(d_object_name + "::U_mean_nc", NDIM, false);
        d_velocity_product_mean_node_var =
            new NodeVariable<NDIM, double>(d_object_name + "::UU_mean_nc", SYM_TENSOR_DEPTH, false);
        d_velocity_node_idx = var_db->registerVariableAndContext(d_velocity_node_var, ctx, IntVector<NDIM>(0));
        d_velocity_product_node_idx =
            var_db->registerVariableAndContext(d_velocity_product_node_var, ctx, IntVector<NDIM>(0));
        d_velocity_mean_node_idx =
            var_db->registerVariableAndContext(d_velocity_mean_node_var, ctx, IntVector<NDIM>(0));
        d_velocity_product_mean_node_idx =
            var_db->registerVariableAndContext(d_velocity_product_mean_node_var, ctx, IntVector<NDIM>(0));
        if (d_output_reynolds_stress)
        {
            d_reynolds_stress_node_var =
                new NodeVariable<NDIM, double>(d_object_name + "::ReynoldsStress_nc", NDIM * NDIM, false);
            d_reynolds_stress_plot_idx =
                var_db->registerVariableAndContext(d_reynolds_stress_node_var, plot_ctx, IntVector<NDIM>(0));
        }
        if (d_output_tke)
        {
            d_tke_node_var = new NodeVariable<NDIM, double>(d_object_name + "::TKE_nc", 1, false);
            d_tke_plot_idx = var_db->registerVariableAndContext(d_tke_node_var, plot_ctx, IntVector<NDIM>(0));
        }
        d_velocity_average_manager = new IBTK::HierarchyAveragedDataManager(d_object_name + "::VelocityAveraging",
                                                                            d_velocity_node_var,
                                                                            velocity_averaging_db,
                                                                            grid_geom,
                                                                            register_for_restart);
        d_velocity_product_average_manager =
            new IBTK::HierarchyAveragedDataManager(d_object_name + "::VelocityProductAveraging",
                                                   d_velocity_product_node_var,
                                                   velocity_product_averaging_db,
                                                   grid_geom,
                                                   register_for_restart);
        break;
    default:
        TBOX_ERROR(d_object_name << ": unsupported data centering enum value " << static_cast<int>(d_data_centering)
                                 << "\n");
    }

    if (d_velocity_average_manager->getSnapshotTimePoints() !=
        d_velocity_product_average_manager->getSnapshotTimePoints())
    {
        TBOX_ERROR(d_object_name
                   << ": velocity and velocity-product averaging configurations must use the same snapshot times.\n");
    }
}

bool
INSAveragingTurbulenceStatistics::updateStatistics(const int U_idx,
                                                   const Pointer<SideVariable<NDIM, double>> /*U_var*/,
                                                   const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                                                   const double data_time,
                                                   const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                   const Pointer<IBTK::HierarchyMathOps> hier_math_ops)
{
    if (!shouldUpdateStatistics(data_time)) return false;

    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
    int velocity_idx = IBTK::invalid_index;
    int velocity_product_idx = IBTK::invalid_index;
    int wgt_idx = IBTK::invalid_index;
    switch (d_data_centering)
    {
    case DataCentering::CELL:
        velocity_idx = d_velocity_cell_idx;
        velocity_product_idx = d_velocity_product_cell_idx;
        wgt_idx = wgt_cc_idx;
        break;
    case DataCentering::NODE:
        velocity_idx = d_velocity_node_idx;
        velocity_product_idx = d_velocity_product_node_idx;
        break;
    default:
        TBOX_ERROR(d_object_name << ": unsupported data centering enum value " << static_cast<int>(d_data_centering)
                                 << "\n");
    }

    IBTK::allocate_patch_data(
        { d_velocity_side_scratch_idx, velocity_idx, velocity_product_idx }, data_time, hierarchy);

    Pointer<PatchHierarchy<NDIM>> hierarchy_nc = hierarchy;
    auto* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double>> hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(d_velocity_side_scratch_var, hierarchy_nc, true);
    hier_sc_data_ops->copyData(d_velocity_side_scratch_idx, U_idx);

    using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    IBTK::HierarchyGhostCellInterpolation ghost_fill;
    ghost_fill.initializeOperatorState(ITC(d_velocity_side_scratch_idx,
                                           "CONSERVATIVE_LINEAR_REFINE",
                                           true,
                                           "CONSERVATIVE_COARSEN",
                                           "LINEAR",
                                           false,
                                           velocity_bc_coefs),
                                       hierarchy);
    ghost_fill.fillData(data_time);

    switch (d_data_centering)
    {
    case DataCentering::CELL:
        hier_math_ops->interp(d_velocity_cell_idx,
                              d_velocity_cell_var,
                              d_velocity_side_scratch_idx,
                              d_velocity_side_scratch_var,
                              nullptr,
                              data_time,
                              true);
        break;
    case DataCentering::NODE:
        hier_math_ops->interp(d_velocity_node_idx,
                              d_velocity_node_var,
                              true,
                              d_velocity_side_scratch_idx,
                              d_velocity_side_scratch_var,
                              nullptr,
                              data_time,
                              true);
        break;
    default:
        TBOX_ERROR(d_object_name << ": unsupported data centering enum value " << static_cast<int>(d_data_centering)
                                 << "\n");
    }

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            switch (d_data_centering)
            {
            case DataCentering::CELL:
            {
                Pointer<CellData<NDIM, double>> velocity_data = patch->getPatchData(d_velocity_cell_idx);
                Pointer<CellData<NDIM, double>> velocity_product_data =
                    patch->getPatchData(d_velocity_product_cell_idx);
                for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = *ci;
                    fill_symmetric_second_moment(*velocity_product_data, *velocity_data, idx);
                }
                break;
            }
            case DataCentering::NODE:
            {
                Pointer<NodeData<NDIM, double>> velocity_data = patch->getPatchData(d_velocity_node_idx);
                Pointer<NodeData<NDIM, double>> velocity_product_data =
                    patch->getPatchData(d_velocity_product_node_idx);
                for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
                {
                    const NodeIndex<NDIM>& idx = *ni;
                    fill_symmetric_second_moment(*velocity_product_data, *velocity_data, idx);
                }
                break;
            }
            default:
                TBOX_ERROR(d_object_name << ": unsupported data centering enum value "
                                         << static_cast<int>(d_data_centering) << "\n");
            }
        }
    }

    const bool U_steady =
        d_velocity_average_manager->updateTimeAveragedSnapshot(velocity_idx, data_time, hierarchy, wgt_idx);
    const bool UU_steady = d_velocity_product_average_manager->updateTimeAveragedSnapshot(
        velocity_product_idx, data_time, hierarchy, wgt_idx);
    d_has_statistics_samples = true;

    IBTK::deallocate_patch_data({ d_velocity_side_scratch_idx, velocity_idx, velocity_product_idx }, hierarchy);
    return U_steady && UU_steady;
}

bool
INSAveragingTurbulenceStatistics::isAtSteadyState() const
{
    return d_has_statistics_samples && d_velocity_average_manager->isAtPeriodicSteadyState() &&
           d_velocity_product_average_manager->isAtPeriodicSteadyState();
}

void
INSAveragingTurbulenceStatistics::registerVisItDataWriter(Pointer<VisItDataWriter<NDIM>> visit_writer)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(visit_writer);
#endif
    d_visit_writer = visit_writer;
    if (d_plot_quantities_registered) return;

    if (d_output_reynolds_stress)
    {
        d_visit_writer->registerPlotQuantity(d_reynolds_stress_output_name, "TENSOR", d_reynolds_stress_plot_idx);
    }
    if (d_output_tke)
    {
        d_visit_writer->registerPlotQuantity(d_tke_output_name, "SCALAR", d_tke_plot_idx);
    }
    d_plot_quantities_registered = true;
}

void
INSAveragingTurbulenceStatistics::setupPlotData(const double data_time,
                                                const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                const Pointer<IBTK::HierarchyMathOps> hier_math_ops)
{
    if (!d_output_reynolds_stress && !d_output_tke) return;

    if (d_output_reynolds_stress) IBTK::allocate_patch_data(d_reynolds_stress_plot_idx, data_time, hierarchy);
    if (d_output_tke) IBTK::allocate_patch_data(d_tke_plot_idx, data_time, hierarchy);

    if (!haveStoredSnapshots())
    {
        zeroPlotData(hierarchy);
        return;
    }

    switch (d_data_centering)
    {
    case DataCentering::CELL:
        if (d_output_reynolds_stress)
        {
            fillReynoldsStressSnapshot(
                d_reynolds_stress_plot_idx, d_reynolds_stress_cell_var, data_time, hierarchy, hier_math_ops);
        }
        if (d_output_tke)
        {
            fillTurbulentKineticEnergySnapshot(d_tke_plot_idx, d_tke_cell_var, data_time, hierarchy, hier_math_ops);
        }
        break;
    case DataCentering::NODE:
        if (d_output_reynolds_stress)
        {
            fillReynoldsStressSnapshot(
                d_reynolds_stress_plot_idx, d_reynolds_stress_node_var, data_time, hierarchy, hier_math_ops);
        }
        if (d_output_tke)
        {
            fillTurbulentKineticEnergySnapshot(d_tke_plot_idx, d_tke_node_var, data_time, hierarchy, hier_math_ops);
        }
        break;
    default:
        TBOX_ERROR(d_object_name << ": unsupported data centering enum value " << static_cast<int>(d_data_centering)
                                 << "\n");
    }
}

void
INSAveragingTurbulenceStatistics::deallocatePlotData(const Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    if (d_output_reynolds_stress) IBTK::deallocate_patch_data(d_reynolds_stress_plot_idx, hierarchy);
    if (d_output_tke) IBTK::deallocate_patch_data(d_tke_plot_idx, hierarchy);
}

DataCentering
INSAveragingTurbulenceStatistics::getDataCentering() const
{
    return d_data_centering;
}

IBTK::HierarchyAveragedDataManager&
INSAveragingTurbulenceStatistics::getAveragedVelocityManager()
{
    return *d_velocity_average_manager;
}

const IBTK::HierarchyAveragedDataManager&
INSAveragingTurbulenceStatistics::getAveragedVelocityManager() const
{
    return *d_velocity_average_manager;
}

IBTK::HierarchyAveragedDataManager&
INSAveragingTurbulenceStatistics::getAveragedVelocityProductManager()
{
    return *d_velocity_product_average_manager;
}

const IBTK::HierarchyAveragedDataManager&
INSAveragingTurbulenceStatistics::getAveragedVelocityProductManager() const
{
    return *d_velocity_product_average_manager;
}

void
INSAveragingTurbulenceStatistics::fillReynoldsStressSnapshot(const int R_idx,
                                                             const Pointer<CellVariable<NDIM, double>> /*R_var*/,
                                                             const double time,
                                                             const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                             const Pointer<IBTK::HierarchyMathOps> /*hier_math_ops*/,
                                                             const double tol) const
{
    if (d_data_centering != DataCentering::CELL)
    {
        TBOX_ERROR(d_object_name << "::fillReynoldsStressSnapshot(): cell-centered Reynolds stresses requested, but "
                                 << "analysis_centering = "
                                 << IBAMR::enum_to_string<IBAMR::DataCentering>(d_data_centering) << "\n");
    }

    fillAveragedSnapshots(time, hierarchy, tol);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<CellData<NDIM, double>> U_mean_data = patch->getPatchData(d_velocity_mean_cell_idx);
            Pointer<CellData<NDIM, double>> UU_mean_data = patch->getPatchData(d_velocity_product_mean_cell_idx);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = *ci;
                fill_reynolds_tensor(*R_data, *U_mean_data, *UU_mean_data, idx, d_object_name);
            }
        }
    }

    IBTK::deallocate_patch_data({ d_velocity_mean_cell_idx,
                                  d_velocity_cell_idx,
                                  d_velocity_product_mean_cell_idx,
                                  d_velocity_product_cell_idx },
                                hierarchy);
}

void
INSAveragingTurbulenceStatistics::fillReynoldsStressSnapshot(const int R_idx,
                                                             const Pointer<NodeVariable<NDIM, double>> /*R_var*/,
                                                             const double time,
                                                             const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                             const Pointer<IBTK::HierarchyMathOps> /*hier_math_ops*/,
                                                             const double tol) const
{
    if (d_data_centering != DataCentering::NODE)
    {
        TBOX_ERROR(d_object_name << "::fillReynoldsStressSnapshot(): node-centered Reynolds stresses requested, but "
                                 << "analysis_centering = "
                                 << IBAMR::enum_to_string<IBAMR::DataCentering>(d_data_centering) << "\n");
    }

    fillAveragedSnapshots(time, hierarchy, tol);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<NodeData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<NodeData<NDIM, double>> U_mean_data = patch->getPatchData(d_velocity_mean_node_idx);
            Pointer<NodeData<NDIM, double>> UU_mean_data = patch->getPatchData(d_velocity_product_mean_node_idx);

            for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
            {
                const NodeIndex<NDIM>& idx = *ni;
                fill_reynolds_tensor(*R_data, *U_mean_data, *UU_mean_data, idx, d_object_name);
            }
        }
    }

    IBTK::deallocate_patch_data({ d_velocity_mean_node_idx,
                                  d_velocity_node_idx,
                                  d_velocity_product_mean_node_idx,
                                  d_velocity_product_node_idx },
                                hierarchy);
}

void
INSAveragingTurbulenceStatistics::fillReynoldsStressSnapshot(const int R_idx,
                                                             const Pointer<Variable<NDIM>> R_var,
                                                             const double time,
                                                             const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                             const Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                             const double tol) const
{
    Pointer<CellVariable<NDIM, double>> R_cc_var = R_var;
    if (R_cc_var)
    {
        fillReynoldsStressSnapshot(R_idx, R_cc_var, time, hierarchy, hier_math_ops, tol);
        return;
    }
    Pointer<NodeVariable<NDIM, double>> R_nc_var = R_var;
    if (R_nc_var)
    {
        fillReynoldsStressSnapshot(R_idx, R_nc_var, time, hierarchy, hier_math_ops, tol);
        return;
    }
    TBOX_ERROR(d_object_name << "::fillReynoldsStressSnapshot(): unsupported destination variable type\n");
}

void
INSAveragingTurbulenceStatistics::fillTurbulentKineticEnergySnapshot(
    const int k_idx,
    const Pointer<CellVariable<NDIM, double>> /*k_var*/,
    const double time,
    const Pointer<PatchHierarchy<NDIM>> hierarchy,
    const Pointer<IBTK::HierarchyMathOps> /*hier_math_ops*/,
    const double tol) const
{
    if (d_data_centering != DataCentering::CELL)
    {
        TBOX_ERROR(d_object_name << "::fillTurbulentKineticEnergySnapshot(): cell-centered TKE requested, but "
                                 << "analysis_centering = "
                                 << IBAMR::enum_to_string<IBAMR::DataCentering>(d_data_centering) << "\n");
    }

    fillAveragedSnapshots(time, hierarchy, tol);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> k_data = patch->getPatchData(k_idx);
            Pointer<CellData<NDIM, double>> U_data = patch->getPatchData(d_velocity_mean_cell_idx);
            Pointer<CellData<NDIM, double>> UU_data = patch->getPatchData(d_velocity_product_mean_cell_idx);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = *ci;
                double trace = 0.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    trace += (*UU_data)(idx, d) - (*U_data)(idx, d) * (*U_data)(idx, d);
                }
                (*k_data)(idx) = 0.5 * trace;
            }
        }
    }
    IBTK::deallocate_patch_data({ d_velocity_mean_cell_idx,
                                  d_velocity_cell_idx,
                                  d_velocity_product_mean_cell_idx,
                                  d_velocity_product_cell_idx },
                                hierarchy);
}

void
INSAveragingTurbulenceStatistics::fillTurbulentKineticEnergySnapshot(
    const int k_idx,
    const Pointer<NodeVariable<NDIM, double>> /*k_var*/,
    const double time,
    const Pointer<PatchHierarchy<NDIM>> hierarchy,
    const Pointer<IBTK::HierarchyMathOps> /*hier_math_ops*/,
    const double tol) const
{
    if (d_data_centering != DataCentering::NODE)
    {
        TBOX_ERROR(d_object_name << "::fillTurbulentKineticEnergySnapshot(): node-centered TKE requested, but "
                                 << "analysis_centering = "
                                 << IBAMR::enum_to_string<IBAMR::DataCentering>(d_data_centering) << "\n");
    }

    fillAveragedSnapshots(time, hierarchy, tol);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<NodeData<NDIM, double>> k_data = patch->getPatchData(k_idx);
            Pointer<NodeData<NDIM, double>> U_data = patch->getPatchData(d_velocity_mean_node_idx);
            Pointer<NodeData<NDIM, double>> UU_data = patch->getPatchData(d_velocity_product_mean_node_idx);
            for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
            {
                const NodeIndex<NDIM>& idx = *ni;
                double trace = 0.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    trace += (*UU_data)(idx, d) - (*U_data)(idx, d) * (*U_data)(idx, d);
                }
                (*k_data)(idx) = 0.5 * trace;
            }
        }
    }
    IBTK::deallocate_patch_data({ d_velocity_mean_node_idx,
                                  d_velocity_node_idx,
                                  d_velocity_product_mean_node_idx,
                                  d_velocity_product_node_idx },
                                hierarchy);
}

bool
INSAveragingTurbulenceStatistics::haveStoredSnapshots() const
{
    return d_velocity_average_manager->getSnapshotCache().getNumSnapshots() > 0 &&
           d_velocity_product_average_manager->getSnapshotCache().getNumSnapshots() > 0;
}

double
INSAveragingTurbulenceStatistics::mapToStoredTime(double time) const
{
    return map_to_period(d_period_start, d_period_end, time);
}

void
INSAveragingTurbulenceStatistics::fillAveragedSnapshots(const double time,
                                                        const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                        const double tol) const
{
    if (!haveStoredSnapshots())
    {
        TBOX_ERROR(d_object_name << "::fillAveragedSnapshots(): no statistics snapshots have been accumulated yet.\n");
    }

    const double snapshot_time = mapToStoredTime(time);
    switch (d_data_centering)
    {
    case DataCentering::CELL:
        IBTK::allocate_patch_data({ d_velocity_mean_cell_idx,
                                    d_velocity_cell_idx,
                                    d_velocity_product_mean_cell_idx,
                                    d_velocity_product_cell_idx },
                                  snapshot_time,
                                  hierarchy);
        fill_snapshot_at_requested_time(d_velocity_average_manager->getSnapshotCache(),
                                        d_velocity_mean_cell_idx,
                                        d_velocity_cell_idx,
                                        d_velocity_cell_var,
                                        snapshot_time,
                                        hierarchy,
                                        d_refine_type,
                                        d_period_length,
                                        tol);
        fill_snapshot_at_requested_time(d_velocity_product_average_manager->getSnapshotCache(),
                                        d_velocity_product_mean_cell_idx,
                                        d_velocity_product_cell_idx,
                                        d_velocity_product_cell_var,
                                        snapshot_time,
                                        hierarchy,
                                        d_refine_type,
                                        d_period_length,
                                        tol);
        break;
    case DataCentering::NODE:
        IBTK::allocate_patch_data({ d_velocity_mean_node_idx,
                                    d_velocity_node_idx,
                                    d_velocity_product_mean_node_idx,
                                    d_velocity_product_node_idx },
                                  snapshot_time,
                                  hierarchy);
        fill_snapshot_at_requested_time(d_velocity_average_manager->getSnapshotCache(),
                                        d_velocity_mean_node_idx,
                                        d_velocity_node_idx,
                                        d_velocity_node_var,
                                        snapshot_time,
                                        hierarchy,
                                        d_refine_type,
                                        d_period_length,
                                        tol);
        fill_snapshot_at_requested_time(d_velocity_product_average_manager->getSnapshotCache(),
                                        d_velocity_product_mean_node_idx,
                                        d_velocity_product_node_idx,
                                        d_velocity_product_node_var,
                                        snapshot_time,
                                        hierarchy,
                                        d_refine_type,
                                        d_period_length,
                                        tol);
        break;
    default:
        TBOX_ERROR(d_object_name << ": unsupported data centering enum value " << static_cast<int>(d_data_centering)
                                 << "\n");
    }
}

void
INSAveragingTurbulenceStatistics::zeroPlotData(const Pointer<PatchHierarchy<NDIM>> hierarchy) const
{
    Pointer<PatchHierarchy<NDIM>> hierarchy_nc = hierarchy;
    switch (d_data_centering)
    {
    case DataCentering::CELL:
        if (d_output_reynolds_stress)
        {
            Pointer<HierarchyDataOpsReal<NDIM, double>> hier_data_ops =
                HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(
                    d_reynolds_stress_cell_var, hierarchy_nc, true);
            hier_data_ops->resetLevels(0, hierarchy_nc->getFinestLevelNumber());
            hier_data_ops->setToScalar(d_reynolds_stress_plot_idx, 0.0);
        }
        if (d_output_tke)
        {
            Pointer<HierarchyDataOpsReal<NDIM, double>> hier_data_ops =
                HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(d_tke_cell_var, hierarchy_nc, true);
            hier_data_ops->resetLevels(0, hierarchy_nc->getFinestLevelNumber());
            hier_data_ops->setToScalar(d_tke_plot_idx, 0.0);
        }
        break;
    case DataCentering::NODE:
        if (d_output_reynolds_stress)
        {
            Pointer<HierarchyDataOpsReal<NDIM, double>> hier_data_ops =
                HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(
                    d_reynolds_stress_node_var, hierarchy_nc, true);
            hier_data_ops->resetLevels(0, hierarchy_nc->getFinestLevelNumber());
            hier_data_ops->setToScalar(d_reynolds_stress_plot_idx, 0.0);
        }
        if (d_output_tke)
        {
            Pointer<HierarchyDataOpsReal<NDIM, double>> hier_data_ops =
                HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(d_tke_node_var, hierarchy_nc, true);
            hier_data_ops->resetLevels(0, hierarchy_nc->getFinestLevelNumber());
            hier_data_ops->setToScalar(d_tke_plot_idx, 0.0);
        }
        break;
    default:
        TBOX_ERROR(d_object_name << ": unsupported data centering enum value " << static_cast<int>(d_data_centering)
                                 << "\n");
    }
}
} // namespace IBAMR
