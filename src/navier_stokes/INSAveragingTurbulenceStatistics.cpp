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

#include <array>
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
    d_data_centering = string_to_enum<DataCentering>(
        input_db->getStringWithDefault("analysis_centering", enum_to_string<DataCentering>(d_data_centering)));

    d_velocity_side_scratch_var = new SideVariable<NDIM, double>(d_object_name + "::U_sc_scratch");
    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext(d_object_name + "::CONTEXT");
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
        d_velocity_average_manager = new IBTK::HierarchyAveragedDataManager(
            d_object_name + "::VelocityAveraging", d_velocity_cell_var, input_db, grid_geom, register_for_restart);
        d_velocity_product_average_manager =
            new IBTK::HierarchyAveragedDataManager(d_object_name + "::VelocityProductAveraging",
                                                   d_velocity_product_cell_var,
                                                   input_db,
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
        d_velocity_average_manager = new IBTK::HierarchyAveragedDataManager(
            d_object_name + "::VelocityAveraging", d_velocity_node_var, input_db, grid_geom, register_for_restart);
        d_velocity_product_average_manager =
            new IBTK::HierarchyAveragedDataManager(d_object_name + "::VelocityProductAveraging",
                                                   d_velocity_product_node_var,
                                                   input_db,
                                                   grid_geom,
                                                   register_for_restart);
        break;
    default:
        TBOX_ERROR(d_object_name << ": unsupported data centering enum value " << static_cast<int>(d_data_centering)
                                 << "\n");
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

    IBTK::deallocate_patch_data({ d_velocity_side_scratch_idx, velocity_idx, velocity_product_idx }, hierarchy);
    return U_steady && UU_steady;
}

bool
INSAveragingTurbulenceStatistics::isAtSteadyState() const
{
    return d_velocity_average_manager->isAtPeriodicSteadyState() &&
           d_velocity_product_average_manager->isAtPeriodicSteadyState();
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
                                 << "analysis_centering = " << enum_to_string<DataCentering>(d_data_centering) << "\n");
    }

    const double snapshot_time = d_velocity_average_manager->getTimePoint(time, tol);
    IBTK::allocate_patch_data({ d_velocity_mean_cell_idx, d_velocity_product_mean_cell_idx }, snapshot_time, hierarchy);

    IBTK::fill_snapshot_on_hierarchy(d_velocity_average_manager->getSnapshotCache(),
                                     d_velocity_mean_cell_idx,
                                     snapshot_time,
                                     hierarchy,
                                     d_refine_type,
                                     tol);
    IBTK::fill_snapshot_on_hierarchy(d_velocity_product_average_manager->getSnapshotCache(),
                                     d_velocity_product_mean_cell_idx,
                                     snapshot_time,
                                     hierarchy,
                                     d_refine_type,
                                     tol);

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

    IBTK::deallocate_patch_data({ d_velocity_mean_cell_idx, d_velocity_product_mean_cell_idx }, hierarchy);
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
                                 << "analysis_centering = " << enum_to_string<DataCentering>(d_data_centering) << "\n");
    }

    const double snapshot_time = d_velocity_average_manager->getTimePoint(time, tol);
    IBTK::allocate_patch_data({ d_velocity_mean_node_idx, d_velocity_product_mean_node_idx }, snapshot_time, hierarchy);

    IBTK::fill_snapshot_on_hierarchy(d_velocity_average_manager->getSnapshotCache(),
                                     d_velocity_mean_node_idx,
                                     snapshot_time,
                                     hierarchy,
                                     d_refine_type,
                                     tol);
    IBTK::fill_snapshot_on_hierarchy(d_velocity_product_average_manager->getSnapshotCache(),
                                     d_velocity_product_mean_node_idx,
                                     snapshot_time,
                                     hierarchy,
                                     d_refine_type,
                                     tol);

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

    IBTK::deallocate_patch_data({ d_velocity_mean_node_idx, d_velocity_product_mean_node_idx }, hierarchy);
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
                                 << "analysis_centering = " << enum_to_string<DataCentering>(d_data_centering) << "\n");
    }

    const double snapshot_time = d_velocity_average_manager->getTimePoint(time, tol);
    IBTK::allocate_patch_data({ d_velocity_mean_cell_idx, d_velocity_product_mean_cell_idx }, snapshot_time, hierarchy);
    IBTK::fill_snapshot_on_hierarchy(d_velocity_average_manager->getSnapshotCache(),
                                     d_velocity_mean_cell_idx,
                                     snapshot_time,
                                     hierarchy,
                                     d_refine_type,
                                     tol);
    IBTK::fill_snapshot_on_hierarchy(d_velocity_product_average_manager->getSnapshotCache(),
                                     d_velocity_product_mean_cell_idx,
                                     snapshot_time,
                                     hierarchy,
                                     d_refine_type,
                                     tol);

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
    IBTK::deallocate_patch_data({ d_velocity_mean_cell_idx, d_velocity_product_mean_cell_idx }, hierarchy);
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
                                 << "analysis_centering = " << enum_to_string<DataCentering>(d_data_centering) << "\n");
    }

    const double snapshot_time = d_velocity_average_manager->getTimePoint(time, tol);
    IBTK::allocate_patch_data({ d_velocity_mean_node_idx, d_velocity_product_mean_node_idx }, snapshot_time, hierarchy);
    IBTK::fill_snapshot_on_hierarchy(d_velocity_average_manager->getSnapshotCache(),
                                     d_velocity_mean_node_idx,
                                     snapshot_time,
                                     hierarchy,
                                     d_refine_type,
                                     tol);
    IBTK::fill_snapshot_on_hierarchy(d_velocity_product_average_manager->getSnapshotCache(),
                                     d_velocity_product_mean_node_idx,
                                     snapshot_time,
                                     hierarchy,
                                     d_refine_type,
                                     tol);

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
    IBTK::deallocate_patch_data({ d_velocity_mean_node_idx, d_velocity_product_mean_node_idx }, hierarchy);
}
} // namespace IBAMR
