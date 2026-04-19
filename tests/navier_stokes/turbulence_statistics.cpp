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

// Config files
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <HierarchyDataOpsManager.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/INSAveragingTurbulenceStatistics.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/snapshot_utilities.h>

// Set up application namespace declarations
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <set>
#include <string>

#include <ibamr/app_namespaces.h>

namespace
{
constexpr int TENSOR_DEPTH = NDIM * (NDIM + 1) / 2;
constexpr double PI = 3.141592653589793238462643383279502884;

using Vector = std::array<double, NDIM>;
using SymTensor = std::array<double, TENSOR_DEPTH>;
using LinearVectorCoefficients = std::array<std::array<double, NDIM + 1>, NDIM>;

struct ManufacturedVelocityCoefficients
{
    LinearVectorCoefficients mean;
    LinearVectorCoefficients cosine;
    LinearVectorCoefficients sine;
};

struct ManufacturedStatisticsSchedule
{
    double sample_start_time = 0.0;
    double sample_period = 0.0;
    int num_samples = 0;
    int num_periods = 1;
    bool periodic_mode = false;
    double statistics_start_time = 0.0;
    double period_start = 0.0;
    double period_end = 0.0;
};

ManufacturedVelocityCoefficients
get_manufactured_velocity_coefficients()
{
#if (NDIM == 2)
    return { { { { 0.45, 0.10, -0.05 }, { -0.20, 0.07, 0.04 } } },
             { { { 0.70, -0.08, 0.06 }, { -0.35, 0.05, -0.09 } } },
             { { { -0.25, 0.09, 0.03 }, { 0.55, -0.04, 0.08 } } } };
#endif
#if (NDIM == 3)
    return { { { { 0.45, 0.05, -0.03, 0.04 }, { -0.20, 0.02, 0.06, -0.05 }, { 0.15, -0.04, 0.03, 0.05 } } },
             { { { 0.70, -0.06, 0.04, 0.03 }, { -0.35, 0.05, -0.07, 0.02 }, { 0.28, 0.03, 0.05, -0.04 } } },
             { { { -0.25, 0.08, 0.02, -0.03 }, { 0.55, -0.04, 0.06, 0.05 }, { -0.18, 0.05, -0.02, 0.04 } } } };
#endif
}

Vector
evaluate_linear_vector(const LinearVectorCoefficients& coeffs, const Vector& X)
{
    Vector values = {};
    for (int comp = 0; comp < NDIM; ++comp)
    {
        values[comp] = coeffs[comp][0];
        for (int d = 0; d < NDIM; ++d)
        {
            values[comp] += coeffs[comp][d + 1] * X[d];
        }
    }
    return values;
}

Vector
evaluate_velocity(const ManufacturedVelocityCoefficients& coeffs,
                  const Vector& X,
                  const double time,
                  const double period)
{
    const double theta = 2.0 * PI * time / period;
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    const Vector mean = evaluate_linear_vector(coeffs.mean, X);
    const Vector cosine = evaluate_linear_vector(coeffs.cosine, X);
    const Vector sine = evaluate_linear_vector(coeffs.sine, X);

    Vector velocity = {};
    for (int comp = 0; comp < NDIM; ++comp)
    {
        velocity[comp] = mean[comp] + cosine[comp] * c + sine[comp] * s;
    }
    return velocity;
}

double
map_to_period_interval(const double period_start, const double period_end, double time)
{
    const double period_length = period_end - period_start;
    if (IBTK::abs_equal_eps(period_length, 0.0)) return period_start;

    time -= period_start;
    time = std::fmod(time, period_length);
    if (time < 0.0) time += period_length;
    return period_start + time;
}

template <class DataType, class IndexType>
void
fill_symmetric_second_moment(DataType& velocity_product_data, const DataType& velocity_data, const IndexType& idx)
{
#if (NDIM == 2)
    const double u0 = velocity_data(idx, 0);
    const double u1 = velocity_data(idx, 1);
    velocity_product_data(idx, 0) = u0 * u0;
    velocity_product_data(idx, 1) = u1 * u1;
    velocity_product_data(idx, 2) = u0 * u1;
#endif
#if (NDIM == 3)
    const double u0 = velocity_data(idx, 0);
    const double u1 = velocity_data(idx, 1);
    const double u2 = velocity_data(idx, 2);
    velocity_product_data(idx, 0) = u0 * u0;
    velocity_product_data(idx, 1) = u1 * u1;
    velocity_product_data(idx, 2) = u2 * u2;
    velocity_product_data(idx, 3) = u1 * u2;
    velocity_product_data(idx, 4) = u0 * u2;
    velocity_product_data(idx, 5) = u0 * u1;
#endif
}

void
fill_side_velocity(const int U_idx,
                   const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                   const ManufacturedVelocityCoefficients& coeffs,
                   const double data_time,
                   const double period)
{
    IBTK::allocate_patch_data(U_idx, data_time, patch_hierarchy);

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> U_data = patch->getPatchData(U_idx);
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const x_lower = patch_geom->getXLower();
            const double* const dx = patch_geom->getDx();
            const hier::Index<NDIM>& patch_lower = patch->getBox().lower();

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                {
                    const SideIndex<NDIM> side_idx = si();
                    Vector X = {};
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] +
                               dx[d] * (static_cast<double>(side_idx(d) - patch_lower(d)) + (axis == d ? 0.0 : 0.5));
                    }
                    const Vector velocity = evaluate_velocity(coeffs, X, data_time, period);
                    (*U_data)(side_idx) = velocity[axis];
                }
            }
        }
    }
}

double verify_cell_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                              const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                              const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                              const ManufacturedVelocityCoefficients& coeffs,
                              const ManufacturedStatisticsSchedule& schedule,
                              const int U_idx,
                              const Pointer<SideVariable<NDIM, double>>& U_var,
                              const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                              const double snapshot_time);

double verify_node_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                              const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                              const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                              const ManufacturedVelocityCoefficients& coeffs,
                              const ManufacturedStatisticsSchedule& schedule,
                              const int U_idx,
                              const Pointer<SideVariable<NDIM, double>>& U_var,
                              const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                              const double snapshot_time);

double verify_periodic_phase_cell_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                                             const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                                             const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                                             const ManufacturedVelocityCoefficients& coeffs,
                                             const ManufacturedStatisticsSchedule& schedule,
                                             const int U_idx,
                                             const Pointer<SideVariable<NDIM, double>>& U_var,
                                             const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs);

double verify_periodic_phase_node_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                                             const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                                             const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                                             const ManufacturedVelocityCoefficients& coeffs,
                                             const ManufacturedStatisticsSchedule& schedule,
                                             const int U_idx,
                                             const Pointer<SideVariable<NDIM, double>>& U_var,
                                             const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs);

template <class MeanDataType, class ProductDataType, class IndexType>
SymTensor
compute_reynolds_tensor(const MeanDataType& mean_velocity_data,
                        const ProductDataType& mean_velocity_product_data,
                        const IndexType& idx)
{
    SymTensor R = {};

#if (NDIM == 2)
    const double u0 = mean_velocity_data(idx, 0);
    const double u1 = mean_velocity_data(idx, 1);
    R[0] = mean_velocity_product_data(idx, 0) - u0 * u0;
    R[1] = mean_velocity_product_data(idx, 1) - u1 * u1;
    R[2] = mean_velocity_product_data(idx, 2) - u0 * u1;
#endif

#if (NDIM == 3)
    const double u0 = mean_velocity_data(idx, 0);
    const double u1 = mean_velocity_data(idx, 1);
    const double u2 = mean_velocity_data(idx, 2);
    R[0] = mean_velocity_product_data(idx, 0) - u0 * u0;
    R[1] = mean_velocity_product_data(idx, 1) - u1 * u1;
    R[2] = mean_velocity_product_data(idx, 2) - u2 * u2;
    R[3] = mean_velocity_product_data(idx, 3) - u1 * u2;
    R[4] = mean_velocity_product_data(idx, 4) - u0 * u2;
    R[5] = mean_velocity_product_data(idx, 5) - u0 * u1;
#endif

    return R;
}

double
compute_tke(const SymTensor& reynolds)
{
#if (NDIM == 2)
    return 0.5 * (reynolds[0] + reynolds[1]);
#endif
#if (NDIM == 3)
    return 0.5 * (reynolds[0] + reynolds[1] + reynolds[2]);
#endif
}

void
fill_cell_reference_sample(const int U_idx,
                           const Pointer<SideVariable<NDIM, double>>& U_var,
                           const int side_scratch_idx,
                           const int velocity_sample_idx,
                           const Pointer<CellVariable<NDIM, double>>& velocity_sample_var,
                           const int velocity_product_sample_idx,
                           const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                           const double data_time,
                           const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                           const Pointer<IBTK::HierarchyMathOps>& hier_math_ops)
{
    IBTK::allocate_patch_data(
        { side_scratch_idx, velocity_sample_idx, velocity_product_sample_idx }, data_time, patch_hierarchy);

    Pointer<PatchHierarchy<NDIM>> hierarchy_nc = patch_hierarchy;
    auto* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double>> hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(U_var, hierarchy_nc, true);
    hier_sc_data_ops->copyData(side_scratch_idx, U_idx);

    using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    IBTK::HierarchyGhostCellInterpolation ghost_fill;
    ghost_fill.initializeOperatorState(ITC(side_scratch_idx,
                                           "CONSERVATIVE_LINEAR_REFINE",
                                           true,
                                           "CONSERVATIVE_COARSEN",
                                           "LINEAR",
                                           false,
                                           velocity_bc_coefs),
                                       patch_hierarchy);
    ghost_fill.fillData(data_time);

    hier_math_ops->interp(velocity_sample_idx, velocity_sample_var, side_scratch_idx, U_var, nullptr, data_time, true);

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> velocity_data = patch->getPatchData(velocity_sample_idx);
            Pointer<CellData<NDIM, double>> velocity_product_data = patch->getPatchData(velocity_product_sample_idx);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = *ci;
                fill_symmetric_second_moment(*velocity_product_data, *velocity_data, idx);
            }
        }
    }
}

void
fill_node_reference_sample(const int U_idx,
                           const Pointer<SideVariable<NDIM, double>>& U_var,
                           const int side_scratch_idx,
                           const int velocity_sample_idx,
                           const Pointer<NodeVariable<NDIM, double>>& velocity_sample_var,
                           const int velocity_product_sample_idx,
                           const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                           const double data_time,
                           const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                           const Pointer<IBTK::HierarchyMathOps>& hier_math_ops)
{
    IBTK::allocate_patch_data(
        { side_scratch_idx, velocity_sample_idx, velocity_product_sample_idx }, data_time, patch_hierarchy);

    Pointer<PatchHierarchy<NDIM>> hierarchy_nc = patch_hierarchy;
    auto* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double>> hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(U_var, hierarchy_nc, true);
    hier_sc_data_ops->copyData(side_scratch_idx, U_idx);

    using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    IBTK::HierarchyGhostCellInterpolation ghost_fill;
    ghost_fill.initializeOperatorState(ITC(side_scratch_idx,
                                           "CONSERVATIVE_LINEAR_REFINE",
                                           true,
                                           "CONSERVATIVE_COARSEN",
                                           "LINEAR",
                                           false,
                                           velocity_bc_coefs),
                                       patch_hierarchy);
    ghost_fill.fillData(data_time);

    hier_math_ops->interp(
        velocity_sample_idx, velocity_sample_var, true, side_scratch_idx, U_var, nullptr, data_time, true);

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<NodeData<NDIM, double>> velocity_data = patch->getPatchData(velocity_sample_idx);
            Pointer<NodeData<NDIM, double>> velocity_product_data = patch->getPatchData(velocity_product_sample_idx);
            for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
            {
                const NodeIndex<NDIM>& idx = *ni;
                fill_symmetric_second_moment(*velocity_product_data, *velocity_data, idx);
            }
        }
    }
}

void
build_cell_reference_snapshot(const int U_idx,
                              const Pointer<SideVariable<NDIM, double>>& U_var,
                              const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                              const Pointer<CellVariable<NDIM, double>>& velocity_sample_var,
                              const Pointer<CellVariable<NDIM, double>>& velocity_ref_var,
                              const Pointer<CellVariable<NDIM, double>>& velocity_product_ref_var,
                              const int side_scratch_idx,
                              const int velocity_sample_idx,
                              const int velocity_ref_idx,
                              const int velocity_product_sample_idx,
                              const int velocity_product_ref_idx,
                              const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                              const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                              const ManufacturedVelocityCoefficients& coeffs,
                              const ManufacturedStatisticsSchedule& schedule,
                              const double snapshot_time)
{
    IBTK::allocate_patch_data({ velocity_ref_idx, velocity_product_ref_idx }, snapshot_time, patch_hierarchy);

    Pointer<PatchHierarchy<NDIM>> hierarchy_nc = patch_hierarchy;
    auto* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double>> velocity_data_ops =
        hier_ops_manager->getOperationsDouble(velocity_ref_var, hierarchy_nc, true);
    Pointer<HierarchyDataOpsReal<NDIM, double>> velocity_product_data_ops =
        hier_ops_manager->getOperationsDouble(velocity_product_ref_var, hierarchy_nc, true);

    int sample_count = 0;
    for (int period = 0; period < schedule.num_periods; ++period)
    {
        for (int sample = 0; sample < schedule.num_samples; ++sample)
        {
            const double phase_time = schedule.sample_start_time + schedule.sample_period *
                                                                       static_cast<double>(sample) /
                                                                       static_cast<double>(schedule.num_samples);
            const double data_time =
                schedule.periodic_mode ? phase_time + static_cast<double>(period) * schedule.sample_period : phase_time;
            if (data_time < schedule.statistics_start_time) continue;

            if (schedule.periodic_mode)
            {
                const double mapped_time =
                    map_to_period_interval(schedule.period_start, schedule.period_end, data_time);
                if (!IBTK::abs_equal_eps(mapped_time, snapshot_time, 1.0e-8)) continue;
            }

            fill_side_velocity(U_idx, patch_hierarchy, coeffs, data_time, schedule.sample_period);
            fill_cell_reference_sample(U_idx,
                                       U_var,
                                       side_scratch_idx,
                                       velocity_sample_idx,
                                       velocity_sample_var,
                                       velocity_product_sample_idx,
                                       velocity_bc_coefs,
                                       data_time,
                                       patch_hierarchy,
                                       hier_math_ops);

            if (sample_count == 0)
            {
                velocity_data_ops->copyData(velocity_ref_idx, velocity_sample_idx);
                velocity_product_data_ops->copyData(velocity_product_ref_idx, velocity_product_sample_idx);
            }
            else
            {
                const double N = static_cast<double>(sample_count);
                velocity_data_ops->linearSum(
                    velocity_ref_idx, 1.0 / (N + 1.0), velocity_sample_idx, N / (N + 1.0), velocity_ref_idx);
                velocity_product_data_ops->linearSum(velocity_product_ref_idx,
                                                     1.0 / (N + 1.0),
                                                     velocity_product_sample_idx,
                                                     N / (N + 1.0),
                                                     velocity_product_ref_idx);
            }
            ++sample_count;

            IBTK::deallocate_patch_data({ side_scratch_idx, velocity_sample_idx, velocity_product_sample_idx },
                                        patch_hierarchy);
        }
    }

    if (sample_count == 0)
    {
        TBOX_ERROR("No manufactured samples contributed to cell snapshot time " << snapshot_time << "\n");
    }
}

void
build_node_reference_snapshot(const int U_idx,
                              const Pointer<SideVariable<NDIM, double>>& U_var,
                              const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                              const Pointer<NodeVariable<NDIM, double>>& velocity_sample_var,
                              const Pointer<NodeVariable<NDIM, double>>& velocity_ref_var,
                              const Pointer<NodeVariable<NDIM, double>>& velocity_product_ref_var,
                              const int side_scratch_idx,
                              const int velocity_sample_idx,
                              const int velocity_ref_idx,
                              const int velocity_product_sample_idx,
                              const int velocity_product_ref_idx,
                              const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                              const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                              const ManufacturedVelocityCoefficients& coeffs,
                              const ManufacturedStatisticsSchedule& schedule,
                              const double snapshot_time)
{
    IBTK::allocate_patch_data({ velocity_ref_idx, velocity_product_ref_idx }, snapshot_time, patch_hierarchy);

    Pointer<PatchHierarchy<NDIM>> hierarchy_nc = patch_hierarchy;
    auto* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double>> velocity_data_ops =
        hier_ops_manager->getOperationsDouble(velocity_ref_var, hierarchy_nc, true);
    Pointer<HierarchyDataOpsReal<NDIM, double>> velocity_product_data_ops =
        hier_ops_manager->getOperationsDouble(velocity_product_ref_var, hierarchy_nc, true);

    int sample_count = 0;
    for (int period = 0; period < schedule.num_periods; ++period)
    {
        for (int sample = 0; sample < schedule.num_samples; ++sample)
        {
            const double phase_time = schedule.sample_start_time + schedule.sample_period *
                                                                       static_cast<double>(sample) /
                                                                       static_cast<double>(schedule.num_samples);
            const double data_time =
                schedule.periodic_mode ? phase_time + static_cast<double>(period) * schedule.sample_period : phase_time;
            if (data_time < schedule.statistics_start_time) continue;

            if (schedule.periodic_mode)
            {
                const double mapped_time =
                    map_to_period_interval(schedule.period_start, schedule.period_end, data_time);
                if (!IBTK::abs_equal_eps(mapped_time, snapshot_time, 1.0e-8)) continue;
            }

            fill_side_velocity(U_idx, patch_hierarchy, coeffs, data_time, schedule.sample_period);
            fill_node_reference_sample(U_idx,
                                       U_var,
                                       side_scratch_idx,
                                       velocity_sample_idx,
                                       velocity_sample_var,
                                       velocity_product_sample_idx,
                                       velocity_bc_coefs,
                                       data_time,
                                       patch_hierarchy,
                                       hier_math_ops);

            if (sample_count == 0)
            {
                velocity_data_ops->copyData(velocity_ref_idx, velocity_sample_idx);
                velocity_product_data_ops->copyData(velocity_product_ref_idx, velocity_product_sample_idx);
            }
            else
            {
                const double N = static_cast<double>(sample_count);
                velocity_data_ops->linearSum(
                    velocity_ref_idx, 1.0 / (N + 1.0), velocity_sample_idx, N / (N + 1.0), velocity_ref_idx);
                velocity_product_data_ops->linearSum(velocity_product_ref_idx,
                                                     1.0 / (N + 1.0),
                                                     velocity_product_sample_idx,
                                                     N / (N + 1.0),
                                                     velocity_product_ref_idx);
            }
            ++sample_count;

            IBTK::deallocate_patch_data({ side_scratch_idx, velocity_sample_idx, velocity_product_sample_idx },
                                        patch_hierarchy);
        }
    }

    if (sample_count == 0)
    {
        TBOX_ERROR("No manufactured samples contributed to node snapshot time " << snapshot_time << "\n");
    }
}

double
verify_cell_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                       const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                       const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                       const ManufacturedVelocityCoefficients& coeffs,
                       const ManufacturedStatisticsSchedule& schedule,
                       const int U_idx,
                       const Pointer<SideVariable<NDIM, double>>& U_var,
                       const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                       const double snapshot_time)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("EX7::VerificationContext");
    Pointer<CellVariable<NDIM, double>> U_mean_var = new CellVariable<NDIM, double>("U_mean_cc", NDIM);
    Pointer<CellVariable<NDIM, double>> R_var = new CellVariable<NDIM, double>("R_cc", TENSOR_DEPTH);
    Pointer<CellVariable<NDIM, double>> U_ref_var = new CellVariable<NDIM, double>("U_ref_cc", NDIM);
    Pointer<CellVariable<NDIM, double>> U_sample_var = new CellVariable<NDIM, double>("U_sample_cc", NDIM);
    Pointer<CellVariable<NDIM, double>> UU_ref_var = new CellVariable<NDIM, double>("UU_ref_cc", TENSOR_DEPTH);
    Pointer<CellVariable<NDIM, double>> UU_sample_var = new CellVariable<NDIM, double>("UU_sample_cc", TENSOR_DEPTH);
    const int U_mean_idx = var_db->registerVariableAndContext(U_mean_var, ctx, IntVector<NDIM>(1));
    const int R_idx = var_db->registerVariableAndContext(R_var, ctx, IntVector<NDIM>(0));
    const int U_ref_idx = var_db->registerVariableAndContext(U_ref_var, ctx, IntVector<NDIM>(1));
    const int U_sample_idx = var_db->registerVariableAndContext(U_sample_var, ctx, IntVector<NDIM>(1));
    const int UU_ref_idx = var_db->registerVariableAndContext(UU_ref_var, ctx, IntVector<NDIM>(0));
    const int UU_sample_idx = var_db->registerVariableAndContext(UU_sample_var, ctx, IntVector<NDIM>(0));
    const int U_side_scratch_idx = var_db->registerVariableAndContext(U_var, ctx, IntVector<NDIM>(1));

    IBTK::allocate_patch_data({ U_mean_idx, R_idx }, snapshot_time, patch_hierarchy);

    IBTK::fill_snapshot_on_hierarchy(statistics->getAveragedVelocityManager().getSnapshotCache(),
                                     U_mean_idx,
                                     snapshot_time,
                                     patch_hierarchy,
                                     "CONSERVATIVE_LINEAR_REFINE");
    statistics->fillReynoldsStressSnapshot(R_idx, R_var, snapshot_time, patch_hierarchy, hier_math_ops);
    build_cell_reference_snapshot(U_idx,
                                  U_var,
                                  velocity_bc_coefs,
                                  U_sample_var,
                                  U_ref_var,
                                  UU_ref_var,
                                  U_side_scratch_idx,
                                  U_sample_idx,
                                  U_ref_idx,
                                  UU_sample_idx,
                                  UU_ref_idx,
                                  patch_hierarchy,
                                  hier_math_ops,
                                  coeffs,
                                  schedule,
                                  snapshot_time);

    double max_mean_error = 0.0;
    double max_reynolds_error = 0.0;
    double max_tke_error = 0.0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> U_mean_data = patch->getPatchData(U_mean_idx);
            Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<CellData<NDIM, double>> U_ref_data = patch->getPatchData(U_ref_idx);
            Pointer<CellData<NDIM, double>> UU_ref_data = patch->getPatchData(UU_ref_idx);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                const SymTensor R_ref = compute_reynolds_tensor(*U_ref_data, *UU_ref_data, idx);
                const double k_ref = compute_tke(R_ref);

                for (int comp = 0; comp < NDIM; ++comp)
                {
                    max_mean_error =
                        std::max(max_mean_error, std::abs((*U_mean_data)(idx, comp) - (*U_ref_data)(idx, comp)));
                }
                for (int comp = 0; comp < TENSOR_DEPTH; ++comp)
                {
                    max_reynolds_error = std::max(max_reynolds_error, std::abs((*R_data)(idx, comp) - R_ref[comp]));
                }

#if (NDIM == 2)
                const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1));
#endif
#if (NDIM == 3)
                const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1) + (*R_data)(idx, 2));
#endif
                max_tke_error = std::max(max_tke_error, std::abs(k_num - k_ref));
            }
        }
    }

    IBTK::deallocate_patch_data({ U_mean_idx, R_idx, U_ref_idx, UU_ref_idx }, patch_hierarchy);

    pout << "CELL statistics errors: |<U>-exact|_max = " << max_mean_error << ", |R-exact|_max = " << max_reynolds_error
         << ", |k-exact|_max = " << max_tke_error << "\n";
    return std::max(max_mean_error, std::max(max_reynolds_error, max_tke_error));
}

double
verify_node_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                       const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                       const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                       const ManufacturedVelocityCoefficients& coeffs,
                       const ManufacturedStatisticsSchedule& schedule,
                       const int U_idx,
                       const Pointer<SideVariable<NDIM, double>>& U_var,
                       const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                       const double snapshot_time)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("EX7::VerificationContext");
    Pointer<NodeVariable<NDIM, double>> U_mean_var = new NodeVariable<NDIM, double>("U_mean_nc", NDIM, false);
    Pointer<NodeVariable<NDIM, double>> R_var = new NodeVariable<NDIM, double>("R_nc", TENSOR_DEPTH, false);
    Pointer<NodeVariable<NDIM, double>> U_ref_var = new NodeVariable<NDIM, double>("U_ref_nc", NDIM, false);
    Pointer<NodeVariable<NDIM, double>> U_sample_var = new NodeVariable<NDIM, double>("U_sample_nc", NDIM, false);
    Pointer<NodeVariable<NDIM, double>> UU_ref_var = new NodeVariable<NDIM, double>("UU_ref_nc", TENSOR_DEPTH, false);
    Pointer<NodeVariable<NDIM, double>> UU_sample_var =
        new NodeVariable<NDIM, double>("UU_sample_nc", TENSOR_DEPTH, false);
    const int U_mean_idx = var_db->registerVariableAndContext(U_mean_var, ctx, IntVector<NDIM>(1));
    const int R_idx = var_db->registerVariableAndContext(R_var, ctx, IntVector<NDIM>(0));
    const int U_ref_idx = var_db->registerVariableAndContext(U_ref_var, ctx, IntVector<NDIM>(1));
    const int U_sample_idx = var_db->registerVariableAndContext(U_sample_var, ctx, IntVector<NDIM>(1));
    const int UU_ref_idx = var_db->registerVariableAndContext(UU_ref_var, ctx, IntVector<NDIM>(0));
    const int UU_sample_idx = var_db->registerVariableAndContext(UU_sample_var, ctx, IntVector<NDIM>(0));
    const int U_side_scratch_idx = var_db->registerVariableAndContext(U_var, ctx, IntVector<NDIM>(1));

    IBTK::allocate_patch_data({ U_mean_idx, R_idx }, snapshot_time, patch_hierarchy);

    IBTK::fill_snapshot_on_hierarchy(statistics->getAveragedVelocityManager().getSnapshotCache(),
                                     U_mean_idx,
                                     snapshot_time,
                                     patch_hierarchy,
                                     "LINEAR_REFINE");
    statistics->fillReynoldsStressSnapshot(R_idx, R_var, snapshot_time, patch_hierarchy, hier_math_ops);
    build_node_reference_snapshot(U_idx,
                                  U_var,
                                  velocity_bc_coefs,
                                  U_sample_var,
                                  U_ref_var,
                                  UU_ref_var,
                                  U_side_scratch_idx,
                                  U_sample_idx,
                                  U_ref_idx,
                                  UU_sample_idx,
                                  UU_ref_idx,
                                  patch_hierarchy,
                                  hier_math_ops,
                                  coeffs,
                                  schedule,
                                  snapshot_time);

    double max_mean_error = 0.0;
    double max_reynolds_error = 0.0;
    double max_tke_error = 0.0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<NodeData<NDIM, double>> U_mean_data = patch->getPatchData(U_mean_idx);
            Pointer<NodeData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<NodeData<NDIM, double>> U_ref_data = patch->getPatchData(U_ref_idx);
            Pointer<NodeData<NDIM, double>> UU_ref_data = patch->getPatchData(UU_ref_idx);

            for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
            {
                const NodeIndex<NDIM>& idx = ni();
                const SymTensor R_ref = compute_reynolds_tensor(*U_ref_data, *UU_ref_data, idx);
                const double k_ref = compute_tke(R_ref);

                for (int comp = 0; comp < NDIM; ++comp)
                {
                    max_mean_error =
                        std::max(max_mean_error, std::abs((*U_mean_data)(idx, comp) - (*U_ref_data)(idx, comp)));
                }
                for (int comp = 0; comp < TENSOR_DEPTH; ++comp)
                {
                    max_reynolds_error = std::max(max_reynolds_error, std::abs((*R_data)(idx, comp) - R_ref[comp]));
                }

#if (NDIM == 2)
                const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1));
#endif
#if (NDIM == 3)
                const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1) + (*R_data)(idx, 2));
#endif
                max_tke_error = std::max(max_tke_error, std::abs(k_num - k_ref));
            }
        }
    }

    IBTK::deallocate_patch_data({ U_mean_idx, R_idx, U_ref_idx, UU_ref_idx }, patch_hierarchy);

    pout << "NODE statistics errors: |<U>-exact|_max = " << max_mean_error << ", |R-exact|_max = " << max_reynolds_error
         << ", |k-exact|_max = " << max_tke_error << "\n";
    return std::max(max_mean_error, std::max(max_reynolds_error, max_tke_error));
}

double
verify_periodic_phase_cell_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                                      const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                                      const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                                      const ManufacturedVelocityCoefficients& coeffs,
                                      const ManufacturedStatisticsSchedule& schedule,
                                      const int U_idx,
                                      const Pointer<SideVariable<NDIM, double>>& U_var,
                                      const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("EX7::PeriodicVerificationContext");
    Pointer<CellVariable<NDIM, double>> U_phase_var = new CellVariable<NDIM, double>("U_phase_cc", NDIM);
    Pointer<CellVariable<NDIM, double>> R_var = new CellVariable<NDIM, double>("R_phase_cc", TENSOR_DEPTH);
    Pointer<CellVariable<NDIM, double>> U_ref_var = new CellVariable<NDIM, double>("U_ref_phase_cc", NDIM);
    Pointer<CellVariable<NDIM, double>> U_sample_var = new CellVariable<NDIM, double>("U_sample_phase_cc", NDIM);
    Pointer<CellVariable<NDIM, double>> UU_ref_var = new CellVariable<NDIM, double>("UU_ref_phase_cc", TENSOR_DEPTH);
    Pointer<CellVariable<NDIM, double>> UU_sample_var =
        new CellVariable<NDIM, double>("UU_sample_phase_cc", TENSOR_DEPTH);
    const int U_phase_idx = var_db->registerVariableAndContext(U_phase_var, ctx, IntVector<NDIM>(1));
    const int R_idx = var_db->registerVariableAndContext(R_var, ctx, IntVector<NDIM>(0));
    const int U_ref_idx = var_db->registerVariableAndContext(U_ref_var, ctx, IntVector<NDIM>(1));
    const int U_sample_idx = var_db->registerVariableAndContext(U_sample_var, ctx, IntVector<NDIM>(1));
    const int UU_ref_idx = var_db->registerVariableAndContext(UU_ref_var, ctx, IntVector<NDIM>(0));
    const int UU_sample_idx = var_db->registerVariableAndContext(UU_sample_var, ctx, IntVector<NDIM>(0));
    const int U_side_scratch_idx = var_db->registerVariableAndContext(U_var, ctx, IntVector<NDIM>(1));

    double max_phase_error = 0.0;
    double max_reynolds_error = 0.0;
    double max_tke_error = 0.0;
    for (const double snapshot_time : statistics->getAveragedVelocityManager().getSnapshotTimePoints())
    {
        IBTK::allocate_patch_data({ U_phase_idx, R_idx }, snapshot_time, patch_hierarchy);

        IBTK::fill_snapshot_on_hierarchy(statistics->getAveragedVelocityManager().getSnapshotCache(),
                                         U_phase_idx,
                                         snapshot_time,
                                         patch_hierarchy,
                                         "CONSERVATIVE_LINEAR_REFINE");
        statistics->fillReynoldsStressSnapshot(R_idx, R_var, snapshot_time, patch_hierarchy, hier_math_ops);
        build_cell_reference_snapshot(U_idx,
                                      U_var,
                                      velocity_bc_coefs,
                                      U_sample_var,
                                      U_ref_var,
                                      UU_ref_var,
                                      U_side_scratch_idx,
                                      U_sample_idx,
                                      U_ref_idx,
                                      UU_sample_idx,
                                      UU_ref_idx,
                                      patch_hierarchy,
                                      hier_math_ops,
                                      coeffs,
                                      schedule,
                                      snapshot_time);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> U_phase_data = patch->getPatchData(U_phase_idx);
                Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
                Pointer<CellData<NDIM, double>> U_ref_data = patch->getPatchData(U_ref_idx);
                Pointer<CellData<NDIM, double>> UU_ref_data = patch->getPatchData(UU_ref_idx);

                for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    const SymTensor R_ref = compute_reynolds_tensor(*U_ref_data, *UU_ref_data, idx);
                    const double k_ref = compute_tke(R_ref);
                    for (int comp = 0; comp < NDIM; ++comp)
                    {
                        max_phase_error =
                            std::max(max_phase_error, std::abs((*U_phase_data)(idx, comp) - (*U_ref_data)(idx, comp)));
                    }
                    for (int comp = 0; comp < TENSOR_DEPTH; ++comp)
                    {
                        max_reynolds_error = std::max(max_reynolds_error, std::abs((*R_data)(idx, comp) - R_ref[comp]));
                    }
#if (NDIM == 2)
                    const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1));
#endif
#if (NDIM == 3)
                    const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1) + (*R_data)(idx, 2));
#endif
                    max_tke_error = std::max(max_tke_error, std::abs(k_num - k_ref));
                }
            }
        }

        IBTK::deallocate_patch_data({ U_phase_idx, R_idx, U_ref_idx, UU_ref_idx }, patch_hierarchy);
    }

    pout << "CELL periodic-phase errors: |U_phase-exact|_max = " << max_phase_error
         << ", |R|_max = " << max_reynolds_error << ", |k|_max = " << max_tke_error << "\n";
    return std::max(max_phase_error, std::max(max_reynolds_error, max_tke_error));
}

double
verify_periodic_phase_node_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                                      const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                                      const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                                      const ManufacturedVelocityCoefficients& coeffs,
                                      const ManufacturedStatisticsSchedule& schedule,
                                      const int U_idx,
                                      const Pointer<SideVariable<NDIM, double>>& U_var,
                                      const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("EX7::PeriodicVerificationContext");
    Pointer<NodeVariable<NDIM, double>> U_phase_var = new NodeVariable<NDIM, double>("U_phase_nc", NDIM, false);
    Pointer<NodeVariable<NDIM, double>> R_var = new NodeVariable<NDIM, double>("R_phase_nc", TENSOR_DEPTH, false);
    Pointer<NodeVariable<NDIM, double>> U_ref_var = new NodeVariable<NDIM, double>("U_ref_phase_nc", NDIM, false);
    Pointer<NodeVariable<NDIM, double>> U_sample_var = new NodeVariable<NDIM, double>("U_sample_phase_nc", NDIM, false);
    Pointer<NodeVariable<NDIM, double>> UU_ref_var =
        new NodeVariable<NDIM, double>("UU_ref_phase_nc", TENSOR_DEPTH, false);
    Pointer<NodeVariable<NDIM, double>> UU_sample_var =
        new NodeVariable<NDIM, double>("UU_sample_phase_nc", TENSOR_DEPTH, false);
    const int U_phase_idx = var_db->registerVariableAndContext(U_phase_var, ctx, IntVector<NDIM>(1));
    const int R_idx = var_db->registerVariableAndContext(R_var, ctx, IntVector<NDIM>(0));
    const int U_ref_idx = var_db->registerVariableAndContext(U_ref_var, ctx, IntVector<NDIM>(1));
    const int U_sample_idx = var_db->registerVariableAndContext(U_sample_var, ctx, IntVector<NDIM>(1));
    const int UU_ref_idx = var_db->registerVariableAndContext(UU_ref_var, ctx, IntVector<NDIM>(0));
    const int UU_sample_idx = var_db->registerVariableAndContext(UU_sample_var, ctx, IntVector<NDIM>(0));
    const int U_side_scratch_idx = var_db->registerVariableAndContext(U_var, ctx, IntVector<NDIM>(1));

    double max_phase_error = 0.0;
    double max_reynolds_error = 0.0;
    double max_tke_error = 0.0;
    for (const double snapshot_time : statistics->getAveragedVelocityManager().getSnapshotTimePoints())
    {
        IBTK::allocate_patch_data({ U_phase_idx, R_idx }, snapshot_time, patch_hierarchy);

        IBTK::fill_snapshot_on_hierarchy(statistics->getAveragedVelocityManager().getSnapshotCache(),
                                         U_phase_idx,
                                         snapshot_time,
                                         patch_hierarchy,
                                         "LINEAR_REFINE");
        statistics->fillReynoldsStressSnapshot(R_idx, R_var, snapshot_time, patch_hierarchy, hier_math_ops);
        build_node_reference_snapshot(U_idx,
                                      U_var,
                                      velocity_bc_coefs,
                                      U_sample_var,
                                      U_ref_var,
                                      UU_ref_var,
                                      U_side_scratch_idx,
                                      U_sample_idx,
                                      U_ref_idx,
                                      UU_sample_idx,
                                      UU_ref_idx,
                                      patch_hierarchy,
                                      hier_math_ops,
                                      coeffs,
                                      schedule,
                                      snapshot_time);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<NodeData<NDIM, double>> U_phase_data = patch->getPatchData(U_phase_idx);
                Pointer<NodeData<NDIM, double>> R_data = patch->getPatchData(R_idx);
                Pointer<NodeData<NDIM, double>> U_ref_data = patch->getPatchData(U_ref_idx);
                Pointer<NodeData<NDIM, double>> UU_ref_data = patch->getPatchData(UU_ref_idx);

                for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
                {
                    const NodeIndex<NDIM>& idx = ni();
                    const SymTensor R_ref = compute_reynolds_tensor(*U_ref_data, *UU_ref_data, idx);
                    const double k_ref = compute_tke(R_ref);
                    for (int comp = 0; comp < NDIM; ++comp)
                    {
                        max_phase_error =
                            std::max(max_phase_error, std::abs((*U_phase_data)(idx, comp) - (*U_ref_data)(idx, comp)));
                    }
                    for (int comp = 0; comp < TENSOR_DEPTH; ++comp)
                    {
                        max_reynolds_error = std::max(max_reynolds_error, std::abs((*R_data)(idx, comp) - R_ref[comp]));
                    }
#if (NDIM == 2)
                    const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1));
#endif
#if (NDIM == 3)
                    const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1) + (*R_data)(idx, 2));
#endif
                    max_tke_error = std::max(max_tke_error, std::abs(k_num - k_ref));
                }
            }
        }

        IBTK::deallocate_patch_data({ U_phase_idx, R_idx, U_ref_idx, UU_ref_idx }, patch_hierarchy);
    }

    pout << "NODE periodic-phase errors: |U_phase-exact|_max = " << max_phase_error
         << ", |R|_max = " << max_reynolds_error << ", |k|_max = " << max_tke_error << "\n";
    return std::max(max_phase_error, std::max(max_reynolds_error, max_tke_error));
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    Pointer<INSStaggeredHierarchyIntegrator> time_integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

    Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
        "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
    time_integrator->registerVelocityInitialConditions(u_init);
    Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
        "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
    time_integrator->registerPressureInitialConditions(p_init);

    time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

    Pointer<INSAveragingTurbulenceStatistics> statistics = time_integrator->getTurbulenceStatistics();
    if (statistics.isNull())
    {
        TBOX_ERROR(
            "Expected INSStaggeredHierarchyIntegrator to construct an AVERAGED_VELOCITY turbulence statistics "
            "object.\n");
    }

    Pointer<Database> verification_db = app_initializer->getComponentDatabase("ManufacturedStatistics");
    Pointer<Database> turbulence_statistics_db =
        app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator")->getDatabase("TurbulenceStatistics");
    const double sample_start_time = verification_db->getDoubleWithDefault("sample_start_time", 0.0);
    const double sample_period = verification_db->getDouble("sample_period");
    const int num_samples = verification_db->getInteger("num_samples");
    const int num_periods = verification_db->getIntegerWithDefault("num_periods", 1);
    const bool periodic_mode = verification_db->getBoolWithDefault("periodic_mode", false);
    const bool compact_output = verification_db->getBoolWithDefault("compact_output", false);
    const double verification_tol = verification_db->getDoubleWithDefault("verification_tol", 1.0e-10);
    const double statistics_start_time = turbulence_statistics_db->getDoubleWithDefault("statistics_start_time", 0.0);
    const double period_start = turbulence_statistics_db->getDouble("period_start");
    const double period_end = turbulence_statistics_db->getDouble("period_end");

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, double>> U_var = time_integrator->getVelocityVariable();
    const int U_idx = var_db->mapVariableAndContextToIndex(U_var, time_integrator->getCurrentContext());

    Pointer<IBTK::HierarchyMathOps> hier_math_ops = new IBTK::HierarchyMathOps("HierarchyMathOps", patch_hierarchy);
    hier_math_ops->resetLevels(0, patch_hierarchy->getFinestLevelNumber());

    const ManufacturedVelocityCoefficients coeffs = get_manufactured_velocity_coefficients();
    const ManufacturedStatisticsSchedule schedule = { sample_start_time, sample_period, num_samples,
                                                      num_periods,       periodic_mode, statistics_start_time,
                                                      period_start,      period_end };
    for (int period = 0; period < num_periods; ++period)
    {
        for (int sample = 0; sample < num_samples; ++sample)
        {
            const double phase_time =
                sample_start_time + sample_period * static_cast<double>(sample) / static_cast<double>(num_samples);
            const double data_time =
                periodic_mode ? phase_time + static_cast<double>(period) * sample_period : phase_time;
            fill_side_velocity(U_idx, patch_hierarchy, coeffs, data_time, sample_period);
            statistics->updateStatistics(U_idx,
                                         U_var,
                                         time_integrator->getVelocityBoundaryConditions(),
                                         data_time,
                                         patch_hierarchy,
                                         hier_math_ops);
        }
    }

    const DataCentering data_centering = statistics->getDataCentering();
    double max_error = std::numeric_limits<double>::quiet_NaN();
    if (periodic_mode)
    {
        switch (data_centering)
        {
        case DataCentering::CELL:
            max_error = verify_periodic_phase_cell_statistics(statistics,
                                                              patch_hierarchy,
                                                              hier_math_ops,
                                                              coeffs,
                                                              schedule,
                                                              U_idx,
                                                              U_var,
                                                              time_integrator->getVelocityBoundaryConditions());
            break;
        case DataCentering::NODE:
            max_error = verify_periodic_phase_node_statistics(statistics,
                                                              patch_hierarchy,
                                                              hier_math_ops,
                                                              coeffs,
                                                              schedule,
                                                              U_idx,
                                                              U_var,
                                                              time_integrator->getVelocityBoundaryConditions());
            break;
        default:
            TBOX_ERROR("Unsupported analysis centering " << IBAMR::enum_to_string<IBAMR::DataCentering>(data_centering)
                                                         << "\n");
        }
    }
    else
    {
        const double snapshot_time = *(statistics->getAveragedVelocityManager().getSnapshotTimePoints().begin());
        switch (data_centering)
        {
        case DataCentering::CELL:
            max_error = verify_cell_statistics(statistics,
                                               patch_hierarchy,
                                               hier_math_ops,
                                               coeffs,
                                               schedule,
                                               U_idx,
                                               U_var,
                                               time_integrator->getVelocityBoundaryConditions(),
                                               snapshot_time);
            break;
        case DataCentering::NODE:
            max_error = verify_node_statistics(statistics,
                                               patch_hierarchy,
                                               hier_math_ops,
                                               coeffs,
                                               schedule,
                                               U_idx,
                                               U_var,
                                               time_integrator->getVelocityBoundaryConditions(),
                                               snapshot_time);
            break;
        default:
            TBOX_ERROR("Unsupported analysis centering " << IBAMR::enum_to_string<IBAMR::DataCentering>(data_centering)
                                                         << "\n");
        }
    }

    const bool steady_state = periodic_mode ? statistics->isAtSteadyState() : true;
    const bool success = max_error <= verification_tol && steady_state;
    if (compact_output)
    {
        pout << "mode = " << (periodic_mode ? "PERIODIC_PHASE" : "RUNNING_MEAN") << "\n";
        pout << "analysis_centering = " << IBAMR::enum_to_string<IBAMR::DataCentering>(data_centering) << "\n";
        if (periodic_mode) pout << "periodic_steady_state = " << (steady_state ? "TRUE" : "FALSE") << "\n";
        pout << "verification = " << (success ? "PASS" : "FAIL") << "\n";
    }
    else
    {
        if (periodic_mode) pout << "Periodic steady state flag = " << (steady_state ? "TRUE" : "FALSE") << "\n";
        pout << "Maximum verification error = " << max_error << "\n";
        pout << (success ? "Manufactured turbulence statistics verification passed.\n" :
                           "Manufactured turbulence statistics verification failed.\n");
    }
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
