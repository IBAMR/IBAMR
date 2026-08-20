// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software; you can redistribute it and/or modify it under
// the terms of the 3-clause BSD license. The full text of the license is
// in the file COPYRIGHT at the top level of the IBAMR distribution.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>

#include <petscmat.h>
#include <petscvec.h>

#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
struct KernelSpec
{
    std::string name;
    void (*component_fcn)(double, double*);
    int component_stencil;
    void (*transverse_fcn)(double, double*);
    int transverse_stencil;
};

std::vector<KernelSpec>
get_kernel_specs()
{
    using PMU = IBTK::PETScMatUtilities;
    return {
        { "PIECEWISE_LINEAR",
          PMU::piecewise_linear_delta_fcn,
          PMU::piecewise_linear_delta_stencil,
          PMU::piecewise_linear_delta_fcn,
          PMU::piecewise_linear_delta_stencil },
        { "BSPLINE_3",
          PMU::bspline_3_delta_fcn,
          PMU::bspline_3_delta_stencil,
          PMU::bspline_3_delta_fcn,
          PMU::bspline_3_delta_stencil },
        { "BSPLINE_4",
          PMU::bspline_4_delta_fcn,
          PMU::bspline_4_delta_stencil,
          PMU::bspline_4_delta_fcn,
          PMU::bspline_4_delta_stencil },
        { "BSPLINE_5",
          PMU::bspline_5_delta_fcn,
          PMU::bspline_5_delta_stencil,
          PMU::bspline_5_delta_fcn,
          PMU::bspline_5_delta_stencil },
        { "BSPLINE_6",
          PMU::bspline_6_delta_fcn,
          PMU::bspline_6_delta_stencil,
          PMU::bspline_6_delta_fcn,
          PMU::bspline_6_delta_stencil },
        { "COMPOSITE_BSPLINE_32",
          PMU::bspline_3_delta_fcn,
          PMU::bspline_3_delta_stencil,
          PMU::piecewise_linear_delta_fcn,
          PMU::piecewise_linear_delta_stencil },
        { "COMPOSITE_BSPLINE_43",
          PMU::bspline_4_delta_fcn,
          PMU::bspline_4_delta_stencil,
          PMU::bspline_3_delta_fcn,
          PMU::bspline_3_delta_stencil },
        { "COMPOSITE_BSPLINE_54",
          PMU::bspline_5_delta_fcn,
          PMU::bspline_5_delta_stencil,
          PMU::bspline_4_delta_fcn,
          PMU::bspline_4_delta_stencil },
        { "COMPOSITE_BSPLINE_65",
          PMU::bspline_6_delta_fcn,
          PMU::bspline_6_delta_stencil,
          PMU::bspline_5_delta_fcn,
          PMU::bspline_5_delta_stencil },
        { "IB_3", PMU::ib_3_delta_fcn, PMU::ib_3_delta_stencil, PMU::ib_3_delta_fcn, PMU::ib_3_delta_stencil },
        { "IB_4", PMU::ib_4_delta_fcn, PMU::ib_4_delta_stencil, PMU::ib_4_delta_fcn, PMU::ib_4_delta_stencil },
        { "IB_5", PMU::ib_5_delta_fcn, PMU::ib_5_delta_stencil, PMU::ib_5_delta_fcn, PMU::ib_5_delta_stencil },
        { "IB_6", PMU::ib_6_delta_fcn, PMU::ib_6_delta_stencil, PMU::ib_6_delta_fcn, PMU::ib_6_delta_stencil },
    };
}

bool
is_close(const double actual, const double expected)
{
    const double scale = std::max({ 1.0, std::abs(actual), std::abs(expected) });
    // The fifth-degree B-spline formulas contain substantial cancellation.
    // Keep a fixed near-roundoff bound that remains well below any meaningful
    // coefficient or stencil defect while avoiding compiler-dependent flakes.
    return std::abs(actual - expected) <= 1024.0 * std::numeric_limits<double>::epsilon() * scale;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");

    int test_failures = 0;
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("sc_interp_matrix_live_kernel_consistency is a serial focused test\n");
    }

    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app_initializer);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = std::get<0>(hierarchy_tuple);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);
    const int u_idx = std::get<5>(hierarchy_tuple);
    if (level->getNumberOfPatches() != 1)
    {
        TBOX_ERROR("sc_interp_matrix_live_kernel_consistency requires one patch\n");
    }
    Pointer<Patch<NDIM>> patch = level->getPatch(0);
    Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("sc_interp_matrix_live_kernel_consistency_ctx");
    Pointer<SideVariable<NDIM, int>> dof_var =
        new SideVariable<NDIM, int>("sc_interp_matrix_live_kernel_consistency_dof");
    const int dof_idx = var_db->registerVariableAndContext(dof_var, ctx, IntVector<NDIM>(4));
    level->allocatePatchData(dof_idx);

    std::vector<int> num_dofs_per_proc;
    PETScVecUtilities::constructPatchLevelDOFIndices(num_dofs_per_proc, dof_idx, level);
    Pointer<SideData<NDIM, int>> dof_data = patch->getPatchData(dof_idx);

    std::map<PetscInt, SideIndex<NDIM>> column_to_side;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const SideIndex<NDIM> side_idx(b(), axis, SideIndex<NDIM>::Lower);
            column_to_side[(*dof_data)(side_idx)] = side_idx;
        }
    }

    const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
    const double* const x_lower = patch_geom->getXLower();
    const double* const x_upper = patch_geom->getXUpper();
    const double* const dx = patch_geom->getDx();
    constexpr int n_points = 7;
    std::vector<double> positions(n_points * NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        const double center = 0.5 * (x_lower[d] + x_upper[d]);
        positions[0 * NDIM + d] = center + (0.17 + 0.11 * d) * dx[d];
        positions[1 * NDIM + d] = center - (0.23 + 0.07 * d) * dx[d];
        positions[2 * NDIM + d] = center;
        positions[3 * NDIM + d] = center + 0.5 * dx[d];
        positions[4 * NDIM + d] = x_lower[d] + (0.13 + 0.05 * d) * dx[d];
        positions[5 * NDIM + d] = x_upper[d] - (0.19 + 0.04 * d) * dx[d];
        // Exercise the near-gridline coordinate arithmetic used by the CAV
        // parity case in addition to generic and periodic-boundary points.
        const double trace_fraction = d == 0 ? 0.5 : (d == 1 ? 0.3 : 0.4);
        const double trace_position = x_lower[d] + trace_fraction * (x_upper[d] - x_lower[d]);
        positions[6 * NDIM + d] = d == 0 ? std::nextafter(trace_position, x_lower[d]) : trace_position;
    }

    Vec X_vec = nullptr;
    int ierr = VecCreateMPI(PETSC_COMM_WORLD, static_cast<PetscInt>(positions.size()), PETSC_DETERMINE, &X_vec);
    IBTK_CHKERRQ(ierr);
    std::vector<PetscInt> position_indices(positions.size());
    for (PetscInt k = 0; k < static_cast<PetscInt>(positions.size()); ++k) position_indices[k] = k;
    ierr = VecSetValues(
        X_vec, static_cast<PetscInt>(positions.size()), position_indices.data(), positions.data(), INSERT_VALUES);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(X_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(X_vec);
    IBTK_CHKERRQ(ierr);

    double cell_volume = 1.0;
    for (int d = 0; d < NDIM; ++d) cell_volume *= dx[d];

    pout << std::setprecision(17);
    for (const KernelSpec& kernel : get_kernel_specs())
    {
        Mat J = nullptr;
        PETScMatUtilities::constructPatchLevelSCInterpOp(J,
                                                         kernel.component_fcn,
                                                         kernel.component_stencil,
                                                         kernel.transverse_fcn,
                                                         kernel.transverse_stencil,
                                                         X_vec,
                                                         num_dofs_per_proc,
                                                         dof_idx,
                                                         level);

        double max_interp_error = 0.0;
        double max_spread_error = 0.0;
        bool boundary_wrap_exercised = false;
        bool formula_mutation_detected = false;
        bool weight_order_mutation_detected = false;
        for (PetscInt row = 0; row < n_points * NDIM; ++row)
        {
            PetscInt ncols = 0;
            const PetscInt* cols = nullptr;
            const PetscScalar* vals = nullptr;
            ierr = MatGetRow(J, row, &ncols, &cols, &vals);
            IBTK_CHKERRQ(ierr);
            std::unordered_map<PetscInt, double> matrix_row;
            for (PetscInt k = 0; k < ncols; ++k) matrix_row[cols[k]] = PetscRealPart(vals[k]);
            ierr = MatRestoreRow(J, row, &ncols, &cols, &vals);
            IBTK_CHKERRQ(ierr);

            if (row / NDIM >= 4)
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    bool has_lower_index = false;
                    bool has_upper_index = false;
                    for (const auto& entry : matrix_row)
                    {
                        const auto side_it = column_to_side.find(entry.first);
                        if (side_it == column_to_side.end()) continue;
                        has_lower_index = has_lower_index || side_it->second(d) <= patch->getBox().lower()(d) + 1;
                        has_upper_index = has_upper_index || side_it->second(d) >= patch->getBox().upper()(d) - 1;
                    }
                    boundary_wrap_exercised = boundary_wrap_exercised || (has_lower_index && has_upper_index);
                }
            }

            std::unordered_map<PetscInt, double> live_interp_row;
            for (const auto& entry : matrix_row)
            {
                if (column_to_side.find(entry.first) == column_to_side.end())
                {
                    ++test_failures;
                    continue;
                }
                u_data->fillAll(0.0);
                // Populate every periodic representative of this global DOF
                // so the live kernel sees the same canonical basis vector as
                // the assembled matrix, including across a boundary.
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    const Box<NDIM> ghost_side_box = SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis);
                    for (Box<NDIM>::Iterator b(ghost_side_box); b; b++)
                    {
                        const SideIndex<NDIM> side_idx(b(), axis, SideIndex<NDIM>::Lower);
                        if ((*dof_data)(side_idx) == entry.first) (*u_data)(side_idx) = 1.0;
                    }
                }
                std::vector<double> live_values(n_points * NDIM, std::numeric_limits<double>::quiet_NaN());
                LEInteractor::interpolate(
                    live_values, NDIM, positions, NDIM, u_data, patch, patch->getBox(), kernel.name);
                live_interp_row[entry.first] = live_values[row];
                const double error = std::abs(live_values[row] - entry.second);
                max_interp_error = std::max(max_interp_error, error);
                if (!is_close(live_values[row], entry.second)) ++test_failures;
            }

            // Activate representative formula and stencil-order mutations
            // against the live interpolation row. These controls ensure that
            // the comparison fails for the two defects it is intended to
            // guard against instead of only checking the unmutated path.
            for (const auto& entry : matrix_row)
            {
                if (std::abs(entry.second) > 1.0e-14)
                {
                    const double mutated = entry.second + 1.0e-8 * std::max(1.0, std::abs(entry.second));
                    formula_mutation_detected =
                        formula_mutation_detected || !is_close(mutated, live_interp_row.at(entry.first));
                    break;
                }
            }
            for (auto first = matrix_row.begin(); first != matrix_row.end() && !weight_order_mutation_detected; ++first)
            {
                for (auto second = std::next(first); second != matrix_row.end(); ++second)
                {
                    if (!is_close(first->second, second->second))
                    {
                        const bool first_mismatch = !is_close(second->second, live_interp_row.at(first->first));
                        const bool second_mismatch = !is_close(first->second, live_interp_row.at(second->first));
                        weight_order_mutation_detected = first_mismatch || second_mismatch;
                        break;
                    }
                }
            }

            u_data->fillAll(0.0);
            std::vector<double> lagrangian_values(n_points * NDIM, 0.0);
            lagrangian_values[row] = 1.0;
            LEInteractor::spread(u_data, lagrangian_values, NDIM, positions, NDIM, patch, patch->getBox(), kernel.name);
            std::unordered_map<PetscInt, double> live_spread_row;
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> ghost_side_box = SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis);
                for (Box<NDIM>::Iterator b(ghost_side_box); b; b++)
                {
                    const SideIndex<NDIM> side_idx(b(), axis, SideIndex<NDIM>::Lower);
                    const PetscInt dof = (*dof_data)(side_idx);
                    const double value = cell_volume * (*u_data)(side_idx);
                    if (dof >= 0)
                    {
                        live_spread_row[dof] += value;
                    }
                    else if (!is_close(value, 0.0))
                    {
                        ++test_failures;
                    }
                }
            }
            for (const auto& side_entry : column_to_side)
            {
                const auto matrix_it = matrix_row.find(side_entry.first);
                const double expected = matrix_it == matrix_row.end() ? 0.0 : matrix_it->second;
                const double actual = live_spread_row[side_entry.first];
                const double error = std::abs(actual - expected);
                max_spread_error = std::max(max_spread_error, error);
                if (!is_close(actual, expected)) ++test_failures;
            }
        }
        if (!boundary_wrap_exercised) ++test_failures;
        if (!formula_mutation_detected) ++test_failures;
        if (!weight_order_mutation_detected) ++test_failures;
        pout << kernel.name << " interpolation_max_abs = " << max_interp_error
             << " spreading_max_abs = " << max_spread_error << " boundary_wrap_exercised = " << boundary_wrap_exercised
             << " formula_mutation_detected = " << formula_mutation_detected
             << " weight_order_mutation_detected = " << weight_order_mutation_detected << '\n';
        ierr = MatDestroy(&J);
        IBTK_CHKERRQ(ierr);
    }

    ierr = VecDestroy(&X_vec);
    IBTK_CHKERRQ(ierr);
    level->deallocatePatchData(dof_idx);
    pout << "test_failures = " << test_failures << '\n';
    return test_failures > 0 ? 1 : 0;
}
