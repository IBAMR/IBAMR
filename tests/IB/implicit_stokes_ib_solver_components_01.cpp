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

// Direct interpolation and implicit-operator contracts with a live IBMethod.
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
#include <ibtk/IBKernelEvaluators.h>
#include <ibtk/IBKernelTensorProductEvaluator.h>
#include <ibtk/IBOperatorRegistry.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScMFFDJacobianOperator.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>

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
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

namespace
{
constexpr double epsilon = 1.0 / 1024.0;
// Grid coordinates just below, at, and above cell- and side-centering ties.
const std::array<double, 6> probe = { 4.0 - epsilon, 4.0, 4.0 + epsilon, 4.5 - epsilon, 4.5, 4.5 + epsilon };
// Explicit expected indices, independent of the production rounding formula.
const std::array<int, 6> side_nearest = { 4, 4, 4, 4, 5, 5 };
const std::array<int, 6> cell_nearest = { 3, 4, 4, 4, 4, 4 };
const std::array<int, 6> side_even_lower = { 3, 4, 4, 4, 4, 4 };
const std::array<int, 6> cell_even_lower = { 3, 3, 3, 3, 3, 4 };

// Deliberately asymmetric, distance-dependent weights expose reversal and
// incorrect lower-stencil coordinates as well as misplaced columns.
// Application-defined, move-only evaluator state must survive registration.
struct ProbeEvaluator
{
    std::unique_ptr<double> slope = std::make_unique<double>(1.0 / 32.0);

    std::array<double, 3> operator()(double r) const
    {
        return { 0.125 + r * *slope, 0.375 - r * *slope, 0.5 };
    }
};

template <class Evaluator, std::size_t N>
double
sample_error(const Evaluator& evaluator, double r, const std::array<double, N>& expected)
{
    const std::array<double, N> weights = evaluator(r);
    static_assert(std::tuple_size<decltype(weights)>::value == N, "Natural stencil size changed");
    double error = 0.0;
    for (std::size_t i = 0; i < N; ++i)
    {
        const double entry_error = std::abs(weights[i] - expected[i]);
        if (!(entry_error <= 1.0e-12))
        {
            TBOX_ERROR("Kernel sample error = " << entry_error << '\n');
        }
        error = std::max(error, entry_error);
    }
    return error;
}

template <class Evaluator>
double
moment_error(const Evaluator& evaluator)
{
    constexpr int width = std::tuple_size<decltype(evaluator(0.0))>::value;
    double error = 0.0;
    for (int k = 0; k < 64; ++k)
    {
        const double r = 0.5 * width - 1.0 + k / 64.0;
        const std::array<double, width> w = evaluator(r);
        double sum = 0.0, moment = 0.0;
        for (int i = 0; i < width; ++i)
        {
            sum += w[i];
            moment += i * w[i];
        }
        const double sum_error = std::abs(sum - 1.0), first_moment_error = std::abs(moment - r);
        if (!(sum_error <= 1.0e-12 && first_moment_error <= 1.0e-12))
        {
            TBOX_ERROR("Kernel moment errors = " << sum_error << ", " << first_moment_error << '\n');
        }
        error = std::max({ error, sum_error, first_moment_error });
    }
    return error;
}

int
check_kernels()
{
    const double a = (2.0 - std::sqrt(2.0)) / 8.0, b = (2.0 + std::sqrt(2.0)) / 8.0;
    const double K6 = (59.0 - std::sqrt(261.0)) / 60.0;
    double error = std::max(
        { sample_error(IBKernelEvaluatorBSpline<1>{}, -0.25, std::array<double, 1>{ 1.0 }),
          sample_error(IBKernelEvaluatorBSpline<2>{}, 0.25, std::array<double, 2>{ 0.75, 0.25 }),
          sample_error(IBKernelEvaluatorBSpline<3>{}, 1.0, std::array<double, 3>{ 0.125, 0.75, 0.125 }),
          sample_error(
              IBKernelEvaluatorBSpline<4>{}, 1.5, std::array<double, 4>{ 1.0 / 48, 23.0 / 48, 23.0 / 48, 1.0 / 48 }),
          sample_error(
              IBKernelEvaluatorBSpline<5>{}, 1.5, std::array<double, 5>{ 1.0 / 24, 11.0 / 24, 11.0 / 24, 1.0 / 24, 0 }),
          sample_error(IBKernelEvaluatorBSpline<6>{},
                       2.5,
                       std::array<double, 6>{
                           1.0 / 3840, 237.0 / 3840, 1682.0 / 3840, 1682.0 / 3840, 237.0 / 3840, 1.0 / 3840 }),
          sample_error(IBKernelEvaluatorIB3{}, 1.0, std::array<double, 3>{ 1.0 / 6, 2.0 / 3, 1.0 / 6 }),
          sample_error(IBKernelEvaluatorIB4{}, 1.5, std::array<double, 4>{ a, b, b, a }),
          sample_error(
              IBKernelEvaluatorIB5{},
              1.5,
              std::array<double, 5>{
                  0.0612224005711746881, 0.438777599428825312, 0.438777599428825312, 0.0612224005711746881, 0 }),
          sample_error(
              IBKernelEvaluatorIB6{},
              3.0,
              std::array<double, 6>{ 0, -1.0 / 16 + K6 / 8, 0.25, 5.0 / 8 - K6 / 4, 0.25, -1.0 / 16 + K6 / 8 }) });
    // Exact rational values from the truncated-power definition at r = 11/4.
    error = std::max(error,
                     sample_error(IBKernelEvaluatorBSpline<7>{},
                                  2.75,
                                  std::array<double, 7>{ 729.0 / 2949120,
                                                         112546.0 / 2949120,
                                                         963327.0 / 2949120,
                                                         1434812.0 / 2949120,
                                                         422087.0 / 2949120,
                                                         15618.0 / 2949120,
                                                         1.0 / 2949120 }));
    // Exact rational values from the truncated-power definition at r = 7/2.
    error = std::max(error,
                     sample_error(IBKernelEvaluatorBSpline<8>{},
                                  3.5,
                                  std::array<double, 8>{ 1.0 / 645120,
                                                         2179.0 / 645120,
                                                         60657.0 / 645120,
                                                         259723.0 / 645120,
                                                         259723.0 / 645120,
                                                         60657.0 / 645120,
                                                         2179.0 / 645120,
                                                         1.0 / 645120 }));
    // Independently evaluated Fortran definitions, with natural odd-width
    // coordinates on either side of the nearest-center change.
    error = std::max({ error,
                       sample_error(IBKernelEvaluatorIB5{},
                                    2.25,
                                    std::array<double, 5>{ 0.000539644595320609716,
                                                           0.128737522475479593,
                                                           0.514244366143986938,
                                                           0.333140121904304905,
                                                           0.0233383448809079538 }),
                       sample_error(IBKernelEvaluatorIB5{},
                                    1.75,
                                    std::array<double, 5>{ 0.0233383448809079538,
                                                           0.333140121904304905,
                                                           0.514244366143986938,
                                                           0.128737522475479593,
                                                           0.000539644595320609716 }),
                       sample_error(IBKernelEvaluatorIB6{},
                                    2.25,
                                    std::array<double, 6>{ 0.00965617417165844278,
                                                           0.174648694040214713,
                                                           0.431221688477088836,
                                                           0.325168575099164853,
                                                           0.0591221373512527211,
                                                           0.000182730860620434541 }),
                       sample_error(IBKernelEvaluatorIB6{},
                                    2.75,
                                    std::array<double, 6>{ 0.000182730860620434541,
                                                           0.0591221373512527211,
                                                           0.325168575099164853,
                                                           0.431221688477088836,
                                                           0.174648694040214713,
                                                           0.00965617417165844278 }) });
    const double moments = std::max({ moment_error(IBKernelEvaluatorBSpline<2>{}),
                                      moment_error(IBKernelEvaluatorBSpline<3>{}),
                                      moment_error(IBKernelEvaluatorBSpline<4>{}),
                                      moment_error(IBKernelEvaluatorBSpline<5>{}),
                                      moment_error(IBKernelEvaluatorBSpline<6>{}),
                                      moment_error(IBKernelEvaluatorBSpline<7>{}),
                                      moment_error(IBKernelEvaluatorBSpline<8>{}),
                                      moment_error(IBKernelEvaluatorIB3{}),
                                      moment_error(IBKernelEvaluatorIB4{}),
                                      moment_error(IBKernelEvaluatorIB5{}),
                                      moment_error(IBKernelEvaluatorIB6{}) });

    pout << "kernel_sample_max_error = " << error << '\n';
    pout << "kernel_moment_max_error = " << moments << '\n';
    return error > 1.0e-12 || moments > 1.0e-12;
}

void
register_probe_kernels()
{
    const IBKernel probe_kernel("PROBE");
    const auto register_pair = [&](const IBKernel& other, auto evaluator)
    {
        using E = decltype(evaluator);
        IBOperatorRegistry::register_interpolation_matrix_sc({ probe_kernel, other },
                                                             IBKernelTensorProductEvaluator{ ProbeEvaluator{}, E{} });
        IBOperatorRegistry::register_interpolation_matrix_sc({ other, probe_kernel },
                                                             IBKernelTensorProductEvaluator{ E{}, ProbeEvaluator{} });
    };
    register_pair(IBKernel::BSPLINE_1, IBKernelEvaluatorBSpline<1>{});
    register_pair(IBKernel::BSPLINE_2, IBKernelEvaluatorBSpline<2>{});
    register_pair(IBKernel::IB_4, IBKernelEvaluatorIB4{});
    IBOperatorRegistry::register_interpolation_matrix_sc(
        probe_kernel, IBKernelTensorProductEvaluator{ ProbeEvaluator{}, ProbeEvaluator{} });
}

void
generate_probes(const unsigned int& structure,
                const int& level,
                int& count,
                std::vector<IBTK::Point>& positions,
                void* ctx)
{
    if (ctx && *static_cast<bool*>(ctx))
    {
        count = 1;
        positions.resize(1);
        for (int d = 0; d < NDIM; ++d)
        {
            positions[0](d) = 0.25 / 16.0;
        }
        return;
    }
    TBOX_ASSERT(structure == 0 && level == 0);
    count = probe.size() * probe.size();
    positions.resize(count);
    for (unsigned int j = 0; j < probe.size(); ++j)
    {
        for (unsigned int i = 0; i < probe.size(); ++i)
        {
            positions[j * probe.size() + i](0) = probe[i] / 16.0;
            positions[j * probe.size() + i](1) = probe[j] / 16.0;
        }
    }
}

double
expected_weight(int width, int offset, double distance, bool bspline = false)
{
    if (bspline && width > 2)
    {
        // Independent truncated-power definition of the centered cardinal
        // B-spline, rather than the production knot-interval recurrence.
        double result = 0.0, binomial = 1.0, factorial = 1.0;
        for (int i = 2; i < width; ++i)
        {
            factorial *= i;
        }
        for (int k = 0; k <= width; ++k)
        {
            const double x = std::max(0.0, distance + 0.5 * width - k);
            result += (k % 2 ? -1.0 : 1.0) * binomial * std::pow(x, width - 1) / factorial;
            binomial *= static_cast<double>(width - k) / (k + 1);
        }
        return result;
    }
    if (width == 1)
    {
        return 1.0;
    }
    if (width == 2)
    {
        return std::max(0.0, 1.0 - std::abs(distance));
    }
    if (width == 3)
    {
        const double r_lower = distance + offset;
        return offset == 0 ? 0.125 + r_lower / 32.0 : (offset == 1 ? 0.375 - r_lower / 32.0 : 0.5);
    }
    // The symmetric radial definition of the original even-width IB_4 kernel.
    const double r = std::abs(distance);
    if (r <= 1.0)
    {
        return (3.0 - 2.0 * r + std::sqrt(1.0 + 4.0 * r - 4.0 * r * r)) / 8.0;
    }
    return (5.0 - 2.0 * r - std::sqrt(-7.0 + 12.0 * r - 4.0 * r * r)) / 8.0;
}

bool
check_matrix(Mat matrix,
             Vec positions,
             Pointer<SideData<NDIM, int>> dofs,
             int component_width,
             int transverse_width,
             bool bspline = false)
{
    PetscErrorCode ierr;
    PetscInt begin, end;
    ierr = MatGetOwnershipRange(matrix, &begin, &end);
    IBTK_CHKERRQ(ierr);
    if (begin != 0 || end != static_cast<PetscInt>(NDIM * probe.size() * probe.size()))
    {
        return false;
    }
    const PetscScalar* coordinates;
    ierr = VecGetArrayRead(positions, &coordinates);
    IBTK_CHKERRQ(ierr);
    bool valid = true;
    for (PetscInt row = begin; row < end; ++row)
    {
        const int axis = row % NDIM;
        std::array<int, NDIM> lower, width, sample;
        for (int d = 0; d < NDIM; ++d)
        {
            const double q = 16.0 * PetscRealPart(coordinates[(row - begin) / NDIM * NDIM + d]);
            sample[d] = 0;
            while (sample[d] < static_cast<int>(probe.size()) && q != probe[sample[d]])
            {
                ++sample[d];
            }
            if (sample[d] == static_cast<int>(probe.size()))
            {
                TBOX_ERROR("Unexpected interpolation probe\n");
            }
            width[d] = d == axis ? component_width : transverse_width;
            if (width[d] % 2)
            {
                lower[d] = (d == axis ? side_nearest : cell_nearest)[sample[d]] - width[d] / 2;
            }
            else
            {
                lower[d] = (d == axis ? side_even_lower : cell_even_lower)[sample[d]] - (width[d] / 2 - 1);
            }
        }
        std::map<PetscInt, double> expected;
        for (int j = 0; j < width[1]; ++j)
        {
            for (int i = 0; i < width[0]; ++i)
            {
                SAMRAI::hier::Index<NDIM> index;
                index(0) = lower[0] + i;
                index(1) = lower[1] + j;
                const SideIndex<NDIM> side(index, axis, SideIndex<NDIM>::Lower);
                const double x = probe[sample[0]] - (index(0) + (axis == 0 ? 0.0 : 0.5));
                const double y = probe[sample[1]] - (index(1) + (axis == 1 ? 0.0 : 0.5));
                const int column = (*dofs)(side);
                valid = valid && column >= 0;
                expected[column] = expected_weight(width[0], i, x, bspline) * expected_weight(width[1], j, y, bspline);
            }
        }
        PetscInt count;
        const PetscInt* columns;
        const PetscScalar* values;
        ierr = MatGetRow(matrix, row, &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
        valid = valid && count == static_cast<PetscInt>(expected.size());
        for (PetscInt k = 0; k < count; ++k)
        {
            const auto found = expected.find(columns[k]);
            valid = valid && found != expected.end() && std::abs(PetscRealPart(values[k]) - found->second) < 1.0e-12;
        }
        ierr = MatRestoreRow(matrix, row, &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(positions, &coordinates);
    IBTK_CHKERRQ(ierr);
    return valid;
}

int
check_periodic_interpolation(const std::vector<int>& counts, int dof, Pointer<PatchLevel<NDIM>> level)
{
    // Use only interior DOF entries for the reference; assembly reads ghost entries.
    Pointer<SideData<NDIM, int>> dofs = level->getPatch(0)->getPatchData(dof);
    constexpr int cells = 16, points = 4;
    const std::array<double, 2> locations = { 0.25, 15.75 };
    Vec X = nullptr;
    PetscErrorCode ierr = VecCreateSeq(PETSC_COMM_SELF, NDIM * points, &X);
    IBTK_CHKERRQ(ierr);
    PetscScalar* coordinates;
    ierr = VecGetArray(X, &coordinates);
    IBTK_CHKERRQ(ierr);
    for (int point = 0; point < points; ++point)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            coordinates[NDIM * point + d] = locations[(point >> d) & 1] / cells;
        }
    }
    ierr = VecRestoreArray(X, &coordinates);
    IBTK_CHKERRQ(ierr);
    double weight_error = 0.0, action_error = 0.0;
    int column_mismatches = 0;
    for (bool normal_odd : { true, false })
    {
        Mat matrix = nullptr;
        if (normal_odd)
        {
            PETScMatUtilities::constructPatchLevelSCInterpOp(
                matrix,
                IBKernelTensorProductEvaluator{ IBKernelEvaluatorBSpline<3>{}, IBKernelEvaluatorBSpline<2>{} },
                X,
                counts,
                dof,
                level);
        }
        else
        {
            PETScMatUtilities::constructPatchLevelSCInterpOp(
                matrix,
                IBKernelTensorProductEvaluator{ IBKernelEvaluatorBSpline<2>{}, IBKernelEvaluatorBSpline<3>{} },
                X,
                counts,
                dof,
                level);
        }
        Vec field = nullptr, result = nullptr;
        ierr = MatCreateVecs(matrix, &field, &result);
        IBTK_CHKERRQ(ierr);
        ierr = VecSet(field, 0.0);
        IBTK_CHKERRQ(ierr);
        const auto field_value = [](int axis, int i, int j) { return 1.0 + axis + 0.125 * i + 0.03125 * j * j; };
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (int j = 0; j < cells; ++j)
            {
                for (int i = 0; i < cells; ++i)
                {
                    hier::Index<NDIM> index;
                    index(0) = i;
                    index(1) = j;
                    const int column = (*dofs)(SideIndex<NDIM>(index, axis, SideIndex<NDIM>::Lower));
                    if (column < 0)
                    {
                        TBOX_ERROR("Missing interior periodic DOF\n");
                    }
                    ierr = VecSetValue(field, column, field_value(axis, i, j), INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }
        }
        ierr = VecAssemblyBegin(field);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(field);
        IBTK_CHKERRQ(ierr);
        ierr = MatMult(matrix, field, result);
        IBTK_CHKERRQ(ierr);
        const PetscScalar* actual;
        ierr = VecGetArrayRead(result, &actual);
        IBTK_CHKERRQ(ierr);
        for (int point = 0; point < points; ++point)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                std::array<int, NDIM> widths, lower;
                std::array<double, NDIM> q;
                for (int d = 0; d < NDIM; ++d)
                {
                    widths[d] = (d == axis) == normal_odd ? 3 : 2;
                    q[d] = locations[(point >> d) & 1] - (d == axis ? 0.0 : 0.5);
                    lower[d] = widths[d] == 3 ? static_cast<int>(std::floor(q[d] + 0.5)) - 1 :
                                                static_cast<int>(std::floor(q[d]));
                }
                std::map<PetscInt, double> expected;
                double expected_action = 0.0;
                for (int j = 0; j < widths[1]; ++j)
                {
                    for (int i = 0; i < widths[0]; ++i)
                    {
                        const int ix = lower[0] + i, iy = lower[1] + j;
                        hier::Index<NDIM> wrapped;
                        wrapped(0) = (ix % cells + cells) % cells;
                        wrapped(1) = (iy % cells + cells) % cells;
                        const int column = (*dofs)(SideIndex<NDIM>(wrapped, axis, SideIndex<NDIM>::Lower));
                        const double weight = expected_weight(widths[0], i, q[0] - ix, true) *
                                              expected_weight(widths[1], j, q[1] - iy, true);
                        expected[column] += weight;
                        expected_action += weight * field_value(axis, wrapped(0), wrapped(1));
                    }
                }
                const int row = NDIM * point + axis;
                PetscInt count;
                const PetscInt* columns;
                const PetscScalar* values;
                ierr = MatGetRow(matrix, row, &count, &columns, &values);
                IBTK_CHKERRQ(ierr);
                column_mismatches += count != static_cast<PetscInt>(expected.size());
                for (PetscInt j = 0; j < count; ++j)
                {
                    const auto found = expected.find(columns[j]);
                    if (found == expected.end())
                    {
                        ++column_mismatches;
                    }
                    else
                    {
                        const double entry_error = std::abs(PetscRealPart(values[j]) - found->second);
                        if (!(entry_error <= 1.0e-12))
                        {
                            TBOX_ERROR("Periodic weight error = " << entry_error << '\n');
                        }
                        weight_error = std::max(weight_error, entry_error);
                    }
                }
                ierr = MatRestoreRow(matrix, row, &count, &columns, &values);
                IBTK_CHKERRQ(ierr);
                const double row_error = std::abs(PetscRealPart(actual[row]) - expected_action);
                if (!(row_error <= 1.0e-12))
                {
                    TBOX_ERROR("Periodic interpolation error = " << row_error << '\n');
                }
                action_error = std::max(action_error, row_error);
            }
        }
        ierr = VecRestoreArrayRead(result, &actual);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&field);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&result);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&matrix);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecDestroy(&X);
    IBTK_CHKERRQ(ierr);
    pout << "periodic_column_mismatches = " << column_mismatches << '\n'
         << "periodic_weight_max_error = " << weight_error << '\n'
         << "periodic_action_max_error = " << action_error << '\n';
    return column_mismatches != 0 || weight_error > 1.0e-12 || action_error > 1.0e-12;
}
} // namespace

int
run_interpolation(Pointer<AppInitializer> app, const std::string& input_file)
{
    PetscErrorCode ierr;
    int failures = 0;
    {
        const std::string setup_test = app->getInputDatabase()->getStringWithDefault("setup_test", "");
        bool setup_probe = !setup_test.empty();
        const bool late_fixed = setup_test == "late_fixed";
        if (setup_probe)
        {
            Pointer<Logger::Appender> appender = new TestAppender();
            Logger::getInstance()->setAbortAppender(appender);
        }
        constexpr int max_bspline_order = IBTK_MAX_BSPLINE_ORDER;
        const bool unsupported = input_file.find("registration.unsupported") != std::string::npos;
        if (!unsupported && !setup_probe)
        {
            failures += check_kernels();
        }
        Pointer<IBMethod> method = new IBMethod("IBMethod", app->getComponentDatabase("IBMethod"));
        method->setUseFixedLEOperators(!late_fixed);
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
        initializer->setStructureNamesOnLevel(0, { "probes" });
        initializer->registerInitStructureFunction(generate_probes, &setup_probe);
        method->registerLInitStrategy(initializer);
        gridding->makeCoarsestLevel(hierarchy, 0.0);
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        // Interior probes check centering; setup probes exercise periodic ghost support.
        if (IBTK_MPI::getNodes() != 1 || level->getNumberOfPatches() != 1)
        {
            TBOX_ERROR("Interpolation fixture requires one patch on one rank\n");
        }
        VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = variables->getContext("interpolation");
        Pointer<SideVariable<NDIM, double>> velocity = new SideVariable<NDIM, double>("velocity");
        Pointer<SideVariable<NDIM, int>> indices = new SideVariable<NDIM, int>("indices");
        const int u = variables->registerVariableAndContext(velocity, context, method->getMinimumGhostCellWidth());
        // Match the implicit integrator's DOF allocation, including IBMethod's input minimum.
        const int dof = variables->registerVariableAndContext(indices, context, method->getMinimumGhostCellWidth());
        level->allocatePatchData(u, 0.0);
        level->allocatePatchData(dof, 0.0);
        Pointer<SideData<NDIM, double>> u_data = level->getPatch(0)->getPatchData(u);
        u_data->fillAll(0.0);
        std::vector<int> counts;
        PETScVecUtilities::constructPatchLevelDOFIndices(counts, dof, level);
        Pointer<SideData<NDIM, int>> dofs = level->getPatch(0)->getPatchData(dof);
        std::vector<Pointer<CoarsenSchedule<NDIM>>> synch(1);
        std::vector<Pointer<RefineSchedule<NDIM>>> ghost_fill(1);
        method->initializePatchHierarchy(hierarchy, gridding, u, synch, ghost_fill, 0, 0.0, true);
        method->freeLInitStrategy();
        initializer.setNull();

        if (setup_probe)
        {
            method->preprocessIntegrateData(0.0, 0.125, 1);
            if (late_fixed)
            {
                method->setUseFixedLEOperators(true);
            }
            method->updateFixedLEOperators();
            if (setup_test == "wide_kernel")
            {
                // Applications may still register this evaluator with a smaller built-in catalog.
                if (max_bspline_order < 8)
                {
                    IBOperatorRegistry::register_interpolation_matrix_sc(
                        IBKernel("BSPLINE_8"), IBKernelTensorProductEvaluator{ IBKernelEvaluatorBSpline<8>{} });
                }
                Mat matrix = nullptr;
                method->constructInterpOp(matrix, IBKernel("BSPLINE_8"), counts, dof, 0.125);
                Vec ones = nullptr, residual = nullptr;
                ierr = MatCreateVecs(matrix, &ones, &residual);
                IBTK_CHKERRQ(ierr);
                ierr = VecSet(ones, 1.0);
                IBTK_CHKERRQ(ierr);
                ierr = MatMult(matrix, ones, residual);
                IBTK_CHKERRQ(ierr);
                ierr = VecShift(residual, -1.0);
                IBTK_CHKERRQ(ierr);
                PetscReal error = 0.0;
                ierr = VecNorm(residual, NORM_INFINITY, &error);
                IBTK_CHKERRQ(ierr);
                failures += !(error <= 1.0e-12);
                pout << "dof_ghost_width = " << method->getMinimumGhostCellWidth()(0) << '\n'
                     << "constant_interpolation_error = " << error << '\n';
                ierr = VecDestroy(&ones);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&residual);
                IBTK_CHKERRQ(ierr);
                ierr = MatDestroy(&matrix);
                IBTK_CHKERRQ(ierr);
            }
            method->postprocessIntegrateData(0.0, 0.125, 1);
            method->postprocessData();
            level->deallocatePatchData(u);
            level->deallocatePatchData(dof);
            variables->removePatchDataIndex(u);
            variables->removePatchDataIndex(dof);
            pout << "test_failures = " << failures << std::endl;
            return failures;
        }

        if (unsupported)
        {
            method->preprocessIntegrateData(0.0, 0.125, 1);
            method->updateFixedLEOperators();
            Mat matrix = nullptr;
            method->constructInterpOp(
                matrix, IBKernel(app->getInputDatabase()->getString("matrix_kernel")), counts, dof, 0.125);
            ierr = MatDestroy(&matrix);
            IBTK_CHKERRQ(ierr);
            method->postprocessIntegrateData(0.0, 0.125, 1);
            method->postprocessData();
        }

        // A configured bank does not reserve higher-order kernel specifications.
        IBOperatorRegistry::register_interpolation_matrix_sc(
            IBKernel("BSPLINE_" + std::to_string(IBTK_MAX_BSPLINE_ORDER + 1)),
            IBKernelTensorProductEvaluator{ IBKernelEvaluatorBSpline<IBTK_MAX_BSPLINE_ORDER + 1>{} });
        failures += check_periodic_interpolation(counts, dof, level);

        const std::array<IBKernel, 4> kernel = {
            IBKernel::BSPLINE_1, IBKernel::BSPLINE_2, IBKernel("PROBE"), IBKernel::IB_4
        };
        bool placement = true, scalar_equivalence = true, lifecycle = true;
        // Repeat setup/use/cleanup on the same IBMethod to exercise scratch state
        // invalidation, with fixed-operator updates at both new and half times.
        for (int step = 0; step < 2; ++step)
        {
            const double current_time = step * 0.125, new_time = current_time + 0.125;
            method->preprocessIntegrateData(current_time, new_time, 1);
            method->updateFixedLEOperators();
            Vec X = method->getLDataManager()->getLData("X", 0)->getVec();
            // Direct construction precedes registration and uses distinct functions
            // with the same widths. Both paths must retain the actual weights.
            Mat direct = nullptr, builtin = nullptr, registered = nullptr;
            PETScMatUtilities::constructPatchLevelSCInterpOp(
                direct,
                IBKernelTensorProductEvaluator{ ProbeEvaluator{}, IBKernelEvaluatorBSpline<2>{} },
                X,
                counts,
                dof,
                level);
            PETScMatUtilities::constructPatchLevelSCInterpOp(
                builtin,
                IBKernelTensorProductEvaluator{ IBKernelEvaluatorBSpline<3>{}, IBKernelEvaluatorBSpline<2>{} },
                X,
                counts,
                dof,
                level);
            placement = check_matrix(direct, X, dofs, 3, 2) && placement;
            placement = check_matrix(builtin, X, dofs, 3, 2, true) && placement;
            PetscBool equal;
            ierr = MatEqual(direct, builtin, &equal);
            IBTK_CHKERRQ(ierr);
            if (equal)
            {
                ++failures;
            }
            if (step == 0)
            {
                register_probe_kernels();
            }
            method->constructInterpOp(registered, { IBKernel("PROBE"), IBKernel::BSPLINE_2 }, counts, dof, new_time);
            ierr = MatEqual(direct, registered, &equal);
            IBTK_CHKERRQ(ierr);
            if (!equal)
            {
                ++failures;
            }
            if (max_bspline_order >= 3)
            {
                method->constructInterpOp(
                    registered, { IBKernel::BSPLINE_3, IBKernel::BSPLINE_2 }, counts, dof, new_time);
                ierr = MatEqual(builtin, registered, &equal);
                IBTK_CHKERRQ(ierr);
                if (!equal)
                {
                    ++failures;
                }
            }

            // Direct evaluation is independent of the configured registrations.
            // Supply this combination when it is outside the automatic catalog.
            PETScMatUtilities::constructPatchLevelSCInterpOp(
                direct,
                IBKernelTensorProductEvaluator{ IBKernelEvaluatorBSpline<7>{}, IBKernelEvaluatorBSpline<2>{} },
                X,
                counts,
                dof,
                level);
            placement = check_matrix(direct, X, dofs, 7, 2, true) && placement;
            if (step == 0 && max_bspline_order < 7)
            {
                IBOperatorRegistry::register_interpolation_matrix_sc(
                    { IBKernel("BSPLINE_7"), IBKernel::BSPLINE_2 },
                    IBKernelTensorProductEvaluator{ IBKernelEvaluatorBSpline<7>{}, IBKernelEvaluatorBSpline<2>{} });
            }
            method->constructInterpOp(
                registered, { IBKernel("BSPLINE_7"), IBKernel::BSPLINE_2 }, counts, dof, new_time);
            ierr = MatEqual(direct, registered, &equal);
            IBTK_CHKERRQ(ierr);
            if (!equal)
            {
                ++failures;
            }

            // The B-spline limit does not restrict the other supplied kernels.
            PETScMatUtilities::constructPatchLevelSCInterpOp(
                direct, IBKernelTensorProductEvaluator{ IBKernelEvaluatorIB6{} }, X, counts, dof, level);
            method->constructInterpOp(registered, IBKernel::IB_6, counts, dof, new_time);
            ierr = MatEqual(direct, registered, &equal);
            IBTK_CHKERRQ(ierr);
            if (!equal)
            {
                ++failures;
            }
            ierr = MatDestroy(&direct);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&builtin);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&registered);
            IBTK_CHKERRQ(ierr);
            for (int cw = 1; cw <= 4; ++cw)
            {
                for (int tw = 1; tw <= 4; ++tw)
                {
                    if (max_bspline_order < 2 && (cw == 2 || tw == 2))
                    {
                        continue;
                    }
                    Mat matrix = nullptr;
                    IBImplicitStrategy& strategy = *method;
                    strategy.constructInterpOp(matrix, { kernel[cw - 1], kernel[tw - 1] }, counts, dof, new_time);
                    placement = check_matrix(matrix, X, dofs, cw, tw) && placement;
                    if (cw == tw)
                    {
                        Mat scalar = nullptr;
                        strategy.constructInterpOp(scalar, kernel[cw - 1], counts, dof, new_time);
                        PetscBool equal;
                        ierr = MatEqual(matrix, scalar, &equal);
                        IBTK_CHKERRQ(ierr);
                        scalar_equivalence = scalar_equivalence && equal;
                        IBOperatorRegistry::construct_interpolation_matrix_sc(
                            scalar, kernel[cw - 1], X, counts, dof, level);
                        ierr = MatEqual(matrix, scalar, &equal);
                        IBTK_CHKERRQ(ierr);
                        scalar_equivalence = scalar_equivalence && equal;
                        ierr = MatDestroy(&scalar);
                        IBTK_CHKERRQ(ierr);
                    }
                    ierr = MatDestroy(&matrix);
                    IBTK_CHKERRQ(ierr);
                }
            }
            // Keep the highest-order sample inside this fixed patch geometry.
            const int highest_order = std::min(max_bspline_order, 8);
            const int transverse_order = std::min(max_bspline_order, 2);
            for (bool normal : { true, false })
            {
                const int cw = normal ? highest_order : transverse_order;
                const int tw = normal ? transverse_order : highest_order;
                Mat matrix = nullptr;
                method->constructInterpOp(
                    matrix,
                    { IBKernel("BSPLINE_" + std::to_string(cw)), IBKernel("BSPLINE_" + std::to_string(tw)) },
                    counts,
                    dof,
                    new_time);
                placement = check_matrix(matrix, X, dofs, cw, tw, true) && placement;
                ierr = MatDestroy(&matrix);
                IBTK_CHKERRQ(ierr);
            }
            // Natural named odd widths in both orientations, including mixed
            // parity, checked against independent scalar B-spline values.
            for (int cw : { 2, 3, 5 })
            {
                for (int tw : { 2, 3, 5 })
                {
                    if (cw > max_bspline_order || tw > max_bspline_order || (cw == 2 && tw == 2))
                    {
                        continue;
                    }
                    Mat matrix = nullptr;
                    method->constructInterpOp(
                        matrix,
                        { IBKernel("BSPLINE_" + std::to_string(cw)), IBKernel("BSPLINE_" + std::to_string(tw)) },
                        counts,
                        dof,
                        new_time);
                    placement = check_matrix(matrix, X, dofs, cw, tw, true) && placement;
                    ierr = MatDestroy(&matrix);
                    IBTK_CHKERRQ(ierr);
                }
            }
            Mat half = nullptr;
            method->constructInterpOp(half, kernel[0], counts, dof, current_time + 0.0625);
            lifecycle = check_matrix(half, X, dofs, 1, 1) && lifecycle;
            ierr = MatDestroy(&half);
            IBTK_CHKERRQ(ierr);
            method->postprocessIntegrateData(current_time, new_time, 1);
            method->postprocessData();
        }
        failures += !placement + !scalar_equivalence + !lifecycle;
        pout << "stencil_columns_and_weights_valid = " << (placement ? "true" : "false") << '\n';
        pout << "isotropic_equivalence_valid = " << (scalar_equivalence ? "true" : "false") << '\n';
        pout << "standalone_lifecycle_valid = " << (lifecycle ? "true" : "false") << '\n';
        level->deallocatePatchData(u);
        level->deallocatePatchData(dof);
        variables->removePatchDataIndex(u);
        variables->removePatchDataIndex(dof);
        pout << "test_failures = " << failures << std::endl;
    }
    return failures;
}

namespace
{
using HierarchyVector = SAMRAIVectorReal<NDIM, double>;

// Count evaluations and boundary dispatch without changing the Stokes operations.
class BoundaryCheckedStokesOperator : public StaggeredStokesOperator
{
public:
    BoundaryCheckedStokesOperator() : StaggeredStokesOperator("operator_test::stokes")
    {
    }
    void apply(HierarchyVector& x, HierarchyVector& y) override
    {
        ++evaluations;
        StaggeredStokesOperator::apply(x, y);
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
        if (edge.first > edge.second) std::swap(edge.first, edge.second);
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

int
run_operators(Pointer<AppInitializer> app)
{
    PetscErrorCode ierr;
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
    initializer->setStructureNamesOnLevel(0, { "curve" });
    initializer->registerInitStructureFunction(generate_operator_structure, &physical_boundary);
    initializer->registerInitSpringDataFunction(generate_operator_springs);
    method->registerLInitStrategy(initializer);
    gridding->makeCoarsestLevel(hierarchy, current);
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    if (IBTK_MPI::getNodes() != 1 || level->getNumberOfPatches() != 1)
        TBOX_ERROR("Operator fixture requires one patch on one rank\n");

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
    std::vector<Pointer<CoarsenSchedule<NDIM>>> synch(1);
    std::vector<Pointer<RefineSchedule<NDIM>>> fill(1), prolong(1);
    std::vector<Pointer<LocationIndexRobinBcCoefs<NDIM>>> physical_coefs(NDIM);
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM, nullptr);
    Pointer<CartSideRobinPhysBdryOp> physical_bc;
    if (physical_boundary)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            physical_coefs[axis] = new LocationIndexRobinBcCoefs<NDIM>("velocity_bc", nullptr);
            for (int face = 0; face < 2 * NDIM; ++face)
                physical_coefs[axis]->setBoundaryValue(face, axis == 1 && face < 2 ? 0.3 : 0.0);
            bc_coefs[axis] = physical_coefs[axis];
        }
        physical_bc = new CartSideRobinPhysBdryOp(u_current, bc_coefs, false);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            physical_bc->setPhysicalBoundaryConditions(*level->getPatch(p()), current, ghosts);
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
        fill[0] = ghost_fill.createSchedule(level, physical_bc.getPointer());
    }

    Pointer<HierarchySideDataOpsReal<NDIM, double>> side_ops =
        new HierarchySideDataOpsReal<NDIM, double>(hierarchy, 0, 0);
    Pointer<HierarchyCellDataOpsReal<NDIM, double>> cell_ops =
        new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, 0, 0);
    HierarchyMathOps math_ops("operator_test::math", hierarchy);
    const int u = variables->registerVariableAndContext(u_var, context, ghosts);
    const int p = variables->registerVariableAndContext(p_var, context, IntVector<NDIM>(1));
    Pointer<HierarchyVector> base = new HierarchyVector("base", hierarchy, 0, 0);
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
        Vec X0 = method->getLDataManager()->getLData("X", 0)->getVec();
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
                ierr = VecCopy((*positions)[0]->getVec(), physical_position);
                IBTK_CHKERRQ(ierr);
                // A zero-boundary control proves the marker interpolation sees the wall data.
                for (int face = 0; face < 2; ++face) physical_coefs[1]->setBoundaryValue(face, 0.0);
                nonlinear.apply(*base, *expected);
                ierr = VecAXPY(physical_position, -1.0, (*positions)[0]->getVec());
                IBTK_CHKERRQ(ierr);
                PetscReal position_change = 0.0;
                ierr = VecNorm(physical_position, NORM_INFINITY, &position_change);
                IBTK_CHKERRQ(ierr);
                boundary_position_change = std::max(boundary_position_change, position_change);
                for (int face = 0; face < 2; ++face) physical_coefs[1]->setBoundaryValue(face, 0.3);
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
        for (Pointer<HierarchyVector>& vector : vectors) free_vector_components(*vector);
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
    method->constructInterpOp(J, IBKernel::IB_4, counts, u_dof, force_time);
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
    Vec X0 = method->getLDataManager()->getLData("X", 0)->getVec();
    Pointer<LData> expected_force = new LData("expected_force", 16, NDIM);
    const double cell_volume = geometry->getDx()[0] * geometry->getDx()[1];
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
    ctx.patch_level = level;
    ctx.u_idx = scratch;
    ctx.f_idx = f_scratch;
    ctx.u_current_idx = u_current;
    ctx.u_dof_index_idx = u_dof;
    ctx.p_dof_index_idx = p_dof;
    ctx.time_stepping_type = type;
    StaggeredStokesIBOperator nonlinear("nonlinear");
    StaggeredStokesIBJacobianOperator jacobian("jacobian");
    PETScMFFDJacobianOperator mffd("mffd");
    nonlinear.setOperatorContext(ctx);
    jacobian.setOperatorContext(ctx);
    mffd.setOperator(stokes);
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
            pout << label << ": error = " << std::setprecision(17) << error << ", bound = " << bound
                 << std::setprecision(6) << '\n';
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
        mffd.initializeOperatorState(*mffd_base, *residual);
        method->constructInterpOp(J, IBKernel::IB_4, counts, u_dof, force_time);
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
                ierr = VecAXPY(position, dt / 2, (*U_current_data)[0]->getVec());
                IBTK_CHKERRQ(ierr);
            }
            ierr = VecCopy(position, expected_position);
            IBTK_CHKERRQ(ierr);
            std::vector<Pointer<LData>>*X_data, *U_data;
            bool* X_ghost;
            const TimePoint force_point = midpoint ? TimePoint::HALF_TIME : TimePoint::NEW_TIME;
            method->getPositionData(&X_data, &X_ghost, force_point);
            method->getVelocityData(&U_data, force_point);
            ierr = VecAXPY(position, -1.0, (*X_data)[0]->getVec());
            IBTK_CHKERRQ(ierr);
            PetscReal position_error;
            ierr = VecNorm(position, NORM_INFINITY, &position_error);
            IBTK_CHKERRQ(ierr);
            time_valid = time_valid && std::isfinite(position_error) && position_error < 1.0e-12;
            stokes->apply(*base, *expected);
            ierr = VecSet(expected_force->getVec(), 0.0);
            IBTK_CHKERRQ(ierr);
            force->computeLagrangianForce(
                expected_force, (*X_data)[0], (*U_data)[0], hierarchy, 0, force_time, method->getLDataManager());
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
                                                  (*X_data)[0],
                                                  0.0,
                                                  nullptr,
                                                  hierarchy,
                                                  0,
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
                ierr = VecWAXPY(position, -1.0, expected_position, (*X_data)[0]->getVec());
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
                first_action->copyVector(action);
            else
                nontrivial = nontrivial && !close(action, first_action, 1.0e-7);
            jacobian.setIBCouplingJacobian(coupling);
            // Release the creator's reference, then reinstall the borrowed handle.
            // The operator must retain its own reference throughout replacement.
            Mat coupling_alias = coupling;
            ierr = MatDestroy(&coupling);
            IBTK_CHKERRQ(ierr);
            jacobian.setIBCouplingJacobian(coupling_alias);
            jacobian.apply(*direction, *expected);
            assembled_valid = close(action, expected, JACOBIAN_TOL, "assembled Jacobian") && assembled_valid;
            Mat no_coupling = nullptr;
            jacobian.setIBCouplingJacobian(no_coupling);

            const double h = 1.0e-5;
            plus->linearSum(1.0, base, h, direction);
            minus->linearSum(1.0, base, -h, direction);
            nonlinear.apply(*plus, *expected);
            nonlinear.apply(*minus, *work);
            finite_difference->linearSum(0.5 / h, expected, -0.5 / h, work);
            derivative_valid =
                close(action, finite_difference, FD_TOL, "centered finite difference") && derivative_valid;
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
            mffd_valid = close(expected, work, FD_TOL, "MFFD Stokes action") && work->maxNorm() > 1.0e-12 && mffd_valid;

            nonlinear.applyAdd(*base, *direction, *expected);
            work->add(residual, direction);
            nonlinear_add_valid = close(expected, work, RESIDUAL_TOL, "nonlinear applyAdd") && nonlinear_add_valid;
            jacobian.applyAdd(*direction, *base, *expected);
            work->add(action, base);
            jacobian_add_valid = close(expected, work, JACOBIAN_TOL, "Jacobian applyAdd") && jacobian_add_valid;
        }
        for (GeneralOperator* op :
             { static_cast<GeneralOperator*>(&nonlinear), static_cast<GeneralOperator*>(&jacobian) })
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
        const std::array<int, 4> base_indices = { jacobian.getBaseVector()->getComponentDescriptorIndex(0),
                                                  jacobian.getBaseVector()->getComponentDescriptorIndex(1),
                                                  mffd.getBaseVector()->getComponentDescriptorIndex(0),
                                                  mffd.getBaseVector()->getComponentDescriptorIndex(1) };
        mffd.deallocateOperatorState();
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
    Mat no_coupling = nullptr;
    jacobian.setIBCouplingJacobian(no_coupling);
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
    for (const std::tuple<const char*, bool, double>& check :
         { std::make_tuple("nonlinear_residual", residual_valid, RESIDUAL_TOL),
           std::make_tuple("assembled_jacobian", assembled_valid, JACOBIAN_TOL),
           std::make_tuple("centered_fd", derivative_valid, FD_TOL),
           std::make_tuple("nonlinear_apply_add", nonlinear_add_valid, RESIDUAL_TOL),
           std::make_tuple("jacobian_apply_add", jacobian_add_valid, JACOBIAN_TOL),
           std::make_tuple("mffd_stokes_action", mffd_valid, FD_TOL),
           std::make_tuple("initial_supplied_action", initial_supplied_valid, INITIALIZATION_TOL),
           std::make_tuple("outer_initialization_strategy", strategy_valid, INITIALIZATION_TOL),
           std::make_tuple("post_initialization_supplied", supplied_valid, INITIALIZATION_TOL) })
    {
        pout << std::get<0>(check) << " = " << (std::get<1>(check) ? "true" : "false")
             << ", tolerance = " << std::get<2>(check) << '\n';
        failures += !std::get<1>(check);
    }
    for (const std::pair<std::string, bool>& check :
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
    for (Pointer<HierarchyVector>& vector : vectors) free_vector_components(*vector);
    for (int idx : allocated)
    {
        level->deallocatePatchData(idx);
        variables->removePatchDataIndex(idx);
    }
    pout << "test_failures = " << failures << std::endl;
    return failures;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    // Keep optional visualization warnings out of the compared test output.
    Logger::getInstance()->setWarning(false);
    const std::string input_file = argc > 1 ? argv[1] : "";
    const bool duplicate_custom = input_file.find("registration.duplicate_custom") != std::string::npos;
    if (input_file.find("registration.") != std::string::npos)
    {
        Pointer<Logger::Appender> appender = new TestAppender();
        Logger::getInstance()->setAbortAppender(appender);
        PIO::logOnlyNodeZero("output");
        if (input_file.find("duplicate") != std::string::npos ||
            input_file.find("registration.unknown") != std::string::npos)
        {
            std::ifstream input(input_file);
            std::string kernel_name, transverse_name;
            input >> kernel_name;
            const IBKernelTensorProduct kernel =
                input >> transverse_name ? IBKernelTensorProduct{ IBKernel(kernel_name), IBKernel(transverse_name) } :
                                           IBKernelTensorProduct{ IBKernel(kernel_name) };
            if (duplicate_custom)
            {
                IBOperatorRegistry::register_interpolation_matrix_sc(
                    kernel, IBKernelTensorProductEvaluator{ IBKernelEvaluatorIB4{} });
            }
            IBOperatorRegistry::register_interpolation_matrix_sc(
                kernel, IBKernelTensorProductEvaluator{ IBKernelEvaluatorIB4{} });
            return 0;
        }
    }
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "components.log");
    const std::string test_case = app->getInputDatabase()->getStringWithDefault("test_case", "interpolation");
    if (test_case == "interpolation") return run_interpolation(app, input_file);
    if (test_case == "operators") return run_operators(app);
    TBOX_ERROR("Unknown component test case: " << test_case << '\n');
    return 1;
}
