// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/IBOperatorBuilder.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/ib_kernel_dispatch.h>
#include <ibtk/ib_kernel_evaluators.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/interpolation_utilities.h>

#include <tbox/Logger.h>
#include <tbox/Utilities.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <array>
#include <concepts>
#include <iomanip>
#include <limits>
#include <map>
#include <string>
#include <utility>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

struct LinearIBKernel
{
    static constexpr std::size_t get_stencil_width()
    {
        return 2;
    }
    template <class Output, std::floating_point Input>
    requires IBTK::detail::IBKernelWritableWeights<Output, 2> Output evaluate(const Input& r) const
    {
        using Coefficient = ib_kernel_weights_value_t<Output>;
        const Coefficient x = r;
        Output weights{};
        weights[0] = 1 - x;
        weights[1] = x;
        return weights;
    }
    template <class Output, class Input>
    Output evaluate(Input&) const = delete;
    template <class Output, class Input>
    Output evaluate(const Input&&) const = delete;
};

static_assert(IBKernelEvaluatorScalar<LinearIBKernel>);

struct ExplicitMoveEvaluator
{
    ExplicitMoveEvaluator() = default;
    explicit ExplicitMoveEvaluator(ExplicitMoveEvaluator&&) = default;
    ExplicitMoveEvaluator(const ExplicitMoveEvaluator&) = delete;

    template <int Axis>
    static constexpr std::array<std::size_t, NDIM> get_stencil_widths()
    {
        std::array<std::size_t, NDIM> widths;
        widths.fill(1);
        return widths;
    }
    template <int Axis, IBKernelWeights Output, std::floating_point Input>
    requires(Axis >= 0 && Axis < NDIM && ib_kernel_weights_extent_v<Output> == 1) Output
        evaluate(const std::array<Input, NDIM>&) const
    {
        return Output{ 1 };
    }
};

struct ImmovableEvaluator : ExplicitMoveEvaluator
{
    ImmovableEvaluator() = default;
    ImmovableEvaluator(const ImmovableEvaluator&) = delete;
    ImmovableEvaluator(ImmovableEvaluator&&) = delete;
};

struct CopyOnlyEvaluator : ExplicitMoveEvaluator
{
    CopyOnlyEvaluator() = default;
    CopyOnlyEvaluator(const CopyOnlyEvaluator&) : ExplicitMoveEvaluator()
    {
    }
    CopyOnlyEvaluator(CopyOnlyEvaluator&&) = delete;
};

static void
construct_operator_builders()
{
    CopyOnlyEvaluator copy_only;
    ExplicitMoveEvaluator explicit_move;
    IBKernelEvaluatorTensorProduct<IBKernelEvaluators::IB4> lvalue{ IBKernelEvaluators::IB4{} };
    const IBOperatorBuilder from_copy_only(copy_only);
    const IBOperatorBuilder from_explicit_move(std::move(explicit_move));
    const IBOperatorBuilder from_prvalue{ ExplicitMoveEvaluator{} };
    const IBOperatorBuilder from_lvalue(lvalue);
}

static_assert(IBKernelEvaluatorCartesian<ExplicitMoveEvaluator>);
static_assert(std::constructible_from<IBOperatorBuilder, ExplicitMoveEvaluator>);
static_assert(!std::constructible_from<IBOperatorBuilder, ExplicitMoveEvaluator&>);
static_assert(IBKernelEvaluatorCartesian<CopyOnlyEvaluator>);
static_assert(std::constructible_from<IBOperatorBuilder, CopyOnlyEvaluator&>);
static_assert(IBKernelEvaluatorCartesian<ImmovableEvaluator>);
static_assert(!std::constructible_from<IBOperatorBuilder, ImmovableEvaluator>);
static_assert(std::constructible_from<IBOperatorBuilder, IBKernelEvaluatorTensorProduct<IBKernelEvaluators::IB4>&>);

double
exact_fcn(const VectorNd& x)
{
    double ret = 1.0;
    for (int d = 0; d < NDIM; ++d) ret += x[d] * static_cast<double>(d + 1);
    return ret;
}

namespace
{
// std::max keeps its first argument when the second is NaN, so a NaN
// difference must make the accumulated error infinite instead.
void
accumulate_error(double& maximum, const double error)
{
    maximum = std::isfinite(error) ? std::max(maximum, error) : std::numeric_limits<double>::infinity();
}

// A point in the neighboring rank's first cell is assigned to the nearest local patch. For a stencil of width w, its
// stencil then extends w / 2 + 1 cells past that patch, for even and for odd w, which is the ghost width that
// IBOperatorBuilder::getMinimumGhostWidth() reports. One ghost layer less must fail.
template <class Evaluator>
int
check_worst_case_ghost_width(Pointer<PatchLevel<NDIM>> level, const bool remove_layer)
{
    if (IBTK_MPI::getNodes() != 2)
    {
        TBOX_ERROR("The worst-case ghost width fixture requires two ranks.\n");
    }
    const IBKernelEvaluatorTensorProduct kernel{ Evaluator{} };
    const int ghost_width = IBOperatorBuilder(kernel).getMinimumGhostWidth() - (remove_layer ? 1 : 0);
    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, int>> indices = new SideVariable<NDIM, int>("worst_case_indices");
    const int dof = variables->registerVariableAndContext(indices, variables->getContext("worst_case"), ghost_width);
    level->allocatePatchData(dof);
    std::vector<int> counts;
    PETScVecUtilities::constructPatchLevelDOFIndices(counts, dof, level);
    PatchLevel<NDIM>::Iterator local_patch(level);
    TBOX_ASSERT(local_patch);
    Pointer<CartesianPatchGeometry<NDIM>> patch_geometry = level->getPatch(local_patch())->getPatchGeometry();
    std::array<double, NDIM> position;
    for (int d = 0; d < NDIM; ++d)
    {
        position[d] = 0.5 * (patch_geometry->getXLower()[d] + patch_geometry->getXUpper()[d]);
    }
    // The domain is split across x = 0.5. Each point lies in the other rank's cell next to the split, on the side of
    // that cell's center that makes the stencil extend farthest.
    const bool left = patch_geometry->getXLower()[0] < 0.5;
    position[0] = 0.5 + (left ? 0.9 : -0.9) * patch_geometry->getDx()[0];
    Vec X = nullptr;
    PetscErrorCode ierr = VecCreateMPI(PETSC_COMM_WORLD, NDIM, PETSC_DECIDE, &X);
    IBTK_CHKERRQ(ierr);
    PetscScalar* coordinates;
    ierr = VecGetArray(X, &coordinates);
    IBTK_CHKERRQ(ierr);
    for (int d = 0; d < NDIM; ++d)
    {
        coordinates[d] = position[d];
    }
    ierr = VecRestoreArray(X, &coordinates);
    IBTK_CHKERRQ(ierr);
    if (remove_layer)
    {
        Pointer<SAMRAI::tbox::Logger::Appender> abort_appender = new TestAppender();
        SAMRAI::tbox::Logger::getInstance()->setAbortAppender(abort_appender);
    }
    Mat matrix = nullptr;
    PETScMatUtilities::constructPatchLevelSCInterpOp(matrix, kernel, X, counts, dof, level);
    Vec field = nullptr, result = nullptr;
    ierr = MatCreateVecs(matrix, &field, &result);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(field, 1.0);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(matrix, field, result);
    IBTK_CHKERRQ(ierr);
    ierr = VecShift(result, -1.0);
    IBTK_CHKERRQ(ierr);
    PetscReal norm;
    ierr = VecNorm(result, NORM_INFINITY, &norm);
    IBTK_CHKERRQ(ierr);
    if (!std::isfinite(norm) || norm > 1.0e-12)
    {
        TBOX_ERROR("The worst-case point does not reproduce a constant: error = " << norm << ".\n");
    }
    plog << "worst-case ghost width = " << ghost_width << '\n' << "worst-case constant error = 0\n";
    ierr = VecDestroy(&result);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&field);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&matrix);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&X);
    IBTK_CHKERRQ(ierr);
    level->deallocatePatchData(dof);
    variables->removePatchDataIndex(dof);
    return 0;
}

int
check_matrix_assembly(Pointer<PatchLevel<NDIM>> level, Pointer<CartesianGridGeometry<NDIM>> geometry, bool periodic)
{
    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, int>> indices = new SideVariable<NDIM, int>("matrix_indices");
    const int dof = variables->registerVariableAndContext(indices, variables->getContext("matrix"), periodic ? 4 : 2);
    level->allocatePatchData(dof);
    std::vector<int> counts;
    PETScVecUtilities::constructPatchLevelDOFIndices(counts, dof, level);
    PatchLevel<NDIM>::Iterator local_patch(level);
    TBOX_ASSERT(local_patch);
    Pointer<Patch<NDIM>> patch = level->getPatch(local_patch());
    Pointer<CartesianPatchGeometry<NDIM>> patch_geometry = patch->getPatchGeometry();
    Pointer<SideData<NDIM, int>> dofs = patch->getPatchData(dof);
    std::array<double, NDIM> position;
    for (int d = 0; d < NDIM; ++d)
    {
        position[d] =
            0.5 * (patch_geometry->getXLower()[d] + patch_geometry->getXUpper()[d]) + 0.13 * patch_geometry->getDx()[d];
    }
    // The distributed fixture splits the domain across x = 0.5.
    if (IBTK_MPI::getNodes() > 1)
    {
        position[0] = 0.5 + (patch_geometry->getXLower()[0] < 0.5 ? -0.2 : 0.2) * patch_geometry->getDx()[0];
    }
    Vec X = nullptr;
    PetscErrorCode ierr = VecCreateMPI(PETSC_COMM_WORLD, NDIM, PETSC_DECIDE, &X);
    IBTK_CHKERRQ(ierr);
    PetscScalar* coordinates;
    ierr = VecGetArray(X, &coordinates);
    IBTK_CHKERRQ(ierr);
    for (int d = 0; d < NDIM; ++d)
    {
        coordinates[d] = position[d];
    }
    ierr = VecRestoreArray(X, &coordinates);
    IBTK_CHKERRQ(ierr);
    Mat matrix = nullptr;
    if (periodic)
    {
        // A width-eight stencil wraps around the four-cell periodic domain.
        // Contributions mapping to the same column must accumulate.
        const IBKernelEvaluatorTensorProduct kernel{ IBKernelEvaluators::BSpline<8>{} };
        double error = 0.0;
        for (int cycle = 0; cycle < 2; ++cycle)
        {
            PETScMatUtilities::constructPatchLevelSCInterpOp(matrix, kernel, X, counts, dof, level);
            Vec field = nullptr, result = nullptr;
            ierr = MatCreateVecs(matrix, &field, &result);
            IBTK_CHKERRQ(ierr);
            ierr = VecSet(field, 1.0);
            IBTK_CHKERRQ(ierr);
            ierr = MatMult(matrix, field, result);
            IBTK_CHKERRQ(ierr);
            ierr = VecShift(result, -1.0);
            IBTK_CHKERRQ(ierr);
            PetscReal norm;
            ierr = VecNorm(result, NORM_INFINITY, &norm);
            IBTK_CHKERRQ(ierr);
            TBOX_ASSERT(std::isfinite(norm));
            accumulate_error(error, static_cast<double>(norm));
            ierr = VecDestroy(&result);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&field);
            IBTK_CHKERRQ(ierr);
        }
        TBOX_ASSERT(std::isfinite(error) && error <= 1.0e-12);
        plog << "periodic constant error = " << error << '\n';
        ierr = MatDestroy(&matrix);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&X);
        IBTK_CHKERRQ(ierr);
        level->deallocatePatchData(dof);
        variables->removePatchDataIndex(dof);
        return 0;
    }
    ierr = VecLockReadPush(X);
    IBTK_CHKERRQ(ierr);
    PETScMatUtilities::constructPatchLevelSCInterpOp(
        matrix,
        IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::BSpline<3>{}, LinearIBKernel{} },
        X,
        counts,
        dof,
        level);
    ierr = VecLockReadPop(X);
    IBTK_CHKERRQ(ierr);
    Vec field = nullptr, result = nullptr;
    ierr = MatCreateVecs(matrix, &field, &result);
    IBTK_CHKERRQ(ierr);
    PetscInt column_begin, column_end, row_begin, row_end;
    ierr = VecGetOwnershipRange(field, &column_begin, &column_end);
    IBTK_CHKERRQ(ierr);
    PetscScalar* field_values;
    ierr = VecGetArray(field, &field_values);
    IBTK_CHKERRQ(ierr);
    for (PetscInt column = column_begin; column < column_end; ++column)
    {
        field_values[column - column_begin] = 1.0 + 0.001 * column;
    }
    ierr = VecRestoreArray(field, &field_values);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(matrix, field, result);
    IBTK_CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(matrix, &row_begin, &row_end);
    IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(row_end - row_begin == NDIM);
    const PetscScalar* actual;
    ierr = VecGetArrayRead(result, &actual);
    IBTK_CHKERRQ(ierr);
    int mismatches = 0, off_process = 0;
    double weight_error = 0.0, action_error = 0.0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        std::map<PetscInt, double> expected;
        // Scan grid points and use the explicit quadratic/linear spline formulas,
        // independently of the assembly stencil and evaluator implementation.
        for (SideIterator<NDIM> side(dofs->getGhostBox(), axis); side; side++)
        {
            double weight = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                const double grid_x =
                    geometry->getXLower()[d] + (side()(d) + (d == axis ? 0.0 : 0.5)) * patch_geometry->getDx()[d];
                const double r = std::abs((position[d] - grid_x) / patch_geometry->getDx()[d]);
                weight *= d == axis ? (r < 0.5 ? 0.75 - r * r : (r < 1.5 ? 0.5 * (1.5 - r) * (1.5 - r) : 0.0)) :
                                      std::max(0.0, 1.0 - r);
            }
            if (weight > 0.0)
            {
                TBOX_ASSERT((*dofs)(side()) >= 0);
                expected[(*dofs)(side())] = weight;
            }
        }
        PetscInt n;
        const PetscInt* columns;
        const PetscScalar* values;
        ierr = MatGetRow(matrix, row_begin + axis, &n, &columns, &values);
        IBTK_CHKERRQ(ierr);
        mismatches += n != static_cast<PetscInt>(expected.size());
        for (PetscInt k = 0; k < n; ++k)
        {
            const auto found = expected.find(columns[k]);
            if (found == expected.end())
            {
                ++mismatches;
            }
            else
            {
                accumulate_error(weight_error, std::abs(PetscRealPart(values[k]) - found->second));
            }
            off_process += columns[k] < column_begin || columns[k] >= column_end;
        }
        double reference = 0.0;
        for (const std::pair<const PetscInt, double>& entry : expected)
        {
            reference += entry.second * (1.0 + 0.001 * entry.first);
        }
        accumulate_error(action_error, std::abs(PetscRealPart(actual[axis]) - reference));
        ierr = MatRestoreRow(matrix, row_begin + axis, &n, &columns, &values);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(result, &actual);
    IBTK_CHKERRQ(ierr);
    mismatches = IBTK_MPI::sumReduction(mismatches);
    off_process = IBTK_MPI::sumReduction(off_process);
    weight_error = IBTK_MPI::maxReduction(weight_error);
    action_error = IBTK_MPI::maxReduction(action_error);
    plog << std::setprecision(12) << "column mismatches = " << mismatches << '\n'
         << "weight error = " << weight_error << '\n'
         << "action error = " << action_error << '\n';
    const bool valid = mismatches == 0 && weight_error <= 1.0e-12 && action_error <= 1.0e-12 &&
                       (IBTK_MPI::getNodes() == 1 || off_process > 0);
    ierr = VecDestroy(&result);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&field);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&matrix);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&X);
    IBTK_CHKERRQ(ierr);
    level->deallocatePatchData(dof);
    variables->removePatchDataIndex(dof);
    if (!valid)
    {
        TBOX_ERROR("Matrix assembly column, weight, action, or off-process coverage failed.\n");
    }
    return 0;
}
// Compare the matrix of the builder with the matrix of the evaluator that the kernel name should select. This checks
// the dispatch from a name to its evaluator; the kernel weights are checked in ib_kernel. Build the matrix with the DOF
// index ghost width that the builder requires, so that a stencil the builder under-reports fails inside the
// construction.
template <class Evaluator>
void
compare_operator_builder(const IBOperatorBuilder& builder,
                         const Evaluator& evaluator,
                         Vec X,
                         Pointer<PatchLevel<NDIM>> level,
                         int& compared,
                         int& matrix_mismatches,
                         int& row_length_mismatches)
{
    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, int>> indices =
        new SideVariable<NDIM, int>("builder_indices_" + std::to_string(compared));
    const int dof = variables->registerVariableAndContext(
        indices, variables->getContext("builder"), builder.getMinimumGhostWidth());
    level->allocatePatchData(dof);
    std::vector<int> counts;
    PETScVecUtilities::constructPatchLevelDOFIndices(counts, dof, level);
    const IBOperatorBuilder copy = builder;
    Mat from_builder = nullptr, direct = nullptr;
    copy.constructInterpolationMatrixSide(from_builder, X, counts, dof, level);
    PETScMatUtilities::constructPatchLevelSCInterpOp(direct, evaluator, X, counts, dof, level);
    PetscBool equal;
    PetscErrorCode ierr = MatEqual(from_builder, direct, &equal);
    IBTK_CHKERRQ(ierr);
    PetscInt n;
    const PetscInt* columns;
    const PetscScalar* values;
    ierr = MatGetRow(from_builder, 0, &n, &columns, &values);
    IBTK_CHKERRQ(ierr);
    row_length_mismatches += n != static_cast<PetscInt>(detail::ib_kernel_stencil_size<Evaluator, 0>());
    ierr = MatRestoreRow(from_builder, 0, &n, &columns, &values);
    IBTK_CHKERRQ(ierr);
    ++compared;
    matrix_mismatches += !equal;
    ierr = MatDestroy(&from_builder);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&direct);
    IBTK_CHKERRQ(ierr);
    level->deallocatePatchData(dof);
    variables->removePatchDataIndex(dof);
}

template <class Evaluator>
void
compare_operator_builder(const std::string& name,
                         const Evaluator& evaluator,
                         Vec X,
                         Pointer<PatchLevel<NDIM>> level,
                         int& compared,
                         int& matrix_mismatches,
                         int& row_length_mismatches)
{
    TBOX_ASSERT(IBOperatorBuilder::is_built_in(IBKernelTensorProduct(name)));
    compare_operator_builder(IBOperatorBuilder{ IBKernelTensorProduct(name) },
                             evaluator,
                             X,
                             level,
                             compared,
                             matrix_mismatches,
                             row_length_mismatches);
}

template <class... Arguments, std::size_t... Index>
void
compare_bspline_kernels(std::index_sequence<Index...>, Arguments&... arguments)
{
    (compare_operator_builder("BSPLINE_" + std::to_string(Index + 1),
                              IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::BSpline<Index + 1>{} },
                              arguments...),
     ...);
}

template <class... Arguments, std::size_t... Index>
void
compare_composite_bspline_kernels(std::index_sequence<Index...>, Arguments&... arguments)
{
    (compare_operator_builder("COMPOSITE_BSPLINE_" + std::to_string(Index + 1) + "_" + std::to_string(Index + 2),
                              IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::BSpline<Index + 1>{},
                                                              IBKernelEvaluators::BSpline<Index + 2>{} },
                              arguments...),
     ...);
    (compare_operator_builder("COMPOSITE_BSPLINE_" + std::to_string(Index + 2) + "_" + std::to_string(Index + 1),
                              IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::BSpline<Index + 2>{},
                                                              IBKernelEvaluators::BSpline<Index + 1>{} },
                              arguments...),
     ...);
}

int
check_operator_builder(Pointer<PatchLevel<NDIM>> level)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    TBOX_ASSERT(level->getNumberOfPatches() > 1);
    Pointer<CartesianGridGeometry<NDIM>> geometry = level->getGridGeometry();
    Vec X = nullptr;
    PetscErrorCode ierr = VecCreateMPI(PETSC_COMM_WORLD, NDIM, PETSC_DECIDE, &X);
    IBTK_CHKERRQ(ierr);
    PetscScalar* coordinates;
    ierr = VecGetArray(X, &coordinates);
    IBTK_CHKERRQ(ierr);
    for (int d = 0; d < NDIM; ++d)
    {
        // Just below the corner shared by the lowest patches, so that every
        // stencil extends into the neighboring patches.
        coordinates[d] = 0.5 * (geometry->getXLower()[d] + geometry->getXUpper()[d]) - 0.13 * geometry->getDx()[d];
    }
    ierr = VecRestoreArray(X, &coordinates);
    IBTK_CHKERRQ(ierr);

    int compared = 0, matrix_mismatches = 0, row_length_mismatches = 0;
    compare_bspline_kernels(std::make_index_sequence<MAX_BUILT_IN_BSPLINE_ORDER>{},
                            X,
                            level,
                            compared,
                            matrix_mismatches,
                            row_length_mismatches);
    compare_composite_bspline_kernels(std::make_index_sequence<MAX_BUILT_IN_BSPLINE_ORDER - 1>{},
                                      X,
                                      level,
                                      compared,
                                      matrix_mismatches,
                                      row_length_mismatches);
    compare_operator_builder("IB_3",
                             IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::IB3{} },
                             X,
                             level,
                             compared,
                             matrix_mismatches,
                             row_length_mismatches);
    compare_operator_builder("IB_4",
                             IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::IB4{} },
                             X,
                             level,
                             compared,
                             matrix_mismatches,
                             row_length_mismatches);
    compare_operator_builder("IB_5",
                             IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::IB5{} },
                             X,
                             level,
                             compared,
                             matrix_mismatches,
                             row_length_mismatches);
    compare_operator_builder("IB_6",
                             IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::IB6{} },
                             X,
                             level,
                             compared,
                             matrix_mismatches,
                             row_length_mismatches);
    compare_operator_builder(
        "DISCONTINUOUS_LINEAR",
        IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::BSpline<2>{}, IBKernelEvaluators::BSpline<1>{} },
        X,
        level,
        compared,
        matrix_mismatches,
        row_length_mismatches);

    // A kernel the library does not define, supplied as an evaluator.
    int application_mismatches = 0;
    {
        const IBKernelEvaluatorTensorProduct evaluator{ IBKernelEvaluators::BSpline<3>{}, LinearIBKernel{} };
        int matrix = 0, row_length = 0;
        compare_operator_builder(IBOperatorBuilder(evaluator), evaluator, X, level, compared, matrix, row_length);
        application_mismatches = matrix + row_length;
    }

    int query_errors = 0;
    for (const char* name : { "BSPLINE_9",
                              "IB_4_W8",
                              "PIECEWISE_CUBIC",
                              "COMPOSITE_BSPLINE_1_3",
                              "COMPOSITE_BSPLINE_8_9",
                              "APPLICATION_KERNEL" })
    {
        query_errors += IBOperatorBuilder::is_built_in(IBKernelTensorProduct(name));
    }
    // The builder asks for the ghost width that LEInteractor asks for when the stencil size is even, and for one layer
    // less when it is odd, because LEInteractor rounds the three- and five-point kernels up to the next even size.
    int ghost_width_mismatches = 0;
    std::string ghost_width_report;
    for (const char* name : { "BSPLINE_3",
                              "BSPLINE_4",
                              "BSPLINE_6",
                              "IB_3",
                              "IB_4",
                              "IB_5",
                              "IB_6",
                              "COMPOSITE_BSPLINE_2_3",
                              "DISCONTINUOUS_LINEAR" })
    {
        const IBKernelTensorProduct kernel(name);
        const int builder_width = IBOperatorBuilder(kernel).getMinimumGhostWidth();
        const int interactor_width = LEInteractor::getMinimumGhostWidth(kernel);
        const bool odd = std::string(name) == "BSPLINE_3" || std::string(name) == "IB_3" ||
                         std::string(name) == "IB_5" || std::string(name) == "COMPOSITE_BSPLINE_2_3";
        if (builder_width != interactor_width - (odd ? 1 : 0))
        {
            ++ghost_width_mismatches;
            ghost_width_report += std::string(" ") + name + ": builder " + std::to_string(builder_width) +
                                  ", LEInteractor " + std::to_string(interactor_width) + ";";
        }
    }
    ierr = VecDestroy(&X);
    IBTK_CHKERRQ(ierr);
    plog << "kernels compared = " << compared << '\n';
    plog << "matrix mismatches = " << matrix_mismatches << '\n';
    plog << "row length mismatches = " << row_length_mismatches << '\n';
    plog << "application kernel mismatches = " << application_mismatches << '\n';
    plog << "unsupported kernels reported built in = " << query_errors << '\n';
    if (ghost_width_mismatches != 0)
    {
        TBOX_ERROR("The operator builder and LEInteractor report different ghost widths for "
                   << ghost_width_mismatches << " kernels:" << ghost_width_report << "\n");
    }
    if (matrix_mismatches + row_length_mismatches + application_mismatches + query_errors != 0)
    {
        TBOX_ERROR("Operator builder comparison failed; see the printed mismatch counts.\n");
    }
    return 0;
}
} // namespace

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    construct_operator_builders();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create cell-centered data and extrapolate that data at physical
        // boundaries to obtain ghost cell values.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = var_db->getContext("CONTEXT");
        Pointer<CellVariable<NDIM, double>> cc_var = new CellVariable<NDIM, double>("cc");
        Pointer<SideVariable<NDIM, double>> sc_var = new SideVariable<NDIM, double>("sc");
        Pointer<NodeVariable<NDIM, double>> nc_var = new NodeVariable<NDIM, double>("nc", NDIM);
        const int gcw = 4;
        const int cc_idx = var_db->registerVariableAndContext(cc_var, context, gcw);
        const int sc_idx = var_db->registerVariableAndContext(sc_var, context, gcw);
        const int nc_idx = var_db->registerVariableAndContext(nc_var, context, gcw);

        // Initialize the AMR patch hierarchy.
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        if (input_db->getBoolWithDefault("operator_builder", false))
        {
            if (input_db->keyExists("unsupported_kernel"))
            {
                SAMRAI::tbox::Pointer<SAMRAI::tbox::Logger::Appender> abort_appender = new TestAppender();
                SAMRAI::tbox::Logger::getInstance()->setAbortAppender(abort_appender);
                const IBOperatorBuilder builder{ IBKernelTensorProduct(input_db->getString("unsupported_kernel")) };
                return 0;
            }
            return check_operator_builder(patch_hierarchy->getPatchLevel(0));
        }

        if (input_db->getBoolWithDefault("worst_case_ghost_width", false))
        {
            const bool remove_layer = input_db->getBoolWithDefault("remove_ghost_layer", false);
            const int width = input_db->getIntegerWithDefault("worst_case_kernel_width", 4);
            if (width == 3)
            {
                return check_worst_case_ghost_width<IBKernelEvaluators::BSpline<3>>(patch_hierarchy->getPatchLevel(0),
                                                                                    remove_layer);
            }
            return check_worst_case_ghost_width<IBKernelEvaluators::BSpline<4>>(patch_hierarchy->getPatchLevel(0),
                                                                                remove_layer);
        }

        if (input_db->getBoolWithDefault("matrix_assembly", false))
        {
            return check_matrix_assembly(patch_hierarchy->getPatchLevel(0),
                                         grid_geometry,
                                         input_db->getBoolWithDefault("periodic_matrix", false));
        }

        // Allocate and fill in patch data
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(cc_idx);
            level->allocatePatchData(nc_idx);
            level->allocatePatchData(sc_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> cc_data = patch->getPatchData(cc_idx);
                Pointer<SideData<NDIM, double>> sc_data = patch->getPatchData(sc_idx);
                Pointer<NodeData<NDIM, double>> nc_data = patch->getPatchData(nc_idx);
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    VectorNd x;
                    for (int d = 0; d < NDIM; ++d)
                        x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                    (*cc_data)(idx) = exact_fcn(x);
                }

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                    {
                        const SideIndex<NDIM>& idx = si();
                        VectorNd x;
                        for (int d = 0; d < NDIM; ++d)
                            x[d] =
                                xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                        (*sc_data)(idx) = exact_fcn(x);
                    }
                }

                for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
                {
                    const NodeIndex<NDIM>& idx = ni();
                    VectorNd x;
                    for (int d = 0; d < NDIM; ++d) x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)));
                    for (int d = 0; d < NDIM; ++d) (*nc_data)(idx, d) = exact_fcn(x);
                }
            }
        }

        // Now fill ghost cells
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comps{ ITC(cc_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE"),
                                           ITC(sc_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE"),
                                           ITC(nc_idx, "LINEAR_REFINE", false, "NONE") };
        HierarchyGhostCellInterpolation ghost_cell_fill;
        ghost_cell_fill.initializeOperatorState(ghost_cell_comps, patch_hierarchy, coarsest_ln, finest_ln);
        ghost_cell_fill.fillData(0.0);

        // Now interpolate to the specified point.
        std::vector<VectorNd> x_pt(2);
        for (int d = 0; d < NDIM; ++d) x_pt[0][d] = 0.7;
        for (int d = 0; d < NDIM; ++d) x_pt[1][d] = 0.2;

        // Cell centered
        pout << "Interpolating cell centered values\n";
        std::vector<double> interped_val = interpolate(x_pt, cc_idx, cc_var, 1, patch_hierarchy, "IB_4");
        const IBKernelTensorProduct resolved({ IBKernel::IB_4 });
        TBOX_ASSERT(interped_val == interpolate(x_pt, cc_idx, cc_var, 1, patch_hierarchy, resolved));
        TBOX_ASSERT(interpolate(x_pt[0], cc_idx, cc_var, 1, patch_hierarchy, resolved)[0] == interped_val[0]);
        for (int i = 0; i < 2; ++i)
        {
            bool correct = std::abs(interped_val[i] - exact_fcn(x_pt[i])) < 1.0e-12;
            correct = IBTK_MPI::maxReduction(correct ? 0 : 1) == 0;
            if (!correct)
            {
                plog << "Interpolant number " << i << " was not exact!\n";
                plog << "Expected " << exact_fcn(x_pt[i]) << " and got " << interped_val[i] << "\n";
                plog << "Error: " << interped_val[i] - exact_fcn(x_pt[i]) << "\n";
            }
        }

        // Side centered
        pout << "Interpolating side centered values\n";
        interped_val = interpolate(x_pt, sc_idx, sc_var, 1, patch_hierarchy, "IB_4");
        TBOX_ASSERT(interped_val == interpolate(x_pt, sc_idx, sc_var, 1, patch_hierarchy, resolved));
        for (int i = 0; i < 2; ++i)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                bool correct = std::abs(interped_val[i * NDIM + d] - exact_fcn(x_pt[i])) < 1.0e-12;
                correct = IBTK_MPI::maxReduction(correct ? 0 : 1) == 0;
                if (!correct)
                {
                    plog << "Interpolant number " << i << " and depth " << d << " was not exact!\n";
                    plog << "Expected " << exact_fcn(x_pt[i]) << " and got " << interped_val[i] << "\n";
                    plog << "Error: " << interped_val[i] - exact_fcn(x_pt[i]) << "\n";
                }
            }
        }

        // Node centered
        pout << "Interpolating node centered values\n";
        interped_val = interpolate(x_pt, nc_idx, nc_var, NDIM, patch_hierarchy, "IB_4");
        TBOX_ASSERT(interped_val == interpolate(x_pt, nc_idx, nc_var, NDIM, patch_hierarchy, resolved));
        for (int i = 0; i < 2; ++i)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                bool correct = std::abs(interped_val[i * NDIM + d] - exact_fcn(x_pt[i])) < 1.0e-12;
                correct = IBTK_MPI::maxReduction(correct ? 0 : 1) == 0;
                if (!correct)
                {
                    plog << "Interpolant number " << i << " and depth " << d << " was not exact!\n";
                    plog << "Expected " << exact_fcn(x_pt[i]) << " and got " << interped_val[i] << "\n";
                    plog << "Error: " << interped_val[i] - exact_fcn(x_pt[i]) << "\n";
                }
            }
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
