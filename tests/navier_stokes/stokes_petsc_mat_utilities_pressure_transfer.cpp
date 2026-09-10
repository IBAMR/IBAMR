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
 * This test checks linear cell-centered pressure prolongation and restriction.
 * For each transfer it compares the PETSc matrix with the corresponding SAMRAI
 * schedule on an analytic profile. It also checks the matrix and the schedule
 * against the exact result for a constant profile and for an affine profile away
 * from the domain boundary, which verifies constant preservation and interior
 * affine reproduction independently of that comparison. Restriction is also
 * checked for the cell-volume adjoint of prolongation. The domain may be
 * periodic or bounded and the levels may be partitioned into patches and
 * distributed over MPI ranks.
 */

#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartCellDoubleLinearCoarsen.h>
#include <ibtk/CartCellDoubleLinearRefine.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/ibtk_utilities.h>

#include <petscmat.h>

#include <boost/math/constants/constants.hpp>

#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellGeometry.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <HierarchyCellDataOpsReal.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefineSchedule.h>
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
#include <limits>
#include <numeric>
#include <string>
#include <utility>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
enum ProfileKind
{
    TEST_PROFILE,
    CONSTANT_PROFILE,
    AFFINE_PROFILE
};

// Prolongation and restriction use different test profiles: the restriction
// coefficients are read with the prefix "restriction_".
double
profile_value(const ProfileKind kind, const VectorNd& X, Pointer<Database> test_db, const std::string& prefix)
{
    if (kind == CONSTANT_PROFILE)
    {
        return 1.3;
    }
    if (kind == AFFINE_PROFILE)
    {
        double value = 0.3;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            value += (0.7 - 1.1 * axis) * X[axis];
        }
        return value;
    }
    const bool nonlinear = test_db->getString("profile_type") == "nonlinear";
    const double pi = boost::math::constants::pi<double>();
    double value = test_db->getDoubleWithDefault(prefix + "profile_constant", 0.7);
    double phase = 0.0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const double coefficient =
            test_db->getDoubleWithDefault(prefix + "profile_coefficient_" + std::to_string(axis), 1.0);
        value += coefficient * std::sin(2.0 * pi * (axis + 1.0) * X[axis]);
        phase += (axis + 1.0) * X[axis];
    }
    if (nonlinear)
    {
        value += 0.35 * std::sin(2.0 * pi * phase) + 0.2 * std::cos(4.0 * pi * phase);
    }
    return value;
}

// Each depth is a distinct affine function of the profile.
double
depth_value(const double value, const int depth)
{
    return (depth + 1.0) * value + 0.1 * depth;
}

void
set_pressure_profile(Pointer<CellData<NDIM, double>> data,
                     Pointer<Patch<NDIM>> patch,
                     const ProfileKind kind,
                     Pointer<Database> test_db,
                     const std::string& prefix)
{
    for (Box<NDIM>::Iterator b(data->getGhostBox()); b; b++)
    {
        const double value = profile_value(kind, IndexUtilities::getCellCenter(*patch, b()), test_db, prefix);
        for (int depth = 0; depth < data->getDepth(); ++depth)
        {
            (*data)(b(), depth) = depth_value(value, depth);
        }
    }
}

// Cells farther than a coarse cell width from every domain boundary, where
// linear interpolation reproduces affine data exactly, whether the boundary is
// physical or periodic.
bool
is_interior(const VectorNd& X, const double* x_lower, const double* x_upper, const double* coarse_dx)
{
    for (int axis = 0; axis < NDIM; ++axis)
    {
        if (X[axis] < x_lower[axis] + coarse_dx[axis] || X[axis] > x_upper[axis] - coarse_dx[axis])
        {
            return false;
        }
    }
    return true;
}

// std::max keeps its first argument when the second is NaN, so a nonfinite
// difference must make the accumulated error infinite instead.
double
accumulate_error(const double error, const double difference)
{
    return std::isfinite(difference) ? std::max(error, difference) : std::numeric_limits<double>::infinity();
}

// The pressure rows of the prolongation matrix, whose fine-level DOFs are numbered by p_dof_index_idx.
// Return the number of stored entries in those rows that lie in columns owned by another rank, and the
// smallest number of local pressure rows on any rank.
std::pair<PetscInt, PetscInt>
count_off_process_pressure_entries(Mat mat, const int p_dof_index_idx, Pointer<PatchLevel<NDIM>> fine_level)
{
    PetscInt column_first = 0, column_last = 0;
    int ierr = MatGetOwnershipRangeColumn(mat, &column_first, &column_last);
    IBTK_CHKERRQ(ierr);
    PetscInt off_process_entries = 0, local_rows = 0;
    for (PatchLevel<NDIM>::Iterator p(fine_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = fine_level->getPatch(p());
        Pointer<CellData<NDIM, int>> dof_data = patch->getPatchData(p_dof_index_idx);
        for (Box<NDIM>::Iterator b(patch->getBox()); b; b++)
        {
            for (int depth = 0; depth < dof_data->getDepth(); ++depth)
            {
                PetscInt num_columns = 0;
                const PetscInt* columns = nullptr;
                ierr = MatGetRow(mat, (*dof_data)(b(), depth), &num_columns, &columns, nullptr);
                IBTK_CHKERRQ(ierr);
                for (PetscInt k = 0; k < num_columns; ++k)
                {
                    off_process_entries += columns[k] < column_first || columns[k] >= column_last;
                }
                ierr = MatRestoreRow(mat, (*dof_data)(b(), depth), &num_columns, &columns, nullptr);
                IBTK_CHKERRQ(ierr);
                ++local_rows;
            }
        }
    }
    return { IBTK_MPI::sumReduction(off_process_entries), IBTK_MPI::minReduction(local_rows) };
}

struct TransferResult
{
    double input_norm = 0.0;
    double matrix_norm = 0.0;
    double schedule_norm = 0.0;
    double transfer_error = 0.0;
    double exact_error = 0.0;
};

struct TransferFixture
{
    Pointer<Database> test_db;
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy;
    Pointer<PatchLevel<NDIM>> coarse_level, fine_level;
    IntVector<NDIM> level_ratio;
    bool cell_only = false;
    int depth = 1;
    int u_idx = -1, p_idx = -1, result_idx = -1;
    int u_dof_index_idx = -1, p_dof_index_idx = -1;
    Mat prolong_mat = nullptr;
    const double* x_lower = nullptr;
    const double* x_upper = nullptr;
    const double* coarse_dx = nullptr;
};

// Transfer one profile with the matrix and with the SAMRAI schedule.
TransferResult
transfer_profile(TransferFixture& fixture, const bool restriction, const ProfileKind kind, Mat transfer)
{
    Pointer<PatchLevel<NDIM>> input_level = restriction ? fixture.fine_level : fixture.coarse_level;
    Pointer<PatchLevel<NDIM>> output_level = restriction ? fixture.coarse_level : fixture.fine_level;
    const std::string prefix = restriction ? "restriction_" : "";
    for (PatchLevel<NDIM>::Iterator p(input_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = input_level->getPatch(p());
        Pointer<CellData<NDIM, double>> data = patch->getPatchData(fixture.p_idx);
        set_pressure_profile(data, patch, kind, fixture.test_db, prefix);
    }
    for (PatchLevel<NDIM>::Iterator p(output_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = output_level->getPatch(p());
        Pointer<CellData<NDIM, double>> data = patch->getPatchData(fixture.result_idx);
        data->fillAll(0.0);
    }

    if (restriction)
    {
        const int required_width = std::max(1, fixture.level_ratio.max() / 2);
        Pointer<CoarsenAlgorithm<NDIM>> algorithm = new CoarsenAlgorithm<NDIM>();
        algorithm->registerCoarsen(fixture.result_idx,
                                   fixture.p_idx,
                                   new CartCellDoubleLinearCoarsen(IntVector<NDIM>(
                                       fixture.test_db->getIntegerWithDefault("declared_width", required_width))),
                                   IntVector<NDIM>(0));
        Pointer<CoarsenSchedule<NDIM>> schedule = algorithm->createSchedule(fixture.coarse_level, fixture.fine_level);
        schedule->coarsenData();
    }
    else
    {
        Pointer<RefineAlgorithm<NDIM>> algorithm = new RefineAlgorithm<NDIM>();
        algorithm->registerRefine(
            fixture.result_idx, fixture.p_idx, fixture.result_idx, new CartCellDoubleLinearRefine(), nullptr);
        Pointer<RefineSchedule<NDIM>> schedule = algorithm->createSchedule(
            fixture.fine_level, Pointer<PatchLevel<NDIM>>(), 0, fixture.patch_hierarchy, nullptr);
        schedule->fillData(0.0);
    }

    // Constant profiles are reproduced everywhere for prolongation, and away from the boundary for
    // restriction; affine profiles are reproduced away from the boundary.
    const auto exact_error = [&]()
    {
        double error = 0.0;
        for (PatchLevel<NDIM>::Iterator p(output_level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = output_level->getPatch(p());
            Pointer<CellData<NDIM, double>> data = patch->getPatchData(fixture.result_idx);
            for (Box<NDIM>::Iterator b(patch->getBox()); b; b++)
            {
                const VectorNd X = IndexUtilities::getCellCenter(*patch, b());
                if ((kind == AFFINE_PROFILE || restriction) &&
                    !is_interior(X, fixture.x_lower, fixture.x_upper, fixture.coarse_dx))
                {
                    continue;
                }
                const double expected = profile_value(kind, X, fixture.test_db, prefix);
                for (int d = 0; d < fixture.depth; ++d)
                {
                    error = accumulate_error(error, std::abs((*data)(b(), d) - depth_value(expected, d)));
                }
            }
        }
        return IBTK_MPI::maxReduction(error);
    };
    TransferResult result;
    if (kind != TEST_PROFILE)
    {
        result.exact_error = exact_error();
    }

    Vec x = nullptr, matrix_result = nullptr, schedule_result = nullptr;
    int ierr = MatCreateVecs(transfer, &x, &matrix_result);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(matrix_result, &schedule_result);
    IBTK_CHKERRQ(ierr);
    if (fixture.cell_only)
    {
        PETScVecUtilities::copyToPatchLevelVec(x, fixture.p_idx, fixture.p_dof_index_idx, input_level);
        PETScVecUtilities::copyToPatchLevelVec(
            schedule_result, fixture.result_idx, fixture.p_dof_index_idx, output_level);
    }
    else
    {
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            x, fixture.u_idx, fixture.u_dof_index_idx, fixture.p_idx, fixture.p_dof_index_idx, input_level);
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(schedule_result,
                                                                     fixture.u_idx,
                                                                     fixture.u_dof_index_idx,
                                                                     fixture.result_idx,
                                                                     fixture.p_dof_index_idx,
                                                                     output_level);
    }
    ierr = MatMult(transfer, x, matrix_result);
    IBTK_CHKERRQ(ierr);
    ierr = VecNorm(x, NORM_INFINITY, &result.input_norm);
    IBTK_CHKERRQ(ierr);
    ierr = VecNorm(matrix_result, NORM_INFINITY, &result.matrix_norm);
    IBTK_CHKERRQ(ierr);
    ierr = VecNorm(schedule_result, NORM_INFINITY, &result.schedule_norm);
    IBTK_CHKERRQ(ierr);
    if (kind == TEST_PROFILE)
    {
        ierr = VecAXPY(matrix_result, -1.0, schedule_result);
        IBTK_CHKERRQ(ierr);
        ierr = VecNorm(matrix_result, NORM_INFINITY, &result.transfer_error);
        IBTK_CHKERRQ(ierr);
    }
    else
    {
        // Check the matrix result against the same exact values as the schedule, retaining the
        // larger error.
        if (fixture.cell_only)
        {
            PETScVecUtilities::copyFromPatchLevelVec(
                matrix_result, fixture.result_idx, fixture.p_dof_index_idx, output_level, nullptr, nullptr);
        }
        else
        {
            IBAMR::StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(matrix_result,
                                                                           fixture.u_idx,
                                                                           fixture.u_dof_index_idx,
                                                                           fixture.result_idx,
                                                                           fixture.p_dof_index_idx,
                                                                           output_level,
                                                                           nullptr,
                                                                           nullptr);
        }
        result.exact_error = accumulate_error(result.exact_error, exact_error());
    }
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&matrix_result);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&schedule_result);
    IBTK_CHKERRQ(ierr);
    return result;
}

// A direct sub-box call must reproduce the schedule and preserve all other cells.
double
check_write_range(TransferFixture& fixture)
{
    Pointer<Patch<NDIM>> input_patch = fixture.fine_level->getPatch(0);
    Pointer<Patch<NDIM>> output_patch = fixture.coarse_level->getPatch(0);
    Pointer<CellData<NDIM, double>> p_output = output_patch->getPatchData(fixture.result_idx);
    Pointer<CellData<NDIM, double>> saved =
        new CellData<NDIM, double>(p_output->getBox(), fixture.depth, p_output->getGhostCellWidth());
    saved->copy(*p_output);
    const int required_width = std::max(1, fixture.level_ratio.max() / 2);
    double write_range_error = 0.0;
    for (int side = 0; side < 2; ++side)
    {
        Box<NDIM> sub_box = output_patch->getBox();
        for (int axis = 0; axis < NDIM; ++axis)
        {
            if (side == 0)
            {
                sub_box.upper()(axis) = std::min(sub_box.upper()(axis), sub_box.lower()(axis) + 1);
            }
            else
            {
                sub_box.lower()(axis) = std::max(sub_box.lower()(axis), sub_box.upper()(axis) - 1);
            }
        }
        p_output->fillAll(-12345.0);
        CartCellDoubleLinearCoarsen coarsen_op{ IntVector<NDIM>(required_width) };
        coarsen_op.coarsen(
            *output_patch, *input_patch, fixture.result_idx, fixture.p_idx, sub_box, fixture.level_ratio);
        for (Box<NDIM>::Iterator b(p_output->getGhostBox()); b; b++)
        {
            for (int d = 0; d < fixture.depth; ++d)
            {
                const double expected = sub_box.contains(b()) ? (*saved)(b(), d) : -12345.0;
                write_range_error = accumulate_error(write_range_error, std::abs((*p_output)(b(), d) - expected));
            }
        }
    }
    p_output->copy(*saved);
    return write_range_error;
}

// Restriction is the cell-volume adjoint of prolongation: (P u_c, v_f) V_f = (u_c, R v_f) V_c.
double
check_volume_adjoint(TransferFixture& fixture, Vec fine_input, Vec coarse_output)
{
    Vec coarse_probe = nullptr, prolonged_probe = nullptr;
    int ierr = MatCreateVecs(fixture.prolong_mat, &coarse_probe, &prolonged_probe);
    IBTK_CHKERRQ(ierr);
    PetscScalar* values = nullptr;
    PetscInt first = 0, last = 0;
    ierr = VecGetOwnershipRange(coarse_probe, &first, &last);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(coarse_probe, &values);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = first; i < last; ++i)
    {
        values[i - first] = std::cos(0.17 * (i + 1));
    }
    ierr = VecRestoreArray(coarse_probe, &values);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(fixture.prolong_mat, coarse_probe, prolonged_probe);
    IBTK_CHKERRQ(ierr);
    PetscScalar fine_product = 0.0, coarse_product = 0.0;
    ierr = VecDot(prolonged_probe, fine_input, &fine_product);
    IBTK_CHKERRQ(ierr);
    ierr = VecDot(coarse_probe, coarse_output, &coarse_product);
    IBTK_CHKERRQ(ierr);
    for (int axis = 0; axis < NDIM; ++axis)
    {
        fine_product *= fixture.coarse_dx[axis] / fixture.level_ratio(axis);
        coarse_product *= fixture.coarse_dx[axis];
    }
    ierr = VecDestroy(&coarse_probe);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&prolonged_probe);
    IBTK_CHKERRQ(ierr);
    return std::abs(fine_product - coarse_product);
}

bool
transfer_ok(const TransferResult& r, const bool require_signal)
{
    const bool finite = std::isfinite(r.input_norm) && std::isfinite(r.matrix_norm) && std::isfinite(r.schedule_norm) &&
                        std::isfinite(r.transfer_error) && std::isfinite(r.exact_error);
    return finite && r.transfer_error < 1.0e-12 && r.exact_error < 1.0e-12 &&
           (!require_signal || (r.input_norm > 1.0e-8 && r.matrix_norm > 1.0e-8 && r.schedule_norm > 1.0e-8));
}

void
print_transfer(const std::string& name, const TransferResult& r, const bool print_exact)
{
    plog << name << " input norm = " << r.input_norm << '\n'
         << name << " matrix norm = " << r.matrix_norm << '\n'
         << name << " schedule norm = " << r.schedule_norm << '\n'
         << name << " transfer error = " << r.transfer_error << '\n';
    if (print_exact)
    {
        plog << name << " exact error = " << r.exact_error << '\n';
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Logger::Appender> abort_appender = new TestAppender();
    Logger::getInstance()->setAbortAppender(abort_appender);
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    Pointer<Database> test_db = input_db->getDatabase("test");
    if (test_db->keyExists("invalid_width"))
    {
        CartCellDoubleLinearCoarsen coarsen_op(IntVector<NDIM>(test_db->getInteger("invalid_width")));
        return 0;
    }

    TransferFixture fixture;
    fixture.test_db = test_db;
    fixture.depth = test_db->getIntegerWithDefault("depth", 1);
    fixture.cell_only = fixture.depth > 1;
    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app_initializer);
    fixture.patch_hierarchy = std::get<0>(hierarchy_tuple);
    if (fixture.patch_hierarchy->getFinestLevelNumber() != 1)
    {
        TBOX_ERROR("Pressure transfer test requires two levels.\n");
    }
    fixture.coarse_level = fixture.patch_hierarchy->getPatchLevel(0);
    fixture.fine_level = fixture.patch_hierarchy->getPatchLevel(1);
    fixture.level_ratio = fixture.fine_level->getRatio() / fixture.coarse_level->getRatio();
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = fixture.patch_hierarchy->getGridGeometry();
    fixture.x_lower = grid_geometry->getXLower();
    fixture.x_upper = grid_geometry->getXUpper();
    fixture.coarse_dx = grid_geometry->getDx();
    const bool partitioned = test_db->getBoolWithDefault("partitioned", false);
    if (partitioned && fixture.fine_level->getNumberOfPatches() <= 1)
    {
        TBOX_ERROR("The partitioned transfer case requires multiple fine patches.\n");
    }
    const bool single_patch = fixture.coarse_level->getNumberOfPatches() == 1 &&
                              fixture.fine_level->getNumberOfPatches() == 1 && IBTK_MPI::getNodes() == 1;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("stokes_petsc_mat_utilities_pressure_transfer_ctx");
    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("pressure_transfer_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var =
        new CellVariable<NDIM, int>("pressure_transfer_p_dof", fixture.depth);
    fixture.u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, ctx, IntVector<NDIM>(1));
    fixture.p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, ctx, IntVector<NDIM>(1));
    for (Pointer<PatchLevel<NDIM>> level : { fixture.coarse_level, fixture.fine_level })
    {
        level->allocatePatchData(fixture.u_dof_index_idx);
        level->allocatePatchData(fixture.p_dof_index_idx);
    }

    std::vector<int> num_coarse_dofs_per_proc, num_fine_dofs_per_proc;
    AO coarse_level_ao = nullptr;
    int u_coarse_ao_offset = 0, p_coarse_ao_offset = 0;
    if (fixture.cell_only)
    {
        PETScVecUtilities::constructPatchLevelDOFIndices(
            num_coarse_dofs_per_proc, fixture.p_dof_index_idx, fixture.coarse_level);
        PETScVecUtilities::constructPatchLevelDOFIndices(
            num_fine_dofs_per_proc, fixture.p_dof_index_idx, fixture.fine_level);
        PETScVecUtilities::constructPatchLevelAO(coarse_level_ao,
                                                 num_coarse_dofs_per_proc,
                                                 fixture.p_dof_index_idx,
                                                 fixture.coarse_level,
                                                 p_coarse_ao_offset);
        PETScMatUtilities::constructProlongationOp(fixture.prolong_mat,
                                                   "LINEAR",
                                                   fixture.p_dof_index_idx,
                                                   num_fine_dofs_per_proc,
                                                   num_coarse_dofs_per_proc,
                                                   fixture.fine_level,
                                                   fixture.coarse_level,
                                                   coarse_level_ao,
                                                   p_coarse_ao_offset);
    }
    else
    {
        IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
            num_coarse_dofs_per_proc, fixture.u_dof_index_idx, fixture.p_dof_index_idx, fixture.coarse_level);
        IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
            num_fine_dofs_per_proc, fixture.u_dof_index_idx, fixture.p_dof_index_idx, fixture.fine_level);
        IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelAO(coarse_level_ao,
                                                                       num_coarse_dofs_per_proc,
                                                                       fixture.u_dof_index_idx,
                                                                       fixture.p_dof_index_idx,
                                                                       fixture.coarse_level,
                                                                       u_coarse_ao_offset,
                                                                       p_coarse_ao_offset);
        IBAMR::StaggeredStokesPETScMatUtilities::constructProlongationOp(fixture.prolong_mat,
                                                                         "LINEAR",
                                                                         "LINEAR",
                                                                         fixture.u_dof_index_idx,
                                                                         fixture.p_dof_index_idx,
                                                                         num_fine_dofs_per_proc,
                                                                         num_coarse_dofs_per_proc,
                                                                         fixture.fine_level,
                                                                         fixture.coarse_level,
                                                                         coarse_level_ao,
                                                                         u_coarse_ao_offset,
                                                                         p_coarse_ao_offset);
    }

    Pointer<SideVariable<NDIM, double>> velocity = new SideVariable<NDIM, double>("velocity");
    Pointer<CellVariable<NDIM, double>> pressure = new CellVariable<NDIM, double>("pressure", fixture.depth);
    Pointer<CellVariable<NDIM, double>> transferred = new CellVariable<NDIM, double>("transferred", fixture.depth);
    const int required_width = std::max(1, fixture.level_ratio.max() / 2);
    fixture.u_idx = var_db->registerVariableAndContext(velocity, ctx, IntVector<NDIM>(1));
    fixture.p_idx = var_db->registerVariableAndContext(pressure, ctx, IntVector<NDIM>(required_width));
    fixture.result_idx = var_db->registerVariableAndContext(transferred, ctx, IntVector<NDIM>(1));
    for (Pointer<PatchLevel<NDIM>> level : { fixture.coarse_level, fixture.fine_level })
    {
        level->allocatePatchData(fixture.u_idx);
        level->allocatePatchData(fixture.p_idx);
        level->allocatePatchData(fixture.result_idx);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<SideData<NDIM, double>> u = level->getPatch(p())->getPatchData(fixture.u_idx);
            u->fillAll(0.0);
        }
    }

    // Prolongation.
    bool success = true;
    const TransferResult prolongation_test = transfer_profile(fixture, false, TEST_PROFILE, fixture.prolong_mat);
    const TransferResult prolongation_constant =
        transfer_profile(fixture, false, CONSTANT_PROFILE, fixture.prolong_mat);
    const TransferResult prolongation_affine = transfer_profile(fixture, false, AFFINE_PROFILE, fixture.prolong_mat);
    print_transfer("prolongation", prolongation_test, false);
    plog << "prolongation constant error = " << prolongation_constant.exact_error << '\n'
         << "prolongation affine interior error = " << prolongation_affine.exact_error << '\n';
    success = success && transfer_ok(prolongation_test, true) && transfer_ok(prolongation_constant, false) &&
              transfer_ok(prolongation_affine, false);

    // Restriction is the cell-volume-scaled transpose of prolongation.
    Mat restriction_mat = nullptr;
    int ierr = MatTranspose(fixture.prolong_mat, MAT_INITIAL_MATRIX, &restriction_mat);
    IBTK_CHKERRQ(ierr);
    double volume_ratio = 1.0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        volume_ratio /= fixture.level_ratio(axis);
    }
    ierr = MatScale(restriction_mat, volume_ratio);
    IBTK_CHKERRQ(ierr);
    const TransferResult restriction_constant = transfer_profile(fixture, true, CONSTANT_PROFILE, restriction_mat);
    const TransferResult restriction_affine = transfer_profile(fixture, true, AFFINE_PROFILE, restriction_mat);
    const TransferResult restriction_test = transfer_profile(fixture, true, TEST_PROFILE, restriction_mat);
    double write_range_error = 0.0;
    if (single_patch)
    {
        write_range_error = check_write_range(fixture);
    }

    // Recompute the test-profile vectors for the adjoint check.
    Vec fine_input = nullptr, coarse_output = nullptr;
    ierr = MatCreateVecs(restriction_mat, &fine_input, &coarse_output);
    IBTK_CHKERRQ(ierr);
    if (fixture.cell_only)
    {
        PETScVecUtilities::copyToPatchLevelVec(fine_input, fixture.p_idx, fixture.p_dof_index_idx, fixture.fine_level);
        PETScVecUtilities::copyToPatchLevelVec(
            coarse_output, fixture.result_idx, fixture.p_dof_index_idx, fixture.coarse_level);
    }
    else
    {
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(fine_input,
                                                                     fixture.u_idx,
                                                                     fixture.u_dof_index_idx,
                                                                     fixture.p_idx,
                                                                     fixture.p_dof_index_idx,
                                                                     fixture.fine_level);
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(coarse_output,
                                                                     fixture.u_idx,
                                                                     fixture.u_dof_index_idx,
                                                                     fixture.result_idx,
                                                                     fixture.p_dof_index_idx,
                                                                     fixture.coarse_level);
    }
    const double adjoint_error = check_volume_adjoint(fixture, fine_input, coarse_output);
    ierr = VecDestroy(&fine_input);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&coarse_output);
    IBTK_CHKERRQ(ierr);

    print_transfer("restriction", restriction_test, false);
    plog << "restriction constant error = " << restriction_constant.exact_error << '\n'
         << "restriction affine interior error = " << restriction_affine.exact_error << '\n';
    if (single_patch)
    {
        plog << "restriction write range error = " << write_range_error << '\n';
    }
    plog << "restriction volume adjoint error = " << adjoint_error << '\n';
    success = success && transfer_ok(restriction_test, true) && transfer_ok(restriction_constant, false) &&
              transfer_ok(restriction_affine, false) && std::isfinite(adjoint_error) && adjoint_error < 1.0e-12 &&
              std::isfinite(write_range_error) && write_range_error < 1.0e-12;

    // With several ranks, every rank must own pressure rows and some pressure stencil entries must
    // refer to columns owned by another rank.
    if (IBTK_MPI::getNodes() > 1)
    {
        const std::pair<PetscInt, PetscInt> counts =
            count_off_process_pressure_entries(fixture.prolong_mat, fixture.p_dof_index_idx, fixture.fine_level);
        plog << "off-process pressure matrix entries = " << counts.first << '\n';
        success = success && counts.first > 0 && counts.second > 0;
    }

    ierr = MatDestroy(&restriction_mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&fixture.prolong_mat);
    IBTK_CHKERRQ(ierr);
    ierr = AODestroy(&coarse_level_ao);
    IBTK_CHKERRQ(ierr);
    return success ? 0 : 1;
}
