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
 * This test checks linear cell-centered pressure prolongation and restriction
 * parity between the PETSc transfer operators and the corresponding SAMRAI
 * transfer schedules on analytic profiles, including physical boundaries.
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
#include <numeric>
#include <string>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
void
set_pressure_profile(Pointer<CellData<NDIM, double>> data, Pointer<Patch<NDIM>> patch, Pointer<Database> test_db)
{
    const bool nonlinear = test_db->getString("profile_type") == "nonlinear";
    const double pi = boost::math::constants::pi<double>();
    for (Box<NDIM>::Iterator b(data->getGhostBox()); b; b++)
    {
        const VectorNd X = IndexUtilities::getCellCenter(*patch, b());
        double value = test_db->getDoubleWithDefault("profile_constant", 0.7);
        double phase = 0.0;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const double coefficient =
                test_db->getDoubleWithDefault("profile_coefficient_" + std::to_string(axis), 1.0);
            value += coefficient * std::sin(2.0 * pi * (axis + 1.0) * X[axis]);
            phase += (axis + 1.0) * X[axis];
        }
        if (nonlinear)
        {
            value += 0.35 * std::sin(2.0 * pi * phase) + 0.2 * std::cos(4.0 * pi * phase);
        }
        for (int depth = 0; depth < data->getDepth(); ++depth)
        {
            (*data)(b(), depth) = (depth + 1.0) * value + 0.1 * depth;
        }
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
    const int depth = test_db->getIntegerWithDefault("depth", 1);
    const bool cell_only = depth > 1;
    const bool restriction = test_db->getString("test_mode") == "restriction";
    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app_initializer);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = std::get<0>(hierarchy_tuple);
    if (patch_hierarchy->getFinestLevelNumber() != 1 || IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("Pressure transfer test requires two levels on one MPI rank.\n");
    }
    const int coarse_ln = 0;
    const int fine_ln = 1;
    Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(coarse_ln);
    Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(fine_ln);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("stokes_petsc_mat_utilities_pressure_transfer_ctx");

    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("pressure_transfer_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("pressure_transfer_p_dof", depth);
    const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, ctx, IntVector<NDIM>(1));
    const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, ctx, IntVector<NDIM>(1));

    coarse_level->allocatePatchData(u_dof_index_idx);
    coarse_level->allocatePatchData(p_dof_index_idx);
    fine_level->allocatePatchData(u_dof_index_idx);
    fine_level->allocatePatchData(p_dof_index_idx);

    std::vector<int> num_coarse_dofs_per_proc, num_fine_dofs_per_proc;
    if (cell_only)
    {
        PETScVecUtilities::constructPatchLevelDOFIndices(num_coarse_dofs_per_proc, p_dof_index_idx, coarse_level);
        PETScVecUtilities::constructPatchLevelDOFIndices(num_fine_dofs_per_proc, p_dof_index_idx, fine_level);
    }
    else
    {
        IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
            num_coarse_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, coarse_level);
        IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
            num_fine_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, fine_level);
    }

    AO coarse_level_ao = nullptr;
    int u_coarse_ao_offset = 0, p_coarse_ao_offset = 0;
    if (cell_only)
    {
        PETScVecUtilities::constructPatchLevelAO(
            coarse_level_ao, num_coarse_dofs_per_proc, p_dof_index_idx, coarse_level, p_coarse_ao_offset);
    }
    else
    {
        IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelAO(coarse_level_ao,
                                                                       num_coarse_dofs_per_proc,
                                                                       u_dof_index_idx,
                                                                       p_dof_index_idx,
                                                                       coarse_level,
                                                                       u_coarse_ao_offset,
                                                                       p_coarse_ao_offset);
    }

    Mat prolong_mat = nullptr;
    if (cell_only)
    {
        PETScMatUtilities::constructProlongationOp(prolong_mat,
                                                   "LINEAR",
                                                   p_dof_index_idx,
                                                   num_fine_dofs_per_proc,
                                                   num_coarse_dofs_per_proc,
                                                   fine_level,
                                                   coarse_level,
                                                   coarse_level_ao,
                                                   p_coarse_ao_offset);
    }
    else
    {
        IBAMR::StaggeredStokesPETScMatUtilities::constructProlongationOp(prolong_mat,
                                                                         "LINEAR",
                                                                         "LINEAR",
                                                                         u_dof_index_idx,
                                                                         p_dof_index_idx,
                                                                         num_fine_dofs_per_proc,
                                                                         num_coarse_dofs_per_proc,
                                                                         fine_level,
                                                                         coarse_level,
                                                                         coarse_level_ao,
                                                                         u_coarse_ao_offset,
                                                                         p_coarse_ao_offset);
    }

    Pointer<PatchLevel<NDIM>> input_level = restriction ? fine_level : coarse_level;
    Pointer<PatchLevel<NDIM>> output_level = restriction ? coarse_level : fine_level;
    const bool single_patch = input_level->getNumberOfPatches() == 1 && output_level->getNumberOfPatches() == 1;
    if (test_db->getBoolWithDefault("partitioned", false) && fine_level->getNumberOfPatches() <= 1)
    {
        TBOX_ERROR("The partitioned transfer case requires multiple fine patches.\n");
    }
    Pointer<SideVariable<NDIM, double>> velocity = new SideVariable<NDIM, double>("velocity");
    Pointer<CellVariable<NDIM, double>> pressure = new CellVariable<NDIM, double>("pressure", depth);
    Pointer<CellVariable<NDIM, double>> transferred = new CellVariable<NDIM, double>("transferred", depth);
    const IntVector<NDIM> level_ratio = fine_level->getRatio() / coarse_level->getRatio();
    const int required_width = std::max(1, level_ratio.max() / 2);
    const int u_idx = var_db->registerVariableAndContext(velocity, ctx, IntVector<NDIM>(1));
    const int p_idx = var_db->registerVariableAndContext(pressure, ctx, IntVector<NDIM>(required_width));
    const int result_idx = var_db->registerVariableAndContext(transferred, ctx, IntVector<NDIM>(1));
    input_level->allocatePatchData(u_idx);
    input_level->allocatePatchData(p_idx);
    output_level->allocatePatchData(u_idx);
    output_level->allocatePatchData(result_idx);
    Pointer<Patch<NDIM>> input_patch = input_level->getPatch(0);
    Pointer<Patch<NDIM>> output_patch = output_level->getPatch(0);
    for (PatchLevel<NDIM>::Iterator p(input_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = input_level->getPatch(p());
        Pointer<CellData<NDIM, double>> data = patch->getPatchData(p_idx);
        Pointer<SideData<NDIM, double>> u = patch->getPatchData(u_idx);
        u->fillAll(0.0);
        set_pressure_profile(data, patch, test_db);
    }
    for (PatchLevel<NDIM>::Iterator p(output_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = output_level->getPatch(p());
        Pointer<CellData<NDIM, double>> data = patch->getPatchData(result_idx);
        Pointer<SideData<NDIM, double>> u = patch->getPatchData(u_idx);
        u->fillAll(0.0);
        data->fillAll(0.0);
    }

    double write_range_error = 0.0;
    Mat transfer = prolong_mat;
    int ierr = 0;
    if (restriction)
    {
        ierr = MatTranspose(prolong_mat, MAT_INITIAL_MATRIX, &transfer);
        IBTK_CHKERRQ(ierr);
        double volume_ratio = 1.0;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            volume_ratio /= level_ratio(axis);
        }
        ierr = MatScale(transfer, volume_ratio);
        IBTK_CHKERRQ(ierr);
        Pointer<CoarsenAlgorithm<NDIM>> algorithm = new CoarsenAlgorithm<NDIM>();
        algorithm->registerCoarsen(result_idx,
                                   p_idx,
                                   new CartCellDoubleLinearCoarsen(IntVector<NDIM>(
                                       test_db->getIntegerWithDefault("declared_width", required_width))),
                                   IntVector<NDIM>(0));
        Pointer<CoarsenSchedule<NDIM>> schedule = algorithm->createSchedule(coarse_level, fine_level);
        schedule->coarsenData();
        if (single_patch)
        {
            Pointer<CellData<NDIM, double>> p_output = output_patch->getPatchData(result_idx);
            // A direct sub-box call must reproduce the schedule and preserve all other cells.
            Pointer<CellData<NDIM, double>> saved =
                new CellData<NDIM, double>(p_output->getBox(), depth, p_output->getGhostCellWidth());
            saved->copy(*p_output);
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
                coarsen_op.coarsen(*output_patch, *input_patch, result_idx, p_idx, sub_box, level_ratio);
                for (Box<NDIM>::Iterator b(p_output->getGhostBox()); b; b++)
                {
                    for (int d = 0; d < depth; ++d)
                    {
                        const double expected = sub_box.contains(b()) ? (*saved)(b(), d) : -12345.0;
                        write_range_error = std::max(write_range_error, std::abs((*p_output)(b(), d) - expected));
                    }
                }
            }
            p_output->copy(*saved);
        }
    }
    else
    {
        Pointer<RefineAlgorithm<NDIM>> algorithm = new RefineAlgorithm<NDIM>();
        algorithm->registerRefine(result_idx, p_idx, result_idx, new CartCellDoubleLinearRefine(), nullptr);
        Pointer<RefineSchedule<NDIM>> schedule =
            algorithm->createSchedule(fine_level, Pointer<PatchLevel<NDIM>>(), coarse_ln, patch_hierarchy, nullptr);
        schedule->fillData(0.0);
    }

    Vec x = nullptr, matrix_result = nullptr, schedule_result = nullptr;
    ierr = MatCreateVecs(transfer, &x, &matrix_result);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(matrix_result, &schedule_result);
    IBTK_CHKERRQ(ierr);
    if (cell_only)
    {
        PETScVecUtilities::copyToPatchLevelVec(x, p_idx, p_dof_index_idx, input_level);
        PETScVecUtilities::copyToPatchLevelVec(schedule_result, result_idx, p_dof_index_idx, output_level);
    }
    else
    {
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            x, u_idx, u_dof_index_idx, p_idx, p_dof_index_idx, input_level);
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            schedule_result, u_idx, u_dof_index_idx, result_idx, p_dof_index_idx, output_level);
    }
    double adjoint_error = 0.0;
    if (restriction)
    {
        Vec coarse_probe = nullptr, prolonged_probe = nullptr;
        ierr = MatCreateVecs(prolong_mat, &coarse_probe, &prolonged_probe);
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
        ierr = MatMult(prolong_mat, coarse_probe, prolonged_probe);
        IBTK_CHKERRQ(ierr);
        PetscScalar fine_product = 0.0, coarse_product = 0.0;
        ierr = VecDot(prolonged_probe, x, &fine_product);
        IBTK_CHKERRQ(ierr);
        ierr = VecDot(coarse_probe, schedule_result, &coarse_product);
        IBTK_CHKERRQ(ierr);
        Pointer<CartesianPatchGeometry<NDIM>> fine_geom = input_patch->getPatchGeometry();
        Pointer<CartesianPatchGeometry<NDIM>> coarse_geom = output_patch->getPatchGeometry();
        for (int axis = 0; axis < NDIM; ++axis)
        {
            fine_product *= fine_geom->getDx()[axis];
            coarse_product *= coarse_geom->getDx()[axis];
        }
        adjoint_error = std::abs(fine_product - coarse_product);
        ierr = VecDestroy(&coarse_probe);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&prolonged_probe);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatMult(transfer, x, matrix_result);
    IBTK_CHKERRQ(ierr);
    double input_norm = 0.0, matrix_norm = 0.0, schedule_norm = 0.0, error = 0.0;
    ierr = VecNorm(x, NORM_INFINITY, &input_norm);
    IBTK_CHKERRQ(ierr);
    ierr = VecNorm(matrix_result, NORM_INFINITY, &matrix_norm);
    IBTK_CHKERRQ(ierr);
    ierr = VecNorm(schedule_result, NORM_INFINITY, &schedule_norm);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(matrix_result, -1.0, schedule_result);
    IBTK_CHKERRQ(ierr);
    ierr = VecNorm(matrix_result, NORM_INFINITY, &error);
    IBTK_CHKERRQ(ierr);
    plog << std::setprecision(12) << "input norm = " << input_norm << '\n'
         << "matrix norm = " << matrix_norm << '\n'
         << "schedule norm = " << schedule_norm << '\n'
         << "transfer error = " << error << '\n';
    if (restriction)
    {
        if (single_patch)
        {
            plog << "write range error = " << write_range_error << '\n';
        }
        plog << "volume adjoint error = " << adjoint_error << '\n';
    }
    const bool success = std::isfinite(adjoint_error) && adjoint_error < 1.0e-12 && std::isfinite(write_range_error) &&
                         write_range_error < 1.0e-12 && std::isfinite(error) && error < 1.0e-12 &&
                         input_norm > 1.0e-8 && matrix_norm > 1.0e-8 && schedule_norm > 1.0e-8;
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&matrix_result);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&schedule_result);
    IBTK_CHKERRQ(ierr);
    if (restriction)
    {
        ierr = MatDestroy(&transfer);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatDestroy(&prolong_mat);
    IBTK_CHKERRQ(ierr);
    ierr = AODestroy(&coarse_level_ao);
    IBTK_CHKERRQ(ierr);
    return success ? 0 : 1;
}
