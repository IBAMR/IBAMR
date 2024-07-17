// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/StaggeredStokesBoxRelaxationFACOperator.h"
#include "ibamr/StaggeredStokesFACPreconditionerStrategy.h"

#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/IBTK_CHKERRQ.h"

#include "ArrayData.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "BoxList.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "ProcessorMapping.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscsys.h"
#include "petscvec.h"

#include <algorithm>
#include <array>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int GHOSTS = 1;

inline int
compute_side_index(const hier::IndexNd& i, const BoxNd& box, const unsigned int axis)
{
    const BoxNd side_box = SideGeometryNd::toSideBox(box, axis);
    if (!side_box.contains(i)) return -1;
    int offset = 0;
    for (unsigned int d = 0; d < axis; ++d)
    {
        offset += SideGeometryNd::toSideBox(box, d).size();
    }
    return offset + side_box.offset(i);
} // compute_side_index

inline int
compute_cell_index(const hier::IndexNd& i, const BoxNd& box)
{
    if (!box.contains(i)) return -1;
    int offset = 0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        offset += SideGeometryNd::toSideBox(box, axis).size();
    }
    return box.offset(i) + offset;
} // compute_cell_index

void
buildBoxOperator(Mat& A,
                 const PoissonSpecifications& U_problem_coefs,
                 const BoxNd& box,
                 const BoxNd& ghost_box,
                 const std::array<double, NDIM>& dx)
{
    int ierr;

    const double C = U_problem_coefs.getCConstant();
    const double D = U_problem_coefs.getDConstant();

    // Allocate a PETSc matrix for the box operator.
    std::array<BoxNd, NDIM> side_boxes;
    std::array<BoxListNd, NDIM> side_ghost_boxes;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        side_boxes[axis] = SideGeometryNd::toSideBox(box, axis);
        side_ghost_boxes[axis] = SideGeometryNd::toSideBox(ghost_box, axis);
        side_ghost_boxes[axis].removeIntersections(side_boxes[axis]);
    }
    BoxListNd cell_ghost_boxes(ghost_box);
    cell_ghost_boxes.removeIntersections(box);

    int size = 0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        size += SideGeometryNd::toSideBox(ghost_box, axis).size();
    }
    size += ghost_box.size();

    static const int U_stencil_sz = 2 * NDIM + 3;
    static const int P_stencil_sz = 2 * NDIM + 1;
    std::vector<int> nnz(size, 0);

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (BoxNd::Iterator b(side_boxes[axis]); b; b++)
        {
            nnz[compute_side_index(b(), ghost_box, axis)] = std::min(size, U_stencil_sz);
        }
        for (BoxListNd::Iterator bl(side_ghost_boxes[axis]); bl; bl++)
        {
            for (BoxNd::Iterator b(bl()); b; b++)
            {
                nnz[compute_side_index(b(), ghost_box, axis)] = 1;
            }
        }
    }

    for (BoxNd::Iterator b(box); b; b++)
    {
        nnz[compute_cell_index(b(), ghost_box)] = std::min(size, P_stencil_sz);
    }
    for (BoxListNd::Iterator bl(cell_ghost_boxes); bl; bl++)
    {
        for (BoxNd::Iterator b(bl()); b; b++)
        {
            nnz[compute_cell_index(b(), ghost_box)] = 1;
        }
    }

    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, 0, size ? &nnz[0] : nullptr, &A);
    IBTK_CHKERRQ(ierr);

// Set some general matrix options.
#if !defined(NDEBUG)
    ierr = MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
#endif

    // Set the matrix coefficients to correspond to the standard finite
    // difference approximation to the time-dependent incompressible Stokes
    // operator.
    //
    // Note that boundary conditions at both physical boundaries and at
    // coarse-fine interfaces are implicitly treated by setting ghost cell
    // values appropriately.  Thus the matrix coefficients are independent of
    // any boundary conditions.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (BoxNd::Iterator b(side_boxes[axis]); b; b++)
        {
            hier::IndexNd i = b();
            const int mat_row = compute_side_index(i, ghost_box, axis);

            std::vector<int> mat_cols(U_stencil_sz, -1);
            std::vector<double> mat_vals(U_stencil_sz, 0.0);

            mat_cols[0] = mat_row;
            mat_vals[0] = C;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                hier::IndexNd shift = 0;
                shift(d) = 1;
                const hier::IndexNd u_left = i - shift;
                const hier::IndexNd u_rght = i + shift;
                mat_cols[2 * d + 1] = compute_side_index(u_left, ghost_box, axis);
                mat_cols[2 * d + 2] = compute_side_index(u_rght, ghost_box, axis);

                mat_vals[0] -= 2.0 * D / (dx[d] * dx[d]);
                mat_vals[2 * d + 1] = D / (dx[d] * dx[d]);
                mat_vals[2 * d + 2] = D / (dx[d] * dx[d]);
            }

            hier::IndexNd shift = 0;
            shift(axis) = 1;
            const hier::IndexNd p_left = i - shift;
            const hier::IndexNd p_rght = i;
            mat_cols[2 * NDIM + 1] = compute_cell_index(p_left, ghost_box);
            mat_cols[2 * NDIM + 2] = compute_cell_index(p_rght, ghost_box);

            mat_vals[2 * NDIM + 1] = -1.0 / dx[axis];
            mat_vals[2 * NDIM + 2] = 1.0 / dx[axis];

            static const int m = 1;
            static const int n = U_stencil_sz;
            ierr = MatSetValues(A, m, &mat_row, n, &mat_cols[0], &mat_vals[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    for (BoxNd::Iterator b(box); b; b++)
    {
        hier::IndexNd i = b();
        const int mat_row = compute_cell_index(i, ghost_box);

        std::vector<int> mat_cols(P_stencil_sz, -1);
        std::vector<double> mat_vals(P_stencil_sz, 0.0);

        mat_cols[0] = mat_row;
        mat_vals[0] = 0.0;

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            hier::IndexNd shift = 0;
            shift(axis) = 1;
            const hier::IndexNd u_left = i;
            const hier::IndexNd u_rght = i + shift;
            mat_cols[2 * axis + 1] = compute_side_index(u_left, ghost_box, axis);
            mat_cols[2 * axis + 2] = compute_side_index(u_rght, ghost_box, axis);

            mat_vals[2 * axis + 1] = 1.0 / dx[axis];
            mat_vals[2 * axis + 2] = -1.0 / dx[axis];
        }

        static const int m = 1;
        static const int n = P_stencil_sz;
        ierr = MatSetValues(A, m, &mat_row, n, &mat_cols[0], &mat_vals[0], INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }

    // Set the entries in the ghost cell region so that ghost cell values are
    // not modified by the smoother.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (BoxListNd::Iterator bl(side_ghost_boxes[axis]); bl; bl++)
        {
            for (BoxNd::Iterator b(bl()); b; b++)
            {
                const int i = compute_side_index(b(), ghost_box, axis);
                ierr = MatSetValue(A, i, i, 1.0, INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    for (BoxListNd::Iterator bl(cell_ghost_boxes); bl; bl++)
    {
        for (BoxNd::Iterator b(bl()); b; b++)
        {
            const int i = compute_cell_index(b(), ghost_box);
            ierr = MatSetValue(A, i, i, 1.0, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrices.
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return;
} // buildBoxOperator

void
modifyRhsForBcs(Vec& v,
                const SideDataNd<double>& U_data,
                const CellDataNd<double>& P_data,
                const PoissonSpecifications& U_problem_coefs,
                const BoxNd& box,
                const BoxNd& ghost_box,
                const double* const dx)
{
    int ierr;

    const double D = U_problem_coefs.getDConstant();
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const BoxNd side_box = SideGeometryNd::toSideBox(box, axis);
        for (BoxNd::Iterator b(side_box); b; b++)
        {
            hier::IndexNd i = b();
            const int idx = compute_side_index(i, ghost_box, axis);

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                hier::IndexNd shift = 0;
                shift(d) = 1;
                const hier::IndexNd u_left = i - shift;
                const hier::IndexNd u_rght = i + shift;
                if (!side_box.contains(u_left))
                {
                    ierr = VecSetValue(v,
                                       idx,
                                       +D * U_data(SideIndexNd(u_left, axis, SideIndexNd::Lower)) / (dx[d] * dx[d]),
                                       ADD_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
                if (!side_box.contains(u_rght))
                {
                    ierr = VecSetValue(v,
                                       idx,
                                       +D * U_data(SideIndexNd(u_rght, axis, SideIndexNd::Lower)) / (dx[d] * dx[d]),
                                       ADD_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }

            hier::IndexNd shift = 0;
            shift(axis) = 1;
            const hier::IndexNd p_left = i - shift;
            const hier::IndexNd p_rght = i;
            if (!box.contains(p_left))
            {
                ierr = VecSetValue(v, idx, +P_data(p_left) / dx[axis], ADD_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            if (!box.contains(p_rght))
            {
                ierr = VecSetValue(v, idx, -P_data(p_rght) / dx[axis], ADD_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    ierr = VecAssemblyBegin(v);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);
    IBTK_CHKERRQ(ierr);
    return;
} // modifyRhsForBcs

inline void
copyToVec(Vec& v,
          const SideDataNd<double>& U_data,
          const CellDataNd<double>& P_data,
          const BoxNd& box,
          const BoxNd& ghost_box)
{
    int ierr;

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const BoxNd side_box = SideGeometryNd::toSideBox(box, axis);
        for (BoxNd::Iterator b(side_box); b; b++)
        {
            const hier::IndexNd& i = b();
            const SideIndexNd s_i(i, axis, 0);
            const int idx = compute_side_index(i, ghost_box, axis);
            ierr = VecSetValue(v, idx, U_data(s_i), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    for (BoxNd::Iterator b(box); b; b++)
    {
        const hier::IndexNd& i = b();
        const int idx = compute_cell_index(i, ghost_box);
        ierr = VecSetValue(v, idx, P_data(i), INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(v);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);
    IBTK_CHKERRQ(ierr);
    return;
} // copyToVec

inline void
copyFromVec(Vec& v, SideDataNd<double>& U_data, CellDataNd<double>& P_data, const BoxNd& box, const BoxNd& ghost_box)
{
    int ierr;

    const double omega = 0.65;

    double U;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const BoxNd side_box = SideGeometryNd::toSideBox(box, axis);
        for (BoxNd::Iterator b(side_box); b; b++)
        {
            const hier::IndexNd& i = b();
            const SideIndexNd s_i(i, axis, SideIndexNd::Lower);
            const int idx = compute_side_index(i, ghost_box, axis);
            ierr = VecGetValues(v, 1, &idx, &U);
            IBTK_CHKERRQ(ierr);
            U_data(s_i) = (1.0 - omega) * U_data(s_i) + omega * U;
        }
    }

    double P;
    for (BoxNd::Iterator b(box); b; b++)
    {
        const hier::IndexNd& i = b();
        const int idx = compute_cell_index(i, ghost_box);
        ierr = VecGetValues(v, 1, &idx, &P);
        IBTK_CHKERRQ(ierr);
        P_data(i) = (1.0 - omega) * P_data(i) + omega * P;
    }
    return;
} // copyFromVec
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesBoxRelaxationFACOperator::StaggeredStokesBoxRelaxationFACOperator(
    const std::string& object_name,
    const Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : StaggeredStokesFACPreconditionerStrategy(object_name, GHOSTS, input_db, default_options_prefix)
{
    // intentionally blank
    return;
} // StaggeredStokesBoxRelaxationFACOperator

StaggeredStokesBoxRelaxationFACOperator::~StaggeredStokesBoxRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~StaggeredStokesBoxRelaxationFACOperator

void
StaggeredStokesBoxRelaxationFACOperator::smoothError(SAMRAIVectorRealNd<double>& error,
                                                     const SAMRAIVectorRealNd<double>& residual,
                                                     int level_num,
                                                     int num_sweeps,
                                                     bool /*performing_pre_sweeps*/,
                                                     bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    int ierr;
    Pointer<PatchLevelNd> level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);
    const int U_scratch_idx = d_side_scratch_idx;
    const int P_scratch_idx = d_cell_scratch_idx;

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevelNd::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<PatchNd> patch = level->getPatch(p());

            Pointer<SideDataNd<double> > U_error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideDataNd<double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
#if !defined(NDEBUG)
            const BoxNd& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
            TBOX_ASSERT(U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                U_scratch_data->getArrayData(axis).copy(U_error_data->getArrayData(axis),
                                                        d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                        IntVectorNd(0));
            }

            Pointer<CellDataNd<double> > P_error_data = error.getComponentPatchData(1, *patch);
            Pointer<CellDataNd<double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
#if !defined(NDEBUG)
            const BoxNd& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
            TBOX_ASSERT(P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            P_scratch_data->getArrayData().copy(
                P_error_data->getArrayData(), d_patch_cell_bc_box_overlap[level_num][patch_counter], IntVectorNd(0));
        }
    }

    // Smooth the error by the specified number of sweeps.
    for (int isweep = 0; isweep < num_sweeps; ++isweep)
    {
        // Re-fill ghost cell data as needed.
        if (level_num > d_coarsest_ln)
        {
            if (isweep > 0)
            {
                // Copy the coarse-fine interface ghost cell values which are
                // cached in the scratch data into the error data.
                int patch_counter = 0;
                for (PatchLevelNd::Iterator p(level); p; p++, ++patch_counter)
                {
                    Pointer<PatchNd> patch = level->getPatch(p());

                    Pointer<SideDataNd<double> > U_error_data = error.getComponentPatchData(0, *patch);
                    Pointer<SideDataNd<double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
#if !defined(NDEBUG)
                    const BoxNd& U_ghost_box = U_error_data->getGhostBox();
                    TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
                    TBOX_ASSERT(U_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        U_error_data->getArrayData(axis).copy(
                            U_scratch_data->getArrayData(axis),
                            d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                            IntVectorNd(0));
                    }

                    Pointer<CellDataNd<double> > P_error_data = error.getComponentPatchData(1, *patch);
                    Pointer<CellDataNd<double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
#if !defined(NDEBUG)
                    const BoxNd& P_ghost_box = P_error_data->getGhostBox();
                    TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
                    TBOX_ASSERT(P_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    P_error_data->getArrayData().copy(P_scratch_data->getArrayData(),
                                                      d_patch_cell_bc_box_overlap[level_num][patch_counter],
                                                      IntVectorNd(0));
                }

                // Fill the non-coarse-fine interface ghost cell values.
                const std::pair<int, int> error_idxs = std::make_pair(U_error_idx, P_error_idx);
                xeqScheduleGhostFillNoCoarse(error_idxs, level_num);
            }

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension.
            d_U_cf_bdry_op->setPatchDataIndex(U_error_idx);
            d_P_cf_bdry_op->setPatchDataIndex(P_error_idx);
            const IntVectorNd& ratio = level->getRatioToCoarserLevel();
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                const IntVectorNd& ghost_width_to_fill = d_gcw;
                d_U_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
                d_P_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
            }
        }
        else if (isweep > 0)
        {
            const std::pair<int, int> error_idxs = std::make_pair(U_error_idx, P_error_idx);
            xeqScheduleGhostFillNoCoarse(error_idxs, level_num);
        }

        // Smooth the error on the patches.
        Vec& e = d_box_e[level_num];
        Vec& r = d_box_r[level_num];
        KSP& ksp = d_box_ksp[level_num];
        int patch_counter = 0;
        for (PatchLevelNd::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<PatchNd> patch = level->getPatch(p());
            Pointer<SideDataNd<double> > U_error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideDataNd<double> > U_residual_data = residual.getComponentPatchData(0, *patch);
#if !defined(NDEBUG)
            const BoxNd& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_residual_data->getGhostBox());
            TBOX_ASSERT(U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_residual_data->getGhostCellWidth() == d_gcw);
#endif
            Pointer<CellDataNd<double> > P_error_data = error.getComponentPatchData(1, *patch);
            Pointer<CellDataNd<double> > P_residual_data = residual.getComponentPatchData(1, *patch);
#if !defined(NDEBUG)
            const BoxNd& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_residual_data->getGhostBox());
            TBOX_ASSERT(P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_residual_data->getGhostCellWidth() == d_gcw);
#endif
            // Smooth the error on the patch.
            const BoxNd& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometryNd> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            for (BoxNd::Iterator b(patch_box); b; b++)
            {
                const hier::IndexNd& i = b();
                const BoxNd box(i, i);
                copyToVec(e, *U_error_data, *P_error_data, box, box);
                copyToVec(r, *U_residual_data, *P_residual_data, box, box);
                modifyRhsForBcs(r, *U_error_data, *P_error_data, d_U_problem_coefs, box, box, dx);
                ierr = KSPSolve(ksp, r, e);
                IBTK_CHKERRQ(ierr);
                copyFromVec(e, *U_error_data, *P_error_data, box, box);
            }
        }
    }

    // Synchronize data along patch boundaries.
    xeqScheduleDataSynch(U_error_idx, level_num);
    return;
} // smoothError

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredStokesBoxRelaxationFACOperator::initializeOperatorStateSpecialized(const SAMRAIVectorRealNd<double>&
                                                                            /*solution*/,
                                                                            const SAMRAIVectorRealNd<double>&
                                                                            /*rhs*/,
                                                                            const int coarsest_reset_ln,
                                                                            const int finest_reset_ln)
{
    // Initialize the box relaxation data on each level of the patch hierarchy.
    d_box_op.resize(d_finest_ln + 1);
    d_box_e.resize(d_finest_ln + 1);
    d_box_r.resize(d_finest_ln + 1);
    d_box_ksp.resize(d_finest_ln + 1);
    const BoxNd box(hier::IndexNd(0), hier::IndexNd(0));
    Pointer<CartesianGridGeometryNd> geometry = d_hierarchy->getGridGeometry();
    const double* const dx_coarsest = geometry->getDx();
    std::array<double, NDIM> dx;
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        const IntVectorNd& ratio = d_hierarchy->getPatchLevel(ln)->getRatio();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx_coarsest[d] / static_cast<double>(ratio(d));
        }
        buildBoxOperator(d_box_op[ln], d_U_problem_coefs, box, box, dx);
        int ierr;
        ierr = MatCreateVecs(d_box_op[ln], &d_box_e[ln], &d_box_r[ln]);
        IBTK_CHKERRQ(ierr);
        ierr = KSPCreate(PETSC_COMM_SELF, &d_box_ksp[ln]);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetOperators(d_box_ksp[ln], d_box_op[ln], d_box_op[ln]);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetReusePreconditioner(d_box_ksp[ln], PETSC_TRUE);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetType(d_box_ksp[ln], KSPPREONLY);
        IBTK_CHKERRQ(ierr);
        PC box_pc;
        ierr = KSPGetPC(d_box_ksp[ln], &box_pc);
        IBTK_CHKERRQ(ierr);
        ierr = PCSetType(box_pc, PCLU);
        IBTK_CHKERRQ(ierr);
        ierr = PCFactorReorderForNonzeroDiagonal(box_pc, std::numeric_limits<double>::epsilon());
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetUp(d_box_ksp[ln]);
        IBTK_CHKERRQ(ierr);
    }

    // Get overlap information for setting patch boundary conditions.
    d_patch_side_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevelNd> level = d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_side_bc_box_overlap[ln].resize(num_local_patches);
        int patch_counter = 0;
        for (PatchLevelNd::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<PatchNd> patch = level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const BoxNd side_box = SideGeometryNd::toSideBox(patch_box, axis);
                const BoxNd side_ghost_box = BoxNd::grow(side_box, 1);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis] = BoxListNd(side_ghost_box);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis].removeIntersections(side_box);
            }
        }
    }

    d_patch_cell_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevelNd> level = d_hierarchy->getPatchLevel(ln);

        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_cell_bc_box_overlap[ln].resize(num_local_patches);

        int patch_counter = 0;
        for (PatchLevelNd::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<PatchNd> patch = level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();
            const BoxNd& ghost_box = BoxNd::grow(patch_box, 1);

            d_patch_cell_bc_box_overlap[ln][patch_counter] = BoxListNd(ghost_box);
            d_patch_cell_bc_box_overlap[ln][patch_counter].removeIntersections(patch_box);
        }
    }
    return;
} // initializeOperatorStateSpecialized

void
StaggeredStokesBoxRelaxationFACOperator::deallocateOperatorStateSpecialized(const int coarsest_reset_ln,
                                                                            const int finest_reset_ln)
{
    if (!d_is_initialized) return;
    for (int ln = coarsest_reset_ln; ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        int ierr;
        ierr = MatDestroy(&d_box_op[ln]);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&d_box_e[ln]);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&d_box_r[ln]);
        IBTK_CHKERRQ(ierr);
        ierr = KSPDestroy(&d_box_ksp[ln]);
        IBTK_CHKERRQ(ierr);
        d_patch_side_bc_box_overlap[ln].resize(0);
        d_patch_cell_bc_box_overlap[ln].resize(0);
    }
    return;
} // deallocateOperatorStateSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
