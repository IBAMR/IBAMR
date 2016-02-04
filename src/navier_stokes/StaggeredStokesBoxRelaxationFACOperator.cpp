// Filename: StaggeredStokesBoxRelaxationFACOperator.cpp
// Created on 11 Jun 2010 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <algorithm>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

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
#include "boost/array.hpp"
#include "ibamr/StaggeredStokesBoxRelaxationFACOperator.h"
#include "ibamr/StaggeredStokesFACPreconditionerStrategy.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int GHOSTS = 1;

inline int
compute_side_index(const Index<NDIM>& i, const Box<NDIM>& box, const unsigned int axis)
{
    const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(box, axis);
    if (!side_box.contains(i)) return -1;
    int offset = 0;
    for (unsigned int d = 0; d < axis; ++d)
    {
        offset += SideGeometry<NDIM>::toSideBox(box, d).size();
    }
    return offset + side_box.offset(i);
} // compute_side_index

inline int
compute_cell_index(const Index<NDIM>& i, const Box<NDIM>& box)
{
    if (!box.contains(i)) return -1;
    int offset = 0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        offset += SideGeometry<NDIM>::toSideBox(box, axis).size();
    }
    return box.offset(i) + offset;
} // compute_cell_index

void
buildBoxOperator(Mat& A,
                 const PoissonSpecifications& U_problem_coefs,
                 const Box<NDIM>& box,
                 const Box<NDIM>& ghost_box,
                 const boost::array<double, NDIM>& dx)
{
    int ierr;

    const double C = U_problem_coefs.getCConstant();
    const double D = U_problem_coefs.getDConstant();

    // Allocate a PETSc matrix for the box operator.
    boost::array<Box<NDIM>, NDIM> side_boxes;
    boost::array<BoxList<NDIM>, NDIM> side_ghost_boxes;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        side_boxes[axis] = SideGeometry<NDIM>::toSideBox(box, axis);
        side_ghost_boxes[axis] = SideGeometry<NDIM>::toSideBox(ghost_box, axis);
        side_ghost_boxes[axis].removeIntersections(side_boxes[axis]);
    }
    BoxList<NDIM> cell_ghost_boxes(ghost_box);
    cell_ghost_boxes.removeIntersections(box);

    int size = 0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        size += SideGeometry<NDIM>::toSideBox(ghost_box, axis).size();
    }
    size += ghost_box.size();

    static const int U_stencil_sz = 2 * NDIM + 3;
    static const int P_stencil_sz = 2 * NDIM + 1;
    std::vector<int> nnz(size, 0);

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (Box<NDIM>::Iterator b(side_boxes[axis]); b; b++)
        {
            nnz[compute_side_index(b(), ghost_box, axis)] = std::min(size, U_stencil_sz);
        }
        for (BoxList<NDIM>::Iterator bl(side_ghost_boxes[axis]); bl; bl++)
        {
            for (Box<NDIM>::Iterator b(bl()); b; b++)
            {
                nnz[compute_side_index(b(), ghost_box, axis)] = 1;
            }
        }
    }

    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        nnz[compute_cell_index(b(), ghost_box)] = std::min(size, P_stencil_sz);
    }
    for (BoxList<NDIM>::Iterator bl(cell_ghost_boxes); bl; bl++)
    {
        for (Box<NDIM>::Iterator b(bl()); b; b++)
        {
            nnz[compute_cell_index(b(), ghost_box)] = 1;
        }
    }

    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, 0, size ? &nnz[0] : NULL, &A);
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
        for (Box<NDIM>::Iterator b(side_boxes[axis]); b; b++)
        {
            Index<NDIM> i = b();
            const int mat_row = compute_side_index(i, ghost_box, axis);

            std::vector<int> mat_cols(U_stencil_sz, -1);
            std::vector<double> mat_vals(U_stencil_sz, 0.0);

            mat_cols[0] = mat_row;
            mat_vals[0] = C;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                Index<NDIM> shift = 0;
                shift(d) = 1;
                const Index<NDIM> u_left = i - shift;
                const Index<NDIM> u_rght = i + shift;
                mat_cols[2 * d + 1] = compute_side_index(u_left, ghost_box, axis);
                mat_cols[2 * d + 2] = compute_side_index(u_rght, ghost_box, axis);

                mat_vals[0] -= 2.0 * D / (dx[d] * dx[d]);
                mat_vals[2 * d + 1] = D / (dx[d] * dx[d]);
                mat_vals[2 * d + 2] = D / (dx[d] * dx[d]);
            }

            Index<NDIM> shift = 0;
            shift(axis) = 1;
            const Index<NDIM> p_left = i - shift;
            const Index<NDIM> p_rght = i;
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

    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        Index<NDIM> i = b();
        const int mat_row = compute_cell_index(i, ghost_box);

        std::vector<int> mat_cols(P_stencil_sz, -1);
        std::vector<double> mat_vals(P_stencil_sz, 0.0);

        mat_cols[0] = mat_row;
        mat_vals[0] = 0.0;

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            Index<NDIM> shift = 0;
            shift(axis) = 1;
            const Index<NDIM> u_left = i;
            const Index<NDIM> u_rght = i + shift;
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
        for (BoxList<NDIM>::Iterator bl(side_ghost_boxes[axis]); bl; bl++)
        {
            for (Box<NDIM>::Iterator b(bl()); b; b++)
            {
                const int i = compute_side_index(b(), ghost_box, axis);
                ierr = MatSetValue(A, i, i, 1.0, INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    for (BoxList<NDIM>::Iterator bl(cell_ghost_boxes); bl; bl++)
    {
        for (Box<NDIM>::Iterator b(bl()); b; b++)
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
                const SideData<NDIM, double>& U_data,
                const CellData<NDIM, double>& P_data,
                const PoissonSpecifications& U_problem_coefs,
                const Box<NDIM>& box,
                const Box<NDIM>& ghost_box,
                const double* const dx)
{
    int ierr;

    const double D = U_problem_coefs.getDConstant();
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            Index<NDIM> i = b();
            const int idx = compute_side_index(i, ghost_box, axis);

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                Index<NDIM> shift = 0;
                shift(d) = 1;
                const Index<NDIM> u_left = i - shift;
                const Index<NDIM> u_rght = i + shift;
                if (!side_box.contains(u_left))
                {
                    ierr = VecSetValue(v,
                                       idx,
                                       +D * U_data(SideIndex<NDIM>(u_left, axis, SideIndex<NDIM>::Lower)) /
                                           (dx[d] * dx[d]),
                                       ADD_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
                if (!side_box.contains(u_rght))
                {
                    ierr = VecSetValue(v,
                                       idx,
                                       +D * U_data(SideIndex<NDIM>(u_rght, axis, SideIndex<NDIM>::Lower)) /
                                           (dx[d] * dx[d]),
                                       ADD_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }

            Index<NDIM> shift = 0;
            shift(axis) = 1;
            const Index<NDIM> p_left = i - shift;
            const Index<NDIM> p_rght = i;
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
          const SideData<NDIM, double>& U_data,
          const CellData<NDIM, double>& P_data,
          const Box<NDIM>& box,
          const Box<NDIM>& ghost_box)
{
    int ierr;

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const Index<NDIM>& i = b();
            const SideIndex<NDIM> s_i(i, axis, 0);
            const int idx = compute_side_index(i, ghost_box, axis);
            ierr = VecSetValue(v, idx, U_data(s_i), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        const Index<NDIM>& i = b();
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
copyFromVec(Vec& v,
            SideData<NDIM, double>& U_data,
            CellData<NDIM, double>& P_data,
            const Box<NDIM>& box,
            const Box<NDIM>& ghost_box)
{
    int ierr;

    const double omega = 0.65;

    double U;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const Index<NDIM>& i = b();
            const SideIndex<NDIM> s_i(i, axis, SideIndex<NDIM>::Lower);
            const int idx = compute_side_index(i, ghost_box, axis);
            ierr = VecGetValues(v, 1, &idx, &U);
            IBTK_CHKERRQ(ierr);
            U_data(s_i) = (1.0 - omega) * U_data(s_i) + omega * U;
        }
    }

    double P;
    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        const Index<NDIM>& i = b();
        const int idx = compute_cell_index(i, ghost_box);
        ierr = VecGetValues(v, 1, &idx, &P);
        IBTK_CHKERRQ(ierr);
        P_data(i) = (1.0 - omega) * P_data(i) + omega * P;
    }
    return;
} // copyFromVec
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesBoxRelaxationFACOperator::StaggeredStokesBoxRelaxationFACOperator(
    const std::string& object_name,
    const Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : StaggeredStokesFACPreconditionerStrategy(object_name, GHOSTS, input_db, default_options_prefix),
      d_box_op(),
      d_box_e(),
      d_box_r(),
      d_box_ksp(),
      d_patch_side_bc_box_overlap(),
      d_patch_cell_bc_box_overlap()
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
StaggeredStokesBoxRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
                                                     const SAMRAIVectorReal<NDIM, double>& residual,
                                                     int level_num,
                                                     int num_sweeps,
                                                     bool /*performing_pre_sweeps*/,
                                                     bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    int ierr;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);
    const int U_scratch_idx = d_side_scratch_idx;
    const int P_scratch_idx = d_cell_scratch_idx;

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > U_error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM, double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
#if !defined(NDEBUG)
            const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
            TBOX_ASSERT(U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                U_scratch_data->getArrayData(axis).copy(U_error_data->getArrayData(axis),
                                                        d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                        IntVector<NDIM>(0));
            }

            Pointer<CellData<NDIM, double> > P_error_data = error.getComponentPatchData(1, *patch);
            Pointer<CellData<NDIM, double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
#if !defined(NDEBUG)
            const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
            TBOX_ASSERT(P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            P_scratch_data->getArrayData().copy(P_error_data->getArrayData(),
                                                d_patch_cell_bc_box_overlap[level_num][patch_counter],
                                                IntVector<NDIM>(0));
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
                for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());

                    Pointer<SideData<NDIM, double> > U_error_data = error.getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM, double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
#if !defined(NDEBUG)
                    const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
                    TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
                    TBOX_ASSERT(U_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        U_error_data->getArrayData(axis)
                            .copy(U_scratch_data->getArrayData(axis),
                                  d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                  IntVector<NDIM>(0));
                    }

                    Pointer<CellData<NDIM, double> > P_error_data = error.getComponentPatchData(1, *patch);
                    Pointer<CellData<NDIM, double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
#if !defined(NDEBUG)
                    const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
                    TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
                    TBOX_ASSERT(P_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    P_error_data->getArrayData().copy(P_scratch_data->getArrayData(),
                                                      d_patch_cell_bc_box_overlap[level_num][patch_counter],
                                                      IntVector<NDIM>(0));
                }

                // Fill the non-coarse-fine interface ghost cell values.
                const std::pair<int, int> error_idxs = std::make_pair(U_error_idx, P_error_idx);
                xeqScheduleGhostFillNoCoarse(error_idxs, level_num);
            }

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension.
            d_U_cf_bdry_op->setPatchDataIndex(U_error_idx);
            d_P_cf_bdry_op->setPatchDataIndex(P_error_idx);
            const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const IntVector<NDIM>& ghost_width_to_fill = d_gcw;
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
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > U_error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM, double> > U_residual_data = residual.getComponentPatchData(0, *patch);
#if !defined(NDEBUG)
            const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_residual_data->getGhostBox());
            TBOX_ASSERT(U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_residual_data->getGhostCellWidth() == d_gcw);
#endif
            Pointer<CellData<NDIM, double> > P_error_data = error.getComponentPatchData(1, *patch);
            Pointer<CellData<NDIM, double> > P_residual_data = residual.getComponentPatchData(1, *patch);
#if !defined(NDEBUG)
            const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_residual_data->getGhostBox());
            TBOX_ASSERT(P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_residual_data->getGhostCellWidth() == d_gcw);
#endif
            // Smooth the error on the patch.
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            for (Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const Index<NDIM>& i = b();
                const Box<NDIM> box(i, i);
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
StaggeredStokesBoxRelaxationFACOperator::initializeOperatorStateSpecialized(const SAMRAIVectorReal<NDIM, double>&
                                                                            /*solution*/,
                                                                            const SAMRAIVectorReal<NDIM, double>&
                                                                            /*rhs*/,
                                                                            const int coarsest_reset_ln,
                                                                            const int finest_reset_ln)
{
    // Initialize the box relaxation data on each level of the patch hierarchy.
    d_box_op.resize(d_finest_ln + 1);
    d_box_e.resize(d_finest_ln + 1);
    d_box_r.resize(d_finest_ln + 1);
    d_box_ksp.resize(d_finest_ln + 1);
    const Box<NDIM> box(Index<NDIM>(0), Index<NDIM>(0));
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    const double* const dx_coarsest = geometry->getDx();
    boost::array<double, NDIM> dx;
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        const IntVector<NDIM>& ratio = d_hierarchy->getPatchLevel(ln)->getRatio();
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
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_side_bc_box_overlap[ln].resize(num_local_patches);
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                const Box<NDIM> side_ghost_box = Box<NDIM>::grow(side_box, 1);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis] = BoxList<NDIM>(side_ghost_box);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis].removeIntersections(side_box);
            }
        }
    }

    d_patch_cell_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_cell_bc_box_overlap[ln].resize(num_local_patches);

        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Box<NDIM>& ghost_box = Box<NDIM>::grow(patch_box, 1);

            d_patch_cell_bc_box_overlap[ln][patch_counter] = BoxList<NDIM>(ghost_box);
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

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
