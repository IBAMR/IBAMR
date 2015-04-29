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

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "boost/array.hpp"
#include "ibamr/StaggeredStokesBoxRelaxationFACOperator.h"
#include "ibamr/StaggeredStokesFACPreconditionerStrategy.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/ibtk_utilities.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscsys.h"
#include "petscvec.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"

#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int GHOSTS = 1;

inline int compute_side_index(const Index& i, const Box& box, const unsigned int axis)
{
    const Box side_box = SideGeometry::toSideBox(box, axis);
    if (!side_box.contains(i)) return -1;
    int offset = 0;
    for (unsigned int d = 0; d < axis; ++d)
    {
        offset += SideGeometry::toSideBox(box, d).size();
    }
    return offset + static_cast<unsigned int>(side_box.offset(i));
}

inline int compute_cell_index(const Index& i, const Box& box)
{
    if (!box.contains(i)) return -1;
    int offset = 0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        offset += SideGeometry::toSideBox(box, axis).size();
    }
    return offset + static_cast<unsigned int>(box.offset(i));
}

void buildBoxOperator(Mat& A,
                      const PoissonSpecifications& U_problem_coefs,
                      const Box& box,
                      const Box& ghost_box,
                      const boost::array<double, NDIM>& dx)
{
    int ierr;

    const double C = U_problem_coefs.getCConstant();
    const double D = U_problem_coefs.getDConstant();

    // Allocate a PETSc matrix for the box operator.
    std::vector<Box> side_boxes(NDIM, Box(DIM));
    boost::array<BoxContainer, NDIM> side_ghost_boxes;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        side_boxes[axis] = SideGeometry::toSideBox(box, axis);
        side_ghost_boxes[axis] = BoxContainer(SideGeometry::toSideBox(ghost_box, axis));
        side_ghost_boxes[axis].removeIntersections(side_boxes[axis]);
    }
    BoxContainer cell_ghost_boxes(ghost_box);
    cell_ghost_boxes.removeIntersections(box);

    int size = 0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        size += SideGeometry::toSideBox(ghost_box, axis).size();
    }
    size += ghost_box.size();

    static const int U_stencil_sz = 2 * NDIM + 3;
    static const int P_stencil_sz = 2 * NDIM + 1;
    std::vector<int> nnz(size, 0);

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (auto b = side_boxes[axis].begin(), e = side_boxes[axis].end(); b != e; ++b)
        {
            nnz[compute_side_index(*b, ghost_box, axis)] = std::min(size, U_stencil_sz);
        }
        for (auto bl = side_ghost_boxes[axis].begin(), e = side_ghost_boxes[axis].end(); bl != e; ++bl)
        {
            const Box& box = *bl;
            for (auto b = box.begin(), e = box.end(); b != e; ++b)
            {
                nnz[compute_side_index(*b, ghost_box, axis)] = 1;
            }
        }
    }

    for (auto b = box.begin(), e = box.end(); b != e; ++b)
    {
        nnz[compute_cell_index(*b, ghost_box)] = std::min(size, P_stencil_sz);
    }
    for (auto bl = cell_ghost_boxes.begin(), e = cell_ghost_boxes.end(); bl != e; ++bl)
    {
        const Box& box = *bl;
        for (auto b = box.begin(), e = box.end(); b != e; ++b)
        {
            nnz[compute_cell_index(*b, ghost_box)] = 1;
        }
    }

    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, PETSC_DEFAULT, &nnz[0], &A);
    IBTK_CHKERRQ(ierr);

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
        for (auto b = side_boxes[axis].begin(), e = side_boxes[axis].end(); b != e; ++b)
        {
            const Index& i = *b;
            const int mat_row = compute_side_index(i, ghost_box, axis);

            std::vector<int> mat_cols(U_stencil_sz, -1);
            std::vector<double> mat_vals(U_stencil_sz, 0.0);

            mat_cols[0] = mat_row;
            mat_vals[0] = C;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                Index shift = Index::getZeroIndex(DIM);
                shift(d) = 1;
                const Index u_left = i - shift;
                const Index u_rght = i + shift;
                mat_cols[2 * d + 1] = compute_side_index(u_left, ghost_box, axis);
                mat_cols[2 * d + 2] = compute_side_index(u_rght, ghost_box, axis);

                mat_vals[0] += 2.0 * D / (dx[d] * dx[d]);
                mat_vals[2 * d + 1] = -D / (dx[d] * dx[d]);
                mat_vals[2 * d + 2] = -D / (dx[d] * dx[d]);
            }

            Index shift = Index::getZeroIndex(DIM);
            shift(axis) = 1;
            const Index p_left = i - shift;
            const Index p_rght = i;
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

    for (auto b = box.begin(), e = box.end(); b != e; ++b)
    {
        const Index& i = *b;
        const int mat_row = compute_cell_index(i, ghost_box);

        std::vector<int> mat_cols(P_stencil_sz, -1);
        std::vector<double> mat_vals(P_stencil_sz, 0.0);

        mat_cols[0] = mat_row;
        mat_vals[0] = 0.0;

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            Index shift = Index::getZeroIndex(DIM);
            shift(axis) = 1;
            const Index u_left = i;
            const Index u_rght = i + shift;
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
        for (auto bl = side_ghost_boxes[axis].begin(), e = side_ghost_boxes[axis].end(); bl != e; ++bl)
        {
            const Box& box = *bl;
            for (auto b = box.begin(), e = box.end(); b != e; ++b)
            {
                const int i = compute_side_index(*b, ghost_box, axis);
                ierr = MatSetValue(A, i, i, 1.0, INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    for (auto bl = cell_ghost_boxes.begin(), e = cell_ghost_boxes.end(); bl != e; ++bl)
    {
        const Box& box = *bl;
        for (auto b = box.begin(), e = box.end(); b != e; ++b)
        {
            const int i = compute_cell_index(*b, ghost_box);
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
}

void modifyRhsForBcs(Vec& v,
                     const SideData<double>& U_data,
                     const CellData<double>& P_data,
                     const PoissonSpecifications& U_problem_coefs,
                     const Box& box,
                     const Box& ghost_box,
                     const double* const dx)
{
    int ierr;

    const double D = U_problem_coefs.getDConstant();
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const Box side_box = SideGeometry::toSideBox(box, axis);
        for (auto b = side_box.begin(), e = side_box.end(); b != e; ++b)
        {
            const Index& i = *b;
            const int idx = compute_side_index(i, ghost_box, axis);

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                Index shift = Index::getZeroIndex(DIM);
                shift(d) = 1;
                const Index u_left = i - shift;
                const Index u_rght = i + shift;
                if (!side_box.contains(u_left))
                {
                    ierr = VecSetValue(v, idx, +D * U_data(SideIndex(u_left, axis, SideIndex::Lower)) / (dx[d] * dx[d]),
                                       ADD_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
                if (!side_box.contains(u_rght))
                {
                    ierr = VecSetValue(v, idx, +D * U_data(SideIndex(u_rght, axis, SideIndex::Lower)) / (dx[d] * dx[d]),
                                       ADD_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }

            Index shift = Index::getZeroIndex(DIM);
            shift(axis) = 1;
            const CellIndex p_left(i - shift);
            const CellIndex p_rght(i);
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
}

inline void
copyToVec(Vec& v, const SideData<double>& U_data, const CellData<double>& P_data, const Box& box, const Box& ghost_box)
{
    int ierr;

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const Box side_box = SideGeometry::toSideBox(box, axis);
        for (auto b = side_box.begin(), e = side_box.end(); b != e; ++b)
        {
            const Index& i = *b;
            const SideIndex s_i(i, axis, 0);
            const int idx = compute_side_index(i, ghost_box, axis);
            ierr = VecSetValue(v, idx, U_data(s_i), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    for (auto b = CellGeometry::begin(box), e = CellGeometry::end(box); b != e; ++b)
    {
        const CellIndex& i = *b;
        const int idx = compute_cell_index(i, ghost_box);
        ierr = VecSetValue(v, idx, P_data(i), INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(v);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);
    IBTK_CHKERRQ(ierr);
    return;
}

inline void
copyFromVec(Vec& v, SideData<double>& U_data, CellData<double>& P_data, const Box& box, const Box& ghost_box)
{
    int ierr;

    const double omega = 0.65;

    double U;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const Box side_box = SideGeometry::toSideBox(box, axis);
        for (auto b = side_box.begin(), e = side_box.end(); b != e; ++b)
        {
            const Index& i = *b;
            const SideIndex s_i(i, axis, SideIndex::Lower);
            const int idx = compute_side_index(i, ghost_box, axis);
            ierr = VecGetValues(v, 1, &idx, &U);
            IBTK_CHKERRQ(ierr);
            U_data(s_i) = (1.0 - omega) * U_data(s_i) + omega * U;
        }
    }

    double P;
    for (auto b = CellGeometry::begin(box), e = CellGeometry::end(box); b != e; ++b)
    {
        const CellIndex& i = *b;
        const int idx = compute_cell_index(i, ghost_box);
        ierr = VecGetValues(v, 1, &idx, &P);
        IBTK_CHKERRQ(ierr);
        P_data(i) = (1.0 - omega) * P_data(i) + omega * P;
    }
    return;
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesBoxRelaxationFACOperator::StaggeredStokesBoxRelaxationFACOperator(
    const std::string& object_name,
    const boost::shared_ptr<Database>& input_db,
    const std::string& default_options_prefix)
    : StaggeredStokesFACPreconditionerStrategy(object_name, GHOSTS, input_db, default_options_prefix), d_box_op(),
      d_box_e(), d_box_r(), d_box_ksp(), d_patch_side_bc_box_overlap(), d_patch_cell_bc_box_overlap()
{
    // intentionally blank
    return;
}

StaggeredStokesBoxRelaxationFACOperator::~StaggeredStokesBoxRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}

void StaggeredStokesBoxRelaxationFACOperator::smoothError(SAMRAIVectorReal<double>& error,
                                                          const SAMRAIVectorReal<double>& residual,
                                                          int level_num,
                                                          int num_sweeps,
                                                          bool /*performing_pre_sweeps*/,
                                                          bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    int ierr;
    auto level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);
    const int U_scratch_idx = d_side_scratch_idx;
    const int P_scratch_idx = d_cell_scratch_idx;

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (auto p = level->begin(); p != level->end(); ++p, ++patch_counter)
        {
            auto patch = *p;

            auto U_error_data = BOOST_CAST<SideData<double> >(error.getComponentPatchData(0, *patch));
            auto U_scratch_data = BOOST_CAST<SideData<double> >(patch->getPatchData(U_scratch_idx));
            const Box& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box.isSpatiallyEqual(U_scratch_data->getGhostBox()));
            TBOX_ASSERT(U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                U_scratch_data->getArrayData(axis).copy(U_error_data->getArrayData(axis),
                                                        d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                        IntVector::getZero(DIM));
            }

            auto P_error_data = BOOST_CAST<CellData<double> >(error.getComponentPatchData(1, *patch));
            auto P_scratch_data = BOOST_CAST<CellData<double> >(patch->getPatchData(P_scratch_idx));
            const Box& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box.isSpatiallyEqual(P_scratch_data->getGhostBox()));
            TBOX_ASSERT(P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
            P_scratch_data->getArrayData().copy(P_error_data->getArrayData(),
                                                d_patch_cell_bc_box_overlap[level_num][patch_counter],
                                                IntVector::getZero(DIM));
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
                for (auto p = level->begin(); p != level->end(); ++p, ++patch_counter)
                {
                    auto patch = *p;

                    auto U_error_data = BOOST_CAST<SideData<double> >(error.getComponentPatchData(0, *patch));
                    auto U_scratch_data = BOOST_CAST<SideData<double> >(patch->getPatchData(U_scratch_idx));
                    const Box& U_ghost_box = U_error_data->getGhostBox();
                    TBOX_ASSERT(U_ghost_box.isSpatiallyEqual(U_scratch_data->getGhostBox()));
                    TBOX_ASSERT(U_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        U_error_data->getArrayData(axis)
                            .copy(U_scratch_data->getArrayData(axis),
                                  d_patch_side_bc_box_overlap[level_num][patch_counter][axis], IntVector::getZero(DIM));
                    }

                    auto P_error_data = BOOST_CAST<CellData<double> >(error.getComponentPatchData(1, *patch));
                    auto P_scratch_data = BOOST_CAST<CellData<double> >(patch->getPatchData(P_scratch_idx));
                    const Box& P_ghost_box = P_error_data->getGhostBox();
                    TBOX_ASSERT(P_ghost_box.isSpatiallyEqual(P_scratch_data->getGhostBox()));
                    TBOX_ASSERT(P_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
                    P_error_data->getArrayData().copy(P_scratch_data->getArrayData(),
                                                      d_patch_cell_bc_box_overlap[level_num][patch_counter],
                                                      IntVector::getZero(DIM));
                }

                // Fill the non-coarse-fine interface ghost cell values.
                const std::pair<int, int> error_idxs = std::make_pair(U_error_idx, P_error_idx);
                xeqScheduleGhostFillNoCoarse(error_idxs, level_num);
            }

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension.
            d_U_cf_bdry_op->setPatchDataIndex(U_error_idx);
            d_P_cf_bdry_op->setPatchDataIndex(P_error_idx);
            const IntVector& ratio = level->getRatioToCoarserLevel();
            for (auto p = level->begin(); p != level->end(); ++p)
            {
                auto patch = *p;
                const IntVector ghost_width_to_fill(DIM, d_gcw);
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
        Vec& err = d_box_e[level_num];
        Vec& res = d_box_r[level_num];
        KSP& ksp = d_box_ksp[level_num];
        int patch_counter = 0;
        for (auto p = level->begin(); p != level->end(); ++p, ++patch_counter)
        {
            auto patch = *p;
            auto U_error_data = BOOST_CAST<SideData<double> >(error.getComponentPatchData(0, *patch));
            auto U_residual_data = BOOST_CAST<SideData<double> >(residual.getComponentPatchData(0, *patch));
            const Box& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box.isSpatiallyEqual(U_residual_data->getGhostBox()));
            TBOX_ASSERT(U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_residual_data->getGhostCellWidth() == d_gcw);
            auto P_error_data = BOOST_CAST<CellData<double> >(error.getComponentPatchData(1, *patch));
            auto P_residual_data = BOOST_CAST<CellData<double> >(residual.getComponentPatchData(1, *patch));
            const Box& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box.isSpatiallyEqual(P_residual_data->getGhostBox()));
            TBOX_ASSERT(P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_residual_data->getGhostCellWidth() == d_gcw);

            // Smooth the error on the patch.
            const Box& patch_box = patch->getBox();
            const auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
            const double* const dx = pgeom->getDx();
            for (auto b = patch_box.begin(), e = patch_box.end(); b != e; ++b)
            {
                const Index& i = *b;
                const Box box(i, i, BlockId::invalidId());
                copyToVec(err, *U_error_data, *P_error_data, box, box);
                copyToVec(res, *U_residual_data, *P_residual_data, box, box);
                modifyRhsForBcs(res, *U_error_data, *P_error_data, d_U_problem_coefs, box, box, dx);
                ierr = KSPSolve(ksp, res, err);
                IBTK_CHKERRQ(ierr);
                copyFromVec(err, *U_error_data, *P_error_data, box, box);
            }
        }
    }

    // Synchronize data along patch boundaries.
    xeqScheduleDataSynch(U_error_idx, level_num);
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void StaggeredStokesBoxRelaxationFACOperator::initializeOperatorStateSpecialized(const SAMRAIVectorReal<double>&
                                                                                 /*solution*/,
                                                                                 const SAMRAIVectorReal<double>&
                                                                                 /*rhs*/,
                                                                                 const int coarsest_reset_ln,
                                                                                 const int finest_reset_ln)
{
    // Initialize the box relaxation data on each level of the patch hierarchy.
    d_box_op.resize(d_finest_ln + 1);
    d_box_e.resize(d_finest_ln + 1);
    d_box_r.resize(d_finest_ln + 1);
    d_box_ksp.resize(d_finest_ln + 1);
    const Box box(Index::getZeroIndex(DIM), Index::getZeroIndex(DIM), BlockId::invalidId());
    auto geometry = BOOST_CAST<CartesianGridGeometry>(d_hierarchy->getGridGeometry());
    const double* const dx_coarsest = geometry->getDx();
    boost::array<double, NDIM> dx;
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        const IntVector& ratio = d_hierarchy->getPatchLevel(ln)->getRatioToLevelZero();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx_coarsest[d] / static_cast<double>(ratio(d));
        }
        buildBoxOperator(d_box_op[ln], d_U_problem_coefs, box, box, dx);
        int ierr;
        ierr = MatGetVecs(d_box_op[ln], &d_box_e[ln], &d_box_r[ln]);
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
        auto level = d_hierarchy->getPatchLevel(ln);
        const auto num_local_patches = level->getProcessorMapping().getLocalIndices().size();
        d_patch_side_bc_box_overlap[ln].resize(num_local_patches);
        int patch_counter = 0;
        for (auto p = level->begin(); p != level->end(); ++p, ++patch_counter)
        {
            auto patch = *p;
            const Box& patch_box = patch->getBox();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box side_box = SideGeometry::toSideBox(patch_box, axis);
                const Box side_ghost_box = Box::grow(side_box, IntVector::getOne(DIM));
                d_patch_side_bc_box_overlap[ln][patch_counter][axis] = BoxContainer(side_ghost_box);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis].removeIntersections(side_box);
            }
        }
    }

    d_patch_cell_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);

        const auto num_local_patches = level->getProcessorMapping().getLocalIndices().size();
        d_patch_cell_bc_box_overlap[ln].resize(num_local_patches);

        int patch_counter = 0;
        for (auto p = level->begin(); p != level->end(); ++p, ++patch_counter)
        {
            auto patch = *p;
            const Box& patch_box = patch->getBox();
            const Box& ghost_box = Box::grow(patch_box, IntVector::getOne(DIM));

            d_patch_cell_bc_box_overlap[ln][patch_counter] = BoxContainer(ghost_box);
            d_patch_cell_bc_box_overlap[ln][patch_counter].removeIntersections(patch_box);
        }
    }
    return;
}

void StaggeredStokesBoxRelaxationFACOperator::deallocateOperatorStateSpecialized(const int coarsest_reset_ln,
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
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
