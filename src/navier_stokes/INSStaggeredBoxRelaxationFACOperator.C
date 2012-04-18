// Filename: INSStaggeredBoxRelaxationFACOperator.C
// Created on 11 Jun 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "INSStaggeredBoxRelaxationFACOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CartCellDoubleQuadraticCFInterpolation.h>
#include <ibtk/CartSideDoubleQuadraticCFInterpolation.h>
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/SideSynchCopyFillPattern.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int GHOSTS = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values; used only to evaluate composite grid
// residuals.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries; used only to evaluate
// composite grid residuals.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells; used only to evaluate composite grid residuals.
static const bool CONSISTENT_TYPE_2_BDRY = false;

inline int
compute_side_index(
    const Index<NDIM>& i,
    const Box<NDIM>& box,
    const unsigned int axis)
{
    const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(box,axis);
    if (!side_box.contains(i)) return -1;
    int offset = 0;
    for (unsigned int d = 0; d < axis; ++d)
    {
        offset += SideGeometry<NDIM>::toSideBox(box,d).size();
    }
    return offset + side_box.offset(i);
}// compute_side_index

inline int
compute_cell_index(
    const Index<NDIM>& i,
    const Box<NDIM>& box)
{
    if (!box.contains(i)) return -1;
    int offset = 0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        offset += SideGeometry<NDIM>::toSideBox(box,axis).size();
    }
    return box.offset(i) + offset;
}// compute_cell_index

void
buildBoxOperator(
    Mat& A,
    const INSProblemCoefs& problem_coefs,
    const double dt,
    const Box<NDIM>& box,
    const Box<NDIM>& ghost_box,
    const blitz::TinyVector<double,NDIM>& dx)
{
    int ierr;

    const double rho = problem_coefs.getRho();
    const double mu = problem_coefs.getMu();
    const double lambda = problem_coefs.getLambda();

    // Allocate a PETSc matrix for the box operator.
    blitz::TinyVector<Box<NDIM>,NDIM> side_boxes;
    blitz::TinyVector<BoxList<NDIM>,NDIM> side_ghost_boxes;
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

    static const int U_stencil_sz = 2*NDIM+3;
    static const int P_stencil_sz = 2*NDIM+1;
    std::vector<int> nnz(size, 0);

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (Box<NDIM>::Iterator b(side_boxes[axis]); b; b++)
        {
            nnz[compute_side_index(b(), ghost_box, axis)] = std::min(size,U_stencil_sz);
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
        nnz[compute_cell_index(b(), ghost_box)] = std::min(size,P_stencil_sz);
    }
    for (BoxList<NDIM>::Iterator bl(cell_ghost_boxes); bl; bl++)
    {
        for (Box<NDIM>::Iterator b(bl()); b; b++)
        {
            nnz[compute_cell_index(b(), ghost_box)] = 1;
        }
    }

    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, PETSC_DEFAULT, &nnz[0], &A);  IBTK_CHKERRQ(ierr);

    // Set some general matrix options.
#ifdef DEBUG_CHECK_ASSERTIONS
    ierr = MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);  IBTK_CHKERRQ(ierr);
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

            std::vector<int> mat_cols(U_stencil_sz,-1);
            std::vector<double> mat_vals(U_stencil_sz,0.0);

            mat_cols[0] = mat_row;
            mat_vals[0] = (rho/dt) + 0.5*lambda;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                Index<NDIM> shift = 0;
                shift(d) = 1;
                const Index<NDIM> u_left = i - shift;
                const Index<NDIM> u_rght = i + shift;
                mat_cols[2*d+1] = compute_side_index(u_left, ghost_box, axis);
                mat_cols[2*d+2] = compute_side_index(u_rght, ghost_box, axis);

                mat_vals[    0] +=      mu/(dx[d]*dx[d]);
                mat_vals[2*d+1]  = -0.5*mu/(dx[d]*dx[d]);
                mat_vals[2*d+2]  = -0.5*mu/(dx[d]*dx[d]);
            }

            Index<NDIM> shift = 0;
            shift(axis) = 1;
            const Index<NDIM> p_left = i - shift;
            const Index<NDIM> p_rght = i;
            mat_cols[2*NDIM+1] = compute_cell_index(p_left, ghost_box);
            mat_cols[2*NDIM+2] = compute_cell_index(p_rght, ghost_box);

            mat_vals[2*NDIM+1] = -1.0/dx[axis];
            mat_vals[2*NDIM+2] =  1.0/dx[axis];

            static const int m = 1;
            static const int n = U_stencil_sz;
            ierr = MatSetValues(A, m, &mat_row, n, &mat_cols[0], &mat_vals[0], INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        Index<NDIM> i = b();
        const int mat_row = compute_cell_index(i, ghost_box);

        std::vector<int> mat_cols(P_stencil_sz,-1);
        std::vector<double> mat_vals(P_stencil_sz,0.0);

        mat_cols[0] = mat_row;
        mat_vals[0] = 0.0;

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            Index<NDIM> shift = 0;
            shift(axis) = 1;
            const Index<NDIM> u_left = i;
            const Index<NDIM> u_rght = i + shift;
            mat_cols[2*axis+1] = compute_side_index(u_left, ghost_box, axis);
            mat_cols[2*axis+2] = compute_side_index(u_rght, ghost_box, axis);

            mat_vals[2*axis+1] =  1.0/dx[axis];
            mat_vals[2*axis+2] = -1.0/dx[axis];
        }

        static const int m = 1;
        static const int n = P_stencil_sz;
        ierr = MatSetValues(A, m, &mat_row, n, &mat_cols[0], &mat_vals[0], INSERT_VALUES);  IBTK_CHKERRQ(ierr);
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
                ierr = MatSetValue(A, i, i, 1.0, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    for (BoxList<NDIM>::Iterator bl(cell_ghost_boxes); bl; bl++)
    {
        for (Box<NDIM>::Iterator b(bl()); b; b++)
        {
            const int i = compute_cell_index(b(), ghost_box);
            ierr = MatSetValue(A, i, i, 1.0, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrices.
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    return;
}// buildBoxOperator

void
modifyRhsForBcs(
    Vec& v,
    const SideData<NDIM,double>& U_data,
    const CellData<NDIM,double>& P_data,
    const INSProblemCoefs& problem_coefs,
    const double /*dt*/,
    const Box<NDIM>& box,
    const Box<NDIM>& ghost_box,
    const double* const dx)
{
    int ierr;

    const double mu = problem_coefs.getMu();
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
                    ierr = VecSetValue(v, idx, +0.5*mu*U_data(SideIndex<NDIM>(u_left, axis, SideIndex<NDIM>::Lower))/(dx[d]*dx[d]), ADD_VALUES); IBTK_CHKERRQ(ierr);
                }
                if (!side_box.contains(u_rght))
                {
                    ierr = VecSetValue(v, idx, +0.5*mu*U_data(SideIndex<NDIM>(u_rght, axis, SideIndex<NDIM>::Lower))/(dx[d]*dx[d]), ADD_VALUES); IBTK_CHKERRQ(ierr);
                }
            }

            Index<NDIM> shift = 0;
            shift(axis) = 1;
            const Index<NDIM> p_left = i - shift;
            const Index<NDIM> p_rght = i;
            if (!box.contains(p_left))
            {
                ierr = VecSetValue(v, idx, +P_data(p_left)/dx[axis], ADD_VALUES); IBTK_CHKERRQ(ierr);
            }
            if (!box.contains(p_rght))
            {
                ierr = VecSetValue(v, idx, -P_data(p_rght)/dx[axis], ADD_VALUES); IBTK_CHKERRQ(ierr);
            }
        }
    }

    ierr = VecAssemblyBegin(v);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);  IBTK_CHKERRQ(ierr);
    return;
}// modifyRhsForBcs

inline void
copyToVec(
    Vec& v,
    const SideData<NDIM,double>& U_data,
    const CellData<NDIM,double>& P_data,
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
            ierr = VecSetValue(v, idx, U_data(s_i), INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        const Index<NDIM>& i = b();
        const int idx = compute_cell_index(i, ghost_box);
        ierr = VecSetValue(v, idx, P_data(i), INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(v);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);  IBTK_CHKERRQ(ierr);
    return;
}// copyToVec

inline void
copyFromVec(
    Vec& v,
    SideData<NDIM,double>& U_data,
    CellData<NDIM,double>& P_data,
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
            ierr = VecGetValues(v, 1, &idx, &U);  IBTK_CHKERRQ(ierr);
            U_data(s_i) = (1.0-omega)*U_data(s_i) + omega*U;
        }
    }

    double P;
    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        const Index<NDIM>& i = b();
        const int idx = compute_cell_index(i, ghost_box);
        ierr = VecGetValues(v, 1, &idx, &P);  IBTK_CHKERRQ(ierr);
        P_data(i) = (1.0-omega)*P_data(i) + omega*P;
    }
    return;
}// copyFromVec
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredBoxRelaxationFACOperator::INSStaggeredBoxRelaxationFACOperator(
    const std::string& object_name,
    const Pointer<Database> input_db)
    : INSStaggeredFACPreconditionerStrategy(object_name, GHOSTS, input_db),
      d_problem_coefs(),
      d_default_U_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_U_bc_coef", Pointer<Database>(NULL))),
      d_U_bc_coefs(blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>(d_default_U_bc_coef)),
      d_default_P_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_P_bc_coef", Pointer<Database>(NULL))),
      d_P_bc_coef(d_default_P_bc_coef),
      d_box_op(),
      d_box_e(),
      d_box_r(),
      d_box_ksp(),
      d_patch_side_bc_box_overlap(),
      d_patch_cell_bc_box_overlap()
{
    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions for the velocity and homogeneous Neumann
    // boundary conditions for the pressure.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_U_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_U_bc_coef->setBoundaryValue(2*d+1,0.0);
        d_default_P_bc_coef->setBoundarySlope(2*d  ,0.0);
        d_default_P_bc_coef->setBoundarySlope(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setPhysicalBcCoefs(blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>(d_default_U_bc_coef),d_default_P_bc_coef);
    return;
}// INSStaggeredBoxRelaxationFACOperator

INSStaggeredBoxRelaxationFACOperator::~INSStaggeredBoxRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    delete d_default_U_bc_coef;
    delete d_default_P_bc_coef;
    return;
}// ~INSStaggeredBoxRelaxationFACOperator

void
INSStaggeredBoxRelaxationFACOperator::setProblemCoefficients(
    const INSProblemCoefs& problem_coefs)
{
    d_problem_coefs = problem_coefs;
    return;
}// setProblemCoefficients

void
INSStaggeredBoxRelaxationFACOperator::setPhysicalBcCoefs(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& U_bc_coefs,
    RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (U_bc_coefs[d] != NULL)
        {
            d_U_bc_coefs[d] = U_bc_coefs[d];
        }
        else
        {
            d_U_bc_coefs[d] = d_default_U_bc_coef;
        }
    }

    if (P_bc_coef != NULL)
    {
        d_P_bc_coef = P_bc_coef;
    }
    else
    {
        d_P_bc_coef = d_default_P_bc_coef;
    }
    return;
}// setPhysicalBcCoefs

void
INSStaggeredBoxRelaxationFACOperator::smoothError(
    SAMRAIVectorReal<NDIM,double>& error,
    const SAMRAIVectorReal<NDIM,double>& residual,
    int level_num,
    int num_sweeps,
    bool /*performing_pre_sweeps*/,
    bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    int ierr;
    const double dt = d_new_time-d_current_time;
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

            Pointer<SideData<NDIM,double> >   U_error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM,double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
            TBOX_ASSERT(  U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                U_scratch_data->getArrayData(axis).copy(
                    U_error_data->getArrayData(axis),
                    d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                    IntVector<NDIM>(0));
            }

            Pointer<CellData<NDIM,double> >   P_error_data = error.getComponentPatchData(1, *patch);
            Pointer<CellData<NDIM,double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
            TBOX_ASSERT(  P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            P_scratch_data->getArrayData().copy(
                P_error_data->getArrayData(),
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

                    Pointer<SideData<NDIM,double> >   U_error_data = error.getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM,double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                    const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
                    TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
                    TBOX_ASSERT(  U_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        U_error_data->getArrayData(axis).copy(
                            U_scratch_data->getArrayData(axis),
                            d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                            IntVector<NDIM>(0));
                    }

                    Pointer<CellData<NDIM,double> >   P_error_data = error.getComponentPatchData(1, *patch);
                    Pointer<CellData<NDIM,double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                    const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
                    TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
                    TBOX_ASSERT(  P_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    P_error_data->getArrayData().copy(
                        P_scratch_data->getArrayData(),
                        d_patch_cell_bc_box_overlap[level_num][patch_counter],
                        IntVector<NDIM>(0));
                }

                // Fill the non-coarse-fine interface ghost cell values.
                const std::pair<int,int> error_idxs = std::make_pair(U_error_idx,P_error_idx);
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
            const std::pair<int,int> error_idxs = std::make_pair(U_error_idx, P_error_idx);
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
            Pointer<SideData<NDIM,double> >    U_error_data = error   .getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM,double> > U_residual_data = residual.getComponentPatchData(0, *patch);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_residual_data->getGhostBox());
            TBOX_ASSERT(   U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_residual_data->getGhostCellWidth() == d_gcw);
#endif
            Pointer<CellData<NDIM,double> >    P_error_data = error   .getComponentPatchData(1, *patch);
            Pointer<CellData<NDIM,double> > P_residual_data = residual.getComponentPatchData(1, *patch);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_residual_data->getGhostBox());
            TBOX_ASSERT(   P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_residual_data->getGhostCellWidth() == d_gcw);
#endif
            // Smooth the error on the patch.
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            for (Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const Index<NDIM>& i = b();
                const Box<NDIM> box(i,i);
                copyToVec(e, *U_error_data, *P_error_data, box, box);
                copyToVec(r, *U_residual_data, *P_residual_data, box, box);
                modifyRhsForBcs(r, *U_error_data, *P_error_data, d_problem_coefs, dt, box, box, dx);
                ierr = KSPSolve(ksp, r, e);  IBTK_CHKERRQ(ierr);
                copyFromVec(e, *U_error_data, *P_error_data, box, box);
            }
        }
    }

    // Synchronize data along patch boundaries.
    xeqScheduleSideDataSynch(U_error_idx, level_num);
    return;
}// smoothError

bool
INSStaggeredBoxRelaxationFACOperator::solveCoarsestLevel(
    SAMRAIVectorReal<NDIM,double>& error,
    const SAMRAIVectorReal<NDIM,double>& residual,
    int coarsest_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
#endif
    smoothError(error, residual, coarsest_ln, d_coarse_solver_max_its, false, false);
    return true;
}// solveCoarsestLevel

void
INSStaggeredBoxRelaxationFACOperator::computeResidual(
    SAMRAIVectorReal<NDIM,double>& residual,
    const SAMRAIVectorReal<NDIM,double>& solution,
    const SAMRAIVectorReal<NDIM,double>& rhs,
    int coarsest_level_num,
    int finest_level_num)
{
    const int U_res_idx = residual.getComponentDescriptorIndex(0);
    const int U_sol_idx = solution.getComponentDescriptorIndex(0);
    const int U_rhs_idx = rhs.getComponentDescriptorIndex(0);

    const Pointer<SideVariable<NDIM,double> > U_res_sc_var = residual.getComponentVariable(0);
    const Pointer<SideVariable<NDIM,double> > U_sol_sc_var = solution.getComponentVariable(0);
    const Pointer<SideVariable<NDIM,double> > U_rhs_sc_var = rhs.getComponentVariable(0);

    const int P_res_idx = residual.getComponentDescriptorIndex(1);
    const int P_sol_idx = solution.getComponentDescriptorIndex(1);
    const int P_rhs_idx = rhs.getComponentDescriptorIndex(1);

    const Pointer<CellVariable<NDIM,double> > P_res_cc_var = residual.getComponentVariable(1);
    const Pointer<CellVariable<NDIM,double> > P_sol_cc_var = solution.getComponentVariable(1);
    const Pointer<CellVariable<NDIM,double> > P_rhs_cc_var = rhs.getComponentVariable(1);

    // Fill ghost-cell values.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    Pointer<VariableFillPattern<NDIM> > sc_fill_pattern = new SideNoCornersFillPattern(GHOSTS, false, false, true);
    Pointer<VariableFillPattern<NDIM> > cc_fill_pattern = new CellNoCornersFillPattern(GHOSTS, false, false, true);
    InterpolationTransactionComponent U_scratch_component(U_sol_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs, sc_fill_pattern);
    InterpolationTransactionComponent P_scratch_component(P_sol_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef , cc_fill_pattern);
    std::vector<InterpolationTransactionComponent> U_P_components(2);
    U_P_components[0] = U_scratch_component;
    U_P_components[1] = P_scratch_component;
    if (d_hier_bdry_fill_ops[finest_level_num].isNull())
    {
        d_hier_bdry_fill_ops[finest_level_num] = new HierarchyGhostCellInterpolation();
        d_hier_bdry_fill_ops[finest_level_num]->initializeOperatorState(U_P_components, d_hierarchy, coarsest_level_num, finest_level_num);
    }
    else
    {
        d_hier_bdry_fill_ops[finest_level_num]->resetTransactionComponents(U_P_components);
    }
    d_hier_bdry_fill_ops[finest_level_num]->setHomogeneousBc(true);
    d_hier_bdry_fill_ops[finest_level_num]->fillData(d_new_time);

    // Compute the residual, r = f - A*u.
    if (d_hier_math_ops[finest_level_num].isNull())
    {
        std::ostringstream stream;
        stream << d_object_name << "::hier_math_ops_" << finest_level_num;
        d_hier_math_ops[finest_level_num] = new HierarchyMathOps(stream.str(), d_hierarchy, coarsest_level_num, finest_level_num);
    }
    d_hier_math_ops[finest_level_num]->grad(U_res_idx, U_res_sc_var, /*cf_bdry_synch*/ true, 1.0, P_sol_idx, P_sol_cc_var, NULL, d_new_time);
    const double rho = d_problem_coefs.getRho();
    const double mu = d_problem_coefs.getMu();
    const double lambda = d_problem_coefs.getLambda();
    const double dt = d_new_time-d_current_time;
    PoissonSpecifications helmholtz_spec("");
    helmholtz_spec.setCConstant((rho/dt)+0.5*lambda);
    helmholtz_spec.setDConstant(        -0.5*mu    );
    d_hier_math_ops[finest_level_num]->laplace(U_res_idx, U_res_sc_var, helmholtz_spec, U_sol_idx, U_sol_sc_var, NULL, d_new_time, 1.0, U_res_idx, U_res_sc_var);
    HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    hier_sc_data_ops.axpy(U_res_idx, -1.0, U_res_idx, U_rhs_idx, false);
    d_hier_math_ops[finest_level_num]->div(P_res_idx, P_res_cc_var, -1.0, U_sol_idx, U_sol_sc_var, NULL, d_new_time, /*cf_bdry_synch*/ true);
    HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    hier_cc_data_ops.axpy(P_res_idx, -1.0, P_res_idx, P_rhs_idx, false);
    return;
}// computeResidual

/////////////////////////////// PROTECTED ////////////////////////////////////

void
INSStaggeredBoxRelaxationFACOperator::initializeOperatorStateSpecialized(
    const SAMRAIVectorReal<NDIM,double>& /*solution*/,
    const SAMRAIVectorReal<NDIM,double>& /*rhs*/,
    const int coarsest_reset_ln,
    const int finest_reset_ln)
{
    // Setup boundary condition handling objects.
    d_U_bc_op = new CartSideRobinPhysBdryOp(d_side_scratch_idx, d_U_bc_coefs, false);
    d_P_bc_op = new CartCellRobinPhysBdryOp(d_cell_scratch_idx, d_P_bc_coef , false);
    d_U_cf_bdry_op = new CartSideDoubleQuadraticCFInterpolation();
    d_P_cf_bdry_op = new CartCellDoubleQuadraticCFInterpolation();
    d_U_op_stencil_fill_pattern = new SideNoCornersFillPattern(GHOSTS, false, false, false);
    d_P_op_stencil_fill_pattern = new CellNoCornersFillPattern(GHOSTS, false, false, false);
    d_U_synch_fill_pattern = new SideSynchCopyFillPattern();

    // Initialize the box relaxation data on each level of the patch hierarchy.
    d_box_op.resize(d_finest_ln+1);
    d_box_e.resize(d_finest_ln+1);
    d_box_r.resize(d_finest_ln+1);
    d_box_ksp.resize(d_finest_ln+1);
    const Box<NDIM> box(Index<NDIM>(0),Index<NDIM>(0));
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    const double* const dx_coarsest = geometry->getDx();
    blitz::TinyVector<double,NDIM> dx;
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        const IntVector<NDIM>& ratio = d_hierarchy->getPatchLevel(ln)->getRatio();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx_coarsest[d]/static_cast<double>(ratio(d));
        }
        const double dt = d_new_time-d_current_time;
        buildBoxOperator(d_box_op[ln], d_problem_coefs, dt, box, box, dx);
        int ierr;
        ierr = MatGetVecs(d_box_op[ln], &d_box_e[ln], &d_box_r[ln]);  IBTK_CHKERRQ(ierr);
        ierr = KSPCreate(PETSC_COMM_SELF, &d_box_ksp[ln]);  IBTK_CHKERRQ(ierr);
        ierr = KSPSetOperators(d_box_ksp[ln], d_box_op[ln], d_box_op[ln], SAME_PRECONDITIONER);  IBTK_CHKERRQ(ierr);
        ierr = KSPSetType(d_box_ksp[ln], KSPPREONLY);  IBTK_CHKERRQ(ierr);
        PC box_pc;
        ierr = KSPGetPC(d_box_ksp[ln], &box_pc);  IBTK_CHKERRQ(ierr);
        ierr = PCSetType(box_pc, PCLU);  IBTK_CHKERRQ(ierr);
        ierr = PCFactorReorderForNonzeroDiagonal(box_pc, std::numeric_limits<double>::epsilon());  IBTK_CHKERRQ(ierr);
        ierr = KSPSetUp(d_box_ksp[ln]);  IBTK_CHKERRQ(ierr);
    }

    // Get overlap information for setting patch boundary conditions.
    d_patch_side_bc_box_overlap.resize(d_finest_ln+1);
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
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box,axis);
                const Box<NDIM> side_ghost_box = Box<NDIM>::grow(side_box, 1);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis] = BoxList<NDIM>(side_ghost_box);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis].removeIntersections(side_box);
            }
        }
    }

    d_patch_cell_bc_box_overlap.resize(d_finest_ln+1);
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
}// initializeOperatorStateSpecialized

void
INSStaggeredBoxRelaxationFACOperator::deallocateOperatorStateSpecialized(
    const int coarsest_reset_ln,
    const int finest_reset_ln)
{
    if (!d_is_initialized) return;
    for (int ln = coarsest_reset_ln; ln <= std::min(d_finest_ln,finest_reset_ln); ++ln)
    {
        int ierr;
        ierr = MatDestroy(&d_box_op [ln]);  IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&d_box_e  [ln]);  IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&d_box_r  [ln]);  IBTK_CHKERRQ(ierr);
        ierr = KSPDestroy(&d_box_ksp[ln]);  IBTK_CHKERRQ(ierr);
        d_patch_side_bc_box_overlap[ln].resize(0);
        d_patch_cell_bc_box_overlap[ln].resize(0);
    }
    return;
}// deallocateOperatorStateSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
