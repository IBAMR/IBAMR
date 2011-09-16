// Filename: CartSideDoubleDivPreservingRefine.C
// Created on 09 Nov 2008 by Boyce Griffith
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

#include "CartSideDoubleDivPreservingRefine.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <SideData.h>
#include <SideVariable.h>
#include <tbox/MathUtilities.h>

// BLITZ++ INCLUDES
#include <blitz/tinyvec.h>

// C++ STDLIB INCLUDES
#include <limits>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
void
buildBoxOperator(
    Mat& A,
    const Box<NDIM>& box,
    const double* const dx)
{
    int ierr;

    // Allocate a PETSc matrix for the box operator.
    const int size = box.size();
    std::vector<int> nnz(size);
    for (int j = 0; j < size; ++j)
    {
        if (j == 0) nnz[j] = size;
        else        nnz[j] = std::min(size,2*NDIM+1);
    }
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, PETSC_DEFAULT, &nnz[0], &A);  IBTK_CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);

    // Set the matrix coefficients to correspond to the standard finite
    // difference approximation to the Laplacian with homogeneous Neumann
    // boundary conditions.  We explicitly force the mean value of the solution
    // to equal zero.
    for (int j = 0; j < size; ++j)
    {
        ierr = MatSetValue(A, 0, j, 1.0, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }
    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        const Index<NDIM>& i = b();
        const int mat_row = box.offset(i);
        if (mat_row == 0) continue;  // don't reset values in the first row
        std::vector<double> mat_vals(2*NDIM+1,0.0);
        std::vector<int   > mat_cols(2*NDIM+1, -1);
        mat_vals[0] = 0.0;
        mat_cols[0] = mat_row;
        for (int axis = 0, s = 1; axis < NDIM; ++axis)
        {
            for (int shift = -1; shift <= 1; shift += 2, ++s)
            {
                Index<NDIM> i_shift = i;
                i_shift(axis) += shift;
                if (!box.contains(i_shift)) continue;
                const int mat_col = box.offset(i_shift);
                mat_vals[0] -= 1.0/(dx[axis]*dx[axis]);  //     diagonal entry
                mat_vals[s] += 1.0/(dx[axis]*dx[axis]);  // off-diagonal entry
                mat_cols[s]  = mat_col;
            }
        }
        static const int m = 1;
        static const int n = mat_vals.size();
        ierr = MatSetValues(A, m, &mat_row, n, &mat_cols[0], &mat_vals[0], INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);

    Vec M_sum;
    VecCreateSeq(PETSC_COMM_SELF, size, &M_sum);
    MatGetRowSum(A,M_sum);
    for (int k = 0; k < size; ++k)
    {
        double sum;
        VecGetValues(M_sum, 1, &k, &sum);
    }
    return;
}// buildBoxOperator

inline void
formRHS(
    Vec& f_vec,
    const SideData<NDIM,double>& U_data,
    const Box<NDIM>& box,
    const double* const dx)
{
    int ierr;

    // u* = u + Grad Phi ===> Div Grad Phi = f = Div u* - Div u.
    //
    // Because we wish to maintain the discrete divergence of the overlying
    // coarse grid cell, we impose Div u = (Div u)_coarse instead of Div = 0.
    double div_u_coarse = 0.0;
    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        const Index<NDIM>& i = b();
        const int idx = box.offset(i);
        double div_u_star = 0.0;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const SideIndex<NDIM> s_i_upper(i, axis, SideIndex<NDIM>::Upper);
            const SideIndex<NDIM> s_i_lower(i, axis, SideIndex<NDIM>::Lower);
            div_u_star += (U_data(s_i_upper)-U_data(s_i_lower))/dx[axis];
        }
        div_u_coarse += div_u_star;
        ierr = VecSetValue(f_vec, idx, div_u_star, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }
    div_u_coarse /= static_cast<double>(box.size());
    ierr = VecAssemblyBegin(f_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd  (f_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecShift(f_vec, -div_u_coarse);  IBTK_CHKERRQ(ierr);
    ierr = VecSetValue(f_vec, 0, 0.0, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(f_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd  (f_vec);  IBTK_CHKERRQ(ierr);
    return;
}// formRHS

inline void
correctVelocity(
    Vec& phi_vec,
    SideData<NDIM,double>& U_data,
    const Box<NDIM>& box,
    const double* const dx)
{
    int ierr;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Box<NDIM> side_box = box;
        side_box.lower()(axis) += 1;
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const Index<NDIM>& i_upper = b();
            Index<NDIM>        i_lower = i_upper;
            i_lower(axis) -= 1;
            int idxs[2] = { box.offset(i_lower) , box.offset(i_upper) };
            double phi_vals[2];
            ierr = VecGetValues(phi_vec, 2, idxs, phi_vals);  IBTK_CHKERRQ(ierr);
            const SideIndex<NDIM> s_i(i_upper, axis, SideIndex<NDIM>::Lower);
            U_data(s_i) -= (phi_vals[1]-phi_vals[0])/dx[axis];
        }
    }
    return;
}// correctVelocity
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideDoubleDivPreservingRefine::CartSideDoubleDivPreservingRefine(
    const int u_dst_idx,
    const int u_src_idx,
    const int indicator_idx)
    : d_u_dst_idx(u_dst_idx),
      d_u_src_idx(u_src_idx),
      d_indicator_idx(indicator_idx)
{
    // intentionally blank
    return;
}// CartSideDoubleDivPreservingRefine

CartSideDoubleDivPreservingRefine::~CartSideDoubleDivPreservingRefine()
{
    // intentionally blank
    return;
}// ~CartSideDoubleDivPreservingRefine

void
CartSideDoubleDivPreservingRefine::setPhysicalBoundaryConditions(
    Patch<NDIM>& /*patch*/,
    const double /*fill_time*/,
    const IntVector<NDIM>& /*ghost_width_to_fill*/)
{
    // intentionally blank
    return;
}// setPhysicalBoundaryConditions

IntVector<NDIM>
CartSideDoubleDivPreservingRefine::getRefineOpStencilWidth() const
{
    return 0;
}// getRefineOpStencilWidth

void
CartSideDoubleDivPreservingRefine::preprocessRefine(
    Patch<NDIM>& /*fine*/,
    const Patch<NDIM>& /*coarse*/,
    const Box<NDIM>& /*fine_box*/,
    const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
    return;
}// preprocessRefine

void
CartSideDoubleDivPreservingRefine::postprocessRefine(
    Patch<NDIM>& fine,
    const Patch<NDIM>& /*coarse*/,
    const Box<NDIM>& fine_box,
    const IntVector<NDIM>& ratio)
{
    Pointer<SideData<NDIM,double> > fdata = fine.getPatchData(d_u_dst_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!fdata.isNull());
#endif
    const int fdata_depth = fdata->getDepth();

    // Reset the values of any fine grid values for which the indicator data is
    // set to "1".
    if (fine.checkAllocated(d_u_src_idx) && fine.checkAllocated(d_indicator_idx))
    {
        Pointer<SideData<NDIM,double> >     u_src_data = fine.getPatchData(    d_u_src_idx);
        Pointer<SideData<NDIM,double> > indicator_data = fine.getPatchData(d_indicator_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(fine_box,axis)); b; b++)
            {
                const Index<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i,axis,0);
                if (std::abs((*indicator_data)(i_s)-1.0) < 1.0e-12)
                {
                    for (int depth = 0; depth < fdata_depth; ++depth)
                    {
                        (*fdata)(i_s,depth) = (*u_src_data)(i_s,depth);
                    }
                }
            }
        }
    }

    // Determine the box on which we need to compute the divergence-preserving
    // correction.
    const Box<NDIM> correction_box = Box<NDIM>::refine(Box<NDIM>::coarsen(fine_box,ratio),ratio);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(fdata->getGhostBox().contains(correction_box));
#endif

    // Setup a linear solver to compute the local projection.
    int ierr;

    Box<NDIM> box;
    Index<NDIM>& box_lower = box.lower();
    Index<NDIM>& box_upper = box.upper();
    box_lower = Index<NDIM>(0);
    box_upper = ratio-Index<NDIM>(1);

    Mat L_mat;
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = fine.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    buildBoxOperator(L_mat, box, dx);

    Vec phi_vec, f_vec;
    ierr = MatGetVecs(L_mat, &phi_vec, &f_vec);  IBTK_CHKERRQ(ierr);

    KSP ksp;
    ierr = KSPCreate(PETSC_COMM_SELF, &ksp);  IBTK_CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, L_mat, L_mat, SAME_PRECONDITIONER);  IBTK_CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPPREONLY);  IBTK_CHKERRQ(ierr);

    PC pc;
    ierr = KSPGetPC(ksp, &pc);  IBTK_CHKERRQ(ierr);
    ierr = PCSetType(pc, PCLU);  IBTK_CHKERRQ(ierr);

    // Perform local projections to impose a discrete divergence condition.
#if (NDIM == 3)
    for (int n2 = fine_box.lower()(2); n2 < fine_box.upper()(2); n2 += ratio(2))
    {
        box_lower(2) = n2;
        box_upper(2) = n2+ratio(2)-1;
#endif
        for (int n1 = fine_box.lower()(1); n1 < fine_box.upper()(1); n1 += ratio(1))
        {
            box_lower(1) = n1;
            box_upper(1) = n1+ratio(1)-1;
            for (int n0 = fine_box.lower()(0); n0 < fine_box.upper()(0); n0 += ratio(0))
            {
                box_lower(0) = n0;
                box_upper(0) = n0+ratio(0)-1;

                formRHS(f_vec, *fdata, box, dx);
                ierr = KSPSolve(ksp, f_vec, phi_vec);  IBTK_CHKERRQ(ierr);
                correctVelocity(phi_vec, *fdata, box, dx);
            }
        }
#if (NDIM == 3)
    }
#endif

    // Clean-up allocated data.
    ierr = VecDestroy(phi_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(f_vec);  IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(L_mat);  IBTK_CHKERRQ(ierr);
    ierr = KSPDestroy(ksp);  IBTK_CHKERRQ(ierr);
    return;
}// postprocessRefine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
