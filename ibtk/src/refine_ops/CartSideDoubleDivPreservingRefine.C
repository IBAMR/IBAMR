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

// PETSC INCLUDES
#include <petsc.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <SideData.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline int
compute_cell_index(
    const Index<NDIM>& i,
    const Box<NDIM>& box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(box.contains(i));
#endif
    return box.offset(i);
}// compute_cell_index

void
buildBoxOperator(
    Mat& A,
    const Box<NDIM>& box,
    const double* const dx)
{
    int ierr;

    // Allocate a PETSc matrix for the box operator.
    const int size = box.size();
    ierr = MatCreateSeqDense(PETSC_COMM_SELF, size, size, PETSC_NULL, &A);  IBTK_CHKERRQ(ierr);

    // Set the matrix coefficients to correspond to the standard finite
    // difference approximation to the Laplacian with homogeneous Neumann
    // boundary conditions.  We clamp the (0,0) value to equal zero in order to
    // make the linear system nonsingular.
    ierr = MatSetValue(A, 0, 0, 1.0, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    const Index<NDIM>& box_lower = box.lower();
    const Index<NDIM>& box_upper = box.upper();
    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        Index<NDIM> i = b();
        Index<NDIM> i_shift;
        const int mat_row = compute_cell_index(i,box);
        if (mat_row == 0) continue;
        std::vector<int> mat_cols;
        std::vector<double> mat_vals;
        mat_cols.reserve(2*NDIM+1);
        mat_vals.reserve(2*NDIM+1);
        mat_cols.push_back(mat_row);
        mat_vals.push_back(0.0);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            if (i(axis) > box_lower(axis))
            {
                i_shift = i;
                i_shift(axis) -= 1;
                const int i_off_diagonal = compute_cell_index(i_shift,box);
                if (i_off_diagonal == 0) continue;
                mat_cols.push_back(i_off_diagonal);
                mat_vals[0] +=     -1.0/(dx[axis]*dx[axis]) ;
                mat_vals.push_back(+1.0/(dx[axis]*dx[axis]));
            }
            if (i(axis) < box_upper(axis))
            {
                i_shift = i;
                i_shift(axis) += 1;
                const int i_off_diagonal = compute_cell_index(i_shift,box);
                if (i_off_diagonal == 0) continue;
                mat_cols.push_back(i_off_diagonal);
                mat_vals[0] +=     -1.0/(dx[axis]*dx[axis]) ;
                mat_vals.push_back(+1.0/(dx[axis]*dx[axis]));
            }
        }
        static const int m = 1;
        static const int n = mat_vals.size();
        ierr = MatSetValues(A, m, &mat_row, n, &mat_cols[0], &mat_vals[0], INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    return;
}// buildBoxOperator

inline void
formRHS(
    Vec& div_u_star_vec,
    const SideData<NDIM,double>& U_data,
    const Box<NDIM>& box,
    const double* const dx)
{
    int ierr;
    double div_u_coarse = 0.0;

    // Div Grad Phi = f = Div u* - (Div u)_coarse.
    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        const Index<NDIM>& i = b();
        const int idx = compute_cell_index(i, box);
        double div_u_star = 0.0;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const SideIndex<NDIM> s_i_upper(i, axis, SideIndex<NDIM>::Upper);
            const SideIndex<NDIM> s_i_lower(i, axis, SideIndex<NDIM>::Lower);
            div_u_star += (U_data(s_i_upper)-U_data(s_i_lower))/dx[axis];
        }
        div_u_coarse += div_u_star;
        ierr = VecSetValue(div_u_star_vec, idx, div_u_star, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(div_u_star_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(div_u_star_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecShift(div_u_star_vec, -div_u_coarse);  IBTK_CHKERRQ(ierr);
    ierr = VecSetValue(div_u_star_vec, 0, 0.0, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(div_u_star_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(div_u_star_vec);  IBTK_CHKERRQ(ierr);
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
        pout << "axis = " << axis << "\n";
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const Index<NDIM>& i_upper = b();
            Index<NDIM>        i_lower = i_upper;
            i_lower(axis) -= 1;
            int idxs[2] = { compute_cell_index(i_lower,box) , compute_cell_index(i_upper,box) };
            double phi_vals[2];
            pout << "correcting index " << i_upper << "\n";
            ierr = VecGetValues(phi_vec, 2, idxs, phi_vals);  IBTK_CHKERRQ(ierr);
            const SideIndex<NDIM> s_i(b(), axis, SideIndex<NDIM>::Lower);
            U_data(s_i) -= (phi_vals[1]-phi_vals[0])/dx[axis];
        }
    }
    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        const Index<NDIM>& i = b();
        double div_u_star = 0.0;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const SideIndex<NDIM> s_i_upper(i, axis, SideIndex<NDIM>::Upper);
            const SideIndex<NDIM> s_i_lower(i, axis, SideIndex<NDIM>::Lower);
            div_u_star += (U_data(s_i_upper)-U_data(s_i_lower))/dx[axis];
        }
        pout << "i = " << i << " div u = " << div_u_star << "\n";
    }
    return;
}// correctVelocity
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideDoubleDivPreservingRefine::CartSideDoubleDivPreservingRefine(
    const int u_idx)
    : d_u_idx(u_idx)
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
    const Box<NDIM>& unrestricted_fine_box,
    const IntVector<NDIM>& ratio)
{
    int ierr;

    // NOTE: We can only perform the corrections on groups of interior fine
    // boxes that are completely covered by overlying coarse grid boxes.
    Pointer<SideData<NDIM,double> > U_fine_data = fine.getPatchData(d_u_idx);
    const Box<NDIM>& ghost_box = U_fine_data->getGhostBox();
    Box<NDIM> fine_box = Box<NDIM>::refine(Box<NDIM>::coarsen(unrestricted_fine_box,ratio),ratio);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        while (fine_box.lower()(axis) < ghost_box.lower()(axis)) fine_box.lower(axis) += ratio(axis);
        while (fine_box.upper()(axis) > ghost_box.upper()(axis)) fine_box.upper(axis) -= ratio(axis);
    }
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = fine.getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Setup a linear solver to compute the local projection.
    Box<NDIM> box;
    Index<NDIM>& box_lower = box.lower();
    Index<NDIM>& box_upper = box.upper();
    box_lower = Index<NDIM>(0);
    box_upper = ratio-Index<NDIM>(1);

    Mat A;
    buildBoxOperator(A, box, dx);

    Vec phi_vec, div_u_vec;
    ierr = MatGetVecs(A, &phi_vec, &div_u_vec);  IBTK_CHKERRQ(ierr);

    KSP ksp;
    ierr = KSPCreate(PETSC_COMM_SELF, &ksp);  IBTK_CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A, SAME_PRECONDITIONER);  IBTK_CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPPREONLY);  IBTK_CHKERRQ(ierr);

    PC pc;
    ierr = KSPGetPC(ksp, &pc);  IBTK_CHKERRQ(ierr);
    ierr = PCSetType(pc, PCLU);  IBTK_CHKERRQ(ierr);

    // Perform local projections.
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

                // Impose the discrete divergence-free condition.
                formRHS(div_u_vec, *U_fine_data, box, dx);
                ierr = KSPSolve(ksp, div_u_vec, phi_vec);  IBTK_CHKERRQ(ierr);
                correctVelocity(phi_vec, *U_fine_data, box, dx);
            }
        }
#if (NDIM == 3)
    }
#endif
    return;
}// postprocessRefine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
