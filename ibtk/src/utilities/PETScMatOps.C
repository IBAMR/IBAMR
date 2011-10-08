// Filename: PETScMatOps.C
// Created on 27 Sep 2010 by Boyce Griffith
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

#include "PETScMatOps.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#if PETSC_CLANGUAGE_C
extern "C"
#endif
{
#include "../src/mat/impls/aij/mpi/mpiaij.h"
#include "petscblaslapack.h"
}

// C++ STDLIB INCLUDES
#include <cstring>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

#undef __FUNCT__
#define __FUNCT__ "MatAXPY"
PetscErrorCode
PETScMatOps::MatAXPY(
    Mat Y,
    PetscScalar a,
    Mat X)
{
    PetscErrorCode ierr;
    const MatType X_type;
    const MatType Y_type;

    PetscFunctionBegin;

    ierr = MatGetType(X,&X_type);CHKERRQ(ierr);
    ierr = MatGetType(Y,&Y_type);CHKERRQ(ierr);
    if      ((strcmp(X_type,MATSEQAIJ) == 0) && (strcmp(Y_type,MATSEQAIJ) == 0))
    {
        ierr = MatAXPY_SeqAIJ(Y,a,X);CHKERRQ(ierr);
    }
    else if ((strcmp(X_type,MATMPIAIJ) == 0) && (strcmp(Y_type,MATMPIAIJ) == 0))
    {
        ierr = MatAXPY_MPIAIJ(Y,a,X);CHKERRQ(ierr);
    }
    else
    {
        SETERRQ2(PETSC_ERR_ARG_WRONG, "unsupported matrix types; X type = %s, Y type = %s", X_type, Y_type);
    }

    PetscFunctionReturn(0);
}// MatAXPY

#undef __FUNCT__
#define __FUNCT__ "MatAXPY_MPIAIJ"
PetscErrorCode
PETScMatOps::MatAXPY_MPIAIJ(
    Mat Y,
    PetscScalar a,
    Mat X)
{
    PetscErrorCode ierr;
    Mat_MPIAIJ *X_aij = (Mat_MPIAIJ*)X->data;
    Mat_MPIAIJ *Y_aij = (Mat_MPIAIJ*)Y->data;

    PetscFunctionBegin;

    /*
     * Perform the AXPY operation on the diagonal and off-diagonal submatrices.
     */
    ierr = MatAXPY_SeqAIJ(Y_aij->A,a,X_aij->A);CHKERRQ(ierr);
    ierr = MatAXPY_SeqAIJ(Y_aij->B,a,X_aij->B);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}// MatAXPY_MPIAIJ

#undef __FUNCT__
#define __FUNCT__ "MatAXPY_SeqAIJ"
PetscErrorCode
PETScMatOps::MatAXPY_SeqAIJ(
    Mat Y,
    PetscScalar a,
    Mat X)
{
    PetscErrorCode ierr;
    Mat W;
    Mat_SeqAIJ *X_aij = (Mat_SeqAIJ*)X->data;
    Mat_SeqAIJ *Y_aij = (Mat_SeqAIJ*)Y->data;
    PetscInt m, m_X, m_Y, n, n_X, n_Y;
    PetscInt *nnz, row;
    PetscInt *X_ai = X_aij->i, *X_aj = X_aij->j, *X_ailen = X_aij->ilen;
    PetscInt *Y_ai = Y_aij->i, *Y_aj = Y_aij->j, *Y_ailen = Y_aij->ilen;
    PetscInt *j_X, *j_Y, *j_X_next, *j_Y_next;
    PetscInt X_len, X_len_max;
    MPI_Comm comm, comm_X, comm_Y;
    PetscScalar *v;
    PetscBLASInt one = 1, bn;

    PetscFunctionBegin;

    /*
     * Ensure matrices are row-oriented.
     */
    if (!Y_aij->roworiented || !X_aij->roworiented)
    {
        SETERRQ(PETSC_ERR_ARG_WRONG, "X and Y both must be row-oriented matrices");
    }

    /*
     * Ensure matrices have conforming sizes.
     */
    ierr = MatGetSize(X, &m_X, &n_X);CHKERRQ(ierr);
    ierr = MatGetSize(Y, &m_Y, &n_Y);CHKERRQ(ierr);
    if (m_X != m_Y)
    {
        SETERRQ(PETSC_ERR_ARG_SIZ, "X and Y must have the same number of rows");
    }
    if (n_X != n_Y)
    {
        SETERRQ(PETSC_ERR_ARG_SIZ, "X and Y must have the same number of cols");
    }
    m = m_X;
    n = n_X;

    /*
     * Ensure matrices share the same MPI communicator.
     */
    ierr = PetscObjectGetComm((PetscObject)X,&comm_X);CHKERRQ(ierr);
    ierr = PetscObjectGetComm((PetscObject)Y,&comm_Y);CHKERRQ(ierr);
    if (comm_X != comm_Y)
    {
        SETERRQ(PETSC_ERR_ARG_WRONG, "X and Y both must share the same MPI communicator");
    }
    comm = comm_X;

    /*
     * Determine non-zero structure of a*X + Y.
     *
     * NOTE: Assumes that column indices are sorted in increasing order.
     */
    ierr = PetscMalloc(m*sizeof(PetscInt), &nnz);CHKERRQ(ierr);
    X_len_max = 0;
    for (row = 0; row < m; ++row)
    {
        nnz[row] = 0;

        j_X = X_aj + X_ai[row];
        j_Y = Y_aj + Y_ai[row];

        j_X_next = X_aj + X_ai[row+1];
        j_Y_next = Y_aj + Y_ai[row+1];

        X_len = X_ai[row+1] - X_ai[row];
        if (X_len > X_len_max) X_len_max = X_len;

        while (j_X != j_X_next || j_Y != j_Y_next)
        {
            nnz[row] += 1;
            if (j_X != j_X_next && j_Y != j_Y_next)
            {
                if      (*j_X < *j_Y)  /* X col index is lower than Y col index */
                {
                    ++j_X;
                }
                else if (*j_X > *j_Y)  /* Y col index is lower than X col index */
                {
                    ++j_Y;
                }
                else                   /* X and Y col indices are equal         */
                {
                    ++j_X;
                    ++j_Y;
                }
            }
            else if (j_X != j_X_next)
            {
                ++j_X;
            }
            else if (j_Y != j_Y_next)
            {
                ++j_Y;
            }
        }
    }

    /*
     * Create the new matrix.
     */
    ierr = MatCreateSeqAIJ(comm,m,n,PETSC_DEFAULT,nnz,&W);CHKERRQ(ierr);
    ierr = PetscFree(nnz);CHKERRQ(ierr);

#ifdef DEBUG_CHECK_ASSERTIONS
    ierr = MatSetOption(W, MAT_NEW_NONZERO_LOCATION_ERR  , PETSC_TRUE); IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(W, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); IBTK_CHKERRQ(ierr);
#endif

    /*
     * Accumulate values in the new matrix and assemble it.
     */
    if (a != 0.0 && a != 1.0)
    {
        ierr = PetscMalloc(X_len_max*sizeof(PetscScalar), &v);CHKERRQ(ierr);
    }
    for (row = 0; row < m; ++row)
    {
        if (a == 1.0)
        {
            ierr = MatSetValues(W,
                                1, &row,
                                X_ailen[row], X_aj + X_ai[row],
                                X_aij->a + X_ai[row],
                                ADD_VALUES);CHKERRQ(ierr);
        }
        else if (a != 0.0)
        {
            ierr = PetscMemcpy(X_aij->a + X_ai[row], v, X_ailen[row]*sizeof(PetscScalar));CHKERRQ(ierr);
            bn = PetscBLASIntCast(X_ailen[row]);
            BLASscal_(&bn,&a,v,&one);
            ierr = MatSetValues(W,
                                1, &row,
                                X_ailen[row], X_aj + X_ai[row],
                                v,
                                ADD_VALUES);CHKERRQ(ierr);
        }
        ierr = MatSetValues(W,
                            1, &row,
                            Y_ailen[row], Y_aj + Y_ai[row],
                            Y_aij->a + Y_ai[row],
                            ADD_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(W,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    if (a != 0.0 && a != 1.0)
    {
        ierr = PetscFree(v);CHKERRQ(ierr);
    }
    ierr = MatAssemblyEnd(W,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatHeaderReplace(Y,W);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// MatAXPY_SeqAIJ

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
