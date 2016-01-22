// Filename: PETScMatLOWrapper.cpp
// Created on 26 Dec 2003 by Boyce Griffith
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

#include <stddef.h>
#include <string>

#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/PETScMatLOWrapper.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "mpi.h"
#include "petscmat.h"
#include "petscsys.h"
#include "tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScMatLOWrapper::PETScMatLOWrapper(const std::string& object_name, const Mat& petsc_mat)
    : LinearOperator(object_name),
      d_petsc_mat(petsc_mat),
      d_x(NULL),
      d_y(NULL),
      d_z(NULL),
      d_petsc_x(NULL),
      d_petsc_y(NULL),
      d_petsc_z(NULL)
{
    // intentionally blank
    return;
} // PETScMatLOWrapper()

PETScMatLOWrapper::~PETScMatLOWrapper()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~PETScMatLOWrapper()

const Mat&
PETScMatLOWrapper::getPETScMat() const
{
    return d_petsc_mat;
} // getPETScMat

void
PETScMatLOWrapper::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    if (!d_is_initialized) initializeOperatorState(x, y);

    // Update the PETSc Vec wrappers.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_y, Pointer<SAMRAIVectorReal<NDIM, double> >(&y, false));

    // Apply the operator.
    int ierr = MatMult(d_petsc_mat, d_petsc_x, d_petsc_y);
    IBTK_CHKERRQ(ierr);
    return;
} // apply

void
PETScMatLOWrapper::applyAdd(SAMRAIVectorReal<NDIM, double>& x,
                            SAMRAIVectorReal<NDIM, double>& y,
                            SAMRAIVectorReal<NDIM, double>& z)
{
    if (!d_is_initialized) initializeOperatorState(x, y);

    // Update the PETSc Vec wrappers.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_y, Pointer<SAMRAIVectorReal<NDIM, double> >(&y, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_z, Pointer<SAMRAIVectorReal<NDIM, double> >(&z, false));

    // Apply the operator.
    int ierr = MatMultAdd(d_petsc_mat, d_petsc_x, d_petsc_y, d_petsc_z);
    IBTK_CHKERRQ(ierr);
    return;
} // applyAdd

void
PETScMatLOWrapper::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                           const SAMRAIVectorReal<NDIM, double>& out)
{
    if (d_is_initialized) deallocateOperatorState();
    d_x = in.cloneVector("");
    d_y = out.cloneVector("");
    d_z = out.cloneVector("");
    MPI_Comm comm;
    int ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(d_petsc_mat), &comm);
    IBTK_CHKERRQ(ierr);
    d_petsc_x = PETScSAMRAIVectorReal::createPETScVector(d_x, comm);
    d_petsc_y = PETScSAMRAIVectorReal::createPETScVector(d_y, comm);
    d_petsc_z = PETScSAMRAIVectorReal::createPETScVector(d_z, comm);
    d_is_initialized = true;
    return;
} // initializeOperatorState

void
PETScMatLOWrapper::deallocateOperatorState()
{
    if (!d_is_initialized) return;
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_x);
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_y);
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_z);
    d_x->freeVectorComponents();
    d_x.setNull();
    d_y->freeVectorComponents();
    d_y.setNull();
    d_z->freeVectorComponents();
    d_z.setNull();
    d_is_initialized = false;
    return;
} // deallocateOperatorState

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
