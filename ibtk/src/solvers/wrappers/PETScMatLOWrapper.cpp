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

#include "SAMRAI/hier/IntVector.h"

#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/PETScMatLOWrapper.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "mpi.h"
#include "petscmat.h"
#include "petscsys.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScMatLOWrapper::PETScMatLOWrapper(const std::string& object_name, const Mat& petsc_mat)
    : LinearOperator(object_name), d_petsc_mat(petsc_mat), d_x(NULL), d_y(NULL), d_z(NULL), d_petsc_x(NULL),
      d_petsc_y(NULL), d_petsc_z(NULL)
{
    // intentionally blank
    return;
}

PETScMatLOWrapper::~PETScMatLOWrapper()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}

const Mat& PETScMatLOWrapper::getPETScMat() const
{
    return d_petsc_mat;
}

void PETScMatLOWrapper::apply(SAMRAIVectorReal<double>& x, SAMRAIVectorReal<double>& y)
{
    if (!d_is_initialized) initializeOperatorState(x, y);

    // Update the PETSc Vec wrappers.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x,
                                               boost::shared_ptr<SAMRAIVectorReal<double> >(&x, NullDeleter()));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_y,
                                               boost::shared_ptr<SAMRAIVectorReal<double> >(&y, NullDeleter()));

    // Apply the operator.
    int ierr = MatMult(d_petsc_mat, d_petsc_x, d_petsc_y);
    IBTK_CHKERRQ(ierr);
    return;
}

void PETScMatLOWrapper::applyAdd(SAMRAIVectorReal<double>& x, SAMRAIVectorReal<double>& y, SAMRAIVectorReal<double>& z)
{
    if (!d_is_initialized) initializeOperatorState(x, y);

    // Update the PETSc Vec wrappers.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x,
                                               boost::shared_ptr<SAMRAIVectorReal<double> >(&x, NullDeleter()));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_y,
                                               boost::shared_ptr<SAMRAIVectorReal<double> >(&y, NullDeleter()));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_z,
                                               boost::shared_ptr<SAMRAIVectorReal<double> >(&z, NullDeleter()));

    // Apply the operator.
    int ierr = MatMultAdd(d_petsc_mat, d_petsc_x, d_petsc_y, d_petsc_z);
    IBTK_CHKERRQ(ierr);
    return;
}

void PETScMatLOWrapper::initializeOperatorState(const SAMRAIVectorReal<double>& in, const SAMRAIVectorReal<double>& out)
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
}

void PETScMatLOWrapper::deallocateOperatorState()
{
    if (!d_is_initialized) return;
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_x);
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_y);
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_z);
    d_x->freeVectorComponents();
    d_x.reset();
    d_y->freeVectorComponents();
    d_y.reset();
    d_z->freeVectorComponents();
    d_z.reset();
    d_is_initialized = false;
    return;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
