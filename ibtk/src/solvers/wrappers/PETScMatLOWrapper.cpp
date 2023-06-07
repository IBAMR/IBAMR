// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/PETScMatLOWrapper.h"
#include "ibtk/PETScSAMRAIVectorReal.h"

#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAIVectorReal.h"

#include "petscmat.h"

#include <mpi.h>

#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScMatLOWrapper::PETScMatLOWrapper(std::string object_name, Mat petsc_mat)
    : LinearOperator(std::move(object_name)), d_petsc_mat(std::move(petsc_mat))
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
    free_vector_components(*d_x);
    d_x.setNull();
    free_vector_components(*d_y);
    d_y.setNull();
    free_vector_components(*d_z);
    d_z.setNull();
    d_is_initialized = false;
    return;
} // deallocateOperatorState

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
