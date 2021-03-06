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

#include "ibtk/GeneralOperator.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/PETScSNESFunctionGOWrapper.h"

#include "Box.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#include "petscvec.h"

#include <mpi.h>

#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScSNESFunctionGOWrapper::PETScSNESFunctionGOWrapper(
    std::string object_name,
    SNES petsc_snes,
    PetscErrorCode (*const petsc_snes_form_func)(SNES, Vec, Vec, void*),
    void* const petsc_snes_func_ctx)
    : GeneralOperator(std::move(object_name)),
      d_petsc_snes(std::move(petsc_snes)),
      d_petsc_snes_form_func(petsc_snes_form_func),
      d_petsc_snes_func_ctx(petsc_snes_func_ctx)
{
    // intentionally blank
    return;
} // PETScSNESFunctionGOWrapper()

PETScSNESFunctionGOWrapper::~PETScSNESFunctionGOWrapper()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~PETScSNESFunctionGOWrapper()

const SNES&
PETScSNESFunctionGOWrapper::getPETScSNES() const
{
    return d_petsc_snes;
} // getPETScSNES

PetscErrorCode (*PETScSNESFunctionGOWrapper::getPETScSNESFormFunction())(SNES, Vec, Vec, void*)
{
    return d_petsc_snes_form_func;
} // getPETScSNESFormFunction

void*
PETScSNESFunctionGOWrapper::getPETScSNESFunctionContext() const
{
    return d_petsc_snes_func_ctx;
} // getPETScSNESFunctionContext

void
PETScSNESFunctionGOWrapper::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    if (!d_is_initialized) initializeOperatorState(x, y);

    // Update the PETSc Vec wrappers.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_y, Pointer<SAMRAIVectorReal<NDIM, double> >(&y, false));

    // Apply the operator.
    int ierr = d_petsc_snes_form_func(d_petsc_snes, d_petsc_x, d_petsc_y, d_petsc_snes_func_ctx);
    IBTK_CHKERRQ(ierr);
    return;
} // apply

void
PETScSNESFunctionGOWrapper::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                    const SAMRAIVectorReal<NDIM, double>& out)
{
    if (d_is_initialized) deallocateOperatorState();
    d_x = in.cloneVector("");
    d_y = out.cloneVector("");
    MPI_Comm comm;
    int ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(d_petsc_snes), &comm);
    IBTK_CHKERRQ(ierr);
    d_petsc_x = PETScSAMRAIVectorReal::createPETScVector(d_x, comm);
    d_petsc_y = PETScSAMRAIVectorReal::createPETScVector(d_y, comm);
    d_is_initialized = true;
    return;
} // initializeOperatorState

void
PETScSNESFunctionGOWrapper::deallocateOperatorState()
{
    if (!d_is_initialized) return;
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_x);
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_y);
    d_x->freeVectorComponents();
    d_x.setNull();
    d_y->freeVectorComponents();
    d_y.setNull();
    d_is_initialized = false;
    return;
} // deallocateOperatorState

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
