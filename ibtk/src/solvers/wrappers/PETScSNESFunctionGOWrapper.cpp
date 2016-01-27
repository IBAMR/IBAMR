// Filename: PETScSNESFunctionGOWrapper.cpp
// Created on 27 Dec 2003 by Boyce Griffith
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
#include "ibtk/GeneralOperator.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/PETScSNESFunctionGOWrapper.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "mpi.h"
#include "petscsnes.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScSNESFunctionGOWrapper::PETScSNESFunctionGOWrapper(
    const std::string& object_name,
    const SNES& petsc_snes,
    PetscErrorCode (*const petsc_snes_form_func)(SNES, Vec, Vec, void*),
    void* const petsc_snes_func_ctx)
    : GeneralOperator(object_name),
      d_petsc_snes(petsc_snes),
      d_petsc_snes_form_func(petsc_snes_form_func),
      d_petsc_snes_func_ctx(petsc_snes_func_ctx),
      d_x(NULL),
      d_y(NULL),
      d_petsc_x(NULL),
      d_petsc_y(NULL)
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
