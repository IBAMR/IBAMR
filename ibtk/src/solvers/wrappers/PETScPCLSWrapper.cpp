// Filename: PETScPCLSWrapper.cpp
// Created on 11 Nov 2004 by Boyce Griffith
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
#include <ostream>
#include <string>

#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScPCLSWrapper.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "mpi.h"
#include "petscpc.h"
#include "petscsys.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScPCLSWrapper::PETScPCLSWrapper(const std::string& object_name, const PC& petsc_pc)
    : d_petsc_pc(petsc_pc), d_x(NULL), d_b(NULL), d_petsc_x(NULL), d_petsc_b(NULL)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ true);
    return;
} // PETScPCLSWrapper()

PETScPCLSWrapper::~PETScPCLSWrapper()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~PETScPCLSWrapper()

const PC&
PETScPCLSWrapper::getPETScPC() const
{
    return d_petsc_pc;
} // getPETScPC

bool
PETScPCLSWrapper::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    if (!d_is_initialized) initializeSolverState(x, b);

    // Update the PETSc Vec wrappers.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_b, Pointer<SAMRAIVectorReal<NDIM, double> >(&b, false));

    // Apply the preconditioner.
    int ierr = PCApply(d_petsc_pc, d_petsc_x, d_petsc_b);
    IBTK_CHKERRQ(ierr);
    return true;
} // solveSystem

void
PETScPCLSWrapper::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                        const SAMRAIVectorReal<NDIM, double>& b)
{
    if (d_is_initialized) deallocateSolverState();
    d_x = x.cloneVector("");
    d_b = b.cloneVector("");
    MPI_Comm comm;
    int ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(d_petsc_pc), &comm);
    IBTK_CHKERRQ(ierr);
    d_petsc_x = PETScSAMRAIVectorReal::createPETScVector(d_x, comm);
    d_petsc_b = PETScSAMRAIVectorReal::createPETScVector(d_b, comm);
    d_is_initialized = true;
    return;
} // initializeSolverState

void
PETScPCLSWrapper::deallocateSolverState()
{
    if (!d_is_initialized) return;
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_x);
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_b);
    d_x->freeVectorComponents();
    d_x.setNull();
    d_b->freeVectorComponents();
    d_b.setNull();
    d_is_initialized = false;
    return;
} // deallocateSolverState

void
PETScPCLSWrapper::setInitialGuessNonzero(bool /*initial_guess_nonzero*/)
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::setInitialGuessNonzero() not supported" << std::endl););
    return;
} // setInitialGuessNonzero

bool
PETScPCLSWrapper::getInitialGuessNonzero() const
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::getInitialGuessNonzero() not supported" << std::endl););
    return true;
} // getInitialGuessNonzero

void
PETScPCLSWrapper::setMaxIterations(int /*max_iterations*/)
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::setMaxIterations() not supported" << std::endl););
    return;
} // setMaxIterations

int
PETScPCLSWrapper::getMaxIterations() const
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::getMaxIterations() not supported" << std::endl););
    return 0;
} // getMaxIterations

void
PETScPCLSWrapper::setAbsoluteTolerance(double /*abs_residual_tol*/)
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::setAbsoluteTolerance() not supported" << std::endl););
    return;
} // setAbsoluteTolerance

double
PETScPCLSWrapper::getAbsoluteTolerance() const
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::getAbsoluteTolerance() not supported" << std::endl););
    return 0.0;
} // getAbsoluteTolerance

void
PETScPCLSWrapper::setRelativeTolerance(double /*rel_residual_tol*/)
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::setRelativeTolerance() not supported" << std::endl););
    return;
} // setRelativeTolerance

double
PETScPCLSWrapper::getRelativeTolerance() const
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::getRelativeTolerance() not supported" << std::endl););
    return 0.0;
} // getRelativeTolerance

int
PETScPCLSWrapper::getNumIterations() const
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::getNumIteratons() not supported" << std::endl););
    return 0;
} // getNumIterations

double
PETScPCLSWrapper::getResidualNorm() const
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::getResidualNorm() not supported" << std::endl););
    return 0.0;
} // getResidualNorm

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
