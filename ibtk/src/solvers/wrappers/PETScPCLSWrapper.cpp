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

#include "ibtk/GeneralSolver.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/PETScPCLSWrapper.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#include "petscpc.h"
#include "petscpctypes.h"
#include <petsclog.h>

#include <mpi.h>

#include <ostream>
#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScPCLSWrapper::PETScPCLSWrapper(std::string object_name, PC petsc_pc) : d_petsc_pc(std::move(petsc_pc))
{
    GeneralSolver::init(std::move(object_name), /*homogeneous_bc*/ true);
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
