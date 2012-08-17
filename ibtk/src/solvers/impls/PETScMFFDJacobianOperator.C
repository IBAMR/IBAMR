// Filename: PETScMFFDJacobianOperator.C
// Created on 27 Aug 2010 by Boyce Griffith
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

#include "PETScMFFDJacobianOperator.h"

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
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScMFFDJacobianOperator::PETScMFFDJacobianOperator(
    const std::string& options_prefix)
    : JacobianOperator(),
      d_F(NULL),
      d_nonlinear_solver(NULL),
      d_petsc_jac(PETSC_NULL),
      d_op_u(NULL),
      d_op_x(NULL),
      d_op_y(NULL),
      d_petsc_u(PETSC_NULL),
      d_petsc_x(PETSC_NULL),
      d_petsc_y(PETSC_NULL),
      d_options_prefix(options_prefix),
      d_is_initialized(false)
{
    // intentionally blank
    return;
}// PETScMFFDJacobianOperator()

PETScMFFDJacobianOperator::~PETScMFFDJacobianOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}// ~PETScMFFDJacobianOperator()

void
PETScMFFDJacobianOperator::setOperator(
    Pointer<GeneralOperator> F)
{
    d_F = F;
    return;
}// setOperator

void
PETScMFFDJacobianOperator::setNewtonKrylovSolver(
    Pointer<PETScNewtonKrylovSolver> nonlinear_solver)
{
    d_nonlinear_solver = nonlinear_solver;
    return;
}// setNewtonKrylovSolver

void
PETScMFFDJacobianOperator::formJacobian(
    SAMRAIVectorReal<NDIM,double>& u)
{
    int ierr;
    if (d_nonlinear_solver.isNull())
    {
        d_op_u->copyVector(Pointer<SAMRAIVectorReal<NDIM,double> >(&u,false), false);
        ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(d_petsc_u)); IBTK_CHKERRQ(ierr);
        ierr = MatMFFDSetBase(d_petsc_jac, d_petsc_u, PETSC_NULL); IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyBegin(d_petsc_jac, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(d_petsc_jac, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    }
    else
    {
        SNES snes = d_nonlinear_solver->getPETScSNES();
        Vec u, f;
        ierr = SNESGetSolution(snes, &u); IBTK_CHKERRQ(ierr);
        ierr = SNESGetFunction(snes, &f, PETSC_NULL, PETSC_NULL); IBTK_CHKERRQ(ierr);
        ierr = MatMFFDSetBase(d_petsc_jac, u, f); IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyBegin(d_petsc_jac, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(d_petsc_jac, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    }
    return;
}// formJacobian

Pointer<SAMRAIVectorReal<NDIM,double> >
PETScMFFDJacobianOperator::getBaseVector() const
{
    if (d_nonlinear_solver.isNull())
    {
        return d_op_u;
    }
    else
    {
        SNES snes = d_nonlinear_solver->getPETScSNES();
        Vec u;
        int ierr = SNESGetSolution(snes, &u); IBTK_CHKERRQ(ierr);
        return PETScSAMRAIVectorReal::getSAMRAIVector(u);
    }
    return Pointer<SAMRAIVectorReal<NDIM,double> >(NULL);
}// getBaseVector

void
PETScMFFDJacobianOperator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    // Compute the action of the operator.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, Pointer<SAMRAIVectorReal<NDIM,PetscScalar> >(&x,false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_y, Pointer<SAMRAIVectorReal<NDIM,PetscScalar> >(&y,false));
    int ierr = MatMult(d_petsc_jac, d_petsc_x, d_petsc_y); IBTK_CHKERRQ(ierr);
    return;
}// apply

void
PETScMFFDJacobianOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    if (d_is_initialized) deallocateOperatorState();

    int ierr;
    MPI_Comm comm = PETSC_COMM_WORLD;

    // Setup the matrix-free matrix.
    ierr = MatCreateMFFD(comm, 0, 0, 0, 0, &d_petsc_jac); IBTK_CHKERRQ(ierr);
    ierr = MatMFFDSetFunction(d_petsc_jac, reinterpret_cast<PetscErrorCode(*)(void*, Vec, Vec)>(FormFunction_SAMRAI), this); IBTK_CHKERRQ(ierr);
    if (!d_options_prefix.empty())
    {
        ierr = MatSetOptionsPrefix(d_petsc_jac, d_options_prefix.c_str()); IBTK_CHKERRQ(ierr);
    }
    ierr = MatSetFromOptions(d_petsc_jac); IBTK_CHKERRQ(ierr);

    // Setup solution and rhs vectors.
    d_op_u = in.cloneVector(in.getName());
    d_op_u->allocateVectorData();
    d_petsc_u = PETScSAMRAIVectorReal::createPETScVector(d_op_u, comm);

    d_op_x = in.cloneVector(in.getName());
    d_petsc_x = PETScSAMRAIVectorReal::createPETScVector(d_op_x, comm);

    d_op_y = out.cloneVector(out.getName());
    d_petsc_y = PETScSAMRAIVectorReal::createPETScVector(d_op_y, comm);

    // Indicate that the operator is initialized.
    d_is_initialized = true;
    return;
}// initializeOperatorState

void
PETScMFFDJacobianOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_u);
    d_petsc_u = PETSC_NULL;
    d_op_u->resetLevels(d_op_u->getCoarsestLevelNumber(), std::min(d_op_u->getFinestLevelNumber(),d_op_u->getPatchHierarchy()->getFinestLevelNumber()));
    d_op_u->freeVectorComponents();
    d_op_u.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_x);
    d_petsc_x = PETSC_NULL;
    d_op_x->resetLevels(d_op_x->getCoarsestLevelNumber(), std::min(d_op_x->getFinestLevelNumber(),d_op_x->getPatchHierarchy()->getFinestLevelNumber()));
    d_op_x->freeVectorComponents();
    d_op_x.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_y);
    d_petsc_y = PETSC_NULL;
    d_op_y->resetLevels(d_op_y->getCoarsestLevelNumber(), std::min(d_op_y->getFinestLevelNumber(),d_op_y->getPatchHierarchy()->getFinestLevelNumber()));
    d_op_y->freeVectorComponents();
    d_op_y.setNull();

    int ierr = MatDestroy(&d_petsc_jac); IBTK_CHKERRQ(ierr);
    d_petsc_jac = PETSC_NULL;

    d_is_initialized = false;
    return;
}// deallocateOperatorState

void
PETScMFFDJacobianOperator::enableLogging(
    bool /*enabled*/)
{
    // intentionally blank
    return;
}// enableLogging

/////////////////////////////// PRIVATE //////////////////////////////////////

PetscErrorCode
PETScMFFDJacobianOperator::FormFunction_SAMRAI(
    void* p_ctx,
    Vec x,
    Vec f)
{
    PETScMFFDJacobianOperator* jac_op = static_cast<PETScMFFDJacobianOperator*>(p_ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(jac_op != NULL);
    TBOX_ASSERT(!jac_op->d_F.isNull());
#endif
    int ierr;
    jac_op->d_F->apply(*PETScSAMRAIVectorReal::getSAMRAIVector(x), *PETScSAMRAIVectorReal::getSAMRAIVector(f));
    if (!jac_op->d_nonlinear_solver.isNull())
    {
        SNES snes = jac_op->d_nonlinear_solver->getPETScSNES();
        Vec rhs;
        ierr = SNESGetRhs(snes, &rhs); IBTK_CHKERRQ(ierr);
        if (rhs != PETSC_NULL)
        {
            VecAXPY(f, -1.0, rhs); IBTK_CHKERRQ(ierr);
        }
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f)); IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// FormFunction_SAMRAI

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
