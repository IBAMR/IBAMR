// Filename: PETScMFFDJacobianOperator.cpp
// Created on 27 Aug 2010 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#include "ibtk/GeneralOperator.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/JacobianOperator.h"
#include "ibtk/PETScMFFDJacobianOperator.h"
#include "ibtk/PETScNewtonKrylovSolver.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include "petscerror.h"
#include "petscmat.h"
#include "petscmath.h"
#include "petscsnes.h"
#include "petscsys.h"
#include "petscvec.h"

#include <mpi.h>

#include <algorithm>
#include <ostream>
#include <string>
#include <utility>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScMFFDJacobianOperator::PETScMFFDJacobianOperator(std::string object_name, std::string options_prefix)
    : JacobianOperator(std::move(object_name)), d_options_prefix(std::move(options_prefix))
{
    // intentionally blank
}

PETScMFFDJacobianOperator::~PETScMFFDJacobianOperator()
{
    if (d_is_initialized) deallocateOperatorState();
}

void
PETScMFFDJacobianOperator::setOperator(Pointer<GeneralOperator> F)
{
    d_F = F;
}

void
PETScMFFDJacobianOperator::setNewtonKrylovSolver(Pointer<PETScNewtonKrylovSolver> nonlinear_solver)
{
    d_nonlinear_solver = nonlinear_solver;
}

void
PETScMFFDJacobianOperator::formJacobian(SAMRAIVectorReal<NDIM, double>& u)
{
    int ierr;
    if (d_nonlinear_solver)
    {
        SNES snes = d_nonlinear_solver->getPETScSNES();
        Vec u, f;
        ierr = SNESGetSolution(snes, &u);
        IBTK_CHKERRQ(ierr);
        ierr = SNESGetFunction(snes, &f, nullptr, nullptr);
        IBTK_CHKERRQ(ierr);
        ierr = MatMFFDSetBase(d_petsc_jac, u, f);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyBegin(d_petsc_jac, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(d_petsc_jac, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
    }
    else
    {
        d_op_u->allocateVectorData();
        d_op_u->copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&u, false), false);
        PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_u, d_op_u);
        ierr = MatMFFDSetBase(d_petsc_jac, d_petsc_u, nullptr);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyBegin(d_petsc_jac, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(d_petsc_jac, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        d_op_u->deallocateVectorData();
    }
}

Pointer<SAMRAIVectorReal<NDIM, double> >
PETScMFFDJacobianOperator::getBaseVector() const
{
    if (d_nonlinear_solver)
    {
        SNES snes = d_nonlinear_solver->getPETScSNES();
        Vec u;
        int ierr = SNESGetSolution(snes, &u);
        IBTK_CHKERRQ(ierr);
        Pointer<SAMRAIVectorReal<NDIM, double> > samrai_u;
        PETScSAMRAIVectorReal::getSAMRAIVector(u, &samrai_u);
        Pointer<SAMRAIVectorReal<NDIM, double> > samrai_u_ptr = samrai_u;
        PETScSAMRAIVectorReal::restoreSAMRAIVector(u, &samrai_u);
        return samrai_u_ptr;
    }
    else
    {
        return d_op_u;
    }
    return Pointer<SAMRAIVectorReal<NDIM, double> >(nullptr);
}

void
PETScMFFDJacobianOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    // Compute the action of the operator.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, Pointer<SAMRAIVectorReal<NDIM, PetscScalar> >(&x, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_y, Pointer<SAMRAIVectorReal<NDIM, PetscScalar> >(&y, false));
    int ierr = MatMult(d_petsc_jac, d_petsc_x, d_petsc_y);
    IBTK_CHKERRQ(ierr);
}

void
PETScMFFDJacobianOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                   const SAMRAIVectorReal<NDIM, double>& out)
{
    if (d_is_initialized) deallocateOperatorState();

    int ierr;
    MPI_Comm comm = PETSC_COMM_WORLD;

    // Setup the matrix-free matrix.
    ierr = MatCreateMFFD(comm, 1, 1, PETSC_DETERMINE, PETSC_DETERMINE, &d_petsc_jac);
    IBTK_CHKERRQ(ierr);
    ierr = MatMFFDSetFunction(
        d_petsc_jac, reinterpret_cast<PetscErrorCode (*)(void*, Vec, Vec)>(FormFunction_SAMRAI), this);
    IBTK_CHKERRQ(ierr);
    if (!d_options_prefix.empty())
    {
        ierr = MatSetOptionsPrefix(d_petsc_jac, d_options_prefix.c_str());
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatSetFromOptions(d_petsc_jac);
    IBTK_CHKERRQ(ierr);

    // Setup solution and rhs vectors.
    d_op_u = in.cloneVector(in.getName());
    d_petsc_u = PETScSAMRAIVectorReal::createPETScVector(d_op_u, comm);

    d_op_x = in.cloneVector(in.getName());
    d_petsc_x = PETScSAMRAIVectorReal::createPETScVector(d_op_x, comm);

    d_op_y = out.cloneVector(out.getName());
    d_petsc_y = PETScSAMRAIVectorReal::createPETScVector(d_op_y, comm);

    // Indicate that the operator is initialized.
    d_is_initialized = true;
}

void
PETScMFFDJacobianOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_u);
    d_petsc_u = nullptr;
    d_op_u->resetLevels(0,
                        std::min(d_op_u->getFinestLevelNumber(), d_op_u->getPatchHierarchy()->getFinestLevelNumber()));
    d_op_u->freeVectorComponents();
    d_op_u.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_x);
    d_petsc_x = nullptr;
    d_op_x->freeVectorComponents();
    d_op_x.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_y);
    d_petsc_y = nullptr;
    d_op_y->freeVectorComponents();
    d_op_y.setNull();

    int ierr = MatDestroy(&d_petsc_jac);
    IBTK_CHKERRQ(ierr);
    d_petsc_jac = nullptr;

    d_is_initialized = false;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

PetscErrorCode
PETScMFFDJacobianOperator::FormFunction_SAMRAI(void* p_ctx, Vec x, Vec f)
{
    auto jac_op = static_cast<PETScMFFDJacobianOperator*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(jac_op);
    TBOX_ASSERT(jac_op->d_F);
#endif
    int ierr;
    Pointer<SAMRAIVectorReal<NDIM, double> > samrai_x, samrai_f;
    PETScSAMRAIVectorReal::getSAMRAIVectorRead(x, &samrai_x);
    PETScSAMRAIVectorReal::getSAMRAIVector(f, &samrai_f);
    jac_op->d_F->apply(*samrai_x, *samrai_f);
    PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(x, &samrai_x);
    PETScSAMRAIVectorReal::restoreSAMRAIVector(f, &samrai_f);
    if (jac_op->d_nonlinear_solver)
    {
        SNES snes = jac_op->d_nonlinear_solver->getPETScSNES();
        Vec rhs;
        ierr = SNESGetRhs(snes, &rhs);
        CHKERRQ(ierr);
        if (rhs)
        {
            ierr = VecAXPY(f, -1.0, rhs);
            CHKERRQ(ierr);
        }
    }
    PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
