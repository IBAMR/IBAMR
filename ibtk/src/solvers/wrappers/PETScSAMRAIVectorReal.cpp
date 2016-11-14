// Filename: PETScSAMRAIVectorReal.cpp
// Created on 10 Nov 2004 by Boyce Griffith
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

#include <math.h>
#include <algorithm>
#include <ostream>

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/NormOps.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "mpi.h"
#include "petscerror.h"
#include "petscis.h"
#include "petscmath.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_vec_duplicate;
static Timer* t_vec_dot;
static Timer* t_vec_m_dot;
static Timer* t_vec_norm;
static Timer* t_vec_t_dot;
static Timer* t_vec_m_t_dot;
static Timer* t_vec_scale;
static Timer* t_vec_copy;
static Timer* t_vec_set;
static Timer* t_vec_swap;
static Timer* t_vec_axpy;
static Timer* t_vec_axpby;
static Timer* t_vec_maxpy;
static Timer* t_vec_aypx;
static Timer* t_vec_waxpy;
static Timer* t_vec_axpbypcz;
static Timer* t_vec_pointwise_mult;
static Timer* t_vec_pointwise_divide;
static Timer* t_vec_get_size;
static Timer* t_vec_get_local_size;
static Timer* t_vec_max;
static Timer* t_vec_min;
static Timer* t_vec_set_random;
static Timer* t_vec_destroy;
static Timer* t_vec_dot_local;
static Timer* t_vec_t_dot_local;
static Timer* t_vec_norm_local;
static Timer* t_vec_m_dot_local;
static Timer* t_vec_m_t_dot_local;
static Timer* t_vec_max_pointwise_divide;
static Timer* t_vec_dot_norm2;

#define PSVR_CAST1(v) (static_cast<PETScSAMRAIVectorReal*>(v->data))
#define PSVR_CAST2(v) (static_cast<PETScSAMRAIVectorReal*>(v->data)->d_samrai_vector)

#if !defined(NDEBUG)
#define PSVR_CHECK1(v)                                                                                                 \
    TBOX_ASSERT((v));                                                                                                  \
    TBOX_ASSERT(!PSVR_CAST1((v))->d_vector_checked_out_read);
#define PSVR_CHECK2(v1, v2)                                                                                            \
    PSVR_CHECK1((v1));                                                                                                 \
    PSVR_CHECK1((v2));
#define PSVR_CHECK3(v1, v2, v3)                                                                                        \
    PSVR_CHECK1((v1));                                                                                                 \
    PSVR_CHECK1((v2));                                                                                                 \
    PSVR_CHECK1((v3));
#define PSVR_CHECKN(v, N)                                                                                              \
    for (int i = 0; i < static_cast<int>(N); ++i)                                                                      \
    {                                                                                                                  \
        PSVR_CHECK1((v)[i]);                                                                                           \
    }
#else
#define PSVR_CHECK1(v)
#define PSVR_CHECK2(v1, v2)
#define PSVR_CHECKN(v, N)
#endif
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

PETScSAMRAIVectorReal::PETScSAMRAIVectorReal(Pointer<SAMRAIVectorReal<NDIM, PetscScalar> > samrai_vector,
                                             bool vector_created_via_duplicate,
                                             MPI_Comm comm)
    : d_samrai_vector(samrai_vector),
      d_vector_created_via_duplicate(vector_created_via_duplicate),
      d_vector_checked_out_read_write(false),
      d_vector_checked_out_read(false)
{
    // Setup Timers.
    IBTK_DO_ONCE(
        t_vec_duplicate = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDuplicate()");
        t_vec_dot = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDot()");
        t_vec_m_dot = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMDot()");
        t_vec_norm = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecNorm()");
        t_vec_t_dot = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecTDot()");
        t_vec_m_t_dot = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMTDot()");
        t_vec_scale = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecScale()");
        t_vec_copy = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecCopy()");
        t_vec_set = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecSet()");
        t_vec_swap = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecSwap()");
        t_vec_axpy = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAXPY()");
        t_vec_axpby = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAXPBY()");
        t_vec_maxpy = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMAXPY()");
        t_vec_aypx = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAYPX()");
        t_vec_waxpy = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecWAXPY()");
        t_vec_axpbypcz = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAXPBYPCZ()");
        t_vec_pointwise_mult = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecPointwiseMult()");
        t_vec_pointwise_divide =
            TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecPointwiseDivide()");
        t_vec_get_size = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecGetSize()");
        t_vec_get_local_size = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecGetLocalSize()");
        t_vec_max = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMax()");
        t_vec_min = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMin()");
        t_vec_set_random = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecSetRandom()");
        t_vec_destroy = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDestroy()");
        t_vec_dot_local = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDot_local()");
        t_vec_t_dot_local = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecTDot_local()");
        t_vec_norm_local = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecNorm_local()");
        t_vec_m_dot_local = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMDot_local()");
        t_vec_m_t_dot_local = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMTDot_local()");
        t_vec_max_pointwise_divide =
            TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMaxPointwiseDivide()");
        t_vec_dot_norm2 = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDotNorm2()"););

    int ierr;
    ierr = VecCreate(comm, &d_petsc_vector);
    IBTK_CHKERRQ(ierr);

    // Assign vector operations to PETSc vector object.
    static struct _VecOps DvOps;
    IBTK_DO_ONCE(DvOps.duplicate = PETScSAMRAIVectorReal::VecDuplicate_SAMRAI;
                 DvOps.duplicatevecs = VecDuplicateVecs_SAMRAI;
                 DvOps.destroyvecs = VecDestroyVecs_SAMRAI;
                 DvOps.dot = VecDot_SAMRAI;
                 DvOps.mdot = VecMDot_SAMRAI;
                 DvOps.norm = VecNorm_SAMRAI;
                 DvOps.tdot = VecTDot_SAMRAI;
                 DvOps.mtdot = VecMTDot_SAMRAI;
                 DvOps.scale = VecScale_SAMRAI;
                 DvOps.copy = VecCopy_SAMRAI;
                 DvOps.set = VecSet_SAMRAI;
                 DvOps.swap = VecSwap_SAMRAI;
                 DvOps.axpy = VecAXPY_SAMRAI;
                 DvOps.axpby = VecAXPBY_SAMRAI;
                 DvOps.maxpy = VecMAXPY_SAMRAI;
                 DvOps.aypx = VecAYPX_SAMRAI;
                 DvOps.waxpy = VecWAXPY_SAMRAI;
                 DvOps.axpbypcz = VecAXPBYPCZ_SAMRAI;
                 DvOps.pointwisemult = VecPointwiseMult_SAMRAI;
                 DvOps.pointwisedivide = VecPointwiseDivide_SAMRAI;
                 DvOps.getsize = VecGetSize_SAMRAI;
                 DvOps.getlocalsize = VecGetLocalSize_SAMRAI;
                 DvOps.max = VecMax_SAMRAI;
                 DvOps.min = VecMin_SAMRAI;
                 DvOps.setrandom = VecSetRandom_SAMRAI;
                 DvOps.destroy = PETScSAMRAIVectorReal::VecDestroy_SAMRAI;
                 DvOps.dot_local = VecDot_local_SAMRAI;
                 DvOps.tdot_local = VecTDot_local_SAMRAI;
                 DvOps.norm_local = VecNorm_local_SAMRAI;
                 DvOps.mdot_local = VecMDot_local_SAMRAI;
                 DvOps.mtdot_local = VecMTDot_local_SAMRAI;
                 DvOps.maxpointwisedivide = VecMaxPointwiseDivide_SAMRAI;
                 DvOps.dotnorm2 = VecDotNorm2_SAMRAI;);
    ierr = PetscMemcpy(d_petsc_vector->ops, &DvOps, sizeof(DvOps));
    IBTK_CHKERRQ(ierr);

    // Set PETSc vector data.
    d_petsc_vector->data = this;
    d_petsc_vector->petscnative = PETSC_FALSE;
    int size;
    MPI_Comm_size(comm, &size);
    d_petsc_vector->map->n = 1;    // NOTE: Here we are giving a bogus local  size.
    d_petsc_vector->map->N = size; // NOTE: Here we are giving a bogus global size.
    d_petsc_vector->map->bs = 1;   // NOTE: Here we are giving a bogus block  size.

    // Set the PETSc vector type name.
    ierr = PetscObjectChangeTypeName(reinterpret_cast<PetscObject>(d_petsc_vector), "Vec_SAMRAI");
    IBTK_CHKERRQ(ierr);

    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(d_petsc_vector));
    IBTK_CHKERRQ(ierr);
}

PETScSAMRAIVectorReal::~PETScSAMRAIVectorReal()
{
    if (!d_vector_created_via_duplicate)
    {
        d_petsc_vector->ops->destroy = 0;
        int ierr = VecDestroy(&d_petsc_vector);
        IBTK_CHKERRQ(ierr);
    }
}

/////////////////////////////// PRIVATE //////////////////////////////////////

PetscErrorCode
PETScSAMRAIVectorReal::VecDuplicate_SAMRAI(Vec v, Vec* newv)
{
    IBTK_TIMER_START(t_vec_duplicate);
    PSVR_CHECK1(v);
    PetscErrorCode ierr;
    Pointer<SAMRAIVectorReal<NDIM, PetscScalar> > samrai_vec = PSVR_CAST2(v)->cloneVector(PSVR_CAST2(v)->getName());
    samrai_vec->allocateVectorData();
    static const bool vector_created_via_duplicate = true;
    MPI_Comm comm;
    ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(v), &comm);
    CHKERRQ(ierr);
    PETScSAMRAIVectorReal* new_psv = new PETScSAMRAIVectorReal(samrai_vec, vector_created_via_duplicate, comm);
    *newv = new_psv->d_petsc_vector;
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*newv));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_duplicate);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDestroy_SAMRAI(Vec v)
{
    IBTK_TIMER_START(t_vec_destroy);
    PSVR_CHECK1(v);
    if (PSVR_CAST1(v)->d_vector_created_via_duplicate)
    {
        PSVR_CAST2(v)->resetLevels(0,
                                   std::min(PSVR_CAST2(v)->getFinestLevelNumber(),
                                            PSVR_CAST2(v)->getPatchHierarchy()->getFinestLevelNumber()));
        PSVR_CAST2(v)->deallocateVectorData();
        PSVR_CAST2(v)->freeVectorComponents();
        PSVR_CAST2(v).setNull();
        destroyPETScVector(PSVR_CAST1(v)->d_petsc_vector);
    }
    IBTK_TIMER_STOP(t_vec_destroy);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDuplicateVecs_SAMRAI(Vec v, PetscInt m, Vec* V[])
{
    PSVR_CHECK1(v);
    PetscErrorCode ierr;
    ierr = PetscMalloc1(m, V);
    CHKERRQ(ierr);
    for (PetscInt i = 0; i < m; ++i)
    {
        ierr = VecDuplicate(v, *V + i);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDestroyVecs_SAMRAI(PetscInt m, Vec vv[])
{
    PSVR_CHECKN(vv, m);
    PetscErrorCode ierr;
    for (PetscInt i = 0; i < m; ++i)
    {
        ierr = VecDestroy(&vv[i]);
        CHKERRQ(ierr);
    }
    ierr = PetscFree(vv);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDot_SAMRAI(Vec x, Vec y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_dot);
    PSVR_CHECK2(x, y);
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y));
    IBTK_TIMER_STOP(t_vec_dot);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMDot_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_m_dot);
    PSVR_CHECK1(x);
    PSVR_CHECKN(y, nv);
    static const bool local_only = true;
    for (PetscInt i = 0; i < nv; ++i)
    {
        val[i] = PSVR_CAST2(x)->dot(PSVR_CAST2(y[i]), local_only);
    }
    SAMRAI_MPI::sumReduction(val, nv);
    IBTK_TIMER_STOP(t_vec_m_dot);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecNorm_SAMRAI(Vec x, NormType type, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_norm);
    PSVR_CHECK1(x);
    if (type == NORM_1)
    {
        *val = NormOps::L1Norm(PSVR_CAST2(x));
    }
    else if (type == NORM_2)
    {
        *val = NormOps::L2Norm(PSVR_CAST2(x));
    }
    else if (type == NORM_INFINITY)
    {
        *val = NormOps::maxNorm(PSVR_CAST2(x));
    }
    else if (type == NORM_1_AND_2)
    {
        static const bool local_only = true;
        val[0] = NormOps::L1Norm(PSVR_CAST2(x), local_only);
        val[1] = NormOps::L2Norm(PSVR_CAST2(x), local_only);
        val[1] = val[1] * val[1];
        SAMRAI_MPI::sumReduction(val, 2);
        val[1] = std::sqrt(val[1]);
    }
    else
    {
        TBOX_ERROR("PETScSAMRAIVectorReal::norm()\n"
                   << "  vector norm type "
                   << static_cast<int>(type)
                   << " unsupported"
                   << std::endl);
    }
    IBTK_TIMER_STOP(t_vec_norm);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecTDot_SAMRAI(Vec x, Vec y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_t_dot);
    PSVR_CHECK2(x, y);
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y));
    IBTK_TIMER_STOP(t_vec_t_dot);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMTDot_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_m_t_dot);
    PSVR_CHECK1(x);
    PSVR_CHECKN(y, nv);
    static const bool local_only = true;
    for (PetscInt i = 0; i < nv; ++i)
    {
        val[i] = PSVR_CAST2(x)->dot(PSVR_CAST2(y[i]), local_only);
    }
    SAMRAI_MPI::sumReduction(val, nv);
    IBTK_TIMER_STOP(t_vec_m_t_dot);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecScale_SAMRAI(Vec x, PetscScalar alpha)
{
    IBTK_TIMER_START(t_vec_scale);
    PSVR_CHECK1(x);
    static const bool interior_only = false;
    PSVR_CAST2(x)->scale(alpha, PSVR_CAST2(x), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_scale);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecCopy_SAMRAI(Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_copy);
    PSVR_CHECK2(x, y);
    static const bool interior_only = false;
    PSVR_CAST2(y)->copyVector(PSVR_CAST2(x), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_copy);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecSet_SAMRAI(Vec x, PetscScalar alpha)
{
    IBTK_TIMER_START(t_vec_set);
    PSVR_CHECK1(x);
    static const bool interior_only = false;
    PSVR_CAST2(x)->setToScalar(alpha, interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_set);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecSwap_SAMRAI(Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_swap);
    PSVR_CHECK2(x, y);
    PSVR_CAST2(x)->swapVectors(PSVR_CAST2(y));
    int ierr;
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_swap);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecAXPY_SAMRAI(Vec y, PetscScalar alpha, Vec x)
{
    IBTK_TIMER_START(t_vec_axpy);
    PSVR_CHECK2(x, y);
    static const bool interior_only = false;
    if (MathUtilities<double>::equalEps(alpha, 1.0))
    {
        PSVR_CAST2(y)->add(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else if (MathUtilities<double>::equalEps(alpha, -1.0))
    {
        PSVR_CAST2(y)->subtract(PSVR_CAST2(y), PSVR_CAST2(x), interior_only);
    }
    else
    {
        PSVR_CAST2(y)->axpy(alpha, PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_axpy);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecAXPBY_SAMRAI(Vec y, PetscScalar alpha, PetscScalar beta, Vec x)
{
    IBTK_TIMER_START(t_vec_axpby);
    PSVR_CHECK2(x, y);
    static const bool interior_only = false;
    if (MathUtilities<double>::equalEps(alpha, 1.0) && MathUtilities<double>::equalEps(beta, 1.0))
    {
        PSVR_CAST2(y)->add(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else if (MathUtilities<double>::equalEps(beta, 1.0))
    {
        PSVR_CAST2(y)->axpy(alpha, PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else if (MathUtilities<double>::equalEps(alpha, 1.0))
    {
        PSVR_CAST2(y)->axpy(beta, PSVR_CAST2(y), PSVR_CAST2(x), interior_only);
    }
    else
    {
        PSVR_CAST2(y)->linearSum(alpha, PSVR_CAST2(x), beta, PSVR_CAST2(y), interior_only);
    }
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_axpby);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMAXPY_SAMRAI(Vec y, PetscInt nv, const PetscScalar* alpha, Vec* x)
{
    IBTK_TIMER_START(t_vec_maxpy);
    PSVR_CHECK1(y);
    PSVR_CHECKN(x, nv);
    static const bool interior_only = false;
    for (PetscInt i = 0; i < nv; ++i)
    {
        if (MathUtilities<double>::equalEps(alpha[i], 1.0))
        {
            PSVR_CAST2(y)->add(PSVR_CAST2(x[i]), PSVR_CAST2(y), interior_only);
        }
        else if (MathUtilities<double>::equalEps(alpha[i], -1.0))
        {
            PSVR_CAST2(y)->subtract(PSVR_CAST2(y), PSVR_CAST2(x[i]), interior_only);
        }
        else
        {
            PSVR_CAST2(y)->axpy(alpha[i], PSVR_CAST2(x[i]), PSVR_CAST2(y), interior_only);
        }
    }
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_maxpy);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecAYPX_SAMRAI(Vec y, const PetscScalar alpha, Vec x)
{
    IBTK_TIMER_START(t_vec_aypx);
    PSVR_CHECK2(x, y);
    static const bool interior_only = false;
    if (MathUtilities<double>::equalEps(alpha, 1.0))
    {
        PSVR_CAST2(y)->add(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else if (MathUtilities<double>::equalEps(alpha, -1.0))
    {
        PSVR_CAST2(y)->subtract(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else
    {
        PSVR_CAST2(y)->axpy(alpha, PSVR_CAST2(y), PSVR_CAST2(x), interior_only);
    }
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_aypx);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecWAXPY_SAMRAI(Vec w, PetscScalar alpha, Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_waxpy);
    PSVR_CHECK3(w, x, y);
    static const bool interior_only = false;
    if (MathUtilities<double>::equalEps(alpha, 1.0))
    {
        PSVR_CAST2(w)->add(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else if (MathUtilities<double>::equalEps(alpha, -1.0))
    {
        PSVR_CAST2(w)->subtract(PSVR_CAST2(y), PSVR_CAST2(x), interior_only);
    }
    else
    {
        PSVR_CAST2(w)->axpy(alpha, PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_waxpy);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecAXPBYPCZ_SAMRAI(Vec z, PetscScalar alpha, PetscScalar beta, PetscScalar gamma, Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_axpbypcz);
    PSVR_CHECK3(x, y, z);
    static const bool interior_only = false;
    PSVR_CAST2(z)->linearSum(alpha, PSVR_CAST2(x), gamma, PSVR_CAST2(z), interior_only);
    PSVR_CAST2(z)->axpy(beta, PSVR_CAST2(y), PSVR_CAST2(z), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(z));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_axpbypcz);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecPointwiseMult_SAMRAI(Vec w, Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_pointwise_mult);
    PSVR_CHECK3(w, x, y);
    static const bool interior_only = false;
    PSVR_CAST2(w)->multiply(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_pointwise_mult);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecPointwiseDivide_SAMRAI(Vec w, Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_pointwise_divide);
    PSVR_CHECK3(w, x, y);
    static const bool interior_only = false;
    PSVR_CAST2(w)->divide(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_pointwise_divide);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecGetSize_SAMRAI(Vec v, PetscInt* size)
{
    IBTK_TIMER_START(t_vec_get_size);
    PSVR_CHECK1(v);
    *size = v->map->N;
    IBTK_TIMER_STOP(t_vec_get_size);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecGetLocalSize_SAMRAI(Vec v, PetscInt* size)
{
    IBTK_TIMER_START(t_vec_get_local_size);
    PSVR_CHECK1(v);
    *size = v->map->n;
    IBTK_TIMER_STOP(t_vec_get_local_size);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMax_SAMRAI(Vec x, PetscInt* p, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_max);
    PSVR_CHECK1(x);
    *p = -1;
    *val = PSVR_CAST2(x)->max();
    IBTK_TIMER_STOP(t_vec_max);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMin_SAMRAI(Vec x, PetscInt* p, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_min);
    PSVR_CHECK1(x);
    *p = -1;
    *val = PSVR_CAST2(x)->min();
    IBTK_TIMER_STOP(t_vec_min);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecSetRandom_SAMRAI(Vec x, PetscRandom rctx)
{
    IBTK_TIMER_START(t_vec_set_random);
    PSVR_CHECK1(x);
    PetscScalar lo, hi;
    int ierr;
    ierr = PetscRandomGetInterval(rctx, &lo, &hi);
    CHKERRQ(ierr);
    PSVR_CAST2(x)->setRandomValues(hi - lo, lo);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_set_random);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDot_local_SAMRAI(Vec x, Vec y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_dot_local);
    PSVR_CHECK2(x, y);
    static const bool local_only = true;
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y), local_only);
    IBTK_TIMER_STOP(t_vec_dot_local);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecTDot_local_SAMRAI(Vec x, Vec y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_t_dot_local);
    PSVR_CHECK2(x, y);
    static const bool local_only = true;
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y), local_only);
    IBTK_TIMER_STOP(t_vec_t_dot_local);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecNorm_local_SAMRAI(Vec x, NormType type, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_norm_local);
    PSVR_CHECK1(x);
    static const bool local_only = true;
    if (type == NORM_1)
    {
        *val = NormOps::L1Norm(PSVR_CAST2(x), local_only);
    }
    else if (type == NORM_2)
    {
        *val = NormOps::L2Norm(PSVR_CAST2(x), local_only);
    }
    else if (type == NORM_INFINITY)
    {
        *val = NormOps::maxNorm(PSVR_CAST2(x), local_only);
    }
    else if (type == NORM_1_AND_2)
    {
        val[0] = NormOps::L1Norm(PSVR_CAST2(x), local_only);
        val[1] = NormOps::L2Norm(PSVR_CAST2(x), local_only);
    }
    else
    {
        TBOX_ERROR("PETScSAMRAIVectorReal::norm()\n"
                   << "  vector norm type "
                   << static_cast<int>(type)
                   << " unsupported"
                   << std::endl);
    }
    IBTK_TIMER_STOP(t_vec_norm_local);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMDot_local_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_m_dot_local);
    PSVR_CHECK1(x);
    PSVR_CHECKN(y, nv);
    static const bool local_only = true;
    for (PetscInt i = 0; i < nv; ++i)
    {
        val[i] = PSVR_CAST2(x)->dot(PSVR_CAST2(y[i]), local_only);
    }
    IBTK_TIMER_STOP(t_vec_m_dot_local);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMTDot_local_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_m_t_dot_local);
    PSVR_CHECK1(x);
    PSVR_CHECKN(y, nv);
    static const bool local_only = true;
    for (PetscInt i = 0; i < nv; ++i)
    {
        val[i] = PSVR_CAST2(x)->dot(PSVR_CAST2(y[i]), local_only);
    }
    IBTK_TIMER_STOP(t_vec_m_t_dot_local);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMaxPointwiseDivide_SAMRAI(Vec x, Vec y, PetscScalar* max)
{
    IBTK_TIMER_START(t_vec_max_pointwise_divide);
    PSVR_CHECK2(x, y);
    *max = PSVR_CAST2(x)->maxPointwiseDivide(PSVR_CAST2(y));
    IBTK_TIMER_STOP(t_vec_max_pointwise_divide);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDotNorm2_SAMRAI(Vec s, Vec t, PetscScalar* dp, PetscScalar* nm)
{
    IBTK_TIMER_START(t_vec_dot_norm2);
    PSVR_CHECK2(s, t);
    *dp = PSVR_CAST2(s)->dot(PSVR_CAST2(t));
    *nm = PSVR_CAST2(t)->dot(PSVR_CAST2(t));
    IBTK_TIMER_STOP(t_vec_dot_norm2);
    PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
