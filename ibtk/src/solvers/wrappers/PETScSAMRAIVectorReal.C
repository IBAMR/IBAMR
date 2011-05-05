// Filename: PETScSAMRAIVectorReal.C
// Created on 10 Nov 2004 by Boyce Griffith
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

#include "PETScSAMRAIVectorReal.h"

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
#include <ibtk/NormOps.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <tbox/Pointer.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <ostream>

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

// Static functions for linkage with PETSc solver package routines.  These
// functions are intended to match those in the PETSc _VecOps structure.

#define PSVR_CAST1(v) (static_cast<PETScSAMRAIVectorReal*>(v->data))
#define PSVR_CAST2(v) (PETScSAMRAIVectorReal::getSAMRAIVector(v))

PetscErrorCode
VecDot_SAMRAI(
    Vec x,
    Vec y,
    PetscScalar* val)
{
    t_vec_dot->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y));
    t_vec_dot->stop();
    PetscFunctionReturn(0);
}// VecDot

PetscErrorCode
VecMDot_SAMRAI(
    Vec x,
    PetscInt nv,
    const Vec* y,
    PetscScalar* val)
{
    t_vec_m_dot->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    for (PetscInt i = 0; i < nv; ++i)
    {
        TBOX_ASSERT(y[i] != static_cast<Vec>(NULL));
    }
#endif
    for (PetscInt i = 0; i < nv; ++i)
    {
        val[i] = PSVR_CAST2(x)->dot(PSVR_CAST2(y[i]));
    }
    t_vec_m_dot->stop();
    PetscFunctionReturn(0);
}// VecMDot

PetscErrorCode
VecNorm_SAMRAI(
    Vec x,
    NormType type,
    PetscScalar* val)
{
    t_vec_norm->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
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
        val[0] = NormOps::L1Norm(PSVR_CAST2(x));
        val[1] = NormOps::L2Norm(PSVR_CAST2(x));
    }
    else
    {
        TBOX_ERROR("PETScSAMRAIVectorReal::norm()\n" <<
                   "  vector norm type " << int(type) << " unsupported" << std::endl);
    }
    t_vec_norm->stop();
    PetscFunctionReturn(0);
}// VecNorm

PetscErrorCode
VecTDot_SAMRAI(
    Vec x,
    Vec y,
    PetscScalar* val)
{
    t_vec_t_dot->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y));
    t_vec_t_dot->stop();
    PetscFunctionReturn(0);
}// VecTDot

PetscErrorCode
VecMTDot_SAMRAI(
    Vec x,
    PetscInt nv,
    const Vec* y,
    PetscScalar* val)
{
    t_vec_m_t_dot->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    for (PetscInt i = 0; i < nv; ++i)
    {
        TBOX_ASSERT(y[i] != static_cast<Vec>(NULL));
    }
#endif
    for (PetscInt i = 0; i < nv; ++i)
    {
        val[i] = PSVR_CAST2(x)->dot(PSVR_CAST2(y[i]));
    }
    t_vec_m_t_dot->stop();
    PetscFunctionReturn(0);
}// VecMTDot

PetscErrorCode
VecScale_SAMRAI(
    Vec x,
    PetscScalar alpha)
{
    t_vec_scale->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
    static const bool interior_only = false;
    PSVR_CAST2(x)->scale(alpha,PSVR_CAST2(x), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); IBTK_CHKERRQ(ierr);
    t_vec_scale->stop();
    PetscFunctionReturn(0);
}// VecScale

PetscErrorCode
VecCopy_SAMRAI(
    Vec x,
    Vec y)
{
    t_vec_copy->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
    static const bool interior_only = false;
    PSVR_CAST2(y)->copyVector(PSVR_CAST2(x), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);
    t_vec_copy->stop();
    PetscFunctionReturn(0);
}// VecCopy

PetscErrorCode
VecSet_SAMRAI(
    Vec x,
    PetscScalar alpha)
{
    t_vec_set->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
    static const bool interior_only = false;
    PSVR_CAST2(x)->setToScalar(alpha, interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); IBTK_CHKERRQ(ierr);
    t_vec_set->stop();
    PetscFunctionReturn(0);
}// VecSet

PetscErrorCode
VecSwap_SAMRAI(
    Vec x,
    Vec y)
{
    t_vec_swap->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
    PSVR_CAST2(x)->swapVectors(PSVR_CAST2(y));
    int ierr;
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);
    t_vec_swap->stop();
    PetscFunctionReturn(0);
}// VecSwap

PetscErrorCode
VecAXPY_SAMRAI(
    Vec y,
    PetscScalar alpha,
    Vec x)
{
    t_vec_axpy->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
    static const bool interior_only = false;
    if (MathUtilities<double>::equalEps(alpha,1.0))
    {
        PSVR_CAST2(y)->add(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else if (MathUtilities<double>::equalEps(alpha,-1.0))
    {
        PSVR_CAST2(y)->subtract(PSVR_CAST2(y), PSVR_CAST2(x), interior_only);
    }
    else
    {
        PSVR_CAST2(y)->axpy(alpha, PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);
    t_vec_axpy->stop();
    PetscFunctionReturn(0);
}// VecAXPY

PetscErrorCode
VecAXPBY_SAMRAI(
    Vec y,
    PetscScalar alpha,
    PetscScalar beta,
    Vec x)
{
    t_vec_axpby->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
    static const bool interior_only = false;
    if (MathUtilities<double>::equalEps(alpha,1.0) && MathUtilities<double>::equalEps(beta,1.0))
    {
        PSVR_CAST2(y)->add(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else if (MathUtilities<double>::equalEps(beta,1.0))
    {
        PSVR_CAST2(y)->axpy(alpha, PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else if (MathUtilities<double>::equalEps(alpha,1.0))
    {
        PSVR_CAST2(y)->axpy(beta, PSVR_CAST2(y), PSVR_CAST2(x), interior_only);
    }
    else
    {
        PSVR_CAST2(y)->linearSum(alpha, PSVR_CAST2(x), beta, PSVR_CAST2(y), interior_only);
    }
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);
    t_vec_axpby->stop();
    PetscFunctionReturn(0);
}// VecAXPBY

PetscErrorCode
VecMAXPY_SAMRAI(
    Vec y,
    PetscInt nv,
    const PetscScalar* alpha,
    Vec* x)
{
    t_vec_maxpy->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
    for (PetscInt i = 0; i < nv; ++i)
    {
        TBOX_ASSERT(x[i] != static_cast<Vec>(NULL));
    }
#endif
    static const bool interior_only = false;
    for (PetscInt i = 0; i < nv; ++i)
    {
        if (MathUtilities<double>::equalEps(alpha[i],1.0))
        {
            PSVR_CAST2(y)->add(PSVR_CAST2(x[i]), PSVR_CAST2(y), interior_only);
        }
        else if (MathUtilities<double>::equalEps(alpha[i],-1.0))
        {
            PSVR_CAST2(y)->subtract(PSVR_CAST2(y), PSVR_CAST2(x[i]), interior_only);
        }
        else
        {
            PSVR_CAST2(y)->axpy(alpha[i], PSVR_CAST2(x[i]), PSVR_CAST2(y), interior_only);
        }
    }
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);
    t_vec_maxpy->stop();
    PetscFunctionReturn(0);
}// VecMAXPY

PetscErrorCode
VecAYPX_SAMRAI(
    Vec y,
    const PetscScalar alpha,
    Vec x)
{
    t_vec_aypx->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
    static const bool interior_only = false;
    if (MathUtilities<double>::equalEps(alpha,1.0))
    {
        PSVR_CAST2(y)->add(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else if (MathUtilities<double>::equalEps(alpha,-1.0))
    {
        PSVR_CAST2(y)->subtract(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else
    {
        PSVR_CAST2(y)->axpy(alpha, PSVR_CAST2(y), PSVR_CAST2(x), interior_only);
    }
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);
    t_vec_aypx->stop();
    PetscFunctionReturn(0);
}// VecAYPX

PetscErrorCode
VecWAXPY_SAMRAI(
    Vec w,
    PetscScalar alpha,
    Vec x,
    Vec y)
{
    t_vec_waxpy->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
    TBOX_ASSERT(w != static_cast<Vec>(NULL));
#endif
    static const bool interior_only = false;
    if (MathUtilities<double>::equalEps(alpha,1.0))
    {
        PSVR_CAST2(w)->add(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    else if (MathUtilities<double>::equalEps(alpha,-1.0))
    {
        PSVR_CAST2(w)->subtract(PSVR_CAST2(y), PSVR_CAST2(x), interior_only);
    }
    else
    {
        PSVR_CAST2(w)->axpy(alpha, PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    }
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w)); IBTK_CHKERRQ(ierr);
    t_vec_waxpy->stop();
    PetscFunctionReturn(0);
}// VecWAXPY

PetscErrorCode
VecAXPBYPCZ_SAMRAI(
    Vec z,
    PetscScalar alpha,
    PetscScalar beta,
    PetscScalar gamma,
    Vec x,
    Vec y)
{
    t_vec_axpbypcz->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
    TBOX_ASSERT(z != static_cast<Vec>(NULL));
#endif
    static const bool interior_only = false;
    PSVR_CAST2(z)->linearSum(alpha, PSVR_CAST2(x), gamma, PSVR_CAST2(z), interior_only);
    PSVR_CAST2(z)->axpy(beta, PSVR_CAST2(y), PSVR_CAST2(z), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(z)); IBTK_CHKERRQ(ierr);
    t_vec_axpbypcz->stop();
    PetscFunctionReturn(0);
}// VecAXPBYPCZ

PetscErrorCode
VecPointwiseMult_SAMRAI(
    Vec w,
    Vec x,
    Vec y)
{
    t_vec_pointwise_mult->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
    TBOX_ASSERT(w != static_cast<Vec>(NULL));
#endif
    static const bool interior_only = false;
    PSVR_CAST2(w)->multiply(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w)); IBTK_CHKERRQ(ierr);
    t_vec_pointwise_mult->stop();
    PetscFunctionReturn(0);
}// VecPointwiseMult

PetscErrorCode
VecPointwiseDivide_SAMRAI(
    Vec w,
    Vec x,
    Vec y)
{
    t_vec_pointwise_divide->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
    TBOX_ASSERT(w != static_cast<Vec>(NULL));
#endif
    static const bool interior_only = false;
    PSVR_CAST2(w)->divide(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w)); IBTK_CHKERRQ(ierr);
    t_vec_pointwise_divide->stop();
    PetscFunctionReturn(0);
}// VecPointwiseDivide

PetscErrorCode
VecGetSize_SAMRAI(
    Vec v,
    PetscInt* n)
{
    t_vec_get_size->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(v != static_cast<Vec>(NULL));
#endif
    *n = 0;
    t_vec_get_size->stop();
    PetscFunctionReturn(0);
}// VecGetSize

PetscErrorCode
VecGetLocalSize_SAMRAI(
    Vec v,
    PetscInt* n)
{
    t_vec_get_local_size->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(v != static_cast<Vec>(NULL));
#endif
    *n = 0;
    t_vec_get_local_size->stop();
    PetscFunctionReturn(0);
}// VecGetLocalSize

PetscErrorCode
VecMax_SAMRAI(
    Vec x,
    PetscInt* p,
    PetscScalar* val)
{
    t_vec_max->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
    *p = -1;
    *val = PSVR_CAST2(x)->max();
    t_vec_max->stop();
    PetscFunctionReturn(0);
}// VecMax

PetscErrorCode
VecMin_SAMRAI(
    Vec x,
    PetscInt* p,
    PetscScalar* val)
{
    t_vec_min->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
    *p = -1;
    *val = PSVR_CAST2(x)->min();
    t_vec_min->stop();
    PetscFunctionReturn(0);
}// VecMin

PetscErrorCode
VecSetRandom_SAMRAI(
    Vec x,
    PetscRandom rctx)
{
    t_vec_set_random->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
    PetscScalar lo, hi;
    int ierr;
    ierr = PetscRandomGetInterval(rctx, &lo, &hi); IBTK_CHKERRQ(ierr);
    PSVR_CAST2(x)->setRandomValues(hi-lo, lo);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); IBTK_CHKERRQ(ierr);
    t_vec_set_random->stop();
    PetscFunctionReturn(0);
}// VecSetRandom

PetscErrorCode
VecDot_local_SAMRAI(
    Vec x,
    Vec y,
    PetscScalar* val)
{
    t_vec_dot_local->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
    static const bool local_only = true;
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y), local_only);
    t_vec_dot_local->stop();
    PetscFunctionReturn(0);
}// VecDot_local

PetscErrorCode
VecTDot_local_SAMRAI(
    Vec x,
    Vec y,
    PetscScalar* val)
{
    t_vec_t_dot_local->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
    static const bool local_only = true;
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y), local_only);
    t_vec_t_dot_local->stop();
    PetscFunctionReturn(0);
}// VecTDot_local

PetscErrorCode
VecNorm_local_SAMRAI(
    Vec x,
    NormType type,
    PetscScalar* val)
{
    t_vec_norm_local->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
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
        TBOX_ERROR("PETScSAMRAIVectorReal::norm()\n" <<
                   "  vector norm type " << int(type) << " unsupported" << std::endl);
    }
    t_vec_norm_local->stop();
    PetscFunctionReturn(0);
}// VecNorm_local

PetscErrorCode
VecMDot_local_SAMRAI(
    Vec x,
    PetscInt nv,
    const Vec* y,
    PetscScalar* val)
{
    t_vec_m_dot_local->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    for (PetscInt i = 0; i < nv; ++i)
    {
        TBOX_ASSERT(y[i] != static_cast<Vec>(NULL));
    }
#endif
    static const bool local_only = true;
    for (PetscInt i = 0; i < nv; ++i)
    {
        val[i] = PSVR_CAST2(x)->dot(PSVR_CAST2(y[i]), local_only);
    }
    t_vec_m_dot_local->stop();
    PetscFunctionReturn(0);
}// VecMDot_local

PetscErrorCode
VecMTDot_local_SAMRAI(
    Vec x,
    PetscInt nv,
    const Vec* y,
    PetscScalar* val)
{
    t_vec_m_t_dot_local->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    for (PetscInt i = 0; i < nv; ++i)
    {
        TBOX_ASSERT(y[i] != static_cast<Vec>(NULL));
    }
#endif
    static const bool local_only = true;
    for (PetscInt i = 0; i < nv; ++i)
    {
        val[i] = PSVR_CAST2(x)->dot(PSVR_CAST2(y[i]), local_only);
    }
    t_vec_m_t_dot_local->stop();
    PetscFunctionReturn(0);
}// VecMTDot_local

PetscErrorCode
VecMaxPointwiseDivide_SAMRAI(
    Vec x,
    Vec y,
    PetscScalar* max)
{
    t_vec_max_pointwise_divide->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x != static_cast<Vec>(NULL));
    TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
    *max = PSVR_CAST2(x)->maxPointwiseDivide(PSVR_CAST2(y));
    t_vec_max_pointwise_divide->stop();
    PetscFunctionReturn(0);
}// VecMaxPointwiseDivide

PetscErrorCode
VecDotNorm2_SAMRAI(
    Vec s,
    Vec t,
    PetscScalar* dp,
    PetscScalar* nm)
{
    t_vec_dot_norm2->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(s != static_cast<Vec>(NULL));
    TBOX_ASSERT(t != static_cast<Vec>(NULL));
#endif
    *dp = PSVR_CAST2(s)->dot(PSVR_CAST2(t));
    *nm = PSVR_CAST2(t)->dot(PSVR_CAST2(t));
    t_vec_dot_norm2->stop();
    PetscFunctionReturn(0);
}// VecDotNorm2_SAMRAI
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

PETScSAMRAIVectorReal::PETScSAMRAIVectorReal(
    Pointer<SAMRAIVectorReal<NDIM,PetscScalar> > samrai_vector,
    bool vector_created_via_duplicate,
    MPI_Comm comm)
    : d_samrai_vector(samrai_vector),
      d_vector_created_via_duplicate(vector_created_via_duplicate)
{
    // Setup Timers.
    IBTK_DO_ONCE(
        t_vec_duplicate            = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDuplicate()");
        t_vec_dot                  = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDot()");
        t_vec_m_dot                = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMDot()");
        t_vec_norm                 = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecNorm()");
        t_vec_t_dot                = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecTDot()");
        t_vec_m_t_dot              = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMTDot()");
        t_vec_scale                = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecScale()");
        t_vec_copy                 = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecCopy()");
        t_vec_set                  = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecSet()");
        t_vec_swap                 = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecSwap()");
        t_vec_axpy                 = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAXPY()");
        t_vec_axpby                = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAXPBY()");
        t_vec_maxpy                = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMAXPY()");
        t_vec_aypx                 = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAYPX()");
        t_vec_waxpy                = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecWAXPY()");
        t_vec_axpbypcz             = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAXPBYPCZ()");
        t_vec_pointwise_mult       = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecPointwiseMult()");
        t_vec_pointwise_divide     = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecPointwiseDivide()");
        t_vec_get_size             = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecGetSize()");
        t_vec_get_local_size       = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecGetLocalSize()");
        t_vec_max                  = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMax()");
        t_vec_min                  = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMin()");
        t_vec_set_random           = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecSetRandom()");
        t_vec_destroy              = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDestroy()");
        t_vec_dot_local            = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDot_local()");
        t_vec_t_dot_local          = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecTDot_local()");
        t_vec_norm_local           = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecNorm_local()");
        t_vec_m_dot_local          = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMDot_local()");
        t_vec_m_t_dot_local        = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMTDot_local()");
        t_vec_max_pointwise_divide = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMaxPointwiseDivide()");
        t_vec_dot_norm2            = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDotNorm2()");
                 );

    int ierr;
    ierr = VecCreate(comm, &d_petsc_vector); IBTK_CHKERRQ(ierr);

    // Set PETSc vector data to this abstract vector object
    d_petsc_vector->petscnative                  = PETSC_FALSE;
    d_petsc_vector->map->n                       = 0;
    d_petsc_vector->map->N                       = 0;
    d_petsc_vector->map->bs                      = 1;

    // Assign vector operations to PETSc vector object.
    d_petsc_vector->ops->duplicate               = PETScSAMRAIVectorReal::VecDuplicate_SAMRAI;
    d_petsc_vector->ops->duplicatevecs           = VecDuplicateVecs_Default;
    d_petsc_vector->ops->destroyvecs             = VecDestroyVecs_Default;
    d_petsc_vector->ops->dot                     = VecDot_SAMRAI;
    d_petsc_vector->ops->mdot                    = VecMDot_SAMRAI;
    d_petsc_vector->ops->norm                    = VecNorm_SAMRAI;
    d_petsc_vector->ops->tdot                    = VecTDot_SAMRAI;
    d_petsc_vector->ops->mtdot                   = VecMTDot_SAMRAI;
    d_petsc_vector->ops->scale                   = VecScale_SAMRAI;
    d_petsc_vector->ops->copy                    = VecCopy_SAMRAI;
    d_petsc_vector->ops->set                     = VecSet_SAMRAI;
    d_petsc_vector->ops->swap                    = VecSwap_SAMRAI;
    d_petsc_vector->ops->axpy                    = VecAXPY_SAMRAI;
    d_petsc_vector->ops->axpby                   = VecAXPBY_SAMRAI;
    d_petsc_vector->ops->maxpy                   = VecMAXPY_SAMRAI;
    d_petsc_vector->ops->aypx                    = VecAYPX_SAMRAI;
    d_petsc_vector->ops->waxpy                   = VecWAXPY_SAMRAI;
    d_petsc_vector->ops->axpbypcz                = VecAXPBYPCZ_SAMRAI;
    d_petsc_vector->ops->pointwisemult           = VecPointwiseMult_SAMRAI;
    d_petsc_vector->ops->pointwisedivide         = VecPointwiseDivide_SAMRAI;
    d_petsc_vector->ops->setvalues               = PETSC_NULL;
    d_petsc_vector->ops->assemblybegin           = PETSC_NULL;
    d_petsc_vector->ops->assemblyend             = PETSC_NULL;
    d_petsc_vector->ops->getarray                = PETSC_NULL;
    d_petsc_vector->ops->getsize                 = VecGetSize_SAMRAI;
    d_petsc_vector->ops->getlocalsize            = VecGetLocalSize_SAMRAI;
    d_petsc_vector->ops->restorearray            = PETSC_NULL;
    d_petsc_vector->ops->max                     = VecMax_SAMRAI;
    d_petsc_vector->ops->min                     = VecMin_SAMRAI;
    d_petsc_vector->ops->setrandom               = VecSetRandom_SAMRAI;
    d_petsc_vector->ops->setoption               = PETSC_NULL;
    d_petsc_vector->ops->setvaluesblocked        = PETSC_NULL;
    d_petsc_vector->ops->destroy                 = PETScSAMRAIVectorReal::VecDestroy_SAMRAI;
    d_petsc_vector->ops->view                    = PETSC_NULL;
    d_petsc_vector->ops->placearray              = PETSC_NULL;
    d_petsc_vector->ops->replacearray            = PETSC_NULL;
    d_petsc_vector->ops->dot_local               = VecDot_local_SAMRAI;
    d_petsc_vector->ops->tdot_local              = VecTDot_local_SAMRAI;
    d_petsc_vector->ops->norm_local              = VecNorm_local_SAMRAI;
    d_petsc_vector->ops->mdot_local              = VecMDot_local_SAMRAI;
    d_petsc_vector->ops->mtdot_local             = VecMTDot_local_SAMRAI;
    d_petsc_vector->ops->loadintovector          = PETSC_NULL;
    d_petsc_vector->ops->loadintovectornative    = PETSC_NULL;
    d_petsc_vector->ops->reciprocal              = PETSC_NULL;
    d_petsc_vector->ops->viewnative              = PETSC_NULL;
    d_petsc_vector->ops->conjugate               = PETSC_NULL;
    d_petsc_vector->ops->setlocaltoglobalmapping = PETSC_NULL;
    d_petsc_vector->ops->setvalueslocal          = PETSC_NULL;
    d_petsc_vector->ops->resetarray              = PETSC_NULL;
    d_petsc_vector->ops->setfromoptions          = PETSC_NULL;
    d_petsc_vector->ops->maxpointwisedivide      = VecMaxPointwiseDivide_SAMRAI;
    d_petsc_vector->ops->load                    = PETSC_NULL;
    d_petsc_vector->ops->pointwisemax            = PETSC_NULL;
    d_petsc_vector->ops->pointwisemaxabs         = PETSC_NULL;
    d_petsc_vector->ops->pointwisemin            = PETSC_NULL;
    d_petsc_vector->ops->getvalues               = PETSC_NULL;
    d_petsc_vector->ops->sqrt                    = PETSC_NULL;
    d_petsc_vector->ops->abs                     = PETSC_NULL;
    d_petsc_vector->ops->exp                     = PETSC_NULL;
    d_petsc_vector->ops->log                     = PETSC_NULL;
    d_petsc_vector->ops->shift                   = PETSC_NULL;
    d_petsc_vector->ops->create                  = PETSC_NULL;
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1)
    d_petsc_vector->ops->dotnorm2                = VecDotNorm2_SAMRAI;
#endif

    ierr = PetscObjectChangeTypeName(reinterpret_cast<PetscObject>(d_petsc_vector),"Vec_SAMRAI");  IBTK_CHKERRQ(ierr);

    d_petsc_vector->data = this;
    return;
}// PETScSAMRAIVectorReal

PETScSAMRAIVectorReal::~PETScSAMRAIVectorReal()
{
    if (!d_vector_created_via_duplicate)
    {
        d_petsc_vector->ops->destroy = 0;
        int ierr = VecDestroy(d_petsc_vector); IBTK_CHKERRQ(ierr);
    }
    return;
}// ~PETScSAMRAIVectorReal

/////////////////////////////// PRIVATE //////////////////////////////////////

PetscErrorCode
PETScSAMRAIVectorReal::VecDuplicate_SAMRAI(
    Vec v,
    Vec* newv)
{
    t_vec_duplicate->start();
    PetscErrorCode ierr;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(v != static_cast<Vec>(NULL));
#endif
    Pointer<SAMRAIVectorReal<NDIM,PetscScalar> > samrai_vec = PSVR_CAST2(v)->cloneVector(PSVR_CAST2(v)->getName());
    samrai_vec->allocateVectorData();
    static const bool vector_created_via_duplicate = true;
    MPI_Comm comm;
    ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(v), &comm); IBTK_CHKERRQ(ierr);
    PETScSAMRAIVectorReal* new_psv = new PETScSAMRAIVectorReal(samrai_vec, vector_created_via_duplicate, comm);
    *newv = new_psv->d_petsc_vector;
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*newv)); IBTK_CHKERRQ(ierr);
    t_vec_duplicate->stop();
    PetscFunctionReturn(0);
}// VecDuplicate

PetscErrorCode
PETScSAMRAIVectorReal::VecDestroy_SAMRAI(
    Vec v)
{
    t_vec_destroy->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(v != static_cast<Vec>(NULL));
#endif
    if (PSVR_CAST1(v)->d_vector_created_via_duplicate)
    {
        PSVR_CAST2(v)->resetLevels(0, std::min(PSVR_CAST2(v)->getFinestLevelNumber(),
                                               PSVR_CAST2(v)->getPatchHierarchy()->getFinestLevelNumber()));
        PSVR_CAST2(v)->freeVectorComponents();
        PSVR_CAST2(v).setNull();
        destroyPETScVector(PSVR_CAST1(v)->d_petsc_vector);
    }
    t_vec_destroy->stop();
    PetscFunctionReturn(0);
}// VecDestroy

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

//////////////////////////////////////////////////////////////////////////////
