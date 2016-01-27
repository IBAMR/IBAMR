// Filename: PETScMultiVec.cpp
// Created on 13 Mar 2008 by Boyce Griffith
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
#include <stddef.h>
#include <algorithm>
#include <ostream>

#include "ibtk/PETScMultiVec.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "mpi.h"
#include "petsc/private/vecimpl.h" // IWYU pragma: keep
#include "petscerror.h"
#include "petscis.h"
#include "petscmath.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

struct Vec_MultiVec
{
    PetscInt n;           // the number of component vector objects
    Vec* array;           // the array  of component vector objects
    Vec* array_allocated; // if the array was allocated by PetscMalloc, this is its pointer
};

#undef __FUNCT__
#define __FUNCT__ "VecDuplicate_MultiVec"
PetscErrorCode
VecDuplicate_MultiVec(Vec v, Vec* newv)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(v);
#endif
    Vec_MultiVec* mv = static_cast<Vec_MultiVec*>(v->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mv);
#endif
    PetscErrorCode ierr;
    Vec* newvarray;
    ierr = PetscMalloc(mv->n * sizeof(Vec), &newvarray);
    CHKERRQ(ierr);
    for (PetscInt k = 0; k < mv->n; ++k)
    {
        ierr = VecDuplicate(mv->array[k], &newvarray[k]);
        CHKERRQ(ierr);
    }
    MPI_Comm comm;
    ierr = PetscObjectGetComm((PetscObject)v, &comm);
    CHKERRQ(ierr);
    ierr = VecCreateMultiVec(comm, mv->n, newvarray, newv);
    CHKERRQ(ierr);
    Vec_MultiVec* mnewv = static_cast<Vec_MultiVec*>((*newv)->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mnewv);
#endif
    mnewv->array_allocated = newvarray;
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*newv));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecDuplicate_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecDot_MultiVec"
PetscErrorCode
VecDot_MultiVec(Vec x, Vec y, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
    TBOX_ASSERT(y);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my);
    TBOX_ASSERT(mx->n == my->n);
#endif
    PetscErrorCode ierr;
    *val = 0.0;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
        PetscScalar component_val;
        ierr = VecDot(mx->array[k], my->array[k], &component_val);
        CHKERRQ(ierr);
        *val += component_val;
    }
    PetscFunctionReturn(0);
} // VecDot_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecMDot_MultiVec"
PetscErrorCode
VecMDot_MultiVec(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
    for (PetscInt i = 0; i < nv; ++i)
    {
        TBOX_ASSERT(y[i]);
    }
#endif
    PetscErrorCode ierr;
    for (PetscInt i = 0; i < nv; ++i)
    {
        ierr = VecDot_MultiVec(x, y[i], &val[i]);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} // VecMDot_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecNorm_MultiVec"
PetscErrorCode
VecNorm_MultiVec(Vec x, NormType type, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
#endif
    PetscErrorCode ierr;
    if (type == NORM_1)
    {
        *val = 0.0;
        for (PetscInt k = 0; k < mx->n; ++k)
        {
            PetscScalar component_val;
            ierr = VecNorm(mx->array[k], type, &component_val);
            CHKERRQ(ierr);
            *val += component_val;
        }
    }
    else if (type == NORM_2)
    {
        *val = 0.0;
        for (PetscInt k = 0; k < mx->n; ++k)
        {
            PetscScalar component_val;
            ierr = VecNorm(mx->array[k], type, &component_val);
            CHKERRQ(ierr);
            *val += component_val * component_val;
        }
        *val = sqrt(*val);
    }
    else if (type == NORM_INFINITY)
    {
        *val = PETSC_MIN_REAL;
        for (PetscInt k = 0; k < mx->n; ++k)
        {
            PetscScalar component_val;
            ierr = VecNorm(mx->array[k], type, &component_val);
            CHKERRQ(ierr);
            *val = std::max(*val, component_val);
        }
    }
    else if (type == NORM_1_AND_2)
    {
        ierr = VecNorm_MultiVec(x, NORM_1, &val[0]);
        CHKERRQ(ierr);
        ierr = VecNorm_MultiVec(x, NORM_2, &val[1]);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} // VecNorm_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecTDot_MultiVec"
PetscErrorCode
VecTDot_MultiVec(Vec x, Vec y, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
    TBOX_ASSERT(y);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my);
    TBOX_ASSERT(mx->n == my->n);
#endif
    PetscErrorCode ierr;
    *val = 0.0;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
        PetscScalar component_val;
        ierr = VecTDot(mx->array[k], my->array[k], &component_val);
        CHKERRQ(ierr);
        *val += component_val;
    }
    PetscFunctionReturn(0);
} // VecTDot_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecMTDot_MultiVec"
PetscErrorCode
VecMTDot_MultiVec(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
    for (PetscInt i = 0; i < nv; ++i)
    {
        TBOX_ASSERT(y[i]);
    }
#endif
    PetscErrorCode ierr;
    for (PetscInt i = 0; i < nv; ++i)
    {
        ierr = VecTDot_MultiVec(x, y[i], &val[i]);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} // VecMTDot_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecScale_MultiVec"
PetscErrorCode
VecScale_MultiVec(Vec x, PetscScalar alpha)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
        ierr = VecScale(mx->array[k], alpha);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecScale_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecCopy_MultiVec"
PetscErrorCode
VecCopy_MultiVec(Vec x, Vec y)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
    TBOX_ASSERT(y);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my);
    TBOX_ASSERT(mx->n == my->n);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
        ierr = VecCopy(mx->array[k], my->array[k]);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecCopy_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecSet_MultiVec"
PetscErrorCode
VecSet_MultiVec(Vec x, PetscScalar alpha)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
        ierr = VecSet(mx->array[k], alpha);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecSet_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecSwap_MultiVec"
PetscErrorCode
VecSwap_MultiVec(Vec x, Vec y)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
    TBOX_ASSERT(y);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my);
    TBOX_ASSERT(mx->n == my->n);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
        ierr = VecSwap(mx->array[k], my->array[k]);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecSwap_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecAXPY_MultiVec"
PetscErrorCode
VecAXPY_MultiVec(Vec y, PetscScalar alpha, Vec x)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(y);
    TBOX_ASSERT(x);
#endif
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(my);
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my->n == mx->n);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < my->n; ++k)
    {
        ierr = VecAXPY(my->array[k], alpha, mx->array[k]);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecAXPY_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecAXPBY_MultiVec"
PetscErrorCode
VecAXPBY_MultiVec(Vec y, PetscScalar alpha, PetscScalar beta, Vec x)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(y);
    TBOX_ASSERT(x);
#endif
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(my);
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my->n == mx->n);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < my->n; ++k)
    {
        ierr = VecAXPBY(my->array[k], alpha, beta, mx->array[k]);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecAXPBY_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecMAXPY_MultiVec"
PetscErrorCode
VecMAXPY_MultiVec(Vec y, PetscInt nv, const PetscScalar* alpha, Vec* x)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(y);
    for (PetscInt i = 0; i < nv; ++i)
    {
        TBOX_ASSERT(x[i]);
    }
#endif
    PetscErrorCode ierr;
    for (PetscInt i = 0; i < nv; ++i)
    {
        ierr = VecAXPY_MultiVec(y, alpha[i], x[i]);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecMAXPY_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecAYPX_MultiVec"
PetscErrorCode
VecAYPX_MultiVec(Vec y, const PetscScalar alpha, Vec x)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(y);
    TBOX_ASSERT(x);
#endif
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(my);
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my->n == mx->n);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < my->n; ++k)
    {
        ierr = VecAYPX(my->array[k], alpha, mx->array[k]);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecAYPX_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecWAXPY_MultiVec"
PetscErrorCode
VecWAXPY_MultiVec(Vec w, PetscScalar alpha, Vec x, Vec y)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(w);
    TBOX_ASSERT(x);
    TBOX_ASSERT(y);
#endif
    Vec_MultiVec* mw = static_cast<Vec_MultiVec*>(w->data);
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mw);
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my);
    TBOX_ASSERT(mw->n == mx->n);
    TBOX_ASSERT(mw->n == my->n);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < mw->n; ++k)
    {
        ierr = VecWAXPY(mw->array[k], alpha, mx->array[k], my->array[k]);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecWAXPY_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecAXPBYPCZ_MultiVec"
PetscErrorCode
VecAXPBYPCZ_MultiVec(Vec w, PetscScalar alpha, PetscScalar beta, PetscScalar gamma, Vec x, Vec y)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(w);
    TBOX_ASSERT(x);
    TBOX_ASSERT(y);
#endif
    Vec_MultiVec* mw = static_cast<Vec_MultiVec*>(w->data);
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mw);
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my);
    TBOX_ASSERT(mw->n == mx->n);
    TBOX_ASSERT(mw->n == my->n);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < mw->n; ++k)
    {
        ierr = VecAXPBYPCZ(mw->array[k], alpha, beta, gamma, mx->array[k], my->array[k]);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecAXPBYPCZ_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecPointwiseMult_MultiVec"
PetscErrorCode
VecPointwiseMult_MultiVec(Vec w, Vec x, Vec y)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(w);
    TBOX_ASSERT(x);
    TBOX_ASSERT(y);
#endif
    Vec_MultiVec* mw = static_cast<Vec_MultiVec*>(w->data);
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mw);
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my);
    TBOX_ASSERT(mw->n == mx->n);
    TBOX_ASSERT(mw->n == my->n);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < mw->n; ++k)
    {
        ierr = VecPointwiseMult(mw->array[k], mx->array[k], my->array[k]);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecPointwiseMult_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecPointwiseDivide_MultiVec"
PetscErrorCode
VecPointwiseDivide_MultiVec(Vec w, Vec x, Vec y)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(w);
    TBOX_ASSERT(x);
    TBOX_ASSERT(y);
#endif
    Vec_MultiVec* mw = static_cast<Vec_MultiVec*>(w->data);
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mw);
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my);
    TBOX_ASSERT(mw->n == mx->n);
    TBOX_ASSERT(mw->n == my->n);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < mw->n; ++k)
    {
        ierr = VecPointwiseDivide(mw->array[k], mx->array[k], my->array[k]);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecPointwiseDivide_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecGetSize_MultiVec"
PetscErrorCode
VecGetSize_MultiVec(Vec x, PetscInt* size)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
#endif
    *size = x->map->N;
    PetscFunctionReturn(0);
} // VecGetSize_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecGetLocalSize_MultiVec"
PetscErrorCode
VecGetLocalSize_MultiVec(Vec x, PetscInt* size)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
#endif
    *size = x->map->n;
    PetscFunctionReturn(0);
} // VecGetLocalSize_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecMax_MultiVec"
PetscErrorCode
VecMax_MultiVec(Vec x, PetscInt* p, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
#endif
    PetscErrorCode ierr;
    *p = -1;
    *val = PETSC_MIN_REAL;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
        PetscInt component_p;
        PetscScalar component_val;
        ierr = VecMax(mx->array[k], &component_p, &component_val);
        CHKERRQ(ierr);
        *val = std::max(*val, component_val);
    }
    PetscFunctionReturn(0);
} // VecMax_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecMin_MultiVec"
PetscErrorCode
VecMin_MultiVec(Vec x, PetscInt* p, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
#endif
    PetscErrorCode ierr;
    *p = -1;
    *val = PETSC_MAX_REAL;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
        PetscInt component_p;
        PetscScalar component_val;
        ierr = VecMin(mx->array[k], &component_p, &component_val);
        CHKERRQ(ierr);
        *val = std::min(*val, component_val);
    }
    PetscFunctionReturn(0);
} // VecMin_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecSetRandom_MultiVec"
PetscErrorCode
VecSetRandom_MultiVec(Vec x, PetscRandom rctx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
#endif
    PetscErrorCode ierr;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
        ierr = VecSetRandom(mx->array[k], rctx);
        CHKERRQ(ierr);
    }
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecSetRandom_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecDestroy_MultiVec"
PetscErrorCode
VecDestroy_MultiVec(Vec v)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(v);
#endif
    Vec_MultiVec* mv = static_cast<Vec_MultiVec*>(v->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mv);
#endif
    PetscErrorCode ierr;
    if (mv->array_allocated)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(mv->array_allocated == mv->array);
#endif
        for (PetscInt k = 0; k < mv->n; ++k)
        {
            ierr = VecDestroy(&mv->array_allocated[k]);
            CHKERRQ(ierr);
        }
        ierr = PetscFree(mv->array_allocated);
        CHKERRQ(ierr);
    }
    ierr = PetscFree(mv);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecDestroy_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecDot_local_MultiVec"
PetscErrorCode
VecDot_local_MultiVec(Vec x, Vec y, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
    TBOX_ASSERT(y);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my);
    TBOX_ASSERT(mx->n == my->n);
#endif
    PetscErrorCode ierr;
    *val = 0.0;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(mx->array[k]->ops->dot_local);
        TBOX_ASSERT(mx->array[k]->ops->dot_local == my->array[k]->ops->dot_local);
#endif
        PetscScalar component_val;
        ierr = (*mx->array[k]->ops->dot_local)(mx->array[k], my->array[k], &component_val);
        CHKERRQ(ierr);
        *val += component_val;
    }
    PetscFunctionReturn(0);
} // VecDot_local_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecTDot_local_MultiVec"
PetscErrorCode
VecTDot_local_MultiVec(Vec x, Vec y, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
    TBOX_ASSERT(y);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my);
    TBOX_ASSERT(mx->n == my->n);
#endif
    PetscErrorCode ierr;
    *val = 0.0;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(mx->array[k]->ops->tdot_local);
        TBOX_ASSERT(mx->array[k]->ops->tdot_local == my->array[k]->ops->dot_local);
#endif
        PetscScalar component_val;
        ierr = (*mx->array[k]->ops->tdot_local)(mx->array[k], my->array[k], &component_val);
        CHKERRQ(ierr);
        *val += component_val;
    }
    PetscFunctionReturn(0);
} // VecTDot_local_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecNorm_local_MultiVec"
PetscErrorCode
VecNorm_local_MultiVec(Vec x, NormType type, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
#endif
    PetscErrorCode ierr;
    if (type == NORM_1)
    {
        *val = 0.0;
        for (PetscInt k = 0; k < mx->n; ++k)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(mx->array[k]->ops->norm_local);
#endif
            PetscScalar component_val;
            ierr = (*mx->array[k]->ops->norm_local)(mx->array[k], type, &component_val);
            CHKERRQ(ierr);
            *val += component_val;
        }
    }
    else if (type == NORM_2)
    {
        *val = 0.0;
        for (PetscInt k = 0; k < mx->n; ++k)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(mx->array[k]->ops->norm_local);
#endif
            PetscScalar component_val;
            ierr = (*mx->array[k]->ops->norm_local)(mx->array[k], type, &component_val);
            CHKERRQ(ierr);
            *val += component_val * component_val;
        }
        *val = sqrt(*val);
    }
    else if (type == NORM_INFINITY)
    {
        *val = PETSC_MIN_REAL;
        for (PetscInt k = 0; k < mx->n; ++k)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(mx->array[k]->ops->norm_local);
#endif
            PetscScalar component_val;
            ierr = (*mx->array[k]->ops->norm_local)(mx->array[k], type, &component_val);
            CHKERRQ(ierr);
            *val = std::max(*val, component_val);
        }
    }
    else if (type == NORM_1_AND_2)
    {
        ierr = VecNorm_local_MultiVec(x, NORM_1, &val[0]);
        CHKERRQ(ierr);
        ierr = VecNorm_local_MultiVec(x, NORM_2, &val[1]);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} // VecNorm_local_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecMDot_local_MultiVec"
PetscErrorCode
VecMDot_local_MultiVec(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
    for (PetscInt i = 0; i < nv; ++i)
    {
        TBOX_ASSERT(y[i]);
    }
#endif
    PetscErrorCode ierr;
    for (PetscInt i = 0; i < nv; ++i)
    {
        ierr = VecDot_local_MultiVec(x, y[i], &val[i]);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} // VecMDot_local_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecMTDot_local_MultiVec"
PetscErrorCode
VecMTDot_local_MultiVec(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
    for (PetscInt i = 0; i < nv; ++i)
    {
        TBOX_ASSERT(y[i]);
    }
#endif
    PetscErrorCode ierr;
    for (PetscInt i = 0; i < nv; ++i)
    {
        ierr = VecTDot_local_MultiVec(x, y[i], &val[i]);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} // VecMTDot_local_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecMaxPointwiseDivide_MultiVec"
PetscErrorCode
VecMaxPointwiseDivide_MultiVec(Vec x, Vec y, PetscScalar* max)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(x);
    TBOX_ASSERT(y);
#endif
    Vec_MultiVec* mx = static_cast<Vec_MultiVec*>(x->data);
    Vec_MultiVec* my = static_cast<Vec_MultiVec*>(y->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mx);
    TBOX_ASSERT(my);
    TBOX_ASSERT(mx->n == my->n);
#endif
    PetscErrorCode ierr;
    *max = PETSC_MIN_REAL;
    for (PetscInt k = 0; k < mx->n; ++k)
    {
        PetscScalar component_max;
        ierr = VecMaxPointwiseDivide(mx->array[k], my->array[k], &component_max);
        CHKERRQ(ierr);
        *max = std::max(*max, component_max);
    }
    PetscFunctionReturn(0);
} // VecMaxPointwiseDivide_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecDotNorm2_MultiVec"
PetscErrorCode
VecDotNorm2_MultiVec(Vec s, Vec t, PetscScalar* dp, PetscScalar* nm)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(s);
    TBOX_ASSERT(t);
#endif
    PetscErrorCode ierr;
    ierr = VecDot_MultiVec(s, t, dp);
    CHKERRQ(ierr);
    ierr = VecDot_MultiVec(t, t, nm);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecDotNorm2_MultiVec

#undef __FUNCT__
#define __FUNCT__ "VecCreateMultiVec"
PetscErrorCode
VecCreateMultiVec(MPI_Comm comm, PetscInt n, Vec vv[], Vec* v)
{
    PetscErrorCode ierr;
    ierr = VecCreate(comm, v);
    CHKERRQ(ierr);

    // Assign vector operations to PETSc vector object.
    static struct _VecOps DvOps;
    IBTK_DO_ONCE(DvOps.duplicate = VecDuplicate_MultiVec; DvOps.duplicatevecs = VecDuplicateVecs_Default;
                 DvOps.destroyvecs = VecDestroyVecs_Default;
                 DvOps.dot = VecDot_MultiVec;
                 DvOps.mdot = VecMDot_MultiVec;
                 DvOps.norm = VecNorm_MultiVec;
                 DvOps.tdot = VecTDot_MultiVec;
                 DvOps.mtdot = VecMTDot_MultiVec;
                 DvOps.scale = VecScale_MultiVec;
                 DvOps.copy = VecCopy_MultiVec;
                 DvOps.set = VecSet_MultiVec;
                 DvOps.swap = VecSwap_MultiVec;
                 DvOps.axpy = VecAXPY_MultiVec;
                 DvOps.axpby = VecAXPBY_MultiVec;
                 DvOps.maxpy = VecMAXPY_MultiVec;
                 DvOps.aypx = VecAYPX_MultiVec;
                 DvOps.waxpy = VecWAXPY_MultiVec;
                 DvOps.axpbypcz = VecAXPBYPCZ_MultiVec;
                 DvOps.pointwisemult = VecPointwiseMult_MultiVec;
                 DvOps.pointwisedivide = VecPointwiseDivide_MultiVec;
                 DvOps.getsize = VecGetSize_MultiVec;
                 DvOps.getlocalsize = VecGetLocalSize_MultiVec;
                 DvOps.max = VecMax_MultiVec;
                 DvOps.min = VecMin_MultiVec;
                 DvOps.setrandom = VecSetRandom_MultiVec;
                 DvOps.destroy = VecDestroy_MultiVec;
                 DvOps.dot_local = VecDot_local_MultiVec;
                 DvOps.tdot_local = VecTDot_local_MultiVec;
                 DvOps.norm_local = VecNorm_local_MultiVec;
                 DvOps.mdot_local = VecMDot_local_MultiVec;
                 DvOps.mtdot_local = VecMTDot_local_MultiVec;
                 DvOps.maxpointwisedivide = VecMaxPointwiseDivide_MultiVec;
                 DvOps.dotnorm2 = VecDotNorm2_MultiVec;);
    ierr = PetscMemcpy((*v)->ops, &DvOps, sizeof(DvOps));
    CHKERRQ(ierr);

    // Set PETSc vector data.
    Vec_MultiVec* mv;
    ierr = PetscNew(&mv);
    CHKERRQ(ierr);
    mv->n = n;
    mv->array = vv;
    mv->array_allocated = NULL;
    (*v)->data = mv;
    (*v)->petscnative = PETSC_FALSE;
    (*v)->map->n = 0;
    (*v)->map->N = 0;
    for (PetscInt k = 0; k < mv->n; ++k)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(mv->array[k]);
#endif
        (*v)->map->n += mv->array[k]->map->n;
        (*v)->map->N += mv->array[k]->map->N;
    }
    (*v)->map->bs = 1; // NOTE: Here we are giving a bogus block size.

    // Set the PETSc vector type name.
    ierr = PetscObjectChangeTypeName(reinterpret_cast<PetscObject>(*v), "Vec_MultiVec");
    CHKERRQ(ierr);

    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*v));
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // VecCreateMultiVec

#undef __FUNCT__
#define __FUNCT__ "VecMultiVecGetNumberOfSubVecs"
PetscErrorCode
VecMultiVecGetNumberOfSubVecs(Vec v, PetscInt* n)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(v);
#endif
    Vec_MultiVec* mv = static_cast<Vec_MultiVec*>(v->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mv);
#endif
    *n = mv->n;
    PetscFunctionReturn(0);
} // VecMultiVecGetNumberOfSubVecs

#undef __FUNCT__
#define __FUNCT__ "VecMultiVecGetSubVecs"
PetscErrorCode
VecMultiVecGetSubVecs(Vec v, Vec* vv[])
{
#if !defined(NDEBUG)
    TBOX_ASSERT(v);
#endif
    Vec_MultiVec* mv = static_cast<Vec_MultiVec*>(v->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mv);
#endif
    *vv = mv->array;
    PetscFunctionReturn(0);
} // VecMultiVecGetSubVecs

#undef __FUNCT__
#define __FUNCT__ "VecMultiVecGetSubVec"
PetscErrorCode
VecMultiVecGetSubVec(Vec v, PetscInt idx, Vec* subv)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(v);
#endif
    Vec_MultiVec* mv = static_cast<Vec_MultiVec*>(v->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(mv);
    TBOX_ASSERT(0 <= idx && idx < mv->n);
#endif
    *subv = mv->array[idx];
    PetscFunctionReturn(0);
} // VecMultiVecGetSubVec

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
