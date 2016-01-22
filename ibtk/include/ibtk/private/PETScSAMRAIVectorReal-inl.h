// Filename: PETScSAMRAIVectorReal-inl.h
// Created on 27 Aug 2010 by Boyce Griffith
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

#ifndef included_PETScSAMRAIVectorReal_inl_h
#define included_PETScSAMRAIVectorReal_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "petsc/private/petscimpl.h" // IWYU pragma: keep
#include "petsc/private/vecimpl.h"   // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline Vec
PETScSAMRAIVectorReal::createPETScVector(
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> > samrai_vec,
    MPI_Comm comm)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(samrai_vec);
#endif
    static const bool vector_created_via_duplicate = false;
    PETScSAMRAIVectorReal* psv = new PETScSAMRAIVectorReal(samrai_vec, vector_created_via_duplicate, comm);
    return psv->d_petsc_vector;
} // createPETScVector

inline void
PETScSAMRAIVectorReal::destroyPETScVector(Vec petsc_vec)
{
    if (petsc_vec)
    {
        PETScSAMRAIVectorReal* psv = static_cast<PETScSAMRAIVectorReal*>(petsc_vec->data);
#if !defined(NDEBUG)
        TBOX_ASSERT(psv);
#endif
        delete psv;
    }
    return;
} // destroyPETScVector

inline SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> >
PETScSAMRAIVectorReal::getSAMRAIVector(Vec petsc_vec)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(petsc_vec);
#endif
    PETScSAMRAIVectorReal* psv = static_cast<PETScSAMRAIVectorReal*>(petsc_vec->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(psv);
#endif
    return psv->d_samrai_vector;
} // getSAMRAIVector

inline void
PETScSAMRAIVectorReal::replaceSAMRAIVector(
    Vec petsc_vec,
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> > samrai_vec)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(petsc_vec);
    TBOX_ASSERT(samrai_vec);
#endif
    PETScSAMRAIVectorReal* psv = static_cast<PETScSAMRAIVectorReal*>(petsc_vec->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(psv);
#endif
    psv->d_samrai_vector = samrai_vec;
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_vec));
    IBTK_CHKERRQ(ierr);
    return;
} // getSAMRAIVector

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PETScSAMRAIVectorReal_inl_h
