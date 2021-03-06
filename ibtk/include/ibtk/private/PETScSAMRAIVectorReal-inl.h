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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_PETScSAMRAIVectorReal_inl_h
#define included_IBTK_PETScSAMRAIVectorReal_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

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
}

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
}

inline void
PETScSAMRAIVectorReal::getSAMRAIVector(
    Vec petsc_vec,
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> >* samrai_vec)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(petsc_vec);
#endif
    PETScSAMRAIVectorReal* psv = static_cast<PETScSAMRAIVectorReal*>(petsc_vec->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(psv);
    TBOX_ASSERT(!psv->d_vector_checked_out_read);
#endif
    psv->d_vector_checked_out_read_write = true;
    *samrai_vec = psv->d_samrai_vector;
}

inline void
PETScSAMRAIVectorReal::restoreSAMRAIVector(
    Vec petsc_vec,
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> >* samrai_vec)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(petsc_vec);
#endif
    PETScSAMRAIVectorReal* psv = static_cast<PETScSAMRAIVectorReal*>(petsc_vec->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(psv);
    TBOX_ASSERT(psv->d_vector_checked_out_read_write);
    TBOX_ASSERT(psv->d_samrai_vector.getPointer() == *samrai_vec);
#endif
    psv->d_vector_checked_out_read_write = false;
    *samrai_vec = NULL;
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_vec));
    IBTK_CHKERRQ(ierr);
}

inline void
PETScSAMRAIVectorReal::getSAMRAIVectorRead(
    Vec petsc_vec,
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> >* samrai_vec)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(petsc_vec);
#endif
    PETScSAMRAIVectorReal* psv = static_cast<PETScSAMRAIVectorReal*>(petsc_vec->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(psv);
    TBOX_ASSERT(!psv->d_vector_checked_out_read_write);
#endif
    psv->d_vector_checked_out_read = true;
    *samrai_vec = psv->d_samrai_vector;
}

inline void
PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(
    Vec petsc_vec,
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> >* samrai_vec)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(petsc_vec);
#endif
    PETScSAMRAIVectorReal* psv = static_cast<PETScSAMRAIVectorReal*>(petsc_vec->data);
#if !defined(NDEBUG)
    TBOX_ASSERT(psv);
    TBOX_ASSERT(psv->d_vector_checked_out_read);
    TBOX_ASSERT(psv->d_samrai_vector.getPointer() == *samrai_vec);
#endif
    psv->d_vector_checked_out_read = false;
    *samrai_vec = NULL;
}

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
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PETScSAMRAIVectorReal_inl_h
