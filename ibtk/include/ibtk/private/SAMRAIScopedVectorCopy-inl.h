// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
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

#ifndef included_IBTK_SAMRAIScopedVectorCopy_inl_h
#define included_IBTK_SAMRAIScopedVectorCopy_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/SAMRAIScopedVectorCopy.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <SAMRAIVectorReal.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
template <typename TYPE>
SAMRAIScopedVectorCopy<TYPE>::SAMRAIScopedVectorCopy(
    const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> >& vector,
    const std::string& name)
    : SAMRAIScopedVectorCopy(checked_dereference(vector), name)
{
}

template <typename TYPE>
SAMRAIScopedVectorCopy<TYPE>::SAMRAIScopedVectorCopy(const SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE>& vector,
                                                     const std::string& name)
    : SAMRAIScopedVectorDuplicate<TYPE>(vector, name)
{
    this->d_vector->copyVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> >(
                                   const_cast<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE>*>(&vector), false),
                               false);
}
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_SAMRAIScopedVectorCopy_inl_h
