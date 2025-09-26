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
template <typename INPUT_TYPE>
SAMRAIScopedVectorCopy<TYPE>::SAMRAIScopedVectorCopy(
    const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, INPUT_TYPE> >& vector,
    const std::string& name)
    : SAMRAIScopedVectorCopy(checked_dereference(vector), name)
{
}

template <>
template <>
inline SAMRAIScopedVectorCopy<float>::SAMRAIScopedVectorCopy(const SAMRAI::solv::SAMRAIVectorReal<NDIM, float>& vector,
                                                             const std::string& name)
    : SAMRAIScopedVectorDuplicate<float>(vector, name)
{
    this->d_vector->copyVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, float> >(
                                   const_cast<SAMRAI::solv::SAMRAIVectorReal<NDIM, float>*>(&vector), false),
                               false);
}

template <>
template <>
inline SAMRAIScopedVectorCopy<double>::SAMRAIScopedVectorCopy(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& vector,
    const std::string& name)
    : SAMRAIScopedVectorDuplicate<double>(vector, name)
{
    this->d_vector->copyVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >(
                                   const_cast<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>(&vector), false),
                               false);
}

template <typename TYPE>
template <typename INPUT_TYPE>
SAMRAIScopedVectorCopy<TYPE>::SAMRAIScopedVectorCopy(const SAMRAI::solv::SAMRAIVectorReal<NDIM, INPUT_TYPE>& vector,
                                                     const std::string& name)
    : SAMRAIScopedVectorDuplicate<TYPE>(vector, name)
{
    SAMRAIScopedVectorDuplicate<TYPE>::transformFromVector(vector);
}

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_SAMRAIScopedVectorCopy_inl_h
