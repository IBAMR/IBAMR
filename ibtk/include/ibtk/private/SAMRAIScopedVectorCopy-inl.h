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
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIPointer.h>
#include <SAMRAISAMRAIVectorReal.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
template <typename TYPE>
SAMRAIScopedVectorCopy<TYPE>::SAMRAIScopedVectorCopy(const SAMRAIPointer<SAMRAISAMRAIVectorReal<TYPE>>& vector,
                                                     const std::string& name)
    : SAMRAIScopedVectorCopy(checked_dereference(vector), name)
{
}

template <typename TYPE>
SAMRAIScopedVectorCopy<TYPE>::SAMRAIScopedVectorCopy(const SAMRAISAMRAIVectorReal<TYPE>& vector,
                                                     const std::string& name)
    : SAMRAIScopedVectorDuplicate<TYPE>(vector, name)
{
    this->d_vector->copyVector(
        SAMRAIPointer<SAMRAISAMRAIVectorReal<TYPE>>(const_cast<SAMRAISAMRAIVectorReal<TYPE>*>(&vector), false), false);
}
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_SAMRAIScopedVectorCopy_inl_h
