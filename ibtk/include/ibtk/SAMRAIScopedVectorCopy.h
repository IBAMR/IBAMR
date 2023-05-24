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

#ifndef included_IBTK_SAMRAIScopedVectorCopy
#define included_IBTK_SAMRAIScopedVectorCopy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/SAMRAIScopedVectorDuplicate.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <SAMRAIVectorReal.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/**
 * Wrapper class around a SAMRAIVectorReal with RAII semantics. Creates a new
 * copy of an input vector.
 *
 * @note The name of this class is analogous to the meaning of copy in PETSc's
 * VecCopy() function.
 */
template <typename TYPE>
class SAMRAIScopedVectorCopy : public SAMRAIScopedVectorDuplicate<TYPE>
{
public:
    /*!
     * Constructor. Sets up a vector equivalent to @p vector.
     */
    SAMRAIScopedVectorCopy(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> >& vector,
                           const std::string& name = "");

    /*!
     * Constructor. Sets up a vector equivalent to @p vector.
     */
    SAMRAIScopedVectorCopy(const SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE>& vector, const std::string& name = "");
};
} // namespace IBTK

#include "ibtk/private/SAMRAIScopedVectorCopy-inl.h"

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_SAMRAIScopedVectorCopy
