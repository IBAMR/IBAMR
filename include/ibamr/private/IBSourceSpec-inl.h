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

#ifndef included_IBAMR_IBSourceSpec_inl
#define included_IBAMR_IBSourceSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBSourceSpec.h"

#include "ibtk/StreamableManager.h"

#include "tbox/PIO.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBSourceSpec::getIsRegisteredWithStreamableManager()
{
    return (STREAMABLE_CLASS_ID != IBTK::StreamableManager::getUnregisteredID());
} // getIsRegisteredWithStreamableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline IBSourceSpec::IBSourceSpec(const int master_idx, const int source_idx)
    : d_master_idx(master_idx), d_source_idx(source_idx)
{
#if !defined(NDEBUG)
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBSourceSpec::IBSourceSpec():\n"
                   << "  must call IBSourceSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBSourceSpec objects.\n");
    }
#endif
    return;
} // IBSourceSpec

inline IBSourceSpec::~IBSourceSpec()
{
    // intentionally blank
    return;
} // ~IBSourceSpec

inline const int&
IBSourceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
} // getMasterNodeIndex

inline int&
IBSourceSpec::getMasterNodeIndex()
{
    return d_master_idx;
} // getMasterNodeIndex

inline const int&
IBSourceSpec::getSourceIndex() const
{
    return d_source_idx;
} // getSourceIndex

inline int&
IBSourceSpec::getSourceIndex()
{
    return d_source_idx;
} // getSourceIndex

inline int
IBSourceSpec::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
} // getStreamableClassID

inline size_t
IBSourceSpec::getDataStreamSize() const
{
    return 2 * SAMRAI::tbox::AbstractStream::sizeofInt();
} // getDataStreamSize

inline void
IBSourceSpec::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_master_idx, 1);
    stream.pack(&d_source_idx, 1);
    return;
} // packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// NAMESPACE ////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBSourceSpec_inl
