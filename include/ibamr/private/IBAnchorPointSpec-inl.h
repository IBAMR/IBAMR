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

#ifndef included_IBAMR_IBAnchorPointSpec_inl
#define included_IBAMR_IBAnchorPointSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBAnchorPointSpec.h"

#include "ibtk/StreamableManager.h"

#include "tbox/PIO.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBAnchorPointSpec::getIsRegisteredWithStreamableManager()
{
    return (STREAMABLE_CLASS_ID != IBTK::StreamableManager::getUnregisteredID());
} // getIsRegisteredWithStreamableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline IBAnchorPointSpec::IBAnchorPointSpec(const int node_idx) : d_node_idx(node_idx)
{
#if !defined(NDEBUG)
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBAnchorPointSpec::IBAnchorPointSpec():\n"
                   << "  must call IBAnchorPointSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBAnchorPointSpec objects.\n");
    }
#endif
    return;
} // IBAnchorPointSpec

inline IBAnchorPointSpec::~IBAnchorPointSpec()
{
    // intentionally blank
    return;
} // ~IBAnchorPointSpec

inline const int&
IBAnchorPointSpec::getNodeIndex() const
{
    return d_node_idx;
} // getNodeIndex

inline int&
IBAnchorPointSpec::getNodeIndex()
{
    return d_node_idx;
} // getNodeIndex

inline int
IBAnchorPointSpec::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
} // getStreamableClassID

inline size_t
IBAnchorPointSpec::getDataStreamSize() const
{
    return SAMRAI::tbox::AbstractStream::sizeofInt();
} // getDataStreamSize

inline void
IBAnchorPointSpec::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_node_idx, 1);
    return;
} // packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBAnchorPointSpec_inl
