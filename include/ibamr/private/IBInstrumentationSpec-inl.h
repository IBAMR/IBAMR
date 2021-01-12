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

#ifndef included_IBAMR_IBInstrumentationSpec_inl
#define included_IBAMR_IBInstrumentationSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/StreamableManager.h"

#include "SAMRAI_config.h"
#include "tbox/AbstractStream.h"
#include "tbox/Utilities.h"

#include <string>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBInstrumentationSpec::getIsRegisteredWithStreamableManager()
{
    return (STREAMABLE_CLASS_ID != IBTK::StreamableManager::getUnregisteredID());
} // getIsRegisteredWithStreamableManager

inline const std::vector<std::string>&
IBInstrumentationSpec::getInstrumentNames()
{
    return s_instrument_names;
} // getInstrumentNames

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline IBInstrumentationSpec::IBInstrumentationSpec(const int master_idx, const int meter_idx, const int node_idx)
    : d_master_idx(master_idx), d_meter_idx(meter_idx), d_node_idx(node_idx)
{
#if !defined(NDEBUG)
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBInstrumentationSpec::IBInstrumentationSpec():\n"
                   << "  must call IBInstrumentationSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBInstrumentationSpec objects.\n");
    }
#endif
    return;
} // IBInstrumentationSpec

inline IBInstrumentationSpec::~IBInstrumentationSpec()
{
    // intentionally blank
    return;
} // ~IBInstrumentationSpec

inline const int&
IBInstrumentationSpec::getMasterNodeIndex() const
{
    return d_master_idx;
} // getMasterNodeIndex

inline int&
IBInstrumentationSpec::getMasterNodeIndex()
{
    return d_master_idx;
} // getMasterNodeIndex

inline const int&
IBInstrumentationSpec::getMeterIndex() const
{
    return d_meter_idx;
} // getMeterIndex

inline int&
IBInstrumentationSpec::getMeterIndex()
{
    return d_meter_idx;
} // getMeterIndex

inline const int&
IBInstrumentationSpec::getNodeIndex() const
{
    return d_node_idx;
} // getNodeIndex

inline int&
IBInstrumentationSpec::getNodeIndex()
{
    return d_node_idx;
} // getNodeIndex

inline int
IBInstrumentationSpec::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
} // getStreamableClassID

inline size_t
IBInstrumentationSpec::getDataStreamSize() const
{
    return 3 * SAMRAI::tbox::AbstractStream::sizeofInt();
} // getDataStreamSize

inline void
IBInstrumentationSpec::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_master_idx, 1);
    stream.pack(&d_meter_idx, 1);
    stream.pack(&d_node_idx, 1);
    return;
} // packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// NAMESPACE ////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBInstrumentationSpec_inl
