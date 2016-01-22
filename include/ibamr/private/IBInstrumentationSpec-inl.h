// Filename: IBInstrumentationSpec-inl.h
// Created on 11 Jun 2007 by Boyce Griffith
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

#ifndef included_IBInstrumentationSpec_inl
#define included_IBInstrumentationSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <unistd.h>
#include <ostream>
#include <string>
#include <vector>

#include "SAMRAI_config.h"
#include "ibtk/StreamableManager.h"
#include "tbox/AbstractStream.h"
#include "tbox/Utilities.h"

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

#endif //#ifndef included_IBInstrumentationSpec_inl
