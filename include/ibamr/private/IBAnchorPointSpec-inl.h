// Filename: IBAnchorPointSpec-inl.h
// Created on 18 Aug 2008 by Boyce Griffith
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

#ifndef included_IBAnchorPointSpec_inl
#define included_IBAnchorPointSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

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

#endif //#ifndef included_IBAnchorPointSpec_inl
