// Filename: IBSourceSpec-inl.h
// Created on 28 Apr 2011 by Boyce Griffith
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

#ifndef included_IBSourceSpec_inl
#define included_IBSourceSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

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

#endif //#ifndef included_IBSourceSpec_inl
