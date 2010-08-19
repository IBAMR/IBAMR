// Filename: IBAnchorPointSpec.C
// Created on 18 Aug 2008 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBAnchorPointSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
IBAnchorPointSpec::IBAnchorPointSpec(
    const int node_idx)
    : d_node_idx(node_idx)
{
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBAnchorPointSpec::IBAnchorPointSpec():\n"
                   << "  must call IBAnchorPointSpec::registerWithStashableManager() before\n"
                   << "  creating any IBAnchorPointSpec objects.\n");
    }
    return;
}// IBAnchorPointSpec

inline
IBAnchorPointSpec::~IBAnchorPointSpec()
{
    // intentionally blank
    return;
}// ~IBAnchorPointSpec

inline const int&
IBAnchorPointSpec::getNodeIndex() const
{
    return d_node_idx;
}// getNodeIndex

inline int&
IBAnchorPointSpec::getNodeIndex()
{
    return d_node_idx;
}// getNodeIndex

inline int
IBAnchorPointSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBAnchorPointSpec::getDataStreamSize() const
{
    return SAMRAI::tbox::AbstractStream::sizeofInt();
}// getDataStreamSize

inline void
IBAnchorPointSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_node_idx,1);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
