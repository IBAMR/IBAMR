// Filename: IBInstrumentationSpec.I
// Created on 11 Jun 2007 by Boyce Griffith
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
IBInstrumentationSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

inline const std::vector<std::string>&
IBInstrumentationSpec::getInstrumentNames()
{
    return s_instrument_names;
}// getInstrumentNames

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
IBInstrumentationSpec::IBInstrumentationSpec(
    const int master_idx,
    const int meter_idx,
    const int node_idx)
    : d_master_idx(master_idx),
      d_meter_idx(meter_idx),
      d_node_idx(node_idx)
{
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBInstrumentationSpec::IBInstrumentationSpec():\n"
                   << "  must call IBInstrumentationSpec::registerWithStashableManager() before\n"
                   << "  creating any IBInstrumentationSpec objects.\n");
    }
    return;
}// IBInstrumentationSpec

inline
IBInstrumentationSpec::~IBInstrumentationSpec()
{
    // intentionally blank
    return;
}// ~IBInstrumentationSpec

inline const int&
IBInstrumentationSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

inline int&
IBInstrumentationSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

inline const int&
IBInstrumentationSpec::getMeterIndex() const
{
    return d_meter_idx;
}// getMeterIndex

inline int&
IBInstrumentationSpec::getMeterIndex()
{
    return d_meter_idx;
}// getMeterIndex

inline const int&
IBInstrumentationSpec::getNodeIndex() const
{
    return d_node_idx;
}// getNodeIndex

inline int&
IBInstrumentationSpec::getNodeIndex()
{
    return d_node_idx;
}// getNodeIndex

inline int
IBInstrumentationSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBInstrumentationSpec::getDataStreamSize() const
{
    return 5*SAMRAI::tbox::AbstractStream::sizeofInt();
}// getDataStreamSize

inline void
IBInstrumentationSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_master_idx,1);
    stream.pack(&d_meter_idx,1);
    stream.pack(&d_node_idx,1);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// NAMESPACE ////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
