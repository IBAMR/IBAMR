// Filename: LNode.I
// Created on 05 May 2011 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/StreamableManager.h>
#include <ibtk/ibtk_utilities.h>

// SAMRAI INCLUDES
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
LNode::LNode(
    const int lagrangian_nidx,
    const int global_petsc_nidx,
    const int local_petsc_nidx,
    const SAMRAI::hier::IntVector<NDIM>& periodic_offset,
    const blitz::TinyVector<double,NDIM>& periodic_displacement,
    const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data)
    : LNodeIndex(lagrangian_nidx, global_petsc_nidx, local_petsc_nidx, periodic_offset, periodic_displacement),
      d_node_data_sz(node_data.size()),
      d_node_data()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_node_data_sz <= LNode::MAX_SIZE);
#endif
    for (int k = 0; k < d_node_data_sz; ++k)
    {
        d_node_data[k] = node_data[k];
    }
    return;
}// LNode

inline
LNode::LNode(
    const LNode& from)
    : LNodeIndex(from),
      d_node_data_sz(from.d_node_data_sz),
      d_node_data()
{
    for (int k = 0; k < d_node_data_sz; ++k)
    {
        d_node_data[k] = from.d_node_data[k];
    }
    return;
}// LNode

inline
LNode::LNode(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
    : LNodeIndex(),
      d_node_data_sz(0),
      d_node_data()
{
    unpackStream(stream,offset);
    return;
}// LNode

inline
LNode::~LNode()
{
    // intentionally blank
    return;
}// ~LNode

inline LNode&
LNode::operator=(
    const LNode& that)
{
    if (this != &that)
    {
        LNodeIndex::operator=(that);
        assignThatToThis(that);
    }
    return *this;
}// operator=

inline void
LNode::registerPeriodicShift(
    const SAMRAI::hier::IntVector<NDIM>& offset,
    const blitz::TinyVector<double,NDIM>& displacement)
{
    LNodeIndex::registerPeriodicShift(offset, displacement);
    for (int k = 0; k < d_node_data_sz; ++k)
    {
        d_node_data[k]->registerPeriodicShift(offset, displacement);
    }
    return;
}// registerPeriodicShift

inline std::vector<SAMRAI::tbox::Pointer<Streamable> >
LNode::getNodeData() const
{
    return std::vector<SAMRAI::tbox::Pointer<Streamable> >(d_node_data,d_node_data+d_node_data_sz);
}// getNodeData

inline void
LNode::setNodeData(
    const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data)
{
    const int original_data_sz = d_node_data_sz;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(node_data.size()) <= LNode::MAX_SIZE);
#endif
    d_node_data_sz = node_data.size();
    for (int k = 0; k < d_node_data_sz; ++k)
    {
        d_node_data[k] = node_data[k];
    }
    for (int k = d_node_data_sz; k < original_data_sz; ++k)
    {
        d_node_data[k].setNull();
    }
    return;
}// setNodeData

template<typename T>
inline T*
LNode::getNodeDataItem() const
{
    T* ret_val = NULL;
    for (int k = 0; k < d_node_data_sz; ++k)
    {
        ret_val = dynamic_cast<T*>(d_node_data[k].getPointer());
        if (ret_val != NULL) break;
    }
    return ret_val;
}// getNodeDataItem

inline void
LNode::copySourceItem(
    const SAMRAI::hier::Index<NDIM>& src_index,
    const SAMRAI::hier::IntVector<NDIM>& src_offset,
    const LNodeIndex& src_item)
{
    LNodeIndex::copySourceItem(src_index, src_offset, src_item);
    const LNode* const p_src_item = dynamic_cast<const LNode*>(&src_item);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(p_src_item != NULL);
#endif
    assignThatToThis(*p_src_item);
    return;
}// copySourceItem

inline size_t
LNode::getDataStreamSize() const
{
    return LNodeIndex::getDataStreamSize() + SAMRAI::tbox::AbstractStream::sizeofInt() + StreamableManager::getManager()->getDataStreamSize(d_node_data, d_node_data_sz);
}// getDataStreamSize

inline void
LNode::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    LNodeIndex::packStream(stream);
    stream << d_node_data_sz;
    StreamableManager::getManager()->packStream(stream, d_node_data, d_node_data_sz);
    return;
}// packStream

inline void
LNode::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    LNodeIndex::unpackStream(stream, offset);
    for (int k = 0; k < d_node_data_sz; ++k) d_node_data[k].setNull();
    stream >> d_node_data_sz;
    StreamableManager::getManager()->unpackStream(stream, offset, d_node_data, d_node_data_sz);
    return;
}// unpackStream

/////////////////////////////// PRIVATE //////////////////////////////////////

inline void
LNode::assignThatToThis(
    const LNode& that)
{
    const int original_data_sz = d_node_data_sz;
    d_node_data_sz = that.d_node_data_sz;
    for (int k = 0; k < d_node_data_sz; ++k)
    {
        d_node_data[k] = that.d_node_data[k];
    }
    for (int k = d_node_data_sz; k < original_data_sz; ++k)
    {
        d_node_data[k].setNull();
    }
    return;
}// assignThatToThis

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
