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
      d_node_data(node_data)
{
    // intentionally blank
    return;
}// LNode

inline
LNode::LNode(
    const LNode& from)
    : LNodeIndex(from),
      d_node_data(from.d_node_data)
{
    // intentionally blank
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
    for (std::vector<SAMRAI::tbox::Pointer<Streamable> >::iterator it = d_node_data.begin();
         it != d_node_data.end(); ++it)
    {
        (*it)->registerPeriodicShift(offset, displacement);
    }
    return;
}// registerPeriodicShift

inline const std::vector<SAMRAI::tbox::Pointer<Streamable> >&
LNode::getNodeData() const
{
    return d_node_data;
}// getNodeData

inline std::vector<SAMRAI::tbox::Pointer<Streamable> >&
LNode::getNodeData()
{
    return d_node_data;
}// getNodeData

inline void
LNode::setNodeData(
    const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data)
{
    d_node_data = node_data;
    return;
}// setNodeData

template<typename T>
inline T*
LNode::getNodeDataItem() const
{
    T* ret_val = NULL;
    for (std::vector<SAMRAI::tbox::Pointer<Streamable> >::const_iterator cit = d_node_data.begin();
         cit != d_node_data.end(); ++cit)
    {
        ret_val = dynamic_cast<T*>(cit->getPointer());
        if (ret_val != NULL) break;
    }
    return ret_val;
}// getNodeDataItem

template<typename T>
inline std::vector<T*>
LNode::getNodeDataVector() const
{
    std::vector<T*> ret_val;
    for (std::vector<SAMRAI::tbox::Pointer<Streamable> >::const_iterator cit = d_node_data.begin();
         cit != d_node_data.end(); ++cit)
    {
        T* it_val = dynamic_cast<T*>(cit->getPointer());
        if (it_val != NULL) ret_val.push_back(it_val);
    }
    return ret_val;
}// getNodeDataVector

inline void
LNode::copySourceItem(
    const SAMRAI::hier::Index<NDIM>& src_index,
    const SAMRAI::hier::IntVector<NDIM>& src_offset,
    const LNode& src_item)
{
    LNodeIndex::copySourceItem(src_index, src_offset, src_item);
    assignThatToThis(src_item);
    return;
}// copySourceItem

inline size_t
LNode::getDataStreamSize() const
{
    return LNodeIndex::getDataStreamSize() + StreamableManager::getManager()->getDataStreamSize(d_node_data);
}// getDataStreamSize

inline void
LNode::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    LNodeIndex::packStream(stream);
    StreamableManager::getManager()->packStream(stream, d_node_data);
    return;
}// packStream

inline void
LNode::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    LNodeIndex::unpackStream(stream, offset);
    d_node_data.clear();
    StreamableManager::getManager()->unpackStream(stream, offset, d_node_data);
    return;
}// unpackStream

/////////////////////////////// PRIVATE //////////////////////////////////////

inline void
LNode::assignThatToThis(
    const LNode& that)
{
    d_node_data = that.d_node_data;
    return;
}// assignThatToThis

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
