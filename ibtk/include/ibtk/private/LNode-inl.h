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

#ifndef included_IBTK_LNode_inl_h
#define included_IBTK_LNode_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LNode.h"
#include "ibtk/StreamableManager.h"
#include "ibtk/compiler_hints.h"

#include "tbox/AbstractStream.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline LNode::LNode(const int lagrangian_nidx,
                    const int global_petsc_nidx,
                    const int local_petsc_nidx,
                    const SAMRAI::hier::IntVector<NDIM>& initial_periodic_offset,
                    const SAMRAI::hier::IntVector<NDIM>& current_periodic_offset,
                    const Vector& initial_periodic_displacement,
                    const Vector& current_periodic_displacement,
                    const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data)
    : LNodeIndex(lagrangian_nidx,
                 global_petsc_nidx,
                 local_petsc_nidx,
                 initial_periodic_offset,
                 current_periodic_offset,
                 initial_periodic_displacement,
                 current_periodic_displacement),
      d_node_data(node_data)
{
    setupNodeDataTypeArray();
    return;
} // LNode

inline LNode::LNode(const LNode& from) : LNodeIndex(from), d_node_data(from.d_node_data)
{
    setupNodeDataTypeArray();
    return;
} // LNode

inline LNode::LNode(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& offset)
    : LNodeIndex(), d_node_data()
{
    unpackStream(stream, offset);
    return;
} // LNode

inline LNode::~LNode()
{
    // intentionally blank
    return;
} // ~LNode

inline LNode&
LNode::operator=(const LNode& that)
{
    if (this != &that)
    {
        LNodeIndex::operator=(that);
        assignThatToThis(that);
    }
    return *this;
} // operator=

inline void
LNode::registerPeriodicShift(const SAMRAI::hier::IntVector<NDIM>& offset, const Vector& displacement)
{
    LNodeIndex::registerPeriodicShift(offset, displacement);
    for (const auto& node : d_node_data)
    {
        node->registerPeriodicShift(offset, displacement);
    }
    return;
} // registerPeriodicShift

inline const std::vector<SAMRAI::tbox::Pointer<Streamable> >&
LNode::getNodeData() const
{
    return d_node_data;
} // getNodeData

inline void
LNode::setNodeData(const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data)
{
    d_node_data = node_data;
    setupNodeDataTypeArray();
    return;
} // setNodeData

inline void
LNode::appendNodeDataItem(const SAMRAI::tbox::Pointer<Streamable>& node_data_item)
{
    d_node_data.push_back(node_data_item);
    const int class_id = node_data_item->getStreamableClassID();
    if (class_id < MAX_SIZE && !d_node_data_type_arr[class_id])
    {
        d_node_data_type_arr[class_id] = node_data_item.getPointer();
    }
    return;
} // appendNodeDataItem

inline void
LNode::removeNodeDataItem(const SAMRAI::tbox::Pointer<Streamable>& node_data_item)
{
    std::vector<SAMRAI::tbox::Pointer<Streamable> >::iterator it =
        std::find(d_node_data.begin(), d_node_data.end(), node_data_item);
    if (it != d_node_data.end())
    {
        d_node_data.erase(it);
        setupNodeDataTypeArray();
    }
    return;
} // removeNodeDataItem

template <typename T>
inline T*
LNode::getNodeDataItem() const
{
    if (LIKELY(T::STREAMABLE_CLASS_ID < MAX_SIZE))
    {
        return static_cast<T*>(d_node_data_type_arr[T::STREAMABLE_CLASS_ID]);
    }
    else
    {
        const size_t node_data_sz = d_node_data.size();
        T* ret_val = nullptr;
        Streamable* it_val;
        size_t k;
        for (k = 0; k < node_data_sz && !ret_val; ++k)
        {
            it_val = d_node_data[k];
            if (it_val->getStreamableClassID() == T::STREAMABLE_CLASS_ID)
            {
                ret_val = static_cast<T*>(it_val);
            }
        }
        return ret_val;
    }
} // getNodeDataItem

template <typename T>
inline std::vector<T*>
LNode::getNodeDataVector() const
{
    const size_t node_data_sz = d_node_data.size();
    std::vector<T*> ret_val;
    ret_val.reserve(node_data_sz);
    Streamable* it_val;
    for (size_t k = 0; k < node_data_sz; ++k)
    {
        it_val = d_node_data[k];
        if (it_val->getStreamableClassID() == T::STREAMABLE_CLASS_ID)
        {
            ret_val.push_back(static_cast<T*>(it_val));
        }
    }
    return ret_val;
} // getNodeDataVector

inline void
LNode::copySourceItem(const SAMRAI::hier::Index<NDIM>& src_index,
                      const SAMRAI::hier::IntVector<NDIM>& src_offset,
                      const LNodeIndex& src_item)
{
    LNodeIndex::copySourceItem(src_index, src_offset, src_item);
    const LNode* const p_src_item = dynamic_cast<const LNode*>(&src_item);
#if !defined(NDEBUG)
    TBOX_ASSERT(p_src_item);
#endif
    assignThatToThis(*p_src_item);
    return;
} // copySourceItem

inline size_t
LNode::getDataStreamSize() const
{
    return LNodeIndex::getDataStreamSize() + StreamableManager::getManager()->getDataStreamSize(d_node_data);
} // getDataStreamSize

inline void
LNode::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    LNodeIndex::packStream(stream);
    StreamableManager::getManager()->packStream(stream, d_node_data);
    return;
} // packStream

inline void
LNode::unpackStream(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& offset)
{
    LNodeIndex::unpackStream(stream, offset);
    d_node_data.clear();
    StreamableManager::getManager()->unpackStream(stream, offset, d_node_data);
    setupNodeDataTypeArray();
    return;
} // unpackStream

/////////////////////////////// PRIVATE //////////////////////////////////////

inline void
LNode::assignThatToThis(const LNode& that)
{
    d_node_data = that.d_node_data;
    setupNodeDataTypeArray();
    return;
} // assignThatToThis

inline void
LNode::setupNodeDataTypeArray()
{
    std::fill(d_node_data_type_arr, d_node_data_type_arr + MAX_SIZE, static_cast<Streamable*>(nullptr));
    Streamable* it_val;
    int class_id;
    for (std::vector<SAMRAI::tbox::Pointer<Streamable> >::const_iterator cit = d_node_data.begin();
         cit != d_node_data.end();
         ++cit)
    {
        it_val = *cit;
        class_id = it_val->getStreamableClassID();
        if (LIKELY(class_id < MAX_SIZE && !d_node_data_type_arr[class_id]))
        {
            d_node_data_type_arr[class_id] = it_val;
        }
    }
    return;
} // setupNodeDataTypeArray

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LNode_inl_h
