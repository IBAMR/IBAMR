// Filename: LNodeIndex.I
// Created on 28 Feb 2004 by Boyce Griffith
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
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
LNodeIndex::LNodeIndex(
    const int lagrangian_nidx,
    const int local_petsc_nidx,
    double* const X_ptr,
    const SAMRAI::hier::IntVector<NDIM>& periodic_offset,
    const std::vector<double>& periodic_displacement,
    const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data)
    : d_lagrangian_nidx(lagrangian_nidx),
      d_local_petsc_nidx(local_petsc_nidx),
      d_X_ptr(X_ptr),
      d_offset(periodic_offset),
      d_displacement(periodic_displacement),
      d_node_data(node_data)
{
    // intentionally blank
    return;
}// LNodeIndex

inline
LNodeIndex::LNodeIndex(
    const LNodeIndex& from)
    : d_lagrangian_nidx(from.d_lagrangian_nidx),
      d_local_petsc_nidx(from.d_local_petsc_nidx),
      d_X_ptr(from.d_X_ptr),
      d_offset(from.d_offset),
      d_displacement(from.d_displacement),
      d_node_data(from.d_node_data)
{
    // intentionally blank
    return;
}// LNodeIndex

inline
LNodeIndex::~LNodeIndex()
{
    // intentionally blank
    return;
}// ~LNodeIndex

inline LNodeIndex&
LNodeIndex::operator=(
    const LNodeIndex& that)
{
    if (this != &that)
    {
        assignThatToThis(that);
    }
    return *this;
}// operator=

inline int
LNodeIndex::getLagrangianIndex() const
{
    return d_lagrangian_nidx;
}// getLagrangianIndex

inline void
LNodeIndex::setLagrangianIndex(
    const int lagrangian_nidx)
{
    d_lagrangian_nidx = lagrangian_nidx;
    return;
}// setLagrangianIndex

inline int
LNodeIndex::getLocalPETScIndex() const
{
    return d_local_petsc_nidx;
}// getLocalPETScIndex

inline void
LNodeIndex::setLocalPETScIndex(
    const int local_petsc_nidx)
{
    d_local_petsc_nidx = local_petsc_nidx;
    return;
}// setLocalPETScIndex

inline double*
LNodeIndex::getNodeLocation() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_X_ptr != NULL);
#endif
    return d_X_ptr;
}// getNodeLocation

inline void
LNodeIndex::setNodeLocation(
    double* const X_ptr)
{
    d_X_ptr = X_ptr;
    return;
}// setNodeLocation

inline void
LNodeIndex::registerPeriodicShift(
    const SAMRAI::hier::IntVector<NDIM>& offset,
    const std::vector<double>& displacement)
{
    d_offset += offset;
    for (int d = 0; d < NDIM; ++d)
    {
        d_displacement[d] += displacement[d];
    }
    for (std::vector<SAMRAI::tbox::Pointer<Streamable> >::iterator it = d_node_data.begin();
         it != d_node_data.end(); ++it)
    {
        (*it)->registerPeriodicShift(offset, displacement);
    }
    return;
}// registerPeriodicShift

inline const SAMRAI::hier::IntVector<NDIM>&
LNodeIndex::getPeriodicOffset() const
{
    return d_offset;
}// getPeriodicOffset

inline const std::vector<double>&
LNodeIndex::getPeriodicDisplacement() const
{
    return d_displacement;
}// getPeriodicDisplacement

inline const std::vector<SAMRAI::tbox::Pointer<Streamable> >&
LNodeIndex::getNodeData() const
{
    return d_node_data;
}// getNodeData

inline std::vector<SAMRAI::tbox::Pointer<Streamable> >&
LNodeIndex::getNodeData()
{
    return d_node_data;
}// getNodeData

inline void
LNodeIndex::setNodeData(
    const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data)
{
    d_node_data = node_data;
    return;
}// setNodeData

template<class T>
inline SAMRAI::tbox::Pointer<T>
LNodeIndex::getNodeData() const
{
    SAMRAI::tbox::Pointer<T> ret_val = NULL;
    for (std::vector<SAMRAI::tbox::Pointer<Streamable> >::const_iterator cit = d_node_data.begin();
         cit != d_node_data.end(); ++cit)
    {
        ret_val = (*cit);
        if (!ret_val.isNull()) break;
    }
    return ret_val;
}// getNodeData

template<class T>
inline std::vector<SAMRAI::tbox::Pointer<T> >
LNodeIndex::getNodeDataVector() const
{
    std::vector<SAMRAI::tbox::Pointer<T> > ret_val;
    for (std::vector<SAMRAI::tbox::Pointer<Streamable> >::const_iterator cit = d_node_data.begin();
         cit != d_node_data.end(); ++cit)
    {
        SAMRAI::tbox::Pointer<T> it_val = (*cit);
        if (!it_val.isNull()) ret_val.push_back(it_val);
    }
    return ret_val;
}// getNodeDataVector

inline const std::vector<SAMRAI::tbox::Pointer<Streamable> >&
LNodeIndex::getStashData() const
{
    IBTK_DEPRECATED_MEMBER_FUNCTION2("IBTK::LNodeIndex","getStashData()","getNodeData()");
    return getNodeData();
}// getStashData

inline std::vector<SAMRAI::tbox::Pointer<Streamable> >&
LNodeIndex::getStashData()
{
    IBTK_DEPRECATED_MEMBER_FUNCTION2("IBTK::LNodeIndex","getStashData()","getNodeData()");
    return getNodeData();
}// getStashData

inline void
LNodeIndex::setStashData(
    const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data)
{
    IBTK_DEPRECATED_MEMBER_FUNCTION2("IBTK::LNodeIndex","setStashData()","setNodeData()");
    setNodeData(node_data);
    return;
}// setStashData

template<class T>
inline SAMRAI::tbox::Pointer<T>
LNodeIndex::getStashData() const
{
    IBTK_DEPRECATED_MEMBER_FUNCTION2("IBTK::LNodeIndex","getStashData()","getNodeData()");
    return getNodeData<T>();
}// getStashData

template<class T>
inline std::vector<SAMRAI::tbox::Pointer<T> >
LNodeIndex::getStashDataVector() const
{
    IBTK_DEPRECATED_MEMBER_FUNCTION2("IBTK::LNodeIndex","getStashDataVector()","getNodeDataVector()");
    return getNodeDataVector<T>();
}// getStashDataVector

inline void
LNodeIndex::copySourceItem(
    const SAMRAI::hier::Index<NDIM>& src_index,
    const SAMRAI::hier::IntVector<NDIM>& src_offset,
    const LNodeIndex& src_item)
{
    assignThatToThis(src_item);
    return;
}// copySourceItem

inline size_t
LNodeIndex::getDataStreamSize() const
{
    return (2+NDIM)*SAMRAI::tbox::AbstractStream::sizeofInt() + NDIM*SAMRAI::tbox::AbstractStream::sizeofDouble()
        + StreamableManager::getManager()->getDataStreamSize(d_node_data);
}// getDataStreamSize

inline void
LNodeIndex::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_lagrangian_nidx, 1);
    stream.pack(&d_local_petsc_nidx, 1);
    stream.pack(d_offset, NDIM);
    stream.pack(&d_displacement[0], NDIM);
    StreamableManager::getManager()->packStream(stream, d_node_data);
    return;
}// packStream

inline void
LNodeIndex::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    stream.unpack(&d_lagrangian_nidx, 1);
    stream.unpack(&d_local_petsc_nidx, 1);
    d_X_ptr = NULL;
    stream.unpack(d_offset, NDIM);
    d_displacement.resize(NDIM);
    stream.unpack(&d_displacement[0], NDIM);
    d_node_data.clear();
    StreamableManager::getManager()->unpackStream(stream, offset, d_node_data);
    return;
}// unpackStream

/////////////////////////////// PRIVATE //////////////////////////////////////

inline void
LNodeIndex::assignThatToThis(
    const LNodeIndex& that)
{
    d_lagrangian_nidx = that.d_lagrangian_nidx;
    d_local_petsc_nidx = that.d_local_petsc_nidx;
    d_X_ptr = that.d_X_ptr;
    d_offset = that.d_offset;
    d_displacement = that.d_displacement;
    d_node_data = that.d_node_data;
    return;
}// assignThatToThis

//////////////////////////////////////////////////////////////////////////////

inline bool
operator<(
    const LNodeIndex& lhs,
    const LNodeIndex& rhs)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(lhs.d_X_ptr != NULL);
    TBOX_ASSERT(rhs.d_X_ptr != NULL);
#if ((NDIM>3) || (NDIM<1))
    TBOX_ERROR("operator<(const LNodeIndex&,const LNodeIndex&): not implemented for NDIM=="
               << NDIM << endl);
#endif
#endif
    return (lhs.d_X_ptr[0]<rhs.d_X_ptr[0]) ||
#if (NDIM>1)
        ((lhs.d_X_ptr[0]==rhs.d_X_ptr[0])&&
         (lhs.d_X_ptr[1]< rhs.d_X_ptr[1])) ||
#if (NDIM>2)
        ((lhs.d_X_ptr[0]==rhs.d_X_ptr[0])&&
         (lhs.d_X_ptr[1]==rhs.d_X_ptr[1])&&
         (lhs.d_X_ptr[2]< rhs.d_X_ptr[2])) ||
#endif
#endif
        ((lhs.d_X_ptr[0]==rhs.d_X_ptr[0])&&
#if (NDIM>1)
         (lhs.d_X_ptr[1]==rhs.d_X_ptr[1])&&
#if (NDIM>2)
         (lhs.d_X_ptr[2]==rhs.d_X_ptr[2])&&
#endif
#endif
         (lhs.d_lagrangian_nidx<rhs.d_lagrangian_nidx));
}// operator<

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
