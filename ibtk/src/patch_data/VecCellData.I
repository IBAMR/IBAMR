// Filename: VecCellData.I
// Created on 09 Apr 2010 by Boyce Griffith
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

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

template<class TYPE>
inline unsigned int
VecCellData<TYPE>::getDepth() const
{
    return d_depth;
}// getDepth

template<class TYPE>
inline TYPE*
VecCellData<TYPE>::getPointer()
{
    return d_data.data();
}// getPointer

template<class TYPE>
inline const TYPE*
VecCellData<TYPE>::getPointer() const
{
    return d_data.data();
}// getPointer

template<class TYPE>
inline blitz::Array<TYPE,1>
VecCellData<TYPE>::operator()(
    const SAMRAI::pdat::CellIndex<NDIM>& i)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(getGhostBox().contains(i));
#endif
#if   (NDIM == 2)
    return d_data(blitz::Range::all(),i(0),i(1));
#elif (NDIM == 3)
    return d_data(blitz::Range::all(),i(0),i(1),i(2));
#else
    TBOX_ERROR("VecCellData implemented only for NDIM=2,3\n");
    static blitz::Array<TYPE,1> dummy;
    return dummy;
#endif
}// operator()

template<class TYPE>
inline TYPE&
VecCellData<TYPE>::operator()(
    const SAMRAI::pdat::CellIndex<NDIM>& i,
    const unsigned int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(depth < d_depth);
    TBOX_ASSERT(getGhostBox().contains(i));
#endif
#if   (NDIM == 2)
    return d_data(static_cast<int>(depth),i(0),i(1));
#elif (NDIM == 3)
    return d_data(static_cast<int>(depth),i(0),i(1),i(2));
#else
    TBOX_ERROR("VecCellData implemented only for NDIM=2,3\n");
    static TYPE dummy;
    return dummy;
#endif
}// operator()

template<class TYPE>
inline const TYPE&
VecCellData<TYPE>::operator()(
    const SAMRAI::pdat::CellIndex<NDIM>& i,
    const unsigned int depth) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(depth < d_depth);
    TBOX_ASSERT(getGhostBox().contains(i));
#endif
#if   (NDIM == 2)
    return d_data(static_cast<int>(depth),i(0),i(1));
#elif (NDIM == 3)
    return d_data(static_cast<int>(depth),i(0),i(1),i(2));
#else
    TBOX_ERROR("VecCellData implemented only for NDIM=2,3\n");
    static TYPE dummy;
    return dummy;
#endif
}// operator()

template<class TYPE>
inline blitz::Array<TYPE,NDIM+1>&
VecCellData<TYPE>::getArray()
{
    return d_data;
}// getArray

template<class TYPE>
inline const blitz::Array<TYPE,NDIM+1>&
VecCellData<TYPE>::getArray() const
{
    return d_data;
}// getArray

template<class TYPE>
inline void
VecCellData<TYPE>::copyOnBox(
    const VecCellData<TYPE>& src,
    const SAMRAI::hier::Box<NDIM>& box)
{
    const SAMRAI::hier::Box<NDIM> copy_box = box * getGhostBox() * src.getGhostBox();
    const SAMRAI::hier::Index<NDIM>& lower = copy_box.lower();
    const SAMRAI::hier::Index<NDIM>& upper = copy_box.upper();
#if   (NDIM == 2)
    this->  d_data(blitz::Range::all(),blitz::Range(lower(0),upper(0)),blitz::Range(lower(1),upper(1))) =
        src.d_data(blitz::Range::all(),blitz::Range(lower(0),upper(0)),blitz::Range(lower(1),upper(1)));
#elif (NDIM == 3)
    this->  d_data(blitz::Range::all(),blitz::Range(lower(0),upper(0)),blitz::Range(lower(1),upper(1)),blitz::Range(lower(2),upper(2))) =
        src.d_data(blitz::Range::all(),blitz::Range(lower(0),upper(0)),blitz::Range(lower(1),upper(1)),blitz::Range(lower(2),upper(2)));
#else
    TBOX_ERROR("VecCellData implemented only for NDIM=2,3\n");
#endif
    return;
}// copyOnBox

template<class TYPE>
inline void
VecCellData<TYPE>::copyOnBox(
    const VecCellData<TYPE>& src,
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::hier::IntVector<NDIM>& src_offset)
{
    const SAMRAI::hier::Box<NDIM> copy_box = box * getGhostBox() * SAMRAI::hier::Box<NDIM>::shift(src.getGhostBox(),src_offset);
    const SAMRAI::hier::Index<NDIM>& dst_lower = copy_box.lower();
    const SAMRAI::hier::Index<NDIM>& dst_upper = copy_box.upper();
    const SAMRAI::hier::Index<NDIM>  src_lower = dst_lower - src_offset;
    const SAMRAI::hier::Index<NDIM>  src_upper = dst_upper - src_offset;
#if   (NDIM == 2)
    this->  d_data(blitz::Range::all(),blitz::Range(dst_lower(0),dst_upper(0)),blitz::Range(dst_lower(1),dst_upper(1))) =
        src.d_data(blitz::Range::all(),blitz::Range(src_lower(0),src_upper(0)),blitz::Range(src_lower(1),src_upper(1)));
#elif (NDIM == 3)
    this->  d_data(blitz::Range::all(),blitz::Range(dst_lower(0),dst_upper(0)),blitz::Range(dst_lower(1),dst_upper(1)),blitz::Range(dst_lower(2),dst_upper(2))) =
        src.d_data(blitz::Range::all(),blitz::Range(src_lower(0),src_upper(0)),blitz::Range(src_lower(1),src_upper(1)),blitz::Range(src_lower(2),src_upper(2)));
#else
    TBOX_ERROR("VecCellData implemented only for NDIM=2,3\n");
#endif
    return;
}// copyOnBox

template<class TYPE>
inline void
VecCellData<TYPE>::copyOnBox(
    const SAMRAI::pdat::CellData<NDIM,TYPE>& src,
    const SAMRAI::hier::Box<NDIM>& box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(src.getDepth() >= 0);
    TBOX_ASSERT(d_depth == static_cast<unsigned int>(src.getDepth()));
#endif
    const SAMRAI::hier::Box<NDIM> copy_box = box * getGhostBox() * src.getGhostBox();
    for (Iterator i(copy_box); i; i++)
    {
        const SAMRAI::hier::Index<NDIM>& idx = i();
        for (unsigned int d = 0; d < d_depth; ++d)
        {
            (*this)(idx,d) = src(idx,d);
        }
    }
    return;
}// copyOnBox

template<class TYPE>
inline void
VecCellData<TYPE>::copy2OnBox(
    SAMRAI::pdat::CellData<NDIM,TYPE>& dst,
    const SAMRAI::hier::Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst.getDepth() >= 0);
    TBOX_ASSERT(d_depth == static_cast<unsigned int>(dst.getDepth()));
#endif
    const SAMRAI::hier::Box<NDIM> copy_box = box * getGhostBox() * dst.getGhostBox();
    for (Iterator i(copy_box); i; i++)
    {
        const SAMRAI::hier::Index<NDIM>& idx = i();
        for (unsigned int d = 0; d < d_depth; ++d)
        {
            dst(idx,d) = (*this)(idx,d);
        }
    }
    return;
}// copy2OnBox

template<class TYPE>
inline void
VecCellData<TYPE>::fill(
    const TYPE& t,
    const unsigned int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d < d_depth);
#endif
#if   (NDIM == 2)
    d_data(d,blitz::Range::all(),blitz::Range::all()) = t;
#elif (NDIM == 3)
    d_data(d,blitz::Range::all(),blitz::Range::all(),blitz::Range::all()) = t;
#else
    TBOX_ERROR("VecCellData implemented only for NDIM=2,3\n");
#endif
    return;
}// fill

template<class TYPE>
inline void
VecCellData<TYPE>::fill(
    const TYPE& t,
    const SAMRAI::hier::Box<NDIM>& box,
    const unsigned int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d < d_depth);
#endif
    const SAMRAI::hier::Index<NDIM>& lower = box.lower();
    const SAMRAI::hier::Index<NDIM>& upper = box.upper();
#if   (NDIM == 2)
    d_data(d,blitz::Range(lower(0),upper(0)),blitz::Range(lower(1),upper(1))) = t;
#elif (NDIM == 3)
    d_data(d,blitz::Range(lower(0),upper(0)),blitz::Range(lower(1),upper(1)),blitz::Range(lower(2),upper(2))) = t;
#else
    TBOX_ERROR("VecCellData implemented only for NDIM=2,3\n");
#endif
    return;
}// fill

template<class TYPE>
inline void VecCellData<TYPE>::fillAll(
    const TYPE& t)
{
    d_data = t;
    return;
}// fillAll

template<class TYPE>
inline void
VecCellData<TYPE>::fillAll(
    const TYPE& t,
    const SAMRAI::hier::Box<NDIM>& box)
{
    const SAMRAI::hier::Index<NDIM>& lower = box.lower();
    const SAMRAI::hier::Index<NDIM>& upper = box.upper();
#if   (NDIM == 2)
    d_data(blitz::Range::all(),blitz::Range(lower(0),upper(0)),blitz::Range(lower(1),upper(1))) = t;
#elif (NDIM == 3)
    d_data(blitz::Range::all(),blitz::Range(lower(0),upper(0)),blitz::Range(lower(1),upper(1)),blitz::Range(lower(2),upper(2))) = t;
#else
    TBOX_ERROR("VecCellData implemented only for NDIM=2,3\n");
#endif
    return;
}// fillAll

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
