// Filename: PatchVecCellDataBasicOps.cpp
// Created on 09 Apr 2010 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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

#include <ostream>

#include "Index.h"
#include "IntVector.h"
#include "PatchVecCellDataBasicOps.h"
#include "SAMRAI_config.h"
#include "blitz/array.h"
#include "blitz/array/expr.h"
#include "blitz/array/ops.h"
#include "blitz/array/reduce.h"
#include "blitz/etbase.h"
#include "blitz/range.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {
template <class TYPE> class Pointer;
}  // namespace tbox
}  // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

#if   (NDIM == 2)
#define GET_ARRAY_BOX(arr,box)                          \
    arr(blitz::Range::all(),                            \
        blitz::Range(box.lower()(0),box.upper()(0)),    \
        blitz::Range(box.lower()(1),box.upper()(1)))
#elif (NDIM == 3)
#define GET_ARRAY_BOX(arr,box)                          \
    arr(blitz::Range::all(),                            \
        blitz::Range(box.lower()(0),box.upper()(0)),    \
        blitz::Range(box.lower()(1),box.upper()(1)),    \
        blitz::Range(box.lower()(2),box.upper()(2)))
#endif

/////////////////////////////// PUBLIC ///////////////////////////////////////

template<class TYPE>
PatchVecCellDataBasicOps<TYPE>::PatchVecCellDataBasicOps()
{
    // intentionally blank
    return;
}// PatchVecCellDataBasicOps

template<class TYPE>
PatchVecCellDataBasicOps<TYPE>::~PatchVecCellDataBasicOps()
{
    // intentionally blank
    return;
}// ~PatchVecCellDataBasicOps

template<class TYPE>
void
PatchVecCellDataBasicOps<TYPE>::scale(
    Pointer<VecCellData<TYPE> >& dst,
    const TYPE& alpha,
    const Pointer<VecCellData<TYPE> >& src,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst && src);
#endif
    const Box<NDIM> data_box = box * dst->getGhostBox() * src->getGhostBox();
    blitz::Array<double,NDIM+1> dst_arr(GET_ARRAY_BOX(dst->getArray(),data_box));
    blitz::Array<double,NDIM+1> src_arr(GET_ARRAY_BOX(src->getArray(),data_box));
    dst_arr = alpha*src_arr;
    return;
}// scale

template<class TYPE>
void
PatchVecCellDataBasicOps<TYPE>::addScalar(
    Pointer<VecCellData<TYPE> >& dst,
    const Pointer<VecCellData<TYPE> >& src,
    const TYPE& alpha,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst && src);
#endif
    const Box<NDIM> data_box = box * dst->getGhostBox() * src->getGhostBox();
    blitz::Array<double,NDIM+1> dst_arr(GET_ARRAY_BOX(dst->getArray(),data_box));
    blitz::Array<double,NDIM+1> src_arr(GET_ARRAY_BOX(src->getArray(),data_box));
    dst_arr = src_arr+alpha;
    return;
}// addScalar

template<class TYPE>
void
PatchVecCellDataBasicOps<TYPE>::add(
    Pointer<VecCellData<TYPE> >& dst,
    const Pointer<VecCellData<TYPE> >& src1,
    const Pointer<VecCellData<TYPE> >& src2,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst && src1 && src2);
#endif
    const Box<NDIM> data_box = box * dst->getGhostBox() * src1->getGhostBox() * src2->getGhostBox();
    blitz::Array<double,NDIM+1>  dst_arr(GET_ARRAY_BOX( dst->getArray(),data_box));
    blitz::Array<double,NDIM+1> src1_arr(GET_ARRAY_BOX(src1->getArray(),data_box));
    blitz::Array<double,NDIM+1> src2_arr(GET_ARRAY_BOX(src2->getArray(),data_box));
    dst_arr = src1_arr + src2_arr;
    return;
}// add

template<class TYPE>
void
PatchVecCellDataBasicOps<TYPE>::subtract(
    Pointer<VecCellData<TYPE> >& dst,
    const Pointer<VecCellData<TYPE> >& src1,
    const Pointer<VecCellData<TYPE> >& src2,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst && src1 && src2);
#endif
    const Box<NDIM> data_box = box * dst->getGhostBox() * src1->getGhostBox() * src2->getGhostBox();
    blitz::Array<double,NDIM+1>  dst_arr(GET_ARRAY_BOX( dst->getArray(),data_box));
    blitz::Array<double,NDIM+1> src1_arr(GET_ARRAY_BOX(src1->getArray(),data_box));
    blitz::Array<double,NDIM+1> src2_arr(GET_ARRAY_BOX(src2->getArray(),data_box));
    dst_arr = src1_arr - src2_arr;
    return;
}// subtract

template<class TYPE>
void
PatchVecCellDataBasicOps<TYPE>::multiply(
    Pointer<VecCellData<TYPE> >& dst,
    const Pointer<VecCellData<TYPE> >& src1,
    const Pointer<VecCellData<TYPE> >& src2,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst && src1 && src2);
#endif
    const Box<NDIM> data_box = box * dst->getGhostBox() * src1->getGhostBox() * src2->getGhostBox();
    blitz::Array<double,NDIM+1>  dst_arr(GET_ARRAY_BOX( dst->getArray(),data_box));
    blitz::Array<double,NDIM+1> src1_arr(GET_ARRAY_BOX(src1->getArray(),data_box));
    blitz::Array<double,NDIM+1> src2_arr(GET_ARRAY_BOX(src2->getArray(),data_box));
    dst_arr = src1_arr * src2_arr;
    return;
}// multiply

template<class TYPE>
void
PatchVecCellDataBasicOps<TYPE>::divide(
    Pointer<VecCellData<TYPE> >& dst,
    const Pointer<VecCellData<TYPE> >& src1,
    const Pointer<VecCellData<TYPE> >& src2,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst && src1 && src2);
#endif
    const Box<NDIM> data_box = box * dst->getGhostBox() * src1->getGhostBox() * src2->getGhostBox();
    blitz::Array<double,NDIM+1>  dst_arr(GET_ARRAY_BOX( dst->getArray(),data_box));
    blitz::Array<double,NDIM+1> src1_arr(GET_ARRAY_BOX(src1->getArray(),data_box));
    blitz::Array<double,NDIM+1> src2_arr(GET_ARRAY_BOX(src2->getArray(),data_box));
    dst_arr = src1_arr / src2_arr;
    return;
}// divide

template<class TYPE>
void
PatchVecCellDataBasicOps<TYPE>::reciprocal(
    Pointer<VecCellData<TYPE> >& dst,
    const Pointer<VecCellData<TYPE> >& src,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst && src);
#endif
    const Box<NDIM> data_box = box * dst->getGhostBox() * src->getGhostBox();
    blitz::Array<double,NDIM+1> dst_arr(GET_ARRAY_BOX(dst->getArray(),data_box));
    blitz::Array<double,NDIM+1> src_arr(GET_ARRAY_BOX(src->getArray(),data_box));
    dst_arr = 1.0 / src_arr;
    return;
}// reciprocal

template<class TYPE>
void
PatchVecCellDataBasicOps<TYPE>::linearSum(
    Pointer<VecCellData<TYPE> >& dst,
    const TYPE& alpha,
    const Pointer<VecCellData<TYPE> >& src1,
    const TYPE& beta,
    const Pointer<VecCellData<TYPE> >& src2,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst && src1 && src2);
#endif
    const Box<NDIM> data_box = box * dst->getGhostBox() * src1->getGhostBox() * src2->getGhostBox();
    blitz::Array<double,NDIM+1>  dst_arr(GET_ARRAY_BOX( dst->getArray(),data_box));
    blitz::Array<double,NDIM+1> src1_arr(GET_ARRAY_BOX(src1->getArray(),data_box));
    blitz::Array<double,NDIM+1> src2_arr(GET_ARRAY_BOX(src2->getArray(),data_box));
    dst_arr = (alpha*src1_arr) + (beta*src2_arr);
    return;
}// linearSum

template<class TYPE>
void
PatchVecCellDataBasicOps<TYPE>::axpy(
    Pointer<VecCellData<TYPE> >& dst,
    const TYPE& alpha,
    const Pointer<VecCellData<TYPE> >& src1,
    const Pointer<VecCellData<TYPE> >& src2,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst && src1 && src2);
#endif
    const Box<NDIM> data_box = box * dst->getGhostBox() * src1->getGhostBox() * src2->getGhostBox();
    blitz::Array<double,NDIM+1>  dst_arr(GET_ARRAY_BOX( dst->getArray(),data_box));
    blitz::Array<double,NDIM+1> src1_arr(GET_ARRAY_BOX(src1->getArray(),data_box));
    blitz::Array<double,NDIM+1> src2_arr(GET_ARRAY_BOX(src2->getArray(),data_box));
    dst_arr = (alpha*src1_arr) + src2_arr;
    return;
}// axpy

template<class TYPE>
void
PatchVecCellDataBasicOps<TYPE>::axmy(
    Pointer<VecCellData<TYPE> >& dst,
    const TYPE& alpha,
    const Pointer<VecCellData<TYPE> >& src1,
    const Pointer<VecCellData<TYPE> >& src2,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst && src1 && src2);
#endif
    const Box<NDIM> data_box = box * dst->getGhostBox() * src1->getGhostBox() * src2->getGhostBox();
    blitz::Array<double,NDIM+1>  dst_arr(GET_ARRAY_BOX( dst->getArray(),data_box));
    blitz::Array<double,NDIM+1> src1_arr(GET_ARRAY_BOX(src1->getArray(),data_box));
    blitz::Array<double,NDIM+1> src2_arr(GET_ARRAY_BOX(src2->getArray(),data_box));
    dst_arr = (alpha*src1_arr) - src2_arr;
    return;
}// axmy

template<class TYPE>
TYPE
PatchVecCellDataBasicOps<TYPE>::min(
    const Pointer<VecCellData<TYPE> >& data,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data);
#endif
    const Box<NDIM> data_box = box * data->getGhostBox();
    blitz::Array<double,NDIM+1> data_arr(GET_ARRAY_BOX(data->getArray(),data_box));
    return blitz::min(data_arr);
}// min

template<class TYPE>
TYPE
PatchVecCellDataBasicOps<TYPE>::max(
    const Pointer<VecCellData<TYPE> >& data,
    const Box<NDIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data);
#endif
    const Box<NDIM> data_box = box * data->getGhostBox();
    blitz::Array<double,NDIM+1> data_arr(GET_ARRAY_BOX(data->getArray(),data_box));
    return blitz::max(data_arr);
}// max

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include "ibtk/VecCellData.h"  // IWYU pragma: keep

template class IBTK::PatchVecCellDataBasicOps<double>;

//////////////////////////////////////////////////////////////////////////////
