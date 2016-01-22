// Filename: LMarker-inl.h
// Created on 12 Sep 2007 by Boyce Griffith
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

#ifndef included_LMarker_inl_h
#define included_LMarker_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LMarker.h"
#include "tbox/AbstractStream.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline LMarker::LMarker(const int idx,
                        const Point& X,
                        const Vector& U,
                        const SAMRAI::hier::IntVector<NDIM>& periodic_offset)
    : d_idx(idx), d_X(X), d_U(U), d_offset(periodic_offset)
{
    // intentionally blank
    return;
} // LMarker

inline LMarker::LMarker(const LMarker& from) : d_idx(from.d_idx), d_X(from.d_X), d_U(from.d_U), d_offset(from.d_offset)
{
    // intentionally blank
    return;
} // LMarker

inline LMarker::LMarker(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& offset)
    : d_idx(-1), d_X(Point::Zero()), d_U(Vector::Zero()), d_offset(offset)
{
    unpackStream(stream, offset);
    return;
} // LMarker

inline LMarker::~LMarker()
{
    // intentionally blank
    return;
} // ~LMarker

inline LMarker& LMarker::operator=(const LMarker& that)
{
    if (this == &that) return *this; // check for self-assignment
    d_idx = that.d_idx;
    d_X = that.d_X;
    d_U = that.d_U;
    d_offset = that.d_offset;
    return *this;
} // operator=

inline const int&
LMarker::getIndex() const
{
    return d_idx;
} // getIndex

inline int&
LMarker::getIndex()
{
    return d_idx;
} // getIndex

inline void
LMarker::setIndex(int idx)
{
    d_idx = idx;
    return;
} // setIndex

inline const Point&
LMarker::getPosition() const
{
    return d_X;
} // getPosition

inline Point&
LMarker::getPosition()
{
    return d_X;
} // getPosition

inline void
LMarker::setPosition(const Point& X)
{
    d_X = X;
    return;
} // setPosition

inline const Vector&
LMarker::getVelocity() const
{
    return d_U;
} // getVelocity

inline Vector&
LMarker::getVelocity()
{
    return d_U;
} // getVelocity

inline void
LMarker::setVelocity(const Vector& U)
{
    d_U = U;
    return;
} // setVelocity

inline const SAMRAI::hier::IntVector<NDIM>&
LMarker::getPeriodicOffset() const
{
    return d_offset;
} // getPeriodicOffset

inline void
LMarker::setPeriodicOffset(const SAMRAI::hier::IntVector<NDIM>& offset)
{
    d_offset = offset;
    return;
} // setPeriodicOffset

inline void
LMarker::copySourceItem(const SAMRAI::hier::Index<NDIM>& /*src_index*/,
                        const SAMRAI::hier::IntVector<NDIM>& src_offset,
                        const LMarker& src_item)
{
    d_idx = src_item.d_idx;
    d_X = src_item.d_X;
    d_U = src_item.d_U;
    d_offset = src_offset;
    return;
} // copySourceItem

inline size_t
LMarker::getDataStreamSize() const
{
    return (1 * SAMRAI::tbox::AbstractStream::sizeofInt() + 2 * NDIM * SAMRAI::tbox::AbstractStream::sizeofDouble());
} // getDataStreamSize

inline void
LMarker::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_idx, 1);
    stream.pack(d_X.data(), NDIM);
    stream.pack(d_U.data(), NDIM);
    return;
} // packStream

inline void
LMarker::unpackStream(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& /*offset*/)
{
    stream.unpack(&d_idx, 1);
    stream.unpack(d_X.data(), NDIM);
    stream.unpack(d_U.data(), NDIM);
    return;
} // unpackStream

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LMarker_inl_h
