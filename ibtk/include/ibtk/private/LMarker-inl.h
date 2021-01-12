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

#ifndef included_IBTK_LMarker_inl_h
#define included_IBTK_LMarker_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

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

inline LMarker&
LMarker::operator=(const LMarker& that)
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

#endif //#ifndef included_IBTK_LMarker_inl_h
