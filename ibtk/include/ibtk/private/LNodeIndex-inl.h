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

#ifndef included_IBTK_LNodeIndex_inl_h
#define included_IBTK_LNodeIndex_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LNodeIndex.h"

#include "tbox/AbstractStream.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline LNodeIndex::LNodeIndex(const int lagrangian_nidx,
                              const int global_petsc_nidx,
                              const int local_petsc_nidx,
                              const SAMRAI::hier::IntVector<NDIM>& initial_periodic_offset,
                              const SAMRAI::hier::IntVector<NDIM>& current_periodic_offset,
                              const Vector& initial_periodic_displacement,
                              const Vector& current_periodic_displacement)
    : d_lagrangian_nidx(lagrangian_nidx),
      d_global_petsc_nidx(global_petsc_nidx),
      d_local_petsc_nidx(local_petsc_nidx),
      d_offset_0(initial_periodic_offset),
      d_offset(current_periodic_offset),
      d_displacement_0(initial_periodic_displacement),
      d_displacement(current_periodic_displacement)
{
    // intentionally blank
    return;
} // LNodeIndex

inline LNodeIndex::LNodeIndex(const LNodeIndex& from)
    : d_lagrangian_nidx(from.d_lagrangian_nidx),
      d_global_petsc_nidx(from.d_global_petsc_nidx),
      d_local_petsc_nidx(from.d_local_petsc_nidx),
      d_offset_0(from.d_offset_0),
      d_offset(from.d_offset),
      d_displacement_0(from.d_displacement_0),
      d_displacement(from.d_displacement)
{
    // intentionally blank
    return;
} // LNodeIndex

inline LNodeIndex::LNodeIndex(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& offset)
    : d_lagrangian_nidx(-1),
      d_global_petsc_nidx(-1),
      d_local_petsc_nidx(-1),
      d_offset_0(0),
      d_offset(0),
      d_displacement_0(Vector::Zero()),
      d_displacement(Vector::Zero())
{
    unpackStream(stream, offset);
    return;
} // LNodeIndex

inline LNodeIndex::~LNodeIndex()
{
    // intentionally blank
    return;
} // ~LNodeIndex

inline LNodeIndex&
LNodeIndex::operator=(const LNodeIndex& that)
{
    if (this != &that)
    {
        assignThatToThis(that);
    }
    return *this;
} // operator=

inline int
LNodeIndex::getLagrangianIndex() const
{
    return d_lagrangian_nidx;
} // getLagrangianIndex

inline void
LNodeIndex::setLagrangianIndex(const int lagrangian_nidx)
{
    d_lagrangian_nidx = lagrangian_nidx;
    return;
} // setLagrangianIndex

inline int
LNodeIndex::getGlobalPETScIndex() const
{
    return d_global_petsc_nidx;
} // getGlobalPETScIndex

inline void
LNodeIndex::setGlobalPETScIndex(const int global_petsc_nidx)
{
    d_global_petsc_nidx = global_petsc_nidx;
    return;
} // setGlobalPETScIndex

inline int
LNodeIndex::getLocalPETScIndex() const
{
    return d_local_petsc_nidx;
} // getLocalPETScIndex

inline void
LNodeIndex::setLocalPETScIndex(const int local_petsc_nidx)
{
    d_local_petsc_nidx = local_petsc_nidx;
    return;
} // setLocalPETScIndex

inline void
LNodeIndex::registerPeriodicShift(const SAMRAI::hier::IntVector<NDIM>& offset, const Vector& displacement)
{
    d_offset += offset;
    d_displacement += displacement;
    return;
} // registerPeriodicShift

inline const SAMRAI::hier::IntVector<NDIM>&
LNodeIndex::getInitialPeriodicOffset() const
{
    return d_offset_0;
} // getInitialPeriodicOffset

inline const SAMRAI::hier::IntVector<NDIM>&
LNodeIndex::getPeriodicOffset() const
{
    return d_offset;
} // getPeriodicOffset

inline const Vector&
LNodeIndex::getInitialPeriodicDisplacement() const
{
    return d_displacement_0;
} // getInitialPeriodicDisplacement

inline const Vector&
LNodeIndex::getPeriodicDisplacement() const
{
    return d_displacement;
} // getPeriodicDisplacement

inline void
LNodeIndex::copySourceItem(const SAMRAI::hier::Index<NDIM>& /*src_index*/,
                           const SAMRAI::hier::IntVector<NDIM>& /*src_offset*/,
                           const LNodeIndex& src_item)
{
    assignThatToThis(src_item);
    return;
} // copySourceItem

inline size_t
LNodeIndex::getDataStreamSize() const
{
    return (3 + 2 * NDIM) * SAMRAI::tbox::AbstractStream::sizeofInt() +
           (2 * NDIM) * SAMRAI::tbox::AbstractStream::sizeofDouble();
} // getDataStreamSize

inline void
LNodeIndex::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_lagrangian_nidx, 1);
    stream.pack(&d_global_petsc_nidx, 1);
    stream.pack(&d_local_petsc_nidx, 1);
    stream.pack(d_offset_0, NDIM);
    stream.pack(d_offset, NDIM);
    stream.pack(d_displacement_0.data(), NDIM);
    stream.pack(d_displacement.data(), NDIM);
    return;
} // packStream

inline void
LNodeIndex::unpackStream(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& /*offset*/)
{
    stream.unpack(&d_lagrangian_nidx, 1);
    stream.unpack(&d_global_petsc_nidx, 1);
    stream.unpack(&d_local_petsc_nidx, 1);
    stream.unpack(d_offset_0, NDIM);
    stream.unpack(d_offset, NDIM);
    stream.unpack(d_displacement_0.data(), NDIM);
    stream.unpack(d_displacement.data(), NDIM);
    return;
} // unpackStream

/////////////////////////////// PRIVATE //////////////////////////////////////

inline void
LNodeIndex::assignThatToThis(const LNodeIndex& that)
{
    d_lagrangian_nidx = that.d_lagrangian_nidx;
    d_global_petsc_nidx = that.d_global_petsc_nidx;
    d_local_petsc_nidx = that.d_local_petsc_nidx;
    d_offset_0 = that.d_offset_0;
    d_offset = that.d_offset;
    d_displacement_0 = that.d_displacement_0;
    d_displacement = that.d_displacement;
    return;
} // assignThatToThis

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LNodeIndex_inl_h
