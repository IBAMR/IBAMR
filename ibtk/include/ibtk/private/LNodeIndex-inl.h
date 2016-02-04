// Filename: LNodeIndex-inl.h
// Created on 28 Feb 2004 by Boyce Griffith
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

#ifndef included_LNodeIndex_inl_h
#define included_LNodeIndex_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LNodeIndex.h"
#include "tbox/AbstractStream.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline LNodeIndex::LNodeIndex(const int lagrangian_nidx,
                              const int global_petsc_nidx,
                              const int local_petsc_nidx,
                              const SAMRAI::hier::IntVector<NDIM>& periodic_offset,
                              const Vector& periodic_displacement)
    : d_lagrangian_nidx(lagrangian_nidx),
      d_global_petsc_nidx(global_petsc_nidx),
      d_local_petsc_nidx(local_petsc_nidx),
      d_offset(periodic_offset),
      d_displacement(periodic_displacement)
{
    // intentionally blank
    return;
} // LNodeIndex

inline LNodeIndex::LNodeIndex(const LNodeIndex& from)
    : d_lagrangian_nidx(from.d_lagrangian_nidx),
      d_global_petsc_nidx(from.d_global_petsc_nidx),
      d_local_petsc_nidx(from.d_local_petsc_nidx),
      d_offset(from.d_offset),
      d_displacement(from.d_displacement)
{
    // intentionally blank
    return;
} // LNodeIndex

inline LNodeIndex::LNodeIndex(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& offset)
    : d_lagrangian_nidx(-1),
      d_global_petsc_nidx(-1),
      d_local_petsc_nidx(-1),
      d_offset(0),
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

inline LNodeIndex& LNodeIndex::operator=(const LNodeIndex& that)
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
LNodeIndex::getPeriodicOffset() const
{
    return d_offset;
} // getPeriodicOffset

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
    return (3 + NDIM) * SAMRAI::tbox::AbstractStream::sizeofInt() + NDIM * SAMRAI::tbox::AbstractStream::sizeofDouble();
} // getDataStreamSize

inline void
LNodeIndex::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_lagrangian_nidx, 1);
    stream.pack(&d_global_petsc_nidx, 1);
    stream.pack(&d_local_petsc_nidx, 1);
    stream.pack(d_offset, NDIM);
    stream.pack(d_displacement.data(), NDIM);
    return;
} // packStream

inline void
LNodeIndex::unpackStream(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& /*offset*/)
{
    stream.unpack(&d_lagrangian_nidx, 1);
    stream.unpack(&d_global_petsc_nidx, 1);
    stream.unpack(&d_local_petsc_nidx, 1);
    stream.unpack(d_offset, NDIM);
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
    d_offset = that.d_offset;
    d_displacement = that.d_displacement;
    return;
} // assignThatToThis

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndex_inl_h
