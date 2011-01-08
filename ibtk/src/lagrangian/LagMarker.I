// Filename: LagMarker.I
// Created on 12 Sep 2007 by Boyce Griffith
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

inline
LagMarker::LagMarker()
    : d_X(),
      d_U(),
      d_idx(),
      d_offset(0)
{
    // intentionally blank
    return;
}// LagMarker

inline
LagMarker::LagMarker(
    const LagMarker& from)
    : d_X(from.d_X),
      d_U(from.d_U),
      d_idx(from.d_idx),
      d_offset(from.d_offset)
{
    // intentionally blank
    return;
}// LagMarker

inline
LagMarker::~LagMarker()
{
    // intentionally blank
    return;
}// ~LagMarker

inline LagMarker&
LagMarker::operator=(
    const LagMarker& that)
{
    if (this == &that) return *this;  // check for self-assignment

    d_X      = that.d_X;
    d_U      = that.d_U;
    d_idx    = that.d_idx;
    d_offset = that.d_offset;
    return *this;
}// operator=

inline int
LagMarker::getNumberOfMarkers() const
{
    return d_X.size()/NDIM;
}// getNumberOfMarkers

inline const std::vector<double>&
LagMarker::getPositions() const
{
    return d_X;
}// getPositions

inline std::vector<double>&
LagMarker::getPositions()
{
    return d_X;
}// getPositions

inline void
LagMarker::setPositions(
    const std::vector<double>& X)
{
    d_X = X;
    return;
}// setPositions

inline const std::vector<double>&
LagMarker::getVelocities() const
{
    return d_U;
}// getVelocities

inline std::vector<double>&
LagMarker::getVelocities()
{
    return d_U;
}// getVelocities

inline void
LagMarker::setVelocities(
    const std::vector<double>& U)
{
    d_U = U;
    return;
}// setVelocities

inline const std::vector<int>&
LagMarker::getIndices() const
{
    return d_idx;
}// getIndices

inline std::vector<int>&
LagMarker::getIndices()
{
    return d_idx;
}// getIndices

inline void
LagMarker::setIndices(
    const std::vector<int>& idx)
{
    d_idx = idx;
    return;
}// setIndices

inline void
LagMarker::addMarker(
    const LagMarker& that)
{
    this->d_X.insert(this->d_X.end(),that.d_X.begin(),that.d_X.end());
    this->d_U.insert(this->d_U.end(),that.d_U.begin(),that.d_U.end());
    this->d_idx.insert(this->d_idx.end(),that.d_idx.begin(),that.d_idx.end());
    return;
}// addMarker

inline const SAMRAI::hier::IntVector<NDIM>&
LagMarker::getPeriodicOffset() const
{
    return d_offset;
}// getPeriodicOffset

inline void
LagMarker::setPeriodicOffset(
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    d_offset = offset;
    return;
}// setPeriodicOffset

inline void
LagMarker::copySourceItem(
    const SAMRAI::hier::Index<NDIM>& src_index,
    const SAMRAI::hier::IntVector<NDIM>& src_offset,
    const LagMarker& src_item)
{
    d_X      = src_item.d_X;
    d_U      = src_item.d_U;
    d_idx    = src_item.d_idx;
    d_offset = src_offset;
    return;
}// copySourceItem

inline size_t
LagMarker::getDataStreamSize() const
{
    return (SAMRAI::tbox::AbstractStream::sizeofDouble()*(d_X.size()+d_U.size())+
            SAMRAI::tbox::AbstractStream::sizeofInt()*(d_idx.size()+3));
}// getDataStreamSize

inline void
LagMarker::packStream(
    SAMRAI::tbox::AbstractStream& stream) const
{
    stream << int(d_X.size());
    stream.pack(&d_X[0],d_X.size());

    stream << int(d_U.size());
    stream.pack(&d_U[0],d_U.size());

    stream << int(d_idx.size());
    stream.pack(&d_idx[0],d_idx.size());
    return;
}// packStream

inline void
LagMarker::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int X_size;
    stream >> X_size;
    d_X.resize(X_size);
    stream.unpack(&d_X[0],X_size);

    int U_size;
    stream >> U_size;
    d_U.resize(U_size);
    stream.unpack(&d_U[0],U_size);

    int idx_size;
    stream >> idx_size;
    d_idx.resize(idx_size);
    stream.unpack(&d_idx[0],idx_size);
    return;
}// unpackStream

inline void
LagMarker::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database) const
{
    database->putInteger("d_X.size()", d_X.size());
    database->putDoubleArray("d_X", &d_X[0], d_X.size());

    database->putInteger("d_U.size()", d_U.size());
    database->putDoubleArray("d_U", &d_U[0], d_U.size());

    database->putInteger("d_idx.size()", d_idx.size());
    database->putIntegerArray("d_idx", &d_idx[0], d_idx.size());

    database->putIntegerArray("d_offset", d_offset, NDIM);
    return;
}// putToDatabase

inline void
LagMarker::getFromDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database)
{
    const int X_size = database->getInteger("d_X.size()");
    d_X.resize(X_size);
    database->getDoubleArray("d_X", &d_X[0], d_X.size());

    const int U_size = database->getInteger("d_U.size()");
    d_U.resize(U_size);
    database->getDoubleArray("d_U", &d_U[0], d_U.size());

    const int idx_size = database->getInteger("d_idx.size()");
    d_idx.resize(idx_size);
    database->getIntegerArray("d_idx", &d_idx[0], d_idx.size());

    database->getIntegerArray("d_offset", d_offset, NDIM);
    return;
}// getFromDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
