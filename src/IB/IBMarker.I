// Filename: IBMarker.I
// Last modified: <13.Sep.2007 01:20:24 griffith@box221.cims.nyu.edu>
// Created on 12 Sep 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
IBMarker::IBMarker()
    : d_X(),
      d_U(),
      d_idx(),
      d_offset(0)
{
    // intentionally blank
    return;
}// IBMarker

inline
IBMarker::IBMarker(
    const IBMarker& from)
    : d_X(from.d_X),
      d_U(from.d_U),
      d_idx(from.d_idx),
      d_offset(from.d_offset)
{
    // intentionally blank
    return;
}// IBMarker

inline
IBMarker::~IBMarker()
{
    // intentionally blank
    return;
}// ~IBMarker

inline IBMarker&
IBMarker::operator=(
    const IBMarker& that)
{
    if (this == &that) return *this;  // check for self-assignment

    d_X      = that.d_X;
    d_U      = that.d_U;
    d_idx    = that.d_idx;
    d_offset = that.d_offset;
    return *this;
}// operator=

inline int
IBMarker::getNumberOfMarkers() const
{
    return d_X.size()/NDIM;
}// getNumberOfMarkers

inline const std::vector<double>&
IBMarker::getPositions() const
{
    return d_X;
}// getPositions

inline std::vector<double>&
IBMarker::getPositions()
{
    return d_X;
}// getPositions

inline void
IBMarker::setPositions(
    const std::vector<double>& X)
{
    d_X = X;
    return;
}// setPositions

inline const std::vector<double>&
IBMarker::getVelocities() const
{
    return d_U;
}// getVelocities

inline std::vector<double>&
IBMarker::getVelocities()
{
    return d_U;
}// getVelocities

inline void
IBMarker::setVelocities(
    const std::vector<double>& U)
{
    d_U = U;
    return;
}// setVelocities

inline const std::vector<int>&
IBMarker::getIndices() const
{
    return d_idx;
}// getIndices

inline std::vector<int>&
IBMarker::getIndices()
{
    return d_idx;
}// getIndices

inline void
IBMarker::setIndices(
    const std::vector<int>& idx)
{
    d_idx = idx;
    return;
}// setIndices

inline const SAMRAI::hier::IntVector<NDIM>&
IBMarker::getPeriodicOffset() const
{
    return d_offset;
}// getPeriodicOffset

inline void
IBMarker::setPeriodicOffset(
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    d_offset = offset;
    return;
}// setPeriodicOffset

inline void
IBMarker::copySourceItem(
    const SAMRAI::hier::Index<NDIM>& src_index,
    const SAMRAI::hier::IntVector<NDIM>& src_offset,
    const IBMarker& src_item)
{
    d_X      = src_item.d_X;
    d_U      = src_item.d_U;
    d_idx    = src_item.d_idx;
    d_offset = src_offset;
    return;
}// copySourceItem

inline size_t
IBMarker::getDataStreamSize() const
{
    return (SAMRAI::tbox::AbstractStream::sizeofDouble()*2*d_X.size() +
            SAMRAI::tbox::AbstractStream::sizeofInt()*(d_idx.size()+3));
}// getDataStreamSize

inline void
IBMarker::packStream(
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
IBMarker::unpackStream(
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
IBMarker::putToDatabase(
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
IBMarker::getFromDatabase(
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

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
