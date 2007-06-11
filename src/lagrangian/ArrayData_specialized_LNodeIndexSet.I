#ifndef included_ArrayData_specialized_LNodeIndexSet
#define included_ArrayData_specialized_LNodeIndexSet

// Filename: ArrayData_specialized_LNodeIndexSet.I
// Last modified: <10.Jun.2007 22:34:44 boyce@bigboy.nyconnect.com>
// Created on 04 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include <ArrayData.h>
#include <CellData.h>
#include <ibamr/LNodeIndexSet.h>
#include <cassert>

namespace SAMRAI
{
namespace pdat
{
template<>
inline size_t
ArrayData<NDIM,IBAMR::LNodeIndexSet>::getSizeOfData(
    const hier::Box<NDIM>& box,
    const int depth)
{
    return tbox::Arena::align(box.size()*depth*sizeof(IBAMR::LNodeIndexSet));
}// getSizeOfData

template<>
inline bool
ArrayData<NDIM,IBAMR::LNodeIndexSet>::isStandardType()
{
    return false;
}// isStandardType

template<>
inline void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::copy(
    const ArrayData<NDIM,IBAMR::LNodeIndexSet>& src,
    const hier::Box<NDIM>& box)
{
    const hier::Box<NDIM> dst_box(box*d_box*src.d_box);

    // Clear out the destination region.
    fillAll(IBAMR::LNodeIndexSet(),dst_box);

    // Copy items from the source data.
    for (hier::Box<NDIM>::Iterator it(dst_box); it; it++)
    {
        const hier::Index<NDIM>& i_dst = it();
        const hier::Index<NDIM>& i_src = i_dst;
        IBAMR::LNodeIndexSet& dst_node_set = (*this)(i_dst,0);
        const IBAMR::LNodeIndexSet& src_node_set = src(i_src,0);
        if (!src_node_set.empty())
        {
            dst_node_set = src_node_set;
        }
    }
    return;
}// copy

template<>
inline void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::copy(
    const ArrayData<NDIM,IBAMR::LNodeIndexSet>& src,
    const hier::Box<NDIM>& box,
    const hier::IntVector<NDIM>& offset)
{
    const hier::Box<NDIM> dst_box(box*d_box*hier::Box<NDIM>::shift(src.d_box,offset));

    // Clear out the destination region.
    fillAll(IBAMR::LNodeIndexSet(),dst_box);

    // Copy items from the source data.
    for (hier::Box<NDIM>::Iterator it(dst_box); it; it++)
    {
        const hier::Index<NDIM>& i_dst = it();
        const hier::Index<NDIM>& i_src = i_dst-offset;
        IBAMR::LNodeIndexSet& dst_node_set = (*this)(i_dst,0);
        const IBAMR::LNodeIndexSet& src_node_set = src(i_src,0);
        if (!src_node_set.empty())
        {
            dst_node_set.copySourceItem(i_src,offset,src_node_set);
        }
    }
    return;
}// copy

template<>
inline void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::copy(
    const ArrayData<NDIM,IBAMR::LNodeIndexSet>& src,
    const hier::BoxList<NDIM>& boxes,
    const hier::IntVector<NDIM>& offset)
{
    for (hier::BoxList<NDIM>::Iterator b(boxes); b; b++)
    {
        copy(src, b(), offset);
    }
    return;
}// copy

template<>
inline void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::copyDepth(
    int dst_depth,
    const ArrayData<NDIM,IBAMR::LNodeIndexSet>& src,
    int src_depth,
    const hier::Box<NDIM>& box)
{
    TBOX_ERROR("CellData<DIM,TYPE>::copyDepth() is not implemented for TYPE=IBAMR::LNodeIndexSet" << std::endl);
    return;
}// copyDepth

template<>
inline bool
ArrayData<NDIM,IBAMR::LNodeIndexSet>::canEstimateStreamSizeFromBox()
{
    return false;
}// canEstimateStreamSizeFromBox

template<>
inline int
ArrayData<NDIM,IBAMR::LNodeIndexSet>::getDataStreamSize(
    const hier::BoxList<NDIM>& dest_boxes,
    const hier::IntVector<NDIM>& source_offset) const
{
    size_t bytes = tbox::AbstractStream::sizeofInt();
    for (hier::BoxList<NDIM>::Iterator b(dest_boxes); b; b++)
    {
        const hier::Box<NDIM> src_box(d_box*hier::Box<NDIM>::shift(b(), -source_offset));
        for (hier::Box<NDIM>::Iterator it(src_box); it; it++)
        {
            const hier::Index<NDIM>& i = it();
            const IBAMR::LNodeIndexSet& node_set = (*this)(i,0);
            if (!node_set.empty())
            {
                bytes += NDIM*tbox::AbstractStream::sizeofInt() + node_set.getDataStreamSize();
            }
        }
    }
    return bytes;
}// getDataStreamSize

template<>
inline void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::packStream(
    tbox::AbstractStream& stream,
    const hier::Box<NDIM>& dest_box,
    const hier::IntVector<NDIM>& source_offset) const
{
    const hier::Box<NDIM> src_box(d_box*hier::Box<NDIM>::shift(dest_box, -source_offset));

    // Count the number of non-empty index set locations.
    int num_items = 0;
    for (hier::Box<NDIM>::Iterator it(src_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        const IBAMR::LNodeIndexSet& node_set = (*this)(i,0);
        if (!node_set.empty()) ++num_items;
    }
    stream << num_items;

    // Pack the non-empty index sets into the stream.
    for (hier::Box<NDIM>::Iterator it(src_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        const IBAMR::LNodeIndexSet& node_set = (*this)(i,0);
        if (!node_set.empty())
        {
            stream.pack(i, NDIM);
            node_set.packStream(stream);
        }
    }
    return;
}// packStream

template<>
inline void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::packStream(
    tbox::AbstractStream& stream,
    const hier::BoxList<NDIM>& dest_boxes,
    const hier::IntVector<NDIM>& source_offset) const
{
    // Count the number of non-empty index set locations.
    int num_items = 0;
    for (hier::BoxList<NDIM>::Iterator b(dest_boxes); b; b++)
    {
        const hier::Box<NDIM> src_box(d_box*hier::Box<NDIM>::shift(b(), -source_offset));
        for (hier::Box<NDIM>::Iterator it(src_box); it; it++)
        {
            const hier::Index<NDIM>& i = it();
            const IBAMR::LNodeIndexSet& node_set = (*this)(i,0);
            if (!node_set.empty()) ++num_items;
        }
    }
    stream << num_items;

    // Pack the non-empty index sets into the stream.
    for (hier::BoxList<NDIM>::Iterator b(dest_boxes); b; b++)
    {
        const hier::Box<NDIM> src_box(d_box*hier::Box<NDIM>::shift(b(), -source_offset));
        for (hier::Box<NDIM>::Iterator it(src_box); it; it++)
        {
            const hier::Index<NDIM>& i = it();
            const IBAMR::LNodeIndexSet& node_set = (*this)(i,0);
            if (!node_set.empty())
            {
                stream.pack(i, NDIM);
                node_set.packStream(stream);
            }
        }
    }
    return;
}// packStream

template<>
inline void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::unpackStream(
    tbox::AbstractStream& stream,
    const hier::Box<NDIM>& dest_box,
    const hier::IntVector<NDIM>& source_offset)
{
    // Clear out the destination region.
    fillAll(IBAMR::LNodeIndexSet(),dest_box*d_box);

    // Get the number of non-empty index set locations.
    int num_items;
    stream >> num_items;

    // Unpack the individual items.
    hier::Index<NDIM> i;
    for (int k = 0; k < num_items; ++k)
    {
        stream.unpack(i, NDIM);
        (*this)(i,0).unpackStream(stream,source_offset);
    }
    return;
}// unpackStream

template<>
inline void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::unpackStream(
    tbox::AbstractStream& stream,
    const hier::BoxList<NDIM>& dest_boxes,
    const hier::IntVector<NDIM>& source_offset)
{
    // Clear out the destination region.
    for (hier::BoxList<NDIM>::Iterator b(dest_boxes); b; b++)
    {
        const hier::Box<NDIM>& dest_box(b());
        fillAll(IBAMR::LNodeIndexSet(),dest_box*d_box);
    }

    // Get the number of non-empty index set locations.
    int num_items;
    stream >> num_items;

    // Unpack the individual items.
    hier::Index<NDIM> i;
    for (int k = 0; k < num_items; ++k)
    {
        stream.unpack(i, NDIM);
        (*this)(i,0).unpackStream(stream,source_offset);
    }
    return;
}// unpackStream

template <>
inline void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::putSpecializedToDatabase(
    tbox::Pointer<tbox::Database> database)
{
    const hier::Box<NDIM> src_box(d_box);

    // Count the number of non-empty index set locations.
    int num_items = 0;
    for (hier::Box<NDIM>::Iterator it(src_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        const IBAMR::LNodeIndexSet& node_set = (*this)(i,0);
        if (!node_set.empty()) ++num_items;
    }
    database->putInteger("num_items",num_items);

    // Pack the non-empty index sets into the database.
    int k = 0;
    for (hier::Box<NDIM>::Iterator it(src_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        const IBAMR::LNodeIndexSet& node_set = (*this)(i,0);
        if (!node_set.empty())
        {
            std::ostringstream idx_stream;
            idx_stream << "idx_" << k;
            database->putIntegerArray(idx_stream.str(),i,NDIM);

            std::ostringstream db_stream;
            db_stream << "db_" << k;
            SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> sub_db = database->putDatabase(db_stream.str());
            node_set.putToDatabase(sub_db);

            ++k;
        }
    }
    return;
}// putSpecializedToDatabase

template <>
inline void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::getSpecializedFromDatabase(
    tbox::Pointer<tbox::Database> database)
{
    const hier::Box<NDIM> dst_box(d_box);

    // Clear out the destination region.
    fillAll(IBAMR::LNodeIndexSet(),dst_box);

    // Get the number of non-empty index set locations.
    const int num_items = database->getInteger("num_items");

    // Unpack the individual items.
    hier::Index<NDIM> i;
    for (int k = 0; k < num_items; ++k)
    {
        std::ostringstream idx_stream;
        idx_stream << "idx_" << k;
        database->getIntegerArray(idx_stream.str(),i,NDIM);

        std::ostringstream db_stream;
        db_stream << "db_" << k;
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> sub_db = database->getDatabase(db_stream.str());
        IBAMR::LNodeIndexSet& node_set = (*this)(i,0);
        node_set.getFromDatabase(sub_db);
    }
    return;
}// getSpecializedFromDatabase

template <>
inline void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::undefineData()
{
    fillAll(IBAMR::LNodeIndexSet());
    return;
}// undefineData

template <>
inline void
CellData<NDIM,IBAMR::LNodeIndexSet>::print(
    const hier::Box<NDIM>& box,
    std::ostream& os,
    const int prec) const
{
    TBOX_ERROR("CellData<DIM,TYPE>::print() is not implemented for TYPE=IBAMR::LNodeIndexSet" << std::endl);
    return;
}// print

template<>
inline void
CellData<NDIM,IBAMR::LNodeIndexSet>::print(
    const hier::Box<NDIM>& box,
    const int d,
    std::ostream& os,
    int prec) const
{
    TBOX_ERROR("CellData<DIM,TYPE>::print() is not implemented for TYPE=IBAMR::LNodeIndexSet" << std::endl);
    return;
}// print
}
}

#endif //#ifndef included_ArrayData_specialized_LNodeIndexSet
