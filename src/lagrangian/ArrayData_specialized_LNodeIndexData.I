// Filename: ArrayData_specialized_LNodeIndexSet.I
// Last modified: <04.Jun.2007 16:26:29 griffith@box221.cims.nyu.edu>
// Created on 04 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include <ArrayData.C>
#include <ibamr/LNodeIndexSet.h>
#include <cassert>

namespace SAMRAI
{
namespace pdat
{
template<>
size_t
ArrayData<NDIM,IBAMR::LNodeIndexSet>::getSizeOfData(
    const hier::Box<NDIM>& box,
    const int depth)
{
    return tbox::Arena::align(box.size()*depth*sizeof(IBAMR::LNodeIndexSet));
}// getSizeOfData

template<>
bool
ArrayData<NDIM,IBAMR::LNodeIndexSet>::isStandardType()
{
    return false;
}// isStandardType

template<>
bool
ArrayData<NDIM,IBAMR::LNodeIndexSet>::canEstimateStreamSizeFromBox()
{
    return false;
}// canEstimateStreamSizeFromBox

template<>
int
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
void
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
void
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
void
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
void
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
void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::putSpecializedToDatabase(
    tbox::Pointer<tbox::Database> database)
{
    assert(false);
    return;
}// putSpecializedToDatabase

template <>
void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::getSpecializedFromDatabase(
    tbox::Pointer<tbox::Database> database)
{
    assert(false);
    return;
}// getSpecializedFromDatabase

template <>
void
ArrayData<NDIM,IBAMR::LNodeIndexSet>::undefineData()
{
    fillAll(IBAMR::LNodeIndexSet());
    return;
}// undefineData

template <>
void
CellData<NDIM,IBAMR::LNodeIndexSet>::print(
    const hier::Box<NDIM>& box,
    std::ostream& os,
    const int prec) const
{
    TBOX_ERROR("CellData<DIM,TYPE>::print() is not implemented for TYPE=IBAMR::LNodeIndexSet" << std::endl);
    return;
}// print

template<>
void
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
