// Filename: LNodeIndexSet.C
// Last modified: <04.Jun.2007 16:16:44 griffith@box221.cims.nyu.edu>
// Created on 29 Feb 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "LNodeIndexSet.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <ArrayData.C>
#include <CellData.C>
#include <CellDataFactory.C>
#include <CellVariable.C>
#include <IndexData.C>
#include <IndexDataFactory.C>
#include <IndexVariable.C>
#include <tbox/Array.C>
#include <tbox/List.C>
#include <tbox/Pointer.C>

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

template class SAMRAI::tbox::Pointer<IBAMR::LNodeIndexSet>;

#if (NDIM == 1)
template class SAMRAI::pdat::ArrayData<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellData<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellDataFactory<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellVariable<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexData<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataFactory<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexIterator<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexVariable<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<1,IBAMR::LNodeIndexSet> >;
#endif

#if (NDIM == 2)
template class SAMRAI::pdat::ArrayData<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellData<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellDataFactory<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellVariable<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexData<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataFactory<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexIterator<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexVariable<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<2,IBAMR::LNodeIndexSet> >;
#endif

#if (NDIM == 3)
template class SAMRAI::pdat::ArrayData<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellData<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellDataFactory<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellVariable<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexData<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataFactory<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexIterator<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexVariable<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<3,IBAMR::LNodeIndexSet> >;
#endif

//////////////////////////////////////////////////////////////////////////////
