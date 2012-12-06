// Filename: VecCellData.C
// Created on 09 Apr 2010 by Boyce Griffith
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

#include "VecCellData.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/namespaces.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int VEC_CELL_DATA_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

template<class TYPE>
size_t
VecCellData<TYPE>::getSizeOfData(
    const Box<NDIM>& box,
    const unsigned int depth,
    const IntVector<NDIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(depth > 0);
#endif
    const Box<NDIM> ghost_box = Box<NDIM>::grow(box, ghosts);
    return sizeof(TYPE)*depth*ghost_box.size();
}// getSizeOfData

template<class TYPE>
VecCellData<TYPE>::VecCellData(
    const Box<NDIM>& box,
    const unsigned int depth,
    const IntVector<NDIM>& ghosts)
    : PatchData<NDIM>(box, ghosts),
      d_depth(depth),
#if   (NDIM == 2)
      d_data(blitz::Range(0,depth-1),
             blitz::Range(box.lower()(0)-ghosts(0),box.upper()(0)+ghosts(0)),
             blitz::Range(box.lower()(1)-ghosts(1),box.upper()(1)+ghosts(1)),
             blitz::ColumnMajorArray<NDIM+1>())
#elif (NDIM == 3)
      d_data(blitz::Range(0,depth-1),
             blitz::Range(box.lower()(0)-ghosts(0),box.upper()(0)+ghosts(0)),
             blitz::Range(box.lower()(1)-ghosts(1),box.upper()(1)+ghosts(1)),
             blitz::Range(box.lower()(2)-ghosts(2),box.upper()(2)+ghosts(2)),
             blitz::ColumnMajorArray<NDIM+1>())
#endif
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(depth > 0);
    TBOX_ASSERT(ghosts.min() >= 0);
#endif
#if (NDIM != 2 && NDIM != 3)
    TBOX_ERROR("VecCellData implemented only for NDIM=2,3\n");
#endif
#ifdef DEBUG_INITIALIZE_UNDEFINED
    fillAll(MathUtilities<TYPE>::getSignalingNaN());
#endif
    return;
}// VecCellData

template<class TYPE>
VecCellData<TYPE>::~VecCellData()
{
    // intentionally blank
    return;
}// ~VecCellData

template<class TYPE>
void
VecCellData<TYPE>::copy(
    const PatchData<NDIM>& src)
{
    const VecCellData<TYPE>* const t_vcd_src = dynamic_cast<const VecCellData<TYPE>*>(&src);
    const CellData<NDIM,TYPE>* const t_cd_src = dynamic_cast<const CellData<NDIM,TYPE>*>(&src);
    if (t_vcd_src != NULL)
    {
        const Box<NDIM> box = getGhostBox() * t_vcd_src->getGhostBox();
        copyOnBox(*t_vcd_src, box);
    }
    else if (t_cd_src != NULL)
    {
        const Box<NDIM> box = getGhostBox() * t_cd_src->getGhostBox();
        copyOnBox(*t_cd_src, box);
    }
    else
    {
        src.copy2(*this);
    }
    return;
}// copy

template<class TYPE>
void
VecCellData<TYPE>::copy2(
    PatchData<NDIM>& dst) const
{
    VecCellData<TYPE>* const t_vcd_dst = dynamic_cast<VecCellData<TYPE>*>(&dst);
    CellData<NDIM,TYPE>* const t_cd_dst = dynamic_cast<CellData<NDIM,TYPE>*>(&dst);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(t_vcd_dst != NULL || t_cd_dst != NULL);
#endif
    if (t_vcd_dst != NULL)
    {
        const Box<NDIM> box = getGhostBox() * t_vcd_dst->getGhostBox();
        t_vcd_dst->copyOnBox(*this, box);
    }
    else if (t_cd_dst != NULL)
    {
        const Box<NDIM> box = getGhostBox() * t_cd_dst->getGhostBox();
        copy2OnBox(*t_cd_dst, box);
    }
    return;
}// copy2

template<class TYPE>
void
VecCellData<TYPE>::copy(
    const PatchData<NDIM>& src,
    const BoxOverlap<NDIM>& overlap)
{
    const VecCellData<TYPE>* const t_src = dynamic_cast<const VecCellData<TYPE>*>(&src);
    const CellOverlap<NDIM>* const t_overlap = dynamic_cast<const CellOverlap<NDIM>*>(&overlap);
    if ((t_src == NULL) || (t_overlap == NULL))
    {
        src.copy2(*this, overlap);
    }
    else
    {
        const BoxList<NDIM>& dst_box_list = t_overlap->getDestinationBoxList();
        const IntVector<NDIM>& src_offset = t_overlap->getSourceOffset();
        for (BoxList<NDIM>::Iterator b(dst_box_list); b; b++)
        {
            const Box<NDIM> box = b() * getGhostBox() * Box<NDIM>::shift(t_src->getGhostBox(),src_offset);
            copyOnBox(*t_src, box, src_offset);
        }
    }
    return;
}// copy

template<class TYPE>
void VecCellData<TYPE>::copy2(
    PatchData<NDIM>& dst, const BoxOverlap<NDIM>& overlap) const
{
    VecCellData<TYPE>* const t_dst = dynamic_cast<VecCellData<TYPE>*>(&dst);
    const CellOverlap<NDIM>* const t_overlap = dynamic_cast<const CellOverlap<NDIM>*>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(t_dst != NULL);
    TBOX_ASSERT(t_overlap != NULL);
#endif
    const BoxList<NDIM>& dst_box_list = t_overlap->getDestinationBoxList();
    const IntVector<NDIM>& src_offset = t_overlap->getSourceOffset();
    for (BoxList<NDIM>::Iterator b(dst_box_list); b; b++)
    {
        const Box<NDIM> box = b() * t_dst->getGhostBox() * Box<NDIM>::shift(getGhostBox(),src_offset);
        t_dst->copyOnBox(*this, box, src_offset);
    }
    return;
}// copy2

template<>
bool
VecCellData<bool>::canEstimateStreamSizeFromBox() const
{
    return true;
}// canEstimateStreamSizeFromBox

template<>
bool
VecCellData<char>::canEstimateStreamSizeFromBox() const
{
    return true;
}// canEstimateStreamSizeFromBox

template<>
bool
VecCellData<double>::canEstimateStreamSizeFromBox() const
{
    return true;
}// canEstimateStreamSizeFromBox

template<>
bool
VecCellData<float>::canEstimateStreamSizeFromBox() const
{
    return true;
}// canEstimateStreamSizeFromBox

template<>
bool
VecCellData<int>::canEstimateStreamSizeFromBox() const
{
    return true;
}// canEstimateStreamSizeFromBox

template<>
bool
VecCellData<dcomplex>::canEstimateStreamSizeFromBox() const
{
    return true;
}// canEstimateStreamSizeFromBox

template<>
int
VecCellData<bool>::getDataStreamSize(
    const BoxOverlap<NDIM>& overlap) const
{
    const CellOverlap<NDIM>* const t_overlap = dynamic_cast<const CellOverlap<NDIM>*>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(t_overlap != NULL);
#endif
    int size = 0;
    const BoxList<NDIM>& dst_box_list = t_overlap->getDestinationBoxList();
    for (BoxList<NDIM>::Iterator b(dst_box_list); b; b++)
    {
        const Box<NDIM>& box = b();
        size += sizeof(bool)*d_depth*box.size();
    }
    return size;
}// getDataStreamSize

template<>
int
VecCellData<char>::getDataStreamSize(
    const BoxOverlap<NDIM>& overlap) const
{
    const CellOverlap<NDIM>* const t_overlap = dynamic_cast<const CellOverlap<NDIM>*>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(t_overlap != NULL);
#endif
    int size = 0;
    const BoxList<NDIM>& dst_box_list = t_overlap->getDestinationBoxList();
    for (BoxList<NDIM>::Iterator b(dst_box_list); b; b++)
    {
        const Box<NDIM>& box = b();
        size += sizeof(char)*d_depth*box.size();
    }
    return size;
}// getDataStreamSize

template<>
int
VecCellData<double>::getDataStreamSize(
    const BoxOverlap<NDIM>& overlap) const
{
    const CellOverlap<NDIM>* const t_overlap = dynamic_cast<const CellOverlap<NDIM>*>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(t_overlap != NULL);
#endif
    int size = 0;
    const BoxList<NDIM>& dst_box_list = t_overlap->getDestinationBoxList();
    for (BoxList<NDIM>::Iterator b(dst_box_list); b; b++)
    {
        const Box<NDIM>& box = b();
        size += sizeof(double)*d_depth*box.size();
    }
    return size;
}// getDataStreamSize

template<>
int
VecCellData<float>::getDataStreamSize(
    const BoxOverlap<NDIM>& overlap) const
{
    const CellOverlap<NDIM>* const t_overlap = dynamic_cast<const CellOverlap<NDIM>*>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(t_overlap != NULL);
#endif
    int size = 0;
    const BoxList<NDIM>& dst_box_list = t_overlap->getDestinationBoxList();
    for (BoxList<NDIM>::Iterator b(dst_box_list); b; b++)
    {
        const Box<NDIM>& box = b();
        size += sizeof(float)*d_depth*box.size();
    }
    return size;
}// getDataStreamSize

template<>
int
VecCellData<dcomplex>::getDataStreamSize(
    const BoxOverlap<NDIM>& overlap) const
{
    const CellOverlap<NDIM>* const t_overlap = dynamic_cast<const CellOverlap<NDIM>*>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(t_overlap != NULL);
#endif
    int size = 0;
    const BoxList<NDIM>& dst_box_list = t_overlap->getDestinationBoxList();
    for (BoxList<NDIM>::Iterator b(dst_box_list); b; b++)
    {
        const Box<NDIM>& box = b();
        size += sizeof(dcomplex)*d_depth*box.size();
    }
    return size;
}// getDataStreamSize

template<class TYPE>
void
VecCellData<TYPE>::packStream(
    AbstractStream& stream,
    const BoxOverlap<NDIM>& overlap) const
{
    const CellOverlap<NDIM>* const t_overlap = dynamic_cast<const CellOverlap<NDIM>*>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(t_overlap != NULL);
#endif
    const BoxList<NDIM>& dst_box_list = t_overlap->getDestinationBoxList();
    const IntVector<NDIM>& src_offset = t_overlap->getSourceOffset();
    for (BoxList<NDIM>::Iterator b(dst_box_list); b; b++)
    {
        const Box<NDIM> box = Box<NDIM>::shift(b(),-src_offset);
        const int size = d_depth*box.size();
        std::vector<TYPE> buffer;
        buffer.reserve(size);
        for (Iterator i(box); i; i++)
        {
            const Index<NDIM>& idx = i();
            for (unsigned int d = 0; d < d_depth; ++d)
            {
                buffer.push_back((*this)(idx,d));
            }
        }
        if (!buffer.empty()) stream.pack(&buffer[0], size);
    }
    return;
}// packStream

template<class TYPE>
void
VecCellData<TYPE>::unpackStream(
    AbstractStream& stream,
    const BoxOverlap<NDIM>& overlap)
{
    const CellOverlap<NDIM>* const t_overlap = dynamic_cast<const CellOverlap<NDIM>*>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(t_overlap != NULL);
#endif
    const BoxList<NDIM>& dst_box_list = t_overlap->getDestinationBoxList();
    for (BoxList<NDIM>::Iterator b(dst_box_list); b; b++)
    {
        const Box<NDIM>& box = b();
        const int size = d_depth*box.size();
        std::vector<TYPE> buffer(size);
        if (!buffer.empty()) stream.unpack(&buffer[0], size);
        int k = 0;
        for (Iterator i(box); i; i++)
        {
            const Index<NDIM>& idx = i();
            for (unsigned int d = 0; d < d_depth; ++d, ++k)
            {
                (*this)(idx,d) = buffer[k];
            }
        }
    }
    return;
}// unpackStream

template<class TYPE>
void
VecCellData<TYPE>::print(
    const Box<NDIM>& box,
    std::ostream& os,
    const int prec) const
{
    for (unsigned int d = 0; d < d_depth; ++d)
    {
        os << "Array depth = " << d << std::endl;
        print(box, d, os, prec);
    }
    return;
}// print

template<class TYPE>
void
VecCellData<TYPE>::print(
    const Box<NDIM>& box,
    const unsigned int depth,
    std::ostream& os,
    const int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(depth < d_depth);
#endif
    os.precision(prec);
    for (Iterator i(box); i; i++)
    {
        const Index<NDIM>& idx = i();
        os << "array" << i() << " = " << (*this)(idx,depth) << std::endl << std::flush;
        os << std::flush;
    }
    return;
}// print

template<class TYPE>
void VecCellData<TYPE>::getSpecializedFromDatabase(
    Pointer<Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(database);
#endif
    const int ver = database->getInteger("VEC_CELL_DATA_VERSION");
    if (ver != VEC_CELL_DATA_VERSION)
    {
        TBOX_ERROR("VecCellData<NDIM>::getSpecializedFromDatabase error...\n"
                   << "Restart file version different than class version" << std::endl);
    }
    d_depth = database->getInteger("d_depth");
    const Box<NDIM>& ghost_box = getGhostBox();
    const Index<NDIM>& lower = ghost_box.lower();
    const Index<NDIM>& upper = ghost_box.upper();
#if   (NDIM == 2)
    d_data.resize(blitz::Range(0,d_depth-1),
                  blitz::Range(lower(0),upper(0)),
                  blitz::Range(lower(1),upper(1)));
#elif (NDIM == 3)
    d_data.resize(blitz::Range(0,d_depth-1),
                  blitz::Range(lower(0),upper(0)),
                  blitz::Range(lower(1),upper(1)),
                  blitz::Range(lower(2),upper(2)));
#else
    TBOX_ERROR("VecCellData implemented only for NDIM=2,3\n");
#endif
    getArrayFromDatabase(database);
    return;
}// getSpecializedFromDatabase

template<class TYPE>
void
VecCellData<TYPE>::putSpecializedToDatabase(
    Pointer<Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(database);
#endif
    database->putInteger("VEC_CELL_DATA_VERSION", VEC_CELL_DATA_VERSION);
    database->putInteger("d_depth", d_depth);
    putArrayToDatabase(database);
    return;
}// putSpecializedToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

template<>
void
VecCellData<float>::getArrayFromDatabase(
    Pointer<Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(database);
#endif
    database->getFloatArray("d_data", d_data.data(), d_data.size());
    return;
}// getArrayFromDatabase

template<>
void
VecCellData<double>::getArrayFromDatabase(
    Pointer<Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(database);
#endif
    database->getDoubleArray("d_data", d_data.data(), d_data.size());
    return;
}// getArrayFromDatabase

template<>
void
VecCellData<float>::putArrayToDatabase(
    Pointer<Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(database);
#endif
    database->putFloatArray("d_data", d_data.data(), d_data.size());
    return;
}// putArrayToDatabase

template<>
void
VecCellData<double>::putArrayToDatabase(
    Pointer<Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(database);
#endif
    database->putDoubleArray("d_data", d_data.data(), d_data.size());
    return;
}// putArrayToDatabase

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class IBTK::VecCellData<double>;
template class Pointer<IBTK::VecCellData<double> >;

//////////////////////////////////////////////////////////////////////////////
