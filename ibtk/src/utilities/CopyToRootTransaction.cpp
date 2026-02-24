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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/CopyToRootTransaction.h"
#include "ibtk/samrai_compatibility_names.h"

#include "BoxArray.h"
#include "SAMRAIBox.h"
#include "SAMRAIBoxGeometry.h"
#include "SAMRAIBoxOverlap.h"
#include "SAMRAIGridGeometry.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchData.h"
#include "SAMRAIPatchDataFactory.h"
#include "SAMRAIPatchDescriptor.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"

#include <ostream>

#include "ibtk/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CopyToRootTransaction::CopyToRootTransaction(const int src_proc,
                                             const int dst_proc,
                                             SAMRAIPointer<SAMRAIPatchLevel> patch_level,
                                             const int src_patch_data_idx,
                                             SAMRAIPointer<SAMRAIPatchData> dst_patch_data)
    : d_src_proc(src_proc),
      d_dst_proc(dst_proc),
      d_patch_level(patch_level),
      d_src_patch_data_idx(src_patch_data_idx),
      d_dst_patch_data(dst_patch_data)
{
    // intentionally blank
    return;
} // CopyToRootTransaction

SAMRAIPointer<SAMRAIPatchData>
CopyToRootTransaction::getRootPatchData() const
{
    return d_dst_patch_data;
} // getRootPatchData

bool
CopyToRootTransaction::canEstimateIncomingMessageSize()
{
    return false;
} // canEstimateIncomingMessageSize

int
CopyToRootTransaction::computeIncomingMessageSize()
{
    return 0;
} // computeIncomingMessageSize

int
CopyToRootTransaction::computeOutgoingMessageSize()
{
    SAMRAIPointer<SAMRAIPatchDataFactory> pdat_factory =
        d_patch_level->getPatchDescriptor()->getPatchDataFactory(d_src_patch_data_idx);

    SAMRAIPointer<SAMRAIGridGeometry> grid_geom = d_patch_level->getGridGeometry();
#if !defined(NDEBUG)
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
#endif
    const SAMRAIBox& dst_box = grid_geom->getPhysicalDomain()[0];
    SAMRAIPointer<SAMRAIBoxGeometry> dst_box_geometry = pdat_factory->getBoxGeometry(dst_box);

    int size = AbstractStream::sizeofInt();
    for (SAMRAIPatchLevel::Iterator p(d_patch_level); p; p++)
    {
        const int src_patch_num = p();
        size += AbstractStream::sizeofInt();
        SAMRAIPointer<SAMRAIPatch> patch = d_patch_level->getPatch(src_patch_num);
        const SAMRAIBox& src_box = patch->getBox();
        SAMRAIPointer<SAMRAIBoxGeometry> src_box_geometry = pdat_factory->getBoxGeometry(src_box);
        const SAMRAIBox& src_mask = dst_box;
        const bool overwrite_interior = true;
        const SAMRAIIntVector src_shift = 0;
        SAMRAIPointer<SAMRAIBoxOverlap> box_overlap =
            dst_box_geometry->calculateOverlap(*src_box_geometry, src_mask, overwrite_interior, src_shift);
        size += patch->getPatchData(d_src_patch_data_idx)->getDataStreamSize(*box_overlap);
    }
    return size;
} // computeOutgoingMessageSize

int
CopyToRootTransaction::getSourceProcessor()
{
    return d_src_proc;
} // getSourceProcessor

int
CopyToRootTransaction::getDestinationProcessor()
{
    return d_dst_proc;
} // getDestinationProcessor

void
CopyToRootTransaction::packStream(AbstractStream& stream)
{
    SAMRAIPointer<SAMRAIPatchDataFactory> pdat_factory =
        d_patch_level->getPatchDescriptor()->getPatchDataFactory(d_src_patch_data_idx);

    SAMRAIPointer<SAMRAIGridGeometry> grid_geom = d_patch_level->getGridGeometry();
#if !defined(NDEBUG)
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
#endif
    const SAMRAIBox& dst_box = grid_geom->getPhysicalDomain()[0];
    SAMRAIPointer<SAMRAIBoxGeometry> dst_box_geometry = pdat_factory->getBoxGeometry(dst_box);

    int src_patch_count = 0;
    for (SAMRAIPatchLevel::Iterator p(d_patch_level); p; p++)
    {
        ++src_patch_count;
    }
    stream << src_patch_count;

    for (SAMRAIPatchLevel::Iterator p(d_patch_level); p; p++)
    {
        const int src_patch_num = p();
        stream << src_patch_num;
        SAMRAIPointer<SAMRAIPatch> patch = d_patch_level->getPatch(src_patch_num);
        const SAMRAIBox& src_box = patch->getBox();
        SAMRAIPointer<SAMRAIBoxGeometry> src_box_geometry = pdat_factory->getBoxGeometry(src_box);
        const SAMRAIBox& src_mask = dst_box;
        const bool overwrite_interior = true;
        const SAMRAIIntVector src_shift = 0;
        SAMRAIPointer<SAMRAIBoxOverlap> box_overlap =
            dst_box_geometry->calculateOverlap(*src_box_geometry, src_mask, overwrite_interior, src_shift);
        patch->getPatchData(d_src_patch_data_idx)->packStream(stream, *box_overlap);
    }
    return;
} // packStream

void
CopyToRootTransaction::unpackStream(AbstractStream& stream)
{
    SAMRAIPointer<SAMRAIPatchDataFactory> pdat_factory =
        d_patch_level->getPatchDescriptor()->getPatchDataFactory(d_src_patch_data_idx);

    SAMRAIPointer<SAMRAIGridGeometry> grid_geom = d_patch_level->getGridGeometry();
#if !defined(NDEBUG)
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
#endif
    const SAMRAIBox& dst_box = grid_geom->getPhysicalDomain()[0];
    SAMRAIPointer<SAMRAIBoxGeometry> dst_box_geometry = pdat_factory->getBoxGeometry(dst_box);

    int src_patch_count;
    stream >> src_patch_count;
    for (int p = 0; p < src_patch_count; ++p)
    {
        int src_patch_num;
        stream >> src_patch_num;
        const SAMRAIBox& src_box = d_patch_level->getBoxes()[src_patch_num];
        SAMRAIPointer<SAMRAIBoxGeometry> src_box_geometry = pdat_factory->getBoxGeometry(src_box);
        const SAMRAIBox& src_mask = dst_box;
        const bool overwrite_interior = true;
        const SAMRAIIntVector src_shift = 0;
        SAMRAIPointer<SAMRAIBoxOverlap> box_overlap =
            dst_box_geometry->calculateOverlap(*src_box_geometry, src_mask, overwrite_interior, src_shift);
        d_dst_patch_data->unpackStream(stream, *box_overlap);
    }
    return;
} // unpackStream

void
CopyToRootTransaction::copyLocalData()
{
    SAMRAIPointer<SAMRAIPatchDataFactory> pdat_factory =
        d_patch_level->getPatchDescriptor()->getPatchDataFactory(d_src_patch_data_idx);

    SAMRAIPointer<SAMRAIGridGeometry> grid_geom = d_patch_level->getGridGeometry();
#if !defined(NDEBUG)
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
#endif
    const SAMRAIBox& dst_box = grid_geom->getPhysicalDomain()[0];
    SAMRAIPointer<SAMRAIBoxGeometry> dst_box_geometry = pdat_factory->getBoxGeometry(dst_box);

    for (SAMRAIPatchLevel::Iterator p(d_patch_level); p; p++)
    {
        int src_patch_num = p();
        SAMRAIPointer<SAMRAIPatch> patch = d_patch_level->getPatch(src_patch_num);
        const SAMRAIBox& src_box = patch->getBox();
        SAMRAIPointer<SAMRAIBoxGeometry> src_box_geometry = pdat_factory->getBoxGeometry(src_box);
        const SAMRAIBox& src_mask = dst_box;
        const bool overwrite_interior = true;
        const SAMRAIIntVector src_shift = 0;
        SAMRAIPointer<SAMRAIBoxOverlap> box_overlap =
            dst_box_geometry->calculateOverlap(*src_box_geometry, src_mask, overwrite_interior, src_shift);
        d_dst_patch_data->copy(*patch->getPatchData(d_src_patch_data_idx), *box_overlap);
    }
    return;
} // copyLocalData

void
CopyToRootTransaction::printClassData(std::ostream& stream) const
{
    stream << "CopyToRootTransaction::printClassData() is not implemented\n";
    return;
} // printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
