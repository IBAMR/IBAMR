// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/private/ArrayDataTransaction.h>

#include <tbox/AbstractStream.h>
#include <tbox/Utilities.h>

#include <BoxGeometry.h>
#include <BoxGeometryFillPattern.h>
#include <BoxOverlap.h>
#include <CellOverlap.h>
#include <PatchDataFactory.h>
#include <PatchDescriptor.h>
#include <SideOverlap.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <ostream>
#include <typeinfo>
#include <utility>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

ArrayDataTransaction::ArrayDataTransaction(const int src_proc,
                                           const int dst_proc,
                                           const Box<NDIM>& dst_box,
                                           const IntVector<NDIM>& src_shift,
                                           const int depth,
                                           ArrayAccessor src_array,
                                           ArrayAccessor dst_array)
    : d_src_proc(src_proc),
      d_dst_proc(dst_proc),
      d_dst_box(dst_box),
      d_src_shift(src_shift),
      d_depth(depth),
      d_src_array(std::move(src_array)),
      d_dst_array(std::move(dst_array))
{
    // intentionally blank
} // ArrayDataTransaction

bool
ArrayDataTransaction::canEstimateIncomingMessageSize()
{
    return true;
} // canEstimateIncomingMessageSize

int
ArrayDataTransaction::computeIncomingMessageSize()
{
    return AbstractStream::sizeofDouble(d_dst_box.size() * d_depth);
} // computeIncomingMessageSize

int
ArrayDataTransaction::computeOutgoingMessageSize()
{
    return AbstractStream::sizeofDouble(d_dst_box.size() * d_depth);
} // computeOutgoingMessageSize

int
ArrayDataTransaction::getSourceProcessor()
{
    return d_src_proc;
} // getSourceProcessor

int
ArrayDataTransaction::getDestinationProcessor()
{
    return d_dst_proc;
} // getDestinationProcessor

void
ArrayDataTransaction::packStream(AbstractStream& stream)
{
    d_src_array().packStream(stream, d_dst_box, d_src_shift);
} // packStream

void
ArrayDataTransaction::unpackStream(AbstractStream& stream)
{
    d_dst_array().unpackStream(stream, d_dst_box, d_src_shift);
} // unpackStream

void
ArrayDataTransaction::copyLocalData()
{
    d_dst_array().copy(d_src_array(), d_dst_box, d_src_shift);
} // copyLocalData

void
ArrayDataTransaction::printClassData(std::ostream& stream) const
{
    stream << "ArrayDataTransaction: " << d_src_proc << " -> " << d_dst_proc << ", box " << d_dst_box << ", shift "
           << d_src_shift << ", depth " << d_depth << "\n";
} // printClassData

bool
FillPatternShape::operator==(const FillPatternShape& that) const
{
    return kind == that.kind && width == that.width && overwrite_interior == that.overwrite_interior &&
           other == that.other;
} // operator==

bool
FillPatternShape::operator!=(const FillPatternShape& that) const
{
    return !(*this == that);
} // operator!=

FillPatternShape
fill_pattern_shape(Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern, const IntVector<NDIM>& ghost_width)
{
    FillPatternShape shape;
    shape.width = ghost_width;
    if (!fill_pattern)
    {
        return shape;
    }
    const SAMRAI::xfer::VariableFillPattern<NDIM>& pattern_object = *fill_pattern;
    if (typeid(pattern_object) == typeid(SAMRAI::xfer::BoxGeometryFillPattern<NDIM>))
    {
        return shape;
    }
    IntVector<NDIM> stencil_width;
    if (typeid(pattern_object) == typeid(CellNoCornersFillPattern))
    {
        Pointer<CellNoCornersFillPattern> pattern = fill_pattern;
        shape.kind = FillPatternShape::Kind::CELL_NO_CORNERS;
        stencil_width = pattern->getStencilWidth();
        shape.overwrite_interior = pattern->getOverwriteInterior();
    }
    else if (typeid(pattern_object) == typeid(SideNoCornersFillPattern))
    {
        Pointer<SideNoCornersFillPattern> pattern = fill_pattern;
        shape.kind = FillPatternShape::Kind::SIDE_NO_CORNERS;
        stencil_width = pattern->getStencilWidth();
        shape.overwrite_interior = pattern->getOverwriteInterior();
    }
    else
    {
        shape.kind = FillPatternShape::Kind::OTHER;
        shape.other = fill_pattern.getPointer();
        return shape;
    }
    for (int d = 0; d < NDIM; ++d)
    {
        shape.width(d) = std::min(ghost_width(d), stencil_width(d));
    }
    return shape;
} // fill_pattern_shape

std::vector<BoxList<NDIM>>
fill_pattern_indices(const int data_idx,
                     Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern,
                     const Box<NDIM>& patch_box,
                     const int ln)
{
    Pointer<PatchDataFactory<NDIM>> factory =
        VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor()->getPatchDataFactory(data_idx);
    Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> pattern = fill_pattern;
    if (!pattern)
    {
        pattern = new SAMRAI::xfer::BoxGeometryFillPattern<NDIM>();
    }
    // A source that covers the ghost box, as a refine schedule to level ln computes the overlap.
    const Box<NDIM> ghost_box = Box<NDIM>::grow(patch_box, factory->getGhostCellWidth());
    pattern->setTargetPatchLevelNumber(ln);
    Pointer<BoxOverlap<NDIM>> overlap = pattern->calculateOverlapOnLevel(*factory->getBoxGeometry(patch_box),
                                                                         *factory->getBoxGeometry(ghost_box),
                                                                         patch_box,
                                                                         ghost_box,
                                                                         /*overwrite_interior*/ true,
                                                                         IntVector<NDIM>(0),
                                                                         ln,
                                                                         ln);
    pattern->setTargetPatchLevelNumber(-1);
    if (Pointer<CellOverlap<NDIM>> cell_overlap = overlap)
    {
        return { cell_overlap->getDestinationBoxList() };
    }
    Pointer<SideOverlap<NDIM>> side_overlap = overlap;
    TBOX_ASSERT(side_overlap);
    std::vector<BoxList<NDIM>> indices(NDIM);
    for (int axis = 0; axis < NDIM; ++axis)
    {
        indices[axis] = side_overlap->getDestinationBoxList(axis);
    }
    return indices;
} // fill_pattern_indices

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
