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

#include <ibtk/private/data_kind.h>

#include <tbox/Pointer.h>

#include <CellDataFactory.h>
#include <EdgeDataFactory.h>
#include <FaceDataFactory.h>
#include <NodeDataFactory.h>
#include <PatchDataFactory.h>
#include <PatchDescriptor.h>
#include <SideDataFactory.h>
#include <VariableDatabase.h>

#include <array>
#include <string>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
DataKind::operator==(const DataKind& other) const
{
    return centering == other.centering && depth == other.depth;
} // operator==

bool
DataKind::operator!=(const DataKind& other) const
{
    return !(*this == other);
} // operator!=

DataKind
get_data_kind(const int data_idx)
{
    Pointer<PatchDataFactory<NDIM>> factory =
        VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor()->getPatchDataFactory(data_idx);
    if (Pointer<CellDataFactory<NDIM, double>> f = factory)
    {
        return { DataCentering::CELL, f->getDefaultDepth() };
    }
    if (Pointer<NodeDataFactory<NDIM, double>> f = factory)
    {
        return { DataCentering::NODE, f->getDefaultDepth() };
    }
    if (Pointer<SideDataFactory<NDIM, double>> f = factory)
    {
        return { DataCentering::SIDE, f->getDefaultDepth() };
    }
    if (Pointer<EdgeDataFactory<NDIM, double>> f = factory)
    {
        return { DataCentering::EDGE, f->getDefaultDepth() };
    }
    if (Pointer<FaceDataFactory<NDIM, double>> f = factory)
    {
        return { DataCentering::FACE, f->getDefaultDepth() };
    }
    return { DataCentering::OTHER, 0 };
} // get_data_kind

std::string
data_kind_to_string(const DataKind& kind)
{
    static const std::array<std::string, 6> names = { "cell", "node", "side", "edge", "face", "other" };
    return names[static_cast<std::size_t>(kind.centering)] + "-centered data of depth " + std::to_string(kind.depth);
} // data_kind_to_string

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
