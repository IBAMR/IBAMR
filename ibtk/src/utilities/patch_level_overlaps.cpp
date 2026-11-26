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

#include <ibtk/private/patch_level_overlaps.h>

#include <tbox/Array.h>

#include <BoxArray.h>
#include <BoxTree.h>

#include <ibtk/app_namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
std::optional<std::pair<int, int>>
find_overlapping_patches(Pointer<PatchLevel<NDIM>> level)
{
    const BoxArray<NDIM>& boxes = level->getBoxes();
    Pointer<BoxTree<NDIM>> box_tree = level->getBoxTree();
    tbox::Array<int> indices;
    for (int i = 0; i < boxes.getNumberOfBoxes(); ++i)
    {
        // The indices include i itself.
        box_tree->findOverlapIndices(indices, boxes[i]);
        int j = boxes.getNumberOfBoxes();
        for (int k = 0; k < indices.getSize(); ++k)
        {
            if (indices[k] > i && indices[k] < j) j = indices[k];
        }
        if (j < boxes.getNumberOfBoxes()) return std::make_pair(i, j);
    }
    return std::nullopt;
}
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
