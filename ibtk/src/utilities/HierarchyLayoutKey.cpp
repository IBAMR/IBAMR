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

#include <ibtk/private/HierarchyLayoutKey.h>

#include <tbox/Pointer.h>

#include <GridGeometry.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>

#include <utility>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
bool
equal_boxes(const BoxArray<NDIM>& a, const BoxArray<NDIM>& b)
{
    if (a.size() != b.size())
    {
        return false;
    }
    for (int i = 0; i < a.size(); ++i)
    {
        if (a[i] != b[i])
        {
            return false;
        }
    }
    return true;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
HierarchyLayoutKey::changed(const PatchHierarchy<NDIM>& hierarchy, const int coarsest_ln, const int finest_ln)
{
    HierarchyLayoutKey key(hierarchy, coarsest_ln, finest_ln);
    if (key == *this)
    {
        return false;
    }
    *this = std::move(key);
    return true;
} // changed

/////////////////////////////// PRIVATE //////////////////////////////////////

HierarchyLayoutKey::HierarchyLayoutKey(const PatchHierarchy<NDIM>& hierarchy,
                                       const int coarsest_ln,
                                       const int finest_ln)
    : d_coarsest_ln(coarsest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        const Pointer<PatchLevel<NDIM>> level = hierarchy.getPatchLevel(ln);
        Level key_level;
        key_level.boxes = level->getBoxes();
        const ProcessorMapping& mapping = level->getProcessorMapping();
        key_level.mapping.resize(mapping.getSizeOfMappingArray());
        for (int p = 0; p < mapping.getSizeOfMappingArray(); ++p)
        {
            key_level.mapping[p] = mapping.getProcessorAssignment(p);
        }
        key_level.ratio = level->getRatio();
        key_level.periodic_shift = level->getGridGeometry()->getPeriodicShift(key_level.ratio);
        d_levels.push_back(std::move(key_level));
    }
} // HierarchyLayoutKey

bool
HierarchyLayoutKey::operator==(const HierarchyLayoutKey& other) const
{
    return d_coarsest_ln == other.d_coarsest_ln && d_levels == other.d_levels;
} // operator==

bool
HierarchyLayoutKey::Level::operator==(const Level& other) const
{
    return equal_boxes(boxes, other.boxes) && mapping == other.mapping && ratio == other.ratio &&
           periodic_shift == other.periodic_shift;
} // operator==

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
