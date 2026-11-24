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

#include <ibtk/private/HierarchyLevelsIdentity.h>

#include <tbox/Pointer.h>

#include <PatchLevel.h>

#include <utility>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
HierarchyLevelsIdentity::changed(const PatchHierarchy<NDIM>& hierarchy, const int coarsest_ln, const int finest_ln)
{
    std::vector<Pointer<PatchLevel<NDIM>>> levels;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        levels.push_back(hierarchy.getPatchLevel(ln));
    }
    bool same = d_hierarchy == &hierarchy && d_coarsest_ln == coarsest_ln && d_levels.size() == levels.size();
    for (std::size_t k = 0; same && k < levels.size(); ++k)
    {
        same = d_levels[k].getPointer() == levels[k].getPointer();
    }
    if (same)
    {
        return false;
    }
    d_hierarchy = &hierarchy;
    d_coarsest_ln = coarsest_ln;
    d_levels = std::move(levels);
    return true;
} // changed

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
