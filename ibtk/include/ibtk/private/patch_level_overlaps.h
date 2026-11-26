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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_private_patch_level_overlaps
#define included_IBTK_private_patch_level_overlaps

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <PatchLevel.h>

#include <optional>
#include <utility>

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBTK
{
namespace detail
{
/**
 * Return the lowest-numbered pair of overlapping patches of @p level, or nothing. Every processor returns the same pair
 * without communication.
 */
std::optional<std::pair<int, int>>
find_overlapping_patches(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level);
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_private_patch_level_overlaps
