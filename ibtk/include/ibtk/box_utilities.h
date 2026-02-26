// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2024 by the IBAMR developers
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

#ifndef included_IBTK_box_utilities
#define included_IBTK_box_utilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIBox.h>

#include <vector>

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBTK
{

/**
 * Given a set of boxes (which describe a region in index space), return
 * another set of boxes whose union covers the same index space. The
 * returned set of boxes is formed by merging boxes in @p boxes along
 * their longest edges.
 */
std::vector<SAMRAIBox> merge_boxes_by_longest_edge(const std::vector<SAMRAIBox>& boxes);

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_ibtk_utilities
