// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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

#include <Box.h>

#include <vector>

namespace IBTK
{
/**
 * Given a set of boxes (which describe a region in index space), return
 * another set of boxes whose union covers the same index space. The
 * returned set of boxes is formed by merging boxes in @p boxes along
 * their longest edges.
 */
std::vector<SAMRAI::hier::Box<NDIM> > merge_boxes_by_longest_edge(const std::vector<SAMRAI::hier::Box<NDIM> >& boxes);
} // namespace IBTK

#endif
