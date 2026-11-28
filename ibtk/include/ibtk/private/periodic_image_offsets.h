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

#ifndef included_IBTK_periodic_image_offsets
#define included_IBTK_periodic_image_offsets

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <IntVector.h>

#include <vector>

namespace IBTK
{
namespace detail
{
/*!
 * \brief Return the offsets of the periodic images of a box, including the zero offset: the combinations of -1, 0, and
 * 1 times periodic_shift along the axes with a nonzero shift.
 */
std::vector<SAMRAI::hier::IntVector<NDIM>> periodic_image_offsets(const SAMRAI::hier::IntVector<NDIM>& periodic_shift);
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_periodic_image_offsets
