// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_CubicCoarsenLayer
#define included_IBTK_CubicCoarsenLayer

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ArrayData.h>
#include <Box.h>
#include <IntVector.h>

/////////////////////////////// FUNCTIONS ////////////////////////////////////

namespace IBTK
{
/*!
 * Set the values of \a data outside \a box to quiet NaN.
 */
void mark_cubic_coarsen_layer(SAMRAI::pdat::ArrayData<NDIM, double>& data, const SAMRAI::hier::Box<NDIM>& box);

/*!
 * Replace each NaN value in the layer of width(d) positions beyond \a box along each axis d (width(d) is 0 or 1) by
 * the cubic extrapolation from the values of \a data nearest to the position, along the smallest set of axes that
 * makes all of those values valid ones. A value of \a data that is not NaN is valid; so is every value in \a box.
 * Along an axis where \a box has fewer than four positions, the extrapolation has the highest order the positions
 * allow. A NaN inside \a box is not replaced.
 */
void extrapolate_cubic_coarsen_layer(SAMRAI::pdat::ArrayData<NDIM, double>& data,
                                     const SAMRAI::hier::Box<NDIM>& box,
                                     const SAMRAI::hier::IntVector<NDIM>& width);
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CubicCoarsenLayer
