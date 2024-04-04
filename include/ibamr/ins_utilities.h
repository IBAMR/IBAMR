// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2024 by the IBAMR developers
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

#ifndef included_IBAMR_ins_utilities
#define included_IBAMR_ins_utilities

#include <ibtk/config.h>

#include "tbox/Pointer.h"

namespace IBTK
{
class HierarchyMathOps;
}

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////
namespace IBAMR
{
/*!
 * \return Integral of inflow \f$ -\vec{u} \cdot \vec{n}\f$ at a physical boundary.
 *
 * \param location_idx of the boundary at which the integral is evaluated
 */
double computeNetInflowPhysicalBoundary(SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                        int u_idx,
                                        int bdry_loc_idx);
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_utilities
