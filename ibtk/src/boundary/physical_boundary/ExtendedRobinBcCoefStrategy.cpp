// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include "ibtk/ExtendedRobinBcCoefStrategy.h"

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
ExtendedRobinBcCoefStrategy::setTargetPatchDataIndex(int target_data_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(target_data_idx >= 0);
#endif
    d_target_data_idx = target_data_idx;
    return;
} // setTargetPatchDataIndex

void
ExtendedRobinBcCoefStrategy::clearTargetPatchDataIndex()
{
    d_target_data_idx = invalid_index;
    return;
} // clearTargetPatchDataIndex

void
ExtendedRobinBcCoefStrategy::setHomogeneousBc(bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
} // setHomogeneousBc

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
