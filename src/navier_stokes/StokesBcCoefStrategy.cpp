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

#include "ibamr/StokesBcCoefStrategy.h"
#include "ibamr/ibamr_enums.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

namespace IBAMR
{
class StokesSpecifications;
} // namespace IBAMR

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
StokesBcCoefStrategy::setStokesSpecifications(const StokesSpecifications* const problem_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(problem_coefs);
#endif
    d_problem_coefs = problem_coefs;
    return;
} // setStokesSpecifications

void
StokesBcCoefStrategy::setTargetVelocityPatchDataIndex(int u_target_data_idx)
{
    d_u_target_data_idx = u_target_data_idx;
    return;
} // setTargetVelocityPatchDataIndex

void
StokesBcCoefStrategy::clearTargetVelocityPatchDataIndex()
{
    d_u_target_data_idx = invalid_index;
    return;
} // clearTargetVelocityPatchDataIndex

void
StokesBcCoefStrategy::setTargetPressurePatchDataIndex(int p_target_data_idx)
{
    d_p_target_data_idx = p_target_data_idx;
    return;
} // setPressurePatchDataIndex

void
StokesBcCoefStrategy::clearTargetPressurePatchDataIndex()
{
    d_p_target_data_idx = invalid_index;
    return;
} // clearPressurePatchDataIndex

void
StokesBcCoefStrategy::setTractionBcType(TractionBcType bc_type)
{
    d_traction_bc_type = bc_type;
    return;
} // setTractionBcType

TractionBcType
StokesBcCoefStrategy::getTractionBcType() const
{
    return d_traction_bc_type;
} // getTractionBcType

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
