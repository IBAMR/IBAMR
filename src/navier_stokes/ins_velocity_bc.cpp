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

#include <ibamr/private/ins_velocity_bc.h>

#include <tbox/Utilities.h>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
namespace detail
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityBcType
classify_velocity_bc(const double a, const double b)
{
    const bool is_velocity = (a == 1.0 && b == 0.0);
    const bool is_traction = (a == 0.0 && b == 1.0);
    if (!is_velocity && !is_traction)
    {
        TBOX_ERROR("classify_velocity_bc():\n"
                   << "  unsupported velocity boundary condition coefficients (a, b) = (" << a << ", " << b << ")\n"
                   << "  use (a, b) = (1, 0) to prescribe the velocity or (0, 1) to prescribe the traction\n");
    }
    return is_velocity ? VelocityBcType::VELOCITY : VelocityBcType::TRACTION;
} // classify_velocity_bc

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
