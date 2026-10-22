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

#include <tbox/Utilities.h>

#include "./ins_velocity_bc.h"

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
namespace velocity_bc
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

Type
get_type(const double a, const double b, const char* const function)
{
    const bool is_velocity = (a == 1.0 && b == 0.0);
    const bool is_traction = (a == 0.0 && b == 1.0);
    if (!is_velocity && !is_traction)
    {
        TBOX_ERROR(function << ":\n"
                            << "  unsupported velocity boundary condition coefficients (a, b) = (" << a << ", " << b
                            << ")\n"
                            << "  use (a, b) = (1, 0) to prescribe the velocity or (0, 1) to prescribe the traction\n");
    }
    return is_velocity ? Type::VELOCITY : Type::TRACTION;
} // get_type

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace velocity_bc
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
