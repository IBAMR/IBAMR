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

#ifndef included_IBAMR_ins_velocity_bc
#define included_IBAMR_ins_velocity_bc

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

// Classification of the velocity boundary conditions of the incompressible Navier-Stokes classes. This header is not
// installed.
namespace IBAMR
{
namespace velocity_bc
{
enum class Type
{
    VELOCITY,
    TRACTION
};

/*!
 * Return VELOCITY if the Robin coefficients (a, b) are exactly (1, 0) and TRACTION if they are exactly (0, 1). Any
 * other coefficients stop the program with an error that names \a function.
 */
Type get_type(double a, double b, const char* function);
} // namespace velocity_bc
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_velocity_bc
