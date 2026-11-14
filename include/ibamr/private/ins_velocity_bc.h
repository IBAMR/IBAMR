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

namespace IBAMR
{
namespace detail
{
/*!
 * \brief The condition that a pair of Robin coefficients imposes on a velocity component.
 */
enum class VelocityBcType
{
    /*! The velocity is prescribed. */
    VELOCITY,
    /*! The traction is prescribed. */
    TRACTION
};

/*!
 * \brief Return VELOCITY if the Robin coefficients (a, b) are exactly (1, 0) and TRACTION if they are exactly (0, 1).
 * Stop with an error for any other coefficients.
 */
VelocityBcType classify_velocity_bc(double a, double b);
} // namespace detail
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_velocity_bc
