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

#ifndef included_IBAMR_BoundaryTreatment
#define included_IBAMR_BoundaryTreatment

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <string>

/////////////////////////////// STRUCT DEFINITION ////////////////////////////

namespace IBAMR
{
/*!
 * \brief Struct BoundaryTreatment specifies the discretization of the velocity and pressure at physical boundaries for
 * a fluid integrator and the objects it creates or serves.
 */
struct BoundaryTreatment
{
    /*!
     * The interpolation type of the ghost values at physical boundaries at which boundary conditions are present,
     * "LINEAR" or "QUADRATIC".
     */
    std::string bdry_interp_type = "LINEAR";

    /*!
     * The extrapolation type of coarse data beyond physical boundaries for coarse-fine interpolation.
     */
    std::string bdry_extrap_type = "LINEAR";
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_BoundaryTreatment
