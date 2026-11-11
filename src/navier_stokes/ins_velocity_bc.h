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

#include <BoundaryBox.h>
#include <Patch.h>
#include <RobinBcCoefStrategy.h>

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
 * Return VELOCITY if the Robin coefficients (a, b) are exactly (1, 0) and TRACTION if they are exactly (0, 1). Stop
 * with an error that names \a caller for any other coefficients.
 */
Type get_type(double a, double b, const char* caller);

/*!
 * Return whether the TRACTION conditions of the velocity component with the Robin coefficient object \a comp_bc_coef
 * make its boundary data on \a bdry_box of \a patch at \a fill_time depend on the velocity: the component is tangential
 * to the box and a face that the fill sets has TRACTION coefficients.
 */
bool traction_data_depend_on_velocity(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* comp_bc_coef,
                                      unsigned int comp_idx,
                                      const SAMRAI::hier::Patch<NDIM>& patch,
                                      const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                                      double fill_time);
} // namespace velocity_bc
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_velocity_bc
