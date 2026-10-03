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

#ifndef included_IBAMR_private_StaggeredStokesIBTimeSteppingUtilities_h
#define included_IBAMR_private_StaggeredStokesIBTimeSteppingUtilities_h

#include <ibamr/config.h>

#include <ibamr/IBImplicitStrategy.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/RobinPhysBdryPatchStrategy.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <IntVector.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>

#include <limits>
#include <string>

namespace IBAMR
{
// Time-stepping parameters and helpers shared by the nonlinear and Jacobian operators.
namespace detail
{
struct StaggeredStokesIBTimeStepParameters
{
    double evaluation_time = std::numeric_limits<double>::quiet_NaN();
    double nonlinear_force_scale = std::numeric_limits<double>::quiet_NaN();
    double jacobian_force_scale = std::numeric_limits<double>::quiet_NaN();
    // Location of force evaluation between current and updated positions.
    double force_position_fraction = std::numeric_limits<double>::quiet_NaN();
};

/*! \brief Require a time-stepping rule supported by the Stokes-IB operators. */
inline void
require_staggered_stokes_ib_time_stepping_type(const TimeSteppingType time_stepping_type, const std::string& caller)
{
    switch (time_stepping_type)
    {
    case BACKWARD_EULER:
    case TRAPEZOIDAL_RULE:
    case MIDPOINT_RULE:
        return;
    default:
        TBOX_ERROR(
            caller << ":\n  unsupported time stepping type; use BACKWARD_EULER, TRAPEZOIDAL_RULE, or MIDPOINT_RULE.\n");
    }
}

/*! \brief Select the coupling time and force scaling for the time-stepping rule. */
inline StaggeredStokesIBTimeStepParameters
get_staggered_stokes_ib_time_step_parameters(const TimeSteppingType time_stepping_type,
                                             const double current_time,
                                             const double new_time,
                                             const std::string& caller)
{
    StaggeredStokesIBTimeStepParameters parameters;
    switch (time_stepping_type)
    {
    case BACKWARD_EULER:
        parameters.evaluation_time = new_time;
        parameters.nonlinear_force_scale = 1.0;
        parameters.jacobian_force_scale = 1.0;
        parameters.force_position_fraction = 1.0;
        break;
    case TRAPEZOIDAL_RULE:
        parameters.evaluation_time = new_time;
        parameters.nonlinear_force_scale = 0.5;
        parameters.jacobian_force_scale = 0.5;
        parameters.force_position_fraction = 1.0;
        break;
    case MIDPOINT_RULE:
        parameters.evaluation_time = current_time + 0.5 * (new_time - current_time);
        parameters.nonlinear_force_scale = 1.0;
        parameters.jacobian_force_scale = 0.5;
        parameters.force_position_fraction = 0.5;
        break;
    default:
        TBOX_ERROR(caller << ":\n  unsupported time stepping type.\n");
    }
    return parameters;
}

/*! \brief Advance the strategy's positions with the selected time-stepping rule. */
inline void
advance_staggered_stokes_ib_strategy(IBImplicitStrategy& ib_implicit_ops,
                                     const TimeSteppingType time_stepping_type,
                                     const double current_time,
                                     const double new_time,
                                     const std::string& caller)
{
    switch (time_stepping_type)
    {
    case BACKWARD_EULER:
        ib_implicit_ops.backwardEulerStep(current_time, new_time);
        break;
    case TRAPEZOIDAL_RULE:
        ib_implicit_ops.trapezoidalStep(current_time, new_time);
        break;
    case MIDPOINT_RULE:
        ib_implicit_ops.midpointStep(current_time, new_time);
        break;
    default:
        TBOX_ERROR(caller << ":\n  unsupported time stepping type.\n");
    }
    return;
}

/*!
 * \brief Zero the spread force at f_idx where the velocity on the physical
 * boundary is prescribed, on levels coarsest_ln through finest_ln.
 *
 * The Stokes rows of those velocities impose the boundary condition, so the
 * IB force does not enter them. Homogeneous boundary filling with
 * f_phys_bdry_op, the strategy used for velocity interpolation and force
 * spreading, sets exactly those values to zero.
 */
inline void
zero_force_at_prescribed_boundary_velocity(IBTK::RobinPhysBdryPatchStrategy& f_phys_bdry_op,
                                           const int f_idx,
                                           SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                           const int coarsest_ln,
                                           const int finest_ln,
                                           const double data_time)
{
    f_phys_bdry_op.setPatchDataIndex(f_idx);
    f_phys_bdry_op.setHomogeneousBc(true);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            f_phys_bdry_op.setPhysicalBoundaryConditions(
                *level->getPatch(p()), data_time, SAMRAI::hier::IntVector<NDIM>(1));
        }
    }
    return;
}
} // namespace detail
} // namespace IBAMR

#endif // #ifndef included_IBAMR_private_StaggeredStokesIBTimeSteppingUtilities_h
