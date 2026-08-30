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

#ifndef included_IBAMR_private_StaggeredStokesIBTimeSteppingUtilities_inl_h
#define included_IBAMR_private_StaggeredStokesIBTimeSteppingUtilities_inl_h

#include <ibamr/IBImplicitStrategy.h>
#include <ibamr/ibamr_enums.h>

#include <tbox/Utilities.h>

#include <limits>
#include <string>

namespace IBAMR
{
enum class StaggeredStokesIBVelocityState
{
    NEW,
    MIDPOINT_AVERAGE
};

struct StaggeredStokesIBTimeStepParameters
{
    double half_time = std::numeric_limits<double>::quiet_NaN();
    double velocity_time = std::numeric_limits<double>::quiet_NaN();
    double force_time = std::numeric_limits<double>::quiet_NaN();
    StaggeredStokesIBVelocityState velocity_state = StaggeredStokesIBVelocityState::NEW;
    double nonlinear_force_scale = std::numeric_limits<double>::quiet_NaN();
    double jacobian_force_scale = std::numeric_limits<double>::quiet_NaN();
    // Location of force evaluation between current and updated positions.
    double force_position_fraction = std::numeric_limits<double>::quiet_NaN();
};

inline StaggeredStokesIBTimeStepParameters
get_staggered_stokes_ib_time_step_parameters(const TimeSteppingType time_stepping_type,
                                             const double current_time,
                                             const double new_time,
                                             const std::string& caller)
{
    StaggeredStokesIBTimeStepParameters parameters;
    parameters.half_time = current_time + 0.5 * (new_time - current_time);
    switch (time_stepping_type)
    {
    case BACKWARD_EULER:
        parameters.velocity_time = new_time;
        parameters.force_time = new_time;
        parameters.velocity_state = StaggeredStokesIBVelocityState::NEW;
        parameters.nonlinear_force_scale = 1.0;
        parameters.jacobian_force_scale = 1.0;
        parameters.force_position_fraction = 1.0;
        break;
    case TRAPEZOIDAL_RULE:
        parameters.velocity_time = new_time;
        parameters.force_time = new_time;
        parameters.velocity_state = StaggeredStokesIBVelocityState::NEW;
        parameters.nonlinear_force_scale = 0.5;
        parameters.jacobian_force_scale = 0.5;
        parameters.force_position_fraction = 1.0;
        break;
    case MIDPOINT_RULE:
        parameters.velocity_time = parameters.half_time;
        parameters.force_time = parameters.half_time;
        parameters.velocity_state = StaggeredStokesIBVelocityState::MIDPOINT_AVERAGE;
        parameters.nonlinear_force_scale = 1.0;
        parameters.jacobian_force_scale = 0.5;
        parameters.force_position_fraction = 0.5;
        break;
    default:
        TBOX_ERROR(caller << ": unsupported time stepping type\n");
    }
    return parameters;
}

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
        TBOX_ERROR(caller << ": unsupported time stepping type\n");
    }
    return;
}
} // namespace IBAMR

#endif // #ifndef included_IBAMR_private_StaggeredStokesIBTimeSteppingUtilities_inl_h
