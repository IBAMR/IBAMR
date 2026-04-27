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

#ifndef included_IBAMR_private_StaggeredStokesIBOperatorUtilities_inl_h
#define included_IBAMR_private_StaggeredStokesIBOperatorUtilities_inl_h

#include <ibamr/StaggeredStokesIBOperator.h>
#include <ibamr/private/StaggeredStokesIBTimeSteppingUtilities-inl.h>

#include <ibtk/RobinPhysBdryPatchStrategy.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <HierarchyDataOpsReal.h>
#include <SAMRAIVectorReal.h>

#include <string>
#include <utility>

namespace IBAMR
{
inline void
set_staggered_stokes_ib_operator_times(const SAMRAI::tbox::Pointer<StaggeredStokesOperator>& stokes_op,
                                       const std::pair<double, double>& time_interval,
                                       const double solution_time)
{
    if (!stokes_op) return;
    stokes_op->setTimeInterval(time_interval.first, time_interval.second);
    stokes_op->setSolutionTime(solution_time);
    return;
}

inline void
apply_staggered_stokes_ib_stokes_operator(const StaggeredStokesIBOperator::Context& ctx,
                                          const std::pair<double, double>& time_interval,
                                          const double solution_time,
                                          SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                          SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y)
{
    set_staggered_stokes_ib_operator_times(ctx.stokes_op, time_interval, solution_time);
    ctx.stokes_op->setHomogeneousBc(true);
    ctx.stokes_op->apply(x, y);
    return;
}

inline void
copy_staggered_stokes_ib_velocity_state(const StaggeredStokesIBOperator::Context& ctx,
                                        const StaggeredStokesIBTimeSteppingSchedule& schedule,
                                        const int u_new_idx,
                                        const int u_current_idx,
                                        const std::string& caller)
{
    switch (schedule.velocity_state)
    {
    case StaggeredStokesIBVelocityState::NEW:
        ctx.hier_velocity_data_ops->copyData(ctx.u_idx, u_new_idx);
        break;
    case StaggeredStokesIBVelocityState::MIDPOINT_AVERAGE:
        ctx.hier_velocity_data_ops->linearSum(ctx.u_idx, 0.5, u_new_idx, 0.5, u_current_idx);
        break;
    default:
        TBOX_ERROR(caller << ": unsupported velocity state\n");
    }
    return;
}

inline void
set_staggered_stokes_ib_velocity_bdry_state(const StaggeredStokesIBOperator::Context& ctx,
                                            const int patch_data_idx,
                                            const bool homogeneous_bc)
{
    if (!ctx.u_phys_bdry_op) return;
    ctx.u_phys_bdry_op->setPatchDataIndex(patch_data_idx);
    ctx.u_phys_bdry_op->setHomogeneousBc(homogeneous_bc);
    return;
}

inline void
prepare_staggered_stokes_ib_bc_forwarding(const StaggeredStokesIBOperator::Context& ctx,
                                          const std::pair<double, double>& time_interval,
                                          const double solution_time,
                                          const bool homogeneous_bc,
                                          const std::string& caller)
{
    if (!ctx.stokes_op)
    {
        TBOX_ERROR(caller << ": missing Stokes operator\n");
    }
    set_staggered_stokes_ib_operator_times(ctx.stokes_op, time_interval, solution_time);
    ctx.stokes_op->setHomogeneousBc(homogeneous_bc);
    return;
}
} // namespace IBAMR

#endif // #ifndef included_IBAMR_private_StaggeredStokesIBOperatorUtilities_inl_h
