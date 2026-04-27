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

#include <ibamr/IBImplicitStrategy.h>
#include <ibamr/StaggeredStokesIBOperator.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/private/StaggeredStokesIBOperatorUtilities-inl.h>
#include <ibamr/private/StaggeredStokesIBTimeSteppingUtilities-inl.h>

#include <ibtk/RobinPhysBdryPatchStrategy.h>

#include <tbox/Pointer.h>

#include <PatchHierarchy.h>
#include <SAMRAIVectorReal.h>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesIBOperator::StaggeredStokesIBOperator(const std::string& object_name, const bool homogeneous_bc)
    : GeneralOperator(object_name, homogeneous_bc)
{
    // intentionally blank
    return;
} // StaggeredStokesIBOperator

StaggeredStokesIBOperator::~StaggeredStokesIBOperator()
{
    deallocateOperatorState();
    return;
} // ~StaggeredStokesIBOperator

void
StaggeredStokesIBOperator::setOperatorContext(const StaggeredStokesIBOperator::Context& ctx)
{
    d_ctx = ctx;
    return;
} // setOperatorContext

void
StaggeredStokesIBOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(getIsInitialized());
#endif
    if (!d_ctx.ib_implicit_ops || !d_ctx.stokes_op || !d_ctx.hier_velocity_data_ops)
    {
        TBOX_ERROR(d_object_name << "::apply(): incomplete operator context\n");
    }
    if (d_ctx.u_idx == IBTK::invalid_index || d_ctx.f_idx == IBTK::invalid_index)
    {
        TBOX_ERROR(d_object_name << "::apply(): invalid scratch data indices\n");
    }
    if (d_ctx.u_current_idx == IBTK::invalid_index)
    {
        TBOX_ERROR(d_object_name << "::apply(): invalid current velocity data index\n");
    }

    const double current_time = getTimeInterval().first;
    const double new_time = getTimeInterval().second;
    const auto schedule = get_staggered_stokes_ib_time_stepping_schedule(
        d_ctx.time_stepping_type, current_time, new_time, d_object_name + "::apply()");

    const int u_new_idx = x.getComponentDescriptorIndex(0);
    const int f_u_idx = y.getComponentDescriptorIndex(0);

    apply_staggered_stokes_ib_stokes_operator(d_ctx, getTimeInterval(), getSolutionTime(), x, y);
    copy_staggered_stokes_ib_velocity_state(
        d_ctx, schedule, u_new_idx, d_ctx.u_current_idx, d_object_name + "::apply()");

    set_staggered_stokes_ib_velocity_bdry_state(d_ctx, d_ctx.u_idx, false);
    d_ctx.ib_implicit_ops->interpolateVelocity(
        d_ctx.u_idx, d_ctx.u_synch_scheds, d_ctx.u_ghost_fill_scheds, schedule.velocity_time);

    advance_staggered_stokes_ib_strategy(
        *d_ctx.ib_implicit_ops, d_ctx.time_stepping_type, current_time, new_time, d_object_name + "::apply()");

    d_ctx.ib_implicit_ops->computeLagrangianForce(schedule.force_time);
    d_ctx.hier_velocity_data_ops->setToScalar(d_ctx.f_idx, 0.0, /*interior_only*/ false);
    set_staggered_stokes_ib_velocity_bdry_state(d_ctx, d_ctx.f_idx, true);
    d_ctx.ib_implicit_ops->spreadForce(
        d_ctx.f_idx, d_ctx.u_phys_bdry_op, d_ctx.f_prolongation_scheds, schedule.force_time);
    d_ctx.hier_velocity_data_ops->axpy(f_u_idx, -schedule.nonlinear_force_scale, d_ctx.f_idx, f_u_idx);

    return;
} // apply

void
StaggeredStokesIBOperator::applyAdd(SAMRAIVectorReal<NDIM, double>& x,
                                    SAMRAIVectorReal<NDIM, double>& y,
                                    SAMRAIVectorReal<NDIM, double>& z)
{
    GeneralOperator::applyAdd(x, y, z);
    return;
} // applyAdd

void
StaggeredStokesIBOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                   const SAMRAIVectorReal<NDIM, double>& out)
{
    if (getIsInitialized()) deallocateOperatorState();

    if (!d_ctx.ib_implicit_ops || !d_ctx.stokes_op)
    {
        TBOX_ERROR(d_object_name << "::initializeOperatorState(): missing operator dependencies\n");
    }

    d_ctx.ib_implicit_ops->setUseFixedLEOperators(d_ctx.use_fixed_le_operators);
    d_ctx.stokes_op->initializeOperatorState(in, out);
    GeneralOperator::initializeOperatorState(in, out);
    return;
} // initializeOperatorState

void
StaggeredStokesIBOperator::deallocateOperatorState()
{
    if (d_ctx.stokes_op) d_ctx.stokes_op->deallocateOperatorState();
    GeneralOperator::deallocateOperatorState();
    return;
} // deallocateOperatorState

void
StaggeredStokesIBOperator::modifyRhsForBcs(SAMRAIVectorReal<NDIM, double>& y)
{
    prepare_staggered_stokes_ib_bc_forwarding(
        d_ctx, getTimeInterval(), getSolutionTime(), getHomogeneousBc(), d_object_name + "::modifyRhsForBcs()");
    d_ctx.stokes_op->modifyRhsForBcs(y);
    return;
} // modifyRhsForBcs

void
StaggeredStokesIBOperator::imposeSolBcs(SAMRAIVectorReal<NDIM, double>& u)
{
    prepare_staggered_stokes_ib_bc_forwarding(
        d_ctx, getTimeInterval(), getSolutionTime(), getHomogeneousBc(), d_object_name + "::imposeSolBcs()");
    d_ctx.stokes_op->imposeSolBcs(u);
    return;
} // imposeSolBcs

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
