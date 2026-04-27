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

#include <ibtk/RobinPhysBdryPatchStrategy.h>

#include <tbox/Pointer.h>

#include <PatchHierarchy.h>
#include <SAMRAIVectorReal.h>

#include <limits>
#include <utility>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
void
set_stokes_operator_times(const Pointer<StaggeredStokesOperator>& stokes_op,
                          const std::pair<double, double>& time_interval,
                          const double solution_time)
{
    if (!stokes_op) return;
    stokes_op->setTimeInterval(time_interval.first, time_interval.second);
    stokes_op->setSolutionTime(solution_time);
}
} // namespace

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
    const double half_time = current_time + 0.5 * getDt();

    const int u_current_idx = d_ctx.u_current_idx;
    const int u_new_idx = x.getComponentDescriptorIndex(0);
    const int f_u_idx = y.getComponentDescriptorIndex(0);

    set_stokes_operator_times(d_ctx.stokes_op, getTimeInterval(), getSolutionTime());
    d_ctx.stokes_op->setHomogeneousBc(true);
    d_ctx.stokes_op->apply(x, y);

    double velocity_time = std::numeric_limits<double>::quiet_NaN();
    switch (d_ctx.time_stepping_type)
    {
    case BACKWARD_EULER:
    case TRAPEZOIDAL_RULE:
        d_ctx.hier_velocity_data_ops->copyData(d_ctx.u_idx, u_new_idx);
        velocity_time = new_time;
        break;
    case MIDPOINT_RULE:
        d_ctx.hier_velocity_data_ops->linearSum(d_ctx.u_idx, 0.5, u_new_idx, 0.5, u_current_idx);
        velocity_time = half_time;
        break;
    default:
        TBOX_ERROR(d_object_name << "::apply(): unsupported time stepping type\n");
    }

    if (d_ctx.u_phys_bdry_op)
    {
        d_ctx.u_phys_bdry_op->setPatchDataIndex(d_ctx.u_idx);
        d_ctx.u_phys_bdry_op->setHomogeneousBc(false);
    }
    d_ctx.ib_implicit_ops->interpolateVelocity(
        d_ctx.u_idx, d_ctx.u_synch_scheds, d_ctx.u_ghost_fill_scheds, velocity_time);

    switch (d_ctx.time_stepping_type)
    {
    case BACKWARD_EULER:
        d_ctx.ib_implicit_ops->backwardEulerStep(current_time, new_time);
        break;
    case TRAPEZOIDAL_RULE:
        d_ctx.ib_implicit_ops->trapezoidalStep(current_time, new_time);
        break;
    case MIDPOINT_RULE:
        d_ctx.ib_implicit_ops->midpointStep(current_time, new_time);
        break;
    default:
        TBOX_ERROR(d_object_name << "::apply(): unsupported time stepping type\n");
    }
    double force_time = std::numeric_limits<double>::quiet_NaN();
    double kappa = std::numeric_limits<double>::quiet_NaN();
    switch (d_ctx.time_stepping_type)
    {
    case BACKWARD_EULER:
        force_time = new_time;
        kappa = 1.0;
        break;
    case TRAPEZOIDAL_RULE:
        force_time = new_time;
        kappa = 0.5;
        break;
    case MIDPOINT_RULE:
        force_time = half_time;
        kappa = 1.0;
        break;
    default:
        TBOX_ERROR(d_object_name << "::apply(): unsupported time stepping type\n");
    }

    d_ctx.ib_implicit_ops->computeLagrangianForce(force_time);
    d_ctx.hier_velocity_data_ops->setToScalar(d_ctx.f_idx, 0.0, /*interior_only*/ false);
    if (d_ctx.u_phys_bdry_op)
    {
        d_ctx.u_phys_bdry_op->setPatchDataIndex(d_ctx.f_idx);
        d_ctx.u_phys_bdry_op->setHomogeneousBc(true);
    }
    d_ctx.ib_implicit_ops->spreadForce(d_ctx.f_idx, d_ctx.u_phys_bdry_op, d_ctx.f_prolongation_scheds, force_time);
    d_ctx.hier_velocity_data_ops->axpy(f_u_idx, -kappa, d_ctx.f_idx, f_u_idx);

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
    if (!d_ctx.stokes_op)
    {
        TBOX_ERROR(d_object_name << "::modifyRhsForBcs(): missing Stokes operator\n");
    }
    set_stokes_operator_times(d_ctx.stokes_op, getTimeInterval(), getSolutionTime());
    d_ctx.stokes_op->setHomogeneousBc(getHomogeneousBc());
    d_ctx.stokes_op->modifyRhsForBcs(y);
    return;
} // modifyRhsForBcs

void
StaggeredStokesIBOperator::imposeSolBcs(SAMRAIVectorReal<NDIM, double>& u)
{
    if (!d_ctx.stokes_op)
    {
        TBOX_ERROR(d_object_name << "::imposeSolBcs(): missing Stokes operator\n");
    }
    set_stokes_operator_times(d_ctx.stokes_op, getTimeInterval(), getSolutionTime());
    d_ctx.stokes_op->setHomogeneousBc(getHomogeneousBc());
    d_ctx.stokes_op->imposeSolBcs(u);
    return;
} // imposeSolBcs

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
