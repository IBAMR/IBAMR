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
#include <ibamr/StaggeredStokesIBJacobianOperator.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/private/StaggeredStokesIBTimeSteppingUtilities.h>

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SAMRAIVectorReal.h>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesIBJacobianOperator::StaggeredStokesIBJacobianOperator(const std::string& object_name)
    : JacobianOperator(object_name)
{
    // intentionally blank
    return;
} // StaggeredStokesIBJacobianOperator

StaggeredStokesIBJacobianOperator::~StaggeredStokesIBJacobianOperator()
{
    deallocateOperatorState();
    if (d_SAJ_mat)
    {
        PetscErrorCode ierr = MatDestroy(&d_SAJ_mat);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // ~StaggeredStokesIBJacobianOperator

void
StaggeredStokesIBJacobianOperator::setOperatorContext(const StaggeredStokesIBOperator::Context& ctx)
{
    if (getIsInitialized())
    {
        TBOX_ERROR(d_object_name << "::setOperatorContext():\n"
                                 << "  deallocate operator state before replacing the Context.\n");
    }
    d_ctx = ctx;
    return;
} // setOperatorContext

void
StaggeredStokesIBJacobianOperator::setIBCouplingJacobian(Mat SAJ_mat)
{
    if (d_SAJ_mat == SAJ_mat)
    {
        return;
    }
    if (getIsInitialized())
    {
        validateContext(SAJ_mat != nullptr);
    }
    if (d_SAJ_mat)
    {
        PetscErrorCode ierr = MatDestroy(&d_SAJ_mat);
        IBTK_CHKERRQ(ierr);
    }
    // The work vectors were created from the previous matrix.
    if (d_input_vec)
    {
        PetscErrorCode ierr = VecDestroy(&d_input_vec);
        IBTK_CHKERRQ(ierr);
    }
    if (d_output_vec)
    {
        PetscErrorCode ierr = VecDestroy(&d_output_vec);
        IBTK_CHKERRQ(ierr);
    }
    d_SAJ_mat = SAJ_mat;
    if (d_SAJ_mat)
    {
        PetscErrorCode ierr = PetscObjectReference(reinterpret_cast<PetscObject>(d_SAJ_mat));
        IBTK_CHKERRQ(ierr);
    }
    return;
} // setIBCouplingJacobian

void
StaggeredStokesIBJacobianOperator::formJacobian(SAMRAIVectorReal<NDIM, double>& x)
{
    if (!getIsInitialized())
    {
        TBOX_ERROR(d_object_name << "::formJacobian():\n  the operator state must be initialized.\n");
    }
    if (!d_base_vector)
    {
        d_base_vector = x.cloneVector(d_object_name + "::base_vector");
        d_base_vector->allocateVectorData();
    }
    d_base_vector->copyVector(Pointer<SAMRAIVectorReal<NDIM, double>>(&x, false));
    if (d_SAJ_mat)
    {
        // The supplied matrix already defines the coupling. The base vector above may differ from the one the
        // strategy last linearized at, so a later formJacobian() with no matrix installed must redo that work.
        d_strategy_linearization_formed = false;
        return;
    }

    const double current_time = getTimeInterval().first;
    const double new_time = getTimeInterval().second;
    const detail::StaggeredStokesIBTimeStepParameters step_parameters =
        detail::get_staggered_stokes_ib_time_step_parameters(
            d_ctx.time_stepping_type, current_time, new_time, d_object_name + "::formJacobian()");

    const int u_new_idx = x.getComponentDescriptorIndex(0);

    if (d_ctx.time_stepping_type == MIDPOINT_RULE)
    {
        d_ctx.hier_velocity_data_ops->linearSum(d_ctx.u_idx, 0.5, u_new_idx, 0.5, d_ctx.u_current_idx);
    }
    else
    {
        d_ctx.hier_velocity_data_ops->copyData(d_ctx.u_idx, u_new_idx);
    }

    if (!d_solver_X || !d_solver_X0)
    {
        d_ctx.ib_implicit_ops->createSolverVecs(&d_solver_X, &d_solver_X0);
    }
    d_ctx.ib_implicit_ops->setupSolverVecs(&d_solver_X0, nullptr);

    if (d_ctx.u_phys_bdry_op)
    {
        d_ctx.u_phys_bdry_op->setPatchDataIndex(d_ctx.u_idx);
        d_ctx.u_phys_bdry_op->setHomogeneousBc(false);
    }

    if (d_ctx.time_stepping_type == BACKWARD_EULER)
    {
        // Backward Euler evaluates the force purely at the new position (force_position_fraction is
        // trivially 1), so it uses the strategy's own dedicated linearized-residual interface
        // (interpolateLinearizedVelocity()/computeLinearizedResidual()) instead of the general residual
        // interface the else branch below uses, which the other rules need for their force-position blend.
        //
        // Interpolate the physical base, including inhomogeneous boundary data.
        // Negating u before interpolation would not negate the boundary contribution.
        d_ctx.ib_implicit_ops->interpolateLinearizedVelocity(
            d_ctx.u_idx, d_ctx.u_synch_scheds, d_ctx.u_ghost_fill_scheds, step_parameters.evaluation_time);
        d_ctx.ib_implicit_ops->computeLinearizedResidual(d_solver_X0, d_solver_X);
        // The position residual is X_current - dt*U; reflect it to obtain X_new.
        PetscErrorCode ierr = VecAXPBY(d_solver_X, 2.0, -1.0, d_solver_X0);
        IBTK_CHKERRQ(ierr);
    }
    else
    {
        // Populate the same nonlinear velocity state as apply() before changing
        // positions. In particular, trapezoidal stepping retains the stored
        // Lagrangian current velocity, not a new interpolation at the endpoint.
        d_ctx.ib_implicit_ops->interpolateVelocity(
            d_ctx.u_idx, d_ctx.u_synch_scheds, d_ctx.u_ghost_fill_scheds, step_parameters.evaluation_time);
        // The position residual at X_current is -dt*U_half. Recover and restore
        // the endpoint, then select the force evaluation position separately.
        d_ctx.ib_implicit_ops->setUpdatedPosition(d_solver_X0);
        d_ctx.ib_implicit_ops->computeResidual(d_solver_X);
        PetscErrorCode ierr = VecAYPX(d_solver_X, -1.0, d_solver_X0);
        IBTK_CHKERRQ(ierr);
        d_ctx.ib_implicit_ops->setUpdatedPosition(d_solver_X);
        const double fraction = step_parameters.force_position_fraction;
        ierr = VecAXPBY(d_solver_X, 1.0 - fraction, fraction, d_solver_X0);
        IBTK_CHKERRQ(ierr);
    }
    d_ctx.ib_implicit_ops->setLinearizedPosition(d_solver_X, step_parameters.evaluation_time);
    d_strategy_linearization_formed = true;
    return;
} // formJacobian

Pointer<SAMRAIVectorReal<NDIM, double>>
StaggeredStokesIBJacobianOperator::getBaseVector() const
{
    return d_base_vector;
} // getBaseVector

void
StaggeredStokesIBJacobianOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    if (!getIsInitialized())
    {
        TBOX_ERROR(d_object_name << "::apply():\n  the operator state must be initialized.\n");
    }
    if (d_SAJ_mat)
    {
        if (x.getCoarsestLevelNumber() != x.getFinestLevelNumber())
        {
            TBOX_ERROR(d_object_name << "::apply():\n  a supplied coupling matrix requires single-level vectors.\n");
        }

        Pointer<PatchLevel<NDIM>> level = x.getPatchHierarchy()->getPatchLevel(x.getCoarsestLevelNumber());
        PetscErrorCode ierr = 0;
        if (!d_input_vec || !d_output_vec)
        {
            ierr = MatCreateVecs(d_SAJ_mat, &d_input_vec, &d_output_vec);
            IBTK_CHKERRQ(ierr);
        }

        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(d_input_vec,
                                                              x.getComponentDescriptorIndex(0),
                                                              d_ctx.u_dof_index_idx,
                                                              x.getComponentDescriptorIndex(1),
                                                              d_ctx.p_dof_index_idx,
                                                              level);
        d_ctx.stokes_op->setTimeInterval(getTimeInterval().first, getTimeInterval().second);
        d_ctx.stokes_op->setSolutionTime(getSolutionTime());
        d_ctx.stokes_op->setHomogeneousBc(true);
        d_ctx.stokes_op->apply(x, y);
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(d_output_vec,
                                                              y.getComponentDescriptorIndex(0),
                                                              d_ctx.u_dof_index_idx,
                                                              y.getComponentDescriptorIndex(1),
                                                              d_ctx.p_dof_index_idx,
                                                              level);
        // The supplied coupling matrix already includes its sign and time-step
        // factors: add it directly to the Stokes action in coupled DOF ordering.
        ierr = MatMultAdd(d_SAJ_mat, d_input_vec, d_output_vec, d_output_vec);
        IBTK_CHKERRQ(ierr);
        StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(d_output_vec,
                                                                y.getComponentDescriptorIndex(0),
                                                                d_ctx.u_dof_index_idx,
                                                                y.getComponentDescriptorIndex(1),
                                                                d_ctx.p_dof_index_idx,
                                                                level,
                                                                nullptr,
                                                                nullptr);
        return;
    }

    if (!d_strategy_linearization_formed)
    {
        TBOX_ERROR(d_object_name << "::apply():\n"
                                 << "  formJacobian() must be called before the strategy action is applied.\n");
    }
    const double current_time = getTimeInterval().first;
    const double new_time = getTimeInterval().second;
    const detail::StaggeredStokesIBTimeStepParameters step_parameters =
        detail::get_staggered_stokes_ib_time_step_parameters(
            d_ctx.time_stepping_type, current_time, new_time, d_object_name + "::apply()");

    const int u_idx = x.getComponentDescriptorIndex(0);
    const int f_u_idx = y.getComponentDescriptorIndex(0);

    d_ctx.stokes_op->setTimeInterval(current_time, new_time);
    d_ctx.stokes_op->setSolutionTime(getSolutionTime());
    d_ctx.stokes_op->setHomogeneousBc(true);
    d_ctx.stokes_op->apply(x, y);

    if (!d_solver_X || !d_solver_X0)
    {
        d_ctx.ib_implicit_ops->createSolverVecs(&d_solver_X, &d_solver_X0);
    }
    // With zero position input, the linearized position residual is -dt*J*u.
    // Interpolating -beta*u therefore produces dt*beta*J*u below, where beta
    // is 1 for backward Euler and 1/2 for trapezoidal/midpoint.
    d_ctx.ib_implicit_ops->setupSolverVecs(nullptr, &d_solver_X0);

    d_ctx.hier_velocity_data_ops->scale(d_ctx.u_idx, -step_parameters.jacobian_force_scale, u_idx);
    if (d_ctx.u_phys_bdry_op)
    {
        d_ctx.u_phys_bdry_op->setPatchDataIndex(d_ctx.u_idx);
        d_ctx.u_phys_bdry_op->setHomogeneousBc(true);
    }
    d_ctx.ib_implicit_ops->interpolateLinearizedVelocity(
        d_ctx.u_idx, d_ctx.u_synch_scheds, d_ctx.u_ghost_fill_scheds, step_parameters.evaluation_time);
    d_ctx.ib_implicit_ops->computeLinearizedResidual(d_solver_X0, d_solver_X);

    // Apply the force derivative K at the position cached by formJacobian(),
    // spread it, and subtract beta times the result. With fixed interpolation
    // J and spreading S, the added momentum block is -dt*beta^2*S*K*J.
    // For trapezoidal the two halves are the position and force weights; for
    // midpoint their product accounts for both velocity and position averaging.
    // The Stokes pressure/divergence action is unchanged.
    d_ctx.ib_implicit_ops->computeLinearizedLagrangianForce(d_solver_X, step_parameters.evaluation_time);
    d_ctx.hier_velocity_data_ops->setToScalar(d_ctx.f_idx, 0.0, /*interior_only*/ false);
    if (d_ctx.u_phys_bdry_op)
    {
        d_ctx.u_phys_bdry_op->setPatchDataIndex(d_ctx.f_idx);
        d_ctx.u_phys_bdry_op->setHomogeneousBc(true);
    }
    d_ctx.ib_implicit_ops->spreadLinearizedForce(
        d_ctx.f_idx, d_ctx.u_phys_bdry_op, d_ctx.f_prolongation_scheds, step_parameters.evaluation_time);
    d_ctx.hier_velocity_data_ops->axpy(f_u_idx, -step_parameters.jacobian_force_scale, d_ctx.f_idx, f_u_idx);

    return;
} // apply

void
StaggeredStokesIBJacobianOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                           const SAMRAIVectorReal<NDIM, double>& out)
{
    if (getIsInitialized())
    {
        deallocateOperatorState();
    }

    validateContext(d_SAJ_mat != nullptr);
    if (!d_ctx.stokes_op->getIsInitialized())
    {
        TBOX_ERROR(d_object_name << "::initializeOperatorState():\n"
                                 << "  the owner of the Stokes operator must initialize it first.\n");
    }
    JacobianOperator::initializeOperatorState(in, out);
    return;
} // initializeOperatorState

void
StaggeredStokesIBJacobianOperator::deallocateOperatorState()
{
    d_strategy_linearization_formed = false;
    if (d_base_vector)
    {
        free_vector_components(*d_base_vector);
        d_base_vector.setNull();
    }
    if (d_input_vec)
    {
        PetscErrorCode ierr = VecDestroy(&d_input_vec);
        IBTK_CHKERRQ(ierr);
    }
    if (d_output_vec)
    {
        PetscErrorCode ierr = VecDestroy(&d_output_vec);
        IBTK_CHKERRQ(ierr);
    }
    if (d_solver_X)
    {
        PetscErrorCode ierr = VecDestroy(&d_solver_X);
        IBTK_CHKERRQ(ierr);
    }
    if (d_solver_X0)
    {
        PetscErrorCode ierr = VecDestroy(&d_solver_X0);
        IBTK_CHKERRQ(ierr);
    }
    JacobianOperator::deallocateOperatorState();
    return;
} // deallocateOperatorState

void
StaggeredStokesIBJacobianOperator::modifyRhsForBcs(SAMRAIVectorReal<NDIM, double>& y)
{
    if (!d_ctx.stokes_op)
    {
        TBOX_ERROR(d_object_name << "::modifyRhsForBcs():\n  missing Stokes operator.\n");
    }
    d_ctx.stokes_op->setTimeInterval(getTimeInterval().first, getTimeInterval().second);
    d_ctx.stokes_op->setSolutionTime(getSolutionTime());
    d_ctx.stokes_op->setHomogeneousBc(getHomogeneousBc());
    d_ctx.stokes_op->modifyRhsForBcs(y);
    return;
} // modifyRhsForBcs

void
StaggeredStokesIBJacobianOperator::imposeSolBcs(SAMRAIVectorReal<NDIM, double>& u)
{
    if (!d_ctx.stokes_op)
    {
        TBOX_ERROR(d_object_name << "::imposeSolBcs():\n  missing Stokes operator.\n");
    }
    d_ctx.stokes_op->setTimeInterval(getTimeInterval().first, getTimeInterval().second);
    d_ctx.stokes_op->setSolutionTime(getSolutionTime());
    d_ctx.stokes_op->setHomogeneousBc(getHomogeneousBc());
    d_ctx.stokes_op->imposeSolBcs(u);
    return;
} // imposeSolBcs

/////////////////////////////// PRIVATE //////////////////////////////////////

void
StaggeredStokesIBJacobianOperator::validateContext(bool supplied_matrix) const
{
    if (!d_ctx.stokes_op)
    {
        TBOX_ERROR(d_object_name << "::validateContext():\n  the operator Context requires a Stokes operator.\n");
    }
    if (supplied_matrix)
    {
        if (d_ctx.u_dof_index_idx < 0 || d_ctx.p_dof_index_idx < 0)
        {
            TBOX_ERROR(
                d_object_name << "::validateContext():\n"
                              << "  a supplied coupling matrix requires velocity and pressure DOF-index data.\n");
        }
    }
    else
    {
        if (!d_ctx.ib_implicit_ops || !d_ctx.hier_velocity_data_ops)
        {
            TBOX_ERROR(
                d_object_name << "::validateContext():\n"
                              << "  the strategy action requires an IB strategy and velocity data operations.\n");
        }
        if (d_ctx.u_idx < 0 || d_ctx.f_idx < 0 || d_ctx.u_current_idx < 0)
        {
            TBOX_ERROR(
                d_object_name << "::validateContext():\n"
                              << "  the strategy action requires velocity, force and current-velocity data indices.\n");
        }
        detail::require_staggered_stokes_ib_time_stepping_type(d_ctx.time_stepping_type,
                                                               d_object_name + "::validateContext()");
    }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
