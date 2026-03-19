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

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>

#include <tbox/Pointer.h>

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

StaggeredStokesIBJacobianOperator::StaggeredStokesIBJacobianOperator(const std::string& object_name)
    : JacobianOperator(object_name)
{
    // intentionally blank
    return;
} // StaggeredStokesIBJacobianOperator

StaggeredStokesIBJacobianOperator::~StaggeredStokesIBJacobianOperator()
{
    deallocateOperatorState();
    return;
} // ~StaggeredStokesIBJacobianOperator

void
StaggeredStokesIBJacobianOperator::setOperatorContext(const StaggeredStokesIBOperatorContext& ctx)
{
    d_ctx = ctx;
    return;
} // setOperatorContext

void
StaggeredStokesIBJacobianOperator::setIBCouplingJacobian(Mat& SAJ_mat)
{
    PetscErrorCode ierr = MatDestroy(&d_SAJ_mat);
    IBTK_CHKERRQ(ierr);
    d_SAJ_mat = SAJ_mat;
    if (d_SAJ_mat)
    {
        ierr = PetscObjectReference(reinterpret_cast<PetscObject>(d_SAJ_mat));
        IBTK_CHKERRQ(ierr);
    }
    return;
} // setIBCouplingJacobian

void
StaggeredStokesIBJacobianOperator::formJacobian(SAMRAIVectorReal<NDIM, double>& x)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(getIsInitialized());
#endif
    if (!d_ctx.ib_implicit_ops || !d_ctx.stokes_op || !d_ctx.hier_velocity_data_ops)
    {
        TBOX_ERROR(d_object_name << "::formJacobian(): incomplete operator context\n");
    }
    if (d_ctx.u_idx == IBTK::invalid_index)
    {
        TBOX_ERROR(d_object_name << "::formJacobian(): invalid velocity scratch index\n");
    }
    if (d_ctx.u_current_idx == IBTK::invalid_index)
    {
        TBOX_ERROR(d_object_name << "::formJacobian(): invalid current velocity data index\n");
    }

    if (!d_base_vector)
    {
        d_base_vector = x.cloneVector(d_object_name + "::base_vector");
        d_base_vector->allocateVectorData();
    }
    d_base_vector->copyVector(Pointer<SAMRAIVectorReal<NDIM, double>>(&x, false));

    const double current_time = getTimeInterval().first;
    const double new_time = getTimeInterval().second;
    const double half_time = current_time + 0.5 * getDt();

    const int u_current_idx = d_ctx.u_current_idx;
    const int u_new_idx = x.getComponentDescriptorIndex(0);

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
        TBOX_ERROR(d_object_name << "::formJacobian(): unsupported time stepping type\n");
    }

    Vec X_new = nullptr;
    Vec X0 = nullptr;
    d_ctx.ib_implicit_ops->createSolverVecs(&X_new, &X0);
    d_ctx.ib_implicit_ops->setupSolverVecs(&X0, nullptr);

    d_ctx.hier_velocity_data_ops->scale(d_ctx.u_idx, -1.0, d_ctx.u_idx);
    if (d_ctx.u_phys_bdry_op)
    {
        d_ctx.u_phys_bdry_op->setPatchDataIndex(d_ctx.u_idx);
        d_ctx.u_phys_bdry_op->setHomogeneousBc(true);
    }
    d_ctx.ib_implicit_ops->interpolateLinearizedVelocity(
        d_ctx.u_idx, d_ctx.u_synch_scheds, d_ctx.u_ghost_fill_scheds, velocity_time);
    d_ctx.ib_implicit_ops->computeLinearizedResidual(X0, X_new);
    d_ctx.ib_implicit_ops->setLinearizedPosition(X_new, velocity_time);

    PetscErrorCode ierr = VecDestroy(&X_new);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&X0);
    IBTK_CHKERRQ(ierr);
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
#if !defined(NDEBUG)
    TBOX_ASSERT(getIsInitialized());
#endif
    if (d_SAJ_mat)
    {
        if (!d_ctx.patch_level)
        {
            TBOX_ERROR(d_object_name << "::apply(): SAJ apply path requires a valid patch level\n");
        }
        if (d_ctx.u_dof_index_idx == IBTK::invalid_index || d_ctx.p_dof_index_idx == IBTK::invalid_index)
        {
            TBOX_ERROR(d_object_name << "::apply(): SAJ apply path requires valid DOF-index patch data indices\n");
        }
        if (x.getCoarsestLevelNumber() != x.getFinestLevelNumber())
        {
            TBOX_ERROR(d_object_name << "::apply(): SAJ apply path supports single-level vectors only\n");
        }
        if (!d_ctx.stokes_op)
        {
            TBOX_ERROR(d_object_name << "::apply(): SAJ apply path requires Stokes operator\n");
        }

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
                                                              d_ctx.patch_level);
        set_stokes_operator_times(d_ctx.stokes_op, getTimeInterval(), getSolutionTime());
        d_ctx.stokes_op->setHomogeneousBc(true);
        d_ctx.stokes_op->apply(x, y);
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(d_output_vec,
                                                              y.getComponentDescriptorIndex(0),
                                                              d_ctx.u_dof_index_idx,
                                                              y.getComponentDescriptorIndex(1),
                                                              d_ctx.p_dof_index_idx,
                                                              d_ctx.patch_level);
        ierr = MatMultAdd(d_SAJ_mat, d_input_vec, d_output_vec, d_output_vec);
        IBTK_CHKERRQ(ierr);
        StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(d_output_vec,
                                                                y.getComponentDescriptorIndex(0),
                                                                d_ctx.u_dof_index_idx,
                                                                y.getComponentDescriptorIndex(1),
                                                                d_ctx.p_dof_index_idx,
                                                                d_ctx.patch_level,
                                                                nullptr,
                                                                nullptr);
        return;
    }

    if (!d_ctx.ib_implicit_ops || !d_ctx.stokes_op || !d_ctx.hier_velocity_data_ops)
    {
        TBOX_ERROR(d_object_name << "::apply(): incomplete operator context\n");
    }
    if (d_ctx.u_idx == IBTK::invalid_index || d_ctx.f_idx == IBTK::invalid_index)
    {
        TBOX_ERROR(d_object_name << "::apply(): invalid scratch data indices\n");
    }

    const double current_time = getTimeInterval().first;
    const double new_time = getTimeInterval().second;
    const double half_time = current_time + 0.5 * getDt();

    const int u_idx = x.getComponentDescriptorIndex(0);
    const int f_u_idx = y.getComponentDescriptorIndex(0);

    set_stokes_operator_times(d_ctx.stokes_op, getTimeInterval(), getSolutionTime());
    d_ctx.stokes_op->setHomogeneousBc(true);
    d_ctx.stokes_op->apply(x, y);

    double force_time = std::numeric_limits<double>::quiet_NaN();
    double velocity_time = std::numeric_limits<double>::quiet_NaN();
    double kappa = std::numeric_limits<double>::quiet_NaN();
    switch (d_ctx.time_stepping_type)
    {
    case BACKWARD_EULER:
        force_time = new_time;
        velocity_time = new_time;
        kappa = 1.0;
        break;
    case TRAPEZOIDAL_RULE:
        force_time = new_time;
        velocity_time = new_time;
        kappa = 0.5;
        break;
    case MIDPOINT_RULE:
        force_time = half_time;
        velocity_time = half_time;
        kappa = 0.5;
        break;
    default:
        TBOX_ERROR(d_object_name << "::apply(): unsupported time stepping type\n");
    }

    Vec X = nullptr;
    Vec X0 = nullptr;
    d_ctx.ib_implicit_ops->createSolverVecs(&X, &X0);
    d_ctx.ib_implicit_ops->setupSolverVecs(nullptr, &X0);

    d_ctx.hier_velocity_data_ops->scale(d_ctx.u_idx, -kappa, u_idx);
    if (d_ctx.u_phys_bdry_op)
    {
        d_ctx.u_phys_bdry_op->setPatchDataIndex(d_ctx.u_idx);
        d_ctx.u_phys_bdry_op->setHomogeneousBc(true);
    }
    d_ctx.ib_implicit_ops->interpolateLinearizedVelocity(
        d_ctx.u_idx, d_ctx.u_synch_scheds, d_ctx.u_ghost_fill_scheds, velocity_time);
    d_ctx.ib_implicit_ops->computeLinearizedResidual(X0, X);

    d_ctx.ib_implicit_ops->computeLinearizedLagrangianForce(X, force_time);
    d_ctx.hier_velocity_data_ops->setToScalar(d_ctx.f_idx, 0.0, /*interior_only*/ false);
    if (d_ctx.u_phys_bdry_op)
    {
        d_ctx.u_phys_bdry_op->setPatchDataIndex(d_ctx.f_idx);
        d_ctx.u_phys_bdry_op->setHomogeneousBc(true);
    }
    d_ctx.ib_implicit_ops->spreadLinearizedForce(
        d_ctx.f_idx, d_ctx.u_phys_bdry_op, d_ctx.f_prolongation_scheds, force_time);
    d_ctx.hier_velocity_data_ops->axpy(f_u_idx, -kappa, d_ctx.f_idx, f_u_idx);

    PetscErrorCode ierr = VecDestroy(&X);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&X0);
    IBTK_CHKERRQ(ierr);
    return;
} // apply

void
StaggeredStokesIBJacobianOperator::applyAdd(SAMRAIVectorReal<NDIM, double>& x,
                                            SAMRAIVectorReal<NDIM, double>& y,
                                            SAMRAIVectorReal<NDIM, double>& z)
{
    GeneralOperator::applyAdd(x, y, z);
    return;
} // applyAdd

void
StaggeredStokesIBJacobianOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                           const SAMRAIVectorReal<NDIM, double>& out)
{
    if (getIsInitialized()) deallocateOperatorState();

    if (!d_ctx.ib_implicit_ops || !d_ctx.stokes_op)
    {
        TBOX_ERROR(d_object_name << "::initializeOperatorState(): missing operator dependencies\n");
    }

    d_ctx.ib_implicit_ops->setUseFixedLEOperators(true);
    d_ctx.ib_implicit_ops->updateFixedLEOperators();
    d_ctx.stokes_op->initializeOperatorState(in, out);
    JacobianOperator::initializeOperatorState(in, out);
    return;
} // initializeOperatorState

void
StaggeredStokesIBJacobianOperator::deallocateOperatorState()
{
    if (d_base_vector)
    {
        d_base_vector->deallocateVectorData();
        d_base_vector.setNull();
    }
    if (d_ctx.stokes_op) d_ctx.stokes_op->deallocateOperatorState();
    PetscErrorCode ierr = VecDestroy(&d_input_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&d_output_vec);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&d_SAJ_mat);
    IBTK_CHKERRQ(ierr);
    JacobianOperator::deallocateOperatorState();
    return;
} // deallocateOperatorState

void
StaggeredStokesIBJacobianOperator::modifyRhsForBcs(SAMRAIVectorReal<NDIM, double>& y)
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
StaggeredStokesIBJacobianOperator::imposeSolBcs(SAMRAIVectorReal<NDIM, double>& u)
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
