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
#include <ibamr/StaggeredStokesIBJacobianFACPreconditioner.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/FACPreconditionerStrategy.h>
#include <ibtk/SAMRAIScopedVectorCopy.h>
#include <ibtk/SAMRAIScopedVectorDuplicate.h>
#include <ibtk/ibtk_enums.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesIBJacobianFACPreconditioner::StaggeredStokesIBJacobianFACPreconditioner(
    const std::string& object_name,
    Pointer<IBTK::FACPreconditionerStrategy> fac_strategy,
    Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : FACPreconditioner(object_name, fac_strategy, input_db, default_options_prefix)
{
    // intentionally blank
    return;
} // StaggeredStokesIBJacobianFACPreconditioner

void
StaggeredStokesIBJacobianFACPreconditioner::setVelocityPoissonSpecifications(
    const PoissonSpecifications& U_problem_coefs)
{
    StaggeredStokesSolver::setVelocityPoissonSpecifications(U_problem_coefs);
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setVelocityPoissonSpecifications(U_problem_coefs);
    return;
} // setVelocityPoissonSpecifications

void
StaggeredStokesIBJacobianFACPreconditioner::setComponentsHaveNullSpace(const bool has_velocity_nullspace,
                                                                       const bool has_pressure_nullspace)
{
    StaggeredStokesSolver::setComponentsHaveNullSpace(has_velocity_nullspace, has_pressure_nullspace);
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy)
    {
        fac_strategy->setComponentsHaveNullSpace(d_has_velocity_nullspace, d_has_pressure_nullspace);
    }
    return;
} // setComponentsHaveNullSpace

void
StaggeredStokesIBJacobianFACPreconditioner::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
    RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
    StaggeredStokesSolver::setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    return;
} // setPhysicalBcCoefs

void
StaggeredStokesIBJacobianFACPreconditioner::setPhysicalBoundaryHelper(
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
    StaggeredStokesSolver::setPhysicalBoundaryHelper(bc_helper);
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setPhysicalBoundaryHelper(d_bc_helper);
    return;
} // setPhysicalBoundaryHelper

void
StaggeredStokesIBJacobianFACPreconditioner::setIBTimeSteppingType(const TimeSteppingType time_stepping_type)
{
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setIBTimeSteppingType(time_stepping_type);
    return;
} // setIBTimeSteppingType

void
StaggeredStokesIBJacobianFACPreconditioner::setIBForceJacobian(Mat& A_mat)
{
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setIBForceJacobian(A_mat);
    return;
} // setIBForceJacobian

void
StaggeredStokesIBJacobianFACPreconditioner::setIBInterpOp(Mat& J_mat)
{
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setIBInterpOp(J_mat);
    return;
} // setIBInterpOp

void
StaggeredStokesIBJacobianFACPreconditioner::setIBImplicitStrategy(Pointer<IBImplicitStrategy> ib_implicit_ops)
{
    d_ib_implicit_ops = ib_implicit_ops;
    return;
} // setIBImplicitStrategy

void
StaggeredStokesIBJacobianFACPreconditioner::initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                                                  const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b)
{
    if (d_ib_implicit_ops)
    {
        d_ib_implicit_ops->setUseFixedLEOperators(true);
        d_ib_implicit_ops->updateFixedLEOperators();
    }
    FACPreconditioner::initializeSolverState(x, b);
    return;
} // initializeSolverState

bool
StaggeredStokesIBJacobianFACPreconditioner::solveSystemWithCoarseLevelDiagnostics(
    SAMRAIVectorReal<NDIM, double>& x,
    SAMRAIVectorReal<NDIM, double>& b,
    SAMRAIVectorReal<NDIM, double>* coarse_correction_out,
    SAMRAIVectorReal<NDIM, double>* coarse_rhs_out,
    int diagnostic_level_num,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* pre_smooth_input_out_by_level,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* pre_smooth_solution_out_by_level,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* post_smooth_input_out_by_level,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* post_smooth_solution_out_by_level)
{
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    x.setToScalar(0.0, /*interior_only*/ false);
    if (coarse_correction_out) coarse_correction_out->setToScalar(0.0, /*interior_only*/ false);
    if (coarse_rhs_out) coarse_rhs_out->setToScalar(0.0, /*interior_only*/ false);
    if (pre_smooth_input_out_by_level)
    {
        for (auto* pre_smooth_input_out : *pre_smooth_input_out_by_level)
        {
            if (pre_smooth_input_out) pre_smooth_input_out->setToScalar(0.0, /*interior_only*/ false);
        }
    }
    if (pre_smooth_solution_out_by_level)
    {
        for (auto* pre_smooth_solution_out : *pre_smooth_solution_out_by_level)
        {
            if (pre_smooth_solution_out) pre_smooth_solution_out->setToScalar(0.0, /*interior_only*/ false);
        }
    }
    if (post_smooth_input_out_by_level)
    {
        for (auto* post_smooth_input_out : *post_smooth_input_out_by_level)
        {
            if (post_smooth_input_out) post_smooth_input_out->setToScalar(0.0, /*interior_only*/ false);
        }
    }
    if (post_smooth_solution_out_by_level)
    {
        for (auto* post_smooth_solution_out : *post_smooth_solution_out_by_level)
        {
            if (post_smooth_solution_out) post_smooth_solution_out->setToScalar(0.0, /*interior_only*/ false);
        }
    }

    const int diag_level = diagnostic_level_num >= 0 ? diagnostic_level_num : d_coarsest_ln;
    bool captured_diagnostics = false;
    if (d_cycle_type == V_CYCLE && d_num_pre_sweeps == 0)
    {
        FACVCycleNoPreSmoothingWithDiagnostics(
            x,
            b,
            d_finest_ln,
            diag_level,
            coarse_correction_out,
            coarse_rhs_out,
            captured_diagnostics,
            pre_smooth_input_out_by_level,
            pre_smooth_solution_out_by_level,
            post_smooth_input_out_by_level,
            post_smooth_solution_out_by_level);
    }
    else
    {
        SAMRAIScopedVectorCopy<double> f(b);
        SAMRAIScopedVectorDuplicate<double> r(b);
        switch (d_cycle_type)
        {
        case F_CYCLE:
            TBOX_ERROR(
                d_object_name << "::solveSystemWithCoarseLevelDiagnostics(): F_CYCLE diagnostics are not implemented."
                              << std::endl);
        case FMG_CYCLE:
            TBOX_ERROR(
                d_object_name << "::solveSystemWithCoarseLevelDiagnostics(): FMG_CYCLE diagnostics are not implemented."
                              << std::endl);
        case V_CYCLE:
            muCycleWithDiagnostics(
                x,
                f,
                r,
                d_finest_ln,
                1,
                diag_level,
                coarse_correction_out,
                coarse_rhs_out,
                captured_diagnostics,
                pre_smooth_input_out_by_level,
                pre_smooth_solution_out_by_level,
                post_smooth_input_out_by_level,
                post_smooth_solution_out_by_level);
            break;
        case W_CYCLE:
            muCycleWithDiagnostics(
                x,
                f,
                r,
                d_finest_ln,
                2,
                diag_level,
                coarse_correction_out,
                coarse_rhs_out,
                captured_diagnostics,
                pre_smooth_input_out_by_level,
                pre_smooth_solution_out_by_level,
                post_smooth_input_out_by_level,
                post_smooth_solution_out_by_level);
            break;
        default:
            TBOX_ERROR(d_object_name << "::solveSystemWithCoarseLevelDiagnostics():\n"
                                     << "  unsupported FAC cycle type: "
                                     << IBTK::enum_to_string<IBTK::MGCycleType>(d_cycle_type) << "." << std::endl);
        }
    }

    if (!captured_diagnostics && (coarse_correction_out || coarse_rhs_out))
    {
        TBOX_WARNING(d_object_name << "::solveSystemWithCoarseLevelDiagnostics():\n"
                                   << "  requested diagnostics were not captured at level " << diag_level << ".\n");
    }

    if (deallocate_after_solve) deallocateSolverState();
    return true;
} // solveSystemWithCoarseLevelDiagnostics

Pointer<StaggeredStokesIBLevelRelaxationFACOperator>
StaggeredStokesIBJacobianFACPreconditioner::getIBFACPreconditionerStrategy() const
{
    return d_fac_strategy;
} // getIBFACPreconditionerStrategy

/////////////////////////////// PRIVATE //////////////////////////////////////

void
StaggeredStokesIBJacobianFACPreconditioner::captureCoarseLevelDiagnostics(
    const SAMRAIVectorReal<NDIM, double>& u,
    const SAMRAIVectorReal<NDIM, double>& f,
    const int level_num,
    const int diagnostic_level_num,
    SAMRAIVectorReal<NDIM, double>* coarse_correction_out,
    SAMRAIVectorReal<NDIM, double>* coarse_rhs_out,
    bool& captured_diagnostics,
    const bool set_captured_diagnostics) const
{
    if (captured_diagnostics || level_num != diagnostic_level_num) return;
    if (coarse_rhs_out)
    {
        coarse_rhs_out->copyVector(
            Pointer<SAMRAIVectorReal<NDIM, double>>(const_cast<SAMRAIVectorReal<NDIM, double>*>(&f), false));
    }
    if (coarse_correction_out)
    {
        coarse_correction_out->copyVector(
            Pointer<SAMRAIVectorReal<NDIM, double>>(const_cast<SAMRAIVectorReal<NDIM, double>*>(&u), false));
    }
    if (set_captured_diagnostics) captured_diagnostics = true;
    return;
} // captureCoarseLevelDiagnostics

void
StaggeredStokesIBJacobianFACPreconditioner::captureSmootherLevelSolution(
    const SAMRAIVectorReal<NDIM, double>& u,
    const int level_num,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* smooth_solution_out_by_level) const
{
    if (!smooth_solution_out_by_level) return;
    if (level_num < 0) return;
    if (static_cast<std::size_t>(level_num) >= smooth_solution_out_by_level->size()) return;

    SAMRAIVectorReal<NDIM, double>* level_solution_out = (*smooth_solution_out_by_level)[static_cast<std::size_t>(level_num)];
    if (!level_solution_out) return;

    level_solution_out->copyVector(Pointer<SAMRAIVectorReal<NDIM, double>>(
        const_cast<SAMRAIVectorReal<NDIM, double>*>(&u), /*managed*/ false));
    return;
} // captureSmootherLevelSolution

void
StaggeredStokesIBJacobianFACPreconditioner::FACVCycleNoPreSmoothingWithDiagnostics(
    SAMRAIVectorReal<NDIM, double>& u,
    SAMRAIVectorReal<NDIM, double>& f,
    const int level_num,
    const int diagnostic_level_num,
    SAMRAIVectorReal<NDIM, double>* coarse_correction_out,
    SAMRAIVectorReal<NDIM, double>* coarse_rhs_out,
    bool& captured_diagnostics,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* pre_smooth_input_out_by_level,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* pre_smooth_solution_out_by_level,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* post_smooth_input_out_by_level,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* post_smooth_solution_out_by_level)
{
    if (level_num == d_coarsest_ln)
    {
        captureCoarseLevelDiagnostics(
            u, f, level_num, diagnostic_level_num, nullptr, coarse_rhs_out, captured_diagnostics, false);
        d_fac_strategy->solveCoarsestLevel(u, f, level_num);
        captureCoarseLevelDiagnostics(
            u, f, level_num, diagnostic_level_num, coarse_correction_out, nullptr, captured_diagnostics, true);
    }
    else
    {
        d_fac_strategy->restrictResidual(f, f, level_num - 1);
        FACVCycleNoPreSmoothingWithDiagnostics(
            u,
            f,
            level_num - 1,
            diagnostic_level_num,
            coarse_correction_out,
            coarse_rhs_out,
            captured_diagnostics,
            pre_smooth_input_out_by_level,
            pre_smooth_solution_out_by_level,
            post_smooth_input_out_by_level,
            post_smooth_solution_out_by_level);
        d_fac_strategy->prolongError(u, u, level_num);
        if (d_num_post_sweeps > 0)
        {
            captureSmootherLevelSolution(u, level_num, post_smooth_input_out_by_level);
            d_fac_strategy->smoothError(u, f, level_num, d_num_post_sweeps, false, true);
            captureSmootherLevelSolution(u, level_num, post_smooth_solution_out_by_level);
        }
    }
    return;
} // FACVCycleNoPreSmoothingWithDiagnostics

void
StaggeredStokesIBJacobianFACPreconditioner::muCycleWithDiagnostics(
    SAMRAIVectorReal<NDIM, double>& u,
    SAMRAIVectorReal<NDIM, double>& f,
    SAMRAIVectorReal<NDIM, double>& r,
    const int level_num,
    const int mu,
    const int diagnostic_level_num,
    SAMRAIVectorReal<NDIM, double>* coarse_correction_out,
    SAMRAIVectorReal<NDIM, double>* coarse_rhs_out,
    bool& captured_diagnostics,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* pre_smooth_input_out_by_level,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* pre_smooth_solution_out_by_level,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* post_smooth_input_out_by_level,
    const std::vector<SAMRAIVectorReal<NDIM, double>*>* post_smooth_solution_out_by_level)
{
    if (level_num == d_coarsest_ln)
    {
        captureCoarseLevelDiagnostics(
            u, f, level_num, diagnostic_level_num, nullptr, coarse_rhs_out, captured_diagnostics, false);
        d_fac_strategy->solveCoarsestLevel(u, f, level_num);
        captureCoarseLevelDiagnostics(
            u, f, level_num, diagnostic_level_num, coarse_correction_out, nullptr, captured_diagnostics, true);
    }
    else
    {
        if (d_num_pre_sweeps > 0)
        {
            captureSmootherLevelSolution(u, level_num, pre_smooth_input_out_by_level);
            d_fac_strategy->smoothError(u, f, level_num, d_num_pre_sweeps, true, false);
            captureSmootherLevelSolution(u, level_num, pre_smooth_solution_out_by_level);
        }
        d_fac_strategy->computeResidual(r, u, f, level_num - 1, level_num);
        d_fac_strategy->restrictResidual(r, f, level_num - 1);
        d_fac_strategy->setToZero(u, level_num - 1);
        for (int k = 0; k < mu; ++k)
        {
            muCycleWithDiagnostics(u,
                                   f,
                                   r,
                                   level_num - 1,
                                   mu,
                                   diagnostic_level_num,
                                   coarse_correction_out,
                                   coarse_rhs_out,
                                   captured_diagnostics,
                                   pre_smooth_input_out_by_level,
                                   pre_smooth_solution_out_by_level,
                                   post_smooth_input_out_by_level,
                                   post_smooth_solution_out_by_level);
        }
        d_fac_strategy->prolongErrorAndCorrect(u, u, level_num);
        if (d_num_post_sweeps > 0)
        {
            captureSmootherLevelSolution(u, level_num, post_smooth_input_out_by_level);
            d_fac_strategy->smoothError(u, f, level_num, d_num_post_sweeps, false, true);
            captureSmootherLevelSolution(u, level_num, post_smooth_solution_out_by_level);
        }
    }
    return;
} // muCycleWithDiagnostics

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
