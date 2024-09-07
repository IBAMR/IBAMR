// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
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

#include "ibtk/FACPreconditioner.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/SAMRAIScopedVectorCopy.h"
#include "ibtk/SAMRAIScopedVectorDuplicate.h"
#include "ibtk/ibtk_enums.h"

#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <ostream>
#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
FACPreconditioner<T>::FACPreconditioner(std::string object_name,
                                        Pointer<FACPreconditionerStrategy<T> > fac_strategy,
                                        tbox::Pointer<tbox::Database> input_db,
                                        const std::string& /*default_options_prefix*/)
    : d_fac_strategy(fac_strategy)
{
    // Setup default options.
    GeneralSolver::init(std::move(object_name), /*homogeneous_bc*/ true);
    d_initial_guess_nonzero = false;
    d_rel_residual_tol = 1.0e-5;
    d_abs_residual_tol = 1.0e-50;
    d_max_iterations = 1;

    // Register this class with the FACPreconditionerStrategy object.
    d_fac_strategy->setFACPreconditioner(Pointer<FACPreconditioner>(this, false));

    // Initialize object with data read from input database.
    if (input_db)
    {
        getFromInput(input_db);
    }
    return;
} // FACPreconditioner

template <class T>
FACPreconditioner<T>::~FACPreconditioner()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~FACPreconditioner

template <class T>
void
FACPreconditioner<T>::setHomogeneousBc(const bool homogeneous_bc)
{
    LinearSolver::setHomogeneousBc(homogeneous_bc);
    d_fac_strategy->setHomogeneousBc(homogeneous_bc);
    return;
} // setHomogeneousBc

template <class T>
void
FACPreconditioner<T>::setSolutionTime(const double solution_time)
{
    LinearSolver::setSolutionTime(solution_time);
    d_fac_strategy->setSolutionTime(solution_time);
    return;
} // setSolutionTime

template <class T>
void
FACPreconditioner<T>::setTimeInterval(const double current_time, const double new_time)
{
    LinearSolver::setTimeInterval(current_time, new_time);
    d_fac_strategy->setTimeInterval(current_time, new_time);
    return;
} // setTimeInterval

template <>
void
FACPreconditioner<double>::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& solution,
                                                 const SAMRAIVectorReal<NDIM, double>& rhs)
{
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        deallocateSolverState();
    }

    // Setup operator state.
    d_hierarchy = solution.getPatchHierarchy();
    d_coarsest_ln = solution.getCoarsestLevelNumber();
    d_finest_ln = solution.getFinestLevelNumber();

#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == rhs.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == rhs.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == rhs.getFinestLevelNumber());
#endif
    d_fac_strategy->initializeOperatorState(solution, rhs);

    // Allocate scratch data.
    d_fac_strategy->allocateScratchData();

    // Indicate the operator is initialized.
    d_is_initialized = true;
    return;
} // initializeSolverState

template <class T>
void
FACPreconditioner<T>::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& solution_full_precision,
                                            const SAMRAIVectorReal<NDIM, double>& rhs_full_precision)
{
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        deallocateSolverState();
    }

    // Setup operator state.
    d_hierarchy = solution_full_precision.getPatchHierarchy();
    d_coarsest_ln = solution_full_precision.getCoarsestLevelNumber();
    d_finest_ln = solution_full_precision.getFinestLevelNumber();

#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == rhs_full_precision.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == rhs_full_precision.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == rhs_full_precision.getFinestLevelNumber());
#endif
    SAMRAIScopedVectorDuplicate<T> solution_reduced_precision(solution_full_precision);
    SAMRAIScopedVectorDuplicate<T> rhs_reduced_precision(rhs_full_precision);
    d_fac_strategy->initializeOperatorState(solution_reduced_precision, rhs_reduced_precision);

    // Allocate scratch data.
    d_fac_strategy->allocateScratchData();

    // Indicate the operator is initialized.
    d_is_initialized = true;
    return;
} // initializeSolverState

template <>
bool
FACPreconditioner<double>::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    // Set the initial guess to equal zero.
    x.setToScalar(0.0, /*interior_only*/ false);

    // Apply a single FAC cycle.
    if (d_cycle_type == V_CYCLE && d_num_pre_sweeps == 0)
    {
        // V-cycle MG without presmoothing keeps the residual equal to the
        // initial right-hand-side vector f, so we can simply use that vector
        // for the residual in the FAC algorithm.
        FACVCycleNoPreSmoothing(x, b, d_finest_ln);
    }
    else
    {
        // Clone the right-hand-side vector to avoid modifying it during the
        // preconditioning operation.
        SAMRAIScopedVectorCopy<double> f(b);
        SAMRAIScopedVectorDuplicate<double> r(b);

        switch (d_cycle_type)
        {
        case F_CYCLE:
            FCycle(x, f, r, d_finest_ln);
            break;
        case FMG_CYCLE:
            FMGCycle(x, f, r, d_finest_ln, 1);
            break;
        case V_CYCLE:
            muCycle(x, f, r, d_finest_ln, 1);
            break;
        case W_CYCLE:
            muCycle(x, f, r, d_finest_ln, 2);
            break;
        default:
            TBOX_ERROR(d_object_name << "::solveSystem():\n"
                                     << "  unsupported FAC cycle type: " << enum_to_string<MGCycleType>(d_cycle_type)
                                     << "." << std::endl);
        }
    }

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();
    return true;
} // solveSystem

template <class T>
bool
FACPreconditioner<T>::solveSystem(SAMRAIVectorReal<NDIM, double>& x_full_precision,
                                  SAMRAIVectorReal<NDIM, double>& b_full_precision)
{
    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x_full_precision, b_full_precision);

    // Set the initial guess to equal zero.
    x_full_precision.setToScalar(0.0, /*interior_only*/ false);

    // Clone the right-hand-side vector to avoid modifying it during the
    // preconditioning operation.
    SAMRAIScopedVectorCopy<T> x_reduced_precision(x_full_precision);
    SAMRAIScopedVectorCopy<T> f_reduced_precision(b_full_precision);
    SAMRAIScopedVectorDuplicate<T> r_reduced_precision(b_full_precision);

    // Apply a single FAC cycle.
    if (d_cycle_type == V_CYCLE && d_num_pre_sweeps == 0)
    {
        // V-cycle MG without presmoothing keeps the residual equal to the
        // initial right-hand-side vector f, so we can simply use that vector
        // for the residual in the FAC algorithm.
        FACVCycleNoPreSmoothing(x_reduced_precision, f_reduced_precision, d_finest_ln);
    }
    else
    {
        switch (d_cycle_type)
        {
        case F_CYCLE:
            FCycle(x_reduced_precision, f_reduced_precision, r_reduced_precision, d_finest_ln);
            break;
        case FMG_CYCLE:
            FMGCycle(x_reduced_precision, f_reduced_precision, r_reduced_precision, d_finest_ln, 1);
            break;
        case V_CYCLE:
            muCycle(x_reduced_precision, f_reduced_precision, r_reduced_precision, d_finest_ln, 1);
            break;
        case W_CYCLE:
            muCycle(x_reduced_precision, f_reduced_precision, r_reduced_precision, d_finest_ln, 2);
            break;
        default:
            TBOX_ERROR(d_object_name << "::solveSystem():\n"
                                     << "  unsupported FAC cycle type: " << enum_to_string<MGCycleType>(d_cycle_type)
                                     << "." << std::endl);
        }
    }

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    // Restore the solution vector.
    x_reduced_precision.transformToVector(x_full_precision);
    return true;
} // solveSystem

template <class T>
void
FACPreconditioner<T>::deallocateSolverState()
{
    if (!d_is_initialized) return;

    // Deallocate scratch data.
    d_fac_strategy->deallocateScratchData();

    // Deallocate operator state.
    d_fac_strategy->deallocateOperatorState();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;
    return;
} // deallocateSolverState

template <class T>
void
FACPreconditioner<T>::setInitialGuessNonzero(bool initial_guess_nonzero)
{
    if (initial_guess_nonzero)
    {
        TBOX_ERROR(d_object_name << "::setInitialGuessNonzero()\n"
                                 << "  class IBTK::FACPreconditioner requires a zero initial guess" << std::endl);
    }
    return;
} // setInitialGuessNonzero

template <class T>
void
FACPreconditioner<T>::setMaxIterations(int max_iterations)
{
    if (max_iterations != 1)
    {
        TBOX_ERROR(d_object_name << "::setMaxIterations()\n"
                                 << "  class IBTK::FACPreconditioner only performs a single iteration" << std::endl);
    }
    return;
} // setMaxIterations

template <class T>
void
FACPreconditioner<T>::setMGCycleType(MGCycleType cycle_type)
{
    d_cycle_type = cycle_type;
    return;
} // setMGCycleType

template <class T>
MGCycleType
FACPreconditioner<T>::getMGCycleType() const
{
    return d_cycle_type;
} // getMGCycleType

template <class T>
void
FACPreconditioner<T>::setNumPreSmoothingSweeps(int num_pre_sweeps)
{
    d_num_pre_sweeps = num_pre_sweeps;
    return;
} // setNumPreSmoothingSweeps

template <class T>
int
FACPreconditioner<T>::getNumPreSmoothingSweeps() const
{
    return d_num_pre_sweeps;
} // getNumPreSmoothingSweeps

template <class T>
void
FACPreconditioner<T>::setNumPostSmoothingSweeps(int num_post_sweeps)
{
    d_num_post_sweeps = num_post_sweeps;
    return;
} // setNumPostSmoothingSweeps

template <class T>
int
FACPreconditioner<T>::getNumPostSmoothingSweeps() const
{
    return d_num_post_sweeps;
} // getNumPostSmoothingSweeps

template <class T>
Pointer<FACPreconditionerStrategy<T> >
FACPreconditioner<T>::getFACPreconditionerStrategy() const
{
    return d_fac_strategy;
} // getFACPreconditionerStrategy

/////////////////////////////// PROTECTED ////////////////////////////////////

template <class T>
void
FACPreconditioner<T>::FACVCycleNoPreSmoothing(SAMRAIVectorReal<NDIM, T>& u, SAMRAIVectorReal<NDIM, T>& f, int level_num)
{
    if (level_num == d_coarsest_ln)
    {
        // Solve Au = f on the coarsest level.
        d_fac_strategy->solveCoarsestLevel(u, f, level_num);
    }
    else
    {
        // Restrict the residual to the next coarser level.
        d_fac_strategy->restrictResidual(f, f, level_num - 1);

        // Recursively call the FAC algorithm.
        FACVCycleNoPreSmoothing(u, f, level_num - 1);

        // Prolong the error from the next coarser level.  Because we did not
        // perform any presmoothing, we do not need to correct the solution on
        // the current level.
        d_fac_strategy->prolongError(u, u, level_num);

        // Smooth error on the current level.
        if (d_num_post_sweeps > 0)
        {
            d_fac_strategy->smoothError(u, f, level_num, d_num_post_sweeps, false, true);
        }
    }
    return;
} // FACVCycleNoPreSmoothing

template <class T>
void
FACPreconditioner<T>::muCycle(SAMRAIVectorReal<NDIM, T>& u,
                              SAMRAIVectorReal<NDIM, T>& f,
                              SAMRAIVectorReal<NDIM, T>& r,
                              int level_num,
                              int mu)
{
    if (level_num == d_coarsest_ln)
    {
        d_fac_strategy->solveCoarsestLevel(u, f, level_num);
    }
    else
    {
        if (d_num_pre_sweeps > 0)
        {
            d_fac_strategy->smoothError(u, f, level_num, d_num_pre_sweeps, true, false);
        }
        d_fac_strategy->computeResidual(r, u, f, level_num - 1, level_num);
        d_fac_strategy->restrictResidual(r, f, level_num - 1);
        d_fac_strategy->setToZero(u, level_num - 1);
        for (int k = 0; k < mu; ++k) muCycle(u, f, r, level_num - 1, mu);
        d_fac_strategy->prolongErrorAndCorrect(u, u, level_num);
        if (d_num_post_sweeps > 0)
        {
            d_fac_strategy->smoothError(u, f, level_num, d_num_post_sweeps, false, true);
        }
    }
    return;
} // muCycle

template <class T>
void
FACPreconditioner<T>::FCycle(SAMRAIVectorReal<NDIM, T>& u,
                             SAMRAIVectorReal<NDIM, T>& f,
                             SAMRAIVectorReal<NDIM, T>& r,
                             int level_num)
{
    if (level_num == d_coarsest_ln)
    {
        d_fac_strategy->solveCoarsestLevel(u, f, level_num);
    }
    else
    {
        if (d_num_pre_sweeps > 0)
        {
            d_fac_strategy->smoothError(u, f, level_num, d_num_pre_sweeps, true, false);
        }
        d_fac_strategy->computeResidual(r, u, f, level_num - 1, level_num);
        d_fac_strategy->restrictResidual(r, f, level_num - 1);
        d_fac_strategy->setToZero(u, level_num - 1);
        muCycle(u, f, r, level_num - 1, 2);
        muCycle(u, f, r, level_num - 1, 1);
        d_fac_strategy->prolongErrorAndCorrect(u, u, level_num);
        if (d_num_post_sweeps > 0)
        {
            d_fac_strategy->smoothError(u, f, level_num, d_num_post_sweeps, false, true);
        }
    }
    return;
} // FCycle

template <class T>
void
FACPreconditioner<T>::FMGCycle(SAMRAIVectorReal<NDIM, T>& u,
                               SAMRAIVectorReal<NDIM, T>& f,
                               SAMRAIVectorReal<NDIM, T>& r,
                               int level_num,
                               int mu)
{
    if (level_num == d_coarsest_ln)
    {
        d_fac_strategy->setToZero(u, level_num);
    }
    else
    {
        d_fac_strategy->restrictResidual(f, f, level_num - 1);
        FMGCycle(u, f, r, level_num - 1, mu);
        d_fac_strategy->prolongErrorAndCorrect(u, u, level_num);
    }
    muCycle(u, f, r, level_num, mu);
    return;
} // FMGCycle

/////////////////////////////// PRIVATE //////////////////////////////////////

template <class T>
void
FACPreconditioner<T>::getFromInput(tbox::Pointer<tbox::Database> db)
{
    if (!db) return;
    if (db->keyExists("cycle_type")) setMGCycleType(string_to_enum<MGCycleType>(db->getString("cycle_type")));
    if (db->keyExists("num_pre_sweeps")) setNumPreSmoothingSweeps(db->getInteger("num_pre_sweeps"));
    if (db->keyExists("num_post_sweeps")) setNumPostSmoothingSweeps(db->getInteger("num_post_sweeps"));
    if (db->keyExists("enable_logging")) setLoggingEnabled(db->getBool("enable_logging"));
    return;
} // getFromInput

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template class FACPreconditioner<float>;
template class FACPreconditioner<double>;

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
