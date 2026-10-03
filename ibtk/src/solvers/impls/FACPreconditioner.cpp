// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <ibtk/FACPreconditioner.h>
#include <ibtk/FACPreconditionerStrategy.h>
#include <ibtk/LinearSolver.h>
#include <ibtk/SAMRAIScopedVectorCopy.h>
#include <ibtk/SAMRAIScopedVectorDuplicate.h>
#include <ibtk/ibtk_enums.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <MultiblockDataTranslator.h>
#include <PatchHierarchy.h>
#include <SAMRAIVectorReal.h>

#include <ostream>
#include <string>
#include <utility>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Return a vector on levels coarsest_ln through finest_ln that refers to the patch data of vector.
Pointer<SAMRAIVectorReal<NDIM, double>>
get_range_vector(const SAMRAIVectorReal<NDIM, double>& vector, const int coarsest_ln, const int finest_ln)
{
    Pointer<SAMRAIVectorReal<NDIM, double>> view = new SAMRAIVectorReal<NDIM, double>(
        vector.getName() + "::range", vector.getPatchHierarchy(), coarsest_ln, finest_ln);
    for (int comp = 0; comp < vector.getNumberOfComponents(); ++comp)
    {
        view->addComponent(vector.getComponentVariable(comp),
                           vector.getComponentDescriptorIndex(comp),
                           vector.getControlVolumeIndex(comp));
    }
    return view;
} // get_range_vector

// Allocate a zero vector with the components of vector on levels up to finest_ln.
Pointer<SAMRAIVectorReal<NDIM, double>>
allocate_range_vector(const SAMRAIVectorReal<NDIM, double>& vector, const std::string& name, const int finest_ln)
{
    Pointer<SAMRAIVectorReal<NDIM, double>> scratch =
        get_range_vector(vector, vector.getCoarsestLevelNumber(), finest_ln)->cloneVector(name);
    scratch->allocateVectorData();
    scratch->setToScalar(0.0, /*interior_only*/ false);
    return scratch;
} // allocate_range_vector
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

FACPreconditioner::FACPreconditioner(std::string object_name,
                                     Pointer<FACPreconditionerStrategy> fac_strategy,
                                     tbox::Pointer<tbox::Database> input_db,
                                     const std::string& /*default_options_prefix*/)
    : d_fac_strategy(fac_strategy)
{
    // Setup default options. FACPreconditioner implements only the FAC correction
    // scheme, which requires homogeneous boundary conditions; see setHomogeneousBc().
    GeneralSolver::init(std::move(object_name), FACPreconditionerStrategy::ALWAYS_HOMOGENEOUS_BC);
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

FACPreconditioner::~FACPreconditioner()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~FACPreconditioner

void
FACPreconditioner::setHomogeneousBc(const bool homogeneous_bc)
{
    if (homogeneous_bc != FACPreconditionerStrategy::ALWAYS_HOMOGENEOUS_BC)
    {
        TBOX_ERROR(d_object_name << "::setHomogeneousBc():\n"
                                 << "  FACPreconditioner implements only the FAC correction scheme and\n"
                                 << "  so requires homogeneous boundary conditions." << std::endl);
    }
    LinearSolver::setHomogeneousBc(homogeneous_bc);
    return;
} // setHomogeneousBc

void
FACPreconditioner::setSolutionTime(const double solution_time)
{
    LinearSolver::setSolutionTime(solution_time);
    d_fac_strategy->setSolutionTime(solution_time);
    return;
} // setSolutionTime

void
FACPreconditioner::setTimeInterval(const double current_time, const double new_time)
{
    LinearSolver::setTimeInterval(current_time, new_time);
    d_fac_strategy->setTimeInterval(current_time, new_time);
    return;
} // setTimeInterval

bool
FACPreconditioner::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

#if !defined(NDEBUG)
    TBOX_ASSERT(x.getPatchHierarchy() == d_hierarchy && b.getPatchHierarchy() == d_hierarchy);
    TBOX_ASSERT(x.getCoarsestLevelNumber() == d_coarsest_ln && b.getCoarsestLevelNumber() == d_coarsest_ln);
    TBOX_ASSERT(x.getFinestLevelNumber() == d_finest_ln && b.getFinestLevelNumber() == d_finest_ln);
#endif

    // Allocate the scratch data that the current cycle options require.
    allocateCycleScratchData(x, b);

    // Set the initial guess to equal zero.
    x.setToScalar(0.0, /*interior_only*/ false);

    // Apply a single FAC cycle.
    if (d_cycle_type == V_CYCLE || d_cycle_type == W_CYCLE)
    {
        zeroStartCycle(x, b, d_finest_ln, d_cycle_type);
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

void
FACPreconditioner::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& solution,
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
    if (d_coarsest_ln > 0 && d_finest_ln > d_coarsest_ln)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                 << "  vectors that span more than one level must start at level zero." << std::endl);
    }
    d_fac_strategy->initializeOperatorState(solution, rhs);

    // Allocate scratch data.
    d_fac_strategy->allocateScratchData();
    d_residual_vectors.resize(d_finest_ln + 1);
    d_rhs_vectors.resize(d_finest_ln + 1);
    d_correction_vectors.resize(d_finest_ln + 1);

    // Indicate the operator is initialized.
    d_is_initialized = true;
    return;
} // initializeSolverState

void
FACPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    // Free the patch data and their descriptor indices.
    for (auto* vectors : { &d_residual_vectors, &d_rhs_vectors, &d_correction_vectors })
    {
        for (auto& vector : *vectors)
        {
            if (vector)
            {
                free_vector_components(*vector);
            }
        }
        vectors->clear();
    }

    // Deallocate scratch data.
    d_fac_strategy->deallocateScratchData();

    // Deallocate operator state.
    d_fac_strategy->deallocateOperatorState();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;
    return;
} // deallocateSolverState

void
FACPreconditioner::setInitialGuessNonzero(bool initial_guess_nonzero)
{
    if (initial_guess_nonzero)
    {
        TBOX_ERROR(d_object_name << "::setInitialGuessNonzero()\n"
                                 << "  class IBTK::FACPreconditioner requires a zero initial guess" << std::endl);
    }
    return;
} // setInitialGuessNonzero

void
FACPreconditioner::setMaxIterations(int max_iterations)
{
    if (max_iterations != 1)
    {
        TBOX_ERROR(d_object_name << "::setMaxIterations()\n"
                                 << "  class IBTK::FACPreconditioner only performs a single iteration" << std::endl);
    }
    return;
} // setMaxIterations

void
FACPreconditioner::setMGCycleType(MGCycleType cycle_type)
{
    d_cycle_type = cycle_type;
    return;
} // setMGCycleType

MGCycleType
FACPreconditioner::getMGCycleType() const
{
    return d_cycle_type;
} // getMGCycleType

void
FACPreconditioner::setNumPreSmoothingSweeps(int num_pre_sweeps)
{
    d_num_pre_sweeps = num_pre_sweeps;
    return;
} // setNumPreSmoothingSweeps

int
FACPreconditioner::getNumPreSmoothingSweeps() const
{
    return d_num_pre_sweeps;
} // getNumPreSmoothingSweeps

void
FACPreconditioner::setNumPostSmoothingSweeps(int num_post_sweeps)
{
    d_num_post_sweeps = num_post_sweeps;
    return;
} // setNumPostSmoothingSweeps

int
FACPreconditioner::getNumPostSmoothingSweeps() const
{
    return d_num_post_sweeps;
} // getNumPostSmoothingSweeps

Pointer<FACPreconditionerStrategy>
FACPreconditioner::getFACPreconditionerStrategy() const
{
    return d_fac_strategy;
} // getFACPreconditionerStrategy

/////////////////////////////// PROTECTED ////////////////////////////////////

void
FACPreconditioner::muCycle(SAMRAIVectorReal<NDIM, double>& u,
                           SAMRAIVectorReal<NDIM, double>& f,
                           SAMRAIVectorReal<NDIM, double>& r,
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

void
FACPreconditioner::FCycle(SAMRAIVectorReal<NDIM, double>& u,
                          SAMRAIVectorReal<NDIM, double>& f,
                          SAMRAIVectorReal<NDIM, double>& r,
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

void
FACPreconditioner::FMGCycle(SAMRAIVectorReal<NDIM, double>& u,
                            SAMRAIVectorReal<NDIM, double>& f,
                            SAMRAIVectorReal<NDIM, double>& r,
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

void
FACPreconditioner::zeroStartCycle(SAMRAIVectorReal<NDIM, double>& u,
                                  SAMRAIVectorReal<NDIM, double>& f,
                                  const int level_num,
                                  const MGCycleType cycle_type)
{
    if (level_num == d_coarsest_ln)
    {
        d_fac_strategy->solveCoarsestLevel(u, f, level_num);
    }
    else
    {
        // Without presmoothing the residual is f.
        SAMRAIVectorReal<NDIM, double>& residual = d_num_pre_sweeps > 0 ? *d_residual_vectors[level_num] : f;
        if (d_num_pre_sweeps > 0)
        {
            d_fac_strategy->smoothError(u, f, level_num, d_num_pre_sweeps, true, false);

            // u is nonzero only on level_num, but the composite-grid residual
            // differs from f on the coarser levels as well.
            d_fac_strategy->computeResidual(residual, u, f, d_coarsest_ln, level_num);
            d_fac_strategy->restrictResidual(residual, residual, level_num - 1);

            // computeResidual() overwrote u on covered coarse cells; zero the
            // coarser levels before using them for the coarse correction.
            for (int ln = d_coarsest_ln; ln < level_num; ++ln)
            {
                d_fac_strategy->setToZero(u, ln);
            }
        }
        else
        {
            d_fac_strategy->restrictResidual(f, f, level_num - 1);
        }

        zeroStartCycle(u, residual, level_num - 1, cycle_type);

        // A W-cycle visits the coarser levels a second time with another W-cycle.
        if (cycle_type == W_CYCLE)
        {
            improveCycle(u, residual, level_num - 1, W_CYCLE);
        }

        if (d_num_pre_sweeps > 0)
        {
            // Presmoothing made u nonzero on level_num; add the prolongation of the
            // coarse correction to it.
            d_fac_strategy->prolongErrorAndCorrect(u, u, level_num);
        }
        else
        {
            // u is zero on level_num, including its ghost values, so the prolongation
            // of the coarse correction is the correction.
            d_fac_strategy->prolongError(u, u, level_num);
        }
        if (d_num_post_sweeps > 0)
        {
            d_fac_strategy->smoothError(u, f, level_num, d_num_post_sweeps, false, true);
        }
    }
    return;
} // zeroStartCycle

void
FACPreconditioner::improveCycle(SAMRAIVectorReal<NDIM, double>& u,
                                const SAMRAIVectorReal<NDIM, double>& f,
                                const int level_num,
                                const MGCycleType cycle_type)
{
    Pointer<SAMRAIVectorReal<NDIM, double>> solution = get_range_vector(u, d_coarsest_ln, level_num);
    Pointer<SAMRAIVectorReal<NDIM, double>> rhs = d_rhs_vectors[level_num];
    Pointer<SAMRAIVectorReal<NDIM, double>> correction = d_correction_vectors[level_num];
    d_fac_strategy->computeResidual(*rhs, *solution, f, d_coarsest_ln, level_num);
    correction->setToScalar(0.0, /*interior_only*/ false);
    zeroStartCycle(*correction, *rhs, level_num, cycle_type);
    solution->add(solution, correction, /*interior_only*/ false);
    return;
} // improveCycle

void
FACPreconditioner::allocateCycleScratchData(const SAMRAIVectorReal<NDIM, double>& solution,
                                            const SAMRAIVectorReal<NDIM, double>& rhs)
{
    // Only V- and W-cycles on more than one level use scratch data.
    if (d_coarsest_ln == d_finest_ln || d_cycle_type == F_CYCLE || d_cycle_type == FMG_CYCLE)
    {
        return;
    }
    if (d_num_pre_sweeps > 0)
    {
        for (int ln = d_coarsest_ln + 1; ln <= d_finest_ln; ++ln)
        {
            if (!d_residual_vectors[ln])
            {
                d_residual_vectors[ln] =
                    allocate_range_vector(rhs, d_object_name + "::residual::level_" + std::to_string(ln), ln);
            }
        }
    }
    if (d_cycle_type == W_CYCLE)
    {
        for (int ln = d_coarsest_ln; ln < d_finest_ln; ++ln)
        {
            const std::string suffix = "::level_" + std::to_string(ln);
            if (!d_rhs_vectors[ln])
            {
                d_rhs_vectors[ln] = allocate_range_vector(rhs, d_object_name + "::rhs" + suffix, ln);
            }
            if (!d_correction_vectors[ln])
            {
                d_correction_vectors[ln] = allocate_range_vector(solution, d_object_name + "::correction" + suffix, ln);
            }
        }
    }
    return;
} // allocateCycleScratchData

void
FACPreconditioner::getFromInput(tbox::Pointer<tbox::Database> db)
{
    if (!db) return;
    if (db->keyExists("cycle_type")) setMGCycleType(string_to_enum<MGCycleType>(db->getString("cycle_type")));
    if (db->keyExists("num_pre_sweeps")) setNumPreSmoothingSweeps(db->getInteger("num_pre_sweeps"));
    if (db->keyExists("num_post_sweeps")) setNumPostSmoothingSweeps(db->getInteger("num_post_sweeps"));
    if (db->keyExists("enable_logging")) setLoggingEnabled(db->getBool("enable_logging"));
    return;
} // getFromInput

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
