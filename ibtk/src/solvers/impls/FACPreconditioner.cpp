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
#include <ibtk/ibtk_enums.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <HierarchyDataOpsManager.h>
#include <PatchHierarchy.h>
#include <SAMRAIVectorReal.h>
#include <Variable.h>
#include <VariableDatabase.h>

#include <ostream>
#include <string>
#include <utility>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

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
    // Parameters may have changed since initialization. Previously allocated
    // workspaces remain reusable; allocate only any newly required storage.
    allocateCycleScratchData(x, b);
    x.setToScalar(0.0, /*interior_only*/ false);

    if (d_cycle_type == V_CYCLE && d_num_pre_sweeps == 0)
    {
        FACVCycleNoPreSmoothing(x, b, d_finest_ln);
    }
    else if (d_cycle_type == FMG_CYCLE)
    {
        FMGCycle(x, b, d_finest_ln);
    }
    else
    {
        zeroStartCycle(x, b, d_finest_ln, d_cycle_type);
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
    d_fac_strategy->initializeOperatorState(solution, rhs);

    // Allocate scratch data.
    d_fac_strategy->allocateScratchData();
    d_residual_vectors.resize(d_finest_ln + 1);
    d_rhs_vectors.resize(d_finest_ln + 1);
    d_correction_vectors.resize(d_finest_ln + 1);
    d_fmg_rhs_vectors.resize(d_finest_ln + 1);
    allocateCycleScratchData(solution, rhs);

    // Indicate the operator is initialized.
    d_is_initialized = true;
    return;
} // initializeSolverState

void
FACPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    // Free patch indices as well as data, including after hierarchy changes.
    for (std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>>* vectors :
         { &d_residual_vectors, &d_rhs_vectors, &d_correction_vectors, &d_fmg_rhs_vectors })
    {
        for (Pointer<SAMRAIVectorReal<NDIM, double>>& vector : *vectors)
        {
            if (vector)
            {
                free_vector_components(*vector);
            }
        }
        vectors->clear();
    }
    if (d_evaluation_vector)
    {
        free_vector_components(*d_evaluation_vector);
    }
    d_evaluation_vector.setNull();
    d_vector_data_ops.clear();
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
    switch (cycle_type)
    {
    case V_CYCLE:
    case W_CYCLE:
    case F_CYCLE:
    case FMG_CYCLE:
        break;
    default:
        TBOX_ERROR(d_object_name << "::setMGCycleType(): unsupported FAC cycle type." << std::endl);
    }
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
FACPreconditioner::FACVCycleNoPreSmoothing(SAMRAIVectorReal<NDIM, double>& u,
                                           SAMRAIVectorReal<NDIM, double>& f,
                                           int level_num)
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
        Pointer<SAMRAIVectorReal<NDIM, double>> residual = d_residual_vectors[level_num];
        if (d_num_pre_sweeps > 0)
        {
            d_fac_strategy->smoothError(u, f, level_num, d_num_pre_sweeps, true, false);

            // Only the top level has been smoothed. This invocation owns u, so
            // residual evaluation can synchronize its covered coarse values.
            // The child equation needs the residual over the whole parent range.
            d_fac_strategy->computeResidual(*residual, u, f, d_coarsest_ln, level_num);
            d_fac_strategy->restrictResidual(*residual, *residual, level_num - 1);
        }
        else
        {
            // A(0) = 0. Preserve the borrowed RHS: restriction itself copies
            // the child top level, so only its lower prefix needs an explicit copy.
            if (level_num > d_coarsest_ln + 1)
            {
                Pointer<SAMRAIVectorReal<NDIM, double>> prefix =
                    getRangeVector(*residual, d_coarsest_ln, level_num - 2);
                prefix->copyVector(getRangeVector(f, d_coarsest_ln, level_num - 2), /*interior_only*/ false);
            }
            d_fac_strategy->restrictResidual(f, *residual, level_num - 1);
        }
        Pointer<SAMRAIVectorReal<NDIM, double>> child_u = getRangeVector(u, d_coarsest_ln, level_num - 1);
        Pointer<SAMRAIVectorReal<NDIM, double>> child_rhs = getRangeVector(*residual, d_coarsest_ln, level_num - 1);
        if (d_num_pre_sweeps > 0)
        {
            // Remove every covered value written by residual synchronization,
            // including ghosts, before reusing this span as a child correction.
            child_u->setToScalar(0.0, /*interior_only*/ false);
        }

        // The first child starts at zero. Each additional visit improves the
        // accumulated child correction against the same, independently owned RHS.
        if (cycle_type == F_CYCLE)
        {
            zeroStartCycle(*child_u, *child_rhs, level_num - 1, F_CYCLE);
            improveCycle(*child_u, *child_rhs, level_num - 1, V_CYCLE);
        }
        else
        {
            int multiplicity = 1;
            switch (cycle_type)
            {
            case V_CYCLE:
                break;
            case W_CYCLE:
                multiplicity = 2;
                break;
            default:
                TBOX_ERROR(d_object_name << "::zeroStartCycle(): unsupported FAC cycle type." << std::endl);
            }
            zeroStartCycle(*child_u, *child_rhs, level_num - 1, cycle_type);
            for (int visit = 1; visit < multiplicity; ++visit)
            {
                improveCycle(*child_u, *child_rhs, level_num - 1, cycle_type);
            }
        }

        // The child already occupies u's complete lower span. The same-vector
        // transfer adds only the fine correction and refreshes its coarse-fine seeds.
        d_fac_strategy->prolongErrorAndCorrect(u, u, level_num);
        if (d_num_post_sweeps > 0)
        {
            d_fac_strategy->smoothError(u, f, level_num, d_num_post_sweeps, false, true);
        }
    }
    d_fac_strategy->fillGhostCellsNoCoarse(u, level_num);
    return;
} // zeroStartCycle

void
FACPreconditioner::improveCycle(SAMRAIVectorReal<NDIM, double>& u,
                                SAMRAIVectorReal<NDIM, double>& f,
                                const int level_num,
                                const MGCycleType cycle_type)
{
    Pointer<SAMRAIVectorReal<NDIM, double>> evaluation = getRangeVector(*d_evaluation_vector, d_coarsest_ln, level_num);
    evaluation->copyVector(getRangeVector(u, d_coarsest_ln, level_num), /*interior_only*/ false);
    Pointer<SAMRAIVectorReal<NDIM, double>> rhs = d_rhs_vectors[level_num];
    Pointer<SAMRAIVectorReal<NDIM, double>> correction = d_correction_vectors[level_num];
    d_fac_strategy->computeResidual(*rhs, *evaluation, f, d_coarsest_ln, level_num);
    correction->setToScalar(0.0, /*interior_only*/ false);
    zeroStartCycle(*correction, *rhs, level_num, cycle_type);
    Pointer<SAMRAIVectorReal<NDIM, double>> solution = getRangeVector(u, d_coarsest_ln, level_num);
    solution->add(solution, correction, /*interior_only*/ false);
    return;
} // improveCycle

void
FACPreconditioner::FMGCycle(SAMRAIVectorReal<NDIM, double>& u, SAMRAIVectorReal<NDIM, double>& f, const int level_num)
{
    if (level_num == d_coarsest_ln)
    {
        zeroStartCycle(u, f, level_num, V_CYCLE);
    }
    else
    {
        // Keep the parent RHS intact while constructing the next nested equation.
        Pointer<SAMRAIVectorReal<NDIM, double>> child_rhs = d_fmg_rhs_vectors[level_num - 1];
        if (level_num > d_coarsest_ln + 1)
        {
            Pointer<SAMRAIVectorReal<NDIM, double>> prefix = getRangeVector(*child_rhs, d_coarsest_ln, level_num - 2);
            prefix->copyVector(getRangeVector(f, d_coarsest_ln, level_num - 2), /*interior_only*/ false);
        }
        d_fac_strategy->restrictResidual(f, *child_rhs, level_num - 1);
        Pointer<SAMRAIVectorReal<NDIM, double>> child_u = getRangeVector(u, d_coarsest_ln, level_num - 1);
        FMGCycle(*child_u, *child_rhs, level_num - 1);
        d_fac_strategy->prolongError(u, u, level_num);
        improveCycle(u, f, level_num, V_CYCLE);
    }
    return;
} // FMGCycle

/////////////////////////////// PRIVATE //////////////////////////////////////

Pointer<SAMRAIVectorReal<NDIM, double>>
FACPreconditioner::getRangeVector(const SAMRAIVectorReal<NDIM, double>& vector,
                                  const int coarsest_ln,
                                  const int finest_ln) const
{
    Pointer<SAMRAIVectorReal<NDIM, double>> view = new SAMRAIVectorReal<NDIM, double>(
        vector.getName() + "::range", vector.getPatchHierarchy(), coarsest_ln, finest_ln);
    for (int comp = 0; comp < vector.getNumberOfComponents(); ++comp)
    {
        view->addComponent(vector.getComponentVariable(comp),
                           vector.getComponentDescriptorIndex(comp),
                           vector.getControlVolumeIndex(comp),
                           d_vector_data_ops[comp]);
    }
    return view;
} // getRangeVector

Pointer<SAMRAIVectorReal<NDIM, double>>
FACPreconditioner::allocateRangeVector(const SAMRAIVectorReal<NDIM, double>& vector,
                                       const std::string& name,
                                       const int finest_ln) const
{
    Pointer<SAMRAIVectorReal<NDIM, double>> scratch = new SAMRAIVectorReal<NDIM, double>(
        name, vector.getPatchHierarchy(), vector.getCoarsestLevelNumber(), finest_ln);
    VariableDatabase<NDIM>* variable_db = VariableDatabase<NDIM>::getDatabase();
    for (int comp = 0; comp < vector.getNumberOfComponents(); ++comp)
    {
        const Pointer<Variable<NDIM>> variable = vector.getComponentVariable(comp);
        const int index = variable_db->registerClonedPatchDataIndex(variable, vector.getComponentDescriptorIndex(comp));
        scratch->addComponent(variable, index, vector.getControlVolumeIndex(comp), d_vector_data_ops[comp]);
    }
    scratch->allocateVectorData();
    scratch->setToScalar(0.0, /*interior_only*/ false);
    return scratch;
} // allocateRangeVector

void
FACPreconditioner::allocateCycleScratchData(const SAMRAIVectorReal<NDIM, double>& solution,
                                            const SAMRAIVectorReal<NDIM, double>& rhs)
{
    // Preserve the inexpensive default V-cycle path, including its existing
    // in-place restriction of covered RHS data.
    if (d_cycle_type == V_CYCLE && d_num_pre_sweeps == 0)
    {
        return;
    }
    if (d_coarsest_ln != 0)
    {
        TBOX_ERROR(d_object_name << "::allocateCycleScratchData():\n"
                                 << "  this FAC cycle requires coarsest level zero: the current strategy residual\n"
                                 << "  ghost-fill operators can access data below a nonzero coarsest level."
                                 << std::endl);
    }
    if (d_vector_data_ops.empty())
    {
        // cloneVector() shares operation objects with its source. Restricting
        // their level range would also restrict a caller's subsequent direct
        // hierarchy operations, so neither scratch nor views may borrow them.
        HierarchyDataOpsManager<NDIM>* manager = HierarchyDataOpsManager<NDIM>::getManager();
        for (int comp = 0; comp < solution.getNumberOfComponents(); ++comp)
        {
            d_vector_data_ops.push_back(
                manager->getOperationsDouble(solution.getComponentVariable(comp), d_hierarchy, /*get_unique*/ true));
        }
    }
    const bool repeated = d_cycle_type == W_CYCLE || d_cycle_type == F_CYCLE;
    const bool fmg = d_cycle_type == FMG_CYCLE;
    const int finest_warm_ln = fmg ? d_finest_ln : d_finest_ln - 1;
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        const std::string suffix = "::level_" + std::to_string(ln);
        if (ln > d_coarsest_ln)
        {
            const int residual_finest_ln = d_num_pre_sweeps > 0 ? ln : ln - 1;
            Pointer<SAMRAIVectorReal<NDIM, double>>& residual = d_residual_vectors[ln];
            if (residual && residual->getFinestLevelNumber() < residual_finest_ln)
            {
                free_vector_components(*residual);
                residual.setNull();
            }
            if (!residual)
            {
                residual = allocateRangeVector(rhs, d_object_name + "::residual" + suffix, residual_finest_ln);
            }
        }
        if ((repeated || (fmg && ln > d_coarsest_ln)) && ln <= finest_warm_ln)
        {
            if (!d_rhs_vectors[ln])
            {
                d_rhs_vectors[ln] = allocateRangeVector(rhs, d_object_name + "::rhs" + suffix, ln);
            }
            if (!d_correction_vectors[ln])
            {
                d_correction_vectors[ln] = allocateRangeVector(solution, d_object_name + "::correction" + suffix, ln);
            }
        }
        if (fmg && ln < d_finest_ln && !d_fmg_rhs_vectors[ln])
        {
            d_fmg_rhs_vectors[ln] = allocateRangeVector(rhs, d_object_name + "::fmg_rhs" + suffix, ln);
        }
    }
    if ((repeated || fmg) && finest_warm_ln >= d_coarsest_ln && d_finest_ln > d_coarsest_ln)
    {
        if (d_evaluation_vector && d_evaluation_vector->getFinestLevelNumber() < finest_warm_ln)
        {
            free_vector_components(*d_evaluation_vector);
            d_evaluation_vector.setNull();
        }
        if (!d_evaluation_vector)
        {
            d_evaluation_vector = allocateRangeVector(solution, d_object_name + "::evaluation", finest_warm_ln);
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
