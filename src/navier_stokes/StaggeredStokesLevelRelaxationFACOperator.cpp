// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2026 by the IBAMR developers
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

#include <ibamr/StaggeredStokesLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/StaggeredStokesSolverManager.h>
#include <ibamr/ibamr_utilities.h>

#include <ibtk/LinearSolver.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/solver_utilities.h>

#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/MemoryDatabase.h>
#include <tbox/Pointer.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

#include <petscksp.h>

#include <ArrayData.h>
#include <Box.h>
#include <CellData.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>

#include <algorithm>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_smooth_error;
static const int GHOST_CELL_WIDTH = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesLevelRelaxationFACOperator::StaggeredStokesLevelRelaxationFACOperator(
    const std::string& object_name,
    const Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : StaggeredStokesFACPreconditionerStrategy(object_name, GHOST_CELL_WIDTH, input_db, default_options_prefix),
      d_level_solver_default_options_prefix(default_options_prefix + "level_"),
      d_level_solver_db(new MemoryDatabase(object_name + "::level_solver_db"))
{
    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("level_solver_type")) d_level_solver_type = input_db->getString("level_solver_type");
        if (input_db->keyExists("level_solver_rel_residual_tol"))
            d_level_solver_rel_residual_tol = input_db->getDouble("level_solver_rel_residual_tol");
        if (input_db->keyExists("level_solver_abs_residual_tol"))
            d_level_solver_abs_residual_tol = input_db->getDouble("level_solver_abs_residual_tol");
        if (input_db->keyExists("level_solver_max_iterations"))
            d_level_solver_max_iterations = input_db->getInteger("level_solver_max_iterations");
        if (input_db->isDatabase("level_solver_db"))
        {
            d_level_solver_db = input_db->getDatabase("level_solver_db");
        }
    }

    // Setup Timers.
    IBAMR_DO_ONCE(t_smooth_error = TimerManager::getManager()->getTimer(
                      "IBAMR::StaggeredStokesLevelRelaxationFACOperator::smoothError()"););
    return;
} // StaggeredStokesLevelRelaxationFACOperator

StaggeredStokesLevelRelaxationFACOperator::~StaggeredStokesLevelRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~StaggeredStokesLevelRelaxationFACOperator

void
StaggeredStokesLevelRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
                                                       SAMRAIVectorReal<NDIM, double>& residual,
                                                       int level_num,
                                                       int num_sweeps,
                                                       bool performing_pre_sweeps,
                                                       bool performing_post_sweeps)
{
    if (num_sweeps == 0) return;

    IBAMR_TIMER_START(t_smooth_error);

    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);

    Pointer<SAMRAIVectorReal<NDIM, double>> e_level = getLevelSAMRAIVectorReal(error, level_num);
    Pointer<SAMRAIVectorReal<NDIM, double>> r_level = getLevelSAMRAIVectorReal(residual, level_num);

    prepareCoarseFineValues(U_error_idx, P_error_idx, level_num, performing_pre_sweeps, performing_post_sweeps);

    // Smooth the error by the specified number of sweeps.
    for (int isweep = 0; isweep < num_sweeps; ++isweep)
    {
        fillCoarseFineValues(level_num);

        // Smooth the error on the level.
        Pointer<StaggeredStokesSolver> level_solver = d_level_solvers[level_num];
        level_solver->setSolutionTime(d_solution_time);
        level_solver->setTimeInterval(d_current_time, d_new_time);
        level_solver->setComponentsHaveNullSpace(d_has_velocity_nullspace, d_has_pressure_nullspace);
        auto p_level_solver = dynamic_cast<LinearSolver*>(level_solver.getPointer());
        if (p_level_solver)
        {
            bool initial_guess_nonzero = true;

            auto p_petsc_solver = dynamic_cast<PETScKrylovLinearSolver*>(p_level_solver);
            auto p_petsc_level_solver = dynamic_cast<PETScLevelSolver*>(p_level_solver);

            if (p_petsc_solver || p_petsc_level_solver)
            {
                const KSP& petsc_ksp =
                    p_petsc_solver ? p_petsc_solver->getPETScKSP() : p_petsc_level_solver->getPETScKSP();
                KSPType ksp_type;
                KSPGetType(petsc_ksp, &ksp_type);
                if (std::string_view(ksp_type) == KSPPREONLY)
                {
                    // The solve starts from zero, which would discard the correction that post-smoothing keeps.
                    if (performing_post_sweeps)
                    {
                        TBOX_ERROR(d_object_name << "::smoothError():\n"
                                                 << "  post-smoothing with a level solver of KSP type preonly discards "
                                                    "the coarse-grid correction\n"
                                                 << "  use KSP type richardson" << std::endl);
                    }
                    initial_guess_nonzero = false;
                }
            }
            p_level_solver->setInitialGuessNonzero(initial_guess_nonzero);
        }

        level_solver->solveSystem(*e_level, *r_level);
    }

    finishCoarseFineValues(level_num);
    IBAMR_TIMER_STOP(t_smooth_error);
    return;
} // smoothError

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredStokesLevelRelaxationFACOperator::initializeOperatorStateSpecialized(
    const SAMRAIVectorReal<NDIM, double>& /*solution*/,
    const SAMRAIVectorReal<NDIM, double>& /*rhs*/,
    const int coarsest_reset_ln,
    const int finest_reset_ln)
{
    // Initialize the fine level solvers when needed.
    d_level_solvers.resize(d_finest_ln + 1);
    for (int ln = std::max(0, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<StaggeredStokesSolver>& level_solver = d_level_solvers[ln];
        if (!level_solver)
        {
            level_solver = StaggeredStokesSolverManager::getManager()->allocateSolver(
                d_level_solver_type,
                d_object_name + "::level_solver",
                d_level_solver_db,
                d_level_solver_default_options_prefix + std::to_string(ln) + "_");
            set_fixed_iteration_ksp_defaults(level_solver.getPointer(), d_level_solver_db);
        }

        level_solver->setSolutionTime(d_solution_time);
        level_solver->setTimeInterval(d_current_time, d_new_time);
        level_solver->setVelocityPoissonSpecifications(d_U_problem_coefs);
        level_solver->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
        if (d_bc_helper)
        {
            level_solver->setPhysicalBoundaryHelper(d_bc_helper);
        }
        level_solver->setMaxIterations(d_level_solver_max_iterations);
        level_solver->setAbsoluteTolerance(d_level_solver_abs_residual_tol);
        level_solver->setRelativeTolerance(d_level_solver_rel_residual_tol);
        level_solver->setHomogeneousBc(ALWAYS_HOMOGENEOUS_BC);
        level_solver->setComponentsHaveNullSpace(d_has_velocity_nullspace, d_has_pressure_nullspace);
        level_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, ln),
                                            *getLevelSAMRAIVectorReal(*d_rhs, ln));
    }

    // Nullify any fill pattern spec objects which maybe set by the base class.
    // However, we are NOT doing that with this class.
    return;
} // initializeOperatorStateSpecialized

void
StaggeredStokesLevelRelaxationFACOperator::deallocateOperatorStateSpecialized(const int coarsest_reset_ln,
                                                                              const int finest_reset_ln)
{
    if (!d_is_initialized) return;
    for (int ln = coarsest_reset_ln; ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        if (d_level_solvers[ln]) d_level_solvers[ln]->deallocateSolverState();
    }
    return;
} // deallocateOperatorStateSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
