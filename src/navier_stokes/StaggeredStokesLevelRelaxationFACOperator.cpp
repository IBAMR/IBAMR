// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2020 by the IBAMR developers
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

#include "ibamr/StaggeredStokesLevelRelaxationFACOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include "ibtk/PETScLevelSolver.h"

#include "ArrayData.h"
#include "Box.h"
#include "BoxList.h"
#include "CellData.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "ProcessorMapping.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

#include "petscksp.h"

#include <algorithm>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

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

    // Configure the coarse level solver.
    setCoarseSolverType(d_coarse_solver_type);

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
StaggeredStokesLevelRelaxationFACOperator::setSmootherType(const std::string& level_solver_type)
{
    StaggeredStokesFACPreconditionerStrategy::setSmootherType(level_solver_type);

    if (d_level_solver_type != level_solver_type)
    {
        d_level_solvers.clear();
    }
    d_level_solver_type = level_solver_type;

    return;
} // setSmootherType

void
StaggeredStokesLevelRelaxationFACOperator::smoothError(SAMRAIVectorRealNd<double>& error,
                                                       const SAMRAIVectorRealNd<double>& residual,
                                                       int level_num,
                                                       int num_sweeps,
                                                       bool /*performing_pre_sweeps*/,
                                                       bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBAMR_TIMER_START(t_smooth_error);

    Pointer<PatchLevelNd> level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);
    const int U_scratch_idx = d_side_scratch_idx;
    const int P_scratch_idx = d_cell_scratch_idx;

    Pointer<SAMRAIVectorRealNd<double> > e_level = getLevelSAMRAIVectorReal(error, level_num);
    Pointer<SAMRAIVectorRealNd<double> > r_level = getLevelSAMRAIVectorReal(residual, level_num);

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevelNd::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<PatchNd> patch = level->getPatch(p());

            Pointer<SideDataNd<double> > U_error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideDataNd<double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
#if !defined(NDEBUG)
            const BoxNd& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
            TBOX_ASSERT(U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                U_scratch_data->getArrayData(axis).copy(U_error_data->getArrayData(axis),
                                                        d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                        IntVectorNd(0));
            }

            Pointer<CellDataNd<double> > P_error_data = error.getComponentPatchData(1, *patch);
            Pointer<CellDataNd<double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
#if !defined(NDEBUG)
            const BoxNd& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
            TBOX_ASSERT(P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            P_scratch_data->getArrayData().copy(
                P_error_data->getArrayData(), d_patch_cell_bc_box_overlap[level_num][patch_counter], IntVectorNd(0));
        }
    }

    // Smooth the error by the specified number of sweeps.
    for (int isweep = 0; isweep < num_sweeps; ++isweep)
    {
        // Re-fill ghost cell data as needed.
        if (level_num > d_coarsest_ln)
        {
            if (isweep > 0)
            {
                // Copy the coarse-fine interface ghost cell values which are
                // cached in the scratch data into the error data.
                int patch_counter = 0;
                for (PatchLevelNd::Iterator p(level); p; p++, ++patch_counter)
                {
                    Pointer<PatchNd> patch = level->getPatch(p());

                    Pointer<SideDataNd<double> > U_error_data = error.getComponentPatchData(0, *patch);
                    Pointer<SideDataNd<double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        U_error_data->getArrayData(axis).copy(
                            U_scratch_data->getArrayData(axis),
                            d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                            IntVectorNd(0));
                    }

                    Pointer<CellDataNd<double> > P_error_data = error.getComponentPatchData(1, *patch);
                    Pointer<CellDataNd<double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
                    P_error_data->getArrayData().copy(P_scratch_data->getArrayData(),
                                                      d_patch_cell_bc_box_overlap[level_num][patch_counter],
                                                      IntVectorNd(0));
                }
            }

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension. This uses latest interior values on fine patch
            // but cached tangential values at the fine locations of the
            // coarse-fine interface
            d_U_cf_bdry_op->setPatchDataIndex(U_error_idx);
            d_P_cf_bdry_op->setPatchDataIndex(P_error_idx);
            const IntVectorNd& ratio = level->getRatioToCoarserLevel();
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                const IntVectorNd& ghost_width_to_fill = d_gcw;
                d_U_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
                d_P_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
            }
        }

        // Smooth the error on the level.
        Pointer<StaggeredStokesSolver> level_solver = d_level_solvers[level_num];
        level_solver->setSolutionTime(d_solution_time);
        level_solver->setTimeInterval(d_current_time, d_new_time);
        level_solver->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);
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
                if (!std::strcmp(ksp_type, "preonly")) initial_guess_nonzero = false;
            }
            p_level_solver->setInitialGuessNonzero(initial_guess_nonzero);
        }

        level_solver->solveSystem(*e_level, *r_level);
    }

    IBAMR_TIMER_STOP(t_smooth_error);
    return;
} // smoothError

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredStokesLevelRelaxationFACOperator::initializeOperatorStateSpecialized(
    const SAMRAIVectorRealNd<double>& /*solution*/,
    const SAMRAIVectorRealNd<double>& /*rhs*/,
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
            level_solver =
                StaggeredStokesSolverManager::getManager()->allocateSolver(d_level_solver_type,
                                                                           d_object_name + "::level_solver",
                                                                           d_level_solver_db,
                                                                           d_level_solver_default_options_prefix);
        }

        level_solver->setSolutionTime(d_solution_time);
        level_solver->setTimeInterval(d_current_time, d_new_time);
        level_solver->setVelocityPoissonSpecifications(d_U_problem_coefs);
        level_solver->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
        level_solver->setPhysicalBoundaryHelper(d_bc_helper);
        level_solver->setMaxIterations(d_level_solver_max_iterations);
        level_solver->setAbsoluteTolerance(d_level_solver_abs_residual_tol);
        level_solver->setRelativeTolerance(d_level_solver_rel_residual_tol);
        level_solver->setHomogeneousBc(true);
        level_solver->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);
        level_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, ln),
                                            *getLevelSAMRAIVectorReal(*d_rhs, ln));
    }

    // Nullify any fill pattern spec objects which maybe set by the base class.
    // However, we are NOT doing that with this class.

    // Get overlap information for setting patch boundary conditions.
    d_patch_side_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevelNd> level = d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_side_bc_box_overlap[ln].resize(num_local_patches);
        int patch_counter = 0;
        for (PatchLevelNd::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<PatchNd> patch = level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const BoxNd side_box = SideGeometryNd::toSideBox(patch_box, axis);
                const BoxNd side_ghost_box = BoxNd::grow(side_box, 1);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis] = BoxListNd(side_ghost_box);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis].removeIntersections(side_box);
            }
        }
    }

    d_patch_cell_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevelNd> level = d_hierarchy->getPatchLevel(ln);

        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_cell_bc_box_overlap[ln].resize(num_local_patches);

        int patch_counter = 0;
        for (PatchLevelNd::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<PatchNd> patch = level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();
            const BoxNd& ghost_box = BoxNd::grow(patch_box, 1);

            d_patch_cell_bc_box_overlap[ln][patch_counter] = BoxListNd(ghost_box);
            d_patch_cell_bc_box_overlap[ln][patch_counter].removeIntersections(patch_box);
        }
    }
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
        d_patch_side_bc_box_overlap[ln].resize(0);
        d_patch_cell_bc_box_overlap[ln].resize(0);
    }
    return;
} // deallocateOperatorStateSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
