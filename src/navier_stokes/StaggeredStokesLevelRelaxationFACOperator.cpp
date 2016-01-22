// Filename: StaggeredStokesLevelRelaxationFACOperator.cpp
// Created on 16 Mar 2015 by Amneet Bhalla
//
// Copyright (c) 2002-2015, Amneet Bhalla and Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <algorithm>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ArrayData.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "BoxList.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "ProcessorMapping.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/StaggeredStokesLevelRelaxationFACOperator.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include "ibtk/PETScLevelSolver.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscsys.h"
#include "petscvec.h"
#include "boost/array.hpp"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_smooth_error;
static const int GHOST_CELL_WIDTH = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesLevelRelaxationFACOperator::StaggeredStokesLevelRelaxationFACOperator(
    const std::string& object_name,
    const Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : StaggeredStokesFACPreconditionerStrategy(object_name, GHOST_CELL_WIDTH, input_db, default_options_prefix)
{
    // Set some default values for the level solvers.
    d_level_solver_type = StaggeredStokesSolverManager::PETSC_LEVEL_SOLVER;
    d_level_solver_default_options_prefix = default_options_prefix + "level_";
    d_level_solver_rel_residual_tol = 1.0e-5;
    d_level_solver_abs_residual_tol = 1.0e-50;
    d_level_solver_max_iterations = 1;
    d_level_solver_db = new MemoryDatabase(object_name + "::level_solver_db");

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
StaggeredStokesLevelRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
                                                       const SAMRAIVectorReal<NDIM, double>& residual,
                                                       int level_num,
                                                       int num_sweeps,
                                                       bool /*performing_pre_sweeps*/,
                                                       bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBAMR_TIMER_START(t_smooth_error);

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);
    const int U_scratch_idx = d_side_scratch_idx;
    const int P_scratch_idx = d_cell_scratch_idx;

    Pointer<SAMRAIVectorReal<NDIM, double> > e_level = getLevelSAMRAIVectorReal(error, level_num);
    Pointer<SAMRAIVectorReal<NDIM, double> > r_level = getLevelSAMRAIVectorReal(residual, level_num);

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > U_error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM, double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
#if !defined(NDEBUG)
            const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
            TBOX_ASSERT(U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                U_scratch_data->getArrayData(axis).copy(U_error_data->getArrayData(axis),
                                                        d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                        IntVector<NDIM>(0));
            }

            Pointer<CellData<NDIM, double> > P_error_data = error.getComponentPatchData(1, *patch);
            Pointer<CellData<NDIM, double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
#if !defined(NDEBUG)
            const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
            TBOX_ASSERT(P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            P_scratch_data->getArrayData().copy(P_error_data->getArrayData(),
                                                d_patch_cell_bc_box_overlap[level_num][patch_counter],
                                                IntVector<NDIM>(0));
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
                for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());

                    Pointer<SideData<NDIM, double> > U_error_data = error.getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM, double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        U_error_data->getArrayData(axis)
                            .copy(U_scratch_data->getArrayData(axis),
                                  d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                  IntVector<NDIM>(0));
                    }

                    Pointer<CellData<NDIM, double> > P_error_data = error.getComponentPatchData(1, *patch);
                    Pointer<CellData<NDIM, double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
                    P_error_data->getArrayData().copy(P_scratch_data->getArrayData(),
                                                      d_patch_cell_bc_box_overlap[level_num][patch_counter],
                                                      IntVector<NDIM>(0));
                }
            }

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension. This uses latest interior values on fine patch
            // but cached tangential values at the fine locations of the
            // coarse-fine interface
            d_U_cf_bdry_op->setPatchDataIndex(U_error_idx);
            d_P_cf_bdry_op->setPatchDataIndex(P_error_idx);
            const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const IntVector<NDIM>& ghost_width_to_fill = d_gcw;
                d_U_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
                d_P_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
            }
        }

        // Smooth the error on the level.
        Pointer<StaggeredStokesSolver> level_solver = d_level_solvers[level_num];
        level_solver->setSolutionTime(d_solution_time);
        level_solver->setTimeInterval(d_current_time, d_new_time);
        level_solver->setMaxIterations(d_level_solver_max_iterations);
        level_solver->setAbsoluteTolerance(d_level_solver_abs_residual_tol);
        level_solver->setRelativeTolerance(d_level_solver_rel_residual_tol);
        level_solver->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);
        LinearSolver* p_level_solver = dynamic_cast<LinearSolver*>(level_solver.getPointer());
        if (p_level_solver)
        {
            bool initial_guess_nonzero = true;

            PETScKrylovLinearSolver* p_petsc_solver = dynamic_cast<PETScKrylovLinearSolver*>(p_level_solver);
            PETScLevelSolver* p_petsc_level_solver = dynamic_cast<PETScLevelSolver*>(p_level_solver);

            if (p_petsc_solver || p_petsc_level_solver)
            {
                const KSP& petsc_ksp =
                    p_petsc_solver ? p_petsc_solver->getPETScKSP() : p_petsc_level_solver->getPETScKSP();
                KSPType ksp_type;
                KSPGetType(petsc_ksp, &ksp_type);
                if (!strcmp(ksp_type, "preonly")) initial_guess_nonzero = false;
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
    const SAMRAIVectorReal<NDIM, double>& /*solution*/,
    const SAMRAIVectorReal<NDIM, double>& /*rhs*/,
    const int coarsest_reset_ln,
    const int finest_reset_ln)
{
    // Initialize the fine level solvers when needed.
    d_level_solvers.resize(d_finest_ln + 1);
    for (int ln = std::max(1, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<StaggeredStokesSolver>& level_solver = d_level_solvers[ln];
        if (!level_solver)
        {
            std::ostringstream level_solver_stream;
            level_solver_stream << d_level_solver_default_options_prefix;
            level_solver = StaggeredStokesSolverManager::getManager()->allocateSolver(
                d_level_solver_type, d_object_name + "::level_solver", d_level_solver_db, level_solver_stream.str());
        }

        level_solver->setSolutionTime(d_solution_time);
        level_solver->setTimeInterval(d_current_time, d_new_time);
        level_solver->setVelocityPoissonSpecifications(d_U_problem_coefs);
        level_solver->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
        level_solver->setPhysicalBoundaryHelper(d_bc_helper);
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
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_side_bc_box_overlap[ln].resize(num_local_patches);
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                const Box<NDIM> side_ghost_box = Box<NDIM>::grow(side_box, 1);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis] = BoxList<NDIM>(side_ghost_box);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis].removeIntersections(side_box);
            }
        }
    }

    d_patch_cell_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_cell_bc_box_overlap[ln].resize(num_local_patches);

        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Box<NDIM>& ghost_box = Box<NDIM>::grow(patch_box, 1);

            d_patch_cell_bc_box_overlap[ln][patch_counter] = BoxList<NDIM>(ghost_box);
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
