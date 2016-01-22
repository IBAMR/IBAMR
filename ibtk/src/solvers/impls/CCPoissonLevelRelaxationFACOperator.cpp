// Filename: CCPoissonLevelRelaxationFACOperator.cpp
// Created on 8 Mar 2015 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#include <stddef.h>
#include <algorithm>
#include <functional>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ArrayData.h"
#include "Box.h"
#include "BoxList.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "CoarsenOperator.h"
#include "HierarchyCellDataOpsReal.h"
#include "IBTK_config.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchDescriptor.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "ProcessorMapping.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideIndex.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "VariableFillPattern.h"
#include "boost/array.hpp"
#include "ibtk/CCPoissonLevelRelaxationFACOperator.h"
#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartCellDoubleCubicCoarsen.h"
#include "ibtk/CartCellDoubleQuadraticCFInterpolation.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include "ibtk/PETScLevelSolver.h"
#include "ibtk/PoissonFACPreconditionerStrategy.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "petscmat.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_smooth_error;
static Timer* t_solve_coarsest_level;
static Timer* t_compute_residual;

// Default data depth.
static const int DEFAULT_DATA_DEPTH = 1;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries; used only to evaluate
// composite grid residuals.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells; used only to evaluate composite grid residuals.
static const bool CONSISTENT_TYPE_2_BDRY = false;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCPoissonLevelRelaxationFACOperator::CCPoissonLevelRelaxationFACOperator(const std::string& object_name,
                                                                         const Pointer<Database> input_db,
                                                                         const std::string& default_options_prefix)
    : PoissonFACPreconditionerStrategy(
          object_name,
          new CellVariable<NDIM, double>(object_name + "::cell_scratch", DEFAULT_DATA_DEPTH),
          CELLG,
          input_db,
          default_options_prefix),
      d_level_solver_type(),
      d_level_solvers(),
      d_level_solver_db(),
      d_coarse_solver(NULL),
      d_coarse_solver_db()
{
    // Set some default values.
    d_prolongation_method = "LINEAR_REFINE";
    d_restriction_method = "CONSERVATIVE_COARSEN";
    d_level_solver_type = CCPoissonSolverManager::PETSC_LEVEL_SOLVER;
    d_level_solver_default_options_prefix = default_options_prefix + "level_";
    d_level_solver_rel_residual_tol = 1.0e-5;
    d_level_solver_abs_residual_tol = 1.0e-50;
    d_level_solver_max_iterations = 1;
    d_level_solver_db = new MemoryDatabase(object_name + "::level_solver_db");
    d_coarse_solver_type = CCPoissonSolverManager::PETSC_LEVEL_SOLVER;
    d_coarse_solver_default_options_prefix = default_options_prefix + "level_0_";
    d_coarse_solver_rel_residual_tol = 1.0e-5;
    d_coarse_solver_abs_residual_tol = 1.0e-50;
    d_coarse_solver_max_iterations = 1;
    d_coarse_solver_db = new MemoryDatabase(object_name + "::coarse_solver_db");

    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("smoother_type")) d_smoother_type = input_db->getString("smoother_type");
        if (input_db->keyExists("prolongation_method"))
            d_prolongation_method = input_db->getString("prolongation_method");
        if (input_db->keyExists("restriction_method")) d_restriction_method = input_db->getString("restriction_method");
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
        if (input_db->keyExists("coarse_solver_type")) d_coarse_solver_type = input_db->getString("coarse_solver_type");
        if (input_db->keyExists("coarse_solver_rel_residual_tol"))
            d_coarse_solver_rel_residual_tol = input_db->getDouble("coarse_solver_rel_residual_tol");
        if (input_db->keyExists("coarse_solver_abs_residual_tol"))
            d_coarse_solver_abs_residual_tol = input_db->getDouble("coarse_solver_abs_residual_tol");
        if (input_db->keyExists("coarse_solver_max_iterations"))
            d_coarse_solver_max_iterations = input_db->getInteger("coarse_solver_max_iterations");
        if (input_db->isDatabase("coarse_solver_db"))
        {
            d_coarse_solver_db = input_db->getDatabase("coarse_solver_db");
        }
    }

    // Configure the coarse level solver.
    setCoarseSolverType(d_coarse_solver_type);

    // Setup Timers.
    IBTK_DO_ONCE(t_smooth_error =
                     TimerManager::getManager()->getTimer("IBTK::CCPoissonLevelRelaxationFACOperator::smoothError()");
                 t_solve_coarsest_level = TimerManager::getManager()->getTimer(
                     "IBTK::CCPoissonLevelRelaxationFACOperator::solveCoarsestLevel()");
                 t_compute_residual = TimerManager::getManager()->getTimer(
                     "IBTK::CCPoissonLevelRelaxationFACOperator::computeResidual()"););
    return;
}

CCPoissonLevelRelaxationFACOperator::~CCPoissonLevelRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}

void
CCPoissonLevelRelaxationFACOperator::setSmootherType(const std::string& level_solver_type)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSmootherType():\n"
                                 << "  cannot be called while operator state is initialized"
                                 << std::endl);
    }
    if (d_level_solver_type != level_solver_type)
    {
        d_level_solvers.clear();
    }
    d_level_solver_type = level_solver_type;
    return;
}

void
CCPoissonLevelRelaxationFACOperator::setCoarseSolverType(const std::string& coarse_solver_type)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setCoarseSolverType():\n"
                                 << "  cannot be called while operator state is initialized"
                                 << std::endl);
    }
    if (d_coarse_solver_type != coarse_solver_type) d_coarse_solver.setNull();
    d_coarse_solver_type = coarse_solver_type;
    if (!d_coarse_solver)
    {
        d_coarse_solver = CCPoissonSolverManager::getManager()->allocateSolver(d_coarse_solver_type,
                                                                               d_object_name + "::coarse_solver",
                                                                               d_coarse_solver_db,
                                                                               d_coarse_solver_default_options_prefix);
    }
    return;
}

void
CCPoissonLevelRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
                                                 const SAMRAIVectorReal<NDIM, double>& residual,
                                                 int level_num,
                                                 int num_sweeps,
                                                 bool /*performing_pre_sweeps*/,
                                                 bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBTK_TIMER_START(t_smooth_error);

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int error_idx = error.getComponentDescriptorIndex(0);
    const int scratch_idx = d_scratch_idx;

    Pointer<SAMRAIVectorReal<NDIM, double> > e_level = getLevelSAMRAIVectorReal(error, level_num);
    Pointer<SAMRAIVectorReal<NDIM, double> > r_level = getLevelSAMRAIVectorReal(residual, level_num);

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > error_data = error.getComponentPatchData(0, *patch);
            Pointer<CellData<NDIM, double> > scratch_data = patch->getPatchData(scratch_idx);
#if !defined(NDEBUG)
            TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
            scratch_data->getArrayData().copy(
                error_data->getArrayData(), d_patch_bc_box_overlap[level_num][patch_counter], IntVector<NDIM>(0));
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
                    Pointer<CellData<NDIM, double> > error_data = error.getComponentPatchData(0, *patch);
                    Pointer<CellData<NDIM, double> > scratch_data = patch->getPatchData(scratch_idx);
#if !defined(NDEBUG)
                    TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    error_data->getArrayData().copy(scratch_data->getArrayData(),
                                                    d_patch_bc_box_overlap[level_num][patch_counter],
                                                    IntVector<NDIM>(0));
                }
            }

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension.
            d_cf_bdry_op->setPatchDataIndex(error_idx);
            const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const IntVector<NDIM>& ghost_width_to_fill = d_gcw;
                d_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
            }
        }

        // Smooth the error.
        Pointer<PoissonSolver> level_solver = d_level_solvers[level_num];
        level_solver->setSolutionTime(d_solution_time);
        level_solver->setTimeInterval(d_current_time, d_new_time);
        level_solver->setMaxIterations(d_level_solver_max_iterations);
        level_solver->setAbsoluteTolerance(d_level_solver_abs_residual_tol);
        level_solver->setRelativeTolerance(d_level_solver_rel_residual_tol);
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
    IBTK_TIMER_STOP(t_smooth_error);
    return;
}

bool
CCPoissonLevelRelaxationFACOperator::solveCoarsestLevel(SAMRAIVectorReal<NDIM, double>& error,
                                                        const SAMRAIVectorReal<NDIM, double>& residual,
                                                        int coarsest_ln)
{
    IBTK_TIMER_START(t_solve_coarsest_level);
#if !defined(NDEBUG)
    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
    TBOX_ASSERT(d_coarse_solver);
#endif
    Pointer<SAMRAIVectorReal<NDIM, double> > e_level = getLevelSAMRAIVectorReal(error, coarsest_ln);
    Pointer<SAMRAIVectorReal<NDIM, double> > r_level = getLevelSAMRAIVectorReal(residual, coarsest_ln);
    d_coarse_solver->setSolutionTime(d_solution_time);
    d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
    d_coarse_solver->setMaxIterations(d_coarse_solver_max_iterations);
    d_coarse_solver->setAbsoluteTolerance(d_coarse_solver_abs_residual_tol);
    d_coarse_solver->setRelativeTolerance(d_coarse_solver_rel_residual_tol);
    LinearSolver* p_coarse_solver = dynamic_cast<LinearSolver*>(d_coarse_solver.getPointer());
    if (p_coarse_solver)
    {
        bool initial_guess_nonzero = true;
        PETScKrylovLinearSolver* p_petsc_solver = dynamic_cast<PETScKrylovLinearSolver*>(p_coarse_solver);
        PETScLevelSolver* p_petsc_level_solver = dynamic_cast<PETScLevelSolver*>(p_coarse_solver);
        if (p_petsc_solver || p_petsc_level_solver)
        {
            const KSP& petsc_ksp = p_petsc_solver ? p_petsc_solver->getPETScKSP() : p_petsc_level_solver->getPETScKSP();
            KSPType ksp_type;
            KSPGetType(petsc_ksp, &ksp_type);

            if (!strcmp(ksp_type, "preonly")) initial_guess_nonzero = false;
        }
        p_coarse_solver->setInitialGuessNonzero(initial_guess_nonzero);
    }

    d_coarse_solver->solveSystem(*e_level, *r_level);
    IBTK_TIMER_STOP(t_solve_coarsest_level);
    return true;
}

void
CCPoissonLevelRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
                                                     const SAMRAIVectorReal<NDIM, double>& solution,
                                                     const SAMRAIVectorReal<NDIM, double>& rhs,
                                                     int coarsest_level_num,
                                                     int finest_level_num)
{
    IBTK_TIMER_START(t_compute_residual);

    const int res_idx = residual.getComponentDescriptorIndex(0);
    const int sol_idx = solution.getComponentDescriptorIndex(0);
    const int rhs_idx = rhs.getComponentDescriptorIndex(0);

    const Pointer<CellVariable<NDIM, double> > res_var = residual.getComponentVariable(0);
    const Pointer<CellVariable<NDIM, double> > sol_var = solution.getComponentVariable(0);
    const Pointer<CellVariable<NDIM, double> > rhs_var = rhs.getComponentVariable(0);

    // Fill ghost-cell values.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    Pointer<CellNoCornersFillPattern> fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    InterpolationTransactionComponent transaction_comp(sol_idx,
                                                       DATA_REFINE_TYPE,
                                                       USE_CF_INTERPOLATION,
                                                       DATA_COARSEN_TYPE,
                                                       BDRY_EXTRAP_TYPE,
                                                       CONSISTENT_TYPE_2_BDRY,
                                                       d_bc_coefs,
                                                       fill_pattern);
    if (d_level_bdry_fill_ops[finest_level_num])
    {
        d_level_bdry_fill_ops[finest_level_num]->resetTransactionComponent(transaction_comp);
    }
    else
    {
        d_level_bdry_fill_ops[finest_level_num] = new HierarchyGhostCellInterpolation();
        d_level_bdry_fill_ops[finest_level_num]->initializeOperatorState(
            transaction_comp, d_hierarchy, coarsest_level_num, finest_level_num);
    }
    d_level_bdry_fill_ops[finest_level_num]->setHomogeneousBc(true);
    d_level_bdry_fill_ops[finest_level_num]->fillData(d_solution_time);
    InterpolationTransactionComponent default_transaction_comp(d_solution->getComponentDescriptorIndex(0),
                                                               DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_bc_coefs,
                                                               fill_pattern);
    d_level_bdry_fill_ops[finest_level_num]->resetTransactionComponent(default_transaction_comp);

    // Compute the residual, r = f - A*u.
    if (!d_level_math_ops[finest_level_num])
    {
        std::ostringstream stream;
        stream << d_object_name << "::hier_math_ops_" << finest_level_num;
        d_level_math_ops[finest_level_num] =
            new HierarchyMathOps(stream.str(), d_hierarchy, coarsest_level_num, finest_level_num);
    }
    d_level_math_ops[finest_level_num]->laplace(
        res_idx, res_var, d_poisson_spec, sol_idx, sol_var, NULL, d_solution_time);
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    hier_cc_data_ops.axpy(res_idx, -1.0, res_idx, rhs_idx, false);

    IBTK_TIMER_STOP(t_compute_residual);
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void
CCPoissonLevelRelaxationFACOperator::initializeOperatorStateSpecialized(const SAMRAIVectorReal<NDIM, double>& solution,
                                                                        const SAMRAIVectorReal<NDIM, double>& rhs,
                                                                        const int coarsest_reset_ln,
                                                                        const int finest_reset_ln)
{
    // Setup solution and rhs vectors.
    Pointer<CellVariable<NDIM, double> > solution_var = solution.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > rhs_var = rhs.getComponentVariable(0);

    Pointer<CellDataFactory<NDIM, double> > solution_pdat_fac = solution_var->getPatchDataFactory();
    Pointer<CellDataFactory<NDIM, double> > rhs_pdat_fac = rhs_var->getPatchDataFactory();

#if !defined(NDEBUG)
    TBOX_ASSERT(solution_var);
    TBOX_ASSERT(rhs_var);
    TBOX_ASSERT(solution_pdat_fac);
    TBOX_ASSERT(rhs_pdat_fac);
#endif

    if (solution_pdat_fac->getDefaultDepth() != rhs_pdat_fac->getDefaultDepth())
    {
        TBOX_ERROR("CCPoissonLevelRelaxationFACOperator::initializeOperatorState()\n"
                   << "  solution and rhs vectors must have the same data depths\n"
                   << "  solution data depth = "
                   << solution_pdat_fac->getDefaultDepth()
                   << "\n"
                   << "  rhs      data depth = "
                   << rhs_pdat_fac->getDefaultDepth()
                   << std::endl);
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellDataFactory<NDIM, double> > scratch_pdat_fac =
        var_db->getPatchDescriptor()->getPatchDataFactory(d_scratch_idx);
    scratch_pdat_fac->setDefaultDepth(solution_pdat_fac->getDefaultDepth());

    // Initialize the fine level solvers when needed.
    d_level_solvers.resize(d_finest_ln + 1);
    for (int ln = std::max(1, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<PoissonSolver>& level_solver = d_level_solvers[ln];
        if (!level_solver)
        {
            std::ostringstream level_solver_stream;
            level_solver_stream << d_level_solver_default_options_prefix << ln << "_";
            level_solver = CCPoissonSolverManager::getManager()->allocateSolver(
                d_level_solver_type, d_object_name + "::level_solver", d_level_solver_db, level_solver_stream.str());
        }

        level_solver->setSolutionTime(d_solution_time);
        level_solver->setTimeInterval(d_current_time, d_new_time);
        level_solver->setPoissonSpecifications(d_poisson_spec);
        level_solver->setPhysicalBcCoefs(d_bc_coefs);
        level_solver->setHomogeneousBc(true);
        level_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, ln),
                                            *getLevelSAMRAIVectorReal(*d_rhs, ln));
    }

    // Initialize the coarse level solvers when needed.
    if (coarsest_reset_ln == d_coarsest_ln)
    {
        // Note that since the coarse level solver is solving for the error, it
        // must always employ homogeneous boundary conditions.
        d_coarse_solver->setSolutionTime(d_solution_time);
        d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
        d_coarse_solver->setPoissonSpecifications(d_poisson_spec);
        d_coarse_solver->setPhysicalBcCoefs(d_bc_coefs);
        d_coarse_solver->setHomogeneousBc(true);
        d_coarse_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, d_coarsest_ln),
                                               *getLevelSAMRAIVectorReal(*d_rhs, d_coarsest_ln));
    }

    PetscViewer matlab_viewer;
    const int debug_ln = d_coarsest_ln;
    const KSP& level_ksp = debug_ln == 0 ?
                               dynamic_cast<PETScLevelSolver*>(d_coarse_solver.getPointer())->getPETScKSP() :
                               dynamic_cast<PETScLevelSolver*>(d_level_solvers[debug_ln].getPointer())->getPETScKSP();
    KSPSetUp(level_ksp);
    Mat level_mat;
    KSPGetOperators(level_ksp, &level_mat, NULL);
    std::ostringstream level_mat_filename;
    level_mat_filename << "level_" << debug_ln << "_mat";
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, level_mat_filename.str().c_str(), FILE_MODE_WRITE, &matlab_viewer);
    PetscViewerSetFormat(matlab_viewer, PETSC_VIEWER_NATIVE);
    MatView(level_mat, matlab_viewer);

    PC pc;
    KSPGetPC(level_ksp, &pc);
    KSP* subksp;
    PetscInt nlocal, first;
    PCASMGetSubKSP(pc, &nlocal, &first, &subksp);
    for (int i = 0; i < nlocal; i++)
    {
        Mat sub_mat;
        KSPGetOperators(subksp[i], &sub_mat, PETSC_NULL);

        std::ostringstream filename;
        filename << "level_" << debug_ln << "_sub_mat_" << i;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.str().c_str(), FILE_MODE_WRITE, &matlab_viewer);
        PetscViewerSetFormat(matlab_viewer, PETSC_VIEWER_NATIVE);
        MatView(sub_mat, matlab_viewer);
    }

    // Destroy Petsc reader
    PetscViewerDestroy(&matlab_viewer);

    // Setup specialized transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    IBTK_DO_ONCE(geometry->addSpatialCoarsenOperator(new CartCellDoubleCubicCoarsen()););

    // Setup coarse-fine interface and physical boundary operators.
    d_cf_bdry_op = new CartCellDoubleQuadraticCFInterpolation();
    d_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_cf_bdry_op->setPatchDataIndex(d_scratch_idx);
    d_cf_bdry_op->setPatchHierarchy(d_hierarchy);
    d_bc_op = new CartCellRobinPhysBdryOp(d_scratch_idx, d_bc_coefs, false);

    // Setup fill pattern spec objects.
    if (d_poisson_spec.dIsConstant())
    {
        d_op_stencil_fill_pattern = new CellNoCornersFillPattern(CELLG, true, false, false);
    }
    else
    {
        d_op_stencil_fill_pattern = NULL;
    }

    // Get overlap information for setting patch boundary conditions.
    d_patch_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_bc_box_overlap[ln].resize(num_local_patches);
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Box<NDIM>& ghost_box = Box<NDIM>::grow(patch_box, 1);
            d_patch_bc_box_overlap[ln][patch_counter] = BoxList<NDIM>(ghost_box);
            d_patch_bc_box_overlap[ln][patch_counter].removeIntersections(patch_box);
        }
    }
    return;
}

void
CCPoissonLevelRelaxationFACOperator::deallocateOperatorStateSpecialized(const int /*coarsest_reset_ln*/,
                                                                        const int /*finest_reset_ln*/)
{
    if (!d_is_initialized) return;
    if (!d_in_initialize_operator_state)
    {
        d_patch_bc_box_overlap.clear();
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            if (d_level_solvers[ln]) d_level_solvers[ln]->deallocateSolverState();
        }
        if (d_coarse_solver) d_coarse_solver->deallocateSolverState();
    }
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
