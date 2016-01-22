// Filename: CCPoissonPointRelaxationFACOperator.cpp
// Created on 10 Feb 2005 by Boyce Griffith
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
#include "ibtk/CCPoissonPointRelaxationFACOperator.h"
#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartCellDoubleCubicCoarsen.h"
#include "ibtk/CartCellDoubleQuadraticCFInterpolation.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LinearSolver.h"
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

// FORTRAN ROUTINES
#if (NDIM == 2)
#define GS_SMOOTH_FC IBTK_FC_FUNC(gssmooth2d, GSSMOOTH2D)
#define RB_GS_SMOOTH_FC IBTK_FC_FUNC(rbgssmooth2d, RBGSSMOOTH2D)
#endif
#if (NDIM == 3)
#define GS_SMOOTH_FC IBTK_FC_FUNC(gssmooth3d, GSSMOOTH3D)
#define RB_GS_SMOOTH_FC IBTK_FC_FUNC(rbgssmooth3d, RBGSSMOOTH3D)
#endif

// Function interfaces
extern "C" {
void GS_SMOOTH_FC(double* U,
                  const int& U_gcw,
                  const double& alpha,
                  const double& beta,
                  const double* F,
                  const int& F_gcw,
                  const int& ilower0,
                  const int& iupper0,
                  const int& ilower1,
                  const int& iupper1,
#if (NDIM == 3)
                  const int& ilower2,
                  const int& iupper2,
#endif
                  const double* dx);

void RB_GS_SMOOTH_FC(double* U,
                     const int& U_gcw,
                     const double& alpha,
                     const double& beta,
                     const double* F,
                     const int& F_gcw,
                     const int& ilower0,
                     const int& iupper0,
                     const int& ilower1,
                     const int& iupper1,
#if (NDIM == 3)
                     const int& ilower2,
                     const int& iupper2,
#endif
                     const double* dx,
                     const int& red_or_black);
}

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

struct IndexComp : std::binary_function<Index<NDIM>, Index<NDIM>, bool>
{
    inline bool operator()(const Index<NDIM>& lhs, const Index<NDIM>& rhs) const
    {
        return ((lhs(0) < rhs(0))
#if (NDIM > 1)
                ||
                (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                ||
                (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
                    );
    } // operator()
};

enum SmootherType
{
    PATCH_GAUSS_SEIDEL,
    PROCESSOR_GAUSS_SEIDEL,
    RED_BLACK_GAUSS_SEIDEL,
    UNKNOWN = -1
};

inline SmootherType
get_smoother_type(const std::string& smoother_type_string)
{
    if (smoother_type_string == "PATCH_GAUSS_SEIDEL") return PATCH_GAUSS_SEIDEL;
    if (smoother_type_string == "PROCESSOR_GAUSS_SEIDEL") return PROCESSOR_GAUSS_SEIDEL;
    if (smoother_type_string == "RED_BLACK_GAUSS_SEIDEL")
        return RED_BLACK_GAUSS_SEIDEL;
    else
        return UNKNOWN;
} // get_smoother_type

inline bool
use_red_black_ordering(SmootherType smoother_type)
{
    if (smoother_type == RED_BLACK_GAUSS_SEIDEL)
    {
        return true;
    }
    else
    {
        return false;
    }
} // use_red_black_ordering

inline bool
do_local_data_update(SmootherType smoother_type)
{
    if (smoother_type == PROCESSOR_GAUSS_SEIDEL || smoother_type == RED_BLACK_GAUSS_SEIDEL)
    {
        return true;
    }
    else
    {
        return false;
    }
} // do_local_data_update
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCPoissonPointRelaxationFACOperator::CCPoissonPointRelaxationFACOperator(const std::string& object_name,
                                                                         const Pointer<Database> input_db,
                                                                         const std::string& default_options_prefix)
    : PoissonFACPreconditionerStrategy(
          object_name,
          new CellVariable<NDIM, double>(object_name + "::cell_scratch", DEFAULT_DATA_DEPTH),
          CELLG,
          input_db,
          default_options_prefix),
      d_coarse_solver(NULL),
      d_coarse_solver_db(),
      d_patch_bc_box_overlap(),
      d_patch_neighbor_overlap()
{
    // Set some default values.
    d_smoother_type = "PATCH_GAUSS_SEIDEL";
    d_prolongation_method = "LINEAR_REFINE";
    d_restriction_method = "CONSERVATIVE_COARSEN";
    d_coarse_solver_type = CCPoissonSolverManager::HYPRE_LEVEL_SOLVER;
    d_coarse_solver_rel_residual_tol = 1.0e-5;
    d_coarse_solver_abs_residual_tol = 1.0e-50;
    d_coarse_solver_max_iterations = 1;
    d_coarse_solver_db = new MemoryDatabase(object_name + "::coarse_solver_db");
    d_coarse_solver_db->putString("solver_type", "PFMG");
    d_coarse_solver_db->putInteger("num_pre_relax_steps", 0);
    d_coarse_solver_db->putInteger("num_post_relax_steps", 2);

    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("smoother_type")) d_smoother_type = input_db->getString("smoother_type");
        if (input_db->keyExists("prolongation_method"))
            d_prolongation_method = input_db->getString("prolongation_method");
        if (input_db->keyExists("restriction_method")) d_restriction_method = input_db->getString("restriction_method");
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
        if (input_db->isDatabase("bottom_solver"))
        {
            tbox::pout << "WARNING: ``bottom_solver'' input entry is no longer used by class "
                          "CCPoissonPointRelaxationFACOperator.\n"
                       << "         use ``coarse_solver_db'' input entry instead.\n";
        }
        if (input_db->isDatabase("hypre_solver"))
        {
            tbox::pout << "WARNING: ``hypre_solver'' input entry is no longer used by class "
                          "CCPoissonPointRelaxationFACOperator.\n"
                       << "         use ``coarse_solver_db'' input entry instead.\n";
        }
        if (input_db->isDatabase("petsc_solver"))
        {
            tbox::pout << "WARNING: ``petsc_solver'' input entry is no longer used by class "
                          "CCPoissonPointRelaxationFACOperator.\n"
                       << "         use ``coarse_solver_db'' input entry instead.\n";
        }
    }

    // Configure the coarse level solver.
    setCoarseSolverType(d_coarse_solver_type);

    // Setup Timers.
    IBTK_DO_ONCE(t_smooth_error =
                     TimerManager::getManager()->getTimer("IBTK::CCPoissonPointRelaxationFACOperator::smoothError()");
                 t_solve_coarsest_level = TimerManager::getManager()->getTimer(
                     "IBTK::CCPoissonPointRelaxationFACOperator::solveCoarsestLevel()");
                 t_compute_residual = TimerManager::getManager()->getTimer(
                     "IBTK::CCPoissonPointRelaxationFACOperator::computeResidual()"););
    return;
} // CCPoissonPointRelaxationFACOperator

CCPoissonPointRelaxationFACOperator::~CCPoissonPointRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~CCPoissonPointRelaxationFACOperator

void
CCPoissonPointRelaxationFACOperator::setSmootherType(const std::string& smoother_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(get_smoother_type(smoother_type) != UNKNOWN);
#endif
    d_smoother_type = smoother_type;
    return;
} // setSmootherType

void
CCPoissonPointRelaxationFACOperator::setCoarseSolverType(const std::string& coarse_solver_type)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setCoarseSolverType():\n"
                                 << "  cannot be called while operator state is initialized"
                                 << std::endl);
    }
    if (d_coarse_solver_type != coarse_solver_type) d_coarse_solver.setNull();
    d_coarse_solver_type = coarse_solver_type;
    if (get_smoother_type(d_coarse_solver_type) == UNKNOWN && !d_coarse_solver)
    {
        d_coarse_solver = CCPoissonSolverManager::getManager()->allocateSolver(d_coarse_solver_type,
                                                                               d_object_name + "::coarse_solver",
                                                                               d_coarse_solver_db,
                                                                               d_coarse_solver_default_options_prefix);
    }
    return;
} // setCoarseSolverType

void
CCPoissonPointRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
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

    // Determine the smoother type.
    const std::string& smoother_type_string = (level_num == d_coarsest_ln ? d_coarse_solver_type : d_smoother_type);
    const SmootherType smoother_type = get_smoother_type(smoother_type_string);
#if !defined(NDEBUG)
    TBOX_ASSERT(smoother_type != UNKNOWN);
#endif
    const bool red_black_ordering = use_red_black_ordering(smoother_type);
    const bool update_local_data = do_local_data_update(smoother_type);

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
            const Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
            TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
            scratch_data->getArrayData().copy(
                error_data->getArrayData(), d_patch_bc_box_overlap[level_num][patch_counter], IntVector<NDIM>(0));
        }
    }

    // Smooth the error by the specified number of sweeps.
    if (red_black_ordering) num_sweeps *= 2;
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
                    const Box<NDIM>& ghost_box = error_data->getGhostBox();
                    TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
                    TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    error_data->getArrayData().copy(scratch_data->getArrayData(),
                                                    d_patch_bc_box_overlap[level_num][patch_counter],
                                                    IntVector<NDIM>(0));
                }

                // Fill the non-coarse-fine interface ghost cell values.
                xeqScheduleGhostFillNoCoarse(error_idx, level_num);
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
        else if (isweep > 0)
        {
            xeqScheduleGhostFillNoCoarse(error_idx, level_num);
        }

        // Smooth the error on the patches.
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > error_data = error.getComponentPatchData(0, *patch);
            Pointer<CellData<NDIM, double> > residual_data = residual.getComponentPatchData(0, *patch);
#if !defined(NDEBUG)
            const Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == residual_data->getGhostBox());
            TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(residual_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(error_data->getDepth() == residual_data->getDepth());
#endif
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            // Copy updated values from neighboring local patches.
            if (update_local_data)
            {
                const std::map<int, Box<NDIM> > neighbor_overlap = d_patch_neighbor_overlap[level_num][patch_counter];
                for (std::map<int, Box<NDIM> >::const_iterator cit = neighbor_overlap.begin();
                     cit != neighbor_overlap.end();
                     ++cit)
                {
                    const int src_patch_num = cit->first;
                    const Box<NDIM>& overlap = cit->second;
                    Pointer<Patch<NDIM> > src_patch = level->getPatch(src_patch_num);
                    Pointer<CellData<NDIM, double> > src_error_data = error.getComponentPatchData(0, *src_patch);
                    error_data->getArrayData().copy(src_error_data->getArrayData(), overlap, IntVector<NDIM>(0));
                }
            }

            // Smooth the error for each data depth.
            //
            // NOTE: Since the boundary conditions are handled "implicitly" by
            // setting ghost cell values, we can re-use the same patch operators
            // for each data depth even if different boundary conditions are
            // imposed on different components of the vector-valued solution
            // data.
            const double& alpha = d_poisson_spec.getDConstant();
            const double& beta = d_poisson_spec.cIsZero() ? 0.0 : d_poisson_spec.getCConstant();
            for (int depth = 0; depth < error_data->getDepth(); ++depth)
            {
                double* const U = error_data->getPointer(depth);
                const int U_ghosts = (error_data->getGhostCellWidth()).max();
                const double* const F = residual_data->getPointer(depth);
                const int F_ghosts = (residual_data->getGhostCellWidth()).max();
                if (red_black_ordering)
                {
                    int red_or_black = isweep % 2; // "red" = 0, "black" = 1
                    RB_GS_SMOOTH_FC(U,
                                    U_ghosts,
                                    alpha,
                                    beta,
                                    F,
                                    F_ghosts,
                                    patch_box.lower(0),
                                    patch_box.upper(0),
                                    patch_box.lower(1),
                                    patch_box.upper(1),
#if (NDIM == 3)
                                    patch_box.lower(2),
                                    patch_box.upper(2),
#endif
                                    dx,
                                    red_or_black);
                }
                else
                {
                    GS_SMOOTH_FC(U,
                                 U_ghosts,
                                 alpha,
                                 beta,
                                 F,
                                 F_ghosts,
                                 patch_box.lower(0),
                                 patch_box.upper(0),
                                 patch_box.lower(1),
                                 patch_box.upper(1),
#if (NDIM == 3)
                                 patch_box.lower(2),
                                 patch_box.upper(2),
#endif
                                 dx);
                }
            }
        }
    }
    IBTK_TIMER_STOP(t_smooth_error);
    return;
} // smoothError

bool
CCPoissonPointRelaxationFACOperator::solveCoarsestLevel(SAMRAIVectorReal<NDIM, double>& error,
                                                        const SAMRAIVectorReal<NDIM, double>& residual,
                                                        int coarsest_ln)
{
    IBTK_TIMER_START(t_solve_coarsest_level);
#if !defined(NDEBUG)
    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
#endif
    if (d_coarse_solver)
    {
        d_coarse_solver->setSolutionTime(d_solution_time);
        d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
        d_coarse_solver->setMaxIterations(d_coarse_solver_max_iterations);
        d_coarse_solver->setAbsoluteTolerance(d_coarse_solver_abs_residual_tol);
        d_coarse_solver->setRelativeTolerance(d_coarse_solver_rel_residual_tol);
        LinearSolver* p_coarse_solver = dynamic_cast<LinearSolver*>(d_coarse_solver.getPointer());
        if (p_coarse_solver) p_coarse_solver->setInitialGuessNonzero(true);
        d_coarse_solver->solveSystem(*getLevelSAMRAIVectorReal(error, d_coarsest_ln),
                                     *getLevelSAMRAIVectorReal(residual, d_coarsest_ln));
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(get_smoother_type(d_coarse_solver_type) != UNKNOWN);
#endif
        smoothError(error, residual, coarsest_ln, d_coarse_solver_max_iterations, false, false);
    }
    IBTK_TIMER_STOP(t_solve_coarsest_level);
    return true;
} // solveCoarsestLevel

void
CCPoissonPointRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
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
} // computeResidual

/////////////////////////////// PROTECTED ////////////////////////////////////

void
CCPoissonPointRelaxationFACOperator::initializeOperatorStateSpecialized(const SAMRAIVectorReal<NDIM, double>& solution,
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
        TBOX_ERROR("CCPoissonPointRelaxationFACOperator::initializeOperatorState()\n"
                   << "  solution and rhs vectors must have the same data depths\n"
                   << "  solution data depth = "
                   << solution_pdat_fac->getDefaultDepth()
                   << "\n"
                   << "  rhs      data depth = "
                   << rhs_pdat_fac->getDefaultDepth()
                   << std::endl);
    }

    const bool constant_coefficients =
        (d_poisson_spec.cIsZero() || d_poisson_spec.cIsConstant()) && d_poisson_spec.dIsConstant();
    if (!constant_coefficients)
    {
        TBOX_ERROR(d_object_name << "::initializeOperatorState():\n"
                                 << "  requires constant coefficients"
                                 << std::endl);
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellDataFactory<NDIM, double> > scratch_pdat_fac =
        var_db->getPatchDescriptor()->getPatchDataFactory(d_scratch_idx);
    scratch_pdat_fac->setDefaultDepth(solution_pdat_fac->getDefaultDepth());

    // Initialize the coarse level solvers when needed.
    if (coarsest_reset_ln == d_coarsest_ln && d_coarse_solver)
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
    d_op_stencil_fill_pattern = new CellNoCornersFillPattern(CELLG, true, false, false);

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

    // Get overlap information for re-setting patch boundary conditions during
    // smoothing.
    d_patch_neighbor_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_neighbor_overlap[ln].resize(num_local_patches);
        int patch_counter1 = 0;
        for (PatchLevel<NDIM>::Iterator p1(level); p1; p1++, ++patch_counter1)
        {
            d_patch_neighbor_overlap[ln][patch_counter1].clear();
            Pointer<Patch<NDIM> > dst_patch = level->getPatch(p1());
            const Box<NDIM>& dst_patch_box = dst_patch->getBox();
            const Box<NDIM>& dst_ghost_box = Box<NDIM>::grow(dst_patch_box, 1);
            int patch_counter2 = 0;
            for (PatchLevel<NDIM>::Iterator p2(level); patch_counter2 < patch_counter1; p2++, ++patch_counter2)
            {
                Pointer<Patch<NDIM> > src_patch = level->getPatch(p2());
                const Box<NDIM>& src_patch_box = src_patch->getBox();
                const Box<NDIM> overlap = dst_ghost_box * src_patch_box;
                if (!overlap.empty())
                {
                    d_patch_neighbor_overlap[ln][patch_counter1][p2()] = overlap;
                }
            }
        }
    }
    return;
} // initializeOperatorStateSpecialized

void
CCPoissonPointRelaxationFACOperator::deallocateOperatorStateSpecialized(const int /*coarsest_reset_ln*/,
                                                                        const int /*finest_reset_ln*/)
{
    if (!d_is_initialized) return;

    if (!d_in_initialize_operator_state)
    {
        d_patch_bc_box_overlap.clear();
        d_patch_neighbor_overlap.clear();
        if (d_coarse_solver) d_coarse_solver->deallocateSolverState();
    }
    return;
} // deallocateOperatorStateSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
