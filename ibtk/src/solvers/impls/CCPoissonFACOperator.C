// Filename: CCPoissonFACOperator.C
// Created on 10 Feb 2005 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "CCPoissonFACOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/CartCellDoubleCubicCoarsen.h>
#include <ibtk/CartCellDoubleQuadraticCFInterpolation.h>
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/NormOps.h>
#include <ibtk/RefinePatchStrategySet.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>

// BLITZ++ INCLUDES
#include <blitz/tinyvec.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <map>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define GS_SMOOTH_FC FC_FUNC(gssmooth2d,GSSMOOTH2D)
#define RB_GS_SMOOTH_FC FC_FUNC(rbgssmooth2d,RBGSSMOOTH2D)
#endif
#if (NDIM == 3)
#define GS_SMOOTH_FC FC_FUNC(gssmooth3d,GSSMOOTH3D)
#define RB_GS_SMOOTH_FC FC_FUNC(rbgssmooth3d,RBGSSMOOTH3D)
#endif

// Function interfaces
extern "C"
{
    void
    GS_SMOOTH_FC(
        double* U, const int& U_gcw,
        const double& alpha, const double& beta,
        const double* F, const int& F_gcw,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
#if (NDIM == 3)
        const int& ilower2, const int& iupper2,
#endif
        const double* dx, const int& sweeps);

    void
    RB_GS_SMOOTH_FC(
        double* U, const int& U_gcw,
        const double& alpha, const double& beta,
        const double* F, const int& F_gcw,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
#if (NDIM == 3)
        const int& ilower2, const int& iupper2,
#endif
        const double* dx, const int& sweeps);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_restrict_residual;
static Timer* t_prolong_error;
static Timer* t_prolong_error_and_correct;
static Timer* t_smooth_error;
static Timer* t_solve_coarsest_level;
static Timer* t_compute_residual;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values; used only to evaluate composite grid
// residuals.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries; used only to evaluate
// composite grid residuals.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells; used only to evaluate composite grid residuals.
static const bool CONSISTENT_TYPE_2_BDRY = false;

struct IndexComp
    : std::binary_function<Index<NDIM>,Index<NDIM>,bool>
{
    inline bool
    operator()(
        const Index<NDIM>& lhs,
        const Index<NDIM>& rhs) const
        {
            return ((lhs(0) < rhs(0))
#if (NDIM > 1)
                    || (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                    || (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
                    );
        }// operator()
};
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCPoissonFACOperator::CCPoissonFACOperator(
    const std::string& object_name,
    const Pointer<Database> input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_solution(NULL),
      d_rhs(NULL),
      d_depth(1),
      d_using_petsc_smoothers(true),
      d_gcw(CELLG),
      d_patch_vec_e(),
      d_patch_vec_f(),
      d_patch_mat(),
      d_patch_bc_box_overlap(),
      d_patch_smoother_bc_boxes(),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_hier_bdry_fill_ops(),
      d_hier_math_ops(),
      d_in_initialize_operator_state(false),
      d_coarsest_reset_ln(-1),
      d_finest_reset_ln(-1),
      d_poisson_spec(d_object_name+"::poisson_spec"),
      d_smoother_choice("additive"),
      d_prolongation_method("LINEAR_REFINE"),
      d_restriction_method("CONSERVATIVE_COARSEN"),
      d_preconditioner(NULL),
      d_coarse_solver_choice("block_jacobi"),
      d_coarse_solver_tol(1.0e-6),
      d_coarse_solver_max_its(10),
      d_using_hypre(d_coarse_solver_choice == "hypre"),
      d_hypre_solvers(d_depth),
      d_hypre_db(),
      d_using_petsc(d_coarse_solver_choice == "petsc"),
      d_petsc_solver(NULL),
      d_petsc_db(),
      d_context(NULL),
      d_cell_scratch_idx(-1),
      d_bc_op(NULL),
      d_default_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(
                            d_object_name+"::default_bc_coef", Pointer<Database>(NULL))),
      d_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(1,d_default_bc_coef)),
      d_apply_time(0.0),
      d_cf_bdry_op(),
      d_op_stencil_fill_pattern(),
      d_prolongation_refine_operator(),
      d_prolongation_refine_patch_strategy(),
      d_prolongation_refine_algorithm(),
      d_prolongation_refine_schedules(),
      d_restriction_coarsen_operator(),
      d_restriction_coarsen_algorithm(),
      d_restriction_coarsen_schedules(),
      d_ghostfill_nocoarse_refine_operator(),
      d_ghostfill_nocoarse_refine_algorithm(),
      d_ghostfill_nocoarse_refine_schedules()
{
    // Get values from the input database.
    if (!input_db.isNull())
    {
        d_smoother_choice = input_db->getStringWithDefault("smoother_choice", d_smoother_choice);
        d_prolongation_method = input_db->getStringWithDefault("prolongation_method", d_prolongation_method);
        d_restriction_method = input_db->getStringWithDefault("restriction_method", d_restriction_method);
        d_coarse_solver_choice = input_db->getStringWithDefault("coarse_solver_choice", d_coarse_solver_choice);
        d_coarse_solver_tol = input_db->getDoubleWithDefault("coarse_solver_tolerance", d_coarse_solver_tol);
        d_coarse_solver_max_its = input_db->getIntegerWithDefault("coarse_solver_max_iterations", d_coarse_solver_max_its);
        if (input_db->isDatabase("hypre_solver"))
        {
            d_hypre_db = input_db->getDatabase("hypre_solver");
        }
        if (input_db->isDatabase("petsc_solver"))
        {
            d_petsc_db = input_db->getDatabase("petsc_solver");
        }
    }
    sanityCheck();

    // Configure the bottom solver.
    setCoarsestLevelSolverChoice(d_coarse_solver_choice);

    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_bc_coef->setBoundaryValue(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setPhysicalBcCoef(d_default_bc_coef);

    // Setup scratch variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name+"::CONTEXT");

    const IntVector<NDIM> cell_ghosts = d_gcw;

    Pointer<CellVariable<NDIM,double> > cell_scratch_var =
        new CellVariable<NDIM,double>(d_object_name+"::cell_scratch", d_depth);
    if (var_db->checkVariableExists(cell_scratch_var->getName()))
    {
        cell_scratch_var = var_db->getVariable(cell_scratch_var->getName());
        d_cell_scratch_idx = var_db->mapVariableAndContextToIndex(cell_scratch_var, d_context);
        var_db->removePatchDataIndex(d_cell_scratch_idx);
    }
    d_cell_scratch_idx = var_db->registerVariableAndContext(cell_scratch_var, d_context, cell_ghosts);

    // Setup Timers.
    IBTK_DO_ONCE(
        t_restrict_residual         = TimerManager::getManager()->getTimer("IBTK::CCPoissonFACOperator::restrictResidual()");
        t_prolong_error             = TimerManager::getManager()->getTimer("IBTK::CCPoissonFACOperator::prolongError()");
        t_prolong_error_and_correct = TimerManager::getManager()->getTimer("IBTK::CCPoissonFACOperator::prolongErrorAndCorrect()");
        t_smooth_error              = TimerManager::getManager()->getTimer("IBTK::CCPoissonFACOperator::smoothError()");
        t_solve_coarsest_level      = TimerManager::getManager()->getTimer("IBTK::CCPoissonFACOperator::solveCoarsestLevel()");
        t_compute_residual          = TimerManager::getManager()->getTimer("IBTK::CCPoissonFACOperator::computeResidual()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBTK::CCPoissonFACOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBTK::CCPoissonFACOperator::deallocateOperatorState()");
                 );
    return;
}// CCPoissonFACOperator

CCPoissonFACOperator::~CCPoissonFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    delete d_default_bc_coef;
    return;
}// ~CCPoissonFACOperator

void
CCPoissonFACOperator::setPoissonSpecifications(
    const PoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
}// setPoissonSpecifications

void
CCPoissonFACOperator::setPhysicalBcCoef(
    RobinBcCoefStrategy<NDIM>* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(1,bc_coef));
    return;
}// setPhysicalBcCoef

void
CCPoissonFACOperator::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    d_bc_coefs.resize(bc_coefs.size());
    for (unsigned int l = 0; l < bc_coefs.size(); ++l)
    {
        if (bc_coefs[l] != NULL)
        {
            d_bc_coefs[l] = bc_coefs[l];
        }
        else
        {
            d_bc_coefs[l] = d_default_bc_coef;
        }
    }
    return;
}// setPhysicalBcCoefs

void
CCPoissonFACOperator::setPhysicalBcCoefs(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(&bc_coefs[0],&bc_coefs[0]+NDIM));
    return;
}// setPhysicalBcCoefs

void
CCPoissonFACOperator::setTime(
    const double time)
{
    d_apply_time = time;
    for (unsigned int k = 0; k < d_hypre_solvers.size(); ++k)
    {
        if (!d_hypre_solvers[k].isNull()) d_hypre_solvers[k]->setTime(d_apply_time);
    }
    if (!d_petsc_solver.isNull()) d_petsc_solver->setTime(d_apply_time);
    return;
}// setTime

void
CCPoissonFACOperator::setResetLevels(
    const int coarsest_ln,
    const int finest_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((coarsest_ln == -1 && finest_ln == -1) ||
                (coarsest_ln >=  0 && finest_ln >= coarsest_ln));
#endif
    if (d_is_initialized)
    {
        d_coarsest_reset_ln = coarsest_ln;
        d_finest_reset_ln = finest_ln;
    }
    return;
}// setResetLevels

void
CCPoissonFACOperator::setGhostCellWidth(
    const IntVector<NDIM>& ghost_cell_width)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setGhostCellWidth()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_gcw = ghost_cell_width;
    sanityCheck();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;
    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    Pointer<CellVariable<NDIM,double> > cell_scratch_var = var;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!cell_scratch_var.isNull());
#endif
    var_db->removePatchDataIndex(d_cell_scratch_idx);
    const IntVector<NDIM> cell_ghosts = d_gcw;
    d_cell_scratch_idx = var_db->registerVariableAndContext(cell_scratch_var, d_context, cell_ghosts);
#ifdef DEBUG_CHECK_ASSERTIONS
    Pointer<CellDataFactory<NDIM,double> > cell_scratch_pdat_fac = var_db->getPatchDescriptor()->getPatchDataFactory(d_cell_scratch_idx);
    TBOX_ASSERT(!cell_scratch_pdat_fac.isNull());
    TBOX_ASSERT(cell_scratch_pdat_fac->getGhostCellWidth() == d_gcw);
#endif
    return;
}// setGhostCellWidth

void
CCPoissonFACOperator::setSmootherChoice(
    const std::string& smoother_choice)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSmootherChoice()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_smoother_choice = smoother_choice;
    sanityCheck();
    return;
}// setSmootherChoice

void
CCPoissonFACOperator::setCoarsestLevelSolverChoice(
    const std::string& coarse_solver_choice)
{
    d_coarse_solver_choice = coarse_solver_choice;

    if (d_coarse_solver_choice == "hypre")
    {
        d_using_hypre = true;
        d_hypre_solvers.resize(d_depth);
        for (int depth = 0; depth < d_depth; ++depth)
        {
            std::ostringstream stream;
            stream << depth;
            d_hypre_solvers[depth] = new CCPoissonHypreLevelSolver(d_object_name+"::hypre_solver_"+stream.str(), d_hypre_db);
            d_hypre_solvers[depth]->setDataDepth(depth);
        }
        if (d_is_initialized)
        {
            initializeHypreLevelSolvers();
        }
    }
    else
    {
        d_using_hypre = false;
        d_hypre_solvers.clear();
    }

    if (d_coarse_solver_choice == "petsc")
    {
        d_using_petsc = true;
        if (d_is_initialized)
        {
            initializePETScLevelSolver();
        }
    }
    else
    {
        d_using_petsc = false;
        d_petsc_solver.setNull();
    }

    sanityCheck();
    return;
}// setCoarsestLevelSolverChoice

void
CCPoissonFACOperator::setCoarsestLevelSolverTolerance(
    double coarse_solver_tol)
{
    d_coarse_solver_tol = coarse_solver_tol;
    sanityCheck();
    return;
}// setCoarsestLevelSolverTolerance

void
CCPoissonFACOperator::setCoarsestLevelSolverMaxIterations(
    int coarse_solver_max_its)
{
    d_coarse_solver_max_its = coarse_solver_max_its;
    sanityCheck();
    return;
}// setCoarsestLevelSolverMaxIterations

void
CCPoissonFACOperator::setProlongationMethod(
    const std::string& prolongation_method)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setProlongationMethod()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_prolongation_method = prolongation_method;
    sanityCheck();
    return;
}// setProlongationMethod

void
CCPoissonFACOperator::setRestrictionMethod(
    const std::string& restriction_method)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setRestrictionMethod()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_restriction_method = restriction_method;
    sanityCheck();
    return;
}// setRestrictionMethod

///
///  The following routines:
///
///      setFACPreconditioner(),
///      restrictResidual(),
///      prolongError(),
///      prolongErrorAndCorrect(),
///      smoothError(),
///      solveCoarsestLevel(),
///      computeResidual(),
///      initializeOperatorState(),
///      deallocateOperatorState()
///
///  are concrete implementations of functions declared in the
///  FACPreconditionerStrategy abstract base class.
///

void
CCPoissonFACOperator::setFACPreconditioner(
    ConstPointer<FACPreconditioner> preconditioner)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setFACPreconditioner()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_preconditioner = preconditioner;
    sanityCheck();
    return;
}// setFACPreconditioner

void
CCPoissonFACOperator::restrictResidual(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    IBTK_TIMER_START(t_restrict_residual);

    const int src_idx = src.getComponentDescriptorIndex(0);
    const int dst_idx = dst.getComponentDescriptorIndex(0);

    if (src_idx != dst_idx)
    {
        HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(d_hierarchy, dst_ln, dst_ln);
        static const bool interior_only = false;
        hier_cc_data_ops.copyData(dst_idx, src_idx, interior_only);
    }
    xeqScheduleRestriction(dst_idx, src_idx, dst_ln);

    IBTK_TIMER_STOP(t_restrict_residual);
    return;
}// restrictResidual

void
CCPoissonFACOperator::prolongError(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    IBTK_TIMER_START(t_prolong_error);

    const int dst_idx = dst.getComponentDescriptorIndex(0);
    const int src_idx = src.getComponentDescriptorIndex(0);

    // Prolong the correction from the coarse src level data into the fine level
    // dst data.
    xeqScheduleProlongation(dst_idx, src_idx, dst_ln);

    IBTK_TIMER_STOP(t_prolong_error);
    return;
}// prolongError

void
CCPoissonFACOperator::prolongErrorAndCorrect(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    IBTK_TIMER_START(t_prolong_error_and_correct);

    const int dst_idx = dst.getComponentDescriptorIndex(0);
    const int src_idx = src.getComponentDescriptorIndex(0);

    // Prolong the correction from the coarse level src data into the fine level
    // scratch data and then correct the fine level dst data.
    static const bool interior_only = false;
    if (src_idx != dst_idx)
    {
        HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops_coarse(d_hierarchy, dst_ln-1, dst_ln-1);
        hier_cc_data_ops_coarse.add(dst_idx, dst_idx, src_idx, interior_only);
    }
    xeqScheduleProlongation(d_cell_scratch_idx, src_idx, dst_ln);
    HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops_fine(d_hierarchy, dst_ln, dst_ln);
    hier_cc_data_ops_fine.add(dst_idx, dst_idx, d_cell_scratch_idx, interior_only);

    IBTK_TIMER_STOP(t_prolong_error_and_correct);
    return;
}// prolongErrorAndCorrect

void
CCPoissonFACOperator::smoothError(
    SAMRAIVectorReal<NDIM,double>& error,
    const SAMRAIVectorReal<NDIM,double>& residual,
    int level_num,
    int num_sweeps,
    bool /*performing_pre_sweeps*/,
    bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBTK_TIMER_START(t_smooth_error);

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int error_idx = error.getComponentDescriptorIndex(0);
    const int scratch_idx = d_cell_scratch_idx;

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM,double> >   error_data = error.getComponentPatchData(0, *patch);
            Pointer<CellData<NDIM,double> > scratch_data = patch->getPatchData(scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
            TBOX_ASSERT(  error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
            scratch_data->getArrayData().copy(
                error_data->getArrayData(),
                d_patch_bc_box_overlap[level_num][patch_counter],
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
                    Pointer<CellData<NDIM,double> >   error_data = error.getComponentPatchData(0, *patch);
                    Pointer<CellData<NDIM,double> > scratch_data = patch->getPatchData(scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                    const Box<NDIM>& ghost_box = error_data->getGhostBox();
                    TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
                    TBOX_ASSERT(  error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    error_data->getArrayData().copy(
                        scratch_data->getArrayData(),
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
            Pointer<CellData<NDIM,double> >    error_data = error   .getComponentPatchData(0, *patch);
            Pointer<CellData<NDIM,double> > residual_data = residual.getComponentPatchData(0, *patch);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == residual_data->getGhostBox());
            TBOX_ASSERT(   error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(residual_data->getGhostCellWidth() == d_gcw);
#endif
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            // Copy updated values from other local patches.
            if (d_smoother_choice == "multiplicative")
            {
                const std::map<int,Box<NDIM> > smoother_bc_boxes = d_patch_smoother_bc_boxes[level_num][patch_counter];
                for (std::map<int,Box<NDIM> >::const_iterator cit = smoother_bc_boxes.begin();
                     cit != smoother_bc_boxes.end(); ++cit)
                {
                    const int src_patch_num = cit->first;
                    const Box<NDIM>& overlap = cit->second;
                    Pointer<Patch<NDIM> > src_patch = level->getPatch(src_patch_num);
                    Pointer<CellData<NDIM,double> > src_error_data = error.getComponentPatchData(0, *src_patch);
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
            if (d_using_petsc_smoothers)
            {
                // Reset ghost cell values in the residual data so that patch
                // boundary conditions are properly handled.
                residual_data->getArrayData().copy(
                    error_data->getArrayData(),
                    d_patch_bc_box_overlap[level_num][patch_counter],
                    IntVector<NDIM>(0));

                for (int depth = 0; depth < d_depth; ++depth)
                {
                    // Setup the PETSc Vec wrappers for the given patch data and
                    // data depth.
                    int ierr;

                    Vec& e = d_patch_vec_e[level_num][patch_counter];
                    Vec& f = d_patch_vec_f[level_num][patch_counter];

                    ierr = VecPlaceArray(e,    error_data->getPointer(depth));  IBTK_CHKERRQ(ierr);
                    ierr = VecPlaceArray(f, residual_data->getPointer(depth));  IBTK_CHKERRQ(ierr);

                    // Smooth the error on the patch using PETSc.  Here, we are
                    // approximately solving
                    //
                    //     Ae = f
                    //
                    // using an iteration of the form
                    //
                    //     e <- e + PC(f - Ae) = e + PC(r) = e + x.
                    //
                    // Presently, we simply employ symmetric Gauss-Seidel as the
                    // patch smoother.
                    static const double omega = 1.0;
                    static const double shift = 0.0;
                    static const int its = 1;
                    Mat& A = d_patch_mat[level_num][patch_counter];
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0)
                    ierr = MatRelax(A, f, omega, SOR_SYMMETRIC_SWEEP, shift, its, its, e);  IBTK_CHKERRQ(ierr);
#endif
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1)
                    ierr = MatSOR(A, f, omega, SOR_SYMMETRIC_SWEEP, shift, its, its, e);  IBTK_CHKERRQ(ierr);
#endif
                    // Reset the PETSc Vec wrappers.
                    ierr = VecResetArray(e);  IBTK_CHKERRQ(ierr);
                    ierr = VecResetArray(f);  IBTK_CHKERRQ(ierr);
                }
            }
            else
            {
                // Smooth the error via red-black Gauss-Seidel.
                const double& alpha = d_poisson_spec.getDConstant();
                const double& beta = d_poisson_spec.cIsZero() ? 0.0 : d_poisson_spec.getCConstant();
                for (int depth = 0; depth < d_depth; ++depth)
                {
                    double* const U = error_data->getPointer(depth);
                    const int U_ghosts = (error_data->getGhostCellWidth()).max();
                    const double* const F = residual_data->getPointer(depth);
                    const int F_ghosts = (residual_data->getGhostCellWidth()).max();
                    static const int its = 1;
                    RB_GS_SMOOTH_FC(
                        U, U_ghosts,
                        alpha, beta,
                        F, F_ghosts,
                        patch_box.lower(0), patch_box.upper(0),
                        patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                        patch_box.lower(2), patch_box.upper(2),
#endif
                        dx, its);
                }
            }
        }
    }

    IBTK_TIMER_STOP(t_smooth_error);
    return;
}// smoothError

bool
CCPoissonFACOperator::solveCoarsestLevel(
    SAMRAIVectorReal<NDIM,double>& error,
    const SAMRAIVectorReal<NDIM,double>& residual,
    int coarsest_ln)
{
    IBTK_TIMER_START(t_solve_coarsest_level);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
#endif
    if (!d_using_hypre && !d_using_petsc)
    {
        smoothError(error, residual, coarsest_ln, d_coarse_solver_max_its, false, false);
    }
    else
    {
        SAMRAIVectorReal<NDIM,double> error_level(error.getName()+"::level", error.getPatchHierarchy(), coarsest_ln, coarsest_ln);
        for (int comp = 0; comp < error.getNumberOfComponents(); ++comp)
        {
            error_level.addComponent(error.getComponentVariable(comp), error.getComponentDescriptorIndex(comp), error.getControlVolumeIndex(comp));
        }

        SAMRAIVectorReal<NDIM,double> residual_level(residual.getName()+"::level", residual.getPatchHierarchy(), coarsest_ln, coarsest_ln);
        for (int comp = 0; comp < residual.getNumberOfComponents(); ++comp)
        {
            residual_level.addComponent(residual.getComponentVariable(comp), residual.getComponentDescriptorIndex(comp), residual.getControlVolumeIndex(comp));
        }

        if (d_using_hypre)
        {
            for (int depth = 0; depth < d_depth; ++depth)
            {
                d_hypre_solvers[depth]->setInitialGuessNonzero(true);
                d_hypre_solvers[depth]->setMaxIterations(d_coarse_solver_max_its);
                d_hypre_solvers[depth]->setRelativeTolerance(d_coarse_solver_tol);
                d_hypre_solvers[depth]->solveSystem(error_level,residual_level);
            }
        }
        else if (d_using_petsc)
        {
            d_petsc_solver->setInitialGuessNonzero(true);
            d_petsc_solver->setMaxIterations(d_coarse_solver_max_its);
            d_petsc_solver->setRelativeTolerance(d_coarse_solver_tol);
            d_petsc_solver->solveSystem(error_level,residual_level);
        }
    }

    IBTK_TIMER_STOP(t_solve_coarsest_level);
    return true;
}// solveCoarsestLevel

void
CCPoissonFACOperator::computeResidual(
    SAMRAIVectorReal<NDIM,double>& residual,
    const SAMRAIVectorReal<NDIM,double>& solution,
    const SAMRAIVectorReal<NDIM,double>& rhs,
    int coarsest_level_num,
    int finest_level_num)
{
    IBTK_TIMER_START(t_compute_residual);

    const int res_idx = residual.getComponentDescriptorIndex(0);
    const int sol_idx = solution.getComponentDescriptorIndex(0);
    const int rhs_idx = rhs.getComponentDescriptorIndex(0);

    const Pointer<CellVariable<NDIM,double> > res_var = residual.getComponentVariable(0);
    const Pointer<CellVariable<NDIM,double> > sol_var = solution.getComponentVariable(0);
    const Pointer<CellVariable<NDIM,double> > rhs_var = rhs.getComponentVariable(0);

    // Fill ghost-cell values.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    Pointer<CellNoCornersFillPattern> fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    InterpolationTransactionComponent transaction_comp(sol_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_bc_coefs, fill_pattern);
    if (d_hier_bdry_fill_ops[finest_level_num].isNull())
    {
        d_hier_bdry_fill_ops[finest_level_num] = new HierarchyGhostCellInterpolation();
        d_hier_bdry_fill_ops[finest_level_num]->initializeOperatorState(transaction_comp, d_hierarchy, coarsest_level_num, finest_level_num);
    }
    else
    {
        d_hier_bdry_fill_ops[finest_level_num]->resetTransactionComponent(transaction_comp);
    }
    d_hier_bdry_fill_ops[finest_level_num]->setHomogeneousBc(true);
    d_hier_bdry_fill_ops[finest_level_num]->fillData(d_apply_time);

    // Compute the residual, r = f - A*u.
    if (d_hier_math_ops[finest_level_num].isNull())
    {
        std::ostringstream stream;
        stream << d_object_name << "::hier_math_ops_" << finest_level_num;
        d_hier_math_ops[finest_level_num] = new HierarchyMathOps(stream.str(), d_hierarchy, coarsest_level_num, finest_level_num);
    }
    d_hier_math_ops[finest_level_num]->laplace(res_idx, res_var, d_poisson_spec, sol_idx, sol_var, NULL, d_apply_time);
    HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    hier_cc_data_ops.axpy(res_idx, -1.0, res_idx, rhs_idx, false);

    IBTK_TIMER_STOP(t_compute_residual);
    return;
}// computeResidual

void
CCPoissonFACOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& solution,
    const SAMRAIVectorReal<NDIM,double>& rhs)
{
    IBTK_TIMER_START(t_initialize_operator_state);

    d_in_initialize_operator_state = true;

    // Cache the level range to be reset.
    //
    // NOTE: We cannot use d_coarsest_reset_ln and d_finest_reset_ln since those
    // values are reset by deallocateOperatorState().
    const int coarsest_reset_ln =
        (d_coarsest_reset_ln != -1 && d_finest_reset_ln != -1
         ? d_coarsest_reset_ln
         : solution.getCoarsestLevelNumber());
    const int finest_reset_ln =
        (d_coarsest_reset_ln != -1 && d_finest_reset_ln != -1
         ? d_finest_reset_ln
         : solution.getFinestLevelNumber());

    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_solution = solution.cloneVector(solution.getName());
    d_solution->allocateVectorData();

    d_rhs = rhs.cloneVector(rhs.getName());
    d_rhs->allocateVectorData();

    Pointer<CellVariable<NDIM,double> > solution_var = solution.getComponentVariable(0);
    Pointer<CellVariable<NDIM,double> >      rhs_var =      rhs.getComponentVariable(0);

    Pointer<CellDataFactory<NDIM,double> > solution_pdat_fac = solution_var->getPatchDataFactory();
    Pointer<CellDataFactory<NDIM,double> >      rhs_pdat_fac =      rhs_var->getPatchDataFactory();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!solution_var.isNull());
    TBOX_ASSERT(!rhs_var.isNull());
    TBOX_ASSERT(!solution_pdat_fac.isNull());
    TBOX_ASSERT(!rhs_pdat_fac.isNull());
#endif

    if (solution_pdat_fac->getDefaultDepth() != rhs_pdat_fac->getDefaultDepth())
    {
        TBOX_ERROR("CCPoissonFACOperator::initializeOperatorState()\n"
                   << "  solution and rhs vectors must have the same data depths\n"
                   << "  solution data depth = " << solution_pdat_fac->getDefaultDepth() << "\n"
                   << "  rhs      data depth = " << rhs_pdat_fac     ->getDefaultDepth() << std::endl);
    }

    const int old_depth = d_depth;
    d_depth = solution_pdat_fac->getDefaultDepth();

    if (d_depth != old_depth)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<CellDataFactory<NDIM,double> > cell_scratch_pdat_fac =
            var_db->getPatchDescriptor()->getPatchDataFactory(d_cell_scratch_idx);
        cell_scratch_pdat_fac->setDefaultDepth(d_depth);
    }

    // Reset the hierarchy configuration.
    d_hierarchy   = solution.getPatchHierarchy();
    d_coarsest_ln = solution.getCoarsestLevelNumber();
    d_finest_ln   = solution.getFinestLevelNumber();

    // Setup level operators.
    d_hier_bdry_fill_ops.resize(d_finest_ln+1, NULL);
    d_hier_math_ops.resize(d_finest_ln+1, NULL);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        d_hier_bdry_fill_ops[ln].setNull();
        d_hier_math_ops[ln].setNull();
    }

    // Allocate scratch data.
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_cell_scratch_idx)) level->allocatePatchData(d_cell_scratch_idx);
    }

    // Initialize the hypre solver when needed.
    if (d_using_hypre && (coarsest_reset_ln == d_coarsest_ln))
    {
        initializeHypreLevelSolvers();
    }

    // Initialize the petsc solver when needed.
    if (d_using_petsc && (coarsest_reset_ln == d_coarsest_ln))
    {
        initializePETScLevelSolver();
    }

    // Get the transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    IBTK_DO_ONCE(
        geometry->addSpatialCoarsenOperator(new CartCellDoubleCubicCoarsen());
                 );

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    d_prolongation_refine_operator = geometry->lookupRefineOperator(var, d_prolongation_method);
    d_cf_bdry_op = new CartCellDoubleQuadraticCFInterpolation();
    d_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_cf_bdry_op->setPatchDataIndex(d_cell_scratch_idx);
    d_cf_bdry_op->setPatchHierarchy(d_hierarchy);

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    d_restriction_coarsen_operator = geometry->lookupCoarsenOperator(var, d_restriction_method);
    d_ghostfill_nocoarse_refine_operator = NULL;

    // Make space for saving communication schedules.  There is no need to
    // delete the old schedules first because we have deallocated the solver
    // state above.
    d_bc_op = new CartCellRobinPhysBdryOp(d_cell_scratch_idx, d_bc_coefs, false);

    if (d_poisson_spec.dIsConstant())
    {
        d_op_stencil_fill_pattern = new CellNoCornersFillPattern(CELLG, true, false, false);
    }
    else
    {
        d_op_stencil_fill_pattern = NULL;
    }

    std::vector<RefinePatchStrategy<NDIM>*> prolongation_refine_patch_strategies;
    prolongation_refine_patch_strategies.push_back(d_cf_bdry_op);
    prolongation_refine_patch_strategies.push_back(d_bc_op);
    d_prolongation_refine_patch_strategy = new RefinePatchStrategySet(
        prolongation_refine_patch_strategies.begin(), prolongation_refine_patch_strategies.end(), false);

    d_prolongation_refine_schedules.resize(d_finest_ln+1);
    d_restriction_coarsen_schedules.resize(d_finest_ln);
    d_ghostfill_nocoarse_refine_schedules.resize(d_finest_ln+1);

    d_prolongation_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_restriction_coarsen_algorithm = new CoarsenAlgorithm<NDIM>();
    d_ghostfill_nocoarse_refine_algorithm = new RefineAlgorithm<NDIM>();

    d_prolongation_refine_algorithm->registerRefine(
        d_cell_scratch_idx,
        solution.getComponentDescriptorIndex(0),
        d_cell_scratch_idx,
        d_prolongation_refine_operator,
        d_op_stencil_fill_pattern);
    d_restriction_coarsen_algorithm->registerCoarsen(
        d_cell_scratch_idx,
        rhs.getComponentDescriptorIndex(0),
        d_restriction_coarsen_operator);
    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        d_ghostfill_nocoarse_refine_operator,
        d_op_stencil_fill_pattern);

    for (int dst_ln = d_coarsest_ln+1; dst_ln <= d_finest_ln; ++dst_ln)
    {
        d_prolongation_refine_schedules[dst_ln] =
            d_prolongation_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln),
                Pointer<PatchLevel<NDIM> >(),
                dst_ln-1, d_hierarchy, d_prolongation_refine_patch_strategy.getPointer());

        d_ghostfill_nocoarse_refine_schedules[dst_ln] =
            d_ghostfill_nocoarse_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln), d_bc_op.getPointer());
    }

    d_ghostfill_nocoarse_refine_schedules[d_coarsest_ln] =
        d_ghostfill_nocoarse_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_ln), d_bc_op.getPointer());

    for (int dst_ln = d_coarsest_ln; dst_ln < d_finest_ln; ++dst_ln)
    {
        d_restriction_coarsen_schedules[dst_ln] =
            d_restriction_coarsen_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln  ),
                d_hierarchy->getPatchLevel(dst_ln+1));
    }

    // Initialize all cached PETSc data.
    d_using_petsc_smoothers = !(d_poisson_spec.cIsZero() || d_poisson_spec.cIsConstant()) || !d_poisson_spec.dIsConstant();
    if (d_using_petsc_smoothers)
    {
        d_patch_vec_e.resize(d_finest_ln+1);
        d_patch_vec_f.resize(d_finest_ln+1);
        d_patch_mat.resize(d_finest_ln+1);
        for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
            d_patch_mat  [ln].resize(num_local_patches);
            d_patch_vec_e[ln].resize(num_local_patches);
            d_patch_vec_f[ln].resize(num_local_patches);

            int patch_counter = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Box<NDIM>& ghost_box = Box<NDIM>::grow(patch_box, d_gcw);
                const int size = ghost_box.size();

                int ierr;

                Vec& e = d_patch_vec_e[ln][patch_counter];
                Vec& f = d_patch_vec_f[ln][patch_counter];
                ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, size, PETSC_NULL, &e);  IBTK_CHKERRQ(ierr);
                ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, size, PETSC_NULL, &f);  IBTK_CHKERRQ(ierr);

                Mat& A = d_patch_mat[ln][patch_counter];
                buildPatchLaplaceOperator(A, d_poisson_spec, patch, d_gcw);
            }
        }
    }

    // Get overlap information for setting patch boundary conditions.
    d_patch_bc_box_overlap.resize(d_finest_ln+1);
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
    // multiplicative smoothing.
    if (d_smoother_choice == "multiplicative")
    {
        d_patch_smoother_bc_boxes.resize(d_finest_ln+1);
        for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
            d_patch_smoother_bc_boxes[ln].resize(num_local_patches);

            int patch_counter1 = 0;
            for (PatchLevel<NDIM>::Iterator p1(level); p1; p1++, ++patch_counter1)
            {
                d_patch_smoother_bc_boxes[ln][patch_counter1].clear();

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
                        d_patch_smoother_bc_boxes[ln][patch_counter1][p2()] = overlap;
                    }
                }
            }
        }
    }
    else
    {
        d_patch_smoother_bc_boxes.clear();
    }

    // Indicate that the operator is initialized.
    d_is_initialized = true;
    d_in_initialize_operator_state = false;

    IBTK_TIMER_STOP(t_initialize_operator_state);
    return;
}// initializeOperatorState

void
CCPoissonFACOperator::deallocateOperatorState()
{
    if (d_is_initialized && !d_in_initialize_operator_state &&
        (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
    {
        return;
    }

    IBTK_TIMER_START(t_deallocate_operator_state);

    if (d_is_initialized)
    {
        const int coarsest_reset_ln =
            (d_in_initialize_operator_state &&
             (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
            ? d_coarsest_reset_ln : d_coarsest_ln;
        const int finest_reset_ln =
            (d_in_initialize_operator_state &&
             (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
            ? d_finest_reset_ln : d_finest_ln;

        for (int ln = coarsest_reset_ln; ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
        {
            if (ln <= d_hierarchy->getFinestLevelNumber())
            {
                Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(d_cell_scratch_idx)) level->deallocatePatchData(d_cell_scratch_idx);
            }

            if (d_using_petsc_smoothers)
            {
                int ierr;

                for (std::vector<Vec>::iterator it = d_patch_vec_e[ln].begin();
                     it != d_patch_vec_e[ln].end(); ++it)
                {
                    Vec& e = *it;
                    ierr = VecDestroy(e);  IBTK_CHKERRQ(ierr);
                }
                d_patch_vec_e[ln].clear();
                for (std::vector<Vec>::iterator it = d_patch_vec_f[ln].begin();
                     it != d_patch_vec_f[ln].end(); ++it)
                {
                    Vec& f = *it;
                    ierr = VecDestroy(f);  IBTK_CHKERRQ(ierr);
                }
                d_patch_vec_f[ln].clear();
                for (std::vector<Mat>::iterator it = d_patch_mat[ln].begin();
                     it != d_patch_mat[ln].end(); ++it)
                {
                    Mat& A = *it;
                    ierr = MatDestroy(A);  IBTK_CHKERRQ(ierr);
                }
                d_patch_mat[ln].clear();
            }
        }

        if (d_using_hypre && (coarsest_reset_ln == d_coarsest_ln))
        {
            for (int depth = 0; depth < d_depth; ++depth)
            {
                d_hypre_solvers[depth]->deallocateSolverState();
            }
        }

        if (d_using_petsc && (coarsest_reset_ln == d_coarsest_ln))
        {
            d_petsc_solver->deallocateSolverState();
        }

        // Delete the solution and rhs vectors.
        d_solution->resetLevels(d_solution->getCoarsestLevelNumber(), std::min(d_solution->getFinestLevelNumber(),d_hierarchy->getFinestLevelNumber()));
        d_solution->freeVectorComponents();
        d_solution.setNull();

        d_rhs->resetLevels(d_rhs->getCoarsestLevelNumber(), std::min(d_rhs->getFinestLevelNumber(),d_hierarchy->getFinestLevelNumber()));
        d_rhs->freeVectorComponents();
        d_rhs.setNull();

        if (!d_in_initialize_operator_state ||
            (d_coarsest_reset_ln == -1) || (d_finest_reset_ln == -1))
        {
            d_patch_vec_e.resize(0);
            d_patch_vec_f.resize(0);
            d_patch_mat.resize(0);
            d_patch_bc_box_overlap.resize(0);
            d_patch_smoother_bc_boxes.resize(0);

            d_hierarchy.setNull();
            d_coarsest_ln = -1;
            d_finest_ln   = -1;

            d_hier_bdry_fill_ops.clear();
            d_hier_math_ops.clear();

            d_prolongation_refine_operator      .setNull();
            d_cf_bdry_op                        .setNull();

            d_prolongation_refine_patch_strategy.setNull();
            d_prolongation_refine_algorithm     .setNull();
            d_prolongation_refine_schedules     .resize(0);

            d_restriction_coarsen_operator .setNull();
            d_restriction_coarsen_algorithm.setNull();
            d_restriction_coarsen_schedules.resize(0);

            d_ghostfill_nocoarse_refine_operator .setNull();
            d_ghostfill_nocoarse_refine_algorithm.setNull();
            d_ghostfill_nocoarse_refine_schedules.resize(0);
        }
    }

    // Clear the "reset level" range.
    d_coarsest_reset_ln = -1;
    d_finest_reset_ln   = -1;

    // Indicate that the operator is not initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_operator_state);
    return;
}// deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CCPoissonFACOperator::xeqScheduleProlongation(
    const int dst_idx,
    const int src_idx,
    const int dst_ln)
{
    d_bc_op->setPatchDataIndex(dst_idx);
    d_bc_op->setPhysicalBcCoefs(d_bc_coefs);
    d_bc_op->setHomogeneousBc(true);
    d_cf_bdry_op->setPatchDataIndex(dst_idx);

    RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(dst_idx, src_idx, dst_idx, d_prolongation_refine_operator, d_op_stencil_fill_pattern);
    refiner.resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    d_prolongation_refine_schedules[dst_ln]->fillData(d_apply_time);
    d_prolongation_refine_algorithm->resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    return;
}// xeqScheduleProlongation

void
CCPoissonFACOperator::xeqScheduleRestriction(
    const int dst_idx,
    const int src_idx,
    const int dst_ln)
{
    CoarsenAlgorithm<NDIM> coarsener;
    coarsener.registerCoarsen(dst_idx, src_idx, d_restriction_coarsen_operator);
    coarsener.resetSchedule(d_restriction_coarsen_schedules[dst_ln]);
    d_restriction_coarsen_schedules[dst_ln]->coarsenData();
    d_restriction_coarsen_algorithm->resetSchedule(d_restriction_coarsen_schedules[dst_ln]);
    return;
}// xeqScheduleRestriction

void
CCPoissonFACOperator::xeqScheduleGhostFillNoCoarse(
    const int dst_idx,
    const int dst_ln)
{
    d_bc_op->setPatchDataIndex(dst_idx);
    d_bc_op->setPhysicalBcCoefs(d_bc_coefs);
    d_bc_op->setHomogeneousBc(true);

    RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(dst_idx, dst_idx, dst_idx, d_ghostfill_nocoarse_refine_operator, d_op_stencil_fill_pattern);
    refiner.resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    d_ghostfill_nocoarse_refine_schedules[dst_ln]->fillData(d_apply_time);
    d_ghostfill_nocoarse_refine_algorithm->resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    return;
}// xeqScheduleGhostFillNoCoarse

void
CCPoissonFACOperator::initializeHypreLevelSolvers()
{
    SAMRAIVectorReal<NDIM,double> solution_level(d_solution->getName()+"::level", d_solution->getPatchHierarchy(), d_coarsest_ln, d_coarsest_ln);
    for (int comp = 0; comp < d_solution->getNumberOfComponents(); ++comp)
    {
        solution_level.addComponent(d_solution->getComponentVariable(comp), d_solution->getComponentDescriptorIndex(comp), d_solution->getControlVolumeIndex(comp));
    }

    SAMRAIVectorReal<NDIM,double> rhs_level(d_rhs->getName()+"::level", d_rhs->getPatchHierarchy(), d_coarsest_ln, d_coarsest_ln);
    for (int comp = 0; comp < d_rhs->getNumberOfComponents(); ++comp)
    {
        rhs_level.addComponent(d_rhs->getComponentVariable(comp), d_rhs->getComponentDescriptorIndex(comp), d_rhs->getControlVolumeIndex(comp));
    }

    // Note that since the bottom solver is solving for the error, it must
    // always employ homogeneous boundary conditions.
    if (d_hypre_solvers.size() != static_cast<unsigned int>(d_depth))
    {
        d_hypre_solvers.resize(d_depth);
        for (int depth = 0; depth < d_depth; ++depth)
        {
            std::ostringstream stream;
            stream << depth;
            d_hypre_solvers[depth] = new CCPoissonHypreLevelSolver(d_object_name+"::hypre_solver_"+stream.str(), d_hypre_db);
            d_hypre_solvers[depth]->setDataDepth(depth);
        }
    }
    for (int depth = 0; depth < d_depth; ++depth)
    {
        d_hypre_solvers[depth]->setPoissonSpecifications(d_poisson_spec);
        d_hypre_solvers[depth]->setPhysicalBcCoef(d_bc_coefs[depth]);
        d_hypre_solvers[depth]->setTime(d_apply_time);
        d_hypre_solvers[depth]->setHomogeneousBc(true);
        d_hypre_solvers[depth]->initializeSolverState(solution_level, rhs_level);
    }
    return;
}// initializeHypreLevelSolvers

void
CCPoissonFACOperator::initializePETScLevelSolver()
{
    d_petsc_solver = new CCPoissonPETScLevelSolver(d_object_name+"::petsc_solver", d_petsc_db);
    d_petsc_solver->setTime(d_apply_time);

    SAMRAIVectorReal<NDIM,double> solution_level(d_solution->getName()+"::level", d_solution->getPatchHierarchy(), d_coarsest_ln, d_coarsest_ln);
    for (int comp = 0; comp < d_solution->getNumberOfComponents(); ++comp)
    {
        solution_level.addComponent(d_solution->getComponentVariable(comp), d_solution->getComponentDescriptorIndex(comp), d_solution->getControlVolumeIndex(comp));
    }

    SAMRAIVectorReal<NDIM,double> rhs_level(d_rhs->getName()+"::level", d_rhs->getPatchHierarchy(), d_coarsest_ln, d_coarsest_ln);
    for (int comp = 0; comp < d_rhs->getNumberOfComponents(); ++comp)
    {
        rhs_level.addComponent(d_rhs->getComponentVariable(comp), d_rhs->getComponentDescriptorIndex(comp), d_rhs->getControlVolumeIndex(comp));
    }

    // Note that since the bottom solver is solving for the error, it must
    // always employ homogeneous boundary conditions.
    for (int depth = 0; depth < d_depth; ++depth)
    {
        d_petsc_solver->setPoissonSpecifications(d_poisson_spec);
        d_petsc_solver->setPhysicalBcCoefs(d_bc_coefs);
        d_petsc_solver->setHomogeneousBc(true);
        d_petsc_solver->initializeSolverState(solution_level, rhs_level);
    }
    return;
}// initializePETScLevelSolver

void
CCPoissonFACOperator::buildPatchLaplaceOperator(
    Mat& A,
    const PoissonSpecifications& poisson_spec,
    const Pointer<Patch<NDIM> > patch,
    const IntVector<NDIM>& ghost_cell_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    if (ghost_cell_width.min() == 0)
    {
        TBOX_ERROR("CCPoissonFACOperator::buildPatchLaplaceOperator():\n"
                   << "  ghost cells are required in all directions" << std::endl);
    }
#endif

    // Get the Poisson problem coefficients.
    const Box<NDIM>& patch_box = patch->getBox();
    static const IntVector<NDIM> no_ghosts = 0;

    Pointer<CellData<NDIM,double> > C_data;
    if (!poisson_spec.cIsZero() && !poisson_spec.cIsConstant())
    {
        C_data = patch->getPatchData(poisson_spec.getCPatchDataId());
        if (C_data.isNull())
        {
            TBOX_ERROR("CCPoissonFACOperator::buildPatchLaplaceOperator()\n"
                       << "  to solve (C u + div D grad u) = f with non-constant C,\n"
                       << "  C must be cell-centered double precision data" << std::endl);
        }
    }
    else
    {
        C_data = new CellData<NDIM,double>(patch_box, 1, no_ghosts);
        if (poisson_spec.cIsZero()) C_data->fill(0.0);
        else C_data->fill(poisson_spec.getCConstant());
    }

    Pointer<SideData<NDIM,double> > D_data;
    if (!poisson_spec.dIsConstant())
    {
        D_data = patch->getPatchData(poisson_spec.getDPatchDataId());
        if (D_data.isNull())
        {
            TBOX_ERROR("CCPoissonFACOperator::buildPatchLaplaceOperator()\n"
                       << "  to solve C u + div D grad u = f with non-constant D,\n"
                       << "  D must be side-centered double precision data" << std::endl);
        }
    }
    else
    {
        D_data = new SideData<NDIM,double>(patch_box, 1, no_ghosts);
        D_data->fill(poisson_spec.getDConstant());
    }

    // Build the patch operator.
    if (D_data->getDepth() == 1)
    {
        // Isotropic diffusion or grid aligned anisotropy.
        buildPatchLaplaceOperator_aligned(A, C_data, D_data, patch, ghost_cell_width);
    }
    else if (D_data->getDepth() == NDIM)
    {
        // Non-grid aligned anisotropy.
        buildPatchLaplaceOperator_nonaligned(A, C_data, D_data, patch, ghost_cell_width);
    }
    else
    {
        TBOX_ERROR("CCPoissonFACOperator::buildPatchLaplaceOperator()\n"
                   << "  D must be side-centered patch data with either 1 or NDIM components" << std::endl);
    }
    return;
}// buildPatchLaplaceOperator

void
CCPoissonFACOperator::buildPatchLaplaceOperator_aligned(
    Mat& A,
    const Pointer<CellData<NDIM,double> > C_data,
    const Pointer<SideData<NDIM,double> > D_data,
    const Pointer<Patch<NDIM> > patch,
    const IntVector<NDIM>& ghost_cell_width)
{
    int ierr;

    // Allocate a PETSc matrix for the patch operator.
    const Box<NDIM>& patch_box = patch->getBox();
    const Box<NDIM>& ghost_box = Box<NDIM>::grow(patch_box, ghost_cell_width);
    const int size = ghost_box.size();

    static const int stencil_sz = 2*NDIM+1;

    BoxList<NDIM> ghost_boxes(ghost_box);
    ghost_boxes.removeIntersections(patch_box);
    std::vector<int> nnz(size, stencil_sz);
    for (BoxList<NDIM>::Iterator bl(ghost_boxes); bl; bl++)
    {
        for (Box<NDIM>::Iterator b(bl()); b; b++)
        {
            nnz[ghost_box.offset(b())] = 1;
        }
    }
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, PETSC_DEFAULT, &nnz[0], &A);  IBTK_CHKERRQ(ierr);

    // Set some general matrix options.
#ifdef DEBUG_CHECK_ASSERTIONS
    ierr = MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);  IBTK_CHKERRQ(ierr);
#endif

    // Setup the finite difference stencil.  The stencil order is chosen to
    // optimize performance when setting the matrix coefficients.
    blitz::TinyVector<int,NDIM> num_cells;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        num_cells[d] = ghost_box.numberCells(d);
    }
    std::vector<int> mat_stencil(stencil_sz);
#if (NDIM == 2)
    static const int x_axis = 0;
    mat_stencil[0] = -num_cells[x_axis]; // ylower
    mat_stencil[1] = -1;                 // xlower
    mat_stencil[2] = 0;
    mat_stencil[3] = +1;                 // xupper
    mat_stencil[4] = +num_cells[x_axis]; // yupper
#endif
#if (NDIM == 3)
    static const int x_axis = 0;
    static const int y_axis = 1;
    mat_stencil[0] = -num_cells[x_axis]*num_cells[y_axis]; // zlower
    mat_stencil[1] = -num_cells[x_axis];                   // ylower
    mat_stencil[2] = -1;                                   // xlower
    mat_stencil[3] = 0;
    mat_stencil[4] = +1;                                   // xupper
    mat_stencil[5] = +num_cells[x_axis];                   // yupper
    mat_stencil[6] = +num_cells[x_axis]*num_cells[y_axis]; // zupper
#endif

    // Set the matrix coefficients to correspond to the standard finite
    // difference stencil for the Laplace operator.
    //
    // Note that boundary conditions at both physical boundaries and at
    // coarse-fine interfaces are implicitly treated by setting ghost cell
    // values appropriately.  Thus the matrix coefficients are independent of
    // any boundary conditions.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    for (Box<NDIM>::Iterator b(patch_box); b; b++)
    {
        const Index<NDIM>& i = b();

        std::vector<double> mat_vals(stencil_sz,0.0);
        mat_vals[NDIM] = (*C_data)(i);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const double& h = dx[axis];
            {
                const SideIndex<NDIM> ilower(i, axis, SideIndex<NDIM>::Lower);
                const double& D_lower = (*D_data)(ilower);
                mat_vals[NDIM-axis-1] += D_lower/(h*h);
                mat_vals[NDIM       ] -= D_lower/(h*h);
            }
            {
                const SideIndex<NDIM> iupper(i, axis, SideIndex<NDIM>::Upper);
                const double& D_upper = (*D_data)(iupper);
                mat_vals[NDIM+axis+1] += D_upper/(h*h);
                mat_vals[NDIM       ] -= D_upper/(h*h);
            }
        }

        static const int m = 1;
        static const int n = stencil_sz;
        std::vector<int> idxn(stencil_sz);
        const int idxm = ghost_box.offset(i);

        std::transform(mat_stencil.begin(), mat_stencil.end(),
                       idxn.begin(), std::bind2nd(std::plus<int>(), idxm));
        ierr = MatSetValues(A, m, &idxm, n, &idxn[0], &mat_vals[0],
                            INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Set the entries in the ghost cell region so that ghost cell values are
    // not modified by the smoother.
    for (BoxList<NDIM>::Iterator bl(ghost_boxes); bl; bl++)
    {
        for (Box<NDIM>::Iterator b(bl()); b; b++)
        {
            const int i = ghost_box.offset(b());
            ierr = MatSetValue(A, i, i, 1.0, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    return;
}// buildPatchLaplaceOperator_aligned

void
CCPoissonFACOperator::buildPatchLaplaceOperator_nonaligned(
    Mat& A,
    const Pointer<CellData<NDIM,double> > C_data,
    const Pointer<SideData<NDIM,double> > D_data,
    const Pointer<Patch<NDIM> > patch,
    const IntVector<NDIM>& ghost_cell_width)
{
    int ierr;

    // Allocate a PETSc matrix for the patch operator.
    const Box<NDIM>& patch_box = patch->getBox();
    const Box<NDIM>& ghost_box = Box<NDIM>::grow(patch_box, ghost_cell_width);
    const int size = ghost_box.size();

    static const int stencil_sz = (NDIM == 2 ? 9 : 19);

    BoxList<NDIM> ghost_boxes(ghost_box);
    ghost_boxes.removeIntersections(patch_box);
    std::vector<int> nnz(size, stencil_sz);
    for (BoxList<NDIM>::Iterator bl(ghost_boxes); bl; bl++)
    {
        for (Box<NDIM>::Iterator b(bl()); b; b++)
        {
            nnz[ghost_box.offset(b())] = 1;
        }
    }
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, PETSC_DEFAULT, &nnz[0], &A);  IBTK_CHKERRQ(ierr);

    // Set some general matrix options.
#ifdef DEBUG_CHECK_ASSERTIONS
    ierr = MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);  IBTK_CHKERRQ(ierr);
#endif

    // Setup the finite difference stencil.
    blitz::TinyVector<int,NDIM> num_cells;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        num_cells[d] = ghost_box.numberCells(d);
    }
    std::vector<int> mat_stencil(stencil_sz);
    std::map<Index<NDIM>,int,IndexComp> stencil_indices;
    int stencil_index = 0;
    static const int x_axis = 0;
#if (NDIM == 3)
    static const int y_axis = 1;
    for (int z_offset = -1; z_offset <= 1; ++z_offset)
#endif
    {
        for (int y_offset = -1; y_offset <= 1; ++y_offset)
        {
            for (int x_offset = -1; x_offset <= 1; ++x_offset)
            {
#if (NDIM == 3)
                // No full-corner coupling in 3D.
                if (x_offset == 0 || y_offset == 0 || z_offset == 0)
                {
#endif
#if (NDIM == 2)
                    const Index<NDIM> i(x_offset,y_offset);
#endif
#if (NDIM == 3)
                    const Index<NDIM> i(x_offset,y_offset,z_offset);
#endif
                    stencil_indices[i] = stencil_index++;
                    mat_stencil[stencil_indices[i]] = x_offset + y_offset*num_cells[x_axis]
#if (NDIM == 3)
                        + z_offset*num_cells[x_axis]*num_cells[y_axis]
#endif
                        ;
#if (NDIM == 3)
                }
#endif
            }
        }
    }

    // Set the matrix coefficients to correspond to a second-order accurate
    // finite difference stencil for the Laplace operator.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    for (Box<NDIM>::Iterator b(patch_box); b; b++)
    {
        const Index<NDIM>& i = b();
        static const Index<NDIM> i_stencil_center(0);
        const int stencil_center = stencil_indices[i_stencil_center];

        std::vector<double> mat_vals(stencil_sz,0.0);
        mat_vals[stencil_center] = (*C_data)(i);

        // The grid aligned part of the stencil (normal derivatives).
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const double& h = dx[axis];

            // Lower side normal flux.
            {
                Index<NDIM> i_stencil_lower(0);
                --i_stencil_lower[axis];
                const int stencil_lower = stencil_indices[i_stencil_lower];

                const SideIndex<NDIM> ilower(i, axis, SideIndex<NDIM>::Lower);
                const double& D_lower = (*D_data)(ilower,axis);
                mat_vals[stencil_lower ] += D_lower/(h*h);
                mat_vals[stencil_center] -= D_lower/(h*h);
            }

            // Upper side normal flux.
            {
                Index<NDIM> i_stencil_upper(0);
                ++i_stencil_upper[axis];
                const int stencil_upper = stencil_indices[i_stencil_upper];

                const SideIndex<NDIM> iupper(i, axis, SideIndex<NDIM>::Upper);
                const double& D_upper = (*D_data)(iupper,axis);
                mat_vals[stencil_upper ] += D_upper/(h*h);
                mat_vals[stencil_center] -= D_upper/(h*h);
            }
        }

        // The non-grid aligned part of the stencil (transverse derivatives).
        for (unsigned int norm_axis = 0; norm_axis < NDIM; ++norm_axis)
        {
            const double& norm_h = dx[norm_axis];
            for (unsigned int trans_axis = 0; trans_axis < NDIM; ++trans_axis)
            {
                if (norm_axis == trans_axis) break;
                const double& trans_h = dx[trans_axis];

                // Lower side transverse flux.
                {
                    const SideIndex<NDIM> ilower(i, norm_axis, SideIndex<NDIM>::Lower);
                    for (int norm_shift = -1; norm_shift <= 0; ++norm_shift)
                    {
                        for (int trans_shift = -1; trans_shift <= 1; trans_shift += 2)
                        {
                            Index<NDIM> i_stencil(0);
                            i_stencil[ norm_axis] +=  norm_shift;
                            i_stencil[trans_axis] += trans_shift;
                            const int stencil_index = stencil_indices[i_stencil];
                            if (trans_shift == 1)
                            {
                                mat_vals[stencil_index] -= 0.25*(*D_data)(ilower,norm_axis)/(norm_h*trans_h);
                            }
                            else
                            {
                                mat_vals[stencil_index] += 0.25*(*D_data)(ilower,norm_axis)/(norm_h*trans_h);
                            }
                        }
                    }
                }

                // Upper side transverse flux.
                {
                    const SideIndex<NDIM> iupper(i, norm_axis, SideIndex<NDIM>::Upper);
                    for (int norm_shift = 0; norm_shift <= 1; ++norm_shift)
                    {
                        for (int trans_shift = -1; trans_shift <= 1; trans_shift += 2)
                        {
                            Index<NDIM> i_stencil(0);
                            i_stencil[ norm_axis] +=  norm_shift;
                            i_stencil[trans_axis] += trans_shift;
                            const int stencil_index = stencil_indices[i_stencil];
                            if (trans_shift == 1)
                            {
                                mat_vals[stencil_index] += 0.25*(*D_data)(iupper,norm_axis)/(norm_h*trans_h);
                            }
                            else
                            {
                                mat_vals[stencil_index] -= 0.25*(*D_data)(iupper,norm_axis)/(norm_h*trans_h);
                            }
                        }
                    }
                }
            }
        }

        static const int m = 1;
        static const int n = stencil_sz;
        std::vector<int> idxn(stencil_sz);
        const int idxm = ghost_box.offset(i);

        std::transform(mat_stencil.begin(), mat_stencil.end(),
                       idxn.begin(), std::bind2nd(std::plus<int>(), idxm));
        ierr = MatSetValues(A, m, &idxm, n, &idxn[0], &mat_vals[0],
                            INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Set the entries in the ghost cell region so that ghost cell values are
    // not modified by the smoother.
    for (BoxList<NDIM>::Iterator bl(ghost_boxes); bl; bl++)
    {
        for (Box<NDIM>::Iterator b(bl()); b; b++)
        {
            const int i = ghost_box.offset(b());
            ierr = MatSetValue(A, i, i, 1.0, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    return;
}// buildPatchLaplaceOperator_nonaligned

void
CCPoissonFACOperator::sanityCheck()
{
    if (d_gcw.min() <= 0)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  ghost_cell_width.min() must be greater than zero" << std::endl);
    }

    if (d_smoother_choice != "additive" &&
        d_smoother_choice != "multiplicative")
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  unknown smoother type: " << d_smoother_choice << "\n"
                   << "  valid choices are: additive, multiplicative" << std::endl);
    }

    if (d_coarse_solver_choice != "block_jacobi" &&
        d_coarse_solver_choice != "hypre" &&
        d_coarse_solver_choice != "petsc")
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  unknown coarse solver type: " << d_coarse_solver_choice << "\n"
                   << "  valid choices are: block_jacobi, hypre, petsc" << std::endl);
    }

    if (d_coarse_solver_tol < 0.0)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  invalid coarse solver tolerance: " << d_coarse_solver_tol << std::endl);
    }

    if (d_coarse_solver_max_its <= 0)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  invalid coarse solver maximum iterations: " << d_coarse_solver_max_its << std::endl);
    }

    for (unsigned int l = 0; l < d_bc_coefs.size(); ++l)
    {
        if (d_bc_coefs[l] == NULL)
        {
            TBOX_ERROR(d_object_name << ":\n"
                       << "  invalid physical bc object at depth = " << l << std::endl);
        }
    }
    return;
}// sanityCheck

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
