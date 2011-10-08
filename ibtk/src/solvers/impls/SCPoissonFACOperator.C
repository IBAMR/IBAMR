// Filename: SCPoissonFACOperator.C
// Created on 13 Nov 2008 by Boyce Griffith
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

#include "SCPoissonFACOperator.h"

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
#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/CartSideDoubleQuadraticCFInterpolation.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/NormOps.h>
#include <ibtk/RefinePatchStrategySet.h>
#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/SideSynchCopyFillPattern.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>

// C++ STDLIB INCLUDES
#include <algorithm>

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
static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

SCPoissonFACOperator::SCPoissonFACOperator(
    const std::string& object_name,
    const Pointer<Database> input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_solution(NULL),
      d_rhs(NULL),
      d_depth(1),
      d_using_petsc_smoothers(true),
      d_gcw(SIDEG),
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
      d_prolongation_method("CONSTANT_REFINE"),
      d_restriction_method("CONSERVATIVE_COARSEN"),
      d_preconditioner(NULL),
      d_coarse_solver_choice("block_jacobi"),
      d_coarse_solver_tol(1.0e-6),
      d_coarse_solver_max_its(10),
      d_using_hypre(d_coarse_solver_choice == "hypre"),
      d_hypre_solver(NULL),
      d_hypre_db(),
      d_using_petsc(d_coarse_solver_choice == "petsc"),
      d_petsc_solver(NULL),
      d_petsc_db(),
      d_context(NULL),
      d_side_scratch_idx(-1),
      d_bc_op(NULL),
      d_default_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_bc_coef", Pointer<Database>(NULL))),
      d_bc_coefs(blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>(d_default_bc_coef)),
      d_apply_time(0.0),
      d_cf_bdry_op(),
      d_op_stencil_fill_pattern(),
      d_side_synch_fill_pattern(),
      d_prolongation_refine_operator(),
      d_prolongation_refine_patch_strategy(),
      d_prolongation_refine_algorithm(),
      d_prolongation_refine_schedules(),
      d_restriction_coarsen_operator(),
      d_restriction_coarsen_algorithm(),
      d_restriction_coarsen_schedules(),
      d_ghostfill_nocoarse_refine_operator(),
      d_ghostfill_nocoarse_refine_algorithm(),
      d_ghostfill_nocoarse_refine_schedules(),
      d_side_synch_refine_operator(),
      d_side_synch_refine_algorithm(),
      d_side_synch_refine_schedules()
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

    // Create the hypre solver, if needed.
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
    setPhysicalBcCoefs(blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>(d_default_bc_coef));

    // Setup scratch variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name+"::CONTEXT");

    const IntVector<NDIM> side_ghosts = d_gcw;

    Pointer<SideVariable<NDIM,double> > side_scratch_var =
        new SideVariable<NDIM,double>(d_object_name+"::side_scratch", d_depth, true);
    if (var_db->checkVariableExists(side_scratch_var->getName()))
    {
        side_scratch_var = var_db->getVariable(side_scratch_var->getName());
        d_side_scratch_idx = var_db->mapVariableAndContextToIndex(side_scratch_var, d_context);
        var_db->removePatchDataIndex(d_side_scratch_idx);
    }
    d_side_scratch_idx = var_db->registerVariableAndContext(side_scratch_var, d_context, side_ghosts);

    // Setup Timers.
    IBTK_DO_ONCE(
        t_restrict_residual         = TimerManager::getManager()->getTimer("IBTK::SCPoissonFACOperator::restrictResidual()");
        t_prolong_error             = TimerManager::getManager()->getTimer("IBTK::SCPoissonFACOperator::prolongError()");
        t_prolong_error_and_correct = TimerManager::getManager()->getTimer("IBTK::SCPoissonFACOperator::prolongErrorAndCorrect()");
        t_smooth_error              = TimerManager::getManager()->getTimer("IBTK::SCPoissonFACOperator::smoothError()");
        t_solve_coarsest_level      = TimerManager::getManager()->getTimer("IBTK::SCPoissonFACOperator::solveCoarsestLevel()");
        t_compute_residual          = TimerManager::getManager()->getTimer("IBTK::SCPoissonFACOperator::computeResidual()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBTK::SCPoissonFACOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBTK::SCPoissonFACOperator::deallocateOperatorState()");
                 );
    return;
}// SCPoissonFACOperator

SCPoissonFACOperator::~SCPoissonFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    delete d_default_bc_coef;
    return;
}// ~SCPoissonFACOperator

void
SCPoissonFACOperator::setPoissonSpecifications(
    const PoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
}// setPoissonSpecifications

void
SCPoissonFACOperator::setPhysicalBcCoefs(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs)
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (bc_coefs[d] != NULL)
        {
            d_bc_coefs[d] = bc_coefs[d];
        }
        else
        {
            d_bc_coefs[d] = d_default_bc_coef;
        }
    }
    return;
}// setPhysicalBcCoefs

void
SCPoissonFACOperator::setTime(
    const double time)
{
    d_apply_time = time;
    if (!d_hypre_solver.isNull()) d_hypre_solver->setTime(d_apply_time);
    if (!d_petsc_solver.isNull()) d_petsc_solver->setTime(d_apply_time);
    return;
}// setTime

void
SCPoissonFACOperator::setResetLevels(
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
SCPoissonFACOperator::setGhostCellWidth(
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
    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    Pointer<SideVariable<NDIM,double> > side_scratch_var = var;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!side_scratch_var.isNull());
#endif
    var_db->removePatchDataIndex(d_side_scratch_idx);
    const IntVector<NDIM> side_ghosts = d_gcw;
    d_side_scratch_idx = var_db->registerVariableAndContext(side_scratch_var, d_context, side_ghosts);
#ifdef DEBUG_CHECK_ASSERTIONS
    Pointer<SideDataFactory<NDIM,double> > side_scratch_pdat_fac = var_db->getPatchDescriptor()->getPatchDataFactory(d_side_scratch_idx);
    TBOX_ASSERT(!side_scratch_pdat_fac.isNull());
    TBOX_ASSERT(side_scratch_pdat_fac->getGhostCellWidth() == d_gcw);
#endif
    return;
}// setGhostCellWidth

void
SCPoissonFACOperator::setSmootherChoice(
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
SCPoissonFACOperator::setCoarsestLevelSolverChoice(
    const std::string& coarse_solver_choice)
{
    d_coarse_solver_choice = coarse_solver_choice;

    if (d_coarse_solver_choice == "hypre")
    {
        d_using_hypre = true;
        d_hypre_solver = new SCPoissonHypreLevelSolver(d_object_name+"::hypre_solver", d_hypre_db);
        if (d_is_initialized)
        {
            initializeHypreLevelSolver();
        }
    }
    else
    {
        d_using_hypre = false;
        d_hypre_solver.setNull();
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
SCPoissonFACOperator::setCoarsestLevelSolverTolerance(
    double coarse_solver_tol)
{
    d_coarse_solver_tol = coarse_solver_tol;
    sanityCheck();
    return;
}// setCoarsestLevelSolverTolerance

void
SCPoissonFACOperator::setCoarsestLevelSolverMaxIterations(
    int coarse_solver_max_its)
{
    d_coarse_solver_max_its = coarse_solver_max_its;
    sanityCheck();
    return;
}// setCoarsestLevelSolverMaxIterations

void
SCPoissonFACOperator::setProlongationMethod(
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
SCPoissonFACOperator::setRestrictionMethod(
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
///  FACOperatorStrategy abstract base class.
///

void
SCPoissonFACOperator::setFACPreconditioner(
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
SCPoissonFACOperator::restrictResidual(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    IBTK_TIMER_START(t_restrict_residual);

    const int src_idx = src.getComponentDescriptorIndex(0);
    const int dst_idx = dst.getComponentDescriptorIndex(0);

    if (src_idx != dst_idx)
    {
        HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops(d_hierarchy, dst_ln, dst_ln);
        static const bool interior_only = false;
        hier_sc_data_ops.copyData(dst_idx, src_idx, interior_only);
    }
    xeqScheduleRestriction(dst_idx, src_idx, dst_ln);

    IBTK_TIMER_STOP(t_restrict_residual);
    return;
}// restrictResidual

void
SCPoissonFACOperator::prolongError(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    IBTK_TIMER_START(t_prolong_error);

    const int dst_idx = dst.getComponentDescriptorIndex(0);
    const int src_idx = src.getComponentDescriptorIndex(0);

    // Refine the correction from the coarse level src data into the fine
    // level error.
    xeqScheduleProlongation(dst_idx, src_idx, dst_ln);

    IBTK_TIMER_STOP(t_prolong_error);
    return;
}// prolongError

void
SCPoissonFACOperator::prolongErrorAndCorrect(
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
        HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops_coarse(d_hierarchy, dst_ln-1, dst_ln-1);
        hier_sc_data_ops_coarse.add(dst_idx, dst_idx, src_idx, interior_only);
    }
    xeqScheduleProlongation(d_side_scratch_idx, src_idx, dst_ln);
    HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops_fine(d_hierarchy, dst_ln, dst_ln);
    hier_sc_data_ops_fine.add(dst_idx, dst_idx, d_side_scratch_idx, interior_only);

    IBTK_TIMER_STOP(t_prolong_error_and_correct);
    return;
}// prolongErrorAndCorrect

void
SCPoissonFACOperator::smoothError(
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
    const int scratch_idx = d_side_scratch_idx;

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM,double> >   error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM,double> > scratch_data = patch->getPatchData(scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
            TBOX_ASSERT(  error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                scratch_data->getArrayData(axis).copy(
                    error_data->getArrayData(axis),
                    d_patch_bc_box_overlap[level_num][patch_counter][axis],
                    IntVector<NDIM>(0));
            }
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
                    Pointer<SideData<NDIM,double> >   error_data = error.getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM,double> > scratch_data = patch->getPatchData(scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                    const Box<NDIM>& ghost_box = error_data->getGhostBox();
                    TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
                    TBOX_ASSERT(  error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        error_data->getArrayData(axis).copy(
                            scratch_data->getArrayData(axis),
                            d_patch_bc_box_overlap[level_num][patch_counter][axis],
                            IntVector<NDIM>(0));
                    }
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
            Pointer<SideData<NDIM,double> >    error_data = error   .getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM,double> > residual_data = residual.getComponentPatchData(0, *patch);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == residual_data->getGhostBox());
            TBOX_ASSERT(   error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(residual_data->getGhostCellWidth() == d_gcw);
#endif
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                // Copy updated values from other local patches.
                if (d_smoother_choice == "multiplicative")
                {
                    const std::map<int,Box<NDIM> > smoother_bc_boxes = d_patch_smoother_bc_boxes[level_num][patch_counter][axis];
                    for (std::map<int,Box<NDIM> >::const_iterator cit = smoother_bc_boxes.begin();
                         cit != smoother_bc_boxes.end(); ++cit)
                    {
                        const int src_patch_num = cit->first;
                        const Box<NDIM>& overlap = cit->second;
                        Pointer<Patch<NDIM> > src_patch = level->getPatch(src_patch_num);
                        Pointer<SideData<NDIM,double> > src_error_data = error.getComponentPatchData(0, *src_patch);
                        error_data->getArrayData(axis).copy(src_error_data->getArrayData(axis), overlap, IntVector<NDIM>(0));
                    }
                }

                const Box<NDIM> side_patch_box = SideGeometry<NDIM>::toSideBox(patch_box,axis);

                // Smooth the error for each data depth.
                //
                // NOTE: Since the boundary conditions are handled "implicitly"
                // by setting ghost cell values, we can re-use the same patch
                // operators for each data depth even if different boundary
                // conditions are imposed on different components of the
                // vector-valued solution data.
                if (d_using_petsc_smoothers)
                {
                    // Reset ghost cell values in the copy of the residual data
                    // so that patch boundary conditions are properly handled.
                    residual_data->getArrayData(axis).copy(
                        error_data->getArrayData(axis),
                        d_patch_bc_box_overlap[level_num][patch_counter][axis],
                        IntVector<NDIM>(0));

                    for (int depth = 0; depth < d_depth; ++depth)
                    {
                        // Setup the PETSc Vec wrappers for the given patch
                        // data, axis, and data depth.
                        int ierr;

                        Vec& e = d_patch_vec_e[level_num][patch_counter][axis];
                        Vec& f = d_patch_vec_f[level_num][patch_counter][axis];

                        ierr = VecPlaceArray(e,    error_data->getPointer(axis,depth));  IBTK_CHKERRQ(ierr);
                        ierr = VecPlaceArray(f, residual_data->getPointer(axis,depth));  IBTK_CHKERRQ(ierr);

                        // Smooth the error on the patch using PETSc.  Here, we
                        // are approximately solving
                        //
                        //     Ae = f
                        //
                        // using an iteration of the form
                        //
                        //     e <- e + PC(f - Ae) = e + PC(r) = e + x.
                        //
                        // Presently, we simply employ symmetric Gauss-Seidel as
                        // the patch smoother.
                        static const double omega = 1.0;
                        static const double shift = 0.0;
                        static const int its = 1;
                        Mat& A = d_patch_mat[level_num][patch_counter][axis];
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
                    // Smooth the error using red-black Gauss-Seidel.
                    const double& alpha = d_poisson_spec.getDConstant();
                    const double& beta = d_poisson_spec.cIsZero() ? 0.0 : d_poisson_spec.getCConstant();
                    for (int depth = 0; depth < d_depth; ++depth)
                    {
                        double* const U = error_data->getPointer(axis,depth);
                        const int U_ghosts = (error_data->getGhostCellWidth()).max();
                        const double* const F = residual_data->getPointer(axis,depth);
                        const int F_ghosts = (residual_data->getGhostCellWidth()).max();
                        static const int its = 1;
                        RB_GS_SMOOTH_FC(
                            U, U_ghosts,
                            alpha, beta,
                            F, F_ghosts,
                            side_patch_box.lower(0), side_patch_box.upper(0),
                            side_patch_box.lower(1), side_patch_box.upper(1),
#if (NDIM == 3)
                            side_patch_box.lower(2), side_patch_box.upper(2),
#endif
                            dx, its);
                    }
                }
            }
        }
    }

    // Synchronize data along patch boundaries.
    xeqScheduleSideDataSynch(error_idx, level_num);

    IBTK_TIMER_STOP(t_smooth_error);
    return;
}// smoothError

bool
SCPoissonFACOperator::solveCoarsestLevel(
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
            if (d_depth > 1)
            {
                TBOX_ERROR("SCPoissonFACOperator::solveCoarsestLevel()\n"
                           << "  hypre level solver does not support non-scalar-valued data" << std::endl);
            }
            d_hypre_solver->setInitialGuessNonzero(true);
            d_hypre_solver->setMaxIterations(d_coarse_solver_max_its);
            d_hypre_solver->setRelativeTolerance(d_coarse_solver_tol);
            d_hypre_solver->solveSystem(error_level,residual_level);
        }
        else if (d_using_petsc)
        {
            d_petsc_solver->setInitialGuessNonzero(true);
            d_petsc_solver->setMaxIterations(d_coarse_solver_max_its);
            d_petsc_solver->setRelativeTolerance(d_coarse_solver_tol);
            d_petsc_solver->solveSystem(error_level,residual_level);
        }

        // Synchronize data along patch boundaries.
        const int error_idx = error.getComponentDescriptorIndex(0);
        xeqScheduleSideDataSynch(error_idx, coarsest_ln);
    }

    IBTK_TIMER_STOP(t_solve_coarsest_level);
    return true;
}// solveCoarsestLevel

void
SCPoissonFACOperator::computeResidual(
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

    const Pointer<SideVariable<NDIM,double> > res_var = residual.getComponentVariable(0);
    const Pointer<SideVariable<NDIM,double> > sol_var = solution.getComponentVariable(0);
    const Pointer<SideVariable<NDIM,double> > rhs_var = rhs.getComponentVariable(0);

    // Fill ghost-cell values.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    Pointer<SideNoCornersFillPattern> fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);
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
    HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    hier_sc_data_ops.axpy(res_idx, -1.0, res_idx, rhs_idx, false);

    IBTK_TIMER_STOP(t_compute_residual);
    return;
}// computeResidual

void
SCPoissonFACOperator::initializeOperatorState(
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

    Pointer<SideVariable<NDIM,double> > solution_var = solution.getComponentVariable(0);
    Pointer<SideVariable<NDIM,double> >      rhs_var =      rhs.getComponentVariable(0);

    Pointer<SideDataFactory<NDIM,double> > solution_pdat_fac = solution_var->getPatchDataFactory();
    Pointer<SideDataFactory<NDIM,double> >      rhs_pdat_fac =      rhs_var->getPatchDataFactory();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!solution_var.isNull());
    TBOX_ASSERT(!rhs_var.isNull());
    TBOX_ASSERT(!solution_pdat_fac.isNull());
    TBOX_ASSERT(!rhs_pdat_fac.isNull());
#endif

    if (solution_pdat_fac->getDefaultDepth() != rhs_pdat_fac->getDefaultDepth())
    {
        TBOX_ERROR("SCPoissonFACOperator::initializeOperatorState()\n"
                   << "  solution and rhs vectors must have the same data depths\n"
                   << "  solution data depth = " << solution_pdat_fac->getDefaultDepth() << "\n"
                   << "  rhs      data depth = " << rhs_pdat_fac     ->getDefaultDepth() << std::endl);
    }

    const int old_depth = d_depth;
    d_depth = solution_pdat_fac->getDefaultDepth();

    if (d_depth != old_depth)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<SideDataFactory<NDIM,double> > side_scratch_pdat_fac =
            var_db->getPatchDescriptor()->getPatchDataFactory(d_side_scratch_idx);
        side_scratch_pdat_fac->setDefaultDepth(d_depth);
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
        if (!level->checkAllocated(d_side_scratch_idx)) level->allocatePatchData(d_side_scratch_idx);
    }

    // Initialize the hypre solver when needed.
    if (d_using_hypre && (coarsest_reset_ln == d_coarsest_ln))
    {
        initializeHypreLevelSolver();
    }

    // Initialize the petsc solver when needed.
    if (d_using_petsc && (coarsest_reset_ln == d_coarsest_ln))
    {
        initializePETScLevelSolver();
    }

    // Get the transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    IBTK_DO_ONCE(
        geometry->addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen());
                 );

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    d_prolongation_refine_operator = geometry->lookupRefineOperator(var, d_prolongation_method);
    d_cf_bdry_op = new CartSideDoubleQuadraticCFInterpolation();
    d_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_cf_bdry_op->setPatchDataIndex(d_side_scratch_idx);
    d_cf_bdry_op->setPatchHierarchy(d_hierarchy);

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    d_restriction_coarsen_operator = geometry->lookupCoarsenOperator(var, d_restriction_method);
    d_ghostfill_nocoarse_refine_operator = NULL;
    d_side_synch_refine_operator = NULL;

    // Make space for saving communication schedules.  There is no need to
    // delete the old schedules first because we have deallocated the solver
    // state above.
    d_bc_op = new CartSideRobinPhysBdryOp(d_side_scratch_idx, d_bc_coefs, false);

    if (d_poisson_spec.dIsConstant())
    {
        d_op_stencil_fill_pattern = new SideNoCornersFillPattern(SIDEG, true, false, false);
    }
    else
    {
        d_op_stencil_fill_pattern = NULL;
    }
    d_side_synch_fill_pattern = new SideSynchCopyFillPattern();

    std::vector<RefinePatchStrategy<NDIM>*> prolongation_refine_patch_strategies;
    prolongation_refine_patch_strategies.push_back(d_cf_bdry_op);
    prolongation_refine_patch_strategies.push_back(d_bc_op);
    d_prolongation_refine_patch_strategy = new RefinePatchStrategySet(
        prolongation_refine_patch_strategies.begin(), prolongation_refine_patch_strategies.end(), false);

    d_prolongation_refine_schedules.resize(d_finest_ln+1);
    d_restriction_coarsen_schedules.resize(d_finest_ln+1);
    d_ghostfill_nocoarse_refine_schedules.resize(d_finest_ln+1);
    d_side_synch_refine_schedules.resize(d_finest_ln+1);

    d_prolongation_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_restriction_coarsen_algorithm = new CoarsenAlgorithm<NDIM>();
    d_ghostfill_nocoarse_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_side_synch_refine_algorithm = new RefineAlgorithm<NDIM>();

    d_prolongation_refine_algorithm->registerRefine(
        d_side_scratch_idx,
        solution.getComponentDescriptorIndex(0),
        d_side_scratch_idx,
        d_prolongation_refine_operator,
        d_op_stencil_fill_pattern);
    d_restriction_coarsen_algorithm->registerCoarsen(
        d_side_scratch_idx,
        rhs.getComponentDescriptorIndex(0),
        d_restriction_coarsen_operator);
    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        d_ghostfill_nocoarse_refine_operator,
        d_op_stencil_fill_pattern);
    d_side_synch_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        d_side_synch_refine_operator,
        d_side_synch_fill_pattern);

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

        d_side_synch_refine_schedules[dst_ln] =
            d_side_synch_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln));
    }

    d_ghostfill_nocoarse_refine_schedules[d_coarsest_ln] =
        d_ghostfill_nocoarse_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_ln), d_bc_op.getPointer());

    d_side_synch_refine_schedules[d_coarsest_ln] =
        d_side_synch_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_ln));

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
                const Box<NDIM>  ghost_box = Box<NDIM>::grow(patch_box, d_gcw);

                blitz::TinyVector<Vec,NDIM>& e = d_patch_vec_e[ln][patch_counter];
                blitz::TinyVector<Vec,NDIM>& f = d_patch_vec_f[ln][patch_counter];
                blitz::TinyVector<Mat,NDIM>& A = d_patch_mat  [ln][patch_counter];
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    int ierr;

                    const Box<NDIM> axis_ghost_box = SideGeometry<NDIM>::toSideBox(ghost_box,axis);
                    const int size = axis_ghost_box.size();

                    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, size, PETSC_NULL, &e[axis]);  IBTK_CHKERRQ(ierr);
                    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, size, PETSC_NULL, &f[axis]);  IBTK_CHKERRQ(ierr);

                    buildPatchLaplaceOperator(A[axis], d_poisson_spec, patch, axis, d_gcw);
                }
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
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box,axis);
                const Box<NDIM> side_ghost_box = Box<NDIM>::grow(side_box, 1);
                d_patch_bc_box_overlap[ln][patch_counter][axis] = BoxList<NDIM>(side_ghost_box);
                d_patch_bc_box_overlap[ln][patch_counter][axis].removeIntersections(side_box);
            }
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
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    d_patch_smoother_bc_boxes[ln][patch_counter1][axis].clear();
                }

                Pointer<Patch<NDIM> > dst_patch = level->getPatch(p1());
                const Box<NDIM>& dst_patch_box = dst_patch->getBox();
                const Box<NDIM>& dst_ghost_box = Box<NDIM>::grow(dst_patch_box, 1);

                int patch_counter2 = 0;
                for (PatchLevel<NDIM>::Iterator p2(level); patch_counter2 < patch_counter1; p2++, ++patch_counter2)
                {
                    Pointer<Patch<NDIM> > src_patch = level->getPatch(p2());
                    const Box<NDIM>& src_patch_box = src_patch->getBox();

                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        const Box<NDIM> overlap =
                            SideGeometry<NDIM>::toSideBox(dst_ghost_box,axis) *
                            SideGeometry<NDIM>::toSideBox(src_patch_box,axis);
                        if (!overlap.empty())
                        {
                            d_patch_smoother_bc_boxes[ln][patch_counter1][axis][p2()] = overlap;
                        }
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
SCPoissonFACOperator::deallocateOperatorState()
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
                if (level->checkAllocated(d_side_scratch_idx)) level->deallocatePatchData(d_side_scratch_idx);
            }

            if (d_using_petsc_smoothers)
            {
                int ierr;

                for (std::vector<blitz::TinyVector<Vec,NDIM> >::iterator it = d_patch_vec_e[ln].begin();
                     it != d_patch_vec_e[ln].end(); ++it)
                {
                    blitz::TinyVector<Vec,NDIM>& e = *it;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        ierr = VecDestroy(e[d]);  IBTK_CHKERRQ(ierr);
                    }
                }
                d_patch_vec_e[ln].clear();
                for (std::vector<blitz::TinyVector<Vec,NDIM> >::iterator it = d_patch_vec_f[ln].begin();
                     it != d_patch_vec_f[ln].end(); ++it)
                {
                    blitz::TinyVector<Vec,NDIM>& f = *it;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        ierr = VecDestroy(f[d]);  IBTK_CHKERRQ(ierr);
                    }
                }
                d_patch_vec_f[ln].clear();
                for (std::vector<blitz::TinyVector<Mat,NDIM> >::iterator it = d_patch_mat[ln].begin();
                     it != d_patch_mat[ln].end(); ++it)
                {
                    blitz::TinyVector<Mat,NDIM>& A = *it;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        ierr = MatDestroy(A[d]);  IBTK_CHKERRQ(ierr);
                    }
                }
                d_patch_mat[ln].clear();
            }
        }

        if (d_using_hypre && (coarsest_reset_ln == d_coarsest_ln))
        {
            d_hypre_solver->deallocateSolverState();
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

            d_side_synch_refine_operator .setNull();
            d_side_synch_refine_algorithm.setNull();
            d_side_synch_refine_schedules.resize(0);
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
SCPoissonFACOperator::xeqScheduleProlongation(
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
SCPoissonFACOperator::xeqScheduleRestriction(
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
SCPoissonFACOperator::xeqScheduleGhostFillNoCoarse(
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
SCPoissonFACOperator::xeqScheduleSideDataSynch(
    const int dst_idx,
    const int dst_ln)
{
    RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(dst_idx, dst_idx, dst_idx, d_side_synch_refine_operator, d_side_synch_fill_pattern);
    refiner.resetSchedule(d_side_synch_refine_schedules[dst_ln]);
    d_side_synch_refine_schedules[dst_ln]->fillData(d_apply_time);
    d_side_synch_refine_algorithm->resetSchedule(d_side_synch_refine_schedules[dst_ln]);
    return;
}// xeqScheduleSideDataSynch

void
SCPoissonFACOperator::initializeHypreLevelSolver()
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
    d_hypre_solver->setPoissonSpecifications(d_poisson_spec);
    d_hypre_solver->setPhysicalBcCoefs(d_bc_coefs);
    d_hypre_solver->setTime(d_apply_time);
    d_hypre_solver->setHomogeneousBc(true);
    d_hypre_solver->initializeSolverState(solution_level, rhs_level);
    return;
}// initializeHypreLevelSolver

void
SCPoissonFACOperator::initializePETScLevelSolver()
{
    d_petsc_solver = new SCPoissonPETScLevelSolver(d_object_name+"::petsc_solver", d_petsc_db);
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
SCPoissonFACOperator::buildPatchLaplaceOperator(
    Mat& A,
    const PoissonSpecifications& poisson_spec,
    const Pointer<Patch<NDIM> > patch,
    const int component_axis,
    const IntVector<NDIM>& ghost_cell_width)
{
    const double C = poisson_spec.cIsZero() ? 0.0 : poisson_spec.getCConstant();
    const double D = poisson_spec.getDConstant();

    int ierr;

    // Allocate a PETSc matrix for the patch operator.
    const Box<NDIM>& patch_box = patch->getBox();
    const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, component_axis);
    const Box<NDIM> ghost_box = Box<NDIM>::grow(side_box, ghost_cell_width);
    const int size = ghost_box.size();
    IntVector<NDIM> component_axis_shift = 0;
    component_axis_shift(component_axis) = -1;

    static const int stencil_sz = 2*NDIM+1;

    BoxList<NDIM> ghost_boxes(ghost_box);
    ghost_boxes.removeIntersections(side_box);
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

    for (Box<NDIM>::Iterator b(side_box); b; b++)
    {
        const Index<NDIM>& i = b();

        std::vector<double> mat_vals(stencil_sz,0.0);
        mat_vals[NDIM] = C;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const double& h = dx[axis];
            {
                mat_vals[NDIM-axis-1] += D/(h*h);
                mat_vals[NDIM       ] -= D/(h*h);
            }
            {
                mat_vals[NDIM+axis+1] += D/(h*h);
                mat_vals[NDIM       ] -= D/(h*h);
            }
        }

        static const int m = 1;
        static const int n = stencil_sz;
        std::vector<int> idxn(stencil_sz);
        const int idxm = ghost_box.offset(i);

        std::transform(mat_stencil.begin(), mat_stencil.end(), idxn.begin(), std::bind2nd(std::plus<int>(), idxm));
        ierr = MatSetValues(A, m, &idxm, n, &idxn[0], &mat_vals[0], INSERT_VALUES);  IBTK_CHKERRQ(ierr);
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

    // Assemble the matrices.
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    return;
}// buildPatchLaplaceOperator

void
SCPoissonFACOperator::sanityCheck()
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

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d_bc_coefs[d] == NULL)
        {
            TBOX_ERROR(d_object_name << ":\n"
                       << "  invalid physical bc object at depth = " << d << std::endl);
        }
    }
    return;
}// sanityCheck

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
