// Filename: SCPoissonPointRelaxationFACOperator.C
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

#include "SCPoissonPointRelaxationFACOperator.h"

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
#include <ibtk/RefinePatchStrategySet.h>
#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/SideSynchCopyFillPattern.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define GS_SMOOTH_FC FC_FUNC(gssmooth2d,GSSMOOTH2D)
#endif
#if (NDIM == 3)
#define GS_SMOOTH_FC FC_FUNC(gssmooth3d,GSSMOOTH3D)
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

SCPoissonPointRelaxationFACOperator::SCPoissonPointRelaxationFACOperator(
    const std::string& object_name,
    const Pointer<Database> input_db)
    : PoissonFACPreconditionerStrategy(object_name, PoissonSpecifications(object_name+"::poisson_spec"), new LocationIndexRobinBcCoefs<NDIM>(object_name+"::default_bc_coef", Pointer<Database>(NULL)), std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM,NULL), new SideVariable<NDIM,double>(object_name+"::side_scratch", DEFAULT_DATA_DEPTH), SIDEG, input_db),
      d_depth(DEFAULT_DATA_DEPTH),
      d_using_hypre(false),
      d_hypre_solver(NULL),
      d_hypre_db(),
      d_using_petsc(false),
      d_petsc_solver(NULL),
      d_petsc_db(),
      d_patch_bc_box_overlap(),
      d_patch_smoother_bc_boxes()
{
    // Get values from the input database.
    if (!input_db.isNull())
    {
        if (input_db->isDatabase("hypre_solver"))
        {
            d_hypre_db = input_db->getDatabase("hypre_solver");
        }
        if (input_db->isDatabase("petsc_solver"))
        {
            d_petsc_db = input_db->getDatabase("petsc_solver");
        }
    }

    // Create the hypre solver, if needed.
    setCoarsestLevelSolverChoice(d_coarse_solver_choice);

    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        LocationIndexRobinBcCoefs<NDIM>* p_default_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_bc_coef);
        p_default_bc_coef->setBoundaryValue(2*d  ,0.0);
        p_default_bc_coef->setBoundaryValue(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setPhysicalBcCoefs(blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>(d_default_bc_coef));

    // Setup Timers.
    IBTK_DO_ONCE(
        t_smooth_error         = TimerManager::getManager()->getTimer("IBTK::SCPoissonPointRelaxationFACOperator::smoothError()");
        t_solve_coarsest_level = TimerManager::getManager()->getTimer("IBTK::SCPoissonPointRelaxationFACOperator::solveCoarsestLevel()");
        t_compute_residual     = TimerManager::getManager()->getTimer("IBTK::SCPoissonPointRelaxationFACOperator::computeResidual()");
                 );
    return;
}// SCPoissonPointRelaxationFACOperator

SCPoissonPointRelaxationFACOperator::~SCPoissonPointRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}// ~SCPoissonPointRelaxationFACOperator

void
SCPoissonPointRelaxationFACOperator::setSmootherChoice(
    const std::string& smoother_choice)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSmootherChoice():\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    if (smoother_choice != "additive" && smoother_choice != "multiplicative")
    {
        TBOX_ERROR(d_object_name << "::setSmootherChoice():\n"
                   << "  unknown smoother type: " << smoother_choice << "\n"
                   << "  valid choices are: additive, multiplicative" << std::endl);
    }
    d_smoother_choice = smoother_choice;
    return;
}// setSmootherChoice

void
SCPoissonPointRelaxationFACOperator::setCoarsestLevelSolverChoice(
    const std::string& coarse_solver_choice)
{
     if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setCoarsestLevelSolverChoice():\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    if (coarse_solver_choice != "block_jacobi" && coarse_solver_choice != "hypre" && coarse_solver_choice != "petsc")
    {
        TBOX_ERROR(d_object_name << "::setCoarsestLevelSolverChoice():\n"
                   << "  unknown coarse solver type: " << d_coarse_solver_choice << "\n"
                   << "  valid choices are: block_jacobi, hypre, petsc" << std::endl);
    }
   d_coarse_solver_choice = coarse_solver_choice;

    if (d_coarse_solver_choice == "hypre")
    {
        d_using_hypre = true;
    }
    else
    {
        d_using_hypre = false;
        d_hypre_solver.setNull();
    }

    if (d_coarse_solver_choice == "petsc")
    {
        d_using_petsc = true;
    }
    else
    {
        d_using_petsc = false;
        d_petsc_solver.setNull();
    }
    return;
}// setCoarsestLevelSolverChoice

void
SCPoissonPointRelaxationFACOperator::smoothError(
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
    const int scratch_idx = d_scratch_idx;

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
                scratch_data->getArrayData(axis).copy(error_data->getArrayData(axis), d_patch_bc_box_overlap[level_num][patch_counter][axis], IntVector<NDIM>(0));
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
                        error_data->getArrayData(axis).copy(scratch_data->getArrayData(axis), d_patch_bc_box_overlap[level_num][patch_counter][axis], IntVector<NDIM>(0));
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
                    for (std::map<int,Box<NDIM> >::const_iterator cit = smoother_bc_boxes.begin(); cit != smoother_bc_boxes.end(); ++cit)
                    {
                        const int src_patch_num = cit->first;
                        const Box<NDIM>& overlap = cit->second;
                        Pointer<Patch<NDIM> > src_patch = level->getPatch(src_patch_num);
                        Pointer<SideData<NDIM,double> > src_error_data = error.getComponentPatchData(0, *src_patch);
                        error_data->getArrayData(axis).copy(src_error_data->getArrayData(axis), overlap, IntVector<NDIM>(0));
                    }
                }

                // Smooth the error using Gauss-Seidel.
                const double& alpha = d_poisson_spec.getDConstant();
                const double& beta = d_poisson_spec.cIsZero() ? 0.0 : d_poisson_spec.getCConstant();
                const Box<NDIM> side_patch_box = SideGeometry<NDIM>::toSideBox(patch_box,axis);
                for (int depth = 0; depth < d_depth; ++depth)
                {
                    double* const U = error_data->getPointer(axis,depth);
                    const int U_ghosts = (error_data->getGhostCellWidth()).max();
                    const double* const F = residual_data->getPointer(axis,depth);
                    const int F_ghosts = (residual_data->getGhostCellWidth()).max();
                    static const int its = 1;
                    GS_SMOOTH_FC(
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

    // Synchronize data along patch boundaries.
    xeqScheduleDataSynch(error_idx, level_num);
    IBTK_TIMER_STOP(t_smooth_error);
    return;
}// smoothError

bool
SCPoissonPointRelaxationFACOperator::solveCoarsestLevel(
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
                TBOX_ERROR("SCPoissonPointRelaxationFACOperator::solveCoarsestLevel()\n"
                           << "  hypre level solver does not support non-scalar-valued data" << std::endl);
            }
            d_hypre_solver->setSolutionTime(d_solution_time);
            d_hypre_solver->setTimeInterval(d_current_time, d_new_time);
            d_hypre_solver->setInitialGuessNonzero(true);
            d_hypre_solver->setMaxIterations(d_coarse_solver_max_its);
            d_hypre_solver->setRelativeTolerance(d_coarse_solver_tol);
            d_hypre_solver->solveSystem(error_level,residual_level);
        }
        else if (d_using_petsc)
        {
            d_petsc_solver->setSolutionTime(d_solution_time);
            d_petsc_solver->setTimeInterval(d_current_time, d_new_time);
            d_petsc_solver->setInitialGuessNonzero(true);
            d_petsc_solver->setMaxIterations(d_coarse_solver_max_its);
            d_petsc_solver->setRelativeTolerance(d_coarse_solver_tol);
            d_petsc_solver->solveSystem(error_level,residual_level);
        }

        // Synchronize data along patch boundaries.
        const int error_idx = error.getComponentDescriptorIndex(0);
        xeqScheduleDataSynch(error_idx, coarsest_ln);
    }
    IBTK_TIMER_STOP(t_solve_coarsest_level);
    return true;
}// solveCoarsestLevel

void
SCPoissonPointRelaxationFACOperator::computeResidual(
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
    d_hier_bdry_fill_ops[finest_level_num]->fillData(d_solution_time);

    // Compute the residual, r = f - A*u.
    if (d_hier_math_ops[finest_level_num].isNull())
    {
        std::ostringstream stream;
        stream << d_object_name << "::hier_math_ops_" << finest_level_num;
        d_hier_math_ops[finest_level_num] = new HierarchyMathOps(stream.str(), d_hierarchy, coarsest_level_num, finest_level_num);
    }
    d_hier_math_ops[finest_level_num]->laplace(res_idx, res_var, d_poisson_spec, sol_idx, sol_var, NULL, d_solution_time);
    HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    hier_sc_data_ops.axpy(res_idx, -1.0, res_idx, rhs_idx, false);

    IBTK_TIMER_STOP(t_compute_residual);
    return;
}// computeResidual

/////////////////////////////// PROTECTED ////////////////////////////////////


void
SCPoissonPointRelaxationFACOperator::initializeOperatorStateSpecialized(
    const SAMRAIVectorReal<NDIM,double>& solution,
    const SAMRAIVectorReal<NDIM,double>& rhs,
    const int coarsest_reset_ln,
    const int finest_reset_ln)
{
    // Setup solution and rhs vectors.
    Pointer<SideVariable<NDIM,double> > solution_var = solution.getComponentVariable(0);
    Pointer<SideVariable<NDIM,double> >      rhs_var =      rhs.getComponentVariable(0);

    Pointer<SideDataFactory<NDIM,double> > solution_pdat_fac = solution_var->getPatchDataFactory();
    Pointer<SideDataFactory<NDIM,double> >      rhs_pdat_fac =      rhs_var->getPatchDataFactory();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!solution_var.isNull());
    TBOX_ASSERT(!     rhs_var.isNull());
    TBOX_ASSERT(!solution_pdat_fac.isNull());
    TBOX_ASSERT(!     rhs_pdat_fac.isNull());
#endif

    if (solution_pdat_fac->getDefaultDepth() != rhs_pdat_fac->getDefaultDepth())
    {
        TBOX_ERROR("SCPoissonPointRelaxationFACOperator::initializeOperatorState()\n"
                   << "  solution and rhs vectors must have the same data depths\n"
                   << "  solution data depth = " << solution_pdat_fac->getDefaultDepth() << "\n"
                   << "  rhs      data depth = " << rhs_pdat_fac     ->getDefaultDepth() << std::endl);
    }

    const int old_depth = d_depth;
    d_depth = solution_pdat_fac->getDefaultDepth();

    if (d_depth != old_depth)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<SideDataFactory<NDIM,double> > scratch_pdat_fac = var_db->getPatchDescriptor()->getPatchDataFactory(d_scratch_idx);
        scratch_pdat_fac->setDefaultDepth(d_depth);
    }

    // Initialize the bottom solvers when needed.
    if (coarsest_reset_ln == d_coarsest_ln)
    {
        if (d_using_hypre) initializeHypreLevelSolver();
        if (d_using_petsc) initializePETScLevelSolver();
    }

    // Setup specialized transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    IBTK_DO_ONCE(
        geometry->addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen());
                 );

    // Setup coarse-fine interface and physical boundary operators.
    d_cf_bdry_op = new CartSideDoubleQuadraticCFInterpolation();
    d_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_cf_bdry_op->setPatchDataIndex(d_scratch_idx);
    d_cf_bdry_op->setPatchHierarchy(d_hierarchy);
    d_bc_op = new CartSideRobinPhysBdryOp(d_scratch_idx, d_bc_coefs, false);

    // Setup fill pattern spec objects.
    if (d_poisson_spec.dIsConstant())
    {
        d_op_stencil_fill_pattern = new SideNoCornersFillPattern(SIDEG, true, false, false);
    }
    else
    {
        d_op_stencil_fill_pattern = NULL;
    }
    d_synch_fill_pattern = new SideSynchCopyFillPattern();

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
                        const Box<NDIM> overlap = SideGeometry<NDIM>::toSideBox(dst_ghost_box,axis) * SideGeometry<NDIM>::toSideBox(src_patch_box,axis);
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
    return;
}// initializeOperatorStateSpecialized

void
SCPoissonPointRelaxationFACOperator::deallocateOperatorStateSpecialized(
    const int /*coarsest_reset_ln*/,
    const int /*finest_reset_ln*/)
{
    if (!d_is_initialized) return;

    if (!d_in_initialize_operator_state)
    {
        d_patch_bc_box_overlap.clear();
        d_patch_smoother_bc_boxes.clear();
        if (d_using_hypre) d_hypre_solver->deallocateSolverState();
        if (d_using_petsc) d_petsc_solver->deallocateSolverState();
    }
    return;
}// deallocateOperatorStateSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SCPoissonPointRelaxationFACOperator::initializeHypreLevelSolver()
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
    d_hypre_solver = new SCPoissonHypreLevelSolver(d_object_name+"::hypre_solver", d_hypre_db);
    d_hypre_solver->setSolutionTime(d_solution_time);
    d_hypre_solver->setTimeInterval(d_current_time, d_new_time);
    d_hypre_solver->setPoissonSpecifications(d_poisson_spec);
    d_hypre_solver->setPhysicalBcCoefs(d_bc_coefs);
    d_hypre_solver->setHomogeneousBc(true);
    d_hypre_solver->initializeSolverState(solution_level, rhs_level);
    return;
}// initializeHypreLevelSolver

void
SCPoissonPointRelaxationFACOperator::initializePETScLevelSolver()
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
    d_petsc_solver = new SCPoissonPETScLevelSolver(d_object_name+"::petsc_solver", d_petsc_db);
    d_petsc_solver->setSolutionTime(d_solution_time);
    d_petsc_solver->setTimeInterval(d_current_time, d_new_time);
    d_petsc_solver->setPoissonSpecifications(d_poisson_spec);
    d_petsc_solver->setPhysicalBcCoefs(d_bc_coefs);
    d_petsc_solver->setHomogeneousBc(true);
    d_petsc_solver->initializeSolverState(solution_level, rhs_level);
    return;
}// initializePETScLevelSolver

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
