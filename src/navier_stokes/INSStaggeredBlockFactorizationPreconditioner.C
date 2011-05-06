// Filename: INSStaggeredBlockFactorizationPreconditioner.C
// Created on 22 Sep 2008 by Boyce Griffith
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

#include "INSStaggeredBlockFactorizationPreconditioner.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CellNoCornersFillPattern.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_solve_system;
static Timer* t_initialize_solver_state;
static Timer* t_deallocate_solver_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredBlockFactorizationPreconditioner::INSStaggeredBlockFactorizationPreconditioner(
    const INSCoefs& problem_coefs,
    RobinBcCoefStrategy<NDIM>* P_bc_coef,
    const bool normalize_pressure,
    Pointer<LinearSolver> velocity_helmholtz_solver,
    Pointer<LinearSolver> pressure_poisson_solver,
    Pointer<HierarchyCellDataOpsReal<NDIM,double> > hier_cc_data_ops,
    Pointer<HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops,
    Pointer<HierarchyMathOps> hier_math_ops)
    : d_do_log(false),
      d_is_initialized(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_problem_coefs(problem_coefs),
      d_pressure_helmholtz_spec("INSStaggeredBlockFactorizationPreconditioner::pressure_helmholtz_spec"),
      d_normalize_pressure(normalize_pressure),
      d_velocity_helmholtz_solver(velocity_helmholtz_solver),
      d_pressure_poisson_solver(pressure_poisson_solver),
      d_hier_cc_data_ops(hier_cc_data_ops),
      d_hier_sc_data_ops(hier_sc_data_ops),
      d_hier_math_ops(hier_math_ops),
      d_wgt_cc_var(NULL),
      d_wgt_sc_var(NULL),
      d_wgt_cc_idx(-1),
      d_wgt_sc_idx(-1),
      d_volume(std::numeric_limits<double>::quiet_NaN()),
      d_P_bc_coef(P_bc_coef),
      d_P_bdry_fill_op(Pointer<HierarchyGhostCellInterpolation>(NULL)),
      d_no_fill_op(NULL),
      d_hierarchy(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_U_var(NULL),
      d_U_scratch_idx(-1),
      d_P_var(NULL),
      d_P_scratch_idx(-1)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("INSStaggeredBlockFactorizationPreconditionerOperator::CONTEXT");

    const std::string U_var_name = "INSStaggeredBlockFactorizationPreconditioner::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var.isNull())
    {
        d_U_var = new SideVariable<NDIM,double>(U_var_name);
        d_U_scratch_idx = var_db->registerVariableAndContext(d_U_var, context, IntVector<NDIM>(SIDEG));
    }
    else
    {
        d_U_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_U_scratch_idx >= 0);
#endif
    const std::string P_var_name = "INSStaggeredBlockFactorizationPreconditioner::P";
    d_P_var = var_db->getVariable(P_var_name);
    if (d_P_var.isNull())
    {
        d_P_var = new CellVariable<NDIM,double>(P_var_name);
        d_P_scratch_idx = var_db->registerVariableAndContext(d_P_var, context, IntVector<NDIM>(CELLG));
    }
    else
    {
        d_P_scratch_idx = var_db->mapVariableAndContextToIndex(d_P_var, context);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_P_scratch_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_solve_system            = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredBlockFactorizationPreconditioner::solveSystem()");
        t_initialize_solver_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredBlockFactorizationPreconditioner::initializeSolverState()");
        t_deallocate_solver_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredBlockFactorizationPreconditioner::deallocateSolverState()");
                  );
    return;
}// INSStaggeredBlockFactorizationPreconditioner

INSStaggeredBlockFactorizationPreconditioner::~INSStaggeredBlockFactorizationPreconditioner()
{
    deallocateSolverState();
    return;
}// ~INSStaggeredBlockFactorizationPreconditioner

void
INSStaggeredBlockFactorizationPreconditioner::setTimeInterval(
    const double current_time,
    const double new_time)
{
    const double rho    = d_problem_coefs.getRho();
    const double mu     = d_problem_coefs.getMu();
    const double lambda = d_problem_coefs.getLambda();
    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = d_new_time-d_current_time;
    d_pressure_helmholtz_spec.setCConstant(-(rho/d_dt+0.5*lambda));
    d_pressure_helmholtz_spec.setDConstant(-(        -0.5*mu    ));
    return;
}// setTimeInterval

bool
INSStaggeredBlockFactorizationPreconditioner::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    t_solve_system->start();

    // Initialize the solver (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeSolverState(x,b);

    // Get the vector components.
    const int U_in_idx = b.getComponentDescriptorIndex(0);
    const int P_in_idx = b.getComponentDescriptorIndex(1);

    const Pointer<Variable<NDIM> >& U_in_var = b.getComponentVariable(0);
    const Pointer<Variable<NDIM> >& P_in_var = b.getComponentVariable(1);

    Pointer<SideVariable<NDIM,double> > U_in_sc_var = U_in_var;
    Pointer<CellVariable<NDIM,double> > P_in_cc_var = P_in_var;

    const int U_out_idx = x.getComponentDescriptorIndex(0);
    const int P_out_idx = x.getComponentDescriptorIndex(1);

    const Pointer<Variable<NDIM> >& U_out_var = x.getComponentVariable(0);
    const Pointer<Variable<NDIM> >& P_out_var = x.getComponentVariable(1);

    Pointer<SideVariable<NDIM,double> > U_out_sc_var = U_out_var;
    Pointer<CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

    // Setup the component solver vectors.
    Pointer<SAMRAIVectorReal<NDIM,double> > U_scratch_vec;
    U_scratch_vec = new SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredBlockFactorizationPreconditioner::U_scratch", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_scratch_vec->addComponent(d_U_var, d_U_scratch_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > P_scratch_vec;
    P_scratch_vec = new SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredBlockFactorizationPreconditioner::P_scratch", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_scratch_vec->addComponent(d_P_var, d_P_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > U_out_vec;
    U_out_vec = new SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredBlockFactorizationPreconditioner::U_out", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_out_vec->addComponent(U_out_sc_var, U_out_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > P_in_vec;
    P_in_vec = new SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredBlockFactorizationPreconditioner::P_in", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_in_vec->addComponent(P_in_cc_var, P_in_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > P_out_vec;
    P_out_vec = new SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredBlockFactorizationPreconditioner::P_out", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_out_vec->addComponent(P_out_cc_var, P_out_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    // Setup the interpolation transaction information.
    Pointer<VariableFillPattern<NDIM> > fill_pattern = new CellNoCornersFillPattern(CELLG, false, true);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent     P_out_transaction_comp(      P_out_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef, fill_pattern);
    InterpolationTransactionComponent P_scratch_transaction_comp(d_P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef, fill_pattern);

    // Solve the pressure sub-problem.
    //
    // P_out := -((rho/dt)*I-0.5*mu*L) * (-L)^{-1} * P_in
    d_pressure_poisson_solver->solveSystem(*P_scratch_vec,*P_in_vec);
    d_hier_math_ops->laplace(
        P_out_idx, P_out_cc_var,   // dst
        d_pressure_helmholtz_spec, // Poisson spec
        d_P_scratch_idx, d_P_var,  // src
        d_P_bdry_fill_op,          // src_bdry_fill
        d_current_time+0.5*d_dt);  // src_bdry_fill_time
    if (d_normalize_pressure)
    {
        const double P_mean = (1.0/d_volume)*d_hier_cc_data_ops->integral(P_out_idx, d_wgt_cc_idx);
        d_hier_cc_data_ops->addScalar(P_out_idx, P_out_idx, -P_mean);
    }

    // Compute the right-hand-side for the Helmholtz solve.
    d_P_bdry_fill_op->resetTransactionComponent(P_out_transaction_comp);
    static const bool cf_bdry_synch = true;
    d_hier_math_ops->grad(
        d_U_scratch_idx, d_U_var, // dst
        cf_bdry_synch,            // dst_cf_bdry_synch
        -1.0,                     // alpha
        P_out_idx, P_out_cc_var,  // src1
        d_P_bdry_fill_op,         // src1_bdry_fill
        d_current_time+0.5*d_dt,  // src1_bdry_fill_time
        1.0,                      // beta
        U_in_idx, U_in_sc_var);   // src2
    d_P_bdry_fill_op->resetTransactionComponent(P_scratch_transaction_comp);

    // Solve the velocity sub-problem.
    //
    // U_out := (rho/dt)*I-0.5*mu*L)^{-1} * [U_in + Grad * ((rho/dt)*I-0.5*mu*L) * (-L)^{-1} * P_in]
    //        = (rho/dt)*I-0.5*mu*L)^{-1} * [U_in - Grad * P_out]
    d_velocity_helmholtz_solver->solveSystem(*U_out_vec,*U_scratch_vec);

    // Deallocate the solver (if necessary).
    if (deallocate_at_completion) deallocateSolverState();

    t_solve_system->stop();
    return true;
}// solveSystem

void
INSStaggeredBlockFactorizationPreconditioner::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& b)
{
    t_initialize_solver_state->start();

    if (d_is_initialized) deallocateSolverState();

    // Get the hierarchy configuration.
    d_hierarchy = x.getPatchHierarchy();
    d_coarsest_ln = x.getCoarsestLevelNumber();
    d_finest_ln = x.getFinestLevelNumber();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy == b.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == b.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == b.getFinestLevelNumber());
#endif
    d_wgt_cc_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_sc_var = d_hier_math_ops->getSideWeightVariable();
    d_wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    d_wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    Pointer<VariableFillPattern<NDIM> > fill_pattern = new CellNoCornersFillPattern(CELLG, false, true);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent P_scratch_component(d_P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef, fill_pattern);
    d_P_bdry_fill_op = new HierarchyGhostCellInterpolation();
    d_P_bdry_fill_op->initializeOperatorState(P_scratch_component, d_hierarchy);

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_U_scratch_idx))
        {
            level->allocatePatchData(d_U_scratch_idx);
        }
        if (!level->checkAllocated(d_P_scratch_idx))
        {
            level->allocatePatchData(d_P_scratch_idx);
        }
    }
    d_is_initialized = true;

    t_initialize_solver_state->stop();
    return;
}// initializeSolverState

void
INSStaggeredBlockFactorizationPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    t_deallocate_solver_state->start();

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_U_scratch_idx))
        {
            level->deallocatePatchData(d_U_scratch_idx);
        }
        if (level->checkAllocated(d_P_scratch_idx))
        {
            level->deallocatePatchData(d_P_scratch_idx);
        }
    }
    d_is_initialized = false;

    t_deallocate_solver_state->stop();
    return;
}// deallocateSolverState

void
INSStaggeredBlockFactorizationPreconditioner::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    // intentionally blank
    return;
}// setInitialGuessNonzero

bool
INSStaggeredBlockFactorizationPreconditioner::getInitialGuessNonzero() const
{
    // intentionally blank
    return true;
}// getInitialGuessNonzero

void
INSStaggeredBlockFactorizationPreconditioner::setMaxIterations(
    int max_iterations)
{
    // intentionally blank
    return;
}// setMaxIterations

int
INSStaggeredBlockFactorizationPreconditioner::getMaxIterations() const
{
    // intentionally blank
    return 1;
}// getMaxIterations

void
INSStaggeredBlockFactorizationPreconditioner::setAbsoluteTolerance(
    double abs_residual_tol)
{
    // intentionally blank
    return;
}// setAbsoluteTolerance

double
INSStaggeredBlockFactorizationPreconditioner::getAbsoluteTolerance() const
{
    // intentionally blank
    return 0.0;
}// getAbsoluteTolerance

void
INSStaggeredBlockFactorizationPreconditioner::setRelativeTolerance(
    double rel_residual_tol)
{
    // intentionally blank
    return;
}// setRelativeTolerance

double
INSStaggeredBlockFactorizationPreconditioner::getRelativeTolerance() const
{
    // intentionally blank
    return 0.0;
}// getRelativeTolerance

int
INSStaggeredBlockFactorizationPreconditioner::getNumIterations() const
{
    // intentionally blank
    return 0;
}// getNumIterations

double
INSStaggeredBlockFactorizationPreconditioner::getResidualNorm() const
{
    return 0.0;
}// getResidualNorm

void
INSStaggeredBlockFactorizationPreconditioner::enableLogging(
    bool enabled)
{
    d_do_log = enabled;
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::INSStaggeredBlockFactorizationPreconditioner>;

//////////////////////////////////////////////////////////////////////////////
