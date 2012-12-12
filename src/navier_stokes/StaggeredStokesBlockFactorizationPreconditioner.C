// Filename: StaggeredStokesBlockFactorizationPreconditioner.C
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

#include "StaggeredStokesBlockFactorizationPreconditioner.h"

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

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE     = "NONE";
static const bool        USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE    = "CUBIC_COARSEN";

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

StaggeredStokesBlockFactorizationPreconditioner::StaggeredStokesBlockFactorizationPreconditioner(
    const std::string& object_name,
    Pointer<Database> /*input_db*/,
    const std::string& /*default_options_prefix*/)
    : StaggeredStokesBlockPreconditioner(/*needs_velocity_solver*/ true, /*needs_pressure_solver*/ true),
      d_hierarchy(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_velocity_data_ops(NULL),
      d_pressure_data_ops(NULL),
      d_velocity_wgt_idx(-1),
      d_pressure_wgt_idx(-1),
      d_hier_math_ops(NULL),
      d_P_bdry_fill_op(NULL),
      d_no_fill_op(NULL),
      d_U_var(NULL),
      d_F_U_mod_idx(-1),
      d_P_var(NULL),
      d_P_scratch_idx(-1)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ true);

    // Present implementation requires zero initial guess and can perform only
    // one iteration.
    d_initial_guess_nonzero = false;
    d_max_iterations = 1;

    // Setup variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name+"::CONTEXT");

    const std::string U_var_name = d_object_name+"::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var)
    {
        d_F_U_mod_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
    else
    {
        d_U_var = new SideVariable<NDIM,double>(U_var_name);
        d_F_U_mod_idx = var_db->registerVariableAndContext(d_U_var, context, IntVector<NDIM>(SIDEG));
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_F_U_mod_idx >= 0);
#endif
    const std::string P_var_name = d_object_name+"::P";
    d_P_var = var_db->getVariable(P_var_name);
    if (d_P_var)
    {
        d_P_scratch_idx = var_db->mapVariableAndContextToIndex(d_P_var, context);
    }
    else
    {
        d_P_var = new CellVariable<NDIM,double>(P_var_name);
        d_P_scratch_idx = var_db->registerVariableAndContext(d_P_var, context, IntVector<NDIM>(CELLG));
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_P_scratch_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_solve_system            = TimerManager::getManager()->getTimer("IBAMR::StaggeredStokesBlockFactorizationPreconditioner::solveSystem()");
        t_initialize_solver_state = TimerManager::getManager()->getTimer("IBAMR::StaggeredStokesBlockFactorizationPreconditioner::initializeSolverState()");
        t_deallocate_solver_state = TimerManager::getManager()->getTimer("IBAMR::StaggeredStokesBlockFactorizationPreconditioner::deallocateSolverState()");
                  );
    return;
}// StaggeredStokesBlockFactorizationPreconditioner

StaggeredStokesBlockFactorizationPreconditioner::~StaggeredStokesBlockFactorizationPreconditioner()
{
    deallocateSolverState();
    return;
}// ~StaggeredStokesBlockFactorizationPreconditioner

bool
StaggeredStokesBlockFactorizationPreconditioner::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    IBAMR_TIMER_START(t_solve_system);

    // Initialize the solver (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeSolverState(x,b);

    // Set the initial guess to equal zero.
    x.setToScalar(0.0,false);

    // Get the vector components.
    const int F_U_idx = b.getComponentDescriptorIndex(0);
    const int F_P_idx = b.getComponentDescriptorIndex(1);

    const Pointer<Variable<NDIM> >& F_U_var = b.getComponentVariable(0);
    const Pointer<Variable<NDIM> >& F_P_var = b.getComponentVariable(1);

    Pointer<SideVariable<NDIM,double> > F_U_sc_var = F_U_var;
    Pointer<CellVariable<NDIM,double> > F_P_cc_var = F_P_var;

    const int U_idx = x.getComponentDescriptorIndex(0);
    const int P_idx = x.getComponentDescriptorIndex(1);

    const Pointer<Variable<NDIM> >& U_var = x.getComponentVariable(0);
    const Pointer<Variable<NDIM> >& P_var = x.getComponentVariable(1);

    Pointer<SideVariable<NDIM,double> > U_sc_var = U_var;
    Pointer<CellVariable<NDIM,double> > P_cc_var = P_var;

    // Setup the component solver vectors.
    Pointer<SAMRAIVectorReal<NDIM,double> > F_U_mod_vec;
    F_U_mod_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::F_U_mod", d_hierarchy, d_coarsest_ln, d_finest_ln);
    F_U_mod_vec->addComponent(d_U_var, d_F_U_mod_idx, d_velocity_wgt_idx, d_velocity_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > U_vec;
    U_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::U", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_vec->addComponent(U_sc_var, U_idx, d_velocity_wgt_idx, d_velocity_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > P_scratch_vec;
    P_scratch_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::P_scratch", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_scratch_vec->addComponent(d_P_var, d_P_scratch_idx, d_pressure_wgt_idx, d_pressure_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > F_P_vec;
    F_P_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::F_P", d_hierarchy, d_coarsest_ln, d_finest_ln);
    F_P_vec->addComponent(F_P_cc_var, F_P_idx, d_pressure_wgt_idx, d_pressure_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > P_vec;
    P_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::P", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_vec->addComponent(P_cc_var, P_idx, d_pressure_wgt_idx, d_pressure_data_ops);

    // Setup the interpolation transaction information.
    Pointer<VariableFillPattern<NDIM> > fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent         P_transaction_comp(          P_idx, DATA_REFINE_TYPE, USE_CF_INTERPOLATION, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef, fill_pattern);
    InterpolationTransactionComponent P_scratch_transaction_comp(d_P_scratch_idx, DATA_REFINE_TYPE, USE_CF_INTERPOLATION, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef, fill_pattern);

    // (1) Solve the pressure sub-problem.
    //
    // P := -(C*I+D*L) * (-L)^{-1} * F_P
    //
    // NOTE: When C=0, then
    //
    // P := -(D*L) * (-L)^{-1} * F_P
    //    = D*F_P
    if (d_U_problem_coefs.cIsZero() || MathUtilities<double>::equalEps(d_U_problem_coefs.getCConstant(),0.0))
    {
        d_pressure_data_ops->scale(P_idx, d_U_problem_coefs.getDConstant(), F_P_idx);
    }
    else
    {
        d_pressure_solver->setHomogeneousBc(true);
        LinearSolver* p_pressure_solver = dynamic_cast<LinearSolver*>(d_pressure_solver.getPointer());
        if (p_pressure_solver) p_pressure_solver->setInitialGuessNonzero(false);
        d_pressure_solver->solveSystem(*P_scratch_vec,*F_P_vec);
        d_hier_math_ops->laplace(P_idx, P_cc_var, d_U_problem_coefs, d_P_scratch_idx, d_P_var, d_P_bdry_fill_op, d_pressure_solver->getSolutionTime());
        d_pressure_data_ops->scale(P_idx, -1.0, P_idx);
    }
    d_P_bdry_fill_op->resetTransactionComponent(P_transaction_comp);
    d_P_bdry_fill_op->fillData(d_pressure_solver->getSolutionTime());
    d_P_bdry_fill_op->resetTransactionComponent(P_scratch_transaction_comp);

    // (2) Solve the velocity sub-problem.
    //
    // U := (C*I+D*L)^{-1} * [F_U + Grad * (C*I+D*L) * (-L)^{-1} * F_P]
    //    = (C*I+D*L)^{-1} * [F_U - Grad * P]
    static const bool cf_bdry_synch = true;
    d_hier_math_ops->grad(d_F_U_mod_idx, d_U_var, cf_bdry_synch, -1.0, P_idx, P_cc_var, d_no_fill_op, d_pressure_solver->getSolutionTime(), 1.0, F_U_idx, F_U_sc_var);
    d_velocity_solver->setHomogeneousBc(true);
    LinearSolver* p_velocity_solver = dynamic_cast<LinearSolver*>(d_velocity_solver.getPointer());
    if (p_velocity_solver) p_velocity_solver->setInitialGuessNonzero(false);
    d_velocity_solver->solveSystem(*U_vec,*F_U_mod_vec);

    // (3) Account for any nullspace vectors.
    if (p_velocity_solver)
    {
        const std::vector<Pointer<SAMRAIVectorReal<NDIM,double> > >& U_nul_vecs = p_velocity_solver->getNullspaceBasisVectors();
        if (!U_nul_vecs.empty())
        {
            for (unsigned int k = 0; k < U_nul_vecs.size(); ++k)
            {
                const double alpha = U_vec->dot(U_nul_vecs[k])/U_nul_vecs[k]->dot(U_nul_vecs[k]);
                U_vec->axpy(-alpha, U_nul_vecs[k], U_vec);
            }
        }
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!p_velocity_solver->getNullspaceContainsConstantVector());
#endif
    }
    LinearSolver* p_pressure_solver = dynamic_cast<LinearSolver*>(d_pressure_solver.getPointer());
    if (p_pressure_solver)
    {
        if (p_pressure_solver->getNullspaceContainsConstantVector())
        {
            const double volume = d_hier_math_ops->getVolumeOfPhysicalDomain();
            const double P_mean = (1.0/volume)*d_pressure_data_ops->integral(P_idx, d_pressure_wgt_idx);
            d_pressure_data_ops->addScalar(P_idx, P_idx, -P_mean);
        }
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(p_pressure_solver->getNullspaceBasisVectors().empty());
#endif
    }

    // Deallocate the solver (if necessary).
    if (deallocate_at_completion) deallocateSolverState();

    IBAMR_TIMER_STOP(t_solve_system);
    return true;
}// solveSystem

void
StaggeredStokesBlockFactorizationPreconditioner::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& b)
{
    IBAMR_TIMER_START(t_initialize_solver_state);

    if (d_is_initialized) deallocateSolverState();

    // Get the hierarchy configuration.
    d_hierarchy = x.getPatchHierarchy();
    d_coarsest_ln = x.getCoarsestLevelNumber();
    d_finest_ln = x.getFinestLevelNumber();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy == b.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == b.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == b.getFinestLevelNumber());
#else
    NULL_USE(b);
#endif

    // Setup hierarchy operators.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();

    d_velocity_data_ops = hier_ops_manager->getOperationsDouble(x.getComponentVariable(0), d_hierarchy, true);
    d_velocity_data_ops->setPatchHierarchy(d_hierarchy);
    d_velocity_data_ops->resetLevels(d_coarsest_ln, d_finest_ln);
    d_velocity_wgt_idx = x.getControlVolumeIndex(0);

    d_pressure_data_ops = hier_ops_manager->getOperationsDouble(x.getComponentVariable(1), d_hierarchy, true);
    d_pressure_data_ops->setPatchHierarchy(d_hierarchy);
    d_pressure_data_ops->resetLevels(d_coarsest_ln, d_finest_ln);
    d_pressure_wgt_idx = x.getControlVolumeIndex(1);

    d_hier_math_ops = new HierarchyMathOps(d_object_name+"::HierarchyMathOps", d_hierarchy, d_coarsest_ln, d_finest_ln);

    Pointer<VariableFillPattern<NDIM> > fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent P_scratch_component(d_P_scratch_idx, DATA_REFINE_TYPE, USE_CF_INTERPOLATION, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef, fill_pattern);
    d_P_bdry_fill_op = new HierarchyGhostCellInterpolation();
    d_P_bdry_fill_op->setHomogeneousBc(true);
    d_P_bdry_fill_op->initializeOperatorState(P_scratch_component, d_hierarchy);

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_F_U_mod_idx)) level->allocatePatchData(d_F_U_mod_idx);
        if (!level->checkAllocated(d_P_scratch_idx)) level->allocatePatchData(d_P_scratch_idx);
    }
    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_solver_state);
    return;
}// initializeSolverState

void
StaggeredStokesBlockFactorizationPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_solver_state);

    d_velocity_data_ops.setNull();
    d_pressure_data_ops.setNull();
    d_velocity_wgt_idx = -1;
    d_pressure_wgt_idx = -1;
    d_hier_math_ops.setNull();
    d_P_bdry_fill_op.setNull();
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_F_U_mod_idx)) level->deallocatePatchData(d_F_U_mod_idx);
        if (level->checkAllocated(d_P_scratch_idx)) level->deallocatePatchData(d_P_scratch_idx);
    }
    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_solver_state);
    return;
}// deallocateSolverState

void
StaggeredStokesBlockFactorizationPreconditioner::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    if (initial_guess_nonzero)
    {
        TBOX_ERROR(d_object_name+"::setInitialGuessNonzero()\n"
                   << "  class IBAMR::StaggeredStokesBlockFactorizationPreconditioner requires a zero initial guess" << std::endl);
    }
    return;
}// setInitialGuessNonzero

void
StaggeredStokesBlockFactorizationPreconditioner::setMaxIterations(
    int max_iterations)
{
    if (max_iterations != 1)
    {
        TBOX_ERROR(d_object_name+"::setMaxIterations()\n"
                   << "  class IBAMR::StaggeredStokesBlockFactorizationPreconditioner only performs a single iteration" << std::endl);
    }
    return;
}// setMaxIterations

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
