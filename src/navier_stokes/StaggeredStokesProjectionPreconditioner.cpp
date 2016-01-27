// Filename: StaggeredStokesProjectionPreconditioner.cpp
// Created on 29 Apr 2008 by Boyce Griffith
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
#include <ostream>
#include <string>

#include "CellVariable.h"
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VariableFillPattern.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesProjectionPreconditioner.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PoissonSolver.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
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

StaggeredStokesProjectionPreconditioner::StaggeredStokesProjectionPreconditioner(
    const std::string& object_name,
    Pointer<Database> /*input_db*/,
    const std::string& /*default_options_prefix*/)
    : StaggeredStokesBlockPreconditioner(/*needs_velocity_solver*/ true,
                                         /*needs_pressure_solver*/ true),
      d_Phi_bdry_fill_op(NULL),
      d_no_fill_op(NULL),
      d_Phi_var(NULL),
      d_F_Phi_var(NULL),
      d_Phi_scratch_idx(-1),
      d_F_Phi_idx(-1)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ true);

    // Present implementation requires zero initial guess and can perform only
    // one iteration.
    d_initial_guess_nonzero = false;
    d_max_iterations = 1;

    // Setup variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONTEXT");

    const std::string Phi_var_name = d_object_name + "::Phi";
    d_Phi_var = var_db->getVariable(Phi_var_name);
    if (d_Phi_var)
    {
        d_Phi_scratch_idx = var_db->mapVariableAndContextToIndex(d_Phi_var, context);
    }
    else
    {
        d_Phi_var = new CellVariable<NDIM, double>(Phi_var_name);
        d_Phi_scratch_idx = var_db->registerVariableAndContext(d_Phi_var, context, IntVector<NDIM>(CELLG));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_Phi_scratch_idx >= 0);
#endif
    const std::string F_var_name = d_object_name + "::F";
    d_F_Phi_var = var_db->getVariable(F_var_name);
    if (d_F_Phi_var)
    {
        d_F_Phi_idx = var_db->mapVariableAndContextToIndex(d_F_Phi_var, context);
    }
    else
    {
        d_F_Phi_var = new CellVariable<NDIM, double>(F_var_name);
        d_F_Phi_idx = var_db->registerVariableAndContext(d_F_Phi_var, context, IntVector<NDIM>(CELLG));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_F_Phi_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(t_solve_system = TimerManager::getManager()->getTimer(
                      "IBAMR::StaggeredStokesProjectionPreconditioner::solveSystem()");
                  t_initialize_solver_state = TimerManager::getManager()->getTimer(
                      "IBAMR::StaggeredStokesProjectionPreconditioner::initializeSolverState()");
                  t_deallocate_solver_state = TimerManager::getManager()->getTimer(
                      "IBAMR::StaggeredStokesProjectionPreconditioner::deallocateSolverState()"););
    return;
} // StaggeredStokesProjectionPreconditioner

StaggeredStokesProjectionPreconditioner::~StaggeredStokesProjectionPreconditioner()
{
    deallocateSolverState();
    return;
} // ~StaggeredStokesProjectionPreconditioner

bool
StaggeredStokesProjectionPreconditioner::solveSystem(SAMRAIVectorReal<NDIM, double>& x,
                                                     SAMRAIVectorReal<NDIM, double>& b)
{
    IBAMR_TIMER_START(t_solve_system);

    // Initialize the solver (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeSolverState(x, b);

    // Determine whether we are solving a steady-state problem.
    const bool steady_state =
        d_U_problem_coefs.cIsZero() ||
        (d_U_problem_coefs.cIsConstant() && MathUtilities<double>::equalEps(d_U_problem_coefs.getCConstant(), 0.0));

    // Get the vector components.
    const int F_U_idx = b.getComponentDescriptorIndex(0);
    const int F_P_idx = b.getComponentDescriptorIndex(1);

    const Pointer<Variable<NDIM> >& F_U_var = b.getComponentVariable(0);
    const Pointer<Variable<NDIM> >& F_P_var = b.getComponentVariable(1);

    Pointer<SideVariable<NDIM, double> > F_U_sc_var = F_U_var;
    Pointer<CellVariable<NDIM, double> > F_P_cc_var = F_P_var;

    const int U_idx = x.getComponentDescriptorIndex(0);
    const int P_idx = x.getComponentDescriptorIndex(1);

    const Pointer<Variable<NDIM> >& U_var = x.getComponentVariable(0);
    const Pointer<Variable<NDIM> >& P_var = x.getComponentVariable(1);

    Pointer<SideVariable<NDIM, double> > U_sc_var = U_var;
    Pointer<CellVariable<NDIM, double> > P_cc_var = P_var;

    // Setup the component solver vectors.
    Pointer<SAMRAIVectorReal<NDIM, double> > F_U_vec;
    F_U_vec = new SAMRAIVectorReal<NDIM, double>(d_object_name + "::F_U", d_hierarchy, d_coarsest_ln, d_finest_ln);
    F_U_vec->addComponent(F_U_sc_var, F_U_idx, d_velocity_wgt_idx, d_velocity_data_ops);

    Pointer<SAMRAIVectorReal<NDIM, double> > U_vec;
    U_vec = new SAMRAIVectorReal<NDIM, double>(d_object_name + "::U", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_vec->addComponent(U_sc_var, U_idx, d_velocity_wgt_idx, d_velocity_data_ops);

    Pointer<SAMRAIVectorReal<NDIM, double> > Phi_scratch_vec;
    Phi_scratch_vec =
        new SAMRAIVectorReal<NDIM, double>(d_object_name + "::Phi_scratch", d_hierarchy, d_coarsest_ln, d_finest_ln);
    Phi_scratch_vec->addComponent(d_Phi_var, d_Phi_scratch_idx, d_pressure_wgt_idx, d_pressure_data_ops);

    Pointer<SAMRAIVectorReal<NDIM, double> > F_Phi_vec;
    F_Phi_vec = new SAMRAIVectorReal<NDIM, double>(d_object_name + "::F_Phi", d_hierarchy, d_coarsest_ln, d_finest_ln);
    F_Phi_vec->addComponent(d_F_Phi_var, d_F_Phi_idx, d_pressure_wgt_idx, d_pressure_data_ops);

    Pointer<SAMRAIVectorReal<NDIM, double> > P_vec;
    P_vec = new SAMRAIVectorReal<NDIM, double>(d_object_name + "::P", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_vec->addComponent(P_cc_var, P_idx, d_pressure_wgt_idx, d_pressure_data_ops);

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_Phi_scratch_idx);
        level->allocatePatchData(d_F_Phi_idx);
    }

    // (1) Solve the velocity sub-problem for an initial approximation to U.
    //
    // U^* := inv(rho/dt - K*mu*L) F_U
    //
    // An approximate Helmholtz solver is used.
    d_velocity_solver->setHomogeneousBc(true);
    LinearSolver* p_velocity_solver = dynamic_cast<LinearSolver*>(d_velocity_solver.getPointer());
    if (p_velocity_solver) p_velocity_solver->setInitialGuessNonzero(false);
    d_velocity_solver->solveSystem(*U_vec, *F_U_vec);

    // (2) Solve the pressure sub-problem.
    //
    // We treat two cases:
    //
    // (i) rho/dt = 0.
    //
    // In this case,
    //
    //    U - U^* + G Phi = 0
    //    -D U = F_P
    //
    // so that
    //
    //    Phi := inv(-L_p) * F_Phi = inv(-L_p) * (-F_P - D U^*)
    //    P   := -K*mu*F_Phi
    //
    // in which L_p = D*G.
    //
    // (ii) rho/dt != 0.
    //
    // In this case,
    //
    //    rho (U - U^*) + G Phi = 0
    //    -D U = F_P
    //
    // so that
    //
    //    Phi := inv(-L_rho) * F_phi = inv(-L_rho) * (-F_P - D U^*)
    //    P   := (1/dt - K*mu*L_rho)*Phi = (1/dt) Phi - K*mu*F_phi
    //
    // in which L_rho = D*(1/rho)*G.
    //
    // Approximate Poisson solvers are used in both cases.
    d_hier_math_ops->div(d_F_Phi_idx,
                         d_F_Phi_var,
                         -1.0,
                         U_idx,
                         U_sc_var,
                         d_no_fill_op,
                         d_new_time,
                         /*cf_bdry_synch*/ true,
                         -1.0,
                         F_P_idx,
                         F_P_cc_var);
    d_pressure_solver->setHomogeneousBc(true);
    LinearSolver* p_pressure_solver = dynamic_cast<LinearSolver*>(d_pressure_solver.getPointer());
    p_pressure_solver->setInitialGuessNonzero(false);
    d_pressure_solver->solveSystem(*Phi_scratch_vec, *F_Phi_vec);
    if (steady_state)
    {
        d_pressure_data_ops->scale(P_idx, -d_U_problem_coefs.getDConstant(), d_F_Phi_idx);
    }
    else
    {
        d_pressure_data_ops->linearSum(
            P_idx, 1.0 / getDt(), d_Phi_scratch_idx, -d_U_problem_coefs.getDConstant(), d_F_Phi_idx);
    }

    // (3) Evaluate U in terms of U^* and Phi.
    //
    // We treat two cases:
    //
    // (i) rho = 0.  In this case,
    //
    //    U = U^* - G Phi
    //
    // (ii) rho != 0.  In this case,
    //
    //    U = U^* - (1.0/rho) G Phi
    double coef;
    if (steady_state)
    {
        coef = -1.0;
    }
    else
    {
        coef = d_P_problem_coefs.getDConstant();
    }
    d_hier_math_ops->grad(U_idx,
                          U_sc_var,
                          /*cf_bdry_synch*/ true,
                          coef,
                          d_Phi_scratch_idx,
                          d_Phi_var,
                          d_Phi_bdry_fill_op,
                          d_pressure_solver->getSolutionTime(),
                          1.0,
                          U_idx,
                          U_sc_var);

    // Account for nullspace vectors.
    correctNullspace(U_vec, P_vec);

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_Phi_scratch_idx);
        level->deallocatePatchData(d_F_Phi_idx);
    }

    // Deallocate the solver (if necessary).
    if (deallocate_at_completion) deallocateSolverState();

    IBAMR_TIMER_STOP(t_solve_system);
    return true;
} // solveSystem

void
StaggeredStokesProjectionPreconditioner::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                                               const SAMRAIVectorReal<NDIM, double>& b)
{
    IBAMR_TIMER_START(t_initialize_solver_state);

    if (d_is_initialized) deallocateSolverState();

    // Parent class initialization.
    StaggeredStokesBlockPreconditioner::initializeSolverState(x, b);

    // Setup hierarchy operators.
    Pointer<VariableFillPattern<NDIM> > fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent P_scratch_component(d_Phi_scratch_idx,
                                                          DATA_REFINE_TYPE,
                                                          USE_CF_INTERPOLATION,
                                                          DATA_COARSEN_TYPE,
                                                          BDRY_EXTRAP_TYPE,
                                                          CONSISTENT_TYPE_2_BDRY,
                                                          d_P_bc_coef,
                                                          fill_pattern);
    d_Phi_bdry_fill_op = new HierarchyGhostCellInterpolation();
    d_Phi_bdry_fill_op->setHomogeneousBc(true);
    d_Phi_bdry_fill_op->initializeOperatorState(P_scratch_component, d_hierarchy);

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_solver_state);
    return;
} // initializeSolverState

void
StaggeredStokesProjectionPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_solver_state);

    // Parent class deallocation.
    StaggeredStokesBlockPreconditioner::deallocateSolverState();

    // Deallocate hierarchy operators.
    d_Phi_bdry_fill_op.setNull();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

void
StaggeredStokesProjectionPreconditioner::setInitialGuessNonzero(bool initial_guess_nonzero)
{
    if (initial_guess_nonzero)
    {
        TBOX_ERROR(d_object_name + "::setInitialGuessNonzero()\n"
                   << "  class IBAMR::StaggeredStokesProjectionPreconditioner requires a zero "
                      "initial guess"
                   << std::endl);
    }
    return;
} // setInitialGuessNonzero

void
StaggeredStokesProjectionPreconditioner::setMaxIterations(int max_iterations)
{
    if (max_iterations != 1)
    {
        TBOX_ERROR(d_object_name + "::setMaxIterations()\n"
                   << "  class IBAMR::StaggeredStokesProjectionPreconditioner only performs a "
                      "single iteration"
                   << std::endl);
    }
    return;
} // setMaxIterations

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
