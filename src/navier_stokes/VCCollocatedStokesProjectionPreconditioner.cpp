// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
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

#include "ibamr/VCCollocatedStokesProjectionPreconditioner.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PoissonSolver.h"

#include "CellVariable.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VariableFillPattern.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include <ostream>
#include <string>

#include "ibamr/namespaces.h" // IWYU pragma: keep

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
static const std::string DATA_REFINE_TYPE = "LINEAR_REFINE";
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
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

VCCollocatedStokesProjectionPreconditioner::VCCollocatedStokesProjectionPreconditioner(const std::string& object_name,
                                                                                       Pointer<Database> input_db)
    : LinearSolver(),
      d_Phi_bdry_fill_op(nullptr),
      d_Phi_var(nullptr),
      d_F_Phi_var(nullptr),
      d_Phi_scratch_idx(-1),
      d_F_Phi_idx(-1)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ true);

    if (input_db)
    {
        d_eps = input_db->getDoubleWithDefault("pressure_stabilization_coef", d_eps);
        d_alpha = input_db->getDoubleWithDefault("velocity_correction_coef", d_alpha);
    }

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
                      "IBAMR::VCCollocatedStokesProjectionPreconditioner::solveSystem()");
                  t_initialize_solver_state =
                      TimerManager::getManager()->getTimer("IBAMR::VCCollocatedStokesProjectionPreconditioner::"
                                                           "initializeSolverState()");
                  t_deallocate_solver_state =
                      TimerManager::getManager()->getTimer("IBAMR::VCCollocatedStokesProjectionPreconditioner::"
                                                           "deallocateSolverState()"););
    return;
} // VCCollocatedStokesProjectionPreconditioner

VCCollocatedStokesProjectionPreconditioner::~VCCollocatedStokesProjectionPreconditioner()
{
    deallocateSolverState();
    return;
} // ~VCCollocatedStokesProjectionPreconditioner

void
VCCollocatedStokesProjectionPreconditioner::setVelocitySubdomainSolver(Pointer<PoissonSolver> velocity_solver)
{
    d_velocity_solver = velocity_solver;
    return;
} // setVelocitySubdomainSolver

void
VCCollocatedStokesProjectionPreconditioner::setPressureSubdomainSolver(Pointer<PoissonSolver> pressure_solver)
{
    d_pressure_solver = pressure_solver;
    return;
} // setPressureSubdomainSolver

void
VCCollocatedStokesProjectionPreconditioner::setPressureStabilizationCoeffient(double eps)
{
    d_eps = eps;
    return;
} // setPressureStabilizationCoeffient

void
VCCollocatedStokesProjectionPreconditioner::setVelocityCorrectionCoeffient(double alpha)
{
    d_alpha = alpha;
    return;
} // setVelocityCorrectionCoeffient

void
VCCollocatedStokesProjectionPreconditioner::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
    RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(U_bc_coefs.size() == NDIM);
#endif
    d_U_bc_coefs = U_bc_coefs;
    d_P_bc_coef = P_bc_coef;
    if (d_velocity_solver) d_velocity_solver->setPhysicalBcCoefs(d_U_bc_coefs);
    if (d_pressure_solver) d_pressure_solver->setPhysicalBcCoef(d_P_bc_coef);
    return;
} // setPhysicalBcCoefs

bool
VCCollocatedStokesProjectionPreconditioner::solveSystem(SAMRAIVectorReal<NDIM, double>& x,
                                                        SAMRAIVectorReal<NDIM, double>& b)
{
    IBAMR_TIMER_START(t_solve_system);

    // Initialize the solver (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeSolverState(x, b);

    // Get the vector components.
    const int F_U_idx = b.getComponentDescriptorIndex(0);
    const int F_P_idx = b.getComponentDescriptorIndex(1);

    const Pointer<Variable<NDIM> >& F_U_var = b.getComponentVariable(0);
    const Pointer<Variable<NDIM> >& F_P_var = b.getComponentVariable(1);

    Pointer<CellVariable<NDIM, double> > F_U_cc_var = F_U_var;
    Pointer<CellVariable<NDIM, double> > F_P_cc_var = F_P_var;

    const int U_idx = x.getComponentDescriptorIndex(0);
    const int P_idx = x.getComponentDescriptorIndex(1);

    const Pointer<Variable<NDIM> >& U_var = x.getComponentVariable(0);
    const Pointer<Variable<NDIM> >& P_var = x.getComponentVariable(1);

    Pointer<CellVariable<NDIM, double> > U_cc_var = U_var;
    Pointer<CellVariable<NDIM, double> > P_cc_var = P_var;

    // Setup the component solver vectors.
    Pointer<SAMRAIVectorReal<NDIM, double> > F_U_vec;
    F_U_vec = new SAMRAIVectorReal<NDIM, double>(d_object_name + "::F_U", d_hierarchy, d_coarsest_ln, d_finest_ln);
    F_U_vec->addComponent(F_U_cc_var, F_U_idx, d_velocity_wgt_idx, d_velocity_data_ops);

    Pointer<SAMRAIVectorReal<NDIM, double> > U_vec;
    U_vec = new SAMRAIVectorReal<NDIM, double>(d_object_name + "::U", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_vec->addComponent(U_cc_var, U_idx, d_velocity_wgt_idx, d_velocity_data_ops);

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
    // U^* := A^-1 (F_U)
    //
    // An approximate Helmholtz solver is used.
    d_velocity_solver->setHomogeneousBc(true);
    auto p_velocity_solver = dynamic_cast<LinearSolver*>(d_velocity_solver.getPointer());
    if (p_velocity_solver) p_velocity_solver->setInitialGuessNonzero(false);
    d_velocity_solver->solveSystem(*U_vec, *F_U_vec);

    // (2) Solve the pressure sub-problem.
    //
    // In this case,
    //
    //    U - U^* + alpha G Phi = 0
    //    -D U + eps Lp Phi = F_P

    //    D U - DU^* + alpha DG Phi = eps Lp (Phi) - F_p - DU^* + alpha Lp Phi = 0
    // => (alpha + eps)Lp Phi = F_p + DU^*
    // so that
    //
    //    Phi := inv(-Lp) * (-F_P - D U^*) / (alpha + eps)
    //    P   := Phi
    //
    //  Note that here we allow for a generalized version of the divergence
    //  operator, i.e., D(.) = div (coef .). In this case Lp should be interpreted
    // as Lp = D (coef.) G
    //
    const double fac = 1.0 / (d_alpha + d_eps);
    d_hier_math_ops->div(
        d_F_Phi_idx, d_F_Phi_var, -fac, U_idx, U_cc_var, d_no_fill_op, d_new_time, -fac, F_P_idx, F_P_cc_var);

    d_pressure_solver->setHomogeneousBc(true);
    auto p_pressure_solver = dynamic_cast<LinearSolver*>(d_pressure_solver.getPointer());
    p_pressure_solver->setInitialGuessNonzero(false);
    d_pressure_solver->solveSystem(*Phi_scratch_vec, *F_Phi_vec);

    // P = Phi
    d_pressure_data_ops->copyData(P_idx, d_Phi_scratch_idx);

    // (3) Evaluate U in terms of U^* and Phi.

    //
    //    U = U^* - alpha G Phi
    //
    d_hier_math_ops->grad(U_idx,
                          U_cc_var,
                          -d_alpha,
                          d_Phi_scratch_idx,
                          d_Phi_var,
                          d_Phi_bdry_fill_op,
                          d_pressure_solver->getSolutionTime(),
                          1.0,
                          U_idx,
                          U_cc_var);

    // Account for nullspace vectors.
    correctNullSpace(U_vec, P_vec);

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
VCCollocatedStokesProjectionPreconditioner::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                                                  const SAMRAIVectorReal<NDIM, double>& b)
{
    IBAMR_TIMER_START(t_initialize_solver_state);

    if (d_is_initialized) deallocateSolverState();

    // Get the hierarchy configuration.
    d_hierarchy = x.getPatchHierarchy();
    d_coarsest_ln = x.getCoarsestLevelNumber();
    d_finest_ln = x.getFinestLevelNumber();
#if !defined(NDEBUG)
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

    d_hier_math_ops =
        new HierarchyMathOps(d_object_name + "::HierarchyMathOps", d_hierarchy, d_coarsest_ln, d_finest_ln);

    // Setup hierarchy operators.
    Pointer<VariableFillPattern<NDIM> > fill_pattern = nullptr;
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
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
VCCollocatedStokesProjectionPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_solver_state);

    d_velocity_data_ops.setNull();
    d_pressure_data_ops.setNull();
    d_velocity_wgt_idx = invalid_index;
    d_pressure_wgt_idx = invalid_index;
    d_hier_math_ops.setNull();
    d_Phi_bdry_fill_op.setNull();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

void
VCCollocatedStokesProjectionPreconditioner::setInitialGuessNonzero(bool initial_guess_nonzero)
{
    if (initial_guess_nonzero)
    {
        TBOX_ERROR(d_object_name + "::setInitialGuessNonzero()\n"
                   << "  class IBAMR::VCCollocatedStokesProjectionPreconditioner "
                      "requires a zero "
                      "initial guess"
                   << std::endl);
    }
    return;
} // setInitialGuessNonzero

void
VCCollocatedStokesProjectionPreconditioner::setMaxIterations(int max_iterations)
{
    if (max_iterations != 1)
    {
        TBOX_ERROR(d_object_name + "::setMaxIterations()\n"
                   << "  class IBAMR::VCCollocatedStokesProjectionPreconditioner "
                      "only performs a "
                      "single iteration"
                   << std::endl);
    }
    return;
} // setMaxIterations

/////////////////////////////// PROTECTED ////////////////////////////////////

void
VCCollocatedStokesProjectionPreconditioner::correctNullSpace(Pointer<SAMRAIVectorReal<NDIM, double> > U_vec,
                                                             Pointer<SAMRAIVectorReal<NDIM, double> > P_vec)
{
    auto p_velocity_solver = dynamic_cast<LinearSolver*>(d_velocity_solver.getPointer());
    if (p_velocity_solver)
    {
        const std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > >& U_nul_vecs =
            p_velocity_solver->getNullSpaceBasisVectors();
        if (!U_nul_vecs.empty())
        {
            for (const auto& U_nul_vec : U_nul_vecs)
            {
                const double alpha = U_vec->dot(U_nul_vec) / U_nul_vec->dot(U_nul_vec);
                U_vec->axpy(-alpha, U_nul_vec, U_vec);
            }
        }
#if !defined(NDEBUG)
        TBOX_ASSERT(!p_velocity_solver->getNullSpaceContainsConstantVector());
#endif
    }

    auto p_pressure_solver = dynamic_cast<LinearSolver*>(d_pressure_solver.getPointer());
    if (p_pressure_solver)
    {
        if (p_pressure_solver->getNullSpaceContainsConstantVector())
        {
            const int P_idx = P_vec->getComponentDescriptorIndex(0);
            const double volume = d_hier_math_ops->getVolumeOfPhysicalDomain();
            const double P_mean = (1.0 / volume) * d_pressure_data_ops->integral(P_idx, d_pressure_wgt_idx);
            d_pressure_data_ops->addScalar(P_idx, P_idx, -P_mean);
        }
#if !defined(NDEBUG)
        TBOX_ASSERT(p_pressure_solver->getNullSpaceBasisVectors().empty());
#endif
    }
    return;
} // correctNullSpace

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
