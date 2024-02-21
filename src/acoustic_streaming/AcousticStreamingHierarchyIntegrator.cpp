// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#include "ibamr/AcousticStreamingHierarchyIntegrator.h"
#include "ibamr/BrinkmanPenalizationStrategy.h"
#include "ibamr/FOAcousticStreamingPETScLevelSolver.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/INSProjectionBcCoef.h"
#include "ibamr/INSVCStaggeredPressureBcCoef.h"
#include "ibamr/INSVCStaggeredVelocityBcCoef.h"
#include "ibamr/PETScKrylovStaggeredStokesSolver.h"
#include "ibamr/VCStaggeredStokesOperator.h"
#include "ibamr/VCStaggeredStokesProjectionPreconditioner.h"

#include "ibtk/CartGridFunctionSet.h"
#include "ibtk/CartSideDoubleDivPreservingRefine.h"
#include "ibtk/CartSideDoubleRT0Refine.h"
#include "ibtk/CartSideDoubleSpecializedLinearRefine.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/VCSCViscousOpPointRelaxationFACOperator.h"
#include "ibtk/VCSCViscousOperator.h"
#include "ibtk/VCSCViscousPETScLevelSolver.h"

#include "HierarchyDataOpsManager.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

// FORTRAN ROUTINES

extern "C"
{
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{

// Version of AcousticStreamingHierarchyIntegrator restart file data.
static const int ACOUSTIC_STREAMING_HIERARCHY_INTEGRATOR_VERSION = 1;

// Timers.
static Timer* t_setup_plot_data_specialized;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;
static const int NODEG = 1;
static const int EDGEG = 1;
static const int MUCELLG = 2;

// Real and imaginary component naming
static const int REAL = 0;
static const int IMAG = 1;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Default solvers and preconditioners to register with the solver managers
static const std::string DEFAULT_VC_STAGGERED_STOKES_SOLVER = "VC_STAGGERED_STOKES_PETSC_KRYLOV_SOLVER";
static const std::string DEFAULT_VC_STAGGERED_STOKES_PRECOND = "VC_STAGGERED_STOKES_PROJECTION_PRECONDITIONER";
static const std::string DEFAULT_VC_VELOCITY_SOLVER = "VC_VELOCITY_PETSC_KRYLOV_SOLVER";
static const std::string DEFAULT_VC_VELOCITY_PRECOND = "VC_VELOCITY_POINT_RELAXATION_FAC_PRECONDITIONER";
static const std::string DEFAULT_VC_VELOCITY_LEVEL_SOLVER = "VC_VELOCITY_PETSC_LEVEL_SOLVER";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

Pointer<StaggeredStokesSolver>
allocate_vc_stokes_krylov_solver(const std::string& solver_object_name,
                                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
                                 const std::string& solver_default_options_prefix)
{
    Pointer<PETScKrylovStaggeredStokesSolver> krylov_solver =
        new PETScKrylovStaggeredStokesSolver(solver_object_name, solver_input_db, solver_default_options_prefix);
    krylov_solver->setOperator(new VCStaggeredStokesOperator(solver_object_name + "::vc_staggered_stokes_operator"));
    return krylov_solver;
}

Pointer<PoissonSolver>
allocate_vc_velocity_krylov_solver(const std::string& solver_object_name,
                                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
                                   const std::string& solver_default_options_prefix)
{
    Pointer<PETScKrylovPoissonSolver> krylov_solver =
        new PETScKrylovPoissonSolver(solver_object_name, solver_input_db, solver_default_options_prefix);
    krylov_solver->setOperator(new VCSCViscousOperator(solver_object_name + "::vc_viscous_operator"));
    return krylov_solver;
}

} // namespace

AcousticStreamingHierarchyIntegrator::AcousticStreamingHierarchyIntegrator(std::string object_name,
                                                                           Pointer<Database> input_db,
                                                                           bool register_for_restart)
    : HierarchyIntegrator(std::move(object_name), input_db, register_for_restart),
      d_default_so_bc_coefs(d_object_name + "::default_so_bc_coefs", Pointer<Database>(nullptr)),
      d_so_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(nullptr)),
      d_vc_stokes_op_spec(),
      d_vc_projection_pc_spec()
{
    d_U1_var = new SideVariable<NDIM, double>(object_name + "::U1", /*depth*/ 2);
    d_U2_var = new SideVariable<NDIM, double>(object_name + "::U2", /*depth*/ 1);

    d_P1_var = new CellVariable<NDIM, double>(object_name + "::P1", /*depth*/ 2);
    d_P2_var = new CellVariable<NDIM, double>(object_name + "::P2", /*depth*/ 1);

    d_F1_var = new SideVariable<NDIM, double>(object_name + "::F1", /*depth*/ 2);
    d_F2_var = new SideVariable<NDIM, double>(object_name + "::F2", /*depth*/ 1);

    d_Q1_var = new CellVariable<NDIM, double>(object_name + "::Q1", /*depth*/ 2);
    d_Q2_var = new CellVariable<NDIM, double>(object_name + "::Q2", /*depth*/ 1);

    d_U1_plot_var = new CellVariable<NDIM, double>(object_name + "::U1_plot", /*depth*/ 2 * NDIM);
    d_U2_plot_var = new CellVariable<NDIM, double>(object_name + "::U2_plot", /*depth*/ NDIM);

    d_rho_plot_var = new CellVariable<NDIM, double>(object_name + "::rho_plot", /*depth*/ NDIM);

#if (NDIM == 2)
    d_Omega2_var = new CellVariable<NDIM, double>(d_object_name + "::Omega2");
#endif
#if (NDIM == 3)
    d_Omega2_var = new CellVariable<NDIM, double>(d_object_name + "::Omega2", NDIM);
#endif
#if (NDIM == 3)
    d_Omega2_Norm_var = new CellVariable<NDIM, double>(d_object_name + "::|Omega2|_2");
#endif

    d_U2_regrid_var = new SideVariable<NDIM, double>(d_object_name + "::U2_regrid");
    d_U2_src_var = new SideVariable<NDIM, double>(d_object_name + "::U2_src");
    d_indicator2_var = new SideVariable<NDIM, double>(d_object_name + "::indicator2");

    auto set_timer = [&](const char* name) { return tbox::TimerManager::getManager()->getTimer(name); };
    IBTK_DO_ONCE(t_setup_plot_data_specialized =
                     set_timer("IBAMR::AcousticStreamingHierarchyIntegrator::setupPlotDataSpecialized()"););

    // Setup physical boundary conditions objects for the first-order system.
    for (int comp = 0; comp < 2; ++comp)
    {
        d_U1_bc_coefs[comp].resize(NDIM, nullptr);
    }

    // Setup physical boundary conditions objects for the second-order system.

    // First, setup default boundary condition objects that specify homogeneous
    // Dirichlet (solid-wall) boundary conditions for the velocity.
    // These BCs will be replaced further downstream during the problem setup (e.g., in main.cpp)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_so_bc_coefs.setBoundaryValue(2 * d, 0.0);
        d_default_so_bc_coefs.setBoundaryValue(2 * d + 1, 0.0);
    }
    registerSecondOrderPhysicalBoundaryConditions(
        std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, &d_default_so_bc_coefs));

    d_bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    d_U2_bc_coefs.resize(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_U2_bc_coefs[d] = new INSVCStaggeredVelocityBcCoef(d, this, d_so_bc_coefs, IBAMR::TRACTION);
    }
    d_P2_bc_coef = new INSVCStaggeredPressureBcCoef(this, d_so_bc_coefs, IBAMR::TRACTION);

    d_U2_star_bc_coefs.resize(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_U2_star_bc_coefs[d] = new INSIntermediateVelocityBcCoef(d, d_so_bc_coefs);
    }
    d_Phi_bc_coef = new INSProjectionBcCoef(d_so_bc_coefs);

    // Register solver factory functions for variable coefficient Stokes and
    // viscous solvers for the second order system.
    StaggeredStokesSolverManager::getManager()->registerSolverFactoryFunction(DEFAULT_VC_STAGGERED_STOKES_SOLVER,
                                                                              allocate_vc_stokes_krylov_solver);
    StaggeredStokesSolverManager::getManager()->registerSolverFactoryFunction(
        DEFAULT_VC_STAGGERED_STOKES_PRECOND, VCStaggeredStokesProjectionPreconditioner::allocate_solver);
    SCPoissonSolverManager::getManager()->registerSolverFactoryFunction(DEFAULT_VC_VELOCITY_SOLVER,
                                                                        allocate_vc_velocity_krylov_solver);
    SCPoissonSolverManager::getManager()->registerSolverFactoryFunction(
        DEFAULT_VC_VELOCITY_PRECOND, VCSCViscousOpPointRelaxationFACOperator::allocate_solver);
    SCPoissonSolverManager::getManager()->registerSolverFactoryFunction(DEFAULT_VC_VELOCITY_LEVEL_SOLVER,
                                                                        VCSCViscousPETScLevelSolver::allocate_solver);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    return;
} // AcousticStreamingHierarchyIntegrator

AcousticStreamingHierarchyIntegrator::~AcousticStreamingHierarchyIntegrator()
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        delete d_U2_bc_coefs[d];
        delete d_U2_star_bc_coefs[d];

        d_U2_bc_coefs[d] = nullptr;
        d_U2_star_bc_coefs[d] = nullptr;
    }
    delete d_Phi_bc_coef;
    d_Phi_bc_coef = nullptr;

    return;
} // ~AcousticStreamingHierarchyIntegrator

void
AcousticStreamingHierarchyIntegrator::registerFirstOrderPhysicalBoundaryConditions(
    const std::array<std::vector<RobinBcCoefStrategy<NDIM>*>, 2>& U_bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(U_bc_coefs[REAL].size() == NDIM);
    TBOX_ASSERT(U_bc_coefs[IMAG].size() == NDIM);
#endif

    d_U1_bc_coefs = U_bc_coefs;
    return;
} // registerFirstOrderPhysical BoundaryConditions

void
AcousticStreamingHierarchyIntegrator::registerSecondOrderPhysicalBoundaryConditions(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    d_so_bc_coefs = bc_coefs;
    return;
} // registerSecondOrderPhysicalBoundaryConditions

void
AcousticStreamingHierarchyIntegrator::registerFirstOrderVelocityInitialConditions(Pointer<CartGridFunction> U1_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_U1_init = U1_init;
    return;
} // registerFirstOrderVelocityInitialConditions

void
AcousticStreamingHierarchyIntegrator::registerFirstOrderPressureInitialConditions(Pointer<CartGridFunction> P1_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_P1_init = P1_init;
    return;
} // registerFirstOrderPressureInitialConditions

void
AcousticStreamingHierarchyIntegrator::registerSecondOrderVelocityInitialConditions(Pointer<CartGridFunction> U2_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_U2_init = U2_init;
    return;
} // registerSecondOrderVelocityInitialConditions

void
AcousticStreamingHierarchyIntegrator::registerSecondOrderPressureInitialConditions(Pointer<CartGridFunction> P2_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_P2_init = P2_init;
    return;
} // registerSecondOrderPressureInitialConditions

void
AcousticStreamingHierarchyIntegrator::registerFirstOrderBodyForceFunction(Pointer<CartGridFunction> F1_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    if (d_F1_fcn)
    {
        Pointer<CartGridFunctionSet> p_F1_fcn = d_F1_fcn;
        if (!p_F1_fcn)
        {
            pout << d_object_name << "::registerFirstOrderBodyForceFunction(): WARNING:\n"
                 << "  body force function has already been set.\n"
                 << "  functions will be evaluated in the order in which they were "
                    "registered "
                    "with "
                    "the solver\n"
                 << "  when evaluating the body force term value.\n";
            p_F1_fcn = new CartGridFunctionSet(d_object_name + "::fo_body_force_function_set");
            p_F1_fcn->addFunction(d_F1_fcn);
        }
        p_F1_fcn->addFunction(F1_fcn);
        d_F1_fcn = p_F1_fcn;
    }
    else
    {
        d_F1_fcn = F1_fcn;
    }
    return;
} // registerFirstOrderBodyForceFunction

void
AcousticStreamingHierarchyIntegrator::registerSecondOrderBodyForceFunction(Pointer<CartGridFunction> F2_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    if (d_F2_fcn)
    {
        Pointer<CartGridFunctionSet> p_F2_fcn = d_F2_fcn;
        if (!p_F2_fcn)
        {
            pout << d_object_name << "::registerSecondOrderBodyForceFunction(): WARNING:\n"
                 << "  body force function has already been set.\n"
                 << "  functions will be evaluated in the order in which they were "
                    "registered "
                    "with "
                    "the solver\n"
                 << "  when evaluating the body force term value.\n";
            p_F2_fcn = new CartGridFunctionSet(d_object_name + "::so_body_force_function_set");
            p_F2_fcn->addFunction(d_F2_fcn);
        }
        p_F2_fcn->addFunction(F2_fcn);
        d_F2_fcn = p_F2_fcn;
    }
    else
    {
        d_F2_fcn = F2_fcn;
    }
    return;
} // registerSecondOrderBodyForceFunction

void
AcousticStreamingHierarchyIntegrator::registerFirstOrderVelocityDivergenceFunction(Pointer<CartGridFunction> Q1_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    if (d_Q1_fcn)
    {
        Pointer<CartGridFunctionSet> p_Q1_fcn = d_Q1_fcn;
        if (!p_Q1_fcn)
        {
            pout << d_object_name << "::registerFirstOrderVelocityDivergenceFunction(): WARNING:\n"
                 << "  velocity divergence function has already been set.\n"
                 << "  functions will be evaluated in the order in which they were "
                    "registered with the solver\n"
                 << "  when evaluating the fluid source term value.\n";
            p_Q1_fcn = new CartGridFunctionSet(d_object_name + "::fo_velocity_divergence_function_set");
            p_Q1_fcn->addFunction(d_Q1_fcn);
        }
        p_Q1_fcn->addFunction(Q1_fcn);
    }
    else
    {
        d_Q1_fcn = Q1_fcn;
    }
    return;
} // registerFirstOrderVelocityDivergenceFunction

void
AcousticStreamingHierarchyIntegrator::registerSecondOrderVelocityDivergenceFunction(Pointer<CartGridFunction> Q2_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    if (d_Q2_fcn)
    {
        Pointer<CartGridFunctionSet> p_Q2_fcn = d_Q2_fcn;
        if (!p_Q2_fcn)
        {
            pout << d_object_name << "::registerSecondOrderVelocityDivergenceFunction(): WARNING:\n"
                 << "  velocity divergence function has already been set.\n"
                 << "  functions will be evaluated in the order in which they were "
                    "registered with the solver\n"
                 << "  when evaluating the fluid source term value.\n";
            p_Q2_fcn = new CartGridFunctionSet(d_object_name + "::so_velocity_divergence_function_set");
            p_Q2_fcn->addFunction(d_Q2_fcn);
        }
        p_Q2_fcn->addFunction(Q2_fcn);
    }
    else
    {
        d_Q2_fcn = Q2_fcn;
    }
    return;
} // registerSecondOrderVelocityDivergenceFunction

Pointer<PoissonSolver>
AcousticStreamingHierarchyIntegrator::getSecondOrderVelocitySubdomainSolver()
{
    return d_velocity_solver;
} // getSecondOrderVelocitySubdomainSolver

Pointer<PoissonSolver>
AcousticStreamingHierarchyIntegrator::getSecondOrderPressureSubdomainSolver()
{
    return d_pressure_solver;
} // getSecondOrderPressureSubdomainSolver

void
AcousticStreamingHierarchyIntegrator::registerMassDensityVariable(Pointer<Variable<NDIM> > rho_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_rho_var);
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_rho_var = rho_var;
    return;
} // registerMassDensityVariable

Pointer<Variable<NDIM> >
AcousticStreamingHierarchyIntegrator::getMassDensityVariable() const
{
    return d_rho_var;
} // getMassDensityVariable

void
AcousticStreamingHierarchyIntegrator::registerShearViscosityVariable(Pointer<Variable<NDIM> > mu_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_mu_var);
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_mu_var = mu_var;
    return;
} // registerShearViscosityVariable

Pointer<Variable<NDIM> >
AcousticStreamingHierarchyIntegrator::getShearViscosityVariable() const
{
    return d_mu_var;
} // getShearViscosityVariable

void
AcousticStreamingHierarchyIntegrator::registerBulkViscosityVariable(Pointer<Variable<NDIM> > lambda_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_lambda_var);
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_lambda_var = lambda_var;
    return;
} // registerBulkViscosityVariable

Pointer<Variable<NDIM> >
AcousticStreamingHierarchyIntegrator::getBulkViscosityVariable() const
{
    return d_lambda_var;
} // getBulkViscosityVariable

void
AcousticStreamingHierarchyIntegrator::setDensityVCInterpolationType(const IBTK::VCInterpType vc_interp_type)
{
    d_rho_vc_interp_type = vc_interp_type;
    return;
} // setDensityVCInterpolationType

void
AcousticStreamingHierarchyIntegrator::setShearViscosityVCInterpolationType(const IBTK::VCInterpType vc_interp_type)
{
    d_mu_vc_interp_type = vc_interp_type;
    return;
} // setShearViscosityVCInterpolationType

void
AcousticStreamingHierarchyIntegrator::registerResetFluidDensityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx)
{
    d_reset_rho_fcns.push_back(callback);
    d_reset_rho_fcns_ctx.push_back(ctx);
    return;
} // registerResetFluidDensityFcn

void
AcousticStreamingHierarchyIntegrator::registerResetFluidShearViscosityFcn(ResetFluidPropertiesFcnPtr callback,
                                                                          void* ctx)
{
    d_reset_mu_fcns.push_back(callback);
    d_reset_mu_fcns_ctx.push_back(ctx);
    return;
} // registerResetFluidShearViscosityFcn

void
AcousticStreamingHierarchyIntegrator::registerResetFluidBulkViscosityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx)
{
    d_reset_lambda_fcns.push_back(callback);
    d_reset_lambda_fcns_ctx.push_back(ctx);
    return;
} // registerResetFluidBulkViscosityFcn

void
AcousticStreamingHierarchyIntegrator::registerMassDensityInitialConditions(const Pointer<CartGridFunction> rho_init_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_rho_init_fcn = rho_init_fcn;
    return;
} // registerMassDensityInitialConditions

void
AcousticStreamingHierarchyIntegrator::registerShearViscosityInitialConditions(
    const Pointer<CartGridFunction> mu_init_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_mu_init_fcn = mu_init_fcn;
    return;
} // registerShearViscosityInitialConditions

void
AcousticStreamingHierarchyIntegrator::registerBulkViscosityInitialConditions(
    const Pointer<CartGridFunction> lambda_init_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_lambda_init_fcn = lambda_init_fcn;
    return;
} // registerBulkViscosityInitialConditions

void
AcousticStreamingHierarchyIntegrator::registerMassDensityBoundaryConditions(RobinBcCoefStrategy<NDIM>* rho_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(rho_bc_coef);
#endif
    d_rho_bc_coef = rho_bc_coef;
    return;
} // registerMassDensityBoundaryConditions

void
AcousticStreamingHierarchyIntegrator::registerShearViscosityBoundaryConditions(
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* mu_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(mu_bc_coef);
#endif
    d_mu_bc_coef = mu_bc_coef;
    return;
} // registerShearViscosityBoundaryCondition

void
AcousticStreamingHierarchyIntegrator::registerBulkViscosityBoundaryConditions(
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* lambda_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(lambda_bc_coef);
#endif
    d_lambda_bc_coef = lambda_bc_coef;
    return;
} // registerBulkViscosityBoundaryConditions

void
AcousticStreamingHierarchyIntegrator::registerBrinkmanPenalizationStrategy(
    Pointer<BrinkmanPenalizationStrategy> brinkman_force)
{
    d_brinkman_force.push_back(brinkman_force);
    return;
} // registerBrinkmanPenalizationStrategy

void
AcousticStreamingHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                                    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;
    const int max_levels = gridding_alg->getMaxLevels();

    // Setup solver for the first-order system
    d_first_order_solver = new FOAcousticStreamingPETScLevelSolver(
        d_object_name + "::first_order_solver", d_first_order_solver_db, "fo_solver_");
    d_first_order_solver_needs_init = true;

    // Setup solvers for the second order system.
    d_stokes_solver_type = DEFAULT_VC_STAGGERED_STOKES_SOLVER;
    d_stokes_precond_type = DEFAULT_VC_STAGGERED_STOKES_PRECOND;
    d_stokes_solver = StaggeredStokesSolverManager::getManager()->allocateSolver(d_stokes_solver_type,
                                                                                 d_object_name + "::so_stokes_solver",
                                                                                 d_stokes_solver_db,
                                                                                 "so_solver_",
                                                                                 d_stokes_precond_type,
                                                                                 d_object_name + "::stokes_precond",
                                                                                 d_stokes_precond_db,
                                                                                 "so_solver_pc_");
    d_stokes_solver_needs_init = true;

    d_velocity_solver_type = DEFAULT_VC_VELOCITY_SOLVER;
    if (max_levels == 1)
    {
        d_velocity_precond_type = DEFAULT_VC_VELOCITY_LEVEL_SOLVER;
    }
    else
    {
        d_velocity_precond_type = DEFAULT_VC_VELOCITY_PRECOND;
    }
    d_velocity_solver = SCPoissonSolverManager::getManager()->allocateSolver(d_velocity_solver_type,
                                                                             d_object_name + "::so_velocity_solver",
                                                                             d_velocity_solver_db,
                                                                             "so_velocity_",
                                                                             d_velocity_precond_type,
                                                                             d_object_name + "::so_velocity_precond",
                                                                             d_velocity_precond_db,
                                                                             "so_velocity_pc_");
    d_velocity_solver_needs_init = true;

    d_pressure_solver_type = CCPoissonSolverManager::PETSC_KRYLOV_SOLVER;
    if (max_levels == 1)
    {
        d_pressure_precond_type = CCPoissonSolverManager::DEFAULT_LEVEL_SOLVER;
    }
    else
    {
        d_pressure_precond_type = CCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
    }
    d_pressure_solver = CCPoissonSolverManager::getManager()->allocateSolver(d_pressure_solver_type,
                                                                             d_object_name + "::so_pressure_solver",
                                                                             d_pressure_solver_db,
                                                                             "so_pressure_",
                                                                             d_pressure_precond_type,
                                                                             d_object_name + "::so_pressure_precond",
                                                                             d_pressure_precond_db,
                                                                             "so_pressure_pc_");
    d_pressure_solver_needs_init = true;

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_cc_data_ops =
        hier_ops_manager->getOperationsDouble(new CellVariable<NDIM, double>("cc_var"), hierarchy, true);
    d_hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(new SideVariable<NDIM, double>("sc_var"), hierarchy, true);
    d_hier_nc_data_ops =
        hier_ops_manager->getOperationsDouble(new NodeVariable<NDIM, double>("nc_var"), hierarchy, true);
    d_hier_ec_data_ops =
        hier_ops_manager->getOperationsDouble(new EdgeVariable<NDIM, double>("ec_var"), hierarchy, true);
    d_hier_math_ops = buildHierarchyMathOps(d_hierarchy);

    // Register state variables that are maintained by the
    // INSVCStaggeredHierarchyIntegrator.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    grid_geom->addSpatialRefineOperator(new CartSideDoubleRT0Refine());
    grid_geom->addSpatialRefineOperator(new CartSideDoubleSpecializedLinearRefine());

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> node_ghosts = NODEG;
    const IntVector<NDIM> edge_ghosts = EDGEG;
    const IntVector<NDIM> no_ghosts = 0;
    const IntVector<NDIM> mu_cell_ghosts = MUCELLG;
    const IntVector<NDIM> lambda_cell_ghosts = CELLG;

    registerVariable(d_U1_current_idx,
                     d_U1_new_idx,
                     d_U1_scratch_idx,
                     d_U1_var,
                     side_ghosts,
                     d_U_coarsen_type,
                     d_U_refine_type,
                     d_U1_init);

    registerVariable(d_P1_current_idx,
                     d_P1_new_idx,
                     d_P1_scratch_idx,
                     d_P1_var,
                     cell_ghosts,
                     d_P_coarsen_type,
                     d_P_refine_type,
                     d_P1_init);

    if (d_F1_fcn)
    {
        registerVariable(d_F1_current_idx,
                         d_F1_new_idx,
                         d_F1_scratch_idx,
                         d_F1_var,
                         side_ghosts,
                         d_F_coarsen_type,
                         d_F_refine_type,
                         d_F1_fcn);
    }
    else
    {
        d_F1_current_idx = invalid_index;
        d_F1_new_idx = invalid_index;
        d_F1_scratch_idx = invalid_index;
    }

    if (d_Q1_fcn)
    {
        registerVariable(d_Q1_current_idx,
                         d_Q1_new_idx,
                         d_Q1_scratch_idx,
                         d_Q1_var,
                         cell_ghosts,
                         d_Q_coarsen_type,
                         d_Q_refine_type,
                         d_Q1_fcn);
    }
    else
    {
        d_Q1_current_idx = invalid_index;
        d_Q1_new_idx = invalid_index;
        d_Q1_scratch_idx = invalid_index;
    }

    registerVariable(d_U2_current_idx,
                     d_U2_new_idx,
                     d_U2_scratch_idx,
                     d_U2_var,
                     side_ghosts,
                     d_U_coarsen_type,
                     d_U_refine_type,
                     d_U2_init);

    registerVariable(d_P2_current_idx,
                     d_P2_new_idx,
                     d_P2_scratch_idx,
                     d_P2_var,
                     cell_ghosts,
                     d_P_coarsen_type,
                     d_P_refine_type,
                     d_P2_init);

    if (d_F2_fcn)
    {
        registerVariable(d_F2_current_idx,
                         d_F2_new_idx,
                         d_F2_scratch_idx,
                         d_F2_var,
                         side_ghosts,
                         d_F_coarsen_type,
                         d_F_refine_type,
                         d_F2_fcn);
    }
    else
    {
        d_F2_current_idx = invalid_index;
        d_F2_new_idx = invalid_index;
        d_F2_scratch_idx = invalid_index;
    }

    if (d_Q2_fcn)
    {
        registerVariable(d_Q2_current_idx,
                         d_Q2_new_idx,
                         d_Q2_scratch_idx,
                         d_Q2_var,
                         cell_ghosts,
                         d_Q_coarsen_type,
                         d_Q_refine_type,
                         d_Q2_fcn);
    }
    else
    {
        d_Q2_current_idx = invalid_index;
        d_Q2_new_idx = invalid_index;
        d_Q2_scratch_idx = invalid_index;
    }

#if !defined(NDEBUG)
    // AcousticStreamingHierarchyIntegrator should initialize the shear viscosity variable.
    TBOX_ASSERT(d_mu_init_fcn || d_reset_mu_fcns.size() > 0);
#endif
    registerVariable(d_mu_current_idx,
                     d_mu_new_idx,
                     d_mu_scratch_idx,
                     d_mu_var,
                     mu_cell_ghosts,
                     d_mu_coarsen_type,
                     d_mu_refine_type,
                     d_mu_init_fcn);

    if (d_lambda_var)
    {
#if !defined(NDEBUG)
        // AcousticStreamingHierarchyIntegrator should initialize the bulk viscosity
        // variable.
        TBOX_ASSERT(d_lambda_init_fcn || d_reset_lambda_fcns.size() > 0);
#endif
        registerVariable(d_lambda_current_idx,
                         d_lambda_new_idx,
                         d_lambda_scratch_idx,
                         d_lambda_var,
                         lambda_cell_ghosts,
                         d_lambda_coarsen_type,
                         d_lambda_refine_type,
                         d_lambda_init_fcn);
    }
    else
    {
        d_lambda_current_idx = invalid_index;
        d_lambda_new_idx = invalid_index;
        d_lambda_scratch_idx = invalid_index;
    }

#if !defined(NDEBUG)
    // AcousticStreamingHierarchyIntegrator should initialize the density variable.
    TBOX_ASSERT(d_rho_init_fcn || d_reset_rho_fcns.size() > 0);
#endif
    registerVariable(d_rho_current_idx,
                     d_rho_new_idx,
                     d_rho_scratch_idx,
                     d_rho_var,
                     side_ghosts,
                     d_rho_coarsen_type,
                     d_rho_refine_type,
                     d_rho_init_fcn);

    // Register plot variables that are maintained by the
    // AcousticStreamingHierarchyIntegrator.
    registerVariable(d_U1_plot_idx, d_U1_plot_var, no_ghosts, getCurrentContext());
    registerVariable(d_U2_plot_idx, d_U2_plot_var, no_ghosts, getCurrentContext());
    registerVariable(d_rho_plot_idx, d_rho_plot_var, no_ghosts, getCurrentContext());
    registerVariable(d_Omega2_idx, d_Omega2_var, no_ghosts, getCurrentContext());

    // Register scratch variables that are maintained by the
    // AcousticStreamingHierarchyIntegrator.
#if (NDIM == 3)
    registerVariable(d_Omega2_Norm_idx, d_Omega2_Norm_var, no_ghosts);
#endif
    registerVariable(d_U2_regrid_idx, d_U2_regrid_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);
    registerVariable(d_U2_src_idx, d_U2_src_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);
    registerVariable(d_indicator2_idx, d_indicator2_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);

    // Register variables for plotting.
    if (d_visit_writer)
    {
        if (d_output_U1)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0)
                {
                    d_visit_writer->registerPlotQuantity("U1_real_x", "SCALAR", d_U1_plot_idx, d);
                    d_visit_writer->registerPlotQuantity("U1_imag_x", "SCALAR", d_U1_plot_idx, NDIM + d);
                }
                if (d == 1)
                {
                    d_visit_writer->registerPlotQuantity("U1_real_y", "SCALAR", d_U1_plot_idx, d);
                    d_visit_writer->registerPlotQuantity("U1_imag_y", "SCALAR", d_U1_plot_idx, NDIM + d);
                }
                if (d == 2)
                {
                    d_visit_writer->registerPlotQuantity("U1_real_z", "SCALAR", d_U1_plot_idx, d);
                    d_visit_writer->registerPlotQuantity("U1_imag_z", "SCALAR", d_U1_plot_idx, NDIM + d);
                }
            }
        }

        if (d_output_rho && d_rho_var)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0)
                {
                    d_visit_writer->registerPlotQuantity("rho_acoustic_x", "SCALAR", d_rho_plot_idx, d);
                }
                if (d == 1)
                {
                    d_visit_writer->registerPlotQuantity("rho_acoustic_y", "SCALAR", d_rho_plot_idx, d);
                }
                if (d == 2)
                {
                    d_visit_writer->registerPlotQuantity("rho_acoustic_z", "SCALAR", d_rho_plot_idx, d);
                }
            }
        }

        if (d_output_P1)
        {
            d_visit_writer->registerPlotQuantity("P1_real", "SCALAR", d_P1_current_idx, REAL);
            d_visit_writer->registerPlotQuantity("P1_imag", "SCALAR", d_P1_current_idx, IMAG);
        }

        if (d_output_U2)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0) d_visit_writer->registerPlotQuantity("U2_x", "SCALAR", d_U2_plot_idx, d);
                if (d == 1) d_visit_writer->registerPlotQuantity("U2_y", "SCALAR", d_U2_plot_idx, d);
                if (d == 2) d_visit_writer->registerPlotQuantity("U2_z", "SCALAR", d_U2_plot_idx, d);
            }
        }

        if (d_output_P2)
        {
            d_visit_writer->registerPlotQuantity("P2", "SCALAR", d_P2_current_idx, 0);
        }

        if (d_output_mu && d_mu_var)
        {
            d_visit_writer->registerPlotQuantity("mu_acoustic", "SCALAR", d_mu_current_idx, 0);
        }

        if (d_output_lambda && d_lambda_var)
        {
            d_visit_writer->registerPlotQuantity("lambda_acoustic", "SCALAR", d_lambda_current_idx, 0);
        }

        if (d_output_Omega2)
        {
#if (NDIM == 2)
            d_visit_writer->registerPlotQuantity("Omega2", "SCALAR", d_Omega2_idx, 0);
#endif
#if (NDIM == 3)
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0) d_visit_writer->registerPlotQuantity("Omega2_x", "SCALAR", d_Omega2_idx, d);
                if (d == 1) d_visit_writer->registerPlotQuantity("Omega2_y", "SCALAR", d_Omega2_idx, d);
                if (d == 2) d_visit_writer->registerPlotQuantity("Omega2_z", "SCALAR", d_Omega2_idx, d);
            }
#endif
        }
    }

    // Initialize and register variables used to compute second-order solver coefficients
    d_pressure_D_var = new SideVariable<NDIM, double>(d_object_name + "::pressure_D");
    d_pressure_D_idx = var_db->registerVariableAndContext(d_pressure_D_var, getCurrentContext(), side_ghosts);

    d_velocity_C_var = new SideVariable<NDIM, double>(d_object_name + "::velocity_C");
    d_velocity_C_idx = var_db->registerVariableAndContext(d_velocity_C_var, getCurrentContext(), no_ghosts);

    d_velocity_L_var = new SideVariable<NDIM, double>(d_object_name + "::velocity_L");
    d_velocity_L_idx = var_db->registerVariableAndContext(d_velocity_L_var, getCurrentContext(), side_ghosts);

    d_projection_D_var = new SideVariable<NDIM, double>(d_object_name + "::projection_D");
    d_projection_D_idx = var_db->registerVariableAndContext(d_projection_D_var, getCurrentContext(), side_ghosts);

#if (NDIM == 2)
    d_velocity_D_var = new NodeVariable<NDIM, double>(d_object_name + "::velocity_D");
#elif (NDIM == 3)
    d_velocity_D_var = new EdgeVariable<NDIM, double>(d_object_name + "::velocity_D");
#endif
    d_velocity_D_idx = var_db->registerVariableAndContext(
        d_velocity_D_var, getCurrentContext(), (NDIM == 2 ? node_ghosts : edge_ghosts));

    d_velocity_D_cc_var = new CellVariable<NDIM, double>(d_object_name + "::velocity_D_cc");
    d_velocity_D_cc_idx = var_db->registerVariableAndContext(d_velocity_D_cc_var, getCurrentContext(), no_ghosts);

#if (NDIM == 2)
    d_mu_interp_var = new NodeVariable<NDIM, double>(d_object_name + "::mu_interp");
#elif (NDIM == 3)
    d_mu_interp_var = new EdgeVariable<NDIM, double>(d_object_name + "::mu_interp");
#endif
    d_mu_interp_idx =
        var_db->registerVariableAndContext(d_mu_interp_var, getCurrentContext(), NDIM == 2 ? node_ghosts : edge_ghosts);

    // Register persistent variables to be used for boundary conditions and other
    // applications.
    // Note: these will not be deallocated.
    Pointer<CellVariable<NDIM, double> > mu_cc_linear_op_var =
        new CellVariable<NDIM, double>(d_object_name + "_mu_cc_linear_op_var",
                                       /*depth*/ 1);
    d_mu_linear_op_idx = var_db->registerVariableAndContext(
        mu_cc_linear_op_var, var_db->getContext(d_object_name + "::mu_linear_op"), mu_cell_ghosts);
    d_mu_interp_linear_op_idx =
        var_db->registerVariableAndContext(d_mu_interp_var,
                                           var_db->getContext(d_object_name + "::mu_interp_linear_op"),
                                           NDIM == 2 ? node_ghosts : edge_ghosts);
    Pointer<SideVariable<NDIM, double> > rho_sc_linear_op_var =
        new SideVariable<NDIM, double>(d_object_name + "_rho_sc_linear_op_var",
                                       /*depth*/ 1);
    d_rho_linear_op_idx = var_db->registerVariableAndContext(
        rho_sc_linear_op_var, var_db->getContext(d_object_name + "::rho_linear_op_var"), side_ghosts);

    // Setup a specialized coarsen algorithm.
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
    Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(d_U2_var, d_U_coarsen_type);
    coarsen_alg->registerCoarsen(d_U2_scratch_idx, d_U2_scratch_idx, coarsen_op);
    registerCoarsenAlgorithm(d_object_name + "::U2_COARSEN_OP", coarsen_alg);

    // Setup the first order solver
    {
        d_first_order_solver->setAcousticAngularFrequency(d_acoustic_freq);
        d_first_order_solver->setSoundSpeed(d_sound_speed);
        d_first_order_solver->setBoundaryConditionCoefficients(d_U1_bc_coefs);
        d_first_order_solver->setHomogeneousBc(false);
    }

    // Setup the Stokes solver for the second order system.
    {
        Pointer<KrylovLinearSolver> p_stokes_krylov_solver = d_stokes_solver;
        Pointer<VCStaggeredStokesOperator> p_vc_stokes_op = p_stokes_krylov_solver->getOperator();
        p_vc_stokes_op->setDPatchDataInterpolationType(d_mu_vc_interp_type);

        Pointer<StaggeredStokesBlockPreconditioner> p_stokes_block_pc = p_stokes_krylov_solver->getPreconditioner();
        p_stokes_block_pc->setVelocitySubdomainSolver(d_velocity_solver);
        p_stokes_block_pc->setPressureSubdomainSolver(d_pressure_solver);

        // Set the velocity subdomain solver interpolation type
        auto p_velocity_solver = dynamic_cast<IBTK::PETScKrylovLinearSolver*>(d_velocity_solver.getPointer());
        Pointer<PoissonFACPreconditioner> p_poisson_fac_pc = p_velocity_solver->getPreconditioner();
        if (p_poisson_fac_pc)
        {
            Pointer<VCSCViscousOpPointRelaxationFACOperator> p_vc_point_fac_op =
                p_poisson_fac_pc->getFACPreconditionerStrategy();
            p_vc_point_fac_op->setDPatchDataInterpolationType(d_mu_vc_interp_type);
            Pointer<VCSCViscousPETScLevelSolver> p_vc_level_solver = p_vc_point_fac_op->getCoarseSolver();
            p_vc_level_solver->setViscosityInterpolationType(d_mu_vc_interp_type);
        }
        Pointer<VCSCViscousPETScLevelSolver> p_vc_level_solver = p_velocity_solver->getPreconditioner();
        if (p_vc_level_solver)
        {
            p_vc_level_solver->setViscosityInterpolationType(d_mu_vc_interp_type);
        }

        Pointer<VCSCViscousOperator> p_velocity_op = p_velocity_solver->getOperator();
        p_velocity_op->setDPatchDataInterpolationType(d_mu_vc_interp_type);
    }

    // Setup a boundary op to set velocity boundary conditions on regrid.
    //(NOT SURE WHAT TO DO HERE?)
    // d_fill_after_regrid_phys_bdry_bc_op.reset(new CartSideRobinPhysBdryOp(d_U_scratch_idx, d_U_bc_coefs, false));

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
AcousticStreamingHierarchyIntegrator::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                               Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    HierarchyIntegrator::initializePatchHierarchy(hierarchy, gridding_alg);
    return;
} // initializePatchHierarchy

void
AcousticStreamingHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                   const double new_time,
                                                                   const int num_cycles)
{
    HierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
        level->allocatePatchData(d_velocity_C_idx, current_time);
        level->allocatePatchData(d_velocity_L_idx, current_time);
        level->allocatePatchData(d_velocity_D_idx, current_time);
        level->allocatePatchData(d_velocity_D_cc_idx, current_time);
        level->allocatePatchData(d_pressure_D_idx, current_time);
        level->allocatePatchData(d_projection_D_idx, current_time);
        level->allocatePatchData(d_mu_interp_idx, current_time);
        if (!level->checkAllocated(d_mu_linear_op_idx)) level->allocatePatchData(d_mu_linear_op_idx, current_time);
        if (!level->checkAllocated(d_mu_interp_linear_op_idx))
            level->allocatePatchData(d_mu_interp_linear_op_idx, current_time);
        if (!level->checkAllocated(d_rho_linear_op_idx)) level->allocatePatchData(d_rho_linear_op_idx, current_time);
    }

    // Preprocess the operators and solvers
    preprocessOperatorsAndSolvers(current_time, new_time);

    // Preprocess Brinkman penalization objects.
    for (auto& brinkman_force : d_brinkman_force)
    {
        brinkman_force->setTimeInterval(current_time, new_time);
        brinkman_force->preprocessComputeBrinkmanPenalization(current_time, new_time, num_cycles);
    }

    // Allocate memory for RHS vectors of the two solvers.
    d_rhs1_vec->allocateVectorData(current_time);
    d_rhs1_vec->setToScalar(0.0);
    d_U2_rhs_vec->allocateVectorData(current_time);
    d_U2_rhs_vec->setToScalar(0.0);
    d_P2_rhs_vec->allocateVectorData(current_time);
    d_P2_rhs_vec->setToScalar(0.0);

    // Copy the current solution to new variables (to be used as guess solution)
    d_hier_sc_data_ops->copyData(d_U1_new_idx, d_U1_current_idx);
    d_hier_sc_data_ops->copyData(d_U2_new_idx, d_U2_current_idx);
    d_hier_cc_data_ops->copyData(d_P1_new_idx, d_P1_current_idx);
    d_hier_cc_data_ops->copyData(d_P2_new_idx, d_P2_current_idx);

    // Cache BC data.
    d_bc_helper->cacheBcCoefData(d_so_bc_coefs, new_time, d_hierarchy);

    // Set up inhomogeneous BCs.
    d_first_order_solver->setHomogeneousBc(false);
    d_stokes_solver->setHomogeneousBc(false);

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);

    return;
} // preprocessIntegrateHierarchy

void
AcousticStreamingHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                         const double new_time,
                                                         const int cycle_num)
{
    HierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);

    // Get the coarsest and finest level numbers.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Update density
    const double apply_time = new_time;
    for (unsigned k = 0; k < d_reset_rho_fcns.size(); ++k)
    {
        d_reset_rho_fcns[k](d_rho_new_idx,
                            d_rho_var,
                            d_hier_math_ops,
                            cycle_num,
                            apply_time,
                            current_time,
                            new_time,
                            d_reset_rho_fcns_ctx[k]);
    }
    d_hier_sc_data_ops->copyData(d_rho_scratch_idx,
                                 d_rho_new_idx,
                                 /*interior_only*/ true);
    d_rho_bdry_bc_fill_op->fillData(new_time);

    // Update shear viscosity.
    for (unsigned k = 0; k < d_reset_mu_fcns.size(); ++k)
    {
        d_reset_mu_fcns[k](d_mu_new_idx,
                           d_mu_var,
                           d_hier_math_ops,
                           cycle_num,
                           apply_time,
                           current_time,
                           new_time,
                           d_reset_mu_fcns_ctx[k]);
    }
    d_hier_cc_data_ops->copyData(d_mu_scratch_idx,
                                 d_mu_new_idx,
                                 /*interior_only*/ true);
    d_mu_bdry_bc_fill_op->fillData(new_time);

    // Interpolate onto node or edge centers
    if (d_mu_vc_interp_type == VC_AVERAGE_INTERP)
    {
        d_hier_math_ops->interp_ghosted(
            d_mu_interp_idx, d_mu_interp_var, d_mu_scratch_idx, d_mu_var, d_no_fill_op, new_time);
    }
    else if (d_mu_vc_interp_type == VC_HARMONIC_INTERP)
    {
        d_hier_math_ops->harmonic_interp_ghosted(
            d_mu_interp_idx, d_mu_interp_var, d_mu_scratch_idx, d_mu_var, d_no_fill_op, new_time);
    }
    else
    {
        TBOX_ERROR("this statement should not be reached");
    }

    // Store viscosity for later use
    d_hier_cc_data_ops->copyData(d_mu_linear_op_idx,
                                 d_mu_scratch_idx,
                                 /*interior_only*/ false);
#if (NDIM == 2)
    d_hier_nc_data_ops->copyData(d_mu_interp_linear_op_idx,
                                 d_mu_interp_idx,
                                 /*interior_only*/ false);
#elif (NDIM == 3)
    d_hier_ec_data_ops->copyData(d_mu_interp_linear_op_idx,
                                 d_mu_interp_idx,
                                 /*interior_only*/ false);
#endif

    // Update bulk viscosity.
    if (d_lambda_var)
    {
        for (unsigned k = 0; k < d_reset_lambda_fcns.size(); ++k)
        {
            d_reset_lambda_fcns[k](d_lambda_new_idx,
                                   d_lambda_var,
                                   d_hier_math_ops,
                                   cycle_num,
                                   apply_time,
                                   current_time,
                                   new_time,
                                   d_reset_lambda_fcns_ctx[k]);
        }
        d_hier_cc_data_ops->copyData(d_lambda_scratch_idx,
                                     d_lambda_new_idx,
                                     /*interior_only*/ true);
        d_lambda_bdry_bc_fill_op->fillData(new_time);
    }

    // Synchronize the newest density
    using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
    SynchronizationTransactionComponent rho_scratch_synch_transaction =
        SynchronizationTransactionComponent(d_rho_scratch_idx, d_rho_coarsen_type);
    d_side_synch2_op->resetTransactionComponent(rho_scratch_synch_transaction);
    d_side_synch2_op->synchronizeData(d_integrator_time);
    SynchronizationTransactionComponent default_synch_transaction =
        SynchronizationTransactionComponent(d_U2_scratch_idx, d_U_coarsen_type);
    d_side_synch2_op->resetTransactionComponent(default_synch_transaction);

    // Store density for later use (e.g., apps defining gravity body force)
    d_hier_sc_data_ops->copyData(d_rho_linear_op_idx,
                                 d_rho_scratch_idx,
                                 /*interior_only*/ true);

    // Update the solvers and operators to take into account new state variables
    updateOperatorsAndSolvers(current_time, new_time, cycle_num);

    // Setup the solution and right-hand-side vectors.
    setupSolverVectors(d_sol1_vec, d_rhs1_vec, d_sol2_vec, d_rhs2_vec, current_time, new_time, cycle_num);

    // Solve for u1(n+1), p1(n+1)
    d_first_order_solver->solveSystem(*d_sol1_vec, *d_rhs1_vec);

    // Solve for u2(n+1), p2(n+1).
    d_stokes_solver->solveSystem(*d_sol2_vec, *d_rhs2_vec);

    if (d_enable_logging && d_enable_logging_solver_iterations)
        plog << d_object_name
             << "::integrateHierarchy(): stokes solve number of iterations = " << d_stokes_solver->getNumIterations()
             << "\n";
    if (d_enable_logging)
        plog << d_object_name
             << "::integrateHierarchy(): stokes solve residual norm        = " << d_stokes_solver->getResidualNorm()
             << "\n";
    if (d_explicitly_remove_so_nullspace) removeSecondOrderNullSpace(d_sol2_vec);

    // Reset the solution and right-hand-side vectors.
    resetSolverVectors(d_sol1_vec, d_rhs1_vec, d_sol2_vec, d_rhs2_vec, current_time, new_time, cycle_num);

    // Re-update viscosity if it is maintained by the integrator
    // using the newest available data from INS and advection-diffusion solvers
    if (d_current_num_cycles == cycle_num + 1)
    {
        for (unsigned k = 0; k < d_reset_mu_fcns.size(); ++k)
        {
            d_reset_mu_fcns[k](d_mu_new_idx,
                               d_mu_var,
                               d_hier_math_ops,
                               cycle_num,
                               apply_time,
                               current_time,
                               new_time,
                               d_reset_mu_fcns_ctx[k]);
        }
    }

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy

void
AcousticStreamingHierarchyIntegrator::postprocessIntegrateHierarchy(double current_time,
                                                                    double new_time,
                                                                    bool skip_synchronize_new_state_data,
                                                                    int num_cycles)
{
    HierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_enable_logging)
            plog << d_object_name << "::postprocessIntegrateHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(NEW_DATA);
    }

    // Determine the CFL number.
    if (!d_parent_integrator)
    {
        double cfl_max = 0.0;
        PatchSideDataOpsReal<NDIM, double> patch_sc_ops;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double dx_min = *(std::min_element(dx, dx + NDIM));
                Pointer<SideData<NDIM, double> > u2_sc_new_data = patch->getPatchData(d_U2_new_idx);
                double u2_max = 0.0;
                u2_max = patch_sc_ops.maxNorm(u2_sc_new_data, patch_box);
                cfl_max = std::max(cfl_max, u2_max * dt / dx_min);
            }
        }
        cfl_max = IBTK_MPI::maxReduction(cfl_max);
        if (d_enable_logging)
            plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl_max << "\n";
    }

    // Compute max |Omega|_2.
    if (d_using_vorticity_tagging)
    {
        d_hier_sc_data_ops->copyData(d_U2_scratch_idx, d_U2_new_idx);
        d_hier_math_ops->curl(d_Omega2_idx, d_Omega2_var, d_U2_scratch_idx, d_U2_var, d_U2_bdry_bc_fill_op, new_time);

#if (NDIM == 3)
        d_hier_math_ops->pointwiseL2Norm(d_Omega2_Norm_idx, d_Omega2_Norm_var, d_Omega2_idx, d_Omega2_var);
#endif
        const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
#if (NDIM == 2)
        d_Omega2_max = d_hier_cc_data_ops->maxNorm(d_Omega2_idx, wgt_cc_idx);
#endif
#if (NDIM == 3)
        d_Omega2_max = d_hier_cc_data_ops->max(d_Omega2_Norm_idx, wgt_cc_idx);
#endif
    }

    // Deallocate scratch data.
    deallocate_vector_data(*d_U2_rhs_vec);
    deallocate_vector_data(*d_P2_rhs_vec);

    // Deallocate any temporary data used to compute coefficients
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_velocity_C_idx);
        level->deallocatePatchData(d_velocity_L_idx);
        level->deallocatePatchData(d_velocity_D_idx);
        level->deallocatePatchData(d_velocity_D_cc_idx);
        level->deallocatePatchData(d_pressure_D_idx);
        level->deallocatePatchData(d_projection_D_idx);
        level->deallocatePatchData(d_mu_interp_idx);
    }

    // Postprocess Brinkman penalization objects.
    for (auto& brinkman_force : d_brinkman_force)
    {
        brinkman_force->postprocessComputeBrinkmanPenalization(current_time, new_time, num_cycles);
    }

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;

} // postprocessIntegrateHierarchy

void
AcousticStreamingHierarchyIntegrator::removeSecondOrderNullSpace(
    const Pointer<SAMRAIVectorReal<NDIM, double> >& sol2_vec)
{
    if (d_nul2_vecs.empty()) return;
    for (const auto& nul2_vec : d_nul2_vecs)
    {
        const double sol2_dot_nul2 = sol2_vec->dot(nul2_vec);
        const double nul2_L2_norm_square = nul2_vec->dot(nul2_vec);
        sol2_vec->axpy(-sol2_dot_nul2 / nul2_L2_norm_square, nul2_vec, sol2_vec);
    }
    return;
} // removeSecondOrderNullSpace

/////////////////////////////// PROTECTED //////////////////////////////////////
void
AcousticStreamingHierarchyIntegrator::regridHierarchyBeginSpecialized()
{
    // intentionally left blank
    return;
} // regridHierarchyBeginSpecialized

void
AcousticStreamingHierarchyIntegrator::regridHierarchyEndSpecialized()
{
    // intentionally left blank
    return;
} // regridHierarchyEndSpecialized

void
AcousticStreamingHierarchyIntegrator::initializeCompositeHierarchyDataSpecialized(double /*init_data_time*/,
                                                                                  bool /*initial_time*/)
{
    // intentionally left blank
    return;
} // initializeCompositeHierarchyDataSpecialized

void
AcousticStreamingHierarchyIntegrator::initializeLevelDataSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool /*can_be_refined*/,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > base_old_level,
    const bool /*allocate_data*/)
{
    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const Pointer<PatchLevel<NDIM> > old_level = base_old_level;
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (old_level)
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Correct the divergence of the interpolated velocity data.
    if (!initial_time && level_number > 0)
    {
        // Allocate scratch data.
        ComponentSelector scratch_data;
        scratch_data.setFlag(d_U2_regrid_idx);
        scratch_data.setFlag(d_U2_src_idx);
        scratch_data.setFlag(d_indicator2_idx);
        level->allocatePatchData(scratch_data, init_data_time);
        if (old_level) old_level->allocatePatchData(scratch_data, init_data_time);

        // Set the indicator data to equal "0" in each patch of the new patch
        // level, and initialize values of U to cause floating point errors if
        // we fail to re-initialize it properly.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > indicator2_data = patch->getPatchData(d_indicator2_idx);
            indicator2_data->fillAll(0.0);

            Pointer<SideData<NDIM, double> > U2_current_data = patch->getPatchData(d_U2_current_idx);
            Pointer<SideData<NDIM, double> > U2_regrid_data = patch->getPatchData(d_U2_regrid_idx);
            Pointer<SideData<NDIM, double> > U2_src_data = patch->getPatchData(d_U2_src_idx);
            U2_current_data->fillAll(std::numeric_limits<double>::quiet_NaN());
            U2_regrid_data->fillAll(std::numeric_limits<double>::quiet_NaN());
            U2_src_data->fillAll(std::numeric_limits<double>::quiet_NaN());
        }

        if (old_level)
        {
            // Set the indicator data to equal "1" on each patch of the old
            // patch level and reset U.
            for (PatchLevel<NDIM>::Iterator p(old_level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = old_level->getPatch(p());

                Pointer<SideData<NDIM, double> > indicator2_data = patch->getPatchData(d_indicator2_idx);
                indicator2_data->fillAll(1.0);

                Pointer<SideData<NDIM, double> > U2_current_data = patch->getPatchData(d_U2_current_idx);
                Pointer<SideData<NDIM, double> > U2_regrid_data = patch->getPatchData(d_U2_regrid_idx);
                Pointer<SideData<NDIM, double> > U2_src_data = patch->getPatchData(d_U2_src_idx);
                U2_regrid_data->copy(*U2_current_data);
                U2_src_data->copy(*U2_current_data);
            }

            // Create a communications schedule to copy data from the old patch
            // level to the new patch level.
            //
            // Note that this will set the indicator data to equal "1" at each
            // location in the new patch level that is a copy of a location from
            // the old patch level.
            RefineAlgorithm<NDIM> copy_data;
            copy_data.registerRefine(d_U2_regrid_idx, d_U2_regrid_idx, d_U2_regrid_idx, nullptr);
            copy_data.registerRefine(d_U2_src_idx, d_U2_src_idx, d_U2_src_idx, nullptr);
            copy_data.registerRefine(d_indicator2_idx, d_indicator2_idx, d_indicator2_idx, nullptr);
            ComponentSelector bc_fill_data;
            bc_fill_data.setFlag(d_U2_regrid_idx);
            bc_fill_data.setFlag(d_U2_src_idx);
            CartSideRobinPhysBdryOp phys_bdry_bc_op(bc_fill_data, d_U2_bc_coefs, false);
            copy_data.createSchedule(level, old_level, &phys_bdry_bc_op)->fillData(init_data_time);
        }

        // Setup the divergence- and curl-preserving prolongation refine
        // algorithm and refine the velocity data.
        RefineAlgorithm<NDIM> fill_div_free_prolongation;
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        fill_div_free_prolongation.registerRefine(d_U2_current_idx, d_U2_current_idx, d_U2_regrid_idx, nullptr);
        Pointer<RefineOperator<NDIM> > refine_op = grid_geom->lookupRefineOperator(d_U2_var, d_U_refine_type);
        Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(d_U2_var, d_U_coarsen_type);
        CartSideRobinPhysBdryOp phys_bdry_bc_op(d_U2_regrid_idx, d_U2_bc_coefs, false);
        CartSideDoubleDivPreservingRefine div_preserving_op(
            d_U2_regrid_idx, d_U2_src_idx, d_indicator2_idx, refine_op, coarsen_op, init_data_time, &phys_bdry_bc_op);
        fill_div_free_prolongation.createSchedule(level, old_level, level_number - 1, hierarchy, &div_preserving_op)
            ->fillData(init_data_time);

        // Free scratch data.
        level->deallocatePatchData(scratch_data);
        if (old_level) old_level->deallocatePatchData(scratch_data);
    }

    // Initialize level data at the initial time.
    if (initial_time)
    {
        // Initialize the maximum value of |Omega|_2 on the grid.
        if (d_using_vorticity_tagging)
        {
            if (level_number == 0) d_Omega2_max = 0.0;

            // Allocate scratch data.
            for (int ln = 0; ln <= level_number; ++ln)
            {
                hierarchy->getPatchLevel(ln)->allocatePatchData(d_U2_scratch_idx, init_data_time);
#if (NDIM == 3)
                hierarchy->getPatchLevel(ln)->allocatePatchData(d_Omega2_Norm_idx, init_data_time);
#endif
            }

            // Fill ghost cells.
            HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
            Pointer<HierarchyCellDataOpsReal<NDIM, double> > hier_cc_data_ops =
                hier_ops_manager->getOperationsDouble(d_U2_plot_var, d_hierarchy, true);
            Pointer<HierarchySideDataOpsReal<NDIM, double> > hier_sc_data_ops =
                hier_ops_manager->getOperationsDouble(d_U2_var, d_hierarchy, true);
            hier_sc_data_ops->resetLevels(0, level_number);
            hier_sc_data_ops->copyData(d_U2_scratch_idx, d_U2_current_idx);
            using InterpolationTransactionComponent =
                HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            InterpolationTransactionComponent U2_bc_component(d_U2_scratch_idx,
                                                              DATA_REFINE_TYPE,
                                                              USE_CF_INTERPOLATION,
                                                              DATA_COARSEN_TYPE,
                                                              d_bdry_extrap_type, // TODO: update variable name
                                                              CONSISTENT_TYPE_2_BDRY,
                                                              d_U2_bc_coefs);
            HierarchyGhostCellInterpolation U2_bdry_bc_fill_op;
            U2_bdry_bc_fill_op.initializeOperatorState(U2_bc_component, d_hierarchy, 0, level_number);
            U2_bdry_bc_fill_op.fillData(init_data_time);

            // Compute max |Omega|_2.
            HierarchyMathOps hier_math_ops(d_object_name + "::HierarchyLevelMathOps", d_hierarchy, 0, level_number);
            hier_math_ops.curl(
                d_Omega2_idx, d_Omega2_var, d_U2_scratch_idx, d_U2_var, d_U2_bdry_bc_fill_op, init_data_time);
#if (NDIM == 3)
            hier_math_ops.pointwiseL2Norm(d_Omega2_Norm_idx, d_Omega2_Norm_var, d_Omega2_idx, d_Omega2_var);
#endif
            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
#if (NDIM == 2)
            d_Omega2_max = hier_cc_data_ops->maxNorm(d_Omega2_idx, wgt_cc_idx);
#endif
#if (NDIM == 3)
            d_Omega2_max = hier_cc_data_ops->max(d_Omega2_Norm_idx, wgt_cc_idx);
#endif
            // Deallocate scratch data.
            for (int ln = 0; ln <= level_number; ++ln)
            {
                hierarchy->getPatchLevel(ln)->deallocatePatchData(d_U2_scratch_idx);
#if (NDIM == 3)
                hierarchy->getPatchLevel(ln)->deallocatePatchData(d_Omega2_Norm_idx);
#endif
            }
        }
    }
    return;
} // initializeLevelDataSpecialized

void
AcousticStreamingHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((coarsest_level >= 0) && (coarsest_level <= finest_level) &&
                (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(hierarchy->getPatchLevel(ln));
    }
#else
    NULL_USE(coarsest_level);
    NULL_USE(finest_level);
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset the hierarchy operations objects for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_nc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_nc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_ec_data_ops->setPatchHierarchy(hierarchy);
    d_hier_ec_data_ops->resetLevels(0, finest_hier_level);

    // Setup the patch boundary filling objects.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent U2_bc_component(d_U2_scratch_idx,
                                                      DATA_REFINE_TYPE,
                                                      USE_CF_INTERPOLATION,
                                                      DATA_COARSEN_TYPE,
                                                      d_bdry_extrap_type, // TODO: update variable name
                                                      CONSISTENT_TYPE_2_BDRY,
                                                      d_U2_bc_coefs);
    d_U2_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_U2_bdry_bc_fill_op->initializeOperatorState(U2_bc_component, d_hierarchy);

    InterpolationTransactionComponent P2_bc_component(d_P2_scratch_idx,
                                                      DATA_REFINE_TYPE,
                                                      USE_CF_INTERPOLATION,
                                                      DATA_COARSEN_TYPE,
                                                      d_bdry_extrap_type, // TODO: update variable name
                                                      CONSISTENT_TYPE_2_BDRY,
                                                      d_P2_bc_coef);
    d_P2_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_P2_bdry_bc_fill_op->initializeOperatorState(P2_bc_component, d_hierarchy);

    InterpolationTransactionComponent mu_bc_component(
        d_mu_scratch_idx, d_mu_refine_type, false, d_mu_coarsen_type, d_mu_bdry_extrap_type, false, d_mu_bc_coef);
    d_mu_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_mu_bdry_bc_fill_op->initializeOperatorState(mu_bc_component, d_hierarchy);

    if (d_lambda_var)
    {
        InterpolationTransactionComponent lambda_bc_component(d_lambda_scratch_idx,
                                                              d_lambda_refine_type,
                                                              false,
                                                              d_lambda_coarsen_type,
                                                              d_lambda_bdry_extrap_type,
                                                              false,
                                                              d_lambda_bc_coef);
        d_lambda_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
        d_lambda_bdry_bc_fill_op->initializeOperatorState(lambda_bc_component, d_hierarchy);
    }

    InterpolationTransactionComponent rho_bc_component(
        d_rho_scratch_idx, d_rho_refine_type, false, d_rho_coarsen_type, d_rho_bdry_extrap_type, false, d_rho_bc_coef);
    d_rho_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_rho_bdry_bc_fill_op->initializeOperatorState(rho_bc_component, d_hierarchy);

    // Setup the patch boundary synchronization objects.
    using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;

    SynchronizationTransactionComponent synch1_transaction =
        SynchronizationTransactionComponent(d_U1_scratch_idx, d_U_coarsen_type);
    d_side_synch1_op = new SideDataSynchronization();
    d_side_synch1_op->initializeOperatorState(synch1_transaction, d_hierarchy);

    SynchronizationTransactionComponent synch2_transaction =
        SynchronizationTransactionComponent(d_U2_scratch_idx, d_U_coarsen_type);
    d_side_synch2_op = new SideDataSynchronization();
    d_side_synch2_op->initializeOperatorState(synch2_transaction, d_hierarchy);

    // Indicate that vectors and solvers need to be re-initialized.
    d_coarsest_reset_ln = coarsest_level;
    d_finest_reset_ln = finest_level;
    d_vectors_need_init = true;
    d_velocity_solver_needs_init = true;
    d_pressure_solver_needs_init = true;
    d_stokes_solver_needs_init = true;
    d_first_order_solver_needs_init = true;
    return;
} // resetHierarchyConfigurationSpecialized

void
AcousticStreamingHierarchyIntegrator::applyGradientDetectorSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double /*error_data_time*/,
    const int tag_index,
    const bool /*initial_time*/,
    const bool /*uses_richardson_extrapolation_too*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Tag cells based on the magnitude of the vorticity.
    //
    // Note that if either the relative or absolute threshold is zero for a
    // particular level, no tagging is performed on that level.
    if (d_using_vorticity_tagging)
    {
        double Omega2_rel_thresh = 0.0;
        if (d_Omega2_rel_thresh.size() > 0)
        {
            Omega2_rel_thresh =
                d_Omega2_rel_thresh[std::max(std::min(level_number, d_Omega2_rel_thresh.size() - 1), 0)];
        }
        double Omega2_abs_thresh = 0.0;
        if (d_Omega2_abs_thresh.size() > 0)
        {
            Omega2_abs_thresh =
                d_Omega2_abs_thresh[std::max(std::min(level_number, d_Omega2_abs_thresh.size() - 1), 0)];
        }
        if (Omega2_rel_thresh > 0.0 || Omega2_abs_thresh > 0.0)
        {
            double thresh = std::numeric_limits<double>::max();
            if (Omega2_rel_thresh > 0.0) thresh = std::min(thresh, Omega2_rel_thresh * d_Omega2_max);
            if (Omega2_abs_thresh > 0.0) thresh = std::min(thresh, Omega2_abs_thresh);
            thresh += std::sqrt(std::numeric_limits<double>::epsilon());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, int> > tags_data = patch->getPatchData(tag_index);
                Pointer<CellData<NDIM, double> > Omega2_data = patch->getPatchData(d_Omega2_idx);
                for (CellIterator<NDIM> ic(patch_box); ic; ic++)
                {
                    const hier::Index<NDIM>& i = ic();
#if (NDIM == 2)
                    if (std::abs((*Omega2_data)(i)) > thresh)
                    {
                        (*tags_data)(i) = 1;
                    }
#endif
#if (NDIM == 3)
                    double norm_Omega2_sq = 0.0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        norm_Omega2_sq += (*Omega2_data)(i, d) * (*Omega2_data)(i, d);
                    }
                    const double norm_Omega2 = std::sqrt(norm_Omega2_sq);
                    if (norm_Omega2 > thresh)
                    {
                        (*tags_data)(i) = 1;
                    }
#endif
                }
            }
        }
    }
    return;
} // applyGradientDetectorSpecialized

void
AcousticStreamingHierarchyIntegrator::setupPlotDataSpecialized()
{
    Pointer<VariableContext> ctx = getCurrentContext();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    static const bool synch_cf_interface = true;

    // Interpolate u to cell centers.
    if (d_output_U1)
    {
        const int U1_sc_idx = var_db->mapVariableAndContextToIndex(d_U1_var, ctx);
        const int U1_cc_plot_idx = var_db->mapVariableAndContextToIndex(d_U1_plot_var, ctx);
        d_hier_math_ops->interp(
            U1_cc_plot_idx, d_U1_plot_var, U1_sc_idx, d_U1_var, d_no_fill_op, d_integrator_time, synch_cf_interface);
    }
    if (d_output_U2)
    {
        const int U2_sc_idx = var_db->mapVariableAndContextToIndex(d_U2_var, ctx);
        const int U2_cc_plot_idx = var_db->mapVariableAndContextToIndex(d_U2_plot_var, ctx);
        d_hier_math_ops->interp(
            U2_cc_plot_idx, d_U2_plot_var, U2_sc_idx, d_U2_var, d_no_fill_op, d_integrator_time, synch_cf_interface);
    }
    if (d_output_rho)
    {
        const int rho_sc_idx = var_db->mapVariableAndContextToIndex(d_rho_var, ctx);
        const int rho_cc_plot_idx = var_db->mapVariableAndContextToIndex(d_rho_plot_var, ctx);
        Pointer<SideVariable<NDIM, double> > rho_sc_var = d_rho_var;
#if !defined(NDEBUG)
        TBOX_ASSERT(rho_sc_var);
#endif
        d_hier_math_ops->interp(rho_cc_plot_idx,
                                d_rho_plot_var,
                                rho_sc_idx,
                                rho_sc_var,
                                d_no_fill_op,
                                d_integrator_time,
                                synch_cf_interface);
    }

    // Compute Omega = curl U.
    if (d_output_Omega2)
    {
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(d_U2_scratch_idx, d_integrator_time);
        }
        d_hier_sc_data_ops->copyData(d_U2_scratch_idx, d_U2_current_idx);
        d_U2_bdry_bc_fill_op->fillData(d_integrator_time);
        d_hier_math_ops->curl(d_Omega2_idx, d_Omega2_var, d_U2_scratch_idx, d_U2_var, d_no_fill_op, d_integrator_time);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(d_U2_scratch_idx);
        }
    }

    return;
} // setupPlotDataSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
AcousticStreamingHierarchyIntegrator::getFromInput(Pointer<Database> input_db, const bool is_from_restart)
{
    if (input_db->keyExists("acoustic_frequency"))
        d_acoustic_freq = 2.0 * M_PI * input_db->getDouble("acoustic_frequency");
    if (input_db->keyExists("frequency")) d_acoustic_freq = 2.0 * M_PI * input_db->getDouble("frequency");
    if (input_db->keyExists("acoustic_angular_frequency"))
        d_acoustic_freq = input_db->getDouble("acoustic_angular_frequency");
    if (input_db->keyExists("angular_frequency")) d_acoustic_freq = input_db->getDouble("angular_frequency");
    if (input_db->keyExists("sound_speed")) d_sound_speed = input_db->getDouble("sound_speed");

    // Get the interpolation type for the material properties
    if (input_db->keyExists("vc_interpolation_type"))
    {
        d_rho_vc_interp_type = IBTK::string_to_enum<VCInterpType>(input_db->getString("vc_interpolation_type"));
        d_mu_vc_interp_type = IBTK::string_to_enum<VCInterpType>(input_db->getString("vc_interpolation_type"));
    }
    if (input_db->keyExists("rho_vc_interpolation_type"))
    {
        d_rho_vc_interp_type = IBTK::string_to_enum<VCInterpType>(input_db->getString("rho_vc_interpolation_type"));
    }
    if (input_db->keyExists("mu_vc_interpolation_type"))
    {
        d_mu_vc_interp_type = IBTK::string_to_enum<VCInterpType>(input_db->getString("mu_vc_interpolation_type"));
    }
    switch (d_rho_vc_interp_type)
    {
    case VC_HARMONIC_INTERP:
    case VC_AVERAGE_INTERP:
        break;
    default:
        TBOX_ERROR(d_object_name << "::getFromInput():\n"
                                 << "  unsupported density interpolation type: "
                                 << IBTK::enum_to_string<VCInterpType>(d_rho_vc_interp_type) << " \n"
                                 << "  valid choices are: VC_HARMONIC_INTERP, VC_AVERAGE_INTERP\n");
    }
    switch (d_mu_vc_interp_type)
    {
    case VC_HARMONIC_INTERP:
    case VC_AVERAGE_INTERP:
        break;
    default:
        TBOX_ERROR(d_object_name << "::getFromInput():\n"
                                 << "  unsupported viscosity interpolation type: "
                                 << IBTK::enum_to_string<VCInterpType>(d_mu_vc_interp_type) << " \n"
                                 << "  valid choices are: VC_HARMONIC_INTERP, VC_AVERAGE_INTERP\n");
    }

    if (input_db->keyExists("first_order_solver_db"))
    {
        d_first_order_solver_db = input_db->getDatabase("first_order_solver_db");
    }
    else
    {
        d_first_order_solver_db = new MemoryDatabase("first_order_solver_db");
        d_first_order_solver_db->putString("ksp_type", "preonly");
        d_first_order_solver_db->putString("pc_type", "lu");
        d_first_order_solver_db->putBool("initial_guess_nonzero", false);
    }

    if (input_db->keyExists("stokes_solver_db"))
    {
        d_stokes_solver_db = input_db->getDatabase("stokes_solver_db");
    }
    else
    {
        d_stokes_solver_db = new MemoryDatabase("stokes_solver_db");
        d_stokes_solver_db->putString("ksp_type", "fgmres");
    }

    if (input_db->keyExists("precond_reinit_interval"))
        d_precond_reinit_interval = input_db->getInteger("precond_reinit_interval");
    if (input_db->keyExists("stokes_precond_db"))
    {
        d_stokes_precond_db = input_db->getDatabase("stokes_precond_db");
    }
    else
    {
        d_stokes_precond_db = new MemoryDatabase("stokes_precond_db");
        d_stokes_precond_db->putInteger("max_iterations", 1);
    }

    if (input_db->keyExists("velocity_solver_db"))
    {
        d_velocity_solver_db = input_db->getDatabase("velocity_solver_db");
    }
    else
    {
        d_velocity_solver_db = new MemoryDatabase("velocity_solver_db");
        d_velocity_solver_db->putString("ksp_type", "richardson");
        d_velocity_solver_db->putInteger("max_iterations", 10);
        d_velocity_solver_db->putDouble("rel_residual_tol", 1.0e-1);
    }

    if (input_db->keyExists("velocity_precond_db"))
    {
        d_velocity_precond_db = input_db->getDatabase("velocity_precond_db");
    }
    else
    {
        d_velocity_precond_db = new MemoryDatabase("velocity_precond_db");
        d_velocity_precond_db->putString("coarse_solver_type", DEFAULT_VC_VELOCITY_LEVEL_SOLVER);
        d_velocity_precond_db->putInteger("max_iterations", 1);
    }

    if (input_db->keyExists("pressure_solver_db"))
    {
        d_pressure_solver_db = input_db->getDatabase("pressure_solver_db");
    }
    else
    {
        d_pressure_solver_db = new MemoryDatabase("pressure_solver_db");
        d_pressure_solver_db->putString("ksp_type", "richardson");
        d_pressure_solver_db->putInteger("max_iterations", 10);
        d_pressure_solver_db->putDouble("rel_residual_tol", 1.0e-1);
    }

    if (input_db->keyExists("pressure_precond_db"))
    {
        d_pressure_precond_db = input_db->getDatabase("pressure_precond_db");
    }
    else
    {
        d_pressure_precond_db = new MemoryDatabase("pressure_precond_db");
        d_pressure_precond_db->putInteger("max_iterations", 1);
    }

    if (input_db->keyExists("num_cycles")) d_num_cycles = input_db->getInteger("num_cycles");

    if (input_db->keyExists("cfl"))
        d_cfl_max = input_db->getDouble("cfl");
    else if (input_db->keyExists("cfl_max"))
        d_cfl_max = input_db->getDouble("cfl_max");
    else if (input_db->keyExists("CFL"))
        d_cfl_max = input_db->getDouble("CFL");
    else if (input_db->keyExists("CFL_max"))
        d_cfl_max = input_db->getDouble("CFL_max");

    if (input_db->keyExists("normalize_pressure")) d_normalize_pressure = input_db->getBool("normalize_pressure");
    if (input_db->keyExists("normalize_velocity")) d_normalize_velocity = input_db->getBool("normalize_velocity");
    if (input_db->keyExists("explicitly_remove_nullspace"))
        d_explicitly_remove_so_nullspace = input_db->getBool("explicitly_remove_nullspace");
    else if (input_db->keyExists("explicitly_remove_so_nullspace"))
        d_explicitly_remove_so_nullspace = input_db->getBool("explicitly_remove_so_nullspace");
    else if (input_db->keyExists("explicitly_remove_second_order_nullspace"))
        d_explicitly_remove_so_nullspace = input_db->getBool("explicitly_remove_second_order_nullspace");

    if (input_db->keyExists("output_U1")) d_output_U1 = input_db->getBool("output_U1");
    if (input_db->keyExists("output_U2")) d_output_U2 = input_db->getBool("output_U2");
    if (input_db->keyExists("output_P1")) d_output_P1 = input_db->getBool("output_P1");
    if (input_db->keyExists("output_P2")) d_output_P2 = input_db->getBool("output_P2");
    if (input_db->keyExists("output_Omega2")) d_output_Omega2 = input_db->getBool("output_Omega2");
    if (input_db->keyExists("output_rho")) d_output_rho = input_db->getBool("output_rho");
    if (input_db->keyExists("output_mu")) d_output_mu = input_db->getBool("output_mu");
    if (input_db->keyExists("output_lambda")) d_output_lambda = input_db->getBool("output_lambda");

    return;
} // getFromInput

void
AcousticStreamingHierarchyIntegrator::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("ACOUSTIC_STREAMING_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != ACOUSTIC_STREAMING_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }

    d_acoustic_freq = db->getDouble("acoustic_angular_frequency");
    d_sound_speed = db->getDouble("sound_speed");

    d_num_cycles = db->getInteger("num_cycles");
    d_cfl_max = db->getDouble("cfl_max");

    d_normalize_pressure = db->getBool("normalize_pressure");
    d_normalize_velocity = db->getBool("normalize_velocity");
    return;
} // getFromRestart

void
AcousticStreamingHierarchyIntegrator::preprocessOperatorsAndSolvers(const double current_time, const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();

    // Setup solver vectors for the second order system.
    const bool has_velocity_nullspace = d_normalize_velocity;
    const bool has_pressure_nullspace = d_normalize_pressure;
    if (d_vectors_need_init)
    {
        // Vectors for the first-order system
        d_sol1_vec =
            new SAMRAIVectorReal<NDIM, double>(d_object_name + "::sol1_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_sol1_vec->addComponent(d_U1_var, d_U1_scratch_idx, wgt_sc_idx, d_hier_sc_data_ops);
        d_sol1_vec->addComponent(d_P1_var, d_P1_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

        if (d_rhs1_vec) free_vector_components(*d_rhs1_vec);
        d_rhs1_vec = d_sol1_vec->cloneVector(d_object_name + "::rhs1_vec");

        // Vectors for the second-order system
        d_U2_scratch_vec =
            new SAMRAIVectorReal<NDIM, double>(d_object_name + "::U2_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_U2_scratch_vec->addComponent(d_U2_var, d_U2_scratch_idx, wgt_sc_idx, d_hier_sc_data_ops);

        d_P2_scratch_vec =
            new SAMRAIVectorReal<NDIM, double>(d_object_name + "::P2_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_P2_scratch_vec->addComponent(d_P2_var, d_P2_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

        if (d_U2_rhs_vec) free_vector_components(*d_U2_rhs_vec);
        if (d_P2_rhs_vec) free_vector_components(*d_P2_rhs_vec);

        d_U2_rhs_vec = d_U2_scratch_vec->cloneVector(d_object_name + "::U2_rhs_vec");
        d_P2_rhs_vec = d_P2_scratch_vec->cloneVector(d_object_name + "::P2_rhs_vec");

        d_sol2_vec =
            new SAMRAIVectorReal<NDIM, double>(d_object_name + "::sol2_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_sol2_vec->addComponent(d_U2_var, d_U2_scratch_idx, wgt_sc_idx, d_hier_sc_data_ops);
        d_sol2_vec->addComponent(d_P2_var, d_P2_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

        d_rhs2_vec =
            new SAMRAIVectorReal<NDIM, double>(d_object_name + "::rhs2_vec", d_hierarchy, coarsest_ln, finest_ln);
        const int U2_rhs_idx = d_U2_rhs_vec->getComponentDescriptorIndex(0);
        d_rhs2_vec->addComponent(d_U2_var, U2_rhs_idx, wgt_sc_idx, d_hier_sc_data_ops);
        const int P2_rhs_idx = d_P2_rhs_vec->getComponentDescriptorIndex(0);
        d_rhs2_vec->addComponent(d_P2_var, P2_rhs_idx, wgt_cc_idx, d_hier_cc_data_ops);

        for (const auto& nul2_vec : d_nul2_vecs)
        {
            if (nul2_vec) free_vector_components(*nul2_vec);
        }
        const int n_nul2_vecs = (has_pressure_nullspace ? 1 : 0) + (has_velocity_nullspace ? NDIM : 0);
        d_nul2_vecs.resize(n_nul2_vecs);

        for (const auto& U2_nul_vec : d_U2_nul_vecs)
        {
            if (U2_nul_vec) free_vector_components(*U2_nul_vec);
        }
        const int n_U2_nul_vecs = (has_velocity_nullspace ? NDIM : 0);
        d_U2_nul_vecs.resize(n_U2_nul_vecs);

        if (has_velocity_nullspace)
        {
            for (unsigned int k = 0; k < NDIM; ++k)
            {
                d_nul2_vecs[k] = d_sol2_vec->cloneVector(d_object_name + "::nul2_vec_U_" + std::to_string(k));
                d_nul2_vecs[k]->allocateVectorData(current_time);
                d_nul2_vecs[k]->setToScalar(0.0);
                d_U2_nul_vecs[k] = d_U2_scratch_vec->cloneVector(d_object_name + "::U2_nul_vec_U_" + std::to_string(k));
                d_U2_nul_vecs[k]->allocateVectorData(current_time);
                d_U2_nul_vecs[k]->setToScalar(0.0);
                for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<SideData<NDIM, double> > nul2_data =
                            patch->getPatchData(d_nul2_vecs[k]->getComponentDescriptorIndex(0));
                        nul2_data->getArrayData(k).fillAll(1.0);
                        Pointer<SideData<NDIM, double> > U2_nul_data =
                            patch->getPatchData(d_U2_nul_vecs[k]->getComponentDescriptorIndex(0));
                        U2_nul_data->getArrayData(k).fillAll(1.0);
                    }
                }
            }
        }

        if (has_pressure_nullspace)
        {
            d_nul2_vecs.back() = d_sol2_vec->cloneVector(d_object_name + "::nul2_vec_p");
            d_nul2_vecs.back()->allocateVectorData(current_time);
            d_hier_sc_data_ops->setToScalar(d_nul2_vecs.back()->getComponentDescriptorIndex(0), 0.0);
            d_hier_cc_data_ops->setToScalar(d_nul2_vecs.back()->getComponentDescriptorIndex(1), 1.0);
        }

        d_vectors_need_init = false;
    }

    // Setup boundary conditions objects for the second order system.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto U2_bc_coef = dynamic_cast<INSVCStaggeredVelocityBcCoef*>(d_U2_bc_coefs[d]);
        // U2_bc_coef->setStokesSpecifications(nullptr); // intentionally set to a nullptr to throw an error if used
        U2_bc_coef->setPhysicalBcCoefs(d_so_bc_coefs);
        U2_bc_coef->setSolutionTime(new_time);
        U2_bc_coef->setTimeInterval(current_time, new_time);
    }
    auto P2_bc_coef = dynamic_cast<INSVCStaggeredPressureBcCoef*>(d_P2_bc_coef);
    // P2_bc_coef->setStokesSpecifications(nullptr); // intentionally set to a nullptr to throw an error if used
    P2_bc_coef->setPhysicalBcCoefs(d_so_bc_coefs);
    P2_bc_coef->setSolutionTime(new_time);
    P2_bc_coef->setTimeInterval(current_time, new_time);
    P2_bc_coef->setViscosityInterpolationType(d_mu_vc_interp_type);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto U2_star_bc_coef = dynamic_cast<INSIntermediateVelocityBcCoef*>(d_U2_star_bc_coefs[d]);
        U2_star_bc_coef->setPhysicalBcCoefs(d_so_bc_coefs);
        U2_star_bc_coef->setSolutionTime(new_time);
        U2_star_bc_coef->setTimeInterval(current_time, new_time);
    }
    auto Phi_bc_coef = dynamic_cast<INSProjectionBcCoef*>(d_Phi_bc_coef);
    Phi_bc_coef->setPhysicalBcCoefs(d_so_bc_coefs);
    Phi_bc_coef->setSolutionTime(new_time);
    Phi_bc_coef->setTimeInterval(current_time, new_time);

    return;
} // preprocessOperatorsAndSolvers

void
AcousticStreamingHierarchyIntegrator::updateOperatorsAndSolvers(const double current_time,
                                                                const double new_time,
                                                                int cycle_num)
{
    // Compute the Brinkman contribution to the second order velocity operator.
    const bool has_brinkman_force = d_brinkman_force.size();
    if (has_brinkman_force)
    {
        d_hier_sc_data_ops->setToScalar(d_velocity_L_idx, 0.0);

        for (auto& brinkman_force : d_brinkman_force)
        {
            brinkman_force->demarcateBrinkmanZone(d_velocity_L_idx, new_time, cycle_num);
        }

        // Synchronize Brinkman coefficients.
        using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
        SynchronizationTransactionComponent L_synch_transaction =
            SynchronizationTransactionComponent(d_velocity_L_idx, "CONSERVATIVE_COARSEN");
        d_side_synch2_op->resetTransactionComponent(L_synch_transaction);
        d_side_synch2_op->synchronizeData(d_integrator_time);
        SynchronizationTransactionComponent default_synch2_transaction =
            SynchronizationTransactionComponent(d_U2_scratch_idx, d_U_coarsen_type);
        d_side_synch2_op->resetTransactionComponent(default_synch2_transaction);
    }

    // Modify the viscosity coefficients based on the backward Euler scheme
    const double K = 1.0;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        // Operate on a single level at a time
        d_hier_cc_data_ops->resetLevels(ln, ln);
        d_hier_sc_data_ops->resetLevels(ln, ln);
#if (NDIM == 2)
        d_hier_nc_data_ops->resetLevels(ln, ln);
#elif (NDIM == 3)
        d_hier_ec_data_ops->resetLevels(ln, ln);
#endif

        // C_sc = L(x,t^n+1)
        if (has_brinkman_force)
        {
            d_hier_sc_data_ops->copyData(d_velocity_C_idx, d_velocity_L_idx, /*interior_only*/ false);
        }
        // D_{ec,nc} = -K * mu
        // Lambda{cc} = -K * Lambda
#if (NDIM == 2)
        d_hier_nc_data_ops->scale(d_velocity_D_idx, -K, d_mu_interp_idx, /*interior_only*/ false);
#elif (NDIM == 3)
        d_hier_ec_data_ops->scale(d_velocity_D_idx, -K, d_mu_interp_idx, /*interior_only*/ false);
#endif
        d_hier_cc_data_ops->scale(d_velocity_D_cc_idx, -K, d_mu_scratch_idx, /*interior_only*/ false);

        if (d_lambda_var)
        {
            d_hier_cc_data_ops->scale(d_lambda_scratch_idx, -K, d_lambda_scratch_idx, /*interior_only*/ false);
        }
    }

    // Ensure that hierarchy operators operate on all levels from now on.
    d_hier_cc_data_ops->resetLevels(coarsest_ln, finest_ln);
    d_hier_sc_data_ops->resetLevels(coarsest_ln, finest_ln);
#if (NDIM == 2)
    d_hier_nc_data_ops->resetLevels(coarsest_ln, finest_ln);
#elif (NDIM == 3)
    d_hier_ec_data_ops->resetLevels(coarsest_ln, finest_ln);
#endif

    PoissonSpecifications U2_problem_coefs(d_object_name + "::U2_problem_coefs");
    PoissonSpecifications P2_problem_coefs(d_object_name + "::P2_problem_coefs");
    d_vc_stokes_op_spec.reset();
    d_vc_projection_pc_spec.reset();

    if (has_brinkman_force)
    {
        U2_problem_coefs.setCPatchDataId(d_velocity_C_idx);
        d_vc_stokes_op_spec.d_C_idx = d_velocity_C_idx;
        d_vc_stokes_op_spec.d_C_is_const = false;
        d_vc_projection_pc_spec.d_steady_state = false;
        d_vc_projection_pc_spec.d_theta = 1.0;
    }
    else
    {
        U2_problem_coefs.setCZero(); // pure Stokes problem
        d_vc_stokes_op_spec.d_C_const = 0.0;
        d_vc_stokes_op_spec.d_C_is_const = true;
        d_vc_projection_pc_spec.d_steady_state = true;
        d_vc_projection_pc_spec.d_theta = 0.0;
    }

    U2_problem_coefs.setDPatchDataId(d_velocity_D_idx);
    d_vc_stokes_op_spec.d_D_idx = d_velocity_D_idx;
    d_vc_stokes_op_spec.d_D_is_const = false;

    if (d_lambda_var)
    {
        d_vc_stokes_op_spec.d_L_idx = d_lambda_scratch_idx;
        d_vc_stokes_op_spec.d_L_is_const = false;
    }

    // Set the coefficients for the pressure solver within the projection preconditioner
    // D_sc = -rho/chi inside Brinkman zone and -rho outside.
    d_hier_sc_data_ops->setToScalar(d_velocity_L_idx, 1.0, /*interior_only*/ false);
    for (auto& brinkman_force : d_brinkman_force)
    {
        brinkman_force->demarcateBrinkmanZone(d_velocity_L_idx, new_time, cycle_num);
    }
    d_hier_sc_data_ops->reciprocal(d_projection_D_idx,
                                   d_velocity_L_idx,
                                   /*interior_only*/ false);
    d_hier_sc_data_ops->scale(d_projection_D_idx,
                              -1.0,
                              d_projection_D_idx,
                              /*interior_only*/ false);
    d_hier_sc_data_ops->multiply(d_pressure_D_idx, d_rho_scratch_idx, d_projection_D_idx, /*interior_only*/ false);

    // Synchronize pressure patch data coefficient
    using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
    SynchronizationTransactionComponent p_coef_synch_transaction =
        SynchronizationTransactionComponent(d_pressure_D_idx, "CONSERVATIVE_COARSEN");
    d_side_synch2_op->resetTransactionComponent(p_coef_synch_transaction);
    d_side_synch2_op->synchronizeData(d_integrator_time);
    SynchronizationTransactionComponent default_synch2_transaction =
        SynchronizationTransactionComponent(d_U2_scratch_idx, d_U_coarsen_type);
    d_side_synch2_op->resetTransactionComponent(default_synch2_transaction);

    P2_problem_coefs.setCZero();
    P2_problem_coefs.setDPatchDataId(d_pressure_D_idx);

    // Set the coefficients used to compute density weighted divergence of velocity for the RHS of pressure problem
    // as well as to correct velocity after pressure is computed within the preconditioner.
    d_vc_projection_pc_spec.d_div_coef_idx = d_rho_linear_op_idx;
    d_vc_projection_pc_spec.d_D_idx = d_projection_D_idx;
    d_vc_projection_pc_spec.d_D_is_const = false;

    // Set the cell centered "local viscosity" patch data index to be used within the projection preconditioner.
    d_vc_projection_pc_spec.d_mu_cc_idx = d_velocity_D_cc_idx;

    // Set the coefficient to compute density weighted divergence of velocity in the VCStaggeredStokesOperator.
    d_vc_stokes_op_spec.d_div_coef_idx = d_rho_linear_op_idx;

    // Ensure that solver components are appropriately reinitialized at the
    // correct intervals
    const bool precond_reinit = d_integrator_step % d_precond_reinit_interval == 0;
    if (precond_reinit)
    {
        d_velocity_solver_needs_init = true;
        d_pressure_solver_needs_init = true;
        d_stokes_solver_needs_init = true;
    }

    // Setup subdomain solvers.
    const bool has_velocity_nullspace = d_normalize_velocity;
    const bool has_pressure_nullspace = d_normalize_pressure;

    d_velocity_solver->setPoissonSpecifications(U2_problem_coefs);
    d_velocity_solver->setPhysicalBcCoefs(d_U2_star_bc_coefs);
    d_velocity_solver->setSolutionTime(new_time);
    d_velocity_solver->setTimeInterval(current_time, new_time);
    if (d_velocity_solver_needs_init)
    {
        if (d_enable_logging)
            plog << d_object_name
                 << "::updateOperatorsAndSolvers`(): initializing "
                    "velocity subdomain solver"
                 << std::endl;
        auto p_velocity_solver = dynamic_cast<LinearSolver*>(d_velocity_solver.getPointer());
        if (p_velocity_solver)
        {
            p_velocity_solver->setInitialGuessNonzero(false);
            if (has_velocity_nullspace) p_velocity_solver->setNullspace(false, d_U2_nul_vecs);
        }
        d_velocity_solver->initializeSolverState(*d_U2_scratch_vec, *d_U2_rhs_vec);
        d_velocity_solver_needs_init = false;
    }

    d_pressure_solver->setPoissonSpecifications(P2_problem_coefs);
    d_pressure_solver->setPhysicalBcCoef(d_Phi_bc_coef);
    d_pressure_solver->setSolutionTime(new_time);
    d_pressure_solver->setTimeInterval(current_time, new_time);
    if (d_pressure_solver_needs_init)
    {
        if (d_enable_logging)
            plog << d_object_name
                 << "::updateOperatorsAndSolvers(): initializing "
                    "pressure subdomain solver"
                 << std::endl;
        auto p_pressure_solver = dynamic_cast<LinearSolver*>(d_pressure_solver.getPointer());
        if (p_pressure_solver)
        {
            p_pressure_solver->setInitialGuessNonzero(false);
            if (has_pressure_nullspace) p_pressure_solver->setNullspace(true);
        }
        d_pressure_solver->initializeSolverState(*d_P2_scratch_vec, *d_P2_rhs_vec);
        d_pressure_solver_needs_init = false;
    }

    // Setup Stokes solver.
    d_stokes_solver->setPhysicalBcCoefs(d_U2_bc_coefs, d_P2_bc_coef);
    d_stokes_solver->setPhysicalBoundaryHelper(d_bc_helper);
    d_stokes_solver->setSolutionTime(new_time);
    d_stokes_solver->setTimeInterval(current_time, new_time);
    d_stokes_solver->setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);

    auto p_stokes_linear_solver = dynamic_cast<LinearSolver*>(d_stokes_solver.getPointer());
    auto p_stokes_krylov_solver = dynamic_cast<KrylovLinearSolver*>(p_stokes_linear_solver);
    Pointer<VCStaggeredStokesOperator> p_vc_stokes_op = p_stokes_krylov_solver->getOperator();
    p_vc_stokes_op->setProblemSpecification(&d_vc_stokes_op_spec);

    auto p_stokes_block_pc =
        dynamic_cast<StaggeredStokesBlockPreconditioner*>(p_stokes_krylov_solver->getPreconditioner().getPointer());
    auto p_vc_stokes_proj_pc = dynamic_cast<VCStaggeredStokesProjectionPreconditioner*>(
        p_stokes_krylov_solver->getPreconditioner().getPointer());

    p_stokes_block_pc->setPhysicalBcCoefs(d_U2_star_bc_coefs, d_Phi_bc_coef);
    p_stokes_block_pc->setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);

    p_vc_stokes_proj_pc->setProblemSpecification(&d_vc_projection_pc_spec);

    if (d_stokes_solver_needs_init)
    {
        if (d_enable_logging)
            plog << d_object_name
                 << "::updateOperatorsAndSolvers(): initializing "
                    "incompressible Stokes solver"
                 << std::endl;

        p_stokes_linear_solver->setInitialGuessNonzero(true);
        if (has_velocity_nullspace || has_pressure_nullspace) p_stokes_linear_solver->setNullspace(false, d_nul2_vecs);

        d_stokes_solver->initializeSolverState(*d_sol2_vec, *d_rhs2_vec);
        d_stokes_solver_needs_init = false;
    }

    // Setup the first-order solver. Because the coefficients (viscosity etc.) have
    // changed we need to rebuild the matrix.
    d_first_order_solver->deallocateSolverState();

    d_first_order_solver->setSoundSpeed(d_sound_speed);
    d_first_order_solver->setAcousticAngularFrequency(d_acoustic_freq);
    d_first_order_solver->setBoundaryConditionCoefficients(d_U1_bc_coefs);
    d_first_order_solver->setMassDensityPatchDataIndex(d_rho_scratch_idx);
    d_first_order_solver->setShearViscosityPatchDataIndex(d_velocity_D_idx);
    d_first_order_solver->setViscosityInterpolationType(d_mu_vc_interp_type);
    d_first_order_solver->setBulkViscosityPatchDataIndex(d_lambda_scratch_idx);

    d_first_order_solver->initializeSolverState(*d_sol1_vec, *d_rhs1_vec);

} // updateOperatorsAndSolvers

void
AcousticStreamingHierarchyIntegrator::setupSolverVectors(const Pointer<SAMRAIVectorReal<NDIM, double> >& sol1_vec,
                                                         const Pointer<SAMRAIVectorReal<NDIM, double> >& rhs1_vec,
                                                         const Pointer<SAMRAIVectorReal<NDIM, double> >& sol2_vec,
                                                         const Pointer<SAMRAIVectorReal<NDIM, double> >& rhs2_vec,
                                                         const double current_time,
                                                         const double new_time,
                                                         const int cycle_num)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Account for body forcing terms.
    if (d_F1_fcn)
    {
        d_F1_fcn->setDataOnPatchHierarchy(d_F1_scratch_idx, d_F1_var, d_hierarchy, new_time);
        d_hier_sc_data_ops->add(
            rhs1_vec->getComponentDescriptorIndex(0), rhs1_vec->getComponentDescriptorIndex(0), d_F1_scratch_idx);
    }
    if (d_F2_fcn)
    {
        d_F2_fcn->setDataOnPatchHierarchy(d_F2_scratch_idx, d_F2_var, d_hierarchy, new_time);
        d_hier_sc_data_ops->add(
            rhs2_vec->getComponentDescriptorIndex(0), rhs2_vec->getComponentDescriptorIndex(0), d_F2_scratch_idx);
    }

    // Account for mass source/sink terms.
    if (d_Q1_fcn)
    {
        d_Q1_fcn->setDataOnPatchHierarchy(d_Q1_new_idx, d_Q1_var, d_hierarchy, new_time);
        d_hier_cc_data_ops->add(
            rhs1_vec->getComponentDescriptorIndex(1), rhs1_vec->getComponentDescriptorIndex(1), d_Q1_new_idx);
    }
    if (d_Q2_fcn)
    {
        d_Q2_fcn->setDataOnPatchHierarchy(d_Q2_new_idx, d_Q2_var, d_hierarchy, new_time);
        d_hier_cc_data_ops->add(
            rhs2_vec->getComponentDescriptorIndex(1), rhs2_vec->getComponentDescriptorIndex(1), d_Q2_new_idx);
    }

    // Compute and add the Brinkman term to the RHS vector.
    const bool has_brinkman = d_brinkman_force.size();
    if (has_brinkman)
    {
        d_hier_sc_data_ops->setToScalar(d_velocity_L_idx, 0.0);
        for (auto& brinkman_force : d_brinkman_force)
        {
            brinkman_force->computeBrinkmanVelocity(d_velocity_L_idx, new_time, cycle_num);
        }
        d_hier_sc_data_ops->add(
            rhs2_vec->getComponentDescriptorIndex(0), rhs2_vec->getComponentDescriptorIndex(0), d_velocity_L_idx);
    }

    // Set solution components to equal most recent approximations to u(n+1) and
    // p(n+1).
    d_hier_sc_data_ops->copyData(sol1_vec->getComponentDescriptorIndex(0), d_U1_new_idx);
    d_hier_cc_data_ops->copyData(sol1_vec->getComponentDescriptorIndex(1), d_P1_new_idx);
    d_hier_sc_data_ops->copyData(sol2_vec->getComponentDescriptorIndex(0), d_U2_new_idx);
    d_hier_cc_data_ops->copyData(sol2_vec->getComponentDescriptorIndex(1), d_P2_new_idx);

    // Synchronize solution and right-hand-side data before solve for the second order system.
    using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;

    SynchronizationTransactionComponent sol1_synch_transaction =
        SynchronizationTransactionComponent(sol1_vec->getComponentDescriptorIndex(0), d_U_coarsen_type);
    d_side_synch1_op->resetTransactionComponent(sol1_synch_transaction);
    d_side_synch1_op->synchronizeData(current_time);
    SynchronizationTransactionComponent rhs1_synch_transaction =
        SynchronizationTransactionComponent(rhs1_vec->getComponentDescriptorIndex(0), d_F_coarsen_type);
    d_side_synch1_op->resetTransactionComponent(rhs1_synch_transaction);
    d_side_synch1_op->synchronizeData(current_time);
    SynchronizationTransactionComponent default_synch1_transaction =
        SynchronizationTransactionComponent(d_U1_scratch_idx, d_U_coarsen_type);
    d_side_synch1_op->resetTransactionComponent(default_synch1_transaction);

    SynchronizationTransactionComponent sol2_synch_transaction =
        SynchronizationTransactionComponent(sol2_vec->getComponentDescriptorIndex(0), d_U_coarsen_type);
    d_side_synch2_op->resetTransactionComponent(sol2_synch_transaction);
    d_side_synch2_op->synchronizeData(current_time);
    SynchronizationTransactionComponent rhs2_synch_transaction =
        SynchronizationTransactionComponent(rhs2_vec->getComponentDescriptorIndex(0), d_F_coarsen_type);
    d_side_synch2_op->resetTransactionComponent(rhs2_synch_transaction);
    d_side_synch2_op->synchronizeData(current_time);
    SynchronizationTransactionComponent default_synch2_transaction =
        SynchronizationTransactionComponent(d_U2_scratch_idx, d_U_coarsen_type);
    d_side_synch2_op->resetTransactionComponent(default_synch2_transaction);

    return;
} // setupSolverVectors

void
AcousticStreamingHierarchyIntegrator::resetSolverVectors(const Pointer<SAMRAIVectorReal<NDIM, double> >& sol1_vec,
                                                         const Pointer<SAMRAIVectorReal<NDIM, double> >& rhs1_vec,
                                                         const Pointer<SAMRAIVectorReal<NDIM, double> >& sol2_vec,
                                                         const Pointer<SAMRAIVectorReal<NDIM, double> >& rhs2_vec,
                                                         const double current_time,
                                                         const double /*new_time*/,
                                                         const int /*cycle_num*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Synchronize solution data after solve.
    using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
    SynchronizationTransactionComponent sol1_synch_transaction =
        SynchronizationTransactionComponent(sol1_vec->getComponentDescriptorIndex(0), d_U_coarsen_type);
    d_side_synch1_op->synchronizeData(current_time);
    SynchronizationTransactionComponent default_synch1_transaction =
        SynchronizationTransactionComponent(d_U1_scratch_idx, d_U_coarsen_type);
    d_side_synch1_op->resetTransactionComponent(default_synch1_transaction);

    SynchronizationTransactionComponent sol2_synch_transaction =
        SynchronizationTransactionComponent(sol2_vec->getComponentDescriptorIndex(0), d_U_coarsen_type);
    d_side_synch2_op->synchronizeData(current_time);
    SynchronizationTransactionComponent default_synch2_transaction =
        SynchronizationTransactionComponent(d_U2_scratch_idx, d_U_coarsen_type);
    d_side_synch2_op->resetTransactionComponent(default_synch2_transaction);

    // Pull out solution components.
    d_hier_sc_data_ops->copyData(d_U1_new_idx, sol1_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_P1_new_idx, sol1_vec->getComponentDescriptorIndex(1));
    d_hier_sc_data_ops->copyData(d_U2_new_idx, sol2_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_P2_new_idx, sol2_vec->getComponentDescriptorIndex(1));

    // Remove body force terms
    if (d_F1_fcn)
    {
        d_hier_sc_data_ops->subtract(
            rhs1_vec->getComponentDescriptorIndex(0), rhs1_vec->getComponentDescriptorIndex(0), d_F1_scratch_idx);
        d_hier_sc_data_ops->copyData(d_F1_new_idx, d_F1_scratch_idx);
    }
    d_hier_sc_data_ops->subtract(
        rhs2_vec->getComponentDescriptorIndex(0), rhs2_vec->getComponentDescriptorIndex(0), d_velocity_L_idx);
    if (d_F2_fcn)
    {
        d_hier_sc_data_ops->subtract(
            rhs2_vec->getComponentDescriptorIndex(0), rhs2_vec->getComponentDescriptorIndex(0), d_F2_scratch_idx);
        d_hier_sc_data_ops->copyData(d_F2_new_idx, d_F2_scratch_idx);
    }

    // Remove mass source/sink terms.
    if (d_Q1_fcn)
    {
        d_hier_cc_data_ops->subtract(
            rhs1_vec->getComponentDescriptorIndex(1), rhs1_vec->getComponentDescriptorIndex(1), d_Q1_new_idx);
    }
    if (d_Q2_fcn)
    {
        d_hier_cc_data_ops->subtract(
            rhs2_vec->getComponentDescriptorIndex(1), rhs2_vec->getComponentDescriptorIndex(1), d_Q2_new_idx);
    }

    return;
} // resetSolverVectors

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////