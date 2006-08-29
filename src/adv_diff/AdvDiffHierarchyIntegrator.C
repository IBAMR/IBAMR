//
// AdvectionDiffusionHierarchyIntegrator.C
//
// Created on 17 Mar 2004
//         by Boyce Griffith (boyce@bigboy.speakeasy.net).
//
// Last modified: <30.Jan.2006 22:23:40 boyce@boyce.cims.nyu.edu>
//

#include "AdvectionDiffusionHierarchyIntegrator.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

// C++ STDLIB INCLUDES
//
#include <limits>

// STL INCLUDES
//
#include <set>

// SAMRAI-tools INCLUDES
//
#include "AbstractKrylovLinearSolver.h"
#include "CCPoissonAMRHypreSolver.h"
#include "FACPreconditionerALSWrapper.h"
#include "PETScKrylovLinearSolver.h"

// SAMRAI INCLUDES
//
#include "CartesianGridGeometry.h"
#include "CellDataFactory.h"
#include "CoarsenOperator.h"
#include "HierarchyDataOpsManager.h"
#include "RefineOperator.h"
#include "SimpleCellRobinBcCoefs.h"
#include "VariableDatabase.h"
#include "tbox/NullDatabase.h"
#include "tbox/RestartManager.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// INLINE ///////////////////////////////////////

//#ifdef DEBUG_NO_INLINE
//#include "AdvectionDiffusionHierarchyIntegrator.I"
//#endif

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
    // Timers.
    static tbox::Pointer<tbox::Timer> t_initialize_hierarchy_integrator;
    static tbox::Pointer<tbox::Timer> t_initialize_hierarchy;
    static tbox::Pointer<tbox::Timer> t_advance_hierarchy;
    static tbox::Pointer<tbox::Timer> t_rebalance_coarsest_level;
    static tbox::Pointer<tbox::Timer> t_regrid_hierarchy;
    static tbox::Pointer<tbox::Timer> t_integrate_hierarchy;
    static tbox::Pointer<tbox::Timer> t_synchronize_hierarchy;
    static tbox::Pointer<tbox::Timer> t_synchronize_new_levels;
    static tbox::Pointer<tbox::Timer> t_reset_time_dependent_data;
    static tbox::Pointer<tbox::Timer> t_reset_data_to_preadvance_state;
    static tbox::Pointer<tbox::Timer> t_initialize_level_data;
    static tbox::Pointer<tbox::Timer> t_reset_hierarchy_configuration;
    static tbox::Pointer<tbox::Timer> t_apply_gradient_detector;
    static tbox::Pointer<tbox::Timer> t_put_to_database;

    // Number of ghosts cells used for each variable quantity.
    static const int CELLG = 1;

    // Version of AdvectionDiffusionHierarchyIntegrator restart file data
    static const int ADVECTION_DIFFUSION_HIERARCHY_INTEGRATOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvectionDiffusionHierarchyIntegrator::AdvectionDiffusionHierarchyIntegrator(
    const string& object_name,
    tbox::Pointer<tbox::Database> input_db,
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    tbox::Pointer<GodunovAdvector> explicit_predictor,
    bool register_for_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
    assert(!hierarchy.isNull());
    assert(!explicit_predictor.isNull());
#endif
    d_object_name = object_name;
    d_registered_for_restart = register_for_restart;

    d_hierarchy = hierarchy;
    
    if (d_registered_for_restart)
    {
        tbox::RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }
    
    // Set some default values.
    d_start_time = 0.0;
    d_end_time = numeric_limits<double>::max();
    d_grow_dt = 2.0;
    d_max_integrator_steps = numeric_limits<int>::max();
    d_regrid_interval = 1;
    d_old_dt = -1.0;
    d_integrator_time = numeric_limits<double>::quiet_NaN();
    d_integrator_step = numeric_limits<int>::max();

    d_solver_package = "PETSc";
    d_max_iterations = 25;
    d_abs_residual_tol = 1.0e-30;
    d_rel_residual_tol = 1.0e-8;
    
    d_is_initialized = false;

    d_do_log = false;

    d_maintain_helmholtz3_solvers = false;
    d_maintain_helmholtz4_solvers = false;
    
    // Initialize object with data read from the input and restart
    // databases.
    bool from_restart = tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Initialize the algs::HyperbolicPatchStrategy<NDIM> and
    // algs::HyperbolicLevelIntegrator<NDIM> objects used to provide the numerical
    // routines for explicitly integrating the advective terms.
    tbox::Pointer<tbox::Database> hyp_patch_strategy_input_db;
    if (input_db->keyExists("AdvDiffHypPatchStrategy"))
    {
        hyp_patch_strategy_input_db = input_db->
            getDatabase("AdvDiffHypPatchStrategy");
    }
    else
    {
        hyp_patch_strategy_input_db = new tbox::NullDatabase();
    }
    
    d_hyp_patch_strategy =
        new AdvDiffHypPatchStrategy(object_name+"::AdvDiffHypPatchStrategy",
                                    hyp_patch_strategy_input_db,
                                    explicit_predictor,
                                    d_hierarchy->getGridGeometry(),
                                    d_registered_for_restart);
    
    tbox::Pointer<tbox::Database> hyp_level_integrator_input_db;
    if (input_db->keyExists("HyperbolicLevelIntegrator"))
    {
        hyp_level_integrator_input_db = input_db->
            getDatabase("HyperbolicLevelIntegrator");
    }
    else
    {
        hyp_level_integrator_input_db = new tbox::NullDatabase();
    }

    const bool use_time_refinement = false;
    d_hyp_level_integrator = new algs::HyperbolicLevelIntegrator<NDIM>(
        object_name+"::HyperbolicLevelIntegrator",
        hyp_level_integrator_input_db,
        d_hyp_patch_strategy,
        register_for_restart,
        use_time_refinement);

    // Get initialization data for the FAC ops and FAC preconditioners.
    if (input_db->keyExists("FACOps"))
    {
        d_fac_ops_db = input_db->getDatabase("FACOps");
    }
    if (input_db->keyExists("FACPreconditioners"))
    {
        d_fac_pcs_db = input_db->getDatabase("FACPreconditioners");
    }
    if (input_db->keyExists("CCPoissonAMRHypreSolver"))
    {
        d_hypre_solver_db = input_db->getDatabase("CCPoissonAMRHypreSolver");
    }
    
    // Obtain the Hierarchy data operations objects.    
    math::HierarchyDataOpsManager<NDIM>* hier_ops_manager =
        math::HierarchyDataOpsManager<NDIM>::getManager();
    
    tbox::Pointer<pdat::CellVariable<NDIM,double> > cc_var = new pdat::CellVariable<NDIM,double>(
        "cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(
        cc_var, hierarchy);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::advanceHierarchy()");
        t_rebalance_coarsest_level = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::rebalanceCoarsestLevel()");
        t_regrid_hierarchy = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::regridHierarchy()");
        t_integrate_hierarchy = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::integrateHierarchy()");
        t_synchronize_hierarchy = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::resetTimeDependentData()");
        t_reset_data_to_preadvance_state = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::resetDataToPreadvanceState()");
        t_initialize_level_data = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::AdvectionDiffusionHierarchyIntegrator::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// AdvectionDiffusionHierarchyIntegrator

AdvectionDiffusionHierarchyIntegrator::~AdvectionDiffusionHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    // Deallocate all solver components.
    d_helmholtz1_solvers.clear();
    d_helmholtz1_fac_pcs.clear();
    d_helmholtz1_fac_ops.clear();
    d_helmholtz1_ops.clear();
    d_helmholtz1_specs.clear();
    d_helmholtz1_bc_coefs.clear();
    
    d_helmholtz2_solvers.clear();
    d_helmholtz2_fac_pcs.clear();
    d_helmholtz2_fac_ops.clear();
    d_helmholtz2_ops.clear();
    d_helmholtz2_specs.clear();
    d_helmholtz2_bc_coefs.clear();

    d_helmholtz3_solvers.clear();
    d_helmholtz3_fac_pcs.clear();
    d_helmholtz3_fac_ops.clear();
    d_helmholtz3_ops.clear();
    d_helmholtz3_specs.clear();
    d_helmholtz3_bc_coefs.clear();

    d_helmholtz4_solvers.clear();
    d_helmholtz4_fac_pcs.clear();
    d_helmholtz4_fac_ops.clear();
    d_helmholtz4_ops.clear();
    d_helmholtz4_specs.clear();
    d_helmholtz4_bc_coefs.clear();

    return;
}// ~AdvectionDiffusionHierarchyIntegrator

///
///  The following routines:
///
///      registerAdvectedAndDiffusedQuantity(),
///      registerAdvectedAndDiffusedQuantityWithSourceTerm(),
///      registerAdvectionVelocity(),
///      registerVisItDataWriter()
///
///  allow the specification of quantities to be advected and
///  diffused.
///

void AdvectionDiffusionHierarchyIntegrator::registerAdvectedAndDiffusedQuantity(
    tbox::Pointer<pdat::CellVariable<NDIM,double> > Q_var,
    const double Q_mu,
    const bool conservation_form,
    tbox::Pointer<SetDataStrategy> Q_init,
    tbox::Pointer<PhysicalBCDataStrategy> Q_bc,
    tbox::Pointer<pdat::FaceVariable<NDIM,double> > grad_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!Q_var.isNull());
    assert(Q_mu >= 0.0);
#endif
    d_Q_vars .push_back(Q_var);
    d_Q_inits.push_back(Q_init);
    d_Q_bcs  .push_back(Q_bc);
    d_Q_mus  .push_back(Q_mu);

    d_grad_vars.push_back(grad_var);
    
    d_F_vars.push_back(NULL);
    d_F_sets.push_back(NULL);
    
    tbox::Pointer<pdat::CellDataFactory<NDIM,double> > Q_factory =
        Q_var->getPatchDataFactory();
    const int Q_depth = Q_factory->getDefaultDepth();
    tbox::Pointer<pdat::CellVariable<NDIM,double> > Psi_var =
        new pdat::CellVariable<NDIM,double>(Q_var->getName()+"::Psi",Q_depth);
    
    d_Psi_vars.push_back(Psi_var);
    
    d_hyp_patch_strategy->registerAdvectedQuantityWithSourceTerm(
        Q_var, Psi_var, conservation_form, Q_init, Q_bc,
        tbox::Pointer<SetDataStrategy>(NULL), grad_var);
    
    return;
}// registerAdvectedAndDiffusedQuantity

void AdvectionDiffusionHierarchyIntegrator::registerAdvectedAndDiffusedQuantityWithSourceTerm(
    tbox::Pointer<pdat::CellVariable<NDIM,double> > Q_var,
    const double Q_mu,
    tbox::Pointer<pdat::CellVariable<NDIM,double> > F_var,
    const bool conservation_form,
    tbox::Pointer<SetDataStrategy> Q_init,
    tbox::Pointer<PhysicalBCDataStrategy> Q_bc,
    tbox::Pointer<SetDataStrategy> F_set,
    tbox::Pointer<pdat::FaceVariable<NDIM,double> > grad_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!Q_var.isNull());
    assert(Q_mu >= 0.0);
    assert(!F_var.isNull());
#endif
    d_Q_vars .push_back(Q_var);
    d_Q_inits.push_back(Q_init);
    d_Q_bcs  .push_back(Q_bc);
    d_Q_mus  .push_back(Q_mu);
    
    d_grad_vars.push_back(grad_var);

    d_F_vars.push_back(F_var);
    d_F_sets.push_back(F_set);
    
    tbox::Pointer<pdat::CellDataFactory<NDIM,double> > Q_factory =
        Q_var->getPatchDataFactory();
    const int Q_depth = Q_factory->getDefaultDepth();

#ifdef DEBUG_CHECK_ASSERTIONS
    tbox::Pointer<pdat::CellDataFactory<NDIM,double> > F_factory =
        F_var->getPatchDataFactory();
    assert(Q_depth == F_factory->getDefaultDepth());
#endif

    tbox::Pointer<pdat::CellVariable<NDIM,double> > Psi_var =
        new pdat::CellVariable<NDIM,double>(Q_var->getName()+"::Psi",Q_depth);

    d_Psi_vars.push_back(Psi_var);
    
    d_hyp_patch_strategy->registerAdvectedQuantityWithSourceTerm(
        Q_var, Psi_var, conservation_form, Q_init, Q_bc,
        tbox::Pointer<SetDataStrategy>(NULL), grad_var);
    
    return;
}// registerAdvectedAndDiffusedQuantityWithSourceTerm

void AdvectionDiffusionHierarchyIntegrator::registerAdvectionVelocity(
    tbox::Pointer<pdat::FaceVariable<NDIM,double> > u_var,
    const bool u_is_div_free,
    tbox::Pointer<SetDataStrategy> u_set)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!u_var.isNull());
#endif
    d_hyp_patch_strategy->registerAdvectionVelocity(
        u_var,u_is_div_free,u_set);
    d_u_var = u_var;
    d_u_set = u_set;
    d_u_is_div_free = u_is_div_free;
    return;
}// registerAdvectionVelocity

void AdvectionDiffusionHierarchyIntegrator::registerConvergenceMonitor(
    tbox::Pointer<ConvergenceMonitor> monitor)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!(monitor.isNull()));
#endif
    d_convergence_monitor = monitor;
    return;
}// registerConvergenceMonitor

void AdvectionDiffusionHierarchyIntegrator::registerVisItDataWriter(
    tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!(visit_writer.isNull()));
#endif
    d_hyp_patch_strategy->registerVisItDataWriter(visit_writer);
    return;
}// registerVisItDataWriter

///
///  The following routines:
///
///      getHierarchyMathOps(),
///      setHierarchyMathOps(),
///      isManagingHierarchyMathOps()
///
///  allow for the sharing of a single HierarchyMathOps object between
///  mutiple HierarchyIntegrator objects.
///

tbox::Pointer<HierarchyMathOps> AdvectionDiffusionHierarchyIntegrator::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!d_hier_math_ops.isNull());
#endif
    return(d_hier_math_ops);
}// getHierarchyMathOps

void AdvectionDiffusionHierarchyIntegrator::setHierarchyMathOps(
    tbox::Pointer<HierarchyMathOps> hier_math_ops,
    const bool manage_ops)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hier_math_ops.isNull());
#endif
    d_hier_math_ops = hier_math_ops;
    d_is_managing_hier_math_ops = manage_ops;
    return;
}// setHierarchyMathOps

bool AdvectionDiffusionHierarchyIntegrator::isManagingHierarchyMathOps() const
{
    return(d_is_managing_hier_math_ops);
}// isManagingHierarchyMathOps

///
///  The following routines:
///
///      getHelmholtzSpecs(),
///      getHelmholtzBcCoefs(),
///      getHelmholtzSolvers(),
///      maintainExtraSolvers()
///
///  allow other objects to access the Helmholtz solvers and related
///  data used by this integrator.
///  

vector<const solv::PoissonSpecifications*> AdvectionDiffusionHierarchyIntegrator::getHelmholtzSpecs(
    const double mu)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_helmholtz1_specs.count(mu) > 0);
    assert(d_helmholtz2_specs.count(mu) > 0);
    assert(d_helmholtz3_specs.count(mu) > 0);
    assert(d_helmholtz4_specs.count(mu) > 0);
#endif
    vector<const solv::PoissonSpecifications*> return_vec;
    return_vec.push_back(d_helmholtz1_specs[mu]);
    return_vec.push_back(d_helmholtz2_specs[mu]);
    return_vec.push_back(d_helmholtz3_specs[mu]);
    return_vec.push_back(d_helmholtz4_specs[mu]);

    return(return_vec);
}// getHelmholtzSpecs

vector<const solv::RobinBcCoefStrategy<NDIM>*> AdvectionDiffusionHierarchyIntegrator::getHelmholtzBcCoefs(
    const double mu)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_helmholtz1_bc_coefs.count(mu) > 0);
    assert(d_helmholtz2_bc_coefs.count(mu) > 0);
    assert(d_helmholtz3_bc_coefs.count(mu) > 0);
    assert(d_helmholtz4_bc_coefs.count(mu) > 0);
#endif
    vector<const solv::RobinBcCoefStrategy<NDIM>*> return_vec;
    return_vec.push_back(d_helmholtz1_bc_coefs[mu]);
    return_vec.push_back(d_helmholtz2_bc_coefs[mu]);
    return_vec.push_back(d_helmholtz3_bc_coefs[mu]);
    return_vec.push_back(d_helmholtz4_bc_coefs[mu]);
    
    return(return_vec);
}// getHelmholtzBcCoefs

vector<tbox::Pointer<AbstractLinearSolver> > AdvectionDiffusionHierarchyIntegrator::getHelmholtzSolvers(
    const double mu)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_helmholtz1_solvers.count(mu) > 0);
    assert(d_helmholtz2_solvers.count(mu) > 0);
#endif
    vector<tbox::Pointer<AbstractLinearSolver> > return_vec;
    return_vec.push_back(d_helmholtz1_solvers[mu]);
    return_vec.push_back(d_helmholtz2_solvers[mu]);
    return_vec.push_back(d_helmholtz3_solvers[mu]);
    return_vec.push_back(d_helmholtz4_solvers[mu]);
    
    return(return_vec);
}// getHelmholtzSolvers

void AdvectionDiffusionHierarchyIntegrator::maintainExtraSolvers(
    const int coeff)
{
    if (coeff == 3)
    {
        d_maintain_helmholtz3_solvers = true;
    }
    else if (coeff == 4)
    {
        d_maintain_helmholtz4_solvers = true;
    }
    else
    {
        TBOX_ERROR(d_object_name << "::maintainExtraSolvers()\n"
                   << "  unknown solver coefficient: " << coeff << "\n");
    }
    return;
}// maintainExtraSolvers

///
///  The following routines:
///
///      initializeHierarchyIntegrator(),
///      initializeHierarchy(),
///      advanceHierarchy(),
///      atRegridPoint(),
///      getIntegratorTime(),
///      getStartTime(),
///      getEndTime(),
///      getIntegratorStep(),
///      getMaxIntegratorSteps(),
///      stepsRemaining(),
///      getPatchHierarchy(),
///      getGriddingAlgorithm(),
///      getHyperbolicLevelIntegrator(),
///      getHyperbolicPatchStrategy()
///
///  allow the AdvectionDiffusionHierarchyIntegrator to be used as a
///  hierarchy integrator.
///

void AdvectionDiffusionHierarchyIntegrator::initializeHierarchyIntegrator(
    tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();
    
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Setup the tag buffer.
    if (d_using_default_tag_buffer)
    {
        d_tag_buffer.resizeArray(d_gridding_alg->getMaxLevels());
        for (int i = 0; i < d_gridding_alg->getMaxLevels(); ++i)
        {
            d_tag_buffer[i] = d_regrid_interval;
        }
    }
    else
    {
        if (d_tag_buffer.getSize() < d_gridding_alg->getMaxLevels())
        {
            int tsize = d_tag_buffer.getSize();
            d_tag_buffer.resizeArray(d_gridding_alg->getMaxLevels());
            for (int i = tsize; i < d_gridding_alg->getMaxLevels(); ++i)
            {
                d_tag_buffer[i] = d_tag_buffer[tsize-1];
            }
        }
    }
    
    // Register forcing term data.
    const hier::IntVector<NDIM> cell_ghosts = CELLG;
    
    typedef vector<tbox::Pointer<pdat::CellVariable<NDIM,double> > > CellVariableVector;
    for (CellVariableVector::iterator it = d_F_vars.begin();
         it != d_F_vars.end(); ++it)
    {
        // Advected and diffused quantities do not necessarily have
        // source terms.
        tbox::Pointer<pdat::CellVariable<NDIM,double> > F_var = *it;

        if (!F_var.isNull())
        {
            d_hyp_level_integrator->
                registerVariable(F_var, cell_ghosts,
                                 algs::HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                                 d_hierarchy->getGridGeometry(),
                                 "CONSERVATIVE_COARSEN",
                                 "CONSERVATIVE_LINEAR_REFINE");
        }
    }
    
    // Register rhs and sol data used by the linear solvers.    
    d_sol_var = new pdat::CellVariable<NDIM,double>(
        "AdvectionDiffusionHierarchyIntegrator::sol",1);
    d_hyp_level_integrator->
        registerVariable(d_sol_var, cell_ghosts,
                         algs::HyperbolicLevelIntegrator<NDIM>::NO_FILL,
                         d_hierarchy->getGridGeometry(),
                         "NO_COARSEN", "NO_REFINE");
    
    d_sol_idx = hier::VariableDatabase<NDIM>::getDatabase()->
        mapVariableAndContextToIndex(d_sol_var, getCurrentContext());

    d_rhs_var = new pdat::CellVariable<NDIM,double>(
        "AdvectionDiffusionHierarchyIntegrator::rhs",1);
    d_hyp_level_integrator->
        registerVariable(d_rhs_var, cell_ghosts,
                         algs::HyperbolicLevelIntegrator<NDIM>::NO_FILL,
                         d_hierarchy->getGridGeometry(),
                         "NO_COARSEN", "NO_REFINE");
                         
    d_rhs_idx = hier::VariableDatabase<NDIM>::getDatabase()->
        mapVariableAndContextToIndex(d_rhs_var, getCurrentContext());

    d_tmp_var = new pdat::CellVariable<NDIM,double>(
        "AdvectionDiffusionHierarchyIntegrator::tmp",1);
    d_hyp_level_integrator->
        registerVariable(d_tmp_var, cell_ghosts,
                         algs::HyperbolicLevelIntegrator<NDIM>::NO_FILL,
                         d_hierarchy->getGridGeometry(),
                         "NO_COARSEN", "NO_REFINE");
    
    d_tmp_idx = hier::VariableDatabase<NDIM>::getDatabase()->
        mapVariableAndContextToIndex(d_tmp_var, getCurrentContext());
    
    // Initialize the algs::HyperbolicLevelIntegrator<NDIM>.
    // NOTE: This must be done after all variables are registered.
    d_hyp_level_integrator->initializeLevelIntegrator(d_gridding_alg);
    
    // Create several communications schedules used in setting up the
    // rhs terms and in synchronizing solution data.
    hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    tbox::Pointer<xfer::RefineOperator<NDIM> > refine_operator;
    tbox::Pointer<xfer::CoarsenOperator<NDIM> > coarsen_operator;
    
    d_ralgs["sol->sol::::CONSTANT_REFINE"] = new xfer::RefineAlgorithm<NDIM>();
    
    refine_operator = grid_geom->lookupRefineOperator(
        d_sol_var, "CONSTANT_REFINE");
    d_ralgs["sol->sol::::CONSTANT_REFINE"]->
        registerRefine(d_sol_idx, // destination
                       d_sol_idx, // source
                       d_tmp_idx, // temporary work space
                       refine_operator);

    d_calgs["SYNCH_NEW_STATE_DATA"] = new xfer::CoarsenAlgorithm<NDIM>();
    
    typedef vector<tbox::Pointer<pdat::CellVariable<NDIM,double> > > CellVariableVector;
    for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
    {
        tbox::Pointer<pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
        coarsen_operator = grid_geom->lookupCoarsenOperator(
            Q_var, "CONSERVATIVE_COARSEN");
        const int Q_new_idx = var_db->
            mapVariableAndContextToIndex(Q_var, getNewContext());
        
        d_calgs["SYNCH_NEW_STATE_DATA"]->
            registerCoarsen(Q_new_idx, // destination
                            Q_new_idx, // source
                            coarsen_operator);
    }

    // Setup the Hierarchy math operations object.
    if (d_hier_math_ops.isNull())
    {
        d_hier_math_ops = new HierarchyMathOps(
            d_object_name+"::HierarchyMathOps", d_hierarchy);
        d_is_managing_hier_math_ops = true;
    }
    
    // Setup the FAC preconditioners.
    set<double> initialized_mus;
    
    for (vector<double>::const_iterator it = d_Q_mus.begin();
         it != d_Q_mus.end(); ++it)
    {
        const double mu = *it;
        
        if (initialized_mus.count(mu) == 0)
        {
            ostringstream stream;
            stream << mu;
            const string& mu_name = stream.str();

            // Helmholtz problem #1 data
            d_helmholtz1_specs[mu] = new solv::PoissonSpecifications(
                d_object_name+"::Helmholtz Specs 1::"+mu_name);
            d_helmholtz1_bc_coefs[mu] = new solv::SimpleCellRobinBcCoefs<NDIM>();

            d_helmholtz1_fac_ops[mu] = new CCPoissonFACOperator(
                d_object_name+"::FAC Ops 1::"+mu_name, "adv_diff_fac_", d_fac_ops_db);
            d_helmholtz1_fac_ops[mu]->setPoissonSpecifications(
                *d_helmholtz1_specs[mu]);
            d_helmholtz1_fac_ops[mu]->setPhysicalBcCoefObject(
                d_helmholtz1_bc_coefs[mu]);
                
            d_helmholtz1_fac_pcs[mu] = new solv::FACPreconditioner<NDIM>(
                d_object_name+"::FAC Preconditioner 1::"+mu_name,
                *d_helmholtz1_fac_ops[mu], d_fac_pcs_db);
            d_helmholtz1_fac_ops[mu]->setPreconditioner(
                d_helmholtz1_fac_pcs[mu]);
            
            // Helmholtz problem #2 data
            d_helmholtz2_specs[mu] = new solv::PoissonSpecifications(
                d_object_name+"::Helmholtz Specs 2::"+mu_name);
            d_helmholtz2_bc_coefs[mu] = new solv::SimpleCellRobinBcCoefs<NDIM>();
            
            d_helmholtz2_fac_ops[mu] = new CCPoissonFACOperator(
                d_object_name+"::FAC Ops 2::"+mu_name,  "adv_diff_fac_", d_fac_ops_db);
            d_helmholtz2_fac_ops[mu]->setPoissonSpecifications(
                *d_helmholtz2_specs[mu]);
            d_helmholtz2_fac_ops[mu]->setPhysicalBcCoefObject(
                d_helmholtz2_bc_coefs[mu]);
                
            d_helmholtz2_fac_pcs[mu] = new solv::FACPreconditioner<NDIM>(
                d_object_name+"::FAC Preconditioner 2::"+mu_name,
                *d_helmholtz2_fac_ops[mu], d_fac_pcs_db);
            d_helmholtz2_fac_ops[mu]->setPreconditioner(
                d_helmholtz2_fac_pcs[mu]);
            
            // Helmholtz problem #3 data
            d_helmholtz3_specs[mu] = new solv::PoissonSpecifications(
                d_object_name+"::Helmholtz Specs 3::"+mu_name);
            d_helmholtz3_bc_coefs[mu] = new solv::SimpleCellRobinBcCoefs<NDIM>();

            if (d_maintain_helmholtz3_solvers)
            {
                d_helmholtz3_fac_ops[mu] = new CCPoissonFACOperator(
                    d_object_name+"::FAC Ops 3::"+mu_name,  "adv_diff_fac_", d_fac_ops_db);
                d_helmholtz3_fac_ops[mu]->setPoissonSpecifications(
                    *d_helmholtz3_specs[mu]);
                d_helmholtz3_fac_ops[mu]->setPhysicalBcCoefObject(
                    d_helmholtz3_bc_coefs[mu]);
                    
                d_helmholtz3_fac_pcs[mu] = new solv::FACPreconditioner<NDIM>(
                    d_object_name+"::FAC Preconditioner 3::"+mu_name,
                    *d_helmholtz3_fac_ops[mu], d_fac_pcs_db);
                d_helmholtz3_fac_ops[mu]->setPreconditioner(
                    d_helmholtz3_fac_pcs[mu]);
            }
            
            // Helmholtz problem #4 data
            d_helmholtz4_specs[mu] = new solv::PoissonSpecifications(
                d_object_name+"::Helmholtz Specs 4::"+mu_name);
            d_helmholtz4_bc_coefs[mu] = new solv::SimpleCellRobinBcCoefs<NDIM>();
            
            if (d_maintain_helmholtz4_solvers)
            {
                d_helmholtz4_fac_ops[mu] = new CCPoissonFACOperator(
                    d_object_name+"::FAC Ops 4::"+mu_name,  "adv_diff_fac_", d_fac_ops_db);
                d_helmholtz4_fac_ops[mu]->setPoissonSpecifications(
                    *d_helmholtz4_specs[mu]);
                d_helmholtz4_fac_ops[mu]->setPhysicalBcCoefObject(
                    d_helmholtz4_bc_coefs[mu]);
                    
                d_helmholtz4_fac_pcs[mu] = new solv::FACPreconditioner<NDIM>(
                    d_object_name+"::FAC Preconditioner 4::"+mu_name,
                    *d_helmholtz4_fac_ops[mu], d_fac_pcs_db);
                d_helmholtz4_fac_ops[mu]->setPreconditioner(
                    d_helmholtz4_fac_pcs[mu]);
            }
            
            // Initialize the Helmholtz solvers.
            if (d_solver_package == "PETSc")
            {
                d_helmholtz1_solvers[mu] = new PETScKrylovLinearSolver(
                    d_object_name+"::PETSc Krylov solver 1::"+mu_name, "adv_diff_");
                d_helmholtz2_solvers[mu] = new PETScKrylovLinearSolver(
                    d_object_name+"::PETSc Krylov solver 2::"+mu_name, "adv_diff_");

                if (d_maintain_helmholtz3_solvers)
                {
                    d_helmholtz3_solvers[mu] = new PETScKrylovLinearSolver(
                        d_object_name+"::PETSc Krylov solver 3::"+mu_name, "adv_diff_");
                }
                if (d_maintain_helmholtz4_solvers)
                {
                    d_helmholtz4_solvers[mu] = new PETScKrylovLinearSolver(
                        d_object_name+"::PETSc Krylov solver 4::"+mu_name, "adv_diff_");
                }
                
                d_using_ksp_method = true;
            }
            else if (d_solver_package == "SAMRAI")
            {
                d_helmholtz1_solvers[mu] = new FACPreconditionerALSWrapper(d_helmholtz1_fac_pcs[mu]);
                d_helmholtz2_solvers[mu] = new FACPreconditionerALSWrapper(d_helmholtz2_fac_pcs[mu]);

                if (d_maintain_helmholtz3_solvers)
                {
                    d_helmholtz3_solvers[mu] = new FACPreconditionerALSWrapper(d_helmholtz3_fac_pcs[mu]);
                }
                if (d_maintain_helmholtz4_solvers)
                {
                    d_helmholtz4_solvers[mu] = new FACPreconditionerALSWrapper(d_helmholtz4_fac_pcs[mu]);
                }
                
                d_using_ksp_method = false;
            }
            else if (d_solver_package == "hypre")
            {
                d_helmholtz1_solvers[mu] = new CCPoissonAMRHypreSolver(
                    d_object_name+"::hypre AMR solver 1::"+mu_name, d_hypre_solver_db);
                d_helmholtz2_solvers[mu] = new CCPoissonAMRHypreSolver(
                    d_object_name+"::hypre AMR solver 2::"+mu_name, d_hypre_solver_db);

                if (d_maintain_helmholtz3_solvers)
                {
                    d_helmholtz3_solvers[mu] = new CCPoissonAMRHypreSolver(
                        d_object_name+"::hypre AMR solver 3::"+mu_name, d_hypre_solver_db);
                }
                if (d_maintain_helmholtz4_solvers)
                {
                    d_helmholtz4_solvers[mu] = new CCPoissonAMRHypreSolver(
                        d_object_name+"::hypre AMR solver 4::"+mu_name, d_hypre_solver_db);
                }

                d_using_ksp_method = false;
            }
            else
            {
                TBOX_ERROR(d_object_name << "::AdvectionDiffusionHierarchyIntegrator():\n" <<
                           "  unknown linear solver package " << d_solver_package << "\n");
            }
            
            d_helmholtz1_solvers[mu]->setMaxIterations(d_max_iterations);
            d_helmholtz1_solvers[mu]->setAbsoluteTolerance(d_abs_residual_tol);
            d_helmholtz1_solvers[mu]->setRelativeTolerance(d_rel_residual_tol);

            d_helmholtz2_solvers[mu]->setMaxIterations(d_max_iterations);
            d_helmholtz2_solvers[mu]->setAbsoluteTolerance(d_abs_residual_tol);
            d_helmholtz2_solvers[mu]->setRelativeTolerance(d_rel_residual_tol);

            if (d_maintain_helmholtz3_solvers)
            {
                d_helmholtz3_solvers[mu]->setMaxIterations(d_max_iterations);
                d_helmholtz3_solvers[mu]->setAbsoluteTolerance(d_abs_residual_tol);
                d_helmholtz3_solvers[mu]->setRelativeTolerance(d_rel_residual_tol);
            }
            if (d_maintain_helmholtz4_solvers)
            {
                d_helmholtz4_solvers[mu]->setMaxIterations(d_max_iterations);
                d_helmholtz4_solvers[mu]->setAbsoluteTolerance(d_abs_residual_tol);
                d_helmholtz4_solvers[mu]->setRelativeTolerance(d_rel_residual_tol);
            }
            
            // Initilize the Helmholtz operators.
            if (d_using_ksp_method)
            {    
                d_helmholtz1_ops[mu] = new CCLaplaceOperator(d_object_name+"::Helmholtz Operator 1::"+mu_name);
                d_helmholtz2_ops[mu] = new CCLaplaceOperator(d_object_name+"::Helmholtz Operator 2::"+mu_name);
                
                tbox::Pointer<AbstractKrylovLinearSolver> krylov1_solver = d_helmholtz1_solvers[mu];
                tbox::Pointer<AbstractKrylovLinearSolver> krylov2_solver = d_helmholtz2_solvers[mu];
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(!krylov1_solver.isNull());
                assert(!krylov2_solver.isNull());
#endif
                krylov1_solver->setInitialGuessNonzero(true);
                krylov2_solver->setInitialGuessNonzero(true);
                
                krylov1_solver->setOperator(d_helmholtz1_ops[mu]);
                krylov2_solver->setOperator(d_helmholtz2_ops[mu]);
                
                krylov1_solver->setPreconditioner(
                    new FACPreconditionerALSWrapper(d_helmholtz1_fac_pcs[mu]));
                krylov2_solver->setPreconditioner(
                    new FACPreconditionerALSWrapper(d_helmholtz2_fac_pcs[mu]));

                if (d_maintain_helmholtz3_solvers)
                {
                    d_helmholtz3_ops[mu] = new CCLaplaceOperator(d_object_name+"::Helmholtz Operator 3::"+mu_name);
                    tbox::Pointer<AbstractKrylovLinearSolver> krylov3_solver = d_helmholtz3_solvers[mu];
#ifdef DEBUG_CHECK_ASSERTIONS
                    assert(!krylov3_solver.isNull());
#endif                
                    krylov3_solver->setInitialGuessNonzero(true);
                    krylov3_solver->setOperator(d_helmholtz3_ops[mu]);
                    krylov3_solver->setPreconditioner(
                        new FACPreconditionerALSWrapper(d_helmholtz3_fac_pcs[mu]));
                }
                if (d_maintain_helmholtz4_solvers)
                {
                    d_helmholtz4_ops[mu] = new CCLaplaceOperator(d_object_name+"::Helmholtz Operator 4::"+mu_name);
                    tbox::Pointer<AbstractKrylovLinearSolver> krylov4_solver = d_helmholtz4_solvers[mu];
#ifdef DEBUG_CHECK_ASSERTIONS
                    assert(!krylov4_solver.isNull());
#endif                
                    krylov4_solver->setInitialGuessNonzero(true);
                    krylov4_solver->setOperator(d_helmholtz4_ops[mu]);
                    krylov4_solver->setPreconditioner(
                        new FACPreconditionerALSWrapper(d_helmholtz4_fac_pcs[mu]));
                }
            }
            
            d_helmholtz_solvers_need_init[mu] = true;
            
            initialized_mus.insert(mu);
        }
    }
    
    // Set the current integration time.
    if (!tbox::RestartManager::getManager()->isFromRestart())
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }
    
    // Indicate that the integrator has been initialized.
    d_is_initialized = true;
    
    t_initialize_hierarchy_integrator->stop();
    return;
}// initializeHierarchyIntegrator

double AdvectionDiffusionHierarchyIntegrator::initializeHierarchy()
{
    t_initialize_hierarchy->start();
    
    if (!d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchy()\n" <<
                   "  must initialize the integrator prior to call to initializeHierarchy()." << endl);
    }
    
    // Initialize the patch hierarchy.
    const bool initial_time = !tbox::RestartManager::getManager()->isFromRestart();
    
    if (!initial_time)
    {
        d_hierarchy->getFromRestart(d_gridding_alg->getMaxLevels());
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        d_gridding_alg->getTagAndInitializeStrategy()->
            resetHierarchyConfiguration(d_hierarchy, coarsest_ln, finest_ln);
    }
    else
    {
        d_gridding_alg->makeCoarsestLevel(d_hierarchy,d_start_time);

        int level_number = 0;
        bool done = false;
        while (!done && (d_gridding_alg->levelCanBeRefined(level_number)))
        {
            d_gridding_alg->
                makeFinerLevel(d_hierarchy,
                               d_integrator_time, initial_time,
                               d_tag_buffer[level_number]);
            
            done = !d_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }
        
        // After data on each level is initialized at simulation start
        // time, coarser levels are synchronized with finer levels
        // that didn't exist when the coarser level initial data was
        // set.
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        
        if (finest_ln > 0)
        {
            synchronizeNewLevels(d_hierarchy, coarsest_ln, finest_ln,
                                 d_start_time, initial_time);
        }
    }

    // The timestep is given by the minimum allowable timestep over
    // all levels in the patch hierarchy.
    double dt_next = numeric_limits<double>::max();
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt_next = tbox::Utilities::
            dmin(dt_next, d_hyp_level_integrator->
                 getLevelDt(level, d_integrator_time, initial_time));
    }
    if (d_integrator_time+dt_next > d_end_time)
    {
        dt_next = d_end_time - d_integrator_time;
    }

    t_initialize_hierarchy->stop();
    return(dt_next);
}// initializeHierarchy

double AdvectionDiffusionHierarchyIntegrator::advanceHierarchy(
    const double dt,
    const bool rebalance_coarsest)
{
    t_advance_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_end_time >= d_integrator_time+dt);
#endif

    // First: rebalance the coarsest level (when requested).
    if (rebalance_coarsest)
    {
        rebalanceCoarsestLevel();
    }
    
    // Second: integrate the data and synchronize the hierarchy.
    const double dt_next =
        integrateHierarchy(d_integrator_time, d_integrator_time+dt);

    synchronizeHierarchy();
    
    // Third: reset all time dependent data.
    resetTimeDependentData(d_integrator_time+dt);

    // Fourth: regrid (when appropriate).
    const bool do_regrid = (d_regrid_interval == 0
                            ? false
                            : (d_integrator_step % d_regrid_interval == 0));
    if (do_regrid)
    {
        regridHierarchy();
    }
    
    t_advance_hierarchy->stop();
    return(dt_next);
}// advanceHierarchy

bool AdvectionDiffusionHierarchyIntegrator::atRegridPoint() const
{
    const int level_number = 0;
    
    return((d_integrator_step>0) 
           && d_gridding_alg->levelCanBeRefined(level_number)
           && (d_regrid_interval == 0
               ? false
               : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double AdvectionDiffusionHierarchyIntegrator::getIntegratorTime() const
{
    return(d_integrator_time);
}// getIntegratorTime
 
double AdvectionDiffusionHierarchyIntegrator::getStartTime() const
{
    return(d_start_time);
}// getStartTime
 
double AdvectionDiffusionHierarchyIntegrator::getEndTime() const
{
    return(d_end_time);
}// getEndTime
 
int AdvectionDiffusionHierarchyIntegrator::getIntegratorStep() const
{
    return(d_integrator_step);
}// getIntegratorStep
 
int AdvectionDiffusionHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return(d_max_integrator_steps);
}// getMaxIntegratorSteps
 
bool AdvectionDiffusionHierarchyIntegrator::stepsRemaining() const
{
    return(d_integrator_step < d_max_integrator_steps);
}// stepsRemaining

const tbox::Pointer<hier::PatchHierarchy<NDIM> >
AdvectionDiffusionHierarchyIntegrator::getPatchHierarchy() const
{
    return(d_hierarchy);
}// getPatchHierarchy
 
tbox::Pointer<mesh::GriddingAlgorithm<NDIM> >
AdvectionDiffusionHierarchyIntegrator::getGriddingAlgorithm() const
{
    return(d_gridding_alg);
}// getGriddingAlgorithm

tbox::Pointer<algs::HyperbolicLevelIntegrator<NDIM> >
AdvectionDiffusionHierarchyIntegrator::getHyperbolicLevelIntegrator() const
{
    return(d_hyp_level_integrator);
}// getHyperbolicLevelIntegrator

tbox::Pointer<AdvDiffHypPatchStrategy>
AdvectionDiffusionHierarchyIntegrator::getHyperbolicPatchStrategy() const
{
    return(d_hyp_patch_strategy);
}// getHyperbolicPatchStrategy

///
///  The following routines:
///
///      rebalanceCoarsestLevel(),
///      regridHierarchy(),
///      integrateHierarchy(),
///      synchronizeHierarchy(),
///      synchronizeNewLevels(),
///      resetTimeDependentData(),
///      resetDataToPreadvanceState()
///
///  allow the AdvectionDiffusionHierarchyIntegrator to provide data
///  management for a time integrator which making use of this class.
///

void AdvectionDiffusionHierarchyIntegrator::rebalanceCoarsestLevel()
{
    t_rebalance_coarsest_level->start();
    
    d_gridding_alg->makeCoarsestLevel(d_hierarchy,d_integrator_time);

    t_rebalance_coarsest_level->stop();
    return;
}// rebalanceCoarsestLevel

void AdvectionDiffusionHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();
    
    const int coarsest_ln = 0;
    d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln,
                                         d_integrator_time, d_tag_buffer);

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy
    
double AdvectionDiffusionHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time)
{
    t_integrate_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(current_time <= new_time);
    assert(d_end_time > d_integrator_time);
    assert(tbox::Utilities::deq(d_integrator_time,current_time));
#endif

    const double dt = new_time - current_time;
    
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    
    // First: Predict the advective terms and synchronize them across
    // all levels of the patch hierarchy.
    
    hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
    
    typedef vector<tbox::Pointer<pdat::CellVariable<NDIM,double> > > CellVariableVector;
    for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
    {
        tbox::Pointer<pdat::CellVariable<NDIM,double> > Q_var   = d_Q_vars[l];
        tbox::Pointer<pdat::CellVariable<NDIM,double> > F_var   = d_F_vars[l];
        tbox::Pointer<pdat::CellVariable<NDIM,double> > Psi_var = d_Psi_vars[l];
        const double mu = d_Q_mus[l];
        tbox::Pointer<SetDataStrategy> F_set = d_F_sets[l];
        
        // Setup the right hand side for the advective flux
        // prediction.
        solv::PoissonSpecifications mu_spec("mu_spec");
        mu_spec.setCConstant(0.0);
        mu_spec.setDConstant(mu);

        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        
        const int Q_current_idx = var_db->
            mapVariableAndContextToIndex(Q_var,  getCurrentContext());
        const int F_current_idx = (F_var.isNull())
            ? -1
            : var_db->mapVariableAndContextToIndex(
                F_var,getCurrentContext());
        const int Psi_current_idx = var_db->
            mapVariableAndContextToIndex(Psi_var,getCurrentContext());
        
        tbox::Pointer<pdat::CellDataFactory<NDIM,double> > Q_factory =
            Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();
        
        for (int depth = 0; depth < Q_depth; ++depth)
        {
            // Setup F and Psi.
            if (!F_set.isNull() && F_set->isTimeDependent())
            {
                F_set->setDataOnPatchHierarchy(
                    F_var, getCurrentContext(),
                    d_hierarchy, d_integrator_time);
            }
            
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                
                for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                    
                    tbox::Pointer<pdat::CellData<NDIM,double> > sol_data = 
                        patch->getPatchData(d_sol_idx);
                    tbox::Pointer<pdat::CellData<NDIM,double> > Q_current_data = 
                        patch->getPatchData(Q_current_idx);
                    
                    sol_data->copyDepth(0,               // dst_depth
                                        *Q_current_data, // src
                                        depth);          // src_depth
                }
            }
            
            d_hier_math_ops->
                laplace(Psi_current_idx, Psi_var  ,  // dst
                        mu_spec,                     // Poisson spec
                        d_sol_idx      , d_sol_var,  // src1
                        d_rscheds["sol->sol::::CONSTANT_REFINE"],
                        d_integrator_time,           // src1_bdry_fill_time
                        1.0,                         // beta
                        F_current_idx, F_var,        // src2
                        depth, 0, depth);
        }
    }
    
    double dt_next = numeric_limits<double>::max();
    
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        
        const bool first_step = true;
        const bool last_step = false;
        
        double level_dt_next = d_hyp_level_integrator->advanceLevel(
            level, d_hierarchy, current_time, new_time,
            first_step, last_step);
        
        dt_next = tbox::Utilities::dmin(dt_next,level_dt_next);
    }
    
    dt_next = tbox::Utilities::dmin(dt_next,d_grow_dt*dt);
    
    if (new_time+dt_next >= d_end_time)
    {
        dt_next = d_end_time - new_time;
    }

    if (finest_ln > 0)
    {
        d_hyp_level_integrator->
            standardLevelSynchronization(d_hierarchy,
                                         coarsest_ln, finest_ln,
                                         new_time, current_time);
    }
    
    // Second: Solve the appropriate Helmholtz equation(s)
    // corresponding to the diffusive timestepping scheme.  Currently,
    // we use the TGA scheme (see below).
    
    // Indicate that all solvers need to be reinitialized if the
    // current timestep is sufficiently different from the previous
    // one.
    if (!tbox::Utilities::deq(dt,d_old_dt))
    {
        for (map<double,bool>::iterator it = d_helmholtz_solvers_need_init.begin();
             it != d_helmholtz_solvers_need_init.end(); ++it)
        {
            (*it).second = true;
        }
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
    }

    typedef vector<tbox::Pointer<pdat::CellVariable<NDIM,double> > > CellVariableVector;
    for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
    {
        tbox::Pointer<pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
        tbox::Pointer<pdat::CellVariable<NDIM,double> > F_var = d_F_vars[l];
        const double mu = d_Q_mus[l];
        tbox::Pointer<SetDataStrategy> F_set = d_F_sets[l];
        
        tbox::Pointer<solv::PoissonSpecifications> helmholtz1_spec      = d_helmholtz1_specs  [mu];
        tbox::Pointer<CCPoissonFACOperator> helmholtz1_fac_op           = d_helmholtz1_fac_ops[mu];
        tbox::Pointer<solv::FACPreconditioner<NDIM> > helmholtz1_fac_pc = d_helmholtz1_fac_pcs[mu];

        tbox::Pointer<AbstractLinearSolver> helmholtz1_solver = d_helmholtz1_solvers[mu];
        tbox::Pointer<CCLaplaceOperator>    helmholtz1_op     = d_helmholtz1_ops    [mu];
        
        tbox::Pointer<solv::PoissonSpecifications> helmholtz2_spec      = d_helmholtz2_specs  [mu];
        tbox::Pointer<CCPoissonFACOperator> helmholtz2_fac_op           = d_helmholtz2_fac_ops[mu];
        tbox::Pointer<solv::FACPreconditioner<NDIM> > helmholtz2_fac_pc = d_helmholtz2_fac_pcs[mu];
        
        tbox::Pointer<AbstractLinearSolver> helmholtz2_solver = d_helmholtz2_solvers[mu];
        tbox::Pointer<CCLaplaceOperator>    helmholtz2_op     = d_helmholtz2_ops    [mu];
        
        tbox::Pointer<solv::PoissonSpecifications> helmholtz3_spec      = d_helmholtz3_specs  [mu];
        tbox::Pointer<CCPoissonFACOperator> helmholtz3_fac_op           = d_helmholtz3_fac_ops[mu];
        tbox::Pointer<solv::FACPreconditioner<NDIM> > helmholtz3_fac_pc = d_helmholtz3_fac_pcs[mu];
        
        tbox::Pointer<AbstractLinearSolver> helmholtz3_solver = d_helmholtz3_solvers[mu];
        tbox::Pointer<CCLaplaceOperator>    helmholtz3_op     = d_helmholtz3_ops    [mu];
        
        tbox::Pointer<solv::PoissonSpecifications> helmholtz4_spec      = d_helmholtz4_specs  [mu];
        tbox::Pointer<CCPoissonFACOperator> helmholtz4_fac_op           = d_helmholtz4_fac_ops[mu];
        tbox::Pointer<solv::FACPreconditioner<NDIM> > helmholtz4_fac_pc = d_helmholtz4_fac_pcs[mu];
        
        tbox::Pointer<AbstractLinearSolver> helmholtz4_solver = d_helmholtz4_solvers[mu];
        tbox::Pointer<CCLaplaceOperator>    helmholtz4_op     = d_helmholtz4_ops    [mu];
        
        // Coefficients corresponding to the TGA discretization.
        //
        // Ref: McCorquodale, Colella, Johansen.  "A Cartesian grid
        // embedded boundary method for the heat equation on irregular
        // domains." JCP 173, pp. 620-635 (2001)
        
        const double a = 2.0 - sqrt(2.0) - numeric_limits<double>::epsilon();
        const double mu1 = 0.5*(a-sqrt(a*a-4.0*a+2.0))*dt;
        const double mu2 = 0.5*(a+sqrt(a*a-4.0*a+2.0))*dt;
        const double mu3 = (1.0-a)*dt;
        const double mu4 = (0.5-a)*dt;
        
        // Initialize the linear solver.
        if (d_helmholtz_solvers_need_init[mu])
        {
            if (d_do_log) tbox::plog << d_object_name << ": "
                                     << "Initializing Helmholtz solvers for mu = " << mu
                                     << ", dt = " << dt << "\n";
            
            helmholtz1_spec->setCConstant(1.0);
            helmholtz1_spec->setDConstant(-mu*mu1);

            helmholtz2_spec->setCConstant(1.0);
            helmholtz2_spec->setDConstant(-mu*mu2);

            helmholtz3_spec->setCConstant(1.0);
            helmholtz3_spec->setDConstant(mu*mu3);

            helmholtz4_spec->setCConstant(1.0);
            helmholtz4_spec->setDConstant(mu*mu4);

            if (d_solver_package == "hypre")
            {
                tbox::Pointer<CCPoissonAMRHypreSolver> hypre_solver;

                hypre_solver = d_helmholtz1_solvers[mu];
                hypre_solver->setPoissonSpecifications(*helmholtz1_spec);
                hypre_solver->setResetLevels(d_coarsest_reset_ln, d_finest_reset_ln);
                
                hypre_solver = d_helmholtz2_solvers[mu];
                hypre_solver->setPoissonSpecifications(*helmholtz2_spec);
                hypre_solver->setResetLevels(d_coarsest_reset_ln, d_finest_reset_ln);
                
                if (d_maintain_helmholtz3_solvers)
                {
                    hypre_solver = d_helmholtz3_solvers[mu];
                    hypre_solver->setPoissonSpecifications(*helmholtz3_spec);
                    hypre_solver->setResetLevels(d_coarsest_reset_ln, d_finest_reset_ln);
                }
                
                if (d_maintain_helmholtz4_solvers)
                {
                    hypre_solver = d_helmholtz4_solvers[mu];
                    hypre_solver->setPoissonSpecifications(*helmholtz4_spec);
                    hypre_solver->setResetLevels(d_coarsest_reset_ln, d_finest_reset_ln);
                }
            }
            else
            {
                helmholtz1_fac_op->setPoissonSpecifications(*helmholtz1_spec);
                helmholtz1_fac_op->setResetLevels(d_coarsest_reset_ln, d_finest_reset_ln);
            
                helmholtz2_fac_op->setPoissonSpecifications(*helmholtz2_spec);
                helmholtz2_fac_op->setResetLevels(d_coarsest_reset_ln, d_finest_reset_ln);
                
                if (d_maintain_helmholtz3_solvers)
                {
                    helmholtz3_fac_op->setPoissonSpecifications(*helmholtz3_spec);
                    helmholtz3_fac_op->setResetLevels(d_coarsest_reset_ln, d_finest_reset_ln);
                }
                
                if (d_maintain_helmholtz4_solvers)
                {
                    helmholtz4_fac_op->setPoissonSpecifications(*helmholtz4_spec);
                    helmholtz4_fac_op->setResetLevels(d_coarsest_reset_ln, d_finest_reset_ln);
                }
            }
            
            if (d_using_ksp_method)
            {
                helmholtz1_op->setPoissonSpecifications(helmholtz1_spec);
                helmholtz1_op->setHierarchyMathOps(d_hier_math_ops);

                helmholtz2_op->setPoissonSpecifications(helmholtz2_spec);
                helmholtz2_op->setHierarchyMathOps(d_hier_math_ops);

                if (d_maintain_helmholtz3_solvers)
                {
                    helmholtz3_op->setPoissonSpecifications(helmholtz3_spec);
                    helmholtz3_op->setHierarchyMathOps(d_hier_math_ops);
                }
                if (d_maintain_helmholtz4_solvers)
                {
                    helmholtz4_op->setPoissonSpecifications(helmholtz4_spec);
                    helmholtz4_op->setHierarchyMathOps(d_hier_math_ops);
                }
            }

            helmholtz1_solver->initializeSolverState(*d_sol_vec,*d_rhs_vec);
            helmholtz2_solver->initializeSolverState(*d_sol_vec,*d_rhs_vec);
            if (d_maintain_helmholtz3_solvers) helmholtz3_solver->initializeSolverState(*d_sol_vec,*d_rhs_vec);
            if (d_maintain_helmholtz4_solvers) helmholtz4_solver->initializeSolverState(*d_sol_vec,*d_rhs_vec);
            
            d_helmholtz_solvers_need_init[mu] = false;
        }
        
        if (d_using_ksp_method)
        {
            helmholtz1_op->setTime(d_integrator_time);
            helmholtz2_op->setTime(d_integrator_time);
            if (d_maintain_helmholtz3_solvers) helmholtz3_op->setTime(d_integrator_time);
            if (d_maintain_helmholtz4_solvers) helmholtz4_op->setTime(d_integrator_time);
        }
        
        // Setup the rhs terms and solve the systems.
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        
        const int Q_current_idx = var_db->
            mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int F_current_idx = (F_var.isNull())
            ? -1
            : var_db->mapVariableAndContextToIndex(
                F_var,getCurrentContext());
        const int Q_new_idx = var_db->
            mapVariableAndContextToIndex(Q_var, getNewContext());
        
        tbox::Pointer<pdat::CellDataFactory<NDIM,double> > Q_factory =
            Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();
        
        if (!F_set.isNull() && F_set->isTimeDependent())
        {
            F_set->setDataOnPatchHierarchy(
                F_var, getCurrentContext(),
                d_hierarchy, d_integrator_time+0.5*dt);
        }
        
        if (!F_var.isNull())
        {
            d_hier_cc_data_ops->add(Q_new_idx, F_current_idx, Q_new_idx);
        }
        
        for (int depth = 0; depth < Q_depth; ++depth)
        {
            // Setup the rhs data part 1:
            // 
            //      rhs := (I+mu4*L)*F(n+1/2)
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                
                for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                    
                    tbox::Pointer<pdat::CellData<NDIM,double> > sol_data = 
                        patch->getPatchData(d_sol_idx);
                    tbox::Pointer<pdat::CellData<NDIM,double> > Q_new_data = 
                        patch->getPatchData(Q_new_idx);
                    
                    sol_data->copyDepth(0,           // dst_depth
                                        *Q_new_data, // src
                                        depth);      // src_depth
                }
            }
            
            d_hier_math_ops->
                laplace(d_rhs_idx, d_rhs_var,  // dst
                        *helmholtz4_spec,      // Poisson spec
                        d_sol_idx, d_sol_var,  // src
                        d_rscheds["sol->sol::::CONSTANT_REFINE"],
                        d_integrator_time);    // src_bdry_fill_time
            
            // Setup the rhs data part 2:
            // 
            //      rhs := (I+mu3*L)*Q(n) + (I+mu4*L)*F(n+1/2)*dt
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                
                for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                    
                    tbox::Pointer<pdat::CellData<NDIM,double> > sol_data = 
                        patch->getPatchData(d_sol_idx);
                    tbox::Pointer<pdat::CellData<NDIM,double> > Q_current_data = 
                        patch->getPatchData(Q_current_idx);
                    
                    sol_data->copyDepth(0,               // dst_depth
                                        *Q_current_data, // src
                                        depth);          // src_depth
                }
            }
            
            d_hier_math_ops->
                laplace(d_rhs_idx, d_rhs_var,   // dst
                        *helmholtz3_spec,       // Poisson spec
                        d_sol_idx, d_sol_var,   // src1
                        d_rscheds["sol->sol::::CONSTANT_REFINE"],
                        d_integrator_time,      // src1_bdry_fill_time
                        dt,                     // gamma
                        d_rhs_idx, d_rhs_var);  // src2

            // Solve the linear systems.
            helmholtz2_solver->solveSystem(*d_sol_vec,*d_rhs_vec);
            d_rhs_vec->copyVector(d_sol_vec);
            helmholtz1_solver->solveSystem(*d_sol_vec,*d_rhs_vec);
            if (d_do_log) tbox::plog << "AdvectionDiffusionHierarchyIntegrator::integrateHierarchy(): solve 1 number of iterations = " << helmholtz2_solver->getNumIterations() << "\n";
            if (d_do_log) tbox::plog << "AdvectionDiffusionHierarchyIntegrator::integrateHierarchy(): solve 1 residual norm        = " << helmholtz2_solver->getResidualNorm()  << "\n";
            if (d_do_log) tbox::plog << "AdvectionDiffusionHierarchyIntegrator::integrateHierarchy(): solve 2 number of iterations = " << helmholtz1_solver->getNumIterations() << "\n";
            if (d_do_log) tbox::plog << "AdvectionDiffusionHierarchyIntegrator::integrateHierarchy(): solve 2 residual norm        = " << helmholtz1_solver->getResidualNorm()  << "\n";
            
            // Put the sol data into Q_new.
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                
                for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                    
                    tbox::Pointer<pdat::CellData<NDIM,double> > sol_data = 
                        patch->getPatchData(d_sol_idx);
                    tbox::Pointer<pdat::CellData<NDIM,double> > Q_new_data = 
                        patch->getPatchData(Q_new_idx);
                    
                    Q_new_data->copyDepth(depth,     // dst_depth,
                                          *sol_data, // src
                                          0);        // src_depth
                }
            }
        }
    }

    t_integrate_hierarchy->stop();
    return(dt_next);
}// integrateHierarchy

void AdvectionDiffusionHierarchyIntegrator::synchronizeHierarchy()
{
    t_synchronize_hierarchy->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_NEW_STATE_DATA"][ln]->coarsenData();
    }
    
    t_synchronize_hierarchy->stop();
    return;
}// synchronizeHierarchy

void AdvectionDiffusionHierarchyIntegrator::synchronizeNewLevels(
    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const bool initial_time)
{
    t_synchronize_new_levels->start();
    
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((coarsest_level >= 0)
           && (coarsest_level < finest_level)
           && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        assert(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    // We use the HyperbolicLevelIntegrator to handle as much data
    // management as possible.
    d_hyp_level_integrator->
        synchronizeNewLevels(hierarchy, coarsest_level, finest_level,
                             sync_time, initial_time);

    t_synchronize_new_levels->stop();
    return;
}// synchronizeNewLevels

void AdvectionDiffusionHierarchyIntegrator::resetTimeDependentData(
    const double new_time)
{
    t_reset_time_dependent_data->start();

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // Reset the time dependent data.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hyp_level_integrator->
            resetTimeDependentData(d_hierarchy->getPatchLevel(ln),
                                   d_integrator_time,
                                   d_gridding_alg->levelCanBeRefined(ln));
    }

    t_reset_time_dependent_data->stop();
    return;
}// resetTimeDependentData

void AdvectionDiffusionHierarchyIntegrator::resetDataToPreadvanceState()
{
    t_reset_data_to_preadvance_state->start();

    // We use the HyperbolicLevelIntegrator to handle as much data
    // management as possible.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hyp_level_integrator->
            resetDataToPreadvanceState(d_hierarchy->getPatchLevel(ln));
    }

    t_reset_data_to_preadvance_state->stop();
    return;
}// resetDataToPreadvanceState

///
///  The following routines:
///
///      initializeLevelData(),
///      resetHierarchyConfiguration(),
///      applyGradientDetector()
///
///  are concrete implementations of functions declared in the
///  mesh::StandardTagAndInitStrategy<NDIM> abstract base class.
///

void AdvectionDiffusionHierarchyIntegrator::initializeLevelData(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > base_old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const tbox::Pointer<hier::PatchLevel<NDIM> > old_level = base_old_level;
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((level_number >= 0) 
           && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!(old_level.isNull())) {
        assert(level_number == old_level->getLevelNumber());
    }
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    // We use the HyperbolicLevelIntegrator to handle as much data
    // management as possible.
    d_hyp_level_integrator->
        initializeLevelData(hierarchy, level_number, init_data_time,
                            can_be_refined, initial_time, old_level,
                            allocate_data);

    // If a ConvergenceMonitor object is registered with the
    // integrator, initialize any level data associated with it.
    if (!d_convergence_monitor.isNull())
    {
        d_convergence_monitor->
            initializeLevelData(hierarchy, level_number, init_data_time,
                                can_be_refined, initial_time, old_level,
                                allocate_data);
    }
    
    // Set the initial values of any forcing terms.
    if (initial_time)
    {
        hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
        tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
        typedef vector<tbox::Pointer<pdat::CellVariable<NDIM,double> > > CellVariableVector;
        for (CellVariableVector::size_type l = 0; l < d_F_vars.size(); ++l)
        {
            tbox::Pointer<pdat::CellVariable<NDIM,double> > F_var = d_F_vars[l];
        
            if (!F_var.isNull())
            {
                const int F_idx = var_db->mapVariableAndContextToIndex(
                    F_var, getCurrentContext());
                tbox::Pointer<SetDataStrategy> F_set = d_F_sets [l];
                
                if (!F_set.isNull())
                {
                    F_set->setDataOnPatchLevel(F_idx, F_var, level,
                                               init_data_time, initial_time);
                }
                else
                {
                    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                        
                        tbox::Pointer<pdat::CellData<NDIM,double> > F_data =
                            patch->getPatchData(F_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(!F_data.isNull());
#endif
                        F_data->fillAll(0.0);
                    }
                }
            }
        }
    }
    
    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void AdvectionDiffusionHierarchyIntegrator::resetHierarchyConfiguration(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((coarsest_level >= 0)
           && (coarsest_level <= finest_level) 
           && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln0 = 0; ln0 <= finest_level; ++ln0)
    {
        assert(!(hierarchy->getPatchLevel(ln0)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // We use the HyperbolicLevelIntegrator to handle as much data
    // management as possible.
    d_hyp_level_integrator->
        resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    
    // If a ConvergenceMonitor object is registered with the
    // integrator, reset any data associated with it.
    if (!d_convergence_monitor.isNull())
    {
        d_convergence_monitor->resetHierarchyConfiguration(
            hierarchy, coarsest_level, finest_level);
    }
    
    // Reset the Hierarchy data operations for the new hierarchy
    // configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);
    
    // Reset the Hierarchy math operations for the new configuration.
    if (d_is_managing_hier_math_ops)
    {
        d_hier_math_ops->setPatchHierarchy(hierarchy);
        d_hier_math_ops->resetLevels(0, finest_hier_level);
    }
    
    // Get the cell weights data.
    d_wgt_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    
    // Reset the solution and rhs vectors.
    d_sol_vec = new solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::sol_vec", d_hierarchy, 0, finest_hier_level);
    d_sol_vec->addComponent(d_sol_var,d_sol_idx,d_wgt_idx,d_hier_cc_data_ops);
    
    d_rhs_vec = new solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::rhs_vec", d_hierarchy, 0, finest_hier_level);
    d_rhs_vec->addComponent(d_rhs_var,d_rhs_idx,d_wgt_idx,d_hier_cc_data_ops);
    
    // Indicate that all linear solvers must be re-initialized.
    for (map<double,bool>::iterator it =
             d_helmholtz_solvers_need_init.begin();
         it != d_helmholtz_solvers_need_init.end();
         ++it)
    {
        (*it).second = true;
    }

    d_coarsest_reset_ln = coarsest_level;
    d_finest_reset_ln = finest_level;
    
    // If we have added or removed a level, resize the schedule
    // vectors.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it != d_ralgs.end(); ++it)
    {
        d_rscheds[(*it).first].resize(finest_hier_level+1);
    }
    
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it != d_calgs.end(); ++it)
    {
        d_cscheds[(*it).first].resize(finest_hier_level+1);
    }
    
    // (Re)build refine communication schedules.  These are created
    // for all levels in the hierarchy.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it != d_ralgs.end(); ++it)
    {
        for (int ln = coarsest_level; ln <= finest_hier_level; ++ln)
        {
            tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            
            d_rscheds[(*it).first][ln] = (*it).second->
                createSchedule(level, ln-1, hierarchy);
        }
    }
    
    // (Re)build coarsen communication schedules.  These are set only
    // for levels >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it != d_calgs.end(); ++it)
    {
        for (int ln = tbox::Utilities::imax(coarsest_level,1);
             ln <= finest_hier_level; ++ln)
        {
            tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            tbox::Pointer<hier::PatchLevel<NDIM> > coarser_level =
                hierarchy->getPatchLevel(ln-1);
            
            d_cscheds[(*it).first][ln] = (*it).second->
                createSchedule(coarser_level, level);
        }
    }
    
    // Reset the "empty" schedule vectors.
    d_rscheds["NONE"]=vector<tbox::Pointer<xfer::RefineSchedule<NDIM> > >(finest_hier_level+1);
    d_cscheds["NONE"]=vector<tbox::Pointer<xfer::CoarsenSchedule<NDIM> > >(finest_hier_level+1);
    
    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration
    
void AdvectionDiffusionHierarchyIntegrator::applyGradientDetector(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy, 
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    t_apply_gradient_detector->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((level_number >= 0)
           && (level_number <= hierarchy->getFinestLevelNumber()));
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    
    // Due to tag buffers, it is necessary to untag all cells prior to
    // tagging.
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
        tbox::Pointer<pdat::CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
        tags_data->fillAll(0);
    }
    
    // Tag cells for refinement according to the criteria specified by
    // the criteria specified by the algs::HyperbolicLevelIntegrator<NDIM>.
    d_hyp_level_integrator->
        applyGradientDetector(hierarchy, level_number, error_data_time,
                              tag_index, initial_time,
                              uses_richardson_extrapolation_too);
    
    t_apply_gradient_detector->stop();
    return;
}// applyGradientDetector
    
///
///  The following routines:
///
///      getCurrentContext(),
///      getNewContext(),
///      getOldContext(),
///      getScratchContext(),
///      getPlotContext()
///
///  allow access to the various variable contexts maintained by the
///  integrator.
///

///
/// We simply reuse the hier::VariableContext objects defined in the
/// algs::HyperbolicLevelIntegrator<NDIM> object.
///

tbox::Pointer<hier::VariableContext> AdvectionDiffusionHierarchyIntegrator::getCurrentContext() const
{
    return(d_hyp_level_integrator->getCurrentContext());
}// getCurrentContext

tbox::Pointer<hier::VariableContext> AdvectionDiffusionHierarchyIntegrator::getNewContext() const
{
    return(d_hyp_level_integrator->getNewContext());
}// getNewContext

tbox::Pointer<hier::VariableContext> AdvectionDiffusionHierarchyIntegrator::getOldContext() const
{
    return(d_hyp_level_integrator->getOldContext());
}// getOldContext

tbox::Pointer<hier::VariableContext> AdvectionDiffusionHierarchyIntegrator::getScratchContext() const
{
    return(d_hyp_level_integrator->getScratchContext());
}// getScratchContext

tbox::Pointer<hier::VariableContext> AdvectionDiffusionHierarchyIntegrator::getPlotContext() const
{
    return(d_hyp_level_integrator->getPlotContext());
}// getPlotContext

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  tbox::Serializable abstract base class.
///
    
void AdvectionDiffusionHierarchyIntegrator::putToDatabase(
    tbox::Pointer<tbox::Database> db)
{
    t_put_to_database->start();
    
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif
    
    db->putInteger("ADVECTION_DIFFUSION_HIERARCHY_INTEGRATOR_VERSION",
                   ADVECTION_DIFFUSION_HIERARCHY_INTEGRATOR_VERSION);
    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putBool("d_using_default_tag_buffer", d_using_default_tag_buffer);
    db->putIntegerArray("d_tag_buffer", d_tag_buffer);
    db->putDouble("d_integrator_time", d_integrator_time);
    db->putInteger("d_integrator_step", d_integrator_step);

    t_put_to_database->stop();
    return;
}// putToDatabase
    
///
///  The following routines:
///
///      printClassData()
///
///  are provided for your viewing pleasure.
///
    
void AdvectionDiffusionHierarchyIntegrator::printClassData(
    ostream& os) const
{
    os << "\nAdvectionDiffusionHierarchyIntegrator::printClassData..." << endl;
    os << "\nAdvectionDiffusionHierarchyIntegrator: this = "
       << (AdvectionDiffusionHierarchyIntegrator*)this << endl;
    os << "d_object_name = " << d_object_name << endl;
    os << "d_integrator_time = " << d_integrator_time << "\n"
       << "d_start_time = " << d_start_time << "\n"
       << "d_end_time = " << d_end_time << "\n"
       << "d_integrator_step = " << d_integrator_step << "\n"
       << "d_grow_dt = " << d_grow_dt << endl;
    os << "I AM INCOMPLETE!!!!!!!!!" << endl;
    return;
}// printClassData

/////////////////////////////// PRIVATE //////////////////////////////////////

void AdvectionDiffusionHierarchyIntegrator::getFromInput(
    tbox::Pointer<tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif
    // Read in data members from input database.
    d_end_time = db->getDoubleWithDefault("end_time", d_end_time);
    d_grow_dt = db->getDoubleWithDefault("grow_dt", d_grow_dt);
    
    d_max_integrator_steps = db->getIntegerWithDefault(
        "max_integrator_steps", d_max_integrator_steps);
    
    d_regrid_interval = db->getIntegerWithDefault(
        "regrid_interval", d_regrid_interval);
    
    if (db->keyExists("tag_buffer"))
    {
        d_tag_buffer = db->getIntegerArray("tag_buffer");
        d_using_default_tag_buffer = false;
    }
    else
    {
        d_using_default_tag_buffer = true;
        TBOX_WARNING(d_object_name << ":  "
                     << "Key data `tag_buffer' not found in input.  "
                     << "Default values used.  See class header for details.");
    }
    
    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    d_solver_package = db->getStringWithDefault("solver_package", d_solver_package);
    d_max_iterations = db->getIntegerWithDefault("max_iterations", d_max_iterations);
    d_abs_residual_tol = db->getDoubleWithDefault("abs_residual_tol", d_abs_residual_tol);
    d_rel_residual_tol = db->getDoubleWithDefault("rel_residual_tol", d_rel_residual_tol);

    if (!is_from_restart)
    {
        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);
    }

    return;
}// getFromInput

void AdvectionDiffusionHierarchyIntegrator::getFromRestart()
{
    tbox::Pointer<tbox::Database> restart_db =
        tbox::RestartManager::getManager()->getRootDatabase();
    
    tbox::Pointer<tbox::Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }
    
    int ver = db->getInteger("ADVECTION_DIFFUSION_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != ADVECTION_DIFFUSION_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }
    
    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    d_grow_dt = db->getDouble("d_grow_dt");
    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");
    d_regrid_interval = db->getInteger("d_regrid_interval");
    d_using_default_tag_buffer = db->getBool("d_using_default_tag_buffer");
    d_tag_buffer = db->getIntegerArray("d_tag_buffer");
    d_integrator_time = db->getDouble("d_integrator_time");
    d_integrator_step = db->getInteger("d_integrator_step");

    return;
}// getFromRestart
    
/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION

#include "tbox/Pointer.C"

//////////////////////////////////////////////////////////////////////
///
/// These declarations are required to use the
/// AdvectionDiffusionHierarchyIntegrator class.
///
//////////////////////////////////////////////////////////////////////

template class tbox::Pointer<AdvectionDiffusionHierarchyIntegrator>;

#endif

//////////////////////////////////////////////////////////////////////////////
