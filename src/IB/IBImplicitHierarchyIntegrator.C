// Filename: IBImplicitHierarchyIntegrator.C
// Last modified: <12.Feb.2008 21:15:18 griffith@box221.cims.nyu.edu>
// Created on 30 Mar 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBImplicitHierarchyIntegrator.h"

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
#include <ibamr/LEInteractor.h>
#include <ibamr/LNodeIndexData.h>

// STOOLS INCLUDES
#include <stools/CartExtrapPhysBdryOp.h>
#include <stools/PETSC_SAMRAI_ERROR.h>
#include <stools/STOOLS_Utilities.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CellIndex.h>
#include <CoarsenOperator.h>
#include <HierarchyDataOpsManager.h>
#include <Index.h>
#include <Patch.h>
#include <PatchCellDataOpsReal.h>
#include <RefineOperator.h>
#include <VariableDatabase.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <limits>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy_integrator;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_advance_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_rebalance_coarsest_level;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_regrid_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_integrate_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_new_levels;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_time_dependent_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_data_to_preadvance_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_hierarchy_configuration;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_apply_gradient_detector;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;

// Whether to account for periodic shifts in spreading and
// interpolating.
static const bool ENFORCE_PERIODIC_BCS = true;

// Version of IBImplicitHierarchyIntegrator restart file data.
static const int IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION = 1;

// Global variables are evil.
double dt;
std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > X_data;
std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > F_data;
std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > U_data;
std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > X_new_data;
std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > F_new_data;
std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > U_new_data;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitHierarchyIntegrator::IBImplicitHierarchyIntegrator(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    SAMRAI::tbox::Pointer<INSHierarchyIntegrator> ins_hier_integrator,
    SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> force_strategy,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_delta_fcn("IB_4"),
      d_ghosts(-1),
      d_hierarchy(hierarchy),
      d_gridding_alg(NULL),
      d_visit_writer(NULL),
      d_silo_writer(NULL),
      d_load_balancer(NULL),
      d_ins_hier_integrator(ins_hier_integrator),
      d_lag_data_manager(NULL),
      d_lag_init(NULL),
      d_eulerian_force_set(NULL),
      d_force_strategy(force_strategy),
      d_force_strategy_needs_init(true),
      d_start_time(0.0),
      d_end_time(std::numeric_limits<double>::max()),
      d_grow_dt(2.0),
      d_max_integrator_steps(std::numeric_limits<int>::max()),
      d_regrid_interval(1),
      d_old_dt(-1.0),
      d_integrator_time(std::numeric_limits<double>::quiet_NaN()),
      d_integrator_step(std::numeric_limits<int>::max()),
      d_dt_max(std::numeric_limits<double>::max()),
      d_dt_max_time_max(std::numeric_limits<double>::max()),
      d_dt_max_time_min(-(d_dt_max_time_max-std::numeric_limits<double>::epsilon())),
      d_is_initialized(false),
      d_do_log(false),
      d_ralgs(),
      d_rstrategies(),
      d_rscheds(),
      d_calgs(),
      d_cstrategies(),
      d_cscheds(),
      d_force_ralg(),
      d_force_rstrategy(),
      d_force_rscheds(),
      d_hier_cc_data_ops(),
      d_V_var(NULL),
      d_W_var(NULL),
      d_F_var(NULL),
      d_Q_var(NULL),
      d_current(NULL),
      d_scratch(NULL),
      d_V_idx(-1),
      d_W_idx(-1),
      d_F_idx(-1),
      d_F_scratch1_idx(-1),
      d_F_scratch2_idx(-1),
      d_Q_idx(-1),
      d_Q_scratch_idx(-1)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT(!ins_hier_integrator.isNull());
    TBOX_ASSERT(!force_strategy.isNull());
#endif

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the input and restart
    // databases.
    const bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Determine the ghost cell width required for cell-centered
    // spreading and interpolating.
    const int stencil_size = LEInteractor::getStencilSize(d_delta_fcn);
    d_ghosts = static_cast<int>(floor(0.5*static_cast<double>(stencil_size))+1);

    // Get the Lagrangian Data Manager.
    d_lag_data_manager = LDataManager::getManager(
        d_object_name+"::LDataManager", d_ghosts, d_registered_for_restart);

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager =
        SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var = new SAMRAI::pdat::CellVariable<NDIM,double>(
        "cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(
        cc_var, hierarchy);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::advanceHierarchy()");
        t_rebalance_coarsest_level = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::rebalanceCoarsestLevel()");
        t_regrid_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::regridHierarchy()");
        t_integrate_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::integrateHierarchy()");
        t_synchronize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBImplicitHierarchyIntegrator::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// IBImplicitHierarchyIntegrator

IBImplicitHierarchyIntegrator::~IBImplicitHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    for (RefinePatchStrategyMap::iterator it = d_rstrategies.begin();
         it != d_rstrategies.end(); ++it)
    {
        delete (*it).second;
    }

    for (CoarsenPatchStrategyMap::iterator it = d_cstrategies.begin();
         it != d_cstrategies.end(); ++it)
    {
        delete (*it).second;
    }

    return;
}// ~IBImplicitHierarchyIntegrator

const std::string&
IBImplicitHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
IBImplicitHierarchyIntegrator::registerLNodeInitStrategy(
    SAMRAI::tbox::Pointer<LNodeInitStrategy> lag_init)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!lag_init.isNull());
#endif
    d_lag_init = lag_init;
    d_lag_data_manager->registerLNodeInitStrategy(d_lag_init);
    return;
}// registerLNodeInitStrategy

void
IBImplicitHierarchyIntegrator::freeLNodeInitStrategy()
{
    d_lag_init.setNull();
    d_lag_data_manager->freeLNodeInitStrategy();
    return;
}// freeLNodeInitStrategy

void
IBImplicitHierarchyIntegrator::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    d_ins_hier_integrator->registerVisItDataWriter(d_visit_writer);
    d_lag_data_manager->registerVisItDataWriter(d_visit_writer);
    return;
}// registerVisItDataWriter

void
IBImplicitHierarchyIntegrator::registerLagSiloDataWriter(
    SAMRAI::tbox::Pointer<LagSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!silo_writer.isNull());
#endif
    d_silo_writer = silo_writer;
    d_lag_data_manager->registerLagSiloDataWriter(d_silo_writer);
    return;
}// registerLagSiloDataWriter

void
IBImplicitHierarchyIntegrator::registerLoadBalancer(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!load_balancer.isNull());
#endif
    d_load_balancer = load_balancer;
    d_lag_data_manager->registerLoadBalancer(d_load_balancer);
    return;
}// registerLoadBalancer

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
///      getLDataManager()
///
///  allow the IBImplicitHierarchyIntegrator to be used as a hierarchy
///  integrator.
///

void
IBImplicitHierarchyIntegrator::initializeHierarchyIntegrator(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Initialize all variables.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    d_current = var_db->getContext(d_object_name+"::CURRENT");
    d_scratch = var_db->getContext(d_object_name+"::SCRATCH");
    const SAMRAI::hier::IntVector<NDIM> ghosts = d_ghosts;
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    d_V_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::V",NDIM);
    d_V_idx = var_db->registerVariableAndContext(d_V_var, d_current, ghosts);

    d_W_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::W",NDIM);
    d_W_idx = var_db->registerVariableAndContext(d_W_var, d_current, ghosts);

    d_F_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::F",NDIM);
    d_F_idx = var_db->registerVariableAndContext(d_F_var, d_current, no_ghosts);
    d_F_scratch1_idx = var_db->registerVariableAndContext(d_F_var, d_scratch, no_ghosts);
    d_F_scratch2_idx = var_db->registerClonedPatchDataIndex(d_F_var, d_F_scratch1_idx);

    // Initialize the objects used to manage Lagragian-Eulerian
    // interaction.
    d_eulerian_force_set = new IBEulerianForceSetter(
        d_object_name+"::IBEulerianForceSetter", -1, d_F_idx, -1);
    d_ins_hier_integrator->registerForceSpecification(d_eulerian_force_set);

    // Initialize the INSHierarchyIntegrator.
    //
    // NOTE: This must be done after all variables created by the
    // IBImplicitHierarchyIntegrator class have been registered.
    d_ins_hier_integrator->initializeHierarchyIntegrator(d_gridding_alg);

    // Create several communications algorithms, used in filling ghost
    // cell data and synchronizing data on the patch hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_operator;

    const int U_current_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getCurrentContext());

    d_ralgs["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"] =
        new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getVelocityVar(),
        "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"]->
        registerRefine(d_V_idx,       // destination
                       U_current_idx, // source
                       d_V_idx,       // temporary work space
                       refine_operator);
    d_rstrategies["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"] =
        new STOOLS::CartExtrapPhysBdryOp(d_V_idx, "LINEAR");

    const int U_new_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getNewContext());

    d_ralgs["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"] =
        new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getVelocityVar(),
        "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"]->
        registerRefine(d_W_idx,   // destination
                       U_new_idx, // source
                       d_W_idx,   // temporary work space
                       refine_operator);
    d_rstrategies["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"] =
        new STOOLS::CartExtrapPhysBdryOp(d_W_idx, "LINEAR");

    // NOTE: When using conservative averaging to coarsen the velocity
    // from finer levels to coarser levels, the appropriate
    // prolongation operator for the force is constant refinement.
    // This choice results in IB spreading and interpolation being
    // adjoint.
    const int F_current_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getForceVar(),
        d_ins_hier_integrator->getCurrentContext());

    d_force_ralg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getForceVar(), "CONSTANT_REFINE");
    d_force_ralg->registerRefine(F_current_idx,     // destination
                                 F_current_idx,     // source
                                 d_F_scratch1_idx,  // temporary work space
                                 refine_operator);
    refine_operator = grid_geom->lookupRefineOperator(
        d_F_var, "CONSTANT_REFINE");
    d_force_ralg->registerRefine(d_F_idx,           // destination
                                 d_F_idx,           // source
                                 d_F_scratch2_idx,  // temporary work space
                                 refine_operator);
    SAMRAI::hier::ComponentSelector F_scratch_idxs;
    F_scratch_idxs.setFlag(d_F_scratch1_idx);
    F_scratch_idxs.setFlag(d_F_scratch2_idx);
    d_force_rstrategy = new STOOLS::CartExtrapPhysBdryOp(
        F_scratch_idxs, "LINEAR");

    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"] =
        new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_V_var, "CONSERVATIVE_COARSEN");
    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"]->
        registerCoarsen(U_current_idx, // destination
                        U_current_idx, // source
                        coarsen_operator);

    d_calgs["U->U::N->N::CONSERVATIVE_COARSEN"] =
        new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_W_var, "CONSERVATIVE_COARSEN");
    d_calgs["U->U::N->N::CONSERVATIVE_COARSEN"]->
        registerCoarsen(U_new_idx, // destination
                        U_new_idx, // source
                        coarsen_operator);

    d_calgs["F->F::CONSERVATIVE_COARSEN"] =
        new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_ins_hier_integrator->getForceVar(), "CONSERVATIVE_COARSEN");
    d_calgs["F->F::CONSERVATIVE_COARSEN"]->
        registerCoarsen(F_current_idx, // destination
                        F_current_idx, // source
                        coarsen_operator);
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_F_var, "CONSERVATIVE_COARSEN");
    d_calgs["F->F::CONSERVATIVE_COARSEN"]->
        registerCoarsen(d_F_idx, // destination
                        d_F_idx, // source
                        coarsen_operator);

    // Set the current integration time.
    if (!SAMRAI::tbox::RestartManager::getManager()->isFromRestart())
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }

    // Indicate that the integrator has been initialized.
    d_is_initialized = true;

    t_initialize_hierarchy_integrator->stop();
    return;
}// initializeHierarchyIntegrator

double
IBImplicitHierarchyIntegrator::initializeHierarchy()
{
    t_initialize_hierarchy->start();

    // Use the INSHierarchyIntegrator to initialize the patch
    // hierarchy.
    double dt_next = d_ins_hier_integrator->initializeHierarchy();

    if (d_integrator_time >= d_dt_max_time_min &&
        d_integrator_time <= d_dt_max_time_max)
    {
        dt_next = SAMRAI::tbox::Utilities::dmin(dt_next, d_dt_max);
    }

    t_initialize_hierarchy->stop();
    return dt_next;
}// initializeHierarchy

double
IBImplicitHierarchyIntegrator::advanceHierarchy(
    const double _dt,
    const bool rebalance_coarsest)
{
    t_advance_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
#endif

    dt = _dt;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln == finest_ln);

    const double current_time = d_integrator_time;
    const double new_time     = d_integrator_time+dt;
    const bool initial_time   = SAMRAI::tbox::Utilities::deq(current_time,d_start_time);

    // Set the current time interval in the force specification objects.
    d_eulerian_force_set->setTimeInterval(current_time, new_time);
    d_force_strategy->setTimeInterval(current_time, new_time);

    // Rebalance the coarsest level (when requested).
    if (rebalance_coarsest)
    {
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): rebalancing coarsest level\n";
        rebalanceCoarsestLevel();
    }

    // (Re)initialize the force strategy.
    if (d_force_strategy_needs_init)
    {
        resetLagrangianForceStrategy(current_time, initial_time);
        d_force_strategy_needs_init = false;
    }

    // Initialize the Lagrangian data vectors.
    X_data.resize(finest_ln+1);
    F_data.resize(finest_ln+1);
    U_data.resize(finest_ln+1);
    X_new_data.resize(finest_ln+1);
    F_new_data.resize(finest_ln+1);
    U_new_data.resize(finest_ln+1);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln] = d_lag_data_manager->getLNodeLevelData("X",ln);
            F_data[ln] = d_lag_data_manager->createLNodeLevelData("F",ln,NDIM);
            U_data[ln] = d_lag_data_manager->createLNodeLevelData("U",ln,NDIM);

            X_data[ln]->restoreLocalFormVec();
            F_data[ln]->restoreLocalFormVec();
            U_data[ln]->restoreLocalFormVec();

            X_new_data[ln] = d_lag_data_manager->createLNodeLevelData("X_new",ln,NDIM);
            F_new_data[ln] = d_lag_data_manager->createLNodeLevelData("F_new",ln,NDIM);
            U_new_data[ln] = d_lag_data_manager->createLNodeLevelData("U_new",ln,NDIM);

            X_new_data[ln]->restoreLocalFormVec();
            F_new_data[ln]->restoreLocalFormVec();
            U_new_data[ln]->restoreLocalFormVec();
        }
    }

    // Allocate Cartesian grid patch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_V_idx, current_time);
        level->allocatePatchData(d_W_idx, current_time);
        level->allocatePatchData(d_F_idx, current_time);
        level->allocatePatchData(d_F_scratch1_idx, current_time);
        level->allocatePatchData(d_F_scratch2_idx, current_time);
    }

    // Compute F(n) and F(n+1).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // On the coarsest level of the patch hierarchy, simply
        // initialize the Cartesian force density to equal zero.
        //
        // For each of the finer levels in the patch hierarchy,
        // conservatively interpolate the Cartesian force density from
        // coarser levels in the patch hierarchy.
        if (ln == coarsest_ln)
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_current_data = patch->
                    getPatchData(d_ins_hier_integrator->getForceVar(),
                                 d_ins_hier_integrator->getCurrentContext());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_new_data = patch->
                    getPatchData(d_F_idx);
                f_current_data->fillAll(0.0);
                f_new_data->fillAll(0.0);
            }
        }
        else
        {
            // Interpolate the Cartesian force from the next coarser
            // level in the hierarchy.
            //
            // Note that these refine schedules initialize both the
            // force data maintained by the IBHierarchyIntegrator and
            // the current context of the force data maintained by the
            // Navier-Stokes solver.
            d_force_rscheds[ln]->fillData(current_time);
        }

        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            int ierr;

            // 1. Compute F(n+1) = F(X(n+1),n), the Lagrangian force
            //    corresponding to configuration X(n+1) at time t_{n+1}.
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing F(n+1) on level number " << ln << "\n";

            Vec F_new_vec = F_new_data[ln]->getGlobalVec();
            ierr = VecSet(F_new_vec, 0.0);  PETSC_SAMRAI_ERROR(ierr);
            d_force_strategy->computeLagrangianForce(
                F_new_data[ln], X_data[ln], NULL,
                d_hierarchy, ln, new_time, d_lag_data_manager);
            F_new_vec = F_new_data[ln]->getGlobalVec();
            ierr = VecScale(F_new_vec, 2.0);  PETSC_SAMRAI_ERROR(ierr);

            // 2. Communicate ghost data.
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): communicating ghost node data on level number " << ln << "\n";

            X_data[ln]->beginGhostUpdate();
            F_new_data[ln]->beginGhostUpdate();
            X_data[ln]->endGhostUpdate();
            F_new_data[ln]->endGhostUpdate();

            // 3. Spread F(n+1) from the Lagrangian mesh onto the Cartesian
            //    grid.
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): spreading F(n+1) on level number " << ln << "\n";

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_new_data = patch->
                    getPatchData(d_F_idx);
                const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data = patch->
                    getPatchData(d_lag_data_manager->
                                 getLNodeIndexPatchDescriptorIndex());
                LEInteractor::spread(
                    f_new_data, F_new_data[ln], X_data[ln], idx_data,
                    patch, SAMRAI::hier::Box<NDIM>::grow(patch_box,d_ghosts),
                    d_delta_fcn, ENFORCE_PERIODIC_BCS);
            }
        }
    }

    // Solve for u(n+1) and p(n+1/2).
    d_ins_hier_integrator->predictAdvectionVelocity(current_time, new_time);
    d_ins_hier_integrator->integrateAdvDiff(current_time, new_time);
    d_ins_hier_integrator->projectVelocity(current_time, new_time);
    d_ins_hier_integrator->updatePressure(current_time, new_time);
    d_ins_hier_integrator->synchronizeHierarchy();

    // Synchronize the Cartesian grid velocity u(n+1) on the patch hierarchy.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["U->U::N->N::CONSERVATIVE_COARSEN"][ln]->coarsenData();
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_rscheds["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(current_time);
    }

    // In the following loop, we perform the following operations:
    //
    // 1. Interpolate u(n+1) from the Cartesian grid onto the Lagrangian mesh.
    //
    // 2. Compute X(n+1), the final structure configuration at time t_{n+1}.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            // 1. Interpolate u(n+1) from the Cartesian grid onto the
            //    Lagrangian mesh.
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): interpolating u(n+1) to U(n+1) on level number " << ln << "\n";

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > w_data =
                    patch->getPatchData(d_W_idx);
                const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data = patch->getPatchData(
                    d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());

                LEInteractor::interpolate(
                    U_new_data[ln], X_data[ln], idx_data, w_data,
                    patch, patch_box,
                    d_delta_fcn, ENFORCE_PERIODIC_BCS);
            }
        }
    }

    // Solve for X(n+1), u(n+1), and p(n+1/2).
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing X(n+1), u(n+1), p(n+1/2)\n";

    int ierr;
    const int ln = 0;
    Vec X_vec = X_data[ln]->getGlobalVec();

    Vec X_soln_vec, R_soln_vec;
    ierr = VecDuplicate(X_vec, &X_soln_vec);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecDuplicate(X_vec, &R_soln_vec);  PETSC_SAMRAI_ERROR(ierr);

    ierr = VecCopy(X_vec, X_soln_vec);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecSet(R_soln_vec, 0.0);  PETSC_SAMRAI_ERROR(ierr);

    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes);  PETSC_SAMRAI_ERROR(ierr);
    ierr = SNESSetFunction(snes, R_soln_vec, IBImplicitHierarchyIntegrator::FormFunction_PETSC, static_cast<void*>(this));  PETSC_SAMRAI_ERROR(ierr);
    ierr = SNESSetFromOptions(snes);  PETSC_SAMRAI_ERROR(ierr);
    ierr = SNESSolve(snes, PETSC_NULL, X_soln_vec);  PETSC_SAMRAI_ERROR(ierr);
    ierr = SNESDestroy(snes);  PETSC_SAMRAI_ERROR(ierr);

    X_vec = X_data[ln]->getGlobalVec();
    ierr = VecCopy(X_soln_vec, X_vec);  PETSC_SAMRAI_ERROR(ierr);

    ierr = VecDestroy(X_soln_vec);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecDestroy(R_soln_vec);  PETSC_SAMRAI_ERROR(ierr);

    // Deallocate Cartesian grid patch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_V_idx);
        level->deallocatePatchData(d_W_idx);
        level->deallocatePatchData(d_F_idx);
        level->deallocatePatchData(d_F_scratch1_idx);
        level->deallocatePatchData(d_F_scratch2_idx);
    }

    // Reset all time dependent data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierData(new_time);

    // Regrid (when appropriate).
    const bool do_regrid = (d_regrid_interval == 0
                            ? false
                            : (d_integrator_step % d_regrid_interval == 0));
    if (do_regrid)
    {
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): regridding\n";
        regridHierarchy();
    }

    // Determine the next stable timestep.
    double dt_next = d_ins_hier_integrator->getStableTimestep();

    if (d_integrator_time+dt_next >= d_end_time)
    {
        dt_next = d_end_time - d_integrator_time;
    }

    if (d_integrator_time >= d_dt_max_time_min &&
        d_integrator_time <= d_dt_max_time_max)
    {
        dt_next = SAMRAI::tbox::Utilities::dmin(d_dt_max,dt_next);
    }

    if (!initial_time)
    {
        dt_next = SAMRAI::tbox::Utilities::dmin(dt_next,d_grow_dt*d_old_dt);
    }

    t_advance_hierarchy->stop();
    return dt_next;
}// advanceHierarchy

bool
IBImplicitHierarchyIntegrator::atRegridPoint() const
{
    const int level_number = 0;

    return ((d_integrator_step>0)
            && d_gridding_alg->levelCanBeRefined(level_number)
            && (d_regrid_interval == 0
                ? false
                : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double
IBImplicitHierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}// getIntegratorTime

double
IBImplicitHierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}// getStartTime

double
IBImplicitHierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}// getEndTime

int
IBImplicitHierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}// getIntegratorStep

int
IBImplicitHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}// getMaxIntegratorSteps

bool
IBImplicitHierarchyIntegrator::stepsRemaining() const
{
    return d_integrator_step < d_max_integrator_steps;
}// stepsRemaining

const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
IBImplicitHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
IBImplicitHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

LDataManager*
IBImplicitHierarchyIntegrator::getLDataManager() const
{
    return d_lag_data_manager;
}// getLDataManager

///
///  The following routines:
///
///      rebalanceCoarsestLevel(),
///      regridHierarchy(),
///      synchronizeHierarchy(),
///      synchronizeNewLevels(),
///      resetTimeDependentHierData(),
///      resetHierDataToPreadvanceState()
///
///  allow the IBImplicitHierarchyIntegrator to provide data management
///  for a time integrator which making use of this class.
///

void
IBImplicitHierarchyIntegrator::rebalanceCoarsestLevel()
{
    t_rebalance_coarsest_level->start();

    // Update the workload for regridding.
    d_lag_data_manager->updateWorkloadData(
        0,d_hierarchy->getFinestLevelNumber());

    // Before rebalancing, begin Lagrangian data movement.
    d_lag_data_manager->beginDataRedistribution();

    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->rebalanceCoarsestLevel();

    // After rebalancing, finish Lagrangian data movement.
    d_lag_data_manager->endDataRedistribution();

    // Update the workload post-regridding.
    d_lag_data_manager->updateWorkloadData(
        0,d_hierarchy->getFinestLevelNumber());

    t_rebalance_coarsest_level->stop();
    return;
}// rebalanceCoarsestLevel

void
IBImplicitHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    // Update the workload pre-regridding.
    d_lag_data_manager->updateWorkloadData(
        0,d_hierarchy->getFinestLevelNumber());

    // Before regriding, begin Lagrangian data movement.
    d_lag_data_manager->beginDataRedistribution();

    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->regridHierarchy();

    // After regridding, finish Lagrangian data movement.
    d_lag_data_manager->endDataRedistribution();

    // Update the workload post-regridding.
    d_lag_data_manager->updateWorkloadData(
        0,d_hierarchy->getFinestLevelNumber());

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

void
IBImplicitHierarchyIntegrator::synchronizeHierarchy()
{
    t_synchronize_hierarchy->start();

    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->synchronizeHierarchy();

    t_synchronize_hierarchy->stop();
    return;
}// synchronizeHierarchy

void
IBImplicitHierarchyIntegrator::synchronizeNewLevels(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const bool initial_time)
{
    t_synchronize_new_levels->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level < finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->
        synchronizeNewLevels(hierarchy, coarsest_level, finest_level,
                             sync_time, initial_time);

    t_synchronize_new_levels->stop();
    return;
}// synchronizeNewLevels

void
IBImplicitHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    t_reset_time_dependent_data->start();

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->resetTimeDependentHierData(new_time);

    t_reset_time_dependent_data->stop();
    return;
}// resetTimeDependentHierData

void
IBImplicitHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    t_reset_data_to_preadvance_state->start();

    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->resetHierDataToPreadvanceState();

    t_reset_data_to_preadvance_state->stop();
    return;
}// resetHierDataToPreadvanceState

///
///  The following routines:
///
///      initializeLevelData(),
///      resetHierarchyConfiguration(),
///      applyGradientDetector()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::mesh::StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
IBImplicitHierarchyIntegrator::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->
        initializeLevelData(hierarchy, level_number, init_data_time,
                            can_be_refined, initial_time, old_level,
                            allocate_data);

    // We use the LDataManager to handle as much unstructured data
    // management as possible.
    d_lag_data_manager->setPatchHierarchy(d_hierarchy);
    d_lag_data_manager->resetLevels(0,hierarchy->getFinestLevelNumber());

    d_lag_data_manager->
        initializeLevelData(hierarchy, level_number, init_data_time,
                            can_be_refined, initial_time, old_level,
                            allocate_data);

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBImplicitHierarchyIntegrator::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level <= finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->
        resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // We use the LDataManager to handle as much unstructured data
    // management as possible.
    d_lag_data_manager->
        resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Reset the Hierarchy data operations for the new hierarchy
    // configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    // If we have added or removed a level, resize the schedule
    // vectors.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it!= d_ralgs.end(); ++it)
    {
        d_rscheds[(*it).first].resize(finest_hier_level+1);
    }

    d_force_rscheds.resize(finest_hier_level+1);

    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it!= d_calgs.end(); ++it)
    {
        d_cscheds[(*it).first].resize(finest_hier_level+1);
    }

    // (Re)build generic refine communication schedules.  These are
    // created for all levels in the hierarchy.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it!= d_ralgs.end(); ++it)
    {
        for (int ln = coarsest_level; ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

            d_rscheds[(*it).first][ln] = (*it).second->
                createSchedule(level, ln-1, hierarchy, d_rstrategies[(*it).first]);
        }
    }

    // (Re)build specialized refine communication schedules used to
    // compute the Cartesian force density.  These are set only for
    // levels >= 1.
    for (int ln = SAMRAI::tbox::Utilities::imax(1,coarsest_level);
         ln <= finest_hier_level; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

        d_force_rscheds[ln] = d_force_ralg->
            createSchedule(
                level, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(),
                ln-1, hierarchy, d_force_rstrategy);
    }

    // (Re)build coarsen communication schedules.  These are set only
    // for levels >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it!= d_calgs.end(); ++it)
    {
        for (int ln = SAMRAI::tbox::Utilities::imax(coarsest_level,1);
             ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level =
                hierarchy->getPatchLevel(ln-1);

            d_cscheds[(*it).first][ln] = (*it).second->
                createSchedule(coarser_level, level, d_cstrategies[(*it).first]);
        }
    }

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
IBImplicitHierarchyIntegrator::applyGradientDetector(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    t_apply_gradient_detector->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // It is necessary to untag all cells prior to tagging.
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
        tags_data->fillAll(0);
    }

    // Tag cells for refinement according to the criteria specified by
    // the INSHierarchyIntegrator.
    d_ins_hier_integrator->
        applyGradientDetector(hierarchy, level_number, error_data_time,
                              tag_index, initial_time,
                              uses_richardson_extrapolation_too);

    // Tag cells which contain Lagrangian nodes.
    d_lag_data_manager->
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
/// We simply reuse the SAMRAI::hier::VariableContext objects defined in the
/// INSHierarchyIntegrator object.
///

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBImplicitHierarchyIntegrator::getCurrentContext() const
{
    return d_ins_hier_integrator->getCurrentContext();
}// getCurrentContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBImplicitHierarchyIntegrator::getNewContext() const
{
    return d_ins_hier_integrator->getNewContext();
}// getNewContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBImplicitHierarchyIntegrator::getOldContext() const
{
    return d_ins_hier_integrator->getOldContext();
}// getOldContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBImplicitHierarchyIntegrator::getScratchContext() const
{
    return d_ins_hier_integrator->getScratchContext();
}// getScratchContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBImplicitHierarchyIntegrator::getPlotContext() const
{
    return d_ins_hier_integrator->getPlotContext();
}// getPlotContext

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::tbox::Serializable abstract base class.
///

void
IBImplicitHierarchyIntegrator::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION",
                   IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION);

    db->putString("d_delta_fcn", d_delta_fcn);
    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putDouble("d_old_dt", d_old_dt);
    db->putDouble("d_integrator_time", d_integrator_time);
    db->putInteger("d_integrator_step", d_integrator_step);
    db->putDouble("d_dt_max", d_dt_max);
    db->putDouble("d_dt_max_time_max", d_dt_max_time_max);
    db->putDouble("d_dt_max_time_min", d_dt_max_time_min);

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

void
IBImplicitHierarchyIntegrator::printClassData(
    std::ostream& os) const
{
    os << "\nIBImplicitHierarchyIntegrator::printClassData..." << "\n";
    os << "\nIBImplicitHierarchyIntegrator: this = "
       << const_cast<IBImplicitHierarchyIntegrator*>(this) << "\n";
    os << "d_object_name = " << d_object_name << "\n";
    os << "d_integrator_time = " << d_integrator_time << "\n"
       << "d_start_time = " << d_start_time << "\n"
       << "d_end_time = " << d_end_time << "\n"
       << "d_integrator_step = " << d_integrator_step << "\n"
       << "d_grow_dt = " << d_grow_dt << "\n";
    os << "I AM INCOMPLETE!!!!!!!!!" << "\n";
    return;
}// printClassData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBImplicitHierarchyIntegrator::resetLagrangianForceStrategy(
    const double init_data_time,
    const bool initial_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            d_force_strategy->initializeLevelData(
                d_hierarchy, ln, init_data_time, initial_time,
                d_lag_data_manager);
        }
    }

    return;
}// resetLagrangianForceStrategy

void
IBImplicitHierarchyIntegrator::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    // Read data members from input database.
    d_end_time = db->getDoubleWithDefault("end_time", d_end_time);
    d_grow_dt = db->getDoubleWithDefault("grow_dt", d_grow_dt);

    d_max_integrator_steps = db->getIntegerWithDefault(
        "max_integrator_steps", d_max_integrator_steps);

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault(
        "dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault(
        "dt_max_time_min", d_dt_max_time_min);

    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    if (!is_from_restart)
    {
        d_delta_fcn = db->getStringWithDefault("delta_fcn", d_delta_fcn);

        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);
    }

    return;
}// getFromInput

void
IBImplicitHierarchyIntegrator::getFromRestart()
{
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> restart_db =
        SAMRAI::tbox::RestartManager::getManager()->getRootDatabase();

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    d_delta_fcn = db->getString("d_delta_fcn");
    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    d_grow_dt = db->getDouble("d_grow_dt");
    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");
    d_regrid_interval = db->getInteger("d_regrid_interval");
    d_old_dt = db->getDouble("d_old_dt");
    d_integrator_time = db->getDouble("d_integrator_time");
    d_integrator_step = db->getInteger("d_integrator_step");
    d_dt_max = db->getDouble("d_dt_max");
    d_dt_max_time_max = db->getDouble("d_dt_max_time_max");
    d_dt_max_time_min = db->getDouble("d_dt_max_time_min");

    return;
}// getFromRestart

PetscErrorCode
IBImplicitHierarchyIntegrator::FormFunction_PETSC(
    SNES snes,
    Vec X_soln_vec,
    Vec R_soln_vec,
    void* p_ctx)
{
    (void) snes;
    IBImplicitHierarchyIntegrator* solver = static_cast<IBImplicitHierarchyIntegrator*>(p_ctx);
    PetscFunctionReturn(solver->FormFunction(X_soln_vec, R_soln_vec));
}// FormFunction_PETSC

PetscErrorCode
IBImplicitHierarchyIntegrator::FormFunction(
    Vec X_soln_vec,
    Vec R_soln_vec)
{
    int ierr;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln == finest_ln);
    const int ln = 0;

    // Compute X~(n+1).
    {
        Vec X_vec = X_data[ln]->getGlobalVec();
        Vec X_new_vec = X_new_data[ln]->getGlobalVec();
        Vec U_new_vec = U_new_data[ln]->getGlobalVec();

        ierr = VecCopy(X_vec, X_new_vec);  PETSC_SAMRAI_ERROR(ierr);
        ierr = VecAXPY(X_new_vec, +dt, U_new_vec);  PETSC_SAMRAI_ERROR(ierr);
    }

    // Set R = X(n+1)-X~(n+1).
    {
        Vec X_new_vec = X_new_data[ln]->getGlobalVec();
        ierr = VecWAXPY(R_soln_vec, -1.0, X_soln_vec, X_new_vec);  PETSC_SAMRAI_ERROR(ierr);
    }

    PetscFunctionReturn(0);
}// FormFunction

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBImplicitHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
