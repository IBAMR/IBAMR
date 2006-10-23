// Filename: IBHierarchyIntegrator.C
// Created on 12 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <23.Oct.2006 17:35:56 boyce@bigboy.nyconnect.com>

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "IBHierarchyIntegrator.h"

// IBAMR INCLUDES
#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#endif

#include <ibamr/LEInteractor.h>
#include <ibamr/LNodeIndexData.h>

// STOOLS INCLUDES
#include <stools/PETSC_SAMRAI_ERROR.h>
#include <stools/STOOLS_Utilities.h>

// SAMRAI INCLUDES
#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CellIndex.h>
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
#include <cassert>
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

// The delta function to use when interpolating point pressure
// values from the grid (these are required to compute source/sink
// flow rates).
static const std::string PRES_DELTA_FCN = "IB_4";

// The delta function to use when spreading point source/sink flow
// rates onto the grid.
static const std::string SOURCE_DELTA_FCN = "WIDE_IB_4";

// Name of normalization output file.
static const std::string NORM_DATA_FILE_NAME = "source_norm_data.txt";

// Whether to account for periodic shifts in spreading and
// interpolating.
static const bool ENFORCE_PERIODIC_BCS = true;

// Version of IBHierarchyIntegrator restart file data.
static const int IB_HIERARCHY_INTEGRATOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBHierarchyIntegrator::IBHierarchyIntegrator(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    SAMRAI::tbox::Pointer<INSHierarchyIntegrator> ins_hier_integrator,
    SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> force_strategy,
    SAMRAI::tbox::Pointer<IBLagrangianSourceStrategy> source_strategy,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_hierarchy(hierarchy),
      d_gridding_alg(NULL),
      d_visit_writer(NULL),
      d_silo_writer(NULL),
      d_load_balancer(NULL),
      d_ins_hier_integrator(ins_hier_integrator),
      d_lag_data_manager(NULL),
      d_lag_posn_init(NULL),
      d_body_force_set(NULL),
      d_eulerian_force_set(NULL),
      d_force_strategy(force_strategy),
      d_force_strategy_needs_init(true),
      d_eulerian_source_set(NULL),
      d_source_strategy(source_strategy),
      d_source_strategy_needs_init(true),
      d_delta_fcn("IB_4"),
      d_ghosts(-1),
      d_pres_ghosts(-1),
      d_source_ghosts(-1),
      d_start_time(0.0),
      d_end_time(numeric_limits<double>::max()),
      d_grow_dt(2.0),
      d_max_integrator_steps(numeric_limits<int>::max()),
      d_num_cycles(1),
      d_num_init_cycles(5),
      d_regrid_interval(1),
      d_timestepping_order(2),
      d_old_dt(-1.0),
      d_integrator_time(numeric_limits<double>::quiet_NaN()),
      d_integrator_step(numeric_limits<int>::max()),
      d_dt_max(numeric_limits<double>::max()),
      d_dt_max_time_max(numeric_limits<double>::max()),
      d_dt_max_time_min(-(d_dt_max_time_max-numeric_limits<double>::epsilon())),
      d_is_initialized(false),
      d_do_log(false),
      d_ralgs(),
      d_rscheds(),
      d_V_var(NULL),
      d_W_var(NULL),
      d_F_var(NULL),
      d_P_var(NULL),
      d_Q_var(NULL),
      d_context(NULL),
      d_V_idx(-1),
      d_W_idx(-1),
      d_F_idx(-1),
      d_P_idx(-1),
      d_Q_idx(-1)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
    assert(!hierarchy.isNull());
    assert(!ins_hier_integrator.isNull());
    assert(!force_strategy.isNull());
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

    // Set the ghost cell widths required for spreading and
    // interpolation.
    const int stencil_size = LEInteractor::getStencilSize(d_delta_fcn);
    d_ghosts = static_cast<int>(floor(0.5*static_cast<double>(stencil_size))+1);

    const int pres_stencil_size = LEInteractor::getStencilSize(PRES_DELTA_FCN);
    d_pres_ghosts = static_cast<int>(floor(0.5*static_cast<double>(pres_stencil_size))+1);

    const int source_stencil_size = LEInteractor::getStencilSize(SOURCE_DELTA_FCN);
    d_source_ghosts = static_cast<int>(floor(0.5*static_cast<double>(source_stencil_size))+1);

    // Get the Lagrangian Data Manager.
    //
    // XXXX: Make sure this ghost cell width is correct!!!!
    d_lag_data_manager = LDataManager::getManager(
        d_object_name+"::LDataManager", d_ghosts, d_registered_for_restart);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::advanceHierarchy()");
        t_rebalance_coarsest_level = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::rebalanceCoarsestLevel()");
        t_regrid_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::regridHierarchy()");
        t_integrate_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::integrateHierarchy()");
        t_synchronize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::IBHierarchyIntegrator::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// IBHierarchyIntegrator

IBHierarchyIntegrator::~IBHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
    return;
}// ~IBHierarchyIntegrator

void
IBHierarchyIntegrator::registerBodyForceSpecification(
    SAMRAI::tbox::Pointer<SetDataStrategy> body_force_set)
{
    d_body_force_set = body_force_set;
    if (!d_eulerian_force_set.isNull())
    {
        d_eulerian_force_set->registerBodyForce(d_body_force_set);
    }
    return;
}// registerBodyForceSpecification

void
IBHierarchyIntegrator::registerLNodePosnInitStrategy(
    SAMRAI::tbox::Pointer<LNodePosnInitStrategy> lag_posn_init)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!lag_posn_init.isNull());
#endif
    d_lag_posn_init = lag_posn_init;
    d_lag_data_manager->registerLNodePosnInitStrategy(d_lag_posn_init);
    return;
}// registerLNodePosnInitStrategy

void
IBHierarchyIntegrator::freeLNodePosnInitStrategy()
{
    d_lag_posn_init.setNull();
    d_lag_data_manager->freeLNodePosnInitStrategy();
    return;
}// freeLNodePosnInitStrategy

void
IBHierarchyIntegrator::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    d_ins_hier_integrator->registerVisItDataWriter(d_visit_writer);
    d_lag_data_manager->registerVisItDataWriter(d_visit_writer);
    return;
}// registerVisItDataWriter

void
IBHierarchyIntegrator::registerLagSiloDataWriter(
    SAMRAI::tbox::Pointer<LagSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!silo_writer.isNull());
#endif
    d_silo_writer = silo_writer;
    d_lag_data_manager->registerLagSiloDataWriter(d_silo_writer);
    return;
}// registerLagSiloDataWriter

void
IBHierarchyIntegrator::registerLoadBalancer(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!load_balancer.isNull());
#endif
    d_load_balancer = load_balancer;
    d_lag_data_manager->registerLoadBalancer(d_load_balancer);
    return;
}// registerLoadBalancer

void
IBHierarchyIntegrator::gatherAllData(
    const int mpi_root,
    double* const X_structure,
    const int struct_sz,
    double* const X_marker,
    const int marker_sz,
    const int level_num_in)
{
    // XXXX: This should be cleaned up a lot!

    const int level_num = (level_num_in == -1 ? d_hierarchy->getFinestLevelNumber() : level_num_in);

    int ierr;

    // Get needed MPI info.
    const int mpi_rank = SAMRAI::tbox::MPI::getRank();

    // Setup the IS data used to generate the VecScatters that
    // redistribute the distributed data into a single local vector on
    // the specified MPI process
    const int total_sz = (mpi_rank == mpi_root ? struct_sz + marker_sz : 0);
    std::vector<int> ref_is_idxs(total_sz);
    for (int k = 0; k < total_sz; ++k)
    {
        ref_is_idxs[k] = k;
    }
    std::vector<int> dst_is_idxs = ref_is_idxs;

    // Map Lagrangian indices to PETSc indices.
    d_lag_data_manager->mapLagrangianToPETSc(ref_is_idxs, level_num);

    // Setup IS indices for the appropriate data depth.
    std::vector<int> src_is_idxs(total_sz);
    transform(ref_is_idxs.begin(), ref_is_idxs.end(),
              src_is_idxs.begin(),
              bind2nd(multiplies<int>(),NDIM));
    transform(dst_is_idxs.begin(), dst_is_idxs.end(),
              dst_is_idxs.begin(),
              bind2nd(multiplies<int>(),NDIM));

    SAMRAI::tbox::Pointer<LNodeLevelData> X_data = d_lag_data_manager->getLNodeLevelData("X", level_num);
    Vec& src_vec = X_data->getGlobalVec();

    // Create the VecScatter to scatter data from the global
    // PETSc Vec.
    IS         src_is, dst_is;
    Vec        dst_vec;
    VecScatter vec_scatter;
    ierr = ISCreateBlock(PETSC_COMM_WORLD, NDIM, total_sz,
                         &src_is_idxs[0], &src_is);  PETSC_SAMRAI_ERROR(ierr);
    ierr = ISCreateBlock(PETSC_COMM_WORLD, NDIM, total_sz,
                         &dst_is_idxs[0], &dst_is);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, NDIM*total_sz,
                        PETSC_DETERMINE, &dst_vec);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecSetBlockSize(dst_vec, NDIM);           PETSC_SAMRAI_ERROR(ierr);
    ierr = VecScatterCreate(src_vec, src_is, dst_vec, dst_is,
                            &vec_scatter);           PETSC_SAMRAI_ERROR(ierr);

    // Scatter the data.
    ierr = VecScatterBegin(src_vec, dst_vec, INSERT_VALUES, SCATTER_FORWARD,
                           vec_scatter);             PETSC_SAMRAI_ERROR(ierr);
    ierr = VecScatterEnd(  src_vec, dst_vec, INSERT_VALUES, SCATTER_FORWARD,
                           vec_scatter);             PETSC_SAMRAI_ERROR(ierr);

    // Copy the data into the provided arrays.
    double* X_arr;
    ierr = VecGetArray(    dst_vec, &X_arr);         PETSC_SAMRAI_ERROR(ierr);
    if (mpi_rank == mpi_root)
    {
        memcpy(static_cast<void*>(&X_structure[0]), static_cast<void*>(&X_arr[             0]), NDIM*struct_sz*sizeof(double));
        memcpy(static_cast<void*>(&X_marker   [0]), static_cast<void*>(&X_arr[NDIM*struct_sz]), NDIM*marker_sz*sizeof(double));
    }
    ierr = VecRestoreArray(dst_vec, &X_arr);         PETSC_SAMRAI_ERROR(ierr);

    // Destroy all workspace data.
    ierr = ISDestroy(src_is);               PETSC_SAMRAI_ERROR(ierr);
    ierr = ISDestroy(dst_is);               PETSC_SAMRAI_ERROR(ierr);
    ierr = VecDestroy(dst_vec);             PETSC_SAMRAI_ERROR(ierr);
    ierr = VecScatterDestroy(vec_scatter);  PETSC_SAMRAI_ERROR(ierr);

    return;
}// gatherAllData

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
///      getGriddingAlgorithm()
///
///  allow the IBHierarchyIntegrator to be used as a hierarchy
///  integrator.
///

void
IBHierarchyIntegrator::initializeHierarchyIntegrator(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Initialize all variables.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name+"::CONTEXT");
    const SAMRAI::hier::IntVector<NDIM> ghosts = d_ghosts;
    const SAMRAI::hier::IntVector<NDIM> pres_ghosts = d_pres_ghosts;
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    d_V_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::V",NDIM);
    d_V_idx = var_db->registerVariableAndContext(d_V_var, d_context, ghosts);

    d_W_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::W",NDIM);
    d_W_idx = var_db->registerVariableAndContext(d_W_var, d_context, ghosts);

    d_F_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::F",NDIM);
    d_F_idx = var_db->registerVariableAndContext(d_F_var, d_context, no_ghosts);

    if (!d_source_strategy.isNull())
    {
        d_P_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::P",1);
        d_P_idx = var_db->registerVariableAndContext(d_P_var, d_context, pres_ghosts);

        d_Q_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Q",1);
        d_Q_idx = var_db->registerVariableAndContext(d_Q_var, d_context, no_ghosts);
    }

    // Initialize the INSHierarchyIntegrator.
    //
    // NOTE: This must be done after all variables are registered.
    d_eulerian_force_set = new IBEulerianForceSetter(
        d_object_name+"::IBEulerianForceSetter",d_F_idx);
    d_eulerian_force_set->registerBodyForce(d_body_force_set);
    d_ins_hier_integrator->registerForceSpecification(d_eulerian_force_set);

    if (!d_source_strategy.isNull())
    {
        d_eulerian_source_set = new IBEulerianSourceSetter(
            d_object_name+"::IBEulerianSourceSetter",d_Q_idx);
        d_ins_hier_integrator->registerDivergenceSpecification(
            d_eulerian_source_set);
    }

    d_ins_hier_integrator->initializeHierarchyIntegrator(d_gridding_alg);

    // Create several refinement communications algorithms, used in
    // filling ghost cell data.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;

    d_ralgs["U_current->V::CONSERVATIVE_LINEAR_REFINE"] =
        new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getVelocityVar(),
        "CONSERVATIVE_LINEAR_REFINE");

    int U_current_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getCurrentContext());

    d_ralgs["U_current->V::CONSERVATIVE_LINEAR_REFINE"]->
        registerRefine(d_V_idx,       // destination
                       U_current_idx, // source
                       d_V_idx,       // temporary work space
                       refine_operator);

    d_ralgs["U_new->W::CONSERVATIVE_LINEAR_REFINE"] =
        new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getVelocityVar(),
        "CONSERVATIVE_LINEAR_REFINE");

    int U_new_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getNewContext());

    d_ralgs["U_new->W::CONSERVATIVE_LINEAR_REFINE"]->
        registerRefine(d_W_idx,   // destination
                       U_new_idx, // source
                       d_W_idx,   // temporary work space
                       refine_operator);

    if (!d_source_strategy.isNull())
    {
        d_ralgs["P_current->P::LINEAR_REFINE"] =
            new SAMRAI::xfer::RefineAlgorithm<NDIM>();

        refine_operator = grid_geom->lookupRefineOperator(
            d_ins_hier_integrator->getPressureVar(),
            "LINEAR_REFINE");

        int P_current_idx = var_db->mapVariableAndContextToIndex(
            d_ins_hier_integrator->getPressureVar(),
            d_ins_hier_integrator->getCurrentContext());

        d_ralgs["P_current->P::LINEAR_REFINE"]->
            registerRefine(d_P_idx,       // destination
                           P_current_idx, // source
                           d_P_idx,       // temporary work space
                           refine_operator);
    }

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
IBHierarchyIntegrator::initializeHierarchy()
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
IBHierarchyIntegrator::advanceHierarchy(
    const double dt,
    const bool rebalance_coarsest)
{
    t_advance_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_end_time >= d_integrator_time+dt);
#endif

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Rebalance the coarsest level (when requested).
    if (rebalance_coarsest)
    {
        rebalanceCoarsestLevel();
    }

    const double current_time = d_integrator_time;
    const double new_time     = d_integrator_time+dt;
    const bool initial_time   = SAMRAI::tbox::Utilities::deq(current_time,d_start_time);

    // (Re)initialize the force and source strategies.
    if (d_force_strategy_needs_init)
    {
        resetLagrangianForceStrategy(current_time, initial_time);
    }
    if (d_source_strategy_needs_init && !d_source_strategy.isNull())
    {
        resetLagrangianSourceStrategy(current_time, initial_time);
    }

    d_force_strategy_needs_init = false;
    d_source_strategy_needs_init = false;

    std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > X_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > F_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > U_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > X_new_data(finest_ln+1);

    // When necessary, determine the average pressure near the
    // periodic boundary at z=z_min and z=z_max.
    double p_norm = 0.0;
    if (!d_source_strategy.isNull())
    {
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        if (!grid_geom->getDomainIsSingleBox()) TBOX_ERROR("physical domain must be a single box...\n");

        const SAMRAI::hier::Box<NDIM>& domain_box = grid_geom->getPhysicalDomain()(0);
        SAMRAI::hier::Box<NDIM> lower_box = domain_box;
        SAMRAI::hier::Box<NDIM> upper_box = domain_box;
        lower_box.upper(NDIM-1) = domain_box.lower(NDIM-1);
        upper_box.lower(NDIM-1) = domain_box.upper(NDIM-1);

        const double* const ref_xLower = grid_geom->getXLower();
        const double* const ref_xUpper = grid_geom->getXUpper();
        const double* const ref_dx = grid_geom->getDx();
        double vol = 1.0;
        for (int d = 0; d < NDIM-1; ++d)
        {
            vol *= (ref_xUpper[d] - ref_xLower[d]);
        }
        vol *= 2.0*ref_dx[NDIM-1];

        // Compute the discrete integral of the pressure near the
        // periodic boundary.
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        const int P_current_idx = var_db->mapVariableAndContextToIndex(
            d_ins_hier_integrator->getPressureVar(),
            d_ins_hier_integrator->getCurrentContext());
        const int wgt_idx = d_ins_hier_integrator->
            getHierarchyMathOps()->getCellWeightPatchDescriptorIndex();
        SAMRAI::math::PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
        double p_integral = 0.0;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            const SAMRAI::hier::IntVector<NDIM>& ratio = level->getRatio();
            const SAMRAI::hier::Box<NDIM> refined_lower_box =
                SAMRAI::hier::Box<NDIM>::refine(lower_box,ratio);
            const SAMRAI::hier::Box<NDIM> refined_upper_box =
                SAMRAI::hier::Box<NDIM>::refine(upper_box,ratio);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > p_data = patch->
                    getPatchData(P_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > wgt_data = patch->
                    getPatchData(wgt_idx);
                p_integral += patch_cc_data_ops.integral(
                    p_data, patch_box*refined_lower_box, wgt_data);
                p_integral += patch_cc_data_ops.integral(
                    p_data, patch_box*refined_upper_box, wgt_data);
            }
        }

        // Compute the normalization pressure.
        p_norm = SAMRAI::tbox::MPI::sumReduction(p_integral)/vol;
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): p_norm = " << p_norm << "\n";
    }

    // Compute F(n) = F(X(n),n), the preliminary structure
    // configuration X~(n+1), and F~(n+1) = F(X~(n+1),n+1).
    //
    // Also, when specified, compute Q(n), the source/sink strength.
    double q_total = 0.0;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_V_idx, current_time);
        level->allocatePatchData(d_W_idx, current_time);
        level->allocatePatchData(d_F_idx, current_time);

        if (!d_source_strategy.isNull())
        {
            level->allocatePatchData(d_P_idx, current_time);
            level->allocatePatchData(d_Q_idx, current_time);
            d_rscheds["P_current->P::LINEAR_REFINE"][ln]->fillData(current_time);
        }

        if (!d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            // Set the Cartesian force and source/sink data to equal
            // zero.
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_data = patch->
                    getPatchData(d_F_idx);
                f_data->fillAll(0.0);

                if (!d_source_strategy.isNull())
                {
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data = patch->
                        getPatchData(d_Q_idx);
                    q_data->fillAll(0.0);
                }
            }
        }
        else
        {
            if (ln < finest_ln)
            {
                TBOX_ERROR("Forces presently are not correctly computed except on the finest level of the patch hierarchy.\n");
            }

            X_data[ln]     = d_lag_data_manager->getLNodeLevelData("X",ln);
            F_data[ln]     = d_lag_data_manager->createLNodeLevelData("F",ln,NDIM);
            U_data[ln]     = d_lag_data_manager->createLNodeLevelData("U",ln,NDIM);
            X_new_data[ln] = d_lag_data_manager->createLNodeLevelData("X_new",ln,NDIM);

            X_data[ln]    ->restoreLocalFormVec();
            F_data[ln]    ->restoreLocalFormVec();
            U_data[ln]    ->restoreLocalFormVec();
            X_new_data[ln]->restoreLocalFormVec();

            // Compute Q(n) = Q(X(n),n) and spread the source/sink
            // strengths onto the Cartesian grid.
            if (!d_source_strategy.isNull())
            {
                // Get the present source locations.
                //
                // IMPORTANT NOTE: We assume that each MPI process is
                // assigned the source locations.
                if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing source locations\n";

                const int nsrc = d_source_strategy->getNumSources();
                std::vector<double> P_src(nsrc,0.0);
                std::vector<double> Q_src(nsrc,0.0);
                std::vector<std::vector<double> > X_src(nsrc, std::vector<double>(NDIM,0.0));
                d_source_strategy->getSourceLocations(
                    X_src, X_data[ln],
                    d_hierarchy, ln, current_time, d_lag_data_manager);

                // Interpolate the Cartesian pressure to the point
                // source locations.
                //
                // For each point source/sink, we only interpolate the
                // pressures on a *single* process, namely the process
                // that owns the patch containing the location of the
                // point source/sink.
                //
                // IMPORTANT NOTE: The following code assumes that the
                // number of point source/sinks is O(1) per process.
                if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing source pressures\n";

                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

                    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                    const double* const xLower = pgeom->getXLower();
                    const double* const xUpper = pgeom->getXUpper();

                    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > p_data =
                        patch->getPatchData(d_P_idx);

                    for (int n = 0; n < nsrc; ++n)
                    {
                        bool patch_owns_src = true;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            patch_owns_src = patch_owns_src &&
                                (xLower[d] <= X_src[n][d] &&
                                 xUpper[d] >  X_src[n][d]);
                        }

                        if (patch_owns_src)
                        {
                            static const int n_vals = 1;
                            LEInteractor::interpolate(
                                &P_src[n], 1, &X_src[n][0], NDIM, n_vals, p_data,
                                patch, patch_box,
                                PRES_DELTA_FCN);
                            P_src[n] -= p_norm;
                        }
                    }
                }

                // Perform an MPI reduction to properly set the values
                // of the pressures on all processes.
                //
                // NOTE: This only works because we have interpolated
                // and normalized the pressure at a given source/sink
                // position on exactly one MPI process.
                SAMRAI::tbox::MPI::sumReduction(&P_src[0], P_src.size());

                // Compute the source strengths.
                //
                // IMPORTANT NOTE: We assume that each MPI process is
                // assigned the same source strengths.
                if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing source strengths\n";

                d_source_strategy->computeSourceStrengths(
                    Q_src, P_src, X_src, ln, current_time);
                q_total = accumulate(Q_src.begin(), Q_src.end(), q_total);

                // Spread the sources onto the Cartesian grid.
                //
                // IMPORTANT NOTE: This operation assumes that the
                // number of point sources is O(1) per process.
                if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): spreading source strengths to Cartesian grid\n";

                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

                    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                    const double* const xLower = pgeom->getXLower();
                    const double* const xUpper = pgeom->getXUpper();
                    const double* const dx = pgeom->getDx();

                    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data =
                        patch->getPatchData(d_Q_idx);
                    q_data->fillAll(0.0);

                    for (int n = 0; n < nsrc; ++n)
                    {
                        const SAMRAI::hier::IntVector<NDIM>& ghosts = d_source_ghosts;
                        bool spread_source = true;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            spread_source = spread_source &&
                                (xLower[d]-static_cast<double>(ghosts(d))*dx[d] <= X_src[n][d] &&
                                 xUpper[d]+static_cast<double>(ghosts(d))*dx[d] >= X_src[n][d]);
                        }

                        if (spread_source)
                        {
                            static const int n_vals = 1;
                            LEInteractor::spread(
                                q_data, &Q_src[n], 1, &X_src[n][0], NDIM, n_vals,
                                patch, SAMRAI::hier::Box<NDIM>::grow(patch_box,ghosts),
                                SOURCE_DELTA_FCN);
                        }
                    }
                }
            }

            // Compute F(n) = F(X(n),n) and spread the force onto the
            // Cartesian grid.
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing F(n)\n";

            {// XXXX
                Vec F_level_vec = F_data[ln]->getGlobalVec();
                int ierr = VecSet(F_level_vec, 0.0);  PETSC_SAMRAI_ERROR(ierr);
            }

            d_force_strategy->computeLagrangianForce(
                F_data[ln], X_data[ln],
                d_hierarchy, ln, current_time, d_lag_data_manager);
            X_data[ln]->beginGhostUpdate();
            F_data[ln]->beginGhostUpdate();
            X_data[ln]->endGhostUpdate();
            F_data[ln]->endGhostUpdate();

            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): spreading F(n) to f(n)\n";

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_data = patch->
                    getPatchData(d_ins_hier_integrator->getForceVar(),
                                 d_ins_hier_integrator->getCurrentContext());
                const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data = patch->
                    getPatchData(d_lag_data_manager->
                                 getLNodeIndexPatchDescriptorIndex());
                f_data->fillAll(0.0);
                LEInteractor::spread(
                    f_data, F_data[ln], X_data[ln], idx_data,
                    patch, SAMRAI::hier::Box<NDIM>::grow(patch_box,d_ghosts),
                    d_delta_fcn, ENFORCE_PERIODIC_BCS);
            }

            // Interpolate the Cartesian velocity onto the Lagrangian
            // mesh, U(X(n),n) = u(X(n),n).
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): interpolating u(n) to U(n)\n";

            d_rscheds["U_current->V::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(current_time);

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > v_data =
                    patch->getPatchData(d_V_idx);
                const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data = patch->getPatchData(
                    d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());

                LEInteractor::interpolate(
                    U_data[ln], X_data[ln], idx_data, v_data,
                    patch, patch_box,
                    d_delta_fcn, ENFORCE_PERIODIC_BCS);
            }

            // Set X~(n+1) = X(n) + dt*U(X(n),n).
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing X~(n+1)\n";

            int ierr;
            ierr = VecCopy(X_data[ln]->getGlobalVec(), X_new_data[ln]->getGlobalVec());
            PETSC_SAMRAI_ERROR(ierr);
            ierr = VecAXPY(X_new_data[ln]->getGlobalVec(),
                           dt, U_data[ln]->getGlobalVec());
            PETSC_SAMRAI_ERROR(ierr);

            // Compute F~(n+1) = F(X~(n+1),n+1) and spread the force
            // onto the Cartesian grid.
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing F~(n+1)\n";

            {// XXXX
                Vec F_level_vec = F_data[ln]->getGlobalVec();
                int ierr = VecSet(F_level_vec, 0.0);  PETSC_SAMRAI_ERROR(ierr);
            }

            d_force_strategy->computeLagrangianForce(
                F_data[ln], X_new_data[ln],
                d_hierarchy, ln, new_time, d_lag_data_manager);
            X_new_data[ln]->beginGhostUpdate();
            F_data    [ln]->beginGhostUpdate();
            X_new_data[ln]->endGhostUpdate();
            F_data    [ln]->endGhostUpdate();

            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): spreading F~(n+1) to f~(n+1)\n";

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_data = patch->
                    getPatchData(d_F_idx);
                const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data = patch->
                    getPatchData(d_lag_data_manager->
                                 getLNodeIndexPatchDescriptorIndex());
                f_data->fillAll(0.0);
                LEInteractor::spread(
                    f_data, F_data[ln], X_new_data[ln], idx_data,
                    patch,SAMRAI::hier::Box<NDIM>::grow(patch_box,d_ghosts),
                    d_delta_fcn, ENFORCE_PERIODIC_BCS);
            }
        }
    }

    // When necessary, normalize the Cartesian source/sink field.
    //
    // NOTE: We assume here that the point sources and sinks are
    // sufficiently far away from the periodic boundary that their
    // values are not modeified by the following operation.
    if (!d_source_strategy.isNull())
    {
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): q_total = " << q_total << "\n";

        // Place source/sink planes near the periodic boundary at
        // z=z_min and z=z_max so that the mean (discrete integral) of
        // the divergence field is zero.
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        if (!grid_geom->getDomainIsSingleBox()) TBOX_ERROR("physical domain must be a single box...\n");

        const SAMRAI::hier::Box<NDIM>& domain_box = grid_geom->getPhysicalDomain()(0);
        SAMRAI::hier::Box<NDIM> lower_box = domain_box;
        SAMRAI::hier::Box<NDIM> upper_box = domain_box;
        lower_box.upper(NDIM-1) = domain_box.lower(NDIM-1);
        upper_box.lower(NDIM-1) = domain_box.upper(NDIM-1);

        const double* const ref_xLower = grid_geom->getXLower();
        const double* const ref_xUpper = grid_geom->getXUpper();
        const double* const ref_dx = grid_geom->getDx();
        double vol = 1.0;
        for (int d = 0; d < NDIM-1; ++d)
        {
            vol *= (ref_xUpper[d] - ref_xLower[d]);
        }
        vol *= 2.0*ref_dx[NDIM-1];

        // Fill the appropriate portions of the domain.
        const double q_norm = -q_total/vol;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            const SAMRAI::hier::IntVector<NDIM>& ratio = level->getRatio();
            const SAMRAI::hier::Box<NDIM> refined_lower_box =
                SAMRAI::hier::Box<NDIM>::refine(lower_box,ratio);
            const SAMRAI::hier::Box<NDIM> refined_upper_box =
                SAMRAI::hier::Box<NDIM>::refine(upper_box,ratio);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data = patch->
                    getPatchData(d_Q_idx);
                q_data->fill(q_norm, patch_box*refined_lower_box);
                q_data->fill(q_norm, patch_box*refined_upper_box);
            }
        }

        // Write out the normalization data.
        static const int mpi_root = 0;
        if (SAMRAI::tbox::MPI::getRank() == mpi_root)
        {
            static bool file_initialized = false;
            if (!file_initialized)
            {
                ifstream fin(NORM_DATA_FILE_NAME.c_str(),ios::in);
                file_initialized = fin.is_open();
            }

            if (!file_initialized)
            {
                ofstream fout(NORM_DATA_FILE_NAME.c_str(),ios::out);
#if (NDIM == 2)
                fout << "mean pressure and normalizing source strength near y-periodic boundary" << endl;
#endif
#if (NDIM == 3)
                fout << "mean pressure and normalizing source strength near z-periodic boundary" << endl;
#endif
                fout << "note that -q_norm*vol equals the sum of the point source/sink strengths" << endl;

                fout.width(13);
                fout << "time";
                fout.width(14);
                fout << "volume";
                fout.width(14);
                fout << "p_norm";
                fout.width(14);
                fout << "q_norm";

                fout << endl;

                file_initialized = true;
            }

            ofstream fout(NORM_DATA_FILE_NAME.c_str(),ios::app);

            fout.setf(ios_base::showpoint);
            fout.width(13); fout.precision(10);
            fout << current_time;
            fout.setf(ios_base::scientific);
            fout.setf(ios_base::showpos);
            fout.width(14); fout.precision(5);
            fout << vol;
            fout.width(14); fout.precision(5);
            fout << p_norm;
            fout.width(14); fout.precision(5);
            fout << q_norm;
            fout.unsetf(ios_base::scientific);
            fout.unsetf(ios_base::showpos);

            fout << endl;
        }
    }

    // The pressure at start_time is not an initial value for the
    // incompressible Navier-Stokes equations, so we solve for it by
    // cycling the solution for the first timestep:
    //
    // Solve the Navier-Stokes equations with initial guess P(n=0)=0,
    // solve for P(n=1/2), discard U(n=1) and u(n=1), rinse, wash,
    // repeat.
    //
    // For all other timesteps, we just use the previous value of P as
    // the guess for P(n+1/2).
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing u(n+1), p(n+1/2)\n";

    int cycle = 0;
    const bool performing_init_cycles = initial_time;

    const int num_cycles = performing_init_cycles ?
        d_num_init_cycles : d_num_cycles;

    for (cycle = 0; cycle < num_cycles; ++cycle)
    {
        if (performing_init_cycles)
        {
            if (d_do_log) SAMRAI::tbox::plog << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++\n";
            if (d_do_log) SAMRAI::tbox::plog << "+\n";
            if (d_do_log) SAMRAI::tbox::plog << "+ Performing cycle " << cycle+1 << " of "
                                             << d_num_init_cycles << " to initialize P(n=1/2)\n";
            if (d_do_log) SAMRAI::tbox::plog << "+\n";
            if (d_do_log) SAMRAI::tbox::plog << "++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
        }

        // Solve the Navier-Stokes equations for U(n+1), u(n+1),
        // P(n+1/2).  For better or worse, each of the major
        // algorithmic steps of the projection method is separated
        // into its own member function.
        d_ins_hier_integrator->predictAdvectionVelocity(
            current_time, new_time);

        d_ins_hier_integrator->integrateAdvDiff(current_time, new_time);

        d_ins_hier_integrator->projectVelocity(current_time, new_time);

        d_ins_hier_integrator->updatePressure(current_time, new_time, true);

        d_ins_hier_integrator->synchronizeHierarchy();

        if (cycle < (num_cycles-1))
        {
            // Reset data to the preadvance state.
            d_ins_hier_integrator->resetHierDataToPreadvanceState();

            // Reset time of all current data.
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                level->setTime(current_time);
            }
        }
        else
        {
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                level->setTime(new_time);
            }
        }
    }

    // Compute the final structure configuration, X(n+1).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            // XXXX
#if 0
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;

            refine_operator = grid_geom->lookupRefineOperator(
                d_ins_hier_integrator->getVelocityVar(),
                "CONSERVATIVE_LINEAR_REFINE");

            int U_new_idx = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase()->
                mapVariableAndContextToIndex(
                    d_ins_hier_integrator->getVelocityVar(),
                    d_ins_hier_integrator->getNewContext());

            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > bdry_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
            bdry_fill_alg->registerRefine(d_W_idx,   // destination
                                          U_new_idx, // source
                                          d_W_idx,   // temporary work space
                                          refine_operator);

            bdry_fill_alg->resetSchedule(d_rscheds["U_current->V::CONSERVATIVE_LINEAR_REFINE"][ln]);
            d_rscheds["U_current->V::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(new_time);
            d_ralgs["U_current->V::CONSERVATIVE_LINEAR_REFINE"]->resetSchedule(
                d_rscheds["U_current->V::CONSERVATIVE_LINEAR_REFINE"][ln]);
#endif
            // XXXX

            d_rscheds["U_new->W::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(new_time);

            // Interpolate U(X~(n+1),n+1) = u(X~(n+1),n+1).
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): interpolating u(n+1) to U(n+1)\n";

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > w_data =
                    patch->getPatchData(d_W_idx);
                const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data = patch->getPatchData(
                    d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());

                LEInteractor::interpolate(
                    U_data[ln], X_new_data[ln], idx_data, w_data,
                    patch, patch_box,
                    d_delta_fcn, ENFORCE_PERIODIC_BCS);
            }

            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing X(n+1)\n";

            switch (d_timestepping_order)
            {
                case 2:
                {
                    // Set X(n+1) = 0.5[X(n) + X~(n+1)] + dt*U(X~(n+1),n+1)).
                    int ierr;
                    ierr = VecAXPY(X_new_data[ln]->getGlobalVec(),
                                   dt, U_data[ln]->getGlobalVec());
                    PETSC_SAMRAI_ERROR(ierr);

                    ierr = VecAXPBY(X_data[ln]->getGlobalVec(),
                                    0.5, 0.5, X_new_data[ln]->getGlobalVec());
                    PETSC_SAMRAI_ERROR(ierr);

                    break;
                }
                case 3:
                {
                    // Set X~(n+1/2) = (3/4)*X(n) + (1/4)*X~(n+2)
                    int ierr;
                    ierr = VecAXPY(X_new_data[ln]->getGlobalVec(),
                                   dt, U_data[ln]->getGlobalVec());
                    PETSC_SAMRAI_ERROR(ierr);

                    ierr = VecAXPBY(X_new_data[ln]->getGlobalVec(),
                                    0.75, 0.25, X_data[ln]->getGlobalVec());
                    PETSC_SAMRAI_ERROR(ierr);

                    // Interpolate the Cartesian velocity onto the Lagrangian
                    // mesh, U(X~(n+1/2),n+1/2) = u(X~(n+1/2),n+1/2).
                    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > v_data =
                            patch->getPatchData(d_V_idx);
                        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > w_data =
                            patch->getPatchData(d_W_idx);
                        const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data = patch->getPatchData(
                            d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());

                        SAMRAI::math::PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;

                        patch_cc_data_ops.linearSum(
                            w_data, 0.5, v_data, 0.5, w_data,
                            v_data->getGhostBox() * w_data->getGhostBox());

                        LEInteractor::interpolate(
                            U_data[ln], X_new_data[ln], idx_data, w_data,
                            patch, patch_box,
                            d_delta_fcn, ENFORCE_PERIODIC_BCS);
                    }

                    // Set X~(n+1) = (1/3)*X(n) + (2/3)*X~(n+3/2)
                    ierr = VecAXPY(X_new_data[ln]->getGlobalVec(),
                                   dt, U_data[ln]->getGlobalVec());
                    PETSC_SAMRAI_ERROR(ierr);

                    ierr = VecAXPBY(X_data[ln]->getGlobalVec(),
                                    2.0/3.0, 1.0/3.0, X_new_data[ln]->getGlobalVec());
                    PETSC_SAMRAI_ERROR(ierr);
                    break;
                }
                default:
                {
                    TBOX_ERROR("Unsupported order for the Lagrangian timestepping scheme, " << d_timestepping_order << endl);
                    break;
                }
            }
        }
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_V_idx);
        level->deallocatePatchData(d_W_idx);
        level->deallocatePatchData(d_F_idx);
        if (!d_source_strategy.isNull())
        {
            level->deallocatePatchData(d_P_idx);
            level->deallocatePatchData(d_Q_idx);
        }
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
IBHierarchyIntegrator::atRegridPoint() const
{
    const int level_number = 0;

    return ((d_integrator_step>0)
            && d_gridding_alg->levelCanBeRefined(level_number)
            && (d_regrid_interval == 0
                ? false
                : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double
IBHierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}// getIntegratorTime

double
IBHierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}// getStartTime

double
IBHierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}// getEndTime

int
IBHierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}// getIntegratorStep

int
IBHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}// getMaxIntegratorSteps

bool
IBHierarchyIntegrator::stepsRemaining() const
{
    return d_integrator_step < d_max_integrator_steps;
}// stepsRemaining

const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
IBHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
IBHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

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
///  allow the IBHierarchyIntegrator to provide data management
///  for a time integrator which making use of this class.
///

void
IBHierarchyIntegrator::rebalanceCoarsestLevel()
{
    t_rebalance_coarsest_level->start();

    // Update the workload for regridding.
    d_lag_data_manager->updateWorkloadAndNodeCount(
        0,d_hierarchy->getFinestLevelNumber());

    // Before rebalancing, begin Lagrangian data movement.
    d_lag_data_manager->beginDataRedistribution();

    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->rebalanceCoarsestLevel();

    // After rebalancing, finish Lagrangian data movement.
    d_lag_data_manager->endDataRedistribution();

    // Update the workload post-regridding.
    d_lag_data_manager->updateWorkloadAndNodeCount(
        0,d_hierarchy->getFinestLevelNumber());

    // Indicate that the force and source strategies need to be
    // re-initialized.
    d_force_strategy_needs_init  = true;
    d_source_strategy_needs_init = true;

    t_rebalance_coarsest_level->stop();
    return;
}// rebalanceCoarsestLevel

void
IBHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    // Update the workload pre-regridding.
    d_lag_data_manager->updateWorkloadAndNodeCount(
        0,d_hierarchy->getFinestLevelNumber());

    // Before regriding, begin Lagrangian data movement.
    d_lag_data_manager->beginDataRedistribution();

    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->regridHierarchy();

    // After regridding, finish Lagrangian data movement.
    d_lag_data_manager->endDataRedistribution();

    // Update the workload post-regridding.
    d_lag_data_manager->updateWorkloadAndNodeCount(
        0,d_hierarchy->getFinestLevelNumber());

    // Indicate that the force and source strategies need to be
    // re-initialized.
    d_force_strategy_needs_init  = true;
    d_source_strategy_needs_init = true;

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

void
IBHierarchyIntegrator::synchronizeHierarchy()
{
    t_synchronize_hierarchy->start();

    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->synchronizeHierarchy();

    t_synchronize_hierarchy->stop();
    return;
}// synchronizeHierarchy

void
IBHierarchyIntegrator::synchronizeNewLevels(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
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
    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->
        synchronizeNewLevels(hierarchy, coarsest_level, finest_level,
                             sync_time, initial_time);

    t_synchronize_new_levels->stop();
    return;
}// synchronizeNewLevels

void
IBHierarchyIntegrator::resetTimeDependentHierData(
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
IBHierarchyIntegrator::resetHierDataToPreadvanceState()
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
IBHierarchyIntegrator::initializeLevelData(
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
    assert(!hierarchy.isNull());
    assert((level_number >= 0)
           && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!(old_level.isNull())) {
        assert(level_number == old_level->getLevelNumber());
    }
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
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
IBHierarchyIntegrator::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

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

    // We use the INSHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->
        resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // We use the LDataManager to handle as much unstructured data
    // management as possible.
    d_lag_data_manager->
        resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // If we have added or removed a level, resize the schedule
    // vectors.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it!= d_ralgs.end(); ++it)
    {
        d_rscheds[(*it).first].resize(finest_hier_level+1);
    }

    // (Re)build refine communication schedules.  These are created
    // for all levels in the hierarchy.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it!= d_ralgs.end(); ++it)
    {
        for (int ln = coarsest_level; ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

            d_rscheds[(*it).first][ln] = (*it).second->
                createSchedule(level, ln-1, hierarchy);
        }
    }

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
IBHierarchyIntegrator::applyGradientDetector(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
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

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#if 0
    // Tag cells which contain point sources and sinks.
    if (!d_source_strategy.isNull())
    {
        // Get the present source locations.
        //
        // IMPORTANT NOTE: We assume that each MPI process is assigned
        // the source locations.
        static bool warned = false;
        if (!warned) TBOX_WARNING("STUPID STUPID STUPID TAGGING STRATEGY!\n");
        warned = true;

        const int nsrc = d_source_strategy->getNumSources();
        std::vector<std::vector<double> > X_src(nsrc, std::vector<double>(NDIM,0.0));
        d_source_strategy->getSourceLocations(
            X_src, NULL, // XXXX This will break the heart!
            d_hierarchy, level_number, error_data_time, d_lag_data_manager);

        bool tag_box_set = false;
        std::vector<SAMRAI::hier::Box<NDIM> > tag_box(nsrc);

        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);

            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
            const SAMRAI::hier::Index<NDIM>& patch_upper = patch_box.upper();

            const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const xLower = pgeom->getXLower();
            const double* const xUpper = pgeom->getXUpper();
            const double* const dx     = pgeom->getDx();

            if (!tag_box_set)
            {
                for (int n = 0; n < nsrc; ++n)
                {
                    const SAMRAI::hier::Index<NDIM> idx =
                        STOOLS::STOOLS_Utilities::getCellIndex(
                            &X_src[n][0],xLower,xUpper,dx,patch_lower,patch_upper);
                    tag_box[n] = SAMRAI::hier::Box<NDIM>(idx,idx);
#if 0
                    for (int d = 0; d < NDIM; ++d)
                    {
                        if (X_src[n][d] <= xLower[d]+(static_cast<double>(idx(d))+0.5)*dx[d])
                        {
                            tag_box[n].growLower(d,1);
                        }
                        else
                        {
                            tag_box[n].growUpper(d,1);
                        }
                    }
#endif
                }
                tag_box_set = true;
            }

            for (int n = 0; n < nsrc; ++n)
            {
                if (patch_box.intersects(tag_box[n]))
                {
                    tags_data->fillAll(1,patch_box*tag_box[n]);
                }
            }
        }
    }
#endif
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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
IBHierarchyIntegrator::getCurrentContext() const
{
    return d_ins_hier_integrator->getCurrentContext();
}// getCurrentContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBHierarchyIntegrator::getNewContext() const
{
    return d_ins_hier_integrator->getNewContext();
}// getNewContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBHierarchyIntegrator::getOldContext() const
{
    return d_ins_hier_integrator->getOldContext();
}// getOldContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBHierarchyIntegrator::getScratchContext() const
{
    return d_ins_hier_integrator->getScratchContext();
}// getScratchContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBHierarchyIntegrator::getPlotContext() const
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
IBHierarchyIntegrator::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif
    db->putInteger("IB_HIERARCHY_INTEGRATOR_VERSION",
                   IB_HIERARCHY_INTEGRATOR_VERSION);
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
IBHierarchyIntegrator::printClassData(
    std::ostream& os) const
{
    os << "\nIBHierarchyIntegrator::printClassData..." << endl;
    os << "\nIBHierarchyIntegrator: this = "
       << const_cast<IBHierarchyIntegrator*>(this) << endl;
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

void
IBHierarchyIntegrator::resetLagrangianForceStrategy(
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
IBHierarchyIntegrator::resetLagrangianSourceStrategy(
    const double init_data_time,
    const bool initial_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            d_source_strategy->initializeLevelData(
                d_hierarchy, ln, init_data_time, initial_time,
                d_lag_data_manager);
        }
    }

    return;
}// resetLagrangianSourceStrategy

void
IBHierarchyIntegrator::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif
    // Read data members from input database.
    d_delta_fcn = db->getStringWithDefault("delta_fcn", d_delta_fcn);

    d_end_time = db->getDoubleWithDefault("end_time", d_end_time);
    d_grow_dt = db->getDoubleWithDefault("grow_dt", d_grow_dt);

    d_max_integrator_steps = db->getIntegerWithDefault(
        "max_integrator_steps", d_max_integrator_steps);

    d_num_cycles = db->getIntegerWithDefault("num_cycles", d_num_cycles);

    d_regrid_interval = db->getIntegerWithDefault(
        "regrid_interval", d_regrid_interval);

    d_timestepping_order = db->getIntegerWithDefault(
        "timestepping_order", d_timestepping_order);

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault(
        "dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault(
        "dt_max_time_min", d_dt_max_time_min);

    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    if (!is_from_restart)
    {
        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);

        d_num_init_cycles = db->getIntegerWithDefault(
            "num_init_cycles", d_num_init_cycles);
    }

    return;
}// getFromInput

void
IBHierarchyIntegrator::getFromRestart()
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

    int ver = db->getInteger("IB_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_HIERARCHY_INTEGRATOR_VERSION)
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

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
