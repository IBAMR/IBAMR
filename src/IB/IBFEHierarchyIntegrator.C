// Filename: IBFEHierarchyIntegrator.C
// Last modified: <28.Apr.2010 11:47:48 griffith@boyce-griffiths-mac-pro.local>
// Created on 27 Jul 2009 by Boyce Griffith (griffith@griffith-macbook-pro.local)

#include "IBFEHierarchyIntegrator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// C++ STDLIB INCLUDES
#include <algorithm>
#include <iterator>
#include <limits>
#include <numeric>

// IBAMR INCLUDES
#include <ibamr/IBInstrumentationSpec.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>

// LIBMESH INCLUDES
#include <boundary_info.h>
#include <dense_vector.h>
#include <dof_map.h>
#include <fe.h>
#include <explicit_system.h>
#include <numeric_vector.h>
#include <petsc_vector.h>
#include <quadrature_trap.h>
#include <quadrature_gauss.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CoarsenOperator.h>
#include <HierarchyDataOpsManager.h>
#include <Index.h>
#include <IndexData.h>
#include <Patch.h>
#include <VariableDatabase.h>
#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

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
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_regrid_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_new_levels;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_time_dependent_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_data_to_preadvance_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_hierarchy_configuration;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_apply_gradient_detector;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;

// Version of IBFEHierarchyIntegrator restart file data.
static const int IB_FE_HIERARCHY_INTEGRATOR_VERSION = 1;
}

const std::string IBFEHierarchyIntegrator::       COORDINATES_SYSTEM_NAME = "IB coordinates system";
const std::string IBFEHierarchyIntegrator::COORDINATE_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string IBFEHierarchyIntegrator::             FORCE_SYSTEM_NAME = "IB force system";
const std::string IBFEHierarchyIntegrator::          VELOCITY_SYSTEM_NAME = "IB velocity system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEHierarchyIntegrator::IBFEHierarchyIntegrator(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    SAMRAI::tbox::Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
    IBTK::FEDataManager* fe_data_manager,
    IBTK::LDataManager* lag_data_manager,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_fe_data_manager(fe_data_manager),
      d_split_interior_and_bdry_forces(false),
      d_coordinate_mapping_function(NULL),
      d_coordinate_mapping_function_ctx(NULL),
      d_PK1_stress_function(NULL),
      d_PK1_stress_function_ctx(NULL),
      d_lag_data_manager(lag_data_manager),
      d_hierarchy(hierarchy),
      d_gridding_alg(NULL),
      d_visit_writer(NULL),
      d_load_balancer(NULL),
      d_ins_hier_integrator(ins_hier_integrator),
      d_body_force_fcn(NULL),
      d_eulerian_force_fcn(NULL),
      d_start_time(0.0),
      d_end_time(std::numeric_limits<double>::max()),
      d_grow_dt(2.0),
      d_max_integrator_steps(std::numeric_limits<int>::max()),
      d_num_cycles(1),
      d_regrid_interval(1),
      d_old_dt(-1.0),
      d_integrator_time(std::numeric_limits<double>::quiet_NaN()),
      d_integrator_step(std::numeric_limits<int>::max()),
      d_dt_max(std::numeric_limits<double>::max()),
      d_dt_max_time_max(std::numeric_limits<double>::max()),
      d_dt_max_time_min(-(d_dt_max_time_max-std::numeric_limits<double>::epsilon())),
      d_is_initialized(false),
      d_do_log(false),
      d_hier_sc_data_ops(),
      d_ralgs(),
      d_rstrategies(),
      d_rscheds(),
      d_calgs(),
      d_cstrategies(),
      d_cscheds(),
      d_V_var(NULL),
      d_F_var(NULL),
      d_current(NULL),
      d_scratch(NULL),
      d_V_idx(-1),
      d_F_idx(-1)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT(!ins_hier_integrator.isNull());
#endif

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager = SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > sc_var = new SAMRAI::pdat::SideVariable<NDIM,double>("sc_var");
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, hierarchy);

    // Initialize all variable contexts.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    d_current = var_db->getContext(d_object_name+"::CURRENT");
    d_scratch = var_db->getContext(d_object_name+"::SCRATCH");

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy            = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy               = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()");
        t_regrid_hierarchy                = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::regridHierarchy()");
        t_synchronize_hierarchy           = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels          = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data       = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state  = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data           = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration   = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector         = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database                 = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// IBFEHierarchyIntegrator

IBFEHierarchyIntegrator::~IBFEHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    for (RefinePatchStrategyMap::iterator it = d_rstrategies.begin(); it != d_rstrategies.end(); ++it)
    {
        delete (*it).second;
    }

    for (CoarsenPatchStrategyMap::iterator it = d_cstrategies.begin(); it != d_cstrategies.end(); ++it)
    {
        delete (*it).second;
    }
    return;
}// ~IBFEHierarchyIntegrator

const std::string&
IBFEHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
IBFEHierarchyIntegrator::setInitialCoordinateMappingFunction(
    Point (*coordinate_mapping_function)(const Point& s, void* ctx),
    void* coordinate_mapping_function_ctx)
{
    d_coordinate_mapping_function = coordinate_mapping_function;
    d_coordinate_mapping_function_ctx = coordinate_mapping_function_ctx;
    return;
}// setInitialCoordinateMappingFunction

void
IBFEHierarchyIntegrator::setPK1StressTensorFunction(
    TensorValue<double> (*PK1_stress_function)(const TensorValue<double>& dX_ds, const Point& X, const Point& s, Elem* const elem, const double& time, void* ctx),
    void* PK1_stress_function_ctx)
{
    d_PK1_stress_function = PK1_stress_function;
    d_PK1_stress_function_ctx = PK1_stress_function_ctx;
    return;
}// setPK1StressTensorFunction

void
IBFEHierarchyIntegrator::setLagrangianForceStrategy(
    SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> lag_force_strategy)
{
    d_lag_force_strategy = lag_force_strategy;
    d_lag_force_strategy_needs_init = true;
    return;
}// setLagrangianForceStrategy

void
IBFEHierarchyIntegrator::registerVelocityInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> U_init)
{
    d_ins_hier_integrator->registerVelocityInitialConditions(U_init);
    return;
}// registerVelocityInitialConditions

void
IBFEHierarchyIntegrator::registerVelocityPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::registerVelocityPhysicalBcCoefs()\n"
                   << "  velocity boundary conditions must be registered prior to initialization\n"
                   << "  of the hierarchy integrator object.\n");
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned l = 0; l < U_bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(U_bc_coefs[l] != NULL);
    }
#endif
    d_ins_hier_integrator->registerVelocityPhysicalBcCoefs(U_bc_coefs);
    return;
}// registerVelocityPhysicalBcCoefs

void
IBFEHierarchyIntegrator::registerBodyForceSpecification(
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> body_force_fcn)
{
    d_body_force_fcn = body_force_fcn;
    return;
}// registerBodyForceSpecification

void
IBFEHierarchyIntegrator::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    d_ins_hier_integrator->registerVisItDataWriter(d_visit_writer);
    return;
}// registerVisItDataWriter

void
IBFEHierarchyIntegrator::registerLoadBalancer(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!load_balancer.isNull());
#endif
    d_load_balancer = load_balancer;
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
///      getGriddingAlgorithm()
///
///  allow the IBFEHierarchyIntegrator to be used as a hierarchy integrator.
///

void
IBFEHierarchyIntegrator::initializeHierarchyIntegrator(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Initialize FE system data.
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();

    d_fe_data_manager->COORDINATES_SYSTEM_NAME = COORDINATES_SYSTEM_NAME;
    ExplicitSystem& coords_system = equation_systems->add_system<ExplicitSystem>(COORDINATES_SYSTEM_NAME);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        std::ostringstream os;
        os << "X_" << d;
        coords_system.add_variable(os.str(), FIRST, LAGRANGE);
    }

    ExplicitSystem& coords_mapping_system = equation_systems->add_system<ExplicitSystem>(COORDINATE_MAPPING_SYSTEM_NAME);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        std::ostringstream os;
        os << "dX_" << d;
        coords_mapping_system.add_variable(os.str(), FIRST, LAGRANGE);
    }

    ExplicitSystem& force_system = equation_systems->add_system<ExplicitSystem>(FORCE_SYSTEM_NAME);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        std::ostringstream os;
        os << "F_" << d;
        force_system.add_variable(os.str(), FIRST, LAGRANGE);
    }

    ExplicitSystem& velocity_system = equation_systems->add_system<ExplicitSystem>(VELOCITY_SYSTEM_NAME);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        std::ostringstream os;
        os << "U_" << d;
        velocity_system.add_variable(os.str(), FIRST, LAGRANGE);
    }

    // Initialize all variables.
    const SAMRAI::hier::IntVector<NDIM> ghosts = d_fe_data_manager->getGhostCellWidth();
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    d_V_var = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::V");
    d_V_idx = var_db->registerVariableAndContext(d_V_var, d_scratch, ghosts);

    d_F_var = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::F");
    d_F_idx = var_db->registerVariableAndContext(d_F_var, d_scratch, no_ghosts);

    // Initialize the objects used to manage Lagragian-Eulerian interaction.
    d_eulerian_force_fcn = new IBEulerianForceFunction(d_object_name+"::IBEulerianForceFunction", -1, -1, d_F_idx);
    d_ins_hier_integrator->registerBodyForceSpecification(d_eulerian_force_fcn);

    // Initialize the INSStaggeredHierarchyIntegrator.
    //
    // NOTE: This must be done after all variables created by the
    // IBFEHierarchyIntegrator class have been registered.
    d_ins_hier_integrator->initializeHierarchyIntegrator(d_gridding_alg);

    // Create several communications algorithms, used in filling ghost cell data
    // and synchronizing data on the patch hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_operator;

    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getCurrentContext());

    d_ralgs["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getVelocityVar(), "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"]->registerRefine(d_V_idx, U_current_idx, d_V_idx, refine_operator);
    d_rstrategies["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"] = new IBTK::CartExtrapPhysBdryOp(d_V_idx, "QUADRATIC");

    d_ralgs["V->V::S->S::CONSERVATIVE_LINEAR_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getVelocityVar(), "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["V->V::S->S::CONSERVATIVE_LINEAR_REFINE"]->registerRefine(d_V_idx, d_V_idx, d_V_idx, refine_operator);
    d_rstrategies["V->V::S->S::CONSERVATIVE_LINEAR_REFINE"] = new IBTK::CartExtrapPhysBdryOp(d_V_idx, "QUADRATIC");

    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(d_V_var, "CONSERVATIVE_COARSEN");
    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"]->registerCoarsen(U_current_idx, U_current_idx, coarsen_operator);

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
IBFEHierarchyIntegrator::initializeHierarchy()
{
    t_initialize_hierarchy->start();

    // Initialize the data structures and state variables for the FE equation
    // systems.
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    equation_systems->init();
    initializeCoordinates();
    updateCoordinateMapping();

    System& coords_system = equation_systems->get_system<System>(COORDINATES_SYSTEM_NAME);
    coords_system.assemble_before_solve = false;
    coords_system.assemble();

    System& coords_mapping_system = equation_systems->get_system<System>(COORDINATE_MAPPING_SYSTEM_NAME);
    coords_mapping_system.assemble_before_solve = false;
    coords_mapping_system.assemble();

    // Use the INSStaggeredHierarchyIntegrator to initialize the patch
    // hierarchy.
    double dt_next = d_ins_hier_integrator->initializeHierarchy();
    if (d_integrator_time >= d_dt_max_time_min && d_integrator_time <= d_dt_max_time_max)
    {
        dt_next = std::min(dt_next, d_dt_max);
    }

    // Initialize the FE data manager.
    d_fe_data_manager->reinitElementMappings();

    // Reset the Lagrangian data manager.
    if (d_lag_data_manager != NULL)
    {
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        d_lag_data_manager->beginDataRedistribution();
        d_lag_data_manager->endDataRedistribution();
        d_lag_data_manager->updateWorkloadData(coarsest_ln, finest_ln);
    }

    // Indicate that the force strategy needs to be re-initialized.
    d_lag_force_strategy_needs_init = true;

    t_initialize_hierarchy->stop();
    return dt_next;
}// initializeHierarchy

double
IBFEHierarchyIntegrator::advanceHierarchy(
    const double dt)
{
    t_advance_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
#endif

    PetscErrorCode ierr;

    const double current_time = d_integrator_time;
    const double new_time     = d_integrator_time+dt;
    const bool initial_time   = SAMRAI::tbox::MathUtilities<double>::equalEps(current_time,d_start_time);

    // Regrid the patch hierarchy.
    const bool do_regrid = (d_regrid_interval == 0 ? false : (d_integrator_step % d_regrid_interval == 0));
    if (do_regrid)
    {
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): regridding prior to timestep " << d_integrator_step << "\n";
        regridHierarchy();
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Set the current time interval in the force specification objects.
    d_eulerian_force_fcn->registerBodyForceSpecification(d_body_force_fcn);
    d_eulerian_force_fcn->setTimeInterval(current_time, new_time);
    if (!d_lag_force_strategy.isNull()) d_lag_force_strategy->setTimeInterval(current_time, new_time);

    // (Re)initialize the force strategy object.
    if (d_lag_force_strategy_needs_init)
    {
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            if (d_lag_data_manager != NULL && d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                d_lag_force_strategy->initializeLevelData(d_hierarchy, ln, current_time, initial_time, d_lag_data_manager);
            }
        }
        d_lag_force_strategy_needs_init = false;
    }

    // Extract the FE vectors.
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    System& coords_system = equation_systems->get_system<System>(COORDINATES_SYSTEM_NAME);
    System& force_system = equation_systems->get_system<System>(FORCE_SYSTEM_NAME);
    System& velocity_system = equation_systems->get_system<System>(VELOCITY_SYSTEM_NAME);

    NumericVector<double>& X_current = *(coords_system.current_local_solution);
    coords_system.get_dof_map().enforce_constraints_exactly(coords_system, &X_current);

    AutoPtr<NumericVector<double> > X_new_ptr = X_current.clone();
    NumericVector<double>& X_new = *X_new_ptr;

    AutoPtr<NumericVector<double> > X_half_ptr = X_current.clone();
    NumericVector<double>& X_half = *X_half_ptr;

    NumericVector<double>& F_half = *(   force_system.current_local_solution);
    NumericVector<double>& U_half = *(velocity_system.current_local_solution);

    NumericVector<double>* X_half_IB_ghost_ptr = d_fe_data_manager->getGhostedCoordsVector();
    NumericVector<double>& X_half_IB_ghost = *X_half_IB_ghost_ptr;

    NumericVector<double>* F_half_IB_ghost_ptr = d_fe_data_manager->getGhostedSolutionVector(FORCE_SYSTEM_NAME);
    NumericVector<double>& F_half_IB_ghost = *F_half_IB_ghost_ptr;

    // Initialize the various LNodeLevelData objects on each level of the patch hierarchy.
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_current_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_new_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_half_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > U_half_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_half_data(finest_ln+1);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager != NULL && d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_current_data[ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,ln);
            X_new_data[ln] = d_lag_data_manager->createLNodeLevelData("X_new",ln,NDIM);
            X_half_data[ln] = d_lag_data_manager->createLNodeLevelData("X_half",ln,NDIM);
            U_half_data[ln] = d_lag_data_manager->createLNodeLevelData("U_half",ln,NDIM);
            F_half_data[ln] = d_lag_data_manager->createLNodeLevelData("F_half",ln,NDIM);

            X_current_data[ln]->restoreLocalFormVec();
            X_new_data[ln]->restoreLocalFormVec();
            X_half_data[ln]->restoreLocalFormVec();
            U_half_data[ln]->restoreLocalFormVec();
            F_half_data[ln]->restoreLocalFormVec();
        }
    }

    // Get patch data descriptors for the current and new velocity data.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getCurrentContext());
    const int U_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getNewContext());

    // Allocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_V_idx, current_time);
        level->allocatePatchData(d_F_idx, current_time);
    }

    // Synchronize the Cartesian grid velocity u(n) on the patch hierarchy.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["U->U::C->C::CONSERVATIVE_COARSEN"][ln]->coarsenData();
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_rscheds["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(current_time);
    }

    // Initialize X(n+1) to equal X(n).
    ierr = VecCopy(dynamic_cast<PetscVector<double>*>(&X_current)->vec(),
                   dynamic_cast<PetscVector<double>*>(&X_new)->vec()); IBTK_CHKERRQ(ierr);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager != NULL && d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_current_vec = X_current_data[ln]->getGlobalVec();
            Vec X_new_vec = X_new_data[ln]->getGlobalVec();
            ierr = VecCopy(X_current_vec, X_new_vec); IBTK_CHKERRQ(ierr);
        }
    }

    // Perform one or more cycles to compute the updated configuration of the
    // coupled fluid-structure system.
    d_ins_hier_integrator->integrateHierarchy_initialize(current_time, new_time);
    for (int cycle = 0; cycle < d_num_cycles; ++cycle)
    {
        // Set X(n+1/2) = 0.5*(X(n) + X(n+1)).
        ierr = VecAXPBYPCZ(dynamic_cast<PetscVector<double>*>(&X_half)->vec(),
                           0.5, 0.5, 0.0,
                           dynamic_cast<PetscVector<double>*>(&X_current)->vec(),
                           dynamic_cast<PetscVector<double>*>(&X_new)->vec()); IBTK_CHKERRQ(ierr);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_lag_data_manager != NULL && d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                Vec X_current_vec = X_current_data[ln]->getGlobalVec();
                Vec X_new_vec = X_new_data[ln]->getGlobalVec();
                Vec X_half_vec = X_half_data[ln]->getGlobalVec();
                ierr = VecAXPBYPCZ(X_half_vec, 0.5, 0.5, 0.0, X_current_vec, X_new_vec); IBTK_CHKERRQ(ierr);
            }
        }

        // Compute F(n+1/2) = F(X(n+1/2),t(n+1/2)).
        computeInteriorForceDensity(F_half, X_half, current_time+0.5*dt);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_lag_data_manager != NULL && d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                Vec F_half_vec = F_half_data[ln]->getGlobalVec();
                ierr = VecSet(F_half_vec, 0.0); IBTK_CHKERRQ(ierr);
                d_lag_force_strategy->computeLagrangianForce(F_half_data[ln], X_half_data[ln], U_half_data[ln], d_hierarchy, ln, current_time+0.5*dt, d_lag_data_manager);
            }
        }

        // Copy data into the "IB ghosted" vectors.
        ierr = VecCopy(dynamic_cast<PetscVector<double>*>(&X_half)->vec(), dynamic_cast<PetscVector<double>*>(&X_half_IB_ghost)->vec()); IBTK_CHKERRQ(ierr);
        ierr = VecCopy(dynamic_cast<PetscVector<double>*>(&F_half)->vec(), dynamic_cast<PetscVector<double>*>(&F_half_IB_ghost)->vec()); IBTK_CHKERRQ(ierr);

        // Spread F(n+1/2) to f(n+1/2).
        d_hier_sc_data_ops->setToScalar(d_F_idx, 0.0);
        d_fe_data_manager->spread(d_F_idx, F_half_IB_ghost, X_half_IB_ghost, FORCE_SYSTEM_NAME);
        if (d_split_interior_and_bdry_forces)
        {
            spreadBoundaryForceDensity(d_F_idx, X_half_IB_ghost, current_time+0.5*dt);
        }
        if (d_lag_data_manager != NULL)
        {
            d_lag_data_manager->spread(d_F_idx, F_half_data, X_half_data);
        }

        // Solve the incompressible Navier-Stokes equations.
        d_ins_hier_integrator->integrateHierarchy(current_time, new_time);

        // Set u(n+1/2) = 0.5*(u(n) + u(n+1)).
        d_hier_sc_data_ops->linearSum(d_V_idx, 0.5, U_current_idx, 0.5, U_new_idx);

        // Interpolate u(n+1/2) to U(n+1/2).
        d_fe_data_manager->interpolate(d_V_idx, U_half, X_half_IB_ghost, VELOCITY_SYSTEM_NAME, d_rscheds["V->V::S->S::CONSERVATIVE_LINEAR_REFINE"], current_time);
        if (d_lag_data_manager != NULL)
        {
            d_lag_data_manager->interpolate(d_V_idx, U_half_data, X_half_data, std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(), current_time, false);
        }

        // Set X(n+1) = X(n) + dt*U(n+1/2).
        ierr = VecWAXPY(dynamic_cast<PetscVector<double>*>(&X_new)->vec(),
                        dt,
                        dynamic_cast<PetscVector<double>*>(&U_half)->vec(),
                        dynamic_cast<PetscVector<double>*>(&X_current)->vec()); IBTK_CHKERRQ(ierr);
        coords_system.get_dof_map().enforce_constraints_exactly(coords_system, &X_new);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_lag_data_manager != NULL && d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                Vec X_current_vec = X_current_data[ln]->getGlobalVec();
                Vec X_new_vec = X_new_data[ln]->getGlobalVec();
                Vec U_half_vec = U_half_data[ln]->getGlobalVec();
                ierr = VecWAXPY(X_new_vec, dt, U_half_vec, X_current_vec); IBTK_CHKERRQ(ierr);
            }
        }
    }
    d_ins_hier_integrator->integrateHierarchy_finalize(current_time, new_time);

    // Reset X_current to equal X_new.
    ierr = VecCopy(dynamic_cast<PetscVector<double>*>(&X_new)->vec(),
                   dynamic_cast<PetscVector<double>*>(&X_current)->vec()); IBTK_CHKERRQ(ierr);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager != NULL && d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_current_vec = X_current_data[ln]->getGlobalVec();
            Vec X_new_vec = X_new_data[ln]->getGlobalVec();
            ierr = VecCopy(X_new_vec, X_current_vec); IBTK_CHKERRQ(ierr);
        }
    }

    // Copy the ghosted data into the non-ghosted solution vectors.
    ierr = VecCopy(dynamic_cast<PetscVector<double>*>(  coords_system.current_local_solution.get())->vec(),
                   dynamic_cast<PetscVector<double>*>(  coords_system.              solution.get())->vec()); IBTK_CHKERRQ(ierr);
    ierr = VecCopy(dynamic_cast<PetscVector<double>*>(   force_system.current_local_solution.get())->vec(),
                   dynamic_cast<PetscVector<double>*>(   force_system.              solution.get())->vec()); IBTK_CHKERRQ(ierr);
    ierr = VecCopy(dynamic_cast<PetscVector<double>*>(velocity_system.current_local_solution.get())->vec(),
                   dynamic_cast<PetscVector<double>*>(velocity_system.              solution.get())->vec()); IBTK_CHKERRQ(ierr);

    // Update the coordinate mapping dX = X - s.
    updateCoordinateMapping();

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_V_idx);
        level->deallocatePatchData(d_F_idx);
    }

    // Synchronize all data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): synchronizing data\n";
    synchronizeHierarchy();

    // Reset all time dependent data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierData(new_time);

    // Determine the next stable timestep.
    double dt_next = d_ins_hier_integrator->getStableTimestep(getCurrentContext());

    if (d_integrator_time+dt_next >= d_end_time)
    {
        dt_next = d_end_time - d_integrator_time;
    }

    if (d_integrator_time >= d_dt_max_time_min && d_integrator_time <= d_dt_max_time_max)
    {
        dt_next = std::min(d_dt_max,dt_next);
    }

    if (!initial_time)
    {
        dt_next = std::min(dt_next,d_grow_dt*d_old_dt);
    }

    t_advance_hierarchy->stop();
    return dt_next;
}// advanceHierarchy

bool
IBFEHierarchyIntegrator::atRegridPoint() const
{
    static const int level_number = 0;
    return ((d_integrator_step>0)
            && d_gridding_alg->levelCanBeRefined(level_number)
            && (d_regrid_interval == 0
                ? false
                : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double
IBFEHierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}// getIntegratorTime

double
IBFEHierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}// getStartTime

double
IBFEHierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}// getEndTime

int
IBFEHierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}// getIntegratorStep

int
IBFEHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}// getMaxIntegratorSteps

bool
IBFEHierarchyIntegrator::stepsRemaining() const
{
    return d_integrator_step < d_max_integrator_steps;
}// stepsRemaining

const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
IBFEHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
IBFEHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

///
///  The following routines:
///
///      regridHierarchy(),
///      synchronizeHierarchy(),
///      synchronizeNewLevels(),
///      resetTimeDependentHierData(),
///      resetHierDataToPreadvanceState()
///
///  allow the IBFEHierarchyIntegrator to provide data management for a time
///  integrator which making use of this class.
///

void
IBFEHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    // Begin Lagrangian data redistribution.
    if (d_lag_data_manager != NULL)
    {
        d_lag_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());
        d_lag_data_manager->beginDataRedistribution();
    }

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->regridHierarchy();

    // Complete Lagrangian data redistribution.
    if (d_lag_data_manager != NULL)
    {
        d_lag_data_manager->endDataRedistribution();
        d_lag_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());
    }
    d_lag_force_strategy_needs_init = true;

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

void
IBFEHierarchyIntegrator::synchronizeHierarchy()
{
    t_synchronize_hierarchy->start();

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->synchronizeHierarchy();

    t_synchronize_hierarchy->stop();
    return;
}// synchronizeHierarchy

void
IBFEHierarchyIntegrator::synchronizeNewLevels(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const bool initial_time)
{
    t_synchronize_new_levels->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0) && (coarsest_level < finest_level) && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->synchronizeNewLevels(hierarchy, coarsest_level, finest_level, sync_time, initial_time);

    t_synchronize_new_levels->stop();
    return;
}// synchronizeNewLevels

void
IBFEHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    t_reset_time_dependent_data->start();

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->resetTimeDependentHierData(new_time);

    t_reset_time_dependent_data->stop();
    return;
}// resetTimeDependentHierData

void
IBFEHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    t_reset_data_to_preadvance_state->start();

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
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
IBFEHierarchyIntegrator::initializeLevelData(
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
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    // Initialize the FE data manager.
    d_fe_data_manager->setPatchHierarchy(hierarchy);
    d_fe_data_manager->resetLevels(0,hierarchy->getFinestLevelNumber());
    d_fe_data_manager->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    // Initialize the Lagrangian data manager.
    if (d_lag_data_manager != NULL)
    {
        d_lag_data_manager->setPatchHierarchy(hierarchy);
        d_lag_data_manager->resetLevels(0,hierarchy->getFinestLevelNumber());
        d_lag_data_manager->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    }

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBFEHierarchyIntegrator::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0) && (coarsest_level <= finest_level) && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset the FE data manager.
    d_fe_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_hier_level);

    // Reset the Lagrangian data manager.
    if (d_lag_data_manager != NULL)
    {
        d_lag_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_hier_level);
    }

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->resetLevels(0, finest_hier_level);

    // If we have added or removed a level, resize the schedule vectors.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin(); it != d_ralgs.end(); ++it)
    {
        d_rscheds[(*it).first].resize(finest_hier_level+1);
    }

    for (CoarsenAlgMap::const_iterator it = d_calgs.begin(); it != d_calgs.end(); ++it)
    {
        d_cscheds[(*it).first].resize(finest_hier_level+1);
    }

    // (Re)build generic refine communication schedules.  These are created for
    // all levels in the hierarchy.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin(); it != d_ralgs.end(); ++it)
    {
        for (int ln = coarsest_level; ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            d_rscheds[(*it).first][ln] = (*it).second->createSchedule(level, ln-1, hierarchy, d_rstrategies[(*it).first]);
        }
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin(); it!= d_calgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(ln-1);
            d_cscheds[(*it).first][ln] = (*it).second->createSchedule(coarser_level, level, d_cstrategies[(*it).first]);
        }
    }

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
IBFEHierarchyIntegrator::applyGradientDetector(
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
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
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

    // Tag cells for refinement according to the criteria specified by the
    // INSStaggeredHierarchyIntegrator.
    d_ins_hier_integrator->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);

    // Tag cells which contain Lagragian nodes.
    d_fe_data_manager->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    if (d_lag_data_manager != NULL)
    {
        d_lag_data_manager->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }

    t_apply_gradient_detector->stop();
    return;
}// applyGradientDetector

///
///  The following routines:
///
///      getCurrentContext(),
///      getNewContext(),
///      getScratchContext()
///
///  allow access to the various variable contexts maintained by the integrator.
///

///
/// We simply reuse the SAMRAI::hier::VariableContext objects defined in the
/// INSStaggeredHierarchyIntegrator object.
///

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBFEHierarchyIntegrator::getCurrentContext() const
{
    return d_ins_hier_integrator->getCurrentContext();
}// getCurrentContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBFEHierarchyIntegrator::getNewContext() const
{
    return d_ins_hier_integrator->getNewContext();
}// getNewContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBFEHierarchyIntegrator::getScratchContext() const
{
    return d_ins_hier_integrator->getScratchContext();
}// getScratchContext

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::tbox::Serializable abstract base class.
///

void
IBFEHierarchyIntegrator::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("IB_FE_HIERARCHY_INTEGRATOR_VERSION", IB_FE_HIERARCHY_INTEGRATOR_VERSION);

    db->putBool("d_split_interior_and_bdry_forces", d_split_interior_and_bdry_forces);
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

    const std::vector<std::string> instrument_names = IBInstrumentationSpec::getInstrumentNames();
    if (!instrument_names.empty())
    {
        db->putInteger("instrument_names_sz", instrument_names.size());
        db->putStringArray("instrument_names", &instrument_names[0], instrument_names.size());
    }

    t_put_to_database->stop();
    return;
}// putToDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBFEHierarchyIntegrator::computeInteriorForceDensity(
    NumericVector<double>& G,
    NumericVector<double>& X,
    const double& time)
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();

    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    QGauss qrule(dim, FIFTH);
    QGauss qrule_face(dim-1, FIFTH);

    System& system = equation_systems->get_system<System>(FORCE_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dof_map.variable_type(0) == dof_map.variable_type(d));
    }
#endif

    AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    fe->attach_quadrature_rule(&qrule);
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi = fe->get_dphi();

    AutoPtr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
    fe_face->attach_quadrature_rule(&qrule_face);
    const std::vector<double>& JxW_face = fe_face->get_JxW();
    const std::vector<std::vector<double> >& phi_face = fe_face->get_phi();
    const std::vector<Point>& normal_face = fe_face->get_normals();

    System& coords_system = equation_systems->get_system<System>(COORDINATES_SYSTEM_NAME);
    const DofMap& coords_dof_map = coords_system.get_dof_map();
    std::vector<std::vector<unsigned int> > coords_dof_indices(NDIM);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(coords_dof_map.variable_type(0) == coords_dof_map.variable_type(d));
        TBOX_ASSERT(coords_dof_map.variable_type(d) ==        dof_map.variable_type(d));
    }
#endif

    AutoPtr<FEBase> coords_fe(FEBase::build(dim, coords_dof_map.variable_type(0)));
    coords_fe->attach_quadrature_rule(&qrule);
    const std::vector<Point>& coords_q_point = coords_fe->get_xyz();
    const std::vector<std::vector<double> >& coords_phi = coords_fe->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& coords_dphi = coords_fe->get_dphi();

    AutoPtr<FEBase> coords_fe_face(FEBase::build(dim, coords_dof_map.variable_type(0)));
    coords_fe_face->attach_quadrature_rule(&qrule_face);
    const std::vector<Point>& coords_q_point_face = coords_fe_face->get_xyz();
    const std::vector<std::vector<double> >& coords_phi_face = coords_fe_face->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& coords_dphi_face = coords_fe_face->get_dphi();

    // Loop over the elements to accumulate the interior forces at the nodes of
    // the mesh.  These are computed via
    //
    //    F_k = -int{P(s,t) grad phi_k(s)}ds + int{P(s,t) N(s,t) phi_k(s)}dA(s)
    //
    // This right-hand side vector is used to solve for the nodal values of the
    // interior elastic force density.
    AutoPtr<NumericVector<double> > F = G.clone();
    F->zero();
    std::vector<DenseVector<double> > F_e(NDIM);
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el = mesh.active_local_elements_begin(); el != end_el; ++el)
    {
        Elem* const elem = *el;

        fe->reinit(elem);
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            dof_map.dof_indices(elem, dof_indices[i], i);
            F_e[i].resize(dof_indices[i].size());
        }

        coords_fe->reinit(elem);
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            coords_dof_map.dof_indices(elem, coords_dof_indices[i], i);
        }

        // Loop over interior quadrature points.
        for (unsigned int qp = 0; qp < qrule.n_points(); ++qp)
        {
            // Compute the value of the first Piola-Kirchoff stress tensor at
            // the quadrature point.
            const Point& s_qp = coords_q_point[qp];
            const Point& X_qp = IBTK::compute_coordinate(qp,X,coords_phi,coords_dof_indices);
            const TensorValue<double> dX_ds = IBTK::compute_coordinate_mapping_jacobian(qp,X,coords_dphi,coords_dof_indices);
            const TensorValue<double> P = d_PK1_stress_function(dX_ds,X_qp,s_qp,elem,time,d_PK1_stress_function_ctx);

            // Accumulate the nodal forces.
            for (unsigned int k = 0; k < phi.size(); ++k)
            {
                const VectorValue<double> F_qp = -P*(dphi[k][qp])*JxW[qp];
                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    F_e[i](k) += F_qp(i);
                }
            }
        }

        if (d_split_interior_and_bdry_forces)
        {
            // Loop over the element boundaries.
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries and physical boundaries with
                // constraints.
                const short int boundary_id = mesh.boundary_info->boundary_id(elem,side);
                const bool at_physical_bdry = elem->neighbor(side) == NULL && !dof_map.is_periodic_boundary(boundary_id);
                const bool normal_dirichlet_bdry = boundary_id == NORMAL_DIRICHLET_BOUNDARY_ID || boundary_id == NORMAL_AND_TANGENTIAL_DIRICHLET_BOUNDARY_ID;
                if (at_physical_bdry && !normal_dirichlet_bdry)
                {
                    fe_face->reinit(elem, side);
                    coords_fe_face->reinit(elem, side);
                    for (unsigned int qp = 0; qp < qrule_face.n_points(); ++qp)
                    {
                        // Compute the value of the first Piola-Kirchoff stress
                        // tensor at the quadrature point.
                        const Point& s_qp = coords_q_point_face[qp];
                        const Point& X_qp = IBTK::compute_coordinate(qp,X,coords_phi_face,coords_dof_indices);
                        const TensorValue<double> dX_ds = IBTK::compute_coordinate_mapping_jacobian(qp,X,coords_dphi_face,coords_dof_indices);
                        const TensorValue<double> P = d_PK1_stress_function(dX_ds,X_qp,s_qp,elem,time,d_PK1_stress_function_ctx);
                        const VectorValue<double> F = P*normal_face[qp]*JxW_face[qp];

                        const TensorValue<double> dX_ds_inv_trans = IBTK::tensor_inverse_transpose(dX_ds, NDIM);
                        VectorValue<double> n = dX_ds_inv_trans*normal_face[qp];
                        n /= n.size();
                        const VectorValue<double> F_n = (F*n)*n;

                        // Split off the normal component of the force.
                        for (unsigned int k = 0; k < phi_face.size(); ++k)
                        {
                            VectorValue<double> F_qp;
                            F_qp += F_n*phi_face[k][qp];
                            for (unsigned int i = 0; i < NDIM; ++i)
                            {
                                F_e[i](k) += F_qp(i);
                            }
                        }
                    }
                }
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions) and
        // add the elemental contributions to the global vector.
        for (int i = 0; i < NDIM; ++i)
        {
            dof_map.constrain_element_vector(F_e[i], dof_indices[i]);
            F->add_vector(F_e[i], dof_indices[i]);
        }
    }

    // Assemble the right hand side vector.
    F->close();
    dof_map.enforce_constraints_exactly(system, F.get());

    // Solve for G, the nodal interior elastic force density.
    std::pair<LinearSolver<double>*,SparseMatrix<double>*> proj_solver_components = d_fe_data_manager->getL2ProjectionSolver(FORCE_SYSTEM_NAME);
    LinearSolver<double>* solver = proj_solver_components.first;
    SparseMatrix<double>* M = proj_solver_components.second;
    const double tol = 1.0e-10;
    const unsigned int max_its = 100;
    solver->solve(*M, *M, G, *F, tol, max_its);
    dof_map.enforce_constraints_exactly(system, &G);
    return;
}// computeInteriorForceDensity

void
IBFEHierarchyIntegrator::spreadBoundaryForceDensity(
    const int f_data_idx,
    NumericVector<double>& X_ghost,
    const double& time)
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();

    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    QBase* const qrule = d_fe_data_manager->getQuadratureRule();
    AutoPtr<QBase> qrule_face = QBase::build(qrule->type(), dim-1, qrule->get_order());

    System& system = equation_systems->get_system<System>(COORDINATES_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int d = 1; d < NDIM; ++d)
    {
        TBOX_ASSERT(dof_map.variable_type(0) == dof_map.variable_type(d));
    }
#endif

    AutoPtr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
    fe_face->attach_quadrature_rule(qrule_face.get());
    const std::vector<double>& JxW_face = fe_face->get_JxW();
    const std::vector<Point>& q_point_face = fe_face->get_xyz();
    const std::vector<std::vector<double> >& phi_face = fe_face->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi_face = fe_face->get_dphi();
    const std::vector<Point>& normal_face = fe_face->get_normals();

    // Loop over the patches to spread the nodal values onto the grid.
    const int level_num = d_fe_data_manager->getLevelNumber();
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > f_data = patch->getPatchData(f_data_idx);
        const SAMRAI::hier::Box<NDIM>& box = f_data->getGhostBox();

        // The relevant range of elements.
        const std::vector<Elem*>& active_patch_elems = d_fe_data_manager->getActivePatchGhostBoxElems(level_num, patch_num);
        const std::vector<Elem*>::const_iterator el_begin = active_patch_elems.begin();
        const std::vector<Elem*>::const_iterator el_end   = active_patch_elems.end();

        // Loop over the elements and compute the values to be spread and the
        // positions of the quadrature points.
        int qp_offset = 0;
        std::vector<double> T_bdry, X_bdry;
        for (std::vector<Elem*>::const_iterator cit = el_begin; cit != el_end; ++cit)
        {
            Elem* const elem = *cit;

            // Loop over the element boundaries.
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries and physical boundaries with
                // constraints.
                const short int boundary_id = mesh.boundary_info->boundary_id(elem,side);
                const bool at_physical_bdry = elem->neighbor(side) == NULL && !dof_map.is_periodic_boundary(boundary_id);
                const bool normal_dirichlet_bdry = boundary_id == NORMAL_DIRICHLET_BOUNDARY_ID || boundary_id == NORMAL_AND_TANGENTIAL_DIRICHLET_BOUNDARY_ID;
                if (at_physical_bdry && !normal_dirichlet_bdry)
                {
                    fe_face->reinit(elem, side);

                    // Determine the degrees of freedom for the current element.
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        dof_map.dof_indices(elem, dof_indices[i], i);
                    }

                    // Loop over boundary quadrature points.
                    T_bdry.resize(T_bdry.size()+NDIM*qrule_face->n_points(),0.0);
                    X_bdry.resize(X_bdry.size()+NDIM*qrule_face->n_points(),0.0);
                    for (unsigned int qp = 0; qp < qrule_face->n_points(); ++qp)
                    {
                        // Compute the value of the first Piola-Kirchoff stress
                        // tensor at the quadrature point and evaluate the
                        // transmission force density.
                        const Point& s_qp = q_point_face[qp];
                        const Point& X_qp = IBTK::compute_coordinate(qp,X_ghost,phi_face,dof_indices);
                        const TensorValue<double> dX_ds = IBTK::compute_coordinate_mapping_jacobian(qp,X_ghost,dphi_face,dof_indices);
                        const TensorValue<double> P = d_PK1_stress_function(dX_ds,X_qp,s_qp,elem,time,d_PK1_stress_function_ctx);
                        const VectorValue<double> F = P*normal_face[qp]*JxW_face[qp];

                        const TensorValue<double> dX_ds_inv_trans = IBTK::tensor_inverse_transpose(dX_ds, NDIM);
                        VectorValue<double> n = dX_ds_inv_trans*normal_face[qp];
                        n /= n.size();
                        const VectorValue<double> F_n = (F*n)*n;

                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            T_bdry[NDIM*(qp+qp_offset)+i] = -F_n(i);
                            X_bdry[NDIM*(qp+qp_offset)+i] = X_qp(i);
                        }
                    }

                    qp_offset += qrule_face->n_points();
                }
            }
        }

        // Spread the boundary forces to the grid.
        if (qp_offset > 0)
        {
            IBTK::LEInteractor::spread(f_data, T_bdry, NDIM, X_bdry, NDIM, patch, box, d_fe_data_manager->getSpreadWeightingFunction());
        }
    }
    return;
}// spreadBoundaryForceDensity

void
IBFEHierarchyIntegrator::initializeCoordinates()
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& coords_system = equation_systems->get_system<System>(COORDINATES_SYSTEM_NAME);
    const unsigned int coords_system_number = coords_system.number();
    NumericVector<double>& X_coords = *coords_system.solution;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(coords_system_number) > 0)
        {
            libmesh_assert(n->n_vars(coords_system_number) == NDIM);
            const Point& s = *n;
            const Point X = d_coordinate_mapping_function == NULL ? s : d_coordinate_mapping_function(s, d_coordinate_mapping_function_ctx);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int dof_index = n->dof_number(coords_system_number,d,0);
                X_coords.set(dof_index,X(d));
            }
        }
    }
    X_coords.close();
    X_coords.localize(*coords_system.current_local_solution);
    return;
}// initializeCoordinates

void
IBFEHierarchyIntegrator::updateCoordinateMapping()
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& coords_system = equation_systems->get_system<System>(COORDINATES_SYSTEM_NAME);
    const unsigned int coords_system_number = coords_system.number();
    NumericVector<double>& X_coords = *coords_system.solution;
    System& coords_mapping_system = equation_systems->get_system<System>(COORDINATE_MAPPING_SYSTEM_NAME);
    const unsigned int coords_mapping_system_number = coords_mapping_system.number();
    NumericVector<double>& dX_coords = *coords_mapping_system.solution;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(coords_system_number) > 0)
        {
            libmesh_assert(n->n_vars(coords_system_number) == NDIM);
            libmesh_assert(n->n_vars(coords_mapping_system_number) == NDIM);
            const Point& s = *n;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int X_dof_index = n->dof_number(coords_system_number,d,0);
                const int dX_dof_index = n->dof_number(coords_mapping_system_number,d,0);
                dX_coords.set(dX_dof_index,X_coords(X_dof_index)-s(d));
            }
        }
    }
    dX_coords.close();
    dX_coords.localize(*coords_mapping_system.current_local_solution);
    return;
}// updateCoordinateMapping

void
IBFEHierarchyIntegrator::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    // Read data members from input database.
    d_end_time = db->getDoubleWithDefault("end_time", d_end_time);
    d_grow_dt = db->getDoubleWithDefault("grow_dt", d_grow_dt);

    d_max_integrator_steps = db->getIntegerWithDefault("max_integrator_steps", d_max_integrator_steps);

    d_num_cycles = db->getIntegerWithDefault("num_cycles", d_num_cycles);

    d_regrid_interval = db->getIntegerWithDefault("regrid_interval", d_regrid_interval);

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault("dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault("dt_max_time_min", d_dt_max_time_min);

    d_split_interior_and_bdry_forces = db->getBoolWithDefault("split_interior_and_bdry_forces", d_split_interior_and_bdry_forces);

    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    if (!is_from_restart)
    {
        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);
    }
    return;
}// getFromInput

void
IBFEHierarchyIntegrator::getFromRestart()
{
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> restart_db = SAMRAI::tbox::RestartManager::getManager()->getRootDatabase();
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

    int ver = db->getInteger("IB_FE_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_FE_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    d_split_interior_and_bdry_forces = db->getBool("d_split_interior_and_bdry_forces");
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

    if (db->keyExists("instrument_names"))
    {
        const int sz = db->getInteger("instrument_names_sz");
        std::vector<std::string> instrument_names(sz);
        db->getStringArray("instrument_names", &instrument_names[0], sz);
        IBInstrumentationSpec::setInstrumentNames(instrument_names);
    }

    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBFEHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
