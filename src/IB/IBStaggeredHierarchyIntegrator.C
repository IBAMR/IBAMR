// Filename: IBStaggeredHierarchyIntegrator.C
// Last modified: <23.Sep.2008 14:06:02 griffith@box230.cims.nyu.edu>
// Created on 08 May 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "IBStaggeredHierarchyIntegrator.h"

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
#include <ibamr/IBAnchorPointSpec.h>
#include <ibamr/INSStaggeredIntermediateVelocityBcCoef.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/INSStaggeredProjectionBcCoef.h>
#include <ibamr/INSStaggeredVelocityBcCoef.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LNodeIndexData2.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScMultiVec.h>
#include <ibtk/PETScSAMRAIVectorReal.h>

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
#include <IndexData.h>
#include <Patch.h>
#include <PatchCellDataOpsReal.h>
#include <VariableDatabase.h>
#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <iterator>
#include <limits>
#include <numeric>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_SC_STABLEDT_F77 F77_FUNC_(navier_stokes_sc_stabledt2d, NAVIER_STOKES_SC_STABLEDT2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_SC_STABLEDT_F77 F77_FUNC_(navier_stokes_sc_stabledt3d, NAVIER_STOKES_SC_STABLEDT3D)
#endif

extern "C"
{
    void
    NAVIER_STOKES_SC_STABLEDT_F77(
        const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
#endif
        double&
                                  );
}

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
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_get_stable_timestep;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_regrid_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_new_levels;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_time_dependent_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_data_to_preadvance_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_hierarchy_configuration;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_apply_gradient_detector;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CONSERVATIVE_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Version of IBStaggeredHierarchyIntegrator restart file data.
static const int IB_STAGGERED_HIERARCHY_INTEGRATOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStaggeredHierarchyIntegrator::IBStaggeredHierarchyIntegrator(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    SAMRAI::tbox::Pointer<HierarchyProjector> hier_projector,
    SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> force_strategy,
    bool register_for_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT(!force_strategy.isNull());
#endif
    d_object_name = object_name;
    d_registered_for_restart = register_for_restart;
    d_hierarchy = hierarchy;
    d_hier_projector = hier_projector;
    d_force_strategy = force_strategy;

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }

    // Set some default values.
    d_u_scale = 1.0;
    d_p_scale = 1.0;
    d_f_scale = 1.0;

    d_use_exact_jacobian = false;

    d_start_time = 0.0;
    d_end_time = std::numeric_limits<double>::max();
    d_grow_dt = 2.0;
    d_max_integrator_steps = std::numeric_limits<int>::max();

    d_using_vorticity_tagging = false;
    d_omega_max = 0.0;

    d_normalize_pressure = false;

    d_regrid_interval = 1;
    d_old_dt = -1.0;
    d_integrator_time = std::numeric_limits<double>::quiet_NaN();
    d_integrator_step = std::numeric_limits<int>::max();

    d_conservation_form = false;

    d_output_u = false;
    d_output_p = false;
    d_output_f = false;
    d_output_omega = false;
    d_output_div_u = false;

    d_rho    = std::numeric_limits<double>::quiet_NaN();
    d_mu     = std::numeric_limits<double>::quiet_NaN();
    d_lambda = std::numeric_limits<double>::quiet_NaN();

    d_cfl = 0.9;

    d_dt_max = std::numeric_limits<double>::max();
    d_dt_max_time_max = std::numeric_limits<double>::max();
    d_dt_max_time_min = -(d_dt_max_time_max-std::numeric_limits<double>::epsilon());

    d_is_initialized = false;

    d_do_log = false;

    // Setup default boundary condition objects that specify homogeneous
    // Dirichlet boundary conditions for the velocity.
    d_default_u_bc_coef = new SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>(
        d_object_name+"::default_u_bc_coef", SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL));
    for (int d = 0; d < NDIM; ++d)
    {
        d_default_u_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_u_bc_coef->setBoundaryValue(2*d+1,0.0);
    }

    registerVelocityPhysicalBcCoefs(
        std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(
            NDIM,d_default_u_bc_coef));

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Set the problem coefs.
    d_problem_coefs = new INSCoefs(d_rho, d_mu, d_lambda);

    // Determine the ghost cell width required for cell-centered spreading and
    // interpolating.
    const int stencil_size = IBTK::LEInteractor::getStencilSize(d_delta_fcn);
    d_ghosts = 8; // XXXX int(floor(0.5*double(stencil_size)))+1;
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";
    SAMRAI::tbox::pout << "please fix me!!!!\n";

    // Get the Lagrangian Data Manager.
    d_lag_data_manager = IBTK::LDataManager::getManager(
        d_object_name+"::LDataManager", d_ghosts, d_registered_for_restart);

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager =
        SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var =
        new SAMRAI::pdat::CellVariable<NDIM,double>("cc_var");
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > sc_var =
        new SAMRAI::pdat::SideVariable<NDIM,double>("sc_var");

    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy);
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, hierarchy);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::advanceHierarchy()");
        t_get_stable_timestep = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::getStableTimestep()");
        t_regrid_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::regridHierarchy()");
        t_synchronize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBStaggeredHierarchyIntegrator::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// IBStaggeredHierarchyIntegrator

IBStaggeredHierarchyIntegrator::~IBStaggeredHierarchyIntegrator()
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

    delete d_helmholtz_spec;
    delete d_default_u_bc_coef;
    if (!d_u_bc_coefs.empty())
    {
        for (int d = 0; d < NDIM; ++d)
        {
            delete d_u_bc_coefs[d];
            delete d_u_star_bc_coefs[d];
        }
        delete d_p_bc_coef;
        delete d_phi_bc_coef;
    }
    return;
}// ~IBStaggeredHierarchyIntegrator

const std::string&
IBStaggeredHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
IBStaggeredHierarchyIntegrator::registerVelocityInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> u_init)
{
    d_u_init = u_init;
    return;
}// registerVelocityInitialConditions

void
IBStaggeredHierarchyIntegrator::registerVelocityPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::registerVelocityPhysicalBcCoefs():\n"
                   << "  velocity boundary conditions must be registered prior to initialization\n"
                   << "  of the hierarchy integrator object." << std::endl);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned l = 0; l < u_bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(u_bc_coefs[l] != NULL);
    }
#endif
    if (!d_u_bc_coefs.empty())
    {
        for (int d = 0; d < NDIM; ++d)
        {
            delete d_u_bc_coefs[d];
            delete d_u_star_bc_coefs[d];
        }
        delete d_p_bc_coef;
        delete d_phi_bc_coef;
    }
    d_u_bc_coefs.clear();
    d_u_bc_coefs.resize(NDIM,NULL);
    for (int d = 0; d < NDIM; ++d)
    {
        d_u_bc_coefs[d] = new INSStaggeredVelocityBcCoef(d,*d_problem_coefs,u_bc_coefs);
    }
    d_u_star_bc_coefs.clear();
    d_u_star_bc_coefs.resize(NDIM,NULL);
    for (int d = 0; d < NDIM; ++d)
    {
        d_u_star_bc_coefs[d] = new INSStaggeredIntermediateVelocityBcCoef(d,u_bc_coefs);
    }
    d_p_bc_coef = new INSStaggeredPressureBcCoef(*d_problem_coefs,u_bc_coefs);
    d_phi_bc_coef = new INSStaggeredProjectionBcCoef(u_bc_coefs);
    d_hier_projector->setVelocityPhysicalBcCoefs(d_u_star_bc_coefs);
    d_hier_projector->setPressurePhysicalBcCoef(d_phi_bc_coef);
    return;
}// registerVelocityPhysicalBcCoefs

void
IBStaggeredHierarchyIntegrator::registerPressureInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> p_init)
{
    d_p_init = p_init;
    return;
}// registerPressureInitialConditions

void
IBStaggeredHierarchyIntegrator::registerBodyForceSpecification(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> body_force_set)
{
    d_body_force_set = body_force_set;
    return;
}// registerBodyForceSpecification

void
IBStaggeredHierarchyIntegrator::registerLNodeInitStrategy(
    SAMRAI::tbox::Pointer<IBTK::LNodeInitStrategy> lag_init)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!lag_init.isNull());
#endif
    d_lag_init = lag_init;
    d_lag_data_manager->registerLNodeInitStrategy(d_lag_init);
    return;
}// registerLNodeInitStrategy

void
IBStaggeredHierarchyIntegrator::freeLNodeInitStrategy()
{
    d_lag_init.setNull();
    d_lag_data_manager->freeLNodeInitStrategy();
    return;
}// freeLNodeInitStrategy

void
IBStaggeredHierarchyIntegrator::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    d_lag_data_manager->registerVisItDataWriter(d_visit_writer);
    return;
}// registerVisItDataWriter

void
IBStaggeredHierarchyIntegrator::registerLagSiloDataWriter(
    SAMRAI::tbox::Pointer<IBTK::LagSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!silo_writer.isNull());
#endif
    d_silo_writer = silo_writer;
    d_lag_data_manager->registerLagSiloDataWriter(d_silo_writer);
    return;
}// registerLagSiloDataWriter

#if (NDIM == 3)
void
IBStaggeredHierarchyIntegrator::registerLagM3DDataWriter(
    SAMRAI::tbox::Pointer<IBTK::LagM3DDataWriter> m3D_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!m3D_writer.isNull());
#endif
    d_m3D_writer = m3D_writer;
    d_lag_data_manager->registerLagM3DDataWriter(d_m3D_writer);
    return;
}// registerLagM3DDataWriter
#endif

void
IBStaggeredHierarchyIntegrator::registerLoadBalancer(
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
///  allow the IBStaggeredHierarchyIntegrator to be used as a hierarchy integrator.
///

void
IBStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!gridding_alg.isNull());
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

    // Initialize all variable contexts.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    d_current_context = var_db->getContext(d_object_name+"::CURRENT");
    d_new_context     = var_db->getContext(d_object_name+"::NEW"    );
    d_scratch_context = var_db->getContext(d_object_name+"::SCRATCH");

    // Initialize all variables.
    d_u_var          = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::u"          );
    d_u_cc_var       = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::u_cc",  NDIM);
    d_p_var          = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::p"          );
    d_f_var          = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::f"          );
    d_f_cc_var       = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::f_cc",  NDIM);
#if ( NDIM == 2)
    d_omega_var      = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::omega"      );
#endif
#if ( NDIM == 3)
    d_omega_var      = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::omega", NDIM);
    d_omega_norm_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::||omega||_2");
#endif
    d_div_u_var      = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::div_u"      );
    d_u_interp_var   = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::u_interp"   );

    // Create the default communication algorithms.
    d_fill_after_regrid = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    d_calgs["SYNCH_CURRENT_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    d_calgs["SYNCH_NEW_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();

    // Register state variables that are maintained by the
    // INSStaggeredHierarchyIntegrator.

    const SAMRAI::hier::IntVector<NDIM> cell_ghosts = CELLG;
    const SAMRAI::hier::IntVector<NDIM> side_ghosts = SIDEG;
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    registerVariable(d_u_current_idx, d_u_new_idx, d_u_scratch_idx,
                     d_u_var, side_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_u_cc_current_idx, d_u_cc_new_idx, d_u_cc_scratch_idx,
                     d_u_cc_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_p_current_idx, d_p_new_idx, d_p_scratch_idx,
                     d_p_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");

    registerVariable(d_f_current_idx, d_f_new_idx, d_f_scratch_idx,
                     d_f_var, side_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_f_cc_current_idx, d_f_cc_new_idx, d_f_cc_scratch_idx,
                     d_f_cc_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_omega_current_idx, d_omega_new_idx, d_omega_scratch_idx,
                     d_omega_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");
#if (NDIM == 3)
    registerVariable(d_omega_norm_current_idx, d_omega_norm_new_idx, d_omega_norm_scratch_idx,
                     d_omega_norm_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");
#endif
    registerVariable(d_div_u_current_idx, d_div_u_new_idx, d_div_u_scratch_idx,
                     d_div_u_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    // Register scratch variables.

    registerVariable(d_u_interp_idx, d_u_interp_var, d_ghosts);

    d_u_rhs_idx  = var_db->registerClonedPatchDataIndex(d_u_var, d_u_scratch_idx);
    d_u_nul_idx  = var_db->registerClonedPatchDataIndex(d_u_var, d_u_scratch_idx);
    d_u_half_idx = var_db->registerClonedPatchDataIndex(d_u_var, d_u_scratch_idx);
    d_n_idx      = var_db->registerClonedPatchDataIndex(d_u_var, d_u_scratch_idx);
    d_p_rhs_idx  = var_db->registerClonedPatchDataIndex(d_p_var, d_p_scratch_idx);
    d_p_nul_idx  = var_db->registerClonedPatchDataIndex(d_p_var, d_p_scratch_idx);

    d_scratch_data.setFlag(d_u_rhs_idx);
    d_scratch_data.setFlag(d_u_nul_idx);
    d_scratch_data.setFlag(d_u_half_idx);
    d_scratch_data.setFlag(d_n_idx);
    d_scratch_data.setFlag(d_p_rhs_idx);
    d_scratch_data.setFlag(d_p_nul_idx);

    // Register variables for plotting.
    if (!d_visit_writer.isNull())
    {
        if (d_output_u)
        {
            d_visit_writer->registerPlotQuantity(
                "u", "VECTOR", d_u_cc_current_idx, 0, d_u_scale);

            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    "u_"+stream.str(), "SCALAR", d_u_cc_current_idx, d, d_u_scale);
            }
        }

        if (d_output_p)
        {
            d_visit_writer->registerPlotQuantity(
                "p", "SCALAR", d_p_current_idx, 0, d_p_scale);
        }

        if (d_output_f)
        {
            d_visit_writer->registerPlotQuantity(
                "f", "VECTOR", d_f_cc_current_idx, 0, d_f_scale);

            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    "f_"+stream.str(), "SCALAR", d_f_cc_current_idx, d, d_f_scale);
            }
        }

        if (d_output_omega)
        {
            d_visit_writer->registerPlotQuantity(
                "omega", (NDIM == 2) ? "SCALAR" : "VECTOR", d_omega_current_idx);
#if (NDIM == 3)
            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    "omega_"+stream.str(), "SCALAR", d_omega_current_idx, d);
            }
#endif
        }

        if (d_output_div_u)
        {
            d_visit_writer->registerPlotQuantity(
                "div u", "SCALAR", d_div_u_current_idx);
        }
    }

    // Create several refinement communications algorithms, used in filling
    // ghost cell data.
    std::string ralg_name;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;

    d_ralgs["u_interp->u_interp::S->S::CONSERVATIVE_LINEAR_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_u_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["u_interp->u_interp::S->S::CONSERVATIVE_LINEAR_REFINE"]->
        registerRefine(d_u_interp_idx, // destination
                       d_u_interp_idx, // source
                       d_u_interp_idx, // temporary work space
                       refine_operator);
    d_rstrategies["u_interp->u_interp::S->S::CONSERVATIVE_LINEAR_REFINE"] =
        new IBTK::CartExtrapPhysBdryOp(d_u_interp_idx, BDRY_EXTRAP_TYPE);

    // Setup the Hierarchy math operations object.
    d_hier_math_ops = new IBTK::HierarchyMathOps(d_object_name+"::HierarchyMathOps", d_hierarchy);
    d_is_managing_hier_math_ops = true;
    d_hier_projector->setHierarchyMathOps(d_hier_math_ops);

    // Set the current integration time.
    if (!SAMRAI::tbox::RestartManager::getManager()->isFromRestart())
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }

    // Setup physical boundary conditions helper.
    d_u_bc_helper = new INSStaggeredPhysicalBoundaryHelper();

    // Setup the Stokes operator.
    d_stokes_op = new INSStaggeredStokesOperator(
        *d_problem_coefs,
        d_u_bc_coefs, d_p_bc_coef,
        d_hier_math_ops);

    // Setup the convective operator.
    d_convective_op = new INSStaggeredPPMConvectiveOperator(
        *d_problem_coefs,
        d_conservation_form);

    // Setup the Helmholtz solver.
    d_helmholtz_spec = new SAMRAI::solv::PoissonSpecifications(d_object_name+"::helmholtz_spec");
    d_helmholtz_op = new IBTK::SCLaplaceOperator(
        d_object_name+"::Helmholtz Operator", *d_helmholtz_spec, d_u_star_bc_coefs, true);
    d_helmholtz_op->setHierarchyMathOps(d_hier_math_ops);

    d_helmholtz_solver_needs_init = true;
    d_helmholtz_solver = new IBTK::PETScKrylovLinearSolver(
        d_object_name+"::PETSc Krylov solver", "helmholtz_");
    d_helmholtz_solver->setInitialGuessNonzero(false);
    d_helmholtz_solver->setOperator(d_helmholtz_op);

    // Setup the projection preconditioner.
    d_projection_pc_needs_init = true;
#if 0
    d_projection_pc = new INSStaggeredProjectionPreconditioner(
        *d_problem_coefs,
        d_normalize_pressure,
        d_helmholtz_solver, d_hier_projector,
        d_hier_cc_data_ops, d_hier_sc_data_ops, d_hier_math_ops);
#else
    TBOX_ASSERT(false);  // XXXX
#endif

    // Indicate that the integrator has been initialized.
    d_is_initialized = true;

    t_initialize_hierarchy_integrator->stop();
    return;
}// initializeHierarchyIntegrator

double
IBStaggeredHierarchyIntegrator::initializeHierarchy()
{
    t_initialize_hierarchy->start();

    if (!d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchy():\n" <<
                   "  must call initializeHierarchyIntegrator() prior to calling initializeHierarchy()." << std::endl);
    }

    // Initialize the patch hierarchy.
    const bool initial_time = !SAMRAI::tbox::RestartManager::getManager()->isFromRestart();

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

        // After data on each level is initialized at simulation start time,
        // coarser levels are synchronized with finer levels that didn't exist
        // when the coarser level initial data was set.
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();

        if (finest_ln > 0)
        {
            synchronizeNewLevels(d_hierarchy, coarsest_ln, finest_ln,
                                 d_start_time, initial_time);
        }
    }

    // Reset the Lagrangian data distribution.
    d_lag_data_manager->beginDataRedistribution();
    d_lag_data_manager->endDataRedistribution();

    // Update the workload.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    d_lag_data_manager->updateWorkloadData(
        coarsest_ln, finest_ln);

    // Indicate that the force strategy needs to be re-initialized.
    d_force_strategy_needs_init = true;

    // The next timestep is given by the minimum allowable timestep over all
    // levels in the patch hierarchy.
    const double dt_next = getStableTimestep(getCurrentContext());

    t_initialize_hierarchy->stop();
    return dt_next;
}// initializeHierarchy

double
IBStaggeredHierarchyIntegrator::advanceHierarchy(
    const double dt)
{
    t_advance_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
#endif

    PetscErrorCode ierr;

    const double current_time = d_integrator_time;
    const double new_time = d_integrator_time+dt;
    const bool initial_time = SAMRAI::tbox::MathUtilities<double>::equalEps(d_integrator_time,d_start_time);

    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = dt;

    // Set the guess for the initial pressure to zero.
    if (initial_time)
    {
        d_hier_cc_data_ops->setToScalar(d_p_current_idx, 0.0);
    }

    // Regrid the patch hierarchy.
    const bool do_regrid = ((d_regrid_interval == 0)
                            ? false
                            : (d_integrator_step % d_regrid_interval == 0));
    if (do_regrid)
    {
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): regridding prior to timestep " << d_integrator_step << "\n";
        regridHierarchy();
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Reset the various Lagrangian objects.  XXXX
    d_X_data     .resize(finest_ln+1);
    d_X_mid_data .resize(finest_ln+1);
    d_X_new_data .resize(finest_ln+1);
    d_X_half_data.resize(finest_ln+1);
    d_U_half_data.resize(finest_ln+1);
    d_F_half_data.resize(finest_ln+1);
    d_J_mat.resize(finest_ln+1,static_cast<Mat>(NULL));
    d_strct_mat.resize(finest_ln+1,static_cast<Mat>(NULL));
    d_strct_pc_mat.resize(finest_ln+1,static_cast<Mat>(NULL));
    d_strct_ksp.resize(finest_ln+1,static_cast<KSP>(NULL));
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(ln == finest_ln);
#endif
            d_X_data     [ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,ln);
            d_X_mid_data [ln] = d_lag_data_manager->createLNodeLevelData("X_mid" ,ln,NDIM);
            d_X_new_data [ln] = d_lag_data_manager->createLNodeLevelData("X_new" ,ln,NDIM);
            d_X_half_data[ln] = d_lag_data_manager->createLNodeLevelData("X_half",ln,NDIM);
            d_U_half_data[ln] = d_lag_data_manager->createLNodeLevelData("U_half",ln,NDIM);
            d_F_half_data[ln] = d_lag_data_manager->createLNodeLevelData("F_half",ln,NDIM);
        }
    }

    // Setup the boundary conditions specifications.
    for (int d = 0; d < NDIM; ++d)
    {
        INSStaggeredVelocityBcCoef* u_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(d_u_bc_coefs[d]);
        u_bc_coef->setTimeInterval(current_time,new_time);
    }
    INSStaggeredPressureBcCoef* p_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(d_p_bc_coef);
    p_bc_coef->setTimeInterval(current_time,new_time);
    p_bc_coef->setVelocityCurrentPatchDataIndex(d_u_current_idx);
    p_bc_coef->setVelocityNewPatchDataIndex(d_u_new_idx);

    // Set the current time interval.
    d_stokes_op->setTimeInterval(current_time,new_time);
    d_projection_pc->setTimeInterval(current_time,new_time);
    d_force_strategy->setTimeInterval(current_time, new_time);

    // (Re)initialize the force strategy.
    if (d_force_strategy_needs_init)
    {
        resetLagrangianForceStrategy(current_time, initial_time);
        d_force_strategy_needs_init = false;
    }

    // Synchronize current state data.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_CURRENT_STATE_DATA"][ln]->coarsenData();
    }

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data    ,     new_time);
    }

    d_hier_sc_data_ops->setToScalar(d_u_scratch_idx,0.0);
    d_hier_sc_data_ops->setToScalar(d_u_rhs_idx,0.0);
    d_hier_sc_data_ops->setToScalar(d_u_nul_idx,0.0);
    d_hier_sc_data_ops->setToScalar(d_u_half_idx,0.0);
    d_hier_sc_data_ops->setToScalar(d_n_idx,0.0);

    d_hier_cc_data_ops->setToScalar(d_p_scratch_idx,0.0);
    d_hier_cc_data_ops->setToScalar(d_p_rhs_idx,0.0);
    d_hier_cc_data_ops->setToScalar(d_p_nul_idx,1.0);

    // Setup fluid vectors.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > u_sol_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::u_sol_vec", d_hierarchy, 0, finest_ln);
    u_sol_vec->addComponent(d_u_var,d_u_scratch_idx,d_wgt_sc_idx,d_hier_sc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > u_rhs_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::u_rhs_vec", d_hierarchy, 0, finest_ln);
    u_rhs_vec->addComponent(d_u_var,d_u_rhs_idx,d_wgt_sc_idx,d_hier_sc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_sol_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::fluid_sol_vec", d_hierarchy, 0, finest_ln);
    fluid_sol_vec->addComponent(d_u_var,d_u_scratch_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    fluid_sol_vec->addComponent(d_p_var,d_p_scratch_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_rhs_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::fluid_rhs_vec", d_hierarchy, 0, finest_ln);
    fluid_rhs_vec->addComponent(d_u_var,d_u_rhs_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    fluid_rhs_vec->addComponent(d_p_var,d_p_rhs_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_nul_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::fluid_nul_vec", d_hierarchy, 0, finest_ln);
    fluid_nul_vec->addComponent(d_u_var,d_u_nul_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    fluid_nul_vec->addComponent(d_p_var,d_p_nul_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    Vec petsc_fluid_sol_vec = IBTK::PETScSAMRAIVectorReal<double>::createPETScVector(fluid_sol_vec, PETSC_COMM_WORLD);
    Vec petsc_fluid_rhs_vec = IBTK::PETScSAMRAIVectorReal<double>::createPETScVector(fluid_rhs_vec, PETSC_COMM_WORLD);
    Vec petsc_fluid_nul_vec = IBTK::PETScSAMRAIVectorReal<double>::createPETScVector(fluid_nul_vec, PETSC_COMM_WORLD);

    // Setup the structure vectors.
    Vec petsc_structure_sol_vec;
    ierr = VecDuplicate(d_X_data[finest_ln]->getGlobalVec(),&petsc_structure_sol_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSet(petsc_structure_sol_vec,0.0);

    Vec petsc_structure_rhs_vec;
    ierr = VecDuplicate(d_X_data[finest_ln]->getGlobalVec(),&petsc_structure_rhs_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSet(petsc_structure_rhs_vec,0.0);

    Vec petsc_structure_nul_vec;
    ierr = VecDuplicate(d_X_data[finest_ln]->getGlobalVec(),&petsc_structure_nul_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSet(petsc_structure_nul_vec,0.0);

    // Set the initial guess for u(n+1) to equal u(n).
    d_hier_sc_data_ops->copyData(d_u_new_idx, d_u_current_idx);
    d_hier_sc_data_ops->copyData(d_u_scratch_idx, d_u_current_idx);

    // Set the initial guess for p(n+1/2) to equal p(n-1/2).
    d_hier_cc_data_ops->copyData(d_p_new_idx, d_p_current_idx);
    d_hier_cc_data_ops->copyData(d_p_scratch_idx, d_p_current_idx);

    // Set the initial guess for X(n+1) to equal X(n).
    Vec X_vec     = d_X_data    [finest_ln]->getGlobalVec();
    Vec X_new_vec = d_X_new_data[finest_ln]->getGlobalVec();
    ierr = VecCopy(X_vec,X_new_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecCopy(X_vec,petsc_structure_sol_vec);  IBTK_CHKERRQ(ierr);

    // Setup solvers.
    d_helmholtz_spec->setCConstant((d_rho/dt)+0.5*d_lambda);
    d_helmholtz_spec->setDConstant(          -0.5*d_mu    );

    d_helmholtz_op->setPoissonSpecifications(*d_helmholtz_spec);
    d_helmholtz_op->setPhysicalBcCoefs(d_u_star_bc_coefs);
    d_helmholtz_op->setHomogeneousBc(true);
    d_helmholtz_op->setTime(new_time);
    d_helmholtz_op->setHierarchyMathOps(d_hier_math_ops);

    d_helmholtz_solver->setInitialGuessNonzero(false);
    d_helmholtz_solver->setOperator(d_helmholtz_op);
    if (d_helmholtz_solver_needs_init || !SAMRAI::tbox::MathUtilities<double>::equalEps(dt,d_old_dt))
    {
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << ": Initializing Helmholtz solver; dt = " << dt << "\n";
        d_helmholtz_solver->initializeSolverState(*u_sol_vec,*u_rhs_vec);
        d_helmholtz_solver_needs_init = false;
    }

    if (d_projection_pc_needs_init)
    {
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << ": Initializing projection preconditioner; dt = " << dt << "\n";
        d_projection_pc->initializeSolverState(*fluid_sol_vec,*fluid_rhs_vec);
    }

    d_helmholtz_solver->enableLogging(false);
    d_projection_pc->enableLogging(false);

    // Hierarchy ghost cell transaction components.
    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent u_bc_component(d_u_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_u_bc_coefs);

    // Compute the right-hand side terms.
    SAMRAI::solv::PoissonSpecifications rhs_spec(d_object_name+"::rhs_spec");
    rhs_spec.setCConstant((d_rho/dt)-0.5*d_lambda);
    rhs_spec.setDConstant(          +0.5*d_mu    );
    d_hier_sc_data_ops->copyData(d_u_scratch_idx, d_u_current_idx);
    d_hier_math_ops->laplace(
        d_u_rhs_idx, d_u_var,
        rhs_spec,
        d_u_scratch_idx, d_u_var,
        d_u_bdry_bc_fill_op /* XXXX d_u_bdry_extrap_fill_op XXXX */, current_time);
    if (!d_body_force_set.isNull())
    {
        d_body_force_set->setDataOnPatchHierarchy(
            d_f_scratch_idx, d_f_var, d_hierarchy, current_time+0.5*dt);
        d_hier_sc_data_ops->add(d_u_rhs_idx, d_f_scratch_idx, d_u_rhs_idx);
    }

    // Setup the multivectors.
    Vec petsc_sol_multivec;
    Vec petsc_sol_multivecs[2] = { petsc_fluid_sol_vec , petsc_structure_sol_vec };
    ierr = IBTK::VecCreateMultiVec(PETSC_COMM_WORLD, 2, petsc_sol_multivecs, &petsc_sol_multivec);  IBTK_CHKERRQ(ierr);
    d_petsc_x_vec = petsc_sol_multivec;

    Vec petsc_rhs_multivec;
    Vec petsc_rhs_multivecs[2] = { petsc_fluid_rhs_vec , petsc_structure_rhs_vec };
    ierr = IBTK::VecCreateMultiVec(PETSC_COMM_WORLD, 2, petsc_rhs_multivecs, &petsc_rhs_multivec);  IBTK_CHKERRQ(ierr);
    d_petsc_f_vec = petsc_sol_multivec;

    Vec petsc_nul_multivec;
    Vec petsc_nul_multivecs[2] = { petsc_fluid_nul_vec , petsc_structure_nul_vec };
    ierr = IBTK::VecCreateMultiVec(PETSC_COMM_WORLD, 2, petsc_nul_multivecs, &petsc_nul_multivec);  IBTK_CHKERRQ(ierr);

    Vec petsc_res_multivec;
    ierr = VecDuplicate(petsc_rhs_multivec, &petsc_res_multivec);  IBTK_CHKERRQ(ierr);

    int local_sz;
    ierr = VecGetLocalSize(petsc_sol_multivec, &local_sz);  IBTK_CHKERRQ(ierr);

    // Setup the nullspace object.
    MatNullSpace petsc_nullsp;
    Vec vecs[1] = { petsc_nul_multivec };
    static const PetscTruth has_cnst = PETSC_FALSE;
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, has_cnst, 1, vecs, &petsc_nullsp);  IBTK_CHKERRQ(ierr);

    // Setup the PETSc solver.
    static const std::string options_prefix = "ib_";
    ierr = SNESCreate(PETSC_COMM_WORLD, &d_petsc_snes);  IBTK_CHKERRQ(ierr);
    ierr = SNESSetFunction(d_petsc_snes, petsc_res_multivec, IBStaggeredHierarchyIntegrator::FormFunction_SAMRAI, static_cast<void*>(this));  IBTK_CHKERRQ(ierr);

    Mat petsc_jac;
    if (d_use_exact_jacobian)
    {
        ierr = MatCreateShell(PETSC_COMM_WORLD, local_sz, local_sz, PETSC_DETERMINE, PETSC_DETERMINE, static_cast<void*>(this), &petsc_jac);  IBTK_CHKERRQ(ierr);
        ierr = MatShellSetOperation(petsc_jac, MATOP_MULT, reinterpret_cast<void(*)(void)>(IBStaggeredHierarchyIntegrator::MatVecMult_SAMRAI));  IBTK_CHKERRQ(ierr);
    }
    else
    {
        ierr = MatCreateSNESMF(d_petsc_snes, &petsc_jac);  IBTK_CHKERRQ(ierr);
    }

    PetscContainer petsc_jac_ctx;
    ierr = PetscContainerCreate(PETSC_COMM_WORLD, &petsc_jac_ctx);  IBTK_CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(petsc_jac_ctx, static_cast<void*>(this));  IBTK_CHKERRQ(ierr);
    ierr = PetscObjectCompose(reinterpret_cast<PetscObject>(petsc_jac), "petsc_jac_ctx", reinterpret_cast<PetscObject>(petsc_jac_ctx));  IBTK_CHKERRQ(ierr);
    ierr = MatShellSetOperation(petsc_jac, MATOP_GET_VECS, reinterpret_cast<void(*)(void)>(IBStaggeredHierarchyIntegrator::MatGetVecs_SAMRAI)); IBTK_CHKERRQ(ierr);

    ierr = SNESSetJacobian(d_petsc_snes, petsc_jac, petsc_jac, IBStaggeredHierarchyIntegrator::FormJacobian_SAMRAI, static_cast<void*>(this));  IBTK_CHKERRQ(ierr);

    KSP petsc_ksp;
    ierr = SNESGetKSP(d_petsc_snes, &petsc_ksp);  IBTK_CHKERRQ(ierr);
    if (d_normalize_pressure)
    {
        ierr = KSPSetNullSpace(petsc_ksp, petsc_nullsp);  IBTK_CHKERRQ(ierr);
    }

    PC petsc_pc;
    ierr = KSPGetPC(petsc_ksp, &petsc_pc);  IBTK_CHKERRQ(ierr);
    ierr = PCSetType(petsc_pc, PCCOMPOSITE);  IBTK_CHKERRQ(ierr);
    int composite_pc_counter = 0;

    PC petsc_fluid_pc;
    static const int fluid_pc_number = composite_pc_counter++;
    ierr = PCCompositeAddPC(petsc_pc, PCSHELL);  IBTK_CHKERRQ(ierr);
    ierr = PCCompositeGetPC(petsc_pc, fluid_pc_number, &petsc_fluid_pc);  IBTK_CHKERRQ(ierr);
    ierr = PCShellSetContext(petsc_fluid_pc, static_cast<void*>(this));  IBTK_CHKERRQ(ierr);
    ierr = PCShellSetApply(petsc_fluid_pc, IBStaggeredHierarchyIntegrator::PCApplyFluid_SAMRAI);  IBTK_CHKERRQ(ierr);
    ierr = PCShellSetName(petsc_fluid_pc, "IBStaggeredHierarchyIntegrator Fluid PC");  IBTK_CHKERRQ(ierr);

    PC petsc_strct_pc;
    static const int strct_pc_number = composite_pc_counter++;
    ierr = PCCompositeAddPC(petsc_pc, PCSHELL);  IBTK_CHKERRQ(ierr);
    ierr = PCCompositeGetPC(petsc_pc, strct_pc_number, &petsc_strct_pc);  IBTK_CHKERRQ(ierr);
    ierr = PCShellSetContext(petsc_strct_pc, static_cast<void*>(this));  IBTK_CHKERRQ(ierr);
    ierr = PCShellSetApply(petsc_strct_pc, IBStaggeredHierarchyIntegrator::PCApplyStrct_SAMRAI);  IBTK_CHKERRQ(ierr);
    ierr = PCShellSetName(petsc_strct_pc, "IBStaggeredHierarchyIntegrator Strct PC");  IBTK_CHKERRQ(ierr);

    ierr = SNESSetOptionsPrefix(d_petsc_snes, options_prefix.c_str());  IBTK_CHKERRQ(ierr);
    ierr = SNESSetFromOptions(d_petsc_snes);  IBTK_CHKERRQ(ierr);

    // Setup inhomogeneous boundary conditions.
    d_u_bc_helper->cacheBcCoefData(d_u_var, d_u_bc_coefs, new_time, SAMRAI::hier::IntVector<NDIM>(SIDEG), d_hierarchy);
    d_u_bc_helper->zeroValuesAtDirichletBoundaries(fluid_rhs_vec->getComponentDescriptorIndex(0));
    TBOX_ASSERT(fluid_sol_vec->getComponentDescriptorIndex(0) == d_u_scratch_idx);
    d_u_bdry_bc_fill_op->fillData(new_time);

    // Solve for u(n+1), p(n+1/2), and X(n+1) using (truncated) fixed point
    // iteration.
    static const int num_cycles = 3;
    for (int cycle = 0; cycle < num_cycles; ++cycle)
    {
        // Compute u_half := 0.5*(u(n)+u(n+1)).
        d_hier_sc_data_ops->linearSum(d_u_half_idx, 0.5, d_u_current_idx, 0.5, d_u_new_idx);

        // Compute X_mid := 0.5*(X(n)+X(n+1)).
        Vec X_vec     = d_X_data    [finest_ln]->getGlobalVec();
        Vec X_new_vec = d_X_new_data[finest_ln]->getGlobalVec();
        Vec X_mid_vec = d_X_mid_data[finest_ln]->getGlobalVec();
        ierr = VecCopy(X_vec,X_mid_vec);  IBTK_CHKERRQ(ierr);
        ierr = VecAXPBY(X_mid_vec,0.5,0.5,X_new_vec);  IBTK_CHKERRQ(ierr);
        d_X_mid_data[finest_ln]->beginGhostUpdate();
        d_X_mid_data[finest_ln]->endGhostUpdate();

        // Compute (u_half*grad)u_half.
        d_convective_op->applyConvectiveOperator(d_u_half_idx, d_n_idx);
        d_u_bc_helper->zeroValuesAtDirichletBoundaries(d_n_idx);

        // Setup the right-hand side vector.
        d_hier_sc_data_ops->axpy(fluid_rhs_vec->getComponentDescriptorIndex(0), -d_rho, d_n_idx, fluid_rhs_vec->getComponentDescriptorIndex(0));

        // Solve for u(n+1), p(n+1/2), and X(n+1).
        ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_sol_multivec));  IBTK_CHKERRQ(ierr);
        ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_rhs_multivec));  IBTK_CHKERRQ(ierr);
        ierr = SNESSolve(d_petsc_snes, petsc_rhs_multivec, petsc_sol_multivec);  IBTK_CHKERRQ(ierr);

        // Pull out solution components.
        d_hier_sc_data_ops->copyData(d_u_new_idx, fluid_sol_vec->getComponentDescriptorIndex(0));
        d_hier_cc_data_ops->copyData(d_p_new_idx, fluid_sol_vec->getComponentDescriptorIndex(1));
        ierr = VecCopy(petsc_structure_sol_vec, d_X_new_data[finest_ln]->getGlobalVec());  IBTK_CHKERRQ(ierr);

        // Reset the right-hand side vector.
        d_hier_sc_data_ops->axpy(fluid_rhs_vec->getComponentDescriptorIndex(0), +d_rho, d_n_idx, fluid_rhs_vec->getComponentDescriptorIndex(0));
    }

    /*
     * Compute the drag and lift coefficients by integrating the components of
     * the Lagrangian force field over the computational domain.
     */
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            std::ofstream drag_stream, lift_stream;
            if (SAMRAI::tbox::SAMRAI_MPI::getRank() == 0)
            {
                drag_stream.open("C_D.curve", SAMRAI::tbox::RestartManager::getManager()->isFromRestart() || !initial_time ? std::ios::app : std::ios::out);
                lift_stream.open("C_L.curve", SAMRAI::tbox::RestartManager::getManager()->isFromRestart() || !initial_time ? std::ios::app : std::ios::out);

                drag_stream.setf(std::ios_base::scientific);
                drag_stream.setf(std::ios_base::showpos);
                drag_stream.setf(std::ios_base::showpoint);
                drag_stream.width(16); drag_stream.precision(15);

                lift_stream.setf(std::ios_base::scientific);
                lift_stream.setf(std::ios_base::showpos);
                lift_stream.setf(std::ios_base::showpoint);
                lift_stream.width(16); lift_stream.precision(15);
            }

            Vec X_vec      = d_X_data     [ln]->getGlobalVec();
            Vec X_new_vec  = d_X_new_data [ln]->getGlobalVec();
            Vec X_half_vec = d_X_half_data[ln]->getGlobalVec();
            ierr = VecCopy(X_vec,X_half_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(X_half_vec,0.5,0.5,X_new_vec);  IBTK_CHKERRQ(ierr);
            Vec U_half_vec = d_U_half_data[ln]->getGlobalVec();
            ierr = VecCopy(X_vec,U_half_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(U_half_vec,1.0/d_dt,-1.0/d_dt,X_new_vec);  IBTK_CHKERRQ(ierr);

            d_force_strategy->computeLagrangianForce(
                d_F_half_data[ln], d_X_half_data[ln], d_U_half_data[ln],
                d_hierarchy, ln, d_current_time+0.5*d_dt, d_lag_data_manager);

            double F_D = 0.0;
            double F_L = 0.0;
            for (int i = 0; i < d_F_half_data[ln]->getLocalNodeCount(); ++i)
            {
                F_D -= (*d_F_half_data[ln])(i,0);
                F_L -= (*d_F_half_data[ln])(i,1);
            }

            F_D = SAMRAI::tbox::SAMRAI_MPI::sumReduction(F_D);
            F_L = SAMRAI::tbox::SAMRAI_MPI::sumReduction(F_L);

            static const double radius = 0.5;
            if (SAMRAI::tbox::SAMRAI_MPI::getRank() == 0)
            {
                drag_stream << d_current_time+0.5*d_dt << " " << F_D/radius << std::endl;
                lift_stream << d_current_time+0.5*d_dt << " " << F_L/radius << std::endl;
            }

            if (SAMRAI::tbox::SAMRAI_MPI::getRank() == 0)
            {
                drag_stream.close();
                lift_stream.close();
            }
        }
    }

    // Compute the updated force for visualization purposes.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_new_vec     = d_X_new_data [ln]->getGlobalVec();
            Vec X_current_vec = d_X_data     [ln]->getGlobalVec();
            Vec X_half_vec    = d_X_half_data[ln]->getGlobalVec();
            ierr = VecCopy(X_current_vec,X_half_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(X_half_vec,0.5,0.5,X_new_vec);  IBTK_CHKERRQ(ierr);
            Vec U_half_vec = d_U_half_data[ln]->getGlobalVec();
            ierr = VecCopy(X_current_vec,U_half_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(U_half_vec,1.0/d_dt,-1.0/d_dt,X_new_vec);  IBTK_CHKERRQ(ierr);

            Vec F_half_vec = d_F_half_data[ln]->getGlobalVec();
            ierr = VecSet(F_half_vec, 0.0);  IBTK_CHKERRQ(ierr);
            d_force_strategy->computeLagrangianForce(
                d_F_half_data[ln], d_X_half_data[ln], d_U_half_data[ln],
                d_hierarchy, ln, d_current_time+0.5*d_dt, d_lag_data_manager);
            resetAnchorPointValues(d_F_half_data, ln, ln);
        }
    }
    spread(d_f_new_idx, d_F_half_data, true, d_X_mid_data, false);

    // Deallocate the force Jacobian objects and preconditioner solvers.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_strct_mat[ln] != static_cast<Mat>(NULL))
        {
            ierr = MatDestroy(d_strct_mat[ln]);
            d_strct_mat[ln] = static_cast<Mat>(NULL);
        }
        if (d_strct_pc_mat[ln] != static_cast<Mat>(NULL))
        {
            ierr = MatDestroy(d_strct_pc_mat[ln]);
            d_strct_pc_mat[ln] = static_cast<Mat>(NULL);
        }
        if (d_strct_ksp[ln] != static_cast<KSP>(NULL))
        {
            ierr = KSPDestroy(d_strct_ksp[ln]);
            d_strct_ksp[ln] = static_cast<KSP>(NULL);
        }
        if (d_J_mat[ln] != static_cast<Mat>(NULL))
        {
            ierr = MatDestroy(d_J_mat[ln]);
            d_J_mat[ln] = static_cast<Mat>(NULL);
        }
    }

    // Deallocate the Jacobian object.
    ierr = PetscContainerDestroy(petsc_jac_ctx);  IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(petsc_jac);  IBTK_CHKERRQ(ierr);

    // Deallocate the PETSc solver.
    ierr = SNESDestroy(d_petsc_snes);  IBTK_CHKERRQ(ierr);

    // Deallocate the nullspace object.
    ierr = MatNullSpaceDestroy(petsc_nullsp);  IBTK_CHKERRQ(ierr);

    // Deallocate the multivectors.
    ierr = VecDestroy(petsc_sol_multivec);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(petsc_rhs_multivec);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(petsc_nul_multivec);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(petsc_res_multivec);  IBTK_CHKERRQ(ierr);

    // Deallocate the structure vectors.
    ierr = VecDestroy(petsc_structure_sol_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(petsc_structure_rhs_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(petsc_structure_nul_vec);  IBTK_CHKERRQ(ierr);

    // Deallocate the fluid vectors.
    IBTK::PETScSAMRAIVectorReal<double>::destroyPETScVector(petsc_fluid_sol_vec);
    IBTK::PETScSAMRAIVectorReal<double>::destroyPETScVector(petsc_fluid_rhs_vec);
    IBTK::PETScSAMRAIVectorReal<double>::destroyPETScVector(petsc_fluid_nul_vec);

    // Synchronize new state data.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_NEW_STATE_DATA"][ln]->coarsenData();
    }

    // Compute the cell-centered approximation to u^{n+1} (used for
    // visualization only).
    reinterpolateVelocity(getNewContext());

    // Compute the cell-centered approximation to f^{n+1} (used for
    // visualization only).
    reinterpolateForce(getNewContext());

    // Setup u_scratch to allow for ghost cell filling.
    d_hier_sc_data_ops->copyData(d_u_scratch_idx, d_u_new_idx);
    d_u_bdry_bc_fill_op->resetTransactionComponent(u_bc_component);
    d_u_bdry_bc_fill_op->fillData(new_time);

    // Compute omega = curl U.
    d_hier_math_ops->curl(
        d_omega_new_idx, d_omega_var,
        d_u_scratch_idx, d_u_var,
        d_no_fill_op, new_time);
#if (NDIM == 3)
    d_hier_math_ops->pointwise_L2Norm(
        d_omega_norm_new_idx, d_omega_norm_var,
        d_omega_new_idx, d_omega_var);
#endif

    // Compute max ||omega||_2.
#if (NDIM == 2)
    d_omega_max = std::max(+d_hier_cc_data_ops->max(d_omega_new_idx),
                           -d_hier_cc_data_ops->min(d_omega_new_idx));
#endif
#if (NDIM == 3)
    d_omega_max = d_hier_cc_data_ops->max(d_omega_norm_new_idx);
#endif
    // Compute Div U.
    d_hier_math_ops->div(
        d_div_u_new_idx, d_div_u_var,
        1.0, d_u_scratch_idx, d_u_var,
        d_no_fill_op, new_time, false);

    // Synchronize all data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): synchronizing data\n";
    synchronizeHierarchy();

    // Reset all time dependent data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierData(new_time);

    t_advance_hierarchy->stop();
    return getStableTimestep(getCurrentContext());
}// advanceHierarchy

double
IBStaggeredHierarchyIntegrator::getStableTimestep(
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx)
{
    t_get_stable_timestep->start();

    const bool initial_time = SAMRAI::tbox::MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    double dt_next = std::numeric_limits<double>::max();

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt_next = std::min(d_cfl*getLevelDt(level,ctx),dt_next);
    }

    if (d_integrator_time+dt_next >= d_end_time)
    {
        dt_next = d_end_time - d_integrator_time;
    }

    if (d_integrator_time >= d_dt_max_time_min &&
        d_integrator_time <= d_dt_max_time_max)
    {
        dt_next = std::min(d_dt_max,dt_next);
    }

    if (!initial_time)
    {
        dt_next = std::min(dt_next,d_grow_dt*d_old_dt);
    }

    t_get_stable_timestep->stop();
    return dt_next;
}// getStableTimestep()

bool
IBStaggeredHierarchyIntegrator::atRegridPoint() const
{
    const int level_number = 0;
    return ((d_integrator_step > 0)
            && d_gridding_alg->levelCanBeRefined(level_number)
            && (d_regrid_interval == 0
                ? false
                : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double
IBStaggeredHierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}// getIntegratorTime

double
IBStaggeredHierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}// getStartTime

double
IBStaggeredHierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}// getEndTime

int
IBStaggeredHierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}// getIntegratorStep

int
IBStaggeredHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}// getMaxIntegratorSteps

bool
IBStaggeredHierarchyIntegrator::stepsRemaining() const
{
    return (d_integrator_step < d_max_integrator_steps);
}// stepsRemaining

const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
IBStaggeredHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
IBStaggeredHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

SAMRAI::tbox::Pointer<HierarchyProjector>
IBStaggeredHierarchyIntegrator::getHierarchyProjector() const
{
    return d_hier_projector;
}// getHierarchyProjector

IBTK::LDataManager*
IBStaggeredHierarchyIntegrator::getLDataManager() const
{
    return d_lag_data_manager;
}// getLDataManager

///
///  The following routines:
///
///      regridHierarchy(),
///      synchronizeHierarchy(),
///      synchronizeNewLevels(),
///      resetTimeDependentHierData(),
///      resetHierDataToPreadvanceState()
///
///  allow the IBStaggeredHierarchyIntegrator to provide data management for a time
///  integrator which making use of this class.
///

void
IBStaggeredHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    // Update the workload pre-regridding.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_lag_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());

    // Before regriding, begin Lagrangian data movement.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): starting Lagrangian data movement.\n";
    d_lag_data_manager->beginDataRedistribution();

    // Regrid the hierarchy.
    const int coarsest_ln = 0;
    d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln, d_integrator_time, d_tag_buffer);
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Synchronize the state data on the patch hierarchy.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_CURRENT_STATE_DATA"][ln]->coarsenData();
    }

    // After regridding, finish Lagrangian data movement.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): finishing Lagrangian data movement.\n";
    d_lag_data_manager->endDataRedistribution();

    // Reset the various Lagrangian objects.  XXXX
    d_X_data     .resize(finest_ln+1);
    d_X_mid_data .resize(finest_ln+1);
    d_X_new_data .resize(finest_ln+1);
    d_X_half_data.resize(finest_ln+1);
    d_U_half_data.resize(finest_ln+1);
    d_F_half_data.resize(finest_ln+1);
    d_J_mat.resize(finest_ln+1,static_cast<Mat>(NULL));
    d_strct_mat.resize(finest_ln+1,static_cast<Mat>(NULL));
    d_strct_pc_mat.resize(finest_ln+1,static_cast<Mat>(NULL));
    d_strct_ksp.resize(finest_ln+1,static_cast<KSP>(NULL));
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(ln == finest_ln);
#endif
            d_X_data     [ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,ln);
            d_X_mid_data [ln] = d_lag_data_manager->createLNodeLevelData("X_mid" ,ln,NDIM);
            d_X_new_data [ln] = d_lag_data_manager->createLNodeLevelData("X_new" ,ln,NDIM);
            d_X_half_data[ln] = d_lag_data_manager->createLNodeLevelData("X_half",ln,NDIM);
            d_U_half_data[ln] = d_lag_data_manager->createLNodeLevelData("U_half",ln,NDIM);
            d_F_half_data[ln] = d_lag_data_manager->createLNodeLevelData("F_half",ln,NDIM);
        }
    }

    // Update the workload post-regridding.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_lag_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());

    // Indicate that the force strategy needs to be re-initialized.
    d_force_strategy_needs_init  = true;

    // Compute the set of local anchor points.
    static const double eps = 2.0*sqrt(std::numeric_limits<double>::epsilon());
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const grid_xLower = grid_geom->getXLower();
    const double* const grid_xUpper = grid_geom->getXUpper();
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        d_anchor_point_local_idxs[ln].clear();
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(lag_node_index_idx);
                for (IBTK::LNodeIndexData2::Iterator it(patch_box); it; it++)
                {
                    const SAMRAI::pdat::CellIndex<NDIM>& i = *it;
                    const IBTK::LNodeIndexSet& node_set = (*idx_data)(i);
                    for (IBTK::LNodeIndexSet::const_iterator n = node_set.begin();
                         n != node_set.end(); ++n)
                    {
                        const IBTK::LNodeIndexSet::value_type& node_idx = *n;
                        const std::vector<SAMRAI::tbox::Pointer<IBTK::Stashable> >& stash_data = node_idx->getStashData();
                        for (unsigned l = 0; l < stash_data.size(); ++l)
                        {
                            SAMRAI::tbox::Pointer<IBAnchorPointSpec> anchor_point_spec = stash_data[l];
                            if (!anchor_point_spec.isNull())
                            {
                                d_anchor_point_local_idxs[ln].insert(node_idx->getLocalPETScIndex());
                            }
                        }
                    }
                }
            }

            SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,ln);
            for (int i = 0; i < X_data->getLocalNodeCount(); ++i)
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    if ((periodic_shift[d] == 0) &&
                        ((*X_data)(i,d) - grid_xLower[d] <= eps ||
                         grid_xUpper[d] - (*X_data)(i,d) <= eps))
                    {
                        d_anchor_point_local_idxs[ln].insert(i);
                        break;
                    }
                }
            }
        }
    }

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

void
IBStaggeredHierarchyIntegrator::synchronizeHierarchy()
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

void
IBStaggeredHierarchyIntegrator::synchronizeNewLevels(
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
    // Synchronize initial data on the hierarchy.
    if (finest_level > 0)
    {
        for (int ln = finest_level; ln > coarsest_level; --ln)
        {
            d_cscheds["SYNCH_CURRENT_STATE_DATA"][ln]->coarsenData();
        }
    }

    t_synchronize_new_levels->stop();
    return;
}// synchronizeNewLevels

void
IBStaggeredHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    t_reset_time_dependent_data->start();

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // Swap SAMRAI::hier::PatchData<NDIM> pointers between the current and new contexts.
    for (std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > >::const_iterator sv =
             d_state_variables.begin();
         sv != d_state_variables.end(); ++sv)
    {
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& v = *sv;

        const int src_idx = var_db->mapVariableAndContextToIndex(
            v, getNewContext());
        const int dst_idx = var_db->mapVariableAndContextToIndex(
            v, getCurrentContext());

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > src_data =
                    patch->getPatchData(src_idx);
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > dst_data =
                    patch->getPatchData(dst_idx);

#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(src_data->getBox() == dst_data->getBox());
                TBOX_ASSERT(src_data->getGhostCellWidth() ==
                       dst_data->getGhostCellWidth());
#endif

                patch->setPatchData(dst_idx, src_data);
                patch->setPatchData(src_idx, dst_data);
            }

            if (d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                PetscErrorCode ierr;
                Vec X_vec     = d_X_data    [ln]->getGlobalVec();
                Vec X_new_vec = d_X_new_data[ln]->getGlobalVec();
                ierr = VecCopy(X_new_vec,X_vec);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    // Deallocate the scratch and new data and reset the time of the current
    // data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->setTime(new_time, d_current_data);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }

    t_reset_time_dependent_data->stop();
    return;
}// resetTimeDependentHierData

void
IBStaggeredHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    t_reset_data_to_preadvance_state->start();

    // Deallocate the scratch and new data and reset the time of the current
    // data.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->setTime(d_integrator_time, d_current_data);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }

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
IBStaggeredHierarchyIntegrator::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > base_old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > old_level = base_old_level;
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
    // Allocate storage needed to initialize the level and fill data from
    // coarser levels in AMR hierarchy, if any.
    //
    // Since time gets set when we allocate data, re-stamp it to current time if
    // we don't need to allocate.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    if (allocate_data)
    {
        level->allocatePatchData(d_current_data, init_data_time);
    }
    else
    {
        level->setTime(init_data_time, d_current_data);
    }

    // Fill data from coarser levels in AMR hierarchy.
    if (level_number > 0 || !old_level.isNull())
    {
        level->allocatePatchData(d_scratch_data, init_data_time);

        IBTK::CartExtrapPhysBdryOp fill_after_regrid_bc_op(d_fill_after_regrid_bc_idxs, BDRY_EXTRAP_TYPE);
        d_fill_after_regrid->createSchedule(level,
                                            old_level,
                                            level_number-1,
                                            hierarchy,
                                            &fill_after_regrid_bc_op)->fillData(init_data_time);
        level->deallocatePatchData(d_scratch_data);
    }

    // Initialize level data at the initial time.
    if (initial_time)
    {
        // If no initialization object is provided, initialize the velocity,
        // divergance, and vorticity to zero.  Otherwise, use the initialization
        // object to set the velocity to some specified value and compute the
        // divergance and vorticity corresponding to the initial velocity.
        if (d_u_init.isNull())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > u_current_data =
                    patch->getPatchData(d_u_current_idx);
                u_current_data->fillAll(0.0);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > u_cc_current_data =
                    patch->getPatchData(d_u_cc_current_idx);
                u_cc_current_data->fillAll(0.0);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > omega_current_data =
                    patch->getPatchData(d_omega_current_idx);
                omega_current_data->fillAll(0.0);
#if (NDIM == 3)
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > omega_norm_current_data =
                    patch->getPatchData(d_omega_norm_current_idx);
                omega_norm_current_data->fillAll(0.0);
#endif
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > div_u_current_data =
                    patch->getPatchData(d_div_u_current_idx);
                div_u_current_data->fillAll(0.0);
            }
        }
        else
        {
            level->allocatePatchData(d_u_scratch_idx, init_data_time);

            // Initialize U.
            d_u_init->setDataOnPatchLevel(
                d_u_current_idx, d_u_var, level,
                init_data_time, initial_time);
            IBTK::PatchMathOps patch_math_ops;
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > u_current_data =
                    patch->getPatchData(d_u_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > u_cc_current_data =
                    patch->getPatchData(d_u_cc_current_idx);

                patch_math_ops.interp(u_cc_current_data, u_current_data, patch);
            }

            // Fill in U boundary data from coarser levels.
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom =
                d_hierarchy->getGridGeometry();
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > ralg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator = grid_geom->lookupRefineOperator(
                d_u_var, "CONSERVATIVE_LINEAR_REFINE");
            ralg->registerRefine(d_u_scratch_idx, // destination
                                 d_u_current_idx, // source
                                 d_u_scratch_idx, // temporary work space
                                 refine_operator);
            IBTK::CartExtrapPhysBdryOp bc_op(d_u_scratch_idx, BDRY_EXTRAP_TYPE);
            ralg->createSchedule(level, level_number-1, hierarchy, &bc_op)->fillData(init_data_time);

            // Initialize quantities derived from the initial value of U.
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > u_scratch_data =
                    patch->getPatchData(d_u_scratch_idx);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > omega_current_data =
                    patch->getPatchData(d_omega_current_idx);
                patch_math_ops.curl(omega_current_data, u_scratch_data, patch);
#if (NDIM == 3)
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > omega_norm_current_data =
                    patch->getPatchData(d_omega_norm_current_idx);
                patch_math_ops.pointwise_L2Norm(omega_norm_current_data, omega_current_data, patch);
#endif
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > div_u_current_data =
                    patch->getPatchData(d_div_u_current_idx);
                patch_math_ops.div(div_u_current_data, 1.0, u_scratch_data, 0.0, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> >(NULL), patch);
            }

            level->deallocatePatchData(d_u_scratch_idx);
        }

        // Initialize the maximum value of ||omega||_2 on the grid.
        if (level_number == 0)
        {
            d_omega_max = 0.0;
        }

        SAMRAI::math::PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
#if (NDIM == 2)
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > omega_current_data =
                patch->getPatchData(d_omega_current_idx);
            d_omega_max = std::max(d_omega_max, +patch_cc_data_ops.max(omega_current_data, patch_box));
            d_omega_max = std::max(d_omega_max, -patch_cc_data_ops.min(omega_current_data, patch_box));
#endif
#if (NDIM == 3)
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > omega_norm_current_data =
                patch->getPatchData(d_omega_norm_current_idx);
            d_omega_max = std::max(d_omega_max, patch_cc_data_ops.max(omega_norm_current_data, patch_box));
#endif
        }

        // If no initialization object is provided, initialize the pressure to
        // zero.  Otherwise, use the initialization object to set the pressure
        // to some specified value.
        //
        // NOTE: This initial value for the pressure IS NOT USED by the time
        // integrator and is only specified for purposes of visualization.
        if (d_p_init.isNull())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > p_current_data =
                    patch->getPatchData(d_p_current_idx);
                p_current_data->fillAll(0.0);
            }
        }
        else
        {
            // Initialize P.
            d_p_init->setDataOnPatchLevel(
                d_p_current_idx, d_p_var, level,
                init_data_time, initial_time);
        }

        // If no initialization object is provided, initialize the body force to
        // zero.  Otherwise, use the initialization object to set the body force
        // to some specified value.
        if (d_f_set.isNull())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > f_current_data =
                    patch->getPatchData(d_f_current_idx);
                f_current_data->fillAll(0.0);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_cc_current_data =
                    patch->getPatchData(d_f_cc_current_idx);
                f_cc_current_data->fillAll(0.0);
            }
        }
        else
        {
            // Initialize F.
            d_f_set->setDataOnPatchLevel(
                d_f_current_idx, d_f_var, level,
                init_data_time, initial_time);
            IBTK::PatchMathOps patch_math_ops;
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > f_current_data =
                    patch->getPatchData(d_f_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_cc_current_data =
                    patch->getPatchData(d_f_cc_current_idx);

                patch_math_ops.interp(f_cc_current_data, f_current_data, patch);
            }

        }
    }

    // We use the LDataManager to handle unstructured data management.
    d_lag_data_manager->setPatchHierarchy(hierarchy);
    d_lag_data_manager->resetLevels(0,hierarchy->getFinestLevelNumber());
    d_lag_data_manager->initializeLevelData(
        hierarchy, level_number, init_data_time,
        can_be_refined, initial_time, old_level,
        allocate_data);

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBStaggeredHierarchyIntegrator::resetHierarchyConfiguration(
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

    // We use the LDataManager to handle unstructured data management.
    d_lag_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Reset the HierarchyProjector object.
    d_hier_projector->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->resetLevels(0, finest_hier_level);

    // Reset the Hierarchy math operations for the new configuration.
    if (d_is_managing_hier_math_ops)
    {
        d_hier_math_ops->setPatchHierarchy(hierarchy);
        d_hier_math_ops->resetLevels(0, finest_hier_level);
    }

    // Get the control volume weight variables and patch data descriptor
    // indices.
    d_wgt_cc_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    d_wgt_sc_var = d_hier_math_ops->getSideWeightVariable();
    d_wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();

    // Get the volume of the physical domain.
    d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    // Setup the patch boundary filling objects.
    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;

    InterpolationTransactionComponent u_bc_component(d_u_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_u_bc_coefs);
    d_u_bdry_bc_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_u_bdry_bc_fill_op->initializeOperatorState(u_bc_component, d_hierarchy);

    InterpolationTransactionComponent u_extrap_component(d_u_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);
    d_u_bdry_extrap_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_u_bdry_extrap_fill_op->initializeOperatorState(u_extrap_component, d_hierarchy);

    // Reset operators and solvers.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_sol_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::fluid_sol_vec", d_hierarchy, 0, finest_hier_level);
    fluid_sol_vec->addComponent(d_u_var,d_u_scratch_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    fluid_sol_vec->addComponent(d_p_var,d_p_scratch_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_rhs_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::fluid_rhs_vec", d_hierarchy, 0, finest_hier_level);
    fluid_rhs_vec->addComponent(d_u_var,d_u_rhs_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    fluid_rhs_vec->addComponent(d_p_var,d_p_rhs_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    d_convective_op->initializeOperatorState(*fluid_sol_vec,*fluid_rhs_vec);
    d_stokes_op->initializeOperatorState(*fluid_sol_vec,*fluid_rhs_vec);

    // If we have added or removed a level, resize the anchor point vectors.
    d_anchor_point_local_idxs.clear();
    d_anchor_point_local_idxs.resize(finest_hier_level+1);

    // If we have added or removed a level, resize the schedule vectors.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it!= d_ralgs.end(); ++it)
    {
        d_rscheds[(*it).first].resize(finest_hier_level+1);
    }

    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it!= d_calgs.end(); ++it)
    {
        d_cscheds[(*it).first].resize(finest_hier_level+1);
    }

    // (Re)build refine communication schedules.  These are created for all
    // levels in the hierarchy.
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

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it!= d_calgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level =
                hierarchy->getPatchLevel(ln-1);

            d_cscheds[(*it).first][ln] = (*it).second->
                createSchedule(coarser_level, level, d_cstrategies[(*it).first]);
        }
    }

    // Indicate that solvers need to be re-initialized.
    d_helmholtz_solver_needs_init = true;
    d_projection_pc_needs_init = true;

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
IBStaggeredHierarchyIntegrator::applyGradientDetector(
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

    // Tag cells which contain Lagrangian nodes.
    d_lag_data_manager->
        applyGradientDetector(hierarchy, level_number, error_data_time,
                              tag_index, initial_time,
                              uses_richardson_extrapolation_too);

    // Tag cells based on the magnatude of the vorticity.
    //
    // Note that if either the relative or absolute threshold is zero for a
    // particular level, no tagging is performed on that level.
    if (d_using_vorticity_tagging)
    {
        const double omega_rel_thresh =
            (level_number >= 0 && level_number < d_omega_rel_thresh.getSize()
             ? d_omega_rel_thresh[level_number]
             : (level_number < 0
                ? d_omega_rel_thresh[0]
                : d_omega_rel_thresh[d_omega_rel_thresh.size()-1]));
        const double omega_abs_thresh =
            (level_number >= 0 && level_number < d_omega_abs_thresh.getSize()
             ? d_omega_abs_thresh[level_number]
             : (level_number < 0
                ? d_omega_abs_thresh[0]
                : d_omega_abs_thresh[d_omega_abs_thresh.size()-1]));
        if (omega_rel_thresh > 0.0 && omega_abs_thresh > 0.0)
        {
            const double thresh = sqrt(std::numeric_limits<double>::epsilon()) +
                std::min(omega_rel_thresh*d_omega_max, omega_abs_thresh);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > omega_current_data = patch->getPatchData(d_omega_current_idx);
                for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = ic();
#if (NDIM == 2)
                    if (std::abs((*omega_current_data)(i)) > thresh)
                    {
                        (*tags_data)(i) = 1;
                    }
#endif
#if (NDIM == 3)
                    double norm_omega_sq = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        norm_omega_sq += (*omega_current_data)(i,d)*(*omega_current_data)(i,d);
                    }
                    const double norm_omega = sqrt(norm_omega_sq);
                    if (norm_omega > thresh)
                    {
                        (*tags_data)(i) = 1;
                    }
#endif
                }
            }
        }
    }

    t_apply_gradient_detector->stop();
    return;
}// applyGradientDetector

///
///  The following routines:
///
///      getVelocityVar(),
///      getPressureVar(),
///      getForceVar()
///
///  allows access to the various state variables maintained by the integrator.
///

SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >
IBStaggeredHierarchyIntegrator::getVelocityVar()
{
    return d_u_var;
}// getVelocityVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
IBStaggeredHierarchyIntegrator::getPressureVar()
{
    return d_p_var;
}// getPressureVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >
IBStaggeredHierarchyIntegrator::getForceVar()
{
    return d_f_var;
}// getForceVar

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
/// AdvDiffIntegrator object.
///

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBStaggeredHierarchyIntegrator::getCurrentContext() const
{
    return d_current_context;
}// getCurrentContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBStaggeredHierarchyIntegrator::getNewContext() const
{
    return d_new_context;
}// getNewContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBStaggeredHierarchyIntegrator::getScratchContext() const
{
    return d_scratch_context;
}// getScratchContext

///
/// The following routines:
///
///      reinterpolateVelocity(),
///      reinterpolateForce()
///
/// are miscelaneous utility functions.

void
IBStaggeredHierarchyIntegrator::reinterpolateVelocity(
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx)
{
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int    u_idx = var_db->mapVariableAndContextToIndex(   d_u_var, ctx);
    const int u_cc_idx = var_db->mapVariableAndContextToIndex(d_u_cc_var, ctx);
    static const bool synch_cf_interface = true;
    d_hier_math_ops->interp(
        u_cc_idx, d_u_cc_var,
        u_idx   , d_u_var   ,
        d_no_fill_op, d_integrator_time, synch_cf_interface);
    return;
}// reinterpolateVelocity

void
IBStaggeredHierarchyIntegrator::reinterpolateForce(
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx)
{
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int    f_idx = var_db->mapVariableAndContextToIndex(   d_f_var, ctx);
    const int f_cc_idx = var_db->mapVariableAndContextToIndex(d_f_cc_var, ctx);
    static const bool synch_cf_interface = true;
    d_hier_math_ops->interp(
        f_cc_idx, d_f_cc_var,
        f_idx   , d_f_var   ,
        d_no_fill_op, d_integrator_time, synch_cf_interface);
    return;
}// reinterpolateForce

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::tbox::Serializable abstract base class.
///

void
IBStaggeredHierarchyIntegrator::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    TBOX_ASSERT(false);

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
IBStaggeredHierarchyIntegrator::printClassData(
    std::ostream& os) const
{
    TBOX_ASSERT(false);
    return;
}// printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBStaggeredHierarchyIntegrator::registerVariable(
    int& current_idx,
    int& new_idx,
    int& scratch_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
    const SAMRAI::hier::IntVector<NDIM>& scratch_ghosts,
    const std::string& coarsen_name,
    const std::string& refine_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!variable.isNull());
#endif
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    current_idx = -1; // insure that uninitialized variable patch data
    new_idx     = -1; // descriptor indices cause errors
    scratch_idx = -1;

    d_state_variables.push_back(variable);

    // Setup the current context.
    current_idx = var_db->registerVariableAndContext(
        variable, getCurrentContext(), no_ghosts);
    d_current_data.setFlag(current_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(current_idx);
    }

    // Setup the new context.
    new_idx = var_db->registerVariableAndContext(
        variable, getNewContext(), no_ghosts);
    d_new_data.setFlag(new_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(new_idx);
    }

    // Setup the scratch context.
    scratch_idx = var_db->registerVariableAndContext(
        variable, getScratchContext(), scratch_ghosts);
    d_scratch_data.setFlag(scratch_idx);

    // Setup the communication algorithms.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator =
        grid_geom->lookupRefineOperator(variable, refine_name);
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_operator =
        grid_geom->lookupCoarsenOperator(variable, coarsen_name);

    // Setup the refine algorithm used to fill data in new or modified patch
    // levels following a regrid operation.
    if (!refine_operator.isNull())
    {
        d_fill_after_regrid->
            registerRefine(current_idx, // destination
                           current_idx, // source
                           scratch_idx, // temporary work space
                           refine_operator);

        // Keep track of the cell-centered scratch data indices, for use in the
        // refinement schedule used to fill data in new or modified patch levels
        // following a regrid operation.
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var = variable;
        if (!cc_var.isNull())
        {
            d_fill_after_regrid_bc_idxs.setFlag(scratch_idx);
        }
    }

    // Setup the SYNCH_CURRENT_STATE_DATA and SYNCH_NEW_STATE_DATA algorithms,
    // used to synchronize the data on the hierarchy.
    if (!coarsen_operator.isNull())
    {
        d_calgs["SYNCH_CURRENT_STATE_DATA"]->
            registerCoarsen(current_idx, // destination
                            current_idx, // source
                            coarsen_operator);

        d_calgs["SYNCH_NEW_STATE_DATA"]->
            registerCoarsen(new_idx, // destination
                            new_idx, // source
                            coarsen_operator);
    }
    return;
}// registerVariable

void
IBStaggeredHierarchyIntegrator::registerVariable(
    int& scratch_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
    const SAMRAI::hier::IntVector<NDIM>& scratch_ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!variable.isNull());
#endif

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    scratch_idx = -1; // insure that uninitialized variable patch data
                      // descriptor indices cause errors

    d_scratch_variables.push_back(variable);

    // Setup the scratch context.
    scratch_idx = var_db->registerVariableAndContext(
        variable, getScratchContext(), scratch_ghosts);
    d_scratch_data.setFlag(scratch_idx);
    return;
}// registerVariable

/////////////////////////////// PRIVATE //////////////////////////////////////

#undef __FUNCT__
#define __FUNCT__ "IBStaggeredHierarchyIntegrator::FormFunction_SAMRAI"
PetscErrorCode
IBStaggeredHierarchyIntegrator::FormFunction_SAMRAI(
    SNES snes,
    Vec x,
    Vec f,
    void* p_ctx)
{
    PetscErrorCode ierr;
    IBStaggeredHierarchyIntegrator* hier_integrator = static_cast<IBStaggeredHierarchyIntegrator*>(p_ctx);
#if DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(hier_integrator != NULL);
#endif
    hier_integrator->FormFunction(x,f);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f));  IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// FormFunction_SAMRAI

void
IBStaggeredHierarchyIntegrator::FormFunction(
    Vec x,
    Vec f)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    PetscErrorCode ierr;

#ifdef DEBUG_CHECK_ASSERTIONS
    int N_x;
    ierr = IBTK::VecMultiVecGetNumberOfVecs(x,&N_x);  IBTK_CHKERRQ(ierr);
    int N_f;
    ierr = IBTK::VecMultiVecGetNumberOfVecs(f,&N_f);  IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(N_x == 2);
    TBOX_ASSERT(N_f == 2);
#endif
    Vec* petsc_x_multivecs;
    ierr = IBTK::VecMultiVecGetVecs(x, &petsc_x_multivecs);  IBTK_CHKERRQ(ierr);
    Vec petsc_fluid_x_vec = petsc_x_multivecs[0];
    Vec petsc_strct_x_vec = petsc_x_multivecs[1];

    Vec* petsc_f_multivecs;
    ierr = IBTK::VecMultiVecGetVecs(f, &petsc_f_multivecs);  IBTK_CHKERRQ(ierr);
    Vec petsc_fluid_f_vec = petsc_f_multivecs[0];
    Vec petsc_strct_f_vec = petsc_f_multivecs[1];

    // The fluid function evaluation.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_x_vec = IBTK::PETScSAMRAIVectorReal<double>::getSAMRAIVector(petsc_fluid_x_vec);
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_f_vec = IBTK::PETScSAMRAIVectorReal<double>::getSAMRAIVector(petsc_fluid_f_vec);
    static const bool homogeneous_bc = false;
    d_stokes_op->apply(homogeneous_bc, *fluid_x_vec, *fluid_f_vec);

    // The structure function evaluation.
    Vec X_new_vec     = petsc_strct_x_vec;
    Vec X_current_vec = d_X_data[finest_ln]->getGlobalVec();
    Vec X_half_vec    = d_X_half_data[finest_ln]->getGlobalVec();
    ierr = VecCopy(X_current_vec, X_half_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAXPBY(X_half_vec, 0.5, 0.5, X_new_vec);  IBTK_CHKERRQ(ierr);
    Vec U_half_vec = d_U_half_data[finest_ln]->getGlobalVec();
    ierr = VecCopy(X_current_vec,U_half_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAXPBY(U_half_vec,1.0/d_dt,-1.0/d_dt,X_new_vec);  IBTK_CHKERRQ(ierr);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            Vec F_half_vec = d_F_half_data[ln]->getGlobalVec();
            PetscErrorCode ierr = VecSet(F_half_vec, 0.0);  IBTK_CHKERRQ(ierr);
            d_force_strategy->computeLagrangianForce(
                d_F_half_data[ln], d_X_half_data[ln], d_U_half_data[ln],
                d_hierarchy, ln, d_current_time+0.5*d_dt, d_lag_data_manager);
            resetAnchorPointValues(d_F_half_data, ln, ln);
        }
    }
    spread(d_f_scratch_idx, d_F_half_data, true, d_X_mid_data, false);

    d_hier_sc_data_ops->axpy(fluid_f_vec->getComponentDescriptorIndex(0), -1.0, d_f_scratch_idx, fluid_f_vec->getComponentDescriptorIndex(0));
    d_u_bc_helper->zeroValuesAtDirichletBoundaries(fluid_f_vec->getComponentDescriptorIndex(0));

    d_hier_sc_data_ops->linearSum(d_u_interp_idx, 0.5, fluid_x_vec->getComponentDescriptorIndex(0), 0.5, d_u_current_idx);
    interp(d_U_half_data, d_u_interp_idx, true, d_X_mid_data, false);
    resetAnchorPointValues(d_U_half_data, coarsest_ln, finest_ln);

    Vec R_vec = petsc_strct_f_vec;
    ierr = VecWAXPY(R_vec,-1.0,X_current_vec,X_new_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAXPBY(R_vec,-1.0,(1.0/d_dt),U_half_vec);  IBTK_CHKERRQ(ierr);

    // Flush any cached data associated with the modified PETSc vectors.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_fluid_f_vec));  IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_strct_f_vec));  IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f));  IBTK_CHKERRQ(ierr);
    return;
}// FormFunction

#undef __FUNCT__
#define __FUNCT__ "IBStaggeredHierarchyIntegrator::FormJacobian_SAMRAI"
PetscErrorCode
IBStaggeredHierarchyIntegrator::FormJacobian_SAMRAI(
    SNES snes,
    Vec x,
    Mat* A,
    Mat* B,
    MatStructure* mat_structure,
    void* p_ctx)
{
    IBStaggeredHierarchyIntegrator* hier_integrator = static_cast<IBStaggeredHierarchyIntegrator*>(p_ctx);
#if DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(hier_integrator != NULL);
#endif
    hier_integrator->FormJacobian(x);
    if (!hier_integrator->d_use_exact_jacobian)
    {
        PetscErrorCode ierr = MatMFFDComputeJacobian(snes,x,A,B,mat_structure,p_ctx);  CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}// FormJacobian_SAMRAI

void
IBStaggeredHierarchyIntegrator::FormJacobian(
    Vec x)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    PetscErrorCode ierr;

    Vec* petsc_x_multivecs;
#ifdef DEBUG_CHECK_ASSERTIONS
    int N_x;
    ierr = IBTK::VecMultiVecGetNumberOfVecs(x,&N_x);  IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(N_x == 2);
#endif
    ierr = IBTK::VecMultiVecGetVecs(x, &petsc_x_multivecs);  IBTK_CHKERRQ(ierr);
    Vec petsc_fluid_x_vec = petsc_x_multivecs[0];
    Vec petsc_strct_x_vec = petsc_x_multivecs[1];
    (void) petsc_fluid_x_vec;

    // Destroy any old Jacobian data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_strct_mat[ln] != static_cast<Mat>(NULL))
        {
            ierr = MatDestroy(d_strct_mat[ln]);
            d_strct_mat[ln] = static_cast<Mat>(NULL);
        }
        if (d_strct_pc_mat[ln] != static_cast<Mat>(NULL))
        {
            ierr = MatDestroy(d_strct_pc_mat[ln]);
            d_strct_pc_mat[ln] = static_cast<Mat>(NULL);
        }
        if (d_strct_ksp[ln] != static_cast<KSP>(NULL))
        {
            ierr = KSPDestroy(d_strct_ksp[ln]);
            d_strct_ksp[ln] = static_cast<KSP>(NULL);
        }
        if (d_J_mat[ln] != static_cast<Mat>(NULL))
        {
            ierr = MatDestroy(d_J_mat[ln]);
            d_J_mat[ln] = static_cast<Mat>(NULL);
        }
    }

    // Compute J = dF/dX at the current approximation to X_half and U_half.
    Vec X_new_vec     = petsc_strct_x_vec;
    Vec X_current_vec = d_X_data[finest_ln]->getGlobalVec();
    Vec X_half_vec    = d_X_half_data[finest_ln]->getGlobalVec();
    ierr = VecCopy(X_current_vec, X_half_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAXPBY(X_half_vec, 0.5, 0.5, X_new_vec);  IBTK_CHKERRQ(ierr);
    Vec U_half_vec = d_U_half_data[finest_ln]->getGlobalVec();
    ierr = VecCopy(X_current_vec,U_half_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAXPBY(U_half_vec,1.0/d_dt,-1.0/d_dt,X_new_vec);  IBTK_CHKERRQ(ierr);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            // Determine the non-zero structure of the force Jacobian matrix and
            // allocate a block AIJ matrix.
            const int num_local_nodes = d_lag_data_manager->getNumberOfLocalNodes(ln);
            std::vector<int> d_nnz(num_local_nodes), o_nnz(num_local_nodes);
            d_force_strategy->computeLagrangianForceJacobianNonzeroStructure(
                d_nnz, o_nnz, d_hierarchy, ln, d_current_time+0.5*d_dt, d_lag_data_manager);
            ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                                    NDIM, NDIM*num_local_nodes, NDIM*num_local_nodes,
                                    PETSC_DETERMINE, PETSC_DETERMINE,
                                    PETSC_DEFAULT, &d_nnz[0],
                                    PETSC_DEFAULT, &o_nnz[0],
                                    &d_J_mat[ln]);  IBTK_CHKERRQ(ierr);

            // Compute the Jacobian of the force.
            d_force_strategy->computeLagrangianForceJacobian(
                d_J_mat[ln], MAT_FINAL_ASSEMBLY, 0.5, d_X_half_data[ln], 1.0/d_dt, d_U_half_data[ln],
                d_hierarchy, ln, d_current_time+0.5*d_dt, d_lag_data_manager);

            // Reset any rows corresponding to anchor nodes.
            const int num_zeroed_rows = d_anchor_point_local_idxs[ln].size();
            std::vector<int> zeroed_row_idxs;
            const int global_node_offset = d_lag_data_manager->getGlobalNodeOffset(ln);
            for (std::set<int>::const_iterator cit = d_anchor_point_local_idxs[ln].begin();
                 cit != d_anchor_point_local_idxs[ln].end(); ++cit)
            {
                const int& i = *cit;
                zeroed_row_idxs.push_back(i+global_node_offset);
            }
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(num_zeroed_rows == int(zeroed_row_idxs.size()));
#endif
            ierr = MatZeroRows(d_J_mat[ln], num_zeroed_rows, &zeroed_row_idxs[0], 0.0);  IBTK_CHKERRQ(ierr);
        }
    }

    // Setup a shell matrix for the approximate structure Jacobian matrix, i.e.,
    // I - (dt^2)/(2*rho) S S^{*} dF/dX(n+1), and form a preconditioner matrix
    // obtained by replacing S S^{*} with diag(S S^{*}).
    int local_sz;
    ierr = VecGetLocalSize(petsc_strct_x_vec, &local_sz);  IBTK_CHKERRQ(ierr);
    static const double C = pow(IBTK::LEInteractor::getC(d_delta_fcn),NDIM);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            ierr = MatCreateShell(PETSC_COMM_WORLD, local_sz, local_sz, PETSC_DETERMINE, PETSC_DETERMINE, static_cast<void*>(this), &d_strct_mat[ln]);  IBTK_CHKERRQ(ierr);
            ierr = MatShellSetOperation(d_strct_mat[ln], MATOP_MULT, reinterpret_cast<void(*)(void)>(IBStaggeredHierarchyIntegrator::MatVecMult_structure_SAMRAI));  IBTK_CHKERRQ(ierr);

            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            const double* const dx_coarsest = grid_geom->getDx();
            const SAMRAI::hier::IntVector<NDIM> ratio = level->getRatio();
            double volume_element = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                volume_element *= dx_coarsest[d]/double(ratio[d]);
            }
            ierr = MatDuplicate(d_J_mat[ln], MAT_COPY_VALUES, &d_strct_pc_mat[ln]);  IBTK_CHKERRQ(ierr);
//          ierr = MatConvert(d_J_mat[ln], MATMPIAIJ, MAT_INITIAL_MATRIX, &d_strct_pc_mat[ln]);  IBTK_CHKERRQ(ierr);
            ierr = MatScale(d_strct_pc_mat[ln], -0.5*d_dt*d_dt*(C/volume_element)/d_rho);  IBTK_CHKERRQ(ierr);
            ierr = MatShift(d_strct_pc_mat[ln], 1.0);  IBTK_CHKERRQ(ierr);

            static const std::string options_prefix = "strct_";
            ierr = KSPCreate(PETSC_COMM_WORLD, &d_strct_ksp[ln]);  IBTK_CHKERRQ(ierr);
            ierr = KSPSetOperators(d_strct_ksp[ln], d_strct_mat[ln], d_strct_pc_mat[ln], SAME_PRECONDITIONER);  IBTK_CHKERRQ(ierr);
            ierr = KSPSetOptionsPrefix(d_strct_ksp[ln], options_prefix.c_str());  IBTK_CHKERRQ(ierr);
            ierr = KSPSetFromOptions(d_strct_ksp[ln]);  IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// FormJacobian

#undef __FUNCT__
#define __FUNCT__ "IBStaggeredHierarchyIntegrator::MatVecMult_SAMRAI"
PetscErrorCode
IBStaggeredHierarchyIntegrator::MatVecMult_SAMRAI(
    Mat A,
    Vec x,
    Vec y)
{
    PetscErrorCode ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx);  IBTK_CHKERRQ(ierr);
    IBStaggeredHierarchyIntegrator* hier_integrator = static_cast<IBStaggeredHierarchyIntegrator*>(p_ctx);
#if DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(hier_integrator != NULL);
#endif
    hier_integrator->MatVecMult(x,y);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));  IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// MatVecMult_SAMRAI

void
IBStaggeredHierarchyIntegrator::MatVecMult(
    Vec x,
    Vec y)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    PetscErrorCode ierr;
#ifdef DEBUG_CHECK_ASSERTIONS
    int N_x;
    ierr = IBTK::VecMultiVecGetNumberOfVecs(x,&N_x);  IBTK_CHKERRQ(ierr);
    int N_y;
    ierr = IBTK::VecMultiVecGetNumberOfVecs(y,&N_y);  IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(N_x == 2);
    TBOX_ASSERT(N_y == 2);
#endif
    Vec* petsc_x_multivecs;
    ierr = IBTK::VecMultiVecGetVecs(x, &petsc_x_multivecs);  IBTK_CHKERRQ(ierr);
    Vec petsc_fluid_x_vec = petsc_x_multivecs[0];
    Vec petsc_strct_x_vec = petsc_x_multivecs[1];

    Vec* petsc_y_multivecs;
    ierr = IBTK::VecMultiVecGetVecs(y, &petsc_y_multivecs);  IBTK_CHKERRQ(ierr);
    Vec petsc_fluid_y_vec = petsc_y_multivecs[0];
    Vec petsc_strct_y_vec = petsc_y_multivecs[1];

    // The fluid Jacobian application.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_x_vec = IBTK::PETScSAMRAIVectorReal<double>::getSAMRAIVector(petsc_fluid_x_vec);
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_y_vec = IBTK::PETScSAMRAIVectorReal<double>::getSAMRAIVector(petsc_fluid_y_vec);
    d_stokes_op->apply(*fluid_x_vec, *fluid_y_vec);

    // The structure Jacobian application.
    d_hier_sc_data_ops->scale(d_u_interp_idx, 0.5, fluid_x_vec->getComponentDescriptorIndex(0));
    d_u_bc_helper->zeroValuesAtDirichletBoundaries(d_u_interp_idx);
    interp(d_U_half_data, d_u_interp_idx, true, d_X_mid_data, false);
    resetAnchorPointValues(d_U_half_data, coarsest_ln, finest_ln);

    Vec X_new_vec  = petsc_strct_x_vec;
    Vec F_half_vec = d_F_half_data[finest_ln]->getGlobalVec();
    ierr = MatMult(d_J_mat[finest_ln], X_new_vec, F_half_vec);  IBTK_CHKERRQ(ierr);
    spread(d_f_scratch_idx, d_F_half_data, true, d_X_mid_data, false);
    d_hier_sc_data_ops->axpy(fluid_y_vec->getComponentDescriptorIndex(0), -1.0, d_f_scratch_idx, fluid_y_vec->getComponentDescriptorIndex(0));
    d_u_bc_helper->zeroValuesAtDirichletBoundaries(fluid_y_vec->getComponentDescriptorIndex(0));

    Vec U_half_vec = d_U_half_data[finest_ln]->getGlobalVec();
    Vec R_vec      = petsc_strct_y_vec;
    ierr = VecCopy(X_new_vec,R_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAXPBY(R_vec,-1.0,(1.0/d_dt),U_half_vec);  IBTK_CHKERRQ(ierr);

    // Flush any cached data associated with the modified PETSc vectors.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_fluid_y_vec));  IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_strct_y_vec));  IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));  IBTK_CHKERRQ(ierr);
    return;
}// MatVecMult

#undef __FUNCT__
#define __FUNCT__ "IBStaggeredHierarchyIntegrator::MatGetVecs_SAMRAI"
PetscErrorCode
IBStaggeredHierarchyIntegrator::MatGetVecs_SAMRAI(
    Mat A,
    Vec* right,
    Vec* left)
{
    PetscErrorCode ierr;

    PetscContainer petsc_jac_ctx;
    ierr = PetscObjectQuery(reinterpret_cast<PetscObject>(A), "petsc_jac_ctx", reinterpret_cast<PetscObject*>(&petsc_jac_ctx));  IBTK_CHKERRQ(ierr);
    void* p_ctx;
    ierr = PetscContainerGetPointer(petsc_jac_ctx, &p_ctx);  IBTK_CHKERRQ(ierr);
    IBStaggeredHierarchyIntegrator* hier_integrator = static_cast<IBStaggeredHierarchyIntegrator*>(p_ctx);
#if DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(hier_integrator != NULL);
#endif
    if (right != PETSC_NULL)
    {
        // vector that the matrix can be multiplied against
        ierr = VecDuplicate(hier_integrator->d_petsc_x_vec, right); IBTK_CHKERRQ(ierr);
        ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*right)); IBTK_CHKERRQ(ierr);
    }
    if (left != PETSC_NULL)
    {
        // vector that the matrix vector product can be stored in
        ierr = VecDuplicate(hier_integrator->d_petsc_f_vec, left); IBTK_CHKERRQ(ierr);
        ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*left)); IBTK_CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}// MatGetVecs_SAMRAI

#undef __FUNCT__
#define __FUNCT__ "IBStaggeredHierarchyIntegrator::PCApplyFluid_SAMRAI"
PetscErrorCode
IBStaggeredHierarchyIntegrator::PCApplyFluid_SAMRAI(
    void* p_ctx,
    Vec x,
    Vec y)
{
    PetscErrorCode ierr;
    IBStaggeredHierarchyIntegrator* hier_integrator = static_cast<IBStaggeredHierarchyIntegrator*>(p_ctx);
#if DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(hier_integrator != NULL);
#endif
    hier_integrator->PCApplyFluid(x,y);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));  IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// PCApplyFluid_SAMRAI

void
IBStaggeredHierarchyIntegrator::PCApplyFluid(
    Vec x,
    Vec y)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    PetscErrorCode ierr;
#ifdef DEBUG_CHECK_ASSERTIONS
    int N_x;
    ierr = IBTK::VecMultiVecGetNumberOfVecs(x,&N_x);  IBTK_CHKERRQ(ierr);
    int N_y;
    ierr = IBTK::VecMultiVecGetNumberOfVecs(y,&N_y);  IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(N_x == 2);
    TBOX_ASSERT(N_y == 2);
#endif
    Vec* petsc_x_multivecs;
    ierr = IBTK::VecMultiVecGetVecs(x, &petsc_x_multivecs);  IBTK_CHKERRQ(ierr);
    Vec petsc_fluid_x_vec = petsc_x_multivecs[0];
    Vec petsc_strct_x_vec = petsc_x_multivecs[1];

    Vec* petsc_y_multivecs;
    ierr = IBTK::VecMultiVecGetVecs(y, &petsc_y_multivecs);  IBTK_CHKERRQ(ierr);
    Vec petsc_fluid_y_vec = petsc_y_multivecs[0];
    Vec petsc_strct_y_vec = petsc_y_multivecs[1];

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_x_vec = IBTK::PETScSAMRAIVectorReal<double>::getSAMRAIVector(petsc_fluid_x_vec);
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_y_vec = IBTK::PETScSAMRAIVectorReal<double>::getSAMRAIVector(petsc_fluid_y_vec);

    // Compute the action of the fluid preconditioner.
    d_projection_pc->solveSystem(*fluid_y_vec, *fluid_x_vec);
    d_u_bc_helper->zeroValuesAtDirichletBoundaries(fluid_y_vec->getComponentDescriptorIndex(0));

    // Update the structure data appropriately.
    d_hier_sc_data_ops->copyData(d_u_interp_idx, fluid_y_vec->getComponentDescriptorIndex(0));
    interp(d_U_half_data, d_u_interp_idx, true, d_X_mid_data, false);
    resetAnchorPointValues(d_U_half_data, coarsest_ln, finest_ln);
    ierr = VecWAXPY(petsc_strct_y_vec, 0.5, d_U_half_data[finest_ln]->getGlobalVec(), petsc_strct_x_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecScale(petsc_strct_y_vec, d_dt);  IBTK_CHKERRQ(ierr);

    // Flush any cached data associated with the modified PETSc vectors.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_fluid_y_vec));  IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_strct_y_vec));  IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));  IBTK_CHKERRQ(ierr);
    return;
}// PCApplyFluid

int
IBStaggeredHierarchyIntegrator::PCApplyStrct_SAMRAI(
    void* p_ctx,
    Vec x,
    Vec y)
{
    PetscErrorCode ierr;
    IBStaggeredHierarchyIntegrator* hier_integrator = static_cast<IBStaggeredHierarchyIntegrator*>(p_ctx);
#if DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(hier_integrator != NULL);
#endif
    hier_integrator->PCApplyStrct(x,y);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));  IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// PCApplyStrct_SAMRAI

void
IBStaggeredHierarchyIntegrator::PCApplyStrct(
    Vec x,
    Vec y)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    PetscErrorCode ierr;
#ifdef DEBUG_CHECK_ASSERTIONS
    int N_x;
    ierr = IBTK::VecMultiVecGetNumberOfVecs(x,&N_x);  IBTK_CHKERRQ(ierr);
    int N_y;
    ierr = IBTK::VecMultiVecGetNumberOfVecs(y,&N_y);  IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(N_x == 2);
    TBOX_ASSERT(N_y == 2);
#endif
    Vec* petsc_x_multivecs;
    ierr = IBTK::VecMultiVecGetVecs(x, &petsc_x_multivecs);  IBTK_CHKERRQ(ierr);
    Vec petsc_fluid_x_vec = petsc_x_multivecs[0];
    Vec petsc_strct_x_vec = petsc_x_multivecs[1];

    Vec* petsc_y_multivecs;
    ierr = IBTK::VecMultiVecGetVecs(y, &petsc_y_multivecs);  IBTK_CHKERRQ(ierr);
    Vec petsc_fluid_y_vec = petsc_y_multivecs[0];
    Vec petsc_strct_y_vec = petsc_y_multivecs[1];

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_x_vec = IBTK::PETScSAMRAIVectorReal<double>::getSAMRAIVector(petsc_fluid_x_vec);
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > fluid_y_vec = IBTK::PETScSAMRAIVectorReal<double>::getSAMRAIVector(petsc_fluid_y_vec);

    // Initialze the result to equal zero.
    ierr = VecSet(y, 0.0);  IBTK_CHKERRQ(ierr);

    // Allocate temporary data.
    Vec petsc_strct_f_vec;
    ierr = VecDuplicate(petsc_strct_x_vec, &petsc_strct_f_vec);  IBTK_CHKERRQ(ierr);

    // Compute the right hand side term.
    ierr = VecCopy(petsc_strct_x_vec, petsc_strct_f_vec);  IBTK_CHKERRQ(ierr);
    d_hier_sc_data_ops->copyData(d_u_interp_idx, fluid_x_vec->getComponentDescriptorIndex(0));
    d_u_bc_helper->zeroValuesAtDirichletBoundaries(d_u_interp_idx);
    interp(d_U_half_data, d_u_interp_idx, true, d_X_mid_data, false);
    resetAnchorPointValues(d_U_half_data, coarsest_ln, finest_ln);
    ierr = VecAXPBY(petsc_strct_f_vec, 0.5*d_dt*d_dt/d_rho, d_dt, d_U_half_data[finest_ln]->getGlobalVec());  IBTK_CHKERRQ(ierr);

    // Compute the application of the structure preconditioner.
    ierr = KSPSolve(d_strct_ksp[finest_ln], petsc_strct_f_vec, petsc_strct_y_vec);  IBTK_CHKERRQ(ierr);

    // Compute the fluid component.
    ierr = MatMult(d_J_mat[finest_ln], petsc_strct_y_vec, d_F_half_data[finest_ln]->getGlobalVec());  IBTK_CHKERRQ(ierr);
    spread(d_f_scratch_idx, d_F_half_data, true, d_X_mid_data, false);
    d_hier_sc_data_ops->linearSum(fluid_y_vec->getComponentDescriptorIndex(0), d_dt/d_rho, fluid_x_vec->getComponentDescriptorIndex(0), d_dt/d_rho, d_f_scratch_idx);
    d_u_bc_helper->zeroValuesAtDirichletBoundaries(fluid_y_vec->getComponentDescriptorIndex(0));

    // Deallocate temporary data.
    ierr = VecDestroy(petsc_strct_f_vec);  IBTK_CHKERRQ(ierr);

    // Flush any cached data associated with the modified PETSc vectors.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_fluid_y_vec));  IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(petsc_strct_y_vec));  IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));  IBTK_CHKERRQ(ierr);
    return;
}// PCApplyStrct

#undef __FUNCT__
#define __FUNCT__ "IBStaggeredHierarchyIntegrator::MatVecMult_structure_SAMRAI"
PetscErrorCode
IBStaggeredHierarchyIntegrator::MatVecMult_structure_SAMRAI(
    Mat A,
    Vec x,
    Vec y)
{
    PetscErrorCode ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx);  IBTK_CHKERRQ(ierr);
    IBStaggeredHierarchyIntegrator* hier_integrator = static_cast<IBStaggeredHierarchyIntegrator*>(p_ctx);
#if DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(hier_integrator != NULL);
#endif
    hier_integrator->MatVecMult_structure(x,y);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));  IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// MatVecMult_structure_SAMRAI

void
IBStaggeredHierarchyIntegrator::MatVecMult_structure(
    Vec x,
    Vec y)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    PetscErrorCode ierr;
    ierr = MatMult(d_J_mat[finest_ln], x, d_F_half_data[finest_ln]->getGlobalVec());  IBTK_CHKERRQ(ierr);
    resetAnchorPointValues(d_F_half_data, coarsest_ln, finest_ln);
    spread(d_u_interp_idx, d_F_half_data, true, d_X_mid_data, false);
    d_u_bc_helper->zeroValuesAtDirichletBoundaries(d_u_interp_idx);
    interp(d_U_half_data, d_u_interp_idx, true, d_X_mid_data, false);
    ierr = VecWAXPY(y, -0.5*d_dt*d_dt/d_rho, d_U_half_data[finest_ln]->getGlobalVec(), x);  IBTK_CHKERRQ(ierr);
    return;
}// MatVecMult_structure

void
IBStaggeredHierarchyIntegrator::spread(
    const int f_data_idx,
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_data,
    const bool F_data_ghost_node_update,
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data,
    const bool X_data_ghost_node_update,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? d_hierarchy->getFinestLevelNumber() : finest_ln_in);
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            if (F_data_ghost_node_update) F_data[ln]->beginGhostUpdate();
            if (X_data_ghost_node_update) X_data[ln]->beginGhostUpdate();
        }
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            if (F_data_ghost_node_update) F_data[ln]->endGhostUpdate();
            if (X_data_ghost_node_update) X_data[ln]->endGhostUpdate();
        }
    }

    d_hier_sc_data_ops->setToScalar(f_data_idx, 0.0);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > f_data = patch->getPatchData(f_data_idx);
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());
                const SAMRAI::hier::Box<NDIM>& box = idx_data->getGhostBox();
                IBTK::LEInteractor::spread(
                    f_data, F_data[ln], X_data[ln], idx_data,
                    patch, box, periodic_shift,
                    d_delta_fcn);
            }
        }
    }
    return;
}// spread

void
IBStaggeredHierarchyIntegrator::interp(
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > U_data,
    const int u_data_idx,
    const bool u_data_ghost_cell_update,
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data,
    const bool X_data_ghost_node_update,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? d_hierarchy->getFinestLevelNumber() : finest_ln_in);
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    if (u_data_ghost_cell_update)
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var;
        var_db->mapIndexToVariable(u_data_idx, var);
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_op = grid_geom->lookupRefineOperator(
            var, "CONSERVATIVE_LINEAR_REFINE");
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > refine_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
        refine_alg->registerRefine(u_data_idx, // destination
                                   u_data_idx, // source
                                   u_data_idx, // temporary work space
                                   refine_op);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > refine_sched = d_rscheds["u_interp->u_interp::S->S::CONSERVATIVE_LINEAR_REFINE"][ln];
            refine_alg->resetSchedule(refine_sched);
            refine_sched->fillData(d_integrator_time);
            d_ralgs["u_interp->u_interp::S->S::CONSERVATIVE_LINEAR_REFINE"]->resetSchedule(refine_sched);
        }
    }

    if (X_data_ghost_node_update)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                X_data[ln]->beginGhostUpdate();
            }
        }
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                X_data[ln]->endGhostUpdate();
            }
        }
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > u_data = patch->getPatchData(u_data_idx);
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());
                const SAMRAI::hier::Box<NDIM>& box = idx_data->getBox();
                IBTK::LEInteractor::interpolate(
                    U_data[ln], X_data[ln], idx_data, u_data,
                    patch, box, periodic_shift,
                    d_delta_fcn);
            }
        }
    }
    return;
}// interp

double
IBStaggeredHierarchyIntegrator::getLevelDt(
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx) const
{
    double stable_dt = std::numeric_limits<double>::max();
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        stable_dt = std::min(stable_dt,getPatchDt(patch,ctx));
    }
    stable_dt = SAMRAI::tbox::SAMRAI_MPI::minReduction(stable_dt);
    return stable_dt;
}// getLevelDt

double
IBStaggeredHierarchyIntegrator::getPatchDt(
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx) const
{
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
        patch->getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const SAMRAI::hier::Index<NDIM>& ilower = patch->getBox().lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch->getBox().upper();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > u_data =
        patch->getPatchData(d_u_var, ctx);
    const SAMRAI::hier::IntVector<NDIM>& u_ghost_cells = u_data->getGhostCellWidth();

    double stable_dt = std::numeric_limits<double>::max();
#if (NDIM == 2)
    NAVIER_STOKES_SC_STABLEDT_F77(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),
        u_ghost_cells(0),u_ghost_cells(1),
        u_data->getPointer(0),u_data->getPointer(1),
        stable_dt);
#endif
#if (NDIM == 3)
    NAVIER_STOKES_SC_STABLEDT_F77(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
        u_ghost_cells(0),u_ghost_cells(1),u_ghost_cells(2),
        u_data->getPointer(0),u_data->getPointer(1),u_data->getPointer(2),
        stable_dt);
#endif
    return stable_dt;
}// getPatchDt

void
IBStaggeredHierarchyIntegrator::resetLagrangianForceStrategy(
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
IBStaggeredHierarchyIntegrator::resetAnchorPointValues(
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > V_data,
    const int coarsest_ln,
    const int finest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            const int depth = V_data[ln]->getDepth();
            Vec V_vec = V_data[ln]->getGlobalVec();
            double* V_arr;
            int ierr = VecGetArray(V_vec, &V_arr);  IBTK_CHKERRQ(ierr);
            for (std::set<int>::const_iterator cit = d_anchor_point_local_idxs[ln].begin();
                 cit != d_anchor_point_local_idxs[ln].end(); ++cit)
            {
                const int& i = *cit;
                for (int d = 0; d < depth; ++d)
                {
                    V_arr[depth*i+d] = 0.0;
                }
            }
            ierr = VecRestoreArray(V_vec, &V_arr);  IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// resetAnchorPointValues

void
IBStaggeredHierarchyIntegrator::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
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

    d_conservation_form = db->getBoolWithDefault(
        "conservation_form", d_conservation_form);

    d_using_vorticity_tagging = db->getBoolWithDefault(
        "using_vorticity_tagging", d_using_vorticity_tagging);

    if (d_using_vorticity_tagging)
    {
        if (db->keyExists("vorticity_rel_thresh"))
        {
            d_omega_rel_thresh = db->getDoubleArray("vorticity_rel_thresh");
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n"
                         << "  Vorticity tagging is enabled but key data `vorticity_rel_thresh' not found in input.\n"
                         << "  Using default values for all levels in the locally refined grid.\n");
            d_omega_rel_thresh.resizeArray(1);
            d_omega_rel_thresh[0] = 1.0;
        }

        for (int i = 0; i < d_omega_rel_thresh.getSize(); ++i)
        {
            if (d_omega_rel_thresh[i] < 0.0 || d_omega_rel_thresh[i] > 1.0)
            {
                TBOX_ERROR(d_object_name << ":  "
                           << "relative vorticity thresholds for each level must lie in the interval [0,1].\n");
            }
        }

        if (db->keyExists("vorticity_abs_thresh"))
        {
            d_omega_abs_thresh = db->getDoubleArray("vorticity_abs_thresh");
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n"
                         << "  Vorticity tagging is enabled but key data `vorticity_abs_thresh' not found in input.\n"
                         << "  Using default values for all levels in the locally refined grid.\n");
            d_omega_abs_thresh.resizeArray(1);
            d_omega_abs_thresh[0] = std::numeric_limits<double>::max();
        }

        for (int i = 0; i < d_omega_abs_thresh.getSize(); ++i)
        {
            if (d_omega_abs_thresh[i] < 0.0)
            {
                TBOX_ERROR(d_object_name << ":  "
                           << "absolute vorticity thresholds for each level must be nonnegative.\n");
            }
        }
    }

    d_output_u      = db->getBoolWithDefault("output_u"     , d_output_u    );
    d_output_p      = db->getBoolWithDefault("output_p"     , d_output_p    );
    d_output_f      = db->getBoolWithDefault("output_f"     , d_output_f    );
    d_output_omega  = db->getBoolWithDefault("output_omega" , d_output_omega);
    d_output_div_u  = db->getBoolWithDefault("output_div_u" , d_output_div_u);

    d_u_scale = db->getDoubleWithDefault("u_scale", d_u_scale);
    d_p_scale = db->getDoubleWithDefault("p_scale", d_p_scale);
    d_f_scale = db->getDoubleWithDefault("f_scale", d_f_scale);

    d_use_exact_jacobian = db->getBoolWithDefault("use_exact_jacobian", d_use_exact_jacobian);

    d_cfl = db->getDoubleWithDefault("cfl",d_cfl);

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault(
        "dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault(
        "dt_max_time_min", d_dt_max_time_min);

    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    if (!is_from_restart)
    {
        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);

        d_normalize_pressure = db->getBoolWithDefault(
            "normalize_pressure", d_normalize_pressure);

        if (db->keyExists("rho"))
        {
            d_rho = db->getDouble("rho");
        }
        else
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `rho' not found in input.");
        }

        if (db->keyExists("mu"))
        {
            d_mu = db->getDouble("mu");
        }
        else
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `mu' not found in input.");
        }

        if (db->keyExists("lambda"))
        {
            d_lambda = db->getDouble("lambda");
        }
        else
        {
            d_lambda = 0.0;
        }

        d_delta_fcn = db->getStringWithDefault("delta_fcn", d_delta_fcn);
    }
    return;
}// getFromInput

void
IBStaggeredHierarchyIntegrator::getFromRestart()
{
    TBOX_ASSERT(false);
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBStaggeredHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
