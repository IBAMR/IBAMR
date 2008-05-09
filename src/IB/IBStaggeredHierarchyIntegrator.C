// Filename: IBStaggeredHierarchyIntegrator.C
// Last modified: <08.May.2008 21:46:57 griffith@box230.cims.nyu.edu>
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
#include <ibamr/INSStaggeredConvectiveOperator.h>
#include <ibamr/INSStaggeredNavierStokesOperator.h>
#include <ibamr/INSStaggeredProjectionPreconditioner.h>
#include <ibamr/INSStaggeredStokesOperator.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LNodeIndexData2.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScNewtonKrylovSolver.h>
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/SCLaplaceOperator.h>

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

// NOTE: The number of ghost cells required by the Godunov advection scheme
// depends on the order of the reconstruction.  These values were chosen to work
// with xsPPM7 (the modified piecewise parabolic method of Rider, Greenough, and
// Kamm).
static const int GADVECTG = 4;

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

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }

    // Set some default values.
    d_U_scale = 1.0;
    d_P_scale = 1.0;
    d_F_scale = 1.0;

    d_start_time = 0.0;
    d_end_time = std::numeric_limits<double>::max();
    d_grow_dt = 2.0;
    d_max_integrator_steps = std::numeric_limits<int>::max();

    d_using_vorticity_tagging = false;
    d_Omega_max = 0.0;

    d_second_order_pressure_update = false;
    d_normalize_pressure = false;

    d_num_init_cycles = 5;

    d_regrid_interval = 1;
    d_old_dt = -1.0;
    d_integrator_time = std::numeric_limits<double>::quiet_NaN();
    d_integrator_step = std::numeric_limits<int>::max();

    d_conservation_form = false;

    d_output_U = false;
    d_output_P = false;
    d_output_F = false;
    d_output_Omega = false;
    d_output_Div_U = false;

    d_rho    = std::numeric_limits<double>::quiet_NaN();
    d_mu     = std::numeric_limits<double>::quiet_NaN();
    d_nu     = std::numeric_limits<double>::quiet_NaN();
    d_lambda = std::numeric_limits<double>::quiet_NaN();

    d_cfl = 0.9;

    d_dt_max = std::numeric_limits<double>::max();
    d_dt_max_time_max = std::numeric_limits<double>::max();
    d_dt_max_time_min = -(d_dt_max_time_max-std::numeric_limits<double>::epsilon());

    d_is_initialized = false;

    d_do_log = false;

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Determine the ghost cell width required for cell-centered spreading and
    // interpolating.
    const int stencil_size = IBTK::LEInteractor::getStencilSize(d_delta_fcn);
    d_ghosts = int(floor(0.5*double(stencil_size)))+1;

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
    return;
}// ~IBStaggeredHierarchyIntegrator

const std::string&
IBStaggeredHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
IBStaggeredHierarchyIntegrator::registerVelocityInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> U_init)
{
    d_U_init = U_init;
    return;
}// registerVelocityInitialConditions

void
IBStaggeredHierarchyIntegrator::registerPressureInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> P_init)
{
    d_P_init = P_init;
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
    d_U_var          = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::U"          );
    d_U_cc_var       = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::U_cc",  NDIM);
    d_P_var          = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::P"          );
    d_F_var          = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::F"          );
    d_F_cc_var       = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::F_cc",  NDIM);
#if ( NDIM == 2)
    d_Omega_var      = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Omega"      );
#endif
#if ( NDIM == 3)
    d_Omega_var      = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Omega", NDIM);
    d_Omega_Norm_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::||Omega||_2");
#endif
    d_Div_U_var      = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Div_U"      );
    d_gadvect_U_var  = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::gadvect_U"  );

    // Create the default communication algorithms.
    d_fill_after_regrid = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    d_calgs["SYNCH_CURRENT_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    d_calgs["SYNCH_NEW_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();

    // Register state variables that are maintained by the
    // INSStaggeredHierarchyIntegrator.

    const SAMRAI::hier::IntVector<NDIM> cell_ghosts = CELLG;
    const SAMRAI::hier::IntVector<NDIM> side_ghosts = SIDEG;
    const SAMRAI::hier::IntVector<NDIM> gadvect_velocity_ghosts = std::max(GADVECTG,SIDEG);
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    registerVariable(d_U_current_idx, d_U_new_idx, d_U_scratch_idx,
                     d_U_var, side_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_U_cc_current_idx, d_U_cc_new_idx, d_U_cc_scratch_idx,
                     d_U_cc_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_P_current_idx, d_P_new_idx, d_P_scratch_idx,
                     d_P_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");

    registerVariable(d_F_current_idx, d_F_new_idx, d_F_scratch_idx,
                     d_F_var, side_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_F_cc_current_idx, d_F_cc_new_idx, d_F_cc_scratch_idx,
                     d_F_cc_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_Omega_current_idx, d_Omega_new_idx, d_Omega_scratch_idx,
                     d_Omega_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");
#if (NDIM == 3)
    registerVariable(d_Omega_Norm_current_idx, d_Omega_Norm_new_idx, d_Omega_Norm_scratch_idx,
                     d_Omega_Norm_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");
#endif
    registerVariable(d_Div_U_current_idx, d_Div_U_new_idx, d_Div_U_scratch_idx,
                     d_Div_U_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    // Register other scratch variables.
    d_gadvect_U_scratch_idx = var_db->registerVariableAndContext(
        d_gadvect_U_var, getScratchContext(), gadvect_velocity_ghosts);

    // Register variables for plotting.
    if (!d_visit_writer.isNull())
    {
        if (d_output_U)
        {
            d_visit_writer->registerPlotQuantity(
                d_U_var->getName(), "VECTOR", d_U_cc_current_idx, 0, d_U_scale);

            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    d_U_var->getName()+stream.str(), "SCALAR", d_U_cc_current_idx, d, d_U_scale);
            }
        }

        if (d_output_P)
        {
            d_visit_writer->registerPlotQuantity(
                d_P_var->getName(), "SCALAR", d_P_current_idx, 0, d_P_scale);
        }

        if (d_output_F)
        {
            d_visit_writer->registerPlotQuantity(
                d_F_var->getName(), "VECTOR", d_F_cc_current_idx, 0, d_F_scale);

            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    d_F_var->getName()+stream.str(), "SCALAR", d_F_cc_current_idx, d, d_F_scale);
            }
        }

        if (d_output_Omega)
        {
            d_visit_writer->registerPlotQuantity(
                d_Omega_var->getName(), (NDIM == 2) ? "SCALAR" : "VECTOR", d_Omega_current_idx);
#if (NDIM == 3)
            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    d_Omega_var->getName()+stream.str(), "SCALAR", d_Omega_current_idx, d);
            }
#endif
        }

        if (d_output_Div_U)
        {
            d_visit_writer->registerPlotQuantity(
                d_Div_U_var->getName(), "SCALAR", d_Div_U_current_idx);
        }
    }

    // Create several refinement communications algorithms, used in filling
    // ghost cell data.
    std::string ralg_name;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_operator;

    ralg_name = "gadvect_U_scratch_bdry_fill";
    d_ralgs[ralg_name] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(d_gadvect_U_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs[ralg_name]->registerRefine(d_gadvect_U_scratch_idx, // destination
                                       d_U_scratch_idx,         // source
                                       d_gadvect_U_scratch_idx, // temporary work space
                                       refine_operator);
    d_rstrategies[ralg_name] = new IBTK::CartExtrapPhysBdryOp(d_gadvect_U_scratch_idx, BDRY_EXTRAP_TYPE);

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
#if 0
    t_advance_hierarchy->start();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
#endif

    const double current_time = d_integrator_time;
    const double new_time     = d_integrator_time+dt;
    const bool initial_time   = SAMRAI::tbox::MathUtilities<double>::equalEps(current_time,d_start_time);

    // Regrid the patch hierarchy.
    const bool do_regrid = (d_regrid_interval == 0
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

    // Set the current time interval in the force specification object.
    d_force_strategy->setTimeInterval(current_time, new_time);

    // (Re)initialize the force strategy.
    if (d_force_strategy_needs_init)
    {
        resetLagrangianForceStrategy(current_time, initial_time);
        d_force_strategy_needs_init = false;
    }

    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > U_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_new_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_new_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > U_new_data(finest_ln+1);

    // Synchronize the Cartesian grid velocity u(n) on the patch hierarchy.
    //
    // NOTE: Since we are maintaining the Lagrangian velocity data, this step is
    // skipped for each timestep following the initial one, except for those
    // timesteps that immediately follow a regridding.
    if (initial_time || d_reinterpolate_after_regrid)
    {
        for (int ln = finest_ln; ln > coarsest_ln; --ln)
        {
            d_cscheds["U->U::C->C::CONSERVATIVE_COARSEN"][ln]->coarsenData();
        }
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(d_V_idx, current_time);
            d_rscheds["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(current_time);
        }
    }

    // Initialize the various LNodeLevelData objects, including X_data, F_data,
    // U_data, and X_new_data, on each level of the patch hierarchy.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,ln);
            U_data[ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::VEL_DATA_NAME,ln);
            F_data[ln] = d_lag_data_manager->createLNodeLevelData("F",ln,NDIM);

            X_data[ln]->restoreLocalFormVec();
            U_data[ln]->restoreLocalFormVec();
            F_data[ln]->restoreLocalFormVec();

            X_new_data[ln] = d_lag_data_manager->createLNodeLevelData("X_new",ln,NDIM);
            U_new_data[ln] = d_lag_data_manager->createLNodeLevelData("U_new",ln,NDIM);
            F_new_data[ln] = d_lag_data_manager->createLNodeLevelData("F_new",ln,NDIM);

            X_new_data[ln]->restoreLocalFormVec();
            U_new_data[ln]->restoreLocalFormVec();
            F_new_data[ln]->restoreLocalFormVec();
        }
    }

    // 1. Interpolate u(n) from the Cartesian grid onto the Lagrangian mesh.
    //
    // NOTE: Since we are maintaining the Lagrangian velocity data, this step is
    // skipped for each timestep following the initial one execpt immediately
    // following a regridding operation.
    if (initial_time || d_reinterpolate_after_regrid)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
            if (d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): interpolating u(n) to U(n) on level number " << ln << "\n";
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > v_data =
                        patch->getPatchData(d_V_idx);
                    const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(
                        d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());
                    IBTK::LEInteractor::interpolate(
                        U_data[ln], X_data[ln], idx_data, v_data,
                        patch, patch_box, periodic_shift,
                        d_delta_fcn);
                }
            }
        }
    }

    // 2. Compute F(n) = F(X(n),n), the Lagrangian force corresponding to
    //    configuration X(n) at time t_{n}.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing F(n) on level number " << ln << "\n";
            Vec F_vec = F_data[ln]->getGlobalVec();
            int ierr = VecSet(F_vec, 0.0);  IBTK_CHKERRQ(ierr);
            d_force_strategy->computeLagrangianForce(
                F_data[ln], X_data[ln],
                d_hierarchy, ln, current_time, d_lag_data_manager);
        }
    }

    // 3. Spread F(n) from the Lagrangian mesh onto the Cartesian grid.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
        level->allocatePatchData(d_F_scratch1_idx, current_time);

        // On the coarsest level of the patch hierarchy, simply initialize the
        // Cartesian force density to equal zero.
        //
        // For each of the finer levels in the patch hierarchy, conservatively
        // interpolate the Cartesian force density from coarser levels in the
        // patch hierarchy.
        if (ln == coarsest_ln)
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_current_data = patch->
                    getPatchData(d_ins_hier_integrator->getForceVar(),
                                 d_ins_hier_integrator->getCurrentContext());
                f_current_data->fillAll(0.0);
            }
        }
        else
        {
            // Interpolate the Cartesian force from the next coarser level in
            // the hierarchy.
            //
            // Note that these refine schedules initialize both the force data
            // maintained by the IBStaggeredHierarchyIntegrator and the current context
            // of the force data maintained by the Navier-Stokes solver.
            d_force_current_rscheds[ln]->fillData(current_time);
        }

        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln]->beginGhostUpdate();
            F_data[ln]->beginGhostUpdate();
            X_data[ln]->endGhostUpdate();
            F_data[ln]->endGhostUpdate();
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): spreading F(n) to f(n) on level number " << ln << "\n";
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_current_data = patch->
                    getPatchData(d_ins_hier_integrator->getForceVar(),
                                 d_ins_hier_integrator->getCurrentContext());
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->
                    getPatchData(d_lag_data_manager->
                                 getLNodeIndexPatchDescriptorIndex());
                IBTK::LEInteractor::spread(
                    f_current_data, F_data[ln], X_data[ln], idx_data,
                    patch, SAMRAI::hier::Box<NDIM>::grow(patch_box,d_ghosts), periodic_shift,
                    d_delta_fcn);
            }
        }
    }

    // 4. Compute X~(n+1), the preliminary structure configuration at time
    //    t_{n+1}.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing X~(n+1) on level number " << ln << "\n";
            Vec U_vec = U_data[ln]->getGlobalVec();
            Vec X_vec = X_data[ln]->getGlobalVec();
            Vec X_new_vec = X_new_data[ln]->getGlobalVec();
            int ierr = VecWAXPY(X_new_vec, dt, U_vec, X_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // 5. Compute F~(n+1) = F(X~(n+1),n+1), the Lagrangian force corresponding
    //    to configuration X~(n+1) at time t_{n+1}.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing F~(n+1) on level number " << ln << "\n";
            Vec F_new_vec = F_new_data[ln]->getGlobalVec();
            int ierr = VecSet(F_new_vec, 0.0);  IBTK_CHKERRQ(ierr);
            d_force_strategy->computeLagrangianForce(
                F_new_data[ln], X_new_data[ln],
                d_hierarchy, ln, new_time, d_lag_data_manager);
        }
    }

    // 6. Spread F~(n+1) from the Lagrangian mesh onto the Cartesian grid.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
        level->allocatePatchData(d_F_idx, current_time);
        level->allocatePatchData(d_F_scratch2_idx, current_time);

        // On the coarsest level of the patch hierarchy, simply initialize the
        // Cartesian force density to equal zero.
        //
        // For each of the finer levels in the patch hierarchy, conservatively
        // interpolate the Cartesian force density from coarser levels in the
        // patch hierarchy.
        if (ln == coarsest_ln)
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_new_data = patch->
                    getPatchData(d_F_idx);
                f_new_data->fillAll(0.0);
            }
        }
        else
        {
            // Interpolate the Cartesian force from the next coarser level in
            // the hierarchy.
            //
            // Note that these refine schedules initialize both the force data
            // maintained by the IBStaggeredHierarchyIntegrator and the current context
            // of the force data maintained by the Navier-Stokes solver.
            d_force_new_rscheds[ln]->fillData(current_time);
        }

        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_new_data[ln]->beginGhostUpdate();
            F_new_data[ln]->beginGhostUpdate();
            X_new_data[ln]->endGhostUpdate();
            F_new_data[ln]->endGhostUpdate();
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): spreading F~(n+1) to f~(n+1) on level number " << ln << "\n";
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_new_data = patch->
                    getPatchData(d_F_idx);
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->
                    getPatchData(d_lag_data_manager->
                                 getLNodeIndexPatchDescriptorIndex());
                IBTK::LEInteractor::spread(
                    f_new_data, F_new_data[ln], X_new_data[ln], idx_data,
                    patch, SAMRAI::hier::Box<NDIM>::grow(patch_box,d_ghosts), periodic_shift,
                    d_delta_fcn);
            }
        }
    }

    // 7. If an additional body force specification object is provided, compute
    //    the body force at the beginning and end of the time step, and add
    //    those values to f(n) and f(n+1).
    if (!d_body_force_set.isNull())
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        const int F_current_idx = var_db->mapVariableAndContextToIndex(
            d_ins_hier_integrator->getForceVar(),
            d_ins_hier_integrator->getCurrentContext());

        d_body_force_set->setDataOnPatchHierarchy(
            d_F_scratch1_idx, d_F_var, d_hierarchy,
            current_time, false, coarsest_ln, finest_ln);

        d_hier_cc_data_ops->add(F_current_idx, F_current_idx, d_F_scratch1_idx);

        d_body_force_set->setDataOnPatchHierarchy(
            d_F_scratch2_idx, d_F_var, d_hierarchy,
            current_time+dt, false, coarsest_ln, finest_ln);

        d_hier_cc_data_ops->add(d_F_idx, d_F_idx, d_F_scratch2_idx);
    }

    // Synchronize the Cartesian grid force field on the patch hierarchy.
    //
    // Note that these coarsen schedules synchronize both the force data
    // maintained by the IBStaggeredHierarchyIntegrator and the current context of the
    // force data maintained by the Navier-Stokes solver.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["F->F::CONSERVATIVE_COARSEN"][ln]->coarsenData();
    }

    // Deallocate unneeded scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (initial_time || d_reinterpolate_after_regrid) level->deallocatePatchData(d_V_idx);
        level->deallocatePatchData(d_F_scratch1_idx);
        level->deallocatePatchData(d_F_scratch2_idx);
    }

    // Reset the reinterpolation flag.
    d_reinterpolate_after_regrid = false;

    // Solve the incompressible Navier-Stokes equations for U(n+1) and P(n+1/2).
    //
    // Note that the pressure at start_time is not an initial value for the
    // incompressible Navier-Stokes equations, so we must solve for it by
    // cycling the solution for the first timestep:
    //
    // Solve the Navier-Stokes equations with initial guess P(n=0)=0, solve for
    // P(n=1/2), discard U(n=1) and u(n=1), rinse, wash, repeat.
    //
    // For all other timesteps, we just use the previous value of P as the guess
    // for P(n+1/2).
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing u(n+1), p(n+1/2)\n";

    int cycle = 0;
    const bool performing_init_cycles = initial_time;
    const int num_cycles = (performing_init_cycles
                            ? d_num_init_cycles
                            : d_num_cycles);
    for (cycle = 0; cycle < num_cycles; ++cycle)
    {
        if (d_do_log && performing_init_cycles)
        {
            SAMRAI::tbox::plog << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++\n";
            SAMRAI::tbox::plog << "+\n";
            SAMRAI::tbox::plog << "+ Performing cycle " << cycle+1 << " of " << d_num_init_cycles << " to initialize P(n=1/2)\n";
            SAMRAI::tbox::plog << "+\n";
            SAMRAI::tbox::plog << "++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
        }

        // Solve the Navier-Stokes equations for U(n+1), u(n+1), P(n+1/2).  For
        // better or worse, each of the major algorithmic steps of the
        // projection method is separated into its own member function.
        d_ins_hier_integrator->predictAdvectionVelocity(current_time, new_time);

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

    // Synchronize the Cartesian grid velocity u(n+1) on the patch hierarchy.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["U->U::N->N::CONSERVATIVE_COARSEN"][ln]->coarsenData();
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_W_idx, current_time);
        d_rscheds["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(current_time);
    }

    // In the following loop, we perform the following operations:
    //
    // 1. Interpolate u(n+1) from the Cartesian grid onto the Lagrangian mesh
    //    using X~(n+1).
    //
    // 2. Compute X(n+1), the final structure configuration at time t_{n+1}.
    //
    // 3. Interpolate u(n+1) from the Cartesian grid onto the Lagrangian mesh
    //    using X(n+1).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(
            level->getRatio());

        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            int ierr;

            // 1. Interpolate u(n+1) from the Cartesian grid onto the Lagrangian
            //    mesh.
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): interpolating u(n+1) to U(n+1) on level number " << ln << "\n";

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > w_data =
                    patch->getPatchData(d_W_idx);
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(
                    d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());

                IBTK::LEInteractor::interpolate(
                    U_new_data[ln], X_new_data[ln], idx_data, w_data,
                    patch, patch_box, periodic_shift,
                    d_delta_fcn);
            }

            // 2. Compute X(n+1), the final structure configuration at time
            //    t_{n+1}.
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing X(n+1) on level number " << ln << "\n";

            Vec U_new_vec = U_new_data[ln]->getGlobalVec();
            Vec X_vec = X_data[ln]->getGlobalVec();
            Vec X_new_vec = X_new_data[ln]->getGlobalVec();
            ierr = VecAXPY(X_new_vec, dt, U_new_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(X_vec, 0.5, 0.5, X_new_vec);  IBTK_CHKERRQ(ierr);

            // 3. Interpolate u(n+1) from the Cartesian grid onto the Lagrangian
            //    mesh using X(n+1).
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): interpolating u(n+1) to U(n+1) on level number " << ln << "\n";

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > w_data =
                    patch->getPatchData(d_W_idx);
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(
                    d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());

                IBTK::LEInteractor::interpolate(
                    U_data[ln], X_data[ln], idx_data, w_data,
                    patch, patch_box, periodic_shift,
                    d_delta_fcn);
            }
        }
    }

    // Deallocate any remaining scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_W_idx);
        level->deallocatePatchData(d_F_idx);
    }

    // Reset all time dependent data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierData(new_time);

    // Determine the next stable timestep.
    double dt_next = d_ins_hier_integrator->getStableTimestep();

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
    t_advance_hierarchy->stop();
    return dt_next;
#endif
    t_advance_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
#endif

    const double current_time = d_integrator_time;
    const double new_time = d_integrator_time+dt;
    const bool initial_time = SAMRAI::tbox::MathUtilities<double>::equalEps(d_integrator_time,d_start_time);

    // Set the guess for the initial pressure to zero.
    if (initial_time)
    {
        d_hier_cc_data_ops->setToScalar(d_P_current_idx, 0.0);
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

    // Hierarchy ghost cell transaction components.
    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_scratch_component(d_U_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);
    InterpolationTransactionComponent P_scratch_component(d_P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);

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

    // Setup the solver vectors.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_scratch_vec, U_rhs_vec;
    U_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::U_scratch_vec", d_hierarchy, 0, finest_ln);
    U_scratch_vec->addComponent(d_U_var, d_U_scratch_idx, d_wgt_sc_idx, d_hier_sc_data_ops);
    U_rhs_vec = U_scratch_vec->cloneVector(d_object_name+"::U_rhs_vec");
    U_rhs_vec->allocateVectorData(current_time);
    const int U_rhs_idx = U_rhs_vec->getComponentDescriptorIndex(0);
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_rhs_var = U_rhs_vec->getComponentVariable(0);
    d_hier_sc_data_ops->setToScalar(U_rhs_idx,0.0);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > P_scratch_vec, P_rhs_vec;
    P_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::P_scratch_vec", d_hierarchy, 0, finest_ln);
    P_scratch_vec->addComponent(d_P_var, d_P_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);
    P_rhs_vec = P_scratch_vec->cloneVector(d_object_name+"::P_rhs_vec");
    P_rhs_vec->allocateVectorData(current_time);
    const int P_rhs_idx = P_rhs_vec->getComponentDescriptorIndex(0);
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_rhs_var = P_rhs_vec->getComponentVariable(0);
    d_hier_cc_data_ops->setToScalar(P_rhs_idx,0.0);

    const int num_cycles = initial_time ? d_num_init_cycles : d_num_cycles;
    for (int cycle = 0; cycle < num_cycles; ++cycle)
    {
        if (d_do_log && initial_time)
        {
            SAMRAI::tbox::plog << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++\n";
            SAMRAI::tbox::plog << "+\n";
            SAMRAI::tbox::plog << "+ Performing cycle " << cycle+1 << " of " << d_num_init_cycles << " to initialize the pressure at time t=dt/2\n";
            SAMRAI::tbox::plog << "+\n";
            SAMRAI::tbox::plog << "++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
        }

        // Initialize the right-hand side terms.
        SAMRAI::solv::PoissonSpecifications rhs_spec(d_object_name+"::rhs_spec");
        rhs_spec.setCConstant((d_rho/dt)-0.5*d_lambda);
        rhs_spec.setDConstant(          +0.5*d_mu    );
        d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
        d_U_scratch_bdry_fill_op->resetTransactionComponent(U_scratch_component);
        d_hier_math_ops->laplace(
            U_rhs_idx, U_rhs_var,
            rhs_spec,
            d_U_scratch_idx, d_U_var,
            d_U_scratch_bdry_fill_op, current_time);
        if (!d_F_set.isNull())
        {
            if (d_F_set->isTimeDependent())
            {
                d_F_set->setDataOnPatchHierarchy(
                    d_F_new_idx, d_F_var, d_hierarchy, new_time);
            }
            d_hier_sc_data_ops->axpy(U_rhs_idx, +0.5, d_F_current_idx, U_rhs_idx);
            d_hier_sc_data_ops->axpy(U_rhs_idx, +0.5, d_F_new_idx    , U_rhs_idx);
        }
        else
        {
            d_hier_sc_data_ops->setToScalar(d_F_new_idx, 0.0);
        }

        // Setup the Helmholtz solver.
        SAMRAI::solv::PoissonSpecifications helmholtz_spec(d_object_name+"::helmholtz_spec");
        helmholtz_spec.setCConstant((d_rho/dt)+0.5*d_lambda);
        helmholtz_spec.setDConstant(          -0.5*d_mu    );
        SAMRAI::tbox::Pointer<IBTK::SCLaplaceOperator> helmholtz_operator = new IBTK::SCLaplaceOperator(
            d_object_name+"::Helmholtz Operator", helmholtz_spec, NULL);  // XXXX

        helmholtz_operator->setPoissonSpecifications(helmholtz_spec);
        //helmholtz_operator->setPhysicalBcCoefs(U_bc_coefs);  // XXXX
        helmholtz_operator->setHomogeneousBc(false);  // XXXX
        helmholtz_operator->setTime(new_time);
        helmholtz_operator->setHierarchyMathOps(d_hier_math_ops);

        SAMRAI::tbox::Pointer<IBTK::KrylovLinearSolver> helmholtz_solver = new IBTK::PETScKrylovLinearSolver(
            d_object_name+"::PETSc Krylov solver", "adv_diff_");
        helmholtz_solver->setInitialGuessNonzero(false);
        helmholtz_solver->setOperator(helmholtz_operator);

        // Reset the solution, rhs, and nullspace vectors.
        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > sol_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
            d_object_name+"::sol_vec", d_hierarchy, 0, finest_ln);
        sol_vec->addComponent(d_U_var,d_U_scratch_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
        sol_vec->addComponent(d_P_var,d_P_scratch_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > rhs_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
            d_object_name+"::rhs_vec", d_hierarchy, 0, finest_ln);
        rhs_vec->addComponent(d_U_var,U_rhs_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
        rhs_vec->addComponent(d_P_var,P_rhs_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > nul_vec = sol_vec->cloneVector(d_object_name+"::nul_vec");
        nul_vec->allocateVectorData(current_time);
        d_hier_sc_data_ops->setToScalar(nul_vec->getComponentDescriptorIndex(0), 0.0);
        d_hier_cc_data_ops->setToScalar(nul_vec->getComponentDescriptorIndex(1), 1.0);

        // Setup the nullspace object.
        PetscErrorCode ierr;
        MatNullSpace petsc_nullsp;
        Vec petsc_nullsp_vec = IBTK::PETScSAMRAIVectorReal<double>::createPETScVector(nul_vec, PETSC_COMM_SELF);
        Vec vecs[] = {petsc_nullsp_vec};
        static const PetscTruth has_cnst = PETSC_FALSE;
        ierr = MatNullSpaceCreate(PETSC_COMM_SELF, has_cnst, 1, vecs, &petsc_nullsp); IBTK_CHKERRQ(ierr);

        // Setup the Stokes operator.
        SAMRAI::tbox::Pointer<INSStaggeredStokesOperator> stokes_op = new INSStaggeredStokesOperator(
            d_rho, d_mu, d_lambda,
            d_hier_math_ops,
            d_U_P_scratch_bdry_fill_op);
        stokes_op->setCurrentTimeInterval(current_time,new_time);

        // Setup the convective operator.
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator = grid_geom->lookupRefineOperator(d_gadvect_U_var, "CONSERVATIVE_LINEAR_REFINE");
        SAMRAI::tbox::Pointer<INSStaggeredConvectiveOperator> convective_op = new INSStaggeredConvectiveOperator(
            d_rho, d_mu, d_lambda,
            d_conservation_form,
            d_ralgs["gadvect_U_scratch_bdry_fill"], refine_operator, d_rscheds["gadvect_U_scratch_bdry_fill"]);
        convective_op->setVelocityScratchPatchDataIndex(d_gadvect_U_scratch_idx);

        // Setup the Navier-Stokes operator.
        SAMRAI::tbox::Pointer<INSStaggeredNavierStokesOperator> navier_stokes_op = new INSStaggeredNavierStokesOperator(
            d_rho, d_mu, d_lambda,
            stokes_op, convective_op,
            d_hier_sc_data_ops);
        navier_stokes_op->setVelocityCurrentPatchDataIndex(d_U_current_idx);

        // Setup the projection preconditioner.
        static const std::string projection_type = "pressure_increment";
        SAMRAI::tbox::Pointer<INSStaggeredProjectionPreconditioner> projection_pc = new INSStaggeredProjectionPreconditioner(
            d_rho, d_mu, d_lambda,
            projection_type,
            d_normalize_pressure,
            helmholtz_solver, d_hier_projector,
            d_hier_cc_data_ops, d_hier_sc_data_ops, d_hier_math_ops,
            d_P_scratch_bdry_fill_op);
        projection_pc->setCurrentTimeInterval(current_time,new_time);
        projection_pc->initializeSolverState(*sol_vec,*rhs_vec);

        // Setup the nonlinear solver.
        SAMRAI::tbox::Pointer<IBTK::PETScNewtonKrylovSolver> nonlinear_solver = new IBTK::PETScNewtonKrylovSolver(
            d_object_name+"::nonlinear_solver", "ins_");
        nonlinear_solver->setOperator(navier_stokes_op);
        nonlinear_solver->initializeSolverState(*sol_vec,*rhs_vec);
        if (d_normalize_pressure)
        {
            SNES petsc_snes = nonlinear_solver->getPETScSNES();
            KSP petsc_ksp;
            ierr = SNESGetKSP(petsc_snes, &petsc_ksp); IBTK_CHKERRQ(ierr);
            ierr = KSPSetNullSpace(petsc_ksp, petsc_nullsp); IBTK_CHKERRQ(ierr);
        }
        nonlinear_solver->getLinearSolver()->setPreconditioner(projection_pc);

        // Solve system.
        d_hier_sc_data_ops->copyData(sol_vec->getComponentDescriptorIndex(0), (cycle == 0 ? d_U_current_idx : d_U_new_idx));
        d_hier_cc_data_ops->copyData(sol_vec->getComponentDescriptorIndex(1), (cycle == 0 ? d_P_current_idx : d_P_new_idx));
        nonlinear_solver->solveSystem(*sol_vec,*rhs_vec);

        // Pull out solution components.
        d_hier_sc_data_ops->copyData(d_U_new_idx, sol_vec->getComponentDescriptorIndex(0));
        d_hier_cc_data_ops->copyData(d_P_new_idx, sol_vec->getComponentDescriptorIndex(1));

        // Deallocate the nullspace object.
        if (petsc_nullsp != static_cast<MatNullSpace>(NULL))
        {
            ierr = MatNullSpaceDestroy(petsc_nullsp); IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(petsc_nullsp_vec); IBTK_CHKERRQ(ierr);
        }

        // Deallocate scratch data.
        nul_vec->deallocateVectorData();

        // Reset the value the current estimate of the pressure when performing
        // multiple cycles.
        if (cycle < num_cycles-1)
        {
            d_hier_cc_data_ops->copyData(d_P_current_idx, d_P_new_idx);
        }

        // Reset the transaction components.
        typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        InterpolationTransactionComponent U_scratch_component(d_U_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX
        InterpolationTransactionComponent P_scratch_component(d_P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX

        std::vector<InterpolationTransactionComponent> U_P_scratch_components(2);
        U_P_scratch_components[0] = U_scratch_component;
        U_P_scratch_components[1] = P_scratch_component;

        d_U_P_scratch_bdry_fill_op->resetTransactionComponents(U_P_scratch_components);
        d_U_scratch_bdry_fill_op->resetTransactionComponent(U_scratch_component);
        d_P_scratch_bdry_fill_op->resetTransactionComponent(P_scratch_component);
    }

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

    // Setup U_scratch to allow for ghost cell filling.
    d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_new_idx);
    d_U_scratch_bdry_fill_op->resetTransactionComponent(U_scratch_component);
    d_U_scratch_bdry_fill_op->fillData(new_time);

    // Compute Omega = curl U.
    d_hier_math_ops->curl(
        d_Omega_new_idx, d_Omega_var,
        d_U_scratch_idx, d_U_var,
        d_no_fill_op, new_time);
#if (NDIM == 3)
    d_hier_math_ops->pointwise_L2Norm(
        d_Omega_Norm_new_idx, d_Omega_Norm_var,
        d_Omega_new_idx, d_Omega_var);
#endif

    // Compute max ||Omega||_2.
#if (NDIM == 2)
    d_Omega_max = std::max(+d_hier_cc_data_ops->max(d_Omega_new_idx),
                           -d_hier_cc_data_ops->min(d_Omega_new_idx));
#endif
#if (NDIM == 3)
    d_Omega_max = d_hier_cc_data_ops->max(d_Omega_Norm_new_idx);
#endif
    // Compute Div U.
    d_hier_math_ops->div(
        d_Div_U_new_idx, d_Div_U_var,
        1.0, d_U_scratch_idx, d_U_var,
        d_no_fill_op, new_time, false);

    // Deallocate scratch data.
    U_rhs_vec->deallocateVectorData();
    P_rhs_vec->deallocateVectorData();

    // Synchronize all data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): synchronizing data\n";
    synchronizeHierarchy();

    // Reset all time dependent data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierData(new_time);

    // Determine the next stable timestep from u(n+1).
    const double dt_next = getStableTimestep(getCurrentContext());

    t_advance_hierarchy->stop();
    return dt_next;
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

    // Update the workload post-regridding.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_lag_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());

    // Indicate that the force strategy needs to be re-initialized.
    d_force_strategy_needs_init  = true;

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
        if (d_U_init.isNull())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_current_data =
                    patch->getPatchData(d_U_current_idx);
                U_current_data->fillAll(0.0);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_cc_current_data =
                    patch->getPatchData(d_U_cc_current_idx);
                U_cc_current_data->fillAll(0.0);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_current_data =
                    patch->getPatchData(d_Omega_current_idx);
                Omega_current_data->fillAll(0.0);
#if (NDIM == 3)
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_Norm_current_data =
                    patch->getPatchData(d_Omega_Norm_current_idx);
                Omega_Norm_current_data->fillAll(0.0);
#endif
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Div_U_current_data =
                    patch->getPatchData(d_Div_U_current_idx);
                Div_U_current_data->fillAll(0.0);
            }
        }
        else
        {
            level->allocatePatchData(d_U_scratch_idx, init_data_time);

            // Initialize U.
            d_U_init->setDataOnPatchLevel(
                d_U_current_idx, d_U_var, level,
                init_data_time, initial_time);
            IBTK::PatchMathOps patch_math_ops;
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_current_data =
                    patch->getPatchData(d_U_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_cc_current_data =
                    patch->getPatchData(d_U_cc_current_idx);

                patch_math_ops.interp(U_cc_current_data, U_current_data, patch);
            }

            // Fill in U boundary data from coarser levels.
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom =
                d_hierarchy->getGridGeometry();
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > ralg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator = grid_geom->lookupRefineOperator(
                d_U_var, "CONSERVATIVE_LINEAR_REFINE");
            ralg->registerRefine(d_U_scratch_idx, // destination
                                 d_U_current_idx, // source
                                 d_U_scratch_idx, // temporary work space
                                 refine_operator);
            IBTK::CartExtrapPhysBdryOp bc_op(d_U_scratch_idx, BDRY_EXTRAP_TYPE);
            ralg->createSchedule(level, level_number-1, hierarchy, &bc_op)->fillData(init_data_time);

            // Initialize quantities derived from the initial value of U.
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_scratch_data =
                    patch->getPatchData(d_U_scratch_idx);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_current_data =
                    patch->getPatchData(d_Omega_current_idx);
                patch_math_ops.curl(Omega_current_data, U_scratch_data, patch);
#if (NDIM == 3)
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_Norm_current_data =
                    patch->getPatchData(d_Omega_Norm_current_idx);
                patch_math_ops.pointwise_L2Norm(Omega_Norm_current_data, Omega_current_data, patch);
#endif
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Div_U_current_data =
                    patch->getPatchData(d_Div_U_current_idx);
                patch_math_ops.div(Div_U_current_data, 1.0, U_scratch_data, 0.0, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> >(NULL), patch);
            }

            level->deallocatePatchData(d_U_scratch_idx);
        }

        // Initialize the maximum value of ||Omega||_2 on the grid.
        if (level_number == 0)
        {
            d_Omega_max = 0.0;
        }

        SAMRAI::math::PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
#if (NDIM == 2)
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_current_data =
                patch->getPatchData(d_Omega_current_idx);
            d_Omega_max = std::max(d_Omega_max, +patch_cc_data_ops.max(Omega_current_data, patch_box));
            d_Omega_max = std::max(d_Omega_max, -patch_cc_data_ops.min(Omega_current_data, patch_box));
#endif
#if (NDIM == 3)
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_Norm_current_data =
                patch->getPatchData(d_Omega_Norm_current_idx);
            d_Omega_max = std::max(d_Omega_max, patch_cc_data_ops.max(Omega_Norm_current_data, patch_box));
#endif
        }

        // If no initialization object is provided, initialize the pressure to
        // zero.  Otherwise, use the initialization object to set the pressure
        // to some specified value.
        //
        // NOTE: This initial value for the pressure IS NOT USED by the time
        // integrator and is only specified for purposes of visualization.
        if (d_P_init.isNull())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > P_current_data =
                    patch->getPatchData(d_P_current_idx);
                P_current_data->fillAll(0.0);
            }
        }
        else
        {
            // Initialize P.
            d_P_init->setDataOnPatchLevel(
                d_P_current_idx, d_P_var, level,
                init_data_time, initial_time);
        }

        // If no initialization object is provided, initialize the body force to
        // zero.  Otherwise, use the initialization object to set the body force
        // to some specified value.
        if (d_F_set.isNull())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > F_current_data =
                    patch->getPatchData(d_F_current_idx);
                F_current_data->fillAll(0.0);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_cc_current_data =
                    patch->getPatchData(d_F_cc_current_idx);
                F_cc_current_data->fillAll(0.0);
            }
        }
        else
        {
            // Initialize F.
            d_F_set->setDataOnPatchLevel(
                d_F_current_idx, d_F_var, level,
                init_data_time, initial_time);
            IBTK::PatchMathOps patch_math_ops;
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > F_current_data =
                    patch->getPatchData(d_F_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_cc_current_data =
                    patch->getPatchData(d_F_cc_current_idx);

                patch_math_ops.interp(F_cc_current_data, F_current_data, patch);
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
    InterpolationTransactionComponent U_scratch_component(d_U_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX
    InterpolationTransactionComponent P_scratch_component(d_P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX

    std::vector<InterpolationTransactionComponent> U_P_scratch_components(2);
    U_P_scratch_components[0] = U_scratch_component;
    U_P_scratch_components[1] = P_scratch_component;

    d_U_P_scratch_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_U_P_scratch_bdry_fill_op->initializeOperatorState(U_P_scratch_components, d_hierarchy);

    d_U_scratch_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_U_scratch_bdry_fill_op->initializeOperatorState(U_scratch_component, d_hierarchy);

    d_P_scratch_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_P_scratch_bdry_fill_op->initializeOperatorState(P_scratch_component, d_hierarchy);

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
        const double Omega_rel_thresh =
            (level_number >= 0 && level_number < d_Omega_rel_thresh.getSize()
             ? d_Omega_rel_thresh[level_number]
             : (level_number < 0
                ? d_Omega_rel_thresh[0]
                : d_Omega_rel_thresh[d_Omega_rel_thresh.size()-1]));
        const double Omega_abs_thresh =
            (level_number >= 0 && level_number < d_Omega_abs_thresh.getSize()
             ? d_Omega_abs_thresh[level_number]
             : (level_number < 0
                ? d_Omega_abs_thresh[0]
                : d_Omega_abs_thresh[d_Omega_abs_thresh.size()-1]));
        if (Omega_rel_thresh > 0.0 && Omega_abs_thresh > 0.0)
        {
            const double thresh = sqrt(std::numeric_limits<double>::epsilon()) +
                std::min(Omega_rel_thresh*d_Omega_max, Omega_abs_thresh);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags_data = patch->
                    getPatchData(tag_index);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_current_data = patch->
                    getPatchData(d_Omega_current_idx);

                for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = ic();
#if (NDIM == 2)
                    if (std::abs((*Omega_current_data)(i)) > thresh)
                    {
                        (*tags_data)(i) = 1;
                    }
#endif
#if (NDIM == 3)
                    double norm_Omega_sq = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        norm_Omega_sq += (*Omega_current_data)(i,d)*(*Omega_current_data)(i,d);
                    }
                    const double norm_Omega = sqrt(norm_Omega_sq);
                    if (norm_Omega > thresh)
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
    return d_U_var;
}// getVelocityVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
IBStaggeredHierarchyIntegrator::getPressureVar()
{
    return d_P_var;
}// getPressureVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >
IBStaggeredHierarchyIntegrator::getForceVar()
{
    return d_F_var;
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
    const int    U_idx = var_db->mapVariableAndContextToIndex(   d_U_var, ctx);
    const int U_cc_idx = var_db->mapVariableAndContextToIndex(d_U_cc_var, ctx);
    static const bool synch_cf_interface = true;
    d_hier_math_ops->interp(
        U_cc_idx, d_U_cc_var,
        U_idx   , d_U_var   ,
        d_no_fill_op, d_integrator_time, synch_cf_interface);
    return;
}// reinterpolateVelocity

void
IBStaggeredHierarchyIntegrator::reinterpolateForce(
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx)
{
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int    F_idx = var_db->mapVariableAndContextToIndex(   d_F_var, ctx);
    const int F_cc_idx = var_db->mapVariableAndContextToIndex(d_F_cc_var, ctx);
    static const bool synch_cf_interface = true;
    d_hier_math_ops->interp(
        F_cc_idx, d_F_cc_var,
        F_idx   , d_F_var   ,
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

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_data =
        patch->getPatchData(d_U_var, ctx);
    const SAMRAI::hier::IntVector<NDIM>& U_ghost_cells = U_data->getGhostCellWidth();

    double stable_dt = std::numeric_limits<double>::max();
#if (NDIM == 2)
    NAVIER_STOKES_SC_STABLEDT_F77(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),
        U_ghost_cells(0),U_ghost_cells(1),
        U_data->getPointer(0),U_data->getPointer(1),
        stable_dt);
#endif
#if (NDIM == 3)
    NAVIER_STOKES_SC_STABLEDT_F77(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
        U_ghost_cells(0),U_ghost_cells(1),U_ghost_cells(2),
        U_data->getPointer(0),U_data->getPointer(1),U_data->getPointer(2),
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

    d_num_cycles = db->getIntegerWithDefault("num_cycles", d_num_cycles);

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
            d_Omega_rel_thresh = db->getDoubleArray("vorticity_rel_thresh");
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n"
                         << "  Vorticity tagging is enabled but key data `vorticity_rel_thresh' not found in input.\n"
                         << "  Using default values for all levels in the locally refined grid.\n");
            d_Omega_rel_thresh.resizeArray(1);
            d_Omega_rel_thresh[0] = 1.0;
        }

        for (int i = 0; i < d_Omega_rel_thresh.getSize(); ++i)
        {
            if (d_Omega_rel_thresh[i] < 0.0 || d_Omega_rel_thresh[i] > 1.0)
            {
                TBOX_ERROR(d_object_name << ":  "
                           << "relative vorticity thresholds for each level must lie in the interval [0,1].\n");
            }
        }

        if (db->keyExists("vorticity_abs_thresh"))
        {
            d_Omega_abs_thresh = db->getDoubleArray("vorticity_abs_thresh");
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n"
                         << "  Vorticity tagging is enabled but key data `vorticity_abs_thresh' not found in input.\n"
                         << "  Using default values for all levels in the locally refined grid.\n");
            d_Omega_abs_thresh.resizeArray(1);
            d_Omega_abs_thresh[0] = std::numeric_limits<double>::max();
        }

        for (int i = 0; i < d_Omega_abs_thresh.getSize(); ++i)
        {
            if (d_Omega_abs_thresh[i] < 0.0)
            {
                TBOX_ERROR(d_object_name << ":  "
                           << "absolute vorticity thresholds for each level must be nonnegative.\n");
            }
        }
    }

    d_output_U      = db->getBoolWithDefault("output_U"     , d_output_U     );
    d_output_P      = db->getBoolWithDefault("output_P"     , d_output_P     );
    d_output_F      = db->getBoolWithDefault("output_F"     , d_output_F     );
    d_output_Omega  = db->getBoolWithDefault("output_Omega" , d_output_Omega );
    d_output_Div_U  = db->getBoolWithDefault("output_Div_U" , d_output_Div_U );

    d_U_scale = db->getDoubleWithDefault("U_scale", d_U_scale);
    d_P_scale = db->getDoubleWithDefault("P_scale", d_P_scale);
    d_F_scale = db->getDoubleWithDefault("F_scale", d_F_scale);

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

        d_num_init_cycles = db->getIntegerWithDefault(
            "num_init_cycles", d_num_init_cycles);

        d_second_order_pressure_update = db->getBoolWithDefault(
            "second_order_pressure_update", d_second_order_pressure_update);

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

        d_nu = d_mu/d_rho;  // the kinematic viscosity

        if (db->keyExists("lambda"))
        {
            d_lambda = db->getDouble("lambda");
        }
        else
        {
            d_lambda = 0.0;
        }

        d_lambda = d_lambda/d_rho;

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
