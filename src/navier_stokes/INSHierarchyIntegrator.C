// Filename: INSHierarchyIntegrator.C
// Last modified: <18.Feb.2007 22:31:01 boyce@boyce-griffiths-powerbook-g4-15.local>
// Created on 02 Apr 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "INSHierarchyIntegrator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// STOOLS INCLUDES
#include <stools/PatchMathOps.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellDataFactory.h>
#include <CellIterator.h>
#include <CoarseFineBoundary.h>
#include <CoarsenOperator.h>
#include <HierarchyDataOpsManager.h>
#include <Index.h>
#include <IntVector.h>
#include <PatchCellDataOpsReal.h>
#include <PatchData.h>
#include <RefineOperator.h>
#include <SimpleCellRobinBcCoefs.h>
#include <VariableDatabase.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_ADVECTIVE_DIVSOURCE_F77 F77_FUNC_(navier_stokes_advective_divsource2d, NAVIER_STOKES_ADVECTIVE_DIVSOURCE2D)
#define NAVIER_STOKES_CONSERVATIVE_DIVSOURCE_F77 F77_FUNC_(navier_stokes_conservative_divsource2d, NAVIER_STOKES_CONSERVATIVE_DIVSOURCE2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_ADVECTIVE_DIVSOURCE_F77 F77_FUNC_(navier_stokes_advective_divsource3d, NAVIER_STOKES_ADVECTIVE_DIVSOURCE3D)
#define NAVIER_STOKES_CONSERVATIVE_DIVSOURCE_F77 F77_FUNC_(navier_stokes_conservative_divsource3d, NAVIER_STOKES_CONSERVATIVE_DIVSOURCE3D)
#endif

// Function interfaces
extern "C" {
    void NAVIER_STOKES_ADVECTIVE_DIVSOURCE_F77(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
#endif
        const double* ,
        double*);

    void NAVIER_STOKES_CONSERVATIVE_DIVSOURCE_F77(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
#endif
        const double* ,
        double*);
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
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_rebalance_coarsest_level;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_regrid_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_predict_advection_velocity;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_integrate_adv_diff;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_project_velocity;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_update_pressure;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_new_levels;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_time_dependent_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_data_to_preadvance_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_hierarchy_configuration;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_apply_gradient_detector;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_div_source_term;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int FACEG = 1;

// Version of INSHierarchyIntegrator restart file data.
static const int INS_HIERARCHY_INTEGRATOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSHierarchyIntegrator::INSHierarchyIntegrator(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    SAMRAI::tbox::Pointer<GodunovAdvector> explicit_predictor,
    SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator,
    SAMRAI::tbox::Pointer<HierarchyProjector> hier_projector,
    bool register_for_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
    assert(!hierarchy.isNull());
    assert(!explicit_predictor.isNull());
    assert(!adv_diff_hier_integrator.isNull());
    assert(!hier_projector.isNull());
#endif
    d_object_name = object_name;
    d_registered_for_restart = register_for_restart;

    d_hierarchy = hierarchy;

    d_explicit_predictor = explicit_predictor;
    d_adv_diff_hier_integrator = adv_diff_hier_integrator;
    d_hier_projector = hier_projector;

    d_hyp_level_integrator = d_adv_diff_hier_integrator->
        getHyperbolicLevelIntegrator();

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }

    // Set some default values.
    d_P_scale = 1.0;
    d_F_scale = 1.0;
    d_Q_scale = 1.0;

    d_start_time = 0.0;
    d_end_time = std::numeric_limits<double>::max();
    d_grow_dt = 2.0;
    d_max_integrator_steps = std::numeric_limits<int>::max();

    d_using_synch_projection = true;

    d_using_vorticity_tagging = false;
    d_Omega_max = 0.0;

    d_project_predicted_flux = false;

    d_second_order_pressure_update = true;

    d_num_cycles = 1;
    d_num_init_cycles = 5;

    d_cycle = 0;
    d_performing_init_cycles = false;

    d_regrid_interval = 1;
    d_old_dt = -1.0;
    d_integrator_time = std::numeric_limits<double>::quiet_NaN();
    d_integrator_step = std::numeric_limits<int>::max();

    d_conservation_form = false;

    d_output_P = false;
    d_output_F = false;
    d_output_Q = false;

    d_output_Omega = false;

    d_output_Div_U = false;
    d_output_Div_u = false;
    d_output_Div_u_adv = false;

    d_rho = std::numeric_limits<double>::quiet_NaN();
    d_mu  = std::numeric_limits<double>::quiet_NaN();
    d_nu  = std::numeric_limits<double>::quiet_NaN();

    d_dt_max = std::numeric_limits<double>::max();
    d_dt_max_time_max = std::numeric_limits<double>::max();
    d_dt_max_time_min = -(d_dt_max_time_max-numeric_limits<double>::epsilon());

    d_is_initialized = false;

    d_do_log = false;

    d_reproject_after_regrid = true;

    // Initialize object with data read from the input and restart
    // databases.
    const bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager =
        SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var = new SAMRAI::pdat::CellVariable<NDIM,double>(
        "cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(
        cc_var, hierarchy);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > fc_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(
        "fc_var");
    d_hier_fc_data_ops = hier_ops_manager->getOperationsDouble(
        fc_var, hierarchy);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::advanceHierarchy()");
        t_get_stable_timestep = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::getStableTimestep()");
        t_rebalance_coarsest_level = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::rebalanceCoarsestLevel()");
        t_regrid_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::regridHierarchy()");
        t_predict_advection_velocity = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::predictAdvectionVelocity()");
        t_integrate_adv_diff = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::integrateAdvDiff()");
        t_project_velocity = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::projectVelocity()");
        t_update_pressure = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::updatePressure()");
        t_synchronize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::putToDatabase()");
        t_compute_div_source_term = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSHierarchyIntegrator::computeDivSourceTerm()");
        timers_need_init = false;
    }
    return;
}// INSHierarchyIntegrator

INSHierarchyIntegrator::~INSHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
    return;
}// ~INSHierarchyIntegrator

void
INSHierarchyIntegrator::registerVelocityInitialConditions(
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> U_init)
{
    d_U_init = U_init;
    return;
}// registerVelocityInitialConditions

void
INSHierarchyIntegrator::registerVelocityPhysicalBcCoef(
    const std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs)
{
    d_U_bc_coefs = U_bc_coefs;
    return;
}// registerVelocityPhysicalBcCoef

void
INSHierarchyIntegrator::registerPressureInitialConditions(
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> P_init)
{
    d_P_init = P_init;
    return;
}// registerPressureInitialConditions

void
INSHierarchyIntegrator::registerPressurePhysicalBcCoef(
    const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const P_bc_coef)
{
    d_P_bc_coef = P_bc_coef;
    return;
}// registerPressurePhysicalBcCoef

void
INSHierarchyIntegrator::registerForceSpecification(
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> F_set)
{
    d_F_set = F_set;
    return;
}// registerForceSpecification

void
INSHierarchyIntegrator::registerDivergenceSpecification(
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> Q_set)
{
    d_Q_set = Q_set;
    return;
}// registerDivergenceSpecification

void
INSHierarchyIntegrator::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!visit_writer.isNull());
#endif
    d_adv_diff_hier_integrator->registerVisItDataWriter(visit_writer);
    d_visit_writer = visit_writer;
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

SAMRAI::tbox::Pointer<STOOLS::HierarchyMathOps>
INSHierarchyIntegrator::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!d_hier_math_ops.isNull());
#endif
    return d_hier_math_ops;
}// getHierarchyMathOps

void
INSHierarchyIntegrator::setHierarchyMathOps(
    SAMRAI::tbox::Pointer<STOOLS::HierarchyMathOps> hier_math_ops,
    const bool manage_ops)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hier_math_ops.isNull());
#endif
    d_hier_math_ops = hier_math_ops;
    d_is_managing_hier_math_ops = manage_ops;
    return;
}// setHierarchyMathOps

bool
INSHierarchyIntegrator::isManagingHierarchyMathOps() const
{
    return d_is_managing_hier_math_ops;
}// isManagingHierarchyMathOps

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
///      getGodunovAdvector(),
///      getAdvDiffHierarchyIntegrator(),
///      getHierarchyProjector()
///
///  allow the INSHierarchyIntegrator to be used as a hierarchy
///  integrator.
///

void
INSHierarchyIntegrator::initializeHierarchyIntegrator(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg)
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

    // Initialize all variables.
    //
    // Data corresponding to U, Grad_P, and u_adv is maintained by the
    // AdvDiffHierarchyIntegrator.
    //
    // All other data is managed by the INSHierarchyIntegrator.
    d_U_var = new SAMRAI::pdat::CellVariable<NDIM,double>(/*d_object_name+"::U"*/"u",NDIM);
    d_u_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(d_object_name+"::u");
    d_u_adv_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(d_object_name+"::u_adv");

    d_P_var = new SAMRAI::pdat::CellVariable<NDIM,double>(/*d_object_name+"::P"*/"p");
    d_Grad_P_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Grad P",NDIM);

    d_Phi_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Phi");
    d_Grad_Phi_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Grad_Phi",NDIM);
    d_grad_Phi_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(d_object_name+"::grad_Phi");

    d_G_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::G",NDIM);
    d_H_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::H",NDIM);

    d_V_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::V",NDIM);

    if (!d_F_set.isNull())
    {
        d_F_var = new SAMRAI::pdat::CellVariable<NDIM,double>(/*d_object_name+"::F"*/"f",NDIM);
    }
    if (!d_Q_set.isNull())
    {
        d_Q_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Q");
        d_F_div_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::F_div",NDIM);
    }
    if (d_output_Omega || d_using_vorticity_tagging)
    {
        const int depth = (NDIM == 2) ? 1 : NDIM;
        d_Omega_var = new SAMRAI::pdat::CellVariable<NDIM,double>(/*d_object_name+"::Omega"*/"omega",depth);
#if (NDIM == 3)
        d_Omega_Norm_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::|Omega|_2");
#endif
    }
    if (d_output_Div_U)
    {
        d_Div_U_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Div U");
    }
    if (d_output_Div_u)
    {
        d_Div_u_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Div u");
    }
    if (d_output_Div_u_adv)
    {
        d_Div_u_adv_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Div u_adv");
    }

    // Create the default communication algorithms.
    d_fill_after_regrid = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    d_calgs["SYNCH_CURRENT_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    d_calgs["SYNCH_NEW_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();

    d_rscheds["NONE"] = std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >();
    d_cscheds["NONE"] = std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >();

    // Register state variables that are maintained by the
    // INSHierarchyIntegrator.

    const SAMRAI::hier::IntVector<NDIM> cell_ghosts = CELLG;
    const SAMRAI::hier::IntVector<NDIM> face_ghosts = CELLG;
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    registerVariable(d_u_current_idx, d_u_new_idx, d_u_scratch_idx,
                     d_u_var, face_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_P_current_idx, d_P_new_idx, d_P_scratch_idx,
                     d_P_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");

    if (!d_F_var.isNull())
    {
        registerVariable(d_F_current_idx, d_F_new_idx, d_F_scratch_idx,
                         d_F_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
    }
    else
    {
        d_F_current_idx = -1;
        d_F_scratch_idx = -1;
        d_F_new_idx = -1;
    }

    if (!d_Q_var.isNull())
    {
        registerVariable(d_Q_current_idx, d_Q_new_idx, d_Q_scratch_idx,
                         d_Q_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
    }
    else
    {
        d_Q_current_idx = -1;
        d_Q_scratch_idx = -1;
        d_Q_new_idx = -1;
    }

    if (!d_F_div_var.isNull())
    {
        registerVariable(d_F_div_current_idx, d_F_div_new_idx,
                         d_F_div_scratch_idx,
                         d_F_div_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
    }
    else
    {
        d_F_div_current_idx = -1;
        d_F_div_scratch_idx = -1;
        d_F_div_new_idx = -1;
    }

    if (!d_Omega_var.isNull())
    {
        registerVariable(d_Omega_current_idx, d_Omega_new_idx,
                         d_Omega_scratch_idx,
                         d_Omega_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
#if (NDIM == 3)
        registerVariable(d_Omega_Norm_idx, d_Omega_Norm_var, no_ghosts);
#endif
    }
    else
    {
        d_Omega_current_idx = -1;
        d_Omega_scratch_idx = -1;
        d_Omega_new_idx = -1;
#if (NDIM == 3)
        d_Omega_Norm_idx = -1;
#endif
    }

    if (!d_Div_U_var.isNull())
    {
        registerVariable(d_Div_U_current_idx, d_Div_U_new_idx,
                         d_Div_U_scratch_idx,
                         d_Div_U_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
    }
    else
    {
        d_Div_U_current_idx = -1;
        d_Div_U_scratch_idx = -1;
        d_Div_U_new_idx = -1;
    }

    if (!d_Div_u_var.isNull())
    {
        registerVariable(d_Div_u_current_idx, d_Div_u_new_idx,
                         d_Div_u_scratch_idx,
                         d_Div_u_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
    }
    else
    {
        d_Div_u_current_idx = -1;
        d_Div_u_scratch_idx = -1;
        d_Div_u_new_idx = -1;
    }

    if (!d_Div_u_adv_var.isNull())
    {
        registerVariable(d_Div_u_adv_current_idx, d_Div_u_adv_new_idx,
                         d_Div_u_adv_scratch_idx,
                         d_Div_u_adv_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
    }
    else
    {
        d_Div_u_adv_current_idx = -1;
        d_Div_u_adv_scratch_idx = -1;
        d_Div_u_adv_new_idx = -1;
    }

    // Register scratch variables that are maintained by the
    // INSHierarchyIntegrator.

    registerVariable(d_Phi_idx, d_Phi_var, cell_ghosts);

    registerVariable(d_Grad_Phi_idx, d_Grad_Phi_var, no_ghosts);
    registerVariable(d_grad_Phi_idx, d_grad_Phi_var,
                     d_project_predicted_flux ? face_ghosts : no_ghosts);

    registerVariable(d_G_idx, d_G_var, no_ghosts);
    registerVariable(d_H_idx, d_H_var, cell_ghosts);

    registerVariable(d_V_idx, d_V_var, cell_ghosts);

    // Register state variables that are maintained by the
    // AdvDiffHierarchyIntegrator.

    const bool u_adv_is_div_free = d_Q_set.isNull();
    d_adv_diff_hier_integrator->registerAdvectionVelocity(
        d_u_adv_var, u_adv_is_div_free);

    d_adv_diff_hier_integrator->
        registerAdvectedAndDiffusedQuantityWithSourceTerm(
            d_U_var, d_nu, d_Grad_P_var, d_conservation_form,
            d_U_init, d_U_bc_coefs,
            SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy>(NULL),
            d_project_predicted_flux ? d_grad_Phi_var : SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >(NULL));

    // Initialize the AdvDiffHierarchyIntegrator.
    //
    // NOTE: This must be done after all variables are registered.
    d_adv_diff_hier_integrator->initializeHierarchyIntegrator(d_gridding_alg);

    // Obtain the patch data descriptor indices for all variables
    // registered with the AdvDiffHierarchyIntegrator.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    d_u_adv_current_idx = var_db->mapVariableAndContextToIndex(
        d_u_adv_var, getCurrentContext());
    d_u_adv_scratch_idx = var_db->mapVariableAndContextToIndex(
        d_u_adv_var, getScratchContext());
    d_u_adv_new_idx = var_db->mapVariableAndContextToIndex(
        d_u_adv_var, getNewContext());

    d_U_current_idx = var_db->mapVariableAndContextToIndex(
        d_U_var, getCurrentContext());
    d_U_scratch_idx = var_db->mapVariableAndContextToIndex(
        d_U_var, getScratchContext());
    d_U_new_idx = var_db->mapVariableAndContextToIndex(
        d_U_var, getNewContext());

    d_Grad_P_current_idx = var_db->mapVariableAndContextToIndex(
        d_Grad_P_var, getCurrentContext());
    d_Grad_P_scratch_idx = var_db->mapVariableAndContextToIndex(
        d_Grad_P_var, getScratchContext());
    d_Grad_P_new_idx = var_db->mapVariableAndContextToIndex(
        d_Grad_P_var, getNewContext());

    // Obtain the sol, rhs, and tmp variables from the
    // AdvDiffHierarchyIntegrator.
    //
    // XXXX: Kludge!
    d_sol_var = var_db->getVariable("AdvDiffHierarchyIntegrator::sol");
    d_rhs_var = var_db->getVariable("AdvDiffHierarchyIntegrator::rhs");
    d_tmp_var = var_db->getVariable("AdvDiffHierarchyIntegrator::tmp");

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!d_sol_var.isNull());
    assert(!d_rhs_var.isNull());
    assert(!d_tmp_var.isNull());
#endif

    d_sol_idx = var_db->mapVariableAndContextToIndex(
        d_sol_var, getCurrentContext());
    d_rhs_idx = var_db->mapVariableAndContextToIndex(
        d_rhs_var, getCurrentContext());
    d_tmp_idx = var_db->mapVariableAndContextToIndex(
        d_tmp_var, getCurrentContext());

    // Register variables for plotting.
    //
    // Note that U is automatically registered for plotting by the
    // advection-diffusion integrator.
    if (!d_visit_writer.isNull())
    {
        if (d_output_P)
        {
            d_visit_writer->registerPlotQuantity(
                d_P_var->getName(), "SCALAR", d_P_current_idx, 0, d_P_scale);
        }

        if (!d_F_var.isNull() && d_output_F)
        {
            d_visit_writer->registerPlotQuantity(
                d_F_var->getName(), "VECTOR", d_F_current_idx, 0, d_F_scale);

            for (int d = 0; d < NDIM; ++d)
            {
                ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    d_F_var->getName()+stream.str(), "SCALAR",
                    d_F_current_idx, d, d_F_scale);
            }
        }

        if (!d_Q_var.isNull() && d_output_Q)
        {
            d_visit_writer->registerPlotQuantity(
                d_Q_var->getName(), "SCALAR", d_Q_current_idx, 0, d_Q_scale);
        }

        if (!d_Omega_var.isNull() && d_output_Omega)
        {
            d_visit_writer->registerPlotQuantity(
                d_Omega_var->getName(), (NDIM == 2) ? "SCALAR" : "VECTOR",
                d_Omega_current_idx);
#if (NDIM == 3)
            for (int d = 0; d < NDIM; ++d)
            {
                ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    d_Omega_var->getName()+stream.str(), "SCALAR",
                    d_Omega_current_idx, d);
            }
#endif
        }

        if (!d_Div_U_var.isNull() && d_output_Div_U)
        {
            d_visit_writer->registerPlotQuantity(
                d_Div_U_var->getName(), "SCALAR", d_Div_U_current_idx);
        }

        if (!d_Div_u_var.isNull() && d_output_Div_u)
        {
            d_visit_writer->registerPlotQuantity(
                d_Div_u_var->getName(), "SCALAR", d_Div_u_current_idx);
        }

        if (!d_Div_u_adv_var.isNull() && d_output_Div_u_adv)
        {
            d_visit_writer->registerPlotQuantity(
                d_Div_u_adv_var->getName(), "SCALAR", d_Div_u_adv_current_idx);
        }
    }

    // Create several refinement communications algorithms, used in
    // filling ghost cell data.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;

    d_ralgs["U->V::C->S::CONSTANT_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    refine_operator = grid_geom->lookupRefineOperator(
        d_U_var, "CONSTANT_REFINE");
    d_ralgs["U->V::C->S::CONSTANT_REFINE"]->
        registerRefine(d_V_idx,         // destination
                       d_U_current_idx, // source
                       d_V_idx,         // temporary work space
                       refine_operator);

    d_ralgs["P->P::C->S::CONSTANT_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    refine_operator = grid_geom->lookupRefineOperator(
        d_P_var, "CONSTANT_REFINE");
    d_ralgs["P->P::C->S::CONSTANT_REFINE"]->
        registerRefine(d_P_scratch_idx, // destination
                       d_P_current_idx, // source
                       d_P_scratch_idx, // temporary work space
                       refine_operator);

    d_ralgs["Phi->Phi::CONSTANT_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    refine_operator = grid_geom->lookupRefineOperator(
        d_Phi_var, "CONSTANT_REFINE");
    d_ralgs["Phi->Phi::CONSTANT_REFINE"]->
        registerRefine(d_Phi_idx, // destination
                       d_Phi_idx, // source
                       d_tmp_idx, // temporary work space
                       refine_operator);

    d_ralgs["grad_Phi->grad_Phi::CONSERVATIVE_LINEAR_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    refine_operator = grid_geom->lookupRefineOperator(
        d_grad_Phi_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["grad_Phi->grad_Phi::CONSERVATIVE_LINEAR_REFINE"]->
        registerRefine(d_grad_Phi_idx,      // destination
                       d_grad_Phi_idx,      // source
                       d_u_adv_scratch_idx, // temporary work space
                       refine_operator);

    d_ralgs["predictAdvectionVelocity"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    refine_operator = grid_geom->lookupRefineOperator(
        d_u_adv_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["predictAdvectionVelocity"]->
        registerRefine(d_u_adv_scratch_idx, // destination
                       d_u_adv_current_idx, // source
                       d_u_adv_scratch_idx, // temporary work space
                       refine_operator);

    refine_operator = grid_geom->lookupRefineOperator(
        d_U_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["predictAdvectionVelocity"]->
        registerRefine(d_U_scratch_idx, // destination
                       d_U_current_idx, // source
                       d_U_scratch_idx, // temporary work space
                       refine_operator);

    refine_operator = grid_geom->lookupRefineOperator(
        d_H_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["predictAdvectionVelocity"]->
        registerRefine(d_H_idx, // destination
                       d_G_idx, // source
                       d_H_idx, // temporary work space
                       refine_operator);

    // Create several coarsening communications algorithms, used in
    // synchronizing refined regions of coarse data with the
    // underlying fine data.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_operator;

    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_U_var, "CONSERVATIVE_COARSEN");
    d_calgs["SYNCH_CURRENT_STATE_DATA"]->
        registerCoarsen(d_U_current_idx, // destination
                        d_U_current_idx, // source
                        coarsen_operator);

    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_u_adv_var, "CONSERVATIVE_COARSEN");
    d_calgs["SYNCH_CURRENT_STATE_DATA"]->
        registerCoarsen(d_u_adv_current_idx, // destination
                        d_u_adv_current_idx, // source
                        coarsen_operator);

    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_U_var, "CONSERVATIVE_COARSEN");
    d_calgs["SYNCH_NEW_STATE_DATA"]->
        registerCoarsen(d_U_new_idx, // destination
                        d_U_new_idx, // source
                        coarsen_operator);

    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_u_adv_var, "CONSERVATIVE_COARSEN");
    d_calgs["SYNCH_NEW_STATE_DATA"]->
        registerCoarsen(d_u_adv_new_idx, // destination
                        d_u_adv_new_idx, // source
                        coarsen_operator);

    // Setup the Hierarchy math operations object.
    setHierarchyMathOps(d_adv_diff_hier_integrator->getHierarchyMathOps());
    d_hier_projector->setHierarchyMathOps(d_hier_math_ops);

    // Setup the FAC preconditioner and obtain the Poisson
    // specifications used for a second order update to the pressure.
    if (d_second_order_pressure_update)
    {
        assert(false); // XXXX
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
INSHierarchyIntegrator::initializeHierarchy()
{
    t_initialize_hierarchy->start();

    if (!d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchy()\n" <<
                   "  must call initializeHierarchyIntegrator() prior to calling initializeHierarchy()." << endl);
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

    // The next timestep is given by the minimum allowable timestep
    // over all levels in the patch hierarchy.
    double dt_next = getStableTimestep();

    t_initialize_hierarchy->stop();
    return dt_next;
}// initializeHierarchy

double
INSHierarchyIntegrator::advanceHierarchy(
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
    if (rebalance_coarsest) rebalanceCoarsestLevel();

    const double current_time = d_integrator_time;
    const double new_time = d_integrator_time+dt;
    const bool initial_time = SAMRAI::tbox::Utilities::deq(d_integrator_time,d_start_time);

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

    d_cycle = 0;
    d_performing_init_cycles = initial_time;

    // Set the guess for the initial pressure to zero.
    if (d_performing_init_cycles)
    {
        d_hier_cc_data_ops->setToScalar(d_P_current_idx, // data
                                        0.0);            // alpha
    }

    const int num_cycles = d_performing_init_cycles ?
        d_num_init_cycles : d_num_cycles;

    for (d_cycle = 0; d_cycle < num_cycles; ++d_cycle)
    {
        if (d_performing_init_cycles)
        {
            if (d_do_log) SAMRAI::tbox::plog << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++\n";
            if (d_do_log) SAMRAI::tbox::plog << "+\n";
            if (d_do_log) SAMRAI::tbox::plog << "+ Performing cycle " << d_cycle+1 << " of "
                                             << d_num_init_cycles << " to initialize P(n=1/2)\n";
            if (d_do_log) SAMRAI::tbox::plog << "+\n";
            if (d_do_log) SAMRAI::tbox::plog << "++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
        }

        // Solve the Navier-Stokes equations for U(n+1), u(n+1),
        // P(n+1/2).  Each of the major algorithmic steps is separated
        // into its own member function.
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): predicting advection velocity\n";
        predictAdvectionVelocity(current_time, new_time);

        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): integrating advection-diffusion equation\n";
        integrateAdvDiff(current_time, new_time);

        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): projecting intermediate velocity\n";
        projectVelocity(current_time, new_time);

        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): updating pressure\n";
        updatePressure(current_time, new_time, true);

        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): synchronizing hierarhcy\n";
        synchronizeHierarchy();

        if (d_cycle < (num_cycles-1))
        {
            // Reset data to the preadvance state.
            resetHierDataToPreadvanceState();

            // Reset time of all current data.
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                level->setTime(current_time);
            }
        }
    }

    // Reset all time dependent data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierData(new_time);

    // Regrid (when appropriate).
    const bool do_regrid = ((d_regrid_interval == 0)
                            ? false
                            : (d_integrator_step % d_regrid_interval == 0));
    if (do_regrid)
    {
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): regridding\n";
        regridHierarchy();
    }

    // Determine the next stable timestep from u(n+1).
    const double dt_next = getStableTimestep();

    t_advance_hierarchy->stop();
    return dt_next;
}// advanceHierarchy

double
INSHierarchyIntegrator::getStableTimestep()
{
    t_get_stable_timestep->start();

    const bool initial_time =
        SAMRAI::tbox::Utilities::deq(d_integrator_time, d_start_time);
    double dt_next = std::numeric_limits<double>::max();

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt_next = SAMRAI::tbox::Utilities::
            dmin(dt_next, d_hyp_level_integrator->
                 getLevelDt(level, d_integrator_time, initial_time));
    }

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

    t_get_stable_timestep->stop();
    return dt_next;
}// getStableTimestep()

bool
INSHierarchyIntegrator::atRegridPoint() const
{
    const int level_number = 0;
    return ((d_integrator_step > 0)
            && d_gridding_alg->levelCanBeRefined(level_number)
            && (d_regrid_interval == 0
                ? false
                : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double
INSHierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}// getIntegratorTime

double
INSHierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}// getStartTime

double
INSHierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}// getEndTime

int
INSHierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}// getIntegratorStep

int
INSHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}// getMaxIntegratorSteps

bool
INSHierarchyIntegrator::stepsRemaining() const
{
    return (d_integrator_step < d_max_integrator_steps);
}// stepsRemaining

const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
INSHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
INSHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

SAMRAI::tbox::Pointer<GodunovAdvector>
INSHierarchyIntegrator::getGodunovAdvector() const
{
    return d_explicit_predictor;
}// getGodunovAdvector

SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator>
INSHierarchyIntegrator::getAdvDiffHierarchyIntegrator() const
{
    return d_adv_diff_hier_integrator;
}// getAdvDiffHierarchyIntegrator

SAMRAI::tbox::Pointer<HierarchyProjector>
INSHierarchyIntegrator::getHierarchyProjector() const
{
    return d_hier_projector;
}// getHierarchyProjector

///
///  The following routines:
///
///      rebalanceCoarsestLevel(),
///      regridHierarchy(),
///      predictAdvectionVelocity(),
///      integrateAdvDiff(),
///      projectHierarchy(),
///      updatePressure(),
///      synchronizeHierarchy(),
///      synchronizeNewLevels(),
///      resetTimeDependentHierData(),
///      resetHierDataToPreadvanceState()
///
///  allow the INSHierarchyIntegrator to provide data management
///  for a time integrator which making use of this class.
///

void
INSHierarchyIntegrator::rebalanceCoarsestLevel()
{
    t_rebalance_coarsest_level->start();

    d_gridding_alg->makeCoarsestLevel(d_hierarchy,d_integrator_time);

    t_rebalance_coarsest_level->stop();
    return;
}// rebalanceCoarsestLevel

void
INSHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    const int coarsest_ln = 0;
    d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln,
                                         d_integrator_time, d_tag_buffer);

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

void
INSHierarchyIntegrator::predictAdvectionVelocity(
    const double current_time,
    const double new_time)
{
    t_predict_advection_velocity->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(current_time <= new_time);
    assert(d_end_time > d_integrator_time);
    assert(SAMRAI::tbox::Utilities::deq(d_integrator_time,current_time));
#endif

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double dt = new_time - current_time;

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
    }

    // Immediately following a regrid, we re-project the velocity.  We
    // also project the velocity at the start time.
    if ((d_using_synch_projection && d_reproject_after_regrid) ||
        (d_performing_init_cycles && d_cycle == 0))
    {
        // Interpolate U(*)->u(*).
        static const bool u_scratch_cf_bdry_synch = true;
        d_hier_math_ops->interp(
            d_u_scratch_idx, d_u_var, // u(*)
            u_scratch_cf_bdry_synch,  // synch u(*) coarse-fine bdry
            d_V_idx        , d_V_var, // U(*)
            d_rscheds["U->V::C->S::CONSTANT_REFINE"],
            current_time);            // data time

        // Project u^(*)->u(N) and re-use Phi to project U^(*)->U^(n).
        d_hier_cc_data_ops->setToScalar(d_Phi_idx, // data
                                        0.0);      // alpha

        d_hier_projector->projectHierarchy(
            d_u_current_idx, d_u_var       , // u(n)
            d_Phi_idx      , d_Phi_var     , // Phi
            d_rscheds["Phi->Phi::CONSTANT_REFINE"],
            current_time,                    // Phi data time
            d_grad_Phi_idx , d_grad_Phi_var, // grad Phi
            d_u_scratch_idx, d_u_var       , // u(n,*)
            d_rscheds["NONE"],               // no fill needed
            current_time,                    // u_scratch data time
            false,                           // no synch needed
            d_Q_current_idx, d_Q_var);       // div u(n)

        static const bool grad_Phi_cf_bdry_synch = false;
        d_hier_math_ops->interp(
            d_Grad_Phi_idx, d_Grad_Phi_var, // dst
            d_grad_Phi_idx, d_grad_Phi_var, // src
            d_rscheds["NONE"],              // no fill needed
            current_time,                   // data time
            grad_Phi_cf_bdry_synch);        // src_cf_bdry_synch

        d_hier_cc_data_ops->subtract(d_U_current_idx, // dst
                                     d_U_current_idx, // src1
                                     d_Grad_Phi_idx); // src2

        for (int ln = finest_ln; ln > coarsest_ln; --ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom =
                d_hierarchy->getGridGeometry();
            SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > coarsen_alg = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
            SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_op;

            coarsen_op = grid_geom->lookupCoarsenOperator(
                d_u_var, "CONSERVATIVE_COARSEN");
            coarsen_alg->registerCoarsen(d_u_current_idx, // destination
                                         d_u_current_idx, // source
                                         coarsen_op);

            coarsen_op = grid_geom->lookupCoarsenOperator(
                d_U_var, "CONSERVATIVE_COARSEN");
            coarsen_alg->registerCoarsen(d_U_current_idx, // destination
                                         d_U_current_idx, // source
                                         coarsen_op);

            coarsen_alg->createSchedule(
                d_hierarchy->getPatchLevel(ln-1),
                d_hierarchy->getPatchLevel(ln))->coarsenData();
        }

        d_reproject_after_regrid = false;
    }

    // Initialize the advection velocity to equal u(n).
    d_hier_fc_data_ops->copyData(d_u_adv_current_idx, d_u_current_idx);

    // Setup the forcing terms for velocity prediction.
    d_hier_math_ops->grad(
        d_Grad_P_current_idx, d_Grad_P_var, // dst
        -1.0/d_rho,                         // alpha
        d_P_scratch_idx     , d_P_var     , // src
        d_rscheds["P->P::C->S::CONSTANT_REFINE"],
        current_time);                      // data time

    if (!d_F_var.isNull())
    {
        d_hier_cc_data_ops->axpy(d_G_idx,               // dst
                                 1.0/d_rho,             // alpha
                                 d_F_current_idx,       // src1
                                 d_Grad_P_current_idx); // src2
    }
    else
    {
        d_hier_cc_data_ops->copyData(d_G_idx,               // dst
                                     d_Grad_P_current_idx); // src
    }

    if (!d_Q_var.isNull())
    {
        computeDivSourceTerm(
            d_F_div_current_idx, d_Q_current_idx, d_u_current_idx,
            coarsest_ln, finest_ln);

        d_hier_cc_data_ops->add(d_G_idx,              // dst
                                d_G_idx,              // src1
                                d_F_div_current_idx); // src2
    }

    SAMRAI::solv::PoissonSpecifications spec("spec");
    spec.setCConstant(0.0);
    spec.setDConstant(d_nu);

    for (int d = 0; d < NDIM; ++d)
    {
        std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& V_bdry_fill =
            (d == 0) ?
            d_rscheds["U->V::C->S::CONSTANT_REFINE"] :
            d_rscheds["NONE"];

        d_hier_math_ops->laplace(
            d_G_idx, d_G_var, // dst
            spec,             // Poisson spec
            d_V_idx, d_V_var, // src1
            V_bdry_fill,      // src1 bdry fill
            current_time,     // data time
            1.0,              // beta
            d_G_idx, d_G_var, // src2
            d, d, d);         // depths
    }

    // Predict the time centered advection velocity.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_U_scratch_idx, current_time);
        level->allocatePatchData(d_u_adv_scratch_idx, current_time);
        d_rscheds["predictAdvectionVelocity"][ln]->fillData(current_time);

        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_adv_current_data =
                patch->getPatchData(d_u_adv_current_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_adv_scratch_data =
                patch->getPatchData(d_u_adv_scratch_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_scratch_data =
                patch->getPatchData(d_U_scratch_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > H_data =
                patch->getPatchData(d_H_idx);

            d_explicit_predictor->predictNormalVelocityWithSourceTerm(
                *u_adv_current_data, *u_adv_scratch_data,
                *U_scratch_data, *H_data,
                *patch, dt);
        }

        level->deallocatePatchData(d_U_scratch_idx);
        level->deallocatePatchData(d_u_adv_scratch_idx);
    }

    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_CURRENT_STATE_DATA"][ln]->coarsenData();
    }

    // Project the advection velocity.
    if (!d_Q_set.isNull())
    {
        d_Q_set->setDataOnPatchHierarchy(
            d_Q_new_idx, d_Q_var, d_hierarchy, current_time+0.5*dt);
    }

    d_hier_cc_data_ops->setToScalar(d_Phi_idx, // data
                                    0.0);      // alpha

    d_hier_projector->projectHierarchy(
        d_u_adv_current_idx, d_u_adv_var   , // u(n+1/2)
        d_Phi_idx          , d_Phi_var     , // Phi
        d_rscheds["Phi->Phi::CONSTANT_REFINE"],
        current_time,                        // Phi data time
        d_grad_Phi_idx     , d_grad_Phi_var, // grad Phi
        d_u_adv_current_idx, d_u_adv_var   , // u(n+1/2,*)
        d_rscheds["NONE"],                   // no fill needed
        current_time,                        // data time
        true,                                // synch u(n+1/2,*)
        d_Q_new_idx, d_Q_var);               // div u(n+1/2)

    if (!d_Q_set.isNull())
    {
        computeDivSourceTerm(
            d_F_div_new_idx, d_Q_new_idx, d_u_adv_current_idx,
            coarsest_ln, finest_ln);
    }

    if (!d_Div_u_adv_var.isNull())
    {
        static const bool u_adv_current_cf_bdry_synch = true;
        d_hier_math_ops->div(
            d_Div_u_adv_new_idx, d_Div_u_adv_var, // dst
            1.0,                                  // alpha
            d_u_adv_current_idx, d_u_adv_var    , // src
            d_rscheds["NONE"],                    // don't need to fill u(n+1/2) data
            current_time,                         // data time
            u_adv_current_cf_bdry_synch);         // don't re-synch u(n+1/2) coarse-fine bdry

        if (d_do_log) SAMRAI::tbox::plog << "||Div u_adv||_oo = "
                                         << d_hier_cc_data_ops->maxNorm(d_Div_u_adv_new_idx, // data
                                                                        d_wgt_idx)           // vol
                                         << endl;
    }

    t_predict_advection_velocity->stop();
    return;
}// predictAdvectionVelocity

void
INSHierarchyIntegrator::integrateAdvDiff(
    const double current_time,
    const double new_time)
{
    t_integrate_adv_diff->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(current_time <= new_time);
    assert(d_end_time > d_integrator_time);
    assert(SAMRAI::tbox::Utilities::deq(d_integrator_time,current_time));
#endif

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double dt = new_time - current_time;

    if (d_project_predicted_flux)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            d_rscheds["grad_Phi->grad_Phi::CONSERVATIVE_LINEAR_REFINE"][ln]->
                fillData(current_time);
        }
    }

    // Setup time centered forcing terms.
    if (!d_F_set.isNull())
    {
        if (d_F_set->isTimeDependent())
        {
            d_F_set->setDataOnPatchHierarchy(
                d_F_new_idx, d_F_var, d_hierarchy, current_time+dt);
        }

        d_hier_cc_data_ops->axpy(d_Grad_P_current_idx,  // dst
                                 0.5/d_rho,             // alpha
                                 d_F_current_idx,       // src1
                                 d_Grad_P_current_idx); // src2

        d_hier_cc_data_ops->axpy(d_Grad_P_current_idx,  // dst
                                 0.5/d_rho,             // alpha
                                 d_F_new_idx,           // src1
                                 d_Grad_P_current_idx); // src2
    }

    if (!d_Q_set.isNull())
    {
        d_hier_cc_data_ops->add(d_Grad_P_current_idx, // dst
                                d_Grad_P_current_idx, // src1
                                d_F_div_new_idx);     // src2
    }

    // Solve the advection-diffusion equation for U^(*).
    d_adv_diff_hier_integrator->integrateHierarchy(current_time, new_time);
    d_adv_diff_hier_integrator->synchronizeHierarchy();

    // Fix current data corresponding to Grad P.
    if (!d_F_set.isNull())
    {
        d_hier_cc_data_ops->axpy(d_Grad_P_current_idx,  // dst
                                 -0.5/d_rho,            // alpha
                                 d_F_current_idx,       // src1
                                 d_Grad_P_current_idx); // src2

        d_hier_cc_data_ops->axpy(d_Grad_P_current_idx,  // dst
                                 -0.5/d_rho,            // alpha
                                 d_F_new_idx,           // src1
                                 d_Grad_P_current_idx); // src2
    }

    if (!d_Q_set.isNull())
    {
        d_hier_cc_data_ops->subtract(d_Grad_P_current_idx, // dst
                                     d_Grad_P_current_idx, // src1
                                     d_F_div_new_idx);     // src2
    }

    t_integrate_adv_diff->stop();
    return;
}// integrateAdvDiff

void
INSHierarchyIntegrator::projectVelocity(
    const double current_time,
    const double new_time)
{
    t_project_velocity->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(current_time <= new_time);
    assert(d_end_time > d_integrator_time);
    assert(SAMRAI::tbox::Utilities::deq(d_integrator_time,current_time));
#endif

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double dt = new_time - current_time;

    // Compute Q = div u(n+1).
    if (!d_Q_set.isNull())
    {
        d_Q_set->setDataOnPatchHierarchy(
            d_Q_new_idx, d_Q_var, d_hierarchy, new_time);
    }

    // Set U^{*} = U^{*} + (dt/rho) Grad P(n-1/2).
    //
    // NOTE: d_Grad_P_current_idx = -(1/rho)*Grad P(n-1/2).
    d_hier_cc_data_ops->axpy(d_U_new_idx,          // dst
                             -dt,                  // alpha
                             d_Grad_P_current_idx, // src1
                             d_U_new_idx);         // src2

    // Compute u^(*) from U^(*).
    {
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > bdry_fill_op = grid_geom->
            lookupRefineOperator(d_U_var, "CONSTANT_REFINE");
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > bdry_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
        bdry_fill_alg->registerRefine(d_V_idx,     // destination
                                      d_U_new_idx, // source
                                      d_V_idx,     // temporary work space
                                      bdry_fill_op);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            bdry_fill_alg->resetSchedule(
                d_rscheds["U->V::C->S::CONSTANT_REFINE"][ln]);
        }

        static const bool u_scratch_cf_bdry_synch = true;
        d_hier_math_ops->interp(
            d_u_scratch_idx, d_u_var, // u(*)
            u_scratch_cf_bdry_synch,  // synch u(*) coarse-fine bdry
            d_V_idx        , d_V_var, // U(*)
            d_rscheds["U->V::C->S::CONSTANT_REFINE"],
            new_time);                // data time

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            d_ralgs["U->V::C->S::CONSTANT_REFINE"]->resetSchedule(
                d_rscheds["U->V::C->S::CONSTANT_REFINE"][ln]);
        }

        // Project u^(*)->u(N+1) and re-use Phi to project U^(*)->U^(n+1).
        d_hier_cc_data_ops->scale(d_Phi_idx,        // dst
                                  dt/d_rho,         // alpha
                                  d_P_current_idx); // src

        d_hier_projector->projectHierarchy(
            d_u_new_idx    , d_u_var       , // u(n+1)
            d_Phi_idx      , d_Phi_var     , // Phi
            d_rscheds["Phi->Phi::CONSTANT_REFINE"],
            current_time,                    // Phi data time
            d_grad_Phi_idx , d_grad_Phi_var, // grad Phi
            d_u_scratch_idx, d_u_var       , // u(n+1,*)
            d_rscheds["NONE"],               // no fill needed
            current_time,                    // u_scratch data time
            false,                           // no synch needed
            d_Q_new_idx, d_Q_var);           // div u(n+1)

        static const bool grad_Phi_cf_bdry_synch = false;
        d_hier_math_ops->interp(
            d_Grad_Phi_idx, d_Grad_Phi_var, // dst
            d_grad_Phi_idx, d_grad_Phi_var, // src
            d_rscheds["NONE"],              // no fill needed
            current_time,                   // data time
            grad_Phi_cf_bdry_synch);        // src_cf_bdry_synch

        d_hier_cc_data_ops->subtract(d_U_new_idx,     // dst
                                     d_U_new_idx,     // src1
                                     d_Grad_Phi_idx); // src2
    }

    // Optionally compute some quantities.
    if (!d_Omega_var.isNull() || !d_Div_U_var.isNull())
    {
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > bdry_fill_op = grid_geom->
            lookupRefineOperator(d_U_var, "CONSTANT_REFINE");
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > bdry_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
        bdry_fill_alg->registerRefine(d_V_idx,     // destination
                                      d_U_new_idx, // source
                                      d_V_idx,     // temporary work space
                                      bdry_fill_op);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            bdry_fill_alg->resetSchedule(
                d_rscheds["U->V::C->S::CONSTANT_REFINE"][ln]);
        }

        if (!d_Omega_var.isNull() && !d_Div_U_var.isNull())
        {
            d_hier_math_ops->curl(
                d_Omega_new_idx, d_Omega_var, // dst
                d_V_idx        , d_V_var    , // src
                d_rscheds["U->V::C->S::CONSTANT_REFINE"],
                new_time);                    // data time

            d_hier_math_ops->div(
                d_Div_U_new_idx, d_Div_U_var, // dst
                1.0,                          // alpha
                d_V_idx        , d_V_var    , // src
                d_rscheds["NONE"],            // no data fill needed
                new_time);                    // data time
        }
        else if (!d_Omega_var.isNull())
        {
            d_hier_math_ops->curl(
                d_Omega_new_idx, d_Omega_var, // dst
                d_V_idx        , d_V_var    , // src
                d_rscheds["U->V::C->S::CONSTANT_REFINE"],
                new_time);                    // data time
        }
        else if (!d_Div_U_var.isNull())
        {
            d_hier_math_ops->div(
                d_Div_U_new_idx, d_Div_U_var, // dst
                1.0,                          // alpha
                d_V_idx        , d_V_var    , // src
                d_rscheds["U->V::C->S::CONSTANT_REFINE"],
                new_time);                    // data time
        }

        if (!d_Omega_var.isNull())
        {
#if (NDIM == 2)
            d_Omega_max = d_hier_cc_data_ops->maxNorm(d_Omega_new_idx);
#endif
#if (NDIM == 3)
            d_hier_math_ops->pointwise_L2Norm(
                d_Omega_Norm_idx, d_Omega_Norm_var, // dst
                d_Omega_new_idx , d_Omega_var);     // src

            d_Omega_max = d_hier_cc_data_ops->maxNorm(d_Omega_Norm_idx);
#endif
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            d_ralgs["U->V::C->S::CONSTANT_REFINE"]->resetSchedule(
                d_rscheds["U->V::C->S::CONSTANT_REFINE"][ln]);
        }
    }
    if (!d_Div_U_var.isNull())
    {
        if (d_do_log) SAMRAI::tbox::plog << "||Div U||_1  = "
                                         << d_hier_cc_data_ops->L1Norm(d_Div_U_new_idx, // data
                                                                       d_wgt_idx)       // vol
                                         << endl;
        if (d_do_log) SAMRAI::tbox::plog << "||Div U||_2  = "
                                         << d_hier_cc_data_ops->L2Norm(d_Div_U_new_idx, // data
                                                                       d_wgt_idx)       // vol
                                         << endl;
        if (d_do_log) SAMRAI::tbox::plog << "||Div U||_oo = "
                                         << d_hier_cc_data_ops->maxNorm(d_Div_U_new_idx, // data
                                                                        d_wgt_idx)       // vol
                                         << endl;
    }
    if (!d_Div_u_var.isNull())
    {
        static const bool u_new_cf_bdry_synch = true;
        d_hier_math_ops->div(
            d_Div_u_new_idx, d_Div_u_var, // dst
            1.0,                          // alpha
            d_u_new_idx    , d_u_var    , // src
            d_rscheds["NONE"],            // no fill needed
            new_time,                     // data time
            u_new_cf_bdry_synch);         // don't synch src c-f bdry

        if (d_do_log) SAMRAI::tbox::plog << "||Div u||_oo = "
                                         << d_hier_cc_data_ops->maxNorm(d_Div_u_new_idx, // data
                                                                        d_wgt_idx)       // vol
                                         << endl;
    }

    t_project_velocity->stop();
    return;
}// projectVelocity

void
INSHierarchyIntegrator::updatePressure(
    const double current_time,
    const double new_time,
    const bool override_current_pressure)
{
    t_update_pressure->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(current_time <= new_time);
    assert(d_end_time > d_integrator_time);
    assert(SAMRAI::tbox::Utilities::deq(d_integrator_time,current_time));
#endif

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double dt = new_time - current_time;

    // We have to solve an additional Helmholtz problem in order to
    // obtain a formally consistent second order timestepping scheme
    // for the pressure.
    //
    // Setup the rhs data to solve for the second order update to the
    // pressure.
    if (d_second_order_pressure_update)
    {
        d_hier_cc_data_ops->linearSum(d_Phi_idx,        // dst
                                      d_rho/dt,         // alpha
                                      d_Phi_idx,        // src1
                                      -1.0,             // beta
                                      d_P_current_idx); // src2

        d_hier_math_ops->laplace(
            d_sol_idx, d_sol_var, // dst
            *d_helmholtz1_spec,   // Poisson spec
            d_Phi_idx, d_Phi_var, // src
            d_rscheds["Phi->Phi::CONSTANT_REFINE"],
            current_time);        // data time

        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom =
            d_hierarchy->getGridGeometry();

        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > bdry_fill_op = grid_geom->
            lookupRefineOperator(d_sol_var, "CONSTANT_REFINE");
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > bdry_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
        bdry_fill_alg->registerRefine(d_sol_idx, // destination
                                      d_sol_idx, // source
                                      d_tmp_idx, // temporary work space
                                      bdry_fill_op);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            bdry_fill_alg->resetSchedule(
                d_rscheds["Phi->Phi::CONSTANT_REFINE"][ln]);
        }

        d_hier_math_ops->laplace(
            d_rhs_idx, d_rhs_var, // dst
            *d_helmholtz2_spec,   // Poisson spec
            d_sol_idx, d_sol_var, // src
            d_rscheds["Phi->Phi::CONSTANT_REFINE"],
            current_time);        // data time

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            d_ralgs["Phi->Phi::CONSTANT_REFINE"]->resetSchedule(
                d_rscheds["Phi->Phi::CONSTANT_REFINE"][ln]);
        }

        // Solve the linear system.
        d_helmholtz4_solver->solveSystem(*d_sol_vec,*d_rhs_vec);

        // Update the pressure.
        d_hier_cc_data_ops->add(d_P_new_idx,      // dst
                                d_sol_idx,        // src1
                                d_P_current_idx); // src2
    }
    else
    {
        d_hier_cc_data_ops->scale(d_P_new_idx, // dst
                                  d_rho/dt,    // alpha
                                  d_Phi_idx);  // src
    }

    // Normalize P(n+1/2) to have mean (discrete integral) zero.
    const double P_mean = (1.0/d_volume)*
        d_hier_cc_data_ops->integral(d_P_new_idx, // data
                                     d_wgt_idx);  // vol

    d_hier_cc_data_ops->addScalar(d_P_new_idx, // dst
                                  d_P_new_idx, // src
                                  -P_mean);    // alpha

    if (override_current_pressure)
    {
        d_hier_cc_data_ops->copyData(d_P_current_idx, // dst
                                     d_P_new_idx);    // src
    }

    t_update_pressure->stop();
    return;
}// updatePressure

void
INSHierarchyIntegrator::synchronizeHierarchy()
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
INSHierarchyIntegrator::synchronizeNewLevels(
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
    // We use the AdvDiffHierarchyIntegrator to handle as
    // much data management as possible.
    d_adv_diff_hier_integrator->
        synchronizeNewLevels(hierarchy, coarsest_level, finest_level,
                             sync_time, initial_time);

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
INSHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    t_reset_time_dependent_data->start();

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // We use the AdvDiffHierarchyIntegrator to handle as
    // much data management as possible.
    d_adv_diff_hier_integrator->resetTimeDependentHierData(new_time);

    // Reset the time dependent data.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Swap SAMRAI::hier::PatchData<NDIM> pointers between the current and new contexts.
    for (list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > >::const_iterator sv =
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
                assert(src_data->getBox() == dst_data->getBox());
                assert(src_data->getGhostCellWidth() ==
                       dst_data->getGhostCellWidth());
#endif

                patch->setPatchData(dst_idx, src_data);
                patch->setPatchData(src_idx, dst_data);
            }
        }
    }

    // Deallocate the scratch and new data and reset the time of the
    // current data.
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
INSHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    t_reset_data_to_preadvance_state->start();

    // We use the AdvDiffHierarchyIntegrator to handle as
    // much data management as possible.
    d_adv_diff_hier_integrator->resetHierDataToPreadvanceState();

    // Deallocate the scratch and new data and reset the time of the
    // current data.
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
///  mesh::StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
INSHierarchyIntegrator::initializeLevelData(
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
    assert(!hierarchy.isNull());
    assert((level_number >= 0)
           && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        assert(level_number == old_level->getLevelNumber());
    }
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    // We use the AdvDiffHierarchyIntegrator and
    // HierarchyProjector objects to handle as much data management as
    // possible.
    d_adv_diff_hier_integrator->
        initializeLevelData(hierarchy, level_number, init_data_time,
                            can_be_refined, initial_time, old_level,
                            allocate_data);
    d_hier_projector->
        initializeLevelData(hierarchy, level_number, init_data_time,
                            can_be_refined, initial_time, old_level,
                            allocate_data);

    // Allocate storage needed to initialize the level and fill
    // data from coarser levels in AMR hierarchy, if any.
    //
    // Since time gets set when we allocate data, re-stamp it to
    // current time if we don't need to allocate.
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

        d_fill_after_regrid->
            createSchedule(level,
                           old_level,
                           level_number-1,
                           hierarchy)->fillData(init_data_time);

        level->deallocatePatchData(d_scratch_data);
    }

    // Initialize level data at the initial time.
    if (initial_time)
    {
        // If no initialization object is provided, initialize the
        // velocity, divergance, and vorticity to zero.  Otherwise,
        // use the initialization object to set the velocity to some
        // specified value and compute the divergance and vorticity
        // corresponding to the initial velocity.
        if (d_U_init.isNull())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_current_data = patch->
                    getPatchData(d_u_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_current_data = patch->
                    getPatchData(d_U_current_idx);

                u_current_data->fillAll(0.0);
                U_current_data->fillAll(0.0);
            }
            if (!d_Omega_var.isNull())
            {
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_current_data =
                        patch->getPatchData(d_Omega_current_idx);
                    Omega_current_data->fillAll(0.0);
                }
            }
            if (!d_Div_U_var.isNull())
            {
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Div_U_current_data =
                        patch->getPatchData(d_Div_U_current_idx);
                    Div_U_current_data->fillAll(0.0);
                }
            }
            if (!d_Div_u_var.isNull())
            {
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Div_u_current_data =
                        patch->getPatchData(d_Div_u_current_idx);
                    Div_u_current_data->fillAll(0.0);
                }
            }
            if (!d_Div_u_adv_var.isNull())
            {
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Div_u_adv_current_data =
                        patch->getPatchData(d_Div_u_adv_current_idx);
                    Div_u_adv_current_data->fillAll(0.0);
                }
            }
        }
        else
        {
            d_U_init->setDataOnPatchLevel(
                d_U_current_idx, d_U_var, level,
                init_data_time, initial_time);

            level->allocatePatchData(d_U_scratch_idx, init_data_time);

            // Fill in U boundary data from coarser levels.
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom =
                d_hierarchy->getGridGeometry();
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator =
                grid_geom->lookupRefineOperator(
                    d_U_var, "CONSERVATIVE_LINEAR_REFINE");
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > ralg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
            ralg->registerRefine(d_U_scratch_idx, // destination
                                 d_U_current_idx, // source
                                 d_U_scratch_idx, // temporary work space
                                 refine_operator);
            ralg->createSchedule(level, level_number-1, hierarchy)->
                fillData(init_data_time);

            STOOLS::PatchMathOps patch_math_ops;

            const int nghost = 1;
            std::vector<SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM> > > cf_boundary;
            cf_boundary.clear();
            cf_boundary.resize(level_number+1);

            for (int ln = 0; ln <= level_number; ++ln)
            {
                cf_boundary[ln] = new SAMRAI::hier::CoarseFineBoundary<NDIM>(*hierarchy,ln,nghost);
            }

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_current_data =
                    patch->getPatchData(d_u_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_adv_current_data =
                    patch->getPatchData(d_u_adv_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_scratch_data =
                    patch->getPatchData(d_U_scratch_idx);

                patch_math_ops.interp(
                    u_current_data, U_scratch_data,
                    patch, cf_boundary[level_number], level->getRatioToCoarserLevel());

                u_adv_current_data->copy(*u_current_data);
            }

            if (!d_Omega_var.isNull())
            {
                SAMRAI::math::PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;

                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_current_data =
                        patch->getPatchData(d_Omega_current_idx);
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_scratch_data =
                        patch->getPatchData(d_U_scratch_idx);

                    patch_math_ops.curl(
                        Omega_current_data, U_scratch_data,
                        patch, cf_boundary[level_number], level->getRatioToCoarserLevel());

#if (NDIM == 2)
                    d_Omega_max = SAMRAI::tbox::Utilities::dmax(
                        d_Omega_max,
                        patch_cc_data_ops.max(Omega_current_data, patch->getBox()));

                    d_Omega_max = SAMRAI::tbox::Utilities::dmax(
                        d_Omega_max,
                        -patch_cc_data_ops.min(Omega_current_data, patch->getBox()));
#endif
#if (NDIM == 3)
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_Norm_data =
                        new SAMRAI::pdat::CellData<NDIM,double>(patch->getBox(), 1, 0);

                    patch_math_ops.pointwise_L2Norm(
                        Omega_Norm_data, Omega_current_data,
                        patch);

                    d_Omega_max = SAMRAI::tbox::Utilities::dmax(
                        d_Omega_max,
                        patch_cc_data_ops.max(Omega_Norm_data, patch->getBox()));
#endif
                }
            }
            if (!d_Div_U_var.isNull())
            {
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Div_U_current_data =
                        patch->getPatchData(d_Div_U_current_idx);
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_scratch_data =
                        patch->getPatchData(d_U_scratch_idx);

                    patch_math_ops.div(
                        Div_U_current_data,
                        1.0, U_scratch_data,
                        0.0, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >(NULL),
                        patch, cf_boundary[level_number], level->getRatioToCoarserLevel());
                }
            }
            if (!d_Div_u_var.isNull())
            {
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Div_u_current_data =
                        patch->getPatchData(d_Div_u_current_idx);
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_current_data =
                        patch->getPatchData(d_u_current_idx);

                    patch_math_ops.div(
                        Div_u_current_data,
                        1.0, u_current_data,
                        0.0, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >(NULL),
                        patch, cf_boundary[level_number], level->getRatioToCoarserLevel());
                }
            }
            if (!d_Div_u_adv_var.isNull())
            {
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Div_u_adv_current_data =
                        patch->getPatchData(d_Div_u_adv_current_idx);
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_adv_current_data =
                        patch->getPatchData(d_u_adv_current_idx);

                    patch_math_ops.div(
                        Div_u_adv_current_data,
                        1.0, u_adv_current_data,
                        0.0, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >(NULL),
                        patch, cf_boundary[level_number], level->getRatioToCoarserLevel());
                }
            }

            level->deallocatePatchData(d_U_scratch_idx);
        }

        // If no initialization object is provided, initialize the
        // pressure to zero.  Otherwise, use the initialization object
        // to set the velocity to some specified value.
        //
        // NOTE: This initial value for the pressure IS NOT USED by
        // the time integrator and is only specified for visualization
        // purposes.  The computed pressure is initialized by
        // iterating the solution method over the initial timestep.
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
            d_P_init->setDataOnPatchLevel(
                d_P_current_idx, d_P_var, level,
                init_data_time, initial_time);
        }

        // If an initialization object is provided, initialize the
        // applied force.
        if (!d_F_set.isNull())
        {
            d_F_set->setDataOnPatchLevel(
                d_F_current_idx, d_F_var, level,
                init_data_time, initial_time);
        }

        // If an initialization object is provided, initialize the
        // divergence.
        if (!d_Q_set.isNull())
        {
            d_Q_set->setDataOnPatchLevel(
                d_Q_current_idx, d_Q_var, level,
                init_data_time, initial_time);

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_div_current_data =
                    patch->getPatchData(d_F_div_current_idx);
                F_div_current_data->fillAll(0.0);
            }
        }

        // Initialize Grad P to equal zero.
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Grad_P_current_data =
                patch->getPatchData(d_Grad_P_current_idx);
            Grad_P_current_data->fillAll(0.0);
        }
    }

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
INSHierarchyIntegrator::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((coarsest_level >= 0)
           && (coarsest_level <= finest_level)
           && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        assert(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // We use the AdvDiffHierarchyIntegrator and
    // HierarchyProjector objects to handle as much data management as
    // possible.
    d_adv_diff_hier_integrator->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    d_hier_projector          ->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Indicate that the velocity field needs to be re-projected (but
    // only in the multi-level case).
    d_reproject_after_regrid = d_reproject_after_regrid || (finest_level>0);

    // Reset the Hierarchy data operations for the new hierarchy
    // configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_fc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_fc_data_ops->resetLevels(0, finest_hier_level);

    // Reset the Hierarchy math operations for the new configuration.
    if (d_is_managing_hier_math_ops)
    {
        d_hier_math_ops->setPatchHierarchy(hierarchy);
        d_hier_math_ops->resetLevels(0, finest_hier_level);
    }

    // Get the cell weight variable and patch data descriptor index.
    d_wgt_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Get the volume of the physical domain.
    d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    // Reset the solution and rhs vectors.
    d_sol_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::sol_vec", d_hierarchy, 0, finest_hier_level);
    d_sol_vec->addComponent(d_sol_var,d_sol_idx,d_wgt_idx,d_hier_cc_data_ops);

    d_rhs_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::rhs_vec", d_hierarchy, 0, finest_hier_level);
    d_rhs_vec->addComponent(d_rhs_var,d_rhs_idx,d_wgt_idx,d_hier_cc_data_ops);

    // If we have added or removed a level, resize the schedule
    // vectors.
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
                createSchedule(coarser_level, level);
        }
    }

    // Reset the "empty" schedule vectors.
    d_rscheds["NONE"]=vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(finest_hier_level+1);
    d_cscheds["NONE"]=vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >(finest_hier_level+1);

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
INSHierarchyIntegrator::applyGradientDetector(
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
    // the AdvDiffHierarchyIntegrator.
    d_adv_diff_hier_integrator->
        applyGradientDetector(hierarchy, level_number, error_data_time,
                              tag_index, initial_time,
                              uses_richardson_extrapolation_too);

    // Tag cells based on the magnatude of the vorticity.
    if (d_using_vorticity_tagging)
    {
        const double Omega_eps =
            (level_number >= 0 && level_number < d_Omega_eps.getSize()
             ? d_Omega_eps[level_number]
             : (level_number < 0
                ? d_Omega_eps[0]
                : d_Omega_eps[d_Omega_eps.size()-1]));
        if (Omega_eps > 0.0)
        {
            const double tol = Omega_eps*d_Omega_max + sqrt(std::numeric_limits<double>::epsilon());
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
                    if (SAMRAI::tbox::Utilities::dabs((*Omega_current_data)(i)) > tol)
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
                    if (norm_Omega > tol)
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
///      getAdvectionVelocityVar(),
///      getForceVar(),
///      getDivergenceVar()
///
///  allows access to the various state variables maintained by
///  the integrator.
///

SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
INSHierarchyIntegrator::getVelocityVar()
{
    return d_U_var;
}// getVelocityVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
INSHierarchyIntegrator::getPressureVar()
{
    return d_P_var;
}// getPressureVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >
INSHierarchyIntegrator::getAdvectionVelocityVar()
{
    return d_u_adv_var;
}// getAdvectionVelocityVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
INSHierarchyIntegrator::getForceVar()
{
    return d_F_var;
}// getForceVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
INSHierarchyIntegrator::getDivergenceVar()
{
    return d_Q_var;
}// getDivergenceVar

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
/// We simply reuse the SAMRAI::hier::VariableContext objects defined
/// in the AdvDiffIntegrator object.
///

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSHierarchyIntegrator::getCurrentContext() const
{
    return d_adv_diff_hier_integrator->getCurrentContext();
}// getCurrentContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSHierarchyIntegrator::getNewContext() const
{
    return d_adv_diff_hier_integrator->getNewContext();
}// getNewContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSHierarchyIntegrator::getOldContext() const
{
    return d_adv_diff_hier_integrator->getOldContext();
}// getOldContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSHierarchyIntegrator::getScratchContext() const
{
    return d_adv_diff_hier_integrator->getScratchContext();
}// getScratchContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSHierarchyIntegrator::getPlotContext() const
{
    return d_adv_diff_hier_integrator->getPlotContext();
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
INSHierarchyIntegrator::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif
    db->putInteger("INS_HIERARCHY_INTEGRATOR_VERSION",
                   INS_HIERARCHY_INTEGRATOR_VERSION);
    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_num_cycles", d_num_cycles);
    db->putInteger("d_num_init_cycles", d_num_init_cycles);
    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putBool("d_using_default_tag_buffer", d_using_default_tag_buffer);
    db->putIntegerArray("d_tag_buffer", d_tag_buffer);
    db->putBool("d_using_synch_projection", d_using_synch_projection);
    db->putBool("d_conservation_form", d_conservation_form);
    db->putBool("d_using_vorticity_tagging", d_using_vorticity_tagging);
    if (d_using_vorticity_tagging) db->putDoubleArray("d_Omega_eps", d_Omega_eps);
    db->putDouble("d_Omega_max", d_Omega_max);
    db->putBool("d_project_predicted_flux", d_project_predicted_flux);
    db->putBool("d_second_order_pressure_update", d_second_order_pressure_update);
    db->putDouble("d_old_dt", d_old_dt);
    db->putDouble("d_integrator_time", d_integrator_time);
    db->putInteger("d_integrator_step", d_integrator_step);
    db->putDouble("d_dt_max", d_dt_max);
    db->putDouble("d_dt_max_time_max", d_dt_max_time_max);
    db->putDouble("d_dt_max_time_min", d_dt_max_time_min);
    db->putDouble("d_rho", d_rho);
    db->putDouble("d_mu", d_mu);
    db->putDouble("d_nu", d_nu);

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
INSHierarchyIntegrator::printClassData(
    std::ostream& os) const
{
    os << "\nINSHierarchyIntegrator::printClassData..." << endl;
    os << "\nINSHierarchyIntegrator: this = " << const_cast<INSHierarchyIntegrator*>(this) << endl;
    os << "d_object_name = " << d_object_name << endl;
    os << "d_integrator_time = " << d_integrator_time << "\n"
       << "d_start_time = " << d_start_time << "\n"
       << "d_end_time = " << d_end_time << "\n"
       << "d_integrator_step = " << d_integrator_step << "\n"
       << "d_grow_dt = " << d_grow_dt << endl;
    os << "d_rho = " << d_rho << "\n"
       << "d_mu = " << d_mu << "\n"
       << "d_nu = " << d_nu << endl;
    os << "I AM INCOMPLETE!!!!!!!!!" << endl;
    return;
}// printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

void
INSHierarchyIntegrator::registerVariable(
    int& current_idx,
    int& new_idx,
    int& scratch_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
    const SAMRAI::hier::IntVector<NDIM>& scratch_ghosts,
    const std::string& coarsen_name,
    const std::string& refine_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!variable.isNull());
#endif
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    current_idx = -1;  // insure that uninitialized variable
    new_idx     = -1;  // patch data descriptor indices
    scratch_idx = -1;  // cause errors

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

    // Setup the refine algorithm used to fill data in new or modified
    // patch levels following a regrid operation.
    if (!refine_operator.isNull())
    {
        d_fill_after_regrid->
            registerRefine(current_idx, // destination
                           current_idx, // source
                           scratch_idx, // temporary work space
                           refine_operator);
    }

    // Setup the SYNCH_CURRENT_STATE_DATA and SYNCH_NEW_STATE_DATA
    // algorithms, used to synchronize the data on the hierarchy.
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
INSHierarchyIntegrator::registerVariable(
    int& scratch_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
    const SAMRAI::hier::IntVector<NDIM>& scratch_ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!variable.isNull());
#endif

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    scratch_idx = -1; // insure that uninitialized variable
                      // patch data descriptor indices
                      // cause errors

    d_scratch_variables.push_back(variable);

    // Setup the scratch context.
    scratch_idx = var_db->registerVariableAndContext(
        variable, getScratchContext(), scratch_ghosts);
    d_scratch_data.setFlag(scratch_idx);

    return;
}// registerVariable

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSHierarchyIntegrator::computeDivSourceTerm(
    const int F_idx,
    const int Q_idx,
    const int u_idx,
    const int coarsest_ln,
    const int finest_ln)
{
    t_compute_div_source_term->start();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

            const SAMRAI::hier::Index<NDIM>& ilower = patch->getBox().lower();
            const SAMRAI::hier::Index<NDIM>& iupper = patch->getBox().upper();

            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_data = patch->getPatchData(u_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Q_data = patch->getPatchData(Q_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_data = patch->getPatchData(F_idx);

            const SAMRAI::hier::IntVector<NDIM>& u_data_gc = u_data->getGhostCellWidth();
            const SAMRAI::hier::IntVector<NDIM>& Q_data_gc = Q_data->getGhostCellWidth();
            const SAMRAI::hier::IntVector<NDIM>& F_data_gc = F_data->getGhostCellWidth();

            if (d_conservation_form)
            {
                NAVIER_STOKES_ADVECTIVE_DIVSOURCE_F77(
#if (NDIM == 2)
                    ilower(0),iupper(0),ilower(1),iupper(1),
                    u_data_gc(0),u_data_gc(1),
                    Q_data_gc(0),Q_data_gc(1),
                    F_data_gc(0),F_data_gc(1),
                    u_data->getPointer(0),u_data->getPointer(1),
#endif
#if (NDIM == 3)
                    ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                    u_data_gc(0),u_data_gc(1),u_data_gc(2),
                    Q_data_gc(0),Q_data_gc(1),Q_data_gc(2),
                    F_data_gc(0),F_data_gc(1),F_data_gc(2),
                    u_data->getPointer(0),u_data->getPointer(1),u_data->getPointer(2),
#endif
                    Q_data->getPointer(),
                    F_data->getPointer());
            }
            else
            {
                NAVIER_STOKES_CONSERVATIVE_DIVSOURCE_F77(
#if (NDIM == 2)
                    ilower(0),iupper(0),ilower(1),iupper(1),
                    u_data_gc(0),u_data_gc(1),
                    Q_data_gc(0),Q_data_gc(1),
                    F_data_gc(0),F_data_gc(1),
                    u_data->getPointer(0),u_data->getPointer(1),
#endif
#if (NDIM == 3)
                    ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                    u_data_gc(0),u_data_gc(1),u_data_gc(2),
                    Q_data_gc(0),Q_data_gc(1),Q_data_gc(2),
                    F_data_gc(0),F_data_gc(1),F_data_gc(2),
                    u_data->getPointer(0),u_data->getPointer(1),u_data->getPointer(2),
#endif
                    Q_data->getPointer(),
                    F_data->getPointer());
            }
        }
    }

    t_compute_div_source_term->stop();
    return;
}// computeDivSourceTerm

void
INSHierarchyIntegrator::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
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

    d_using_synch_projection = db->getBoolWithDefault(
        "using_synch_projection", d_using_synch_projection);

    d_conservation_form = db->getBoolWithDefault(
        "conservation_form", d_conservation_form);

    d_using_vorticity_tagging = db->getBoolWithDefault(
        "using_vorticity_tagging", d_using_vorticity_tagging);

    if (d_using_vorticity_tagging)
    {
        if (db->keyExists("vorticity_threshold"))
        {
            d_Omega_eps = db->getDoubleArray("vorticity_threshold");
        }
        else
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `vorticity_threshold' not found in input.");
        }

        for (int i = 0; i < d_Omega_eps.getSize(); ++i)
        {
            if (d_Omega_eps[i] < 0.0 || d_Omega_eps[i] > 1.0)
            {
                TBOX_ERROR(d_object_name << ":  "
                           << "Vorticity thesholds for each level must lie in the interval [0,1].\n");
            }
        }
    }

    d_output_P = db->getBoolWithDefault("output_P", d_output_P);
    d_output_F = db->getBoolWithDefault("output_F", d_output_F);
    d_output_Q = db->getBoolWithDefault("output_Q", d_output_Q);

    d_P_scale = db->getDoubleWithDefault("P_scale", d_P_scale);
    d_F_scale = db->getDoubleWithDefault("F_scale", d_F_scale);
    d_Q_scale = db->getDoubleWithDefault("Q_scale", d_Q_scale);

    d_output_Omega = db->getBoolWithDefault("output_Omega", d_output_Omega);

    d_output_Div_U = db->getBoolWithDefault("output_Div_U", d_output_Div_U);
    d_output_Div_u = db->getBoolWithDefault("output_Div_u", d_output_Div_u);
    d_output_Div_u_adv = db->getBoolWithDefault("output_Div_u_adv", d_output_Div_u_adv);

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

        d_project_predicted_flux = db->getBoolWithDefault(
            "project_predicted_flux", d_project_predicted_flux);

        d_second_order_pressure_update = db->getBoolWithDefault(
            "second_order_pressure_update", d_second_order_pressure_update);

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
    }

    return;
}// getFromInput

void
INSHierarchyIntegrator::getFromRestart()
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

    int ver = db->getInteger("INS_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != INS_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    d_grow_dt = db->getDouble("d_grow_dt");
    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");
    d_num_cycles = db->getInteger("d_num_cycles");
    d_num_init_cycles = db->getInteger("d_num_init_cycles");
    d_regrid_interval = db->getInteger("d_regrid_interval");
    d_using_default_tag_buffer = db->getBool("d_using_default_tag_buffer");
    d_tag_buffer = db->getIntegerArray("d_tag_buffer");
    d_using_synch_projection = db->getBool("d_using_synch_projection");
    d_conservation_form = db->getBool("d_conservation_form");
    d_using_vorticity_tagging = db->getBool("d_using_vorticity_tagging");
    if (d_using_vorticity_tagging) d_Omega_eps = db->getDoubleArray("d_Omega_eps");
    d_Omega_max = db->getDouble("d_Omega_max");
    d_project_predicted_flux = db->getBool("d_project_predicted_flux");
    d_second_order_pressure_update = db->getBool("d_second_order_pressure_update");
    d_old_dt = db->getDouble("d_old_dt");
    d_integrator_time = db->getDouble("d_integrator_time");
    d_integrator_step = db->getInteger("d_integrator_step");
    d_dt_max = db->getDouble("d_dt_max");
    d_dt_max_time_max = db->getDouble("d_dt_max_time_max");
    d_dt_max_time_min = db->getDouble("d_dt_max_time_min");
    d_rho = db->getDouble("d_rho");
    d_mu = db->getDouble("d_mu");
    d_nu = db->getDouble("d_nu");

    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
