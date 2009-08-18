// Filename: INSStaggeredHierarchyIntegrator.C
// Last modified: <17.Aug.2009 16:20:56 griffith@boyce-griffiths-mac-pro.local>
// Created on 20 Mar 2008 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "INSStaggeredHierarchyIntegrator.h"

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
#include <ibamr/INSStaggeredIntermediateVelocityBcCoef.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/INSStaggeredProjectionBcCoef.h>
#include <ibamr/INSStaggeredVelocityBcCoef.h>

// IBTK INCLUDES
#include <ibtk/CartSideDoubleDivPreservingRefine.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/FACPreconditionerLSWrapper.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/RefinePatchStrategySet.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <HierarchyDataOpsManager.h>
#include <PatchCellDataOpsReal.h>
#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <iterator>
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_SC_REGRID_COPY_FC FC_FUNC_(navier_stokes_sc_regrid_copy2d,NAVIER_STOKES_SC_REGRID_COPY2D)
#define NAVIER_STOKES_SC_REGRID_APPLY_CORRECTION_FC FC_FUNC_(navier_stokes_sc_regrid_apply_correction2d,NAVIER_STOKES_SC_REGRID_APPLY_CORRECTION2D)
#define NAVIER_STOKES_SC_STABLEDT_FC FC_FUNC_(navier_stokes_sc_stabledt2d, NAVIER_STOKES_SC_STABLEDT2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_SC_REGRID_COPY_FC FC_FUNC_(navier_stokes_sc_regrid_copy3d,NAVIER_STOKES_SC_REGRID_COPY3D)
#define NAVIER_STOKES_SC_REGRID_APPLY_CORRECTION_FC FC_FUNC_(navier_stokes_sc_regrid_apply_correction3d,NAVIER_STOKES_SC_REGRID_APPLY_CORRECTION3D)
#define NAVIER_STOKES_SC_STABLEDT_FC FC_FUNC_(navier_stokes_sc_stabledt3d, NAVIER_STOKES_SC_STABLEDT3D)
#endif

extern "C"
{
    void
    NAVIER_STOKES_SC_REGRID_COPY_FC(
        double* u_dst0, double* u_dst1,
#if (NDIM == 3)
        double* u_dst2,
#endif
        const int& u_dst_gcw,
        const double* u_src0, const double* u_src1,
#if (NDIM == 3)
        const double* u_src2,
#endif
        const int& u_src_gcw,
        const int* indicator,
        const int& indicator_gcw,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1
#if (NDIM == 3)
        ,const int& ilower2,const int& iupper2
#endif
                                     );

    void
    NAVIER_STOKES_SC_REGRID_APPLY_CORRECTION_FC(
        double* u0, double* u1,
#if (NDIM == 3)
        double* u2,
#endif
        const int& u_gcw,
        const int* indicator,
        const int& indicator_gcw,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
#if (NDIM == 3)
        const int& ilower2, const int& iupper2,
#endif
        const double* dx_fine);

    void
    NAVIER_STOKES_SC_STABLEDT_FC(
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
int
register_stokes_solver_options(
    const std::string& stokes_prefix,
    double& div_u_abstol,
    const double & default_div_u_abstol)
{
    PetscErrorCode ierr;
    ierr = PetscOptionsBegin(PETSC_COMM_WORLD, stokes_prefix.c_str(), "additional options for incompressible Stokes solver", "");  CHKERRQ(ierr);
    ierr = PetscOptionsReal("-div_u_atol", "absolute solver congergence tolerance for the value of ||div u||_oo", "", default_div_u_abstol, &div_u_abstol, PETSC_NULL);  CHKERRQ(ierr);
    ierr = PetscOptionsEnd();  CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// register_stokes_solver_options

// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy_integrator;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_advance_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_get_stable_timestep;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_regrid_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_integrate_hierarchy_initialize;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_integrate_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_integrate_hierarchy_finalize;
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
static const std::string CELL_DATA_COARSEN_TYPE = "CUBIC_COARSEN";
static const std::string SIDE_DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Version of INSStaggeredHierarchyIntegrator restart file data.
static const int INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION = 1;
}

#define OUTPUT_REGRID_DIV_U 0

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredHierarchyIntegrator::INSStaggeredHierarchyIntegrator(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    bool register_for_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!hierarchy.isNull());
#endif
    d_object_name = object_name;
    d_registered_for_restart = register_for_restart;

    d_hierarchy = hierarchy;

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Set some default values.
    d_U_scale = 1.0;
    d_P_scale = 1.0;
    d_F_scale = 1.0;

    d_start_time = 0.0;
    d_end_time = std::numeric_limits<double>::max();
    d_grow_dt = 2.0;
    d_max_integrator_steps = std::numeric_limits<int>::max();

    d_num_cycles = 3;

    d_using_vorticity_tagging = false;
    d_Omega_max = 0.0;

    d_normalize_pressure = false;

    d_regrid_interval = 1;
    d_old_dt = -1.0;
    d_op_and_solver_init_dt = -1.0;
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
    d_lambda = std::numeric_limits<double>::quiet_NaN();

    d_cfl = 0.9;

    d_dt_max = std::numeric_limits<double>::max();
    d_dt_max_time_max = std::numeric_limits<double>::max();
    d_dt_max_time_min = -(d_dt_max_time_max-std::numeric_limits<double>::epsilon());

    d_is_initialized = false;

    d_do_log = false;

    d_regrid_max_div_growth_factor = 1.1;

    // Setup default boundary condition objects that specify homogeneous
    // Dirichlet boundary conditions for the velocity.
    d_default_U_bc_coef = new SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>(
        d_object_name+"::default_U_bc_coef",
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL));
    for (int d = 0; d < NDIM; ++d)
    {
        d_default_U_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_U_bc_coef->setBoundaryValue(2*d+1,0.0);
    }
    registerVelocityPhysicalBcCoefs(std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(
                                        NDIM,d_default_U_bc_coef));

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Set the problem coefs.
    d_problem_coefs = new INSCoefs(d_rho, d_mu, d_lambda);

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager = SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var = new SAMRAI::pdat::CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > sc_var = new SAMRAI::pdat::SideVariable<NDIM,double>("sc_var");
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, hierarchy);

    // Initialize all variable contexts.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    d_current_context = var_db->getContext(d_object_name+"::CURRENT");
    d_new_context     = var_db->getContext(d_object_name+"::NEW"    );
    d_scratch_context = var_db->getContext(d_object_name+"::SCRATCH");

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::advanceHierarchy()");
        t_get_stable_timestep = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::getStableTimestep()");
        t_regrid_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::regridHierarchy()");
        t_integrate_hierarchy_initialize = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::integrateHierarchy_initialize()");
        t_integrate_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::integrateHierarchy()");
        t_integrate_hierarchy_finalize = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::integrateHierarchy_finalize()");
        t_synchronize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// INSStaggeredHierarchyIntegrator

INSStaggeredHierarchyIntegrator::~INSStaggeredHierarchyIntegrator()
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

    if (d_helmholtz_spec != NULL) delete d_helmholtz_spec;
    if (d_poisson_spec != NULL) delete d_poisson_spec;
    delete d_default_U_bc_coef;
    if (!d_U_bc_coefs.empty())
    {
        for (int d = 0; d < NDIM; ++d)
        {
            delete d_U_bc_coefs[d];
            delete d_U_star_bc_coefs[d];
        }
        delete d_P_bc_coef;
        delete d_Phi_bc_coef;
    }
    delete d_problem_coefs;
    if (d_regrid_projection_spec != NULL) delete d_regrid_projection_spec;
    return;
}// ~INSStaggeredHierarchyIntegrator

const std::string&
INSStaggeredHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
INSStaggeredHierarchyIntegrator::registerVelocityInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> U_init)
{
    d_U_init = U_init;
    return;
}// registerVelocityInitialConditions

void
INSStaggeredHierarchyIntegrator::registerVelocityPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::registerVelocityPhysicalBcCoefs():\n"
                   << "  velocity boundary conditions must be registered prior to initialization\n"
                   << "  of the hierarchy integrator object." << std::endl);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned l = 0; l < U_bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(U_bc_coefs[l] != NULL);
    }
#endif
    if (!d_U_bc_coefs.empty())
    {
        for (int d = 0; d < NDIM; ++d)
        {
            delete d_U_bc_coefs[d];
            delete d_U_star_bc_coefs[d];
        }
        delete d_P_bc_coef;
        delete d_Phi_bc_coef;
    }
    d_U_bc_coefs.clear();
    d_U_bc_coefs.resize(NDIM,NULL);
    for (int d = 0; d < NDIM; ++d)
    {
        d_U_bc_coefs[d] = new INSStaggeredVelocityBcCoef(d,*d_problem_coefs,U_bc_coefs);
    }
    d_U_star_bc_coefs.clear();
    d_U_star_bc_coefs.resize(NDIM,NULL);
    for (int d = 0; d < NDIM; ++d)
    {
        d_U_star_bc_coefs[d] = new INSStaggeredIntermediateVelocityBcCoef(d,U_bc_coefs);
    }
    d_P_bc_coef = new INSStaggeredPressureBcCoef(*d_problem_coefs,U_bc_coefs);
    d_Phi_bc_coef = new INSStaggeredProjectionBcCoef(U_bc_coefs);
    return;
}// registerVelocityPhysicalBcCoefs

void
INSStaggeredHierarchyIntegrator::registerPressureInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> P_init)
{
    d_P_init = P_init;
    return;
}// registerPressureInitialConditions

void
INSStaggeredHierarchyIntegrator::registerBodyForceSpecification(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> F_set)
{
    d_F_set = F_set;
    return;
}// registerBodyForceSpecification

void
INSStaggeredHierarchyIntegrator::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    return;
}// registerVisItDataWriter

void
INSStaggeredHierarchyIntegrator::registerRegridHierarchyCallback(
    void (*callback)(const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy, const double regrid_data_time, const bool initial_time, void* ctx),
    void* ctx)
{
    d_regrid_hierarchy_callbacks.push_back(callback);
    d_regrid_hierarchy_callback_ctxs.push_back(ctx);
    return;
}// registerRegridHierarchyCallback

void
INSStaggeredHierarchyIntegrator::registerApplyGradientDetectorCallback(
    void (*callback)(const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy, const int level_number, const double error_data_time, const int tag_index, const bool initial_time, const bool uses_richardson_extrapolation_too, void* ctx),
    void* ctx)
{
    d_apply_gradient_detector_callbacks.push_back(callback);
    d_apply_gradient_detector_callback_ctxs.push_back(ctx);
    return;
}// registerApplyGradientDetectorCallback

///
///  The following routines:
///
///      getHierarchyMathOps(),
///      setHierarchyMathOps(),
///      isManagingHierarchyMathOps()
///
///  allow for the sharing of a single HierarchyMathOps object between mutiple
///  HierarchyIntegrator objects.
///

SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps>
INSStaggeredHierarchyIntegrator::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hier_math_ops.isNull());
#endif
    return d_hier_math_ops;
}// getHierarchyMathOps

void
INSStaggeredHierarchyIntegrator::setHierarchyMathOps(
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
    const bool manage_ops)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hier_math_ops.isNull());
#endif
    d_hier_math_ops = hier_math_ops;
    d_is_managing_hier_math_ops = manage_ops;
    return;
}// setHierarchyMathOps

bool
INSStaggeredHierarchyIntegrator::isManagingHierarchyMathOps() const
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
///      getGriddingAlgorithm()
///
///  allow the INSStaggeredHierarchyIntegrator to be used as a hierarchy integrator.
///

void
INSStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

    int ierr;
    PetscTruth flg;

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
    d_Phi_var        = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Phi"        );
    d_indicator_var  = new SAMRAI::pdat::CellVariable<NDIM,int   >(d_object_name+"::indicator"  );

    // Create the default communication algorithms.
    d_fill_after_regrid = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    d_calgs["SYNCH_CURRENT_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    d_calgs["SYNCH_NEW_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();

    // Register state variables that are maintained by the
    // INSStaggeredHierarchyIntegrator.

    const SAMRAI::hier::IntVector<NDIM> cell_ghosts = CELLG;
    const SAMRAI::hier::IntVector<NDIM> side_ghosts = SIDEG;
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

    // Register scratch variables that are maintained by the
    // INSStaggeredHierarchyIntegrator.

    registerVariable(d_Phi_idx, d_Phi_var, cell_ghosts);
    registerVariable(d_indicator_idx, d_indicator_var, side_ghosts);

    // Setup regridding refine algorithm for resetting velocity values along the
    // coarse-fine interface.
    d_fill_cf_interface_after_regrid = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    d_fill_cf_interface_after_regrid->registerRefine(d_U_current_idx, // destination
                                                     d_U_current_idx, // source
                                                     d_U_scratch_idx, // temporary work space
                                                     NULL);
    d_fill_cf_interface_after_regrid->registerRefine(d_indicator_idx, // destination
                                                     d_indicator_idx, // source
                                                     d_indicator_idx, // temporary work space
                                                     NULL);

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

    // Set the current integration time.
    if (!SAMRAI::tbox::RestartManager::getManager()->isFromRestart())
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }

    // Setup the Hierarchy math operations object.
    d_hier_math_ops = new IBTK::HierarchyMathOps(d_object_name+"::HierarchyMathOps", d_hierarchy);
    d_is_managing_hier_math_ops = true;

    // Setup physical boundary conditions helper.
    d_U_bc_helper = new INSStaggeredPhysicalBoundaryHelper();

    // Setup the Stokes operator.
    d_stokes_op = new INSStaggeredStokesOperator(
        *d_problem_coefs,
        d_U_bc_coefs, d_U_bc_helper, d_P_bc_coef,
        d_hier_math_ops);

    // Setup the convective operator.
    d_convective_op_needs_init = true;
    d_convective_op = new INSStaggeredPPMConvectiveOperator(
        *d_problem_coefs,
        d_conservation_form);

    // Setup the linear solver.
    const std::string stokes_prefix = "stokes_";
    d_stokes_solver_needs_init = true;
    d_stokes_solver = new IBTK::PETScKrylovLinearSolver(d_object_name+"::stokes_solver", stokes_prefix);
    d_stokes_solver->setInitialGuessNonzero(true);
    d_stokes_solver->setOperator(d_stokes_op);
    d_stokes_solver->setKSPType("fgmres");
    d_div_u_abstol = 1.0e-5;
    ierr = register_stokes_solver_options(stokes_prefix, d_div_u_abstol, d_div_u_abstol);  IBTK_CHKERRQ(ierr);

    // Setup the preconditioner and preconditioner sub-solvers.
    std::vector<std::string> pc_shell_types(3);
    pc_shell_types[0] = "projection";
    pc_shell_types[1] = "block_factorization";
    pc_shell_types[2] = "none";
    d_stokes_solver->setValidPCShellTypes(pc_shell_types);

    size_t len = 255;
    char stokes_pc_type_str[len];
    ierr = PetscOptionsGetString("stokes_", "-pc_type", stokes_pc_type_str, len, &flg);  IBTK_CHKERRQ(ierr);
    std::string stokes_pc_type = "shell";
    if (flg)
    {
        stokes_pc_type = std::string(stokes_pc_type_str);
    }

    if (!(stokes_pc_type == "none" || stokes_pc_type == "shell"))
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                   "  invalid stokes preconditioner type: " << stokes_pc_type << "\n"
                   "  valid stokes preconditioner types: shell, none" << std::endl);
    }

    std::string stokes_pc_shell_type;
    if (stokes_pc_type == "shell")
    {
        char stokes_pc_shell_type_str[len];
        ierr = PetscOptionsGetString("stokes_", "-pc_shell_type", stokes_pc_shell_type_str, len, &flg);  IBTK_CHKERRQ(ierr);
        stokes_pc_shell_type = "projection";
        if (flg)
        {
            stokes_pc_shell_type = std::string(stokes_pc_shell_type_str);
        }

        if (!(stokes_pc_shell_type == "none" || stokes_pc_shell_type == "projection" || stokes_pc_shell_type == "block_factorization"))
        {
            TBOX_ERROR(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                       "  invalid stokes shell preconditioner type: " << stokes_pc_shell_type << "\n"
                       "  valid stokes shell preconditioner types: projection, block_factorization, none" << std::endl);
        }
    }
    else
    {
        stokes_pc_shell_type = "none";
    }

    // Setup the velocity subdomain solver.
    const bool needs_helmholtz_solver = stokes_pc_type == "shell" && (stokes_pc_shell_type == "projection" || stokes_pc_shell_type == "block_factorization");
    if (needs_helmholtz_solver)
    {
        const std::string helmholtz_prefix = "helmholtz_";

        // Setup the various solver components.
        d_helmholtz_spec = new SAMRAI::solv::PoissonSpecifications(d_object_name+"::helmholtz_spec");
        d_helmholtz_op = new IBTK::SCLaplaceOperator(d_object_name+"::Helmholtz Operator", *d_helmholtz_spec, d_U_star_bc_coefs, true);
        d_helmholtz_op->setHierarchyMathOps(d_hier_math_ops);

        d_helmholtz_solver_needs_init = true;
        d_helmholtz_solver = new IBTK::PETScKrylovLinearSolver(d_object_name+"::Helmholtz Krylov Solver", helmholtz_prefix);
        d_helmholtz_solver->setInitialGuessNonzero(false);
        d_helmholtz_solver->setOperator(d_helmholtz_op);

        if (d_gridding_alg->getMaxLevels() == 1)
        {
            if (d_helmholtz_hypre_pc_db.isNull())
            {
                TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                             "  helmholtz hypre pc solver database is null." << std::endl);
            }
            d_helmholtz_hypre_pc = new IBTK::SCPoissonHypreLevelSolver(d_object_name+"::Helmholtz Preconditioner", d_helmholtz_hypre_pc_db);
            d_helmholtz_hypre_pc->setPoissonSpecifications(*d_helmholtz_spec);

            d_helmholtz_solver->setPreconditioner(d_helmholtz_hypre_pc);
        }
        else
        {
            if (d_helmholtz_fac_pc_db.isNull())
            {
                TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                             "  helmholtz fac pc solver database is null." << std::endl);
            }
            d_helmholtz_fac_op = new IBTK::SCPoissonFACOperator(d_object_name+"::Helmholtz FAC Operator", d_helmholtz_fac_pc_db);
            d_helmholtz_fac_op->setPoissonSpecifications(*d_helmholtz_spec);

            d_helmholtz_fac_pc = new SAMRAI::solv::FACPreconditioner<NDIM>(d_object_name+"::Helmholtz Preconditioner", *d_helmholtz_fac_op, d_helmholtz_fac_pc_db);
            d_helmholtz_fac_op->setPreconditioner(d_helmholtz_fac_pc);

            d_helmholtz_solver->setPreconditioner(new IBTK::FACPreconditionerLSWrapper(d_helmholtz_fac_pc, d_helmholtz_fac_pc_db));
        }

        // Set some default options.
        d_helmholtz_solver->setKSPType(d_gridding_alg->getMaxLevels() == 1 ? "preonly" : "gmres");
        d_helmholtz_solver->setAbsoluteTolerance(1.0e-30);
        d_helmholtz_solver->setRelativeTolerance(1.0e-02);
        d_helmholtz_solver->setMaxIterations(25);
    }
    else
    {
        d_helmholtz_spec = NULL;
        d_helmholtz_op = NULL;
        d_helmholtz_hypre_pc = NULL;
        d_helmholtz_fac_op = NULL;
        d_helmholtz_fac_pc = NULL;
        d_helmholtz_solver = NULL;
    }

    // Setup the pressure subdomain solver.
    const bool needs_poisson_solver = stokes_pc_type == "shell" && (stokes_pc_shell_type == "projection" || stokes_pc_shell_type == "block_factorization");
    if (needs_poisson_solver)
    {
        const std::string poisson_prefix = "poisson_";

        // Setup the various solver components.
        d_poisson_spec = new SAMRAI::solv::PoissonSpecifications(d_object_name+"::poisson_spec");
        d_poisson_op = new IBTK::CCLaplaceOperator(d_object_name+"::Poisson Operator", *d_poisson_spec, d_U_star_bc_coefs, true);
        d_poisson_op->setHierarchyMathOps(d_hier_math_ops);

        d_poisson_solver_needs_init = true;
        d_poisson_solver = new IBTK::PETScKrylovLinearSolver(d_object_name+"::Poisson Krylov Solver", poisson_prefix);
        d_poisson_solver->setInitialGuessNonzero(false);
        d_poisson_solver->setOperator(d_poisson_op);

        if (d_gridding_alg->getMaxLevels() == 1)
        {
            if (d_poisson_hypre_pc_db.isNull())
            {
                TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                             "  poisson hypre pc solver database is null." << std::endl);
            }
            d_poisson_hypre_pc = new IBTK::CCPoissonHypreLevelSolver(d_object_name+"::Poisson Preconditioner", d_poisson_hypre_pc_db);
            d_poisson_hypre_pc->setPoissonSpecifications(*d_poisson_spec);

            d_poisson_solver->setPreconditioner(d_poisson_hypre_pc);
        }
        else
        {
            if (d_poisson_fac_pc_db.isNull())
            {
                TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                             "  poisson fac pc solver database is null." << std::endl);
            }
            d_poisson_fac_op = new IBTK::CCPoissonFACOperator(d_object_name+"::Poisson FAC Operator", d_poisson_fac_pc_db);
            d_poisson_fac_op->setPoissonSpecifications(*d_poisson_spec);

            d_poisson_fac_pc = new SAMRAI::solv::FACPreconditioner<NDIM>(d_object_name+"::Poisson Preconditioner", *d_poisson_fac_op, d_poisson_fac_pc_db);
            d_poisson_fac_op->setPreconditioner(d_poisson_fac_pc);

            d_poisson_solver->setPreconditioner(new IBTK::FACPreconditionerLSWrapper(d_poisson_fac_pc, d_poisson_fac_pc_db));
        }

        // Set some default options.
        d_poisson_solver->setKSPType(d_gridding_alg->getMaxLevels() == 1 ? "preonly" : "gmres");
        d_poisson_solver->setAbsoluteTolerance(1.0e-30);
        d_poisson_solver->setRelativeTolerance(1.0e-02);
        d_poisson_solver->setMaxIterations(25);
        const bool constant_null_space = d_normalize_pressure;
        if (constant_null_space)
        {
            std::string iname = std::string("-") + poisson_prefix + std::string("ksp_constant_null_space");
            ierr = PetscOptionsSetValue(iname.c_str(), PETSC_NULL);  IBTK_CHKERRQ(ierr);
        }
    }
    else
    {
        d_poisson_spec = NULL;
        d_poisson_op = NULL;
        d_poisson_hypre_pc = NULL;
        d_poisson_fac_op = NULL;
        d_poisson_fac_pc = NULL;
        d_poisson_solver = NULL;
    }

    // Setup the Stokes preconditioner.
    if (stokes_pc_type == "shell")
    {
        if (stokes_pc_shell_type == "projection")
        {
            d_projection_pc_needs_init = true;
            d_projection_pc = new INSStaggeredProjectionPreconditioner(*d_problem_coefs, d_Phi_bc_coef, d_normalize_pressure, d_helmholtz_solver, d_poisson_solver, d_hier_cc_data_ops, d_hier_sc_data_ops, d_hier_math_ops);
            d_stokes_solver->setPreconditioner(d_projection_pc);
        }
        else if (stokes_pc_shell_type == "block_factorization")
        {
            d_block_pc_needs_init = true;
            d_block_pc = new INSStaggeredBlockFactorizationPreconditioner(*d_problem_coefs, d_Phi_bc_coef, d_normalize_pressure, d_helmholtz_solver, d_poisson_solver, d_hier_cc_data_ops, d_hier_sc_data_ops, d_hier_math_ops);
            d_stokes_solver->setPreconditioner(d_block_pc);
        }
    }

    // Setup the regrid projection Poisson solver.
    const bool needs_regrid_projection_solver = (d_gridding_alg->getMaxLevels() > 1);
    if (needs_regrid_projection_solver)
    {
        const std::string regrid_projection_prefix = "regrid_projection_";

        // Setup the various solver components.
        for (int d = 0; d < NDIM; ++d)
        {
            d_regrid_projection_bc_coef.setBoundarySlope(2*d  ,0.0);
            d_regrid_projection_bc_coef.setBoundarySlope(2*d+1,0.0);
        }

        d_regrid_projection_spec = new SAMRAI::solv::PoissonSpecifications(d_object_name+"::regrid_projection_spec");
        d_regrid_projection_op = new IBTK::CCLaplaceOperator(d_object_name+"::Regrid Projection Poisson Operator", *d_regrid_projection_spec, &d_regrid_projection_bc_coef, true);
        d_regrid_projection_op->setHierarchyMathOps(d_hier_math_ops);

        d_regrid_projection_solver = new IBTK::PETScKrylovLinearSolver(d_object_name+"::Regrid Projection Poisson Krylov Solver", regrid_projection_prefix);
        d_regrid_projection_solver->setInitialGuessNonzero(false);
        d_regrid_projection_solver->setOperator(d_regrid_projection_op);

        TBOX_ASSERT(d_gridding_alg->getMaxLevels() > 1);

        if (d_regrid_projection_fac_pc_db.isNull())
        {
            TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                         "  regrid projection poisson fac pc solver database is null." << std::endl);
        }
        d_regrid_projection_fac_op = new IBTK::CCPoissonFACOperator(d_object_name+"::Regrid Projection Poisson FAC Operator", d_regrid_projection_fac_pc_db);
        d_regrid_projection_fac_op->setPoissonSpecifications(*d_regrid_projection_spec);

        d_regrid_projection_fac_pc = new SAMRAI::solv::FACPreconditioner<NDIM>(d_object_name+"::Regrid Projection Poisson Preconditioner", *d_regrid_projection_fac_op, d_regrid_projection_fac_pc_db);
        d_regrid_projection_fac_op->setPreconditioner(d_regrid_projection_fac_pc);

        d_regrid_projection_solver->setPreconditioner(new IBTK::FACPreconditionerLSWrapper(d_regrid_projection_fac_pc, d_regrid_projection_fac_pc_db));

        // Set some default options.
        d_regrid_projection_solver->setKSPType("gmres");
        d_regrid_projection_solver->setAbsoluteTolerance(1.0e-12);
        d_regrid_projection_solver->setRelativeTolerance(1.0e-08);
        d_regrid_projection_solver->setMaxIterations(25);

        // NOTE: We always use homogeneous Neumann boundary conditions for the
        // regrid projection Poisson solver.
        static const bool constant_null_space = true;
        if (constant_null_space)
        {
            std::string iname = std::string("-") + regrid_projection_prefix + std::string("ksp_constant_null_space");
            ierr = PetscOptionsSetValue(iname.c_str(), PETSC_NULL);  IBTK_CHKERRQ(ierr);
        }
    }
    else
    {
        d_regrid_projection_spec = NULL;
        d_regrid_projection_op = NULL;
        d_regrid_projection_fac_op = NULL;
        d_regrid_projection_fac_pc = NULL;
        d_regrid_projection_solver = NULL;
    }

    // Indicate that the integrator has been initialized.
    d_is_initialized = true;

    t_initialize_hierarchy_integrator->stop();
    return;
}// initializeHierarchyIntegrator

double
INSStaggeredHierarchyIntegrator::initializeHierarchy()
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
        d_gridding_alg->getTagAndInitializeStrategy()->resetHierarchyConfiguration(
            d_hierarchy, coarsest_ln, finest_ln);
    }
    else
    {
        d_gridding_alg->makeCoarsestLevel(d_hierarchy,d_start_time);

        int level_number = 0;
        bool done = false;
        while (!done && (d_gridding_alg->levelCanBeRefined(level_number)))
        {
            d_gridding_alg->makeFinerLevel(
                d_hierarchy, d_integrator_time, initial_time, d_tag_buffer[level_number]);

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
            synchronizeNewLevels(d_hierarchy, coarsest_ln, finest_ln, d_start_time, initial_time);
        }
    }

    // The next timestep is given by the minimum allowable timestep over all
    // levels in the patch hierarchy.
    const double dt_next = getStableTimestep(getCurrentContext());

    // Initialize the operators and solvers.
    const double current_time = initial_time ? d_start_time : d_integrator_time;
    const double new_time = current_time + dt_next;
    initializeOperatorsAndSolvers(current_time, new_time);

    t_initialize_hierarchy->stop();
    return dt_next;
}// initializeHierarchy

double
INSStaggeredHierarchyIntegrator::advanceHierarchy(
    const double dt)
{
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

    // Integrate all time dependent data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): integrating time dependent data\n";
    integrateHierarchy_initialize(current_time, new_time);
    for (int cycle = 0; cycle < d_num_cycles; ++cycle)
    {
        integrateHierarchy(current_time, new_time);
    }
    integrateHierarchy_finalize(current_time, new_time);

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
INSStaggeredHierarchyIntegrator::getStableTimestep(
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

    if (d_integrator_time >= d_dt_max_time_min && d_integrator_time <= d_dt_max_time_max)
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
INSStaggeredHierarchyIntegrator::atRegridPoint() const
{
    const int level_number = 0;
    return ((d_integrator_step > 0)
            && d_gridding_alg->levelCanBeRefined(level_number)
            && (d_regrid_interval == 0
                ? false
                : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double
INSStaggeredHierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}// getIntegratorTime

double
INSStaggeredHierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}// getStartTime

double
INSStaggeredHierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}// getEndTime

int
INSStaggeredHierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}// getIntegratorStep

int
INSStaggeredHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}// getMaxIntegratorSteps

bool
INSStaggeredHierarchyIntegrator::stepsRemaining() const
{
    return (d_integrator_step < d_max_integrator_steps);
}// stepsRemaining

const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
INSStaggeredHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
INSStaggeredHierarchyIntegrator::getGriddingAlgorithm() const
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
///  allow the INSStaggeredHierarchyIntegrator to provide data management for a time
///  integrator which making use of this class.
///

void
INSStaggeredHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    const bool initial_time = SAMRAI::tbox::MathUtilities<double>::equalEps(d_integrator_time,d_start_time);

    const int coarsest_ln = 0;

    // Determine the divergence of the velocity field before regridding.
    d_hier_math_ops->div(d_Div_U_current_idx, d_Div_U_var, 1.0, d_U_current_idx, d_U_var, d_no_fill_op, d_integrator_time, false);
    const double Div_U_norm_1_pre  = d_hier_cc_data_ops->L1Norm( d_Div_U_current_idx, d_wgt_cc_idx);
    const double Div_U_norm_2_pre  = d_hier_cc_data_ops->L2Norm( d_Div_U_current_idx, d_wgt_cc_idx);
    const double Div_U_norm_oo_pre = d_hier_cc_data_ops->maxNorm(d_Div_U_current_idx, d_wgt_cc_idx);

    // Copy current data to regrid data.
    const int finest_ln_before_regrid = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln_before_regrid; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_regrid_data, d_integrator_time);
        for (std::map<int,int>::const_iterator cit = d_regrid_current_idx_map.begin();
             cit != d_regrid_current_idx_map.end(); ++cit)
        {
            const int current_idx = (*cit).first;
            const int regrid_current_idx = (*cit).second;
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(current_idx != regrid_current_idx);
#endif
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > current_data = patch->getPatchData(current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > regrid_current_data = patch->getPatchData(regrid_current_idx);
                regrid_current_data->copy(*current_data);
            }
        }
    }

    // Regrid the hierarchy
    d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln, d_integrator_time, d_tag_buffer);

    // Swap current data with regrid data.
    const int finest_ln_after_regrid = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln_after_regrid; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (std::map<int,int>::const_iterator cit = d_regrid_current_idx_map.begin();
             cit != d_regrid_current_idx_map.end(); ++cit)
        {
            const int current_idx = (*cit).first;
            const int regrid_current_idx = (*cit).second;
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(current_idx != regrid_current_idx);
#endif
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > current_data = patch->getPatchData(current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > regrid_current_data = patch->getPatchData(regrid_current_idx);
                patch->setPatchData(current_idx, regrid_current_data);
                patch->setPatchData(regrid_current_idx, current_data);
            }
        }
        level->deallocatePatchData(d_regrid_data);
    }

    // Determine the divergence of the velocity field after regridding.
    d_hier_math_ops->div(d_Div_U_current_idx, d_Div_U_var, 1.0, d_U_current_idx, d_U_var, d_no_fill_op, d_integrator_time, true);
    const double Div_U_norm_1_post  = d_hier_cc_data_ops->L1Norm( d_Div_U_current_idx, d_wgt_cc_idx);
    const double Div_U_norm_2_post  = d_hier_cc_data_ops->L2Norm( d_Div_U_current_idx, d_wgt_cc_idx);
    const double Div_U_norm_oo_post = d_hier_cc_data_ops->maxNorm(d_Div_U_current_idx, d_wgt_cc_idx);

    // Project the interpolated velocity if needed.
    if (d_needs_regrid_projection && (Div_U_norm_1_post  > d_regrid_max_div_growth_factor*Div_U_norm_1_pre ||
                                      Div_U_norm_2_post  > d_regrid_max_div_growth_factor*Div_U_norm_2_pre ||
                                      Div_U_norm_oo_post > d_regrid_max_div_growth_factor*Div_U_norm_oo_pre))
    {
        regridProjection();
    }
    d_needs_regrid_projection = false;

    // Synchronize the state data on the patch hierarchy.
    for (int ln = finest_ln_after_regrid; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_CURRENT_STATE_DATA"][ln]->coarsenData();
    }

    // Execute any registered callback functions.
    for (size_t i = 0; i < d_regrid_hierarchy_callbacks.size(); ++i)
    {
        (*d_regrid_hierarchy_callbacks[i])(d_hierarchy, d_integrator_time, initial_time, d_regrid_hierarchy_callback_ctxs[i]);
    }

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

void
INSStaggeredHierarchyIntegrator::integrateHierarchy_initialize(
    const double current_time,
    const double new_time)
{
    t_integrate_hierarchy_initialize->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::integrateHierarchy_initialize(): current_time = " << current_time << ", new_time = " << new_time << ", dt = " << dt << "\n";

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
    d_U_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(d_object_name+"::U_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
    d_U_scratch_vec->addComponent(d_U_var, d_U_scratch_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    d_U_rhs_vec = d_U_scratch_vec->cloneVector(d_object_name+"::U_rhs_vec");
    d_U_rhs_vec->allocateVectorData(current_time);
    const int U_rhs_idx = d_U_rhs_vec->getComponentDescriptorIndex(0);
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_rhs_var = d_U_rhs_vec->getComponentVariable(0);
    d_hier_sc_data_ops->setToScalar(U_rhs_idx,0.0);

    d_U_half_vec = d_U_scratch_vec->cloneVector(d_object_name+"::U_half_vec");
    d_U_half_vec->allocateVectorData(current_time);
    const int U_half_idx = d_U_half_vec->getComponentDescriptorIndex(0);
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_half_var = d_U_half_vec->getComponentVariable(0);
    d_hier_sc_data_ops->setToScalar(U_half_idx,0.0);

    d_N_vec = d_U_scratch_vec->cloneVector(d_object_name+"::N_vec");
    d_N_vec->allocateVectorData(current_time);
    const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > N_var = d_N_vec->getComponentVariable(0);
    d_hier_sc_data_ops->setToScalar(N_idx,0.0);

    d_P_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(d_object_name+"::P_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
    d_P_scratch_vec->addComponent(d_P_var, d_P_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    d_P_rhs_vec = d_P_scratch_vec->cloneVector(d_object_name+"::P_rhs_vec");
    d_P_rhs_vec->allocateVectorData(current_time);
    const int P_rhs_idx = d_P_rhs_vec->getComponentDescriptorIndex(0);
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_rhs_var = d_P_rhs_vec->getComponentVariable(0);
    d_hier_cc_data_ops->setToScalar(P_rhs_idx,0.0);

    // Initialize the right-hand side terms.
    SAMRAI::solv::PoissonSpecifications rhs_spec(d_object_name+"::rhs_spec");
    rhs_spec.setCConstant((d_rho/dt)-0.5*d_lambda);
    rhs_spec.setDConstant(          +0.5*d_mu    );
    d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
    d_hier_math_ops->laplace(
        U_rhs_idx, U_rhs_var,
        rhs_spec,
        d_U_scratch_idx, d_U_var,
        d_U_bdry_bc_fill_op, current_time);

    // Reset the solution, rhs, and nullspace vectors.
    d_sol_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    d_sol_vec->addComponent(d_U_var,d_U_scratch_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    d_sol_vec->addComponent(d_P_var,d_P_scratch_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    d_rhs_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    d_rhs_vec->addComponent(d_U_var,U_rhs_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    d_rhs_vec->addComponent(d_P_var,P_rhs_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    // Setup the operators and solvers.
    initializeOperatorsAndSolvers(current_time, new_time);

    // Setup the convergence test.
    PetscErrorCode ierr;
    KSP petsc_ksp = d_stokes_solver->getPETScKSP();
    ierr = KSPDefaultConvergedCreate(&d_default_conv_ctx);  IBTK_CHKERRQ(ierr);
    ierr = KSPSetConvergenceTest(petsc_ksp, INSStaggeredHierarchyIntegrator::KSPDivUConvergenceTest, static_cast<void*>(this), PETSC_NULL);  IBTK_CHKERRQ(ierr);

    // Setup the nullspace object.
    if (d_normalize_pressure)
    {
        d_nul_vec = d_sol_vec->cloneVector(d_object_name+"::nul_vec");
        d_nul_vec->allocateVectorData(current_time);
        d_hier_sc_data_ops->setToScalar(d_nul_vec->getComponentDescriptorIndex(0), 0.0);
        d_hier_cc_data_ops->setToScalar(d_nul_vec->getComponentDescriptorIndex(1), 1.0);

        MatNullSpace petsc_nullsp;
        Vec petsc_nullsp_vec = IBTK::PETScSAMRAIVectorReal<double>::createPETScVector(d_nul_vec, PETSC_COMM_WORLD);
        Vec vecs[] = {petsc_nullsp_vec};
        static const PetscTruth has_cnst = PETSC_FALSE;
        ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, has_cnst, 1, vecs, &petsc_nullsp); IBTK_CHKERRQ(ierr);
        KSP petsc_ksp = d_stokes_solver->getPETScKSP();
        ierr = KSPSetNullSpace(petsc_ksp, petsc_nullsp); IBTK_CHKERRQ(ierr);
    }

    // Set the initial guess.
    d_hier_sc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(0), d_U_current_idx);
    d_hier_cc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(1), d_P_current_idx);
    d_hier_sc_data_ops->copyData(d_U_new_idx, d_sol_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_P_new_idx, d_sol_vec->getComponentDescriptorIndex(1));

    // Setup inhomogeneous boundary conditions.
    d_U_bc_helper->clearBcCoefData();
    d_U_bc_helper->cacheBcCoefData(d_U_scratch_idx, d_U_var, d_U_bc_coefs, new_time, SAMRAI::hier::IntVector<NDIM>(SIDEG), d_hierarchy);

    d_stokes_op->setHomogeneousBc(false);
    d_stokes_op->modifyRhsForInhomogeneousBc(*d_rhs_vec);
    d_stokes_op->setHomogeneousBc(true);

    t_integrate_hierarchy_initialize->start();
    return;
}// integrateHierarchy_initialize

void
INSStaggeredHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time)
{
    t_integrate_hierarchy->start();

    const double dt = new_time-current_time;

    // Perform a single step of fixed point iteration.

    // Compute U_half := 0.5*(u(n)+u(n+1)).
    const int U_half_idx = d_U_half_vec->getComponentDescriptorIndex(0);
    d_hier_sc_data_ops->linearSum(U_half_idx, 0.5, d_U_current_idx, 0.5, d_U_new_idx);

    // Compute (u_half*grad)u_half.
    const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
    d_convective_op->applyConvectiveOperator(U_half_idx, N_idx);

    // Setup the right-hand side vector.
    d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), -d_rho, N_idx, d_rhs_vec->getComponentDescriptorIndex(0));
    if (!d_F_set.isNull())
    {
        d_F_set->setDataOnPatchHierarchy(d_F_scratch_idx, d_F_var, d_hierarchy, current_time+0.5*dt);
        d_hier_sc_data_ops->add(d_rhs_vec->getComponentDescriptorIndex(0), d_rhs_vec->getComponentDescriptorIndex(0), d_F_scratch_idx);
    }

    // Ensure there is no forcing at Dirichlet boundaries (the Dirichlet
    // boundary condition takes precedence).
    d_U_bc_helper->zeroValuesAtDirichletBoundaries(d_rhs_vec->getComponentDescriptorIndex(0));

    // Solve for u(n+1), p(n+1/2).
    d_stokes_solver->solveSystem(*d_sol_vec,*d_rhs_vec);

    // Enforce Dirichlet boundary conditions.
    d_U_bc_helper->resetValuesAtDirichletBoundaries(d_sol_vec->getComponentDescriptorIndex(0));

    // Pull out solution components.
    d_hier_sc_data_ops->copyData(d_U_new_idx, d_sol_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_P_new_idx, d_sol_vec->getComponentDescriptorIndex(1));

    // Reset the right-hand side vector.
    d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), +d_rho, N_idx, d_rhs_vec->getComponentDescriptorIndex(0));
    if (!d_F_set.isNull())
    {
        d_hier_sc_data_ops->subtract(d_rhs_vec->getComponentDescriptorIndex(0), d_rhs_vec->getComponentDescriptorIndex(0), d_F_scratch_idx);
        d_hier_sc_data_ops->copyData(d_F_new_idx, d_F_scratch_idx);
    }

    t_integrate_hierarchy->start();
    return;
}// integrateHierarchy

void
INSStaggeredHierarchyIntegrator::integrateHierarchy_finalize(
    const double current_time,
    const double new_time)
{
    t_integrate_hierarchy_finalize->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

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

    // Compute Omega = curl U.
    //
    // NOTE: Re-filling ghost cell data here overrides the conservative
    // coarsening of U from fine levels to coarse levels.  In particular, this
    // means that we need to re-coarsen the data associated with patch data
    // descriptor index d_U_scratch_idx if we wish to compute, e.g., the
    // divergence of the velocity field.  This operation
    d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_new_idx);
    d_hier_math_ops->curl(
        d_Omega_new_idx, d_Omega_var,
        d_U_scratch_idx, d_U_var,
        d_U_bdry_bc_fill_op, new_time);
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
#if OUTPUT_REGRID_DIV_U
    d_hier_cc_data_ops->copyData(d_Div_U_new_idx,d_Div_U_current_idx);
#else
    d_hier_math_ops->div(
        d_Div_U_new_idx, d_Div_U_var,
        1.0, d_U_new_idx, d_U_var,
        d_no_fill_op, new_time, false);
#endif

    // Deallocate the nullspace object.
    if (d_normalize_pressure)
    {
        PetscErrorCode ierr;
        MatNullSpace petsc_nullsp;
        KSP petsc_ksp = d_stokes_solver->getPETScKSP();
        ierr = KSPGetNullSpace(petsc_ksp, &petsc_nullsp); IBTK_CHKERRQ(ierr);
        ierr = MatNullSpaceDestroy(petsc_nullsp); IBTK_CHKERRQ(ierr);
    }

    // Deallocate the convergence test context.
    PetscErrorCode ierr = KSPDefaultConvergedDestroy(d_default_conv_ctx); IBTK_CHKERRQ(ierr);

    // Deallocate scratch data.
    d_U_rhs_vec->freeVectorComponents();
    d_U_half_vec->freeVectorComponents();
    d_N_vec->freeVectorComponents();
    d_P_rhs_vec->freeVectorComponents();
    if (d_normalize_pressure)
    {
        d_nul_vec->freeVectorComponents();
    }

    // Deallocate solver vectors.
    d_U_scratch_vec.setNull();
    d_U_rhs_vec.setNull();
    d_U_half_vec.setNull();
    d_N_vec.setNull();
    d_P_scratch_vec.setNull();
    d_P_rhs_vec.setNull();
    d_sol_vec.setNull();
    d_rhs_vec.setNull();
    d_nul_vec.setNull();

    t_integrate_hierarchy_finalize->stop();
    return;
}// integrateHierarchy_finalize

void
INSStaggeredHierarchyIntegrator::synchronizeHierarchy()
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
INSStaggeredHierarchyIntegrator::synchronizeNewLevels(
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
INSStaggeredHierarchyIntegrator::resetTimeDependentHierData(
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
INSStaggeredHierarchyIntegrator::resetHierDataToPreadvanceState()
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
///  mesh::StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
INSStaggeredHierarchyIntegrator::initializeLevelData(
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
    if (!initial_time && (level_number > 0 || !old_level.isNull()))
    {
        level->allocatePatchData(d_regrid_data, init_data_time);
        level->allocatePatchData(d_scratch_data, init_data_time);

        IBTK::CartExtrapPhysBdryOp fill_after_regrid_extrap_bc_op(d_fill_after_regrid_bc_idxs, BDRY_EXTRAP_TYPE);
        IBTK::CartSideRobinPhysBdryOp fill_after_regrid_phys_bdry_bc_op(d_regrid_scratch_idx_map[d_U_scratch_idx], d_U_bc_coefs, false);
        IBTK::CartSideDoubleDivPreservingRefine fill_after_regrid_div_preserving_op(d_regrid_scratch_idx_map[d_U_scratch_idx]);
        std::vector<SAMRAI::xfer::RefinePatchStrategy<NDIM>*> refine_patch_strategies(3);
        refine_patch_strategies[0] = &fill_after_regrid_extrap_bc_op;
        refine_patch_strategies[1] = &fill_after_regrid_phys_bdry_bc_op;
        refine_patch_strategies[2] = &fill_after_regrid_div_preserving_op;
        IBTK::RefinePatchStrategySet fill_after_regrid_patch_strategy_set(refine_patch_strategies.begin(), refine_patch_strategies.end(), false);
        d_fill_after_regrid->createSchedule(level,
                                            old_level,
                                            level_number-1,
                                            hierarchy,
                                            &fill_after_regrid_patch_strategy_set)->fillData(init_data_time);

        if (!old_level.isNull() && level_number == hierarchy->getFinestLevelNumber())
        {
            old_level->allocatePatchData(d_indicator_idx, init_data_time);

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(old_level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = old_level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > indicator_data = patch->getPatchData(d_indicator_idx);
                indicator_data->fillAll(1);
            }

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > indicator_data = patch->getPatchData(d_indicator_idx);
                indicator_data->fillAll(0);
            }

            IBTK::CartSideRobinPhysBdryOp fill_cf_interface_after_regrid_bc_op(d_U_scratch_idx, d_U_bc_coefs, false);
            d_fill_cf_interface_after_regrid->createSchedule(level, old_level, &fill_cf_interface_after_regrid_bc_op)->fillData(init_data_time);

            const SAMRAI::hier::IntVector<NDIM>& ratio_to_coarser_level = level->getRatioToCoarserLevel();
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const dx = patch_geom->getDx();

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_dst_data = patch->getPatchData(d_regrid_current_idx_map[d_U_current_idx]);
                const int U_dst_ghosts = U_dst_data->getGhostCellWidth().max();

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_src_data = patch->getPatchData(d_U_current_idx);
                const int U_src_ghosts = U_src_data->getGhostCellWidth().max();

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > indicator_data = patch->getPatchData(d_indicator_idx);
                const int indicator_ghosts = indicator_data->getGhostCellWidth().max();

                NAVIER_STOKES_SC_REGRID_COPY_FC(
                    U_dst_data->getPointer(0),
                    U_dst_data->getPointer(1),
#if (NDIM == 3)
                    U_dst_data->getPointer(2),
#endif
                    U_dst_ghosts,
                    U_src_data->getPointer(0),
                    U_src_data->getPointer(1),
#if (NDIM == 3)
                    U_src_data->getPointer(2),
#endif
                    U_src_ghosts,
                    indicator_data->getPointer(0),
                    indicator_ghosts,
                    patch_box.lower()(0), patch_box.upper()(0),
                    patch_box.lower()(1), patch_box.upper()(1)
#if (NDIM == 3)
                    ,patch_box.lower()(2),patch_box.upper()(2)
#endif
                                                 );

                if (ratio_to_coarser_level == SAMRAI::hier::IntVector<NDIM>(2))
                {
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(U_dst_ghosts%2 == 0);
                    TBOX_ASSERT(U_dst_ghosts <= indicator_ghosts);
#endif
                    NAVIER_STOKES_SC_REGRID_APPLY_CORRECTION_FC(
                        U_dst_data->getPointer(0),
                        U_dst_data->getPointer(1),
#if (NDIM == 3)
                        U_dst_data->getPointer(2),
#endif
                        U_dst_ghosts,
                        indicator_data->getPointer(0),
                        indicator_ghosts,
                        patch_box.lower()(0), patch_box.upper()(0),
                        patch_box.lower()(1), patch_box.upper()(1),
#if (NDIM == 3)
                        patch_box.lower()(2), patch_box.upper()(2),
#endif
                        dx                                      );
                }
            }
            old_level->deallocatePatchData(d_indicator_idx);
        }
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

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
INSStaggeredHierarchyIntegrator::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
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

    InterpolationTransactionComponent U_bc_component(d_U_scratch_idx, SIDE_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs);
    d_U_bdry_bc_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_U_bdry_bc_fill_op->initializeOperatorState(U_bc_component, d_hierarchy);

//  InterpolationTransactionComponent U_extrap_component(d_U_scratch_idx, SIDE_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);
//  d_U_bdry_extrap_fill_op = new IBTK::HierarchyGhostCellInterpolation();
//  d_U_bdry_extrap_fill_op->initializeOperatorState(U_extrap_component, d_hierarchy);

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
    d_convective_op_needs_init = true;
    d_helmholtz_solver_needs_init = true;
    d_poisson_solver_needs_init = true;
    d_projection_pc_needs_init = true;
    d_block_pc_needs_init = true;
    d_stokes_solver_needs_init = true;

    // Indicate that we need to perform a regrid projection.
    d_needs_regrid_projection = true;

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
INSStaggeredHierarchyIntegrator::applyGradientDetector(
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
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_current_data = patch->getPatchData(d_Omega_current_idx);
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

    // Allow callback functions to tag cells for refinement.
    for (size_t i = 0; i < d_apply_gradient_detector_callbacks.size(); ++i)
    {
        (*d_apply_gradient_detector_callbacks[i])(
            hierarchy, level_number, error_data_time,
            tag_index, initial_time,
            uses_richardson_extrapolation_too, d_apply_gradient_detector_callback_ctxs[i]);
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
INSStaggeredHierarchyIntegrator::getVelocityVar()
{
    return d_U_var;
}// getVelocityVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
INSStaggeredHierarchyIntegrator::getPressureVar()
{
    return d_P_var;
}// getPressureVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >
INSStaggeredHierarchyIntegrator::getForceVar()
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
INSStaggeredHierarchyIntegrator::getCurrentContext() const
{
    return d_current_context;
}// getCurrentContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSStaggeredHierarchyIntegrator::getNewContext() const
{
    return d_new_context;
}// getNewContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSStaggeredHierarchyIntegrator::getScratchContext() const
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
INSStaggeredHierarchyIntegrator::reinterpolateVelocity(
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
INSStaggeredHierarchyIntegrator::reinterpolateForce(
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
INSStaggeredHierarchyIntegrator::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION",
                   INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION);

    db->putDouble("d_U_scale", d_U_scale);
    db->putDouble("d_P_scale", d_P_scale);
    db->putDouble("d_F_scale", d_F_scale);
    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_num_cycles", d_num_cycles);
    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putBool("d_using_default_tag_buffer", d_using_default_tag_buffer);
    db->putIntegerArray("d_tag_buffer", d_tag_buffer);
    db->putBool("d_conservation_form", d_conservation_form);
    db->putBool("d_using_vorticity_tagging", d_using_vorticity_tagging);
    db->putDoubleArray("d_Omega_rel_thresh", d_Omega_rel_thresh);
    db->putDoubleArray("d_Omega_abs_thresh", d_Omega_abs_thresh);
    db->putDouble("d_Omega_max", d_Omega_max);
    db->putBool("d_normalize_pressure", d_normalize_pressure);
    db->putBool("d_output_U", d_output_U);
    db->putBool("d_output_P", d_output_P);
    db->putBool("d_output_F", d_output_F);
    db->putBool("d_output_Omega", d_output_Omega);
    db->putBool("d_output_Div_U", d_output_Div_U);
    db->putDouble("d_old_dt", d_old_dt);
    db->putDouble("d_op_and_solver_init_dt", d_op_and_solver_init_dt);
    db->putDouble("d_integrator_time", d_integrator_time);
    db->putInteger("d_integrator_step", d_integrator_step);
    db->putDouble("d_cfl", d_cfl);
    db->putDouble("d_dt_max", d_dt_max);
    db->putDouble("d_dt_max_time_max", d_dt_max_time_max);
    db->putDouble("d_dt_max_time_min", d_dt_max_time_min);
    db->putBool("d_do_log", d_do_log);
    db->putDouble("d_rho", d_rho);
    db->putDouble("d_mu", d_mu);
    db->putDouble("d_lambda", d_lambda);
    db->putDouble("d_div_u_abstol", d_div_u_abstol);
    db->putDouble("d_regrid_max_div_growth_factor", d_regrid_max_div_growth_factor);

    t_put_to_database->stop();
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

void
INSStaggeredHierarchyIntegrator::registerVariable(
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
    const bool data_lives_on_patch_border = variable->dataLivesOnPatchBorder();
    const bool fine_boundary_represents_var = variable->fineBoundaryRepresentsVariable();

    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    current_idx = -1; // insure that uninitialized variable patch data
    new_idx     = -1; // descriptor indices cause errors
    scratch_idx = -1;

    d_state_variables.push_back(variable);

     // Setup the current context.
    current_idx = var_db->registerVariableAndContext(variable, getCurrentContext(), no_ghosts);
    d_current_data.setFlag(current_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(current_idx);
    }

    // Setup the new context.
    new_idx = var_db->registerVariableAndContext(variable, getNewContext(), no_ghosts);
    d_new_data.setFlag(new_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(new_idx);
    }

    // Setup the scratch context.
    scratch_idx = var_db->registerVariableAndContext(variable, getScratchContext(), scratch_ghosts);
    d_scratch_data.setFlag(scratch_idx);

    // Setup the regrid data.
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > regrid_variable = NULL;
    int regrid_current_idx = -1;
    int regrid_scratch_idx = -1;
    const SAMRAI::hier::IntVector<NDIM> regrid_scratch_ghosts = 2;
    if ((!data_lives_on_patch_border || !fine_boundary_represents_var) && scratch_ghosts >= regrid_scratch_ghosts)
    {
        regrid_variable = variable;
        regrid_current_idx = current_idx;
        regrid_scratch_idx = scratch_idx;
    }
    else
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var = variable;
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > fc_var = variable;
        SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> > nc_var = variable;
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > sc_var = variable;
        if (!cc_var.isNull())
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > cc_pd_factory = cc_var->getPatchDataFactory();
            regrid_variable = new SAMRAI::pdat::CellVariable<NDIM,double>(cc_var->getName() + "::REGRID", cc_pd_factory->getDefaultDepth());
        }
        else if (!fc_var.isNull())
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceDataFactory<NDIM,double> > sc_pd_factory = sc_var->getPatchDataFactory();
            static const bool fine_bdry_represents_var = false;
            regrid_variable = new SAMRAI::pdat::FaceVariable<NDIM,double>(sc_var->getName() + "::REGRID", sc_pd_factory->getDefaultDepth(), fine_bdry_represents_var);
        }
        else if (!nc_var.isNull())
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeDataFactory<NDIM,double> > nc_pd_factory = nc_var->getPatchDataFactory();
            regrid_variable = new SAMRAI::pdat::NodeVariable<NDIM,double>(nc_var->getName() + "::REGRID", nc_pd_factory->getDefaultDepth());
        }
        else if (!sc_var.isNull())
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideDataFactory<NDIM,double> > sc_pd_factory = sc_var->getPatchDataFactory();
            static const bool fine_bdry_represents_var = false;
            regrid_variable = new SAMRAI::pdat::SideVariable<NDIM,double>(sc_var->getName() + "::REGRID", sc_pd_factory->getDefaultDepth(), fine_bdry_represents_var);
        }
        else
        {
            TBOX_ERROR(d_object_name << "::registerVariable():\n"
                       << "  unsupported variable type." << std::endl);
        }

        regrid_current_idx = var_db->registerVariableAndContext(regrid_variable, getCurrentContext(), no_ghosts);
        d_regrid_data.setFlag(regrid_current_idx);
        d_regrid_current_idx_map[current_idx] = regrid_current_idx;

        regrid_scratch_idx = var_db->registerVariableAndContext(regrid_variable, getScratchContext(), regrid_scratch_ghosts);
        d_regrid_data.setFlag(regrid_scratch_idx);
        d_regrid_scratch_idx_map[scratch_idx] = regrid_scratch_idx;
   }

    // Setup the refine algorithm used to fill data in new or modified patch
    // levels following a regrid operation.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > regrid_refine_operator =
        grid_geom->lookupRefineOperator(regrid_variable, refine_name);
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > regrid_coarsen_operator =
        grid_geom->lookupCoarsenOperator(regrid_variable, coarsen_name);
    if (!regrid_refine_operator.isNull())
    {
        d_fill_after_regrid_bc_idxs.setFlag(regrid_scratch_idx);
        d_fill_after_regrid->registerRefine(regrid_current_idx, // destination
                                            regrid_current_idx, // source
                                            regrid_scratch_idx, // temporary work space
                                            regrid_refine_operator);
    }

    // Setup the SYNCH_CURRENT_STATE_DATA and SYNCH_NEW_STATE_DATA algorithms,
    // used to synchronize the data on the hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator =
        grid_geom->lookupRefineOperator(variable, refine_name);
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_operator =
        grid_geom->lookupCoarsenOperator(variable, coarsen_name);
    if (!coarsen_operator.isNull())
    {
        d_calgs["SYNCH_CURRENT_STATE_DATA"]->registerCoarsen(current_idx, // destination
                                                             current_idx, // source
                                                             coarsen_operator);

        d_calgs["SYNCH_NEW_STATE_DATA"]->registerCoarsen(new_idx, // destination
                                                         new_idx, // source
                                                         coarsen_operator);
    }
    return;
}// registerVariable

void
INSStaggeredHierarchyIntegrator::registerVariable(
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

PetscErrorCode
INSStaggeredHierarchyIntegrator::KSPDivUConvergenceTest(
    KSP ksp,
    PetscInt n,
    PetscReal rnorm,
    KSPConvergedReason* reason,
    void* convergence_test_ctx)
{
    INSStaggeredHierarchyIntegrator* hier_integrator = static_cast<INSStaggeredHierarchyIntegrator*>(convergence_test_ctx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(hier_integrator != NULL);
#endif
    PetscErrorCode ierr;
    ierr = KSPDefaultConverged(ksp, n, rnorm, reason, hier_integrator->d_default_conv_ctx);  IBTK_CHKERRQ(ierr);

    // Whenever the default convergence test is satisfied, compute the discrete
    // divergence of the solution vector and ensure that it satisfies the
    // relevant convergence tolerance.
    //
    // NOTE: Here, we assume that the flow is incompressible with no internal
    // fluid sources or sinks (i.e., div u = 0 holds throughout the
    // computational domain).  This convergence test will require some minor
    // modifications to support internal fluid sources and sinks.
    const double& div_u_abstol = hier_integrator->d_div_u_abstol;
    if (int(*reason) > 0 && div_u_abstol > 0.0)
    {
        SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > hier_cc_data_ops = hier_integrator->d_hier_cc_data_ops;
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops = hier_integrator->d_hier_math_ops;
        const int wgt_cc_idx = hier_integrator->d_wgt_cc_idx;

        const int Div_U_idx = hier_integrator->d_Div_U_scratch_idx;
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Div_U_cc_var = hier_integrator->d_Div_U_var;
        const SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>& no_fill_op = hier_integrator->d_no_fill_op;
        const double& integrator_time = hier_integrator->d_integrator_time;

        Vec petsc_sol_vec;
        ierr = KSPBuildSolution(ksp, PETSC_NULL, &petsc_sol_vec);  IBTK_CHKERRQ(ierr);
        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > sol_vec = IBTK::PETScSAMRAIVectorReal<double>::getSAMRAIVector(petsc_sol_vec);
        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > sol_vec_clone = sol_vec->cloneVector("");
        sol_vec_clone->allocateVectorData(integrator_time);
        sol_vec_clone->copyVector(sol_vec);

        const int U_idx = sol_vec_clone->getComponentDescriptorIndex(0);
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_var = sol_vec_clone->getComponentVariable(0);
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_sc_var = U_var;

        hier_integrator->d_U_bc_helper->resetValuesAtDirichletBoundaries(U_idx);

        static const bool U_cf_bdry_synch = true;
        hier_math_ops->div(
            Div_U_idx, Div_U_cc_var, // dst
            +1.0,                    // alpha
            U_idx, U_sc_var,         // src
            no_fill_op,              // src_bdry_fill
            integrator_time,         // src_bdry_fill_time
            U_cf_bdry_synch);        // src_cf_bdry_synch

        sol_vec_clone->freeVectorComponents();
        sol_vec_clone->deallocateVectorData();

        const double Div_U_oo = hier_cc_data_ops->maxNorm(Div_U_idx, wgt_cc_idx);
        if (Div_U_oo > div_u_abstol)
        {
            PetscInfo3(ksp,"Linear solver has converged according to KSPDefaultConverged, but solution does not yet satisfy the absolute convergence tolerance on div U. Divergence max-norm %G is greater than the absolute tolerance %G at iteration %D\n", Div_U_oo, div_u_abstol, n);
            *reason = KSP_CONVERGED_ITERATING;
        }
        else
        {
            PetscInfo3(ksp,"Linear solver has converged according to KSPDefaultConverged, and solution satisfies the absolute convergence tolerance on div U. Divergence max-norm %G is less than or equal to the absolute tolerance %G at iteration %D\n", Div_U_oo, div_u_abstol, n);
        }
    }
    PetscFunctionReturn(0);
}// KSPDivUConvergenceTest

void
INSStaggeredHierarchyIntegrator::regridProjection()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate temporary data.
    SAMRAI::hier::ComponentSelector scratch_idxs;
    scratch_idxs.setFlag(d_U_scratch_idx);
    scratch_idxs.setFlag(d_Phi_idx);
    scratch_idxs.setFlag(d_Div_U_scratch_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(scratch_idxs, d_integrator_time);
    }

    // Compute div U before applying the projection operator.
    const bool U_current_cf_bdry_synch = true;
    d_hier_math_ops->div(
        d_Div_U_scratch_idx, d_Div_U_var, // dst
        +1.0,                             // alpha
        d_U_current_idx, d_U_var,         // src
        d_no_fill_op,                     // src_bdry_fill
        d_integrator_time,                // src_bdry_fill_time
        U_current_cf_bdry_synch);         // src_cf_bdry_synch
    if (d_do_log)
    {
        const double Div_U_norm_1  = d_hier_cc_data_ops->L1Norm( d_Div_U_scratch_idx, d_wgt_cc_idx);
        const double Div_U_norm_2  = d_hier_cc_data_ops->L2Norm( d_Div_U_scratch_idx, d_wgt_cc_idx);
        const double Div_U_norm_oo = d_hier_cc_data_ops->maxNorm(d_Div_U_scratch_idx, d_wgt_cc_idx);
        SAMRAI::tbox::plog << d_object_name << "::regridProjection():\n"
                           << "  performing regrid projection\n"
                           << "  before projection:\n"
                           << "    ||Div U||_1  = " << Div_U_norm_1  << "\n"
                           << "    ||Div U||_2  = " << Div_U_norm_2  << "\n"
                           << "    ||Div U||_oo = " << Div_U_norm_oo << "\n";
    }

    // Setup the solver vectors.
    d_hier_cc_data_ops->setToScalar(d_Phi_idx, 0.0, false);
    d_hier_cc_data_ops->scale(d_Div_U_scratch_idx, -1.0, d_Div_U_scratch_idx);
    const double Div_U_mean = (1.0/d_volume)*d_hier_cc_data_ops->integral(d_Div_U_scratch_idx, d_wgt_cc_idx);
    d_hier_cc_data_ops->addScalar(d_Div_U_scratch_idx, d_Div_U_scratch_idx, -Div_U_mean);

    SAMRAI::solv::SAMRAIVectorReal<NDIM,double> sol_vec(d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(d_Phi_var, d_Phi_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAI::solv::SAMRAIVectorReal<NDIM,double> rhs_vec(d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_Div_U_var, d_Div_U_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    // Setup the Poisson solver.
    d_regrid_projection_spec->setCZero();
    d_regrid_projection_spec->setDConstant(-1.0);

    d_regrid_projection_op->setPoissonSpecifications(*d_regrid_projection_spec);
    d_regrid_projection_op->setPhysicalBcCoef(&d_regrid_projection_bc_coef);
    d_regrid_projection_op->setHomogeneousBc(true);
    d_regrid_projection_op->setTime(d_integrator_time);
    d_regrid_projection_op->setHierarchyMathOps(d_hier_math_ops);

    d_regrid_projection_fac_op->setPoissonSpecifications(*d_regrid_projection_spec);
    d_regrid_projection_fac_op->setPhysicalBcCoef(&d_regrid_projection_bc_coef);
    d_regrid_projection_fac_op->setHomogeneousBc(true);
    d_regrid_projection_fac_op->setTime(d_integrator_time);

    d_regrid_projection_solver->setInitialGuessNonzero(false);
    d_regrid_projection_solver->setOperator(d_regrid_projection_op);

    // Solve the projection Poisson problem.
    d_regrid_projection_solver->initializeSolverState(sol_vec,rhs_vec);
    d_regrid_projection_solver->solveSystem(sol_vec,rhs_vec);
    d_regrid_projection_solver->deallocateSolverState();

    // NOTE: We always use homogeneous Neumann boundary conditions for the
    // regrid projection Poisson solver.
    static const bool constant_null_space = true;
    if (constant_null_space)
    {
        const double Phi_mean = (1.0/d_volume)*d_hier_cc_data_ops->integral(d_Phi_idx, d_wgt_cc_idx);
        d_hier_cc_data_ops->addScalar(d_Phi_idx, d_Phi_idx, -Phi_mean);
    }

    // Setup the interpolation transaction information.
    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent Phi_bc_component(d_Phi_idx, CELL_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, &d_regrid_projection_bc_coef);
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> Phi_bdry_bc_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    Phi_bdry_bc_fill_op->initializeOperatorState(Phi_bc_component, d_hierarchy);

    // Fill the physical boundary conditions for Phi.
    Phi_bdry_bc_fill_op->setHomogeneousBc(true);
    Phi_bdry_bc_fill_op->fillData(d_integrator_time);

    // Set U := U - grad Phi.
    const bool U_scratch_cf_bdry_synch = true;
    d_hier_math_ops->grad(
        d_U_scratch_idx, d_U_var,  // dst
        U_scratch_cf_bdry_synch,   // dst_cf_bdry_synch
        1.0,                       // alpha
        d_Phi_idx, d_Phi_var,      // src
        d_no_fill_op,              // src_bdry_fill
        d_integrator_time);        // src_bdry_fill_time
    d_hier_sc_data_ops->axpy(d_U_current_idx, -1.0, d_U_scratch_idx, d_U_current_idx);

    // Compute div U after applying the projection operator
    if (d_do_log)
    {
        const bool U_current_cf_bdry_synch = true;
        d_hier_math_ops->div(
            d_Div_U_scratch_idx, d_Div_U_var, // dst
            +1.0,                             // alpha
            d_U_current_idx, d_U_var,         // src
            d_no_fill_op,                     // src_bdry_fill
            d_integrator_time,                // src_bdry_fill_time
            U_current_cf_bdry_synch);         // src_cf_bdry_synch
        const double Div_U_norm_1  = d_hier_cc_data_ops->L1Norm( d_Div_U_scratch_idx, d_wgt_cc_idx);
        const double Div_U_norm_2  = d_hier_cc_data_ops->L2Norm( d_Div_U_scratch_idx, d_wgt_cc_idx);
        const double Div_U_norm_oo = d_hier_cc_data_ops->maxNorm(d_Div_U_scratch_idx, d_wgt_cc_idx);
        SAMRAI::tbox::plog << "  after projection:\n"
                           << "    ||Div U||_1  = " << Div_U_norm_1  << "\n"
                           << "    ||Div U||_2  = " << Div_U_norm_2  << "\n"
                           << "    ||Div U||_oo = " << Div_U_norm_oo << "\n";
    }

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }
    return;
}// regridProjection

void
INSStaggeredHierarchyIntegrator::initializeOperatorsAndSolvers(
    const double current_time,
    const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::U_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
    U_scratch_vec->addComponent(d_U_var,d_U_scratch_idx,d_wgt_sc_idx,d_hier_sc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_rhs_vec = U_scratch_vec->cloneVector(d_object_name+"::U_rhs_vec");
    const int U_rhs_idx = U_rhs_vec->getComponentDescriptorIndex(0);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > P_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::P_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
    P_scratch_vec->addComponent(d_P_var,d_P_scratch_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > P_rhs_vec = P_scratch_vec->cloneVector(d_object_name+"::P_rhs_vec");
    const int P_rhs_idx = P_rhs_vec->getComponentDescriptorIndex(0);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > sol_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec->addComponent(d_U_var,d_U_scratch_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    sol_vec->addComponent(d_P_var,d_P_scratch_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > rhs_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec->addComponent(d_U_var,U_rhs_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    rhs_vec->addComponent(d_P_var,P_rhs_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    for (int d = 0; d < NDIM; ++d)
    {
        INSStaggeredVelocityBcCoef* U_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(d_U_bc_coefs[d]);
        U_bc_coef->setTimeInterval(current_time,new_time);
    }
    INSStaggeredPressureBcCoef* P_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(d_P_bc_coef);
    P_bc_coef->setTimeInterval(current_time,new_time);
    P_bc_coef->setVelocityCurrentPatchDataIndex(d_U_current_idx);
    P_bc_coef->setVelocityNewPatchDataIndex(d_U_new_idx);

    if (!d_helmholtz_solver.isNull())
    {
        d_helmholtz_spec->setCConstant((d_rho/dt)+0.5*d_lambda);
        d_helmholtz_spec->setDConstant(          -0.5*d_mu    );

        d_helmholtz_op->setPoissonSpecifications(*d_helmholtz_spec);
        d_helmholtz_op->setPhysicalBcCoefs(d_U_star_bc_coefs);
        d_helmholtz_op->setHomogeneousBc(true);
        d_helmholtz_op->setTime(new_time);
        d_helmholtz_op->setHierarchyMathOps(d_hier_math_ops);

        if (!d_helmholtz_hypre_pc.isNull())
        {
            d_helmholtz_hypre_pc->setPoissonSpecifications(*d_helmholtz_spec);
            d_helmholtz_hypre_pc->setPhysicalBcCoefs(d_U_star_bc_coefs);
            d_helmholtz_hypre_pc->setHomogeneousBc(true);
            d_helmholtz_hypre_pc->setTime(new_time);
        }
        else if (!d_helmholtz_fac_op.isNull())
        {
            d_helmholtz_fac_op->setPoissonSpecifications(*d_helmholtz_spec);
            d_helmholtz_fac_op->setPhysicalBcCoefs(d_U_star_bc_coefs);
            d_helmholtz_fac_op->setHomogeneousBc(true);
            d_helmholtz_fac_op->setTime(new_time);
        }

        d_helmholtz_solver->setInitialGuessNonzero(false);
        d_helmholtz_solver->setOperator(d_helmholtz_op);
        if (d_helmholtz_solver_needs_init || !SAMRAI::tbox::MathUtilities<double>::equalEps(dt,d_op_and_solver_init_dt))
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::integrateHierarchy(): Initializing Helmholtz solver" << std::endl;
            d_helmholtz_solver->initializeSolverState(*U_scratch_vec,*U_rhs_vec);
        }
        d_helmholtz_solver_needs_init = false;
    }

    if (!d_poisson_solver.isNull())
    {
        d_poisson_spec->setCZero();
        d_poisson_spec->setDConstant(-1.0);

        d_poisson_op->setPoissonSpecifications(*d_poisson_spec);
        d_poisson_op->setPhysicalBcCoef(d_Phi_bc_coef);
        d_poisson_op->setHomogeneousBc(true);
        d_poisson_op->setTime(current_time+0.5*dt);
        d_poisson_op->setHierarchyMathOps(d_hier_math_ops);

        if (!d_poisson_hypre_pc.isNull())
        {
            d_poisson_hypre_pc->setPoissonSpecifications(*d_poisson_spec);
            d_poisson_hypre_pc->setPhysicalBcCoef(d_Phi_bc_coef);
            d_poisson_hypre_pc->setHomogeneousBc(true);
            d_poisson_hypre_pc->setTime(current_time+0.5*dt);
        }
        else if (!d_poisson_fac_op.isNull())
        {
            d_poisson_fac_op->setPoissonSpecifications(*d_poisson_spec);
            d_poisson_fac_op->setPhysicalBcCoef(d_Phi_bc_coef);
            d_poisson_fac_op->setHomogeneousBc(true);
            d_poisson_fac_op->setTime(current_time+0.5*dt);
        }

        d_poisson_solver->setInitialGuessNonzero(false);
        d_poisson_solver->setOperator(d_poisson_op);
        if (d_poisson_solver_needs_init)
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::integrateHierarchy(): Initializing Poisson solver" << std::endl;
            d_poisson_solver->initializeSolverState(*P_scratch_vec,*P_rhs_vec);
        }
        d_poisson_solver_needs_init = false;
    }

    if (!d_projection_pc.isNull())
    {
        d_projection_pc->setTimeInterval(current_time,new_time,dt);
        if (d_projection_pc_needs_init && !d_stokes_solver_needs_init)
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::integrateHierarchy(): Initializing projection preconditioner" << std::endl;
            d_projection_pc->initializeSolverState(*sol_vec,*rhs_vec);
        }
        d_projection_pc_needs_init = false;
    }

    if (!d_block_pc.isNull())
    {
        d_block_pc->setTimeInterval(current_time,new_time);
        if (d_block_pc_needs_init && !d_stokes_solver_needs_init)
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::integrateHierarchy(): Initializing block-factorization preconditioner" << std::endl;
            d_block_pc->initializeSolverState(*sol_vec,*rhs_vec);
        }
        d_block_pc_needs_init = false;
    }

    if (!d_stokes_solver.isNull())
    {
        d_stokes_op->setTimeInterval(current_time,new_time);
        d_stokes_solver->setOperator(d_stokes_op);
        if (d_stokes_solver_needs_init)
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::integrateHierarchy(): Initializing incompressible Stokes solver" << std::endl;
            d_stokes_solver->initializeSolverState(*sol_vec,*rhs_vec);
        }
        d_stokes_solver_needs_init = false;
    }

    if (!d_convective_op.isNull())
    {
        if (d_convective_op_needs_init)
        {
            d_convective_op->initializeOperatorState(*U_scratch_vec,*U_rhs_vec);
        }
        d_convective_op_needs_init = false;
    }

    U_rhs_vec->freeVectorComponents();
    P_rhs_vec->freeVectorComponents();

    // Keep track of the timestep size to avoid unnecessary re-initialization.
    d_op_and_solver_init_dt = dt;
    return;
}// initializeOperatorsAndSolvers

double
INSStaggeredHierarchyIntegrator::getLevelDt(
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
INSStaggeredHierarchyIntegrator::getPatchDt(
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
    NAVIER_STOKES_SC_STABLEDT_FC(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),
        U_ghost_cells(0),U_ghost_cells(1),
        U_data->getPointer(0),U_data->getPointer(1),
        stable_dt);
#endif
#if (NDIM == 3)
    NAVIER_STOKES_SC_STABLEDT_FC(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
        U_ghost_cells(0),U_ghost_cells(1),U_ghost_cells(2),
        U_data->getPointer(0),U_data->getPointer(1),U_data->getPointer(2),
        stable_dt);
#endif
    return stable_dt;
}// getPatchDt

void
INSStaggeredHierarchyIntegrator::getFromInput(
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

    d_num_cycles = db->getIntegerWithDefault(
        "num_cycles", d_num_cycles);

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

    d_output_U     = db->getBoolWithDefault("output_U"    , d_output_U     );
    d_output_P     = db->getBoolWithDefault("output_P"    , d_output_P     );
    d_output_F     = db->getBoolWithDefault("output_F"    , d_output_F     );
    d_output_Omega = db->getBoolWithDefault("output_Omega", d_output_Omega );
    d_output_Div_U = db->getBoolWithDefault("output_Div_U", d_output_Div_U );

    d_U_scale = db->getDoubleWithDefault("U_scale", d_U_scale);
    d_P_scale = db->getDoubleWithDefault("P_scale", d_P_scale);
    d_F_scale = db->getDoubleWithDefault("F_scale", d_F_scale);

    d_cfl = db->getDoubleWithDefault("cfl",d_cfl);

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault(
        "dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault(
        "dt_max_time_min", d_dt_max_time_min);

    d_helmholtz_hypre_pc_db       = db->isDatabase("HelmholtzHypreSolver") ? db->getDatabase("HelmholtzHypreSolver") : SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL);
    d_helmholtz_fac_pc_db         = db->isDatabase("HelmholtzFACSolver"  ) ? db->getDatabase("HelmholtzFACSolver"  ) : SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL);
    d_poisson_hypre_pc_db         = db->isDatabase("PoissonHypreSolver"  ) ? db->getDatabase("PoissonHypreSolver"  ) : SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL);
    d_poisson_fac_pc_db           = db->isDatabase("PoissonFACSolver"    ) ? db->getDatabase("PoissonFACSolver"    ) : SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL);
    d_regrid_projection_fac_pc_db = db->isDatabase("PoissonFACSolver"    ) ? db->getDatabase("PoissonFACSolver"    ) : SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL);

    d_regrid_max_div_growth_factor = db->getDoubleWithDefault("regrid_max_div_growth_factor", d_regrid_max_div_growth_factor);

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
    }
    return;
}// getFromInput

void
INSStaggeredHierarchyIntegrator::getFromRestart()
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

    int ver = db->getInteger("INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    d_U_scale = db->getDouble("d_U_scale");
    d_P_scale = db->getDouble("d_P_scale");
    d_F_scale = db->getDouble("d_F_scale");
    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    d_grow_dt = db->getDouble("d_grow_dt");
    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");
    d_num_cycles = db->getInteger("d_num_cycles");
    d_regrid_interval = db->getInteger("d_regrid_interval");
    d_using_default_tag_buffer = db->getBool("d_using_default_tag_buffer");
    d_tag_buffer = db->getIntegerArray("d_tag_buffer");
    d_conservation_form = db->getBool("d_conservation_form");
    d_using_vorticity_tagging = db->getBool("d_using_vorticity_tagging");
    d_Omega_rel_thresh = db->getDoubleArray("d_Omega_rel_thresh");
    d_Omega_abs_thresh = db->getDoubleArray("d_Omega_abs_thresh");
    d_Omega_max = db->getDouble("d_Omega_max");
    d_normalize_pressure = db->getBool("d_normalize_pressure");
    d_output_U = db->getBool("d_output_U");
    d_output_P = db->getBool("d_output_P");
    d_output_F = db->getBool("d_output_F");
    d_output_Omega = db->getBool("d_output_Omega");
    d_output_Div_U = db->getBool("d_output_Div_U");
    d_old_dt = db->getDouble("d_old_dt");
    d_op_and_solver_init_dt = db->getDouble("d_op_and_solver_init_dt");
    d_integrator_time = db->getDouble("d_integrator_time");
    d_integrator_step = db->getInteger("d_integrator_step");
    d_cfl = db->getDouble("d_cfl");
    d_dt_max = db->getDouble("d_dt_max");
    d_dt_max_time_max = db->getDouble("d_dt_max_time_max");
    d_dt_max_time_min = db->getDouble("d_dt_max_time_min");
    d_do_log = db->getBool("d_do_log");
    d_rho = db->getDouble("d_rho");
    d_mu = db->getDouble("d_mu");
    d_lambda = db->getDouble("d_lambda");
    d_div_u_abstol = db->getDouble("d_div_u_abstol");
    d_regrid_max_div_growth_factor = db->getDouble("d_regrid_max_div_growth_factor");
    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
