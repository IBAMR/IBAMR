// Filename: IBImplicitHierarchyIntegrator.C
// Created on 08 May 2008 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

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
#include <ibamr/IBAnchorPointSpec.h>
#include <ibamr/INSStaggeredIntermediateVelocityBcCoef.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/INSStaggeredProjectionBcCoef.h>
#include <ibamr/INSStaggeredVelocityBcCoef.h>
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CartSideDoubleDivPreservingRefine.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LNodeIndexData.h>
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
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
#define NAVIER_STOKES_SC_STABLEDT_FC FC_FUNC_(navier_stokes_sc_stabledt2d, NAVIER_STOKES_SC_STABLEDT2D)
#define NAVIER_STOKES_SIDE_TO_fACE_FC FC_FUNC_(navier_stokes_side_to_face2d, NAVIER_STOKES_SIDE_TO_fACE2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_SC_STABLEDT_FC FC_FUNC_(navier_stokes_sc_stabledt3d, NAVIER_STOKES_SC_STABLEDT3D)
#define NAVIER_STOKES_SIDE_TO_fACE_FC FC_FUNC_(navier_stokes_side_to_face3d, NAVIER_STOKES_SIDE_TO_fACE3D)
#endif

extern "C"
{
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

    void
    NAVIER_STOKES_SIDE_TO_fACE_FC(
#if (NDIM == 2)
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* , const int& ,
        double* , double* , const int&
#endif
#if (NDIM == 3)
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* , const int& ,
        double* , double* , double* , const int&
#endif
                                  );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Pointer<Timer> t_initialize_hierarchy_integrator;
static Pointer<Timer> t_initialize_hierarchy;
static Pointer<Timer> t_advance_hierarchy;
static Pointer<Timer> t_get_stable_timestep;
static Pointer<Timer> t_regrid_hierarchy;
static Pointer<Timer> t_integrate_hierarchy_initialize;
static Pointer<Timer> t_integrate_hierarchy;
static Pointer<Timer> t_integrate_hierarchy_finalize;
static Pointer<Timer> t_synchronize_hierarchy;
static Pointer<Timer> t_synchronize_new_levels;
static Pointer<Timer> t_reset_time_dependent_data;
static Pointer<Timer> t_reset_data_to_preadvance_state;
static Pointer<Timer> t_initialize_level_data;
static Pointer<Timer> t_reset_hierarchy_configuration;
static Pointer<Timer> t_apply_gradient_detector;
static Pointer<Timer> t_put_to_database;

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

// Version of IBImplicitHierarchyIntegrator restart file data.
static const int IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitHierarchyIntegrator::IBImplicitHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
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
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Set some default values.
    d_u_scale = 1.0;
    d_p_scale = 1.0;
    d_f_scale = 1.0;
    d_q_scale = 1.0;

    d_start_time = 0.0;
    d_end_time = std::numeric_limits<double>::max();
    d_grow_dt = 2.0;
    d_max_integrator_steps = std::numeric_limits<int>::max();

    d_num_cycles = 3;

    d_using_vorticity_tagging = false;
    d_omega_max = 0.0;

    d_normalize_pressure = false;
    d_convective_difference_form = ADVECTIVE;
    d_creeping_flow = false;

    d_regrid_interval = 1;
    d_regrid_mode = STANDARD;
    d_old_dt = -1.0;
    d_op_and_solver_init_dt = -1.0;
    d_integrator_time = std::numeric_limits<double>::quiet_NaN();
    d_integrator_step = std::numeric_limits<int>::max();

    d_output_u = false;
    d_output_p = false;
    d_output_f = false;
    d_output_q = false;
    d_output_omega = false;
    d_output_div_u = false;

    d_rho    = std::numeric_limits<double>::quiet_NaN();
    d_mu     = std::numeric_limits<double>::quiet_NaN();
    d_lambda = std::numeric_limits<double>::quiet_NaN();

    d_cfl = 0.9;

    d_dt_max = std::numeric_limits<double>::max();
    d_dt_max_time_max = std::numeric_limits<double>::max();
    d_dt_max_time_min = -(d_dt_max_time_max-std::numeric_limits<double>::epsilon());
    d_dt_init = std::numeric_limits<double>::max();

    d_is_initialized = false;

    d_do_log = false;

    d_regrid_max_div_growth_factor = 1.1;

    d_delta_fcn = "IB_4";
    d_ghosts = -1;

    // Setup default boundary condition objects that specify homogeneous
    // Dirichlet boundary conditions for the velocity.
    d_default_u_bc_coef = new LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_u_bc_coef", Pointer<Database>(NULL));
    for (int d = 0; d < NDIM; ++d)
    {
        d_default_u_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_u_bc_coef->setBoundaryValue(2*d+1,0.0);
    }
    registerVelocityPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM,d_default_u_bc_coef));

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Set the problem coefs.
    d_problem_coefs = new INSCoefs(d_rho, d_mu, d_lambda);

    // Get the Lagrangian Data Manager.
    d_lag_data_manager = LDataManager::getManager(d_object_name+"::LDataManager",
                                                  d_delta_fcn, d_delta_fcn,
                                                  d_ghosts,
                                                  d_registered_for_restart);
    d_ghosts = d_lag_data_manager->getGhostCellWidth();

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();

    Pointer<CellVariable<NDIM,double> > cc_var = new CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy, true);

    Pointer<SideVariable<NDIM,double> > sc_var = new SideVariable<NDIM,double>("sc_var");
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, hierarchy, true);

    // Initialize all variable contexts.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_current_context = var_db->getContext(d_object_name+"::CURRENT");
    d_new_context     = var_db->getContext(d_object_name+"::NEW"    );
    d_scratch_context = var_db->getContext(d_object_name+"::SCRATCH");

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_initialize_hierarchy_integrator = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy            = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy               = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::advanceHierarchy()");
        t_get_stable_timestep             = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::getStableTimestep()");
        t_regrid_hierarchy                = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::regridHierarchy()");
        t_integrate_hierarchy_initialize  = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::integrateHierarchy_initialize()");
        t_integrate_hierarchy             = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::integrateHierarchy()");
        t_integrate_hierarchy_finalize    = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::integrateHierarchy_finalize()");
        t_synchronize_hierarchy           = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels          = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data       = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state  = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data           = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration   = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector         = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database                 = TimerManager::getManager()->getTimer("IBAMR::IBImplicitHierarchyIntegrator::putToDatabase()");
                  );
    return;
}// IBImplicitHierarchyIntegrator

IBImplicitHierarchyIntegrator::~IBImplicitHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
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
    delete d_problem_coefs;
    if (d_regrid_projection_spec != NULL) delete d_regrid_projection_spec;
    return;
}// ~IBImplicitHierarchyIntegrator

const std::string&
IBImplicitHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
IBImplicitHierarchyIntegrator::registerVelocityInitialConditions(
    Pointer<CartGridFunction> u_init)
{
    d_u_init = u_init;
    return;
}// registerVelocityInitialConditions

void
IBImplicitHierarchyIntegrator::registerVelocityPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
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
    return;
}// registerVelocityPhysicalBcCoefs

void
IBImplicitHierarchyIntegrator::registerPressureInitialConditions(
    Pointer<CartGridFunction> p_init)
{
    d_p_init = p_init;
    return;
}// registerPressureInitialConditions

void
IBImplicitHierarchyIntegrator::registerBodyForceSpecification(
    Pointer<CartGridFunction> f_fcn)
{
    d_f_fcn = f_fcn;
    return;
}// registerBodyForceSpecification

void
IBImplicitHierarchyIntegrator::registerSourceSpecification(
    Pointer<CartGridFunction> q_fcn)
{
    d_q_fcn = q_fcn;
    pout << "\n"
         << "WARNING: There is an extra term which should be added to the momentum equation\n"
         << "         in the case that div u != 0.\n"
         << "         At the present time, this term has NOT been incorporated into the\n"
         << "         staggered-grid incompressible Navier-Stokes solver.\n"
         << "\n";
    return;
}// registerSourceSpecification

void
IBImplicitHierarchyIntegrator::registerAdvDiffHierarchyIntegrator(
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator)
{
    d_adv_diff_hier_integrator = adv_diff_hier_integrator;
    return;
}// registerAdvDiffHierarchyIntegrator

void
IBImplicitHierarchyIntegrator::registerLNodeInitStrategy(
    Pointer<LNodeInitStrategy> lag_init_strategy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!lag_init_strategy.isNull());
#endif
    d_lag_init_strategy = lag_init_strategy;
    d_lag_data_manager->registerLNodeInitStrategy(d_lag_init_strategy);
    return;
}// registerLNodeInitStrategy

void
IBImplicitHierarchyIntegrator::freeLNodeInitStrategy()
{
    d_lag_init_strategy.setNull();
    d_lag_data_manager->freeLNodeInitStrategy();
    return;
}// freeLNodeInitStrategy

void
IBImplicitHierarchyIntegrator::registerIBLagrangianForceStrategy(
    Pointer<IBLagrangianForceStrategy> lag_force_strategy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!lag_force_strategy.isNull());
#endif
    d_lag_force_strategy = lag_force_strategy;
    return;
}// registerIBLagrangianForceStrategy

void
IBImplicitHierarchyIntegrator::registerVisItDataWriter(
    Pointer<VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    if (!d_adv_diff_hier_integrator.isNull())
    {
        d_adv_diff_hier_integrator->registerVisItDataWriter(visit_writer);
    }
    return;
}// registerVisItDataWriter

void
IBImplicitHierarchyIntegrator::registerLagSiloDataWriter(
    Pointer<LagSiloDataWriter> silo_writer)
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
IBImplicitHierarchyIntegrator::registerLagM3DDataWriter(
    Pointer<LagM3DDataWriter> m3D_writer)
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
IBImplicitHierarchyIntegrator::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > load_balancer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!load_balancer.isNull());
#endif
    d_load_balancer = load_balancer;
    d_lag_data_manager->registerLoadBalancer(d_load_balancer);
    return;
}// registerLoadBalancer

void
IBImplicitHierarchyIntegrator::registerRegridHierarchyCallback(
    void (*callback)(const Pointer<BasePatchHierarchy<NDIM> > hierarchy, const double regrid_data_time, const bool initial_time, void* ctx),
    void* ctx)
{
    d_regrid_hierarchy_callbacks.push_back(callback);
    d_regrid_hierarchy_callback_ctxs.push_back(ctx);
    return;
}// registerRegridHierarchyCallback

void
IBImplicitHierarchyIntegrator::registerApplyGradientDetectorCallback(
    void (*callback)(const Pointer<BasePatchHierarchy<NDIM> > hierarchy, const int level_number, const double error_data_time, const int tag_index, const bool initial_time, const bool uses_richardson_extrapolation_too, void* ctx),
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
///  allow for the sharing of a single HierarchyMathOps object between multiple
///  HierarchyIntegrator objects.
///

Pointer<HierarchyMathOps>
IBImplicitHierarchyIntegrator::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hier_math_ops.isNull());
#endif
    return d_hier_math_ops;
}// getHierarchyMathOps

void
IBImplicitHierarchyIntegrator::setHierarchyMathOps(
    Pointer<HierarchyMathOps> hier_math_ops,
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
IBImplicitHierarchyIntegrator::isManagingHierarchyMathOps() const
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
///  allow the IBImplicitHierarchyIntegrator to be used as a hierarchy integrator.
///

void
IBImplicitHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
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

    // Initialize all variables.
    d_u_var          = new SideVariable<NDIM,double>(d_object_name+"::u"          );
    d_u_cc_var       = new CellVariable<NDIM,double>(d_object_name+"::u_cc",  NDIM);
    d_p_var          = new CellVariable<NDIM,double>(d_object_name+"::p"          );
    d_p_extrap_var   = new CellVariable<NDIM,double>(d_object_name+"::p_extrap"   );
    if (!d_f_fcn.isNull())
    {
        d_f_var      = new SideVariable<NDIM,double>(d_object_name+"::f"          );
        d_f_cc_var   = new CellVariable<NDIM,double>(d_object_name+"::f_cc",  NDIM);
    }
    if (!d_q_fcn.isNull())
    {
        d_q_var      = new CellVariable<NDIM,double>(d_object_name+"::q"          );
    }
#if ( NDIM == 2)
    d_omega_var      = new CellVariable<NDIM,double>(d_object_name+"::omega"      );
#endif
#if ( NDIM == 3)
    d_omega_var      = new CellVariable<NDIM,double>(d_object_name+"::omega", NDIM);
    d_omega_norm_var = new CellVariable<NDIM,double>(d_object_name+"::||omega||_2");
#endif
    d_div_u_var      = new CellVariable<NDIM,double>(d_object_name+"::div_u"      );
    d_phi_var        = new CellVariable<NDIM,double>(d_object_name+"::phi"        );
    d_ib_dof_var     = new SideVariable<NDIM,int   >(d_object_name+"::ib_dof"     );
    d_u_half_ib_var  = new SideVariable<NDIM,double>(d_object_name+"::u_half_ib"  );
    d_u_regrid_var   = new SideVariable<NDIM,double>(d_object_name+"::u_regrid"   );
    d_u_src_var      = new SideVariable<NDIM,double>(d_object_name+"::u_src"      );
    d_indicator_var  = new SideVariable<NDIM,double>(d_object_name+"::indicator"  );

    // Create the default communication algorithms.
    d_fill_after_regrid = new RefineAlgorithm<NDIM>();

    d_calgs["SYNCH_CURRENT_STATE_DATA"] = new CoarsenAlgorithm<NDIM>();
    d_calgs["SYNCH_NEW_STATE_DATA"] = new CoarsenAlgorithm<NDIM>();

    // Register state variables that are maintained by the
    // IBImplicitHierarchyIntegrator.

    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> no_ghosts = 0;

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

    registerVariable(d_p_extrap_current_idx, d_p_extrap_new_idx, d_p_extrap_scratch_idx,
                     d_p_extrap_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");

    if (!d_f_fcn.isNull())
    {
        registerVariable(d_f_current_idx, d_f_new_idx, d_f_scratch_idx,
                         d_f_var, side_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

        registerVariable(d_f_cc_current_idx, d_f_cc_new_idx, d_f_cc_scratch_idx,
                         d_f_cc_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
    }
    else
    {
        d_f_current_idx = -1;
        d_f_new_idx     = -1;
        d_f_scratch_idx = -1;
        d_f_cc_current_idx = -1;
        d_f_cc_new_idx     = -1;
        d_f_cc_scratch_idx = -1;
    }

    if (!d_q_fcn.isNull())
    {
        registerVariable(d_q_current_idx, d_q_new_idx, d_q_scratch_idx,
                         d_q_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSTANT_REFINE");
    }
    else
    {
        d_q_current_idx = -1;
        d_q_new_idx     = -1;
        d_q_scratch_idx = -1;
    }

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

    // Register scratch variables that are maintained by the
    // IBImplicitHierarchyIntegrator.

    registerVariable(      d_phi_idx,       d_phi_var, cell_ghosts);
    registerVariable(   d_ib_dof_idx,    d_ib_dof_var,    d_ghosts);
    registerVariable(d_u_half_ib_idx, d_u_half_ib_var, side_ghosts);  // XXXX
    registerVariable( d_u_regrid_idx,  d_u_regrid_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);
    registerVariable(    d_u_src_idx,     d_u_src_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);
    registerVariable(d_indicator_idx, d_indicator_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);

    d_scratch_data.clrFlag(d_ib_dof_idx);  // handle data allocation manually

    // Register variables for plotting.
    if (!d_visit_writer.isNull())
    {
        if (d_output_u)
        {
            d_visit_writer->registerPlotQuantity(d_u_var->getName(), "VECTOR", d_u_cc_current_idx, 0, d_u_scale);
            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(d_u_var->getName()+stream.str(), "SCALAR", d_u_cc_current_idx, d, d_u_scale);
            }
        }

        if (d_output_p)
        {
            d_visit_writer->registerPlotQuantity(d_p_var->getName(), "SCALAR", d_p_current_idx, 0, d_p_scale);
        }

        if (!d_f_fcn.isNull() && d_output_f)
        {
            d_visit_writer->registerPlotQuantity(d_f_var->getName(), "VECTOR", d_f_cc_current_idx, 0, d_f_scale);
            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(d_f_var->getName()+stream.str(), "SCALAR", d_f_cc_current_idx, d, d_f_scale);
            }
        }

        if (!d_q_fcn.isNull() && d_output_q)
        {
            d_visit_writer->registerPlotQuantity(d_q_var->getName(), "SCALAR", d_q_current_idx, 0, d_q_scale);
        }

        if (d_output_omega)
        {
            d_visit_writer->registerPlotQuantity(d_omega_var->getName(), (NDIM == 2) ? "SCALAR" : "VECTOR", d_omega_current_idx);
#if (NDIM == 3)
            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(d_omega_var->getName()+stream.str(), "SCALAR", d_omega_current_idx, d);
            }
#endif
        }

        if (d_output_div_u)
        {
            d_visit_writer->registerPlotQuantity(d_div_u_var->getName(), "SCALAR", d_div_u_current_idx);
        }
    }

    // Setup (optional) advection-diffusion variables.
    if (!d_adv_diff_hier_integrator.isNull())
    {
        // NOTE: Memory management for the face-centered advection velocity is
        // handled by the advection-diffusion hierarchy integrator class.
        d_u_fc_var = new FaceVariable<NDIM,double>(d_object_name+"::u_fc");
        d_adv_diff_hier_integrator->registerAdvectionVelocity(d_u_fc_var, d_q_fcn.isNull());
        d_adv_diff_hier_integrator->initializeHierarchyIntegrator(d_gridding_alg);
    }

    // Set the current integration time.
    if (!RestartManager::getManager()->isFromRestart())
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }

    // Setup the Hierarchy math operations object.
    d_hier_math_ops = new HierarchyMathOps(d_object_name+"::HierarchyMathOps", d_hierarchy);
    d_is_managing_hier_math_ops = true;

    // Setup physical boundary conditions helper.
    d_u_bc_helper = new INSStaggeredPhysicalBoundaryHelper();

    // Setup the Stokes operator.
    d_stokes_op_needs_init = true;
    d_stokes_op = new INSStaggeredStokesOperator(
        *d_problem_coefs,
        d_u_bc_coefs, d_u_bc_helper, d_p_bc_coef,
        d_hier_math_ops);

    // Setup the convective operator.
    d_convective_op_needs_init = true;
    d_convective_op = new INSStaggeredPPMConvectiveOperator(
        d_convective_difference_form);

    // Setup the IB operator.
    d_ib_op_needs_init = true;
    d_ib_SFR_op = new IBImplicitSFROperator(this);
    d_ib_op = new IBImplicitOperator(d_stokes_op, d_ib_SFR_op);

    // Setup the nonlinear solver.
    const std::string ib_prefix = "ib_";
    d_ib_solver_needs_init = true;
    d_ib_solver = new PETScNewtonKrylovSolver(d_object_name+"::ib_solver", ib_prefix);
    d_ib_solver->setOperator(d_ib_op);
    d_ib_SJR_op = new PETScMFFDJacobianOperator(ib_prefix);  // XXXX
    Pointer<PETScMFFDJacobianOperator> ib_SJR_op = d_ib_SJR_op;
    ib_SJR_op->setOperator(d_ib_SFR_op);
//  d_ib_SJR_op = new IBImplicitSJROperator(this);
    d_ib_mod_helmholtz_pc = new IBImplicitModHelmholtzPETScLevelSolver("IBImplicitModHelmholtzPETScLevelSolver", d_helmholtz_petsc_pc_db);
    d_ib_jac_op = new IBImplicitJacobian(d_stokes_op, d_ib_SJR_op, d_ib_mod_helmholtz_pc);
    d_ib_solver->setJacobian(d_ib_jac_op);

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

        d_regrid_projection_spec = new PoissonSpecifications(d_object_name+"::regrid_projection_spec");
        d_regrid_projection_op = new CCLaplaceOperator(d_object_name+"::Regrid Projection Poisson Operator", *d_regrid_projection_spec, &d_regrid_projection_bc_coef, true);
        d_regrid_projection_op->setHierarchyMathOps(d_hier_math_ops);

        d_regrid_projection_solver = new PETScKrylovLinearSolver(d_object_name+"::Regrid Projection Poisson Krylov Solver", regrid_projection_prefix);
        d_regrid_projection_solver->setInitialGuessNonzero(false);
        d_regrid_projection_solver->setOperator(d_regrid_projection_op);

        TBOX_ASSERT(d_gridding_alg->getMaxLevels() > 1);

        if (d_regrid_projection_fac_pc_db.isNull())
        {
            TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                         "  regrid projection poisson fac pc solver database is null." << std::endl);
        }
        d_regrid_projection_fac_op = new CCPoissonFACOperator(d_object_name+"::Regrid Projection Poisson FAC Operator", d_regrid_projection_fac_pc_db);
        d_regrid_projection_fac_op->setPoissonSpecifications(*d_regrid_projection_spec);
        d_regrid_projection_fac_pc = new IBTK::FACPreconditioner(d_object_name+"::Regrid Projection Poisson Preconditioner", *d_regrid_projection_fac_op, d_regrid_projection_fac_pc_db);
        d_regrid_projection_solver->setPreconditioner(d_regrid_projection_fac_pc);

        // Set some default options.
        d_regrid_projection_solver->setKSPType("gmres");
        d_regrid_projection_solver->setAbsoluteTolerance(1.0e-12);
        d_regrid_projection_solver->setRelativeTolerance(1.0e-08);
        d_regrid_projection_solver->setMaxIterations(25);

        // NOTE: We always use homogeneous Neumann boundary conditions for the
        // regrid projection Poisson solver.
        d_regrid_projection_solver->setNullspace(true, NULL);
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
IBImplicitHierarchyIntegrator::initializeHierarchy()
{
    t_initialize_hierarchy->start();

    if (!d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchy():\n" <<
                   "  must call initializeHierarchyIntegrator() prior to calling initializeHierarchy()." << std::endl);
    }

    // Initialize the patch hierarchy.
    const bool initial_time = !RestartManager::getManager()->isFromRestart();

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

    // Reset the Lagrangian data distribution.
    d_lag_data_manager->beginDataRedistribution();
    d_lag_data_manager->endDataRedistribution();

    // Update the workload.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    d_lag_data_manager->updateWorkloadData(coarsest_ln, finest_ln);

    // Indicate that the force strategy needs to be re-initialized.
    d_lag_force_strategy_needs_init = true;

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
IBImplicitHierarchyIntegrator::advanceHierarchy(
    const double dt)
{
    t_advance_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
#endif

    const double current_time = d_integrator_time;
    const double new_time = d_integrator_time+dt;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time,d_start_time);

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
        if (d_do_log) plog << d_object_name << "::advanceHierarchy(): regridding prior to timestep " << d_integrator_step << "\n";
        regridHierarchy();
    }

    // Integrate all time dependent data.
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): integrating time dependent data\n";
    integrateHierarchy_initialize(current_time, new_time);
    for (int cycle = 0; cycle < d_num_cycles; ++cycle)
    {
        integrateHierarchy(current_time, new_time);
    }
    integrateHierarchy_finalize(current_time, new_time);

    // Synchronize all data.
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): synchronizing data\n";
    synchronizeHierarchy();

    // Reset all time dependent data.
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierData(new_time);

    // Determine the next stable timestep from u(n+1).
    const double dt_next = getStableTimestep(getCurrentContext());

    t_advance_hierarchy->stop();
    return dt_next;
}// advanceHierarchy

double
IBImplicitHierarchyIntegrator::getStableTimestep(
    Pointer<VariableContext> ctx)
{
    t_get_stable_timestep->start();

    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    double dt_next = std::numeric_limits<double>::max();

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
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

    if (initial_time)
    {
        dt_next = std::min(dt_next,d_dt_init);
    }
    else
    {
        dt_next = std::min(dt_next,d_grow_dt*d_old_dt);
    }

    t_get_stable_timestep->stop();
    return dt_next;
}// getStableTimestep()

bool
IBImplicitHierarchyIntegrator::atRegridPoint() const
{
    const int level_number = 0;
    return ((d_integrator_step > 0)
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
    return (d_integrator_step < d_max_integrator_steps);
}// stepsRemaining

const Pointer<PatchHierarchy<NDIM> >
IBImplicitHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

Pointer<GriddingAlgorithm<NDIM> >
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
///      regridHierarchy(),
///      synchronizeHierarchy(),
///      synchronizeNewLevels(),
///      resetTimeDependentHierData(),
///      resetHierDataToPreadvanceState()
///
///  allow the IBImplicitHierarchyIntegrator to provide data management for a time
///  integrator which making use of this class.
///

void
IBImplicitHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time,d_start_time);

    const int coarsest_ln = 0;

    // Determine the divergence of the velocity field before regridding.
    d_hier_math_ops->div(d_div_u_current_idx, d_div_u_var, 1.0, d_u_current_idx, d_u_var, d_no_fill_op, d_integrator_time, false);
    const double div_u_norm_1_pre  = d_hier_cc_data_ops->L1Norm( d_div_u_current_idx, d_wgt_cc_idx);
    const double div_u_norm_2_pre  = d_hier_cc_data_ops->L2Norm( d_div_u_current_idx, d_wgt_cc_idx);
    const double div_u_norm_oo_pre = d_hier_cc_data_ops->maxNorm(d_div_u_current_idx, d_wgt_cc_idx);

    // Update the workload pre-regridding.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_lag_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());

    // Before regridding, begin Lagrangian data movement.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): starting Lagrangian data movement.\n";
    d_lag_data_manager->beginDataRedistribution();

    // Regrid the hierarchy.
    switch (d_regrid_mode)
    {
        case STANDARD:
            d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln, d_integrator_time, d_tag_buffer);
            break;
        case AGGRESSIVE:
            for (int k = 0; k < std::max(1,d_hierarchy->getFinestLevelNumber()); ++k)
            {
                d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln, d_integrator_time, d_tag_buffer);
            }
            break;
        default:
            TBOX_ERROR(d_object_name << "::regridHierarchy():\n"
                       << "  unrecognized regrid mode: " << enum_to_string<RegridMode>(d_regrid_mode) << "." << std::endl);
    }

    // Determine the divergence of the velocity field after regridding.
    d_hier_math_ops->div(d_div_u_current_idx, d_div_u_var, 1.0, d_u_current_idx, d_u_var, d_no_fill_op, d_integrator_time, true);
    const double div_u_norm_1_post  = d_hier_cc_data_ops->L1Norm( d_div_u_current_idx, d_wgt_cc_idx);
    const double div_u_norm_2_post  = d_hier_cc_data_ops->L2Norm( d_div_u_current_idx, d_wgt_cc_idx);
    const double div_u_norm_oo_post = d_hier_cc_data_ops->maxNorm(d_div_u_current_idx, d_wgt_cc_idx);

    // Project the interpolated velocity if needed.
    if (d_needs_regrid_projection && (div_u_norm_1_post  > d_regrid_max_div_growth_factor*div_u_norm_1_pre ||
                                      div_u_norm_2_post  > d_regrid_max_div_growth_factor*div_u_norm_2_pre ||
                                      div_u_norm_oo_post > d_regrid_max_div_growth_factor*div_u_norm_oo_pre))
    {
        regridProjection();
    }
    d_needs_regrid_projection = false;

    // Synchronize the state data on the patch hierarchy.
    const int finest_ln_after_regrid = d_hierarchy->getFinestLevelNumber();
    for (int ln = finest_ln_after_regrid; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_CURRENT_STATE_DATA"][ln]->coarsenData();
    }

    // After regridding, finish Lagrangian data movement.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): finishing Lagrangian data movement.\n";
    d_lag_data_manager->endDataRedistribution();

    // Reset the various Lagrangian objects.
    d_X_data     .resize(finest_ln_after_regrid+1);
    d_X_mid_data .resize(finest_ln_after_regrid+1);
    d_X_new_data .resize(finest_ln_after_regrid+1);
    d_X_half_data.resize(finest_ln_after_regrid+1);
    d_U_half_data.resize(finest_ln_after_regrid+1);
    d_F_half_data.resize(finest_ln_after_regrid+1);
    for (int ln = coarsest_ln; ln <= finest_ln_after_regrid; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(ln == finest_ln_after_regrid);
#endif
            d_X_data     [ln] = d_lag_data_manager->getLNodeLevelData(LDataManager::POSN_DATA_NAME,ln);
            d_X_mid_data [ln] = d_lag_data_manager->createLNodeLevelData("X_mid" ,ln,NDIM);
            d_X_new_data [ln] = d_lag_data_manager->createLNodeLevelData("X_new" ,ln,NDIM);
            d_X_half_data[ln] = d_lag_data_manager->createLNodeLevelData("X_half",ln,NDIM);
            d_U_half_data[ln] = d_lag_data_manager->createLNodeLevelData("U_half",ln,NDIM);
            d_F_half_data[ln] = d_lag_data_manager->createLNodeLevelData("F_half",ln,NDIM);
        }
    }

    // Update the workload post-regridding.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_lag_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());

    // Indicate that the force strategy needs to be re-initialized.
    d_lag_force_strategy_needs_init  = true;

    // Compute the set of local anchor points.
    static const double eps = 2.0*sqrt(std::numeric_limits<double>::epsilon());
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const grid_xLower = grid_geom->getXLower();
    const double* const grid_xUpper = grid_geom->getXUpper();
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();
    for (int ln = coarsest_ln; ln <= finest_ln_after_regrid; ++ln)
    {
        d_anchor_point_local_idxs[ln].clear();
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
                const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
                for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
                     it != idx_data->lnode_index_end(); ++it)
                {
                    const LNodeIndex& node_idx = *it;
                    Pointer<IBAnchorPointSpec> anchor_point_spec = node_idx.getNodeData<IBAnchorPointSpec>();
                    if (!anchor_point_spec.isNull())
                    {
                        d_anchor_point_local_idxs[ln].insert(node_idx.getLocalPETScIndex());
                    }
                }
            }

            Pointer<LNodeLevelData> X_data = d_lag_data_manager->getLNodeLevelData(LDataManager::POSN_DATA_NAME,ln);
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

    // Execute any registered callback functions.
    for (size_t i = 0; i < d_regrid_hierarchy_callbacks.size(); ++i)
    {
        (*d_regrid_hierarchy_callbacks[i])(d_hierarchy, d_integrator_time, initial_time, d_regrid_hierarchy_callback_ctxs[i]);
    }

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

void
IBImplicitHierarchyIntegrator::integrateHierarchy_initialize(
    const double current_time,
    const double new_time)
{
    t_integrate_hierarchy_initialize->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    if (d_do_log) plog << d_object_name << "::integrateHierarchy_initialize(): current_time = " << current_time << ", new_time = " << new_time << ", dt = " << dt << "\n";

    // Reset the various Lagrangian objects.
    d_X_data     .resize(finest_ln+1);
    d_X_mid_data .resize(finest_ln+1);
    d_X_new_data .resize(finest_ln+1);
    d_X_half_data.resize(finest_ln+1);
    d_U_half_data.resize(finest_ln+1);
    d_F_half_data.resize(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(ln == finest_ln);
#endif
            d_X_data     [ln] = d_lag_data_manager->getLNodeLevelData(LDataManager::POSN_DATA_NAME,ln);
            d_X_mid_data [ln] = d_lag_data_manager->createLNodeLevelData("X_mid" ,ln,NDIM);
            d_X_new_data [ln] = d_lag_data_manager->createLNodeLevelData("X_new" ,ln,NDIM);
            d_X_half_data[ln] = d_lag_data_manager->createLNodeLevelData("X_half",ln,NDIM);
            d_U_half_data[ln] = d_lag_data_manager->createLNodeLevelData("U_half",ln,NDIM);
            d_F_half_data[ln] = d_lag_data_manager->createLNodeLevelData("F_half",ln,NDIM);
        }
    }

    // Synchronize current state data.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_CURRENT_STATE_DATA"][ln]->coarsenData();
    }

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data    ,     new_time);
    }

    // Setup the solver vectors.
    d_u_scratch_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::u_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
    d_u_scratch_vec->addComponent(d_u_var, d_u_scratch_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    d_u_rhs_vec = d_u_scratch_vec->cloneVector(d_object_name+"::u_rhs_vec");
    d_u_rhs_vec->allocateVectorData(current_time);
    const int u_rhs_idx = d_u_rhs_vec->getComponentDescriptorIndex(0);
    const Pointer<SideVariable<NDIM,double> > u_rhs_var = d_u_rhs_vec->getComponentVariable(0);
    d_hier_sc_data_ops->setToScalar(u_rhs_idx,0.0);

    d_u_half_vec = d_u_scratch_vec->cloneVector(d_object_name+"::u_half_vec");
    d_u_half_vec->allocateVectorData(current_time);
    const int u_half_idx = d_u_half_vec->getComponentDescriptorIndex(0);
    const Pointer<SideVariable<NDIM,double> > u_half_var = d_u_half_vec->getComponentVariable(0);
    d_hier_sc_data_ops->setToScalar(u_half_idx,0.0);

    d_n_vec = d_u_scratch_vec->cloneVector(d_object_name+"::n_vec");
    d_n_vec->allocateVectorData(current_time);
    const int n_idx = d_n_vec->getComponentDescriptorIndex(0);
    const Pointer<SideVariable<NDIM,double> > n_var = d_n_vec->getComponentVariable(0);
    d_hier_sc_data_ops->setToScalar(n_idx,0.0);

    d_p_scratch_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::p_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
    d_p_scratch_vec->addComponent(d_p_var, d_p_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    d_p_rhs_vec = d_p_scratch_vec->cloneVector(d_object_name+"::p_rhs_vec");
    d_p_rhs_vec->allocateVectorData(current_time);
    const int p_rhs_idx = d_p_rhs_vec->getComponentDescriptorIndex(0);
    const Pointer<CellVariable<NDIM,double> > p_rhs_var = d_p_rhs_vec->getComponentVariable(0);
    d_hier_cc_data_ops->setToScalar(p_rhs_idx,0.0);

    // Initialize the right-hand side terms.
    PoissonSpecifications rhs_spec(d_object_name+"::rhs_spec");
    rhs_spec.setCConstant((d_rho/dt)-0.5*d_lambda);
    rhs_spec.setDConstant(          +0.5*d_mu    );
    d_hier_sc_data_ops->copyData(d_u_scratch_idx, d_u_current_idx);
    d_hier_math_ops->laplace(
        u_rhs_idx, u_rhs_var,
        rhs_spec,
        d_u_scratch_idx, d_u_var,
        d_u_bdry_bc_fill_op, current_time);

    // Reset the solution, rhs, and nullspace vectors.
    d_sol_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    d_sol_vec->addComponent(d_u_var,d_u_scratch_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    d_sol_vec->addComponent(d_p_var,d_p_scratch_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    d_rhs_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    d_rhs_vec->addComponent(d_u_var,u_rhs_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    d_rhs_vec->addComponent(d_p_var,p_rhs_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    // Setup the operators and solvers.
    initializeOperatorsAndSolvers(current_time, new_time);

    // Setup the nullspace object.
    if (d_normalize_pressure)
    {
        d_nul_vec = d_sol_vec->cloneVector(d_object_name+"::nul_vec");
        d_nul_vec->allocateVectorData(current_time);
        d_hier_sc_data_ops->setToScalar(d_nul_vec->getComponentDescriptorIndex(0), 0.0);
        d_hier_cc_data_ops->setToScalar(d_nul_vec->getComponentDescriptorIndex(1), 1.0);
        d_ib_solver->getLinearSolver()->setNullspace(false, d_nul_vec);
        if (!d_poisson_solver.isNull()) d_poisson_solver->setNullspace(true, NULL);
    }

    // Set the initial guess.
    d_hier_sc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(0), d_u_current_idx);
    d_hier_cc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(1), d_p_current_idx);
    d_hier_sc_data_ops->copyData(d_u_new_idx, d_sol_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_p_new_idx, d_sol_vec->getComponentDescriptorIndex(1));

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            PetscErrorCode ierr;
            Vec X_vec     = d_X_data    [ln]->getGlobalVec();
            Vec X_new_vec = d_X_new_data[ln]->getGlobalVec();
            ierr = VecCopy(X_vec,X_new_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Setup inhomogeneous boundary conditions.
    d_u_bc_helper->clearBcCoefData();
    d_u_bc_helper->cacheBcCoefData(d_u_scratch_idx, d_u_var, d_u_bc_coefs, new_time, IntVector<NDIM>(SIDEG), d_hierarchy);

    t_integrate_hierarchy_initialize->stop();
    return;
}// integrateHierarchy_initialize

void
IBImplicitHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time)
{
    t_integrate_hierarchy->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double dt = new_time-current_time;

    // Perform a single step of fixed point iteration.

    // Compute u_half := 0.5*(u(n)+u(n+1)).
    const int u_half_idx = d_u_half_vec->getComponentDescriptorIndex(0);
    d_hier_sc_data_ops->linearSum(u_half_idx, 0.5, d_u_current_idx, 0.5, d_u_new_idx);

    // Compute (u_half*grad)u_half.
    const int n_idx = d_n_vec->getComponentDescriptorIndex(0);
    if (d_creeping_flow)
    {
        d_hier_sc_data_ops->setToScalar(n_idx, 0.0);
    }
    else
    {
        d_convective_op->applyConvectiveOperator(u_half_idx, n_idx);
    }

    // Setup the right-hand side vector.
    d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), -d_rho, n_idx, d_rhs_vec->getComponentDescriptorIndex(0));
    if (!d_f_fcn.isNull())
    {
        d_f_fcn->setDataOnPatchHierarchy(d_f_scratch_idx, d_f_var, d_hierarchy, current_time+0.5*dt);
        d_hier_sc_data_ops->add(d_rhs_vec->getComponentDescriptorIndex(0), d_rhs_vec->getComponentDescriptorIndex(0), d_f_scratch_idx);
    }
    if (!d_q_fcn.isNull())
    {
        d_q_fcn->setDataOnPatchHierarchy(d_q_scratch_idx, d_q_var, d_hierarchy, current_time+0.5*dt);
        d_hier_sc_data_ops->subtract(d_rhs_vec->getComponentDescriptorIndex(1), d_rhs_vec->getComponentDescriptorIndex(1), d_q_scratch_idx);
    }

    // Compute X_mid := 0.5*(X(n)+X(n+1)).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            PetscErrorCode ierr;
            Vec X_vec     = d_X_data    [ln]->getGlobalVec();
            Vec X_new_vec = d_X_new_data[ln]->getGlobalVec();
            Vec X_mid_vec = d_X_mid_data[ln]->getGlobalVec();
            ierr = VecCopy(X_vec,X_mid_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(X_mid_vec,0.5,0.5,X_new_vec);  IBTK_CHKERRQ(ierr);
            d_X_mid_data[ln]->beginGhostUpdate();
            d_X_mid_data[ln]->endGhostUpdate();
        }
    }

    // Construct the interpolation and spreading operators.
    d_R_mats.resize(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            int ierr;
            Pointer<PatchLevel<NDIM> > patch_level = d_hierarchy->getPatchLevel(ln);
            Vec X_mid_vec = d_X_mid_data[ln]->getGlobalVec();
            PETScMatUtilities::constructPatchLevelInterpOp(
                d_R_mats[ln], X_mid_vec,
                d_u_half_ib_idx, d_u_half_ib_var, d_ib_dof_idx, d_ib_dof_var,
                patch_level);
            const int global_node_offset = d_lag_data_manager->getGlobalNodeOffset(ln);
            std::vector<int> anchored_rows;
            anchored_rows.reserve(NDIM*d_anchor_point_local_idxs.size());
            for (std::set<int>::const_iterator cit = d_anchor_point_local_idxs[ln].begin();
                 cit != d_anchor_point_local_idxs[ln].end(); ++cit)
            {
                const int local_node_idx = *cit;
                for (int d = 0; d < NDIM; ++d)
                {
                    anchored_rows.push_back(NDIM*(local_node_idx+global_node_offset)+d);
                }
            }
            ierr = MatZeroRows(d_R_mats[ln], anchored_rows.size(), anchored_rows.empty() ? PETSC_NULL : &anchored_rows[0], 0.0); IBTK_CHKERRQ(ierr);
        }
    }

    // Ensure there is no forcing at Dirichlet boundaries (the Dirichlet
    // boundary condition takes precedence).
    d_u_bc_helper->zeroValuesAtDirichletBoundaries(d_rhs_vec->getComponentDescriptorIndex(0));

    // Synchronize solution and right-hand-side data before solve.
    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
    SynchronizationTransactionComponent sol_synch_transaction =
        SynchronizationTransactionComponent(d_sol_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");
    SynchronizationTransactionComponent rhs_synch_transaction =
        SynchronizationTransactionComponent(d_rhs_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");

    d_side_synch_op->resetTransactionComponent(sol_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);
    d_side_synch_op->resetTransactionComponent(rhs_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);

    // Solve for u(n+1), p(n+1/2), and X(n+1).
    d_ib_solver->solveSystem(*d_sol_vec,*d_rhs_vec);

    // Synchronize solution data after solve.
    d_side_synch_op->resetTransactionComponent(sol_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);

    // Update the values of any advected-and-diffused quantities registered with
    // the optional advection-diffusion solver.
    if (!d_adv_diff_hier_integrator.isNull())
    {
        d_adv_diff_hier_integrator->resetHierDataToPreadvanceState();

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int u_adv_idx = var_db->mapVariableAndContextToIndex(d_u_fc_var, d_adv_diff_hier_integrator->getCurrentContext());
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Index<NDIM>& patch_lower = patch_box.lower();
                const Index<NDIM>& patch_upper = patch_box.upper();
                Pointer<SideData<NDIM,double> > u_half_data = patch->getPatchData(u_half_idx);
                const int u_half_data_gcw = u_half_data->getGhostCellWidth().min();
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(u_half_data_gcw == u_half_data->getGhostCellWidth().max());
#endif
                Pointer<FaceData<NDIM,double> > u_adv_data = patch->getPatchData(u_adv_idx);
                const int u_adv_data_gcw = u_adv_data->getGhostCellWidth().min();
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(u_adv_data_gcw == u_adv_data->getGhostCellWidth().max());
#endif
                NAVIER_STOKES_SIDE_TO_fACE_FC(
                    patch_lower(0), patch_upper(0),
                    patch_lower(1), patch_upper(1),
#if (NDIM == 3)
                    patch_lower(2), patch_upper(2),
#endif
                    u_half_data->getPointer(0),
                    u_half_data->getPointer(1),
#if (NDIM == 3)
                    u_half_data->getPointer(2),
#endif
                    u_half_data_gcw,
                    u_adv_data->getPointer(0),
                    u_adv_data->getPointer(1),
#if (NDIM == 3)
                    u_adv_data->getPointer(2),
#endif
                    u_adv_data_gcw);
            }
        }
        d_adv_diff_hier_integrator->integrateHierarchy(current_time, new_time);
        d_adv_diff_hier_integrator->synchronizeHierarchy();
    }

    // Enforce Dirichlet boundary conditions.
    d_u_bc_helper->resetValuesAtDirichletBoundaries(d_sol_vec->getComponentDescriptorIndex(0));

    // Pull out solution components.
    d_hier_sc_data_ops->copyData(d_u_new_idx, d_sol_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_p_new_idx, d_sol_vec->getComponentDescriptorIndex(1));

    // Reset the right-hand side vector.
    d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), +d_rho, n_idx, d_rhs_vec->getComponentDescriptorIndex(0));
    if (!d_f_fcn.isNull())
    {
        d_hier_sc_data_ops->subtract(d_rhs_vec->getComponentDescriptorIndex(0), d_rhs_vec->getComponentDescriptorIndex(0), d_f_scratch_idx);
        d_hier_sc_data_ops->copyData(d_f_new_idx, d_f_scratch_idx);
    }
    if (!d_q_fcn.isNull())
    {
        d_hier_sc_data_ops->add(d_rhs_vec->getComponentDescriptorIndex(1), d_rhs_vec->getComponentDescriptorIndex(1), d_q_scratch_idx);
        d_hier_sc_data_ops->copyData(d_q_new_idx, d_q_scratch_idx);
    }

    // Set u(n+1/2) = 0.5*(u(n) + u(n+1)).
    d_hier_sc_data_ops->linearSum(d_u_half_ib_idx, 0.5, d_u_current_idx, 0.5, d_u_new_idx);

    // Interpolate u(n+1/2) to U(n+1/2) and set X(n+1) = X(n) + dt*U(n+1/2).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            int ierr;
            Pointer<PatchLevel<NDIM> > patch_level = d_hierarchy->getPatchLevel(ln);
            Vec U_half_vec = d_U_half_data[ln]->getGlobalVec();
            Vec u_half_ib_vec = static_cast<Vec>(NULL);
            PETScVecUtilities::constructPatchLevelVec(u_half_ib_vec, d_u_half_ib_idx, d_u_half_ib_var, patch_level);
            PETScVecUtilities::copyToPatchLevelVec(u_half_ib_vec, d_u_half_ib_idx, d_u_half_ib_var, patch_level);
            ierr = MatMult(d_R_mats[ln], u_half_ib_vec, U_half_vec); IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(u_half_ib_vec); IBTK_CHKERRQ(ierr);
            Vec X_vec = d_X_data[ln]->getGlobalVec();
            Vec X_new_vec = d_X_new_data[ln]->getGlobalVec();
            ierr = VecWAXPY(X_new_vec, dt, U_half_vec, X_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Destroy the restriction matrices.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            int ierr;
            ierr = MatDestroy(d_R_mats[ln]);  IBTK_CHKERRQ(ierr);
        }
    }
    d_R_mats.clear();

    t_integrate_hierarchy->stop();
    return;
}// integrateHierarchy

void
IBImplicitHierarchyIntegrator::integrateHierarchy_finalize(
    const double current_time,
    const double new_time)
{
    t_integrate_hierarchy_finalize->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;

    // Synchronize new state data.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_NEW_STATE_DATA"][ln]->coarsenData();
    }

    // Extrapolate the pressure forward in time to obtain an approximation to
    // p^{n+1}.
    if (d_old_dt > 0.0)
    {
        d_hier_cc_data_ops->linearSum(d_p_extrap_new_idx, (2.0*dt+d_old_dt)/(dt+d_old_dt), d_p_new_idx, -dt/(dt+d_old_dt), d_p_current_idx);
    }

    // Compute the cell-centered approximation to u^{n+1} (used for
    // visualization only).
    reinterpolateVelocity(getNewContext());

    // Compute the cell-centered approximation to f^{n+1} (used for
    // visualization only).
    if (!d_f_fcn.isNull()) reinterpolateForce(getNewContext());

    // Compute omega = curl u.
    //
    // NOTE: Re-filling ghost cell data here overrides the conservative
    // coarsening of u from fine levels to coarse levels.  In particular, this
    // means that we need to re-coarsen the data associated with patch data
    // descriptor index d_u_scratch_idx if we wish to compute, e.g., the
    // divergence of the velocity field.  This operation
    d_hier_sc_data_ops->copyData(d_u_scratch_idx, d_u_new_idx);
    d_hier_math_ops->curl(
        d_omega_new_idx, d_omega_var,
        d_u_scratch_idx, d_u_var,
        d_u_bdry_bc_fill_op, new_time);
#if (NDIM == 3)
    d_hier_math_ops->pointwiseL2Norm(
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

    // Compute div u.
    d_hier_math_ops->div(
        d_div_u_new_idx, d_div_u_var,
        1.0, d_u_new_idx, d_u_var,
        d_no_fill_op, new_time, false);

    // Deallocate the nullspace object.
    if (d_normalize_pressure)
    {
        PetscErrorCode ierr;
        MatNullSpace petsc_nullsp;
        Pointer<PETScKrylovLinearSolver> ksp_solver = d_ib_solver->getLinearSolver();
        KSP petsc_ksp = ksp_solver->getPETScKSP();
        ierr = KSPGetNullSpace(petsc_ksp, &petsc_nullsp); IBTK_CHKERRQ(ierr);
        ierr = MatNullSpaceDestroy(petsc_nullsp); IBTK_CHKERRQ(ierr);
    }

    // Deallocate scratch data.
    d_u_rhs_vec->freeVectorComponents();
    d_u_half_vec->freeVectorComponents();
    d_n_vec->freeVectorComponents();
    d_p_rhs_vec->freeVectorComponents();
    if (d_normalize_pressure)
    {
        d_nul_vec->freeVectorComponents();
    }

    // Deallocate solver vectors.
    d_u_scratch_vec.setNull();
    d_u_rhs_vec.setNull();
    d_u_half_vec.setNull();
    d_n_vec.setNull();
    d_p_scratch_vec.setNull();
    d_p_rhs_vec.setNull();
    d_sol_vec.setNull();
    d_rhs_vec.setNull();
    d_nul_vec.setNull();

    t_integrate_hierarchy_finalize->stop();
    return;
}// integrateHierarchy_finalize

void
IBImplicitHierarchyIntegrator::synchronizeHierarchy()
{
    t_synchronize_hierarchy->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_NEW_STATE_DATA"][ln]->coarsenData();
    }

    if (!d_adv_diff_hier_integrator.isNull())
    {
        d_adv_diff_hier_integrator->synchronizeHierarchy();
    }

    t_synchronize_hierarchy->stop();
    return;
}// synchronizeHierarchy

void
IBImplicitHierarchyIntegrator::synchronizeNewLevels(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
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

    if (!d_adv_diff_hier_integrator.isNull())
    {
        d_adv_diff_hier_integrator->synchronizeNewLevels(hierarchy, coarsest_level, finest_level,
                                                         sync_time, initial_time);
    }

    t_synchronize_new_levels->stop();
    return;
}// synchronizeNewLevels

void
IBImplicitHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    t_reset_time_dependent_data->start();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // Swap PatchData<NDIM> pointers between the current and new contexts.
    for (std::list<Pointer<Variable<NDIM> > >::const_iterator sv = d_state_variables.begin();
         sv != d_state_variables.end(); ++sv)
    {
        const Pointer<Variable<NDIM> >& v = *sv;
        const int src_idx = var_db->mapVariableAndContextToIndex(v, getNewContext());
        const int dst_idx = var_db->mapVariableAndContextToIndex(v, getCurrentContext());
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());

                Pointer<PatchData<NDIM> > src_data = patch->getPatchData(src_idx);
                Pointer<PatchData<NDIM> > dst_data = patch->getPatchData(dst_idx);

#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(src_data->getBox() == dst_data->getBox());
                TBOX_ASSERT(src_data->getGhostCellWidth() == dst_data->getGhostCellWidth());
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
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->setTime(new_time, d_current_data);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }

    if (!d_adv_diff_hier_integrator.isNull())
    {
        d_adv_diff_hier_integrator->resetTimeDependentHierData(new_time);
    }

    t_reset_time_dependent_data->stop();
    return;
}// resetTimeDependentHierData

void
IBImplicitHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    t_reset_data_to_preadvance_state->start();

    // Deallocate the scratch and new data and reset the time of the current
    // data.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->setTime(d_integrator_time, d_current_data);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }

    if (!d_adv_diff_hier_integrator.isNull())
    {
        d_adv_diff_hier_integrator->resetHierDataToPreadvanceState();
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
IBImplicitHierarchyIntegrator::initializeLevelData(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > base_old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const Pointer<PatchLevel<NDIM> > old_level = base_old_level;
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
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    if (allocate_data)
    {
        level->allocatePatchData(d_current_data, init_data_time);
        level->allocatePatchData(d_ib_dof_idx, init_data_time);
    }
    else
    {
        level->setTime(init_data_time, d_current_data);
        level->setTime(d_ib_dof_idx, d_current_data);
    }

    // Initialize DOF index data.
    level->allocatePatchData(d_u_half_ib_idx, init_data_time);
    PETScVecUtilities::constructPatchLevelDOFIndices(d_ib_dof_idx, d_ib_dof_var, d_u_half_ib_idx, d_u_half_ib_var, level);
    level->deallocatePatchData(d_u_half_ib_idx);

    // Fill data from coarser levels in AMR hierarchy.
    if (!initial_time && (level_number > 0 || !old_level.isNull()))
    {
        level->allocatePatchData(d_scratch_data, init_data_time);

        CartExtrapPhysBdryOp extrap_bc_op(d_fill_after_regrid_bc_idxs, BDRY_EXTRAP_TYPE);
        CartSideRobinPhysBdryOp phys_bdry_bc_op(d_u_scratch_idx, d_u_bc_coefs, false);
        std::vector<RefinePatchStrategy<NDIM>*> refine_patch_strategies(2);
        refine_patch_strategies[0] = &extrap_bc_op;
        refine_patch_strategies[1] = &phys_bdry_bc_op;
        RefinePatchStrategySet patch_strategy_set(refine_patch_strategies.begin(), refine_patch_strategies.end(), false);
        d_fill_after_regrid->createSchedule(level, old_level, level_number-1, hierarchy, &patch_strategy_set)->fillData(init_data_time);

        if (level_number > 0)
        {
            // Use divergence- and curl-preserving prolongation to re-set the
            // velocity data.
            ComponentSelector div_free_prolongation_scratch_data;
            div_free_prolongation_scratch_data.setFlag( d_u_regrid_idx);
            div_free_prolongation_scratch_data.setFlag(    d_u_src_idx);
            div_free_prolongation_scratch_data.setFlag(d_indicator_idx);

            // Set the indicator data to equal "0" in each patch of the new
            // patch level and initialize values of u to cause a floating point
            // error if used incorrectly.
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());

                Pointer<SideData<NDIM,double> > indicator_data = patch->getPatchData(d_indicator_idx);
                indicator_data->fillAll(0.0);

                Pointer<SideData<NDIM,double> > u_current_data = patch->getPatchData(d_u_current_idx);
                Pointer<SideData<NDIM,double> >  u_regrid_data = patch->getPatchData( d_u_regrid_idx);
                Pointer<SideData<NDIM,double> >     u_src_data = patch->getPatchData(    d_u_src_idx);
                u_current_data->fillAll(std::numeric_limits<double>::quiet_NaN());
                u_regrid_data ->fillAll(std::numeric_limits<double>::quiet_NaN());
                u_src_data    ->fillAll(std::numeric_limits<double>::quiet_NaN());
            }

            if (!old_level.isNull())
            {
                // Set the indicator data to equal "1" on each patch of the old
                // patch level and reset u.
                old_level->allocatePatchData(div_free_prolongation_scratch_data, init_data_time);
                for (PatchLevel<NDIM>::Iterator p(old_level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = old_level->getPatch(p());

                    Pointer<SideData<NDIM,double> > indicator_data = patch->getPatchData(d_indicator_idx);
                    indicator_data->fillAll(1.0);

                    Pointer<SideData<NDIM,double> > u_current_data = patch->getPatchData(d_u_current_idx);
                    Pointer<SideData<NDIM,double> >  u_regrid_data = patch->getPatchData( d_u_regrid_idx);
                    Pointer<SideData<NDIM,double> >     u_src_data = patch->getPatchData(    d_u_src_idx);
                    u_regrid_data->copy(*u_current_data);
                    u_src_data   ->copy(*u_current_data);
                }

                // Create a communications schedule to copy data from the old patch
                // level to the new patch level.
                //
                // Note that this will set the indicator data to equal "1" at
                // each location in the new patch level which is a copy of a
                // location from the old patch level.
                RefineAlgorithm<NDIM> copy_data;
                copy_data.registerRefine( d_u_regrid_idx,  d_u_regrid_idx,  d_u_regrid_idx, NULL);
                copy_data.registerRefine(    d_u_src_idx,     d_u_src_idx,     d_u_src_idx, NULL);
                copy_data.registerRefine(d_indicator_idx, d_indicator_idx, d_indicator_idx, NULL);
                ComponentSelector bc_fill_data;
                bc_fill_data.setFlag(d_u_regrid_idx);
                bc_fill_data.setFlag(   d_u_src_idx);
                CartSideRobinPhysBdryOp phys_bdry_bc_op(bc_fill_data, d_u_bc_coefs, false);
                copy_data.createSchedule(level, old_level, &phys_bdry_bc_op)->fillData(init_data_time);
            }

            // Setup the divergence- and curl-preserving prolongation refine
            // algorithm and refine the velocity data.
            RefineAlgorithm<NDIM> fill_div_free_prolongation;
            fill_div_free_prolongation.registerRefine(d_u_current_idx, d_u_current_idx, d_u_regrid_idx, NULL);
            CartSideRobinPhysBdryOp phys_bdry_bc_op(d_u_regrid_idx, d_u_bc_coefs, false);
            CartSideDoubleDivPreservingRefine div_preserving_op(d_u_regrid_idx, d_u_src_idx, d_indicator_idx, init_data_time, Pointer<RefinePatchStrategy<NDIM> >(&phys_bdry_bc_op,false));
            fill_div_free_prolongation.createSchedule(level, old_level, level_number-1, hierarchy, &div_preserving_op)->fillData(init_data_time);

            // Free scratch data.
            if (!old_level.isNull()) old_level->deallocatePatchData(div_free_prolongation_scratch_data);
        }
        level->deallocatePatchData(d_scratch_data);
    }

    // Initialize level data at the initial time.
    if (initial_time)
    {
        // If no initialization object is provided, initialize the velocity,
        // divergence, and vorticity to zero.  Otherwise, use the initialization
        // object to set the velocity to some specified value and compute the
        // divergence and vorticity corresponding to the initial velocity.
        if (d_u_init.isNull())
        {
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());

                Pointer<SideData<NDIM,double> > u_current_data = patch->getPatchData(d_u_current_idx);
                u_current_data->fillAll(0.0);

                Pointer<CellData<NDIM,double> > u_cc_current_data = patch->getPatchData(d_u_cc_current_idx);
                u_cc_current_data->fillAll(0.0);

                Pointer<CellData<NDIM,double> > omega_current_data = patch->getPatchData(d_omega_current_idx);
                omega_current_data->fillAll(0.0);
#if (NDIM == 3)
                Pointer<CellData<NDIM,double> > omega_norm_current_data = patch->getPatchData(d_omega_norm_current_idx);
                omega_norm_current_data->fillAll(0.0);
#endif
                Pointer<CellData<NDIM,double> > div_u_current_data = patch->getPatchData(d_div_u_current_idx);
                div_u_current_data->fillAll(0.0);
            }
        }
        else
        {
            level->allocatePatchData(d_u_scratch_idx, init_data_time);

            // Initialize u.
            d_u_init->setDataOnPatchLevel(d_u_current_idx, d_u_var, level, init_data_time, initial_time);
            PatchMathOps patch_math_ops;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());

                Pointer<SideData<NDIM,double> > u_current_data = patch->getPatchData(d_u_current_idx);
                Pointer<CellData<NDIM,double> > u_cc_current_data = patch->getPatchData(d_u_cc_current_idx);

                patch_math_ops.interp(u_cc_current_data, u_current_data, patch);
            }

            // Fill in u boundary data from coarser levels.
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
            Pointer<RefineAlgorithm<NDIM> > ralg = new RefineAlgorithm<NDIM>();
            Pointer<RefineOperator<NDIM> > refine_operator = grid_geom->lookupRefineOperator(
                d_u_var, "CONSERVATIVE_LINEAR_REFINE");
            ralg->registerRefine(d_u_scratch_idx, // destination
                                 d_u_current_idx, // source
                                 d_u_scratch_idx, // temporary work space
                                 refine_operator);
            CartExtrapPhysBdryOp bc_op(d_u_scratch_idx, BDRY_EXTRAP_TYPE);
            ralg->createSchedule(level, level_number-1, hierarchy, &bc_op)->fillData(init_data_time);

            // Initialize quantities derived from the initial value of u.
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());

                Pointer<SideData<NDIM,double> > u_scratch_data = patch->getPatchData(d_u_scratch_idx);

                Pointer<CellData<NDIM,double> > omega_current_data = patch->getPatchData(d_omega_current_idx);
                patch_math_ops.curl(omega_current_data, u_scratch_data, patch);
#if (NDIM == 3)
                Pointer<CellData<NDIM,double> > omega_norm_current_data = patch->getPatchData(d_omega_norm_current_idx);
                patch_math_ops.pointwiseL2Norm(omega_norm_current_data, omega_current_data, patch);
#endif
                Pointer<CellData<NDIM,double> > div_u_current_data = patch->getPatchData(d_div_u_current_idx);
                patch_math_ops.div(div_u_current_data, 1.0, u_scratch_data, 0.0, Pointer<CellData<NDIM,double> >(NULL), patch);
            }

            level->deallocatePatchData(d_u_scratch_idx);
        }

        // Initialize the maximum value of ||omega||_2 on the grid.
        if (level_number == 0)
        {
            d_omega_max = 0.0;
        }

        PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
#if (NDIM == 2)
            Pointer<CellData<NDIM,double> > omega_current_data = patch->getPatchData(d_omega_current_idx);
            d_omega_max = std::max(d_omega_max, +patch_cc_data_ops.max(omega_current_data, patch_box));
            d_omega_max = std::max(d_omega_max, -patch_cc_data_ops.min(omega_current_data, patch_box));
#endif
#if (NDIM == 3)
            Pointer<CellData<NDIM,double> > omega_norm_current_data = patch->getPatchData(d_omega_norm_current_idx);
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
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());

                Pointer<CellData<NDIM,double> > p_current_data = patch->getPatchData(d_p_current_idx);
                p_current_data->fillAll(0.0);

                Pointer<CellData<NDIM,double> > p_extrap_current_data = patch->getPatchData(d_p_extrap_current_idx);
                p_extrap_current_data->fillAll(0.0);
            }
        }
        else
        {
            // Initialize P.
            d_p_init->setDataOnPatchLevel(d_p_current_idx, d_p_var, level, init_data_time, initial_time);
            d_p_init->setDataOnPatchLevel(d_p_extrap_current_idx, d_p_extrap_var, level, init_data_time, initial_time);
        }

        // Use the initialization object to set the body force to some specified
        // value.
        if (!d_f_fcn.isNull())
        {
            d_f_fcn->setDataOnPatchLevel(d_f_current_idx, d_f_var, level, init_data_time, initial_time);
            PatchMathOps patch_math_ops;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM,double> > f_current_data = patch->getPatchData(d_f_current_idx);
                Pointer<CellData<NDIM,double> > f_cc_current_data = patch->getPatchData(d_f_cc_current_idx);
                patch_math_ops.interp(f_cc_current_data, f_current_data, patch);
            }

        }

        // Use the initialization object to set the source/sink strength to some
        // specified value.
        if (!d_q_fcn.isNull())
        {
            d_q_fcn->setDataOnPatchLevel(d_q_current_idx, d_q_var, level, init_data_time, initial_time);
        }
    }

    if (!d_adv_diff_hier_integrator.isNull())
    {
        d_adv_diff_hier_integrator->
            initializeLevelData(hierarchy, level_number, init_data_time,
                                can_be_refined, initial_time, old_level,
                                allocate_data);
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
IBImplicitHierarchyIntegrator::resetHierarchyConfiguration(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
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
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent u_bc_component(d_u_scratch_idx, SIDE_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_u_bc_coefs);
    d_u_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_u_bdry_bc_fill_op->initializeOperatorState(u_bc_component, d_hierarchy);

    // Setup the patch boundary synchronization objects.
    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
    SynchronizationTransactionComponent synch_transaction = SynchronizationTransactionComponent(d_u_scratch_idx, "CONSERVATIVE_COARSEN");
    d_side_synch_op = new SideDataSynchronization();
    d_side_synch_op->initializeOperatorState(synch_transaction, d_hierarchy);

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
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            d_rscheds[(*it).first][ln] = (*it).second->createSchedule(level, ln-1, hierarchy, d_rstrategies[(*it).first]);
        }
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it!= d_calgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(ln-1);
            d_cscheds[(*it).first][ln] = (*it).second->createSchedule(coarser_level, level, d_cstrategies[(*it).first]);
        }
    }

    // If we have added or removed a level, resize the anchor point vectors.
    d_anchor_point_local_idxs.clear();
    d_anchor_point_local_idxs.resize(finest_hier_level+1);

    // Indicate that solvers need to be re-initialized.
    d_stokes_op_needs_init = true;
    d_convective_op_needs_init = true;
    d_ib_op_needs_init = true;
    d_helmholtz_solver_needs_init = true;
    d_poisson_solver_needs_init = true;
    d_projection_pc_needs_init = true;
    d_vanka_pc_needs_init = true;
    d_block_pc_needs_init = true;
    d_ib_solver_needs_init = true;

    // Indicate that we need to perform a regrid projection.
    d_needs_regrid_projection = true;

    // Reset the hierarchy configuration for the advection-diffusion solver.
    if (!d_adv_diff_hier_integrator.isNull())
    {
        d_adv_diff_hier_integrator->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    }

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
IBImplicitHierarchyIntegrator::applyGradientDetector(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
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

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // It is necessary to untag all cells prior to tagging.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
        tags_data->fillAll(0);
    }

    if (!d_adv_diff_hier_integrator.isNull())
    {
        d_adv_diff_hier_integrator->applyGradientDetector(hierarchy, level_number, error_data_time,
                                                          tag_index, initial_time,
                                                          uses_richardson_extrapolation_too);
    }

    // Tag cells which contain Lagrangian nodes.
    d_lag_data_manager->applyGradientDetector(hierarchy, level_number, error_data_time,
                                              tag_index, initial_time,
                                              uses_richardson_extrapolation_too);

    // Tag cells based on the magnitude of the vorticity.
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
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
                Pointer<CellData<NDIM,double> > omega_current_data = patch->getPatchData(d_omega_current_idx);
                for (CellIterator<NDIM> ic(patch_box); ic; ic++)
                {
                    const Index<NDIM>& i = ic();
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
///      getExtrapolatedPressureVar(),
///      getForceVar(),
///      getSourceVar()
///
///  allows access to the various state variables maintained by the integrator.
///

Pointer<SideVariable<NDIM,double> >
IBImplicitHierarchyIntegrator::getVelocityVar()
{
    return d_u_var;
}// getVelocityVar

Pointer<CellVariable<NDIM,double> >
IBImplicitHierarchyIntegrator::getPressureVar()
{
    return d_p_var;
}// getPressureVar

Pointer<CellVariable<NDIM,double> >
IBImplicitHierarchyIntegrator::getExtrapolatedPressureVar()
{
    return d_p_extrap_var;
}// getExtrapolatedPressureVar

Pointer<SideVariable<NDIM,double> >
IBImplicitHierarchyIntegrator::getForceVar()
{
    return d_f_var;
}// getForceVar

Pointer<CellVariable<NDIM,double> >
IBImplicitHierarchyIntegrator::getSourceVar()
{
    return d_q_var;
}// getSourceVar

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
/// We simply reuse the VariableContext objects defined in the
/// AdvDiffHierarchyIntegrator object.
///

Pointer<VariableContext>
IBImplicitHierarchyIntegrator::getCurrentContext() const
{
    return d_current_context;
}// getCurrentContext

Pointer<VariableContext>
IBImplicitHierarchyIntegrator::getNewContext() const
{
    return d_new_context;
}// getNewContext

Pointer<VariableContext>
IBImplicitHierarchyIntegrator::getScratchContext() const
{
    return d_scratch_context;
}// getScratchContext

///
/// The following routines:
///
///      reinterpolateVelocity(),
///      reinterpolateForce()
///
/// are miscellaneous utility functions.

void
IBImplicitHierarchyIntegrator::reinterpolateVelocity(
    Pointer<VariableContext> ctx)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
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
IBImplicitHierarchyIntegrator::reinterpolateForce(
    Pointer<VariableContext> ctx)
{
    if (!d_f_fcn.isNull())
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int    f_idx = var_db->mapVariableAndContextToIndex(   d_f_var, ctx);
        const int f_cc_idx = var_db->mapVariableAndContextToIndex(d_f_cc_var, ctx);
        static const bool synch_cf_interface = true;
        d_hier_math_ops->interp(
            f_cc_idx, d_f_cc_var,
            f_idx   , d_f_var   ,
            d_no_fill_op, d_integrator_time, synch_cf_interface);
    }
    return;
}// reinterpolateForce

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  Serializable abstract base class.
///

void
IBImplicitHierarchyIntegrator::putToDatabase(
    Pointer<Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION",
                   IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION);

    db->putDouble("d_u_scale", d_u_scale);
    db->putDouble("d_p_scale", d_p_scale);
    db->putDouble("d_f_scale", d_f_scale);
    db->putDouble("d_q_scale", d_q_scale);
    db->putBool("d_output_u", d_output_u);
    db->putBool("d_output_p", d_output_p);
    db->putBool("d_output_f", d_output_f);
    db->putBool("d_output_q", d_output_q);
    db->putBool("d_output_omega", d_output_omega);
    db->putBool("d_output_div_u", d_output_div_u);

    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);

    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);

    db->putInteger("d_num_cycles", d_num_cycles);

    db->putInteger("d_regrid_interval", d_regrid_interval);

    db->putBool("d_using_default_tag_buffer", d_using_default_tag_buffer);
    db->putIntegerArray("d_tag_buffer", d_tag_buffer);

    db->putBool("d_using_vorticity_tagging", d_using_vorticity_tagging);
    db->putDoubleArray("d_omega_rel_thresh", d_omega_rel_thresh);
    db->putDoubleArray("d_omega_abs_thresh", d_omega_abs_thresh);
    db->putDouble("d_omega_max", d_omega_max);

    db->putBool("d_normalize_pressure", d_normalize_pressure);
    db->putString("d_convective_difference_form", enum_to_string<ConvectiveDifferencingType>(d_convective_difference_form));
    db->putBool("d_creeping_flow", d_creeping_flow);

    db->putDouble("d_old_dt", d_old_dt);
    db->putDouble("d_op_and_solver_init_dt", d_op_and_solver_init_dt);
    db->putDouble("d_integrator_time", d_integrator_time);
    db->putInteger("d_integrator_step", d_integrator_step);

    db->putDouble("d_cfl", d_cfl);

    db->putDouble("d_dt_max", d_dt_max);
    db->putDouble("d_dt_max_time_max", d_dt_max_time_max);
    db->putDouble("d_dt_max_time_min", d_dt_max_time_min);
    db->putDouble("d_dt_init", d_dt_init);

    db->putBool("d_do_log", d_do_log);

    db->putDouble("d_rho", d_rho);
    db->putDouble("d_mu", d_mu);
    db->putDouble("d_lambda", d_lambda);

    db->putDouble("d_regrid_max_div_growth_factor", d_regrid_max_div_growth_factor);

    t_put_to_database->stop();
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBImplicitHierarchyIntegrator::registerVariable(
    int& current_idx,
    int& new_idx,
    int& scratch_idx,
    const Pointer<Variable<NDIM> > variable,
    const IntVector<NDIM>& scratch_ghosts,
    const std::string& coarsen_name,
    const std::string& refine_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!variable.isNull());
#endif
    const IntVector<NDIM> no_ghosts = 0;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

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

    // Get the data transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineOperator<NDIM> > refine_operator = grid_geom->lookupRefineOperator(variable, refine_name);
    Pointer<CoarsenOperator<NDIM> > coarsen_operator = grid_geom->lookupCoarsenOperator(variable, coarsen_name);

    // Setup the refine algorithm used to fill data in new or modified patch
    // levels following a regrid operation.
    if (!refine_operator.isNull())
    {
        d_fill_after_regrid_bc_idxs.setFlag(scratch_idx);
        d_fill_after_regrid->registerRefine(current_idx, // destination
                                            current_idx, // source
                                            scratch_idx, // temporary work space
                                            refine_operator);
    }

    // Setup the SYNCH_CURRENT_STATE_DATA and SYNCH_NEW_STATE_DATA algorithms,
    // used to synchronize the data on the hierarchy.
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
IBImplicitHierarchyIntegrator::registerVariable(
    int& scratch_idx,
    const Pointer<Variable<NDIM> > variable,
    const IntVector<NDIM>& scratch_ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!variable.isNull());
#endif

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    scratch_idx = -1; // insure that uninitialized variable patch data
                      // descriptor indices cause errors

    d_scratch_variables.push_back(variable);

    // Setup the scratch context.
    scratch_idx = var_db->registerVariableAndContext(variable, getScratchContext(), scratch_ghosts);
    d_scratch_data.setFlag(scratch_idx);
    return;
}// registerVariable

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBImplicitHierarchyIntegrator::regridProjection()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate temporary data.
    ComponentSelector scratch_idxs;
    scratch_idxs.setFlag(d_u_scratch_idx);
    scratch_idxs.setFlag(d_phi_idx);
    scratch_idxs.setFlag(d_div_u_scratch_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(scratch_idxs, d_integrator_time);
    }

    // Compute div u before applying the projection operator.
    const bool u_current_cf_bdry_synch = true;
    d_hier_math_ops->div(
        d_div_u_scratch_idx, d_div_u_var, // dst
        +1.0,                             // alpha
        d_u_current_idx, d_u_var,         // src
        d_no_fill_op,                     // src_bdry_fill
        d_integrator_time,                // src_bdry_fill_time
        u_current_cf_bdry_synch);         // src_cf_bdry_synch
    if (d_do_log)
    {
        const double div_u_norm_1  = d_hier_cc_data_ops->L1Norm( d_div_u_scratch_idx, d_wgt_cc_idx);
        const double div_u_norm_2  = d_hier_cc_data_ops->L2Norm( d_div_u_scratch_idx, d_wgt_cc_idx);
        const double div_u_norm_oo = d_hier_cc_data_ops->maxNorm(d_div_u_scratch_idx, d_wgt_cc_idx);
        plog << d_object_name << "::regridProjection():\n"
             << "  performing regrid projection\n"
             << "  before projection:\n"
             << "    ||div u||_1  = " << div_u_norm_1  << "\n"
             << "    ||div u||_2  = " << div_u_norm_2  << "\n"
             << "    ||div u||_oo = " << div_u_norm_oo << "\n";
    }

    // Setup the solver vectors.
    d_hier_cc_data_ops->setToScalar(d_phi_idx, 0.0, false);
    d_hier_cc_data_ops->scale(d_div_u_scratch_idx, -1.0, d_div_u_scratch_idx);
    const double div_u_mean = (1.0/d_volume)*d_hier_cc_data_ops->integral(d_div_u_scratch_idx, d_wgt_cc_idx);
    d_hier_cc_data_ops->addScalar(d_div_u_scratch_idx, d_div_u_scratch_idx, -div_u_mean);

    SAMRAIVectorReal<NDIM,double> sol_vec(d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(d_phi_var, d_phi_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM,double> rhs_vec(d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_div_u_var, d_div_u_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

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
    d_regrid_projection_fac_op->setTime(d_integrator_time);

    d_regrid_projection_solver->setInitialGuessNonzero(false);
    d_regrid_projection_solver->setOperator(d_regrid_projection_op);

    // Solve the projection Poisson problem.
    d_regrid_projection_solver->initializeSolverState(sol_vec,rhs_vec);
    d_regrid_projection_solver->solveSystem(sol_vec,rhs_vec);
    d_regrid_projection_solver->deallocateSolverState();

    // NOTE: We always use homogeneous Neumann boundary conditions for the
    // regrid projection Poisson solver.
    d_regrid_projection_solver->setNullspace(true, NULL);

    // Setup the interpolation transaction information.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent phi_bc_component(d_phi_idx, CELL_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, &d_regrid_projection_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> phi_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    phi_bdry_bc_fill_op->initializeOperatorState(phi_bc_component, d_hierarchy);

    // Fill the physical boundary conditions for phi.
    phi_bdry_bc_fill_op->setHomogeneousBc(true);
    phi_bdry_bc_fill_op->fillData(d_integrator_time);

    // Set u := u - grad phi.
    const bool u_scratch_cf_bdry_synch = true;
    d_hier_math_ops->grad(
        d_u_scratch_idx, d_u_var,  // dst
        u_scratch_cf_bdry_synch,   // dst_cf_bdry_synch
        1.0,                       // alpha
        d_phi_idx, d_phi_var,      // src
        d_no_fill_op,              // src_bdry_fill
        d_integrator_time);        // src_bdry_fill_time
    d_hier_sc_data_ops->axpy(d_u_current_idx, -1.0, d_u_scratch_idx, d_u_current_idx);

    // Compute div u after applying the projection operator
    if (d_do_log)
    {
        const bool u_current_cf_bdry_synch = true;
        d_hier_math_ops->div(
            d_div_u_scratch_idx, d_div_u_var, // dst
            +1.0,                             // alpha
            d_u_current_idx, d_u_var,         // src
            d_no_fill_op,                     // src_bdry_fill
            d_integrator_time,                // src_bdry_fill_time
            u_current_cf_bdry_synch);         // src_cf_bdry_synch
        const double div_u_norm_1  = d_hier_cc_data_ops->L1Norm( d_div_u_scratch_idx, d_wgt_cc_idx);
        const double div_u_norm_2  = d_hier_cc_data_ops->L2Norm( d_div_u_scratch_idx, d_wgt_cc_idx);
        const double div_u_norm_oo = d_hier_cc_data_ops->maxNorm(d_div_u_scratch_idx, d_wgt_cc_idx);
        plog << "  after projection:\n"
             << "    ||div u||_1  = " << div_u_norm_1  << "\n"
             << "    ||div u||_2  = " << div_u_norm_2  << "\n"
             << "    ||div u||_oo = " << div_u_norm_oo << "\n";
    }

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }
    return;
}// regridProjection

void
IBImplicitHierarchyIntegrator::initializeOperatorsAndSolvers(
    const double current_time,
    const double new_time)
{
    const bool initial_time = MathUtilities<double>::equalEps(current_time,d_start_time);
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;

    if (d_lag_force_strategy_needs_init)
    {
        resetLagrangianForceStrategy(current_time, initial_time);
        d_lag_force_strategy_needs_init = false;
    }

    Pointer<SAMRAIVectorReal<NDIM,double> > u_scratch_vec = new SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::u_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
    u_scratch_vec->addComponent(d_u_var,d_u_scratch_idx,d_wgt_sc_idx,d_hier_sc_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > u_rhs_vec = u_scratch_vec->cloneVector(d_object_name+"::u_rhs_vec");
    const int u_rhs_idx = u_rhs_vec->getComponentDescriptorIndex(0);

    Pointer<SAMRAIVectorReal<NDIM,double> > p_scratch_vec = new SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::p_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
    p_scratch_vec->addComponent(d_p_var,d_p_scratch_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > p_rhs_vec = p_scratch_vec->cloneVector(d_object_name+"::p_rhs_vec");
    const int p_rhs_idx = p_rhs_vec->getComponentDescriptorIndex(0);

    Pointer<SAMRAIVectorReal<NDIM,double> > sol_vec = new SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec->addComponent(d_u_var,d_u_scratch_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    sol_vec->addComponent(d_p_var,d_p_scratch_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > rhs_vec = new SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec->addComponent(d_u_var,u_rhs_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
    rhs_vec->addComponent(d_p_var,p_rhs_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

    for (int d = 0; d < NDIM; ++d)
    {
        INSStaggeredVelocityBcCoef* u_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(d_u_bc_coefs[d]);
        u_bc_coef->setTimeInterval(current_time,new_time);
    }
    INSStaggeredPressureBcCoef* p_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(d_p_bc_coef);
    p_bc_coef->setTimeInterval(current_time,new_time);
    p_bc_coef->setVelocityCurrentPatchDataIndex(d_u_current_idx);
    p_bc_coef->setVelocityNewPatchDataIndex(d_u_new_idx);

    if (!d_helmholtz_solver.isNull())
    {
        d_helmholtz_spec->setCConstant((d_rho/dt)+0.5*d_lambda);
        d_helmholtz_spec->setDConstant(          -0.5*d_mu    );

        d_helmholtz_op->setPoissonSpecifications(*d_helmholtz_spec);
        d_helmholtz_op->setPhysicalBcCoefs(d_u_star_bc_coefs);
        d_helmholtz_op->setHomogeneousBc(true);
        d_helmholtz_op->setTime(new_time);
        d_helmholtz_op->setHierarchyMathOps(d_hier_math_ops);

        d_ib_mod_helmholtz_pc->setPoissonSpecifications(*d_helmholtz_spec);
        d_ib_mod_helmholtz_pc->setPhysicalBcCoefs(d_u_star_bc_coefs);
        d_ib_mod_helmholtz_pc->setHomogeneousBc(true);
        d_ib_mod_helmholtz_pc->setTime(new_time);

        if (!d_helmholtz_hypre_pc.isNull())
        {
            d_helmholtz_hypre_pc->setPoissonSpecifications(*d_helmholtz_spec);
            d_helmholtz_hypre_pc->setPhysicalBcCoefs(d_u_star_bc_coefs);
            d_helmholtz_hypre_pc->setHomogeneousBc(true);
            d_helmholtz_hypre_pc->setTime(new_time);
        }
        else if (!d_helmholtz_fac_op.isNull())
        {
            d_helmholtz_fac_op->setPoissonSpecifications(*d_helmholtz_spec);
            d_helmholtz_fac_op->setPhysicalBcCoefs(d_u_star_bc_coefs);
            d_helmholtz_fac_op->setTime(new_time);
        }

        d_helmholtz_solver->setInitialGuessNonzero(false);
        d_helmholtz_solver->setOperator(d_mod_helmholtz_op);
        if (d_helmholtz_solver_needs_init || !MathUtilities<double>::equalEps(dt,d_op_and_solver_init_dt))
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing Helmholtz solver" << std::endl;
            d_helmholtz_solver->initializeSolverState(*u_scratch_vec,*u_rhs_vec);
        }
        d_helmholtz_solver_needs_init = false;
    }

    if (!d_poisson_solver.isNull())
    {
        d_poisson_spec->setCZero();
        d_poisson_spec->setDConstant(-1.0);

        d_poisson_op->setPoissonSpecifications(*d_poisson_spec);
        d_poisson_op->setPhysicalBcCoef(d_phi_bc_coef);
        d_poisson_op->setHomogeneousBc(true);
        d_poisson_op->setTime(current_time+0.5*dt);
        d_poisson_op->setHierarchyMathOps(d_hier_math_ops);

        if (!d_poisson_hypre_pc.isNull())
        {
            d_poisson_hypre_pc->setPoissonSpecifications(*d_poisson_spec);
            d_poisson_hypre_pc->setPhysicalBcCoef(d_phi_bc_coef);
            d_poisson_hypre_pc->setHomogeneousBc(true);
            d_poisson_hypre_pc->setTime(current_time+0.5*dt);
        }
        else if (!d_poisson_fac_op.isNull())
        {
            d_poisson_fac_op->setPoissonSpecifications(*d_poisson_spec);
            d_poisson_fac_op->setPhysicalBcCoef(d_phi_bc_coef);
            d_poisson_fac_op->setTime(current_time+0.5*dt);
        }

        d_poisson_solver->setInitialGuessNonzero(false);
        d_poisson_solver->setOperator(d_poisson_op);
        if (d_poisson_solver_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing Poisson solver" << std::endl;
            d_poisson_solver->initializeSolverState(*p_scratch_vec,*p_rhs_vec);
        }
        d_poisson_solver_needs_init = false;
    }

    if (!d_projection_pc.isNull())
    {
        d_projection_pc->setTimeInterval(current_time,new_time,dt);
        if (d_projection_pc_needs_init && !d_ib_solver_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing projection preconditioner" << std::endl;
            d_projection_pc->initializeSolverState(*sol_vec,*rhs_vec);
        }
        d_projection_pc_needs_init = false;
    }

    if (!d_vanka_fac_pc.isNull())
    {
        d_vanka_fac_op->setProblemCoefficients(*d_problem_coefs,dt);
        d_vanka_fac_op->setTimeInterval(current_time,new_time);
        d_vanka_fac_op->setPhysicalBcCoefs(d_u_star_bc_coefs,d_phi_bc_coef);
        if (d_vanka_pc_needs_init && !d_ib_solver_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing Vanka preconditioner" << std::endl;
            d_vanka_fac_pc->initializeSolverState(*sol_vec,*rhs_vec);
        }
        d_vanka_pc_needs_init = false;
    }

    if (!d_block_pc.isNull())
    {
        d_block_pc->setTimeInterval(current_time,new_time);
        if (d_block_pc_needs_init && !d_ib_solver_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing block-factorization preconditioner" << std::endl;
            d_block_pc->initializeSolverState(*sol_vec,*rhs_vec);
        }
        d_block_pc_needs_init = false;
    }

    if (!d_stokes_op.isNull())
    {
        if (d_stokes_op_needs_init && !d_ib_solver_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing incompressible Stokes operator" << std::endl;
            d_stokes_op->initializeOperatorState(*u_scratch_vec,*u_rhs_vec);
        }
        d_stokes_op_needs_init = false;
    }

    if (!d_ib_op.isNull())
    {
        if (d_ib_op_needs_init && !d_ib_solver_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing implicit IB operator" << std::endl;
            d_ib_op->initializeOperatorState(*u_scratch_vec,*u_rhs_vec);
        }
        d_ib_op_needs_init = false;
    }

    if (!d_ib_solver.isNull())
    {
        d_ib_op->setTimeInterval(current_time,new_time);
        Pointer<IBImplicitSJROperator> ib_SJR_op = d_ib_SJR_op;  // XXXX
        if (!ib_SJR_op.isNull()) ib_SJR_op->setTimeInterval(current_time,new_time);
        d_ib_solver->setOperator(d_ib_op);
        d_ib_solver->setJacobian(d_ib_jac_op);
        if (d_ib_solver_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing implicit IB solver" << std::endl;
            d_ib_solver->initializeSolverState(*sol_vec,*rhs_vec);
            Pointer<PETScKrylovLinearSolver> ksp_solver = d_ib_solver->getLinearSolver();

            int ierr;
            PetscTruth flg;

            // Setup the preconditioner and preconditioner sub-solvers.
            std::vector<std::string> pc_shell_types(4);
            pc_shell_types[0] = "projection";
            pc_shell_types[1] = "vanka";
            pc_shell_types[2] = "block_factorization";
            pc_shell_types[3] = "none";
            ksp_solver->setValidPCShellTypes(pc_shell_types);

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

            char poisson_pc_type_str[len];
            ierr = PetscOptionsGetString("poisson_", "-pc_type", poisson_pc_type_str, len, &flg);  IBTK_CHKERRQ(ierr);
            std::string poisson_pc_type = "shell";
            if (flg)
            {
                poisson_pc_type = std::string(poisson_pc_type_str);
            }

            std::string stokes_pc_shell_type = "none";
            if (stokes_pc_type == "shell")
            {
                char stokes_pc_shell_type_str[len];
                ierr = PetscOptionsGetString("stokes_", "-pc_shell_type", stokes_pc_shell_type_str, len, &flg);  IBTK_CHKERRQ(ierr);
                stokes_pc_shell_type = "projection";
                if (flg)
                {
                    stokes_pc_shell_type = std::string(stokes_pc_shell_type_str);
                }

                if (!(stokes_pc_shell_type == "none" || stokes_pc_shell_type == "projection" || stokes_pc_shell_type == "vanka" || stokes_pc_shell_type == "block_factorization"))
                {
                    TBOX_ERROR(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                               "  invalid stokes shell preconditioner type: " << stokes_pc_shell_type << "\n"
                               "  valid stokes shell preconditioner types: projection, vanka, block_factorization, none" << std::endl);
                }
            }

            // Setup the velocity subdomain solver.
            const bool needs_helmholtz_solver = stokes_pc_type == "shell" && (stokes_pc_shell_type == "projection" || stokes_pc_shell_type == "block_factorization");
            if (needs_helmholtz_solver)
            {
                const std::string helmholtz_prefix = "helmholtz_";

                // Setup the various solver components.
                d_helmholtz_spec = new PoissonSpecifications(d_object_name+"::helmholtz_spec");
                d_helmholtz_op = new SCLaplaceOperator(d_object_name+"::Helmholtz Operator", *d_helmholtz_spec, d_u_star_bc_coefs, true);
                d_helmholtz_op->setHierarchyMathOps(d_hier_math_ops);

                d_mod_helmholtz_op = new IBImplicitModHelmholtzOperator(d_helmholtz_op, d_ib_SJR_op);

                d_helmholtz_solver_needs_init = true;
                d_helmholtz_solver = new PETScKrylovLinearSolver(d_object_name+"::Helmholtz Krylov Solver", helmholtz_prefix);
                d_helmholtz_solver->setInitialGuessNonzero(false);
                d_helmholtz_solver->setOperator(d_mod_helmholtz_op);
                d_helmholtz_solver->setPreconditioner(d_ib_mod_helmholtz_pc);

                // Set some default options.
                d_helmholtz_solver->setKSPType("fgmres");
                d_helmholtz_solver->setAbsoluteTolerance(1.0e-30);
                d_helmholtz_solver->setRelativeTolerance(1.0e-02);
                d_helmholtz_solver->setMaxIterations(25);
            }
            else
            {
                d_helmholtz_spec = NULL;
                d_helmholtz_op = NULL;
                d_mod_helmholtz_op = NULL;
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
                d_poisson_spec = new PoissonSpecifications(d_object_name+"::poisson_spec");
                d_poisson_op = new CCLaplaceOperator(d_object_name+"::Poisson Operator", *d_poisson_spec, d_phi_bc_coef, true);
                d_poisson_op->setHierarchyMathOps(d_hier_math_ops);

                d_poisson_solver_needs_init = true;
                d_poisson_solver = new PETScKrylovLinearSolver(d_object_name+"::Poisson Krylov Solver", poisson_prefix);
                d_poisson_solver->setInitialGuessNonzero(false);
                d_poisson_solver->setOperator(d_poisson_op);

                if (poisson_pc_type != "none")
                {
                    if (d_gridding_alg->getMaxLevels() == 1)
                    {
                        if (d_poisson_hypre_pc_db.isNull())
                        {
                            TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                                         "  Poisson hypre PC solver database is null." << std::endl);
                        }
                        d_poisson_hypre_pc = new CCPoissonHypreLevelSolver(d_object_name+"::Poisson Preconditioner", d_poisson_hypre_pc_db);
                        d_poisson_hypre_pc->setPoissonSpecifications(*d_poisson_spec);

                        d_poisson_solver->setPreconditioner(d_poisson_hypre_pc);
                    }
                    else
                    {
                        if (d_poisson_fac_pc_db.isNull())
                        {
                            TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                                         "  Poisson FAC PC solver database is null." << std::endl);
                        }
                        d_poisson_fac_op = new CCPoissonFACOperator(d_object_name+"::Poisson FAC Operator", d_poisson_fac_pc_db);
                        d_poisson_fac_op->setPoissonSpecifications(*d_poisson_spec);
                        d_poisson_fac_pc = new IBTK::FACPreconditioner(d_object_name+"::Poisson Preconditioner", *d_poisson_fac_op, d_poisson_fac_pc_db);
                        d_poisson_solver->setPreconditioner(d_poisson_fac_pc);
                    }
                }

                // Set some default options.
                d_poisson_solver->setKSPType(d_gridding_alg->getMaxLevels() == 1 ? "preonly" : "gmres");
                d_poisson_solver->setAbsoluteTolerance(1.0e-30);
                d_poisson_solver->setRelativeTolerance(1.0e-02);
                d_poisson_solver->setMaxIterations(25);
                d_poisson_solver->setNullspace(d_normalize_pressure, NULL);
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
                    d_projection_pc = new INSStaggeredProjectionPreconditioner(*d_problem_coefs, d_phi_bc_coef, d_normalize_pressure, d_helmholtz_solver, d_poisson_solver, d_hier_cc_data_ops, d_hier_sc_data_ops, d_hier_math_ops);
                    ksp_solver->setPreconditioner(d_projection_pc);
                }
                else if (stokes_pc_shell_type == "vanka")
                {
                    d_vanka_pc_needs_init = true;

                    if (d_vanka_fac_pc_db.isNull())
                    {
                        TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                                     "  Vanka FAC PC solver database is null." << std::endl);
                    }
                    d_vanka_fac_op = new INSStaggeredBoxRelaxationFACOperator(d_object_name+"::Vanka FAC Operator", d_vanka_fac_pc_db);
                    d_vanka_fac_pc = new IBTK::FACPreconditioner(d_object_name+"::Vanka Preconditioner", *d_vanka_fac_op, d_vanka_fac_pc_db);
                    ksp_solver->setPreconditioner(d_vanka_fac_pc);
                }
                else if (stokes_pc_shell_type == "block_factorization")
                {
                    d_block_pc_needs_init = true;
                    d_block_pc = new INSStaggeredBlockFactorizationPreconditioner(*d_problem_coefs, d_phi_bc_coef, d_normalize_pressure, d_helmholtz_solver, d_poisson_solver, d_hier_cc_data_ops, d_hier_sc_data_ops, d_hier_math_ops);
                    ksp_solver->setPreconditioner(d_block_pc);
                }
            }
        }
        d_ib_solver_needs_init = false;
    }

    if (!d_convective_op.isNull())
    {
        if (d_convective_op_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing convective operator" << std::endl;
            d_convective_op->initializeOperatorState(*u_scratch_vec,*u_rhs_vec);
        }
        d_convective_op_needs_init = false;
    }

    u_rhs_vec->freeVectorComponents();
    p_rhs_vec->freeVectorComponents();

    // Keep track of the timestep size to avoid unnecessary re-initialization.
    d_op_and_solver_init_dt = dt;
    return;
}// initializeOperatorsAndSolvers

double
IBImplicitHierarchyIntegrator::getLevelDt(
    Pointer<PatchLevel<NDIM> > level,
    Pointer<VariableContext> ctx) const
{
    double stable_dt = std::numeric_limits<double>::max();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        stable_dt = std::min(stable_dt,getPatchDt(patch,ctx));
    }
    stable_dt = SAMRAI_MPI::minReduction(stable_dt);
    return stable_dt;
}// getLevelDt

double
IBImplicitHierarchyIntegrator::getPatchDt(
    Pointer<Patch<NDIM> > patch,
    Pointer<VariableContext> ctx) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const Index<NDIM>& ilower = patch->getBox().lower();
    const Index<NDIM>& iupper = patch->getBox().upper();

    Pointer<SideData<NDIM,double> > u_data = patch->getPatchData(d_u_var, ctx);
    const IntVector<NDIM>& u_ghost_cells = u_data->getGhostCellWidth();

    double stable_dt = std::numeric_limits<double>::max();
#if (NDIM == 2)
    NAVIER_STOKES_SC_STABLEDT_FC(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),
        u_ghost_cells(0),u_ghost_cells(1),
        u_data->getPointer(0),u_data->getPointer(1),
        stable_dt);
#endif
#if (NDIM == 3)
    NAVIER_STOKES_SC_STABLEDT_FC(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
        u_ghost_cells(0),u_ghost_cells(1),u_ghost_cells(2),
        u_data->getPointer(0),u_data->getPointer(1),u_data->getPointer(2),
        stable_dt);
#endif
    return stable_dt;
}// getPatchDt

void
IBImplicitHierarchyIntegrator::resetLagrangianForceStrategy(
    const double init_data_time,
    const bool initial_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            d_lag_force_strategy->initializeLevelData(
                d_hierarchy, ln, init_data_time, initial_time,
                d_lag_data_manager);
        }
    }
    return;
}// resetLagrangianForceStrategy

void
IBImplicitHierarchyIntegrator::resetAnchorPointValues(
    std::vector<Pointer<LNodeLevelData> > V_data,
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
IBImplicitHierarchyIntegrator::getFromInput(
    Pointer<Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif

    // Read in data members from input database.
    d_end_time = db->getDoubleWithDefault("end_time", d_end_time);
    d_grow_dt = db->getDoubleWithDefault("grow_dt", d_grow_dt);

    d_max_integrator_steps = db->getIntegerWithDefault("max_integrator_steps", d_max_integrator_steps);

    d_num_cycles = db->getIntegerWithDefault("num_cycles", d_num_cycles);

    d_regrid_interval = db->getIntegerWithDefault("regrid_interval", d_regrid_interval);
    d_regrid_mode = string_to_enum<RegridMode>(db->getStringWithDefault("regrid_mode", enum_to_string<RegridMode>(d_regrid_mode)));

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

    d_using_vorticity_tagging = db->getBoolWithDefault("using_vorticity_tagging", d_using_vorticity_tagging);

    d_omega_rel_thresh.resizeArray(1);
    d_omega_rel_thresh[0] = 0.3;

    d_omega_abs_thresh.resizeArray(1);
    d_omega_abs_thresh[0] = std::numeric_limits<double>::max();

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
        }

        for (int i = 0; i < d_omega_abs_thresh.getSize(); ++i)
        {
            if (d_omega_abs_thresh[i] < 0.0)
            {
                TBOX_ERROR(d_object_name << ":  "
                           << "absolute vorticity thresholds for each level must be non-negative.\n");
            }
        }
    }

    d_output_u     = db->getBoolWithDefault("output_u"    , d_output_u    );
    d_output_p     = db->getBoolWithDefault("output_p"    , d_output_p    );
    d_output_f     = db->getBoolWithDefault("output_f"    , d_output_f    );
    d_output_q     = db->getBoolWithDefault("output_q"    , d_output_q    );
    d_output_omega = db->getBoolWithDefault("output_omega", d_output_omega);
    d_output_div_u = db->getBoolWithDefault("output_div_u", d_output_div_u);

    d_u_scale = db->getDoubleWithDefault("u_scale", d_u_scale);
    d_p_scale = db->getDoubleWithDefault("p_scale", d_p_scale);
    d_f_scale = db->getDoubleWithDefault("f_scale", d_f_scale);
    d_q_scale = db->getDoubleWithDefault("q_scale", d_q_scale);

    d_cfl = db->getDoubleWithDefault("cfl",d_cfl);

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault("dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault("dt_max_time_min", d_dt_max_time_min);
    d_dt_init = db->getDoubleWithDefault("dt_init",d_dt_init);

    d_helmholtz_hypre_pc_db       = db->isDatabase("HelmholtzHypreSolver") ? db->getDatabase("HelmholtzHypreSolver") : Pointer<Database>(NULL);
    d_helmholtz_petsc_pc_db       = db->isDatabase("HelmholtzPETScSolver") ? db->getDatabase("HelmholtzPETScSolver") : Pointer<Database>(NULL);
    d_helmholtz_fac_pc_db         = db->isDatabase("HelmholtzFACSolver"  ) ? db->getDatabase("HelmholtzFACSolver"  ) : Pointer<Database>(NULL);
    d_poisson_hypre_pc_db         = db->isDatabase("PoissonHypreSolver"  ) ? db->getDatabase("PoissonHypreSolver"  ) : Pointer<Database>(NULL);
    d_poisson_fac_pc_db           = db->isDatabase("PoissonFACSolver"    ) ? db->getDatabase("PoissonFACSolver"    ) : Pointer<Database>(NULL);
    d_vanka_fac_pc_db             = db->isDatabase("VankaFACSolver"      ) ? db->getDatabase("VankaFACSolver"      ) : Pointer<Database>(NULL);
    d_regrid_projection_fac_pc_db = db->isDatabase("PoissonFACSolver"    ) ? db->getDatabase("PoissonFACSolver"    ) : Pointer<Database>(NULL);

    d_regrid_max_div_growth_factor = db->getDoubleWithDefault("regrid_max_div_growth_factor", d_regrid_max_div_growth_factor);

    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    if (!is_from_restart)
    {
        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);

        d_normalize_pressure = db->getBoolWithDefault("normalize_pressure", d_normalize_pressure);
        d_convective_difference_form = string_to_enum<ConvectiveDifferencingType>(
            db->getStringWithDefault("convective_difference_form", enum_to_string<ConvectiveDifferencingType>(d_convective_difference_form)));
        d_creeping_flow = db->getBoolWithDefault("creeping_flow", d_creeping_flow);

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
        if (db->isInteger("min_ghost_cell_width"))
        {
            d_ghosts = db->getInteger("min_ghost_cell_width");
        }
        else if (db->isDouble("min_ghost_cell_width"))
        {
            d_ghosts = int(std::ceil(db->getDouble("min_ghost_cell_width")));
        }
    }
    return;
}// getFromInput

void
IBImplicitHierarchyIntegrator::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
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

    d_u_scale = db->getDouble("d_u_scale");
    d_p_scale = db->getDouble("d_p_scale");
    d_f_scale = db->getDouble("d_f_scale");
    d_q_scale = db->getDouble("d_q_scale");
    d_output_u = db->getBool("d_output_u");
    d_output_p = db->getBool("d_output_p");
    d_output_f = db->getBool("d_output_f");
    d_output_q = db->getBool("d_output_q");
    d_output_omega = db->getBool("d_output_omega");
    d_output_div_u = db->getBool("d_output_div_u");

    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    d_grow_dt = db->getDouble("d_grow_dt");

    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");

    d_num_cycles = db->getInteger("d_num_cycles");

    d_regrid_interval = db->getInteger("d_regrid_interval");

    d_using_default_tag_buffer = db->getBool("d_using_default_tag_buffer");
    d_tag_buffer = db->getIntegerArray("d_tag_buffer");

    d_using_vorticity_tagging = db->getBool("d_using_vorticity_tagging");
    d_omega_rel_thresh = db->getDoubleArray("d_omega_rel_thresh");
    d_omega_abs_thresh = db->getDoubleArray("d_omega_abs_thresh");
    d_omega_max = db->getDouble("d_omega_max");

    d_normalize_pressure = db->getBool("d_normalize_pressure");
    d_convective_difference_form = string_to_enum<ConvectiveDifferencingType>(db->getString("d_convective_difference_form"));
    d_creeping_flow = db->getBool("d_creeping_flow");

    d_old_dt = db->getDouble("d_old_dt");
    d_op_and_solver_init_dt = db->getDouble("d_op_and_solver_init_dt");
    d_integrator_time = db->getDouble("d_integrator_time");
    d_integrator_step = db->getInteger("d_integrator_step");

    d_cfl = db->getDouble("d_cfl");

    d_dt_max = db->getDouble("d_dt_max");
    d_dt_max_time_max = db->getDouble("d_dt_max_time_max");
    d_dt_max_time_min = db->getDouble("d_dt_max_time_min");
    d_dt_init = db->getDouble("d_dt_init");

    d_do_log = db->getBool("d_do_log");

    d_rho = db->getDouble("d_rho");
    d_mu = db->getDouble("d_mu");
    d_lambda = db->getDouble("d_lambda");

    d_regrid_max_div_growth_factor = db->getDouble("d_regrid_max_div_growth_factor");

    d_delta_fcn = db->getString("d_delta_fcn");
    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBImplicitHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
