// Filename: IBHierarchyIntegrator.C
// Created on 12 Jul 2004 by Boyce Griffith
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

#include "IBHierarchyIntegrator.h"

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
#include <ibamr/IBInstrumentationSpec.h>
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IndexUtilities.h>

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

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_initialize_hierarchy_integrator;
static Timer* t_initialize_hierarchy;
static Timer* t_advance_hierarchy;
static Timer* t_regrid_hierarchy;
static Timer* t_synchronize_hierarchy;
static Timer* t_synchronize_new_levels;
static Timer* t_reset_time_dependent_data;
static Timer* t_reset_data_to_preadvance_state;
static Timer* t_initialize_level_data;
static Timer* t_reset_hierarchy_configuration;
static Timer* t_apply_gradient_detector;
static Timer* t_put_to_database;

// Version of IBHierarchyIntegrator restart file data.
static const int IB_HIERARCHY_INTEGRATOR_VERSION = 1;

inline std::string
discard_comments(
    const std::string& input_string)
{
    // Create a copy of the input string, but without any text following a '!',
    // '#', or '%' character.
    std::string output_string = input_string;
    std::istringstream string_stream;

    // Discard any text following a '!' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '!');
    string_stream.clear();

    // Discard any text following a '#' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '#');
    string_stream.clear();

    // Discard any text following a '%' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '%');
    string_stream.clear();
    return output_string;
}// discard_comments

inline double
cos_delta(
    const double x,
    const double eps)
{
    if (std::abs(x) > eps)
    {
        return 0.0;
    }
    else
    {
        return 0.5*(1.0+cos(M_PI*x/eps))/eps;
    }
}// cos_delta
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBHierarchyIntegrator::IBHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<INSHierarchyIntegrator> ins_hier_integrator,
    Pointer<IBLagrangianForceStrategy> force_strategy,
    Pointer<IBPostProcessStrategy> post_processor,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_interp_delta_fcn("IB_4"),
      d_spread_delta_fcn("IB_4"),
      d_ghosts(-1),
      d_hierarchy(hierarchy),
      d_gridding_alg(NULL),
      d_visit_writer(NULL),
      d_silo_writer(NULL),
      d_load_balancer(NULL),
      d_ins_hier_integrator(ins_hier_integrator),
      d_l_data_manager(NULL),
      d_instrument_panel(NULL),
      d_total_flow_volume(),
      d_U_init(NULL),
      d_P_init(NULL),
      d_U_bc_coefs(),
      d_P_bc_coef(NULL),
      d_lag_init(NULL),
      d_body_force_fcn(NULL),
      d_eulerian_force_fcn(NULL),
      d_force_strategy(force_strategy),
      d_force_strategy_needs_init(true),
      d_post_processor(post_processor),
      d_post_processor_needs_init(true),
      d_start_time(0.0),
      d_end_time(std::numeric_limits<double>::max()),
      d_grow_dt(2.0),
      d_max_integrator_steps(std::numeric_limits<int>::max()),
      d_num_cycles(1),
      d_num_init_cycles(5),
      d_regrid_interval(1),
      d_old_dt(-1.0),
      d_integrator_time(std::numeric_limits<double>::quiet_NaN()),
      d_integrator_step(std::numeric_limits<int>::max()),
      d_dt_max(std::numeric_limits<double>::max()),
      d_dt_max_time_max(std::numeric_limits<double>::max()),
      d_dt_max_time_min(-(d_dt_max_time_max-std::numeric_limits<double>::epsilon())),
      d_is_initialized(false),
      d_do_log(false),
      d_reinterpolate_after_regrid(false),
      d_hier_cc_data_ops(),
      d_ralgs(),
      d_rstrategies(),
      d_rscheds(),
      d_calgs(),
      d_cstrategies(),
      d_cscheds(),
      d_force_current_ralg(),
      d_force_new_ralg(),
      d_force_current_rstrategy(),
      d_force_new_rstrategy(),
      d_force_current_rscheds(),
      d_force_new_rscheds(),
      d_V_var(NULL),
      d_W_var(NULL),
      d_F_var(NULL),
      d_current(NULL),
      d_scratch(NULL),
      d_V_idx(-1),
      d_W_idx(-1),
      d_F_idx(-1),
      d_F_scratch1_idx(-1),
      d_F_scratch2_idx(-1)
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
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Check the choices for the delta function.
    if (d_interp_delta_fcn != d_spread_delta_fcn)
    {
        pout << "WARNING: different delta functions are being used for velocity interpolation and force spreading.\n"
             << "         recommended usage is to employ the same delta functions for both interpolation and spreading.\n";
    }

    // Get the Lagrangian Data Manager.
    d_l_data_manager = LDataManager::getManager(d_object_name+"::LDataManager",
                                                  d_interp_delta_fcn, d_spread_delta_fcn,
                                                  d_ghosts,
                                                  d_registered_for_restart);
    d_ghosts = d_l_data_manager->getGhostCellWidth();

    // Create the instrument panel object.
    d_instrument_panel = new IBInstrumentPanel(d_object_name+"::IBInstrumentPanel", (input_db->isDatabase("IBInstrumentPanel") ? input_db->getDatabase("IBInstrumentPanel") : Pointer<Database>(NULL)));

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<CellVariable<NDIM,double> > cc_var = new CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy, true);

    // Initialize all variable contexts.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_current = var_db->getContext(d_object_name+"::CURRENT");
    d_scratch = var_db->getContext(d_object_name+"::SCRATCH");

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_initialize_hierarchy_integrator = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy            = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy               = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::advanceHierarchy()");
        t_regrid_hierarchy                = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::regridHierarchy()");
        t_synchronize_hierarchy           = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels          = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data       = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state  = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data           = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration   = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector         = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database                 = TimerManager::getManager()->getTimer("IBAMR::IBHierarchyIntegrator::putToDatabase()");
                  );
    return;
}// IBHierarchyIntegrator

IBHierarchyIntegrator::~IBHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    for (RefinePatchStrategyMap::iterator it = d_rstrategies.begin(); it != d_rstrategies.end(); ++it)
    {
        delete it->second;
    }

    for (CoarsenPatchStrategyMap::iterator it = d_cstrategies.begin(); it != d_cstrategies.end(); ++it)
    {
        delete it->second;
    }
    return;
}// ~IBHierarchyIntegrator

const std::string&
IBHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
IBHierarchyIntegrator::registerVelocityInitialConditions(
    Pointer<CartGridFunction> U_init)
{
    d_U_init = U_init;
    d_ins_hier_integrator->registerVelocityInitialConditions(d_U_init);
    return;
}// registerVelocityInitialConditions

void
IBHierarchyIntegrator::registerVelocityPhysicalBcCoefs(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& U_bc_coefs)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::registerVelocityPhysicalBcCoefs()\n"
                   << "  velocity boundary conditions must be registered prior to initialization\n"
                   << "  of the hierarchy integrator object.\n");
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(U_bc_coefs[d] != NULL);
    }
#endif
    d_U_bc_coefs = U_bc_coefs;
    d_ins_hier_integrator->registerVelocityPhysicalBcCoefs(d_U_bc_coefs);
    return;
}// registerVelocityPhysicalBcCoefs

void
IBHierarchyIntegrator::registerPressureInitialConditions(
    Pointer<CartGridFunction> P_init)
{
    d_P_init = P_init;
    d_ins_hier_integrator->registerPressureInitialConditions(d_P_init);
    return;
}// registerPressureInitialConditions

void
IBHierarchyIntegrator::registerPressurePhysicalBcCoef(
    RobinBcCoefStrategy<NDIM>* const P_bc_coef)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::registerPressurePhysicalBcCoef():\n"
                   << "  pressure boundary conditions must be registered prior to initialization\n"
                   << "  of the hierarchy integrator object." << std::endl);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(P_bc_coef != NULL);
#endif
    d_P_bc_coef = P_bc_coef;
    d_ins_hier_integrator->registerPressurePhysicalBcCoef(d_P_bc_coef);
    return;
}// registerPressurePhysicalBcCoef

void
IBHierarchyIntegrator::registerBodyForceSpecification(
    Pointer<CartGridFunction> body_force_fcn)
{
    d_body_force_fcn = body_force_fcn;
    return;
}// registerBodyForceSpecification

void
IBHierarchyIntegrator::registerLInitStrategy(
    Pointer<LInitStrategy> lag_init)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!lag_init.isNull());
#endif
    d_lag_init = lag_init;
    d_l_data_manager->registerLInitStrategy(d_lag_init);
    return;
}// registerLInitStrategy

void
IBHierarchyIntegrator::freeLInitStrategy()
{
    d_lag_init.setNull();
    d_l_data_manager->freeLInitStrategy();
    return;
}// freeLInitStrategy

void
IBHierarchyIntegrator::registerVisItDataWriter(
    Pointer<VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    d_ins_hier_integrator->registerVisItDataWriter(d_visit_writer);
    d_l_data_manager->registerVisItDataWriter(d_visit_writer);
    return;
}// registerVisItDataWriter

void
IBHierarchyIntegrator::registerLSiloDataWriter(
    Pointer<LSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!silo_writer.isNull());
#endif
    d_silo_writer = silo_writer;
    d_l_data_manager->registerLSiloDataWriter(d_silo_writer);
    return;
}// registerLSiloDataWriter

#if (NDIM == 3)
void
IBHierarchyIntegrator::registerLM3DDataWriter(
    Pointer<LM3DDataWriter> m3D_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!m3D_writer.isNull());
#endif
    d_m3D_writer = m3D_writer;
    d_l_data_manager->registerLM3DDataWriter(d_m3D_writer);
    return;
}// registerLM3DDataWriter
#endif

void
IBHierarchyIntegrator::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > load_balancer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!load_balancer.isNull());
#endif
    d_load_balancer = load_balancer;
    d_l_data_manager->registerLoadBalancer(d_load_balancer);
    return;
}// registerLoadBalancer

///
///  The following routines:
///
///      initializeHierarchyIntegrator(),
///      initializeHierarchy(),
///      advanceHierarchy(),
///      postProcessData(),
///      atRegridPoint(),
///      getIntegratorTime(),
///      getStartTime(),
///      getEndTime(),
///      getIntegratorStep(),
///      getMaxIntegratorSteps(),
///      stepsRemaining(),
///      getPatchHierarchy(),
///      getGriddingAlgorithm(),
///      getLDataManager(),
///      getIBInstrumentPanel()
///
///  allow the IBHierarchyIntegrator to be used as a hierarchy integrator.
///

void
IBHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    IBAMR_TIMER_START(t_initialize_hierarchy_integrator);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Initialize all variables.
    const IntVector<NDIM> ghosts = d_ghosts;
    const IntVector<NDIM> no_ghosts = 0;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    d_V_var = new CellVariable<NDIM,double>(d_object_name+"::V",NDIM);
    d_V_idx = var_db->registerVariableAndContext(d_V_var, d_current, ghosts);

    d_W_var = new CellVariable<NDIM,double>(d_object_name+"::W",NDIM);
    d_W_idx = var_db->registerVariableAndContext(d_W_var, d_current, ghosts);

    d_F_var = new CellVariable<NDIM,double>(d_object_name+"::F",NDIM);
    d_F_idx = var_db->registerVariableAndContext(d_F_var, d_current, no_ghosts);
    d_F_scratch1_idx = var_db->registerVariableAndContext(d_F_var, d_scratch, no_ghosts);
    d_F_scratch2_idx = var_db->registerClonedPatchDataIndex(d_F_var, d_F_scratch1_idx);

    // Initialize the objects used to manage Lagrangian-Eulerian interaction.
    //
    // NOTE: The IBEulerianForceFunction only has to set the new Cartesian grid
    // force.  The current Cartesian grid force is set manually by
    // IBHierarchyIntegrator::advanceHierarchy().
    d_eulerian_force_fcn = new IBEulerianForceFunction(d_object_name+"::IBEulerianForceFunction", -1, d_F_idx, -1);
    d_ins_hier_integrator->registerBodyForceSpecification(d_eulerian_force_fcn);

    // Initialize the INSHierarchyIntegrator.
    //
    // NOTE: This must be done after all variables created by the
    // IBHierarchyIntegrator class have been registered.
    d_ins_hier_integrator->initializeHierarchyIntegrator(d_gridding_alg);

    // Create several communications algorithms, used in filling ghost cell data
    // and synchronizing data on the patch hierarchy.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineOperator<NDIM> > refine_operator;
    Pointer<CoarsenOperator<NDIM> > coarsen_operator;

    const int U_current_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getCurrentContext());

    d_ralgs["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"] = new RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getVelocityVar(),
        "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"]->registerRefine(
        d_V_idx,       // destination
        U_current_idx, // source
        d_V_idx,       // temporary work space
        refine_operator);
    d_rstrategies["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"] = new CartExtrapPhysBdryOp(d_V_idx, "LINEAR");

    const int U_new_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getNewContext());

    d_ralgs["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"] = new RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getVelocityVar(),
        "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"]->registerRefine(
        d_W_idx,   // destination
        U_new_idx, // source
        d_W_idx,   // temporary work space
        refine_operator);
    d_rstrategies["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"] = new CartExtrapPhysBdryOp(d_W_idx, "LINEAR");

    const int U_scratch_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getScratchContext());
    const int P_new_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getPressureVar(),
        d_ins_hier_integrator->getNewContext());
    const int P_scratch_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getPressureVar(),
        d_ins_hier_integrator->getScratchContext());

    d_ralgs["INSTRUMENTATION_DATA_FILL"] = new RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getVelocityVar(),
        "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["INSTRUMENTATION_DATA_FILL"]->registerRefine(
        U_scratch_idx,  // destination
        U_new_idx,      // source
        U_scratch_idx,  // temporary work space
        refine_operator);

    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getPressureVar(),
        "LINEAR_REFINE");
    d_ralgs["INSTRUMENTATION_DATA_FILL"]->registerRefine(
        P_scratch_idx,  // destination
        P_new_idx,      // source
        P_scratch_idx,  // temporary work space
        refine_operator);

    ComponentSelector instrumentation_data_fill_bc_idxs;
    instrumentation_data_fill_bc_idxs.setFlag(U_scratch_idx);
    instrumentation_data_fill_bc_idxs.setFlag(P_scratch_idx);
    d_rstrategies["INSTRUMENTATION_DATA_FILL"] = new CartExtrapPhysBdryOp(instrumentation_data_fill_bc_idxs, "LINEAR");

    // NOTE: When using conservative averaging to coarsen the velocity from
    // finer levels to coarser levels, the appropriate prolongation operator for
    // the force is constant refinement.  This choice results in IB spreading
    // and interpolation being adjoint.
    const int F_current_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getForceVar(),
        d_ins_hier_integrator->getCurrentContext());

    d_force_current_ralg = new RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getForceVar(), "CONSTANT_REFINE");
    d_force_current_ralg->registerRefine(F_current_idx,     // destination
                                         F_current_idx,     // source
                                         d_F_scratch1_idx,  // temporary work space
                                         refine_operator);
    d_force_current_rstrategy = new CartExtrapPhysBdryOp(d_F_scratch1_idx, "LINEAR");

    d_force_new_ralg = new RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getForceVar(), "CONSTANT_REFINE");
    d_force_new_ralg->registerRefine(d_F_idx,           // destination
                                     d_F_idx,           // source
                                     d_F_scratch2_idx,  // temporary work space
                                     refine_operator);
    d_force_new_rstrategy = new CartExtrapPhysBdryOp(d_F_scratch2_idx, "LINEAR");

    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"] = new CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_V_var, "CONSERVATIVE_COARSEN");
    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"]->registerCoarsen(
        U_current_idx, // destination
        U_current_idx, // source
        coarsen_operator);

    d_calgs["U->U::N->N::CONSERVATIVE_COARSEN"] = new CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_W_var, "CONSERVATIVE_COARSEN");
    d_calgs["U->U::N->N::CONSERVATIVE_COARSEN"]->registerCoarsen(
        U_new_idx, // destination
        U_new_idx, // source
        coarsen_operator);

    d_calgs["F->F::CONSERVATIVE_COARSEN"] = new CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_ins_hier_integrator->getForceVar(), "CONSERVATIVE_COARSEN");
    d_calgs["F->F::CONSERVATIVE_COARSEN"]->registerCoarsen(
        F_current_idx, // destination
        F_current_idx, // source
        coarsen_operator);
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_F_var, "CONSERVATIVE_COARSEN");
    d_calgs["F->F::CONSERVATIVE_COARSEN"]->registerCoarsen(
        d_F_idx, // destination
        d_F_idx, // source
        coarsen_operator);

    // Set the current integration time.
    if (!RestartManager::getManager()->isFromRestart())
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }

    // Indicate that the integrator has been initialized.
    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_hierarchy_integrator);
    return;
}// initializeHierarchyIntegrator

double
IBHierarchyIntegrator::initializeHierarchy()
{
    IBAMR_TIMER_START(t_initialize_hierarchy);

    // Use the INSHierarchyIntegrator to initialize the patch hierarchy.
    double dt_next = d_ins_hier_integrator->initializeHierarchy();

    if (d_integrator_time >= d_dt_max_time_min && d_integrator_time <= d_dt_max_time_max)
    {
        dt_next = std::min(dt_next, d_dt_max);
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Begin Lagrangian data movement.
    d_l_data_manager->beginDataRedistribution();

    // Finish Lagrangian data movement.
    d_l_data_manager->endDataRedistribution();

    // Update the workload.
    d_l_data_manager->updateWorkloadData(coarsest_ln, finest_ln);

    // Initialize the instrumentation data.
    d_instrument_panel->initializeHierarchyIndependentData(d_hierarchy, d_l_data_manager);
    if (d_instrument_panel->isInstrumented())
    {
        d_instrument_panel->initializeHierarchyDependentData(d_hierarchy, d_l_data_manager, d_integrator_step, d_integrator_time);
        if (d_total_flow_volume.empty())
        {
            d_total_flow_volume.resize(d_instrument_panel->getFlowValues().size(),0.0);
        }
    }

    // Indicate that the force strategy and the post processor need to be
    // re-initialized.
    d_force_strategy_needs_init  = true;
    d_post_processor_needs_init  = true;

    IBAMR_TIMER_STOP(t_initialize_hierarchy);
    return dt_next;
}// initializeHierarchy

double
IBHierarchyIntegrator::advanceHierarchy(
    const double dt)
{
    IBAMR_TIMER_START(t_advance_hierarchy);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
#endif

    const double current_time = d_integrator_time;
    const double new_time     = d_integrator_time+dt;
    const bool initial_time   = MathUtilities<double>::equalEps(current_time,d_start_time);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Regrid the patch hierarchy.
    const bool do_regrid = (d_regrid_interval == 0 ? false : (d_integrator_step % d_regrid_interval == 0));
    if (do_regrid)
    {
        if (d_do_log) plog << d_object_name << "::advanceHierarchy(): regridding prior to timestep " << d_integrator_step << "\n";
        regridHierarchy();
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Set the current time interval in the force specification objects.
    d_eulerian_force_fcn->setTimeInterval(current_time, new_time);
    d_force_strategy->setTimeInterval(current_time, new_time);

    // (Re)initialize the force strategy and the post-processor.
    if (d_force_strategy_needs_init)
    {
        resetLagrangianForceStrategy(current_time, initial_time);
        d_force_strategy_needs_init = false;
    }
    if (d_post_processor_needs_init && !d_post_processor.isNull())
    {
        resetPostProcessor(current_time, initial_time);
        d_post_processor_needs_init = false;
    }

    std::vector<Pointer<LData> > X_data(finest_ln+1);
    std::vector<Pointer<LData> > F_data(finest_ln+1);
    std::vector<Pointer<LData> > U_data(finest_ln+1);
    std::vector<Pointer<LData> > X_new_data(finest_ln+1);
    std::vector<Pointer<LData> > F_new_data(finest_ln+1);
    std::vector<Pointer<LData> > U_new_data(finest_ln+1);

    std::vector<Pointer<LData> > K_data(finest_ln+1);
    std::vector<Pointer<LData> > M_data(finest_ln+1);
    std::vector<Pointer<LData> > Y_data(finest_ln+1);
    std::vector<Pointer<LData> > dY_dt_data(finest_ln+1);
    std::vector<Pointer<LData> > F_K_data(finest_ln+1);
    std::vector<Pointer<LData> > Y_new_data(finest_ln+1);
    std::vector<Pointer<LData> > dY_dt_new_data(finest_ln+1);
    std::vector<Pointer<LData> > F_K_new_data(finest_ln+1);

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
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(d_V_idx, current_time);
            d_rscheds["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(current_time);
        }
    }

    // Initialize the various LData objects, including X_data, F_data,
    // U_data, and X_new_data, on each level of the patch hierarchy.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            X_data    [ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
            U_data    [ln] = d_l_data_manager->getLData(LDataManager::VEL_DATA_NAME,ln);
            F_data    [ln] = d_l_data_manager->createLData("F",ln,NDIM);

            X_new_data[ln] = d_l_data_manager->createLData("X_new",ln,NDIM);
            U_new_data[ln] = d_l_data_manager->createLData("U_new",ln,NDIM);
            F_new_data[ln] = d_l_data_manager->createLData("F_new",ln,NDIM);
        }
    }

    // 1. Interpolate u(n) from the Cartesian grid onto the Lagrangian mesh.
    //
    // NOTE: Since we are maintaining the Lagrangian velocity data, this step is
    // skipped for each timestep following the initial one except immediately
    // following a regridding operation.
    if (initial_time || d_reinterpolate_after_regrid)
    {
        if (d_do_log) plog << d_object_name << "::advanceHierarchy(): interpolating u(n) to U(n)\n";
        d_l_data_manager->interp(d_V_idx, U_data, X_data);
    }

    resetAnchorPointValues(U_data, coarsest_ln, finest_ln);

    // 2. Compute F(n) = F(X(n),n), the Lagrangian force corresponding to
    //    configuration X(n) at time t_{n}.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            if (d_do_log) plog << d_object_name << "::advanceHierarchy(): computing F(n) on level number " << ln << "\n";
            Vec F_vec = F_data[ln]->getVec();
            int ierr = VecSet(F_vec, 0.0);  IBTK_CHKERRQ(ierr);
            d_force_strategy->computeLagrangianForce(
                F_data[ln], X_data[ln], U_data[ln],
                d_hierarchy, ln, current_time, d_l_data_manager);
        }
    }

    resetAnchorPointValues(F_data, coarsest_ln, finest_ln);

    // 3. Spread F(n) from the Lagrangian mesh onto the Cartesian grid.
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): spreading F(n) to f(n)\n";
    const int F_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getForceVar(), d_ins_hier_integrator->getCurrentContext());
    d_l_data_manager->spread(F_current_idx, F_data, X_data);

    // 4. Compute X~(n+1), the preliminary structure configuration at time
    //    t_{n+1}.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            if (d_do_log) plog << d_object_name << "::advanceHierarchy(): computing X~(n+1) on level number " << ln << "\n";
            Vec U_vec = U_data[ln]->getVec();
            Vec X_vec = X_data[ln]->getVec();
            Vec X_new_vec = X_new_data[ln]->getVec();
            int ierr = VecWAXPY(X_new_vec, dt, U_vec, X_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // 5. Compute F~(n+1) = F(X~(n+1),n+1), the Lagrangian force corresponding
    //    to configuration X~(n+1) at time t_{n+1}.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            if (d_do_log) plog << d_object_name << "::advanceHierarchy(): computing F~(n+1) on level number " << ln << "\n";
            Vec F_new_vec = F_new_data[ln]->getVec();
            int ierr = VecSet(F_new_vec, 0.0);  IBTK_CHKERRQ(ierr);
            d_force_strategy->computeLagrangianForce(
                F_new_data[ln], X_new_data[ln], U_data[ln],
                d_hierarchy, ln, new_time, d_l_data_manager);
        }
    }

    resetAnchorPointValues(F_new_data, coarsest_ln, finest_ln);

    // 6. Spread F~(n+1) from the Lagrangian mesh onto the Cartesian grid.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
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
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM,double> > f_new_data = patch->getPatchData(d_F_idx);
                f_new_data->fillAll(0.0);
            }
        }
        else
        {
            // Interpolate the Cartesian force from the next coarser level in
            // the hierarchy.
            //
            // Note that these refine schedules initialize both the force data
            // maintained by the IBHierarchyIntegrator and the current context
            // of the force data maintained by the Navier-Stokes solver.
            d_force_new_rscheds[ln]->fillData(current_time);
        }

        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            if (d_do_log) plog << d_object_name << "::advanceHierarchy(): spreading F~(n+1) to f~(n+1) on level number " << ln << "\n";
            d_l_data_manager->spread(d_F_idx, F_new_data, X_new_data);
        }
    }

    // 7. If an additional body force specification object is provided, compute
    //    the body force at the beginning and end of the time step, and add
    //    those values to f(n) and f(n+1).
    if (!d_body_force_fcn.isNull())
    {
        d_body_force_fcn->setDataOnPatchHierarchy(
            d_F_scratch1_idx, d_F_var, d_hierarchy,
            current_time, false, coarsest_ln, finest_ln);

        d_hier_cc_data_ops->add(F_current_idx, F_current_idx, d_F_scratch1_idx);

        d_body_force_fcn->setDataOnPatchHierarchy(
            d_F_scratch2_idx, d_F_var, d_hierarchy,
            current_time+dt, false, coarsest_ln, finest_ln);

        d_hier_cc_data_ops->add(d_F_idx, d_F_idx, d_F_scratch2_idx);
    }

    // Synchronize the Cartesian grid force field on the patch hierarchy.
    //
    // Note that these coarsen schedules synchronize both the force data
    // maintained by the IBHierarchyIntegrator and the current context of the
    // force data maintained by the Navier-Stokes solver.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["F->F::CONSERVATIVE_COARSEN"][ln]->coarsenData();
    }

    // Deallocate unneeded scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
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
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): computing u(n+1), p(n+1/2)\n";

    int cycle = 0;
    const bool performing_init_cycles = initial_time;
    const int num_cycles = (performing_init_cycles
                            ? d_num_init_cycles
                            : d_num_cycles);
    for (cycle = 0; cycle < num_cycles; ++cycle)
    {
        if (d_do_log && performing_init_cycles)
        {
            plog << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++\n";
            plog << "+\n";
            plog << "+ Performing cycle " << cycle+1 << " of " << d_num_init_cycles << " to initialize P(n=1/2)\n";
            plog << "+\n";
            plog << "++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
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
                Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                level->setTime(current_time);
            }
        }
        else
        {
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
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
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_W_idx, current_time);
        d_rscheds["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(current_time);
    }

    // Interpolate u(n+1) from the Cartesian grid onto the Lagrangian mesh using
    // X~(n+1).
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): interpolating u(n+1) to U(n+1)\n";
    d_l_data_manager->interp(d_W_idx, U_new_data, X_new_data);
    resetAnchorPointValues(U_new_data, coarsest_ln, finest_ln);

    // Compute X(n+1), the final structure configuration at time t_{n+1}.
    //
    // Also, advance the positions and velocities of the massive ghost particles
    // forward in time via 2nd order SSP Runge-Kutta.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            int ierr;

            if (d_do_log) plog << d_object_name << "::advanceHierarchy(): computing X(n+1) on level number " << ln << "\n";

            Vec U_new_vec = U_new_data[ln]->getVec();
            Vec X_vec = X_data[ln]->getVec();
            Vec X_new_vec = X_new_data[ln]->getVec();
            ierr = VecAXPY(X_new_vec, dt, U_new_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(X_vec, 0.5, 0.5, X_new_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Interpolate u(n+1) from the Cartesian grid onto the Lagrangian mesh using
    // X(n+1).
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): interpolating u(n+1) to U(n+1)\n";
    d_l_data_manager->interp(d_W_idx, U_data, X_data);

    // Update the instrumentation data.
    updateIBInstrumentationData(d_integrator_step+1,new_time);
    if (d_instrument_panel->isInstrumented())
    {
        const std::vector<std::string>& instrument_name = d_instrument_panel->getInstrumentNames();
        const std::vector<double>& flow_data = d_instrument_panel->getFlowValues();
        for (unsigned int m = 0; m < flow_data.size(); ++m)
        {
            // NOTE: Flow volume is calculated in default units.
            d_total_flow_volume[m] += flow_data[m]*dt;
            plog << "flow volume through " << instrument_name[m] << ":\t " << d_total_flow_volume[m] << "\n";
        }
    }

    // Deallocate any remaining scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_W_idx);
        level->deallocatePatchData(d_F_idx);
    }

    // Reset all time dependent data.
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
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

    IBAMR_TIMER_STOP(t_advance_hierarchy);
    return dt_next;
}// advanceHierarchy

void
IBHierarchyIntegrator::postProcessData()
{
    if (d_post_processor.isNull()) return;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getCurrentContext());
    const int P_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVar(), d_ins_hier_integrator->getCurrentContext());
    const int F_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getForceVar(), d_ins_hier_integrator->getCurrentContext());

    const double current_time = d_integrator_time;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Initialize data on each level of the patch hierarchy.
    std::vector<Pointer<LData> > X_data(finest_ln+1);
    std::vector<Pointer<LData> > F_data(finest_ln+1);
    std::vector<Pointer<LData> > U_data(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_V_idx, current_time);
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
            U_data[ln] = d_l_data_manager->getLData(LDataManager::VEL_DATA_NAME,ln);
            F_data[ln] = d_l_data_manager->createLData("F",ln,NDIM);
        }
    }

    // Interpolate u(n) from the Cartesian grid onto the Lagrangian mesh.
    d_l_data_manager->interp(d_V_idx, U_data, X_data);
    resetAnchorPointValues(U_data, coarsest_ln, finest_ln);

    // Compute F(n) = F(X(n),n), the Lagrangian force corresponding to
    // configuration X(n) at time t_{n}.
    if (d_force_strategy_needs_init)
    {
        const bool initial_time = MathUtilities<double>::equalEps(current_time,d_start_time);
        resetLagrangianForceStrategy(current_time, initial_time);
        d_force_strategy_needs_init = false;
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec F_vec = F_data[ln]->getVec();
            int ierr = VecSet(F_vec, 0.0);  IBTK_CHKERRQ(ierr);
            d_force_strategy->computeLagrangianForce(
                F_data[ln], X_data[ln], U_data[ln],
                d_hierarchy, ln, current_time, d_l_data_manager);
        }
    }
    resetAnchorPointValues(F_data, coarsest_ln, finest_ln);

    // Perform the user-defined post-processing.
    d_post_processor->postProcessData(U_current_idx, P_current_idx, F_current_idx, F_data, X_data, U_data, d_hierarchy, coarsest_ln, finest_ln, current_time, d_l_data_manager);

    // Deallocate data on each level of the patch hierarchy.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_V_idx);
    }
    return;
}// postProcessData

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

const Pointer<PatchHierarchy<NDIM> >
IBHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

Pointer<GriddingAlgorithm<NDIM> >
IBHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

LDataManager*
IBHierarchyIntegrator::getLDataManager() const
{
    return d_l_data_manager;
}// getLDataManager

Pointer<IBInstrumentPanel>
IBHierarchyIntegrator::getIBInstrumentPanel() const
{
    return d_instrument_panel;
}// getIBInstrumentPanel

///
///  The following routines:
///
///      regridHierarchy(),
///      synchronizeHierarchy(),
///      synchronizeNewLevels(),
///      resetTimeDependentHierData(),
///      resetHierDataToPreadvanceState()
///
///  allow the IBHierarchyIntegrator to provide data management for a time
///  integrator which making use of this class.
///

void
IBHierarchyIntegrator::regridHierarchy()
{
    IBAMR_TIMER_START(t_regrid_hierarchy);

    // Determine the current range of hierarchy levels.
    const int coarsest_ln_before_regrid = 0;
    const int finest_ln_before_regrid = d_hierarchy->getFinestLevelNumber();

    // Update the workload pre-regridding.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_l_data_manager->updateWorkloadData(coarsest_ln_before_regrid,finest_ln_before_regrid);

    // Before regridding, begin Lagrangian data movement.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): starting Lagrangian data movement.\n";
    d_l_data_manager->beginDataRedistribution(coarsest_ln_before_regrid,finest_ln_before_regrid);

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): calling INSHierarchyIntegrator::regridHierarchy().\n";
    d_ins_hier_integrator->regridHierarchy();

    // After regridding, finish Lagrangian data movement.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): finishing Lagrangian data movement.\n";
    d_l_data_manager->endDataRedistribution(coarsest_ln_before_regrid,finest_ln_before_regrid);

    // Update the workload post-regridding.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_l_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());

    // Indicate that the force strategy and post-processor need to be
    // re-initialized.
    d_force_strategy_needs_init = true;
    d_post_processor_needs_init = true;

    // Look up the re-distributed Lagrangian position data.
    std::vector<Pointer<LData> > X_data(d_hierarchy->getFinestLevelNumber()+1);
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
        }
    }

    // Compute the set of local anchor points.
    static const double eps = 2.0*sqrt(std::numeric_limits<double>::epsilon());
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const grid_xLower = grid_geom->getXLower();
    const double* const grid_xUpper = grid_geom->getXUpper();
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        d_anchor_point_local_idxs[ln].clear();
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
            const std::vector<LNode*>& local_nodes = mesh->getNodes();
            for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
                 cit != local_nodes.end(); ++cit)
            {
                const LNode* const node_idx = *cit;
                const IBAnchorPointSpec* const anchor_point_spec = node_idx->getNodeDataItem<IBAnchorPointSpec>();
                if (anchor_point_spec != NULL)
                {
                    d_anchor_point_local_idxs[ln].insert(node_idx->getLocalPETScIndex());
                }
            }

            const blitz::Array<double,2>& X_array = *X_data[ln]->getLocalFormVecArray();
            for (int i = 0; i < static_cast<int>(X_data[ln]->getLocalNodeCount()); ++i)
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    if ((periodic_shift[d] == 0) &&
                        (X_array(i,d) <= grid_xLower[d]+eps || X_array(i,d) >= grid_xUpper[d]-eps))
                    {
                        d_anchor_point_local_idxs[ln].insert(i);
                        break;
                    }
                }
            }
            X_data[ln]->restoreArrays();
        }
    }

    IBAMR_TIMER_STOP(t_regrid_hierarchy);
    return;
}// regridHierarchy

void
IBHierarchyIntegrator::synchronizeHierarchy()
{
    IBAMR_TIMER_START(t_synchronize_hierarchy);

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->synchronizeHierarchy();

    IBAMR_TIMER_STOP(t_synchronize_hierarchy);
    return;
}// synchronizeHierarchy

void
IBHierarchyIntegrator::synchronizeNewLevels(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const bool initial_time)
{
    IBAMR_TIMER_START(t_synchronize_new_levels);

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
    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->synchronizeNewLevels(
        hierarchy, coarsest_level, finest_level,
        sync_time, initial_time);

    IBAMR_TIMER_STOP(t_synchronize_new_levels);
    return;
}// synchronizeNewLevels

void
IBHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    IBAMR_TIMER_START(t_reset_time_dependent_data);

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->resetTimeDependentHierData(new_time);

    IBAMR_TIMER_STOP(t_reset_time_dependent_data);
    return;
}// resetTimeDependentHierData

void
IBHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    IBAMR_TIMER_START(t_reset_data_to_preadvance_state);

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->resetHierDataToPreadvanceState();

    IBAMR_TIMER_STOP(t_reset_data_to_preadvance_state);
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
///  StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
IBHierarchyIntegrator::initializeLevelData(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
    IBAMR_TIMER_START(t_initialize_level_data);

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

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->initializeLevelData(
        hierarchy, level_number, init_data_time,
        can_be_refined, initial_time, old_level,
        allocate_data);

    // We use the LDataManager to handle as much unstructured data management as
    // possible.
    d_l_data_manager->setPatchHierarchy(hierarchy);
    d_l_data_manager->resetLevels(0,hierarchy->getFinestLevelNumber());

    d_l_data_manager->initializeLevelData(
        hierarchy, level_number, init_data_time,
        can_be_refined, initial_time, old_level,
        allocate_data);

    IBAMR_TIMER_STOP(t_initialize_level_data);
    return;
}// initializeLevelData

void
IBHierarchyIntegrator::resetHierarchyConfiguration(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    IBAMR_TIMER_START(t_reset_hierarchy_configuration);

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

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // We use the LDataManager to handle as much unstructured data management as
    // possible.
    d_l_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Indicate that the velocity field needs to be re-interpolated (but only in
    // the multi-level case).
    d_reinterpolate_after_regrid = d_reinterpolate_after_regrid || (finest_level>0);

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    // If we have added or removed a level, resize the anchor point vectors.
    d_anchor_point_local_idxs.clear();
    d_anchor_point_local_idxs.resize(finest_hier_level+1);

    // If we have added or removed a level, resize the schedule vectors.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin(); it != d_ralgs.end(); ++it)
    {
        d_rscheds[it->first].resize(finest_hier_level+1);
    }

    d_force_current_rscheds.resize(finest_hier_level+1);
    d_force_new_rscheds    .resize(finest_hier_level+1);

    for (CoarsenAlgMap::const_iterator it = d_calgs.begin(); it != d_calgs.end(); ++it)
    {
        d_cscheds[it->first].resize(finest_hier_level+1);
    }

    // (Re)build generic refine communication schedules.  These are created for
    // all levels in the hierarchy.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin(); it != d_ralgs.end(); ++it)
    {
        for (int ln = coarsest_level; ln <= finest_hier_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

            d_rscheds[it->first][ln] = it->second->createSchedule(
                level, ln-1, hierarchy, d_rstrategies[it->first]);
        }
    }

    // (Re)build specialized refine communication schedules used to compute the
    // Cartesian force density.  These are set only for levels >= 1.
    for (int ln = std::max(1,coarsest_level); ln <= finest_hier_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

        d_force_current_rscheds[ln] = d_force_current_ralg->createSchedule(
            level, Pointer<PatchLevel<NDIM> >(),
            ln-1, hierarchy, d_force_current_rstrategy);

        d_force_new_rscheds[ln] = d_force_new_ralg->createSchedule(
            level, Pointer<PatchLevel<NDIM> >(),
            ln-1, hierarchy, d_force_new_rstrategy);
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin(); it != d_calgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(ln-1);
            d_cscheds[it->first][ln] = it->second->createSchedule(
                coarser_level, level, d_cstrategies[it->first]);
        }
    }

    IBAMR_TIMER_STOP(t_reset_hierarchy_configuration);
    return;
}// resetHierarchyConfiguration

void
IBHierarchyIntegrator::applyGradientDetector(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    IBAMR_TIMER_START(t_apply_gradient_detector);

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

    // Tag cells for refinement according to the criteria specified by the
    // INSHierarchyIntegrator.
    d_ins_hier_integrator->applyGradientDetector(
        hierarchy, level_number, error_data_time,
        tag_index, initial_time,
        uses_richardson_extrapolation_too);

    // Tag cells which contain Lagrangian nodes.
    d_l_data_manager->applyGradientDetector(
        hierarchy, level_number, error_data_time,
        tag_index, initial_time,
        uses_richardson_extrapolation_too);

    IBAMR_TIMER_STOP(t_apply_gradient_detector);
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
/// We simply reuse the VariableContext objects defined in the
/// INSHierarchyIntegrator object.
///

Pointer<VariableContext>
IBHierarchyIntegrator::getCurrentContext() const
{
    return d_ins_hier_integrator->getCurrentContext();
}// getCurrentContext

Pointer<VariableContext>
IBHierarchyIntegrator::getNewContext() const
{
    return d_ins_hier_integrator->getNewContext();
}// getNewContext

Pointer<VariableContext>
IBHierarchyIntegrator::getScratchContext() const
{
    return d_ins_hier_integrator->getScratchContext();
}// getScratchContext

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  Serializable abstract base class.
///

void
IBHierarchyIntegrator::putToDatabase(
    Pointer<Database> db)
{
    IBAMR_TIMER_START(t_put_to_database);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("IB_HIERARCHY_INTEGRATOR_VERSION",
                   IB_HIERARCHY_INTEGRATOR_VERSION);

    db->putString("d_interp_delta_fcn", d_interp_delta_fcn);
    db->putString("d_spread_delta_fcn", d_spread_delta_fcn);
    db->putInteger("d_total_flow_volume_sz", d_total_flow_volume.size());
    if (!d_total_flow_volume.empty())
    {
        db->putDoubleArray("d_total_flow_volume", &d_total_flow_volume[0], d_total_flow_volume.size());
    }
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

    IBAMR_TIMER_STOP(t_put_to_database);
    return;
}// putToDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBHierarchyIntegrator::resetLagrangianForceStrategy(
    const double init_data_time,
    const bool initial_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            d_force_strategy->initializeLevelData(
                d_hierarchy, ln, init_data_time, initial_time,
                d_l_data_manager);
        }
    }
    return;
}// resetLagrangianForceStrategy

void
IBHierarchyIntegrator::resetPostProcessor(
    const double init_data_time,
    const bool initial_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            d_post_processor->initializeLevelData(
                d_hierarchy, ln, init_data_time, initial_time,
                d_l_data_manager);
        }
    }
    return;
}// resetPostProcessor

void
IBHierarchyIntegrator::updateIBInstrumentationData(
    const int timestep_num,
    const double data_time)
{
    if (!d_instrument_panel->isInstrumented()) return;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Compute the positions of the flow meter nets.
    d_instrument_panel->initializeHierarchyDependentData(
        d_hierarchy, d_l_data_manager, timestep_num, data_time);

    // Compute the flow rates and pressures.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int U_scratch_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getScratchContext());
    const int P_scratch_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getPressureVar(),
        d_ins_hier_integrator->getScratchContext());

    std::vector<bool> deallocate_U_scratch_data(finest_ln+1,false);
    std::vector<bool> deallocate_P_scratch_data(finest_ln+1,false);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(U_scratch_idx))
        {
            deallocate_U_scratch_data[ln] = true;
            level->allocatePatchData(U_scratch_idx, data_time);
        }
        if (!level->checkAllocated(P_scratch_idx))
        {
            deallocate_P_scratch_data[ln] = true;
            level->allocatePatchData(P_scratch_idx, data_time);
        }
        d_rscheds["INSTRUMENTATION_DATA_FILL"][ln]->fillData(data_time);
    }

    d_instrument_panel->readInstrumentData(
        U_scratch_idx, P_scratch_idx,
        d_hierarchy, d_l_data_manager,
        timestep_num, data_time);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (deallocate_U_scratch_data[ln]) level->deallocatePatchData(U_scratch_idx);
        if (deallocate_P_scratch_data[ln]) level->deallocatePatchData(P_scratch_idx);
    }
    return;
}// updateIBInstrumentationData

void
IBHierarchyIntegrator::resetAnchorPointValues(
    std::vector<Pointer<LData> > V_data,
    const int coarsest_ln,
    const int finest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            const int depth = V_data[ln]->getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(depth == NDIM);
#endif
            Vec V_vec = V_data[ln]->getVec();
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
IBHierarchyIntegrator::getFromInput(
    Pointer<Database> db,
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

    d_num_cycles = db->getIntegerWithDefault("num_cycles", d_num_cycles);

    d_regrid_interval = db->getIntegerWithDefault(
        "regrid_interval", d_regrid_interval);

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault(
        "dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault(
        "dt_max_time_min", d_dt_max_time_min);

    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    if (!is_from_restart)
    {
        if (db->isString("interp_delta_fcn") && db->isString("spread_delta_fcn"))
        {
            d_interp_delta_fcn = db->getStringWithDefault("interp_delta_fcn", d_interp_delta_fcn);
            d_spread_delta_fcn = db->getStringWithDefault("spread_delta_fcn", d_spread_delta_fcn);
        }
        else
        {
            d_interp_delta_fcn = db->getStringWithDefault("delta_fcn", d_interp_delta_fcn);
            d_spread_delta_fcn = db->getStringWithDefault("delta_fcn", d_spread_delta_fcn);
        }

        if (db->isInteger("min_ghost_cell_width"))
        {
            d_ghosts = db->getInteger("min_ghost_cell_width");
        }
        else if (db->isDouble("min_ghost_cell_width"))
        {
            d_ghosts = static_cast<int>(std::ceil(db->getDouble("min_ghost_cell_width")));
        }

        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);

        d_num_init_cycles = db->getIntegerWithDefault(
            "num_init_cycles", d_num_init_cycles);
    }
    return;
}// getFromInput

void
IBHierarchyIntegrator::getFromRestart()
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

    int ver = db->getInteger("IB_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    if (db->isString("d_interp_delta_fcn") && db->isString("d_spread_delta_fcn"))
    {
        d_interp_delta_fcn = db->getString("d_interp_delta_fcn");
        d_spread_delta_fcn = db->getString("d_spread_delta_fcn");
    }
    else if (db->isString("d_delta_fcn"))
    {
        d_interp_delta_fcn = db->getString("d_delta_fcn");
        d_spread_delta_fcn = db->getString("d_delta_fcn");
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to "
                   << d_object_name << " does not contain keys ``d_interp_delta_fcn'', ``d_spread_delta_fcn'', or ``d_delta_fcn''.");
    }

    const int total_flow_volume_sz = db->getInteger("d_total_flow_volume_sz");
    d_total_flow_volume.resize(total_flow_volume_sz, std::numeric_limits<double>::quiet_NaN());
    if (!d_total_flow_volume.empty())
    {
        db->getDoubleArray("d_total_flow_volume", &d_total_flow_volume[0], d_total_flow_volume.size());
    }
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

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
