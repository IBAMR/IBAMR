// Filename: IBStaggeredHierarchyIntegrator.C
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
#include <ibamr/IBInstrumentationSpec.h>
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CartSideDoubleSpecializedLinearRefine.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IndexUtilities.h>
#if (NDIM == 3)
#include <ibtk/LM3DDataWriter.h>
#endif
#include <ibtk/LMarkerSetData.h>
#include <ibtk/LMarkerUtilities.h>
#include <ibtk/LSiloDataWriter.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CoarsenOperator.h>
#include <HierarchyDataOpsManager.h>
#include <Index.h>
#include <Patch.h>
#include <VariableDatabase.h>
#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// BLITZ++ INCLUDES
#include <blitz/array.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <iterator>
#include <limits>
#include <numeric>

// FORTRAN ROUTINES
#define DSQRTM_FC FC_FUNC(dsqrtm,DSQRTM)
extern "C"
{
    void
    DSQRTM_FC(
        double* X, const double* A);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
void
interpolate_directors(
    const std::vector<double*>& D_half,
    const std::vector<double*>& D,
    const std::vector<double*>& D_new)
{
    blitz::Array<double,2> A(3,3,blitz::ColumnMajorArray<2>());
    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < 3; ++i)
        {
            A(i,j) = D_new[0][i]*D[0][j] + D_new[1][i]*D[1][j] + D_new[2][i]*D[2][j];
        }
    }

    blitz::Array<double,2> sqrt_A(3,3,blitz::ColumnMajorArray<2>());
    DSQRTM_FC(sqrt_A.data(), A.data());

    for (int alpha = 0; alpha < 3; ++alpha)
    {
        blitz::Array<double,1> D_half_arr(D_half[alpha], blitz::shape(3), blitz::neverDeleteData);
        blitz::Array<double,1> D_arr(D[alpha], blitz::shape(3), blitz::neverDeleteData);

        blitz::firstIndex i;
        blitz::secondIndex j;
        D_half_arr = blitz::sum(sqrt_A(i,j)*D_arr(j),j);
    }
    return;
}// interpolate_directors

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

// Version of IBStaggeredHierarchyIntegrator restart file data.
static const int IB_STAGGERED_HIERARCHY_INTEGRATOR_VERSION = 1;

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

IBStaggeredHierarchyIntegrator::IBStaggeredHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
    Pointer<IBLagrangianForceStrategy> force_strategy,
    Pointer<IBLagrangianSourceStrategy> source_strategy,
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
      d_lag_init(NULL),
      d_body_force_fcn(NULL),
      d_eulerian_force_fcn(NULL),
      d_force_strategy(force_strategy),
      d_force_strategy_needs_init(true),
      d_eulerian_source_fcn(NULL),
      d_source_strategy(source_strategy),
      d_source_strategy_needs_init(true),
      d_post_processor(post_processor),
      d_post_processor_needs_init(true),
      d_using_pIB_method(false),
      d_gravitational_acceleration(0.0),
      d_using_orthonormal_directors(false),
      d_start_time(0.0),
      d_end_time(std::numeric_limits<double>::max()),
      d_grow_dt(2.0),
      d_max_integrator_steps(std::numeric_limits<int>::max()),
      d_num_cycles(2),
      d_regrid_interval(1),
      d_regrid_cfl_interval(0.0),
      d_regrid_cfl_estimate(0.0),
      d_old_dt(-1.0),
      d_integrator_time(std::numeric_limits<double>::quiet_NaN()),
      d_integrator_step(std::numeric_limits<int>::max()),
      d_dt_max(std::numeric_limits<double>::max()),
      d_dt_max_time_max(std::numeric_limits<double>::max()),
      d_dt_max_time_min(-(d_dt_max_time_max-std::numeric_limits<double>::epsilon())),
      d_is_initialized(false),
      d_do_log(false),
      d_mark_input_file_name(""),
      d_mark_init_posns(),
      d_hier_cc_data_ops(),
      d_hier_sc_data_ops(),
      d_ralgs(),
      d_rstrategies(),
      d_rscheds(),
      d_palgs(),
      d_pstrategies(),
      d_pscheds(),
      d_calgs(),
      d_cstrategies(),
      d_cscheds(),
      d_V_var(NULL),
      d_W_var(NULL),
      d_F_var(NULL),
      d_N_var(NULL),
      d_mark_var(NULL),
      d_current(NULL),
      d_scratch(NULL),
      d_V_idx(-1),
      d_W_idx(-1),
      d_F_idx(-1),
      d_N_idx(-1),
      d_mark_current_idx(-1),
      d_mark_scratch_idx(-1)
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

    // Read in the marker initial positions.
    if (!from_restart) LMarkerUtilities::readMarkerPositions(d_mark_init_posns, d_mark_input_file_name, d_hierarchy);

    // Create the instrument panel object.
    d_instrument_panel = new IBInstrumentPanel(d_object_name+"::IBInstrumentPanel", (input_db->isDatabase("IBInstrumentPanel") ? input_db->getDatabase("IBInstrumentPanel") : Pointer<Database>(NULL)));

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<CellVariable<NDIM,double> > cc_var = new CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy, true);
    Pointer<SideVariable<NDIM,double> > sc_var = new SideVariable<NDIM,double>("sc_var");
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, hierarchy, true);

    // Initialize all variable contexts.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_current = var_db->getContext(d_object_name+"::CURRENT");
    d_scratch = var_db->getContext(d_object_name+"::SCRATCH");

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_initialize_hierarchy_integrator = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy            = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy               = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::advanceHierarchy()");
        t_regrid_hierarchy                = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::regridHierarchy()");
        t_synchronize_hierarchy           = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels          = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data       = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state  = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data           = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration   = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector         = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database                 = TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::putToDatabase()");
                  );
    return;
}// IBStaggeredHierarchyIntegrator

IBStaggeredHierarchyIntegrator::~IBStaggeredHierarchyIntegrator()
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
}// ~IBStaggeredHierarchyIntegrator

const std::string&
IBStaggeredHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
IBStaggeredHierarchyIntegrator::registerVelocityInitialConditions(
    Pointer<CartGridFunction> U_init)
{
    d_ins_hier_integrator->registerVelocityInitialConditions(U_init);
    return;
}// registerVelocityInitialConditions

void
IBStaggeredHierarchyIntegrator::registerVelocityPhysicalBcCoefs(
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
    d_ins_hier_integrator->registerVelocityPhysicalBcCoefs(U_bc_coefs);
    return;
}// registerVelocityPhysicalBcCoefs

void
IBStaggeredHierarchyIntegrator::registerBodyForceSpecification(
    Pointer<CartGridFunction> body_force_fcn)
{
    d_body_force_fcn = body_force_fcn;
    return;
}// registerBodyForceSpecification

void
IBStaggeredHierarchyIntegrator::registerLInitStrategy(
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
IBStaggeredHierarchyIntegrator::freeLInitStrategy()
{
    d_lag_init.setNull();
    d_l_data_manager->freeLInitStrategy();
    return;
}// freeLInitStrategy

void
IBStaggeredHierarchyIntegrator::registerVisItDataWriter(
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
IBStaggeredHierarchyIntegrator::registerLSiloDataWriter(
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
IBStaggeredHierarchyIntegrator::registerLM3DDataWriter(
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
IBStaggeredHierarchyIntegrator::registerLoadBalancer(
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
///  allow the IBStaggeredHierarchyIntegrator to be used as a hierarchy integrator.
///

void
IBStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Initialize all variables.
    const IntVector<NDIM> ghosts = d_ghosts;
    const IntVector<NDIM> no_ghosts = 0;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    d_V_var = new SideVariable<NDIM,double>(d_object_name+"::V");
    d_V_idx = var_db->registerVariableAndContext(d_V_var, d_scratch, ghosts);

    d_F_var = new SideVariable<NDIM,double>(d_object_name+"::F");
    d_F_idx = var_db->registerVariableAndContext(d_F_var, d_scratch, SIDEG);

    if (d_using_orthonormal_directors)
    {
        d_W_var = new SideVariable<NDIM,double>(d_object_name+"::W");
        d_W_idx = var_db->registerVariableAndContext(d_W_var, d_scratch, ghosts);

        d_N_var = new SideVariable<NDIM,double>(d_object_name+"::N");
        d_N_idx = var_db->registerVariableAndContext(d_N_var, d_scratch, SIDEG);
    }

    if (!d_source_strategy.isNull())
    {
        d_Q_var = new CellVariable<NDIM,double>(d_object_name+"::Q",1);
        d_Q_idx = var_db->registerVariableAndContext(d_Q_var, d_scratch, no_ghosts);
    }

    d_mark_var = new LMarkerSetVariable(d_object_name+"::mark");
    d_mark_current_idx = var_db->registerVariableAndContext(d_mark_var, getCurrentContext(), ghosts);
    d_mark_scratch_idx = var_db->registerVariableAndContext(d_mark_var, getScratchContext(), ghosts);
    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_mark_current_idx);
    }

    // Initialize the objects used to manage Lagrangian-Eulerian interaction.
    d_eulerian_force_fcn = new IBEulerianForceFunction(d_object_name+"::IBEulerianForceFunction", -1, -1, d_F_idx);
    d_ins_hier_integrator->registerBodyForceSpecification(d_eulerian_force_fcn);

    if (!d_source_strategy.isNull())
    {
        d_eulerian_source_fcn = new IBEulerianSourceFunction(d_object_name+"::IBEulerianSourceFunction", d_Q_idx, d_Q_idx, d_Q_idx);
        d_ins_hier_integrator->registerSourceSpecification(d_eulerian_source_fcn);
    }

    // Initialize the INSStaggeredHierarchyIntegrator.
    //
    // NOTE: This must be done after all variables created by the
    // IBStaggeredHierarchyIntegrator class have been registered.
    d_ins_hier_integrator->initializeHierarchyIntegrator(d_gridding_alg);

    // Create several communications algorithms, used in filling ghost cell data
    // and synchronizing data on the patch hierarchy.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    IBAMR_DO_ONCE(
        grid_geom->addSpatialRefineOperator(new CartSideDoubleSpecializedLinearRefine());
                  );

    Pointer<RefineOperator<NDIM> > refine_operator;
    Pointer<CoarsenOperator<NDIM> > coarsen_operator;

    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getCurrentContext());
    const int U_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getNewContext());
    const int U_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getScratchContext());
    const int P_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVar(), d_ins_hier_integrator->getNewContext());
    const int P_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVar(), d_ins_hier_integrator->getScratchContext());

    d_ralgs["U->V::C->S::GHOST_FILL"] = new RefineAlgorithm<NDIM>();
    refine_operator = NULL;
    d_ralgs["U->V::C->S::GHOST_FILL"]->registerRefine(d_V_idx, U_current_idx, d_V_idx, refine_operator);

    d_ralgs["V->V::S->S::GHOST_FILL"] = new RefineAlgorithm<NDIM>();
    refine_operator = NULL;
    d_ralgs["V->V::S->S::GHOST_FILL"]->registerRefine(d_V_idx, d_V_idx, d_V_idx, refine_operator);

    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"] = new CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(d_V_var, "CONSERVATIVE_COARSEN");
    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"]->registerCoarsen(U_current_idx, U_current_idx, coarsen_operator);

    d_calgs["V->V::S->S::CONSERVATIVE_COARSEN"] = new CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(d_V_var, "CONSERVATIVE_COARSEN");
    d_calgs["V->V::S->S::CONSERVATIVE_COARSEN"]->registerCoarsen(d_V_idx, d_V_idx, coarsen_operator);

    d_palgs["F->F::S->S::PROLONGATION"] = new RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(d_F_var, "SPECIALIZED_LINEAR_REFINE");
    d_palgs["F->F::S->S::PROLONGATION"]->registerRefine(d_F_idx, d_F_idx, d_F_idx, refine_operator);

    d_ralgs["INSTRUMENTATION_DATA_FILL"] = new RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getVelocityVar(), "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["INSTRUMENTATION_DATA_FILL"]->registerRefine(U_scratch_idx, U_new_idx, U_scratch_idx, refine_operator);
    refine_operator = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getPressureVar(), "LINEAR_REFINE");
    d_ralgs["INSTRUMENTATION_DATA_FILL"]->registerRefine(P_scratch_idx, P_new_idx, P_scratch_idx, refine_operator);
    ComponentSelector instrumentation_data_fill_bc_idxs;
    instrumentation_data_fill_bc_idxs.setFlag(U_scratch_idx);
    instrumentation_data_fill_bc_idxs.setFlag(P_scratch_idx);
    d_rstrategies["INSTRUMENTATION_DATA_FILL"] = new CartExtrapPhysBdryOp(instrumentation_data_fill_bc_idxs, "LINEAR");

    if (d_using_orthonormal_directors)
    {
        d_ralgs["N->N::S->S::CONSERVATIVE_LINEAR_REFINE"] = new RefineAlgorithm<NDIM>();
        refine_operator = grid_geom->lookupRefineOperator(d_N_var, "CONSERVATIVE_LINEAR_REFINE");
        d_ralgs["N->N::S->S::CONSERVATIVE_LINEAR_REFINE"]->registerRefine(d_N_idx, d_N_idx, d_N_idx, refine_operator);

        d_palgs["N->N::S->S::PROLONGATION"] = new RefineAlgorithm<NDIM>();
        refine_operator = grid_geom->lookupRefineOperator(d_N_var, "SPECIALIZED_LINEAR_REFINE");
        d_palgs["N->N::S->S::PROLONGATION"]->registerRefine(d_N_idx, d_N_idx, d_N_idx, refine_operator);

        d_ralgs["W->W::S->S::GHOST_FILL"] = new RefineAlgorithm<NDIM>();
        refine_operator = NULL;
        d_ralgs["W->W::S->S::GHOST_FILL"]->registerRefine(d_W_idx, d_W_idx, d_W_idx, refine_operator);

        d_calgs["W->W::S->S::CONSERVATIVE_COARSEN"] = new CoarsenAlgorithm<NDIM>();
        coarsen_operator = grid_geom->lookupCoarsenOperator(d_W_var, "CONSERVATIVE_COARSEN");
        d_calgs["W->W::S->S::CONSERVATIVE_COARSEN"]->registerCoarsen(d_W_idx, d_W_idx, coarsen_operator);
    }

    // Set the current integration time.
    if (!RestartManager::getManager()->isFromRestart())
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

    // Use the INSStaggeredHierarchyIntegrator to initialize the patch
    // hierarchy.
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

    // Initialize various Lagrangian data objects.
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time,d_start_time);
    const double init_data_time = d_integrator_time;
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        const bool can_be_refined = level_number < finest_ln || d_gridding_alg->levelCanBeRefined(level_number);

        // Initialize marker data.
        LMarkerUtilities::initializeMarkersOnLevel(d_mark_current_idx, d_mark_init_posns, d_hierarchy, level_number, initial_time, Pointer<PatchLevel<NDIM> >(NULL));

        // Initialize source/sink data.
        if (!d_source_strategy.isNull())
        {
            d_source_strategy->initializeLevelData(d_hierarchy, level_number, init_data_time, initial_time, d_l_data_manager);
            d_n_src[level_number] = d_source_strategy->getNumSources(d_hierarchy, level_number, d_integrator_time, d_l_data_manager);
            d_X_src[level_number].resize(d_n_src[level_number], blitz::TinyVector<double,NDIM>(std::numeric_limits<double>::quiet_NaN()));
            d_r_src[level_number].resize(d_n_src[level_number], std::numeric_limits<double>::quiet_NaN());
            d_P_src[level_number].resize(d_n_src[level_number], std::numeric_limits<double>::quiet_NaN());
            d_Q_src[level_number].resize(d_n_src[level_number], std::numeric_limits<double>::quiet_NaN());
            Pointer<LData> X_data = NULL;
            if (d_l_data_manager->levelContainsLagrangianData(level_number)) X_data = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,level_number);
            d_source_strategy->getSourceLocations(d_X_src[level_number], d_r_src[level_number], X_data, d_hierarchy, level_number, init_data_time, d_l_data_manager);
        }

        // Remaining initialization is done only for levels that contain
        // Lagrangian data.
        if (!d_l_data_manager->levelContainsLagrangianData(level_number)) continue;

        // Setup data needed for pIB method.
        if (initial_time && d_using_pIB_method)
        {
            static const bool manage_data = true;
            Pointer<LData> M_data = d_l_data_manager->createLData("M",level_number,1,manage_data);
            Pointer<LData> K_data = d_l_data_manager->createLData("K",level_number,1,manage_data);
            Pointer<LData> Y_data = d_l_data_manager->createLData("Y",level_number,NDIM,manage_data);
            Pointer<LData> dY_dt_data = d_l_data_manager->createLData("dY_dt",level_number,NDIM,manage_data);
            static const int global_index_offset = 0;
            static const int local_index_offset = 0;
            d_lag_init->initializeMassDataOnPatchLevel(global_index_offset, local_index_offset, M_data, K_data, d_hierarchy, level_number, init_data_time, can_be_refined, initial_time, d_l_data_manager);
            if (!d_silo_writer.isNull())
            {
                d_silo_writer->registerVariableData("M", M_data, level_number);
                d_silo_writer->registerVariableData("Y", Y_data, level_number);
            }
        }

        // Setup data needed for gIB method.
        if (initial_time && d_using_orthonormal_directors)
        {
            static const bool manage_data = true;
            Pointer<LData> D_data = d_l_data_manager->createLData("D",level_number,3*3,manage_data);
            static const int global_index_offset = 0;
            static const int local_index_offset = 0;
            d_lag_init->initializeDirectorDataOnPatchLevel(global_index_offset, local_index_offset, D_data, d_hierarchy, level_number, init_data_time, can_be_refined, initial_time, d_l_data_manager);
            if (!d_silo_writer.isNull())
            {
                d_silo_writer->registerVariableData("D1", D_data, 0, 3, level_number);
                d_silo_writer->registerVariableData("D2", D_data, 3, 3, level_number);
                d_silo_writer->registerVariableData("D3", D_data, 6, 3, level_number);
            }
        }
    }

    // Prune duplicate markers following initialization.
    LMarkerUtilities::pruneDuplicateMarkers(d_mark_current_idx, d_hierarchy);

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

    // Indicate that the force and source strategies and the post processor need
    // to be re-initialized.
    d_force_strategy_needs_init = true;
    d_source_strategy_needs_init = true;
    d_post_processor_needs_init = true;

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

    const double current_time = d_integrator_time;
    const double new_time     = d_integrator_time+dt;
    const bool initial_time   = MathUtilities<double>::equalEps(current_time,d_start_time);

    Pointer<HierarchyMathOps> hier_math_ops = d_ins_hier_integrator->getHierarchyMathOps();

    // Regrid the patch hierarchy.
    if (atRegridPoint())
    {
        if (d_do_log) plog << d_object_name << "::advanceHierarchy(): regridding prior to timestep " << d_integrator_step << "\n";
        regridHierarchy();
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Set the current time interval in the force and (optional) source
    // specification objects.
    d_eulerian_force_fcn->registerBodyForceSpecification(d_body_force_fcn);
    d_eulerian_force_fcn->setTimeInterval(current_time, new_time);
    d_force_strategy->setTimeInterval(current_time, new_time);
    if (!d_source_strategy.isNull())
    {
        d_eulerian_source_fcn->setTimeInterval(current_time, new_time);
        d_source_strategy->setTimeInterval(current_time, new_time);
    }

    // (Re)initialize the force and (optional) source strategies and the
    // post-processor.
    if (d_force_strategy_needs_init)
    {
        resetLagrangianForceStrategy(current_time, initial_time);
        d_force_strategy_needs_init = false;
    }
    if (d_source_strategy_needs_init && !d_source_strategy.isNull())
    {
        resetLagrangianSourceStrategy(current_time, initial_time);
        d_source_strategy_needs_init = false;
    }
    if (d_post_processor_needs_init && !d_post_processor.isNull())
    {
        resetPostProcessor(current_time, initial_time);
        d_post_processor_needs_init = false;
    }

    // Allocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_V_idx, current_time);
        level->allocatePatchData(d_F_idx, current_time);
        if (d_using_orthonormal_directors)
        {
            level->allocatePatchData(d_W_idx, current_time);
            level->allocatePatchData(d_N_idx, current_time);
        }
        if (!d_source_strategy.isNull())
        {
            level->allocatePatchData(d_Q_idx, current_time);
        }
    }

    // Initialize the various LData objects on each level of the patch hierarchy.
    std::vector<Pointer<LData> > X_data(finest_ln+1);
    std::vector<Pointer<LData> > U_data(finest_ln+1);
    std::vector<Pointer<LData> > X_new_data(finest_ln+1);
    std::vector<Pointer<LData> > X_half_data(finest_ln+1);
    std::vector<Pointer<LData> > U_half_data(finest_ln+1);
    std::vector<Pointer<LData> > F_half_data(finest_ln+1);

    std::vector<Pointer<LData> > K_data(finest_ln+1);
    std::vector<Pointer<LData> > M_data(finest_ln+1);
    std::vector<Pointer<LData> > Y_data(finest_ln+1);
    std::vector<Pointer<LData> > dY_dt_data(finest_ln+1);
    std::vector<Pointer<LData> > Y_new_data(finest_ln+1);
    std::vector<Pointer<LData> > dY_dt_new_data(finest_ln+1);
    std::vector<Pointer<LData> > F_K_half_data(finest_ln+1);

    std::vector<Pointer<LData> > D_data(finest_ln+1);
    std::vector<Pointer<LData> > D_new_data(finest_ln+1);
    std::vector<Pointer<LData> > D_half_data(finest_ln+1);
    std::vector<Pointer<LData> > N_half_data(finest_ln+1);
    std::vector<Pointer<LData> > W_half_data(finest_ln+1);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
            U_data[ln] = d_l_data_manager->getLData(LDataManager::VEL_DATA_NAME,ln);
            X_new_data[ln] = d_l_data_manager->createLData("X_new",ln,NDIM);
            X_half_data[ln] = d_l_data_manager->createLData("X_half",ln,NDIM);
            U_half_data[ln] = d_l_data_manager->createLData("U_half",ln,NDIM);
            F_half_data[ln] = d_l_data_manager->createLData("F_half",ln,NDIM);

            if (d_using_pIB_method)
            {
                K_data[ln]         = d_l_data_manager->getLData("K",ln);
                M_data[ln]         = d_l_data_manager->getLData("M",ln);
                Y_data[ln]         = d_l_data_manager->getLData("Y",ln);
                dY_dt_data[ln]     = d_l_data_manager->getLData("dY_dt",ln);
                Y_new_data[ln]     = d_l_data_manager->createLData("Y_new",ln,NDIM);
                dY_dt_new_data[ln] = d_l_data_manager->createLData("dY_dt_new",ln,NDIM);
                F_K_half_data[ln]  = d_l_data_manager->createLData("F_K_half",ln,NDIM);
            }

            if (d_using_orthonormal_directors)
            {
                D_data[ln]      = d_l_data_manager->getLData("D",ln);
                D_new_data[ln]  = d_l_data_manager->createLData("D_new",ln,3*3);
                D_half_data[ln] = d_l_data_manager->createLData("D_half",ln,3*3);
                N_half_data[ln] = d_l_data_manager->createLData("N_half",ln,3);
                W_half_data[ln] = d_l_data_manager->createLData("W_half",ln,3);
            }
        }
    }

    // Get patch data descriptors for the current and new velocity data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getCurrentContext());
    const int U_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getNewContext());

    // Synchronize the Cartesian grid velocity u(n) on the patch hierarchy, then
    // interpolate the velocity field.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["U->U::C->C::CONSERVATIVE_COARSEN"][ln]->coarsenData();
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_rscheds["U->V::C->S::GHOST_FILL"][ln]->fillData(current_time);
    }
    d_l_data_manager->interp(d_V_idx, U_data, X_data);
    resetAnchorPointValues(U_data, coarsest_ln, finest_ln);

    // Initialize X(n+1/2) to equal X(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_vec = X_data[ln]->getVec();
            Vec X_half_vec = X_half_data[ln]->getVec();
            int ierr = VecCopy(X_vec, X_half_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Initialize X(n+1) to equal X(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_vec = X_data[ln]->getVec();
            Vec X_new_vec = X_new_data[ln]->getVec();
            int ierr = VecCopy(X_vec, X_new_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Initialize X_mark(n+1) to equal X_mark(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_mark_scratch_idx, d_integrator_time);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<LMarkerSetData> mark_data         = patch->getPatchData(d_mark_current_idx);
            Pointer<LMarkerSetData> mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);
            mark_scratch_data->copy(*mark_data);
        }
    }

    // Initialize U(n+1/2) to equal U(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec U_vec = U_data[ln]->getVec();
            Vec U_half_vec = U_half_data[ln]->getVec();
            int ierr = VecCopy(U_vec, U_half_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    if (d_using_pIB_method)
    {
        // Set the initial values of Y and dY/dt.
        if (initial_time)
        {
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                if (d_l_data_manager->levelContainsLagrangianData(ln))
                {
                    int ierr;
                    Vec X_vec = X_data[ln]->getVec();
                    Vec U_vec = U_data[ln]->getVec();
                    Vec Y_vec = Y_data[ln]->getVec();
                    Vec dY_dt_vec = dY_dt_data[ln]->getVec();
                    ierr = VecCopy(X_vec, Y_vec);  IBTK_CHKERRQ(ierr);
                    ierr = VecCopy(U_vec, dY_dt_vec);  IBTK_CHKERRQ(ierr);
                }
            }
        }

        // Initialize Y(n+1) to equal Y(n).
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                Vec Y_vec = Y_data[ln]->getVec();
                Vec Y_new_vec = Y_new_data[ln]->getVec();
                int ierr = VecCopy(Y_vec, Y_new_vec);  IBTK_CHKERRQ(ierr);
            }
        }

        // Initialize dY/dt(n+1) to equal dY/dt(n).
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                Vec dY_dt_vec = dY_dt_data[ln]->getVec();
                Vec dY_dt_new_vec = dY_dt_new_data[ln]->getVec();
                int ierr = VecCopy(dY_dt_vec, dY_dt_new_vec);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    if (d_using_orthonormal_directors)
    {
        // Initialize D(n+1/2) to equal D(n).
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                Vec D_vec = D_data[ln]->getVec();
                Vec D_half_vec = D_half_data[ln]->getVec();
                int ierr = VecCopy(D_vec, D_half_vec);  IBTK_CHKERRQ(ierr);
            }
        }

        // Initialize D(n+1) to equal D(n).
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                Vec D_vec = D_data[ln]->getVec();
                Vec D_new_vec = D_new_data[ln]->getVec();
                int ierr = VecCopy(D_vec, D_new_vec);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    // Perform one or more cycles of fixed point iteration to compute the
    // updated configuration of the coupled fluid-structure system.
    d_ins_hier_integrator->integrateHierarchy_initialize(current_time, new_time);
    for (int cycle = 0; cycle < d_num_cycles; ++cycle)
    {
        // Compute F(n+1/2) = F(X(n+1/2),U(n+1/2),t_{n+1/2}).
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                Vec F_half_vec = F_half_data[ln]->getVec();
                int ierr = VecSet(F_half_vec, 0.0);  IBTK_CHKERRQ(ierr);
                if (!d_using_orthonormal_directors)
                {
                    d_force_strategy->computeLagrangianForce(F_half_data[ln], X_half_data[ln], U_half_data[ln], d_hierarchy, ln, current_time+0.5*dt, d_l_data_manager);
                }
                else
                {
                    Vec N_half_vec = N_half_data[ln]->getVec();
                    int ierr = VecSet(N_half_vec, 0.0);  IBTK_CHKERRQ(ierr);
                    d_force_strategy->computeLagrangianForceAndTorque(F_half_data[ln], N_half_data[ln], X_half_data[ln], D_half_data[ln], U_half_data[ln], d_hierarchy, ln, current_time+0.5*dt, d_l_data_manager);
                }
            }
        }

        if (d_using_pIB_method)
        {
            // Compute pIB-related penalty forces, F_K = K*(Y-X), and update the
            // pIB-related state variables, Y and dY/dt.
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                if (d_l_data_manager->levelContainsLagrangianData(ln))
                {
                    int ierr;

                    Vec K_vec = K_data[ln]->getVec();
                    Vec M_vec = M_data[ln]->getVec();
                    Vec X_half_vec = X_half_data[ln]->getVec();
                    Vec Y_vec = Y_data[ln]->getVec();
                    Vec Y_new_vec = Y_new_data[ln]->getVec();
                    Vec dY_dt_vec = dY_dt_data[ln]->getVec();
                    Vec dY_dt_new_vec = dY_dt_new_data[ln]->getVec();
                    Vec F_K_half_vec = F_K_half_data[ln]->getVec();

                    int n_local = 0;
                    ierr = VecGetLocalSize(M_vec, &n_local);  IBTK_CHKERRQ(ierr);

                    double* K_arr, * M_arr, * X_half_arr, * Y_arr, * Y_new_arr, * dY_dt_arr, * dY_dt_new_arr, * F_K_half_arr;
                    ierr = VecGetArray(K_vec, &K_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecGetArray(M_vec, &M_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecGetArray(X_half_vec, &X_half_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecGetArray(Y_vec, &Y_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecGetArray(Y_new_vec, &Y_new_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecGetArray(dY_dt_vec, &dY_dt_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecGetArray(dY_dt_new_vec, &dY_dt_new_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecGetArray(F_K_half_vec, &F_K_half_arr);  IBTK_CHKERRQ(ierr);

                    static double max_displacement = 0.0;
                    double max_config_displacement = 0.0;
                    for (int i = 0; i < n_local; ++i)
                    {
                        const double& K = K_arr[i];
                        const double& M = M_arr[i];
                        const double* const X_half = &X_half_arr[NDIM*i];
                        const double* const Y = &Y_arr[NDIM*i];
                        const double* const dY_dt = &dY_dt_arr[NDIM*i];
                        double* const Y_new = &Y_new_arr[NDIM*i];
                        double* const dY_dt_new = &dY_dt_new_arr[NDIM*i];
                        double* const F_K_half = &F_K_half_arr[NDIM*i];

                        double displacement = 0.0;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            double Y_minus_X = 0.5*(Y_new[d] + Y[d]) - X_half[d];
                            displacement += Y_minus_X*Y_minus_X;
                            F_K_half[d] = K*Y_minus_X;
                            Y_new[d] = Y[d] + 0.5*dt*(dY_dt_new[d]+dY_dt[d]);
                            dY_dt_new[d] = dY_dt[d] - (dt/M)*F_K_half[d] + dt*d_gravitational_acceleration[d];
                        }
                        displacement = sqrt(displacement);
                        max_config_displacement = std::max(displacement,max_config_displacement);
                    }
                    max_config_displacement = SAMRAI_MPI::maxReduction(max_config_displacement);

                    max_displacement = std::max(max_config_displacement,max_displacement);
                    if (d_do_log)
                    {
                        plog << d_object_name << "::advanceHierarchy():\n";
                        plog << "  maximum massive boundary point displacement [present configuration] = " << max_config_displacement << "\n";
                        plog << "  maximum massive boundary point displacement [entire simulation]     = " << max_displacement << "\n";
                    }

                    ierr = VecRestoreArray(K_vec, &K_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecRestoreArray(M_vec, &M_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecRestoreArray(X_half_vec, &X_half_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecRestoreArray(Y_vec, &Y_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecRestoreArray(Y_new_vec, &Y_new_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecRestoreArray(dY_dt_vec, &dY_dt_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecRestoreArray(dY_dt_new_vec, &dY_dt_new_arr);  IBTK_CHKERRQ(ierr);
                    ierr = VecRestoreArray(F_K_half_vec, &F_K_half_arr);  IBTK_CHKERRQ(ierr);

                    Vec F_half_vec = F_half_data[ln]->getVec();
                    ierr = VecAXPY(F_half_vec,1.0,F_K_half_vec);  IBTK_CHKERRQ(ierr);
                }
            }
        }

        // Spread F(n+1/2) (and possibly N(n+1/2)) to f(n+1/2).
        d_hier_sc_data_ops->setToScalar(d_F_idx, 0.0);
        if (d_using_orthonormal_directors)
        {
            d_hier_sc_data_ops->setToScalar(d_N_idx, 0.0, false);
            resetAnchorPointValues(N_half_data, coarsest_ln, finest_ln);
            d_l_data_manager->spread(d_N_idx, N_half_data, X_half_data, d_pscheds["N->N::S->S::PROLONGATION"], true, true);
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                d_rscheds["N->N::S->S::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(current_time);
            }
            hier_math_ops->curl(d_F_idx, d_F_var, d_N_idx, d_N_var, NULL, current_time);
            d_hier_sc_data_ops->scale(d_F_idx, 0.5, d_F_idx);
        }
        resetAnchorPointValues(F_half_data, coarsest_ln, finest_ln);
        d_l_data_manager->spread(d_F_idx, F_half_data, X_half_data, d_pscheds["F->F::S->S::PROLONGATION"], true, true);

        // Compute the source/sink strengths corresponding to any distributed
        // internal fluid sources or sinks.
        if (!d_source_strategy.isNull())
        {
            computeSourceStrengths(coarsest_ln, finest_ln, current_time+0.5*dt, X_half_data);
        }

        // Solve the incompressible Navier-Stokes equations.
        d_ins_hier_integrator->integrateHierarchy(current_time, new_time, cycle);

        // Set u(n+1/2) = 0.5*(u(n) + u(n+1)).
        d_hier_sc_data_ops->linearSum(d_V_idx, 0.5, U_current_idx, 0.5, U_new_idx);

        // Interpolate u(n+1/2) to U(n+1/2).
        d_l_data_manager->interp(d_V_idx, U_half_data, X_half_data, d_cscheds["V->V::S->S::CONSERVATIVE_COARSEN"], d_rscheds["V->V::S->S::GHOST_FILL"], current_time);
        resetAnchorPointValues(U_half_data, coarsest_ln, finest_ln);

        // Set X(n+1) = X(n) + dt*U(n+1/2).
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                int ierr;
                Vec X_vec = X_data[ln]->getVec();
                Vec X_new_vec = X_new_data[ln]->getVec();
                Vec U_half_vec = U_half_data[ln]->getVec();
                ierr = VecWAXPY(X_new_vec, dt, U_half_vec, X_vec);  IBTK_CHKERRQ(ierr);
            }
        }

        // Set X_mark(n+1) = X_mark(n) + dt*U(n+1/2).
        LMarkerUtilities::advectMarkers(d_mark_current_idx, d_mark_scratch_idx, d_V_idx, dt, d_interp_delta_fcn, d_hierarchy);

        // Set X(n+1/2) = 0.5*(X(n) + X(n+1)).
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                int ierr;
                Vec X_vec = X_data[ln]->getVec();
                Vec X_new_vec = X_new_data[ln]->getVec();
                Vec X_half_vec = X_half_data[ln]->getVec();
                ierr = VecCopy(X_vec, X_half_vec);  IBTK_CHKERRQ(ierr);
                ierr = VecAXPBY(X_half_vec, 0.5, 0.5, X_new_vec);  IBTK_CHKERRQ(ierr);
            }
        }

        if (d_using_orthonormal_directors)
        {
            // Set w(n+1/2) = curl u(n+1/2).
            hier_math_ops->curl(d_W_idx, d_W_var, d_V_idx, d_V_var, NULL, current_time);

            // Interpolate w(n+1/2) to W(n+1/2).
            d_l_data_manager->interp(d_W_idx, W_half_data, X_half_data, d_cscheds["W->W::S->S::CONSERVATIVE_COARSEN"], d_rscheds["W->W::S->S::GHOST_FILL"], current_time);
            resetAnchorPointValues(W_half_data, coarsest_ln, finest_ln);

            // Compute D(n+1) from D(n) and W(n+1/2), and compute D(n+1/2) from
            // D(n) and D(n+1).
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                if (d_l_data_manager->levelContainsLagrangianData(ln))
                {
                    blitz::Array<double,2>&      D_array = *     D_data[ln]->getLocalFormVecArray();
                    blitz::Array<double,2>& D_half_array = *D_half_data[ln]->getLocalFormVecArray();
                    blitz::Array<double,2>&  D_new_array = * D_new_data[ln]->getLocalFormVecArray();
                    blitz::Array<double,2>& W_half_array = *W_half_data[ln]->getLocalFormVecArray();
                    const int n_local = d_l_data_manager->getNumberOfLocalNodes(ln);
                    for (int l = 0; l < n_local; ++l)
                    {
                        std::vector<double*> D(3);
                        for (int alpha = 0; alpha < 3; ++alpha)
                        {
                            D[alpha] = &D_array(l,3*alpha);
                        }

                        std::vector<double*> D_new(3);
                        for (int alpha = 0; alpha < 3; ++alpha)
                        {
                            D_new[alpha] = &D_new_array(l,3*alpha);
                        }

                        double* W_half = &W_half_array(l,0);

                        blitz::Array<double,1> e(W_half, blitz::shape(3), blitz::duplicateData);
                        blitz::firstIndex i;
                        const double norm_e = sqrt(blitz::sum(e(i)*e(i)));

                        if (norm_e > std::numeric_limits<double>::epsilon())
                        {
                            e /= norm_e;
                            const double theta = 0.5*dt*norm_e;

                            const double c_theta = cos(theta);
                            const double s_theta = sin(theta);
                            blitz::Array<double,2> R(3,3);
                            R(0,0) =  c_theta     ;  R(0,1) = -s_theta*e(2);  R(0,2) =  s_theta*e(1);
                            R(1,0) =  s_theta*e(2);  R(1,1) =  c_theta     ;  R(1,2) = -s_theta*e(0);
                            R(2,0) = -s_theta*e(1);  R(2,1) =  s_theta*e(0);  R(2,2) =  c_theta     ;
                            for (int j = 0; j < 3; ++j)
                            {
                                for (int i = 0; i < 3; ++i)
                                {
                                    R(i,j) += (1.0-c_theta)*e(i)*e(j);
                                }
                            }

                            for (int alpha = 0; alpha < 3; ++alpha)
                            {
                                blitz::Array<double,1> D_arr(D[alpha], blitz::shape(3), blitz::neverDeleteData);
                                blitz::Array<double,1> D_new_arr(D_new[alpha], blitz::shape(3), blitz::neverDeleteData);

                                blitz::firstIndex i;
                                blitz::secondIndex j;
                                D_new_arr = blitz::sum(R(i,j)*D_arr(j),j);
                            }
                        }
                        else
                        {
                            for (int alpha = 0; alpha < 3; ++alpha)
                            {
                                blitz::Array<double,1> D_arr(D[alpha], blitz::shape(3), blitz::neverDeleteData);
                                blitz::Array<double,1> D_new_arr(D_new[alpha], blitz::shape(3), blitz::neverDeleteData);
                                D_new_arr = D_arr;
                            }
                        }

                        std::vector<double*> D_half(3);
                        for (int alpha = 0; alpha < 3; ++alpha)
                        {
                            D_half[alpha] = &D_half_array(l,3*alpha);
                        }

                        interpolate_directors(D_half, D, D_new);
                    }
                    D_data     [ln]->restoreArrays();
                    D_half_data[ln]->restoreArrays();
                    D_new_data [ln]->restoreArrays();
                    W_half_data[ln]->restoreArrays();
                }
            }
        }

        // Compute the pressure at the updated locations of any distributed
        // internal fluid sources or sinks.
        if (!d_source_strategy.isNull())
        {
            computeSourcePressures(coarsest_ln, finest_ln, current_time+0.5*dt, X_half_data);
        }
    }
    d_ins_hier_integrator->integrateHierarchy_finalize(current_time, new_time);

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

    // Reset X to equal X_new.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_vec = X_data[ln]->getVec();
            Vec X_new_vec = X_new_data[ln]->getVec();
            int ierr = VecCopy(X_new_vec, X_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Reset X_mark to equal X_mark_new.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<LMarkerSetData> mark_data = patch->getPatchData(d_mark_current_idx);
            Pointer<LMarkerSetData> mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);
            mark_data->copy(*mark_scratch_data);
        }
        level->deallocatePatchData(d_mark_scratch_idx);
    }

    if (d_using_pIB_method)
    {
        // Reset Y to equal Y_new.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                Vec Y_vec = Y_data[ln]->getVec();
                Vec Y_new_vec = Y_new_data[ln]->getVec();
                int ierr = VecCopy(Y_new_vec, Y_vec);  IBTK_CHKERRQ(ierr);
            }
        }

        // Reset dY_dt to equal dY_dt_new.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                Vec dY_dt_vec = dY_dt_data[ln]->getVec();
                Vec dY_dt_new_vec = dY_dt_new_data[ln]->getVec();
                int ierr = VecCopy(dY_dt_new_vec, dY_dt_vec);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    if (d_using_orthonormal_directors)
    {
        // Reset D to equal D_new.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                Vec D_vec = D_data[ln]->getVec();
                Vec D_new_vec = D_new_data[ln]->getVec();
                int ierr = VecCopy(D_new_vec, D_vec);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_V_idx);
        level->deallocatePatchData(d_F_idx);
        if (d_using_orthonormal_directors)
        {
            level->deallocatePatchData(d_W_idx);
            level->deallocatePatchData(d_N_idx);
        }
        if (!d_source_strategy.isNull())
        {
            level->deallocatePatchData(d_Q_idx);
        }
    }

    // Synchronize all data.
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): synchronizing data\n";
    synchronizeHierarchy();

    // Reset all time dependent data.
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierData(new_time);

    // Determine the current flow CFL number.
    double cfl_max = 0.0;
    PatchSideDataOpsReal<NDIM,double> patch_ops;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double dx_min = *(std::min_element(dx,dx+NDIM));
            Pointer<SideData<NDIM,double> > U_current_data = patch->getPatchData(U_current_idx);
            const double U_max = std::max(+patch_ops.max(U_current_data, patch_box),
                                          -patch_ops.min(U_current_data, patch_box));
            cfl_max = std::max(cfl_max, U_max*dt/dx_min);
        }
    }
    cfl_max = SAMRAI_MPI::maxReduction(cfl_max);
    d_regrid_cfl_estimate += cfl_max;

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

void
IBStaggeredHierarchyIntegrator::postProcessData()
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
    d_hier_sc_data_ops->copyData(d_V_idx, U_current_idx);
    d_l_data_manager->interp(d_V_idx, U_data, X_data, d_cscheds["V->V::S->S::CONSERVATIVE_COARSEN"], d_rscheds["V->V::S->S::GHOST_FILL"], d_integrator_time);
    resetAnchorPointValues(U_data, coarsest_ln, finest_ln);

    // Compute F(n) = F(X(n),U(n),n), the Lagrangian force corresponding to
    // configuration X(n) and velocity U(n) at time t_{n}.
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
            d_force_strategy->computeLagrangianForce(F_data[ln], X_data[ln], U_data[ln], d_hierarchy, ln, current_time, d_l_data_manager);
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
IBStaggeredHierarchyIntegrator::atRegridPoint() const
{
    return d_regrid_cfl_interval > 0.0 ? (d_regrid_cfl_estimate >= d_regrid_cfl_interval) : (d_regrid_interval == 0 ? false : (d_integrator_step % d_regrid_interval == 0));
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
    return d_integrator_step < d_max_integrator_steps;
}// stepsRemaining

const Pointer<PatchHierarchy<NDIM> >
IBStaggeredHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

Pointer<GriddingAlgorithm<NDIM> >
IBStaggeredHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

LDataManager*
IBStaggeredHierarchyIntegrator::getLDataManager() const
{
    return d_l_data_manager;
}// getLDataManager

Pointer<IBInstrumentPanel>
IBStaggeredHierarchyIntegrator::getIBInstrumentPanel() const
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
///  allow the IBStaggeredHierarchyIntegrator to provide data management for a time
///  integrator which making use of this class.
///

void
IBStaggeredHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    // Determine the current range of hierarchy levels.
    const int coarsest_ln_before_regrid = 0;
    const int finest_ln_before_regrid = d_hierarchy->getFinestLevelNumber();

    // Update the marker data.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): resetting markers particles.\n";
    LMarkerUtilities::collectMarkersOnPatchHierarchy(d_mark_current_idx, d_hierarchy);

    // Update the workload pre-regridding.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_l_data_manager->updateWorkloadData(coarsest_ln_before_regrid,finest_ln_before_regrid);

    // Before regridding, begin Lagrangian data movement.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): starting Lagrangian data movement.\n";
    d_l_data_manager->beginDataRedistribution(coarsest_ln_before_regrid,finest_ln_before_regrid);

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): calling INSStaggeredHierarchyIntegrator::regridHierarchy().\n";
    d_ins_hier_integrator->regridHierarchy();

    // After regridding, finish Lagrangian data movement.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): finishing Lagrangian data movement.\n";
    d_l_data_manager->endDataRedistribution(coarsest_ln_before_regrid,finest_ln_before_regrid);

    // Update the workload post-regridding.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_l_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());

    // Prune duplicate markers following regridding.
    LMarkerUtilities::pruneDuplicateMarkers(d_mark_current_idx, d_hierarchy);

    // Indicate that the force and (optional) source strategies and
    // post-processor need to be re-initialized.
    d_force_strategy_needs_init  = true;
    d_source_strategy_needs_init = true;
    d_post_processor_needs_init  = true;

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

    // Reset the regrid CFL estimate.
    d_regrid_cfl_estimate = 0.0;

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

void
IBStaggeredHierarchyIntegrator::synchronizeHierarchy()
{
    t_synchronize_hierarchy->start();

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->synchronizeHierarchy();

    t_synchronize_hierarchy->stop();
    return;
}// synchronizeHierarchy

void
IBStaggeredHierarchyIntegrator::synchronizeNewLevels(
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
    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->synchronizeNewLevels(hierarchy, coarsest_level, finest_level, sync_time, initial_time);

    t_synchronize_new_levels->stop();
    return;
}// synchronizeNewLevels

void
IBStaggeredHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    t_reset_time_dependent_data->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->resetTimeDependentHierData(new_time);

    // Reset the time of the current data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->setTime(new_time, d_mark_current_idx);
    }

    t_reset_time_dependent_data->stop();
    return;
}// resetTimeDependentHierData

void
IBStaggeredHierarchyIntegrator::resetHierDataToPreadvanceState()
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
///  StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
IBStaggeredHierarchyIntegrator::initializeLevelData(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > old_level,
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

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    // We use the LDataManager to handle as much unstructured data management as
    // possible.
    d_l_data_manager->setPatchHierarchy(hierarchy);
    d_l_data_manager->resetLevels(0,hierarchy->getFinestLevelNumber());
    d_l_data_manager->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    // Allocate storage needed to initialize the level and fill data from
    // coarser levels in AMR hierarchy, if any.
    //
    // Since time gets set when we allocate data, re-stamp it to current time if
    // we don't need to allocate.
    if (allocate_data)
    {
        level->allocatePatchData(d_mark_current_idx, init_data_time);
    }
    else
    {
        level->setTime(init_data_time, d_mark_current_idx);
    }

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBStaggeredHierarchyIntegrator::resetHierarchyConfiguration(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
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

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // We use the LDataManager to handle as much unstructured data management as
    // possible.
    d_l_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);
    d_hier_sc_data_ops->resetLevels(0, finest_hier_level);

    // If we have added or removed a level, resize the anchor point vectors.
    d_anchor_point_local_idxs.clear();
    d_anchor_point_local_idxs.resize(finest_hier_level+1);

    // If we have added or removed a level, resize the source/sink data vectors.
    d_X_src.resize(finest_hier_level+1);
    d_r_src.resize(finest_hier_level+1);
    d_P_src.resize(finest_hier_level+1);
    d_Q_src.resize(finest_hier_level+1);
    d_n_src.resize(finest_hier_level+1,0);

    // If we have added or removed a level, resize the schedule vectors.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin(); it != d_ralgs.end(); ++it)
    {
        d_rscheds[it->first].resize(finest_hier_level+1);
    }

    for (RefineAlgMap::const_iterator it = d_palgs.begin(); it != d_palgs.end(); ++it)
    {
        d_pscheds[it->first].resize(finest_hier_level+1);
    }

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
            d_rscheds[it->first][ln] = it->second->createSchedule(level, ln-1, hierarchy, d_rstrategies[it->first]);
        }
    }

    // (Re)build generic prolongation communication schedules.  These are set only for levels
    // >= 1.
    for (RefineAlgMap::const_iterator it = d_palgs.begin(); it != d_palgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            d_pscheds[it->first][ln] = it->second->createSchedule(level, Pointer<PatchLevel<NDIM> >(), ln-1, d_hierarchy, d_pstrategies[it->first]);
        }
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin(); it != d_calgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(ln-1);
            d_cscheds[it->first][ln] = it->second->createSchedule(coarser_level, level, d_cstrategies[it->first]);
        }
    }

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
IBStaggeredHierarchyIntegrator::applyGradientDetector(
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

    // Tag cells for refinement according to the criteria specified by the
    // INSStaggeredHierarchyIntegrator.
    d_ins_hier_integrator->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);

    // Tag cells which contain Lagrangian nodes.
    d_l_data_manager->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);

    // Tag cells for refinement where the Cartesian source/sink strength is
    // nonzero.
    if (!d_source_strategy.isNull() && !initial_time && hierarchy->finerLevelExists(level_number))
    {
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        if (!grid_geom->getDomainIsSingleBox()) TBOX_ERROR("physical domain must be a single box...\n");

        const Index<NDIM>& lower = grid_geom->getPhysicalDomain()[0].lower();
        const Index<NDIM>& upper = grid_geom->getPhysicalDomain()[0].upper();
        const double* const xLower = grid_geom->getXLower();
        const double* const xUpper = grid_geom->getXUpper();
        const double* const dx = grid_geom->getDx();

        const int finer_level_number = level_number+1;
        Pointer<PatchLevel<NDIM> > finer_level = hierarchy->getPatchLevel(finer_level_number);
        for (int n = 0; n < d_n_src[finer_level_number]; ++n)
        {
            blitz::TinyVector<double,NDIM> dx_finer;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dx_finer[d] = dx[d]/static_cast<double>(finer_level->getRatio()(d));
            }

            // The source radius must be an integer multiple of the grid
            // spacing.
            blitz::TinyVector<double,NDIM> r;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                r[d] = floor(d_r_src[finer_level_number][n]/dx_finer[d])*dx_finer[d];
                r[d] = std::max(r[d],2.0*dx_finer[d]);
            }

            // Determine the approximate source stencil box.
            const Index<NDIM> i_center = IndexUtilities::getCellIndex(d_X_src[finer_level_number][n], xLower, xUpper, dx_finer.data(), lower, upper);
            Box<NDIM> stencil_box(i_center,i_center);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                stencil_box.grow(d, static_cast<int>(ceil(r[d]/dx_finer[d])));
            }
            const Box<NDIM> coarsened_stencil_box = Box<NDIM>::coarsen(stencil_box, finer_level->getRatioToCoarserLevel());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
                tags_data->fillAll(1, coarsened_stencil_box);
            }
        }
    }

    t_apply_gradient_detector->stop();
    return;
}// applyGradientDetector

///
///  The following routines:
///
///      getLMarkerSetVar()
///
///  allows access to the various state variables maintained by the integrator.
///

Pointer<LMarkerSetVariable>
IBStaggeredHierarchyIntegrator::getLMarkerSetVar() const
{
    return d_mark_var;
}// getLMarkerSetVar

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
/// INSStaggeredHierarchyIntegrator object.
///

Pointer<VariableContext>
IBStaggeredHierarchyIntegrator::getCurrentContext() const
{
    return d_ins_hier_integrator->getCurrentContext();
}// getCurrentContext

Pointer<VariableContext>
IBStaggeredHierarchyIntegrator::getNewContext() const
{
    return d_ins_hier_integrator->getNewContext();
}// getNewContext

Pointer<VariableContext>
IBStaggeredHierarchyIntegrator::getScratchContext() const
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
IBStaggeredHierarchyIntegrator::putToDatabase(
    Pointer<Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("IB_STAGGERED_HIERARCHY_INTEGRATOR_VERSION",
                   IB_STAGGERED_HIERARCHY_INTEGRATOR_VERSION);

    db->putString("d_interp_delta_fcn", d_interp_delta_fcn);
    db->putString("d_spread_delta_fcn", d_spread_delta_fcn);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);

    const std::vector<std::string>& instrument_names = IBInstrumentationSpec::getInstrumentNames();
    if (!instrument_names.empty())
    {
        db->putInteger("instrument_names_sz", instrument_names.size());
        db->putStringArray("instrument_names", &instrument_names[0], instrument_names.size());
    }
    db->putInteger("d_total_flow_volume_sz", d_total_flow_volume.size());
    if (!d_total_flow_volume.empty())
    {
        db->putDoubleArray("d_total_flow_volume", &d_total_flow_volume[0], d_total_flow_volume.size());
    }

    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    db->putInteger("finest_hier_level", finest_hier_level);
    db->putIntegerArray("d_n_src", &d_n_src[0], finest_hier_level+1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        for (int n = 0; n < d_n_src[ln]; ++n)
        {
            std::ostringstream id_stream;
            id_stream << ln << "_" << n;
            const std::string id_string = id_stream.str();
            db->putDoubleArray("d_X_src_"+id_string, &d_X_src[ln][n][0], NDIM);
            db->putDouble("d_r_src_"+id_string, d_r_src[ln][n]);
            db->putDouble("d_P_src_"+id_string, d_P_src[ln][n]);
            db->putDouble("d_Q_src_"+id_string, d_Q_src[ln][n]);
        }
    }

    db->putBool("d_using_pIB_method", d_using_pIB_method);
    if (d_using_pIB_method)
    {
        db->putDoubleArray("d_gravitational_acceleration", &d_gravitational_acceleration[0], NDIM);
    }

    db->putBool("d_using_orthonormal_directors", d_using_orthonormal_directors);

    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);

    db->putInteger("d_num_cycles", d_num_cycles);

    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putDouble("d_regrid_cfl_interval", d_regrid_cfl_interval);
    db->putDouble("d_regrid_cfl_estimate", d_regrid_cfl_estimate);

    db->putDouble("d_old_dt", d_old_dt);
    db->putDouble("d_integrator_time", d_integrator_time);
    db->putInteger("d_integrator_step", d_integrator_step);

    db->putDouble("d_dt_max", d_dt_max);
    db->putDouble("d_dt_max_time_max", d_dt_max_time_max);
    db->putDouble("d_dt_max_time_min", d_dt_max_time_min);

    t_put_to_database->stop();
    return;
}// putToDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBStaggeredHierarchyIntegrator::resetLagrangianForceStrategy(
    const double init_data_time,
    const bool initial_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            d_force_strategy->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
        }
    }
    return;
}// resetLagrangianForceStrategy

void
IBStaggeredHierarchyIntegrator::resetLagrangianSourceStrategy(
    const double init_data_time,
    const bool initial_time)
{
    if (!d_source_strategy.isNull())
    {
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                d_source_strategy->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
            }
        }
    }
    return;
}// resetLagrangianSourceStrategy

void
IBStaggeredHierarchyIntegrator::resetPostProcessor(
    const double init_data_time,
    const bool initial_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            d_post_processor->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
        }
    }
    return;
}// resetPostProcessor

void
IBStaggeredHierarchyIntegrator::updateIBInstrumentationData(
    const int timestep_num,
    const double data_time)
{
    if (!d_instrument_panel->isInstrumented()) return;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Compute the positions of the flow meter nets.
    d_instrument_panel->initializeHierarchyDependentData(d_hierarchy, d_l_data_manager, timestep_num, data_time);

    // Compute the flow rates and pressures.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int U_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getScratchContext());
    const int P_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVar(), d_ins_hier_integrator->getScratchContext());

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

    d_instrument_panel->readInstrumentData(U_scratch_idx, P_scratch_idx, d_hierarchy, d_l_data_manager, timestep_num, data_time);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (deallocate_U_scratch_data[ln]) level->deallocatePatchData(U_scratch_idx);
        if (deallocate_P_scratch_data[ln]) level->deallocatePatchData(P_scratch_idx);
    }
    return;
}// updateIBInstrumentationData

void
IBStaggeredHierarchyIntegrator::resetAnchorPointValues(
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
            for (std::set<int>::const_iterator cit = d_anchor_point_local_idxs[ln].begin(); cit != d_anchor_point_local_idxs[ln].end(); ++cit)
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
IBStaggeredHierarchyIntegrator::computeSourceStrengths(
    const int coarsest_level,
    const int finest_level,
    const double data_time,
    const std::vector<Pointer<LData> >& X_data)
{
    if (d_source_strategy.isNull()) return;

    // Reset the values of Q_src.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        d_Q_src[ln] = std::vector<double>(d_n_src[ln],0.0);
    }

    // Get the present source locations.
    //
    // IMPORTANT NOTE: Here, we require that each MPI process is assigned the
    // same source locations and radii.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        if (d_n_src[ln] > 0)
        {
            if (d_do_log) plog << d_object_name << "::advanceHierarchy(): computing source locations on level number " << ln << "\n";
            d_source_strategy->getSourceLocations(d_X_src[ln], d_r_src[ln], X_data[ln], d_hierarchy, ln, data_time, d_l_data_manager);
        }
    }

    // Compute the source strengths for each of the sources/sinks.
    //
    // IMPORTANT NOTE: Here, we require that each MPI process is assigned the
    // same source strengths.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        if (d_n_src[ln] > 0)
        {
            if (d_do_log) plog << d_object_name << "::advanceHierarchy(): computing fluid source strengths on level number " << ln << "\n";
            d_source_strategy->computeSourceStrengths(d_Q_src[ln], d_hierarchy, ln, data_time, d_l_data_manager);
        }
    }

    // Spread the sources/sinks onto the Cartesian grid.
    d_hier_cc_data_ops->setToScalar(d_Q_idx, 0.0);
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        if (d_n_src[ln] > 0)
        {
            TBOX_ASSERT(ln == d_hierarchy->getFinestLevelNumber());
            if (d_do_log) plog << d_object_name << "::advanceHierarchy(): spreading fluid source strengths to the Cartesian grid on level number " << ln << "\n";
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Index<NDIM>& patch_lower = patch_box.lower();
                const Index<NDIM>& patch_upper = patch_box.upper();

                const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

                const double* const xLower = pgeom->getXLower();
                const double* const xUpper = pgeom->getXUpper();
                const double* const dx = pgeom->getDx();

                const Pointer<CellData<NDIM,double> > q_data = patch->getPatchData(d_Q_idx);
                for (int n = 0; n < d_n_src[ln]; ++n)
                {
                    // The source radius must be an integer multiple of the grid
                    // spacing.
                    blitz::TinyVector<double,NDIM> r;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        r[d] = floor(d_r_src[ln][n]/dx[d])*dx[d];
                        r[d] = std::max(r[d],2.0*dx[d]);
                    }

                    // Determine the approximate source stencil box.
                    const Index<NDIM> i_center = IndexUtilities::getCellIndex(d_X_src[ln][n], xLower, xUpper, dx, patch_lower, patch_upper);
                    Box<NDIM> stencil_box(i_center,i_center);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        stencil_box.grow(d, static_cast<int>(ceil(r[d]/dx[d])));
                    }

                    // Spread the source strength onto the Cartesian grid.
                    for (Box<NDIM>::Iterator b(patch_box*stencil_box); b; b++)
                    {
                        const Index<NDIM>& i = b();

                        double wgt = 1.0;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            const double X_center = xLower[d] + dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
                            wgt *= cos_delta(X_center - d_X_src[ln][n][d], r[d]);
                        }
                        (*q_data)(i) += d_Q_src[ln][n]*wgt;
                    }
                }
            }
        }
    }

    // Compute the net inflow into the computational domain.
    const int wgt_idx = d_ins_hier_integrator->getHierarchyMathOps()->getCellWeightPatchDescriptorIndex();
    PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
    double Q_sum = 0.0;
    double Q_max = 0.0;
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        Q_sum = std::accumulate(d_Q_src[ln].begin(), d_Q_src[ln].end(), Q_sum);
        for (unsigned int k = 0; k < d_Q_src[ln].size(); ++k)
        {
            Q_max = std::max(Q_max,std::abs(d_Q_src[ln][k]));
        }
    }
    const double q_total = d_hier_cc_data_ops->integral(d_Q_idx, wgt_idx);

    if (d_do_log)
    {
        plog << d_object_name << "::advanceHierarchy():\n";
#if (NDIM == 2)
        plog << "    Sum_{i,j} q_{i,j} h^2     = " << q_total << "\n"
             << "    Sum_{l=1,...,n_src} Q_{l} = " << Q_sum << "\n";
#endif
#if (NDIM == 3)
        plog << "    Sum_{i,j,k} q_{i,j,k} h^3 = " << q_total << "\n"
             << "    Sum_{l=1,...,n_src} Q_{l} = " << Q_sum <<  "\n";
#endif
    }

    if (std::abs(q_total-Q_sum)                     > 1.0e-12 &&
        std::abs(q_total-Q_sum)/std::max(Q_max,1.0) > 1.0e-12)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Lagrangian and Eulerian source/sink strengths are inconsistent.");
    }

    // When necessary, balance the net inflow/outflow with outflow/inflow along
    // the upper/lower boundaries of the computational domain.
    if (std::abs(q_total) > 1.0e-12)
    {
        plog << "    adding ``external'' source/sink to offset net inflow/outflow into domain.\n";
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
        const Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];
        const double* const dx_coarsest = grid_geom->getDx();

        Box<NDIM> interior_box = domain_box;
        for (unsigned int d = 0; d < NDIM-1; ++d)
        {
            interior_box.grow(d,-1);
        }

        BoxList<NDIM> bdry_boxes;
        bdry_boxes.removeIntersections(domain_box,interior_box);
        double vol = static_cast<double>(bdry_boxes.getTotalSizeOfBoxes());
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            vol *= dx_coarsest[d];
        }

        const double q_norm = -q_total/vol;
        for (int ln = coarsest_level; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            BoxList<NDIM> level_bdry_boxes(bdry_boxes);
            level_bdry_boxes.refine(level->getRatio());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM,double> > q_data = patch->getPatchData(d_Q_idx);
                for (BoxList<NDIM>::Iterator blist(level_bdry_boxes); blist; blist++)
                {
                    for (Box<NDIM>::Iterator b(blist()*patch_box); b; b++)
                    {
                        (*q_data)(b()) += q_norm;
                    }
                }
            }
        }

        if (std::abs(d_hier_cc_data_ops->integral(d_Q_idx, wgt_idx) > 1.0e-12))
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "``external'' source/sink does not correctly offset net inflow/outflow into domain.\n"
                       << "integral{q} = " << d_hier_cc_data_ops->integral(d_Q_idx, wgt_idx) << " != 0.\n");
        }
    }
    return;
}// computeSourceStrengths

void
IBStaggeredHierarchyIntegrator::computeSourcePressures(
    const int coarsest_level,
    const int finest_level,
    const double data_time,
    const std::vector<Pointer<LData> >& X_data)
{
    if (d_source_strategy.isNull()) return;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int P_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVar(), d_ins_hier_integrator->getNewContext());
    const int wgt_idx = d_ins_hier_integrator->getHierarchyMathOps()->getCellWeightPatchDescriptorIndex();

    // Compute the normalization pressure.
    double p_norm = 0.0;
    if (true)
    {
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
        const Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];

        Box<NDIM> interior_box = domain_box;
        for (unsigned int d = 0; d < NDIM-1; ++d)
        {
            interior_box.grow(d,-1);
        }

        BoxList<NDIM> bdry_boxes;
        bdry_boxes.removeIntersections(domain_box,interior_box);

        double vol = 0.0;
        for (int ln = coarsest_level; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            BoxList<NDIM> level_bdry_boxes(bdry_boxes);
            level_bdry_boxes.refine(level->getRatio());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM,double> > p_data = patch->getPatchData(P_new_idx);
                const Pointer<CellData<NDIM,double> > wgt_data = patch->getPatchData(wgt_idx);
                for (BoxList<NDIM>::Iterator blist(level_bdry_boxes); blist; blist++)
                {
                    for (Box<NDIM>::Iterator b(blist()*patch_box); b; b++)
                    {
                        const Index<NDIM>& i = b();
                        p_norm += (*p_data)(i)*(*wgt_data)(i);
                        vol += (*wgt_data)(i);
                    }
                }
            }
        }

        SAMRAI_MPI::sumReduction(&p_norm,1);
        SAMRAI_MPI::sumReduction(&vol,1);

        p_norm /= vol;
    }

    // Reset the values of P_src.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        d_P_src[ln] = std::vector<double>(d_n_src[ln],0.0);
    }

    // Get the present source locations.
    //
    // IMPORTANT NOTE: Here, we require that each MPI process is assigned the
    // same source locations and radii.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        if (d_n_src[ln] > 0)
        {
            if (d_do_log) plog << d_object_name << "::advanceHierarchy(): computing source locations on level number " << ln << "\n";

            d_source_strategy->getSourceLocations(d_X_src[ln], d_r_src[ln], X_data[ln], d_hierarchy, ln, data_time, d_l_data_manager);
        }
    }

    // Compute the mean pressure at the sources/sinks associated with each level
    // of the Cartesian grid.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        if (d_n_src[ln] > 0)
        {
            if (d_do_log) plog << d_object_name << "::advanceHierarchy(): computing source pressures on level number " << ln << "\n";

            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Index<NDIM>& patch_lower = patch_box.lower();
                const Index<NDIM>& patch_upper = patch_box.upper();

                const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

                const double* const xLower = pgeom->getXLower();
                const double* const xUpper = pgeom->getXUpper();
                const double* const dx = pgeom->getDx();

                const Pointer<CellData<NDIM,double> > p_data = patch->getPatchData(P_new_idx);
                for (int n = 0; n < d_n_src[ln]; ++n)
                {
                    // The source radius must be an integer multiple of the grid
                    // spacing.
                    blitz::TinyVector<double,NDIM> r;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        r[d] = floor(d_r_src[ln][n]/dx[d])*dx[d];
                        r[d] = std::max(r[d],2.0*dx[d]);
                    }

                    // Determine the approximate source stencil box.
                    const Index<NDIM> i_center = IndexUtilities::getCellIndex(d_X_src[ln][n], xLower, xUpper, dx, patch_lower, patch_upper);
                    Box<NDIM> stencil_box(i_center,i_center);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        stencil_box.grow(d, static_cast<int>(ceil(r[d]/dx[d])));
                    }

                    // Interpolate the pressure from the Cartesian grid.
                    for (Box<NDIM>::Iterator b(patch_box*stencil_box); b; b++)
                    {
                        const Index<NDIM>& i = b();

                        double wgt = 1.0;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            const double X_center = xLower[d] + dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
                            wgt *= cos_delta(X_center - d_X_src[ln][n][d], r[d])*dx[d];
                        }
                        d_P_src[ln][n] += (*p_data)(i)*wgt;
                    }
                }
            }

            SAMRAI_MPI::sumReduction(&d_P_src[ln][0], d_P_src[ln].size());
            std::transform(d_P_src[ln].begin(), d_P_src[ln].end(), d_P_src[ln].begin(),
                           std::bind2nd(std::plus<double>(),-p_norm));

            // Set the pressures for the Lagrangian source strategy.
            d_source_strategy->setSourcePressures(d_P_src[ln], d_hierarchy, ln, data_time, d_l_data_manager);
        }
    }
    return;
}// computeSourcePressures

void
IBStaggeredHierarchyIntegrator::getFromInput(
    Pointer<Database> db,
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
    d_regrid_cfl_interval = db->getDoubleWithDefault("regrid_cfl_interval", d_regrid_cfl_interval);

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault("dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault("dt_max_time_min", d_dt_max_time_min);

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
        d_using_pIB_method = db->getBoolWithDefault("using_pIB_method", d_using_pIB_method);
        if (d_using_pIB_method)
        {
            if (db->keyExists("gravitational_acceleration"))
            {
                db->getDoubleArray("gravitational_acceleration", &d_gravitational_acceleration[0], NDIM);
            }
            else
            {
                TBOX_WARNING(d_object_name << ":  "
                             << "Using penalty-IB method but key data `gravitational_acceleration' not found in input.");
            }
        }
        d_using_orthonormal_directors = db->getBoolWithDefault("using_orthonormal_directors", d_using_orthonormal_directors);
        if (d_using_orthonormal_directors && NDIM != 3)
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "Directors are currently only supported for three-dimensional problems." << std::endl);
        }
        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);
        d_mark_input_file_name = db->getStringWithDefault("marker_input_file_name", d_mark_input_file_name);
    }
    return;
}// getFromInput

void
IBStaggeredHierarchyIntegrator::getFromRestart()
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

    int ver = db->getInteger("IB_STAGGERED_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_STAGGERED_HIERARCHY_INTEGRATOR_VERSION)
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

    db->getIntegerArray("d_ghosts", d_ghosts, NDIM);

    if (db->keyExists("instrument_names"))
    {
        const int sz = db->getInteger("instrument_names_sz");
        std::vector<std::string> instrument_names(sz);
        db->getStringArray("instrument_names", &instrument_names[0], sz);
        IBInstrumentationSpec::setInstrumentNames(instrument_names);
    }
    const int total_flow_volume_sz = db->getInteger("d_total_flow_volume_sz");
    d_total_flow_volume.resize(total_flow_volume_sz, std::numeric_limits<double>::quiet_NaN());
    if (!d_total_flow_volume.empty())
    {
        db->getDoubleArray("d_total_flow_volume", &d_total_flow_volume[0], d_total_flow_volume.size());
    }

    const int finest_hier_level = db->getInteger("finest_hier_level");
    d_X_src.resize(finest_hier_level+1);
    d_r_src.resize(finest_hier_level+1);
    d_P_src.resize(finest_hier_level+1);
    d_Q_src.resize(finest_hier_level+1);
    d_n_src.resize(finest_hier_level+1,0);
    db->getIntegerArray("d_n_src", &d_n_src[0], finest_hier_level+1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        d_X_src[ln].resize(d_n_src[ln],blitz::TinyVector<double,NDIM>(std::numeric_limits<double>::quiet_NaN()));
        d_r_src[ln].resize(d_n_src[ln],std::numeric_limits<double>::quiet_NaN());
        d_P_src[ln].resize(d_n_src[ln],std::numeric_limits<double>::quiet_NaN());
        d_Q_src[ln].resize(d_n_src[ln],std::numeric_limits<double>::quiet_NaN());
        for (int n = 0; n < d_n_src[ln]; ++n)
        {
            std::ostringstream id_stream;
            id_stream << ln << "_" << n;
            const std::string id_string = id_stream.str();
            db->getDoubleArray("d_X_src_"+id_string, &d_X_src[ln][n][0], NDIM);
            d_r_src[ln][n] = db->getDouble("d_r_src_"+id_string);
            d_P_src[ln][n] = db->getDouble("d_P_src_"+id_string);
            d_Q_src[ln][n] = db->getDouble("d_Q_src_"+id_string);
        }
    }

    d_using_pIB_method = db->getBoolWithDefault("d_using_pIB_method",false);
    if (d_using_pIB_method)
    {
        db->getDoubleArray("d_gravitational_acceleration", &d_gravitational_acceleration[0], NDIM);
    }

    d_using_orthonormal_directors = db->getBoolWithDefault("d_using_orthonormal_directors",false);

    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    d_grow_dt = db->getDouble("d_grow_dt");
    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");

    d_num_cycles = db->getInteger("d_num_cycles");

    d_regrid_interval = db->getInteger("d_regrid_interval");
    d_regrid_cfl_interval = db->getDouble("d_regrid_cfl_interval");
    d_regrid_cfl_estimate = db->getDouble("d_regrid_cfl_estimate");

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

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBStaggeredHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
