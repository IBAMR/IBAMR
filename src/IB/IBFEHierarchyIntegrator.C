// Filename: IBFEHierarchyIntegrator.C
// Created on 27 Jul 2009 by Boyce Griffith
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
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LagMarkerUtilities.h>
#include <ibtk/libmesh_utilities.h>

// LIBMESH INCLUDES
#include <boundary_info.h>
#include <dense_vector.h>
#include <dof_map.h>
#include <explicit_system.h>
#include <fe.h>
#include <mesh_refinement.h>
#include <numeric_vector.h>
#include <petsc_vector.h>
#include <quadrature_trap.h>
#include <quadrature_gauss.h>
#include <string_to_enum.h>
using namespace libMesh;

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
static Pointer<Timer> t_initialize_hierarchy_integrator;
static Pointer<Timer> t_initialize_hierarchy;
static Pointer<Timer> t_advance_hierarchy;
static Pointer<Timer> t_advance_hierarchy_init;
static Pointer<Timer> t_advance_hierarchy_phase1;
static Pointer<Timer> t_advance_hierarchy_phase2;
static Pointer<Timer> t_advance_hierarchy_phase3;
static Pointer<Timer> t_advance_hierarchy_phase4;
static Pointer<Timer> t_advance_hierarchy_phase5;
static Pointer<Timer> t_advance_hierarchy_phase6;
static Pointer<Timer> t_advance_hierarchy_phase7;
static Pointer<Timer> t_advance_hierarchy_phase8;
static Pointer<Timer> t_advance_hierarchy_finalize;
static Pointer<Timer> t_regrid_hierarchy;
static Pointer<Timer> t_synchronize_hierarchy;
static Pointer<Timer> t_synchronize_new_levels;
static Pointer<Timer> t_reset_time_dependent_data;
static Pointer<Timer> t_reset_data_to_preadvance_state;
static Pointer<Timer> t_initialize_level_data;
static Pointer<Timer> t_reset_hierarchy_configuration;
static Pointer<Timer> t_apply_gradient_detector;
static Pointer<Timer> t_put_to_database;

// Version of IBFEHierarchyIntegrator restart file data.
static const int IB_FE_HIERARCHY_INTEGRATOR_VERSION = 1;
}

const std::string IBFEHierarchyIntegrator::       COORDS_SYSTEM_NAME = "IB coordinates system";
const std::string IBFEHierarchyIntegrator::COORD_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string IBFEHierarchyIntegrator::        FORCE_SYSTEM_NAME = "IB force system";
const std::string IBFEHierarchyIntegrator::     VELOCITY_SYSTEM_NAME = "IB velocity system";
const std::string IBFEHierarchyIntegrator::  PROJ_STRAIN_SYSTEM_NAME = "IB projected dilataional strain system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEHierarchyIntegrator::IBFEHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
    FEDataManager* fe_data_manager,
    LDataManager* lag_data_manager,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_fe_data_manager(fe_data_manager),
      d_fe_order(FIRST),
      d_fe_family(LAGRANGE),
      d_split_interior_and_bdry_forces(false),
      d_use_consistent_mass_matrix(true),
      d_use_fbar_projection(false),
      d_J_bar_fe_order(CONSTANT),
      d_J_bar_fe_family(MONOMIAL),
      d_coordinate_mapping_fcn(NULL),
      d_coordinate_mapping_fcn_ctx(NULL),
      d_PK1_stress_fcn(NULL),
      d_PK1_stress_fcn_ctx(NULL),
      d_lag_body_force_fcn(NULL),
      d_lag_body_force_fcn_ctx(NULL),
      d_lag_pressure_fcn(NULL),
      d_lag_pressure_fcn_ctx(NULL),
      d_lag_data_manager(lag_data_manager),
      d_hierarchy(hierarchy),
      d_gridding_alg(NULL),
      d_visit_writer(NULL),
      d_load_balancer(NULL),
      d_ins_hier_integrator(ins_hier_integrator),
      d_eul_body_force_fcn(NULL),
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
      d_mark_input_file_name(""),
      d_num_mark(0),
      d_mark_init_posns(),
      d_hier_sc_data_ops(),
      d_ralgs(),
      d_rstrategies(),
      d_rscheds(),
      d_calgs(),
      d_cstrategies(),
      d_cscheds(),
      d_V_var(NULL),
      d_F_var(NULL),
      d_mark_var(NULL),
      d_current(NULL),
      d_scratch(NULL),
      d_V_idx(-1),
      d_F_idx(-1),
      d_mark_current_idx(-1),
      d_mark_scratch_idx(-1)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT(!ins_hier_integrator.isNull());
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

    // Read in the marker initial positions.
    if (!from_restart) d_num_mark = LagMarkerUtilities::readMarkerPositions(d_mark_init_posns, d_mark_input_file_name, d_hierarchy);

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<SideVariable<NDIM,double> > sc_var = new SideVariable<NDIM,double>("sc_var");
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, hierarchy, true);

    // Initialize all variable contexts.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_current = var_db->getContext(d_object_name+"::CURRENT");
    d_scratch = var_db->getContext(d_object_name+"::SCRATCH");

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_initialize_hierarchy_integrator = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy            = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy               = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()");
        t_advance_hierarchy_init          = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()_init");
        t_advance_hierarchy_phase1        = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()_phase1");
        t_advance_hierarchy_phase2        = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()_phase2");
        t_advance_hierarchy_phase3        = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()_phase3");
        t_advance_hierarchy_phase4        = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()_phase4");
        t_advance_hierarchy_phase5        = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()_phase5");
        t_advance_hierarchy_phase6        = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()_phase6");
        t_advance_hierarchy_phase7        = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()_phase7");
        t_advance_hierarchy_phase8        = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()_phase8");
        t_advance_hierarchy_finalize      = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::advanceHierarchy()_finalize");
        t_regrid_hierarchy                = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::regridHierarchy()");
        t_synchronize_hierarchy           = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels          = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data       = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state  = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data           = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration   = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector         = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database                 = TimerManager::getManager()->getTimer("IBAMR::IBFEHierarchyIntegrator::putToDatabase()");
                  );
    return;
}// IBFEHierarchyIntegrator

IBFEHierarchyIntegrator::~IBFEHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
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
IBFEHierarchyIntegrator::registerInitialCoordinateMappingFunction(
    Point (*coordinate_mapping_fcn)(const Point& s, void* ctx),
    void* coordinate_mapping_fcn_ctx)
{
    d_coordinate_mapping_fcn = coordinate_mapping_fcn;
    d_coordinate_mapping_fcn_ctx = coordinate_mapping_fcn_ctx;
    return;
}// registerInitialCoordinateMappingFunction

void
IBFEHierarchyIntegrator::registerPK1StressTensorFunction(
    TensorValue<double> (*PK1_stress_fcn)(const TensorValue<double>& dX_ds, const Point& X, const Point& s, Elem* const elem, const double& time, void* ctx),
    void* PK1_stress_fcn_ctx)
{
    d_PK1_stress_fcn = PK1_stress_fcn;
    d_PK1_stress_fcn_ctx = PK1_stress_fcn_ctx;
    return;
}// registerPK1StressTensorFunction

void
IBFEHierarchyIntegrator::registerLagBodyForceFunction(
    VectorValue<double> (*lag_body_force_fcn)(const TensorValue<double>& dX_ds, const Point& X, const Point& s, Elem* const elem, const double& time, void* ctx),
    void* lag_body_force_fcn_ctx)
{
    d_lag_body_force_fcn = lag_body_force_fcn;
    d_lag_body_force_fcn_ctx = lag_body_force_fcn_ctx;
    return;
}// registerLagBodyForceFunction

void
IBFEHierarchyIntegrator::registerLagPressureFunction(
    double (*lag_pressure_fcn)(const TensorValue<double>& dX_ds, const Point& X, const Point& s, Elem* const elem, const unsigned short int side, const double& time, void* ctx),
    void* lag_pressure_fcn_ctx)
{
    d_lag_pressure_fcn = lag_pressure_fcn;
    d_lag_pressure_fcn_ctx = lag_pressure_fcn_ctx;
    return;
}// registerLagPressureFunction

void
IBFEHierarchyIntegrator::registerIBLagrangianForceStrategy(
    Pointer<IBLagrangianForceStrategy> ib_lag_force_strategy)
{
    d_ib_lag_force_strategy = ib_lag_force_strategy;
    d_ib_lag_force_strategy_needs_init = true;
    return;
}// registerIBLagrangianForceStrategy

void
IBFEHierarchyIntegrator::registerVelocityInitialConditions(
    Pointer<CartGridFunction> U_init)
{
    d_ins_hier_integrator->registerVelocityInitialConditions(U_init);
    return;
}// registerVelocityInitialConditions

void
IBFEHierarchyIntegrator::registerVelocityPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs)
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
IBFEHierarchyIntegrator::registerEulBodyForceFunction(
    Pointer<CartGridFunction> eul_body_force_fcn)
{
    d_eul_body_force_fcn = eul_body_force_fcn;
    return;
}// registerEulBodyForceFunction

void
IBFEHierarchyIntegrator::registerVisItDataWriter(
    Pointer<VisItDataWriter<NDIM> > visit_writer)
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
    Pointer<LoadBalancer<NDIM> > load_balancer)
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
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Initialize FE system data.
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();

    d_fe_data_manager->COORDINATES_SYSTEM_NAME = COORDS_SYSTEM_NAME;
    ExplicitSystem& X_system = equation_systems->add_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        std::ostringstream os;
        os << "X_" << d;
        X_system.add_variable(os.str(), d_fe_order, d_fe_family);
    }

    ExplicitSystem& X_mapping_system = equation_systems->add_system<ExplicitSystem>(COORD_MAPPING_SYSTEM_NAME);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        std::ostringstream os;
        os << "dX_" << d;
        X_mapping_system.add_variable(os.str(), d_fe_order, d_fe_family);
    }

    ExplicitSystem& F_system = equation_systems->add_system<ExplicitSystem>(FORCE_SYSTEM_NAME);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        std::ostringstream os;
        os << "F_" << d;
        F_system.add_variable(os.str(), d_fe_order, d_fe_family);
    }

    ExplicitSystem& U_system = equation_systems->add_system<ExplicitSystem>(VELOCITY_SYSTEM_NAME);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        std::ostringstream os;
        os << "U_" << d;
        U_system.add_variable(os.str(), d_fe_order, d_fe_family);
    }

    if (d_use_fbar_projection)
    {
        ExplicitSystem& J_bar_system = equation_systems->add_system<ExplicitSystem>(PROJ_STRAIN_SYSTEM_NAME);
        J_bar_system.add_variable("J_bar", d_J_bar_fe_order, d_J_bar_fe_family);
    }

    // Initialize all variables.
    const IntVector<NDIM> ghosts = d_fe_data_manager->getGhostCellWidth();
    const IntVector<NDIM> no_ghosts = 0;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    d_V_var = new SideVariable<NDIM,double>(d_object_name+"::V");
    d_V_idx = var_db->registerVariableAndContext(d_V_var, d_scratch, ghosts);

    d_F_var = new SideVariable<NDIM,double>(d_object_name+"::F");
    d_F_idx = var_db->registerVariableAndContext(d_F_var, d_scratch, no_ghosts);

    d_mark_var = new IndexVariable<NDIM,LagMarker,CellGeometry<NDIM> >(d_object_name+"::mark");
    d_mark_current_idx = var_db->registerVariableAndContext(d_mark_var, getCurrentContext(), ghosts);
    d_mark_scratch_idx = var_db->registerVariableAndContext(d_mark_var, getScratchContext(), ghosts);
    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_mark_current_idx);
    }

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
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineOperator<NDIM> > refine_operator;
    Pointer<CoarsenOperator<NDIM> > coarsen_operator;

    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getCurrentContext());

    d_ralgs["U->V::C->S::GHOST_FILL"] = new RefineAlgorithm<NDIM>();
    refine_operator = NULL;
    d_ralgs["U->V::C->S::GHOST_FILL"]->registerRefine(d_V_idx, U_current_idx, d_V_idx, refine_operator);

    d_ralgs["V->V::S->S::GHOST_FILL"] = new RefineAlgorithm<NDIM>();
    refine_operator = NULL;
    d_ralgs["V->V::S->S::GHOST_FILL"]->registerRefine(d_V_idx, d_V_idx, d_V_idx, refine_operator);

    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"] = new CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(d_V_var, "CONSERVATIVE_COARSEN");
    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"]->registerCoarsen(U_current_idx, U_current_idx, coarsen_operator);

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
IBFEHierarchyIntegrator::initializeHierarchy()
{
    t_initialize_hierarchy->start();

    // Initialize the data structures and state variables for the FE equation
    // systems.
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    equation_systems->init();
    initializeCoordinates();
    updateCoordinateMapping();

    // Set up boundary conditions.  Specifically, add appropriate boundary IDs
    // to the BoundaryInfo object associated with the mesh, and add DOF
    // constraints for the nodal forces and velocities.
    System& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
    X_system.assemble_before_solve = false;
    X_system.assemble();

    System& X_mapping_system = equation_systems->get_system<System>(COORD_MAPPING_SYSTEM_NAME);
    X_mapping_system.assemble_before_solve = false;
    X_mapping_system.assemble();

    ExplicitSystem& F_system = equation_systems->get_system<ExplicitSystem>(   FORCE_SYSTEM_NAME);
    ExplicitSystem& U_system = equation_systems->get_system<ExplicitSystem>(VELOCITY_SYSTEM_NAME);

    const MeshBase& mesh = equation_systems->get_mesh();
    DofMap& F_dof_map = F_system.get_dof_map();
    DofMap& U_dof_map = U_system.get_dof_map();
    const unsigned int F_sys_num = F_system.number();
    const unsigned int U_sys_num = U_system.number();

    MeshBase::const_element_iterator       el_it  = mesh.elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.elements_end();
    for ( ; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        for (unsigned int side = 0; side < elem->n_sides(); ++side)
        {
            const bool at_mesh_bdry = elem->neighbor(side) == NULL;
            if (at_mesh_bdry)
            {
                const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
                const bool at_dirichlet_bdry = std::find(bdry_ids.begin(), bdry_ids.end(), FEDataManager::DIRICHLET_BDRY_ID) != bdry_ids.end();
                if (at_dirichlet_bdry)
                {
                    for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                    {
                        if (elem->is_node_on_side(n, side))
                        {
                            Node* node = elem->get_node(n);
                            mesh.boundary_info->add_node(node, FEDataManager::DIRICHLET_BDRY_ID);

                            if (node->n_dofs(F_sys_num) > 0)
                            {
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    const int F_dof_index = node->dof_number(F_sys_num,d,0);
                                    DofConstraintRow F_constraint_row;
                                    F_constraint_row[F_dof_index] = 1.0;
                                    F_dof_map.add_constraint_row(F_dof_index, F_constraint_row, false);
                                }
                            }

                            if (node->n_dofs(U_sys_num) > 0)
                            {
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    const int U_dof_index = node->dof_number(U_sys_num,d,0);
                                    DofConstraintRow U_constraint_row;
                                    U_constraint_row[U_dof_index] = 1.0;
                                    U_dof_map.add_constraint_row(U_dof_index, U_constraint_row, false);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Use the INSStaggeredHierarchyIntegrator to initialize the patch
    // hierarchy.
    double dt_next = d_ins_hier_integrator->initializeHierarchy();
    if (d_integrator_time >= d_dt_max_time_min && d_integrator_time <= d_dt_max_time_max)
    {
        dt_next = std::min(dt_next, d_dt_max);
    }

    // Reset the Lagrangian data manager.
    if (d_lag_data_manager != NULL)
    {
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        d_lag_data_manager->beginDataRedistribution();
        d_lag_data_manager->endDataRedistribution();
        d_lag_data_manager->updateWorkloadData(coarsest_ln, finest_ln);
    }

    // Initialize the FE data manager.
    d_fe_data_manager->reinitElementMappings();

    // Compute cached FE data.
    computeCachedProjectedDilatationalStrainFEData();
    computeCachedInteriorForceDensityFEData();
    computeCachedTransmissionForceDensityFEData();

    // Prune duplicate markers following initialization.
    LagMarkerUtilities::pruneDuplicateMarkers(d_mark_current_idx, d_hierarchy);

    // Indicate that the force strategy needs to be re-initialized.
    d_ib_lag_force_strategy_needs_init = true;

    t_initialize_hierarchy->stop();
    return dt_next;
}// initializeHierarchy

double
IBFEHierarchyIntegrator::advanceHierarchy(
    const double dt)
{
    t_advance_hierarchy->start();
    t_advance_hierarchy_init->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
#endif

    PetscErrorCode ierr;

    const double current_time = d_integrator_time;
    const double new_time     = d_integrator_time+dt;
    const bool initial_time   = MathUtilities<double>::equalEps(current_time,d_start_time);

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
    d_eulerian_force_fcn->registerBodyForceSpecification(d_eul_body_force_fcn);
    d_eulerian_force_fcn->setTimeInterval(current_time, new_time);
    if (!d_ib_lag_force_strategy.isNull()) d_ib_lag_force_strategy->setTimeInterval(current_time, new_time);

    // (Re)initialize the force strategy object.
    if (d_ib_lag_force_strategy_needs_init)
    {
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            if (d_lag_data_manager != NULL && d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                d_ib_lag_force_strategy->initializeLevelData(d_hierarchy, ln, current_time, initial_time, d_lag_data_manager);
            }
        }
        d_ib_lag_force_strategy_needs_init = false;
    }

    // Extract the FE vectors.
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    System& X_system = equation_systems->get_system<System>(  COORDS_SYSTEM_NAME);
    System& F_system = equation_systems->get_system<System>(   FORCE_SYSTEM_NAME);
    System& U_system = equation_systems->get_system<System>(VELOCITY_SYSTEM_NAME);

    NumericVector<double>& X_current = *(X_system.solution);
    X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_current);

    AutoPtr<NumericVector<double> > X_new_ptr = X_current.clone();
    NumericVector<double>& X_new = *X_new_ptr;

    NumericVector<double>& X_half = *(X_system.current_local_solution);
    NumericVector<double>& U_half = *(U_system.              solution);
    NumericVector<double>& F_half = *(F_system.              solution);

    NumericVector<double>* X_half_IB_ghost_ptr = d_fe_data_manager->getGhostedCoordsVector();
    NumericVector<double>& X_half_IB_ghost = *X_half_IB_ghost_ptr;

    NumericVector<double>* F_half_IB_ghost_ptr = d_fe_data_manager->getGhostedSolutionVector(FORCE_SYSTEM_NAME);
    NumericVector<double>& F_half_IB_ghost = *F_half_IB_ghost_ptr;

    NumericVector<double>* J_bar_half = NULL;
    NumericVector<double>* J_bar_half_IB_ghost = NULL;
    if (d_use_fbar_projection)
    {
        System& J_bar_system = equation_systems->get_system<System>(PROJ_STRAIN_SYSTEM_NAME);
        J_bar_half = J_bar_system.solution.get();
        J_bar_half_IB_ghost = d_fe_data_manager->getGhostedSolutionVector(PROJ_STRAIN_SYSTEM_NAME);
    }

    // Initialize the various LNodeLevelData objects on each level of the patch hierarchy.
    std::vector<Pointer<LNodeLevelData> > X_current_data(finest_ln+1);
    std::vector<Pointer<LNodeLevelData> > X_new_data(finest_ln+1);
    std::vector<Pointer<LNodeLevelData> > X_half_data(finest_ln+1);
    std::vector<Pointer<LNodeLevelData> > U_half_data(finest_ln+1);
    std::vector<Pointer<LNodeLevelData> > F_half_data(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager != NULL && d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_current_data[ln] = d_lag_data_manager->getLNodeLevelData(LDataManager::POSN_DATA_NAME,ln);
            X_new_data    [ln] = d_lag_data_manager->createLNodeLevelData("X_new" ,ln,NDIM);
            X_half_data   [ln] = d_lag_data_manager->createLNodeLevelData("X_half",ln,NDIM);
            U_half_data   [ln] = d_lag_data_manager->createLNodeLevelData("U_half",ln,NDIM);
            F_half_data   [ln] = d_lag_data_manager->createLNodeLevelData("F_half",ln,NDIM);

            X_current_data[ln]->restoreLocalFormVec();
            X_new_data    [ln]->restoreLocalFormVec();
            X_half_data   [ln]->restoreLocalFormVec();
            U_half_data   [ln]->restoreLocalFormVec();
            F_half_data   [ln]->restoreLocalFormVec();
        }
    }

    // Get patch data descriptors for the current and new velocity data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getCurrentContext());
    const int U_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getNewContext());

    // Allocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
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
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_rscheds["U->V::C->S::GHOST_FILL"][ln]->fillData(current_time);
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

    // Initialize X_mark(n+1) to equal X_mark(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_mark_scratch_idx, d_integrator_time);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<IndexData<NDIM,LagMarker,CellGeometry<NDIM> > > mark_data         = patch->getPatchData(d_mark_current_idx);
            Pointer<IndexData<NDIM,LagMarker,CellGeometry<NDIM> > > mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);
            mark_scratch_data->copy(*mark_data);
        }
    }

    t_advance_hierarchy_init->stop();

    // Perform one or more cycles to compute the updated configuration of the
    // coupled fluid-structure system.
    d_ins_hier_integrator->integrateHierarchy_initialize(current_time, new_time);
    for (int cycle = 0; cycle < d_num_cycles; ++cycle)
    {
        // Set X(n+1/2) = 0.5*(X(n) + X(n+1)).
        t_advance_hierarchy_phase1->start();
        ierr = VecAXPBYPCZ(dynamic_cast<PetscVector<double>*>(&X_half)->vec(),
                           0.5, 0.5, 0.0,
                           dynamic_cast<PetscVector<double>*>(&X_current)->vec(),
                           dynamic_cast<PetscVector<double>*>(&X_new)->vec()); IBTK_CHKERRQ(ierr);
        X_half.close();
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
        t_advance_hierarchy_phase1->stop();

        // Project the deformation, if necessary.
        t_advance_hierarchy_phase2->start();
        if (d_use_fbar_projection) computeProjectedDilatationalStrain(*J_bar_half, X_half);
        t_advance_hierarchy_phase2->stop();

        // Compute F(n+1/2) = F(X(n+1/2),t(n+1/2)).
        t_advance_hierarchy_phase3->start();
        computeInteriorForceDensity(F_half, X_half, J_bar_half, current_time+0.5*dt);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_lag_data_manager != NULL && d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                Vec F_half_vec = F_half_data[ln]->getGlobalVec();
                ierr = VecSet(F_half_vec, 0.0); IBTK_CHKERRQ(ierr);
                d_ib_lag_force_strategy->computeLagrangianForce(F_half_data[ln], X_half_data[ln], U_half_data[ln], d_hierarchy, ln, current_time+0.5*dt, d_lag_data_manager);
            }
        }
        t_advance_hierarchy_phase3->stop();

        // Copy data into the "IB ghosted" vectors.
        t_advance_hierarchy_phase4->start();
        ierr = VecCopy(dynamic_cast<PetscVector<double>*>(&X_half)->vec(), dynamic_cast<PetscVector<double>*>(&X_half_IB_ghost)->vec()); IBTK_CHKERRQ(ierr);
        ierr = VecCopy(dynamic_cast<PetscVector<double>*>(&F_half)->vec(), dynamic_cast<PetscVector<double>*>(&F_half_IB_ghost)->vec()); IBTK_CHKERRQ(ierr);
        if (d_use_fbar_projection)
        {
            ierr = VecCopy(dynamic_cast<PetscVector<double>*>(J_bar_half)->vec(), dynamic_cast<PetscVector<double>*>(J_bar_half_IB_ghost)->vec()); IBTK_CHKERRQ(ierr);
        }
        t_advance_hierarchy_phase4->stop();

        // Spread F(n+1/2) to f(n+1/2).
        t_advance_hierarchy_phase5->start();
        d_hier_sc_data_ops->setToScalar(d_F_idx, 0.0);
        d_fe_data_manager->spread(d_F_idx, F_half_IB_ghost, X_half_IB_ghost, FORCE_SYSTEM_NAME, true, true);
        if (d_split_interior_and_bdry_forces)
        {
            spreadTransmissionForceDensity(d_F_idx, X_half_IB_ghost, J_bar_half_IB_ghost, current_time+0.5*dt);
        }
        if (d_lag_data_manager != NULL)
        {
            d_lag_data_manager->spread(d_F_idx, F_half_data, X_half_data);
        }
        t_advance_hierarchy_phase5->stop();

        // Solve the incompressible Navier-Stokes equations.
        t_advance_hierarchy_phase6->start();
        d_ins_hier_integrator->integrateHierarchy(current_time, new_time, cycle);
        t_advance_hierarchy_phase6->stop();

        // Set u(n+1/2) = 0.5*(u(n) + u(n+1)) and interpolate u(n+1/2) to
        // U(n+1/2).
        t_advance_hierarchy_phase7->start();
        d_hier_sc_data_ops->linearSum(d_V_idx, 0.5, U_current_idx, 0.5, U_new_idx);
        d_fe_data_manager->interp(d_V_idx, U_half, X_half_IB_ghost, VELOCITY_SYSTEM_NAME, d_rscheds["V->V::S->S::GHOST_FILL"], current_time, false);
        if (d_lag_data_manager != NULL)
        {
            d_lag_data_manager->interp(d_V_idx, U_half_data, X_half_data, std::vector<Pointer<CoarsenSchedule<NDIM> > >(), std::vector<Pointer<RefineSchedule<NDIM> > >(), current_time, false);
        }
        t_advance_hierarchy_phase7->stop();

        // Set X(n+1) = X(n) + dt*U(n+1/2).
        t_advance_hierarchy_phase8->start();
        ierr = VecWAXPY(dynamic_cast<PetscVector<double>*>(&X_new)->vec(),
                        dt,
                        dynamic_cast<PetscVector<double>*>(&U_half)->vec(),
                        dynamic_cast<PetscVector<double>*>(&X_current)->vec()); IBTK_CHKERRQ(ierr);
        X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_new);
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
        LagMarkerUtilities::advectMarkers(d_mark_current_idx, d_mark_scratch_idx, d_V_idx, dt, d_fe_data_manager->getInterpWeightingFunction(), d_hierarchy);
        t_advance_hierarchy_phase8->stop();
    }
    d_ins_hier_integrator->integrateHierarchy_finalize(current_time, new_time);

    t_advance_hierarchy_finalize->start();

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

    // Update the coordinate mapping dX = X - s.
    updateCoordinateMapping();

    // Reset X_mark to equal X_mark_new.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<IndexData<NDIM,LagMarker,CellGeometry<NDIM> > > mark_data = patch->getPatchData(d_mark_current_idx);
            Pointer<IndexData<NDIM,LagMarker,CellGeometry<NDIM> > > mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);
            mark_data->copy(*mark_scratch_data);
        }
        level->deallocatePatchData(d_mark_scratch_idx);
    }

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_V_idx);
        level->deallocatePatchData(d_F_idx);
    }

    // Synchronize all data.
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): synchronizing data\n";
    synchronizeHierarchy();

    // Reset all time dependent data.
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
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

    t_advance_hierarchy_finalize->stop();
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

const Pointer<PatchHierarchy<NDIM> >
IBFEHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

Pointer<GriddingAlgorithm<NDIM> >
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

    // Update the marker data.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): resetting markers particles.\n";
    LagMarkerUtilities::collectMarkersOnPatchHierarchy(d_mark_current_idx, d_hierarchy);

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
    d_ib_lag_force_strategy_needs_init = true;

    // Reinitialize the FE data manager.
    d_fe_data_manager->reinitElementMappings();

    // Recompute cached FE data.
    computeCachedTransmissionForceDensityFEData();

    // Prune duplicate markers following regridding.
    LagMarkerUtilities::pruneDuplicateMarkers(d_mark_current_idx, d_hierarchy);

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
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
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
///  StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
IBFEHierarchyIntegrator::initializeLevelData(
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
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

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

    // Initialize marker data.
    LagMarkerUtilities::initializeMarkersOnLevel(d_mark_current_idx, d_mark_init_posns, hierarchy, level_number, initial_time, old_level);

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBFEHierarchyIntegrator::resetHierarchyConfiguration(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
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
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            d_rscheds[(*it).first][ln] = (*it).second->createSchedule(level, ln-1, hierarchy, d_rstrategies[(*it).first]);
        }
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin(); it!= d_calgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(ln-1);
            d_cscheds[(*it).first][ln] = (*it).second->createSchedule(coarser_level, level, d_cstrategies[(*it).first]);
        }
    }

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
IBFEHierarchyIntegrator::applyGradientDetector(
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
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
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
///      getLagMarkerVar()
///
///  allows access to the various state variables maintained by the integrator.
///

Pointer<IndexVariable<NDIM,LagMarker,CellGeometry<NDIM> > >
IBFEHierarchyIntegrator::getLagMarkerVar() const
{
    return d_mark_var;
}// getLagMarkerVar

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
IBFEHierarchyIntegrator::getCurrentContext() const
{
    return d_ins_hier_integrator->getCurrentContext();
}// getCurrentContext

Pointer<VariableContext>
IBFEHierarchyIntegrator::getNewContext() const
{
    return d_ins_hier_integrator->getNewContext();
}// getNewContext

Pointer<VariableContext>
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
///  Serializable abstract base class.
///

void
IBFEHierarchyIntegrator::putToDatabase(
    Pointer<Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("IB_FE_HIERARCHY_INTEGRATOR_VERSION", IB_FE_HIERARCHY_INTEGRATOR_VERSION);

    db->putBool("d_split_interior_and_bdry_forces", d_split_interior_and_bdry_forces);
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
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

    db->putString("d_mark_input_file_name", d_mark_input_file_name);
    db->putInteger("d_num_mark", d_num_mark);
    if (d_num_mark > 0)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(NDIM*d_num_mark == int(d_mark_init_posns.size()));
#endif
        db->putDoubleArray("d_mark_init_posns", &d_mark_init_posns[0], d_mark_init_posns.size());
    }

    t_put_to_database->stop();
    return;
}// putToDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBFEHierarchyIntegrator::computeProjectedDilatationalStrain(
    NumericVector<double>& J_bar,
    NumericVector<double>& X)
{
    if (!d_use_fbar_projection) return;

    // Extract the relevant systems and DOF maps.
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    System& J_bar_system = equation_systems->get_system<System>(PROJ_STRAIN_SYSTEM_NAME);
    const DofMap& J_bar_dof_map = J_bar_system.get_dof_map();

    // Setup global and elemental right-hand-side vectors.
    AutoPtr<NumericVector<double> > F = J_bar.clone();
    F->zero();
    DenseVector<double> F_e;

    // Loop over the elements to compute the right-hand side vector.
    const MeshBase& mesh = equation_systems->get_mesh();
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
    for (int e = 0; e < std::distance(el_begin,el_end); ++e)
    {
        // Lookup cached data.
        std::vector<unsigned int> J_bar_dof_indices = d_proj_strain_J_bar_dof_indices(e);
        const blitz::Array<double,2>& J_bar_phi_JxW = d_proj_strain_J_bar_phi_JxW(e);

        const std::vector<std::vector<unsigned int> >& X_dof_indices = d_proj_strain_X_dof_indices(e);
        const blitz::Array<VectorValue<double>,2>& X_dphi = d_proj_strain_X_dphi(e);

        // Compute the elemental right-hand-side entries.
        F_e.resize(J_bar_dof_indices.size());
        const int n_qp    = J_bar_phi_JxW.extent(blitz::firstDim);
        const int n_basis = J_bar_phi_JxW.extent(blitz::secondDim);
        for (int qp = 0; qp < n_qp; ++qp)
        {
            const double J = compute_coordinate_mapping_jacobian_det(qp,X,X_dphi,X_dof_indices);
            for (int k = 0; k < n_basis; ++k)
            {
                F_e(k) += J*J_bar_phi_JxW(qp,k);
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions) and
        // add the elemental contributions to the global vector.
        J_bar_dof_map.constrain_element_vector(F_e, J_bar_dof_indices);
        F->add_vector(F_e, J_bar_dof_indices);
    }

    // Solve for J_bar.
    d_fe_data_manager->computeL2Projection(J_bar, *F, PROJ_STRAIN_SYSTEM_NAME, d_use_consistent_mass_matrix);
    return;
}// computeProjectedDilatationalStrain

void
IBFEHierarchyIntegrator::computeInteriorForceDensity(
    NumericVector<double>& G,
    NumericVector<double>& X,
    NumericVector<double>* J_bar,
    const double& time)
{
    // Extract the relevant systems and DOF maps.
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    System& F_system = equation_systems->get_system<System>(FORCE_SYSTEM_NAME);
    const DofMap& F_dof_map = F_system.get_dof_map();

    // Setup global and elemental right-hand-side vectors.
    AutoPtr<NumericVector<double> > F = G.clone();
    F->zero();
    DenseVector<double> F_e[NDIM];

    // Setup data structures related to (optional) F-bar projection.
    const DofMap* J_bar_dof_map = NULL;
    if (d_use_fbar_projection)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(J_bar != NULL);
#endif
        System& J_bar_system = equation_systems->get_system<System>(PROJ_STRAIN_SYSTEM_NAME);
        J_bar_dof_map = &J_bar_system.get_dof_map();
    }

    // Loop over the elements to accumulate the interior forces at the nodes of
    // the mesh.  These are computed via
    //
    //    F_k = -int{PP(s,t) grad phi_k(s)}ds + int{PP(s,t) N(s,t) phi_k(s)}dA(s)
    //
    // This right-hand side vector is used to solve for the nodal values of the
    // interior elastic force density.
    const MeshBase& mesh = equation_systems->get_mesh();
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        const unsigned int e = std::distance(el_begin,el_it);

        // Lookup cached data.
        const blitz::Array<Point,1>& q_point = d_interior_q_point(e);

        std::vector<std::vector<unsigned int> > F_dof_indices = d_interior_F_dof_indices(e);
        const blitz::Array<double,2>& F_phi_JxW = d_interior_F_phi_JxW(e);
        const blitz::Array<VectorValue<double>,2>& F_dphi_JxW = d_interior_F_dphi_JxW(e);

        const std::vector<std::vector<unsigned int> >& X_dof_indices = d_interior_X_dof_indices(e);
        const blitz::Array<double,2>& X_phi = d_interior_X_phi(e);
        const blitz::Array<VectorValue<double>,2>& X_dphi = d_interior_X_dphi(e);

        const std::vector<unsigned int>* const J_bar_dof_indices = (d_use_fbar_projection ? &d_interior_J_bar_dof_indices(e) : NULL);
        const blitz::Array<double,2>* const J_bar_phi = (d_use_fbar_projection ? &d_interior_J_bar_phi(e) : NULL);

        // Compute the elemental right-hand-side entries.
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            F_e[i].resize(F_dof_indices[i].size());
        }
        const int n_qp    = F_phi_JxW.extent(blitz::firstDim);
        const int n_basis = F_phi_JxW.extent(blitz::secondDim);
        for (int qp = 0; qp < n_qp; ++qp)
        {
            const Point& s_qp = q_point(qp);
            const Point& X_qp = compute_coordinate(qp,X,X_phi,X_dof_indices);
            const TensorValue<double> dX_ds = compute_coordinate_mapping_jacobian(qp,X,X_dphi,X_dof_indices,J_bar,J_bar_phi,J_bar_dof_indices);

            // Compute the value of the first Piola-Kirchhoff stress tensor at
            // the quadrature point and add the corresponding forces to the
            // right-hand-side vector.
            const TensorValue<double> PP = d_PK1_stress_fcn(dX_ds,X_qp,s_qp,elem,time,d_PK1_stress_fcn_ctx);
            for (int k = 0; k < n_basis; ++k)
            {
                const VectorValue<double> F_qp = -PP*F_dphi_JxW(qp,k);
                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    F_e[i](k) += F_qp(i);
                }
            }

            if (d_lag_body_force_fcn != NULL)
            {
                // Compute the value of the body force at the quadrature point
                // and add the corresponding forces to the right-hand-side
                // vector.
                const VectorValue<double> F_b = d_lag_body_force_fcn(dX_ds,X_qp,s_qp,elem,time,d_lag_body_force_fcn_ctx);
                for (int k = 0; k < n_basis; ++k)
                {
                    const VectorValue<double> F_qp = F_phi_JxW(qp,k)*F_b;
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_e[i](k) += F_qp(i);
                    }
                }
            }
        }

        // Loop over the element boundaries.
        for (unsigned short int side = 0; side < elem->n_sides(); ++side)
        {
            // Skip non-physical boundaries.
            const bool at_physical_bdry = d_interior_elem_side_at_physical_bdry(e)(side);
            if (!at_physical_bdry) continue;

            // Determine whether we need to compute surface forces along this
            // part of the physical boundary; if not, skip the present side.
            const bool at_dirichlet_bdry = d_interior_elem_side_at_dirichlet_bdry(e)(side);
            const bool compute_transmission_force = (( d_split_interior_and_bdry_forces && !at_dirichlet_bdry) ||
                                                     (!d_split_interior_and_bdry_forces &&  at_dirichlet_bdry));
            const bool compute_pressure = (d_split_interior_and_bdry_forces && d_lag_pressure_fcn != NULL);
            if (!(compute_transmission_force || compute_pressure)) continue;

            // Lookup cached data.
            const blitz::Array<Point,1>& q_point_face = d_interior_q_point_face(e)(side);

            const blitz::Array<VectorValue<double>,1>& F_normal_face = d_interior_F_normal_face(e)(side);
            const blitz::Array<double,2>& F_phi_JxW_face = d_interior_F_phi_JxW_face(e)(side);

            const blitz::Array<double,2>& X_phi_face = d_interior_X_phi_face(e)(side);
            const blitz::Array<VectorValue<double>,2>& X_dphi_face = d_interior_X_dphi_face(e)(side);

            const blitz::Array<double,2>* const J_bar_phi_face = (d_use_fbar_projection ? &d_interior_J_bar_phi_face(e)(side) : NULL);

            // Compute the elemental right-hand-side entries.
            const int n_qp    = F_phi_JxW_face.extent(blitz::firstDim);
            const int n_basis = F_phi_JxW_face.extent(blitz::secondDim);
            for (int qp = 0; qp < n_qp; ++qp)
            {
                const Point& s_qp = q_point_face(qp);
                const Point& X_qp = compute_coordinate(qp,X,X_phi_face,X_dof_indices);
                const TensorValue<double> dX_ds = compute_coordinate_mapping_jacobian(qp,X,X_dphi_face,X_dof_indices,J_bar,J_bar_phi_face,J_bar_dof_indices);

                VectorValue<double> F;

                if (compute_transmission_force)
                {
                    // Compute the value of the first Piola-Kirchhoff stress
                    // tensor at the quadrature point and add the corresponding
                    // force to the right-hand-side vector.
                    const TensorValue<double> PP = d_PK1_stress_fcn(dX_ds,X_qp,s_qp,elem,time,d_PK1_stress_fcn_ctx);
                    F += PP*F_normal_face(qp);
                }

                if (compute_pressure && !at_dirichlet_bdry)
                {
                    // Compute the value of the pressure at the quadrature point
                    // and add the corresponding force to the right-hand-side
                    // vector.
                    const double P = d_lag_pressure_fcn(dX_ds,X_qp,s_qp,elem,side,time,d_lag_pressure_fcn_ctx);
                    const double J = dX_ds.det();
                    F -= P*J*tensor_inverse_transpose(dX_ds,NDIM)*F_normal_face(qp);
                }

                for (int k = 0; k < n_basis; ++k)
                {
                    const VectorValue<double> F_qp = F_phi_JxW_face(qp,k)*F;
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_e[i](k) += F_qp(i);
                    }
                }
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions) and
        // add the elemental contributions to the global vector.
        for (int i = 0; i < NDIM; ++i)
        {
            F_dof_map.constrain_element_vector(F_e[i], F_dof_indices[i]);
            F->add_vector(F_e[i], F_dof_indices[i]);
        }
    }

    // Solve for G.
    d_fe_data_manager->computeL2Projection(G, *F, FORCE_SYSTEM_NAME, d_use_consistent_mass_matrix);
    return;
}// computeInteriorForceDensity

void
IBFEHierarchyIntegrator::spreadTransmissionForceDensity(
    const int f_data_idx,
    NumericVector<double>& X_ghost,
    NumericVector<double>* J_bar_ghost,
    const double& time)
{
    if (!d_split_interior_and_bdry_forces) return;

    // Loop over the patches to spread the transmission elastic force density
    // onto the grid.
    const int level_num = d_fe_data_manager->getLevelNumber();
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(f_data_idx);
        const Box<NDIM>& box = f_data->getGhostBox();

        // Loop over the elements and compute the values to be spread and the
        // positions of the quadrature points.
        const std::vector<Elem*>& active_patch_elems = d_fe_data_manager->getActivePatchElements()[local_patch_num];
        if (active_patch_elems.empty()) continue;
        const std::vector<Elem*>::const_iterator el_begin = active_patch_elems.begin();
        const std::vector<Elem*>::const_iterator el_end   = active_patch_elems.end();

        const blitz::Array<std::vector<std::vector<unsigned int> >,1>& patch_X_dof_indices = d_transmission_X_dof_indices(local_patch_num);
        const blitz::Array<std::vector<unsigned int>,1>* const patch_J_bar_dof_indices = (d_use_fbar_projection ? &d_transmission_J_bar_dof_indices(local_patch_num) : NULL);

        const blitz::Array<blitz::Array<bool,1>,1>& patch_elem_side_at_physical_bdry = d_transmission_elem_side_at_physical_bdry(local_patch_num);
        const blitz::Array<blitz::Array<bool,1>,1>& patch_elem_side_at_dirichlet_bdry = d_transmission_elem_side_at_dirichlet_bdry(local_patch_num);

        const blitz::Array<blitz::Array<blitz::Array<libMesh::Point,1>,1>,1>& patch_q_point_face = d_transmission_q_point_face(local_patch_num);

        const blitz::Array<blitz::Array<blitz::Array<libMesh::VectorValue<double>,1>,1>,1>& patch_X_normal_JxW_face = d_transmission_X_normal_JxW_face(local_patch_num);
        const blitz::Array<blitz::Array<blitz::Array<double,2>,1>,1>& patch_X_phi_face = d_transmission_X_phi_face(local_patch_num);
        const blitz::Array<blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1>,1>& patch_X_dphi_face = d_transmission_X_dphi_face(local_patch_num);

        const blitz::Array<blitz::Array<blitz::Array<double,2>,1>,1>* const patch_J_bar_phi_face = (d_use_fbar_projection ? &d_transmission_J_bar_phi_face(local_patch_num) : NULL);

        int qp_offset = 0;
        std::vector<double> T_bdry, X_bdry;
        for (std::vector<Elem*>::const_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            const unsigned int e = std::distance(el_begin,el_it);

            // Loop over the element boundaries.
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries.
                const bool at_physical_bdry = patch_elem_side_at_physical_bdry(e)(side);
                if (!at_physical_bdry) continue;

                // Determine whether we need to compute surface forces along
                // this part of the physical boundary; if not, skip the present
                // side.
                const bool at_dirichlet_bdry = patch_elem_side_at_dirichlet_bdry(e)(side);
                const bool compute_transmission_force = d_split_interior_and_bdry_forces && !at_dirichlet_bdry;
                const bool compute_pressure = (d_split_interior_and_bdry_forces && d_lag_pressure_fcn != NULL);
                if (!(compute_transmission_force || compute_pressure)) continue;

                // Lookup cached data.
                const blitz::Array<Point,1>& q_point_face = patch_q_point_face(e)(side);

                const std::vector<std::vector<unsigned int> >& X_dof_indices = patch_X_dof_indices(e);
                const blitz::Array<VectorValue<double>,1>& X_normal_JxW_face = patch_X_normal_JxW_face(e)(side);
                const blitz::Array<double,2>& X_phi_face = patch_X_phi_face(e)(side);
                const blitz::Array<VectorValue<double>,2>& X_dphi_face = patch_X_dphi_face(e)(side);

                const std::vector<unsigned int>* const J_bar_dof_indices = (d_use_fbar_projection ? &(*patch_J_bar_dof_indices)(e) : NULL);
                const blitz::Array<double,2>* const J_bar_phi_face = (d_use_fbar_projection ? &(*patch_J_bar_phi_face)(e)(side) : NULL);

                // Loop over boundary quadrature points.
                const int n_qp = q_point_face.size();
                T_bdry.resize(T_bdry.size()+NDIM*n_qp,0.0);
                X_bdry.resize(X_bdry.size()+NDIM*n_qp,0.0);
                for (int qp = 0; qp < n_qp; ++qp, ++qp_offset)
                {
                    const Point& s_qp = q_point_face(qp);
                    const Point& X_qp = compute_coordinate(qp,X_ghost,X_phi_face,X_dof_indices);
                    const TensorValue<double> dX_ds = compute_coordinate_mapping_jacobian(qp,X_ghost,X_dphi_face,X_dof_indices,J_bar_ghost,J_bar_phi_face,J_bar_dof_indices);

                    VectorValue<double> F;
                    if (compute_transmission_force)
                    {
                        // Compute the value of the first Piola-Kirchhoff stress
                        // tensor at the quadrature point and compute the
                        // corresponding force.
                        const TensorValue<double> PP = d_PK1_stress_fcn(dX_ds,X_qp,s_qp,elem,time,d_PK1_stress_fcn_ctx);
                        F -= PP*X_normal_JxW_face(qp);
                    }

                    if (compute_pressure)
                    {
                        // Compute the value of the pressure at the quadrature
                        // point and compute the corresponding force.
                        const double P = d_lag_pressure_fcn(dX_ds,X_qp,s_qp,elem,side,time,d_lag_pressure_fcn_ctx);
                        const double J = dX_ds.det();
                        F -= P*J*tensor_inverse_transpose(dX_ds,NDIM)*X_normal_JxW_face(qp);
                    }

                    const int idx = NDIM*qp_offset;
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        T_bdry[idx+i] = F(i);
                    }

                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        X_bdry[idx+i] = X_qp(i);
                    }
                }
            }
        }

        // Spread the boundary forces to the grid.
        if (qp_offset > 0)
        {
            LEInteractor::spread(f_data, T_bdry, NDIM, X_bdry, NDIM, patch, box, d_fe_data_manager->getSpreadWeightingFunction());
        }
    }
    return;
}// spreadTransmissionForceDensity

void
IBFEHierarchyIntegrator::initializeCoordinates()
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    const bool identity_mapping = d_coordinate_mapping_fcn == NULL;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num) > 0)
        {
            libmesh_assert(n->n_vars(X_sys_num) == NDIM);
            const Point& s = *n;
            const Point X = identity_mapping ? s : d_coordinate_mapping_fcn(s, d_coordinate_mapping_fcn_ctx);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int dof_index = n->dof_number(X_sys_num,d,0);
                X_coords.set(dof_index,X(d));
            }
        }
    }
    X_coords.close();
    return;
}// initializeCoordinates

void
IBFEHierarchyIntegrator::updateCoordinateMapping()
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    System& X_mapping_system = equation_systems->get_system<System>(COORD_MAPPING_SYSTEM_NAME);
    const unsigned int X_mapping_sys_num = X_mapping_system.number();
    NumericVector<double>& dX_coords = *X_mapping_system.solution;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num) > 0)
        {
            libmesh_assert(n->n_vars(X_sys_num) == NDIM);
            libmesh_assert(n->n_vars(X_mapping_sys_num) == NDIM);
            const Point& s = *n;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int X_dof_index = n->dof_number(X_sys_num,d,0);
                const int dX_dof_index = n->dof_number(X_mapping_sys_num,d,0);
                dX_coords.set(dX_dof_index,X_coords(X_dof_index)-s(d));
            }
        }
    }
    dX_coords.close();
    return;
}// updateCoordinateMapping

void
IBFEHierarchyIntegrator::computeCachedProjectedDilatationalStrainFEData()
{
    if (!d_use_fbar_projection) return;

    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();

    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    Order quad_order = CONSTANT;
    switch (d_fe_order)
    {
        case FIRST:
            quad_order = THIRD;
            break;
        case SECOND:
            quad_order = FIFTH;
            break;
        default:
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ERROR("this statement should not be reached\n");
#endif
    }
    QGauss qrule(dim, quad_order);
    QGauss qrule_face(dim-1, quad_order);

    System& J_bar_system = equation_systems->get_system<System>(PROJ_STRAIN_SYSTEM_NAME);
    const DofMap& J_bar_dof_map = J_bar_system.get_dof_map();
    std::vector<unsigned int> J_bar_dof_indices;

    AutoPtr<FEBase> J_bar_fe(FEBase::build(dim, J_bar_dof_map.variable_type(0)));
    J_bar_fe->attach_quadrature_rule(&qrule);
    const std::vector<double>& J_bar_JxW = J_bar_fe->get_JxW();
    const std::vector<std::vector<double> >& J_bar_phi = J_bar_fe->get_phi();

    System& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(0) == X_dof_map.variable_type(d));
    }
#endif

    AutoPtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
    X_fe->attach_quadrature_rule(&qrule);
    const std::vector<std::vector<VectorValue<double> > >& X_dphi = X_fe->get_dphi();

    // Compute cached values.
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();

    const unsigned int sz = std::distance(el_begin,el_end);

    d_proj_strain_J_bar_dof_indices.resize(sz);
    d_proj_strain_J_bar_phi_JxW.resize(sz);

    d_proj_strain_X_dof_indices.resize(sz);
    d_proj_strain_X_dphi.resize(sz);

    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        const unsigned int e = std::distance(el_begin,el_it);

        // Compute cached data.
        J_bar_fe->reinit(elem);
        J_bar_dof_map.dof_indices(elem, J_bar_dof_indices);
        d_proj_strain_J_bar_dof_indices(e) = J_bar_dof_indices;

        X_fe->reinit(elem);
        d_proj_strain_X_dof_indices(e).resize(NDIM);
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[i], i);
            d_proj_strain_X_dof_indices(e)[i] = X_dof_indices[i];
        }

        const int n_qp = qrule.n_points();
        d_proj_strain_J_bar_phi_JxW(e).resize(n_qp,J_bar_phi.size());
        d_proj_strain_X_dphi(e).resize(n_qp,X_dphi.size());
        for (int qp = 0; qp < n_qp; ++qp)
        {
            for (int k = 0; k < int(J_bar_phi.size()); ++k)
            {
                d_proj_strain_J_bar_phi_JxW(e)(qp,k) = J_bar_JxW[qp]*J_bar_phi[k][qp];
            }
            for (int k = 0; k < int(X_dphi.size()); ++k)
            {
                d_proj_strain_X_dphi(e)(qp,k) = X_dphi[k][qp];
            }
        }
    }
    return;
}// computeCachedProjectedDilatationalStrainFEData

void
IBFEHierarchyIntegrator::computeCachedInteriorForceDensityFEData()
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();

    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    Order quad_order = CONSTANT;
    switch (d_fe_order)
    {
        case FIRST:
            quad_order = THIRD;
            break;
        case SECOND:
            quad_order = FIFTH;
            break;
        default:
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ERROR("this statement should not be reached\n");
#endif
    }
    QGauss qrule(dim, quad_order);
    QGauss qrule_face(dim-1, quad_order);

    System& F_system = equation_systems->get_system<System>(FORCE_SYSTEM_NAME);
    const DofMap& F_dof_map = F_system.get_dof_map();
    std::vector<std::vector<unsigned int> > F_dof_indices(NDIM);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(F_dof_map.variable_type(0) == F_dof_map.variable_type(d));
    }
#endif

    AutoPtr<FEBase> F_fe(FEBase::build(dim, F_dof_map.variable_type(0)));
    F_fe->attach_quadrature_rule(&qrule);
    const std::vector<Point>& q_point = F_fe->get_xyz();
    const std::vector<double>& F_JxW = F_fe->get_JxW();
    const std::vector<std::vector<double> >& F_phi = F_fe->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& F_dphi = F_fe->get_dphi();

    AutoPtr<FEBase> F_fe_face(FEBase::build(dim, F_dof_map.variable_type(0)));
    F_fe_face->attach_quadrature_rule(&qrule_face);
    const std::vector<Point>& q_point_face = F_fe_face->get_xyz();
    const std::vector<double>& F_JxW_face = F_fe_face->get_JxW();
    const std::vector<std::vector<double> >& F_phi_face = F_fe_face->get_phi();
    const std::vector<Point>& F_normal_face = F_fe_face->get_normals();

    System& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(0) == X_dof_map.variable_type(d));
    }
#endif

    AutoPtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
    X_fe->attach_quadrature_rule(&qrule);
    const std::vector<std::vector<double> >& X_phi = X_fe->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& X_dphi = X_fe->get_dphi();

    AutoPtr<FEBase> X_fe_face(FEBase::build(dim, X_dof_map.variable_type(0)));
    X_fe_face->attach_quadrature_rule(&qrule_face);
    const std::vector<std::vector<double> >& X_phi_face = X_fe_face->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& X_dphi_face = X_fe_face->get_dphi();

    const DofMap* J_bar_dof_map = NULL;
    std::vector<unsigned int> J_bar_dof_indices;
    AutoPtr<FEBase> J_bar_fe;
    const std::vector<std::vector<double> >* J_bar_phi = NULL;
    AutoPtr<FEBase> J_bar_fe_face;
    const std::vector<std::vector<double> >* J_bar_phi_face = NULL;
    if (d_use_fbar_projection)
    {
        System& J_bar_system = equation_systems->get_system<System>(PROJ_STRAIN_SYSTEM_NAME);
        J_bar_dof_map = &J_bar_system.get_dof_map();

        J_bar_fe = FEBase::build(dim, J_bar_dof_map->variable_type(0));
        J_bar_fe->attach_quadrature_rule(&qrule);
        J_bar_phi = &(J_bar_fe->get_phi());

        J_bar_fe_face = FEBase::build(dim, J_bar_dof_map->variable_type(0));
        J_bar_fe_face->attach_quadrature_rule(&qrule_face);
        J_bar_phi_face = &(J_bar_fe_face->get_phi());
    }

    // Compute cached values.
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();

    const unsigned int sz = std::distance(el_begin,el_end);

    d_interior_q_point.resize(sz);

    d_interior_F_dof_indices.resize(sz);
    d_interior_F_phi_JxW.resize(sz);
    d_interior_F_dphi_JxW.resize(sz);

    d_interior_X_dof_indices.resize(sz);
    d_interior_X_phi.resize(sz);
    d_interior_X_dphi.resize(sz);

    d_interior_q_point_face.resize(sz);

    d_interior_F_normal_face.resize(sz);
    d_interior_F_phi_JxW_face.resize(sz);

    d_interior_X_phi_face.resize(sz);
    d_interior_X_dphi_face.resize(sz);

    d_interior_elem_side_at_physical_bdry.resize(sz);
    d_interior_elem_side_at_dirichlet_bdry.resize(sz);

    if (d_use_fbar_projection)
    {
        d_interior_J_bar_dof_indices.resize(sz);
        d_interior_J_bar_phi.resize(sz);

        d_interior_J_bar_phi_face.resize(sz);
    }

    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        const unsigned int e = std::distance(el_begin,el_it);

        // Compute cached data.
        F_fe->reinit(elem);
        d_interior_F_dof_indices(e).resize(NDIM);
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            F_dof_map.dof_indices(elem, F_dof_indices[i], i);
            d_interior_F_dof_indices(e)[i] = F_dof_indices[i];
        }

        X_fe->reinit(elem);
        d_interior_X_dof_indices(e).resize(NDIM);
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[i], i);
            d_interior_X_dof_indices(e)[i] = X_dof_indices[i];
        }

        if (d_use_fbar_projection)
        {
            J_bar_fe->reinit(elem);
            J_bar_dof_map->dof_indices(elem, J_bar_dof_indices);
            d_interior_J_bar_dof_indices(e) = J_bar_dof_indices;
        }

        const int n_qp = qrule.n_points();
        d_interior_q_point(e).resize(n_qp);
        d_interior_F_phi_JxW(e).resize(n_qp,F_phi.size());
        d_interior_F_dphi_JxW(e).resize(n_qp,F_dphi.size());
        d_interior_X_phi(e).resize(n_qp,X_phi.size());
        d_interior_X_dphi(e).resize(n_qp,X_dphi.size());
        if (d_use_fbar_projection) d_interior_J_bar_phi(e).resize(n_qp,J_bar_phi->size());
        for (int qp = 0; qp < n_qp; ++qp)
        {
            d_interior_q_point(e)(qp) = q_point[qp];
            for (int k = 0; k < int(F_phi.size()); ++k)
            {
                d_interior_F_phi_JxW(e)(qp,k) = F_JxW[qp]*F_phi[k][qp];
            }
            for (int k = 0; k < int(F_dphi.size()); ++k)
            {
                d_interior_F_dphi_JxW(e)(qp,k) = F_dphi[k][qp]*F_JxW[qp];
            }
            for (int k = 0; k < int(X_phi.size()); ++k)
            {
                d_interior_X_phi(e)(qp,k) = X_phi[k][qp];
            }
            for (int k = 0; k < int(X_dphi.size()); ++k)
            {
                d_interior_X_dphi(e)(qp,k) = X_dphi[k][qp];
            }
            if (d_use_fbar_projection)
            {
                for (int k = 0; k < int(J_bar_phi->size()); ++k)
                {
                    d_interior_J_bar_phi(e)(qp,k) = (*J_bar_phi)[k][qp];
                }
            }
        }

        // Loop over the element boundaries.
        const unsigned int n_sides = elem->n_sides();
        d_interior_elem_side_at_physical_bdry(e).resize(n_sides);
        d_interior_elem_side_at_dirichlet_bdry(e).resize(n_sides);
        d_interior_q_point_face(e).resize(n_sides);
        d_interior_F_normal_face(e).resize(n_sides);
        d_interior_F_phi_JxW_face(e).resize(n_sides);
        d_interior_X_phi_face(e).resize(n_sides);
        d_interior_X_dphi_face(e).resize(n_sides);
        if (d_use_fbar_projection) d_interior_J_bar_phi_face(e).resize(n_sides);
        for (unsigned short int side = 0; side < n_sides; ++side)
        {
            // Skip non-physical boundaries.
            const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
            bool at_physical_bdry  = false;
            bool at_dirichlet_bdry = false;
            for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
            {
                const short int bdry_id = *cit;
                at_physical_bdry  = at_physical_bdry  || (elem->neighbor(side) == NULL && !X_dof_map.is_periodic_boundary(bdry_id));
                at_dirichlet_bdry = at_dirichlet_bdry || (bdry_id == FEDataManager::DIRICHLET_BDRY_ID);
            }
            d_interior_elem_side_at_physical_bdry (e)(side) = at_physical_bdry ;
            d_interior_elem_side_at_dirichlet_bdry(e)(side) = at_dirichlet_bdry;
            if (!at_physical_bdry) continue;

            // Determine whether we need to compute surface forces along the
            // physical boundary.
            const bool compute_transmission_force = (( d_split_interior_and_bdry_forces && !at_dirichlet_bdry) ||
                                                     (!d_split_interior_and_bdry_forces &&  at_dirichlet_bdry));
            const bool compute_pressure = (d_split_interior_and_bdry_forces && d_lag_pressure_fcn != NULL);
            if (!(compute_transmission_force || compute_pressure)) continue;

            F_fe_face->reinit(elem, side);
            X_fe_face->reinit(elem, side);
            if (d_use_fbar_projection) J_bar_fe_face->reinit(elem, side);

            // Compute cached data.
            const int n_qp_face = qrule_face.n_points();
            d_interior_q_point_face(e)(side).resize(n_qp_face);
            d_interior_F_normal_face(e)(side).resize(n_qp_face);
            d_interior_F_phi_JxW_face(e)(side).resize(n_qp_face,F_phi_face.size());
            d_interior_X_phi_face(e)(side).resize(n_qp_face,X_phi_face.size());
            d_interior_X_dphi_face(e)(side).resize(n_qp_face,X_dphi_face.size());
            if (d_use_fbar_projection) d_interior_J_bar_phi_face(e)(side).resize(n_qp_face,J_bar_phi_face->size());
            for (int qp = 0; qp < n_qp_face; ++qp)
            {
                d_interior_q_point_face(e)(side)(qp) = q_point_face[qp];
                d_interior_F_normal_face(e)(side)(qp) = F_normal_face[qp];
                for (int k = 0; k < int(F_phi_face.size()); ++k)
                {
                    d_interior_F_phi_JxW_face(e)(side)(qp,k) = F_phi_face[k][qp]*F_JxW_face[qp];
                }
                for (int k = 0; k < int(X_phi_face.size()); ++k)
                {
                    d_interior_X_phi_face(e)(side)(qp,k) = X_phi_face[k][qp];
                }
                for (int k = 0; k < int(X_dphi_face.size()); ++k)
                {
                    d_interior_X_dphi_face(e)(side)(qp,k) = X_dphi_face[k][qp];
                }
                if (d_use_fbar_projection)
                {
                    for (int k = 0; k < int(J_bar_phi_face->size()); ++k)
                    {
                        d_interior_J_bar_phi_face(e)(side)(qp,k) = (*J_bar_phi_face)[k][qp];
                    }
                }
            }
        }
    }
    return;
}// computeCachedInteriorForceDensityFEData

void
IBFEHierarchyIntegrator::computeCachedTransmissionForceDensityFEData()
{
    if (!d_split_interior_and_bdry_forces) return;

    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();

    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    QBase* qrule = d_fe_data_manager->getQuadratureRule();
    AutoPtr<QBase> qrule_face = QBase::build(qrule->type(), dim-1, qrule->get_order());

    System& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(0) == X_dof_map.variable_type(d));
    }
#endif

    AutoPtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
    X_fe->attach_quadrature_rule(qrule);

    AutoPtr<FEBase> X_fe_face(FEBase::build(dim, X_dof_map.variable_type(0)));
    X_fe_face->attach_quadrature_rule(qrule_face.get());
    const std::vector<Point>& q_point_face = X_fe_face->get_xyz();
    const std::vector<double>& X_JxW_face = X_fe_face->get_JxW();
    const std::vector<std::vector<double> >& X_phi_face = X_fe_face->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& X_dphi_face = X_fe_face->get_dphi();
    const std::vector<Point>& X_normal_face = X_fe_face->get_normals();

    const DofMap* J_bar_dof_map = NULL;
    std::vector<unsigned int> J_bar_dof_indices;
    AutoPtr<FEBase> J_bar_fe;
    AutoPtr<FEBase> J_bar_fe_face;
    const std::vector<std::vector<double> >* J_bar_phi_face = NULL;
    if (d_use_fbar_projection)
    {
        System& J_bar_system = equation_systems->get_system<System>(PROJ_STRAIN_SYSTEM_NAME);
        J_bar_dof_map = &J_bar_system.get_dof_map();

        J_bar_fe = FEBase::build(dim, J_bar_dof_map->variable_type(0));
        J_bar_fe->attach_quadrature_rule(qrule);

        J_bar_fe_face = FEBase::build(dim, J_bar_dof_map->variable_type(0));
        J_bar_fe_face->attach_quadrature_rule(qrule_face.get());
        J_bar_phi_face = &(J_bar_fe_face->get_phi());
    }

    // Compute cached values.
    const int level_num = d_fe_data_manager->getLevelNumber();
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);

    const int n_local_patches = level->getProcessorMapping().getNumberOfLocalIndices();
    d_transmission_X_dof_indices.resize(n_local_patches);
    if (d_use_fbar_projection) d_transmission_J_bar_dof_indices.resize(n_local_patches);
    d_transmission_elem_side_at_physical_bdry.resize(n_local_patches);
    d_transmission_elem_side_at_dirichlet_bdry.resize(n_local_patches);
    d_transmission_q_point_face.resize(n_local_patches);
    d_transmission_X_normal_JxW_face.resize(n_local_patches);
    d_transmission_X_phi_face.resize(n_local_patches);
    d_transmission_X_dphi_face.resize(n_local_patches);
    if (d_use_fbar_projection) d_transmission_J_bar_phi_face.resize(n_local_patches);

    // Loop over the patches.
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // Loop over the elements.
        const std::vector<Elem*>& active_patch_elems = d_fe_data_manager->getActivePatchElements()[local_patch_num];
        if (active_patch_elems.empty()) continue;
        const std::vector<Elem*>::const_iterator el_begin = active_patch_elems.begin();
        const std::vector<Elem*>::const_iterator el_end   = active_patch_elems.end();

        blitz::Array<std::vector<std::vector<unsigned int> >,1>& patch_X_dof_indices = d_transmission_X_dof_indices(local_patch_num);
        blitz::Array<std::vector<unsigned int>,1>* const patch_J_bar_dof_indices = (d_use_fbar_projection ? &d_transmission_J_bar_dof_indices(local_patch_num) : NULL);

        blitz::Array<blitz::Array<bool,1>,1>& patch_elem_side_at_physical_bdry = d_transmission_elem_side_at_physical_bdry(local_patch_num);
        blitz::Array<blitz::Array<bool,1>,1>& patch_elem_side_at_dirichlet_bdry = d_transmission_elem_side_at_dirichlet_bdry(local_patch_num);

        blitz::Array<blitz::Array<blitz::Array<libMesh::Point,1>,1>,1>& patch_q_point_face = d_transmission_q_point_face(local_patch_num);

        blitz::Array<blitz::Array<blitz::Array<libMesh::VectorValue<double>,1>,1>,1>& patch_X_normal_JxW_face = d_transmission_X_normal_JxW_face(local_patch_num);
        blitz::Array<blitz::Array<blitz::Array<double,2>,1>,1>& patch_X_phi_face = d_transmission_X_phi_face(local_patch_num);
        blitz::Array<blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1>,1>& patch_X_dphi_face = d_transmission_X_dphi_face(local_patch_num);

        blitz::Array<blitz::Array<blitz::Array<double,2>,1>,1>* const patch_J_bar_phi_face = (d_use_fbar_projection ? &d_transmission_J_bar_phi_face(local_patch_num) : NULL);

        const unsigned int sz = active_patch_elems.size();
        patch_X_dof_indices.resize(sz);
        if (d_use_fbar_projection) patch_J_bar_dof_indices->resize(sz);
        patch_elem_side_at_physical_bdry.resize(sz);
        patch_elem_side_at_dirichlet_bdry.resize(sz);
        patch_q_point_face.resize(sz);
        patch_X_normal_JxW_face.resize(sz);
        patch_X_phi_face.resize(sz);
        patch_X_dphi_face.resize(sz);
        if (d_use_fbar_projection) patch_J_bar_phi_face->resize(sz);

        for (std::vector<Elem*>::const_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            const unsigned int e = std::distance(el_begin,el_it);

            // Compute cached data.
            X_fe->reinit(elem);
            patch_X_dof_indices(e).resize(NDIM);
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                X_dof_map.dof_indices(elem, X_dof_indices[i], i);
                patch_X_dof_indices(e)[i] = X_dof_indices[i];
            }

            if (d_use_fbar_projection)
            {
                J_bar_fe->reinit(elem);
                J_bar_dof_map->dof_indices(elem, J_bar_dof_indices);
                (*patch_J_bar_dof_indices)(e) = J_bar_dof_indices;
            }

            // Loop over the element boundaries.
            const unsigned int n_sides = elem->n_sides();
            patch_elem_side_at_physical_bdry(e).resize(n_sides);
            patch_elem_side_at_dirichlet_bdry(e).resize(n_sides);
            patch_q_point_face(e).resize(n_sides);
            patch_X_normal_JxW_face(e).resize(n_sides);
            patch_X_phi_face(e).resize(n_sides);
            patch_X_dphi_face(e).resize(n_sides);
            if (d_use_fbar_projection) (*patch_J_bar_phi_face)(e).resize(n_sides);
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries.
                const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
                bool at_physical_bdry  = false;
                bool at_dirichlet_bdry = false;
                for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
                {
                    const short int bdry_id = *cit;
                    at_physical_bdry  = at_physical_bdry  || (elem->neighbor(side) == NULL && !X_dof_map.is_periodic_boundary(bdry_id));
                    at_dirichlet_bdry = at_dirichlet_bdry || (bdry_id == FEDataManager::DIRICHLET_BDRY_ID);
                }
                patch_elem_side_at_physical_bdry (e)(side) = at_physical_bdry ;
                patch_elem_side_at_dirichlet_bdry(e)(side) = at_dirichlet_bdry;
                if (!at_physical_bdry) continue;

                // Determine whether we need to compute surface forces along
                // this part of the physical boundary; if not, skip the present
                // side.
                const bool compute_transmission_force = d_split_interior_and_bdry_forces && !at_dirichlet_bdry;
                const bool compute_pressure = (d_split_interior_and_bdry_forces && d_lag_pressure_fcn != NULL);
                if (!(compute_transmission_force || compute_pressure)) continue;

                X_fe_face->reinit(elem, side);
                if (d_use_fbar_projection) J_bar_fe_face->reinit(elem, side);

                // Compute cached data.
                const int n_qp_face = qrule_face->n_points();
                patch_q_point_face(e)(side).resize(n_qp_face);
                patch_X_normal_JxW_face(e)(side).resize(n_qp_face);
                patch_X_phi_face(e)(side).resize(n_qp_face,X_phi_face.size());
                patch_X_dphi_face(e)(side).resize(n_qp_face,X_dphi_face.size());
                if (d_use_fbar_projection) (*patch_J_bar_phi_face)(e)(side).resize(n_qp_face,J_bar_phi_face->size());
                for (int qp = 0; qp < n_qp_face; ++qp)
                {
                    patch_q_point_face(e)(side)(qp) = q_point_face[qp];
                    patch_X_normal_JxW_face(e)(side)(qp) = X_normal_face[qp]*X_JxW_face[qp];
                    for (int k = 0; k < int(X_phi_face.size()); ++k)
                    {
                        patch_X_phi_face(e)(side)(qp,k) = X_phi_face[k][qp];
                    }
                    for (int k = 0; k < int(X_dphi_face.size()); ++k)
                    {
                        patch_X_dphi_face(e)(side)(qp,k) = X_dphi_face[k][qp];
                    }
                    if (d_use_fbar_projection)
                    {
                        for (int k = 0; k < int(J_bar_phi_face->size()); ++k)
                        {
                            (*patch_J_bar_phi_face)(e)(side)(qp,k) = (*J_bar_phi_face)[k][qp];
                        }
                    }
                }
            }
        }
    }
    return;
}// computeCachedTransmissionForceDensityFEData

void
IBFEHierarchyIntegrator::getFromInput(
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

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault("dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault("dt_max_time_min", d_dt_max_time_min);

    d_fe_order = Utility::string_to_enum<Order>(db->getStringWithDefault("fe_order", Utility::enum_to_string<Order>(d_fe_order)));
    d_fe_family = Utility::string_to_enum<FEFamily>(db->getStringWithDefault("fe_family", Utility::enum_to_string<FEFamily>(d_fe_family)));
    d_split_interior_and_bdry_forces = db->getBoolWithDefault("split_interior_and_bdry_forces", d_split_interior_and_bdry_forces);
    d_use_consistent_mass_matrix = db->getBoolWithDefault("use_consistent_mass_matrix", d_use_consistent_mass_matrix);

    d_use_fbar_projection = db->getBoolWithDefault("use_fbar_projection", d_use_fbar_projection);
    if (d_use_fbar_projection)
    {
        d_J_bar_fe_order = Utility::string_to_enum<Order>(db->getStringWithDefault("J_bar_fe_order", Utility::enum_to_string<Order>(d_J_bar_fe_order)));
        d_J_bar_fe_family = Utility::string_to_enum<FEFamily>(db->getStringWithDefault("J_bar_fe_family", Utility::enum_to_string<FEFamily>(d_J_bar_fe_family)));
    }

    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    if (!is_from_restart)
    {
        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);
        d_mark_input_file_name = db->getStringWithDefault("marker_input_file_name", d_mark_input_file_name);
    }
    return;
}// getFromInput

void
IBFEHierarchyIntegrator::getFromRestart()
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

    int ver = db->getInteger("IB_FE_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_FE_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    d_split_interior_and_bdry_forces = db->getBool("d_split_interior_and_bdry_forces");
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
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

    d_mark_input_file_name = db->getString("d_mark_input_file_name");
    d_num_mark = db->getInteger("d_num_mark");
    d_mark_init_posns.resize(NDIM*d_num_mark);
    if (d_num_mark > 0)
    {
        db->getDoubleArray("d_mark_init_posns", &d_mark_init_posns[0], d_mark_init_posns.size());
    }
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBFEHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
