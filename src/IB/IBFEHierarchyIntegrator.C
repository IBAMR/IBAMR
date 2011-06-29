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
#include <ibtk/LMarkerSetData.h>
#include <ibtk/LMarkerUtilities.h>
#include <ibtk/libmesh_utilities.h>

// LIBMESH INCLUDES
#include <boundary_info.h>
#include <dense_vector.h>
#include <dof_map.h>
#include <explicit_system.h>
#include <fe.h>
#include <fe_interface.h>
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
static Timer* t_initialize_hierarchy_integrator;
static Timer* t_initialize_hierarchy;
static Timer* t_advance_hierarchy;
static Timer* t_advance_hierarchy_init;
static Timer* t_advance_hierarchy_phase1;
static Timer* t_advance_hierarchy_phase2;
static Timer* t_advance_hierarchy_phase3;
static Timer* t_advance_hierarchy_phase4;
static Timer* t_advance_hierarchy_phase5;
static Timer* t_advance_hierarchy_phase6;
static Timer* t_advance_hierarchy_phase7;
static Timer* t_advance_hierarchy_phase8;
static Timer* t_advance_hierarchy_finalize;
static Timer* t_regrid_hierarchy;
static Timer* t_synchronize_hierarchy;
static Timer* t_synchronize_new_levels;
static Timer* t_reset_time_dependent_data;
static Timer* t_reset_data_to_preadvance_state;
static Timer* t_initialize_level_data;
static Timer* t_reset_hierarchy_configuration;
static Timer* t_apply_gradient_detector;
static Timer* t_put_to_database;

// Version of IBFEHierarchyIntegrator restart file data.
static const int IB_FE_HIERARCHY_INTEGRATOR_VERSION = 1;
}

const std::string IBFEHierarchyIntegrator::       COORDS_SYSTEM_NAME = "IB coordinates system";
const std::string IBFEHierarchyIntegrator::COORD_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string IBFEHierarchyIntegrator::        FORCE_SYSTEM_NAME = "IB force system";
const std::string IBFEHierarchyIntegrator::     VELOCITY_SYSTEM_NAME = "IB velocity system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEHierarchyIntegrator::IBFEHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
    FEDataManager* fe_data_manager,
    LDataManager* l_data_manager,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_num_parts(1),
      d_fe_data_managers(1,fe_data_manager),
      d_use_IB_spreading_operator(true),
      d_use_IB_interpolation_operator(true),
      d_split_interior_and_bdry_forces(false),
      d_use_jump_conditions(false),
      d_use_consistent_mass_matrix(true),
      d_fe_family(LAGRANGE),
      d_fe_order(INVALID_ORDER),
      d_quad_type(QGAUSS),
      d_quad_order(FIFTH),
      d_coordinate_mapping_fcns(d_num_parts,NULL),
      d_coordinate_mapping_fcn_ctxs(d_num_parts,NULL),
      d_PK1_stress_fcns(d_num_parts,NULL),
      d_PK1_stress_fcn_systems(d_num_parts),
      d_PK1_stress_fcn_ctxs(d_num_parts,NULL),
      d_lag_body_force_fcns(d_num_parts,NULL),
      d_lag_body_force_fcn_systems(d_num_parts),
      d_lag_body_force_fcn_ctxs(d_num_parts,NULL),
      d_lag_pressure_fcns(d_num_parts,NULL),
      d_lag_pressure_fcn_systems(d_num_parts),
      d_lag_pressure_fcn_ctxs(d_num_parts,NULL),
      d_lag_surface_force_fcns(d_num_parts,NULL),
      d_lag_surface_force_fcn_systems(d_num_parts),
      d_lag_surface_force_fcn_ctxs(d_num_parts,NULL),
      d_l_data_manager(l_data_manager),
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
    commonConstructor(input_db);
    return;
}// IBFEHierarchyIntegrator

IBFEHierarchyIntegrator::IBFEHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
    std::vector<FEDataManager*> fe_data_managers,
    LDataManager* l_data_manager,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_num_parts(fe_data_managers.size()),
      d_fe_data_managers(fe_data_managers),
      d_use_IB_spreading_operator(true),
      d_use_IB_interpolation_operator(true),
      d_split_interior_and_bdry_forces(false),
      d_use_jump_conditions(false),
      d_use_consistent_mass_matrix(true),
      d_fe_family(LAGRANGE),
      d_fe_order(INVALID_ORDER),
      d_quad_type(QGAUSS),
      d_quad_order(FIFTH),
      d_coordinate_mapping_fcns(d_num_parts,NULL),
      d_coordinate_mapping_fcn_ctxs(d_num_parts,NULL),
      d_PK1_stress_fcns(d_num_parts,NULL),
      d_PK1_stress_fcn_systems(d_num_parts),
      d_PK1_stress_fcn_ctxs(d_num_parts,NULL),
      d_lag_body_force_fcns(d_num_parts,NULL),
      d_lag_body_force_fcn_systems(d_num_parts),
      d_lag_body_force_fcn_ctxs(d_num_parts,NULL),
      d_lag_pressure_fcns(d_num_parts,NULL),
      d_lag_pressure_fcn_systems(d_num_parts),
      d_lag_pressure_fcn_ctxs(d_num_parts,NULL),
      d_lag_surface_force_fcns(d_num_parts,NULL),
      d_lag_surface_force_fcn_systems(d_num_parts),
      d_lag_surface_force_fcn_ctxs(d_num_parts,NULL),
      d_l_data_manager(l_data_manager),
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
    commonConstructor(input_db);
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
        delete it->second;
    }

    for (CoarsenPatchStrategyMap::iterator it = d_cstrategies.begin(); it != d_cstrategies.end(); ++it)
    {
        delete it->second;
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
    CoordinateMappingFcnPtr coordinate_mapping_fcn,
    void* coordinate_mapping_fcn_ctx,
    const unsigned int mesh_number)
{
    d_coordinate_mapping_fcns    [mesh_number] = coordinate_mapping_fcn;
    d_coordinate_mapping_fcn_ctxs[mesh_number] = coordinate_mapping_fcn_ctx;
    return;
}// registerInitialCoordinateMappingFunction

void
IBFEHierarchyIntegrator::registerPK1StressTensorFunction(
    PK1StressFcnPtr PK1_stress_fcn,
    std::vector<unsigned int> PK1_stress_fcn_systems,
    void* PK1_stress_fcn_ctx,
    const unsigned int mesh_number)
{
    d_PK1_stress_fcns       [mesh_number] = PK1_stress_fcn;
    d_PK1_stress_fcn_systems[mesh_number] = PK1_stress_fcn_systems;
    d_PK1_stress_fcn_ctxs   [mesh_number] = PK1_stress_fcn_ctx;
    return;
}// registerPK1StressTensorFunction

void
IBFEHierarchyIntegrator::registerLagBodyForceFunction(
    LagBodyForceFcnPtr lag_body_force_fcn,
    std::vector<unsigned int> lag_body_force_fcn_systems,
    void* lag_body_force_fcn_ctx,
    const unsigned int mesh_number)
{
    d_lag_body_force_fcns       [mesh_number] = lag_body_force_fcn;
    d_lag_body_force_fcn_systems[mesh_number] = lag_body_force_fcn_systems;
    d_lag_body_force_fcn_ctxs   [mesh_number] = lag_body_force_fcn_ctx;
    return;
}// registerLagBodyForceFunction

void
IBFEHierarchyIntegrator::registerLagPressureFunction(
    LagPressureFcnPtr lag_pressure_fcn,
    std::vector<unsigned int> lag_pressure_fcn_systems,
    void* lag_pressure_fcn_ctx,
    const unsigned int mesh_number)
{
    d_lag_pressure_fcns       [mesh_number] = lag_pressure_fcn;
    d_lag_pressure_fcn_systems[mesh_number] = lag_pressure_fcn_systems;
    d_lag_pressure_fcn_ctxs   [mesh_number] = lag_pressure_fcn_ctx;
    return;
}// registerLagPressureFunction

void
IBFEHierarchyIntegrator::registerLagSurfaceForceFunction(
    LagSurfaceForceFcnPtr lag_surface_force_fcn,
    std::vector<unsigned int> lag_surface_force_fcn_systems,
    void* lag_surface_force_fcn_ctx,
    const unsigned int mesh_number)
{
    d_lag_surface_force_fcns       [mesh_number] = lag_surface_force_fcn;
    d_lag_surface_force_fcn_systems[mesh_number] = lag_surface_force_fcn_systems;
    d_lag_surface_force_fcn_ctxs   [mesh_number] = lag_surface_force_fcn_ctx;
    return;
}// registerLagSurfaceForceFunction

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
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& U_bc_coefs)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::registerVelocityPhysicalBcCoefs()\n"
                   << "  velocity boundary conditions must be registered prior to initialization\n"
                   << "  of the hierarchy integrator object.\n");
    }
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
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->registerLoadBalancer(d_load_balancer);
    }
    if (d_l_data_manager != NULL) d_l_data_manager->registerLoadBalancer(d_load_balancer);
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
    IBAMR_TIMER_START(t_initialize_hierarchy_integrator);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Initialize FE system data.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();

        d_fe_data_managers[part]->COORDINATES_SYSTEM_NAME = COORDS_SYSTEM_NAME;
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

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
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

    IBAMR_TIMER_STOP(t_initialize_hierarchy_integrator);
    return;
}// initializeHierarchyIntegrator

double
IBFEHierarchyIntegrator::initializeHierarchy()
{
    IBAMR_TIMER_START(t_initialize_hierarchy);

    // Initialize the data structures and state variables for the FE equation
    // systems.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
        equation_systems->init();
        initializeCoordinates(part);
        updateCoordinateMapping(part);

        // Set up boundary conditions.  Specifically, add appropriate boundary
        // IDs to the BoundaryInfo object associated with the mesh, and add DOF
        // constraints for the nodal forces and velocities.
        System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
        X_system.assemble_before_solve = false;
        X_system.assemble();

        System& X_mapping_system = equation_systems->get_system(COORD_MAPPING_SYSTEM_NAME);
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
    }

    // Use the INSStaggeredHierarchyIntegrator to initialize the patch
    // hierarchy.
    double dt_next = d_ins_hier_integrator->initializeHierarchy();
    if (d_integrator_time >= d_dt_max_time_min && d_integrator_time <= d_dt_max_time_max)
    {
        dt_next = std::min(dt_next, d_dt_max);
    }

    // Reset the Lagrangian data manager.
    if (d_l_data_manager != NULL)
    {
        d_l_data_manager->beginDataRedistribution();
        d_l_data_manager->endDataRedistribution();
    }

    // Initialize the FE data manager.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->reinitElementMappings();
    }

    // Prune duplicate markers following initialization.
    LMarkerUtilities::pruneDuplicateMarkers(d_mark_current_idx, d_hierarchy);

    // Indicate that the force strategy needs to be re-initialized.
    d_ib_lag_force_strategy_needs_init = true;

    IBAMR_TIMER_STOP(t_initialize_hierarchy);
    return dt_next;
}// initializeHierarchy

double
IBFEHierarchyIntegrator::advanceHierarchy(
    const double dt)
{
    IBAMR_TIMER_START(t_advance_hierarchy);
    IBAMR_TIMER_START(t_advance_hierarchy_init);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
#endif

    PetscErrorCode ierr;

    const double current_time = d_integrator_time;
    const double new_time     = d_integrator_time+dt;
    const bool initial_time   = MathUtilities<double>::equalEps(current_time,d_start_time);

    if (atRegridPoint())
    {
        // Regrid the patch hierarchy.
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
            if (d_l_data_manager != NULL && d_l_data_manager->levelContainsLagrangianData(ln))
            {
                d_ib_lag_force_strategy->initializeLevelData(d_hierarchy, ln, current_time, initial_time, d_l_data_manager);
            }
        }
        d_ib_lag_force_strategy_needs_init = false;
    }

    // Allocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_V_idx, current_time);
        level->allocatePatchData(d_F_idx, current_time);
    }

    // Extract the FE vectors.
    std::vector<EquationSystems*> equation_systems(d_num_parts);
    std::vector<System*> X_systems(d_num_parts), F_systems(d_num_parts), U_systems(d_num_parts);
    std::vector<NumericVector<double>*> X_current_vecs(d_num_parts), X_half_vecs(d_num_parts), X_new_vecs(d_num_parts), U_half_vecs(d_num_parts), F_half_vecs(d_num_parts);
    std::vector<NumericVector<double>*> X_half_IB_ghost_vecs(d_num_parts), F_half_IB_ghost_vecs(d_num_parts);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        equation_systems[part] = d_fe_data_managers[part]->getEquationSystems();

        X_systems[part] = &equation_systems[part]->get_system(  COORDS_SYSTEM_NAME);
        F_systems[part] = &equation_systems[part]->get_system(   FORCE_SYSTEM_NAME);
        U_systems[part] = &equation_systems[part]->get_system(VELOCITY_SYSTEM_NAME);

        X_current_vecs[part] = X_systems[part]->solution.get();
        X_half_vecs   [part] = X_systems[part]->current_local_solution.get();
        X_new_vecs    [part] = X_current_vecs[part]->clone().release();  // NOTE: must be manually deleted
        U_half_vecs   [part] = U_systems[part]->solution.get();
        F_half_vecs   [part] = F_systems[part]->solution.get();

        X_half_IB_ghost_vecs[part] = d_fe_data_managers[part]->buildGhostedCoordsVector();
        F_half_IB_ghost_vecs[part] = d_fe_data_managers[part]->buildGhostedSolutionVector(FORCE_SYSTEM_NAME);
    }

    // Initialize the various LData objects on each level of the patch hierarchy.
    std::vector<Pointer<LData> > X_current_data(finest_ln+1);
    std::vector<Pointer<LData> > U_current_data(finest_ln+1);
    std::vector<Pointer<LData> > X_new_data    (finest_ln+1);
    std::vector<Pointer<LData> > X_half_data   (finest_ln+1);
    std::vector<Pointer<LData> > U_half_data   (finest_ln+1);
    std::vector<Pointer<LData> > F_half_data   (finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager != NULL && d_l_data_manager->levelContainsLagrangianData(ln))
        {
            X_current_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
            U_current_data[ln] = d_l_data_manager->getLData(LDataManager:: VEL_DATA_NAME,ln);
            X_new_data    [ln] = d_l_data_manager->createLData("X_new" ,ln,NDIM);
            X_half_data   [ln] = d_l_data_manager->createLData("X_half",ln,NDIM);
            U_half_data   [ln] = d_l_data_manager->createLData("U_half",ln,NDIM);
            F_half_data   [ln] = d_l_data_manager->createLData("F_half",ln,NDIM);
        }
    }

    // Get patch data descriptors for the current and new velocity data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getCurrentContext());
    const int U_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getNewContext());

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
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        ierr = VecCopy(dynamic_cast<PetscVector<double>*>(X_current_vecs[part])->vec(), dynamic_cast<PetscVector<double>*>(X_new_vecs[part])->vec()); IBTK_CHKERRQ(ierr);
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager != NULL && d_l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_current_vec = X_current_data[ln]->getVec();
            Vec X_new_vec = X_new_data[ln]->getVec();
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
            Pointer<LMarkerSetData> mark_data         = patch->getPatchData(d_mark_current_idx);
            Pointer<LMarkerSetData> mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);
            mark_scratch_data->copy(*mark_data);
        }
    }

    IBAMR_TIMER_STOP(t_advance_hierarchy_init);

    // Perform one or more cycles to compute the updated configuration of the
    // coupled fluid-structure system.
    d_ins_hier_integrator->integrateHierarchy_initialize(current_time, new_time);
    for (int cycle = 0; cycle < d_num_cycles; ++cycle)
    {
        // Set X(n+1/2) = 0.5*(X(n) + X(n+1)).
        IBAMR_TIMER_START(t_advance_hierarchy_phase1);
        for (unsigned part = 0; part < d_num_parts; ++part)
        {
            ierr = VecAXPBYPCZ(dynamic_cast<PetscVector<double>*>(X_half_vecs   [part])->vec(),
                               0.5, 0.5, 0.0,
                               dynamic_cast<PetscVector<double>*>(X_current_vecs[part])->vec(),
                               dynamic_cast<PetscVector<double>*>(X_new_vecs    [part])->vec()); IBTK_CHKERRQ(ierr);
            X_half_vecs[part]->close();
        }
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager != NULL && d_l_data_manager->levelContainsLagrangianData(ln))
            {
                Vec X_current_vec = X_current_data[ln]->getVec();
                Vec X_new_vec = X_new_data[ln]->getVec();
                Vec X_half_vec = X_half_data[ln]->getVec();
                ierr = VecAXPBYPCZ(X_half_vec, 0.5, 0.5, 0.0, X_current_vec, X_new_vec); IBTK_CHKERRQ(ierr);
            }
        }
        IBAMR_TIMER_STOP(t_advance_hierarchy_phase1);

        // Compute F(n+1/2) = F(X(n+1/2),t(n+1/2)).
        IBAMR_TIMER_START(t_advance_hierarchy_phase3);
        for (unsigned part = 0; part < d_num_parts; ++part)
        {
            computeInteriorForceDensity(*F_half_vecs[part], *X_half_vecs[part], current_time+0.5*dt, part);
        }
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager != NULL && d_l_data_manager->levelContainsLagrangianData(ln))
            {
                Vec F_half_vec = F_half_data[ln]->getVec();
                ierr = VecSet(F_half_vec, 0.0); IBTK_CHKERRQ(ierr);
                d_ib_lag_force_strategy->computeLagrangianForce(F_half_data[ln], X_half_data[ln], U_half_data[ln], d_hierarchy, ln, current_time+0.5*dt, d_l_data_manager);
            }
        }
        IBAMR_TIMER_STOP(t_advance_hierarchy_phase3);

        // Copy data into the "IB ghosted" vectors.
        IBAMR_TIMER_START(t_advance_hierarchy_phase4);
        for (unsigned int part = 0; part < d_num_parts; ++part)
        {
            ierr = VecCopy(dynamic_cast<PetscVector<double>*>(X_half_vecs[part])->vec(), dynamic_cast<PetscVector<double>*>(X_half_IB_ghost_vecs[part])->vec()); IBTK_CHKERRQ(ierr);
            ierr = VecCopy(dynamic_cast<PetscVector<double>*>(F_half_vecs[part])->vec(), dynamic_cast<PetscVector<double>*>(F_half_IB_ghost_vecs[part])->vec()); IBTK_CHKERRQ(ierr);
            X_half_IB_ghost_vecs[part]->close();
            F_half_IB_ghost_vecs[part]->close();
        }
        IBAMR_TIMER_STOP(t_advance_hierarchy_phase4);

        // Spread F(n+1/2) to f(n+1/2).
        IBAMR_TIMER_START(t_advance_hierarchy_phase5);
        d_hier_sc_data_ops->setToScalar(d_F_idx, 0.0);
        for (unsigned int part = 0; part < d_num_parts; ++part)
        {
            NumericVector<double>& X_half_IB_ghost = *X_half_IB_ghost_vecs[part];
            NumericVector<double>& F_half_IB_ghost = *F_half_IB_ghost_vecs[part];
            if (d_use_IB_spreading_operator)
            {
                d_fe_data_managers[part]->spread(d_F_idx, F_half_IB_ghost, X_half_IB_ghost, FORCE_SYSTEM_NAME, false, false);
            }
            else
            {
                d_fe_data_managers[part]->prolongDensity(d_F_idx, F_half_IB_ghost, X_half_IB_ghost, FORCE_SYSTEM_NAME, false, false);
            }
            if (d_split_interior_and_bdry_forces)
            {
                if (d_use_jump_conditions)
                {
                    imposeJumpConditions(d_F_idx, F_half_IB_ghost, X_half_IB_ghost, current_time+0.5*dt, part);
                }
                else
                {
                    spreadTransmissionForceDensity(d_F_idx, X_half_IB_ghost, current_time+0.5*dt, part);
                }
            }
        }
        if (d_l_data_manager != NULL)
        {
            d_l_data_manager->spread(d_F_idx, F_half_data, X_half_data);
        }
        IBAMR_TIMER_STOP(t_advance_hierarchy_phase5);

        // Solve the incompressible Navier-Stokes equations.
        IBAMR_TIMER_START(t_advance_hierarchy_phase6);
        d_ins_hier_integrator->integrateHierarchy(current_time, new_time, cycle);
        IBAMR_TIMER_STOP(t_advance_hierarchy_phase6);

        // Set u(n+1/2) = 0.5*(u(n) + u(n+1)) and interpolate u(n+1/2) to
        // U(n+1/2).
        IBAMR_TIMER_START(t_advance_hierarchy_phase7);
        d_hier_sc_data_ops->linearSum(d_V_idx, 0.5, U_current_idx, 0.5, U_new_idx);
        for (unsigned int part = 0; part < d_num_parts; ++part)
        {
            NumericVector<double>& X_half_IB_ghost = *X_half_IB_ghost_vecs[part];
            NumericVector<double>& U_half = *U_half_vecs[part];
            if (d_use_IB_interpolation_operator)
            {
                d_fe_data_managers[part]->interp(d_V_idx, U_half, X_half_IB_ghost, VELOCITY_SYSTEM_NAME, d_rscheds["V->V::S->S::GHOST_FILL"], current_time, false);
            }
            else
            {
                d_fe_data_managers[part]->restrictValue(d_V_idx, U_half, X_half_IB_ghost, VELOCITY_SYSTEM_NAME, false);
            }
        }
        if (d_l_data_manager != NULL)
        {
            d_l_data_manager->interp(d_V_idx, U_half_data, X_half_data, std::vector<Pointer<CoarsenSchedule<NDIM> > >(), std::vector<Pointer<RefineSchedule<NDIM> > >(), current_time, false);
        }
        IBAMR_TIMER_STOP(t_advance_hierarchy_phase7);

        // Set X(n+1) = X(n) + dt*U(n+1/2).
        IBAMR_TIMER_START(t_advance_hierarchy_phase8);
        for (unsigned int part = 0; part < d_num_parts; ++part)
        {
            ierr = VecWAXPY(dynamic_cast<PetscVector<double>*>(X_new_vecs    [part])->vec(),
                            dt,
                            dynamic_cast<PetscVector<double>*>(U_half_vecs   [part])->vec(),
                            dynamic_cast<PetscVector<double>*>(X_current_vecs[part])->vec()); IBTK_CHKERRQ(ierr);
        }
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_l_data_manager != NULL && d_l_data_manager->levelContainsLagrangianData(ln))
            {
                Vec X_current_vec = X_current_data[ln]->getVec();
                Vec X_new_vec = X_new_data[ln]->getVec();
                Vec U_half_vec = U_half_data[ln]->getVec();
                ierr = VecWAXPY(X_new_vec, dt, U_half_vec, X_current_vec); IBTK_CHKERRQ(ierr);
            }
        }
        LMarkerUtilities::advectMarkers(d_mark_current_idx, d_mark_scratch_idx, d_V_idx, dt, d_fe_data_managers[0]->getInterpWeightingFunction(), d_hierarchy);
        IBAMR_TIMER_STOP(t_advance_hierarchy_phase8);
    }
    d_ins_hier_integrator->integrateHierarchy_finalize(current_time, new_time);

    IBAMR_TIMER_START(t_advance_hierarchy_finalize);

    // Reset X_current to equal X_new.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecCopy(dynamic_cast<PetscVector<double>*>(X_new_vecs[part])->vec(), dynamic_cast<PetscVector<double>*>(X_current_vecs[part])->vec()); IBTK_CHKERRQ(ierr);
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager != NULL && d_l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_current_vec = X_current_data[ln]->getVec();
            Vec X_new_vec = X_new_data[ln]->getVec();
            ierr = VecCopy(X_new_vec, X_current_vec); IBTK_CHKERRQ(ierr);
        }
    }

    // Update the coordinate mapping dX = X - s.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        updateCoordinateMapping(part);
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

    // Update the interpolated velocity to corresond to the updated Eulerian
    // velocity field
    d_hier_sc_data_ops->copyData(d_V_idx, U_new_idx);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecCopy(dynamic_cast<PetscVector<double>*>(X_new_vecs[part])->vec(), dynamic_cast<PetscVector<double>*>(X_half_IB_ghost_vecs[part])->vec()); IBTK_CHKERRQ(ierr);
        NumericVector<double>& X_new_IB_ghost = *X_half_IB_ghost_vecs[part];
        NumericVector<double>& U_new = *U_half_vecs[part];
        if (d_use_IB_interpolation_operator)
        {
            d_fe_data_managers[part]->interp(d_V_idx, U_new, X_new_IB_ghost, VELOCITY_SYSTEM_NAME, d_rscheds["V->V::S->S::GHOST_FILL"], current_time, true);
        }
        else
        {
            d_fe_data_managers[part]->restrictValue(d_V_idx, U_new, X_new_IB_ghost, VELOCITY_SYSTEM_NAME, true);
        }
    }

    // Deallocate scratch data.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        delete X_new_vecs[part];
    }
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
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): current CFL number = " << cfl_max << "\n";
    if (d_do_log) plog << d_object_name << "::advanceHierarchy(): estimated upper bound on IB point displacement since last regrid = " << d_regrid_cfl_estimate << "\n";

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

    IBAMR_TIMER_STOP(t_advance_hierarchy_finalize);
    IBAMR_TIMER_STOP(t_advance_hierarchy);
    return dt_next;
}// advanceHierarchy

bool
IBFEHierarchyIntegrator::atRegridPoint() const
{
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time,d_start_time);
    if (initial_time) return true;
    return d_regrid_cfl_interval > 0.0 ? (d_regrid_cfl_estimate >= d_regrid_cfl_interval) : (d_regrid_interval == 0 ? false : (d_integrator_step % d_regrid_interval == 0));
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
    IBAMR_TIMER_START(t_regrid_hierarchy);

    // Update the marker data.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): resetting marker particles.\n";
    LMarkerUtilities::collectMarkersOnPatchHierarchy(d_mark_current_idx, d_hierarchy);

    // Update the workload data.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());
    }
    if (d_l_data_manager != NULL)
    {
        d_l_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());
    }

    // Begin Lagrangian data redistribution.
    if (d_l_data_manager != NULL)
    {
        d_l_data_manager->beginDataRedistribution();
    }

    // Use the INSStaggeredHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->regridHierarchy();

    // Complete Lagrangian data redistribution.
    if (d_l_data_manager != NULL)
    {
        d_l_data_manager->endDataRedistribution();
    }
    d_ib_lag_force_strategy_needs_init = true;

    // Reinitialize the FE data manager.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->reinitElementMappings();
    }

    // Prune duplicate markers following regridding.
    LMarkerUtilities::pruneDuplicateMarkers(d_mark_current_idx, d_hierarchy);

    // Reset the regrid CFL estimate.
    d_regrid_cfl_estimate = 0.0;

    IBAMR_TIMER_STOP(t_regrid_hierarchy);
    return;
}// regridHierarchy

void
IBFEHierarchyIntegrator::synchronizeHierarchy()
{
    IBAMR_TIMER_START(t_synchronize_hierarchy);

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->synchronizeHierarchy();

    IBAMR_TIMER_STOP(t_synchronize_hierarchy);
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
    IBAMR_TIMER_START(t_synchronize_new_levels);

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

    IBAMR_TIMER_STOP(t_synchronize_new_levels);
    return;
}// synchronizeNewLevels

void
IBFEHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    IBAMR_TIMER_START(t_reset_time_dependent_data);

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->resetTimeDependentHierData(new_time);

    IBAMR_TIMER_STOP(t_reset_time_dependent_data);
    return;
}// resetTimeDependentHierData

void
IBFEHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    IBAMR_TIMER_START(t_reset_data_to_preadvance_state);

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
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
IBFEHierarchyIntegrator::initializeLevelData(
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
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Initialize the FE data manager.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->setPatchHierarchy(hierarchy);
        d_fe_data_managers[part]->resetLevels(0,hierarchy->getFinestLevelNumber());
        d_fe_data_managers[part]->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    }

    // Initialize the Lagrangian data manager.
    if (d_l_data_manager != NULL)
    {
        d_l_data_manager->setPatchHierarchy(hierarchy);
        d_l_data_manager->resetLevels(0,hierarchy->getFinestLevelNumber());
        d_l_data_manager->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
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
    LMarkerUtilities::initializeMarkersOnLevel(d_mark_current_idx, d_mark_init_posns, hierarchy, level_number, initial_time, old_level);

    IBAMR_TIMER_STOP(t_initialize_level_data);
    return;
}// initializeLevelData

void
IBFEHierarchyIntegrator::resetHierarchyConfiguration(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    IBAMR_TIMER_START(t_reset_hierarchy_configuration);

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
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_hier_level);
    }

    // Reset the Lagrangian data manager.
    if (d_l_data_manager != NULL)
    {
        d_l_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_hier_level);
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
        d_rscheds[it->first].resize(finest_hier_level+1);
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

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin(); it!= d_calgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(ln-1);
            d_cscheds[it->first][ln] = it->second->createSchedule(coarser_level, level, d_cstrategies[it->first]);
        }
    }

    IBAMR_TIMER_STOP(t_reset_hierarchy_configuration);
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
    IBAMR_TIMER_START(t_apply_gradient_detector);

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
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }
    if (d_l_data_manager != NULL)
    {
        d_l_data_manager->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }

    IBAMR_TIMER_STOP(t_apply_gradient_detector);
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
IBFEHierarchyIntegrator::getLMarkerSetVar() const
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
    IBAMR_TIMER_START(t_put_to_database);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("IB_FE_HIERARCHY_INTEGRATOR_VERSION", IB_FE_HIERARCHY_INTEGRATOR_VERSION);

    db->putBool("d_use_IB_spreading_operator", d_use_IB_spreading_operator);
    db->putBool("d_use_IB_interpolation_operator", d_use_IB_interpolation_operator);
    db->putBool("d_split_interior_and_bdry_forces", d_split_interior_and_bdry_forces);
    db->putBool("d_use_jump_conditions", d_use_jump_conditions);
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putDouble("d_regrid_cfl_interval", d_regrid_cfl_interval);
    db->putDouble("d_regrid_cfl_estimate", d_regrid_cfl_estimate);
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
IBFEHierarchyIntegrator::commonConstructor(
    Pointer<Database> input_db)
{

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Determine whether we should use first-order or second-order shape
    // functions for each part of the structure.
    bool mesh_has_first_order_elems  = false;
    bool mesh_has_second_order_elems = false;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
        const MeshBase& mesh = equation_systems->get_mesh();
        MeshBase::const_element_iterator       el_it  = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for ( ; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            mesh_has_first_order_elems  = mesh_has_first_order_elems  || elem->default_order() == FIRST ;
            mesh_has_second_order_elems = mesh_has_second_order_elems || elem->default_order() == SECOND;
        }
    }
    mesh_has_first_order_elems  = SAMRAI_MPI::maxReduction(mesh_has_first_order_elems );
    mesh_has_second_order_elems = SAMRAI_MPI::maxReduction(mesh_has_second_order_elems);
    if (( mesh_has_first_order_elems &&  mesh_has_second_order_elems) ||
        (!mesh_has_first_order_elems && !mesh_has_second_order_elems))
    {
        TBOX_ERROR(d_object_name << "::IBFEHierarchyIntegrator():\n"
                   << "  a parts of FE mesh must contain only FIRST order elements or only SECOND order elements" << std::endl);
    }
    if (mesh_has_first_order_elems )
    {
        d_fe_order = FIRST;
        d_quad_order = THIRD;
    }
    if (mesh_has_second_order_elems)
    {
        d_fe_order = SECOND;
        d_quad_order = FIFTH;
    }
    pout << d_object_name << ": using " << Utility::enum_to_string<Order>(d_fe_order) << " order " << Utility::enum_to_string<FEFamily>(d_fe_family) << " finite elements.\n\n";

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    getFromInput(input_db, from_restart);

    // Read in the marker initial positions.
    if (!from_restart) d_num_mark = LMarkerUtilities::readMarkerPositions(d_mark_init_posns, d_mark_input_file_name, d_hierarchy);

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<SideVariable<NDIM,double> > sc_var = new SideVariable<NDIM,double>("sc_var");
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, d_hierarchy, true);

    // Initialize all variable contexts.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_current = var_db->getContext(d_object_name+"::CURRENT");
    d_scratch = var_db->getContext(d_object_name+"::SCRATCH");

    // Initialize all variables.
    IntVector<NDIM> ghosts = 0;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ghosts = IntVector<NDIM>::max(ghosts,d_fe_data_managers[part]->getGhostCellWidth());
    }
    const IntVector<NDIM> no_ghosts = 0;

    d_V_var = new SideVariable<NDIM,double>(d_object_name+"::V");
    d_V_idx = var_db->registerVariableAndContext(d_V_var, d_scratch, ghosts);

    d_F_var = new SideVariable<NDIM,double>(d_object_name+"::F");
    d_F_idx = var_db->registerVariableAndContext(d_F_var, d_scratch, no_ghosts);

    d_mark_var = new LMarkerSetVariable(d_object_name+"::mark");
    d_mark_current_idx = var_db->registerVariableAndContext(d_mark_var, getCurrentContext(), ghosts);
    d_mark_scratch_idx = var_db->registerVariableAndContext(d_mark_var, getScratchContext(), ghosts);
    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_mark_current_idx);
    }

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
}// commonConstructor

void
IBFEHierarchyIntegrator::computeInteriorForceDensity(
    NumericVector<double>& G_vec,
    NumericVector<double>& X_vec,
    const double& time,
    const unsigned int part)
{
    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const int dim = mesh.mesh_dimension();
    AutoPtr<QBase> qrule = QBase::build(d_quad_type, dim, d_quad_order);
    AutoPtr<QBase> qrule_face = QBase::build(d_quad_type, dim-1, d_quad_order);

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& system = equation_systems->get_system(FORCE_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d) TBOX_ASSERT(dof_map.variable_type(d) == dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> dof_indices(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) dof_indices(d).reserve(27);
    AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<Point>& q_point = fe->get_xyz();
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi = fe->get_dphi();
    AutoPtr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
    fe_face->attach_quadrature_rule(qrule_face.get());
    const std::vector<Point>& q_point_face = fe_face->get_xyz();
    const std::vector<double>& JxW_face = fe_face->get_JxW();
    const std::vector<Point>& normal_face = fe_face->get_normals();
    const std::vector<std::vector<double> >& phi_face = fe_face->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi_face = fe_face->get_dphi();

#ifdef DEBUG_CHECK_ASSERTIONS
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    for (unsigned int d = 0; d < NDIM; ++d) TBOX_ASSERT(dof_map.variable_type(d) == X_dof_map.variable_type(d));
#endif

    // Setup extra data needed to compute stresses/forces.
    std::vector<NumericVector<double>*> PK1_stress_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_PK1_stress_fcn_systems[part].begin();
         cit != d_PK1_stress_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        system.update();
        PK1_stress_fcn_data.push_back(system.current_local_solution.get());
    }

    std::vector<NumericVector<double>*> lag_body_force_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_lag_body_force_fcn_systems[part].begin();
         cit != d_lag_body_force_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        system.update();
        lag_body_force_fcn_data.push_back(system.current_local_solution.get());
    }

    std::vector<NumericVector<double>*> lag_pressure_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_lag_pressure_fcn_systems[part].begin();
         cit != d_lag_pressure_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        system.update();
        lag_pressure_fcn_data.push_back(system.current_local_solution.get());
    }

    std::vector<NumericVector<double>*> lag_surface_force_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_lag_surface_force_fcn_systems[part].begin();
         cit != d_lag_surface_force_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        system.update();
        lag_surface_force_fcn_data.push_back(system.current_local_solution.get());
    }

    // Setup global and elemental right-hand-side vectors.
    AutoPtr<NumericVector<double> > G_rhs_vec = G_vec.zero_clone();
    DenseVector<double> G_rhs_e[NDIM];

    // Loop over the elements to compute the right-hand side vector.  This is
    // computed via
    //
    //    rhs_k = -int{PP(s,t) grad phi_k(s)}ds + int{PP(s,t) N(s,t) phi_k(s)}dA(s)
    //
    // This right-hand side vector is used to solve for the nodal values of the
    // interior elastic force density.
    TensorValue<double> PP, dX_ds, dX_ds_inv_trans;
    VectorValue<double> F, F_b, F_s, F_qp, n;
    Point X_qp;
    double P;
    blitz::Array<double,2> X_node;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;

        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices(d), d);
            if (G_rhs_e[d].size() != dof_indices(d).size())
            {
                G_rhs_e[d].resize(dof_indices(d).size());  // NOTE: DenseVector::resize() automatically zeroes the vector contents.
            }
            else
            {
                G_rhs_e[d].zero();
            }
        }

        const unsigned int n_qp = qrule->n_points();
        const unsigned int n_basis = dof_indices(0).size();

        get_values_for_interpolation(X_node, X_vec, dof_indices);
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            const Point& s_qp = q_point[qp];
            interpolate(X_qp,qp,X_node,phi);
            jacobian(dX_ds,qp,X_node,dphi);

            // Compute the value of the first Piola-Kirchhoff stress tensor at
            // the quadrature point and add the corresponding forces to the
            // right-hand-side vector.
            d_PK1_stress_fcns[part](PP,dX_ds,X_qp,s_qp,elem,X_vec,PK1_stress_fcn_data,time,d_PK1_stress_fcn_ctxs[part]);
            for (unsigned int k = 0; k < n_basis; ++k)
            {
                F_qp = -PP*dphi[k][qp]*JxW[qp];
                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    G_rhs_e[i](k) += F_qp(i);
                }
            }

            if (d_lag_body_force_fcns[part] != NULL)
            {
                // Compute the value of the body force at the quadrature point
                // and add the corresponding forces to the right-hand-side
                // vector.
                d_lag_body_force_fcns[part](F_b,dX_ds,X_qp,s_qp,elem,X_vec,lag_body_force_fcn_data,time,d_lag_body_force_fcn_ctxs[part]);
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = phi[k][qp]*JxW[k]*F_b;
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        G_rhs_e[i](k) += F_qp(i);
                    }
                }
            }
        }

        // Loop over the element boundaries.
        for (unsigned short int side = 0; side < elem->n_sides(); ++side)
        {
            // Determine whether we are at a physical boundary and, if so,
            // whether it is a Dirichlet boundary.
            bool at_physical_bdry = elem->neighbor(side) == NULL;
            bool at_dirichlet_bdry = false;
            const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
            for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
            {
                const short int bdry_id = *cit;
                at_physical_bdry  = at_physical_bdry  && !dof_map.is_periodic_boundary(bdry_id);
                at_dirichlet_bdry = at_dirichlet_bdry || (bdry_id == FEDataManager::DIRICHLET_BDRY_ID);
            }

            // Skip non-physical boundaries.
            if (!at_physical_bdry) continue;

            // Determine whether we need to compute surface forces along this
            // part of the physical boundary; if not, skip the present side.
            const bool compute_transmission_force = (( d_split_interior_and_bdry_forces && !at_dirichlet_bdry) ||
                                                     (!d_split_interior_and_bdry_forces &&  at_dirichlet_bdry));
            const bool compute_pressure = (!d_split_interior_and_bdry_forces && d_lag_pressure_fcns[part] != NULL && !at_dirichlet_bdry);
            const bool compute_surface_force = (!d_split_interior_and_bdry_forces && d_lag_surface_force_fcns[part] != NULL && !at_dirichlet_bdry);
            if (!(compute_transmission_force || compute_pressure || compute_surface_force)) continue;

            fe_face->reinit(elem, side);

            const unsigned int n_qp = qrule_face->n_points();
            const unsigned int n_basis = dof_indices(0).size();

            get_values_for_interpolation(X_node, X_vec, dof_indices);
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const Point& s_qp = q_point_face[qp];
                interpolate(X_qp,qp,X_node,phi_face);
                jacobian(dX_ds,qp,X_node,dphi_face);
                const double J = std::abs(dX_ds.det());
                tensor_inverse_transpose(dX_ds_inv_trans,dX_ds,NDIM);
                F.zero();
                if (compute_transmission_force)
                {
                    // Compute the value of the first Piola-Kirchhoff stress
                    // tensor at the quadrature point and add the corresponding
                    // force to the right-hand-side vector.
                    d_PK1_stress_fcns[part](PP,dX_ds,X_qp,s_qp,elem,X_vec,PK1_stress_fcn_data,time,d_PK1_stress_fcn_ctxs[part]);
                    F += PP*normal_face[qp];
                }
                if (compute_pressure)
                {
                    // Compute the value of the pressure at the quadrature point
                    // and add the corresponding force to the right-hand-side
                    // vector.
                    d_lag_pressure_fcns[part](P,dX_ds,X_qp,s_qp,elem,side,X_vec,lag_pressure_fcn_data,time,d_lag_pressure_fcn_ctxs[part]);
                    F -= P*J*dX_ds_inv_trans*normal_face[qp];
                }
                if (compute_surface_force)
                {
                    // Compute the value of the surface force at the quadrature
                    // point and add the corresponding force to the
                    // right-hand-side vector.
                    d_lag_surface_force_fcns[part](F_s,dX_ds,X_qp,s_qp,elem,side,X_vec,lag_surface_force_fcn_data,time,d_lag_surface_force_fcn_ctxs[part]);
                    F += F_s;
                }

                // If we are imposing jump conditions, then we keep only the
                // normal part of the force.  This has the effect of projecting
                // the tangential part of the surface force onto the interior
                // force density.
                if (d_split_interior_and_bdry_forces && d_use_jump_conditions && !at_dirichlet_bdry)
                {
                    n = (dX_ds_inv_trans*normal_face[qp]).unit();
                    F = (F*n)*n;
                }

                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = phi_face[k][qp]*JxW_face[qp]*F;
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        G_rhs_e[i](k) += F_qp(i);
                    }
                }
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions) and
        // add the elemental contributions to the global vector.
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            dof_map.constrain_element_vector(G_rhs_e[i], dof_indices(i));
            G_rhs_vec->add_vector(G_rhs_e[i], dof_indices(i));
        }
    }

    // Solve for G.
    G_rhs_vec->close();
    d_fe_data_managers[part]->computeL2Projection(G_vec, *G_rhs_vec, FORCE_SYSTEM_NAME, d_use_consistent_mass_matrix);
    return;
}// computeInteriorForceDensity

void
IBFEHierarchyIntegrator::spreadTransmissionForceDensity(
    const int f_data_idx,
    NumericVector<double>& X_ghost_vec,
    const double& time,
    const unsigned int part)
{
    if (!d_split_interior_and_bdry_forces) return;

    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const int dim = mesh.mesh_dimension();
    QAdaptiveGauss* ib_qrule = d_fe_data_managers[part]->getQuadratureRule();
    QAdaptiveGauss* ib_qrule_face = d_fe_data_managers[part]->getQuadratureRuleFace();

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& system = equation_systems->get_system(FORCE_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d) TBOX_ASSERT(dof_map.variable_type(d) == dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> dof_indices(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) dof_indices(d).reserve(27);
    AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    fe->attach_quadrature_rule(ib_qrule);
    fe->get_xyz();  // prevents FE::reinit() from rebuilding all FE data items...
    AutoPtr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
    fe_face->attach_quadrature_rule(ib_qrule_face);
    const std::vector<Point>& q_point_face = fe_face->get_xyz();
    const std::vector<double>& JxW_face = fe_face->get_JxW();
    const std::vector<Point>& normal_face = fe_face->get_normals();
    const std::vector<std::vector<double> >& phi_face = fe_face->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi_face = fe_face->get_dphi();

    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
#ifdef DEBUG_CHECK_ASSERTIONS
    const DofMap& X_dof_map = X_system.get_dof_map();
    for (unsigned int d = 0; d < NDIM; ++d) TBOX_ASSERT(dof_map.variable_type(d) == X_dof_map.variable_type(d));
#endif

    // Setup extra data needed to compute stresses/forces.
    std::vector<NumericVector<double>*> PK1_stress_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_PK1_stress_fcn_systems[part].begin();
         cit != d_PK1_stress_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        PK1_stress_fcn_data.push_back(d_fe_data_managers[part]->buildGhostedSolutionVector(system.name()));
    }

    std::vector<NumericVector<double>*> lag_pressure_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_lag_pressure_fcn_systems[part].begin();
         cit != d_lag_pressure_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        lag_pressure_fcn_data.push_back(d_fe_data_managers[part]->buildGhostedSolutionVector(system.name()));
    }

    std::vector<NumericVector<double>*> lag_surface_force_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_lag_surface_force_fcn_systems[part].begin();
         cit != d_lag_surface_force_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        lag_surface_force_fcn_data.push_back(d_fe_data_managers[part]->buildGhostedSolutionVector(system.name()));
    }

    // Loop over the patches to spread the transmission elastic force density
    // onto the grid.
    const blitz::Array<blitz::Array<Elem*,1>,1>& active_patch_element_map = d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getLevelNumber();
    TensorValue<double> PP, dX_ds, dX_ds_inv_trans;
    VectorValue<double> F, F_s;
    Point X_qp;
    double P;
    std::vector<Point> elem_X;
    blitz::Array<double,2> X_node;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const blitz::Array<Elem*,1>& patch_elems = active_patch_element_map(local_patch_num);
        const int num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();

        // Setup vectors to store the values of T and X at the quadrature
        // points.  We compute a conservative upper bound on the number of
        // quadrature points to try to avoid unnecessary reallocations.
        static const unsigned int n_qp_estimate = (NDIM == 2 ? 16 : 16*16);
        std::vector<double> T_bdry;
        T_bdry.reserve(NDIM*n_qp_estimate*num_active_patch_elems);
        std::vector<double> X_bdry;
        X_bdry.reserve(NDIM*n_qp_estimate*num_active_patch_elems);

        // Loop over the elements and compute the values to be spread and the
        // positions of the quadrature points.
        int qp_offset = 0;
        for (int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems(e_idx);

            bool has_physical_boundaries = false;
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                bool at_physical_bdry = elem->neighbor(side) == NULL;
                const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
                for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
                {
                    const short int bdry_id = *cit;
                    at_physical_bdry  = at_physical_bdry  && !dof_map.is_periodic_boundary(bdry_id);
                }
                has_physical_boundaries = has_physical_boundaries || at_physical_bdry;
            }
            if (!has_physical_boundaries) continue;

            get_nodal_positions(elem_X, elem, X_ghost_vec, X_system.number());
            ib_qrule->set_elem_data(&elem_X, patch_dx);
            fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dof_map.dof_indices(elem, dof_indices(d), d);
            }

            // Loop over the element boundaries.
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                // Determine whether we are at a physical boundary and, if so,
                // whether it is a Dirichlet boundary.
                bool at_physical_bdry = elem->neighbor(side) == NULL;
                bool at_dirichlet_bdry = false;
                const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
                for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
                {
                    const short int bdry_id = *cit;
                    at_physical_bdry  = at_physical_bdry  && !dof_map.is_periodic_boundary(bdry_id);
                    at_dirichlet_bdry = at_dirichlet_bdry || (bdry_id == FEDataManager::DIRICHLET_BDRY_ID);
                }

                // Skip non-physical boundaries.
                if (!at_physical_bdry) continue;

                // Determine whether we need to compute surface forces along
                // this part of the physical boundary; if not, skip the present
                // side.
                const bool compute_transmission_force = d_split_interior_and_bdry_forces && !at_dirichlet_bdry;
                const bool compute_pressure = (d_split_interior_and_bdry_forces && d_lag_pressure_fcns[part] != NULL);
                const bool compute_surface_force = (d_split_interior_and_bdry_forces && d_lag_surface_force_fcns[part] != NULL);
                if (!(compute_transmission_force || compute_pressure || compute_surface_force)) continue;

                AutoPtr<Elem> side_elem = elem->build_side(side);
                get_nodal_positions(elem_X, side_elem.get(), X_ghost_vec, X_system.number());
                ib_qrule_face->set_elem_data(&elem_X, patch_dx);
                fe_face->reinit(elem, side);

                const unsigned int n_qp = ib_qrule_face->n_points();

                T_bdry.resize(T_bdry.size()+NDIM*n_qp);
                X_bdry.resize(X_bdry.size()+NDIM*n_qp);

                get_values_for_interpolation(X_node, X_ghost_vec, dof_indices);
                for (unsigned int qp = 0; qp < n_qp; ++qp, ++qp_offset)
                {
                    const Point& s_qp = q_point_face[qp];
                    interpolate(X_qp,qp,X_node,phi_face);
                    jacobian(dX_ds,qp,X_node,dphi_face);
                    const double J = std::abs(dX_ds.det());
                    tensor_inverse_transpose(dX_ds_inv_trans,dX_ds,NDIM);
                    F.zero();
                    if (compute_transmission_force)
                    {
                        // Compute the value of the first Piola-Kirchhoff stress
                        // tensor at the quadrature point and compute the
                        // corresponding force.
                        d_PK1_stress_fcns[part](PP,dX_ds,X_qp,s_qp,elem,X_ghost_vec,PK1_stress_fcn_data,time,d_PK1_stress_fcn_ctxs[part]);
                        F -= PP*normal_face[qp]*JxW_face[qp];
                    }
                    if (compute_pressure)
                    {
                        // Compute the value of the pressure at the quadrature
                        // point and compute the corresponding force.
                        d_lag_pressure_fcns[part](P,dX_ds,X_qp,s_qp,elem,side,X_ghost_vec,lag_pressure_fcn_data,time,d_lag_pressure_fcn_ctxs[part]);
                        F -= P*J*dX_ds_inv_trans*normal_face[qp]*JxW_face[qp];
                    }
                    if (compute_surface_force)
                    {
                        // Compute the value of the surface force at the
                        // quadrature point and compute the corresponding force.
                        d_lag_surface_force_fcns[part](F_s,dX_ds,X_qp,s_qp,elem,side,X_ghost_vec,lag_surface_force_fcn_data,time,d_lag_surface_force_fcn_ctxs[part]);
                        F += F_s;
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

        if (qp_offset == 0) continue;

        // Spread the boundary forces to the grid.
        const std::string& spread_weighting_fcn = d_fe_data_managers[part]->getSpreadWeightingFunction();
        const hier::IntVector<NDIM>& ghost_width = d_fe_data_managers[part]->getGhostCellWidth();
        const Box<NDIM> spread_box = Box<NDIM>::grow(patch->getBox(), ghost_width);
        Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(f_data_idx);
        LEInteractor::spread(f_data, T_bdry, NDIM, X_bdry, NDIM, patch, spread_box, spread_weighting_fcn);
    }
    return;
}// spreadTransmissionForceDensity

void
IBFEHierarchyIntegrator::imposeJumpConditions(
    const int f_data_idx,
    NumericVector<double>& F_ghost_vec,
    NumericVector<double>& X_ghost_vec,
    const double& time,
    const unsigned int part)
{
    if (!d_split_interior_and_bdry_forces) return;

    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& system = equation_systems->get_system(FORCE_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d) TBOX_ASSERT(dof_map.variable_type(d) == dof_map.variable_type(0));
#endif
    FEType fe_type = dof_map.variable_type(0);
    blitz::Array<std::vector<unsigned int>,1> dof_indices(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) dof_indices(d).reserve(27);
    blitz::Array<std::vector<unsigned int>,1> side_dof_indices(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) side_dof_indices(d).reserve(9);
    AutoPtr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
    const std::vector<Point>& q_point_face = fe_face->get_xyz();
    const std::vector<Point>& normal_face = fe_face->get_normals();
    const std::vector<std::vector<double> >& phi_face = fe_face->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi_face = fe_face->get_dphi();

#ifdef DEBUG_CHECK_ASSERTIONS
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    for (unsigned int d = 0; d < NDIM; ++d) TBOX_ASSERT(dof_map.variable_type(d) == X_dof_map.variable_type(d));
#endif

    // Setup extra data needed to compute stresses/forces.
    std::vector<NumericVector<double>*> PK1_stress_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_PK1_stress_fcn_systems[part].begin();
         cit != d_PK1_stress_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        PK1_stress_fcn_data.push_back(d_fe_data_managers[part]->buildGhostedSolutionVector(system.name()));
    }

    std::vector<NumericVector<double>*> lag_pressure_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_lag_pressure_fcn_systems[part].begin();
         cit != d_lag_pressure_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        lag_pressure_fcn_data.push_back(d_fe_data_managers[part]->buildGhostedSolutionVector(system.name()));
    }

    std::vector<NumericVector<double>*> lag_surface_force_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_lag_surface_force_fcn_systems[part].begin();
         cit != d_lag_surface_force_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        lag_surface_force_fcn_data.push_back(d_fe_data_managers[part]->buildGhostedSolutionVector(system.name()));
    }

    // Loop over the patches to impose jump conditions on the Eulerian grid that
    // are determined from the interior and transmission elastic force
    // densities.
    const blitz::Array<blitz::Array<Elem*,1>,1>& active_patch_element_map = d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getLevelNumber();
    TensorValue<double> PP, dX_ds, dX_ds_inv_trans;
    VectorValue<double> F, F_s, F_qp, n;
    Point X_qp;
    double P;
    blitz::Array<double,2> F_node, X_node;
    static const unsigned int MAX_NODES = (NDIM == 2 ? 9 : 27);
    Point s_node_cache[MAX_NODES], X_node_cache[MAX_NODES];
    blitz::TinyVector<double,NDIM> X_min, X_max;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const blitz::Array<Elem*,1>& patch_elems = active_patch_element_map(local_patch_num);
        const int num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(f_data_idx);
        const Box<NDIM>& patch_box = patch->getBox();
        const CellIndex<NDIM>& patch_lower = patch_box.lower();
        const CellIndex<NDIM>& patch_upper = patch_box.upper();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const x_lower = patch_geom->getXLower();
        const double* const x_upper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();

        blitz::TinyVector<Box<NDIM>,NDIM> side_boxes;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box,axis);
        }

        // Loop over the elements.
        for (int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems(e_idx);

            bool has_physical_boundaries = false;
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                bool at_physical_bdry = elem->neighbor(side) == NULL;
                const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
                for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
                {
                    const short int bdry_id = *cit;
                    at_physical_bdry  = at_physical_bdry  && !dof_map.is_periodic_boundary(bdry_id);
                }
                has_physical_boundaries = has_physical_boundaries || at_physical_bdry;
            }
            if (!has_physical_boundaries) continue;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dof_map.dof_indices(elem, dof_indices(d), d);
            }

            // Loop over the element boundaries.
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                // Determine whether we are at a physical boundary and, if so,
                // whether it is a Dirichlet boundary.
                bool at_physical_bdry = elem->neighbor(side) == NULL;
                bool at_dirichlet_bdry = false;
                const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
                for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
                {
                    const short int bdry_id = *cit;
                    at_physical_bdry  = at_physical_bdry  && !dof_map.is_periodic_boundary(bdry_id);
                    at_dirichlet_bdry = at_dirichlet_bdry || (bdry_id == FEDataManager::DIRICHLET_BDRY_ID);
                }

                // Skip non-physical boundaries.
                if (!at_physical_bdry) continue;

                // Determine whether we need to compute surface forces along
                // this part of the physical boundary; if not, skip the present
                // side.
                const bool compute_transmission_force = d_split_interior_and_bdry_forces && !at_dirichlet_bdry;
                const bool compute_pressure = (d_split_interior_and_bdry_forces && d_lag_pressure_fcns[part] != NULL);
                const bool compute_surface_force = (d_split_interior_and_bdry_forces && d_lag_surface_force_fcns[part] != NULL);
                if (!(compute_transmission_force || compute_pressure || compute_surface_force)) continue;

                // Construct a side element.
                AutoPtr<Elem> side_elem = elem->build_side(side);
                const unsigned int n_node_side = side_elem->n_nodes();
                for (int d = 0; d < NDIM; ++d)
                {
                    dof_map.dof_indices(side_elem.get(), side_dof_indices(d), d);
                }

                // Cache the nodal and physical coordinates of the side element,
                // determine the bounding box of the current configuration of
                // the side element, and set the nodal coordinates to correspond
                // to the physical coordinates.
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(n_node_side <= MAX_NODES);
#endif
                X_min =  0.5*std::numeric_limits<double>::max();
                X_max = -0.5*std::numeric_limits<double>::max();
                for (unsigned int k = 0; k < n_node_side; ++k)
                {
                    s_node_cache[k] = side_elem->point(k);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X_node_cache[k](d) = X_ghost_vec(side_dof_indices(d)[k]);
                        X_min[d] = std::min(X_min[d],X_node_cache[k](d));
                        X_max[d] = std::max(X_max[d],X_node_cache[k](d));
                    }
                    side_elem->point(k) = X_node_cache[k];
                }

                // Loop over coordinate directions and look for intersections
                // with the background fluid grid.
                std::vector<Point> intersection_master_coords;
                std::vector<int>   intersection_axes;
                static const int estimated_max_size = (NDIM == 2 ? 64 : 512);
                intersection_master_coords.reserve(estimated_max_size);
                intersection_axes         .reserve(estimated_max_size);
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    // Setup a unit vector pointing in the appropriate
                    // coordinate direction.
                    VectorValue<double> q;
                    q(axis) = 1.0;

                    // Loop over the relevant range of indices.
                    blitz::TinyVector<int,NDIM> i_begin, i_end, ic;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            i_begin[d] = 0;
                            i_end  [d] = 1;
                        }
                        else
                        {
                            i_begin[d] = std::ceil((X_min[d]-x_lower[d])/dx[d] - 0.5) + patch_lower[d];
                            i_end  [d] = std::ceil((X_max[d]-x_lower[d])/dx[d] - 0.5) + patch_lower[d];
                        }
                    }
#if (NDIM == 3)
                    for (ic[2] = i_begin[2]; ic[2] < i_end[2]; ++ic[2])
                    {
#endif
                        for (ic[1] = i_begin[1]; ic[1] < i_end[1]; ++ic[1])
                        {
                            for (ic[0] = i_begin[0]; ic[0] < i_end[0]; ++ic[0])
                            {
                                Point r;
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    r(d) = (d == axis ? 0.0 : x_lower[d] + dx[d]*(static_cast<double>(ic[d]-patch_lower[d])+0.5));
                                }
#if (NDIM == 2)
                                std::vector<std::pair<double,Point> > intersections = intersect_line_with_edge(dynamic_cast<Edge*>(side_elem.get()), r, q);
#endif
#if (NDIM == 3)
                                std::vector<std::pair<double,Point> > intersections = intersect_line_with_face(dynamic_cast<Face*>(side_elem.get()), r, q);
#endif
                                for (unsigned int k = 0; k < intersections.size(); ++k)
                                {
                                    intersection_master_coords.push_back(intersections[k].second);
                                    intersection_axes.push_back(axis);
                                }
                            }
                        }
#if (NDIM == 3)
                    }
#endif
                }

                // Restore the element coordinates.
                for (unsigned int k = 0; k < n_node_side; ++k)
                {
                    side_elem->point(k) = s_node_cache[k];
                }

                // If there are no intersection points, then continue to the
                // next side.
                if (intersection_master_coords.empty()) continue;

                // Evaluate the jump conditions and apply them to the Eulerian
                // grid.
                fe_face->reinit(elem, side, TOLERANCE, &intersection_master_coords);

                if (!d_use_IB_spreading_operator)
                {
                    get_values_for_interpolation(F_node, F_ghost_vec, dof_indices);
                }
                get_values_for_interpolation(X_node, X_ghost_vec, dof_indices);
                for (unsigned int qp = 0; qp < intersection_master_coords.size(); ++qp)
                {
                    const unsigned int axis = intersection_axes[qp];
                    interpolate(X_qp,qp,X_node,phi_face);
                    Index<NDIM> i = IndexUtilities::getCellIndex(&X_qp(0), x_lower, x_upper, dx, patch_lower, patch_upper);
                    if (X_qp(axis) > x_lower[axis] + static_cast<double>(i(axis)-patch_lower[axis]+0.5)*dx[axis])
                    {
                        ++i(axis);
                    }
#ifdef DEBUG_CHECK_ASSERTIONS
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            TBOX_ASSERT(x_lower[d]+(static_cast<double>(i(d)-patch_lower[d])-0.5)*dx[d] <= X_qp(d) && x_lower[d]+(static_cast<double>(i(d)-patch_lower[d])+0.5)*dx[d] >= X_qp(d));
                        }
                        else
                        {
                            TBOX_ASSERT(std::abs(x_lower[d]+(static_cast<double>(i(d)-patch_lower[d])+0.5)*dx[d] - X_qp(d)) < sqrt(std::numeric_limits<double>::epsilon()));
                        }
                    }
#endif
                    const int d = axis;
                    const SideIndex<NDIM> s_i(i,d,0);
                    if (!side_boxes[d].contains(s_i)) continue;
                    const Point& s_qp = q_point_face[qp];
                    interpolate(X_qp,qp,X_node,phi_face);
                    jacobian(dX_ds,qp,X_node,dphi_face);
                    const double J = std::abs(dX_ds.det());
                    tensor_inverse_transpose(dX_ds_inv_trans,dX_ds,NDIM);
                    F.zero();
                    if (compute_transmission_force)
                    {
                        // Compute the value of the first Piola-Kirchhoff stress
                        // tensor at the quadrature point and compute the
                        // corresponding force.
                        d_PK1_stress_fcns[part](PP,dX_ds,X_qp,s_qp,elem,X_ghost_vec,PK1_stress_fcn_data,time,d_PK1_stress_fcn_ctxs[part]);
                        F -= PP*normal_face[qp];
                    }
                    if (compute_pressure)
                    {
                        // Compute the value of the pressure at the quadrature
                        // point and compute the corresponding force.
                        d_lag_pressure_fcns[part](P,dX_ds,X_qp,s_qp,elem,side,X_ghost_vec,lag_pressure_fcn_data,time,d_lag_pressure_fcn_ctxs[part]);
                        F -= P*J*dX_ds_inv_trans*normal_face[qp];
                    }
                    if (compute_surface_force)
                    {
                        // Compute the value of the surface force at the
                        // quadrature point and compute the corresponding force.
                        d_lag_surface_force_fcns[part](F_s,dX_ds,X_qp,s_qp,elem,side,X_ghost_vec,lag_surface_force_fcn_data,time,d_lag_surface_force_fcn_ctxs[part]);
                        F += F_s;
                    }

                    // Use Nanson's formula (n da = J F^{-T} N dA) to convert
                    // force per unit area in the reference configuration into
                    // force per unit area in the current configuration.  This
                    // value determines the discontinuity in the pressure at the
                    // fluid-structure interface.
                    n = (dX_ds_inv_trans*normal_face[qp]).unit();
                    const double dA_da = 1.0/(J*(dX_ds_inv_trans*normal_face[qp])*n);
                    F *= dA_da;

                    // Determine the value of the interior force density at the
                    // boundary, and convert it to force per unit volume in the
                    // current configuration.  This value determines the
                    // discontinuity in the normal derivative of the pressure at
                    // the fluid-structure interface.
                    //
                    // NOTE: This additional correction appears to be
                    // ineffective when we use "diffuse" force spreading; hence,
                    // we compute it only when we do NOT use the IB/FE version
                    // of the IB force spreading operator.
                    if (d_use_IB_spreading_operator)
                    {
                        F_qp.zero();
                    }
                    else
                    {
                        interpolate(F_qp,qp,F_node,phi_face);
                        F_qp /= J;
                    }

                    // Impose the jump conditions.
                    const double X = X_qp(d);
                    const double x_cell_bdry = x_lower[d]+static_cast<double>(i(d)-patch_lower[d])*dx[d];
                    const double h = x_cell_bdry + (X > x_cell_bdry ? +0.5 : -0.5)*dx[d] - X;
                    const double C_p = F*n - h*F_qp(d);
                    (*f_data)(s_i) += (n(d) > 0.0 ? +1.0 : -1.0)*(C_p/dx[d]);
                }
            }
        }
    }
    return;
}// imposeJumpConditions

void
IBFEHierarchyIntegrator::initializeCoordinates(
    const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    const bool identity_mapping = d_coordinate_mapping_fcns[part] == NULL;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num) > 0)
        {
            libmesh_assert(n->n_vars(X_sys_num) == NDIM);
            const Point& s = *n;
            Point X = s;
            if (!identity_mapping)
            {
                d_coordinate_mapping_fcns[part](X, s, d_coordinate_mapping_fcn_ctxs[part]);
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int dof_index = n->dof_number(X_sys_num,d,0);
                X_coords.set(dof_index,X(d));
            }
        }
    }
    X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
    X_coords.close();
    return;
}// initializeCoordinates

void
IBFEHierarchyIntegrator::updateCoordinateMapping(
    const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    System& X_mapping_system = equation_systems->get_system(COORD_MAPPING_SYSTEM_NAME);
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
    d_regrid_cfl_interval = db->getDoubleWithDefault("regrid_cfl_interval", d_regrid_cfl_interval);

    d_dt_max = db->getDoubleWithDefault("dt_max", d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault("dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault("dt_max_time_min", d_dt_max_time_min);

    d_use_IB_spreading_operator = db->getBoolWithDefault("use_IB_spreading_operator", d_use_IB_spreading_operator);
    d_use_IB_interpolation_operator = db->getBoolWithDefault("use_IB_interpolation_operator", d_use_IB_interpolation_operator);
    d_split_interior_and_bdry_forces = db->getBoolWithDefault("split_interior_and_bdry_forces", d_split_interior_and_bdry_forces);
    d_use_jump_conditions = db->getBoolWithDefault("use_jump_conditions", d_use_jump_conditions);
    d_use_consistent_mass_matrix = db->getBoolWithDefault("use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    d_quad_type = Utility::string_to_enum<QuadratureType>(db->getStringWithDefault("quad_type", Utility::enum_to_string<QuadratureType>(d_quad_type)));
    d_quad_order = Utility::string_to_enum<Order>(db->getStringWithDefault("quad_order", Utility::enum_to_string<Order>(d_quad_order)));

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

    d_use_IB_spreading_operator = db->getBool("d_use_IB_spreading_operator");
    d_use_IB_interpolation_operator = db->getBool("d_use_IB_interpolation_operator");
    d_split_interior_and_bdry_forces = db->getBool("d_split_interior_and_bdry_forces");
    d_use_jump_conditions = db->getBool("d_use_jump_conditions");
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    d_grow_dt = db->getDouble("d_grow_dt");
    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");
    d_regrid_interval = db->getInteger("d_regrid_interval");
    d_regrid_cfl_interval = db->getDouble("d_regrid_cfl_interval");
    d_regrid_cfl_estimate = db->getDouble("d_regrid_cfl_estimate");
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
template class Pointer<IBAMR::IBFEHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
