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
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IndexUtilities.h>
#include <ibtk/LMarkerSetData.h>
#include <ibtk/LMarkerUtilities.h>
#if (NDIM == 3)
#include <ibtk/LM3DDataWriter.h>
#endif
#include <ibtk/LSiloDataWriter.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBHierarchyIntegrator restart file data.
static const int IB_HIERARCHY_INTEGRATOR_VERSION = 2;

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
    Pointer<INSHierarchyIntegrator> ins_hier_integrator,
    bool register_for_restart)
    : HierarchyIntegrator(object_name, input_db, register_for_restart)
{
    // Register the fluid solver as a child integrator of this integrator
    // object.
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!ins_hier_integrator.isNull());
#endif
    d_ins_hier_integrator = ins_hier_integrator;
    registerChildHierarchyIntegrator(d_ins_hier_integrator);

    // Ensure all pointers to helper objects are NULL.
    d_silo_writer = NULL;
#if (NDIM == 3)
    d_m3D_writer = NULL;
#endif
    d_load_balancer = NULL;
    d_l_initializer = NULL;
    d_body_force_fcn = NULL;
    d_ib_force_fcn = NULL;
    d_ib_force_fcn_needs_init = true;
    d_ib_source_fcn = NULL;
    d_ib_source_fcn_needs_init = true;
    d_normalize_source_strength = false;
    d_post_processor = NULL;
    d_post_processor_needs_init = true;

    // Set some default values.
    d_integrator_is_initialized = false;
    d_num_cycles = 3;
    d_interp_delta_fcn = "IB_4";
    d_spread_delta_fcn = "IB_4";
    d_ghosts = -1;
    d_using_pIB_method = false;
    d_gravitational_acceleration = 0.0;
    d_regrid_interval = 1;
    d_regrid_cfl_interval = 0.0;
    d_regrid_cfl_estimate = 0.0;

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Check the choices for the delta function.
    if (d_interp_delta_fcn != d_spread_delta_fcn)
    {
        pout << "WARNING: different delta functions are being used for velocity interpolation and force spreading.\n"
             << "         recommended usage is to employ the same delta functions for both interpolation and spreading.\n";
    }

    // Get the Lagrangian Data Manager.
    d_l_data_manager = LDataManager::getManager(d_object_name+"::LDataManager", d_interp_delta_fcn, d_spread_delta_fcn, d_ghosts, d_registered_for_restart);
    d_ghosts = d_l_data_manager->getGhostCellWidth();

    // Read in the marker initial positions.
    if (!from_restart) LMarkerUtilities::readMarkerPositions(d_mark_init_posns, d_mark_input_file_name, d_hierarchy);

    // Create the instrument panel object.
    d_instrument_panel = new IBInstrumentPanel(d_object_name+"::IBInstrumentPanel", (input_db->isDatabase("IBInstrumentPanel") ? input_db->getDatabase("IBInstrumentPanel") : Pointer<Database>(NULL)));
    return;
}// IBHierarchyIntegrator

IBHierarchyIntegrator::~IBHierarchyIntegrator()
{
    // intentionally blank.
    return;
}// ~IBHierarchyIntegrator

void
IBHierarchyIntegrator::registerBodyForceFunction(
    Pointer<CartGridFunction> F_fcn)
{
    d_body_force_fcn = F_fcn;
    return;
}// registerBodyForceFunction

void
IBHierarchyIntegrator::registerIBLagrangianForceFunction(
    Pointer<IBLagrangianForceStrategy> ib_force_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!ib_force_fcn.isNull());
#endif
    d_ib_force_fcn = ib_force_fcn;
    return;
}// registerIBLagrangianForceFunction

void
IBHierarchyIntegrator::registerIBLagrangianSourceFunction(
    Pointer<IBLagrangianSourceStrategy> ib_source_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!ib_source_fcn.isNull());
#endif
    d_ib_source_fcn = ib_source_fcn;
    return;
}// registerIBLagrangianSourceFunction

void
IBHierarchyIntegrator::registerLInitStrategy(
    Pointer<LInitStrategy> l_initializer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!l_initializer.isNull());
#endif
    d_l_initializer = l_initializer;
    d_l_data_manager->registerLInitStrategy(d_l_initializer);
    return;
}// registerLInitStrategy

void
IBHierarchyIntegrator::freeLInitStrategy()
{
    d_l_initializer.setNull();
    d_l_data_manager->freeLInitStrategy();
    return;
}// freeLInitStrategy

LDataManager*
IBHierarchyIntegrator::getLDataManager() const
{
    return d_l_data_manager;
}// getLDataManager

Pointer<LMarkerSetVariable>
IBHierarchyIntegrator::getLMarkerSetVariable() const
{
    return d_mark_var;
}// getLMarkerSetVariable

Pointer<IBInstrumentPanel>
IBHierarchyIntegrator::getIBInstrumentPanel() const
{
    return d_instrument_panel;
}// getIBInstrumentPanel

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

void
IBHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_velocity_data_ops = hier_ops_manager->getOperationsDouble(d_ins_hier_integrator->getVelocityVariable(), hierarchy, true);
    d_hier_pressure_data_ops = hier_ops_manager->getOperationsDouble(d_ins_hier_integrator->getPressureVariable(), hierarchy, true);

    // Initialize all variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const IntVector<NDIM>    ghosts = d_ghosts;
    const IntVector<NDIM> no_ghosts = 0;

    Pointer<CellVariable<NDIM,double> > U_cc_var = d_ins_hier_integrator->getVelocityVariable();
    Pointer<SideVariable<NDIM,double> > U_sc_var = d_ins_hier_integrator->getVelocityVariable();
    if (!U_cc_var.isNull())
    {
        d_V_var = new CellVariable<NDIM,double>(d_object_name+"::V", NDIM);
        d_F_var = new CellVariable<NDIM,double>(d_object_name+"::F", NDIM);
    }
    else if (!U_sc_var.isNull())
    {
        d_V_var = new SideVariable<NDIM,double>(d_object_name+"::V");
        d_F_var = new SideVariable<NDIM,double>(d_object_name+"::F");
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchyIntegrator():\n"
                   << "  unsupported velocity data centering" << std::endl);
    }
    d_V_idx = var_db->registerVariableAndContext(d_V_var, getScratchContext(),    ghosts);
    d_F_idx = var_db->registerVariableAndContext(d_F_var, getScratchContext(), no_ghosts);

    Pointer<CellVariable<NDIM,double> > P_cc_var = d_ins_hier_integrator->getPressureVariable();
    if (!P_cc_var.isNull())
    {
        if (!d_ib_source_fcn.isNull()) d_Q_var = new CellVariable<NDIM,double>(d_object_name+"::Q");
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchyIntegrator():\n"
                   << "  unsupported pressure data centering" << std::endl);
    }
    if (!d_ib_source_fcn.isNull()) d_Q_idx = var_db->registerVariableAndContext(d_Q_var, getScratchContext(), no_ghosts);

    d_mark_var = new LMarkerSetVariable(d_object_name+"::mark");
    d_mark_current_idx = var_db->registerVariableAndContext(d_mark_var, getCurrentContext(), ghosts);
    d_mark_scratch_idx = var_db->registerVariableAndContext(d_mark_var, getScratchContext(), ghosts);
    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_mark_current_idx);
    }

    // Initialize the objects used to manage Lagrangian-Eulerian interaction.
    d_eulerian_force_fcn = new IBEulerianForceFunction(d_object_name+"::IBEulerianForceFunction", -1, -1, d_F_idx);
    d_ins_hier_integrator->registerBodyForceFunction(d_eulerian_force_fcn);
    if (!d_ib_source_fcn.isNull())
    {
        d_eulerian_source_fcn = new IBEulerianSourceFunction(d_object_name+"::IBEulerianSourceFunction", d_Q_idx, d_Q_idx, d_Q_idx);
        d_ins_hier_integrator->registerFluidSourceFunction(d_eulerian_source_fcn);
    }

    // Initialize the fluid solver.  It is necessary to do this after
    // registering a body force function or fluid source distribution function.
    d_ins_hier_integrator->initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Create several communications algorithms, used in filling ghost cell data
    // and synchronizing data on the patch hierarchy.
    Pointer<Geometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineOperator<NDIM> > refine_operator;
    Pointer<CoarsenOperator<NDIM> > coarsen_operator;

    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getCurrentContext());
    const int U_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getNewContext());
    const int U_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getScratchContext());
    const int P_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVariable(), d_ins_hier_integrator->getNewContext());
    const int P_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVariable(), d_ins_hier_integrator->getScratchContext());

    d_refine_algs["U->V::C->S::GHOST_FILL"] = new RefineAlgorithm<NDIM>();
    refine_operator = NULL;
    d_refine_algs["U->V::C->S::GHOST_FILL"]->registerRefine(d_V_idx, U_current_idx, d_V_idx, refine_operator);

    d_refine_algs["V->V::S->S::GHOST_FILL"] = new RefineAlgorithm<NDIM>();
    refine_operator = NULL;
    d_refine_algs["V->V::S->S::GHOST_FILL"]->registerRefine(d_V_idx, d_V_idx, d_V_idx, refine_operator);

    d_coarsen_algs["U->U::C->C::CONSERVATIVE_COARSEN"] = new CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(d_V_var, "CONSERVATIVE_COARSEN");
    d_coarsen_algs["U->U::C->C::CONSERVATIVE_COARSEN"]->registerCoarsen(U_current_idx, U_current_idx, coarsen_operator);

    d_coarsen_algs["V->V::S->S::CONSERVATIVE_COARSEN"] = new CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(d_V_var, "CONSERVATIVE_COARSEN");
    d_coarsen_algs["V->V::S->S::CONSERVATIVE_COARSEN"]->registerCoarsen(d_V_idx, d_V_idx, coarsen_operator);

    d_refine_algs["INSTRUMENTATION_DATA_FILL"] = new RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getVelocityVariable(), "CONSERVATIVE_LINEAR_REFINE");
    d_refine_algs["INSTRUMENTATION_DATA_FILL"]->registerRefine(U_scratch_idx, U_new_idx, U_scratch_idx, refine_operator);
    refine_operator = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getPressureVariable(), "LINEAR_REFINE");
    d_refine_algs["INSTRUMENTATION_DATA_FILL"]->registerRefine(P_scratch_idx, P_new_idx, P_scratch_idx, refine_operator);
    ComponentSelector instrumentation_data_fill_bc_idxs;
    instrumentation_data_fill_bc_idxs.setFlag(U_scratch_idx);
    instrumentation_data_fill_bc_idxs.setFlag(P_scratch_idx);
    d_refine_strategies["INSTRUMENTATION_DATA_FILL"] = new CartExtrapPhysBdryOp(instrumentation_data_fill_bc_idxs, "LINEAR");

    // Set the current integration time.
    if (!RestartManager::getManager()->isFromRestart())
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
}// initializeHierarchyIntegrator

void
IBHierarchyIntegrator::initializePatchHierarchy(
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_hierarchy_is_initialized) return;

    // Initialize Eulerian data.
    HierarchyIntegrator::initializePatchHierarchy(hierarchy, gridding_alg);

    // Begin Lagrangian data movement.
    d_l_data_manager->beginDataRedistribution();

    // Finish Lagrangian data movement.
    d_l_data_manager->endDataRedistribution();

    // Update the workload estimates.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    d_l_data_manager->updateWorkloadData(coarsest_ln, finest_ln);

    // Initialize various Lagrangian data objects.
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time,d_start_time);
    if (initial_time)
    {
        const double init_data_time = d_integrator_time;

        // Initialize the interpolated velocity field.
        std::vector<Pointer<LData> > X_data(finest_ln+1);
        std::vector<Pointer<LData> > U_data(finest_ln+1);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(d_V_idx, init_data_time);
            if (d_l_data_manager->levelContainsLagrangianData(ln))
            {
                X_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
                U_data[ln] = d_l_data_manager->getLData(LDataManager:: VEL_DATA_NAME,ln);
            }
        }
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getCurrentContext());
        d_hier_velocity_data_ops->copyData(d_V_idx, U_current_idx);
        d_l_data_manager->interp(d_V_idx, U_data, X_data, d_coarsen_scheds["V->V::S->S::CONSERVATIVE_COARSEN"], d_refine_scheds["V->V::S->S::GHOST_FILL"], init_data_time);
        resetAnchorPointValues(U_data, coarsest_ln, finest_ln);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(d_V_idx);
        }

        for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
        {
            const bool can_be_refined = level_number < finest_ln || d_gridding_alg->levelCanBeRefined(level_number);

            // Initialize marker data.
            LMarkerUtilities::initializeMarkersOnLevel(d_mark_current_idx, d_mark_init_posns, d_hierarchy, level_number, initial_time, Pointer<PatchLevel<NDIM> >(NULL));

            // Initialize source/sink data.
            if (!d_ib_source_fcn.isNull())
            {
                d_ib_source_fcn->initializeLevelData(d_hierarchy, level_number, init_data_time, initial_time, d_l_data_manager);
                d_n_src[level_number] = d_ib_source_fcn->getNumSources(d_hierarchy, level_number, d_integrator_time, d_l_data_manager);
                d_X_src[level_number].resize(d_n_src[level_number], blitz::TinyVector<double,NDIM>(std::numeric_limits<double>::quiet_NaN()));
                d_r_src[level_number].resize(d_n_src[level_number], std::numeric_limits<double>::quiet_NaN());
                d_P_src[level_number].resize(d_n_src[level_number], std::numeric_limits<double>::quiet_NaN());
                d_Q_src[level_number].resize(d_n_src[level_number], std::numeric_limits<double>::quiet_NaN());
                d_ib_source_fcn->getSourceLocations(d_X_src[level_number], d_r_src[level_number], X_data[level_number], d_hierarchy, level_number, init_data_time, d_l_data_manager);
            }

            // Setup data needed for pIB method.
            if (d_l_data_manager->levelContainsLagrangianData(level_number) && d_using_pIB_method)
            {
                static const bool manage_data = true;
                Pointer<LData> M_data = d_l_data_manager->createLData("M",level_number,1,manage_data);
                Pointer<LData> K_data = d_l_data_manager->createLData("K",level_number,1,manage_data);
                Pointer<LData> Y_data = d_l_data_manager->createLData("Y",level_number,NDIM,manage_data);
                Pointer<LData> dY_dt_data = d_l_data_manager->createLData("dY_dt",level_number,NDIM,manage_data);
                static const int global_index_offset = 0;
                static const int local_index_offset = 0;
                d_l_initializer->initializeMassDataOnPatchLevel(global_index_offset, local_index_offset, M_data, K_data, d_hierarchy, level_number, init_data_time, can_be_refined, initial_time, d_l_data_manager);
                if (!d_silo_writer.isNull())
                {
                    d_silo_writer->registerVariableData("M", M_data, level_number);
                    d_silo_writer->registerVariableData("Y", Y_data, level_number);
                }

                // Set initial conditions.
                Vec X_vec = X_data[level_number]->getVec();
                Vec U_vec = U_data[level_number]->getVec();
                Vec Y_vec = Y_data->getVec();
                Vec dY_dt_vec = dY_dt_data->getVec();
                PetscErrorCode ierr;
                ierr = VecCopy(X_vec, Y_vec);  IBTK_CHKERRQ(ierr);
                ierr = VecCopy(U_vec, dY_dt_vec);  IBTK_CHKERRQ(ierr);
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
    d_ib_force_fcn_needs_init = true;
    d_ib_source_fcn_needs_init = true;
    d_post_processor_needs_init = true;

    // Indicate that the hierarchy is initialized.
    d_hierarchy_is_initialized = true;
    return;
}// initializePatchHierarchy

void
IBHierarchyIntegrator::preprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const int num_cycles)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const bool initial_time = MathUtilities<double>::equalEps(current_time,d_start_time);

    // Setup the Eulerian body force and fluid source functions.
    d_eulerian_force_fcn->registerBodyForceFunction(d_body_force_fcn);
    d_eulerian_force_fcn->setTimeInterval(current_time, new_time);
    if (!d_ib_force_fcn.isNull())
    {
        d_ib_force_fcn->setTimeInterval(current_time, new_time);
    }
    if (!d_ib_source_fcn.isNull())
    {
        d_eulerian_source_fcn->setTimeInterval(current_time, new_time);
        d_ib_source_fcn->setTimeInterval(current_time, new_time);
    }

    // Reset Lagrangian helper objects.
    if (d_ib_force_fcn_needs_init && !d_ib_force_fcn.isNull())
    {
        resetLagrangianForceFunction(current_time, initial_time);
        d_ib_force_fcn_needs_init = false;
    }
    if (d_ib_source_fcn_needs_init && !d_ib_source_fcn.isNull())
    {
        resetLagrangianSourceFunction(current_time, initial_time);
        d_ib_source_fcn_needs_init = false;
    }
    if (d_post_processor_needs_init && !d_post_processor.isNull())
    {
        resetPostProcessor(current_time, initial_time);
        d_post_processor_needs_init = false;
    }

    // Allocate Eulerian scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_V_idx, current_time);
        level->allocatePatchData(d_F_idx, current_time);
        if (!d_ib_source_fcn.isNull())
        {
            level->allocatePatchData(d_Q_idx, current_time);
        }
    }

    // Allocate Lagangian scratch data.
    d_X_current_data.resize(finest_ln+1);
    d_U_current_data.resize(finest_ln+1);
    d_X_new_data    .resize(finest_ln+1);
    d_X_half_data   .resize(finest_ln+1);
    d_U_half_data   .resize(finest_ln+1);
    d_F_half_data   .resize(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_X_current_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
        d_U_current_data[ln] = d_l_data_manager->getLData(LDataManager:: VEL_DATA_NAME,ln);
        d_X_new_data    [ln] = d_l_data_manager->createLData("X_new" ,ln,NDIM);
        d_X_half_data   [ln] = d_l_data_manager->createLData("X_half",ln,NDIM);
        d_U_half_data   [ln] = d_l_data_manager->createLData("U_half",ln,NDIM);
        d_F_half_data   [ln] = d_l_data_manager->createLData("F_half",ln,NDIM);
    }

    if (d_using_pIB_method)
    {
        d_K_data            .resize(finest_ln+1);
        d_M_data            .resize(finest_ln+1);
        d_Y_current_data    .resize(finest_ln+1);
        d_dY_dt_current_data.resize(finest_ln+1);
        d_Y_new_data        .resize(finest_ln+1);
        d_dY_dt_new_data    .resize(finest_ln+1);
        d_F_K_half_data     .resize(finest_ln+1);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            d_K_data            [ln] = d_l_data_manager->getLData("K"    ,ln);
            d_M_data            [ln] = d_l_data_manager->getLData("M"    ,ln);
            d_Y_current_data    [ln] = d_l_data_manager->getLData("Y"    ,ln);
            d_dY_dt_current_data[ln] = d_l_data_manager->getLData("dY_dt",ln);
            d_Y_new_data        [ln] = d_l_data_manager->createLData("Y_new"    ,ln,NDIM);
            d_dY_dt_new_data    [ln] = d_l_data_manager->createLData("dY_dt_new",ln,NDIM);
            d_F_K_half_data     [ln] = d_l_data_manager->createLData("F_K_half" ,ln,NDIM);
        }
    }

    // Synchronize the Cartesian grid velocity u(n) on the patch hierarchy.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_coarsen_scheds["U->U::C->C::CONSERVATIVE_COARSEN"][ln]->coarsenData();
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_refine_scheds["U->V::C->S::GHOST_FILL"][ln]->fillData(current_time);
    }

    // Initialize X(n+1) to equal X(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_current_vec = d_X_current_data[ln]->getVec();
            Vec X_new_vec = d_X_new_data[ln]->getVec();
            ierr = VecCopy(X_current_vec, X_new_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Initialize U(n+1/2) to equal U(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec U_current_vec = d_U_current_data[ln]->getVec();
            Vec U_half_vec = d_U_half_data[ln]->getVec();
            ierr = VecCopy(U_current_vec, U_half_vec);  IBTK_CHKERRQ(ierr);
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

    // Initialize Y(n+1) to equal Y(n), and initialize dY/dt(n+1) to equal
    // dY/dt(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln) && d_using_pIB_method)
        {
            Vec Y_current_vec = d_Y_current_data[ln]->getVec();
            Vec Y_new_vec = d_Y_new_data[ln]->getVec();
            ierr = VecCopy(Y_current_vec, Y_new_vec);  IBTK_CHKERRQ(ierr);

            Vec dY_dt_current_vec = d_dY_dt_current_data[ln]->getVec();
            Vec dY_dt_new_vec = d_dY_dt_new_data[ln]->getVec();
            ierr = VecCopy(dY_dt_current_vec, dY_dt_new_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Initialize the fluid solver.
    d_ins_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, num_cycles);
    return;
}// preprocessIntegrateHierarchy

void
IBHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time,
    const int cycle_num)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getCurrentContext());
    const int U_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getNewContext());

    // Set X(n+1/2) = 0.5*(X(n) + X(n+1)).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        Vec X_current_vec = d_X_current_data[ln]->getVec();
        Vec X_new_vec = d_X_new_data[ln]->getVec();
        Vec X_half_vec = d_X_half_data[ln]->getVec();
        ierr = VecCopy(X_current_vec, X_half_vec);  IBTK_CHKERRQ(ierr);
        ierr = VecAXPBY(X_half_vec, 0.5, 0.5, X_new_vec);  IBTK_CHKERRQ(ierr);
    }

    // Compute F(n+1/2) = F(X(n+1/2),U(n+1/2),t_{n+1/2}).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        Vec F_half_vec = d_F_half_data[ln]->getVec();
        ierr = VecSet(F_half_vec, 0.0);  IBTK_CHKERRQ(ierr);
        if (!d_ib_force_fcn.isNull())
        {
            d_ib_force_fcn->computeLagrangianForce(d_F_half_data[ln], d_X_half_data[ln], d_U_half_data[ln], d_hierarchy, ln, current_time+0.5*dt, d_l_data_manager);
        }
    }

    // Compute pIB penalty forces, F_K = K*(Y-X), and update the pIB state
    // variables, Y and dY/dt.
    if (d_using_pIB_method)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

            Vec K_vec = d_K_data[ln]->getVec();
            Vec M_vec = d_M_data[ln]->getVec();
            Vec X_half_vec = d_X_half_data[ln]->getVec();
            Vec Y_current_vec = d_Y_current_data[ln]->getVec();
            Vec Y_new_vec = d_Y_new_data[ln]->getVec();
            Vec dY_dt_current_vec = d_dY_dt_current_data[ln]->getVec();
            Vec dY_dt_new_vec = d_dY_dt_new_data[ln]->getVec();
            Vec F_K_half_vec = d_F_K_half_data[ln]->getVec();

            int n_local = 0;
            ierr = VecGetLocalSize(M_vec, &n_local);  IBTK_CHKERRQ(ierr);

            double* K_arr, * M_arr, * X_half_arr, * Y_arr, * Y_new_arr, * dY_dt_arr, * dY_dt_new_arr, * F_K_half_arr;
            ierr = VecGetArray(K_vec, &K_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecGetArray(M_vec, &M_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecGetArray(X_half_vec, &X_half_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecGetArray(Y_current_vec, &Y_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecGetArray(Y_new_vec, &Y_new_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecGetArray(dY_dt_current_vec, &dY_dt_arr);  IBTK_CHKERRQ(ierr);
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
                plog << d_object_name << "::integrateHierarchy():\n";
                plog << "  maximum massive boundary point displacement [present configuration] = " << max_config_displacement << "\n";
                plog << "  maximum massive boundary point displacement [entire simulation]     = " << max_displacement << "\n";
            }

            ierr = VecRestoreArray(K_vec, &K_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecRestoreArray(M_vec, &M_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecRestoreArray(X_half_vec, &X_half_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecRestoreArray(Y_current_vec, &Y_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecRestoreArray(Y_new_vec, &Y_new_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecRestoreArray(dY_dt_current_vec, &dY_dt_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecRestoreArray(dY_dt_new_vec, &dY_dt_new_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecRestoreArray(F_K_half_vec, &F_K_half_arr);  IBTK_CHKERRQ(ierr);

            Vec F_half_vec = d_F_half_data[ln]->getVec();
            ierr = VecAXPY(F_half_vec,1.0,F_K_half_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Spread F(n+1/2) to f(n+1/2).
    d_hier_velocity_data_ops->setToScalar(d_F_idx, 0.0);
    resetAnchorPointValues(d_F_half_data, coarsest_ln, finest_ln);
    d_l_data_manager->spread(d_F_idx, d_F_half_data, d_X_half_data, d_refine_scheds["NONE"], true, true);

    // Compute the source/sink strengths corresponding to any distributed
    // internal fluid sources or sinks.
    if (!d_ib_source_fcn.isNull())
    {
        computeSourceStrengths(coarsest_ln, finest_ln, current_time+0.5*dt, d_X_half_data);
    }

    // Solve the incompressible Navier-Stokes equations.
    d_ins_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);

    // Set u(n+1/2) = 0.5*(u(n) + u(n+1)) and interpolate u(n+1/2) to U(n+1/2).
    d_hier_velocity_data_ops->linearSum(d_V_idx, 0.5, U_current_idx, 0.5, U_new_idx);
    d_l_data_manager->interp(d_V_idx, d_U_half_data, d_X_half_data, d_coarsen_scheds["V->V::S->S::CONSERVATIVE_COARSEN"], d_refine_scheds["V->V::S->S::GHOST_FILL"], current_time);
    resetAnchorPointValues(d_U_half_data, coarsest_ln, finest_ln);

    // Set X(n+1) = X(n) + dt*U(n+1/2).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        Vec X_current_vec = d_X_current_data[ln]->getVec();
        Vec X_new_vec = d_X_new_data[ln]->getVec();
        Vec U_half_vec = d_U_half_data[ln]->getVec();
        ierr = VecWAXPY(X_new_vec, dt, U_half_vec, X_current_vec);  IBTK_CHKERRQ(ierr);
    }
    LMarkerUtilities::advectMarkers(d_mark_current_idx, d_mark_scratch_idx, d_V_idx, dt, d_interp_delta_fcn, d_hierarchy);

    // Compute the pressure at the updated locations of any distributed internal
    // fluid sources or sinks.
    if (!d_ib_source_fcn.isNull())
    {
        computeSourcePressures(coarsest_ln, finest_ln, current_time+0.5*dt, d_X_half_data);
    }
    return;
}// integrateHierarchy

void
IBHierarchyIntegrator::postprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const bool skip_synchronize_new_state_data,
    const int num_cycles)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getCurrentContext());
    const int U_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getNewContext());

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
            if (d_do_log) plog << "flow volume through " << instrument_name[m] << ":\t " << d_total_flow_volume[m] << "\n";
        }
    }

    // Reset X_current to equal X_new.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_current_vec = d_X_current_data[ln]->getVec();
            Vec X_new_vec = d_X_new_data[ln]->getVec();
            ierr = VecCopy(X_new_vec, X_current_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Interpolate u(n+1) to U(n+1).
    d_hier_velocity_data_ops->copyData(d_V_idx, U_new_idx);
    d_l_data_manager->interp(d_V_idx, d_U_current_data, d_X_current_data, d_coarsen_scheds["V->V::S->S::CONSERVATIVE_COARSEN"], d_refine_scheds["V->V::S->S::GHOST_FILL"], d_integrator_time);
    resetAnchorPointValues(d_U_current_data, coarsest_ln, finest_ln);

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

    // Reset Y to equal Y_new, and reset dY_dt to equal dY_dt_new.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln) && d_using_pIB_method)
        {
            Vec Y_current_vec = d_Y_current_data[ln]->getVec();
            Vec Y_new_vec = d_Y_new_data[ln]->getVec();
            ierr = VecCopy(Y_new_vec, Y_current_vec);  IBTK_CHKERRQ(ierr);

            Vec dY_dt_current_vec = d_dY_dt_current_data[ln]->getVec();
            Vec dY_dt_new_vec = d_dY_dt_new_data[ln]->getVec();
            ierr = VecCopy(dY_dt_new_vec, dY_dt_current_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Deallocate Eulerian scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_V_idx);
        level->deallocatePatchData(d_F_idx);
        if (!d_ib_source_fcn.isNull())
        {
            level->deallocatePatchData(d_Q_idx);
        }
    }

    // Deallocate Lagrangian scratch data.
    d_X_current_data.clear();
    d_U_current_data.clear();
    d_X_new_data    .clear();
    d_X_half_data   .clear();
    d_U_half_data   .clear();
    d_F_half_data   .clear();
    if (d_using_pIB_method)
    {
        d_K_data            .clear();
        d_M_data            .clear();
        d_Y_current_data    .clear();
        d_dY_dt_current_data.clear();
        d_Y_new_data        .clear();
        d_dY_dt_new_data    .clear();
        d_F_K_half_data     .clear();
    }

    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_do_log) plog << d_object_name << "::postprocessIntegrateHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(CURRENT_DATA);
    }

    // Determine the current flow CFL number.
    double cfl_max = 0.0;
    PatchCellDataOpsReal<NDIM,double> patch_cc_ops;
    PatchSideDataOpsReal<NDIM,double> patch_sc_ops;
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
            Pointer<CellData<NDIM,double> > U_cc_current_data = patch->getPatchData(U_current_idx);
            Pointer<SideData<NDIM,double> > U_sc_current_data = patch->getPatchData(U_current_idx);
            double U_max = 0.0;
            if (!U_cc_current_data.isNull()) U_max = patch_cc_ops.maxNorm(U_cc_current_data, patch_box);
            if (!U_sc_current_data.isNull()) U_max = patch_sc_ops.maxNorm(U_sc_current_data, patch_box);
            cfl_max = std::max(cfl_max, U_max*dt/dx_min);
        }
    }
    cfl_max = SAMRAI_MPI::maxReduction(cfl_max);
    d_regrid_cfl_estimate += cfl_max;
    if (d_do_log) plog << d_object_name << "::postprocessIntegrateHierarchy(): current CFL number = " << cfl_max << "\n";
    if (d_do_log) plog << d_object_name << "::postprocessIntegrateHierarchy(): estimated upper bound on IB point displacement since last regrid = " << d_regrid_cfl_estimate << "\n";

    // Deallocate the fluid solver.
d_ins_hier_integrator->postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
}// postprocessIntegrateHierarchy

void
IBHierarchyIntegrator::regridHierarchy()
{
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

    // Use the INSHierarchyIntegrator to handle Eulerian data management.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): calling INSHierarchyIntegrator::regridHierarchy().\n";
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
    d_ib_force_fcn_needs_init  = true;
    d_ib_source_fcn_needs_init = true;
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
    return;
}// regridHierarchy

void
IBHierarchyIntegrator::postProcessData()
{
    if (d_post_processor.isNull()) return;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getCurrentContext());
    const int P_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVariable(), d_ins_hier_integrator->getCurrentContext());
    const int F_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getForceVariable(), d_ins_hier_integrator->getCurrentContext());

    PetscErrorCode ierr;

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
    d_hier_velocity_data_ops->copyData(d_V_idx, U_current_idx);
    d_l_data_manager->interp(d_V_idx, U_data, X_data, d_coarsen_scheds["V->V::S->S::CONSERVATIVE_COARSEN"], d_refine_scheds["V->V::S->S::GHOST_FILL"], d_integrator_time);
    resetAnchorPointValues(U_data, coarsest_ln, finest_ln);

    // Compute F(n) = F(X(n),U(n),n), the Lagrangian force corresponding to
    // configuration X(n) and velocity U(n) at time t_{n}.
    if (d_ib_force_fcn_needs_init)
    {
        const bool initial_time = MathUtilities<double>::equalEps(current_time,d_start_time);
        resetLagrangianForceFunction(current_time, initial_time);
        d_ib_force_fcn_needs_init = false;
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec F_vec = F_data[ln]->getVec();
            ierr = VecSet(F_vec, 0.0);  IBTK_CHKERRQ(ierr);
            if (!d_ib_force_fcn.isNull())
            {
                d_ib_force_fcn->computeLagrangianForce(F_data[ln], X_data[ln], U_data[ln], d_hierarchy, ln, current_time, d_l_data_manager);
            }
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

/////////////////////////////// PROTECTED ////////////////////////////////////

bool
IBHierarchyIntegrator::atRegridPointSpecialized() const
{
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time,d_start_time);
    if (initial_time) return true;
    return d_regrid_cfl_interval > 0.0 ? (d_regrid_cfl_estimate >= d_regrid_cfl_interval) : (d_regrid_interval == 0 ? false : (d_integrator_step % d_regrid_interval == 0));
}// atRegridPointSpecialized

void
IBHierarchyIntegrator::initializeLevelDataSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > base_old_level,
    const bool allocate_data)
{
    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const Pointer<PatchLevel<NDIM> > old_level = base_old_level;
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

    // Use the LDataManager to handle unstructured data management.
    d_l_data_manager->setPatchHierarchy(hierarchy);
    d_l_data_manager->resetLevels(0,hierarchy->getFinestLevelNumber());
    d_l_data_manager->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    return;
}// initializeLevelDataSpecialized

void
IBHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0) && (coarsest_level <= finest_level) && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset the LDataManager.
    d_l_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_velocity_data_ops->setPatchHierarchy(hierarchy);
    d_hier_pressure_data_ops->setPatchHierarchy(hierarchy);
    d_hier_velocity_data_ops->resetLevels(0, finest_hier_level);
    d_hier_pressure_data_ops->resetLevels(0, finest_hier_level);

    // If we have added or removed a level, resize the anchor point vectors.
    d_anchor_point_local_idxs.clear();
    d_anchor_point_local_idxs.resize(finest_hier_level+1);

    // If we have added or removed a level, resize the source/sink data vectors.
    d_X_src.resize(finest_hier_level+1);
    d_r_src.resize(finest_hier_level+1);
    d_P_src.resize(finest_hier_level+1);
    d_Q_src.resize(finest_hier_level+1);
    d_n_src.resize(finest_hier_level+1,0);
    return;
}// resetHierarchyConfigurationSpecialized

void
IBHierarchyIntegrator::applyGradientDetectorSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Tag cells that contain Lagrangian nodes.
    d_l_data_manager->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);

    // Tag cells where the Cartesian source/sink strength is nonzero.
    if (!d_ib_source_fcn.isNull() && !initial_time && hierarchy->finerLevelExists(level_number))
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
                r[d] = std::max(std::floor(d_r_src[finer_level_number][n]/dx_finer[d]+0.5),2.0)*dx_finer[d];
            }

            // Determine the approximate source stencil box.
            const Index<NDIM> i_center = IndexUtilities::getCellIndex(d_X_src[finer_level_number][n], xLower, xUpper, dx_finer.data(), lower, upper);
            Box<NDIM> stencil_box(i_center,i_center);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                stencil_box.grow(d, static_cast<int>(r[d]/dx_finer[d])+1);
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
    return;
}// applyGradientDetectorSpecialized

void
IBHierarchyIntegrator::putToDatabaseSpecialized(
    Pointer<Database> db)
{
    db->putInteger("IB_HIERARCHY_INTEGRATOR_VERSION",IB_HIERARCHY_INTEGRATOR_VERSION);
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
    db->putBool("d_normalize_source_strength", d_normalize_source_strength);
    db->putBool("d_using_pIB_method", d_using_pIB_method);
    if (d_using_pIB_method)
    {
        db->putDoubleArray("d_gravitational_acceleration", &d_gravitational_acceleration[0], NDIM);
    }
    db->putDouble("d_regrid_cfl_interval", d_regrid_cfl_interval);
    db->putDouble("d_regrid_cfl_estimate", d_regrid_cfl_estimate);
    return;
}// putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBHierarchyIntegrator::resetLagrangianForceFunction(
    const double init_data_time,
    const bool initial_time)
{
    if (d_ib_force_fcn.isNull()) return;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            d_ib_force_fcn->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
        }
    }
    return;
}// resetLagrangianForceFunction

void
IBHierarchyIntegrator::resetLagrangianSourceFunction(
    const double init_data_time,
    const bool initial_time)
{
    if (d_ib_source_fcn.isNull()) return;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_l_data_manager->levelContainsLagrangianData(ln))
        {
            d_ib_source_fcn->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
        }
    }
    return;
}// resetLagrangianSourceFunction

void
IBHierarchyIntegrator::resetPostProcessor(
    const double init_data_time,
    const bool initial_time)
{
    if (d_post_processor.isNull()) return;
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
IBHierarchyIntegrator::updateIBInstrumentationData(
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
    const int U_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getScratchContext());
    const int P_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVariable(), d_ins_hier_integrator->getScratchContext());

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
        d_refine_scheds["INSTRUMENTATION_DATA_FILL"][ln]->fillData(data_time);
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
IBHierarchyIntegrator::resetAnchorPointValues(
    std::vector<Pointer<LData> > V_data,
    const int coarsest_ln,
    const int finest_ln)
{
    PetscErrorCode ierr;
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
            ierr = VecGetArray(V_vec, &V_arr);  IBTK_CHKERRQ(ierr);
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
IBHierarchyIntegrator::computeSourceStrengths(
    const int coarsest_level,
    const int finest_level,
    const double data_time,
    const std::vector<Pointer<LData> >& X_data)
{
    if (d_ib_source_fcn.isNull()) return;

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
            d_ib_source_fcn->getSourceLocations(d_X_src[ln], d_r_src[ln], X_data[ln], d_hierarchy, ln, data_time, d_l_data_manager);
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
            d_ib_source_fcn->computeSourceStrengths(d_Q_src[ln], d_hierarchy, ln, data_time, d_l_data_manager);
        }
    }

    // Spread the sources/sinks onto the Cartesian grid.
    d_hier_pressure_data_ops->setToScalar(d_Q_idx, 0.0);
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
                        r[d] = std::max(std::floor(d_r_src[ln][n]/dx[d]+0.5),2.0)*dx[d];
                    }

                    // Determine the approximate source stencil box.
                    const Index<NDIM> i_center = IndexUtilities::getCellIndex(d_X_src[ln][n], xLower, xUpper, dx, patch_lower, patch_upper);
                    Box<NDIM> stencil_box(i_center,i_center);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        stencil_box.grow(d, static_cast<int>(r[d]/dx[d])+1);
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
    const int wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
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
    const double q_total = d_hier_pressure_data_ops->integral(d_Q_idx, wgt_idx);

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
        TBOX_ERROR(d_object_name << "::integrateHierarchy()  Lagrangian and Eulerian source/sink strengths are inconsistent.");
    }

    // Balance the net inflow/outflow with outflow/inflow along the upper/lower
    // boundaries of the computational domain (if needed).
    if (d_normalize_source_strength)
    {
        if (d_do_log) plog << "    adding ``external'' source/sink to offset net inflow/outflow into domain.\n";
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

        if (std::abs(d_hier_pressure_data_ops->integral(d_Q_idx, wgt_idx) > 1.0e-12))
        {
            TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                       << "  ``external'' source/sink does not correctly offset net inflow/outflow into domain.\n"
                       << "  integral{q} = " << d_hier_pressure_data_ops->integral(d_Q_idx, wgt_idx) << " != 0.\n");
        }
    }
    return;
}// computeSourceStrengths

void
IBHierarchyIntegrator::computeSourcePressures(
    const int coarsest_level,
    const int finest_level,
    const double data_time,
    const std::vector<Pointer<LData> >& X_data)
{
    if (d_ib_source_fcn.isNull()) return;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int P_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVariable(), d_ins_hier_integrator->getNewContext());
    const int wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Compute the normalization pressure (if needed).
    double p_norm = 0.0;
    if (d_normalize_source_strength)
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

            d_ib_source_fcn->getSourceLocations(d_X_src[ln], d_r_src[ln], X_data[ln], d_hierarchy, ln, data_time, d_l_data_manager);
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
                        r[d] = std::max(std::floor(d_r_src[ln][n]/dx[d]+0.5),2.0)*dx[d];
                    }

                    // Determine the approximate source stencil box.
                    const Index<NDIM> i_center = IndexUtilities::getCellIndex(d_X_src[ln][n], xLower, xUpper, dx, patch_lower, patch_upper);
                    Box<NDIM> stencil_box(i_center,i_center);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        stencil_box.grow(d, static_cast<int>(r[d]/dx[d])+1);
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
            d_ib_source_fcn->setSourcePressures(d_P_src[ln], d_hierarchy, ln, data_time, d_l_data_manager);
        }
    }
    return;
}// computeSourcePressures

void
IBHierarchyIntegrator::getFromInput(
    Pointer<Database> db,
    bool is_from_restart)
{
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
        d_normalize_source_strength = db->getBoolWithDefault("normalize_source_strength", d_normalize_source_strength);
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
        d_mark_input_file_name = db->getStringWithDefault("marker_input_file_name", d_mark_input_file_name);
    }
    if (db->keyExists("num_cycles")) d_num_cycles = db->getInteger("num_cycles");
    if (db->keyExists("regrid_cfl_interval")) d_regrid_cfl_interval = db->getDouble("regrid_cfl_interval");
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
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("IB_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
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
    d_normalize_source_strength = db->getBool("d_normalize_source_strength");
    d_using_pIB_method = db->getBoolWithDefault("d_using_pIB_method",false);
    if (d_using_pIB_method)
    {
        db->getDoubleArray("d_gravitational_acceleration", &d_gravitational_acceleration[0], NDIM);
    }
    d_regrid_cfl_interval = db->getDouble("d_regrid_cfl_interval");
    d_regrid_cfl_estimate = db->getDouble("d_regrid_cfl_estimate");
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
