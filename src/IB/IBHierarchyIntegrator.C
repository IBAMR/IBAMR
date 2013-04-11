// Filename: IBHierarchyIntegrator.C
// Created on 12 Jul 2004 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/LMarkerUtilities.h>

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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBHierarchyIntegrator::~IBHierarchyIntegrator()
{
    // intentionally blank
    return;
}// ~IBHierarchyIntegrator

Pointer<IBStrategy>
IBHierarchyIntegrator::getIBStrategy() const
{
    return d_ib_method_ops;
}// getIBStrategy

void
IBHierarchyIntegrator::registerBodyForceFunction(
    Pointer<CartGridFunction> f_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    if (d_body_force_fcn)
    {
        Pointer<CartGridFunctionSet> p_body_force_fcn = d_body_force_fcn;
        if (!p_body_force_fcn)
        {
            pout << d_object_name << "::registerBodyForceFunction(): WARNING:\n"
                 << "  body force function has already been set.\n"
                 << "  functions will be evaluated in the order in which they were registered with the solver\n"
                 << "  when evaluating the body force term value.\n";
            p_body_force_fcn = new CartGridFunctionSet(d_object_name+"::body_force_function_set");
            p_body_force_fcn->addFunction(d_body_force_fcn);
        }
        p_body_force_fcn->addFunction(f_fcn);
    }
    else
    {
        d_body_force_fcn = f_fcn;
    }
    return;
}// registerBodyForceFunction

void
IBHierarchyIntegrator::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > load_balancer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(load_balancer);
#endif
    d_load_balancer = load_balancer;
    if (d_workload_idx == -1)
    {
        d_workload_var = new CellVariable<NDIM,double>(d_object_name+"::workload");
        registerVariable(d_workload_idx, d_workload_var, 0, getCurrentContext());
    }
    d_ib_method_ops->registerLoadBalancer(load_balancer, d_workload_idx);
    return;
}// registerLoadBalancer

Pointer<Variable<NDIM> >
IBHierarchyIntegrator::getVelocityVariable() const
{
    return d_u_var;
}// getVelocityVariable

Pointer<Variable<NDIM> >
IBHierarchyIntegrator::getPressureVariable() const
{
    return d_p_var;
}// getPressureVariable

Pointer<Variable<NDIM> >
IBHierarchyIntegrator::getBodyForceVariable() const
{
    return d_f_var;
}// getBodyForceVariable

Pointer<Variable<NDIM> >
IBHierarchyIntegrator::getFluidSourceVariable() const
{
    return d_q_var;
}// getFluidSourceVariable

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
    d_hier_velocity_data_ops = hier_ops_manager->getOperationsDouble(getVelocityVariable(), hierarchy, true);
    d_hier_pressure_data_ops = hier_ops_manager->getOperationsDouble(getPressureVariable(), hierarchy, true);
    d_hier_cc_data_ops       = hier_ops_manager->getOperationsDouble(new CellVariable<NDIM,double>("cc_var"), hierarchy, true);

    // Initialize all variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const IntVector<NDIM> ib_ghosts = d_ib_method_ops->getMinimumGhostCellWidth();
    const IntVector<NDIM>    ghosts = 1;

    d_u_idx = var_db->registerVariableAndContext(d_u_var, d_ib_context, ib_ghosts);
    d_f_idx = var_db->registerVariableAndContext(d_f_var, d_ib_context,    ghosts);
    if (d_time_stepping_type == TRAPEZOIDAL_RULE)
    {
        d_f_current_idx = var_db->registerClonedPatchDataIndex(d_f_var, d_f_idx);
    }
    else
    {
        d_f_current_idx = -1;
    }

    if (d_ib_method_ops->hasFluidSources())
    {
        d_p_idx = var_db->registerVariableAndContext(d_p_var, d_ib_context, ib_ghosts);
        d_q_idx = var_db->registerVariableAndContext(d_q_var, d_ib_context,    ghosts);
    }
    else
    {
        d_q_var = NULL;
        d_q_idx = -1;
    }

    if (!d_mark_file_name.empty())
    {
        d_mark_var = new LMarkerSetVariable(d_object_name+"::markers");
        registerVariable(d_mark_current_idx, d_mark_new_idx, d_mark_scratch_idx, d_mark_var, ghosts);
    }

    // Initialize the fluid solver.
    if (d_ib_method_ops->hasFluidSources())
    {
        d_ins_hier_integrator->registerFluidSourceFunction(new IBEulerianSourceFunction(this));
    }
    d_ins_hier_integrator->initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Have the IB method ops object register any additional Eulerian variables
    // and communications algorithms that it requires.
    d_ib_method_ops->registerEulerianVariables();
    d_ib_method_ops->registerEulerianCommunicationAlgorithms();

    // Create several communications algorithms, used in filling ghost cell data
    // and synchronizing data on the patch hierarchy.
    Pointer<Geometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineAlgorithm<NDIM> > refine_alg;
    Pointer<RefineOperator<NDIM> > refine_op;
    RefinePatchStrategy<NDIM>* refine_patch_strategy;
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg;
    Pointer<CoarsenOperator<NDIM> > coarsen_op;

    const int u_new_idx     = var_db->mapVariableAndContextToIndex(getVelocityVariable(), getNewContext());
    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(getVelocityVariable(), getScratchContext());
    const int p_new_idx     = var_db->mapVariableAndContextToIndex(getPressureVariable(), getNewContext());
    const int p_scratch_idx = var_db->mapVariableAndContextToIndex(getPressureVariable(), getScratchContext());

    refine_alg = new RefineAlgorithm<NDIM>();
    refine_op = NULL;
    refine_alg->registerRefine(d_u_idx, d_u_idx, d_u_idx, refine_op);
    registerGhostfillRefineAlgorithm(d_object_name+"::u", refine_alg);

    coarsen_alg = new CoarsenAlgorithm<NDIM>();
    coarsen_op = grid_geom->lookupCoarsenOperator(d_u_var, "CONSERVATIVE_COARSEN");
    coarsen_alg->registerCoarsen(d_u_idx, d_u_idx, coarsen_op);
    registerCoarsenAlgorithm(d_object_name+"::u::CONSERVATIVE_COARSEN", coarsen_alg);

    refine_alg = new RefineAlgorithm<NDIM>();
    refine_op = grid_geom->lookupRefineOperator(d_f_var, "CONSERVATIVE_LINEAR_REFINE");
    refine_alg->registerRefine(d_f_idx, d_f_idx, d_f_idx, refine_op);
    registerProlongRefineAlgorithm(d_object_name+"::f", refine_alg);

    refine_alg = new RefineAlgorithm<NDIM>();
    refine_op = grid_geom->lookupRefineOperator(getVelocityVariable(), "CONSERVATIVE_LINEAR_REFINE");
    refine_alg->registerRefine(u_scratch_idx, u_new_idx, u_scratch_idx, refine_op);
    refine_op = grid_geom->lookupRefineOperator(getPressureVariable(), "LINEAR_REFINE");
    refine_alg->registerRefine(p_scratch_idx, p_new_idx, p_scratch_idx, refine_op);
    ComponentSelector instrumentation_data_fill_bc_idxs;
    instrumentation_data_fill_bc_idxs.setFlag(u_scratch_idx);
    instrumentation_data_fill_bc_idxs.setFlag(p_scratch_idx);
    refine_patch_strategy = new CartExtrapPhysBdryOp(instrumentation_data_fill_bc_idxs, "LINEAR");
    registerGhostfillRefineAlgorithm(d_object_name+"::INSTRUMENTATION_DATA_FILL", refine_alg, refine_patch_strategy);

    if (d_ib_method_ops->hasFluidSources())
    {
        refine_alg = new RefineAlgorithm<NDIM>();
        refine_op = NULL;
        refine_alg->registerRefine(d_p_idx, d_p_idx, d_p_idx, refine_op);
        registerGhostfillRefineAlgorithm(d_object_name+"::p", refine_alg);

        coarsen_alg = new CoarsenAlgorithm<NDIM>();
        coarsen_op = grid_geom->lookupCoarsenOperator(d_p_var, "CONSERVATIVE_COARSEN");
        coarsen_alg->registerCoarsen(d_p_idx, d_p_idx, coarsen_op);
        registerCoarsenAlgorithm(d_object_name+"::p::CONSERVATIVE_COARSEN", coarsen_alg);

        refine_alg = new RefineAlgorithm<NDIM>();
        refine_op = grid_geom->lookupRefineOperator(d_q_var, "CONSERVATIVE_LINEAR_REFINE");
        refine_alg->registerRefine(d_q_idx, d_q_idx, d_q_idx, refine_op);
        registerProlongRefineAlgorithm(d_object_name+"::q", refine_alg);
    }

    // Read in initial marker positions.
    if (!d_mark_file_name.empty())
    {
        LMarkerUtilities::readMarkerPositions(d_mark_init_posns, d_mark_file_name, hierarchy->getGridGeometry());
    }

    // Setup the tag buffer.
    const int finest_hier_ln = gridding_alg->getMaxLevels()-1;
    const int tsize = d_tag_buffer.size();
    d_tag_buffer.resizeArray(finest_hier_ln);
    for (int i = tsize; i < finest_hier_ln; ++i) d_tag_buffer[i] = 0;
    for (int i = std::max(tsize,1); i < d_tag_buffer.size(); ++i)
    {
        d_tag_buffer[i] = d_tag_buffer[i-1];
    }
    d_ib_method_ops->setupTagBuffer(d_tag_buffer, d_gridding_alg);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
}// initializeHierarchyIntegrator

void
IBHierarchyIntegrator::initializePatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_hierarchy_is_initialized) return;

    // Initialize Eulerian data.
    HierarchyIntegrator::initializePatchHierarchy(hierarchy, gridding_alg);

    // Begin Lagrangian data movement.
    d_ib_method_ops->beginDataRedistribution(hierarchy, gridding_alg);

    // Finish Lagrangian data movement.
    d_ib_method_ops->endDataRedistribution(hierarchy, gridding_alg);

    // Initialize Lagrangian data on the patch hierarchy.
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_u_idx, d_integrator_time);
        level->allocatePatchData(d_scratch_data, d_integrator_time);
    }
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(getVelocityVariable(), getCurrentContext());
    d_hier_velocity_data_ops->copyData(d_u_idx, u_current_idx);
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time,d_start_time);
    d_ib_method_ops->initializePatchHierarchy(hierarchy, gridding_alg, d_u_idx, getCoarsenSchedules(d_object_name+"::u::CONSERVATIVE_COARSEN"), getGhostfillRefineSchedules(d_object_name+"::u"), d_integrator_step, d_integrator_time, initial_time);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_u_idx);
        level->deallocatePatchData(d_scratch_data);
    }

    // Indicate that the hierarchy is initialized.
    d_hierarchy_is_initialized = true;
    return;
}// initializePatchHierarchy

void
IBHierarchyIntegrator::regridHierarchy()
{
    // Update the workload pre-regridding.
    if (d_load_balancer)
    {
        if (d_enable_logging) plog << d_object_name << "::regridHierarchy(): updating workload estimates\n";
        d_hier_cc_data_ops->setToScalar(d_workload_idx, 1.0);
        d_ib_method_ops->updateWorkloadEstimates(d_hierarchy, d_workload_idx);
    }

    // Collect the marker particles to level 0 of the patch hierarchy.
    if (d_mark_var)
    {
        LMarkerUtilities::collectMarkersOnPatchHierarchy(d_mark_current_idx, d_hierarchy);
    }

    // Before regridding, begin Lagrangian data movement.
    if (d_enable_logging) plog << d_object_name << "::regridHierarchy(): starting Lagrangian data movement\n";
    d_ib_method_ops->beginDataRedistribution(d_hierarchy, d_gridding_alg);

    // Use the INSHierarchyIntegrator to handle Eulerian data management.
    if (d_enable_logging) plog << d_object_name << "::regridHierarchy(): regridding the patch hierarchy\n";
    HierarchyIntegrator::regridHierarchy();

    // After regridding, finish Lagrangian data movement.
    if (d_enable_logging) plog << d_object_name << "::regridHierarchy(): finishing Lagrangian data movement\n";
    d_ib_method_ops->endDataRedistribution(d_hierarchy, d_gridding_alg);

    // Prune any duplicated markers located in the "invalid" regions of coarser
    // levels of the patch hierarchy.
    if (d_mark_var)
    {
        LMarkerUtilities::pruneInvalidMarkers(d_mark_current_idx, d_hierarchy);
    }

    // Reset the regrid CFL estimate.
    d_regrid_cfl_estimate = 0.0;
    return;
}// regridHierarchy

/////////////////////////////// PROTECTED ////////////////////////////////////

IBHierarchyIntegrator::IBHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<IBStrategy> ib_method_ops,
    Pointer<INSHierarchyIntegrator> ins_hier_integrator,
    bool register_for_restart)
    : HierarchyIntegrator(object_name, input_db, register_for_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(ib_method_ops);
    TBOX_ASSERT(ins_hier_integrator);
#endif

    // Set the IB method operations objects.
    d_ib_method_ops = ib_method_ops;
    d_ib_method_ops->registerIBHierarchyIntegrator(this);

    // Register the fluid solver as a child integrator of this integrator object
    // and reuse the variables and variable contexts of the INS solver.
    d_ins_hier_integrator = ins_hier_integrator;
    registerChildHierarchyIntegrator(d_ins_hier_integrator);
    d_u_var = d_ins_hier_integrator->getVelocityVariable();
    d_p_var = d_ins_hier_integrator->getPressureVariable();
    d_f_var = d_ins_hier_integrator->getBodyForceVariable();
    d_q_var = d_ins_hier_integrator->getFluidSourceVariable();
    d_current_context = d_ins_hier_integrator->getCurrentContext();
    d_scratch_context = d_ins_hier_integrator->getScratchContext();
    d_new_context = d_ins_hier_integrator->getNewContext();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_ib_context = var_db->getContext(d_object_name+"::IB");

    // Set some default values.
    d_integrator_is_initialized = false;
    d_time_stepping_type = MIDPOINT_RULE;
    d_regrid_cfl_interval = 0.0;
    d_regrid_cfl_estimate = 0.0;
    d_error_on_dt_change = true;
    d_warn_on_dt_change = false;

    // Do not allocate a workload variable by default.
    d_workload_var.setNull();
    d_workload_idx = -1;

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);
    return;
}// IBHierarchyIntegrator

bool
IBHierarchyIntegrator::atRegridPointSpecialized() const
{
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time,d_start_time);
    if (initial_time) return true;
    if (d_regrid_cfl_interval > 0.0)
    {
        return (d_regrid_cfl_estimate >= d_regrid_cfl_interval);
    }
    else if (d_regrid_interval != 0)
    {
        return (d_integrator_step % d_regrid_interval == 0);
    }
    return false;
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
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (old_level)
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif

    // Initialize workload data.
    if (d_workload_idx != -1)
    {
        HierarchyCellDataOpsReal<NDIM,double> level_cc_data_ops(hierarchy,level_number,level_number);
        level_cc_data_ops.setToScalar(d_workload_idx, 1.0);
        d_load_balancer->setUniformWorkload(level_number);
    }

    // Initialize marker data
    if (d_mark_var)
    {
        LMarkerUtilities::initializeMarkersOnLevel(d_mark_current_idx, d_mark_init_posns, hierarchy, level_number, initial_time, old_level);
    }

    // Initialize IB data.
    d_ib_method_ops->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
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
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((coarsest_level >= 0) && (coarsest_level <= finest_level) && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(hierarchy->getPatchLevel(ln));
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset IB data.
    d_ib_method_ops->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_velocity_data_ops->setPatchHierarchy(hierarchy);
    d_hier_pressure_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops      ->setPatchHierarchy(hierarchy);
    d_hier_velocity_data_ops->resetLevels(0, finest_hier_level);
    d_hier_pressure_data_ops->resetLevels(0, finest_hier_level);
    d_hier_cc_data_ops      ->resetLevels(0, finest_hier_level);
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
    // Tag cells for refinement.
    d_ib_method_ops->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    return;
}// applyGradientDetectorSpecialized

void
IBHierarchyIntegrator::putToDatabaseSpecialized(
    Pointer<Database> db)
{
    db->putInteger("IB_HIERARCHY_INTEGRATOR_VERSION",IB_HIERARCHY_INTEGRATOR_VERSION);
    db->putString("d_time_stepping_type", enum_to_string<TimeSteppingType>(d_time_stepping_type));
    db->putDouble("d_regrid_cfl_interval", d_regrid_cfl_interval);
    db->putDouble("d_regrid_cfl_estimate", d_regrid_cfl_estimate);
    return;
}// putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBHierarchyIntegrator::getFromInput(
    Pointer<Database> db,
    bool /*is_from_restart*/)
{
    if      (db->keyExists("regrid_cfl_interval")) d_regrid_cfl_interval = db->getDouble("regrid_cfl_interval");
    if      (db->keyExists("error_on_dt_change")       ) d_error_on_dt_change = db->getBool("error_on_dt_change");
    else if (db->keyExists("error_on_timestep_change") ) d_error_on_dt_change = db->getBool("error_on_timestep_change");
    else if (db->keyExists("error_on_time_step_change")) d_error_on_dt_change = db->getBool("error_on_time_step_change");
    if      (db->keyExists("warn_on_dt_change")       ) d_warn_on_dt_change = db->getBool("warn_on_dt_change");
    else if (db->keyExists("warn_on_time_step_change")) d_warn_on_dt_change = db->getBool("warn_on_time_step_change");
    else if (db->keyExists("warn_on_timestep_change") ) d_warn_on_dt_change = db->getBool("warn_on_timestep_change");
    if      (db->keyExists("time_stepping_type")) d_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("time_stepping_type"));
    else if (db->keyExists("timestepping_type") ) d_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("timestepping_type" ));
    if      (db->keyExists("marker_file_name")) d_mark_file_name = db->getString("marker_file_name");
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
    d_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("d_time_stepping_type"));
    d_regrid_cfl_interval = db->getDouble("d_regrid_cfl_interval");
    d_regrid_cfl_estimate = db->getDouble("d_regrid_cfl_estimate");
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
