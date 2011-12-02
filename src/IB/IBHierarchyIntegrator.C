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
#include <ibtk/IndexUtilities.h>
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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBHierarchyIntegrator::IBHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<IBStrategy> ib_method_ops,
    Pointer<INSHierarchyIntegrator> ins_hier_integrator,
    bool register_for_restart)
    : HierarchyIntegrator(object_name, input_db, register_for_restart)
{
    // Set the IB method operations objects.
    d_ib_method_ops = ib_method_ops;
    d_ib_method_ops->registerIBHierarchyIntegrator(this);

    // Register the fluid solver as a child integrator of this integrator
    // object.
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!ins_hier_integrator.isNull());
#endif
    d_ins_hier_integrator = ins_hier_integrator;
    registerChildHierarchyIntegrator(d_ins_hier_integrator);

    // Set some default values.
    d_integrator_is_initialized = false;
    d_regrid_cfl_interval = 0.0;
    d_regrid_cfl_estimate = 0.0;
    d_error_on_dt_change = false;
    d_warn_on_dt_change = true;

    // Do not allocate a workload variable by default.
    d_workload_var.setNull();
    d_workload_idx = -1;

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);
    return;
}// IBHierarchyIntegrator

IBHierarchyIntegrator::~IBHierarchyIntegrator()
{
    // intentionally blank.
    return;
}// ~IBHierarchyIntegrator

void
IBHierarchyIntegrator::registerBodyForceFunction(
    Pointer<CartGridFunction> f_fcn)
{
    d_body_force_fcn = f_fcn;
    return;
}// registerBodyForceFunction

void
IBHierarchyIntegrator::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > load_balancer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!load_balancer.isNull());
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
    d_hier_velocity_data_ops    = hier_ops_manager->getOperationsDouble(d_ins_hier_integrator->getVelocityVariable(), hierarchy, true);
    d_hier_pressure_cc_data_ops = hier_ops_manager->getOperationsDouble(d_ins_hier_integrator->getPressureVariable(), hierarchy, true);
    d_hier_cc_data_ops          = hier_ops_manager->getOperationsDouble(new CellVariable<NDIM,double>("cc_var"), hierarchy, true);

    // Initialize all variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const IntVector<NDIM> ib_ghosts = d_ib_method_ops->getMinimumGhostCellWidth();
    const IntVector<NDIM>    ghosts = 1;

    Pointer<CellVariable<NDIM,double> > u_cc_var = d_ins_hier_integrator->getVelocityVariable();
    Pointer<SideVariable<NDIM,double> > u_sc_var = d_ins_hier_integrator->getVelocityVariable();
    if (!u_cc_var.isNull())
    {
        d_u_var = new CellVariable<NDIM,double>(d_object_name+"::u", NDIM);
        d_f_var = new CellVariable<NDIM,double>(d_object_name+"::f", NDIM);
    }
    else if (!u_sc_var.isNull())
    {
        d_u_var = new SideVariable<NDIM,double>(d_object_name+"::u");
        d_f_var = new SideVariable<NDIM,double>(d_object_name+"::f");
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchyIntegrator():\n"
                   << "  unsupported velocity data centering" << std::endl);
    }
    d_u_idx = var_db->registerVariableAndContext(d_u_var, getScratchContext(), ib_ghosts);
    d_f_idx = var_db->registerVariableAndContext(d_f_var, getScratchContext(),    ghosts);

    if (d_ib_method_ops->hasFluidSources())
    {
        Pointer<CellVariable<NDIM,double> > p_cc_var = d_ins_hier_integrator->getPressureVariable();
        if (!p_cc_var.isNull())
        {
            d_p_var = new CellVariable<NDIM,double>(d_object_name+"::p");
            d_q_var = new CellVariable<NDIM,double>(d_object_name+"::q");
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeHierarchyIntegrator():\n"
                       << "  unsupported pressure data centering" << std::endl);
        }
        d_p_idx = var_db->registerVariableAndContext(d_p_var, getScratchContext(), ib_ghosts);
        d_q_idx = var_db->registerVariableAndContext(d_q_var, getScratchContext(),    ghosts);
    }
    else
    {
        d_q_var = NULL;
        d_q_idx = -1;
    }

    // Initialize the objects used to manage Lagrangian-Eulerian interaction.
    d_eulerian_force_fcn = new IBEulerianForceFunction(d_object_name+"::IBEulerianForceFunction", d_f_idx, d_f_idx, d_f_idx);
    d_ins_hier_integrator->registerBodyForceFunction(d_eulerian_force_fcn);
    if (d_ib_method_ops->hasFluidSources())
    {
        d_eulerian_source_fcn = new IBEulerianSourceFunction(d_object_name+"::IBEulerianSourceFunction", d_q_idx, d_q_idx, d_q_idx);
        d_ins_hier_integrator->registerFluidSourceFunction(d_eulerian_source_fcn);
    }

    // Initialize the fluid solver.  It is necessary to do this after
    // registering a body force function or fluid source distribution function.
    d_ins_hier_integrator->initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Create several communications algorithms, used in filling ghost cell data
    // and synchronizing data on the patch hierarchy.
    Pointer<Geometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineAlgorithm<NDIM> > refine_alg;
    Pointer<RefineOperator<NDIM> > refine_op;
    RefinePatchStrategy<NDIM>* refine_patch_strategy;
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg;
    Pointer<CoarsenOperator<NDIM> > coarsen_op;

    const int u_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getNewContext());
    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getScratchContext());
    const int p_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVariable(), d_ins_hier_integrator->getNewContext());
    const int p_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVariable(), d_ins_hier_integrator->getScratchContext());

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
    refine_op = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getVelocityVariable(), "CONSERVATIVE_LINEAR_REFINE");
    refine_alg->registerRefine(u_scratch_idx, u_new_idx, u_scratch_idx, refine_op);
    refine_op = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getPressureVariable(), "LINEAR_REFINE");
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
    }
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                                   d_ins_hier_integrator->getCurrentContext());
    d_hier_velocity_data_ops->copyData(d_u_idx, u_current_idx);
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time,d_start_time);
    d_ib_method_ops->initializePatchHierarchy(hierarchy, gridding_alg, d_u_idx, getCoarsenSchedules(d_object_name+"::u::CONSERVATIVE_COARSEN"), getGhostfillRefineSchedules(d_object_name+"::u"), d_integrator_step, d_integrator_time, initial_time);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_u_idx);
    }

    // Indicate that the hierarchy is initialized.
    d_hierarchy_is_initialized = true;
    return;
}// initializePatchHierarchy

int
IBHierarchyIntegrator::getNumberOfCycles()
{
    return std::max(d_num_cycles,d_ins_hier_integrator->getNumberOfCycles());
}// getNumberOfCycles

void
IBHierarchyIntegrator::preprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const int num_cycles)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time,d_start_time);

    // Determine whether there has been a time step size change.
    if ((d_error_on_dt_change || d_warn_on_dt_change) &&
        (!initial_time && !MathUtilities<double>::equalEps(dt, d_dt_previous[0])))
    {
        if (d_error_on_dt_change)
        {
            TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():  Time step size change encountered.\n"
                       << "Aborting." << std::endl);
        }
        if (d_warn_on_dt_change)
        {
            pout << d_object_name << "::preprocessIntegrateHierarchy():  WARNING: Time step size change encountered.\n"
                 << "Suggest reducing maximum time step size in input file." << std::endl;
        }
    }

    // Setup the Eulerian body force and fluid source/sink functions.
    d_eulerian_force_fcn->registerBodyForceFunction(d_body_force_fcn);
    d_eulerian_force_fcn ->setTimeInterval(current_time, new_time);
    if (d_ib_method_ops->hasFluidSources())
    {
        d_eulerian_source_fcn->setTimeInterval(current_time, new_time);
    }

    // Allocate Eulerian scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_u_idx, current_time);
        level->allocatePatchData(d_f_idx, current_time);
        if (d_ib_method_ops->hasFluidSources())
        {
            level->allocatePatchData(d_q_idx, current_time);
        }
    }

    // Initialize the fluid solver.
    d_ins_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    // Initialize IB data.
    d_ib_method_ops->preprocessIntegrateData(current_time, new_time, num_cycles);

    if (num_cycles == 1)
    {
        // Compute an initial prediction of the updated positions of the
        // Lagrangian structure.
        //
        // NOTE: The velocity should already have been interpolated to the
        // curvilinear mesh and should not need to be re-interpolated.
        if (d_do_log) plog << d_object_name << "::preprocessIntegrateHierarchy(): performing Lagrangian forward Euler step\n";
        d_ib_method_ops->eulerStep(current_time, new_time);
    }
    return;
}// preprocessIntegrateHierarchy

void
IBHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time,
    const int cycle_num)
{
    const double half_time = current_time+0.5*(new_time-current_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                                   d_ins_hier_integrator->getCurrentContext());
    const int u_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                                   d_ins_hier_integrator->getNewContext());
    const int p_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVariable(),
                                                                   d_ins_hier_integrator->getNewContext());

    // Compute the Lagrangian forces and spread them to the Eulerian grid.
    if (d_do_log) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
    d_ib_method_ops->computeLagrangianForce(half_time);
    if (d_do_log) plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
    d_ib_method_ops->spreadForce(d_f_idx, getProlongRefineSchedules(d_object_name+"::f"), half_time);

    // Compute the Lagrangian source/sink strengths and spread them to the
    // Eulerian grid.
    if (d_ib_method_ops->hasFluidSources())
    {
        d_ib_method_ops->computeLagrangianFluidSource(half_time);
        d_hier_pressure_cc_data_ops->setToScalar(d_q_idx, 0.0);
        d_ib_method_ops->spreadFluidSource(d_q_idx, getProlongRefineSchedules(d_object_name+"::q"), half_time);
    }

    // Solve the incompressible Navier-Stokes equations.
    d_ib_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    if (d_do_log) plog << d_object_name << "::integrateHierarchy(): solving the incompressible Navier-Stokes equations\n";
    d_ins_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);
    d_ib_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
    if (d_do_log) plog << d_object_name << "::integrateHierarchy(): interpolating Eulerian velocity to the Lagrangian mesh\n";
    d_ib_method_ops->interpolateVelocity(d_u_idx, getCoarsenSchedules(d_object_name+"::u::CONSERVATIVE_COARSEN"), getGhostfillRefineSchedules(d_object_name+"::u"), half_time);

    // Compute an updated prediction of the updated positions of the Lagrangian
    // structure.
    if (d_current_num_cycles > 1 && d_current_cycle_num == 0)
    {
        if (d_do_log) plog << d_object_name << "::integrateHierarchy(): performing Lagrangian forward Euler step\n";
        d_ib_method_ops->eulerStep(current_time, new_time);
    }
    else
    {
        if (d_do_log) plog << d_object_name << "::integrateHierarchy(): performing Lagrangian midpoint step\n";
        d_ib_method_ops->midpointStep(current_time, new_time);
    }

    // Compute the pressure at the updated locations of any distributed internal
    // fluid sources or sinks.
    if (d_ib_method_ops->hasFluidSources())
    {
        d_hier_pressure_cc_data_ops->copyData(d_p_idx, p_new_idx);
        d_ib_method_ops->interpolatePressure(d_p_idx, getCoarsenSchedules(d_object_name+"::p::CONSERVATIVE_COARSEN"), getGhostfillRefineSchedules(d_object_name+"::p"), half_time);
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
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
    if (d_do_log) plog << d_object_name << "::postprocessIntegrateHierarchy(): interpolating Eulerian velocity to the Lagrangian mesh\n";
    d_ib_method_ops->interpolateVelocity(d_u_idx, getCoarsenSchedules(d_object_name+"::u::CONSERVATIVE_COARSEN"), getGhostfillRefineSchedules(d_object_name+"::u"), new_time);

    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_do_log) plog << d_object_name << "::postprocessIntegrateHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(NEW_DATA);
    }

    // Determine the CFL number.
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
            Pointer<CellData<NDIM,double> > u_cc_new_data = patch->getPatchData(u_new_idx);
            Pointer<SideData<NDIM,double> > u_sc_new_data = patch->getPatchData(u_new_idx);
            double u_max = 0.0;
            if (!u_cc_new_data.isNull()) u_max = patch_cc_ops.maxNorm(u_cc_new_data, patch_box);
            if (!u_sc_new_data.isNull()) u_max = patch_sc_ops.maxNorm(u_sc_new_data, patch_box);
            cfl_max = std::max(cfl_max, u_max*dt/dx_min);
        }
    }
    cfl_max = SAMRAI_MPI::maxReduction(cfl_max);
    d_regrid_cfl_estimate += cfl_max;
    if (d_do_log) plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl_max << "\n";
    if (d_do_log) plog << d_object_name << "::postprocessIntegrateHierarchy(): estimated upper bound on IB point displacement since last regrid = " << d_regrid_cfl_estimate << "\n";

    // Deallocate the fluid solver.
    d_ins_hier_integrator->postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Deallocate IB data.
    d_ib_method_ops->postprocessIntegrateData(current_time, new_time, num_cycles);

    // Deallocate Eulerian scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_u_idx);
        level->deallocatePatchData(d_f_idx);
        if (d_ib_method_ops->hasFluidSources())
        {
            level->deallocatePatchData(d_q_idx);
        }
    }
    return;
}// postprocessIntegrateHierarchy

void
IBHierarchyIntegrator::regridHierarchy()
{
    // Update the workload pre-regridding.
    if (!d_load_balancer.isNull())
    {
        if (d_do_log) plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
        d_hier_cc_data_ops->setToScalar(d_workload_idx, 1.0);
        d_ib_method_ops->updateWorkloadEstimates(d_hierarchy, d_workload_idx);
    }

    // Before regridding, begin Lagrangian data movement.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): starting Lagrangian data movement.\n";
    d_ib_method_ops->beginDataRedistribution(d_hierarchy, d_gridding_alg);

    // Use the INSHierarchyIntegrator to handle Eulerian data management.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): calling " << d_ins_hier_integrator->getName() << "::regridHierarchy().\n";
    d_ins_hier_integrator->regridHierarchy();

    // After regridding, finish Lagrangian data movement.
    if (d_do_log) plog << d_object_name << "::regridHierarchy(): finishing Lagrangian data movement.\n";
    d_ib_method_ops->endDataRedistribution(d_hierarchy, d_gridding_alg);

    // Reset the regrid CFL estimate.
    d_regrid_cfl_estimate = 0.0;
    return;
}// regridHierarchy

/////////////////////////////// PROTECTED ////////////////////////////////////

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
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    // Initialize workload data.
    if (d_workload_idx != -1)
    {
        HierarchyCellDataOpsReal<NDIM,double> level_cc_data_ops(hierarchy,level_number,level_number);
        level_cc_data_ops.setToScalar(d_workload_idx, 1.0);
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
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0) && (coarsest_level <= finest_level) && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset IB data.
    d_ib_method_ops->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_velocity_data_ops   ->setPatchHierarchy(hierarchy);
    d_hier_pressure_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops         ->setPatchHierarchy(hierarchy);
    d_hier_velocity_data_ops   ->resetLevels(0, finest_hier_level);
    d_hier_pressure_cc_data_ops->resetLevels(0, finest_hier_level);
    d_hier_cc_data_ops         ->resetLevels(0, finest_hier_level);
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
    if (db->keyExists("regrid_cfl_interval")) d_regrid_cfl_interval = db->getDouble("regrid_cfl_interval");
    if      (db->keyExists("error_on_dt_change")       ) d_error_on_dt_change = db->getBool("error_on_dt_change");
    else if (db->keyExists("error_on_timestep_change") ) d_error_on_dt_change = db->getBool("error_on_timestep_change");
    else if (db->keyExists("error_on_time_step_change")) d_error_on_dt_change = db->getBool("error_on_time_step_change");
    if      (db->keyExists("warn_on_dt_change")       ) d_warn_on_dt_change = db->getBool("warn_on_dt_change");
    else if (db->keyExists("warn_on_timestep_change") ) d_warn_on_dt_change = db->getBool("warn_on_timestep_change");
    else if (db->keyExists("warn_on_time_step_change")) d_warn_on_dt_change = db->getBool("warn_on_time_step_change");
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
    d_regrid_cfl_interval = db->getDouble("d_regrid_cfl_interval");
    d_regrid_cfl_estimate = db->getDouble("d_regrid_cfl_estimate");
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
