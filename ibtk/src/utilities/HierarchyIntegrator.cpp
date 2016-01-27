// Filename: HierarchyIntegrator.cpp
// Created on 10 Aug 2011 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <algorithm>
#include <deque>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CellData.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenPatchStrategy.h"
#include "CoarsenSchedule.h"
#include "ComponentSelector.h"
#include "EdgeData.h"
#include "FaceData.h"
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "NodeData.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "TagAndInitializeStrategy.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"
#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/RefinePatchStrategySet.h"
#include "ibtk/ibtk_enums.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of HierarchyIntegrator restart file data.
static const int HIERARCHY_INTEGRATOR_VERSION = 1;
}

const std::string HierarchyIntegrator::SYNCH_CURRENT_DATA_ALG = "SYNCH_CURRENT_DATA";
const std::string HierarchyIntegrator::SYNCH_NEW_DATA_ALG = "SYNCH_NEW_DATA";

/////////////////////////////// PUBLIC ///////////////////////////////////////

HierarchyIntegrator::HierarchyIntegrator(const std::string& object_name,
                                         Pointer<Database> input_db,
                                         bool register_for_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
#endif
    d_object_name = object_name;
    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
    }

    // Initialize state variables.
    d_hierarchy = NULL;
    d_gridding_alg = NULL;
    d_hierarchy_is_initialized = false;
    d_regridding_hierarchy = false;
    d_at_regrid_time_step = false;
    d_visit_writer = NULL;
    d_parent_integrator = NULL;
    d_current_num_cycles = -1;
    d_current_cycle_num = -1;
    d_current_dt = std::numeric_limits<double>::quiet_NaN();

    // Set default values.
    d_integrator_time = std::numeric_limits<double>::quiet_NaN();
    d_start_time = 0.0;
    d_end_time = std::numeric_limits<double>::max();
    d_dt_init = std::numeric_limits<double>::max();
    d_dt_min = 0.0;
    d_dt_max = std::numeric_limits<double>::max();
    d_dt_growth_factor = 2.0;
    d_integrator_step = 0;
    d_num_cycles = 1;
    d_max_integrator_steps = std::numeric_limits<int>::max();
    d_regrid_interval = 1;
    d_regrid_mode = STANDARD;
    d_enable_logging = false;
    d_bdry_extrap_type = "LINEAR";
    d_manage_hier_math_ops = true;
    d_tag_buffer.resizeArray(1);
    d_tag_buffer[0] = 0;

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    else
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }
    if (input_db) getFromInput(input_db, from_restart);

    // Initialize all variable contexts.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_current_context = var_db->getContext(d_object_name + "::CURRENT");
    d_new_context = var_db->getContext(d_object_name + "::NEW");
    d_scratch_context = var_db->getContext(d_object_name + "::SCRATCH");

    // Create default communications algorithms.
    d_coarsen_algs[SYNCH_CURRENT_DATA_ALG] = new CoarsenAlgorithm<NDIM>();
    d_coarsen_algs[SYNCH_NEW_DATA_ALG] = new CoarsenAlgorithm<NDIM>();
    d_fill_after_regrid_phys_bdry_bc_op = NULL;
    return;
} // HierarchyIntegrator

HierarchyIntegrator::~HierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }

    for (RefinePatchStrategyMap::iterator it = d_ghostfill_strategies.begin(); it != d_ghostfill_strategies.end(); ++it)
    {
        delete it->second;
    }

    for (RefinePatchStrategyMap::iterator it = d_prolong_strategies.begin(); it != d_prolong_strategies.end(); ++it)
    {
        delete it->second;
    }

    for (CoarsenPatchStrategyMap::iterator it = d_coarsen_strategies.begin(); it != d_coarsen_strategies.end(); ++it)
    {
        delete it->second;
    }
    return;
} // ~HierarchyIntegrator

const std::string&
HierarchyIntegrator::getName() const
{
    return d_object_name;
} // getName

void
HierarchyIntegrator::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                              Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_hierarchy_is_initialized || d_parent_integrator) return;

#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT(gridding_alg);
#endif
    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Initialize the hierarchy integrator and any child hierarchy integrators
    // prior to initializing the patch hierarchy.
    std::deque<HierarchyIntegrator*> hier_integrators(1, this);
    while (!hier_integrators.empty())
    {
        HierarchyIntegrator* integrator = hier_integrators.front();
        integrator->initializeHierarchyIntegrator(d_hierarchy, d_gridding_alg);
        integrator->setupTagBuffer(d_gridding_alg);
        for (int i = 0; i < std::min(d_tag_buffer.size(), integrator->d_tag_buffer.size()); ++i)
        {
            d_tag_buffer[i] = std::max(d_tag_buffer[i], integrator->d_tag_buffer[i]);
        }
        hier_integrators.pop_front();
        hier_integrators.insert(
            hier_integrators.end(), integrator->d_child_integrators.begin(), integrator->d_child_integrators.end());
    }
    plog << d_object_name << "::initializePatchHierarchy(): "
         << "tag_buffer =";
    for (int i = 0; i < d_tag_buffer.size(); ++i) plog << " " << d_tag_buffer[i];
    plog << "\n";

    // Initialize the patch hierarchy.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        d_hierarchy->getFromRestart(d_gridding_alg->getMaxLevels());
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        d_gridding_alg->getTagAndInitializeStrategy()->resetHierarchyConfiguration(d_hierarchy, coarsest_ln, finest_ln);
    }
    else
    {
        const bool initial_time = true;
        d_gridding_alg->makeCoarsestLevel(d_hierarchy, d_start_time);
        int level_number = 0;
        bool done = false;
        while (!done && (d_gridding_alg->levelCanBeRefined(level_number)))
        {
            d_gridding_alg->makeFinerLevel(d_hierarchy, d_integrator_time, initial_time, d_tag_buffer[level_number]);
            done = !d_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }
        synchronizeHierarchyData(CURRENT_DATA);
    }

    // Indicate that the hierarchy is initialized.
    hier_integrators.push_back(this);
    while (!hier_integrators.empty())
    {
        HierarchyIntegrator* integrator = hier_integrators.front();
        integrator->d_hierarchy_is_initialized = true;
        hier_integrators.pop_front();
        hier_integrators.insert(
            hier_integrators.end(), integrator->d_child_integrators.begin(), integrator->d_child_integrators.end());
    }
    return;
} // initializePatchHierarchy

void
HierarchyIntegrator::advanceHierarchy(double dt)
{
    const double dt_min = getMinimumTimeStepSize();
    const double dt_max = getMaximumTimeStepSize();
    if (dt < dt_min || dt > dt_max)
    {
        TBOX_ERROR(d_object_name << "::advanceHierarchy():\n"
                                 << "  at time = "
                                 << d_integrator_time
                                 << ": time step size dt = "
                                 << dt
                                 << "\n"
                                 << "  minimum time step size = "
                                 << dt_min
                                 << "\n"
                                 << "  maximum time step size = "
                                 << dt_max
                                 << "\n");
    }

    if (d_integrator_time + dt > d_end_time)
    {
        pout << "WARNING: at time = " << d_integrator_time << ": reducing dt so that current_time+dt = end_time.\n";
        dt = d_end_time - d_integrator_time;
    }
    const double current_time = d_integrator_time;
    const double new_time = d_integrator_time + dt;
    if (dt < 0.0)
    {
        TBOX_ERROR(d_object_name << "::advanceHierarchy():\n"
                                 << "  at time = "
                                 << d_integrator_time
                                 << ": time step size dt = "
                                 << dt
                                 << ".\n");
    }
    else if (dt == 0.0)
    {
        TBOX_ERROR(d_object_name << "::advanceHierarchy():\n"
                                 << "  at time = "
                                 << d_integrator_time
                                 << ": time step size dt = "
                                 << dt
                                 << ".\n");
    }
    else if (current_time == new_time || MathUtilities<double>::equalEps(current_time, new_time))
    {
        TBOX_ERROR(d_object_name << "::advanceHierarchy():\n"
                                 << "  at time = "
                                 << d_integrator_time
                                 << ": time step size dt = "
                                 << dt
                                 << " is zero to machine precision.\n");
    }
    if (d_enable_logging)
        plog << d_object_name << "::advanceHierarchy(): time interval = [" << current_time << "," << new_time
             << "], dt = " << dt << "\n";

    // Regrid the patch hierarchy.
    if (atRegridPoint())
    {
        if (d_enable_logging)
            plog << d_object_name << "::advanceHierarchy(): regridding prior to timestep " << d_integrator_step << "\n";
        d_regridding_hierarchy = true;
        regridHierarchy();
        d_regridding_hierarchy = false;
        d_at_regrid_time_step = true;
    }

    // Determine the number of cycles and the time step size.
    d_current_num_cycles = getNumberOfCycles();
    d_current_dt = new_time - current_time;

    // Execute the preprocessing method of the parent integrator, and
    // recursively execute all preprocessing callbacks registered with the
    // parent and child integrators.
    preprocessIntegrateHierarchy(current_time, new_time, d_current_num_cycles);

    // Perform one or more cycles.  In each cycle, execute the integration
    // method of the parent integrator, and recursively execute all integration
    // callbacks registered with the parent and child integrators.
    if (d_enable_logging) plog << d_object_name << "::advanceHierarchy(): integrating hierarchy\n";
    for (int cycle_num = 0; cycle_num < d_current_num_cycles; ++cycle_num)
    {
        if (d_enable_logging && d_current_num_cycles != 1)
        {
            if (d_enable_logging)
                plog << d_object_name << "::advanceHierarchy(): executing cycle " << cycle_num + 1 << " of "
                     << d_current_num_cycles << "\n";
        }
        integrateHierarchy(current_time, new_time, cycle_num);
    }

    // Execute the postprocessing method of the parent integrator, and
    // recursively execute all postprocessing callbacks registered with the
    // parent and child integrators.
    static const bool skip_synchronize_new_state_data = true;
    postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data, d_current_num_cycles);

    // Ensure that the current values of num_cycles, cycle_num, and dt are
    // reset.
    std::deque<HierarchyIntegrator*> hier_integrators(1, this);
    while (!hier_integrators.empty())
    {
        HierarchyIntegrator* integrator = hier_integrators.front();
        integrator->d_current_num_cycles = -1;
        integrator->d_current_cycle_num = -1;
        integrator->d_current_dt = std::numeric_limits<double>::quiet_NaN();
        hier_integrators.pop_front();
        hier_integrators.insert(
            hier_integrators.end(), integrator->d_child_integrators.begin(), integrator->d_child_integrators.end());
    }

    // Synchronize the updated data.
    if (d_enable_logging) plog << d_object_name << "::advanceHierarchy(): synchronizing updated data\n";
    synchronizeHierarchyData(NEW_DATA);

    // Reset all time dependent data.
    if (d_enable_logging) plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierarchyData(new_time);

    // Reset the regrid indicator.
    d_at_regrid_time_step = false;
    return;
} // advanceHierarchy

double
HierarchyIntegrator::getMinimumTimeStepSize()
{
    double dt = getMinimumTimeStepSizeSpecialized();
    for (std::set<HierarchyIntegrator*>::iterator it = d_child_integrators.begin(); it != d_child_integrators.end();
         ++it)
    {
        dt = std::max(dt, (*it)->getMinimumTimeStepSize());
    }
    return std::min(dt, d_end_time - d_integrator_time);
} // getMinimumTimeStepSize

double
HierarchyIntegrator::getMaximumTimeStepSize()
{
    double dt = getMaximumTimeStepSizeSpecialized();
    for (std::set<HierarchyIntegrator*>::iterator it = d_child_integrators.begin(); it != d_child_integrators.end();
         ++it)
    {
        dt = std::min(dt, (*it)->getMaximumTimeStepSize());
    }
    return std::min(dt, d_end_time - d_integrator_time);
} // getMaximumTimeStepSize

void
HierarchyIntegrator::synchronizeHierarchyData(VariableContextType ctx_type)
{
    synchronizeHierarchyDataSpecialized(ctx_type);
    for (std::set<HierarchyIntegrator*>::iterator it = d_child_integrators.begin(); it != d_child_integrators.end();
         ++it)
    {
        (*it)->synchronizeHierarchyData(ctx_type);
    }
    return;
} // synchronizeHierarchyData

void
HierarchyIntegrator::resetTimeDependentHierarchyData(const double new_time)
{
    resetTimeDependentHierarchyDataSpecialized(new_time);
    for (std::set<HierarchyIntegrator*>::iterator it = d_child_integrators.begin(); it != d_child_integrators.end();
         ++it)
    {
        (*it)->resetTimeDependentHierarchyData(new_time);
    }
    return;
} // resetTimeDependentHierarchyData

void
HierarchyIntegrator::resetIntegratorToPreadvanceState()
{
    resetIntegratorToPreadvanceStateSpecialized();
    for (std::set<HierarchyIntegrator*>::iterator it = d_child_integrators.begin(); it != d_child_integrators.end();
         ++it)
    {
        (*it)->resetIntegratorToPreadvanceState();
    }
    return;
} // resetIntegratorToPreadvanceState

void
HierarchyIntegrator::regridHierarchy()
{
    const int coarsest_ln = 0;

    bool check_volume_change = !d_parent_integrator && d_hierarchy_is_initialized;
    const double old_volume = check_volume_change ? d_hier_math_ops->getVolumeOfPhysicalDomain() : 0.0;

    // Regrid the hierarchy.
    switch (d_regrid_mode)
    {
    case STANDARD:
        d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln, d_integrator_time, d_tag_buffer);
        break;
    case AGGRESSIVE:
        for (int k = 0; k < d_gridding_alg->getMaxLevels(); ++k)
        {
            d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln, d_integrator_time, d_tag_buffer);
        }
        break;
    default:
        TBOX_ERROR(d_object_name << "::regridHierarchy():\n"
                                 << "  unrecognized regrid mode: "
                                 << enum_to_string<RegridMode>(d_regrid_mode)
                                 << "."
                                 << std::endl);
    }

    const double new_volume = check_volume_change ? d_hier_math_ops->getVolumeOfPhysicalDomain() : 0.0;

    if (check_volume_change && !MathUtilities<double>::equalEps(old_volume, new_volume))
    {
        TBOX_WARNING(
            d_object_name << "::regridHierarchy():\n"
                          << "  change in domain volume detected (volume is computed by summing control volumes)\n"
                          << "    old volume = "
                          << old_volume
                          << "\n"
                          << "    new volume = "
                          << new_volume
                          << "\n"
                          << "  this may indicate overlapping patches in the AMR grid hierarchy.");
    }

    // Synchronize the state data on the patch hierarchy.
    synchronizeHierarchyData(CURRENT_DATA);
    return;
} // regridHierarchy

bool
HierarchyIntegrator::atRegridPoint() const
{
    bool regrid_hierarchy = atRegridPointSpecialized();
    for (std::set<HierarchyIntegrator*>::iterator it = d_child_integrators.begin();
         it != d_child_integrators.end() && !regrid_hierarchy;
         ++it)
    {
        regrid_hierarchy = regrid_hierarchy || (*it)->atRegridPoint();
    }
    return regrid_hierarchy;
} // atRegridPoint

double
HierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
} // getIntegratorTime

double
HierarchyIntegrator::getStartTime() const
{
    return d_start_time;
} // getStartTime

double
HierarchyIntegrator::getEndTime() const
{
    return d_end_time;
} // getEndTime

int
HierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
} // getIntegratorStep

int
HierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
} // getMaxIntegratorSteps

bool
HierarchyIntegrator::stepsRemaining() const
{
    return ((d_integrator_step < d_max_integrator_steps) && (d_integrator_time < d_end_time) &&
            !MathUtilities<double>::equalEps(d_integrator_time, d_end_time));
} // stepsRemaining

Pointer<PatchHierarchy<NDIM> >
HierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
} // getPatchHierarchy

Pointer<GriddingAlgorithm<NDIM> >
HierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
} // getGriddingAlgorithm

void
HierarchyIntegrator::registerVisItDataWriter(Pointer<VisItDataWriter<NDIM> > visit_writer)
{
    d_visit_writer = visit_writer;
    for (std::set<HierarchyIntegrator*>::iterator it = d_child_integrators.begin(); it != d_child_integrators.end();
         ++it)
    {
        (*it)->registerVisItDataWriter(visit_writer);
    }
    return;
} // registerVisItDataWriter

Pointer<VisItDataWriter<NDIM> >
HierarchyIntegrator::getVisItDataWriter() const
{
    return d_visit_writer;
}

void
HierarchyIntegrator::setupPlotData()
{
    setupPlotDataSpecialized();
    for (std::set<HierarchyIntegrator*>::iterator it = d_child_integrators.begin(); it != d_child_integrators.end();
         ++it)
    {
        (*it)->setupPlotData();
    }
    return;
} // setupPlotData

int
HierarchyIntegrator::getNumberOfCycles() const
{
    return d_num_cycles;
} // getNumberOfCycles

int
HierarchyIntegrator::getCurrentCycleNumber() const
{
    return d_current_cycle_num;
} // getCurrentCycleNumber

double
HierarchyIntegrator::getCurrentTimeStepSize() const
{
    return d_current_dt;
} // getCurrentTimeStepSize

void
HierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                  const double new_time,
                                                  const int num_cycles)
{
    d_current_num_cycles = num_cycles;
    d_current_cycle_num = -1;
    d_current_dt = new_time - current_time;
#if !defined(NDEBUG)
    TBOX_ASSERT(d_current_num_cycles > 0);
    TBOX_ASSERT(d_current_dt > 0.0);
#endif
    return;
} // preprocessIntegrateHierarchy

void
HierarchyIntegrator::integrateHierarchy(const double current_time, const double new_time, const int cycle_num)
{
    ++d_current_cycle_num;
#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(d_current_dt, new_time - current_time));
    TBOX_ASSERT(d_current_cycle_num == cycle_num);
    TBOX_ASSERT(d_current_cycle_num < d_current_num_cycles);
#else
    NULL_USE(current_time);
    NULL_USE(new_time);
    NULL_USE(cycle_num);
#endif
    return;
} // integrateHierarchy

void
HierarchyIntegrator::skipCycle(const double current_time, const double new_time, const int cycle_num)
{
    ++d_current_cycle_num;
#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(d_current_dt, new_time - current_time));
    TBOX_ASSERT(d_current_cycle_num == cycle_num);
    TBOX_ASSERT(d_current_cycle_num < d_current_num_cycles);
#else
    NULL_USE(current_time);
    NULL_USE(new_time);
    NULL_USE(cycle_num);
#endif
    return;
} // skipCycle

void
HierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                   const double new_time,
                                                   const bool /*skip_synchronize_new_state_data*/,
                                                   const int num_cycles)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(d_current_dt, new_time - current_time));
    TBOX_ASSERT(num_cycles == d_current_num_cycles);
    TBOX_ASSERT(d_current_cycle_num + 1 == d_current_num_cycles);
#else
    NULL_USE(current_time);
    NULL_USE(new_time);
    NULL_USE(num_cycles);
#endif
    d_current_num_cycles = -1;
    d_current_cycle_num = -1;
    d_current_dt = std::numeric_limits<double>::quiet_NaN();
    return;
} // postprocessIntegrateHierarchy

void
HierarchyIntegrator::registerPreprocessIntegrateHierarchyCallback(PreprocessIntegrateHierarchyCallbackFcnPtr callback,
                                                                  void* ctx)
{
    d_preprocess_integrate_hierarchy_callbacks.push_back(callback);
    d_preprocess_integrate_hierarchy_callback_ctxs.push_back(ctx);
    return;
} // registerPreprocessIntegrateHierarchyCallback

void
HierarchyIntegrator::registerIntegrateHierarchyCallback(IntegrateHierarchyCallbackFcnPtr callback, void* ctx)
{
    d_integrate_hierarchy_callbacks.push_back(callback);
    d_integrate_hierarchy_callback_ctxs.push_back(ctx);
    return;
} // registerIntegrateHierarchyCallback

void
HierarchyIntegrator::registerPostprocessIntegrateHierarchyCallback(PostprocessIntegrateHierarchyCallbackFcnPtr callback,
                                                                   void* ctx)
{
    d_postprocess_integrate_hierarchy_callbacks.push_back(callback);
    d_postprocess_integrate_hierarchy_callback_ctxs.push_back(ctx);
    return;
} // registerPostprocessIntegrateHierarchyCallback

void
HierarchyIntegrator::registerApplyGradientDetectorCallback(ApplyGradientDetectorCallbackFcnPtr callback, void* ctx)
{
    d_apply_gradient_detector_callbacks.push_back(callback);
    d_apply_gradient_detector_callback_ctxs.push_back(ctx);
    return;
} // registerApplyGradientDetectorCallback

void
HierarchyIntegrator::initializeLevelData(const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
                                         const int level_number,
                                         const double init_data_time,
                                         const bool can_be_refined,
                                         const bool initial_time,
                                         const Pointer<BasePatchLevel<NDIM> > base_old_level,
                                         const bool allocate_data)
{
    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const Pointer<PatchLevel<NDIM> > old_level = base_old_level;
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (old_level)
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
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
    }
    else
    {
        level->setTime(init_data_time, d_current_data);
    }

    // Fill data from coarser levels in AMR hierarchy.
    if (!initial_time && (level_number > 0 || old_level))
    {
        level->allocatePatchData(d_scratch_data, init_data_time);
        std::vector<RefinePatchStrategy<NDIM>*> fill_after_regrid_prolong_patch_strategies;
        CartExtrapPhysBdryOp fill_after_regrid_extrap_bc_op(d_fill_after_regrid_bc_idxs, d_bdry_extrap_type);
        fill_after_regrid_prolong_patch_strategies.push_back(&fill_after_regrid_extrap_bc_op);
        if (d_fill_after_regrid_phys_bdry_bc_op)
        {
            fill_after_regrid_prolong_patch_strategies.push_back(d_fill_after_regrid_phys_bdry_bc_op);
        }
        RefinePatchStrategySet fill_after_regrid_patch_strategy_set(fill_after_regrid_prolong_patch_strategies.begin(),
                                                                    fill_after_regrid_prolong_patch_strategies.end(),
                                                                    false);
        d_fill_after_regrid_prolong_alg
            .createSchedule(level, old_level, level_number - 1, hierarchy, &fill_after_regrid_patch_strategy_set)
            ->fillData(init_data_time);
        level->deallocatePatchData(d_scratch_data);
    }

    // Initialize level data at the initial time.
    if (initial_time)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        for (std::list<Pointer<Variable<NDIM> > >::const_iterator cit = d_state_variables.begin();
             cit != d_state_variables.end();
             ++cit)
        {
            Pointer<Variable<NDIM> > var = *cit;
            const int var_current_idx = var_db->mapVariableAndContextToIndex(var, getCurrentContext());
            Pointer<CartGridFunction> var_init = d_state_var_init_fcns[var];
            if (var_init)
            {
                var_init->setDataOnPatchLevel(var_current_idx, var, level, init_data_time, initial_time);
            }
            else
            {
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double> > var_current_cc_data = patch->getPatchData(var_current_idx);
                    Pointer<EdgeData<NDIM, double> > var_current_ec_data = patch->getPatchData(var_current_idx);
                    Pointer<FaceData<NDIM, double> > var_current_fc_data = patch->getPatchData(var_current_idx);
                    Pointer<NodeData<NDIM, double> > var_current_nc_data = patch->getPatchData(var_current_idx);
                    Pointer<SideData<NDIM, double> > var_current_sc_data = patch->getPatchData(var_current_idx);
                    if (var_current_cc_data)
                        var_current_cc_data->fillAll(0.0);
                    else if (var_current_ec_data)
                        var_current_ec_data->fillAll(0.0);
                    else if (var_current_fc_data)
                        var_current_fc_data->fillAll(0.0);
                    else if (var_current_nc_data)
                        var_current_nc_data->fillAll(0.0);
                    else if (var_current_sc_data)
                        var_current_sc_data->fillAll(0.0);
                }
            }
        }
    }

    // Perform specialized data initialization.
    initializeLevelDataSpecialized(
        base_hierarchy, level_number, init_data_time, can_be_refined, initial_time, base_old_level, allocate_data);

    // Initialize data associated with any child integrators.
    for (std::set<HierarchyIntegrator*>::iterator it = d_child_integrators.begin(); it != d_child_integrators.end();
         ++it)
    {
        (*it)->initializeLevelData(
            base_hierarchy, level_number, init_data_time, can_be_refined, initial_time, base_old_level, allocate_data);
    }
    return;
} // initializeLevelData

void
HierarchyIntegrator::resetHierarchyConfiguration(const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
                                                 const int coarsest_level,
                                                 const int finest_level)
{
    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((coarsest_level >= 0) && (coarsest_level <= finest_level) &&
                (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(hierarchy->getPatchLevel(ln));
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Initialize or reset the hierarchy math operations object.
    d_hier_math_ops = buildHierarchyMathOps(hierarchy);
    if (d_manage_hier_math_ops)
    {
        d_hier_math_ops->setPatchHierarchy(hierarchy);
        d_hier_math_ops->resetLevels(0, finest_hier_level);
    }

    // If we have added or removed a level, resize the communication schedule
    // vectors.
    for (RefineAlgorithmMap::const_iterator it = d_ghostfill_algs.begin(); it != d_ghostfill_algs.end(); ++it)
    {
        d_ghostfill_scheds[it->first].resize(finest_hier_level + 1);
    }

    for (RefineAlgorithmMap::const_iterator it = d_prolong_algs.begin(); it != d_prolong_algs.end(); ++it)
    {
        d_prolong_scheds[it->first].resize(finest_hier_level + 1);
    }

    for (CoarsenAlgorithmMap::const_iterator it = d_coarsen_algs.begin(); it != d_coarsen_algs.end(); ++it)
    {
        d_coarsen_scheds[it->first].resize(finest_hier_level + 1);
    }

    // (Re)build ghost cell filling communication schedules.  These are created
    // for all levels in the hierarchy.
    for (RefineAlgorithmMap::const_iterator it = d_ghostfill_algs.begin(); it != d_ghostfill_algs.end(); ++it)
    {
        for (int ln = coarsest_level; ln <= std::min(finest_level + 1, finest_hier_level); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            d_ghostfill_scheds[it->first][ln] =
                it->second->createSchedule(level, ln - 1, hierarchy, d_ghostfill_strategies[it->first]);
        }
    }

    // (Re)build data prolongation communication schedules.  These are set only for levels
    // >= 1.
    for (RefineAlgorithmMap::const_iterator it = d_prolong_algs.begin(); it != d_prolong_algs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level, 1); ln <= std::min(finest_level + 1, finest_level); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            d_prolong_scheds[it->first][ln] = it->second->createSchedule(
                level, Pointer<PatchLevel<NDIM> >(), ln - 1, hierarchy, d_prolong_strategies[it->first]);
        }
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgorithmMap::const_iterator it = d_coarsen_algs.begin(); it != d_coarsen_algs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level, 1); ln <= std::min(finest_level + 1, finest_level); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(ln - 1);
            d_coarsen_scheds[it->first][ln] =
                it->second->createSchedule(coarser_level, level, d_coarsen_strategies[it->first]);
        }
    }

    // Perform specialized reset operations.
    resetHierarchyConfigurationSpecialized(base_hierarchy, coarsest_level, finest_level);

    // Reset data associated with any child integrators.
    for (std::set<HierarchyIntegrator*>::iterator it = d_child_integrators.begin(); it != d_child_integrators.end();
         ++it)
    {
        (*it)->resetHierarchyConfiguration(base_hierarchy, coarsest_level, finest_level);
    }
    return;
} // resetHierarchyConfiguration

void
HierarchyIntegrator::applyGradientDetector(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                           const int level_number,
                                           const double error_data_time,
                                           const int tag_index,
                                           const bool initial_time,
                                           const bool uses_richardson_extrapolation_too)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // First untag all cells.
    if (!d_parent_integrator)
    {
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, int> > tags_data = patch->getPatchData(tag_index);
            tags_data->fillAll(0);
        }
    }

    // Tag cells.
    applyGradientDetectorSpecialized(
        hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    executeApplyGradientDetectorCallbackFcns(
        hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);

    // ALlow child integrators to tag cells for refinement.
    for (std::set<HierarchyIntegrator*>::iterator it = d_child_integrators.begin(); it != d_child_integrators.end();
         ++it)
    {
        (*it)->applyGradientDetector(
            hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }
    return;
} // applyGradientDetector

Pointer<VariableContext>
HierarchyIntegrator::getCurrentContext() const
{
    return d_current_context;
} // getCurrentContext

Pointer<VariableContext>
HierarchyIntegrator::getNewContext() const
{
    return d_new_context;
} // getNewContext

Pointer<VariableContext>
HierarchyIntegrator::getScratchContext() const
{
    return d_scratch_context;
} // getScratchContext

bool
HierarchyIntegrator::isAllocatedPatchData(const int data_idx, int coarsest_ln, int finest_ln) const
{
    if (data_idx < 0) return false;
    if (coarsest_ln == -1) coarsest_ln = 0;
    if (finest_ln == -1) finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(data_idx)) return false;
    }
    return true;
} // isAllocatedPatchData

void
HierarchyIntegrator::allocatePatchData(const int data_idx, const double data_time, int coarsest_ln, int finest_ln) const
{
    if (data_idx < 0) return;
    if (coarsest_ln == -1) coarsest_ln = 0;
    if (finest_ln == -1) finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(data_idx)) level->allocatePatchData(data_idx, data_time);
    }
    return;
} // allocatePatchData

void
HierarchyIntegrator::deallocatePatchData(const int data_idx, int coarsest_ln, int finest_ln) const
{
    if (data_idx < 0) return;
    if (coarsest_ln == -1) coarsest_ln = 0;
    if (finest_ln == -1) finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(data_idx)) level->deallocatePatchData(data_idx);
    }
    return;
} // deallocatePatchData

Pointer<HierarchyMathOps>
HierarchyIntegrator::getHierarchyMathOps() const
{
    return d_hier_math_ops;
} // HierarchyMathOps

void
HierarchyIntegrator::putToDatabase(Pointer<Database> db)
{
    db->putInteger("HIERARCHY_INTEGRATOR_VERSION", HIERARCHY_INTEGRATOR_VERSION);
    db->putDouble("d_integrator_time", d_integrator_time);
    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    const int dt_previous_size = static_cast<int>(d_dt_previous.size());
    db->putInteger("d_dt_previous_size", dt_previous_size);
    if (dt_previous_size > 0)
    {
        const std::vector<double> dt_previous_vec(d_dt_previous.begin(), d_dt_previous.end());
        db->putDoubleArray("d_dt_previous_vec", &dt_previous_vec[0], dt_previous_size);
    }
    db->putDouble("d_dt_init", d_dt_init);
    db->putDouble("d_dt_min", d_dt_min);
    db->putDouble("d_dt_max", d_dt_max);
    db->putDouble("d_dt_growth_factor", d_dt_growth_factor);
    db->putInteger("d_integrator_step", d_integrator_step);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_num_cycles", d_num_cycles);
    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putString("d_regrid_mode", enum_to_string<RegridMode>(d_regrid_mode));
    db->putBool("d_enable_logging", d_enable_logging);
    db->putIntegerArray("d_tag_buffer", d_tag_buffer);
    db->putString("d_bdry_extrap_type", d_bdry_extrap_type);
    putToDatabaseSpecialized(db);
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

double
HierarchyIntegrator::getMinimumTimeStepSizeSpecialized()
{
    return d_dt_min;
} // getMinimumTimeStepSizeSpecialized

double
HierarchyIntegrator::getMaximumTimeStepSizeSpecialized()
{
    double dt = d_dt_max;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    if (initial_time)
    {
        dt = std::min(d_dt_init, d_dt_max);
    }
    else
    {
        dt = std::min(dt, std::max(1.0, d_dt_growth_factor) * d_dt_previous[0]);
    }
    return dt;
} // getMaximumTimeStepSizeSpecialized

void
HierarchyIntegrator::synchronizeHierarchyDataSpecialized(VariableContextType ctx_type)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const bool synch_current_data = ctx_type == CURRENT_DATA;
    const bool synch_new_data = ctx_type == NEW_DATA;
#if !defined(NDEBUG)
    TBOX_ASSERT(synch_current_data || synch_new_data);
#endif
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        if (synch_current_data) d_coarsen_scheds[SYNCH_CURRENT_DATA_ALG][ln]->coarsenData();
        if (synch_new_data) d_coarsen_scheds[SYNCH_NEW_DATA_ALG][ln]->coarsenData();
    }
    return;
} // synchronizeHierarchyDataSpecialized

void
HierarchyIntegrator::resetTimeDependentHierarchyDataSpecialized(const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Advance the simulation time.
    d_dt_previous.push_front(new_time - d_integrator_time);
    static const unsigned int MAX_DT_PREVIOUS_SIZE = 31;
    if (d_dt_previous.size() > MAX_DT_PREVIOUS_SIZE) d_dt_previous.pop_back();
    d_integrator_time = new_time;
    ++d_integrator_step;

    // Swap PatchData<NDIM> pointers between the current and new contexts.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    for (std::list<Pointer<Variable<NDIM> > >::const_iterator sv = d_state_variables.begin();
         sv != d_state_variables.end();
         ++sv)
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
#if !defined(NDEBUG)
                TBOX_ASSERT(src_data->getBox() == dst_data->getBox());
                TBOX_ASSERT(src_data->getGhostCellWidth() == dst_data->getGhostCellWidth());
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
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->setTime(d_integrator_time, d_current_data);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }
    return;
} // resetTimeDependentHierarchyDataSpecialized

void
HierarchyIntegrator::resetIntegratorToPreadvanceStateSpecialized()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Deallocate the scratch and new data and reset the time of the current
    // data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
        level->setTime(d_integrator_time, d_current_data);
    }
    return;
} // resetIntegratorToPreadvanceStateSpecialized

bool
HierarchyIntegrator::atRegridPointSpecialized() const
{
    if (d_parent_integrator)
    {
        return false;
    }
    else
    {
        return (d_integrator_step > 0) && (d_regrid_interval != 0) && (d_integrator_step % d_regrid_interval == 0);
    }
} // atRegridPointSpecialized

void
HierarchyIntegrator::setupPlotDataSpecialized()
{
    // intentionally blank
    return;
} // setupPlotDataSpecialized

void
HierarchyIntegrator::initializeLevelDataSpecialized(const Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
                                                    const int /*level_number*/,
                                                    const double /*init_data_time*/,
                                                    const bool /*can_be_refined*/,
                                                    const bool /*initial_time*/,
                                                    const Pointer<BasePatchLevel<NDIM> > /*old_level*/,
                                                    const bool /*allocate_data*/)
{
    // intentionally blank
    return;
} // initializeLevelDataSpecialized

void
HierarchyIntegrator::resetHierarchyConfigurationSpecialized(const Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
                                                            const int /*coarsest_level*/,
                                                            const int /*finest_level*/)
{
    // intentionally blank
    return;
} // resetHierarchyConfigurationSpecialized

void
HierarchyIntegrator::applyGradientDetectorSpecialized(const Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
                                                      const int /*level_number*/,
                                                      const double /*error_data_time*/,
                                                      const int /*tag_index*/,
                                                      const bool /*initial_time*/,
                                                      const bool /*uses_richardson_extrapolation_too*/)
{
    // intentionally blank
    return;
} // applyGradientDetectorSpecialized

void HierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> /*db*/)
{
    // intentionally blank
    return;
} // putToDatabaseSpecialized

void
HierarchyIntegrator::executePreprocessIntegrateHierarchyCallbackFcns(double current_time,
                                                                     double new_time,
                                                                     int num_cycles)
{
    std::vector<PreprocessIntegrateHierarchyCallbackFcnPtr>& callbacks = d_preprocess_integrate_hierarchy_callbacks;
    std::vector<void*>& ctxs = d_preprocess_integrate_hierarchy_callback_ctxs;
    for (unsigned int k = 0; k < callbacks.size(); ++k)
    {
        (*callbacks[k])(current_time, new_time, num_cycles, ctxs[k]);
    }
    return;
} // executePreprocessIntegrateHierarchyCallbackFcns

void
HierarchyIntegrator::executeIntegrateHierarchyCallbackFcns(double current_time, double new_time, int cycle_num)
{
    std::vector<IntegrateHierarchyCallbackFcnPtr>& callbacks = d_integrate_hierarchy_callbacks;
    std::vector<void*>& ctxs = d_integrate_hierarchy_callback_ctxs;
    for (unsigned int k = 0; k < callbacks.size(); ++k)
    {
        (*callbacks[k])(current_time, new_time, cycle_num, ctxs[k]);
    }
    return;
} // executeIntegrateHierarchyCallbackFcns

void
HierarchyIntegrator::executePostprocessIntegrateHierarchyCallbackFcns(double current_time,
                                                                      double new_time,
                                                                      bool skip_synchronize_new_state_data,
                                                                      int num_cycles)
{
    std::vector<PostprocessIntegrateHierarchyCallbackFcnPtr>& callbacks = d_postprocess_integrate_hierarchy_callbacks;
    std::vector<void*>& ctxs = d_postprocess_integrate_hierarchy_callback_ctxs;
    for (unsigned int k = 0; k < callbacks.size(); ++k)
    {
        (*callbacks[k])(current_time, new_time, skip_synchronize_new_state_data, num_cycles, ctxs[k]);
    }
    return;
} // executePostprocessIntegrateHierarchyCallbackFcns

void
HierarchyIntegrator::executeApplyGradientDetectorCallbackFcns(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                                              const int level_number,
                                                              const double error_data_time,
                                                              const int tag_index,
                                                              const bool initial_time,
                                                              const bool uses_richardson_extrapolation_too)
{
    std::vector<ApplyGradientDetectorCallbackFcnPtr>& callbacks = d_apply_gradient_detector_callbacks;
    std::vector<void*>& ctxs = d_apply_gradient_detector_callback_ctxs;
    for (unsigned int k = 0; k < callbacks.size(); ++k)
    {
        (*callbacks[k])(hierarchy,
                        level_number,
                        error_data_time,
                        tag_index,
                        initial_time,
                        uses_richardson_extrapolation_too,
                        ctxs[k]);
    }
    return;
} // executeApplyGradientDetectorCallbackFcns

void
HierarchyIntegrator::registerVariable(int& current_idx,
                                      int& new_idx,
                                      int& scratch_idx,
                                      const Pointer<Variable<NDIM> > variable,
                                      const IntVector<NDIM>& scratch_ghosts,
                                      const std::string& coarsen_name,
                                      const std::string& refine_name,
                                      Pointer<CartGridFunction> init_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(variable);
#endif
    d_state_var_init_fcns[variable] = init_fcn;

    const IntVector<NDIM> no_ghosts = 0;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    current_idx = -1; // insure that uninitialized variable patch data
    new_idx = -1;     // descriptor indices cause errors
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

    // Setup the scratch context.
    scratch_idx = var_db->registerVariableAndContext(variable, getScratchContext(), scratch_ghosts);
    d_scratch_data.setFlag(scratch_idx);

    // Get the data transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineOperator<NDIM> > refine_operator = grid_geom->lookupRefineOperator(variable, refine_name);
    Pointer<CoarsenOperator<NDIM> > coarsen_operator = grid_geom->lookupCoarsenOperator(variable, coarsen_name);

    // Setup the refine algorithm used to fill data in new or modified patch
    // levels following a regrid operation.
    if (refine_operator)
    {
        d_fill_after_regrid_bc_idxs.setFlag(scratch_idx);
        d_fill_after_regrid_prolong_alg.registerRefine(current_idx, current_idx, scratch_idx, refine_operator);
    }

    // Setup the SYNCH_CURRENT_DATA and SYNCH_NEW_DATA algorithms, used to
    // synchronize the data on the hierarchy.
    if (coarsen_operator)
    {
        d_coarsen_algs[SYNCH_CURRENT_DATA_ALG]->registerCoarsen(current_idx, current_idx, coarsen_operator);
        d_coarsen_algs[SYNCH_NEW_DATA_ALG]->registerCoarsen(new_idx, new_idx, coarsen_operator);
    }
    return;
} // registerVariable

void
HierarchyIntegrator::registerVariable(int& idx,
                                      const Pointer<Variable<NDIM> > variable,
                                      const IntVector<NDIM>& ghosts,
                                      Pointer<VariableContext> ctx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(variable);
#endif
    if (!ctx) ctx = getScratchContext();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    idx = -1; // insure that uninitialized variable patch data descriptor indices cause errors

    d_scratch_variables.push_back(variable);

    // Setup the scratch context.
    idx = var_db->registerVariableAndContext(variable, ctx, ghosts);
    if (*ctx == *getCurrentContext())
    {
        d_current_data.setFlag(idx);
        if (d_registered_for_restart)
        {
            var_db->registerPatchDataForRestart(idx);
        }
    }
    else if (*ctx == *getScratchContext())
        d_scratch_data.setFlag(idx);
    else if (*ctx == *getNewContext())
        d_new_data.setFlag(idx);
    else
    {
        TBOX_ERROR(d_object_name << "::registerVariable():\n"
                                 << "  unrecognized variable context: "
                                 << ctx->getName()
                                 << "\n"
                                 << "  variable context should be one of:\n"
                                 << "    "
                                 << getCurrentContext()->getName()
                                 << ", "
                                 << getNewContext()->getName()
                                 << ", or "
                                 << getScratchContext()->getName()
                                 << std::endl);
    }
    return;
} // registerVariable

void
HierarchyIntegrator::registerGhostfillRefineAlgorithm(const std::string& name,
                                                      Pointer<RefineAlgorithm<NDIM> > ghostfill_alg,
                                                      RefinePatchStrategy<NDIM>* ghostfill_patch_strategy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ghostfill_algs.find(name) == d_ghostfill_algs.end());
#endif
    d_ghostfill_algs[name] = ghostfill_alg;
    d_ghostfill_strategies[name] = ghostfill_patch_strategy;
} // registerGhostfillRefineAlgorithm

void
HierarchyIntegrator::registerProlongRefineAlgorithm(const std::string& name,
                                                    Pointer<RefineAlgorithm<NDIM> > prolong_alg,
                                                    RefinePatchStrategy<NDIM>* prolong_patch_strategy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_prolong_algs.find(name) == d_prolong_algs.end());
#endif
    d_prolong_algs[name] = prolong_alg;
    d_prolong_strategies[name] = prolong_patch_strategy;
} // registerProlongRefineAlgorithm

void
HierarchyIntegrator::registerCoarsenAlgorithm(const std::string& name,
                                              Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg,
                                              CoarsenPatchStrategy<NDIM>* coarsen_patch_strategy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_coarsen_algs.find(name) == d_coarsen_algs.end());
#endif
    d_coarsen_algs[name] = coarsen_alg;
    d_coarsen_strategies[name] = coarsen_patch_strategy;
} // registerCoarsenAlgorithm

Pointer<RefineAlgorithm<NDIM> >
HierarchyIntegrator::getGhostfillRefineAlgorithm(const std::string& name) const
{
    RefineAlgorithmMap::const_iterator alg_it = d_ghostfill_algs.find(name);
#if !defined(NDEBUG)
    TBOX_ASSERT(alg_it != d_ghostfill_algs.end());
#endif
    return alg_it->second;
} // getGhostfillRefineAlgorithm

Pointer<RefineAlgorithm<NDIM> >
HierarchyIntegrator::getProlongRefineAlgorithm(const std::string& name) const
{
    RefineAlgorithmMap::const_iterator alg_it = d_prolong_algs.find(name);
#if !defined(NDEBUG)
    TBOX_ASSERT(alg_it != d_prolong_algs.end());
#endif
    return alg_it->second;
} // getProlongRefineAlgorithm

Pointer<CoarsenAlgorithm<NDIM> >
HierarchyIntegrator::getCoarsenAlgorithm(const std::string& name) const
{
    CoarsenAlgorithmMap::const_iterator alg_it = d_coarsen_algs.find(name);
#if !defined(NDEBUG)
    TBOX_ASSERT(alg_it != d_coarsen_algs.end());
#endif
    return alg_it->second;
} // getCoarsenAlgorithm

const std::vector<Pointer<RefineSchedule<NDIM> > >&
HierarchyIntegrator::getGhostfillRefineSchedules(const std::string& name) const
{
    RefineScheduleMap::const_iterator sched_it = d_ghostfill_scheds.find(name);
#if !defined(NDEBUG)
    TBOX_ASSERT(sched_it != d_ghostfill_scheds.end());
#endif
    return sched_it->second;
} // getGhostfillRefineSchedules

const std::vector<Pointer<RefineSchedule<NDIM> > >&
HierarchyIntegrator::getProlongRefineSchedules(const std::string& name) const
{
    RefineScheduleMap::const_iterator sched_it = d_prolong_scheds.find(name);
#if !defined(NDEBUG)
    TBOX_ASSERT(sched_it != d_prolong_scheds.end());
#endif
    return sched_it->second;
} // getProlongRefineSchedules

const std::vector<Pointer<CoarsenSchedule<NDIM> > >&
HierarchyIntegrator::getCoarsenSchedules(const std::string& name) const
{
    CoarsenScheduleMap::const_iterator sched_it = d_coarsen_scheds.find(name);
#if !defined(NDEBUG)
    TBOX_ASSERT(sched_it != d_coarsen_scheds.end());
#endif
    return sched_it->second;
} // getCoarsenSchedules

void
HierarchyIntegrator::registerChildHierarchyIntegrator(HierarchyIntegrator* child_integrator)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(child_integrator != this);
#endif
    child_integrator->d_parent_integrator = this;
    d_child_integrators.insert(child_integrator);
    return;
} // registerChildHierarchyIntegrator

void
HierarchyIntegrator::registerParentHierarchyIntegrator(HierarchyIntegrator* parent_integrator)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(parent_integrator != this);
#endif
    d_parent_integrator = parent_integrator;
    parent_integrator->d_child_integrators.insert(this);
    d_manage_hier_math_ops = false;
    return;
} // registerParentHierarchyIntegrator

Pointer<HierarchyMathOps>
HierarchyIntegrator::buildHierarchyMathOps(Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    if (!d_parent_integrator)
    {
        if (!d_hier_math_ops)
        {
            d_hier_math_ops = new HierarchyMathOps(d_object_name + "::HierarchyMathOps", hierarchy);
        }
        d_manage_hier_math_ops = true;
    }
    else
    {
        d_hier_math_ops = d_parent_integrator->buildHierarchyMathOps(hierarchy);
        d_manage_hier_math_ops = false;
    }
    return d_hier_math_ops;
} // buildHierarchyMathOps

void
HierarchyIntegrator::setupTagBuffer(Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    const int finest_hier_ln = gridding_alg->getMaxLevels() - 1;
    Array<int> new_tag_buffer(std::max(finest_hier_ln, 1));
    new_tag_buffer[0] = 0;
    for (int i = 0; i < finest_hier_ln; ++i)
    {
        if (i < d_tag_buffer.size())
            new_tag_buffer[i] = d_tag_buffer[i];
        else if (i > 0)
            new_tag_buffer[i] = new_tag_buffer[i - 1];
    }
    d_tag_buffer = new_tag_buffer;
    return;
} // setupTagBuffer

/////////////////////////////// PRIVATE //////////////////////////////////////

void
HierarchyIntegrator::getFromInput(Pointer<Database> db, bool is_from_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    // Read in data members from input database.
    if (!is_from_restart && db->keyExists("start_time")) d_start_time = db->getDouble("start_time");
    if (db->keyExists("end_time")) d_end_time = db->getDouble("end_time");
    if (db->keyExists("dt_init")) d_dt_init = db->getDouble("dt_init");
    if (db->keyExists("dt_min")) d_dt_min = db->getDouble("dt_min");
    if (db->keyExists("dt_max")) d_dt_max = db->getDouble("dt_max");
    if (db->keyExists("grow_dt"))
        d_dt_growth_factor = db->getDouble("grow_dt");
    else if (db->keyExists("dt_growth_factor"))
        d_dt_growth_factor = db->getDouble("dt_growth_factor");
    if (db->keyExists("max_integrator_steps")) d_max_integrator_steps = db->getInteger("max_integrator_steps");
    if (db->keyExists("num_cycles")) d_num_cycles = db->getInteger("num_cycles");
    if (db->keyExists("regrid_interval")) d_regrid_interval = db->getInteger("regrid_interval");
    if (db->keyExists("regrid_mode")) d_regrid_mode = string_to_enum<RegridMode>(db->getString("regrid_mode"));
    if (db->keyExists("enable_logging")) d_enable_logging = db->getBool("enable_logging");
    if (db->keyExists("bdry_extrap_type")) d_bdry_extrap_type = db->getString("bdry_extrap_type");
    if (db->keyExists("tag_buffer")) d_tag_buffer = db->getIntegerArray("tag_buffer");
    return;
} // getFromInput

void
HierarchyIntegrator::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  restart database corresponding to " << d_object_name
                                 << " not found in restart file."
                                 << std::endl);
    }
    int ver = db->getInteger("HIERARCHY_INTEGRATOR_VERSION");
    if (ver != HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  restart file version different than class version." << std::endl);
    }
    d_integrator_time = db->getDouble("d_integrator_time");
    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    const int dt_previous_size = db->getInteger("d_dt_previous_size");
    if (dt_previous_size > 0)
    {
        std::vector<double> dt_previous_vec(dt_previous_size);
        db->getDoubleArray("d_dt_previous_vec", &dt_previous_vec[0], dt_previous_size);
        d_dt_previous = std::deque<double>(dt_previous_vec.begin(), dt_previous_vec.end());
    }
    else
    {
        d_dt_previous.clear();
    }
    d_dt_min = db->getDoubleWithDefault("d_dt_min", 0.0);
    d_dt_max = db->getDouble("d_dt_max");
    d_dt_growth_factor = db->getDouble("d_dt_growth_factor");
    d_integrator_step = db->getInteger("d_integrator_step");
    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");
    d_num_cycles = db->getInteger("d_num_cycles");
    d_regrid_interval = db->getInteger("d_regrid_interval");
    d_regrid_mode = string_to_enum<RegridMode>(db->getString("d_regrid_mode"));
    d_enable_logging = db->getBool("d_enable_logging");
    d_bdry_extrap_type = db->getString("d_bdry_extrap_type");
    d_tag_buffer = db->getIntegerArray("d_tag_buffer");
    return;
} // getFromRestart

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
