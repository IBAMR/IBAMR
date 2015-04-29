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

#include "SAMRAI/hier/PatchDataRestartManager.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/mesh/TagAndInitializeStrategy.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/RefinePatchStrategySet.h"
#include "ibtk/ibtk_enums.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"

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
                                         const boost::shared_ptr<Database>& input_db,
                                         bool register_for_restart)
    : d_fill_after_regrid_prolong_alg()
{
    TBOX_ASSERT(!object_name.empty());
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
    d_tag_buffer.resize(1);
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
    auto var_db = VariableDatabase::getDatabase();
    d_current_context = var_db->getContext(d_object_name + "::CURRENT");
    d_new_context = var_db->getContext(d_object_name + "::NEW");
    d_scratch_context = var_db->getContext(d_object_name + "::SCRATCH");

    // Create default communications algorithms.
    d_coarsen_algs[SYNCH_CURRENT_DATA_ALG] = boost::make_shared<CoarsenAlgorithm>(DIM);
    d_coarsen_algs[SYNCH_NEW_DATA_ALG] = boost::make_shared<CoarsenAlgorithm>(DIM);
    d_fill_after_regrid_phys_bdry_bc_op = NULL;
    return;
}

HierarchyIntegrator::~HierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }

    for (auto it = d_ghostfill_strategies.begin(); it != d_ghostfill_strategies.end(); ++it)
    {
        delete it->second;
    }

    for (auto it = d_prolong_strategies.begin(); it != d_prolong_strategies.end(); ++it)
    {
        delete it->second;
    }

    for (auto it = d_coarsen_strategies.begin(); it != d_coarsen_strategies.end(); ++it)
    {
        delete it->second;
    }
    return;
}

const std::string& HierarchyIntegrator::getName() const
{
    return d_object_name;
}

void HierarchyIntegrator::initializePatchHierarchy(const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                                   const boost::shared_ptr<GriddingAlgorithm>& gridding_alg)
{
    if (d_hierarchy_is_initialized || d_parent_integrator) return;

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT(gridding_alg);

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Initialize the hierarchy integrator and any child hierarchy integrators
    // prior to initializing the patch hierarchy.
    std::deque<HierarchyIntegrator*> hier_integrators(1, this);
    while (!hier_integrators.empty())
    {
        HierarchyIntegrator* integrator = hier_integrators.front();
        integrator->initializeHierarchyIntegrator(d_hierarchy, d_gridding_alg);
        integrator->setupTagBuffer();
        for (unsigned int i = 0; i < std::min(d_tag_buffer.size(), integrator->d_tag_buffer.size()); ++i)
        {
            d_tag_buffer[i] = std::max(d_tag_buffer[i], integrator->d_tag_buffer[i]);
        }
        hier_integrators.pop_front();
        hier_integrators.insert(hier_integrators.end(), integrator->d_child_integrators.begin(),
                                integrator->d_child_integrators.end());
    }
    plog << d_object_name << "::initializePatchHierarchy(): "
         << "tag_buffer =";
    for (unsigned int i = 0; i < d_tag_buffer.size(); ++i) plog << " " << d_tag_buffer[i];
    plog << "\n";

    // Initialize the patch hierarchy.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        d_gridding_alg->getTagAndInitializeStrategy()->resetHierarchyConfiguration(d_hierarchy, coarsest_ln, finest_ln);
    }
    else
    {
        const bool initial_time = true;
        d_gridding_alg->makeCoarsestLevel(d_start_time);
        int level_number = 0;
        bool done = false;
        while (!done && (d_hierarchy->levelCanBeRefined(level_number)))
        {
            d_gridding_alg->makeFinerLevel(d_tag_buffer[level_number], initial_time, d_integrator_step, d_integrator_time);
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
        hier_integrators.insert(hier_integrators.end(), integrator->d_child_integrators.begin(),
                                integrator->d_child_integrators.end());
    }
    return;
}

void HierarchyIntegrator::advanceHierarchy(double dt)
{
    const double dt_min = getMinimumTimeStepSize();
    const double dt_max = getMaximumTimeStepSize();
    if (dt < dt_min || dt > dt_max)
    {
        TBOX_ERROR(d_object_name << "::advanceHierarchy():\n"
                                 << "  at time = " << d_integrator_time << ": time step size dt = " << dt << "\n"
                                 << "  minimum time step size = " << dt_min << "\n"
                                 << "  maximum time step size = " << dt_max << "\n");
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
                                 << "  at time = " << d_integrator_time << ": time step size dt = " << dt << ".\n");
    }
    else if (dt == 0.0)
    {
        TBOX_ERROR(d_object_name << "::advanceHierarchy():\n"
                                 << "  at time = " << d_integrator_time << ": time step size dt = " << dt << ".\n");
    }
    else if (current_time == new_time || MathUtilities<double>::equalEps(current_time, new_time))
    {
        TBOX_ERROR(d_object_name << "::advanceHierarchy():\n"
                                 << "  at time = " << d_integrator_time << ": time step size dt = " << dt
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
        hier_integrators.insert(hier_integrators.end(), integrator->d_child_integrators.begin(),
                                integrator->d_child_integrators.end());
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
}

double HierarchyIntegrator::getMinimumTimeStepSize()
{
    double dt = getMinimumTimeStepSizeSpecialized();
    for (auto it = d_child_integrators.begin(); it != d_child_integrators.end(); ++it)
    {
        dt = std::max(dt, (*it)->getMinimumTimeStepSize());
    }
    return std::min(dt, d_end_time - d_integrator_time);
}

double HierarchyIntegrator::getMaximumTimeStepSize()
{
    double dt = getMaximumTimeStepSizeSpecialized();
    for (auto it = d_child_integrators.begin(); it != d_child_integrators.end(); ++it)
    {
        dt = std::min(dt, (*it)->getMaximumTimeStepSize());
    }
    return std::min(dt, d_end_time - d_integrator_time);
}

void HierarchyIntegrator::synchronizeHierarchyData(VariableContextType ctx_type)
{
    synchronizeHierarchyDataSpecialized(ctx_type);
    for (auto it = d_child_integrators.begin(); it != d_child_integrators.end(); ++it)
    {
        (*it)->synchronizeHierarchyData(ctx_type);
    }
    return;
}

void HierarchyIntegrator::resetTimeDependentHierarchyData(const double new_time)
{
    resetTimeDependentHierarchyDataSpecialized(new_time);
    for (auto it = d_child_integrators.begin(); it != d_child_integrators.end(); ++it)
    {
        (*it)->resetTimeDependentHierarchyData(new_time);
    }
    return;
}

void HierarchyIntegrator::resetIntegratorToPreadvanceState()
{
    resetIntegratorToPreadvanceStateSpecialized();
    for (auto it = d_child_integrators.begin(); it != d_child_integrators.end(); ++it)
    {
        (*it)->resetIntegratorToPreadvanceState();
    }
    return;
}

void HierarchyIntegrator::regridHierarchy()
{
    const int coarsest_ln = 0;

    // Regrid the hierarchy.
    switch (d_regrid_mode)
    {
    case STANDARD:
        d_gridding_alg->regridAllFinerLevels(coarsest_ln, d_tag_buffer, d_integrator_step, d_integrator_time);
        break;
    case AGGRESSIVE:
        for (int k = 0; k < d_hierarchy->getMaxNumberOfLevels(); ++k)
        {
            d_gridding_alg->regridAllFinerLevels(coarsest_ln, d_tag_buffer, d_integrator_step, d_integrator_time);
        }
        break;
    default:
        TBOX_ERROR(d_object_name << "::regridHierarchy():\n"
                                 << "  unrecognized regrid mode: " << enum_to_string<RegridMode>(d_regrid_mode) << "."
                                 << std::endl);
    }

    // Synchronize the state data on the patch hierarchy.
    synchronizeHierarchyData(CURRENT_DATA);
    return;
}

bool HierarchyIntegrator::atRegridPoint() const
{
    bool regrid_hierarchy = atRegridPointSpecialized();
    for (auto it = d_child_integrators.begin(); it != d_child_integrators.end() && !regrid_hierarchy; ++it)
    {
        regrid_hierarchy = regrid_hierarchy || (*it)->atRegridPoint();
    }
    return regrid_hierarchy;
}

double HierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}

double HierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}

double HierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}

int HierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}

int HierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}

bool HierarchyIntegrator::stepsRemaining() const
{
    return ((d_integrator_step < d_max_integrator_steps) && (d_integrator_time < d_end_time) &&
            !MathUtilities<double>::equalEps(d_integrator_time, d_end_time));
}

boost::shared_ptr<PatchHierarchy> HierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}

boost::shared_ptr<GriddingAlgorithm> HierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}

void HierarchyIntegrator::registerVisItDataWriter(const boost::shared_ptr<VisItDataWriter>& visit_writer)
{
    d_visit_writer = visit_writer;
    for (auto it = d_child_integrators.begin(); it != d_child_integrators.end(); ++it)
    {
        (*it)->registerVisItDataWriter(visit_writer);
    }
    return;
}

void HierarchyIntegrator::setupPlotData()
{
    setupPlotDataSpecialized();
    for (auto it = d_child_integrators.begin(); it != d_child_integrators.end(); ++it)
    {
        (*it)->setupPlotData();
    }
    return;
}

int HierarchyIntegrator::getNumberOfCycles() const
{
    return d_num_cycles;
}

int HierarchyIntegrator::getCurrentCycleNumber() const
{
    return d_current_cycle_num;
}

double HierarchyIntegrator::getCurrentTimeStepSize() const
{
    return d_current_dt;
}

void HierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                       const double new_time,
                                                       const int num_cycles)
{
    d_current_num_cycles = num_cycles;
    d_current_cycle_num = -1;
    d_current_dt = new_time - current_time;
    TBOX_ASSERT(d_current_num_cycles > 0);
    TBOX_ASSERT(d_current_dt > 0.0);
    return;
}

void HierarchyIntegrator::integrateHierarchy(const double current_time, const double new_time, const int cycle_num)
{
    ++d_current_cycle_num;
    TBOX_ASSERT(MathUtilities<double>::equalEps(d_current_dt, new_time - current_time));
    TBOX_ASSERT(d_current_cycle_num == cycle_num);
    TBOX_ASSERT(d_current_cycle_num < d_current_num_cycles);
    return;
}

void HierarchyIntegrator::skipCycle(const double current_time, const double new_time, const int cycle_num)
{
    ++d_current_cycle_num;
    TBOX_ASSERT(MathUtilities<double>::equalEps(d_current_dt, new_time - current_time));
    TBOX_ASSERT(d_current_cycle_num == cycle_num);
    TBOX_ASSERT(d_current_cycle_num < d_current_num_cycles);
    return;
}

void HierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                        const double new_time,
                                                        const bool /*skip_synchronize_new_state_data*/,
                                                        const int num_cycles)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(d_current_dt, new_time - current_time));
    TBOX_ASSERT(num_cycles == d_current_num_cycles);
    TBOX_ASSERT(d_current_cycle_num + 1 == d_current_num_cycles);
    d_current_num_cycles = -1;
    d_current_cycle_num = -1;
    d_current_dt = std::numeric_limits<double>::quiet_NaN();
    return;
}

void
HierarchyIntegrator::registerPreprocessIntegrateHierarchyCallback(PreprocessIntegrateHierarchyCallbackFcnPtr callback,
                                                                  void* ctx)
{
    d_preprocess_integrate_hierarchy_callbacks.push_back(callback);
    d_preprocess_integrate_hierarchy_callback_ctxs.push_back(ctx);
    return;
}

void HierarchyIntegrator::registerIntegrateHierarchyCallback(IntegrateHierarchyCallbackFcnPtr callback, void* ctx)
{
    d_integrate_hierarchy_callbacks.push_back(callback);
    d_integrate_hierarchy_callback_ctxs.push_back(ctx);
    return;
}

void
HierarchyIntegrator::registerPostprocessIntegrateHierarchyCallback(PostprocessIntegrateHierarchyCallbackFcnPtr callback,
                                                                   void* ctx)
{
    d_postprocess_integrate_hierarchy_callbacks.push_back(callback);
    d_postprocess_integrate_hierarchy_callback_ctxs.push_back(ctx);
    return;
}

void HierarchyIntegrator::registerApplyGradientDetectorCallback(ApplyGradientDetectorCallbackFcnPtr callback, void* ctx)
{
    d_apply_gradient_detector_callbacks.push_back(callback);
    d_apply_gradient_detector_callback_ctxs.push_back(ctx);
    return;
}

void HierarchyIntegrator::initializeLevelData(const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                              const int level_number,
                                              const double init_data_time,
                                              const bool can_be_refined,
                                              const bool initial_time,
                                              const boost::shared_ptr<PatchLevel>& old_level,
                                              const bool allocate_data)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (old_level)
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

    // Allocate storage needed to initialize the level and fill data from
    // coarser levels in AMR hierarchy, if any.
    //
    // Since time gets set when we allocate data, re-stamp it to current time if
    // we don't need to allocate.
    auto level = hierarchy->getPatchLevel(level_number);
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
        std::vector<boost::shared_ptr<RefinePatchStrategy> > fill_after_regrid_prolong_patch_strategies;
        auto fill_after_regrid_extrap_bc_op = boost::make_shared<CartExtrapPhysBdryOp>(d_fill_after_regrid_bc_idxs, d_bdry_extrap_type);
        fill_after_regrid_prolong_patch_strategies.push_back(fill_after_regrid_extrap_bc_op);
        if (d_fill_after_regrid_phys_bdry_bc_op)
        {
            fill_after_regrid_prolong_patch_strategies.push_back(d_fill_after_regrid_phys_bdry_bc_op);
        }
        RefinePatchStrategySet fill_after_regrid_patch_strategy_set(fill_after_regrid_prolong_patch_strategies.begin(),
                                                                    fill_after_regrid_prolong_patch_strategies.end());
        d_fill_after_regrid_prolong_alg.createSchedule(level, old_level, level_number - 1, hierarchy,
                                                       &fill_after_regrid_patch_strategy_set)->fillData(init_data_time);
        level->deallocatePatchData(d_scratch_data);
    }

    // Initialize level data at the initial time.
    if (initial_time)
    {
        auto var_db = VariableDatabase::getDatabase();
        for (auto cit = d_state_variables.begin(); cit != d_state_variables.end(); ++cit)
        {
            auto var = *cit;
            const int var_current_idx = var_db->mapVariableAndContextToIndex(var, getCurrentContext());
            auto var_init = d_state_var_init_fcns[var.get()];
            if (var_init)
            {
                var_init->setDataOnPatchLevel(var_current_idx, var, level, init_data_time, initial_time);
            }
            else
            {
                for (auto p = level->begin(); p != level->end(); ++p)
                {
                    auto patch = *p;
                    auto var_current_data = patch->getPatchData(var_current_idx);
                    auto var_current_cc_data = boost::dynamic_pointer_cast<CellData<double> >(var_current_data);
                    auto var_current_ec_data = boost::dynamic_pointer_cast<EdgeData<double> >(var_current_data);
                    auto var_current_fc_data = boost::dynamic_pointer_cast<FaceData<double> >(var_current_data);
                    auto var_current_nc_data = boost::dynamic_pointer_cast<NodeData<double> >(var_current_data);
                    auto var_current_sc_data = boost::dynamic_pointer_cast<SideData<double> >(var_current_data);
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
    initializeLevelDataSpecialized(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level,
                                   allocate_data);

    // Initialize data associated with any child integrators.
    for (auto it = d_child_integrators.begin(); it != d_child_integrators.end(); ++it)
    {
        (*it)->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level,
                                   allocate_data);
    }
    return;
}

void HierarchyIntegrator::resetHierarchyConfiguration(const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                                      const int coarsest_level,
                                                      const int finest_level)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((coarsest_level >= 0) && (coarsest_level <= finest_level) &&
                (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(hierarchy->getPatchLevel(ln));
    }
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
    for (auto it = d_ghostfill_algs.begin(); it != d_ghostfill_algs.end(); ++it)
    {
        d_ghostfill_scheds[it->first].resize(finest_hier_level + 1);
    }

    for (auto it = d_prolong_algs.begin(); it != d_prolong_algs.end(); ++it)
    {
        d_prolong_scheds[it->first].resize(finest_hier_level + 1);
    }

    for (auto it = d_coarsen_algs.begin(); it != d_coarsen_algs.end(); ++it)
    {
        d_coarsen_scheds[it->first].resize(finest_hier_level + 1);
    }

    // (Re)build ghost cell filling communication schedules.  These are created
    // for all levels in the hierarchy.
    for (auto it = d_ghostfill_algs.begin(); it != d_ghostfill_algs.end(); ++it)
    {
        for (int ln = coarsest_level; ln <= std::min(finest_level + 1, finest_hier_level); ++ln)
        {
            auto level = hierarchy->getPatchLevel(ln);
            d_ghostfill_scheds[it->first][ln] =
                it->second->createSchedule(level, ln - 1, hierarchy, d_ghostfill_strategies[it->first]);
        }
    }

    // (Re)build data prolongation communication schedules.  These are set only for levels
    // >= 1.
    for (auto it = d_prolong_algs.begin(); it != d_prolong_algs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level, 1); ln <= std::min(finest_level + 1, finest_level); ++ln)
        {
            auto level = hierarchy->getPatchLevel(ln);
            d_prolong_scheds[it->first][ln] =
                it->second->createSchedule(level, NULL, ln - 1, hierarchy, d_prolong_strategies[it->first]);
        }
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (auto it = d_coarsen_algs.begin(); it != d_coarsen_algs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level, 1); ln <= std::min(finest_level + 1, finest_level); ++ln)
        {
            auto level = hierarchy->getPatchLevel(ln);
            auto coarser_level = hierarchy->getPatchLevel(ln - 1);
            d_coarsen_scheds[it->first][ln] =
                it->second->createSchedule(coarser_level, level, d_coarsen_strategies[it->first]);
        }
    }

    // Perform specialized reset operations.
    resetHierarchyConfigurationSpecialized(hierarchy, coarsest_level, finest_level);

    // Reset data associated with any child integrators.
    for (auto it = d_child_integrators.begin(); it != d_child_integrators.end(); ++it)
    {
        (*it)->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    }
    return;
}

void HierarchyIntegrator::applyGradientDetector(const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                                const int level_number,
                                                const double error_data_time,
                                                const int tag_index,
                                                const bool initial_time,
                                                const bool uses_richardson_extrapolation_too)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    auto level = hierarchy->getPatchLevel(level_number);

    // First untag all cells.
    if (!d_parent_integrator)
    {
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch = *p;
            auto tags_data = BOOST_CAST<CellData<int> >(patch->getPatchData(tag_index));
            tags_data->fillAll(0);
        }
    }

    // Tag cells.
    applyGradientDetectorSpecialized(hierarchy, level_number, error_data_time, tag_index, initial_time,
                                     uses_richardson_extrapolation_too);
    executeApplyGradientDetectorCallbackFcns(hierarchy, level_number, error_data_time, tag_index, initial_time,
                                             uses_richardson_extrapolation_too);

    // ALlow child integrators to tag cells for refinement.
    for (auto it = d_child_integrators.begin(); it != d_child_integrators.end(); ++it)
    {
        (*it)->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time,
                                     uses_richardson_extrapolation_too);
    }
    return;
}

boost::shared_ptr<VariableContext> HierarchyIntegrator::getCurrentContext() const
{
    return d_current_context;
}

boost::shared_ptr<VariableContext> HierarchyIntegrator::getNewContext() const
{
    return d_new_context;
}

boost::shared_ptr<VariableContext> HierarchyIntegrator::getScratchContext() const
{
    return d_scratch_context;
}

bool HierarchyIntegrator::isAllocatedPatchData(const int data_idx, int coarsest_ln, int finest_ln) const
{
    if (data_idx < 0) return false;
    if (coarsest_ln == -1) coarsest_ln = 0;
    if (finest_ln == -1) finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(data_idx)) return false;
    }
    return true;
}

void
HierarchyIntegrator::allocatePatchData(const int data_idx, const double data_time, int coarsest_ln, int finest_ln) const
{
    if (data_idx < 0) return;
    if (coarsest_ln == -1) coarsest_ln = 0;
    if (finest_ln == -1) finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(data_idx)) level->allocatePatchData(data_idx, data_time);
    }
    return;
}

void HierarchyIntegrator::deallocatePatchData(const int data_idx, int coarsest_ln, int finest_ln) const
{
    if (data_idx < 0) return;
    if (coarsest_ln == -1) coarsest_ln = 0;
    if (finest_ln == -1) finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(data_idx)) level->deallocatePatchData(data_idx);
    }
    return;
}

boost::shared_ptr<HierarchyMathOps> HierarchyIntegrator::getHierarchyMathOps() const
{
    return d_hier_math_ops;
}

void HierarchyIntegrator::putToRestart(const boost::shared_ptr<Database>& db) const
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
        db->putDoubleVector("d_dt_previous_vec", dt_previous_vec);
    }
    db->putDouble("d_dt_min", d_dt_min);
    db->putDouble("d_dt_max", d_dt_max);
    db->putDouble("d_dt_growth_factor", d_dt_growth_factor);
    db->putInteger("d_integrator_step", d_integrator_step);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_num_cycles", d_num_cycles);
    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putString("d_regrid_mode", enum_to_string<RegridMode>(d_regrid_mode));
    db->putBool("d_enable_logging", d_enable_logging);
    db->putIntegerVector("d_tag_buffer", d_tag_buffer);
    db->putString("d_bdry_extrap_type", d_bdry_extrap_type);
    putToRestartSpecialized(db);
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

double HierarchyIntegrator::getMinimumTimeStepSizeSpecialized()
{
    return d_dt_min;
}

double HierarchyIntegrator::getMaximumTimeStepSizeSpecialized()
{
    double dt = d_dt_max;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    if (!initial_time && d_dt_growth_factor >= 1.0)
    {
        dt = std::min(dt, d_dt_growth_factor * d_dt_previous[0]);
    }
    return dt;
}

void HierarchyIntegrator::synchronizeHierarchyDataSpecialized(VariableContextType ctx_type)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const bool synch_current_data = ctx_type == CURRENT_DATA;
    const bool synch_new_data = ctx_type == NEW_DATA;
    TBOX_ASSERT(synch_current_data || synch_new_data);
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        if (synch_current_data) d_coarsen_scheds[SYNCH_CURRENT_DATA_ALG][ln]->coarsenData();
        if (synch_new_data) d_coarsen_scheds[SYNCH_NEW_DATA_ALG][ln]->coarsenData();
    }
    return;
}

void HierarchyIntegrator::resetTimeDependentHierarchyDataSpecialized(const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Advance the simulation time.
    d_dt_previous.push_front(new_time - d_integrator_time);
    static const unsigned int MAX_DT_PREVIOUS_SIZE = 31;
    if (d_dt_previous.size() > MAX_DT_PREVIOUS_SIZE) d_dt_previous.pop_back();
    d_integrator_time = new_time;
    ++d_integrator_step;

    // Swap PatchData pointers between the current and new contexts.
    auto var_db = VariableDatabase::getDatabase();
    for (auto sv = d_state_variables.begin(); sv != d_state_variables.end(); ++sv)
    {
        auto v = *sv;
        const int src_idx = var_db->mapVariableAndContextToIndex(v, getNewContext());
        const int dst_idx = var_db->mapVariableAndContextToIndex(v, getCurrentContext());
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            auto level = d_hierarchy->getPatchLevel(ln);
            for (auto p = level->begin(); p != level->end(); ++p)
            {
                auto patch = *p;
                auto src_data = patch->getPatchData(src_idx);
                auto dst_data = patch->getPatchData(dst_idx);
                TBOX_ASSERT(src_data->getBox().isSpatiallyEqual(dst_data->getBox()));
                TBOX_ASSERT(src_data->getGhostCellWidth() == dst_data->getGhostCellWidth());
                patch->setPatchData(dst_idx, src_data);
                patch->setPatchData(src_idx, dst_data);
            }
        }
    }

    // Deallocate the scratch and new data and reset the time of the current
    // data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        level->setTime(d_integrator_time, d_current_data);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }
    return;
}

void HierarchyIntegrator::resetIntegratorToPreadvanceStateSpecialized()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Deallocate the scratch and new data and reset the time of the current
    // data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
        level->setTime(d_integrator_time, d_current_data);
    }
    return;
}

bool HierarchyIntegrator::atRegridPointSpecialized() const
{
    if (d_parent_integrator)
    {
        return false;
    }
    else
    {
        return (d_integrator_step > 0) && (d_regrid_interval != 0) && (d_integrator_step % d_regrid_interval == 0);
    }
}

void HierarchyIntegrator::setupPlotDataSpecialized()
{
    // intentionally blank
    return;
}

void HierarchyIntegrator::initializeLevelDataSpecialized(const boost::shared_ptr<PatchHierarchy>& /*hierarchy*/,
                                                         const int /*level_number*/,
                                                         const double /*init_data_time*/,
                                                         const bool /*can_be_refined*/,
                                                         const bool /*initial_time*/,
                                                         const boost::shared_ptr<PatchLevel>& /*old_level*/,
                                                         const bool /*allocate_data*/)
{
    // intentionally blank
    return;
}

void HierarchyIntegrator::resetHierarchyConfigurationSpecialized(const boost::shared_ptr<PatchHierarchy>& /*hierarchy*/,
                                                                 const int /*coarsest_level*/,
                                                                 const int /*finest_level*/)
{
    // intentionally blank
    return;
}

void HierarchyIntegrator::applyGradientDetectorSpecialized(const boost::shared_ptr<PatchHierarchy>& /*hierarchy*/,
                                                           const int /*level_number*/,
                                                           const double /*error_data_time*/,
                                                           const int /*tag_index*/,
                                                           const bool /*initial_time*/,
                                                           const bool /*uses_richardson_extrapolation_too*/)
{
    // intentionally blank
    return;
}

void HierarchyIntegrator::putToRestartSpecialized(const boost::shared_ptr<Database>& /*db*/) const
{
    // intentionally blank
    return;
}

void HierarchyIntegrator::executePreprocessIntegrateHierarchyCallbackFcns(double current_time,
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
}

void HierarchyIntegrator::executeIntegrateHierarchyCallbackFcns(double current_time, double new_time, int cycle_num)
{
    std::vector<IntegrateHierarchyCallbackFcnPtr>& callbacks = d_integrate_hierarchy_callbacks;
    std::vector<void*>& ctxs = d_integrate_hierarchy_callback_ctxs;
    for (unsigned int k = 0; k < callbacks.size(); ++k)
    {
        (*callbacks[k])(current_time, new_time, cycle_num, ctxs[k]);
    }
    return;
}

void HierarchyIntegrator::executePostprocessIntegrateHierarchyCallbackFcns(double current_time,
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
}

void HierarchyIntegrator::executeApplyGradientDetectorCallbackFcns(const boost::shared_ptr<PatchHierarchy>& hierarchy,
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
        (*callbacks[k])(hierarchy, level_number, error_data_time, tag_index, initial_time,
                        uses_richardson_extrapolation_too, ctxs[k]);
    }
    return;
}

void HierarchyIntegrator::registerVariable(int& current_idx,
                                           int& new_idx,
                                           int& scratch_idx,
                                           const boost::shared_ptr<Variable>& variable,
                                           const IntVector& scratch_ghosts,
                                           const std::string& coarsen_name,
                                           const std::string& refine_name,
                                           const boost::shared_ptr<CartGridFunction>& init_fcn)
{
    TBOX_ASSERT(variable);
    d_state_var_init_fcns[variable.get()] = init_fcn;

    const IntVector no_ghosts = IntVector::getZero(DIM);

    auto var_db = VariableDatabase::getDatabase();

    current_idx = -1; // insure that uninitialized variable patch data
    new_idx = -1;     // descriptor indices cause errors
    scratch_idx = -1;

    d_state_variables.push_back(variable);

    // Setup the current context.
    current_idx = var_db->registerVariableAndContext(variable, getCurrentContext(), no_ghosts);
    d_current_data.setFlag(current_idx);
    if (d_registered_for_restart)
    {
        auto patch_data_restart_manager = PatchDataRestartManager::getManager();
        patch_data_restart_manager->registerPatchDataForRestart(current_idx);
    }

    // Setup the new context.
    new_idx = var_db->registerVariableAndContext(variable, getNewContext(), no_ghosts);
    d_new_data.setFlag(new_idx);

    // Setup the scratch context.
    scratch_idx = var_db->registerVariableAndContext(variable, getScratchContext(), scratch_ghosts);
    d_scratch_data.setFlag(scratch_idx);

    // Get the data transfer operators.
    auto grid_geom = BOOST_CAST<CartesianGridGeometry>(d_hierarchy->getGridGeometry());
    auto refine_operator = grid_geom->lookupRefineOperator(variable, refine_name);
    auto coarsen_operator = grid_geom->lookupCoarsenOperator(variable, coarsen_name);

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
}

void HierarchyIntegrator::registerVariable(int& idx,
                                           const boost::shared_ptr<Variable> &variable,
                                           const IntVector& ghosts,
                                           const boost::shared_ptr<VariableContext>& ctx)
{
    TBOX_ASSERT(variable);
    auto var_ctx = ctx ? ctx : getScratchContext();
    auto var_db = VariableDatabase::getDatabase();

    idx = -1; // insure that uninitialized variable patch data descriptor indices cause errors

    d_scratch_variables.push_back(variable);

    // Setup the scratch context.
    idx = var_db->registerVariableAndContext(variable, var_ctx, ghosts);
    if (*var_ctx == *getCurrentContext())
    {
        d_current_data.setFlag(idx);
        if (d_registered_for_restart)
        {
            auto patch_data_restart_manager = PatchDataRestartManager::getManager();
            patch_data_restart_manager->registerPatchDataForRestart(idx);
        }
    }
    else if (*var_ctx == *getScratchContext())
        d_scratch_data.setFlag(idx);
    else if (*var_ctx == *getNewContext())
        d_new_data.setFlag(idx);
    else
    {
        TBOX_ERROR(d_object_name << "::registerVariable():\n"
                                 << "  unrecognized variable context: " << var_ctx->getName() << "\n"
                                 << "  variable context should be one of:\n"
                                 << "    " << getCurrentContext()->getName() << ", " << getNewContext()->getName()
                                 << ", or " << getScratchContext()->getName() << std::endl);
    }
    return;
}

void HierarchyIntegrator::registerGhostfillRefineAlgorithm(const std::string& name,
                                                           const boost::shared_ptr<RefineAlgorithm>& ghostfill_alg,
                                                           RefinePatchStrategy* ghostfill_patch_strategy)
{
    TBOX_ASSERT(d_ghostfill_algs.find(name) == d_ghostfill_algs.end());
    d_ghostfill_algs[name] = ghostfill_alg;
    d_ghostfill_strategies[name] = ghostfill_patch_strategy;
}

void HierarchyIntegrator::registerProlongRefineAlgorithm(const std::string& name,
                                                         const boost::shared_ptr<RefineAlgorithm>& prolong_alg,
                                                         RefinePatchStrategy* prolong_patch_strategy)
{
    TBOX_ASSERT(d_prolong_algs.find(name) == d_prolong_algs.end());
    d_prolong_algs[name] = prolong_alg;
    d_prolong_strategies[name] = prolong_patch_strategy;
}

void HierarchyIntegrator::registerCoarsenAlgorithm(const std::string& name,
                                                   const boost::shared_ptr<CoarsenAlgorithm>& coarsen_alg,
                                                   CoarsenPatchStrategy* coarsen_patch_strategy)
{
    TBOX_ASSERT(d_coarsen_algs.find(name) == d_coarsen_algs.end());
    d_coarsen_algs[name] = coarsen_alg;
    d_coarsen_strategies[name] = coarsen_patch_strategy;
}

boost::shared_ptr<RefineAlgorithm> HierarchyIntegrator::getGhostfillRefineAlgorithm(const std::string& name) const
{
    RefineAlgorithmMap::const_iterator alg_it = d_ghostfill_algs.find(name);
    TBOX_ASSERT(alg_it != d_ghostfill_algs.end());
    return alg_it->second;
}

boost::shared_ptr<RefineAlgorithm> HierarchyIntegrator::getProlongRefineAlgorithm(const std::string& name) const
{
    RefineAlgorithmMap::const_iterator alg_it = d_prolong_algs.find(name);
    TBOX_ASSERT(alg_it != d_prolong_algs.end());
    return alg_it->second;
}

boost::shared_ptr<CoarsenAlgorithm> HierarchyIntegrator::getCoarsenAlgorithm(const std::string& name) const
{
    CoarsenAlgorithmMap::const_iterator alg_it = d_coarsen_algs.find(name);
    TBOX_ASSERT(alg_it != d_coarsen_algs.end());
    return alg_it->second;
}

const std::vector<boost::shared_ptr<RefineSchedule> >&
HierarchyIntegrator::getGhostfillRefineSchedules(const std::string& name) const
{
    RefineScheduleMap::const_iterator sched_it = d_ghostfill_scheds.find(name);
    TBOX_ASSERT(sched_it != d_ghostfill_scheds.end());
    return sched_it->second;
}

const std::vector<boost::shared_ptr<RefineSchedule> >&
HierarchyIntegrator::getProlongRefineSchedules(const std::string& name) const
{
    RefineScheduleMap::const_iterator sched_it = d_prolong_scheds.find(name);
    TBOX_ASSERT(sched_it != d_prolong_scheds.end());
    return sched_it->second;
}

const std::vector<boost::shared_ptr<CoarsenSchedule> >&
HierarchyIntegrator::getCoarsenSchedules(const std::string& name) const
{
    CoarsenScheduleMap::const_iterator sched_it = d_coarsen_scheds.find(name);
    TBOX_ASSERT(sched_it != d_coarsen_scheds.end());
    return sched_it->second;
}

void HierarchyIntegrator::registerChildHierarchyIntegrator(HierarchyIntegrator* child_integrator)
{
    TBOX_ASSERT(child_integrator != this);
    child_integrator->d_parent_integrator = this;
    d_child_integrators.insert(child_integrator);
    return;
}

void HierarchyIntegrator::registerParentHierarchyIntegrator(HierarchyIntegrator* parent_integrator)
{
    TBOX_ASSERT(parent_integrator != this);
    d_parent_integrator = parent_integrator;
    parent_integrator->d_child_integrators.insert(this);
    d_manage_hier_math_ops = false;
    return;
}

boost::shared_ptr<HierarchyMathOps>
HierarchyIntegrator::buildHierarchyMathOps(const boost::shared_ptr<PatchHierarchy>& hierarchy)
{
    if (!d_parent_integrator)
    {
        if (!d_hier_math_ops)
        {
            d_hier_math_ops = boost::make_shared<HierarchyMathOps>(d_object_name + "::HierarchyMathOps", hierarchy);
        }
        d_manage_hier_math_ops = true;
    }
    else
    {
        d_hier_math_ops = d_parent_integrator->buildHierarchyMathOps(hierarchy);
        d_manage_hier_math_ops = false;
    }
    return d_hier_math_ops;
}

void HierarchyIntegrator::setupTagBuffer()
{
    const int finest_hier_ln = d_hierarchy->getMaxNumberOfLevels() - 1;
    std::vector<int> new_tag_buffer(std::max(finest_hier_ln, 1));
    new_tag_buffer[0] = 0;
    for (unsigned int i = 0; static_cast<int>(i) < finest_hier_ln; ++i)
    {
        if (i < d_tag_buffer.size())
            new_tag_buffer[i] = d_tag_buffer[i];
        else if (i > 0)
            new_tag_buffer[i] = new_tag_buffer[i - 1];
    }
    d_tag_buffer = new_tag_buffer;
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void HierarchyIntegrator::getFromInput(const boost::shared_ptr<Database>& db, bool is_from_restart)
{
    TBOX_ASSERT(db);

    // Read in data members from input database.
    if (!is_from_restart && db->keyExists("start_time")) d_start_time = db->getDouble("start_time");
    if (db->keyExists("end_time")) d_end_time = db->getDouble("end_time");
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
    if (db->keyExists("tag_buffer")) d_tag_buffer = db->getIntegerVector("tag_buffer");
    return;
}

void HierarchyIntegrator::getFromRestart()
{
    auto restart_db = RestartManager::getManager()->getRootDatabase();
    boost::shared_ptr<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
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
        const std::vector<double>& dt_previous_vec = db->getDoubleVector("d_dt_previous_vec");
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
    d_tag_buffer = db->getIntegerVector("d_tag_buffer");
    return;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
