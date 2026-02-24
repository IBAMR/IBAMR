// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/MergingLoadBalancer.h"
#include "ibtk/SecondaryHierarchy.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIArray.h"
#include "SAMRAIBasePatchHierarchy.h"
#include "SAMRAIBasePatchLevel.h"
#include "SAMRAIBergerRigoutsos.h"
#include "SAMRAIBox.h"
#include "SAMRAIBoxGeneratorStrategy.h"
#include "SAMRAICellData.h"
#include "SAMRAIDatabase.h"
#include "SAMRAIGriddingAlgorithm.h"
#include "SAMRAIHierarchyCellDataOpsInteger.h"
#include "SAMRAIHierarchyCellDataOpsReal.h"
#include "SAMRAIIntVector.h"
#include "SAMRAILoadBalancer.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRefineAlgorithm.h"
#include "SAMRAIRefineOperator.h"
#include "SAMRAIRefinePatchStrategy.h"
#include "SAMRAIStandardTagAndInitStrategy.h"
#include "SAMRAIStandardTagAndInitialize.h"
#include "SAMRAITimerManager.h"

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ibtk/app_namespaces.h"

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
// Timers.
static Timer* t_reinit;
static Timer* t_primary_to_secondary;
static Timer* t_secondary_to_primary;

/**
 * The secondary hierarchy is special in that we already know which cells it
 * should contain - i.e., we don't want to tag anything different than what is
 * already present. This class just tags cells that are already refined.
 */
class CopyRefinementTags : public SAMRAIStandardTagAndInitStrategy
{
    virtual void initializeLevelData(const SAMRAIPointer<SAMRAIBasePatchHierarchy> /*hierarchy*/,
                                     const int /*level_number*/,
                                     const double /*init_data_time*/,
                                     const bool /*can_be_refined*/,
                                     const bool /*initial_time*/,
                                     const SAMRAIPointer<SAMRAIBasePatchLevel> /*old_level*/ = nullptr,
                                     const bool /*allocate_data*/ = true) override
    {
    }

    virtual void resetHierarchyConfiguration(const SAMRAIPointer<SAMRAIBasePatchHierarchy> /*hierarchy*/,
                                             const int /*coarsest_level*/,
                                             const int /*finest_level*/) override
    {
    }

    virtual void applyGradientDetector(const SAMRAIPointer<SAMRAIBasePatchHierarchy> hierarchy,
                                       const int level_number,
                                       const double /*error_data_time*/,
                                       const int tag_index,
                                       const bool /*initial_time*/,
                                       const bool /*uses_richardson_extrapolation_too*/) override
    {
        if (level_number == hierarchy->getFinestLevelNumber()) return;
        SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = hierarchy;
        TBOX_ASSERT(patch_hierarchy);
        SAMRAIPointer<SAMRAIPatchLevel> level = hierarchy->getPatchLevel(level_number);

        if (!level->checkAllocated(tag_index)) level->allocatePatchData(tag_index, 0.0);
        SAMRAIHierarchyCellDataOpsInteger hier_cc_data_ops(hierarchy, level_number, level_number);
        // this is the only time we will tag things so set the tag values to zero
        hier_cc_data_ops.setToScalar(tag_index, 0);

        // OK - now we need to figure out which values need to be refined. Since
        // the hierarchy already satisfies a proper nesting condition we only
        // need to check the next finer level to figure out which cells should
        // be tagged on this level, since we don't care about tagging in the
        // secondary hierarchy - just load balancing (but samrai requires we do
        // both simultaneously)
        SAMRAIPointer<SAMRAIPatchLevel> finer_level = hierarchy->getPatchLevel(level_number + 1);
        const auto& boxes = finer_level->getBoxes();
        for (int i = 0; i < boxes.getNumberOfBoxes(); ++i)
        {
            const SAMRAIBox finer_box = boxes[i];
            auto coarsened_box = finer_box;
            coarsened_box.coarsen(finer_level->getRatioToCoarserLevel());
            {
                // Sanity check that we really do have proper nesting
                auto check_box = coarsened_box;
                check_box.refine(finer_level->getRatioToCoarserLevel());
                TBOX_ASSERT(check_box == finer_box);
            }

            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                const SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                const SAMRAIBox box = patch->getBox();
                SAMRAIPointer<SAMRAICellData<int>> tag_data = patch->getPatchData(tag_index);

                // SAMRAI's intersection code is wrong and returns nonsense
                // answers when the intersection should be empty
                if (box.intersects(coarsened_box))
                {
                    const auto intersection = box * coarsened_box;
                    tag_data->fill(1, intersection);
                }
            }
        }
    }
};
} // namespace

SecondaryHierarchy::SecondaryHierarchy(std::string name,
                                       SAMRAIPointer<SAMRAIDatabase> gridding_algorithm_db,
                                       SAMRAIPointer<SAMRAIDatabase> load_balancer_db)
    : d_object_name(name),
      d_tag_strategy(std::make_unique<CopyRefinementTags>()),
      d_eulerian_data_cache(std::make_shared<SAMRAIDataCache>())
{
    auto set_timer = [&](const char* name) { return SAMRAITimerManager::getManager()->getTimer(name); };

    IBTK_DO_ONCE(t_reinit = set_timer("IBTK::SecondaryHierarchy::reinit()");
                 t_primary_to_secondary = set_timer("IBTK::SecondaryHierarchy::transferPrimaryToSecondary()");
                 t_secondary_to_primary = set_timer("IBTK::SecondaryHierarchy::transferSecondaryToPrimary()");)

    // IBAMR consistently uses applyGradientDetector for everything so do that too
    SAMRAIPointer<InputDatabase> database(new InputDatabase(d_object_name + ":: tag_db"));
    database->putString("tagging_method", "GRADIENT_DETECTOR");
    d_error_detector = new SAMRAIStandardTagAndInitialize(d_object_name + "::tag", d_tag_strategy.get(), database);
    d_box_generator = new SAMRAIBergerRigoutsos();
    const std::string load_balancer_type = load_balancer_db->getStringWithDefault("type", "MERGING");
    if (load_balancer_type == "DEFAULT")
        d_load_balancer = new SAMRAILoadBalancer(load_balancer_db);
    else if (load_balancer_type == "MERGING")
        d_load_balancer = new MergingLoadBalancer(load_balancer_db);
    else
        TBOX_ERROR(d_object_name << "::SecondaryHierarchy():\n"
                                 << "unimplemented load balancer type " << load_balancer_type << std::endl);

    d_gridding_algorithm = new SAMRAIGriddingAlgorithm(d_object_name + "::gridding_alg",
                                                       gridding_algorithm_db,
                                                       d_error_detector,
                                                       d_box_generator,
                                                       d_load_balancer,
                                                       /*due to a bug in SAMRAI this *has* to be true*/ true);
}

void
SecondaryHierarchy::reinit(int coarsest_patch_level_number,
                           int finest_patch_level_number,
                           SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy)
{
    IBTK_TIMER_START(t_reinit);
    d_coarsest_patch_level_number = coarsest_patch_level_number;
    d_finest_patch_level_number = finest_patch_level_number;
    d_primary_hierarchy = patch_hierarchy;

    d_transfer_forward_schedules.clear();
    d_transfer_backward_schedules.clear();

    d_secondary_hierarchy = d_primary_hierarchy->makeRefinedPatchHierarchy(d_object_name + "::secondary_hierarchy",
                                                                           SAMRAIIntVector(1),
                                                                           /*register_for_restart*/ false);
    d_eulerian_data_cache->setPatchHierarchy(d_secondary_hierarchy);
    d_eulerian_data_cache->resetLevels(0, finest_patch_level_number);
    IBTK_TIMER_STOP(t_reinit);
}

void
SecondaryHierarchy::reinit(int coarsest_patch_level_number,
                           int finest_patch_level_number,
                           SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy,
                           int workload_idx)
{
    IBTK_TIMER_START(t_reinit);
    d_coarsest_patch_level_number = coarsest_patch_level_number;
    d_finest_patch_level_number = finest_patch_level_number;
    d_primary_hierarchy = patch_hierarchy;

    d_transfer_forward_schedules.clear();
    d_transfer_backward_schedules.clear();

    d_load_balancer->setWorkloadPatchDataIndex(workload_idx);
    auto new_secondary_hierarchy =
        d_primary_hierarchy->makeRefinedPatchHierarchy(d_object_name + "::secondary_hierarchy",
                                                       SAMRAIIntVector(1),
                                                       /*register_for_restart*/ false);

    // transfer workload data. Since we have cell-centered data we don't have to
    // worry about boundary stuff
    for (int ln = 0; ln <= d_finest_patch_level_number; ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> old_level = patch_hierarchy->getPatchLevel(ln);
        TBOX_ASSERT(old_level->checkAllocated(workload_idx));
        SAMRAIPointer<SAMRAIPatchLevel> new_level = new_secondary_hierarchy->getPatchLevel(ln);
        new_level->allocatePatchData(workload_idx);

        SAMRAIPointer<SAMRAIRefineAlgorithm> refine_algorithm = new SAMRAIRefineAlgorithm();
        SAMRAIPointer<SAMRAIRefineOperator> refine_op = nullptr;
        refine_algorithm->registerRefine(workload_idx, workload_idx, workload_idx, refine_op);
        auto schedule = refine_algorithm->createSchedule("DEFAULT_FILL", new_level, old_level);

        schedule->fillData(0.0);
    }

    // case where everything is on the coarsest level:
    if (d_finest_patch_level_number == 0)
    {
        d_gridding_algorithm->makeCoarsestLevel(new_secondary_hierarchy, 0.0 /*d_current_time*/);
    }
    else
    {
        // We don't need to buffer the tagging since we are just copying the refinement levels
        SAMRAIArray<int> tag_buffer(patch_hierarchy->getFinestLevelNumber() + 1);
        std::fill(tag_buffer.getPointer(), tag_buffer.getPointer() + tag_buffer.getSize(), 0);
        d_gridding_algorithm->regridAllFinerLevels(new_secondary_hierarchy,
                                                   std::max(d_coarsest_patch_level_number - 1, 0),
                                                   0.0 /*d_current_time*/,
                                                   tag_buffer);
    }

    d_secondary_hierarchy = new_secondary_hierarchy;
    d_eulerian_data_cache->setPatchHierarchy(new_secondary_hierarchy);
    d_eulerian_data_cache->resetLevels(0, d_finest_patch_level_number);

    // for analysis purposes, also transfer the workload data to the new partitioning
    for (int ln = 0; ln <= d_finest_patch_level_number; ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> old_level = patch_hierarchy->getPatchLevel(ln);
        TBOX_ASSERT(old_level->checkAllocated(workload_idx));
        SAMRAIPointer<SAMRAIPatchLevel> new_level = d_secondary_hierarchy->getPatchLevel(ln);
        if (!new_level->checkAllocated(workload_idx)) new_level->allocatePatchData(workload_idx);
        SAMRAIHierarchyCellDataOpsReal<double> hier_cc_data_ops(d_secondary_hierarchy, ln, ln);
        hier_cc_data_ops.setToScalar(workload_idx, 0.0, false);

        SAMRAIPointer<SAMRAIRefineAlgorithm> refine_algorithm = new SAMRAIRefineAlgorithm();
        SAMRAIPointer<SAMRAIRefineOperator> refine_op = nullptr;
        refine_algorithm->registerRefine(workload_idx, workload_idx, workload_idx, refine_op);
        auto schedule = refine_algorithm->createSchedule("DEFAULT_FILL", new_level, old_level);
        schedule->fillData(0.0);
    }
    IBTK_TIMER_STOP(t_reinit);
}

void
SecondaryHierarchy::transferPrimaryToSecondary(const int level_number,
                                               const int primary_data_idx,
                                               const int scratch_data_idx,
                                               const double data_time,
                                               SAMRAIRefinePatchStrategy* patch_strategy)
{
    TBOX_ASSERT(d_secondary_hierarchy);
    IBTK_TIMER_START(t_primary_to_secondary);
    const auto key = std::make_pair(level_number, std::make_pair(primary_data_idx, scratch_data_idx));
    if (d_transfer_forward_schedules.count(key) == 0)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = d_primary_hierarchy->getPatchLevel(level_number);
        SAMRAIPointer<SAMRAIPatchLevel> scratch_level = d_secondary_hierarchy->getPatchLevel(level_number);
        if (!scratch_level->checkAllocated(scratch_data_idx)) scratch_level->allocatePatchData(scratch_data_idx, 0.0);
        SAMRAIPointer<SAMRAIRefineAlgorithm> refine_algorithm = new SAMRAIRefineAlgorithm();
        SAMRAIPointer<SAMRAIRefineOperator> refine_op_f = nullptr;
        refine_algorithm->registerRefine(scratch_data_idx, primary_data_idx, scratch_data_idx, refine_op_f);
        d_transfer_forward_schedules[key] =
            refine_algorithm->createSchedule("DEFAULT_FILL", scratch_level, level, patch_strategy);
    }
    d_transfer_forward_schedules[key]->fillData(data_time);
    IBTK_TIMER_STOP(t_primary_to_secondary);
} // transferPrimaryToSecondary

void
SecondaryHierarchy::transferSecondaryToPrimary(const int level_number,
                                               const int primary_data_idx,
                                               const int scratch_data_idx,
                                               const double data_time,
                                               SAMRAIRefinePatchStrategy* patch_strategy)
{
    TBOX_ASSERT(d_secondary_hierarchy);
    IBTK_TIMER_START(t_secondary_to_primary);
    const auto key = std::make_pair(level_number, std::make_pair(primary_data_idx, scratch_data_idx));
    if (d_transfer_backward_schedules.count(key) == 0)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = d_primary_hierarchy->getPatchLevel(level_number);
        SAMRAIPointer<SAMRAIPatchLevel> scratch_level = d_secondary_hierarchy->getPatchLevel(level_number);
        SAMRAIPointer<SAMRAIRefineAlgorithm> refine_algorithm = new SAMRAIRefineAlgorithm();
        SAMRAIPointer<SAMRAIRefineOperator> refine_op_b = nullptr;
        refine_algorithm->registerRefine(primary_data_idx, scratch_data_idx, primary_data_idx, refine_op_b);
        d_transfer_backward_schedules[key] =
            refine_algorithm->createSchedule("DEFAULT_FILL", level, scratch_level, patch_strategy);
    }
    d_transfer_backward_schedules[key]->fillData(data_time);
    IBTK_TIMER_STOP(t_secondary_to_primary);
} // transferSecondaryToPrimary

std::shared_ptr<SAMRAIDataCache>
SecondaryHierarchy::getSAMRAIDataCache()
{
    return d_eulerian_data_cache;
}

SAMRAIPointer<SAMRAIPatchHierarchy>
SecondaryHierarchy::getPrimaryHierarchy()
{
    return d_primary_hierarchy;
}

SAMRAIPointer<SAMRAIPatchHierarchy>
SecondaryHierarchy::getSecondaryHierarchy()
{
    return d_secondary_hierarchy;
}

} // namespace IBTK
