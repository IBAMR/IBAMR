// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2021 by the IBAMR developers
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

#include "BergerRigoutsos.h"
#include "BoxGeneratorStrategy.h"
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "RefinePatchStrategy.h"
#include "StandardTagAndInitialize.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ibtk/app_namespaces.h"

namespace IBTK
{
SecondaryHierarchy::SecondaryHierarchy(std::string name,
                                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> gridding_algorithm_db,
                                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> load_balancer_db,
                                       SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>* tag_strategy)
    : d_object_name(name),
      d_coarsest_patch_level_number(-1),
      d_finest_patch_level_number(-1),
      d_eulerian_data_cache(std::make_shared<SAMRAIDataCache>())
{
    // IBFEMethod implements applyGradientDetector so we have to turn that on
    Pointer<InputDatabase> database(new InputDatabase(d_object_name + ":: tag_db"));
    database->putString("tagging_method", "GRADIENT_DETECTOR");
    d_error_detector = new StandardTagAndInitialize<NDIM>(d_object_name + "::tag", tag_strategy, database);
    d_box_generator = new BergerRigoutsos<NDIM>();
    const std::string load_balancer_type = load_balancer_db->getStringWithDefault("type", "MERGING");
    if (load_balancer_type == "DEFAULT")
        d_load_balancer = new LoadBalancer<NDIM>(load_balancer_db);
    else if (load_balancer_type == "MERGING")
        d_load_balancer = new MergingLoadBalancer(load_balancer_db);
    else
        TBOX_ERROR(d_object_name << "::SecondaryHierarchy():\n"
                                 << "unimplemented load balancer type " << load_balancer_type << std::endl);

    d_gridding_algorithm = new GriddingAlgorithm<NDIM>(d_object_name + "::gridding_alg",
                                                       gridding_algorithm_db,
                                                       d_error_detector,
                                                       d_box_generator,
                                                       d_load_balancer,
                                                       /*due to a bug in SAMRAI this *has* to be true*/ true);
}

void
SecondaryHierarchy::reinit(int coarsest_patch_level_number,
                           int finest_patch_level_number,
                           Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
{
    d_coarsest_patch_level_number = coarsest_patch_level_number;
    d_finest_patch_level_number = finest_patch_level_number;
    d_primary_hierarchy = patch_hierarchy;

    d_transfer_forward_schedules.clear();
    d_transfer_backward_schedules.clear();

    d_secondary_hierarchy = d_primary_hierarchy->makeRefinedPatchHierarchy(d_object_name + "::secondary_hierarchy",
                                                                           IntVector<NDIM>(1),
                                                                           /*register_for_restart*/ false);
    d_eulerian_data_cache->setPatchHierarchy(d_secondary_hierarchy);
    d_eulerian_data_cache->resetLevels(0, finest_patch_level_number);
}

void
SecondaryHierarchy::reinit(int coarsest_patch_level_number,
                           int finest_patch_level_number,
                           Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                           const tbox::Array<int>& tag_buffer)
{
    d_coarsest_patch_level_number = coarsest_patch_level_number;
    d_finest_patch_level_number = finest_patch_level_number;
    d_primary_hierarchy = patch_hierarchy;

    d_transfer_forward_schedules.clear();
    d_transfer_backward_schedules.clear();

    // case where everything is on the coarsest level:
    if (d_finest_patch_level_number == 0)
    {
        d_gridding_algorithm->makeCoarsestLevel(d_secondary_hierarchy, 0.0 /*d_current_time*/);
    }
    else
    {
        if (d_coarsest_patch_level_number == 0)
        {
            d_gridding_algorithm->makeCoarsestLevel(d_secondary_hierarchy, 0.0 /*d_current_time*/);
        }
        d_gridding_algorithm->regridAllFinerLevels(
            d_secondary_hierarchy, std::max(d_coarsest_patch_level_number - 1, 0), 0.0 /*d_current_time*/, tag_buffer);
    }

    d_eulerian_data_cache->setPatchHierarchy(d_secondary_hierarchy);
    d_eulerian_data_cache->resetLevels(0, finest_patch_level_number);
}

SAMRAI::xfer::RefineSchedule<NDIM>&
SecondaryHierarchy::getPrimaryToScratchSchedule(const int level_number,
                                                const int primary_data_idx,
                                                const int scratch_data_idx,
                                                SAMRAI::xfer::RefinePatchStrategy<NDIM>* patch_strategy)
{
    TBOX_ASSERT(d_secondary_hierarchy);
    const auto key = std::make_pair(level_number, std::make_pair(primary_data_idx, scratch_data_idx));
    if (d_transfer_forward_schedules.count(key) == 0)
    {
        Pointer<PatchLevel<NDIM> > level = d_primary_hierarchy->getPatchLevel(level_number);
        Pointer<PatchLevel<NDIM> > scratch_level = d_secondary_hierarchy->getPatchLevel(level_number);
        if (!scratch_level->checkAllocated(scratch_data_idx)) scratch_level->allocatePatchData(scratch_data_idx, 0.0);
        Pointer<RefineAlgorithm<NDIM> > refine_algorithm = new RefineAlgorithm<NDIM>();
        Pointer<RefineOperator<NDIM> > refine_op_f = nullptr;
        refine_algorithm->registerRefine(scratch_data_idx, primary_data_idx, scratch_data_idx, refine_op_f);
        d_transfer_forward_schedules[key] =
            refine_algorithm->createSchedule("DEFAULT_FILL", scratch_level, level, patch_strategy, false, nullptr);
    }
    return *d_transfer_forward_schedules[key];
} // getPrimaryToScratchSchedule

SAMRAI::xfer::RefineSchedule<NDIM>&
SecondaryHierarchy::getScratchToPrimarySchedule(const int level_number,
                                                const int primary_data_idx,
                                                const int scratch_data_idx,
                                                SAMRAI::xfer::RefinePatchStrategy<NDIM>* patch_strategy)
{
    TBOX_ASSERT(d_secondary_hierarchy);
    const auto key = std::make_pair(level_number, std::make_pair(primary_data_idx, scratch_data_idx));
    if (d_transfer_backward_schedules.count(key) == 0)
    {
        Pointer<PatchLevel<NDIM> > level = d_primary_hierarchy->getPatchLevel(level_number);
        Pointer<PatchLevel<NDIM> > scratch_level = d_secondary_hierarchy->getPatchLevel(level_number);
        Pointer<RefineAlgorithm<NDIM> > refine_algorithm = new RefineAlgorithm<NDIM>();
        Pointer<RefineOperator<NDIM> > refine_op_b = nullptr;
        refine_algorithm->registerRefine(primary_data_idx, scratch_data_idx, primary_data_idx, refine_op_b);
        d_transfer_backward_schedules[key] =
            refine_algorithm->createSchedule("DEFAULT_FILL", level, scratch_level, patch_strategy, false, nullptr);
    }
    return *d_transfer_backward_schedules[key];
} // getScratchToPrimarySchedule

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
SecondaryHierarchy::getGriddingAlgorithm()
{
    return d_gridding_algorithm;
}

std::shared_ptr<IBTK::SAMRAIDataCache>
SecondaryHierarchy::getSAMRAIDataCache()
{
    return d_eulerian_data_cache;
}

} // namespace IBTK
