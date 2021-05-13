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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_SecondaryHierarchy
#define included_IBTK_SecondaryHierarchy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/SAMRAIDataCache.h"

#include "BoxGeneratorStrategy.h"
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "RefinePatchStrategy.h"
#include "StandardTagAndInitStrategy.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace IBTK
{
/**
 * Class encapsulating the idea of a secondary hierarchy. In practice, this
 * is an equivalent (same index space) hierarchy to a given PatchHierarchy,
 * but is partitioned in parallel in a different way.
 */
struct SecondaryHierarchy
{
public:
    /**
     * Constructor - requires that the object be reinitialized. At the
     * present time this object still requires manual setup of the various
     * gridding classes.
     */
    SecondaryHierarchy(std::string name,
                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> gridding_algorithm_db,
                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> load_balancer_db,
                       SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>* tag_strategy);

    /**
     * Reinitialize the secondary hierarchy by performing a copy of the
     * input hierarchy. This is useful during staggered initialization when
     * we may not yet want to actually perform hierarchy-specific load
     * balancing.
     */
    void reinit(int coarsest_patch_level_number,
                int finest_patch_level_number,
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /**
     * Reinitialize the secondary hierarchy based on a new patch hierarchy.
     */
    void reinit(int coarsest_patch_level_number,
                int finest_patch_level_number,
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                const SAMRAI::tbox::Array<int>& tag_buffer,
                int workload_idx);

    /*!
     * Get the transfer schedule from the primary hierarchy to the scratch
     * hierarchy associated with the given level and index. If necessary the
     * schedule is created and stored in a map.
     *
     * If needed, a SAMRAI::xfer::RefinePatchStrategy object can be provided
     * for filling ghost data at physical boundaries.
     */
    SAMRAI::xfer::RefineSchedule<NDIM>&
    getPrimaryToScratchSchedule(int level_number,
                                int primary_data_idx,
                                int scratch_data_idx,
                                SAMRAI::xfer::RefinePatchStrategy<NDIM>* patch_strategy = nullptr);

    /*!
     * Get the transfer schedule from the scratch hierarchy to the primary
     * hierarchy associated with the given level and index. If necessary the
     * schedule is created and stored in a map.
     *
     * If needed, a SAMRAI::xfer::RefinePatchStrategy object can be provided
     * for filling ghost data at physical boundaries.
     */
    SAMRAI::xfer::RefineSchedule<NDIM>&
    getScratchToPrimarySchedule(int level_number,
                                int primary_data_idx,
                                int scratch_data_idx,
                                SAMRAI::xfer::RefinePatchStrategy<NDIM>* patch_strategy = nullptr);

    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > getGriddingAlgorithm();

    /*!
     * Get a copy of the pointer to the scratch patch object.
     */
    std::shared_ptr<IBTK::SAMRAIDataCache> getSAMRAIDataCache();

    /*!
     * Pointer to the primary patch hierarchy (i.e., the one not by this class).
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_primary_hierarchy;

    /*!
     * Pointer to the secondary patch hierarchy (i.e., the one managed by this class).
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_secondary_hierarchy;

protected:
    /**
     * Object name.
     */
    std::string d_object_name;

    /**
     * Coarsest level on which there are patches with elements (i.e.,
     * patches which should be considered by this object).
     */
    int d_coarsest_patch_level_number;

    /**
     * Finest level on which there are patches with elements (i.e.,
     * patches which should be considered by this object).
     */
    int d_finest_patch_level_number;

    /**
     * Error detector.
     *
     * @note this object has to be persistent since d_gridding_alg
     * requires it: see the note for that member object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::TagAndInitializeStrategy<NDIM> > d_error_detector;

    /**
     * Box generator.
     *
     * @note this object has to be persistent since d_gridding_alg
     * requires it: see the note for that member object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::BoxGeneratorStrategy<NDIM> > d_box_generator;

    /**
     * Load balancer.
     *
     * @note this object has to be persistent since d_scratch_gridding_alg
     * requires it: see the note for that member object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;

    /**
     * Gridding algorithm.
     *
     * @note this object has to be persistent because, due to a bug in SAMRAI,
     * it is impossible to create a SAMRAI::mesh::GriddingAlgorithm object in
     * a restarted simulation without a corresponding entry in the restart
     * database.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_algorithm;

    /*!
     * Refinement schedules for transferring data from the primary hierarchy to
     * the secondary (i.e., this) hierarchy. The key type is the level number
     * and a pair of indices (the primary and scratch, in that order).
     *
     * @note this function assumes that only data on the finest level needs to
     * be transferred.
     */
    std::map<std::pair<int, std::pair<int, int> >, SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >
        d_transfer_forward_schedules;

    /*!
     * Refinement schedules for transferring data from the secondary hierarchy
     * (i.e., the one managed by this object) to the primary hierarchy. The key
     * type is the level number and a pair of indices (the primary and scratch,
     * in that order).
     *
     * @note this function assumes that only data on the finest level needs to
     * be transferred.
     */
    std::map<std::pair<int, std::pair<int, int> >, SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >
        d_transfer_backward_schedules;

    /**
     * Scratch index data object for the patch hierarchy managed by this class.
     */
    std::shared_ptr<IBTK::SAMRAIDataCache> d_eulerian_data_cache;
};
} // namespace IBTK

#endif
