// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2026 by the IBAMR developers
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

#ifndef included_IBTK_SAMRAIGhostDataAccumulator
#define included_IBTK_SAMRAIGhostDataAccumulator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <BasePatchHierarchy.h>
#include <IntVector.h>
#include <Variable.h>

#include <memory>
#include <vector>

namespace IBTK
{
class RobinPhysBdryPatchStrategy;

/*!
 * \brief Class SAMRAIGhostDataAccumulator sums the ghost-region values of patches into the patches that own them.
 *
 * Each level in the range given to the constructor is summed independently, including its periodic images. Ghost cells
 * that no patch of the level owns, such as those at a coarse-fine interface, are not summed. Cell-centered and
 * side-centered double precision data are supported. Afterward the copies of a location in the patch interiors agree to
 * round-off. The values in ghost regions are unspecified.
 *
 * The accumulation stops with an error unless:
 * - patches are at least as wide as the ghost cell width in every direction in which the level has adjacent patches
 *   (within the ghost cell width of each other) or is periodic;
 * - each patch touches a physical boundary or is at least as far from it as the ghost cell width (for cell-centered
 *   data, also a periodic boundary);
 * - the levels are those that the constructor received.
 *
 * Destroy the object before regridding or redistributing the levels.
 */
class SAMRAIGhostDataAccumulator
{
public:
    /*!
     * Constructor.
     */
    SAMRAIGhostDataAccumulator(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> patch_hierarchy,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                               const SAMRAI::hier::IntVector<NDIM> gcw,
                               const int coarsest_ln,
                               const int finest_ln);

    /*!
     * Destructor.
     */
    ~SAMRAIGhostDataAccumulator();

    /*!
     * Accumulate data without applying the transpose of a physical boundary fill.
     *
     * @deprecated Use the overload with a boundary operator, and pass a null operator.
     */
    IBTK_DEPRECATED(
        "Use accumulateGhostData(int, RobinPhysBdryPatchStrategy*, double, const "
        "std::vector<std::vector<int>>*) with a null boundary operator instead.")
    void accumulateGhostData(const int idx);

    /*!
     * Accumulate the data with patch data index \a idx in place, including the transpose of the physical boundary fill
     * of \a bdry_op. All processes must call this function.
     *
     * \param idx  Patch data index of the data, with the ghost cell width of the constructor.
     * \param bdry_op  Boundary operator whose patch data index is set to \a idx, or null for no transpose.
     * \param fill_time  Time at which the boundary operator is evaluated.
     * \param active_patch_nums  If not null, one vector per level, coarsest first, of the numbers of the patches owned
     * by this process whose data may be nonzero. The function skips the other patches, whose data must be zero.
     */
    void accumulateGhostData(int idx,
                             RobinPhysBdryPatchStrategy* bdry_op,
                             double fill_time,
                             const std::vector<std::vector<int>>* active_patch_nums = nullptr);

private:
    SAMRAIGhostDataAccumulator(const SAMRAIGhostDataAccumulator&) = delete;
    SAMRAIGhostDataAccumulator& operator=(const SAMRAIGhostDataAccumulator&) = delete;

    /*!
     * The schedules associated with one level and their state, built when the level is first used.
     */
    struct LevelSchedules;

    /*!
     * Check the level and build its schedule, unless this has been done.
     */
    void initializeLevel(int ln);

    /*!
     * Pointer to the patch hierarchy under consideration.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> d_hierarchy;

    /*!
     * Pointer to the variable whose data layout we copied.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> d_var;

    /*!
     * Ghost cell width.
     */
    const SAMRAI::hier::IntVector<NDIM> d_gcw;

    /*!
     * Coarsest level of the patch hierarchy on which we work.
     */
    const int d_coarsest_ln = invalid_level_number;

    /*!
     * Finest level of the patch hierarchy on which we work.
     */
    const int d_finest_ln = invalid_level_number;

    /*!
     * Boolean indicating whether or not we have cell-centered data.
     */
    bool d_cc_data = true;

    /*!
     * Patch data index of the copy of the values that are summed.
     */
    int d_scratch_idx = IBTK::invalid_index;

    /*!
     * The schedules associated with each level, indexed by level number.
     */
    std::vector<std::unique_ptr<LevelSchedules>> d_level_schedules;
};
} // namespace IBTK

#endif
