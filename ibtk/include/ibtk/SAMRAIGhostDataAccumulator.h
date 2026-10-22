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
 * \brief Class that sums data stored in the ghost regions of the patches of a patch hierarchy into the patches that own
 * them.
 *
 * Each level in the range given to the constructor is summed independently, including its periodic images. Ghost cells
 * that no patch of the level owns, such as those at a coarse-fine interface, are not summed. Afterwards the interior
 * copies of a location agree and the values in ghost regions are unspecified. The patches of a level must be at least
 * as wide as the ghost cell width in every direction in which the level has adjacent patches or is periodic; the
 * accumulation emits an error otherwise. The caller destroys the object before the levels are regridded or
 * redistributed. Cell-centered and side-centered data with double precision values are supported.
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
     * Accumulate data by summing values in ghost regions into the owning patch, without applying the transpose of a
     * physical boundary fill.
     *
     * \deprecated Use accumulateGhostData(int, RobinPhysBdryPatchStrategy*,
     * double, const std::vector<std::vector<int>>*) with a null boundary
     * operator.
     */
    void accumulateGhostData(const int idx);

    /*!
     * Accumulate the data with patch data index \a idx, including the transpose of the physical boundary fill of
     * \a bdry_op (skipped if null). All processes must call this function.
     *
     * \param idx  Patch data index of the data, accumulated in place, with the ghost cell width of the constructor.
     * \param bdry_op  Boundary operator whose transpose of the physical boundary fill the function applies, or null.
     * \param active_patch_nums  If not null, one vector per level, coarsest first, of the numbers of the patches owned
     * by this process whose data may be nonzero. The data of the other patches must be zero and the function applies no
     * boundary transpose to them, which is exact only if the transpose of their zero data is zero (a homogeneous
     * boundary operator).
     */
    void accumulateGhostData(int idx,
                             RobinPhysBdryPatchStrategy* bdry_op,
                             double fill_time,
                             const std::vector<std::vector<int>>* active_patch_nums = nullptr);

private:
    SAMRAIGhostDataAccumulator(const SAMRAIGhostDataAccumulator&) = delete;
    SAMRAIGhostDataAccumulator& operator=(const SAMRAIGhostDataAccumulator&) = delete;

    /*!
     * Data of one level, built when the level is first used.
     */
    struct LevelInfo;

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
     * Index of the patch data that holds the copy of the values that are
     * summed.
     */
    int d_scratch_idx = IBTK::invalid_index;

    /*!
     * Level information, indexed by level number.
     */
    std::vector<std::unique_ptr<LevelInfo>> d_level_info;
};
} // namespace IBTK

#endif
