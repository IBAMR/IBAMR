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
 * \brief Class that sums data stored in the ghost regions of the patches of a
 * patch hierarchy into the patches that own them.
 *
 * Each level in the range given to the constructor is treated on its own; no
 * data are moved between levels. On a level, the data of a patch p hold values
 * in their ghost data box, which is the ghost cell box for cell-centered data
 * and the side boxes of the ghost cell box for side-centered data. A location in
 * the ghost data box of p that lies in the interior data box of another patch q
 * of the same level is owned by q. The interior data box of q is the patch box
 * of q for cell-centered data and the side boxes of the patch box of q for
 * side-centered data, which include the sides on the boundary of q. For
 * side-centered data, the sides on the outermost layer of the ghost side box of
 * p can be on the boundary of a patch q that has no cell in common with the
 * ghost cell box of p; such a side is owned by q, too. Periodic
 * boundaries are taken into account: a location is identified with its periodic
 * images, and p and q may be the same patch if its ghost region overlaps its own
 * periodic image.
 *
 * Call the values that the patches of a level and their periodic images store
 * for one location, in patch interiors or in ghost regions, the copies of the
 * location. The accumulation replaces the value that each patch q stores for
 * each location in its interior data box by the sum of all copies of the
 * location. A side on the boundary between two patches is in the interior data
 * box of both: afterwards the interior copies of the two patches both hold the
 * sum of all copies, so all interior copies of a location agree. The values in
 * ghost regions are unspecified after an accumulation. Values at locations that
 * are in the interior data box of no patch of the level are not summed
 * anywhere: these are ghost cells at a coarse-fine interface that no patch of
 * the same level covers and locations outside the physical domain.
 *
 * The patches of a level must be at least as wide as the ghost cell width in
 * every direction in which the level has patches next to each other or is
 * periodic. SAMRAI gives periodic images only to patches that touch a periodic
 * boundary, and the transpose of the physical boundary fill acts only on patches
 * that touch the physical boundary, so a ghost region that extends past a narrow
 * patch to patches that do not touch the boundary would lose values. The
 * accumulation emits an error if a level does not satisfy this requirement.
 *
 * The communication is done with a SAMRAI schedule for each level, which does
 * not depend on the patch data index of the data that are accumulated. A
 * schedule is built the first time that its level is used and is kept until the
 * object is destroyed. The object holds a copy of the data that are being
 * accumulated. It is valid for one configuration of the patches of the levels: it
 * must be destroyed whenever the levels are regridded or their patches are
 * redistributed, and before the levels themselves are destroyed.
 *
 * Cell-centered and side-centered data with double precision values are
 * supported.
 */
class SAMRAIGhostDataAccumulator
{
public:
    /*!
     * Constructor.
     *
     * \param patch_hierarchy  Patch hierarchy that holds the data.
     * \param var              Variable of the data that are accumulated.
     * \param gcw              Ghost cell width of the data that are accumulated.
     * \param coarsest_ln      Coarsest level on which data are accumulated.
     * \param finest_ln        Finest level on which data are accumulated.
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
     * Accumulate data by summing values in ghost positions into the entry on
     * the owning patch, without applying the transpose of a physical boundary
     * fill.
     *
     * \deprecated Use accumulateGhostData(int, RobinPhysBdryPatchStrategy*,
     * double, const std::vector<std::vector<int>>*), which includes the
     * transpose of the physical boundary fill. With a null boundary operator it
     * does what this function does.
     */
    void accumulateGhostData(const int idx);

    /*!
     * Accumulate the data with patch data index \a idx, including the transpose
     * of the physical boundary fill. This function must be called by all
     * processes.
     *
     * On every patch of every level of this object the function first applies
     * RobinPhysBdryPatchStrategy::accumulateFromPhysicalBoundaryData() of
     * \a bdry_op, after setting the patch data index of \a bdry_op to \a idx and
     * with the ghost cell width of this object. This step is skipped if
     * \a bdry_op is null. It then sums the values in ghost cells or ghost sides
     * that are in the interior data box of another patch of the same level,
     * including periodic images, into that patch. After the call the values in
     * patch interiors, including the values on sides shared between patches,
     * hold the sums of the copies of the locations as described in the class
     * documentation, and all copies of a shared value agree. The values in ghost
     * regions are unspecified. No data move between levels.
     *
     * \param idx         Patch data index of the data. The data are accumulated
     *                    in place. The ghost cell width of the data must be the
     *                    one given to the constructor.
     * \param bdry_op     Boundary operator whose transpose of the physical
     *                    boundary fill is applied, or null.
     * \param fill_time   Time that is passed to the boundary operator.
     * \param active_patch_nums  If not null, one vector for each level of this
     *                    object, starting with the coarsest, that holds the
     *                    patch numbers (see SAMRAI::hier::Patch::getPatchNumber())
     *                    of the patches owned by this process whose data may be
     *                    nonzero. The data of all other patches must be zero. For
     *                    the other patches, the transpose of the physical boundary
     *                    fill is not applied, and the values in their ghost
     *                    regions are not sent to the patches that own them. A
     *                    level without any such patch on any process is skipped.
     *                    The result is the same as for a call without
     *                    \a active_patch_nums only if the transpose of the
     *                    physical boundary fill, applied to the data of an
     *                    inactive patch, would leave them zero. This is not
     *                    always the case: for side-centered data the transpose
     *                    of the fill of CartSideRobinPhysBdryOp assigns g/a to a
     *                    side that is normal to the physical boundary where the
     *                    condition is a Dirichlet condition with coefficients a
     *                    and g. It does leave the data zero if g is zero or the
     *                    boundary operator is homogeneous. The call also
     *                    performs one reduction over all processes.
     */
    void accumulateGhostData(int idx,
                             RobinPhysBdryPatchStrategy* bdry_op,
                             double fill_time,
                             const std::vector<std::vector<int>>* active_patch_nums = nullptr);

private:
    SAMRAIGhostDataAccumulator(const SAMRAIGhostDataAccumulator&) = delete;
    SAMRAIGhostDataAccumulator& operator=(const SAMRAIGhostDataAccumulator&) = delete;

    /*!
     * Data of one level: the level that holds the scratch patch data, and the
     * data that are built when the level is first used.
     */
    struct LevelInfo;

    /*!
     * Check the level and build its communication schedule, unless this has
     * been done.
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
