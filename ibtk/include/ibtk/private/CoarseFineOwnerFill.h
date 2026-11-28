// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_CoarseFineOwnerFill
#define included_IBTK_CoarseFineOwnerFill

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>
#include <tbox/Schedule.h>

#include <ArrayData.h>
#include <Box.h>
#include <CellData.h>
#include <IntVector.h>
#include <PatchData.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefinePatchStrategy.h>
#include <SideData.h>
#include <VariableFillPattern.h>

#include <array>
#include <map>
#include <memory>
#include <vector>

namespace IBTK
{
namespace detail
{
/*!
 * \brief Class CoarseFineOwnerFill sets the coarse-fine interface ghost values of cell-centered or side-centered double
 * precision data on the levels of a patch hierarchy, one value for each ghost cell (or face), whatever the patch
 * layout, the number of processes, and the periodic boundaries. Along each direction (d, s), a ghost cell (face) c
 * outside the fine level gets the quadratic interpolation of the coarse-fine kernels if the first cell (face)
 * c + k s e_d of the fine level is at the distance k, at most the ratio to the coarser level and the fill width, and
 * if c + (k + 1) s e_d is also in the fine level. The value of a cell (face) with several directions is the average of
 * the values; all other ghost cells are unchanged. The fill width is the smaller of the ghost width of the data and
 * the width of the fill pattern, and it must not exceed the ratio to the coarser level.
 *
 * Use postprocessRefine() as a refine hook of the schedule that fills the destination. After that schedule, call
 * computeValues(), run a schedule that has the transactions that appendTransactions() adds, and call copyValues().
 */
class CoarseFineOwnerFill : public SAMRAI::xfer::RefinePatchStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor. Sets up the levels coarsest_ln to finest_ln of the hierarchy (level 0 has no coarse-fine
     * interface).
     */
    CoarseFineOwnerFill(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                        int coarsest_ln,
                        int finest_ln,
                        int dst_data_idx,
                        int src_data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern);

    /*!
     * \brief Use other patch data indices and fill pattern. The centering, depth, and fill width of the data and the
     * shape of the fill pattern must not change; a change stops with an error.
     */
    void resetPatchDataIndices(int dst_data_idx,
                               int src_data_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern);

    /*!
     * \return True if this process has coarse-fine ghost values to compute or to set on level ln.
     */
    bool hasValues(int ln) const;

    /*!
     * \brief Complete the values of level ln that the local patches own from the tangential values that
     * postprocessRefine() computed. Call it after the refine schedule of the level.
     */
    void computeValues(int ln);

    /*!
     * \brief Append to a schedule the transactions that send the values of level ln to the patches that get them.
     * Every process must append the transactions of the same objects in the same order.
     */
    void appendTransactions(int ln, SAMRAI::tbox::Schedule& schedule);

    /*!
     * \brief Set the ghost cells of the destination data of level ln to the averages of the values that the schedule
     * brought.
     */
    void copyValues(int ln);

    /*!
     * \name Implementation of the RefinePatchStrategy interface.
     */
    //\{
    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch<NDIM>& patch,
                                       double fill_time,
                                       const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) override;

    SAMRAI::hier::IntVector<NDIM> getRefineOpStencilWidth() const override;

    void preprocessRefine(SAMRAI::hier::Patch<NDIM>& fine,
                          const SAMRAI::hier::Patch<NDIM>& coarse,
                          const SAMRAI::hier::Box<NDIM>& fine_box,
                          const SAMRAI::hier::IntVector<NDIM>& ratio) override;

    void postprocessRefine(SAMRAI::hier::Patch<NDIM>& fine,
                           const SAMRAI::hier::Patch<NDIM>& coarse,
                           const SAMRAI::hier::Box<NDIM>& fine_box,
                           const SAMRAI::hier::IntVector<NDIM>& ratio) override;
    //\}

private:
    CoarseFineOwnerFill() = delete;
    CoarseFineOwnerFill(const CoarseFineOwnerFill& from) = delete;
    CoarseFineOwnerFill& operator=(const CoarseFineOwnerFill& that) = delete;

    /*!
     * \brief The value along (d, s) on the cells or faces index_box of component comp, which a patch computes into
     * owned.
     */
    struct Task
    {
        //! The component: -1 for cell data, else the axis of the faces.
        int comp = -1;
        //! The direction (d, s): the axis d and the sign s of the step toward the fine level.
        int d = 0;
        int s = 0;
        //! The cells (faces) that get values and the box that the kernels fill.
        SAMRAI::hier::Box<NDIM> index_box;
        SAMRAI::hier::Box<NDIM> kernel_box;
        //! The values that the kernels compute.
        std::unique_ptr<SAMRAI::pdat::ArrayData<NDIM, double>> owned;
    };

    /*!
     * \brief The tasks and the work data of a local patch.
     */
    struct PatchWork
    {
        //! The box of the patch.
        SAMRAI::hier::Box<NDIM> box;
        //! The tasks of the patch.
        std::vector<Task> tasks;
        //! The tangential values of each direction, indexed by the location index.
        std::array<SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM>>, 2 * NDIM> tangential;
        //! The data that the normal kernels read and write.
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM>> work, result;
        //! The side-centered indicator that the kernels use.
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, int>> indicator;
    };

    /*!
     * \brief A box of ghost cells (faces) that gets the average of the values of the directions in directions.
     */
    struct Piece
    {
        //! The cells (faces) of the piece.
        SAMRAI::hier::Box<NDIM> box;
        //! The directions that contribute, indexed by the location index.
        std::array<bool, 2 * NDIM> directions = {};
    };

    /*!
     * \brief The values that a local patch gets: for each component, the values of each direction and the pieces.
     */
    struct Receiver
    {
        //! The values of each direction, indexed by the location index.
        std::vector<std::array<std::unique_ptr<SAMRAI::pdat::ArrayData<NDIM, double>>, 2 * NDIM>> slots;
        //! The pieces to average.
        std::vector<std::vector<Piece>> pieces;
    };

    /*!
     * \brief The values of task task of the patch src_patch that go to the ghost cells (faces) box of the patch
     * dst_patch, shifted by offset.
     */
    struct Transfer
    {
        //! The patch numbers of the destination and the source, and the index of the task in the source.
        int dst_patch = -1;
        int src_patch = -1;
        int task = -1;
        //! The index of the component and the location index of the direction.
        std::size_t ci = 0;
        int code = 0;
        //! The periodic shift of the source.
        SAMRAI::hier::IntVector<NDIM> offset = SAMRAI::hier::IntVector<NDIM>(0);
        //! The cells (faces) of the destination.
        SAMRAI::hier::Box<NDIM> box;
    };

    /*!
     * \brief The data of one level.
     */
    struct Level
    {
        //! Whether this process has values to compute or to set on the level.
        bool active = false;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level;
        //! The work of the local patches, by patch number.
        std::map<int, PatchWork> work;
        //! The values that the local patches get, by patch number.
        std::map<int, Receiver> receivers;
        //! The transfers that involve a local patch, in the same order on every process.
        std::vector<Transfer> transfers;
    };

    /*!
     * \brief Set the ghost width, the fill width, and the width of the work data from the destination data and the
     * fill pattern.
     */
    void setGhostWidth();

    /*!
     * \brief Stop with an error if the fill width exceeds the ratio of level ln to the coarser level.
     */
    void checkFillWidth(int ln) const;

    /*!
     * \brief Find the tasks and the transfers of level ln.
     */
    void setUpLevel(int ln);

    /*!
     * \brief Allocate the values and find the pieces that the local patches of level ln get.
     */
    void setUpReceivers(int ln);

    /*!
     * \brief Allocate the work data of the local patches of level ln.
     */
    void setUpWork(int ln);

    /*!
     * \brief Compute the tangential values of a patch from the coarse data of the refine schedule.
     */
    void computeTangential(int ln,
                           PatchWork& pw,
                           const SAMRAI::hier::Patch<NDIM>& coarse,
                           const SAMRAI::hier::Box<NDIM>& fine_box,
                           const SAMRAI::hier::IntVector<NDIM>& ratio);

    /*!
     * \brief Compute the values of a patch from its tangential values and its fine data.
     */
    void computeNormal(PatchWork& pw, SAMRAI::hier::Patch<NDIM>& patch, const SAMRAI::hier::IntVector<NDIM>& ratio);

    //! The patch hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;

    //! The levels that this object sets up.
    int d_coarsest_ln = 0;
    int d_finest_ln = 0;

    //! The patch data indices of the destination and of the source.
    int d_dst_idx = invalid_index;
    int d_src_idx = invalid_index;

    //! Whether the data are cell-centered (otherwise side-centered).
    bool d_is_cell = true;

    //! The ghost width of the destination data and the width of the fill.
    SAMRAI::hier::IntVector<NDIM> d_ghost_width = SAMRAI::hier::IntVector<NDIM>(0);
    SAMRAI::hier::IntVector<NDIM> d_fill_width = SAMRAI::hier::IntVector<NDIM>(0);

    //! The ghost width of the work data.
    int d_work_gcw = 2;

    //! The depth of the data.
    int d_depth = 1;

    //! The fill pattern of the destination.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> d_fill_pattern;

    //! The components: -1 for cell data, else the axes.
    std::vector<int> d_comps;

    //! The data of each level.
    std::vector<Level> d_levels;
};
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CoarseFineOwnerFill
