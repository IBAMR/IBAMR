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

#ifndef included_IBTK_CFOwnerFill
#define included_IBTK_CFOwnerFill

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <BoxList.h>
#include <CellData.h>
#include <IntVector.h>
#include <PatchData.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefinePatchStrategy.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <VariableFillPattern.h>

#include <array>
#include <map>
#include <vector>

namespace IBTK
{
/*!
 * \brief Class CFOwnerFill sets the coarse-fine interface ghost values of cell-centered or side-centered double
 * precision data on the levels of a patch hierarchy, one value for each ghost cell (or face), whatever the patch
 * layout. Along each direction (d, s), a ghost cell (face) c outside the fine level gets the quadratic interpolation of
 * the coarse-fine kernels if the first cell (face) c + k s e_d of the fine level is at the distance k, at most the
 * ratio to the coarser level and the fill width, and if c + (k + 1) s e_d is also in the fine level. The value of a
 * cell (face) with several directions is the average of the values; all other ghost cells are unchanged. The fill
 * width is the smaller of the ghost width of the data and the width of the fill pattern, and it must not exceed the
 * ratio to the coarser level.
 *
 * The patch that contains c + k s e_d computes the value along (d, s) from its own interior and from the coarse data
 * that its ordinary fill brings to postprocessRefine(). The values then travel through temporary levels of ghost
 * width zero: one level of owned cells (faces) for each direction and component, and one level of the cells (faces)
 * that get values, where they are averaged and from which they are copied into the ghost cells of every patch with
 * the fill pattern of the destination. The temporary levels do not change the maximum ghost width of the registered
 * patch data.
 *
 * Use postprocessRefine() as a refine hook of the schedule that fills the destination, then computeValues() and
 * copyValues() after the schedule.
 */
class CFOwnerFill : public SAMRAI::xfer::RefinePatchStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor. Sets up the levels coarsest_ln to finest_ln of the hierarchy (level 0 has no coarse-fine
     * interface).
     */
    CFOwnerFill(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                int coarsest_ln,
                int finest_ln,
                int dst_data_idx,
                int src_data_idx,
                SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern);

    ~CFOwnerFill();

    /*!
     * \brief Use other patch data indices and fill pattern. The centering, depth, and fill width of the data must not
     * change; a change stops with an error.
     */
    void resetPatchDataIndices(int dst_data_idx,
                               int src_data_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern);

    /*!
     * \return True if the level has coarse-fine ghost values to set.
     */
    bool hasValues(int ln) const;

    /*!
     * \brief Complete the values of level ln from the tangential values that postprocessRefine() computed, and
     * average them. Call it after the refine schedule of the level.
     */
    void computeValues(int ln);

    /*!
     * \brief Copy the values into the ghost cells of the destination data of level ln.
     */
    void copyValues(int ln, double fill_time);

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
    CFOwnerFill() = delete;
    CFOwnerFill(const CFOwnerFill& from) = delete;
    CFOwnerFill& operator=(const CFOwnerFill& that) = delete;

    // The value along (d, s) on the cells or faces index_box of component comp, which a patch computes.
    struct Task
    {
        int comp;
        int d;
        int s;
        SAMRAI::hier::Box<NDIM> index_box;
        SAMRAI::hier::Box<NDIM> kernel_box;
        int owned_level; // The level of owned values that receives the result, and its patch.
        int owned_patch;
        SAMRAI::pdat::ArrayData<NDIM, double>* owned_array; // The data of that patch.
    };

    // The tasks and the work data of a patch.
    struct PatchWork
    {
        SAMRAI::hier::Box<NDIM> box;
        std::vector<Task> tasks;
        std::array<SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM>>, 2 * NDIM> tangential;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM>> work, result;
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, int>> indicator;
    };

    // The cells (faces) of a component on a level of owned values.
    struct OwnedLevel
    {
        int code;
        int comp;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level;
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM>> alg;
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>> sched;
    };

    // A box of the level of values that get the average of the values of the directions in mask.
    struct Piece
    {
        SAMRAI::hier::Box<NDIM> box;
        int mask;
    };

    // The pieces and the data of a local patch of the level of values, for each component.
    struct ValuePatch
    {
        std::vector<std::vector<Piece>> pieces;
        std::vector<SAMRAI::pdat::ArrayData<NDIM, double>*> result;
        std::vector<std::array<SAMRAI::pdat::ArrayData<NDIM, double>*, 2 * NDIM>> slots;
    };

    // A task of any patch of the level, to find the values that go to a patch of values on the same process.
    struct TaskInfo
    {
        std::size_t ci;
        int code;
        int owned_level;
        int owned_patch;
        int rank;
        SAMRAI::hier::BoxList<NDIM> local_parts, wrapped_parts;
    };

    // A copy of values from a patch of owned values to a patch of values on the same process.
    struct DirectCopy
    {
        SAMRAI::pdat::ArrayData<NDIM, double>* dst;
        const SAMRAI::pdat::ArrayData<NDIM, double>* src;
        SAMRAI::hier::Box<NDIM> box;
    };

    struct Level
    {
        bool active = false;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level, value_level;
        // The cells (faces) of each component with a value along each (d, s), and all of them with periodic images.
        std::vector<std::array<SAMRAI::hier::BoxList<NDIM>, 2 * NDIM>> regions;
        std::vector<SAMRAI::hier::BoxList<NDIM>> all_regions;
        std::map<int, PatchWork> work;
        std::vector<OwnedLevel> owned;
        std::vector<ValuePatch> value_patches; // The local patches of the value level.
        std::vector<TaskInfo> task_info;
        std::vector<DirectCopy> direct;
        // The cells (faces) of each component whose values go to the level of values through a schedule.
        std::vector<std::array<SAMRAI::hier::BoxList<NDIM>, 2 * NDIM>> cross_regions;
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM>> copy_alg;
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>> copy_sched;
    };

    void setGhostWidth();
    void clearLevels();
    void checkFillWidth(int ln) const;
    void setUpLevel(int ln);
    void setUpSchedules(int ln);
    void computeGeometry(int ln, SAMRAI::hier::BoxList<NDIM>& value_boxes);
    void setUpValueLevel(int ln, const SAMRAI::hier::BoxList<NDIM>& value_boxes);
    void setUpDirect(int ln);
    void setUpPieces(int ln);
    void setUpWork(int ln);
    void computeTangential(int ln,
                           PatchWork& pw,
                           const SAMRAI::hier::Patch<NDIM>& coarse,
                           const SAMRAI::hier::Box<NDIM>& fine_box,
                           const SAMRAI::hier::IntVector<NDIM>& ratio);
    void computeNormal(PatchWork& pw, SAMRAI::hier::Patch<NDIM>& patch, const SAMRAI::hier::IntVector<NDIM>& ratio);

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    int d_coarsest_ln, d_finest_ln;
    int d_dst_idx, d_src_idx;
    int d_owned_idx = IBTK::invalid_index, d_result_idx = IBTK::invalid_index;
    std::array<int, 2 * NDIM> d_slot_idx;
    bool d_is_cell = true;
    SAMRAI::hier::IntVector<NDIM> d_ghost_width, d_fill_width;
    int d_work_gcw = 2;
    int d_depth = 1;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> d_fill_pattern;
    std::vector<int> d_comps; // -1 for cell data, else the axes.
    std::vector<Level> d_levels;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CFOwnerFill
