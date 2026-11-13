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

#ifndef included_IBTK_RegionFillPattern
#define included_IBTK_RegionFillPattern

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <BoxList.h>
#include <IntVector.h>
#include <VariableFillPattern.h>

#include <string>
#include <vector>

namespace IBTK
{
/*!
 * \brief Class RegionFillPattern limits the destination boxes of another fill pattern to given regions of cells (cell
 * data) or faces (side data).
 *
 * For side data, the regions are one list of face index boxes for each axis. For cell data, the first list contains
 * the cells.
 */
class RegionFillPattern : public SAMRAI::xfer::VariableFillPattern<NDIM>
{
public:
    /*!
     * \brief Constructor.
     */
    RegionFillPattern(SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> inner,
                      std::vector<SAMRAI::hier::BoxList<NDIM>> regions);

    /*!
     * \name Implementation of the VariableFillPattern interface.
     */
    //\{
    SAMRAI::tbox::Pointer<SAMRAI::hier::BoxOverlap<NDIM>>
    calculateOverlap(const SAMRAI::hier::BoxGeometry<NDIM>& dst_geometry,
                     const SAMRAI::hier::BoxGeometry<NDIM>& src_geometry,
                     const SAMRAI::hier::Box<NDIM>& dst_patch_box,
                     const SAMRAI::hier::Box<NDIM>& src_mask,
                     bool overwrite_interior,
                     const SAMRAI::hier::IntVector<NDIM>& src_offset) const override;

    SAMRAI::tbox::Pointer<SAMRAI::hier::BoxOverlap<NDIM>>
    calculateOverlapOnLevel(const SAMRAI::hier::BoxGeometry<NDIM>& dst_geometry,
                            const SAMRAI::hier::BoxGeometry<NDIM>& src_geometry,
                            const SAMRAI::hier::Box<NDIM>& dst_patch_box,
                            const SAMRAI::hier::Box<NDIM>& src_mask,
                            bool overwrite_interior,
                            const SAMRAI::hier::IntVector<NDIM>& src_offset,
                            int dst_level_num,
                            int src_level_num) const override;

    void setTargetPatchLevelNumber(int level_num) override;

    SAMRAI::hier::IntVector<NDIM>& getStencilWidth() override;

    const std::string& getPatternName() const override;
    //\}

private:
    SAMRAI::tbox::Pointer<SAMRAI::hier::BoxOverlap<NDIM>>
    restrictOverlap(SAMRAI::tbox::Pointer<SAMRAI::hier::BoxOverlap<NDIM>> overlap,
                    const SAMRAI::hier::IntVector<NDIM>& src_offset) const;

    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> d_inner;
    std::vector<SAMRAI::hier::BoxList<NDIM>> d_regions;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_RegionFillPattern
