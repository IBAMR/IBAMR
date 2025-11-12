// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2025 by the IBAMR developers
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

#ifndef included_IBTK_CellNoCornersFillPattern
#define included_IBTK_CellNoCornersFillPattern

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "IntVector.h"
#include "VariableFillPattern.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoxGeometry;

template <int DIM>
class BoxOverlap;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CellCellNoCornersFillPattern is a concrete implementation of the
 * abstract base class SAMRAI::xfer::VariableFillPattern.  It is used to
 * calculate overlaps according to a pattern that limits overlaps to the
 * cell-centered ghost region surrounding a patch on the target level,
 * excluding all corners (and, in 3D, patch edges).
 *
 * On levels other than the target level (or in cases in which the target level
 * cannot be determined), the overlap pattern defaults to that provided by class
 * SAMRAI::pdat::CellOverlap.
 *
 */
class CellNoCornersFillPattern : public SAMRAI::xfer::VariableFillPattern<NDIM>
{
public:
    /*!
     * \brief Constructor.
     *
     * \param stencil_width        the width to fill
     * \param overwrite_interior   whether to include the patch interior
     *
     * \note The parameter overwrite_interior takes precedence over the value
     * passed in to the function calculateOverlap on the target patch level.
     */
    CellNoCornersFillPattern(int stencil_width, bool overwrite_interior);

    /*!
     * \brief Destructor.
     */
    ~CellNoCornersFillPattern() = default;

    /*!
     * Calculate overlaps between the destination and source geometries according
     * to the desired pattern.  This will return the portion of the intersection
     * of the geometries that lies in the ghost region of the specified width
     * surrounding the patch, excluding all edges and corners.  The patch is
     * identified by the argument dst_patch_box.
     *
     * \param dst_geometry        geometry object for destination box
     * \param src_geometry        geometry object for source box
     * \param dst_patch_box       box for the destination patch
     * \param src_mask            the source mask, the box resulting from shifting the source
     *box
     * \param overwrite_interior  controls whether or not to include the destination box
     *interior in
     *the overlap
     * \param src_offset          the offset between source and destination index space (src +
     *src_offset = dst)
     *
     * \return                    pointer to the calculated overlap object
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::BoxOverlap<NDIM> >
    calculateOverlap(const SAMRAI::hier::BoxGeometry<NDIM>& dst_geometry,
                     const SAMRAI::hier::BoxGeometry<NDIM>& src_geometry,
                     const SAMRAI::hier::Box<NDIM>& dst_patch_box,
                     const SAMRAI::hier::Box<NDIM>& src_mask,
                     bool overwrite_interior,
                     const SAMRAI::hier::IntVector<NDIM>& src_offset) const override;

    /*!
     * Calculate overlaps between the destination and source geometries according
     * to the desired pattern.  This will return the portion of the intersection
     * of the geometries that lies in the ghost region of the specified width
     * surrounding the patch, excluding all edges and corners.  The patch is
     * identified by the argument dst_patch_box.
     *
     * \param dst_geometry        geometry object for destination box
     * \param src_geometry        geometry object for source box
     * \param dst_patch_box       box for the destination patch
     * \param src_mask            the source mask, the box resulting from shifting the source
     *box
     * \param overwrite_interior  controls whether or not to include the destination box
     *interior in
     *the overlap
     * \param src_offset          the offset between source and destination index space (src +
     *src_offset = dst)
     * \param dst_level_num       the level of the patch hierarchy on which the dst boxes are
     *located
     * \param src_level_num       the level of the patch hierarchy on which the src boxes are
     *located
     *
     * \return                    pointer to the calculated overlap object
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::BoxOverlap<NDIM> >
    calculateOverlapOnLevel(const SAMRAI::hier::BoxGeometry<NDIM>& dst_geometry,
                            const SAMRAI::hier::BoxGeometry<NDIM>& src_geometry,
                            const SAMRAI::hier::Box<NDIM>& dst_patch_box,
                            const SAMRAI::hier::Box<NDIM>& src_mask,
                            bool overwrite_interior,
                            const SAMRAI::hier::IntVector<NDIM>& src_offset,
                            int dst_level_num,
                            int src_level_num) const override;

    /*!
     * Set the target patch level number for the variable fill pattern.
     */
    void setTargetPatchLevelNumber(int level_num) override;

    /*!
     * Returns the stencil width.
     */
    SAMRAI::hier::IntVector<NDIM>& getStencilWidth() override;

    /*!
     * Returns a string name identifier "CELL_NO_CORNERS_FILL_PATTERN".
     */
    const std::string& getPatternName() const override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CellNoCornersFillPattern() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CellNoCornersFillPattern(const CellNoCornersFillPattern& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CellNoCornersFillPattern& operator=(const CellNoCornersFillPattern& that) = delete;

    SAMRAI::hier::IntVector<NDIM> d_stencil_width;
    const bool d_overwrite_interior;
    int d_target_level_num = IBTK::invalid_level_number;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CellNoCornersFillPattern
