// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_SideSynchCopyFillPattern
#define included_IBTK_SideSynchCopyFillPattern

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "Box.h"
#include "IntVector.h"
#include "VariableFillPattern.h"
#include "tbox/Pointer.h"

#include <string>

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
 * \brief Class SideCellSynchCopyFillPattern is a concrete implementation of the
 * abstract base class SAMRAI::xfer::VariableFillPattern.  It is used to
 * calculate overlaps according to a pattern which limits overlaps to the
 * side-centered ghost region surrounding a patch appropriate for
 * "synchronizing" side-centered values at patch boundaries.
 */
class SideSynchCopyFillPattern : public SAMRAI::xfer::VariableFillPattern<NDIM>
{
public:
    /*!
     * \brief Default constructor
     */
    SideSynchCopyFillPattern() = default;

    /*!
     * \brief Destructor
     */
    ~SideSynchCopyFillPattern() = default;

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
     * Returns the stencil width.
     */
    SAMRAI::hier::IntVector<NDIM>& getStencilWidth() override;

    /*!
     * Returns a string name identifier "SIDE_SYNCH_COPY_FILL_PATTERN".
     */
    const std::string& getPatternName() const override;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SideSynchCopyFillPattern(const SideSynchCopyFillPattern& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SideSynchCopyFillPattern& operator=(const SideSynchCopyFillPattern& that) = delete;

    SAMRAI::hier::IntVector<NDIM> d_stencil_width = 1;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_SideSynchCopyFillPattern
