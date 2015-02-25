// Filename: SideNoCornersFillPattern.h
// Created on 09 Mar 2010 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_SideNoCornersFillPattern
#define included_SideNoCornersFillPattern

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/xfer/VariableFillPattern.h"

namespace SAMRAI
{
namespace hier
{

class BoxGeometry;

class BoxOverlap;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class SideCellNoCornersFillPattern is a concrete implementation of the
 * abstract base class SAMRAI::xfer::VariableFillPattern.  It is used to
 * calculate overlaps according to a pattern which limits overlaps to the
 * cell-centered ghost region surrounding a patch, excluding all corners.  In
 * 3D, it is also possible to configure this fill pattern object also to exclude
 * all edges.
 */
class SideNoCornersFillPattern : public SAMRAI::xfer::VariableFillPattern
{
public:
    /*!
     * \brief Constructor.
     *
     * \note Parameters include_edges_on_dst_level and
     * include_edges_on_src_level have no effect for 2D problems.
     */
    SideNoCornersFillPattern(int stencil_width,
                             bool include_dst_patch_box,
                             bool include_edges_on_dst_level,
                             bool include_edges_on_src_level);

    /*!
     * \brief Destructor
     */
    ~SideNoCornersFillPattern();

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
     * \param src_mask            the source mask, the box resulting from shifting the source box
     * \param fill_box            the box to be filled
     * \param overwrite_interior  controls whether or not to include the destination box interior in the overlap
     * \param transformation      the transformation between source and destination index space
     *
     * \return                    pointer to the calculated overlap object
     */
    boost::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(const SAMRAI::hier::BoxGeometry& dst_geometry,
                     const SAMRAI::hier::BoxGeometry& src_geometry,
                     const SAMRAI::hier::Box& dst_patch_box,
                     const SAMRAI::hier::Box& src_mask,
                     const SAMRAI::hier::Box& fill_box,
                     bool overwrite_interior,
                     const SAMRAI::hier::Transformation& transformation) const;

    /*!
     * Compute overlaps that define the space to be filled by a refinement operation.
     *
     * \param fill_boxes          list representing all of the space on a patch or its ghost region that may be filled
     *                            by a refine operator (cell-centered represtentation)
     * \param node_fill_boxes     node-centered representation of fill_boxes
     * \param patch_box           box representing the patch where a refine operator will fill data (cell-centered
     *                            representation)
     * \param data_box            box representing the full extent of the region covered by a patch data object,
     *                            including all ghosts (cell-centered representation)
     * \param patch_data_factory  patch data factory for the data that is to be filled
     *
     * \return                    pointer to the calculated overlap object
     */
    boost::shared_ptr<SAMRAI::hier::BoxOverlap>
    computeFillBoxesOverlap(const SAMRAI::hier::BoxContainer& fill_boxes,
                            const SAMRAI::hier::BoxContainer& node_fill_boxes,
                            const SAMRAI::hier::Box& patch_box,
                            const SAMRAI::hier::Box& data_box,
                            const SAMRAI::hier::PatchDataFactory& patch_data_factory) const;

    /*!
     * Returns the stencil width.
     */
    SAMRAI::hier::IntVector& getStencilWidth();

    /*!
     * Returns a string name identifier "SIDE_NO_CORNERS_FILL_PATTERN".
     */
    const std::string& getPatternName() const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SideNoCornersFillPattern();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SideNoCornersFillPattern(const SideNoCornersFillPattern& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SideNoCornersFillPattern& operator=(const SideNoCornersFillPattern& that);

    SAMRAI::hier::IntVector d_stencil_width;
    const bool d_include_dst_patch_box;
    const bool d_include_edges_on_dst_level;
    const bool d_include_edges_on_src_level;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_SideNoCornersFillPattern
