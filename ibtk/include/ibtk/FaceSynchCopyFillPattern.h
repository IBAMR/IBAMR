// Filename: FaceSynchCopyFillPattern.h
// Created on 03 Feb 2011 by Boyce Griffith
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

#ifndef included_FaceSynchCopyFillPattern
#define included_FaceSynchCopyFillPattern

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

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
 * \brief Class FaceCellSynchCopyFillPattern is a concrete implementation of the
 * abstract base class SAMRAI::xfer::VariableFillPattern.  It is used to
 * calculate overlaps according to a pattern which limits overlaps to the
 * face-centered ghost region surrounding a patch appropriate for
 * "synchronizing" face-centered values at patch boundaries.
 */
class FaceSynchCopyFillPattern : public SAMRAI::xfer::VariableFillPattern<NDIM>
{
public:
    /*!
     * \brief Default constructor
     */
    FaceSynchCopyFillPattern();

    /*!
     * \brief Destructor
     */
    ~FaceSynchCopyFillPattern();

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
                     const SAMRAI::hier::IntVector<NDIM>& src_offset) const;

    /*!
     * Returns the stencil width.
     */
    SAMRAI::hier::IntVector<NDIM>& getStencilWidth();

    /*!
     * Returns a string name identifier "FACE_SYNCH_COPY_FILL_PATTERN".
     */
    const std::string& getPatternName() const;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FaceSynchCopyFillPattern(const FaceSynchCopyFillPattern& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FaceSynchCopyFillPattern& operator=(const FaceSynchCopyFillPattern& that);

    SAMRAI::hier::IntVector<NDIM> d_stencil_width;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FaceSynchCopyFillPattern
