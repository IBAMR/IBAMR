// Filename: PhysicalBoundaryUtilities.h
// Created on 30 Sep 2006 by Boyce Griffith
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

#ifndef included_PhysicalBoundaryUtilities
#define included_PhysicalBoundaryUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "BoundaryBox.h"
#include "Box.h"
#include "tbox/Array.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PhysicalBoundaryUtilities is a utility class to organize
 * functions related to setting values at physical boundaries.
 */
class PhysicalBoundaryUtilities
{
public:
    /*!
     * \brief Indicate whether the given boundary information indicates a lower
     * boundary region (i.e., the associated box region contains higher values
     * along the axis in the coordinate direction than the boundary region).
     */
    static bool isLower(int loc, int codim, int direction);

    /*!
     * \brief Indicate whether the given boundary information indicates a upper
     * boundary region (i.e., the associated box region contains lower values
     * along the axis in the coordinate direction than the boundary region).
     */
    static bool isUpper(int loc, int codim, int direction);

    /*!
     * \brief Return the co-dimension 1 boundary boxes corresponding to the
     * physical boundaries of the supplied patch.
     */
    static SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> >
    getPhysicalBoundaryCodim1Boxes(const SAMRAI::hier::Patch<NDIM>& patch);

    /*!
     * \brief Return the co-dimension 2 boundary boxes corresponding to the
     * physical boundaries of the supplied patch.
     */
    static SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> >
    getPhysicalBoundaryCodim2Boxes(const SAMRAI::hier::Patch<NDIM>& patch);

    /*!
     * \brief Return the co-dimension 3 boundary boxes corresponding to the
     * physical boundaries of the supplied patch.
     */
    static SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> >
    getPhysicalBoundaryCodim3Boxes(const SAMRAI::hier::Patch<NDIM>& patch);

    /*!
     * \brief Trim a co-dimension 1 boundary box so that it does not stick out
     * past the patch domain in directions transverse to the boundary normal.
     *
     * \note The supplied boundary box must be of type 1.
     *
     * \see SAMRAI::hier::BoundaryBox::getBoundaryType
     */
    static SAMRAI::hier::BoundaryBox<NDIM> trimBoundaryCodim1Box(const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                                                                 const SAMRAI::hier::Patch<NDIM>& patch);

    /*!
     * \brief Return box describing the side index space of surfaces defined by
     * a boundary box.
     *
     * \note The supplied boundary box must be of type 1.
     *
     * \see SAMRAI::hier::BoundaryBox::getBoundaryType
     */
    static SAMRAI::hier::Box<NDIM> makeSideBoundaryCodim1Box(const SAMRAI::hier::BoundaryBox<NDIM>& boundary_box);

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PhysicalBoundaryUtilities();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PhysicalBoundaryUtilities(const PhysicalBoundaryUtilities& from);

    /*!
     * \brief Destructor.
     */
    ~PhysicalBoundaryUtilities();

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PhysicalBoundaryUtilities& operator=(const PhysicalBoundaryUtilities& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PhysicalBoundaryUtilities
