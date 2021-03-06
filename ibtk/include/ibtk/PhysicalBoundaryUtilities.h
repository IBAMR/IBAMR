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

#ifndef included_IBTK_PhysicalBoundaryUtilities
#define included_IBTK_PhysicalBoundaryUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

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
    PhysicalBoundaryUtilities() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PhysicalBoundaryUtilities(const PhysicalBoundaryUtilities& from) = delete;

    /*!
     * \brief Destructor.
     */
    ~PhysicalBoundaryUtilities() = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PhysicalBoundaryUtilities& operator=(const PhysicalBoundaryUtilities& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PhysicalBoundaryUtilities
