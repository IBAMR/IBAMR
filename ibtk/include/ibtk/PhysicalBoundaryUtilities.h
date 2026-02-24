// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2023 by the IBAMR developers
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

#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIArray.h>
#include <SAMRAIBoundaryBox.h>
#include <SAMRAIBox.h>
#include <SAMRAIPatch.h>

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
    static SAMRAIArray<SAMRAIBoundaryBox> getPhysicalBoundaryCodim1Boxes(const SAMRAIPatch& patch);

    /*!
     * \brief Return the co-dimension 2 boundary boxes corresponding to the
     * physical boundaries of the supplied patch.
     */
    static SAMRAIArray<SAMRAIBoundaryBox> getPhysicalBoundaryCodim2Boxes(const SAMRAIPatch& patch);

    /*!
     * \brief Return the co-dimension 3 boundary boxes corresponding to the
     * physical boundaries of the supplied patch.
     */
    static SAMRAIArray<SAMRAIBoundaryBox> getPhysicalBoundaryCodim3Boxes(const SAMRAIPatch& patch);

    /*!
     * \brief Trim a co-dimension 1 boundary box so that it does not stick out
     * past the patch domain in directions transverse to the boundary normal.
     *
     * \note The supplied boundary box must be of type 1.
     *
     * \see SAMRAI::hier::BoundaryBox::getBoundaryType
     */
    static SAMRAIBoundaryBox trimBoundaryCodim1Box(const SAMRAIBoundaryBox& bdry_box, const SAMRAIPatch& patch);

    /*!
     * \brief Return box describing the side index space of surfaces defined by
     * a boundary box.
     *
     * \note The supplied boundary box must be of type 1.
     *
     * \see SAMRAI::hier::BoundaryBox::getBoundaryType
     */
    static SAMRAIBox makeSideBoundaryCodim1Box(const SAMRAIBoundaryBox& boundary_box);

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

#endif // #ifndef included_IBTK_PhysicalBoundaryUtilities
