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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_ShiftedPatchGeometry
#define included_IBTK_ShiftedPatchGeometry

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <Patch.h>
#include <PatchGeometry.h>

#include <array>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class ShiftedPatchGeometry is the Cartesian patch geometry of a patch, displaced by a shift.
 *
 * Boundary condition objects locate their coefficients with the geometry of the patch they receive. Replacing the
 * geometry of the patch with a ShiftedPatchGeometry makes them evaluate the coefficients at the locations of data with
 * another centering, such as the sides normal to an axis.
 */
class ShiftedPatchGeometry
{
public:
    /*!
     * \brief Constructor. The geometry has the grid spacing, refinement ratio, and boundary flags of the geometry of \a
     * patch, and its lower and upper corners are displaced by \a shift.
     */
    ShiftedPatchGeometry(const SAMRAI::hier::Patch<NDIM>& patch, const std::array<double, NDIM>& shift);

    /*!
     * \brief Class Scope replaces the geometry of a patch by a ShiftedPatchGeometry until the scope ends.
     *
     * \note SAMRAI::solv::RobinBcCoefStrategy::setBcCoefs() takes a const patch, so the constructor sets the geometry
     * of a const patch. The destructor restores the original geometry.
     */
    class Scope
    {
    public:
        /*!
         * \brief Constructor.
         */
        Scope(const SAMRAI::hier::Patch<NDIM>& patch, const ShiftedPatchGeometry& shifted);

        /*!
         * \brief Destructor.
         */
        ~Scope();

        Scope(const Scope& from) = delete;
        Scope& operator=(const Scope& that) = delete;

    private:
        SAMRAI::hier::Patch<NDIM>& d_patch;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchGeometry<NDIM>> d_original;
    };

private:
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchGeometry<NDIM>> d_geometry;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_ShiftedPatchGeometry
