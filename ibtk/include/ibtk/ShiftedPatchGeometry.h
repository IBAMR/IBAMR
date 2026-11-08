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
 * \brief Class ShiftedPatchGeometry is the Cartesian patch geometry of a patch displaced by a shift.
 *
 * Use a ScopedGeometryReplacement to replace the geometry of a patch, so that boundary condition coefficients are
 * evaluated at the locations of data with another centering, such as the sides normal to an axis.
 */
class ShiftedPatchGeometry
{
public:
    /*!
     * \brief Constructor. The geometry has the grid spacing, refinement ratio, and boundary flags of \a patch.
     * Its lower and upper corners are displaced by \a shift.
     */
    ShiftedPatchGeometry(SAMRAI::hier::Patch<NDIM>& patch, const std::array<double, NDIM>& shift);

    /*!
     * \brief Class ScopedGeometryReplacement replaces the geometry of a patch by a ShiftedPatchGeometry until the scope
     * ends.
     */
    class ScopedGeometryReplacement
    {
    public:
        /*!
         * \brief Constructor.
         */
        ScopedGeometryReplacement(SAMRAI::hier::Patch<NDIM>& patch, const ShiftedPatchGeometry& shifted_geometry);

        /*!
         * \brief Destructor.
         */
        ~ScopedGeometryReplacement();

        ScopedGeometryReplacement(const ScopedGeometryReplacement& from) = delete;
        ScopedGeometryReplacement& operator=(const ScopedGeometryReplacement& that) = delete;

    private:
        /*! The patch and its original geometry, restored by the destructor. */
        SAMRAI::hier::Patch<NDIM>& d_patch;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchGeometry<NDIM>> d_original_geometry;
    };

private:
    /*! The shifted geometry. */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchGeometry<NDIM>> d_shifted_geometry;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_ShiftedPatchGeometry
