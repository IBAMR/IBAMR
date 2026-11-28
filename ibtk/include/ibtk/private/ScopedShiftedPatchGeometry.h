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

#ifndef included_IBTK_ScopedShiftedPatchGeometry
#define included_IBTK_ScopedShiftedPatchGeometry

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <Patch.h>
#include <PatchGeometry.h>

#include <array>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace detail
{
/*!
 * \brief Class ScopedShiftedPatchGeometry replaces the geometry of a patch by the Cartesian geometry displaced by a
 * shift until the scope ends, so that boundary condition coefficients are evaluated at the locations of data with
 * another centering, such as the sides normal to an axis.
 */
class ScopedShiftedPatchGeometry
{
public:
    /*!
     * \brief Constructor. The shifted geometry has the grid spacing, refinement ratio, and boundary flags of \a patch.
     * Its lower and upper corners are displaced by \a shift.
     */
    ScopedShiftedPatchGeometry(SAMRAI::hier::Patch<NDIM>& patch, const std::array<double, NDIM>& shift);

    /*!
     * \brief Destructor. Restores the original geometry of the patch.
     */
    ~ScopedShiftedPatchGeometry();

    ScopedShiftedPatchGeometry(const ScopedShiftedPatchGeometry& from) = delete;
    ScopedShiftedPatchGeometry& operator=(const ScopedShiftedPatchGeometry& that) = delete;

private:
    /*! The patch and its original geometry, restored by the destructor. */
    SAMRAI::hier::Patch<NDIM>& d_patch;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchGeometry<NDIM>> d_original_geometry;
};
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_ScopedShiftedPatchGeometry
