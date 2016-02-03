// Filename: CartSideDoubleDivPreservingRefine.cpp
// Created on 09 Nov 2008 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <cmath>
#include <limits>
#include <ostream>

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CoarsenOperator.h"
#include "IBTK_config.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchDescriptor.h"
#include "PatchGeometry.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "ibtk/CartSideDoubleDivPreservingRefine.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Array.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define DIV_PRESERVING_CORRECTION_FC IBTK_FC_FUNC_(div_preserving_correction2d, DIV_PRESERVING_CORRECTION2D)
#endif

#if (NDIM == 3)
#define DIV_PRESERVING_CORRECTION_FC IBTK_FC_FUNC_(div_preserving_correction3d, DIV_PRESERVING_CORRECTION3D)
#endif

extern "C" {
void DIV_PRESERVING_CORRECTION_FC(double* u0,
                                  double* u1,
#if (NDIM == 3)
                                  double* u2,
#endif
                                  const int& u_gcw,
                                  const int& ilower0,
                                  const int& iupper0,
                                  const int& ilower1,
                                  const int& iupper1,
#if (NDIM == 3)
                                  const int& ilower2,
                                  const int& iupper2,
#endif
                                  const int& correction_box_ilower0,
                                  const int& correction_box_iupper0,
                                  const int& correction_box_ilower1,
                                  const int& correction_box_iupper1,
#if (NDIM == 3)
                                  const int& correction_box_ilower2,
                                  const int& correction_box_iupper2,
#endif
                                  const int* ratio,
                                  const double* dx_fine);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideDoubleDivPreservingRefine::CartSideDoubleDivPreservingRefine(const int u_dst_idx,
                                                                     const int u_src_idx,
                                                                     const int indicator_idx,
                                                                     Pointer<RefineOperator<NDIM> > refine_op,
                                                                     Pointer<CoarsenOperator<NDIM> > coarsen_op,
                                                                     const double fill_time,
                                                                     RefinePatchStrategy<NDIM>* const phys_bdry_op)
    : d_u_dst_idx(u_dst_idx),
      d_u_src_idx(u_src_idx),
      d_indicator_idx(indicator_idx),
      d_fill_time(fill_time),
      d_phys_bdry_op(phys_bdry_op),
      d_refine_op(refine_op),
      d_coarsen_op(coarsen_op)
{
    // intentionally blank
    return;
} // CartSideDoubleDivPreservingRefine

CartSideDoubleDivPreservingRefine::~CartSideDoubleDivPreservingRefine()
{
    // intentionally blank
    return;
} // ~CartSideDoubleDivPreservingRefine

void
CartSideDoubleDivPreservingRefine::setPhysicalBoundaryConditions(Patch<NDIM>& patch,
                                                                 const double fill_time,
                                                                 const IntVector<NDIM>& ghost_width_to_fill)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(fill_time, d_fill_time));
#endif
    if (d_phys_bdry_op) d_phys_bdry_op->setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    return;
} // setPhysicalBoundaryConditions

IntVector<NDIM>
CartSideDoubleDivPreservingRefine::getRefineOpStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
} // getRefineOpStencilWidth

void
CartSideDoubleDivPreservingRefine::preprocessRefine(Patch<NDIM>& /*fine*/,
                                                    const Patch<NDIM>& /*coarse*/,
                                                    const Box<NDIM>& /*fine_box*/,
                                                    const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
CartSideDoubleDivPreservingRefine::postprocessRefine(Patch<NDIM>& fine,
                                                     const Patch<NDIM>& coarse,
                                                     const Box<NDIM>& unrestricted_fine_box,
                                                     const IntVector<NDIM>& ratio)
{
    // NOTE: This operator cannot fill the full ghost cell width of the
    // destination data.  We instead restrict the size of the fine box to ensure
    // that we have adequate data to apply the divergence- and curl-preserving
    // corrections.
    const Box<NDIM> fine_box = unrestricted_fine_box * Box<NDIM>::grow(fine.getBox(), 2);

#if !defined(NDEBUG)
    for (int d = 0; d < NDIM; ++d)
    {
        if (ratio(d) % 2 != 0)
        {
            TBOX_ERROR("CartSideDoubleDivPreservingRefine::postprocessRefine():\n"
                       << "  refinement ratio must be a power of 2 for divergence- and "
                          "curl-preserving refinement operator."
                       << std::endl);
        }
    }
#endif

    Pointer<SideData<NDIM, double> > fdata = fine.getPatchData(d_u_dst_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(fdata);
#endif
    const int fdata_ghosts = fdata->getGhostCellWidth().max();
#if !defined(NDEBUG)
    TBOX_ASSERT(fdata_ghosts == fdata->getGhostCellWidth().min());
#endif
    const int fdata_depth = fdata->getDepth();

    Pointer<SideData<NDIM, double> > cdata = coarse.getPatchData(d_u_dst_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(cdata);
    const int cdata_ghosts = cdata->getGhostCellWidth().max();
    TBOX_ASSERT(cdata_ghosts == cdata->getGhostCellWidth().min());
    const int cdata_depth = cdata->getDepth();
    TBOX_ASSERT(cdata_depth == fdata_depth);
#endif

    if (ratio == IntVector<NDIM>(2))
    {
        // Perform (limited) conservative prolongation of the coarse grid data.
        d_refine_op->refine(fine, coarse, d_u_dst_idx, d_u_dst_idx, fine_box, ratio);

        Pointer<SideData<NDIM, double> > u_src_data = fine.getPatchData(d_u_src_idx);
        Pointer<SideData<NDIM, double> > indicator_data = fine.getPatchData(d_indicator_idx);

        // Ensure that we do not modify any of the data from the old level by
        // setting the value of the fine grid data to equal u_src wherever the
        // indicator data equals "1".
        if (u_src_data && indicator_data)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(fine_box, axis)); b; b++)
                {
                    const Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, axis, 0);
                    if (std::abs((*indicator_data)(i_s)-1.0) < 1.0e-12)
                    {
                        for (int depth = 0; depth < fdata_depth; ++depth)
                        {
                            (*fdata)(i_s, depth) = (*u_src_data)(i_s, depth);
                        }
                    }
                }
            }
        }

        // Reinterpolate data in the normal direction in the newly refined part
        // of the level wherever the indicator data does NOT equal "1".  Notice
        // that this loop actually modifies only data that is NOT covered by an
        // overlying coarse grid cell face.
        if (indicator_data)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(fine_box, axis)); b; b++)
                {
                    const Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, axis, 0);
                    if (!(std::abs((*indicator_data)(i_s)-1.0) < 1.0e-12))
                    {
                        const Index<NDIM> i_coarse_lower = IndexUtilities::coarsen(i, ratio);
                        const Index<NDIM> i_lower = IndexUtilities::refine(i_coarse_lower, ratio);
                        if (i(axis) == i_lower(axis)) continue;

                        Index<NDIM> i_coarse_upper = i_coarse_lower;
                        i_coarse_upper(axis) += 1;
                        const Index<NDIM> i_upper = IndexUtilities::refine(i_coarse_upper, ratio);

                        const double w1 =
                            static_cast<double>(i(axis) - i_lower(axis)) / static_cast<double>(ratio(axis));
                        const double w0 = 1.0 - w1;

                        const SideIndex<NDIM> i_s_lower(i_lower, axis, 0);
                        const SideIndex<NDIM> i_s_upper(i_upper, axis, 0);
                        for (int depth = 0; depth < fdata_depth; ++depth)
                        {
                            (*fdata)(i_s, depth) = w0 * (*fdata)(i_s_lower, depth) + w1* (*fdata)(i_s_upper, depth);
                        }
                    }
                }
            }
        }

        // Determine the box on which we need to compute the divergence- and
        // curl-preserving correction.
        const Box<NDIM> correction_box = Box<NDIM>::refine(Box<NDIM>::coarsen(fine_box, 2), 2);
#if !defined(NDEBUG)
        TBOX_ASSERT(fdata->getGhostBox().contains(correction_box));
#endif
        // Apply the divergence- and curl-preserving correction to the fine grid
        // data.
        Pointer<CartesianPatchGeometry<NDIM> > pgeom_fine = fine.getPatchGeometry();
        const double* const dx_fine = pgeom_fine->getDx();
        for (int d = 0; d < fdata_depth; ++d)
        {
            DIV_PRESERVING_CORRECTION_FC(fdata->getPointer(0, d),
                                         fdata->getPointer(1, d),
#if (NDIM == 3)
                                         fdata->getPointer(2, d),
#endif
                                         fdata_ghosts,
                                         fdata->getBox().lower()(0),
                                         fdata->getBox().upper()(0),
                                         fdata->getBox().lower()(1),
                                         fdata->getBox().upper()(1),
#if (NDIM == 3)
                                         fdata->getBox().lower()(2),
                                         fdata->getBox().upper()(2),
#endif
                                         correction_box.lower()(0),
                                         correction_box.upper()(0),
                                         correction_box.lower()(1),
                                         correction_box.upper()(1),
#if (NDIM == 3)
                                         correction_box.lower()(2),
                                         correction_box.upper()(2),
#endif
                                         ratio,
                                         dx_fine);
        }
    }
    else
    {
        // Setup an intermediate patch.
        const Box<NDIM> intermediate_patch_box = Box<NDIM>::refine(coarse.getBox(), 2);
        Patch<NDIM> intermediate(intermediate_patch_box, coarse.getPatchDescriptor());
        intermediate.allocatePatchData(d_u_dst_idx);

        // Setup a patch geometry object for the intermediate patch.
        Pointer<CartesianPatchGeometry<NDIM> > pgeom_coarse = coarse.getPatchGeometry();
        const IntVector<NDIM>& ratio_to_level_zero_coarse = pgeom_coarse->getRatio();
        Array<Array<bool> > touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            touches_regular_bdry[axis].resizeArray(2);
            touches_periodic_bdry[axis].resizeArray(2);
            for (int upperlower = 0; upperlower < 2; ++upperlower)
            {
                touches_regular_bdry[axis][upperlower] = pgeom_coarse->getTouchesRegularBoundary(axis, upperlower);
                touches_periodic_bdry[axis][upperlower] = pgeom_coarse->getTouchesPeriodicBoundary(axis, upperlower);
            }
        }
        const double* const dx_coarse = pgeom_coarse->getDx();

        const IntVector<NDIM> ratio_to_level_zero_intermediate = ratio_to_level_zero_coarse * 2;
        double dx_intermediate[NDIM], x_lower_intermediate[NDIM], x_upper_intermediate[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            dx_intermediate[d] = 0.5 * dx_coarse[d];
            x_lower_intermediate[d] = pgeom_coarse->getXLower()[d];
            x_upper_intermediate[d] = pgeom_coarse->getXUpper()[d];
        }
        intermediate.setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero_intermediate,
                                                                       touches_regular_bdry,
                                                                       touches_periodic_bdry,
                                                                       dx_intermediate,
                                                                       x_lower_intermediate,
                                                                       x_upper_intermediate));

        // The intermediate box where we need to fill data must be large enough
        // to provide ghost cell values for the fine fill box.
        const Box<NDIM> intermediate_box = Box<NDIM>::grow(Box<NDIM>::coarsen(fine_box, ratio / 2), 2);

        // Setup the original velocity and indicator data.
        if (fine.checkAllocated(d_u_src_idx) && fine.checkAllocated(d_indicator_idx))
        {
            intermediate.allocatePatchData(d_u_src_idx);
            intermediate.allocatePatchData(d_indicator_idx);
            Pointer<SideData<NDIM, double> > u_src_idata = intermediate.getPatchData(d_u_src_idx);
            Pointer<SideData<NDIM, double> > indicator_idata = intermediate.getPatchData(d_indicator_idx);
            u_src_idata->fillAll(std::numeric_limits<double>::quiet_NaN());
            indicator_idata->fillAll(-1.0);
#if !defined(NDEBUG)
            Pointer<SideData<NDIM, double> > u_src_fdata = fine.getPatchData(d_u_src_idx);
            Pointer<SideData<NDIM, double> > indicator_fdata = fine.getPatchData(d_indicator_idx);
            TBOX_ASSERT(u_src_fdata->getGhostBox().contains(Box<NDIM>::refine(intermediate_box, ratio / 2)));
            TBOX_ASSERT(indicator_fdata->getGhostBox().contains(Box<NDIM>::refine(intermediate_box, ratio / 2)));
#endif
            d_coarsen_op->coarsen(intermediate, fine, d_u_src_idx, d_u_src_idx, intermediate_box, ratio / 2);
            d_coarsen_op->coarsen(intermediate, fine, d_indicator_idx, d_indicator_idx, intermediate_box, ratio / 2);
        }

        // Recursively refine from the coarse patch to the fine patch.
        postprocessRefine(intermediate, coarse, intermediate_box, 2);
        postprocessRefine(fine, intermediate, fine_box, ratio / 2);

        // Deallocate any allocated patch data.
        intermediate.deallocatePatchData(d_u_dst_idx);
        if (fine.checkAllocated(d_u_src_idx) && fine.checkAllocated(d_indicator_idx))
        {
            intermediate.deallocatePatchData(d_u_src_idx);
            intermediate.deallocatePatchData(d_indicator_idx);
        }
    }
    return;
} // postprocessRefine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
