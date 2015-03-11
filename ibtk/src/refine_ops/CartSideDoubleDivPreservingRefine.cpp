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

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "IBTK_config.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/pdat/CellOverlap.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "ibtk/CartSideDoubleDivPreservingRefine.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include "SAMRAI/tbox/Utilities.h"

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
                                                                     boost::shared_ptr<RefineOperator> refine_op,
                                                                     boost::shared_ptr<CoarsenOperator> coarsen_op,
                                                                     const double fill_time,
                                                                     RefinePatchStrategy* const phys_bdry_op)
    : RefinePatchStrategy(DIM), d_u_dst_idx(u_dst_idx), d_u_src_idx(u_src_idx), d_indicator_idx(indicator_idx),
      d_fill_time(fill_time), d_phys_bdry_op(phys_bdry_op), d_refine_op(refine_op), d_coarsen_op(coarsen_op)
{
    // intentionally blank
    return;
}

CartSideDoubleDivPreservingRefine::~CartSideDoubleDivPreservingRefine()
{
    // intentionally blank
    return;
}

void CartSideDoubleDivPreservingRefine::setPhysicalBoundaryConditions(Patch& patch,
                                                                      const double fill_time,
                                                                      const IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(fill_time, d_fill_time));
    if (d_phys_bdry_op) d_phys_bdry_op->setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    return;
}

IntVector CartSideDoubleDivPreservingRefine::getRefineOpStencilWidth() const
{
    return IntVector(DIM, REFINE_OP_STENCIL_WIDTH);
}

void CartSideDoubleDivPreservingRefine::preprocessRefine(Patch& /*fine*/,
                                                         const Patch& /*coarse*/,
                                                         const Box& /*fine_box*/,
                                                         const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
}

void CartSideDoubleDivPreservingRefine::postprocessRefine(Patch& fine,
                                                          const Patch& coarse,
                                                          const Box& unrestricted_fine_box,
                                                          const IntVector& ratio)
{
    // NOTE: This operator cannot fill the full ghost cell width of the
    // destination data.  We instead restrict the size of the fine box to ensure
    // that we have adequate data to apply the divergence- and curl-preserving
    // corrections.
    const Box fine_box = unrestricted_fine_box * Box::grow(fine.getBox(), IntVector(DIM, 2));

    for (int d = 0; d < NDIM; ++d)
    {
        if (ratio(d) % 2 != 0)
        {
            TBOX_ERROR("CartSideDoubleDivPreservingRefine::postprocessRefine():\n"
                       << "  refinement ratio must be a power of 2 for divergence- and "
                          "curl-preserving refinement operator." << std::endl);
        }
    }

    auto fdata = BOOST_CAST<SideData<double> >(fine.getPatchData(d_u_dst_idx));
    TBOX_ASSERT(fdata);
    const int fdata_ghosts = fdata->getGhostCellWidth().max();
    TBOX_ASSERT(fdata_ghosts == fdata->getGhostCellWidth().min());
    const int fdata_depth = fdata->getDepth();

    auto cdata = BOOST_CAST<SideData<double> >(coarse.getPatchData(d_u_dst_idx));
    TBOX_ASSERT(cdata);
    const int cdata_ghosts = cdata->getGhostCellWidth().max();
    TBOX_ASSERT(cdata_ghosts == cdata->getGhostCellWidth().min());
    const int cdata_depth = cdata->getDepth();
    TBOX_ASSERT(cdata_depth == fdata_depth);

    if (ratio == IntVector(DIM, 2))
    {
        // Perform (limited) conservative prolongation of the coarse grid data.
        CellOverlap fine_overlap(BoxContainer(fine_box), IntVector::getZero(DIM)); // should this be SideOverlap?
        d_refine_op->refine(fine, coarse, d_u_dst_idx, d_u_dst_idx, fine_overlap, ratio);

        auto u_src_data = BOOST_CAST<SideData<double> >(fine.getPatchData(d_u_src_idx));
        auto indicator_data = BOOST_CAST<SideData<double> >(fine.getPatchData(d_indicator_idx));

        // Ensure that we do not modify any of the data from the old level by
        // setting the value of the fine grid data to equal u_src wherever the
        // indicator data equals "1".
        if (u_src_data && indicator_data)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator b(fine_box, axis); b; b++)
                {
                    const SideIndex& i_s = b();
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
                for (SideIterator b(fine_box, axis); b; b++)
                {
                    const SideIndex& i_s = b();
                    const CellIndex i(i_s.toCell(SideIndex::Upper));
                    if (!(std::abs((*indicator_data)(i_s)-1.0) < 1.0e-12))
                    {
                        const Index i_coarse_lower = IndexUtilities::coarsen(i, ratio);
                        const Index i_lower = IndexUtilities::refine(i_coarse_lower, ratio);
                        if (i(axis) == i_lower(axis)) continue;

                        Index i_coarse_upper = i_coarse_lower;
                        i_coarse_upper(axis) += 1;
                        const Index i_upper = IndexUtilities::refine(i_coarse_upper, ratio);

                        const double w1 =
                            static_cast<double>(i(axis) - i_lower(axis)) / static_cast<double>(ratio(axis));
                        const double w0 = 1.0 - w1;

                        const SideIndex i_s_lower(i_lower, axis, 0);
                        const SideIndex i_s_upper(i_upper, axis, 0);
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
        const Box correction_box = Box::refine(Box::coarsen(fine_box, IntVector(DIM, 2)), IntVector(DIM, 2));
        TBOX_ASSERT(fdata->getGhostBox().contains(correction_box));

        // Apply the divergence- and curl-preserving correction to the fine grid
        // data.
        auto pgeom_fine = BOOST_CAST<CartesianPatchGeometry>(fine.getPatchGeometry());
        const double* const dx_fine = pgeom_fine->getDx();
        for (int d = 0; d < fdata_depth; ++d)
        {
            DIV_PRESERVING_CORRECTION_FC(fdata->getPointer(0, d), fdata->getPointer(1, d),
#if (NDIM == 3)
                                         fdata->getPointer(2, d),
#endif
                                         fdata_ghosts, fdata->getBox().lower()(0), fdata->getBox().upper()(0),
                                         fdata->getBox().lower()(1), fdata->getBox().upper()(1),
#if (NDIM == 3)
                                         fdata->getBox().lower()(2), fdata->getBox().upper()(2),
#endif
                                         correction_box.lower()(0), correction_box.upper()(0),
                                         correction_box.lower()(1), correction_box.upper()(1),
#if (NDIM == 3)
                                         correction_box.lower()(2), correction_box.upper()(2),
#endif
                                         &ratio(0), dx_fine);
        }
    }
    else
    {
        // Setup an intermediate patch.
        Box intermediate_patch_box = Box::refine(coarse.getBox(), IntVector(DIM, 2));
        MappedBox intermediate_mapped_box(DIM);
        const GlobalId& global_id = coarse.getGlobalId();
        const LocalId& local_id = global_id.getLocalId();
        const int owner_rank = global_id.getOwnerRank();
        intermediate_mapped_box.initialize(intermediate_patch_box, local_id, owner_rank);
        Patch intermediate(intermediate_mapped_box, coarse.getPatchDescriptor());
        intermediate.allocatePatchData(d_u_dst_idx);

        // Setup a patch geometry object for the intermediate patch.
        auto pgeom_crse = BOOST_CAST<CartesianPatchGeometry>(coarse.getPatchGeometry());
        const IntVector& ratio_to_level_zero_coarse = pgeom_crse->getRatio();
        PatchGeometry::TwoDimBool touches_regular_bdry(DIM), touches_periodic_bdry(DIM);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (int side = 0; side < 2; ++side)
            {
                touches_regular_bdry(axis, side) = pgeom_crse->getTouchesRegularBoundary(axis, side);
                touches_periodic_bdry(axis, side) = pgeom_crse->getTouchesPeriodicBoundary(axis, side);
            }
        }
        const double* const dx_coarse = pgeom_crse->getDx();

        const IntVector ratio_to_level_zero_intermediate = ratio_to_level_zero_coarse * 2;
        double dx_intermediate[NDIM], x_lower_intermediate[NDIM], x_upper_intermediate[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            dx_intermediate[d] = 0.5 * dx_coarse[d];
            x_lower_intermediate[d] = pgeom_crse->getXLower()[d];
            x_upper_intermediate[d] = pgeom_crse->getXUpper()[d];
        }
        intermediate.setPatchGeometry(boost::make_shared<CartesianPatchGeometry>(
            ratio_to_level_zero_intermediate, touches_regular_bdry, touches_periodic_bdry, dx_intermediate,
            x_lower_intermediate, x_upper_intermediate));

        // The intermediate box where we need to fill data must be large enough
        // to provide ghost cell values for the fine fill box.
        const Box intermediate_box = Box::grow(Box::coarsen(fine_box, ratio / 2), IntVector(DIM, 2));

        // Setup the original velocity and indicator data.
        if (fine.checkAllocated(d_u_src_idx) && fine.checkAllocated(d_indicator_idx))
        {
            intermediate.allocatePatchData(d_u_src_idx);
            intermediate.allocatePatchData(d_indicator_idx);
            auto u_src_idata = BOOST_CAST<SideData<double> >(intermediate.getPatchData(d_u_src_idx));
            auto indicator_idata = BOOST_CAST<SideData<double> >(intermediate.getPatchData(d_indicator_idx));
            u_src_idata->fillAll(std::numeric_limits<double>::quiet_NaN());
            indicator_idata->fillAll(-1.0);
            auto u_src_fdata = BOOST_CAST<SideData<double> >(fine.getPatchData(d_u_src_idx));
            auto indicator_fdata = BOOST_CAST<SideData<double> >(fine.getPatchData(d_indicator_idx));
            TBOX_ASSERT(u_src_fdata->getGhostBox().contains(Box::refine(intermediate_box, ratio / 2)));
            TBOX_ASSERT(indicator_fdata->getGhostBox().contains(Box::refine(intermediate_box, ratio / 2)));
            d_coarsen_op->coarsen(intermediate, fine, d_u_src_idx, d_u_src_idx, intermediate_box, ratio / 2);
            d_coarsen_op->coarsen(intermediate, fine, d_indicator_idx, d_indicator_idx, intermediate_box, ratio / 2);
        }

        // Recursively refine from the coarse patch to the fine patch.
        postprocessRefine(intermediate, coarse, intermediate_box, IntVector(DIM, 2));
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
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
