// Filename: CartCellDoubleBoundsPreservingConservativeLinearRefine.cpp
// Created on 06 Jul 2010 by Boyce Griffith
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

#include <algorithm>
#include <limits>
#include <ostream>
#include <string>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "ibtk/CartCellDoubleBoundsPreservingConservativeLinearRefine.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{
class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string CartCellDoubleBoundsPreservingConservativeLinearRefine::s_op_name =
    "BOUNDS_PRESERVING_CONSERVATIVE_LINEAR_REFINE";

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartCellDoubleBoundsPreservingConservativeLinearRefine::CartCellDoubleBoundsPreservingConservativeLinearRefine()
    : RefineOperator(s_op_name), d_conservative_linear_refine_op(), d_constant_refine_op()
{
    // intentionally blank
    return;
}

CartCellDoubleBoundsPreservingConservativeLinearRefine::~CartCellDoubleBoundsPreservingConservativeLinearRefine()
{
    // intentionally blank
    return;
}

int CartCellDoubleBoundsPreservingConservativeLinearRefine::getOperatorPriority() const
{
    return d_conservative_linear_refine_op.getOperatorPriority();
}

IntVector CartCellDoubleBoundsPreservingConservativeLinearRefine::getStencilWidth(const Dimension& dim) const
{
    return d_conservative_linear_refine_op.getStencilWidth(dim);
}

void CartCellDoubleBoundsPreservingConservativeLinearRefine::refine(Patch& fine,
                                                                    const Patch& coarse,
                                                                    const int dst_component,
                                                                    const int src_component,
                                                                    const BoxOverlap& fine_overlap,
                                                                    const IntVector& ratio) const
{
    auto fine_cell_overlap = CPP_CAST<const CellOverlap*>(&fine_overlap);
    TBOX_ASSERT(fine_cell_overlap);
    const BoxContainer& fine_boxes = fine_cell_overlap->getDestinationBoxList();
    for (auto bl = fine_boxes.begin(), e = fine_boxes.end(); bl != e; ++bl)
    {
        const Box& fine_box = bl();
        // Determine the box over which we can apply the bounds-preserving
        // correction, and construct a list of boxes that will not be corrected.
        bool empty_correction_box = false;
        Box correction_box = Box::refine(Box::coarsen(fine_box, ratio), ratio);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            int& lower = correction_box.lower(axis);
            while (lower < fine_box.lower(axis))
            {
                lower += ratio(axis);
            }

            int& upper = correction_box.upper(axis);
            while (upper > fine_box.upper(axis))
            {
                upper -= ratio(axis);
            }

            if (lower >= upper)
            {
                empty_correction_box = true;
            }
        }
        const Box coarse_correction_box = Box::coarsen(correction_box, ratio);
        BoxContainer uncorrected_boxes(fine_box);
        if (!empty_correction_box)
        {
            uncorrected_boxes.removeIntersections(correction_box);
        }

        // Employ limited conservative interpolation to prolong data on the
        // correction box.
        d_conservative_linear_refine_op.refine(fine, coarse, dst_component, src_component, correction_box, ratio);

        // Employ constant interpolation to prolong data on the rest of the fine
        // box.
        for (auto b(uncorrected_boxes); b; b++)
        {
            d_constant_refine_op.refine(fine, coarse, dst_component, src_component, b(), ratio);
        }

        // There is nothing left to do if the correction box is empty.
        if (empty_correction_box) return;

        // Correct the data within the correction box.
        auto fdata = BOOST_CAST<CellData<double> >(fine.getPatchData(dst_component));
        auto cdata = BOOST_CAST<CellData<double> >(coarse.getPatchData(src_component));
        TBOX_ASSERT(fdata);
        TBOX_ASSERT(cdata);
        TBOX_ASSERT(fdata->getDepth() == cdata->getDepth());
        const int data_depth = fdata->getDepth();
        const Box& patch_box_crse = coarse.getBox();
        const Index& patch_lower_crse = patch_box_crse.lower();
        const Index& patch_upper_crse = patch_box_crse.upper();
        auto pgeom_crse = BOOST_CAST<CartesianPatchGeometry>(coarse.getPatchGeometry());
        for (int depth = 0; depth < data_depth; ++depth)
        {
            for (auto b(coarse_correction_box); b; b++)
            {
                const Index& i_crse = b();
                const Index i_fine = i_crse * ratio;

                // Determine the lower/upper bounds.
                Box stencil_box_crse(i_crse, i_crse);
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    if (i_crse(axis) > patch_lower_crse(axis) || !pgeom_crse->getTouchesRegularBoundary(axis, 0))
                    {
                        stencil_box_crse.growLower(axis, 1);
                    }
                    if (i_crse(axis) < patch_upper_crse(axis) || !pgeom_crse->getTouchesRegularBoundary(axis, 1))
                    {
                        stencil_box_crse.growUpper(axis, 1);
                    }
                }

                double l = std::numeric_limits<double>::max();
                double u = -(l - std::numeric_limits<double>::epsilon());
                for (auto b = CellGeometry::begin(stencil_box_crse), e = CellGeometry::end(stencil_box_crse); b != e;
                     ++b)
                {
                    const double& m = (*cdata)(b(), depth);
                    l = std::min(l, m);
                    u = std::max(u, m);
                }

                // Force all refined data to lie within the bounds, accumulating the
                // discrepancy.
                Box stencil_box_fine(i_fine, i_fine);
                stencil_box_fine.growUpper(ratio - IntVector::getOne(DIM));
                double Delta = 0.0;
                for (auto b = CellGeometry::begin(stencil_box_fine), e = CellGeometry::end(stencil_box_fine); b != e;
                     ++b)
                {
                    double& m = (*fdata)(b(), depth);
                    Delta += std::max(0.0, m - u) - std::max(0.0, l - m);
                    m = std::max(std::min(m, u), l);
                }

                // Distribute the discrepancy to maintain conservation.
                if (Delta >= std::numeric_limits<double>::epsilon())
                {
                    double K = 0.0;
                    for (auto b = CellGeometry::begin(stencil_box_fine), e = CellGeometry::end(stencil_box_fine);
                         b != e; ++b)
                    {
                        const double& m = (*fdata)(b(), depth);
                        double k = u - m;
                        K += k;
                    }
                    for (auto b = CellGeometry::begin(stencil_box_fine), e = CellGeometry::end(stencil_box_fine);
                         b != e; ++b)
                    {
                        double& m = (*fdata)(b(), depth);
                        double k = u - m;
                        m += Delta * k / K;
                    }
                }
                else if (Delta <= -std::numeric_limits<double>::epsilon())
                {
                    double K = 0.0;
                    for (auto b = CellGeometry::begin(stencil_box_fine), e = CellGeometry::end(stencil_box_fine);
                         b != e; ++b)
                    {
                        const double& m = (*fdata)(b(), depth);
                        double k = m - l;
                        K += k;
                    }
                    for (auto b = CellGeometry::begin(stencil_box_fine), e = CellGeometry::end(stencil_box_fine);
                         b != e; ++b)
                    {
                        double& m = (*fdata)(b(), depth);
                        double k = m - l;
                        m += Delta * k / K;
                    }
                }
            }
        }
    }
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
