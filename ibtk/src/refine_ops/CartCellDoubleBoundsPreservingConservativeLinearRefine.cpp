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

#include "Box.h"
#include "BoxList.h"
#include "CartesianCellDoubleConservativeLinearRefine.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellDoubleConstantRefine.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "ibtk/CartCellDoubleBoundsPreservingConservativeLinearRefine.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
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
    : d_conservative_linear_refine_op(), d_constant_refine_op()
{
    // intentionally blank
    return;
} // CartCellDoubleBoundsPreservingConservativeLinearRefine

CartCellDoubleBoundsPreservingConservativeLinearRefine::~CartCellDoubleBoundsPreservingConservativeLinearRefine()
{
    // intentionally blank
    return;
} // ~CartCellDoubleBoundsPreservingConservativeLinearRefine

bool
CartCellDoubleBoundsPreservingConservativeLinearRefine::findRefineOperator(const Pointer<Variable<NDIM> >& var,
                                                                           const std::string& op_name) const
{
    const Pointer<CellVariable<NDIM, double> > cc_var = var;
    return (cc_var && op_name == s_op_name);
} // findRefineOperator

const std::string&
CartCellDoubleBoundsPreservingConservativeLinearRefine::getOperatorName() const
{
    return s_op_name;
} // getOperatorName

int
CartCellDoubleBoundsPreservingConservativeLinearRefine::getOperatorPriority() const
{
    return d_conservative_linear_refine_op.getOperatorPriority();
} // getOperatorPriority

IntVector<NDIM>
CartCellDoubleBoundsPreservingConservativeLinearRefine::getStencilWidth() const
{
    return d_conservative_linear_refine_op.getStencilWidth();
} // getStencilWidth

void
CartCellDoubleBoundsPreservingConservativeLinearRefine::refine(Patch<NDIM>& fine,
                                                               const Patch<NDIM>& coarse,
                                                               const int dst_component,
                                                               const int src_component,
                                                               const Box<NDIM>& fine_box,
                                                               const IntVector<NDIM>& ratio) const
{
    // Determine the box over which we can apply the bounds-preserving
    // correction, and construct a list of boxes that will not be corrected.
    bool empty_correction_box = false;
    Box<NDIM> correction_box = Box<NDIM>::refine(Box<NDIM>::coarsen(fine_box, ratio), ratio);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        int& lower = correction_box.lower()(axis);
        while (lower < fine_box.lower()(axis))
        {
            lower += ratio(axis);
        }

        int& upper = correction_box.upper()(axis);
        while (upper > fine_box.upper()(axis))
        {
            upper -= ratio(axis);
        }

        if (lower >= upper)
        {
            empty_correction_box = true;
        }
    }
    const Box<NDIM> coarse_correction_box = Box<NDIM>::coarsen(correction_box, ratio);

    BoxList<NDIM> uncorrected_boxes(fine_box);
    if (!empty_correction_box)
    {
        uncorrected_boxes.removeIntersections(correction_box);
    }

    // Employ limited conservative interpolation to prolong data on the
    // correction box.
    d_conservative_linear_refine_op.refine(fine, coarse, dst_component, src_component, correction_box, ratio);

    // Employ constant interpolation to prolong data on the rest of the fine
    // box.
    for (BoxList<NDIM>::Iterator b(uncorrected_boxes); b; b++)
    {
        d_constant_refine_op.refine(fine, coarse, dst_component, src_component, b(), ratio);
    }

    // There is nothing left to do if the correction box is empty.
    if (empty_correction_box) return;

    // Correct the data within the correction box.
    Pointer<CellData<NDIM, double> > fdata = fine.getPatchData(dst_component);
    Pointer<CellData<NDIM, double> > cdata = coarse.getPatchData(src_component);
#if !defined(NDEBUG)
    TBOX_ASSERT(fdata);
    TBOX_ASSERT(cdata);
    TBOX_ASSERT(fdata->getDepth() == cdata->getDepth());
#endif
    const int data_depth = fdata->getDepth();
    const Box<NDIM>& patch_box_crse = coarse.getBox();
    const Index<NDIM>& patch_lower_crse = patch_box_crse.lower();
    const Index<NDIM>& patch_upper_crse = patch_box_crse.upper();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom_crse = coarse.getPatchGeometry();
    for (int depth = 0; depth < data_depth; ++depth)
    {
        for (Box<NDIM>::Iterator b(coarse_correction_box); b; b++)
        {
            const Index<NDIM>& i_crse = b();
            const Index<NDIM> i_fine = i_crse * ratio;

            // Determine the lower/upper bounds.
            Box<NDIM> stencil_box_crse(i_crse, i_crse);
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
            for (Box<NDIM>::Iterator b(stencil_box_crse); b; b++)
            {
                const double& m = (*cdata)(b(), depth);
                l = std::min(l, m);
                u = std::max(u, m);
            }

            // Force all refined data to lie within the bounds, accumulating the
            // discrepancy.
            Box<NDIM> stencil_box_fine(i_fine, i_fine);
            stencil_box_fine.growUpper(ratio - IntVector<NDIM>(1));
            double Delta = 0.0;
            for (Box<NDIM>::Iterator b(stencil_box_fine); b; b++)
            {
                double& m = (*fdata)(b(), depth);
                Delta += std::max(0.0, m - u) - std::max(0.0, l - m);
                m = std::max(std::min(m, u), l);
            }

            // Distribute the discrepancy to maintain conservation.
            if (Delta >= std::numeric_limits<double>::epsilon())
            {
                double K = 0.0;
                for (Box<NDIM>::Iterator b(stencil_box_fine); b; b++)
                {
                    const double& m = (*fdata)(b(), depth);
                    double k = u - m;
                    K += k;
                }
                for (Box<NDIM>::Iterator b(stencil_box_fine); b; b++)
                {
                    double& m = (*fdata)(b(), depth);
                    double k = u - m;
                    m += Delta * k / K;
                }
            }
            else if (Delta <= -std::numeric_limits<double>::epsilon())
            {
                double K = 0.0;
                for (Box<NDIM>::Iterator b(stencil_box_fine); b; b++)
                {
                    const double& m = (*fdata)(b(), depth);
                    double k = m - l;
                    K += k;
                }
                for (Box<NDIM>::Iterator b(stencil_box_fine); b; b++)
                {
                    double& m = (*fdata)(b(), depth);
                    double k = m - l;
                    m += Delta * k / K;
                }
            }
        }
    }
    return;
} // refine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
