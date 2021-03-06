// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/CartCellDoubleBoundsPreservingConservativeLinearRefine.h"

#include "Box.h"
#include "BoxList.h"
#include "CartesianPatchGeometry.h"
#include "CellVariable.h"
#include "Patch.h"

#include <algorithm>
#include <limits>
#include <string>

#include "ibtk/namespaces.h" // IWYU pragma: keep

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
    const hier::Index<NDIM>& patch_lower_crse = patch_box_crse.lower();
    const hier::Index<NDIM>& patch_upper_crse = patch_box_crse.upper();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom_crse = coarse.getPatchGeometry();
    for (int depth = 0; depth < data_depth; ++depth)
    {
        for (Box<NDIM>::Iterator b(coarse_correction_box); b; b++)
        {
            const hier::Index<NDIM>& i_crse = b();
            const hier::Index<NDIM> i_fine = i_crse * ratio;

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
            double u = -std::numeric_limits<double>::max();
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
