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

#include "ibtk/CartCellDoubleQuadraticRefine.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellVariable.h"
#include "Patch.h"

#include <array>
#include <string>
#include <vector>

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

const std::string CartCellDoubleQuadraticRefine::s_op_name = "QUADRATIC_REFINE";

namespace
{
static const int REFINE_OP_PRIORITY = 0;
static const int REFINE_OP_STENCIL_WIDTH = 1;

inline int
coarsen(const int& index, const int& ratio)
{
    return (index < 0 ? (index + 1) / ratio - 1 : index / ratio);
} // coarsen

inline hier::Index<NDIM>
coarsen(const hier::Index<NDIM>& index, const IntVector<NDIM>& ratio)
{
    hier::Index<NDIM> coarse_index;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        coarse_index(d) = coarsen(index(d), ratio(d));
    }
    return coarse_index;
} // coarsen
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
CartCellDoubleQuadraticRefine::findRefineOperator(const Pointer<Variable<NDIM> >& var, const std::string& op_name) const
{
    const Pointer<CellVariable<NDIM, double> > cc_var = var;
    return (cc_var && op_name == s_op_name);
} // findRefineOperator

const std::string&
CartCellDoubleQuadraticRefine::getOperatorName() const
{
    return s_op_name;
} // getOperatorName

int
CartCellDoubleQuadraticRefine::getOperatorPriority() const
{
    return REFINE_OP_PRIORITY;
} // getOperatorPriority

IntVector<NDIM>
CartCellDoubleQuadraticRefine::getStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
} // getStencilWidth

void
CartCellDoubleQuadraticRefine::refine(Patch<NDIM>& fine,
                                      const Patch<NDIM>& coarse,
                                      const int dst_component,
                                      const int src_component,
                                      const Box<NDIM>& fine_box,
                                      const IntVector<NDIM>& ratio) const
{
    // Get the patch data.
    Pointer<CellData<NDIM, double> > fdata = fine.getPatchData(dst_component);
    Pointer<CellData<NDIM, double> > cdata = coarse.getPatchData(src_component);
#if !defined(NDEBUG)
    TBOX_ASSERT(fdata);
    TBOX_ASSERT(cdata);
    TBOX_ASSERT(fdata->getDepth() == cdata->getDepth());
#endif
    const int data_depth = fdata->getDepth();

    const Box<NDIM>& patch_box_fine = fine.getBox();
    const hier::Index<NDIM>& patch_lower_fine = patch_box_fine.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom_fine = fine.getPatchGeometry();
    const double* const XLower_fine = pgeom_fine->getXLower();
    const double* const dx_fine = pgeom_fine->getDx();

    const Box<NDIM>& patch_box_crse = coarse.getBox();
    const hier::Index<NDIM>& patch_lower_crse = patch_box_crse.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom_crse = coarse.getPatchGeometry();
    const double* const XLower_crse = pgeom_crse->getXLower();
    const double* const dx_crse = pgeom_crse->getDx();

    // Set all values in the fine box via quadratic interpolation from the
    // overlying coarse grid data.
    for (Box<NDIM>::Iterator b(fine_box); b; b++)
    {
        const hier::Index<NDIM>& i_fine = b();
        const hier::Index<NDIM> i_crse = coarsen(i_fine, ratio);

        // Determine the interpolation stencil in the coarse index space.
        Box<NDIM> stencil_box_crse(i_crse, i_crse);
        stencil_box_crse.grow(IntVector<NDIM>(1));

        // Determine the interpolation weights.
        static const int degree = 2;
        std::array<std::array<double, degree + 1>, NDIM> wgts(array_constant<std::array<double, degree + 1>, NDIM>(
            std::array<double, degree + 1>(array_constant<double, degree + 1>(0.0))));
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const double X =
                XLower_fine[axis] + dx_fine[axis] * (static_cast<double>(i_fine(axis) - patch_lower_fine(axis)) + 0.5);
            std::vector<double> X_crse(degree + 1, 0.0);
            for (int i_crse = stencil_box_crse.lower()(axis), k = 0; i_crse <= stencil_box_crse.upper()(axis);
                 ++i_crse, ++k)
            {
                X_crse[k] =
                    XLower_crse[axis] + dx_crse[axis] * (static_cast<double>(i_crse - patch_lower_crse(axis)) + 0.5);
            }
            wgts[axis][0] = ((X - X_crse[1]) * (X - X_crse[2])) / ((X_crse[0] - X_crse[1]) * (X_crse[0] - X_crse[2]));
            wgts[axis][1] = ((X - X_crse[0]) * (X - X_crse[2])) / ((X_crse[1] - X_crse[0]) * (X_crse[1] - X_crse[2]));
            wgts[axis][2] = ((X - X_crse[0]) * (X - X_crse[1])) / ((X_crse[2] - X_crse[0]) * (X_crse[2] - X_crse[1]));
        }

        // Interpolate from the coarse grid to the fine grid.
        hier::Index<NDIM> i_intrp;
        for (int d = 0; d < data_depth; ++d)
        {
            (*fdata)(i_fine, d) = 0.0;
#if (NDIM > 2)
            for (int i2 = 0; i2 <= degree; ++i2)
            {
                const double& wgt2 = wgts[2][i2];
                i_intrp(2) = stencil_box_crse.lower()(2) + i2;
#else
            const double wgt2 = 1.0;
#endif
#if (NDIM > 1)
                for (int i1 = 0; i1 <= degree; ++i1)
                {
                    const double& wgt1 = wgts[1][i1];
                    i_intrp(1) = stencil_box_crse.lower()(1) + i1;
#else
            const double wgt1 = 1.0;
#endif
                    for (int i0 = 0; i0 <= degree; ++i0)
                    {
                        const double& wgt0 = wgts[0][i0];
                        i_intrp(0) = stencil_box_crse.lower()(0) + i0;

                        (*fdata)(i_fine, d) += wgt0 * wgt1 * wgt2 * (*cdata)(i_intrp, d);
                    }
#if (NDIM > 1)
                }
#endif
#if (NDIM > 2)
            }
#endif
        }
    }
    return;
} // refine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
