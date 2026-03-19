// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/CartCellDoubleLinearRefine.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellVariable.h>
#include <Patch.h>

#include <string>

#include <ibtk/namespaces.h> // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

#if (NDIM == 2)
#define CC_LINEAR_REFINE_FC IBTK_FC_FUNC_(cart_cell_linear_refine2d, CART_CELL_LINEAR_REFINE2D)
#endif
#if (NDIM == 3)
#define CC_LINEAR_REFINE_FC IBTK_FC_FUNC_(cart_cell_linear_refine3d, CART_CELL_LINEAR_REFINE3D)
#endif

extern "C"
{
    void CC_LINEAR_REFINE_FC(double*,
                             const int&,
                             const int&,
                             const int&,
                             const int&,
                             const int&,
#if (NDIM == 3)
                             const int&,
                             const int&,
#endif
                             const double*,
                             const int&,
                             const int&,
                             const int&,
                             const int&,
                             const int&,
#if (NDIM == 3)
                             const int&,
                             const int&,
#endif
                             const int&,
                             const int&,
                             const int&,
                             const int&,
#if (NDIM == 3)
                             const int&,
                             const int&,
#endif
                             const int*,
                             const int*);
}

namespace IBTK
{
const std::string CartCellDoubleLinearRefine::s_op_name = "IBTK_LINEAR_REFINE";

namespace
{
static const int REFINE_OP_PRIORITY = 0;
static const int REFINE_OP_STENCIL_WIDTH = 1;
} // namespace

bool
CartCellDoubleLinearRefine::findRefineOperator(const Pointer<Variable<NDIM>>& var, const std::string& op_name) const
{
    const Pointer<CellVariable<NDIM, double>> cc_var = var;
    return cc_var && op_name == s_op_name;
}

const std::string&
CartCellDoubleLinearRefine::getOperatorName() const
{
    return s_op_name;
}

int
CartCellDoubleLinearRefine::getOperatorPriority() const
{
    return REFINE_OP_PRIORITY;
}

IntVector<NDIM>
CartCellDoubleLinearRefine::getStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
}

void
CartCellDoubleLinearRefine::refine(Patch<NDIM>& fine,
                                   const Patch<NDIM>& coarse,
                                   const int dst_component,
                                   const int src_component,
                                   const Box<NDIM>& fine_box,
                                   const IntVector<NDIM>& ratio) const
{
    Pointer<CellData<NDIM, double>> fdata = fine.getPatchData(dst_component);
    Pointer<CellData<NDIM, double>> cdata = coarse.getPatchData(src_component);
#if !defined(NDEBUG)
    TBOX_ASSERT(fdata);
    TBOX_ASSERT(cdata);
    TBOX_ASSERT(fdata->getDepth() == cdata->getDepth());
#endif

    const int data_depth = fdata->getDepth();
    const Box<NDIM>& fdata_box = fdata->getBox();
    const int fdata_gcw = fdata->getGhostCellWidth().max();
    const Box<NDIM>& cdata_box = cdata->getBox();
    const int cdata_gcw = cdata->getGhostCellWidth().max();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = coarse.getPatchGeometry();
    int touches_regular_bdry[2 * NDIM];
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        touches_regular_bdry[2 * axis + 0] = pgeom->getTouchesRegularBoundary(axis, 0);
        touches_regular_bdry[2 * axis + 1] = pgeom->getTouchesRegularBoundary(axis, 1);
    }

    for (int depth = 0; depth < data_depth; ++depth)
    {
        CC_LINEAR_REFINE_FC(fdata->getPointer(depth),
                            fdata_gcw,
                            fdata_box.lower()(0),
                            fdata_box.upper()(0),
                            fdata_box.lower()(1),
                            fdata_box.upper()(1),
#if (NDIM == 3)
                            fdata_box.lower()(2),
                            fdata_box.upper()(2),
#endif
                            cdata->getPointer(depth),
                            cdata_gcw,
                            cdata_box.lower()(0),
                            cdata_box.upper()(0),
                            cdata_box.lower()(1),
                            cdata_box.upper()(1),
#if (NDIM == 3)
                            cdata_box.lower()(2),
                            cdata_box.upper()(2),
#endif
                            fine_box.lower()(0),
                            fine_box.upper()(0),
                            fine_box.lower()(1),
                            fine_box.upper()(1),
#if (NDIM == 3)
                            fine_box.lower()(2),
                            fine_box.upper()(2),
#endif
                            ratio,
                            touches_regular_bdry);
    }
}
} // namespace IBTK
