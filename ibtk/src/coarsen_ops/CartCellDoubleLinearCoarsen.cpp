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

#include <ibtk/CartCellDoubleLinearCoarsen.h>
#include <ibtk/PhysicalBoundaryUtilities.h>

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
#define CC_LINEAR_COARSEN_FC IBTK_FC_FUNC(cart_cell_linear_coarsen2d, CART_CELL_LINEAR_COARSEN2D)
#define CC_LINEAR_COARSEN_BDRY_FC IBTK_FC_FUNC(cart_cell_linear_coarsen_bdry2d, CART_CELL_LINEAR_COARSEN_BDRY2D)
#endif
#if (NDIM == 3)
#define CC_LINEAR_COARSEN_FC IBTK_FC_FUNC(cart_cell_linear_coarsen3d, CART_CELL_LINEAR_COARSEN3D)
#define CC_LINEAR_COARSEN_BDRY_FC IBTK_FC_FUNC(cart_cell_linear_coarsen_bdry3d, CART_CELL_LINEAR_COARSEN_BDRY3D)
#endif

extern "C"
{
    void CC_LINEAR_COARSEN_FC(double*,
                              const int&,
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
                              const int&,
                              const int&,
                              const int&,
                              const int&,
#if (NDIM == 3)
                              const int&,
                              const int&,
#endif
                              const int*);
    void CC_LINEAR_COARSEN_BDRY_FC(double*,
                                   const int&,
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
                                   const int&,
                                   const int&,
                                   const int&,
                                   const int&,
#if (NDIM == 3)
                                   const int&,
                                   const int&,
#endif
                                   const int*,
                                   const int*,
                                   const int&,
                                   const int&,
                                   const int&,
                                   const int&,
#if (NDIM == 3)
                                   const int&,
                                   const int&,
#endif
                                   const int&);
}

namespace IBTK
{
const std::string CartCellDoubleLinearCoarsen::s_op_name = "IBTK_LINEAR_COARSEN";

namespace
{
static const int COARSEN_OP_PRIORITY = 0;
static const int COARSEN_OP_STENCIL_WIDTH = 1;
} // namespace

bool
CartCellDoubleLinearCoarsen::findCoarsenOperator(const Pointer<Variable<NDIM>>& var, const std::string& op_name) const
{
    const Pointer<CellVariable<NDIM, double>> cc_var = var;
    return cc_var && op_name == s_op_name;
}

const std::string&
CartCellDoubleLinearCoarsen::getOperatorName() const
{
    return s_op_name;
}

int
CartCellDoubleLinearCoarsen::getOperatorPriority() const
{
    return COARSEN_OP_PRIORITY;
}

IntVector<NDIM>
CartCellDoubleLinearCoarsen::getStencilWidth() const
{
    return COARSEN_OP_STENCIL_WIDTH;
}

void
CartCellDoubleLinearCoarsen::coarsen(Patch<NDIM>& coarse,
                                     const Patch<NDIM>& fine,
                                     const int dst_component,
                                     const int src_component,
                                     const Box<NDIM>& coarse_box,
                                     const IntVector<NDIM>& ratio) const
{
    Pointer<CellData<NDIM, double>> cdata = coarse.getPatchData(dst_component);
    Pointer<CellData<NDIM, double>> fdata = fine.getPatchData(src_component);
#if !defined(NDEBUG)
    TBOX_ASSERT(cdata);
    TBOX_ASSERT(fdata);
    TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif

    const int data_depth = cdata->getDepth();
    const Box<NDIM>& cdata_box = cdata->getBox();
    const int cdata_gcw = cdata->getGhostCellWidth().max();
    const Box<NDIM>& fdata_box = fdata->getBox();
    const int fdata_gcw = fdata->getGhostCellWidth().max();
    const Box<NDIM>& patch_box_crse = coarse.getBox();
    const Box<NDIM>& patch_box_fine = fine.getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = coarse.getPatchGeometry();
    Array<BoundaryBox<NDIM>> bboxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(coarse);
    int touches_regular_bdry[2 * NDIM];
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        touches_regular_bdry[2 * axis + 0] = pgeom->getTouchesRegularBoundary(axis, 0);
        touches_regular_bdry[2 * axis + 1] = pgeom->getTouchesRegularBoundary(axis, 1);
    }

    for (int depth = 0; depth < data_depth; ++depth)
    {
        CC_LINEAR_COARSEN_FC(cdata->getPointer(depth),
                             cdata_gcw,
                             fdata->getPointer(depth),
                             fdata_gcw,
                             cdata_box.lower()(0),
                             cdata_box.upper()(0),
                             cdata_box.lower()(1),
                             cdata_box.upper()(1),
#if (NDIM == 3)
                             cdata_box.lower()(2),
                             cdata_box.upper()(2),
#endif
                             fdata_box.lower()(0),
                             fdata_box.upper()(0),
                             fdata_box.lower()(1),
                             fdata_box.upper()(1),
#if (NDIM == 3)
                             fdata_box.lower()(2),
                             fdata_box.upper()(2),
#endif
                             coarse_box.lower()(0),
                             coarse_box.upper()(0),
                             coarse_box.lower()(1),
                             coarse_box.upper()(1),
#if (NDIM == 3)
                             coarse_box.lower()(2),
                             coarse_box.upper()(2),
#endif
                             ratio);

        for (int k = 0; k < bboxes.getSize(); ++k)
        {
            const auto trimmed_bbox = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bboxes[k], coarse);
            const Box<NDIM>& cell_bdry_box = trimmed_bbox.getBox();
            const int bdry_normal_axis = static_cast<int>(trimmed_bbox.getLocationIndex() / 2);
            const int bdry_lower_side = (trimmed_bbox.getLocationIndex() % 2) == 0 ? 1 : 0;
            CC_LINEAR_COARSEN_BDRY_FC(cdata->getPointer(depth),
                                      cdata_gcw,
                                      fdata->getPointer(depth),
                                      fdata_gcw,
                                      cdata_box.lower()(0),
                                      cdata_box.upper()(0),
                                      cdata_box.lower()(1),
                                      cdata_box.upper()(1),
#if (NDIM == 3)
                                      cdata_box.lower()(2),
                                      cdata_box.upper()(2),
#endif
                                      fdata_box.lower()(0),
                                      fdata_box.upper()(0),
                                      fdata_box.lower()(1),
                                      fdata_box.upper()(1),
#if (NDIM == 3)
                                      fdata_box.lower()(2),
                                      fdata_box.upper()(2),
#endif
                                      cell_bdry_box.lower()(0),
                                      cell_bdry_box.upper()(0),
                                      cell_bdry_box.lower()(1),
                                      cell_bdry_box.upper()(1),
#if (NDIM == 3)
                                      cell_bdry_box.lower()(2),
                                      cell_bdry_box.upper()(2),
#endif
                                      ratio,
                                      touches_regular_bdry,
                                      patch_box_crse.lower()(0),
                                      patch_box_crse.upper()(0),
                                      patch_box_crse.lower()(1),
                                      patch_box_crse.upper()(1),
#if (NDIM == 3)
                                      patch_box_crse.lower()(2),
                                      patch_box_crse.upper()(2),
#endif
                                      bdry_normal_axis);
        }
    }
}
} // namespace IBTK
