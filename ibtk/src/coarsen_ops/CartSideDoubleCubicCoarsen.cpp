// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <Box.h>
#include <SideData.h>
#include <SideIndex.h>
#include <SideIterator.h>
#include <SideVariable.h>

#include <array>
#include <cmath>
#include <ostream>
#include <string>

#include "./CubicCoarsenStencils.h"

#include <ibtk/namespaces.h> // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

// FORTRAN ROUTINES
#if (NDIM == 2)
#define SC_CUBIC_COARSEN_FC IBTK_FC_FUNC(sccubiccoarsen2d, SCCUBICCOARSEN2D)
#endif
#if (NDIM == 3)
#define SC_CUBIC_COARSEN_FC IBTK_FC_FUNC(sccubiccoarsen3d, SCCUBICCOARSEN3D)
#endif

// Function interfaces
extern "C"
{
    void SC_CUBIC_COARSEN_FC(double* U_coarse0,
                             double* U_coarse1,
#if (NDIM == 3)
                             double* U_coarse2,
#endif
                             const int& U_crse_gcw,
                             const double* U_fine0,
                             const double* U_fine1,
#if (NDIM == 3)
                             const double* U_fine2,
#endif
                             const int& U_fine_gcw,
                             const int& ilowerc0,
                             const int& iupperc0,
                             const int& ilowerc1,
                             const int& iupperc1,
#if (NDIM == 3)
                             const int& ilowerc2,
                             const int& iupperc2,
#endif
                             const int& ilowerf0,
                             const int& iupperf0,
                             const int& ilowerf1,
                             const int& iupperf1,
#if (NDIM == 3)
                             const int& ilowerf2,
                             const int& iupperf2,
#endif
                             const int* ratio_to_coarser,
                             const int* fblower,
                             const int* fbupper);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string CartSideDoubleCubicCoarsen::s_op_name = "CUBIC_COARSEN";

namespace
{
static const int COARSEN_OP_PRIORITY = 0;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideDoubleCubicCoarsen::CartSideDoubleCubicCoarsen(const bool use_source_ghost_layer)
    : d_use_source_ghost_layer(use_source_ghost_layer)
{
    // intentionally blank
    return;
} // CartSideDoubleCubicCoarsen

void
CartSideDoubleCubicCoarsen::setFineLevel(Pointer<PatchLevel<NDIM>> fine_level)
{
    d_fine_level = fine_level;
    return;
} // setFineLevel

bool
CartSideDoubleCubicCoarsen::findCoarsenOperator(const Pointer<Variable<NDIM>>& var, const std::string& op_name) const
{
    Pointer<SideVariable<NDIM, double>> sc_var = var;
    return (sc_var && op_name == s_op_name);
} // findCoarsenOperator

const std::string&
CartSideDoubleCubicCoarsen::getOperatorName() const
{
    return s_op_name;
} // getOperatorName

int
CartSideDoubleCubicCoarsen::getOperatorPriority() const
{
    return COARSEN_OP_PRIORITY;
} // getOperatorPriority

IntVector<NDIM>
CartSideDoubleCubicCoarsen::getStencilWidth() const
{
    return d_use_source_ghost_layer ? IntVector<NDIM>(1) : d_weighted_average_coarsen_op.getStencilWidth();
} // getStencilWidth

void
CartSideDoubleCubicCoarsen::coarsen(Patch<NDIM>& coarse,
                                    const Patch<NDIM>& fine,
                                    const int dst_component,
                                    const int src_component,
                                    const Box<NDIM>& coarse_box,
                                    const IntVector<NDIM>& ratio) const
{
    if (ratio.min() < 4)
    {
        bool even_ratios = true;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            even_ratios = even_ratios && ratio(d) % 2 == 0;
        }
        if (!even_ratios || !d_use_source_ghost_layer)
        {
            IBTK_DO_ONCE(
                TBOX_WARNING("CartSideDoubleCubicCoarsen::coarsen():\n"
                             << "  cubic coarsening requires refinement ratios of 4 or larger, or even ratios and a "
                                "source ghost layer\n"
                             << "  reverting to weighted averaging\n"););
            d_weighted_average_coarsen_op.coarsen(coarse, fine, dst_component, src_component, coarse_box, ratio);
            return;
        }
        if (!d_fine_level)
        {
            TBOX_ERROR("CartSideDoubleCubicCoarsen::coarsen():\n"
                       << "   the fine level is not set\n"
                       << "   call setFineLevel() before coarsening at a refinement ratio less than 4\n");
        }
    }
    Pointer<SideData<NDIM, double>> cdata = coarse.getPatchData(dst_component);
    Pointer<SideData<NDIM, double>> fdata = fine.getPatchData(src_component);
    const int U_fine_ghosts = (fdata->getGhostCellWidth()).max();
    const int U_crse_ghosts = (cdata->getGhostCellWidth()).max();
    if (U_fine_ghosts != (fdata->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartSideDoubleCubicCoarsen::coarsen():\n"
                   << "   fine patch data does not have uniform ghost cell widths" << std::endl);
    }
    if (U_crse_ghosts != (cdata->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartSideDoubleCubicCoarsen::coarsen():\n"
                   << "   coarse patch data does not have uniform ghost cell widths" << std::endl);
    }
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (ratio(d) % 2 == 1)
        {
            TBOX_ERROR("CartSideDoubleCubicCoarsen::coarsen():\n"
                       << "   refinement ratio between coarse and fine index spaces is odd" << std::endl);
        }
    }
    const int data_depth = cdata->getDepth();
    if (data_depth != fdata->getDepth())
    {
        TBOX_ERROR("CartSideDoubleCubicCoarsen::coarsen():\n"
                   << "   coarse and fine patch data have different depths\n");
    }
    const Box<NDIM>& patch_box_fine = fine.getBox();
    const Box<NDIM>& patch_box_crse = coarse.getBox();
    for (int depth = 0; depth < data_depth; ++depth)
    {
        double* const U_crse0 = cdata->getPointer(0, depth);
        double* const U_crse1 = cdata->getPointer(1, depth);
#if (NDIM == 3)
        double* const U_crse2 = cdata->getPointer(2, depth);
#endif
        const double* const U_fine0 = fdata->getPointer(0, depth);
        const double* const U_fine1 = fdata->getPointer(1, depth);
#if (NDIM == 3)
        const double* const U_fine2 = fdata->getPointer(2, depth);
#endif
        SC_CUBIC_COARSEN_FC(U_crse0,
                            U_crse1,
#if (NDIM == 3)
                            U_crse2,
#endif
                            U_crse_ghosts,
                            U_fine0,
                            U_fine1,
#if (NDIM == 3)
                            U_fine2,
#endif
                            U_fine_ghosts,
                            patch_box_crse.lower(0),
                            patch_box_crse.upper(0),
                            patch_box_crse.lower(1),
                            patch_box_crse.upper(1),
#if (NDIM == 3)
                            patch_box_crse.lower(2),
                            patch_box_crse.upper(2),
#endif
                            patch_box_fine.lower(0),
                            patch_box_fine.upper(0),
                            patch_box_fine.lower(1),
                            patch_box_fine.upper(1),
#if (NDIM == 3)
                            patch_box_fine.lower(2),
                            patch_box_fine.upper(2),
#endif
                            ratio,
                            coarse_box.lower(),
                            coarse_box.upper());
    }

    if (ratio.min() >= 4)
    {
        return;
    }

    // Coarsen the sides whose stencil involves a fine side outside the level (NaN in the source data) with stencils
    // that use only covered sides. A NaN in the source data at a covered side propagates.
    const LevelCoverage coverage(d_fine_level);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        std::array<CubicCoarsenAxis, NDIM> axes;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d == axis)
            {
                axes[d] = CubicCoarsenAxis::NORMAL;
            }
            else
            {
                axes[d] = ratio(d) == 2 ? CubicCoarsenAxis::RATIO_2 : CubicCoarsenAxis::INTERIOR;
            }
        }
        for (SideIterator<NDIM> is(coarse_box, axis); is; is++)
        {
            const SideIndex<NDIM>& i_c = is();
            if (!std::isnan((*cdata)(i_c, 0)))
            {
                continue;
            }
            Index<NDIM> base;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                base(d) = d == axis ? ratio(d) * i_c(d) : ratio(d) * i_c(d) + ratio(d) / 2 - 1;
            }
            // A coarse neighbor is covered if the level covers a cell on either side of the fine side of the cubic
            // stencil in it.
            const CubicCoarsenStencil stencil(axes,
                                              [&](const IntVector<NDIM>& offset)
                                              {
                                                  Index<NDIM> fine_cell = base;
                                                  for (unsigned int d = 0; d < NDIM; ++d)
                                                  {
                                                      fine_cell(d) += offset(d) == 1 ? 2 : offset(d);
                                                  }
                                                  if (coverage.covers(fine_cell))
                                                  {
                                                      return true;
                                                  }
                                                  fine_cell(axis) -= 1;
                                                  return coverage.covers(fine_cell);
                                              });
            for (int depth = 0; depth < data_depth; ++depth)
            {
                (*cdata)(i_c, depth) = stencil.apply(fdata->getArrayData(axis), base, depth);
            }
        }
    }
    return;
} // coarsen

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
