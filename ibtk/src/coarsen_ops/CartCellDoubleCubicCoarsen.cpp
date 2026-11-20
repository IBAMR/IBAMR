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

#include <ibtk/CartCellDoubleCubicCoarsen.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <Box.h>
#include <CartesianCellDoubleWeightedAverage.h>
#include <CellData.h>
#include <CellIterator.h>

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
#define CC_CUBIC_COARSEN_FC IBTK_FC_FUNC(cccubiccoarsen2d, CCCUBICCOARSEN2D)
#endif
#if (NDIM == 3)
#define CC_CUBIC_COARSEN_FC IBTK_FC_FUNC(cccubiccoarsen3d, CCCUBICCOARSEN3D)
#endif

// Function interfaces
extern "C"
{
    void CC_CUBIC_COARSEN_FC(double* U_coarse,
                             const int& U_crse_gcw,
                             const double* U_fine0,
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

const std::string CartCellDoubleCubicCoarsen::s_op_name = "CUBIC_COARSEN";

namespace
{
static const int COARSEN_OP_PRIORITY = 0;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartCellDoubleCubicCoarsen::CartCellDoubleCubicCoarsen(const bool use_source_ghost_layer)
    : d_use_source_ghost_layer(use_source_ghost_layer)
{
    // intentionally blank
    return;
} // CartCellDoubleCubicCoarsen

void
CartCellDoubleCubicCoarsen::setFineLevel(Pointer<PatchLevel<NDIM>> fine_level)
{
    d_fine_level = fine_level;
    return;
} // setFineLevel

bool
CartCellDoubleCubicCoarsen::findCoarsenOperator(const Pointer<Variable<NDIM>>& var, const std::string& op_name) const
{
    Pointer<CellVariable<NDIM, double>> cc_var = var;
    return (cc_var && op_name == s_op_name);
} // findCoarsenOperator

const std::string&
CartCellDoubleCubicCoarsen::getOperatorName() const
{
    return s_op_name;
} // getOperatorName

int
CartCellDoubleCubicCoarsen::getOperatorPriority() const
{
    return COARSEN_OP_PRIORITY;
} // getOperatorPriority

IntVector<NDIM>
CartCellDoubleCubicCoarsen::getStencilWidth() const
{
    return d_use_source_ghost_layer ? IntVector<NDIM>(1) : d_weighted_average_coarsen_op.getStencilWidth();
} // getStencilWidth

void
CartCellDoubleCubicCoarsen::coarsen(Patch<NDIM>& coarse,
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
                TBOX_WARNING("CartCellDoubleCubicCoarsen::coarsen():\n"
                             << "  cubic coarsening requires refinement ratios of 4 or larger, or even ratios and a "
                                "source ghost layer\n"
                             << "  reverting to weighted averaging\n"););
            d_weighted_average_coarsen_op.coarsen(coarse, fine, dst_component, src_component, coarse_box, ratio);
            return;
        }
        if (!d_fine_level)
        {
            TBOX_ERROR("CartCellDoubleCubicCoarsen::coarsen():\n"
                       << "   the fine level is not set\n"
                       << "   call setFineLevel() before coarsening at a refinement ratio less than 4\n");
        }
    }
    Pointer<CellData<NDIM, double>> cdata = coarse.getPatchData(dst_component);
    Pointer<CellData<NDIM, double>> fdata = fine.getPatchData(src_component);
    const int U_fine_ghosts = (fdata->getGhostCellWidth()).max();
    const int U_crse_ghosts = (cdata->getGhostCellWidth()).max();
    if (U_fine_ghosts != (fdata->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartCellDoubleCubicCoarsen::coarsen():\n"
                   << "   fine patch data does not have uniform ghost cell widths" << std::endl);
    }
    if (U_crse_ghosts != (cdata->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartCellDoubleCubicCoarsen::coarsen():\n"
                   << "   coarse patch data does not have uniform ghost cell widths" << std::endl);
    }
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (ratio(d) % 2 == 1)
        {
            TBOX_ERROR("CartCellDoubleCubicCoarsen::coarsen():\n"
                       << "   refinement ratio between coarse and fine index spaces is odd" << std::endl);
        }
    }
    const int data_depth = cdata->getDepth();
    if (data_depth != fdata->getDepth())
    {
        TBOX_ERROR("CartCellDoubleCubicCoarsen::coarsen():\n"
                   << "   coarse and fine patch data have different depths\n");
    }
    const Box<NDIM>& patch_box_fine = fine.getBox();
    const Box<NDIM>& patch_box_crse = coarse.getBox();
    for (int depth = 0; depth < data_depth; ++depth)
    {
        double* const U_crse = cdata->getPointer(depth);
        const double* const U_fine = fdata->getPointer(depth);
        CC_CUBIC_COARSEN_FC(U_crse,
                            U_crse_ghosts,
                            U_fine,
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

    // Coarsen the cells whose stencil involves a fine cell outside the level (NaN in the source data) with stencils
    // that use only covered cells. A NaN in the source data at a covered cell propagates.
    const LevelCoverage coverage(d_fine_level);
    std::array<CubicCoarsenAxis, NDIM> axes;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        axes[d] = ratio(d) == 2 ? CubicCoarsenAxis::RATIO_2 : CubicCoarsenAxis::INTERIOR;
    }
    for (CellIterator<NDIM> ic(coarse_box); ic; ic++)
    {
        const Index<NDIM>& i_c = ic();
        if (!std::isnan((*cdata)(i_c, 0)))
        {
            continue;
        }
        Index<NDIM> base;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            base(d) = ratio(d) * i_c(d) + ratio(d) / 2 - 1;
        }
        // A coarse neighbor is covered if the fine cell of the cubic stencil in it is.
        const CubicCoarsenStencil stencil(axes,
                                          [&](const IntVector<NDIM>& offset)
                                          {
                                              Index<NDIM> fine_cell = base;
                                              for (unsigned int d = 0; d < NDIM; ++d)
                                              {
                                                  fine_cell(d) += offset(d) == 1 ? 2 : offset(d);
                                              }
                                              return coverage.covers(fine_cell);
                                          });
        for (int depth = 0; depth < data_depth; ++depth)
        {
            (*cdata)(i_c, depth) = stencil.apply(fdata->getArrayData(), base, depth);
        }
    }
    return;
} // coarsen

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
