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

#include "ibtk/CartCellDoubleCubicCoarsen.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianCellDoubleWeightedAverage.h"
#include "CellData.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <ostream>
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

bool
CartCellDoubleCubicCoarsen::findCoarsenOperator(const Pointer<Variable<NDIM> >& var, const std::string& op_name) const
{
    Pointer<CellVariable<NDIM, double> > cc_var = var;
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
    return d_weighted_average_coarsen_op.getStencilWidth();
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
        IBTK_DO_ONCE(TBOX_WARNING("CartCellDoubleCubicCoarsen::coarsen():\n"
                                  << "  cubic coarsening requires a refinement ratio of 4 or larger.\n"
                                  << "  reverting to weighted averaging." << std::endl););
        d_weighted_average_coarsen_op.coarsen(coarse, fine, dst_component, src_component, coarse_box, ratio);
        return;
    }
    Pointer<CellData<NDIM, double> > cdata = coarse.getPatchData(dst_component);
    Pointer<CellData<NDIM, double> > fdata = fine.getPatchData(src_component);
    const int U_fine_ghosts = (fdata->getGhostCellWidth()).max();
    const int U_crse_ghosts = (cdata->getGhostCellWidth()).max();
#if !defined(NDEBUG)
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
#endif
    const int data_depth = cdata->getDepth();
#if !defined(NDEBUG)
    TBOX_ASSERT(data_depth == fdata->getDepth());
#endif
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
    return;
} // coarsen

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
