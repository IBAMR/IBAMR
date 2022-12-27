// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2022 by the IBAMR developers
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

#include "ibtk/IndexUtilities.h"

#include "CellVariable.h"
#include "TotalAmountRefineAndCoarsen.h"

#include "ibtk/app_namespaces.h"

#include "ibtk/namespaces.h"

namespace IBTK
{
std::string TotalAmountRefine::s_object_name = "AMOUNT_CONSTANT_REFINE";
// The priority determines the order that refinement operations occur. We assign an arbitrary number here since we don't
// care about the order.
static const int REFINE_OP_PRIORITY = 0;

bool
TotalAmountRefine::findRefineOperator(const Pointer<hier::Variable<NDIM> >& var, const std::string& op_name) const
{
    // This operation is only valid on CellVariable's
    Pointer<CellVariable<NDIM, double> > cast_var = var;
    if (!cast_var.isNull() && op_name == s_object_name)
        return true;
    else
        return false;
} // TotalAmountRefine::findRefineOperator

const std::string&
TotalAmountRefine::getOperatorName() const
{
    return s_object_name;
} // TotalAmountRefine::getOperatorName

int
TotalAmountRefine::getOperatorPriority() const
{
    return REFINE_OP_PRIORITY;
} // TotalAmountRefine::getOperatorPriority

IntVector<NDIM>
TotalAmountRefine::getStencilWidth() const
{
    return IntVector<NDIM>(0);
} // TotalAmountRefine::getStencilWidth

void
TotalAmountRefine::refine(Patch<NDIM>& fine,
                          const Patch<NDIM>& coarse,
                          const int dst_component,
                          const int src_component,
                          const Box<NDIM>& fine_box,
                          const IntVector<NDIM>& ratio) const
{
    // Fill in fine_box with refined values from the coarse patch. We assume a constant profile across the coarse cell.
    // Therefore, the fine cells will be the amount from the coarse cell divided by the total number of fine cells.
    Pointer<CellData<NDIM, double> > fine_data = fine.getPatchData(dst_component);
    Pointer<CellData<NDIM, double> > coarse_data = coarse.getPatchData(src_component);
    int fine_cells_per_coarse = 1;
    for (int d = 0; d < NDIM; ++d) fine_cells_per_coarse *= ratio(d);
    // Loop over the fine box.
    for (CellIterator<NDIM> ci(fine_box); ci; ci++)
    {
        const CellIndex<NDIM>& fine_idx = ci();
        // Get the corresponding coarse index.
        const CellIndex<NDIM>& coarse_idx = IndexUtilities::coarsen(fine_idx, ratio);
        (*fine_data)(fine_idx) = (*coarse_data)(coarse_idx) / static_cast<double>(fine_cells_per_coarse);
    }
} // TotalAmountRefine::refine

std::string TotalAmountCoarsen::s_object_name = "AMOUNT_CONSTANT_COARSEN";
// The priority determines the order that coarsening operations occur. We assign an arbitrary number here since we don't
// care about the order.
static const int COARSEN_OP_PRIORITY = 0;

bool
TotalAmountCoarsen::findCoarsenOperator(const Pointer<hier::Variable<NDIM> >& var, const std::string& op_name) const
{
    // This operation is only valid on CellVariable's
    Pointer<CellVariable<NDIM, double> > cast_var = var;
    if (!cast_var.isNull() && op_name == s_object_name)
        return true;
    else
        return false;
} // TotalAmountCoarsen::findCoarsenOperator

const std::string&
TotalAmountCoarsen::getOperatorName() const
{
    return s_object_name;
} // TotalAmountCoarsen::getOperatorName

int
TotalAmountCoarsen::getOperatorPriority() const
{
    return COARSEN_OP_PRIORITY;
} // TotalAmountCoarsen::getOperatorPriority

IntVector<NDIM>
TotalAmountCoarsen::getStencilWidth() const
{
    return IntVector<NDIM>(0);
} // TotalAmountCoarsen::getStencilWidth

void
TotalAmountCoarsen::coarsen(Patch<NDIM>& coarse,
                            const Patch<NDIM>& fine,
                            const int dst_component,
                            const int src_component,
                            const Box<NDIM>& coarse_box,
                            const IntVector<NDIM>& ratio) const
{
    // Fill in fine_box with refined values from the coarse patch. The coarse data will simply be the sum of all refined
    // cells inside the coarse cell.
    Pointer<CellData<NDIM, double> > fine_data = fine.getPatchData(src_component);
    Pointer<CellData<NDIM, double> > coarse_data = coarse.getPatchData(dst_component);
    // Loop over the coarse box.
    for (CellIterator<NDIM> ci(coarse_box); ci; ci++)
    {
        const CellIndex<NDIM>& coarse_idx = ci();
        (*coarse_data)(coarse_idx) = 0.0;
        // Get the index on the fine box.
        const CellIndex<NDIM>& fine_idx = IndexUtilities::refine(coarse_idx, ratio);
        for (int i = 0; i < ratio(0); ++i)
        {
            for (int j = 0; j < ratio(1); ++j)
            {
#if (NDIM == 3)
                for (int k = 0; k < ratio(2); ++k)
                {
#endif
                    (*coarse_data)(coarse_idx) += (*fine_data)(fine_idx + IntVector<NDIM>(i,
                                                                                          j
#if (NDIM == 3)
                                                                                          ,
                                                                                          k
#endif
                                                                                          ));
#if (NDIM == 3)
                }
#endif
            }
        }
    }
} // TotalAmountCoarsen::coarsen
} // namespace IBTK
