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

#include <ibtk/CCScaledLaplaceOperator.h>
#include <ibtk/HierarchyMathOps.h>

#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

CCScaledLaplaceOperator::CCScaledLaplaceOperator(std::string object_name,
                                                 Pointer<Database> input_db,
                                                 const bool homogeneous_bc)
    : CCLaplaceOperator(std::move(object_name), input_db, homogeneous_bc)
{
    return;
} // CCScaledLaplaceOperator

void
CCScaledLaplaceOperator::setScalePatchDataIndex(const int scale_idx)
{
#if !defined(NDEBUG)
    if (scale_idx < 0)
    {
        TBOX_ERROR(d_object_name << "::setScalePatchDataIndex():\n"
                                 << "  scale_idx must be a valid patch data index.\n");
    }
#endif
    d_scale_idx = scale_idx;

    return;
} // setScalePatchDataIndex

/////////////////////////////// PROTECTED ////////////////////////////////////

void
CCScaledLaplaceOperator::applySpecialized(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    // CCLaplaceOperator::apply() has already filled all input ghost cells.
    for (int comp = 0; comp < x.getNumberOfComponents(); ++comp)
    {
        Pointer<CellVariable<NDIM, double>> x_var = x.getComponentVariable(comp);
        Pointer<CellVariable<NDIM, double>> y_var = y.getComponentVariable(comp);

        const int x_idx = x.getComponentDescriptorIndex(comp);
        const int y_idx = y.getComponentDescriptorIndex(comp);

        for (unsigned int depth = 0; depth < d_bc_coefs.size(); ++depth)
        {
            d_hier_math_ops->scaled_laplace(y_idx,
                                            y_var,
                                            d_poisson_spec,
                                            d_scale_idx,
                                            x_idx,
                                            x_var,
                                            nullptr,
                                            d_solution_time,
                                            0.0,
                                            invalid_index,
                                            nullptr,
                                            depth,
                                            depth);
        }
    }

    return;
} // applySpecialized

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK
