// ---------------------------------------------------------------------
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
// This file is part of IBAMR and is distributed under the 3-clause BSD
// license. The full text of the license can be found in COPYRIGHT.
// ---------------------------------------------------------------------

#ifndef included_IBAMR_IBImplicitStaggeredHierarchyIntegrator_inl
#define included_IBAMR_IBImplicitStaggeredHierarchyIntegrator_inl

#include <ibamr/config.h>

#include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>

#include <ibtk/PETScMatUtilities.h>

#include <tbox/Utilities.h>

#include <memory>
#include <utility>

namespace IBAMR
{
template <IBTK::IBKernelEvaluatorCartesian Evaluator>
inline void
IBImplicitStaggeredHierarchyIntegrator::setJacobianInterpolationKernel(Evaluator evaluator)
{
    if (d_integrator_is_initialized)
    {
        TBOX_ERROR("The interpolation kernel must be configured before initialization.\n");
    }
    d_interp_matrix_builder = make_matrix_builder(std::move(evaluator));
}

template <IBTK::IBKernelEvaluatorCartesian Evaluator>
inline IBImplicitStaggeredHierarchyIntegrator::InterpolationMatrixBuilder
IBImplicitStaggeredHierarchyIntegrator::make_matrix_builder(Evaluator evaluator)
{
    // std::function requires a copyable capture even when the evaluator is move-only.
    auto kernel = std::make_shared<const Evaluator>(std::move(evaluator));
    return [kernel](Mat& J,
                    Vec X,
                    const std::vector<int>& counts,
                    int dof,
                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level)
    { IBTK::PETScMatUtilities::constructPatchLevelSCInterpOp(J, *kernel, X, counts, dof, level); };
}
} // namespace IBAMR
#endif
