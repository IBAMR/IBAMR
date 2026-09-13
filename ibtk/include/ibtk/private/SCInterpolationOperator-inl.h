// ---------------------------------------------------------------------
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
// This file is part of IBAMR and is distributed under the 3-clause BSD
// license. The full text of the license can be found in COPYRIGHT.
// ---------------------------------------------------------------------

#ifndef included_IBTK_SCInterpolationOperator_inl
#define included_IBTK_SCInterpolationOperator_inl
#include <ibtk/config.h>

#include <ibtk/PETScMatUtilities.h>
#include <ibtk/SCInterpolationOperator.h>

#include <memory>
#include <utility>

namespace IBTK
{
template <TensorKernel Evaluator>
SCInterpolationOperator::SCInterpolationOperator(IBKernelTensorProduct kernel, Evaluator evaluator)
    : d_kernel(std::move(kernel)),
      d_builder([owned = std::make_shared<const Evaluator>(std::move(evaluator))](
                    Mat& mat,
                    Vec X,
                    const std::vector<int>& num_dofs,
                    int dof_index,
                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level)
                { PETScMatUtilities::constructPatchLevelSCInterpOp(mat, *owned, X, num_dofs, dof_index, level); })
{
}
} // namespace IBTK
#endif
