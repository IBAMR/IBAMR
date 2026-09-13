// ---------------------------------------------------------------------
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
// This file is part of IBAMR and is distributed under the 3-clause BSD
// license. The full text of the license can be found in COPYRIGHT.
// ---------------------------------------------------------------------

#include <ibtk/IBKernelTensorProductEvaluator.h>
#include <ibtk/SCInterpolationOperator.h>
#include <ibtk/kernels.h>

#include <tbox/Utilities.h>

#include <cstdlib>
#include <string>

#include <ibtk/namespaces.h>

namespace
{
template <std::size_t N = 1, class Select>
IBTK::SCInterpolationOperator
select_kernel(const IBTK::IBKernel& kernel, const Select& select)
{
    using namespace IBTK;
    if (kernel == IBKernel("BSPLINE_" + std::to_string(N)))
    {
        return select(Kernels::BSpline<N>{});
    }
    if constexpr (N < IBTK_MAX_BSPLINE_ORDER)
    {
        return select_kernel<N + 1>(kernel, select);
    }
    else
    {
        if (kernel == IBKernel::IB_3)
        {
            return select(Kernels::IB3{});
        }
        if (kernel == IBKernel::IB_4)
        {
            return select(Kernels::IB4{});
        }
        if (kernel == IBKernel::IB_5)
        {
            return select(Kernels::IB5{});
        }
        if (kernel == IBKernel::IB_6)
        {
            return select(Kernels::IB6{});
        }
        TBOX_ERROR("SCInterpolationOperator: no supplied evaluator for kernel " << kernel.getName() << '\n');
        std::abort();
    }
}
} // namespace

namespace IBTK
{
SCInterpolationOperator::SCInterpolationOperator(const IBKernelTensorProduct& kernel)
    : SCInterpolationOperator(select_kernel(
          kernel[0],
          [&](auto normal)
          {
              return select_kernel(
                  kernel[kernel.size() - 1],
                  [&](auto tangential) {
                      return SCInterpolationOperator(kernel, IBKernelTensorProductEvaluator{ normal, tangential });
                  });
          }))
{
}

const IBKernelTensorProduct&
SCInterpolationOperator::getKernel() const
{
    return d_kernel;
}

void
SCInterpolationOperator::constructMatrix(Mat& mat,
                                         Vec X,
                                         const std::vector<int>& num_dofs,
                                         int dof_index,
                                         Pointer<PatchLevel<NDIM>> level) const
{
    d_builder(mat, X, num_dofs, dof_index, level);
}
} // namespace IBTK
