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

#include <ibtk/IBKernelEvaluators.h>
#include <ibtk/IBKernelTensorProductEvaluator.h>
#include <ibtk/SCInterpOpRegistry.h>

#include <tuple>
#include <utility>

#include <ibtk/namespaces.h>

namespace IBTK
{
std::map<IBKernelTensorProduct, SCInterpOpRegistry::Builder>&
SCInterpOpRegistry::get_builders()
{
    static auto builders = []
    {
        std::map<IBKernelTensorProduct, Builder> result;
        const auto kernels = std::make_tuple(std::make_pair(IBKernel::BSPLINE_1, IBKernelEvaluatorBSpline1{}),
                                             std::make_pair(IBKernel::BSPLINE_2, IBKernelEvaluatorBSpline2{}),
                                             std::make_pair(IBKernel::BSPLINE_3, IBKernelEvaluatorBSpline3{}),
                                             std::make_pair(IBKernel::BSPLINE_4, IBKernelEvaluatorBSpline4{}),
                                             std::make_pair(IBKernel::BSPLINE_5, IBKernelEvaluatorBSpline5{}),
                                             std::make_pair(IBKernel::BSPLINE_6, IBKernelEvaluatorBSpline6{}),
                                             std::make_pair(IBKernel::IB_3, IBKernelEvaluatorIB3{}),
                                             std::make_pair(IBKernel::IB_4, IBKernelEvaluatorIB4{}),
                                             std::make_pair(IBKernel::IB_5, IBKernelEvaluatorIB5{}),
                                             std::make_pair(IBKernel::IB_6, IBKernelEvaluatorIB6{}));
        const auto add_normal = [&](const auto& normal)
        {
            std::apply(
                [&](const auto&... tangential)
                {
                    (result.emplace(IBKernelTensorProduct{ normal.first, tangential.first },
                                    make_builder(IBKernelTensorProductEvaluator{ normal.second, tangential.second })),
                     ...);
                },
                kernels);
        };
        std::apply([&](const auto&... normal) { (add_normal(normal), ...); }, kernels);
        return result;
    }();
    return builders;
}

void
SCInterpOpRegistry::construct(Mat& mat,
                              const IBKernelTensorProduct& kernel,
                              Vec& X_vec,
                              const std::vector<int>& num_dofs_per_proc,
                              int dof_index_idx,
                              Pointer<PatchLevel<NDIM>> patch_level)
{
    const auto& builders = get_builders();
    const auto builder = builders.find(kernel);
    if (builder == builders.end())
    {
        TBOX_ERROR("SCInterpOpRegistry::construct(): no registered evaluator for kernel " << kernel << "\n");
    }
    builder->second(mat, X_vec, num_dofs_per_proc, dof_index_idx, patch_level);
}

} // namespace IBTK
