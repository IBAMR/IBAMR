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
#include <ibtk/IBOperatorRegistry.h>

#include <string>
#include <tuple>
#include <utility>

#include <ibtk/namespaces.h>

namespace
{
template <std::size_t... I>
auto
make_bspline_kernels(std::index_sequence<I...>)
{
    return std::make_tuple(
        std::make_pair(IBTK::IBKernel("BSPLINE_" + std::to_string(I + 1)), IBTK::IBKernelEvaluatorBSpline<I + 1>{})...);
}
} // namespace

namespace IBTK
{
std::map<IBKernelTensorProduct, IBOperatorRegistry::Builder>&
IBOperatorRegistry::get_builders()
{
    static auto builders = []
    {
        std::map<IBKernelTensorProduct, Builder> result;
        const auto kernels = std::tuple_cat(make_bspline_kernels(std::make_index_sequence<IBTK_MAX_BSPLINE_ORDER>{}),
                                            std::make_tuple(std::make_pair(IBKernel::IB_3, IBKernelEvaluatorIB3{}),
                                                            std::make_pair(IBKernel::IB_4, IBKernelEvaluatorIB4{}),
                                                            std::make_pair(IBKernel::IB_5, IBKernelEvaluatorIB5{}),
                                                            std::make_pair(IBKernel::IB_6, IBKernelEvaluatorIB6{})));
        const auto add_normal = [&](const auto& normal)
        {
            const auto add_tangential = [&](const auto& tangential)
            {
                result.emplace(IBKernelTensorProduct{ normal.first, tangential.first },
                               make_builder(IBKernelTensorProductEvaluator{ normal.second, tangential.second }));
            };
            std::apply([&](const auto&... tangential) { (add_tangential(tangential), ...); }, kernels);
        };
        std::apply([&](const auto&... normal) { (add_normal(normal), ...); }, kernels);
        return result;
    }();
    return builders;
}

void
IBOperatorRegistry::construct_interpolation_matrix_sc(Mat& mat,
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
        TBOX_ERROR("IBOperatorRegistry::construct_interpolation_matrix_sc(): no registered evaluator for kernel "
                   << kernel << "\n");
    }
    builder->second(mat, X_vec, num_dofs_per_proc, dof_index_idx, patch_level);
}

} // namespace IBTK
