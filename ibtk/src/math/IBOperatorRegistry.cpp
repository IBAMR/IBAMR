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
const auto&
supplied_kernels()
{
    static const auto kernels =
        std::tuple_cat(make_bspline_kernels(std::make_index_sequence<IBTK_MAX_BSPLINE_ORDER>{}),
                       std::make_tuple(std::make_pair(IBTK::IBKernel::IB_3, IBTK::IBKernelEvaluatorIB3{}),
                                       std::make_pair(IBTK::IBKernel::IB_4, IBTK::IBKernelEvaluatorIB4{}),
                                       std::make_pair(IBTK::IBKernel::IB_5, IBTK::IBKernelEvaluatorIB5{}),
                                       std::make_pair(IBTK::IBKernel::IB_6, IBTK::IBKernelEvaluatorIB6{})));
    return kernels;
}
} // namespace

namespace IBTK
{
std::map<IBKernelTensorProduct, IBOperatorRegistry::Builder>&
IBOperatorRegistry::get_builders()
{
    static std::map<IBKernelTensorProduct, Builder> builders;
    return builders;
}

bool
IBOperatorRegistry::is_supplied_kernel(const IBKernelTensorProduct& kernel)
{
    const auto& kernels = supplied_kernels();
    for (std::size_t d = 0; d < kernel.size(); ++d)
    {
        const bool supplied =
            std::apply([&](const auto&... scalar) { return ((kernel[d] == scalar.first) || ...); }, kernels);
        if (!supplied)
        {
            return false;
        }
    }
    return true;
}

IBOperatorRegistry::Builder
IBOperatorRegistry::make_supplied_builder(const IBKernelTensorProduct& kernel)
{
    const auto& kernels = supplied_kernels();
    Builder builder;
    const auto select_normal = [&](const auto& normal)
    {
        if (normal.first != kernel[0])
        {
            return;
        }
        const auto select_tangential = [&](const auto& tangential)
        {
            if (tangential.first == kernel[kernel.size() - 1])
            {
                builder = make_builder(IBKernelTensorProductEvaluator{ normal.second, tangential.second });
            }
        };
        std::apply([&](const auto&... tangential) { (select_tangential(tangential), ...); }, kernels);
    };
    std::apply([&](const auto&... normal) { (select_normal(normal), ...); }, kernels);
    return builder;
}

void
IBOperatorRegistry::construct_interpolation_matrix_sc(Mat& mat,
                                                      const IBKernelTensorProduct& kernel,
                                                      Vec X_vec,
                                                      const std::vector<int>& num_dofs_per_proc,
                                                      int dof_index_idx,
                                                      Pointer<PatchLevel<NDIM>> patch_level)
{
    for (std::size_t d = 0; d < kernel.size(); ++d)
    {
        if (kernel[d] == IBKernel::UNKNOWN)
        {
            TBOX_ERROR("IBOperatorRegistry::construct_interpolation_matrix_sc(): unspecified kernel " << kernel
                                                                                                      << '\n');
        }
    }
    std::map<IBKernelTensorProduct, Builder>& builders = get_builders();
    auto builder = builders.find(kernel);
    if (builder == builders.end())
    {
        Builder supplied = make_supplied_builder(kernel);
        if (supplied)
        {
            builder = builders.emplace(kernel, std::move(supplied)).first;
        }
    }
    if (builder == builders.end())
    {
        TBOX_ERROR("IBOperatorRegistry::construct_interpolation_matrix_sc(): no registered evaluator for kernel "
                   << kernel << "\n");
    }
    builder->second(mat, X_vec, num_dofs_per_proc, dof_index_idx, patch_level);
}

} // namespace IBTK
