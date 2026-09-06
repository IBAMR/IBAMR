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

#include <tuple>
#include <utility>

#include <ibtk/namespaces.h>

namespace
{
constexpr int MAX_SUPPLIED_BSPLINE_ORDER = 6;

template <class Evaluator>
bool
is_within_bspline_limit(const Evaluator&, int)
{
    return true;
}

template <std::size_t N>
bool
is_within_bspline_limit(const IBTK::IBKernelEvaluatorBSpline<N>&, const int max_order)
{
    return N <= static_cast<std::size_t>(max_order);
}
} // namespace

namespace IBTK
{
int IBOperatorRegistry::s_max_bspline_order = MAX_SUPPLIED_BSPLINE_ORDER;
bool IBOperatorRegistry::s_builders_initialized = false;

void
IBOperatorRegistry::set_max_bspline_order(const int max_order)
{
    if (max_order < 1 || max_order > MAX_SUPPLIED_BSPLINE_ORDER)
    {
        TBOX_ERROR("IBOperatorRegistry::set_max_bspline_order(): order must be between 1 and "
                   << MAX_SUPPLIED_BSPLINE_ORDER << "\n");
    }
    if (s_builders_initialized)
    {
        TBOX_ERROR("IBOperatorRegistry::set_max_bspline_order(): registry is already initialized\n");
    }
    s_max_bspline_order = max_order;
}

std::map<IBKernelTensorProduct, IBOperatorRegistry::Builder>&
IBOperatorRegistry::get_builders()
{
    static auto builders = []
    {
        std::map<IBKernelTensorProduct, Builder> result;
        const auto kernels = std::make_tuple(std::make_pair(IBKernel::BSPLINE_1, IBKernelEvaluatorBSpline<1>{}),
                                             std::make_pair(IBKernel::BSPLINE_2, IBKernelEvaluatorBSpline<2>{}),
                                             std::make_pair(IBKernel::BSPLINE_3, IBKernelEvaluatorBSpline<3>{}),
                                             std::make_pair(IBKernel::BSPLINE_4, IBKernelEvaluatorBSpline<4>{}),
                                             std::make_pair(IBKernel::BSPLINE_5, IBKernelEvaluatorBSpline<5>{}),
                                             std::make_pair(IBKernel::BSPLINE_6, IBKernelEvaluatorBSpline<6>{}),
                                             std::make_pair(IBKernel::IB_3, IBKernelEvaluatorIB3{}),
                                             std::make_pair(IBKernel::IB_4, IBKernelEvaluatorIB4{}),
                                             std::make_pair(IBKernel::IB_5, IBKernelEvaluatorIB5{}),
                                             std::make_pair(IBKernel::IB_6, IBKernelEvaluatorIB6{}));
        const auto add_normal = [&](const auto& normal)
        {
            if (!is_within_bspline_limit(normal.second, s_max_bspline_order)) return;
            const auto add_tangential = [&](const auto& tangential)
            {
                if (is_within_bspline_limit(tangential.second, s_max_bspline_order))
                    result.emplace(IBKernelTensorProduct{ normal.first, tangential.first },
                                   make_builder(IBKernelTensorProductEvaluator{ normal.second, tangential.second }));
            };
            std::apply([&](const auto&... tangential) { (add_tangential(tangential), ...); }, kernels);
        };
        std::apply([&](const auto&... normal) { (add_normal(normal), ...); }, kernels);
        s_builders_initialized = true;
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
