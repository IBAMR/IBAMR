// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_ib_kernel_dispatch_inl
#define included_IBTK_ib_kernel_dispatch_inl

#include <ibtk/config.h>

#include <ibtk/IBKernel.h>
#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/ib_kernel_dispatch.h>
#include <ibtk/ib_kernel_evaluators.h>

#include <string>
#include <utility>

namespace IBTK
{
namespace detail
{
inline IBKernel
ib_kernel_bspline(const std::size_t order)
{
    return IBKernel("BSPLINE_" + std::to_string(order));
}

template <class Visitor, class... Factors>
bool
call_if_ib_kernel_equals(const bool matches, Visitor& visitor, Factors... factors)
{
    if (!matches) return false;
    visitor(IBKernelEvaluatorTensorProduct{ std::move(factors)... });
    return true;
}

template <class Visitor, std::size_t... Index>
bool
dispatch_isotropic_ib_kernel_evaluator(const IBKernel& factor, Visitor& visitor, std::index_sequence<Index...>)
{
    return (call_if_ib_kernel_equals(
                factor == ib_kernel_bspline(Index + 1), visitor, IBKernelEvaluators::BSpline<Index + 1>{}) ||
            ...) ||
           call_if_ib_kernel_equals(factor == IBKernel::IB_3, visitor, IBKernelEvaluators::IB3{}) ||
           call_if_ib_kernel_equals(factor == IBKernel::IB_4, visitor, IBKernelEvaluators::IB4{}) ||
           call_if_ib_kernel_equals(factor == IBKernel::IB_5, visitor, IBKernelEvaluators::IB5{}) ||
           call_if_ib_kernel_equals(factor == IBKernel::IB_6, visitor, IBKernelEvaluators::IB6{});
}

template <class Visitor, std::size_t... Index>
bool
dispatch_composite_ib_kernel_evaluator(const IBKernel& normal,
                                       const IBKernel& transverse,
                                       Visitor& visitor,
                                       std::index_sequence<Index...>)
{
    return (
        (call_if_ib_kernel_equals(normal == ib_kernel_bspline(Index + 1) && transverse == ib_kernel_bspline(Index + 2),
                                  visitor,
                                  IBKernelEvaluators::BSpline<Index + 1>{},
                                  IBKernelEvaluators::BSpline<Index + 2>{}) ||
         call_if_ib_kernel_equals(normal == ib_kernel_bspline(Index + 2) && transverse == ib_kernel_bspline(Index + 1),
                                  visitor,
                                  IBKernelEvaluators::BSpline<Index + 2>{},
                                  IBKernelEvaluators::BSpline<Index + 1>{})) ||
        ...);
}
} // namespace detail

template <class Visitor>
bool
dispatch_ib_kernel_evaluator(const IBKernelTensorProduct& kernel, Visitor&& visitor)
{
    if (kernel.isIsotropic())
    {
        return detail::dispatch_isotropic_ib_kernel_evaluator(
            kernel[0], visitor, std::make_index_sequence<MAX_BUILT_IN_BSPLINE_ORDER>{});
    }
    return detail::dispatch_composite_ib_kernel_evaluator(
        kernel[0], kernel[1], visitor, std::make_index_sequence<MAX_BUILT_IN_BSPLINE_ORDER - 1>{});
}
} // namespace IBTK
#endif
