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

#ifndef included_IBTK_ib_kernel_dispatch
#define included_IBTK_ib_kernel_dispatch

#include <ibtk/config.h>

#include <ibtk/IBKernelTensorProduct.h>

#include <cstddef>

namespace IBTK
{
/*!
 * \brief Highest B-spline order that dispatch_ib_kernel_evaluator() selects by name.
 *
 * This bounds only the selection of kernels by name. IBKernelEvaluators::BSpline<N> can
 * be instantiated for any positive N, and an application can pass such an
 * evaluator, for example IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::BSpline<12>{} },
 * to the code that would otherwise select a built-in kernel by name.
 */
inline constexpr std::size_t MAX_BUILT_IN_BSPLINE_ORDER = 8;

/*!
 * \brief Call visitor with the evaluator for a built-in IB kernel.
 *
 * The built-in kernels are BSPLINE_1 through MAX_BUILT_IN_BSPLINE_ORDER, IB_3
 * through IB_6, and composite B-splines whose normal and transverse orders
 * differ by one (for example COMPOSITE_BSPLINE_2_3 and COMPOSITE_BSPLINE_3_2),
 * up to MAX_BUILT_IN_BSPLINE_ORDER. The visitor is called once, with an
 * IBKernelEvaluatorTensorProduct prvalue, so the kernel is a compile-time
 * type inside the visitor and its evaluation can be inlined.
 *
 * \return true if the visitor was called, and false if kernel is not built in.
 */
template <class Visitor>
bool dispatch_ib_kernel_evaluator(const IBKernelTensorProduct& kernel, Visitor&& visitor);
} // namespace IBTK

#include <ibtk/private/ib_kernel_dispatch-inl.h>
#endif
