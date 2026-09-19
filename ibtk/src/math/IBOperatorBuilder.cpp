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

#include <ibtk/IBOperatorBuilder.h>
#include <ibtk/ib_kernel_dispatch.h>

#include <tbox/Utilities.h>

#include <memory>
#include <type_traits>
#include <utility>

namespace IBTK
{
IBOperatorBuilder::IBOperatorBuilder(const IBKernelTensorProduct& kernel)
{
    const bool built_in = dispatch_ib_kernel_evaluator(
        kernel,
        [this](auto evaluator)
        { d_operations = std::make_shared<Model<std::remove_cvref_t<decltype(evaluator)>>>(std::move(evaluator)); });
    if (!built_in)
    {
        TBOX_ERROR("IBOperatorBuilder::IBOperatorBuilder():\n"
                   << "  no built-in evaluator for kernel " << kernel << ".\n"
                   << "  Built-in kernels are BSPLINE_1 through BSPLINE_" << MAX_BUILT_IN_BSPLINE_ORDER
                   << ", IB_3 through IB_6, and composite B-splines whose orders differ by one.\n"
                   << "  Construct the builder from an evaluator to use any other kernel.");
    }
}

bool
IBOperatorBuilder::is_built_in(const IBKernelTensorProduct& kernel)
{
    return dispatch_ib_kernel_evaluator(kernel, [](const auto&) {});
}
} // namespace IBTK
