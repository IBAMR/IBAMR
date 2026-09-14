// ---------------------------------------------------------------------
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
// This file is part of IBAMR and is distributed under the 3-clause BSD
// license. The full text of the license can be found in COPYRIGHT.
// ---------------------------------------------------------------------
#ifndef included_IBTK_IBKernelConcepts_inl
#define included_IBTK_IBKernelConcepts_inl

#include <ibtk/config.h>

#include <ibtk/IBKernelConcepts.h>

#include <limits>
#include <utility>

namespace IBTK::detail
{
template <class T, int Axis>
constexpr bool
ib_kernel_cartesian_shape()
{
    using Weights =
        decltype(std::declval<const T&>().template evaluate<Axis>(std::declval<const std::array<double, NDIM>&>()));
    constexpr std::array<std::size_t, NDIM> widths = T::template get_stencil_widths<Axis>();
    std::size_t count = 1;
    for (std::size_t width : widths)
    {
        if (width == 0 || width > std::numeric_limits<std::size_t>::max() / count)
        {
            return false;
        }
        count *= width;
    }
    return count == IBKernelWeightsTraits<Weights>::extent;
}
} // namespace IBTK::detail
#endif
