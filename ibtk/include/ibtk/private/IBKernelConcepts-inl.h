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

namespace IBTK::detail
{
template <class T, int Axis>
constexpr std::size_t
ib_kernel_stencil_size()
{
    constexpr std::array<std::size_t, NDIM> widths = T::template get_stencil_widths<Axis>();
    std::size_t count = 1;
    for (std::size_t width : widths)
    {
        count *= width;
    }
    return count;
}
} // namespace IBTK::detail
#endif
