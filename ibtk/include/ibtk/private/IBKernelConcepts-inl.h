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
#ifndef included_IBTK_IBKernelConcepts_inl
#define included_IBTK_IBKernelConcepts_inl

#include <ibtk/config.h>

#include <ibtk/IBKernelConcepts.h>

namespace IBTK::detail
{
constexpr std::size_t
ib_kernel_width_product(const std::array<std::size_t, NDIM>& widths)
{
    std::size_t count = 1;
    for (const std::size_t width : widths)
    {
        count *= width;
    }
    return count;
}

template <class T, int Axis>
constexpr std::size_t
ib_kernel_stencil_size()
{
    return ib_kernel_width_product(T::template get_stencil_widths<Axis>());
}
} // namespace IBTK::detail
#endif
