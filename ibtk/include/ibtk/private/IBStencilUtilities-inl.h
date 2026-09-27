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

#ifndef included_IBTK_IBStencilUtilities_inl
#define included_IBTK_IBStencilUtilities_inl
#include <ibtk/config.h>

#include <ibtk/private/IBStencilUtilities.h>

#include <cmath>
namespace IBTK::detail
{
inline int
ib_stencil_lower(const double x,
                 const double x_lower,
                 const double dx,
                 const int domain_lower,
                 const int cell,
                 const int width,
                 const double offset)
{
    if (width % 2 != 0)
    {
        const double q = (x - x_lower) / dx + static_cast<double>(domain_lower) - offset;
        return static_cast<int>(std::floor(q + 0.5)) - width / 2;
    }
    const double center = (static_cast<double>(cell - domain_lower) + 0.5) * dx + x_lower;
    return cell - width / 2 + (offset == 0.0 || x > center ? 1 : 0);
}
inline double
ib_stencil_coordinate(const double x,
                      const double x_lower,
                      const double dx,
                      const int domain_lower,
                      const int lower,
                      const double offset)
{
    const double first = (static_cast<double>(lower - domain_lower) + offset) * dx + x_lower;
    return (x - first) / dx;
}
} // namespace IBTK::detail
#endif
