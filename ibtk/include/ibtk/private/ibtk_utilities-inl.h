// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_ibtk_utilities_inl_h
#define included_IBTK_ibtk_utilities_inl_h

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>

#include <algorithm>
#include <cmath>

namespace IBTK
{
inline double
smooth_heaviside(const double& phi, const double& alpha)
{
    if (phi <= -alpha)
    {
        return 0.0;
    }
    if (phi >= alpha)
    {
        return 1.0;
    }
    if (phi == 0.0)
    {
        return 0.5;
    }

    const double z = (alpha - std::abs(phi)) / alpha;
    double fraction;
    if (z <= 0.5)
    {
        // For x = pi*z <= pi/2, use (x - sin(x))/(2*pi) through x^21.
        // The alternating-series remainder relative to the result is bounded
        // by 6*x^20/(23!*(1-x^2/20)) < 2.3e-18. Above this switch the direct
        // subtraction has condition number at most (pi+2)/(pi-2) < 4.51.
        const double x = M_PI * z;
        const double t = x * x;
        double series = -1.0 / 51090942171709440000.0;
        series = 1.0 / 121645100408832000.0 + t * series;
        series = -1.0 / 355687428096000.0 + t * series;
        series = 1.0 / 1307674368000.0 + t * series;
        series = -1.0 / 6227020800.0 + t * series;
        series = 1.0 / 39916800.0 + t * series;
        series = -1.0 / 362880.0 + t * series;
        series = 1.0 / 5040.0 + t * series;
        series = -1.0 / 120.0 + t * series;
        series = 1.0 / 6.0 + t * series;
        fraction = 0.5 * z * t * series;
    }
    else
    {
        fraction = 0.5 * (z - std::sin(M_PI * z) / M_PI);
    }
    // Guard the closed range against final rounding, preserving the small tail.
    fraction = std::max(0.0, std::min(0.5, fraction));
    return phi < 0.0 ? fraction : 1.0 - fraction;
}

inline double
smooth_delta(const double& phi, const double& alpha)
{
    if (std::abs(phi) >= alpha)
    {
        return 0.0;
    }
    const double z = (alpha - std::abs(phi)) / alpha;
    const double s = std::sin(0.5 * M_PI * z);
    return (s * s) / alpha;
}
} // namespace IBTK

#endif // #ifndef included_IBTK_ibtk_utilities_inl_h
