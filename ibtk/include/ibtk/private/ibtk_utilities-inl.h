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

    // Evaluate the smaller of H(phi) and 1 - H(phi) from the distance to the nearer
    // cutoff. This avoids cancellation in 0.5 + 0.5*phi/alpha near phi = -alpha;
    // the clamp keeps final rounding within [0, 1/2].
    const double z = (alpha - std::abs(phi)) / alpha;
    const double fraction = std::max(0.0, std::min(0.5, 0.5 * (z - std::sin(M_PI * z) / M_PI)));
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
