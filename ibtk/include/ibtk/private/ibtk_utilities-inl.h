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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_ibtk_utilities_inl_h
#define included_IBTK_ibtk_utilities_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

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

inline std::string
format_samrai_output_filename(const int iteration_num, const std::string& data_dump_dirname, const std::string& prefix)
{
    std::ostringstream oss;
    oss << data_dump_dirname << "/" << prefix << ".";
    oss << std::setw(5) << std::setfill('0') << iteration_num << ".samrai." << std::setw(5) << IBTK_MPI::getRank();
    return oss.str();
}

inline std::string
format_iteration_output_filename(const int iteration_num,
                                 const std::string& data_dump_dirname,
                                 const std::string& prefix)
{
    std::ostringstream oss;
    oss << data_dump_dirname << "/" << prefix << ".";
    oss << std::setw(5) << std::setfill('0') << iteration_num;
    return oss.str();
}
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_ibtk_utilities_inl_h
