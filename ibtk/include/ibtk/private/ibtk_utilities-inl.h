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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_ibtk_utilities_inl_h
#define included_IBTK_ibtk_utilities_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <iomanip>
#include <sstream>
#include <string>

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBTK
{
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
