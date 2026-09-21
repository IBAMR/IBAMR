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

#ifndef included_IBTK_string_utilities
#define included_IBTK_string_utilities

#include <ibtk/config.h>

#include <string_view>

namespace IBTK
{
/*!
 * \brief Whether two strings are equal when the ASCII letters are compared without regard to case.
 *
 * The comparison does not depend on the locale, so it is suitable for the names of options, kernels and
 * solvers, which are ASCII.
 */
constexpr bool equals_ignore_case(std::string_view a, std::string_view b);
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/private/string_utilities-inl.h> // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_string_utilities
