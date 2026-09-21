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

#ifndef included_IBTK_string_utilities_inl_h
#define included_IBTK_string_utilities_inl_h

#include <ibtk/config.h>

#include <ibtk/string_utilities.h>

namespace IBTK
{
constexpr bool
equals_ignore_case(const std::string_view a, const std::string_view b)
{
    if (a.size() != b.size())
    {
        return false;
    }
    const auto lower = [](const char c) { return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c; };
    for (std::string_view::size_type i = 0; i < a.size(); ++i)
    {
        if (lower(a[i]) != lower(b[i]))
        {
            return false;
        }
    }
    return true;
}
} // namespace IBTK

#endif
