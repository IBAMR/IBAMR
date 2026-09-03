// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_IBKernel_inl
#define included_IBTK_IBKernel_inl

namespace IBTK
{
constexpr std::array<std::uint64_t, IBKernel::NAME_WORD_COUNT>
IBKernel::encode_name(std::string_view name)
{
    // Zero pads each word; 38^12 < 2^63. Call only with canonical names.
    std::array<std::uint64_t, NAME_WORD_COUNT> words{};
    for (std::size_t i = 0; i < NAME_WORD_COUNT * DIGITS_PER_WORD; ++i)
    {
        const char c = i < name.size() ? name[i] : '\0';
        const unsigned int digit = c == '\0' ? 0 : c == '_' ? 37 : c >= 'A' ? c - 'A' + 1 : c - '0' + 27;
        words[i / DIGITS_PER_WORD] = ENCODING_BASE * words[i / DIGITS_PER_WORD] + digit;
    }
    return words;
}

inline bool
IBKernel::operator==(const IBKernel& other) const
{
    return d_name == other.d_name;
}

inline bool
IBKernel::operator!=(const IBKernel& other) const
{
    return !(*this == other);
}
} // namespace IBTK
#endif
