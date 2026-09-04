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
constexpr std::array<std::uint64_t, IBKernel::NAME_BLOCK_COUNT>
IBKernel::encode_name(std::string_view name)
{
    // Zero pads each block; the chosen radix and digit count fit in a signed 64-bit value.
    std::array<std::uint64_t, NAME_BLOCK_COUNT> blocks{};
    for (std::size_t i = 0; i < NAME_BLOCK_COUNT * DIGITS_PER_BLOCK; ++i)
    {
        const char c = i < name.size() ? name[i] : '\0';
        const unsigned int digit = c == '\0' ? 0 :
                                   c == '_'  ? UNDERSCORE_DIGIT :
                                   c >= 'A'  ? c - 'A' + LETTER_DIGIT_OFFSET :
                                               c - '0' + NUMBER_DIGIT_OFFSET;
        blocks[i / DIGITS_PER_BLOCK] = ENCODING_BASE * blocks[i / DIGITS_PER_BLOCK] + digit;
    }
    return blocks;
}

constexpr IBKernel
IBKernel::from_canonical_name(std::string_view name)
{
    return IBKernel(encode_name(name));
}

constexpr IBKernel::IBKernel(const std::array<std::uint64_t, NAME_BLOCK_COUNT>& name) : d_name(name)
{
}

inline constexpr IBKernel IBKernel::BSPLINE_1(IBKernel::from_canonical_name("BSPLINE_1"));
inline constexpr IBKernel IBKernel::BSPLINE_2(IBKernel::from_canonical_name("BSPLINE_2"));
inline constexpr IBKernel IBKernel::BSPLINE_3(IBKernel::from_canonical_name("BSPLINE_3"));
inline constexpr IBKernel IBKernel::BSPLINE_4(IBKernel::from_canonical_name("BSPLINE_4"));
inline constexpr IBKernel IBKernel::BSPLINE_5(IBKernel::from_canonical_name("BSPLINE_5"));
inline constexpr IBKernel IBKernel::BSPLINE_6(IBKernel::from_canonical_name("BSPLINE_6"));
inline constexpr IBKernel IBKernel::IB_3(IBKernel::from_canonical_name("IB_3"));
inline constexpr IBKernel IBKernel::IB_4(IBKernel::from_canonical_name("IB_4"));
inline constexpr IBKernel IBKernel::IB_4_W8(IBKernel::from_canonical_name("IB_4_W8"));
inline constexpr IBKernel IBKernel::IB_5(IBKernel::from_canonical_name("IB_5"));
inline constexpr IBKernel IBKernel::IB_6(IBKernel::from_canonical_name("IB_6"));
inline constexpr IBKernel IBKernel::PIECEWISE_CUBIC(IBKernel::from_canonical_name("PIECEWISE_CUBIC"));

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

inline bool
IBKernel::operator<(const IBKernel& other) const
{
    return d_name < other.d_name;
}
} // namespace IBTK
#endif
