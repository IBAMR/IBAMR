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

#include <ibtk/IBKernel.h>

#include <tbox/Utilities.h>

#include <algorithm>
#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace IBTK
{
bool
IBKernel::isValidName(const std::string& name)
{
    std::array<std::uint64_t, NAME_BLOCK_COUNT> encoded_name;
    return try_encode_name(name, encoded_name);
}

bool
IBKernel::try_encode_name(std::string_view name, std::array<std::uint64_t, NAME_BLOCK_COUNT>& encoded_name)
{
    if (name.empty() || name.size() > MAX_NAME_LENGTH) return false;

    // Kernel names use ASCII, so case conversion does not require a locale.
    std::array<char, MAX_NAME_LENGTH> canonical_name{};
    for (std::size_t i = 0; i < name.size(); ++i)
    {
        const char input = name[i];
        const char canonical = input >= 'a' && input <= 'z' ? static_cast<char>(input - 'a' + 'A') : input;
        if (!((canonical >= 'A' && canonical <= 'Z') || (canonical >= '0' && canonical <= '9') || canonical == '_'))
            return false;
        canonical_name[i] = canonical;
    }

    const std::string_view canonical(canonical_name.data(), name.size());
    if (canonical == "PIECEWISE_CONSTANT")
        encoded_name = encode_name("BSPLINE_1");
    else if (canonical == "PIECEWISE_LINEAR")
        encoded_name = encode_name("BSPLINE_2");
    else
    {
        encoded_name = {};
        for (std::size_t i = 0; i < canonical.size(); ++i)
        {
            const char c = canonical[i];
            const unsigned int digit = c == '_'             ? UNDERSCORE_DIGIT :
                                       c >= 'A' && c <= 'Z' ? c - 'A' + LETTER_DIGIT_OFFSET :
                                                              c - '0' + NUMBER_DIGIT_OFFSET;
            encoded_name[i / DIGITS_PER_BLOCK] = ENCODING_BASE * encoded_name[i / DIGITS_PER_BLOCK] + digit;
        }
        static constexpr auto encoding_base_powers = []()
        {
            std::array<std::uint64_t, DIGITS_PER_BLOCK + 1> powers{};
            powers[0] = 1;
            for (std::size_t i = 1; i < powers.size(); ++i) powers[i] = ENCODING_BASE * powers[i - 1];
            return powers;
        }();
        for (std::size_t block = 0; block < NAME_BLOCK_COUNT; ++block)
        {
            const std::size_t block_start = block * DIGITS_PER_BLOCK;
            const std::size_t active_digits =
                canonical.size() <= block_start ? 0 : std::min(DIGITS_PER_BLOCK, canonical.size() - block_start);
            encoded_name[block] *= encoding_base_powers[DIGITS_PER_BLOCK - active_digits];
        }
    }
    return true;
}

IBKernel::IBKernel(const std::string& name)
{
    if (!try_encode_name(name, d_name))
        TBOX_ERROR("IBKernel requires 1 to " << MAX_NAME_LENGTH << " ASCII letters, digits, or underscores: " << name
                                             << '\n');
}

IBKernel::IBKernel(const char* name)
{
    if (!name || !try_encode_name(name, d_name))
        TBOX_ERROR("IBKernel requires a nonnull name with 1 to " << MAX_NAME_LENGTH
                                                                 << " ASCII letters, digits, or underscores\n");
}

const std::vector<IBKernel>&
IBKernel::getStandardKernels()
{
    static const std::vector<IBKernel> kernels = { BSPLINE_1, BSPLINE_2, BSPLINE_3, BSPLINE_4,
                                                   BSPLINE_5, BSPLINE_6, IB_3,      IB_4,
                                                   IB_4_W8,   IB_5,      IB_6,      PIECEWISE_CUBIC };
    return kernels;
}

std::string
IBKernel::getName() const
{
    std::string name(NAME_BLOCK_COUNT * DIGITS_PER_BLOCK, '\0');
    for (std::size_t block = 0; block < NAME_BLOCK_COUNT; ++block)
    {
        std::uint64_t encoded_block = d_name[block];
        for (std::size_t i = DIGITS_PER_BLOCK; i > 0; --i)
        {
            const unsigned int digit = encoded_block % ENCODING_BASE;
            encoded_block /= ENCODING_BASE;
            name[DIGITS_PER_BLOCK * block + i - 1] = digit == 0                  ? '\0' :
                                                     digit < NUMBER_DIGIT_OFFSET ? 'A' + digit - LETTER_DIGIT_OFFSET :
                                                     digit < UNDERSCORE_DIGIT    ? '0' + digit - NUMBER_DIGIT_OFFSET :
                                                                                   '_';
        }
    }
    const auto padding = name.find('\0');
    if (padding != std::string::npos) name.resize(padding);
    return name;
}

} // namespace IBTK
