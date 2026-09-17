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
IBKernel::is_valid_name(const std::string& name)
{
    std::array<std::uint64_t, NAME_BLOCK_COUNT> encoded_name;
    return try_encode_name(name, encoded_name);
}

bool
IBKernel::try_encode_name(std::string_view name, std::array<std::uint64_t, NAME_BLOCK_COUNT>& encoded_name)
{
    if (name.empty() || name.size() > MAX_NAME_LENGTH)
    {
        return false;
    }

    std::array<char, MAX_NAME_LENGTH> canonical_name{};
    for (std::size_t i = 0; i < name.size(); ++i)
    {
        const char input = name[i];
        const char canonical = input >= 'a' && input <= 'z' ? static_cast<char>(input - 'a' + 'A') : input;
        if (!((canonical >= 'A' && canonical <= 'Z') || (canonical >= '0' && canonical <= '9') || canonical == '_'))
        {
            return false;
        }
        canonical_name[i] = canonical;
    }

    const std::string_view canonical(canonical_name.data(), name.size());
    if (canonical == "PIECEWISE_CONSTANT")
    {
        encoded_name = encode_name("BSPLINE_1");
    }
    else if (canonical == "PIECEWISE_LINEAR")
    {
        encoded_name = encode_name("BSPLINE_2");
    }
    else
    {
        encoded_name = encode_name(canonical);
    }
    return true;
}

IBKernel::IBKernel(const std::string& name)
{
    if (!try_encode_name(name, d_name))
    {
        TBOX_ERROR("IBKernel requires 1 to " << MAX_NAME_LENGTH << " ASCII letters, digits, or underscores: " << name
                                             << '\n');
    }
}

IBKernel::IBKernel(const char* name)
{
    if (!name || !try_encode_name(name, d_name))
    {
        TBOX_ERROR("IBKernel requires a nonnull name with 1 to " << MAX_NAME_LENGTH
                                                                 << " ASCII letters, digits, or underscores\n");
    }
}

const std::vector<IBKernel>&
IBKernel::get_standard_kernels()
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
    const std::string::size_type padding = name.find('\0');
    if (padding != std::string::npos)
    {
        name.resize(padding);
    }
    return name;
}

} // namespace IBTK
