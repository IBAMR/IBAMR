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

#include <locale>
#include <string>
#include <vector>

namespace IBTK
{
bool
IBKernel::isValidName(const std::string& name)
{
    if (name.empty()) return false;
    std::string canonical = name;
    std::use_facet<std::ctype<char>>(std::locale::classic())
        .toupper(canonical.data(), canonical.data() + canonical.size());
    if (canonical == "PIECEWISE_CONSTANT")
        canonical = "BSPLINE_1";
    else if (canonical == "PIECEWISE_LINEAR")
        canonical = "BSPLINE_2";
    if (canonical.size() > MAX_NAME_LENGTH) return false;
    for (const char c : canonical)
        if (!((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_')) return false;
    return true;
}

IBKernel::IBKernel(const std::string& name)
{
    if (!isValidName(name))
        TBOX_ERROR("IBKernel requires 1 to " << MAX_NAME_LENGTH << " ASCII letters, digits, or underscores: " << name
                                             << '\n');
    std::string canonical = name;
    std::use_facet<std::ctype<char>>(std::locale::classic())
        .toupper(canonical.data(), canonical.data() + canonical.size());
    if (canonical == "PIECEWISE_CONSTANT")
        canonical = "BSPLINE_1";
    else if (canonical == "PIECEWISE_LINEAR")
        canonical = "BSPLINE_2";
    d_name = encode_name(canonical);
}

IBKernel::IBKernel(const char* name) : IBKernel(name ? std::string(name) : std::string())
{
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
            name[DIGITS_PER_BLOCK * block + i - 1] = digit == 0  ? '\0' :
                                                     digit <= 26 ? 'A' + digit - 1 :
                                                     digit <= 36 ? '0' + digit - 27 :
                                                                   '_';
        }
    }
    const auto padding = name.find('\0');
    if (padding != std::string::npos) name.resize(padding);
    return name;
}

} // namespace IBTK
