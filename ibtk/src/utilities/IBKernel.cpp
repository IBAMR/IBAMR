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
#include <ibtk/IBKernelTensorProduct.h>

#include <locale>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{

// Zero means this is not a recognized composite spelling.
int
composite_orders(const std::string& name)
{
    if (name == "DISCONTINUOUS_LINEAR") return 21;
    const int pairs[] = { 12, 21, 23, 32, 34, 43, 45, 54, 56, 65 };
    for (const int pair : pairs)
        if (name == "COMPOSITE_BSPLINE_" + std::to_string(pair)) return pair;
    return 0;
}
} // namespace

namespace IBTK
{
IBKernel::IBKernel(Builtin kernel)
{
    // Same order as Builtin. Constant initialization avoids repeated packing.
    static constexpr std::array<std::array<std::uint64_t, NAME_WORD_COUNT>, 13> builtin_keys{
        { encode_name("BSPLINE_1"),
          encode_name("BSPLINE_2"),
          encode_name("BSPLINE_3"),
          encode_name("BSPLINE_4"),
          encode_name("BSPLINE_5"),
          encode_name("BSPLINE_6"),
          encode_name("PIECEWISE_CUBIC"),
          encode_name("IB_3"),
          encode_name("IB_4"),
          encode_name("IB_4_W8"),
          encode_name("IB_5"),
          encode_name("IB_6"),
          encode_name("USER_DEFINED") }
    };
    d_name = builtin_keys.at(static_cast<unsigned int>(kernel));
}

IBKernel::IBKernel(const std::string& name)
{
    if (name.empty()) throw std::invalid_argument("IBKernel requires a nonempty scalar name");
    std::string canonical = name;
    std::use_facet<std::ctype<char>>(std::locale::classic())
        .toupper(canonical.data(), canonical.data() + canonical.size());
    if (composite_orders(canonical))
        throw std::invalid_argument("IBKernel requires a scalar name, not a composite kernel name");
    if (canonical == "PIECEWISE_CONSTANT")
        canonical = "BSPLINE_1";
    else if (canonical == "PIECEWISE_LINEAR")
        canonical = "BSPLINE_2";
    if (canonical.size() > MAX_NAME_LENGTH)
        throw std::invalid_argument("IBKernel scalar name exceeds " + std::to_string(MAX_NAME_LENGTH) + " characters");
    for (const char c : canonical)
        if (!((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_'))
            throw std::invalid_argument("IBKernel scalar name requires ASCII letters, digits, or underscore");
    d_name = encode_name(canonical);
}

std::string
IBKernel::getName() const
{
    std::string name(NAME_WORD_COUNT * DIGITS_PER_WORD, '\0');
    for (std::size_t chunk = 0; chunk < NAME_WORD_COUNT; ++chunk)
    {
        std::uint64_t word = d_name[chunk];
        for (std::size_t i = DIGITS_PER_WORD; i > 0; --i)
        {
            const unsigned int digit = word % ENCODING_BASE;
            word /= ENCODING_BASE;
            name[DIGITS_PER_WORD * chunk + i - 1] = digit == 0  ? '\0' :
                                                    digit <= 26 ? 'A' + digit - 1 :
                                                    digit <= 36 ? '0' + digit - 27 :
                                                                  '_';
        }
    }
    const auto padding = name.find('\0');
    if (padding != std::string::npos) name.resize(padding);
    return name;
}

IBKernelTensorProduct::IBKernelTensorProduct(std::vector<IBKernel> factors) : d_factors(std::move(factors))
{
    if (d_factors.empty()) throw std::invalid_argument("IBKernelTensorProduct requires at least one factor");
}

std::size_t
IBKernelTensorProduct::getNumberOfKernels() const
{
    return d_factors.size();
}

const IBKernel&
IBKernelTensorProduct::getKernel(std::size_t slot) const
{
    if (slot >= d_factors.size()) throw std::out_of_range("IB kernel slot out of range");
    return d_factors[slot];
}

bool
IBKernelTensorProduct::isIsotropic() const
{
    for (std::size_t d = 1; d < d_factors.size(); ++d)
        if (d_factors[d] != d_factors[0]) return false;
    return true;
}

bool
IBKernelTensorProduct::operator==(const IBKernelTensorProduct& other) const
{
    return d_factors == other.d_factors;
}

bool
IBKernelTensorProduct::operator!=(const IBKernelTensorProduct& other) const
{
    return !(*this == other);
}

IBKernelTensorProduct
IBKernelTensorProduct::from_name(const std::string& name)
{
    std::string canonical = name;
    std::use_facet<std::ctype<char>>(std::locale::classic())
        .toupper(canonical.data(), canonical.data() + canonical.size());
    const int orders = composite_orders(canonical);
    if (orders)
        return IBKernelTensorProduct(
            { IBKernel("BSPLINE_" + std::to_string(orders / 10)), IBKernel("BSPLINE_" + std::to_string(orders % 10)) });
    return IBKernelTensorProduct({ IBKernel(canonical) });
}

std::ostream&
operator<<(std::ostream& stream, const IBKernelTensorProduct& kernel)
{
    stream << '(';
    for (std::size_t d = 0; d < kernel.getNumberOfKernels(); ++d)
    {
        if (d) stream << ',';
        stream << kernel.getKernel(d).getName();
    }
    return stream << ')';
}
} // namespace IBTK
