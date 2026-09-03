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

#include <ostream>
#include <stdexcept>
#include <string>

namespace
{
// Each word holds twelve base-38 digits. Zero is reserved for right padding;
// 38^12 < 2^63, so neither chunk overflows. Call only with canonical names.
constexpr std::array<std::uint64_t, 2>
encode_name(const char* name, std::size_t length)
{
    std::array<std::uint64_t, 2> words{};
    for (std::size_t i = 0; i < 24; ++i)
    {
        const char c = i < length ? name[i] : '\0';
        const unsigned int digit = c == '\0' ? 0 : c == '_' ? 37 : c >= 'A' ? c - 'A' + 1 : c - '0' + 27;
        words[i / 12] = 38 * words[i / 12] + digit;
    }
    return words;
}

template <std::size_t N>
constexpr std::array<std::uint64_t, 2>
encode_literal(const char (&name)[N])
{
    return encode_name(name, N - 1);
}

// Same order as IBKernel::Builtin. Constant initialization avoids packing names
// whenever a built-in is constructed during numerical dispatch.
constexpr std::array<std::array<std::uint64_t, 2>, 13> builtin_keys{ { encode_literal("BSPLINE_1"),
                                                                       encode_literal("BSPLINE_2"),
                                                                       encode_literal("BSPLINE_3"),
                                                                       encode_literal("BSPLINE_4"),
                                                                       encode_literal("BSPLINE_5"),
                                                                       encode_literal("BSPLINE_6"),
                                                                       encode_literal("PIECEWISE_CUBIC"),
                                                                       encode_literal("IB_3"),
                                                                       encode_literal("IB_4"),
                                                                       encode_literal("IB_4_W8"),
                                                                       encode_literal("IB_5"),
                                                                       encode_literal("IB_6"),
                                                                       encode_literal("USER_DEFINED") } };

std::string
uppercase_name(std::string name)
{
    for (char& c : name)
        if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    return name;
}

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
IBKernel::IBKernel(Builtin kernel) : d_name(builtin_keys.at(static_cast<unsigned int>(kernel)))
{
}

IBKernel::IBKernel(const std::string& name)
{
    if (name.empty()) throw std::invalid_argument("IBKernel requires a nonempty scalar name");
    std::string canonical = uppercase_name(name);
    if (composite_orders(canonical))
        throw std::invalid_argument("IBKernel requires a scalar name, not a composite kernel name");
    if (canonical == "PIECEWISE_CONSTANT")
        canonical = "BSPLINE_1";
    else if (canonical == "PIECEWISE_LINEAR")
        canonical = "BSPLINE_2";
    if (canonical.size() > 24) throw std::invalid_argument("IBKernel scalar name exceeds 24 characters");
    for (const char c : canonical)
        if (!((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_'))
            throw std::invalid_argument("IBKernel scalar name requires ASCII letters, digits, or underscore");
    d_name = encode_name(canonical.data(), canonical.size());
}

std::string
IBKernel::getName() const
{
    std::string name(24, '\0');
    for (std::size_t chunk = 0; chunk < 2; ++chunk)
    {
        std::uint64_t word = d_name[chunk];
        for (std::size_t i = 12; i > 0; --i)
        {
            const unsigned int digit = word % 38;
            word /= 38;
            name[12 * chunk + i - 1] = digit == 0  ? '\0' :
                                       digit <= 26 ? 'A' + digit - 1 :
                                       digit <= 36 ? '0' + digit - 27 :
                                                     '_';
        }
    }
    const auto padding = name.find('\0');
    if (padding != std::string::npos) name.resize(padding);
    return name;
}

IBKernelTensorProduct::IBKernelTensorProduct(const IBKernel& kernel) : IBKernelTensorProduct(kernel, kernel)
{
}

IBKernelTensorProduct::IBKernelTensorProduct(const IBKernel& normal, const IBKernel& tangential)
    : d_factors{ { normal,
                   tangential
#if (NDIM == 3)
                   ,
                   tangential
#endif
      } },
      d_component_relative(normal != tangential)
{
}

IBKernelTensorProduct::IBKernelTensorProduct(const std::array<IBKernel, NDIM>& factors)
    : d_factors(factors), d_component_relative(false)
{
}

const IBKernel&
IBKernelTensorProduct::getKernel(unsigned int axis, unsigned int component_axis) const
{
    if (axis >= NDIM || component_axis >= NDIM) throw std::out_of_range("IB kernel direction out of range");
    return d_factors[d_component_relative ? (axis == component_axis ? 0 : 1) : axis];
}

bool
IBKernelTensorProduct::isIsotropic() const
{
    for (unsigned int d = 1; d < NDIM; ++d)
        if (d_factors[d] != d_factors[0]) return false;
    return true;
}

bool
IBKernelTensorProduct::isComponentRelative() const
{
    return d_component_relative;
}

const IBKernel&
IBKernelTensorProduct::getNormalKernel() const
{
    if (!d_component_relative && !isIsotropic()) throw std::logic_error("Cartesian kernel has no normal factor");
    return d_factors[0];
}

const IBKernel&
IBKernelTensorProduct::getTangentialKernel() const
{
    if (!d_component_relative && !isIsotropic()) throw std::logic_error("Cartesian kernel has no tangential factor");
    return d_factors[1];
}

bool
IBKernelTensorProduct::operator==(const IBKernelTensorProduct& other) const
{
    return d_component_relative == other.d_component_relative && d_factors == other.d_factors;
}

bool
IBKernelTensorProduct::operator!=(const IBKernelTensorProduct& other) const
{
    return !(*this == other);
}

IBKernelTensorProduct
IBKernelTensorProduct::from_name(const std::string& name)
{
    const int orders = composite_orders(uppercase_name(name));
    if (orders)
        return IBKernelTensorProduct(IBKernel("BSPLINE_" + std::to_string(orders / 10)),
                                     IBKernel("BSPLINE_" + std::to_string(orders % 10)));
    return IBKernelTensorProduct(IBKernel(name));
}

std::ostream&
operator<<(std::ostream& stream, const IBKernelTensorProduct& kernel)
{
    stream << (kernel.isComponentRelative() ? "component-relative(" : "Cartesian(");
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d) stream << ',';
        stream << kernel.getKernel(d, 0).getName();
    }
    return stream << ')';
}
} // namespace IBTK
