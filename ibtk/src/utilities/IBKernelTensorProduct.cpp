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

#include <ibtk/IBKernelTensorProduct.h>

#include <tbox/Utilities.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

namespace
{
bool
equal_ignoring_case(std::string_view lhs, std::string_view rhs)
{
    // Kernel names use ASCII, so case conversion does not require a locale.
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        const char lhs_char = lhs[i] >= 'a' && lhs[i] <= 'z' ? static_cast<char>(lhs[i] - 'a' + 'A') : lhs[i];
        if (lhs_char != rhs[i]) return false;
    }
    return true;
}

std::optional<std::array<unsigned int, 2>>
composite_bspline_orders(std::string_view name)
{
    if (equal_ignoring_case(name, "DISCONTINUOUS_LINEAR")) return std::array<unsigned int, 2>{ { 2, 1 } };

    constexpr std::string_view prefix = "COMPOSITE_BSPLINE_";
    if (name.size() <= prefix.size() || !equal_ignoring_case(name.substr(0, prefix.size()), prefix))
        return std::nullopt;

    const std::string_view orders = name.substr(prefix.size());
    if (orders.size() == 2 && orders[0] >= '1' && orders[0] <= '9' && orders[1] >= '1' && orders[1] <= '9')
        return std::array<unsigned int, 2>{ { static_cast<unsigned int>(orders[0] - '0'),
                                              static_cast<unsigned int>(orders[1] - '0') } };

    const std::size_t separator = orders.find('_');
    if (separator == std::string_view::npos) return std::nullopt;
    unsigned int first = 0, second = 0;
    const char* const begin = orders.data();
    const char* const end = begin + orders.size();
    const auto first_result = std::from_chars(begin, begin + separator, first);
    const auto second_result = std::from_chars(begin + separator + 1, end, second);
    if (first_result.ec != std::errc{} || first_result.ptr != begin + separator || second_result.ec != std::errc{} ||
        second_result.ptr != end || first == 0 || second == 0)
        return std::nullopt;
    return std::array<unsigned int, 2>{ { first, second } };
}

IBTK::IBKernel bspline_kernel(unsigned int order);

std::optional<IBTK::IBKernel>
standard_scalar_kernel(std::string_view name)
{
    constexpr std::string_view bspline_prefix = "BSPLINE_";
    if (name.size() == bspline_prefix.size() + 1 &&
        equal_ignoring_case(name.substr(0, bspline_prefix.size()), bspline_prefix) && name.back() >= '1' &&
        name.back() <= '6')
        return bspline_kernel(name.back() - '0');

    if (name.size() == 4 && equal_ignoring_case(name.substr(0, 3), "IB_") && name.back() >= '3' && name.back() <= '6')
    {
        switch (name.back())
        {
        case '3':
            return IBTK::IBKernel::IB_3;
        case '4':
            return IBTK::IBKernel::IB_4;
        case '5':
            return IBTK::IBKernel::IB_5;
        case '6':
            return IBTK::IBKernel::IB_6;
        }
    }
    if (equal_ignoring_case(name, "IB_4_W8")) return IBTK::IBKernel::IB_4_W8;
    if (equal_ignoring_case(name, "PIECEWISE_CONSTANT")) return IBTK::IBKernel::BSPLINE_1;
    if (equal_ignoring_case(name, "PIECEWISE_LINEAR")) return IBTK::IBKernel::BSPLINE_2;
    if (equal_ignoring_case(name, "PIECEWISE_CUBIC")) return IBTK::IBKernel::PIECEWISE_CUBIC;
    return std::nullopt;
}

IBTK::IBKernel
bspline_kernel(unsigned int order)
{
    if (order == 0) TBOX_ERROR("B-spline order must be positive\n");
    switch (order)
    {
    case 1:
        return IBTK::IBKernel::BSPLINE_1;
    case 2:
        return IBTK::IBKernel::BSPLINE_2;
    case 3:
        return IBTK::IBKernel::BSPLINE_3;
    case 4:
        return IBTK::IBKernel::BSPLINE_4;
    case 5:
        return IBTK::IBKernel::BSPLINE_5;
    case 6:
        return IBTK::IBKernel::BSPLINE_6;
    default:
        return IBTK::IBKernel("BSPLINE_" + std::to_string(order));
    }
}

} // namespace

namespace IBTK
{
IBKernelTensorProduct::IBKernelTensorProduct(const IBKernel& factor) : IBKernelTensorProduct({ factor })
{
}

IBKernelTensorProduct::IBKernelTensorProduct(std::initializer_list<IBKernel> factors)
    : IBKernelTensorProduct(canonicalize(factors))
{
}

IBKernelTensorProduct::IBKernelTensorProduct(const std::string& name) : IBKernelTensorProduct(parse_name(name))
{
}

IBKernelTensorProduct::IBKernelTensorProduct(const char* name) : IBKernelTensorProduct(parse_name(name))
{
}

bool
IBKernelTensorProduct::isValidName(const std::string& name)
{
    return standard_scalar_kernel(name) || composite_bspline_orders(name).has_value() || IBKernel::isValidName(name);
}

IBKernelTensorProduct::IBKernelTensorProduct(CanonicalFactors factors)
    : d_factors(std::move(factors.factors)), d_size(factors.size)
{
}

IBKernelTensorProduct::CanonicalFactors
IBKernelTensorProduct::canonicalize(std::initializer_list<IBKernel> factors)
{
    if (factors.size() < MIN_ACTIVE_FACTORS || factors.size() > MAX_ACTIVE_FACTORS)
        TBOX_ERROR("IBKernelTensorProduct requires one or two factors\n");
    const IBKernel first = *factors.begin();
    const IBKernel second = factors.size() == MAX_ACTIVE_FACTORS ? *(factors.begin() + 1) : first;
    return { { { first, second } },
             static_cast<std::size_t>(first == second ? MIN_ACTIVE_FACTORS : MAX_ACTIVE_FACTORS) };
}

IBKernelTensorProduct::CanonicalFactors
IBKernelTensorProduct::parse_name(const std::string& name)
{
    if (const auto scalar = standard_scalar_kernel(name)) return IBKernelTensorProduct::canonicalize({ *scalar });
    const auto orders = composite_bspline_orders(name);
    if (orders)
        return IBKernelTensorProduct::canonicalize({ bspline_kernel((*orders)[0]), bspline_kernel((*orders)[1]) });
    return IBKernelTensorProduct::canonicalize({ IBKernel(name) });
}

IBKernelTensorProduct::CanonicalFactors
IBKernelTensorProduct::parse_name(const char* name)
{
    if (!name) TBOX_ERROR("Invalid null IB kernel tensor-product name\n");
    if (const auto scalar = standard_scalar_kernel(name)) return IBKernelTensorProduct::canonicalize({ *scalar });
    const auto orders = composite_bspline_orders(name);
    if (orders)
        return IBKernelTensorProduct::canonicalize({ bspline_kernel((*orders)[0]), bspline_kernel((*orders)[1]) });
    return IBKernelTensorProduct::canonicalize({ IBKernel(name) });
}

std::size_t
IBKernelTensorProduct::size() const
{
    return d_size;
}

const IBKernel&
IBKernelTensorProduct::operator[](std::size_t slot) const
{
    return d_factors[slot];
}

bool
IBKernelTensorProduct::isIsotropic() const
{
    return d_size == MIN_ACTIVE_FACTORS;
}

bool
IBKernelTensorProduct::operator==(const IBKernelTensorProduct& other) const
{
    return d_size == other.d_size && std::equal(d_factors.begin(), d_factors.begin() + d_size, other.d_factors.begin());
}

bool
IBKernelTensorProduct::operator!=(const IBKernelTensorProduct& other) const
{
    return !(*this == other);
}

bool
IBKernelTensorProduct::operator<(const IBKernelTensorProduct& other) const
{
    return std::lexicographical_compare(
        d_factors.begin(), d_factors.begin() + d_size, other.d_factors.begin(), other.d_factors.begin() + other.d_size);
}

bool
operator==(const IBKernelTensorProduct& product, const IBKernel& kernel)
{
    return product.isIsotropic() && product[0] == kernel;
}

bool
operator==(const IBKernel& kernel, const IBKernelTensorProduct& product)
{
    return product == kernel;
}

bool
operator!=(const IBKernelTensorProduct& product, const IBKernel& kernel)
{
    return !(product == kernel);
}

bool
operator!=(const IBKernel& kernel, const IBKernelTensorProduct& product)
{
    return !(kernel == product);
}

std::ostream&
operator<<(std::ostream& stream, const IBKernelTensorProduct& kernel)
{
    stream << '(';
    for (std::size_t d = 0; d < kernel.size(); ++d)
    {
        if (d) stream << ',';
        stream << kernel[d].getName();
    }
    return stream << ')';
}
} // namespace IBTK
