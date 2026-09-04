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
#include <locale>
#include <optional>
#include <ostream>
#include <string>
#include <utility>

namespace
{
std::optional<std::array<unsigned int, 2>>
composite_bspline_orders(const std::string& name)
{
    static const std::pair<const char*, std::array<unsigned int, 2>> composites[] = {
        { "COMPOSITE_BSPLINE_12", { { 1, 2 } } }, { "COMPOSITE_BSPLINE_21", { { 2, 1 } } },
        { "COMPOSITE_BSPLINE_23", { { 2, 3 } } }, { "COMPOSITE_BSPLINE_32", { { 3, 2 } } },
        { "COMPOSITE_BSPLINE_34", { { 3, 4 } } }, { "COMPOSITE_BSPLINE_43", { { 4, 3 } } },
        { "COMPOSITE_BSPLINE_45", { { 4, 5 } } }, { "COMPOSITE_BSPLINE_54", { { 5, 4 } } },
        { "COMPOSITE_BSPLINE_56", { { 5, 6 } } }, { "COMPOSITE_BSPLINE_65", { { 6, 5 } } },
        { "DISCONTINUOUS_LINEAR", { { 2, 1 } } }
    };
    for (const auto& composite : composites)
        if (name == composite.first) return composite.second;
    return std::nullopt;
}

std::string
canonicalize_name(std::string name)
{
    std::use_facet<std::ctype<char>>(std::locale::classic()).toupper(name.data(), name.data() + name.size());
    return name;
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

IBKernelTensorProduct::IBKernelTensorProduct(const char* name)
    : IBKernelTensorProduct(name ? std::string(name) : std::string())
{
}

bool
IBKernelTensorProduct::isValidName(const std::string& name)
{
    const std::string canonical = canonicalize_name(name);
    return composite_bspline_orders(canonical).has_value() || IBKernel::isValidName(canonical);
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
    return { { { first, second, second } },
             static_cast<std::size_t>(first == second ? MIN_ACTIVE_FACTORS : MAX_ACTIVE_FACTORS) };
}

IBKernelTensorProduct::CanonicalFactors
IBKernelTensorProduct::parse_name(const std::string& name)
{
    if (!isValidName(name)) TBOX_ERROR("Invalid IB kernel tensor-product name: " << name << '\n');
    const std::string canonical = canonicalize_name(name);
    const auto orders = composite_bspline_orders(canonical);
    if (orders)
        return IBKernelTensorProduct::canonicalize({ IBKernel("BSPLINE_" + std::to_string((*orders)[0])),
                                                     IBKernel("BSPLINE_" + std::to_string((*orders)[1])) });
    return IBKernelTensorProduct::canonicalize({ IBKernel(canonical) });
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
