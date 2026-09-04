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
#include <locale>
#include <ostream>
#include <string>
#include <utility>

namespace
{
int
composite_orders(const std::string& name)
{
    if (name == "DISCONTINUOUS_LINEAR") return 21;
    const int pairs[] = { 12, 21, 23, 32, 34, 43, 45, 54, 56, 65 };
    for (const int pair : pairs)
        if (name == "COMPOSITE_BSPLINE_" + std::to_string(pair)) return pair;
    return 0;
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
    return composite_orders(canonical) != 0 || IBKernel::isValidName(canonical);
}

IBKernelTensorProduct::IBKernelTensorProduct(CanonicalFactors factors)
    : d_factors(std::move(factors.factors)), d_size(factors.size)
{
}

IBKernelTensorProduct::CanonicalFactors
IBKernelTensorProduct::canonicalize(std::initializer_list<IBKernel> factors)
{
    if (factors.size() == 0 || factors.size() > 2) TBOX_ERROR("IBKernelTensorProduct requires one or two factors\n");
    const IBKernel first = *factors.begin();
    const IBKernel second = factors.size() == 2 ? *(factors.begin() + 1) : first;
    return { { { first, second, second } }, static_cast<std::size_t>(first == second ? 1 : 2) };
}

IBKernelTensorProduct::CanonicalFactors
IBKernelTensorProduct::parse_name(const std::string& name)
{
    if (!isValidName(name)) TBOX_ERROR("Invalid IB kernel tensor-product name: " << name << '\n');
    const std::string canonical = canonicalize_name(name);
    const int orders = composite_orders(canonical);
    if (orders)
        return IBKernelTensorProduct::canonicalize(
            { IBKernel("BSPLINE_" + std::to_string(orders / 10)), IBKernel("BSPLINE_" + std::to_string(orders % 10)) });
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
    return d_size == 1;
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
