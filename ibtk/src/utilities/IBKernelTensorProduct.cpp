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
canonicalize(std::string name)
{
    std::use_facet<std::ctype<char>>(std::locale::classic()).toupper(name.data(), name.data() + name.size());
    return name;
}
} // namespace

namespace IBTK
{
IBKernelTensorProduct::IBKernelTensorProduct(std::vector<IBKernel> factors) : d_factors(std::move(factors))
{
    if (d_factors.empty()) TBOX_ERROR("IBKernelTensorProduct requires at least one factor\n");
}

IBKernelTensorProduct::IBKernelTensorProduct(std::initializer_list<IBKernel> factors) : d_factors(factors)
{
    if (d_factors.empty()) TBOX_ERROR("IBKernelTensorProduct requires at least one factor\n");
}

IBKernelTensorProduct::IBKernelTensorProduct(const std::string& name)
{
    if (!isValidName(name)) TBOX_ERROR("Invalid IB kernel tensor-product name: " << name << '\n');
    const std::string canonical = canonicalize(name);
    const int orders = composite_orders(canonical);
    if (orders)
    {
        d_factors = { IBKernel("BSPLINE_" + std::to_string(orders / 10)),
                      IBKernel("BSPLINE_" + std::to_string(orders % 10)) };
    }
    else
    {
        d_factors = { IBKernel(canonical) };
    }
}

bool
IBKernelTensorProduct::isValidName(const std::string& name)
{
    const std::string canonical = canonicalize(name);
    return composite_orders(canonical) != 0 || IBKernel::isValidName(canonical);
}

std::size_t
IBKernelTensorProduct::size() const
{
    return d_factors.size();
}

const IBKernel&
IBKernelTensorProduct::operator[](std::size_t slot) const
{
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
