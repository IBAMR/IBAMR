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

#include <stdexcept>
#include <string>

namespace
{
// Only built-in spellings are case-insensitive. Never modify a custom identity.
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
IBKernel::IBKernel(const std::string& name) : d_name(name)
{
    if (name.empty()) throw std::invalid_argument("IBKernel requires a nonempty scalar name");
    const std::string upper = uppercase_name(name);
    if (composite_orders(upper))
        throw std::invalid_argument("IBKernel requires a scalar name, not a composite kernel name");
    if (upper == "PIECEWISE_CONSTANT")
        d_name = "BSPLINE_1";
    else if (upper == "PIECEWISE_LINEAR")
        d_name = "BSPLINE_2";
    else
    {
        const char* scalar_names[] = {
            "BSPLINE_1", "BSPLINE_2", "BSPLINE_3", "BSPLINE_4", "BSPLINE_5", "BSPLINE_6",   "PIECEWISE_CUBIC",
            "IB_3",      "IB_4",      "IB_4_W8",   "IB_5",      "IB_6",      "USER_DEFINED"
        };
        for (const char* scalar_name : scalar_names)
            if (upper == scalar_name)
            {
                d_name = scalar_name;
                break;
            }
    }
}

const std::string&
IBKernel::getName() const
{
    return d_name;
}

bool
IBKernel::operator==(const IBKernel& other) const
{
    return d_name == other.d_name;
}

bool
IBKernel::operator!=(const IBKernel& other) const
{
    return !(*this == other);
}

IBKernelTensorProduct::IBKernelTensorProduct(const IBKernel& kernel) : IBKernelTensorProduct(kernel, kernel)
{
}

IBKernelTensorProduct::IBKernelTensorProduct(const IBKernel& normal, const IBKernel& tangential)
    : d_normal(normal), d_tangential(tangential)
{
}

const IBKernel&
IBKernelTensorProduct::getNormalKernel() const
{
    return d_normal;
}

const IBKernel&
IBKernelTensorProduct::getTangentialKernel() const
{
    return d_tangential;
}

bool
IBKernelTensorProduct::operator==(const IBKernelTensorProduct& other) const
{
    return d_normal == other.d_normal && d_tangential == other.d_tangential;
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
} // namespace IBTK
