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

#include <mutex>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>

namespace
{
const std::array<std::string, 13>&
builtin_names()
{
    static const std::array<std::string, 13> names{ { "BSPLINE_1",
                                                      "BSPLINE_2",
                                                      "BSPLINE_3",
                                                      "BSPLINE_4",
                                                      "BSPLINE_5",
                                                      "BSPLINE_6",
                                                      "PIECEWISE_CUBIC",
                                                      "IB_3",
                                                      "IB_4",
                                                      "IB_4_W8",
                                                      "IB_5",
                                                      "IB_6",
                                                      "USER_DEFINED" } };
    return names;
}

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
IBKernel::IBKernel(Builtin kernel) : d_name(&builtin_names().at(static_cast<unsigned int>(kernel)))
{
}

IBKernel::IBKernel(const std::string& name)
{
    if (name.empty()) throw std::invalid_argument("IBKernel requires a nonempty scalar name");
    const std::string upper = uppercase_name(name);
    if (composite_orders(upper))
        throw std::invalid_argument("IBKernel requires a scalar name, not a composite kernel name");
    std::string canonical = name;
    if (upper == "PIECEWISE_CONSTANT")
        canonical = "BSPLINE_1";
    else if (upper == "PIECEWISE_LINEAR")
        canonical = "BSPLINE_2";
    for (const auto& builtin : builtin_names())
        if (upper == builtin || canonical == builtin)
        {
            d_name = &builtin;
            return;
        }
    // set insertion preserves existing elements and their addresses. This is
    // only name interning: it has no evaluator registration or dispatch role.
    static std::set<std::string> names;
    static std::mutex mutex;
    const std::lock_guard<std::mutex> lock(mutex);
    d_name = &*names.insert(canonical).first;
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
