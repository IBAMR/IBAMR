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

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using IBTK::IBKernel;
using IBTK::IBKernelTensorProduct;

namespace
{
void
require(bool valid, const std::string& message)
{
    if (!valid) throw std::runtime_error(message);
}

template <class F>
void
require_invalid(F operation)
{
    bool rejected = false;
    try
    {
        operation();
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "invalid scalar/composite description was accepted");
}

std::vector<std::string>
spellings(const std::string& name)
{
    std::string lower = name, mixed = name;
    for (std::size_t i = 0; i < name.size(); ++i)
        if (name[i] >= 'A' && name[i] <= 'Z')
        {
            lower[i] = static_cast<char>(name[i] - 'A' + 'a');
            if (i % 2 == 0) mixed[i] = lower[i];
        }
    return { name, lower, mixed };
}
} // namespace

int
main()
{
    static_assert(std::is_trivially_copyable<IBKernel>::value, "identity copies must be trivial");
    static_assert(sizeof(IBKernel) <= sizeof(void*), "identity must be compact");
    try
    {
        const IBKernel builtins[] = { IBKernel::BSPLINE_1,       IBKernel::BSPLINE_2, IBKernel::BSPLINE_3,
                                      IBKernel::BSPLINE_4,       IBKernel::BSPLINE_5, IBKernel::BSPLINE_6,
                                      IBKernel::PIECEWISE_CUBIC, IBKernel::IB_3,      IBKernel::IB_4,
                                      IBKernel::IB_4_W8,         IBKernel::IB_5,      IBKernel::IB_6,
                                      IBKernel::USER_DEFINED };
        for (const auto& kernel : builtins) require(kernel == IBKernel(kernel.getName()), "public built-in identity");
        const char* scalar_names[] = {
            "BSPLINE_1", "BSPLINE_2", "BSPLINE_3", "BSPLINE_4", "BSPLINE_5", "BSPLINE_6",   "PIECEWISE_CUBIC",
            "IB_3",      "IB_4",      "IB_4_W8",   "IB_5",      "IB_6",      "USER_DEFINED"
        };
        for (const char* name : scalar_names)
            for (const auto& spelling : spellings(name))
            {
                const IBKernel kernel(spelling);
                require(kernel.getName() == name, "scalar canonical spelling");
                require(kernel == IBKernel(name), "scalar equality");
                require(IBKernelTensorProduct::from_name(spelling) == IBKernelTensorProduct(kernel, kernel),
                        "isotropic legacy name");
                require(IBKernelTensorProduct(kernel) == IBKernelTensorProduct(kernel, kernel),
                        "isotropic constructor");
            }

        const std::pair<const char*, const char*> aliases[] = { { "PIECEWISE_CONSTANT", "BSPLINE_1" },
                                                                { "PIECEWISE_LINEAR", "BSPLINE_2" } };
        for (const auto& alias : aliases)
            for (const auto& spelling : spellings(alias.first))
            {
                require(IBKernel(spelling).getName() == alias.second, "scalar alias");
                require(IBKernelTensorProduct::from_name(spelling) == IBKernelTensorProduct::from_name(alias.second),
                        "isotropic alias");
            }

        struct Composite
        {
            const char* name;
            const char* normal;
            const char* tangential;
        };
        const Composite composites[] = {
            { "COMPOSITE_BSPLINE_12", "BSPLINE_1", "BSPLINE_2" }, { "COMPOSITE_BSPLINE_21", "BSPLINE_2", "BSPLINE_1" },
            { "COMPOSITE_BSPLINE_23", "BSPLINE_2", "BSPLINE_3" }, { "COMPOSITE_BSPLINE_32", "BSPLINE_3", "BSPLINE_2" },
            { "COMPOSITE_BSPLINE_34", "BSPLINE_3", "BSPLINE_4" }, { "COMPOSITE_BSPLINE_43", "BSPLINE_4", "BSPLINE_3" },
            { "COMPOSITE_BSPLINE_45", "BSPLINE_4", "BSPLINE_5" }, { "COMPOSITE_BSPLINE_54", "BSPLINE_5", "BSPLINE_4" },
            { "COMPOSITE_BSPLINE_56", "BSPLINE_5", "BSPLINE_6" }, { "COMPOSITE_BSPLINE_65", "BSPLINE_6", "BSPLINE_5" },
            { "DISCONTINUOUS_LINEAR", "BSPLINE_2", "BSPLINE_1" }
        };
        for (const auto& expected : composites)
            for (const auto& spelling : spellings(expected.name))
            {
                const auto product = IBKernelTensorProduct::from_name(spelling);
                require(product.getNormalKernel().getName() == expected.normal, "normal factor");
                require(product.getTangentialKernel().getName() == expected.tangential, "tangential factor");
                require(product == IBKernelTensorProduct(IBKernel(expected.normal), IBKernel(expected.tangential)),
                        "ordered product");
                require(product != IBKernelTensorProduct(IBKernel(expected.tangential), IBKernel(expected.normal)),
                        "reversed product");
                require_invalid([&] { IBKernel scalar(spelling); });
            }

        std::string name = "ApplicationKernel";
        const IBKernel custom(name);
        const std::string& retained_name = custom.getName();
        for (int i = 0; i < 1000; ++i) IBKernel("OtherApplicationKernel" + std::to_string(i));
        require(&retained_name == &IBKernel(name).getName() && custom == IBKernel(name), "interned identity lifetime");
        name = "changed";
        require(custom.getName() == "ApplicationKernel", "scalar owns its name");
        require(custom != IBKernel("applicationkernel"), "custom case sensitivity");
        require(IBKernel("IB_4_custom").getName() == "IB_4_custom", "custom built-in prefix");
        require(IBKernel("composite_bspline_custom").getName() == "composite_bspline_custom",
                "custom composite-like prefix");
        require(IBKernelTensorProduct::from_name(custom.getName()) == IBKernelTensorProduct(custom),
                "custom isotropic interpretation");
        IBKernel normal("ApplicationKernel"), tangential("AnotherKernel");
        const IBKernelTensorProduct product(normal, tangential);
        normal = IBKernel("changed");
        tangential = IBKernel("changed");
        require(product.getNormalKernel() == custom && product.getTangentialKernel() == IBKernel("AnotherKernel"),
                "product owns its factors");
        require(IBKernelTensorProduct(custom) != IBKernelTensorProduct(IBKernel("applicationkernel")),
                "custom product case sensitivity");
        require(IBKernelTensorProduct(IBKernel("IB_4"), IBKernel("IB_3")) !=
                    IBKernelTensorProduct(IBKernel("IB_3"), IBKernel("IB_4")),
                "non-B-spline product");
        require(IBKernelTensorProduct(IBKernel("piecewise_constant"), IBKernel("BSPLINE_1")) ==
                    IBKernelTensorProduct(IBKernel("BSPLINE_1")),
                "normalized factor equality");

        require_invalid([] { IBKernel empty(""); });
        require_invalid([] { IBKernelTensorProduct::from_name(""); });

        const std::array<IBKernel, NDIM> axes{ { IBKernel::IB_4,
                                                 IBKernel::IB_3
#if (NDIM == 3)
                                                 ,
                                                 IBKernel::IB_6
#endif
        } };
        const IBKernelTensorProduct cartesian(axes), relative(IBKernel::IB_4, IBKernel::IB_3);
        std::array<IBKernel, NDIM> equal_axes = axes;
        equal_axes.fill(IBKernel::IB_4);
        require(IBKernelTensorProduct(equal_axes) == IBKernelTensorProduct(IBKernel::IB_4), "equal Cartesian mapping");
        require(cartesian != relative, "Cartesian and relative mappings differ");
        for (unsigned int component = 0; component < NDIM; ++component)
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                require(cartesian.getKernel(axis, component) == axes[axis], "Cartesian x/y/z mapping");
                require(relative.getKernel(axis, component) ==
                            IBKernel(component == axis ? IBKernel::IB_4 : IBKernel::IB_3),
                        "component-relative mapping");
            }

        std::ofstream out("output");
        out << "Scalar names and aliases: PASS\n"
            << "Ordered and isotropic products: PASS\n"
            << "Custom identities and value ownership: PASS\n"
            << "Public compact identities and Cartesian/component direction mapping: PASS\n"
            << "Scalar/composite separation and empty-name rejection: PASS\n";
        if (!out) return 1;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return 1;
    }
    return 0;
}
