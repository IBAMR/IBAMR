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
#include <ibtk/LEInteractor.h>

#include <mpi.h>

#include <cstring>
#include <fstream>
#include <iostream>
#include <locale>
#include <sstream>
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
    std::use_facet<std::ctype<char>>(std::locale::classic()).tolower(lower.data(), lower.data() + lower.size());
    for (std::size_t i = 0; i < name.size(); i += 2) mixed[i] = lower[i];
    return { name, lower, mixed };
}

// Test-only inspection of the two-word layout, without adding a public key API.
std::array<std::uint64_t, 2>
numeric_key(const IBKernel& kernel)
{
    std::array<std::uint64_t, 2> words{};
    static_assert(sizeof(kernel) == sizeof(words), "identity must contain exactly two words");
    std::memcpy(words.data(), &kernel, sizeof(kernel));
    return words;
}
} // namespace

int
main(int argc, char* argv[])
{
    static_assert(std::is_trivially_copyable<IBKernel>::value, "identity copies must be trivial");
    static_assert(std::is_standard_layout<IBKernel>::value, "test inspects the two-word layout");
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int status = 0;
    try
    {
        // Ranks construct shared names in opposite orders with disjoint extras.
        // Compare the actual integers via typed MPI transport, not decoded names
        // or host-order bytes. The serial fixtures check the same value path.
        const char* common_names[] = { "IB_4", "APPLICATION_KERNEL", "ABCDEFGHIJKLMNOPQRSTUVWX" };
        std::array<std::uint64_t, 6> common_words{};
        for (int j = 0; j < 3; ++j)
        {
            const int i = rank % 2 ? 2 - j : j;
            const IBKernel extra((rank % 2 ? "ODD_EXTRA_" : "EVEN_EXTRA_") + std::to_string(j));
            const IBKernel common(common_names[i]);
            const auto words = numeric_key(common);
            common_words[2 * i] = words[0];
            common_words[2 * i + 1] = words[1];
        }
        auto root_words = common_words;
        MPI_Bcast(root_words.data(), 6, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        require(common_words == root_words, "rank order/subset changed numeric identity");

        // Independently tabulated positional values: A=1,...,Z=26,0=27,...,
        // 9=36,_=37, padded with zeros to two twelve-digit base-38 chunks.
        struct Encoding
        {
            const char* name;
            std::array<std::uint64_t, 2> words;
        };
        const Encoding encodings[] = {
            { "A", { { UINT64_C(238572050223552512), 0 } } },
            { "_", { { UINT64_C(8827165858271442944), 0 } } },
            { "ABCDEFGHIJKL", { { UINT64_C(251642104107238734), 0 } } },
            { "ABCDEFGHIJKLM", { { UINT64_C(251642104107238734), UINT64_C(3101436652906182656) } } },
            { "ABCDEFGHIJKLMNOPQRSTUVWX", { { UINT64_C(251642104107238734), UINT64_C(3191881425781291314) } } },
            { "ABCDEFGHIJKLMNOPQRSTUVW_", { { UINT64_C(251642104107238734), UINT64_C(3191881425781291327) } } },
            { "USER_DEFINED", { { UINT64_C(5130207666400688530), 0 } } },
            { "IB_4", { { UINT64_C(2165952653010967808), 0 } } }
        };
        for (const auto& expected : encodings)
        {
            const IBKernel kernel(expected.name), copy(kernel);
            require(numeric_key(kernel) == expected.words, "exact numeric encoding");
            require(copy == kernel && numeric_key(copy) == expected.words, "two-word value copy");
            require(copy.getName() == expected.name, "exact name round trip");
        }
        require(IBKernel("ABCDEFGHIJKL") != IBKernel("ABCDEFGHIJKLM"), "second chunk affects equality");
        require(IBKernel("ABCDEFGHIJKLMNOPQRSTUVWX") != IBKernel("ABCDEFGHIJKLMNOPQRSTUVW_"),
                "final character affects equality");
        require(IBKernel("A") != IBKernel("A_") && IBKernel("A") != IBKernel("AA"), "padding is not a name character");
        require_invalid([] { IBKernel kernel("ABCDEFGHIJKLMNOPQRSTUVWXY"); });
        for (const std::string invalid :
             { std::string("A B"), std::string("A-B"), std::string("A\0B", 3), std::string("A\x80", 2) })
            require_invalid([&] { IBKernel kernel(invalid); });

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
        for (std::size_t i = 0; i < sizeof(builtins) / sizeof(builtins[0]); ++i)
            require(builtins[i] == IBKernel(scalar_names[i]), "precomputed built-in key matches its public name");
        for (const char* name : scalar_names)
            for (const auto& spelling : spellings(name))
            {
                const IBKernel kernel(spelling);
                require(kernel.getName() == name, "scalar canonical spelling");
                require(kernel == IBKernel(name), "scalar equality");
                require(IBKernelTensorProduct::from_name(spelling) == IBKernelTensorProduct({ kernel }),
                        "isotropic legacy name");
                require(IBKernelTensorProduct::from_name(spelling).getNumberOfKernels() == 1,
                        "scalar name supplies one factor");
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
                require(product.getNumberOfKernels() == 2, "composite name supplies two factors");
                require(product.getKernel(0).getName() == expected.normal, "first factor");
                require(product.getKernel(1).getName() == expected.tangential, "second factor");
                require(product == IBKernelTensorProduct({ IBKernel(expected.normal), IBKernel(expected.tangential) }),
                        "ordered product");
                require(product != IBKernelTensorProduct({ IBKernel(expected.tangential), IBKernel(expected.normal) }),
                        "reversed product");
                require_invalid([&] { IBKernel scalar(spelling); });
            }

        std::string name = "ApplicationKernel";
        const IBKernel custom(name);
        name = "changed";
        require(custom.getName() == "APPLICATIONKERNEL", "scalar owns its value");
        for (const auto& spelling : spellings("APPLICATIONKERNEL"))
            require(custom == IBKernel(spelling), "custom ASCII uppercase canonicalization");
        require(custom != IBKernel::USER_DEFINED, "custom is not the legacy callback selector");
        require(IBKernel("IB_4_custom").getName() == "IB_4_CUSTOM", "custom built-in prefix");
        require(IBKernel("composite_bspline_custom").getName() == "COMPOSITE_BSPLINE_CUSTOM",
                "custom composite-like prefix");
        require(IBKernelTensorProduct::from_name(custom.getName()) == IBKernelTensorProduct({ custom }),
                "custom isotropic interpretation");
        IBKernel normal("ApplicationKernel"), tangential("AnotherKernel");
        const IBKernelTensorProduct product({ normal, tangential });
        normal = IBKernel("changed");
        tangential = IBKernel("changed");
        require(product.getKernel(0) == custom && product.getKernel(1) == IBKernel("AnotherKernel"),
                "product owns its factors");
        require(IBKernelTensorProduct({ custom }) == IBKernelTensorProduct({ IBKernel("applicationkernel") }),
                "custom product canonicalization");
        require(IBKernelTensorProduct({ IBKernel("IB_4"), IBKernel("IB_3") }) !=
                    IBKernelTensorProduct({ IBKernel("IB_3"), IBKernel("IB_4") }),
                "non-B-spline product");
        require(IBKernelTensorProduct({ IBKernel("piecewise_constant"), IBKernel("BSPLINE_1") }) ==
                    IBKernelTensorProduct({ IBKernel("BSPLINE_1"), IBKernel("BSPLINE_1") }),
                "normalized factor equality");

        require_invalid([] { IBKernel empty(""); });
        require_invalid([] { IBKernelTensorProduct::from_name(""); });
        require_invalid([] { IBKernelTensorProduct empty({}); });

        const std::vector<IBKernel> factors
        {
            IBKernel::IB_4, IBKernel::IB_3
#if (NDIM == 3)
                ,
                IBKernel::IB_6
#endif
        };
        const IBKernelTensorProduct ordered(factors), pair({ IBKernel::IB_4, IBKernel::IB_3 });
        const IBKernelTensorProduct single({ IBKernel::IB_4 }), repeated(std::vector<IBKernel>(NDIM, IBKernel::IB_4));
        const IBKernelTensorProduct triple({ IBKernel::IB_4, IBKernel::IB_3, IBKernel::IB_3 });
        const IBKernelTensorProduct four({ IBKernel::IB_4, IBKernel::IB_4, IBKernel::IB_4, IBKernel::IB_4 });
        require(single.getNumberOfKernels() == 1 && pair.getNumberOfKernels() == 2 &&
                    triple.getNumberOfKernels() == 3 && four.getNumberOfKernels() == 4,
                "supplied factor count");
        require(single != repeated && pair != triple, "length-sensitive sequence equality");
        require(single.isIsotropic() && repeated.isIsotropic() && four.isIsotropic() && !pair.isIsotropic(),
                "stored-factor isotropy");
        require(!IBTK::LEInteractor::isKnownKernel(four), "overlong isotropic description is unsupported");
        require(four.getKernel(3) == IBKernel::IB_4, "access beyond spatial dimension");
        bool bounds_rejected = false;
        try
        {
            four.getKernel(four.getNumberOfKernels());
        }
        catch (const std::out_of_range&)
        {
            bounds_rejected = true;
        }
        require(bounds_rejected, "stored-count bounds check");
        std::ostringstream diagnostic;
        diagnostic << single << ' ' << four;
        require(diagnostic.str() == "(IB_4) (IB_4,IB_4,IB_4,IB_4)", "diagnostic preserves exact sequence");
#if (NDIM == 3)
        require(ordered == IBKernelTensorProduct({ IBKernel::IB_4, IBKernel::IB_3, IBKernel::IB_6 }),
                "three-factor construction");
#else
        require(!IBTK::LEInteractor::isKnownKernel(triple), "three factors are unsupported in 2D");
#endif
        for (unsigned int slot = 0; slot < NDIM; ++slot)
            require(ordered.getKernel(slot) == factors[slot], "ordered slot access");
        // Rows are distinguished component axes; columns are Cartesian axes.
        const unsigned int expected_slots[NDIM][NDIM] = {
#if (NDIM == 2)
            { 0, 1 },
            { 1, 0 }
#else
            { 0, 1, 2 },
            { 1, 0, 2 },
            { 1, 2, 0 }
#endif
        };
        for (unsigned int component = 0; component < NDIM; ++component)
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                require(IBTK::LEInteractor::get_kernel_factor(ordered, axis, component) ==
                            factors[expected_slots[component][axis]],
                        "shared consumer directional mapping");
                require(IBTK::LEInteractor::get_kernel_factor(pair, axis, component) ==
                            IBKernel(component == axis ? IBKernel::IB_4 : IBKernel::IB_3),
                        "two-factor consumer mapping");
                require(IBTK::LEInteractor::get_kernel_factor(single, axis, component) == IBKernel::IB_4 &&
                            IBTK::LEInteractor::get_kernel_factor(repeated, axis, component) == IBKernel::IB_4,
                        "one-factor and repeated-factor consumer mappings agree");
#if (NDIM == 3)
                require(IBTK::LEInteractor::get_kernel_factor(pair, axis, component) ==
                            IBTK::LEInteractor::get_kernel_factor(triple, axis, component),
                        "two-factor and repeated-transverse consumer mappings agree");
#endif
            }

        if (rank == 0)
        {
            std::ofstream out("output");
            out << "Scalar names and aliases: PASS\n"
                << "Ordered and isotropic products: PASS\n"
                << "Custom canonical identities and two-word value copies: PASS\n"
                << "Exact encoding, 24-character round trip, and rank order/subset independence: PASS\n"
                << "Variable-length sequences and shared consumer expansion/mapping: PASS\n"
                << "Scalar/composite separation and invalid-name rejection: PASS\n";
            require(static_cast<bool>(out), "output write failed");
        }
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        status = 1;
    }
    int global_status = 0;
    MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Finalize();
    return global_status;
}
