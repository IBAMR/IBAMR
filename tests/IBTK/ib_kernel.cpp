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
#include <ibtk/IBTKInit.h>

#include <tbox/Utilities.h>

#include <mpi.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iterator>
#include <locale>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using IBTK::IBKernel;
using IBTK::IBKernelTensorProduct;

bool ib_kernel_static_initialization_valid();

namespace
{
constexpr std::size_t ENCODED_BLOCK_COUNT = 2;

void
require(bool valid, const std::string& message)
{
    if (!valid) TBOX_ERROR(message << '\n');
}

std::vector<std::string>
spellings(const std::string& name)
{
    std::string lower = name, mixed = name;
    std::use_facet<std::ctype<char>>(std::locale::classic()).tolower(lower.data(), lower.data() + lower.size());
    for (std::size_t i = 0; i < name.size(); i += 2) mixed[i] = lower[i];
    return { name, lower, mixed };
}

// Test-only inspection of the compact layout, without adding a public key API.
std::array<std::uint64_t, ENCODED_BLOCK_COUNT>
numeric_blocks(const IBKernel& kernel)
{
    std::array<std::uint64_t, ENCODED_BLOCK_COUNT> blocks{};
    static_assert(sizeof(kernel) == sizeof(blocks), "identity must contain the expected number of blocks");
    std::memcpy(blocks.data(), &kernel, sizeof(kernel));
    return blocks;
}

IBKernelTensorProduct
copy_product(IBKernelTensorProduct kernel)
{
    return kernel;
}
} // namespace

int
main(int argc, char* argv[])
{
    static_assert(std::is_trivially_copyable<IBKernel>::value, "identity copies must be trivial");
    static_assert(std::is_standard_layout<IBKernel>::value, "test inspects the two-block layout");
    static_assert(std::is_convertible<std::string, IBKernelTensorProduct>::value,
                  "legacy strings must convert to tensor products");
    static_assert(std::is_convertible<const char*, IBKernelTensorProduct>::value,
                  "legacy C strings must convert to tensor products");
    static_assert(std::is_convertible<IBKernel, IBKernelTensorProduct>::value,
                  "scalar kernels must convert to tensor products");
    static_assert(!std::is_default_constructible<IBKernelTensorProduct>::value,
                  "tensor products have no unspecified state");
    static_assert(!std::is_constructible<IBKernelTensorProduct, std::vector<IBKernel>>::value,
                  "dynamic factor containers are not a public construction surface");
    static_assert(!std::is_constructible<IBKernelTensorProduct, std::array<IBKernel, NDIM>>::value,
                  "exact arrays are not a public construction surface");
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const std::string input_file = argc > 1 ? argv[1] : "";
    if (input_file.find("invalid_scalar") != std::string::npos)
    {
        std::ofstream("output") << "IBKernel rejects malformed scalar names\n";
        const IBKernel invalid("BAD NAME");
        TBOX_ERROR("IBKernel invalid-scalar test continued unexpectedly: " << invalid.getName() << '\n');
    }
    if (input_file.find("empty_product") != std::string::npos)
    {
        std::ofstream("output") << "IBKernelTensorProduct rejects empty factor sequences\n";
        const IBKernelTensorProduct invalid(std::initializer_list<IBKernel>{});
        TBOX_ERROR("IBKernelTensorProduct empty-factor test continued unexpectedly: " << invalid << '\n');
    }
    if (input_file.find("null_product") != std::string::npos)
    {
        std::ofstream("output") << "IBKernelTensorProduct rejects null names\n";
        const IBKernelTensorProduct invalid(static_cast<const char*>(nullptr));
        TBOX_ERROR("IBKernelTensorProduct null-name test continued unexpectedly: " << invalid << '\n');
    }
    if (input_file.find("three_factor") != std::string::npos)
    {
        std::ofstream("output") << "IBKernelTensorProduct rejects three factors\n";
        const IBKernelTensorProduct invalid({ IBKernel::IB_4, IBKernel::IB_4, IBKernel::IB_4 });
        TBOX_ERROR("IBKernelTensorProduct three-factor test continued unexpectedly: " << invalid << '\n');
    }

    // Ranks construct shared names in opposite orders with disjoint extras.
    // Compare the actual integers via typed MPI transport, not decoded names
    // or host-order bytes. The serial fixtures check the same value path.
    constexpr const char* COMMON_NAMES[] = { "IB_4", "APPLICATION_KERNEL", "ABCDEFGHIJKLMNOPQRSTUVWX" };
    constexpr std::size_t COMMON_NAME_COUNT = std::size(COMMON_NAMES);
    std::array<std::uint64_t, ENCODED_BLOCK_COUNT * COMMON_NAME_COUNT> common_blocks{};
    for (std::size_t j = 0; j < COMMON_NAME_COUNT; ++j)
    {
        const std::size_t i = rank % 2 ? COMMON_NAME_COUNT - 1 - j : j;
        const IBKernel extra((rank % 2 ? "ODD_EXTRA_" : "EVEN_EXTRA_") + std::to_string(j));
        const IBKernel common(COMMON_NAMES[i]);
        const auto blocks = numeric_blocks(common);
        for (std::size_t block = 0; block < ENCODED_BLOCK_COUNT; ++block)
            common_blocks[ENCODED_BLOCK_COUNT * i + block] = blocks[block];
    }
    auto root_blocks = common_blocks;
    MPI_Bcast(root_blocks.data(), static_cast<int>(root_blocks.size()), MPI_UINT64_T, 0, MPI_COMM_WORLD);
    require(common_blocks == root_blocks, "rank order/subset changed numeric identity");

    // Independently tabulated values for the canonical letter, digit, and
    // underscore encoding with fixed-width block padding.
    struct Encoding
    {
        const char* name;
        std::array<std::uint64_t, ENCODED_BLOCK_COUNT> blocks;
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
        require(numeric_blocks(kernel) == expected.blocks, "exact numeric encoding");
        require(copy == kernel && numeric_blocks(copy) == expected.blocks, "two-block value copy");
        require(copy.getName() == expected.name, "exact name round trip");
    }
    require(IBKernel("ABCDEFGHIJKL") != IBKernel("ABCDEFGHIJKLM"), "second chunk affects equality");
    require(IBKernel("ABCDEFGHIJKLMNOPQRSTUVWX") != IBKernel("ABCDEFGHIJKLMNOPQRSTUVW_"),
            "final character affects equality");
    require(IBKernel("A") != IBKernel("A_") && IBKernel("A") != IBKernel("AA"), "padding is not a name character");
    const char* const scalar_c_string = "piecewise_constant";
    require(IBKernel{ scalar_c_string } == IBKernel(IBKernel::BSPLINE_1), "scalar C-string construction");
    require(!IBKernel::isValidName("ABCDEFGHIJKLMNOPQRSTUVWXY"), "over-capacity scalar name accepted");
    for (const std::string invalid :
         { std::string("A B"), std::string("A-B"), std::string("A\0B", 3), std::string("A\x80", 2) })
        require(!IBKernel::isValidName(invalid), "malformed scalar name accepted");

    const IBKernel standard_values[] = { IBKernel::BSPLINE_1, IBKernel::BSPLINE_2, IBKernel::BSPLINE_3,
                                         IBKernel::BSPLINE_4, IBKernel::BSPLINE_5, IBKernel::BSPLINE_6,
                                         IBKernel::IB_3,      IBKernel::IB_4,      IBKernel::IB_4_W8,
                                         IBKernel::IB_5,      IBKernel::IB_6,      IBKernel::PIECEWISE_CUBIC };
    const char* scalar_names[] = { "BSPLINE_1", "BSPLINE_2", "BSPLINE_3", "BSPLINE_4", "BSPLINE_5", "BSPLINE_6",
                                   "IB_3",      "IB_4",      "IB_4_W8",   "IB_5",      "IB_6",      "PIECEWISE_CUBIC" };
    const auto& catalog = IBKernel::getStandardKernels();
    const std::size_t n_standard_values = sizeof(standard_values) / sizeof(standard_values[0]);
    require(catalog.size() == n_standard_values, "standard catalog size");
    require(std::is_sorted(catalog.begin(), catalog.end()), "standard catalog order");
    require(std::set<IBKernel>(catalog.begin(), catalog.end()).size() == catalog.size(),
            "standard kernels are not unique container keys");
    for (std::size_t i = 0; i < n_standard_values; ++i)
    {
        require(standard_values[i] == IBKernel(scalar_names[i]), "constant-initialized standard value");
        require(catalog[i] == standard_values[i], "standard catalog value and order");
        require(catalog[i].getName() != "USER_DEFINED" && catalog[i].getName() != "PIECEWISE_CONSTANT" &&
                    catalog[i].getName() != "PIECEWISE_LINEAR",
                "catalog contains aliases or callback selector");
        for (std::size_t j = 0; j < n_standard_values; ++j)
            require((standard_values[i] < standard_values[j]) == (i < j), "exact standard kernel order");
    }
    require(ib_kernel_static_initialization_valid(), "separate-TU standard value initialization");
    for (const char* name : scalar_names)
        for (const auto& spelling : spellings(name))
        {
            const IBKernel kernel(spelling);
            require(kernel.getName() == name, "scalar canonical spelling");
            require(kernel == IBKernel(name), "scalar equality");
            require(IBKernelTensorProduct(spelling) == IBKernelTensorProduct({ kernel }), "isotropic legacy name");
            require(IBKernelTensorProduct(spelling).size() == 1, "scalar name supplies one factor");
        }

    const std::pair<const char*, const char*> aliases[] = { { "PIECEWISE_CONSTANT", "BSPLINE_1" },
                                                            { "PIECEWISE_LINEAR", "BSPLINE_2" } };
    for (const auto& alias : aliases)
        for (const auto& spelling : spellings(alias.first))
        {
            require(IBKernel(spelling).getName() == alias.second, "scalar alias");
            require(IBKernelTensorProduct(spelling) == IBKernelTensorProduct(alias.second), "isotropic alias");
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
            const auto product = IBKernelTensorProduct(spelling);
            require(product.size() == 2, "composite name supplies two factors");
            require(product[0].getName() == expected.normal, "first factor");
            require(product[1].getName() == expected.tangential, "second factor");
            require(product == IBKernelTensorProduct({ IBKernel(expected.normal), IBKernel(expected.tangential) }),
                    "ordered product");
            require(product != IBKernelTensorProduct({ IBKernel(expected.tangential), IBKernel(expected.normal) }),
                    "reversed product");
            require(IBKernel(spelling).getName() == std::string(expected.name),
                    "scalar identity does not interpret composite names");
        }

    std::string name = "ApplicationKernel";
    const IBKernel custom(name);
    name = "changed";
    require(custom.getName() == "APPLICATIONKERNEL", "scalar owns its value");
    for (const auto& spelling : spellings("APPLICATIONKERNEL"))
        require(custom == IBKernel(spelling), "custom ASCII uppercase canonicalization");
    require(custom != IBKernel("USER_DEFINED"), "custom identity is distinct from legacy callback spelling");
    require(IBKernel("IB_4_custom").getName() == "IB_4_CUSTOM", "custom built-in prefix");
    require(IBKernel("composite_bspline_custom").getName() == "COMPOSITE_BSPLINE_CUSTOM",
            "custom composite-like prefix");
    require(IBKernelTensorProduct(custom.getName()) == IBKernelTensorProduct({ custom }),
            "custom isotropic interpretation");
    IBKernel normal("ApplicationKernel"), tangential("AnotherKernel");
    const IBKernelTensorProduct product({ normal, tangential });
    normal = IBKernel("changed");
    tangential = IBKernel("changed");
    require(product[0] == custom && product[1] == IBKernel("AnotherKernel"), "product owns its factors");
    require(IBKernelTensorProduct({ custom }) == IBKernelTensorProduct({ IBKernel("applicationkernel") }),
            "custom product canonicalization");
    const IBKernelTensorProduct legacy_c_string{ "DISCONTINUOUS_LINEAR" };
    require(legacy_c_string == IBKernelTensorProduct({ IBKernel::BSPLINE_2, IBKernel::BSPLINE_1 }),
            "tensor-product C-string construction");
    require(IBKernelTensorProduct({ IBKernel("IB_4"), IBKernel("IB_3") }) !=
                IBKernelTensorProduct({ IBKernel("IB_3"), IBKernel("IB_4") }),
            "non-B-spline product");
    require(IBKernelTensorProduct({ IBKernel("piecewise_constant"), IBKernel("BSPLINE_1") }) ==
                IBKernelTensorProduct({ IBKernel("BSPLINE_1"), IBKernel("BSPLINE_1") }),
            "normalized factor equality");

    require(!IBKernel::isValidName("") && !IBKernelTensorProduct::isValidName(""), "empty names are not valid");

    const IBKernel runtime_normal("IB_4"), runtime_tangential("IB_3");
    const IBKernelTensorProduct scalar(IBKernel::IB_4), single({ runtime_normal });
    const IBKernelTensorProduct repeated = { IBKernel::IB_4, IBKernel::IB_4 };
    const IBKernelTensorProduct pair({ runtime_normal, runtime_tangential });
    const IBKernelTensorProduct ordered({ runtime_normal, runtime_tangential });
    require(scalar.size() == 1 && single.size() == 1 && repeated.size() == 1 && pair.size() == 2,
            "products expose canonical axis-relative arity");
    require(scalar == single && single == repeated && scalar == repeated, "normalized isotropic transitivity");
    require(scalar == IBKernel::IB_4 && IBKernel::IB_4 == scalar, "symmetric scalar equality");
    require(pair != IBKernel::IB_4 && IBKernel::IB_4 != pair, "symmetric scalar inequality");
    require(single.isIsotropic() && repeated.isIsotropic() && !pair.isIsotropic(), "normalized isotropy");
    const IBKernelTensorProduct ordered_products[] = { IBKernelTensorProduct(IBKernel::BSPLINE_1),
                                                       IBKernelTensorProduct(
                                                           { IBKernel::BSPLINE_1, IBKernel::BSPLINE_2 }),
                                                       IBKernelTensorProduct(IBKernel::BSPLINE_2) };
    require(std::is_sorted(std::begin(ordered_products), std::end(ordered_products)),
            "tensor products have lexicographic container-key order");
    require(std::set<IBKernelTensorProduct>(std::begin(ordered_products), std::end(ordered_products)).size() ==
                std::size(ordered_products),
            "tensor products are not unique container keys");
    require(!(single < repeated) && !(repeated < single), "equal products differ in container-key order");
    require(copy_product(IBKernel::IB_4) == scalar, "implicit scalar function argument");
    require(copy_product(std::string("IB_4")) == scalar, "implicit string function argument");
    require(copy_product("IB_4") == scalar, "implicit C-string function argument");
    std::ostringstream diagnostic;
    diagnostic << single << ' ' << pair;
    require(diagnostic.str() == "(IB_4) (IB_4,IB_3)", "diagnostic describes canonical product");
    for (std::size_t slot = 0; slot < ordered.size(); ++slot)
    {
        const IBKernel expected = slot == 0 ? runtime_normal : runtime_tangential;
        require(ordered[slot] == expected, "ordered slot access");
    }

    if (rank == 0)
    {
        std::ofstream out("output");
        out << "Scalar names and aliases: PASS\n"
            << "Ordered and isotropic products: PASS\n"
            << "Custom canonical identities and two-block value copies: PASS\n"
            << "Exact encoding, maximum-length round trip, and rank order/subset independence: PASS\n"
            << "Canonical axis-relative products, equality, and diagnostics: PASS\n"
            << "Scalar/composite parsing ownership and validity predicates: PASS\n";
        require(static_cast<bool>(out), "output write failed");
    }
    return 0;
}
