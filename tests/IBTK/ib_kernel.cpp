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

#include "../tests.h"

using IBTK::IBKernel;
using IBTK::IBKernelTensorProduct;

bool ib_kernel_static_initialization_valid();

namespace
{
constexpr std::size_t ENCODED_BLOCK_COUNT = 2;

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
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Logger::Appender> abort_appender = new TestAppender();
    SAMRAI::tbox::Logger::getInstance()->setAbortAppender(abort_appender);
    if (input_file.find("invalid_scalar") != std::string::npos)
    {
        SAMRAI::tbox::PIO::logOnlyNodeZero("output");
        const IBKernel invalid("BAD NAME");
        return 0;
    }
    if (input_file.find("empty_product") != std::string::npos)
    {
        SAMRAI::tbox::PIO::logOnlyNodeZero("output");
        const IBKernelTensorProduct invalid(std::initializer_list<IBKernel>{});
        return 0;
    }
    if (input_file.find("null_product") != std::string::npos)
    {
        SAMRAI::tbox::PIO::logOnlyNodeZero("output");
        const IBKernelTensorProduct invalid(static_cast<const char*>(nullptr));
        return 0;
    }
    if (input_file.find("three_factor") != std::string::npos)
    {
        SAMRAI::tbox::PIO::logOnlyNodeZero("output");
        const IBKernelTensorProduct invalid({ IBKernel::IB_4, IBKernel::IB_4, IBKernel::IB_4 });
        return 0;
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
    TBOX_ASSERT(common_blocks == root_blocks);

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
        TBOX_ASSERT(numeric_blocks(kernel) == expected.blocks);
        TBOX_ASSERT(copy == kernel && numeric_blocks(copy) == expected.blocks);
        TBOX_ASSERT(copy.getName() == expected.name);
    }
    TBOX_ASSERT(IBKernel("ABCDEFGHIJKL") != IBKernel("ABCDEFGHIJKLM"));
    TBOX_ASSERT(IBKernel("ABCDEFGHIJKLMNOPQRSTUVWX") != IBKernel("ABCDEFGHIJKLMNOPQRSTUVW_"));
    TBOX_ASSERT(IBKernel("A") != IBKernel("A_") && IBKernel("A") != IBKernel("AA"));
    const char* const scalar_c_string = "piecewise_constant";
    TBOX_ASSERT(IBKernel{ scalar_c_string } == IBKernel(IBKernel::BSPLINE_1));
    TBOX_ASSERT(!IBKernel::isValidName("ABCDEFGHIJKLMNOPQRSTUVWXY"));
    for (const std::string& invalid :
         { std::string("A B"), std::string("A-B"), std::string("A\0B", 3), std::string("A\x80", 2) })
        TBOX_ASSERT(!IBKernel::isValidName(invalid));

    const IBKernel standard_values[] = { IBKernel::BSPLINE_1, IBKernel::BSPLINE_2, IBKernel::BSPLINE_3,
                                         IBKernel::BSPLINE_4, IBKernel::BSPLINE_5, IBKernel::BSPLINE_6,
                                         IBKernel::IB_3,      IBKernel::IB_4,      IBKernel::IB_4_W8,
                                         IBKernel::IB_5,      IBKernel::IB_6,      IBKernel::PIECEWISE_CUBIC };
    const char* scalar_names[] = { "BSPLINE_1", "BSPLINE_2", "BSPLINE_3", "BSPLINE_4", "BSPLINE_5", "BSPLINE_6",
                                   "IB_3",      "IB_4",      "IB_4_W8",   "IB_5",      "IB_6",      "PIECEWISE_CUBIC" };
    const auto& catalog = IBKernel::getStandardKernels();
    const std::size_t n_standard_values = sizeof(standard_values) / sizeof(standard_values[0]);
    TBOX_ASSERT(catalog.size() == n_standard_values);
    TBOX_ASSERT(std::is_sorted(catalog.begin(), catalog.end()));
    TBOX_ASSERT(std::set<IBKernel>(catalog.begin(), catalog.end()).size() == catalog.size());
    for (std::size_t i = 0; i < n_standard_values; ++i)
    {
        TBOX_ASSERT(standard_values[i] == IBKernel(scalar_names[i]));
        TBOX_ASSERT(catalog[i] == standard_values[i]);
        TBOX_ASSERT(catalog[i].getName() != "USER_DEFINED" && catalog[i].getName() != "PIECEWISE_CONSTANT" &&
                    catalog[i].getName() != "PIECEWISE_LINEAR");
        for (std::size_t j = 0; j < n_standard_values; ++j)
            TBOX_ASSERT((standard_values[i] < standard_values[j]) == (i < j));
    }
    TBOX_ASSERT(ib_kernel_static_initialization_valid());
    for (const char* name : scalar_names)
        for (const auto& spelling : spellings(name))
        {
            const IBKernel kernel(spelling);
            TBOX_ASSERT(kernel.getName() == name);
            TBOX_ASSERT(kernel == IBKernel(name));
            TBOX_ASSERT(IBKernelTensorProduct(spelling) == IBKernelTensorProduct({ kernel }));
            TBOX_ASSERT(IBKernelTensorProduct(spelling).size() == 1);
        }

    const std::pair<const char*, const char*> aliases[] = { { "PIECEWISE_CONSTANT", "BSPLINE_1" },
                                                            { "PIECEWISE_LINEAR", "BSPLINE_2" } };
    for (const auto& alias : aliases)
        for (const auto& spelling : spellings(alias.first))
        {
            TBOX_ASSERT(IBKernel(spelling).getName() == alias.second);
            TBOX_ASSERT(IBKernelTensorProduct(spelling) == IBKernelTensorProduct(alias.second));
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
            TBOX_ASSERT(product.size() == 2);
            TBOX_ASSERT(product[0].getName() == expected.normal);
            TBOX_ASSERT(product[1].getName() == expected.tangential);
            TBOX_ASSERT(product == IBKernelTensorProduct({ IBKernel(expected.normal), IBKernel(expected.tangential) }));
            TBOX_ASSERT(product != IBKernelTensorProduct({ IBKernel(expected.tangential), IBKernel(expected.normal) }));
            TBOX_ASSERT(IBKernel(spelling).getName() == std::string(expected.name));
        }

    std::string name = "ApplicationKernel";
    const IBKernel custom(name);
    name = "changed";
    TBOX_ASSERT(custom.getName() == "APPLICATIONKERNEL");
    for (const auto& spelling : spellings("APPLICATIONKERNEL")) TBOX_ASSERT(custom == IBKernel(spelling));
    TBOX_ASSERT(custom != IBKernel("USER_DEFINED"));
    TBOX_ASSERT(IBKernel("IB_4_custom").getName() == "IB_4_CUSTOM");
    TBOX_ASSERT(IBKernel("composite_bspline_custom").getName() == "COMPOSITE_BSPLINE_CUSTOM");
    TBOX_ASSERT(IBKernelTensorProduct(custom.getName()) == IBKernelTensorProduct({ custom }));
    IBKernel normal("ApplicationKernel"), tangential("AnotherKernel");
    const IBKernelTensorProduct product({ normal, tangential });
    normal = IBKernel("changed");
    tangential = IBKernel("changed");
    TBOX_ASSERT(product[0] == custom && product[1] == IBKernel("AnotherKernel"));
    TBOX_ASSERT(IBKernelTensorProduct({ custom }) == IBKernelTensorProduct({ IBKernel("applicationkernel") }));
    const IBKernelTensorProduct legacy_c_string{ "DISCONTINUOUS_LINEAR" };
    TBOX_ASSERT(legacy_c_string == IBKernelTensorProduct({ IBKernel::BSPLINE_2, IBKernel::BSPLINE_1 }));
    TBOX_ASSERT(IBKernelTensorProduct({ IBKernel("IB_4"), IBKernel("IB_3") }) !=
                IBKernelTensorProduct({ IBKernel("IB_3"), IBKernel("IB_4") }));
    TBOX_ASSERT(IBKernelTensorProduct({ IBKernel("piecewise_constant"), IBKernel("BSPLINE_1") }) ==
                IBKernelTensorProduct({ IBKernel("BSPLINE_1"), IBKernel("BSPLINE_1") }));

    TBOX_ASSERT(!IBKernel::isValidName("") && !IBKernelTensorProduct::isValidName(""));

    const IBKernel runtime_normal("IB_4"), runtime_tangential("IB_3");
    const IBKernelTensorProduct scalar(IBKernel::IB_4), single({ runtime_normal });
    const IBKernelTensorProduct repeated = { IBKernel::IB_4, IBKernel::IB_4 };
    const IBKernelTensorProduct pair({ runtime_normal, runtime_tangential });
    const IBKernelTensorProduct ordered({ runtime_normal, runtime_tangential });
    TBOX_ASSERT(scalar.size() == 1 && single.size() == 1 && repeated.size() == 1 && pair.size() == 2);
    TBOX_ASSERT(scalar == single && single == repeated && scalar == repeated);
    TBOX_ASSERT(scalar == IBKernel::IB_4 && IBKernel::IB_4 == scalar);
    TBOX_ASSERT(pair != IBKernel::IB_4 && IBKernel::IB_4 != pair);
    TBOX_ASSERT(single.isIsotropic() && repeated.isIsotropic() && !pair.isIsotropic());
    const IBKernelTensorProduct ordered_products[] = { IBKernelTensorProduct(IBKernel::BSPLINE_1),
                                                       IBKernelTensorProduct(
                                                           { IBKernel::BSPLINE_1, IBKernel::BSPLINE_2 }),
                                                       IBKernelTensorProduct(IBKernel::BSPLINE_2) };
    TBOX_ASSERT(std::is_sorted(std::begin(ordered_products), std::end(ordered_products)));
    TBOX_ASSERT(std::set<IBKernelTensorProduct>(std::begin(ordered_products), std::end(ordered_products)).size() ==
                std::size(ordered_products));
    TBOX_ASSERT(!(single < repeated) && !(repeated < single));
    TBOX_ASSERT(copy_product(IBKernel::IB_4) == scalar);
    TBOX_ASSERT(copy_product(std::string("IB_4")) == scalar);
    TBOX_ASSERT(copy_product("IB_4") == scalar);
    std::ostringstream diagnostic;
    diagnostic << single << ' ' << pair;
    TBOX_ASSERT(diagnostic.str() == "(IB_4) (IB_4,IB_3)");
    for (std::size_t slot = 0; slot < ordered.size(); ++slot)
    {
        const IBKernel expected = slot == 0 ? runtime_normal : runtime_tangential;
        TBOX_ASSERT(ordered[slot] == expected);
    }

    if (rank == 0)
    {
        std::ofstream out("output");
        TBOX_ASSERT(static_cast<bool>(out));
    }
    return 0;
}
