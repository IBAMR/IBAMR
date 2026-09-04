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

#ifndef included_IBTK_IBKernel
#define included_IBTK_IBKernel

#include <ibtk/config.h>

#include <array>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace IBTK
{
/*!
 * \brief Value identifying a scalar one-dimensional IB kernel.
 *
 * Copies and equality use an exact two-block value, without string operations
 * or allocation. The same canonical name has the same value on every process,
 * independently of construction order. Resolve names at configuration
 * boundaries, not in numerical loops. Store names, not internal values, in
 * input and restart data.
 * The compile-time name capacity is currently 24 characters; changing it
 * requires rebuilding both the library and its applications consistently.
 *
 * The public values below are the standard scalar catalog. Standard names
 * are case-insensitive; PIECEWISE_CONSTANT and PIECEWISE_LINEAR alias BSPLINE_1
 * and BSPLINE_2. All names, including custom names, are normalized to uppercase
 * using ASCII rules. Canonical scalar names contain 1 to 24 characters from
 * A-Z, 0-9, and underscore; whitespace and other bytes are not accepted.
 * Naming a custom kernel does not register an evaluator.
 *
 */
class IBKernel
{
public:
    static const IBKernel BSPLINE_1;
    static const IBKernel BSPLINE_2;
    static const IBKernel BSPLINE_3;
    static const IBKernel BSPLINE_4;
    static const IBKernel BSPLINE_5;
    static const IBKernel BSPLINE_6;
    static const IBKernel PIECEWISE_CUBIC;
    static const IBKernel IB_3;
    static const IBKernel IB_4;
    static const IBKernel IB_4_W8;
    static const IBKernel IB_5;
    static const IBKernel IB_6;

    /*!
     * \brief Construct a scalar identity.
     *
     * The name must contain 1 to 24 ASCII letters, digits, or underscores.
     * Scalar aliases are resolved before constructing the identity.
     */
    explicit IBKernel(const std::string& name);

    //! Construct a scalar identity from a nonnull C string.
    explicit IBKernel(const char* name);

    //! Whether name satisfies the scalar kernel-name grammar.
    static bool isValidName(const std::string& name);

    /*!
     * \brief Return the standard scalar kernels.
     *
     * The order is for iteration only: it is not an ordinal, serialization,
     * or communication format.
     */
    static const std::vector<IBKernel>& getStandardKernels();

    //! Reconstruct the canonical name for configuration or diagnostics, not dispatch.
    std::string getName() const;

    bool operator==(const IBKernel& other) const;

    bool operator!=(const IBKernel& other) const;

private:
    static constexpr std::size_t MAX_NAME_LENGTH = 24;
    static constexpr std::size_t DIGITS_PER_BLOCK = 12;
    static constexpr std::uint64_t ENCODING_BASE = 38;
    static constexpr std::size_t NAME_BLOCK_COUNT = (MAX_NAME_LENGTH + DIGITS_PER_BLOCK - 1) / DIGITS_PER_BLOCK;

    static constexpr std::array<std::uint64_t, NAME_BLOCK_COUNT> encode_name(std::string_view name);

    explicit constexpr IBKernel(const std::array<std::uint64_t, NAME_BLOCK_COUNT>& name) : d_name(name)
    {
    }

    std::array<std::uint64_t, NAME_BLOCK_COUNT> d_name;
};
} // namespace IBTK

#include <ibtk/private/IBKernel-inl.h>

#endif
