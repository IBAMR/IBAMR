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
 * Copies, equality, and ordering use an exact compact value, without string
 * operations or allocation. The same canonical name has the same value on
 * every process, independently of construction order. Store names, not this
 * internal value, in input and restart data.
 *
 * Names are case-insensitive and may contain at most MAX_NAME_LENGTH ASCII
 * letters, digits, or underscores. PIECEWISE_CONSTANT and PIECEWISE_LINEAR are
 * aliases for BSPLINE_1 and BSPLINE_2. Naming a custom kernel does not register
 * an evaluator.
 */
class IBKernel
{
public:
    /*! \name Standard scalar kernels */
    //\{
    static const IBKernel BSPLINE_1;
    static const IBKernel BSPLINE_2;
    static const IBKernel BSPLINE_3;
    static const IBKernel BSPLINE_4;
    static const IBKernel BSPLINE_5;
    static const IBKernel BSPLINE_6;
    static const IBKernel IB_3;
    static const IBKernel IB_4;
    static const IBKernel IB_4_W8;
    static const IBKernel IB_5;
    static const IBKernel IB_6;
    static const IBKernel PIECEWISE_CUBIC;
    //\}

    /*!
     * \brief Construct a scalar identity.
     *
     * \param name Name containing at most MAX_NAME_LENGTH ASCII letters,
     * digits, or underscores. Scalar aliases are accepted.
     */
    explicit IBKernel(const std::string& name);

    /*! \brief Construct a scalar identity from a nonnull C string. */
    explicit IBKernel(const char* name);

    /*! \brief Return whether \p name satisfies the scalar kernel-name grammar. */
    static bool isValidName(const std::string& name);

    /*!
     * \brief Return the standard scalar kernels.
     *
     * The order is for iteration only: it is not an ordinal, serialization,
     * or communication format.
     */
    static const std::vector<IBKernel>& getStandardKernels();

    /*! \brief Reconstruct the canonical name for configuration or diagnostics. */
    std::string getName() const;

    /*! \brief Return whether two kernels have the same canonical identity. */
    bool operator==(const IBKernel& other) const;

    /*! \brief Return whether two kernels have different canonical identities. */
    bool operator!=(const IBKernel& other) const;

    /*! \brief Order kernels by canonical identity for use as container keys. */
    bool operator<(const IBKernel& other) const;

private:
    //! Maximum number of characters in a canonical kernel name.
    static constexpr std::size_t MAX_NAME_LENGTH = 24;

    //! Number of encoded base-38 digits stored in each integer block.
    static constexpr std::size_t DIGITS_PER_BLOCK = 12;

    //! Radix used by the compact canonical-name representation.
    static constexpr std::uint64_t ENCODING_BASE = 38;

    //! Number of integer blocks in the compact representation.
    static constexpr std::size_t NAME_BLOCK_COUNT = (MAX_NAME_LENGTH + DIGITS_PER_BLOCK - 1) / DIGITS_PER_BLOCK;

    /*! \brief Encode a validated canonical name into its compact identity. */
    static constexpr std::array<std::uint64_t, NAME_BLOCK_COUNT> encode_name(std::string_view name);

    /*! \brief Construct a standard value from a validated canonical name. */
    static constexpr IBKernel from_canonical_name(std::string_view name);

    /*! \brief Construct directly from a compact canonical identity. */
    explicit constexpr IBKernel(const std::array<std::uint64_t, NAME_BLOCK_COUNT>& name);

    //! Compact canonical kernel identity.
    std::array<std::uint64_t, NAME_BLOCK_COUNT> d_name;
};
} // namespace IBTK

#include <ibtk/private/IBKernel-inl.h>

#endif
