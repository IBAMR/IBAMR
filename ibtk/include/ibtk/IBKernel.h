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
 * \brief Class IBKernel identifies a one-dimensional IB kernel by name.
 *
 * Names are case-insensitive. PIECEWISE_CONSTANT and PIECEWISE_LINEAR are
 * aliases for BSPLINE_1 and BSPLINE_2, respectively. Applications may define
 * additional names, but must provide the corresponding kernel implementations.
 *
 * The same name has the same encoded value on every MPI process. Use getName()
 * to write kernel names to input or restart files.
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
     * \brief Construct a kernel from its name.
     *
     * \param name Nonempty name containing at most MAX_NAME_LENGTH ASCII letters,
     * digits, or underscores. Scalar aliases are accepted.
     */
    explicit IBKernel(const std::string& name);

    /*!
     * \brief Construct a kernel from its name.
     *
     * \param name Nonnull C string containing a nonempty name of at most
     * MAX_NAME_LENGTH ASCII letters, digits, or underscores. Scalar aliases
     * are accepted.
     */
    explicit IBKernel(const char* name);

    /*! \brief Return whether \p name satisfies the scalar kernel-name requirements. */
    static bool isValidName(const std::string& name);

    /*! \brief Return the standard scalar kernels in alphabetical order. */
    static const std::vector<IBKernel>& getStandardKernels();

    /*! \brief Return the uppercase kernel name, with aliases replaced by standard names. */
    std::string getName() const;

    /*! \brief Return whether two kernels have the same name after resolving aliases. */
    bool operator==(const IBKernel& other) const;

    /*! \brief Return whether two kernels differ. */
    bool operator!=(const IBKernel& other) const;

    /*! \brief Compare encoded names to order kernels in associative containers. */
    bool operator<(const IBKernel& other) const;

private:
    //! Maximum number of characters in a canonical kernel name.
    static constexpr std::size_t MAX_NAME_LENGTH = 24;

    //! Number of encoded characters stored in each integer block.
    static constexpr std::size_t DIGITS_PER_BLOCK = 12;

    //! Base used to encode kernel names.
    static constexpr std::uint64_t ENCODING_BASE = 38;

    //! Encoded value of 'A'.
    static constexpr unsigned int LETTER_DIGIT_OFFSET = 1;

    //! Encoded value of '0'.
    static constexpr unsigned int NUMBER_DIGIT_OFFSET = 27;

    //! Encoded value reserved for underscore.
    static constexpr unsigned int UNDERSCORE_DIGIT = 37;

    //! Number of integer blocks needed to store a kernel name.
    static constexpr std::size_t NAME_BLOCK_COUNT = (MAX_NAME_LENGTH + DIGITS_PER_BLOCK - 1) / DIGITS_PER_BLOCK;

    /*! \brief Encode a valid uppercase name with aliases already resolved. */
    static constexpr std::array<std::uint64_t, NAME_BLOCK_COUNT> encode_name(std::string_view name);

    /*! \brief Resolve aliases and encode a name; return false if the name is invalid. */
    static bool try_encode_name(std::string_view name, std::array<std::uint64_t, NAME_BLOCK_COUNT>& encoded_name);

    /*! \brief Construct a kernel from a valid uppercase name with aliases already resolved. */
    static constexpr IBKernel from_canonical_name(std::string_view name);

    /*! \brief Construct a kernel from an encoded name. */
    explicit constexpr IBKernel(const std::array<std::uint64_t, NAME_BLOCK_COUNT>& name);

    //! Encoded kernel name.
    std::array<std::uint64_t, NAME_BLOCK_COUNT> d_name;
};
} // namespace IBTK

#include <ibtk/private/IBKernel-inl.h>

#endif
