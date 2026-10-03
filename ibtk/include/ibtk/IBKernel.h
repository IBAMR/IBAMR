// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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
#include <compare>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace IBTK
{
/*!
 * \brief Class IBKernel identifies a one-dimensional IB kernel by name.
 *
 * The name is stored packed into a fixed number of integers rather than as a
 * std::string or a closed enumeration. A closed enumeration would compare and
 * copy as cheaply as this encoding does, but adding a name to it needs a
 * library change, so a consumer could never accept a name it does not already
 * know about; a std::string can name any kernel, but every comparison, and
 * every dispatch keyed on the kernel (interpolation and spreading, called for
 * every IB point on every patch, every timestep), would then cost a string
 * comparison; master's LEInteractor dispatched on individual characters of the name to avoid exactly that cost. This
 * encoding gives an unrestricted name (like a string) at the comparison and copy cost of a small fixed-size integer
 * array (like an enum), without needing a character-by-character dispatch.
 *
 * IBKernel itself does not restrict which names may be constructed, other
 * than the reserved prefixes below, but recognizing a name is up to each
 * consumer: LEInteractor recognizes a fixed catalog of built-in kernels plus
 * the single name USER_DEFINED, whose implementation an application supplies
 * through LEInteractor::s_kernel_fcn; it does not otherwise accept
 * application-chosen names.
 *
 * Names are case-insensitive. PIECEWISE_CONSTANT and PIECEWISE_LINEAR are
 * aliases for BSPLINE_1 and BSPLINE_2, respectively.
 * Names that begin with COMPOSITE_BSPLINE_, and the name DISCONTINUOUS_LINEAR,
 * are reserved for IBKernelTensorProduct and are not valid scalar kernel names.
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
     * Use is_valid_name() to check whether a name is accepted.
     *
     * \param name Name of 1 to 24 ASCII letters, digits, or underscores. Names are
     * case-insensitive; scalar aliases are accepted.
     */
    explicit IBKernel(const std::string& name);

    /*!
     * \brief Construct a kernel from its name.
     *
     * \param name Nonnull C string satisfying the name requirements of
     * IBKernel(const std::string&). A null pointer is an error.
     */
    explicit IBKernel(const char* name);

    /*! \brief Return whether \p name satisfies the scalar kernel-name requirements. */
    static bool is_valid_name(const std::string& name);

    /*! \brief Return the standard scalar kernels, in the order of their comparison. */
    static const std::vector<IBKernel>& get_standard_kernels();

    /*! \brief Return the uppercase kernel name, with aliases replaced by standard names. */
    std::string getName() const;

    /*! \brief Return whether two kernels have the same name after resolving aliases. */
    bool operator==(const IBKernel& other) const = default;

    /*! \brief Compare encoded names to order kernels in associative containers. */
    std::strong_ordering operator<=>(const IBKernel& other) const = default;

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

    //! Powers used to pad partially filled encoded blocks.
    static const std::array<std::uint64_t, DIGITS_PER_BLOCK + 1> ENCODING_BASE_POWERS;

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
