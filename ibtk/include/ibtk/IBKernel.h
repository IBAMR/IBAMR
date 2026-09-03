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

namespace IBTK
{
/*!
 * \brief Value identifying a scalar one-dimensional IB kernel.
 *
 * Copies and equality use an exact two-word value, without string operations
 * or allocation. The same canonical name has the same value on every process,
 * independently of construction order. Resolve names at configuration
 * boundaries, not in numerical loops. Store names, not internal values, in
 * input and restart data.
 *
 * \anchor ib_kernel_catalog
 * The public constants below are the built-in scalar catalog. Built-in names
 * are case-insensitive; PIECEWISE_CONSTANT and PIECEWISE_LINEAR alias BSPLINE_1
 * and BSPLINE_2. All names, including custom names, are normalized to uppercase
 * using ASCII rules. Canonical scalar names contain 1 to 24 characters from
 * A-Z, 0-9, and underscore; whitespace and other bytes are not accepted.
 * Naming a custom kernel does not register an evaluator. USER_DEFINED retains
 * its separate legacy callback meaning; it is not an alias for custom names.
 *
 * IBKernelTensorProduct::from_name() also recognizes COMPOSITE_BSPLINE_XY for
 * XY = 12, 21, 23, 32, 34, 43, 45, 54, 56, 65, with component-axis BSPLINE_X
 * and transverse BSPLINE_Y factors. DISCONTINUOUS_LINEAR aliases 21. These
 * product names cannot be used as scalar names. Execution support is a
 * separate property of the consuming operation.
 */
class IBKernel
{
public:
    enum Builtin
    {
        BSPLINE_1,
        BSPLINE_2,
        BSPLINE_3,
        BSPLINE_4,
        BSPLINE_5,
        BSPLINE_6,
        PIECEWISE_CUBIC,
        IB_3,
        IB_4,
        IB_4_W8,
        IB_5,
        IB_6,
        USER_DEFINED
    };

    //! Construct a built-in identity without name lookup.
    IBKernel(Builtin kernel);

    /*!
     * \brief Construct a scalar identity.
     *
     * \throws std::invalid_argument if the canonical name is empty, exceeds
     * 24 characters, contains an invalid character, or is a recognized
     * composite name (including DISCONTINUOUS_LINEAR). Aliases are resolved
     * before checking the canonical scalar length. Use
     * IBKernelTensorProduct::from_name() to interpret those composite names.
     */
    explicit IBKernel(const std::string& name);

    //! Reconstruct the canonical name for configuration or diagnostics, not dispatch.
    std::string getName() const;

    bool operator==(const IBKernel& other) const;

    bool operator!=(const IBKernel& other) const;

private:
    std::array<std::uint64_t, 2> d_name;
};
} // namespace IBTK

#include <ibtk/private/IBKernel-inl.h>

#endif
