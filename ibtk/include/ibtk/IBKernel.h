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

#include <string>

namespace IBTK
{
/*!
 * \brief Value identifying a scalar one-dimensional IB kernel.
 *
 * Built-in names are recognized case-insensitively and stored with their
 * canonical spelling. PIECEWISE_CONSTANT and PIECEWISE_LINEAR are aliases
 * for BSPLINE_1 and BSPLINE_2. Other nonempty scalar names are application-defined
 * identities: they are owned verbatim and compared case-sensitively, without
 * registration.
 *
 * This is a description, not an evaluator or a declaration of backend support.
 * A tensor-product composition is described separately by IBKernelTensorProduct.
 */
class IBKernel
{
public:
    /*!
     * \brief Construct a scalar identity.
     *
     * \throws std::invalid_argument if the name is empty or is a recognized
     * composite name (including DISCONTINUOUS_LINEAR). Use
     * IBKernelTensorProduct::from_name() to interpret those composite names.
     */
    explicit IBKernel(const std::string& name);

    //! Return the owned scalar name, with built-in aliases normalized.
    const std::string& getName() const;

    bool operator==(const IBKernel& other) const;

    bool operator!=(const IBKernel& other) const;

private:
    std::string d_name;
};
} // namespace IBTK

#endif
