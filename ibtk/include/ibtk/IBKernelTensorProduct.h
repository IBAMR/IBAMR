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

#ifndef included_IBTK_IBKernelTensorProduct
#define included_IBTK_IBKernelTensorProduct

#include <ibtk/config.h>

#include <ibtk/IBKernel.h>

namespace IBTK
{
/*!
 * \brief Value describing a tensor product of scalar IBKernel factors.
 *
 * One factor denotes isotropic use in every coordinate direction. In the
 * two-factor side-centered form, the first factor applies in the face-normal
 * direction and the second in every face-tangential direction. Equality compares
 * the ordered factors, so K and (K,K) describe the same product.
 *
 * Describing a product does not imply that any particular interpolation or
 * spreading backend can execute it. No evaluation, storage extent, or execution
 * registration is attached to these values.
 */
class IBKernelTensorProduct
{
public:
    //! Construct the isotropic product (kernel,kernel).
    explicit IBKernelTensorProduct(const IBKernel& kernel) : IBKernelTensorProduct(kernel, kernel)
    {
    }

    //! Construct an ordered face-normal/face-tangential product.
    IBKernelTensorProduct(const IBKernel& normal, const IBKernel& tangential)
        : d_normal(normal), d_tangential(tangential)
    {
    }

    /*!
     * \brief Interpret an existing kernel name or a custom scalar identity.
     *
     * Scalar names denote isotropic products. Recognized COMPOSITE_BSPLINE_XY
     * names (12, 21, 23, 32, 34, 43, 45, 54, 56, and 65) denote the ordered
     * BSPLINE_X/BSPLINE_Y factors. DISCONTINUOUS_LINEAR aliases the 21 product.
     * These names are recognized case-insensitively. All other names follow
     * IBKernel's scalar-name rules; no general composite expression is parsed.
     *
     * \throws std::invalid_argument if the name is empty.
     */
    static IBKernelTensorProduct fromName(const std::string& name);

    const IBKernel& getNormalKernel() const
    {
        return d_normal;
    }

    const IBKernel& getTangentialKernel() const
    {
        return d_tangential;
    }

    bool operator==(const IBKernelTensorProduct& other) const
    {
        return d_normal == other.d_normal && d_tangential == other.d_tangential;
    }

    bool operator!=(const IBKernelTensorProduct& other) const
    {
        return !(*this == other);
    }

private:
    IBKernel d_normal, d_tangential;
};
} // namespace IBTK

#endif
