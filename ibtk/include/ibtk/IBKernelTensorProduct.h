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

#include <array>
#include <initializer_list>
#include <iosfwd>
#include <string>

namespace IBTK
{
/*!
 * \brief Value describing a tensor product of scalar IBKernel factors.
 *
 * Owns one or two axis-relative scalar factors. One factor applies in all
 * directions. Two factors apply respectively on a consumer-defined
 * distinguished axis and on every transverse axis. An equal pair is stored as
 * the equivalent one-factor value.
 */
class IBKernelTensorProduct
{
public:
    //! Implicitly promote a scalar kernel to an isotropic tensor product.
    IBKernelTensorProduct(const IBKernel& factor);

    //! Canonicalize one or two ordered factors.
    IBKernelTensorProduct(std::initializer_list<IBKernel> factors);

    //! Implicitly interpret a legacy scalar or composite kernel name.
    IBKernelTensorProduct(const std::string& name);

    //! Implicitly interpret a nonnull legacy scalar or composite C string.
    IBKernelTensorProduct(const char* name);

    //! Whether name can construct a tensor-product description.
    static bool isValidName(const std::string& name);

    //! Number of meaningful factors, one or two.
    std::size_t size() const;

    //! Access a factor. Valid indices satisfy slot < size().
    const IBKernel& operator[](std::size_t slot) const;

    //! Whether the value has one all-directions factor.
    //! This does not imply that a consumer supports the kernel.
    bool isIsotropic() const;

    bool operator==(const IBKernelTensorProduct& other) const;

    bool operator!=(const IBKernelTensorProduct& other) const;

private:
    struct CanonicalFactors
    {
        std::array<IBKernel, 3> factors;
        std::size_t size;
    };

    explicit IBKernelTensorProduct(CanonicalFactors factors);

    static CanonicalFactors canonicalize(std::initializer_list<IBKernel> factors);

    static CanonicalFactors parse_name(const std::string& name);

    std::array<IBKernel, 3> d_factors;
    std::size_t d_size;
};

//! Compare an isotropic product with a scalar kernel.
bool operator==(const IBKernelTensorProduct& product, const IBKernel& kernel);

//! Compare a scalar kernel with an isotropic product.
bool operator==(const IBKernel& kernel, const IBKernelTensorProduct& product);

bool operator!=(const IBKernelTensorProduct& product, const IBKernel& kernel);

bool operator!=(const IBKernel& kernel, const IBKernelTensorProduct& product);

//! Diagnostic output only; numerical dispatch uses identities, not names.
std::ostream& operator<<(std::ostream& stream, const IBKernelTensorProduct& kernel);
} // namespace IBTK

#endif
