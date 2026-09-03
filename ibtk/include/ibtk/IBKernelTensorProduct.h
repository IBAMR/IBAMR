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

#include <initializer_list>
#include <iosfwd>
#include <vector>

namespace IBTK
{
/*!
 * \brief Value describing a tensor product of scalar IBKernel factors.
 *
 * Owns a nonempty ordered sequence of scalar factors. It may be constructed
 * from factors or from a legacy scalar or composite kernel name.
 */
class IBKernelTensorProduct
{
public:
    //! Own the supplied nonempty sequence of factors.
    explicit IBKernelTensorProduct(std::vector<IBKernel> factors);

    //! Own the supplied nonempty sequence of factors.
    IBKernelTensorProduct(std::initializer_list<IBKernel> factors);

    //! Interpret a legacy scalar or composite kernel name.
    explicit IBKernelTensorProduct(const std::string& name);

    //! Whether name can construct a tensor-product description.
    static bool isValidName(const std::string& name);

    //! Number of supplied factors.
    std::size_t size() const;

    //! Access a factor. Valid indices satisfy slot < size().
    const IBKernel& operator[](std::size_t slot) const;

    //! Whether every supplied factor equals the first (true for one factor).
    //! This does not imply that a consumer supports the count or the kernel.
    bool isIsotropic() const;

    bool operator==(const IBKernelTensorProduct& other) const;

    bool operator!=(const IBKernelTensorProduct& other) const;

private:
    std::vector<IBKernel> d_factors;
};

//! Diagnostic output only; numerical dispatch uses identities, not names.
std::ostream& operator<<(std::ostream& stream, const IBKernelTensorProduct& kernel);
} // namespace IBTK

#endif
