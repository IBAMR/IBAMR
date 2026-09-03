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

#include <iosfwd>
#include <vector>

namespace IBTK
{
/*!
 * \brief Value describing a tensor product of scalar IBKernel factors.
 *
 * Stores exactly the supplied nonempty ordered sequence, of any length.
 * Equality compares both count and ordered identities: {K} differs from {K,K}.
 * Construction and copying may allocate. Numerical consumers take descriptions
 * by reference. Expansion and directional meaning belong to the consuming operation;
 * see \ref le_interactor_kernel_mapping for LEInteractor's convention.
 * See \ref ib_kernel_catalog for built-ins and aliases. A description does not
 * imply that a particular operation can execute it.
 */
class IBKernelTensorProduct
{
public:
    //! Own the supplied factors; throws std::invalid_argument for an empty collection.
    explicit IBKernelTensorProduct(std::vector<IBKernel> factors);

    //! Number of supplied factors, without dimensional expansion.
    std::size_t getNumberOfKernels() const;

    //! Access a factor; throws std::out_of_range for slot >= getNumberOfKernels().
    const IBKernel& getKernel(std::size_t slot) const;

    //! Whether every supplied factor equals the first (true for one factor).
    //! This does not imply that a consumer supports the count or the kernel.
    bool isIsotropic() const;

    /*!
     * \brief Interpret an existing kernel name or a custom scalar identity.
     *
     * Scalar names supply one factor. Composite names supply the two ordered
     * factors documented in \ref ib_kernel_catalog. No dimensional expansion occurs.
     *
     * \throws std::invalid_argument if the scalar name violates the canonical
     * name rules in \ref ib_kernel_catalog.
     */
    static IBKernelTensorProduct from_name(const std::string& name);

    bool operator==(const IBKernelTensorProduct& other) const;

    bool operator!=(const IBKernelTensorProduct& other) const;

private:
    std::vector<IBKernel> d_factors;
};

//! Diagnostic output only; numerical dispatch uses identities, not names.
std::ostream& operator<<(std::ostream& stream, const IBKernelTensorProduct& kernel);
} // namespace IBTK

#endif
