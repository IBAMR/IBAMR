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
#include <iosfwd>

namespace IBTK
{
/*!
 * \brief Value describing a tensor product of scalar IBKernel factors.
 *
 * Stores NDIM ordered factors. One factor repeats in every slot; two factors
 * expand to {A,B} in 2D and {A,B,B} in 3D. Equality compares the expanded
 * factors. Their directional meaning belongs to the consuming operation;
 * see \ref le_interactor_kernel_mapping for LEInteractor's convention.
 * See \ref ib_kernel_catalog for built-ins and aliases. A description does not
 * imply that a particular operation can execute it.
 */
class IBKernelTensorProduct
{
public:
    //! Repeat the scalar factor in every slot.
    explicit IBKernelTensorProduct(const IBKernel& kernel);

    //! Store the first factor followed by repetitions of the second.
    IBKernelTensorProduct(const IBKernel& first, const IBKernel& second);

#if (NDIM == 3)
    //! Specify all three ordered factors.
    IBKernelTensorProduct(const IBKernel& first, const IBKernel& second, const IBKernel& third);
#endif

    //! Specify all ordered factors.
    explicit IBKernelTensorProduct(const std::array<IBKernel, NDIM>& factors);

    //! Access a factor by slot; throws std::out_of_range for slot >= NDIM.
    const IBKernel& getKernel(unsigned int slot) const;

    bool isIsotropic() const;

    /*!
     * \brief Interpret an existing kernel name or a custom scalar identity.
     *
     * Scalar names denote isotropic products. Composite names supply the
     * ordered factors documented in \ref ib_kernel_catalog.
     *
     * \throws std::invalid_argument if the scalar name violates the canonical
     * name rules in \ref ib_kernel_catalog.
     */
    static IBKernelTensorProduct from_name(const std::string& name);

    bool operator==(const IBKernelTensorProduct& other) const;

    bool operator!=(const IBKernelTensorProduct& other) const;

private:
    std::array<IBKernel, NDIM> d_factors;
};

//! Diagnostic output only; numerical dispatch uses identities, not names.
std::ostream& operator<<(std::ostream& stream, const IBKernelTensorProduct& kernel);
} // namespace IBTK

#endif
