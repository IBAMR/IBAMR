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
 * An isotropic product repeats one factor. A component-relative product uses
 * one factor along the component axis and another along every transverse axis:
 * for side data these are the face-normal and tangential directions. Cartesian
 * products instead fix factors in x/y/z order, independently of the component.
 * In 3D, Cartesian factors can distinguish both tangential directions.
 * Equality compares these mappings, with equal-factor forms equivalent.
 * See \ref ib_kernel_catalog for built-ins and aliases. A description does not
 * imply that a particular operation can execute it.
 */
class IBKernelTensorProduct
{
public:
    //! Repeat the scalar factor in every spatial direction.
    explicit IBKernelTensorProduct(const IBKernel& kernel);

    //! Construct an ordered face-normal/face-tangential product.
    IBKernelTensorProduct(const IBKernel& normal, const IBKernel& tangential);

    //! Construct factors in Cartesian x/y/z order, independent of component axis.
    explicit IBKernelTensorProduct(const std::array<IBKernel, NDIM>& factors);

    //! Resolve the factor for a Cartesian axis and a vector component axis.
    const IBKernel& getKernel(unsigned int axis, unsigned int component_axis) const;

    bool isIsotropic() const;

    bool isComponentRelative() const;

    /*!
     * \brief Interpret an existing kernel name or a custom scalar identity.
     *
     * Scalar names denote isotropic products. Composite names use the
     * component-relative mapping documented in \ref ib_kernel_catalog.
     *
     * \throws std::invalid_argument if the scalar name violates the canonical
     * name rules in \ref ib_kernel_catalog.
     */
    static IBKernelTensorProduct from_name(const std::string& name);

    //! Component-axis factor; throws std::logic_error for unequal Cartesian factors.
    const IBKernel& getNormalKernel() const;

    //! Transverse factor; throws std::logic_error for unequal Cartesian factors.
    const IBKernel& getTangentialKernel() const;

    bool operator==(const IBKernelTensorProduct& other) const;

    bool operator!=(const IBKernelTensorProduct& other) const;

private:
    std::array<IBKernel, NDIM> d_factors;
    bool d_component_relative;
};

//! Diagnostic output only; numerical dispatch uses identities, not names.
std::ostream& operator<<(std::ostream& stream, const IBKernelTensorProduct& kernel);
} // namespace IBTK

#endif
