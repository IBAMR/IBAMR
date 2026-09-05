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
 * \brief Class IBKernelTensorProduct describes a tensor product of scalar IB kernels.
 *
 * A single factor applies in every coordinate direction. With two factors,
 * the first applies along an axis selected by the calling code and the second
 * applies in the remaining directions. For example, LEInteractor selects the
 * face-normal direction for side-centered data. Two equal factors are stored as a
 * single factor.
 */
class IBKernelTensorProduct
{
public:
    /*! \brief Construct a tensor product using the same kernel in every direction. */
    IBKernelTensorProduct(const IBKernel& factor);

    /*! \brief Construct a tensor product from one or two scalar kernels. */
    IBKernelTensorProduct(std::initializer_list<IBKernel> factors);

    /*! \brief Construct a tensor product from a scalar or composite kernel name. */
    IBKernelTensorProduct(const std::string& name);

    /*! \brief Construct a tensor product from a nonnull scalar or composite kernel name. */
    IBKernelTensorProduct(const char* name);

    /*! \brief Return whether \p name is a valid scalar or composite kernel name. */
    static bool isValidName(const std::string& name);

    /*! \brief Return the number of factors, treating two equal factors as one. */
    std::size_t size() const;

    /*! \brief Access a factor. Valid indices satisfy \p slot < size(). */
    const IBKernel& operator[](std::size_t slot) const;

    /*! \brief Return whether the same scalar kernel applies in every direction. */
    bool isIsotropic() const;

    /*! \brief Return whether two products have the same factors after combining equal pairs. */
    bool operator==(const IBKernelTensorProduct& other) const;

    /*! \brief Return whether two products differ. */
    bool operator!=(const IBKernelTensorProduct& other) const;

    /*! \brief Compare factors lexicographically using IBKernel::operator<(). */
    bool operator<(const IBKernelTensorProduct& other) const;

private:
    //! Minimum number of factors accepted by the constructors.
    static constexpr std::size_t MIN_ACTIVE_FACTORS = 1;

    //! Maximum number of factors accepted by the constructors.
    static constexpr std::size_t MAX_ACTIVE_FACTORS = 2;

    //! Storage capacity, with room for an additional factor.
    static constexpr std::size_t INLINE_FACTOR_CAPACITY = 3;

    //! Factors and their count after combining equal pairs.
    struct CanonicalFactors
    {
        //! Scalar kernels; only the first size entries are used.
        std::array<IBKernel, INLINE_FACTOR_CAPACITY> factors;

        //! Number of factors in use.
        std::size_t size;
    };

    /*! \brief Construct from factors with equal pairs already combined. */
    explicit IBKernelTensorProduct(CanonicalFactors factors);

    /*! \brief Check the factor count and combine two equal factors into one. */
    static CanonicalFactors canonicalize(std::initializer_list<IBKernel> factors);

    /*! \brief Parse a scalar or composite kernel name. */
    static CanonicalFactors parse_name(const std::string& name);

    /*! \brief Parse a scalar or composite kernel name, rejecting null pointers. */
    static CanonicalFactors parse_name(const char* name);

    //! Scalar kernels in the tensor product.
    std::array<IBKernel, INLINE_FACTOR_CAPACITY> d_factors;

    //! Number of active factors in d_factors.
    std::size_t d_size;
};

//! Compare an isotropic product with a scalar kernel.
bool operator==(const IBKernelTensorProduct& product, const IBKernel& kernel);

//! Compare a scalar kernel with an isotropic product.
bool operator==(const IBKernel& kernel, const IBKernelTensorProduct& product);

bool operator!=(const IBKernelTensorProduct& product, const IBKernel& kernel);

bool operator!=(const IBKernel& kernel, const IBKernelTensorProduct& product);

//! Write the scalar kernel names as a parenthesized, comma-separated list.
std::ostream& operator<<(std::ostream& stream, const IBKernelTensorProduct& kernel);
} // namespace IBTK

#endif
