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
 * Owns a canonical description with one or two axis-relative scalar factors.
 * One factor applies in every direction. Two factors apply respectively on a
 * consumer-defined distinguished axis and on every transverse axis. An equal
 * pair is stored as the equivalent one-factor value.
 */
class IBKernelTensorProduct
{
public:
    /*! \brief Implicitly promote a scalar kernel to an isotropic tensor product. */
    IBKernelTensorProduct(const IBKernel& factor);

    /*! \brief Canonicalize one or two ordered factors. */
    IBKernelTensorProduct(std::initializer_list<IBKernel> factors);

    /*! \brief Implicitly interpret a legacy scalar or composite kernel name. */
    IBKernelTensorProduct(const std::string& name);

    /*! \brief Implicitly interpret a nonnull legacy scalar or composite C string. */
    IBKernelTensorProduct(const char* name);

    /*! \brief Return whether \p name can construct a tensor-product description. */
    static bool isValidName(const std::string& name);

    /*! \brief Return the number of factors in the canonical description. */
    std::size_t size() const;

    /*! \brief Access a factor. Valid indices satisfy \p slot < size(). */
    const IBKernel& operator[](std::size_t slot) const;

    /*! \brief Return whether the same scalar kernel applies in every direction. */
    bool isIsotropic() const;

    /*! \brief Return whether two products have the same canonical description. */
    bool operator==(const IBKernelTensorProduct& other) const;

    /*! \brief Return whether two products have different canonical descriptions. */
    bool operator!=(const IBKernelTensorProduct& other) const;

    /*! \brief Lexicographically order canonical descriptions for use as container keys. */
    bool operator<(const IBKernelTensorProduct& other) const;

private:
    //! Minimum number of factors in a valid public description.
    static constexpr std::size_t MIN_ACTIVE_FACTORS = 1;

    //! Maximum number of factors in the current public description.
    static constexpr std::size_t MAX_ACTIVE_FACTORS = 2;

    //! Inline factor capacity, including reserved headroom for a future role.
    static constexpr std::size_t INLINE_FACTOR_CAPACITY = 3;

    //! Canonical factor storage returned by construction helpers.
    struct CanonicalFactors
    {
        //! Inline factor values; entries beyond size are not observable.
        std::array<IBKernel, INLINE_FACTOR_CAPACITY> factors;

        //! Number of active canonical factors.
        std::size_t size;
    };

    /*! \brief Construct directly from canonical factor storage. */
    explicit IBKernelTensorProduct(CanonicalFactors factors);

    /*! \brief Validate and canonicalize one or two supplied factors. */
    static CanonicalFactors canonicalize(std::initializer_list<IBKernel> factors);

    /*! \brief Parse a validated legacy scalar or composite name. */
    static CanonicalFactors parse_name(const std::string& name);

    //! Inline canonical factors.
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

//! Diagnostic output only; numerical dispatch uses identities, not names.
std::ostream& operator<<(std::ostream& stream, const IBKernelTensorProduct& kernel);
} // namespace IBTK

#endif
