// ---------------------------------------------------------------------
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
// This file is part of IBAMR and is distributed under the 3-clause BSD
// license. The full text of the license can be found in COPYRIGHT.
// ---------------------------------------------------------------------

#ifndef included_IBTK_IBKernelConcepts
#define included_IBTK_IBKernelConcepts

#include <ibtk/config.h>

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace IBTK
{
namespace IBKernels
{
/*! \brief Owning storage for the supplied IB kernel coefficients. */
template <class T, std::size_t N>
using Weights = std::array<T, N>;
} // namespace IBKernels

/*!
 * \brief Coefficient type and extent of an owning fixed-size weight container.
 *
 * Specializations supply value_type and a constant extent. weights[i] supplies
 * coefficient i; its storage layout is unspecified. Copies must have independent
 * coefficient storage.
 */
template <class T>
struct IBKernelWeightsTraits;

/*! \brief Coefficient type and extent of std::array. */
template <class T, std::size_t N>
struct IBKernelWeightsTraits<std::array<T, N>>
{
    using value_type = T;
    static constexpr std::size_t extent = N;
};

/*! \brief An owning, nonempty indexed collection of floating-point coefficients. */
template <class T>
concept IBKernelWeights = requires(const T& weights, std::size_t i)
{
    typename IBKernelWeightsTraits<T>::value_type;
    requires std::floating_point<typename IBKernelWeightsTraits<T>::value_type>;
    typename std::integral_constant<std::size_t, IBKernelWeightsTraits<T>::extent>;
    requires(IBKernelWeightsTraits<T>::extent > 0);
    {
        weights[i]
    } -> std::convertible_to<typename IBKernelWeightsTraits<T>::value_type>;
};

/*!
 * \brief A one-dimensional IB kernel evaluator.
 *
 * kernel(r), with r a const double lvalue, returns N weights by value,
 * with entry i equal to phi(r-i).
 * The argument r is the displacement from the first stencil point to the
 * evaluation point, divided by grid spacing. For odd N, (N-2)/2 <= r < N/2;
 * for even N, N/2-1 <= r <= N/2. No grid-spacing factors are included.
 * Evaluation leaves the evaluator unchanged and returns the same coefficients
 * for the same r.
 */
template <class T>
concept IBKernelEvaluatorScalar =
    std::regular_invocable<const T&, const double&> && IBKernelWeights<std::invoke_result_t<const T&, const double&>>;

namespace detail
{
/*! \brief Check positive widths and their product against the returned extent. */
template <class T, int Axis>
constexpr bool ib_kernel_cartesian_shape();

/*! \brief Cartesian evaluator requirements for one coordinate axis. */
template <class T, int Axis>
concept IBKernelEvaluatorCartesianAxis = requires(const T& kernel, const std::array<double, NDIM>& r)
{
    {
        T::template get_stencil_widths<Axis>()
    } -> std::same_as<std::array<std::size_t, NDIM>>;
    typename std::integral_constant<std::array<std::size_t, NDIM>, T::template get_stencil_widths<Axis>()>;
    {
        kernel.template evaluate<Axis>(r)
    } -> IBKernelWeights;
    requires(ib_kernel_cartesian_shape<T, Axis>());
};
} // namespace detail

/*!
 * \brief An IB kernel evaluator on a rectangular Cartesian stencil.
 *
 * get_stencil_widths<Axis>() gives positive compile-time widths.
 * evaluate<Axis>(r) returns coefficients by value, with coordinate zero varying
 * fastest and extent equal to the product of the widths. Axis selects the
 * distinguished coordinate: x (0), y (1), or z (2).
 * Each r[d] follows the IBKernelEvaluatorScalar convention for that width.
 * Evaluation leaves the evaluator unchanged and returns the same coefficients
 * for the same arguments. The coefficients need not form a separable product.
 */
template <class T>
concept IBKernelEvaluatorCartesian =
    detail::IBKernelEvaluatorCartesianAxis<T, 0> && detail::IBKernelEvaluatorCartesianAxis<T, 1>
#if (NDIM == 3)
    && detail::IBKernelEvaluatorCartesianAxis<T, 2>
#endif
    ;
} // namespace IBTK
#include <ibtk/private/IBKernelConcepts-inl.h>
#endif
