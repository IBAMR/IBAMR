// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_IBKernelConcepts
#define included_IBTK_IBKernelConcepts

#include <ibtk/config.h>

#include <array>
#include <concepts>
#include <cstddef>

namespace IBTK
{
namespace IBKernelEvaluators
{
/*!
 * \brief Owning storage for the supplied IB kernel coefficients.
 *
 * This is the default output container of the evaluators in this namespace.
 * Any other container that models IBKernelWeights can be used by specializing
 * IBKernelWeightsTraits for it.
 */
template <class T, std::size_t N>
using Weights = std::array<T, N>;
} // namespace IBKernelEvaluators

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

/*! \brief typename IBKernelWeightsTraits<T>::value_type. */
template <class T>
using ib_kernel_weights_value_t = typename IBKernelWeightsTraits<T>::value_type;

/*! \brief IBKernelWeightsTraits<T>::extent. */
template <class T>
inline constexpr std::size_t ib_kernel_weights_extent_v = IBKernelWeightsTraits<T>::extent;

/*!
 * \brief An owning, nonempty indexed collection of floating-point coefficients.
 *
 * A type T models this concept when:
 * - IBKernelWeightsTraits<T> supplies a floating-point value_type and a
 *   positive constant extent;
 * - weights[i] is readable and convertible to value_type.
 */
template <class T>
concept IBKernelWeights = requires(const T& weights, std::size_t i)
{
    requires std::floating_point<ib_kernel_weights_value_t<T>>;
    requires(ib_kernel_weights_extent_v<T> > 0);
    {
        weights[i]
    } -> std::convertible_to<ib_kernel_weights_value_t<T>>;
};

namespace detail
{
/*! \brief Storage operations used by the supplied evaluators. */
template <class T, std::size_t N>
concept IBKernelWritableWeights =
    IBKernelWeights<T> && std::default_initializable<T> && std::move_constructible<T> &&
    (ib_kernel_weights_extent_v<T> == N) && requires(T & weights, std::size_t i, ib_kernel_weights_value_t<T> value)
{
    weights[i] = value;
};

/*! \brief Return the product of the compile-time Cartesian stencil widths. */
template <class T, int Axis>
constexpr std::size_t ib_kernel_stencil_size();
} // namespace detail

/*!
 * \brief A one-dimensional IB kernel phi, evaluated over its stencil at a
 * single displacement.
 *
 * A type T models this concept when:
 * - T::get_stencil_width() gives a positive compile-time width N;
 * - T::evaluate<Output>(r) returns N independently owned weights, with entry
 *   i equal to phi(r - i), for Output equal to IBKernelEvaluators::Weights<Coefficient, N>;
 * - evaluation leaves r and the evaluator unchanged and initializes every
 *   coefficient deterministically.
 *
 * Coordinate convention: r is the displacement, in grid spacings, from the
 * first stencil point to the evaluation point. For odd N, (N-2)/2 <= r < N/2;
 * for even N, N/2-1 <= r <= N/2. No grid-spacing factors are included.
 *
 * Example: a user-defined kernel in a tensor product with IB4
 * (\see IBKernelEvaluatorTensorProduct):
 * \code
 * struct MyKernel
 * {
 *     static constexpr std::size_t get_stencil_width() { return 2; }
 *
 *     template <IBTK::IBKernelWeights Output, std::floating_point Input>
 *     Output evaluate(Input r) const
 *     {
 *         using Coefficient = IBTK::ib_kernel_weights_value_t<Output>;
 *         Output w{};
 *         w[0] = Coefficient{ 1 } - r;
 *         w[1] = r;
 *         return w;
 *     }
 * };
 * static_assert(IBTK::IBKernelEvaluatorScalar<MyKernel>);
 * const IBTK::IBKernelEvaluatorTensorProduct product{ MyKernel{}, IBTK::IBKernelEvaluators::IB4{} };
 * \endcode
 */
template <class T, class Input = double, class Coefficient = double>
concept IBKernelEvaluatorScalar = std::floating_point<Input> && std::floating_point<Coefficient> &&
                                  (T::get_stencil_width() > 0) && requires(const T& kernel, const Input& r)
{
    {
        kernel.template evaluate<IBKernelEvaluators::Weights<Coefficient, T::get_stencil_width()>>(r)
    } -> std::same_as<IBKernelEvaluators::Weights<Coefficient, T::get_stencil_width()>>;
};

namespace detail
{
/*! \brief Cartesian evaluation with the selected input and coefficient types. */
template <class T, int Axis, class Input, class Coefficient>
concept IBKernelEvaluatorCartesianAxis = (Axis >= 0 && Axis < NDIM) && requires
{
    T::template get_stencil_widths<Axis>();
}
&&(ib_kernel_stencil_size<T, Axis>() > 0) && requires(const T& kernel, const std::array<Input, NDIM>& r)
{
    {
        kernel.template evaluate<Axis, IBKernelEvaluators::Weights<Coefficient, ib_kernel_stencil_size<T, Axis>()>>(r)
    } -> std::same_as<IBKernelEvaluators::Weights<Coefficient, ib_kernel_stencil_size<T, Axis>()>>;
};
} // namespace detail

/*!
 * \brief An IB kernel evaluator on a rectangular Cartesian stencil, not
 * necessarily a separable product of scalar kernels.
 *
 * A type T models this concept when, for Axis equal to each coordinate 0, ...,
 * NDIM - 1:
 * - T::get_stencil_widths<Axis>() gives positive compile-time widths, one per
 *   coordinate;
 * - T::evaluate<Axis, Output>(r) returns the coefficients over that stencil,
 *   independently owned, with coordinate zero varying fastest and extent
 *   equal to the product of the widths, for Output equal to
 *   IBKernelEvaluators::Weights<Coefficient, N>;
 * - each r[d] follows the IBKernelEvaluatorScalar coordinate convention for
 *   that coordinate's width;
 * - evaluation leaves r and the evaluator unchanged and initializes every
 *   coefficient deterministically.
 *
 * \see IBKernelEvaluatorScalar for the coordinate convention.
 */
template <class T, class Input = double, class Coefficient = double>
concept IBKernelEvaluatorCartesian = std::floating_point<Input> && std::floating_point<Coefficient> &&
                                     detail::IBKernelEvaluatorCartesianAxis<T, 0, Input, Coefficient> &&
                                     detail::IBKernelEvaluatorCartesianAxis<T, 1, Input, Coefficient>
#if (NDIM == 3)
                                     && detail::IBKernelEvaluatorCartesianAxis<T, 2, Input, Coefficient>
#endif
    ;
} // namespace IBTK
#include <ibtk/private/IBKernelConcepts-inl.h>
#endif
