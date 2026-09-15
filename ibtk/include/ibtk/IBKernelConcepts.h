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

namespace detail
{
/*! \brief Storage operations used by the supplied evaluators. */
template <class T, std::size_t N>
concept IBKernelWritableWeights =
    IBKernelWeights<T> && std::default_initializable<T> && std::move_constructible<T> &&
    (IBKernelWeightsTraits<T>::extent == N) &&
    requires(T & weights, std::size_t i, typename IBKernelWeightsTraits<T>::value_type value)
{
    weights[i] = value;
};

/*! \brief Positive compile-time width of a scalar evaluator object. */
template <class T>
concept IBKernelScalarShape = std::is_object_v<T> && requires
{
    {
        T::get_stencil_width()
    } -> std::same_as<std::size_t>;
    typename std::integral_constant<std::size_t, T::get_stencil_width()>;
    requires(T::get_stencil_width() > 0);
};

/*! \brief Return the stencil product, or zero for invalid or overflowing widths. */
template <class T, int Axis>
constexpr std::size_t ib_kernel_stencil_size();

/*! \brief Positive compile-time Cartesian widths with a representable product. */
template <class T, int Axis>
concept IBKernelCartesianShape = (Axis >= 0 && Axis < NDIM) && requires
{
    {
        T::template get_stencil_widths<Axis>()
    } -> std::same_as<std::array<std::size_t, NDIM>>;
    typename std::integral_constant<std::array<std::size_t, NDIM>, T::template get_stencil_widths<Axis>()>;
    requires(ib_kernel_stencil_size<T, Axis>() > 0);
};
} // namespace detail

/*!
 * \brief A scalar IB kernel evaluated with Input displacements and Coefficient weights.
 *
 * get_stencil_width() gives the positive compile-time width N.
 * evaluate<Output>(r) returns N independently owned weights, with entry i equal
 * to phi(r-i). Output supplies the coefficient type and storage.
 * The argument r is the displacement from the first stencil point to the
 * evaluation point, divided by grid spacing. For odd N, (N-2)/2 <= r < N/2;
 * for even N, N/2-1 <= r <= N/2. No grid-spacing factors are included.
 * Evaluation leaves the input and evaluator unchanged and initializes every
 * coefficient deterministically. This concept checks the selected types with
 * IBKernels::Weights storage; other output containers are checked at their call.
 */
template <class T, class Input = double, class Coefficient = double>
concept IBKernelEvaluatorScalar = detail::IBKernelScalarShape<T> && std::floating_point<Input> &&
                                  std::floating_point<Coefficient> && requires(const T& kernel, const Input& r)
{
    {
        kernel.template evaluate<IBKernels::Weights<Coefficient, T::get_stencil_width()>>(r)
    } -> std::same_as<IBKernels::Weights<Coefficient, T::get_stencil_width()>>;
};

namespace detail
{
/*! \brief Cartesian evaluation with the selected input and coefficient types. */
template <class T, int Axis, class Input, class Coefficient>
concept IBKernelEvaluatorCartesianAxis =
    IBKernelCartesianShape<T, Axis> && requires(const T& kernel, const std::array<Input, NDIM>& r)
{
    {
        kernel.template evaluate<Axis, IBKernels::Weights<Coefficient, ib_kernel_stencil_size<T, Axis>()>>(r)
    } -> std::same_as<IBKernels::Weights<Coefficient, ib_kernel_stencil_size<T, Axis>()>>;
};
} // namespace detail

/*!
 * \brief An IB kernel evaluator on a rectangular Cartesian stencil.
 *
 * get_stencil_widths<Axis>() gives positive compile-time widths.
 * evaluate<Axis, Output>(r) returns independently owned coefficients, with
 * coordinate zero varying fastest and extent equal to the product of the widths.
 * Axis selects the distinguished coordinate: x (0), y (1), or z (2).
 * Each r[d] follows the IBKernelEvaluatorScalar convention for that width.
 * Evaluation leaves the input and evaluator unchanged and initializes every
 * coefficient deterministically. The coefficients need not form a separable
 * product. This concept checks Input and Coefficient with IBKernels::Weights
 * storage; other output containers are checked at their call.
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
