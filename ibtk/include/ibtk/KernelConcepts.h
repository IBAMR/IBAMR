// ---------------------------------------------------------------------
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
// This file is part of IBAMR and is distributed under the 3-clause BSD
// license. The full text of the license can be found in COPYRIGHT.
// ---------------------------------------------------------------------

#ifndef included_IBTK_KernelConcepts
#define included_IBTK_KernelConcepts

#include <ibtk/config.h>

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace IBTK
{
/*!
 * \brief Describe an owning fixed-size array of kernel coefficients.
 *
 * Specializations supply value_type and a constant extent. The array must own
 * its coefficients: copying it must not alias another array's storage.
 */
template <class T>
struct KernelWeightTraits;

/*! \brief Coefficient type and extent of std::array. */
template <class T, std::size_t N>
struct KernelWeightTraits<std::array<T, N>>
{
    using value_type = T;
    static constexpr std::size_t extent = N;
};

/*! \brief An owning, nonempty array of floating-point coefficients. */
template <class T>
concept KernelWeights = requires(const T& weights, std::size_t i)
{
    typename KernelWeightTraits<T>::value_type;
    requires std::floating_point<typename KernelWeightTraits<T>::value_type>;
    typename std::integral_constant<std::size_t, KernelWeightTraits<T>::extent>;
    requires(KernelWeightTraits<T>::extent > 0);
    {
        weights[i]
    } -> std::convertible_to<typename KernelWeightTraits<T>::value_type>;
};

/*!
 * \brief A const-callable scalar kernel returning its Weights by value.
 *
 * Entry i is phi(r-i), where r is the displacement from the first grid point
 * to the interpolation point, divided by grid spacing. See Kernels for the
 * first-point convention used by matrix construction.
 */
template <class T>
concept ScalarKernel = requires(const T& kernel, double r)
{
    typename T::Weights;
    requires KernelWeights<typename T::Weights>;
    {
        kernel(r)
    } -> std::same_as<typename T::Weights>;
};

namespace detail
{
/*! \brief Check positive widths and their product against the returned extent. */
template <class T, int Axis>
constexpr bool tensor_kernel_shape();
} // namespace detail
} // namespace IBTK

#include <ibtk/private/KernelConcepts-inl.h>

namespace IBTK
{
/*! \brief Tensor evaluator requirements for one coordinate axis. */
template <class T, int Axis>
concept TensorKernelAxis = requires(const T& kernel, const std::array<double, NDIM>& r)
{
    {
        T::template get_stencil_widths<Axis>()
    } -> std::same_as<std::array<std::size_t, NDIM>>;
    typename std::integral_constant<std::array<std::size_t, NDIM>, T::template get_stencil_widths<Axis>()>;
    {
        kernel.template evaluate<Axis>(r)
    } -> KernelWeights;
    requires(detail::tensor_kernel_shape<T, Axis>());
};

/*!
 * \brief A tensor evaluator for every coordinate axis in NDIM dimensions.
 *
 * Widths are positive compile-time constants. evaluate<Axis>(r) returns owned
 * coefficients, coordinate zero varying fastest, with extent equal to the
 * product of the widths. Axis is x (0), y (1), or z (2).
 * Each r[d] is the displacement from the first stencil point in direction d,
 * divided by grid spacing; see Kernels for the first-point convention.
 * For side-centered interpolation, Axis is the velocity component and selects
 * the coordinate with side-normal centering.
 */
template <class T>
concept TensorKernel = TensorKernelAxis<T, 0> && TensorKernelAxis<T, 1>
#if (NDIM == 3)
                       && TensorKernelAxis<T, 2>
#endif
    ;
} // namespace IBTK
#endif
