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

#ifndef included_IBTK_ib_kernels_inl
#define included_IBTK_ib_kernels_inl

#include <ibtk/config.h>

#include <ibtk/ib_kernels.h>

#include <cmath>

namespace IBTK
{
template <std::size_t N>
requires(N > 0) template <detail::IBKernelWritableWeights<N> Output, std::floating_point Input>
inline Output IBKernels::BSpline<N>::evaluate(const Input r) const
{
    using Coefficient = typename IBKernelWeightsTraits<Output>::value_type;
    const Coefficient x = r;
    // Unit-spacing form of the B-spline basis recurrence; see C. de Boor,
    // "On calculating with B-splines", J. Approx. Theory 6 (1972), 50-62,
    // Section 2, Eq. (27), doi:10.1016/0021-9045(72)90080-9.
    //
    // With w_0^(0) = 1 and missing entries zero, the degree-d weights satisfy
    // w_i^(d) = ((i+1-t)*w_i^(d-1) + (t+d-i)*w_(i-1)^(d-1))/d.
    // Here t is in [0,1]; degree N-1 gives the N weights in stencil order.
    // Each old weight is split into nonnegative contributions whose
    // coefficients sum to one. This avoids cancellation between weights
    // and preserves their sum in exact arithmetic.
    const Coefficient width = N;
    const Coefficient t = x - Coefficient{ 0.5 } * (width - 2);
    Output w{};
    w[0] = 1;
    for (std::size_t degree = 1; degree < N; ++degree)
    {
        const Coefficient inverse_degree = Coefficient{ 1 } / degree;
        Coefficient saved = 0;
        for (std::size_t i = 0; i < degree; ++i)
        {
            const Coefficient term = w[i] * inverse_degree;
            w[i] = saved + ((i + 1) - t) * term;
            saved = (t + (degree - i - 1)) * term;
        }
        w[degree] = saved;
    }
    return w;
}

template <detail::IBKernelWritableWeights<3> Output, std::floating_point Input>
inline Output
IBKernels::IB3::evaluate(const Input r) const
{
    using Coefficient = typename IBKernelWeightsTraits<Output>::value_type;
    const Coefficient x = r;
    const Coefficient s = x - 1;
    const Coefficient q = std::sqrt(1 - 3 * s * s);
    Output w{};
    w[0] = (2 - 3 * s - q) / 6;
    w[1] = (1 + q) / 3;
    w[2] = (2 + 3 * s - q) / 6;
    return w;
}

template <detail::IBKernelWritableWeights<4> Output, std::floating_point Input>
inline Output
IBKernels::IB4::evaluate(const Input r) const
{
    using Coefficient = typename IBKernelWeightsTraits<Output>::value_type;
    const Coefficient x = r;
    Output w;
    // Use kernel symmetry and moment conditions to compute all four weights
    // from one square root.
    const Coefficient r0 = x - 1;
    const Coefficient q = std::sqrt(1 + 4 * r0 * (1 - r0));
    w[0] = Coefficient{ 0.125 } * (3 - 2 * r0 - q);
    w[1] = Coefficient{ 0.125 } * (3 - 2 * r0 + q);
    w[2] = Coefficient{ 0.125 } * (1 + 2 * r0 + q);
    w[3] = Coefficient{ 0.125 } * (1 + 2 * r0 - q);
    return w;
}

template <detail::IBKernelWritableWeights<5> Output, std::floating_point Input>
inline Output
IBKernels::IB5::evaluate(const Input r) const
{
    using Coefficient = typename IBKernelWeightsTraits<Output>::value_type;
    const Coefficient x = r;
    static const Coefficient K = (38 - std::sqrt(Coefficient{ 69 })) / 60;
    const Coefficient r0 = x - 2;
    const Coefficient r2 = r0 * r0;
    const Coefficient r3 = r2 * r0;
    const Coefficient r4 = r2 * r2;
    const Coefficient r6 = r4 * r2;
    const Coefficient phi =
        (136 - 40 * K - 40 * r2 +
         std::sqrt(Coefficient{ 2 }) * std::sqrt(3123 - 6840 * K + 3600 * K * K - 12440 * r2 + 25680 * K * r2 -
                                                 12600 * K * K * r2 + 8080 * r4 - 8400 * K * r4 - 1400 * r6)) /
        280;
    Output w{};
    w[0] = (Coefficient{ 1 } / 12) * (-2 + 2 * phi + 2 * K + r0 - 3 * K * r0 + 2 * r2 - r3);
    w[1] = (Coefficient{ 1 } / 6) * (4 - 4 * phi - K - 4 * r0 + 3 * K * r0 - r2 + r3);
    w[2] = phi;
    w[3] = (Coefficient{ 1 } / 6) * (4 - 4 * phi - K + 4 * r0 - 3 * K * r0 - r2 - r3);
    w[4] = (Coefficient{ 1 } / 12) * (-2 + 2 * phi + 2 * K - r0 + 3 * K * r0 + 2 * r2 + r3);
    return w;
}

template <detail::IBKernelWritableWeights<6> Output, std::floating_point Input>
inline Output
IBKernels::IB6::evaluate(const Input r) const
{
    using Coefficient = typename IBKernelWeightsTraits<Output>::value_type;
    const Coefficient x = r;
    Output w;
    const Coefficient rl = 3 - x;
    const Coefficient r2 = rl * rl;
    const Coefficient r3 = r2 * rl;
    const Coefficient r4 = r3 * rl;
    const Coefficient r5 = r4 * rl;
    static const Coefficient K = (Coefficient{ 59 } / 60) * (1 - std::sqrt(1 - (Coefficient{ 3220 } / 3481)));
    static const Coefficient K2 = K * K;
    static const Coefficient alpha = 28;
    const Coefficient beta = (Coefficient{ 9 } / 4) - (Coefficient{ 3 } / 2) * (K + r2) +
                             ((Coefficient{ 22 } / 3) - 7 * K) * rl - (Coefficient{ 7 } / 3) * r3;
    const Coefficient gamma =
        (Coefficient{ 1 } / 4) *
        (((Coefficient{ 161 } / 36) - (Coefficient{ 59 } / 6) * K + 5 * K2) * (Coefficient{ 1 } / 2) * r2 +
         (-(Coefficient{ 109 } / 24) + 5 * K) * (Coefficient{ 1 } / 3) * r4 + (Coefficient{ 5 } / 18) * r5 * rl);
    const Coefficient discr = beta * beta - 4 * alpha * gamma;
    w[0] = (-beta + std::copysign(Coefficient{ 1 }, (Coefficient{ 3 } / 2) - K) * std::sqrt(discr)) / (2 * alpha);
    w[1] = -3 * w[0] - (Coefficient{ 1 } / 16) + (Coefficient{ 1 } / 8) * (K + r2) +
           (Coefficient{ 1 } / 12) * (3 * K - 1) * rl + (Coefficient{ 1 } / 12) * r3;
    w[2] = 2 * w[0] + (Coefficient{ 1 } / 4) + (Coefficient{ 1 } / 6) * (4 - 3 * K) * rl - (Coefficient{ 1 } / 6) * r3;
    w[3] = 2 * w[0] + (Coefficient{ 5 } / 8) - (Coefficient{ 1 } / 4) * (K + r2);
    w[4] = -3 * w[0] + (Coefficient{ 1 } / 4) - (Coefficient{ 1 } / 6) * (4 - 3 * K) * rl + (Coefficient{ 1 } / 6) * r3;
    w[5] = w[0] - (Coefficient{ 1 } / 16) + (Coefficient{ 1 } / 8) * (K + r2) -
           (Coefficient{ 1 } / 12) * (3 * K - 1) * rl - (Coefficient{ 1 } / 12) * r3;
    return w;
}

template <std::size_t N>
requires(N > 0) constexpr std::size_t IBKernels::BSpline<N>::get_stencil_width()
{
    return N;
}

inline constexpr std::size_t
IBKernels::IB3::get_stencil_width()
{
    return 3;
}

inline constexpr std::size_t
IBKernels::IB4::get_stencil_width()
{
    return 4;
}

inline constexpr std::size_t
IBKernels::IB5::get_stencil_width()
{
    return 5;
}

inline constexpr std::size_t
IBKernels::IB6::get_stencil_width()
{
    return 6;
}

} // namespace IBTK

#endif
