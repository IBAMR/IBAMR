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

#ifndef included_IBTK_IBKernelEvaluators_inl
#define included_IBTK_IBKernelEvaluators_inl

#include <ibtk/config.h>

#include <ibtk/IBKernelEvaluators.h>

#include <cmath>

namespace IBTK
{
template <std::size_t N>
inline typename IBKernelEvaluatorBSpline<N>::Weights
IBKernelEvaluatorBSpline<N>::operator()(const double r) const
{
    // On this knot interval, t is in [0, 1]. The recurrence combines
    // nonnegative contributions without subtracting nearly equal weights.
    const double t = r - 0.5 * (static_cast<double>(N) - 2.0);
    Weights w = {};
    w[0] = 1.0;
    for (std::size_t degree = 1; degree < N; ++degree)
    {
        double saved = 0.0;
        for (std::size_t i = 0; i < degree; ++i)
        {
            const double term = w[i] / static_cast<double>(degree);
            w[i] = saved + (static_cast<double>(i + 1) - t) * term;
            saved = (t + static_cast<double>(degree - i - 1)) * term;
        }
        w[degree] = saved;
    }
    return w;
}

inline IBKernelEvaluatorIB3::Weights
IBKernelEvaluatorIB3::operator()(const double r) const
{
    const double s = r - 1.0;
    const double q = std::sqrt(1.0 - 3.0 * s * s);
    return { (2.0 - 3.0 * s - q) / 6.0, (1.0 + q) / 3.0, (2.0 + 3.0 * s - q) / 6.0 };
}

inline IBKernelEvaluatorIB4::Weights
IBKernelEvaluatorIB4::operator()(const double r) const
{
    Weights w;
    // Match the specialized Fortran recurrence: IB4 symmetry and moment
    // conditions provide all four weights from one square root, avoiding four
    // generic pointwise kernel evaluations.
    const double r0 = r - 1.0;
    const double q = std::sqrt(1.0 + 4.0 * r0 * (1.0 - r0));
    w[0] = 0.125 * (3.0 - 2.0 * r0 - q);
    w[1] = 0.125 * (3.0 - 2.0 * r0 + q);
    w[2] = 0.125 * (1.0 + 2.0 * r0 + q);
    w[3] = 0.125 * (1.0 + 2.0 * r0 - q);
    return w;
}

inline IBKernelEvaluatorIB5::Weights
IBKernelEvaluatorIB5::operator()(const double r) const
{
    static const double K = (38.0 - std::sqrt(69.0)) / 60.0;
    const double r0 = r - 2.0;
    const double r2 = r0 * r0;
    const double r3 = r2 * r0;
    const double r4 = r2 * r2;
    const double r6 = r4 * r2;
    const double phi =
        (136.0 - 40.0 * K - 40.0 * r2 +
         std::sqrt(2.0) * std::sqrt(3123.0 - 6840.0 * K + 3600.0 * K * K - 12440.0 * r2 + 25680.0 * K * r2 -
                                    12600.0 * K * K * r2 + 8080.0 * r4 - 8400.0 * K * r4 - 1400.0 * r6)) /
        280.0;
    return {
        (1.0 / 12.0) * (-2.0 + 2.0 * phi + 2.0 * K + r0 - 3.0 * K * r0 + 2.0 * r2 - r3),
        (1.0 / 6.0) * (4.0 - 4.0 * phi - K - 4.0 * r0 + 3.0 * K * r0 - r2 + r3),
        phi,
        (1.0 / 6.0) * (4.0 - 4.0 * phi - K + 4.0 * r0 - 3.0 * K * r0 - r2 - r3),
        (1.0 / 12.0) * (-2.0 + 2.0 * phi + 2.0 * K - r0 + 3.0 * K * r0 + 2.0 * r2 + r3),
    };
}

inline IBKernelEvaluatorIB6::Weights
IBKernelEvaluatorIB6::operator()(const double r) const
{
    Weights w;
    const double rl = 3.0 - r;
    const double r2 = rl * rl;
    const double r3 = r2 * rl;
    const double r4 = r3 * rl;
    const double r5 = r4 * rl;
    static const double K = (59.0 / 60.0) * (1.0 - std::sqrt(1.0 - (3220.0 / 3481.0)));
    static const double K2 = K * K;
    static const double alpha = 28.0;
    const double beta = (9.0 / 4.0) - (3.0 / 2.0) * (K + r2) + ((22.0 / 3.0) - 7.0 * K) * rl - (7.0 / 3.0) * r3;
    const double gamma = (1.0 / 4.0) * (((161.0 / 36.0) - (59.0 / 6.0) * K + 5.0 * K2) * (1.0 / 2.0) * r2 +
                                        (-(109.0 / 24.0) + 5.0 * K) * (1.0 / 3.0) * r4 + (5.0 / 18.0) * r5 * rl);
    const double discr = beta * beta - 4.0 * alpha * gamma;
    w[0] = (-beta + std::copysign(1.0, (3.0 / 2.0) - K) * std::sqrt(discr)) / (2.0 * alpha);
    w[1] =
        -3.0 * w[0] - (1.0 / 16.0) + (1.0 / 8.0) * (K + r2) + (1.0 / 12.0) * (3.0 * K - 1.0) * rl + (1.0 / 12.0) * r3;
    w[2] = 2.0 * w[0] + (1.0 / 4.0) + (1.0 / 6.0) * (4.0 - 3.0 * K) * rl - (1.0 / 6.0) * r3;
    w[3] = 2.0 * w[0] + (5.0 / 8.0) - (1.0 / 4.0) * (K + r2);
    w[4] = -3.0 * w[0] + (1.0 / 4.0) - (1.0 / 6.0) * (4.0 - 3.0 * K) * rl + (1.0 / 6.0) * r3;
    w[5] = w[0] - (1.0 / 16.0) + (1.0 / 8.0) * (K + r2) - (1.0 / 12.0) * (3.0 * K - 1.0) * rl - (1.0 / 12.0) * r3;
    return w;
}

} // namespace IBTK

#endif
