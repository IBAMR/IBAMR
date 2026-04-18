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

#ifndef included_IBTK_PRIVATE_PETScMatUtilities_inl
#define included_IBTK_PRIVATE_PETScMatUtilities_inl

#include <ibtk/config.h>

#include <ibtk/PETScMatUtilities.h>

#include <algorithm>
#include <cmath>

namespace IBTK
{
inline void
PETScMatUtilities::piecewise_linear_delta_fcn(const double r, double* const w)
{
    w[0] = 1.0 - r;
    w[1] = r;
    return;
}

inline void
PETScMatUtilities::bspline_3_delta_fcn(const double r, double* const w)
{
    const double rr[4] = { r, r - 1.0, r - 2.0, r - 3.0 };
    for (int i = 0; i < 4; ++i)
    {
        const double x = std::abs(rr[i]);
        const double rp = x + 1.5;
        if (x <= 0.5)
        {
            w[i] = 0.5 * (-2.0 * rp * rp + 6.0 * rp - 3.0);
        }
        else if (x <= 1.5)
        {
            w[i] = 0.5 * (rp * rp - 6.0 * rp + 9.0);
        }
        else
        {
            w[i] = 0.0;
        }
    }
    return;
}

inline void
PETScMatUtilities::bspline_4_delta_fcn(const double r, double* const w)
{
    const double rr[4] = { r, r - 1.0, r - 2.0, r - 3.0 };
    for (int i = 0; i < 4; ++i)
    {
        const double x = std::abs(rr[i]);
        const double rp = x + 2.0;
        if (x <= 1.0)
        {
            w[i] = (1.0 / 6.0) * (3.0 * std::pow(rp, 3) - 24.0 * std::pow(rp, 2) + 60.0 * rp - 44.0);
        }
        else if (x <= 2.0)
        {
            w[i] = (1.0 / 6.0) * (-std::pow(rp, 3) + 12.0 * std::pow(rp, 2) - 48.0 * rp + 64.0);
        }
        else
        {
            w[i] = 0.0;
        }
    }
    return;
}

inline void
PETScMatUtilities::bspline_5_delta_fcn(const double r, double* const w)
{
    const double rr[6] = { r, r - 1.0, r - 2.0, r - 3.0, r - 4.0, r - 5.0 };
    for (int i = 0; i < 6; ++i)
    {
        const double x = std::abs(rr[i]);
        const double rp = x + 2.5;
        if (x <= 0.5)
        {
            w[i] = (1.0 / 24.0) *
                   (6.0 * std::pow(rp, 4) - 60.0 * std::pow(rp, 3) + 210.0 * std::pow(rp, 2) - 300.0 * rp + 155.0);
        }
        else if (x <= 1.5)
        {
            w[i] = (1.0 / 24.0) *
                   (-4.0 * std::pow(rp, 4) + 60.0 * std::pow(rp, 3) - 330.0 * std::pow(rp, 2) + 780.0 * rp - 655.0);
        }
        else if (x <= 2.5)
        {
            w[i] = (1.0 / 24.0) *
                   (std::pow(rp, 4) - 20.0 * std::pow(rp, 3) + 150.0 * std::pow(rp, 2) - 500.0 * rp + 625.0);
        }
        else
        {
            w[i] = 0.0;
        }
    }
    return;
}

inline void
PETScMatUtilities::bspline_6_delta_fcn(const double r, double* const w)
{
    const double rr[6] = { r, r - 1.0, r - 2.0, r - 3.0, r - 4.0, r - 5.0 };
    for (int i = 0; i < 6; ++i)
    {
        const double x = std::abs(rr[i]);
        const double rp = x + 3.0;
        if (x <= 1.0)
        {
            w[i] = (1.0 / 60.0) * (2193.0 - 3465.0 * rp + 2130.0 * std::pow(rp, 2) - 630.0 * std::pow(rp, 3) +
                                   90.0 * std::pow(rp, 4) - 5.0 * std::pow(rp, 5));
        }
        else if (x <= 2.0)
        {
            w[i] = (1.0 / 120.0) * (-10974.0 + 12270.0 * rp - 5340.0 * std::pow(rp, 2) + 1140.0 * std::pow(rp, 3) -
                                    120.0 * std::pow(rp, 4) + 5.0 * std::pow(rp, 5));
        }
        else if (x <= 3.0)
        {
            w[i] = (1.0 / 120.0) * (7776.0 - 6480.0 * rp + 2160.0 * std::pow(rp, 2) - 360.0 * std::pow(rp, 3) +
                                    30.0 * std::pow(rp, 4) - std::pow(rp, 5));
        }
        else
        {
            w[i] = 0.0;
        }
    }
    return;
}

inline void
PETScMatUtilities::ib_3_delta_fcn(const double r, double* const w)
{
    const double rr[4] = { r, r - 1.0, r - 2.0, r - 3.0 };
    for (int i = 0; i < 4; ++i)
    {
        const double x = std::abs(rr[i]);
        if (x < 0.5)
        {
            w[i] = (1.0 / 3.0) * (1.0 + std::sqrt(1.0 - 3.0 * x * x));
        }
        else if (x < 1.5)
        {
            const double t = 1.0 - x;
            w[i] = (1.0 / 6.0) * (5.0 - 3.0 * x - std::sqrt(1.0 - 3.0 * t * t));
        }
        else
        {
            w[i] = 0.0;
        }
    }
    return;
}

inline void
PETScMatUtilities::ib_4_delta_fcn(const double r, double* const w)
{
    const double q = std::sqrt(-7.0 + 12.0 * r - 4.0 * r * r);
    w[0] = 0.125 * (5.0 - 2.0 * r - q);
    w[1] = 0.125 * (5.0 - 2.0 * r + q);
    w[2] = 0.125 * (-1.0 + 2.0 * r + q);
    w[3] = 0.125 * (-1.0 + 2.0 * r - q);
    return;
}

inline void
PETScMatUtilities::ib_5_delta_fcn(const double r, double* const w)
{
    static const double K = (38.0 - std::sqrt(69.0)) / 60.0;
    const double rr[6] = { r, r - 1.0, r - 2.0, r - 3.0, r - 4.0, r - 5.0 };
    for (int i = 0; i < 6; ++i)
    {
        const double x = std::abs(rr[i]);
        if (x > 2.5)
        {
            w[i] = 0.0;
            continue;
        }
        const double x2 = x * x;
        const double x4 = x2 * x2;
        const double x6 = x4 * x2;
        const double radicand = std::max(0.0,
                                         2.0 * (3123.0 - 6840.0 * K + 3600.0 * K * K - 12440.0 * x2 + 25680.0 * K * x2 -
                                                12600.0 * K * K * x2 + 8080.0 * x4 - 8400.0 * K * x4 - 1400.0 * x6));
        const double phi = (136.0 - 40.0 * K - 40.0 * x2 + std::sqrt(radicand)) / 280.0;
        if (x <= 0.5)
        {
            w[i] = phi;
        }
        else if (x <= 1.5)
        {
            const double xr = x - 1.0;
            w[i] = (4.0 - 4.0 * phi - K - 4.0 * xr + 3.0 * K * xr - xr * xr + xr * xr * xr) / 6.0;
        }
        else
        {
            const double xr = x - 2.0;
            w[i] = (-2.0 + 2.0 * phi + 2.0 * K + xr - 3.0 * K * xr + 2.0 * xr * xr - xr * xr * xr) / 12.0;
        }
    }
    return;
}

inline void
PETScMatUtilities::ib_6_delta_fcn(const double r, double* const w)
{
    const double rl = r - 2.0;
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
    const double discr = std::max(0.0, beta * beta - 4.0 * alpha * gamma);
    w[0] = (-beta + std::copysign(1.0, (3.0 / 2.0) - K) * std::sqrt(discr)) / (2.0 * alpha);
    w[1] =
        -3.0 * w[0] - (1.0 / 16.0) + (1.0 / 8.0) * (K + r2) + (1.0 / 12.0) * (3.0 * K - 1.0) * rl + (1.0 / 12.0) * r3;
    w[2] = 2.0 * w[0] + (1.0 / 4.0) + (1.0 / 6.0) * (4.0 - 3.0 * K) * rl - (1.0 / 6.0) * r3;
    w[3] = 2.0 * w[0] + (5.0 / 8.0) - (1.0 / 4.0) * (K + r2);
    w[4] = -3.0 * w[0] + (1.0 / 4.0) - (1.0 / 6.0) * (4.0 - 3.0 * K) * rl + (1.0 / 6.0) * r3;
    w[5] = w[0] - (1.0 / 16.0) + (1.0 / 8.0) * (K + r2) - (1.0 / 12.0) * (3.0 * K - 1.0) * rl - (1.0 / 12.0) * r3;
    return;
}

} // namespace IBTK

#endif // #ifndef included_IBTK_PRIVATE_PETScMatUtilities_inl
