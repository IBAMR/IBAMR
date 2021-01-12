// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/FifthOrderStokesWaveGenerator.h"

#include "tbox/Database.h"

#include <cmath>
#include <limits>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FifthOrderStokesWaveGenerator::FifthOrderStokesWaveGenerator(const std::string& object_name, Pointer<Database> input_db)
    : StokesWaveGeneratorStrategy(object_name, input_db)
{
    initStokesCoefficients();
    return;
} // FifthOrderStokesWaveGenerator

double
FifthOrderStokesWaveGenerator::getSurfaceElevation(double x, double time) const
{
    // This computes eta for both shallow and deep water.
    const double& k = d_wave_number;
    const double ka = k * d_amplitude;
    const double omega = sqrt(d_gravity * k) * (d_C[0] + ka * ka * d_C[2] + ka * ka * ka * ka * d_C[4]);
    const double phase = k * x - omega * time;

    return (d_eta[0] * cos(phase) + d_eta[1] * cos(2 * phase) + d_eta[2] * cos(3 * phase) + d_eta[3] * cos(4 * phase) +
            d_eta[4] * cos(5 * phase)) /
           k;

} // getSurfaceElevation

double
FifthOrderStokesWaveGenerator::getVelocity(double x, double z_plus_d, double time, int comp_idx) const
{
    const double& k = d_wave_number;
    const double ka = k * d_amplitude;
    const double omega = sqrt(d_gravity * k) * (d_C[0] + ka * ka * d_C[2] + ka * ka * ka * ka * d_C[4]);
    const double phase = k * x - omega * time;

    if (!d_deep_water_limit)
    {
        double kzd = k * z_plus_d;
        if (comp_idx == 0)
        {
            return d_p[0] * cosh(kzd) * cos(phase) * k + d_p[1] * cosh(2 * kzd) * cos(2 * phase) * k * 2 +
                   d_p[2] * cosh(3 * kzd) * cos(3 * phase) * k * 3 + d_p[3] * cosh(4 * kzd) * cos(4 * phase) * k * 4 +
                   d_p[4] * cosh(5 * kzd) * cos(5 * phase) * k * 5;
        }

        if (comp_idx == 1)
        {
#if (NDIM == 2)
            return d_p[0] * sinh(kzd) * sin(phase) * k + d_p[1] * sinh(2 * kzd) * sin(2 * phase) * k * 2 +
                   d_p[2] * sinh(3 * kzd) * sin(3 * phase) * k * 3 + d_p[3] * sinh(4 * kzd) * sin(4 * phase) * k * 4 +
                   d_p[4] * sinh(5 * kzd) * sin(5 * phase) * k * 5;
#elif (NDIM == 3)
            return 0;
#endif
        }
#if (NDIM == 3)
        if (comp_idx == 2)
        {
            return d_p[0] * sinh(kzd) * sin(phase) * k + d_p[1] * sinh(2 * kzd) * sin(2 * phase) * k * 2 +
                   d_p[2] * sinh(3 * kzd) * sin(3 * phase) * k * 3 + d_p[3] * sinh(4 * kzd) * sin(4 * phase) * k * 4 +
                   d_p[4] * sinh(5 * kzd) * sin(5 * phase) * k * 5;
        }
#endif
    }
    else
    {
        double kz = k * (z_plus_d - d_depth);
        if (comp_idx == 0)
        {
            return d_p[0] * exp(kz) * cos(phase) * k + d_p[1] * exp(2 * kz) * cos(2 * phase) * k * 2 +
                   d_p[2] * exp(3 * kz) * cos(3 * phase) * k * 3;
        }

        if (comp_idx == 1)
        {
#if (NDIM == 2)
            return d_p[0] * exp(kz) * sin(phase) * k + d_p[1] * exp(2 * kz) * sin(2 * phase) * k * 2 +
                   d_p[2] * exp(3 * kz) * sin(3 * phase) * k * 3;
#elif (NDIM == 3)

            return 0;
#endif
        }

#if (NDIM == 3)
        if (comp_idx == 2)
        {
            return d_p[0] * exp(kz) * sin(phase) * k + d_p[1] * exp(2 * kz) * sin(2 * phase) * k * 2 +
                   d_p[2] * exp(3 * kz) * sin(3 * phase) * k * 3;
        }
#endif
    }

    return std::numeric_limits<double>::signaling_NaN();
} // getVelocity

/////////////////////////////// PRIVATE //////////////////////////////////////
void
FifthOrderStokesWaveGenerator::initStokesCoefficients()
{
    const double kd = d_wave_number * d_depth;
    const double ka = d_wave_number * d_amplitude;
    const double ka2 = ka * ka;
    const double ka3 = ka2 * ka;
    const double ka4 = ka3 * ka;
    const double ka5 = ka4 * ka;

    if (d_deep_water_limit)
    {
        // Terms to calculate velocity.
        d_A[1][1] = 1;
        d_A[3][1] = -0.5;
        d_A[4][2] = 0.5;
        d_A[5][1] = -37.0 / 24;
        d_A[5][3] = 1.0 / 12;

        // Terms to calculate eta.
        d_B[2][2] = 0.5;
        d_B[3][1] = -3 / 8.0;
        d_B[4][2] = 1 / 3.0;
        d_B[4][4] = 1 / 3.0;
        d_B[5][3] = 99 / 128.0;
        d_B[5][5] = 125 / 384.0;

        // Wave Dispersion
        d_C[0] = 1;
        d_C[2] = 0.5;
        d_C[4] = 1 / 8.0;

        // Velocity coefficients.
        const double alpha = sqrt(d_gravity / pow(d_wave_number, 3));
        d_p[0] = alpha * (ka * d_A[1][1] + ka3 * d_A[3][1] + ka5 * d_A[5][1]);
        d_p[1] = alpha * (ka4 * d_A[4][2]);
        d_p[2] = alpha * (ka5 * d_A[5][3]);
    }
    else
    {
        double S = 1 / cosh(2 * kd);

        // Terms to calculate velocity.
        d_A[1][1] = 1 / sinh(kd);
        d_A[2][2] = 3 * S * S / (2 * (1 - S) * (1 - S));
        d_A[3][1] = (-4 - 20 * S + 10 * S * S - 13 * S * S * S) / (8 * sinh(kd) * pow((1 - S), 3));
        d_A[3][3] = (-2 * S * S + 11 * S * S * S) / (8 * sinh(kd) * pow((1 - S), 3));
        d_A[4][2] = (12 * S - 14 * S * S - 264 * S * S * S - 45 * S * S * S * S - 13 * S * S * S * S * S) /
                    (24 * pow((1 - S), 5));
        d_A[4][4] = (10 * S * S * S - 174 * S * S * S * S + 291 * pow(S, 5) + 278 * pow(S, 6)) /
                    (48 * (3 + 2 * S) * pow((1 - S), 5));
        d_A[5][1] = (-1184 + 32 * S + 13232 * S * S + 21712 * S * S * S + 20940 * S * S * S * S + 12554 * pow(S, 5) -
                     500 * pow(S, 6) - 3341 * pow(S, 7) - 670 * pow(S, 8)) /
                    (64 * sinh(kd) * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));
        d_A[5][3] = (4 * S + 105 * S * S + 198 * S * S * S - 1376 * S * S * S * S - 1302 * pow(S, 5) - 117 * pow(S, 6) +
                     58 * pow(S, 7)) /
                    (32 * sinh(kd) * (3 + 2 * S) * pow((1 - S), 6));
        d_A[5][5] = (-6 * S * S * S + 272 * S * S * S * S - 1552 * pow(S, 5) + 852 * pow(S, 6) + 2029 * pow(S, 7) +
                     430 * pow(S, 8)) /
                    (64 * sinh(kd) * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));

        // Terms to calculate eta.
        d_B[2][2] = 1 / tanh(kd) * (1 + 2 * S) / (2 * (1 - S));
        d_B[3][1] = -3 * (1 + 3 * S + 3 * S * S + 2 * S * S * S) / (8 * pow((1 - S), 3));
        d_B[4][2] = 1 / tanh(kd) * (6 - 26 * S - 182 * S * S - 204 * S * S * S - 25 * pow(S, 4) + 26 * pow(S, 5)) /
                    (6 * (3 + 2 * S) * pow((1 - S), 4));
        d_B[4][4] = 1 / tanh(kd) * (24 + 92 * S + 122 * S * S + 66 * S * S * S + 67 * pow(S, 4) + 34 * pow(S, 5)) /
                    (24 * (3 + 2 * S) * pow((1 - S), 4));
        d_B[5][3] = 9 *
                    (132 + 17 * S - 2216 * S * S - 5897 * S * S * S - 6292 * pow(S, 4) - 2687 * pow(S, 5) +
                     194 * pow(S, 6) + 467 * pow(S, 7) + 82 * pow(S, 8)) /
                    (128 * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));
        d_B[5][5] = 5 *
                    (300 + 1579 * S + 3176 * S * S + 2949 * S * S * S + 1188 * pow(S, 4) + 675 * pow(S, 5) +
                     1326 * pow(S, 6) + 827 * pow(S, 7) + 130 * pow(S, 8)) /
                    (384 * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));

        // Wave dispersion
        d_C[0] = sqrt(tanh(kd));
        d_C[2] = d_C[0] * (2 + 7 * S * S) / (4 * pow((1 - S), 2));
        d_C[4] = d_C[0] * (4 + 32 * S - 116 * S * S - 400 * S * S * S - 71 * pow(S, 4) + 146 * pow(S, 5)) /
                 (32 * pow((1 - S), 5));

        // Velocity coefficients.
        const double alpha = d_C[0] * sqrt(d_gravity / pow(d_wave_number, 3));
        d_p[0] = alpha * (ka * d_A[1][1] + ka3 * d_A[3][1] + ka5 * d_A[5][1]);
        d_p[1] = alpha * (ka2 * d_A[2][2] + ka4 * d_A[4][2]);
        d_p[2] = alpha * (ka3 * d_A[3][3] + ka5 * d_A[5][3]);
        d_p[3] = alpha * (ka4 * d_A[4][4]);
        d_p[4] = alpha * (ka5 * d_A[5][5]);
    }

    // Component wave eta.
    d_eta[0] = ka + ka3 * d_B[3][1] - ka5 * (d_B[5][3] + d_B[5][5]);
    d_eta[1] = ka2 * d_B[2][2] + ka4 * d_B[4][2];
    d_eta[2] = -ka3 * d_B[3][1] + ka5 * d_B[5][3];
    d_eta[3] = ka4 * d_B[4][4];
    d_eta[4] = ka5 * d_B[5][5];

    return;

} // initStokesCoefficients

} // namespace IBAMR
