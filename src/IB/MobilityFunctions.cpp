// Filename: MobilityFunctions.cpp
// Created on 17 Feb 2016 by Bakytzhan Kallemov and Amneet Bhalla
//
// Copyright (c) 2002-2017, Amneet Bhalla and Boyce Griffith.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "ibamr/MobilityFunctions.h"

#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
#ifndef KRON
#define KRON(i, j) ((i == j) ? 1 : 0) // Kronecker symbol
#endif

const double MOB_FIT_FG_TOL = 1.0e-5; // min distance between blobs to apply empirical fitting
const double ZERO_TOL = 1.0e-10;      // tolerance for zero value
double MOB_FIT_FACTOR;                // constant for normalization
int reUse = 0;                        // flag for reuse data

typedef enum _KERNEL_TYPES
{
    IB3,
    IB4,
    IB6,
    UNKNOWN_TYPE = -1
} KERNEL_TYPES;

// Hydro radii for each kernel IB3,IB4,IB6
const double MOB_FIT_HydroRadius[] = { 0.91, 1.255, 1.4685 };

// current hydrodynamic radius
double HRad;

// coefficients of empirical fit for steady stokes
double F_s[7];
double G_s[4];

// coefficients of empirical fit for finite beta
double F_b[10];
double G_b[6];

// coefficients of empirical fit for f(0)
double Z_b[5];
#if (NDIM == 2)
static double Z_s[3];
#endif

KERNEL_TYPES
GetKernelType(const char* IBKernelName)
{
    if (!std::strcmp(IBKernelName, "IB_3")) return IB3;
    if (!std::strcmp(IBKernelName, "IB_4")) return IB4;
    if (!std::strcmp(IBKernelName, "IB_6")) return IB6;

    std::cerr << "IBEmpiricalMobility: Unknown interpolation kernel type.\n";
    std::exit(EXIT_FAILURE);
} // GetKernelType

// Returns Hydrodynamic radius value
double
getHydroRadius(const char* IBKernelName)
{
    KERNEL_TYPES CurrentKernelType = GetKernelType(IBKernelName);
    return MOB_FIT_HydroRadius[CurrentKernelType];
} // getHydroRadius

// Returns squared norm of the vector.
double
get_sqnorm(const double* a_vec)
{
#if (NDIM == 3)
    return a_vec[0] * a_vec[0] + a_vec[1] * a_vec[1] + a_vec[2] * a_vec[2];
#elif (NDIM == 2)
    return a_vec[0] * a_vec[0] + a_vec[1] * a_vec[1];
#endif

    TBOX_ASSERT(false);
    return -1.0;
} // get_sqnorm

// Returns a value based on linear interpolation of array.
double
InterpolateLinear(const double* Xin, const double* Yin, const int N, double X0)
{
    if (fabs(X0 - Xin[0]) < ZERO_TOL) return Yin[0];
    if (X0 < Xin[0]) return Yin[0] + (Yin[1] - Yin[0]) / (Xin[1] - Xin[0]) * (X0 - Xin[0]);
    if (X0 > Xin[N - 1]) return Yin[N - 1] + (Yin[N - 1] - Yin[N - 2]) / (Xin[N - 1] - Xin[N - 2]) * (X0 - Xin[N - 1]);
    // find nearest neighbour index
    int indx;
    for (indx = 0; indx < N - 1; indx++)
    {
        if (X0 <= Xin[indx + 1]) break;
    }
    return Yin[indx] + (Yin[indx + 1] - Yin[indx]) / (Xin[indx + 1] - Xin[indx]) * (X0 - Xin[indx]);
} // InterpolateLinear

void
InterpolateConstants(KERNEL_TYPES MOB_FIT_current, const double beta)
{
    // data for initialization

    // setting hydrodynamic radius
    HRad = MOB_FIT_HydroRadius[MOB_FIT_current];

#if (NDIM == 3)
    // coefficients for steady stokes 3D
    const double F_stokes_coeff[][3] = { { 1.769, 1.263, 1.097 },   { 1.637, 1.535, 1.60533 },
                                         { 1.362, 0.8134, 0.566 },  { 0.6825, 0.1812, 0.0886 },
                                         { 1.074, 0.6436, 0.4432 }, { 0.6814, 0.181, 0.08848 },
                                         { 0.2609, 0.1774, 0.1567 } };
    const double G_stokes_coeff[][3] = { { 4.314, 11.08, 17.54 }, { 0.1106, 0.1759, 0.2365 } };

    // coefficients for fitting f_beta(0) 3D
    //***********now exist only for IB6 kernel****************
    const int num_cases = 8;
    // Number cases is 8, affects array sizes
    double M_betas[] = { 0., 0.1, 0.25, 0.5, 1., 10., 100., 1000.0 };

    // coefficient to empirical fitting selfmobility f_beta(0)
    double F_beta_zero[] = { 1.0 / 0.0230502, -0.131445787110707, 0.287509131179951, 61.267737026060075 };

    double F_beta_coeff[][8] = { { 2.68030156113612,
                                   2.29464251115180,
                                   2.09297672079185,
                                   1.91583152825747,
                                   1.86066266441860,
                                   1.08621320400389,
                                   1.81206666469132,
                                   0.310021323508742 },
                                 { 1.39009017417989,
                                   1.40943006156447,
                                   1.35006023047952,
                                   1.94612277585268,
                                   3.86995647537320,
                                   11.5138612720591,
                                   1.11416207665514,
                                   -0.0372181749138520 },
                                 { -0.884827389490450,
                                   -0.756475504578064,
                                   -0.608404711842046,
                                   -0.564784303490654,
                                   -0.651175884725500,
                                   -0.427005600208834,
                                   0.00756892008838000,
                                   0.104913106939363 },
                                 { 1.0,
                                   1.06262915172038,
                                   1.49440529353864,
                                   2.03210062812245,
                                   2.18657011795234,
                                   3.37689526048115,
                                   3.67794210129942,
                                   3.67794210129942 },
                                 { 2.92685539141722,
                                   2.67245966301464,
                                   2.42487239022736,
                                   2.70289015143191,
                                   3.34961045584841,
                                   3.07595994173387,
                                   0.174756918359867,
                                   0.129697932247353 },
                                 { 29.0233576593274,
                                   26.6542400026965,
                                   24.6035062102207,
                                   24.4821402826182,
                                   28.4998841491606,
                                   32.3675245644131,
                                   17.3333226235307,
                                   2.20851998706445 },
                                 { 1.67431017814765,
                                   1.55636776784166,
                                   1.39866915583103,
                                   1.25760476142372,
                                   1.26454914224257,
                                   1.54512892334494,
                                   1.04695565772442,
                                   0.470157376261480 },
                                 { 1.18077832670187,
                                   0.866186289480926,
                                   0.456603306920829,
                                   0.157941170393203,
                                   0.133374232562159,
                                   2.93244835654164,
                                   1.07958559850185,
                                   0.221651722533943 },
                                 { 0.0314600987678820,
                                   0.0296450899554450,
                                   0.0322750684446850,
                                   0.0296570052703910,
                                   0.0140708796915960,
                                   0.0426413220771180,
                                   0.0232409154706480,
                                   0.00157997335725900 },
                                 { 0.0320036255275870,
                                   0.0210726756307760,
                                   0.0110662759921420,
                                   0.00380530089032800,
                                   0.00265521991315500,
                                   0.00541862366430600,
                                   5.09744623703744e-05,
                                   0 } };

    double G_beta_coeff[][8] = { { 2.56932494798493,
                                   0.884354880187973,
                                   0.113504063191842,
                                   -0.221849887772746,
                                   0.542698003099015,
                                   0.276258520317914,
                                   0.330745197269599,
                                   -7.88722383366528e-05 },
                                 { -0.171344575435234,
                                   0.431984145806452,
                                   0.482467326122353,
                                   0.461388580082004,
                                   -0.0463322271962110,
                                   -0.131688108575531,
                                   -0.114177160773643,
                                   0.0119064417188160 },
                                 { -0.0846371631980100,
                                   -0.0746824104260310,
                                   -0.0822639676502030,
                                   -0.0689510022688410,
                                   0.108601073154871,
                                   0.105609456589694,
                                   0.0941995195204770,
                                   0.0575579089490820 },
                                 { 1.0,
                                   0.356391856987450,
                                   0.495958169969290,
                                   0.651495025037754,
                                   0.608684335056712,
                                   9.63355010424247,
                                   2.60454636502040,
                                   3.17882803024967 },
                                 { 2.21329860848447,
                                   0.894411380099355,
                                   0.698529765627973,
                                   0.552224660456298,
                                   0.723993778966245,
                                   0.423143094120794,
                                   0.359827255706337,
                                   0.267507608759639 },
                                 { 0.,
                                   0.00283198188374900,
                                   0.00380833037508200,
                                   0.00266380633550800,
                                   -0.00766798833782700,
                                   -0.00440457501794600,
                                   -0.00100982221976800,
                                   1.90168352243139e-05 } };

#elif (NDIM == 2)
    // coefficients for fitting f_beta(0) 2D
    const int num_cases = 7;
    const double F_beta_zero[3][5] = { { 0.1245, 0.05237, 1.873, 0.5579, -0.0159 },
                                       { 0.07016, 0.01957, 1.229, 0.2613, -0.00793 },
                                       { 0.05307, 0.01246, 1.015, 0.1882, -0.005932 } };
    const double F_stokes_zero[3][3] = { { 0., 0., 0. }, { 0., 0., 0. }, { 0.1383, 0.07947, 0.7942 } };

    double M_betas[7] = { 0., 0.1, 0.25, 0.5, 1.0, 5.0, 10 };

    // coefficients for steady stokes in 2d
    const double F_stokes_coeff[3][5] = {
        { 0., 0., 0., 0., 0. },
        { 0., 0., 0., 0., 0. },
        { 0.011515464372641, 0.009893387219365, 0.049266457632092, -1.919065487969394, 0.340774530921814 }
    };
    const double G_stokes_coeff[3][4] = { { 0., 0., 0., 0. },
                                          { 0., 0., 0., 0. },
                                          { 0.1658, 0.007702, -0.1309, 0.1752 } };

    // coefficients for finite beta in 2d
    double F_beta_coeff[3][8][7] = { // IB3
                                     { { 1.678, 2.351, 6.257, 8.813, 7.434, 8.154, 17.1 },
                                       { 0.2768, 0.04019, -0.1125, -0.703, -0.2272, -0.3078, 1.097 },
                                       { -0.03659, -0.01597, -0.0005797, -0.01151, -0.034, -0.01632, 0.0001102 },
                                       { 0.6288, 0.4527, 0.3291, 0.2991, 0.2966, 0.1154, 0.04601 },
                                       { -0.3928, -0.2737, -0.2201, -0.1791, -0.1459, -0.02947, -0.01153 },
                                       { 0.6844, 0.2751, -0.01444, 0.1128, 0.504, 0.7547, 0.3005 },
                                       { -0.9516, -0.3658, 0.06593, 0.06458, -0.04301, -0.07149, -0.01518 },
                                       { 0.5129, 0.2532, 0.07412, 0.04351, 0.04404, 0.009016, 0.00175 } },
                                     // IB4
                                     { { 1.966, 2.28, 2.869, 4.689, 8.277, 18.72, 23.45 },
                                       { 3.215, 3.261, 2.174, 0.5451, 0.1269, -0.9027, -0.9177 },
                                       { -0.0478, -0.03822, -0.02346, -0.01224, -0.01957, -0.02191, -0.007607 },
                                       { 0.4169, 0.348, 0.2652, 0.1992, 0.1799, 0.1015, 0.05124 },
                                       { -0.2006, -0.16, -0.1156, -0.08562, -0.07359, -0.02876, -0.01186 },
                                       { 0.801, 0.6064, 0.343, 0.1487, 0.1572, 0.3007, 0.2231 },
                                       { -0.6851, -0.473, -0.2444, -0.06896, -0.01421, -0.01243, -0.008669 },
                                       { 0.2183, 0.1514, 0.08827, 0.04013, 0.01972, 0.003834, 0.001323 } },
                                     // IB6
                                     { { 2.335, 2.535, 3.042, 4.084, 7.063, 8.999, 10.14 },
                                       { 5.29, 5.373, 4.903, 3.417, 1.568, 0, 0 },
                                       { -0.04505, -0.03687, -0.02636, -0.01617, -0.01339, -0.02194, -0.01141 },
                                       { 0.3293, 0.2808, 0.2269, 0.1745, 0.1393, 0.09043, 0.05078 },
                                       { -0.1387, -0.1141, -0.08859, -0.06522, -0.05038, -0.01973, -0.007966 },
                                       { 0.7132, 0.5568, 0.3829, 0.2235, 0.1364, 0.3424, 0.2881 },
                                       { -0.4959, -0.362, -0.2269, -0.1084, -0.03004, -0.03273, -0.02078 },
                                       { 0.1313, 0.09624, 0.0633, 0.0353, 0.01592, 0.005121, 0.001885 } }
    };

    double G_beta_coeff[3][4][7] = { // IB3
                                     { { 2.64, 0.4042, 0.5145, 0.6641, 0.6132, 0.2806, 0.1883 },
                                       { 64.44, 12.91, 15.72, 26.04, 46.2, 101.4, 141.7 },
                                       { -127.2, -4.817, -3.379, -10.26, -27.95, -40.57, -32.57 },
                                       { 120.3, 0.8398, 1.278, 6.632, 15.22, 22.21, 24.42 } },
                                     // IB4
                                     { { 2.255, 1.87, 1.357, 0.8759, 0.6273, 0.2887, 0.1937 },
                                       { 254, 181.9, 92.03, 65.4, 89.45, 198.1, 260.2 },
                                       { -429.8, -249.2, -69.93, -20.31, -36.4, -73.94, -66.77 },
                                       { 299.9, 172.6, 57.18, 17.99, 18.35, 27.33, 28.69 } },
                                     // IB6
                                     { { 1.918, 1.648, 1.321, 0.9552, 0.6402, 0.291, 0.1964 },
                                       { 423.3, 302, 183.1, 119.8, 122.1, 271.9, 357.2 },
                                       { -620.4, -378, -165.5, -56.36, -40.5, -93.04, -90.56 },
                                       { 353.1, 218.5, 104.6, 40.79, 22.38, 30.02, 31.46 } }
    };

#endif

    int cnt;
#if (NDIM == 3)
    //********3D case
    // setting coeficients for steady stokes fitting
    for (cnt = 0; cnt < 7; cnt++) F_s[cnt] = F_stokes_coeff[cnt][MOB_FIT_current];
    for (cnt = 0; cnt < 2; cnt++) G_s[cnt] = G_stokes_coeff[cnt][MOB_FIT_current];

    // setting coeficients for time dependent fitting (curently only for IB6 kernel)
    for (cnt = 0; cnt < 4; cnt++) Z_b[cnt] = F_beta_zero[cnt];
    for (cnt = 0; cnt < 10; cnt++) F_b[cnt] = InterpolateLinear(M_betas, F_beta_coeff[cnt], num_cases, beta);
    for (cnt = 0; cnt < 6; cnt++) G_b[cnt] = InterpolateLinear(M_betas, G_beta_coeff[cnt], num_cases, beta);

#elif (NDIM == 2)
    for (cnt = 0; cnt < 5; cnt++) Z_b[cnt] = F_beta_zero[MOB_FIT_current][cnt];
    for (cnt = 0; cnt < 3; cnt++) Z_s[cnt] = F_stokes_zero[MOB_FIT_current][cnt];
    for (cnt = 0; cnt < 5; cnt++) F_s[cnt] = F_stokes_coeff[MOB_FIT_current][cnt];
    for (cnt = 0; cnt < 4; cnt++) G_s[cnt] = G_stokes_coeff[MOB_FIT_current][cnt];
    for (cnt = 0; cnt < 8; cnt++)
        F_b[cnt] = InterpolateLinear(M_betas, F_beta_coeff[MOB_FIT_current][cnt], num_cases, beta);
    for (cnt = 0; cnt < 4; cnt++)
        G_b[cnt] = InterpolateLinear(M_betas, G_beta_coeff[MOB_FIT_current][cnt], num_cases, beta);
#endif
    return;
} // InterpolateConstants

void
InitializeAllConstants(const char* IBKernelName, const double MU, const double rho, const double Dt, const double DX)
{
    KERNEL_TYPES CurrentKernelType = GetKernelType(IBKernelName);
    double beta;
    // finding beta
    if (MU <= ZERO_TOL)
        beta = 0.0; // invisid case
    else
        beta = MU * Dt / (rho * DX * DX);

#if (NDIM == 3)
    //******3D case
    if ((rho < ZERO_TOL) || (beta >= 1000.1))
        MOB_FIT_FACTOR = 1. / MU / DX; // 3D steady stokes
    else
        MOB_FIT_FACTOR = Dt / (rho * DX * DX * DX);
#elif (NDIM == 2)
    //*******2D case
    if ((rho < ZERO_TOL) || (beta >= 100.1)) // 2D steady stokes
        MOB_FIT_FACTOR = 1. / MU;
    else
        MOB_FIT_FACTOR = Dt / (rho * DX * DX);
#endif

    InterpolateConstants(CurrentKernelType, beta);

} // InitializeAllConstants

double
_F_R_INF(const double rr, const double Dx, const double L_domain)
{
    const double r = rr / Dx;
#if (NDIM == 3)
    (void)L_domain;
    const double factor = 1.0 / (8.0 * M_PI);
    if (r < 0.8)
        return factor / (3.0 / 4.0 * HRad + F_s[6] * r * r);
    else
        return factor * (std::exp(-F_s[0] * r) * F_s[1] / HRad +
                         (F_s[2] * r + F_s[3] * r * r * r) / (1.0 + F_s[4] * r * r + F_s[5] * std::pow(r, 4)));
#elif (NDIM == 2)
    if (L_domain < ZERO_TOL)
    {
        std::cerr << "IBEMpiricalMobility:_F_R_INF()  L_domain must be non specified in 2D!. Abort.\n";
        std::exit(EXIT_FAILURE);
    }
    double f_0 = (Z_s[0] + Z_s[1] * std::log(L_domain / Z_s[2]));
    if (r < 0.1)
        return f_0;
    else
        return (f_0 + (F_s[0] * r * r + F_s[1] * r * r * r + F_s[2] * r * r * r * std::log(r)) /
                          (1.0 + F_s[3] * r + F_s[4] * r * r - F_s[2] * 4 * M_PI * r * r * r));
#endif
} // _F_R_INF

double
_G_R_INF(const double rr, const double Dx)
{
    const double r = rr / Dx;
#if (NDIM == 3)
    const double factor = 1.0 / (8.0 * M_PI);
    if (r < MOB_FIT_FG_TOL)
        return 0.0;
    else
        return (factor * r * r) / (G_s[0] + G_s[1] * r * r + r * r * r);
#elif (NDIM == 2)
    return (G_s[0] * r * r + G_s[1] * r * r * r) / (1.0 + G_s[2] * r + G_s[3] * r * r + G_s[1] * r * r * r) / 4 / M_PI;
#endif
} // _G_R_INF

double
_F_R_BETA(const double rr, const double Dx, const double beta, const double L_domain)
{
    const double r = rr / Dx;

#if (NDIM == 3)
    (void)L_domain;
    double f_0 = (1.0 + Z_b[1] * std::sqrt(beta) + Z_b[2] * beta) /
                 (Z_b[0] + Z_b[3] * beta + Z_b[2] * 6.0 * M_PI * HRad * beta * beta);

    // invisid case
    if (beta < ZERO_TOL)
        return f_0 / 4.0 / M_PI *
               ((4.0 * M_PI - F_b[4] * r * r) /
                    (1.0 + F_b[0] * r + F_b[1] * r * r + F_b[2] * r * r * r + F_b[4] * f_0 * std::pow(r, 5)) +
                (F_b[5] * r * std::exp(-F_b[6] * r) + F_b[7] * r) /
                    (1.0 + F_b[8] * r * r * r + F_b[9] * std::pow(r, 5)));
    else if (beta < 1000.1)
        return f_0 / 4.0 / M_PI *
               ((4.0 * M_PI +
                 F_b[4] * (-r * r + std::pow(r, 4) * std::exp(-F_b[3] * r / std::sqrt(beta)) / (2.0 * beta))) /
                    (1.0 + F_b[0] * r + F_b[1] * r * r + F_b[2] * r * r * r + F_b[4] * f_0 * std::pow(r, 5)) +
                (F_b[5] * r * std::exp(-F_b[6] * r) + F_b[7] * r) /
                    (1.0 + F_b[8] * r * r * r + F_b[9] * std::pow(r, 5)));
    else
        return _F_R_INF(rr, Dx, 0.0);

#elif (NDIM == 2)
    const double factor = 1.0 / M_PI;
    if (beta < 100.1)
    {
        if (r < 0.8) // if distance is less self mobility is used
            return (Z_b[0] + Z_b[1] * beta * beta) /
                   (Z_b[4] * std::pow(beta, 4) + Z_b[3] * beta * beta * beta + Z_b[2] * beta + 1.0) /
                   (1.0 + r * r / (1.7 + 1.5 * std::log(1.0 + beta)));
        else
            return factor * (-0.5 * std::exp(-F_b[0] / r) / (F_b[1] + r * r) +
                             (F_b[2] + F_b[3] * r + F_b[4] * r * r) /
                                 (1.0 + F_b[5] * r * r + F_b[6] * r * r * r + F_b[7] * std::pow(r, 4)) / r);
    }
    else
    {
        return _F_R_INF(rr, Dx, L_domain);
    }

#endif
} // _F_R_BETA

double
_G_R_BETA(const double rr, const double Dx, const double beta)
{
    const double r = rr / Dx;

#if (NDIM == 3)
    double f_0 = (1.0 + Z_b[1] * std::sqrt(beta) + Z_b[2] * beta) /
                 (Z_b[0] + Z_b[3] * beta + Z_b[2] * 6.0 * M_PI * HRad * beta * beta);
    if (beta < ZERO_TOL)
        return f_0 * 3.0 / (4.0 * M_PI) * G_b[4] * r * r /
               (1.0 + G_b[0] * r + G_b[1] * r * r + G_b[2] * r * r * r + G_b[4] * f_0 * std::pow(r, 5));
    else if (beta < 1001.0)
        return f_0 * 3.0 / (4.0 * M_PI) * G_b[4] *
               (r * r + std::pow(r, 4) * std::exp(-G_b[3] * r / std::sqrt(beta)) / (6.0 * beta)) /
               (1.0 + G_b[0] * r + G_b[1] * r * r + G_b[2] * r * r * r + G_b[5] * std::pow(r, 4) +
                G_b[4] * f_0 * std::pow(r, 5));
    else
        return _G_R_INF(rr, Dx);
#elif (NDIM == 2)
    if (beta < 100.1)
    {
        const double factor = 1.0 / M_PI;
        if (r < MOB_FIT_FG_TOL)
            return 0.0;
        else
        {
            return factor * r / (std::exp(-G_b[0] * r) * (G_b[1] + G_b[2] * r + G_b[3] * r * r) + r * r * r);
        }
    }
    else
    {
        return _G_R_INF(rr, Dx);
    }
#endif
} // _G_R_BETA

// Computes Empirical Mobility components f(r) and g(r)
void
getEmpiricalMobilityComponents(const char* IBKernelName,
                               const double MU,
                               const double rho,
                               const double Dt,
                               const double r,
                               const double DX,
                               const int resetAllConstants,
                               const double L_domain,
                               double* F_MobilityValue,
                               double* G_Mobilityvalue)
{
    // Reuse same static constants for efficiency
    if (resetAllConstants) reUse = 0;
    if (!reUse)
    {
        InitializeAllConstants(IBKernelName, MU, rho, Dt, DX);
        reUse = 1;
    }
    double beta;
    // finding beta
    if (MU <= ZERO_TOL)
        beta = 0.0; // invisid case
    else
        beta = MU * Dt / (rho * DX * DX);

    if (rho < ZERO_TOL)
    {
        *F_MobilityValue = MOB_FIT_FACTOR * _F_R_INF(r, DX, L_domain); // steady stokes term for f(r)
        *G_Mobilityvalue = MOB_FIT_FACTOR * _G_R_INF(r, DX);           // steady stokes term for g(r)
    }
    else
    {
        *F_MobilityValue = MOB_FIT_FACTOR * _F_R_BETA(r, DX, beta, L_domain); // time-dependent f(r)
        *G_Mobilityvalue = MOB_FIT_FACTOR * _G_R_BETA(r, DX, beta);           // time-dependent g(r)
    }
    return;
} // getEmpiricalMobilityComponents
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
MobilityFunctions::constructEmpiricalMobilityMatrix(const char* IBKernelName,
                                                    const double MU,
                                                    const double rho,
                                                    const double Dt,
                                                    const double DX,
                                                    const double* X,
                                                    const int N,
                                                    const int resetAllConstants,
                                                    const double /*PERIODIC_CORRECTION*/,
                                                    const double L_domain,
                                                    double* MM)
{
    int row, col;
    for (row = 0; row < N; row++)
        for (col = 0; col <= row; col++)
        {
            double r_vec[NDIM];
            const int size = N * NDIM;
            int cdir;
            for (cdir = 0; cdir < NDIM; cdir++)
            {
                r_vec[cdir] = X[row * NDIM + cdir] - X[col * NDIM + cdir]; // r(i) - r(j)
            }

            const double rsq = get_sqnorm(r_vec);
            const double r = std::sqrt(rsq);
            double F_R, G_R;

            getEmpiricalMobilityComponents(IBKernelName, MU, rho, Dt, r, DX, resetAllConstants, L_domain, &F_R, &G_R);

            int idir, jdir;
            for (idir = 0; idir < NDIM; idir++)
                for (jdir = 0; jdir <= idir; jdir++)
                {
                    const int index = (col * NDIM + jdir) * size + row * NDIM + idir; // column-major for LAPACK
                    MM[index] = F_R * KRON(idir, jdir);
                    if (row != col)
                    {
                        MM[index] += G_R * r_vec[idir] * r_vec[jdir] / rsq;
                        MM[(row * NDIM + idir) * size + col * NDIM + jdir] = MM[index];
                        if (idir != jdir) MM[(row * NDIM + jdir) * size + col * NDIM + idir] = MM[index];
                    }
                    if (idir != jdir) MM[(col * NDIM + idir) * size + row * NDIM + jdir] = MM[index];
                }
        }
    return;
} // constructEmpiricalMobilityMatrix

void
MobilityFunctions::constructRPYMobilityMatrix(const char* IBKernelName,
                                              const double MU,
                                              const double DX,
                                              const double* X,
                                              const int N,
                                              const double PERIODIC_CORRECTION,
                                              double* MM)
{
    HRad = getHydroRadius(IBKernelName) * DX;
    const double mu_tt = 1. / (6.0 * M_PI * MU * HRad);

    double r_vec[NDIM];
    int size = N * NDIM;
    int row, col;
    for (row = 0; row < N; row++)
        for (col = 0; col <= row; col++)
        {
            if (row == col)
            {
                int idir, jdir;
                for (idir = 0; idir < NDIM; idir++)
                    for (jdir = 0; jdir < NDIM; jdir++)
                    {
                        const int index = (col * NDIM + jdir) * size + row * NDIM + idir; // column-major for LAPACK
                        MM[index] = (mu_tt - PERIODIC_CORRECTION) * KRON(idir, jdir);
                    }
            }
            else
            {
                int cdir;
                for (cdir = 0; cdir < NDIM; cdir++)
                {
                    r_vec[cdir] = X[row * NDIM + cdir] - X[col * NDIM + cdir]; // r(i) - r(j)
                }

                const double rsq = get_sqnorm(r_vec);
                const double r = std::sqrt(rsq);
                int idir, jdir;
                for (idir = 0; idir < NDIM; idir++)
                    for (jdir = 0; jdir <= idir; jdir++)
                    {
                        const int index = (col * NDIM + jdir) * size + row * NDIM + idir; // column-major for LAPACK
                        if (r <= 2.0 * HRad)
                        {
                            MM[index] = (mu_tt * (1 - 9.0 / 32.0 * r / HRad) - PERIODIC_CORRECTION) * KRON(idir, jdir) +
                                        mu_tt * r_vec[idir] * r_vec[jdir] / rsq * 3.0 * r / 32. / HRad;
                            MM[(row * NDIM + idir) * size + col * NDIM + jdir] = MM[index];
                            if (idir != jdir)
                            {
                                MM[(col * NDIM + idir) * size + row * NDIM + jdir] = MM[index];
                                MM[(row * NDIM + jdir) * size + col * NDIM + idir] = MM[index];
                            }
                        }
                        else
                        {
                            double cube = HRad * HRad * HRad / r / r / r;
                            MM[index] =
                                (mu_tt * (3.0 / 4.0 * HRad / r + 1.0 / 2.0 * cube) - PERIODIC_CORRECTION) *
                                    KRON(idir, jdir) +
                                mu_tt * r_vec[idir] * r_vec[jdir] / rsq * (3.0 / 4.0 * HRad / r - 3.0 / 2.0 * cube);
                            MM[(row * NDIM + idir) * size + col * NDIM + jdir] = MM[index];
                            if (idir != jdir)
                            {
                                MM[(col * NDIM + idir) * size + row * NDIM + jdir] = MM[index];
                                MM[(row * NDIM + jdir) * size + col * NDIM + idir] = MM[index];
                            }
                        }
                    } // jdir
            }
        } // column loop
    return;
} // constructRPYMobilityMatrix

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
