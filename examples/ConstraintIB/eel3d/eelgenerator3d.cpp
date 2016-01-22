// Filename : eelgenerator3d.cpp
// Created by Amneet Bhalla on 05/19/2012.

// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
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

//////////////////////////////////// INCLUDES ////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <numeric>
#include <vector>
#include <utility>

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_integration.h"

// Set Fish Related Parameters.
static const double PII = 3.1415926535897932384626433832795;
static const double LENGTH_FISH = 1.0;
static const double WIDTH_HEAD = 0.04 * LENGTH_FISH;
static const double LENGTH_TILLHEAD = 0.04 * LENGTH_FISH;
static const double MAJOR_AXIS = 0.51 * LENGTH_FISH;
static const double MINOR_AXIS = 0.08 * LENGTH_FISH;

double
xPosition(double s, void* params)
{
    double* input = static_cast<double*>(params);
    double a0 = input[0];
    double a1 = input[1];
    double a2 = input[2];
    double a3 = input[3];
    double tau = input[4];
    double t = input[5];
    double T = input[6];

    double x =
        std::cos((1 / (8 * std::pow(PII * tau, 4))) * (std::tanh(PII * t / T)) *
                 (2 * PII * tau * (a2 - 2 * a0 * std::pow(PII * tau, 2)) * std::cos((2 * PII * t) / T) +
                  2 * PII * tau * (-a2 - 3 * a3 * s + 2 * PII * PII * (a0 + s * (a1 + s * (a2 + a3 * s))) * tau * tau) *
                      cos((2 * PII * (t - s * T * tau)) / T) +
                  (3 * a3 - 2 * a1 * std::pow(PII * tau, 2)) * sin((2 * PII * t) / T) +
                  (-3 * a3 + 2 * PII * PII * (a1 + 2 * a2 * s + 3 * a3 * s * s) * tau * tau) *
                      std::sin((2 * PII * (t - s * T * tau)) / T)));

    return x;

} // xposition

double
yPosition(double s, void* params)
{
    double* input = static_cast<double*>(params);
    double a0 = input[0];
    double a1 = input[1];
    double a2 = input[2];
    double a3 = input[3];
    double tau = input[4];
    double t = input[5];
    double T = input[6];

    double y =
        std::sin((1 / (8 * std::pow(PII * tau, 4))) * (std::tanh(PII * t / T)) *
                 (2 * PII * tau * (a2 - 2 * a0 * std::pow(PII * tau, 2)) * std::cos((2 * PII * t) / T) +
                  2 * PII * tau * (-a2 - 3 * a3 * s + 2 * PII * PII * (a0 + s * (a1 + s * (a2 + a3 * s))) * tau * tau) *
                      cos((2 * PII * (t - s * T * tau)) / T) +
                  (3 * a3 - 2 * a1 * std::pow(PII * tau, 2)) * sin((2 * PII * t) / T) +
                  (-3 * a3 + 2 * PII * PII * (a1 + 2 * a2 * s + 3 * a3 * s * s) * tau * tau) *
                      std::sin((2 * PII * (t - s * T * tau)) / T)));

    return y;

} // yPosition

int
main()
{
    const double Lx = 8;
    const double Ly = 4;
    const double Lz = 1;
    const int Nx = 32 * 4 * 4 * 2;
    const int Ny = 16 * 4 * 4 * 2;
    const int Nz = 8 * 4 * 4 * 2;

    const double dx = Lx / Nx;
    const double dy = Ly / Ny;
    const double dz = Lz / Nz;

    const double tau_tail = 1.52;
    const double a0 = 1.29, a1 = -22.57, a2 = 78.39, a3 = -52.83;
    const double t = 0.0;
    const double time_period = 1.0;

    double interp_coefs[4];
    interp_coefs[0] = a0;
    interp_coefs[1] = a1;
    interp_coefs[2] = a2;
    interp_coefs[3] = a3;

    // No. of points on the backbone and till head.
    const int headNs = int(ceil(LENGTH_TILLHEAD / dx));
    const int tailNs = int(ceil((LENGTH_FISH - LENGTH_TILLHEAD) / dx));
    const int bodyNs = headNs + tailNs + 1;

    std::vector<std::pair<int, int> > immersedBodyData(bodyNs);
    std::vector<std::pair<double, double> > immersedBodyWidthHeight(bodyNs);

    for (int i = 1; i <= headNs + 1; ++i)
    {
        const double s = (i - 1) * dx;
        const double section = sqrt(2 * WIDTH_HEAD * s - s * s);
        const double height = MINOR_AXIS * std::sqrt(1 - pow((s - MAJOR_AXIS) / MAJOR_AXIS, 2));
        const int numPtsInSection = int(ceil(section / dy));
        const int numPtsInHeight = int(ceil(height / dz));
        immersedBodyData[i - 1] = std::make_pair(numPtsInSection, numPtsInHeight);
        immersedBodyWidthHeight[i - 1] = std::make_pair(section, height);
    }

    for (int i = headNs + 2; i <= bodyNs; ++i)
    {
        const double s = (i - 1) * dx;
        const double section = WIDTH_HEAD * (LENGTH_FISH - s) / (LENGTH_FISH - LENGTH_TILLHEAD);
        const double height = MINOR_AXIS * std::sqrt(1 - pow((s - MAJOR_AXIS) / MAJOR_AXIS, 2));
        const int numPtsInSection = int(ceil(section / dy));
        const int numPtsInHeight = int(ceil(height / dz));
        immersedBodyData[i - 1] = std::make_pair(numPtsInSection, numPtsInHeight);
        immersedBodyWidthHeight[i - 1] = std::make_pair(section, height);
    }

    int total_lag_pts = 0;

    double input[7];
    input[0] = interp_coefs[0];
    input[1] = interp_coefs[1];
    input[2] = interp_coefs[2];
    input[3] = interp_coefs[3];
    input[4] = tau_tail;
    input[5] = t;
    input[6] = time_period;

    gsl_function Fx, Fy;
    Fx.function = xPosition;
    Fx.params = input;
    Fy.function = yPosition;
    Fy.params = input;

    double ybase, xbase, errory, errorx;
    size_t nevalsy, nevalsx;

    // Find the deformed shape. Rotate the shape about center of mass.
    std::vector<std::vector<double> > shape_new(3);
    for (int i = 1; i <= bodyNs; ++i)
    {
        const int numPtsInSection = immersedBodyData[i - 1].first;
        const int numPtsInHeight = immersedBodyData[i - 1].second;
        const double width = immersedBodyWidthHeight[i - 1].first;
        const double depth = immersedBodyWidthHeight[i - 1].second;
        const double s = (i - 1) * dx;

        gsl_integration_qng(&Fx, 0, s, 1e-8, 0.0, &xbase, &errorx, &nevalsx);
        gsl_integration_qng(&Fy, 0, s, 1e-8, 0.0, &ybase, &errory, &nevalsy);

        if (numPtsInSection && numPtsInHeight)
        {
            // Fill the middle line first.
            for (int k = -numPtsInHeight; k <= numPtsInHeight; ++k)
            {
                shape_new[0].push_back(xbase);
                shape_new[1].push_back(ybase);
                shape_new[2].push_back(k * dz);
            } // middle line filled.

            total_lag_pts += 2 * numPtsInHeight + 1;

            // Fill the rest of the cross section next.
            for (int j = 1; j <= numPtsInSection; ++j)
            {
                const double y = j * dy;
                for (int k = -numPtsInHeight; k <= numPtsInHeight; ++k)
                {
                    const double z = k * dz;
                    if ((std::pow(y / width, 2) + std::pow(z / depth, 2)) <= 1) // use elliptical cross sections
                    {
                        shape_new[0].push_back(xbase); // right side.
                        shape_new[1].push_back(ybase + y);
                        shape_new[2].push_back(z);

                        shape_new[0].push_back(xbase); // left side.
                        shape_new[1].push_back(ybase - y);
                        shape_new[2].push_back(z);

                        total_lag_pts += 2;
                    }
                }
            } // cross section filled
        }
    }
    std::fstream eelstream;
    eelstream.open("eel3d.vertex", std::fstream::out);
    assert((unsigned)total_lag_pts == shape_new[0].size());
    eelstream << total_lag_pts << "\n";

    for (int k = 1; k <= total_lag_pts; ++k)
        eelstream << shape_new[0][k - 1] + 5.5 << "\t" << shape_new[1][k - 1] << "\t" << shape_new[2][k - 1] << "\n";

    return 0;
}
