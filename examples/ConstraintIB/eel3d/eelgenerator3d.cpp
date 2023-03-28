// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

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
    const int headNs = static_cast<int>(ceil(LENGTH_TILLHEAD / dx));
    const int tailNs = static_cast<int>(ceil((LENGTH_FISH - LENGTH_TILLHEAD) / dx));
    const int bodyNs = headNs + tailNs + 1;

    std::vector<std::pair<int, int> > immersedBodyData(bodyNs);
    std::vector<std::pair<double, double> > immersedBodyWidthHeight(bodyNs);

    for (int i = 1; i <= headNs + 1; ++i)
    {
        const double s = (i - 1) * dx;
        const double section = sqrt(2 * WIDTH_HEAD * s - s * s);
        const double height = MINOR_AXIS * std::sqrt(1 - pow((s - MAJOR_AXIS) / MAJOR_AXIS, 2));
        const int numPtsInSection = static_cast<int>(ceil(section / dy));
        const int numPtsInHeight = static_cast<int>(ceil(height / dz));
        immersedBodyData[i - 1] = std::make_pair(numPtsInSection, numPtsInHeight);
        immersedBodyWidthHeight[i - 1] = std::make_pair(section, height);
    }

    for (int i = headNs + 2; i <= bodyNs; ++i)
    {
        const double s = (i - 1) * dx;
        const double section = WIDTH_HEAD * (LENGTH_FISH - s) / (LENGTH_FISH - LENGTH_TILLHEAD);
        const double height = MINOR_AXIS * std::sqrt(1 - pow((s - MAJOR_AXIS) / MAJOR_AXIS, 2));
        const int numPtsInSection = static_cast<int>(ceil(section / dy));
        const int numPtsInHeight = static_cast<int>(ceil(height / dz));
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

    // Find the deformed shape. Rotate the shape about center of mass.
    std::vector<std::vector<double> > shape_new(3);
    for (int i = 1; i <= bodyNs; ++i)
    {
        const int numPtsInSection = immersedBodyData[i - 1].first;
        const int numPtsInHeight = immersedBodyData[i - 1].second;
        const double width = immersedBodyWidthHeight[i - 1].first;
        const double depth = immersedBodyWidthHeight[i - 1].second;
        const double s = (i - 1) * dx;

        auto f_x = [&](const double x) { return xPosition(x, input); };

        auto f_y = [&](const double y) { return yPosition(y, input); };

        namespace bmq = boost::math::quadrature;
        const double xbase =
            s == 0.0 ? 0.0 : bmq::gauss_kronrod<double, 15>::integrate(f_x, 0.0, s, 15, 1e-12, nullptr);
        const double ybase =
            s == 0.0 ? 0.0 : bmq::gauss_kronrod<double, 15>::integrate(f_y, 0.0, s, 15, 1e-12, nullptr);

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
    assert(static_cast<std::size_t>(total_lag_pts) == shape_new[0].size());
    eelstream << total_lag_pts << "\n";

    for (int k = 1; k <= total_lag_pts; ++k)
        eelstream << shape_new[0][k - 1] + 5.5 << "\t" << shape_new[1][k - 1] << "\t" << shape_new[2][k - 1] << "\n";

    return 0;
}
