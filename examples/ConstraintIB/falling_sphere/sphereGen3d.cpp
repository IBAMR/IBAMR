// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Written by namu patel on mahavira@mech.northwestern.edu

#include <math.h>

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// Set sphere related parameters.
static const double DIAMETER = 0.625;
static const int NDIM = 3;

int
main()
{
    // Eulerian domain
    const double Lx = 2.0;
    const double Ly = 8.0;
    const double Lz = 2.0;
    const int Nx = 8 * 4 * 2;
    const int Ny = 32 * 4 * 2;
    const int Nz = 8 * 4 * 2;
    const double dx = Lx / static_cast<double>(Nx);
    const double dy = Ly / static_cast<double>(Ny);
    const double dz = Lz / static_cast<double>(Nz);

    // Lagrangian domain
    const double DX = dx;
    const double DY = dy;
    const double DZ = dz;
    const int NX = ceil(DIAMETER / dx);
    const int NY = ceil(DIAMETER / dy);
    const int NZ = ceil(DIAMETER / dz);
    const double X_cm = 0.0;
    const double Y_cm = 0.0;
    const double Z_cm = 0.0;

    int total_lag_pts = 0;
    std::vector<std::vector<double> > sphere_coords(NDIM);
    for (int k = 0; k < NZ; ++k)
    {
        const double Z = Z_cm - 0.5 * DIAMETER + k * DZ;
        for (int i = 0; i < NX; ++i)
        {
            const double X = X_cm - 0.5 * DIAMETER + i * DX;
            for (int j = 0; j < NY; ++j)
            {
                const double Y = Y_cm - 0.5 * DIAMETER + j * DY;

                if (((X - X_cm) * (X - X_cm) + (Y - Y_cm) * (Y - Y_cm) + (Z - Z_cm) * (Z - Z_cm)) <
                    (0.25 * DIAMETER * DIAMETER))
                {
                    sphere_coords[0].push_back(X);
                    sphere_coords[1].push_back(Y);
                    sphere_coords[2].push_back(Z);
                    ++total_lag_pts;
                }
            }
        }
    }

    std::fstream sphere_coord_stream;
    sphere_coord_stream.open("sphere3d.vertex", std::fstream::out);
    assert(static_cast<std::size_t>(total_lag_pts) == sphere_coords[0].size());
    sphere_coord_stream << total_lag_pts << "\n";
    for (int k = 0; k < total_lag_pts; ++k)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            sphere_coord_stream << std::fixed << std::setprecision(7) << std::setfill('0') << sphere_coords[d][k]
                                << "\t";
        }
        sphere_coord_stream << "\n";
    }

    return 0;
}
