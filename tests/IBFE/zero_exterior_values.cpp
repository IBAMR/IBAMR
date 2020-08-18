// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files

#include <SAMRAI_config.h>

// Headers for basic SAMRAI objects
#include <CartesianPatchGeometry.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/FEDataManager.h>

#include <fstream>
#include <vector>

int
main()
{
    std::ofstream out("output");

    SAMRAI::tbox::Array<SAMRAI::tbox::Array<bool> > bdry_data(NDIM);
    for (int i = 0; i < NDIM; ++i)
    {
        SAMRAI::tbox::Array<bool>& entry = bdry_data[i];
        entry.resizeArray(2);
        entry[0] = false;
        entry[1] = false;
    }
    const double dx[NDIM] = { 1.0 / 8.0 };

#if NDIM == 2
    {
        const double x_lo[NDIM] = { 0.0, 1.0 };
        const double x_up[NDIM] = { 1.0, 3.0 };
        SAMRAI::geom::CartesianPatchGeometry<NDIM> patch_geo(
            SAMRAI::hier::IntVector<NDIM>(1), bdry_data, bdry_data, dx, x_lo, x_up);

        std::vector<double> coordinates{
            0.0, 0.0, // outside patch
            0.0, 2.0, // left boundary (x = 0)
            1.0, 1.5, // right boundary (x = 1)
            0.5, 1.0, // bottom boundary (y = 1)
            0.9, 3.0, // top boundary (y = 3)
            0.9, 2.9  // interior
        };

        std::vector<double> values{ 1.0, 1.1, 2.0, 2.1, 3.0, 3.1, 4.0, 4.1, 5.0, 5.1, 6.0, 6.1 };
        IBTK::FEDataManager::zeroExteriorValues(patch_geo, coordinates, values, NDIM);

        for (unsigned int qp_n = 0; qp_n < 6; ++qp_n)
        {
            out << coordinates[qp_n * NDIM + 0] << ", " << coordinates[qp_n * NDIM + 1] << ": "
                << values[qp_n * NDIM + 0] << ", " << values[qp_n * NDIM + 1] << '\n';
        }
    }
#elif NDIM == 3
    {
        const double x_lo[NDIM] = { 0.0, 1.0, 2.0 };
        const double x_up[NDIM] = { 1.0, 3.0, 6.0 };

        SAMRAI::geom::CartesianPatchGeometry<NDIM> patch_geo(
            SAMRAI::hier::IntVector<NDIM>(1), bdry_data, bdry_data, dx, x_lo, x_up);

        std::vector<double> coordinates{
            0.0, 0.0, 0.0, // outside patch
            0.0, 2.0, 5.0, // left boundary (x = 0)
            1.0, 1.5, 4.0, // right boundary (x = 1)
            0.5, 1.0, 2.5, // front boundary (y = 1)
            0.5, 3.0, 2.5, // back boundary (y = 3)
            0.5, 1.5, 2.0, // bottom boundary (z = 0)
            0.9, 2.0, 6.0, // top boundary (z = 6)
            0.9, 2.9, 5.9  // interior
        };

        std::vector<double> values{ 1.0, 1.1, 1.2, 2.0, 2.1, 2.2, 3.0, 3.1, 3.2, 4.0, 4.1, 4.2,
                                    5.0, 5.1, 5.2, 6.0, 6.1, 6.2, 7.0, 7.1, 7.2, 8.0, 8.1, 8.2 };
        IBTK::FEDataManager::zeroExteriorValues(patch_geo, coordinates, values, NDIM);

        for (unsigned int qp_n = 0; qp_n < 8; ++qp_n)
        {
            out << coordinates[qp_n * NDIM + 0] << ", " << coordinates[qp_n * NDIM + 1] << ", "
                << coordinates[qp_n * NDIM + 2] << ": " << values[qp_n * NDIM + 0] << ", " << values[qp_n * NDIM + 1]
                << ", " << values[qp_n * NDIM + 2] << '\n';
        }
    }
#endif
}
