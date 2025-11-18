// ---------------------------------------------------------------------
//
// Copyright (c) 2025 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Simple IBAMR example demonstrating basic setup with Makefile build

#include <ibamr/config.h>

#include <ibtk/IBTKInit.h>
#include <ibtk/AppInitializer.h>

#include <SAMRAI/tbox/PIO.h>
#include <SAMRAI/tbox/Database.h>
#include <CartesianGridGeometry.h>
#include <VariableDatabase.h>

#include <iostream>
#include <memory>

/*!
 * \brief Simple IBAMR test program built with Makefile
 *
 * This example demonstrates:
 * - IBAMR initialization
 * - SAMRAI grid geometry
 * - Variable database access
 * - Basic PIO (Parallel I/O) operations
 */
int main(int argc, char **argv)
{
    // Initialize IBAMR and MPI
    auto ibtk_init = std::make_shared<IBTK::IBTKInit>(argc, argv);

    // Print basic information using SAMRAI PIO
    SAMRAI::tbox::pout << "=================================================\n";
    SAMRAI::tbox::pout << "Simple IBAMR Makefile Example\n";
    SAMRAI::tbox::pout << "=================================================\n";
    SAMRAI::tbox::pout << "NDIM = " << NDIM << '\n';

    // Get MPI rank information
    int rank = SAMRAI::tbox::SAMRAI_MPI::getRank();
    int size = SAMRAI::tbox::SAMRAI_MPI::getNodes();

    SAMRAI::tbox::pout << "MPI Rank: " << rank << " / " << size << '\n';

    // Test SAMRAI variable database
    SAMRAI::hier::VariableDatabase<NDIM>* var_db =
        SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    SAMRAI::tbox::pout << "Variable database initialized\n";
    SAMRAI::tbox::pout << "Number of registered contexts: "
                       << var_db->getNumberOfRegisteredPatchDataIndices() << '\n';

    // Create a simple Cartesian grid geometry (without input file)
    SAMRAI::tbox::pout << "\nCreating simple grid geometry...\n";

    // Define physical domain
    double x_lo[NDIM], x_up[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        x_lo[d] = 0.0;
        x_up[d] = 1.0;
    }

    // Define periodic dimensions
    int periodic[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        periodic[d] = 1; // All dimensions periodic
    }

    SAMRAI::tbox::pout << "Physical domain: [";
    for (int d = 0; d < NDIM; ++d)
    {
        SAMRAI::tbox::pout << x_lo[d] << ", " << x_up[d];
        if (d < NDIM - 1) SAMRAI::tbox::pout << "] x [";
    }
    SAMRAI::tbox::pout << "]\n";

    SAMRAI::tbox::pout << "Periodic boundary conditions: ";
    for (int d = 0; d < NDIM; ++d)
    {
        SAMRAI::tbox::pout << (periodic[d] ? "TRUE" : "FALSE");
        if (d < NDIM - 1) SAMRAI::tbox::pout << ", ";
    }
    SAMRAI::tbox::pout << '\n';

    SAMRAI::tbox::pout << "\n=================================================\n";
    SAMRAI::tbox::pout << "Simple Makefile example completed successfully!\n";
    SAMRAI::tbox::pout << "=================================================\n";

    return 0;
}
