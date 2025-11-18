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

// Simple IBAMR example demonstrating basic setup and initialization

#include <ibamr/config.h>

#include <ibtk/IBTKInit.h>
#include <ibtk/SAMRAIDataCache.h>

#include <SAMRAI/tbox/PIO.h>
#include <VariableDatabase.h>

#include <Eigen/Core>

#ifdef IBTK_HAVE_LIBMESH
#include <libmesh/point.h>
#endif

#include <iostream>
#include <memory>

/*!
 * \brief Simple IBAMR test program
 *
 * This example demonstrates:
 * - IBAMR initialization
 * - Basic SAMRAI variable database access
 * - Eigen library integration
 * - Optional libMesh integration
 */
int main(int argc, char **argv)
{
    // Initialize IBAMR and MPI
    auto ibtk_init = std::make_shared<IBTK::IBTKInit>(argc, argv);

    // Print basic information
    SAMRAI::tbox::pout << "=================================================\n";
    SAMRAI::tbox::pout << "Simple IBAMR CMake Example\n";
    SAMRAI::tbox::pout << "=================================================\n";

    // Test SAMRAI variable database access
    SAMRAI::hier::VariableDatabase<NDIM>* var_db =
        SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    SAMRAI::tbox::pout << "Number of registered patch data indices: "
                       << var_db->getNumberOfRegisteredPatchDataIndices() << '\n';

    // Test Eigen library
    SAMRAI::tbox::pout << "Eigen version: "
                       << EIGEN_WORLD_VERSION << '.'
                       << EIGEN_MAJOR_VERSION << '.'
                       << EIGEN_MINOR_VERSION << '\n';

    // Create a simple Eigen matrix
    Eigen::Matrix3d test_matrix;
    test_matrix << 1, 2, 3,
                   4, 5, 6,
                   7, 8, 9;
    SAMRAI::tbox::pout << "Eigen matrix test: trace = "
                       << test_matrix.trace() << '\n';

#ifdef IBTK_HAVE_LIBMESH
    // Test libMesh if available
    const libMesh::Point point(1.0, 2.0, 3.0);
    SAMRAI::tbox::pout << "libMesh is available\n";
    SAMRAI::tbox::pout << "Test point: ("
                       << point(0) << ", "
                       << point(1) << ", "
                       << point(2) << ")\n";
#else
    SAMRAI::tbox::pout << "libMesh is NOT available\n";
#endif

    SAMRAI::tbox::pout << "=================================================\n";
    SAMRAI::tbox::pout << "Simple CMake example completed successfully!\n";
    SAMRAI::tbox::pout << "=================================================\n";

    return 0;
}
