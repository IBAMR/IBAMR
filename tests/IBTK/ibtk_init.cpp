// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBTK
    IBTKInit ibtkInit(argc, argv, PETSC_COMM_WORLD);

#ifdef IBTK_HAVE_LIBMESH
    const LibMeshInit& libmesh_init = ibtkInit.getLibMeshInit();
    NULL_USE(libmesh_init);
    // Check if libMesh is initialized.
    TBOX_ASSERT(libMesh::initialized());
#endif
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "ibtk_init.output");
    // Check whether MPI is initialized
    int test;
    MPI_Initialized(&test);
    pout << "MPI " << (test ? "is " : "is not ") << "initialized.\n";
} // main
