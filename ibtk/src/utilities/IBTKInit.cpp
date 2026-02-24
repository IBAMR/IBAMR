// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
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

#include "ibtk/samrai_compatibility_names.h"
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include "SAMRAISAMRAIManager.h"
#include "SAMRAISAMRAI_MPI.h"

#include <petscsys.h>

#include <ibtk/app_namespaces.h>

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////
bool IBTKInit::s_initialized = false;

/////////////////////////////// PUBLIC ///////////////////////////////////////
IBTKInit::IBTKInit(int argc, char** argv, MPI_Comm communicator, char* petsc_file, char* petsc_help)
#ifdef IBTK_HAVE_LIBMESH
    : d_libmesh_init(argc, argv, communicator)
#endif
{
    if (s_initialized) TBOX_ERROR("IBAMR has already been initialized.\n");
#ifdef IBTK_HAVE_LIBMESH
    libMesh::ReferenceCounter::disable_print_counter_info();
    NULL_USE(petsc_file);
    NULL_USE(petsc_help);
#else
    // We need to initialize PETSc.
    PetscInitialize(&argc, &argv, petsc_file, petsc_help);
#endif
    SAMRAISAMRAIManager::setMaxNumberPatchDataEntries(2048);
    SAMRAISAMRAI_MPI::setCommunicator(communicator);
    SAMRAISAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAISAMRAIManager::startup();
    IBTK_MPI::setCommunicator(communicator);
    s_initialized = true;
}

IBTKInit::~IBTKInit()
{
    SAMRAISAMRAIManager::shutdown();
#ifndef IBTK_HAVE_LIBMESH
    PetscFinalize();
#endif
    s_initialized = false;
}
} // namespace IBTK
