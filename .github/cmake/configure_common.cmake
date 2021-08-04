## ---------------------------------------------------------------------
##
## Copyright (c) 2021 - 2021 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

set(CTEST_USE_LAUNCHERS "ON" CACHE STRING "")

set(IBAMR_ENABLE_TESTING "ON" CACHE BOOL "")

include("${CMAKE_CURRENT_LIST_DIR}/configure_sccache.cmake")

set(CMAKE_C_FLAGS "-O1" CACHE STRING "C flags")
set(CMAKE_CXX_FLAGS "-O1" CACHE STRING "C++ flags")
set(CMAKE_Fortran_FLAGS "-O3" CACHE STRING "Fortran flags")
set(CMAKE_INSTALL_PREFIX "/ibamr" CACHE PATH "Install destination")
set(SAMRAI_ROOT "/samrai" CACHE PATH "Location of SAMRAI")
set(LIBMESH_ROOT "/libmesh" CACHE PATH "Location of libmesh")
set(LIBMESH_METHOD "OPT" CACHE STRING "Type of libmesh build (OPT or DBG)")
set(PETSC_ROOT "/petsc/x86_64" CACHE PATH "Location of PetSC")
set(HYPRE_ROOT "/petsc/x86_64" CACHE PATH "Location of Hypre")
set(NUMDIFF_ROOT "/numdiff" CACHE PATH "Location of numdiff")
