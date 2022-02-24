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

SET(CTEST_USE_LAUNCHERS "ON" CACHE STRING "")

SET(IBAMR_ENABLE_TESTING "ON" CACHE BOOL "")

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/configure_sccache.cmake")

SET(CMAKE_C_FLAGS "-O1" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O1" CACHE STRING "C++ flags")
SET(CMAKE_Fortran_FLAGS "-O3" CACHE STRING "Fortran flags")
SET(CMAKE_INSTALL_PREFIX "/ibamr" CACHE PATH "Install destination")
SET(SAMRAI_ROOT "/samrai" CACHE PATH "Location of SAMRAI")
SET(PETSC_ROOT "/petsc/x86_64" CACHE PATH "Location of PETSc")
SET(HYPRE_ROOT "/petsc/x86_64" CACHE PATH "Location of Hypre")
SET(NUMDIFF_ROOT "/numdiff" CACHE PATH "Location of numdiff")
