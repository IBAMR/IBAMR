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

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/configure_ccache.cmake")

SET(CMAKE_INSTALL_PREFIX "/ibamr" CACHE PATH "Install destination")
SET(SAMRAI_ROOT "/samrai" CACHE PATH "Location of SAMRAI")
SET(LIBMESH_ROOT "/libmesh" CACHE PATH "Location of libmesh")
SET(LIBMESH_METHOD "OPT" CACHE STRING "Type of libmesh build (OPT or DBG)")
SET(PETSC_ROOT "/petsc/" CACHE PATH "Location of PETSc")
SET(HYPRE_ROOT "/petsc/" CACHE PATH "Location of Hypre")
SET(NUMDIFF_ROOT "/numdiff" CACHE PATH "Location of numdiff")
