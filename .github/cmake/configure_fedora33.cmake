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

SET(CMAKE_C_COMPILER "/usr/lib64/mpich/bin/mpicc" CACHE STRING "C Compiler")
SET(CMAKE_CXX_COMPILER "/usr/lib64/mpich/bin/mpic++" CACHE STRING "C++ Compiler")
SET(CMAKE_Fortran_COMPILER "/usr/lib64/mpich/bin/mpifort" CACHE STRING "Fortran Compiler")
SET(HDF5_ROOT "$ENV{HDF5_DIR}" CACHE PATH "Location of HDF5 binaries")

# Work around a GCC-12 bug with PatchCellDataBasicOps - GCC thinks its member
# variables are not initialized even though they have default constructors.
#
# libMesh uses std::iterator, which is deprecated, so turn off Wdeprecated-declarations
SET(CMAKE_C_FLAGS "-O1 -Wall -Wextra -Wpedantic -Werror -Wno-maybe-uninitialized -Wno-deprecated-declarations"
  CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O1 -Wall -Wextra -Wpedantic -Werror -Wno-maybe-uninitialized -Wno-deprecated-declarations" CACHE STRING "C++ flags")
SET(CMAKE_Fortran_FLAGS "-O3 -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -Wno-compare-reals" CACHE STRING "Fortran flags")

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/configure_common.cmake")
