## ---------------------------------------------------------------------
##
## Copyright (c) 2022 - 2022 by the IBAMR developers
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

# avoid funny problems with recursive ccache calls by providing explicit paths here
SET(CMAKE_C_COMPILER_LAUNCHER "/usr/bin/ccache" CACHE STRING "Use ccache to compile C code.")
SET(CMAKE_CXX_COMPILER_LAUNCHER "/usr/bin/ccache" CACHE STRING "Use ccache to compile C++ code.")

SET(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE STRING "C Compiler")
SET(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE STRING "C++ Compiler")
SET(CMAKE_Fortran_COMPILER "/usr/bin/gfortran" CACHE STRING "Fortran Compiler")

SET(CMAKE_C_FLAGS "-O1" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O1" CACHE STRING "C++ flags")
SET(CMAKE_Fortran_FLAGS "-O3 -Wno-unused-parameter -Wno-compare-reals" CACHE STRING "Fortran flags")
SET(CMAKE_INSTALL_PREFIX "/ibamr" CACHE PATH "Install destination")
SET(SAMRAI_ROOT "/samrai" CACHE PATH "Location of SAMRAI")
SET(PETSC_ROOT "/petsc/" CACHE PATH "Location of PETSc")
# HDF5 is handled by the module system
SET(HYPRE_ROOT "/petsc/" CACHE PATH "Location of Hypre")
SET(NUMDIFF_ROOT "/numdiff" CACHE PATH "Location of numdiff")
SET(LIBMESH_ROOT "/libmesh" CACHE PATH "Location of libmesh")
SET(LIBMESH_METHOD "OPT" CACHE STRING "Type of libmesh build (OPT or DBG)")
