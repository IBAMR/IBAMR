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

SET(IBAMR_ENABLE_TESTING "OFF" CACHE BOOL "")

# avoid funny problems with recursive ccache calls by providing explicit paths here
SET(CMAKE_C_COMPILER_LAUNCHER "/usr/bin/ccache" CACHE STRING "Use ccache to compile C code.")
SET(CMAKE_CXX_COMPILER_LAUNCHER "/usr/bin/ccache" CACHE STRING "Use ccache to compile C++ code.")

# the primary point of this configuration is to verify that we can use a custom
# MPI installation instead of the system one - here PETSc installed a copy of MPICH
SET(CMAKE_C_COMPILER "/petsc/bin/mpicc" CACHE STRING "C Compiler")
SET(CMAKE_CXX_COMPILER "/petsc/bin/mpic++" CACHE STRING "C++ Compiler")
SET(CMAKE_Fortran_COMPILER "/petsc/bin/mpifort" CACHE STRING "Fortran Compiler")

# libMesh uses std::iterator, which is deprecated, so turn off Wdeprecated-declarations
SET(CMAKE_C_FLAGS "-O0 -Wall -Wextra -Wpedantic -Werror -Wno-deprecated-declarations" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O0 -Wall -Wextra -Wpedantic -Werror -Wno-deprecated-declarations" CACHE STRING "C++ flags")
SET(CMAKE_Fortran_FLAGS "-O0 -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -Wno-compare-reals" CACHE STRING "Fortran flags")

SET(CMAKE_INSTALL_PREFIX "/ibamr" CACHE PATH "Install destination")

SET(SILO_ROOT "/petsc" CACHE PATH "Location of SAMRAI")
SET(PETSC_ROOT "/petsc" CACHE PATH "Location of PETSc")
SET(HDF5_ROOT "/petsc" CACHE PATH "Location of PETSc")
SET(HYPRE_ROOT "/petsc" CACHE PATH "Location of Hypre")
SET(SAMRAI_ROOT "/samrai" CACHE PATH "Location of SAMRAI")
SET(NUMDIFF_ROOT "/numdiff" CACHE PATH "Location of numdiff")
SET(LIBMESH_ROOT "/libmesh" CACHE PATH "Location of libmesh")
SET(LIBMESH_METHOD "OPT" CACHE STRING "Type of libmesh build (OPT or DBG)")

# do not use system dependencies - we want to verify that we can override this correctly
SET(IBAMR_FORCE_BUNDLED_Boost ON CACHE BOOL "Use the bundled version of boost")
SET(IBAMR_FORCE_BUNDLED_muParser ON CACHE BOOL "Use the bundled version of muparser")
SET(IBAMR_FORCE_BUNDLED_Eigen3 ON CACHE BOOL "Use the bundled version of eigen")
