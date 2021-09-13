## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2021 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
#            0: Set up some basic information about the environment            #
# ---------------------------------------------------------------------------- #
CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

PROJECT(IBAMR
  DESCRIPTION "Software infrastructure for the IB method with adaptively-refined grids"
  VERSION ${IBAMR_VERSION}
  HOMEPAGE_URL "https://ibamr.github.io"
  # include C so that we can link against C libraries (e.g., MPI::MPI_C) easily
  LANGUAGES C CXX Fortran)

MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")
INCLUDE(GNUInstallDirs)
INCLUDE(CMakePackageConfigHelpers)
include(CTest)

# Version info:
FILE(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" _version LIMIT_COUNT 1)
STRING(REGEX REPLACE "^([0-9]+)\\..*" "\\1" IBTK_VERSION_MAJOR "${_version}")
STRING(REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1" IBTK_VERSION_MINOR "${_version}")
STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" IBTK_VERSION_SUBMINOR "${_version}")
SET(IBTK_VERSION ${IBTK_VERSION_MAJOR}.${IBTK_VERSION_MINOR}.${IBTK_VERSION_SUBMINOR})
SET(IBAMR_VERSION ${IBTK_VERSION})
SET(IBAMR_DEPENDENCY_LINK_ARGUMENTS "")

# Build tests?
option(IBAMR_ENABLE_TESTING "Should tests be compiled and configured to run with ctest?" ON)

# Do we want dynamic or static linking?
OPTION(BUILD_SHARED_LIBS "Whether or not to build shared libraries." ON)
# Build the libraries in one place:
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/")
SET(IBAMR_DIMENSIONS "2" "3")

# Start building up RPATH:
SET(CMAKE_MACOSX_RPATH 1)
SET(CMAKE_BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_MACOSX_RPATH TRUE)

# determine some details about the current Fortran compiler:
INCLUDE(FortranCInterface)
IF("${FortranCInterface_GLOBAL_CASE}" STREQUAL "LOWER")
  SET(_name "name")
ELSE()
  SET(_name "NAME")
ENDIF()
STRING(JOIN " ## " IBTK_FC_FUNC ${FortranCInterface_GLOBAL_PREFIX} ${_name}
  ${FortranCInterface_GLOBAL_SUFFIX})

IF("${FortranCInterface_GLOBAL__CASE}" STREQUAL "LOWER")
  SET(_name "name")
ELSE()
  SET(_name "NAME")
ENDIF()
STRING(JOIN " ## " IBTK_FC_FUNC_ ${FortranCInterface_GLOBAL__PREFIX} ${_name}
  ${FortranCInterface_GLOBAL__SUFFIX})

# also look for a compatible _Pragma implementation:
INCLUDE(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES(
  "
  _Pragma(\"GCC diagnostic push\")
  _Pragma(\"GCC diagnostic ignored \\\\\\\"-Wunknown-pragmas\\\\\\\"\")
  _Pragma(\"GCC diagnostic ignored \\\\\\\"-Wpragmas\\\\\\\"\")
  _Pragma(\"GCC diagnostic ignored \\\\\\\"-Wextra\\\\\\\"\")
  _Pragma(\"GCC diagnostic pop\")
  int main() {}
  "
  IBTK_HAVE_PRAGMA_KEYWORD)

# other CMake configuration:
SET(CMAKE_INSTALL_MESSAGE "LAZY")
# we need this since SAMRAI is usually statically linked
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Permit some names to have lower-case letters (e.g., should users write PETSc
# or PETSC?), but for use in this script convert everything to upper case
SET(_lower_args "Boost_ROOT" "Eigen3_ROOT" "muParser_ROOT" "libMesh_ROOT"
  "PETSc_ROOT" "libMesh_METHOD" "IBAMR_FORCE_BUNDLED_Boost"
  "IBAMR_FORCE_BUNDLED_Eigen3" "IBAMR_FORCE_BUNDLED_muParser")
FOREACH(_lower_arg ${_lower_args})
  STRING(TOUPPER "${_lower_arg}" _upper_arg)
  IF((NOT "${${_lower_arg}}" STREQUAL "") AND (NOT "${${_upper_arg}}" STREQUAL ""))
    MESSAGE(FATAL_ERROR "Both ${_lower_arg} and ${_upper_arg} cannot be set. \
Their values are
      ${_lower_arg} = ${${_lower_arg}}
      ${_upper_arg} = ${${_upper_arg}}")
  ENDIF()
  IF(NOT "${${_lower_arg}}" STREQUAL "")
    MESSAGE(STATUS "Setting ${_upper_arg} = ${_lower_arg} = ${${_lower_arg}}")
    SET(${_upper_arg} "${${_lower_arg}}")
  ENDIF()
ENDFOREACH()

# print out the relevant command line arguments to make debugging easier
#
# With modern CMake, each find module uses the environment variable ending in
# _ROOT to look for an installation.
SET(_required_roots "BOOST_ROOT" "EIGEN3_ROOT" "HDF5_ROOT" "HYPRE_ROOT"
  "MUPARSER_ROOT" "MPI_ROOT" "PETSC_ROOT" "SAMRAI_ROOT")
FOREACH(_root ${_required_roots})
  IF("${${_root}}" STREQUAL "")
    MESSAGE(STATUS "${_root} was not provided to CMake: default search paths will be used.")
  ELSE()
    MESSAGE(STATUS "${_root}=${${_root}}")
  ENDIF()
ENDFOREACH()

SET(_optional_roots "GSL_ROOT" "LIBMESH_ROOT" "NUMDIFF_ROOT" "SILO_ROOT" "ZLIB_ROOT")
FOREACH(_root ${_optional_roots})
  IF("${${_root}}" STREQUAL "")
    MESSAGE(STATUS "${_root} was not provided to CMake: since this is an \
optional dependency IBAMR will be configured without this package.")
  ELSE()
    MESSAGE(STATUS "${_root}=${${_root}}")
  ENDIF()
ENDFOREACH()


# Make sure that other possible command line arguments have values set
IF("${IBAMR_FORCE_BUNDLED_BOOST}" STREQUAL "")
  SET(IBAMR_FORCE_BUNDLED_BOOST "OFF")
ENDIF()
IF("${IBAMR_FORCE_BUNDLED_EIGEN3}" STREQUAL "")
  SET(IBAMR_FORCE_BUNDLED_EIGEN3 "OFF")
ENDIF()
IF("${IBAMR_FORCE_BUNDLED_MUPARSER}" STREQUAL "")
  SET(IBAMR_FORCE_BUNDLED_MUPARSER "OFF")
ENDIF()

# ---------------------------------------------------------------------------- #
#                       1: manage mandatory dependencies                       #
# ---------------------------------------------------------------------------- #

#
# We and our dependencies require MPI so set that up first:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up MPI")
FIND_PACKAGE(MPI REQUIRED)
# If we are using the compiler wrappers then CMake may not set MPI_C_LIBRARIES -
# if its empty then try to add something in anyway
IF("${MPI_C_LIBRARIES}" STREQUAL "")
  FIND_LIBRARY(_mpi_lib NAMES mpi HINTS ${MPI_ROOT}/lib)
  IF(NOT "${_mpi_lib}" STREQUAL "_mpi_lib-NOTFOUND")
    SET(MPI_C_LIBRARIES ${_mpi_lib})
  ENDIF()
ENDIF()
MESSAGE(STATUS "MPI_C_INCLUDE_DIRS: ${MPI_C_INCLUDE_DIRS}")
MESSAGE(STATUS "MPI_C_LIBRARIES: ${MPI_C_LIBRARIES}")

#
# Boost, which may be bundled:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up Boost")
IF(${IBAMR_FORCE_BUNDLED_BOOST})
  SET(IBAMR_USE_BUNDLED_BOOST TRUE)
ELSE()
  # set up the root with the proper package name too
  SET(Boost_ROOT ${BOOST_ROOT})
  FIND_PACKAGE(Boost 1.57 QUIET)
  IF(${Boost_FOUND})
    SET(IBAMR_USE_BUNDLED_BOOST FALSE)
  ELSE()
    SET(IBAMR_USE_BUNDLED_BOOST TRUE)
  ENDIF()
ENDIF()
IF(NOT ${IBAMR_USE_BUNDLED_BOOST})
  MESSAGE(STATUS "Found external boost ${Boost_VERSION} at ${Boost_INCLUDE_DIRS}")
ELSE()
  MESSAGE(STATUS "Setting up boost as a bundled dependency")
  ADD_LIBRARY(BUNDLED_BOOST INTERFACE)
  TARGET_INCLUDE_DIRECTORIES(
    BUNDLED_BOOST
    INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/ibtk/contrib/boost>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/contrib/boost>)
  INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/ibtk/contrib/boost DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/contrib)
  INSTALL(TARGETS BUNDLED_BOOST EXPORT IBAMRTargets)
ENDIF()

#
# Eigen3, which may be bundled:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up Eigen3")
IF (${IBAMR_FORCE_BUNDLED_EIGEN3})
  SET(IBAMR_USE_BUNDLED_EIGEN3 TRUE)
ELSE()
  # set up the root with the proper package name too
  SET(Eigen3_ROOT ${EIGEN3_ROOT})
  FIND_PACKAGE(Eigen3 3.2.5 QUIET)
  IF(${Eigen3_FOUND})
    SET(IBAMR_USE_BUNDLED_EIGEN3 FALSE)
  ELSE()
    SET(IBAMR_USE_BUNDLED_EIGEN3 TRUE)
  ENDIF()
ENDIF()
IF(NOT ${IBAMR_USE_BUNDLED_EIGEN3})
  MESSAGE(STATUS "Found external Eigen3 ${Eigen3_VERSION} at ${EIGEN3_INCLUDE_DIRS}")
ELSE()
  MESSAGE(STATUS "Setting up Eigen3 as a bundled dependency")
  ADD_LIBRARY(BUNDLED_EIGEN3 INTERFACE)
  TARGET_INCLUDE_DIRECTORIES(
    BUNDLED_EIGEN3
    INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/ibtk/contrib/eigen>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/contrib/eigen>)
  INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/ibtk/contrib/eigen/Eigen DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/contrib/eigen)
  INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/ibtk/contrib/eigen/unsupported DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/contrib/eigen)
  INSTALL(TARGETS BUNDLED_EIGEN3 EXPORT IBAMRTargets)
ENDIF()

#
# muParser, which may be bundled:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up muParser")
IF(${IBAMR_FORCE_BUNDLED_MUPARSER})
  SET(IBAMR_USE_BUNDLED_MUPARSER TRUE)
ELSE()
  FIND_PATH(MUPARSER_INCLUDE_DIRS NAMES muParser.h HINTS ${MUPARSER_ROOT}/include /usr/include/)
  FIND_LIBRARY(MUPARSER_LIBRARIES NAMES muparser HINTS ${MUPARSER_ROOT}/lib /usr/lib)
  IF("${MUPARSER_INCLUDE_DIRS}" STREQUAL "MUPARSER_INCLUDE_DIRS-NOTFOUND" OR
      "${MUPARSER_LIBRARIES}" STREQUAL "MUPARSER_LIBRARIES-NOTFOUND")
    SET(IBAMR_USE_BUNDLED_MUPARSER TRUE)
  ELSE()
    SET(IBAMR_USE_BUNDLED_MUPARSER FALSE)
  ENDIF()
ENDIF()
IF(${IBAMR_USE_BUNDLED_MUPARSER})
  MESSAGE(STATUS "Setting up muParser as a bundled dependency")
  ADD_LIBRARY(BUNDLED_MUPARSER)
  TARGET_INCLUDE_DIRECTORIES(
    BUNDLED_MUPARSER
    PUBLIC
    $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/ibtk/contrib/muparser/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/contrib/muparser>)
  SET_PROPERTY(TARGET BUNDLED_MUPARSER PROPERTY CXX_STANDARD 11)
  TARGET_COMPILE_FEATURES(BUNDLED_MUPARSER PUBLIC cxx_std_11)

  TARGET_SOURCES(BUNDLED_MUPARSER
    PRIVATE
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParser.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserBase.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserBytecode.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserCallback.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserDLL.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserError.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserInt.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserTest.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserTokenReader.cpp)

  INSTALL(FILES
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParser.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserBase.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserBytecode.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserCallback.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserDLL.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserDef.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserError.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserFixes.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserInt.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserStack.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserTemplateMagic.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserTest.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserToken.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserTokenReader.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/contrib/muparser)

  INSTALL(TARGETS BUNDLED_MUPARSER EXPORT IBAMRTargets COMPONENT library)
ELSE()
  MESSAGE(STATUS "Found external muParser")
  MESSAGE(STATUS "MUPARSER_INCLUDE_DIRS: ${MUPARSER_INCLUDE_DIRS}")
  MESSAGE(STATUS "MUPARSER_LIBRARIES: ${MUPARSER_LIBRARIES}")
ENDIF()

#
# HDF5:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up HDF5")
SET(HDF5_FIND_DEBUG TRUE)
FIND_PACKAGE(HDF5 REQUIRED COMPONENTS C)
IF("${HDF5_LIBRARIES}" STREQUAL "HDF5_LIBRARIES-NOTFOUND")
  MESSAGE(FATAL_ERROR "Unable to find a valid HDF5 installation.")
ENDIF()

#
# hypre:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up HYPRE")
FIND_PATH(HYPRE_INCLUDE_DIRS REQUIRED NAMES HYPRE.h HINTS ${HYPRE_ROOT}/include /usr/include/hypre)
MESSAGE(STATUS "HYPRE_INCLUDE_DIRS: ${HYPRE_INCLUDE_DIRS}")
FIND_LIBRARY(HYPRE_LIBRARIES REQUIRED NAMES HYPRE HINTS ${HYPRE_ROOT}/lib /usr/lib)
MESSAGE(STATUS "HYPRE_LIBRARIES: ${HYPRE_LIBRARIES}")
IF(${HYPRE_INCLUDE_DIRS} STREQUAL "HYPRE_INCLUDE_DIRS-NOTFOUND" OR
    ${HYPRE_LIBRARIES} STREQUAL "HYPRE_LIBRARIES-NOTFOUND")
  MESSAGE(FATAL_ERROR "Unable to find a valid HYPRE installation.")
ENDIF()

#
# SAMRAI:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up SAMRAI")
SET(_samrai_library_suffixes "algs" "appu" "geom" "hier"
  "math_std" "mesh" "pdat_std" "solv" "xfer")

# samrai version info:
FILE(STRINGS "${SAMRAI_ROOT}/include/SAMRAI_config.h" SAMRAI_VERSION_MAJOR_LINE
  REGEX "#define.*SAMRAI_VERSION_MAJOR")
STRING(REGEX REPLACE "^.*SAMRAI_VERSION_MAJOR.* ([0-9]+).*" "\\1"
  SAMRAI_VERSION_MAJOR "${SAMRAI_VERSION_MAJOR_LINE}"
  )
FILE(STRINGS "${SAMRAI_ROOT}/include/SAMRAI_config.h" SAMRAI_VERSION_MINOR_LINE
  REGEX "#define.*SAMRAI_VERSION_MINOR")
STRING(REGEX REPLACE "^.*SAMRAI_VERSION_MINOR.* ([0-9]+).*" "\\1"
  SAMRAI_VERSION_MINOR "${SAMRAI_VERSION_MINOR_LINE}"
  )
FILE(STRINGS "${SAMRAI_ROOT}/include/SAMRAI_config.h" SAMRAI_VERSION_PATCHLEVEL_LINE
  REGEX "#define.*SAMRAI_VERSION_PATCHLEVEL")
STRING(REGEX REPLACE "^.*SAMRAI_VERSION_PATCHLEVEL.* ([0-9]+).*" "\\1"
  SAMRAI_VERSION_PATCHLEVEL "${SAMRAI_VERSION_PATCHLEVEL_LINE}"
  )
SET(SAMRAI_VERSION_STRING "${SAMRAI_VERSION_MAJOR}.${SAMRAI_VERSION_MINOR}\
.${SAMRAI_VERSION_PATCHLEVEL}")
MESSAGE(STATUS "Found SAMRAI ${SAMRAI_VERSION_STRING} at ${SAMRAI_ROOT}")

# samrai libs:
SET(SAMRAI_INCLUDE_DIRS "${SAMRAI_ROOT}/include")
SET(SAMRAI2d_LIBRARIES)
SET(SAMRAI3d_LIBRARIES)
ADD_LIBRARY(SAMRAI UNKNOWN IMPORTED)
FIND_LIBRARY(SAMRAI_path SAMRAI REQUIRED HINTS "${SAMRAI_ROOT}/lib")
LIST(APPEND SAMRAI2d_LIBRARIES ${SAMRAI_path})
LIST(APPEND SAMRAI3d_LIBRARIES ${SAMRAI_path})

FOREACH(_d ${IBAMR_DIMENSIONS})
  FOREACH(_suffix ${_samrai_library_suffixes})
    SET(_lib_name "SAMRAI${_d}d_${_suffix}")
    ADD_LIBRARY( ${_lib_name} UNKNOWN IMPORTED)

    FIND_LIBRARY("${_lib_name}_path" ${_lib_name} REQUIRED
      HINTS "${SAMRAI_ROOT}/lib")
    IF("${${_lib_name}_path}" STREQUAL "${_lib_name}-NOTFOUND")
      MESSAGE(FATAL_ERROR
        "Unable to find required library ${_lib_name} in directory ${SAMRAI_ROOT}/lib")
    ENDIF()

    LIST(APPEND "SAMRAI${_d}d_LIBRARIES" "${${_lib_name}_path}")
  ENDFOREACH()
ENDFOREACH()

#
# PETSc:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up PETSc")
FIND_FILE(PETSC_VARIABLES_FILE petscvariables HINTS ${PETSC_ROOT}
  PATH_SUFFIXES conf lib/petsc/conf)
IF(${PETSC_VARIABLES_FILE} STREQUAL "PETSC_VARIABLES_FILE-NOTFOUND")
  MESSAGE(FATAL_ERROR "unable to find the petscvariables configuration file")
ELSE()
  FILE(STRINGS ${PETSC_VARIABLES_FILE} _petsc_raw_includes REGEX "^PETSC_CC_INCLUDES =.*")
  IF ("${_petsc_raw_includes}" STREQUAL "")
    MESSAGE(FATAL_ERROR
      "The configuration script was unable to find the list of                 \
PETSc include directories in the file ${PETSC_VARIABLES_FILE}. This usually    \
indicates that PETSC_ROOT was set to PETSC_DIR when it should be set to        \
PETSC_DIR/PETSC_ARCH (i.e., PETSC_ROOT must be set to the complete path to the \
PETSc installation).")
  ENDIF()

  # We do not yet support PETSc with 64-bit integers
  #
  # TODO - we can use REQUIRED in FIND_FILE in CMake 3.18 and newer
  FIND_FILE(PETSC_CONF_FILE petscconf.h HINTS ${PETSC_ROOT}
    PATH_SUFFIXES petsc include include/petsc)
  IF(${PETSC_VARIABLES_FILE} STREQUAL "PETSC_CONF_FILE-NOTFOUND")
    MESSAGE(FATAL_ERROR "unable to find the petscconf.h configuration file")
  ELSE()
    FILE(STRINGS ${PETSC_CONF_FILE} _petsc_have_64bit_indices REGEX
      "^ *# *define.*PETSC_USE_64BIT_INDICES *1")
    IF(NOT ${_petsc_have_64bit_indices} STREQUAL "")
      MESSAGE(FATAL_ERROR "IBAMR does not support using 64 bit indices with \
PETSc at this time. Please recompile PETSc (and any dependencies, such as   \
libMesh, it may have) to NOT use 64 bit indices.")
    ENDIF()
  ENDIF()

  STRING(REGEX REPLACE "^PETSC_CC_INCLUDES =(.*)" "\\1" _petsc_raw_includes ${_petsc_raw_includes})
  SEPARATE_ARGUMENTS(_petsc_raw_includes)
  # Get rid of preceding -Is (CMake wants just directory names):
  FOREACH(_include ${_petsc_raw_includes})
    STRING(REGEX REPLACE "^-I" "" _directory "${_include}")
    LIST(APPEND PETSC_INCLUDE_DIRS ${_directory})
  ENDFOREACH()

  FILE(STRINGS ${PETSC_VARIABLES_FILE} PETSC_LIBRARIES REGEX "^PETSC_WITH_EXTERNAL_LIB =.*")
  STRING(REGEX REPLACE "^PETSC_WITH_EXTERNAL_LIB =(.*)" "\\1" PETSC_LIBRARIES ${PETSC_LIBRARIES})
  # If PETSC is not installed then ${PETSC_DIR}/${PETSC_ARCH} may appear in the list
  STRING(REGEX REPLACE "\\$.PETSC_DIR./\\$.PETSC_ARCH." "${PETSC_ROOT}" PETSC_LIBRARIES ${PETSC_LIBRARIES})
  SEPARATE_ARGUMENTS(PETSC_LIBRARIES)
  FIND_LIBRARY(PETSC_LIBRARY REQUIRED NAMES "petsc" HINTS ${PETSC_ROOT}/lib)
  LIST(PREPEND PETSC_LIBRARIES ${PETSC_LIBRARY})
  MESSAGE(STATUS "PETSC LIBRARIES: ${PETSC_LIBRARIES}")

  # extract the version numbers:
  FIND_FILE(PETSC_VERSION_FILE petscversion.h HINTS ${PETSC_INCLUDE_DIRS})
  FILE(STRINGS "${PETSC_VERSION_FILE}" PETSC_VERSION_MAJOR_STRING
    REGEX "#define.*PETSC_VERSION_MAJOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_MAJOR.* ([0-9]+).*" "\\1"
    PETSC_VERSION_MAJOR "${PETSC_VERSION_MAJOR_STRING}"
    )
  FILE(STRINGS "${PETSC_VERSION_FILE}" PETSC_VERSION_MINOR_STRING
    REGEX "#define.*PETSC_VERSION_MINOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_MINOR.* ([0-9]+).*" "\\1"
    PETSC_VERSION_MINOR "${PETSC_VERSION_MINOR_STRING}"
    )
  FILE(STRINGS "${PETSC_VERSION_FILE}" PETSC_VERSION_SUBMINOR_STRING
    REGEX "#define.*PETSC_VERSION_SUBMINOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_SUBMINOR.* ([0-9]+).*" "\\1"
    PETSC_VERSION_SUBMINOR "${PETSC_VERSION_SUBMINOR_STRING}"
    )
  SET(PETSC_VERSION
    "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}")
  IF(${PETSC_VERSION} VERSION_LESS 3.7.0)
    MESSAGE(FATAL_ERROR
      "IBAMR requires PETSc version 3.7.0 or newer but the version provided at
${PETSC_ROOT} is version ${PETSC_VERSION}")
  ENDIF()
ENDIF()

# ---------------------------------------------------------------------------- #
#                       2: manage optional dependencies                        #
# ---------------------------------------------------------------------------- #
MESSAGE(STATUS "")
SET(IBAMR_HAVE_LIBMESH FALSE)
IF(NOT "${LIBMESH_ROOT}" STREQUAL "")
  MESSAGE(STATUS "Setting up libMesh")
  FIND_PROGRAM(_libmesh_config "libmesh-config" HINTS ${LIBMESH_ROOT}/bin/)
  IF("${_libmesh_config}" STREQUAL "_libmesh_config-NOTFOUND")
    MESSAGE(FATAL_ERROR "\
libMesh (an optional dependency) was specified with LIBMESH_ROOT=${LIBMESH_ROOT} \
but a valid libmesh-config script could not be found in ${LIBMESH_ROOT}/bin/. \
Please check the value of LIBMESH_ROOT and rerun CMake.")
  ELSE()
    SET(IBAMR_HAVE_LIBMESH TRUE)
    # libMesh requires that we set an environment variable to query libmesh-config
    IF("${LIBMESH_METHOD}" STREQUAL "")
      MESSAGE(FATAL_ERROR "If LIBMESH_ROOT is provided to CMake then \
LIBMESH_METHOD must also be provided to CMake.")
    ENDIF()
    STRING(TOLOWER ${LIBMESH_METHOD} _lower_method)
    SET(ENV{METHOD} ${_lower_method})
    EXECUTE_PROCESS(COMMAND "${_libmesh_config}" "--version"
      OUTPUT_VARIABLE LIBMESH_VERSION)
    STRING(REGEX REPLACE "\n$" "" LIBMESH_VERSION ${LIBMESH_VERSION})
    IF(${LIBMESH_VERSION} VERSION_LESS 1.1.0)
      MESSAGE(FATAL_ERROR
        "IBAMR requires libMesh version 1.1.0 or newer but the version provided at
${LIBMESH_ROOT} is version ${LIBMESH_VERSION}")
    ENDIF()

    # Extract libraries and convert to a CMake list:
    EXECUTE_PROCESS(COMMAND "${_libmesh_config}" "--libs" OUTPUT_VARIABLE LIBMESH_LIBRARIES)
    STRING(REGEX REPLACE "\n$" "" LIBMESH_LIBRARIES ${LIBMESH_LIBRARIES})
    SEPARATE_ARGUMENTS(LIBMESH_LIBRARIES)
    FIND_LIBRARY(LIBMESH_LIBRARY REQUIRED NAMES "mesh_${_lower_method}" HINTS ${LIBMESH_ROOT}/lib)
    IF("${LIBMESH_LIBRARY}" STREQUAL "${LIBMESH_LIBRARY}-NOTFOUND")
      MESSAGE(FATAL_ERROR "
Unable to find the libMesh library itself. This usually happens when either \
libMesh is installed incorrectly or the wrong LIBMESH_METHOD was given to \
CMake.")
    ENDIF()
    LIST(PREPEND LIBMESH_LIBRARIES ${LIBMESH_LIBRARY})
    MESSAGE(STATUS "LIBMESH_LIBRARIES: ${LIBMESH_LIBRARIES}")

    EXECUTE_PROCESS(COMMAND "${_libmesh_config}" "--include" OUTPUT_VARIABLE _libmesh_raw_includes)
    STRING(REGEX REPLACE "\n$" "" _libmesh_raw_includes ${_libmesh_raw_includes})
    SEPARATE_ARGUMENTS(_libmesh_raw_includes)
    # Get rid of preceding -Is (CMake wants just directory names):
    SET(LIBMESH_INCLUDE_DIRS)
    FOREACH(_include ${_libmesh_raw_includes})
      STRING(REGEX REPLACE "^-I" "" _directory "${_include}")
      LIST(APPEND LIBMESH_INCLUDE_DIRS ${_directory})
    ENDFOREACH()
    MESSAGE(STATUS "LIBMESH_INCLUDE_DIRS: ${LIBMESH_INCLUDE_DIRS}")

    # we won't use it directly but see if libMesh expects C++14 by checking its C++ flags:
    EXECUTE_PROCESS(COMMAND "${_libmesh_config}" "--cxxflags" OUTPUT_VARIABLE _libmesh_cxxflags)
    STRING(REGEX REPLACE "\n$" " " _libmesh_cxxflags ${_libmesh_cxxflags})
    SEPARATE_ARGUMENTS(_libmesh_cxxflags)
    ## CXX_FLAGS is not comma-separated so pull things out again
    STRING(REPLACE ";" " " ${_libmesh_cxxflags} CMAKE_REQUIRED_FLAGS)
    SET(CMAKE_REQUIRED_INCLUDES ${LIBMESH_INCLUDE_DIRS})

    # figure out which version of C++ libMesh wants to use:
    SET(LIBMESH_CXX_VERSION 11) # we require C++11 anyway so use it as the default
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <libmesh/libmesh_config.h>
      // older versions do not define LIBMESH_HAVE_CXX14
      #ifdef LIBMESH_HAVE_CXX14_MAKE_UNIQUE
      // OK
      #else
      #error
      #endif
      int main() {}
      "
      LIBMESH_WITH_CXX14
      )
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <libmesh/libmesh_config.h>
      #ifdef LIBMESH_HAVE_CXX17
      // OK
      #else
      #error
      #endif
      int main() {}
      "
      LIBMESH_WITH_CXX17
      )

    IF(${LIBMESH_WITH_CXX17})
      SET(LIBMESH_CXX_VERSION 17)
    ELSEIF(${LIBMESH_WITH_CXX14})
      SET(LIBMESH_CXX_VERSION 14)
    ENDIF()

    # We are not compatible with libMesh's bundled version of boost, so make
    # sure that doesn't happen:
    IF(EXISTS ${LIBMESH_ROOT}/include/boost/)
      IF(EXISTS ${LIBMESH_ROOT}/include/boost/multi_array.hpp)
        # this is OK: perhaps libMesh was installed in /usr/
      ELSE()
        MESSAGE(FATAL_ERROR "\
The directory ${LIBMESH_ROOT}/include/boost/ exists and contains a boost \
installation that does not provide all the headers that IBAMR needs. This can \
happen when either libMesh is installed with its own bundled version of boost \
(which is not compatible with IBAMR) or when libMesh is installed into a \
directory which happens to contain a copy of boost (e.g., when libMesh is \
recompiled without boost and installed into the same location as a previous \
copy of libMesh). If you want to use libMesh with boost then both libMesh and \
IBAMR must use the same external copy of boost. The best way to fix this \
problem is to delete ${LIBMESH_ROOT} and reinstall libMesh with either no \
boost support or an external boost library.")
      ENDIF()
    ENDIF()

    # Use slightly different error messages for bundled and non-bundled Eigen
    IF(EXISTS ${LIBMESH_ROOT}/include/Eigen)
      IF("${EIGEN3_ROOT}" STREQUAL "${LIBMESH_ROOT}/include/Eigen")
        # OK
      ELSEIF(${IBAMR_USE_BUNDLED_EIGEN3})
        MESSAGE(FATAL_ERROR "\
libMesh appears to have been compiled with a bundled copy of Eigen3 (i.e., the \
directory ${LIBMESH_ROOT}/include/Eigen exists) but IBAMR was configured to \
use its own bundled copy of Eigen3. IBAMR does not support using multiple copies \
of Eigen3. To fix this, we recommend deleting the current libMesh installation, \
recompiling libMesh without Eigen support, and then reinstalling libMesh. IBAMR \
does not support using libMesh's bundled copy of Eigen3 since that version is \
not installed with Eigen's CMake configuration files. If you need to use Eigen \
with libMesh then you will need to use an external installation for both \
libMesh and IBAMR.")
      ELSE()
        MESSAGE(FATAL_ERROR "\
libMesh appears to have been compiled with a bundled copy of Eigen3 (i.e., the \
directory ${LIBMESH_ROOT}/include/Eigen exists) and this version of Eigen is \
different from the one provided to CMake via EIGEN3_ROOT = ${EIGEN3_ROOT}: i.e., \
there are two versions of Eigen3 in the configuration path, which is not \
supported. To fix this, delete the current libMesh installation, recompile \
libMesh without Eigen, and reinstall libMesh. If you aren't sure what to do, \
the best option is to ensure that libMesh is compiled without support for \
Eigen and then let IBAMR either detect an Eigen installation or use its own \
bundled copy.")
      ENDIF()
    ENDIF()

    # Verify that libMesh uses the same PETSc that we do:
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <libmesh/libmesh_config.h>
      #ifdef LIBMESH_HAVE_PETSC
      // OK
      #else
      #error
      #endif
      int main() {}
      "
      LIBMESH_WITH_PETSC
      )
    IF(NOT "${LIBMESH_WITH_PETSC}")
      MESSAGE(FATAL_ERROR "IBAMR requires that libMesh be compiled with PETSc.")
    ENDIF()
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <libmesh/libmesh_config.h>
      #if LIBMESH_DETECTED_PETSC_VERSION_MAJOR != ${PETSC_VERSION_MAJOR}
      #error
      #endif
      #if LIBMESH_DETECTED_PETSC_VERSION_MINOR != ${PETSC_VERSION_MINOR}
      #error
      #endif
      #if LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR != ${PETSC_VERSION_SUBMINOR}
      #error
      #endif
      int main() {}
      "
      LIBMESH_WITH_SAME_PETSC
      )
    IF(NOT "${LIBMESH_WITH_SAME_PETSC}")
      MESSAGE(FATAL_ERROR "\
The version of PETSc detected by libMesh differs from the version of PETSc \
detected by IBAMR. This is not allowed.")
    ENDIF()

    # Verify that libMesh uses MPI:
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <libmesh/libmesh_config.h>
      #ifdef LIBMESH_HAVE_MPI
      // OK
      #else
      #error
      #endif
      int main() {}
      "
      LIBMESH_WITH_MPI
      )
    IF(NOT "${LIBMESH_WITH_MPI}")
      MESSAGE(FATAL_ERROR "IBAMR requires that libMesh be compiled with MPI.")
    ENDIF()

    # Verify that libMesh uses XDR (we need it for restarts and the test suite):
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <libmesh/libmesh_config.h>
      #ifdef LIBMESH_HAVE_XDR
      // OK
      #else
      #error
      #endif
      int main() {}
      "
      LIBMESH_WITH_XDR
      )
    IF(NOT "${LIBMESH_WITH_XDR}")
      MESSAGE(FATAL_ERROR "IBAMR requires that libMesh be compiled with XDR.")
    ENDIF()

    # Note that we don't technically require ExodusII - even though we have
    # AppInitializer::getExodusIIFilename() we never actually use ExodusII APIs.

    SET(CMAKE_REQUIRED_FLAGS "")
    SET(CMAKE_REQUIRED_INCLUDES "")
  ENDIF()
ELSE()
  MESSAGE(STATUS "LIBMESH_ROOT was not specified so IBAMR will be configured without it.")
ENDIF()

MESSAGE(STATUS "")
SET(IBAMR_HAVE_SILO FALSE)
IF(NOT "${SILO_ROOT}" STREQUAL "")
  MESSAGE(STATUS "Setting up Silo")
  FIND_PATH(SILO_INCLUDE_DIRS NAMES silo.h HINTS ${SILO_ROOT}/include /usr/include/)
  FIND_LIBRARY(SILO_LIBRARIES NAMES silo siloh5 HINTS ${SILO_ROOT}/lib /usr/lib)
  IF(NOT "${SILO_LIBRARIES}" STREQUAL "SILO_LIBRARIES-NOTFOUND" AND
      NOT "${SILO_INCLUDE_DIRS}" STREQUAL "SILO_INCLUDE_DIRS-NOTFOUND")
    MESSAGE(STATUS "SILO_INCLUDE_DIRS: ${SILO_INCLUDE_DIRS}")
    MESSAGE(STATUS "SILO_LIBRARIES: ${SILO_LIBRARIES}")
    SET(IBAMR_HAVE_SILO TRUE)

    # Silo depends on zlib
    FIND_PACKAGE(ZLIB REQUIRED)
  ELSE()
    MESSAGE(FATAL_ERROR "\
Silo (an optional dependency) was specified with SILO_ROOT=${SILO_ROOT} but a \
silo installation could not be found in that location. Please check the value \
of SILO_ROOT and rerun CMake.")
  ENDIF()
ELSE()
  MESSAGE(STATUS "SILO_ROOT was not specified so IBAMR will be configured without it.")
ENDIF()

MESSAGE(STATUS "")
SET(IBAMR_HAVE_GSL FALSE)
IF(NOT "${GSL_ROOT}" STREQUAL "")
  MESSAGE(STATUS "Setting up GSL")
  # GSL has not yet been updated to the modern convention so help it
  SET(GSL_ROOT_DIR GSL_ROOT)
  FIND_PACKAGE(GSL)
  IF(GSL_FOUND)
    SET(IBAMR_HAVE_GSL TRUE)
  ENDIF()
ENDIF()

# ---------------------------------------------------------------------------- #
#                 3: Check for conflicts between dependencies                  #
# ---------------------------------------------------------------------------- #

#
# Helper function that checks that we can compile and link _src with the given
# dependency libraries and dependency includes against the copy of MPI we found
# - i.e., this function verifies that a dependency uses the same MPI as the one
# we just found. As an extra check we use PETSc's MPI consistency checks in
# petscsys.h to verify that we do have the same MPI version at the preprocessor
# step.
#
FUNCTION(IBAMR_CHECK_COMPILATION_WITH_MPI _dependency_name _src _dependency_libraries
    _dependency_includes)
  SET(CMAKE_REQUIRED_INCLUDES)
  LIST(APPEND CMAKE_REQUIRED_INCLUDES "${_dependency_includes}")
  LIST(APPEND CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDE_DIRS}")
  LIST(APPEND CMAKE_REQUIRED_INCLUDES "${MPI_C_INCLUDE_DIRS}")
  SET(CMAKE_REQUIRED_LIBRARIES)
  LIST(APPEND CMAKE_REQUIRED_LIBRARIES "${_dependency_libraries}")
  LIST(APPEND CMAKE_REQUIRED_LIBRARIES "${PETSC_LIBRARIES}")
  LIST(APPEND CMAKE_REQUIRED_LIBRARIES "${MPI_C_LIBRARIES}")
  LIST(APPEND CMAKE_REQUIRED_LIBRARIES "${MPI_CXX_LIBRARIES}")

  SET(_test_name "${_dependency_name}_SAME_MPI_COMPILE_TEST")
  CHECK_CXX_SOURCE_COMPILES("${_src}" ${_test_name})
  IF(NOT "${${_test_name}}")
    MESSAGE(FATAL_ERROR "\
Unable to compile a test program dependent on ${_dependency_name} that links \
against MPI. This usually means that the dependency is misconfigured or is \
using a different version of MPI than the one supplied to IBAMR.")
  ENDIF()
  SET(CMAKE_REQUIRED_INCLUDES)
  SET(CMAKE_REQUIRED_LIBRARIES)

  FIND_PROGRAM(_ldd "ldd")
  FIND_PROGRAM(_readlink "readlink")
  IF(NOT "${_ldd}" STREQUAL "_ldd-NOTFOUND" AND NOT "${_readlink}" STREQUAL "_readlink-NOTFOUND")
    MESSAGE(STATUS "\
Verifying that ${_dependency_name} and IBAMR use the same MPI implementation")
    FOREACH(_lib ${_dependency_libraries})
      IF(EXISTS ${_lib})
        EXECUTE_PROCESS(COMMAND "${_ldd}" "${_lib}"
          OUTPUT_VARIABLE _libs ERROR_QUIET)
        SEPARATE_ARGUMENTS(_libs)
        FOREACH(_dep ${_libs})
          STRING(REGEX MATCH "libmpi\." _has_mpi ${_dep})
          IF("${_has_mpi}" STREQUAL "libmpi\." AND EXISTS "${_dep}")
            # We found an MPI library: first resolve symbolic links and then
            # check for equality
            SET(_found FALSE)
            EXECUTE_PROCESS(COMMAND "${_readlink}" "-f" "${_dep}"
              OUTPUT_VARIABLE _resolved_dep)
            FOREACH(_mpi_lib "${MPI_C_LIBRARIES}")
              EXECUTE_PROCESS(COMMAND "${_readlink}" "-f" "${_mpi_lib}"
                OUTPUT_VARIABLE _resolved_mpi_lib)
              IF("${_resolved_mpi_lib}" STREQUAL "${_resolved_dep}")
                SET(_found TRUE)
                BREAK()
              ENDIF()
            ENDFOREACH()
            IF(NOT ${_found})
              MESSAGE(FATAL_ERROR "\
The library
    ${_lib}
a part of ${_dependency_name}, links against
    ${_dep}
which conflicts with the MPI implementation explicitly provided to IBAMR, which includes
    ${MPI_C_LIBRARIES}
Please recompile ${_dependency} to link against the same version of MPI.")
            ENDIF()
          ENDIF()
        ENDFOREACH()
      ENDIF()
    ENDFOREACH()
    MESSAGE(STATUS "\
Verifying that ${_dependency_name} and IBAMR use the same MPI implementation - Success")
  ENDIF()
ENDFUNCTION()

#
# PETSc:
#
IBAMR_CHECK_COMPILATION_WITH_MPI(PETSC
  "\
#include <mpi.h>
#include <petscvec.h>

int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, NULL, NULL);

    Vec vec;
    double values[2] = {42.0, 42.0};
    VecCreateGhostWithArray(MPI_COMM_WORLD, 2, PETSC_DECIDE, 0, NULL,
                            values, &vec);
    VecDestroy(&vec);

    PetscFinalize();
}"
  "${PETSC_LIBRARIES}"
  "${PETSC_INCLUDE_DIRS}")

#
# HDF5:
#
IBAMR_CHECK_COMPILATION_WITH_MPI(HDF5
  "\
#include <hdf5.h>
#include <mpi.h>
#include <petscsys.h> // this is the preprocessor check

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    hid_t  file_id = H5Fcreate(\"file.h5\", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    (void)file_id;
    herr_t status  = H5Fclose(file_id);
    (void)status;
    MPI_Finalize();
}"
  "${HDF5_LIBRARIES}"
  "${HDF5_INCLUDE_DIRS}")

#
# HYPRE:
#
IBAMR_CHECK_COMPILATION_WITH_MPI(HYPRE
  "\
#include <HYPRE_struct_ls.h>
#include <mpi.h>
#include <petscsys.h>

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    HYPRE_StructGrid grid;
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);
    HYPRE_StructGridAssemble(grid);
    HYPRE_StructGridDestroy(grid);
    MPI_Finalize();
}"
  "${HYPRE_LIBRARIES}"
  "${HYPRE_INCLUDE_DIRS}")

# SAMRAI is statically linked and doesn't export its linkage information in a
# way that's easy to parse so we cannot execute the MPI check

#
# libMesh:
#
IBAMR_CHECK_COMPILATION_WITH_MPI(LIBMESH
  "\
#include <mpi.h>
#include <petscsys.h>
// This test is a little different: since libMesh may use C++11 C++14, or C++17,
// but we but have not yet configured that flag yet, we instead use libMesh's
// petsc dependency to check that we have the same MPI version.
int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, NULL, NULL);
    PetscFinalize();
}"
  "${LIBMESH_LIBRARIES}"
  "${LIBMESH_INCLUDE_DIRS}")

# ---------------------------------------------------------------------------- #
#                       4: IBAMR-specific configuration                        #
# ---------------------------------------------------------------------------- #

# generate the configuration header. Define a few more things for IBTK.
SET(IBTK_HAVE_LIBMESH ${IBAMR_HAVE_LIBMESH})
SET(IBTK_HAVE_SILO ${IBAMR_HAVE_SILO})
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/ibtk/include/ibtk/config.h.in
  ${CMAKE_BINARY_DIR}/ibtk/include/ibtk/config.h)
INSTALL(FILES ${CMAKE_BINARY_DIR}/ibtk/include/ibtk/config.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ibtk/)

INCLUDE(CMakeParseArguments)
#
# Macro to process m4 files to generate Fortran.
#
# This macro includes some ad-hoc logic to attempt to detect M4 dependencies
# (i.e., include() statements in the provided file) so that the build system is
# aware of them. However, since this is done with regular expressions it is
# fragile and may not always pick up file dependencies. For example, writing
#
# input(
# dnl skip this line
# /path/to/file.m4
# )
#
# is perfectly valid M4 but beyond the capacity of this macro to parse. In such
# cases we simply do not add the dependency and print a warning - compilation
# from scratch will still succeed but the build system will not be able to track
# the dependency.
#
# Required arguments:
#
# - NAME: Relative path (from the directory in which the macro is called) to the
#   relevant Fortran file. For example: if we want to preprocess
#   foo/bar/baz.f.m4 then we should pass foo/bar/baz.f as the NAME argument.
MACRO(IBAMR_PROCESS_M4)
  SET(options)
  SET(singleValueArgs NAME)

  CMAKE_PARSE_ARGUMENTS(arg "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})
  GET_FILENAME_COMPONENT(_in_directory "${CMAKE_CURRENT_SOURCE_DIR}/${arg_NAME}.m4" DIRECTORY)
  GET_FILENAME_COMPONENT(_out_directory "${CMAKE_CURRENT_BINARY_DIR}/${arg_NAME}.f" DIRECTORY)
  FILE(MAKE_DIRECTORY ${_out_directory})

  SET(_samrai_fortdir ${SAMRAI_ROOT}/include)
  SET(_current_srcdir ${_in_directory})
  SET(_top_srcdir ${CMAKE_SOURCE_DIR})
  SET(_input "${CMAKE_CURRENT_SOURCE_DIR}/${arg_NAME}.m4")
  SET(_output "${CMAKE_CURRENT_BINARY_DIR}/${arg_NAME}")
  SET(_args "-DSAMRAI_FORTDIR=${_samrai_fortdir}"
    "-DCURRENT_SRCDIR=${_current_srcdir}"
    "-DTOP_SRCDIR=${_top_srcdir}")

  FILE(STRINGS ${_input} _parsed_includes REGEX "include.*DIR")
  SET(_expanded_files)
  FOREACH(_parsed_include ${_parsed_includes})
    STRING(REGEX REPLACE "include\\((.*)\\)(dnl)?$" "\\1" _tmp ${_parsed_include})
    STRING(REGEX REPLACE "SAMRAI_FORTDIR" "${_samrai_fortdir}" _tmp ${_tmp})
    STRING(REGEX REPLACE "CURRENT_SRCDIR" "${_current_srcdir}" _tmp ${_tmp})
    STRING(REGEX REPLACE "TOP_SRCDIR" "${_top_srcdir}" _tmp ${_tmp})
    LIST(APPEND _expanded_files ${_tmp})
  ENDFOREACH()

  SET(_dependencies)
  LIST(APPEND _dependencies ${_input})
  FOREACH(_expanded_file ${_expanded_files})
    IF (EXISTS ${_expanded_file})
      LIST(APPEND _dependencies ${_expanded_file})
    ELSE()
      MESSAGE(WARNING
"Unable to locate file\n    ${_expanded_file}\n which is a dependency of
${_input_file}. Since dependency tracking of M4 files is done manually this is
not a fatal problem but it may result in compilation failures later in the build process.")
    ENDIF()
  ENDFOREACH()

  ADD_CUSTOM_COMMAND(
    OUTPUT ${_output}
    DEPENDS ${_dependencies}
    COMMAND m4 ${_args} ${_input} > ${_output}
    VERBATIM)

  SET_SOURCE_FILES_PROPERTIES(${_output} PROPERTIES GENERATED true)
  SET_SOURCE_FILES_PROPERTIES(${_output} PROPERTIES LANGUAGE Fortran)
ENDMACRO()

# Macro to setup an IBAMR target library with all common features (i.e.,
# everything but the source files)
FUNCTION(IBAMR_SETUP_TARGET_LIBRARY target_library)
  MESSAGE(STATUS "setting up target ${target_library}")

  SET_TARGET_PROPERTIES(${target_library}
    PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH TRUE)

  # Figure out if we are a 2D or a 3D library:
  SET(_2_location = -1)
  SET(_3_location = -1)
  STRING(FIND "${target_library}" "2" _2_location)
  STRING(FIND "${target_library}" "3" _3_location)
  SET(_d "")
  IF(NOT ${_2_location} STREQUAL "-1")
    SET(_d "2")
  ELSE()
    SET(_d "3")
  ENDIF()

  MESSAGE(STATUS "Adding flag -DNDIM=${_d} to target ${target_library}")
  TARGET_COMPILE_OPTIONS(${target_library} PUBLIC -DNDIM=${_d})

  IF(${IBAMR_HAVE_LIBMESH})
    MESSAGE(STATUS "libMesh expects C++${LIBMESH_CXX_VERSION}, so ${target_library} will be built with C++${LIBMESH_CXX_VERSION} as a requirement.")
    SET_PROPERTY(TARGET ${target_library} PROPERTY CXX_STANDARD ${LIBMESH_CXX_VERSION})
    TARGET_COMPILE_FEATURES(${target_library} PUBLIC "cxx_std_${LIBMESH_CXX_VERSION}")
  ELSE()
    SET_PROPERTY(TARGET ${target_library} PROPERTY CXX_STANDARD 11)
    TARGET_COMPILE_FEATURES(${target_library} PUBLIC "cxx_std_11")
  ENDIF()

  # we and our users will use these MPI functions so make the interface public:
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC MPI::MPI_C)
  # If openMPI is built with the old C++ bindings enabled then we need to link
  # against them, since we don't consistently define OMPI_SKIP_MPICXX. We always
  # need this with libMesh since libMesh is underlinked.
  IF(${MPI_MPICXX_FOUND})
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC MPI::MPI_CXX)
  ENDIF()
  # Silo is underlinked and depends on HDF5, so do it first:
  IF(${IBAMR_HAVE_SILO})
    TARGET_LINK_LIBRARIES(${target_library} PRIVATE "${SILO_LIBRARIES}")
    TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${SILO_INCLUDE_DIRS}")
    # Silo is statically linked so we need to export its zlib dependency
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC ZLIB::ZLIB)
  ENDIF()
  # SAMRAI is also underlinked and depends on HDF5:
  TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC ${SAMRAI_INCLUDE_DIRS})
  FOREACH(_lib ${SAMRAI${_d}d_LIBRARIES})
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC ${_lib})
  ENDFOREACH()
  # HDF5:
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC "${HDF5_LIBRARIES}")
  TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${HDF5_INCLUDE_DIRS}")
  # Hypre:
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC "${HYPRE_LIBRARIES}")
  TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${HYPRE_INCLUDE_DIRS}")
  # Boost (we only need headers):
  IF(${IBAMR_USE_BUNDLED_BOOST})
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC BUNDLED_BOOST)
  ELSE()
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC Boost::headers)
  ENDIF()
  # Eigen3:
  IF(${IBAMR_USE_BUNDLED_EIGEN3})
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC BUNDLED_EIGEN3)
  ELSE()
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC Eigen3::Eigen)
  ENDIF()
  # muParser:
  IF(${IBAMR_USE_BUNDLED_MUPARSER})
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC BUNDLED_MUPARSER)
  ELSE()
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC "${MUPARSER_LIBRARIES}")
    TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${MUPARSER_INCLUDE_DIRS}")
  ENDIF()
  # PETSc:
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC "${PETSC_LIBRARIES}")
  TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${PETSC_INCLUDE_DIRS}")
  # libMesh:
  IF(${IBAMR_HAVE_LIBMESH})
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC "${LIBMESH_LIBRARIES}")
    TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${LIBMESH_INCLUDE_DIRS}")
  ENDIF()
ENDFUNCTION()

MESSAGE(STATUS "")
MESSAGE(STATUS "IBAMR dependencies have been successfully set up.")

#
# IBTK and IBAMR are compiled in 2D and 3D, but the headers are only installed
# once - hence set up the headers as separate targets
#
ADD_LIBRARY(IBTKHeaders INTERFACE)
TARGET_INCLUDE_DIRECTORIES(
  IBTKHeaders
  INTERFACE
  $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/ibtk/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/>)

ADD_LIBRARY(IBAMRHeaders INTERFACE)
TARGET_INCLUDE_DIRECTORIES(
  IBAMRHeaders
  INTERFACE
  $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/>)

INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ibamr
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/ibtk/include/ibtk
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

INSTALL(TARGETS IBTKHeaders EXPORT IBAMRTargets)
INSTALL(TARGETS IBAMRHeaders EXPORT IBAMRTargets)

# ---------------------------------------------------------------------------- #
#                               5: Define targets                              #
# ---------------------------------------------------------------------------- #

#
# Set up specific targets for executables and libraries that want to link to
# IBAMR. As noted above, we have to do this since SAMRAI might not be compiled
# with -fPIC: i.e., if we added SAMRAI as a link dependency to libIBTK and
# libIBTK as a link dependency to libFoo, then the linker would encounter errors
# as a result. With executables we can use the full link interface.
#
SET(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/ibamr/)

# Set up actual files containing the export target information:
INSTALL(EXPORT "IBAMRTargets" FILE "IBAMRTargets.cmake"
  NAMESPACE IBAMR::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ibamr/)

# Do the same for in-tree:
EXPORT(EXPORT "IBAMRTargets"
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/IBAMRTargets.cmake"
  NAMESPACE IBAMR::)

WRITE_BASIC_PACKAGE_VERSION_FILE(
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/IBAMRConfigVersion.cmake"
  VERSION ${IBTK_VERSION}
  COMPATIBILITY AnyNewerVersion)

CONFIGURE_PACKAGE_CONFIG_FILE(
  ${CMAKE_SOURCE_DIR}/cmake/IBAMRConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/IBAMRConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CONFIGDIR})

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/IBAMRConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/IBAMRConfigVersion.cmake
  DESTINATION ${INSTALL_CONFIGDIR})

# proceed to compilation units:
ADD_SUBDIRECTORY(ibtk)
ADD_SUBDIRECTORY(src)

if (IBAMR_ENABLE_TESTING)
  enable_testing()
  ADD_SUBDIRECTORY(tests)
endif()
ADD_SUBDIRECTORY(examples)
