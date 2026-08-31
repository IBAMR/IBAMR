## ---------------------------------------------------------------------
##
## Copyright (c) 2026 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/../../cmake/FixPETScHomebrewPaths.cmake")
IF(NOT DEFINED FIXTURE_ROOT)
  MESSAGE(FATAL_ERROR "Pass a FIXTURE_ROOT in the build or temporary directory")
ENDIF()
FILE(MAKE_DIRECTORY "${FIXTURE_ROOT}")
# Exercise Apple path handling on every test host, then check the non-Apple no-op.
SET(APPLE TRUE)
SET(_checks 0)

FUNCTION(CHECK _name _input _expected)
  SET(_flags "${_input}")
  IBAMR_FIX_PETSC_HOMEBREW_PATHS(_flags)
  IF(NOT "${_flags}" STREQUAL "${_expected}")
    MESSAGE(FATAL_ERROR "${_name}\nexpected: ${_expected}\nactual:   ${_flags}")
  ENDIF()
  IBAMR_FIX_PETSC_HOMEBREW_PATHS(_flags)
  IF(NOT "${_flags}" STREQUAL "${_expected}")
    MESSAGE(FATAL_ERROR "${_name}: not idempotent")
  ENDIF()
  MATH(EXPR _checks "${_checks} + 1")
  SET(_checks "${_checks}" PARENT_SCOPE)
  MESSAGE(STATUS "PASS ${_name}")
ENDFUNCTION()

FOREACH(_layout arm64 intel custom)
  IF(_layout STREQUAL "arm64")
    SET(_prefix "${FIXTURE_ROOT}/arm64/opt/homebrew")
  ELSEIF(_layout STREQUAL "intel")
    SET(_prefix "${FIXTURE_ROOT}/intel/usr/local")
  ELSE()
    SET(_prefix "${FIXTURE_ROOT}/custom/brew")
  ENDIF()
  FOREACH(_formula open-mpi gcc)
    FILE(MAKE_DIRECTORY "${_prefix}/Cellar/${_formula}/current/include"
                        "${_prefix}/Cellar/${_formula}/current/lib/gcc/current"
                        "${_prefix}/Cellar/${_formula}/old-but-present/lib"
                        "${_prefix}/opt")
    FILE(CREATE_LINK "${_prefix}/Cellar/${_formula}/current"
                     "${_prefix}/opt/${_formula}" SYMBOLIC RESULT _error)
    IF(NOT _error STREQUAL "0")
      MESSAGE(FATAL_ERROR "Cannot create fixture opt link: ${_error}")
    ENDIF()
  ENDFOREACH()
  SET(_old_mpi "${_prefix}/Cellar/open-mpi/missing")
  SET(_new_mpi "${_prefix}/opt/open-mpi")
  SET(_old_gcc "${_prefix}/Cellar/gcc/missing_1")
  SET(_new_gcc "${_prefix}/opt/gcc")

  CHECK("${_layout}: MPI include" "-I${_old_mpi}/include" "-I${_new_mpi}/include")
  CHECK("${_layout}: GCC library directory" "-L${_old_gcc}/lib/gcc/current" "-L${_new_gcc}/lib/gcc/current")
  CHECK("${_layout}: rpath" "-Wl,-rpath,${_old_gcc}/lib/gcc/current" "-Wl,-rpath,${_new_gcc}/lib/gcc/current")
  CHECK("${_layout}: two rpaths" "-Wl,-rpath,${_old_mpi}/lib,-rpath,${_old_gcc}/lib" "-Wl,-rpath,${_new_mpi}/lib,-rpath,${_new_gcc}/lib")
  CHECK("${_layout}: rpath equals" "-Wl,-rpath=${_old_mpi}/lib" "-Wl,-rpath=${_new_mpi}/lib")
  CHECK("${_layout}: separate option" "-isystem;${_old_mpi}/include" "-isystem;${_new_mpi}/include")
  CHECK("${_layout}: existing version" "-L${_prefix}/Cellar/gcc/old-but-present/lib" "-L${_prefix}/Cellar/gcc/old-but-present/lib")
  CHECK("${_layout}: opt unchanged" "-I${_new_mpi}/include" "-I${_new_mpi}/include")
  CHECK("${_layout}: missing suffix" "-L${_old_gcc}/lib/missing" "-L${_old_gcc}/lib/missing")
  CHECK("${_layout}: nonexistent library" "${_old_mpi}/lib/no-such-library.dylib" "${_old_mpi}/lib/no-such-library.dylib")
  CHECK("${_layout}: prefix collision" "-Wl,-rpath,${_old_gcc}/lib,-rpath,${_old_gcc}/lib_missing" "-Wl,-rpath,${_new_gcc}/lib,-rpath,${_old_gcc}/lib_missing")
  CHECK("${_layout}: duplicates and order" "-L${_old_mpi}/lib;-lmpi;-L${_old_mpi}/lib;-lm" "-L${_new_mpi}/lib;-lmpi;-L${_new_mpi}/lib;-lm")
  CHECK("${_layout}: framework and unrelated flags" "-framework;Accelerate;-pthread;-Wl,-undefined,dynamic_lookup" "-framework;Accelerate;-pthread;-Wl,-undefined,dynamic_lookup")
  CHECK("${_layout}: non-Homebrew paths" "-I${FIXTURE_ROOT}/ordinary/include;-L${FIXTURE_ROOT}/ordinary/lib" "-I${FIXTURE_ROOT}/ordinary/include;-L${FIXTURE_ROOT}/ordinary/lib")
  FILE(MAKE_DIRECTORY "${_prefix}/opt/hdf5/include")
  CHECK("${_layout}: other formula" "-I${_prefix}/Cellar/hdf5/missing/include" "-I${_prefix}/Cellar/hdf5/missing/include")
  FILE(MAKE_DIRECTORY "${_prefix}/opt/gcc@14/include")
  CHECK("${_layout}: versioned formula" "-I${_prefix}/Cellar/gcc@14/missing/include" "-I${_prefix}/Cellar/gcc@14/missing/include")
  CHECK("${_layout}: missing opt" "-L${_prefix}/other/Cellar/gcc/missing/lib" "-L${_prefix}/other/Cellar/gcc/missing/lib")
  SET(APPLE FALSE)
  CHECK("${_layout}: non-Apple no-op" "-I${_old_mpi}/include;-L${_old_gcc}/lib" "-I${_old_mpi}/include;-L${_old_gcc}/lib")
  SET(APPLE TRUE)
ENDFOREACH()

CHECK("empty flags" "" "")
CHECK("ordinary comma flags" "-Wl,,--as-needed," "-Wl,,--as-needed,")
CHECK("leading comma" ",opaque,,flag," ",opaque,,flag,")
CHECK("empty list elements" ";-lm;;-pthread;" ";-lm;;-pthread;")
CHECK("duplicate paths inside one argument" "-Wl,-rpath,${_old_mpi}/lib,-rpath,${_old_mpi}/lib" "-Wl,-rpath,${_new_mpi}/lib,-rpath,${_new_mpi}/lib")
CHECK("preserve empty elements during correction" ";-L${_old_gcc}/lib;;" ";-L${_new_gcc}/lib;;")
CHECK("unrecognized option" "-DROOT=${_old_mpi}/include" "-DROOT=${_old_mpi}/include")
CHECK("library list stays ordered" "libA.a;-L${_old_mpi}/lib;libB.a;libA.a" "libA.a;-L${_new_mpi}/lib;libB.a;libA.a")
CHECK("formula root without suffix" "${_old_mpi}" "${_new_mpi}")
CHECK("leading empty before one include" ";-I${_old_mpi}/include" ";-I${_new_mpi}/include")
MESSAGE(STATUS "${_checks} checks passed, each also checked for idempotence")
