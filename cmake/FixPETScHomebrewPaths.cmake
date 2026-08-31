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

# Homebrew relocates GCC/Open MPI paths in PETSc's pkg-config files, but
# petscvariables can retain obsolete Cellar paths. Repair only missing paths
# with existing opt counterparts; never modify the PETSc installation.
FUNCTION(IBAMR_FIX_PETSC_HOMEBREW_PATHS _variable)
  IF(NOT APPLE)
    RETURN()
  ENDIF()

  SET(_result "")
  SET(_argument_separator "")
  FOREACH(_argument IN LISTS ${_variable})
    # A linker argument can contain several comma-separated paths.
    STRING(REPLACE "," ";" _parts "${_argument}")
    SET(_fixed_argument "")
    SET(_part_separator "")
    FOREACH(_part IN LISTS _parts)
      IF(_part MATCHES "^(-I|-L|-rpath=)?(/[^ \t\r\n]+)(/Cellar/(gcc|open-mpi)/[^/ \t\r\n]+)(/[^ \t\r\n]*)?$")
        SET(_option "${CMAKE_MATCH_1}")
        SET(_path "${CMAKE_MATCH_2}${CMAKE_MATCH_3}${CMAKE_MATCH_5}")
        SET(_replacement "${CMAKE_MATCH_2}/opt/${CMAKE_MATCH_4}${CMAKE_MATCH_5}")
        IF(NOT EXISTS "${_path}" AND EXISTS "${_replacement}")
          SET(_part "${_option}${_replacement}")
        ENDIF()
      ENDIF()
      STRING(APPEND _fixed_argument "${_part_separator}${_part}")
      SET(_part_separator ",")
    ENDFOREACH()
    STRING(APPEND _result "${_argument_separator}${_fixed_argument}")
    SET(_argument_separator ";")
  ENDFOREACH()
  SET(${_variable} "${_result}" PARENT_SCOPE)
ENDFUNCTION()
