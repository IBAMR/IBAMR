## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

ADD_CUSTOM_TARGET(examples)

SET(EXAMPLE_DIRECTORIES CIB ConstraintIB IB IBFE IBLevelSet IMP adv_diff
  advect complex_fluids fe_mechanics level_set multiphase_flow navier_stokes
  vc_navier_stokes wave_tank)

#
# Helper macro for setting up an example. Call as, for example,
#
# IBAMR_ADD_EXAMPLE(
#   TARGET_NAME
#     "IBFE-ex0"
#   OUTPUT_DIRECTORY
#     "${CMAKE_BINARY_DIR}/examples/IBFE/explicit/ex0"
#   OUTPUT_NAME
#     main2d
#   EXAMPLE_GROUP
#     examples-IBFE
#   SOURCES
#     "example.cpp"
#   REQUIRES
#     IBAMR_HAVE_LIBMESH
#   LINK_TARGETS
#     IBAMR2d
#   INPUT_FILES
#     input2d cylinder2d.vertex petsc_options.dat
#   )
#
# Here the target name is the name used by CMake to identify what needs to be
# compiled. The actual executable will be named (in this case) main2d and will
# be placed in the specified OUTPUT_DIRECTORY.
#
# TARGET_NAME, OUTPUT_DIRECTORY, OUTPUT_NAME, EXAMPLE_GROUP, and SOURCES are all
# mandatory arguments. The rest can be skipped.
#
MACRO(IBAMR_ADD_EXAMPLE)
  SET(OPTIONS)

  SET(ONE_VALUE_ARGS TARGET_NAME OUTPUT_DIRECTORY OUTPUT_NAME EXAMPLE_GROUP)
  SET(MULTI_VALUE_ARGS SOURCES REQUIRES LINK_TARGETS INPUT_FILES EXTRA_FILES)

  CMAKE_PARSE_ARGUMENTS(EX "${OPTIONS}" "${ONE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}"
    ${ARGN})

  SET(_failed_requirement)
  SET(_requirements_met TRUE)
  FOREACH(_requirement ${EX_REQUIRES})
    IF(NOT ${${_requirement}})
      SET(_requirements_met FALSE)
      SET(_failed_requirement ${_requirement})
    ENDIF()
  ENDFOREACH()

  IF(${_requirements_met})
    ADD_EXECUTABLE("${EX_TARGET_NAME}" EXCLUDE_FROM_ALL ${EX_SOURCES})
    SET_TARGET_PROPERTIES("${EX_TARGET_NAME}"
      PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY
      "${EX_OUTPUT_DIRECTORY}"
      OUTPUT_NAME
      "${EX_OUTPUT_NAME}")
    TARGET_INCLUDE_DIRECTORIES("${EX_TARGET_NAME}" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

    IF("${EX_LINK_TARGETS}" STREQUAL "")
      # If we have no target link libraries we need to set up the C++ standard ourselves
      TARGET_COMPILE_FEATURES("${EX_TARGET_NAME}" PUBLIC cxx_std_11)
    ELSE()
      TARGET_LINK_LIBRARIES("${EX_TARGET_NAME}" PRIVATE ${EX_LINK_TARGETS})
    ENDIF()
    ADD_DEPENDENCIES("${EX_EXAMPLE_GROUP}" "${EX_TARGET_NAME}")

    FOREACH(_input_file ${EX_INPUT_FILES})
      CONFIGURE_FILE("${_input_file}" "${EX_OUTPUT_DIRECTORY}" COPYONLY)
    ENDFOREACH()

    FOREACH(_extra_file ${EX_EXTRA_FILES})
      CONFIGURE_FILE("${_extra_file}" "${EX_OUTPUT_DIRECTORY}" COPYONLY)
    ENDFOREACH()
  ELSE()
    MESSAGE(STATUS "Example ${EX_TARGET_NAME} requires libMesh - skipping")
  ENDIF()
ENDMACRO()

#
# Add example subdirectories:
#
FOREACH(_dir ${EXAMPLE_DIRECTORIES})
  ADD_CUSTOM_TARGET("examples-${_dir}")
  ADD_DEPENDENCIES(examples "examples-${_dir}")
ENDFOREACH()

FOREACH(_dir ${EXAMPLE_DIRECTORIES})
  ADD_SUBDIRECTORY(${_dir})
ENDFOREACH()
