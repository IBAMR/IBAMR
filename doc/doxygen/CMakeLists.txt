## ---------------------------------------------------------------------
##
## Copyright (c) 2025 - 2025 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

FIND_PACKAGE(Doxygen)
IF(NOT DOXYGEN_FOUND)
  MESSAGE(FATAL_ERROR
    "Could not find doxygen which is required for building the documentation"
    )
ELSE()
  MESSAGE(STATUS "Using doxygen version ${DOXYGEN_VERSION}")
ENDIF()

#
# Use a custom layout file to disable the 'Related Pages' tab
#
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/DoxygenLayout.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml
  )

#
# Doxygen configuration file:
#
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/options.dox.in
  ${CMAKE_CURRENT_BINARY_DIR}/options.dox
  )

#
# Pre-process the Markdown files to get something Doxygen understands
#
FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated-headers")
SET(_generated_headers "")
MACRO(IBAMR_CONVERT_MD_TO_H _md_file _h_file _groupname _grouptitle _brieftext)
  SET(_outfile "${CMAKE_CURRENT_BINARY_DIR}/generated-headers/${_h_file}")
  ADD_CUSTOM_COMMAND(
    DEPENDS
      ${CMAKE_SOURCE_DIR}/doc/${_md_file}
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/filter.pl
    OUTPUT ${_outfile}
    COMMAND
      ${CMAKE_COMMAND} -E echo "/\*\* @defgroup ${_groupname} ${_grouptitle}"
        > ${_outfile} &&
      ${CMAKE_COMMAND} -E echo "@brief ${_brieftext}" >> ${_outfile} &&
      ${CMAKE_COMMAND} -E cat ${CMAKE_SOURCE_DIR}/doc/${_md_file}
        >> ${_outfile} &&
      ${CMAKE_COMMAND} -E echo "\*/" >> ${_outfile}
      VERBATIM
  )
  LIST(APPEND _generated_headers
    ${CMAKE_CURRENT_BINARY_DIR}/generated-headers/${_h_file})
ENDMACRO()

#
# Rather than list all headers use GLOB
#
FILE(GLOB _doxygen_dependencies
  "${SAMRAI_ROOT}/include/*.h"
  "${CMAKE_SOURCE_DIR}/include/ibamr/*.h"
  "${CMAKE_SOURCE_DIR}/ibtk/include/ibtk/*.h"
)

LIST(APPEND _doxygen_dependencies ${CMAKE_CURRENT_SOURCE_DIR}/main.md)
LIST(APPEND _doxygen_dependencies ${_generated_headers})

SET(_doxygen_input
  ${SAMRAI_ROOT}/include
  ${CMAKE_SOURCE_DIR}/include/ibamr
  ${CMAKE_SOURCE_DIR}/ibtk/include/ibtk

  ${CMAKE_BINARY_DIR}/ibtk/include/ibtk/config.h

  ${CMAKE_CURRENT_SOURCE_DIR}/main.md

  ${CMAKE_CURRENT_BINARY_DIR}/generated-headers
  )

MACRO(TO_STRING _variable)
  SET(${_variable} "")
  FOREACH(_var  ${ARGN})
    SET(${_variable} "${${_variable}} ${_var}")
  ENDFOREACH()
  STRING(STRIP "${${_variable}}" ${_variable})
ENDMACRO()

TO_STRING(_doxygen_input_string ${_doxygen_input})

FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/options.dox"
  "
  INPUT = ${_doxygen_input_string}
  "
  )

# If we use a reasonably modern doxygen version, make sure it is run in parallel
IF(NOT (${DOXYGEN_VERSION} VERSION_LESS 1.9))
  FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/options.dox" "
#---------------------------------------------------------------------------
# Configuration options related parallelization
#---------------------------------------------------------------------------
NUM_PROC_THREADS       = 0
DOT_NUM_THREADS        = 0"
    )
  MESSAGE(STATUS "Letting doxygen run with multiple threads")
ENDIF()

#
## Run Doxygen
##
ADD_CUSTOM_COMMAND(
  OUTPUT
    ${CMAKE_BINARY_DIR}/doxygen.log
  COMMAND ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    > ${CMAKE_BINARY_DIR}/doxygen.log 2>&1 # *pssst*
    || ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/doxygen.log ${CMAKE_BINARY_DIR}/doxygen.err
  COMMAND ${CMAKE_COMMAND} -E echo "-- Documentation is available at ${CMAKE_CURRENT_BINARY_DIR}/ibamr/index.html"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS
    ${_doxygen_dependencies}
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/filter.pl
  COMMENT "Generating documentation via doxygen."
  VERBATIM
  )
ADD_CUSTOM_TARGET(doxygen ALL
  DEPENDS ${CMAKE_BINARY_DIR}/doxygen.log
  )
ADD_DEPENDENCIES(documentation doxygen)
