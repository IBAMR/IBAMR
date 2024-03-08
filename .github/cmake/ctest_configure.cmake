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

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/github_ci.cmake")

SET(cmake_args
  -C "${CMAKE_CURRENT_LIST_DIR}/configure_$ENV{CMAKE_CONFIGURATION}.cmake"
)

# Create an entry in CDash.
CTEST_START(Experimental TRACK "${ctest_track}")

# Gather update information.
FIND_PACKAGE(Git)
SET(CTEST_UPDATE_VERSION_ONLY ON)
SET(CTEST_UPDATE_COMMAND "${GIT_EXECUTABLE}")
CTEST_UPDATE()

# Configure the project.
CTEST_CONFIGURE(
  OPTIONS "${cmake_args}"
  RETURN_VALUE configure_result
)

# Read the files from the build directory.
CTEST_READ_CUSTOM_FILES("${CTEST_BINARY_DIRECTORY}")

IF (configure_result)
  MESSAGE(FATAL_ERROR "Failed to configure ${configure_result}")
ENDIF ()
