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

# Read the files from the build directory.
ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
# Uncomment the line below and all ctest_submit_multi lines
# if CTestConfig.cmake defines `drop_sites`.
# include("${CMAKE_CURRENT_LIST_DIR}/ctest_submit_multi.cmake")

# Pick up from where the configure left off.
ctest_start(APPEND)

# Other generators automatically parallelize, but makefiles
# require an explicit thread-count.
IF(CTEST_CMAKE_GENERATOR STREQUAL "Unix Makefiles")
  INCLUDE(ProcessorCount)
  ProcessorCount(nproc)
  SET(CTEST_BUILD_FLAGS "-j${nproc}")
ENDIF()

ctest_build(
  NUMBER_WARNINGS num_warnings
  RETURN_VALUE build_result)
# ctest_submit_multi(PARTS Build)
ctest_submit(PARTS Build)

IF(build_result)
  MESSAGE(FATAL_ERROR "Failed to build")
ENDIF()

IF("$ENV{CTEST_NO_WARNINGS_ALLOWED}" AND num_warnings GREATER 0)
  MESSAGE(FATAL_ERROR "Found ${num_warnings} warnings (treating as fatal).")
ENDIF()
