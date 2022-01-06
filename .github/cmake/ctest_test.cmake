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
CTEST_READ_CUSTOM_FILES("${CTEST_BINARY_DIRECTORY}")
# Uncomment the line below and all ctest_submit_multi lines
# if CTestConfig.cmake defines `drop_sites`.
# include("${CMAKE_CURRENT_LIST_DIR}/ctest_submit_multi.cmake")

# Pick up from where the configure left off.
CTEST_START(APPEND)

INCLUDE(ProcessorCount)
ProcessorCount(nproc)

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/ctest_exclusions.cmake")
CTEST_TEST(
  PARALLEL_LEVEL "${nproc}"
  RETURN_VALUE test_result
  EXCLUDE "${test_exclusions}")

# Below is an example location for packaged binaries.
# Adjust as needed should we want to upload binaries to cdash.
# file(GLOB packages
#   "${CTEST_BINARY_DIRECTORY}/PluginTests/*/build/plugin/build/*.tar.gz")
#
# if (NOT packages STREQUAL "")
#   ctest_upload(FILES ${packages})
# endif()

# ctest_submit_multi(PARTS Test)
CTEST_SUBMIT(PARTS Test)

IF (test_result)
  MESSAGE(FATAL_ERROR "Failed to test")
ENDIF ()
