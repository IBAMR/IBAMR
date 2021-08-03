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

cmake_minimum_required(VERSION 3.8)

include("${CMAKE_CURRENT_LIST_DIR}/github_ci.cmake")

# Read the files from the build directory.
ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
# Uncomment the line below and all ctest_submit_multi lines
# if CTestConfig.cmake defines `drop_sites`.
# include("${CMAKE_CURRENT_LIST_DIR}/ctest_submit_multi.cmake")

# Pick up from where the configure left off.
ctest_start(APPEND)

include(ProcessorCount)
ProcessorCount(nproc)

include("${CMAKE_CURRENT_LIST_DIR}/ctest_exclusions.cmake")
ctest_test(
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
ctest_submit(PARTS Test)

if (test_result)
  message(FATAL_ERROR
    "Failed to test")
endif ()
