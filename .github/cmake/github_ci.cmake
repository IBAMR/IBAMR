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

IF(NOT DEFINED "ENV{CI}")
  MESSAGE(FATAL_ERROR "This script assumes it is run inside a Github CI action.")
ENDIF()

# Set the build metadata.
SET(CTEST_BUILD_NAME "$ENV{GITHUB_REPOSITORY}-$ENV{CMAKE_CONFIGURATION}")
SET(CTEST_SITE "github-ci")

# Set up the source and build directories
SET(CTEST_SOURCE_DIRECTORY "$ENV{GITHUB_WORKSPACE}")
SET(CTEST_BINARY_DIRECTORY "/build/ibamr")

# On fedora, MPI must be loaded using environment modules:
IF("$ENV{CMAKE_CONFIGURATION}" MATCHES "fedora")
  FIND_PACKAGE(EnvModules REQUIRED)
  ENV_MODULE(purge)
  ENV_MODULE(load modules)
  ENV_MODULE(load mpi)
endif()

# Default to Release builds.
IF(NOT "$ENV{CMAKE_BUILD_TYPE}" STREQUAL "")
  SET(CTEST_BUILD_CONFIGURATION "$ENV{CMAKE_BUILD_TYPE}")
ENDIF()
IF(NOT CTEST_BUILD_CONFIGURATION)
  SET(CTEST_BUILD_CONFIGURATION "Release")
ENDIF()

# Default to using Ninja.
IF(NOT "$ENV{CMAKE_GENERATOR}" STREQUAL "")
  SET(CTEST_CMAKE_GENERATOR "$ENV{CMAKE_GENERATOR}")
ENDIF()
IF(NOT CTEST_CMAKE_GENERATOR)
  SET(CTEST_CMAKE_GENERATOR "Ninja")
ENDIF()

# Determine which track to submit to.
SET(ctest_track "experimental")
IF("$ENV{GITHUB_EVENT_NAME}" STREQUAL "pull_request")
  # Put pull request testing in the "merge requests" track
  SET(ctest_track "merge-requests")
ELSEIF("$ENV{GITHUB_EVENT_NAME}" STREQUAL "push")
  # Put all branch tests in a track named the same as the branch.
  # Currently, only master is tested, but eventually branch names
  # such as "release-0.9" may exist and should have their own tracks.
  STRING(REGEX REPLACE "refs\\/heads\\/" "" ctest_track "$ENV{GITHUB_REF}")
ENDIF()
