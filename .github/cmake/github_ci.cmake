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

if (NOT DEFINED "ENV{CI}")
  message(FATAL_ERROR "This script assumes it is run inside a Github CI action.")
endif()

# Set the build metadata.
set(CTEST_BUILD_NAME "$ENV{GITHUB_REPOSITORY}-$ENV{CMAKE_CONFIGURATION}")
set(CTEST_SITE "github-ci")

# Set up the source and build directories
set(CTEST_SOURCE_DIRECTORY "$ENV{GITHUB_WORKSPACE}")
set(CTEST_BINARY_DIRECTORY "/build/ibamr")

# On fedora, MPI must be loaded using environment modules:
if ("$ENV{CMAKE_CONFIGURATION}" MATCHES "fedora")
  find_package(EnvModules REQUIRED)
  env_module(purge)
  env_module(load modules)
  env_module(load mpi)
endif()

# Default to Release builds.
if (NOT "$ENV{CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CTEST_BUILD_CONFIGURATION "$ENV{CMAKE_BUILD_TYPE}")
endif ()
if (NOT CTEST_BUILD_CONFIGURATION)
  set(CTEST_BUILD_CONFIGURATION "Release")
endif ()

# Default to using Ninja.
if (NOT "$ENV{CMAKE_GENERATOR}" STREQUAL "")
  set(CTEST_CMAKE_GENERATOR "$ENV{CMAKE_GENERATOR}")
endif ()
if (NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Ninja")
endif ()

# Determine which track to submit to.
set(ctest_track "experimental")
if ("$ENV{GITHUB_EVENT_NAME}" STREQUAL "pull_request")
  # Put pull request testing in the "merge requests" track
  set(ctest_track "merge-requests")
elseif ("$ENV{GITHUB_EVENT_NAME}" STREQUAL "push")
  # Put all branch tests in a track named the same as the branch.
  # Currently, only master is tested, but eventually branch names
  # such as "release-0.9" may exist and should have their own tracks.
  string(REGEX_REPLACE "refs\\/heads\\/" "" ctest_track "$ENV{GITHUB_REF}")
endif ()
