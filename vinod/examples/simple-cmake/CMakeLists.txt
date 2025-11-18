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

PROJECT(simple-cmake-example)
CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

# Find IBAMR package
FIND_PACKAGE(IBAMR 0.10.0 REQUIRED)

# Source files for the example
SET(SOURCE_FILES
    main.cpp
)

# Create executable
ADD_EXECUTABLE(simple-cmake-example ${SOURCE_FILES})

# Link against IBAMR
TARGET_LINK_LIBRARIES(simple-cmake-example IBAMR::IBAMR2d)

# Set compiler flags from IBAMR
SET(CMAKE_CXX_FLAGS ${IBAMR_CXX_FLAGS})

# Set C++ standard
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
