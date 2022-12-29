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

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "navier_stokes-ex1-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/navier_stokes/ex1"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-navier_stokes
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "navier_stokes-ex1-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/navier_stokes/ex1"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-navier_stokes
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "navier_stokes-convergence-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/navier_stokes/ex1"
  OUTPUT_NAME
    convergence_tester2d
  EXAMPLE_GROUP
    examples-navier_stokes
  SOURCES
    convergence_tester.cpp
  LINK_TARGETS
    IBAMR2d
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "navier_stokes-convergence-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/navier_stokes/ex1"
  OUTPUT_NAME
    convergence_tester3d
  EXAMPLE_GROUP
    examples-navier_stokes
  SOURCES
    convergence_tester.cpp
  LINK_TARGETS
    IBAMR3d
  )
