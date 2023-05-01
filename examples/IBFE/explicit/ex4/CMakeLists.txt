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
    "IBFE-ex4-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IBFE/explicit/ex4"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-IBFE
  SOURCES
    example.cpp
  REQUIRES
    IBAMR_HAVE_LIBMESH
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d input2d.markers
    )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "IBFE-ex4-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IBFE/explicit/ex4"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-IBFE
  SOURCES
    example.cpp
  REQUIRES
    IBAMR_HAVE_LIBMESH
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d
    )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "IBFE-ex4-convergence-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IBFE/explicit/ex4"
  OUTPUT_NAME
    convergence_tester2d
  EXAMPLE_GROUP
    examples-IBFE
  SOURCES
    convergence_tester.cpp
  REQUIRES
    IBAMR_HAVE_LIBMESH
  LINK_TARGETS
    IBAMR2d
    )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "IBFE-ex4-convergence-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IBFE/explicit/ex4"
  OUTPUT_NAME
    convergence_tester3d
  EXAMPLE_GROUP
    examples-IBFE
  SOURCES
    convergence_tester.cpp
  REQUIRES
    IBAMR_HAVE_LIBMESH
  LINK_TARGETS
    IBAMR3d
    )
