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
    "fe_mechanics-ex0-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/fe_mechanics/ex0"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-fe_mechanics
  SOURCES
    example.cpp
  REQUIRES
    IBAMR_HAVE_LIBMESH
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d
    )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "fe_mechanics-ex0-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/fe_mechanics/ex0"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-fe_mechanics
  SOURCES
    example.cpp
  REQUIRES
    IBAMR_HAVE_LIBMESH
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d
    )
