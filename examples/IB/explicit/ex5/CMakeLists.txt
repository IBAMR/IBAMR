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
    "IB-ex5-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IB/explicit/ex5"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-IB
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d points2d.vertex
    )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "IB-ex5-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IB/explicit/ex5"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-IB
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d points3d.vertex
    )
