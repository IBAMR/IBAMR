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
    "IB-ex2-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IB/explicit/ex2"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-IB
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d.cylinder input2d.cylinder_stabilized
    cylinder2d_128.spring cylinder2d_128.vertex
    cylinder2d_256.spring cylinder2d_256.vertex
    cylinder2d_512.spring cylinder2d_512.vertex
    cylinder2d_1024.spring cylinder2d_1024.vertex
    )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "IB-ex2-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IB/explicit/ex2"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-IB
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d.cylinder input3d.cylinder.test
    input3d.sphere
    cylinder3d_128.vertex
    cylinder3d_256.vertex
    cylinder3d_512.vertex
    sphere3d_32.spring sphere3d_32.vertex
    sphere3d_64.spring sphere3d_64.vertex
    sphere3d_128.spring sphere3d_128.vertex
    sphere3d_256.spring sphere3d_256.vertex)
