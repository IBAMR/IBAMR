## ---------------------------------------------------------------------
##
## Copyright (c) 2026 by the IBAMR developers
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
    "navier_stokes-ex8-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/navier_stokes/ex8"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-navier_stokes
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    README.md input2d.uniform input2d.rigid_rotation input2d.simple_shear input2d.trigonometric
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "navier_stokes-ex8-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/navier_stokes/ex8"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-navier_stokes
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    README.md input3d.uniform input3d.rigid_rotation input3d.simple_shear input3d.trigonometric
  )
