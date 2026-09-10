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
    "navier_stokes-ex7-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/navier_stokes/ex7"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-navier_stokes
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    README.md input2d.cell input2d.node input2d.periodic.cell input2d.periodic.node
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "navier_stokes-ex7-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/navier_stokes/ex7"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-navier_stokes
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    README.md input3d.cell input3d.node input3d.periodic.cell input3d.periodic.node
  )
