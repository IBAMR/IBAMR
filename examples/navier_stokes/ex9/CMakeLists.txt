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
    "navier_stokes-ex9-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/navier_stokes/ex9"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-navier_stokes
  SOURCES
    ../ex0/example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    README.md
    input3d.taylor_green_les
    input3d.decaying_hit
    input3d.channel_flow_les
    input3d.channel_flow_retau180_spinup
    input3d.channel_flow_retau180_average
  )
