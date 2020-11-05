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
    "navier_stokes-ex2"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/navier_stokes/ex2"
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
  EXTRA_FILES
    L1norm.m L1normx.m L1normy.m
    L2norm.m L2normx.m L2normy.m
    analyzer.m
  )
