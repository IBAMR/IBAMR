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
    "CIB-ex2"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/CIB/ex2"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-CIB
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d.plate input2d.plate.amr petsc_options.dat plate2d.vertex
  EXTRA_FILES
    plate2d.m)
