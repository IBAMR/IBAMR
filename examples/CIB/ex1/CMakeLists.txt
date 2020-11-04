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
    "CIB-ex1"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/CIB/ex1"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-CIB
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d shell_3d_in.vertex shell_3d_out.vertex petsc_options.dat
  )
