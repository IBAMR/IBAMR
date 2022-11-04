## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2022 by the IBAMR developers
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
    "DLM-ex0"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/DLM/ex0"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-DLM
  SOURCES
    example.cpp
  REQUIRES
    IBAMR_HAVE_LIBMESH
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d input2d.stiff input2d.very_stiff input2d.small_excess_mass input2d.smaller_excess_mass
    )
