## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2021 by the IBAMR developers
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
    "IIM-ex5-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IIM/ex5"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-IIM
  SOURCES
    example.cpp
  REQUIRES
    IBAMR_HAVE_LIBMESH
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d.steel
    input3d.glass
    )
