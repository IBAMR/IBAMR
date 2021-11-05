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
    "ConstraintIB-eel3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/ConstraintIB/eel3d"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-ConstraintIB
  SOURCES
    IBEELKinematics3d.cpp example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d eel3d.vertex)

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "ConstraintIB-eel3d-generator"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/ConstraintIB/eel3d"
  OUTPUT_NAME
    eelgenerator3d
  EXAMPLE_GROUP
    examples-ConstraintIB
  LINK_TARGETS
    BOOST_INTERFACE
  SOURCES
    eelgenerator3d.cpp)
