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
    "ConstraintIB-impulsively_started_cylinder"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/ConstraintIB/impulsively_started_cylinder"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-ConstraintIB
  SOURCES
    RigidBodyKinematics.cpp example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d cylinder2d.vertex)
