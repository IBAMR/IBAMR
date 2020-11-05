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
    "ConstraintIB-filling_sphere"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/ConstraintIB/falling_sphere"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-ConstraintIB
  SOURCES
    CartGridBodyForce.cpp ForceProjector.cpp RigidBodyKinematics.cpp
    example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d sphere3d.vertex
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "ConstraintIB-spheregen"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/ConstraintIB/falling_sphere"
  OUTPUT_NAME
    sphereGen3d
  EXAMPLE_GROUP
    examples-ConstraintIB
  SOURCES
    sphereGen3d.cpp
  )
