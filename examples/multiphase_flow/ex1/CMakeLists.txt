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
    "multiphase_flow-ex1-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/multiphase_flow/ex1"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-multiphase_flow
  SOURCES
    GravityForcing.cpp LSLocateCircularInterface.cpp LSLocateGasInterface.cpp
    RigidBodyKinematics.cpp SetFluidGasSolidDensity.cpp
    SetFluidGasSolidViscosity.cpp SetLSProperties.cpp TagLSRefinementCells.cpp
    example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d.cyl input2d_free_fall.cyl cylinder2d.vertex
    cylinder2d_free_fall.vertex
  )
