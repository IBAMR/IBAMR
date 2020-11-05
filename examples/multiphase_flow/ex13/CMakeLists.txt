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
    "multiphase_flow-ex13"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/multiphase_flow/ex13"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-multiphase_flow
  SOURCES
    FlowGravityForcing.cpp GravityForcing.cpp LSLocateGasInterface.cpp
    LSLocateStructureInterface.cpp LevelSetGasInitialCondition.cpp
    LevelSetSolidInitialCondition.cpp RigidBodyKinematics.cpp
    SetFluidGasSolidDensity.cpp SetFluidGasSolidViscosity.cpp
    SetLSProperties.cpp TagLSRefinementCells.cpp
    example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d rectangle3d.vertex
  EXTRA_FILES
    Rectangle3d.m
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "multiphase_flow-ex13-rectanglegen"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/multiphase_flow/ex13"
  OUTPUT_NAME
    rectangleGen
  EXAMPLE_GROUP
    examples-multiphase_flow
  SOURCES
    rectangleGen.cpp
  )
