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
    "multiphase_flow-ex12-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/multiphase_flow/ex12"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-multiphase_flow
  SOURCES
    GravityForcing.cpp LSLocateGasInterface.cpp LSLocateStructureInterface.cpp
    LevelSetSolidInitialCondition.cpp RigidBodyKinematics.cpp SetFluidGasSolidDensity.cpp
    SetFluidGasSolidViscosity.cpp SetLSProperties.cpp TagLSRefinementCells.cpp
    example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d
  EXTRA_FILES
    Wedge2d.m wedge2d.vertex
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "multiphase_flow-ex12-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/multiphase_flow/ex12"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-multiphase_flow
  SOURCES
    GravityForcing.cpp LSLocateGasInterface.cpp LSLocateStructureInterface.cpp
    LevelSetSolidInitialCondition.cpp RigidBodyKinematics.cpp SetFluidGasSolidDensity.cpp
    SetFluidGasSolidViscosity.cpp SetLSProperties.cpp TagLSRefinementCells.cpp
    example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d
  EXTRA_FILES
    Wedge3d.m
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "multiphase_flow-ex12-wedgegen"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/multiphase_flow/ex12"
  OUTPUT_NAME
    wedgeGen
  EXAMPLE_GROUP
    examples-multiphase_flow
  SOURCES
    wedgeGen.cpp
  )
