## ---------------------------------------------------------------------
##
## Copyright (c) 2021 - 2021 by the IBAMR developers
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
    "IBLevelSet-ex2"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IBLevelSet/ex2"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-IBLevelSet
  SOURCES
    example.cpp
    FlowGravityForcing.cpp
    GravityForcing.cpp
    LSLocateGasInterface.cpp
    LevelSetInitialCondition.cpp
    SetFluidGasSolidDensity.cpp
    SetFluidGasSolidViscosity.cpp
    SetLSProperties.cpp
    TagLSRefinementCells.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d
    )
