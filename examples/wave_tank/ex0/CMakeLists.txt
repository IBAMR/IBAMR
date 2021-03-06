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
    "wave_tank-ex0-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/wave_tank/ex0"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-wave_tank
  SOURCES
    GravityForcing.cpp LSLocateColumnInterface.cpp SetFluidProperties.cpp
    SetLSProperties.cpp TagLSRefinementCells.cpp example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d input2d_5thOrderStokes input2d_CaseA
    input2d_CaseA_IrregularWave input2d_CaseB input2d_CaseC
  EXTRA_FILES
    StokesWaveLength.m compareIrregularWave.m compareStokesFifthOrder.m
    compareStokesFirstOrder.m compareStokesSecondOrder.m
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "wave_tank-ex0-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/wave_tank/ex0"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-wave_tank
  SOURCES
    GravityForcing.cpp LSLocateColumnInterface.cpp SetFluidProperties.cpp
    SetLSProperties.cpp TagLSRefinementCells.cpp example.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d
  )
