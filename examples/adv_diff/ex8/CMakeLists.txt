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
    "adv_diff-ex8-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/adv_diff/ex8"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-adv_diff
  SOURCES
    example.cpp LevelSetInitialCondition.cpp SetFluidProperties.cpp BoussinesqForcing.cpp TagLSRefinementCells.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d input2d.amr
    )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "adv_diff-ex8-write_ascii"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/adv_diff/ex8"
  OUTPUT_NAME
    write_ascii
  EXAMPLE_GROUP
    examples-adv_diff
  LINK_TARGETS
    IBAMR2d
  SOURCES
    write_ascii.cpp)
