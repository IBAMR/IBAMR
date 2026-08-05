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
    "IBPD-Turek_hron-2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IBPD/explicit/Turek_hron"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-IBPD
  SOURCES
    main.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "IBPD-Turek_hron-generate-sheet2d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IBPD/explicit/Turek_hron"
  OUTPUT_NAME
    sheet2d_generator
  EXAMPLE_GROUP
    examples-IBPD
  SOURCES
    generate_sheet2d.cpp
  )
