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
    "IBPD-ex2-3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IBPD/explicit/ex2"
  OUTPUT_NAME
    main3d
  EXAMPLE_GROUP
    examples-IBPD
  SOURCES
    main.cpp
  LINK_TARGETS
    IBAMR3d
  INPUT_FILES
    input3d.failure
  )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "IBPD-ex2-generate-failed-sheet3d"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IBPD/explicit/ex2"
  OUTPUT_NAME
    failed_sheet3d_generator
  EXAMPLE_GROUP
    examples-IBPD
  SOURCES
    generate_failed_sheet3d.cpp
  )
