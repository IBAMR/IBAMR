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
    "IB-ex3"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IB/explicit/ex3"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-IB
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d
    fila_256.beam fila_256.spring fila_256.target fila_256.vertex
    fila_512.beam fila_512.spring fila_512.target fila_512.vertex
    fila_1024.beam fila_1024.spring fila_1024.target fila_1024.vertex
    )
