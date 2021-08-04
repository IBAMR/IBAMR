## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2021 by the IBAMR developers
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
    "complex_fluids-ex2"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/complex_fluids/ex2"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-complex_fluids
  SOURCES
    example.cpp InterpolationUtilities.cpp
  REQUIRES
    IBAMR_HAVE_LIBMESH
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d
  EXTRA_FILES
    plot_sxx.m
    )

IBAMR_ADD_EXAMPLE(
  TARGET_NAME
    "complex_fluids-ex2-convergence"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/complex_fluids/ex2"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-complex_fluids
  SOURCES
    convergence_tester.cpp
  REQUIRES
    IBAMR_HAVE_LIBMESH
  LINK_TARGETS
    IBAMR2d
    )
