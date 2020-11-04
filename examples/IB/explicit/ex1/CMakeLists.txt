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
    "IB-ex1"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/IB/explicit/ex1"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples-IB
  SOURCES
    example.cpp
  LINK_TARGETS
    IBAMR2d
  INPUT_FILES
    input2d input2d.shell input2d.shell_circum_fibers
    curve2d_64.spring curve2d_64.vertex curve2d_128.spring curve2d_128.vertex
    curve2d_256.spring curve2d_256.vertex curve2d_512.spring curve2d_512.vertex
    shell2d_64.spring shell2d_64.vertex shell2d_128.spring shell2d_128.vertex
    shell2d_256.spring shell2d_256.vertex
    shell2d_radial_64.spring shell2d_radial_64.vertex
    shell2d_radial_128.spring shell2d_radial_128.vertex
    shell2d_radial_256.spring shell2d_radial_256.vertex
    sphere3d.spring sphere3d.vertex
  EXTRA_FILES
    generate_curve2d.m generate_shell2d.m
    )
