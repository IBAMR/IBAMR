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

# List tests here whose results should be omitted.
set(test_exclusions
)

# Platform specific exclusions:
if ("$ENV{CMAKE_CONFIGURATION}" MATCHES "fedora")
  list(APPEND test_exclusions
    # Comment on why the test fails.
    # Regular expression matching test: "^RenderMesh$"
    # MPI Tests >2 MPI ranks oversubscribes on Github Actions
    "mpirun=[3-9]"
    # Additional tests disabled for CI
    "explicit_ex1_2d.mpirun=4.input"
    "explicit_ex2_3d.nodal_quadrature.input"
    "explicit_ex4_3d.mpirun=4.input"
    "explicit_ex5_3d.mpirun=2.input"
    "explicit_ex8_2d.input"
    "explicit_ex8_2d.scratch_hier.input"
    "free_falling_cyl_cib|cib_plate"
    "nwt_cylinder|rotating_barge|cib_double_shell"
  )
endif ()

string(REPLACE ";" "|" test_exclusions "${test_exclusions}")
if (test_exclusions)
  set(test_exclusions "(${test_exclusions})")
endif ()
