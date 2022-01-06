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
SET(test_exclusions)

# Platform specific exclusions:
IF ("$ENV{CMAKE_CONFIGURATION}" MATCHES "fedora")
  LIST(APPEND test_exclusions
    # Comment on why the test fails.
    # Regular expression matching test: "^RenderMesh$"
  )
ENDIF ()

STRING(REPLACE ";" "|" test_exclusions "${test_exclusions}")
IF (test_exclusions)
  SET(test_exclusions "(${test_exclusions})")
ENDIF ()
