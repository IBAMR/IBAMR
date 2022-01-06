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

FUNCTION (ctest_submit_multi)
  FOREACH (site IN LISTS drop_sites)
    FOREACH (variable IN ITEMS SITE METHOD LOCATION SITE_CDASH)
      IF (NOT DEFINED "CTEST_DROP_${variable}_${site}")
        MESSAGE(SEND_ERROR
          "The ${site} site is missing the CTEST_DROP_${variable} "
          "setting.")
      ENDIF ()
      SET("CTEST_DROP_${variable}" "${CTEST_DROP_${variable}_${site}}")
    ENDFOREACH ()

    CTEST_SUBMIT(${ARGN})
  ENDFOREACH ()
ENDFUNCTION ()
