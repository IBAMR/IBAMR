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

function (ctest_submit_multi)
  foreach (site IN LISTS drop_sites)
    foreach (variable IN ITEMS SITE METHOD LOCATION SITE_CDASH)
      if (NOT DEFINED "CTEST_DROP_${variable}_${site}")
        message(SEND_ERROR
          "The ${site} site is missing the CTEST_DROP_${variable} "
          "setting.")
      endif ()
      set("CTEST_DROP_${variable}" "${CTEST_DROP_${variable}_${site}}")
    endforeach ()

    ctest_submit(${ARGN})
  endforeach ()
endfunction ()
