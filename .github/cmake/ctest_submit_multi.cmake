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
