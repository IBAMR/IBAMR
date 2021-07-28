# Test launcher

if (NOT EXISTS ${EXEC})
  message(FATAL_ERROR "Could not find program: ${EXEC}")
endif ()

if (NOT EXISTS ${INPUT})
  message(FATAL_ERROR "Could not find input file: ${INPUT}")
endif ()

set(_command ${EXEC} ${INPUT})

if (NOT WORKDIR)
  set(WORKDIR ${CMAKE_CURRENT_SOURCE_DIR})
endif ()

if (LAUNCHER)
  message(STATUS "${LAUNCHER}")
  list(PREPEND _command ${LAUNCHER})
endif ()

execute_process(
  COMMAND ${_command}
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${WORKDIR}
  COMMAND_ECHO STDOUT
)

if (result)
  if (EXPECT_ERROR)
    return ()
  else ()
    message(FATAL_ERROR "Test failed...")
  endif ()
endif ()

if (RESTART_NUMBER GREATER 0)
  if (NOT EXISTS "${WORKDIR}/restart/")
    message(FATAL_ERROR "Could not find restart directory...")
  endif ()

  execute_process(
    COMMAND ${_command} "${WORKDIR}/restart/" ${RESTART_NUMBER}
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${WORKDIR}
    COMMAND_ECHO STDOUT
  )

  if (result)
    if (EXPECT_ERROR)
      return ()
    else ()
      message(FATAL_ERROR "Test failed (restart)")
    endif ()
  endif ()
endif ()
