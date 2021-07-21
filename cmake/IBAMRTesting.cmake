# Test configuration
set(IBAMR_TEST_DIR "${CMAKE_SOURCE_DIR}/tests/" CACHE PATH "Location of the test cases")
set(IBAMR_TEST_TIMEOUT 600 CACHE STRING "Test timeout")

include(ProcessorCount)
ProcessorCount(sys_nproc)
if (sys_nproc EQUAL 0)
  set(sys_nproc 1)
endif ()
set(IBAMR_MAX_PARALLEL ${sys_nproc} CACHE STRING "Maximum parallel level")
set(IBAMR_TEST_CLEANUP_ON_FAILURE TRUE CACHE BOOL "")

# Find numdiff
find_program(NUMDIFF_EXECUTABLE NAMES numdiff HINTS ${NUMDIFF_ROOT} PATH_SUFFIXES bin)
if ("${NUMDIFF_EXECUTABLE}" STREQUAL "NUMDIFF_EXECUTABLE-NOTFOUND")
  message(WARNING "\
The configuration script was not able to locate numdiff. If you want to run \
the test suite you will need to either edit attest.conf, specify the path to \
numdiff to attest, or rerun CMake with the argument NUMDIFF_ROOT specifying \
numdiff's root installation directory.")
  # clear the value so that attest.conf doesn't contain an invalid path
  set(NUMDIFF_EXECUTABLE "")
  set(IBAMR_ENABLE_TESTING FALSE)
endif ()

set(IBAMR_ALLOW_OVERSUBSCRIBED OFF "Allow tests to be added if they are oversubscribed")
mark_as_advanced(
  NUMDIFF_EXECUTABLE
  IBAMR_TEST_DIR
  IBAMR_TEST_TIMEOUT
  IBAMR_TEST_CLEANUP_ON_FAILURE
  IBAMR_MAX_PARALLEL
  IBAMR_ALLOW_OVERSUBSCRIBED)

function (_ibamr_parse_test name)
  set(n_mpirun 1)
  set(n_restart 0)
  set(expect_error FALSE)

  string(REPLACE "." ";" name_list "${name}")
  foreach (override ${name_list})
    string(REPLACE "=" ";" _x ${override})
    list(LENGTH _x _x_len)
    if (${_x_len} EQUAL 2)
      list(GET _x 0 _x_key)
      list(GET _x 1 _x_value)
      if (_x_key STREQUAL "mpirun")
        set(n_mpirun ${_x_value})
      elseif (_x_key STREQUAL "expect_error")
        string(TOUPPER "${_x_value}" value)
        set(expect_error ${value})
      elseif (_x_key STREQUAL "restart")
        set(n_restart ${_x_value})
      else (_x_key STREQUAL "restart")
        message(WARNING "Encountered unrecognized override ${_x_key}")
      endif ()
    endif ()
  endforeach ()

  set(n_mpirun ${n_mpirun} PARENT_SCOPE)
  set(n_restart ${n_restart} PARENT_SCOPE)
  set(expect_error ${expect_error} PARENT_SCOPE)
endfunction()

function (ibamr_add_tests name)
  cmake_parse_arguments(test
    ""
    "DIRECTORY"
    "LIBRARIES;SOURCES" ${ARGN})

  if (NOT test_DIRECTORY)
    message(FATAL_ERROR "Missing DIRECTORY for test")
  endif ()

  if (NOT IBAMR_ENABLE_TESTING)
    return ()
  endif ()

  # Add test source if it exists (may not because of _2d/_3d suffix)
  if (EXISTS "${IBAMR_TEST_DIR}/${test_DIRECTORY}/${name}.cpp")
    list(APPEND test_SOURCES "${IBAMR_TEST_DIR}/${test_DIRECTORY}/${name}.cpp")
  endif ()

  # If there aren't any sources, the this is not going to work...
  if (NOT test_SOURCES)
    message(FATAL_ERROR "No sources detected for test: ${test_DIRECTORY}/${name}")
  endif ()

  # Find all of the test inputs
  file(GLOB input_files CONFIGURE_DEPENDS "${IBAMR_TEST_DIR}/${test_DIRECTORY}/${name}.*.input")
  if (NOT input_files)
    return ()
  endif ()

  ###############################################################################################
  # Create build targets

  # Add a top level test build target
  if (NOT TARGET "ibamr-tests")
    add_custom_target(ibamr-tests ALL)
  endif ()

  # Add a directory build target and make the top level depend on it
  if (NOT TARGET "ibamr-tests-${test_DIRECTORY}")
    add_custom_target("ibamr-tests-${test_DIRECTORY}")
    add_dependencies(ibamr-tests "ibamr-tests-${test_DIRECTORY}")
  endif ()

  # Create the test executable
  set(_target "ibamr-test-${test_DIRECTORY}-${name}")
  add_executable(${_target} EXCLUDE_FROM_ALL "${test_SOURCES}")
  set_target_properties(${_target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      "${CMAKE_BINARY_DIR}/tests/${_dir}"
    OUTPUT_NAME
      ${name}
    )
  target_link_libraries(${_target} PRIVATE ${test_LIBRARIES})
  target_compile_definitions(${_target}
    PRIVATE
      SOURCE_DIR="${IBAMR_TEST_DIR}/${test_DIRECTORY}/")
  # Connect the target to the directory target
  add_dependencies("ibamr-tests-${test_DIRECTORY}" ${_target})

  # Add the tests for this target
  foreach (input ${input_files})
    get_filename_component(input_stripped ${input} NAME_WLE)
    get_filename_component(input_directory ${input} DIRECTORY)
    set(check_output_file "${input_directory}/${input_stripped}.output")

    _ibamr_parse_test(${input_stripped})

    set(case "${test_DIRECTORY}-${input_stripped}")
    set(test_workdir "${CMAKE_BINARY_DIR}/Testing/${case}")
    set(test_labels ${test_DIRECTORY})
    set(test_list)


    set(_launcher)
    set(_extra_args)
    set(load 1)
    # Add MPI flags
    if (n_mpirun GREATER 1)
      set(load ${n_mpirun})
      # Oversubscribe
      set(_mpi_flags ${MPIEXEC_PREFLAGS} --bind-to none)
      if (n_mpirun GREATER IBAMR_MAX_PARALLEL)
        if (IBAMR_ALLOW_OVERSUBSCRIBED)
          message(WARNING "Adding a test that is oversubscribed")
          list(APPEND _mpi_flags --oversubscribe)
          set(load ${IBAMR_MAX_PARALLEL})
        else ()
          # Don't add tests that would require oversubscription
          continue ()
        endif ()
      endif ()
      list(APPEND _launcher
        ${MPIEXEC} ${_mpi_flags} ${MPIEXEC_NUMPROC_FLAG} ${n_mpirun})
      list(APPEND test_labels "MPI")
    endif ()

    add_test(
      NAME "${case}_run"
      COMMAND
        ${_launcher} $<TARGET_FILE:${_target}> ${input} ${_extra_args}
      )
    list(APPEND test_list "${case}_run")

    # Add a restart test if it is active
    if (n_restart GREATER 0)
      list(PREPEND _extra_args "${test_workdir}/restart/${n_restart}")
      add_test(
        NAME "${case}_restart${n_restart}"
        COMMAND
        ${_launcher} $<TARGET_FILE:${_target}> ${input} ${_extra_args})
      set_tests_properties("${case}_restart${n_restart}"
        PROPERTIES
          DEPENDS "${test_list}")
      list(APPEND test_list "${case}_restart${n_restart}")
      list(APPEND test_labels "RESTART")
    endif ()

    # Set test properties for the test cases
    set_tests_properties(${test_list}
      PROPERTIES
        PARALLEL_LEVEL ${load}
        WILL_FAIL ${expect_error}
        REQUIRED_FILES "${input}"
    )

    if (EXISTS ${check_output_file})
      set(test_output_file "${test_workdir}/output")
      add_test(
        NAME "${case}_check"
        COMMAND
          ${NUMDIFF_EXECUTABLE} -r 1e-6 -a 1e-10 -s "' \\t\\n=,:;<>[](){}^'"
          ${check_output_file} ${test_output_file})
      set_tests_properties("${case}_check"
        PROPERTIES
          DEPENDS "${test_list}"
          REQUIRED_FILES "${check_output_file}")
      list(APPEND test_list "${case}_check")
    endif ()

    # Set working directory and timeout for these tests
    set_tests_properties(${test_list}
      PROPERTIES
        WORKING_DIRECTORY "${test_workdir}"
    )
    set_tests_properties(${test_list}
      PROPERTIES
        TIMEOUT ${IBAMR_TEST_TIMEOUT}
    )

    # Create working direcotry and cleanup test
    file(MAKE_DIRECTORY ${test_workdir})
    add_test(
      NAME "${case}_cleanup"
      COMMAND
        cmake -E remove_directory ${test_workdir})
    if (IBAMR_TEST_CLEANUP_ON_FAILURE)
      # Cleanup, even if the tests fail
      set_tests_properties("${case}_cleanup"
        PROPERTIES
          FIXTURES_CLEANUP ${case})
      set_tests_properties(${test_list}
        PROPERTIES
          FIXTURES_REQUIRED ${case})
    else ()
      # Cleanup, but only if the DEPENDS list finish
      set_tests_properties("${case}_cleanup"
        PROPERTIES
          DEPENDS ${test_list})
    endif ()
    list(APPEND test_list "${case}_cleanup")
    set_tests_properties(${test_list}
      PROPERTIES
        LABELS "${test_labels}"
    )
  endforeach ()
endfunction ()

function (ibamr_add_2d_tests name)
  cmake_parse_arguments(test
    ""
    "DIRECTORY"
    "LIBRARIES;SOURCES" ${ARGN})
  if (NOT test_DIRECTORY)
    message(FATAL_ERROR "Missing DIRECTORY for test")
  endif ()

  if (EXISTS "${IBAMR_TEST_DIR}/${test_DIRECTORY}/${name}.cpp")
    list(APPEND test_SOURCES
      "${IBAMR_TEST_DIR}/${test_DIRECTORY}/${name}.cpp")
  endif ()

  list(APPEND test_LIBRARIES
    IBAMR2d)

  ibamr_add_tests(
    ${name}_2d
    DIRECTORY ${test_DIRECTORY}
    SOURCES ${test_SOURCES}
    LIBRARIES ${test_LIBRARIES}
    )
endfunction ()

function (ibamr_add_3d_tests name)
  cmake_parse_arguments(test
    ""
    "DIRECTORY"
    "LIBRARIES;SOURCES" ${ARGN})

  if (NOT test_DIRECTORY)
    message(FATAL_ERROR "Missing DIRECTORY for test")
  endif ()

  if (EXISTS "${IBAMR_TEST_DIR}/${test_DIRECTORY}/${name}.cpp")
    list(APPEND test_SOURCES
      "${IBAMR_TEST_DIR}/${test_DIRECTORY}/${name}.cpp")
  endif ()

  list(APPEND test_LIBRARIES
    IBAMR3d)

  ibamr_add_tests(
    ${name}_3d
    DIRECTORY ${test_DIRECTORY}
    SOURCES "${test_SOURCES}"
    LIBRARIES "${test_LIBRARIES}"
    )
endfunction ()
