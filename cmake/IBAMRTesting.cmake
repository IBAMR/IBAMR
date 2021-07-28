# Test configuration
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
        if (_x_value STREQUAL "true")
          set(expect_error TRUE)
        else ()
          set(expect_error FALSE)
        endif ()

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
    ""
    "LIBRARIES;SOURCES;CASES" ${ARGN})

  # Get the Test Group name
  string(REPLACE "${${PROJECT_NAME}_SOURCE_DIR}/tests/" "" dir_prefix ${CMAKE_CURRENT_SOURCE_DIR})
  string(REPLACE "/" "_" test_group ${dir_prefix})

  if (NOT test_CASES)
    list(APPEND test_CASES "<default>")
  endif ()

  if (NOT IBAMR_ENABLE_TESTING)
    return ()
  endif ()

  # Add test source if it exists (may not because of _2d/_3d suffix)
  if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/${name}.cpp")
    list(APPEND test_SOURCES
      "${CMAKE_CURRENT_LIST_DIR}/${name}.cpp")
  endif ()

  # If there aren't any sources, the this is not going to work...
  if (NOT test_SOURCES)
    message(FATAL_ERROR "No sources detected for test: ${test_group}/${name}")
  endif ()

  # Find all of the test inputs
  set(input_files)
  foreach (case ${test_CASES})
    _ibamr_parse_test(${case})
    # Skip adding tests that cannot be run
    if (n_mpirun GREATER IBAMR_MAX_PARALLEL AND NOT IBAMR_ALLOW_OVERSUBSCRIBED)
      continue ()
    endif ()

    # Get the filename
    if (case STREQUAL "<default>")
      set(case_filename "${CMAKE_CURRENT_LIST_DIR}/${name}")
    else ()
      set(case_filename "${CMAKE_CURRENT_LIST_DIR}/${name}.${case}")
    endif ()

    # Make sure the input file exists!
    if (NOT EXISTS "${case_filename}.input")
      message(FATAL_ERROR
        "Cannot find input for test case: ${test_group}-${name}-${case}"
        "   ${case_filename}.input")
    endif ()

    # Make sure the output file exists.
    if (NOT expect_error AND NOT EXISTS "${case_filename}.output")
      message(WARNING
        "Cannot find output for test case: ${test_group}-${name}-${case}"
        "   ${case_filename}.output")
    endif ()
    list(APPEND input_files "${case_filename}.input")
  endforeach ()

  if (NOT input_files)
    message(WARNING "No tests discovered for ${test_group}-${name}")
    return ()
  endif ()

  ###############################################################################################
  # Create build targets

  # Add a top level test build target
  if (NOT TARGET "ibamr-tests")
    add_custom_target(ibamr-tests ALL)
  endif ()

  # Add a directory build target and make the top level depend on it
  if (NOT TARGET "ibamr-tests-${test_group}")
    add_custom_target("ibamr-tests-${test_group}")
    add_dependencies(ibamr-tests "ibamr-tests-${test_group}")
  endif ()

  # Create the test executable
  set(_target "ibamr-test-${test_group}-${name}")
  add_executable(${_target} EXCLUDE_FROM_ALL "${test_SOURCES}")
  target_link_libraries(${_target} PRIVATE ${test_LIBRARIES})
  target_compile_definitions(${_target}
    PRIVATE
      SOURCE_DIR="${CMAKE_CURRENT_LIST_DIR}/")
  # Connect the target to the directory target
  add_dependencies("ibamr-tests-${test_group}" ${_target})

  # Add the tests for this target
  foreach (input ${input_files})
    get_filename_component(input_stripped ${input} NAME_WLE)
    get_filename_component(input_directory ${input} DIRECTORY)
    set(check_output_file "${input_directory}/${input_stripped}.output")

    _ibamr_parse_test(${input_stripped})

    set(case "${test_group}-${input_stripped}")
    set(test_workdir "${CMAKE_BINARY_DIR}/Testing/${case}")
    set(test_labels ${test_group})
    set(test_list)


    set(_launcher)
    set(_extra_args)
    set(_test_properties)
    set(_test_cpu_load 1)

    # Add MPI flags
    if (n_mpirun GREATER 1)
      set(_test_cpu_load ${n_mpirun})
      # Oversubscribe
      set(_mpi_flags ${MPIEXEC_PREFLAGS} --bind-to none)
      if (n_mpirun GREATER IBAMR_MAX_PARALLEL)
        if (IBAMR_ALLOW_OVERSUBSCRIBED)
          message(WARNING "Adding a test that is oversubscribed")
          list(APPEND _mpi_flags --oversubscribe)
          set(_test_cpu_load ${IBAMR_MAX_PARALLEL})
          list(APPEND _test_properties RUN_SERIAL TRUE)
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
        cmake
          "-DWORKDIR:PATH=${test_workdir}"
          "-DLAUNCHER:STRING=${_launcher}"
          "-DEXPECT_ERROR:BOOL=${expect_error}"
          "-DRESTART_NUMBER:STRING=${n_restart}"
          "-DEXEC:PATH=$<TARGET_FILE:${_target}>"
          "-DINPUT:PATH=${input}"
          -P ${CMAKE_SOURCE_DIR}/cmake/IBAMRTestLauncher.cmake
      )
    list(APPEND test_list "${case}_run")

    # Set test properties for the test cases
    set_tests_properties(${test_list}
      PROPERTIES
        PROCESSORS ${_test_cpu_load}
        PROCESSOR_AFFINITY FALSE
        REQUIRED_FILES "${input}"
        ${_test_properties}
    )

    if (EXISTS ${check_output_file} AND NOT expect_error)
      set(test_output_file "${test_workdir}/output")
      add_test(
        NAME "${case}_check"
        COMMAND
          ${NUMDIFF_EXECUTABLE} -r 1e-6 -a 1e-10 -s "' \\t\\n=,:;<>[](){}^'"
          ${check_output_file} ${test_output_file})
      set_tests_properties("${case}_check"
        PROPERTIES
          DEPENDS "${case}_run"
          REQUIRED_FILES "${check_output_file}")
      list(APPEND test_list "${case}_check")
    endif ()

    # Set working directory and timeout for these tests
    set_tests_properties(${test_list}
      PROPERTIES
        WORKING_DIRECTORY "${test_workdir}"
        TIMEOUT ${IBAMR_TEST_TIMEOUT}
        ENVIRONMENT "OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1"
    )

    # Setup and Cleanup test directories
    add_test(
      NAME "${case}_setup"
      COMMAND
        cmake -E make_directory ${test_workdir})
    add_test(
      NAME "${case}_cleanup"
      COMMAND
        cmake -E remove_directory ${test_workdir})

    # Configure FIXTURES
    set_tests_properties("${case}_setup"
      PROPERTIES
        FIXTURES_SETUP ${case})
    set_tests_properties(${test_list}
      PROPERTIES
        FIXTURES_REQUIRED ${case})
    if (IBAMR_TEST_CLEANUP_ON_FAILURE)
      # Cleanup, even if the tests fail
      set_tests_properties("${case}_cleanup"
        PROPERTIES
          FIXTURES_CLEANUP ${case})
    else ()
      # Cleanup, but only if the DEPENDS list finish
      # This will also allow a "run" step to be called explicity
      # and leave the results for inspection
      set_tests_properties("${case}_cleanup"
        PROPERTIES
          DEPENDS "${test_list}")
      # Add this to the list of tests that will get the
      # LABLES so dependent cleanup is also called when running
      # with labels
      list(APPEND test_list "${case}_cleanup")
    endif ()
    set_tests_properties(${test_list}
      PROPERTIES
        COST ${_test_cpu_load}
        LABELS "${test_labels}"
    )
  endforeach ()
endfunction ()

function (ibamr_add_2d_tests name)
  cmake_parse_arguments(test
    ""
    ""
    "LIBRARIES;SOURCES;CASES" ${ARGN})

  if (NOT test_CASES)
    list(APPEND test_CASES "<default>")
  endif ()

  if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/${name}.cpp")
    list(APPEND test_SOURCES
      "${CMAKE_CURRENT_LIST_DIR}/${name}.cpp")
  endif ()

  list(APPEND test_LIBRARIES
    IBAMR2d)

  ibamr_add_tests(
    ${name}_2d
    SOURCES ${test_SOURCES}
    LIBRARIES ${test_LIBRARIES}
    CASES ${test_CASES}
    )
endfunction ()

function (ibamr_add_3d_tests name)
  cmake_parse_arguments(test
    ""
    ""
    "LIBRARIES;SOURCES;CASES" ${ARGN})

  if (NOT test_CASES)
    list(APPEND test_CASES "<default>")
  endif ()

  if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/${name}.cpp")
    list(APPEND test_SOURCES
      "${CMAKE_CURRENT_LIST_DIR}/${name}.cpp")
  endif ()

  list(APPEND test_LIBRARIES
    IBAMR3d)

  ibamr_add_tests(
    ${name}_3d
    SOURCES "${test_SOURCES}"
    LIBRARIES "${test_LIBRARIES}"
    CASES ${test_CASES}
    )
endfunction ()
