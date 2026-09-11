# This target is deliberately separate from ordinary library tests.
ADD_CUSTOM_TARGET(tests-example_helpers)
FOREACH(_dim 2 3)
  ADD_EXECUTABLE(tests-example_helpers_${_dim}d EXCLUDE_FROM_ALL helper_contracts.cpp)
  SET_TARGET_PROPERTIES(tests-example_helpers_${_dim}d PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/example_helpers"
    OUTPUT_NAME helper_contracts_${_dim}d)
  TARGET_LINK_LIBRARIES(tests-example_helpers_${_dim}d PRIVATE phase_change_helpers_${_dim}d)
  ADD_DEPENDENCIES(tests-example_helpers tests-example_helpers_${_dim}d)
ENDFOREACH()
ADD_DEPENDENCIES(tests-example_helpers
  multiphase_flow-ex0-2d
  multiphase_flow-ex0-3d
  multiphase_flow-ex3-2d
  multiphase_flow-ex3-3d
  multiphase_flow-ex4
  multiphase_flow-ex14-2d
  multiphase_flow-ex14-3d
  phase_change-ex0-2d
  phase_change-ex1-2d
  phase_change-ex2-2d
  phase_change-ex3-2d
  phase_change-ex4-2d
  phase_change-ex5-2d
  phase_change-ex6-2d
  phase_change-ex7-2d
  phase_change-ex7-3d
  phase_change-ex8-2d
  phase_change-ex8-3d
  phase_change-ex9-2d
  phase_change-ex10-2d)
ADD_CUSTOM_COMMAND(TARGET tests-example_helpers POST_BUILD
  COMMAND bash ${CMAKE_SOURCE_DIR}/tests/link-test-files.sh
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/tests/example_helpers
  VERBATIM)
