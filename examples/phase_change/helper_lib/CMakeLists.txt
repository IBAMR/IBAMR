# Helpers are compiled only when an example or its dedicated verification target is requested.
FOREACH(_dim 2 3)
  ADD_LIBRARY(phase_change_helpers_${_dim}d STATIC EXCLUDE_FROM_ALL
    stefan.cpp thermocapillary.cpp)
  TARGET_INCLUDE_DIRECTORIES(phase_change_helpers_${_dim}d PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
  TARGET_LINK_LIBRARIES(phase_change_helpers_${_dim}d PUBLIC multiphase_flow_helpers_${_dim}d)
ENDFOREACH()
