# Helpers are compiled only when an example or its dedicated verification target is requested.
FOREACH(_dim 2 3)
  ADD_LIBRARY(multiphase_flow_helpers_${_dim}d STATIC EXCLUDE_FROM_ALL
    LSLocateInterface.cpp CapillaryForces.cpp PointwiseLevelSet.cpp)
  TARGET_INCLUDE_DIRECTORIES(multiphase_flow_helpers_${_dim}d PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
  TARGET_LINK_LIBRARIES(multiphase_flow_helpers_${_dim}d PUBLIC IBAMR${_dim}d)
ENDFOREACH()
