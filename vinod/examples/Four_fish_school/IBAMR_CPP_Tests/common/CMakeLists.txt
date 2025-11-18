# CMakeLists.txt for common test utilities

cmake_minimum_required(VERSION 3.12)

# Create library for common utilities
add_library(test_utilities STATIC
    src/AnalyticalSolutions.cpp
    src/ErrorCalculator.cpp
    src/TestUtilities.cpp
)

target_include_directories(test_utilities
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Link with IBAMR and its dependencies
target_link_libraries(test_utilities
    PUBLIC
        IBAMR::IBAMR2d
        SAMRAI::SAMRAI
)

# Export the library
install(TARGETS test_utilities
    EXPORT TestUtilitiesTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/
    DESTINATION include
)
