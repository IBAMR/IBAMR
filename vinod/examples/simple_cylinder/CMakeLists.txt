# CMakeLists.txt for simple_cylinder example

cmake_minimum_required(VERSION 3.15.0)

project(simple_cylinder CXX)

# Find IBAMR
find_package(IBAMR REQUIRED)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add executable
add_executable(main2d main.cpp)

# Link IBAMR libraries
target_link_libraries(main2d PRIVATE IBAMR::IBAMR2d)

# Set output directory
set_target_properties(main2d PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Installation (optional)
install(TARGETS main2d DESTINATION bin)
