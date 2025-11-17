# CMakeLists.txt for Vinod IBAMR Extensions Library
# ====================================================

cmake_minimum_required(VERSION 3.15.0)

project(VinodIBAMR
        VERSION 1.0.0
        DESCRIPTION "Custom IBAMR extensions and utilities"
        LANGUAGES CXX)

# =============================================================================
# Build Configuration
# =============================================================================

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# =============================================================================
# Find Dependencies
# =============================================================================

# Find IBAMR (required)
find_package(IBAMR 0.12 REQUIRED)

if(IBAMR_FOUND)
    message(STATUS "Found IBAMR: ${IBAMR_DIR}")
else()
    message(FATAL_ERROR "IBAMR not found. Please set IBAMR_DIR or CMAKE_PREFIX_PATH.")
endif()

# =============================================================================
# Library Target
# =============================================================================

# Source files
set(VINOD_SOURCES
    src/forces/MultiStructureForceTracker.cpp
    src/CustomForceFunction.cpp
)

# Header files
set(VINOD_HEADERS
    include/vinod/forces/MultiStructureForceTracker.h
    src/CustomForceFunction.h
)

# Create library
add_library(vinod_ibamr ${VINOD_SOURCES} ${VINOD_HEADERS})

# Set library properties
set_target_properties(vinod_ibamr PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER "${VINOD_HEADERS}"
)

# Include directories
target_include_directories(vinod_ibamr
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
)

# Link IBAMR
# NDIM should be set by parent project or passed via -DNDIM=2 or -DNDIM=3
if(NOT DEFINED NDIM)
    set(NDIM 2 CACHE STRING "Spatial dimension (2 or 3)")
    message(STATUS "NDIM not set, defaulting to NDIM=${NDIM}")
endif()

target_link_libraries(vinod_ibamr
    PUBLIC
        IBAMR::IBAMR${NDIM}d
)

# =============================================================================
# Examples
# =============================================================================

# Note: Examples have their own standalone CMakeLists.txt files
# Build them separately from their own directories
# Example:
#   cd examples/simple_cylinder && mkdir build && cd build
#   cmake .. -DIBAMR_DIR=/path/to/IBAMR/install
#   make

# =============================================================================
# Installation
# =============================================================================

include(GNUInstallDirs)

# Install library
install(TARGETS vinod_ibamr
    EXPORT VinodIBAMRTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/vinod
)

# Install headers (preserve directory structure)
install(DIRECTORY include/vinod
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h"
)

# Install CMake config files
install(EXPORT VinodIBAMRTargets
    FILE VinodIBAMRTargets.cmake
    NAMESPACE Vinod::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VinodIBAMR
)

# Create config file
include(CMakePackageConfigHelpers)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/VinodIBAMRConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/VinodIBAMRConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VinodIBAMR
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/VinodIBAMRConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/VinodIBAMRConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/VinodIBAMRConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VinodIBAMR
)

# =============================================================================
# Testing
# =============================================================================

option(BUILD_TESTING "Build test programs" OFF)

if(BUILD_TESTING)
    enable_testing()
    # Add test subdirectory when available
    # add_subdirectory(tests)
endif()

# =============================================================================
# Summary
# =============================================================================

message(STATUS "")
message(STATUS "====================================")
message(STATUS "Vinod IBAMR Extensions Configuration")
message(STATUS "====================================")
message(STATUS "Version: ${PROJECT_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Spatial dimension: ${NDIM}D")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Build testing: ${BUILD_TESTING}")
message(STATUS "====================================")
message(STATUS "")
