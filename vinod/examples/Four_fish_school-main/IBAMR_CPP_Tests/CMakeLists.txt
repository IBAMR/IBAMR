# Master CMakeLists.txt for IBAMR C++ Test Suite
# Verification & Validation tests for scalar transport

cmake_minimum_required(VERSION 3.12)
project(IBAMR_Scalar_Transport_Tests VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(BUILD_ALL_TESTS "Build all tests (1-17)" ON)
option(BUILD_TIER1_ONLY "Build only Tier 1 tests (1-6)" OFF)
option(BUILD_TIER2_ONLY "Build only Tier 2 tests (7-10)" OFF)
option(BUILD_LITERATURE_ONLY "Build only Literature Validation tests (15-17)" OFF)

# Find IBAMR
if(NOT DEFINED IBAMR_ROOT)
    if(DEFINED ENV{IBAMR_ROOT})
        set(IBAMR_ROOT $ENV{IBAMR_ROOT})
    else()
        message(FATAL_ERROR "IBAMR_ROOT not set. Please set IBAMR_ROOT to IBAMR installation directory.")
    endif()
endif()

message(STATUS "Looking for IBAMR in: ${IBAMR_ROOT}")

# Add IBAMR to CMake prefix path
list(APPEND CMAKE_PREFIX_PATH ${IBAMR_ROOT})

# Find required packages
find_package(IBAMR REQUIRED)
find_package(SAMRAI REQUIRED)
find_package(PETSc REQUIRED)
find_package(MPI REQUIRED)

# Include directories
include_directories(
    ${IBAMR_INCLUDE_DIRS}
    ${SAMRAI_INCLUDE_DIRS}
    ${PETSC_INCLUDE_DIRS}
    ${MPI_INCLUDE_PATH}
)

# Add common utilities library
add_subdirectory(common)

# Macro to add a test executable
macro(add_ibamr_test TEST_NAME TEST_DIR)
    add_executable(${TEST_NAME} ${TEST_DIR}/main.cpp)

    target_link_libraries(${TEST_NAME}
        PRIVATE
            test_utilities
            IBAMR::IBAMR2d
            SAMRAI::SAMRAI
            PETSc::PETSc
            MPI::MPI_CXX
    )

    target_include_directories(${TEST_NAME}
        PRIVATE
            ${CMAKE_SOURCE_DIR}/common/include
    )

    # Install test executable
    install(TARGETS ${TEST_NAME}
        RUNTIME DESTINATION bin/${TEST_DIR}
    )

    # Install input file
    install(FILES ${TEST_DIR}/input2d
        DESTINATION bin/${TEST_DIR}
    )
endmacro()

# Define test lists
set(TIER1_TESTS
    test01_smoke
    test02_diffusion
    test03_advection
    test04_mms
    test05_discontinuous
    test06_mass_conservation
)

set(TIER2_TESTS
    test07_bcs
    test08_sphere_source
    test09_high_sc
    test10_moving_ib
)

set(TIER3_TESTS
    test11_amr
    test12_timestep
    test13_long_run
    test14_benchmarks
)

set(LITERATURE_VALIDATION_TESTS
    test15_rotating_cylinder
    test16_3d_sphere
    test17_pitch_plunge
)

# Build tests based on options
if(BUILD_TIER1_ONLY)
    set(TESTS_TO_BUILD ${TIER1_TESTS})
elseif(BUILD_TIER2_ONLY)
    set(TESTS_TO_BUILD ${TIER2_TESTS})
elseif(BUILD_LITERATURE_ONLY)
    set(TESTS_TO_BUILD ${LITERATURE_VALIDATION_TESTS})
elseif(BUILD_ALL_TESTS)
    set(TESTS_TO_BUILD ${TIER1_TESTS} ${TIER2_TESTS} ${TIER3_TESTS} ${LITERATURE_VALIDATION_TESTS})
else()
    set(TESTS_TO_BUILD ${TIER1_TESTS})
endif()

# Add test executables
foreach(TEST ${TESTS_TO_BUILD})
    # Map test name to directory
    if(TEST STREQUAL "test01_smoke")
        add_ibamr_test(${TEST} Test01_SmokeTest)
    elseif(TEST STREQUAL "test02_diffusion")
        add_ibamr_test(${TEST} Test02_Diffusion_Analytic)
    elseif(TEST STREQUAL "test03_advection")
        add_ibamr_test(${TEST} Test03_Advection_Analytic)
    elseif(TEST STREQUAL "test04_mms")
        add_ibamr_test(${TEST} Test04_MMS)
    elseif(TEST STREQUAL "test05_discontinuous")
        add_ibamr_test(${TEST} Test05_Discontinuous)
    elseif(TEST STREQUAL "test06_mass_conservation")
        add_ibamr_test(${TEST} Test06_MassConservation)
    elseif(TEST STREQUAL "test07_bcs")
        add_ibamr_test(${TEST} Test07_BCs)
    elseif(TEST STREQUAL "test08_sphere_source")
        add_ibamr_test(${TEST} Test08_SphereSource)
    elseif(TEST STREQUAL "test09_high_sc")
        add_ibamr_test(${TEST} Test09_HighSc)
    elseif(TEST STREQUAL "test10_moving_ib")
        add_ibamr_test(${TEST} Test10_MovingIB)
    elseif(TEST STREQUAL "test11_amr")
        add_ibamr_test(${TEST} Test11_AMR)
    elseif(TEST STREQUAL "test12_timestep")
        add_ibamr_test(${TEST} Test12_TimeStep)
    elseif(TEST STREQUAL "test13_long_run")
        add_ibamr_test(${TEST} Test13_LongRun)
    elseif(TEST STREQUAL "test14_benchmarks")
        add_ibamr_test(${TEST} Test14_Benchmarks)
    elseif(TEST STREQUAL "test15_rotating_cylinder")
        add_ibamr_test(${TEST} Test15_RotatingCylinder)
    elseif(TEST STREQUAL "test16_3d_sphere")
        add_ibamr_test(${TEST} Test16_3DSphere)
    elseif(TEST STREQUAL "test17_pitch_plunge")
        add_ibamr_test(${TEST} Test17_PitchPlunge)
    endif()
endforeach()

# Print build summary
message(STATUS "================================================")
message(STATUS "IBAMR Scalar Transport Test Suite Configuration")
message(STATUS "================================================")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "IBAMR Root: ${IBAMR_ROOT}")
message(STATUS "Tests to build:")
foreach(TEST ${TESTS_TO_BUILD})
    message(STATUS "  - ${TEST}")
endforeach()
message(STATUS "================================================")

# Create install script
configure_file(
    ${CMAKE_SOURCE_DIR}/README.md
    ${CMAKE_BINARY_DIR}/README.md
    COPYONLY
)
