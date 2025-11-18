PROJECT(eel2d)
CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

# Source files for the undulatory foil simulation
# with Reynolds number and thickness effects
SET(SOURCE_FILES example.cpp IBEELKinematics.cpp)
SET(HEADER_FILES IBEELKinematics.h)

ADD_EXECUTABLE(main2d ${SOURCE_FILES} ${HEADER_FILES})

FIND_PACKAGE(IBAMR REQUIRED)
TARGET_LINK_LIBRARIES(main2d IBAMR::IBAMR2d)

# Set C++ standard
TARGET_COMPILE_FEATURES(main2d PRIVATE cxx_std_11)
