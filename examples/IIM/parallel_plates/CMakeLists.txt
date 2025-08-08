PROJECT(mysolver)
CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

ADD_EXECUTABLE(plates2d plates_example.cpp)
ADD_EXECUTABLE(exe read_file.C)
FIND_PACKAGE(IBAMR 0.11.0 REQUIRED)
TARGET_LINK_LIBRARIES(plates2d IBAMR::IBAMR2d)
TARGET_LINK_LIBRARIES(exe IBAMR::IBAMR2d)
# IBAMR saves the flags it used to compile - you can reuse them if you want to
SET(CMAKE_CXX_FLAGS ${IBAMR_CXX_FLAGS})

