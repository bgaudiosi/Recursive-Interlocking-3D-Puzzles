
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(puzzle)

set (CMAKE_CXX_FLAGS "-std=c++11")

file(GLOB_RECURSE HEADER_CODE ${puzzle_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE SRC_CODE ${puzzle_SOURCE_DIR}/src/*.cpp)

ADD_EXECUTABLE(puzzle ${SRC_CODE} ${HEADER_CODE})
