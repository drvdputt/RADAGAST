cmake_minimum_required(VERSION 3.5)

project(gasmodule)

add_definitions(-DREPOROOT=\"${PROJECT_SOURCE_DIR}\")
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_definitions(-DSILENT)
endif()

SET(EIGEN_INCLUDE_DIR /usr/local/include/eigen3 CACHE PATH "Eigen include directory")

add_subdirectory(src/core)
add_subdirectory(src/mains)
