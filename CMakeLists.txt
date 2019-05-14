cmake_minimum_required(VERSION 3.7)

project(gasmodule)

set(CMAKE_BUILD_TYPE Release)

add_definitions(-DREPOROOT=\"${PROJECT_SOURCE_DIR}\")

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_definitions(-DSILENT)
endif()

add_subdirectory(src/core)
add_subdirectory(src/mains)
