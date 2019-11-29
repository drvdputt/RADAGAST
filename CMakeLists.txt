cmake_minimum_required(VERSION 3.1)
project(gasmodule)

set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wshadow -Wall -Wextra -Wno-missing-braces -Wno-sign-compare -Werror=return-type -pedantic -fno-signed-zeros -fassociative-math -freciprocal-math -march=native -mtune=native)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

add_definitions(-DREPOROOT=\"${PROJECT_SOURCE_DIR}\")
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_definitions(-DSILENT)
endif()

add_definitions(-march=native)

SET(EIGEN_POSSIBLE_LOCATIONS /usr/local/include/eigen3 /usr/include/eigen3 $ENV{HOME}/.local/eigen3)
message(STATUS "Looking for Eigen in the list ${EIGEN_POSSIBLE_LOCATIONS}")
foreach(EIGEN_DIR IN LISTS EIGEN_POSSIBLE_LOCATIONS)
  if (EXISTS ${EIGEN_DIR})
    SET(EIGEN_INCLUDE_DIR ${EIGEN_DIR} CACHE PATH "Eigen include directory")
    break()
  endif()
endforeach()
message(STATUS "Using eigen headers from ${EIGEN_DIR}")

add_subdirectory(src/core)
add_subdirectory(src/mains)
