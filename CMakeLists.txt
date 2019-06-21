cmake_minimum_required(VERSION 3.5)
project(gasmodule)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

add_definitions(-DREPOROOT=\"${PROJECT_SOURCE_DIR}\")
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_definitions(-DSILENT)
endif()

SET(EIGEN_POSSIBLE_LOCATIONS /usr/local/include/eigen3 /usr/include/eigen3)
foreach(EIGEN_DIR IN LISTS EIGEN_POSSIBLE_LOCATIONS)
  if (EXISTS ${EIGEN_DIR})
    SET(EIGEN_INCLUDE_DIR ${EIGEN_DIR} CACHE PATH "Eigen include directory")
    break()
  endif()
endforeach()


add_subdirectory(src/core)
add_subdirectory(src/mains)
