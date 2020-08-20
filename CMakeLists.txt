cmake_minimum_required(VERSION 3.1)
project(RADAGAST)

set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wall -Wextra -Wno-missing-braces -Wno-sign-compare -Werror=return-type -pedantic -fno-signed-zeros -freciprocal-math -march=native -mtune=native)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.\n\
Please remove CMakeCache.txt and call cmake from another directory")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

add_definitions(-DREPOROOT=\"${PROJECT_SOURCE_DIR}\")
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_definitions(-DSILENT)
endif()
add_definitions(-march=native)

add_subdirectory(src/core)
add_subdirectory(src/mains)
