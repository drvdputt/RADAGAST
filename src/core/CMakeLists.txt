set(TARGET gasmodule)

file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.h")

add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

include_directories(../../include)
include_directories(SYSTEM ../../eigen3)

find_package(GSL REQUIRED)
target_link_libraries(${TARGET} GSL::gsl GSL::gslcblas)

include("../../CompilerFlags.cmake")
