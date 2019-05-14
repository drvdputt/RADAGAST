file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.h")
find_package(GSL REQUIRED) # populates GSL_INCLUDED_DIRS

set(TARGET gasmodule_objects)
add_library(${TARGET} OBJECT ${SOURCES} ${HEADERS})
target_include_directories(${TARGET} PRIVATE ../../include)
target_include_directories(${TARGET} SYSTEM PRIVATE ../../eigen3 ${GSL_INCLUDE_DIRS})
include("../../CompilerFlags.cmake")
# do not forget to link to GSL when making a binary from these objects

set(TARGET gasmodule)
add_library(gasmodule STATIC $<TARGET_OBJECTS:gasmodule_objects>)
target_link_libraries(${TARGET} GSL::gsl GSL::gslcblas)
include("../../CompilerFlags.cmake")
# here we already linked GSL into the static libary