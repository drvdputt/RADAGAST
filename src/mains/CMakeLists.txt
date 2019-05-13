include_directories(../core ../../include)

set(TARGET main)
add_executable(${TARGET} main.cpp)
target_link_libraries(${TARGET} gasmodule)
include("../../CompilerFlags.cmake")
# the static or shared library already links to GSL

find_package(GSL REQUIRED) # populates GSL_INCLUDED_DIRS
set(TARGET run_tests)
add_executable(${TARGET} test.cpp $<TARGET_OBJECTS:gasmodule_objects>)
			 target_link_libraries(${TARGET} GSL::gsl GSL::gslcblas)
include("../../CompilerFlags.cmake")
# the testing framework needs the object files instead of the library used
# above, so we still need to link to GSL here
