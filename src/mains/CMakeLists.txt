include_directories(../core ../../include ${EIGEN_INCLUDE_DIR})

foreach(TARGET main krome_compare)
  add_executable(${TARGET} ${TARGET}.cpp)
  # this static library already links to GSL
  target_link_libraries(${TARGET} gasmodule)
  include("../../CompilerFlags.cmake")
endforeach()

find_package(GSL REQUIRED) # populates GSL_INCLUDED_DIRS
set(TARGET run_tests)
add_executable(${TARGET} test.cpp $<TARGET_OBJECTS:gasmodule_objects>)
			 target_link_libraries(${TARGET} GSL::gsl GSL::gslcblas)
include("../../CompilerFlags.cmake")
# the testing framework needs the object files instead of the library used
# above, so we still need to link to GSL here
