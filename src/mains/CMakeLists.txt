find_package(GSL REQUIRED)

# configure multiple executables the same way
foreach(TARGET main krome_compare make_wavelengthgrid)
  add_executable(${TARGET} ${TARGET}.cpp)
  target_link_libraries(${TARGET} radagast)
endforeach()

# Configuration for tests executable is special; it won't work with the statically linked
# 'radagast' archive. Instead, it needs the list of objects from the radagast_objects target
# in ../core. Consequently, GSL still needs to be linked here.
add_executable(run_tests test.cpp $<TARGET_OBJECTS:radagast_obj>)
target_include_directories(run_tests PRIVATE ../../include)
target_link_libraries(run_tests GSL::gsl GSL::gslcblas)
