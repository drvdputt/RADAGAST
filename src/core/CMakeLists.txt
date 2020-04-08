find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)

file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# Build objects as a separate target; needed for building tests executable (see
# ../mains/CMakeLists.txt). Linking the GSL and Eigen3 targets takes care of their include
# directories.
add_library(gasmodule_objects OBJECT ${SOURCES} ${HEADERS})
target_include_directories(gasmodule_objects PRIVATE ../../include)
target_link_libraries(gasmodule_objects Eigen3::Eigen GSL::gsl GSL::gslcblas)

# static library for use in other projects
add_library(gasmodule STATIC $<TARGET_OBJECTS:gasmodule_objects>)
target_link_libraries(gasmodule Eigen3::Eigen GSL::gsl GSL::gslcblas)

# install the static libary and headers to the typical directories under CMAKE_INSTALL_PREFIX
install(TARGETS gasmodule ARCHIVE DESTINATION lib)
install(FILES GasInterface.hpp GrainInterface.hpp GasState.hpp GasDiagnostics.hpp DESTINATION include)
