file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

add_library(gasmodule_obj OBJECT ${SOURCES} ${HEADERS})
target_include_directories(gasmodule_obj PUBLIC ../../include INTERFACE ./)

find_package(Eigen3 3.3 NO_MODULE QUIET)
if(TARGET Eigen3::Eigen)
  # I wanted to use target_link_libraries here, but this is not supported for object libaries in
  # cmake 3.5. So instead, extract the INCLUDE_DIRECTORIES property from the Eigen target using
  # the following generator expression, and use PUBLIC to make sure that the includes are
  # propagated to dependents of gasmodule_obj.
  target_include_directories(gasmodule_obj SYSTEM PUBLIC $<TARGET_PROPERTY:Eigen3::Eigen,INCLUDE_DIRECTORIES>)
elseif(EIGEN3_INCLUDE_DIR)
  # if no package was found, use this variable instead
  message("Using eigen at EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}")
  target_include_directories(gasmodule_obj SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})
else()
  # if both fail, make some suggestions and quit
  message(FATAL_ERROR
    "Could not find Eigen automatically. You have two options.\n"
    " 1) Install an eigen3 package that comes with Eigen3Config.cmake, and set Eigen3_DIR to the directory containing that file if necessary.\n"
    " 2) Download and decompress eigen3, and set -DEIGEN3_INCLUDE_DIR to the location where you chose to put 'eigen3/'.")
endif()
  
find_package(GSL REQUIRED)
target_include_directories(gasmodule_obj SYSTEM PUBLIC
  $<TARGET_PROPERTY:GSL::gsl,INCLUDE_DIRECTORIES>
  $<TARGET_PROPERTY:GSL::gslcblas,INCLUDE_DIRECTORIES>)

add_library(gasmodule $<TARGET_OBJECTS:gasmodule_obj>)
target_include_directories(gasmodule INTERFACE $<TARGET_PROPERTY:gasmodule_obj,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(gasmodule GSL::gsl GSL::gslcblas)
  
# install static libary and headers to the typical directories under CMAKE_INSTALL_PREFIX
install(TARGETS gasmodule ARCHIVE DESTINATION lib)
install(FILES GasInterface.hpp GrainInterface.hpp GasState.hpp GasDiagnostics.hpp DESTINATION include)
