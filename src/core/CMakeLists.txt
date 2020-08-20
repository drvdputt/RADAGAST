file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

add_library(radagast_obj OBJECT ${SOURCES} ${HEADERS})
target_include_directories(radagast_obj PUBLIC ../../include INTERFACE ./)

find_package(Eigen3 3.3 NO_MODULE QUIET)
if(TARGET Eigen3::Eigen)
  # I wanted to use target_link_libraries here, but this is not supported for object libaries in
  # cmake 3.5. And EIGEN3_INCLUDE_DIR set by find_package(Eigen3) is buggy (Eigen3Config.cmake
  # in eigen source code states *do not use*). So I use the target property instead.
  message("Found Eigen ${EIGEN3_VERSION_STRING} package in ${EIGEN3_ROOT_DIR}")
  get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
elseif(EIGEN3_INCLUDE_DIR)
  # if no package was found, try using manually set value
  message("EIGEN3_INCLUDE_DIR was set to ${EIGEN3_INCLUDE_DIR} manually")
endif()

# One last thing: on some systems, the path to an eigen module (e.g. Dense) looks like this:
# somewhere/eigen3/Eigen/Dense. On others, its somwhere/Eigen/Dense, omitting the 'eigen3' part.
# We use the common part, Eigen/Dense, in the code, so in the former case, we should add eigen3
# to the include dir.
if (IS_DIRECTORY ${EIGEN3_INCLUDE_DIR}/eigen3/Eigen)
  set(REAL_EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR}/eigen3)
elseif(IS_DIRECTORY ${EIGEN3_INCLUDE_DIR}/Eigen)
  set(REAL_EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR})
else()
  message(FATAL_ERROR
    "Could not find Eigen in EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}. You have two options.\n"
    " 1) Install an eigen3 package that comes with Eigen3Config.cmake, and set Eigen3_DIR to the directory containing that file if necessary (usually somwhere/share/cmake).\n"
    " 2) Download and decompress eigen3, and set -DEIGEN3_INCLUDE_DIR to the parent directory of 'eigen3' or 'Eigen'.")
endif()

target_include_directories(radagast_obj SYSTEM PUBLIC ${REAL_EIGEN3_INCLUDE_DIR})

find_package(GSL REQUIRED)
target_include_directories(radagast_obj SYSTEM PUBLIC
  $<TARGET_PROPERTY:GSL::gsl,INTERFACE_INCLUDE_DIRECTORIES>
  $<TARGET_PROPERTY:GSL::gslcblas,INTERFACE_INCLUDE_DIRECTORIES>)

add_library(radagast $<TARGET_OBJECTS:radagast_obj>)
# make sure that include directories are inherited
target_include_directories(radagast INTERFACE $<TARGET_PROPERTY:radagast_obj,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(radagast GSL::gsl GSL::gslcblas)

# install static libary and headers to the typical directories under CMAKE_INSTALL_PREFIX
install(TARGETS radagast ARCHIVE DESTINATION lib)
install(FILES GasInterface.hpp GrainInterface.hpp GasState.hpp GasDiagnostics.hpp DESTINATION include)
