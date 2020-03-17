set(TARGET gasmodule_objects)
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")
add_library(${TARGET} OBJECT ${SOURCES} ${HEADERS})

find_package(GSL REQUIRED)
# populates GSL_INCLUDE_DIRS

target_include_directories(${TARGET}
  PRIVATE ../../include
  SYSTEM PRIVATE ${EIGEN_INCLUDE_DIR} ${GSL_INCLUDE_DIRS})

add_library(gasmodule STATIC $<TARGET_OBJECTS:gasmodule_objects>)
target_link_libraries(gasmodule GSL::gsl GSL::gslcblas)

install(TARGETS gasmodule ARCHIVE DESTINATION lib)
install(FILES GasInterface.hpp GrainInterface.hpp GasState.hpp GasDiagnostics.hpp DESTINATION include)
