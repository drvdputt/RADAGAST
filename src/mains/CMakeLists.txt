set(TARGET main)

add_executable(${TARGET} main.cpp)

include_directories(../core ../../include)

target_link_libraries(${TARGET} gasmodule)

include("../../CompilerFlags.cmake")
