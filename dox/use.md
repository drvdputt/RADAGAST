# Using the library

## API and code examples

## Compiling and linking your program

### For one-liners
If compiling your code is a very simple one-liner, use the include flag
`-I<CMAKE_INSTALL_PREFIX>/include`, where `<CMAKE_INSTALL_PREFIX>` is the install directory you
specified at configuration. Example:
``` gcc -c main.cpp -o main.o -I"$(HOME)"/RIDGE-0/cmake_release/include ```
For the linking step, add the linker path and library names using the flags @c
-L\<CMAKE_INSTALL_PREFIX\>/lib and @c -lridge. Your binary will also need to be linked to
GSL using @c -lgsl and @c -lgslcblas.
``` gcc main.o -o -L"$(HOME)"/RIDGE-0/cmake_release/lib -lridge -lgsl -gslcblas ```
The order of the `-l` flags is important here for certain linkers (ld). First, ld needs to
process RIDGE-0, so that it can make a list of symbols that it needs from GSL, and only then
should the GSL libaries be processed.

### CMake
If your program is also built using CMake, then you can take inspiration from the following
snippet

```
  set(GASMODULE_INSTALL_DIR "/your/install/dir")
  find_library(GASMODULE_LIBRARIES
    NAMES libgasmodule.a
    PATHS "${GASMODULE_INSTALL_DIR}/lib")
  find_path(GASMODULE_INCLUDE_PATH
    NAMES "GasInterface.hpp"
    PATHS "${GASMODULE_INSTALL_DIR}/include")
  include_directories(${GASMODULE_INCLUDE_PATH})
  target_link_libraries(${TARGET} ${GASMODULE_LIBRARIES})

  find_package(GSL REQUIRED)
  target_link_libraries(${TARGET} GSL::gsl GSL::gslcblas)
```
