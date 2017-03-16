The 'include' directory contains some headers I stole from either SKIRT or
DIRTY. 

The Eigen headers are in a separate directory 'eigen3' so that the warnings
coming from it can be suppressed more easily. Currently, all of eigen is
included, making it easier to update/refresh in case I break something. In the
future it might be better to throw away everything we're not using, once I know
wat I need exactly from Eigen. 

The makefile doesn’t properly include all of the dependencies yet, but it
probably won’t be our actual build tool. Currently, it builds a binary called
‘test’ which executes whatever was written in src/main.cpp. The object files are
dumped at ../obj, and the binary is placed at ../bin.

Instructions for adding this project as a subproject in CMake:
	
	- In the CMakeLists.txt of the parent project, put
	add_subdirectory(gasmodule).

	- Create a subdirectory called 'gasmodule' next to the other
	subdirectories of the parent project. There put a CMakeLists.txt
	containing something analogous to the following:
"""
# ------------------------------------------------------------------
# Builds a library from the GasModule source
# ------------------------------------------------------------------

# set the target name
set(TARGET gasmodule)

# source files and their headers
file(GLOB SOURCES "/Users/drvdputt/GasModule/git/src/*.cpp")
file(GLOB HEADERS "/Users/drvdputt/GasModule/git/src/*.hpp")
# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# Eigen dependencies
include_directories(SYSTEM /Users/drvdputt/GasModule/git/eigen3)

# GasModule header-only dependencies
include_directories(/Users/drvdputt/GasModule/git/include)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
"""
	The SYSTEM flag for the eigen3 include will suppress the  thousands of
	warnings thrown (in our case, because of compiler flags).
	
	- If any of the subprojects make use of the gasmodule, put the following
	in their CMakeLists.txt:
"""
target_link_libraries(gasmodule)
include_directories(/Users/drvdputt/GasModule/git/src)
""" 