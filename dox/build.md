# Build instructions

## Download
Download a zip of the source code or clone the repository. Due to the way the default build
works, it is best to have a dedicated parent directory.
```
mkdir RADAGAST
cd RADAGAST
mkdir git
git clone https://github.com/drvdputt/RADAGAST.git git
```

## Dependencies
- C++14 compiler (tested with of GCC and Clang)
- A recent enough version of CMake (> 3.1)
- <a href="https://www.gnu.org/software/gsl/">Gnu Scientific Library</a>, with development
  headers. As long as it is installed in a standard location, CMake will be able to find it on
  your system.
- <a href="http://eigen.tuxfamily.org/">Eigen3</a>. If Eigen is available via your package
  manager, install it and CMake will find it automatically in most cases. If not, download the
  latest release and put it at your preferred location. Then, see the CMake commands further
  down this page to make sure that the build system will find the correct location (sometimes
  this is necessary even if Eigen was installed via the package manager).

## Configure
Configuring the build is analogous to most projects that use CMake. In-source builds have been
explicitly disabled, so first create a build directory `RADAGAST/release`, as a sister directory
of `RADAGAST/git/`.

Then run CMake from this directory with the `git/` directory (containing the top level
`CMakeLists.txt`) as the main argument. This is also where build options can be provided using the `-D` option.
Once configured, the variables set for the build can be reviewed by using `cmake -L ../git`. 
```
cd release
cmake ../git -DCMAKE_RELEASE_TYPE=Release \
    -DEIGEN3_INCLUDE_DIR=/path/to/include/eigen3/ -DCMAKE_INSTALL_PREFIX="$(pwd)"
```
Specifying the 'Release' build type enables the default optimization flags for the compiler, and
compiles the code with most debugging output turned off. The last two options are useful when the dependencies are installed in non-standard locations.

### Eigen problems
By default, CMake will look for the Eigen3 headers in a handful of locations. If this
fails, provide the path to the `eigen3/` containing the headers by setting the
EIGEN3_INCLUDE_DIR variable, as shown above.

### GSL problems
If you experience header or linking problems related to GSL, it is usually it is installed in a non-standard location (i.e. the module FindGSL.cmake shipped with cmake can't find it).
In that case, use the option `-DGSL_ROOT_DIR=/your_gsl_install_prefix/`.

## Build and Install
The install prefix for RADAGAST can be chosen freely.
It determines where the library and its headers will be installed once the build has completed.
By default, this is `/usr/local/lib/` and `/usr/local/include/`, but you might want to change this to a location in your own home directory, to avoid polluting your system.

To start the build, run `$ make -j <number of jobs>` from the build directory.
If you want to use RADAGAST in another code, you can opt to run `cmake install` to copy some of the files to the
traditional `lib/` and `include/` directories (located in the directory defined by
`CMAKE_INSTALL_PREFIX`).

> Note: The script `build.sh` provides some reasonable defaults for the above actions. It builds and
> installs to `../cmake_release`.
