# Using the library

## API and code examples
Once the process described in [the build instructions](build.md) is successful, take note of the
install directory and take a look at the header files in the generated `include/` directory. Point
your compiler options, development environment, or toolchain to the correct locations, and try
one of these examples.

### Setup
The idea is that your code manages a single `GasInterface` object. When this object is
constructed, all of its members will be set up too, reading in all of the necessary data. While
it is safe for multiple of these objects to exist alongside each other, doing so should only be
done for experimentation purposes.

The arguments for the `GasInterface` constructor currently consist of three frequency grids.
They are called `iFrequencyv`, `oFrequencyv`, and `eFrequencyv`. They specify the frequency grid
for, respectively:
- the input radiation field
- the calculated opacity of the gas
- the calculated emissivity of the gas

In the snippet below, we initialize the gas module using the same frequency grid for all three
of these physical quantities. What the contents of `frequencyv` are depends on the needs of the
user, but generally speaking, the frequency range just above the H ionization threshold (91.2
nm) should be well sampled.

```
#include "GasInterface.hpp"
#include <valarray>

// Set up frequency grid
std::valarray<double> frequencyv;
// < code to to fill in the frequency grid >

// Initialize the gas module by creating a GasInterface object.
GasModule::GasInterface gasInterface(frequencyv, frequencyv, frequencyv);

// Make some place to store a gas state per cell.
std::vector<GasModule::GasState> gasStatev(numCells);
```

The last line of code allocates space to store the solution for a single computational element.
Allocating more than one of these is optional: You can also just use a single, temporary gas
state and extract what you need before discarding it.

### Solve for equilibrium in a cell
The main goal of this code is to calculate the equilibrium temperature and abundances for a
given density, radiation field, and a local set of grain properties.

To pass the grain properties, we need to use another class : `GrainInterface`. Your code should
gather all the relevant grain properties (grain effective radii, number densities and absorption
efficiencies) into such an object. The default constructor creates an empty object, representing
a medium with zero dust present. Grain populations can then be added using
`GrainInterface::addPopulation`. Read the documentation of this function (see the header of
`GrainInterface`) for details on how to pass the grain properties, and check out the example
further down.

In the example below, the radiation and density are obtained from elsewhere in the code, and an
empty `GrainInterface` is created. This info is passed to the `updateGasState` function of
`GrainInterface`, along with a writable reference to a `GasState` object, to which the most
important output will be written. The pseudocode `parallel_for_loop` emphasizes that running
`updateGasState` is thread safe, as long as different `GasState` objects are given.

```
parallel_for_loop(0, numCells, [&] (int cell) {
  // Mean intensity in frequency units (erg s-1 cm-2 sr-1 Hz-1), and the number density of H
  // nuclei in cm-3
  const std::valarray<double>& meanIntensityv = i_nuv(cell);
  double nH = hdensity(cell);

  // No grains
  const GasModule::GrainInterface grainInterface;

  // Solve for T
  gasInterface.updateGasState(gasStatev[cell], nH, Tguess, grainInterface);
}
```

### Dealing with grains
Essentially, `GrainInterface` is just a wrapper around a vector of `GrainInterface::Population`
objects. Populations can be added using `GrainInterface::addPopulation()`. Each population
represents a collection grains of the same type, with a certain grain size distribution and
absorption efficiency.

For every grain population, a flag specifying silicate or graphite needs to be provided. This
influences the grain photoelectric heating recipe, and determines some constants related to H2
formation on the surface of the grains. The absorption efficiency `Qabs(a, nu)` (\f$
(Q_\mathrm{abs}(a, \nu) \f$) is used to calculate the grain temperatures* and the photoelectric
efficiency. It needs to be provided separately, independent from the carbon or silicate flag. It
should be tabulated for every grain size (first index) and for every frequency of the frequency
grid used for the input radiation field (`iFrequencyv`).

> * While the grain temperatures need to be passed, the gas module actually recalculates the
> grain temperatures by itself. This is because the gas and the grains exchange energy through
> collisions and H2 formation. The given temperatures are only used for the initial guess.

```
#include "GrainInterface.hpp"

// Create grain interface
GasModule::GrainInterface grainInterface;

for (pop = 0; pop < numPop; pop++) {
   // Choose one of two supported types
   GasModule::GrainTypeLabel type = GasModule::GrainTypeLabel::CAR;

   // Gather these physical quantities
   std::valarray<double> sizev = grain_sizes(pop);
   std::valarray<double> densityv = grain_number_density_per_size(pop);
   std::valarray<double> temperaturev = grain_temperature_per_size(pop);
   std::vector<std::valarray<double>> qAbsvv = qabs_nu_per_size(pop);

   // Use the interface to add the population
   grainInterface->addPopulation(type, sizev, densityv, temperaturev, iFrequencyv, qAbsvv);
}
```

> Note: The example code above was written from the standpoint that the grains are completely
> different for each individual cell of a simulation. In most simulations, this is not the case.
> For SKIRT specifically, a fixed grid of sizes is used throughout all the cells, and only the
> densities change. The grain number densities of a population contained in an existing
> `GrainInterface` object can be adjusted safely using
> `GrainInterface::changePopulationDensityv()`. This is more efficient than making a new
> `GrainInterface` each time, since some precalculations are performed every time
> addPopulation() is called.
 
## Compiling and linking your program

### For simple scripts
If compiling your code is a very simple one-liner, use the include flag
`-I<CMAKE_INSTALL_PREFIX>/include`, where `<CMAKE_INSTALL_PREFIX>` is the install directory you
specified at configuration. Example:
```
gcc -c main.cpp -o main.o -I"$(HOME)"/RADAGAST/cmake_release/include
```

For the linking step, add the linker path and library names using the flags
`-L<CMAKE_INSTALL_PREFIX>/lib` and `-lradagast`. Your binary will also need to be linked to
GSL using `-lgsl` and `-lgslcblas`.
```
gcc main.o -o main -L"$(HOME)"/RADAGAST/cmake_release/lib -lradagast -lgsl -gslcblas
```
The order of the `-l` flags is important here for certain linkers (ld). First, ld needs to
process RADAGAST, so that it can make a list of symbols that it needs from GSL, and only then
should the GSL libaries be processed.

### For bigger software
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
Otherwise, you should ensure that all include and linker flags are set properly in your build system.
