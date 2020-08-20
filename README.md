# RADAGAST
**RAD**iation and **D**ust **A**ffecting the **GA**s **ST**ate

## Description
RADAGAST is a local (0D) gas model for use in 3D dust radiative transfer simulations. Given a
certain radiation field and dust model, the equilibrium state for various quantities of the gas
are calculated self-consistently. Currently, the module treats the main species consisting of H
(H+, H and H2), and their interactions with the radiation field and the dust grains.

Quantities calculated self-consistently
- Gas equilibrium temperature
- H+, H, and H2 abundance
- Level populations of H and H2
- Grain temperatures and charge distributions

A non-comprehensive list of processes included in this calculation
- Photoionization of H
- Radiative recombination of H+
- Photodissociation of H2
- Formation of H2 on dust grain surfaces
- Heating by the photoelectric effect on dust grains
- Cooling by collisions between gas and dust grains
- Net heating/cooling by collisional transitions of H2 (and H at high temperatures)

From the resulting state, the emissive properties and cross section of the gas can be calculated, so that a radiative transfer code can include these in its model.
- Cross sections for H photoionization, H2 photoionization, and H2 photodissociation
- Continuum emission of H (free-bound, free-free, and two-photon)
- Emission coefficients, opacities, and shapes for the lines of H and H2.

While RADAGAST was mainly developed for use in the 3D Monte Carlo Dust Radiative Transfer code <a
href="http://www.skirt.ugent.be">SKIRT</a>, it is designed as a stand-alone module which,
ideally, does not depend on SKIRT in any way. It is possible to use this module in any C++ code
where the necessary details about the radiation field and the dust can be provided, without
having to install anything related to SKIRT.

This code was developed during my PhD research, and contains many equations and data from other authors to implement the processes listed above. Consult my PhD thesis (to be published) for details, equations, and citations. Alternatively, one can build the documentation pages or read the comments for each function in the source code.

## Third party libraries and data
RADAGAST includes a release of the <a href="https://github.com/onqtam/doctest">doctest</a> C++
testing framework, to power a small set of unit tests.

It also uses <a href="https://eigen.tuxfamily.org">Eigen3</a> for storing vectors and matrices,
and performing linear algebra, and <a href="https://www.gnu.org/software/gsl/">GSL</a> for
solving several non-linear equations. Currently, these libraries are dependencies, but the git
history might contain a copy.

## Getting started
First read how to [build](dox/build.md) the library.

Then try running the main binary produced at `RADAGAST/../cmake_release/src/mains/main`.
Experiment by giving a few numerical command line arguments to this program. The command
``` main 1e5 1e4 1e3 ```
will run an example model with a gas density of 1e5 cm-3, under the effect
of a blackbody source located at a (fixed) distance of 1 pc with a color temperature of 1e4 K and a
bolometric luminosity of 1e3 solar luminosities.

When this binary works without problems, copy the library built at
`RADAGAST/../cmake_release/src/core/libridge.a` to your preferred location (or set
`CMAKE_INSTALL_PREFIX` and run make install). Then read the [instructions](dox/use.md) on how to
use the API, write a simple C++ script making the calls, and try building and linking your own
C++ script to the library.

## Integration in SKIRT
There currently exists a branch in my fork of SKIRT which supports gas and makes use of
RADAGAST. It is currently not up to data with the master branch of SKIRT, because some new
features have introduces a few conflicts. The reader interested in how RADAGAST is used in
practice, can examine this development branch: https://github.com/drvdputt/SKIRT9/tree/gas-dev
