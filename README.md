# RDGST
Radiation and Dust driven equilibrium of interstellar Gas as a STandalone library

## Description
RDGST is a local (0D) gas model for use in 3D dust radiative transfer simulations. Given a certain radiation field and dust model, the equilibrium state for various quantities of the gas are calculated self-consistently. Currently, the module only treats the main species consisting of H, and various interaction with the radiation field and the dust grains.

Quantities calculated self-consistently
- Gas equilibrium temperature
- H+, H, and H2 abundance
- H level populations up to n = 5
- H2 level populations
- Grain temperatures and charge distributions

A non-comprehensive list of processes included in this calculation
- Photoionization of H 
- Radiative recombination of H+ 
- Photodissociation of H2 
- Formation of H2 on dust grain surfaces
- Heating by the photoelectric effect on dust grains
- Cooling by collisions between gas and dust grains
- Cooling by collisional transitions of H2 (and H at high temperatures)

From the resulting state, the emissive properties and cross section of the gas can be calculated, so that a radiative transfer code can include these in its model for the radiative transfer medium.
- Cross sections for H photoionization, H2 photoionization, and H2 photodissociation
- Continuum of H (free-bound, free-free, and two-photon)
- Emission coefficients and opacities for the lines of H and H2.

While <a href="https://github.com/SKIRT/">SKIRT</a> is our main target, we
have made this a standalone module which, ideally, does not depend on SKIRT in any way. It
should be possible to use this module in any C++ code where the necessary details about the radiation field and the dust 
can be provided.

This code was developed during my PhD research, and the thesis contains all the equations and references for the implementation of the processes listed above. Alternatively, one can build the documentation pages or read the comments for each function in the source code.
