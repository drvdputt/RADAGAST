# RIDGE-0
**R**adiation, **I**nterstellar **D**ust and **G**as in **E**quilibrium; a **0**D model.

## Description
RIDGE-0 is a local (0D) gas model for use in 3D dust radiative transfer simulations. Given a certain radiation field and dust model, the equilibrium state for various quantities of the gas are calculated self-consistently. Currently, the module treats the main species consisting of H (H+, H and H2), and their interactions with the radiation field and the dust grains.

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

While [SKIRT]("https://github.com/SKIRT/") is our main target, RIDGE-0 is designed as a stand-alone module which, ideally, does not depend on SKIRT in any way. It is possible to use this module in any C++ code where the necessary details about the radiation field and the dust can be provided, without having to install anything related to SKIRT.

This code was developed during my PhD research, and contains many equations and data from other authors to implement the processes listed above. Consult my PhD thesis (to be published) for details, equations, and citations. Alternatively, one can build the documentation pages or read the comments for each function in the source code.




