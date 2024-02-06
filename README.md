# ZMC
Fortran code for harmonic modelling of crystalline systems

Overview
--------

Conventional crystal structure refinement and solution relies on an analysis of the Bragg peaks, the sharp, well-defined scattering in a diffraction pattern. In doing this, the structure can be considered as an array of identical unit cells, and this reduces the 'solution' of the crystal structure to the determination of the unit cell contents.

In a disordered material, while the average unit cell does indeed obey the space group symmetry of the crystal structure, such that solving the asymmetric unit within the cell 'solves' the structure, this is not the case for a disordered material -- or even for what would generally be considered as an ordered material if considering instantaneous molecular and atomic motions.

Even an ordered crystal will show thermally induced vibrations of the atoms. If we could 'freeze' the crystal and look at the atoms, we would find that the space group symmetry is not obeyed. It is obeyed on average. In a disordered material it could be that the average occupation of a given site is 50% one type of atom and 50% another, in which case no site obeys the average symmetry.

What does this mean for analysis of diffuse scattering? It means that you can no longer consider all unit cells as identical. It means that you have to now work with an ensemble of unit cells big enough that averaging across it recovers the average and big enough that it can contain a statistically valid population of the local (short-range order) structures present in the crystal.

This means that analysis of diffuse scattering requires different tools from analysis of conventional Bragg scattering. ZMC is an attempt at a program to allow relatively ready implementation of a simulation of diffuse scattering from (particularly but not solely) molecular crystals. 

Current status
--------------

I have some old binaries that seem to work on Windows and Linux, and I want to test them before I post them. The code is so old that some changes to gfortran mean a key module no longer compiles. I was hoping to fix that before I uploaded, but I might upload anyway. It compiles on older versions of the compiler (circa 2016), bit no longer.

**More to come**
