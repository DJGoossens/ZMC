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

Some changes to gfortran broke the program due to some non-robust code in managing things like variable length strings and hash tables. I think gfortran just tighened up a bit and some undefined behaviour got defined in a way that did not suit the program. This has been rectified and it now (March 2026) compiles and runs. Having said that, it has not been exhaustively tested, so there may be routines that my tests have not accessed that are still noncompliant. Sorry.

This archve (see Releases) contains Linux X86_64 (AMD64) binaries compiled on Debian 13. The program has successfully been compiled and run using gfortran on:
* Debian 13 on x86_64
* Cygwin on Windows 10 and Windows 11
* MSYS2 on Windows 10
* Debian 13 inside WSL2 on Windows 11
* Haiku OS (Shredder) on x86_64.

I can report that all of MSYS2, Cygwin and WSL2 give very similar performance for the same underlying hardware. The advantage of WSL is that it ought to run the binaries in the Release without recompiling.

I expect it to compile everywhere gfortran runs. There is not a compile script for non-bash-type shells, but converting it to a Windows BAT file would not be hard. See compile_all_Linux_gfortran.sh. There is no Makefile.

It also compiles and runs on 32-bit Debian on x86 (very old Compaq laptop), Debian 'armel' (old Raspberry Pi) and 64-bit Linux on an ancient DEC Alpha 21164 (EV56) (Debian 14 'sid/forky'). These were just tests if new errors arose. So far, none that matter.

Relevant links
--------------

An article from 2015 that goes into significant detail about how to set up a simulation can be found at:

https://doi.org/10.1155/2015/878463

and

https://onlinelibrary.wiley.com/doi/10.1155/2015/878463

The simulation that was included with that publication and all the associated code and a document to explain it can be found at:

https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1155%2F2015%2F878463&file=acmp878463-sup-0002-f2.zip

This ZIP file has also been added to this github archive. See acmp878463-sup-0002-f2.zip.
