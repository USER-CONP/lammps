# Warning

This implementation of the constant potential method is deprecated.
You can find the new version here: 

https://github.com/robeme/lammps/tree/electrode

# Constant Potential on a Mesh (and more)

This version of the LAMMPS software package features an implementation of the
constant potential method (CPM).

The fix charge_update implements the CPM. Charges of groups specified as
group-ID and with the couple keyword are adapted to meet their respective
potential at every time step. An arbitrary number of electrodes can be set but
the respective groups may not overlap. Electrode charges are smeared with a
Gaussian charge distribution with reciprocal width eta. The energy minimization
is achieved via matrix inversion.

``` fix [ID] [group-ID] charge_update [ψ] [η] [optional keyword] [optional value] ...  ```

`ID` = ID of FIX command

`group-ID` = Group for 'first' electrode

`η` = Gaussian width parameter in angstrom<sup>-1</sup>

`ψ` = charges such that group-ID has a potential ψ volts.

## The optional keywords and values allowed are as follows:

`couple [group-ID2] [ψ]` = second electrode coupled to the first that has a
potential of ψ volts.

`write_inv [file name]` = store inverse of electrode-electrode matrix in a file
for later reuse.

`read_inv [file name]` = read inverse of electrode-electrode matrix.

`symm [on/off]` = use symmetric matrix as defined in Scalfi et al. (PCCP, 2020, 22, 10480) to ensure
charge neutrality.

The fix necessitates the use of a long range solver that can provide the matrix
of electrode-electrode interactions and a vector of electrode-electrolyte
interactions. The Kspace styles ewald/conp and pppm/conp are created
specifically for this task.

For systems with non-periodic boundaries in one or two directions dipole
corrections are available with the kspace_modify. For ewald/conp a
two-dimensional Ewald summation can be used by setting “slab ew2d”:

Specific examples on how to use this code be found in the directory
"./examples/USER/conp/".

In order to install the package with cmake use "cmake -D PKG_USER-CONP=yes ."
and make sure the KSPACE package is installed. A traditional make with
"make-yes user-conp" requires that the linalg library is previously build.  The
linalg library is found in lib/linalg. See the readme therein for more details
on how to compile it.

More information on this specific implementation can be found under: 

https://doi.org/10.1063/5.0063381

----------------------------------------------------------------------

This is the LAMMPS software package.

LAMMPS stands for Large-scale Atomic/Molecular Massively Parallel
Simulator.

Copyright (2003) Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software.  This software is distributed under
the GNU General Public License.

----------------------------------------------------------------------

LAMMPS is a classical molecular dynamics simulation code designed to
run efficiently on parallel computers.  It was developed at Sandia
National Laboratories, a US Department of Energy facility, with
funding from the DOE.  It is an open-source code, distributed freely
under the terms of the GNU Public License (GPL).

The primary author of the code is Steve Plimpton, who can be emailed
at sjplimp@sandia.gov.  The LAMMPS WWW Site at lammps.sandia.gov has
more information about the code and its uses.

The LAMMPS distribution includes the following files and directories:

- README                     this file
- LICENSE                    the GNU General Public License (GPL)
- bench                      benchmark problems
- cmake                      CMake build files
- doc                        documentation
- examples                   simple test problems
- fortran                    Fortran wrapper for LAMMPS
- lib                        additional provided or external libraries
- potentials                 interatomic potential files
- python                     Python wrappers for LAMMPS
- src                        source files
- tools                      pre- and post-processing tools

Point your browser at any of these files to get started:

https://lammps.sandia.gov/doc/Manual.html         LAMMPS manual\
https://lammps.sandia.gov/doc/Intro.html          hi-level introduction\
https://lammps.sandia.gov/doc/Build.html          how to build LAMMPS\
https://lammps.sandia.gov/doc/Run_head.html       how to run LAMMPS\
https://lammps.sandia.gov/doc/Commands_all.html   Table of available commands\
https://lammps.sandia.gov/doc/Library.html        LAMMPS library interfaces\
https://lammps.sandia.gov/doc/Modify.html         how to modify and extend LAMMPS\
https://lammps.sandia.gov/doc/Developer.html      LAMMPS developer info\

You can also create these doc pages locally:

% cd doc\
% make html                # creates HTML pages in doc/html\
% make pdf                 # creates Manual.pdf
