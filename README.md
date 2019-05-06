# Contents
[Installation](#installation)

[File structure](#file-structure-of-nyx)

[More file structure](#more-file-structure)

[Output](#output)

# Installation

This installs Nyx with the Intel compilers.

`mkdir nyx`  
`cd nyx`  
`#clone amrex`  
`git clone `<https://github.com/AMReX-Codes/amrex.git>  
`#clone Nyx`  
`git clone `<https://github.com/AMReX-Astro/Nyx.git>  
`#load intel compiler and MPI`  
`module load intel/compiler/64/2017/17.0.2`  
`module load intel/mpi/64/current/2.174`  
`#necessary to force the intel compilers to be used `  
`export I_MPI_CXX=icpc`  
`export I_MPI_CC=icc`  
`export I_MPI_F90=ifort`  
  
`#build amrex`  
`cd amrex`  
`./configure --comp intel --with-omp yes`  
`#note; I think it is technically not necessary to build amrex. Whatever needs to be built will be built when we compile Nyx. `  
`#However, it is not a bad practice to build amrex first, because it makes it easier to  isolate problems with the environment.`  
`make -j12`  
`#build a random problem in Nyx - note that not all predefined problems with actually compile!`  
`cd ../Nyx/Exec/AMR-density`  
`make -j12 `

A typical problem that one encounters is that the compiler is set to gcc
for a given problem. You need to adapt the corresponding GNUMakefile in
this case.

Nyx expects amrex to exist in <nyx repo>/..; if it is somewhere else,
you need to set the AMREX\_HOME variable accordingly. on startup.

# File structure of Nyx

  - Source/ contains the baseline source files
  - the actually compileable problems (executables) live in sub
    directories of Exec/
  - each sub dir in Exec/ usually has
      - some default *inputs* file, defining the input parameters for a
        simulation run
      - some default *probin* file, defining the input parameters
        exclusively used in Fortran code (when for some reason the
        inputs file's content can be passed down)
      - a *GNUMakefile*, defining how this problem is to be compiled;
        most importantly, in this file you can switch on/off OpenMP
        support, the debug mode, and change the compiler used (with the
        COMP parameter) to gnu/intel/..
      - a *Make.package* file, defining which files need to be compiled.
        You need to add here files only if you are not overriding the
        ones in Source, e.g. if there is x.cpp in Source/ that you
        rewrote and put in the problem folder, it will automatically
        replace the one in Source. However, if you add a 'new' file, you
        need to add it to the list in Make.package; for example, to add
        y.f90 to the list, you would add a line

        `f90EXE_sources += y.f90`

  -   - a *Prob\_3d.f90* file; it contains the code to initialize the
        problem that you want to solve, e.g. sets up the density field
        etc. The important routine here is called fort\_initdata.
      - many problems also have a *Tagging\_3d.f90* and a
        *Nyx\_error.cpp* file; the first defines the criteria for
        refining cells, the latter adds these criteria to the existing
        list of such criteria.

# More file structure

TBD

# Output

Nyx outputs certain global diagnostics at each timestep and plot files at regular
intervals, or at user-specified redshifts. Visualization packages
[VisIt](https://wci.llnl.gov/simulation/computer-codes/visit),
[Paraview](https://www.paraview.org/)
and [yt](http://yt-project.org/)
have built-in support for the AMReX file format used by Nyx.

In addition, Nyx interfaces with two post-processing suites, Reeber and Gimlet. Reeber
uses topological methods to construct merge trees of scalar fields, which Nyx in
turn uses to find halos. Gimlet computes a variety of quantities
related to the Lyman-alpha forest science. These suites are fully MPI-parallel and can
be run either "in situ" or "in-transit", or with a combination of both.


# License
Nyx is released under the LBL's modified BSD license, see the [license.txt](license.txt) file for details.


