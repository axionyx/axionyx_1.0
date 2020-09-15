axionyx is a modification of the the [Nyx code](https://github.com/AMReX-Astro/Nyx.git) for dealing with Fuzzy Dark Matter. The AxioNyx paper https://arxiv.org/abs/2007.08256 was produced using the mixed_dm branch.

# Contents
[Documentation](#Documentation)  
[Installation](#installation)  
[File structure](#file-structure-of-nyx)  
[More file structure](#more-file-structure)   
[Call graph of Axionyx](#call-graph-of-axionyx)  
[Output](#output)  

# Documentation 
A detailed documentation is not available yet. But you can find an API reference here:  
[Available input parameters](mddocs/input_parameters.md)  
[Available test problems](mddocs/test_problems.md)  
[Classes](mddocs/classes.md)  
[Functions](mddocs/functions.md)    
[Files](mddocs/files.md)     

# Installation

clone both the axionyx and the amrex4axionyx repo. axionyx expects amrex to live in ../amrex. You will need FFTW3 for axionyx to compile. Try building e.g. Exec/Test_FDM_lin_growth by running make. Adapt the GNUmakefile as needed.

## Installation in 5 minutes

Clone both repositories and switch branches

    git clone https://github.com/axionyx/amrex4axionyx_1.0.git amrex
    git clone https://github.com/axionyx/axionyx_1.0.git
    cd axionyx_1.0
    git checkout mixed_dm

Configure your compiler. You can do so by editing the GNUMakefile in any of the Exec/Test_FDM_* problems. To tell axionyx where to find FFTW3, you might want to set 

    export FFTW_DIR=path/to/fftw/lib
To compile one of the test problems run
   
    cd Exec/Test_FDM_condensation
    make
And to do a quick run, you can use  
  
    mpirun -np 1 ./Nyx.3d.[...] inputs




    


# File structure

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


# License
axionyx is released under the LBL's modified BSD license, see the [license.txt](license.txt) file for details.


