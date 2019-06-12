FFTW_DIR ?= /usr/users/cbehren2/local_libs/lib

# AMREX_HOME defines the directory in which we will find all the BoxLib code
AMREX_HOME ?= ../../../amrex_pure

HPGMG_DIR ?= ../../Util/hpgmg/finite-volume

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = intel #gnu  # intel

USE_MPI = TRUE
USE_OMP = FALSE

# Analysis
REEBER = FALSE

USE_HPGMG = FALSE
USE_FDM = TRUE

PRECISION = DOUBLE
DEBUG     = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE
#USE_CGRAV = TRUE

Bpack := ./Make.package
Blocs := .

#some stuff for using fftw3

LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3_omp -lfftw3
INCLUDE_LOCATIONS       += $(FFTW_DIR)/../include

include $(AMREX_HOME)/Src/Extern/SWFFT/Make.package
INCLUDE_LOCATIONS       += $(AMREX_HOME)/Src/Extern/SWFFT
VPATH_LOCATIONS         += $(AMREX_HOME)/Src/Extern/SWFFT

include $(TOP)/Exec/Make.Nyx
