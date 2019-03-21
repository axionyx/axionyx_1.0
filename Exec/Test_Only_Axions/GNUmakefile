# AMREX_HOME defines the directory in which we will find all the BoxLib code
AMREX_HOME ?= ../../../amrex

HPGMG_DIR ?= ../../Util/hpgmg/finite-volume

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gnu

USE_MPI = TRUE
USE_OMP = FALSE

# Analysis
REEBER = FALSE
TINY_PROFILE = TRUE

USE_HPGMG = FALSE
USE_FDM = TRUE

PRECISION = DOUBLE
DEBUG     = FALSE
#DEBUG     = TRUE

# physics
DIM      = 3
USE_GRAV = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
