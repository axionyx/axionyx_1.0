export I_MPI_CXX="icpc"
export I_MPI_CC="icc"
export I_MPI_F90="ifort"

export FFTW_DIR="/sw/numerics/fftw3/impi/intel/3.3.8/skl/lib"
export LD_LIBRARY_PATH=$FFTW_DIR:$LD_LIBRARY_PATH

module load intel impi fftw3/impi/intel/3.3.8 hdf5-parallel/impi/intel/1.10.5
