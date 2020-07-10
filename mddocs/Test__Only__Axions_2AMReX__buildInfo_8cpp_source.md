
# File AMReX\_buildInfo.cpp

[**File List**](files.md) **>** [**Exec**](dir_43a12cefb7942b6f49b5b628aafd3192.md) **>** [**Test\_Only\_Axions**](dir_eb24725df855cf6c732a19e4912f662a.md) **>** [**AMReX\_buildInfo.cpp**](Test__Only__Axions_2AMReX__buildInfo_8cpp.md)

[Go to the documentation of this file.](Test__Only__Axions_2AMReX__buildInfo_8cpp.md) 


````cpp

namespace amrex {

const char* buildInfoGetBuildDate() {

  static const char BUILD_DATE[] = "2019-07-16 13:48:54.980096";
  return BUILD_DATE;
}

const char* buildInfoGetBuildDir() {

  static const char BUILD_DIR[] = "/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Exec/Test_Only_Axions";
  return BUILD_DIR;
}

const char* buildInfoGetBuildMachine() {

  static const char BUILD_MACHINE[] = "Linux gwdu101 3.10.0-693.11.6.el7.x86_64 #1 SMP Wed Jan 3 18:09:42 CST 2018 x86_64 x86_64 x86_64 GNU/Linux";
  return BUILD_MACHINE;
}

const char* buildInfoGetAMReXDir() {

  static const char AMREX_DIR[] = "../../../amrex";
  return AMREX_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "intel";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "17.0.2";
  return COMP_VERSION;
}

// deprecated
const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "";
  return FCOMP;
}

// deprecated
const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "";
  return FCOMP_VERSION;
}

const char* buildInfoGetCXXName() {

  static const char CXX_comp_name[] = "mpicxx";
  return CXX_comp_name;
}

const char* buildInfoGetFName() {

  static const char F_comp_name[] = "mpif90";
  return F_comp_name;
}

const char* buildInfoGetCXXFlags() {

  static const char CXX_flags[] = " -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec -std=c++14  -qopenmp -pthread -DNDEBUG -DBL_USE_MPI -DAMREX_USE_MPI -DBL_USE_OMP -DAMREX_USE_OMP -DAMREX_GIT_VERSION=\"19.07-153-g69a0b11aef66-dirty\" -DAMREX_LAUNCH= -DAMREX_DEVICE= -DAMREX_CUDA_FORT_GLOBAL= -DAMREX_CUDA_FORT_DEVICE= -DAMREX_CUDA_FORT_HOST= -DBL_SPACEDIM=3 -DAMREX_SPACEDIM=3 -DBL_FORT_USE_UNDERSCORE -DAMREX_FORT_USE_UNDERSCORE -DBL_Linux -DAMREX_Linux -DAMREX_PARTICLES -DCRSEGRNDOMP -DGRAVITY -DCGRAV -DNO_HYDRO -DFDM -DBL_NOLINEVALUES -I. -I/usr/users/cbehren2/local/include -I../../../amrex/Src/Extern/SWFFT -I/usr/users/cbehren2/local_libs/lib/../include -I. -I../../Source -I../../Source/Src_3d -I../../Source/HydroFortran -I../../Source/Tagging -I../../Source/Initialization -I../../Source/EOS -I../../Source/Network -I../../Source/HeatCool -I../../Source/SourceTerms -I../../Source/DerivedQuantities -I../../Source/Monitors -I../../Source/Gravity -I../../Source/FDM -I../../Source/AMReX_axionyx -I../../../amrex/Src/Base -I../../../amrex/Src/AmrCore -I../../../amrex/Src/Amr -I../../../amrex/Src/Boundary -I../../../amrex/Src/Particle -I../../../amrex/Src/Extern/amrdata -I../../../amrex/Src/Base -I../../../amrex/Src/AmrCore -I../../../amrex/Src/Amr -I../../../amrex/Src/Boundary -I../../../amrex/Src/Particle -I../../../amrex/Src/Extern/amrdata -I../../Source/Constants -I../../../amrex/Src/LinearSolvers/C_CellMG -I../../../amrex/Src/LinearSolvers/C_CellMG -I../../../amrex/Src/LinearSolvers/MLMG -I../../../amrex/Src/LinearSolvers/MLMG -I../../Util/VODE -I../../Util/BLAS -I../../../amrex/Tools/C_scripts";
  return CXX_flags;
}

const char* buildInfoGetFFlags() {

  static const char F_flags[] = " -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec -implicitnone  -qopenmp";
  return F_flags;
}

const char* buildInfoGetLinkFlags() {

  static const char link_flags[] = "-L. ";
  return link_flags;
}

const char* buildInfoGetLibraries() {

  static const char libraries[] = "-L/usr/users/cbehren2/local/lib -lfftw3_mpi -lfftw3 -lfftw3_omp -L/usr/users/cbehren2/local_libs/lib  -lifcore -lifcoremt";
  return libraries;
}

const char* buildInfoGetAux(int i) {

  //static const char AUX1[] = "${AUX[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

int buildInfoGetNumModules() {
  // int num_modules = X;
  int num_modules = 0;

  return num_modules;
}

const char* buildInfoGetModuleName(int i) {

  //static const char MNAME1[] = "${MNAME[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetModuleVal(int i) {

  //static const char MVAL1[] = "${MVAL[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetGitHash(int i) {

  //static const char HASH1[] = "${GIT[1]}";
  static const char HASH1[] = "18.05-325-g1014029-dirty";
  static const char HASH2[] = "19.07-153-g69a0b11-dirty";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return HASH1;
    case 2: return HASH2;

    default: return EMPT;
  }
}

const char* buildInfoGetBuildGitHash() {

  //static const char HASH[] = "${GIT}";
  static const char HASH[] = "";


  return HASH;
}

const char* buildInfoGetBuildGitName() {

  //static const char NAME[] = "";
  static const char NAME[] = "";


  return NAME;
}

}
````

