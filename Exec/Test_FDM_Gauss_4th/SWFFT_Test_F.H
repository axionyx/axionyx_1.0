#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

  void fort_ax_fields (amrex_real* axion, const int* lo, const int* hi);

  void fort_kick (const int* lo, const int* hi, amrex_real* real, amrex_real* imag, amrex_real* phi,
		  const amrex_real* hbaroverm, const amrex_real* dt);

#ifdef __cplusplus
}
#endif
