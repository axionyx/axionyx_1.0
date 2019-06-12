
# File SWFFT\_Test\_F.H

[**File List**](files.md) **>** [**Exec**](dir_43a12cefb7942b6f49b5b628aafd3192.md) **>** [**Test\_Only\_Axions**](dir_eb24725df855cf6c732a19e4912f662a.md) **>** [**SWFFT\_Test\_F.H**](SWFFT__Test__F_8H.md)

[Go to the documentation of this file.](SWFFT__Test__F_8H.md) 


````cpp
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
````

