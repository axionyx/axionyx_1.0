
# File Gravity\_nd.f90

[**File List**](files.md) **>** [**Gravity**](dir_fdbf5007869eac89a42b1cd44aeda050.md) **>** [**Gravity\_nd.f90**](Gravity__nd_8f90.md)

[Go to the documentation of this file.](Gravity__nd_8f90.md) 


````cpp

      subroutine fort_get_grav_const(Gconst_out) bind(C,name="fort_get_grav_const")

         use amrex_fort_module, only : rt => amrex_real
         use fundamental_constants_module, only: gconst

         real(rt) :: Gconst_out

         gconst_out = gconst

      end subroutine fort_get_grav_const

````

