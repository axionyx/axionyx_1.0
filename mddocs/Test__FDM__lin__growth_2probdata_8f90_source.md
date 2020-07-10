
# File probdata.f90

[**File List**](files.md) **>** [**Exec**](dir_43a12cefb7942b6f49b5b628aafd3192.md) **>** [**Test\_FDM\_lin\_growth**](dir_97d68d96f71fb742273f5b8113c9a269.md) **>** [**probdata.f90**](Test__FDM__lin__growth_2probdata_8f90.md)

[Go to the documentation of this file.](Test__FDM__lin__growth_2probdata_8f90.md) 


````cpp
module probdata_module

!     Tagging variables
      integer, save :: max_num_part

!     Comoving variables
!      double precision, save :: comoving_OmAx   (now in comoving_module)

!     Residual variables
      double precision, save :: center(3)

!     Needed for fort_prescribe_grav 
      double precision, save :: dcenx,dceny,dcenz,dmconc,dmmass,dmscale
!     
      
end module probdata_module
````

