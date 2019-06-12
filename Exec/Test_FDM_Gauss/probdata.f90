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
