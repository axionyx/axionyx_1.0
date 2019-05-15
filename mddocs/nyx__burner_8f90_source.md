
# File nyx\_burner.f90

[**File List**](files.md) **>** [**Network**](dir_42bb2cb79beb2277fb25f45fdc565a0d.md) **>** [**nyx\_burner.f90**](nyx__burner_8f90.md)

[Go to the documentation of this file.](nyx__burner_8f90.md) 


````cpp
module nyx_burner_module

  use amrex_fort_module, only : rt => amrex_real
  use eos_module
  use network

contains

  subroutine burner(dens, temp, Xin, ein, dt, time, Xout, eout)

    implicit none
    
    real(rt), intent(in) :: dens, temp, Xin(nspec), ein, dt, time
    real(rt), intent(out) :: Xout(nspec), eout
    
    integer  :: n
    real(rt) :: enuc, dX
    
    xout(:) = xin(:)
    
    enuc = 0.0_rt
    do n = 1, nspec
       dx = xout(n)-xin(n) 
       enuc = enuc - ebin(n) * dx
    enddo
  
    eout = ein + enuc
  
  end subroutine burner

end module nyx_burner_module
````

