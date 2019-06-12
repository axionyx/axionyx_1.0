
# File set\_dirichlet\_bcs\_3d.f90

[**File List**](files.md) **>** [**Gravity**](dir_fdbf5007869eac89a42b1cd44aeda050.md) **>** [**set\_dirichlet\_bcs\_3d.f90**](set__dirichlet__bcs__3d_8f90.md)

[Go to the documentation of this file.](set__dirichlet__bcs__3d_8f90.md) 


````cpp

     subroutine fort_set_homog_bcs(lo,hi,domlo,domhi, &
                                   phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,dx);
 
     use amrex_fort_module, only : rt => amrex_real
     implicit none

     integer         ,intent(in   ) :: lo(3),hi(3),domlo(3),domhi(3)
     integer         ,intent(in   ) :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
     real(rt),intent(  out) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
     real(rt),intent(in   ) :: dx(3)

     phi = 0.d0
 
     end subroutine fort_set_homog_bcs

````

