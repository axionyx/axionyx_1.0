subroutine deposit_fdm_particles(particles, np, state_real, &
     lo_real, hi_real, state_imag, lo_imag, hi_imag, plo, dx, a) &
     bind(c,name='deposit_fdm_particles')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use particle_mod      , only: fdm_particle_t
  use axion_params_module, only : hbaroverm, theta_ax, ii

  integer,              intent(in ) :: np, lo_real(3), hi_real(3), lo_imag(3), hi_imag(3)
  type(fdm_particle_t), intent(in ) :: particles(np)
  real(amrex_real),     intent(out) :: state_real(lo_real(1):hi_real(1), lo_real(2):hi_real(2), lo_imag(3):hi_imag(3),1)
  real(amrex_real),     intent(out) :: state_imag(lo_imag(1):hi_imag(1), lo_imag(2):hi_imag(2), lo_imag(3):hi_imag(3),1)
  real(amrex_real),     intent(in ) :: plo(3), dx(3), a
  
  integer i, j, k, n, i1, j1, k1, rad
  real(amrex_real) pos(3), kernelsize
  complex(amrex_real) :: amp, phi

  do n = 1, np

     amp = cmplx(particles(n)%amp(1),particles(n)%amp(2))
     rad = ceiling(theta_ax/sqrt(2.0*particles(n)%width)/dx(1))
     pos = (particles(n)%pos - plo)/dx + 0.5d0

     i1  = floor(pos(1))
     j1  = floor(pos(2))
     k1  = floor(pos(3))

     do k=-rad,rad,1
        if ( (k1+k).ge.lo_real(3) .and. (k1+k).le.hi_real(3) ) then
           do j=-rad,rad,1
              if ( (j1+j).ge.lo_real(2) .and. (j1+j).le.hi_real(2) ) then
                 do i=-rad,rad,1
                    if ( (i1+i).ge.lo_real(1) .and. (i1+i).le.hi_real(1) ) then
                       
                       kernelsize = ((real(i1+i)+1.0-pos(1))*dx(1)*(real(i1+i)+1.0-pos(1))*dx(1)+ &
                            (real(j1+j)+1.0-pos(2))*dx(2)*(real(j1+j)+1.0-pos(2))*dx(2)+ &
                            (real(k1+k)+1.0-pos(3))*dx(3)*(real(k1+k)+1.0-pos(3))*dx(3))* &
                            particles(n)%width
                       
                       if (kernelsize .le. (theta_ax*theta_ax/2.0)) then
                          
                          phi = amp*exp(-kernelsize)*exp(ii*(particles(n)%phase+ &
                               particles(n)%vel(1)*a*(real(i1+i)-1.0-pos(1))*dx(1)+ &
                               particles(n)%vel(2)*a*(real(j1+j)-1.0-pos(2))*dx(2)+ &
                               particles(n)%vel(3)*a*(real(k1+k)-1.0-pos(3))*dx(3) )/hbaroverm)
                          
                          state_real(i1+i,j1+j,k1+k,1) = state_real(i1+i,j1+j,k1+k,1) + real(real(phi))
                          state_imag(i1+i,j1+j,k1+k,1) = state_imag(i1+i,j1+j,k1+k,1) + real(aimag(phi))

                       endif
                       
                    endif
                 enddo
              endif
           enddo
        endif
     enddo
     
  enddo

end subroutine deposit_fdm_particles

subroutine deposit_fdm_particles_wkb(particles, np, state_real, &
     lo_real, hi_real, state_imag, lo_imag, hi_imag, plo, dx, a) &
     bind(c,name='deposit_fdm_particles_wkb')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use particle_mod      , only: fdm_particle_wkb_t
  use axion_params_module, only : hbaroverm, theta_ax, ii

  integer,              intent(in ) :: np, lo_real(3), hi_real(3), lo_imag(3), hi_imag(3)
  type(fdm_particle_wkb_t), intent(in ) :: particles(np)
  real(amrex_real),     intent(out) :: state_real(lo_real(1):hi_real(1), lo_real(2):hi_real(2), lo_imag(3):hi_imag(3),1)
  real(amrex_real),     intent(out) :: state_imag(lo_imag(1):hi_imag(1), lo_imag(2):hi_imag(2), lo_imag(3):hi_imag(3),1)
  real(amrex_real),     intent(in ) :: plo(3), dx(3), a
  
  integer i, j, k, n, i1, j1, k1, rad
  real(amrex_real) pos(3), kernelsize
  complex(amrex_real) :: amp, phi

  do n = 1, np

     amp = cmplx(particles(n)%amp(1),particles(n)%amp(2))
     rad = ceiling(theta_ax/sqrt(2.0*particles(n)%width)/dx(1))
     pos = (particles(n)%pos - plo)/dx + 0.5d0

     i1  = floor(pos(1))
     j1  = floor(pos(2))
     k1  = floor(pos(3))

     do k=-rad,rad,1
        if ( (k1+k).ge.lo_real(3) .and. (k1+k).le.hi_real(3) ) then
           do j=-rad,rad,1
              if ( (j1+j).ge.lo_real(2) .and. (j1+j).le.hi_real(2) ) then
                 do i=-rad,rad,1
                    if ( (i1+i).ge.lo_real(1) .and. (i1+i).le.hi_real(1) ) then
                       
                       kernelsize = ((real(i1+i)+1.0-pos(1))*dx(1)*(real(i1+i)+1.0-pos(1))*dx(1)+ &
                            (real(j1+j)+1.0-pos(2))*dx(2)*(real(j1+j)+1.0-pos(2))*dx(2)+ &
                            (real(k1+k)+1.0-pos(3))*dx(3)*(real(k1+k)+1.0-pos(3))*dx(3))* &
                            particles(n)%width
                       
                       if (kernelsize .le. (theta_ax*theta_ax/2.0)) then
                          
                          phi = amp*exp(-kernelsize)*exp(ii*(particles(n)%phase+ &
                               particles(n)%vel(1)*a*(real(i1+i)-1.0-pos(1))*dx(1)+ &
                               particles(n)%vel(2)*a*(real(j1+j)-1.0-pos(2))*dx(2)+ &
                               particles(n)%vel(3)*a*(real(k1+k)-1.0-pos(3))*dx(3) )/hbaroverm)
                          
                          state_real(i1+i,j1+j,k1+k,1) = state_real(i1+i,j1+j,k1+k,1) + real(real(phi))
                          state_imag(i1+i,j1+j,k1+k,1) = state_imag(i1+i,j1+j,k1+k,1) + real(aimag(phi))

                       endif
                       
                    endif
                 enddo
              endif
           enddo
        endif
     enddo
     
  enddo

end subroutine deposit_fdm_particles_wkb

subroutine fort_fdm_fields(state, state_l1,state_l2,state_l3,state_h1,state_h2,state_h3)
        
  use meth_params_module, only : NAXVAR, UAXDENS, UAXRE, UAXIM

  implicit none
  
  integer          state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision state( state_l1:state_h1, state_l2:state_h2, state_l3:state_h3, NAXVAR)
  
  state(:,:,:,UAXDENS) = state(:,:,:,UAXRE)**2+state(:,:,:,UAXIM)**2

end subroutine fort_fdm_fields
