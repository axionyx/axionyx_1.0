subroutine deposit_fdm_particles(particles, np, ng, state, &!, ghosts, ng, virts, nv
     lo, hi, plo, dx) &
     bind(c,name='deposit_fdm_particles')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use particle_mod      , only: fdm_particle_t
  use fundamental_constants_module, only: pi
  use axion_params_module, only : hbaroverm, theta_ax, sigma_ax, ii
  use meth_params_module, only : NAXVAR, UAXDENS, UAXRE, UAXIM

  integer                           :: np, ng, nv
  integer                           :: lo(3), hi(3)
  type(fdm_particle_t), intent(in ) :: particles(np)
  ! type(fdm_particle_t), intent(in ) :: ghosts(ng)
  ! type(fdm_particle_t), intent(in ) :: virts(nv)
  real(amrex_real),     intent(out) :: state(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),NAXVAR)
  real(amrex_real)                  :: plo(3), dx(3)
  
  integer i, j, k, n, i1, j1, k1, rad
  real(amrex_real) pos(3), kernelsize, limit, inv_dx(3)
  complex(amrex_real) :: A, phi!, ii

  inv_dx = 1.0d0/dx
  limit = hbaroverm*theta_ax*theta_ax

  do n = 1, np

     A   = cmplx(particles(n)%amp(1),particles(n)%amp(2))
     rad = ceiling(theta_ax/sqrt(2.0*particles(n)%width)*inv_dx(1))
     pos = (particles(n)%pos - plo)*inv_dx + 0.5d0! - real(ng)
     
     ! print *, A, particles(n)%phase, inv_dx(1)

     i1  = floor(pos(1))
     j1  = floor(pos(2))
     k1  = floor(pos(3))
     
     do k=-rad,rad,1
        if ( (k1+k).ge.lo(3) .and. (k1+k).le.hi(3) ) then
           do j=-rad,rad,1
              if ( (j1+j).ge.lo(2) .and. (j1+j).le.hi(2) ) then
                 do i=-rad,rad,1
                    if ( (i1+i).ge.lo(1) .and. (i1+i).le.hi(1) ) then
                       
                       kernelsize = ((real(i1+i)-pos(1))*dx(1)*(real(i1+i)-pos(1))*dx(1)+ &
                            (real(j1+j)-pos(2))*dx(2)*(real(j1+j)-pos(2))*dx(2)+ &
                            (real(k1+k)-pos(3))*dx(3)*(real(k1+k)-pos(3))*dx(3))* &
                            particles(n)%width
                       
                       if (kernelsize .le. (theta_ax*theta_ax/2.0)) then
                          
                          phi = A*exp(-kernelsize)*exp(ii*(particles(n)%phase+ &
                               particles(n)%vel(1)*(real(i1+i)-0.5-pos(1))*dx(1)+ &
                               particles(n)%vel(2)*(real(j1+j)-0.5-pos(2))*dx(2)+ &
                               particles(n)%vel(3)*(real(k1+k)-0.5-pos(3))*dx(3) )/hbaroverm)
                          
                          state(i1+i,j1+j,k1+k,UAXRE) = state(i1+i,j1+j,k1+k,UAXRE) + real(real(phi))
                          state(i1+i,j1+j,k1+k,UAXIM) = state(i1+i,j1+j,k1+k,UAXIM) + real(aimag(phi))

                       endif
                       
                    endif
                 enddo
              endif
           enddo
        endif
     enddo
     
  enddo

!  print *, state(:,:,:,2)

  ! do n = 1, ng

  !    A   = cmplx(ghosts(n)%amp(1),ghosts(n)%amp(2))
  !    rad = ceiling(theta_ax/sqrt(2.0*ghosts(n)%width)*inv_dx(1))
  !    pos = ghosts(n)%pos - plo)*inv_dx + 0.5d0
     
  !    i1  = floor(pos(1))
  !    j1  = floor(pos(2))
  !    k1  = floor(pos(3))
     
  !    do k=-rad,rad,1
  !       if ( (k1+k).ge.1 .and. (k1+k).le.dim3 ) then
  !          do j=-rad,rad,1
  !             if ( (j1+j).ge.1 .and. (j1+j).le.dim2 ) then
  !                do i=-rad,rad,1
  !                   if ( (i1+i).ge.1 .and. (i1+i).le.dim1 ) then
                       
  !                      kernelsize = ((real(i1+i)-0.5-pos(1))*(real(i1+i)-0.5-pos(1))+ &
  !                           (real(j1+j)-0.5-pos(2))*(real(j1+j)-0.5-pos(2))+ &
  !                           (real(k1+k)-0.5-pos(3))*(real(k1+k)-0.5-pos(3)))/ &
  !                           sigma_ax/sigma_ax
                       
  !                      if (kernelsize .le. (theta_ax*theta_ax)) then
                          
  !                         phi = amp*exp(-kernelsize/2.0)*exp(ii*(ghosts(n)%phase+ &
  !                              ghosts(n)%vel(1)*(real(i1+i)-0.5-pos(1))*dx(1)+ &
  !                              ghosts(n)%vel(2)*(real(j1+j)-0.5-pos(2))*dx(2)+ &
  !                              ghosts(n)%vel(3)*(real(k1+k)-0.5-pos(3))*dx(3) )/hbaroverm)
                          
  !                         state(i1+i,j1+j,k1+k,1) = state(i1+i,j1+j,k1+k,1) + real(real(phi))
  !                         state(i1+i,j1+j,k1+k,2) = state(i1+i,j1+j,k1+k,2) + real(aimag(phi))
                          
  !                      endif
                       
  !                   endif
  !                enddo
  !             endif
  !          enddo
  !       endif
  !    enddo
     
  ! enddo

  ! do n = 1, nv

  !    A   = cmplx(virts(n)%amp(1),virts(n)%amp(2))
  !    rad = ceiling(theta_ax/sqrt(2.0*virts(n)%width)*inv_dx(1))
  !    pos = virts(n)%pos - plo)*inv_dx + 0.5d0
     
  !    i1  = floor(pos(1))
  !    j1  = floor(pos(2))
  !    k1  = floor(pos(3))
     
  !    do k=-rad,rad,1
  !       if ( (k1+k).ge.1 .and. (k1+k).le.dim3 ) then
  !          do j=-rad,rad,1
  !             if ( (j1+j).ge.1 .and. (j1+j).le.dim2 ) then
  !                do i=-rad,rad,1
  !                   if ( (i1+i).ge.1 .and. (i1+i).le.dim1 ) then
                       
  !                      kernelsize = ((real(i1+i)-0.5-pos(1))*(real(i1+i)-0.5-pos(1))+ &
  !                           (real(j1+j)-0.5-pos(2))*(real(j1+j)-0.5-pos(2))+ &
  !                           (real(k1+k)-0.5-pos(3))*(real(k1+k)-0.5-pos(3)))/ &
  !                           sigma_ax/sigma_ax
                       
  !                      if (kernelsize .le. (theta_ax*theta_ax)) then
                          
  !                         phi = amp*exp(-kernelsize/2.0)*exp(ii*(virts(n)%phase+ &
  !                              virts(n)%vel(1)*(real(i1+i)-0.5-pos(1))*dx(1)+ &
  !                              virts(n)%vel(2)*(real(j1+j)-0.5-pos(2))*dx(2)+ &
  !                              virts(n)%vel(3)*(real(k1+k)-0.5-pos(3))*dx(3) )/hbaroverm)
                          
  !                         state(i1+i,j1+j,k1+k,1) = state(i1+i,j1+j,k1+k,1) + real(real(phi))
  !                         state(i1+i,j1+j,k1+k,2) = state(i1+i,j1+j,k1+k,2) + real(aimag(phi))
                          
  !                      endif
                       
  !                   endif
  !                enddo
  !             endif
  !          enddo
  !       endif
  !    enddo
     
  ! enddo
  
  ! state(:,:,:,3) = state(:,:,:,1)**2+state(:,:,:,2)**2
  
end subroutine deposit_fdm_particles

subroutine fort_fdm_fields( &
     uout, uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3)
        
  use meth_params_module, only : NAXVAR, UAXDENS, UAXRE, UAXIM
  use fundamental_constants_module

  implicit none
  
  integer          uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
  double precision uout( uout_l1:uout_h1, uout_l2:uout_h2, uout_l3:uout_h3, NAXVAR)
  
  uout(:,:,:,UAXDENS) = uout(:,:,:,UAXRE)**2+uout(:,:,:,UAXIM)**2

  ! print *, uout(:,:,:,UAXRE)

end subroutine fort_fdm_fields
