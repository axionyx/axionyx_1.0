
# File ct\_upwind.f90

[**File List**](files.md) **>** [**MHD**](dir_a5db59ee0cc93a408ad0433ba32613c6.md) **>** [**ct\_upwind.f90**](ct__upwind_8f90.md)

[Go to the documentation of this file.](ct__upwind_8f90.md) 


````cpp
module ct_upwind

 use amrex_fort_module, only : rt => amrex_real
 use hlld_solver, only : hlld
! use hll_solver, only : hll
 use meth_params_module

 implicit none 

 private primtocons
 public corner_transport
 
 
interface checkisnan
       module procedure checkisnanmult
       module procedure checkisnans
end interface checkisnan

 contains

subroutine corner_transport( q, qm, qp, q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3, &  
                  flxx, flxx_l1 , flxx_l2 , flxx_l3 , flxx_h1 , flxx_h2 , flxx_h3, &
                  flxy, flxy_l1 , flxy_l2 , flxy_l3 , flxy_h1 , flxy_h2 , flxy_h3, &
                  flxz, flxz_l1 , flxz_l2 , flxz_l3 , flxz_h1 , flxz_h2 , flxz_h3, &
                  Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                  Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                  Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                  dx, dy, dz, dt)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module, only : qvar
 use electric_field
implicit none

  integer, intent(in)   :: q_l1,q_l2,q_l3,q_h1,q_h2,q_h3

  integer, intent(in)   :: ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3
  integer, intent(in)   :: ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3
  integer, intent(in)   :: ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3

  integer, intent(in)   :: flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3
  integer, intent(in)   :: flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3
  integer, intent(in)   :: flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3

    real(rt), intent(in)  :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR) !prim vars at time t^n
    real(rt), intent(in)  :: qm(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
    real(rt), intent(in)  :: qp(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
    real(rt), intent(out) :: flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,QVAR)!Half Step Fluxes
    real(rt), intent(out) :: flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,QVAR)!Half Step Fluxes
    real(rt), intent(out) :: flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,QVAR)!Half Step Fluxes

    real(rt), intent(out)  :: Ex(ex_l1:ex_h1,ex_l2:ex_h2,ex_l3:ex_h3)
    real(rt), intent(out)  :: Ey(ey_l1:ey_h1,ey_l2:ey_h2,ey_l3:ey_h3)
    real(rt), intent(out)  :: Ez(ez_l1:ez_h1,ez_l2:ez_h2,ez_l3:ez_h3)
 
    real(rt)  :: um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3) !PtoC Vars
    real(rt)  :: up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
    real(rt)  :: cons_temp_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2) !2D Temporary Conservative Vars
    real(rt)  :: cons_temp_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2) !2D Temporary Conservative Vars
    real(rt)  :: cons_half_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3) !Flux Corrected Conservative Vars
    real(rt)  :: cons_half_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
    real(rt)  :: q_temp_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2) !2D Temporary Primitive Vars
    real(rt)  :: q_temp_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2) !2D Temporary Primitive Vars
    real(rt)  :: q_half_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3) !Flux Corrected Primitive Vars
    real(rt)  :: q_half_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3) !Flux Corrected Primitive Vars

    real(rt)  :: flxx1D(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,QVAR)
    real(rt)  :: flxy1D(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,QVAR) 
  real(rt)  :: flxz1D(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,QVAR) !Flux1d for all directions

    real(rt)  :: flxx2D(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,QVAR, 2) !Flux2d for all directions 2 perpendicular directions
    real(rt)  :: flxy2D(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,QVAR, 2) !Flux2d for all directions 2 perpendicular directions
    real(rt)  :: flxz2D(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,QVAR, 2) !Flux2d for all directions 2 perpendicular directions

    real(rt)  :: q2D(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
    real(rt)  :: dx, dy, dz, dt

  integer   :: i, work_lo(3), work_hi(3)


  um = 0.d0
  up = 0.d0
  cons_temp_m = 0.d0
  cons_temp_p = 0.d0
  q_temp_m = 0.d0
  q_temp_p = 0.d0
  cons_half_m = 0.d0
  cons_half_p = 0.d0
  q_half_m = 0.d0
  q_half_p = 0.d0
  q2d = 0.d0

  flxx1d = 0.d0
  flxy1d = 0.d0
  flxz1d = 0.d0

  flxx2d = 0.d0
  flxy2d = 0.d0
  flxz2d = 0.d0

  ex = 0.d0
  ey = 0.d0
  ez = 0.d0

!Calculate Flux 1D
!x-dir
  work_lo(1) = flxx_l1
  work_lo(2) = flxx_l2
  work_lo(3) = flxx_l3
  work_hi(1) = flxx_h1
  work_hi(2) = flxx_h2
  work_hi(3) = flxx_h3
        
  call hlld(work_lo, work_hi, qm,qp,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
            flxx1d(:,:,:,:),flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, 1)

  work_lo(1) = flxy_l1
  work_lo(2) = flxy_l2
  work_lo(3) = flxy_l3
  work_hi(1) = flxy_h1
  work_hi(2) = flxy_h2
  work_hi(3) = flxy_h3

!y-dir  
  call hlld(work_lo, work_hi, qm,qp,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
            flxy1d(:,:,:,:),flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, 2)


  work_lo(1) = flxz_l1
  work_lo(2) = flxz_l2
  work_lo(3) = flxz_l3
  work_hi(1) = flxz_h1
  work_hi(2) = flxz_h2
  work_hi(3) = flxz_h3
                  
!z-dir
  call hlld(work_lo, work_hi, qm,qp,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
            flxz1d(:,:,:,:),flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, 3)

!Prim to Cons
do i = 1,3
  call primtocons(qm(:,:,:,:,i), um(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
  call primtocons(qp(:,:,:,:,i), up(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
enddo


!Use "1D" fluxes To interpolate Temporary Edge Centered Electric Fields


  work_lo(1) = ex_l1+1
  work_lo(2) = ex_l2+1
  work_lo(3) = ex_l3+1
  work_hi(1) = ex_h1-1
  work_hi(2) = ex_h2-1
  work_hi(3) = ex_h3-1
  call electric_edge_x(work_lo, work_hi, &
           q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
           ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
           flxy1d, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
           flxz1d, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

  work_lo(1) = ey_l1+1
  work_lo(2) = ey_l2+1
  work_lo(3) = ey_l3+1
  work_hi(1) = ey_h1-1
  work_hi(2) = ey_h2-1
  work_hi(3) = ey_h3-1
  call electric_edge_y(work_lo, work_hi, &
           q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
           ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
           flxx1d, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
           flxz1d, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

  work_lo(1) = ez_l1+1
  work_lo(2) = ez_l2+1
  work_lo(3) = ez_l3+1
  work_hi(1) = ez_h1-1
  work_hi(2) = ez_h2-1
  work_hi(3) = ez_h3-1
  call electric_edge_z(work_lo, work_hi, &
           q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
           ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
           flxx1d, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
           flxy1d, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3)

!Corner Couple
  work_lo(1) = q_l1+2
  work_lo(2) = q_l2+2
  work_lo(3) = q_l3+2
  work_hi(1) = q_h1-2
  work_hi(2) = q_h2-2
  work_hi(3) = q_h3-2

  call corner_couple(work_lo, work_hi, &
       cons_temp_m, cons_temp_p, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
       flxx1d, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
       flxy1d, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
       flxz1d, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
       dx, dy, dz, dt) !Correct Conservative vars using Transverse Fluxes

    call corner_couple_mag(work_lo, work_hi, &
         cons_temp_m, cons_temp_p, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
         ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
         ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
         ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
         dx, dy, dz, dt) !Correct Conservative vars using Transverse Fluxes


!Cons To Prim
do i = 1,3
  call constoprim(q_temp_m(:,:,:,:,i,1), cons_temp_m(:,:,:,:,i,1), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
  call constoprim(q_temp_p(:,:,:,:,i,1), cons_temp_p(:,:,:,:,i,1), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
  call constoprim(q_temp_m(:,:,:,:,i,2), cons_temp_m(:,:,:,:,i,2), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
  call constoprim(q_temp_p(:,:,:,:,i,2), cons_temp_p(:,:,:,:,i,2), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
enddo

!Calculate Flux 2D
do i = 1,2
  work_lo(1) = flxx_l1+1
  work_lo(2) = flxx_l2+1
  work_lo(3) = flxx_l3+1
  work_hi(1) = flxx_h1-1
  work_hi(2) = flxx_h2-1
  work_hi(3) = flxx_h3-1
!x-dir
  call hlld(work_lo, work_hi, q_temp_m(:,:,:,:,:,i),q_temp_p(:,:,:,:,:,i),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
            flxx2d(:,:,:,:,i),flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, 1)

!y-dir  
  work_lo(1) = flxy_l1+1
  work_lo(2) = flxy_l2+1
  work_lo(3) = flxy_l3+1
  work_hi(1) = flxy_h1-1
  work_hi(2) = flxy_h2-1
  work_hi(3) = flxy_h3-1

  call hlld(work_lo, work_hi, q_temp_m(:,:,:,:,:,i),q_temp_p(:,:,:,:,:,i),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
            flxy2d(:,:,:,:,i),flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, 2)

!z-dir
  work_lo(1) = flxz_l1+1
  work_lo(2) = flxz_l2+1
  work_lo(3) = flxz_l3+1
  work_hi(1) = flxz_h1-1
  work_hi(2) = flxz_h2-1
  work_hi(3) = flxz_h3-1

  call hlld(work_lo, work_hi, q_temp_m(:,:,:,:,:,i),q_temp_p(:,:,:,:,:,i),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
            flxz2d(:,:,:,:,i),flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, 3)
enddo

!Use Averaged 2D fluxes to interpolate temporary Edge Centered Electric Fields, reuse "flx1D"
  flxx1d(:,:,:,:) = 0.5d0*(flxx2d(:,:,:,:,1) + flxx2d(:,:,:,:,2))
  flxy1d(:,:,:,:) = 0.5d0*(flxy2d(:,:,:,:,1) + flxy2d(:,:,:,:,2))
  flxz1d(:,:,:,:) = 0.5d0*(flxz2d(:,:,:,:,1) + flxz2d(:,:,:,:,2))


  work_lo(1) = ex_l1+1!+2
  work_lo(2) = ex_l2+1!+2
  work_lo(3) = ex_l3+1!+2
  work_hi(1) = ex_h1-1!-2
  work_hi(2) = ex_h2-1!-2
  work_hi(3) = ex_h3-1!-2
  call electric_edge_x(work_lo, work_hi, &
              q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
              ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
              flxy1d, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
              flxz1d, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

  work_lo(1) = ey_l1+1!+2
  work_lo(2) = ey_l2+1!+2
  work_lo(3) = ey_l3+1!+2
  work_hi(1) = ey_h1-1!-2
  work_hi(2) = ey_h2-1!-2
  work_hi(3) = ey_h3-1!-2
  call electric_edge_y(work_lo, work_hi, &
              q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
              ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
              flxx1d, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
              flxz1d, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

  work_lo(1) = ez_l1+1!+2
  work_lo(2) = ez_l2+1!+2
  work_lo(3) = ez_l3+1!+2
  work_hi(1) = ez_h1-1!-2
  work_hi(2) = ez_h2-1!-2
  work_hi(3) = ez_h3-1!-2
  call electric_edge_z(work_lo, work_hi, &
              q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
              ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
              flxx1d, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
              flxy1d, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3)

!Half Step conservative vars
  work_lo(1) = q_l1+1
  work_lo(2) = q_l2+1
  work_lo(3) = q_l3+1
  work_hi(1) = q_h1-2
  work_hi(2) = q_h2-2
  work_hi(3) = q_h3-2

  call half_step(work_lo, work_hi, &
           cons_half_m, cons_half_p, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
           flxx2d, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
           flxy2d, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
           flxz2d, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
           dx, dy, dz, dt)

  call half_step_mag(work_lo, work_hi, &
            cons_half_m, cons_half_p, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
            ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
            ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
            ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
            dx, dy, dz, dt)
do i = 1,3
  call constoprim(q_half_m(:,:,:,:,i), cons_half_m(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
  call constoprim(q_half_p(:,:,:,:,i), cons_half_p(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
enddo

!Final Fluxes


!x-dir
  work_lo(1) = flxx_l1+2
  work_lo(2) = flxx_l2+2
  work_lo(3) = flxx_l3+2
  work_hi(1) = flxx_h1-2
  work_hi(2) = flxx_h2-2
  work_hi(3) = flxx_h3-2

  call hlld(work_lo, work_hi, q_half_m, q_half_p,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
            flxx,flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, 1)

!y-dir  
  work_lo(1) = flxy_l1+2
  work_lo(2) = flxy_l2+2
  work_lo(3) = flxy_l3+2
  work_hi(1) = flxy_h1-2
  work_hi(2) = flxy_h2-2
  work_hi(3) = flxy_h3-2

  call hlld(work_lo, work_hi, q_half_m,q_half_p,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
            flxy,flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, 2)


!z-dir
  work_lo(1) = flxz_l1+2
  work_lo(2) = flxz_l2+2
  work_lo(3) = flxz_l3+2
  work_hi(1) = flxz_h1-2
  work_hi(2) = flxz_h2-2
  work_hi(3) = flxz_h3-2

  call hlld(work_lo, work_hi, q_half_m,q_half_p,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
            flxz,flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, 3)

!Primitive update
  call prim_half(q2d,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
          flxx1d, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
          flxy1d, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
          flxz1d, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
          dx, dy, dz, dt)

!Final Electric Field Update
  work_lo(1) = ex_l1+2!+3
  work_lo(2) = ex_l2+2!+3
  work_lo(3) = ex_l3+2!+3
  work_hi(1) = ex_h1-2!-3
  work_hi(2) = ex_h2-2!-3
  work_hi(3) = ex_h3-2!-3
  call electric_edge_x(work_lo, work_hi, &
            q2d, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
            ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
            flxy, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
            flxz, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

  work_lo(1) = ey_l1+2!+3
  work_lo(2) = ey_l2+2!+3
  work_lo(3) = ey_l3+2!+3
  work_hi(1) = ey_h1-2!-3
  work_hi(2) = ey_h2-2!-3
  work_hi(3) = ey_h3-2!-3
  call electric_edge_y(work_lo, work_hi, &
            q2d, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
            ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
            flxx, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
            flxz, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

  work_lo(1) = ez_l1+2!+3
  work_lo(2) = ez_l2+2!+3
  work_lo(3) = ez_l3+2!+3
  work_hi(1) = ez_h1-2!-3
  work_hi(2) = ez_h2-2!-3
  work_hi(3) = ez_h3-2!-3
  call electric_edge_z(work_lo, work_hi, &
               q2d, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
               ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
               flxx, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
               flxy, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3)

end subroutine corner_transport

!================================================= Calculate the Conservative Variables ===============================================

subroutine primtocons(q, u, q_l1 ,q_l2 ,q_l3 ,q_h1 ,q_h2 ,q_h3)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none

  integer,  intent(in )  :: q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
    real(rt), intent(in )  :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
    real(rt), intent(out)  :: u(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  integer                :: i ,j ,k

 do k = q_l3,q_h3
   do j = q_l2,q_h2
     do i = q_l1, q_h1
      u(i,j,k,urho)  = q(i,j,k,qrho)
      u(i,j,k,umx)    = q(i,j,k,qrho)*q(i,j,k,qu)
      u(i,j,k,umy)    = q(i,j,k,qrho)*q(i,j,k,qv)
      u(i,j,k,umz)    = q(i,j,k,qrho)*q(i,j,k,qw)
      u(i,j,k,ueint) = q(i,j,k,qpres)/gamma_minus_1
      u(i,j,k,ueden) = u(i,j,k,ueint) + 0.5d0*q(i,j,k,qrho)*dot_product(q(i,j,k,qu:qw),q(i,j,k,qu:qw)) &
                     + 0.5d0*(dot_product(q(i,j,k,qmagx:qmagz),q(i,j,k,qmagx:qmagz)))
      u(i,j,k,qmagx:qmagz) = q(i,j,k,qmagx:qmagz)
     enddo
   enddo
 enddo
end subroutine primtocons

!================================================= Calculate the Primitve Variables ===============================================

subroutine constoprim(q, u, q_l1 ,q_l2 ,q_l3 ,q_h1 ,q_h2 ,q_h3)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none

  integer , intent(in )  ::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
    real(rt), intent(in )  ::u(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
    real(rt), intent(out)  ::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  integer                :: i ,j ,k

  q = u
 do k = q_l3,q_h3
   do j = q_l2,q_h2
     do i = q_l1, q_h1
      q(i,j,k,qrho)  = u(i,j,k,urho)
      q(i,j,k,qu)    = u(i,j,k,umx)/q(i,j,k,qrho)
      q(i,j,k,qv)    = u(i,j,k,umy)/q(i,j,k,qrho)
      q(i,j,k,qw)    = u(i,j,k,umz)/q(i,j,k,qrho)
      q(i,j,k,qreint) = u(i,j,k,ueden) - 0.5d0*q(i,j,k,qrho)*dot_product(q(i,j,k,qu:qw),q(i,j,k,qu:qw)) &
                                      - 0.5d0*dot_product(u(i,j,k,qmagx:qmagz), u(i,j,k,qmagx:qmagz)) !gives rho e                           
      q(i,j,k,qpres) = q(i,j,k,qreint)*gamma_minus_1! + 0.5*dot_product(u(i,j,k,QMAGX:QMAGZ),u(i,j,k,QMAGX:QMAGZ))
      q(i,j,k,qmagx:qmagz) = u(i,j,k,qmagx:qmagz)
     enddo
   enddo
 enddo
end subroutine constoprim

!======================================= Update the Temporary Conservative Variables with Transverse 1D Fluxes ========================
subroutine corner_couple(lo, hi, &
                         uL, uR, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                         flxx, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                         flxy, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                         flxz, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                         dx, dy, dz, dt)

use amrex_fort_module, only : rt => amrex_real
use meth_params_module

implicit none

  integer, intent(in )  :: lo(3), hi(3)
  integer, intent(in )  :: q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
  integer, intent(in )  :: flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3
  integer, intent(in )  :: flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3
  integer, intent(in )  :: flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3
  
    real(rt), intent(in ) ::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
    real(rt), intent(in ) ::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)

    real(rt), intent(in ) :: flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,QVAR)
    real(rt), intent(in ) :: flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,QVAR)
    real(rt), intent(in ) :: flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,QVAR)

    real(rt), intent(out) :: uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2)
    real(rt), intent(out) :: uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2)

    real(rt)        :: dx, dy, dz, dt, u, v, w
  integer         :: i ,j ,k


  ul(:,:,:,:,:,1) = um
  ul(:,:,:,:,:,2) = um
  ur(:,:,:,:,:,1) = up
  ur(:,:,:,:,:,2) = up

  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
!Left Corrected States
        ul(i,j,k,urho:ueden,1,1) = um(i,j,k,urho:ueden,1) - dt/(3.d0*dy)*(flxy(i,j+1,k,urho:ueden) - flxy(i,j,k,urho:ueden))! y corrected x
        ul(i,j,k,urho:ueden,1,2) = um(i,j,k,urho:ueden,1) - dt/(3.d0*dz)*(flxz(i,j,k+1,urho:ueden) - flxz(i,j,k,urho:ueden))! z corrected x
        ul(i,j,k,urho:ueden,2,1) = um(i,j,k,urho:ueden,2) - dt/(3.d0*dx)*(flxx(i+1,j,k,urho:ueden) - flxx(i,j,k,urho:ueden))! x corrected y
        ul(i,j,k,urho:ueden,2,2) = um(i,j,k,urho:ueden,2) - dt/(3.d0*dz)*(flxz(i,j,k+1,urho:ueden) - flxz(i,j,k,urho:ueden))! z corrected y
        ul(i,j,k,urho:ueden,3,1) = um(i,j,k,urho:ueden,3) - dt/(3.d0*dx)*(flxx(i+1,j,k,urho:ueden) - flxx(i,j,k,urho:ueden))! x corrected z
        ul(i,j,k,urho:ueden,3,2) = um(i,j,k,urho:ueden,3) - dt/(3.d0*dy)*(flxy(i,j+1,k,urho:ueden) - flxy(i,j,k,urho:ueden))! y corrected z

        u = ul(i,j,k,umx,1,1)/ul(i,j,k,urho,1,1)
        v = ul(i,j,k,umy,1,1)/ul(i,j,k,urho,1,1)
        w = ul(i,j,k,umz,1,1)/ul(i,j,k,urho,1,1)
        ul(i,j,k,ueint,1,1) = ul(i,j,k,ueden,1,1) - 0.5d0*ul(i,j,k,urho,1,1)*(u**2 + v**2 + w**2)

        u = ul(i,j,k,umx,1,2)/ul(i,j,k,urho,1,2)
        v = ul(i,j,k,umy,1,2)/ul(i,j,k,urho,1,2)
        w = ul(i,j,k,umz,1,2)/ul(i,j,k,urho,1,2)
        ul(i,j,k,ueint,1,2) = ul(i,j,k,ueden,1,2) - 0.5d0*ul(i,j,k,urho,1,2)*(u**2 + v**2 + w**2)

        u = ul(i,j,k,umx,2,1)/ul(i,j,k,urho,2,1)
        v = ul(i,j,k,umy,2,1)/ul(i,j,k,urho,2,1)
        w = ul(i,j,k,umz,2,1)/ul(i,j,k,urho,2,1)
        ul(i,j,k,ueint,2,1) = ul(i,j,k,ueden,2,1) - 0.5d0*ul(i,j,k,urho,2,1)*(u**2 + v**2 + w**2)

        u = ul(i,j,k,umx,2,2)/ul(i,j,k,urho,2,2)
        v = ul(i,j,k,umy,2,2)/ul(i,j,k,urho,2,2)
        w = ul(i,j,k,umz,2,2)/ul(i,j,k,urho,2,2)
        ul(i,j,k,ueint,2,2) = ul(i,j,k,ueden,2,2) - 0.5d0*ul(i,j,k,urho,2,2)*(u**2 + v**2 + w**2)

        u = ul(i,j,k,umx,3,1)/ul(i,j,k,urho,3,1)
        v = ul(i,j,k,umy,3,1)/ul(i,j,k,urho,3,1)
        w = ul(i,j,k,umz,3,1)/ul(i,j,k,urho,3,1)
        ul(i,j,k,ueint,3,1) = ul(i,j,k,ueden,3,1) - 0.5d0*ul(i,j,k,urho,3,1)*(u**2 + v**2 + w**2)

        u = ul(i,j,k,umx,3,2)/ul(i,j,k,urho,3,2)
        v = ul(i,j,k,umy,3,2)/ul(i,j,k,urho,3,2)
        w = ul(i,j,k,umz,3,2)/ul(i,j,k,urho,3,2)
        ul(i,j,k,ueint,3,2) = ul(i,j,k,ueden,3,2) - 0.5d0*ul(i,j,k,urho,3,2)*(u**2 + v**2 + w**2)


!Right Corrected States
        ur(i,j,k,urho:ueden,1,1) = up(i,j,k,urho:ueden,1) - dt/(3.d0*dy)*(flxy(i+1,j+1,k,urho:ueden) - flxy(i+1,j,k,urho:ueden))
        ur(i,j,k,urho:ueden,1,2) = up(i,j,k,urho:ueden,1) - dt/(3.d0*dz)*(flxz(i+1,j,k+1,urho:ueden) - flxz(i+1,j,k,urho:ueden))
        ur(i,j,k,urho:ueden,2,1) = up(i,j,k,urho:ueden,2) - dt/(3.d0*dx)*(flxx(i+1,j+1,k,urho:ueden) - flxx(i,j+1,k,urho:ueden))
        ur(i,j,k,urho:ueden,2,2) = up(i,j,k,urho:ueden,2) - dt/(3.d0*dz)*(flxz(i,j+1,k+1,urho:ueden) - flxz(i,j+1,k,urho:ueden))
        ur(i,j,k,urho:ueden,3,1) = up(i,j,k,urho:ueden,3) - dt/(3.d0*dx)*(flxx(i+1,j,k+1,urho:ueden) - flxx(i,j,k+1,urho:ueden))
        ur(i,j,k,urho:ueden,3,2) = up(i,j,k,urho:ueden,3) - dt/(3.d0*dy)*(flxy(i,j+1,k+1,urho:ueden) - flxy(i,j,k+1,urho:ueden))

!Magnetic Energy will be subtracted off in the magnetic cc 
        u = ur(i,j,k,umx,1,1)/ur(i,j,k,urho,1,1)
        v = ur(i,j,k,umy,1,1)/ur(i,j,k,urho,1,1)
        w = ur(i,j,k,umz,1,1)/ur(i,j,k,urho,1,1)
        ur(i,j,k,ueint,1,1) = ur(i,j,k,ueden,1,1) - 0.5d0*ur(i,j,k,urho,1,1)*(u**2 + v**2 + w**2)

        u = ur(i,j,k,umx,1,2)/ur(i,j,k,urho,1,2)
        v = ur(i,j,k,umy,1,2)/ur(i,j,k,urho,1,2)
        w = ur(i,j,k,umz,1,2)/ur(i,j,k,urho,1,2)
        ur(i,j,k,ueint,1,2) = ur(i,j,k,ueden,1,2) - 0.5d0*ur(i,j,k,urho,1,2)*(u**2 + v**2 + w**2)

        u = ur(i,j,k,umx,2,1)/ur(i,j,k,urho,2,1)
        v = ur(i,j,k,umy,2,1)/ur(i,j,k,urho,2,1)
        w = ur(i,j,k,umz,2,1)/ur(i,j,k,urho,2,1)
        ur(i,j,k,ueint,2,1) = ur(i,j,k,ueden,2,1) - 0.5d0*ur(i,j,k,urho,2,1)*(u**2 + v**2 + w**2)

        u = ur(i,j,k,umx,2,2)/ur(i,j,k,urho,2,2)
        v = ur(i,j,k,umy,2,2)/ur(i,j,k,urho,2,2)
        w = ur(i,j,k,umz,2,2)/ur(i,j,k,urho,2,2)
        ur(i,j,k,ueint,2,2) = ur(i,j,k,ueden,2,2) - 0.5d0*ur(i,j,k,urho,2,2)*(u**2 + v**2 + w**2)

        u = ur(i,j,k,umx,3,1)/ur(i,j,k,urho,3,1)
        v = ur(i,j,k,umy,3,1)/ur(i,j,k,urho,3,1)
        w = ur(i,j,k,umz,3,1)/ur(i,j,k,urho,3,1)
        ur(i,j,k,ueint,3,1) = ur(i,j,k,ueden,3,1) - 0.5d0*ur(i,j,k,urho,3,1)*(u**2 + v**2 + w**2)

        u = ur(i,j,k,umx,3,2)/ur(i,j,k,urho,3,2)
        v = ur(i,j,k,umy,3,2)/ur(i,j,k,urho,3,2)
        w = ur(i,j,k,umz,3,2)/ur(i,j,k,urho,3,2)
        ur(i,j,k,ueint,3,2) = ur(i,j,k,ueden,3,2) - 0.5d0*ur(i,j,k,urho,3,2)*(u**2 + v**2 + w**2)

      enddo
    enddo
  enddo
end subroutine corner_couple

!================================== Use 1D Electric Fields to Transverse correct the Temporary Magnetic Fields ===========================
subroutine corner_couple_mag(lo, hi, &
                             uL, uR, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                             Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                             Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                             Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                             dx, dy, dz, dt)
use amrex_fort_module, only : rt => amrex_real
use meth_params_module, only : qvar, qmagx, qmagy, qmagz

!Correction using Faraday's Law
implicit none

  integer , intent(in   ) :: lo(3), hi(3)
  integer , intent(in   ) :: q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
  integer , intent(in   ) :: ex_l1,ex_l2,ex_l3,ex_h1,ex_h2, ex_h3
  integer , intent(in   ) :: ey_l1,ey_l2,ey_l3,ey_h1,ey_h2, ey_h3
  integer , intent(in   ) :: ez_l1,ez_l2,ez_l3,ez_h1,ez_h2, ez_h3
    real(rt), intent(inout) :: uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2)
    real(rt), intent(inout) :: uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2)
    real(rt), intent(in   ) :: um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
    real(rt), intent(in   ) :: up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)

    real(rt), intent(in   ) :: Ex(ex_l1:ex_h1,ex_l2:ex_h2,ex_l3:ex_h3)
    real(rt), intent(in   ) :: Ey(ey_l1:ey_h1,ey_l2:ey_h2,ey_l3:ey_h3)
    real(rt), intent(in   ) :: Ez(ez_l1:ez_h1,ez_l2:ez_h2,ez_l3:ez_h3)

    real(rt)                :: dx, dy, dz, dt
  integer                 :: i ,j ,k

  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
!Left State
!X-direction 
!-> Affected by Y flux
      ul(i,j,k,qmagx,1,1) = um(i,j,k,qmagx,1) + dt/(3.d0*dx)*(ez(i,j+1,k) - ez(i,j,k))
      ul(i,j,k,qmagz,1,1) = um(i,j,k,qmagz,1) + dt/(6.d0*dx)* &
                            ((ex(i,j+1,k+1) - ex(i,j,k+1)) + &
                             (ex(i,j+1,k  ) - ex(i,j,k)))
!-> Affected by Z flux
      ul(i,j,k,qmagx,1,2) = um(i,j,k,qmagx,1) - dt/(3.d0*dx)*(ey(i,j,k+1) - ey(i,j,k))
      ul(i,j,k,qmagy,1,1) = um(i,j,k,qmagy,1) - dt/(6.d0*dx)*&
                            ((ex(i,j+1,k+1) - ex(i,j+1,k)) + &
                             (ex(i,j  ,k+1) - ex(i,j  ,k)))
!-> Affected by X flux
      ul(i,j,k,qmagy:qmagz,1,2) = um(i,j,k,qmagy:qmagz,1)
!Y-direction
!-> Affected by X flux
      ul(i,j,k,qmagy,2,1) = um(i,j,k,qmagy,2) - dt/(3.d0*dy)*(ez(i+1,j,k) - ez(i,j,k))
      ul(i,j,k,qmagz,2,1) = um(i,j,k,qmagz,2) - dt/(6.d0*dy)*&
                            ((ey(i+1,j,k+1) - ey(i,j,k+1)) + &
                             (ey(i+1,j,k  ) - ey(i,j,k  )))
!-> Affected by Z flux
      ul(i,j,k,qmagy,2,2) = um(i,j,k,qmagy,2) + dt/(3.d0*dy)*(ex(i,j,k+1) - ex(i,j,k))
      ul(i,j,k,qmagx,2,1) = um(i,j,k,qmagx,2) + dt/(6.d0*dy)*&
                              ((ey(i+1,j,k+1) - ey(i+1,j,k)) + &
                               (ey(i  ,j,k+1) - ey(i  ,j,k)))

      ul(i,j,k,qmagx,2,2) = um(i,j,k,qmagx,2)
      ul(i,j,k,qmagz,2,2) = um(i,j,k,qmagz,2)

!Z-Direction
!-> Affected by X flux
      ul(i,j,k,qmagz,3,1) = um(i,j,k,qmagz,3) + dt/(3.d0*dz)*(ey(i+1,j,k) - ey(i,j,k))
      ul(i,j,k,qmagy,3,1) = um(i,j,k,qmagz,3) + dt/(6.d0*dz)*&
                            ((ez(i+1,j+1,k) - ez(i,j+1,k)) + &
                             (ez(i+1,j  ,k) - ez(i,j  ,k)))
!-> Affected by Y flux
      ul(i,j,k,qmagz,3,2) = um(i,j,k,qmagz,3) - dt/(3.d0*dz)*(ex(i,j+1,k) - ex(i,j,k))
      ul(i,j,k,qmagx,3,1) = um(i,j,k,qmagx,3) - dt/(6.d0*dz)*&
                              ((ez(i+1,j+1,k) - ez(i+1,j,k)) + &
                               (ez(i  ,j+1,k) - ez(i  ,j,k)))

      ul(i,j,k,qmagx:qmagy,3,2) = um(i,j,k,qmagy:qmagz,3)
      ul(i,j,k,ueint,1,1) = ul(i,j,k,ueint,1,1) - 0.5d0*dot_product(ul(i,j,k,qmagx:qmagz,1,1),ul(i,j,k,qmagx:qmagz,1,1))

      ul(i,j,k,ueint,1,2) = ul(i,j,k,ueint,1,2) - 0.5d0*dot_product(ul(i,j,k,qmagx:qmagz,1,2),ul(i,j,k,qmagx:qmagz,1,2))

      ul(i,j,k,ueint,2,1) = ul(i,j,k,ueint,2,1) - 0.5d0*dot_product(ul(i,j,k,qmagx:qmagz,2,1),ul(i,j,k,qmagx:qmagz,2,1))

      ul(i,j,k,ueint,2,2) = ul(i,j,k,ueint,2,2) - 0.5d0*dot_product(ul(i,j,k,qmagx:qmagz,2,2),ul(i,j,k,qmagx:qmagz,2,2))

      ul(i,j,k,ueint,3,1) = ul(i,j,k,ueint,3,1) - 0.5d0*dot_product(ul(i,j,k,qmagx:qmagz,3,1),ul(i,j,k,qmagx:qmagz,3,1))

      ul(i,j,k,ueint,3,2) = ul(i,j,k,ueint,3,2) - 0.5d0*dot_product(ul(i,j,k,qmagx:qmagz,3,2),ul(i,j,k,qmagx:qmagz,3,2))
!Right State
!X-direction 
!-> Affected by Y flux
      ur(i,j,k,qmagx,1,1) = up(i,j,k,qmagx,1) + dt/(3.d0*dx)*(ez(i+1,j+1,k) - ez(i+1,j,k))
      ur(i,j,k,qmagz,1,1) = up(i,j,k,qmagz,1) + dt/(6.d0*dx)*&
                              ((ex(i,j+1,k+1) - ex(i,j,k+1)) + &
                               (ex(i,j+1,k  ) - ex(i,j,k  )))
!-> Affected by Z flux
      ur(i,j,k,qmagx,1,2) = up(i,j,k,qmagx,1) - dt/(3.d0*dx)*(ey(i+1,j,k+1) - ey(i+1,j,k))
      ur(i,j,k,qmagy,1,1) = up(i,j,k,qmagy,1) - dt/(6.d0*dx)*&
                              ((ex(i,j+1,k+1) - ex(i,j+1,k)) + &
                               (ex(i,j  ,k+1) - ex(i,j  ,k)))

      ur(i,j,k,qmagy:qmagz,1,2) = up(i,j,k,qmagy:qmagz,1)
!Y-direction
!-> Affected by X flux
      ur(i,j,k,qmagy,2,1) = up(i,j,k,qmagy,2) - dt/(3.d0*dy)*(ez(i+1,j+1,k) - ez(i,j+1,k))
      ur(i,j,k,qmagz,2,1) = up(i,j,k,qmagz,2) - dt/(6.d0*dy)*&
                              ((ey(i+1,j,k+1) - ey(i,j,k+1)) + &
                               (ey(i+1,j,k  ) - ey(i,j,k  )))
!-> Affected by Z flux
      ur(i,j,k,qmagy,2,2) = up(i,j,k,qmagy,2) + dt/(3.d0*dy)*(ex(i,j+1,k+1) - ex(i,j+1,k))
      ur(i,j,k,qmagx,2,1) = up(i,j,k,qmagx,2) + dt/(6.d0*dy)*&
                              ((ey(i+1,j,k+1) - ey(i+1,j,k)) + &
                               (ey(i  ,j,k+1) - ey(i  ,j,k)))

      ur(i,j,k,qmagx,2,2) = up(i,j,k,qmagx,2)
      ur(i,j,k,qmagz,2,2) = up(i,j,k,qmagz,2)

!Z-Direction
!-> Affected by X flux
      ur(i,j,k,qmagz,3,1) = up(i,j,k,qmagz,3) + dt/(3.d0*dz)*(ey(i+1,j,k+1) - ey(i,j,k+1))
      ur(i,j,k,qmagy,3,1) = up(i,j,k,qmagz,3) + dt/(6.d0*dz)*&
                              ((ez(i+1,j+1,k) - ez(i,j+1,k)) + &
                               (ez(i+1,j  ,k) - ez(i,j  ,k)))
!-> Affected by Y flux
      ur(i,j,k,qmagz,3,2) = up(i,j,k,qmagz,3) - dt/(3.d0*dz)*(ex(i,j+1,k+1) - ex(i,j,k+1))
      ur(i,j,k,qmagx,3,1) = up(i,j,k,qmagx,3) - dt/(6.d0*dz)*&
                            ((ez(i+1,j+1,k) - ez(i+1,j,k)) + &
                             (ez(i  ,j+1,k) - ez(i  ,j,k)))

      ur(i,j,k,qmagx:qmagy,3,2) = up(i,j,k,qmagy:qmagz,3)
      ur(i,j,k,ueint,1,1) = ur(i,j,k,ueint,1,1) - 0.5d0*dot_product(ur(i,j,k,qmagx:qmagz,1,1),ur(i,j,k,qmagx:qmagz,1,1))
      ur(i,j,k,ueint,1,2) = ur(i,j,k,ueint,1,2) - 0.5d0*dot_product(ur(i,j,k,qmagx:qmagz,1,2),ur(i,j,k,qmagx:qmagz,1,2))
      ur(i,j,k,ueint,2,1) = ur(i,j,k,ueint,2,1) - 0.5d0*dot_product(ur(i,j,k,qmagx:qmagz,2,1),ur(i,j,k,qmagx:qmagz,2,1))
      ur(i,j,k,ueint,2,2) = ur(i,j,k,ueint,2,2) - 0.5d0*dot_product(ur(i,j,k,qmagx:qmagz,2,2),ur(i,j,k,qmagx:qmagz,2,2))
      ur(i,j,k,ueint,3,1) = ur(i,j,k,ueint,3,1) - 0.5d0*dot_product(ur(i,j,k,qmagx:qmagz,3,1),ur(i,j,k,qmagx:qmagz,3,1))
      ur(i,j,k,ueint,3,2) = ur(i,j,k,ueint,3,2) - 0.5d0*dot_product(ur(i,j,k,qmagx:qmagz,3,2),ur(i,j,k,qmagx:qmagz,3,2))
      enddo
    enddo
  enddo

end subroutine corner_couple_mag

!====================================================== Final Conservative Corrections================================================================
subroutine half_step(lo, hi, &
         uL, uR, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
         flxx, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
         flxy, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
         flxz, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
         dx, dy, dz, dt)
use amrex_fort_module, only : rt => amrex_real
use meth_params_module, only : qvar, urho, ueden

implicit none

  integer, intent(in )  :: lo(3), hi(3),q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
  integer, intent(in )  :: flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3
  integer, intent(in )  :: flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3
  integer, intent(in )  :: flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3

    real(rt), intent(in ) ::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
    real(rt), intent(in ) ::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)

    real(rt), intent(in ) :: flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,QVAR,2)
    real(rt), intent(in ) :: flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,QVAR,2)
    real(rt), intent(in ) :: flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,QVAR,2)

    real(rt), intent(out) ::uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
    real(rt), intent(out) ::uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)

    real(rt)              :: dx, dy, dz, dt, u, v, w
  integer               :: i ,j ,k, n

  ul = um
  ur = up
  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
!left state             
        ul(i,j,k,urho:ueden,1) = um(i,j,k,urho:ueden,1) - 0.5d0*dt/dy*(flxy(i,j+1,k,urho:ueden,2) - flxy(i,j,k,urho:ueden,2)) &
                               - 0.5d0*dt/dz*(flxz(i,j,k+1,urho:ueden,2) - flxz(i,j,k,urho:ueden,2))
        ul(i,j,k,urho:ueden,2) = um(i,j,k,urho:ueden,2) - 0.5d0*dt/dx*(flxx(i+1,j,k,urho:ueden,2) - flxx(i,j,k,urho:ueden,2)) &
                               - 0.5d0*dt/dz*(flxz(i,j,k+1,urho:ueden,1) - flxz(i,j,k,urho:ueden,1))
        ul(i,j,k,urho:ueden,3) = um(i,j,k,urho:ueden,3) - 0.5d0*dt/dx*(flxx(i+1,j,k,urho:ueden,1) - flxx(i,j,k,urho:ueden,1)) &
                               - 0.5d0*dt/dy*(flxy(i,j+1,k,urho:ueden,1) - flxy(i,j,k,urho:ueden,1))

! Magnetic energy is subtracted in half_step_mag
        do n = 1,3
          u = ul(i,j,k,umx,n)/ul(i,j,k,urho,n)
          v = ul(i,j,k,umy,n)/ul(i,j,k,urho,n)
          w = ul(i,j,k,umz,n)/ul(i,j,k,urho,n)
          ul(i,j,k,ueint,n) = ul(i,j,k,ueden,n) - 0.5d0*ul(i,j,k,urho,n)*(u**2 + v**2 + w**2)
        enddo
!right state                
        ur(i,j,k,urho:ueden,1) = up(i,j,k,urho:ueden,1) - 0.5d0*dt/dy*(flxy(i,j+1,k,urho:ueden,2) - flxy(i,j,k,urho:ueden,2)) &
                               - 0.5d0*dt/dz*(flxz(i,j,k+1,urho:ueden,2) - flxz(i,j,k,urho:ueden,2))
        ur(i,j,k,urho:ueden,2) = up(i,j,k,urho:ueden,2) - 0.5d0*dt/dx*(flxx(i+1,j,k,urho:ueden,2) - flxx(i,j,k,urho:ueden,2)) &
                               - 0.5d0*dt/dz*(flxz(i,j,k+1,urho:ueden,1) - flxz(i,j,k,urho:ueden,1))
        ur(i,j,k,urho:ueden,3) = up(i,j,k,urho:ueden,3) - 0.5d0*dt/dx*(flxx(i+1,j,k,urho:ueden,1) - flxx(i,j,k,urho:ueden,1)) &
                               - 0.5d0*dt/dy*(flxy(i,j+1,k,urho:ueden,1) - flxy(i,j,k,urho:ueden,1))
        do n = 1,3
          u = ur(i,j,k,umx,n)/ur(i,j,k,urho,n)
          v = ur(i,j,k,umy,n)/ur(i,j,k,urho,n)
          w = ur(i,j,k,umz,n)/ur(i,j,k,urho,n)
          ur(i,j,k,ueint,n) = ur(i,j,k,ueden,n) - 0.5d0*ur(i,j,k,urho,n)*(u**2 + v**2 + w**2)
        enddo
      enddo
    enddo
  enddo
end subroutine 

!================================================= Final Magnetic Corrections ========================================================================
subroutine half_step_mag(lo, hi, &
                 uL, uR, um, up, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                 Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                 Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                 Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                 dx, dy, dz, dt)
use amrex_fort_module, only : rt => amrex_real
use meth_params_module, only : qvar, qmagx,qmagy,qmagz

!Correction using Faraday's Law
implicit none

    integer, intent(in)   :: lo(3), hi(3), q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
    integer, intent(in)   :: ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3
    integer, intent(in)   :: ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3
    integer, intent(in)   :: ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3

    real(rt), intent(inout)     ::uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
    real(rt), intent(inout)     ::uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
    real(rt), intent(in)        ::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
    real(rt), intent(in)        ::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)  

    real(rt), intent(in)  :: Ex(ex_l1:ex_h1,ex_l2:ex_h2,ex_l3:ex_h3)
    real(rt), intent(in)  :: Ey(ey_l1:ey_h1,ey_l2:ey_h2,ey_l3:ey_h3)
    real(rt), intent(in)  :: Ez(ez_l1:ez_h1,ez_l2:ez_h2,ez_l3:ez_h3)

    real(rt)                    :: dx, dy, dz, dt
    integer                     :: i ,j ,k, n

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
        !---------------------------------------left state-----------------------------------------------------
                !X-Direction
                !Bx
                ul(i,j,k,qmagx,1) = um(i,j,k,qmagx,1) - 0.5d0*dt/dx*((ey(i,j,k+1) - ey(i,j,k)) &
                                                   - (ez(i,j+1,k) - ez(i,j,k)))
                !By
                ul(i,j,k,qmagy,1) = um(i,j,k,qmagy,1) - 0.5d0*dt/dx*((ex(i  ,j+1,k+1) - ex(i,j+1,k)) &
                                               + (ex(i  ,j  ,k+1) - ex(i,j  ,k)) &
                                                                                   - (ez(i+1,j+1,k  ) - ez(i,j+1,k)) &
                                               - (ez(i+1,j  ,k  ) - ez(i,j  ,k)))
                !Bz
                ul(i,j,k,qmagz,1) = um(i,j,k,qmagz,1) + 0.5d0*dt/dx*((ex(i  ,j+1,k+1) - ex(i,j,k+1)) &
                                               + (ex(i  ,j+1,k  ) - ex(i,j,k  )) &
                                                                   - (ey(i+1,j  ,k+1) - ey(i,j,k+1)) &
                                               - (ey(i+1,j  ,k  ) - ey(i,j,k  )))
                !Y-Direction
                !Bx
                ul(i,j,k,qmagx,2) = um(i,j,k,qmagx,2) + 0.5d0*dt/dy*((ey(i+1,j,k+1) - ey(i+1,j,k)) &
                                                   + (ey(i  ,j,k+1) - ey(i  ,j,k)) &
                                                                                   - (ez(i+1,j+1,k) - ez(i+1,j,k)) &
                                               - (ez(i  ,j+1,k) - ez(i  ,j,k)))
                !By
                ul(i,j,k,qmagy,2) = um(i,j,k,qmagy,2) + 0.5d0*dt/dy*((ex(i,j,k+1) - ex(i,j,k)) &
                                                   - (ez(i+1,j,k) - ez(i,j,k)))
                !Bz
                ul(i,j,k,qmagz,2) = um(i,j,k,qmagz,2) - 0.5d0*dt/dy*((ey(i+1,j,k+1) - ey(i,j,k+1)) &
                                                   + (ey(i+1,j,k  ) - ey(i,j,k  )) &
                                                                                   - (ex(i,j+1,k+1) - ex(i,j,k+1)) &
                                               - (ex(i,j+1,k  ) - ex(i,j,k  )))
                !Z-direction
                !Bx
                ul(i,j,k,qmagx,3) = um(i,j,k,qmagx,3) - 0.5d0*dt/dz*((ez(i+1,j+1,k) - ez(i+1,j,k)) &
                                                   + (ez(i  ,j+1,k) - ez(i  ,j,k)) &
                                                                                   - (ey(i+1,j,k+1) - ey(i+1,j,k)) &
                                               - (ey(i  ,j+1,k) - ey(i  ,j,k)))
                !By
                ul(i,j,k,qmagy,3) = um(i,j,k,qmagy,3) + 0.5d0*dt/dz*((ez(i+1,j+1,k  ) - ez(i,j+1,k)) &
                                                    +(ez(i+1,j  ,k  ) - ez(i,j  ,k)) &
                                                                                    -(ex(i  ,j+1,k+1) - ex(i,j+1,k)) &
                                               - (ex(i  ,j  ,k+1) - ex(i,j  ,k)))
                !Bz
                ul(i,j,k,qmagz,3) = um(i,j,k,qmagz,3) - 0.5d0*dt/dz*((ex(i,j+1,k) - ex(i,j,k)) &
                                               - (ey(i+1,j,k) - ey(i,j,k)))
                do n = 1,3
                    ul(i,j,k,ueint,n) = ul(i,j,k,ueint,n) - 0.5d0*dot_product(ul(i,j,k,qmagx:qmagz,n),ul(i,j,k,qmagx:qmagz,n))
                enddo

    !---------------------------------------right state-----------------------------------------------------
                !X-Direction
                !Bx
                ur(i,j,k,qmagx,1) = up(i,j,k,qmagx,1) - 0.5d0*dt/dx*((ey(i+1,j,k+1) - ey(i+1,j,k)) &
                                    -                (ez(i+1,j+1,k) - ez(i+1,j,k)))
                !By
                ur(i,j,k,qmagy,1) = up(i,j,k,qmagy,1) - 0.5d0*dt/dx*((ex(i  ,j+1,k+1) - ex(i,j+1,k)) &
                                               + (ex(i  ,j  ,k+1) - ex(i,j  ,k)) &
                                                                                   - (ez(i+1,j+1,k  ) - ez(i,j+1,k)) &
                                               - (ez(i+1,j  ,k  ) - ez(i,j  ,k)))
                !Bz
                ur(i,j,k,qmagz,1) = up(i,j,k,qmagz,1) + 0.5d0*dt/dx*((ex(i,j+1,k+1) - ex(i,j,k+1)) &
                                               + (ex(i,j+1,k  ) - ex(i,j,k  )) &
                                                                                   - (ey(i+1,j,k+1) - ey(i,j,k+1)) &
                                               - (ey(i+1,j,k  ) - ey(i,j,k  )))             
                !Y-Direction
                !Bx
                ur(i,j,k,qmagx,2) = up(i,j,k,qmagx,2) + 0.5d0*dt/dy*((ey(i+1,j,k+1) - ey(i+1,j,k)) &
                                                   + (ey(i  ,j,k+1) - ey(i  ,j,k)) &
                                                                                   - (ez(i+1,j+1,k) - ez(i+1,j,k)) &
                                               - (ez(i  ,j+1,k) - ez(i,j,k)))
                !By
                ur(i,j,k,qmagy,2) = up(i,j,k,qmagy,2) + 0.5d0*dt/dy*((ex(i,j+1,k+1) - ex(i,j+1,k)) &
                                    - (ez(i+1,j+1,k) - ez(i,j+1,k)))
                !Bz
                ur(i,j,k,qmagz,2) = up(i,j,k,qmagz,2) - 0.5d0*dt/dy*((ey(i+1,j,k+1) - ey(i,j,k+1)) &
                                        + (ey(i+1,j,k) - ey(i,j,k)) &
                                        - (ex(i,j+1,k+1) - ex(i,j,k+1)) &
                                    - (ex(i,j+1,k) - ex(i,j,k)))    
        
                !Z-direction
                !Bx
                ur(i,j,k,qmagx,3) = up(i,j,k,qmagx,3) - 0.5d0*dt/dz*((ez(i+1,j+1,k) - ez(i+1,j,k)) &
                                            + (ez(i,j+1,k) - ez(i,j,k)) &
                                    - (ey(i+1,j,k+1) - ey(i+1,j,k)) &
                                    - (ey(i,j,k+1) - ey(i,j,k)))
                !By
                ur(i,j,k,qmagy,3) = up(i,j,k,qmagy,3) + 0.5d0*dt/dz*((ez(i+1,j+1,k  ) - ez(i,j+1,k)) &
                                                   + (ez(i+1,j  ,k  ) - ez(i,j  ,k)) &
                                                                                   - (ex(i  ,j+1,k+1) - ex(i,j+1,k)) &
                                               - (ex(i  ,j  ,k+1) - ex(i,j,k)))
                !Bz
                ur(i,j,k,qmagz,3) = up(i,j,k,qmagz,3) - 0.5d0*dt/dz*((ex(i,j+1,k+1) - ex(i,j,k+1)) &
                                               - (ey(i+1,j,k+1) - ey(i,j,k+1)))         
                do n = 1,3
                    ur(i,j,k,ueint,n) = ur(i,j,k,ueint,n) - 0.5d0*dot_product(ur(i,j,k,qmagx:qmagz,n),ur(i,j,k,qmagx:qmagz,n))
                enddo   
            enddo
        enddo
    enddo

end subroutine half_step_mag

!================================== Find the 2D corrected primitive variables =======================================
subroutine prim_half(q2D,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                 flxx, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
             flxy, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
             flxz, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
             dx, dy, dz, dt)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module, only : qvar

implicit none

    integer, intent(in)     ::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3

    integer, intent(in)   :: flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3
    integer, intent(in)   :: flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3
    integer, intent(in)   :: flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3

    real(rt), intent(in)    :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
    real(rt), intent(in)    :: flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,QVAR)
    real(rt), intent(in)    :: flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,QVAR)
    real(rt), intent(in)    :: flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,QVAR)

    real(rt), intent(out)   ::q2D(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)

    real(rt)                ::flx_sum(QVAR)
    real(rt)                ::qflx(QVAR)
    real(rt)                :: dx, dy, dz, dt   
    integer                 ::i, j, k
!   q2D = q
    do k = q_l3+1,q_h3-1
        do j = q_l2+1,q_h2-1
            do i = q_l1+1,q_h1-1
                flx_sum = (flxx(i+1,j,k,:) - flxx(i,j,k,:)) + (flxy(i,j+1,k,:) - flxy(i,j,k,:)) + (flxz(i,j,k+1,:) - flxz(i,j,k,:)) 
                call qflux(qflx,flx_sum,q(i,j,k,:))
                q2d(i,j,k,:) = q(i,j,k,:) - 0.5d0*dt/dx*qflx
            enddo
        enddo
    enddo
end subroutine prim_half


!================================= Calculate the C to P Jacobian applied to the fluxes ===================================

subroutine qflux(qflx,flx,q)
 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none

 real(rt), intent(in)       ::flx(QVAR), q(QVAR)
 real(rt), intent(out)      ::qflx(QVAR)
    qflx = 0.d0
    qflx(qrho)  = flx(urho)
    qflx(qu)    = q(qu)/q(qrho)*flx(urho) + 1.d0/q(qrho)*flx(umx)
    qflx(qv)    = q(qv)/q(qrho)*flx(qrho) + 1.d0/q(qrho)*flx(umy)
    qflx(qw)    = q(qw)/q(qrho)*flx(qrho) + 1.d0/q(qrho)*flx(umz)
    qflx(qpres) = 0.5d0*(q(qu)**2+q(qv)**2+q(qw)**2)*flx(urho) - q(qu)*flx(umx) - q(qv)*flx(umy) - q(qw)*flx(umz) - flx(ueden)  &
                  -q(qmagx)*flx(qmagx) - q(qmagy)*flx(qmagy) - q(qmagz)*flx(qmagz)
    qflx(qpres) = qflx(qpres)*gamma_minus_1               
    qflx(qmagx) = flx(qmagx)
    qflx(qmagy) = flx(qmagy)
    qflx(qmagz) = flx(qmagz)

end subroutine qflux


!============================================ Debug code =====================================================
    subroutine checkisnanmult(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, num)
       use amrex_fort_module, only : rt => amrex_real

    implicit none
    integer, intent(in)  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, num
    real(rt), intent(in) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,num)

    integer :: i,j,k,n


    do n = 1,num
        do k = uout_l3,uout_h3
            do j = uout_l2, uout_h2
                do i = uout_l1,uout_h1
                    if(isnan(uout(i,j,k,n)).or.(abs(uout(i,j,k,n)).ge. 1d14)) then
                        write(*,*) "Bad values ",  uout(i,j,k,:)
                        write(*,*) "Failure to converge ", "i, j, k, n = ", i, j, k, n
                        stop
                    endif
                enddo
            enddo
        enddo
    enddo
    end subroutine checkisnanmult
!============ single =====================  

    subroutine checkisnans(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3)
       use amrex_fort_module, only : rt => amrex_real

    implicit none
    integer, intent(in)  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
    real(rt), intent(in) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3)

    integer :: i,j,k

        do k = uout_l3,uout_h3
            do j = uout_l2, uout_h2
                do i = uout_l1,uout_h1
                    if(isnan(uout(i,j,k)).or.(abs(uout(i,j,k)).ge. 1d14)) then
                        write(*,*) "Bad values ",  uout(i,j,k)
                        write(*,*) "Failure to converge ", "i, j, k = ", i, j, k
                        stop
                    endif
                enddo
            enddo
        enddo
    end subroutine checkisnans

!====================================== Density Check ========================================
    subroutine checknegdens(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3)
       use amrex_fort_module, only : rt => amrex_real

    implicit none
    integer, intent(in)  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
    real(rt), intent(in) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3)

    integer :: i,j,k

        do k = uout_l3,uout_h3
            do j = uout_l2, uout_h2
                do i = uout_l1,uout_h1
                    if(uout(i,j,k).le. 0.d0) then
                        write(*,*) "Non-Positive Density ",  uout(i,j,k)
                        write(*,*) "i, j, k = ", i, j, k
                        stop
                    endif
                enddo
            enddo
        enddo
    end subroutine checknegdens
end module ct_upwind
````

