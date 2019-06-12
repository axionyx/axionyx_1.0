
# File mhd\_plm\_3d.f90

[**File List**](files.md) **>** [**MHD**](dir_a5db59ee0cc93a408ad0433ba32613c6.md) **>** [**mhd\_plm\_3d.f90**](mhd__plm__3d_8f90.md)

[Go to the documentation of this file.](mhd__plm__3d_8f90.md) 


````cpp
module mhd_plm_module

!Module that gives a piecewise linear interpolation for the primitive variables 
!They are projected onto the characteristic variables for tracing. 

   use amrex_fort_module, only : rt => amrex_real
   use meth_params_module
   implicit none

   private vanleer, lvecx, lvecy, lvecz, rvecx, rvecy, rvecz, evals

   public plm 

contains 

  !
  ! characteristics based on u
  !
  !===========================================================================
  ! This is called from within threaded loops in advance_mhd_tile so *no* OMP here ...
  !===========================================================================

   subroutine plm(lo, hi, s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,&
          bx, bxl1, bxl2, bxl3, bxh1, bxh2, bxh3, &
          by, byl1, byl2, byl3, byh1, byh2, byh3, &
          bz, bzl1, bzl2, bzl3, bzh1, bzh2, bzh3, &
                  Ip,Im, ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dy,dz,dt,a_old)


    implicit none
    integer , intent(in   ) ::    s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, lo(3), hi(3)
    integer , intent(in   ) ::    ilo1, ilo2, ilo3, ihi1, ihi2, ihi3
    integer , intent(in   ) ::    bxl1, bxl2, bxl3, bxh1, bxh2, bxh3
    integer , intent(in   ) ::    byl1, byl2, byl3, byh1, byh2, byh3
    integer , intent(in   ) ::    bzl1, bzl2, bzl3, bzh1, bzh2, bzh3
 
    real(rt), intent(in   ) ::    s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QVAR) !Primitive Vars
    real(rt), intent(in   ) ::    bx(bxl1:bxh1, bxl2:bxh2, bxl3:bxh3)   !Face Centered Magnetic Fields
    real(rt), intent(in   ) ::    by(byl1:byh1, byl2:byh2, byl3:byh3)
    real(rt), intent(in   ) ::    bz(bzl1:bzh1, bzl2:bzh2, bzl3:bzh3)

    real(rt), intent(out )  ::    Ip(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,QVAR,3)
    real(rt), intent(out )  ::    Im(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,QVAR,3)

    real(rt), intent(in   ) ::    dx,dy,dz,dt,a_old

    real(rt)                ::    dQL(7), dQR(7), dW, dL, dR, leig(7,7), reig(7,7), lam(7), summ(7)
    real(rt)                ::    temp(s_l1-1:s_h1+1,s_l2-1:s_h2+1,s_l3-1:s_h3+1,8), smhd(7)
    real(rt)                ::    tbx(s_l1-1:s_h1+1,s_l2-1:s_h2+1,s_l3-1:s_h3+1)
    real(rt)                ::    tby(s_l1-1:s_h1+1,s_l2-1:s_h2+1,s_l3-1:s_h3+1)
    real(rt)                ::    tbz(s_l1-1:s_h1+1,s_l2-1:s_h2+1,s_l3-1:s_h3+1)
    real(rt)                ::    dt_over_a
    integer                 ::    ii,ibx,iby,ibz, i , j, k

    ibx = 6
    iby = 7
    ibz = 8

    dt_over_a = dt / a_old

!------------------------workspace variables---------------------------------------------
    temp = 0.d0
    temp(:,:,:,1) = small_dens
    temp(s_l1: s_h1, s_l2: s_h2, s_l3: s_h3,1  :  5) = s(:,:,:,qrho:qpres) !Gas vars Cell Centered
    temp(s_l1: s_h1, s_l2: s_h2, s_l3: s_h3,ibx:ibz) = s(:,:,:,qmagx:qmagz) !Mag vars Cell Centered
!-------------------- Fill Boundaries ---------------------------------------------------
    temp(s_l1-1,s_l2-1,s_l3-1,1:5) = s(s_l1,s_l2,s_l3,qrho:qpres)
    temp(s_l1-1,s_l2-1,s_l3-1,6:8) = s(s_l1,s_l2,s_l3,qmagx:qmagz)
    temp(s_l1-1, s_l2: s_h2, s_l3: s_h3,1:5) = s(s_l1,:,:,qrho:qpres)
    temp(s_l1-1, s_l2: s_h2, s_l3: s_h3,6:8) = s(s_l1,:,:,qmagx:qmagz)
    temp(s_l1:s_h1, s_l2-1, s_l3: s_h3,1:5) = s(:,s_l2,:,qrho:qpres)
    temp(s_l1:s_h1, s_l2-1, s_l3: s_h3,6:8) = s(:,s_l2,:,qmagx:qmagz)
    temp(s_l1:s_h1, s_l2:s_h2, s_l3-1,1:5) = s(:,:,s_l3,qrho:qpres)
    temp(s_l1:s_h1, s_l2:s_h2, s_l3-1,6:8) = s(:,:,s_l3,qmagx:qmagz)
    temp(s_h1+1, s_l2: s_h2, s_l3: s_h3,1:5) = s(s_h1,:,:,qrho:qpres)
    temp(s_h1+1, s_l2: s_h2, s_l3: s_h3,6:8) = s(s_h1,:,:,qmagx:qmagz)
    temp(s_l1:s_h1, s_h2+1, s_l3: s_h3,1:5) = s(:,s_h2,:,qrho:qpres)
    temp(s_l1:s_h1, s_h2+1, s_l3: s_h3,6:8) = s(:,s_h2,:,qmagx:qmagz)
    temp(s_l1:s_h1, s_l2:s_h2, s_h3+1,1:5) = s(:,:,s_h3,qrho:qpres)
    temp(s_l1:s_h1, s_l2:s_h2, s_h3+1,6:8) = s(:,:,s_h3,qmagx:qmagz)
    temp(s_h1+1,s_h2+1,s_h3+1,1:5) = s(s_h1,s_h2,s_h3,qrho:qpres)
    temp(s_h1+1,s_h2+1,s_h3+1,6:8) = s(s_h1,s_h2,s_h3,qmagx:qmagz)
! Temp face centered magnetic fields
  do k = s_l3-1, s_h3+1
    do j = s_l2-1, s_h2+1
      do i = s_l1-1, s_h1+1
!---------------------------------------------- set up temporary bx ------------------------------------------------
      if(i.lt.bxl1) then 
        if(j.lt.bxl2) then
          if(k.lt.bxl3) then 
            tbx(i,j,k) = bx(bxl1,bxl2,bxl3)
          elseif(k.lt.bxh3) then
            tbx(i,j,k) = bx(bxl1, bxl2, k)
          else
            tbx(i,j,k) = bx(bxl1, bxl2, bxh3)
          endif
        elseif(j.lt.bxh2) then
          if(k.lt.bxl3) then 
            tbx(i,j,k) = bx(bxl1,j,bxl3)
          elseif(k.lt.bxh3) then
            tbx(i,j,k) = bx(bxl1, j, k)
          else
            tbx(i,j,k) = bx(bxl1, j, bxh3)
          endif
        else
          if(k.lt.bxl3) then 
            tbx(i,j,k) = bx(bxl1,bxh2,bxl3)
          elseif(k.lt.bxh3) then
            tbx(i,j,k) = bx(bxl1, bxh2, k)
          else
            tbx(i,j,k) = bx(bxl1, bxh2, bxh3)
          endif
        endif
      elseif(i.lt.bxh1) then 
        if(j.lt.bxl2) then
          if(k.lt.bxl3) then 
            tbx(i,j,k) = bx(i,bxl2,bxl3)
          elseif(k.lt.bxh3) then
            tbx(i,j,k) = bx(i, bxl2, k)
          else
            tbx(i,j,k) = bx(i, bxl2, bxh3)
          endif
        elseif(j.lt.bxh2) then
          if(k.lt.bxl3) then 
            tbx(i,j,k) = bx(i,j,bxl3)
          elseif(k.lt.bxh3) then
            tbx(i,j,k) = bx(i, j, k)
          else
            tbx(i,j,k) = bx(i, j, bxh3)
          endif
        else
          if(k.lt.bxl3) then 
            tbx(i,j,k) = bx(i,bxh2,bxl3)
          elseif(k.lt.bxh3) then
            tbx(i,j,k) = bx(i, bxh2, k)
          else
            tbx(i,j,k) = bx(i, bxh2, bxh3)
        endif
      endif
    else
      if(j.lt.bxl2) then
        if(k.lt.bxl3) then 
          tbx(i,j,k) = bx(bxh1,bxl2,bxl3)
        elseif(k.lt.bxh3) then
          tbx(i,j,k) = bx(bxh1, bxl2, k)
        else
          tbx(i,j,k) = bx(bxh1, bxl2, bxh3)
        endif
      elseif(j.lt.bxh2) then
        if(k.lt.bxl3) then 
          tbx(i,j,k) = bx(bxh1,j,bxl3)
        elseif(k.lt.bxh3) then
          tbx(i,j,k) = bx(bxh1, j, k)
        else
                            tbx(i,j,k) = bx(bxh1, j, bxh3)
                        endif
                    else
                        if(k.lt.bxl3) then 
                            tbx(i,j,k) = bx(bxh1,bxh2,bxl3)
                        elseif(k.lt.bxh3) then
                            tbx(i,j,k) = bx(bxh1, bxh2, k)
                        else
                            tbx(i,j,k) = bx(bxh1, bxh2, bxh3)
                        endif
                    endif
                endif

    !---------------------------------------------- set up temporary by ------------------------------------------------
                if(i.lt.byl1) then 
                    if(j.lt.byl2) then
                        if(k.lt.byl3) then 
                            tby(i,j,k) = by(byl1,byl2,byl3)
                        elseif(k.lt.byh3) then
                            tby(i,j,k) = by(byl1, byl2, k)
                        else
                            tby(i,j,k) = by(byl1, byl2, byh3)
                        endif
                    elseif(j.lt.byh2) then
                        if(k.lt.byl3) then 
                            tby(i,j,k) = by(byl1,j,byl3)
                        elseif(k.lt.byh3) then
                            tby(i,j,k) = by(byl1, j, k)
                        else
                            tby(i,j,k) = by(byl1, j, byh3)
                        endif
                    else
                        if(k.lt.byl3) then 
                            tby(i,j,k) = by(byl1,byh2,byl3)
                        elseif(k.lt.byh3) then
                            tby(i,j,k) = by(byl1, byh2, k)
                        else
                            tby(i,j,k) = by(byl1, byh2, byh3)
                        endif
                    endif
                elseif(i.lt.byh1) then 
                    if(j.lt.byl2) then
                        if(k.lt.byl3) then 
                            tby(i,j,k) = by(i,byl2,byl3)
                        elseif(k.lt.byh3) then
                            tby(i,j,k) = by(i, byl2, k)
                        else
                            tby(i,j,k) = by(i, byl2, byh3)
                        endif
                    elseif(j.lt.byh2) then
                        if(k.lt.byl3) then 
                            tby(i,j,k) = by(i,j,byl3)
                        elseif(k.lt.byh3) then
                            tby(i,j,k) = by(i, j, k)
                        else
                            tby(i,j,k) = by(i, j, byh3)
                        endif
                    else
                        if(k.lt.byl3) then 
                            tby(i,j,k) = by(i,byh2,byl3)
                        elseif(k.lt.byh3) then
                            tby(i,j,k) = by(i, byh2, k)
                        else
                            tby(i,j,k) = by(i, byh2, byh3)
                        endif
                    endif
                else
                    if(j.lt.byl2) then
                        if(k.lt.byl3) then 
                            tby(i,j,k) = by(byh1,byl2,byl3)
                        elseif(k.lt.byh3) then
                            tby(i,j,k) = by(byh1, byl2, k)
                        else
                            tby(i,j,k) = by(byh1, byl2, byh3)
                        endif
                    elseif(j.lt.byh2) then
                        if(k.lt.byl3) then 
                            tby(i,j,k) = by(byh1,j,byl3)
                        elseif(k.lt.byh3) then
                            tby(i,j,k) = by(byh1, j, k)
                        else
                            tby(i,j,k) = by(byh1, j, byh3)
                        endif
                    else
                        if(k.lt.byl3) then 
                            tby(i,j,k) = by(byh1,byh2,byl3)
                        elseif(k.lt.byh3) then
                            tby(i,j,k) = by(byh1, byh2, k)
                        else
                            tby(i,j,k) = by(byh1, byh2, byh3)
                        endif
                    endif
                endif

    !---------------------------------------------- set up temporary bz ------------------------------------------------
                if(i.lt.bzl1) then 
                    if(j.lt.bzl2) then
                        if(k.lt.bzl3) then 
                            tbz(i,j,k) = bz(bzl1,bzl2,bzl3)
                        elseif(k.lt.bzh3) then
                            tbz(i,j,k) = bz(bzl1, bzl2, k)
                        else
                            tbz(i,j,k) = bz(bzl1, bzl2, bzh3)
                        endif
                    elseif(j.lt.bzh2) then
                        if(k.lt.bzl3) then 
                            tbz(i,j,k) = bz(bzl1,j,bzl3)
                        elseif(k.lt.bzh3) then
                            tbz(i,j,k) = bz(bzl1, j, k)
                        else
                            tbz(i,j,k) = bz(bzl1, j, bzh3)
                        endif
                    else
                        if(k.lt.bzl3) then 
                            tbz(i,j,k) = bz(bzl1,bzh2,bzl3)
                        elseif(k.lt.bzh3) then
                            tbz(i,j,k) = bz(bzl1, bzh2, k)
                        else
                            tbz(i,j,k) = bz(bzl1, bzh2, bzh3)
                        endif
                    endif
                elseif(i.lt.bzh1) then 
                    if(j.lt.bzl2) then
                        if(k.lt.bzl3) then 
                            tbz(i,j,k) = bz(i,bzl2,bzl3)
                        elseif(k.lt.bzh3) then
                            tbz(i,j,k) = bz(i, bzl2, k)
                        else
                            tbz(i,j,k) = bz(i, bzl2, bzh3)
                        endif
                    elseif(j.lt.bzh2) then
                        if(k.lt.bzl3) then 
                            tbz(i,j,k) = bz(i,j,bzl3)
                        elseif(k.lt.bzh3) then
                            tbz(i,j,k) = bz(i, j, k)
                        else
                            tbz(i,j,k) = bz(i, j, bzh3)
                        endif
                    else
                        if(k.lt.bzl3) then 
                            tbz(i,j,k) = bz(i,bzh2,bzl3)
                        elseif(k.lt.bzh3) then
                            tbz(i,j,k) = bz(i, bzh2, k)
                        else
                            tbz(i,j,k) = bz(i, bzh2, bzh3)
                        endif
                    endif
                else
                    if(j.lt.bzl2) then
                        if(k.lt.bzl3) then 
                            tbz(i,j,k) = bz(bzh1,bzl2,bzl3)
                        elseif(k.lt.bzh3) then
                            tbz(i,j,k) = bz(bzh1, bzl2, k)
                        else
                            tbz(i,j,k) = bz(bzh1, bzl2, bzh3)
                        endif
                    elseif(j.lt.bzh2) then
                        if(k.lt.bzl3) then 
                            tbz(i,j,k) = bz(bzh1,j,bzl3)
                        elseif(k.lt.bzh3) then
                            tbz(i,j,k) = bz(bzh1, j, k)
                        else
                            tbz(i,j,k) = bz(bzh1, j, bzh3)
                        endif
                    else
                        if(k.lt.bzl3) then 
                            tbz(i,j,k) = bz(bzh1,bzh2,bzl3)
                        elseif(k.lt.bzh3) then
                            tbz(i,j,k) = bz(bzh1, bzh2, k)
                        else
                            tbz(i,j,k) = bz(bzh1, bzh2, bzh3)
                        endif
                    endif
                endif
            enddo
        enddo
    enddo

    !============================================== PLM ================================================================
    do k = s_l3, s_h3
        do j = s_l2, s_h2
            do i = s_l1, s_h1
    
    !============================================ X Direction ==============================================
            summ = 0.d0
            smhd = 0.d0
            dql = 0.d0
            dqr = 0.d0
            dw = 0.d0
            reig = 0.d0
            leig = 0.d0
            lam = 0.d0
            !Skip Bx
            dql(1:5) =  temp(i,j,k,1:ibx-1) - temp(i-1,j,k,1:ibx-1) !gas
            dql(6:7) =  temp(i,j,k,ibx+1:8) - temp(i-1,j,k,ibx+1:8) !mag            
            dqr(1:5) =  temp(i+1,j,k,1:ibx-1) - temp(i,j,k,1:ibx-1)
            dqr(6:7) =  temp(i+1,j,k,ibx+1:8) - temp(i,j,k,ibx+1:8)             
            call evals(lam, s(i,j,k,:), 1) !!X dir eigenvalues
            call lvecx(leig,s(i,j,k,:))    !!left eigenvectors
            call rvecx(reig,s(i,j,k,:))    !!right eigenvectors
    !MHD Source Terms 
            smhd(2) = temp(i,j,k,ibx)/temp(i,j,k,1)
            smhd(3) = temp(i,j,k,iby)/temp(i,j,k,1)
            smhd(4) = temp(i,j,k,ibz)/temp(i,j,k,1)
            smhd(5) = temp(i,j,k,ibx)*temp(i,j,k,2) + temp(i,j,k,iby)*temp(i,j,k,3) + temp(i,j,k,ibz)*temp(i,j,k,4)
            smhd(6) = temp(i,j,k,3)
            smhd(7) = temp(i,j,k,4)
            smhd    = smhd*(tbx(i+1,j,k) - tbx(i,j,k))/dx !normal magnetic field direction
    !Interpolate
        !Plus
                !!Using HLLD so sum over all eigenvalues
            do ii = 1,7
                    dl = dot_product(leig(ii,:),dql)
                    dr = dot_product(leig(ii,:),dqr)
                    call vanleer(dw,dl,dr)
                    summ(:) = summ(:) + (1.d0 - dt_over_a/dx*lam(ii))*dw*reig(:,ii)
            enddo
            ip(i,j,k,qrho:qpres,1)   = temp(i,j,k,1:ibx-1) + 0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5)
            ip(i,j,k,qmagx,1)        = temp(i+1,j,k,ibx) !! Bx stuff
            ip(i,j,k,qmagy:qmagz,1)  = temp(i,j,k,iby:ibz) + 0.5d0*summ(6:7) + 0.5d0*dt_over_a*smhd(6:7)
        !Minus
            summ = 0.d0
            do ii = 1,7
                    dl = dot_product(leig(ii,:),dql)
                    dr = dot_product(leig(ii,:),dqr)
                    call vanleer(dw,dl,dr)          
                    summ(:) = summ(:) + (-1.d0 - dt_over_a/dx*lam(ii))*dw*reig(:,ii)
            enddo
            im(i,j,k,qrho:qpres,1)   = temp(i,j,k,1:ibx-1) + 0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5)
            im(i,j,k,qmagx,1)        = temp(i-1,j,k,ibx) !! Bx stuff
            im(i,j,k,qmagy:qmagz,1)  = temp(i,j,k,iby:ibz) + 0.5d0*summ(6:7) + 0.5d0*dt_over_a*smhd(6:7)


    
    !========================================= Y Direction ================================================             

            summ = 0.d0
            smhd = 0.d0
            dql = 0.d0
            dqr = 0.d0
            dw = 0.d0
            reig = 0.d0
            leig = 0.d0
            lam = 0.d0
            !Skip By
            dql(1:6) =  temp(i,j,k,1:ibx) - temp(i,j-1,k,1:ibx) !gas + bx
            dql(7) =    temp(i,j,k,8) - temp(i,j-1,k,iby+1)     !bz         
            dqr(1:6) =  temp(i,j+1,k,1:ibx) - temp(i,j,k,1:ibx)
            dqr(7) =    temp(i,j+1,k,ibz) - temp(i,j,k,ibz)             
            call evals(lam, s(i,j,k,:), 2) !!Y dir eigenvalues
            call lvecy(leig,s(i,j,k,:))    !!left eigenvectors
            call rvecy(reig,s(i,j,k,:))    !!right eigenvectors
            !!Using HLLD so sum over all eigenvalues
            do ii = 1,7
                    dl = dot_product(leig(ii,:),dql)
                    dr = dot_product(leig(ii,:),dqr)
                    call vanleer(dw,dl,dr)          
                    summ(:) = summ(:) + (1 - dt_over_a/dx*lam(ii))*dw*reig(:,ii)
            enddo
    !MHD Source Terms 
            smhd(2) = temp(i,j,k,ibx)/temp(i,j,k,1)
            smhd(3) = temp(i,j,k,iby)/temp(i,j,k,1)
            smhd(4) = temp(i,j,k,ibz)/temp(i,j,k,1)
            smhd(5) = temp(i,j,k,ibx)*temp(i,j,k,2) + temp(i,j,k,iby)*temp(i,j,k,3) + temp(i,j,k,ibz)*temp(i,j,k,4)
            smhd(6) = temp(i,j,k,2)
            smhd(7) = temp(i,j,k,4)
            smhd    = smhd*(tby(i,j+1,k) - tby(i,j,k))/dy !cross-talk of normal magnetic field direction
    !Interpolate
            ip(i,j,k,qrho:qpres,2)  = temp(i,j,k,1:ibx-1) +0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5) !!GAS
            ip(i,j,k,qmagx,2)       = temp(i,j,k,ibx) + 0.5d0*summ(6) + 0.5d0*dt_over_a*smhd(6)
            ip(i,j,k,qmagy,2)       = temp(i,j+1,k,iby) !! By stuff
            ip(i,j,k,qmagz,2)       = temp(i,j,k,ibz) + 0.5d0*summ(7) + 0.5d0*dt_over_a*smhd(7)

            summ = 0.d0
            do ii = 1,7
                    dl = dot_product(leig(ii,:),dql)
                    dr = dot_product(leig(ii,:),dqr)
                    call vanleer(dw,dl,dr)          
                    summ(:) = summ(:) + (- 1 - dt_over_a/dx*lam(ii))*dw*reig(:,ii)
            enddo
            im(i,j,k,qrho:qpres,2)  = temp(i,j,k,1:ibx-1) + 0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5) !!GAS
            im(i,j,k,qmagx,2)       = temp(i,j,k,ibx) + 0.5d0*summ(6) + 0.5d0*dt_over_a*smhd(6)
            im(i,j,k,qmagy,2)       = temp(i,j-1,k,iby) !! By stuff
            im(i,j,k,qmagz,2)       = temp(i,j,k,ibz) + 0.5d0*summ(7) + 0.5d0*dt_over_a*smhd(7)
            
    !========================================= Z Direction ================================================             
            summ = 0.d0
            smhd = 0.d0
            dql = 0.d0
            dqr = 0.d0
            dw = 0.d0
            reig = 0.d0
            leig = 0.d0
            lam = 0.d0
            !Skip Bz
            dql(1:7) =  temp(i,j,k,1:iby) - temp(i,j,k-1,1:iby) 
            dqr(1:7) =  temp(i,j,k+1,1:iby) - temp(i,j,k,1:iby)
            call evals(lam, s(i,j,k,:), 3) !!Z dir eigenvalues
            call lvecz(leig,s(i,j,k,:))    !!left eigenvectors
            call rvecz(reig,s(i,j,k,:))    !!right eigenvectors

            !!Characteristic Tracing
            do ii = 1,7
                    dl = dot_product(leig(ii,:),dql)
                    dr = dot_product(leig(ii,:),dqr)
                    call vanleer(dw,dl,dr)          
                summ(:) = summ(:) + (1 - dt_over_a/dx*lam(ii))*dw*reig(:,ii)
            enddo
    !MHD Source Terms 
            smhd(2) = temp(i,j,k,ibx)/temp(i,j,k,1)
            smhd(3) = temp(i,j,k,iby)/temp(i,j,k,1)
            smhd(4) = temp(i,j,k,ibz)/temp(i,j,k,1)
            smhd(5) = temp(i,j,k,ibx)*temp(i,j,k,2) + temp(i,j,k,iby)*temp(i,j,k,3) + temp(i,j,k,ibz)*temp(i,j,k,4)
            smhd(6) = temp(i,j,k,2)
            smhd(7) = temp(i,j,k,3)
            smhd    = smhd*(tbz(i,j,k+1) - tbz(i,j,k))/dz !cross-talk of normal magnetic field direction

    !Interpolate
            ip(i,j,k,qrho:qpres,3)  = temp(i,j,k,1:ibx-1) + 0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5) !!GAS
            ip(i,j,k,qmagx:qmagy,3) = temp(i,j,k,ibx:iby) + 0.5d0*summ(6:7) + 0.5d0*dt_over_a*smhd(6:7)
            ip(i,j,k,qmagz,3)       = temp(i,j,k+1,ibz) !! Bz stuff
            summ = 0.d0
            do ii = 1,7
                    dl = dot_product(leig(ii,:),dql)
                    dr = dot_product(leig(ii,:),dqr)
                    call vanleer(dw,dl,dr)          
                summ(:) = summ(:) + (- 1 - dt_over_a/dx*lam(ii))*dw*reig(:,ii)
            enddo
            im(i,j,k,qrho:qpres,3)  = temp(i,j,k,1:ibx-1) + 0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5) !!GAS
            im(i,j,k,qmagx:qmagy,3) = temp(i,j,k,ibx:iby) + 0.5d0*summ(6:7) + 0.5d0*dt_over_a*smhd(6:7)
            im(i,j,k,qmagz,3)       = temp(i,j,k-1,ibz) !! Bz stuff
            im(i,j,k,qreint,1) = (im(i,j,k,qpres,1))/gamma_minus_1
            im(i,j,k,qreint,2) = (im(i,j,k,qpres,2))/gamma_minus_1
            im(i,j,k,qreint,3) = (im(i,j,k,qpres,3))/gamma_minus_1

            ip(i,j,k,qreint,1) =  (ip(i,j,k,qpres,1))/gamma_minus_1
            ip(i,j,k,qreint,2) =  (ip(i,j,k,qpres,2))/gamma_minus_1
            ip(i,j,k,qreint,3) =  (ip(i,j,k,qpres,3))/gamma_minus_1
        enddo
        enddo
    enddo
!Need to add source terms, heating cooling, gravity, etc.
    end subroutine plm
!======================================== Minmod TVD slope limiter =========================================
    subroutine minmod(dW, WR, WL)
    
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    real(rt), intent(in) :: WR, WL
    real(rt), intent(out) :: dW
    dw = 0.d0

    if(abs(wr).lt.abs(wl).and.wr*wl.gt. 0.d0) then
        dw = wr
    elseif(abs(wl).lt.abs(wr).and.wr*wl.gt.0.d0) then
        dw = wl
    endif

    end subroutine


!========================================= VanLeer TVD slope limiter =======================================
    subroutine vanleer(dW, WR, WL) 
    
    use amrex_fort_module, only : rt => amrex_real
    
    implicit none
    
    real(rt), intent(in )   ::  WR, WL
    real(rt), intent(out)   ::  dW
    dw = 0.d0   
    
    if( wr*wl .gt. 0.d0 ) then 
    dw = 2.d0*wr*wl/(wr + wl)
    endif 

    end subroutine


!=========================================== Evals =========================================================

    subroutine evals(lam, Q, dir)

    use amrex_fort_module, only : rt => amrex_real

    implicit none
    
    real(rt), intent(in)    :: Q(QVAR)
    real(rt), intent(out)   :: lam(7) !7 waves
    integer, intent(in)     :: dir !Choose direction, 1 for x, 2 for y, 3 for z

    !The characteristic speeds of the system 
    real(rt)                :: cfx, cfy, cfz, cax, cay, caz, csx, csy, csz, ca, as

    !Speeeeeeeedssssss
    as = gamma_const * q(qpres)/q(qrho)
    !Alfven
    ca  = (q(qmagx)**2 + q(qmagy)**2 + q(qmagz)**2)/q(qrho)
    cax = (q(qmagx)**2)/q(qrho)
    cay = (q(qmagy)**2)/q(qrho)
    caz = (q(qmagz)**2)/q(qrho)
    !Sloooooooooow
    csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))
    csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
    csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))
    !Fassssst
    cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))
    cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
    cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))
    if(dir.eq.1) then   
        !Ax eigenvalues
        lam(1) = q(qu) - sqrt(cfx)
        lam(2) = q(qu) - sqrt(cax)
        lam(3) = q(qu) - sqrt(csx)
        lam(4) = q(qu)
        lam(5) = q(qu) + sqrt(csx)
        lam(6) = q(qu) + sqrt(cax)
        lam(7) = q(qu) + sqrt(cfx)
    elseif(dir.eq.2) then
        !Ay eigenvalues
        lam(1) = q(qv) - sqrt(cfy)
        lam(2) = q(qv) - sqrt(cay)
        lam(3) = q(qv) - sqrt(csy)
        lam(4) = q(qv)
        lam(5) = q(qv) + sqrt(csy)
        lam(6) = q(qv) + sqrt(cay)
        lam(7) = q(qv) + sqrt(cfy)
    else
        !Az eigenvalues
        lam(1) = q(qw) - sqrt(cfz)
        lam(2) = q(qw) - sqrt(caz)
        lam(3) = q(qw) - sqrt(csz)
        lam(4) = q(qw)
        lam(5) = q(qw) + sqrt(csz)
        lam(6) = q(qw) + sqrt(caz)
        lam(7) = q(qw) + sqrt(cfz)
    endif
    end subroutine evals
    
!====================================== Left Eigenvectors ===============================================

!x direction
  subroutine lvecx(leig, Q) 
  use amrex_fort_module, only : rt => amrex_real
  
  implicit none
  
  !returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ax
    real(rt), intent(in ) ::Q(QVAR)
    real(rt), intent(out) ::leig(7,7)

  !The characteristic speeds of the system 
    real(rt)        :: cfx, cax, csx, ca, as, S, N
    real(rt)        :: cff, css, Qf, Qs, AAf, AAs, alf, als, bety, betz

  !Speeeeeeeedssssss
  as = gamma_const * q(qpres)/q(qrho)
  !Alfven
  ca = (q(qmagx)**2 + q(qmagy)**2 + q(qmagz)**2)/q(qrho)
  cax = (q(qmagx)**2)/q(qrho)
  !Sloooooooooow
  csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))
  !Fassssst
  cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))
  !useful constants
  alf = sqrt((as - csx)/(cfx - csx))
  als = sqrt((cfx - as)/(cfx - csx))
  if(cfx - as .lt. 0.d0) als = 0.d0
  if(as - csx .lt. 0.d0) alf = 0.d0
  if(abs(q(qmagy)).le. 1.d-14 .and.abs(q(qmagz)).le. 1.d-14) then
    bety = 1.d0/sqrt(2.d0)
    betz = bety
  else
    bety = q(qmagy)/(sqrt(q(qmagy)**2 + q(qmagz)**2))
    betz = q(qmagz)/(sqrt(q(qmagy)**2 + q(qmagz)**2))
  endif
  cff = sqrt(cfx)*alf
  css = sqrt(csx)*als
  s = sign(1.0d0, q(qmagx))
  qf = sqrt(cfx)*alf*s
  qs = sqrt(csx)*als*s
  n = 0.5d0/as
  aaf = sqrt(as)*alf*sqrt(q(qrho))
  aas = sqrt(as)*als*sqrt(q(qrho))

  
  leig(1,:) = (/0.d0, -n*cff, n*qs*bety   , n*qs*betz   , n*alf/q(qrho) , n*aas*bety/q(qrho)            , n*aas*betz/q(qrho)            /) !u - cf
  leig(2,:) = (/0.d0,  0.d0 , -0.5d0*betz , 0.5d0*bety  , 0.d0          , -0.5d0*s*betz/(sqrt(q(qrho))) , 0.5d0*bety*s/(sqrt(q(qrho)))  /) !u - cAx
  leig(3,:) = (/0.d0, -n*css, -n*qf*bety  , -n*qf*betz  , n*als/q(qrho) , -n*aaf*bety/q(qrho)           , -n*aaf*betz/q(qrho)           /) !u - cs
  leig(4,:) = (/1.d0,  0.d0 ,  0.d0       , 0.d0        , -1.d0/as      , 0.d0                          , 0.d0                          /) !u 
  leig(5,:) = (/0.d0,  n*css, n*qf*bety   , n*qf*betz   , n*als/q(qrho) , -n*aaf*bety/q(qrho)           , -n*aaf*betz/q(qrho)           /) !u + cs
  leig(6,:) = (/0.d0,  0.d0 , 0.5d0*betz  , -0.5d0*bety , 0.d0          , -0.5d0*betz*s/(sqrt(q(qrho))) , 0.5d0*bety*s/(sqrt(q(qrho)))  /) !u + cAx
  leig(7,:) = (/0.d0,  n*cff, -n*qs*bety  , -n*qs*betz  , n*alf/q(qrho) , n*aas*bety/q(qrho)            , n*aas*betz/q(qrho)            /) !u + cf
 
  end subroutine lvecx

!y direction
  subroutine lvecy(leig, Q) 
  use amrex_fort_module, only : rt => amrex_real

  implicit none
  
!returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ay
    real(rt), intent(in ) ::Q(QVAR)
    real(rt), intent(out) ::leig(7,7)

  !The characteristic speeds of the system 
    real(rt)        :: cfy, cay, csy, ca, as, S, N
    real(rt)        :: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, betz

  !Speeeeeeeedssssss
  as = gamma_const * q(qpres)/q(qrho)
  !Alfven
  ca = (q(qmagx)**2 + q(qmagy)**2 + q(qmagz)**2)/q(qrho)
  cay = (q(qmagy)**2)/q(qrho)
  !Sloooooooooow
  csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
  !Fassssst
  cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
  !useful constants
  alf = sqrt((as - csy)/(cfy - csy))
  if(as - csy .lt. 0.d0) alf = 0.d0
  als = sqrt((cfy - as)/(cfy - csy))
  if(cfy - as .lt. 0.d0) als = 0.d0
  if(abs(q(qmagx)).le. 1.d-14 .and.abs(q(qmagz)).le. 1.d-14) then
    betx = 1.d0/sqrt(2.d0)
    betz = betx
  else
    betx = q(qmagx)/(sqrt(q(qmagx)**2 + q(qmagz)**2))
    betz = q(qmagz)/(sqrt(q(qmagx)**2 + q(qmagz)**2))
  endif
  cff = sqrt(cfy)*alf
  css = sqrt(csy)*als
  s = sign(1.0d0, q(qmagy))
  qf = sqrt(cfy)*alf*s
  qs = sqrt(csy)*als*s
  aaf = sqrt(as)*alf*sqrt(q(qrho))
  aas = sqrt(as)*als*sqrt(q(qrho))
  n = 0.5d0/as

  leig(1,:) = (/0.d0, n*qs*betx  , -n*cff   , n*qs*betz   , n*alf/q(qrho) , n*aas*betx/q(qrho)            , n*aas*betz/q(qrho)            /) ! v - cf
  leig(2,:) = (/0.d0, -0.5d0*betz,  0.d0    , 0.5d0*betx  , 0.d0          , -0.5d0*betz*s/(sqrt(q(qrho))) , 0.5d0*betx*s/(sqrt(q(qrho)))  /) ! v - cAy
  leig(3,:) = (/0.d0, -n*qf*betx , -n*css   , -n*qf*betz  , n*als/q(qrho) , -n*aaf*betx/q(qrho)           , -n*aaf*betz/q(qrho)           /) ! v - cs
  leig(4,:) = (/1.d0,  0.d0      ,  0.d0    , 0.d0        , -1.d0/as      , 0.d0                          , 0.d0                          /) ! v 
  leig(5,:) = (/0.d0, n*qf*betx  ,  n*css   , n*qf*betz   , n*als/q(qrho) , -n*aaf*betx/q(qrho)           , -n*aaf*betz/q(qrho)           /) ! v + cs
  leig(6,:) = (/0.d0, 0.5d0*betz ,  0.d0    , -0.5d0*betx , 0.d0          , -0.5d0*betz*s/(sqrt(q(qrho))) , 0.5d0*betx*s/(sqrt(q(qrho)))  /) ! v + cAy
  leig(7,:) = (/0.d0, -n*qs*betx ,  n*cff   , -n*qs*betz  , n*alf/q(qrho) , n*aas*betx/q(qrho)            , n*aas*betz/q(qrho)            /) ! v + cf


  end subroutine lvecy

!z direction
  subroutine lvecz(leig, Q) 
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  !returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Az
    real(rt), intent(in ) ::Q(QVAR)
    real(rt), intent(out) ::leig(7,7)

  !The characteristic speeds of the system 
    real(rt)              :: cfz, caz, csz, ca, as, S, N
    real(rt)              :: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, bety

  !Speeeeeeeedssssss
  as = gamma_const * q(qpres)/q(qrho)
  !Alfven
  ca = (q(qmagx)**2 + q(qmagy)**2 + q(qmagz)**2)/q(qrho)
  caz = (q(qmagz)**2)/q(qrho)
  !Sloooooooooow
  csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))
  !Fassssst
  cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))
  !useful constants
  alf = sqrt((as - csz)/(cfz - csz))
  als = sqrt((cfz - as)/(cfz - csz))
  if(cfz - as .lt. 0.d0) als = 0.d0
  if(abs(q(qmagx)).le. 1.d-14 .and.abs(q(qmagy)).le. 1.d-14) then
    betx = 1.d0/sqrt(2.d0)
    bety = betx
  else
    betx = q(qmagx)/(sqrt(q(qmagx)**2 + q(qmagy)**2))
    bety = q(qmagy)/(sqrt(q(qmagx)**2 + q(qmagy)**2))
  endif
  cff = sqrt(cfz)*alf
  css = sqrt(csz)*als
  s = sign(1.0d0, q(qmagz))
  qf = sqrt(cfz)*alf*s
  qs = sqrt(csz)*als*s
  aaf = sqrt(as)*alf*sqrt(q(qrho))
  aas = sqrt(as)*als*sqrt(q(qrho))
  n = 0.5d0/as

  leig(1,:) = (/0.d0, n*qs*bety  , n*qs*betx   , -n*cff , n*alf/q(qrho) , n*aas*betx/q(qrho)           , n*aas*bety/q(qrho)          /) !w - cf
  leig(2,:) = (/0.d0, 0.5d0*betx , -0.5d0*bety , 0.d0   , 0.d0          , -0.5d0*s*bety/(sqrt(q(qrho))), 0.5d0*betx*s/(sqrt(q(qrho)))/) !w - cAz
  leig(3,:) = (/0.d0, -n*qf*bety , -n*qf*betx  , -n*css , n*als/q(qrho) , -n*aaf*betx/q(qrho)          , -n*aaf*bety/q(qrho)         /) !w - cs
  leig(4,:) = (/1.d0, 0.d0       ,  0.d0       , 0.d0   , -1.d0/as      , 0.d0                         , 0.d0                        /) !w
  leig(5,:) = (/0.d0, n*qf*bety  , n*qf*betx   , n*css  , n*als/q(qrho) , -n*aaf*betx/q(qrho)          , -n*aaf*bety/q(qrho)         /) !w + cs
  leig(6,:) = (/0.d0, -0.5d0*betx, 0.5d0*bety  , 0.d0   , 0.d0          , -0.5d0*bety*s/(sqrt(q(qrho))), 0.5d0*betx*s/(sqrt(q(qrho)))/) !w + cAz
  leig(7,:) = (/0.d0, n*qs*bety  , -n*qs*betx  ,  n*cff , n*alf/q(qrho) , n*aas*betx/q(qrho)           , n*aas*bety/q(qrho)          /) !w + cf
  end subroutine lvecz

!====================================== Right Eigenvectors ===============================================
!x direction
  subroutine rvecx(reig, Q) 
  use amrex_fort_module, only : rt => amrex_real
  
  implicit none
  
  !returnes reig, where the cols are the right eigenvectors of the characteristic matrix Ax
    real(rt), intent(in ) ::Q(QVAR)
    real(rt), intent(out) ::reig(7,7)

  !The characteristic speeds of the system 
    real(rt)        :: cfx, cax, csx, ca, as, S
    real(rt)        :: cff, css, Qf, Qs, AAf, AAs, alf, als, bety, betz

  !Speeeeeeeedssssss
  as = gamma_const * q(qpres)/q(qrho)
  !Alfven
  ca = (q(qmagx)**2 + q(qmagy)**2 + q(qmagz)**2)/q(qrho)
  cax = (q(qmagx)**2)/q(qrho)
  !Sloooooooooow
  csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))
  !Fassssst
  cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))
  !useful constants
  alf = sqrt((as - csx)/(cfx - csx))
  als = sqrt((cfx - as)/(cfx - csx))
  if(cfx - as .lt. 0.d0) als = 0.d0
  if(as - csx .lt. 0.d0) alf = 0.d0
  if(abs(q(qmagy)).le. 1.d-14 .and.abs(q(qmagz)).le. 1.d-14) then
    bety = 1.d0/sqrt(2.d0)
    betz = bety
  else
    bety = q(qmagy)/(sqrt(q(qmagy)**2 + q(qmagz)**2))
    betz = q(qmagz)/(sqrt(q(qmagy)**2 + q(qmagz)**2))
  endif
  cff = sqrt(cfx)*alf
  css = sqrt(csx)*als
  s = sign(1.0d0, q(qmagx))
  qf = sqrt(cfx)*alf*s
  qs = sqrt(csx)*als*s
  aaf = sqrt(as)*alf*sqrt(q(qrho))
  aas = sqrt(as)*als*sqrt(q(qrho))
  
               !   u - cf               u - Cax                         u - cs              u       u + cs              u + Cax                         u + cf
  reig(1,:) = (/  q(qrho)*alf   ,  0.d0                , q(qrho)*als   , 1.d0 , q(qrho)*als   ,  0.d0                ,  q(qrho)*alf   /)
  reig(2,:) = (/  -cff          ,  0.d0                , -css          , 0.d0 , css           ,  0.d0                ,  cff           /)
  reig(3,:) = (/  qs*bety       , -betz                , -qf*bety      , 0.d0 , qf*bety       ,  betz                , -qs*bety       /)
  reig(4,:) = (/  qs*betz       ,  bety                , -qf*betz      , 0.d0 , qf*betz       , -bety                , -qs*betz       /)
  reig(5,:) = (/  q(qrho)*as*alf,  0.d0                , q(qrho)*as*als, 0.d0 , q(qrho)*as*als,  0.d0                ,  q(qrho)*as*alf/)
  reig(6,:) = (/  aas*bety      , -betz*s*sqrt(q(qrho)), -aaf*bety     , 0.d0 , -aaf*bety     , -betz*s*sqrt(q(qrho)),  aas*bety      /)
  reig(7,:) = (/  aas*betz      ,  bety*s*sqrt(q(qrho)), -aaf*betz     , 0.d0 , -aaf*betz     ,  bety*s*sqrt(q(qrho)),  aas*betz      /)


  end subroutine rvecx

!y direction
  subroutine rvecy(reig, Q) 
  use amrex_fort_module, only : rt => amrex_real

  implicit none
  
!returnes reig, where the cols are the right eigenvectors of the characteristic matrix Ay
    real(rt), intent(in ) ::Q(QVAR)
    real(rt), intent(out) ::reig(7,7)

!The characteristic speeds of the system 
    real(rt)        :: cfy, cay, csy, ca, as, S
    real(rt)        :: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, betz

!Speeeeeeeedssssss
  as = gamma_const * q(qpres)/q(qrho)
!Alfven
  ca = (q(qmagx)**2 + q(qmagy)**2 + q(qmagz)**2)/q(qrho)
  cay = (q(qmagy)**2)/q(qrho)
!Sloooooooooow
  csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
!Fassssst
  cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
!useful constants
  alf = sqrt((as - csy)/(cfy - csy))
  if(as - csy .lt. 0.d0) alf = 0.d0
  als = sqrt((cfy - as)/(cfy - csy))
  if(cfy - as .lt. 0.d0) als = 0.d0
  if(abs(q(qmagx)).le. 1.d-14 .and.abs(q(qmagz)).le. 1.d-14) then
    betx = 1.d0/sqrt(2.d0)
    betz = betx
  else
    betx = q(qmagx)/(sqrt(q(qmagx)**2 + q(qmagz)**2))
    betz = q(qmagz)/(sqrt(q(qmagx)**2 + q(qmagz)**2))
  endif
  cff = sqrt(cfy)*alf
  css = sqrt(csy)*als
  s = sign(1.0d0, q(qmagy))
  qf = sqrt(cfy)*alf*s
  qs = sqrt(csy)*als*s
  aaf = sqrt(as)*alf*sqrt(q(qrho))
  aas = sqrt(as)*als*sqrt(q(qrho))
  
              !   v - cf                   v - Cay                      v - cs              v         v + cs                 v + Cay                          v + cf
  reig(1,:) = (/  q(qrho)*alf   ,  0.d0                ,  q(qrho)*als   , 1.d0  , q(qrho)*als   ,  0.d0                ,  q(qrho)*alf   /)
  reig(2,:) = (/  qs*betx       , -betz                , -qf*betx       , 0.d0  , qf*betx       ,  betz                , -qs*betx       /)
  reig(3,:) = (/  -cff          ,  0.d0                , -css           , 0.d0  , css           ,  0.d0                ,  cff           /)
  reig(4,:) = (/  qs*betz       ,  betx                , -qf*betz       , 0.d0  , qf*betz       , -betx                , -qs*betz       /)
  reig(5,:) = (/  q(qrho)*as*alf,  0.d0                ,  q(qrho)*as*als, 0.d0  , q(qrho)*as*als,  0.d0                ,  q(qrho)*as*alf/)
  reig(6,:) = (/  aas*betx      , -betz*s*sqrt(q(qrho)), -aaf*betx      , 0.d0  , -aaf*betx     , -betz*s*sqrt(q(qrho)),  aas*betx      /)
  reig(7,:) = (/  aas*betz      ,  betx*s*sqrt(q(qrho)), -aaf*betz      , 0.d0  , -aaf*betz     ,  betx*s*sqrt(q(qrho)),  aas*betz      /)


  end subroutine rvecy

!z direction
  subroutine rvecz(reig, Q) 
  use amrex_fort_module, only : rt => amrex_real
  
  implicit none
  
!returnes reig, where the cols are the right eigenvectors of the characteristic matrix Az
    real(rt), intent(in ) ::Q(QVAR)
    real(rt), intent(out) ::reig(7,7)

!The characteristic speeds of the system 
    real(rt)        :: cfz, caz, csz, ca, as, S
    real(rt)        :: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, bety

!Speeeeeeeedssssss
  as = gamma_const * q(qpres)/q(qrho)
!Alfven
  ca = (q(qmagx)**2 + q(qmagy)**2 + q(qmagz)**2)/q(qrho)
  caz = (q(qmagz)**2)/q(qrho)
!Sloooooooooow
  csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))
!Fassssst
  cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))
!useful constants
  alf = sqrt((as - csz)/(cfz - csz))
  als = sqrt((cfz - as)/(cfz - csz))
  if(cfz - as .lt. 0.d0) als = 0.d0
  if(abs(q(qmagx)).le. 1.d-14 .and.abs(q(qmagy)).le. 1.d-14) then
    betx = 1.d0/sqrt(2.d0)
    bety = betx
  else
    betx = q(qmagx)/(sqrt(q(qmagx)**2 + q(qmagy)**2))
    bety = q(qmagy)/(sqrt(q(qmagx)**2 + q(qmagy)**2))
  endif
  cff = sqrt(cfz)*alf
  css = sqrt(csz)*als
  s = sign(1.0d0, q(qmagz))
  qf = sqrt(cfz)*alf*s
  qs = sqrt(csz)*als*s
  aaf = sqrt(as)*alf*sqrt(q(qrho))
  aas = sqrt(as)*als*sqrt(q(qrho))
              !   w - cf                   w - Caz                      w - cs              w        w + cs               w + Caz                          w + cf
  reig(1,:) = (/  q(qrho)*alf   ,  0.d0                ,  q(qrho)*als   , 1.d0,  q(qrho)*als   ,  0.d0                ,  q(qrho)*alf   /)
  reig(2,:) = (/  qs*bety       ,  betx                , -qf*bety       , 0.d0,  qf*bety       , -betx                , -qs*bety       /)
  reig(3,:) = (/  qs*betx       , -bety                , -qf*betx       , 0.d0,  qf*betx       ,  bety                , -qs*betx       /)
  reig(4,:) = (/  -cff          ,  0.d0                , -css           , 0.d0,  css           ,  0.d0                ,  cff           /)
  reig(5,:) = (/  q(qrho)*as*alf,  0.d0                ,  q(qrho)*as*als, 0.d0,  q(qrho)*as*als,  0.d0                ,  q(qrho)*as*alf/)
  reig(6,:) = (/  aas*betx      , -bety*s*sqrt(q(qrho)), -aaf*betx      , 0.d0, -aaf*betx      , -bety*s*sqrt(q(qrho)),  aas*betx      /)
  reig(7,:) = (/  aas*bety      ,  betx*s*sqrt(q(qrho)), -aaf*bety      , 0.d0, -aaf*bety      ,  betx*s*sqrt(q(qrho)),  aas*bety      /)


  end subroutine rvecz
end module mhd_plm_module
````

