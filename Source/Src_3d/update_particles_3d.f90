
  subroutine update_dm_particles(np, particles, accel, accel_lo, accel_hi, &
                                 plo, dx, dt, a_prev, a_cur, do_move) &
       bind(c,name='update_dm_particles')

    use iso_c_binding
    use amrex_error_module
    use amrex_fort_module, only : amrex_real
    use particle_mod      , only: dm_particle_t

    integer,             intent(in   )        :: np
    type(dm_particle_t), intent(inout)        :: particles(np)
    integer,             intent(in   )        :: accel_lo(3), accel_hi(3)
    real(amrex_real),    intent(in   )        :: accel &
         (accel_lo(1):accel_hi(1),accel_lo(2):accel_hi(2),accel_lo(3):accel_hi(3),3)
    real(amrex_real),    intent(in   )        :: plo(3),dx(3),dt,a_prev,a_cur
    integer,             intent(in   )        :: do_move

    integer          :: i, j, k, n, nc
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) half_dt, a_cur_inv, dt_a_cur_inv
    real(amrex_real) inv_dx(3)
    real(amrex_real) part_accel

    inv_dx = 1.0d0/dx

    half_dt       = 0.5d0 * dt
    a_cur_inv    = 1.d0 / a_cur;
    dt_a_cur_inv = dt * a_cur_inv;

    do n = 1, np

       if (particles(n)%id .le. 0) then
          cycle
       end if

       lx = (particles(n)%pos(1) - plo(1))*inv_dx(1) + 0.5d0
       ly = (particles(n)%pos(2) - plo(2))*inv_dx(2) + 0.5d0
       lz = (particles(n)%pos(3) - plo(3))*inv_dx(3) + 0.5d0

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       if (i-1 .lt. accel_lo(1) .or. i .gt. accel_hi(1) .or. &
           j-1 .lt. accel_lo(2) .or. j .gt. accel_hi(2) .or. &
           k-1 .lt. accel_lo(3) .or. k .gt. accel_hi(3)) then
          print *,'PARTICLE ID ', particles(n)%id,' REACHING OUT OF BOUNDS AT (I,J,K) = ',i,j,k
          call amrex_error('Aborting in move_kick_drift')
       end if

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       ! moveKickDrift: update velocity by half dt:  (a u)^half = (a u)^old  + dt/2 grav^old
       ! moveKick    t: update velocity by half dt:  (a u)^new  = (a u)^half + dt/2 grav^new
       do nc = 1, 3

          part_accel = &
                wx_lo*wy_lo*wz_lo*accel(i-1, j-1, k-1, nc) + &
                wx_lo*wy_lo*wz_hi*accel(i-1, j-1, k  , nc) + &
                wx_lo*wy_hi*wz_lo*accel(i-1, j,   k-1, nc) + &
                wx_lo*wy_hi*wz_hi*accel(i-1, j,   k  , nc) + &
                wx_hi*wy_lo*wz_lo*accel(i,   j-1, k-1, nc) + &
                wx_hi*wy_lo*wz_hi*accel(i,   j-1, k  , nc) + &
                wx_hi*wy_hi*wz_lo*accel(i,   j,   k-1, nc) + &
                wx_hi*wy_hi*wz_hi*accel(i,   j,   k  , nc)

          particles(n)%vel(nc) = a_prev*particles(n)%vel(nc) + half_dt * part_accel
          particles(n)%vel(nc) = particles(n)%vel(nc) * a_cur_inv

       end do

       ! moveKickDrift: Update position by full dt: x^new = x^old + dt u^half / a^half
       if (do_move .eq. 1) then
          do nc = 1, 3
             particles(n)%pos(nc) = particles(n)%pos(nc) + dt_a_cur_inv * particles(n)%vel(nc)
          end do
       end if

    end do

  end subroutine update_dm_particles

  subroutine update_agn_particles(np, particles, accel, accel_lo, accel_hi, &
                                  plo, dx, dt, a_prev, a_cur, do_move) &
       bind(c,name='update_agn_particles')

    use iso_c_binding
    use amrex_error_module
    use amrex_fort_module, only : amrex_real
    use particle_mod      , only: agn_particle_t

    integer,             intent(in   )        :: np
    type(agn_particle_t), intent(inout)        :: particles(np)
    integer,             intent(in   )        :: accel_lo(3), accel_hi(3)
    real(amrex_real),    intent(in   )        :: accel &
         (accel_lo(1):accel_hi(1),accel_lo(2):accel_hi(2),accel_lo(3):accel_hi(3),3)
    real(amrex_real),    intent(in   )        :: plo(3),dx(3),dt,a_prev,a_cur
    integer,             intent(in   )        :: do_move

    integer          :: i, j, k, n, nc
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) half_dt, a_cur_inv, dt_a_cur_inv
    real(amrex_real) inv_dx(3)
    real(amrex_real) part_accel

    inv_dx = 1.0d0/dx

    half_dt       = 0.5d0 * dt
    a_cur_inv    = 1.d0 / a_cur;
    dt_a_cur_inv = dt * a_cur_inv;

    do n = 1, np

       if (particles(n)%id .le. 0) then
          cycle
       end if

       lx = (particles(n)%pos(1) - plo(1))*inv_dx(1) + 0.5d0
       ly = (particles(n)%pos(2) - plo(2))*inv_dx(2) + 0.5d0
       lz = (particles(n)%pos(3) - plo(3))*inv_dx(3) + 0.5d0

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       if (i-1 .lt. accel_lo(1) .or. i .gt. accel_hi(1) .or. &
           j-1 .lt. accel_lo(2) .or. j .gt. accel_hi(2) .or. &
           k-1 .lt. accel_lo(3) .or. k .gt. accel_hi(3)) then
          print *,'PARTICLE ID ', particles(n)%id,' REACHING OUT OF BOUNDS AT (I,J,K) = ',i,j,k
          call amrex_error('Aborting in move_kick_drift')
       end if

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       ! moveKickDrift: update velocity by half dt:  (a u)^half = (a u)^old  + dt/2 grav^old
       ! moveKick    t: update velocity by half dt:  (a u)^new  = (a u)^half + dt/2 grav^new
       do nc = 1, 3

          part_accel = &
                wx_lo*wy_lo*wz_lo*accel(i-1, j-1, k-1, nc) + &
                wx_lo*wy_lo*wz_hi*accel(i-1, j-1, k  , nc) + &
                wx_lo*wy_hi*wz_lo*accel(i-1, j,   k-1, nc) + &
                wx_lo*wy_hi*wz_hi*accel(i-1, j,   k  , nc) + &
                wx_hi*wy_lo*wz_lo*accel(i,   j-1, k-1, nc) + &
                wx_hi*wy_lo*wz_hi*accel(i,   j-1, k  , nc) + &
                wx_hi*wy_hi*wz_lo*accel(i,   j,   k-1, nc) + &
                wx_hi*wy_hi*wz_hi*accel(i,   j,   k  , nc)

          particles(n)%vel(nc) = a_prev*particles(n)%vel(nc) + half_dt * part_accel
          particles(n)%vel(nc) = particles(n)%vel(nc) * a_cur_inv

       end do

       ! moveKickDrift: Update position by full dt: x^new = x^old + dt u^half / a^half
       if (do_move .eq. 1) then
          do nc = 1, 3
             particles(n)%pos(nc) = particles(n)%pos(nc) + dt_a_cur_inv * particles(n)%vel(nc)
          end do
       end if

    end do

  end subroutine update_agn_particles

  subroutine update_fdm_particles(np, particles, accel, accel_lo, accel_hi, &
                                 plo, dx, dt, a_prev, a_cur, do_move) &
       bind(c,name='update_fdm_particles')

    use iso_c_binding
    use amrex_error_module
    use amrex_fort_module, only : amrex_real
    use particle_mod      , only: fdm_particle_t

    integer,              intent(in   )        :: np
    type(fdm_particle_t), intent(inout)        :: particles(np)
    integer,              intent(in   )        :: accel_lo(3), accel_hi(3)
    real(amrex_real),     intent(in   )        :: accel &
         (accel_lo(1):accel_hi(1),accel_lo(2):accel_hi(2),accel_lo(3):accel_hi(3),3)
    real(amrex_real),     intent(in   )        :: plo(3),dx(3),dt,a_prev,a_cur
    integer,              intent(in   )        :: do_move

    integer          :: i, j, k, n, nc
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) half_dt, a_cur_inv, dt_a_cur_inv
    real(amrex_real) inv_dx(3)
    real(amrex_real) part_accel

    inv_dx = 1.0d0/dx

    half_dt       = 0.5d0 * dt
    a_cur_inv    = 1.d0 / a_cur;
    dt_a_cur_inv = dt * a_cur_inv;

    do n = 1, np

       if (particles(n)%id .le. 0) then
          cycle
       end if

       lx = (particles(n)%pos(1) - plo(1))*inv_dx(1) + 0.5d0
       ly = (particles(n)%pos(2) - plo(2))*inv_dx(2) + 0.5d0
       lz = (particles(n)%pos(3) - plo(3))*inv_dx(3) + 0.5d0

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       if (i-1 .lt. accel_lo(1) .or. i .gt. accel_hi(1) .or. &
           j-1 .lt. accel_lo(2) .or. j .gt. accel_hi(2) .or. &
           k-1 .lt. accel_lo(3) .or. k .gt. accel_hi(3)) then
          print *,'PARTICLE ID ', particles(n)%id,' REACHING OUT OF BOUNDS AT (I,J,K) (accel)= ',i,j,k
          call amrex_error('Aborting in move_kick_drift')
       end if

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       ! moveKickDrift: update velocity by half dt:  (a u)^half = (a u)^old  + dt/2 grav^old
       ! moveKick    t: update velocity by half dt:  (a u)^new  = (a u)^half + dt/2 grav^new
       do nc = 1, 3

          part_accel = &
                wx_lo*wy_lo*wz_lo*accel(i-1, j-1, k-1, nc) + &
                wx_lo*wy_lo*wz_hi*accel(i-1, j-1, k  , nc) + &
                wx_lo*wy_hi*wz_lo*accel(i-1, j,   k-1, nc) + &
                wx_lo*wy_hi*wz_hi*accel(i-1, j,   k  , nc) + &
                wx_hi*wy_lo*wz_lo*accel(i,   j-1, k-1, nc) + &
                wx_hi*wy_lo*wz_hi*accel(i,   j-1, k  , nc) + &
                wx_hi*wy_hi*wz_lo*accel(i,   j,   k-1, nc) + &
                wx_hi*wy_hi*wz_hi*accel(i,   j,   k  , nc)

          particles(n)%vel(nc) = a_prev*particles(n)%vel(nc) + half_dt * part_accel
          particles(n)%vel(nc) = particles(n)%vel(nc) * a_cur_inv

       end do

       ! moveKickDrift: Update position by full dt: x^new = x^old + dt u^half / a^half
       if (do_move .eq. 1) then
          do nc = 1, 3
             particles(n)%pos(nc) = particles(n)%pos(nc) + dt_a_cur_inv * particles(n)%vel(nc)
          end do
       end if

    end do

  end subroutine update_fdm_particles

  subroutine update_gaussian_beams(np, particles, accel, accel_lo, accel_hi, &
                                   phi, phi_lo, phi_hi, &
                                   plo, dx, dt, a_prev, a_cur, do_move) &
       bind(c,name='update_gaussian_beams')

    use iso_c_binding
    use amrex_error_module
    use amrex_fort_module, only : amrex_real
    use particle_mod      , only: fdm_particle_t
    use fdm_params_module, only : hbaroverm, ii, gamma_fdm!, sigma_fdm

    integer,              intent(in   )        :: np
    type(fdm_particle_t), intent(inout)        :: particles(np)
    integer,              intent(in   )        :: accel_lo(3), accel_hi(3)
    ! real(amrex_real),     intent(inout)        :: accel &
    real(amrex_real),     intent(in   )        :: accel &
         (accel_lo(1):accel_hi(1),accel_lo(2):accel_hi(2),accel_lo(3):accel_hi(3),3)
    integer,              intent(in   )        :: phi_lo(3), phi_hi(3)
    ! real(amrex_real),     intent(inout)        :: phi &
    real(amrex_real),     intent(in   )        :: phi &
         (phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3),1)
    real(amrex_real),     intent(in   )        :: plo(3),dx(3),dt,a_prev,a_cur
    integer,              intent(in   )        :: do_move

    integer          :: i, j, k, n, nc, nn, nn_hi
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) half_dt, a_cur_inv, dt_a_cur_inv
    real(amrex_real) inv_dx(3), Arg1, Arg2
    real(amrex_real) part_accel, part_phi, vel_square
    real(amrex_real) qq(3,3), pq(3,3), qqold(3,3), pqold(3,3), hessian(3,3)
    real(amrex_real) qp(3,3), pp(3,3), qpold(3,3), ppold(3,3)
    complex(amrex_real) Cpq, ZZ(3,3), det, Cpqtemp, Cpqold
    logical :: test

    inv_dx = 1.0d0/dx

    half_dt       = 0.5d0 * dt
    a_cur_inv    = 1.d0 / a_cur
    dt_a_cur_inv = dt * a_cur_inv

    ! phi = 0.d0
    ! accel = 0.d0
    
    ! do k = phi_lo(3), phi_hi(3)
    !    do j = phi_lo(2), phi_hi(2)
    !       do i = phi_lo(1), phi_hi(1)
    !          phi(i,j,k,1) = -( ((i+0.5d0)*dx(1)+plo(1))*((i+0.5d0)*dx(1)+plo(1))+ &
    !                          ((j+0.5d0)*dx(2)+plo(2))*((j+0.5d0)*dx(2)+plo(2))+ &
    !                          ((k+0.5d0)*dx(3)+plo(3))*((k+0.5d0)*dx(3)+plo(3)) ) * &
    !                        2.0*hbaroverm*hbaroverm*1600.0*1600.0
    !       enddo
    !    enddo
    ! enddo

    ! do k = accel_lo(3), accel_hi(3)
    !    do j = accel_lo(2), accel_hi(2)
    !       do i = accel_lo(1), accel_hi(1)
    !          accel(i,j,k,1) = -((i+0.5d0)*dx(1)+plo(1)) * 4.0*hbaroverm*hbaroverm*1600.0*1600.0
    !          accel(i,j,k,2) = -((j+0.5d0)*dx(2)+plo(2)) * 4.0*hbaroverm*hbaroverm*1600.0*1600.0
    !          accel(i,j,k,3) = -((k+0.5d0)*dx(3)+plo(3)) * 4.0*hbaroverm*hbaroverm*1600.0*1600.0
    !       enddo
    !    enddo
    ! enddo

    do n = 1, np

       if (particles(n)%id .le. 0) then
          cycle
       end if

!                                                                                                                                                                                                              
!     Construct Gaussian beam parameters                                                                                                                                                                       
!                                                                                                                                                                                                              
       Cpq   = cmplx(particles(n)%amp(1),particles(n)%amp(2))
       qq    = reshape((/ particles(n)%qq(1) , particles(n)%qq(2) , particles(n)%qq(3) , &
                          particles(n)%qq(4) , particles(n)%qq(5) , particles(n)%qq(6) , &
                          particles(n)%qq(7) , particles(n)%qq(8) , particles(n)%qq(9) /), shape(qq))
       pq    = reshape((/ particles(n)%pq(1) , particles(n)%pq(2) , particles(n)%pq(3) , &
                          particles(n)%pq(4) , particles(n)%pq(5) , particles(n)%pq(6) , &
                          particles(n)%pq(7) , particles(n)%pq(8) , particles(n)%pq(9) /), shape(pq))
       qp    = reshape((/ particles(n)%qp(1) , particles(n)%qp(2) , particles(n)%qp(3) , &
                          particles(n)%qp(4) , particles(n)%qp(5) , particles(n)%qp(6) , &
                          particles(n)%qp(7) , particles(n)%qp(8) , particles(n)%qp(9) /), shape(qp))
       pp    = reshape((/ particles(n)%pp(1) , particles(n)%pp(2) , particles(n)%pp(3) , &
                          particles(n)%pp(4) , particles(n)%pp(5) , particles(n)%pp(6) , &
                          particles(n)%pp(7) , particles(n)%pp(8) , particles(n)%pp(9) /), shape(pp))
!
!      Store old data
!
       qqold = qq
       pqold = pq
       qpold = qp
       ppold = pp

       lx = (particles(n)%pos(1) - plo(1))*inv_dx(1) + 0.5d0
       ly = (particles(n)%pos(2) - plo(2))*inv_dx(2) + 0.5d0
       lz = (particles(n)%pos(3) - plo(3))*inv_dx(3) + 0.5d0

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       if (i-1 .lt. accel_lo(1) .or. i .gt. accel_hi(1) .or. &
           j-1 .lt. accel_lo(2) .or. j .gt. accel_hi(2) .or. &
           k-1 .lt. accel_lo(3) .or. k .gt. accel_hi(3)) then
          print *,'PARTICLE ID ', particles(n)%id,' REACHING OUT OF BOUNDS AT (I,J,K) (accel)= ',i,j,k
          call amrex_error('Aborting in move_kick_drift')
       end if

       if (i-1 .lt. phi_lo(1) .or. i .gt. phi_hi(1) .or. &
           j-1 .lt. phi_lo(2) .or. j .gt. phi_hi(2) .or. &
           k-1 .lt. phi_lo(3) .or. k .gt. phi_hi(3)) then
          print *,'PARTICLE ID ', particles(n)%id,' REACHING OUT OF BOUNDS AT (I,J,K) (phi)= ',i,j,k
          call amrex_error('Aborting in move_kick_drift')
       end if

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       do nc = 1, 3

          part_accel = &
                wx_lo*wy_lo*wz_lo*accel(i-1, j-1, k-1, nc) + &
                wx_lo*wy_lo*wz_hi*accel(i-1, j-1, k  , nc) + &
                wx_lo*wy_hi*wz_lo*accel(i-1, j,   k-1, nc) + &
                wx_lo*wy_hi*wz_hi*accel(i-1, j,   k  , nc) + &
                wx_hi*wy_lo*wz_lo*accel(i,   j-1, k-1, nc) + &
                wx_hi*wy_lo*wz_hi*accel(i,   j-1, k  , nc) + &
                wx_hi*wy_hi*wz_lo*accel(i,   j,   k-1, nc) + &
                wx_hi*wy_hi*wz_hi*accel(i,   j,   k  , nc)

          particles(n)%vel(nc) = a_prev*particles(n)%vel(nc) + half_dt * part_accel
          particles(n)%vel(nc) = particles(n)%vel(nc) * a_cur_inv

       end do

       part_phi = &
            wx_lo*wy_lo*wz_lo*phi(i-1, j-1, k-1, 1) + &
                wx_lo*wy_lo*wz_hi*phi(i-1, j-1, k  , 1) + &
                wx_lo*wy_hi*wz_lo*phi(i-1, j,   k-1, 1) + &
                wx_lo*wy_hi*wz_hi*phi(i-1, j,   k  , 1) + &
                wx_hi*wy_lo*wz_lo*phi(i,   j-1, k-1, 1) + &
                wx_hi*wy_lo*wz_hi*phi(i,   j-1, k  , 1) + &
                wx_hi*wy_hi*wz_lo*phi(i,   j,   k-1, 1) + &
                wx_hi*wy_hi*wz_hi*phi(i,   j,   k  , 1)

       !Remember different sign convention for phi
       particles(n)%phase = particles(n)%phase + half_dt * part_phi

       hessian(1,1) = &
               (wy_hi*wz_hi*accel(i,   j,   k  , 1) + &
                wy_hi*wz_lo*accel(i,   j,   k-1, 1) + &
                wy_lo*wz_hi*accel(i,   j-1, k  , 1) + &
                wy_lo*wz_lo*accel(i,   j-1, k-1, 1) - &
                wy_hi*wz_hi*accel(i-1, j,   k  , 1) - &
                wy_hi*wz_lo*accel(i-1, j,   k-1, 1) - &
                wy_lo*wz_hi*accel(i-1, j-1, k  , 1) - &
                wy_lo*wz_lo*accel(i-1, j-1, k-1, 1)) * inv_dx(1)

       hessian(1,2) = &
               (wy_hi*wz_hi*accel(i,   j,   k  , 2) + &
                wy_hi*wz_lo*accel(i,   j,   k-1, 2) + &
                wy_lo*wz_hi*accel(i,   j-1, k  , 2) + &
                wy_lo*wz_lo*accel(i,   j-1, k-1, 2) - &
                wy_hi*wz_hi*accel(i-1, j,   k  , 2) - &
                wy_hi*wz_lo*accel(i-1, j,   k-1, 2) - &
                wy_lo*wz_hi*accel(i-1, j-1, k  , 2) - &
                wy_lo*wz_lo*accel(i-1, j-1, k-1, 2)) * inv_dx(1)

       hessian(1,3) = &
               (wy_hi*wz_hi*accel(i,   j,   k  , 3) + &
                wy_hi*wz_lo*accel(i,   j,   k-1, 3) + &
                wy_lo*wz_hi*accel(i,   j-1, k  , 3) + &
                wy_lo*wz_lo*accel(i,   j-1, k-1, 3) - &
                wy_hi*wz_hi*accel(i-1, j,   k  , 3) - &
                wy_hi*wz_lo*accel(i-1, j,   k-1, 3) - &
                wy_lo*wz_hi*accel(i-1, j-1, k  , 3) - &
                wy_lo*wz_lo*accel(i-1, j-1, k-1, 3)) * inv_dx(1)
       
       hessian(2,1) = hessian(1,2)
       
       hessian(2,2) = &
               (wx_hi*wz_hi*accel(i,   j,   k  , 2) + &
                wx_hi*wz_lo*accel(i,   j,   k-1, 2) + &
                wx_lo*wz_hi*accel(i-1, j,   k  , 2) + &
                wx_lo*wz_lo*accel(i-1, j,   k-1, 2) - &
                wx_hi*wz_hi*accel(i,   j-1, k  , 2) - &
                wx_hi*wz_lo*accel(i,   j-1, k-1, 2) - &
                wx_lo*wz_hi*accel(i-1, j-1, k  , 2) - &
                wx_lo*wz_lo*accel(i-1, j-1, k-1, 2)) * inv_dx(2)

       hessian(2,3) = &
               (wx_hi*wz_hi*accel(i,   j,   k  , 3) + &
                wx_hi*wz_lo*accel(i,   j,   k-1, 3) + &
                wx_lo*wz_hi*accel(i-1, j,   k  , 3) + &
                wx_lo*wz_lo*accel(i-1, j,   k-1, 3) - &
                wx_hi*wz_hi*accel(i,   j-1, k  , 3) - &
                wx_hi*wz_lo*accel(i,   j-1, k-1, 3) - &
                wx_lo*wz_hi*accel(i-1, j-1, k  , 3) - &
                wx_lo*wz_lo*accel(i-1, j-1, k-1, 3)) * inv_dx(2)

       hessian(3,1) = hessian(1,3)

       hessian(3,2) = hessian(2,3)

       hessian(3,3) = &
               (wx_hi*wy_hi*accel(i,   j,   k  , 3) + &
                wx_hi*wy_lo*accel(i,   j-1, k  , 3) + &
                wx_lo*wy_hi*accel(i-1, j,   k  , 3) + &
                wx_lo*wy_lo*accel(i-1, j-1, k  , 3) - &
                wx_hi*wy_hi*accel(i,   j  , k-1, 3) - &
                wx_hi*wy_lo*accel(i,   j-1, k-1, 3) - &
                wx_lo*wy_hi*accel(i-1, j  , k-1, 3) - &
                wx_lo*wy_lo*accel(i-1, j-1, k-1, 3)) * inv_dx(3)

       
       !Remember different sign convention for phi
       pq = pq + half_dt * matmul(hessian,qq)
       pp = pp + half_dt * matmul(hessian,qp)

       if (do_move .eq. 1) then

          do nc = 1, 3
             particles(n)%pos(nc) = particles(n)%pos(nc) + dt_a_cur_inv * particles(n)%vel(nc)
          end do
          
          vel_square = particles(n)%vel(1)*particles(n)%vel(1) + &
                       particles(n)%vel(2)*particles(n)%vel(2) + &
                       particles(n)%vel(3)*particles(n)%vel(3) 

          particles(n)%phase = particles(n)%phase + half_dt * vel_square

          qq = qq + dt_a_cur_inv * a_cur_inv * pq
          qp = qp + dt_a_cur_inv * a_cur_inv * pp
          
       else

!
!         Update complex beam amplitude Cpq
!          
          ZZ = particles(n)%width*qq-pq/(2.0*ii*hbaroverm) &
               + (pp-2.0*ii*hbaroverm*particles(n)%width*qp)*gamma_fdm

          det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
               - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
               + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))

        Cpqtemp = sqrt(det)
        Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
        Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)

        if( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
           test = .false.
        else
           if ( Arg1 .le. Arg2 ) then
              Cpq = Cpqtemp
           else
              Cpq = -Cpqtemp
           endif
           test = .true.
        endif
!                                                                                                                                                                                                                  
        if (test .eqv. .false.) then
           nn_hi = 10
           Cpqold = Cpq
           do nn=1, nn_hi
!                                                                                                                                                                                                                  
              ZZ = particles(n)%width*((qqold*(nn_hi-nn)+nn*qq)/nn_hi)-((pqold*(nn_hi-nn)+nn*pq)/nn_hi)/(2.0*ii*hbaroverm) &
                   + (((ppold*(nn_hi-nn)+nn*pp)/nn_hi)-2.0*ii*hbaroverm*particles(n)%width*((qpold*(nn_hi-nn)+nn*qp)/nn_hi))*gamma_fdm                    
!                                                                                                                                                                                                      
              det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
                   - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
                   + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))
!                                                                                                                                                                                                           
              Cpqtemp = sqrt(det)
              Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
              Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)
              if ( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
                 test = .false.
                 exit
              else
                 if ( Arg1 .le. Arg2 ) then
                    Cpq = Cpqtemp
                 else
                    Cpq = -Cpqtemp
                 endif
                 test = .true.
              endif
           enddo
        endif
!                                                                                                                                                                                                                
        if (test .eqv. .false.) then
           nn_hi = 100
           Cpq = Cpqold
           do nn=1, nn_hi
!                                                                                                                                                                                                                
              ZZ = particles(n)%width*((qqold*(nn_hi-nn)+nn*qq)/nn_hi)-((pqold*(nn_hi-nn)+nn*pq)/nn_hi)/(2.0*ii*hbaroverm) &
                   + (((ppold*(nn_hi-nn)+nn*pp)/nn_hi)-2.0*ii*hbaroverm*particles(n)%width*((qpold*(nn_hi-nn)+nn*qp)/nn_hi))*gamma_fdm                    
!                                                                                                                                                                                                                
              det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
                   - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
                   + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))
!                                                                                                                                                                                                                
              Cpqtemp = sqrt(det)
              Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
              Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)
              if ( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
                 test = .false.
                 exit
              else
                 if ( Arg1 .le. Arg2 ) then
                    Cpq = Cpqtemp
                 else
                    Cpq = -Cpqtemp
                 endif
                 test = .true.
              endif
           enddo
        endif
!                                                                                                                                                                                                                
        if (test .eqv. .false.) then
           nn_hi = 1000
           Cpq = Cpqold
           do nn=1, nn_hi
!                                                                                                                                                                                                                
              ZZ = particles(n)%width*((qqold*(nn_hi-nn)+nn*qq)/nn_hi)-((pqold*(nn_hi-nn)+nn*pq)/nn_hi)/(2.0*ii*hbaroverm) &
                   + (((ppold*(nn_hi-nn)+nn*pp)/nn_hi)-2.0*ii*hbaroverm*particles(n)%width*((qpold*(nn_hi-nn)+nn*qp)/nn_hi))*gamma_fdm                    
!                                                                                                                                                                                                              
              det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
                   - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
                   + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))
!                                                                                                                                                                                                                
              Cpqtemp = sqrt(det)
              Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
              Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)
              if ( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
                 test = .false.
                 exit
              else
                 if ( Arg1 .le. Arg2 ) then
                    Cpq = Cpqtemp
                 else
                    Cpq = -Cpqtemp
                 endif
                 test = .true.
              endif
           enddo
        endif
!                                                                                                                                                                                                                
        if (test .eqv. .false.) then
           nn_hi = 10000
           Cpq = Cpqold
           do nn=1, nn_hi
!                                                                                                                                                                                                                
              ZZ = particles(n)%width*((qqold*(nn_hi-nn)+nn*qq)/nn_hi)-((pqold*(nn_hi-nn)+nn*pq)/nn_hi)/(2.0*ii*hbaroverm) &
                   + (((ppold*(nn_hi-nn)+nn*pp)/nn_hi)-2.0*ii*hbaroverm*particles(n)%width*((qpold*(nn_hi-nn)+nn*qp)/nn_hi))*gamma_fdm
!                                                                                                                                                                                                              
              det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
                   - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
                   + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))
!                                                                                                                                                                                                               
              Cpqtemp = sqrt(det)
              Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
              Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)
              if ( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
                 test = .false.
                 exit
              else
                 if ( Arg1 .le. Arg2 ) then
                    Cpq = Cpqtemp
                 else
                    Cpq = -Cpqtemp
                 endif
                 test = .true.
              endif
           enddo
        endif
!                                                                                                                                                                                                                
        if (test .eqv. .false.) then
           nn_hi = 100000
           Cpq = Cpqold
           do nn=1, nn_hi
!                                                                                                                                                                                                                
              ZZ = particles(n)%width*((qqold*(nn_hi-nn)+nn*qq)/nn_hi)-((pqold*(nn_hi-nn)+nn*pq)/nn_hi)/(2.0*ii*hbaroverm) &
                   + (((ppold*(nn_hi-nn)+nn*pp)/nn_hi)-2.0*ii*hbaroverm*particles(n)%width*((qpold*(nn_hi-nn)+nn*qp)/nn_hi))*gamma_fdm
!                                                                                                                                                                                                                
              det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
                   - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
                   + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))
!                                                                                                                                                                                                                
              Cpqtemp = sqrt(det)
              Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
              Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)
              if ( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
                 print *, "Error: functional determinant is changing too fast! Arguments are: ",Arg1, Arg2, Cpqtemp, Cpq, nn, nn_hi
                 call amrex_error('Aborting in update_gaussian_beams')
              else
                 if ( Arg1 .le. Arg2 ) then
                    Cpq = Cpqtemp
                 else
                    Cpq = -Cpqtemp
                 endif
              endif
           enddo
        endif

        particles(n)%amp(1) = real(real(Cpq))
        particles(n)%amp(2) = real(aimag(Cpq))

        ! print *, hbaroverm, ii, gamma_fdm, particles(n)%amp(1), particles(n)%amp(2)

       end if

       particles(n)%qq(1) = qq(1,1)
       particles(n)%qq(2) = qq(2,1)
       particles(n)%qq(3) = qq(3,1)
       particles(n)%qq(4) = qq(1,2)
       particles(n)%qq(5) = qq(2,2)
       particles(n)%qq(6) = qq(3,2)
       particles(n)%qq(7) = qq(1,3)
       particles(n)%qq(8) = qq(2,3)
       particles(n)%qq(9) = qq(3,3)

       particles(n)%pq(1) = pq(1,1)
       particles(n)%pq(2) = pq(2,1)
       particles(n)%pq(3) = pq(3,1)
       particles(n)%pq(4) = pq(1,2)
       particles(n)%pq(5) = pq(2,2)
       particles(n)%pq(6) = pq(3,2)
       particles(n)%pq(7) = pq(1,3)
       particles(n)%pq(8) = pq(2,3)
       particles(n)%pq(9) = pq(3,3)

       particles(n)%qp(1) = qp(1,1)
       particles(n)%qp(2) = qp(2,1)
       particles(n)%qp(3) = qp(3,1)
       particles(n)%qp(4) = qp(1,2)
       particles(n)%qp(5) = qp(2,2)
       particles(n)%qp(6) = qp(3,2)
       particles(n)%qp(7) = qp(1,3)
       particles(n)%qp(8) = qp(2,3)
       particles(n)%qp(9) = qp(3,3)

       particles(n)%pp(1) = pp(1,1)
       particles(n)%pp(2) = pp(2,1)
       particles(n)%pp(3) = pp(3,1)
       particles(n)%pp(4) = pp(1,2)
       particles(n)%pp(5) = pp(2,2)
       particles(n)%pp(6) = pp(3,2)
       particles(n)%pp(7) = pp(1,3)
       particles(n)%pp(8) = pp(2,3)
       particles(n)%pp(9) = pp(3,3)

    end do

  end subroutine update_gaussian_beams

  subroutine update_fdm_particles_wkb(np, particles, accel, accel_lo, accel_hi, &
                                      plo, dx, dt, a_prev, a_cur, do_move) &
       bind(c,name='update_fdm_particles_wkb')

    use iso_c_binding
    use amrex_error_module
    use amrex_fort_module, only : amrex_real
    use particle_mod      , only: fdm_particle_wkb_t

    integer,              intent(in   )        :: np
    type(fdm_particle_wkb_t), intent(inout)    :: particles(np)
    integer,              intent(in   )        :: accel_lo(3), accel_hi(3)
    real(amrex_real),     intent(in   )        :: accel &
         (accel_lo(1):accel_hi(1),accel_lo(2):accel_hi(2),accel_lo(3):accel_hi(3),3)
    real(amrex_real),     intent(in   )        :: plo(3),dx(3),dt,a_prev,a_cur
    integer,              intent(in   )        :: do_move

    integer          :: i, j, k, n, nc
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) half_dt, a_cur_inv, dt_a_cur_inv
    real(amrex_real) inv_dx(3)
    real(amrex_real) part_accel

    inv_dx = 1.0d0/dx

    half_dt       = 0.5d0 * dt
    a_cur_inv    = 1.d0 / a_cur;
    dt_a_cur_inv = dt * a_cur_inv;

    do n = 1, np

       if (particles(n)%id .le. 0) then
          cycle
       end if

       lx = (particles(n)%pos(1) - plo(1))*inv_dx(1) + 0.5d0
       ly = (particles(n)%pos(2) - plo(2))*inv_dx(2) + 0.5d0
       lz = (particles(n)%pos(3) - plo(3))*inv_dx(3) + 0.5d0

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       if (i-1 .lt. accel_lo(1) .or. i .gt. accel_hi(1) .or. &
           j-1 .lt. accel_lo(2) .or. j .gt. accel_hi(2) .or. &
           k-1 .lt. accel_lo(3) .or. k .gt. accel_hi(3)) then
          print *,'PARTICLE ID ', particles(n)%id,' REACHING OUT OF BOUNDS AT (I,J,K) (accel)= ',i,j,k
          call amrex_error('Aborting in move_kick_drift')
       end if

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       ! moveKickDrift: update velocity by half dt:  (a u)^half = (a u)^old  + dt/2 grav^old
       ! moveKick    t: update velocity by half dt:  (a u)^new  = (a u)^half + dt/2 grav^new
       do nc = 1, 3

          part_accel = &
                wx_lo*wy_lo*wz_lo*accel(i-1, j-1, k-1, nc) + &
                wx_lo*wy_lo*wz_hi*accel(i-1, j-1, k  , nc) + &
                wx_lo*wy_hi*wz_lo*accel(i-1, j,   k-1, nc) + &
                wx_lo*wy_hi*wz_hi*accel(i-1, j,   k  , nc) + &
                wx_hi*wy_lo*wz_lo*accel(i,   j-1, k-1, nc) + &
                wx_hi*wy_lo*wz_hi*accel(i,   j-1, k  , nc) + &
                wx_hi*wy_hi*wz_lo*accel(i,   j,   k-1, nc) + &
                wx_hi*wy_hi*wz_hi*accel(i,   j,   k  , nc)

          particles(n)%vel(nc) = a_prev*particles(n)%vel(nc) + half_dt * part_accel
          particles(n)%vel(nc) = particles(n)%vel(nc) * a_cur_inv

       end do

       ! moveKickDrift: Update position by full dt: x^new = x^old + dt u^half / a^half
       if (do_move .eq. 1) then
          do nc = 1, 3
             particles(n)%pos(nc) = particles(n)%pos(nc) + dt_a_cur_inv * particles(n)%vel(nc)
          end do
       end if

    end do

  end subroutine update_fdm_particles_wkb

  subroutine update_gaussian_beams_wkb(np, particles, accel, accel_lo, accel_hi, &
                                       phi, phi_lo, phi_hi, &
                                       plo, dx, dt, a_prev, a_cur, do_move) &
       bind(c,name='update_gaussian_beams_wkb')

    use iso_c_binding
    use amrex_error_module
    use amrex_fort_module, only : amrex_real
    use particle_mod      , only: fdm_particle_wkb_t
    use fdm_params_module, only : hbaroverm, ii!, sigma_fdm, gamma_fdm

    integer,              intent(in   )        :: np
    type(fdm_particle_wkb_t), intent(inout)    :: particles(np)
    integer,              intent(in   )        :: accel_lo(3), accel_hi(3)
    real(amrex_real),     intent(in   )        :: accel &
         (accel_lo(1):accel_hi(1),accel_lo(2):accel_hi(2),accel_lo(3):accel_hi(3),3)
    integer,              intent(in   )        :: phi_lo(3), phi_hi(3)
    real(amrex_real),     intent(in   )        :: phi &
         (phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3),1)
    real(amrex_real),     intent(in   )        :: plo(3),dx(3),dt,a_prev,a_cur
    integer,              intent(in   )        :: do_move

    integer          :: i, j, k, n, nc, nn, nn_hi
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) half_dt, a_cur_inv, dt_a_cur_inv
    real(amrex_real) inv_dx(3), Arg1, Arg2
    real(amrex_real) part_accel, part_phi, vel_square
    real(amrex_real) qq(3,3), pq(3,3), qqold(3,3), pqold(3,3), hessian(3,3)
    complex(amrex_real) Cpq, ZZ(3,3), det, Cpqtemp, Cpqold
    logical :: test

    inv_dx = 1.0d0/dx

    half_dt       = 0.5d0 * dt
    a_cur_inv    = 1.d0 / a_cur
    dt_a_cur_inv = dt * a_cur_inv

    do n = 1, np

       if (particles(n)%id .le. 0) then
          cycle
       end if

!                                                                                                                                                                                                              
!     Construct Gaussian beam parameters                                                                                                                                                                       
!                                                                                                                                                                                                              
       Cpq   = cmplx(particles(n)%amp(1),particles(n)%amp(2))
       qq    = reshape((/ particles(n)%qq(1) , particles(n)%qq(2) , particles(n)%qq(3) , &
                          particles(n)%qq(4) , particles(n)%qq(5) , particles(n)%qq(6) , &
                          particles(n)%qq(7) , particles(n)%qq(8) , particles(n)%qq(9) /), shape(qq))
       pq    = reshape((/ particles(n)%pq(1) , particles(n)%pq(2) , particles(n)%pq(3) , &
                          particles(n)%pq(4) , particles(n)%pq(5) , particles(n)%pq(6) , &
                          particles(n)%pq(7) , particles(n)%pq(8) , particles(n)%pq(9) /), shape(pq))
!
!      Store old data
!
       qqold = qq
       pqold = pq

       lx = (particles(n)%pos(1) - plo(1))*inv_dx(1) + 0.5d0
       ly = (particles(n)%pos(2) - plo(2))*inv_dx(2) + 0.5d0
       lz = (particles(n)%pos(3) - plo(3))*inv_dx(3) + 0.5d0

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       if (i-1 .lt. accel_lo(1) .or. i .gt. accel_hi(1) .or. &
           j-1 .lt. accel_lo(2) .or. j .gt. accel_hi(2) .or. &
           k-1 .lt. accel_lo(3) .or. k .gt. accel_hi(3)) then
          print *,'PARTICLE ID ', particles(n)%id,' REACHING OUT OF BOUNDS AT (I,J,K) (accel)= ',i,j,k
          call amrex_error('Aborting in move_kick_drift')
       end if

       if (i-1 .lt. phi_lo(1) .or. i .gt. phi_hi(1) .or. &
           j-1 .lt. phi_lo(2) .or. j .gt. phi_hi(2) .or. &
           k-1 .lt. phi_lo(3) .or. k .gt. phi_hi(3)) then
          print *,'PARTICLE ID ', particles(n)%id,' REACHING OUT OF BOUNDS AT (I,J,K) (phi)= ',i,j,k
          call amrex_error('Aborting in move_kick_drift')
       end if

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       ! moveKickDrift: update velocity by half dt:  (a u)^half = (a u)^old  + dt/2 grav^old
       ! moveKick    t: update velocity by half dt:  (a u)^new  = (a u)^half + dt/2 grav^new
       do nc = 1, 3

          part_accel = &
                wx_lo*wy_lo*wz_lo*accel(i-1, j-1, k-1, nc) + &
                wx_lo*wy_lo*wz_hi*accel(i-1, j-1, k  , nc) + &
                wx_lo*wy_hi*wz_lo*accel(i-1, j,   k-1, nc) + &
                wx_lo*wy_hi*wz_hi*accel(i-1, j,   k  , nc) + &
                wx_hi*wy_lo*wz_lo*accel(i,   j-1, k-1, nc) + &
                wx_hi*wy_lo*wz_hi*accel(i,   j-1, k  , nc) + &
                wx_hi*wy_hi*wz_lo*accel(i,   j,   k-1, nc) + &
                wx_hi*wy_hi*wz_hi*accel(i,   j,   k  , nc)

          particles(n)%vel(nc) = a_prev*particles(n)%vel(nc) + half_dt * part_accel
          particles(n)%vel(nc) = particles(n)%vel(nc) * a_cur_inv

       end do

       part_phi = &
            wx_lo*wy_lo*wz_lo*phi(i-1, j-1, k-1, 1) + &
                wx_lo*wy_lo*wz_hi*phi(i-1, j-1, k  , 1) + &
                wx_lo*wy_hi*wz_lo*phi(i-1, j,   k-1, 1) + &
                wx_lo*wy_hi*wz_hi*phi(i-1, j,   k  , 1) + &
                wx_hi*wy_lo*wz_lo*phi(i,   j-1, k-1, 1) + &
                wx_hi*wy_lo*wz_hi*phi(i,   j-1, k  , 1) + &
                wx_hi*wy_hi*wz_lo*phi(i,   j,   k-1, 1) + &
                wx_hi*wy_hi*wz_hi*phi(i,   j,   k  , 1)

       !Remember different sign convention for phi
       particles(n)%phase = particles(n)%phase + half_dt * part_phi

       hessian(1,1) = &
               (wy_hi*wz_hi*accel(i,   j,   k  , 1) + &
                wy_hi*wz_lo*accel(i,   j,   k-1, 1) + &
                wy_lo*wz_hi*accel(i,   j-1, k  , 1) + &
                wy_lo*wz_lo*accel(i,   j-1, k-1, 1) - &
                wy_hi*wz_hi*accel(i-1, j,   k  , 1) - &
                wy_hi*wz_lo*accel(i-1, j,   k-1, 1) - &
                wy_lo*wz_hi*accel(i-1, j-1, k  , 1) - &
                wy_lo*wz_lo*accel(i-1, j-1, k-1, 1)) * inv_dx(1)

       hessian(1,2) = &
               (wy_hi*wz_hi*accel(i,   j,   k  , 2) + &
                wy_hi*wz_lo*accel(i,   j,   k-1, 2) + &
                wy_lo*wz_hi*accel(i,   j-1, k  , 2) + &
                wy_lo*wz_lo*accel(i,   j-1, k-1, 2) - &
                wy_hi*wz_hi*accel(i-1, j,   k  , 2) - &
                wy_hi*wz_lo*accel(i-1, j,   k-1, 2) - &
                wy_lo*wz_hi*accel(i-1, j-1, k  , 2) - &
                wy_lo*wz_lo*accel(i-1, j-1, k-1, 2)) * inv_dx(1)

       hessian(1,3) = &
               (wy_hi*wz_hi*accel(i,   j,   k  , 3) + &
                wy_hi*wz_lo*accel(i,   j,   k-1, 3) + &
                wy_lo*wz_hi*accel(i,   j-1, k  , 3) + &
                wy_lo*wz_lo*accel(i,   j-1, k-1, 3) - &
                wy_hi*wz_hi*accel(i-1, j,   k  , 3) - &
                wy_hi*wz_lo*accel(i-1, j,   k-1, 3) - &
                wy_lo*wz_hi*accel(i-1, j-1, k  , 3) - &
                wy_lo*wz_lo*accel(i-1, j-1, k-1, 3)) * inv_dx(1)
       
       hessian(2,1) = hessian(1,2)
       
       hessian(2,2) = &
               (wx_hi*wz_hi*accel(i,   j,   k  , 2) + &
                wx_hi*wz_lo*accel(i,   j,   k-1, 2) + &
                wx_lo*wz_hi*accel(i-1, j,   k  , 2) + &
                wx_lo*wz_lo*accel(i-1, j,   k-1, 2) - &
                wx_hi*wz_hi*accel(i,   j-1, k  , 2) - &
                wx_hi*wz_lo*accel(i,   j-1, k-1, 2) - &
                wx_lo*wz_hi*accel(i-1, j-1, k  , 2) - &
                wx_lo*wz_lo*accel(i-1, j-1, k-1, 2)) * inv_dx(2)

       hessian(2,3) = &
               (wx_hi*wz_hi*accel(i,   j,   k  , 3) + &
                wx_hi*wz_lo*accel(i,   j,   k-1, 3) + &
                wx_lo*wz_hi*accel(i-1, j,   k  , 3) + &
                wx_lo*wz_lo*accel(i-1, j,   k-1, 3) - &
                wx_hi*wz_hi*accel(i,   j-1, k  , 3) - &
                wx_hi*wz_lo*accel(i,   j-1, k-1, 3) - &
                wx_lo*wz_hi*accel(i-1, j-1, k  , 3) - &
                wx_lo*wz_lo*accel(i-1, j-1, k-1, 3)) * inv_dx(2)

       hessian(3,1) = hessian(1,3)

       hessian(3,2) = hessian(2,3)

       hessian(3,3) = &
               (wx_hi*wy_hi*accel(i,   j,   k  , 3) + &
                wx_hi*wy_lo*accel(i,   j-1, k  , 3) + &
                wx_lo*wy_hi*accel(i-1, j,   k  , 3) + &
                wx_lo*wy_lo*accel(i-1, j-1, k  , 3) - &
                wx_hi*wy_hi*accel(i,   j  , k-1, 3) - &
                wx_hi*wy_lo*accel(i,   j-1, k-1, 3) - &
                wx_lo*wy_hi*accel(i-1, j  , k-1, 3) - &
                wx_lo*wy_lo*accel(i-1, j-1, k-1, 3)) * inv_dx(3)

       !Remember different sign convention for phi
       pq = pq + half_dt * matmul(hessian,qq)

       ! moveKickDrift: Update position by full dt: x^new = x^old + dt u^half / a^half
       if (do_move .eq. 1) then

          do nc = 1, 3
             particles(n)%pos(nc) = particles(n)%pos(nc) + dt_a_cur_inv * particles(n)%vel(nc)
          end do
          
          vel_square = particles(n)%vel(1)*particles(n)%vel(1) + &
                       particles(n)%vel(2)*particles(n)%vel(2) + &
                       particles(n)%vel(3)*particles(n)%vel(3) 

          particles(n)%phase = particles(n)%phase + half_dt * vel_square

          qq = qq + dt_a_cur_inv * a_cur_inv * pq
          
       else
!
!         Update complex beam amplitude Cpq
!          
          ZZ = particles(n)%width*qq-pq/(2.0*ii*hbaroverm)

          det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
               - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
               + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))

        Cpqtemp = sqrt(det)
        Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
        Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)

        if( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
           test = .false.
        else
           if ( Arg1 .le. Arg2 ) then
              Cpq = Cpqtemp
           else
              Cpq = -Cpqtemp
           endif
           test = .true.
        endif
!                                                                                                                                                                                                               
        if (test .eqv. .false.) then
           nn_hi = 10
           Cpqold = Cpq
           do nn=1, nn_hi
!                                                                                                                                                                                                               
              ZZ = particles(n)%width*((qqold*(nn_hi-nn)+nn*qq)/nn_hi)-((pqold*(nn_hi-nn)+nn*pq)/nn_hi)/(2.0*ii*hbaroverm)
!                                                                                                                                                                                                      
              det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
                   - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
                   + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))
!                                                                                                                                                                                                           
              Cpqtemp = sqrt(det)
              Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
              Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)
              if ( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
                 test = .false.
                 exit
              else
                 if ( Arg1 .le. Arg2 ) then
                    Cpq = Cpqtemp
                 else
                    Cpq = -Cpqtemp
                 endif
                 test = .true.
              endif
           enddo
        endif
!                                                                                                                                                                                                                
        if (test .eqv. .false.) then
           nn_hi = 100
           Cpq = Cpqold
           do nn=1, nn_hi
!                                                                                                                                                                                                                
              ZZ = particles(n)%width*((qqold*(nn_hi-nn)+nn*qq)/nn_hi)-((pqold*(nn_hi-nn)+nn*pq)/nn_hi)/(2.0*ii*hbaroverm)
!                                                                                                                                                                                                                
              det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
                   - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
                   + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))
!                                                                                                                                                                                                                
              Cpqtemp = sqrt(det)
              Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
              Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)
              if ( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
                 test = .false.
                 exit
              else
                 if ( Arg1 .le. Arg2 ) then
                    Cpq = Cpqtemp
                 else
                    Cpq = -Cpqtemp
                 endif
                 test = .true.
              endif
           enddo
        endif
!                                                                                                                                                                                                                
        if (test .eqv. .false.) then
           nn_hi = 1000
           Cpq = Cpqold
           do nn=1, nn_hi
!                                                                                                                                                                                                                
              ZZ = particles(n)%width*((qqold*(nn_hi-nn)+nn*qq)/nn_hi)-((pqold*(nn_hi-nn)+nn*pq)/nn_hi)/(2.0*ii*hbaroverm)
!                                                                                                                                                                                                              
              det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
                   - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
                   + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))
!                                                                                                                                                                                                                
              Cpqtemp = sqrt(det)
              Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
              Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)
              if ( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
                 test = .false.
                 exit
              else
                 if ( Arg1 .le. Arg2 ) then
                    Cpq = Cpqtemp
                 else
                    Cpq = -Cpqtemp
                 endif
                 test = .true.
              endif
           enddo
        endif
!                                                                                                                                                                                                                
        if (test .eqv. .false.) then
           nn_hi = 10000
           Cpq = Cpqold
           do nn=1, nn_hi
!                                                                                                                                                                                                                
              ZZ = particles(n)%width*((qqold*(nn_hi-nn)+nn*qq)/nn_hi)-((pqold*(nn_hi-nn)+nn*pq)/nn_hi)/(2.0*ii*hbaroverm)
!                                                                                                                                                                                                              
              det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
                   - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
                   + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))
!                                                                                                                                                                                                               
              Cpqtemp = sqrt(det)
              Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
              Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)
              if ( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
                 test = .false.
                 exit
              else
                 if ( Arg1 .le. Arg2 ) then
                    Cpq = Cpqtemp
                 else
                    Cpq = -Cpqtemp
                 endif
                 test = .true.
              endif
           enddo
        endif
!                                                                                                                                                                                                                
        if (test .eqv. .false.) then
           nn_hi = 100000
           Cpq = Cpqold
           do nn=1, nn_hi
!                                                                                                                                                                                                                
              ZZ = particles(n)%width*((qqold*(nn_hi-nn)+nn*qq)/nn_hi)-((pqold*(nn_hi-nn)+nn*pq)/nn_hi)/(2.0*ii*hbaroverm)
!                                                                                                                                                                                                                
              det = (ZZ(1,1)*ZZ(2,2)*ZZ(3,3) - ZZ(1,1)*ZZ(2,3)*ZZ(3,2) &
                   - ZZ(1,2)*ZZ(2,1)*ZZ(3,3) + ZZ(1,2)*ZZ(2,3)*ZZ(3,1) &
                   + ZZ(1,3)*ZZ(2,1)*ZZ(3,2) - ZZ(1,3)*ZZ(2,2)*ZZ(3,1))
!                                                                                                                                                                                                                
              Cpqtemp = sqrt(det)
              Arg1 = abs(Cpq-Cpqtemp)/abs(Cpq)
              Arg2 = abs(Cpq+Cpqtemp)/abs(Cpq)
              if ( (Arg1 .gt. 0.1) .and. (Arg2 .gt. 0.1) ) then
                 print *, "Error: functional determinant is changing too fast! Arguments are: ",Arg1, Arg2, Cpqtemp, Cpq, nn, nn_hi
                 call amrex_error('Aborting in update_gaussian_beams')
              else
                 if ( Arg1 .le. Arg2 ) then
                    Cpq = Cpqtemp
                 else
                    Cpq = -Cpqtemp
                 endif
              endif
           enddo
        endif

        particles(n)%amp(1) = real(real(Cpq))
        particles(n)%amp(2) = real(aimag(Cpq))

       end if

       particles(n)%qq(1) = qq(1,1)
       particles(n)%qq(2) = qq(2,1)
       particles(n)%qq(3) = qq(3,1)
       particles(n)%qq(4) = qq(1,2)
       particles(n)%qq(5) = qq(2,2)
       particles(n)%qq(6) = qq(3,2)
       particles(n)%qq(7) = qq(1,3)
       particles(n)%qq(8) = qq(2,3)
       particles(n)%qq(9) = qq(3,3)

       particles(n)%pq(1) = pq(1,1)
       particles(n)%pq(2) = pq(2,1)
       particles(n)%pq(3) = pq(3,1)
       particles(n)%pq(4) = pq(1,2)
       particles(n)%pq(5) = pq(2,2)
       particles(n)%pq(6) = pq(3,2)
       particles(n)%pq(7) = pq(1,3)
       particles(n)%pq(8) = pq(2,3)
       particles(n)%pq(9) = pq(3,3)

    end do

  end subroutine update_gaussian_beams_wkb
