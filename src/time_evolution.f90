subroutine TVD_RK(m, NX, u_p, u1, u2, u12, u32, f_prime, x, dx, dt, nghost, uf, e_der, rho, betamin)

    use iso_fortran_env, only: RK => real64
    use OMP_LIB
    implicit none

    integer                , intent(in) :: NX, nghost
    real(RK)               , intent(in) :: dx, dt, m
    real(RK)               , intent(inout)   :: betamin
    real(RK), dimension(NX), intent(in)      :: x
    real(RK), dimension(NX), intent(inout)   :: f_prime, e_der, rho
    real(RK), dimension(2*NX), intent(in)    :: u_p
    real(RK), dimension(2*NX), intent(inout) :: u1, u2, u12, u32
    real(RK), dimension(2*NX), intent(inout) :: uf

    integer :: i, idx
    real(RK), parameter :: PI=4._RK*DATAN(1._RK)

    logical, parameter :: correction = .true.

    ! 1**) t      --> t +   dt
    call RK_STEP(NX, u_p, f_prime, x, dx, dt, nghost, u1)
    
    ! 2**) t + dt --> t + 2*dt
    call RK_STEP(NX, u1, f_prime, x, dx, dt, nghost, u2)

    ! 3**) t + dt/2
!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i)
    do i = 1, 2*NX
        u12(i) = 3_RK*u_p(i)/4_RK + u2(i)/4_RK
    end do
!$OMP END DO SIMD

!$OMP SINGLE
    call BC(NX, u12)
!$OMP END SINGLE

    ! 4**) t + dt/2 -- > t + 3*dt/2
    call RK_STEP(NX, u12, f_prime, x, dx, dt, nghost, u32)

!$OMP SINGLE
    betamin = 1.0E20
!$OMP END SINGLE
    ! 5**) t + dt
!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i) REDUCTION(MIN:betamin)
    do i = 1, 2*NX
        uf(i) = u_p(i)/3_RK + 2_RK*u32(i)/3_RK
        if ((i.le.NX).and.(x(i).gt.(2_RK*m)**(1._RK/3._RK))) then
            betamin = min(betamin, uf(i) / x(i)**2)
        end if    
    end do
!$OMP END DO SIMD

!$OMP SINGLE
    call BC(NX, uf)
!$OMP END SINGLE

    if (correction) then
        if (betamin.lt.-PI/2_RK) then
            call compRho(NX, dx, dt, uf(1:NX), u_p(1:NX), uf(NX+1:2*NX), x, e_der, rho)
!$OMP SINGLE
            idx = maxloc(rho, DIM=1, mask=(x>(2_RK*m)**(1._RK/3._RK))) - 1
            if (x(idx+1).gt.(2_RK*m)**(1._RK/3._RK)) then
                ! Try moving one point
                uf(NX+idx) = uf(NX+idx-1)
            end if
!$OMP END SINGLE
        end if
    end if

    return
end subroutine TVD_RK

subroutine RK_STEP(NX, u, f_prime, x, dx, dt, nghost, u_step)

    use iso_fortran_env, only: RK => real64
    use OMP_LIB
    implicit none

    integer                , intent(in) :: NX, nghost
    real(RK)               , intent(in) :: dx, dt
    real(RK), dimension(NX), intent(in) :: x
    real(RK), dimension(NX), intent(inout) :: f_prime
    real(RK), dimension(2*NX), intent(in) :: u
    real(RK), dimension(2*NX), intent(out):: u_step
    
    integer :: i
    real(RK) :: L, Bp, Bm, e_k, e_l, vb
    real(RK), parameter :: PI=4._RK*DATAN(1._RK)

!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i)
    do i = 1, NX
        f_prime(i) = vb(u(i), x(i))
    end do
!$OMP END DO SIMD NOWAIT

!$OMP SINGLE
    u_step(2) = u(2)           ! First phys not evolve
    u_step(NX-1) = u(NX-1)     ! Last phys not evolve

    u_step(NX+2) = u(NX+2)         ! First phys not evolve
    u_step(2*NX-1) = u(2*NX-1)     ! Last phys not evolve

    ! Ghosts
    call BC(NX, u_step)
!$OMP END SINGLE

!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i, L, Bp, Bm, e_k, e_l)
    do i = 3, NX-2
        call WENO_RHS(NX, i, u, f_prime, x, dx, nghost, L, Bp, Bm, e_k, e_l)
        u_step(i)    = u(i) + dt * (L + 0.5_RK * u(NX+i))
        u_step(NX+i) = u(NX+i) + dt * ( vb(Bp, x(i)+0.5_RK*dx) * (u(NX+i) - e_k) + vb(Bm, x(i)-0.5_RK*dx) * (e_l - u(NX+i)) ) / dx
    end do
!$OMP END DO SIMD

    return
end subroutine RK_STEP

subroutine CLF(NX, u, x, dx, vel, dt)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer                , intent(in) :: NX
    real(RK)               , intent(in) :: dx
    real(RK), dimension(NX), intent(in) :: u, x
    real(RK), dimension(NX), intent(inout) :: vel
    real(RK)               , intent(out):: dt

    integer             :: i
    real(RK)            :: v_min, v_max, v_abs
    real(RK), parameter :: fact = 0.25_RK

!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i)
    do i = 1, NX
        vel(i) = 0.5_RK * x(i) * sin(2_RK * u(i) / x(i)**2_RK)
    end do
!$OMP END DO SIMD
!$OMP SINGLE PRIVATE(v_min, v_max, v_abs)
    v_min = minval(vel)
    v_max = maxval(vel)
    v_abs = abs(max(-v_min, v_max))

    dt = fact * dx / v_abs
    if (dt .gt. 0.01_RK*dx) then
       dt = 0.01_RK*dx  ! largest timestep allowed
    end if
!$OMP END SINGLE
    return
end subroutine CLF

subroutine BC(NX, arr)
    
    use iso_fortran_env, only: RK => real64
    implicit none

    integer,                 intent(in)    :: NX
    real(RK), dimension(2*NX), intent(inout) :: arr

    arr(1) = arr(2)
    arr(NX) = arr(NX-1)

    arr(NX+1) = arr(NX+2)
    arr(2*NX) = arr(2*NX-1) 

    return
end subroutine BC
