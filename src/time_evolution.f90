subroutine TVD_RK(NX, u_p, x, dx, dt, nghost, uf)

    use iso_fortran_env, only: RK => real64
    use OMP_LIB
    implicit none

    integer                , intent(in) :: NX, nghost
    real(RK)               , intent(in) :: dx, dt
    real(RK), dimension(NX), intent(in) :: u_p, x
    real(RK), dimension(NX), intent(out):: uf

    real(RK), dimension(NX) :: u1, u2, u12, u32
    integer :: i

!$OMP PARALLEL
    ! 1**) t      --> t +   dt
    call RK_STEP(NX, u_p, x, dx, dt, nghost, u1)

    ! 2**) t + dt --> t + 2*dt
    call RK_STEP(NX, u1, x, dx, dt, nghost, u2)

    ! 3**) t + dt/2
!$OMP DO SCHEDULE(STATIC) private(i)
    do i = 1, 2*NX
        u12(i) = 3_RK*u_p(i)/4_RK + u2(i)/4_RK
    end do
!$OMP END DO

!$OMP SINGLE
    call BC(NX, u12)
!$OMP END SINGLE

    ! 4**) t + dt/2 -- > t + 3*dt/2
    call RK_STEP(NX, u12, x, dx, dt, nghost, u32)

    ! 5**) t + dt
!$OMP DO SCHEDULE(STATIC) private(i)
    do i = 1, 2*NX
        uf(i) = u_p(i)/3_RK + 2_RK*u32(i)/3_RK
    end do
!$OMP END DO
!$OMP END PARALLEL 

    call BC(NX, uf)

    return
end subroutine TVD_RK

subroutine RK_STEP(NX, u, x, dx, dt, nghost, u_step)

    use iso_fortran_env, only: RK => real64
    use OMP_LIB
    implicit none

    integer                , intent(in) :: NX, nghost
    real(RK)               , intent(in) :: dx, dt
    real(RK), dimension(NX), intent(in) :: u, x
    real(RK), dimension(NX), intent(out):: u_step

    integer :: i
    real(RK) :: L, Bp, Bm, e_k, e_l, vb
    real(RK), dimension(NX) :: f_prime

    call fprime(NX, u, x, f_prime)

!$OMP SINGLE
    u_step(2) = u(2)           ! First phys not evolve
    u_step(1) = u_step(2)      ! Ghost
    u_step(NX-1) = u(NX-1)     ! Last phys not evolve
    u_step(NX) = u_step(NX-1)  ! Ghost

    u_step(NX+2) = u(NX+2)           ! First phys not evolve
    u_step(NX+1) = u_step(NX+2)      ! Ghost
    u_step(2*NX-1) = u(2*NX-1)     ! Last phys not evolve
    u_step(2*NX) = u_step(2*NX-1)  ! Ghost
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC) PRIVATE(i, L, Bp, Bm, e_k, e_l)
    do i = 3, NX-2
        call WENO_RHS(NX, i, u, f_prime, x, dx, nghost, L, Bp, Bm, e_k, e_l)
        u_step(i)    = u(i) + dt * (L + u(NX + i))
        u_step(NX+i) = u(NX+i) + dt * ( vb(Bp, x(i)+0.5_RK*dx) * (u(NX+i) - e_k) + vb(Bm, x(i)-0.5_RK*dx) * (e_l - u(NX+i)))
    end do
!$OMP END DO
!$OMP BARRIER

    ! Take care of the boundary conditions on the ghosts
    ! call BC(NX, u_step, nghost)

    return
end subroutine RK_STEP

subroutine CLF(NX, u, x, dx, dt)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer                , intent(in) :: NX
    real(RK)               , intent(in) :: dx
    real(RK), dimension(NX), intent(in) :: u, x
    real(RK)               , intent(out):: dt

    integer :: i
    real(RK)                :: v_min, v_max, v_abs
    real(RK), dimension(NX) :: vel
    real(RK), parameter :: fact = 0.25_RK

    do i = 1, NX
        vel(i) = 0.5_RK * x(i) * sin(2_RK * u(i) / x(i)**2_RK)
    end do

    v_min = minval(vel)
    v_max = maxval(vel)
    v_abs = abs(max(-v_min, v_max))

    dt = fact * dx / v_abs
    if (dt .gt. 0.01_RK*dx) then
       dt = 0.01_RK*dx  ! largest timestep allowed
    end if

    return
end subroutine CLF

subroutine BC(NX, arr)
    
    use iso_fortran_env, only: RK => real64
    implicit none

    integer,                 intent(in)    :: NX
    real(RK), dimension(NX), intent(inout) :: arr

    arr(1) = arr(2)
    arr(NX) = arr(NX-1)

    arr(NX+1) = arr(NX+2)
    arr(2*NX) = arr(2*NX-1) 

    return
end subroutine BC
