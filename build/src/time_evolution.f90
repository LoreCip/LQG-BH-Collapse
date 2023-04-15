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

        ! 1**) t      --> t +   dt
        call RK_STEP(NX, u_p, x, dx, dt, nghost, u1)
        ! 2**) t + dt --> t + 2*dt
        call RK_STEP(NX, u1, x, dx, dt, nghost, u2)

        ! 3**) t + dt/2
        do i = 1, NX
            u12(i) = 3_RK*u_p(i)/4_RK + u2(i)/4_RK
        end do
        call BC(NX, u12, nghost)

        ! 4**) t + dt/2 -- > t + 3*dt/2
        call RK_STEP(NX, u12, x, dx, dt, nghost, u32)

        ! 5**) t + dt
        do i = 1, NX
            uf(i) = u_p(i)/3_RK + 2_RK*u32(i)/3_RK
        end do
        call BC(NX, uf, nghost)

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
        real(RK) :: L
        real(RK), dimension(NX) :: f_prime


        call fprime(NX, u, x, f_prime)
        !$OMP PARALLEL DO
        do i = 1, NX
            ! Skip ghost cells
            if ((i.le.nghost).or.(i.gt.NX-nghost)) cycle

            ! If it is the first physical, set it to zero
            if (i.eq.(1 + nghost)) then
                u_step(i) = 0
            ! If it is the last physical, set it to the previous values
            else if (i.eq.(NX - nghost - 2)) then
                u_step(i) = u(i)
                u_step(i+1) = u(i+1)
            else
                call WENO_RHS(NX, i, u, f_prime, x, dx, nghost, L)
                u_step(i) = u(i) + dt * L
            end if
        end do
        !$OMP END PARALLEL DO
        ! Take care of the boundary conditions on the ghosts
        call BC(NX, u_step, nghost)

        return
    end subroutine RK_STEP

    subroutine CLF(NX, u, x, dx, dt)

        use iso_fortran_env, only: RK => real64
        use OMP_LIB
        implicit none

        integer                , intent(in) :: NX
        real(RK)               , intent(in) :: dx
        real(RK), dimension(NX), intent(in) :: u, x
        real(RK)               , intent(out):: dt

        integer :: i
        real(RK)                :: v_min, v_max, v_abs
        real(RK), dimension(NX) :: vel

        real(RK), parameter :: fact = 0.25_RK

        !$OMP PARALLEL DO
        do i = 1, NX
            vel(i) = 0.5_RK * x(i) * sin(2_RK * u(i) / x(i)**2_RK)
        end do
        !$OMP END PARALLEL DO

        v_min = minval(vel)
        v_max = maxval(vel)
        v_abs = abs(max(-v_min, v_max))

        dt = fact * dx / v_abs
        if (dt .gt. 0.01_RK*dx) then
           dt = 0.01_RK*dx  ! largest timestep allowed
        end if

        return
    end subroutine CLF

    subroutine BC(NX, arr, nghost)
        
        use iso_fortran_env, only: RK => real64
        use OMP_LIB
        implicit none

        integer,                 intent(in)    :: NX, nghost
        real(RK), dimension(NX), intent(inout) :: arr

        integer :: i

        do i = 1, nghost
            arr(i) = arr(nghost + 1 - i)
        end do

        return
    end subroutine BC
