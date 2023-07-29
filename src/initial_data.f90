subroutine initial_data(NX, x, dx, m, r0, a0, idx, u)

    use iso_fortran_env, only: RK => real64
    use OMP_LIB
    implicit none

    integer                , intent(in) :: NX, idx
    real(RK)               , intent(in) :: dx, m, r0, a0
    real(RK), dimension(NX), intent(in) :: x
    real(RK), dimension(2*NX), intent(out):: u

    integer :: i, j
    real(RK) :: th, Mass, heaviside, p
    real(RK), dimension(NX) :: rr, Marray

    real(RK), parameter :: PI=4._RK*DATAN(1._RK)

    ! Full dynamics, starting from step function

!$OMP PARALLEL
    if ( idx .eq. 0 ) then

        ! Physical values
!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i)
        do i = 1, NX
            th = heaviside(r0 - x(i))
            Mass = m*x(i)**3_RK / r0**3_RK * th + m * (1_RK - th)
            ! E(x)
            u(NX+i) = - x(i)**2 / a0**2 * th - r0**2 / a0**2 * (1_RK - th)
            ! B(x)
            u(i) = - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Mass / x(i)**3_RK - 2_RK * u(NX+i) / x(i)**2_RK)
        end do
!$OMP END DO SIMD

    ! Post bounce dynamics, starting from step function
    else if ( idx .eq. 1 ) then

        ! Physical values
!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i)
        do i = 1, NX-1

            th = heaviside(r0 - x(i))
            Mass = m*x(i)**3_RK / r0**3_RK * th + m * (1_RK - th)
            ! E(x)
            u(NX+i) = - x(i)**2 / a0**2 * th - r0**2 / a0**2 * (1_RK - th)
            ! B(x)
            u(i) =  - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Mass / x(i)**3_RK - 2_RK * u(NX+i) / x(i)**2_RK) * (1_RK - th)  &
                    + x(i)**2 * ( - PI + asin(sqrt(2_RK * Mass / x(i)**3_RK + u(NX+i) / x(i)**2_RK)) ) * th
        end do
!$OMP END DO SIMD

    ! Post bounce dynamics, starting from delta function
    else if ( idx .eq. 2 ) then

        ! Physical values
!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i)
        do i = 1, NX-1
            th = heaviside(r0 - x(i))
            ! E(x)
            u(NX+i) = - r0**2 / a0**2 * (1_RK - th)
            ! B(x)
            u(i) =  - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * m / x(i)**3_RK - 2_RK * u(NX+i) / x(i)**2_RK) * (1_RK - th)  &
                    - PI * x(i)**2 * th
        end do
!$OMP END DO SIMD

    ! Full dynamics, starting from atan density profile
    else if ( idx .eq. 3 ) then

!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i)
        do i = 1, NX
            rr(i) = 3_RK * m * (PI/2_RK - atan(x(i) - r0) )  / (8_RK * PI * r0**3)
        end do
!$OMP END DO SIMD
        
!$OMP DO PRIVATE(i, j, p)
        do i = 3, NX-1
            p = 0_RK
!$OMP SIMD REDUCTION(+: p)
            do j = 3, i
                p = p + (rr(j-1) * x(j-1)**2 + rr(j) * x(j)**2)
            end do
!$OMP END SIMD
            Marray(i) = 4_RK * PI * p * dx * 0.5_RK
        end do
!$OMP END DO

!$OMP MASTER
        Marray(2) = 0
!$OMP END MASTER

        ! Physical values
!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i)
        do i = 2, NX-1
            th = heaviside(r0 - x(i))
            if ((x(i).le.(r0+1_RK)).and.(x(i).ge.(r0-1_RK))) then
                u(NX+i) = x(i)**3 - (r0 + 2_RK) * x(i)**2 - (r0**2 - 1_RK) * x(i) + r0*(r0 - 1_RK)**2
                u(NX+i) = u(NX+i) / (4_RK*a0**2)
            else
                u(NX+i) = - x(i)**2 / a0**2 * th - r0**2 / a0**2 * (1-th)
            end if
            ! B(x)
            u(i) = - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Marray(i) / x(i)**3_RK - 2_RK * u(NX+i) /  x(i)**2_RK)
        end do
!$OMP END DO SIMD

    ! Full dynamics, flat case
    else if ( idx .eq. 4 ) then

        ! Physical values
!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i)
        do i = 1, NX
            th = heaviside(r0 - x(i))
            Mass = m*x(i)**3_RK / r0**3_RK * th + m * (1_RK - th)
            ! E(x)
            u(NX+i) = 0_RK
            ! B(x)
            u(i) = - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Mass / x(i)**3_RK - 2_RK * u(NX+i) / x(i)**2_RK)
        end do
!$OMP END DO SIMD

    ! Full dynamics, open case
    else if ( idx .eq. 5 ) then

        ! Physical values
!$OMP DO SIMD SCHEDULE(STATIC) PRIVATE(i)
        do i = 1, NX
            th = heaviside(r0 - x(i))
            Mass = m*x(i)**3_RK / r0**3_RK * th + m * (1_RK - th)
            ! E(x)
            u(NX+i) = + x(i)**2 / a0**2 * th + r0**2 / a0**2 * (1_RK - th)
            ! B(x)
            u(i) = - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Mass / x(i)**3_RK - 2_RK * u(NX+i) / x(i)**2_RK)
        end do
!$OMP END DO SIMD
    end if
!$OMP END PARALLEL

    u(2) = 0_RK
    u(NX+2) = 0_RK

    ! Ghosts
    u(1) = u(2)
    u(NX) = u(NX-1)
    u(NX+1) = u(NX+2)
    u(2*NX) = u(2*NX-1)
    
return
end subroutine initial_data