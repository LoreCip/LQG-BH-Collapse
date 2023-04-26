subroutine initial_data(NX, x, dx, m, r0, a0, idx, u)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer                , intent(in) :: NX, idx
    real(RK)               , intent(in) :: dx, m, r0, a0
    real(RK), dimension(NX), intent(in) :: x
    real(RK), dimension(2*NX), intent(out):: u

    integer :: i, j
    real(RK) :: th, Mass, heaviside, p
    real(RK), dimension(NX) :: rr, Marray

    real(RK), parameter :: PI=4._RK*DATAN(1._RK)

    if ( idx .eq. 0 ) then

        ! Physical values
        do i = 1, NX
            th = heaviside(r0 - x(i))
            Mass = m*x(i)**3_RK / r0**3_RK * th + m * (1_RK - th)
            ! E(x)
            u(NX+i) = + x(i)**2 / a0**2 * th + r0**2 / a0**2 * (1_RK - th)
            ! B(x)
            u(i) = - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Mass / x(i)**3_RK - 2_RK * u(NX+i) / x(i)**2_RK)
        end do
        
    else if ( idx .eq. 1 ) then

        do i = 1, NX
            rr(i) = 3_RK * m * (PI/2_RK - atan(x(i) - r0) )  / (8_RK * PI * r0**3)
        end do

        do i = 3, NX-1
            p = 0_RK
            do j = 3, i
                p = p + (rr(j-1) * x(j-1)**2 + rr(j) * x(j)**2)
            end do
            Marray(i) = 4_RK * PI * p * dx * 0.5_RK
        end do
        Marray(2) = 0

        ! Physical values
        do i = 2, NX-1
            th = heaviside(r0 - x(i))
            ! E(x)
            u(NX+i) = x(i)**2 / a0**2 * th + r0**2 / a0**2 * (1_RK - th)
            ! B(x)
            u(i) = - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Marray(i) / x(i)**3_RK - 2_RK * u(NX+i) /  x(i)**2_RK)
        end do

    end if
    
    u(2) = 0_RK
    u(NX+2) = 0_RK

    ! Ghosts
    u(1) = u(2)
    u(NX) = u(NX-1)
    u(NX+1) = u(NX+2)
    u(2*NX) = u(2*NX-1)

return
end subroutine initial_data