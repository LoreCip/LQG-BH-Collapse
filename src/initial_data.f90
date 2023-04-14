    subroutine initial_data(NX, x, m, r0, u)

        use iso_fortran_env, only: RK => real64
        implicit none

        integer                , intent(in) :: NX
        real(RK)               , intent(in) :: m, r0
        real(RK), dimension(NX), intent(in) :: x
        real(RK), dimension(NX), intent(out):: u

        integer :: i
        real(RK) :: th, Mass, heaviside

        !$OMP PARALLEL DO
        do i = 1, NX
            th = heaviside(r0 - x(i))
            Mass = m*x(i)**3_RK / r0**3_RK * th + m * (1_RK - th)
            u(i) = - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Mass / x(i)**3_RK)
        end do
        !$OMP END PARALLEL DO

    return
    end subroutine initial_data
