subroutine fprime(NX, u, x, f_prime)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer,                 intent(in)  :: NX
    real(RK), dimension(NX), intent(in)  :: u, x
    real(RK), dimension(NX), intent(out) :: f_prime

    integer :: i
!$OMP DO SCHEDULE(STATIC) 
    do i = 1, NX
        f_prime(i) = 0.5_RK * x(i) * sin(2_RK * u(i) / x(i)**2)
    end do
!$OMP END DO
    return
end subroutine fprime


function heaviside(x) result(out)

    use iso_fortran_env, only: RK => real64
    implicit none

    real(RK), intent(in) :: x
    real(RK)             :: out

    if ( x.gt.0 ) then
        out = 1_RK
    else if (x.lt.0) then
        out = 0_RK
    else
        out = 0.5_RK
    end if

    return
end function heaviside

function interpolant(NX, i, x_pt, u, x, dx) result(retval)
    
    use iso_fortran_env, only: RK => real64
    implicit none
    
    integer                , intent(in) :: i, NX
    real(RK)               , intent(in) :: dx, x_pt
    real(RK), dimension(NX), intent(in) :: u, x
    real(RK) :: retval
    
    retval = u(i-1) + (u(i) - u(i-1)) * (x_pt - x(i-1)) / dx
    
    return
end function interpolant