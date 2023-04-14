! FILE: utils.f90

subroutine func(NX, u, x, fun)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer,                 intent(in)  :: NX
    real(RK), dimension(NX), intent(in)  :: u, x
    real(RK), dimension(NX), intent(out) :: fun

    integer :: i

    do i = 1, NX
        fun(i) = 0.5_RK * x(i)**3 * sin(u(i) / x(i)**2)**2
    end do

return
end subroutine func

subroutine fprime(NX, u, x, f_prime)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer,                 intent(in)  :: NX
    real(RK), dimension(NX), intent(in)  :: u, x
    real(RK), dimension(NX), intent(out) :: f_prime

    integer :: i

    do i = 1, NX
        f_prime(i) = 0.5_RK * x(i) * sin(2_RK * u(i) / x(i)**2)
    end do

return
end subroutine fprime


function heaviside(x)

    use iso_fortran_env, only: RK => real64
    implicit none

    real(RK), intent(in) :: x
    real(RK)             :: heaviside

    if ( x.gt.0 ) then
        heaviside = 1_RK
    else if (x.lt.0) then
        heaviside = 0_RK
    else
        heaviside = 0.5_RK
    end if

    return
end function heaviside

subroutine flux(a, b, x_surf, out)

    use iso_fortran_env, only: RK => real64
    implicit none

    real(RK), intent(in) :: a, b, x_surf
    real(RK), intent(out):: out

    real(RK) :: ul, ur, FL, FR
    real(RK), parameter :: PI=4._RK*DATAN(1._RK)

    ul = a / x_surf**2_RK
    ur = b / x_surf**2_RK
    FL = 0.5_RK * x_surf**3 * sin(ul)**2 
    FR = 0.5_RK * x_surf**3 * sin(ur)**2

    if (ul.lt.ur) then
        out = min(FL, FR)
    else if( ul.gt.ur) then
        if ((ur.gt.-PI/2_RK).or.(ul.lt.-PI/2_RK)) then
            out = max(FL, FR)
        else
            out = 0.5_RK * x_surf**3
        end if
    end if
return
end subroutine flux


subroutine BC(NX, arr, nghost)
    use iso_fortran_env, only: RK => real64
    implicit none

    integer,                 intent(in)    :: NX, nghost
    real(RK), dimension(NX), intent(inout) :: arr
    
    integer :: i
    
    do i = 1, nghost
        arr(i) = arr(nghost + 1 - i)
    end do

return
end subroutine BC