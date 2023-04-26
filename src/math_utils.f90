function vb(u, x) result(f_prime)

    use iso_fortran_env, only: RK => real64
    implicit none

    real(RK), intent(in)  :: u, x
    real(RK) :: f_prime

    f_prime = 0.5_RK * x * sin(2_RK * u / x**2)

    return
end function vb


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

subroutine compRho(NX, dx, dt, B, BP, E, x, out)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer,                 intent(in)  :: NX
    real(RK),                intent(in)  :: dx, dt
    real(RK), dimension(NX), intent(in)  :: B, BP, E, x
    real(RK), dimension(NX), intent(out) :: out

    integer :: i
    real(RK), parameter :: PI=4._RK*DATAN(1._RK)
    real(RK), dimension(NX) :: e_der

    do i = 3, NX-2
        e_der(i) = ( E(i-2) - 8*E(i-1) + 8*E(i+1) - E(i+2) ) / (12_RK * dx) 
    end do
    e_der(2) = ( -25_RK*E(i) + 48_RK*E(i+1) - 36_RK*E(i+2) + 16_RK*E(i+3) - 3_RK*E(i+4) ) / (12_RK * dx) 
    e_der(NX-1) = ( 25_RK*E(i) - 48_RK*E(i-1) + 36_RK*E(i-2) - 16_RK*E(i-3) + 3_RK*E(i-4) ) / (12_RK * dx)
    e_der(1) = e_der(2)
    e_der(NX) = e_der(NX-1)

    do i = 2, NX-1
        out(i) = - ( (B(i) - BP(i))/dt + x(i) * e_der(i) / 2_RK ) / (4_RK*PI*x(i)**2)
    end do
    out(1)  = out(2)
    out(NX) = 0

    return
end subroutine compRho