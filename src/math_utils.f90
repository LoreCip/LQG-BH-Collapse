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

subroutine compRho(NX, dx, dt, B, BP, E, x, e_der, out)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer,                 intent(in)  :: NX
    real(RK),                intent(in)  :: dx, dt
    real(RK), dimension(NX), intent(in)  :: B, BP, E, x
    real(RK), dimension(NX), intent(inout) :: e_der
    real(RK), dimension(NX), intent(out) :: out

    integer :: i
    
    real(RK), parameter :: PI=4._RK*DATAN(1._RK)

!$OMP DO SCHEDULE(STATIC) PRIVATE(i)
    do i = 3, NX-2
        e_der(i) = ( E(i-2) - 8*E(i-1) + 8*E(i+1) - E(i+2) ) / (12_RK * dx) 
    end do
!$OMP END DO NOWAIT

!$OMP SINGLE
    e_der(2) = ( -25_RK*E(2) + 48_RK*E(3) - 36_RK*E(4) + 16_RK*E(5) - 3_RK*E(6) ) / (12_RK * dx) 
    e_der(NX-1) = ( 25_RK*E(NX-1) - 48_RK*E(NX-2) + 36_RK*E(NX-3) - 16_RK*E(NX-4) + 3_RK*E(NX-5) ) / (12_RK * dx)
    e_der(1) = e_der(2)
    e_der(NX) = e_der(NX-1)
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC) PRIVATE(i)
    do i = 2, NX-1
        out(i) = - ( (B(i) - BP(i))/dt + x(i) * e_der(i) / 2_RK ) / (4_RK*PI*x(i)**2)
    end do
!$OMP END DO
!$OMP SINGLE
    out(1)  = out(2)
    out(NX) = 0
!$OMP END SINGLE
    
    return
end subroutine compRho

subroutine CompExpansion(NX, B, E, x, theta)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer,                 intent(in)  :: NX
    real(RK), dimension(NX), intent(in)  :: B, E, x
    real(RK), dimension(NX), intent(inout) :: theta

    integer :: i
    
!$OMP DO SCHEDULE(STATIC) PRIVATE(i)
    do i = 1, NX
        theta(i) = 1_RK - x(i)**2 * sin(2_RK * B(i) / x(i)**2)**2 / ( 4_RK * (1_RK + E(i)) )
    end do
!$OMP END DO

    return
end subroutine CompExpansion