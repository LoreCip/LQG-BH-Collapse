! FILE: WENO_reconstruction.f90

subroutine WENO_RHS(NX, i, u, f_prime, x, dx, nghost, L)
    
    use iso_fortran_env, only: RK => real64
    use OMP_LIB
    implicit none
    
    integer                , intent(in) :: i, NX, nghost
    real(RK)               , intent(in) :: dx
    real(RK), dimension(NX), intent(in) :: u, x, f_prime
    real(RK),                intent(out):: L
    
    real(RK) :: x_plus, x_minus, R1, R2, R3, R4, o1, o2
    
    x_plus  = x(i) + dx / 2_RK
    x_minus = x(i) - dx / 2_RK
    
    call R(NX, i  , x_plus,  u, f_prime, x, dx, R1)
    call R(NX, i+1, x_plus,  u, f_prime, x, dx, R2)
    call R(NX, i-1, x_minus, u, f_prime, x, dx, R3)
    call R(NX, i  , x_minus, u, f_prime, x, dx, R4)
    
    if (i.eq.nghost) then
        call flux(R1, R2, x_plus, o1) 
        o2 = 0_RK
    else
        call flux( R1, R2, x_plus , o1) 
        call flux( R3, R4, x_minus, o2)
    end if
    L = - ( o1 - o2 ) / dx
    
    return 
end subroutine WENO_RHS


subroutine R(NX, i, x_pt, u, f_prime, x, dx, Rout)
    
    use iso_fortran_env, only: RK => real64
    use OMP_LIB
    implicit none
    
    integer                , intent(in) :: i, NX
    real(RK)               , intent(in) :: dx, x_pt
    real(RK), dimension(NX), intent(in) :: u, x, f_prime
    real(RK),                intent(out):: Rout
    
    real(RK) :: aj0, aj1, w0, w1, interpolant
    
    call alpha(NX, i, 0, u, f_prime, aj0)
    call alpha(NX, i, 1, u, f_prime, aj1)
    
    w0 = aj0 / (aj0 + aj1)
    w1 = aj1 / (aj0 + aj1)
    Rout = w0 * interpolant(NX, i, x_pt, u, x, dx) + w1 * interpolant(NX, i +1, x_pt, u, x, dx)
    
    return
end subroutine R


subroutine alpha(NX, i, idx, u, f_prime, out)
    
    use iso_fortran_env, only: RK => real64
    use OMP_LIB
    implicit none
    
    integer                , intent(in) :: i, idx, NX
    real(RK), dimension(NX), intent(in) :: u, f_prime
    real(RK),                intent(out):: out
    
    real(RK) :: SI
    real(RK), parameter :: epsilon = 0.000001d0
    
    if (f_prime(i) .gt. 0) then
        if (idx.eq.0) then
            out = 1._RK / 2._RK / (epsilon + SI(u(i), u(i-1)))**2
        else if (idx .eq. 1) then
            out = 1._RK  / (epsilon + SI(u(i+1), u(i)) )**2
        end if
    else if (f_prime(i) .le. 0) then
        if (idx .eq. 0) then
            out = 1._RK / (epsilon + SI(u(i), u(i+1)) )**2
        else if (idx .eq. 1) then
            out = 1._RK / 2._RK / (epsilon + SI(u(i+1), u(i+2)) )**2_RK
        end if
    end if
    
    return
end subroutine alpha

function SI(a, b) result(diff)
    
    use iso_fortran_env, only: RK => real64
    use OMP_LIB
    implicit none
    
    real(RK), intent(in) :: a, b
    real(RK) :: diff
    
    diff = (a - b)**2
    
    return
end function SI

subroutine flux(a, b, x_surf, out)
    
    use iso_fortran_env, only: RK => real64
    use OMP_LIB
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