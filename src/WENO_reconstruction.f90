! FILE: WENO_reconstruction.f90

subroutine WENO_RHS(NX, i, u, f_prime, x, dx, nghost, L, Bp, Bm, e_k, e_l)
    
    use iso_fortran_env, only: RK => real64
    implicit none
    
    integer                , intent(in) :: i, NX, nghost
    real(RK)               , intent(in) :: dx
    real(RK), dimension(NX), intent(in) :: x, f_prime
    real(RK), dimension(2*NX), intent(in) :: u
    real(RK),                intent(out):: L, Bp, Bm, e_k, e_l
    
    real(RK) :: x_plus, x_minus, R1, R2, R3, R4, o1, o2
    integer  :: idxA, idxB
    
    x_plus  = x(i) + dx / 2_RK
    x_minus = x(i) - dx / 2_RK
   
    call R(NX, i  , x_plus,  u(1:NX), f_prime, x, dx, R1)
    call R(NX, i+1, x_plus,  u(1:NX), f_prime, x, dx, R2)
    call R(NX, i-1, x_minus, u(1:NX), f_prime, x, dx, R3)
    call R(NX, i  , x_minus, u(1:NX), f_prime, x, dx, R4)

    if ( i.eq.(2 + nghost) ) then
        call flux(R1, R2, x_plus, o1, Bp, idxA) 
        o2 = 0_RK
        Bm = 0
        idxB = -1
    else
        call flux(R1, R2, x_plus , o1, Bp, idxA) 
        call flux(R3, R4, x_minus, o2, Bm, idxB)
    end if
    L = - ( o1 - o2 ) / dx
    
    if (idxA.eq.1) then
        e_k = u(NX + i)
    else if(idxA.eq.-1) then
        e_k = u(NX + i+1)
    else if (idxA.eq.0) then
        e_k = 0
    end if
    
    if (idxB.eq.1) then
        e_l = u(NX + i-1)
    else if(idxB.eq.-1) then
        e_l = u(NX + i)
    else if (idxB.eq.0) then
        e_l = 0
    end if

    return 
end subroutine WENO_RHS


subroutine R(NX, i, x_pt, u, f_prime, x, dx, Rout)
    
    use iso_fortran_env, only: RK => real64
    implicit none
    
    integer                , intent(in) :: i, NX
    real(RK)               , intent(in) :: dx, x_pt
    real(RK), dimension(NX), intent(in) :: u, x, f_prime
    real(RK),                intent(out):: Rout
    
    real(RK) :: aj0, aj1, w0, w1, interpolant
    
    call alpha(NX, i, u, f_prime, aj0, aj1)
    
    w0 = aj0 / (aj0 + aj1)
    w1 = aj1 / (aj0 + aj1)
    Rout = w0 * interpolant(NX, i, x_pt, u, x, dx) + w1 * interpolant(NX, i+1, x_pt, u, x, dx)

    return
end subroutine R


subroutine alpha(NX, i, u, f_prime, out0, out1)
    
    use iso_fortran_env, only: RK => real64
    implicit none
    
    integer                , intent(in) :: i, NX
    real(RK), dimension(NX), intent(in) :: u, f_prime
    real(RK),                intent(out):: out0, out1
    
    real(RK) :: SI
    real(RK), parameter :: epsilon = 1E-7_RK
    
    if (f_prime(i) .gt. 0) then
        out0 = 1._RK / 2._RK / ( epsilon + SI(u(i), u(i-1)) )**2
        out1 = 1._RK         / ( epsilon + SI(u(i+1), u(i)) )**2
    else        
        out0 = 1._RK         / ( epsilon + SI(u(i), u(i-1)) )**2
        out1 = 1._RK / 2._RK / ( epsilon + SI(u(i+1), u(i)) )**2
    end if

    return
end subroutine alpha

function SI(a, b) result(diff)
    
    use iso_fortran_env, only: RK => real64
    implicit none
    
    real(RK), intent(in) :: a, b
    real(RK) :: diff
    
    diff = (a - b)**2
    
    return
end function SI

subroutine flux(a, b, x_surf, out, Bout, idx)
    
    use iso_fortran_env, only: RK => real64
    implicit none
    
    real(RK), intent(in) :: a, b, x_surf
    real(RK), intent(out):: out, Bout
    integer, intent(out) :: idx
    
    real(RK) :: ul, ur, FL, FR
    real(RK), parameter :: PI=4._RK*DATAN(1._RK)

    ul = a / x_surf**2_RK
    ur = b / x_surf**2_RK
    FL = 0.5_RK * x_surf**3 * sin(ul)**2 
    FR = 0.5_RK * x_surf**3 * sin(ur)**2
    
    if ( ul.le.ur ) then
        ! min(FL, FR)
        if ( FL.lt.FR ) then
            out = FL
            Bout = a
            idx = 1
        else
            out = FR
            Bout = b
            idx = -1
        end if
    else if ( ul.gt.ur ) then
        if ( (ur.gt.-PI/2_RK) .or. (ul.lt.-PI/2_RK) ) then
            ! max(FL, FR)
            if ( FL.gt.FR ) then
                out = FL
                Bout = a
                idx = 1
            else
                out = FR
                Bout = b
                idx = -1
            end if
        else
            out = 0.5 * x_surf**3
            Bout = a
            idx = 0
        end if
    end if
    
    return
end subroutine flux