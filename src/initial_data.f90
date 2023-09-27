subroutine initial_data(NX, x, dx, m, r0, a0, idx, u)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer                , intent(in) :: NX, idx
    real(RK)               , intent(in) :: dx, m, r0, a0
    real(RK), dimension(NX), intent(in) :: x
    real(RK), dimension(2*NX), intent(out):: u

    integer :: i
    real(RK) :: th, Mass, heaviside, p, tmp, xmov
    real(RK), dimension(NX) :: rr, Marray

    real(RK), parameter :: PI=4._RK*DATAN(1._RK), one = 1_RK, zero = 0_RK

    ! Full dynamics, OS closed
    if ( idx .eq. 0 ) then

        ! Physical values
        do i = 1, NX
            th = heaviside(r0 - x(i))
            Mass = m*x(i)**3_RK / r0**3_RK * th + m * (1_RK - th)
            ! E(x)
            u(NX+i) = - x(i)**2 / a0**2 * th - r0**2 / a0**2 * (1_RK - th)
            ! B(x)
            u(i) = - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Mass / x(i)**3_RK - 2_RK * u(NX+i) / x(i)**2_RK)
        end do

    ! Full dynamics, OS flat
    else if ( idx .eq. 1 ) then

        ! Physical values
        do i = 1, NX
            th = heaviside(r0 - x(i))
            Mass = m*x(i)**3_RK / r0**3_RK * th + m * (1_RK - th)
            ! E(x)
            u(NX+i) = 0_RK
            ! B(x)
            u(i) = - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Mass / x(i)**3_RK - 2_RK * u(NX+i) / x(i)**2_RK)
        end do

    ! Full dynamics, OS open
    else if ( idx .eq. 2 ) then

        ! Physical values
        do i = 1, NX
            th = heaviside(r0 - x(i))
            Mass = m*x(i)**3_RK / r0**3_RK * th + m * (1_RK - th)
            ! E(x)
            u(NX+i) = + x(i)**2 / a0**2 * th + r0**2 / a0**2 * (1_RK - th)
            ! B(x)
            u(i) = - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Mass / x(i)**3_RK - 2_RK * u(NX+i) / x(i)**2_RK)
        end do

    ! Full dynamics, starting from arbitrary density profile
    else if ( idx .eq. 3 ) then

        ! Unnormalized density function
        do i = 1, NX
            rr(i) = exp(-2_RK*(x(i) - 7.5_RK)**2)/5_RK + exp(-4_RK*(x(i) - 17.5_RK)**2)/10_RK 
        end do
        
        ! Unnormalized mass function
        Marray(1) = 0_RK
        Marray(2) = 0_RK
        p = 0_RK
        do i = 3, NX-1
            p = p + (rr(i-1)*x(i-1)**2 + rr(i)*x(i)**2)
            Marray(i) = 4_RK*PI * p * dx * 0.5_RK
        end do
        
        ! Normalized mass function
        do i = 1, NX
            Marray(i) =  m * Marray(i) / Marray(NX-1)
        end do
        Marray(NX) = Marray(NX-1)

        ! Physical values
        xmov = x(1)

        do i = 2, NX-1
            th = heaviside(r0 - x(i))

            ! Cubic interpolation around r0 to ensure epsilon is a C1 function
            if ((x(i).le.(r0+1_RK)).and.(x(i).ge.(r0-1_RK))) then
                u(NX+i) = r0**2 * x(i)**3 - r0**2 *(2_RK + r0 + 2_RK*xmov) * x(i)**2 - r0**2 * (1_RK+r0) *(-1_RK+r0-4_RK*xmov) * x(i)
                u(NX+i) = u(NX+i) + r0**2 * (r0**3 - 2_RK * r0**2 * (1_RK+xmov) - 2_RK*xmov * (1_RK + 2_RK*xmov) + r0*(1_RK+4_RK*xmov))
                u(NX+i) = u(NX+i) / (4_RK*a0**2) / (r0 - xmov)**2
            else
                u(NX+i) =  - r0**2 / (r0 - xmov)**2 * (x(i) - xmov)**2 / a0**2 * th - r0**2 / a0**2 * (1_RK - th)
            end if

            ! B(x)
            tmp = 1_RK - 4_RK * Marray(i) / x(i)**3_RK - 2_RK * u(NX+i) /  x(i)**2
            if (tmp.gt.1_RK) then
                u(NX+1:NX+i) = 0_RK
                u(:i) = - 0.5_RK*x(:i)**2 * acos(1_RK - 4_RK * Marray(:i) / x(:i)**3_RK)
                xmov = x(i)
            else
                u(i) = - 0.5_RK*x(i)**2 * acos(tmp)
            end if
        end do

        ! Post bounce dynamics, starting from step function
    else if ( idx .eq. 4 ) then

        ! Physical values
        do i = 1, NX-1

            th = heaviside(r0 - x(i))
            Mass = m*x(i)**3_RK / r0**3_RK * th + m * (1_RK - th)
            ! E(x)
            u(NX+i) = - x(i)**2 / a0**2 * th - r0**2 / a0**2 * (1_RK - th)
            ! B(x)
            u(i) =  - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * Mass / x(i)**3_RK - 2_RK * u(NX+i) / x(i)**2_RK) * (1_RK - th)  &
                    + x(i)**2 * ( - PI + asin(sqrt(2_RK * Mass / x(i)**3_RK + u(NX+i) / x(i)**2_RK)) ) * th
        end do

    ! Post bounce dynamics, starting from delta function
    else if ( idx .eq. 5 ) then

        ! Physical values
        do i = 1, NX-1
            th = heaviside(r0 - x(i))
            ! E(x)
            u(NX+i) = - r0**2 / a0**2 * (1_RK - th)
            ! B(x)
            u(i) =  - 0.5_RK*x(i)**2 * acos(1_RK - 4_RK * m / x(i)**3_RK - 2_RK * u(NX+i) / x(i)**2_RK) * (1_RK - th)  &
                    - PI * x(i)**2 * th
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