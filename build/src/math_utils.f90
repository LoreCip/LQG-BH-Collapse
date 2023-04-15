    subroutine func(NX, u, x, fun)
    
        use iso_fortran_env, only: RK => real64
        use OMP_LIB
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
        use OMP_LIB
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
    
    
    function heaviside(x) result(out)
    
        use iso_fortran_env, only: RK => real64
        use OMP_LIB
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