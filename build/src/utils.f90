! FILE: utils.f90

    subroutine inputParser(path, T_final, r0, m, r, xM, h)

        use iso_fortran_env, only: RK => real64
        use OMP_LIB
        implicit none
        
        character(len=100), intent(in) :: path
        real(RK)          , intent(out):: T_final, r0, m, xM, h
        integer           , intent(out):: r

        real(RK), dimension(6) :: inputs
        integer :: ios, n, nu
        character(100) :: blabla
    
        !! inputs contains, in order:
        !! idx   |   var:
        !! --------------
        !!  1    |  T_final
        !!  2    |  r0
        !!  3    |  m
        !!  4    |  r
        !!  5    |  xM
        !!  6    |  h
        open(newunit=nu, file = path, status='old', iostat=ios)
        if ( ios /= 0 ) stop "Error opening file ParameterFile.dat"
        do n = 1, 6
            read(nu, "(F5.3, A100)", iostat=ios)  inputs(n), blabla
            if (ios /= 0) STOP "Error while reading paramters from ParameterFile.dat"
        end do
        close(nu)

        T_final = inputs(1)        
        r0      = inputs(2)   
        m       = inputs(3)  
        r       = int(inputs(4))  
        xM      = inputs(5)   
        h       = inputs(6)

        print*, T_final,r0     ,m      ,r      ,xM     ,h      
    
        return
    end subroutine inputParser

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
