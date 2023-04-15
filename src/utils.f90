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
    
        return
    end subroutine inputParser

    
