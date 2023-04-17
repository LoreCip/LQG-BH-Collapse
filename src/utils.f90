! FILE: utils.f90

subroutine inputParser(path, T_final, r0, m, r, xM, h, N_output, nthreads)

    use iso_fortran_env, only: RK => real64
    implicit none
    
    character(len=100), intent(in) :: path
    real(RK)          , intent(out):: T_final, r0, m, xM, h
    integer           , intent(out):: r, N_output, nthreads

    real(RK), dimension(8) :: inputs
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
    !!  7    |  N_output
    !!  8    |  nthreads

    open(newunit=nu, file = path, status='old', iostat=ios)
    if ( ios /= 0 ) stop "Error opening file ParameterFile.dat"

    do n = 1, 8
        read(nu, *, iostat=ios)  inputs(n), blabla
        if (ios /= 0) STOP "Error while reading paramters from ParameterFile.dat"
    end do
    close(nu)

    T_final = inputs(1)        
    r0      = inputs(2)   
    m       = inputs(3)  
    r       = int(inputs(4))  
    xM      = inputs(5)   
    h       = inputs(6)
    N_output= int(inputs(7))
    nthreads= int(inputs(8))
    
    return
end subroutine inputParser

subroutine saveOutput(path, NX, xs, u, rho, nghost)
    
    use iso_fortran_env, only: RK => real64
    implicit none

    character(len=100)     , intent(in) :: path
    integer                , intent(in) :: NX, nghost
    real(RK), dimension(NX), intent(in) :: u, xs, rho

    real(RK), dimension(NX-2*nghost) :: ut, xst, rhot
    integer :: i, j
    character(len=100) :: fpath
    
    j = 1
    do i = 1+nghost, Nx-nghost
        ut(j) = u(i)
        xst(j) = xs(i)
        rhot(j) = rho(i)

        j = j+1
    end do

    fpath = trim(path) // '/xs.dat'
    call save(fpath, xst, NX-2*nghost)
    fpath = trim(path) // '/B.dat'
    call save(fpath, ut, NX-2*nghost)
    fpath = trim(path) // '/rho.dat'
    call save(fpath, rhot, NX-2*nghost)

    return
end subroutine saveOutput

subroutine save(path, var, NX)

    use iso_fortran_env, only: RK => real64
    implicit none

    character(len=100)     , intent(in) :: path
    integer                , intent(in) :: NX
    real(RK), dimension(NX), intent(in) :: var

    integer :: ufile        ! Output file
    
    open(newunit=ufile, file=path, status='new', form='unformatted')
    write(ufile) var
    close(ufile)

return
end subroutine save