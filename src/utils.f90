! FILE: utils.f90

subroutine inputParser(path, T_final, r0, m, r, xM, h, N_save, N_output, nthreads)

    use iso_fortran_env, only: RK => real64
    implicit none
    
    character(len=100), intent(in) :: path
    real(RK)          , intent(out):: T_final, r0, m, xM, h
    integer           , intent(out):: r, N_save, N_output, nthreads

    real(RK), dimension(9) :: inputs
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
    !!  7    |  N_save
    !!  8    |  N_output
    !!  9    |  nthreads

    open(newunit=nu, file = path, status='old', iostat=ios)
    if ( ios /= 0 ) stop "Error opening file ParameterFile.dat"

    do n = 1, 9
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
    N_save  = int(inputs(7))
    N_output= int(inputs(8))
    nthreads= int(inputs(9))
    
    return
end subroutine inputParser

subroutine saveOutput(path, NX, u, rho, nghost)
    
    use iso_fortran_env, only: RK => real64
    implicit none

    character(len=100)     , intent(in) :: path
    integer                , intent(in) :: NX, nghost
    real(RK), dimension(NX), intent(in) :: u, rho

    character(len=100) :: fpath

    fpath = trim(path) // '/B.dat'
    call save(fpath, u(nghost+1:NX-nghost-1), NX-2*nghost)
    fpath = trim(path) // '/rho.dat'
    call save(fpath, rho(nghost+1:NX-nghost-1), NX-2*nghost)

    return
end subroutine saveOutput

subroutine save(path, var, NX)

    use iso_fortran_env, only: RK => real64
    implicit none

    character(len=100)     , intent(in) :: path
    integer                , intent(in) :: NX
    real(RK), dimension(NX), intent(in) :: var

    integer :: i, ufile        ! Output file
    
    open(newunit=ufile, file=path, status='UNKNOWN', POSITION='append')
    write(ufile, *) (var(i), i = 1, NX)
    close(ufile)

return
end subroutine save