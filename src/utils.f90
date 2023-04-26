! FILE: utils.f90

subroutine inputParser(path, T_final, r0, a0, m, r, xM, h, N_save, N_output, nthreads)

    use iso_fortran_env, only: RK => real64
    implicit none
    
    character(len=100), intent(in) :: path
    real(RK)          , intent(out):: T_final, r0, a0, m, xM, h
    integer           , intent(out):: r, N_save, N_output, nthreads

    real(RK), dimension(10) :: inputs
    integer :: ios, n, nu
    character(1024) :: blabla

    !! inputs contains, in order:
    !! idx   |   var:
    !! --------------
    !!  1    |  T_final
    !!  2    |  r0
    !!  3    |  a0
    !!  4    |  m
    !!  5    |  r
    !!  6    |  xM
    !!  7    |  h
    !!  8    |  N_save
    !!  9    |  N_output
    !!  10   |  nthreads

    open(newunit=nu, file = path, status='old', iostat=ios)
    if ( ios /= 0 ) stop "Error opening file ParameterFile.dat"

    do n = 1, 10
        read(nu, *, iostat=ios)  inputs(n), blabla
        if (ios /= 0) STOP "Error while reading paramters from ParameterFile.dat"
    end do
    close(nu)

    T_final = inputs(1)        
    r0      = inputs(2)
    a0      = inputs(3)
    m       = inputs(4)  
    r       = int(inputs(5))  
    xM      = inputs(6)   
    h       = inputs(7)
    N_save  = int(inputs(8))
    N_output= int(inputs(9))
    nthreads= int(inputs(10))
    
    return
end subroutine inputParser

subroutine saveOutput(path, NX, B, E, rho, nghost)
    
    use iso_fortran_env, only: RK => real64
    implicit none

    character(len=*)     , intent(in) :: path
    integer                , intent(in) :: NX, nghost
    real(RK), dimension(NX), intent(in) :: B, E, rho

    character(len=1024) :: fpath

    fpath = trim(path) // '/B.dat'
    call save(fpath, B(nghost+1:NX-nghost), NX-2*nghost)
    fpath = trim(path) // '/E.dat'
    call save(fpath, E(nghost+1:NX-nghost), NX-2*nghost)
    fpath = trim(path) // '/rho.dat'
    call save(fpath, rho(nghost+1:NX-nghost), NX-2*nghost)

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