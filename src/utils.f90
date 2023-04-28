! FILE: utils.f90

subroutine inputParser(path, T_final, r0, a0, m, r, xM, h, N_save, N_output, nthreads)

    use iso_fortran_env, only: RK => real64
    implicit none
    
    character(len=100), intent(in) :: path
    real(RK)          , intent(out):: T_final, r0, a0, m, xM, h
    integer           , intent(out):: r, N_save, N_output, nthreads

    real(RK), dimension(10) :: inputs
    integer :: ios, n, nu, error_code
    character(1024) :: blabla, error_string

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
    if ( ios /= 0 ) stop "Error opening parameter file."

    do n = 1, 10
        read(nu, *, iostat=ios)  inputs(n), blabla
        if (ios /= 0) STOP "Error while reading parameters from parameter file."
    end do
    close(nu)

    call paramChecker(inputs, error_code, error_string)
    if (error_code /= 0) then
        write(*,*) "Wrong parameters in parameter file."
        write(*,*) trim(error_string)
        STOP error_code
    end if
    
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



subroutine paramChecker(inputs, error_code, error_string)

    use iso_fortran_env, only: RK => real64
    implicit none

    real(RK), dimension(10), intent(in)  :: inputs
    integer                , intent(out) :: error_code
    character(1024)        , intent(out) :: error_string

    real(RK) :: tmp

    error_code = 0

    ! Physical conditions
    ! Check for:
    ! r0 < a0
    if (inputs(3).lt.inputs(2)) then
        error_code = 1
        error_string = "Characteristic radius r0 must be smaller than scale factor a0."
        return
    end if

    ! xmax < 2 * m * a0^2 / r0^2
    tmp = 2_RK * inputs(4) * inputs(3)**2 / inputs(2)**2
    if (tmp.lt.inputs(6)) then
        error_code = 2
        error_string = "Furthest grid point xM must be smaller than 2 * m * a0^2 / r0^2."
        return
    end if

    ! Sanity checks:
    if ( any(inputs.lt.0) ) then
        error_code = 3
        error_string = "All paramteres must be non negative."
        return
    end if
        
    ! r = 2
    if (int(inputs(5)).ne.2) then
        error_code = 4
        error_string = "Only WENO3 (r=2) is implemtented."
        return
    end if

    return
end subroutine paramChecker


subroutine saveOutput(path, NX, B, E, rho, t, dt, BHpresent, nghost)
    
    use iso_fortran_env, only: RK => real64
    implicit none

    character(len=*)       , intent(in) :: path
    integer                , intent(in) :: NX, nghost
    real(RK), dimension(NX), intent(in) :: B, E, rho
    real(RK),                intent(in) :: t, dt
    logical,                 intent(in) :: BHpresent

    character(len=1024) :: fpath
    real(RK) :: logic2dbl
    real(RK), dimension(3) :: vvv

    fpath = trim(path) // '/B.dat'
    call save(fpath, B(nghost+1:NX-nghost), NX-2*nghost)
    
    fpath = trim(path) // '/E.dat'
    call save(fpath, E(nghost+1:NX-nghost), NX-2*nghost)
    
    fpath = trim(path) // '/rho.dat'
    call save(fpath, rho(nghost+1:NX-nghost), NX-2*nghost)
    
    fpath = trim(path) // '/times.dat'
    vvv = (/t, dt, logic2dbl(BHpresent) /)
    call save(fpath, vvv, size(vvv))

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


double precision function logic2dbl(a)
    
    use iso_fortran_env, only: RK => real64
    implicit none

    logical, intent(in) :: a

    if (a) then
      logic2dbl = 1._RK
    else
      logic2dbl = 0._RK
    end if

end function logic2dbl