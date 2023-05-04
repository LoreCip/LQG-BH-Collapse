! FILE: utils.f90

subroutine inputParser(path, T_final, r0, a0, m, r, xM, h, N_save, N_output, nthreads)

    use iso_fortran_env, only: RK => real64
    implicit none
    
    character(len=4096), intent(in) :: path
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


subroutine openOutput(d, ufiles, path)

    implicit none

    integer,               intent(in) :: d
    integer, dimension(d), intent(in) :: ufiles
    character(len=*),      intent(in) :: path

    character(len=4096)              :: fpath
    character(len=100), dimension(4) :: string_array

    integer :: i

    ! Assign values to string array
    string_array(1) = '/B.dat'
    string_array(2) = '/E.dat'
    string_array(3) = '/rho.dat'
    string_array(4) = '/times.dat'

    do i = 1, d
        fpath = path // string_array(i)
        open(unit=ufiles(i), file=fpath, status='new', POSITION='append')
    end do
    
end subroutine openOutput

subroutine saveOutput(ufile, d, var)
    
    use iso_fortran_env, only: RK => real64
    implicit none

    integer ,               intent(in) :: ufile, d
    real(RK), dimension(d), intent(in) :: var

    integer :: i

    write(ufile, *) (var(i), i = 1, d)

    return
end subroutine saveOutput

subroutine closeOutput(d, ufiles)

    implicit none

    integer,               intent(in) :: d
    integer, dimension(d), intent(in) :: ufiles

    integer :: i
    
    do i = 1, d
        close(ufiles(i))
    end do
    
end subroutine closeOutput


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