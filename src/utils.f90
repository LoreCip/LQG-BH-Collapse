! FILE: utils.f90

subroutine interp(n, x, xi, ui, out)
    ! Linear 1D interpolator
    ! ui and xi are 1D arrays with n entries corresponding to the points needed to compute the n-th order interpolating polynomial.
    use, intrinsic :: iso_fortran_env, dp=>real64

    real(dp), parameter :: EPS = 1E-16

    integer               , intent(in) :: n
    real(dp), dimension(n), intent(in) :: ui, xi
    real(dp)              , intent(in) :: x
    real(dp)              , intent(out):: out 

    integer  :: j
    real(dp) :: Lnj, ssum, xx

    xx = x
    do j = 1, n
        if (x == xi(j)) then
            xx = x + EPS
            exit
        end if
    end do

    ssum = 0
    do j = 0, n
        call prod(n, j, xi, xx, Lnj)
        ssum = ssum + ui(j) * Lnj
    end do

    out = ssum
    return
end subroutine

subroutine prod(n, j, xi, x, out)

    use, intrinsic :: iso_fortran_env, dp=>real64

    integer               , intent(in)  :: n, j
    real(dp)              , intent(in)  :: x
    real(dp), dimension(n), intent(in)  :: xi
    real(dp)              , intent(out) :: out

    integer :: k
    real(dp):: p
    
    p = 1
    do k = 0, n
        if (k == j) then; cycle; end if
        p = p * (x - xi(k)) / (xi(j) - xi(k))
    end do

    out = p
    return
end subroutine


subroutine minimizer()

    use, intrinsic :: iso_fortran_env, dp=>real64

    

    return
end subroutine