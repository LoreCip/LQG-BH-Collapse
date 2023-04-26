program test
    use iso_fortran_env, only: RK => real64

    real(RK) :: x, dx, xxx
    integer :: i, ufile

    dx = real(2_RK) / real(10000_RK)
    xxx = -1

    open(newunit=ufile, file='./test.dat', status='UNKNOWN', POSITION='append')
    do i = 1, 10000
        xxx = xxx + dx
        write(ufile, *) acos(xxx)
    end do
    close(ufile)

end program test