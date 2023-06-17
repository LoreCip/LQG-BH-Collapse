program PerformTests
    
    use iso_fortran_env, only: RK => real64
    implicit none

    real(RK) :: T

    integer, parameter  :: Ntot = 3
    write(*, "(A23)") "Running all unit tests:"

    call STARToutput(1)
    call test_mathutilsAll(T)
    call ENDoutput(1, Ntot, T, "math_utils")

    call STARToutput(2)
    call test_utilsAll(T)
    call ENDoutput(2, Ntot, T, "utils")

    call STARToutput(3)
    call test_hdf5utilsAll(T)
    call ENDoutput(3, Ntot, T, "hdf5_utils")

end program PerformTests




subroutine STARToutput(i)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer, intent(in) :: i

    write(*, "(A8, 1X, I1, 1X, A20)") "Starting", i

end subroutine STARToutput

subroutine ENDoutput(i, Ntot, time, name)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer       , intent(in) :: i, Ntot
    real(RK)      , intent(in) :: time
    character(*), intent(in) :: name

    character(6) :: outcome

    outcome = "Passed"
    write(*, "(I2, A1, I2, 1X, A6, 1X, A19, 1X, A10, 1X, A6, ES10.2, A8)")&
            i, '/', Ntot, "Test :", trim(name), "--------->", trim(outcome), time, " seconds"
    
end subroutine ENDoutput


! subroutine test_<name>All(T)

!     use iso_fortran_env, only: RK => real64
!     implicit none

!     real(RK), intent(out) :: T
!     real(RK) :: T1, T2

!     call cpu_time(T1)


!     call cpu_time(T2)
!     T = T2 - T1

! end subroutine test_<name>All