program PerformTests
    
    use iso_fortran_env, only: RK => real64
    implicit none

    logical, parameter :: passed = .true.

    real(RK) :: T
    logical  :: status

    integer, parameter  :: Ntot = 5
    write(*, "(A23)") "Running all unit tests:"

    call STARToutput(1)
    call test_mathutilsAll(T, status)
    call ENDoutput(1, Ntot, T, status, "math_utils")

end program PerformTests

subroutine STARToutput(i)

    use iso_fortran_env, only: RK => real64
    implicit none

    integer       , intent(in) :: i

    write(*, "(A8, 1X, I1, 1X, A20)") "Starting", i

end subroutine STARToutput

subroutine ENDoutput(i, Ntot, time, status, name)

    use iso_fortran_env, only: RK => real64
    implicit none

    logical, parameter :: passed = .true.

    integer       , intent(in) :: i, Ntot
    real(RK)      , intent(in) :: time
    logical       , intent(in) :: status
    character(*), intent(in) :: name

    character(6) :: outcome

    if (status.eqv.passed) then
        outcome = "Passed"
    else
        outcome = "Failed"
    end if

    write(*, "(I2, A1, I2, 1X, A6, 1X, A19, 1X, A10, 1X, A6, F7.3, A7)")&
            i, '/', Ntot, "Test :", trim(name), "--------->", trim(outcome), time, " seconds"
    
end subroutine ENDoutput


! subroutine test_<name>All(T, status)

!     use iso_fortran_env, only: RK => real64
!     implicit none

!     real(RK), intent(out) :: T
!     logical , intent(out) :: status
!     real(RK) :: T1, T2

!     call cpu_time(T1)


!     call cpu_time(T2)
!     T = T2 - T1

! end subroutine test_<name>All