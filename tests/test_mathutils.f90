subroutine test_mathutilsAll(T, status)

    use iso_fortran_env, only: RK => real64
    implicit none

    real(RK), intent(out) :: T
    logical , intent(out) :: status
    real(RK) :: T1, T2

    call cpu_time(T1)


    call cpu_time(T2)
    T = T2 - T1

end subroutine test_mathutilsAll

subroutine test_heaviside(status)

    logical , intent(out) :: status

    call assert(0, heaviside(-5))

end subroutine test_heaviside