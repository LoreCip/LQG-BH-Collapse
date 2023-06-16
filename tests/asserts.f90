module assert_statments
  implicit none
  
contains
  
    subroutine assertArrayEqual(expected, actual, tol, message)
        
        use iso_fortran_env, only: RK => real64
        implicit none
        
        real(RK), dimension(:), intent(in) :: expected, actual
        real(RK), intent(in) :: tol
        character(len=*), intent(in), optional :: message
        
        integer :: i
        
        if (size(expected) /= size(actual)) then
            write(*, '(A)') "Assertion failed: Array sizes do not match!"
            stop 1
        endif
        
        do i = 1, size(expected)
            if (abs(expected(i) - actual(i)) > tol) then
                write(*, '(A)') "Assertion failed: Array elements differ beyond tolerance!"
                if (present(message)) write(*, '(A)') "Additional Info: ", trim(message)
                stop 1
            endif
        end do
      
    end subroutine assertArrayEqual
    
    subroutine assertFloatEqual(expected, actual, tol, message)
    
        use iso_fortran_env, only: RK => real64
        implicit none
        
        real(RK), intent(in) :: expected, actual
        real(RK), intent(in) :: tol
        character(len=*), intent(in), optional :: message
        
        if (abs(expected - actual) > tol) then
            write(*, '(A,ES10.2)') "Assertion failed: Floating-point values differ beyond tolerance ", tol
            write(*, '(A,ES10.2)') "They differ by ", abs(expected - actual)
            write(*, '(A,ES10.2,A,ES10.2)') "Expected: ", expected, "   Actual : ", actual
            if (present(message)) write(*, '(A)') "Additional Info: ", trim(message)
            stop 1
        endif

    end subroutine assertFloatEqual

    subroutine assertIntEqual(expected, actual, message)

        implicit none
        integer, intent(in) :: actual, expected
        character(*), intent(in), optional :: message
         
        if (actual /= expected) then
            write(*, '(A)') "Assertion failed: integer values differ!"
            write(*, *) "They differ by ", abs(expected - actual)
            write(*, *) "Expected: ", expected, "   Actual : ", actual
            if (present(message)) write(*, *) "Additional Info: ", trim(message)
            stop 1
        endif
    end subroutine assertIntEqual

    subroutine assertStrEqual(expected, actual, message)
        implicit none
        character(len=*), intent(in) :: actual, expected
        character(len=*), intent(in), optional :: message
    
        if (actual /= expected) then
            write(*, '(A)') "Assertion failed: strings content differ!"
            write(*, '(A)') "Expected: ", trim(expected)
            write(*, '(A)') "Actual : ", trim(actual)
            if (present(message)) write(*, '(A)') "Additional Info: ", trim(message)
            stop 1
        endif
    end subroutine assertStrEqual

end module assert_statments

