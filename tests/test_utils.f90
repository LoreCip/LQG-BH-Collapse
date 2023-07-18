subroutine test_utilsAll(T)

    use iso_fortran_env, only: RK => real64
    implicit none

    real(RK), intent(out) :: T
    real(RK) :: T1, T2

    call cpu_time(T1)

    call test_paramChecker()

    call cpu_time(T2)
    T = T2 - T1

end subroutine test_utilsAll

subroutine test_logic2dbl()

    use iso_fortran_env, only: RK => real64
    use assert_statments
    implicit none

    real(RK):: logic2dbl, true
    real(RK), parameter :: tol = 1E-10

    true = 1_RK
    call assertFloatEqual(true, logic2dbl(.true.), tol, "Logic to double converter - Test 1")

    true = 0_RK
    call assertFloatEqual(true, logic2dbl(.false.), tol, "Logic to double converter - Test 2")

end subroutine test_logic2dbl


subroutine test_paramChecker()

    use iso_fortran_env, only: RK => real64
    use assert_statments
    implicit none

    integer :: error_code
    real(RK), dimension(11) :: inputs
    character(len=1024) :: error_string, true_error_string
    real(RK), dimension(11), parameter :: true_inputs = (/ 0., 1.0, 15.0, 40.0, 5.0 , 2.0, 70.0, 0.01, 5000.0, 1000., 4. /)

    inputs = true_inputs
    call paramChecker(inputs, error_code, error_string)
    call assertIntEqual(0, error_code, "ParamChecker - Test 0.1")
    true_error_string = ""
    call assertStrEqual(true_error_string, error_string, "ParamChecker - Test 0.2")
    
    inputs = true_inputs
    inputs(1) = 500
    call paramChecker(inputs, error_code, error_string)
    call assertIntEqual(4, error_code, "ParamChecker - Test 1.1")
    true_error_string = "Initial data index not valid. id = 0,1,2,3,4,5"
    call assertStrEqual(true_error_string, error_string, "ParamChecker - Test 1.2")

    inputs = true_inputs
    inputs(2) = -500
    call paramChecker(inputs, error_code, error_string)
    call assertIntEqual(3, error_code, "ParamChecker - Test 2.1")
    true_error_string = "All parameters must be non negative."
    call assertStrEqual(true_error_string, error_string, "ParamChecker - Test 2.2")

    inputs = true_inputs
    inputs(3) = inputs(4) + 10
    call paramChecker(inputs, error_code, error_string)
    call assertIntEqual(1, error_code, "ParamChecker - Test 3.1")
    true_error_string = "Characteristic radius r0 must be smaller than scale factor a0."
    call assertStrEqual(true_error_string, error_string, "ParamChecker - Test 3.2")

    inputs = true_inputs
    inputs(4) = inputs(3) - 10
    call paramChecker(inputs, error_code, error_string)
    call assertIntEqual(1, error_code, "ParamChecker - Test 4.1")
    true_error_string = "Characteristic radius r0 must be smaller than scale factor a0."
    call assertStrEqual(true_error_string, error_string, "ParamChecker - Test 4.2")

    inputs = true_inputs
    inputs(6) = 5
    call paramChecker(inputs, error_code, error_string)
    call assertIntEqual(5, error_code, "ParamChecker - Test 5.1")
    true_error_string = "Only WENO3 (r=2) is implemented."
    call assertStrEqual(true_error_string, error_string, "ParamChecker - Test 5.2")

    inputs = true_inputs
    inputs(7) = 500
    call paramChecker(inputs, error_code, error_string)
    call assertIntEqual(2, error_code, "ParamChecker - Test 6.1")
    true_error_string = "For a closed universe, the furthest grid point xM must be smaller than 2 * m * a0^2 / r0^2."
    call assertStrEqual(true_error_string, error_string, "ParamChecker - Test 6.2")

end subroutine test_paramChecker

