subroutine test_mathutilsAll(T)

    use iso_fortran_env, only: RK => real64
    implicit none

    real(RK), intent(out) :: T
    real(RK) :: T1, T2

    call cpu_time(T1)

    call test_heaviside()
    call test_vb()
    call test_interpolant()

    call cpu_time(T2)
    T = T2 - T1

end subroutine test_mathutilsAll

subroutine test_heaviside()

    use iso_fortran_env, only: RK => real64
    use assert_statments
    implicit none

    real(RK) :: heaviside, x, true
    real(RK), parameter :: tol = 1E-10
    
    x = -5_RK
    true = 0_RK
    call assertFloatEqual(true, heaviside(x), tol, "Heaviside - Test 1")

    x = 0_RK
    true = 0.5_RK
    call assertFloatEqual(true, heaviside(x), tol, "Heaviside - Test 2")

    x = 5_RK
    true = 1_RK
    call assertFloatEqual(true, heaviside(x), tol, "Heaviside - Test 3")

end subroutine test_heaviside

subroutine test_vb()

    use iso_fortran_env, only: RK => real64
    use assert_statments
    implicit none

    real(RK) :: vb, x, u, true
    real(RK), parameter :: tol = 1E-10

    x = 0_RK
    u = 0_RK
    true = 0_RK
    call assertFloatEqual(true, vb(u, x), tol, "Velocity of B field - Test 1")

    x = 1E-12
    u = 1E-24
    true = 0_RK
    call assertFloatEqual(true, vb(u, x), tol, "Velocity of B field - Test 2")

    x = 1_RK
    u = DATAN(1._RK)
    true = 0.5_RK
    call assertFloatEqual(true, vb(u, x), tol, "Velocity of B field - Test 3")

    x = -1_RK
    u = DATAN(1._RK)
    true = -0.5_RK
    call assertFloatEqual(true, vb(u, x), tol, "Velocity of B field - Test 4")

end subroutine test_vb


subroutine test_interpolant()

    use iso_fortran_env, only: RK => real64
    use assert_statments
    implicit none

    real(RK) :: interpolant, true, xpt
    
    real(RK), parameter :: tol = 1E-10, dx = 1._RK
    real(RK), dimension(2), parameter :: x = (/ 0_RK, 1_RK /), u = (/ 1_RK, 2_RK /)

    true = u(1)
    call assertFloatEqual(true, interpolant(2, 2, x(1), u, x, dx), tol, "Interpolant - Test 1")

    true = u(2)
    call assertFloatEqual(true, interpolant(2, 2, x(2), u, x, dx), tol, "Interpolant - Test 2")

    xpt = 0.5_RK
    true = 1.5_RK
    call assertFloatEqual(true, interpolant(2, 2, xpt, u, x, dx), tol, "Interpolant - Test 3")
   
end subroutine test_interpolant