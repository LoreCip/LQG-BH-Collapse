program LQGeq

    use OMP_LIB
    use iso_fortran_env, only: RK => real64

    implicit none

    real(RK) :: T1, T2

    integer :: error_code, VERBOSE

    ! Physical paramters
    real(RK) :: T_final, &   ! Max integration time
                    r0,      &   ! Density cutoff
                    m,       &   ! Total mass
                    t,       &   ! Integrated time
                    dt           ! Time step
    real(RK), dimension(:), allocatable :: &
                    xs,      &   ! Grid
                    u,       &   ! Solution
                    u_p          ! Previous solution

    ! Computational parameters
    real(RK) :: xM,      &   ! Outer boundary of the grid
                    h            ! Grid spacing
    integer  :: r,       &   ! Order of the WENO method
                nghost,  &   ! Number of ghost cells
                NX           ! Number of points in xs
    
    real(RK), parameter :: eps = 0.0001_RK 
    
    ! Iterables
    integer :: i

    call cpu_time(T1)

    VERBOSE = 0
    
    T_final = 2 ! 8.5*mass**2 + 40*mass
    r0 = 15
    m = 5

    r = 2

    xM = 50.0_RK
    h = 0.1_RK
    nghost = 1

    NX = int(xM / h + 1 + 2*nghost)
    allocate(xs(NX), u(NX), u_p(NX), stat=error_code)
    if(error_code /= 0) STOP "Error during array allocations!"

    !$OMP PARALLEL DO
    do i = 1, NX
      xs(i) = (eps + i - nghost) * h
    end do
    !$OMP END PARALLEL DO

    call initial_data(NX, xs, m, r0, u_p)
    t = 0
    do while (t.lt.T_final)

        call CLF(NX, u_p, xs, h, dt)

        if ((t + dt).gt.T_final) then
            dt = T_final - t
            t = T_final
        else
            t = t + dt
        end if

        call TVD_RK(NX, u_p, xs, h, dt, nghost, u)

        if ( t.lt.T_final ) then
            do i = 1, NX
                u_p(i) = u(i)
            end do
        end if

    end do

    call cpu_time(T2)
    print*, "Total time:", T2 - T1, "seconds."
    print*, u

end program LQGeq