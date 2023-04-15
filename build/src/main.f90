program LQGeq
    
    use iso_fortran_env, only: RK => real64
    use OMP_LIB
    implicit none

    real(RK) :: T1, T2
    integer  :: iTimes1,iTimes2, rate

    integer :: error_code

    ! Physical parameters
    real(RK), parameter :: eps = 0.0001_RK, PI=4._RK*DATAN(1._RK)

    real(RK) :: T_final, &   ! Max integration time
                    r0,      &   ! Density cutoff
                    m,       &   ! Total mass
                    t,       &   ! Integrated time
                    dt           ! Time step
    real(RK), dimension(:), allocatable :: &
                    xs,      &   ! Grid
                    u,       &   ! Solution
                    u_p,     &   ! Previous solution
                    rho          ! Density profile

    ! Computational parameters
    real(RK) :: xM,      &   ! Outer boundary of the grid
                h            ! Grid spacing
    integer  :: r,       &   ! Order of the WENO method
                nghost,  &   ! Number of ghost cells
                NX           ! Number of points in xs
    
    integer :: i, num_args
    character(len=100), dimension(2) :: args

    CALL system_clock(count_rate=rate)
    call cpu_time(T1)
    call SYSTEM_CLOCK(iTimes1)

    ! Read configuration parameters
    num_args = command_argument_count()
    if (num_args .gt. 2) STOP "Only two args are contemplated, the path of the configuration file and the path for the output!"
    do i = 1, 2
        call get_command_argument(i,args(i))
        if (args(i) .eq. '') then
            if (i.eq.1) then
                args(i) = 'build/ParameterFile.dat'
            else if (i .eq. 2) then
                args(i) = 'build/outputs'
            end if
        end if
    end do
    call inputParser(args(1), T_final, r0, m, r, xM, h)

    ! Perform consistency checks
    if (r .eq. 2) then
        nghost = 1
    else
        STOP "Only r = 2 is implemented."
    end if

    ! Define numerical grid
    NX = int(xM / h + 1 + 2*nghost)
    allocate(xs(NX), u(NX), u_p(NX), rho(NX), stat=error_code)
    if(error_code /= 0) STOP "Error during array allocations!"

    !$OMP PARALLEL DO
    do i = 1, NX
      xs(i) = (eps + i - nghost) * h
    end do
    !$OMP END PARALLEL DO

    ! Produce initial data
    call initial_data(NX, xs, m, r0, u_p)
    
    ! Time evolution
    t = 0
    do while (t.lt.T_final)

        ! Determine dt from Courant condition
        call CLF(NX, u_p, xs, h, dt)

        if ((t + dt).gt.T_final) then
            dt = T_final - t
            t = T_final
        else
            t = t + dt
        end if

        ! Perform time step
        call TVD_RK(NX, u_p, xs, h, dt, nghost, u)

        ! Update previous step
        if ( t.lt.T_final ) then
            do i = 1, NX
                u_p(i) = u(i)
            end do
        end if

    end do

    do i = 1, NX
        rho(i) = (u_p(i) - u(i)) / (4*PI * xs(i)**2 * dt)
    end do
    
    call cpu_time(T2)
    call SYSTEM_CLOCK(iTimes2)
    print*, "Total CPU time:", T2 - T1, "seconds."
    print*, "Total system time", real(iTimes2-iTimes1)/real(rate)

    call saveOutput(args(2), NX, xs, u, rho, nghost)

    deallocate(xs, u, u_p, rho)


end program LQGeq