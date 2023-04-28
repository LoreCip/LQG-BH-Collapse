program LQGeq
    
    use iso_fortran_env, only: RK => real64
    use OMP_LIB
    implicit none

    real(RK) :: xxx
    integer  :: iTimes1, iTimes2, rate

    integer :: error_code

    ! Physical parameters
    real(RK), parameter :: eps = 0.0001_RK, PI=4._RK*DATAN(1._RK)

    real(RK) :: T_final, &   ! Max integration time
                    r0,      &   ! Density cutoff
                    a0,      &   ! FRW param
                    m,       &   ! Total mass
                    t,       &   ! Integrated time
                    dt,      &   ! Time step
                    Tformation, &! Time of BH formation
                    Texplosion   ! Time of BH explosion

    real(RK), dimension(:), allocatable :: &
                    xs,      &   ! Grid
                    u,       &   ! Solution
                    u_p,     &   ! Previous solution
                    rho,     &   ! Density profile
                    theta,   &   ! Expansion  
                    vvv
    ! Computational parameters
    real(RK) :: xM,      &   ! Outer boundary of the grid
                h,       &   ! Grid spacing
                ssum        
    integer  :: r,       &   ! Order of the WENO method
                nghost,  &   ! Number of ghost cells
                NX,      &   ! Number of points in xs
                N_output,&   ! Print output every
                N_save,  &   ! Save every
                nthreads     ! Number of threads for OpenMP
                
    ! Iterators           
    integer :: i, counter
    logical :: done, BHpresent

    ! Input arguments
    integer :: num_args
    character(len=1024), dimension(2) :: args
    ! Save
    character(len=1024) :: fpath

    CALL system_clock(count_rate=rate)
    call SYSTEM_CLOCK(iTimes1)

    ! Read configuration parameters
    num_args = command_argument_count()
    if (num_args .gt. 2) STOP "Only two args are contemplated, the path of the configuration file and the path for the output!"
    do i = 1, 2
        call get_command_argument(i,args(i))
        if (args(i) .eq. '') then
            if (i.eq.1) then
                args(i) = 'ParameterFile.par'
            else if (i .eq. 2) then
                args(i) = 'outputs'
            end if
        end if
    end do
    
    call inputParser(args(1), T_final, r0, a0, m, r, xM, h, N_save, N_output, nthreads)

    CALL OMP_SET_NUM_THREADS(nthreads)

    ! PRINT SUMMARY AND NICE OUTPUT
    write(*, "(A65)") "-----------------------------------------------------------------"
    write(*, "(A12)") "WENO3 SOLVER"
    write(*, "(A65)") "-----------------------------------------------------------------"
    write(*, "(A7)") "SUMMARY"
    write(*, "(A35, F9.2)") "   - Simulation time      :    ", T_final
    write(*, "(A35, E9.3)") "   - Grid spacing         :    ", h
    write(*, "(A35, F9.2)") "   - Total mass, M0       :    ", m
    write(*, "(A35, F9.2)") "   - Characteristic radius:    ", r0
    write(*, "(A35, F9.2)") "   - Initial scale factor :    ", a0
    write(*, "(A35, I9)")   "   - Number of threads    :    ", nthreads
    write(*, "(A65)") "-----------------------------------------------------------------"
    write(*, "(A48)") "              Starting simulation...            "
    write(*, "(A65)") "-----------------------------------------------------------------"
    write(*, "(A65)") "    Time     ||    Iteration   ||   BH present   ||    M - M0    "

    ! Array allocation 
    NX = int(xM / h + 1 + 2*nghost)
    allocate(xs(NX), u(2*NX), u_p(2*NX), rho(NX), theta(NX), stat=error_code)
    if(error_code /= 0) STOP "Error during array allocations!"

    
    ! Define numerical grid
    do i = 1, NX
        xs(i) = (eps + i - 1 - nghost) * h
    end do
    if (N_save .gt. 0) then
        fpath = trim(args(2)) // '/xs.dat'
        call save(fpath, xs(nghost+1:NX-nghost), NX-2*nghost)
    end if

    ! Produce initial data
    call initial_data(NX, xs, h, m, r0, a0, 0, u_p)
    
    counter = 0
    ! Time evolution
    t = 0_RK
    done = .false.
    BHpresent = .false.
    do while (t.lt.T_final)

        ! Determine dt from Courant condition
        call CLF(NX, u_p(1:NX), xs, h, dt)

        if ((t + dt).gt.T_final) then
            dt = T_final - t
            t = T_final
            done = .true.
        else
            t = t + dt
        end if      

        ! Perform time step
        call TVD_RK(NX, u_p, xs, h, dt, nghost, u)
        call CompExpansion(NX, u(1:NX), u(NX+1:2*NX), xs, theta)

        if ( (BHpresent.eqv..false.) .and. any(theta.lt.0) ) then
            BHpresent = .true.
            Tformation = t
        else if ( (BHpresent.eqv..true.) .and. all(theta.gt.0) ) then
            BHpresent = .false.
            Texplosion = t
        end if

        ! TERMINAL OUTPUT
        if ( (N_output.gt.0) .and. (mod(counter, N_output).eq.0) ) then
            call compRho(NX, h, dt, u(1:NX), u_p(1:NX), u(NX+1:2*NX), xs, rho)
            ! COMPUTE MASS
            ssum = 0_RK
            do i = 3, NX-1
                ssum = ssum + (rho(i-1) * xs(i-1)**2 + rho(i) * xs(i)**2)
            end do
            ssum = 4 * PI * ssum * h * 0.5_RK
            ! print*, ssum
            write(*, "(2x, F8.3, 3x, A2, 3x, I9, 4x, A2, 7x, l1, 8x, A2, 3x, E11.3)") &
                    t, "||",  counter+1,  "||", BHpresent, "||",  ssum - m
        end if

        if ( (N_save.gt.0) .and. ((mod(counter, N_save).eq.0) .or. done) ) then

            call compRho(NX, h, dt, u(1:NX), u_p(1:NX), u(NX+1:2*NX), xs, rho)
            call saveOutput(args(2), NX, u(1:NX), u(NX+1:2*NX), rho, t, dt, BHpresent, nghost)

            ! write(*, "(A65)") "-----------------------------------------------------------------"
            ! write(*, "(A33, 1X, F9.4)") "Output saved at simulation time", t
            ! write(*, "(A65)") "-----------------------------------------------------------------"

        end if
        counter = counter + 1

        ! Update previous step
        if ( t.lt.T_final ) then
            do i = 1, 2*NX
                u_p(i) = u(i)
            end do
        end if
        
    end do

    deallocate(xs, u, u_p, rho, theta)

    if (N_save .gt. 0) then
        allocate(vvv(4))
        fpath = trim(args(2)) // '/details.dat'
        vvv = (/ m, h, Tformation, Texplosion /)
        call save(fpath, vvv, size(vvv))
        deallocate(vvv)
    end if

    write(*, "(A65)") "-----------------------------------------------------------------"
    write(*, "(A48)") "                  All done!                     "
    write(*, "(A65)") "-----------------------------------------------------------------"

    call SYSTEM_CLOCK(iTimes2)
    xxx = real(iTimes2-iTimes1)/real(rate)
    write(*, "(A21,1x, F7.2, A9)") "Total system runtime:", xxx, " seconds."

end program LQGeq