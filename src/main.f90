program LQGeq

    #ifdef HDF
        use hdf5
    #endif
    
        use iso_fortran_env, only: RK => real64
        use OMP_LIB
        implicit none
    
        real(RK) :: systemtime
        integer  :: iTimes1, iTimes2, rate
    
        integer :: error_code
    #ifdef HDF
        integer(hid_t) :: file_id
    #endif
        integer, dimension(4) :: ufiles
    
        ! Physical parameters
        integer, parameter  :: nghost = 1 ! Number of ghost cells
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
                        u1, u2, u12, u32, &
                        f_prime, &
                        vel,     &
                        e_der,   &
                        rho,     &   ! Density profile
                        theta,   &   ! Expansion  
                        vvv
        ! Computational parameters
        real(RK) :: xM,      &   ! Outer boundary of the grid
                    h,       &   ! Grid spacing
                    ssum,    &
                    logic2dbl
    
        real(RK), dimension(3) :: ttt
    
        integer  :: r,       &   ! Order of the WENO method
                    NX,      &   ! Number of points in xs
                    N_output,&   ! Print output every
                    N_save,  &   ! Save every
                    nthreads,&   ! Number of threads for OpenMP
                    id           ! Initial data index
        ! Iterators           
        integer :: i, counter
        logical :: done, BHpresent, printO, saveO
    
        ! Input arguments
        integer :: num_args
        character(len=4096), dimension(2) :: args
        ! Save
        character(len=4096) :: fpath
    
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
        
        call inputParser(args(1), id, T_final, r0, a0, m, r, xM, h, N_save, N_output, nthreads)
    
        !$ call OMP_SET_NUM_THREADS(nthreads)
    
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
        allocate(e_der(NX), vel(NX), u1(2*NX), u2(2*NX), u12(2*NX), u32(2*NX), f_prime(NX), stat=error_code)
        if(error_code /= 0) STOP "Error during array allocations!"
    
        
        ! Define numerical grid
        do i = 1, NX
            xs(i) = (eps + i - 1 - nghost) * h
        end do
        if (N_save .gt. 0) then
    
    #ifdef HDF
            fpath = trim(args(2)) // '/output.h5'
            call create_hdf5_file(fpath, NX-2, xs(2:NX-1), m, h, file_id)
    #else
            ufiles = (/ 100, 101, 102, 103/)
    
            fpath = trim(args(2)) // '/xs.dat'
            open(unit=99, file=fpath, status='new', POSITION='append')
            write(99, *) "# X grid"
            call saveOutput(99, NX-2, xs(2:NX-1))
            close(99)
    
            call openOutput(size(ufiles), ufiles, trim(args(2)))
    #endif
    
        end if
    
        ! Produce initial data
        call initial_data(NX, xs, h, m, r0, a0, id, u_p)
        
        counter = 0
        ! Time evolution
        t = 0_RK
        done = .false.
        BHpresent = .false.
    !$OMP PARALLEL DEFAULT(SHARED)
        do while (.not.done)
    
            ! Determine dt from Courant condition
            call CLF(NX, u_p(1:NX), xs, h, vel, dt)
    
    !$OMP SINGLE
            if ((t + dt).gt.T_final) then
                dt = T_final - t
                t = T_final
                done = .true.
            else
                t = t + dt
            end if      
    !$OMP END SINGLE NOWAIT
            
            ! Perform time step
            call TVD_RK(NX, u_p,  u1, u2, u12, u32, f_prime, xs, h, dt, nghost, u)
            call CompExpansion(NX, u(1:NX), u(NX+1:2*NX), xs, theta)
    
    !$OMP SECTIONS
        !$OMP SECTION
            if ( (.not.BHpresent) .and. any(theta.lt.0) ) then
                BHpresent = .true.
                Tformation = t
            else if ( (BHpresent) .and. all(theta.gt.0) ) then
                BHpresent = .false.
                Texplosion = t
            end if
        !$OMP SECTION
            saveO = (N_save.gt.0) .and. ((mod(counter, N_save).eq.0) .or. done)
        !$OMP SECTION
            printO = (N_output.gt.0) .and. ((mod(counter, N_output).eq.0) .or. done)
    !$OMP END SECTIONS
    
            if ( saveO .or. printO ) then
                call compRho(NX, h, dt, u(1:NX), u_p(1:NX), u(NX+1:2*NX), xs, e_der, rho)
            end if
    
            ! TERMINAL OUTPUT
            if ( printO ) then
                ! COMPUTE MASS
    !$OMP MASTER
                ssum = 0.5_RK * (rho(2) * xs(2)**2 + rho(NX-1) * xs(NX-1)**2)
    !$OMP END MASTER
    !$OMP DO PRIVATE(i) REDUCTION(+:ssum)
                do i = 3, NX-2
                    ssum = ssum + rho(i) * xs(i)**2
                end do
    !$OMP END DO
    !$OMP MASTER
                ssum = 4 * PI * ssum * h
                write(*, "(2x, F8.3, 3x, A2, 3x, I9, 4x, A2, 7x, l1, 8x, A2, 3x, E11.3, 1X, I1)") &
                        t, "||",  counter,  "||", BHpresent, "||",  ssum - m
    !$OMP END MASTER
            end if
    
            if ( saveO ) then
    
    #ifdef HDF
    !$OMP SINGLE
                call write_arrays_to_hdf5(file_id, NX-2, u(2:NX-1), u(NX+2:2*NX-1), rho(2:NX-1), t, dt, counter)
    !$OMP END SINGLE NOWAIT
    #else
    !$OMP SECTIONS
        !$OMP SECTION
                call saveOutput(ufiles(1), NX-2, u(2:NX-1))
        !$OMP SECTION
                call saveOutput(ufiles(2), NX-2, u(NX+2:2*NX-1))
        !$OMP SECTION
                call saveOutput(ufiles(3), NX-2, rho(2:NX-1))
        !$OMP SECTION
                ttt = (/t, dt, logic2dbl(BHpresent) /)
                call saveOutput(ufiles(4), size(ttt), ttt)
    !$OMP END SECTIONS NOWAIT
    #endif
    
            end if
    
    !$OMP MASTER
            counter = counter + 1
    !$OMP END MASTER
    
            ! Update previous step
            if ( .not.done ) then
    !$OMP DO SIMD PRIVATE(i)
                do i = 1, 2*NX
                    u_p(i) = u(i)
                end do
    !$OMP END DO SIMD 
            end if
    
        
        end do
    !$OMP END PARALLEL
    
        deallocate(xs, e_der, vel, u, u_p, u1, u2, u12, u32, f_prime, rho, theta)
        
        if (N_save .gt. 0) then
    
    #ifdef HDF
            call close_hdf5_file(file_id)
    #else
        call closeOutput(size(ufiles), ufiles)
    
        allocate(vvv(4))
        fpath = trim(args(2)) // '/details.dat'
        vvv = (/ m, h, Tformation, Texplosion /)
        open(unit=104, file=fpath, status='new', POSITION='append')
        write(104, *) "# Mass dx BH_formation_time BH_explosion_time"
        call saveOutput(104, size(vvv), vvv)
        close(104)
        deallocate(vvv)
    
    #endif
    
        end if
    
        write(*, "(A65)") "-----------------------------------------------------------------"
        write(*, "(A48)") "                  All done!                     "
        write(*, "(A65)") "-----------------------------------------------------------------"
    
        call SYSTEM_CLOCK(iTimes2)
        systemtime = real(iTimes2-iTimes1)/real(rate)
        write(*, "(A21,1x, F10.2, A9)") "Total system runtime:", systemtime, " seconds."
    
    end program LQGeq