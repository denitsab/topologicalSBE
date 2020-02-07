!========================================================================!            
!============================ Main program  =============================!
!========================================================================!

        ! Main program 
        program main
          use omp_lib
          use math_constants
          use phys_constants
          use math 
          use iostream
          use grid
          use solid
          use bloch 
          use field          
          use info 
          use ode

          implicit none
          ! Input parameters supplied by the user in the lines below
          real*8                                  :: I0_TWcm2, lambda_nm  ! peak intensity [ TW/cm2 ] and carrier wavelength [ nm ]
          integer                                 :: ncyc                 ! number of cycles defining the FWHM of the electric field envelope 
          real*8                                  :: aANG                 ! lattice constant a (rhombohedral lattice) in Angstroms
          real*8                                  :: cANG                 ! lattice constant c (rhombohedral lattice) in Angstroms
          real*8                                  :: T2frac               ! dephasing time in terms of fraction of the laser period 
          real*8                                  :: f0_                  ! regularization parameter f0 supplied in (Bohr)^(-1) (absolute value!)

          ! Internal variables 
          integer                                 :: i ,j                 ! indices 
          real*8                                  :: t, Time1, Time2      ! time variable, start and end time of the calculation (CPU time)
          complex*16, allocatable, dimension(:)   :: Y0                   ! density matrix vector with initial values (2*nkpts)
          complex*16, allocatable, dimension(:,:) :: Y                    ! density matrix vector at output (2*nkpts, ntpts)
          integer                                 :: nvar                 ! number of coupled ODEs (2*nkpts)
          real*8                                  :: tstart, tend         ! not used 
          complex*16                              :: nc_t, pi_t           ! population and coherence as a function of k at a given t 
          integer                                 :: iprint, iunit
          character(len=256)                      :: fid 
          logical                                 :: wfprint

          ! The variables below are required for the ZVODE routine. See the documentation of zvode.f for details
          integer :: neq, itol, itask, istate, lzw, mf, liw, lrw 
          complex*16, allocatable, dimension(:) :: zwork  
          real*8, allocatable, dimension(:)     :: rwork 
          real*8, dimension(1)                  :: rtol, atol
          integer, allocatable, dimension(:)    :: iwork 
          real*4, dimension(1)                  :: rpar 
          integer, dimension(1)                 :: ipar 

          ! ZVODE control variables to be set/modified by the user
          integer                               :: iopt                   ! Indicates optional input in the ZVODE. 0: use defaults, 1: read in  optional input 
          logical                               :: testing                ! TRUE: do not perform integration, just prints out energies/velocities/dipoles/Berry curvature of VB/CB. FALSE: run full program, i.e. integrate SBE
          real*8                                :: Tzvode1, Tzvode2       ! variables for tracking CPU time specnt in ZVODE on each time step 
          real*8                                :: hmin, hmax, h0         ! parameters for integration time step control in ZVODE: hmin: minimum step size, hmax: maximum step size, h0: trial suggestion for the step size for the first run 
          integer                               :: maxiter                ! maximum number of iterations per step 

          ! OpenMP: number of CPUs, master id, number of threads 
          
          integer :: ncpu, myid, nthreads 

          

          ! Test run vs. full calculation 
          testing = .FALSE.  

          ! Set the grid parameters, the laser field parameters and the regularization constant 
          ! currently, these are hard-coded, but will be supplied in an external file in a future version 
          nkx = 150                     ! number of points in kx direction 
          nky = 150                      ! number of points in ky direction 
          nkpts = nkx * nky              ! total number of k-grid points 
          ntpts = 5000                   ! number of time steps 
          ncyc = 14                       ! number of cycles (FWHM electric field envelope)
          lambda_nm = 7500.d0            ! center wavelength 
          sigp = 0.005                   ! local regularization factor: Gaussian width // eps = (f0*A0)**2 * exp(- (kx**2+ky**2) / sigma**2 )
          ! Note: sigp denotes the percent of the BZ to be covered. The actual width is given by: sigma = sigp * 4*pi/alat/3.0
          f0_  = 0.0001!0.01                 ! regularization parameter entering eps = (f0 * A0)**2. Units: (Bohr)^(-1). Absolute value! 
          bzcut = 1.0
                    
          ! Lattice constants         
          aANG = 4.138                   ! Angstroms; lattice constants of the layered triangular lattices 
          cANG = 28.64                   ! Angstroms; longest diagonal  of the rhombhohedral unit cell 

          I0_TWcm2 = 0.0075d0              ! peak intensity in TW/cm2 
          T2frac   = 0.05d0             ! dephasing time T2, in fraction of a laser period, i.e. T2 = T2frac * Tper

          alphaQWP = 0.d0 * pi / 180.d0! pi/4.d0

          f_kmask = 0.2d0
          f_dmask = 0.6d0
          fabs_vmask = 1.d0*0.3d0
          fbrd_vmask = 1.d0*0.5d0
          q_vmask = 2.d0
          q_kmask = 1.d0/8.d0

          wfprint = .TRUE.


          ! Set number of threads 
          ncpu = 8

          ! Set ZVODE integration step controls 
          h0   = 0.05d0;                 ! Initial guess for the step size 
          hmax = 10.3d0;                  ! Maximum step size 
          hmin = 0.0005d0                ! Minimum step size 

          theta = 0.d0 * pi/180.d0

          ! OMP 

          call omp_set_num_threads(ncpu)

!$omp parallel private(nthreads, myid)
          myid = omp_get_thread_num()
          nthreads = omp_get_num_threads()
          if(myid.eq.0) then
            write(*,*) "Message from id:", myid
            write(*,*) "Starting SFA with ",nthreads, " threads"
          endif
!$omp end parallel

          ! Set up field 
          call field_init( I0_TWcm2, lambda_nm, ncyc )
          ! Set up crystal and momentum space 
          call set_crystal(aANG, cANG, T2frac, f0_)

          ! Print information 
          call print_fields()
          call print_bands()


          ! Allocate space for solution vector
          if(.not.allocated(Y0)) allocate( Y0( 2*nkpts ) )
          if(.not.allocated(Y )) allocate( Y ( 2*nkpts, ntpts ) )

          nvar = 2*nkpts      ! total number of equations 

          ! Adjust the ZVODE parameters (see documentation)
          neq = nvar 
          mf = 10 
          itask = 1
          istate=1
          iopt = 1
          lzw = 15*neq
          lrw = 20+neq 
          liw = 30 
          maxiter= 5000

          itol = 1
          rtol = 1.d-9
          atol = 1.d-9

          ipar = 0
          rpar = 0.d0 

          ! Allocate arrays for ZVODE
  
          if(.not.allocated(zwork)) allocate(zwork(lzw))
          if(.not.allocated(rwork)) allocate(rwork(lrw))
          if(.not.allocated(iwork)) allocate(iwork(liw))


          ! Supply the optional input to ZVODE, see documentation for details 
          if(iopt .ne. 0) then 
            rwork(5:10) = 0.d0
            iwork(5:10) = 0
            rwork(5) = h0
            rwork(6) = hmax
            rwork(7) = hmin
            iwork(6) = maxiter
          else
            h0 = 0.d0
            hmin = 0.d0 
            hmax = 0.0d0 
          endif 

          ! Initial condition 
          Y0(1:neq) = czero

          tstart = tgrid(1)
          tend   = tgrid(ntpts)

          !call CPU_TIME(Time1)
          Time1 = omp_get_wtime()      

        
          if(.not. testing) then 

          open(21, FILE=Opath(1:Olast)//'cpopulation'//trim(calcID)//'.dat', status="replace")
          open(23, FILE=Opath(1:Olast)//'cvcoherence'//trim(calcID)//'.dat', status="replace")

          iunit = 0
          iprint = 100

          do i = 1, ntpts
            t = tgrid(i)

            Tzvode1 = omp_get_wtime()

            call zvode(derivs, nvar, y0, t, t+dt, itol, rtol, atol, itask, &
              istate, iopt, zwork, lzw, rwork, lrw, iwork, liw, jex, mf, rpar, ipar)

            Tzvode2 = omp_get_wtime()

            write(*,*) Tzvode2-Tzvode1

            y0( 1         :   nkpts) = y0(        1 :   nkpts) * kmask
            y0( nkpts + 1 : 2*nkpts) = y0(nkpts + 1 : 2*nkpts) * kmask

            y(:,i) = y0        
            nc_t = sum( y0(       1 :   nkpts ) )*dkx*dky 
            pi_t = sum( y0( nkpts+1 : 2*nkpts ) )*dkx*dky
            write(21, 2002) t, dble(nc_t), aimag(nc_t)
            write(23, 2002) t, dble(pi_t), aimag(pi_t)
            write(88,*) istate

            if( mod(i, iprint) .eq.0 ) then
              iunit = iunit + 1
              call print_snapshot(iunit, y0(1:nkpts), y0(nkpts+1:2*nkpts), t)
            endif
          
          enddo 

          close(21)
          close(23)

          ! Deallocate working space for ZVODE
          if(allocated(zwork)) deallocate(zwork )
          if(allocated(rwork)) deallocate(rwork )
          if(allocated(iwork)) deallocate(iwork )

          !call CPU_TIME(Time2)
          Time2 = omp_get_wtime()
          Time2 = Time2 - Time1

          write(*,*) 'Finished ODE solution step ...'
          write(*,*) 'Time spent in integration routine:', Time2

          write(*,*) "Dumping the wavefunction to files"         

          if (wfprint) then 
            do i = 1, ntpts
              write(fid,'(I10)') i
              open(27, FILE=Opath(1:Olast)//'wfdump'//trim(calcID)//'_'//trim(fid)//'.dat', status="replace", form="unformatted")            
              do j = 1, nkpts 
                write(27) tgrid(i), kxy_grid(j,1), kxy_grid(j,2), dble(y(j,i)), aimag(y(j,i)), &
                dble(y(nkpts+j,i)), aimag(y(nkpts+j,i))
              enddo
              close(27)
            enddo
          endif

          write(*,*) 'Calculating currents ...'


          ! Calculate currents 

          call currents(Y, tgrid)
          
          write(*,*) 'Finished current calculation step ...'


          !call CPU_TIME(Time2)
          Time2 = omp_get_wtime()
          Time2 = Time2 - Time1

          ! Print interband and intraband currents 

          open(18, FILE=Opath(1:Olast)//'interband_dipole'//trim(calcID)//'.dat', status="replace")
          do i = 1, ntpts
            write(18, 2004) tgrid(i), dble(Pt(i,1)), aimag(Pt(i,1)), dble(Pt(i,2)), aimag(Pt(i,2))
          enddo
          close(18)

          open(19, FILE=Opath(1:Olast)//'intraband_dipole'//trim(calcID)//'.dat', status="replace")
          do i = 1, ntpts
            write(19, 2004) tgrid(i), dble(Jt(i,1)), aimag(Jt(i,1)), dble(Jt(i,2)), aimag(Jt(i,2))
          enddo
          close(19)

          end if 

          ! Print summary for the calculation 

          call print_info(I0_TWcm2, lambda_nm, Time2, T2frac,hmin,h0,hmax)

          

          ! Deallocate space 
          call field_destroy()
          call crystal_dealloc()
          if(allocated(Y)) deallocate(Y)
          if(allocated(Y0)) deallocate(Y0)

          

 2002     format (1e16.8,2(2x,1e16.8))  
 2003     format (1e16.8,2x,1e16.8)  
 2004     format (1f16.8,4(2x,1e16.8))  


        end program
