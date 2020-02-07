!========================================================================!            
!========================= Field Variables  =============================!
!========================================================================!
        module field 
          use math_constants
          use phys_constants   
          use iostream
          use grid       
          use math
        implicit none
          real*8,        dimension(:  ), allocatable :: tgrid   ! time vector 
          real*8,        dimension(:,:), allocatable :: grad_t  ! derivative operator with respect to time (not used at the moment)
          complex*16,     dimension(:,:), allocatable :: Pt, Jt ! arrays for intraband current and polarization
          real*8 :: tau, E0, w0, dt, Tper, phiCEP               ! pulse duration, el. field amplitude, angular frequency, time step, 
                                                                ! pulse period, CEP phase 
          integer :: ncycles                                    ! number of cycles in the electric field, defines the FWHM
          real*8  :: alphaQWP, theta                            ! angle of the QWP plate; direction of the GK crystal axis wrt. the field 
          save tgrid, Pt, Jt, tau, E0, w0, dt, Tper, phiCEP, ncycles, grad_t, &
          alphaQWP, theta 

          contains

          subroutine field_init( I0_TWcm2, lambda_nm, ncyc )
            ! initializes the time grid and makes the conversion of field-related input to atomic units 
            ! Input arguments:
            ! I0_TWcm2      : maximum field intensity in [TW / cm2] == [ 10^12 W / cm2 ]
            ! lambda_nm     : carrer wavelength in [ nm ]
            ! ncycles       : electric field envelope FWHM in terms of number of cylces 

            implicit none
            integer :: ncyc 
            real*8 :: I0_TWcm2, lambda_nm
            real*8 :: tstart, tfin,t, fwhm, s
            character(8) :: date
            character(10) :: time
            character(5) :: zone 
            integer, dimension(8) :: idvalues
            integer :: i 
            ! Conversion to atomic units and initialization of the common variables w0, E0, phiCEP
            w0 = 2.d0*pi*cvac0/(lambda_nm*1.d-9)
            E0 = sqrt( 2.d0 * I0_TWcm2 * 1.d12 / (0.01d0**2) / (cvac0*eps0) )
            w0 = w0*atu
            E0 = E0/el_field
            phiCEP = 0.d0     ! Note: phiCEP is not used at the moment. In a future version, it will be treated as an input argument of the function


            Tper = 2.d0*pi/w0
            ncycles = ncyc
            fwhm = dble(ncycles)*Tper
            ! Note: the electric field function envelope is implemented as:
            ! f(t) = exp( - ( t/tau )**2 )
            ! From the definition of the Gaussian function:
            ! f(t) = exp( - ( t/(2*s**2) ) )
            ! we have the FWHM:
            ! FWHM = 2*sqrt(2 * log(2))*s 
            ! Then, it follows that: tau = sqrt(2) * s 
            ! The lines below perform the above transformation
            s = fwhm/(2.d0*sqrt(2.d0*log(2.0)))
            tau = sqrt(2.0)*s 
                        
            ! Initial and final point of hte temporal grid. The density matrix is propagated for a couple of periods after the pulse has subsided
            !tstart = -5.d0 *  tau
            !tfin = 5.d0 * tau 
            tstart = -(dble(ncycles+1))/2.d0*Tper 
            tfin = (dble(ncycles+4))/2.d0*Tper 

            ! Distance of the temporal grid 
            dt = (tfin-tstart)/dble(ntpts-1)

            ! Array allocation
            if(.not.allocated(tgrid)) allocate(tgrid(ntpts))
            if(.not.allocated(Pt))   allocate(Pt(ntpts,2))
            if(.not.allocated(Jt))   allocate(Jt(ntpts,2))
            if(.not.allocated(grad_t)) allocate(grad_t(ntpts, ntpts))
            
            ! Generate a unique label from time and date that is appended to the output files (populations, coherences, dipoles)
            ! Helps to avoid accidental deletion 

            call date_and_time(DATE=date,TIME=time,ZONE=zone,VALUES=idvalues)
            calcID = "_"//trim(date)//"_"//trim(time(1:6))

            ! Initialization
            tgrid = dzero
            Pt    = czero
            Jt    = czero

            do i = 1,ntpts
              t = tstart + dble(i-1)*dt 
              tgrid(i) = t               
            enddo

            ! Generate the differential operator 
            grad_t = DMatrixOP(ntpts, dt)
            
          end subroutine field_init

          ! The next two routines are not used in the current code
          ! They perform a numerical integration of an array by means of cummulative summation of the elements 
          ! Useful for future, in case few-cycle pulses implemented 

          function cumsum1d(x, n) result(y)
            ! Cummulative sum of a 1D array
            ! x: input array 
            ! n: array size 
            integer :: n, i 
            real*8, dimension(n) :: x
            real*8, dimension(n) :: y            
            y(1)=x(1)
            do i=2,n
              y(i) = y(i-1) + x(i)
            enddo
          
            return
          end function cumsum1d

          function cumsum2d(x, n) result(y)
            ! Cummulative sum of a 2D array with dimensions (n,2)
            ! i.e. an array holding the electric field in 2D 

            integer :: n, i 
            real*8, dimension(n,2) :: x
            real*8, dimension(n,2) :: y            
            y(1,:)=x(1,:)
            do i=2,n
              y(i,1) = y(i-1,1) + x(i,1)
              y(i,2) = y(i-1,2) + x(i,2)
            enddo
            return
          end function cumsum2d

          real*8 function fenv(t)
            implicit none 
            real*8 ::t
            ! envelope function: cos2
            if(abs(t).lt.dble(ncycles)*Tper/2.d0) then
              fenv = (cos(pi*t/ ( dble(ncycles) * Tper ) ))**2
            else
             fenv = dzero
            endif
            return
          end function fenv


          real*8 function g_envelope(t)
          ! electric field envelopw function 
            implicit none 
            real*8 :: t 
            g_envelope = exp( - (t/tau)**2 ) 
          end function g_envelope

          real*8 function Efun(t) 
          ! electric field at time t (1D)
            implicit none
            real*8 :: t 
            Efun = E0*g_envelope(t) * cos(w0*t + phiCEP)           
            return
          end function Efun

          function Efun2(t) result(Exy)
          ! electric field in 2D (x,y)-plane at time t 
          ! for future use 
            implicit none
            real*8 :: t 
            real*8, dimension(2) :: Exy
            real*8 :: Etx, Ety 

            !Exy(1) = E0*g_envelope(t) * (cos(w0*t + phiCEP) - cos(2.d0*alphaQWP)*sin(w0*t+ phiCEP))
            !Exy(2) = E0*g_envelope(t) * (                     sin(2.d0*alphaQWP)*sin(w0*t+ phiCEP))
            !Etx = E0*fenv(t) * (cos(w0*t + phiCEP) - cos(2.d0*alphaQWP)*sin(w0*t+ phiCEP))
            !Ety = E0*fenv(t) * (                     sin(2.d0*alphaQWP)*sin(w0*t+ phiCEP))

            if(abs(t).lt.dble(ncycles)*Tper/2.d0) then
              Etx = - (( E0*pi*sin( (2.d0*pi*t)/(ncycles*Tper) ) * ( cos(2.d0*alphaQWP)*cos(t*w0) + sin(t*w0) ) ) / &
              ( ncycles*Tper*w0 ) + E0*cos( (pi*t)/(ncycles*Tper) )**2 * ( -cos(t*w0) + cos(2.d0*alphaQWP)*sin(t*w0) ))

              Ety = ( ( E0*sin(2.d0*alphaQWP) * ( pi*cos(t*w0)*sin( (2.d0*pi*t)/(ncycles*Tper) ) + &
                ncycles*Tper*w0*cos( (Pi*t)/(ncycles*Tper) )**2 * sin(t*w0) ) ) / ( ncycles*Tper*w0 ) )
            else
             Etx = dzero
             Ety = dzero
            endif

            Exy(1) = Etx * cos(theta) - Ety * sin(theta)
            Exy(2) = Etx * sin(theta) + Ety * cos(theta)


            return
          end function Efun2

          real*8 function Afun(t) 
          ! multi-cycle vector potential (approximate) at time t, 1D 
            implicit none
            real*8 :: t, A0
            A0 = E0/w0
            Afun = - A0 * g_envelope(t) * sin(w0*t + phiCEP)           
            return
          end function Afun

          function Afun2(t) result(Axy)
          ! multi-cycle vector potential (approximate) at time t, 2D
            implicit none
            real*8 :: t, A0 
            real*8 :: Atx, Aty 
            real*8, dimension(2) :: Axy
            A0 = E0/w0

            !Axy(1) = - A0*g_envelope(t) * (sin(w0*t + phiCEP) + cos(2.d0*alphaQWP)*cos(w0*t+ phiCEP))
            !Axy(2) = - A0*g_envelope(t) * (                   - sin(2.d0*alphaQWP)*cos(w0*t+ phiCEP))
            Atx = - A0*fenv(t) * (sin(w0*t + phiCEP) + cos(2.d0*alphaQWP)*cos(w0*t+ phiCEP))
            Aty = - A0*fenv(t) * (                   - sin(2.d0*alphaQWP)*cos(w0*t+ phiCEP))


            Axy(1) = Atx * cos(theta) - Aty * sin(theta)
            Axy(2) = Atx * sin(theta) + Aty * cos(theta)
            
            return
          end function Afun2

          subroutine field_destroy()
          ! deallcoate arrays to release memory 
          implicit none
            if(allocated(tgrid)) deallocate(tgrid)
            if(allocated(Jt)) deallocate(Jt)
            if(allocated(Pt)) deallocate(Pt)
            if(allocated(grad_t)) deallocate(grad_t)
            return 
          end subroutine field_destroy

          subroutine print_fields()
          ! print At and Et to files for comparison purposes 
            implicit none
            real*8 :: t 
            real*8, dimension(2) :: At, Et 
            integer :: i 

            open(171,file=Opath(1:Olast)//'efield.dat', status = "replace")
            open(173,file=Opath(1:Olast)//'afield.dat', status = "replace")

            do i = 1, ntpts 
              t = tgrid(i)
              At = Afun2(t)
              Et = Efun2(t) 
            !  write(171, 1001) t, Efun(t)
            !  write(173, 1001) t, Afun(t)   
              write(171, 1003) t, Et(1), Et(2)
              write(173, 1003) t, At(1), At(2)
            enddo

            close(171)
            close(173)

 1001     format (1e16.8,   2x,1e16.8)             
 1003     format (1e16.8, 2(2x,1e16.8))  
          end subroutine print_fields

        end module field
