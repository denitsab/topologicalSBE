!========================================================================!            
!======================= Print info module ==============================!
!========================================================================!
        module info 
          use math_constants
          use phys_constants
          use iostream
          use bloch
          use field
          use solid

          implicit none
          ! A module that prints a summary of the calculation with all necessary details 

          contains

          subroutine print_info(I0_TWcm2, lambda_nm, time_s, T2frac, hmin, hmax, h0)
            implicit none
            real*8 :: I0_TWcm2, lambda_nm, time_s, T2frac
            character(8) :: date
            character(10) :: time
            character(5) :: zone 
            integer, dimension(8) :: idvalues
            real*8                :: eV2au, ANG2au
            real*8                :: hmin, hmax, h0 

            eV2au = eV/hartree
            ANG2au = 1.d-10 / aBohr

            call date_and_time(DATE=date,TIME=time,ZONE=zone,VALUES=idvalues)
  
            open(27, file=Opath(1:Olast)//'info'//trim(calcID)//'.dat', status="replace")

            write(27, '(A)') '**********************************************************************************'
            write(27, '(A)') "*****************    Printing info for calculation *******************************"
            write(27, '(A)') calcID
            write(27, '(A, 3(2X, I4, 2X))') "Finished on", idvalues(3), idvalues(2), idvalues(1)
            write(27, '(A, 3(2X, I2, 2X))') "Time now:", idvalues(5), idvalues(6), idvalues(7)
            write(27, '(A, 2X, f10.4, 2x, A)') "Time spent in integration routine: ", time_s, " seconds"
            write(27, '(A)') 'using:'
            write(27, '(i5, 5X, A)') ncycles,  'cycles long pulse'
            write(27, '(f10.4, 5X, A)') T2frac, '*Tper dephasing time '
            write(27, '(i5, 5X, A)') nkx,  'points in kx-dimension'
            write(27, '(i5, 5X, A)') nky,  'points in ky-dimension'
            write(27, '(i5, 5X, A)') ntpts, 'points in t-domain'
            write(27, '(f10.4, 5x, A)') lambda_nm, 'nm center wavelength'
            write(27, '(A, 5x,f10.4, 5x, A)') 'I = ', I0_TWcm2, 'x 10^12 W/cm2'
            write(27, '(A, 5x,f10.4, 5x, A)') 'QWP angle alpha = ', alphaQWP*180.d0/pi, ' deg'
            write(27, '(A, 5x,f10.4, 5x, A)') 'tilt angle theta = ', theta*180.d0/pi, ' deg'
            write(27, '(A)') ' ' 
            write(27, '(A)') ' Bi2Se3 k.p model parameters: '
            write(27, '(A, 2X, f10.4, 2x, A)') "A0 = ", A0 / eV2au / (ANG2au    ), " eV A " 
            write(27, '(A, 2X, f10.4, 2x, A)') "A2 = ", A2 / eV2au / (ANG2au**2 ), " eV A^2 "
            write(27, '(A, 2X, f10.4, 2x, A)') "B0 = ", B0 / eV2au / (ANG2au    ), " eV A "
            write(27, '(A, 2X, f10.4, 2x, A)') "B2 = ", B2 / eV2au / (ANG2au**2 ), " eV A^2 "
            write(27, '(A, 2X, f10.4, 2x, A)') "C0 = ", C0 / eV2au               ,  " eV "
            write(27, '(A, 2X, f10.4, 2x, A)') "C1 = ", C1 / eV2au / (ANG2au**2 ),  " eV A^2 "
            write(27, '(A, 2X, f10.4, 2x, A)') "C2 = ", C2 / eV2au / (ANG2au**2 ), " eV A^2 "
            write(27, '(A, 2X, f10.4, 2x, A)') "M0 = ", M0 / eV2au               , " eV "
            write(27, '(A, 2X, f10.4, 2x, A)') "M1 = ", M1 / eV2au / (ANG2au**2 ),  " eV A^2 "
            write(27, '(A, 2X, f10.4, 2x, A)') "M2 = ", M2 / eV2au / (ANG2au**2 ),  " eV A^2 "
            write(27, '(A, 2X, f10.4, 2x, A)') "R1 = ", R1 / eV2au / (ANG2au**3 ), " eV A^3 "
            write(27, '(A, 2X, f10.4, 2x, A)') "R2 = ", R2 / eV2au / (ANG2au**3 ), " eV A^3 "
            write(27, '(A, 2X, f10.4, 2x, A)') "alpha1 = ", alpha1, " [] "
            write(27, '(A, 2X, f10.4, 2x, A)') "alpha3 = ", alpha3, " [] ."
            write(27, '(A)') ' Numerical control parameters: '
            write(27, '(A, 2X, f10.4, 2x, A)') "f0 = ", f0glob,                      " a0^(-1) "
            write(27, '(A, 2X, f10.4, 2x, A)') "sigma = ", sigp, " * 4 * pi / 3 / a"
            write(27, '(A, 2X, f10.4, 2x, A)') "hmin  = ", hmin, " minimum step in ZVODE"
            write(27, '(A, 2X, f10.4, 2x, A)') "hmax  = ", hmax, " maximum step in ZVODE"
            write(27, '(A, 2X, f10.4, 2x, A)') "h0    = ", h0  , " initial step in ZVODE"
            write(27, '(A, 2X, f10.4, 2x, A)') "Integration over ", bzcut, " * Brillouin zone"
            write(27, '(A, 2X, f10.4, 2x, A)') "f_kmask = ", f_kmask, " [] "
            write(27, '(A, 2X, f10.4, 2x, A)') "q_kmask = ", q_kmask, " [] "
            write(27, '(A, 2X, f10.4, 2x, A)') "fabs_vmask = ", fabs_vmask, " [] "
            write(27, '(A, 2X, f10.4, 2x, A)') "fbrd_vmask = ", fbrd_vmask, " [] "
            write(27, '(A, 2X, f10.4, 2x, A)') "q_vmask = ", q_vmask, " [] "
            write(27, '(A, 2X, f10.4, 2x, A)') "f_dmask = ", f_dmask, " [] ."
            write(27, '(A)') '**********************************************************************************'
            write(27, '(A)') '******************************* Aufwiedersehen ***********************************'
            write(27, '(A)') '**********************************************************************************'
            close(27)

            return
          end subroutine print_info
        end module info
