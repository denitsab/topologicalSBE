!========================================================================!            
!============================= Band structure ===========================!
!========================================================================!
        module solid
          use math_constants
          use phys_constants
          use grid
          use legendre
          use field
          implicit none 

          real*8 :: dkx, dky                                        ! Momentum grid spacings in kx/ky directions
          real*8 :: A0, A2, B0, B2, C0, C1, C2                      ! Parameters for the Bi2Se3 Hamiltonian, from Liu et al. (2010)
          real*8 :: M0, M1, M2, R1, R2                              ! Parameters for the Bi2Se3 Hamiltonian, from Liu et al. (2010)
          real*8 :: f0DM, f0BC                                      ! Scaling factors for regularization parameter for calculating dipole moments / Berry connection                                       
          real*8 :: alat                                            ! Lattice constant a of the rhombohedral lattice 
          real*8 :: Tdephas_2                                       ! Dephaising time T2, in a.u. 
          real*8 :: f_kmask, q_kmask
          real*8, allocatable, dimension(:,:) :: kxy_grid           ! 2D array holding all nkx * nky grid points (nkx * nky, 2)
          real*8, allocatable, dimension(:)   :: kx1d, ky1d         ! 1D arrays of k-points in x and y dimensions
          real*8, allocatable, dimension(:)   :: kmask
          real*8                                    :: kzpoint      ! Value of the longitudinal (wrt propagation direction) momentum to be used. HHG occurs in the plane perpendicular to the kz-axis
          integer                                   :: icomp        ! Selected spinor from the pair
          real*8                                    :: sigp         ! Width of the Gaussian in the local regularization factor. Denotes the fraction of the BZ to be used as a widt, i.e. 
                                                                    ! the actual width sigma is given by: sigma = sigp * 4*pi/3/alat. 4 pi/3/alat is the extent of the BZ in the Gamma-K direction 
          real*8                                    :: bzcut        ! factor controlling the reduction of the BZ integration, i.e. integrate only over bzcut * BZ. Use to speed up 2D calculations, but with caution!
        

          save Tdephas_2, kxy_grid, kx1d, ky1d, kmask, f_kmask, q_kmask, &
          dkx,dky,  f0DM, f0BC, kzpoint, alat, bzcut, &
           A0, A2, B0, B2, C0, C1, C2, M0, M1, M2, R1, R2, sigp 

          contains

          subroutine set_crystal(aANG, cANG, T2frac, facDM, facBC, kz, ic)
          ! Convert input lattice constant to atomic units 
          ! Set up the NN and NNN vectors and matrices 
          ! Allocate momentum grid 
          ! Input: 
          ! aANG    : lattice constant in Angstroms, from main program
          ! T2frac  : dephasing time in fractions of the laser period, i.e. T2  = T2frac * T0 
          ! facDM   : factor to be used in setting the final scaling   factor for the DM reg. parameter 
          ! facBC   : factor to be used in setting the final scaling   factor for the BC reg. parameter. Please set facDM=facBC
          ! kz      : value for the kzpoint

            implicit none
            real*8 :: T2frac, aANG, cANG, kstart, kend, kx, ky 
            real*8 :: kxstart, kystart, kxend, kyend
            real*8 :: facDM, facBC, kz
            integer :: i, j 
            integer ic 
            real*8 :: eV2au, ANG2au

            real*8 :: kbuff, kabs, kpar

            icomp = ic 


            ! Convert from SI to atomic units 
            eV2au  = eV/hartree 
            ANG2au = 1.d-10 / aBohr

            A0 = 3.33d0    * eV2au * ANG2au    ;
            A2 = 0.d0      * eV2au * ANG2au**2 ;
            B0 = 2.26d0    * eV2au * ANG2au    ;
            B2 = 0.d0      * eV2au * ANG2au**2 ;
            C0 = -0.0083d0 * eV2au             ; 
            C1 = 5.74      * eV2au * ANG2au**2 ; 
            C2 = 30.4      * eV2au * ANG2au**2 ;
            M0 = -0.28d0   * eV2au             ;
            M1 = 6.86d0    * eV2au * ANG2au**2 ; 
            M2 = 44.5d0    * eV2au * ANG2au**2 ; 
            R1 = 50.6d0    * eV2au * ANG2au**3 ;
            R2 = -113.3d0  * eV2au * ANG2au**3 ;
            
            Tdephas_2 = Tper*T2frac        
            alat = aANG*1.d-10/aBohr            

        
            nkpts = nkx * nky                     ! Total number of k-points in 2D 
            nmax = 2*nkpts                        
          

            ! Allocate arrays
            if(.not.allocated(kx1d)) allocate(kx1d(nkx))
            if(.not.allocated(ky1d)) allocate(ky1d(nky))
            if(.not.allocated(kmask)) allocate(kmask(nkpts))
            if(.not.allocated(kxy_grid)) allocate(kxy_grid(nkpts,2))
        

            kystart  = -4.d0*pi/alat/3.d0*bzcut
            kyend    =  4.d0*pi/alat/3.d0*bzcut
            kxstart  = -2.d0*pi/alat/sqrt(3.d0)*bzcut
            kxend    =  2.d0*pi/alat/sqrt(3.d0)*bzcut

            dkx = (kxend-kxstart)/dble(nkx-1)

            ! Initialize 1D grids 
            do i = 1, nkx 
               kx1d(i) = dble(i-1)*dkx + kxstart 
            enddo

            ! If nky is set to 1, the program enters the 1D mode, i.e. ky=0

            if(nky .gt. 1) then 
              dky = (kyend-kystart)/dble(nky-1)                      
              do i = 1, nky 
                ky1d(i) = dble(i-1)*dky + kystart 
              enddo
            else
              ky1d(1) = 0.d0
              dky = 1.d0 
            endif

            kbuff = f_kmask * (kxend - 0.d0)
            kabs  = kxend - kbuff 



            ! Construct the kxy grid. This is a nkpts x 2 grid of all possible momenta. 
            ! Mapping: ( kx(i), ky(j) ) --> kxy ( (i-1)*nky + j, 1:2 )
            ! Note that population, coherences etc are stored in this format 
            do i = 1, nkx 
              kx = kx1d(i)
              do j = 1, nky
                ky = ky1d(j)
                kxy_grid( (i-1)*nky + j, 1 ) = kx
                kxy_grid( (i-1)*nky + j, 2 ) = ky

                kpar = sqrt(kx**2 + ky**2)

                if( kpar .lt. kabs ) then 
                  kmask( (i-1)*nky + j ) = 1.d0 
                else if( kpar .lt. kxend ) then
                  kmask( (i-1)*nky + j ) = abs(cos(abs(kabs-kpar) * pi / (2.d0*kbuff)))**q_kmask
                else
                  kmask( (i-1)*nky + j ) = 0.d0
                endif

              enddo
            enddo            

            kzpoint = kz 
            ! There are two regularization factors defined here, however, the distinction is not yet implemented in the other routines. Set both of them equa, do not mix
            f0DM = facDM 
            f0BC = facBC 

            return
          end subroutine set_crystal

          real*8 function Ebulk(kxyz, m)
          ! Eigenenergy of the bulk VB/CB, 2xdegenerate
          ! kxyz :: 3-vector (kx, ky, kz)
          ! m: = -1 --> VB; = +1 --> CB 

            implicit none
            real*8, dimension(3) ::  kxyz 
            integer              ::  m 
            real*8               ::  kx, ky, kz, kpar
            real*8               ::  sgn 

            sgn  = dble(m)
            kx   = kxyz(1)
            ky   = kxyz(2)
            kz   = kxyz(3)
            kpar = sqrt( kx**2 + ky**2 )

            Ebulk = C0 + C2*kpar**2 + C1*kz**2 + Sqrt(A0**2*kpar**2 + 2.d0*A0*A2*kpar**4 + &
              A2**2*kpar**6 + M0**2 + kz**2*((B0 + B2*kz**2)**2 + 2.d0*M0*M1 + kz**2*M1**2) + &
              2.d0*kx**2*M0*M2 + 2.d0*ky**2*M0*M2 + 2.d0*kx**2*kz**2*M1*M2 + 2.d0*ky**2*kz**2*M1*M2 + &
              kx**4.d0*M2**2 + 2.d0*kx**2*ky**2*M2**2 + ky**4*M2**2 + kx**6*R1**2 - 6.d0*kx**4*ky**2*R1**2 + &
              9.d0*kx**2*ky**4*R1**2 - 2.d0*ky*(-3.d0*kx**2 + ky**2)*kz*(B0 + B2*kz**2)*R2 + &
              ky**2*(-3.d0*kx**2 + ky**2)**2*R2**2)*sgn
            return

          end function Ebulk

          function dEbulk_dk(kxyz, m) result(dE_dk)
          ! Momentum derivative of band eigenenergy; group velocity
            implicit none
            real*8, dimension(3) ::   kxyz, dE_dk
            integer              ::   m 
            real*8               ::   sgn, kx, ky, kz, kpar 
            real*8               ::   dE_dkx, dE_dky, dE_dkz
            kx = kxyz(1)
            ky = kxyz(2)
            kz = kxyz(3)

            kpar = sqrt( kx**2 + ky**2 )

            sgn  = dble(m)
            kpar = sqrt( kx**2 + ky**2 )

            dE_dkx = 2.d0*C2*kx + (kx*((A0 + A2*kpar**2)*(A0 + 3.d0*A2*kpar**2) + 2.d0*M0*M2 + &
              2.d0*kz**2*M1*M2 + 2.d0*kx**2*M2**2 + 2.d0*ky**2*M2**2 + 3.d0*kx**4*R1**2 - &
              12.d0*kx**2*ky**2*R1**2 + 9.d0*ky**4*R1**2 + 6.d0*ky*kz*(B0 + B2*kz**2)*R2 - &
              6.d0*ky**2*(-3.d0*kx**2 + ky**2)*R2**2)*sgn)/Sqrt(A0**2*kpar**2 + 2.d0*A0*A2*kpar**4 +&
              A2**2*kpar**6 + M0**2 + kz**2*((B0 + B2*kz**2)**2 + 2.d0*M0*M1 + kz**2*M1**2) + 2.d0*kx**2*M0*M2 + &
              2.d0*ky**2*M0*M2 + 2.d0*kx**2*kz**2*M1*M2 + 2.d0*ky**2*kz**2*M1*M2 + kx**4*M2**2 + &
              2.d0*kx**2*ky**2*M2**2 + ky**4*M2**2 + kx**6*R1**2 - 6.d0*kx**4*ky**2*R1**2 + &
              9.d0*kx**2*ky**4*R1**2 - 2.d0*ky*(-3.d0*kx**2 + ky**2)*kz*(B0 + B2*kz**2)*R2 + ky**2*(-3.d0*kx**2 + ky**2)**2*R2**2)


            dE_dky = 2.d0*C2*ky + ((ky*((A0 + A2*kpar**2)*(A0 + 3.d0*A2*kpar**2) + 2.d0*M2*(M0 + kz**2*M1 + &
              (kx**2 + ky**2)*M2) - 6.d0*kx**2*(kx**2 - 3.d0*ky**2)*R1**2) + &
              3.d0*(kx - ky)*(kx + ky)*kz*(B0 + B2*kz**2)*R2 + &
              3.d0*ky*(3.d0*kx**4 - 4.d0*kx**2*ky**2 + ky**4)*R2**2)*sgn)/&
              Sqrt(kpar**2*(A0 + A2*kpar**2)**2 + B0**2*kz**2 + 2.d0*B0*B2*kz**4 + &
              B2**2*kz**6 + M0**2 + 2.d0*kz**2*M0*M1 + kz**4*M1**2 + 2.d0*kx**2*M0*M2 +&
              2.d0*ky**2*M0*M2 + 2.d0*kx**2*kz**2*M1*M2 + 2.d0*ky**2*kz**2*M1*M2 + kx**4*M2**2 + &
              2.d0*kx**2*ky**2*M2**2 + ky**4*M2**2 + kx**6*R1**2 - 6.d0*kx**4*ky**2*R1**2 + &
              9.d0*kx**2*ky**4*R1**2 - 2.d0*ky*(-3.d0*kx**2 + ky**2)*kz*(B0 + B2*kz**2)*R2 + &
              ky**2*(-3.d0*kx**2 + ky**2)**2*R2**2)

            dE_dkz = 2.d0*C1*kz + ((kz*(B0**2 + 4.d0*B0*B2*kz**2 + 3.d0*B2**2*kz**4 + &
              2.d0*M1*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2)) - &
              ky*(-3.d0*kx**2 + ky**2)*(B0 + 3.d0*B2*kz**2)*R2)*sgn)/&
              Sqrt(kpar**2*(A0 + A2*kpar**2)**2 + B0**2*kz**2 + 2.d0*B0*B2*kz**4 + &
              B2**2*kz**6 + M0**2 + 2.d0*kz**2*M0*M1 + kz**4*M1**2 + 2.d0*kx**2*M0*M2 + &
              2.d0*ky**2*M0*M2 + 2.d0*kx**2*kz**2*M1*M2 + 2.d0*ky**2*kz**2*M1*M2 + &
              kx**4*M2**2 + 2.d0*kx**2*ky**2*M2**2 + ky**4*M2**2 + kx**6*R1**2 - &
              6.d0*kx**4*ky**2*R1**2 + 9.d0*kx**2*ky**4*R1**2 - 2.d0*ky*(-3.d0*kx**2 + ky**2)*kz*(B0 + B2*kz**2)*R2 + &
              ky**2*(-3.d0*kx**2 + ky**2)**2*R2**2)

            dE_dk = (/ dE_dkx, dE_dky, dE_dkz /)
            return

          end function dEbulk_dk

          real*8 function psiNorm(kxyz, m, f0, sigma)
          ! Norm of the bulk spinor, each component is normalized independently 
            implicit none
            real*8, dimension(3) :: kxyz 
            real*8               :: kx, ky, kz, kpar
            integer              :: m
            real*8               :: sgn, NNorm
            real*8               :: f0, eps 
            real*8, optional     :: sigma ! local regularization, optional 
            real*8               :: s 

            kx = kxyz(1)
            ky = kxyz(2)
            kz = kxyz(3)

            kpar = sqrt( kx**2 + ky**2 )

            s = 0.d0 
            if(present(sigma)) s = 1.d0/( ( sigma*4.d0*pi/alat/3.d0 )**2 )


            eps = (f0 * A0)**2 * exp(-kpar**2 * s)

            kpar = sqrt( kx**2 + ky**2 )
            sgn = dble(m)

            NNorm = A0**2*kpar**2 + B0**2*kz**2 + kx**6*R1**2 - 6.d0*kx**4*ky**2*R1**2 + &
             9.d0*kx**2*ky**4*R1**2 + 6.d0*B0*kx**2*ky*kz*R2 - 2.d0*B0*ky**3*kz*R2 + 9.d0*kx**4*ky**2*R2**2 - &
             6.d0*kx**2*ky**4*R2**2 + ky**6*R2**2 + (M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + &
             Sqrt(A0**2*kpar**2 + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
             kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
             6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) & 
             - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2

            psiNorm = 1.d0/(eps + sqrt(NNorm))
            return
          end function psiNorm

          function psiBulk(kxyz, cmp, m, normflag, f0) result(psi)
          ! Normalized (normflag == TRUE) or not normalized (normflag == FALSE) spinor
          ! cmp = {1,2} --> spinor component of the degenerate pair 
          ! m   = {-1, +1} --> VB / CB 
          ! normflag :: Boolean, whether to retrieve the fully normalized spinor (TRUE) or just the vector as obtained from the Hamiltonian diagonalization with Mathematica (FALSE)
            implicit none
            real*8, dimension(3)     :: kxyz 
            complex*16, dimension(4) :: psi 
            logical                  :: normflag
            real*8                   :: kx, ky, kz, kpar 
            integer                  :: cmp, m 
            real*8                   :: sgn 
            complex*16               :: psi1, psi2, psi3, psi4 
            real*8                   :: NN 
            real*8                   :: f0

            kx = kxyz(1)
            ky = kxyz(2)
            kz = kxyz(3)

            kpar = sqrt( kx**2 + ky**2 )

            sgn = dble(m)

            psi1 = czero
            psi2 = czero
            psi3 = czero
            psi4 = czero

            if(cmp .EQ. 1) then 

              psi1 = M0 + kz**2*M1 + kpar**2*M2 + Sqrt(A0**2*kpar**2 + B0**2*kz**2 + &
                (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                 kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn

              psi2 = B0*kz - dcmplx(0.d0,1.d0)*kx**3*R1 + dcmplx(0.d0,3.d0)*kx*ky**2*R1 +&
               3.d0*kx**2*ky*R2 - ky**3*R2

              psi3 = czero 

              psi4 = A0*(kx + iimag*ky)

            else if ( cmp .EQ. 2 ) then 

              psi1 = czero 

              psi2 = A0*(kx - iimag*ky)

              psi3 = M0 + kz**2*M1 + kpar**2*M2 + Sqrt(A0**2*kpar**2 + B0**2*kz**2 + &
                (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                 kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn

              psi4 = -(B0*kz) - dcmplx(0.d0,1.d0)*kx**3*R1 + dcmplx(0.d0,3.d0)*kx*ky**2*R1 &
              -3.d0*kx**2*ky*R2 + ky**3*R2

            end if 

            NN = psiNorm(kxyz, m, f0)
            psi = (/ psi1, psi2, psi3, psi4 /)
            if( normflag ) psi = NN * psi 

            return 
          end function psiBulk

          function dpsiBulk_dk(kxyz, cmp, m, f0, sigma) result(dpsi_dk)
            ! Bulk wavefunction gradient 
            ! arguments same as above :
            ! kxyz : k-vector in 3D 
            ! cmp = {1,2} --> spinor component of the degenerate pair 
            ! m   = {-1, +1} --> VB / CB 
            ! f0: scaling factor for the regularization parameter eps = (f0 A0)**2
            implicit none
            real*8, dimension(3)       :: kxyz 
            integer                    :: cmp, m 
            complex*16, dimension(4,3) :: dpsi_dk 
            complex*16, dimension(4)   :: dpsi_dkx, dpsi_dky, dpsi_dkz
            real*8                     :: sgn 
            complex*16                 :: dpsi_dkx1, dpsi_dkx2, dpsi_dkx3, dpsi_dkx4 
            complex*16                 :: dpsi_dky1, dpsi_dky2, dpsi_dky3, dpsi_dky4 
            complex*16                 :: dpsi_dkz1, dpsi_dkz2, dpsi_dkz3, dpsi_dkz4 
            complex*16, dimension(4)   :: psi 
            complex*16                 :: dN_dkx, dN_dky, dN_dkz
            real*8                     :: NN 
            real*8                     :: f0, eps 
            real*8                     :: kx, ky, kz 
            real*8, optional           :: sigma 
            real*8                     :: s 

            kx = kxyz(1)
            ky = kxyz(2)
            kz = kxyz(3)

            s = 0.d0 
            if(present(sigma)) s = 1.d0/( ( sigma*4.d0*pi/alat/3.d0 )**2 )


            eps = (f0 * A0)**2 * exp(-(kx**2 + ky**2)* s)


            sgn = dble(m) 

            NN = psiNorm(kxyz, m, f0)

            psi = psiBulk(kxyz, cmp, m, .FALSE.,f0)

            dpsi_dkx1 = czero
            dpsi_dkx2 = czero
            dpsi_dkx3 = czero
            dpsi_dkx4 = czero

            dpsi_dky1 = czero
            dpsi_dky2 = czero
            dpsi_dky3 = czero
            dpsi_dky4 = czero

            dpsi_dkz1 = czero
            dpsi_dkz2 = czero
            dpsi_dkz3 = czero
            dpsi_dkz4 = czero


            if (cmp == 1) then 

              dpsi_dkx1 = 2.d0*kx*M2 + (kx*(A0**2 + 2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                3.d0*kx**4*R1**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2 + &
                2.d0*kx**2*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)/ &
                Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + &
                kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
                6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - &
                2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

              dpsi_dkx2 = dcmplx(0.d0,-3.d0)*(kx - ky)*(kx + ky)*R1 + 6.d0*kx*ky*R2

              dpsi_dkx3 = czero

              dpsi_dkx4 = A0 


              dpsi_dky1 = 2.d0*ky*M2 + ((ky*(A0**2 + 2.d0*M2*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                6.d0*kx**2*(kx**2 - 3.d0*ky**2)*R1**2) + 3.d0*B0*(kx - ky)*(kx + ky)*kz*R2 + 3.d0*ky*(3.d0*kx**4 - &
                4.d0*kx**2*ky**2 + ky**4)*R2**2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) &
                - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

              dpsi_dky2 = dcmplx(0.d0,6.d0)*kx*ky*R1 + 3.d0*kx**2*R2 - 3.d0*ky**2*R2

              dpsi_dky3 = czero

              dpsi_dky4 = iimag * A0 

              dpsi_dkz1 = 2.d0*kz*M1 + ((B0**2*kz + 2.d0*kz*M1*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                B0*ky*(-3.d0*kx**2 + ky**2)*R2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

              dpsi_dkz2 = B0 
              dpsi_dkz3 = czero
              dpsi_dkz4 = czero 

            
            else if (cmp == 2) then

              dpsi_dkx1 = czero 

              dpsi_dkx2 = A0 

              dpsi_dkx3 = 2.d0*kx*M2 + (kx*(A0**2 + 2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                3.d0*kx**4*R1**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2 + &
                2.d0*kx**2*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)/ &
                Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 +&
                kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 -&
                6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - &
                2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

              dpsi_dkx4 = dcmplx(0.d0,-3.d0)*(kx - ky)*(kx + ky)*R1 - 6.d0 *kx*ky*R2

              dpsi_dky1  = czero 

              dpsi_dky2  = -iimag*A0 

              dpsi_dky3  = 2.d0*ky*M2 + ((ky*(A0**2 + 2.d0*M2*(M0 + kz**2*M1 + (kx**2 + &
                ky**2)*M2) - 6.d0*kx**2*(kx**2 - 3.d0*ky**2)*R1**2) + 3.d0*B0*(kx - ky)*(kx + ky)*kz*R2 + &
                3.d0*ky*(3.d0*kx**4 - 4.d0*kx**2*ky**2 + ky**4)*R2**2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + &
                B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + & 
                kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

              dpsi_dky4 = dcmplx(0.d0,6.d0)*kx*ky*R1 - 3.d0*kx**2*R2 + 3.d0*ky**2*R2

              dpsi_dkz1 = czero 

              dpsi_dkz2 = czero 

              dpsi_dkz3 = 2.d0*kz*M1 + ((B0**2*kz + 2.d0*kz*M1*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                B0*ky*(-3.d0*kx**2 + ky**2)*R2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

              dpsi_dkz4 = -B0
            end if 

            dN_dkx =  -(2.d0*A0**2*kx + 4.d0*kx**3*(kx**2 - 3.d0*ky**2)*R1**2 + 2.d0*kx*(kx**2 - 3.d0*ky**2)**2*R1**2 + &
                12.d0*kx*ky*R2*(B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2) + &
                2.d0*(2.d0*kx*M2 + (kx*(A0**2 + 2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                3.d0*kx**4*R1**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2 + &
                2.d0*kx**2*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)/ &
                Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
                kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
                6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + &
                ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2))))*&
                (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn))/&
                (2.d0*(A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2 + &
                (B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2)**2 + &
                (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + &
                B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2)**1.5)

            dN_dky =   -(2.d0*A0**2*ky - 12.d0*kx**2*ky*(kx**2 - 3.d0*ky**2)*R1**2 + &
                6.d0*(kx - ky)*(kx + ky)*R2*(B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2) + &
                2.d0*(2.d0*ky*M2 + ((ky*(A0**2 + 2.d0*M2*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                6.d0*kx**2*(kx**2 - 3.d0*ky**2)*R1**2) + 3.d0*B0*(kx - ky)*(kx + ky)*kz*R2 + &
                3.d0*ky*(3.d0*kx**4 - 4.d0*kx**2*ky**2 + ky**4)*R2**2)*sgn)/ &
                Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
                kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
                6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 +&
                ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2))))*&
                (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2*M0 + 2*kz**2*M1 + ky**2*M2) - 2*B0*ky*kz*R2 + ky**4*R2**2) +&
                kx**4*(M2**2 + 3*ky**2*(-2*R1**2 + 3*R2**2)))*sgn))/ &
                (2.d0*(A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2 + &
                (B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2)**2 + &
                (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + &
                B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 +&
                2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2)**1.5)

            dN_dkz =  -(2.d0*B0*(B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2) + 2.d0*(2.d0*kz*M1 + &
                ((B0**2*kz + 2.d0*kz*M1*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - B0*ky*(-3.d0*kx**2 + ky**2)*R2)*sgn)/&
                Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
                kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 +&
                 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2))))*&
                (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 +&
                2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn))/&
                (2.d0*(A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2 + (B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2)**2 + &
                (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + &
                B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2)**1.5)

            dpsi_dkx = dN_dkx * psi + NN * (/ dpsi_dkx1, dpsi_dkx2, dpsi_dkx3, dpsi_dkx4 /)
            dpsi_dky = dN_dky * psi + NN * (/ dpsi_dky1, dpsi_dky2, dpsi_dky3, dpsi_dky4 /)
            dpsi_dkz = dN_dkz * psi + NN * (/ dpsi_dkz1, dpsi_dkz2, dpsi_dkz3, dpsi_dkz4 /)

            dpsi_dk(:,1) = dpsi_dkx
            dpsi_dk(:,2) = dpsi_dky
            dpsi_dk(:,3) = dpsi_dkz

            return
          end function dpsiBulk_dk

          function DvcB(kxyz, cmpv, cmpc, f0, sigma) result(dipole)
            ! dipole moment between VB and CB (interband)
            ! kxyz : k-vector in 3D
            ! cmpv, cmpc = {1,2} --> spinor component of the degenerate pair for VB and CB, respectively 
            ! f0 : scaling factor for the regularization parameter eps = (f0 A0)**2 to treat discontinuity at Gamma point
            implicit none
            real*8, dimension(3)       :: kxyz
            integer                    :: cmpv, cmpc 
            complex*16, dimension(3)   :: dipole 
            complex*16, dimension(4)   :: psiv
            complex*16, dimension(4,3) :: dpsic_dk
            real*8                     :: f0
            real*8, optional           :: sigma 
            
            psiv = psiBulk(kxyz, cmpv, -1, .TRUE.,0.d0)
            if(present(sigma)) then 
              dpsic_dk = dpsiBulk_dk(kxyz, cmpc, 1, f0, sigma)
            else
              dpsic_dk = dpsiBulk_dk(kxyz, cmpc, 1, f0)
            endif

            dipole(1) = iimag * dot_product(psiv, dpsic_dk(:,1))
            dipole(2) = iimag * dot_product(psiv, dpsic_dk(:,2))
            dipole(3) = iimag * dot_product(psiv, dpsic_dk(:,3))
            return 
          end function DvcB

          function xiB(kxyz, cmp, m, f0,sigma) result (xi)
            ! Berry connection of VB or CB 
            ! kxyz : k-vector in 3D 
            ! cmp  : {1,2} - spinor component of the degenerate pair
            ! m    : {-1,1} - index for VB (-1) or CB (+1)
            ! f0   : regularization parameter scaling factor for treating the singularity at the Gamma point for the CB (M0<0) or the VB (M0>0)
            implicit none
            real*8, dimension(3)       :: kxyz 
            integer                    :: cmp, m 
            complex*16, dimension(3)   :: xi 
            complex*16, dimension(4)   :: psi
            complex*16, dimension(4,3) :: dpsi_dk
            real*8                     :: f0
            real*8, optional           :: sigma 
            
          
            if(present(sigma)) then 
              psi     = psiBulk(kxyz, cmp, m, .TRUE.,f0)
              dpsi_dk = dpsiBulk_dk(kxyz, cmp, m, f0, sigma)
            else
              psi = psiBulk(kxyz, cmp, m, .TRUE.,f0)
              dpsi_dk = dpsiBulk_dk(kxyz, cmp, m, f0)
            endif
            xi(1) = iimag * dot_product(psi, dpsi_dk(:,1))
            xi(2) = iimag * dot_product(psi, dpsi_dk(:,2))
            xi(3) = iimag * dot_product(psi, dpsi_dk(:,3))

            return 
          end function xiB

          function berry(kxyz, cmp, m, f0) result (Omega)
            ! Berry curvature for VB/CB 
            ! kxyz : k-vector in 3D 
            ! cmp  : {1,2} - spinor component of the degenerate pair
            ! m    : {-1,1} - index for VB (-1) or CB (+1)
            ! f0   : regularization parameter scaling factor 

            implicit none 
            real*8, dimension(3)       :: kxyz
            integer                    :: cmp, m 
            complex*16, dimension(3)   :: Omega
            complex*16, dimension(4,3) :: dpsi_dk 
            complex*16, dimension(4)   :: dpsi_dkx, dpsi_dky, dpsi_dkz 
            real*8, optional           :: f0

            if(present(f0)) then 
              dpsi_dk  = dpsiBulk_dk(kxyz, cmp, m, f0)
            else
              dpsi_dk  = dpsiBulk_dk(kxyz, cmp, m, 0.d0)
            endif

              
            dpsi_dkx = dpsi_dk(:,1)
            dpsi_dky = dpsi_dk(:,2)
            dpsi_dkz = dpsi_dk(:,3)

            Omega(1) = iimag * ( dot_product(dpsi_dky, dpsi_dkz) - dot_product(dpsi_dkz, dpsi_dky) )
            Omega(2) = iimag * ( dot_product(dpsi_dkz, dpsi_dkx) - dot_product(dpsi_dkx, dpsi_dkz) )
            Omega(3) = iimag * ( dot_product(dpsi_dkx, dpsi_dky) - dot_product(dpsi_dky, dpsi_dkx) )

            return

          end function berry


          ! Custom implementation of the cross product of two vectors lying in the (xy)-plane, the result is aligned along the z-axis

          real*8 function zcross(v1, v2)
            implicit none
            real*8, dimension(2) :: v1, v2

            zcross = v1(1)*v2(2) - v2(1)*v1(2)
            return
          end function zcross


        

          subroutine get_crystal_kxyz(kxyz, egap, xigap, dmat, fBC_v, fBC_c, sigma)
          ! A routine that gives the band gap (egap), the difference in Berry connections (xigap) and the dipole moment between CB ad VB at a given k-point
          ! Used in the derivs function that calculates the derivatives of the density matrix in the SBE
          ! Inputs:
          ! kxyz: 3-vector, 3D-momentum
          ! fBZ_v, fBC_c: regularization factors for the norm of the VB/CB. The VB does not need to be regularized (M0<0)
          ! sigma: parameter for the local regularization, optional
          ! egap (inout): bandgap
          ! xigap (inout): difference in Berry connections
          ! dmat (inout): dipole moment

            implicit none
            real*8, dimension(3)       :: kxyz
            real*8                     :: egap 
            
            complex*16, dimension(2)   :: dmat, xigap , xi_v, xi_c 
            real*8                     :: fBC_v, fBC_c
            real*8                     :: kx, ky, kz, kpar
            integer                    :: m, cmp
            real*8                     :: sgn, NNorm
            real*8                     :: f0, eps, epsDM
            complex*16, dimension(4)   :: psi 
            complex*16                 :: psi1, psi2, psi3, psi4 
            real*8                     :: NN 
            complex*16, dimension(4,3) :: dpsi_dk 
            complex*16, dimension(4)   :: dpsi_dkx, dpsi_dky, dpsi_dkz
            complex*16                 :: dpsi_dkx1, dpsi_dkx2, dpsi_dkx3, dpsi_dkx4 
            complex*16                 :: dpsi_dky1, dpsi_dky2, dpsi_dky3, dpsi_dky4 
            complex*16                 :: dpsi_dkz1, dpsi_dkz2, dpsi_dkz3, dpsi_dkz4 
            complex*16                 :: dN_dkx, dN_dky, dN_dkz
            integer                    :: cmpv, cmpc 

            real*8                     :: eVB, eCB, eband
            complex*16, dimension(4)   :: psiv, psic 
            complex*16, dimension(4,3) :: dpsiv_dk, dpsic_dk
            real*8                     :: s 
            real*8, optional           :: sigma 

            kx = kxyz(1)
            ky = kxyz(2)
            kz = kxyz(3)

            kpar = sqrt( kx**2 + ky**2 )
            
            kpar = sqrt( kx**2 + ky**2 )
            epsDM = (f0DM * A0)**2
            cmp = icomp 

            s = 0.d0
            if(present(sigma)) s = 1.d0/( ( sigma * 4.d0*pi/3.d0/alat )**2 )
            


            do m = -1, 1, 2 

              sgn = dble(m)

              if( m .eq. -1) then
                f0 = fBC_v 
              else
                f0 = fBC_c 
              endif

              eps = (f0 * A0)**2 * exp( - (kx**2 + ky**2) * s )

              eband = C0 + C2*kpar**2 + C1*kz**2 + Sqrt(A0**2*kpar**2 + 2.d0*A0*A2*kpar**4 + &
              A2**2*kpar**6 + M0**2 + kz**2*((B0 + B2*kz**2)**2 + 2.d0*M0*M1 + kz**2*M1**2) + &
              2.d0*kx**2*M0*M2 + 2.d0*ky**2*M0*M2 + 2.d0*kx**2*kz**2*M1*M2 + 2.d0*ky**2*kz**2*M1*M2 + &
              kx**4.d0*M2**2 + 2.d0*kx**2*ky**2*M2**2 + ky**4*M2**2 + kx**6*R1**2 - 6.d0*kx**4*ky**2*R1**2 + &
              9.d0*kx**2*ky**4*R1**2 - 2.d0*ky*(-3.d0*kx**2 + ky**2)*kz*(B0 + B2*kz**2)*R2 + &
              ky**2*(-3.d0*kx**2 + ky**2)**2*R2**2)*sgn

              ! Caluclate norm of the wavefunction 
              NNorm = A0**2*kpar**2 + B0**2*kz**2 + kx**6*R1**2 - 6.d0*kx**4*ky**2*R1**2 + &
              9.d0*kx**2*ky**4*R1**2 + 6.d0*B0*kx**2*ky*kz*R2 - 2.d0*B0*ky**3*kz*R2 + 9.d0*kx**4*ky**2*R2**2 - &
              6.d0*kx**2*ky**4*R2**2 + ky**6*R2**2 + (M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + &
              Sqrt(A0**2*kpar**2 + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
              kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
              6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) & 
              - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2

              NN = 1.d0/(eps + sqrt(NNorm))

              ! Caluclate  wavefunction 

              if(cmp .EQ. 1) then 

                psi1 = M0 + kz**2*M1 + kpar**2*M2 + Sqrt(A0**2*kpar**2 + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                    2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                   kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn

                psi2 = B0*kz - dcmplx(0.d0,1.d0)*kx**3*R1 + dcmplx(0.d0,3.d0)*kx*ky**2*R1 +&
                 3.d0*kx**2*ky*R2 - ky**3*R2

                psi3 = czero 

                psi4 = A0*(kx + iimag*ky)

              else if ( cmp .EQ. 2 ) then 

                psi1 = czero 

                psi2 = A0*(kx - iimag*ky)

                psi3 = M0 + kz**2*M1 + kpar**2*M2 + Sqrt(A0**2*kpar**2 + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                    2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                   kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn

                psi4 = -(B0*kz) - dcmplx(0.d0,1.d0)*kx**3*R1 + dcmplx(0.d0,3.d0)*kx*ky**2*R1 &
                -3.d0*kx**2*ky*R2 + ky**3*R2

              end if 

              psi = (/ psi1, psi2, psi3, psi4 /)

              ! Caluclate gradient of the wavefunction 

              if (cmp == 1) then 

                dpsi_dkx1 = 2.d0*kx*M2 + (kx*(A0**2 + 2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  3.d0*kx**4*R1**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2 + &
                  2.d0*kx**2*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)/ &
                  Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + &
                  kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
                  6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - &
                  2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dkx2 = dcmplx(0.d0,-3.d0)*(kx - ky)*(kx + ky)*R1 + 6.d0*kx*ky*R2

                dpsi_dkx3 = czero

                dpsi_dkx4 = A0 


                dpsi_dky1 = 2.d0*ky*M2 + ((ky*(A0**2 + 2.d0*M2*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                  6.d0*kx**2*(kx**2 - 3.d0*ky**2)*R1**2) + 3.d0*B0*(kx - ky)*(kx + ky)*kz*R2 + 3.d0*ky*(3.d0*kx**4 - &
                  4.d0*kx**2*ky**2 + ky**4)*R2**2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) &
                  - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dky2 = dcmplx(0.d0,6.d0)*kx*ky*R1 + 3.d0*kx**2*R2 - 3.d0*ky**2*R2

                dpsi_dky3 = czero

                dpsi_dky4 = iimag * A0 

                dpsi_dkz1 = 2.d0*kz*M1 + ((B0**2*kz + 2.d0*kz*M1*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                  B0*ky*(-3.d0*kx**2 + ky**2)*R2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dkz2 = B0 
                dpsi_dkz3 = czero
                dpsi_dkz4 = czero 

              
              else if (cmp == 2) then

                dpsi_dkx1 = czero 

                dpsi_dkx2 = A0 

                dpsi_dkx3 = 2.d0*kx*M2 + (kx*(A0**2 + 2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  3.d0*kx**4*R1**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2 + &
                  2.d0*kx**2*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)/ &
                  Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 +&
                  kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 -&
                  6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - &
                  2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dkx4 = dcmplx(0.d0,-3.d0)*(kx - ky)*(kx + ky)*R1 - 6.d0 *kx*ky*R2

                dpsi_dky1  = czero 

                dpsi_dky2  = -iimag*A0 

                dpsi_dky3  = 2.d0*ky*M2 + ((ky*(A0**2 + 2.d0*M2*(M0 + kz**2*M1 + (kx**2 + &
                  ky**2)*M2) - 6.d0*kx**2*(kx**2 - 3.d0*ky**2)*R1**2) + 3.d0*B0*(kx - ky)*(kx + ky)*kz*R2 + &
                  3.d0*ky*(3.d0*kx**4 - 4.d0*kx**2*ky**2 + ky**4)*R2**2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + &
                  B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + & 
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dky4 = dcmplx(0.d0,6.d0)*kx*ky*R1 - 3.d0*kx**2*R2 + 3.d0*ky**2*R2

                dpsi_dkz1 = czero 

                dpsi_dkz2 = czero 

                dpsi_dkz3 = 2.d0*kz*M1 + ((B0**2*kz + 2.d0*kz*M1*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                  B0*ky*(-3.d0*kx**2 + ky**2)*R2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dkz4 = -B0
              end if 

              dN_dkx =  -(2.d0*A0**2*kx + 4.d0*kx**3*(kx**2 - 3.d0*ky**2)*R1**2 + 2.d0*kx*(kx**2 - 3.d0*ky**2)**2*R1**2 + &
                  12.d0*kx*ky*R2*(B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2) + &
                  2.d0*(2.d0*kx*M2 + (kx*(A0**2 + 2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  3.d0*kx**4*R1**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2 + &
                  2.d0*kx**2*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)/ &
                  Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
                  kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
                  6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + &
                  ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2))))*&
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn))/&
                  (2.d0*(A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2 + &
                  (B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2)**2 + &
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + &
                  B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2)**1.5)

              dN_dky =   -(2.d0*A0**2*ky - 12.d0*kx**2*ky*(kx**2 - 3.d0*ky**2)*R1**2 + &
                  6.d0*(kx - ky)*(kx + ky)*R2*(B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2) + &
                  2.d0*(2.d0*ky*M2 + ((ky*(A0**2 + 2.d0*M2*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                  6.d0*kx**2*(kx**2 - 3.d0*ky**2)*R1**2) + 3.d0*B0*(kx - ky)*(kx + ky)*kz*R2 + &
                  3.d0*ky*(3.d0*kx**4 - 4.d0*kx**2*ky**2 + ky**4)*R2**2)*sgn)/ &
                  Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
                  kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
                  6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 +&
                  ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2))))*&
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2*M0 + 2*kz**2*M1 + ky**2*M2) - 2*B0*ky*kz*R2 + ky**4*R2**2) +&
                  kx**4*(M2**2 + 3*ky**2*(-2*R1**2 + 3*R2**2)))*sgn))/ &
                  (2.d0*(A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2 + &
                  (B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2)**2 + &
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + &
                  B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 +&
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2)**1.5)

              dN_dkz =  -(2.d0*B0*(B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2) + 2.d0*(2.d0*kz*M1 + &
                  ((B0**2*kz + 2.d0*kz*M1*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - B0*ky*(-3.d0*kx**2 + ky**2)*R2)*sgn)/&
                  Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
                  kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 +&
                   9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2))))*&
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 +&
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn))/&
                  (2.d0*(A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2 + (B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2)**2 + &
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + &
                  B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2)**1.5)

              dpsi_dkx = dN_dkx * psi + NN * (/ dpsi_dkx1, dpsi_dkx2, dpsi_dkx3, dpsi_dkx4 /)
              dpsi_dky = dN_dky * psi + NN * (/ dpsi_dky1, dpsi_dky2, dpsi_dky3, dpsi_dky4 /)
              dpsi_dkz = dN_dkz * psi + NN * (/ dpsi_dkz1, dpsi_dkz2, dpsi_dkz3, dpsi_dkz4 /)

              dpsi_dk(:,1) = dpsi_dkx
              dpsi_dk(:,2) = dpsi_dky
              dpsi_dk(:,3) = dpsi_dkz

              if( m .eq. -1) then 
                psiv = NN * psi
                dpsiv_dk = dpsi_dk 
                eVB = eband 
              else
                psic = NN * psi
                dpsic_dk = dpsi_dk
                eCB = eband 
              endif

            enddo

            egap = eCB - eVB

            dmat(1) = iimag * dot_product(psiv, dpsic_dk(:,1))
            dmat(2) = iimag * dot_product(psiv, dpsic_dk(:,2))            

            xigap(1) =  iimag * ( dot_product(psic, dpsic_dk(:,1)) - &
                                 dot_product(psiv, dpsiv_dk(:,1)) ) 
            xigap(2) =  iimag * ( dot_product(psic, dpsic_dk(:,2)) - &
                                 dot_product(psiv, dpsiv_dk(:,2)) )

          end subroutine get_crystal_kxyz


          subroutine get_crystal_velocities_kxyz(kxyz, vgr_v, vgr_c, omg_v, omg_c, dmat, fBC_v, fBC_c, sigma)
            ! This subroutine is invoked in the current calculation step and gives the group velocities of the VB and CB as well as the Berry curvatures
            ! Note: The regularization procedure actually introduces a discontinuity in the CB which was not present before
            ! Therefore, call this subroutine with fBC_v = fBC_c = 0
            ! Arguments:
            ! kxyz: 3D momentum
            ! vgr_v, vgr_c : group velocities of VB/CB
            ! omg_v, omg_c : Berry curvatures of VB/CB 
            ! dmat : dipole moment (WARNING: will be discontinuous when fBC_v=fBC_v=0, re-compute with normal routine)
            ! sigma: optional width of the local regularization factor 
            implicit none
            real*8, dimension(3)       :: kxyz            
            complex*16, dimension(2)   :: dmat
            complex*16                 :: omg_v, omg_c 
            real*8, dimension(2)       :: vgr_v, vgr_c 
            real*8                     :: fBC_v, fBC_c
            real*8                     :: kx, ky, kz, kpar
            real*8                     :: dE_dkx, dE_dky, dE_dkz
            integer                    :: m, cmp
            real*8                     :: sgn, NNorm
            real*8                     :: f0, eps, epsDM
            complex*16, dimension(4)   :: psi 
            complex*16                 :: psi1, psi2, psi3, psi4 
            real*8                     :: NN 
            complex*16, dimension(4,3) :: dpsi_dk 
            complex*16, dimension(4)   :: dpsi_dkx, dpsi_dky, dpsi_dkz
            complex*16                 :: dpsi_dkx1, dpsi_dkx2, dpsi_dkx3, dpsi_dkx4 
            complex*16                 :: dpsi_dky1, dpsi_dky2, dpsi_dky3, dpsi_dky4 
            complex*16                 :: dpsi_dkz1, dpsi_dkz2, dpsi_dkz3, dpsi_dkz4 
            complex*16                 :: dN_dkx, dN_dky, dN_dkz
            integer                    :: cmpv, cmpc             
            complex*16, dimension(4)   :: psiv, psic 
            complex*16, dimension(4,3) :: dpsiv_dk, dpsic_dk
            real*8                     :: s 
            real*8, optional           :: sigma 

            kx = kxyz(1)
            ky = kxyz(2)
            kz = kxyz(3)

            kpar = sqrt( kx**2 + ky**2 )
            
            kpar = sqrt( kx**2 + ky**2 )
            epsDM = (f0DM * A0)**2
            cmp = icomp 
            
            s = 0.d0
            if(present(sigma)) s = 1.d0/( ( sigma * 4.d0*pi/3.d0/alat )**2 )


            do m = -1, 1, 2 

              sgn = dble(m)

              if( m .eq. -1) then
                f0 = fBC_v 
              else
                f0 = fBC_c 
              endif

              eps = (f0 * A0)**2 * exp( -(kx**2 + ky**2) * s )

              dE_dkx = 2.d0*C2*kx + (kx*((A0 + A2*kpar**2)*(A0 + 3.d0*A2*kpar**2) + 2.d0*M0*M2 + &
                2.d0*kz**2*M1*M2 + 2.d0*kx**2*M2**2 + 2.d0*ky**2*M2**2 + 3.d0*kx**4*R1**2 - &
                12.d0*kx**2*ky**2*R1**2 + 9.d0*ky**4*R1**2 + 6.d0*ky*kz*(B0 + B2*kz**2)*R2 - &
                6.d0*ky**2*(-3.d0*kx**2 + ky**2)*R2**2)*sgn)/Sqrt(A0**2*kpar**2 + 2.d0*A0*A2*kpar**4 +&
                A2**2*kpar**6 + M0**2 + kz**2*((B0 + B2*kz**2)**2 + 2.d0*M0*M1 + kz**2*M1**2) + 2.d0*kx**2*M0*M2 + &
                2.d0*ky**2*M0*M2 + 2.d0*kx**2*kz**2*M1*M2 + 2.d0*ky**2*kz**2*M1*M2 + kx**4*M2**2 + &
                2.d0*kx**2*ky**2*M2**2 + ky**4*M2**2 + kx**6*R1**2 - 6.d0*kx**4*ky**2*R1**2 + &
                9.d0*kx**2*ky**4*R1**2 - 2.d0*ky*(-3.d0*kx**2 + ky**2)*kz*(B0 + B2*kz**2)*R2 + ky**2*(-3.d0*kx**2 + ky**2)**2*R2**2)


              dE_dky = 2.d0*C2*ky + ((ky*((A0 + A2*kpar**2)*(A0 + 3.d0*A2*kpar**2) + 2.d0*M2*(M0 + kz**2*M1 + &
                (kx**2 + ky**2)*M2) - 6.d0*kx**2*(kx**2 - 3.d0*ky**2)*R1**2) + &
                3.d0*(kx - ky)*(kx + ky)*kz*(B0 + B2*kz**2)*R2 + &
                3.d0*ky*(3.d0*kx**4 - 4.d0*kx**2*ky**2 + ky**4)*R2**2)*sgn)/&
                Sqrt(kpar**2*(A0 + A2*kpar**2)**2 + B0**2*kz**2 + 2.d0*B0*B2*kz**4 + &
                B2**2*kz**6 + M0**2 + 2.d0*kz**2*M0*M1 + kz**4*M1**2 + 2.d0*kx**2*M0*M2 +&
                2.d0*ky**2*M0*M2 + 2.d0*kx**2*kz**2*M1*M2 + 2.d0*ky**2*kz**2*M1*M2 + kx**4*M2**2 + &
                2.d0*kx**2*ky**2*M2**2 + ky**4*M2**2 + kx**6*R1**2 - 6.d0*kx**4*ky**2*R1**2 + &
                9.d0*kx**2*ky**4*R1**2 - 2.d0*ky*(-3.d0*kx**2 + ky**2)*kz*(B0 + B2*kz**2)*R2 + &
                ky**2*(-3.d0*kx**2 + ky**2)**2*R2**2)

              dE_dkz = 2.d0*C1*kz + ((kz*(B0**2 + 4.d0*B0*B2*kz**2 + 3.d0*B2**2*kz**4 + &
                2.d0*M1*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2)) - &
                ky*(-3.d0*kx**2 + ky**2)*(B0 + 3.d0*B2*kz**2)*R2)*sgn)/&
                Sqrt(kpar**2*(A0 + A2*kpar**2)**2 + B0**2*kz**2 + 2.d0*B0*B2*kz**4 + &
                B2**2*kz**6 + M0**2 + 2.d0*kz**2*M0*M1 + kz**4*M1**2 + 2.d0*kx**2*M0*M2 + &
                2.d0*ky**2*M0*M2 + 2.d0*kx**2*kz**2*M1*M2 + 2.d0*ky**2*kz**2*M1*M2 + &
                kx**4*M2**2 + 2.d0*kx**2*ky**2*M2**2 + ky**4*M2**2 + kx**6*R1**2 - &
                6.d0*kx**4*ky**2*R1**2 + 9.d0*kx**2*ky**4*R1**2 - 2.d0*ky*(-3.d0*kx**2 + ky**2)*kz*(B0 + B2*kz**2)*R2 + &
                ky**2*(-3.d0*kx**2 + ky**2)**2*R2**2)
  

              ! Caluclate norm of the wavefunction 
              NNorm = A0**2*kpar**2 + B0**2*kz**2 + kx**6*R1**2 - 6.d0*kx**4*ky**2*R1**2 + &
              9.d0*kx**2*ky**4*R1**2 + 6.d0*B0*kx**2*ky*kz*R2 - 2.d0*B0*ky**3*kz*R2 + 9.d0*kx**4*ky**2*R2**2 - &
              6.d0*kx**2*ky**4*R2**2 + ky**6*R2**2 + (M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + &
              Sqrt(A0**2*kpar**2 + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
              kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
              6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) & 
              - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2

              NN = 1.d0/(eps + sqrt(NNorm))

              ! Caluclate  wavefunction 

              if(cmp .EQ. 1) then 

                psi1 = M0 + kz**2*M1 + kpar**2*M2 + Sqrt(A0**2*kpar**2 + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                    2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                   kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn

                psi2 = B0*kz - dcmplx(0.d0,1.d0)*kx**3*R1 + dcmplx(0.d0,3.d0)*kx*ky**2*R1 +&
                 3.d0*kx**2*ky*R2 - ky**3*R2

                psi3 = czero 

                psi4 = A0*(kx + iimag*ky)

              else if ( cmp .EQ. 2 ) then 

                psi1 = czero 

                psi2 = A0*(kx - iimag*ky)

                psi3 = M0 + kz**2*M1 + kpar**2*M2 + Sqrt(A0**2*kpar**2 + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                    2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                   kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn

                psi4 = -(B0*kz) - dcmplx(0.d0,1.d0)*kx**3*R1 + dcmplx(0.d0,3.d0)*kx*ky**2*R1 &
                -3.d0*kx**2*ky*R2 + ky**3*R2

              end if 

              psi = (/ psi1, psi2, psi3, psi4 /)

              ! Caluclate gradient of the wavefunction 

              if (cmp == 1) then 

                dpsi_dkx1 = 2.d0*kx*M2 + (kx*(A0**2 + 2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  3.d0*kx**4*R1**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2 + &
                  2.d0*kx**2*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)/ &
                  Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + &
                  kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
                  6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - &
                  2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dkx2 = dcmplx(0.d0,-3.d0)*(kx - ky)*(kx + ky)*R1 + 6.d0*kx*ky*R2

                dpsi_dkx3 = czero

                dpsi_dkx4 = A0 


                dpsi_dky1 = 2.d0*ky*M2 + ((ky*(A0**2 + 2.d0*M2*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                  6.d0*kx**2*(kx**2 - 3.d0*ky**2)*R1**2) + 3.d0*B0*(kx - ky)*(kx + ky)*kz*R2 + 3.d0*ky*(3.d0*kx**4 - &
                  4.d0*kx**2*ky**2 + ky**4)*R2**2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) &
                  - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dky2 = dcmplx(0.d0,6.d0)*kx*ky*R1 + 3.d0*kx**2*R2 - 3.d0*ky**2*R2

                dpsi_dky3 = czero

                dpsi_dky4 = iimag * A0 

                dpsi_dkz1 = 2.d0*kz*M1 + ((B0**2*kz + 2.d0*kz*M1*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                  B0*ky*(-3.d0*kx**2 + ky**2)*R2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dkz2 = B0 
                dpsi_dkz3 = czero
                dpsi_dkz4 = czero 

              
              else if (cmp == 2) then

                dpsi_dkx1 = czero 

                dpsi_dkx2 = A0 

                dpsi_dkx3 = 2.d0*kx*M2 + (kx*(A0**2 + 2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  3.d0*kx**4*R1**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2 + &
                  2.d0*kx**2*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)/ &
                  Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 +&
                  kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 -&
                  6.d0*ky**4*R2**2) + ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - &
                  2.d0*B0*ky*kz*R2 + ky**4*R2**2) + kx**4*(M2**2 + 3*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dkx4 = dcmplx(0.d0,-3.d0)*(kx - ky)*(kx + ky)*R1 - 6.d0 *kx*ky*R2

                dpsi_dky1  = czero 

                dpsi_dky2  = -iimag*A0 

                dpsi_dky3  = 2.d0*ky*M2 + ((ky*(A0**2 + 2.d0*M2*(M0 + kz**2*M1 + (kx**2 + &
                  ky**2)*M2) - 6.d0*kx**2*(kx**2 - 3.d0*ky**2)*R1**2) + 3.d0*B0*(kx - ky)*(kx + ky)*kz*R2 + &
                  3.d0*ky*(3.d0*kx**4 - 4.d0*kx**2*ky**2 + ky**4)*R2**2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + &
                  B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + & 
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dky4 = dcmplx(0.d0,6.d0)*kx*ky*R1 - 3.d0*kx**2*R2 + 3.d0*ky**2*R2

                dpsi_dkz1 = czero 

                dpsi_dkz2 = czero 

                dpsi_dkz3 = 2.d0*kz*M1 + ((B0**2*kz + 2.d0*kz*M1*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                  B0*ky*(-3.d0*kx**2 + ky**2)*R2)*sgn)/Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))

                dpsi_dkz4 = -B0
              end if 

              dN_dkx =  -(2.d0*A0**2*kx + 4.d0*kx**3*(kx**2 - 3.d0*ky**2)*R1**2 + 2.d0*kx*(kx**2 - 3.d0*ky**2)**2*R1**2 + &
                  12.d0*kx*ky*R2*(B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2) + &
                  2.d0*(2.d0*kx*M2 + (kx*(A0**2 + 2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  3.d0*kx**4*R1**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2 + &
                  2.d0*kx**2*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)/ &
                  Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
                  kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
                  6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + &
                  ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2))))*&
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn))/&
                  (2.d0*(A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2 + &
                  (B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2)**2 + &
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + &
                  B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2)**1.5)

              dN_dky =   -(2.d0*A0**2*ky - 12.d0*kx**2*ky*(kx**2 - 3.d0*ky**2)*R1**2 + &
                  6.d0*(kx - ky)*(kx + ky)*R2*(B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2) + &
                  2.d0*(2.d0*ky*M2 + ((ky*(A0**2 + 2.d0*M2*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - &
                  6.d0*kx**2*(kx**2 - 3.d0*ky**2)*R1**2) + 3.d0*B0*(kx - ky)*(kx + ky)*kz*R2 + &
                  3.d0*ky*(3.d0*kx**4 - 4.d0*kx**2*ky**2 + ky**4)*R2**2)*sgn)/ &
                  Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
                  kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + &
                  6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 +&
                  ky**4*R2**2) + kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2))))*&
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 + &
                  9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2*M0 + 2*kz**2*M1 + ky**2*M2) - 2*B0*ky*kz*R2 + ky**4*R2**2) +&
                  kx**4*(M2**2 + 3*ky**2*(-2*R1**2 + 3*R2**2)))*sgn))/ &
                  (2.d0*(A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2 + &
                  (B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2)**2 + &
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + &
                  B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 +&
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2)**1.5)

              dN_dkz =  -(2.d0*B0*(B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2) + 2.d0*(2.d0*kz*M1 + &
                  ((B0**2*kz + 2.d0*kz*M1*(M0 + kz**2*M1 + (kx**2 + ky**2)*M2) - B0*ky*(-3.d0*kx**2 + ky**2)*R2)*sgn)/&
                  Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + &
                  kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + 2.d0*ky**2*M2**2 +&
                   9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2))))*&
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + B0**2*kz**2 + &
                  (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 +&
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) +&
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn))/&
                  (2.d0*(A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2 + (B0*kz + 3.d0*kx**2*ky*R2 - ky**3*R2)**2 + &
                  (eps + M0 + kz**2*M1 + (kx**2 + ky**2)*M2 + Sqrt(A0**2*(kx**2 + ky**2) + &
                  B0**2*kz**2 + (M0 + kz**2*M1)**2 + kx**6*R1**2 + kx**2*(2.d0*M0*M2 + 2.d0*kz**2*M1*M2 + &
                  2.d0*ky**2*M2**2 + 9.d0*ky**4*R1**2 + 6.d0*B0*ky*kz*R2 - 6.d0*ky**4*R2**2) + &
                  ky**2*(M2*(2.d0*M0 + 2.d0*kz**2*M1 + ky**2*M2) - 2.d0*B0*ky*kz*R2 + ky**4*R2**2) + &
                  kx**4*(M2**2 + 3.d0*ky**2*(-2.d0*R1**2 + 3.d0*R2**2)))*sgn)**2)**1.5)

              dpsi_dkx = dN_dkx * psi + NN * (/ dpsi_dkx1, dpsi_dkx2, dpsi_dkx3, dpsi_dkx4 /)
              dpsi_dky = dN_dky * psi + NN * (/ dpsi_dky1, dpsi_dky2, dpsi_dky3, dpsi_dky4 /)
              dpsi_dkz = dN_dkz * psi + NN * (/ dpsi_dkz1, dpsi_dkz2, dpsi_dkz3, dpsi_dkz4 /)

              dpsi_dk(:,1) = dpsi_dkx
              dpsi_dk(:,2) = dpsi_dky
              dpsi_dk(:,3) = dpsi_dkz            

              if( m .eq. -1) then 
                psiv = NN * psi
                dpsiv_dk = dpsi_dk 
                vgr_v = (/ dE_dkx, dE_dky /)
                omg_v = iimag * ( dot_product(dpsi_dkx, dpsi_dky) - dot_product(dpsi_dky, dpsi_dkx) )
              else
                psic = NN * psi
                dpsic_dk = dpsi_dk
                vgr_c = (/ dE_dkx, dE_dky /)
                omg_c = iimag * ( dot_product(dpsi_dkx, dpsi_dky) - dot_product(dpsi_dky, dpsi_dkx) )   
              endif

            enddo            

            dmat(1) = iimag * dot_product(psiv, dpsic_dk(:,1))
            dmat(2) = iimag * dot_product(psiv, dpsic_dk(:,2))            

          end subroutine get_crystal_velocities_kxyz
            


          subroutine crystal_dealloc()
          ! Deallocate common arrays to release memory 
            implicit none
            if(allocated(kx1d)) deallocate(kx1d)
            if(allocated(ky1d)) deallocate(ky1d)
            if(allocated(kxy_grid)) deallocate(kxy_grid)
            if(allocated(kmask)) deallocate(kmask)
            
            return
          end subroutine crystal_dealloc

          subroutine print_snapshot(iunit, n_cc, pi_vc, t)
            implicit none
            integer :: iunit
            complex*16, dimension(nkpts) :: n_cc, pi_vc 
            real*8   :: t 
            integer :: i 

            do i = 1, nkpts 
              write(7000+iunit, 2055) t, kxy_grid(i,1), kxy_grid(i,2), &
              dble(n_cc(i)), aimag(n_cc(i))
              write(9000+iunit, 2055) t, kxy_grid(i,1), kxy_grid(i,2), &
              dble(pi_vc(i)), aimag(pi_vc(i))
              if(mod(i, nky).eq.0) then 
                write(7000+iunit,'(A)') ' ' 
                write(9000+iunit,'(A)') ' ' 
              endif
            enddo
            return
 2055     format (1e16.8,4(2x,1e16.8))  
          end subroutine print_snapshot

          subroutine print_bands()
          ! Print all Bi2Se3 specific quantities, i.e. energy spectra of VB/CB, group velocities, dipoles, Berry connections 
          ! The output format is GNUPLOT-friendly (to be used with splot)
            implicit none
            real*8    , allocatable, dimension(:,:)   :: Ebulk_VB, Ebulk_CB
            real*8    , allocatable, dimension(:,:,:) :: vBulk_VB, vBulk_CB
            real*8    , allocatable, dimension(:,:,:) :: xiBulk_VB1, xiBulk_CB1
            real*8    , allocatable, dimension(:,:,:) :: xiBulk_VB2, xiBulk_CB2
            real*8    , allocatable, dimension(:,:,:) :: berry_VB1, berry_CB1
            real*8    , allocatable, dimension(:,:,:) :: berry_VB2, berry_CB2
            complex*16, allocatable, dimension(:,:,:) :: dipoleVC1, dipoleVC2
            integer :: i, j
            real*8, dimension(3) :: kxyz 
            real*8 :: kx, ky, kz  

            if(.not.allocated(Ebulk_VB)) allocate(Ebulk_VB(nkx, nky))
            if(.not.allocated(Ebulk_CB)) allocate(Ebulk_CB(nkx, nky))
            if(.not.allocated(vBulk_VB)) allocate(vBulk_VB(nkx, nky, 3))
            if(.not.allocated(vBulk_CB)) allocate(vBulk_CB(nkx, nky, 3))
            if(.not.allocated(xiBulk_VB1)) allocate(xiBulk_VB1(nkx, nky, 3))
            if(.not.allocated(xiBulk_VB2)) allocate(xiBulk_VB2(nkx, nky, 3))
            if(.not.allocated(xiBulk_CB1)) allocate(xiBulk_CB1(nkx, nky, 3))
            if(.not.allocated(xiBulk_CB2)) allocate(xiBulk_CB2(nkx, nky, 3))
            if(.not.allocated(berry_VB1)) allocate(berry_VB1(nkx, nky, 3))
            if(.not.allocated(berry_VB2)) allocate(berry_VB2(nkx, nky, 3))
            if(.not.allocated(berry_CB1)) allocate(berry_CB1(nkx, nky, 3))
            if(.not.allocated(berry_CB2)) allocate(berry_CB2(nkx, nky, 3))
            if(.not.allocated(dipoleVC1)) allocate(dipoleVC1(nkx, nky, 3))
            if(.not.allocated(dipoleVC2)) allocate(dipoleVC2(nkx, nky, 3))
            
            kz = 0.d0         

            do j = 1, nky 
              ky = ky1d(j)
              do i = 1, nkx 
                kx = kx1d(i)
                kxyz = (/ kx, ky, kz /)

                Ebulk_VB(i,j) = Ebulk(kxyz, -1)
                Ebulk_CB(i,j) = Ebulk(kxyz,  1)

                vBulk_VB(i,j,:) = dEbulk_dk(kxyz, -1)
                vBulk_CB(i,j,:) = dEbulk_dk(kxyz,  1)

                if( M0 .lt. 0.d0) then 
                  xiBulk_VB1(i,j,:) = dble(xiB( kxyz, 1, -1, 0.d0 ))
                  xiBulk_VB2(i,j,:) = dble(xiB( kxyz, 2, -1, 0.d0 ))
                  xiBulk_CB1(i,j,:) = dble(xiB( kxyz, 1,  1, f0BC, sigp ))
                  xiBulk_CB2(i,j,:) = dble(xiB( kxyz, 2,  1, f0BC, sigp ))
                else
                  xiBulk_VB1(i,j,:) = dble(xiB( kxyz, 1, -1, f0BC))
                  xiBulk_VB2(i,j,:) = dble(xiB( kxyz, 2, -1, f0BC))
                  xiBulk_CB1(i,j,:) = dble(xiB( kxyz, 1,  1, 0.d0 ))
                  xiBulk_CB2(i,j,:) = dble(xiB( kxyz, 2,  1, 0.d0 ))
                endif
                
                berry_VB1(i,j,:) = dble(berry(kxyz, 1, -1))
                berry_VB2(i,j,:) = dble(berry(kxyz, 2, -1))
                berry_CB1(i,j,:) = dble(berry(kxyz, 1,  1))
                berry_CB2(i,j,:) = dble(berry(kxyz, 2,  1))

                dipoleVC1(i,j,:) = DvcB( kxyz, 1, 1, f0DM )
                dipoleVC2(i,j,:) = DvcB( kxyz, 2, 2, f0DM )

              enddo
            enddo

            open(351,file=Opath(1:Olast)//'Ebulk_VB.dat', status = "replace")
            open(353,file=Opath(1:Olast)//'Ebulk_CB.dat', status = "replace")          

            open(355,file=Opath(1:Olast)//'vBulk_VB.dat', status = "replace")
            open(357,file=Opath(1:Olast)//'vBulk_CB.dat', status = "replace")

            open(371,file=Opath(1:Olast)//'xiBulk_VB1.dat', status = "replace")
            open(373,file=Opath(1:Olast)//'xiBulk_CB1.dat', status = "replace")
            open(375,file=Opath(1:Olast)//'xiBulk_VB2.dat', status = "replace")
            open(377,file=Opath(1:Olast)//'xiBulk_CB2.dat', status = "replace")

            open(391,file=Opath(1:Olast)//'berry_VB1.dat', status = "replace")
            open(393,file=Opath(1:Olast)//'berry_CB1.dat', status = "replace")
            open(395,file=Opath(1:Olast)//'berry_VB2.dat', status = "replace")
            open(397,file=Opath(1:Olast)//'berry_CB2.dat', status = "replace")

            open(511,file=Opath(1:Olast)//'dipoleVC1.dat', status = "replace")
            open(513,file=Opath(1:Olast)//'dipoleVC2.dat', status = "replace")

            open(171,file=Opath(1:Olast)//'kxgrid.dat', status = "replace")
            open(173,file=Opath(1:Olast)//'kygrid.dat', status = "replace")

            open(175,file=Opath(1:Olast)//'kmask.dat', status = "replace")


            do i = 1, nkx
              write(171, 5001) kx1d(i)
            enddo

            do j = 1, nky               
              write(173, 5001) ky1d(j)
            enddo

            close(171)
            close(173)

            do i = 1, nkx
              kx = kx1d(i)
              do j = 1, nky 
                ky = ky1d(j)

                write(351, 5002) kx, ky, Ebulk_VB(i,j)
                write(353, 5002) kx, ky, Ebulk_CB(i,j)

                write(355, 5005) kx, ky, vBulk_VB(i,j,1), vBulk_VB(i,j,2), vBulk_VB(i,j,3)
                write(357, 5005) kx, ky, vBulk_CB(i,j,1), vBulk_CB(i,j,2), vBulk_CB(i,j,3)

                write(371, 5005) kx, ky, xiBulk_VB1(i,j,1), xiBulk_VB1(i,j,2), xiBulk_VB1(i,j,3)
                write(373, 5005) kx, ky, xiBulk_CB1(i,j,1), xiBulk_CB1(i,j,2), xiBulk_CB1(i,j,3)
                write(375, 5005) kx, ky, xiBulk_VB2(i,j,1), xiBulk_VB2(i,j,2), xiBulk_VB2(i,j,3)
                write(377, 5005) kx, ky, xiBulk_CB2(i,j,1), xiBulk_CB2(i,j,2), xiBulk_CB2(i,j,3)

                write(391, 5005) kx, ky, berry_VB1(i,j,1), berry_VB1(i,j,2), berry_VB1(i,j,3)
                write(393, 5005) kx, ky, berry_CB1(i,j,1), berry_CB1(i,j,2), berry_CB1(i,j,3)
                write(395, 5005) kx, ky, berry_VB2(i,j,1), berry_VB2(i,j,2), berry_VB2(i,j,3)
                write(397, 5005) kx, ky, berry_CB2(i,j,1), berry_CB2(i,j,2), berry_CB2(i,j,3)

                write(511, 5004) kx, ky, dble(dipoleVC1(i,j,1)), dimag(dipoleVC1(i,j,1)), &
                dble(dipoleVC1(i,j,2)), dimag(dipoleVC1(i,j,2)), dble(dipoleVC1(i,j,3)), dimag(dipoleVC1(i,j,3))
                write(513, 5004) kx, ky, dble(dipoleVC2(i,j,1)), dimag(dipoleVC2(i,j,1)), &
                dble(dipoleVC2(i,j,2)), dimag(dipoleVC2(i,j,2)), dble(dipoleVC2(i,j,3)), dimag(dipoleVC2(i,j,3))

                write(175, 5002) kx, ky, kmask( (i-1)*nky + j )

              enddo

              if(nky .gt. 1) then
                write(351,'(A)') ' '
                write(353,'(A)') ' '
                write(355,'(A)') ' '
                write(357,'(A)') ' '
                write(371,'(A)') ' '
                write(373,'(A)') ' '
                write(375,'(A)') ' '
                write(377,'(A)') ' '
                write(391,'(A)') ' '
                write(393,'(A)') ' '
                write(395,'(A)') ' '
                write(397,'(A)') ' '
                write(511,'(A)') ' '
                write(513,'(A)') ' '
                write(175,'(A)') ' '
              endif
            
            enddo

            close(351)
            close(353)
            close(355)
            close(357)
            close(371)
            close(373)
            close(375)
            close(377)
            close(391)
            close(393)
            close(395)
            close(397)
            close(511)
            close(513)

            close(175)


            if( allocated(Ebulk_VB))    deallocate(Ebulk_VB  )
            if( allocated(Ebulk_CB))    deallocate(Ebulk_CB  )
            if( allocated(vBulk_VB))    deallocate(vBulk_VB  )
            if( allocated(vBulk_CB))    deallocate(vBulk_CB  )
            if( allocated(xiBulk_VB1))  deallocate(xiBulk_VB1)
            if( allocated(xiBulk_VB2))  deallocate(xiBulk_VB2)
            if( allocated(xiBulk_CB1))  deallocate(xiBulk_CB1)
            if( allocated(xiBulk_CB2))  deallocate(xiBulk_CB2)
            if( allocated(berry_VB1))   deallocate(berry_VB1 )
            if( allocated(berry_VB2))   deallocate(berry_VB2 )
            if( allocated(berry_CB1))   deallocate(berry_CB1 )
            if( allocated(berry_CB2))   deallocate(berry_CB2 )
            if( allocated(dipoleVC1))   deallocate(dipoleVC1 )
            if( allocated(dipoleVC2))   deallocate(dipoleVC2 )

 5001     format (1e16.8)  
 5002     format (1e16.8,3(2x,1e16.8))  
 5003     format (1e16.8,2x,1e16.8)  
 5004     format (1f16.8,6(2x,1e16.8))  
 5005     format (1f16.8,4(2x,1e16.8))

          end subroutine print_bands

        end module solid
