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
          real*8 :: alpha1, alpha3                                  ! Parameters for the Bi2Se3 Hamiltonian, from Liu et al. (2010)
          real*8 :: C0t, C2t, Atld, Rt                                ! Parameters for the effective 2D surface Bi2Se3 Hamiltonian, from Liu et al. (2010), see Eq. 34
          real*8 :: alat                                            ! Lattice constant a of the rhombohedral lattice 
          real*8 :: Tdephas_2                                       ! Dephaising time T2, in a.u. 
          real*8, allocatable, dimension(:,:) :: kxy_grid           ! 2D array holding all nkx * nky grid points (nkx * nky, 2)
          real*8, allocatable, dimension(:)   :: kx1d, ky1d         ! 1D arrays of k-points in x and y dimensions
          real*8, allocatable, dimension(:)   :: kmask
          real*8                                    :: sigp         ! Width of the Gaussian in the local regularization factor. Denotes the fraction of the BZ to be used as a widt, i.e. 
                                                                    ! the actual width sigma is given by: sigma = sigp * 4*pi/3/alat. 4 pi/3/alat is the extent of the BZ in the Gamma-K direction 
          real*8                              :: bzcut, f_dmask, f_kmask, q_kmask
          real*8                              :: f0glob            

        

          save Tdephas_2, kxy_grid, kx1d, ky1d, &
          dkx,dky, alat, bzcut, f0glob, q_kmask, &
           A0, A2, B0, B2, C0, C1, C2, M0, M1, M2, R1, R2, sigp, &
           alpha1, alpha3, C0t, C2t, Atld, Rt  , kmask, &
           f_dmask, f_kmask 

          contains

          subroutine set_crystal(aANG, cANG, T2frac, f0)
          ! Convert input lattice constant to atomic units 
          ! Set up the NN and NNN vectors and matrices 
          ! Allocate momentum grid 
          ! Input: 
          ! aANG    : lattice constant in Angstroms, from main program
          ! T2frac  : dephasing time in fractions of the laser period, i.e. T2  = T2frac * T0 

            implicit none
            real*8 :: T2frac, aANG, cANG, kstart, kend, kx, ky 
            real*8 :: kxstart, kystart, kxend, kyend
            integer :: i, j 
            real*8 :: eV2au, ANG2au
            real*8 :: f0 
            real*8 :: kbuff, kabs, kpar

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
            alpha1 = 0.99 
            alpha3 = -0.15;
            C0t    = C0 + alpha3 * M0;
            C2t    = C2 + alpha3 * M2;
            Atld   = A0 *  alpha1    ; 
            Rt     = 0.5d0 * R1 * alpha1 ;   
            
            Tdephas_2 = Tper*T2frac        
            alat = aANG*1.d-10/aBohr            

        
            nkpts = nkx * nky                     ! Total number of k-points in 2D 
            nmax = 2*nkpts                        
          

            ! Allocate arrays
            if(.not.allocated(kx1d)) allocate(kx1d(nkx))
            if(.not.allocated(ky1d)) allocate(ky1d(nky))
            if(.not.allocated(kxy_grid)) allocate(kxy_grid(nkpts,2))
            if(.not.allocated(kmask)) allocate(kmask(nkpts))
        

            kystart  = -4.d0*pi/alat/3.d0 * bzcut
            kyend    =  4.d0*pi/alat/3.d0 * bzcut
            kystart  = -2.d0*pi/alat/sqrt(3.d0) * bzcut
            kyend    =  2.d0*pi/alat/sqrt(3.d0) * bzcut
            kxstart  = -2.d0*pi/alat/sqrt(3.d0) * bzcut
            kxend    =  2.d0*pi/alat/sqrt(3.d0) * bzcut


            dkx = (kxend-kxstart)/dble(nkx-1)

            ! Initialize 1D grids 
            do i = 1, nkx 
               kx1d(i) = dble(i-1)*dkx + kxstart 
            enddo

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
            !kabs = 0.85d0*kyend
            !kbuff = 0.15d0*kyend


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
                  !kmask( (i-1)*nky + j ) = 1.d0 + erf( -(kpar-kabs)**2/(kbuff**2) )
                else
                  kmask( (i-1)*nky + j ) = 0.d0
                endif

              enddo
            enddo                    

            f0glob = f0  

            return
          end subroutine set_crystal

          real*8 function E2Dsurf(kxy, m)
          ! Eigenenergy of the TSS VB/CB
          ! kxy :: 2-vector (kx, ky)
          ! m: = -1 --> VB; = +1 --> CB 

            implicit none
            real*8, dimension(2) ::  kxy
            integer              ::  m 
            real*8               ::  kx, ky, kpar
            real*8               ::  sgn 

            sgn  = dble(m)
            kx   = kxy(1)
            ky   = kxy(2)            
            kpar = sqrt( kx**2 + ky**2 )

            E2DSurf = C0t + C2t*(kx**2 + ky**2) + &
            Sqrt(Atld**2*(kx**2 + ky**2) + 4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2)*sgn
            return

          end function E2DSurf

          function vGroup(kxy, m, f0, sigma) result(vgr)
          ! Momentum derivative of band eigenenergy; group velocity
            implicit none
            real*8, dimension(2) ::   kxy, vgr
            integer              ::   m 
            real*8               ::   sgn, kx, ky, kpar 
            real*8               ::   vgr_x, vgr_y
            real*8               :: f0, eps, s 
            real*8, optional     :: sigma 
            kx = kxy(1)
            ky = kxy(2)
            
            kpar = sqrt( kx**2 + ky**2 )

            sgn  = dble(m)
            kpar = sqrt( kx**2 + ky**2 )

            s = 0.d0 
            if(present(sigma)) then 
              if(sigma .gt. 1.d-5) then 
                s = 1.d0/( (sigma * 4.d0*pi/3.d0/alat)**2 )
              else
                s = 0.d0 
              endif
            endif

            eps = (f0 * Atld)**2 * exp(- s * kpar**2)

            vgr_x = 2.d0*C2t*kx + (kx*(Atld**2 + 12.d0*(kx**4 - 4.d0*kx**2*ky**2 + 3.d0*ky**4)*Rt**2)*sgn)/&
             ( Sqrt(eps +Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))

            vgr_y = 2.d0*C2t*ky + (ky*(Atld**2 - 24.d0*kx**2*(kx**2 - 3.d0*ky**2)*Rt**2)*sgn)/&
             ( Sqrt(eps +Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))

            vgr = (/ vgr_x, vgr_y /)
            return

          end function vGroup

          real*8 function psiNorm(kxy, m, f0, sigma)
          ! Norm of the surface spinor
            implicit none
            real*8, dimension(2) :: kxy
            real*8               :: kx, ky, kpar
            integer              :: m
            real*8               :: sgn, NNorm
            real*8               :: f0, eps 
            real*8, optional     :: sigma 
            real*8               :: s 

            kx = kxy(1)
            ky = kxy(2)
            

            kpar = sqrt( kx**2 + ky**2 )

            s = 0.d0 
            if(present(sigma)) then
              if( sigma .gt. 1.d-5) then 
                s = 1.d0/( ( sigma*4.d0*pi/alat/3.d0 )**2 )
              else
                s = 0.d0
              endif
            endif


            eps = (f0 * A0)**2 * exp(-kpar**2 * s)

            kpar = sqrt( kx**2 + ky**2 )
            sgn = dble(m)

            NNorm = Atld**2*(kx**2 + ky**2) + &
             (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*&
                 sgn)**2
            psiNorm = 1.d0/( sqrt(eps +NNorm))
            return
          end function psiNorm

          function psi2DSurf(kxy, m, normflag, f0, sigma ) result(psi)
          ! Normalized (normflag == TRUE) or not normalized (normflag == FALSE) spinor
          ! m   = {-1, +1} --> VB / CB 
          ! normflag :: Boolean, whether to retrieve the fully normalized spinor (TRUE) or just the vector as obtained from the Hamiltonian diagonalization with Mathematica (FALSE)
            implicit none
            real*8, dimension(2)     :: kxy
            complex*16, dimension(2) :: psi 
            logical                  :: normflag
            real*8                   :: kx, ky, kpar 
            integer                  :: m 
            real*8                   :: sgn 
            complex*16               :: psi1, psi2
            real*8                   :: NN 
            real*8                   :: f0
            real*8, optional         :: sigma 

            kx = kxy(1)
            ky = kxy(2)            

            kpar = sqrt( kx**2 + ky**2 )

            sgn = dble(m)

            psi1 = 2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
            Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn
            psi2 = Atld * (-iimag * kx + ky)
            if(present(sigma)) then
              NN = psiNorm(kxy, m, f0, sigma )
            else 
              NN = psiNorm(kxy, m, f0 )
            endif
            psi = (/ psi1, psi2 /)
            if( normflag ) psi = NN * psi 

            return 
          end function psi2DSurf

          function Dpsi2DSurf_dk(kxy,  m, f0, sigma) result(dpsi_dk)
            ! 2D TSS wavefunction gradient 
            ! arguments same as above :
            ! kxy : k-vector in 2D 
            ! m   = {-1, +1} --> VB / CB 
            ! f0: scaling factor for the regularization parameter eps = (f0 A0)**2
            implicit none
            real*8, dimension(2)       :: kxy
            integer                    :: m 
            complex*16, dimension(2,2) :: dpsi_dk 
            complex*16, dimension(2)   :: dpsi_dkx, dpsi_dky
            real*8                     :: sgn 
            complex*16                 :: dpsi_dkx1, dpsi_dkx2
            complex*16                 :: dpsi_dky1, dpsi_dky2
            complex*16, dimension(2)   :: psi 
            complex*16                 :: dN_dkx, dN_dky
            real*8                     :: NN 
            real*8                     :: f0, eps 
            real*8                     :: kx, ky
            real*8, optional           :: sigma 
            real*8                     :: s 

            kx = kxy(1)
            ky = kxy(2)
            

            s = 0.d0 
            if(present(sigma)) then
              if( sigma .gt. 1.d-5) then 
                s = 1.d0/( ( sigma*4.d0*pi/alat/3.d0 )**2 )
                NN = psiNorm(kxy, m, f0, sigma)
                psi = psi2DSurf(kxy, m, .FALSE.,f0,sigma)
              else
                s = 0.d0
                NN = psiNorm(kxy, m, f0)
                psi = psi2DSurf(kxy, m, .FALSE.,f0)
              endif
            else
              NN = psiNorm(kxy, m, f0)
              psi = psi2DSurf(kxy, m, .FALSE.,f0)
            endif


            eps = (f0 * Atld)**2 * exp(-(kx**2 + ky**2)* s)

            sgn = dble(m) 

            

            dN_dkx = -((2.d0*Atld**2*kx + 2.d0*(2*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)*&
                (6.d0*kx**2*Rt - 6.d0*ky**2*Rt + (kx*(Atld**2 + 12.d0*(kx**4 - &
                4.d0*kx**2*ky**2 + 3.d0*ky**4)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))))/&
                ( 2.d0*(eps +Atld**2*(kx**2 + ky**2) + &
                (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)**2)**1.5))

            dN_dky = -((2.d0*Atld**2*ky + 2.d0*(2*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)*&
                (-12.d0*kx*ky*Rt + (ky*(Atld**2 - 24.d0*kx**2*(kx**2 - 3.d0*ky**2)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))))/&
                ( 2.d0*(eps +Atld**2*(kx**2 + ky**2) + &
                (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)**2)**1.5))

            dpsi_dkx1 =  6.d0*kx**2*Rt - 6.d0*ky**2*Rt + &
                (kx*(Atld**2 + 12.d0*(kx**4 - 4.d0*kx**2*ky**2 + 3.d0*ky**4)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))
            dpsi_dkx2 = -iimag * Atld

            dpsi_dky1 = -12.d0*kx*ky*Rt + (ky*(Atld**2 - 24.d0*kx**2*(kx**2 - 3.d0*ky**2)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))
            dpsi_dky2 =  Atld
            

            dpsi_dkx = dN_dkx * psi + NN * (/ dpsi_dkx1, dpsi_dkx2 /)
            dpsi_dky = dN_dky * psi + NN * (/ dpsi_dky1, dpsi_dky2 /)


            dpsi_dk(:,1) = dpsi_dkx
            dpsi_dk(:,2) = dpsi_dky            

            return
          end function Dpsi2DSurf_dk

          function DvcS(kxy, f0, sigma) result(dipole)
            ! dipole moment between VB and CB (interband)
            ! kxy : k-vector in 2D            
            ! f0 : scaling factor for the regularization parameter eps = (f0 A0)**2 to treat discontinuity at Gamma point
            implicit none
            real*8, dimension(2)       :: kxy
            complex*16, dimension(2)   :: dipole 
            complex*16, dimension(2)   :: psiv
            complex*16, dimension(2,2) :: dpsic_dk
            real*8                     :: f0
            real*8, optional           :: sigma 
            
            if(present(sigma)) then
              psiv = psi2DSurf(kxy, -1, .TRUE., f0, sigma)
              dpsic_dk = Dpsi2DSurf_dk(kxy, 1, f0, sigma)
            else
              psiv = psi2DSurf(kxy, -1, .TRUE., f0)
              dpsic_dk = Dpsi2DSurf_dk(kxy, 1, f0)
            endif
            dipole(1) = iimag * dot_product(psiv, dpsic_dk(:,1))
            dipole(2) = iimag * dot_product(psiv, dpsic_dk(:,2))          
            return 
          end function DvcS

          function xiS(kxy, m, f0,sigma) result (xi)
            ! Berry connection of VB or CB 
            ! kxy  : k-vector in 2D             
            ! m    : {-1,1} - index for VB (-1) or CB (+1)
            ! f0   : regularization parameter scaling factor for treating the singularity at the Gamma point for the CB (M0<0) or the VB (M0>0)
            implicit none
            real*8, dimension(2)       :: kxy
            integer                    :: m 
            complex*16, dimension(2)   :: xi 
            complex*16, dimension(2)   :: psi
            complex*16, dimension(2,2) :: dpsi_dk
            real*8                     :: f0
            real*8, optional           :: sigma 
            
          
            if(present(sigma)) then 
              psi     = psi2DSurf(kxy, m, .TRUE.,f0, sigma )
              dpsi_dk = Dpsi2DSurf_dk(kxy, m, f0, sigma)
            else
              psi     = psi2DSurf(kxy, m, .TRUE.,f0 )
              dpsi_dk = Dpsi2DSurf_dk(kxy, m, f0)
            endif
            xi(1) = iimag * dot_product(psi, dpsi_dk(:,1))
            xi(2) = iimag * dot_product(psi, dpsi_dk(:,2))

            return 
          end function xiS

          function berryS(kxy, m, f0, sigma) result (Omega)
            ! Berry curvature for VB/CB 
            ! kxy : k-vector in 2D             
            ! m    : {-1,1} - index for VB (-1) or CB (+1)
            ! f0   : regularization parameter scaling factor 

            implicit none 
            real*8, dimension(2)       :: kxy
            integer                    :: cmp, m 
            complex*16                 :: Omega          
            real*8                     :: f0, s, eps 
            real*8, optional           :: sigma 
            real*8                     :: kx, ky, sgn 

            kx = kxy(1)
            ky = kxy(2)
            

            s = 0.d0 
            if(present(sigma)) then
              if( sigma .gt. 1.d-5) then 
                s = 1.d0/( ( sigma*4.d0*pi/alat/3.d0 )**2 )
                
              else
                s = 0.d0
                
              endif
            
            endif




            eps = (f0 * Atld)**2 * exp(-(kx**2 + ky**2)* s)


            sgn = dble(m) 

              
            Omega = (A0**2*kx*(kx**2 - 3.d0*ky**2)*R1*sgn)/&
            (eps + A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2)**1.5

            return

          end function berryS


          ! Custom implementation of the cross product of two vectors lying in the (xy)-plane, the result is aligned along the z-axis

          real*8 function zcross(v1, v2)
            implicit none
            real*8, dimension(2) :: v1, v2

            zcross = v1(1)*v2(2) - v2(1)*v1(2)
            return
          end function zcross

          subroutine get_crystal_kxy(kxy, egap, xigap, dmat, f0, sigma)
            implicit none
            real*8, dimension(2)       :: kxy
            real*8                     :: egap             
            complex*16, dimension(2)   :: dmat, xigap 
            real*8                     :: f0
            real*8                     :: kx, ky, kpar
            integer                    :: m
            real*8                     :: sgn, NNorm
            real*8                     :: eps
            complex*16, dimension(2)   :: psi 
            complex*16                 :: psi1, psi2
            real*8                     :: NN 
            complex*16, dimension(2,2) :: dpsi_dk 
            complex*16, dimension(2)   :: dpsi_dkx, dpsi_dky
            complex*16                 :: dpsi_dkx1, dpsi_dkx2
            complex*16                 :: dpsi_dky1, dpsi_dky2
            complex*16                 :: dN_dkx, dN_dky

            real*8                     :: eVB, eCB, eband
            complex*16, dimension(2)   :: psiv, psic 
            complex*16, dimension(2,2) :: dpsiv_dk, dpsic_dk
            real*8                     :: s 
            real*8, optional           :: sigma 

            kx = kxy(1)
            ky = kxy(2)
            

            kpar = sqrt( kx**2 + ky**2 )
            
            s = 0.d0 
            if(present(sigma)) then
              if( sigma .gt. 1.d-5) then 
                s = 1.d0/( ( sigma*4.d0*pi/alat/3.d0 )**2 )
                
              else
                s = 0.d0
                
              endif
            
            endif
            eps = (f0 * Atld)**2 * exp( - (kx**2 + ky**2) * s )

            do m = -1, 1, 2 

              sgn = dble(m)
            
              eband = C0t + C2t*(kx**2 + ky**2) + &
               Sqrt(Atld**2*(kx**2 + ky**2) + 4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2)*sgn

              ! Caluclate norm of the wavefunction 
              NNorm = Atld**2*(kx**2 + ky**2) + &
                (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + Sqrt(Atld**2*(kx**2 + ky**2) + &
                4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)**2

              NN = 1.d0/( sqrt(eps +NNorm))

              ! Caluclate  wavefunction 

              psi1 = 2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                 Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - &
                 3.d0*ky**2)**2*Rt**2)*sgn
              psi2 = Atld * (-iimag * kx + ky)
            
              psi = (/ psi1, psi2 /)

              ! Caluclate gradient of the wavefunction 

              dN_dkx = -((2.d0*Atld**2*kx + 2.d0*(2*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)*&
                (6.d0*kx**2*Rt - 6.d0*ky**2*Rt + (kx*(Atld**2 + 12.d0*(kx**4 - &
                4.d0*kx**2*ky**2 + 3.d0*ky**4)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))))/&
                ( 2.d0*(eps +Atld**2*(kx**2 + ky**2) + &
                (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)**2)**1.5))

              dN_dky = -((2.d0*Atld**2*ky + 2.d0*(2*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)*&
                (-12.d0*kx*ky*Rt + (ky*(Atld**2 - 24.d0*kx**2*(kx**2 - 3.d0*ky**2)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))))/&
                (2.d0*(eps + Atld**2*(kx**2 + ky**2) + &
                (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)**2)**1.5))

              dpsi_dkx1 =  6.d0*kx**2*Rt - 6.d0*ky**2*Rt + &
                (kx*(Atld**2 + 12.d0*(kx**4 - 4.d0*kx**2*ky**2 + 3.d0*ky**4)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))
              dpsi_dkx2 = -iimag * Atld

              dpsi_dky1 = -12.d0*kx*ky*Rt + (ky*(Atld**2 - 24.d0*kx**2*(kx**2 - 3.d0*ky**2)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))
              dpsi_dky2 =  Atld
            

              dpsi_dkx = dN_dkx * psi + NN * (/ dpsi_dkx1, dpsi_dkx2 /)
              dpsi_dky = dN_dky * psi + NN * (/ dpsi_dky1, dpsi_dky2 /)
  
  
              dpsi_dk(:,1) = dpsi_dkx
              dpsi_dk(:,2) = dpsi_dky   

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

          end subroutine get_crystal_kxy

          subroutine get_crystal_kxy_compact(kxy, egap, xigap, dmat, f0, sigma)
            implicit none
            real*8, dimension(2)       :: kxy
            real*8                     :: egap             
            complex*16, dimension(2)   :: dmat, xigap 
            real*8                     :: f0
            real*8                     :: kx, ky, kpar
            integer                    :: m
            real*8                     :: sgn
            real*8                     :: eps            
            real*8                     :: s 
            real*8, optional           :: sigma 
            real*8                     :: eps1, eps2, eps3 
            real*8                     :: krot, throt, dipolemask, kbuff, kabs, kyend

            kx = kxy(1)
            ky = kxy(2)
            

            kpar = sqrt( kx**2 + ky**2 )
            
            s = 0.d0 
            if(present(sigma)) then
              if( sigma .gt. 1.d-5) then 
                s = 1.d0/( ( sigma*4.d0*pi/alat/3.d0 )**2 )
                
              else
                s = 0.d0
                
              endif
            
            endif
            !eps = (f0 * Atld)**2 * exp( - (kx**2 + ky**2) * s )

            krot = kx 
            eps1 = exp( - krot**2 * s )
            krot = cos(pi/3.d0)*kx - sin(pi/3.d0)*ky 
            eps2 = exp( - krot**2 * s )
            krot = cos(-pi/3.d0)*kx - sin(-pi/3.d0)*ky 
            eps3 = exp( -krot**2 * s )
            eps = (f0 * Atld)**2 * (eps1+eps2+eps3)/3.d0 
            eps = (f0 * Atld)**2 * exp( -(kx**2 + ky**2) * s )

            !kyend = ky1d(nky)
            !kbuff = f_dmask*kyend 
            !kabs  = kyend - kbuff 
!
!            !kbuff = 0.3d0 * kyend
!            !kabs  = 0.1d0 * kyend 
!
!            !if( kpar .lt. kabs ) then 
!            !  dipolemask = 1.d0!*(eps1+eps2+eps3)/3.d0 
!            !else if( kpar .lt. kyend ) then
!            !  !dipolemask = abs(cos(abs(kabs-kpar) * pi / (2.d0*kbuff)))**2!0.125!*(eps1+eps2+eps3)/3.d0 
!            !  dipolemask = 1.d0 + erf( -(kpar-kabs)**2/(kbuff**2) )
!            !else
!            !  dipolemask = 0.d0
            !endif

            xigap(1) = (-2.d0*ky*(kx**3 - 3.d0*kx*ky**2)*Rt)/&
               (eps + (kx**2 + ky**2)*Sqrt(Atld**2*(kx**2 + ky**2) + &
                4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))

            xigap(2) = (2.d0*kx**2*(kx**2 - 3.d0*ky**2)*Rt)/&
               (eps + (kx**2 + ky**2)*Sqrt(Atld**2*(kx**2 + ky**2) + &
                4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))

            dmat(1) = (Atld**2*(kx - dcmplx(0.d0,1.d0)*ky)*&
              (4.d0*kx*(kx**2 - 3.d0*ky**2)*&
              (dcmplx(0.d0,-2.d0)*kx**3 + 3.d0*kx**2*ky - 3.d0*ky**3)*Rt**2*&
              (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)) + &
              Atld**2*(kx + dcmplx(0.d0,1.d0)*ky)*&
              (dcmplx(0.d0,-4.d0)*kx**4*Rt + 2.d0*kx**3*ky*Rt - &
              dcmplx(0.d0,6.d0)*kx**2*ky**2*Rt - 6.d0*kx*ky**3*Rt + &
              dcmplx(0.d0,6.d0)*ky**4*Rt + &
              ky*Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))))/&
              (eps + 2.d0*(Atld**2*(kx**2 + ky**2) + &
              2.d0*kx*(kx**2 - 3.d0*ky**2)*Rt*&
              (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)))**1.5*&
              Sqrt((Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*&
              (Atld**2*(kx**2 + ky**2) - &
              2.d0*kx*(kx**2 - 3.d0*ky**2)*Rt*&
              (-2.d0*kx**3*Rt + 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)))))

            dmat(2)=-(Atld**2*kx*(kx - dcmplx(0.d0,1.d0)*ky)*&
              (4.d0*kx*(kx**2 - 3.d0*ky**2)*&
              (kx**2 - dcmplx(0.d0,6.d0)*kx*ky + 3.d0*ky**2)*Rt**2*&
              (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)) + &
              Atld**2*(kx + dcmplx(0.d0,1.d0)*ky)*&
              (2.d0*(kx - dcmplx(0.d0,1.d0)*ky)*&
              (kx**2 - dcmplx(0.d0,6.d0)*kx*ky + 3.d0*ky**2)*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))))/&
              (eps + 2.d0*(Atld**2*(kx**2 + ky**2) + &
              2.d0*kx*(kx**2 - 3.d0*ky**2)*Rt*&
              (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)))**1.5*&
              Sqrt((Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*&
              (Atld**2*(kx**2 + ky**2) - &
              2.d0*kx*(kx**2 - 3.d0*ky**2)*Rt*&
              (-2.d0*kx**3*Rt + 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)))))

            !dmat = dipolemask * dmat! * (eps1+eps2+eps3)/3.d0 
            !dmat = dcmplx(dble(dmat), dipolemask * aimag(dmat))

            egap = 2.d0*Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)

            return
                        

          end subroutine get_crystal_kxy_compact

          subroutine get_crystal_velocities_kxy_compact(kxy, vgr_v, vgr_c, omg_v, omg_c, dmat, f0, sigma)
            implicit none
            real*8, dimension(2)       :: kxy            
            complex*16, dimension(2)   :: dmat
            complex*16                 :: omg_v, omg_c 
            real*8, dimension(2)       :: vgr_v, vgr_c 
            real*8                     :: kx, ky, kpar
            real*8                     :: dE_dkx, dE_dky 
            integer                    :: m
            real*8                     :: sgn
            real*8                     :: f0, eps
            real*8                     :: s 
            real*8, optional           :: sigma 
            real*8                     :: krot, eps1, eps2, eps3
            real*8                     :: kbuff, kabs, dipolemask  , kyend

            kx = kxy(1)
            ky = kxy(2)

            kpar = sqrt( kx**2 + ky**2 )          
            
            s = 0.d0 
            if(present(sigma)) then
              if( sigma .gt. 1.d-5) then 
                s = 1.d0/( ( sigma*4.d0*pi/alat/3.d0 )**2 )
                
              else
                s = 0.d0
                
              endif
            
            endif
            !eps = (f0 * Atld)**2 * exp( -(kx**2 + ky**2) * s )
            krot = kx 
            eps1 = exp( - krot**2 * s )
            krot = cos(pi/3.d0)*kx - sin(pi/3.d0)*ky 
            eps2 = exp( - krot**2 * s )
            krot = cos(-pi/3.d0)*kx - sin(-pi/3.d0)*ky 
            eps3 = exp( -krot**2 * s )
            eps = (f0 * Atld)**2 * (eps1+eps2+eps3)/3.d0 
            eps = (f0 * Atld)**2 * exp( -(kx**2 + ky**2) * s )

            !kyend = ky1d(nky)
            !kbuff = f_dmask*kyend 
            !kabs  = kyend - kbuff 
!
!            !kbuff = 0.3d0 * kyend
            !kabs  = 0.1d0 * kyend

            !if( kpar .lt. kabs ) then 
            !  dipolemask = 1.d0 !* (eps1+eps2+eps3)/3.d0 
            !else if( kpar .lt. kyend ) then
            !  !dipolemask = abs(cos(abs(kabs-kpar) * pi / (2.d0*kbuff)))**2!0.125!*(eps1+eps2+eps3)/3.d0 
            !  dipolemask = 1.d0 + erf( -(kpar-kabs)**2/(kbuff**2) )
            !else
            !  dipolemask = 0.d0
            !endif


            do m = -1, 1, 2 

              sgn = dble(m)              

              dE_dkx = 2.d0*C2t*kx + (kx*(Atld**2 + 12.d0*(kx**4 - 4.d0*kx**2*ky**2 + 3.d0*ky**4)*Rt**2)*sgn)/&
             ( Sqrt(eps*1.d0 +Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))

              dE_dky = 2.d0*C2t*ky + (ky*(Atld**2 - 24.d0*kx**2*(kx**2 - 3.d0*ky**2)*Rt**2)*sgn)/&
             ( Sqrt(eps*1.d0 +Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))

                        

              if( m .eq. -1) then 
                
                vgr_v = (/ dE_dkx, dE_dky /)
               
              else
                
                vgr_c = (/ dE_dkx, dE_dky /)
                
              endif

            enddo   

            omg_v =    (-2.d0*Atld**2*kx*(kx**2 - 3.d0*ky**2)*Rt)/&
              (eps + Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)**1.5

            omg_c =       (2.d0*Atld**2*kx*(kx**2 - 3.d0*ky**2)*Rt)/&
              (eps + Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)**1.5         

            dmat(1) = (Atld**2*(kx - dcmplx(0.d0,1.d0)*ky)*&
              (4.d0*kx*(kx**2 - 3.d0*ky**2)*&
              (dcmplx(0.d0,-2.d0)*kx**3 + 3.d0*kx**2*ky - 3.d0*ky**3)*Rt**2*&
              (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)) + &
              Atld**2*(kx + dcmplx(0.d0,1.d0)*ky)*&
              (dcmplx(0.d0,-4.d0)*kx**4*Rt + 2.d0*kx**3*ky*Rt - &
              dcmplx(0.d0,6.d0)*kx**2*ky**2*Rt - 6.d0*kx*ky**3*Rt + &
              dcmplx(0.d0,6.d0)*ky**4*Rt + &
              ky*Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))))/&
              (eps + 2.d0*(Atld**2*(kx**2 + ky**2) + &
              2.d0*kx*(kx**2 - 3.d0*ky**2)*Rt*&
              (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)))**1.5*&
              Sqrt((Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*&
              (Atld**2*(kx**2 + ky**2) - &
              2.d0*kx*(kx**2 - 3.d0*ky**2)*Rt*&
              (-2.d0*kx**3*Rt + 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)))))

            dmat(2)=-(Atld**2*kx*(kx - dcmplx(0.d0,1.d0)*ky)*&
              (4.d0*kx*(kx**2 - 3.d0*ky**2)*&
              (kx**2 - dcmplx(0.d0,6.d0)*kx*ky + 3.d0*ky**2)*Rt**2*&
              (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)) + &
              Atld**2*(kx + dcmplx(0.d0,1.d0)*ky)*&
              (2.d0*(kx - dcmplx(0.d0,1.d0)*ky)*&
              (kx**2 - dcmplx(0.d0,6.d0)*kx*ky + 3.d0*ky**2)*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))))/&
              (eps + 2.d0*(Atld**2*(kx**2 + ky**2) + &
              2.d0*kx*(kx**2 - 3.d0*ky**2)*Rt*&
              (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)))**1.5*&
              Sqrt((Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*&
              (Atld**2*(kx**2 + ky**2) - &
              2.d0*kx*(kx**2 - 3.d0*ky**2)*Rt*&
              (-2.d0*kx**3*Rt + 6.d0*kx*ky**2*Rt + &
              Sqrt(Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)))))    

            !dmat = dipolemask * dmat !* (eps1+eps2+eps3)/3.d0     
            !dmat = dcmplx(dble(dmat), dipolemask * aimag(dmat))

              return

          end subroutine get_crystal_velocities_kxy_compact


          subroutine get_crystal_velocities_kxy(kxy, vgr_v, vgr_c, omg_v, omg_c, dmat, f0, sigma)
            implicit none
            real*8, dimension(2)       :: kxy            
            complex*16, dimension(2)   :: dmat
            complex*16                 :: omg_v, omg_c 
            real*8, dimension(2)       :: vgr_v, vgr_c 
            real*8                     :: kx, ky, kpar
            real*8                     :: dE_dkx, dE_dky 
            integer                    :: m
            real*8                     :: sgn, NNorm
            real*8                     :: f0, eps
            complex*16, dimension(2)   :: psi 
            complex*16                 :: psi1, psi2
            real*8                     :: NN 
            complex*16, dimension(2,2) :: dpsi_dk 
            complex*16, dimension(2)   :: dpsi_dkx, dpsi_dky
            complex*16                 :: dpsi_dkx1, dpsi_dkx2
            complex*16                 :: dpsi_dky1, dpsi_dky2
            complex*16                 :: dN_dkx, dN_dky         
            complex*16, dimension(2)   :: psiv, psic 
            complex*16, dimension(2,2) :: dpsiv_dk, dpsic_dk
            real*8                     :: s 
            real*8, optional           :: sigma 

            kx = kxy(1)
            ky = kxy(2)

            kpar = sqrt( kx**2 + ky**2 )          
            
            s = 0.d0 
            if(present(sigma)) then
              if( sigma .gt. 1.d-5) then 
                s = 1.d0/( ( sigma*4.d0*pi/alat/3.d0 )**2 )
                
              else
                s = 0.d0
                
              endif
            
            endif
            eps = (f0 * Atld)**2 * exp( -(kx**2 + ky**2) * s )


            do m = -1, 1, 2 

              sgn = dble(m)              

              dE_dkx = 2.d0*C2t*kx + (kx*(Atld**2 + 12.d0*(kx**4 - 4.d0*kx**2*ky**2 + 3.d0*ky**4)*Rt**2)*sgn)/&
             ( Sqrt(eps*1.d0 +Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))

              dE_dky = 2.d0*C2t*ky + (ky*(Atld**2 - 24.d0*kx**2*(kx**2 - 3.d0*ky**2)*Rt**2)*sgn)/&
             ( Sqrt(eps*1.d0 +Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))

               ! Caluclate norm of the wavefunction 
              NNorm = Atld**2*(kx**2 + ky**2) + &
                (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + Sqrt(Atld**2*(kx**2 + ky**2) + &
                4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)**2

              NN = 1.d0/( sqrt(eps +NNorm))

              ! Caluclate  wavefunction 

              psi1 = 2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                 Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - &
                 3.d0*ky**2)**2*Rt**2)*sgn
              psi2 = Atld * (-iimag * kx + ky)
            
              psi = (/ psi1, psi2 /)

              ! Caluclate gradient of the wavefunction 

              dN_dkx = -((2.d0*Atld**2*kx + 2.d0*(2*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)*&
                (6.d0*kx**2*Rt - 6.d0*ky**2*Rt + (kx*(Atld**2 + 12.d0*(kx**4 - &
                4.d0*kx**2*ky**2 + 3.d0*ky**4)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))))/&
                (2.d0*(eps + Atld**2*(kx**2 + ky**2) + &
                (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)**2)**1.5))

              dN_dky = -((2.d0*Atld**2*ky + 2.d0*(2*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)*&
                (-12.d0*kx*ky*Rt + (ky*(Atld**2 - 24.d0*kx**2*(kx**2 - 3.d0*ky**2)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))))/&
                ( 2.d0*(eps +Atld**2*(kx**2 + ky**2) + &
                (2.d0*kx**3*Rt - 6.d0*kx*ky**2*Rt + &
                Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)*sgn)**2)**1.5))

              dpsi_dkx1 =  6.d0*kx**2*Rt - 6.d0*ky**2*Rt + &
                (kx*(Atld**2 + 12.d0*(kx**4 - 4.d0*kx**2*ky**2 + 3.d0*ky**4)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))
              dpsi_dkx2 = -iimag * Atld

              dpsi_dky1 = -12.d0*kx*ky*Rt + (ky*(Atld**2 - 24.d0*kx**2*(kx**2 - 3.d0*ky**2)*Rt**2)*sgn)/&
                (eps*0.d0 + Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2))
              dpsi_dky2 =  Atld
            

              dpsi_dkx = dN_dkx * psi + NN * (/ dpsi_dkx1, dpsi_dkx2 /)
              dpsi_dky = dN_dky * psi + NN * (/ dpsi_dky1, dpsi_dky2 /)
  
  
              dpsi_dk(:,1) = dpsi_dkx
              dpsi_dk(:,2) = dpsi_dky             

              if( m .eq. -1) then 
                psiv = NN * psi
                dpsiv_dk = dpsi_dk 
                vgr_v = (/ dE_dkx, dE_dky /)
                omg_v = (A0**2*kx*(kx**2 - 3.d0*ky**2)*R1*sgn)/&
                 (eps + A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2)**1.5
              else
                psic = NN * psi
                dpsic_dk = dpsi_dk
                vgr_c = (/ dE_dkx, dE_dky /)
                omg_c = (A0**2*kx*(kx**2 - 3.d0*ky**2)*R1*sgn)/&
                 (eps + A0**2*(kx**2 + ky**2) + kx**2*(kx**2 - 3.d0*ky**2)**2*R1**2)**1.5
              endif

            enddo            

            dmat(1) = iimag * dot_product(psiv, dpsic_dk(:,1))
            dmat(2) = iimag * dot_product(psiv, dpsic_dk(:,2))            

          end subroutine get_crystal_velocities_kxy
            


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
            real*8    , allocatable, dimension(:,:)   :: Esurf2D_VB, Esurf2D_CB
            real*8    , allocatable, dimension(:,:,:) :: vGroupS_VB, vGroupS_CB
            real*8    , allocatable, dimension(:,:,:) :: xiSurf2D_VB, xiSurf2D_CB
            real*8    , allocatable, dimension(:,:)   :: berrySS_VB, berrySS_CB
            complex*16, allocatable, dimension(:,:,:) :: dipoleSSVC
            real*8    , allocatable, dimension(:,:)   :: psiNorm_VB, psiNorm_CB
            integer :: i, j
            real*8, dimension(2) :: kxy
            real*8 :: kx, ky
            real*8 :: sigma, f0 

            f0 = f0glob

            sigma = sigp

            if(.not.allocated(Esurf2D_VB)) allocate(Esurf2D_VB(nkx, nky))
            if(.not.allocated(Esurf2D_CB)) allocate(Esurf2D_CB(nkx, nky))
            if(.not.allocated(vGroupS_VB)) allocate(vGroupS_VB(nkx, nky, 2))
            if(.not.allocated(vGroupS_CB)) allocate(vGroupS_CB(nkx, nky, 2))
            if(.not.allocated(xiSurf2D_VB)) allocate(xiSurf2D_VB(nkx, nky, 2))
            if(.not.allocated(xiSurf2D_CB)) allocate(xiSurf2D_CB(nkx, nky, 2))
            if(.not.allocated(berrySS_VB)) allocate(berrySS_VB(nkx, nky))
            if(.not.allocated(berrySS_CB)) allocate(berrySS_CB(nkx, nky))
            if(.not.allocated(dipoleSSVC)) allocate(dipoleSSVC(nkx, nky, 2))
            if(.not.allocated(psiNorm_VB)) allocate(psiNorm_VB(nkx,nky))
            if(.not.allocated(psiNorm_CB)) allocate(psiNorm_CB(nkx,nky))
            
            do j = 1, nky 
              ky = ky1d(j)
              do i = 1, nkx 
                kx = kx1d(i)
                kxy = (/ kx, ky /)

                Esurf2D_VB(i,j) = E2Dsurf(kxy, -1)
                Esurf2D_CB(i,j) = E2Dsurf(kxy,  1)

                vGroupS_VB(i,j,:) = vGroup(kxy, -1, f0, sigma)
                vGroupS_CB(i,j,:) = vGroup(kxy,  1, f0, sigma)

                psiNorm_VB(i,j) = psiNorm(kxy, -1, f0, sigma)
                psiNorm_CB(i,j) = psiNorm(kxy,  1, f0, sigma)

              
                xiSurf2D_VB(i,j,:) = dble(xiS( kxy,  -1, f0, sigma ))
                xiSurf2D_CB(i,j,:) = dble(xiS( kxy,   1, f0, sigma ))                
                
                berrySS_VB(i,j) = dble(berryS(kxy, -1, f0, sigma ))
                berrySS_CB(i,j) = dble(berryS(kxy,  1, f0, sigma ))

                dipoleSSVC(i,j,:) = DvcS( kxy, f0, sigma )                

              enddo
            enddo

            open(351,file=Opath(1:Olast)//'Esurf2D_VB.dat', status = "replace")
            open(353,file=Opath(1:Olast)//'Esurf2D_CB.dat', status = "replace")          

            open(355,file=Opath(1:Olast)//'vGroupS_VB.dat', status = "replace")
            open(357,file=Opath(1:Olast)//'vGroupS_CB.dat', status = "replace")

            open(371,file=Opath(1:Olast)//'xiSurf2D_VB.dat', status = "replace")
            open(373,file=Opath(1:Olast)//'xiSurf2D_CB.dat', status = "replace")            

            open(391,file=Opath(1:Olast)//'berrySS_VB.dat', status = "replace")
            open(393,file=Opath(1:Olast)//'berrySS_CB.dat', status = "replace")            

            open(511,file=Opath(1:Olast)//'dipoleSSVC.dat', status = "replace")

            open(711,file=Opath(1:Olast)//'psiNorm_VB.dat', status = "replace")
            open(713,file=Opath(1:Olast)//'psiNorm_CB.dat', status = "replace")  

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

                write(351, 5002) kx, ky, Esurf2D_VB(i,j)
                write(353, 5002) kx, ky, Esurf2D_CB(i,j)

                write(355, 5005) kx, ky, vGroupS_VB(i,j,1), vGroupS_VB(i,j,2)
                write(357, 5005) kx, ky, vGroupS_CB(i,j,1), vGroupS_CB(i,j,2)

                write(371, 5005) kx, ky, xiSurf2D_VB(i,j,1), xiSurf2D_VB(i,j,2)
                write(373, 5005) kx, ky, xiSurf2D_CB(i,j,1), xiSurf2D_CB(i,j,2)                                

                write(391, 5002) kx, ky, berrySS_VB(i,j)
                write(393, 5002) kx, ky, berrySS_CB(i,j)

                write(711, 5002) kx, ky, psiNorm_VB(i,j)
                write(713, 5002) kx, ky, psiNorm_CB(i,j)

                write(511, 5004) kx, ky, dble(dipoleSSVC(i,j,1)), dimag(dipoleSSVC(i,j,1)), &
                dble(dipoleSSVC(i,j,2)), dimag(dipoleSSVC(i,j,2))

                write(175, 5002) kx, ky, kmask( (i-1)*nky + j )
               
              enddo

              if(nky .gt. 1) then
                write(351,'(A)') ' '
                write(353,'(A)') ' '
                write(355,'(A)') ' '
                write(357,'(A)') ' '
                write(371,'(A)') ' '
                write(373,'(A)') ' '

                write(391,'(A)') ' '
                write(393,'(A)') ' '
      
                write(511,'(A)') ' '
                write(711,'(A)') ' ' 
                write(713,'(A)') ' ' 

                write(175, '(A)') ' ' 
              endif
            
            enddo

            close(175)

            close(351)
            close(353)
            close(355)
            close(357)
            close(371)
            close(373)
           
            close(391)
            close(393)
            
            close(511)
            close(711)
            close(713)


            if( allocated(Esurf2D_VB))    deallocate(Esurf2D_VB  )
            if( allocated(Esurf2D_CB))    deallocate(Esurf2D_CB  )
            if( allocated(vGroupS_VB))    deallocate(vGroupS_VB  )
            if( allocated(vGroupS_CB))    deallocate(vGroupS_CB  )
            if( allocated(xiSurf2D_VB))  deallocate(xiSurf2D_VB)
            if( allocated(xiSurf2D_CB))  deallocate(xiSurf2D_CB)
            if( allocated(berrySS_VB))   deallocate(berrySS_VB )
            if( allocated(berrySS_CB))   deallocate(berrySS_CB )
            if( allocated(dipoleSSVC))   deallocate(dipoleSSVC )
            if(allocated(psiNorm_VB))    deallocate(psiNorm_VB)
            if(allocated(psiNorm_CB))    deallocate(psiNorm_CB)

 5001     format (1e16.8)  
 5002     format (1e16.8,3(2x,1e16.8))  
 5003     format (1e16.8,2x,1e16.8)  
 5004     format (1f16.8,5(2x,1e16.8))  
 5005     format (1f16.8,3(2x,1e16.8))

          end subroutine print_bands

        end module solid
