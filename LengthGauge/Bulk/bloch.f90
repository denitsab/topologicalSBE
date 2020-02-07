!========================================================================!            
!========================= SBE equations ================================!
!========================================================================!

        module bloch
          use math_constants
          use math  
         ! use phys_constants
          use solid
          use field 
          use grid 
          implicit none 
          ! The current module contains the subroutines that calculate the derivatives of the density matrix elements 
          ! These routines are called by the corresponding integrator (e.g. ZVODE)
          ! Consequently, they are adapted to the requirements of the corresponding integrator
          ! see the source file of the integrator for more information on the input requirements 
          real*8    :: fabs_vmask, fbrd_vmask, q_vmask
          
          save fabs_vmask, fbrd_vmask, q_vmask


          contains 

          ! Derivatives of the density matrix elements in the moving k-frame 
          subroutine derivs(neq, x, y, dydx, rpar, ipar)
            implicit none
            ! neq : total number of (coupled) ODEs. Here, neq = 2*nkpts, i.e. populations and coherences as a function of kx, ky, for all points of hte grid 
            ! x   : here, x is the time variable 
            ! y   : the solution vector, i.e. the vector storing nc and pi. Size: 2*nkpts 
            ! dydx : derivative, to be computed by the routine. same dimensions as y 
            ! rpar, ipar - dummy parameters to be compatible with the requirements of ZVODE 
            real*8                          ::    x         ! time variable 
            complex*16,dimension(neq)       ::    y,dydx    ! solution vector, derivative vector 
            real*4, dimension(1)            ::    rpar      ! dummy variable 
            integer, dimension(1)           ::    ipar      ! dummy variable 
            integer                         ::    neq       ! number of ODEs, 2*nkpts 
            real*8                          ::    Et, At    ! electric field / vector potential in 1D    
            integer                         ::    i         ! indices 
            complex*16, dimension (nkpts)   ::    wK                     ! this is w(K) = nv(K) - nc(K)
            real*8, dimension(nkpts, 2 )    ::    kshift                 ! in the moving frame, kshift = k + A(t) 
            real*8, dimension(2)            ::    Et2, At2, kxy          ! Et2: 2D el. field vector (for future use), kxy: (kx, ky)
            real*8, dimension(3)            ::    kxyz 
            complex*16, dimension(2)        ::    dmat , Xigap           ! dipole matrix element, difference of Berry connections 
            complex*16, dimension(3)        ::    dmat3
            real*8                          ::    Egap                   ! energy band gap at a given k                         
            real*8                          ::    fCB, fVB

            ! El. field and vector potential at t 
            !Et  = Efun(x)
            !At  = Afun(x) 
            Et2 = Efun2(x)
            At2 = Afun2(x)


            ! Shift the k-grid by At in the x-direction 
            kshift      = kxy_grid
            !kshift(:,1) = kshift(:,1) + At 
            kshift(:,1) = kshift(:,1) + At2(1)
            kshift(:,2) = kshift(:,2) + At2(2)

            ! Storage scheme: 
            ! elements 1 to nkpts of y corresond to nc 
            ! elements nkpts+1 to 2*nkpts correspond to pi 
            ! the same is used for the derivative dydx 

            dydx = czero 

            if( M0 .lt. 0.d0 ) then 
              fCB = f0BC
              fVB = 0.d0 
            else
              fCB = 0.d0 
              fVB = f0BC 
            endif    


!$omp parallel do private(i, kxyz, Egap, Xigap, dmat) &
!$omp& shared(nkpts, kshift, Et, y, dydx, Et2, fVB, fCB, sigp)
!!$omp& shared(nkpts, kshift, Et, ppp, npv, npc, dotnpc, dotnpv, dotppp, wK)
            do i = 1, nkpts ! This loop samples all momentum points, see the definition of kxy_grid

              kxyz = (/ kshift(i,1), kshift(i,2), kzpoint /)

              call get_crystal_kxyz(kxyz, egap, xigap, dmat, fVB, fCB, sigp)


              ! The lines below are the ODEs for the population/coherence derivatives 
              ! Note the Dvc routine retireves d_{vc}, i.e. the complex-conjugate of d_{cv} in the SM, threrefore, there is one cc operation less compared to Eq. (11) in the paper SM

              !dydx(        i) = - dcmplx(2.d0, 0.d0) * Et * aimag( (dmat(1)) * y(nkpts + i) )
              !dydx(nkpts + i) = -iimag * ( Egap + Et*dble(Xigap(1)) - iimag/Tdephas_2 ) * y(nkpts + i) - &
              !imag * Et * conjg(dmat(1)) * ( dcmplx(1.d0, 0.d0) - dcmplx(2.d0,0.d0) * y(i) )

              dydx(        i) = - dcmplx(2.d0, 0.d0) *  aimag( (Et2(1)*(dmat(1)) + Et2(2)*dmat(2) ) & 
                * y(nkpts + i) )
              dydx(nkpts + i) = -iimag * ( Egap + ( Et2(1)*dble( Xigap(1)) + Et2(2)*dble( Xigap(2)) )  &
                - iimag/Tdephas_2 ) * y(nkpts + i) - &
              iimag * (Et2(1)*conjg(dmat(1)) + Et2(2)*conjg(dmat(2)) ) * ( dcmplx(1.d0, 0.d0) - dcmplx(2.d0,0.d0) * y(i) )
            enddo
!$omp end parallel do                      


        
          return
          end subroutine derivs

          subroutine currents(Y, tvec)
            implicit none
          ! A routine calculating the polarization, the inter- and the intraband currents from the propagated density matrix as a function of time 
          ! Moving frame 
            complex*16, dimension(nmax, ntpts)  :: Y                       ! density matrix as a function of k and t. Format: (2*nkpts, ntpts)
            real*8, dimension    (ntpts)        :: tvec                    ! time vector 
            complex*16, dimension(nkpts)        :: npv, npc, ppp           ! auxiliary arrays for band populations and coherence
            real*8, dimension(3)                :: kxyz, vgr               ! 3D vectors for storing (kx, ky) and group velocity 
            real*8, dimension(2)                :: vgr_v, vgr_c 
            real*8, dimension(2)                :: kxy, At2, Et2           ! 2D vectors for storing (kx, ky) and group velocity 
            complex*16, dimension(3)            :: bOmega                  ! matrix element between CB and VB, 3D 
            complex*16                          :: omg_v, omg_c
            complex*16, dimension(2)            :: dmat
            complex*16, dimension(3)            :: dmat3
            complex*16, dimension(nkpts,2)      :: dmatkgrid               ! array for storing dcv as a function of k, moving frame 
            real*8, dimension(nkpts,2)          :: vgrkgrid_c, vgrkgrid_v  ! array for storing the group velocity as a function of k, moving frame
            real*8, dimension(nkpts,2)          :: vankgrid_c, vankgrid_v  ! array for storing the anpmalous velocity as a function of k, moving frame
            integer                             :: i, j, k            
            real*8                              :: t,  At , Et, bOmegaZ    ! time point, vector potential, el. field, berry curvature 
            real*8, dimension(nkpts,2)          :: kshift                  ! shifted momentum grid, kshift = k + A(t) 
            real*8                              :: fBC_v, fBC_c
            real*8                              :: vmask, kabs, kbuff, kborder , kpar

            dmatkgrid (1:nkpts, 1:2 ) = czero
            vgrkgrid_v(1:nkpts, 1:2 ) = dzero
            vgrkgrid_c(1:nkpts, 1:2 ) = dzero

            if(M0 .lt. 0) then
              fBC_v = 0.d0
              fBC_c = f0BC
            else
              fBC_c = 0.d0
              fBC_v = f0BC
            endif

            do i = 1, ntpts ! Loop over all possible time points 
              kshift = kxy_grid
              t = tvec(i)
              !At = Afun(t)
              !Et = Efun(t)
              At2 = Afun2(t)
              Et2 = Efun2(t)
              ! Shift the k-grid in x-direction 
              !kshift(:,1) = kshift(:,1) + At 
              kshift(:,1) = kshift(:,1) + At2(1) 
              kshift(:,2) = kshift(:,2) + At2(2) 

              kabs = fabs_vmask*ky1d(nky)
              kborder = fbrd_vmask*ky1d(nky)
              kbuff = kborder - kabs 

!!$omp parallel do private(j, kxy, kxyz,  dmat, vgr, bOmega, bOmegaZ ) &
!$omp parallel do private(j, kxyz,  dmat,dmat3, vgr_v, vgr_c, omg_v,omg_c, kpar, vmask ) &
!$omp& shared(nkpts, kshift, dmatkgrid, vgrkgrid_v, vgrkgrid_c,vankgrid_v,vankgrid_c,Et,fBC_v,fBC_c, Et2, At2,f0BC,sigp)
              do j = 1, nkpts  ! Loop over momentum grid; calculate the dipoles, the berry curvature, the group and anpmalous velocities at each shifted k 
                !kxy = kshift(j,:)
                kxyz = (/ kshift(j,1), kshift(j,2), kzpoint /)
                call get_crystal_velocities_kxyz(kxyz, vgr_v, vgr_c, omg_v, omg_c, dmat, 0.d0, 0.d0,sigp)
                dmat3 = DvcB( kxyz, icomp, icomp, f0BC )

                kpar = sqrt( kxyz(1) ** 2 + kxyz(2)**2 )

                if (fabs_vmask .gt. 1.d-6 ) then              
                  if( kpar .lt. kabs ) then 
                    vmask = 1.d0 !* (eps1+eps2+eps3)/3.d0 
                  else if( kpar .lt. kborder ) then
                    vmask = abs(cos(abs(kabs-kpar) * pi / (2.d0*kbuff)))**q_vmask!0.125!*(eps1+eps2+eps3)/3.d0 
                    !vmask = 1.d0 + erf( -(kpar-kabs)**2/(kbuff**2) )
                  else
                    vmask = 0.d0
                  endif
                else
                  vmask = 1.d0 
                endif
                
                dmatkgrid(j,:) = vmask*dmat3 (1:2)
                vgrkgrid_v(j,:) = vmask*vgr_v(1:2)
                !vgrkgrid_v(j,:) = ( sum(abs(vgr_v(1:2))**2) )
                !vankgrid_v(j,:) = (/ 0.d0, - Et * dble(omg_v) /)
                vankgrid_v(j,:) = (/ Et2(2)*dble(omg_v), - Et2(1) * dble(omg_v) /)

                vgrkgrid_c(j,:) = vmask*vgr_c(1:2)
                !vgrkgrid_c(j,:) = ( sum(abs(vgr_c(1:2))**2) )
                vankgrid_c(j,:) = (/ Et2(2)*dble(omg_c), - Et2(1) * dble(omg_c) /)

              enddo ! k
!$omp end parallel do

       

              ! Copy elements of the Y-vector to the populations/coherence arrays according to storage scheme
              npc = Y(1:nkpts, i)
              ppp = Y(nkpts+1:2*nkpts, i)
              npv = dcmplx(1.d0, 0.d0) - npc

              do k = 1,2 ! Loop over x, y directions 
               ! Calculate the polarization / intraband current by integrating over the shifted Brillouin zone at this t 
               !Pt(i,k) = - sum( (dmatkgrid(:,k)) * ppp ) * dkx * dky
               Pt(i,k) = - BZint2DTZ( (dmatkgrid(:,k)) * ppp, nkpts, nkx, nky, dkx, dky ) 
               Pt(i,k) = Pt(i,k) + conjg(Pt(i,k))

               !Jt(i,k) = (- sum( (vgrkgrid_v(:,k)+vankgrid_v(:,k)) * npv )  - &
               !             sum( (vgrkgrid_c(:,k)+vankgrid_c(:,k)) * npc )) * &
               !             dkx * dky
               Jt(i,k) = (- BZint2DTZ( (vgrkgrid_v(:,k)+vankgrid_v(:,k)) * npv, nkpts, nkx, nky, dkx, dky )  - &
                            BZint2DTZ( (vgrkgrid_c(:,k)+vankgrid_c(:,k)) * npc , nkpts, nkx, nky, dkx, dky ) )

             enddo ! x,y 
           enddo   ! t

           ! Calculate the interband current as a time derivative of the polarization. THe same array (and symbol) is used to spare memory.
           Pt(:,1) = derivativeC(Pt(:,1), ntpts)/dt 
           Pt(:,2) = derivativeC(Pt(:,2), ntpts)/dt 

           return

         end subroutine currents

         subroutine print_group_velocity(iunit, t)
            implicit none
            integer :: iunit
            real*8   :: t 
            integer :: i 
            real*8, dimension(3)                :: kxyz, vgr               ! 3D vectors for storing (kx, ky) and group velocity 
            real*8, dimension(2)                :: vgr_v, vgr_c 
            real*8, dimension(2)                :: At2          ! 2D vectors for storing (kx, ky) and group velocity 
            complex*16                          :: omg_v, omg_c
            complex*16, dimension(2)            :: dmat
            real*8, dimension(nkpts,2)          :: kshift                  ! shifted momentum grid, kshift = k + A(t) 

            kshift = kxy_grid
              
            At2 = Afun2(t)
             
            ! Shift the k-grid in x-direction 
            !kshift(:,1) = kshift(:,1) + At 
            kshift(:,1) = kshift(:,1) + At2(1) 
            kshift(:,2) = kshift(:,2) + At2(2) 

            do i = 1, nkpts 

              kxyz = (/ kshift(i,1), kshift(i,2), kzpoint /)
              call get_crystal_velocities_kxyz(kxyz, vgr_v, vgr_c, omg_v, omg_c, dmat, 0.d0, 0.d0,sigp)
              write(7500+iunit, 2055) t, kxy_grid(i,1), kxy_grid(i,2), &
              vgr_v(1), vgr_v(2)
              write(7700+iunit, 2055) t, kxy_grid(i,1), kxy_grid(i,2), &
              vgr_c(1), vgr_c(2)
              
              if(mod(i, nky).eq.0) then 
                write(7500+iunit,'(A)') ' ' 
                write(7700+iunit,'(A)') ' ' 
              endif
            enddo
            return
 2055     format (1e16.8,4(2x,1e16.8))  
          end subroutine print_group_velocity

        subroutine jex()
          implicit none
          return
        end subroutine jex
            
        end module bloch
