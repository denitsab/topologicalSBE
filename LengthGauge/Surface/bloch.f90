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
            real*8, dimension(2)            ::    Et2, At2, kxy               ! Et2: 2D el. field vector (for future use), kxy: (kx, ky)
            complex*16, dimension(2)        ::    dmat , Xigap           ! dipole matrix element, difference of Berry connections 
            real*8                          ::    Egap                   ! energy band gap at a given k                         
            real*8                          ::    Etx, Ety, Atx, Aty
            complex*16, dimension(2)        ::    dmatrot, xirot 

            ! El. field and vector potential at t 
            !Et  = Efun(x)
            !At  = Afun(x) 
            Et2 = Efun2(x)
            At2 = Afun2(x)

          
            ! Shift the k-grid by At in the x-direction 
            kshift      = kxy_grid
            !kshift(:,1) = kshift(:,1) + At 
            !kshift(:,1) = kshift(:,1) + At2(1)
            !kshift(:,2) = kshift(:,2) + At2(2)

            kshift(:,1) = kshift(:,1) + At2(1)
            kshift(:,2) = kshift(:,2) + At2(2)

            ! Storage scheme: 
            ! elements 1 to nkpts of y corresond to nc 
            ! elements nkpts+1 to 2*nkpts correspond to pi 
            ! the same is used for the derivative dydx 

            dydx = czero 

           


!$omp parallel do private(i, kxy, Egap, Xigap, dmat ) &
!$omp& shared(nkpts, kshift, Et, y, dydx, f0glob, sigp, Et2, theta)
            do i = 1, nkpts ! This loop samples all momentum points, see the definition of kxy_grid

              kxy = kshift(i,:)

              call get_crystal_kxy_compact(kxy, egap, xigap, dmat, f0glob, sigp)

              !dmatrot(1) =   dmat(1) * cos(theta) + dmat(2) * sin(theta)
              !dmatrot(2) = - dmat(1) * sin(theta) + dmat(2) * cos(theta)
              !xirot(1)   =   xigap(1) * cos(theta) + xigap(2) * sin(theta)
              !xirot(2)   = - xigap(1) * sin(theta) + xigap(2) * cos(theta)


              ! The lines below are the ODEs for the population/coherence derivatives 
              ! Note the Dvc routine retireves d_{vc}, i.e. the complex-conjugate of d_{cv} in the SM, threrefore, there is one cc operation less compared to Eq. (11) in the paper SM

              !dydx(        i) = - dcmplx(2.d0, 0.d0) * Et * aimag( (dmat(1)) * y(nkpts + i) )
              !dydx(nkpts + i) = -iimag * ( Egap + Et*dble(Xigap(1)) - iimag/Tdephas_2 ) * y(nkpts + i) - &
              !iimag * Et * conjg(dmat(1)) * ( dcmplx(1.d0, 0.d0) - dcmplx(2.d0,0.d0) * y(i) )
              dydx(        i) = - dcmplx(2.d0, 0.d0) * aimag( &
                ( Et2(1)*(dmat(1)) + Et2(2)*(dmat(2)) ) * y(nkpts + i) )
              dydx(nkpts + i) = -iimag * ( Egap + (Et2(1)*dble(xigap(1)) + Et2(2)*dble(xigap(2)) )&
               - iimag/Tdephas_2 ) * y(nkpts + i) - &
              iimag * ( Et2(1)*conjg(dmat(1)) + Et2(2)*conjg(dmat(2)) ) &
               * ( dcmplx(1.d0, 0.d0) - dcmplx(2.d0,0.d0) * y(i) )
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
            real*8, dimension(2)                :: kxy, vgr               ! 3D vectors for storing (kx, ky) and group velocity 
            real*8, dimension(2)                :: vgr_v, vgr_c 
            complex*16                          :: omg_v, omg_c
            complex*16, dimension(2)            :: dmat
            complex*16, dimension(nkpts,2)      :: dmatkgrid               ! array for storing dcv as a function of k, moving frame 
            real*8, dimension(nkpts,2)          :: vgrkgrid_c, vgrkgrid_v  ! array for storing the group velocity as a function of k, moving frame
            real*8, dimension(nkpts,2)          :: vankgrid_c, vankgrid_v  ! array for storing the anpmalous velocity as a function of k, moving frame
            integer                             :: i, j, k            
            real*8                              :: t,  At , Et, bOmegaZ    ! time point, vector potential, el. field, berry curvature 
            real*8, dimension(nkpts,2)          :: kshift                  ! shifted momentum grid, kshift = k + A(t) 
            real*8, dimension(2)                :: At2, Et2 
            complex*16, dimension(2)            :: dmatrot 
            real*8, dimension(2)                :: vgr_rot
            real*8                              :: Etx, Ety, Atx, Aty
            complex*16, dimension(ntpts,2)      :: Jtcry, Ptcry
            real*8                              :: vmask, kabs, kbuff, kborder , kpar

            dmatkgrid (1:nkpts, 1:2 ) = czero
            vgrkgrid_v(1:nkpts, 1:2 ) = dzero
            vgrkgrid_c(1:nkpts, 1:2 ) = dzero

            Jtcry = czero 
            Ptcry = czero

            do i = 1, ntpts ! Loop over all possible time points 
              kshift = kxy_grid
              t = tvec(i)
              !At = Afun(t)
              !Et = Efun(t)
              At2 = Afun2(t)
              Et2 = Efun2(t) 

              

              ! Shift the k-grid in x-direction 
              !kshift(:,1) = kshift(:,1) + At 
              !kshift(:,1) = kshift(:,1) + At2(1)
              !kshift(:,2) = kshift(:,2) + At2(2)

              kshift(:,1) = kshift(:,1) + At2(1)
              kshift(:,2) = kshift(:,2) + At2(2)

              kabs = fabs_vmask*ky1d(nky)
              kborder = fbrd_vmask*ky1d(nky)
              kbuff = kborder - kabs 

              

!$omp parallel do private(j, kxy,  dmat, vgr_v, vgr_c, omg_v,omg_c, vgr_rot, kpar, vmask ) &
!$omp& shared(nkpts, kshift, dmatkgrid, vgrkgrid_v, vgrkgrid_c,vankgrid_v,vankgrid_c,Et, Et2, theta )
              do j = 1, nkpts  ! Loop over momentum grid; calculate the dipoles, the berry curvature, the group and anpmalous velocities at each shifted k 
                kxy = kshift(j,:)
                
                call get_crystal_velocities_kxy_compact(kxy, vgr_v, vgr_c, omg_v, omg_c, dmat, f0glob, sigp)

                !dmatrot(1) =   dmat(1) * cos(theta) + dmat(2) * sin(theta)
                !dmatrot(2) = - dmat(1) * sin(theta) + dmat(2) * cos(theta)

                !vgr_rot(1) =   vgr_v(1) * cos(theta) + vgr_v(2) * sin(theta)
                !vgr_rot(2) = - vgr_v(1) * sin(theta) + vgr_v(2) * cos(theta)
                !vgr_rot(1) =   vgr_v(1) * cos(theta) - vgr_v(2) * sin(theta)
                !vgr_rot(2) =   vgr_v(1) * sin(theta) + vgr_v(2) * cos(theta)

                kpar = sqrt( kxy(1) ** 2 + kxy(2)**2 )

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


                dmatkgrid(j,:) = vmask * dmat (1:2)
                !dmatkgrid(j,:) = dmatrot (1:2)
                vgrkgrid_v(j,:) = vmask * vgr_v(1:2)
                !vgrkgrid_v(j,:) = vgr_rot(1:2)
                !vankgrid_v(j,:) = (/ 0.d0, - Et * dble(omg_v) /)
                !vankgrid_v(j,:) = (/ Et2(2) * dble(omg_v), - Et2(1) * dble(omg_v) /)
                vankgrid_v(j,:) = (/ Et2(2) * dble(omg_v), - Et2(1) * dble(omg_v) /)

                !vgr_rot(1) =   vgr_c(1) * cos(theta) - vgr_c(2) * sin(theta)
                !vgr_rot(2) =   vgr_c(1) * sin(theta) + vgr_c(2) * cos(theta)

                !vgr_rot(1) =   vgr_c(1) * cos(theta) + vgr_c(2) * sin(theta)
                !vgr_rot(2) = - vgr_c(1) * sin(theta) + vgr_c(2) * cos(theta)

                vgrkgrid_c(j,:) = vmask * vgr_c(1:2)
                !vgrkgrid_c(j,:) = vgr_rot(1:2)
                !vankgrid_c(j,:) = (/ 0.d0, - Et * dble(omg_c) /)
                !vankgrid_c(j,:) = (/ Et2(2) * dble(omg_c), - Et2(1) * dble(omg_c) /)
                vankgrid_c(j,:) = (/ Et2(2) * dble(omg_c), - Et2(1) * dble(omg_c) /)

              enddo ! k
!$omp end parallel do


              ! Copy elements of the Y-vector to the populations/coherence arrays according to storage scheme
              npc = Y(1:nkpts, i)
              ppp = Y(nkpts+1:2*nkpts, i)
              npv = dcmplx(1.d0, 0.d0) - npc

              do k = 1,2 ! Loop over x, y directions 
               ! Calculate the polarization / intraband current by integrating over the shifted Brillouin zone at this t 
               !Pt(i,k) = - simpsC( (dmatkgrid(:,k)) * ppp , nkpts) * dkx * dky
               Ptcry(i,k) = - BZint2DTZ( (dmatkgrid(:,k)) * ppp, nkpts, nkx, nky, dkx, dky )
               Ptcry(i,k) = Ptcry(i,k) + conjg(Ptcry(i,k))

               !Jt(i,k) = (- simpsC( (vgrkgrid_v(:,k)+vankgrid_v(:,k)) * npv, nkpts )  - &
               !             simpsC( (vgrkgrid_c(:,k)+vankgrid_c(:,k)) * npc, nkpts )) * &
               !             dkx * dky
               Jtcry(i,k) = (- BZint2DTZ( (vgrkgrid_v(:,k)+vankgrid_v(:,k)) * npv, nkpts, nkx, nky, dkx, dky )  - &
                            BZint2DTZ( (vgrkgrid_c(:,k)+vankgrid_c(:,k)) * npc, nkpts, nkx, nky, dkx, dky ))                            

             enddo ! x,y 
           enddo   ! t

           ! Calculate the interband current as a time derivative of the polarization. THe same array (and symbol) is used to spare memory.
           Ptcry(:,1) = derivativeC(Ptcry(:,1), ntpts)/dt 
           Ptcry(:,2) = derivativeC(Ptcry(:,2), ntpts)/dt 

           Pt(:,1) =   Ptcry(:,1) * cos(theta) + Ptcry(:,2) * sin(theta)
           Pt(:,2) = - Ptcry(:,1) * sin(theta) + Ptcry(:,2) * cos(theta)

           Jt(:,1) =   Jtcry(:,1) * cos(theta) + Jtcry(:,2) * sin(theta)
           Jt(:,2) = - Jtcry(:,1) * sin(theta) + Jtcry(:,2) * cos(theta)

           return

         end subroutine currents

        subroutine jex()
          implicit none
          return
        end subroutine jex
            
        end module bloch