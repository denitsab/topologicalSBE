!========================================================================!            
!======================= Gauss-Legendre-Lobatto =========================!
!========================================================================!        
        module legendre
          use math_constants
          implicit none
          contains

!---     Compute Legendre polynomial
          function Pleg(x, n)
             real(8) :: Pleg, x, y, yp, ypp
             integer :: n, i
             if(n.lt.1) then
               y = 1d0
             elseif(n.eq.1) then
               y = x
             else
               ypp = 1d0
               yp  = x
               do i=1,n-1
                 y = (dble(2*i+1)*x*yp - dble(i)*ypp) / dble(i+1)
                 ypp = yp
                 yp = y
               enddo
             endif
             Pleg = y
          end function Pleg

!---    Compute first derivative of Legendre polynomial
          function dPleg(x,n)
            real*8 dPleg, x
            integer n

            dPleg = n/(x**2-1d0)*(x*Pleg(x,n)-Pleg(x,n-1))
          end function dPleg    

          subroutine gleg_roots(N, x)
          ! Gaus-Legendre collocation points       
          ! using the eigenmethod described in J. Shen, T. Tang & L. Wang, "Spectral Methods:Algorithms, Analysis and Applications",
          ! Springer Series in Computational Mathematics, 41, Springer (2001)
            integer :: N 
            real*8 :: di 
            real*8, dimension(1:N) ::x 
            real*8, dimension(1:N) :: dv 
            real*8, dimension(1:N-1) :: subdv 
            integer :: i, info 
            dv = dzero 
            do i = 1,N-1
              di = dble(i)
              subdv(i) = di/dsqrt( (2.d0*di-1.d0)*(2.d0*di+1.d0) )
            enddo
            call dsterf( N, dv, subdv, info )
            if(info.ne.0) then
              write(*,*) ' Error in subrouitne DSTERF '
              write(*,*) ' Exited with code: ', info
              stop
            endif
            x(1:N) = dv(1:N)
            return
          end subroutine gleg_roots 

          subroutine gleg_weights(N, x, w)
            ! Formula: 
            !                     2 (1-x_i^2)                             2
            !        w_i = ---------------------------- =  ------------------------------
            !                (n+1)^2 [P_{n+1}(x_i)]^2       (1-x_i^2) [P_{n}'(x_i)]^2
            integer :: N, i
            real*8, dimension(1:N)  :: w, x
            if(N.le.1) then
              write(*,*) 'ERROR in SUBROUTINE gleg_lob_weights:'
              write(*,*) ' N should be greater than 1 !!! '
              stop
            endif

            do i = 1, N 
              w(i) = 2.d0/( (1.d0-x(i)**2) * dPleg(x(i), N)**2 )
            enddo 

            return
          end subroutine gleg_weights



          subroutine gleg_lob_roots(N, x)
          ! Gauss-Legendre-Lobatto collocation points
          ! using the eigenmethod described in J. Shen, T. Tang & L. Wang, "Spectral Methods:Algorithms, Analysis and Applications",
          ! Springer Series in Computational Mathematics, 41, Springer (2001)
            integer :: N
            real*8  :: di
            real*8, dimension(0:N)  :: x
            real*8, dimension(1:N-1) :: dv
            real*8, dimension(1:N-2) :: subdv
            integer :: i, info

            if(N.le.1) then
              write(*,*) 'ERROR in SUBROUTINE gleg_lob_roots:'
              write(*,*) ' N should be greater than 1 !!! '
              stop
            else if(N.eq.2) then
              x = (/ -1.d0, 1.d0 /)
              return
            else if(N.eq.3) then
              x = (/ -1.d0, 0.d0, 1.d0 /)
              return
            endif

            dv(1:N-1) = 0.d0

            do i = 1, N-2
              di = dble(i)
              subdv(i) = 0.5d0*dsqrt(di*(di+2.d0)/(di+0.5d0)/(di+1.5d0))
            enddo
            call dsterf( N-1, dv, subdv, info )
            if(info.ne.0) then
              write(*,*) ' Error in subrouitne DSTERF '
              write(*,*) ' Exited with code: ', info
              stop
            endif
            x(0) = -1.d0
            x(N) = 1.d0
            x(1:N-1) = dv(1:N-1)

            return

          end subroutine gleg_lob_roots

          subroutine gleg_lob_weights(N, x, w)
            ! Formula: 
            !                     2
            !        w_i = ------------------
            !              n (n+1) P_n(x_i)^2
            integer :: N, i
            real*8, dimension(0:N)  :: w, x
            if(N.le.1) then
              write(*,*) 'ERROR in SUBROUTINE gleg_lob_weights:'
              write(*,*) ' N should be greater than 1 !!! '
              stop
            endif

            w(0) = 2.d0/(dble(N)*dble(N+1))
            w(N) = w(0)

            do i=1,N-1
              if( (x(i).gt.(-1.d0)) .and. (x(i).lt.1.d0) ) then
                w(i) = 2.d0 / ( dble(N)*dble(N+1) * Pleg(x(i),N)**2 )
              else
                w(i) = -1.d0 / ( dble(N) * dble(N+1) )
              endif
            enddo
            
            return
          end subroutine gleg_lob_weights

          subroutine collocation(N, x, w)
            real*8, dimension(N) :: x, w
            integer :: N

            call gleg_lob_roots  ( N-1, x    )
            call gleg_lob_weights( N-1, x, w )

            return
          end subroutine collocation

           
        end module legendre

      