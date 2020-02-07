!========================================================================!            
!============================ Math Constants ============================!
!========================================================================!

        module math_constants
        implicit none
        save
          real*8     :: pi, half, dzero
          complex*16 :: czero, iimag 
          parameter (pi = 4.d0*datan(1.0d0))     ! PI
          parameter (czero = dcmplx(0.d0, 0.d0)) ! complex zero
          parameter (dzero = 0.d0)               ! double zero 
          parameter (iimag = dcmplx(0.d0, 1.d0)) ! imaginary unit 
          parameter (half = 0.5d0)               ! 1/2 
        end module math_constants

!========================================================================!            
!============================== Math  ===================================!
!========================================================================!

        module math 
        implicit none

        contains 

        function derivative(y,n) result(dy) 
        ! differentiation according to a second order central-difference scheme
        ! y: array to be differentiated
        ! n: size of array
          implicit none
          integer :: n, i 
          real*8, dimension(n) :: y 
          real*8, dimension(n) :: dy 

          dy(1) = ( -3.d0*y(1) + 4.d0*y(2) - y(3) ) / 2.d0

          do i = 2, n-1 
            dy(i) = ( y(i+1) - y(i-1) )/2.d0
          enddo

          dy(n) = ( 3.d0*y(n) - 4.d0*y(n-1) + y(n-2) ) / 2.d0

          return
        end function derivative

        function derivativePBC(y,n) result(dy) 
        ! differentiation according to a second order central-difference scheme under periodic boundary conditions
        ! y: array to be differentiated
        ! n: size of array
          implicit none
          integer :: n
          real*8, dimension(n) :: y 
          real*8, dimension(n) :: dy 

          dy(1) = (y(2) - y(n))/1.d0
          dy(2:n-1) = (y(3:n) - y(1:n-2))/2.d0
          dy(n) = (y(1)-y(n-1))/1.d0

          return
        end function derivativePBC

        function derivativeC(y,n) result(dy)
        ! Routine "derivative", re-written to handle complex input arrays 
          implicit none
          integer :: n, i 
          complex*16, dimension(n) :: y 
          complex*16, dimension(n) :: dy 

          dy(1) = ( -3.d0*y(1) + 4.d0*y(2) - y(3) ) / 2.d0

          do i = 2, n-1 
            dy(i) = ( y(i+1) - y(i-1) )/2.d0
          enddo

          dy(n) = ( 3.d0*y(n) - 4.d0*y(n-1) + y(n-2) ) / 2.d0

          return
        end function derivativeC

        function derivativeC_PBC(y,n) result(dy) 
        ! Routine "derivativePBC", re-written to handle complex input arrays 
          implicit none
          integer :: n
          complex*16, dimension(n) :: y 
          complex*16, dimension(n) :: dy 

          dy(1) = (y(2) - y(n))/1.d0
          dy(2:n-1) = (y(3:n) - y(1:n-2))/2.d0
          dy(n) = (y(1)-y(n-1))/1.d0

          return
        end function derivativeC_PBC

        function DMatrixOP(n, dx) result(drvmat)
        ! Return the differential operator in a matrix form given a grid size of n points and a grid spacing of dx
        ! Array differentiation is performed by matrix multiplication with drvmat
          implicit none
          integer :: n,i 
          real*8  :: dx 
          real*8, dimension(n,n) :: drvmat

          drvmat = 0.d0
          do i = 2, n-1 
            drvmat(i, i+1) =  1.d0/(2.d0*dx)
            drvmat(i, i-1) = -1.d0/(2.d0*dx)
          enddo

          drvmat(1,1) = -3.d0/(2.d0*dx)
          drvmat(1,2) =  2.d0/dx
          drvmat(1,3) = -1.d0/(2.d0*dx)

          drvmat(n,n)   =  3.d0/(2.d0*dx)
          drvmat(n,n-1) = -2.d0/dx
          drvmat(n,n-2) =  1.d0/(2.d0*dx)

          return
        end function DMatrixOP


        function DMatrixOP_PBC(n, dx) result (drvmat)
        ! analogous to DMatrixOP, but adapted to handle complex arrays
          implicit none
          integer :: n, i 
          real*8  :: dx 
          real*8, dimension(n,n) :: drvmat

          drvmat = 0.d0
          do i = 2, n-1 
            drvmat(i, i+1) =  1.d0/(2.d0*dx)
            drvmat(i, i-1) = -1.d0/(2.d0*dx)
          enddo

          drvmat(1, n) = -1.d0/(1.d0*dx)
          drvmat(1, 2) =  1.d0/(1.d0*dx)
          drvmat(n, 1) =  1.d0/(1.d0*dx)
          drvmat(n, n-1) = -1.d0/(1.d0*dx)

          return
        end function DMatrixOP_PBC

        complex*16 function simpsC(y, n) 
        ! Implementation of the Simpson integration rule
        ! y: input vector to be differentiated
        ! n: size of vector y
          implicit none
          integer :: n
          complex*16, dimension(n) :: y 
          complex*16 :: c3

          c3 = dcmplx(1.d0/3.d0)

          simpsC = c3 * ( y(1) + y(n) + dcmplx(4.d0,0.d0)*sum(y(2:n-1:2)) + dcmplx(2.d0,0.d0)*sum(y(3:n-2:2) ) )
          return
        end function simpsC

        function BZint2DSimpson(f,n, nx, ny, dx, dy) result(S2D)
          implicit none
          complex*16, dimension(n)     :: f 
          complex*16                   :: S2D
          integer                      :: n, nx, ny 
          real*8                       :: dx, dy 
          integer                      ::  i, j
          complex*16, dimension(nx,ny) :: fmat 

          do j = 1, ny
            do i = 1, nx  
              fmat(i,j) = f( (i-1)*ny + j )
            enddo
          enddo

          S2D = fmat(1,1) + fmat(1,ny) + fmat(nx,1) + fmat(nx,ny)
          S2D = S2D + 4.d0*sum( fmat(1,3:ny-1:2) ) + 2.d0*sum( fmat(1,2:ny-2:2) ) + &
          + 4.d0*sum( fmat(nx,3:ny-1:2) ) + 2.d0*sum( fmat(nx,2:ny-2:2) ) 
          S2D = S2D + 4.d0*sum( fmat(3:nx-1:2,1) ) + 2.d0*sum( fmat(2:nx-2:2,1) ) + &
          + 4.d0*sum( fmat(3:nx-1:2,ny) ) + 2.d0*sum( fmat(2:nx-2:2,ny) ) 
          S2D = S2D + 16.d0*sum( sum( fmat(3:nx-1:2,3:ny-1:2),1 ) ) + &
          8.d0*sum( sum( fmat(3:nx-1:2,2:ny-2:2),1 ) ) + 8.d0*sum( sum( fmat(2:nx-2:2,3:ny-1:2),1 ) ) + &
          4.d0*sum( sum( fmat(2:nx-2:2,2:ny-2:2),1 ) )

          S2D = 1.d0/9.d0 * dx * dy * S2D
          return
        end function  BZint2DSimpson

        function BZint2DTZ(f,n, nx, ny, dx, dy) result(S2D)
          implicit none
          complex*16, dimension(n)     :: f 
          complex*16                   :: S2D
          integer                      :: n, nx, ny 
          real*8                       :: dx, dy 
          integer                      ::  i, j
          complex*16, dimension(nx,ny) :: fmat 

          do j = 1, ny
            do i = 1, nx  
              fmat(i,j) = f( (i-1)*ny + j )
            enddo
          enddo

          S2D = fmat(1,1) + fmat(1,ny) + fmat(nx,1) + fmat(nx,ny)
          S2D = S2D + 2.d0*sum( fmat(2:nx-1,1) ) + 2.d0*sum( fmat(2:nx-1,ny) ) + &
          2.d0*sum( fmat(1, 2:ny-1) ) + 2.d0*sum( fmat(nx, 2:ny-1) ) + &
          4.d0*sum( sum( fmat(2:nx-1,2:ny-1), 1 ) )
          

          S2D = 1.d0/4.d0 * dx * dy * S2D
          return
        end function  BZint2DTZ

        end module math
