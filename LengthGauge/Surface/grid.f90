!========================================================================!            
!========================= Grid  Variables  =============================!
!========================================================================!
        module grid 
          use math_constants
        ! number of momentun and time grid points 
        ! ntpts : number of points on the time grid 
        ! nkx   : number of points in the x momentum directions
        ! nky   : number of points in the y momentum direction 
        ! nkpts : nkx*nky -> total number of momentum points 
        ! nmax  : reserved (for future use)
        implicit none
          integer :: ntpts, nkpts, nkx, nky
          integer :: nmax
          save ntpts, nkpts, nmax, nky, nkx
        end module grid