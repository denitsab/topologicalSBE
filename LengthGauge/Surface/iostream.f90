!========================================================================!            
!============================ IO Variables  =============================!
!========================================================================!

        module iostream 
        implicit none
        ! holds the input / output file paths and calcID, a unique label for each calculation
        save 
          integer Ilast, Olast
          character*200 Ipath, Opath, calcID
        end module iostream