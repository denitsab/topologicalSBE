!========================================================================!            
!============================ Phy Constants =============================!
!========================================================================!

        module phys_constants
        ! contains commonly encountered natural constants and conversion factors SI <--> a.u.
        implicit none
        save 
          real*8 :: atu, Eh, Iau, cvac0, eV, h_planck, h_bar, m_e, el, eps0, &
              aBohr,  el_field, hartree
          parameter (atu = 2.418884326509d-17)     ! [s] atomic time unit
          parameter (cvac0  = 299792458.d0)           ! [m/s] speed of light in SI
          parameter (Iau = 3.51d16)                ! [W/cm2] atomic unit of intensity
          parameter (Eh  = 4.359744650d-18)        ! [J] Hartree energy
          parameter (eV  = 1.602176565d-19)        ! [J] ElectronVolt
          parameter (h_planck = 6.62606957d-34  )  ! [Js     ]
          parameter (h_bar    = 1.054571800d-34 )  ! [Js     ]
          parameter (m_e      = 9.10938291d-31  )  ! [kg     ] 
          parameter (eps0     = 8.854187817d-12 )  ! [Fm^(-1)]
          parameter (el       = 1.602176565d-19 )  ! [C      ]
          parameter (aBohr    = 5.291772192d-11   )! [m      ] 
          parameter (hartree  = 4.35974417d-18    )! [J      ] 
          parameter (el_field = 5.14220652d11     )! [Vm^(-1)] 

        end module phys_constants