xi_v(1) = (ky*(2*kx**3*Rt - 6*kx*ky**2*Rt + 
     -      Sqrt(Atld**2*(kx**2 + ky**2) + 
     -        4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2)))/
     -  (2.*(kx**2 + ky**2)*
     -    Sqrt(Atld**2*(kx**2 + ky**2) + 
     -      4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2))

            xi_v(2) =     -(kx*(2*kx**3*Rt - 6*kx*ky**2*Rt + 
     -       Sqrt(Atld**2*(kx**2 + ky**2) + 
     -         4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2)))/
     -  (2.*(kx**2 + ky**2)*
     -    Sqrt(Atld**2*(kx**2 + ky**2) + 
     -      4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2)) 

            xi_c(1) =         (ky*(-2*kx**3*Rt + 6*kx*ky**2*Rt + 
     -      Sqrt(Atld**2*(kx**2 + ky**2) + 
     -        4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2)))/
     -  (2.*(kx**2 + ky**2)*
     -    Sqrt(Atld**2*(kx**2 + ky**2) + 
     -      4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2))

            xi_c(2) = -(kx*(-2*kx**3*Rt + 6*kx*ky**2*Rt + 
     -       Sqrt(Atld**2*(kx**2 + ky**2) + 
     -         4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2)))/
     -  (2.*(kx**2 + ky**2)*
     -    Sqrt(Atld**2*(kx**2 + ky**2) + 
     -      4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2))

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
              (2.d0*(Atld**2*(kx**2 + ky**2) + &
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
              (2.d0*(Atld**2*(kx**2 + ky**2) + &
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

            omg_v =    (-2.d0*Atld**2*kx*(kx**2 - 3.d0*ky**2)*Rt)/&
              (Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)**1.5

            omg_c =       (2.d0*Atld**2*kx*(kx**2 - 3.d0*ky**2)*Rt)/&
              (Atld**2*(kx**2 + ky**2) + &
              4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)**1.5


            xigap(1) = (-2*ky*(kx**3 - 3*kx*ky**2)*Rt)/((kx**2 + ky**2)*Sqrt(Atld**2*(kx**2 + ky**2) + 4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2))
            xigap(2) = (2*kx**2*(kx**2 - 3*ky**2)*Rt)/((kx**2 + ky**2)*Sqrt(Atld**2*(kx**2 + ky**2) + 4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2))

            e_v = C0t + C2t*(kx**2 + ky**2) - Sqrt(Atld**2*(kx**2 + ky**2) + 4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2)
            e_c = C0t + C2t*(kx**2 + ky**2) + Sqrt(Atld**2*(kx**2 + ky**2) + 4*kx**2*(kx**2 - 3*ky**2)**2*Rt**2)

            egap = 2.d0*Sqrt(Atld**2*(kx**2 + ky**2) + 4.d0*kx**2*(kx**2 - 3.d0*ky**2)**2*Rt**2)

