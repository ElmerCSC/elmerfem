!
! A partir des 6 viscosites eta_i (loi type boehler)
! retourne les 9 viscosites de la matrice en expression de Voigt      
! dans le repere d'orthotropie
!      
      !  s11 = A11 d11 + A12 d22 + A13 d33
      !  s22 = A12 d11 + A22 d22 + A23 d33
      !  s33 = A13 d11 + A23 d22 + A33 d33
      !  s23 = A44 2*d23
      !  s31 = A55 2*d31
      !  s12 = A66 2*d12

!  Retourne etaVO(9) = Xi
!
!              |X1 X2 X3         |
!              |X2 X4 X5         |
!     Aij =    |X3 X5 X6         |
!              |         X7      |
!              |            X8   |
!              |               X9|

      Subroutine viscVoigt(v,etaVO)

      Implicit None
      Real, Dimension(6) :: v
      Real, Dimension(9) :: etaVO
      Real :: a,b

      a =  4. * v(6) - v(3) 

      etaVO(1) = (v(1) + 2.* v(4) + a )/3.0
      etaVO(4) = (v(2) + 2.* v(5) + a )/3.0
      etaVO(6) = 2.* v(6)

      b = 2.* (v(3) - v(6))

      etaVO(5) = - (v(2) + b + 2.*v(5))/3.0
      etaVO(3) = - (v(1) + b + 2.*v(4))/3.0
      etaVO(2) = etaVO(3) + etaVO(5) + v(3)

      etaVO(7) = .5* (v(5) + v(6)) 
      etaVO(8) = .5* (v(4) + v(6)) 
      etaVO(9) = .5* (v(4) + v(5)) 

      Return
      End 
