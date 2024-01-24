!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Olivier Gagliardini and Julien Brondex
! *  Email:   olivier.gagliardini@univ-grenoble-alpes.fr
! *           julien.brondex@univ-grenoble-alpes.fr
! *           
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: January 2024
! *
! ******************************************************************************/
!--------------------------------------------------------------------------------
!>  Module containing utility routines for the Porous Solver
!--------------------------------------------------------------------------------

MODULE PorousMaterialModels

   USE DefUtils

   IMPLICIT NONE

   CONTAINS

           
!-------------------------------------------------------------------------------           
!> Function a(D) and b(D) from Gagliardini and Meyssonier, 1997.
!> modified to fulfill the condition 3x a(D) >= 2x b(D) for D > 0.1
!> for n=3 only
!-------------------------------------------------------------------------------

!> Function a (Olivier Gagliardini)
   FUNCTION ParameterA ( D ) RESULT(a)
      USE types
      USE CoordinateSystems
      USE SolverUtils
      USE ElementDescription
      USE DefUtils
      IMPLICIT NONE
      REAL(KIND=dp) :: a, D, DD     

      IF (D > 1.0_dp - AEPS) THEN
         a = 1.0_dp           

      ELSE IF (D <= 0.81) THEN
         DD = D
         IF (D < 0.4 ) DD = 0.4_dp
         a = EXP( 13.22240_dp - 15.78652_dp * DD )

      ELSE
         a =  1.0_dp  + 2.0/3.0 * (1.0_dp - D) 
         a = a / ( D**1.5 )
      END IF
   END FUNCTION ParameterA 


!> Function b (Olivier Gagliardini)
   FUNCTION ParameterB ( D ) RESULT(b)
      USE types
      USE CoordinateSystems
      USE SolverUtils
      USE ElementDescription
      USE DefUtils
      IMPLICIT NONE
      REAL(KIND=dp) :: b, D, DD 

      IF (D > 1.0_dp - AEPS) THEN
         b = 0.0_dp           

      ELSE IF (D <= 0.81) THEN
         DD = D
         IF (D < 0.4 ) DD = 0.4_dp
         b = EXP( 15.09371_dp - 20.46489_dp * DD )

      ELSE
         b = (1.0_dp - D)**(1.0/3.0) 
         b = 3.0/4.0 * ( b / (3.0 * (1 - b)) )**1.5

      END IF
   END FUNCTION ParameterB 

!------------------------------------------------------------------------------
!> Returns effective viscosity for Porous Solver (Julien Brondex) 
!------------------------------------------------------------------------------
   FUNCTION PorousEffectiveViscosity( Fluidity, Density, SR, Em, Element, &
        Nodes, n, nd, u, v, w, LocalIP ) RESULT(mu)
     !------------------------------------------------------------------------------

     REAL(KIND=dp)  :: Fluidity, Density, u, v, w, eta, Kcp, mu(2)
     REAL(KIND=dp) ::  SR(:,:)
     TYPE(Nodes_t)  :: Nodes
     INTEGER :: n,nd
     INTEGER, OPTIONAL :: LocalIP
     TYPE(Element_t),POINTER :: Element

     !------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: Material
     REAL(KIND=dp) :: Wn, MinSRInvariant, fa, fb, ss, nn, Em
  
     INTEGER :: i, j
     INTEGER :: body_id
     INTEGER :: old_body = -1

     LOGICAL :: GotIt

     SAVE :: old_body, Wn, MinSRInvariant 

     Material => GetMaterial(Element)

!    Get uniform rheological parameters only once per body to save computational time
!-------------------------------------------------------------------------------------
     body_id = Element % BodyId
     IF (body_id /= old_body) Then 
        old_body = body_id
        ! Get flow law exponent
        Wn = GetConstReal( Material , 'Powerlaw Exponent', GotIt )
        IF (.NOT.GotIt) THEN
           CALL INFO('ComputeDevStress', 'Variable  Powerlaw Exponent not found. &
                          & Setting to 1.0', Level = 20)
           Wn = 1.0
        ELSE
           WRITE(Message,'(A,F10.4)') 'Powerlaw Exponent = ',   Wn
           CALL INFO('ComputeDevStress', Message, Level = 20)
        END IF
        ! Get the Minimum value of the Effective Strain rate 
        MinSRInvariant = 100.0*AEPS
        IF ( Wn > 1.0 ) THEN
           MinSRInvariant = GetConstReal( Material, 'Min Second Invariant', GotIt )
           IF (.NOT.GotIt) THEN
              WRITE(Message,'(A)') 'Variable Min Second Invariant not &
                       &found. Setting to 100.0*AEPS )'
              CALL INFO('ComputeDevStress', Message, Level = 20)
           ELSE
              WRITE(Message,'(A,F14.8)') 'Min Second Invariant = ', MinSRInvariant
              CALL INFO('ComputeDevStress', Message, Level = 20)
           END IF
        END IF
     END IF

!    Get non-uniform rheological parameters for Porous at current point
!-------------------------------------------------------------------------------------
!    Evaluate parameterized functions a and b from Density interpolated at current point
     fa = ParameterA(Density)
     fb = ParameterB(Density)

!    Calculate Strain Rate invariant as E_D^2 = gamma_e^2/fa + E_m^2/fb
     ss = 1.0_dp !case linear by default -> ss = E_D^((1-n)/n) -> ss = 1 if n = 1
     IF ( Wn > 1.0 ) THEN !case non-linear
        ss = 0.0_dp !Initialize ss
        DO i = 1, 3
           DO j = 1, 3
              ss = ss + SR(i,j)**2 ! Gamma_e^2/2 = e_ij e_ij
           END DO
        END DO
        ss = 2.0*ss / fa !Gamma_e^2/fa
        IF ( fb > 1.0e-8 ) ss = ss + Em**2 / fb  !Full E_D^2  = gamma_e^2/fa + E_m^2/fb

        nn =(1.0 - Wn)/Wn !nn = (1 -n)/n

        ss = SQRT(ss)! Full E_D
        IF (ss < MinSRInvariant ) ss = MinSRInvariant !Invariant not less than prescribed min
        ss =  ss**nn ! ss = E_D^((1-n)/n)
     END IF

!    Compute effective viscosity at current point
!-------------------------------------------------------------------------------------
     eta = ss / ( fa * Fluidity ** (1.0 / Wn ) )
!    Compute compressibility parameter  at current point
!-------------------------------------------------------------------------------------
     Kcp = fb * Fluidity ** (1.0 / Wn ) / ss
!    Return output array as [eta, Kcp]
!-------------------------------------------------------------------------------------
     mu = [ eta, Kcp ]
!------------------------------------------------------------------------------
   END FUNCTION PorousEffectiveViscosity        
!------------------------------------------------------------------------------

END MODULE PorousMaterialModels
