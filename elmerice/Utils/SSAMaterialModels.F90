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
! *  Authors: fabien Gillet-Chaulet
! *  Email:   fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: May 2022
! *
! ******************************************************************************/
!--------------------------------------------------------------------------------
!>  Module containing utility routines for the SSA
!--------------------------------------------------------------------------------
MODULE SSAMaterialModels
    USE DefUtils

   IMPLICIT NONE
   CONTAINS

!--------------------------------------------------------------------------------
!>  Return the effective friction coefficient
!--------------------------------------------------------------------------------
   FUNCTION SSAEffectiveFriction(Element,n,Basis,ub,SEP,PartlyGrounded,h,rho,rhow,sealevel,SlipDer) RESULT(Slip)
   IMPLICIT NONE
   REAL(KIND=dp) :: Slip ! the effective friction coefficient
   TYPE(Element_t), POINTER :: Element ! the current element
   INTEGER :: n ! number of nodes
   REAL(KIND=dp) :: Basis(:) ! basis functions
   REAL(KIND=dp) :: ub  ! the velocity for non-linear friction laws
   LOGICAL :: SEP ! Sub-Element Parametrisation of the friction
   LOGICAL :: PartlyGrounded ! is the GL within the current element?
   REAL(KIND=dp) :: h ! for SEP: the ice thickness at current location
   REAL(KIND=dp) :: rho,rhow,sealevel ! density, sea-water density, sea-level
   REAL(KIND=dp),OPTIONAL :: SlipDer ! dSlip/du=dSlip/dv if ub=(u^2+v^2)^1/2 ! required to compute the Jacobian

   TYPE(ValueList_t), POINTER :: Material
   TYPE(Variable_t), POINTER :: GMSol,BedrockSol,NSol
   INTEGER, POINTER :: NodeIndexes(:)
   CHARACTER(LEN=MAX_NAME_LEN) :: Friction
   REAL(KIND=dp) :: Slip2
   REAL(KIND=dp) :: fm,fq,MinN,U0
   REAL(KIND=dp) :: alpha,beta,fB
   INTEGER :: iFriction
   INTEGER :: GLnIP

   REAL(KIND=dp),DIMENSION(n) :: NodalBeta, NodalGM, NodalBed, NodalLinVelo,NodalC,NodalN
   REAL(KIND=dp) :: bedrock,Hf,fC,fN,LinVelo

   LOGICAL :: Found

!  Sub - element GL parameterisation
   IF (SEP) THEN
     GMSol => VariableGet( CurrentModel % Variables, 'GroundedMask',UnFoundFatal=.TRUE. )
     CALL GetLocalSolution( NodalGM,UElement=Element,UVariable=GMSol)

     BedrockSol => VariableGet( CurrentModel % Variables, 'bedrock',UnFoundFatal=.TRUE. )
     CALL GetLocalSolution( NodalBed,UElement=Element,UVariable= BedrockSol)
   END IF

! Friction law
   Material => GetMaterial(Element)
   NodeIndexes => Element % NodeIndexes

   Friction = ListGetString(Material, 'SSA Friction Law',Found, UnFoundFatal=.TRUE.)
   SELECT CASE(Friction)
     CASE('linear')
       iFriction = 1
       fm = 1.0_dp
     CASE('weertman')
       iFriction = 2
     CASE('coulomb')
       iFriction = 3
     CASE('regularized coulomb')
       iFriction = 4
     CASE DEFAULT
       CALL FATAL("SSAEffectiveFriction",'Friction should be linear, Weertman, Coulomb or Regularized coulomb')
    END SELECT

    ! for all friction law
    NodalBeta = 0.0_dp
    NodalBeta(1:n) = ListGetReal( &
           Material, 'SSA Friction Parameter', n, NodeIndexes(1:n), Found,&
           UnFoundFatal=.TRUE.)

    ! for Weertman and Coulomb friction
    IF (iFriction > 1) THEN
      fm = ListGetConstReal( Material, 'SSA Friction Exponent', Found , UnFoundFatal=.TRUE.)

      NodalLinVelo = 0.0_dp
      NodalLinVelo(1:n) = ListGetReal( &
           Material, 'SSA Friction Linear Velocity', n, NodeIndexes(1:n), Found,&
           UnFoundFatal=.TRUE.)
    END IF

    ! only for Coulomb friction
    IF (iFriction == 3) THEN
      fq = ListGetConstReal( Material, 'SSA Friction Post-Peak', Found, UnFoundFatal=.TRUE. )

      NodalC = 0.0_dp
      NodalC(1:n) = ListGetReal( &
           Material, 'SSA Friction Maximum Value', n, NodeIndexes(1:n), Found,&
           UnFoundFatal=.TRUE.)

      ! Get the effective pressure
      NSol => VariableGet( CurrentModel % Variables, 'Effective Pressure', UnFoundFatal=.TRUE. )
      CALL GetLocalSolution( NodalN,UElement=Element, UVariable=NSol)

      MinN = ListGetConstReal( Material, 'SSA Min Effective Pressure', Found, UnFoundFatal=.TRUE.)
    END IF

    ! only for Regularized Coulomb friction law we need to get U0 of Joughin law
    IF (iFriction == 4) THEN
      U0 = ListGetConstReal( Material, 'SSA Friction Threshold Velocity', Found, UnFoundFatal=.TRUE.)
    END IF

    Beta=SUM(Basis(1:n)*NodalBeta(1:n))

    IF (SEP) THEN
      ! Floating
      IF (ALL(NodalGM(1:n).LT.0._dp)) THEN
        beta=0._dp
      ELSE IF (PartlyGrounded) THEN
        bedrock = SUM( NodalBed(1:n) * Basis(1:n) )
        Hf= rhow * (sealevel-bedrock) / rho
        if (h.lt.Hf) beta=0._dp
      END IF
   END IF

   Slip2=0.0_dp
   IF (iFriction > 1) THEN
     LinVelo = SUM( NodalLinVelo(1:n) * Basis(1:n) )
     IF ((iFriction == 2).AND.(fm==1.0_dp)) iFriction=1
     Slip2=1.0_dp
     IF (ub < LinVelo) then
       ub = LinVelo
       Slip2=0.0_dp
     ENDIF
   END IF

   IF (iFriction==3) THEN
     fC = SUM( NodalC(1:n) * Basis(1:n) )
     fN = SUM( NodalN(1:n) * Basis(1:n) )
     ! Effective pressure should be >0 (for the friction law)
     fN = MAX(fN, MinN)
   END IF

   SELECT CASE (iFriction)
    CASE(1)
     Slip = beta
     IF (PRESENT(SlipDer)) SlipDer = 0._dp
    CASE(2)
     Slip = beta * ub**(fm-1.0_dp)
     IF (PRESENT(SlipDer)) SlipDer = Slip2*Slip*(fm-1.0_dp)/(ub*ub)
    CASE(3)
     IF (fq.NE.1.0_dp) THEN
       alpha = (fq-1.0_dp)**(fq-1.0_dp) / fq**fq
     ELSE
       alpha = 1.0_dp
     END IF
     fB = alpha * (beta / (fC*fN))**(fq/fm)
     Slip = beta * ub**(fm-1.0_dp) / (1.0_dp + fB * ub**fq)**fm
     IF (PRESENT(SlipDer)) SlipDer  = Slip2 * Slip * ((fm-1.0_dp) / (ub*ub) - &
             fm*fq*fB*ub**(fq-2.0_dp)/(1.0_dp+fB*ub**fq))
    CASE(4)
     Slip = beta * ub**(fm-1.0_dp) / (ub + U0)**fm
     IF (PRESENT(SlipDer)) SlipDer = Slip2 * Slip * ((fm-1.0_dp) / (ub*ub) - &
              fm*ub**(-1.0_dp)/(ub+U0))
   END SELECT

  END FUNCTION SSAEffectiveFriction

  END MODULE SSAMaterialModels

