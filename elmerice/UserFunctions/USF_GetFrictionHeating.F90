!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Module containing a functions for friction heat
! *
! ******************************************************************************
! *
! *  Authors: Thomas Zwinger, Juha Ruokolainen, Hakime Seddik, Joe Todd
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *  Current date: 20 July 2012 (Martina)
! *
! *****************************************************************************/


! contains function for adding friction heat to geothermal heat flux at the basis and exporting it
! ATTENTION: consistent units for capacity and conductivy are needed.
! ATTENTION: The keyword "Friction Heat" does make reference to Strainheating, not Friction heat!
! usage:


! Solver 1
!   Equation = "Normal vector"
!   Variable = "Normal Vector"   
!   ! in 3dimensional simulations we have 3 entries
!   Variable DOFs = 3 
!   !NB: does not need to actually solve a matrix
!   !    hence no BW optimization needed
!   Optimize Bandwidth = Logical False 
!   Procedure = "ComputeNormal" "ComputeNormalSolver"
!   ! if set to True, all boundary normals would be computed by default
!   ComputeAll = Logical False
!End
!Solver 2
!  Equation = String "StressSolver"
!  Procedure =  File "ComputeDevStressNS" "ComputeDevStress"
!  ! this is just a dummy, hence no output is needed
!  !-----------------------------------------------------------------------
!  Variable = -nooutput "Sij"
!  Variable DOFs = 1
!  ! the name of the variable containing the flow solution (U,V,W,Pressure)
!  !-----------------------------------------------------------------------
!  Flow Solver Name = String "Flow Solution"
!  Exported Variable 1 = "Stress" ! [Sxx, Syy, Szz, Sxy] in 2D
!                                 ! [Sxx, Syy, Szz, Sxy, Syz, Szx] in 3D
!  Exported Variable 1 DOFs = 6   ! 4 in 2D, 6 in 3D
!  Linear System Solver = "Iterative"
!  Linear System Iterative Method = "BiCGStab"
!  Linear System Max Iterations = 300
!  Linear System Convergence Tolerance = 1.0E-09
!  Linear System Abort Not Converged = True
!  Linear System Preconditioning = "ILU0"
!  Linear System Residual Output = 1
!End
!Solver 3
!  Equation = "Dummy"
!  Procedure = File "DummySolver" "DummySolver"
!  Variable = String "friction"
!  Variable DOFs = 1
!End

!Material 1
! ...
!   Cauchy = Logical True !needed for ComputeDevStressNS
! ...
!End

!Boundary Condition 1
!  Name = "bed"
!  ComputeNormal = Logical True !needed for ComputeNormal
! ....
!  Temp Flux BC = Logical True
!  friction = Variable Coordinate 1
!       Real Procedure "./USF_GetFrictionHeating" "getFrictionHeat"
!  Temp Heat Flux = Variable friction
!      Real MATC "0.063*(31556926.0)*1.0E-06 +tx" !assuming 63mW/m2 as geothermal heatflux
!End

!-----------
!function to export and define FrictionHeat at bedrock

FUNCTION getFrictionHeat(  Model, Node, DummyInput)RESULT(frictionheat)
  
  USE DefUtils

  IMPLICIT NONE
  
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: DummyInput, frictionHeat
  !----------------------------------------------------------------------------
  
  INTEGER :: DIM, i, j, Ind(3,3)
  REAL(KIND=dp), POINTER :: FlowValues(:),NormalValues(:),StressValues(:)
  REAL(KIND=dp) :: normal(3), velo(3), un, ut, Sig(3,3), Sn(3), snn, snt
  INTEGER, POINTER :: FlowPerm(:),StressPerm(:), NormalPerm(:)
  LOGICAL :: FirstTime=.TRUE.
  TYPE(Variable_t), POINTER :: FlowSol,StressVariable, NormalVar
  CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName
  
  SAVE FirstTime, DIM, FunctionName,Ind
  
  IF (FirstTime) THEN
     WRITE(FunctionName, '(A)') 'getFrictionHeat'
     DIM = CoordinateSystemDimension()
     FirstTime = .FALSE.
     DO i=1, 3
        Ind(i,i) = i
     END DO
     Ind(1,2) = 4
     Ind(2,1) = 4
     Ind(2,3) = 5
     Ind(3,2) = 5
     Ind(3,1) = 6
     Ind(1,3) = 6
  END IF
  
  ! Get the variable velocity
  !---------------------------
  FlowSol => VariableGet( Model % Variables, 'Flow Solution' )
  IF ( ASSOCIATED( FlowSol ) ) THEN
     FlowPerm    => FlowSol % Perm
     FlowValues  => FlowSol % Values
  ELSE
     CALL FATAL(FunctionName, 'Need NS Solver, Flow Solution not found')
  END IF
  
  ! Get the stress variable
  !------------------------
  StressVariable => VariableGet( Model % Variables, 'Stress' )
  IF ( ASSOCIATED( StressVariable ) ) THEN
     StressPerm    => StressVariable % Perm
     StressValues  => StressVariable % Values
  ELSE
     CALL FATAL(FunctionName,'No variable Stress found')   
  END IF
  
  ! Get the variable for normal vector
  NormalVar =>  VariableGet(Model % Variables,'Normal Vector')
  IF ( ASSOCIATED( NormalVar ) ) THEN
     NormalPerm => NormalVar % Perm
     NormalValues => NormalVar % Values
  ELSE
     CALL FATAL(FunctionName, 'Normal Vector variable not found')
  END IF
  
  DO i=1, DIM
     normal(i) = NormalValues(DIM*(NormalPerm(Node)-1) + i)      
     velo(i) = FlowValues( (DIM+1)*(FlowPerm(Node)-1) + i )
  END DO
  
  !Tangential velocity
  un = SUM(velo(1:DIM)*(normal(1:DIM))) 
  ut = SQRT(SUM( (velo(1:DIM))**2.0_dp ) - (un**2.0_dp) )

  !Tangential Stress
  DO i=1, DIM
     DO j= 1, DIM
        Sig(i,j) =  &
             StressValues( 2*DIM *(StressPerm(Node)-1) + Ind(i,j) )
     END DO
  END DO
  DO i=1, DIM
     Sn(i) = SUM(Sig(i,1:DIM)*normal(1:DIM)) 
  END DO
  Snn = SUM( Sn(1:DIM) * normal(1:DIM) )
  Snt = SQRT( SUM(Sn(1:DIM)**2.0_dp) - (Snn**2.0_dp))

  frictionHeat =ut*Snt
END FUNCTION getFrictionHeat



