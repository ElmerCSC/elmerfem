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

!------------------------------------------------------------------------------
!>  Solves the energy release rate for crack propagation.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ReleaseRateSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element

  LOGICAL :: AllocationsDone = .FALSE., Found

  INTEGER :: n, t, istat
  REAL(KIND=dp) :: Norm

  REAL(KIND=dp), ALLOCATABLE :: LocalDisplacement(:,:)
  REAL(KIND=dp), ALLOCATABLE :: LocalStress(:,:)
  REAL(KIND=dp), ALLOCATABLE :: LocalPropagationShape(:,:)

  REAL(KIND=dp) :: LocalGtheta, Gtheta

  SAVE LocalDisplacement, LocalStress, LocalPropagationShape, AllocationsDone
!------------------------------------------------------------------------------

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough
     ALLOCATE(  LocalDisplacement(3,N), LocalStress(6,N), &
          LocalPropagationShape(3,N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'PoissonSolve', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF

  ! Calculate the energy release rate:
  ! ---------------------------------
   Gtheta = 0.0d0

   DO t = 1, Solver % NumberOfActiveElements
      Element => GetActiveElement(t)
      n = GetElementNOFNodes()
      
      CALL GetVectorLocalSolution( LocalDisplacement, 'True Displacement' )
      CALL GetVectorLocalSolution( LocalStress, 'Stress' )
      CALL GetVectorLocalSolution( LocalPropagationShape, 'Shape Displacement')

      CALL LocalReleaseRate( LocalGtheta, LocalPropagationShape, &
           LocalDisplacement, LocalStress, n )

      Gtheta = Gtheta + LocalGtheta
   END DO

   PRINT *
   PRINT *,'*******************************'
   PRINT *,'Calculated Fracture Parameters:'
   PRINT *,'-------------------------------'
   PRINT *,'2G = ',Gtheta
   PRINT *,'*******************************'
   PRINT *

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalReleaseRate( LocalGtheta, LocalPropagationShape, &
       LocalDisplacement, LocalStress, n )
!------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: LocalGtheta
    REAL(KIND=dp) :: LocalPropagationShape(:,:)
    REAL(KIND=dp) :: LocalDisplacement(:,:), LocalStress(:,:)
    INTEGER :: n

    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,LoadAtIP
    INTEGER :: t, i, j
    LOGICAL :: stat

    REAL(KIND=dp) :: Displacement(3), GradDisplacement(3,3)
    REAL(KIND=dp) :: PropagationShape(3), GradPropagationShape(3,3)
    REAL(KIND=dp) :: r_tensor(3,3), s_tensor(3,3), StressTensor(3,3)
    REAL(KIND=dp) :: DivPropagationShape

    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    LocalGtheta = 0.0d0

    IP = GaussPoints( Element )
    
    ! Loop over integration points:
    !------------------------------
    DO t = 1, IP % n

       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx )

       DO i = 1,3
          Displacement(i) = SUM( Basis(1:n) * LocalDisplacement(i,1:n) )
          PropagationShape(i) = SUM( Basis(1:n) * LocalPropagationShape(i,1:n) )
          
          DO j = 1,3
             GradDisplacement(i,j) = SUM( dBasisdx(1:n,j) * LocalDisplacement(i,1:n) )
             GradPropagationShape(i,j) = SUM( dBasisdx(1:n,j) * LocalPropagationShape(i,1:n) )
          END DO
       END DO
       
       DivPropagationShape = 0.0d0
       DO i = 1,3
          DivPropagationShape = DivPropagationShape + GradPropagationShape(i,i)
       END DO

       StressTensor(1,1) = SUM( Basis(1:n) * LocalStress(1,1:n) )
       StressTensor(2,2) = SUM( Basis(1:n) * LocalStress(2,1:n) )
       StressTensor(3,3) = SUM( Basis(1:n) * LocalStress(3,1:n) )
       StressTensor(1,2) = SUM( Basis(1:n) * LocalStress(4,1:n) )
       StressTensor(2,3) = SUM( Basis(1:n) * LocalStress(5,1:n) )
       StressTensor(1,3) = SUM( Basis(1:n) * LocalStress(6,1:n) )
       StressTensor(2,1) = SUM( Basis(1:n) * LocalStress(4,1:n) )
       StressTensor(3,2) = SUM( Basis(1:n) * LocalStress(5,1:n) )
       StressTensor(3,1) = SUM( Basis(1:n) * LocalStress(6,1:n) )

       r_tensor = -MATMUL( GradDisplacement, GradPropagationShape )

       s_tensor = MATMUL( StressTensor, TRANSPOSE( GradPropagationShape ) ) &
            - DivPropagationShape * StressTensor

       LocalGtheta = LocalGtheta + ( SUM( s_tensor * GradDisplacement)  &
            - SUM( r_tensor * StressTensor ) ) * IP % s(t) * detJ
       
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalReleaseRate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE ReleaseRateSolver
!------------------------------------------------------------------------------
