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
! *  Module for solving elementsize into a nodal field.
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Solves a simple equation for the elementsize using the Galerkin method. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ElementSizeSolver( Model,Solver,dt,TransientSimulation )
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

  INTEGER :: n, t, istat, active
  REAL(KIND=dp) :: Norm, Power
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
  LOGICAL :: GotIt, NoWeight
  REAL(KIND=dp) :: ElemMin, ElemMax
  TYPE(ValueList_t), POINTER :: Params


  SAVE STIFF, FORCE, AllocationsDone
!------------------------------------------------------------------------------

  CALL Info('ElementSizeSolver','Computing nodal element size indicator')

  Mesh => GetMesh()
  Params => GetSolverParams()


  ! Allocate some storage
  !--------------------------------------------------------------
  N = Mesh % MaxElementNodes ! just big enough for elemental arrays
  ALLOCATE( FORCE(N), STIFF(N,N), STAT=istat )
  IF ( istat /= 0 ) THEN
    CALL Fatal( 'ElementSizeSolver', 'Memory allocation error.' )
  END IF
  
  ElemMin = HUGE( ElemMin ) 
  ElemMax = -HUGE( ElemMax ) 

  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()
  
  Power = 1.0_dp / ListGetCReal( Solver % Values,'Element Size Exponent',GotIt)
  IF(.NOT. GotIt ) Power = 1.0_dp
  NoWeight = ListGetLogical( Solver % Values,'No Integration Weight',GotIt)
  
  active = GetNOFActive()
  DO t=1,active
    Element => GetActiveElement(t)
    n = GetElementNOFNodes()
    
    !Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(  STIFF, FORCE, Element, n )
    
    !Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  CALL DefaultFinishBulkAssembly()
  
  ! No flux BCs
  CALL DefaultFinishAssembly()

  CALL DefaultDirichletBCs()
  Norm = DefaultSolve()
  
  DEALLOCATE( FORCE, STIFF )

  WRITE(Message,'(A,ES12.4)') 'Minimum Element Size: ',ElemMin
  CALL Info('ElementSizeSolver',Message)
  WRITE(Message,'(A,ES12.4)') 'Maximum Element Size: ',ElemMax
  CALL Info('ElementSizeSolver',Message)
  WRITE(Message,'(A,ES12.4)') 'Element Size Ratio: ',ElemMax / ElemMin
  CALL Info('ElementSizeSolver',Message)
 
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),DetJ,LoadAtIP,Weight
    LOGICAL :: Stat
    INTEGER :: i,j,t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    ! Numerical integration:
    !----------------------
    IP = GaussPoints( Element )

    DO t=1,IP % n
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t),  detJ, Basis )

       ! The source term at the integration point:
       !------------------------------------------
       LoadAtIP = DetJ ** Power
       IF( NoWeight ) THEN
         Weight = IP % s(t) 
       ELSE
         Weight = IP % s(t) * DetJ
       END IF

       ElemMin = MIN( ElemMin, LoadAtIP )
       ElemMax = MAX( ElemMax, LoadAtIP )
       
       ! Finally, the elemental matrix & vector:
       !----------------------------------------
       DO i = 1, n
         DO j = 1, n
           STIFF(i,j) = STIFF(i,j) + Weight * &
               Basis(i) * Basis(j)
         END DO
         FORCE(i) = FORCE(i) + Weight * Basis(i) * LoadAtIp
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------
END SUBROUTINE ElementSizeSolver
!------------------------------------------------------------------------------
