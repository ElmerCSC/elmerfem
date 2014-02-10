!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - Scientific Computing Ltd., Finland
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
!/******************************************************************************
! *
! *  ExportVertically Solver to export vertically a variable computed on 
! *  the upper or lower boundary in the whole domain 
! *
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini
! *  Email:   gagliar@lgge.obs.ujf-grenoble.fr
! *  Web:     http://www.csc.fi/elmer
! *
! *  Original Date: 30 April 2010 
! * 
! *****************************************************************************


! *****************************************************************************
SUBROUTINE ExportVertically( Model,Solver,dt,TransientSimulation )
!DEC$ATTRIBUTES DLLEXPORT :: ExportVertically
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Export vertically the SSABasal Velocity (given as a Dirichlet Boundary condition) 
!  Compute also the vertical velocity and the pressure
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
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
  TYPE(Element_t),POINTER :: CurrentElement, Element
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material
  TYPE(Variable_t), POINTER :: PointerToVariable

  LOGICAL :: AllocationsDone = .FALSE.

  INTEGER :: i, n, m, t, istat, p, Indexes(128) 
  INTEGER, POINTER :: Permutation(:) 

  REAL(KIND=dp), POINTER :: ForceVector(:), VariableValues(:)
  REAL(KIND=dp) :: Norm

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
       

  SAVE STIFF, FORCE, AllocationsDone, SolverName
!------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  WRITE(SolverName, '(A)') 'ExportVertically'


  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, STIFF)

     ALLOCATE( FORCE(N), STIFF(N,N), STAT=istat )
          
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF

     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS


! No non-linear iteration, no time dependency  
  VariableValues = 0.0d0
  Norm = Solver % Variable % Norm


  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()
  ! bulk assembly
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes()


     CALL LocalMatrix (  STIFF, FORCE, Element, n ) 

     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  
  CALL DefaultFinishAssembly()

  ! Dirichlet 
  CALL DefaultDirichletBCs()
  
  !Solve the system
  Norm = DefaultSolve()


CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n ) 
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ 
    LOGICAL :: Stat
    INTEGER :: t, p, q , dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()


    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       DO p=1,n
         DO q=1,n
           STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
         END DO
       END DO

    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE ExportVertically
!------------------------------------------------------------------------------


