!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini, Ga¨el Durand
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!> Computes nodes surface loads from surface pressure
SUBROUTINE GetHydrostaticLoads( Model,Solver,dt,TransientSimulation )

!------------------------------------------------------------------------------
!******************************************************************************
!
!  
!  TODO: switch y and z as DIM=3
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
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Variable_t), POINTER :: PointerToVariable
  TYPE(Solver_t), POINTER :: PointerToSolver
  TYPE(Nodes_t), SAVE :: Nodes
  TYPE(GaussIntegrationPoints_t) :: IP

  LOGICAL :: AllocationsDone = .FALSE., GotIt, stat

  INTEGER :: i, j, n, m, t, p, Nn, istat, DIM
  INTEGER, POINTER :: Permutation(:)

  REAL(KIND=dp), POINTER :: VariableValues(:)
  REAL(KIND=dp) :: Norm, Normal(3), PwVector(3), pwi, s
  REAL(KIND=dp) :: detJ                                    

  REAL(KIND=dp), ALLOCATABLE :: pwt(:), Basis(:), dBasisdx(:,:), &
                                ddBasisddx(:,:,:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
       

  SAVE AllocationsDone, DIM, SolverName, pwt
  SAVE Basis, dBasisdx, ddBasisddx
  
  !------------------------------------------------------------------------------

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------

  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN

    DIM = CoordinateSystemDimension()
    WRITE(SolverName, '(A)') 'GetHydrostaticLoads'
    n = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
    m = Model % Mesh % NumberOfNodes
    IF (AllocationsDone) DEALLOCATE(pwt, Basis, dBasisdx, ddBasisddx)

    ALLOCATE(pwt(n), Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), STAT=istat )
         
    IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
    END IF

    AllocationsDone = .TRUE.
    CALL INFO( SolverName, 'Memory allocation done.',Level=1 )

  END IF  

  !--------------------------------------------
  ! Calculate the water force for each elements
  !--------------------------------------------

  VariableValues = 0.0_dp

  DO t = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
    n = GetElementNOFNodes()

    BC => GetBC( Element ) 
    pwt(1:n) =  -1.0 * ListGetReal(BC, 'External Pressure', n, &
                    Element % NodeIndexes , GotIt)
!
! Integration
! 
    CALL GetElementNodes( Nodes )
    IP = GaussPoints( Element )
    DO p = 1, IP % n

      stat = ElementInfo( Element, Nodes, IP % U(p), IP % V(p), &
      IP % W(p), detJ, Basis, dBasisdx, ddBasisddx, .FALSE.)          
      s = detJ * IP % S(p)                  

      Normal = NormalVector( Element, Nodes, IP % U(p), IP % V(p), .TRUE.)
!
!  Value of pwt at integration point
!
      pwi = SUM(pwt(1:n)*Basis(1:n))
!
! Compute pw_x, pw_y, pw_z
!
      PwVector(1:DIM) = pwi * Normal(1:DIM)

      DO i = 1, n
        Nn = Permutation(Element % NodeIndexes(i))
        DO j = 1, DIM
          VariableValues(DIM*(Nn-1)+j) = VariableValues(DIM*(Nn-1)+j) + PwVector(j) * s * Basis(i)
        END DO
      END DO

    END DO

  END DO

  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues )

  CALL INFO(SolverName, 'End', level=3)

!------------------------------------------------------------------------------
END SUBROUTINE GetHydrostaticLoads
!------------------------------------------------------------------------------



