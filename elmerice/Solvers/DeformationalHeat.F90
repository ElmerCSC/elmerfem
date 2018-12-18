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
! *  Authors: 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!> DOXYGEN INFORMATION TO BE ADDED
SUBROUTINE DeformationalHeatSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
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
  TYPE(ValueList_t), POINTER :: SolverParams, Material
  TYPE(Variable_t), POINTER :: PointerToVariable,FlowSol
  TYPE(Solver_t), POINTER :: PointerToSolver

  LOGICAL :: AllocationsDone = .FALSE., Found,UnFoundFatal=.TRUE.

  INTEGER :: i, j,n, m, t, istat,k
  INTEGER, POINTER :: Permutation(:), FlowPerm(:), NodeIndexes(:)

  REAL(KIND=dp), POINTER :: VariableValues(:), FlowSolution(:) 
  REAL(KIND=dp) :: Norm
  Integer :: STDOFs,NSDOFs,dim

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), Velo(:,:), Viscosity(:)
  
  CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolName,SolverName

  SAVE STIFF, LOAD, FORCE,  Velo, AllocationsDone, Viscosity
!------------------------------------------------------------------------------
  SolverName = 'Deformational Heat Solver'
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  STDOFs = PointerToVariable % DOFs

  dim = CoordinateSystemDimension()
  
  FlowSolName  = GetString( Solver % Values, 'Flow Solver Name', Found )
  IF (.NOT.Found)   WRITE(FlowSolName,'(A)') 'Flow Solution'
  FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName,UnFoundFatal=UnFoundFatal)
  FlowPerm     => FlowSol % Perm
  NSDOFs     =  FlowSol % DOFs
  FlowSolution => FlowSol % Values

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, Viscosity, Velo)

     ALLOCATE( FORCE(2*STDOFs*N), LOAD(2*STDOFs*N), STIFF(2*STDOFs*N,2*STDOFs*N), Viscosity(N), Velo(3,N), &
          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'HessianSolve', 'Memory allocation error.' )
     END IF
     

     AllocationsDone = .TRUE.
  END IF

  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()

  ! bulk assembly
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     NodeIndexes => Element % NodeIndexes

     Material => GetMaterial()
     Viscosity(1:n) = GetReal( Material,'Viscosity', Found )
     IF (.NOT.Found) CALL FATAL(SolverName,'Could not find  >Viscosity<')


     Velo = 0.0d0
     Do i=1,n
        j = NSDOFs*FlowPerm(NodeIndexes(i))
        Do k=1,dim
          Velo(k,i) =  FlowSolution( j-dim+k-1 )
        End do
     End do

     CALL LocalMatrix(  STIFF, FORCE, Element, n, Velo, Viscosity )
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  

  CALL DefaultFinishAssembly()

  Norm = DefaultSolve()

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element,  n, Velo, Viscosity)
!------------------------------------------------------------------------------
    USE MaterialModels

    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), Velo(:,:), Viscosity(:)
    INTEGER :: n,dim
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3),DetJ,LGrad(3,3)
    REAL(KIND=dp) :: mu,VeloIP(3)
    LOGICAL :: Stat
    INTEGER :: t, p,q , i 
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0
    LGrad = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n

       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )
       


        mu = SUM( Viscosity(1:n) * Basis(1:n) )
        mu = EffectiveViscosity( mu, 1.0_dp , Velo(1,:) , Velo(2,:), Velo(3,:), &
                                   Element, Nodes, n, n, IP % U(t), IP % V(t), IP % W(t), LocalIP=t )

        LGrad = MATMUL( Velo(:,1:n), dBasisdx(1:n,:) )
        VeloIP=0.
        VeloIP(1) = SUM( Velo(1,1:n)*Basis(1:n) )
        VeloIP(2) = SUM( Velo(2,1:n)*Basis(1:n) )
        IF ( dim > 2 ) VeloIP(3) = SUM(Velo(3,1:n)*Basis(1:n) )


            DO p=1,n
              DO q=1,n
                  STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * Basis(p) * Basis(q)
               End do
          END DO
         DO p=1,n
            Force(p) = Force(p) + IP % S(t) * detJ * 0.5_dp * mu * SecondInvariant(VeloIP,LGrad)  * Basis(p)
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE DeformationalHeatSolver
!------------------------------------------------------------------------------

