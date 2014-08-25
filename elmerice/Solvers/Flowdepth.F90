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
! *  Authors: Thomas Zwinger, Martina Schäfer
! *  Email:   Thomas.Zwinger@csc.fi
! *  Web:      http://elmerice.elmerfem.org
! *  Address: CSC - Scientific Computing Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 
! *
! ****************************************************************************/
!> FlowDepthSolver: Solver to inquire the vertical distance from a line (2d) 
!>   or surface (3d). Solves a degenrated Poisson-equation
SUBROUTINE FlowdepthSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the Flowdepth equation!
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
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: SolverParams, BC
  TYPE(Variable_t), POINTER :: PointerToVariable, Grad1Sol, Grad2Sol, SurfSol
  TYPE(Solver_t), POINTER :: PointerToSolver

  LOGICAL :: AllocationsDone = .FALSE., Found, CalcFree = .FALSE., GotIt, SkipBoundary

  INTEGER :: i, n, m, t, istat, DIM
  INTEGER, POINTER :: Permutation(:), NumberOfVisits(:),&
       SurfacePerm(:), GradSurface1Perm(:),GradSurface2Perm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:), Surface(:), GradSurface1(:),GradSurface2(:)
  REAL(KIND=dp) :: Norm, Gradient,GradSurface(3)

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, FreeSurfName, FreeSurfGradName
       

  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName, NumberOfVisits
  !------------------------------------------------------------------------------
  PointerToSolver => Solver
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     DIM = CoordinateSystemDimension()
     WRITE(SolverName, '(A)') 'FlowdepthSolver'
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, NumberOfVisits)

     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), NumberOfVisits(M),&
          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF

     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

  !----------------------------------
  ! get Gradient (change of property
  ! with respect to unit-length of
  ! vertical direction
  !----------------------------------
  SolverParams => GetSolverParams()
  Gradient = GetConstReal( SolverParams, &
                      'Gradient',  Found )
  IF (.NOT. Found) THEN
     CALL WARN(SolverName, 'No keyword >Gradient< found in section Solver')
     CALL WARN(SolverName, 'Assuming value of -1')
     Gradient = -1.0D00
  ELSE
      WRITE(Message,'(A,e12.4,A)') 'Gradient of ',Gradient,' applied'
     CALL INFO(SolverName, Message,Level=1)
  END IF
  !----------------------------------------------------------------
  ! Assign Variables for computations fo free surface and gradients
  !----------------------------------------------------------------
  CalcFree = GetLogical(SolverParams, 'Calc Free Surface', Found)
  IF (.NOT.Found) THEN
     CalcFree = .FALSE.
  ELSE
     IF (CalcFree) THEN
               CALL INFO(SolverName, 'Free surface variable will be calculated', Level=1)
        FreeSurfName = GetString( Solver % Values,'Freesurf Name',GotIt)
        IF (.NOT.GotIt) THEN
           CALL FATAL(SolverName,'Keyaword >Calc Free Surface< set to true, but keyword >Freesurf Name< not found.')
        END IF

        SurfSol => VariableGet( Solver % Mesh % Variables, TRIM(FreeSurfname))
        IF (ASSOCIATED(SurfSol)) THEN
           Surface => SurfSol % Values
           SurfacePerm => SurfSol % Perm
        ELSE
           ALLOCATE(Surface(M), STAT=istat )
           SurfacePerm => Permutation
           IF ( istat /= 0 ) THEN
              CALL Fatal( SolverName, 'Memory allocation error.' )
           END IF
           CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, PointerToSolver, &
                TRIM(FreeSurfname), 1, Surface, SurfacePerm)
        END IF

        FreeSurfGradName = TRIM(FreeSurfname) // 'Grad1'
        Grad1Sol => VariableGet( Solver % Mesh % Variables, FreeSurfGradName)
        IF (ASSOCIATED(Grad1Sol)) THEN
           GradSurface1 => Grad1Sol % Values
           GradSurface1Perm => Grad1Sol % Perm
        ELSE
           ALLOCATE(GradSurface1(M), STAT=istat )
           GradSurface1Perm => Permutation
           IF ( istat /= 0 ) THEN
              CALL Fatal( SolverName, 'Memory allocation error.' )
           END IF
           CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, PointerToSolver, &
                TRIM(FreeSurfGradName), 1, GradSurface1, GradSurface1Perm)
        END IF
        FreeSurfGradName = TRIM(FreeSurfname) // 'Grad2'
        Grad2Sol => VariableGet( Solver % Mesh % Variables, FreeSurfGradName)
        IF (ASSOCIATED(Grad2Sol)) THEN
           GradSurface2 => Grad2Sol % Values
           GradSurface2Perm => Grad2Sol % Perm
        ELSE
           ALLOCATE(GradSurface2(M), STAT=istat )
           GradSurface2Perm => Permutation
           IF ( istat /= 0 ) THEN
              CALL Fatal( SolverName, 'Memory allocation error.' )
           END IF
           CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, PointerToSolver, &
                TRIM(FreeSurfGradName), 1, GradSurface2, GradSurface2Perm)
        END IF
        !---------------
        ! initialization
        !---------------
        Surface = 1.0D00
        GradSurface1 = 2.0D00
        GradSurface2 = 3.0D00
     END IF
  END IF
  

  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()
  ! bulk assembly
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes()
     CALL LocalMatrix(  STIFF, FORCE, Element, n)
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  
  ! Neumann conditions
  DO t=1,Solver % Mesh % NUmberOfBoundaryElements
     Element => GetBoundaryElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     BC => GetBC(Element)
     SkipBoundary = GetLogical(BC, &
                      'Flowdepth Skip',  Found )
     IF (.NOT.Found) SkipBoundary = .FALSE.
     IF (SkipBoundary) CYCLE
     n = GetElementNOFNodes()
     STIFF = 0.0D00
     FORCE = 0.0D00
     ! only select elements with defined normals with respect to the dimension of the problem
     IF( GetElementFamily() .GE. DIM) &
          CALL LocalMatrixBC(  STIFF, FORCE, Element, n, Gradient)
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO

  CALL DefaultFinishAssembly()

  ! Dirichlet 
  CALL DefaultDirichletBCs()

  !Solve the system
  Norm = DefaultSolve()

  !--------------------------------------------------------
  ! post-processing steps for free surface and its gradient
  !--------------------------------------------------------
  IF (Calcfree) THEN   
     Surface = 0.0D00
     GradSurface1 = 0.0D00
     GradSurface2 = 0.0D00
     NumberOfVisits = 0
     DO t=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        CALL GetSurfaceValue(Model, Surface, GradSurface1, GradSurface2,&
             VariableValues, Permutation, &
             SurfacePerm, GradSurface1Perm, GradSurface2Perm, &
             NumberOfVisits, Element, n, Gradient )
     END DO
     DO i=1,Model % Mesh % NumberOfNodes
        GradSurface1(GradSurface1Perm(i)) = GradSurface1(GradSurface1Perm(i))/NumberOfVisits(i)
        GradSurface2(GradSurface2Perm(i)) = GradSurface2(GradSurface2Perm(i))/NumberOfVisits(i)
     END DO
  END IF
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE GetSurfaceValue(Model, Surface, GradSurface1, GradSurface2, &
       VariableValues, Permutation, &
       SurfacePerm, GradSurface1Perm, GradSurface2Perm, &
       NumberOfVisits, Element, n, Gradient )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: LocalSurf, GradSurface(3), Depth, Gradient 
    INTEGER :: n
    INTEGER, POINTER :: Permutation(:), NumberOfVisits(:), &
         SurfacePerm(:), GradSurface1Perm(:), GradSurface2Perm(:)
    REAL(KIND=dp), POINTER :: Surface(:), GradSurface1(:), GradSurface2(:), VariableValues(:)
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3),DetJ,&
         z, U, V, W, SqrtElementMetric
    LOGICAL :: Stat
    INTEGER :: i,j,k,dim
    LOGICAL :: FirstTime
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    dim = CoordinateSystemDimension()
    ! loop over all nodes in Element
    DO i=1,N

       j = Element % NodeIndexes(i) ! get number of node in element in physical space
       NumberOfVisits(j) = NumberOfVisits(j) + 1

       ! get local coordinates of the point i inside the element
       U = Element % Type % NodeU(i)
       V = Element % Type % NodeV(i)
       W = Element % Type % NodeW(i)

       ! get local information on test-functions and derivatives of the point i
       stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
            Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE. )  
       IF (DIM == 2) THEN 
          z = Model % Nodes % y(j)
       ELSE IF (DIM == 3) THEN
          z = Model % Nodes % z(j)
       ELSE
          CALL FATAL(SolverName, 'Flow depth for one-dimensional problem not defined!')
       END IF

       IF (NumberOfVisits(j) == 1) &
            Surface(SurfacePerm(j)) = z -VariableValues(Permutation(j))/Gradient  

       GradSurface1(GradSurface1Perm(j)) = GradSurface1(GradSurface1Perm(j)) +&
            SUM(dBasisdx(1:N,1)*VariableValues(Permutation(Element % NodeIndexes(1:N))))/abs(Gradient)

       IF (DIM > 2) &
            GradSurface2(GradSurface1Perm(j)) = GradSurface2(GradSurface1Perm(j)) +&
            SUM(dBasisdx(1:N,2)*VariableValues(Permutation(Element % NodeIndexes(1:N))))/abs(Gradient)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE GetSurfaceValue
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3),DetJ
    LOGICAL :: Stat
    INTEGER :: t, p,q ,dim
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
  SUBROUTINE LocalMatrixBC(  STIFF, FORCE, Element, n, Gradient)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), Gradient
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
                      DetJ,Normal(3)
    LOGICAL :: Stat
    INTEGER :: t, dim
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

      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      FORCE(1:n) = FORCE(1:n) + Gradient * IP % s(t) * DetJ * Normal(dim) * Basis(1:n)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------
END SUBROUTINE FlowdepthSolver
!------------------------------------------------------------------------------



