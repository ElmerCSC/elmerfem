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
! *
! *  Utilities written as solvers to compute the Helmholtz projection P(A)
! *  of a curl-conforming vector field A. The projection can be obtained as 
! *  P(A) = A - W where  W is the curl-conforming field fitted to represent 
! *  grad Phi, with Phi being a H1-regular scalar field.
! *
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: March 20, 2020
! *
!******************************************************************************


!------------------------------------------------------------------------------
SUBROUTINE HelmholtzProjector_Init(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
!------------------------------------------------------------------------------

  SolverParams => GetSolverParams()
  CALL ListAddLogical(SolverParams, 'Linear System Refactorize', .FALSE.)


  
!------------------------------------------------------------------------------
END SUBROUTINE HelmholtzProjector_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute a H1-regular scalar field to obtain the Helmholtz projection P(A)
!> of a curl-conforming vector field A. Given the solution field Phi of this 
!> solver, the projection can be evaluated as P(A) = A - grad Phi.
!------------------------------------------------------------------------------
SUBROUTINE HelmholtzProjector(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: SolverPtr
  TYPE(Element_t), POINTER :: Element

  LOGICAL :: AllocationsDone = .FALSE.
  LOGICAL :: Found
  LOGICAL :: PiolaVersion, SecondOrder
  LOGICAL :: ConstantBulkMatrix, ReadySystemMatrix

  INTEGER :: dim, PotDOFs
  INTEGER :: i, n, n_pot, nd_pot, t
  INTEGER :: istat, active

  REAL(KIND=dp), ALLOCATABLE :: Stiff(:,:), Force(:), PotSol(:,:)
  REAL(KIND=dp) :: Norm

  CHARACTER(LEN=MAX_NAME_LEN) :: PotName

  SAVE Stiff, Force, PotSol, AllocationsDone
!------------------------------------------------------------------------------
  CALL DefaultStart()

  dim = CoordinateSystemDimension()
  SolverParams => GetSolverParams()
  Mesh => GetMesh()

  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  IF (.NOT. AllocationsDone) THEN
    n = Mesh % MaxElementDOFs
    ALLOCATE( Force(n), Stiff(n,n), PotSol(1,n), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'HelmholtzProjection', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
  END IF

  ! Check if accelerated assembly is desired:
  ConstantBulkMatrix = GetLogical(SolverParams, 'Constant Bulk Matrix', Found)
  ReadySystemMatrix = ASSOCIATED(Solver % Matrix % BulkValues) .AND. &
      ConstantBulkMatrix

  !
  ! Find the variable which is projected:
  !
  PotName = GetString(SolverParams, 'Potential Variable', Found)
  IF (.NOT. Found ) PotName = 'av'
  Found = .FALSE.
  DO i=1,Model % NumberOfSolvers
    SolverPtr => Model % Solvers(i)
    IF (PotName == GetVarName(SolverPtr % Variable)) THEN
      Found = .TRUE.
      EXIT
    END IF
  END DO

  IF(Found) THEN
    CALL Info('HelmholtzProjection', 'Solver inherits potential '&
         //TRIM(PotName)//' from solver: '//I2S(i),Level=7)
  ELSE
    CALL Fatal('HelmholtzProjection', 'Solver associated with potential variable > '&
        //TRIM(PotName)//' < not found!')
  END IF
  
  PotDOFs = SolverPtr % Variable % DOFs
  IF (PotDOFs > 1) CALL Fatal('HelmholtzProjection', 'A real-valued potential expected')

  !
  ! Find some parameters to inherit the vector FE basis as defined in 
  ! the primary solver:
  !
  
  CALL EdgeElementStyle(SolverPtr % Values, PiolaVersion, QuadraticApproximation = SecondOrder )
  IF (PiolaVersion) CALL Info('HelmholtzProjection', &
      'Using Piola-transformed finite elements', Level=5)

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize(Solver, ReadySystemMatrix)

  DO t=1,active
    Element => GetActiveElement(t)
    !
    ! This solver relies on getting basis functions by calling a routine
    ! which returns a curl-conforming basis. It is thus assumed that
    ! the background mesh defines the number of Lagrange basis functions.
    !
    n = GetElementNOFNodes()
   
    ! The DOF counts for the potential (target) variable: 
    n_pot = n*SolverPtr % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)
    nd_pot = GetElementNOFDOFs(USolver=SolverPtr)

    CALL GetVectorLocalSolution(PotSol, PotName, USolver=SolverPtr)

    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(Stiff, Force, Element, n, dim, PiolaVersion, &
        SecondOrder, n_pot, nd_pot, PotSol, ReadySystemMatrix)
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    IF (ReadySystemMatrix) THEN
      CALL DefaultUpdateForce(Force)
    ELSE
      CALL DefaultUpdateEquations(Stiff, Force)
    END IF
  END DO

  IF (ConstantBulkMatrix) THEN 
    CALL DefaultFinishBulkAssembly(BulkUpdate = .NOT.ReadySystemMatrix, RHSUpdate = .FALSE.)
  ELSE
    CALL DefaultFinishBulkAssembly()
  END IF

  !
  ! In this case no need to CALL DefaultFinishAssembly()
  !

  CALL DefaultDirichletBCs()
  !
  ! In connection with the prime application of the solver, the RHS
  ! tends to zero, so zero values may be a better initial guess than
  ! the previous solution:
  !
  Solver % Variable % Values = 0.0d0

  Norm = DefaultSolve()  


CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Stiff, Force, Element, n, dim, PiolaVersion, &
      SecondOrder, n_pot, nd_pot, PotSol, ReadySystemMatrix)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Stiff(:,:), Force(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n   ! The number of background element nodes
    INTEGER :: dim
    LOGICAL :: PiolaVersion, SecondOrder
    INTEGER :: n_pot, nd_pot      ! The size parameters of target field
    REAL(KIND=dp) :: PotSol(:,:)  ! The values of target field DOFS
    LOGICAL :: ReadySystemMatrix  ! A flag to suppress the integration of Stiff
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    LOGICAL :: Stat

    INTEGER :: i, j, p, q, t, EdgeBasisDegree 

    REAL(KIND=dp) :: u, v, w, s, DetJ
    REAL(KIND=dp) :: Basis(n), DBasis(n,3) 
    REAL(KIND=dp) :: WBasis(nd_pot-n_pot,3), CurlWBasis(nd_pot-n_pot,3)
    REAL(KIND=dp) :: A(1,3)
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes)

    Stiff = 0.0d0
    Force = 0.0d0

    IF (SecondOrder) THEN
      EdgeBasisDegree = 2  
      IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
          EdgeBasisDegree=EdgeBasisDegree)
    ELSE
      EdgeBasisDegree = 1
      IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion)
    END IF

    IF( dim == 2 .AND. .NOT. PiolaVersion) THEN
      CALL Fatal('HelmholtzProjection', '"Use Piola Transform = True" needed in 2D')
    END IF
    
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, DBasis, EdgeBasis = WBasis, &
          RotBasis = CurlWBasis, USolver = SolverPtr ) 

      A(1,1:dim) = MATMUL(PotSol(1,n_pot+1:nd_pot), WBasis(1:nd_pot-n_pot,1:dim))
      s = detJ * IP % s(t)

      IF (.NOT. ReadySystemMatrix) THEN
        DO p=1,n
          DO q=1,n
            STIFF(p,q) = STIFF(p,q) + SUM(DBasis(q,1:dim) * DBasis(p,1:dim)) * s
          END DO
        END DO
      END IF
      DO p=1,n
        Force(p) = Force(p) + SUM(A(1,1:dim) * DBasis(p,1:dim)) * s
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE HelmholtzProjector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE RemoveKernelComponent_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: Found, PiolaVersion, SecondOrder, SecondKind
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  CALL ListAddLogical(SolverParams, 'Linear System Refactorize', .FALSE.)

  IF (.NOT. ListCheckPresent(SolverParams, "Element")) THEN
    ! We use one place where all the edge element keywords are defined and checked.
    CALL EdgeElementStyle(SolverParams, PiolaVersion, SecondKind, SecondOrder, Check = .TRUE. )

    IF (SecondOrder) THEN
      CALL ListAddString(SolverParams, "Element", &
          "n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2")
    ELSE IF (SecondKind) THEN
      CALL ListAddString(SolverParams, "Element","n:0 e:2")
    ELSE IF (PiolaVersion) THEN
      CALL ListAddString(SolverParams, "Element", &
          "n:0 e:1 -brick b:3 -quad_face b:2")
    ELSE
      CALL ListAddString( SolverParams, "Element", "n:0 e:1")
    END IF
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE RemoveKernelComponent_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  Apply the Helmholtz projection on a curl-conforming vector field A
!>  when the kernel component grad phi of A (with respect to the curl operator)
!>  has been computed by using the subroutine HelmholtzProjector. This solver
!>  generates the representation W of grad phi in terms of the curl-conforming
!>  basis and finally redefines A := A - W, with W = grad phi. 
!------------------------------------------------------------------------------
SUBROUTINE RemoveKernelComponent(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: SolverPtr, KerSolverPtr
  TYPE(Element_t), POINTER :: Element

  LOGICAL :: AllocationsDone = .FALSE.
  LOGICAL :: Found
  LOGICAL :: PiolaVersion, SecondOrder
  LOGICAL :: ConstantBulkMatrix, ReadySystemMatrix

  INTEGER :: dim, PotDOFs
  INTEGER :: i, j, k, n, nd, n_pot, nd_pot, t
  INTEGER :: istat, active

  REAL(KIND=dp), ALLOCATABLE :: Stiff(:,:), Force(:), Sol(:,:), PotSol(:,:), PhiSol(:)
  REAL(KIND=dp) :: Norm

  CHARACTER(LEN=MAX_NAME_LEN) :: PotName, Name

  SAVE Stiff, Force, Sol, PotSol, PhiSol, AllocationsDone
!------------------------------------------------------------------------------
  CALL DefaultStart()

  dim = CoordinateSystemDimension()
  SolverParams => GetSolverParams()
  Mesh => GetMesh()

  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  IF (.NOT. AllocationsDone) THEN
    n = Mesh % MaxElementDOFs
    ALLOCATE( Force(n), Stiff(n,n), Sol(1,n), PotSol(1,n), PhiSol(n), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'RemoveKernelComponent', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
  END IF

  ! Check if accelerated assembly is desired:
  ConstantBulkMatrix = GetLogical(SolverParams, 'Constant Bulk Matrix', Found)
  ReadySystemMatrix = ASSOCIATED(Solver % Matrix % BulkValues) .AND. &
      ConstantBulkMatrix

  !
  ! Find the variable which is projected:
  !
  PotName = GetString(SolverParams, 'Potential Variable', Found)
  IF (.NOT. Found ) PotName = 'av'
  Found = .FALSE.
  DO i=1,Model % NumberOfSolvers
    SolverPtr => Model % Solvers(i)
    IF (PotName == GetVarName(SolverPtr % Variable)) THEN
      Found = .TRUE.
      EXIT
    END IF
  END DO

  IF (.NOT. Found ) THEN
    CALL Fatal('RemoveKernelComponent', 'Solver associated with potential variable > '&
        //TRIM(PotName)//' < not found!')
  END IF

  PotDOFs = SolverPtr % Variable % DOFs
  IF (PotDOFs > 1) CALL Fatal('RemoveKernelComponent', 'A real-valued potential expected')

  !
  ! Find the variable which defines the kernel component:
  !
  Name = GetString(SolverParams, 'Kernel Variable', Found)
  IF (.NOT. Found ) Name = 'phi'
  Found = .FALSE.
  DO i=1,Model % NumberOfSolvers
    KerSolverPtr => Model % Solvers(i)
    IF (Name == GetVarName(KerSolverPtr % Variable)) THEN
      Found = .TRUE.
      EXIT
    END IF
  END DO

  IF (.NOT. Found ) THEN
    CALL Fatal('RemoveKernelComponent', 'Solver associated with kernel variable > '&
        //TRIM(Name)//' < not found!')
  END IF
  
  IF (KerSolverPtr % Variable % DOFs > 1) CALL Fatal('RemoveKernelComponent', &
      'A real-valued potential expected')

  !
  ! Find some parameters to inherit the vector FE basis as defined in the primary solver:
  !
  
  CALL EdgeElementStyle(SolverPtr % Values, PiolaVersion, QuadraticApproximation = SecondOrder )

  IF (PiolaVersion) CALL Info('RemoveKernelComponent', &
      'Using Piola-transformed finite elements', Level=5)

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize(Solver, ReadySystemMatrix)

  DO t=1,active
    Element => GetActiveElement(t)

    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
   
    ! The DOF counts for the potential variable: 
    n_pot = n*SolverPtr % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)
    nd_pot = GetElementNOFDOFs(USolver=SolverPtr)

    IF (nd /= nd_pot-n_pot) CALL Fatal('RemoveKernelComponent', &
        'Potential variable DOFs count /= the solver DOFs count')

    CALL GetVectorLocalSolution(PotSol, PotName, USolver=SolverPtr)
    CALL GetScalarLocalSolution(PhiSol, Name, USolver=KerSolverPtr)

    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(Stiff, Force, Element, n, nd, dim, PiolaVersion, &
        SecondOrder, PhiSol, ReadySystemMatrix)
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    IF (ReadySystemMatrix) THEN
      CALL DefaultUpdateForce(Force)
    ELSE
      CALL DefaultUpdateEquations(Stiff, Force)
    END IF

  END DO

  IF (ConstantBulkMatrix) THEN 
    CALL DefaultFinishBulkAssembly(BulkUpdate = .NOT.ReadySystemMatrix, RHSUpdate = .FALSE.)
  ELSE
    CALL DefaultFinishBulkAssembly()
  END IF

  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  !
  ! In connection with the prime application of the solver, the RHS
  ! tends to zero, so zero values may be a better initial guess than
  ! the previous solution:
  !
  Solver % Variable % Values = 0.0d0

  Norm = DefaultSolve()  

  !
  ! Finally, redefine the potential variable:
  !
  n = SIZE(Solver % Variable % Perm(:))
  IF (n ==  SIZE(SolverPtr % Variable % Perm(:))) THEN
    DO i=1,n
      j = Solver % Variable % Perm(i)
      IF (j < 1) CYCLE
      k = SolverPtr % Variable % Perm(i)
      IF (k < 1) CALL Fatal('RemoveKernelComponent', &
          'The variable and potential permutations are nonmatching?')
      SolverPtr % Variable % Values(k) = SolverPtr % Variable % Values(k) - &
          Solver % Variable % Values(j)
    END DO
  ELSE
    CALL Fatal('RemoveKernelComponent', 'The variable and potential permutations differ')  
  END IF

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Stiff, Force, Element, n, nd, dim, PiolaVersion, &
      SecondOrder, PhiSol, ReadySystemMatrix)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Stiff(:,:), Force(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, dim
    LOGICAL :: PiolaVersion, SecondOrder
    REAL(KIND=dp) :: PhiSol(:)
    LOGICAL :: ReadySystemMatrix  ! A flag to suppress the integration of Stiff
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    LOGICAL :: Stat

    INTEGER :: i, j, p, q, t, EdgeBasisDegree 

    REAL(KIND=dp) :: s, DetJ
    REAL(KIND=dp) :: Basis(n), DBasis(n,3) 
    REAL(KIND=dp) :: WBasis(nd,3), CurlWBasis(nd,3)
    REAL(KIND=dp) :: A(3)
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes)

    Stiff = 0.0d0
    Force = 0.0d0

    IF (SecondOrder) THEN
      EdgeBasisDegree = 2  
    ELSE
      EdgeBasisDegree = 1
    END IF
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree=EdgeBasisDegree)

    IF( dim == 2 .AND. .NOT. PiolaVersion) THEN
      CALL Fatal('RemoveKernelComponent', '"Use Piola Transform = True" needed in 2D')
    END IF

    
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, DBasis, EdgeBasis = WBasis, &
          RotBasis = CurlWBasis, USolver = SolverPtr ) 
      
      s = detJ * IP % s(t)

      A = 0.0d0
      DO i=1,n
        A(1:dim) = A(1:dim) + PhiSol(i) * DBasis(i,1:dim)
      END DO

      IF (.NOT. ReadySystemMatrix) THEN 
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + SUM(WBasis(q,1:dim) * WBasis(p,1:dim)) * s
          END DO
        END DO
      END IF
      DO p=1,nd
        Force(p) = Force(p) + SUM(A(1:dim) * WBasis(p,1:dim)) * s
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE RemoveKernelComponent
!------------------------------------------------------------------------------
