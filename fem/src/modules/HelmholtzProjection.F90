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
!  LOGICAL :: SecondFamily

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
      CALL Fatal( 'HelmholtzProjector', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
  END IF

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
    CALL Fatal('HelmholtzProjector', 'Solver associated with potential variable > '&
        //TRIM(PotName)//' < not found!')
  END IF

  PotDOFs = SolverPtr % Variable % DOFs
  IF (PotDOFs > 1) CALL Fatal('HelmholtzProjector', 'A real-valued potential expected')

  !
  ! Find some parameters to inherit the vector FE basis as defined in 
  ! the primary solver:
  !
  SecondOrder = GetLogical(SolverPtr % Values, 'Quadratic Approximation', Found)  

  IF (SecondOrder) THEN
    PiolaVersion = .TRUE.
  ELSE
    PiolaVersion = GetLogical(SolverPtr % Values, 'Use Piola Transform', Found) 
  END IF

  IF (PiolaVersion) CALL Info('HelmholtzProjector', &
      'Using Piola-transformed finite elements', Level=5)

!  SecondFamily = GetLogical(SolverPtr % Values, 'Second Kind Basis', Found)

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize()

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
        SecondOrder, n_pot, nd_pot, PotSol)
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL DefaultUpdateEquations(Stiff, Force)

  END DO

  CALL DefaultFinishBulkAssembly()

  CALL DefaultFinishAssembly()
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
      SecondOrder, n_pot, nd_pot, PotSol)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Stiff(:,:), Force(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n   ! The number of background element nodes
    INTEGER :: dim
    LOGICAL :: PiolaVersion, SecondOrder
    INTEGER :: n_pot, nd_pot      ! The size parameters of target field
    REAL(KIND=dp) :: PotSol(:,:)  ! The values of target field DOFS
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

    DO t=1,IP % n

      u = IP % U(t)
      v = IP % V(t)
      w = IP % W(t)

      IF (PiolaVersion) THEN
        stat = EdgeElementInfo(Element, Nodes, u, v, w, DetF=DetJ, &
            Basis=Basis, EdgeBasis=WBasis, dBasisdx=DBasis, &
            BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
      ELSE
        stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis, DBasis)
        IF( dim == 3 ) THEN
          CALL GetEdgeBasis(Element, WBasis, CurlWBasis, Basis, DBasis)
        ELSE
          CALL Fatal('HelmholtzProjector', 'Use Piola Transform = True needed in 2D')
        END IF
      END IF

      A(1,1:dim) = MATMUL(PotSol(1,n_pot+1:nd_pot), WBasis(1:nd_pot-n_pot,1:dim))
      s = detJ * IP % s(t)

      DO p=1,n
        DO q=1,n
          STIFF(p,q) = STIFF(p,q) + SUM(DBasis(q,1:dim) * DBasis(p,1:dim)) * s
        END DO

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
  LOGICAL :: Found, PiolaVersion, SecondOrder
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  IF (.NOT. ListCheckPresent(SolverParams, "Element")) THEN
    !
    ! Automatization is not perfect due to the early phase when this 
    ! routine is called; 'Use Piola Transform' and 'Quadratic Approximation'
    ! must be repeated in two solver sections.
    !
    PiolaVersion = GetLogical(SolverParams, 'Use Piola Transform', Found)   
    SecondOrder = GetLogical(SolverParams, 'Quadratic Approximation', Found)
    IF (.NOT. PiolaVersion .AND. SecondOrder) THEN
      CALL Warn("RemoveKernelComponent_Init0", &
           "Quadratic Approximation requested without Use Piola Transform " &
           //"Setting Use Piola Transform = True.")
      PiolaVersion = .TRUE.
      CALL ListAddLogical(SolverParams, 'Use Piola Transform', .TRUE.)
    END IF

    IF (SecondOrder) THEN
      CALL ListAddString(SolverParams, "Element", &
          "n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2")
    ELSE
      IF (PiolaVersion) THEN
        CALL ListAddString(SolverParams, "Element", &
            "n:0 e:1 -brick b:3 -quad_face b:2")
      ELSE
        CALL ListAddString( SolverParams, "Element", "n:0 e:1")
      END IF
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
!  LOGICAL :: SecondFamily

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
  SecondOrder = GetLogical(SolverPtr % Values, 'Quadratic Approximation', Found)  

  IF (SecondOrder) THEN
    PiolaVersion = .TRUE.
  ELSE
    PiolaVersion = GetLogical(SolverPtr % Values, 'Use Piola Transform', Found) 
  END IF

  IF (PiolaVersion) CALL Info('RemoveKernelComponent', &
      'Using Piola-transformed finite elements', Level=5)

!  SecondFamily = GetLogical(SolverPtr % Values, 'Second Kind Basis', Found)

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize()

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
        SecondOrder, PhiSol)
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL DefaultUpdateEquations(Stiff, Force)

  END DO

  CALL DefaultFinishBulkAssembly()

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
      SecondOrder, PhiSol)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Stiff(:,:), Force(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, dim
    LOGICAL :: PiolaVersion, SecondOrder
    REAL(KIND=dp) :: PhiSol(:)
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    LOGICAL :: Stat

    INTEGER :: i, j, p, q, t, EdgeBasisDegree 

    REAL(KIND=dp) :: u, v, w, s, DetJ
    REAL(KIND=dp) :: Basis(n), DBasis(n,3) 
    REAL(KIND=dp) :: WBasis(nd,3), CurlWBasis(nd,3)
    REAL(KIND=dp) :: A(3)
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


    DO t=1,IP % n

      u = IP % U(t)
      v = IP % V(t)
      w = IP % W(t)

      IF (PiolaVersion) THEN
        stat = EdgeElementInfo(Element, Nodes, u, v, w, DetF=DetJ, &
            Basis=Basis, EdgeBasis=WBasis, dBasisdx=DBasis, &
            BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
      ELSE
        stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis, DBasis)
        IF( dim == 3 ) THEN
          CALL GetEdgeBasis(Element, WBasis, CurlWBasis, Basis, DBasis)
        ELSE
          CALL Fatal('RemoveKernelComponent', 'Use Piola Transform = True needed in 2D')
        END IF
      END IF

      s = detJ * IP % s(t)

      A = 0.0d0
      DO i=1,n
        A(1:dim) = A(1:dim) + PhiSol(i) * DBasis(i,1:dim)
      END DO
 
      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + SUM(WBasis(q,1:dim) * WBasis(p,1:dim)) * s
        END DO

        Force(p) = Force(p) + SUM(A(1:dim) * WBasis(p,1:dim)) * s
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE RemoveKernelComponent
!------------------------------------------------------------------------------
