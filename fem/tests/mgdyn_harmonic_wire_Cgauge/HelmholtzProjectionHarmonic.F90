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
! *  This file contains harmonic version of the transormation and also applies the
!    correction to the V field within conducting regions.
! *
! *
! *  Authors: Mika Malinen, Juha Ruokolainen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: March 20, 2020
! *  Last Modified: June 18, 2021, Juha
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
  CALL ListAddNewLogical( SolverParams, 'Linear System Refactorize', .FALSE.)
  CALL ListAddNewLogical(SolverParams,'Variable Output',.FALSE.) 
  
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
  LOGICAL :: ConstantBulkMatrix
!  LOGICAL :: SecondFamily

  INTEGER :: i, j,k,l,n, n_pot, nd_pot, t
  INTEGER :: dim, PotDOFs
  INTEGER :: istat, active

  REAL(KIND=dp), ALLOCATABLE, TARGET :: Stiff(:,:), Force(:,:), PotSol(:,:), F(:,:)
  REAL(KIND=dp) :: Norm, Omega
  REAL(KIND=dp), POINTER :: SaveRHS(:), SOL(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: PotName

  TYPE(Variable_t), POINTER :: v

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

    ALLOCATE( FORCE(2,n), STIFF(n,n), PotSol(2,n), STAT=istat )
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

  n = Solver % Matrix % NumberOfRows
  ALLOCATE(F(n,2)); F=0

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize(Solver)

  SaveRHS => Solver % Matrix % RHS

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
        SecondOrder, n_pot, nd_pot, PotSol )
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    Solver % Matrix % RHS => F(:,1)
    CALL DefaultUpdateForce(FORCE(1,:))

    Solver % Matrix % RHS => F(:,2)
    CALL DefaultUpdateEquations(STIFF, FORCE(2,:))
  END DO

  CALL DefaultFinishBulkAssembly()
  CALL DefaultDirichletBCs()

  v => VariableGet( Mesh % Variables, 'P'  )
  SOL => v % Values

  Solver % Matrix % RHS => F(:,1)
  CALL DefaultDirichletBCs()
  Norm = DefaultSolve()  
  SOL(1::2) = Solver % Variable % Values

  Solver % Matrix % RHS => F(:,2)
  CALL DefaultDirichletBCs()
  Norm = DefaultSolve()  
  SOL(2::2) = Solver % Variable % Values

  Solver % Matrix % RHS => SaveRHS


  omega = GetAngularFrequency()

  !
  ! Finally, redefine the potential variable:
  ! -----------------------------------------
  DO i=1,Solver % Mesh % NumberOfNodes
    j = Solver % Variable % Perm(i)
    IF(j==0) CYCLE

    k = SolverPtr % Variable % Perm(i)
    IF (k == 0) THEN
      CALL Fatal('RemoveKernelComponent', &
        'The variable and potential permutations are nonmatching?')
    END IF

    SolverPtr % Variable % Values(2*k-1) = SolverPtr % Variable % Values(2*k-1) - &
        omega * SOL(2*j)

    SolverPtr % Variable % Values(2*k) = SolverPtr % Variable % Values(2*k) + &
        omega * SOL(2*j-1)
  END DO

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Stiff, Force, Element, n, dim, PiolaVersion, &
               SecondOrder, n_pot, nd_pot, PotSol )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Stiff(:,:), Force(:,:)
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

    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), A(2,3)
    REAL(KIND=dp) :: u, v, w, s, DetJ
    REAL(KIND=dp) :: WBasis(nd_pot-n_pot,3), CurlWBasis(nd_pot-n_pot,3)
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes)

    STIFF = 0.0_dp
    FORCE = 0.0_dp

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
            Basis=Basis, EdgeBasis=WBasis, dBasisdx=dBasisdx, &
            BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
      ELSE
        stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis, dBasisdx)
        IF( dim == 3 ) THEN
          CALL GetEdgeBasis(Element, WBasis, CurlWBasis, Basis, dBasisdx)
        ELSE
          CALL Fatal('HelmholtzProjector', 'Use Piola Transform = True needed in 2D')
        END IF
      END IF
      s = detJ * IP % s(t)

      A = MATMUL(PotSol(:,n_pot+1:nd_pot), WBasis(1:nd_pot-n_pot,:))

      DO p=1,n
        DO q=1,n
          STIFF(p,q) = STIFF(p,q) + SUM(dBasisdx(q,1:dim) * dBasisdx(p,1:dim)) * s
        END DO
      END DO

      DO p=1,n
        FORCE(1,p) = FORCE(1,p) + SUM(A(1,:) * dBasisdx(p,:)) * s
        FORCE(2,p) = FORCE(2,p) + SUM(A(2,:) * dBasisdx(p,:)) * s
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
  CALL ListAddLogical(SolverParams, 'Linear System Refactorize', .FALSE.)

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

  CALL ListAddNewLogical( SolverParams,"Hcurl Basis",.TRUE.)
  CALL ListAddNewLogical(SolverParams,'Variable Output',.FALSE.) 

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
! LOGICAL :: SecondFamily
  LOGICAL :: PiolaVersion, SecondOrder
  LOGICAL :: ConstantBulkMatrix

  INTEGER :: dim, PotDOFs
  INTEGER :: i, j, k, n, nd, n_pot, nd_pot, t
  INTEGER :: istat, active

  REAL(KIND=dp), ALLOCATABLE, TARGET :: Stiff(:,:), Force(:,:), PhiSol(:,:), &
                    SOL(:), F(:,:)
  REAL(KIND=dp) :: Norm
  CHARACTER(LEN=MAX_NAME_LEN) :: PotName, Name

  REAL(KIND=dp), POINTER :: SaveRHS(:)


  TYPE(Variable_t), POINTER :: v

  SAVE Stiff, Force, PhiSol, AllocationsDone
!------------------------------------------------------------------------------
  CALL DefaultStart()

  dim = CoordinateSystemDimension()
  SolverParams => GetSolverParams()
  Mesh => GetMesh()

  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  IF (.NOT. AllocationsDone) THEN
    n = Mesh % MaxElementDOFs
    ALLOCATE(FORCE(2,n),STIFF(n,n),PhiSol(2,n),STAT=istat)
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
  IF (PotDOFs < 2) CALL Fatal('RemoveKernelComponent', 'A complex-valued potential expected')

  !
  ! Find the variable which defines the kernel component:
  !
  Name = GetString(SolverParams, 'Kernel Variable', Found)
  IF (.NOT. Found ) Name = 'phi'
  V => VariableGet( Mesh % Variables, Name )

  Found = ASSOCIATED(v)
   
  IF (.NOT. Found ) THEN
    CALL Fatal('RemoveKernelComponent', 'Solver associated with kernel variable > '&
        //TRIM(Name)//' < not found!')
  END IF
  
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

  n = Solver % Matrix % NumberOfRows
  ALLOCATE(F(n,2)); F=0
  SaveRHS => Solver % Matrix % RHS

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize(Solver)

  DO t=1,active
    Element => GetActiveElement(t)

    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
   
    ! The DOF counts for the potential variable: 
    n_pot = n*SolverPtr % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)
    nd_pot = GetElementNOFDOFs(USolver=SolverPtr)

    IF (nd /= nd_pot-n_pot) CALL Fatal('RemoveKernelComponent', &
     'Potential variable DOFs count /= the solver DOFs count')

    CALL GetLocalSolution(PhiSol, Name)

    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix( STIFF, FORCE, Element, n, nd, dim, PiolaVersion, &
                SecondOrder, PhiSol )
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    Solver % Matrix % RHS => F(:,1)
    CALL DefaultUpdateForce(FORCE(1,:))

    Solver % Matrix % RHS => F(:,2)
    CALL DefaultUpdateEquations(STIFF, FORCE(2,:))
  END DO

  n = Solver % Matrix % NumberOfRows
  ALLOCATE(SOL(2*n))

  Solver % Matrix % RHS => F(:,1)
  CALL DefaultDirichletBCs()
  Norm = DefaultSolve()  
  SOL(1::2) = Solver % Variable % Values

  Solver % Matrix % RHS => F(:,2)
  CALL DefaultDirichletBCs()
  Norm = DefaultSolve()  
  SOL(2::2) = Solver % Variable % Values

  Solver % Matrix % RHS => SaveRHS

  !
  ! Finally, redefine the potential variable:
  !
  n = SIZE(Solver % Variable % Perm(:))
  IF (n ==  SIZE(SolverPtr % Variable % Perm(:))) THEN
    DO i=Solver % Mesh % NumberOfNodes+1,n
      j = Solver % Variable % Perm(i)
      IF (j<=0) CYCLE

      k = SolverPtr % Variable % Perm(i)
      IF (k<=0) THEN
        CALL Fatal('RemoveKernelComponent', &
          'The variable and potential permutations are nonmatching?')
      END IF

      SolverPtr % Variable % Values(2*k-1) = SolverPtr % Variable % Values(2*k-1) - &
          SOL(2*j-1)

      SolverPtr % Variable % Values(2*k) = SolverPtr % Variable % Values(2*k) - &
          SOL(2*j)
    END DO
  ELSE
    CALL Fatal('RemoveKernelComponent', 'The variable and potential permutations differ')  
  END IF

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Stiff, Force, Element, n, nd, dim, PiolaVersion, &
              SecondOrder, PhiSol )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:,:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, dim
    REAL(KIND=dp) :: PhiSol(:,:)
    LOGICAL :: PiolaVersion, SecondOrder
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    LOGICAL :: Stat

    INTEGER :: i, j, p, q, t, EdgeBasisDegree 

    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), A(2,3)
    REAL(KIND=dp) :: u, v, w, s, DetJ
    REAL(KIND=dp) :: WBasis(nd,3), CurlWBasis(nd,3)
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes)

    STIFF = 0.0_dp
    FORCE = 0.0_dp

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
            Basis=Basis, EdgeBasis=WBasis, dBasisdx=dBasisdx, &
            BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
      ELSE
        stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis, dBasisdx)
        IF( dim == 3 ) THEN
          CALL GetEdgeBasis(Element, WBasis, CurlWBasis, Basis, dBasisdx)
        ELSE
          CALL Fatal('RemoveKernelComponent', 'Use Piola Transform = True needed in 2D')
        END IF
      END IF

      s = detJ * IP % s(t)
      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + s * SUM(WBasis(q,:) * WBasis(p,:))
        END DO
      END DO

      A = MATMUL( PhiSol(:,1:n), dBasisdx(1:n,:) )
      DO q=1,nd
        FORCE(1,q) = FORCE(1,q) + s * SUM(A(1,:) * WBasis(q,:))
        FORCE(2,q) = FORCE(2,q) + s * SUM(A(2,:) * WBasis(q,:))
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE RemoveKernelComponent
!------------------------------------------------------------------------------
