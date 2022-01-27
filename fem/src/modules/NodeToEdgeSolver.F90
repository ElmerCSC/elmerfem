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
! *  Module for mapping results from nodal 2D solution to 3D Hcurl solution.
! *
! *  Author:  Peter Råback, Mika Malinen
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: Jan 2022
! *
! *****************************************************************************/
 
!> \ingroup Solvers
!> \{
!------------------------------------------------------------------------------
SUBROUTINE ExtrudedRestart_init0( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found
    
  Params => GetSolverParams()

  CALL ListAddNewLogical( Params,'Mesh Enforce Local Copy',.TRUE.)
  
END SUBROUTINE ExtrudedRestart_Init0


!------------------------------------------------------------------------------
!> Interpolates fields from one mesh toanother. 
!------------------------------------------------------------------------------
SUBROUTINE ExtrudedRestart( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh, ThisMesh, TargetMesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var, pVar
  CHARACTER(LEN=MAX_NAME_LEN) :: Name, VarName
  INTEGER :: i, j, k, l, dofs, n, m, NoVar, SolverInd, layers, maxperm
  INTEGER, POINTER :: pPerm(:)
  REAL(KIND=dp), POINTER :: pVals(:)
  LOGICAL :: Found, CreateVar
  TYPE(Solver_t), POINTER :: pSolver
  CHARACTER(*), PARAMETER :: Caller = 'ExtrudedRestart'
  
  CALL Info(Caller,'Mapping result between meshes')
    
  ThisMesh => Getmesh()
  CALL Info(Caller,'This mesh name is: '//TRIM(ThisMesh % Name),Level=20)   

  Params => GetSolverParams()

  TargetMesh => NULL()
  
  SolverInd = ListGetInteger( Params,'Target Mesh Solver Index',Found ) 
  IF( Found ) THEN
    ! Target mesh solver is explicitly given
    TargetMesh => CurrentModel % Solvers(SolverInd) % Mesh
    IF( .NOT. ASSOCIATED( TargetMesh ) ) THEN
      CALL Fatal(Caller,'No mesh associated to Solver: '//TRIM(I2S(SolverInd)))
    END IF
  ELSE  
    ! Otherwise use the 1st mesh that is not this old data mesh
    DO i=1,CurrentModel % NumberOfSolvers
      TargetMesh => CurrentModel % Solvers(i) % Mesh 
      IF(.NOT. ASSOCIATED(ThisMesh,TargetMesh) ) THEN
        SolverInd = i
        EXIT
      END IF
    END DO
    IF( SolverInd == 0) THEN
      CALL Fatal(Caller,'Could not find target mesh that is not this mesh!')
    END IF
  END IF

  CALL Info(Caller,'Using target mesh from Solver index: '//TRIM(I2S(SolverInd)),Level=8)
  CALL Info(Caller,'Target mesh name is: '//TRIM(TargetMesh % Name),Level=8)   

  pSolver => Model % Solvers(SolverInd)
  
  NoVar = 0
  DO i = 1,100    
    WRITE (Name,'(A,I0)') 'Variable ',i
    IF( .NOT. ListCheckPresent( Params, Name ) ) EXIT
    NoVar = i
  END DO
  CALL Info(Caller,'Number of variables to be mapped: '//TRIM(I2S(NoVar)),Level=8)   
 
  m = TargetMesh % NumberOfNodes
  IF( ASSOCIATED( pSolver % Variable ) ) THEN
    m = MAX( m, SIZE( pSolver % Variable % Perm ) )
  END IF

  layers = TargetMesh % NumberOfNodes / ThisMesh % NumberOfNodes
  CALL Info(Caller,'Number of layers to be mapped: '//TRIM(I2S(layers)),Level=8)
  
  
  DO i = 1,NoVar
    WRITE (Name,'(A,I0)') 'Variable ',i
    VarName = GetString( Params, Name, Found )  

    Var => VariableGet( ThisMesh % Variables, VarName, ThisOnly = .TRUE. )
    IF(.NOT. ASSOCIATED(Var)) THEN
      CALL Fatal(Caller,'Could not find variable: '//TRIM(VarName))
    END IF
    dofs = Var % Dofs
    maxperm = MAXVAL( Var % Perm ) 

    IF( InfoActive( 20 ) ) THEN
      CALL VectorValuesRange(Var % Values,SIZE(Var % Values),TRIM(VarName))
    END IF
 
    pVar => VariableGet( TargetMesh % Variables, VarName, ThisOnly = .TRUE. )
    CreateVar = .NOT. ASSOCIATED(pVar)
    
    NULLIFY(pPerm,pVals)
    IF(CreateVar) THEN
      ALLOCATE(pPerm(m))
      pPerm = 0
      ALLOCATE(pVals(layers*SIZE(Var % Values)))      
    ELSE
      pPerm => pVar % Perm
      pVals => pVar % Values
    END IF
    pVals = 0.0_dp

    ! Here we assume that the fields to be mapped are nodal ones and the mesh is linear on!!
    n = ThisMesh % NumberOfNodes
    DO j=0,layers-1
      DO k=1,ThisMesh % NumberOfNodes
        IF( Var % Perm(k) == 0 ) CYCLE
        
        IF( CreateVar ) THEN
          pPerm(k+j*n) = Var % Perm(k) + j * maxperm        
        END IF
        DO l=1,dofs
          pVals(dofs*(pPerm(k+j*n)-1)+l) = Var % Values(dofs*(Var % perm(k)-1)+l)
        END DO        
      END DO
    END DO    
    
    IF( InfoActive( 20 ) ) THEN
      CALL VectorValuesRange(pVals,SIZE(pVals),TRIM(VarName)//' extruded')       
    END IF
    
    IF( CreateVar ) THEN
      CALL VariableAddVector( TargetMesh % Variables,TargetMesh,Model % Solvers(SolverInd),&
          VarName,Var % Dofs,pVals,pPerm,VarType=Var % Type)  
      pVar => VariableGet( TargetMesh % Variables, VarName, ThisOnly = .TRUE. )      
      IF(.NOT. ASSOCIATED(pVar) ) THEN
        CALL Fatal(Caller,'Failed to create variable: '//TRIM(VarName))
      END IF
    END IF
  END DO
    
  CALL Info(Caller,'Interpolated variables between meshes',Level=7)

  
END SUBROUTINE ExtrudedRestart




!------------------------------------------------------------------------------
SUBROUTINE NodeToEdgeField_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, PiolaVersion, SecondOrder
!------------------------------------------------------------------------------
  Params => GetSolverParams()
  CALL ListAddLogical(Params, 'Linear System Refactorize', .FALSE.)

  IF (.NOT. ListCheckPresent(Params, "Element")) THEN
    !
    ! Automatization is not perfect due to the early phase when this 
    ! routine is called; 'Use Piola Transform' and 'Quadratic Approximation'
    ! must be repeated in two solver sections.
    !
    PiolaVersion = GetLogical(Params, 'Use Piola Transform', Found)   
    SecondOrder = GetLogical(Params, 'Quadratic Approximation', Found)
    IF (.NOT. PiolaVersion .AND. SecondOrder) THEN
      CALL Warn("NodeToEdgeField_Init0", &
           "Quadratic Approximation requested without Use Piola Transform " &
           //"Setting Use Piola Transform = True.")
      PiolaVersion = .TRUE.
      CALL ListAddLogical(Params, 'Use Piola Transform', .TRUE.)
    END IF

    IF (SecondOrder) THEN
      CALL ListAddString(Params, "Element", &
          "n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2")
    ELSE
      IF (PiolaVersion) THEN
        CALL ListAddString(Params, "Element", &
            "n:0 e:1 -brick b:3 -quad_face b:2")
      ELSE
        CALL ListAddString( Params, "Element", "n:0 e:1")
      END IF
    END IF
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE NodeToEdgeField_Init0
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Project a nodal solution defining a vector field to a Hcurl vector field.
!> If only one component is given then use it as one component of the vector
!> field others being zero.
!------------------------------------------------------------------------------
SUBROUTINE NodeToEdgeField(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Element_t), POINTER :: Element

  LOGICAL :: Found
  LOGICAL :: PiolaVersion, SecondOrder
  LOGICAL :: ConstantBulkMatrix, ReadySystemMatrix
  TYPE(Variable_t), POINTER :: NodalVar, EdgeVar, ThisVar
  
  INTEGER :: dim, dofs, i, j, k, l, n, nd, t
  INTEGER :: istat, active, ActiveComp

  REAL(KIND=dp), ALLOCATABLE :: Stiff(:,:), Force(:), Anodal(:,:)
  REAL(KIND=dp) :: Norm

  CHARACTER(LEN=MAX_NAME_LEN) :: Name
  CHARACTER(*), PARAMETER :: Caller = 'NodeToEdgeField'


!------------------------------------------------------------------------------
  CALL DefaultStart()

  dim = CoordinateSystemDimension()
  Params => GetSolverParams()
  Mesh => GetMesh()

  ! Check if accelerated assembly is desired:
  ConstantBulkMatrix = GetLogical(Params, 'Constant Bulk Matrix', Found)
  ReadySystemMatrix = ASSOCIATED(Solver % Matrix % BulkValues) .AND. &
      ConstantBulkMatrix
    
  ThisVar => Solver % Variable
  IF (.NOT. ASSOCIATED(ThisVar) ) THEN
    CALL Fatal(Caller, 'No variable associated to solver!')
  END IF
  PRINT *,'this var name:',TRIM(ThisVar % Name)
  IF (ThisVar % Dofs /= 1) CALL Fatal(Caller, 'A real-valued potential expected')
  
  ! Find the variable which defines the nodal field to be mapped.
  ! If the nodal field is a scalar field we need to know
  ! to which component it relates to!
  !------------------------------------------------------
  Name = GetString(Params, 'Nodal Variable', Found)
  IF (.NOT. Found ) Name = 'az'
  NodalVar => VariableGet( Mesh % Variables, Name )
  IF(.NOT. ASSOCIATED(NodalVar) ) THEN
    CALL Fatal(Caller,'Nodal variable not associated: '//TRIM(Name))
  END IF

  ! This
  dofs = NodalVar % Dofs
  IF( dofs == 1 ) THEN
    ActiveComp = ListGetInteger( Params,'Active Coordinate',Found )
    IF(.NOT. Found) THEN      
      ActiveComp = dim
      CALL Info(Caller,'Assuming active coordinate to be: '//TRIM(I2S(ActiveComp)),Level=10)
    END IF
  ELSE
    ActiveComp = 0
  END IF
  IF( InfoActive(20) ) THEN
    CALL VectorValuesRange(NodalVar % Values,SIZE(NodalVar % Values),TRIM(NodalVar % Name))       
  END IF

  
  ! These should be consistent with the primary solver!!
  !-------------------------------------------------------------------------------------
  SecondOrder = GetLogical(Params, 'Quadratic Approximation', Found)  
  IF (SecondOrder) THEN
    PiolaVersion = .TRUE.
  ELSE
    PiolaVersion = GetLogical(Params, 'Use Piola Transform', Found) 
  END IF
  IF (PiolaVersion) CALL Info(Caller,'Using Piola-transformed finite elements', Level=5)

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize(Solver, ReadySystemMatrix)

  n = Mesh % MaxElementDOFs
  ALLOCATE( Force(n), Stiff(n,n), Anodal(n,3), STAT=istat )
  
  DO t=1,active
    Element => GetActiveElement(t)

    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()

    IF( dofs == 1 ) THEN
      Anodal = 0.0_dp
      Anodal(1:n,ActiveComp) = NodalVar % Values(NodalVar % Perm(Element % NodeIndexes))
    ELSE
      DO i=1,dofs
        Anodal(1:n,i) = NodalVar % Values(dofs*(NodalVar % Perm(Element % NodeIndexes)-1)+i)
      END DO
    END IF

    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(Stiff, Force, Element, n, nd, dim, PiolaVersion, &
        SecondOrder, Anodal, ReadySystemMatrix)
    
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

  Norm = DefaultSolve()  

  IF( InfoActive(20) ) THEN
    CALL VectorValuesRange(ThisVar % Values,SIZE(ThisVar % Values),TRIM(ThisVar % Name))       
  END IF

  ! Finally, redefine the potential variable:
  ! Find the variable which is projected:
  !---------------------------------------------------------
  Name = GetString(Params, 'Edge Variable', Found)
  IF (.NOT. Found ) Name = 'av'  
  
  EdgeVar => VariableGet( Mesh % Variables, Name ) 
  IF(ASSOCIATED( EdgeVar ) ) THEN
    n = SIZE(Solver % Variable % Perm(:))
    IF (n /=  SIZE(EdgeVar % Perm(:))) THEN
      CALL Fatal(Caller, 'The variable and potential permutations differ')  
    END IF
    
    DO i=1,n
      j = Solver % Variable % Perm(i)
      IF (j < 1) CYCLE
      k = EdgeVar % Perm(i)
      IF (k < 1) CALL Fatal(Caller, &
          'The variable and potential permutations are nonmatching?')
      EdgeVar % Values(k) = Solver % Variable % Values(j)
    END DO
  ELSE
    CALL Warn(Caller,'Could not find target variable for projection!')
  END IF

  IF( InfoActive(20) ) THEN
    CALL VectorValuesRange(EdgeVar % Values,SIZE(EdgeVar % Values),TRIM(EdgeVar % Name))       
  END IF
  
  CALL Info(Caller,'Finished projection to edge basis!')

  
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Stiff, Force, Element, n, nd, dim, PiolaVersion, &
      SecondOrder, Anodal, ReadySystemMatrix)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Stiff(:,:), Force(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, dim
    LOGICAL :: PiolaVersion, SecondOrder
    REAL(KIND=dp) :: Anodal(:,:)
    LOGICAL :: ReadySystemMatrix  ! A flag to suppress the integration of Stiff
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    LOGICAL :: Stat

    INTEGER :: i, j, p, q, t, EdgeBasisDegree 

    REAL(KIND=dp) :: u, v, w, s, DetJ
    REAL(KIND=dp) :: Basis(n), DBasis(n,3) 
    REAL(KIND=dp) :: WBasis(nd,3), CurlWBasis(nd,3)
    REAL(KIND=dp) :: Aip(3)
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
          CALL Fatal(Caller, 'Use Piola Transform = True needed in 2D')
        END IF
      END IF

      s = detJ * IP % s(t)

      Aip = 0.0d0
      DO i=1,dim
        Aip(i) = SUM( Anodal(1:n,i) * Basis(1:n) ) 
      END DO
      
      IF (.NOT. ReadySystemMatrix) THEN 
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + SUM(WBasis(q,1:dim) * WBasis(p,1:dim)) * s
          END DO
        END DO
      END IF
      DO p=1,nd
        Force(p) = Force(p) + SUM(Aip(1:dim) * WBasis(p,1:dim)) * s
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE NodeToEdgeField
!------------------------------------------------------------------------------



!> \}
