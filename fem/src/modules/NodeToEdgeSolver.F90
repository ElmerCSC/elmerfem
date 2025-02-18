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
  CALL ListAddNewLogical( Params,'Mesh Output',.FALSE.)
  CALL ListAddNewInteger( Params,'Active Coordinate', 3)
  CALL ListAddNewInteger( Params,'Interpolation Passive Coordinate', 3)
  
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
  CHARACTER(LEN=MAX_NAME_LEN) :: Name, VarName, TargetName
  INTEGER :: i, j, k, l, dofs, n, m, NoVar, SolverInd, layers, maxperm
  INTEGER, POINTER :: pPerm(:)
  REAL(KIND=dp), POINTER :: pVals(:)
  LOGICAL :: Found, CreateVar, LagrangeCopy, DoIt
  TYPE(Solver_t), POINTER :: pSolver
  CHARACTER(*), PARAMETER :: Caller = 'ExtrudedRestart'
    
  CALL Info(Caller,'Mapping result between meshes')
    
  ThisMesh => Getmesh()
  CALL Info(Caller,'This mesh name is: '//TRIM(ThisMesh % Name),Level=20)   
  CALL Info(Caller,'This mesh dimension is: '//I2S(ThisMesh % MeshDim),Level=20)   

  Params => GetSolverParams()

  TargetMesh => NULL()
  
  SolverInd = ListGetInteger( Params,'Target Mesh Solver Index',Found ) 
  IF( Found ) THEN
    ! Target mesh solver is explicitly given
    TargetMesh => CurrentModel % Solvers(SolverInd) % Mesh
    IF( .NOT. ASSOCIATED( TargetMesh ) ) THEN
      CALL Fatal(Caller,'No mesh associated to Solver: '//I2S(SolverInd))
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

  CALL Info(Caller,'Using target mesh from Solver index: '//I2S(SolverInd),Level=8)
  CALL Info(Caller,'Target mesh name is: '//TRIM(TargetMesh % Name),Level=8)   
  CALL Info(Caller,'Target mesh dimension is: '//I2S(TargetMesh % MeshDim),Level=8)   

  pSolver => Model % Solvers(SolverInd)
  
  NoVar = 0
  DO i = 1,100    
    WRITE (Name,'(A,I0)') 'Variable ',i
    IF( .NOT. ListCheckPresent( Params, Name ) ) EXIT
    NoVar = i
  END DO
  CALL Info(Caller,'Number of variables to be mapped: '//I2S(NoVar),Level=8)   
 
  m = TargetMesh % NumberOfNodes
  IF( ASSOCIATED( pSolver % Variable ) ) THEN
    m = MAX( m, SIZE( pSolver % Variable % Perm ) )
  END IF

  layers = TargetMesh % NumberOfNodes / ThisMesh % NumberOfNodes
  CALL Info(Caller,'Number of layers to be mapped: '//I2S(layers),Level=8)
  
  
  DO i = 1,NoVar
    WRITE (Name,'(A,I0)') 'Variable ',i
    VarName = GetString( Params, Name, Found )  

    Var => VariableGet( ThisMesh % Variables, VarName, ThisOnly = .TRUE. )
    IF(.NOT. ASSOCIATED(Var)) THEN
      IF(InfoActive(20)) THEN
        PRINT *,'List of variable in 2D mesh:'
        Var => ThisMesh % Variables
        DO WHILE(ASSOCIATED(Var))
          PRINT *,TRIM(Var % Name), SIZE(Var % Values), ASSOCIATED(Var % Perm)
          Var => Var % Next
        END DO
      END IF
      CALL Fatal(Caller,'Could not find variable: '//TRIM(VarName))
    END IF
    dofs = Var % Dofs
    
    IF( InfoActive( 20 ) ) THEN
      CALL VectorValuesRange(Var % Values,SIZE(Var % Values),TRIM(VarName))
    END IF
    
    WRITE (Name,'(A,I0)') 'Extruded Variable ',i
    TargetName = GetString( Params, Name, Found )  
    IF(.NOT. Found) TargetName = VarName
           
    pVar => VariableGet( TargetMesh % Variables, TargetName, ThisOnly = .TRUE. )
    CreateVar = .NOT. ASSOCIATED(pVar)    

    IF(.NOT. ASSOCIATED( Var % Perm ) ) THEN
      ! One intended use of this module is to extrude data from 2D electrical machine computation
      ! to 3D one. The it is often desirable also to copy the related electrical circuits that may be
      ! found in the Lagrange multiplier values not associated to any permutation. So these are
      ! copied as one-to-one from 2D to 3D mesh. 
      !------------------------------------------------------------------------------------------------
      n = SIZE(Var % Values)    
      IF( CreateVar ) THEN
        NULLIFY(pVals)
        ALLOCATE(pVals(n))
        pVals = 0.0_dp
        CALL VariableAddVector( TargetMesh % Variables,TargetMesh,&
            Model % Solvers(SolverInd),TargetName,Var % Dofs,pVals)        
      ELSE
        pVals => pVar % Values
      END IF
      pVals(1:n) = Var % Values
      CALL Info(Caller,'Copied variable as such from 2D mesh to 3D mesh: '//TRIM(VarName),Level=8)      
    ELSE      
      maxperm = MAXVAL( Var % Perm )       
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
      
      ! Here we assume that the fields to be mapped are nodal ones and the mesh is linear one!
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

      IF( CreateVar ) THEN
        CALL VariableAddVector( TargetMesh % Variables,TargetMesh,Model % Solvers(SolverInd),&
            TargetName,Var % Dofs,pVals,pPerm,VarType=Var % TYPE)  
        pVar => VariableGet( TargetMesh % Variables, VarName, ThisOnly = .TRUE. )      
        IF(ASSOCIATED(pVar) ) THEN
          CALL Info(Caller,'Created variable: '//TRIM(VarName),Level=20)
        ELSE
          CALL Fatal(Caller,'Failed to create variable: '//TRIM(VarName))
        END IF
      END IF

      IF(.NOT. Var % PeriodicFlipActive ) THEN
        IF( COUNT( Var % Perm > 0 ) > SIZE( Var % Values ) / Var % Dofs ) THEN
          CALL Info(Caller,'The variable that we read in seems to have been conforming!')
          Var % PeriodicFlipActive = .TRUE.
        END IF
      END IF
      
      ! Inherit the periodic flips, if any
      pVar % PeriodicFlipActive = Var % PeriodicFlipActive  

      CALL Info(Caller,'Extruded variable from 2D mesh to 3D mesh: '//TRIM(VarName),Level=8)     
    END IF

    IF( InfoActive( 20 ) ) THEN
      CALL VectorValuesRange(pVals,SIZE(pVals),TRIM(TargetName))
    END IF    

  END DO
    
  CALL Info(Caller,'Transferred '//I2S(NoVar)//' variables from 2d to 3d mesh!',Level=7)

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
  LOGICAL :: Found, PiolaVersion, SecondOrder, SecondKind
!------------------------------------------------------------------------------
  Params => GetSolverParams()
  CALL ListAddLogical(Params, 'Linear System Refactorize', .FALSE.)

  IF (.NOT. ListCheckPresent(Params, "Element")) THEN
    ! We use one place where all the edge element keywords are defined and checked.
    CALL EdgeElementStyle(Params, PiolaVersion, SecondKind, SecondOrder, Check = .TRUE. )

    IF (SecondOrder) THEN
      CALL ListAddString(Params, "Element", &
          "n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2")
    ELSE IF(SecondKind) THEN
      CALL ListAddString(Params, "Element", &
          "n:0 e:2")
    ELSE IF (PiolaVersion) THEN
      CALL ListAddString(Params, "Element", &
          "n:0 e:1 -brick b:3 -quad_face b:2")
    ELSE
      CALL ListAddString( Params, "Element", "n:0 e:1")
    END IF
  END IF

  CALL ListAddNewString( Params,'Variable','-nooutput edgefield proj')
  CALL ListAddInteger( Params, 'Time derivative Order', 0)
  CALL ListAddInteger( Params, 'Active Coordinate', 3)
  CALL ListAddLogical( Params, 'Edge Basis',.TRUE.)
  
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
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Element_t), POINTER :: Element

  LOGICAL :: Found
  LOGICAL :: PiolaVersion
  LOGICAL :: ConstantBulkMatrix, ReadySystemMatrix, IsComplex, IsIm
  TYPE(Variable_t), POINTER :: NodalVar, EdgeVar, ThisVar
  
  INTEGER :: dim, dofs, i, j, k, l, n, nd, t, imoffset
  INTEGER :: istat, active, ActiveComp, EdgeBasisDegree
  INTEGER, POINTER :: NodeIndexes(:)  
  REAL(KIND=dp), ALLOCATABLE :: Stiff(:,:), Force(:), Anodal(:,:)
  REAL(KIND=dp) :: Norm
  TYPE(Solver_t), POINTER :: pSolver  
  CHARACTER(LEN=MAX_NAME_LEN) :: Name
  CHARACTER(*), PARAMETER :: Caller = 'NodeToEdgeField'


  !------------------------------------------------------------------------------

  CALL Info(Caller,'Projecting 3D nodal field to Hcurl field')
  
  CALL DefaultStart()

  dim = CoordinateSystemDimension()
  Params => GetSolverParams()
  Mesh => GetMesh()
  pSolver => Solver
  
  ! Check if accelerated assembly is desired:
  ConstantBulkMatrix = GetLogical(Params, 'Constant Bulk Matrix', Found)
  ReadySystemMatrix = ASSOCIATED(Solver % Matrix % BulkValues) .AND. ConstantBulkMatrix
  IF( ReadySystemMatrix ) CALL Info(Caller,'Assuming that the system matrix is already created!',Level=10)
  
  ThisVar => Solver % Variable
  IF (.NOT. ASSOCIATED(ThisVar) ) THEN
    CALL Fatal(Caller, 'No variable associated to solver!')
  END IF
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

  CALL Info(Caller,'Using nodal variable for projection: '//TRIM(Name),Level=10)  
  IF( InfoActive( 20 ) ) THEN
    CALL VectorValuesRange(NodalVar % Values,SIZE(NodalVar % Values),TRIM(NodalVar % Name))
  END IF
  
  dofs = NodalVar % Dofs
  IsComplex = ( dofs == 6 )
  IsIm = .FALSE.
  imoffset = 0
    
  IF( dofs < dim ) THEN
    ActiveComp = ListGetInteger( Params,'Active Coordinate')
    CALL Info(Caller,'Setting active coordinate to be: '//I2S(ActiveComp),Level=10)
  ELSE
    ActiveComp = 0
  END IF

  
  ! Find the variable which is projected:
  !---------------------------------------------------------
  Name = GetString(Params, 'Edge Variable', Found)
  IF (Found ) THEN
    EdgeVar => VariableGet( Mesh % Variables, Name ) 
    IF(.NOT. ASSOCIATED( EdgeVar ) ) THEN
      CALL Fatal(Caller,'Could not find target variable for projection!')
    END IF
  ELSE
    Name = 'av'    
    EdgeVar => VariableGet( Mesh % Variables, Name ) 
    IF(.NOT. ASSOCIATED( EdgeVar ) ) THEN
      CALL Warn(Caller,'Could not find target variable for projection!')
    END IF
  END IF
    
  IF( InfoActive(20) ) THEN
    CALL VectorValuesRange(NodalVar % Values,SIZE(NodalVar % Values),TRIM(NodalVar % Name))       
  END IF
  
  ! These should be consistent with the primary solver!!
  !-------------------------------------------------------------------------------------
  CALL EdgeElementStyle(Params, PiolaVersion, BasisDegree = EdgeBasisDegree ) 
  IF (PiolaVersion) CALL Info(Caller,'Using Piola-transformed finite elements', Level=5)

  !-----------------------
  ! System assembly:
  !----------------------
  n = Mesh % MaxElementDOFs
  ALLOCATE( Force(n), Stiff(n,n), Anodal(3,n), STAT=istat )

1 CALL DefaultInitialize(Solver, ReadySystemMatrix)
    
  active = GetNOFActive()
  DO t=1,active
    Element => GetActiveElement(t)
    
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    NodeIndexes => Element % NodeIndexes
           
    IF( dofs == 1 ) THEN
      Anodal = 0.0_dp
      WHERE(NodalVar % Perm(NodeIndexes(1:n)) > 0 )
        Anodal(ActiveComp,1:n) = NodalVar % Values( NodalVar % Perm(NodeIndexes(1:n)) )
      END WHERE
    ELSE
      DO i=1,dim
        WHERE(NodalVar % Perm(NodeIndexes(1:n)) > 0 )
          Anodal(i,1:n) = NodalVar % Values( dofs*(NodalVar % Perm(NodeIndexes(1:n))-1)+i+imoffset )
        END WHERE
      END DO
    END IF

    IF( NodalVar % PeriodicFlipActive ) THEN
      DO i=1,dofs
        IF(dofs==1) j=ActiveComp
        WHERE( Mesh % PeriodicFlip(NodeIndexes(1:n) ) )
          Anodal(j,1:n) = -Anodal(j,1:n)
        END WHERE
      END DO
    END IF
        
    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(Stiff, Force, Element, n, nd, dim, PiolaVersion, &
        EdgeBasisDegree, Anodal, ReadySystemMatrix)
       
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

  CALL DefaultFinishBoundaryAssembly()
  IF( InfoActive( 20 ) ) THEN
    CALL VectorValuesRange(Solver % Matrix % Values,SIZE(Solver % Matrix % Values),"A21")
  END IF
  CALL DefaultFinishAssembly()

  IF( InfoActive( 20 ) ) THEN
    CALL VectorValuesRange(Solver % Matrix % Values,SIZE(Solver % Matrix % Values),"A22")
  END IF

  
  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()  

  CALL DefaultFinish()

  IF( InfoActive(20) ) THEN
    CALL VectorValuesRange(ThisVar % Values,SIZE(ThisVar % Values),TRIM(ThisVar % Name))       
  END IF

  ! Finally, redefine the potential variable:
  !---------------------------------------------------------
  IF(ASSOCIATED( EdgeVar ) ) THEN
    n = SIZE(Solver % Variable % Perm)
    IF (n /=  SIZE(EdgeVar % Perm)) THEN
      CALL Fatal(Caller, 'The variable and potential permutations differ')  
    END IF
    
    DO i=1,n
      j = Solver % Variable % Perm(i)
      IF (j < 1) CYCLE
      k = EdgeVar % Perm(i)
      IF (k < 1) CALL Fatal(Caller, &
          'The variable and potential permutations are non-matching?')
      IF(IsComplex ) THEN
        k = 2*k
        IF(.NOT. IsIm) k = k-1
      END IF
      EdgeVar % Values(k) = Solver % Variable % Values(j)
    END DO
  END IF

  ! If we are projecting a complex field then redo for the imaginary component
  IF( IsComplex .AND. .NOT. IsIm ) THEN
    CALL Info(Caller,'Now doing the imaginary component')
    IsIm = .TRUE.
    imoffset = dofs / 2
    ReadySystemMatrix = ConstantBulkMatrix
    GOTO 1    
  END IF

  IF( ASSOCIATED( EdgeVar ) ) THEN
    IF( InfoActive(20) ) THEN
      CALL VectorValuesRange(EdgeVar % Values,SIZE(EdgeVar % Values),TRIM(EdgeVar % Name))       
    END IF
  END IF

  ! We should be visiting this routine only once!
! CALL Info(Caller,'Freeing unneeded matrix structures',Level=10)
! CALL FreeMatrix(Solver % Matrix)
  
  CALL Info(Caller,'Finished projection to edge basis!')

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Stiff, Force, Element, n, nd, dim, PiolaVersion, &
      EdgeBasisDegree, Anodal, ReadySystemMatrix)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Stiff(:,:), Force(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, dim
    LOGICAL :: PiolaVersion
    INTEGER :: EdgeBasisDegree
    REAL(KIND=dp) :: Anodal(:,:)
    LOGICAL :: ReadySystemMatrix  ! A flag to suppress the integration of Stiff
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    LOGICAL :: Stat

    INTEGER :: i, j, p, q, t

    REAL(KIND=dp) :: s, DetJ
    REAL(KIND=dp) :: Basis(n), DBasis(n,3) 
    REAL(KIND=dp) :: WBasis(nd,3), CurlWBasis(nd,3)
    REAL(KIND=dp) :: Aip(3)
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes)

    Stiff = 0.0d0
    Force = 0.0d0

    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree=EdgeBasisDegree)
    IF( dim == 2 .AND. .NOT. PiolaVersion ) THEN
      CALL Fatal(Caller, '"Use Piola Transform = True" needed in 2D')
    END IF

    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, DBasis, EdgeBasis = WBasis, &
          RotBasis = CurlWBasis, USolver = pSolver )             
      s = detJ * IP % s(t)

      Aip = 0.0d0
      DO i=1,dim
        Aip(i) = SUM( Anodal(i,1:n) * Basis(1:n) ) 
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
