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
! *  Module for solving the nodal preconditioning equation associated to AV solver.
! *  It is assumed that the matrix is asssembled by the AV solver and here it is
! *  ready to be used.
! *
! *  Authors: Peter Råback, Mika Malinen
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 3.9.2025
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{

!------------------------------------------------------------------------------
SUBROUTINE APrecSolver_Init( Model,Solver,dt,Transient ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, ElectroDynamics
  LOGICAL :: Monolithic
  CHARACTER(*), PARAMETER :: Caller = 'APrecSolver_Init'
  
  Params => GetSolverParams()
  CALL ListAddNewLogical( Params,'Apply Mortar BCs',.TRUE.)
  CALL ListAddNewLogical( Params,'Use Global Mass Matrix',.TRUE.)
  CALL ListAddNewString( Params,'Exec Solver','never')
  CALL ListAddNewLogical( Params,'Skip Compute Nonlinear Change',.TRUE.)

  Monolithic = ListGetLogical( Params,'Monolithic Solver',Found )  
  IF( Monolithic ) THEN
    CALL ListAddNewString( Params,'Variable','-dofs 3 Nodal A' )
  ELSE
    CALL ListAddNewString( Params,'Variable','-nooutput nodal A tmp')
    CALL ListAddString( Params,&
        NextFreeKeyword('Exported Variable', Params),'-dofs 3 nodal A')
  END IF

  CALL ListAddString( Params,&
      NextFreeKeyword('Exported Variable', Params),'-dofs 3 nodal A rhs')
  CALL ListAddString( Params,&
      NextFreeKeyword('Exported Variable', Params),'-dofs 3 nodal A cum')

  IF( ListGetLogicalAnySolver(Model,'Prec Matrix Cylindrical') ) THEN
    CALL ListAddString( Params,&
        NextFreeKeyword('Exported Variable', Params),'-dofs 3 nodal A cyl')
  END IF


  
!------------------------------------------------------------------------------
END SUBROUTINE APrecSolver_Init ! }}}
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solve the magnetic vector potential expressed in terms of a single component.
!> The solver may take into account rotating boundary conditions.
!> Also optionally compute moments and inertia. 
!------------------------------------------------------------------------------
SUBROUTINE APrecSolver( Model,Solver,dt,Transient ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
  USE CircuitUtils
  USE ZirkaUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: Found
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: i,j,k,n, nb, nd, t
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Mesh_t),   POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: Avar, Svar, NodeResVar, EdgeSolVar, EdgeResVar, pVar
  TYPE(Matrix_t), POINTER, SAVE :: Proj => NULL()
  LOGICAL :: Monolithic, SecondFamily, SecondOrder, PiolaVersion
  TYPE(ValueList_t), POINTER :: EdgeSolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
  LOGICAL, SAVE :: Visited = .FALSE., PrecMatCyl, PrecMatNt, SkipFaces
  REAL(KIND=dp), POINTER :: allrhs(:) => NULL(), BulkValues(:) => NULL()
  INTEGER :: comps, compi, dofs
  INTEGER :: NoVisited = 0
  LOGICAL, POINTER, SAVE :: NodeSkip(:)
  
  CHARACTER(*), PARAMETER :: Caller = 'APrecSolver'


  SAVE :: allrhs, BulkValues, NoVisited 
  
!------------------------------------------------------------------------------

  CALL Info( Caller,'-------------------------------------------------------', Level=10 )
  CALL Info( Caller,'Solving preconditioning equation for vector potential', Level=6 )
  CALL Info( Caller,'-------------------------------------------------------', Level=10 )

  Mesh => Solver % Mesh 
  SolverParams => Solver % Values
  SVar => Solver % Variable
  dofs = SVar % dofs
  
  Monolithic = ListGetLogical( SolverParams,'Monolithic Solver', Found )
  NoVisited = NoVisited + 1

  NodeResVar => VariableGet( Mesh % Variables,'nodal a rhs')
      
  IF(Monolithic) THEN
    ! Solve all 3 components at the same time!
    ! Note that the current assembly in AVSolver is not compatible with this!
    comps = 1
    allrhs => Solver % Matrix % rhs
  ELSE   
    ! Solve one component at a time -> faster!
    comps = 3
    allrhs => NodeResVar % Values
    IF(.NOT. ASSOCIATED(Solver % Matrix % BulkValues)) THEN
      ALLOCATE(Solver % Matrix % BulkValues(SIZE(Solver % Matrix % Values)))
    END IF
    Solver % Matrix % BulkValues = Solver % Matrix % Values
    BulkValues => Solver % Matrix % BulkValues
  END IF

    
  IF(Monolithic) THEN
    AVar => SVar
    PrecMatCyl = ListGetLogicalAnySolver(Model,'Prec Matrix Cylindrical')
  ELSE
    Avar => VariableGet( Mesh % Variables,'Nodal A')        
    IF(.NOT. ASSOCIATED(Avar)) THEN
      CALL Fatal(Caller,'Could not find variable "Nodal A"')
    END IF
    IF(SVar % dofs /= 1) THEN
      CALL Fatal(Caller,'Componentwise solver size should be 1!')
    END IF
    PrecMatCyl = .FALSE.
    PrecMatNt = .FALSE.
  END IF
  IF(AVar % dofs /= 3) THEN
    CALL Fatal(Caller,'Full solution size should be 3!')
  END IF

  EdgeSolVar => NULL()
  sname = ListGetString( SolverParams, 'Edge Update Name', Found)
  IF (Found) THEN
    EdgeSolVar => VariableGet(Mesh % Variables, sname)
    IF (.NOT. ASSOCIATED(EdgeSolVar)) CALL Fatal(Caller, 'Could not found field: '//TRIM(sname))
  ELSE
    CALL Fatal(Caller, 'Give "Edge Update Name" to enable the use as a preconditioner')
  END IF

  EdgeResVar => NULL()  
  sname = ListGetString( SolverParams,'Edge Residual Name',Found)
  IF(.NOT. Found) THEN
    CALL Fatal(Caller, 'Give Edge Residual Name to enable the use as a preconditioner')
  END IF

  EdgeResVar => VariableGet( Mesh % Variables, sname )
  IF(.NOT. ASSOCIATED( EdgeResVar ) ) CALL Fatal(Caller,'Could not find field: '//TRIM(sname))

  EdgeSolverParams => GetSolverParams(EdgeResVar % Solver)
  CALL EdgeElementStyle(EdgeSolverParams, PiolaVersion, SecondFamily, SecondOrder, Check = .TRUE.)

  IF (.NOT. ASSOCIATED(Proj)) THEN
    CALL Info(Caller,'Creating projection matrix to map a nodal solution into vector element space', Level=10)
    SkipFaces = ListGetLogical( SolverParams,'Skip Faces in Projection',Found ) 
    CALL NodalToNedelecInterpolation_GlobalMatrix(Mesh, Avar, EdgeSolVar, Proj, cdim=3, SkipFaces = SkipFaces)
  END IF

  ! Now EdgeResVar represents the residual with respect
  ! to the basis for H(curl). We need to apply a transformation so that
  ! we may solve the residual correction equation by using the nodal basis.
  !-----------------------------------------------------------------------------
  CALL Info(Caller,'Using Transposed Projection Matrix: H(curl) -> H1', Level=10)
  CALL CRS_TransposeMatrixVectorMultiply(Proj, EdgeResVar % Values, allrhs )           
  
  ! Potentially create a mask that avoids residual values being applied on the mortar BC. 
  IF(NoVisited == 1 ) THEN
    n = SIZE(Solver % Matrix % rhs)/dofs
    ALLOCATE(NodeSkip(n))    
    NodeSkip = .FALSE.
    CALL CreateNodeSkipMask(NodeSkip,Solver % Variable)
    n = COUNT(NodeSkip)
    IF(n==0) DEALLOCATE(NodeSkip)
  END IF
  
  IF(Monolithic) THEN
    ! If we use N-T coordinate system to make periodic/rotational BC's easier than we must map the
    ! original residual vector into N-T system. 
    
    IF( Solver % NormalTangential % NormalTangentialNOFNodes > 0 ) THEN
      CALL RotateNtVector( allrhs, Solver )
    END IF

    IF(ListCheckPrefixAnyBodyForce( Model,'Test Load') ) THEN
      CALL SetTestRhs()    
      NodeResVar % Values = allrhs
      PRINT *,'rhs 1:',MINVAL(allrhs(1::3)),MAXVAL(allrhs(1::3))
      PRINT *,'rhs 2:',MINVAL(allrhs(2::3)),MAXVAL(allrhs(2::3))
      PRINT *,'rhs 3:',MINVAL(allrhs(3::3)),MAXVAL(allrhs(3::3))
    END IF
    
    ! By construction do not apply any residual to the mortar boundary. 
    IF(ASSOCIATED(NodeSkip)) THEN
      DO i=1,SIZE(NodeSkip)
        IF(NodeSkip(i)) THEN
          allrhs(dofs*(i-1)+1:dofs*i) = 0.0_dp
        END IF
      END DO
    END IF
      
    IF(ALLOCATED(Solver % Matrix % ConstrainedDOF ) ) &
        Solver % Matrix % ConstrainedDOF = .FALSE.
    CALL DefaultDirichletBCs()
    Norm = DefaultSolve()    
  ELSE
    DO compi = 1, comps      
      Solver % Matrix % Values = BulkValues
      Solver % Matrix % rhs = allrhs(compi::comps)

      ! By construction do not apply any residual to the mortar boundary. 
      IF(ASSOCIATED(NodeSkip)) THEN
        WHERE(NodeSkip)
          Solver % Matrix % rhs = 0.0_dp
        END WHERE
      END IF

      ! Different components will have different Dirichlet BC's!
      ! Note that Solver % Variable now points to correct component of Avar so there is no
      ! need to copy values to it!
      Solver % Variable => VariableGet( Model % Variables,TRIM(Avar % name)//' '//I2S(compi))
      IF(.NOT. ASSOCIATED(Solver % Variable)) THEN
        CALL Fatal(Caller,'Could not find variable for component :'//I2S(compi))
      END IF
      
      IF(ALLOCATED(Solver % Matrix % ConstrainedDOF ) ) &
          Solver % Matrix % ConstrainedDOF = .FALSE.
      CALL DefaultDirichletBCs()
      Norm = DefaultSolve()
    END DO
    Solver % Variable => SVar
  END IF

  IF(PrecMatCyl) THEN
    pVar => VariableGet( Mesh % Variables,'nodal a cyl')
    IF(ASSOCIATED(pVar)) THEN
      pVar % Values = pVar % Values + AVar % Values
    END IF
    CALL CylinderToCartesianProject(SVar % Values)
  END IF
    
  pVar => VariableGet( Mesh % Variables,'nodal a cum')
  IF(ASSOCIATED(pVar)) THEN
    pVar % Values = pVar % Values + AVar % Values
  END IF

  
  CALL Info(Caller,'Projecting nodal solution to vector element space', Level=20)
  CALL CRS_MatrixVectorMultiply(Proj, Avar % Values, EdgeSolVar % Values ) 

  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(Avar % Values,SIZE(Avar % Values),'VecPotNodal')       
    CALL VectorValuesRange(EdgeSolVar % Values,SIZE(EdgeSolVar % Values),'VecPotEdge')       
  END IF



  CALL Info(Caller,'Auxiliary space nodal solution finished!',Level=10)

CONTAINS

  SUBROUTINE CylinderToCartesianProject(b)    
    REAL(KIND=dp) :: b(:)
    REAL(KIND=dp) :: ar, aphi, x, y
    INTEGER :: i,j
    TYPE(Mesh_t), POINTER :: Mesh

    CALL Info(Caller,'Mapping cylindrical vector potential to cartesian!',Level=10)
    
    Mesh => Solver % Mesh 
    
    DO i=1,Mesh % NumberOfNodes
      j = Solver % Variable % Perm(i)
      IF(j==0) CYCLE
      IF(ASSOCIATED(Mesh % PeriodicPerm)) THEN
        IF(Mesh % PeriodicPerm(i) > 0) CYCLE
      END IF
      ar = Solver % Variable % Values(3*j-2)
      aphi = Solver % Variable % Values(3*j-1)
      x = Solver % Mesh % Nodes % x(i)
      y = Solver % Mesh % Nodes % y(i)
      Solver % Variable % Values(3*j-2) = x * ar - y * aphi
      Solver % Variable % Values(3*j-1) = y * ar + x * aphi        
      ! r3d component remains the same!
    END DO
              
  END SUBROUTINE CylinderToCartesianProject
    


!------------------------------------------------------------------------------
  SUBROUTINE RotateNtVector( Vector,Solver )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Vector(:)
    TYPE(Solver_t) :: Solver    
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,m,n,dofs,dim
    REAL(KIND=dp) :: s,Q(3),N1(3),T1(3),T2(3),R(3,3)
    REAL(KIND=dp), POINTER :: Normal(:,:), Tangent1(:,:), Tangent2(:,:)
    !------------------------------------------------------------------------------

    n = Solver % NormalTangential % NormalTangentialNOFNodes
    IF(n==0) RETURN

    dofs = Solver % Variable % dofs
    dim = Mesh % MeshDim

    IF(dim /= dofs) CALL Fatal('RotateNtVector','Currently assuming that dim == dofs !')
    
    Normal => Solver % NormalTangential % BoundaryNormals
    Tangent1 => Solver % NormalTangential % BoundaryTangent1
    Tangent2 => Solver % NormalTangential % BoundaryTangent2
    
    DO i=1,Mesh % NumberOfNodes
      j = Solver % NormalTangential % BoundaryReorder(i)
      k = Solver % Variable % Perm(i)
      IF(j==0 .OR. k==0) CYCLE

      ! Eliminate duplicate projection of conforming nodes.
      IF(ASSOCIATED(Mesh % PeriodicPerm)) THEN
        IF(Mesh % PeriodicPerm(i) > 0) CYCLE
      END IF

      ! Create rotation matrix "R" on this node.
      R = 0.0_dp
      N1 = Normal( j,: )
      
      SELECT CASE(DIM)
      CASE (2)
        R(1,1) =  N1(1)
        R(1,2) =  N1(2)
        R(2,1) = -N1(2)
        R(2,2) =  N1(1)
        R(3,3) = 1.0_dp
      CASE (3)
        T1 = Tangent1( j,: )
        T2 = Tangent2( j,: )

        R(1,1:3) = N1(1:3)
        R(2,1:3) = T1(1:3)
        R(3,1:3) = T2(1:3)
      END SELECT

      ! PRINT *,'NT:',i,j,k,N1,T1,T2
      
      ! Rotate the local vector to N-T coordinates. 
      Q = 0.0_dp
      DO l=1,DOFs
        s = 0.0_dp
        DO m=1,DOFs
          s = s + R(l,m) * Vector(Dofs*(k-1)+m)
        END DO
        Q(l) = s
      END DO
      Vector(DOFs*(k-1)+1:Dofs*k) = Q(1:DOFs)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE RotateNtVector
!------------------------------------------------------------------------------


  SUBROUTINE SetTestRhs()

    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: DetJ,LoadAtIP,Weight
    REAL(KIND=dp), ALLOCATABLE :: FORCE(:), LOAD(:,:), Basis(:), dBasisdx(:,:)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,elem
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Visited = .FALSE.
    SAVE Nodes, Visited, Basis, dBasisdx, FORCE, Load
!------------------------------------------------------------------------------

    allrhs = 0.0_dp

    IF(.NOT. Visited) THEN
      n = Mesh % MaxElementNodes
      ALLOCATE( FORCE(3*n), LOAD(3,n), Basis(n),dBasisdx(n,3) )
      Visited = .TRUE.
    END IF
    
    DO elem=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(elem)


      CALL GetElementNodes( Nodes )
      FORCE = 0._dp
      LOAD = 0._dp

      BodyForce => GetBodyForce(Element)
      IF ( .NOT. ASSOCIATED(BodyForce) ) CYCLE

      Load(1,1:n) = GetReal( BodyForce,'Test Load 1', Found )
      Load(2,1:n) = GetReal( BodyForce,'Test Load 2', Found )
      Load(3,1:n) = GetReal( BodyForce,'Test Load 3', Found )

      IP = GaussPoints( Element )
      
      DO t=1,IP % n
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )
        Weight = IP % s(t) * DetJ
        DO i=1,3      
          LoadAtIP = SUM( Basis(1:n) * LOAD(i,1:n) )
          FORCE(i:3:3*n) = FORCE(i:3:3*n) + Weight * LoadAtIP * Basis(1:n)
        END DO
      END DO

      CALL DefaultUpdateForce(FORCE,Element,Solver)

    END DO
    
  END SUBROUTINE SetTestRhs

  
  
END SUBROUTINE APrecSolver
!------------------------------------------------------------------------------

!> \}

