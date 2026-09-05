!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
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
SUBROUTINE AmsVectorSolver_Init( Model,Solver,dt,Transient ) ! {{{
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
  LOGICAL :: Found, IsMonolithic, IsComplex
  CHARACTER(*), PARAMETER :: Caller = 'AmsVectorSolver_Init'

  Params => GetSolverParams()
  CALL ListAddLogical( Params,'AMS Vector Solver',.TRUE.)
  CALL ListAddNewLogical( Params,'Apply Mortar BCs',.TRUE.)
  CALL ListAddNewLogical( Params,'Use Global Mass Matrix',.TRUE.)
  CALL ListAddNewString( Params,'Exec Solver','never')
  CALL ListAddNewLogical( Params,'Skip Compute Nonlinear Change',.TRUE.)
  CALL ListAddNewInteger( Params,'Nonlinear System Max Iterations', 1)

  !CALL ListAddLogical(Params,'Linear System Refactorize',.TRUE.)
  !CALL ListAddLogical(Params,'Mortar BCs Fixed',.FALSE.)

  
  IsMonolithic = ListGetLogical( Params,'Monolithic Solver',Found )  
  IsComplex = ListGetLogical( Params, 'Linear System Complex', Found )
  
  IF( IsMonolithic ) THEN
    ! We solve the equation as monolithic system
    IF( IsComplex ) THEN
      CALL ListAddNewString( Params,'Variable',&
          'AmsVec[amsx re:1 amsx im:1 amsy re:1 amsy im:1 amsz re:1 amsz im:1]') 
    ELSE
      CALL ListAddNewString( Params,'Variable',&
          'AmsVec[amsx:1 amsy:1 amsz:1]') 
    END IF
  ELSE    
    ! We solve the equation component-wise. Hence the primary variable is a temporary one.
    CALL ListAddNewLogical( Params,'Variable Output',.FALSE.)   
    IF( IsComplex ) THEN
      CALL ListAddNewString( Params,'Variable','Amstmp[amst re:1 amst im:1]')
      CALL ListAddString( Params,&
          NextFreeKeyword('Exported Variable', Params), &
          'AmsVec[amsx re:1 amsx im:1 amsy re:1 amsy im:1 amsz re:1 amsz im:1]') 
    ELSE
      CALL ListAddNewString( Params,'Variable','amstmp')
      CALL ListAddString( Params,&
          NextFreeKeyword('Exported Variable', Params), &
          'AmsVec[amsx:1 amsy:1 amsz:1]') 
    END IF
  END IF

  IF( IsComplex ) THEN
    CALL ListAddString( Params,&
        NextFreeKeyword('Exported Variable', Params),'-dofs 6 nodal amsa rhs')
  ELSE
    CALL ListAddString( Params,&
        NextFreeKeyword('Exported Variable', Params),'-dofs 3 nodal amsa rhs')
  END IF
    
!------------------------------------------------------------------------------
END SUBROUTINE AmsVectorSolver_Init ! }}}
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solve the magnetic vector potential expressed in terms of a single component.
!> The solver may take into account rotating boundary conditions.
!> Also optionally compute moments and inertia. 
!------------------------------------------------------------------------------
SUBROUTINE AmsVectorSolver( Model,Solver,dt,Transient ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulatio
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: Found
  REAL(KIND=dp) :: Norm
  INTEGER :: i,j,k,n, nb, nd, t, ns, n0, n1
  TYPE(Mesh_t),   POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: Avar, SVar, NodeResVar, EdgeSolVar, EdgeResVar, pVar
  TYPE(Matrix_t), POINTER, SAVE :: Proj => NULL()
  LOGICAL :: IsMonolithic, SecondFamily, SecondOrder, PiolaVersion, ExtrudedSol
  TYPE(ValueList_t), POINTER :: EdgeSolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
  LOGICAL, SAVE :: Visited = .FALSE., SkipFaces, PrecMatNt
  REAL(KIND=dp), POINTER :: allrhs(:) => NULL()
  INTEGER :: comps, compi, dofs, dof
  LOGICAL, POINTER, SAVE :: NodeSkip(:)
  TYPE(Matrix_t), POINTER :: A
  REAL(KIND=dp), POINTER :: b(:)
  CHARACTER(*), PARAMETER :: Caller = 'AmsVectorSolver'

  
!------------------------------------------------------------------------------

  CALL Info( Caller,'-------------------------------------------------------', Level=10 )
  CALL Info( Caller,'Solving preconditioning equation for AMS vector', Level=6 )
  CALL Info( Caller,'-------------------------------------------------------', Level=10 )

  Mesh => Solver % Mesh
  SolverParams => Solver % Values
  SVar => Solver % Variable

  A => Solver % Matrix
  b => A % Rhs

  IsMonolithic = ListGetLogical( SolverParams,'Monolithic Solver', Found )
  
  NodeResVar => VariableGet( Mesh % Variables,'nodal amsa rhs',&
      ThisOnly=.TRUE.,UnfoundFatal=.TRUE.)
  allrhs => NodeResVar % Values

  ns = 1
  IF( ListGetLogical( SolverParams,'Linear System Complex',Found ) ) ns = 2
  
  IF(IsMonolithic) THEN
    ! Solve all 3 components at the same time!
    ! Note that the current assembly in AVSolver is not compatible with this!
    comps = 1
  ELSE   
    ! Solve one component at a time -> faster!
    comps = 3
    IF(.NOT. ASSOCIATED(A % BulkValues)) THEN
      ALLOCATE(A % BulkValues(SIZE(A % Values)))
    END IF
    A % BulkValues = A % Values
  END IF
  
  IF(IsMonolithic) THEN
    AVar => SVar
  ELSE
    Avar => VariableGet( Mesh % Variables,'AmsVec',ThisOnly=.TRUE.,UnfoundFatal=.TRUE.)        
    IF(SVar % dofs /= ns) THEN
      CALL Fatal(Caller,'Componentwise solver size should be: '//I2S(ns))
    END IF
    PrecMatNt = .FALSE.
  END IF
  dofs = AVar % dofs
  IF(dofs /= 3*ns) THEN
    CALL Fatal(Caller,'Full solution size should be: '//I2S(2*ns))
  END IF
  AVar % Values = 0.0_dp
  
  sname = ListGetString( SolverParams, 'Edge Update Name', Found)
  IF(.NOT. Found) sname = "ams update"
  EdgeSolVar => VariableGet(Mesh % Variables, sname, ThisOnly=.TRUE.,UnfoundFatal=.TRUE.)
  EdgeSolVar % Values = 0.0_dp

  sname = ListGetString( SolverParams,'Edge Residual Name',Found)
  IF(.NOT. Found) sname = "ams res"
  EdgeResVar => VariableGet( Mesh % Variables, sname, ThisOnly=.TRUE.,UnfoundFatal=.TRUE.)

  n0 = MAXVAL(EdgeResVar % Perm(1:Mesh % NumberOfNodes))
  IF(n0 > 0) THEN
    CALL Info(Caller,'Number of nodal dofs: '//I2S(n0),Level=20)
  END IF
  
  EdgeSolverParams => GetSolverParams(EdgeResVar % Solver)
  CALL EdgeElementStyle(EdgeSolverParams, PiolaVersion, SecondFamily, SecondOrder, Check = .TRUE.)
  IF (SecondOrder) CALL Fatal(Caller, 'The lowest-order edge basis must be assumed') 

  IF (.NOT. ASSOCIATED(Proj)) THEN
    CALL Info(Caller,'Creating projection matrix to map a nodal solution into vector element space', Level=10)
    SkipFaces = ListGetLogical( SolverParams,'Skip Faces in Projection',Found ) 
    CALL NodalToNedelecInterpolation_GlobalMatrix(Mesh, Avar, EdgeSolVar, Proj, cdim=3, &
        SkipFaces = SkipFaces, NodalOffset = n0)
    IF(InfoActive(20)) THEN
      CALL VectorValuesRange(Proj % Values,SIZE(Proj % Values),'Proj')       
    END IF       
  END IF

  ExtrudedSol = ListGetLogical( SolverParams,'Extruded Solution',Found ) 

  ! Now EdgeResVar represents the residual with respect
  ! to the basis for H(curl). We need to apply a transformation so that
  ! we may solve the residual correction equation by using the nodal basis.
  !-----------------------------------------------------------------------------
  CALL Info(Caller,'Using Transposed Projection Matrix: H(curl) -> H1', Level=10)
  CALL CRS_TransposeMatrixVectorMultiply(Proj, EdgeResVar % Values, allrhs )           
   
  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(Allrhs,SIZE(Allrhs),'allrhs')       
  END IF
  
  ! Potentially create a mask that avoids residual values being applied on the mortar BC. 
  IF(.NOT. Visited ) THEN
    n = SIZE(SVar % Values) / SVar % dofs
    ALLOCATE(NodeSkip(n))    
    NodeSkip = .FALSE.
    CALL CreateNodeSkipMask(NodeSkip,SVar)
    n = COUNT(NodeSkip)
    IF(n==0) DEALLOCATE(NodeSkip)
  END IF

  ! By construction do not apply any residual to the mortar boundary. 
  IF(ASSOCIATED(NodeSkip)) THEN
    DO i=1,SIZE(NodeSkip)
      IF(NodeSkip(i)) THEN
        allrhs(dofs*(i-1)+1:dofs*i) = 0.0_dp
      END IF
    END DO
  END IF

  IF(IsMonolithic) THEN
    ! If we use N-T coordinate system to make periodic/rotational BC's easier than we must map the
    ! original residual vector into N-T system. 
    
    IF( Solver % NormalTangential % NormalTangentialNOFNodes > 0 ) THEN
      CALL Info(Caller,'Mapping residual vector into normal-tangential system',Level=10)
      CALL RotateNtVector( allrhs, Solver )
    END IF

#if 0 
    IF(ListCheckPrefixAnyBodyForce( Model,'Test Load') ) THEN
      CALL SetTestRhs()    
      NodeResVar % Values = allrhs
    END IF
#endif

    Solver % Matrix % rhs = allrhs
    
    IF(ALLOCATED(Solver % Matrix % ConstrainedDOF ) ) &
        Solver % Matrix % ConstrainedDOF = .FALSE.
    CALL DefaultDirichletBCs()
    
    Norm = DefaultSolve()    
  ELSE
    DO compi = 1, comps
      IF(ExtrudedSol .AND. compi /= comps ) CYCLE
      
      A % Values = A % BulkValues
      IF( ns == 1 ) THEN      
        b = allrhs(compi::comps)
      ELSE
        ! Picking a stride is not so easy when we want to pick a pair of components from a set of six. 
        b(1::2) = allrhs(2*compi-1::2*comps)
        b(2::2) = allrhs(2*compi::2*comps)
      END IF
        
      ! Nullify the previous Dirichlet conditions to be on the safe side.
      IF(ALLOCATED(A % ConstrainedDOF ) ) A % ConstrainedDOF = .FALSE.

      ! Setting Dirichlet conditions is not possible with the Default routine.
      ! We can still pick the name of the component but just need to potentially
      ! cycle over the Re and Im parts.
      DO dof=1,ns
        sname = ComponentName(AVar % Name,ns*(compi-1)+dof)

        pVar => VariableGet( Model % Variables, sname)
        
        CALL SetDirichletBoundaries( CurrentModel, A, b, sname, & 
            dof, ns, SVar % Perm )
      END DO

      CALL EnforceDirichletConditions( Solver, A, b )
      
      Norm = DefaultSolve()

      IF( ns == 1 ) THEN
        pVar % Values = SVar % Values
      ELSE
        AVar % Values(2*compi-1::2*comps) = SVar % Values(1::2)
        AVar % Values(2*compi::2*comps) = SVar % Values(2::2)        
      END IF
    END DO
    A % Values = A % BulkValues
  END IF

  IF(IsMonolithic .OR. ExtrudedSol ) THEN
    CALL ListAddLogical( SolverParams,'Linear System Refactorize',.FALSE.) 
  END IF
  CALL ListAddLogical( SolverParams,'Mortar BCs Fixed',.TRUE.)
  
  CALL Info(Caller,'Projecting nodal solution to vector element space', Level=20)
  CALL CRS_MatrixVectorMultiply(Proj, Avar % Values, EdgeSolVar % Values ) 

  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(Avar % Values,SIZE(Avar % Values),'VecPotNodal')       
    CALL VectorValuesRange(EdgeSolVar % Values,SIZE(EdgeSolVar % Values),'VecPotEdge')       
  END IF

  Visited = .TRUE.
  CALL Info(Caller,'Auxiliary space nodal vector solution finished!',Level=10)

  
CONTAINS

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

#if 0 
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
          FORCE(i:3*n:3) = FORCE(i:3*n:3) + Weight * LoadAtIP * Basis(1:n)
        END DO
      END DO

      CALL DefaultUpdateForce(FORCE,Element,Solver)

    END DO
    
  END SUBROUTINE SetTestRhs
#endif
  
END SUBROUTINE AmsVectorSolver
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE AmsScalarSolver_Init( Model,Solver,dt,Transient ) ! {{{
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
  LOGICAL :: Found
  CHARACTER(*), PARAMETER :: Caller = 'AmsScalarSolver_Init'

  Params => GetSolverParams()
  CALL ListAddLogical( Params,'AMS Scalar Solver',.TRUE.)
  CALL ListAddNewLogical( Params,'Apply Mortar BCs',.TRUE.)
  CALL ListAddNewLogical( Params,'Use Global Mass Matrix',.TRUE.)
  CALL ListAddNewString( Params,'Exec Solver','never')
  CALL ListAddNewLogical( Params,'Skip Compute Nonlinear Change',.TRUE.)
  CALL ListAddNewInteger( Params,'Nonlinear System Max Iterations', 1)
  
  IF( ListGetLogical( Params,'Linear System Complex', Found ) ) THEN
    CALL ListAddNewString( Params,'Variable','amss[amss re:1 amss im:1]' )
  ELSE
    CALL ListAddNewString( Params,'Variable','amss' )
  END IF
  
!------------------------------------------------------------------------------
END SUBROUTINE AmsScalarSolver_Init ! }}}
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Gradient part of the preconditioner.
!------------------------------------------------------------------------------
SUBROUTINE AmsScalarSolver( Model,Solver,dt,Transient ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
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
  REAL(KIND=dp) :: Norm
  INTEGER :: i,j,k,n, nb, nd, t, dof
  TYPE(Mesh_t),   POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: Vvar, EdgeSolVar, EdgeResVar
  TYPE(Matrix_t), POINTER, SAVE :: Proj => NULL()
  LOGICAL :: SecondFamily, SecondOrder, PiolaVersion
  TYPE(ValueList_t), POINTER :: EdgeSolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
  LOGICAL, SAVE :: Visited = .FALSE., SkipFaces, IsComplex
  LOGICAL, POINTER, SAVE :: NodeSkip(:)
  TYPE(Matrix_t), POINTER :: A
  CHARACTER(*), PARAMETER :: Caller = 'AmsScalarSolver'
  
!------------------------------------------------------------------------------
  
  Mesh => Solver % Mesh 
  SolverParams => Solver % Values
  VVar => Solver % Variable
  A => Solver % Matrix
  
  CALL Info( Caller,'-------------------------------------------------------', Level=10 )
  CALL Info( Caller,'Solving preconditioning equation for AMS scalar', Level=6 )
  CALL Info( Caller,'-------------------------------------------------------', Level=10 )

  IsComplex = ListGetLogical( SolverParams, 'Linear System Complex', Found )

  IF(VVar % dofs > 2) CALL Fatal(Caller,'Solution size should be <=2!')
  VVar % Values = 0.0_dp

  sname = ListGetString( SolverParams, 'Edge Update Name', Found)
  IF(.NOT. Found) sname = "ams update"
  EdgeSolVar => VariableGet(Mesh % Variables, sname, ThisOnly=.TRUE.,UnfoundFatal=.TRUE.)
  EdgeSolVar % Values = 0.0_dp
  
  sname = ListGetString( SolverParams,'Edge Residual Name',Found)
  IF(.NOT. Found) sname = "ams res"
  EdgeResVar => VariableGet( Mesh % Variables, sname, ThisOnly=.TRUE.,UnfoundFatal=.TRUE.)
  
  EdgeSolverParams => GetSolverParams(EdgeResVar % Solver)
  CALL EdgeElementStyle(EdgeSolverParams, PiolaVersion, SecondFamily, SecondOrder, Check = .TRUE.)
  IF (SecondOrder) CALL Fatal(Caller, 'The lowest-order edge basis must be assumed') 

  IF (.NOT. ASSOCIATED(Proj)) THEN
    CALL Info(Caller,'Creating projection matrix to map a nodal solution into gradient space', Level=10)
    SkipFaces = ListGetLogical( SolverParams,'Skip Faces in Projection',Found ) 
    CALL NodalGradientToNedelecInterpolation_GlobalMatrix(Mesh, VVar, EdgeResVar, Proj)
  END IF
    
  ! Now EdgeResVar represents the residual with respect
  ! to the basis for H(curl). We need to apply a transformation so that
  ! we may solve the residual correction equation by using the nodal basis.
  !-----------------------------------------------------------------------------
  CALL Info(Caller,'Using Transposed Projection Matrix: H(curl) -> Grad', Level=10)
  CALL CRS_TransposeMatrixVectorMultiply(Proj, EdgeResVar % Values, A % rhs ) 

  ! Potentially create a mask that avoids residual values being applied on the mortar BC. 
  IF(.NOT. Visited  ) THEN
    n = SIZE(A % rhs)
    ALLOCATE(NodeSkip(n))    
    NodeSkip = .FALSE.
    CALL CreateNodeSkipMask(NodeSkip,VVar)
    n = COUNT(NodeSkip)
    IF(n==0) DEALLOCATE(NodeSkip)
  END IF

  ! By construction do not apply any residual to the mortar boundary. 
  IF(ASSOCIATED(NodeSkip)) THEN
    WHERE(NodeSkip)
      A % rhs = 0.0_dp
    END WHERE
  END IF

  IF(ALLOCATED(A % ConstrainedDOF ) ) A % ConstrainedDOF = .FALSE.
  DO dof=1,VVar % dofs
    sname = ComponentName(VVar,dof)    
    CALL SetDirichletBoundaries( CurrentModel, A, A % rhs, sname, & 
        dof, VVar % dofs, VVar % Perm )
  END DO
  CALL EnforceDirichletConditions( Solver, A, A % rhs )

  Norm = DefaultSolve()    

  CALL ListAddLogical( SolverParams,'Mortar BCs Fixed',.TRUE.)
    
  CALL Info(Caller,'Projecting nodal solution to vector element space', Level=10)
  CALL CRS_MatrixVectorMultiply(Proj, VVar % Values, EdgeSolVar % Values ) 

  BLOCK
    LOGICAL, POINTER :: SkipMask(:)
    LOGICAL :: DoMask
    DoMask = ListGetLogical( SolverParams,'Zero Mortar Fix',Found)
    IF(DoMask) THEN
      n = SIZE(EdgeSolVar % Values)
      ALLOCATE(SkipMask(n))
      SkipMask = .FALSE.
      CALL CreateEdgeSkipMask(SkipMask)
      WHERE(SkipMask) EdgeSolVar % Values = 0.0_dp
      DEALLOCATE(SkipMask)
    END IF
  END BLOCK

  
  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(Vvar % Values,SIZE(Vvar % Values),'ScalarPotNodal')       
    CALL VectorValuesRange(EdgeSolVar % Values,SIZE(EdgeSolVar % Values),'VecPotEdge')       
  END IF

  Visited = .TRUE.
  CALL Info(Caller,'Auxiliary space nodal scalar solution finished!',Level=10)
    
  
END SUBROUTINE AmsScalarSolver
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE AmsVSolver_Init( Model,Solver,dt,Transient ) ! {{{
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
  LOGICAL :: Found
  CHARACTER(*), PARAMETER :: Caller = 'AmsScalarSolver_Init'

  Params => GetSolverParams()
  CALL ListAddLogical( Params,'AMS V Solver',.TRUE.)
  CALL ListAddNewLogical( Params,'Apply Mortar BCs',.TRUE.)
  CALL ListAddNewLogical( Params,'Use Global Mass Matrix',.TRUE.)
  CALL ListAddNewString( Params,'Exec Solver','never')
  CALL ListAddNewLogical( Params,'Skip Compute Nonlinear Change',.TRUE.)
  CALL ListAddNewInteger( Params,'Nonlinear System Max Iterations', 1)
  
  IF( ListGetLogical( Params,'Linear System Complex', Found ) ) THEN
    CALL ListAddNewString( Params,'Variable','amsv[amsv re:1 amsv im:1]' )
  ELSE
    CALL ListAddNewString( Params,'Variable','amsv' )
  END IF
  
!------------------------------------------------------------------------------
END SUBROUTINE AmsVSolver_Init ! }}}
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Preconitioner for the V part of the AV equation when AMS is used.
!------------------------------------------------------------------------------
SUBROUTINE AmsVSolver( Model,Solver,dt,Transient ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
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
  REAL(KIND=dp) :: Norm
  INTEGER :: i,j,k,n, nb, nd, t, dof
  TYPE(Mesh_t),   POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: Vvar, EdgeSolVar, EdgeResVar
  TYPE(ValueList_t), POINTER :: EdgeSolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
  LOGICAL, SAVE :: Visited = .FALSE., IsComplex
  LOGICAL, POINTER, SAVE :: NodeSkip(:)
  TYPE(Matrix_t), POINTER :: VMat, AVMat
  CHARACTER(*), PARAMETER :: Caller = 'AmsVSolver'  
!------------------------------------------------------------------------------
  
  Mesh => Solver % Mesh 
  SolverParams => Solver % Values
  VVar => Solver % Variable
  VMat => Solver % Matrix
  
  CALL Info( Caller,'-------------------------------------------------------', Level=10 )
  CALL Info( Caller,'Solving preconditioning equation for AMS V', Level=6 )
  CALL Info( Caller,'-------------------------------------------------------', Level=10 )

  IsComplex = ListGetLogical( SolverParams, 'Linear System Complex', Found )

  IF(VVar % dofs > 2) CALL Fatal(Caller,'Solution size should be <=2!')
  VVar % Values = 0.0_dp

  sname = ListGetString( SolverParams, 'Edge Update Name', Found)
  IF(.NOT. Found) sname = "ams update"
  EdgeSolVar => VariableGet(Mesh % Variables, sname, ThisOnly=.TRUE.,UnfoundFatal=.TRUE.)
  EdgeSolVar % Values = 0.0_dp
  
  sname = ListGetString( SolverParams,'Edge Residual Name',Found)
  IF(.NOT. Found) sname = "ams res"
  EdgeResVar => VariableGet( Mesh % Variables, sname, ThisOnly=.TRUE.,UnfoundFatal=.TRUE.)
  
  EdgeSolverParams => GetSolverParams(EdgeResVar % Solver)
  AVMat => EdgeResVar % Solver % Matrix
  
  CALL PickNodalSubmatrix(AVMat,VMat,EdgeSolVar,VVar)

  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(VMat % Values,SIZE(VMat % Values),'Vmat values')       
  END IF

  
  ! Now EdgeResVar has first residual related to "V" and then to "A".
  ! We pick just the "V" values.
  n = VMat % NumberOfRows
  VMat % rhs(1:n) = EdgeResVar % Values(1:n)

  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(VMat % rhs,SIZE(VMat % Rhs),'Vrhs')       
  END IF

  
  ! Potentially create a mask that avoids residual values being applied on the mortar BC. 
  IF(.NOT. Visited  ) THEN
    ALLOCATE(NodeSkip(n))    
    NodeSkip = .FALSE.
    CALL CreateNodeSkipMask(NodeSkip,VVar)
    n = COUNT(NodeSkip)
    IF(n==0) DEALLOCATE(NodeSkip)
  END IF

  ! By construction do not apply any residual to the mortar boundary. 
  IF(ASSOCIATED(NodeSkip)) THEN
    WHERE(NodeSkip)
      VMat % rhs = 0.0_dp
    END WHERE
  END IF

  ! No Dirichlet conditions etc. need to be set since the whole matrix is inherited from
  ! the large matrix.
  Norm = DefaultSolve()    
  
  CALL ListAddLogical( SolverParams,'Mortar BCs Fixed',.TRUE.)
    
  EdgeSolVar % Values(1:n) = VVar % Values(1:n)
  
  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(Vvar % Values,SIZE(Vvar % Values),'VNodal')       
  END IF

  Visited = .TRUE.
  CALL Info(Caller,'Auxiliary space nodal V solution finished!',Level=10)

CONTAINS

  ! Picks the part related to nodal dofs from AV matrix.
  ! Assumes that dofs have not been reordered!!
  !-----------------------------------------------------
  SUBROUTINE PickNodalSubmatrix(TotMat,SubMat,TotVar,SubVar)
    TYPE(Matrix_t) :: TotMat, SubMat
    TYPE(Variable_t) :: TotVar, SubVar

    INTEGER, ALLOCATABLE :: TotInvPerm(:),SubInvPerm(:)    
    INTEGER :: i,i1,i2,j1,j2,k1,k2,m,k
    
    SubMat % Values = 0.0_dp

    m = SIZE(TotVar % Perm)
    ALLOCATE(TotInvPerm(m),SubInvPerm(m))
    TotInvPerm = 0
    SubInvPerm = 0
    
    m = Mesh % NumberOfNodes
    DO i=1,m
      i1 = SubVar % Perm(i)
      i2 = TotVar % Perm(i)
      IF(i1>0) SubInvPerm(i1) = i
      IF(i2>0) TotInvPerm(i2) = i
    END DO

    k = 0
    DO i=1,Mesh % NumberOfNodes
      i1 = SubVar % Perm(i)
      i2 = TotVar % Perm(i)
      IF(i1==0 .OR. i2==0) CYCLE

      DO j1=SubMat % Rows(i1),SubMat % Rows(i1+1)-1
        k1 = SubInvPerm(SubMat % Cols(j1)) 
        DO j2=TotMat % Rows(i2),TotMat % Rows(i2+1)-1
          k2 = TotInvPerm(TotMat % Cols(j2)) 
          IF(k1==k2) THEN
            SubMat % Values(j1) = TotMat % Values(j2)
            k = k+1
          END IF
        END DO
      END DO
    END DO

    IF(k<SIZE(SubMat % Values)) THEN
      CALL Fatal(Caller,'Not all values found for matrix!')
    END IF
    
  END SUBROUTINE PickNodalSubmatrix
    
  
END SUBROUTINE AmsVSolver
!------------------------------------------------------------------------------



!> \}

