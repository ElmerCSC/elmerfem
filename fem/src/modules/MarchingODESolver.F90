!/******************************************************************************
! *
! *  Module for solving 1dof Ordinary Differential Equation (ODE) solver on moving
! *  coordinates. The mesh is assumed to be extruded and we assume constant draw
! *  velocity. Hence we can march the ODE at material point in time. 
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 4.12.2019
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{
 

!------------------------------------------------------------------------------
!> Initialization for the primary solver
!------------------------------------------------------------------------------
SUBROUTINE MarchingODESolver_init( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  
  TYPE(Model_t) :: Model    
  REAL(KIND=dp) :: dt       
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'MarchingODESolver_init'
  TYPE(ValueList_t), POINTER :: Params

  Params => GetSolverParams()

  ! By construction time derivative order must be one!
  CALL ListAddNewInteger( Params, 'Time derivative order', 1 )

  ! Global mass matrix allows us to use library routines for time integration
  CALL ListAddNewLogical( Params, 'Use Global Mass Matrix', .TRUE. )

  ! We initialize our own matrix structures
  CALL ListAddNewLogical( Params, 'No Matrix', .TRUE. )
    
!------------------------------------------------------------------------------
END SUBROUTINE MarchingODESolver_init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> This solver combines a simple solver for ODE in time with marching along
!> structured mesh.
!------------------------------------------------------------------------------
SUBROUTINE MarchingODESolver( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model            !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt               !< Timestep size for time dependent simulations
  LOGICAL :: Transient              !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'MarchingODESolver'
  LOGICAL :: Found
  REAL(KIND=dp) :: Norm, Change, dz, dtime, velo, NonLinTol  
  INTEGER :: i,j,n,iter,MaxIter,TimeOrder,nsize,nnodes,BotNodes,layer
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t), POINTER :: PSolver
  INTEGER, POINTER :: BotPointer(:), UpPointer(:)
  INTEGER, POINTER :: BotPerm(:),InvPerm(:),PrevInvPerm(:),OdePerm(:),MaskPerm(:)
  INTEGER :: NumberOfLayers, NoBCNodes
  TYPE(Variable_t), POINTER :: ExtVar, Var1D, Var3D
  TYPE(ValueList_t), POINTER :: Material
  LOGICAL :: MaskExist, ParabolicModel, RequireBC
  REAL(KIND=dp), POINTER :: Coord(:), OdeValues(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: TimeMethod, VarName
  LOGICAL, ALLOCATABLE :: BCNode(:)
  TYPE(Matrix_t), POINTER :: Amat
  REAL(KIND=dp), POINTER :: xvec(:),x0vec(:),dxvec(:)
  
  LOGICAL, SAVE :: Initialized = .FALSE.
!------------------------------------------------------------------------------

  SAVE :: BotPointer, UpPointer, &
      BotPerm, InvPerm, PrevInvPerm, OdePerm, MaskPerm, &
      NumberOfLayers, ExtVar, Var1D, OdeValues, BotNodes, &
      TimeMethod, ParabolicModel, RequireBC, x0vec, dxvec
  
  CALL Info(Caller,'-----------------------------------------------------',Level=6)
  CALL Info(Caller,'Solving ODE on moving coordinates in structured mesh',Level=4)

  Params => GetSolverParams()
  Mesh => Solver % Mesh
  PSolver => Solver

  IF( .NOT. Initialized ) THEN
    CALL Info(Caller,'Initializing structured mesh and ODE structures',Level=7)

    CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar = ExtVar, &
        BotNodePointer = BotPointer, UpNodePointer = UpPointer, &
        NumberOfLayers = NumberOfLayers )
        
    CALL Info(Caller,'Number of element layers: '//TRIM(I2S(NumberOfLayers)),Level=7)

    MaskExist = ASSOCIATED( ExtVar % Perm ) 
    IF( MaskExist ) MaskPerm => ExtVar % Perm
    Coord => ExtVar % Values
    nsize = MIN( SIZE( Coord ), Mesh % NumberOfNodes )
    nnodes = Mesh % NumberOfNodes

    ! We may choose only to apply the ODE to BC nodes
    RequireBC = ListGetLogical( Params,'Apply BCs Only',Found )
    IF( RequireBC ) THEN
      CALL MarkBCNodes( Mesh,BCNode,NoBCNodes)
      IF(NoBCNodes == 0 ) RequireBC = .FALSE.
      CALL Info(Caller,'Number of BC nodes: '//TRIM(I2S(NoBCNodes)),Level=7)
    END IF

    ! Create the permutation using the bottom layer
    BotNodes = 0
    ALLOCATE( BotPerm( Mesh % NumberOfNodes ) )
    BotPerm = 0
    DO i=1,Mesh % NumberOfNodes
      IF( RequireBC ) THEN
        IF( .NOT. BcNode(i) ) CYCLE
      END IF
      j = i
      IF( MaskExist ) THEN
        j = MaskPerm(i)
        IF( j == 0 ) CYCLE
      END IF
      IF(BotPointer(j) == i) THEN
        BotNodes = BotNodes + 1
        BotPerm(i) = BotNodes
      END IF
    END DO
    n = BotNodes
    ALLOCATE( InvPerm(n), PrevInvPerm(n) )

    CALL Info(Caller,'Number of bottom nodes: '//TRIM(I2S(n)),Level=7)

    ! Create the matrix for one layer
    TimeOrder = ListGetInteger( Params,'Time Derivative Order',Found) 
    Solver % Matrix => CreateDiagMatrix( Model, Solver, n, TimeOrder )
    ALLOCATE( Solver % Matrix % rhs(n) )

    ! Create the variable for one layer
    ALLOCATE( OdePerm(n), OdeValues(n) )
    DO i=1,n
      OdePerm(i) = i
    END DO
    OdeValues = 0.0_dp    
    CALL VariableAdd( Mesh % Variables,Mesh,Solver,'ODEVar',1,OdeValues,&
        OdePerm,Output=.FALSE.)
    Var1D => VariableGet( Mesh % Variables,'ODEVar')
    ALLOCATE( Var1D % PrevValues(n,1) )
    Var1D % PrevValues = 0.0_dp

    ! Allocate some vectors to study convergence 
    ALLOCATE( x0vec(n), dxvec(n) )
    x0vec = 0.0_dp
    dxvec = 0.0_dp

    ! This is just to provide correct info for timestepping
    Solver % DoneTime = 1
    Solver % Order = 1
    
    Initialized = .TRUE.
    CALL Info(Caller,'Initialization done',Level=10)
  END IF
  
  ! In principle we could use the same solver to advect many fields on
  ! same extruded mesh. Hence these are not saved. 
  Var3D => Solver % Variable
  VarName = TRIM( Var3D % Name ) 

  ! The variable on the layer
  n = BotNodes
  Solver % Variable => Var1D
  Amat => Solver % Matrix
  xvec => Var1D % Values
    
  CALL Info(Caller,'Solving for variable: '//TRIM(VarName))

  TimeMethod = ListGetString( Params, 'Timestepping Method' )

  ParabolicModel = ListGetLogical( Params,'Parabolic Model',Found )
  IF( ParabolicModel ) THEN
    CALL Info(Caller,'Using parabolic growth model',Level=7)
  END IF
  
  velo = ListGetCReal( Params,'Draw Velocity',UnfoundFatal=.TRUE.)
  
  MaxIter = GetInteger( Params,'Nonlinear System Max Iterations',Found )
  IF(.NOT. Found) MaxIter = 1

  NonLinTol = GetCReal( Params,'Nonlinear System Convergence Tolerance',Found )

  Material => FirstExtrudedMaterial()

  !------------------------------------------------------------------------
  DO layer=0,NumberOfLayers

    CALL Info(Caller,'Solving for layer: '//TRIM(I2S(layer)),Level=8)
    
    ! First layer is determined by the initial conditions (=boundary conditions)
    IF( layer == 0 ) THEN
      DO i=1,Mesh % NumberOfNodes
        j = BotPerm(i)
        IF( j > 0 ) THEN
          InvPerm(j) = i
        END IF
      END DO
      IF( ParabolicModel ) THEN
        xvec(1:n) = 0.5_dp * Var3D % Values(Var3D % Perm(InvPerm))**2 
      ELSE
        xvec(1:n) = Var3D % Values(Var3D % Perm(InvPerm))       
      END IF
      Var1D % PrevValues(1:n,1) = xvec(1:n)
      CYCLE
    END IF

    ! Find the next level of nodes, and remember the previous one. 
    PrevInvPerm = InvPerm
    InvPerm = UpPointer(PrevInvPerm)

    ! If we have parabolic model the equation is on-the-fly written for u^2/2 instead of u
    IF( ParabolicModel ) THEN
      x0vec(1:n) = 0.5_dp * Var3D % Values(Var3D % Perm(InvPerm))**2 
    ELSE
      x0vec(1:n) = Var3D % Values(Var3D % Perm(InvPerm))       
    END IF

    ! This sets the timestep assuming that all nodes are extruded equally.
    ! Hence this only applied to cartesian drawing. 
    dz = Mesh % Nodes % z(InvPerm(1)) - Mesh % Nodes % z(PrevInvPerm(1))
    dtime = dz / velo

    ! We may have iteration if the ODE is nonlinear.
    ! However, more typical could be to iterate over coupled systems. 
    DO iter = 1, MaxIter 

      IF( MaxIter > 1 ) THEN
        CALL Info('ODESolver','Nonlinear iteration: '//TRIM(I2S(iter)),Level=20)
        x0vec = xvec
      END IF
      
      CALL InitializeToZero( Amat, Amat % RHS )
      
      CALL OdeAssembly()
          
      CALL Add1stOrderTime_CRS( Amat, Amat % rhs, dt, PSolver )
      
      ! We don't really need any linear solver for this as there is no coupling among dofs.
      xvec = Amat % rhs / Amat % Values 
      
      Norm = SQRT(SUM(xvec*xvec))

      dxvec = xvec-x0vec

      ! For the 1st iteration the error corresponds to the error with respect to previous solution.
      ! For 2nd etc. iteration the error is of the nonlinear iteration. 
      Change = SQRT(SUM(dxvec*dxvec))/Norm

      IF( Change < NonLinTol ) EXIT
    END DO
              
    PRINT *,'Layer:',layer,dtime,Norm,Change

    IF( ParabolicModel ) THEN
      Var3D % Values(Var3D % Perm(InvPerm)) = SQRT( 2 * xvec ) 
    ELSE      
      Var3D % Values(Var3D % Perm(InvPerm)) = xvec
    END IF
    Var1D % PrevValues(1:n,1) = xvec
  END DO

  ! Take back the pointer, this tampering would allow us to use DefaultSolve()
  ! if we would so choose. 
  Solver % Variable => Var3D
  
  CALL Info(Caller,'All done',Level=5)
  CALL Info(Caller,'-----------------------------------------------------',Level=6)


CONTAINS

  ! Finds pointer to the extruded material.
  ! Note that we assume that this is the only material being extruded.
  !--------------------------------------------------------------------------------
  FUNCTION FirstExtrudedMaterial() RESULT ( Material ) 
    TYPE(Valuelist_t), POINTER :: Material
    TYPE(Element_t), POINTER :: Element
    INTEGER :: t
    
    Material => NULL()
    DO t=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)
      IF( MaskExist ) THEN
        IF( ANY( MaskPerm( Element % NodeIndexes ) == 0 ) ) CYCLE
      END IF
      IF( ANY( Var3D % Perm( Element % NodeIndexes ) == 0 ) ) CYCLE
      Material => GetMaterial(Element)
      EXIT
    END DO

    IF(.NOT. ASSOCIATED( Material ) ) THEN
      CALL Fatal(Caller,'Could not set material for extruded domain!')
    END IF
    
  END FUNCTION FirstExtrudedMaterial
  

  ! Creates the local matrix equation for the ODY before time integration.
  ! We use namespace in order to allow the same solver to be used for
  ! several fields. The equation is of type:
  ! c_t*du/dt + c_r*u = s. 
  !-----------------------------------------------------------------------
  SUBROUTINE OdeAssembly()
    
    TYPE(Matrix_t), POINTER :: A
    
    A => Solver % Matrix
    n = Solver % Matrix % NumberOfRows
    
    A % Rhs(1:n) = ListGetReal( Material,TRIM(VarName)//': Source',&
        n,InvPerm,UnFoundFatal=.TRUE. )    
    A % Values(1:n) = ListGetReal( Material,TRIM(VarName)//': Reaction Coefficient',&
        n,InvPerm,Found ) 
    A % MassValues(1:) = ListGetReal( Material,TRIM(VarName)//': Time Derivative Coefficient',&
        n,invPerm,Found )
    IF(.NOT. Found ) THEN
      A % MassValues(1:n) = 1.0_dp
    END IF

  END SUBROUTINE OdeAssembly

!------------------------------------------------------------------------------
END SUBROUTINE MarchingODESolver
!------------------------------------------------------------------------------
!> \}
