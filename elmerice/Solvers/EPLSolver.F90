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
! *  Authors: Basile de Fleurian
! *  Email:   basile.defleurian@uci.edu; gagliar@lgge.obs.ujf-grenoble.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 02 Jully 2013 
! * 
! *****************************************************************************
!> Limited solver to compute the water head in a sediment layer
!> Exported Variables EPLHead (dof=1)
!> The solver is treated on DIM-1, wit DIM-2 boundary conditions
! *****************************************************************************

RECURSIVE SUBROUTINE EPLSolver( Model,Solver,Timestep,TransientSimulation )
  !------------------------------------------------------------------------------

  USE DiffuseConvective
  USE DiffuseConvectiveGeneral
  USE Differentials
  USE MaterialModels
  USE DefUtils
  !------------------------------------------------------------------------------
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  !******************************************************************************
  !
  !  Solve the diffusion equation in a porous layer within a domain given by a 
  !  mask
  !
  !  ARGUMENTS:
  !
  !  TYPE(Model_t) :: Model,  
  !     INPUT: All model information (mesh,materials,BCs,etc...)
  !
  !  TYPE(Solver_t) :: Solver
  !     INPUT: Linear equation solver options
  !
  !  REAL(KIND=dp) :: Timestep
  !     INPUT: Timestep size for time dependent simulations
  !
  !  LOGICAL :: TransientSimulation
  !     INPUT: Steady state or transient simulation
  !
  !******************************************************************************

  !------------------------------------------------------------------------------
  !    External variables
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER :: Systemmatrix
  TYPE(Element_t), POINTER :: Element
  TYPE(Solver_t), POINTER :: PointerToSolver
  TYPE(Variable_t), POINTER :: EPLHeadSol, Piping, IDSRes, VarEPLHead
  TYPE(ValueList_t), POINTER :: Constants, SolverParams, Equation, &
       Material, BodyForce, BC
  TYPE(Nodes_t) :: ElementNodes

  REAL(KIND=dp), POINTER :: EPLHead(:), OpenEPL(:), IDSResInput(:), &
       ForceVector(:), Hwrk(:,:,:), EPLHeadHomologous(:)
  REAL(KIND=dp), ALLOCATABLE :: Transmitivity(:,:,:) 
  REAL(KIND=dp), ALLOCATABLE :: g(:,:), Nochange(:,:), &
       MASS(:,:), STIFF(:,:)
  REAL(KIND=dp), ALLOCATABLE :: UpperLimit(:), EPLComp(:), &
       Porosity(:), Gravity(:), EPLThick(:) , Density(:), &
       StoringCoef(:), Viscosity(:), C0(:), C1(:), Zero(:), &
       VNull(:), Work(:), Pressure(:), TransferCoeff(:), &
       EPLHeadExt(:), OldValues(:), OldRHS(:),StiffVector(:), &
       IDSHead(:), EPLToIDS(:),DummyRealArray(:), &
       FORCE(:),LOAD(:), TimeForce(:) 
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: LinearTol, NonlinearTol, &
       totat, totst, at, at0, st, &
       WatComp, Norm, PrevNorm, RelativeChange
#else
  REAL(KIND=dp) :: LinearTol, NonlinearTol, &
       totat, totst, at, at0, st, CPUTime, RealTime, &
       WatComp, Norm, PrevNorm, RelativeChange
#endif
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, VariableName, &
       IDSHeadName, IDSResName

  INTEGER, POINTER :: EPLHeadPerm(:), PipingPerm(:), &
       IDSResPerm(:), HomolPerm(:)
  INTEGER :: NonlinearIter, LocalNodes, &
       PenIter, k, N, M, L, t, i, j, istat, &
       body_id, material_id, bf_id, bc_id

  LOGICAL, ALLOCATABLE :: ActiveEPL(:)
  LOGICAL :: Found = .FALSE.,UnFoundFatal=.TRUE., &
       Stabilize = .TRUE. ,&
       UseBubbles = .FALSE., &
       AllocationsDone = .FALSE., &
       FluxBC = .FALSE., &
       IsPeriodicBC = .FALSE., &
       ApplyDirichlet = .FALSE., &
       DummyLogical = .FALSE.

  SAVE ElementNodes,    &
       Hwrk,            &
       Transmitivity,   &
       g,               &
       Nochange,        &
       MASS,            &
       STIFF,           &
       UpperLimit,      &
       EPLComp,         &
       Porosity,        &
       Gravity,         &
       EPLThick,        &
       Density,         &
       StoringCoef,     &
       Viscosity,       &
       C0,              &
       C1,              &
       Zero,            &
       VNull,           &
       Work,            &
       Pressure,        &
       TransferCoeff,   &
       EPLHeadExt,      &
       OldValues,       &
       OldRHS,          &
       StiffVector,     &
       IDSHead,         &
       EPLToIDS,        &
       DummyRealArray,  &
       FORCE,           &
       LOAD,            &
       TimeForce,       &
       ActiveEPL,       &
       SolverName,      &
       VariableName,    &
       LinearTol,       &
       NonlinearTol,    &
       AllocationsDone, &
       DummyLogical


  SystemMatrix => Solver % Matrix
  ForceVector => Solver % Matrix % RHS

  !------------------------------------------------------------------------------
  !    Get variables needed for solution
  !------------------------------------------------------------------------------
  SolverName = 'EPLSolver ('// TRIM(Solver % Variable % Name) // ')'
  VariableName = TRIM(Solver % Variable % Name)

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  PointerToSolver => Solver

  EPLHeadSol => Solver % Variable
  EPLHeadPerm  => EPLHeadSol % Perm
  EPLHead => EPLHeadSol % Values

  LocalNodes = COUNT( EPLHeadPerm .GT. 0 )
  IF ( LocalNodes .LE. 0 ) RETURN  

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
     N = Solver % Mesh % MaxElementNodes
     M = Model % Mesh % NumberOfNodes
     K = SIZE( SystemMatrix % Values )
     L = SIZE( SystemMatrix % RHS )

     IF ( AllocationsDone ) THEN  
        DEALLOCATE(            &
             Transmitivity,    &
             g,                &
             Nochange,         &
             MASS,             &
             STIFF,            &
             UpperLimit,       &
             EPLComp,          &
             Porosity,         &
             Gravity,          &
             EPLThick,         &
             Density,          &
             StoringCoef,      &
             Viscosity,        &
             C0,               &
             C1,               &
             Zero,             &
             VNull,            &
             Work,             &
             Pressure,         &
             TransferCoeff,    &
             EPLHeadExt,       &
             OldValues,        &
             OldRHS,           &
             StiffVector,      &
             IDSHead,          &
             EPLToIDS,         &
             DummyRealArray,   &
             FORCE,            &
             LOAD,             &
             TimeForce,        &
             ActiveEPL,        &
             ElementNodes % x, &
             ElementNodes % y, &
             ElementNodes % z)          
     END IF
     ALLOCATE(                   &
          Transmitivity( 3,3,N ),&
          g( 3,N ),              &
          Nochange( 3,N ),       &
          MASS( 2*N,2*N ),       &
          STIFF( 2*N,2*N ),      &
          UpperLimit( M ),       &
          EPLComp( N ),          &
          Porosity( N ),         &
          Gravity( N ),          &
          EPLThick( N ),         &
          Density( N ),          &
          StoringCoef( N ),      &
          Viscosity( N ),        &
          C0( N ),               &
          C1( N ),               &
          Zero( N ),             &
          VNull( N ),            &
          Work( N ),             &
          Pressure( N ),         &
          TransferCoeff( N ),    &
          EPLHeadExt( N ),       &
          OldValues( K ),        &
          OldRHS( L ),           &
          StiffVector( L ),      &
          IDSHead( N ),          &
          EPLToIDS( N ),         &
          DummyRealArray( N ),   &
          FORCE( 2*N ),          &
          LOAD( N ),             &
          TimeForce( 2*N ),      &
          ActiveEPL( M ),        &
          ElementNodes % x( N ), &
          ElementNodes % y( N ), &
          ElementNodes % z( N ), &
          STAT=istat )

     IF ( istat .NE. 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error in EPL Solve' )
     ELSE
        CALL INFO(SolverName, 'Memory allocation done in EPL Solve', level=1 )
     END IF
     AllocationsDone = .TRUE.
     ActiveEPL = .FALSE.
  END IF

  !------------------------------------------------------------------------------
  !    Say hello
  !------------------------------------------------------------------------------
  WRITE(Message,'(A,A)')&
       'Limited diffusion Solver for variable ', VariableName
  CALL INFO(SolverName,Message,Level=1)
  !------------------------------------------------------------------------------
  !    Read Iteration related constants
  !------------------------------------------------------------------------------
  SolverParams => GetSolverParams()

  Stabilize = GetLogical( SolverParams,'Stabilize', Found )
  IF (.NOT. Found) THEN
     Stabilize = .FALSE.
  END IF

  UseBubbles = GetLogical( SolverParams,'Bubbles', Found )
  IF ( .NOT.Found ) THEN
     UseBubbles = .FALSE.
  END IF

  LinearTol = GetConstReal( SolverParams, &
       'Linear System Convergence Tolerance', Found )
  IF ( .NOT.Found ) THEN
     CALL FATAL(SolverName, 'No >Linear System Convergence Tolerance< found')
  END IF

  NonlinearIter = GetInteger( SolverParams, &
       'Nonlinear System Max Iterations', Found )
  IF ( .NOT.Found ) THEN
     NonlinearIter = 1
  END IF

  NonlinearTol = GetConstReal( SolverParams, &
       'Nonlinear System Convergence Tolerance', Found )
  IF ( .NOT.Found ) THEN
     NonlinearTol = 1e-6
  END IF

  !------------------------------------------------------------------------------
  !    Read Physical constants
  !------------------------------------------------------------------------------
  Constants => GetConstants()

  WatComp = GetConstReal( Constants, &
       'Water Compressibility', Found)
  IF ( .NOT.Found ) THEN
     WatComp = 5.04e-4
  END IF

  totat = 0.0d0
  totst = 0.0d0

  !Variable for the sediment Residual
  !----------------------------------
  IDSResName = GetString(SolverParams,'IDS Residual Name', Found)
  IF (.NOT.Found) THEN 
     WRITE(Message,'(a)') 'No Keyword >IDS Residual Name< defined. Using >IDSHead Residual< as default. '
     CALL INFO(SolverName, Message, level=10)
     WRITE(IDSResName,'(A)') 'IDSHead Residual'
  ELSE
     WRITE(Message,'(a,a)') 'Variable Name for IDS Water Head Residual: ', IDSResName
     CALL INFO(SolverName,Message,Level=12)
  END IF
  
  IDSRes => VariableGet( Model % Mesh % Variables,IDSResName,UnFoundFatal=UnFoundFatal)
  IDSResInput => IDSRes % Values
  IDSResPerm => IDSRes % Perm

  VarEPLHead => VariableGet( Model % Mesh % Variables, TRIM(Solver % Variable % Name) // ' Homologous',UnFoundFatal=UnFoundFatal)
  EPLHeadHomologous => VarEPLHead % Values
  HomolPerm => VarEPLHead % Perm

  !Mask Variable
  !----------------------------------------
  Piping => VariableGet( Model % Mesh % Variables,'Open EPL',UnFoundFatal=UnFoundFatal)
  PipingPerm  => Piping % Perm
  OpenEPL => Piping % Values

  !------------------------------------------------------------------------------
  !       Do we use limiters
  !------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  ApplyDirichlet = GetLogical( SolverParams, &
       'Apply Dirichlet', Found)
  IF ( .NOT.Found ) THEN
     ApplyDirichlet = .FALSE.
  END IF
  IF (.NOT.ApplyDirichlet) THEN
     ActiveEPL = .FALSE.
  END IF

  !-----------------------------------------------------------------------------
  ! Get Names for imported variables
  !-----------------------------------------------------------------------------
  IDSHeadName = GetString(SolverParams,'IDS Head Name', Found)
  IF (.NOT.Found) THEN 
     WRITE(Message,'(a)') 'No Keyword >IDS Head Name< defined. Using >IDSHead< as default. '
     CALL INFO(SolverName, Message, level=10)
     WRITE(IDSHeadName,'(A)') 'IDSHead'
  ELSE
     WRITE(Message,'(a,a)') 'Variable Name for IDS Water Head: ', IDSHeadName
     CALL INFO(SolverName,Message,Level=12)
  END IF

  !------------------------------------------------------------------------------
  !       non-linear system iteration loop
  !------------------------------------------------------------------------------
  DO Peniter=1,NonlinearIter
     
     !------------------------------------------------------------------------------
     ! print out some information
     !------------------------------------------------------------------------------
     at  = CPUTime()
     at0 = RealTime()
     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     WRITE( Message,'(A,A,I3,A,I3)') &
          TRIM(Solver % Variable % Name),  ' iteration no.', Peniter,' of ',NonlinearIter
     CALL Info( SolverName, Message, Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, 'Starting Assembly...', Level=4 )

     !------------------------------------------------------------------------------
     ! lets start
     !------------------------------------------------------------------------------
     CALL DefaultInitialize()
     
     !------------------------------------------------------------------------------
     ! write some info on max/min values
     !------------------------------------------------------------------------------
     WRITE(Message,'(a,e13.6,a,e13.6)') &
          'Max/min values of EPL Water Head: ', MAXVAL(EPLHead(:)), &
          '/',MINVAL(EPLHead(:))
     CALL INFO(SolverName,Message,Level=4)

     !------------------------------------------------------------------------------
     body_id = -1
     NULLIFY(Material)

     !------------------------------------------------------------------------------
     ! Check if new nodes are activated and update the force vector accordingly
     !------------------------------------------------------------------------------
      DO i=1,Model % Mesh % NumberofNodes
        IF ((EPLHeadPerm(i).GT.0).AND.&
             (IDSResPerm(i).GT.0).AND.&
             (PipingPerm(i).GT.0)) THEN
           ForceVector(EPLHeadPerm(i)) = -IDSResInput(IDSResPerm(i))

           !Dealing with Mask Values  
           !------------------------
           IF (IDSResInput(IDSResPerm(i)).LT.0.0)THEN
              OpenEPL(PipingPerm(i)) = -1.0
           ELSEIF(ABS(OpenEPL(PipingPerm(i))).LT.1.0)THEN
              OpenEPL(PipingPerm(i)) = 1.0
           END IF
        END IF
     END DO

     !------------------------------------------------------------------------------
     ! Bulk elements Loop
     !------------------------------------------------------------------------------
     DO t=1,Solver % NumberOfActiveElements
        
        !------------------------------------------------------------------------------
        ! Check if this element belongs to a body where scalar equation
        ! should be calculated
        !------------------------------------------------------------------------------
        Element => GetActiveElement(t,Solver)
        N = GetElementNOFNodes(Element)
        CALL GetElementNodes( ElementNodes,Element )

        IF (.NOT.ASSOCIATED(Element)) CYCLE
        
        IF ( Element % BodyId .NE. body_id ) THEN
           Equation => GetEquation()
           IF (.NOT.ASSOCIATED(Equation)) THEN
              WRITE (Message,'(A,I3)') 'No Equation  found for boundary element no. ', t
              CALL FATAL(SolverName,Message)
           END IF
           
           Material => GetMaterial()
           IF (.NOT.ASSOCIATED(Material)) THEN
              WRITE (Message,'(A,I3)') 'No Material found for boundary element no. ', t
              CALL FATAL(SolverName,Message)
           ELSE
              material_id = GetMaterialId( Element, Found)
              IF(.NOT.Found) THEN
                 WRITE (Message,'(A,I3)') 'No Material ID found for boundary element no. ', t
                 CALL FATAL(SolverName,Message)
              END IF
           END IF
        END IF

        !------------------------------------------------------------------------------
        ! Get element material parameters
        !------------------------------------------------------------------------------              
        EPLComp(1:N) = listGetReal( Material,'EPL Compressibility', N, Element % NodeIndexes, Found,&
              UnFoundFatal=UnFoundFatal)
        ! Previous default value: EPLComp(1:N) = 1.0D-2

        Porosity(1:N) = listGetReal( Material,'EPL Porosity', N, Element % NodeIndexes, Found,&
              UnFoundFatal=UnFoundFatal)
        ! Previous default value: Porosity(1:N) = 0.4D00

        BodyForce => GetBodyForce()
        IF ( ASSOCIATED( BodyForce ) ) THEN
           bf_id = GetBodyForceId()
           g = 0.0_dp  
           g(1,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 1',Found)
           g(2,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 2',Found)
           g(3,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 3',Found)
           Gravity(1:N) = SQRT(SUM(g**2.0/N))
        END IF

        EPLThick(1:N) = listGetReal( Material,'EPL Thickness', N, Element % NodeIndexes, Found,&
              UnFoundFatal=UnFoundFatal)
        ! Previous default value: EPLThick(1:N) = 1.0D00

        Density(1:N) = ListGetReal( Material, 'Water Density',  N, Element % NodeIndexes, Found,&
              UnFoundFatal=UnFoundFatal)
        ! Previous default value: Density(1:N) = 1.0055D-18

        !-----------------------------------------------------------------------------
        ! Get Upper limit:
        !-----------------------------------------------------------------------------
        IF (ApplyDirichlet)THEN
           UpperLimit(Element % Nodeindexes(1:n)) = ListGetReal(Material,TRIM(VariableName) // & 
                ' Upper Limit',n,Element % NodeIndexes, Found)
           
           IF (.NOT. Found) THEN
              WRITE(Message,'(a,i10)') 'No upper limit of solution for element no. ', t
              CALL INFO(SolverName, Message, level=10)
           END IF
        END IF
        !------------------------------------------------------------------------------
        ! add water contribution from transfer from EPL to sediment 
        !------------------------------------------------------------------------------
        LOAD = 0.0D00

        IF (ASSOCIATED(BodyForce)) THEN
           bf_id = GetBodyForceId()
           EPLToIDS(1:N) = ListGetReal( BodyForce, 'EPLToIDS Transfer', &
                N, Element % NodeIndexes(1:N),Found,UnFoundFatal=UnFoundFatal)
        ! Previous default value: EPLToIDS(1:N) = 0.0 
        END IF

        LOAD(1:N) = - EPLToIDS(1:N)

        !Computing the Storing coeficient and transmitivity of the equivqlent porous layer 
        !----------------------------------------------------------------------------------
        CALL ListGetRealArray( Material,'EPL Transmitivity',Hwrk,N, Element % NodeIndexes,Found)
        IF(.NOT.Found) THEN
           WRITE (Message,'(A)') 'A EPL Transmitivity is needed'
           CALL FATAL(SolverName,Message)
        ELSE
           Transmitivity = 0.0D0
           IF ( SIZE(Hwrk,1) == 1 ) THEN
              DO i=1,3
                 Transmitivity( i,i,1:N ) = Hwrk( 1,1,1:N)
              END DO
           ELSE
              WRITE(Message,'(a,a,a)') 'Keyword >EPL Transmitivity< should be isotrop '
              CALL INFO(SolverName,Message,Level=4)              
           END IF
        END IF

        StoringCoef(1:N) = EPLThick(1:N) * &
             Gravity(1:N) * Porosity(1:N) * Density(1:N) * &
             (WatComp + EPLComp(1:N)/Porosity(1:N))

        !------------------------------------------------------------------------------
        ! dummy input array for faking   heat capacity, density, temperature, 
        !                                enthalpy and viscosity
        ! Also faking the velocities to compute the Water load
        !------------------------------------------------------------------------------
        Work = 1.0d00
        Zero = 0.0D00
        VNull = 0.0D00
        Nochange = 0.0D00

        !------------------------------------------------------------------------------
        ! no contribution proportional to Water load by default
        ! No convection by default give C1 equal 0
        ! Viscosity is a dummy
        !------------------------------------------------------------------------------
        C0 = 0.0d00
        C1 = 0.0d00
        Viscosity = 0.0D00

        !------------------------------------------------------------------------------
        ! Get element local matrices, and RHS vectors
        !------------------------------------------------------------------------------
        MASS = 0.0d00
        STIFF = 0.0d00
        FORCE = 0.0D00

        ! cartesian coords
        !----------------
        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           CALL DiffuseConvectiveCompose( MASS, STIFF, FORCE, LOAD, &
                StoringCoef, C0, C1(1:N), Transmitivity, &
                .FALSE., Zero, Zero, VNull, VNull, VNull, &
                Nochange(1,1:N),Nochange(2,1:N),Nochange(3,1:N),&
                Viscosity, Density, Pressure, Zero, Zero,&
                .FALSE., Stabilize, UseBubbles, Element, N, ElementNodes )

           ! special coords (account for metric)
           !-----------------------------------
        ELSE
           CALL DiffuseConvectiveGenCompose( &
                MASS, STIFF, FORCE, LOAD, &
                StoringCoef, C0, C1(1:N), Transmitivity, &
                .FALSE., Zero, Zero, VNull, VNull, VNull, &
                Nochange(1,1:N),Nochange(2,1:N),Nochange(3,1:N), Viscosity,&
                Density, Pressure, Zero, Zero,.FALSE.,&
                Stabilize, Element, N, ElementNodes )
        END IF
        
        !------------------------------------------------------------------------------
        ! If time dependent simulation add mass matrix to stiff matrix
        !------------------------------------------------------------------------------
        TimeForce = FORCE
        IF ( TransientSimulation ) THEN
           IF ( UseBubbles ) FORCE = 0.0d0
           CALL Default1stOrderTime( MASS,STIFF,FORCE )
        END IF

        !------------------------------------------------------------------------------
        !  Update global matrices from local matrices
        !------------------------------------------------------------------------------
        IF (  UseBubbles ) THEN
           CALL Condensate( N, STIFF, FORCE, TimeForce )
           IF (TransientSimulation) CALL DefaultUpdateForce( TimeForce )
        END IF
        CALL DefaultUpdateEquations( STIFF, FORCE )
        
     END DO     !  Bulk elements
     CALL DefaultFinishBulkAssembly()

     !------------------------------------------------------------------------------
     ! Neumann & Newton boundary conditions
     !------------------------------------------------------------------------------
     DO t=1, Solver % Mesh % NumberOfBoundaryElements

        ! get element information
        Element => GetBoundaryElement(t)

        IF ( .NOT.ActiveBoundaryElement() ) CYCLE
        n = GetElementNOFNodes()

        BC => GetBC()
        bc_id = GetBCId( Element )
        CALL GetElementNodes( ElementNodes,Element )            

        IF ( ASSOCIATED( BC ) ) THEN
           
           !Initialize and check that we are on the correct boundary part!
           STIFF=0.0D00
           FORCE=0.0D00
           MASS=0.0D00
           LOAD=0.0D00
           TransferCoeff = 0.0D00
           EPLHeadExt = 0.0D00
           FluxBC = .FALSE.
           FluxBC =  GetLogical(BC,TRIM(Solver % Variable % Name) // ' Flux BC', Found)
           
           IF (FluxBC) THEN

              !BC: -k@Hw/@n = a(Hw - HwExt)
              !Check it if you want to use it this would represent a flux through a leaking media
              !Checking of equations, parameters or units have not been done
              !----------------------------------------------------------------------------------
              TransferCoeff(1:n) = GetReal( BC, TRIM(Solver % Variable % Name) //  ' Transfer Coefficient', Found )
              IF ( ANY(TransferCoeff(1:n).NE.0.0d0) ) THEN
                 EPLHeadExt(1:n) = GetReal( BC, TRIM(Solver % Variable % Name) // ' External Value', Found )   
                 DO j=1,n
                    LOAD(j) = LOAD(j) +  TransferCoeff(j) * EPLHeadExt(j)
                 END DO
              END IF
              
              !---------------
              !BC: -k@T/@n = q
              !---------------
              LOAD(1:n)  = LOAD(1:n) + &
                   GetReal( BC, TRIM(Solver % Variable % Name) // ' Water Flux', Found )
              
              !-------------------------------------
              ! set boundary due to coordinate system
              ! -------------------------------------
              IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                 CALL DiffuseConvectiveBoundary( STIFF,FORCE, &
                      LOAD,TransferCoeff,DummyLogical,DummyRealArray,&
                      DummyRealArray,Element,n,ElementNodes )
              ELSE
                 CALL DiffuseConvectiveGenBoundary(STIFF,FORCE,&
                      LOAD,TransferCoeff,Element,n,ElementNodes ) 
              END IF
           END IF
        END IF
        
        !------------------------------------------------------------------------------
        ! Update global matrices from local matrices
        !------------------------------------------------------------------------------
        IF ( TransientSimulation ) THEN
           MASS = 0.d0
           CALL Default1stOrderTime( MASS, STIFF, FORCE )
        END IF
        CALL DefaultUpdateEquations( STIFF, FORCE )
     END DO ! Neumann & Newton BCs

     CALL DefaultFinishAssembly()
     CALL DefaultDirichletBCs()
     
     OldValues = SystemMatrix % Values
     OldRHS = ForceVector 

     !------------------------------------------------------------------------------
     ! Dirichlet method - matrix and force-vector manipulation 
     !------------------------------------------------------------------------------
     IF (ApplyDirichlet) THEN

        !dealing with EPL spreading
        !------------------------------
        DO t=1,Solver % NumberOfActiveElements

           Element => GetActiveElement(t)
           N = GetElementNOFNodes(Element)
           CALL GetElementNodes( ElementNodes,Element )

           !Spreading is needed if  one node of the element have a water head equal 
           ! or above the upper limit and if one node of the element is still closed
           !------------------------------------------------------------------------
           IF ((ANY(ActiveEPL(Element % NodeIndexes(1:N))))&
                .AND.(ANY(OpenEPL(PipingPerm(Element % NodeIndexes(1:N))).LT.0.0)))THEN
              
              !Get the inefficient layer water head
              !------------------------------------
              CALL GetScalarLocalSolution(IDSHead,IDSHeadName,Element)
              DO i=1,N
                 !The node to be activated is the one with the lower inefficient 
                 !water head
                 !---------------------------------------------------------------
                 IF((IDSHead(i).EQ.MINVAL(IDSHead(1:N))).AND.&
                      (OpenEPL(PipingPerm(Element % NodeIndexes(i))).GT.0.0))&
                      OpenEPL(PipingPerm(Element % NodeIndexes(i))) = -1.0
              END DO
           END IF
        END DO
     END IF

     CALL Info( TRIM(SolverName) // ' EPL Solver', 'Assembly done', Level=4 ) 

     !------------------------------------------------------------------------------
     !     Solve the system
     !------------------------------------------------------------------------------
     at = CPUTime() - at
     st = CPUTime()

     PrevNorm = Solver % Variable % Norm
     Norm = DefaultSolve()

     SystemMatrix % Values = OldValues
     ForceVector = OldRHS

     !-------------------------------------------------------------------------------
     ! determine "active" nodes set, the node for which the EPL Water head is above 
     ! the upper limit
     !-------------------------------------------------------------------------------
     DO i=1,Model % Mesh % NumberOfNodes 
        k = HomolPerm(i)
        l = EPLHeadPerm(i)

        IF ( (k.LE.0).OR.(l.LE.0)) CYCLE
        EPLHeadHomologous(k) = EPLHead(l) - UpperLimit(i)
        
        !----------------------------------------------------------
        ! if upper limit is exceeded, flag the node
        !----------------------------------------------------------
        IF ((ApplyDirichlet).AND.(EPLHeadHomologous(k).GE.0.0 )) THEN
           ActiveEPL(i) = .TRUE.           
           EPLHeadHomologous(k) = LinearTol
        END IF
     END DO

     !------------------------------------------
     ! special treatment for periodic boundaries
     !------------------------------------------
     DO t=1, Solver % Mesh % NumberOfBoundaryElements

        ! get element information
        Element => GetBoundaryElement(t)
        IF ( .NOT.ActiveBoundaryElement() ) CYCLE
        n = GetElementNOFNodes()
        BC => GetBC()
        bc_id = GetBCId( Element )

        CALL GetElementNodes( ElementNodes,Element )

        IF ( ASSOCIATED( BC ) ) THEN    
           IsPeriodicBC = GetLogical(BC,'Periodic BC ' // TRIM(Solver % Variable % Name), Found)
           IF (.NOT.Found) IsPeriodicBC = .FALSE.
           IF (IsPeriodicBC) THEN 
              DO i=1,N
                 IF  (ActiveEPL(Element % NodeIndexes(i))) THEN
                    ActiveEPL(Element % NodeIndexes(i)) = .FALSE.
                 END IF
              END DO
           END IF
        END IF
     END DO
     Norm = Solver % Variable % Norm

     !Checking for convergence
     !------------------------
     st = CPUTime()-st
     totat = totat + at
     totst = totst + st
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',Peniter,' Assembly tot: (s)', at, totat
     CALL Info( SolverName, Message, Level=4 )
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',Peniter,' Solve:    (s)', st, totst
     CALL Info( SolverName, Message, Level=4 )

     IF ((PrevNorm + Norm).NE.0.0d0 ) THEN
        RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
     ELSE
        RelativeChange = 0.0d0
     END IF

     WRITE( Message, * ) 'Result Norm   : ',Norm
     CALL Info( SolverName, Message, Level=4 )
     WRITE( Message, * ) 'Relative Change : ',RelativeChange
     CALL Info( SolverName, Message, Level=4 )
     
     IF (RelativeChange.LT.NonlinearTol) THEN
        EXIT
     END IF

     IF((PenIter.EQ.NonlinearIter).AND.(RelativeChange.GT.NonlinearTol))THEN
        Write(*,*)'NOT CONVERGED'
        STOP
     END IF
  END DO
  
END SUBROUTINE EPLSolver
