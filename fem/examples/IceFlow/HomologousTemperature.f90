!/*****************************************************************************
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
! * Module containing a solver for a genral advection diffusion problem
! *
! ******************************************************************************
! *
! *                     Author:       Juha Ruokolainen
! *
! *                    Address: CSC - IT Center for Science Ltd.
! *                                Keilaranta 14, P.O. BOX 405
! *                                  02101 Espoo, Finland
! *                                  Tel. +358 0 457 2723
! *                                Telefax: +358 0 457 2302
! *                              EMail: Juha.Ruokolainen@csc.fi
! *
! *                       Date: 08 Jun 1997
! *
! *                Modified by: Thomas Zwinger
! *
! *       Date of modification: 28 Nov 2005
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE HomologousTempSolver( Model,Solver,Timestep,TransientSimulation )
!DLLEXPORT HomologousTempSolver
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the convection diffusion equation with limiters!
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
!******************************************************************************
     USE DiffuseConvective
     USE DiffuseConvectiveGeneral
     USE Differentials
     USE MaterialModels

!     USE Adaptive
     USE DefUtils

!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------

     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver

     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Solver_t), POINTER :: PointerToSolver
     TYPE(Matrix_t), POINTER :: StiffMatrix
     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: Element,ParentElement
     TYPE(Variable_t), POINTER :: TempSol,FlowSol,CurrentSol, MeshSol
     TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants

     INTEGER :: i,j,k,l,m,n,t,tt,iter,k1,k2,body_id,eq_id,material_id, &
          istat, LocalNodes,bf_id, bc_id, other_body_id, DIM, &
          NSDOFs,NewtonIter,NonlinearIter, NOFControlHeaters, &
          NumberOfBoundaryNodes
     INTEGER, POINTER :: NodeIndexes(:), TempPerm(:),FlowPerm(:),CurrentPerm(:),MeshPerm(:)
     INTEGER, ALLOCATABLE :: BoundaryMap(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: ConvectionFlag, VariableName, SolverName

     LOGICAL :: Stabilize = .TRUE., Bubbles = .TRUE., UseBubbles, &
          NewtonLinearization = .FALSE., Found, FluxBC, Permeable=.TRUE., &
          AllocationsDone = .FALSE.,  SubroutineVisited = .FALSE., FirstTime=.TRUE.,&
          CorrectNegative, LimitSolution
     LOGICAL, ALLOCATABLE :: IsBoundaryNode(:), LimitedSolution(:,:)

     REAL(KIND=dp) :: NonlinearTol,NewtonTol,NonLinearTolMin,NonLinearTolDegrease, &
           OldNonlinearTol, Relax, RelaxIncrease, RelaxMax, &
          SaveRelax,dt,CumulativeTime, RelativeChange, &
          Norm,PrevNorm,S,C, &
          ReferencePressure=0.0d0, UzawaParameter, &
          HeatCapacityGradient(3)
     

     REAL(KIND=dp), POINTER :: Temp(:),FlowSolution(:), &
          ForceVector(:), &
          PrevSolution(:), HC(:), Hwrk(:,:,:), &
          SlipVelocityX(:), SlipVelocityY(:), SlipVelocityZ(:), PointerOnSlipVelocity(:), &
          PointerOnLagrange(:)
     
     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:), HeatConductivity(:,:,:), &
         FORCE(:), Pressure(:),  MeshVelocity(:,:),&
         TempVeloU(:),TempVeloV(:),TempVeloW(:),TimeForce(:), &
         TransferCoeff(:), LocalTemp(:), Work(:), C1(:), C0(:), Zero(:), Viscosity(:),&
         UpperLimit(:),LowerLimit(:),LagrangeMultiplier1(:),LagrangeMultiplier2(:), &
         HeatCapacity(:), PrevHeatCapacity(:), Density(:), TempExt(:)
     REAL(KIND=dp), ALLOCATABLE, TARGET :: VarSlipVelocity(:)
     REAL(KIND=dp) :: at,at0,totat,st,totst,t1,CPUTime,RealTime

     SAVE MeshVelocity,  &
          TempVeloU, TempVeloV, TempVeloW,&
          MASS, STIFF, LOAD, FORCE,&
          ElementNodes, Hwrk, HeatConductivity, TransferCoeff, &
          AllocationsDone, FirstTime, TimeForce, &
          LocalNodes, LocalTemp, Work, Zero, Viscosity, Dim, Norm, &
          C0, C1, Pressure, CorrectNegative, VariableName, SolverName, &
          UpperLimit, LowerLimit,LagrangeMultiplier1, LagrangeMultiplier2, &
          IsBoundaryNode, BoundaryMap, NumberOfBoundaryNodes, &
          LimitSolution, LimitedSolution, M, HeatCapacity, PrevHeatCapacity, Density, &
          TempExt, OldNonlinearTol, NonLinearTolMin,NonLinearTolDegrease

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
     StiffMatrix => Solver % Matrix
     ForceVector => Solver % Matrix % RHS

     PointerToSolver => Solver

     TempSol => Solver % Variable
     TempPerm  => TempSol % Perm
     Temp => TempSol % Values
     
     LocalNodes = COUNT( TempPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN
     
     FlowSol => VariableGet( Solver % Mesh % Variables, 'Flow Solution' )
     IF ( ASSOCIATED( FlowSol ) ) THEN
        FlowPerm     => FlowSol % Perm
        NSDOFs       =  FlowSol % DOFs
        FlowSolution => FlowSol % Values
     END IF

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
        N = Solver % Mesh % MaxElementNodes
        M = Model % Mesh % NumberOfNodes
        DIM = CoordinateSystemDimension()
        SolverName = 'HomologousTemp ('// TRIM(Solver % Variable % Name) // ')'
        VariableName = TRIM(Solver % Variable % Name)

        IF ( AllocationsDone ) THEN
           DEALLOCATE(                 &
                MeshVelocity,          &
                TempVeloU,             &
                TempVeloV,             &
                TempVeloW,             &
                Pressure,              &
                ElementNodes % x,      &
                ElementNodes % y,      &
                ElementNodes % z,      &
                Work,Zero,             &
                Viscosity,             &
                HeatCapacity,          &
                PrevHeatCapacity,      &
                Density,               &
                TempExt,                  &
                C1,                    &
                C0,                    &
                TransferCoeff,         &
                LocalTemp,             &
                HeatConductivity,      &
                MASS,                  &
                STIFF,LOAD,            &
                FORCE,                 &
                TimeForce,             &
                LagrangeMultiplier1,   &
                LagrangeMultiplier2,   &
                LowerLimit,            &
                UpperLimit,            &
                IsBoundaryNode,        &
                LimitedSolution,       &
                BoundaryMap)
        END IF
        
        ALLOCATE(                                  &
             MeshVelocity( 3,N ),                  &
             TempVeloU( N ),                       &
             TempVeloV( N ),                       &
             TempVeloW( N ),                       &
             Pressure( N ),                        &
             ElementNodes % x( N ),                &
             ElementNodes % y( N ),                &
             ElementNodes % z( N ),                &
             Work( N ), Zero( N ),                 &
             Viscosity( N ),                       &
             HeatCapacity( M ),                    &
             PrevHeatCapacity( M ),                &
             Density( N ),                         &
             TempExt( N ),                         &
             C1(N),                                &
             C0(N),                                &
             TransferCoeff( N ),                   &
             LocalTemp( N ),                       &
             HeatConductivity( 3,3,N ),            &
             MASS(  2*N,2*N ),                     &
             STIFF( 2*N,2*N ),LOAD( N ),           &
             FORCE( 2*N ), TimeForce(2*N),         &
             LagrangeMultiplier1(M),               &
             LagrangeMultiplier2(M),               &
             LowerLimit(M),                        &
             UpperLimit(M),                        &
             IsBoundaryNode(M),                    &
             LimitedSolution(2,M),                 &
             BoundaryMap(M),                       &
             STAT=istat )


        IF ( istat /= 0 ) THEN
           CALL Fatal( SolverName, 'Memory allocation error' )
        ELSE
           CALL INFO(SolverName, 'Memory allocation done', level=1 )
        END IF
        
        
        LagrangeMultiplier1 = 0.0d0
        LagrangeMultiplier2 = 0.0d0

       
        ! Assign Variable for Lagrangemultiplier1
        !       PointerOnLagrange => LagrangeMultiplier1(1:LocalNodes)
!        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, PointerToSolver, &
!             TRIM(Solver % Variable % Name) // ' Lower Limit Lagrange', 1, LagrangeMultiplier1(1:LocalNodes), TempPerm)
!       PointerOnLagrange => LagrangeMultiplier2(1:LocalNodes)
!        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, PointerToSolver, &
!             TRIM(Solver % Variable % Name) // ' Upper Limit Lagrange', 1, LagrangeMultiplier2(1:LocalNodes), TempPerm)

        !-----------------------------------
        
        ! Get Information on Boundary Nodes
        !----------------------------------
        IsBoundaryNode = .FALSE.
        BoundaryMap = 0
        
        DO t=1, Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(t,Solver)
           IF (ActiveBoundaryElement(Element)) THEN
              n = GetElementNOFNodes(Element)
              IsBoundaryNode(Element % NodeIndexes(1:N)) = .TRUE.
           ELSE
              IsBoundaryNode = .FALSE.
           END IF
        END DO
       
        ! Get Bulk to Boundary Mapping
        !-----------------------------
        NumberOfBoundaryNodes = 0
        DO i=1,M
           IF (IsBoundaryNode(i)) THEN
              NumberOfBoundaryNodes = NumberOfBoundaryNodes + 1
              BoundaryMap(i) = NumberOfBoundaryNodes
           END IF
        END DO
        
        NULLIFY( Hwrk )
        
        AllocationsDone = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Say hello
!------------------------------------------------------------------------------
     WRITE(Message,'(A,A)')&
          'Homologous Temperature Solver for variable', VariableName
     CALL INFO(SolverName,Message,Level=1)

!------------------------------------------------------------------------------
!    Read physical and numerical constants and initialize 
!------------------------------------------------------------------------------
     Constants => GetConstants()
     SolverParams => GetSolverParams()

     Stabilize = GetLogical( SolverParams,'Stabilize',Found )
     IF (.NOT. Found) Stabilize = .FALSE.
     UseBubbles = GetLogical( SolverParams,'Bubbles',Found )
     IF ( .NOT.Found .AND. (.NOT.Stabilize)) UseBubbles = .TRUE.

     NonlinearIter = GetInteger(   SolverParams, &
                     'Nonlinear System Max Iterations', Found )
     IF ( .NOT.Found ) NonlinearIter = 1

     IF (FirstTime) THEN
        NonlinearTol  = GetConstReal( SolverParams, &
                     'Nonlinear System Convergence Tolerance',    Found )

        OldNonlinearTol = NonlinearTol  
        NonlinearTolDegrease = GetConstReal( SolverParams, &
             'Nonlinear System Convergence Tolerance Degrease',    Found )
        IF (.NOT.Found) NonlinearTolDegrease = 1.0D00

        NonlinearTolMin  = GetConstReal( SolverParams, &
             'Nonlinear System Convergence Tolerance Min',    Found )
     
        IF (.NOT.Found) NonlinearTolDegrease = 1.0D00
        
     END  IF

     NewtonTol     = GetConstReal( SolverParams, &
                      'Nonlinear System Newton After Tolerance',  Found )
     NewtonIter    = GetInteger(   SolverParams, &
                      'Nonlinear System Newton After Iterations', Found )

     Relax = GetConstReal( SolverParams, &
               'Nonlinear System Relaxation Factor',Found )

     IF ( .NOT.Found ) Relax = 1.0D00

     RelaxIncrease = GetConstReal( SolverParams, &
          'Nonlinear System Relaxation Factor Increase',Found )
     IF (.NOT.Found) RelaxIncrease = 1.0D00

     RelaxMax = GetConstReal( SolverParams, &
          'Nonlinear System Relaxation Factor Max',Found )
     IF (.NOT.Found) RelaxMax = Relax     

     UzawaParameter = GetConstReal( SolverParams, &
               'Uzawa Parameter', Found)
     IF (.NOT.Found) THEN                
        LimitSolution = .FALSE.
        CALL INFO(SolverName, 'No Uzawa Parameter found. No limitting of solution', level=1)
     ELSE
        WRITE(Message,'(a,e13.4)') 'Uzawa Parameter:',  UzawaParameter
        CALL INFO(SolverName, Message, level=1)
        LimitSolution = .TRUE.
     END IF

     SaveRelax = Relax
     dt = Timestep
     CumulativeTime = 0.0d0


!------------------------------------------------------------------------------
!   time stepping loop.
!------------------------------------------------------------------------------
     DO WHILE( CumulativeTime < Timestep-1.0d-12 .OR. .NOT. TransientSimulation )

!------------------------------------------------------------------------------
!       The first time around this has been done by the caller...
!------------------------------------------------------------------------------
        IF ( TransientSimulation .AND. .NOT.FirstTime ) THEN
           CALL InitializeTimestep( Solver )
        END IF

!------------------------------------------------------------------------------
!       Save current solution
!------------------------------------------------------------------------------
        ALLOCATE( PrevSolution(LocalNodes) )
        PrevSolution = Temp(1:LocalNodes)

        totat = 0.0d0
        totst = 0.0d0
        IF (.NOT.FirstTime) THEN
           WRITE( Message, * ) 'Nonlinear Tolerance: old:',OldNonLinearTol,&
                ' new:', MAX(NonLinearTolMin,OldNonLinearTol*NonLinearTolDegrease)
           CALL Info( SolverName, Message, Level=4 )
           NonLinearTol = MAX(NonLinearTolMin,OldNonLinearTol*NonLinearTolDegrease)
           OldNonLinearTol = NonLinearTol
        ELSE
           WRITE( Message, * ) 'Nonlinear Tolerance at start:',NonLinearTol
           CALL Info( SolverName, Message, Level=4 )
           OldNonLinearTol = NonLinearTol
        END IF
!------------------------------------------------------------------------------
!       non-linear system iteration loop
!------------------------------------------------------------------------------
        DO iter=1,NonlinearIter

           IF (.NOT.FirstTime) THEN
              WRITE( Message, * ) 'Relaxation Factor: old:',Relax,&
                   ' new:', MIN(RelaxMax,Relax * Relaxincrease) 
              CALL Info( SolverName, Message, Level=4 )
              Relax = MIN(RelaxMax,Relax * Relaxincrease)
           ELSE
              WRITE( Message, * ) 'Relaxation Factor:',Relax
              CALL Info( SolverName, Message, Level=4 )
           END IF

           FirstTime = .FALSE.


           at  = CPUTime()
           at0 = RealTime()

           CALL Info( SolverName, ' ', Level=4 )
           CALL Info( SolverName, ' ', Level=4 )
           CALL Info( SolverName, '-------------------------------------',Level=4 )
           WRITE( Message,'(A,A,I3,A,I3)') &
                TRIM(Solver % Variable % Name),  ' iteration no.', iter,' of ',NonlinearIter
           CALL Info( SolverName, Message, Level=4 )
           CALL Info( SolverName, '-------------------------------------',Level=4 )
           CALL Info( SolverName, ' ', Level=4 )
           CALL Info( SolverName, 'Starting Assembly...', Level=4 )

!------------------------------------------------------------------------------
           CALL DefaultInitialize()

!-----------------------------------------------------------------------------
!      Update Lagrange multiplier:
!-----------------------------------------------------------------------------
           IF (LimitSolution) THEN
              LowerLimit = 0.0d00
              UpperLimit = 0.0d00

              LimitedSolution = .TRUE.

              ! Get lower and Upper limit

              DO t=1,Solver % NumberOfActiveElements
                 Element => GetActiveElement(t)
                 n = GetElementNOFNodes()
                 CALL GetElementNodes( ElementNodes )
                 Material => GetMaterial()
                 ! lower limit
                 !------------
                 LowerLimit(Element % Nodeindexes(1:N)) = ListGetReal(Material,TRIM(Solver % Variable % Name) // & 
                      ' Lower Limit',n,Element % NodeIndexes, Found)
                 IF (.NOT.Found) THEN
                    LimitedSolution(1,Element % Nodeindexes(1:N)) = .FALSE.
                    WRITE(Message,'(a,i)') 'No lower limit of solution for element no. ', t
                    CALL INFO(SolverName, Message, level=9)
                 END IF                 
                 ! upper limit
                 !------------
                 UpperLimit(Element % Nodeindexes(1:N)) = ListGetReal(Material,TRIM(Solver % Variable % Name) // & 
                      ' Upper Limit',n,Element % NodeIndexes, Found)
                 IF (.NOT. Found) THEN
                    LimitedSolution(2,Element % Nodeindexes(1:N)) = .FALSE.
                    WRITE(Message,'(a,i)') 'No upper limit of solution for element no. ', t
                    CALL INFO(SolverName, Message, level=9)
                 END IF
              END DO


              ! get Lagrange Multipliers
              ! ------------------------
              DO i = 1,M
                 ! upper limit
                 !------------
                 IF (LimitedSolution(2,i)) THEN
                    LagrangeMultiplier2(i) = MAX( 0.0d0, LagrangeMultiplier2(i) + UzawaParameter * &
                         (Temp(TempPerm(i)) - UpperLimit(i)))
!                    LagrangeMultiplier2(i) = MAX( 0.0d0,  UzawaParameter * (Temp(TempPerm(i)) - UpperLimit(i)))
                 ELSE 
                    LagrangeMultiplier2(i) = 0.0d00
                 END IF
                 ! lower limit
                 !-----------
                 IF (LimitedSolution(1,i)) THEN
                    LagrangeMultiplier1(i) = MIN( 0.0d0, LagrangeMultiplier1(i) + UzawaParameter * &
                         ( Temp(TempPerm(i)) - LowerLimit(i)))
!                    LagrangeMultiplier2(i) = MIN( 0.0d0,  UzawaParameter * (Temp(TempPerm(i)) - LowerLimit(i)))
                 ELSE
                    LagrangeMultiplier1(i) = 0.0d00
                 END IF
              END DO
           END IF

           CALL INFO(SolverName,'Max/min values Lagrange multiplier ...',Level=4)
           WRITE(Message,'(a,e13.6,a,e13.6)') &
                '... for lower constraint:', MAXVAL( LagrangeMultiplier1(:) ),'/',MINVAL( LagrangeMultiplier1(:) )
           CALL INFO(SolverName,Message,Level=4)
           WRITE(Message,'(a,e13.6,a,e13.6)') &
                '... or upper constraint:', MAXVAL( LagrangeMultiplier2(:) ),'/',MINVAL( LagrangeMultiplier2(:) )
           CALL INFO(SolverName,Message,Level=4)
           WRITE(Message,'(a,e13.6,a,e13.6)') &
                'Max/min values Temperature:', MAXVAL( Temp(:)),'/',MINVAL( Temp(:))
           CALL INFO(SolverName,Message,Level=4)
!------------------------------------------------------------------------------
           body_id = -1
           NULLIFY(Material)

!           PRINT *, 'Lagr.:', Norm, Solver % Variable % Norm 
!------------------------------------------------------------------------------
!      Bulk elements
!------------------------------------------------------------------------------
           DO t=1,Solver % NumberOfActiveElements

              IF ( RealTime() - at0 > 1.0 ) THEN
                 WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
                      (Solver % NumberOfActiveElements-t) / &
                      (1.0*Solver % NumberOfActiveElements)), ' % done'

                 CALL Info( SolverName, Message, Level=5 )

                 at0 = RealTime()
              END IF
!------------------------------------------------------------------------------
!             Check if this element belongs to a body where scalar equation
!             should be calculated
!------------------------------------------------------------------------------
              Element => GetActiveElement(t,Solver)
              IF (.NOT.ActiveBoundaryElement(Element,Solver)) CYCLE
!------------------------------------------------------------------------------
              IF ( Element % BodyId /= body_id ) THEN
!------------------------------------------------------------------------------
                 Equation => GetEquation()
                 IF (.NOT.ASSOCIATED(Equation)) THEN
                    WRITE (Message,'(A,I3)') 'No Equation  found for boundary element no. ', t
                    CALL FATAL(SolverName,Message)
                 END IF

                 ConvectionFlag = GetString( Equation, 'Convection', Found )

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

!------------------------------------------------------------------------------            
              END IF
!------------------------------------------------------------------------------
              
              N = GetElementNOFNodes(Element)
              CALL GetElementNodes( ElementNodes )

!              CALL GetTempLocalSolution( LocalTemp )
!------------------------------------------------------------------------------
!             Get element material parameters
!------------------------------------------------------------------------------
              CALL ListGetRealArray( Material,TRIM(Solver % Variable % Name) // &
                   ' Heat Conductivity',Hwrk,n, Element % NodeIndexes )

              HeatConductivity = 0.0d0
              IF ( SIZE(Hwrk,1) == 1 ) THEN
                 DO i=1,3
                    HeatConductivity( i,i,1:N ) = Hwrk( 1,1,1:N)
                 END DO
              ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
                 DO i=1,MIN(3,SIZE(Hwrk,1))
                    HeatConductivity(i,i,1:N) = Hwrk(i,1,1:N)
                 END DO
              ELSE
                 DO i=1,MIN(3,SIZE(Hwrk,1))
                    DO j=1,MIN(3,SIZE(Hwrk,2))
                       HeatConductivity(i,j,1:N) = Hwrk(i,j,1:N)
                    END DO
                 END DO
              END IF

              HeatCapacity(1:N) =  ListGetReal( Material,  TRIM(Solver % Variable % Name) // &
                   ' Heat Capacity', n, Element % NodeIndexes, Found )
              IF (.NOT.Found) THEN
                 HeatCapacity = 0.0D00
                 WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', TRIM(Solver % Variable % Name) // &
                   ' Heat Capacity', '< not found for element ', t, ' material ', material_id
                 CALL INFO(SolverName,Message,Level=4)
              END IF

               Density(1:N) = ListGetReal( Material, 'Density',  N, Element % NodeIndexes, Found )
              IF (.NOT.Found) THEN
                 Density = 0.0D00
                 WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Density< not found for element ',&
                      t, ' material ', material_id
                 CALL INFO(SolverName,Message,Level=4)
              END IF

              Viscosity(1:N) = ListGetReal( Material, 'Viscosity',  N, Element % NodeIndexes, Found )
              IF (.NOT.Found) THEN
                 Viscosity = 0.0D00
                 WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Viscosity< not found for element ',&
                      t, ' material ', material_id
                 CALL INFO(SolverName,Message,Level=4)
              END IF

!              PRINT *, 'kappa=', HeatConductivity( 1:MIN(3,SIZE(Hwrk,1)),1:MIN(3,SIZE(Hwrk,2)),1:N ) 
!              PRINT *, 'cp=', HeatCapacity(1:N)

!------------------------------------------------------------------------------
!        Get mesh velocity
!------------------------------------------------------------------------------
              MeshVelocity = 0.0d0
              CALL GetVectorLocalSolution( MeshVelocity, 'Mesh Velocity')
!------------------------------------------------------------------------------
!        Get scalar velocity
!------------------------------------------------------------------------------         
              TempVeloU = 0.0d00
              TempVeloV = 0.0d00
              TempVeloW = 0.0d00
              ! asuming convection or ALE mesh contribution by default
              DO i=1,N
                 C1(i) = Density(i) * HeatCapacity(i)
!                 C1(i) = HeatCapacity(i)
!                 PRINT *, 'C1(',i,')=',C1(i)
              END DO

              IF ( ConvectionFlag == 'constant' ) THEN
                 TempVeloU(1:N)= GetReal( Material, 'Convection Velocity 1', Found )
                 TempVeloV(1:N) = GetReal( Material, 'Convection Velocity 2', Found )
                 TempVeloW(1:N) = GetReal( Material, 'Convection Velocity 3', Found )                 
              ELSE IF ( ConvectionFlag == 'computed' ) THEN
                 DO i=1,n
                    k = FlowPerm(Element % NodeIndexes(i))
                    IF ( k > 0 ) THEN
                       Pressure(i) = FlowSolution(NSDOFs*k) + ReferencePressure
!------------------------------------------------------------------------------
                       SELECT CASE( NSDOFs )
                       CASE(3)
                          TempVeloU(i) = FlowSolution( NSDOFs*k-2 )
                          TempVeloV(i) = FlowSolution( NSDOFs*k-1 )
                          TempVeloW(i) = 0.0D0
                       CASE(4)
                          TempVeloU(i) = FlowSolution( NSDOFs*k-3 )
                          TempVeloV(i) = FlowSolution( NSDOFs*k-2 )
                          TempVeloW(i) = FlowSolution( NSDOFs*k-1 )
                       END SELECT
                    ELSE
                       TempVeloU(i) = 0.0D0
                       TempVeloV(i) = 0.0D0
                       TempVeloW(i) = 0.0D0
                    END IF
                 END DO
                 WRITE(Message,'(a,i5, a, i5)') 'Convection in element ', t, &
                      ' material ',  material_id
              ELSE  ! Conduction and ALE contribution only
                 IF (ANY( MeshVelocity /= 0.0d0 )) THEN
                    WRITE(Message,'(a,i5, a, i5)') 'Only mesh deformation in element ', t,&
                         ' material ',  material_id
                 ELSE ! neither convection nor ALE mesh deformation contribution -> all C1(1:N)=0
                    C1 = 0.0D0 
                    WRITE(Message,'(a,i5, a, i5)') 'No convection and mesh deformation in element ', t,&
                         ' material ',  material_id
                 END IF                 
              END IF
              CALL INFO(SolverName,Message,Level=14)
             
              ! no compressibility by default
              C0=0.0d00


!------------------------------------------------------------------------------
!        Add body forces, if any
!------------------------------------------------------------------------------
              LOAD = 0.0D00
              BodyForce => GetBodyForce()
              IF ( ASSOCIATED( BodyForce ) ) THEN
                 bf_id = GetBodyForceId()
                 !Given volume source
                 !-------------------
                 LOAD(1:N) = LOAD(1:N) +   &
                      GetReal( BodyForce, TRIM(Solver % Variable % Name) // ' Volume Source', Found )
!                 IF (.NOT.Found) LOAD = 0.0D00
              END IF
!------------------------------------------------------------------------------
!    Lagrange multipliers
!------------------------------------------------------------------------------

              LOAD(1:N) = LOAD(1:N) &
                   - LagrangeMultiplier1(Element % NodeIndexes(1:N)) &
                   - LagrangeMultiplier2(Element % NodeIndexes(1:N))

!              PRINT *, t, 'h= ', LOAD(1:N)

!if (iter > 1) print*,maxval(LagrangeMultiplier1(Element % NodeIndexes(1:N))), minval(LagrangeMultiplier1(Element % NodeIndexes(1:N)))

!------------------------------------------------------------------------------
!        Get element local matrices, and RHS vectors
!------------------------------------------------------------------------------
         ! dummy input array for faking   heat capacity, density, temperature, 
         !                                enthalpy and viscosity
              Work = 1.0d00
              Zero = 0.0D00


              MASS = 0.0d00
              STIFF = 0.0d00
              FORCE = 0.0D00
!------------------------------------------------------------------------------
              IF ( CurrentCoordinateSystem() == Cartesian ) THEN
!------------------------------------------------------------------------------
                 CALL DiffuseConvectiveCompose( &
                      MASS, STIFF, FORCE, LOAD, &
                      HeatCapacity, C0, C1(1:N), HeatConductivity, &
                      .FALSE., Zero, Zero, TempVeloU, TempVeloV, TempVeloW, &
                      MeshVelocity(1,1:N),MeshVelocity(2,1:N),MeshVelocity(3,1:N),&
                      Viscosity, Density, Pressure, &
                      .FALSE., Stabilize, UseBubbles, Element, n, ElementNodes )
!------------------------------------------------------------------------------
              ELSE
!------------------------------------------------------------------------------
                 CALL DiffuseConvectiveGenCompose( &
                      MASS, STIFF, FORCE, LOAD, &
                      HeatCapacity, C0, C1(1:N), HeatConductivity, &
                      .FALSE., Zero, Zero, TempVeloU, TempVeloV, TempVeloW, &
                      MeshVelocity(1,1:N),MeshVelocity(2,1:N),MeshVelocity(3,1:N), Viscosity,&
                      Density, Pressure,.FALSE.,&
                      Stabilize, Element, n, ElementNodes )
                 
!------------------------------------------------------------------------------
              END IF
!------------------------------------------------------------------------------
              Bubbles = UseBubbles  .AND. &
                   ( ConvectionFlag == 'computed' .OR. ConvectionFlag == 'constant' )         
!------------------------------------------------------------------------------
!           If time dependent simulation add mass matrix to stiff matrix
!------------------------------------------------------------------------------
              TimeForce  = FORCE
              IF ( TransientSimulation ) THEN
                 IF ( Bubbles ) FORCE = 0.0d0
                 CALL Default1stOrderTime( MASS,STIFF,FORCE )
              END IF
!------------------------------------------------------------------------------
!           Update global matrices from local matrices
!------------------------------------------------------------------------------
              IF (  Bubbles ) THEN
                 CALL Condensate( N, STIFF, FORCE, TimeForce )
                 IF (TransientSimulation) CALL DefaultUpdateForce( TimeForce )
              END IF

              CALL DefaultUpdateEquations( STIFF, FORCE )
!           END IF
              
!------------------------------------------------------------------------------
           END DO     !  Bulk elements
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
           DO t=1, Solver % Mesh % NumberOfBoundaryElements

              ! get element information
              Element => GetBoundaryElement(t)
              IF ( .NOT.ActiveBoundaryElement() ) CYCLE
              n = GetElementNOFNodes()
              IF ( GetElementFamily() == 1 ) CYCLE
              BC => GetBC()
              bc_id = GetBCId( Element )
              CALL GetElementNodes( ElementNodes )


              other_body_id = Element % BoundaryInfo % outbody
              IF (other_body_id < 1) THEN ! only one body in calculation
                 ParentElement => Element % BoundaryInfo % Right
                 IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => Element % BoundaryInfo % Left
              ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
                 ParentElement => Element % BoundaryInfo % Right
                 IF (.NOT.ASSOCIATED(ParentElement) ) THEN ! boundary is declared in input file, but does not xist
                    ParentElement => Element % BoundaryInfo % Left 
                 ELSE
                    IF (ParentElement % BodyId == other_body_id) ParentElement => Element % BoundaryInfo % Left
                 END IF
              END IF
     
              IF (.NOT.ASSOCIATED(ParentElement)) THEN 
                 WRITE (Message,'(A,I5,A)')&
                      'No parent element found on either side for boundary element no. ',t,'.'
                 CALL FATAL(SolverName,Message)
!!$              ELSE
!!$                 Material => GetMaterial(ParentElement)
              END IF


              IF ( ASSOCIATED( BC ) ) THEN            
                 ! Check that we are on the correct boundary part!
                 STIFF=0.0D00
                 FORCE=0.0D00
                 MASS=0.0D00
                 LOAD=0.0D00
                 TransferCoeff = 0.0D00
                 TempExt = 0.0D00
                 FluxBC = .FALSE.
                 FluxBC =  GetLogical(BC,TRIM(Solver % Variable % Name) // ' Flux BC', Found)

                 IF (FluxBC) THEN
                    !------------------------------
                    !BC: -k@T/@n = \alpha(T - TempExt)
                    !------------------------------
                    TransferCoeff(1:N) = GetReal( BC, TRIM(Solver % Variable % Name) //  ' Transfer Coefficient',Found )
                    IF ( ANY(TransferCoeff(1:N) /= 0.0d0) ) THEN
                       TempExt(1:N) = GetReal( BC, TRIM(Solver % Variable % Name) // ' External Value',Found )   
                       DO j=1,n
                          LOAD(j) = LOAD(j) +  TransferCoeff(j) * TempExt(j)
                       END DO
                    END IF

                    !---------------
                    !BC: -k@T/@n = q
                    !---------------
                    LOAD(1:N)  = LOAD(1:N) + &
                         GetReal( BC, TRIM(Solver % Variable % Name) // ' Heat Flux', Found )

                    ! -------------------------------------
                    ! set boundary due to coordinate system
                    ! -------------------------------------
                    IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                       CALL DiffuseConvectiveBoundary( STIFF,FORCE, &
                            LOAD,TransferCoeff,Element,n,ElementNodes )
                    ELSE
                       CALL DiffuseConvectiveGenBoundary(STIFF,FORCE,&
                            LOAD,TransferCoeff,Element,n,ElementNodes ) 
                    END IF
                 END IF
              END IF

!------------------------------------------------------------------------------
!         Update global matrices from local matrices
!------------------------------------------------------------------------------
              IF ( TransientSimulation ) THEN
                 MASS = 0.d0
                 CALL Default1stOrderTime( MASS, STIFF, FORCE )
              END IF
          
              CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
           END DO   ! Neumann & Newton BCs
!------------------------------------------------------------------------------

           CALL DefaultFinishAssembly()
           CALL DefaultDirichletBCs()

           CALL Info( SolverName, 'Assembly done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
           at = CPUTime() - at
           st = CPUTime()

           PrevNorm = Solver % Variable % Norm

           Norm = DefaultSolve()

           st = CPUTIme()-st
           totat = totat + at
           totst = totst + st
           WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
           CALL Info( SolverName, Message, Level=4 )
           WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
           CALL Info( SolverName, Message, Level=4 )

!           IF (LimitSolution) THEN
!              PRINT *, 'NNNNNNNNNNNNNNNNNNNNNNNNNOOOOOOOOOOOOOOOOOOOOOOOOOORRRRRRRRRRRRRRRRRRRRRRRRRMMMMMMMMMMMMMMM'
!              Norm = 0.0d00
!              DO i= 1,M
!                 IF ((LagrangeMultiplier1(i)==0.0d00) .AND. &
!                     (LagrangeMultiplier2(i)==0.0d00)) &
!                      Norm = Norm + SQRT(Temp(TempPerm(i))**2 )
!              END DO
!              Norm = Norm/M
!              Solver % Variable % Norm = Norm
!           END IF

!------------------------------------------------------------------------------


           IF ( PrevNorm + Norm /= 0.0d0 ) THEN
              RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
           ELSE
              RelativeChange = 0.0d0
           END IF

           WRITE( Message, * ) 'Result Norm   : ',Norm
           CALL Info( SolverName, Message, Level=4 )
           WRITE( Message, * ) 'Relative Change : ',RelativeChange
           CALL Info( SolverName, Message, Level=4 )

           IF ( RelativeChange < NewtonTol .OR. iter > NewtonIter ) &
                NewtonLinearization = .TRUE.

           IF ( RelativeChange < NonlinearTol ) EXIT
           

           CALL ListAddConstReal(Solver % Values,  &
                'Nonlinear System Relaxation Factor', Relax )

!------------------------------------------------------------------------------
        END DO ! of the nonlinear iteration
!------------------------------------------------------------------------------

        CALL ListAddConstReal(Solver % Values,  &
             'Nonlinear System Tolerance', NonlinearTol)

!------------------------------------------------------------------------------
!   Compute cumulative time done by now and time remaining
!------------------------------------------------------------------------------
        IF ( .NOT. TransientSimulation ) EXIT
        CumulativeTime = CumulativeTime + dt
        dt = Timestep - CumulativeTime
     END DO ! time interval

!------------------------------------------------------------------------------
     CALL  ListAddConstReal( Solver % Values,  &
          'Nonlinear System Relaxation Factor', SaveRelax )
!------------------------------------------------------------------------------

     DEALLOCATE( PrevSolution )

!   IF ( ListGetLogical( Solver % Values, 'Adaptive Mesh Refinement', Found ) ) &
!      CALL RefineMesh( Model,Solver,Temp,TempPerm, &
!            HeatInsideResidual, HeatEdgeResidual, HeatBoundaryResidual )

     SubroutineVisited = .TRUE.


   CONTAINS

SUBROUTINE HomologousTempBoundary(  StiffMatrix, ForceVector, LoadVector, &
       NodalC1,  Element, n, Permeable )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), ForceVector(:), LoadVector(:), NodalC1(:)
    INTEGER :: n
    LOGICAL :: Permeable
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: SqrtElementMetric, s, C1, Load
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3)
    LOGICAL :: Stat
    INTEGER :: i,j,p,q,t,dim,CoordSys
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

    CALL GetElementNodes( Nodes )
    IP = GaussPoints( Element )

    StiffMatrix = 0.0d0
    ForceVector = 0.0d0

    ! Over all  Gauss-points
    
    DO t = 1, IP % n
       stat = ElementInfo( Element, Nodes, IP % u(t), IP % v(t), IP % w(t), &
            SqrtElementMetric, Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )
       s = SqrtElementMetric * IP % s(t)      
       Load = SUM( LoadVector(1:N)*Basis(1:N) ) ! Load at Gauss-point
       PRINT *, t, ' out of total ', IP % n, SqrtElementMetric, IP % s(t), Load 
       DO q = 1,n
          ForceVector(q) = ForceVector(q) + s * Basis(q) * Load
       END DO
    END DO

  END SUBROUTINE HomologousTempBoundary

!!$  FUNCTION GetDivergence(Field, Element, ElementNodes, n, Nodes) RESULT(divergence)
!!$    USE DefUtils
!!$    USE TYPES
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! external variables
!!$    REAL(KIND=dp) :: Field(:,:), divergence
!!$    INTEGER :: n, Nodes
!!$    TYPE(Nodes_t), TARGET    :: ElementNodes
!!$    TYPE(Element_t), TARGET :: Element
!!$
!!$    !internal variables
!!$    REAL(KIND=dp) :: SqrtElementMetric
!!$    REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3),ddBasisddx(n,3,3)
!!$    INTEGER :: DIM, j
!!$    LOGICAL :: stat
!!$    
!!$    DIM = CoordinateSystemDimension()
!!$
!!$    stat = ElementInfo( Element,ElementNodes, &
!!$         Element % Type % NodeU(n), Element % Type % NodeV(n), Element % Type % NodeW(n), &
!!$         SqrtElementMetric, Basis,dBasisdx,ddBasisddx, .FALSE. )
!!$
!!$    divergence = 0.0d00
!!$ 
!!$    DO j=1,DIM
!!$       divergence = divergence + SUM( dBasisdx(1:N,j) * Field(j,1:N) )
!!$    END DO
!!$  
!!$  END FUNCTION GetDivergence

  FUNCTION GetScalarGradient(Scalar, Element, ElementNodes, NodeNumber, Nodes, ElementNo) RESULT(gradient)
    USE DefUtils
    USE TYPES
    USE CoordinateSystems

    IMPLICIT NONE

    ! external variables
    REAL(KIND=dp) :: Scalar(:), gradient(3)
    INTEGER :: NodeNumber, Nodes, ElementNo
    TYPE(Nodes_t), TARGET    :: ElementNodes
    TYPE(Element_t), TARGET :: Element

    !internal variables
    REAL(KIND=dp) :: SqrtElementMetric
    REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3),ddBasisddx(n,3,3)
    INTEGER :: DIM, j
    LOGICAL :: stat, CylindricSymmetry


    DIM = CoordinateSystemDimension()
    CylindricSymmetry = (CurrentCoordinateSystem() == CylindricSymmetric .OR. &
         CurrentCoordinateSystem() == AxisSymmetric)
    IF ( CurrentCoordinateSystem() == Cartesian .OR. CylindricSymmetry) THEN
       stat = ElementInfo( Element,ElementNodes, &
            Element % Type % NodeU(NodeNumber), &
            Element % Type % NodeV(NodeNumber), &
            Element % Type % NodeW(NodeNumber), &
            SqrtElementMetric, Basis, dBasisdx, ddBasisddx, .FALSE. )

       IF (.NOT.stat) THEN
          WRITE(Message,'(a, i5, a)') 'Degenerated element ', ElementNo, ' found'
          CALL FATAL('GetScalarGradient',Message)
       END IF

       gradient = 0.0d00
 

       DO j=1,DIM
          gradient(j) = SUM( dBasisdx(1:NodeNumber,j) * Scalar(1:NodeNumber) )
       END DO
    ELSE
       CALL FATAL('GetScalarGradient', 'Only Cartesian 2D, CArtesian 3D, AxisSummetric and CylindricSymmetric implemented')
    END IF

  END FUNCTION GetScalarGradient


!------------------------------------------------------------------------------
  END SUBROUTINE HomologousTempSolver
!------------------------------------------------------------------------------



