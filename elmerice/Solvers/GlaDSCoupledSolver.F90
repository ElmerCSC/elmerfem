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
! *  Authors: Olivier Gagliardini, Mauro Werder 
! *  Email:   olivier.gagliardini@ujf-grenoble.fr, m_werder@sfu.ca 
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - Scientific Computing Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 15 February 2013
! *
! *****************************************************************************/
!> Solve for the sheet hydraulic Potential, sheet thickness and channels area
!> similutaneously  (GlaDS model) - This solver replace the 3 solvers solving
!> for these 3 variables independently.  
!>  Equations defined in Werder et al., 2013 
   RECURSIVE SUBROUTINE GlaDSCoupledsolver( Model,Solver,Timestep,TransientSimulation )
!******************************************************************************
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
     USE Differentials
     USE MaterialModels
     USE DefUtils
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    External variables
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Solver_t), POINTER :: PointerToSolver
     TYPE(Matrix_t), POINTER :: Systemmatrix
     TYPE(Nodes_t) :: ElementNodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Edge, Face, Bulk
     TYPE(ValueList_t), POINTER :: Equation, Material, SolverParams, BodyForce, BC, Constants
     TYPE(Variable_t), POINTER :: ChannelAreaVar, ChannelFluxVar, SheetThicknessVar, &
          GroundedMaskVar, HydPotVar
     TYPE(Mesh_t), POINTER :: Mesh 
     
     INTEGER :: i, j, k, l, m, n, t, iter, body_id, eq_id, material_id, &
          istat, LocalNodes,bf_id, bc_id,  DIM, dimSheet, iterC, &
          NSDOFs, NonlinearIter, GhostNodes, NonlinearIterMin, Ne, BDForder, &
          CoupledIter, Nel, ierror, ChannelSolver, FluxVariable, ThicknessSolver, ierr

     TYPE(Variable_t), POINTER :: HydPotSol
     TYPE(Variable_t), POINTER :: ThickSol, AreaSol, VSol, WSol, NSol,  &
              PwSol, ZbSol, qSol, hstoreSol, QcSol, QmSol

     INTEGER, POINTER :: NodeIndexes(:), HydPotPerm(:), PwPerm(:), ZbPerm(:), &
             ThickPerm(:), VPerm(:), WPerm(:), NPerm(:), AreaPerm(:), & 
             qPerm(:), hstorePerm(:), QcPerm(:), QmPerm(:),&
             CAPerm(:), CFPerm(:), SHPerm(:)

     REAL(KIND=dp), POINTER :: HydPot(:), HydPotPrev(:,:), ForceVector(:)
     REAL(KIND=dp), POINTER :: ThickSolution(:), ThickPrev(:,:), VSolution(:), WSolution(:), &
            NSolution(:), PwSolution(:), AreaSolution(:), AreaPrev(:,:), ZbSolution(:), &
            qSolution(:), hstoreSolution(:), QcSolution(:), QmSolution(:),&
            CAValues(:), CFValues(:), SHValues(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, SolverName, MaskName
     CHARACTER(LEN=MAX_NAME_LEN) :: SheetThicknessName, ChannelAreaName, ZbName
     CHARACTER(LEN=MAX_NAME_LEN) :: methodSheet, methodChannels 

     LOGICAL :: Found, FluxBC, Channels, Storage, FirstTime = .TRUE., &
          AllocationsDone = .FALSE.,  SubroutineVisited = .FALSE., &
          meltChannels = .TRUE., NeglectH = .TRUE., Calving = .FALSE., &
          CycleElement=.FALSE., MABool = .FALSE., MaxHBool = .FALSE., LimitEffPres=.FALSE., &
          MinHBool=.FALSE., CycleNode=.FALSE.
     LOGICAL, SAVE :: UseGM, AllowSheetAtGL, ZeroSheetWithHP
     LOGICAL, ALLOCATABLE ::  IsGhostNode(:), NoChannel(:), NodalNoChannel(:)

     ! For use in masking GlaDS floating shelves.  "MASK_HP" is for situations where
     ! Hydraulic potential should be set to zero but not the sheet thickness.  This is
     ! to allow non zero sheet outflow across the grounding line.
     INTEGER :: MaskStatus
     INTEGER, PARAMETER :: MASK_ALL = 0, MASK_NONE = 1, MASK_HP = 2
     
     REAL(KIND=dp) :: NonlinearTol, dt, CumulativeTime, RelativeChange, &
          Norm, PrevNorm, S, C, Qc, MaxArea, MaxH, MinH
     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:), SheetConductivity(:), ChannelConductivity(:),&
       FORCE(:),  C1(:), CT(:), OldValues(:), Refq(:)

     REAL(KIND=dp), ALLOCATABLE :: Vvar(:), ublr(:), hr2(:)

     REAL(KIND=dp), ALLOCATABLE :: lr(:), hr(:), Ar(:), Wopen(:), &
         ub(:), Snn(:), Ev(:), &
         ng(:), alphas(:), betas(:), betac(:), Phi0(:), Phim(:), Afactor(:), Bfactor(:), &
         MoulinArea(:), MoulinFlux(:)

     REAL(KIND=dp), ALLOCATABLE :: IceDensity(:), Ac(:), alphac(:), CCt(:), &
         CCw(:), lc(:)
     REAL(KIND=dp) :: ChannelArea, WaterDensity, gravity, Lw, EdgeTangent(3), &
                 ALPHA, BETA, CoupledTol, CoupledNorm, PrevCoupledNorm, &
                 Discharge(3)
   
     REAL(KIND=dp) :: totat, st, totst, t1
     REAL(KIND=dp) :: Np, pw, zb, Wo, he
     
     REAL(KIND=dp) :: at, at0

     TYPE(ValueHandle_t) :: Load_h
     
     SAVE &
          ElementNodes, EdgeNodes,      &
          C1,                    &
          CT,                    &
          FirstTime,             &
          SheetConductivity,      &
          ChannelConductivity,      &
          MASS,                  &
          STIFF,LOAD,            &
          FORCE,                 &
          IsGhostNode,           &
          AllocationsDone, VariableName, SolverName, NonLinearTol, M, &
          lr, hr, Ar, Wopen, Afactor, Bfactor, &
          MoulinArea, MoulinFlux, &
          ub, Snn, Ev, WaterDensity, gravity, &
          ng, alphas, betas, betac, Phi0, Phim, SheetThicknessName, &
          ChannelAreaName, ZbName, IceDensity, Ac, alphac, CCt, &
          CCw, lc, Lw, NoChannel, NodalNoChannel, &
          Channels, meltChannels, NeglectH, BDForder, &
          Vvar, ublr, hr2, Refq, Nel,&
          Calving, Load_h, LimitEffPres, MaskName

      
     totst = 0.0_dp
     totat = 0.0_dp

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     VariableName = TRIM(Solver % Variable % Name)
     SolverName = 'GlaDSCoupledsolver ('// VariableName // ')'

     CALL ListInitElementKeyword( Load_h, 'Body Force', TRIM(Solver % Variable % Name) // ' Volume Source')

     
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
     SystemMatrix => Solver % Matrix
     ForceVector => Solver % Matrix % RHS

     PointerToSolver => Solver

     HydPotSol => Solver % Variable
     HydPotPerm  => HydPotSol % Perm
     HydPot => HydPotSol % Values
     HydPotPrev => HydPotSol % PrevValues
     
     LocalNodes = COUNT( HydPotPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

     !CHANGE (and at all DIMs)
     Mesh => Solver % Mesh
     DIM = Mesh % MeshDim
     M = Mesh % NumberOfNodes

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Mesh % Changed ) THEN
        N = Mesh % MaxElementNodes
        Ne = Mesh % NumberOfEdges
        Nel = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        K = SIZE( SystemMatrix % Values )
        L = SIZE( SystemMatrix % RHS )

        IF ( AllocationsDone ) THEN
           DEALLOCATE(                    &
                ElementNodes % x,         &
                ElementNodes % y,         &
                ElementNodes % z,         &
                EdgeNodes % x,         &
                EdgeNodes % y,         &
                EdgeNodes % z,         &
                C1,                       &
                CT,                       &
                SheetConductivity,         &
                ChannelConductivity,         &
                MASS,                     &
                STIFF,LOAD,               &
                FORCE,                    &
                IsGhostNode,              &
                lr, hr, Ar, Wopen, Afactor, Bfactor, &
                MoulinArea, MoulinFlux, &
                ub, Snn, Ev, &
                ng, alphas, betas, betac, Phi0, Phim, &
                IceDensity, Ac, alphac, CCt, &
                CCw, lc, OldValues, NoChannel, NodalNoChannel, &
                Vvar, ublr, hr2, &
                Refq)

        END IF                           
        
        ALLOCATE(                                  &
             ElementNodes % x( N ),                &
             ElementNodes % y( N ),                &
             ElementNodes % z( N ),                &
             EdgeNodes % x( N ),                &
             EdgeNodes % y( N ),                &
             EdgeNodes % z( N ),                &
             C1( N ),                              &
             CT( N ),                              &
             SheetConductivity( N ),            &
             ChannelConductivity( N ),            &
             MASS( N, N ),                     &
             STIFF( N, N ), LOAD( N ),           &
             FORCE( N ),                         &
             IsGhostNode( M ),                     &
             lr(N), hr(N), Ar(N), Wopen(N), Afactor(N), Bfactor(N), &
             MoulinArea(N), MoulinFlux(N), &
             ub(N), Snn(N), Ev(N), &
             ng(N), alphas(N), betas(N), betac(N), Phi0(N), Phim(N),   &
             IceDensity(N), Ac(N), alphac(N), CCt(N), &
             CCw(N), lc(N), OldValues(K), NoChannel(M), NodalNoChannel(N), &
             Vvar(M), ublr(M), hr2(M), &
             refq(dim*M), &
             STAT=istat)

        IF ( istat /= 0 ) THEN
           CALL FATAL( SolverName, 'Memory allocation error' )
        ELSE
           CALL INFO(SolverName, 'Memory allocation done', level=1 )
        END IF
        
        ! get the ghost nodes of this partition
        IF ( ParEnv % PEs > 1 ) THEN !only if we have a parallel run
           IsGhostNode( 1:M ) = .FALSE.
           GhostNodes = 0;
           DO t=1,Solver % NumberOfActiveElements
              Element => GetActiveElement(t)
              IF (ParEnv % myPe == Element % partIndex) CYCLE
              DO i=1,GetElementNOFNodes(Element)
                 IsGhostNode(Element % NodeIndexes(i)) = .TRUE.
                 GhostNodes = GhostNodes + 1
              END DO
           END DO
        PRINT *, ParEnv % myPe, ':', GhostNodes, ' ghost nodes'
        END IF

        ! Find the nodes for which we have no channel (on the boundary)
        ! Default is False - We allow Channel to growth everywhere
        NoChannel = .False.
        DO t=1, Mesh % NumberOfBoundaryElements
           ! get element information
           Element => GetBoundaryElement(t)
           !IF ( .NOT.ActiveBoundaryElement() ) CYCLE
           IF (ParEnv % PEs > 1) THEN
              IF (ParEnv % myPe /= Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1)) CYCLE
           END IF
           n = GetElementNOFNodes()
           IF ( GetElementFamily() == 1 ) CYCLE
   
           NULLIFY(BC)
           BC => GetBC( Element )
           bc_id = GetBCId( Element )
           NodalNoChannel = .False.
           NodalNoChannel(1:n) =  GetLogical(BC, 'No Channel BC', Found)
           IF (Found) NoChannel(Element % NodeIndexes(1:n)) = NodalNoChannel(1:n)
        END DO

        AllocationsDone = .TRUE.
     END IF


!------------------------------------------------------------------------------
!    Read physical and numerical constants and initialize 
!------------------------------------------------------------------------------
     IfFirstTime: IF (FirstTime) THEN
        FirstTime = .FALSE.
        Constants => GetConstants()

        WaterDensity = ListGetConstReal( Constants, 'Fresh Water Density', Found )
        IF (.NOT. Found) THEN           
           WaterDensity = ListGetConstReal( Constants, 'Water Density', Found )
           IF (Found) THEN
              WRITE(Message,'(A)') 'Constant name >Water Density< has been &
                   replaced with >Fresh Water Density< due to naming conflict'
              CALL WARN(SolverName, Message )
           END IF
           CALL FATAL(SolverName, 'Constant >Fresh Water Density< not found')
        END IF
        
        gravity = ListGetConstReal( Constants, 'Gravity Norm', UnFoundFatal=.TRUE. )
        Lw = ListGetConstReal( Constants, 'Latent Heat', UnFoundFatal=.TRUE. ) 

        ChannelAreaName = GetString( Constants, &
            'Channel Area Variable Name', Found )
        IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Channel Area Variable Name< not found in section Constants')
           CALL WARN(SolverName,'Taking default value >Channel Area<')
           WRITE(ChannelAreaName,'(A)') 'Channel Area'
        END IF

        SheetThicknessName = GetString( Constants, &
            'Sheet Thickness Variable Name', Found )
        IF(.NOT.Found) THEN        
           CALL WARN(SolverName,'Keyword >Sheet Thickness Variable Name< not found in section Constants')
           CALL WARN(SolverName,'Taking default value >Sheet Thickness<')
           WRITE(SheetThicknessName,'(A)') 'Sheet Thickness'
        END IF

        ZbName = GetString( Constants, 'Bedrock Variable Name', Found )
        IF (Found) THEN
           ZbSol => VariableGet( Solver % Mesh % Variables, ZbName, UnFoundFatal=.TRUE. )
        ELSE
           CALL WARN(SolverName,'Keyword >Bedrock Variable Name< not found in section Constants')
           CALL WARN(SolverName,'Taking default value >zb<')
           WRITE(ZbName,'(A)') 'Zb'
        END IF

        ! To get Channel variables added to this solver mesh if doing
        ! calving and hydrology and consequently having many meshes
        Calving = ListGetLogical(Model % Simulation, 'Calving', Found)
        IF(.NOT.Found) Calving = .FALSE.

        ! Default behaviour relating to marine ice sheets and unglaciated grounded areas is to set the
        ! following switches to false. The defaults change to true when using Samuel Cook's "Calving" 
        ! (set in simulation section of sif).  The defaults will be overwritten for each of the switches
        ! that are specified in the solver section of the sif.       
        SolverParams => GetSolverParams()
        
        UseGM = GetLogical( SolverParams,'Use GroundedMask', Found )
        IF (.NOT. Found) THEN
           IF (Calving) THEN              
              UseGM = .TRUE.
           ELSE
              UseGM = .FALSE.
           END IF
        END IF
        IF (UseGM) THEN
           MaskName = GetString( SolverParams, 'Mask Name', Found )
           IF (.NOT. Found) THEN
              MaskName = "GroundedMask"
           END IF
        END IF
        
        AllowSheetAtGL = GetLogical( SolverParams,'Allow Sheet At GL', Found )
        IF (.NOT. Found) THEN
           AllowSheetAtGL = .TRUE.
        END IF
        ZeroSheetWithHP = GetLogical( SolverParams,'Zero Sheet With HP', Found )
        IF (.NOT. Found) THEN
          IF (Calving) THEN              
            ZeroSheetWithHP = .TRUE.
          ELSE
            ZeroSheetWithHP = .FALSE.
          END IF
        END IF

        IfCalving: IF(Calving) THEN
          DO i=1,Model % NumberOfSolvers
            IF(Model % Solvers(i) % Variable % Name == ChannelAreaName) THEN 
              ChannelSolver = i
              EXIT
            END IF
          END DO
          ChannelAreaVar => VariableGet(Model % Solvers(ChannelSolver) % Mesh&
                     % Variables, ChannelAreaName, ThisOnly=.TRUE.)
          ALLOCATE(CAPerm(SIZE(ChannelAreaVar % Perm)), CAValues(SIZE(ChannelAreaVar % Values)))
          CAPerm = ChannelAreaVar % Perm
          CAValues = ChannelAreaVar % Values
          CALL VariableAdd(Mesh % Variables, Mesh, Solver,&
               'Channel Area', 1, CAValues, CAPerm)
          ChannelAreaVar => VariableGet(Mesh % Variables, 'Channel Area',&
                      ThisOnly=.TRUE.)
          ALLOCATE(ChannelAreaVar % PrevValues(SIZE(ChannelAreaVar % Values),MAX(Solver&
                   % Order, Solver % TimeOrder)))
          ChannelAreaVar % PrevValues(:,1) = ChannelAreaVar % Values
          NULLIFY(ChannelAreaVar)
            
          ChannelFluxVar => VariableGet(Model % Solvers(ChannelSolver) % Mesh&
                     % Variables, 'Channel Flux', ThisOnly=.TRUE.)
          ALLOCATE(CFPerm(SIZE(ChannelFluxVar % Perm)), CFValues(SIZE(ChannelFluxVar % Values)))
          CFPerm = ChannelFluxVar % Perm
          CFValues = ChannelFluxVar % Values
          CALL VariableAdd(Mesh % Variables, Mesh, Solver,&
               'Channel Flux', 1, CFValues, CFPerm) 
          ChannelFluxVar => VariableGet(Mesh % Variables, 'Channel Flux',&
                      ThisOnly=.TRUE.)
          ALLOCATE(ChannelFluxVar % PrevValues(SIZE(ChannelFluxVar % Values),MAX(Solver&
                   % Order, Solver % TimeOrder)))
          ChannelFluxVar % PrevValues(:,1) = ChannelFluxVar % Values
          NULLIFY(ChannelFluxVar)

          !The same for sheet thickness
          DO i=1,Model % NumberOfSolvers
            IF(Model % Solvers(i) % Variable % Name == SheetThicknessName) THEN 
              ThicknessSolver = i
              EXIT
            END IF
          END DO
          SheetThicknessVar => VariableGet(Model % Solvers(ThicknessSolver) % Mesh&
                     % Variables, SheetThicknessName, ThisOnly=.TRUE.)
          ALLOCATE(SHPerm(SIZE(SheetThicknessVar % Perm)), SHValues(SIZE(SheetThicknessVar % Values)))
          SHPerm = SheetThicknessVar % Perm
          SHValues = SheetThicknessVar % Values !Needed to reflect initial condition
          CALL VariableAdd(Mesh % Variables, Mesh, Solver,&
               'Sheet Thickness', 1, SHValues, SHPerm)
          SheetThicknessVar => VariableGet(Mesh % Variables, 'Sheet Thickness',&
                      ThisOnly=.TRUE.)
          ALLOCATE(SheetThicknessVar % PrevValues(SIZE(SheetThicknessVar % Values),MAX(Solver&
                   % Order, Solver % TimeOrder)))
          SheetThicknessVar % PrevValues(:,1) = SheetThicknessVar % Values
          !Necessary to ensure initial condition value reflected in PrevValues
          SheetThicknessVar % PrevValues(:,1) = SheetThicknessVar % Values
          NULLIFY(SheetThicknessVar)
        END IF IfCalving

        ! TODO : implement higher order BDF method
        BDForder = GetInteger(GetSimulation(),'BDF Order', Found)
        IF (.NOT.Found) BDForder = 1
        IF (BDForder /= 1) THEN
           WRITE(Message,'(a)') 'Only working for BDF = 1' 
           CALL FATAL(SolverName, Message)
        END IF
     END IF IfFirstTime

     SolverParams => GetSolverParams()

     NeglectH = GetLogical( SolverParams,'Neglect Sheet Thickness in Potential', Found )
     IF ( .NOT.Found ) THEN
        CALL FATAL(SolverName, 'No >Neglect Sheet Thickness in Potential< found')
     END IF

     methodSheet = GetString( SolverParams, 'Sheet Integration Method', Found )
     IF(.NOT.Found) THEN        
        CALL FATAL(SolverName, 'No >Sheet Integration Methods< found')
     ELSE
        IF ((methodSheet /='explicit').AND.(methodSheet /='implicit').AND.(methodSheet /= 'crank-nicolson')) THEN
           CALL FATAL(SolverName, 'Sheet Integration Method: Implicit, Explicit or Crank-Nicolson')
        END IF 
     END IF

     Channels = GetLogical( SolverParams,'Activate Channels', Found )
     IF (.NOT. Found) Channels = .FALSE.

     !CHANGE
     !To pick up channel and sheet size limiters, if specified
     MaxArea  = GetConstReal( SolverParams, &
          'Max Channel Area',    MABool )
     IF ((.NOT. MABool)) CALL WARN(SolverName,'No max channel area specified. &
          Channels may grow very large')

     LimitEffPres = GetLogical( SolverParams, &
          'Limit Negative Effective Pressure', Found)
     IF (.NOT.Found) LimitEffPres= .FALSE.
     
     MaxH  = GetConstReal( SolverParams, &
          'Max Sheet Thickness',    MaxHBool )
     IF ((.NOT. MaxHBool)) CALL WARN(SolverName,'No max sheet thickness specified.&
          Sheet may grow very large')

     MinH  = GetConstReal( SolverParams, &
          'Min Sheet Thickness',    MinHBool )
  

     IF (Channels) THEN
        meltChannels = GetLogical( SolverParams,'Activate Melt From Channels', Found )
        IF ( .NOT.Found ) THEN
           CALL FATAL(SolverName, 'No >Activate Melt From Channels< found')
        END IF

        methodChannels = GetString( SolverParams, 'Channels Integration Method', Found )
        IF(.NOT.Found) THEN        
           CALL FATAL(SolverName, 'No >Channels Integration Methods< found')
        ELSE
           IF ((methodChannels /= 'explicit').AND.(methodChannels /= 'implicit').AND.(methodChannels /= 'crank-nicolson')) THEN
              CALL FATAL(SolverName, 'Channels Integration Method: Implicit, Explicit or Crank-Nicolson')
           END IF 
        END IF
     END IF

     NonlinearIter = GetInteger( SolverParams, &
                     'Nonlinear System Max Iterations', Found )
     IF ( .NOT.Found ) NonlinearIter = 1

     NonlinearIterMin = GetInteger(   SolverParams, &
                     'Nonlinear System Min Iterations', Found )
     IF ( .NOT.Found ) NonlinearIterMin = 1
 
     NonlinearTol  = GetConstReal( SolverParams, &
          'Nonlinear System Convergence Tolerance',    Found )
     IF ((.Not.Found).AND.(NonlinearIter>1)) CALL FATAL(SolverName,'Need >Nonlinear System Convergence Tolerance<')

     CoupledIter = GetInteger( SolverParams, &
                    'Coupled Max Iterations', Found)
     IF ( .NOT.Found ) CoupledIter = 1

     CoupledTol  = GetConstReal( SolverParams, &
          'Coupled Convergence Tolerance',    Found )
     IF ((.Not.Found).AND.(CoupledIter>1)) CALL FATAL(SolverName,'Need >Nonlinear System Convergence Tolerance<')
     
     ThickSol => VariableGet( Mesh % Variables, SheetThicknessName, UnfoundFatal = .TRUE. )
     ThickPerm     => ThickSol % Perm
     ThickSolution => ThickSol % Values
     ThickPrev => ThickSol % PrevValues

     IF (Channels) THEN
        AreaSol => VariableGet( Mesh % Variables, ChannelAreaName, UnfoundFatal = .TRUE. )
        AreaPerm     => AreaSol % Perm
        AreaSolution => AreaSol % Values
        AreaPrev => AreaSol % PrevValues
        
        ! flux in the channels (for output only) - edge type variable
        QcSol => VariableGet( Mesh % Variables, "Channel Flux" )
        IF ( ASSOCIATED( QcSol ) ) THEN
           QcPerm     => QcSol % Perm
           QcSolution => QcSol % Values
        END IF
     END IF

     ! discharge out of the moulins (for output only)
     QmSol => VariableGet( Mesh % Variables, "Flux from Moulins" )
     IF ( ASSOCIATED( QmSol ) ) THEN
        QmPerm     => QmSol % Perm
        QmSolution => QmSol % Values
     END IF

     ZbSol => VariableGet( Mesh % Variables, ZbName )
     IF ( ASSOCIATED( ZbSol ) ) THEN
        ZbPerm     => ZbSol % Perm
        ZbSolution => ZbSol % Values
     END IF

     dt = Timestep

!------------------------------------------------------------------------------
!    Read the other variables needed                  
!------------------------------------------------------------------------------
     VSol => VariableGet( Mesh % Variables, 'Vclose' )
     IF ( ASSOCIATED( VSol ) ) THEN
          VPerm     => VSol % Perm
          VSolution => VSol % Values
     END IF

     WSol => VariableGet( Mesh % Variables, 'Wopen' )
     IF ( ASSOCIATED( WSol ) ) THEN
          WPerm     => WSol % Perm
          WSolution => WSol % Values
     END IF

     NSol => VariableGet( Mesh % Variables, 'Effective Pressure' )
     IF ( ASSOCIATED( NSol ) ) THEN
          NPerm     => NSol % Perm
          NSolution => NSol % Values
     END IF

     PwSol => VariableGet( Mesh % Variables, 'Water Pressure' )
     IF ( ASSOCIATED( PwSol ) ) THEN
          PwPerm     => PwSol % Perm
          PwSolution => PwSol % Values
     END IF

     qSol => VariableGet( Mesh % Variables, 'Sheet Discharge' )
     IF ( ASSOCIATED( qSol ) ) THEN
          qPerm     => qSol % Perm
          qSolution => qSol % Values
     END IF

     hstoreSol => VariableGet( Mesh % Variables, 'Sheet Storage' )
     IF ( ASSOCIATED( hstoreSol ) ) THEN
          hstorePerm     => hstoreSol % Perm
          hstoreSolution => hstoreSol % Values
     END IF

!------------------------------------------------------------------------------
!   Loop for the coupling of the three equations
!------------------------------------------------------------------------------
    ! check on the coupled Convergence is done on the potential solution only
    PrevCoupledNorm = ComputeNorm( Solver, SIZE(HydPot), HydPot ) 

    
    DO iterC = 1, CoupledIter

!------------------------------------------------------------------------------
!       non-linear system iteration loop
!------------------------------------------------------------------------------
        DO iter = 1, NonlinearIter
           !------------------------------------------------------------------------------
           ! print out some information
           !------------------------------------------------------------------------------
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
              
           Ne = Mesh % NumberOfEdges
            
           !------------------------------------------------------------------------------
           ! lets start
           !------------------------------------------------------------------------------
           CALL DefaultInitialize()
           body_id = -1
           NULLIFY(Material)
           !------------------------------------------------------------------------------
           ! Bulk elements (the sheet layer)
           !------------------------------------------------------------------------------
           DO t=1,Solver % NumberOfActiveElements
              !------------------------------------------------------------------------------
              ! write some info on status of assembly
              !------------------------------------------------------------------------------
              IF ( RealTime() - at0 > 1.0 ) THEN
                 WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
                      (Solver % NumberOfActiveElements-t) / &
                      (1.0*Solver % NumberOfActiveElements)), ' % done'

                 CALL Info( SolverName, Message, Level=5 )
                 at0 = RealTime()
              END IF
              !------------------------------------------------------------------------------
              ! Check if this element belongs to a body where scalar equation
              ! should be calculated
              !------------------------------------------------------------------------------
              Element => GetActiveElement(t,Solver)
              
              ! cycle halo elements
              !-------------------
              IF (ParEnv % myPe  /=  Element % partIndex) CYCLE

              IF (.NOT.ASSOCIATED(Element)) CYCLE
              IF ( Element % BodyId /= body_id ) THEN
                 body_id = Element % bodyId
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
              ! Check if Element dimension is 2 (in x, y plane)
              !------------------------------------------------------------------------------
              dimSheet = Element % TYPE % DIMENSION
              IF(dimSheet>2) THEN
                  WRITE (Message,'(A,I3)')' Work only for 1D or 2D elements' 
                  CALL FATAL(SolverName,Message)
              END IF

              !------------------------------------------------------------------------------
              ! Get element material parameters
              !------------------------------------------------------------------------------              
              N = GetElementNOFNodes(Element)
              CALL GetElementNodes( ElementNodes )

!----------------------------------------------------              
! Get parameters to compute the Total Conductivity 
! K = k . f(h) . |grad HydPot |^(beta-2)
! k = SheetConductivity 
! f(h) = h^alphas
!
! And read parameters to evaluate v - w
! w = u_b/l_r (h_r - h) 
! N = p_i + rhow . g . b - HydPot
! v = A . h . | N |^{n-1} . N 
!----------------------------------------------------              

             CALL GetParametersSheet( Element, Material, N, SheetConductivity, alphas, &
                 betas, Ev, ub, Snn, lr, hr, Ar, ng ) 

              CT = 0.0_dp
              Wopen = 0.0_dp
              Phi0 = 0.0_dp
              DO i=1, n 
                 j = Element % NodeIndexes(i)
                 k = ThickPerm(j)
                 IF ( ASSOCIATED( ZbSol )) THEN
                    zb = ZbSolution(ZbPerm(j))
                 ELSE 
                    IF (dimSheet==1) THEN  
                       zb = Mesh % Nodes % y(j)
                    ELSE
                       zb = Mesh % Nodes % z(j)
                    END IF 
                 END IF
                 CT(i) = Ev(i) /( WaterDensity * gravity)
                 Wopen(i) = MAX(ub(i) / lr(i) * (hr(i) - ThickSolution(k)), 0.0)
 
                 Phi0(i) = Snn(i) + gravity*WaterDensity*zb
                 IF (.NOT. NeglectH) THEN 
                   Phi0(i) = Phi0(i) + gravity*WaterDensity*ThickSolution(k)
                 END IF
              END DO

              !------------------------------------------------------------------------------
              ! Add body forces
              !------------------------------------------------------------------------------
              LOAD = 0.0_dp              
              
              BodyForce => GetBodyForce()
              !IF ( ASSOCIATED( BodyForce ) ) THEN
              !   bf_id = GetBodyForceId()
              !   LOAD(1:N) = LOAD(1:N) + &
              !     GetReal( BodyForce, TRIM(Solver % Variable % Name) // ' Volume Source', Found )
              !END IF
              ! f = m - w + v
              ! v is not added here as it will be linearized for the assembly
              LOAD(1:N) = LOAD(1:N) - Wopen(1:N)

              !------------------------------------------------------------------------------
              ! Get element local matrices, and RHS vectors
              !------------------------------------------------------------------------------
              MASS = 0.0_dp
              STIFF = 0.0_dp 
              FORCE = 0.0_dp
              ! cartesian coords
              !----------------
              IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                  CALL SheetCompose( MASS, STIFF, FORCE, LOAD, &
                       ThickSolution(ThickPerm(Element % NodeIndexes(1:n))), &
                       HydPot(HydPotPerm(Element % NodeIndexes(1:N))), &
                       CT(1:N), SheetConductivity(1:n), alphas(1:n), betas(1:n), & 
                       Phi0(1:n), Ar(1:n), ng(1:n), Element, N, ElementNodes ) 
              ELSE
                 WRITE(Message,'(A)')' Work only for cartesian coordinate'
                 CALL FATAL( SolverName, Message)
              END IF              
              !------------------------------------------------------------------------------
              ! If time dependent simulation add mass matrix to stiff matrix
              !------------------------------------------------------------------------------
              IF ( TransientSimulation ) THEN
                 CALL Default1stOrderTime( MASS, STIFF, FORCE )
              END IF

              CALL DefaultUpdateEquations( STIFF, FORCE )
           END DO     !  Bulk elements

           ! TODO: Is this really needed? 
            CALL DefaultFinishBulkAssembly()
        
        !------------------------------------------------------------------------------
        ! Edge element (the channels)
        ! Edge element are created in the SIF file using Element = 'n=1 e:1' 
        !------------------------------------------------------------------------------
        ! Go only if Activate Channels is True
        !------------------------------------------------------------------------------
        IF (Channels) THEN
           body_id = -1
           NULLIFY(Material)
           DO t=1, Mesh % NumberOfEdges 
              Edge => Mesh % Edges(t)
              IF (.NOT.ASSOCIATED(Edge)) CYCLE
              IF (ParEnv % PEs > 1) THEN
                IF (ParEnv % myPe /= Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1)) CYCLE
              END IF
              n = Edge % TYPE % NumberOfNodes

              ! Work only for 202 elements => n=2
              IF (n/=2) CALL FATAL(SolverName, 'Work only for edge element of type 202')
              ! We keep only the edge which belong in the sheet surface
              ! i.e. Both nodes have Perm > 0
              IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE

              ! We check if we are in a boundary where we want No channel
              IF (ALL(NoChannel(Edge % NodeIndexes(1:n)))) CYCLE


              EdgeNodes % x(1:n) = Mesh % Nodes % x(Edge % NodeIndexes(1:n))
              EdgeNodes % y(1:n) = Mesh % Nodes % y(Edge % NodeIndexes(1:n))
              EdgeNodes % z(1:n) = Mesh % Nodes % z(Edge % NodeIndexes(1:n))

              ! Compute the unit tangent vector of the edge
              EdgeTangent(1) = EdgeNodes % x(2) - EdgeNodes % x(1)
              EdgeTangent(2) = EdgeNodes % y(2) - EdgeNodes % y(1)
              EdgeTangent(3) = EdgeNodes % z(2) - EdgeNodes % z(1)
              EdgeTangent = EdgeTangent / SQRT(SUM(EdgeTangent*EdgeTangent))

              ! Read in the Material Section the needed Parameters for the
              ! Channels - 
              ! In case DIM = 2 Faces have a material ASSOCIATED
              ! In case DIM = 3, Bulk elements have material ASSOCIATED
              IF (DIM==2) THEN 
                 Bulk => Edge % BoundaryInfo % Left
                 IF (.Not.ASSOCIATED(Bulk)) Bulk => Edge % BoundaryInfo % Right
              ELSE 
                 Face => Edge % BoundaryInfo % Left
                 IF (.Not.ASSOCIATED(Face)) Face => Edge % BoundaryInfo % Right
                 IF (ASSOCIATED(Face)) THEN
                    Bulk => Face % BoundaryInfo % Left 
                    IF (.Not.ASSOCIATED(Bulk)) Bulk => Face % BoundaryInfo % Right
                 END IF
              END IF
              IF (.Not.ASSOCIATED(Bulk)) THEN
                 WRITE(Message,'(a,i0)')'No face or bulk element associated to edge no:', t
                 CALL FATAL(SolverName, Message)
              END IF

              IF ( Bulk % BodyId /= body_id ) THEN
                 body_id = Bulk % bodyId
                 Material => GetMaterial( Bulk )
                 IF (.NOT.ASSOCIATED(Material)) THEN
                    WRITE (Message,'(A,I3)') 'No Material found for edge no. ', t
                    CALL FATAL(SolverName,Message)
                 ELSE
                    material_id = GetMaterialId( Bulk, Found)
                    IF(.NOT.Found) THEN
                       WRITE (Message,'(A,I3)') 'No Material ID found for edge no. ', t
                       CALL FATAL(SolverName,Message)
                    END IF
                 END IF
              END IF
              
!----------------------------------------------------              
! Get parameters to compute the Total Conductivity 
! Kc = kc . fc(h) . |grad HydPot |^(betac-2)
! kc = ChannelConductivity 
! fc(h) = S^alphac
! and the closure and opening velocity for Channels
!----------------------------------------------------              

              CALL GetParametersChannel( Edge, Material, n, SheetConductivity, &
                 ChannelConductivity, alphac, betac, alphas, betas, IceDensity, &
                 Snn, Ac, ng, CCt, CCw, lc ) 

              Phi0 = 0.0_dp
              Afactor = 0.0_dp
              Bfactor = 0.0_dp
              Phi0 = 0.0_dp
              Phim = 0.0_dp
              DO i=1,N
                 IF ( ASSOCIATED( ZbSol )) THEN
                    zb = ZbSolution(ZbPerm(Edge % NodeIndexes(i)))
                 ELSE 
                    IF (dimSheet==1) THEN
                       zb = Mesh % Nodes % y(Edge % NodeIndexes(i))
                    ELSE
                       zb = Mesh % Nodes % z(Edge % NodeIndexes(i))
                    END IF
                 END If
                 Phim(i) = gravity*WaterDensity*zb
                 Phi0(i) = Snn(i) + Phim(i) 
                 IF (.NOT. NeglectH) THEN
                    k = ThickPerm(Edge % NodeIndexes(i))
                    Phi0(i) = Phi0(i) + gravity*WaterDensity*ThickSolution(k)
                 END IF
                 Afactor(i) = CCt(i) * CCw(i) * WaterDensity
                 Bfactor(i) = 1.0/(Lw * IceDensity(i)) 
                 IF (meltChannels) Bfactor(i) = Bfactor(i) - 1.0/(Lw * WaterDensity)
              END DO

! The variable Channel Area is defined on the edges only
              IF (AreaPerm(M+t) <=   0)  CYCLE 

              !CHANGE
              !To stabilise channels
              IF(MABool) THEN
                IF(AreaSolution(AreaPerm(M+t)) > MaxArea) AreaSolution(AreaPerm(M+t)) = MaxArea
              END IF
              ChannelArea = AreaSolution(AreaPerm(M+t))

              !------------------------------------------------------------------------------
              ! Get element local matrices, and RHS vectors
              !------------------------------------------------------------------------------
              MASS = 0.0_dp
              STIFF = 0.0_dp
              FORCE = 0.0_dp
              ! cartesian coords
              !----------------

              IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                 CALL ChannelCompose( MASS, STIFF, FORCE, &
                      ThickSolution(ThickPerm(Edge % NodeIndexes(1:n))), &
                      HydPot(HydPotPerm(Edge % NodeIndexes(1:N))), ChannelArea, &
                      ChannelConductivity, alphac, betac, Phi0, Phim, Ac, lc, ng, &
                      SheetConductivity, alphas, betas, Afactor, Bfactor, EdgeTangent, &
                      Edge, n, EdgeNodes )
              ELSE
                 WRITE(Message,'(A)')' Work only for cartesian coordinate'
                 CALL FATAL( SolverName, Message)
              END IF              
              !CHANGE
              !To stop weird channel instability where some channels grow
              !exponentially to stupid levels and eventually mess up whole
              !mesh. Usually seems to be channels with limited fluxes that
              !shouldn't do much, so resetting to 0 seems safest
              !TODO Come up with a better way of fixing this
              IF(MABool) THEN
                k = AreaPerm(M+t)
                IF(AreaSolution(k) > MaxArea) AreaSolution(k) = MaxArea
              END IF
              ! This should be not needed as MASS = 0 here
              IF ( TransientSimulation ) THEN
                 CALL Default1stOrderTime( MASS, STIFF, FORCE, Edge )
              END IF

              CALL DefaultUpdateEquations( STIFF, FORCE, Edge )
           END DO ! Edge elements 
        END IF

           
           !------------------------------------------------------------------------------
           ! Neumann & Newton boundary conditions
           !------------------------------------------------------------------------------
           DO t=1, Mesh % NumberOfBoundaryElements
              ! get element information
              Element => GetBoundaryElement(t)
              ! Don't test this as BC on the side are not active for this solver
              !IF ( .NOT.ActiveBoundaryElement() ) CYCLE

              ! cycle halo elements
              !--------------------
              IF (ParEnv % myPe  /=  Element % partIndex) CYCLE
              n = GetElementNOFNodes()

              IF ( GetElementFamily() == 1 ) THEN 
              ! Moulin case (flux at node)
              
                BC => GetBC( Element )
                bc_id = GetBCId( Element )
                CALL GetElementNodes( ElementNodes )
                IF ( ASSOCIATED( BC ) ) THEN            
                  STIFF=0.0_dp
                  FORCE=0.0_dp
                  MASS=0.0_dp

                  Storage =  GetLogical(BC,'Moulin Storage', Found)
                  IF (Storage) THEN
                    MoulinArea(1:N) = ListGetReal( BC, 'Moulin Area',  N, Element % NodeIndexes, Found, &
                         UnfoundFatal = .TRUE. )
                    ! MASS is a scalar here
                    MASS(1,1) = MoulinArea(1)/(WaterDensity*gravity)
                  END IF

                  ! Is there surface input
                  MoulinFlux(1:N) = ListGetReal( BC, 'Moulin Flux',  N, Element % NodeIndexes, Found, &
                         UnfoundFatal = .TRUE. )
                  FORCE(1) = MoulinFlux(1)

                  ! If variable exist, update the Flux from Moulins variable 
                  IF (ASSOCIATED( QmSol )) THEN
                     j = Element % NodeIndexes(1)
                     IF ( HydPotPerm(j) <=  0 ) CYCLE 
                     QmSolution(QmPerm(j)) = MoulinFlux(1)-MASS(1,1)*(HydPot(HydPotPerm(j))-HydPotPrev(HydPotPerm(j),1))/dt 
                  END IF

                  !------------------------------------------------------------------------------
                  ! Update global matrices from local matrices
                  !------------------------------------------------------------------------------
                  ! We need the Default1stOrderUpdate for the Moulin
                  IF ( TransientSimulation ) THEN
                     CALL Default1stOrderTime( MASS, STIFF, FORCE, Element )
                  END IF
                  CALL DefaultUpdateEquations( STIFF, FORCE, Element )
                END IF

              ELSE
              ! flux over the sheet
                BC => GetBC( Element )
                bc_id = GetBCId( Element )
                CALL GetElementNodes( ElementNodes )

                IF ( ASSOCIATED( BC ) ) THEN            
                  ! Check that we are on the correct boundary part!
                  STIFF=0.0_dp 
                  FORCE=0.0_dp
                  MASS=0.0_dp

                  LOAD=0.0_dp
                  FluxBC =  GetLogical(BC,TRIM(Solver % Variable % Name) // ' Flux BC', Found)

                  IF (FluxBC) THEN
                     !---------------
                     !BC: -k@T/@n = q
                     !---------------
                     LOAD(1:N)  = LOAD(1:N) + &
                             GetReal( BC, 'Sheet Water Flux', Found )
                     ! -------------------------------------
                     ! set boundary due to coordinate system
                     ! -------------------------------------
                     IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                        CALL SheetBoundary( STIFF, FORCE, &
                            LOAD, Element, n, ElementNodes )
                     ELSE
                        WRITE(Message,'(A)')' Work only for cartesian coordinate'
                        CALL FATAL( SolverName, Message)
                     END IF
                     !!! TODO : do we need that as MASS = 0 here
                     IF ( TransientSimulation ) THEN
                        CALL Default1stOrderTime( MASS, STIFF, FORCE, Element )
                     END IF

                     CALL DefaultUpdateEquations( STIFF, FORCE, Element )
                  END IF
                END IF
              END IF ! Element 1O1
          
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
           Norm = 0.0_dp

           PrevNorm = Solver % Variable % Norm
           Norm = DefaultSolve()

           st = CPUTime()-st
           totat = totat + at
           totst = totst + st
           WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
           CALL Info( SolverName, Message, Level=4 )
           WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
           CALL Info( SolverName, Message, Level=4 )


           IF ( PrevNorm + Norm /= 0.0d0 ) THEN
              RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
           ELSE
              RelativeChange = 0.0d0
           END IF

           WRITE( Message, * ) 'Result Norm   : ',Norm
           CALL Info( SolverName, Message, Level=4 )
           WRITE( Message, * ) 'Relative Change : ',RelativeChange
           CALL Info( SolverName, Message, Level=4 )

    !      !----------------------
    !      ! check for convergence
    !      !----------------------
           IF ( RelativeChange < NonlinearTol .AND. iter > NonlinearIterMin ) EXIT 
        END DO ! of the nonlinear iteration

!------------------------------------------------------------------------------
!       Update the Sheet Thickness                 
!------------------------------------------------------------------------------
           Elements: DO t=1,Solver % NumberOfActiveElements
              Element => GetActiveElement(t,Solver)
              IF (ParEnv % myPe  /=  Element % partIndex) CYCLE

              IF (.NOT.ASSOCIATED(Element)) CYCLE
              IF ( Element % BodyId /= body_id ) THEN
                 body_id = Element % bodyId
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

              N = GetElementNOFNodes(Element)
              CALL GetElementNodes( ElementNodes )

              
              IF (UseGM) THEN
                ! Cycle elements with ungrounded nodes and zero all hydrology variables
                CycleElement = .FALSE.

                DO i=1, N
                  MaskStatus = ProcessMask(MaskName, AllowSheetAtGL, Element % NodeIndexes(i))
                  SELECT CASE (MaskStatus)
                  CASE (MASK_ALL)
                    CycleElement = .TRUE.
                    WSolution(WPerm(Element % NodeIndexes(i))) = 0.0
                    Vvar(Element % NodeIndexes(i)) = 0.0
                    NSolution(NPerm(Element % NodeIndexes(i))) = 0.0
                    hstoreSolution(hstorePerm(Element % NodeIndexes(i))) = 0.0
                  CASE (MASK_HP)
                    NSolution(NPerm(Element % NodeIndexes(i))) = 0.0
                  CASE (MASK_NONE)
                  CASE DEFAULT
                    WRITE(Message,'(A)') "MaskStatus not recognised"
                    CALL FATAL( SolverName, Message)
                  END SELECT
                END DO
                IF (CycleElement) THEN
                  CYCLE
                END IF
              END IF

              CALL GetParametersSheet( Element, Material, N, SheetConductivity, alphas, &
                  betas, Ev, ub, Snn, lr, hr, Ar, ng ) 

              CT = 0.0_dp
              Phi0 = 0.0_dp
              DO i=1, n 
                 j = Element % NodeIndexes(i)
                 k = ThickPerm(j)
                 IF ( ASSOCIATED( ZbSol )) THEN
                    zb = ZbSolution(ZbPerm(j))
                 ELSE 
                    IF (dimSheet==1) THEN  
                       zb = Mesh % Nodes % y(j)
                    ELSE
                       zb = Mesh % Nodes % z(j)
                    END IF 
                 END IF
                 CT(i) = Ev(i) /( WaterDensity * gravity)
                 Phi0(i) = Snn(i) + gravity*WaterDensity*zb
                 Np = Phi0(i) - HydPot(HydPotPerm(j))  
                 IF (.NOT. NeglectH) THEN
                    Np = Np + WaterDensity*gravity*ThickSolution(k)
                 END IF
                 pw = HydPot(HydPotPerm(j)) - gravity*WaterDensity*zb
                 Vvar(j) = Ar(i)*ABS(Np)**(ng(i)-1.0)*Np
                 Wo = MAX(ub(i) / lr(i) * (hr(i) - ThickSolution(k)), 0.0) 
                 he = Ev(i)*(HydPot(HydPotPerm(j))/(WaterDensity*gravity)-zb)
                 ublr(j) = ub(i)/lr(i)
                 hr2(j) = hr(i)

                 !To stop it working out values for non-ice covered parts of a
                 !hydromesh in a coupled calving-hydro simulation
                 IF ( ZeroSheetWithHP ) THEN
                   IF(Snn(i)==0.0) THEN
                     Np = 0.0
                     he = 0.0
                   END IF
                 END IF 

                 ! Save this for output only
                 IF (ASSOCIATED(WSol)) WSolution(WPerm(j)) = Wo
                 IF (ASSOCIATED(NSol)) NSolution(NPerm(j)) = Np
                 IF (ASSOCIATED(PwSol)) PwSolution(PwPerm(j)) = pw
                 IF (ASSOCIATED(hstoreSol)) hstoreSolution(hstorePerm(j)) = he
              END DO
           END DO Elements     !  Bulk elements

           ! Loop over all nodes to update ThickSolution
           DO j = 1, Mesh % NumberOfNodes
              k = ThickPerm(j)
              IF (k==0) CYCLE

              CycleNode = .FALSE.
              IF (UseGM) THEN
                ! Cycle ungrounded nodes and zero hydrology variables
                MaskStatus = ProcessMask(MaskName, AllowSheetAtGL, Element % NodeIndexes(i))
                SELECT CASE (MaskStatus)
                CASE (MASK_ALL)
                  CycleNode = .TRUE.
                CASE (MASK_HP, MASK_NONE)
                CASE DEFAULT
                  WRITE(Message,'(A)') "MaskStatus not recognised"
                  CALL FATAL( SolverName, Message)
                END SELECT
              END IF
              IF (ZeroSheetWithHP) THEN
                HydPotVar => VariableGet(Mesh % Variables, "hydraulic potential", ThisOnly=.TRUE., UnfoundFatal=.FALSE.)
                IF(HydPotVar % Values( HydPotVar % perm(j) ).EQ.0.0) THEN
                  CycleNode = .TRUE.
                END IF
                NULLIFY(HydPotVar)
              END IF
              IF (CycleNode) THEN
                ThickSolution(k) = 0.0
                ThickPrev(k,1) = 0.0
                CYCLE
              END IF
              
              IF(MaxHBool) THEN
                IF (ThickSolution(k)>MaxH) THEN
                  ThickSolution(k) = MaxH
                  !ThickPrev(k,1) = 0.0
                END IF
              END IF

              IF(MinHBool) THEN
                IF (ThickSolution(k)<MinH) THEN
                  ThickSolution(k) = MinH
                END IF
              END IF
              
              SELECT CASE(methodSheet)
              CASE('implicit') 
                 IF (ThickSolution(k) > hr2(j)) THEN
                    ThickSolution(k) = MAX(ThickPrev(k,1)/(1.0_dp + dt*Vvar(j)) , AEPS)
                 ELSE
                    ThickSolution(k) = MAX((ThickPrev(k,1) + dt*ublr(j)*hr2(j))/(1.0_dp + dt*(Vvar(j)+ublr(j))) , AEPS)
                 END IF
              CASE('explicit')
                 IF (ThickSolution(k) > hr2(j)) THEN 
                    ThickSolution(k) = MAX(ThickPrev(k,1)*(1.0_dp - dt*Vvar(j)) , AEPS)
                 ELSE
                    ThickSolution(k) = MAX(ThickPrev(k,1)*(1.0_dp - dt*(Vvar(j)+ublr(j))) + dt*ublr(j)*hr2(j), AEPS)
                 END IF
              CASE('crank-nicolson')
                 IF (ThickSolution(k) > hr2(j)) THEN 
                    ThickSolution(k) = MAX(ThickPrev(k,1)*(1.0_dp - 0.5*dt*Vvar(j))/(1.0_dp + 0.5*dt*Vvar(j)) , AEPS)
                 ELSE
                    ThickSolution(k) = MAX((ThickPrev(k,1)*(1.0_dp - 0.5*dt*(Vvar(j)+ublr(j))) + dt*ublr(j)*hr2(j)) & 
                                & /(1.0_dp + 0.5*dt*(Vvar(j)+ublr(j))) , AEPS)
                 END IF
              END SELECT 

              ! Update Vvar
              Vvar(j) = Vvar(j) * ThickSolution(k)

           END DO 
!------------------------------------------------------------------------------
!       Update the Channels Area                 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!       non-linear system iteration loop
!------------------------------------------------------------------------------
     IF (Channels) THEN 
        PrevNorm = ChannelAreaNorm()
        
        DO iter = 1, NonlinearIter
              Edges: DO t=1, Mesh % NumberOfEdges 
                 Edge => Mesh % Edges(t)
                 IF (.NOT.ASSOCIATED(Edge)) CYCLE
                 IF (ParEnv % PEs > 1) THEN
                   IF (ParEnv % myPe  /=  Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1)) CYCLE
                 END IF
                 n = Edge % TYPE % NumberOfNodes
                 IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE
                 IF (ALL(NoChannel(Edge % NodeIndexes(1:n)))) CYCLE

                 IF (UseGM) THEN
                   ! Cycle ungrounded nodes and zero hydrology variables
                   CycleElement = .FALSE.
                   DO i=1, n
                     MaskStatus = ProcessMask(MaskName, AllowSheetAtGL, Element % NodeIndexes(i))
                     SELECT CASE (MaskStatus)
                     CASE (MASK_ALL)
                       CycleElement = .TRUE.
                     CASE (MASK_HP, MASK_NONE)
                     CASE DEFAULT
                       WRITE(Message,'(A)') "MaskStatus not recognised"
                       CALL FATAL( SolverName, Message)
                     END SELECT
                   END DO
                   IF(CycleElement) THEN
                     AreaSolution(AreaPerm(M+t)) = 0.0
                     QcSolution(QcPerm(M+t)) = 0.0
                     CYCLE
                   END IF
                 END IF
                 
                 EdgeNodes % x(1:n) = Mesh % Nodes % x(Edge % NodeIndexes(1:n))
                 EdgeNodes % y(1:n) = Mesh % Nodes % y(Edge % NodeIndexes(1:n))
                 EdgeNodes % z(1:n) = Mesh % Nodes % z(Edge % NodeIndexes(1:n))

                 ! Compute the unit tangent vector of the edge
                 EdgeTangent(1) = EdgeNodes % x(2) - EdgeNodes % x(1)
                 EdgeTangent(2) = EdgeNodes % y(2) - EdgeNodes % y(1)
                 EdgeTangent(3) = EdgeNodes % z(2) - EdgeNodes % z(1)
                 EdgeTangent = EdgeTangent / SQRT(SUM(EdgeTangent*EdgeTangent))

                 IF (DIM==2) THEN 
                    Bulk => Edge % BoundaryInfo % Left
                    IF (.Not.ASSOCIATED(Bulk)) Bulk => Edge % BoundaryInfo % Right
                 ELSE 
                    Face => Edge % BoundaryInfo % Left
                    IF (.Not.ASSOCIATED(Face)) Face => Edge % BoundaryInfo % Right
                    IF (ASSOCIATED(Face)) THEN
                       Bulk => Face % BoundaryInfo % Left 
                       IF (.Not.ASSOCIATED(Bulk)) Bulk => Face % BoundaryInfo % Right
                    END IF
                 END IF
                 IF (.Not.ASSOCIATED(Bulk)) THEN
                    WRITE(Message,'(a,i0)')'No face or bulk element associated to edge no:', t
                    CALL FATAL(SolverName, Message)
                 END IF

                 IF ( Bulk % BodyId /= body_id ) THEN
                    body_id = Bulk % bodyId
                    Material => GetMaterial( Bulk )
                    IF (.NOT.ASSOCIATED(Material)) THEN
                       WRITE (Message,'(A,I3)') 'No Material found for edge no. ', t
                       CALL FATAL(SolverName,Message)
                    ELSE
                       material_id = GetMaterialId( Bulk, Found)
                       IF(.NOT.Found) THEN
                          WRITE (Message,'(A,I3)') 'No Material ID found for edge no. ', t
                          CALL FATAL(SolverName,Message)
                       END IF
                    END IF
                 END IF
              
                 CALL GetParametersChannel( Edge, Material, n, SheetConductivity, &
                    ChannelConductivity, alphac, betac, alphas, betas, IceDensity, &
                    Snn, Ac, ng, CCt, CCw, lc ) 

                 Phi0 = 0.0_dp
                 Afactor = 0.0_dp
                 Bfactor = 0.0_dp
                 Phi0 = 0.0_dp
                 Phim = 0.0_dp
                 DO i=1,N
                    IF ( ASSOCIATED( ZbSol )) THEN
                       zb = ZbSolution(ZbPerm(Edge % NodeIndexes(i)))
                    ELSE 
                       IF (dimSheet==1) THEN
                          zb = Mesh % Nodes % y(Edge % NodeIndexes(i))
                       ELSE
                          zb = Mesh % Nodes % z(Edge % NodeIndexes(i))
                       END IF
                    END If
                    Phim(i) = gravity*WaterDensity*zb
                    Phi0(i) = Snn(i) + Phim(i) 
                    IF (.NOT. NeglectH) THEN
                       k = ThickPerm(Edge % NodeIndexes(i))
                       Phi0(i) = Phi0(i) + gravity*WaterDensity*ThickSolution(k)
                    END IF 
                    Afactor(i) = CCt(i) * CCw(i) * WaterDensity
                    Bfactor(i) = 1.0/(Lw * IceDensity(i)) 
                 END DO

                 IF (AreaPerm(M+t) <= 0 ) CYCLE
                 ChannelArea = AreaSolution(AreaPerm(M+t))

! Compute the force term to evolve the channels area
! Equation of the form dS/dt = S x ALPHA + BETA
                 IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                    CALL GetEvolveChannel( ALPHA, BETA, Qc, ChannelArea, &
                         HydPot(HydPotPerm(Edge % NodeIndexes(1:n))), &
                         ThickSolution(ThickPerm(Edge % NodeIndexes(1:n))), &
                         alphac, betac, ChannelConductivity, Phi0, Phim, Ac, lc, ng, &
                         SheetConductivity, alphas, betas, Afactor, Bfactor, &
                         EdgeTangent, Edge, n, EdgeNodes, LimitEffPres)
                 ELSE
                    WRITE(Message,'(A)')' Work only for cartesian coordinate'
                    CALL FATAL( SolverName, Message)
                 END IF

                 k = AreaPerm(M+t)
                 IF ( k <=  0 ) CYCLE  
                 SELECT CASE(methodChannels)
                 CASE('implicit') 
                    AreaSolution(k) = (AreaPrev(k,1) + dt*BETA)/(1.0_dp - dt*ALPHA)
                 CASE('explicit')
                    AreaSolution(k) = AreaPrev(k,1)*(1.0_dp + dt*ALPHA) + dt*BETA 
                 CASE('crank-nicolson')
                    AreaSolution(k) = (AreaPrev(k,1)*(1.0_dp + 0.5*ALPHA*dt) + dt*BETA)/(1.0_dp - 0.5*dt*ALPHA)
                 END SELECT 

                 !CHANGE
                 !To stop weird channel instability where some channels grow
                 !exponentially to stupid levels and eventually mess up whole
                 !mesh. Usually seems to be channels with limited fluxes that
                 !shouldn't do much, so resetting to 0 seems safest
                 !TODO Come up with a better way of fixing this
                 IF(MABool) THEN
                   IF(AreaPrev(k,1)  /=  0.0) THEN
                     IF(AreaSolution(k)>1.0 .AND. (AreaSolution(k)/AreaPrev(k,1))>5.0) THEN
                       AreaSolution(k) = AreaPrev(k,1)
                     END IF
                   END IF
                   IF(AreaSolution(k) > MaxArea) AreaSolution(k) = MaxArea
                   IF(ISNAN(AreaSolution(k))) AreaSolution(k) = 0.0
                 END IF

                 ! Save Qc if variable exists
                 IF (ASSOCIATED(QcSol)) THEN
                    IF ( QcPerm(M+t) <= 0 ) CYCLE
                    QcSolution(QcPerm(M+t)) = Qc
                 END IF
              END DO Edges
              
           Norm = ChannelAreaNorm()              
           t = Mesh % NumberOfEdges 

           IF ( PrevNorm + Norm /= 0.0d0 ) THEN
              RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
           ELSE
              RelativeChange = 0.0d0
           END IF

           PrevNorm = Norm 

           WRITE( Message, * ) 'S (NRM,RELC) : ',iter, Norm, RelativeChange
           CALL Info( SolverName, Message, Level=3 )

    !      !----------------------
    !      ! check for convergence
    !      !----------------------
           IF ( RelativeChange < NonlinearTol .AND. iter > NonlinearIterMin ) EXIT 
        END DO ! of the nonlinear iteration


        ! Make sure Area >= 0
        AreaSolution = MAX(AreaSolution, 0.0_dp )

        !Stop channels from expanding to eleventy-stupid
        IF(MABool) THEN
          AreaSolution = MIN(AreaSolution,MaxArea)
        END IF

     END IF  ! If Channels

      !   Check for convergence                           
      CoupledNorm = ComputeNorm( Solver, SIZE(HydPot), HydPot ) 

      IF ( PrevCoupledNorm + CoupledNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevCoupledNorm-CoupledNorm ) / (PrevCoupledNorm + CoupledNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF
      PrevCoupledNorm = CoupledNorm 

      WRITE( Message, * ) 'COUPLING LOOP (NRM,RELC) : ',iterC, CoupledNorm, RelativeChange
      CALL Info( SolverName, Message, Level=3 )

      IF ((RelativeChange < CoupledTol).AND. (iterC > 1)) EXIT 
   END DO ! iterC

!--------------------------------------------------------------------------------------------
! Output some useful variables
!--------------------------------------------------------------------------------------------
! Already computed variables
   IF (ASSOCIATED(VSol)) THEN 
      DO i=1, M
        IF ( VPerm(i) <=  0 ) CYCLE
        VSolution(VPerm(i)) = Vvar(i)
      ENDDO
   ENDIF 

! Output the sheet discharge (dimension = dimSheet)
   IF (ASSOCIATED(qSol)) THEN
      Refq = 0.0_dp
      qSolution = 0.0_dp

      ! Loop over all elements are we need to compute grad(Phi)
      DO t=1,Solver % NumberOfActiveElements
         !CHANGE - necessary if using a 2D mesh as is otherwise set to 1 as
         !boundary elements are last in first loop where it's set
         dimSheet = Element % TYPE % DIMENSION
         Element => GetActiveElement(t,Solver)
         IF (ParEnv % myPe /= Element % partIndex) CYCLE

         IF (.NOT.ASSOCIATED(Element)) CYCLE
         IF ( Element % BodyId /= body_id ) THEN
            body_id = Element % bodyId
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

         n = GetElementNOFNodes(Element)
         CALL GetElementNodes( ElementNodes )

         IF (UseGM) THEN
           ! Cycle ungrounded nodes and zero hydrology variables
           CycleElement = .FALSE.
           DO i=1, n
             MaskStatus = ProcessMask(MaskName, AllowSheetAtGL, Element % NodeIndexes(i))
             SELECT CASE (MaskStatus)
             CASE (MASK_ALL)
               CycleElement = .TRUE.
               DO j=1,dimSheet
                 k = dimSheet*(qPerm(Element % NodeIndexes(i))-1)+j
                 qSolution(k) = 0.0
                 Refq(k) = 0.0
               END DO
             CASE (MASK_HP, MASK_NONE)
             CASE DEFAULT
               WRITE(Message,'(A)') "MaskStatus not recognised"
               CALL FATAL( SolverName, Message)
             END SELECT
           END DO
           IF(CycleElement) CYCLE
         END IF
 
         ! we need the SheetConductivity, alphas, betas
         CALL GetParametersSheet( Element, Material, n, SheetConductivity, alphas, &
               betas, Ev, ub, Snn, lr, hr, Ar, ng ) 
         
          ! Go for all nodes of the element
          DO i=1,n
             Discharge = 0.0_dp
             CALL SheetDischargeCompute( & 
                   HydPot(HydPotPerm(Element % NodeIndexes(1:n))), &
                   ThickSolution(ThickPerm(Element % NodeIndexes(i))), &
                   SheetConductivity(i), alphas(i), betas(i), & 
                   Discharge, Element, n, ElementNodes, i ) 

             ! One more value for that node          
             DO j=1,dimSheet
                k = dimSheet*(qPerm(Element % NodeIndexes(i))-1)+j
                Refq(k) = Refq(k) + 1.0_dp
                qSolution(k) = qSolution(k) + Discharge(j)
             END DO  
          END DO
      END DO

      ! Mean nodal value
      DO i=1,n
         DO j=1,dimSheet
            k = dimSheet*(qPerm(Element % NodeIndexes(i))-1)+j
            IF ( Refq(k) > 0.0_dp ) THEN 
              qSolution(k) = qSolution(k)/Refq(k) 
            END IF
         END DO  
      END DO

   END IF

   SubroutineVisited = .TRUE.

   !CHANGE - to make sure PrevValues for added variables in calving updated
   IF(Calving) THEN
     SheetThicknessVar => VariableGet(Mesh % Variables, 'Sheet Thickness',ThisOnly=.TRUE.)
     SheetThicknessVar % PrevValues(:,1) = SheetThicknessVar % Values
     ChannelAreaVar => VariableGet(Mesh % Variables, 'Channel Area',ThisOnly=.TRUE.)
     ChannelAreaVar % PrevValues(:,1) = ChannelAreaVar % Values
     NULLIFY(SheetThicknessVar, ChannelAreaVar)
   END IF

CONTAINS    

  ! Use the grounded mask to decide how to mask the current node.
  !----------------------------------------------------------------------------------------------------------
  FUNCTION ProcessMask(MaskName, AllowSheetAtGL, ii) RESULT( MaskStatus_local )

    CHARACTER(LEN=MAX_NAME_LEN), INTENT(IN) :: MaskName
    LOGICAL, INTENT(IN)                     :: AllowSheetAtGL
    INTEGER, INTENT(IN)                     :: ii ! node index

    INTEGER :: MaskStatus_local

    MaskStatus_local = MASK_NONE
    
    GroundedMaskVar => VariableGet(Mesh % Variables, MaskName, ThisOnly=.TRUE., UnfoundFatal=.TRUE.)

    IF (GroundedMaskVar % Values(GroundedMaskVar % Perm(ii)).LT.0.0) THEN 
       MaskStatus_local = MASK_ALL
    ELSEIF (GroundedMaskVar % Values(GroundedMaskVar % Perm(ii)).EQ.0.0) THEN
       IF (AllowSheetAtGL) THEN
          MaskStatus_local = MASK_HP
       ELSE
          MaskStatus_local = MASK_ALL
       END IF
    END IF

    NULLIFY(GroundedMaskVar)
  
  END FUNCTION ProcessMask


  
  ! Compute consistent channel norm only considering the edges that also have hydrology defined on the nodes.
  ! In parallel only consider the edges in the partition where it is active.
  !----------------------------------------------------------------------------------------------------------
  FUNCTION ChannelAreaNorm() RESULT( AreaNorm ) 
    REAL(KIND=dp) :: AreaNorm
    
    INTEGER :: t,i,n,ncount
    REAL(KIND=dp) :: s, ssum, smin, smax

    ncount = 0
    smin = HUGE(smin)
    smax = 0.0_dp
    ssum = 0.0_dp
    
    DO t=1, Mesh % NumberOfEdges 
      i = AreaPerm(Mesh % NumberOfNodes + t)
      IF(i==0) CYCLE
      
      IF (ParEnv % PEs > 1) THEN
        IF (ParEnv % myPe /= Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1)) CYCLE
      END IF
      
      Edge => Mesh % Edges(t)
      n = Edge % TYPE % NumberOfNodes
      IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE

      s = AreaSolution(i)
      
      ssum = ssum + s**2
      smin = MIN(smin,s)
      smax = MAX(smax,s)
      
      ncount = ncount + 1      
    END DO
      
    IF( ParEnv % PEs > 1 ) THEN
      ncount = ParallelReduction(ncount)
      ssum = ParallelReduction(ssum)
      smin = ParallelReduction(smin,1)
      smax = ParallelReduction(smax,2)      
    END IF
    
    WRITE( Message, * ) 'Range S:', smin, smax
    CALL Info( SolverName, Message )
    
    AreaNorm = SQRT(ssum / ncount ) 
    
  END FUNCTION ChannelAreaNorm


  
!------------------------------------------------------------------------------
  SUBROUTINE SheetCompose( MassMatrix, StiffMatrix, ForceVector,  &
       LoadVector, NodalH, NodalHydPot, NodalCT, NodalC2, Nodalalphas, &
       Nodalbetas, NodalPhi0, NodalAr, NodalNg, Element, n, Nodes )
!------------------------------------------------------------------------------
    USE MaterialModels
    USE Integration
    USE Differentials

    IMPLICIT NONE

    REAL(KIND=dp), DIMENSION(:)   :: ForceVector, LoadVector
    REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix, StiffMatrix
    REAL(KIND=dp) :: NodalHydPot(:), Nodalalphas(:), NodalBetas(:)
    REAL(KIND=dp) :: NodalCT(:), NodalC2(:), NodalH(:)
    REAL(KIND=dp) :: NodalPhi0(:), NodalAr(:), NodalNg(:)

    INTEGER :: n

    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element

    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------

    REAL(KIND=dp) :: dBasisdx(n,3), detJ
    REAL(KIND=dp) :: Basis(n)

    REAL(KIND=dp) :: Force

    REAL(KIND=dp) :: A, M
    REAL(KIND=dp) :: Load

    REAL(KIND=dp) :: x, y, z

    INTEGER :: i, j, k, c, p, q, t, dim, N_Integ, NBasis

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    REAL(KIND=dp) :: s, u, v, w

    REAL(KIND=dp) :: CT, C2, Vfactor, Phi0, PhiG, LoadAtIP

    REAL(KIND=dp) :: gradPhi(3), Ngrad, na, nb, ng, hsheet 

    REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ, V_Integ, W_Integ, S_Integ

    LOGICAL :: Found, Transient, stat

    !------------------------------------------------------------------------------

    Transient = GetString(GetSimulation(),'Simulation type',Found)=='transient'

    ! the dimension for this solver is given by the dim of the elements
    ! should work for 2D and 3D general problem
    dim = Element % TYPE % DIMENSION

    ForceVector = 0.0_dp
    StiffMatrix = 0.0_dp
    MassMatrix  = 0.0_dp
    Load = 0.0_dp

    NBasis = n

    !------------------------------------------------------------------------------
    !    Integration stuff
    !------------------------------------------------------------------------------
    IntegStuff = GaussPoints( element )
    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n


    !------------------------------------------------------------------------------
    !    Now we start integrating
    !------------------------------------------------------------------------------
    DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

       !------------------------------------------------------------------------------
       !      Basis function values & derivatives at the integration point
       !------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )

       s = detJ * S_Integ(t)
       !------------------------------------------------------------------------------
       !      Coefficient of derivative terms
       !      at the integration point
       !------------------------------------------------------------------------------
       CT = SUM( NodalCT(1:n) * Basis(1:n) )
       !------------------------------------------------------------------------------
       !      Compute the effective conductivity K = k h |grad Phi|^{beta-2}
       !------------------------------------------------------------------------------
       gradPhi = 0.0_dp
       Do i=1, dim
          gradPhi(i) = SUM(NodalHydPot(1:n)*dBasisdx(1:n,i))
       END DO
       Ngrad = SQRT(SUM(gradPhi(1:dim)*gradPhi(1:dim)))
       IF (Ngrad < AEPS) Ngrad = AEPS
       nb = SUM(NodalBetas(1:n)*Basis(1:n))
       na = SUM(Nodalalphas(1:n)*Basis(1:n))
       hsheet = SUM(NodalH(1:n)*Basis(1:n))

       C2 = SUM( NodalC2(1:n) * Basis(1:n))
       C2 = C2 * hsheet**na 
       C2 = C2 * Ngrad**(nb-2.0)

       ! Vclose velocity is split into Vfactor . (Phi0 - Phi)^(n-1)
       PhiG = SUM(NodalHydPot(1:n)*Basis(1:n))
       Phi0 = SUM(NodalPhi0(1:n)*Basis(1:n))
       Vfactor = SUM(NodalAr(1:n)*Basis(1:n)) * hsheet
       ng = SUM(NodalNg(1:n)*Basis(1:n)) 
       Vfactor = Vfactor * ABS(Phi0-PhiG)**(ng-1.0_dp)

       Force = SUM( LoadVector(1:n)*Basis(1:n) )
       ! contribution from volume source (using handle)
       LoadAtIP = ListGetElementReal( Load_h, Basis, Element, Found, GaussPoint=t)
       IF (Found) THEN
         Force = Force + LoadAtIP
       END IF
       Force = Force + Vfactor * (Phi0 + (ng - 1.0_dp) * PhiG)

       !------------------------------------------------------------------------------
       !       Loop over basis functions of both unknowns and weights
       !------------------------------------------------------------------------------
       DO p=1,NBasis
          DO q=1,NBasis
             !------------------------------------------------------------------------------
             !         The diffusive-convective equation without stabilization
             !------------------------------------------------------------------------------
             M = CT * Basis(q) * Basis(p)
             !------------------------------------------------------------------------------
             !         The diffusion term
             !------------------------------------------------------------------------------
             A = C2 * SUM( dBasisdx(q,1:dim) * dBasisdx(p,1:dim))
             A = A + Vfactor * ng * Basis(q) * Basis(p)

             StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
             MassMatrix(p,q)  = MassMatrix(p,q)  + s * M
          END DO
          ForceVector(p) = ForceVector(p) + s * Force * Basis(p)
       END DO
    END DO

    !------------------------------------------------------------------------------
  END SUBROUTINE SheetCompose
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Return element local matrices and RSH vector for boundary conditions
!>  of diffusion convection equation: 
!------------------------------------------------------------------------------
   SUBROUTINE SheetBoundary( BoundaryMatrix, BoundaryVector, &
               LoadVector, Element, n, Nodes )
!------------------------------------------------------------------------------
  USE MaterialModels
  USE Integration
  USE Differentials

  IMPLICIT NONE
     REAL(KIND=dp) :: BoundaryMatrix(:,:), BoundaryVector(:), LoadVector(:)
     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n

     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3), detJ

     REAL(KIND=dp) :: U, V, W, S
     REAL(KIND=dp) :: Force
     REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)

     INTEGER :: i, t, q, p, N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------
     BoundaryVector = 0.0_dp
     BoundaryMatrix = 0.0_dp
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       U = U_Integ(t)
       V = V_Integ(t)
       W = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
                  Basis, dBasisdx )

       S = detJ * S_Integ(t)
!------------------------------------------------------------------------------
       Force = SUM( LoadVector(1:n)*Basis(1:n) )

       DO q=1,n
         BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
       END DO
     END DO
   END SUBROUTINE SheetBoundary

!------------------------------------------------------------------------------
SUBROUTINE ChannelCompose( MassMatrix, StiffMatrix, ForceVector, &
      NodalH, NodalHydPot, CArea, NodalKc, NodalAlphac, NodalBetac, NodalPhi0, NodalPhim, &
      NodalAc, Nodallc, Nodalng, NodalKs, NodalAlphas, NodalBetas, &  
      NodalAfactor, NodalBfactor, EdgeTangent, Element, n, Nodes )
!------------------------------------------------------------------------------
  USE MaterialModels
  USE Integration
  USE Differentials

  IMPLICIT NONE

     REAL(KIND=dp), DIMENSION(:)   :: ForceVector
     REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix, StiffMatrix
     REAL(KIND=dp) :: NodalHydPot(:), NodalBetas(:), NodalAlphas(:)
     REAL(KIND=dp) :: NodalH(:), CArea, EdgeTangent(3)
     REAL(KIND=dp) :: NodalPhi0(:), NodalPhim(:), NodalNg(:), NodalAfactor(:)
     REAL(KIND=dp) :: NodalKc(:), NodalAlphac(:), NodalAc(:), Nodallc(:)
     REAL(KIND=dp) :: NodalBetac(:), NodalKs(:), NodalBfactor(:)

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: dBasisdx(n,3), detJ
     REAL(KIND=dp) :: Basis(n), dBasisds(n)

     REAL(KIND=dp) :: Force

     REAL(KIND=dp) :: A, M
     REAL(KIND=dp) :: Load

     REAL(KIND=dp) :: x, y, z

     INTEGER :: i, j, k, c, p, q, t, N_Integ, NBasis

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp) :: s, u, v, w

     REAL(KIND=dp) :: Vfactor, Phi0, PhiG, Afactor, Bfactor, GradPhim, dPw, Ffactor

     REAL(KIND=dp) :: GradPhi, Ngrad, nbc, nbs, nas, nac, hsheet, ng, qc, Kc, Ks, lc
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ, V_Integ, W_Integ, S_Integ

     LOGICAL :: Vms, Found, Transient, stat
     TYPE(ValueList_t), POINTER :: BodyForce

!------------------------------------------------------------------------------

     Transient = GetString(GetSimulation(),'Simulation type',Found)=='transient'

     ForceVector = 0.0_dp
     StiffMatrix = 0.0_dp
     MassMatrix  = 0.0_dp
     Load = 0.0_dp

     NBasis = n

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
             Basis, dBasisdx )

       s = detJ * S_Integ(t)

! derivative along the edge
       DO p=1, n 
          dBasisds(p) = SUM(dBasisdx(p,1:3)*EdgeTangent(1:3))
       END DO 
!       
!------------------------------------------------------------------------------
!      Compute the effective conductivity K = k h^alpha |grad Phi|^{beta-2}
!      in the sheet ks and in the Channel Kc
! gradPhi is here a scalar 
!------------------------------------------------------------------------------
       GradPhi = SUM(NodalHydPot(1:n)*dBasisds(1:n))
       Ngrad = ABS(GradPhi)  
       IF (Ngrad < AEPS) Ngrad = AEPS
       nas = SUM(NodalAlphas(1:n)*Basis(1:n))
       nbs = SUM(NodalBetas(1:n)*Basis(1:n))
       nac = SUM(NodalAlphac(1:n)*Basis(1:n))
       nbc = SUM(NodalBetac(1:n)*Basis(1:n))
       hsheet = SUM(NodalH(1:n)*Basis(1:n))

       Ks = SUM( NodalKs(1:n)*Basis(1:n))
       Ks = Ks * hsheet**nas 
       Ks = Ks * Ngrad**(nbs-2.0) 

       Kc = SUM( NodalKc(1:n)*Basis(1:n))
       Kc = Kc * CArea**nac                        
       Kc = Kc * Ngrad**(nbc-2.0)  

! Vclose velocity is split into Ac . S . (Phi0 - Phi)^n       
       PhiG = SUM(NodalHydPot(1:n)*Basis(1:n))
       Phi0 = SUM(NodalPhi0(1:n)*Basis(1:n))

       ng = SUM(NodalNg(1:n)*Basis(1:n))
       Vfactor = CArea*SUM(NodalAc(1:n)*Basis(1:n))
       Vfactor = Vfactor*ABS(Phi0-PhiG)**(ng-1.0_dp)

       Afactor = SUM(NodalAfactor(1:n)*Basis(1:n))
       Bfactor = SUM(NodalBfactor(1:n)*Basis(1:n))
       Afactor = Afactor * Bfactor
       lc = SUM(Nodallc(1:n)*Basis(1:n))

       qc = -Ks * GradPhi 

       GradPhim = SUM(NodalPhim(1:n)*dBasisds(1:n)) 
       dPw = GradPhi - GradPhim

       IF ((CArea>0.0).OR.(qc*dPw>0.0)) THEN 
         Ffactor = Afactor * lc * qc 
       ELSE
         Ffactor = 0.0_dp
       END IF

       Bfactor = Bfactor * SIGN((ABS(-Kc*GradPhi)+ABS(lc*qc)),GradPhi)

       Force = Ffactor * GradPhim                
       Force = Force + Vfactor*(Phi0+(ng-1.0_dp)*PhiG)
       
!------------------------------------------------------------------------------
!       Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
        DO p=1,NBasis

          Load = Force * Basis(p)
          ForceVector(p) = ForceVector(p) + s * Load

        DO q=1,NBasis
!------------------------------------------------------------------------------
!         The diffusive-convective equation without stabilization
!------------------------------------------------------------------------------
          A = 0.0_dp
!------------------------------------------------------------------------------
!         The diffusion term
!------------------------------------------------------------------------------
          A = A + Kc * dBasisds(p) * dBasisds(q)
          A = A - Afactor * Kc * dPw * Basis(p) * dBasisds(q)
          A = A + Ffactor * Basis(p) * dBasisds(q)
          A = A + Vfactor*ng*Basis(p)*Basis(q)
          A = A + Bfactor * Basis(p) * dBasisds(q) 

          StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
        END DO
        END DO

      END DO
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE ChannelCompose
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE GetEvolveChannel(ALPHA, BETA, Qcc, CArea, NodalHydPot, NodalH, &
      NodalAlphac, NodalBetac, NodalKc, NodalPhi0, NodalPhim, NodalAc, Nodallc, Nodalng, &
      NodalKs, NodalAlphas, NodalBetas, NodalAfactor, NodalBfactor, &
      Tangent, Element, n, Nodes, LimitEffPres)
!------------------------------------------------------------------------------
  USE MaterialModels
  USE Integration
  USE Differentials

  IMPLICIT NONE

            
     REAL(KIND=dp) :: ALPHA, BETA, Qcc
     REAL(KIND=dp) :: NodalHydPot(:), NodalAlphas(:), NodalBetas(:)
     REAL(KIND=dp) :: NodalH(:), CArea, Tangent(3)
     REAL(KIND=dp) :: NodalPhi0(:), NodalPhim(:), NodalNg(:), NodalAfactor(:)
     REAL(KIND=dp) :: NodalKc(:), Nodalalphac(:), NodalAc(:), Nodallc(:)
     REAL(KIND=dp) :: NodalBetac(:), NodalKs(:), NodalBfactor(:)

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

     LOGICAL :: LimitEffPres

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: dBasisdx(n,3), detJ
     REAL(KIND=dp) :: Basis(n), dBasisds(n)


     REAL(KIND=dp) :: A, M
     REAL(KIND=dp) :: Load

     REAL(KIND=dp) :: x, y, z

     INTEGER :: i, j, k, c, p, q, t, N_Integ,NBasis

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp) :: s, u, v, w

     REAL(KIND=dp) :: Phi0, PhiG, EffPressatIP, Afactor, Bfactor, GradPhim, dPw, Ffactor

     REAL(KIND=dp) :: GradPhi, Ngrad, nbc, hsheet, nas, nbs, nac, ng, qc, Kc, Ks, lc
     
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ, V_Integ, W_Integ, S_Integ

     REAL(KIND=dp) :: Xi, Pii, Nef, Vc 

     LOGICAL :: Vms, Found, Transient, stat

!------------------------------------------------------------------------------

     Transient = GetString(GetSimulation(),'Simulation type',Found)=='transient'

     NBasis = n

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

! Compute the force at the middle of the Edge element (u=0)
     u = 0.0_dp
     v = 0.0_dp
     w = 0.0_dp

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
             Basis, dBasisdx )


! derivative along the edge
       DO p=1, n 
          dBasisds(p) = SUM(dBasisdx(p,1:3)*EdgeTangent(1:3))
       END DO 

!------------------------------------------------------------------------------
!      Compute the effective conductivity K = k h |grad Phi|^{beta-2}
!      in the sheet ks and in the Channel Kc
! gradPhi is here a scalar 
!------------------------------------------------------------------------------
! 
       gradPhi = SUM(NodalHydPot(1:n)*dBasisds(1:n))
       Ngrad = ABS(GradPhi) 
       IF (Ngrad < AEPS) Ngrad = AEPS

       nas = SUM(NodalAlphas(1:n)*Basis(1:n))
       nbs = SUM(NodalBetas(1:n)*Basis(1:n))
       nbc = SUM(NodalBetac(1:n)*Basis(1:n))
       nac  = SUM(NodalAlphac(1:n)*Basis(1:n))
       hsheet  = SUM(NodalH(1:n)*Basis(1:n))

       Ks = SUM( NodalKs(1:n) * Basis(1:n))
       Ks = Ks * hsheet**nas 
       Ks = Ks * Ngrad**(nbs-2.0_dp) 

       Kc = SUM( NodalKc(1:n) * Basis(1:n))
       Kc = Kc * MAX(CArea,0.0)**(nac-1.0_dp) 
       Kc = Kc * Ngrad**(nbc-2.0_dp)  

       PhiG = SUM(NodalHydPot(1:n)*Basis(1:n))
       Phi0 = SUM(NodalPhi0(1:n)*Basis(1:n))

       ng = SUM(NodalNg(1:n)*Basis(1:n))
       Vc = SUM(NodalAc(1:n)*Basis(1:n))

       IF (LimitEffPres) THEN
         EffPressatIP = MAX(Phi0-PhiG, 0.0_dp)
       ELSE
         EffPressatIP = Phi0-PhiG
       END IF
       
       Vc = Vc*ABS(EffPressatIP)**(ng-1.0_dp)
       Vc = Vc*(EffPressatIP)

       Afactor = SUM(NodalAfactor(1:n)*Basis(1:n))
       Bfactor = SUM(NodalBfactor(1:n)*Basis(1:n))
       lc = SUM(Nodallc(1:n)*Basis(1:n))

       qc = -Ks * GradPhi 

       GradPhim = SUM(NodalPhim(1:n)*dBasisds(1:n)) 
       dPw = GradPhi - GradPhim

       IF ((CArea>0.0).OR.(qc*dPw>0.0)) THEN 
         Ffactor = lc * qc 
       ELSE
         Ffactor = 0.0_dp
       END IF
  

       ! Terms in ALPHA
       Xi = ABS(-Kc*GradPhi*GradPhi)
       Pii = -Afactor*(-Kc*GradPhi)*dPw
       ALPHA = Bfactor*(Xi - Pii) - Vc 

       ! Terms in BETA
       Xi = ABS(lc*qc*GradPhi)
       Pii = -Afactor*(Ffactor)*dPw
       BETA = Bfactor*(Xi - Pii) 

       ! Channel flux for output
       Qcc = ABS(MAX(CArea,0.0)*Kc*GradPhi)

!------------------------------------------------------------------------------
END SUBROUTINE GetEvolveChannel
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE GetParametersSheet( Element, Material, N, SheetConductivity, alphas, &
                 betas, Ev, ub, Snn, lr, hr, Ar, ng ) 
!------------------------------------------------------------------------------
  USE MaterialModels
  USE Integration
  USE Differentials

  IMPLICIT NONE
  REAL(KIND=dp) :: SheetConductivity(:), alphas(:), betas(:), Ev(:), ub(:), &
           Snn(:), lr(:), hr(:), Ar(:), ng(:)
  INTEGER :: N
  LOGICAL :: Found = .FALSE.
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material

!------------------------------------------------------------------------------
  SheetConductivity = 0.0_dp
  SheetConductivity = ListGetReal(Material, 'Sheet Conductivity', n, Element % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. )
                   
  alphas = 0.0_dp
  alphas(1:N) =  ListGetReal( Material, 'Sheet Flow Exponent alpha', n, Element % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. )

  betas = 0.0_dp
  betas(1:N) =  ListGetReal( Material, 'Sheet Flow Exponent beta', n, Element % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. ) 
                   
  Ev = 0.0_dp
  Ev(1:N) =  ListGetReal( Material, 'Englacial Void Ratio', n, Element % NodeIndexes, & 
           Found, UnfoundFatal = .TRUE. ) 
                   
  ub = 0.0_dp
  ub(1:N) = ListGetReal( Material, 'Sliding Velocity',  n, Element % NodeIndexes, &
            Found, UnfoundFatal = .TRUE. )

  Snn = 0.0_dp
  Snn(1:N) = ListGetReal( Material, 'Ice Normal Stress',  N, Element % NodeIndexes, & 
           Found, UnfoundFatal = .TRUE. ) 
  IF (ANY(Snn(1:N) < 0.0)) THEN
     WRITE(Message,'(A)')'The Ice Normal Stress (overburden pressure) must be positive'
     CALL FATAL(SolverName, Message)
  END IF

  lr = 0.0_dp
  lr(1:N) = ListGetReal( Material, 'Bedrock Bump Length',  N, Element % NodeIndexes, & 
           Found, UnfoundFatal = .TRUE. ) 

  hr = 0.0_dp
  hr(1:N) = ListGetReal( Material, 'Bedrock Bump High',  N, Element % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. ) 

  Ar = 0.0_dp
  Ar(1:N) = ListGetReal( Material, 'Sheet Closure Coefficient',  N, Element % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. ) 

  ng = 0.0_dp
  ng(1:N) = ListGetReal( Material, 'Glen Exponent',  N, Element % NodeIndexes, & 
           Found, UnfoundFatal = .TRUE. ) 
!------------------------------------------------------------------------------
END SUBROUTINE GetParametersSheet 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE GetParametersChannel( Edge, Material, N, SheetConductivity, &
                 ChannelConductivity, alphac, betac, alphas, betas, IceDensity, &
                 Snn, Ac, ng, CCt, CCw, lc ) 
!------------------------------------------------------------------------------
  USE MaterialModels
  USE Integration
  USE Differentials

  IMPLICIT NONE
  REAL(KIND=dp) :: SheetConductivity(:), ChannelConductivity(:), alphac(:), &
           betac(:), alphas(:), betas(:), IceDensity(:), &
           Snn(:), Ac(:), ng(:), CCt(:), CCw(:), lc(:) 
  INTEGER :: N
  LOGICAL :: Found = .FALSE.

  TYPE(Element_t), POINTER :: Edge
  TYPE(ValueList_t), POINTER :: Material

!------------------------------------------------------------------------------

  SheetConductivity = 0.0_dp
  SheetConductivity(1:n) = ListGetReal(Material, 'Sheet Conductivity', n, Edge % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. )
                   
  ChannelConductivity = 0.0_dp
  ChannelConductivity(1:n) = ListGetReal(Material, 'Channel Conductivity', n, Edge % NodeIndexes, & 
           Found, UnfoundFatal = .TRUE. ) 
                   
  alphac = 0.0_dp
  alphac(1:n) =  ListGetReal( Material, 'Channel Flow Exponent Alpha', n, Edge % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. ) 
                   
  betac = 0.0_dp
  betac(1:n) =  ListGetReal( Material, 'Channel Flow Exponent Beta', n, Edge % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. ) 
  
  alphas = 0.0_dp
  alphas(1:n) =  ListGetReal( Material, 'Sheet Flow Exponent alpha', n, Edge % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. ) 

  betas = 0.0_dp
  betas(1:n) =  ListGetReal( Material, 'Sheet Flow Exponent beta', n, Edge % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. ) 
                   
  IceDensity = 0.0_dp
  IceDensity(1:n) = ListGetReal( Material, 'Density',  N, Edge % NodeIndexes, & 
           Found, UnfoundFatal = .TRUE. ) 

  Snn = 0.0_dp
  Snn(1:n) = ListGetReal( Material, 'Ice Normal Stress',  N, Edge % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. ) 
  IF (ANY(Snn(1:n) < 0.0)) THEN
     WRITE(Message,'(A)')'The Ice Normal Stress (overburden pressure) must be positive'
     CALL FATAL(SolverName, Message)
  END IF

  Ac = 0.0_dp
  Ac(1:n) = ListGetReal( Material, 'Channel Closure Coefficient',  N, Edge % NodeIndexes, & 
           Found, UnfoundFatal = .TRUE. ) 

  ng = 0.0_dp
  ng(1:n) = ListGetReal( Material, 'Glen Exponent',  N, Edge % NodeIndexes, & 
           Found, UnfoundFatal = .TRUE. ) 

  CCt = 0.0_dp
  CCt(1:n) = ListGetReal( Material, 'Pressure Melting Coefficient',  N, Edge % NodeIndexes, & 
           Found, UnfoundFatal = .TRUE. ) 

  CCw = 0.0_dp
  CCw(1:n) = ListGetReal( Material, 'Water Heat Capacity',  N, Edge % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. ) 

  lc = 0.0_dp
  lc(1:n) = ListGetReal( Material, 'Sheet Width Over Channel',  N, Edge % NodeIndexes, &
           Found, UnfoundFatal = .TRUE. ) 
!------------------------------------------------------------------------------
END SUBROUTINE GetParametersChannel
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SheetDischargeCompute( &
       NodalHydPot, hsheet, ks, na, nb, &
       Discharge, Element, n, Nodes, NodeNumber )
!------------------------------------------------------------------------------
    USE MaterialModels
    USE Integration
    USE Differentials

    IMPLICIT NONE
    REAL(KIND=dp) :: NodalHydPot(:), Discharge(:)
    REAL(KIND=dp) :: na, nb, ks, hsheet 
    INTEGER :: n, NodeNumber
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element

    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------

    REAL(KIND=dp) :: dBasisdx(n,3), detJ
    REAL(KIND=dp) :: Basis(n)

    INTEGER :: i, dim

    REAL(KIND=dp) :: u, v, w
    REAL(KIND=dp) :: gradPhi(3), Ngrad 
    LOGICAL :: stat

    !------------------------------------------------------------------------------
    ! the dimension for this solver is given by the dim of the elements
    ! should work for 2D and 3D general problem
    dim = Element % TYPE % DIMENSION

    u = Element % Type % NodeU(NodeNumber)
    v = Element % Type % NodeV(NodeNumber)
    w = Element % Type % NodeW(NodeNumber)

    stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )

    !------------------------------------------------------------------------------
    !      Compute Sheet Discharge  = k h^alphas |grad Phi|^{beta-2} grad Phi
    !------------------------------------------------------------------------------
    gradPhi = 0.0_dp
    DO i=1, dim
       gradPhi(i) = SUM(NodalHydPot(1:n)*dBasisdx(1:n,i))
    END DO
    Ngrad = SQRT(SUM(gradPhi(1:dim)*gradPhi(1:dim)))
    IF (Ngrad < AEPS) Ngrad = AEPS

    DO i=1,dim
       Discharge(i) = -ks * (hsheet**na) * (Ngrad**(nb-2.0_dp)) * gradPhi(i) 
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE SheetDischargeCompute
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE GlaDSCoupledsolver
!------------------------------------------------------------------------------

! *****************************************************************************/
!> Dummy solver just to declare the Sheet Thickness variable
RECURSIVE SUBROUTINE GlaDSsheetThickDummy( Model,Solver,Timestep,TransientSimulation )
!******************************************************************************
!******************************************************************************
     USE Differentials
     USE MaterialModels
     USE DefUtils

!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    External variables
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: Timestep
     
     RETURN 
END SUBROUTINE GlaDSsheetThickDummy
!------------------------------------------------------------------------------

! ******************************************************************************
! *
! *  Authors: Rupert Gladstone
! *  Email:   RupertGladstone972@gmail.com
! *  Web: 
! *
! *  Original Date: 
! *   2022/03/06
! *****************************************************************************
!> Solver GlaDS_GLflux
!> 
!> Take GlaDS standard output and a grounded mask and calculate the total
!> subglacial outflow across the grounding line on grounding line nodes.
!> 
!> The grounded mask is assumed to exist and to have the following properties:
!>  Variable name is GroundedMask
!>  GroundedMask==1 only on fully grounded nodes
!>  GroundedMask==0 only on grounding line nodes
!> 
!> GlaDS variable names can be given as follows (default to these values if
!> not prescribed):
!>  subglac sheet thickness variable = String "Sheet Thickness"
!>  subglac sheet discharge variable = String "Sheet Discharge"
!>  subglac channel flux variable = String "Channel Flux"
!> 
!> In any case, the above variables need to exist!
!> 
!> Limitations:
!> Note that the code currently calculates the flux at the GL based on the 
!> assumption that the subglacial water is always flowing from grounded to ocean 
!> nodes.  If there is inflow from ocean to the subglacial system the cross-GL
!> flux will be overestimated. This could be verified for the channel flux by
!> checking the hydraulic potential at both ends of the edges that are included 
!> in the calculation.  If the grounded node has a higher value than the GL 
!> node then the flow is from grounded to ocean. 
!> Checking that sheet discharge is flowing from grounded to ocean nodes is 
!> more awkward because we'd need to calculate the direction of the normal to 
!> the grounding line.
!> 
SUBROUTINE GlaDS_GLflux( Model,Solver,dt,TransientSimulation )

  USE DefUtils
  USE SolverUtils
  IMPLICIT NONE

  ! intent in
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp)  :: dt
  LOGICAL        :: TransientSimulation

  ! local variables
  TYPE(ValueList_t), POINTER :: SolverParams

  TYPE(Element_t), POINTER   :: Edge
  TYPE(Variable_t), POINTER  :: gmVar, channelVar, sheetThickVar, sheetDisVar
  LOGICAL                    :: GotIt, ValidEdge
  CHARACTER(LEN=MAX_NAME_LEN):: channelVarName, sheetThickVarName, sheetDisVarName, SolverName
  CHARACTER(LEN=MAX_NAME_LEN):: MaskName
  REAL(KIND=dp), POINTER     :: gmVals(:), channelVals(:), sheetThickVals(:), sheetDisVals(:)
  REAL(KIND=dp), POINTER     :: GLfluxVals(:)
  REAL(KIND=dp)              :: volFluxSheet, volFluxChannel, sheetDisMag
  INTEGER, POINTER           :: gmPerm(:), channelPerm(:), sheetThickPerm(:), sheetDisPerm(:)
  INTEGER, POINTER           :: GLfluxPerm(:)
  INTEGER                    :: nn, ee, numNodes

  TYPE(Variable_t), POINTER  :: cglfVar, sglfVar
  REAL(KIND=dp), POINTER     :: cglfVals(:), sglfVals(:)
  INTEGER, POINTER           :: cglfPerm(:), sglfPerm(:)


  SolverName = "GlaDS_GLflux"

  CALL Info(SolverName,'Starting subglacial outflow calculation',Level=4)

  SolverParams => GetSolverParams()
  
  !--------------------------------------------------------------------------------------------
  ! The solver variable will contain the total subglacial outflow on nodes.
  ! Units (assuming Elmer/Ice defaults) m^3/a
  GLfluxVals => Solver % Variable % Values
  GLfluxPerm => Solver % Variable % Perm

  !--------------------------------------------------------------------------------------------
  ! Variables containing the GlaDS sheet thickness and discharge and channel flux

  channelVarName = GetString( SolverParams , 'subglac channel flux variable', GotIt )
  IF (.NOT.GotIt) THEN
     CALL Info(SolverName,'>subglac channel flux variable< not found, assuming >Channel Flux<',Level=4)
     channelVarName = "Channel Flux"
  END IF
  channelVar =>  VariableGet(Model % mesh % Variables,TRIM(channelVarName),UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(channelVar)) &
       CALL FATAL(SolverName,"Variable "//TRIM(channelVarName)//" not found")
  channelPerm => channelVar % Perm
  channelVals => channelVar % Values
  
  sheetThickVarName = GetString( SolverParams , 'subglac sheet thickness variable', GotIt )
  IF (.NOT.GotIt) THEN
     CALL Info(SolverName,'>subglac sheet thickness variable< not found, assuming >sheet thickness<',Level=4)
     sheetThickVarName = "Sheet thickness"
  END IF
  sheetThickVar =>  VariableGet(Model % mesh % Variables,TRIM(sheetThickVarName),UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(sheetThickVar)) &
       CALL FATAL(SolverName,"Variable "//TRIM(sheetThickVarName)//" not found")
  sheetThickPerm => sheetThickVar % Perm
  sheetThickVals => sheetThickVar % Values

  sheetDisVarName = GetString( SolverParams , 'subglac sheet discharge variable', GotIt )
  IF (.NOT.GotIt) THEN
     CALL Info(SolverName,'>subglac sheet discharge variable< not found, assuming >sheet discharge<',Level=4)
     sheetDisVarName = "sheet discharge"
  END IF
  sheetDisVar =>  VariableGet(Model % mesh % Variables,TRIM(sheetDisVarName),UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(sheetDisVar)) &
       CALL FATAL(SolverName,"Variable "//TRIM(sheetDisVarName)//" not found")
  sheetDisPerm => sheetDisVar % Perm
  sheetDisVals => sheetDisVar % Values

  ! grounded mask  
  MaskName = GetString( SolverParams , 'grounded mask variable', GotIt )
  IF (.NOT.GotIt) THEN
     CALL Info(SolverName,'>grounded mask variable< not found, assuming >GroundedMask<',Level=4)
     MaskName = "GroundedMask"
  END IF
  gmVar =>  VariableGet(Model % mesh % Variables,TRIM(MaskName),UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(gmVar)) &
       CALL FATAL(SolverName,"Variable >GroundedMask< not found")
  gmPerm => gmVar % Perm
  gmVals => gmVar % Values

  ! The two variables that will contain the sheet and channel fluxes on the GL are also
  ! hard coded (well, their names anyway).
  sglfVar =>  VariableGet(Model % mesh % Variables,TRIM("Sheet GL flux"),UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(sglfVar)) &
       CALL FATAL(SolverName,"Variable >Sheet GL flux< not found")
  sglfPerm => sglfVar % Perm
  sglfVals => sglfVar % Values
  cglfVar =>  VariableGet(Model % mesh % Variables,TRIM("Channel GL flux"),UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(cglfVar)) &
       CALL FATAL(SolverName,"Variable >Channel GL flux< not found")
  cglfPerm => cglfVar % Perm
  cglfVals => cglfVar % Values
  
  ! Loop over all nodes
  numNodes = Solver % Mesh % Nodes % NumberOfNodes
  DO nn = 1, numNodes

     ! We're interested in nodes where the grounded mask is both defined (non-zero permutation)
     ! and has value set to zero (the grounding line).
     IF (gmPerm(nn).le.0) CYCLE
     IF (gmVals(gmPerm(nn)).eq.0) THEN

        ! Sheet discharge multiplied by sheet thickness gives the volume flux from the sheet.
        ! We're hard conding the assumption that the sheet discharge is always a 2D vector,
        ! which should be safe so long as we always run GlaDS in 2D.
        sheetDisMag = ( sheetDisVals( 2*(sheetDisPerm(nn)-1)+1 )**2.0 +      &
                        sheetDisVals( 2*(sheetDisPerm(nn)-1)+2 )**2.0  )**0.5
        volFluxSheet = sheetThickVals(sheetThickPerm(nn)) * sheetDisMag
        
        volFluxChannel = 0.0

        ! work out channel flux.
        ! loop over all edges...
        DO ee=1, Solver % Mesh % NumberOfEdges 
           Edge => Solver % Mesh % Edges(ee)
           IF (.NOT.ASSOCIATED(Edge)) CYCLE
           ! ...ignoring edges not entirely on the lower surface...
           IF (ANY(gmPerm(Edge % NodeIndexes(1:2)).EQ.0)) CYCLE
           ! ... and check whether the edge contains the current node. If so, check whether the
           ! other node is grounded. If yes, the edge is valid for calculating GL flux.
           ValidEdge = .FALSE.
           IF (Edge % NodeIndexes(1).EQ.nn) THEN
             IF (gmVals(gmPerm(Edge % NodeIndexes(2))).EQ.1) ValidEdge = .TRUE.
           ELSEIF (Edge % NodeIndexes(2).EQ.nn) THEN
              IF (gmVals(gmPerm(Edge % NodeIndexes(1))).EQ.1) ValidEdge = .TRUE.
           END IF
           ! Sum channel flux over valid edges
           IF (ValidEdge) THEN
              IF (Solver % Mesh % ParallelInfo % EdgeInterface(ee)) THEN 
                 ! halve value for edges at partition boundaries because these will be 
                 ! counted twice 
                 volFluxChannel = volFluxChannel + 0.5*channelVals(channelPerm(numNodes+ee))
              ELSE
                 volFluxChannel = volFluxChannel + channelVals(channelPerm(numNodes+ee))
              END IF
           END IF
        END DO
        
        cglfVals(cglfPerm(nn)) = volFluxChannel
        sglfVals(sglfPerm(nn)) = volFluxSheet

     END IF

  END DO

  ! Sum nodal values for nodes that exist on multiple partitions
  CALL ParallelSumVector(Solver % Matrix, cglfVals)

  GLfluxVals = 0.0
  
  DO nn = 1, numNodes
     IF (gmPerm(nn).le.0) CYCLE
!     IF (gmVals(gmPerm(nn)).eq.0) GLfluxVals(GLfluxPerm(nn)) = volFluxSheet + volFluxChannel
     IF (gmVals(gmPerm(nn)).eq.0) GLfluxVals(GLfluxPerm(nn)) =  & 
          cglfVals(cglfPerm(nn)) + sglfVals(sglfPerm(nn))
  END DO
  
  NULLIFY(SolverParams)
  NULLIFY(GLfluxVals)
  NULLIFY(GLfluxPerm)
  NULLIFY(gmVals)
  NULLIFY(gmPerm)
  NULLIFY(sheetDisVals)
  NULLIFY(sheetDisPerm)
  NULLIFY(sheetThickVals)
  NULLIFY(sheetThickPerm)
  NULLIFY(channelVals)
  NULLIFY(channelPerm)

END SUBROUTINE GlaDS_GLflux

! Different ways of calculating a grounded melt rate to pass to GlaDS as a
! volume source.
!
! Notes when using this with a 3D Stokes setup:
! 
! Convert a nodal heat to a melt rate at the lower surface of an ice body.  
! Uses nodal weights (area weighting) to convert the nodal heat to heat per 
! unit area, then convert this to a melt rate.  This solver should run on the 
! lower surface only. The calculated melt rate is in m/a water equivalent (so 
! if you want to use this as a normal velocity condition on the lower surface 
! of the ice body you need to use rho_i to convert to m/a ice equivalent).
!
! Note that the nodal heat could be the residual from the temperate ice solver
! or it could come from the friction load (though this ignores GHF and heat
! conducted into the ice, which may approximately balance each other out...).
!
! [Edit: CalculateNodalWeights gives partition boundary artefacts, but the 
!  forcetostress solver seems to produce weights without these artefacts]
!
! Example .sif parameters:
! 
! Constants:
!  Latent Heat = 334000.0 ! Joules per kg
!
! solver params:
!  variable = GroundedMeltRate
!  Mode = "NodalHeat"
!  heat variable name = String "Friction Load"
!  Weights variable name = String "Friction heating boundary weights"
!
!

! Heat is Mega Joules per year.
  ! We multiply by 10^6 to convert from Mega Joules to Joules.
  !
  
  ! Different modes of operation.
  ! "heat"     - a variable providing nodal heat (e.g. could be residual from temperate ice solver) is used
  !              to calculate the melt rate.  Weights (based on area) are also needed in this case.
  !
  ! MeltRate = Heat / (area * density * latent_heat)
  !
  ! "friction" - a sliding velocity variable is provided and used by this routine to calculate basal shear
  !              stress, which is then used (along with the effective linear sliding coefficient ("ceff",
  !              see SSASolver.F90), to calculate melt based on friction heat.
!
RECURSIVE SUBROUTINE GroundedMelt( Model,Solver,Timestep,TransientSimulation )

  USE DefUtils
  
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  !    External variables
  !------------------------------------------------------------------------------
  TYPE(Model_t)          :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL                :: TransientSimulation
  REAL(KIND=dp)          :: Timestep

  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER  :: SolverParams, Material
  TYPE(Variable_t), POINTER   :: MeltVar, WeightsVar, HeatVar, GHFVar, Ceffvar, UbVar 
  LOGICAL, SAVE               :: FirstTime = .TRUE., UseGHF = .FALSE.
  LOGICAL                     :: Found
  CHARACTER(LEN=MAX_NAME_LEN) :: MyName = 'Grounded Melt solver', HeatVarName, WeightsVarName, GHFvarName
  CHARACTER(LEN=MAX_NAME_LEN) :: MeltMode, CeffVarName, UbVarName
  REAL(KIND=dp)               :: rho_fw ! density of fresh water
  REAL(KIND=dp),PARAMETER     :: threshold = 0.001_dp ! threshold friction melt rate for including GHF in melt calc
  REAL(KIND=dp), POINTER      :: WtVals(:), HeatVals(:), MeltVals(:), GHFVals(:), Ceffvals(:), UbVals(:)
  REAL(KIND=dp)               :: LatHeat, GHFscaleFactor, Ub(1)
  INTEGER, POINTER            :: WtPerm(:), HeatPerm(:), MeltPerm(:), GHFPerm(:), Ceffperm(:), UbPerm(:)
  INTEGER                     :: nn
  

  rho_fw = ListGetConstReal( Model % Constants, 'Fresh Water Density', Found )
  IF (.NOT.Found) CALL FATAL(MyName, 'Constant >Fresh Water Density< not found')
  LatHeat = ListGetConstReal( Model % Constants, 'Latent Heat', Found)
  IF (.NOT.Found) CALL Fatal(MyName, '>Latent Heat< not found in constants')

  MeltVar    => Solver%Variable
  MeltVals   => MeltVar%Values 
  MeltPerm   => MeltVar%Perm

  SolverParams => GetSolverParams()

  MeltMode = GetString(SolverParams,'Melt mode', Found)
  IF(.NOT.Found) CALL Fatal(MyName, '>Melt mode< not found in solver params')
  
  SELECT CASE (MeltMode)

  CASE ("heat")      
    HeatVarName = GetString(SolverParams,'heat variable name', Found)
    IF(.NOT.Found) CALL Fatal(MyName, '>Heat variable name< not found in solver params')
    WeightsVarName = GetString(SolverParams,'Weights variable name', Found)
    IF(.NOT.Found) CALL Fatal(MyName, '>Weights variable name< not found in solver params')

    HeatVar    => VariableGet(Model % Variables, HeatVarName, ThisOnly = .TRUE., UnfoundFatal = .TRUE.)
    HeatVals   => HeatVar%Values 
    HeatPerm   => HeatVar%Perm

    WeightsVar => VariableGet(Model % Variables, WeightsVarName, ThisOnly = .TRUE., UnfoundFatal = .TRUE.)
    WtVals     => WeightsVar%Values 
    WtPerm     => WeightsVar%Perm

  CASE ("friction")
    UbVarName = GetString(SolverParams,'Ub variable name', Found)
    IF (.NOT.Found) UbVarName = "SSAVelocity"
    CeffVarName = GetString(SolverParams,'Ceff variable name', Found)
    IF (.NOT.Found) CeffVarName = "Ceff"

    CeffVar    => VariableGet(Model % Variables, CeffVarName, ThisOnly = .TRUE., UnfoundFatal = .TRUE.)
    CeffVals   => CeffVar%Values 
    CeffPerm   => CeffVar%Perm

    UbVar      => VariableGet(Model % Variables, UbVarName, ThisOnly = .TRUE., UnfoundFatal = .TRUE.)
    UbVals     => UbVar%Values 
    UbPerm     => UbVar%Perm

    IF (UbVar % DOFS .NE. 2) THEN
      CALL Fatal(MyName, 'Expecting Ub variable to be 2D')
    END IF
    !    Material => GetMaterial() ! get sliding velocity from material

  CASE DEFAULT
    CALL Fatal(MyName, 'MeltMode not recognised')

  END SELECT
   
  GHFvarName = GetString(SolverParams,'GHF variable name', Found)
  IF (Found) THEN
    UseGHF = .TRUE.
    GHFscaleFactor = GetConstReal( Model % Constants, 'GHF scale factor', Found)
    IF(.NOT.Found) GHFscaleFactor = 1.0
  ELSE
    UseGHF = .FALSE.
  END IF
       
  IF (UseGHF) THEN
    GHFVar     => VariableGet(Model % Variables, GHFvarName, ThisOnly = .TRUE., UnfoundFatal = .TRUE.)
    GHFVals    => GHFVar%Values 
    GHFPerm    => GHFVar%Perm
  END IF

  LoopAllNodes: DO nn=1,Solver % Mesh % NumberOfNodes

    IF (MeltPerm(nn).GT.0) THEN

      SELECT CASE (MeltMode)
      CASE ("heat")      
        MeltVals(MeltPerm(nn)) = ABS( 1.0e6 * HeatVals(HeatPerm(nn)) ) / ( WtVals(WtPerm(nn)) * rho_fw * LatHeat )
      CASE ("friction")
        Ub = (UbVals(2*(UbPerm(nn)-1)+1)**2 + UbVals(2*(UbPerm(nn)-1)+2)**2)**0.5
!        Ub(1:1) = ListGetReal( Material, 'Sliding Velocity', 1, [nn], Found, UnfoundFatal = .TRUE. )
        MeltVals(MeltPerm(nn)) = (Ub(1)**2 * CeffVals(CeffPerm(nn)) ) / ( rho_fw * LatHeat )
      END SELECT
      
      IF (UseGHF) THEN
        ! Scaled GHF is in Mega Joules per m^2 per year.
        MeltVals(MeltPerm(nn)) = MeltVals(MeltPerm(nn)) + &
             ( GHFVals(GHFPerm(nn))*GHFscaleFactor*1.0e6 ) / ( rho_fw*LatHeat )
      END IF
    END IF

  END DO LoopAllNodes
  
  SELECT CASE(MeltMode)
  CASE("heat")
    NULLIFY(HeatVar, HeatVals, HeatPerm, WeightsVar, WtVals, WtPerm)
  CASE("friction")
    NULLIFY(CeffVar, CeffVals, CeffPerm)
  END SELECT
  NULLIFY(MeltVar, MeltVals, MeltPerm)
  IF (UseGHF) THEN
     NULLIFY(GHFVar, GHFVals, GHFPerm)
  END IF
  
END SUBROUTINE GroundedMelt
