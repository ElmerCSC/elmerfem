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

     INTEGER :: i, j, k, l, m, n, t, iter, body_id, eq_id, material_id, &
          istat, LocalNodes,bf_id, bc_id,  DIM, dimSheet, iterC, &
          NSDOFs, NonlinearIter, GhostNodes, NonlinearIterMin, Ne, BDForder, &
          CoupledIter

     TYPE(Variable_t), POINTER :: HydPotSol
     TYPE(Variable_t), POINTER :: ThickSol, AreaSol, VSol, WSol, NSol,  &
              PwSol, ZbSol, qSol, hstoreSol, QcSol, QmSol

     INTEGER, POINTER :: NodeIndexes(:), HydPotPerm(:), PwPerm(:), ZbPerm(:), &
             ThickPerm(:), VPerm(:), WPerm(:), NPerm(:), AreaPerm(:), & 
             qPerm(:), hstorePerm(:), QcPerm(:), QmPerm(:)

     REAL(KIND=dp), POINTER :: HydPot(:), HydPotPrev(:,:), ForceVector(:)
     REAL(KIND=dp), POINTER :: ThickSolution(:), ThickPrev(:,:), VSolution(:), WSolution(:), &
            NSolution(:), PwSolution(:), AreaSolution(:), AreaPrev(:,:), ZbSolution(:), &
            qSolution(:), hstoreSolution(:), QcSolution(:), QmSolution(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, SolverName
     CHARACTER(LEN=MAX_NAME_LEN) :: SheetThicknessName, ChannelAreaName, ZbName
     CHARACTER(LEN=MAX_NAME_LEN) :: methodSheet, methodChannels 

     LOGICAL :: Found, FluxBC, Channels, Storage, FirstTime = .TRUE., &
          AllocationsDone = .FALSE.,  SubroutineVisited = .FALSE., &
          meltChannels = .TRUE., NeglectH = .TRUE. 
     LOGICAL, ALLOCATABLE ::  IsGhostNode(:), NoChannel(:), NodalNoChannel(:)

     REAL(KIND=dp) :: NonlinearTol, dt, CumulativeTime, RelativeChange, &
          Norm, PrevNorm, S, C, Qc
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

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, at0
#else
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif

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
          Vvar, ublr, hr2, Refq

      
     totst = 0.0_dp
     totat = 0.0_dp

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     VariableName = TRIM(Solver % Variable % Name)
     SolverName = 'GlaDSCoupledsolver ('// VariableName // ')'

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

     DIM = CoordinateSystemDimension()

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
        N = Solver % Mesh % MaxElementNodes
        M = Model % Mesh % NumberOfNodes
        Ne = Solver % Mesh % NumberOfEdges
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
                Refq )

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
             STAT=istat )

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
              IF (ParEnv % myPe .EQ. Element % partIndex) CYCLE
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
        DO t=1, Solver % Mesh % NumberOfBoundaryElements
           ! get element information
           Element => GetBoundaryElement(t)
           !IF ( .NOT.ActiveBoundaryElement() ) CYCLE
           IF ((ParEnv % PEs > 1) .AND. &
            (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE

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
     IF (FirstTime) THEN
        FirstTime = .FALSE.
        Constants => GetConstants()

        WaterDensity = ListGetConstReal( Constants, 'Water Density', UnFoundFatal=.TRUE. )
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
        IF(.NOT.Found) THEN        
           CALL WARN(SolverName,'Keyword >Bedrock Variable Name< not found in section Constants')
           CALL WARN(SolverName,'Taking default value >zb<')
           WRITE(ZbName,'(A)') 'Zb'
        END IF

        ! TODO : implement higher order BDF method
        BDForder = GetInteger(GetSimulation(),'BDF Order', Found)
        IF (.NOT.Found) BDForder = 1
        IF (BDForder .NE. 1) THEN
           WRITE(Message,'(a)') 'Only working for BDF = 1' 
           CALL FATAL(SolverName, Message)
        END IF
     END IF ! FirstTime

     SolverParams => GetSolverParams()

     NeglectH = GetLogical( SolverParams,'Neglect Sheet Thickness in Potential', Found )
     IF ( .NOT.Found ) THEN
        CALL FATAL(SolverName, 'No >Neglect Sheet Thickness in Potential< found')
     END IF

     methodSheet = GetString( SolverParams, 'Sheet Integration Method', Found )
     IF(.NOT.Found) THEN        
        CALL FATAL(SolverName, 'No >Sheet Integration Methods< found')
     ELSE
        IF ((methodSheet.NE.'explicit').AND.(methodSheet.NE.'implicit').AND.(methodSheet.NE.'crank-nicolson')) THEN
           CALL FATAL(SolverName, 'Sheet Integration Method: Implicit, Explicit or Crank-Nicolson')
        END IF 
     END IF

     Channels = GetLogical( SolverParams,'Activate Channels', Found )
     IF (.NOT. Found) Channels = .FALSE.

     IF (Channels) THEN
        meltChannels = GetLogical( SolverParams,'Activate Melt From Channels', Found )
        IF ( .NOT.Found ) THEN
           CALL FATAL(SolverName, 'No >Activate Melt From Channels< found')
        END IF

        methodChannels = GetString( SolverParams, 'Channels Integration Method', Found )
        IF(.NOT.Found) THEN        
           CALL FATAL(SolverName, 'No >Channels Integration Methods< found')
        ELSE
           IF ((methodChannels.NE.'explicit').AND.(methodChannels.NE.'implicit').AND.(methodChannels.NE.'crank-nicolson')) THEN
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

     ThickSol => VariableGet( Solver % Mesh % Variables, SheetThicknessName, UnfoundFatal = .TRUE. )
     ThickPerm     => ThickSol % Perm
     ThickSolution => ThickSol % Values
     ThickPrev => ThickSol % PrevValues  

     IF (Channels) THEN
        AreaSol => VariableGet( Solver % Mesh % Variables, ChannelAreaName, UnfoundFatal = .TRUE. )
        AreaPerm     => AreaSol % Perm
        AreaSolution => AreaSol % Values
        AreaPrev => AreaSol % PrevValues
        
        ! flux in the channels (for output only) - edge type variable
        QcSol => VariableGet( Solver % Mesh % Variables, "Channel Flux" )
        IF ( ASSOCIATED( QcSol ) ) THEN
           QcPerm     => QcSol % Perm
           QcSolution => QcSol % Values
        END IF
     END IF

     ! discharge out of the moulins (for output only)
     QmSol => VariableGet( Solver % Mesh % Variables, "Flux from Moulins" )
     IF ( ASSOCIATED( QmSol ) ) THEN
        QmPerm     => QmSol % Perm
        QmSolution => QmSol % Values
     END IF

     ZbSol => VariableGet( Solver % Mesh % Variables, ZbName )
     IF ( ASSOCIATED( ZbSol ) ) THEN
        ZbPerm     => ZbSol % Perm
        ZbSolution => ZbSol % Values
     END IF

     dt = Timestep

!------------------------------------------------------------------------------
!    Read the other variables needed                  
!------------------------------------------------------------------------------
     VSol => VariableGet( Solver % Mesh % Variables, 'Vclose' )
     IF ( ASSOCIATED( VSol ) ) THEN
          VPerm     => VSol % Perm
          VSolution => VSol % Values
     END IF

     WSol => VariableGet( Solver % Mesh % Variables, 'Wopen' )
     IF ( ASSOCIATED( WSol ) ) THEN
          WPerm     => WSol % Perm
          WSolution => WSol % Values
     END IF

     NSol => VariableGet( Solver % Mesh % Variables, 'Effective Pressure' )
     IF ( ASSOCIATED( NSol ) ) THEN
          NPerm     => NSol % Perm
          NSolution => NSol % Values
     END IF

     PwSol => VariableGet( Solver % Mesh % Variables, 'Water Pressure' )
     IF ( ASSOCIATED( PwSol ) ) THEN
          PwPerm     => PwSol % Perm
          PwSolution => PwSol % Values
     END IF

     qSol => VariableGet( Solver % Mesh % Variables, 'Sheet Discharge' )
     IF ( ASSOCIATED( qSol ) ) THEN
          qPerm     => qSol % Perm
          qSolution => qSol % Values
     END IF

     hstoreSol => VariableGet( Solver % Mesh % Variables, 'Sheet Storage' )
     IF ( ASSOCIATED( hstoreSol ) ) THEN
          hstorePerm     => hstoreSol % Perm
          hstoreSolution => hstoreSol % Values
     END IF

!------------------------------------------------------------------------------
!   Loop for the coupling of the three equations
!------------------------------------------------------------------------------
    ! check on the coupled Convergence is done on the potential solution only
    M = Model % Mesh % NumberOfNodes
    IF (ParEnv % PEs > 1) THEN
       PrevCoupledNorm = ParallelNorm(M,HydPot(HydPotPerm(1:M))) 
    ELSE            
       PrevCoupledNorm = SQRT(SUM(HydPot(HydPotPerm(1:M))*HydPot(HydPotPerm(1:M))))/M
    END IF

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
              
           M = Model % Mesh % NumberOfNodes
           Ne = Solver % Mesh % NumberOfEdges
            
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
              IF (ParEnv % myPe .NE. Element % partIndex) CYCLE

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
                       zb = Solver % Mesh % Nodes % y(j)
                    ELSE
                       zb = Solver % Mesh % Nodes % z(j)
                    END IF 
                 END IF
                 CT(i) = Ev(i) /( WaterDensity * gravity)
                 Wopen(i) = MAX(ub(i) / lr(i) * (hr(i) - ThickSolution(k)), 0.0) 
                 Phi0(i) = Snn(i) + gravity*WaterDensity*zb
                 IF (.Not.NeglectH) THEN 
                    Phi0(i) = Phi0(i) + gravity*WaterDensity*ThickSolution(k)
                 END IF
              END DO

              !------------------------------------------------------------------------------
              ! Add body forces
              !------------------------------------------------------------------------------
              LOAD = 0.0_dp
              BodyForce => GetBodyForce()
              IF ( ASSOCIATED( BodyForce ) ) THEN
                 bf_id = GetBodyForceId()
                 LOAD(1:N) = LOAD(1:N) + &
                      GetReal( BodyForce, TRIM(Solver % Variable % Name) // ' Volume Source', Found )
              END IF
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
           DO t=1, Solver % Mesh % NumberOfEdges 
              
              Edge => Solver % Mesh % Edges(t)
              IF (.NOT.ASSOCIATED(Edge)) CYCLE
              IF ((ParEnv % PEs > 1) .AND. &
                (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE
              n = Edge % TYPE % NumberOfNodes

              ! Work only for 202 elements => n=2
              IF (n/=2) CALL FATAL(SolverName, 'Work only for edge element of type 202')

              ! We keep only the edge which belong in the sheet surface
              ! i.e. Both nodes have Perm > 0
              IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE

              ! We check if we are in a boundary where we want No channel
              IF (ALL(NoChannel(Edge % NodeIndexes(1:n)))) CYCLE

              EdgeNodes % x(1:n) = Solver % Mesh % Nodes % x(Edge % NodeIndexes(1:n))
              EdgeNodes % y(1:n) = Solver % Mesh % Nodes % y(Edge % NodeIndexes(1:n))
              EdgeNodes % z(1:n) = Solver % Mesh % Nodes % z(Edge % NodeIndexes(1:n))

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
                       zb = Solver % Mesh % Nodes % y(Edge % NodeIndexes(i))
                    ELSE
                       zb = Solver % Mesh % Nodes % z(Edge % NodeIndexes(i))
                    END IF
                 END If
                 Phim(i) = gravity*WaterDensity*zb
                 Phi0(i) = Snn(i) + Phim(i) 
                 IF (.Not.NeglectH) THEN
                    k = ThickPerm(Edge % NodeIndexes(i))
                    Phi0(i) = Phi0(i) + gravity*WaterDensity*ThickSolution(k)
                 END IF
                 Afactor(i) = CCt(i) * CCw(i) * WaterDensity
                 Bfactor(i) = 1.0/(Lw * IceDensity(i)) 
                 IF (meltChannels) Bfactor(i) = Bfactor(i) - 1.0/(Lw * WaterDensity)
              END DO

! The variable Channel Area is defined on the edges only
              M = Model % Mesh % NumberOfNodes
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
           DO t=1, Solver % Mesh % NumberOfBoundaryElements
              ! get element information
              Element => GetBoundaryElement(t)
              ! Don't test this as BC on the side are not active for this solver
              !IF ( .NOT.ActiveBoundaryElement() ) CYCLE

              ! cycle halo elements
              !--------------------
              IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
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

                  Storage = .False.
                  Storage =  GetLogical(BC,'Moulin Storage', Found)
                  IF (Storage) THEN
                    MoulinArea = 0.0_dp
                    MoulinArea(1:N) = ListGetReal( BC, 'Moulin Area',  N, Element % NodeIndexes, Found, &
                         UnfoundFatal = .TRUE. )
                    ! MASS is a scalar here
                    MASS(1,1) = MoulinArea(1)/(WaterDensity*gravity)
                  END IF

                  ! Is there surface input
                  MoulinFlux = 0.0_dp
                  MoulinFlux(1:N) = ListGetReal( BC, 'Moulin Flux',  N, Element % NodeIndexes, Found, &
                         UnfoundFatal = .TRUE. )
                  FORCE(1) = MoulinFlux(1)

                  ! If variable exist, update the Flux from Moulins variable 
                  IF (ASSOCIATED( QmSol )) THEN
                     j = Element % NodeIndexes(1)
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
                  FluxBC = .FALSE.
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

           DO t=1,Solver % NumberOfActiveElements
              Element => GetActiveElement(t,Solver)
              IF (ParEnv % myPe .NE. Element % partIndex) CYCLE

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
                       zb = Solver % Mesh % Nodes % y(j)
                    ELSE
                       zb = Solver % Mesh % Nodes % z(j)
                    END IF 
                 END IF
                 CT(i) = Ev(i) /( WaterDensity * gravity)
                 Phi0(i) = Snn(i) + gravity*WaterDensity*zb
                 Np = Phi0(i) - HydPot(HydPotPerm(j))  
                 IF (.Not.NeglectH) THEN
                    Np = Np + WaterDensity*gravity*ThickSolution(k)
                 END IF
                 pw = HydPot(HydPotPerm(j)) - gravity*WaterDensity*zb
                 Vvar(j) = Ar(i)*ABS(Np)**(ng(i)-1.0)*Np
                 Wo = MAX(ub(i) / lr(i) * (hr(i) - ThickSolution(k)), 0.0) 
                 he = Ev(i)*(HydPot(HydPotPerm(j))/(WaterDensity*gravity)-zb)
                 ublr(j) = ub(i)/lr(i)
                 hr2(j) = hr(i) 

                 ! Save this for output only
                 IF (ASSOCIATED(WSol)) WSolution(WPerm(j)) = Wo
                 IF (ASSOCIATED(NSol)) NSolution(NPerm(j)) = Np
                 IF (ASSOCIATED(PwSol)) PwSolution(PwPerm(j)) = pw
                 IF (ASSOCIATED(hstoreSol)) hstoreSolution(hstorePerm(j)) = he 

              END DO
           END DO     !  Bulk elements

           ! Loop over all nodes to update ThickSolution
           DO j = 1, Model % Mesh % NumberOfNodes
              k = ThickPerm(j)
              IF (k==0) CYCLE
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
        t = Solver % Mesh % NumberOfEdges 
        M = Model % Mesh % NumberOfNodes

        IF (ParEnv % PEs > 1) THEN
           PrevNorm = ParallelNorm(t,AreaSolution(AreaPerm(M+1:M+t))) 
        ELSE            
           PrevNorm = SQRT(SUM(AreaSolution(AreaPerm(M+1:M+t))*AreaSolution(AreaPerm(M+1:M+t))))/t  
        END IF

        DO iter = 1, NonlinearIter
              DO t=1, Solver % Mesh % NumberOfEdges 
                 Edge => Solver % Mesh % Edges(t)
                 IF (.NOT.ASSOCIATED(Edge)) CYCLE
                 IF ((ParEnv % PEs > 1) .AND. &
                   (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE
                 n = Edge % TYPE % NumberOfNodes
                 IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE
                 IF (ALL(NoChannel(Edge % NodeIndexes(1:n)))) CYCLE

                 EdgeNodes % x(1:n) = Solver % Mesh % Nodes % x(Edge % NodeIndexes(1:n))
                 EdgeNodes % y(1:n) = Solver % Mesh % Nodes % y(Edge % NodeIndexes(1:n))
                 EdgeNodes % z(1:n) = Solver % Mesh % Nodes % z(Edge % NodeIndexes(1:n))

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
                          zb = Solver % Mesh % Nodes % y(Edge % NodeIndexes(i))
                       ELSE
                          zb = Solver % Mesh % Nodes % z(Edge % NodeIndexes(i))
                       END IF
                    END If
                    Phim(i) = gravity*WaterDensity*zb
                    Phi0(i) = Snn(i) + Phim(i) 
                    IF (.Not.NeglectH) THEN
                       k = ThickPerm(Edge % NodeIndexes(i))
                       Phi0(i) = Phi0(i) + gravity*WaterDensity*ThickSolution(k)
                    END IF 
                    Afactor(i) = CCt(i) * CCw(i) * WaterDensity
                    Bfactor(i) = 1.0/(Lw * IceDensity(i)) 
                 END DO

                 M = Model % Mesh % NumberOfNodes
                 ChannelArea = AreaSolution(AreaPerm(M+t))

! Compute the force term to evolve the channels area
! Equation of the form dS/dt = S x ALPHA + BETA
                 IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                    CALL GetEvolveChannel( ALPHA, BETA, Qc, ChannelArea, &
                         HydPot(HydPotPerm(Edge % NodeIndexes(1:n))), &
                         ThickSolution(ThickPerm(Edge % NodeIndexes(1:n))), &
                         alphac, betac, ChannelConductivity, Phi0, Phim, Ac, lc, ng, &
                         SheetConductivity, alphas, betas, Afactor, Bfactor, &
                         EdgeTangent, Edge, n, EdgeNodes )
                 ELSE
                    WRITE(Message,'(A)')' Work only for cartesian coordinate'
                    CALL FATAL( SolverName, Message)
                 END IF

                 k = AreaPerm(M+t)
                 SELECT CASE(methodSheet)
                 CASE('implicit') 
                    AreaSolution(k) = (AreaPrev(k,1) + dt*BETA)/(1.0_dp - dt*ALPHA)
                 CASE('explicit')
                    AreaSolution(k) = AreaPrev(k,1)*(1.0_dp + dt*ALPHA) + dt*BETA 
                 CASE('crank-nicolson')
                    AreaSolution(k) = (AreaPrev(k,1)*(1.0_dp + 0.5*ALPHA*dt) + dt*BETA)/(1.0_dp - 0.5*dt*ALPHA)
                 END SELECT 

                 ! Save Qc if variable exists
                 IF (ASSOCIATED(QcSol)) THEN
                    QcSolution(QcPerm(M+t)) = Qc
                 END IF
              END DO

           t = Solver % Mesh % NumberOfEdges 
           IF (ParEnv % PEs > 1) THEN
              Norm = ParallelNorm(t,AreaSolution(AreaPerm(M+1:M+t))) 
           ELSE            
              Norm = SQRT(SUM(AreaSolution(AreaPerm(M+1:M+t))*AreaSolution(AreaPerm(M+1:M+t))))/t  
           END IF

           IF ( PrevNorm + Norm /= 0.0d0 ) THEN
              RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
           ELSE
              RelativeChange = 0.0d0
           END IF

           PrevNorm = Norm 

           WRITE( Message, * ) 'S (NRM,RELC) : ',iter, Norm, RelativeChange
           CALL Info( SolverName, Message, Level=3 )
           WRITE( Message, * ) 'Max S :', MAXVAL(AreaSolution(AreaPerm(M+1:M+t))),MAXLOC(AreaSolution)  
           CALL Info( SolverName, Message, Level=3 )
           WRITE( Message, * ) 'Min S :', MINVAL(AreaSolution(AreaPerm(M+1:M+t))),MINLOC(AreaSolution)   


    !      !----------------------
    !      ! check for convergence
    !      !----------------------
           IF ( RelativeChange < NonlinearTol .AND. iter > NonlinearIterMin ) EXIT 
        END DO ! of the nonlinear iteration

        ! Make sure Area > 0
        AreaSolution(AreaPerm(M+1:M+t)) = MAX(AreaSolution(AreaPerm(M+1:M+t)),0.0_dp)


     END IF  ! If Channels

      !   Check for convergence                           
      IF (ParEnv % PEs > 1) THEN
         CoupledNorm = ParallelNorm(M,HydPot(HydPotPerm(1:M))) 
      ELSE            
         CoupledNorm = SQRT(SUM(HydPot(HydPotPerm(1:M))*HydPot(HydPotPerm(1:M))))/M
      END IF
      
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
   M = Model % Mesh % NumberOfNodes
   IF (ASSOCIATED(VSol)) VSolution(VPerm(1:M)) = Vvar(1:M)

! Output the sheet discharge (dimension = dimSheet)
   IF (ASSOCIATED(qSol)) THEN
      Refq = 0.0_dp
      qSolution = 0.0_dp

      ! Loop over all elements are we need to compute grad(Phi)
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t,Solver)
         IF (ParEnv % myPe .NE. Element % partIndex) CYCLE

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
      WHERE(Refq > 0.0_dp)
         qSolution = qSolution / Refq
      END WHERE

   END IF

   SubroutineVisited = .TRUE.

CONTAINS    

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

    REAL(KIND=dp) :: CT, C2, Vfactor, Phi0, PhiG

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
      Tangent, Element, n, Nodes )
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

     REAL(KIND=dp) :: Phi0, PhiG, Afactor, Bfactor, GradPhim, dPw, Ffactor

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
       Kc = Kc * MAX(CArea,0.0)**nac
       Kc = Kc * Ngrad**(nbc-2.0_dp)  

       PhiG = SUM(NodalHydPot(1:n)*Basis(1:n))
       Phi0 = SUM(NodalPhi0(1:n)*Basis(1:n))

       ng = SUM(NodalNg(1:n)*Basis(1:n))
       Vc = SUM(NodalAc(1:n)*Basis(1:n))
       Vc = Vc*ABS(Phi0-PhiG)**(ng-1.0_dp)
       Vc = Vc*(Phi0-PhiG)

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
       Qcc = ABS(Kc*GradPhi)

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
  ub(1:N) = ListGetReal( Material, 'Sliding Velocity',  N, Element % NodeIndexes, &
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
