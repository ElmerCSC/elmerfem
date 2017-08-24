! !/*****************************************************************************/
! 
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


! *****************************************************************************/
!>  Save output files for the Channels 
!>  To be used with the GlaDSCoupledSolver   
!>  Solver section also used to declare the Channel Area variable
!>  Need to be executed with the "Exec Solver = After Saving" option
   RECURSIVE SUBROUTINE GlaDSchannelOut( Model,Solver,Timestep,TransientSimulation )
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
     TYPE(Variable_t), POINTER :: HydPotSol, CurrentSol
     TYPE(Variable_t), POINTER :: ThickSol, AreaSol, VSol, WSol, & 
              TimeVar, ZbSol, QcSol
     
     TYPE(ValueList_t), POINTER :: Equation, Material, BC, Constants

     INTEGER :: i, j, k, l, m, n, t, iter, body_id, eq_id, material_id, &
          istat, LocalNodes, bf_id, bc_id,  DIM, NodeSheet, EdgeSheet, NbMoulin, &
          NSDOFs, NonlinearIter, GhostNodes, NonlinearIterMin, it, &
          itOut, BdForder
     INTEGER, ALLOCATABLE :: TableNodeSheet(:), TableMoulin(:)     

     INTEGER, POINTER :: NodeIndexes(:), HydPotPerm(:), ZbPerm(:), &
              ThickPerm(:), VPerm(:), WPerm(:), AreaPerm(:), QcPerm(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, SolverName, HydPotName, &
                OutPutFileName, nit, SheetThicknessName, ZbName

     LOGICAL :: Found, FluxBC, Channels, FileVTK = .FALSE., Storage = .FALSE., &
          AllocationsDone = .FALSE.,  SubroutineVisited = .FALSE., FirstTime = .TRUE., &
          FileASCII = .FALSE., FirstVisit = .TRUE., SaveASCII = .FALSE., &
          SaveVTK = .FALSE., FirstVisit2 = .TRUE. 
          
     LOGICAL, ALLOCATABLE ::  IsGhostNode(:), NoChannel(:), NodalNoChannel(:)

     REAL(KIND=dp) :: NonlinearTol, Relax, &
          SaveRelax, dt, Norm, PrevNorm, S, C
          
     REAL(KIND=dp), POINTER :: HydPot(:), HydPotPrev(:,:), ForceVector(:), PrevSolution(:)
          
     REAL(KIND=dp), POINTER :: ThickSolution(:), VSolution(:), WSolution(:), &
                  AreaSolution(:), AreaPrev(:,:), ZbSolution(:), QcSolution(:)
                         
     REAL(KIND=dp), ALLOCATABLE :: SheetConductivity(:), ChannelConductivity(:)
       
     REAL(KIND=dp), ALLOCATABLE :: ub(:), Snn(:), Ev(:), &
         ng(:), alphas(:), betas(:), betac(:), Phi0(:), Phim(:), Afactor(:), Bfactor(:), &
         MoulinArea(:), MoulinFlux(:), Flux(:)

     REAL(KIND=dp), ALLOCATABLE :: IceDensity(:), Ac(:), alphac(:), CCt(:), &
         CCw(:), lc(:)
     REAL(KIND=dp) :: ChannelArea, EdgeTangent(3), Force, xmean, length, WaterDensity, &
         gravity, Lw, ChannelAreaM1, ChannelAreaM2

     REAL(KIND=dp) :: Np, zb, Xi, Pii, Vc, Nef, x, y, z, dPhidt, Qc 

     SAVE &
          ElementNodes, EdgeNodes,      &
          SheetConductivity,      &
          ChannelConductivity,      &
          IsGhostNode,  OutPutFileName, FileVTK, FileASCII, it,  &
          itOut,    &
          AllocationsDone, FirstTime, FirstVisit, VariableName, SolverName, NonLinearTol, M, &
          Afactor, Bfactor, MoulinArea, MoulinFlux,   &
          ub, Snn, Ev, WaterDensity, gravity, &
          ng, alphas, betas, betac, Phi0, Phim, &
          IceDensity, Ac, alphac, CCt, CCw, lc, Lw, &
          HydPotName, SheetThicknessName, ZbName, &
          NodeSheet, TableNodeSheet, EdgeSheet, NbMoulin, TableMoulin, Flux, &
          NoChannel, NodalNoChannel, BdFOrder, FirstVisit2
!------------------------------------------------------------------------------
!  No sense if not transient
!------------------------------------------------------------------------------
     IF (.Not.TransientSimulation) CALL FATAL( SolverName, 'work only in transient' )

!------------------------------------------------------------------------------
!  This solver is not yet working in parallel, this would require to use 
!  vtu format instead of vtk
!  check here that it is not run in parallel
!------------------------------------------------------------------------------
     IF (ParEnv % PEs > 1) CALL FATAL( SolverName, 'sorry, not yet working in parallel' )

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     SolverName = 'GlaDSchannelOut ('// TRIM(Solver % Variable % Name) // ')'
     VariableName = TRIM(Solver % Variable % Name)
     
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) THEN 
       WRITE(*,*)'Matrix is not associated!!! ',TRIM(SolverName)
       RETURN
     END IF 
     SystemMatrix => Solver % Matrix
     ForceVector => Solver % Matrix % RHS

     PointerToSolver => Solver

     AreaSol => Solver % Variable
     AreaPerm  => AreaSol % Perm
     AreaSolution => AreaSol % Values
     AreaPrev => AreaSol % PrevValues
     
     LocalNodes = COUNT( AreaPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

     DIM = CoordinateSystemDimension()

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
        N = Solver % Mesh % MaxElementNodes
        M = Model % Mesh % NumberOfNodes
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
                SheetConductivity,         &
                ChannelConductivity,         &
                IsGhostNode,              &
                Afactor, Bfactor, MoulinArea, MoulinFlux, &
                ub, Snn, Ev, &
                ng, alphas, betas, betac, Phi0, Phim,  &
                IceDensity, Ac, alphac, CCt, &
                CCw, lc, TableNodeSheet, TableMoulin, Flux, &
                NoChannel, NodalNoChannel)

        END IF                           
        
        ALLOCATE(                                  &
             ElementNodes % x( N ),                &
             ElementNodes % y( N ),                &
             ElementNodes % z( N ),                &
             EdgeNodes % x( N ),                &
             EdgeNodes % y( N ),                &
             EdgeNodes % z( N ),                &
             SheetConductivity( N ),            &
             ChannelConductivity( N ),            &
             IsGhostNode( M ),                     &
             Afactor(N), Bfactor(N), &
             MoulinArea(N), MoulinFlux(N), &
             ub(N), Snn(N), Ev(N), &
             ng(N), alphas(N), betas(N), betac(N), Phi0(N), Phim(N),   &
             IceDensity(N), Ac(N), alphac(N), CCt(N), &
             CCw(N), lc(N), TableNodeSheet(M), TableMoulin(M), Flux(M), &
             NoChannel( M ), NodalNoChannel( N ), &
             STAT=istat )

        IF ( istat /= 0 ) THEN
           CALL FATAL( SolverName, 'Memory allocation error' )
        ELSE
           CALL INFO(SolverName, 'Memory allocation done', level=1 )
        END IF

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

     IF (FirstVisit) THEN 
        FirstVisit = .FALSE.
        ! Read the name of the Hydraulic Potential and Sheet Thickness variables
        ! in the Constant section
        Constants => GetConstants()
        HydPotName = GetString( Constants, &
             'Hydraulic Potential Variable Name', Found )
        IF(.NOT.Found) THEN        
           CALL WARN(SolverName,'Keyword >Hydraulic Potential Variable Name< not found in section Constants')
           CALL WARN(SolverName,'Taking default value >Hydraulic Potential<')
           WRITE(HydPotName,'(A)') 'Hydraulic Potential'
        END IF

        SheetThicknessName = GetString( Constants, &
             'Sheet Thickness Variable Name', Found )
        IF(.NOT.Found) THEN        
           CALL WARN(SolverName,'Keyword >Sheet Thickness Variable Name< not found in section Constants')
           CALL WARN(SolverName,'Taking default value >Sheet Thickness<')
           WRITE(SheetThicknessName,'(A)') 'Sheet Thickness'
        END IF

        ZbName = GetString( Constants, &
             'Bedrock Variable Name', Found )
        IF(.NOT.Found) THEN        
           CALL WARN(SolverName,'Keyword >Bedrock Variable Name< not found in section Constants')
           CALL WARN(SolverName,'Taking default value >zb<')
           WRITE(ZbName,'(A)') 'Zb'
        END IF

        WaterDensity = ListGetConstReal( Constants, 'Water Density', UnFoundFatal=.TRUE. )
        gravity = ListGetConstReal( Constants, 'Gravity Norm', UnFoundFatal=.TRUE. )
        Lw = ListGetConstReal( Constants, 'Latent Heat', UnFoundFatal=.TRUE. )
     END IF

     ! Point on the needed variables
     HydPotSol => VariableGet( Solver % Mesh % Variables, HydPotName, UnfoundFatal = .TRUE. )
     HydPotPerm     => HydPotSol % Perm
     HydPot => HydPotSol % Values
     HydPotPrev => HydPotSol % PrevValues

     ThickSol => VariableGet( Solver % Mesh % Variables, SheetThicknessName, UnfoundFatal = .TRUE. )
     ThickPerm     => ThickSol % Perm
     ThickSolution => ThickSol % Values


     QcSol => VariableGet( Solver % Mesh % Variables, "Channel Flux" )
     IF (ASSOCIATED(QcSol)) THEN
        QcPerm     => QcSol % Perm
        QcSolution => QcSol % Values
     END IF

     ZbSol => VariableGet( Solver % Mesh % Variables, ZbName )
     ! If ZbSol is not associated, we will use mesh coordinate to get the 
     ! bedrock elevation
     IF ( ASSOCIATED( ZbSol ) ) THEN
        ZbPerm     => ZbSol % Perm
        ZbSolution => ZbSol % Values
     END IF

     TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
     IF (FirstTime) THEN
          FirstTime = .FALSE. 
          itOut = 0
          it = 0
          FileVTK = GetLogical( Solver % Values, 'VTK OutPutFile' )
          FileASCII = GetLogical( Solver % Values, 'ASCII OutPutFile' )

          IF (FileVTK.OR.FileASCII) THEN 
             OutPutFileName = GetString( Solver % Values, 'Channels OutPut File Name', Found )
             IF (.NOT.Found) THEN
                WRITE(Message,'(a)') 'Keyword > Channels OutPut File Name < not found'
                CALL FATAL(SolverName, Message)
             END IF

          END IF
                                           
          IF (FileASCII) THEN
             OPEN(10, file = 'results/'//TRIM(OutputFileName)//'.table')
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

                EdgeNodes % x(1:n) = Solver % Mesh % Nodes % x(Edge % NodeIndexes(1:n))
                EdgeNodes % y(1:n) = Solver % Mesh % Nodes % y(Edge % NodeIndexes(1:n))
                EdgeNodes % z(1:n) = Solver % Mesh % Nodes % z(Edge % NodeIndexes(1:n))

                WRITE(10,'(i6,6(x,D15.7))')t, EdgeNodes % x(1), EdgeNodes % y(1), EdgeNodes % z(1),& 
                    EdgeNodes % x(2), EdgeNodes % y(2), EdgeNodes % z(2)
             END DO
             CLOSE(10)
          END IF 
            
          ! Save some information regarding the mesh for VTK output
          ! Number of nodes in the sheet layer: NodeSheet
          ! Number of channels (edges) in the sheet layer: EdgeSheet 
          ! Permutation table for the NodeSheet nodes and global nodes 
          IF (FileVTK) THEN
             NodeSheet = 0 
             TableNodeSheet = 0
             DO i = 1, Model % NumberOfNodes
                IF (HydPotPerm(i)>0) THEN 
                   NodeSheet = NodeSheet + 1 
                   TableNodeSheet(i) = NodeSheet 
                END IF
             END DO
               
             EdgeSheet = 0 
             DO t=1, Solver % Mesh % NumberOfEdges 
                Edge => Solver % Mesh % Edges(t)
                n = Edge % TYPE % NumberOfNodes
                IF (.NOT.ASSOCIATED(Edge)) CYCLE
                IF ((ParEnv % PEs > 1) .AND. &
                   (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE
                IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE
                  EdgeSheet = EdgeSheet + 1
             END DO 
             WRITE(Message,'(a,i0)')'Number of Channels (edges): ', EdgeSheet
             CALL INFO(SolverName, Message, level=2 )
             WRITE(Message,'(a,i0)')'Number of nodes in the sheet layer: ', NodeSheet
             CALL INFO(SolverName, Message, level=2 )

             NbMoulin = 0 
             TableMoulin = 0
             ! Moulins are located on 101 Boundary elements
             DO t=1, Solver % Mesh % NumberOfBoundaryElements
                Element => GetBoundaryElement(t)
                ! IF ( .NOT.ActiveBoundaryElement() ) CYCLE
                IF ((ParEnv % PEs > 1) .AND. &
                   (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE
                n = GetElementNOFNodes()

                IF ( GetElementFamily() > 1 ) CYCLE
                NbMoulin = NbMoulin + 1
                j = Element % NodeIndexes(1)
                TableMoulin(j) = NbMoulin
             END DO
             WRITE(Message,'(a,i0)')'Number of Moulins: ', NbMoulin 
             CALL INFO(SolverName, Message, level=2 )
          END IF
     END IF ! FirstTime
     
     ! If we are here, it means we have to save 
     ! assuming Exec Solver = After Saving is set in the solver section
     ! number of file is incremented by one
     itOut = itOut + 1

     SaveASCII = FileASCII 
     SaveVTK = FileVTK 

     dt = Timestep

     CALL Info( SolverName, ' Channels Output will be saved ', Level=4 )

     !------------------------------------------------------------------------------
     ! Edge element (the channels)
     ! Edge element are created in the SIF file using Element = 'n=0 e:1' 
     !------------------------------------------------------------------------------

     IF (SaveASCII) THEN
           WRITE(nit,'(i4.4)')itOut
           nit = ADJUSTL(nit)
           OPEN(10, file='results/'//TRIM(OutputFileName)//'.'//TRIM(nit) )
           WRITE(10,*)'# time', itOut, it, TimeVar % Values(1)
     END IF
           
     body_id = -1
     NULLIFY(Material)
     DO t=1, Solver % Mesh % NumberOfEdges 
              Edge => Solver % Mesh % Edges(t)
              IF (.NOT.ASSOCIATED(Edge)) CYCLE
              IF ((ParEnv % PEs > 1) .AND. &
                 (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE
              n = Edge % TYPE % NumberOfNodes

              EdgeNodes % x(1:n) = Solver % Mesh % Nodes % x(Edge % NodeIndexes(1:n))
              EdgeNodes % y(1:n) = Solver % Mesh % Nodes % y(Edge % NodeIndexes(1:n))
              EdgeNodes % z(1:n) = Solver % Mesh % Nodes % z(Edge % NodeIndexes(1:n))
              
              ! Work only for 202 elements => n=2
              IF (n/=2) CALL FATAL(SolverName, 'Work only for edge element of type 202')

              ! We keep only the edge which belong in the sheet surface
              ! i.e. Both nodes have Perm > 0
              IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE

              xmean = 0.5*(EdgeNodes % x(1) + EdgeNodes % x(2))

              ! Compute the unit tangent vector of the edge
              EdgeTangent(1) = EdgeNodes % x(2) - EdgeNodes % x(1)
              EdgeTangent(2) = EdgeNodes % y(2) - EdgeNodes % y(1)
              EdgeTangent(3) = EdgeNodes % z(2) - EdgeNodes % z(1)

              length = SQRT(SUM(EdgeTangent*EdgeTangent))
              EdgeTangent = EdgeTangent / length 

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
!
              SheetConductivity = 0.0_dp
              SheetConductivity = ListGetReal(Material, 'Sheet Conductivity', n, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 

              ChannelConductivity = 0.0_dp
              ChannelConductivity = ListGetReal(Material, 'Channel Conductivity', n, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 

              betac = 0.0_dp
              betac(1:N) =  ListGetReal( Material, 'Channel Flow Exponent beta', n, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 
                   
              alphac = 0.0_dp
              alphac(1:N) =  ListGetReal( Material, 'Channel Flow Exponent alpha', n, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 
                   
              alphas = 0.0_dp
              alphas(1:N) =  ListGetReal( Material, 'Sheet Flow Exponent alpha', n, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 
                   
              betas = 0.0_dp
              betas(1:N) =  ListGetReal( Material, 'Sheet Flow Exponent beta', n, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 

              IceDensity = 0.0_dp
              IceDensity(1:N) = ListGetReal( Material, 'Density',  N, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 

              Snn = 0.0_dp
              Snn(1:N) = ListGetReal( Material, 'Ice Normal Stress',  N, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 
              IF (ANY(Snn(1:N) < 0.0)) THEN
                 WRITE(Message,'(A)')'The Ice Normal Stress (overburden &
                       & pressure) must be positive'
                 CALL FATAL(SolverName, Message)
              END IF

              Ac = 0.0_dp
              Ac(1:N) = ListGetReal( Material, 'Channel Closure Coefficient',  N, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 

              ng = 0.0_dp
              ng(1:N) = ListGetReal( Material, 'Glen Exponent',  N, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 

              CCt = 0.0_dp
              CCt(1:N) = ListGetReal( Material, 'Pressure Melting Coefficient',  N, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 

              CCw = 0.0_dp
              CCw(1:N) = ListGetReal( Material, 'Water Heat Capacity',  N, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 

              lc = 0.0_dp
              lc(1:N) = ListGetReal( Material, 'Sheet Width Over Channel',  N, Edge % NodeIndexes, &
                   Found, UnfoundFatal = .TRUE. ) 

              Phi0 = 0.0_dp
              Phim = 0.0_dp
              Afactor = 0.0_dp
              Bfactor = 0.0_dp
              DO i=1,N
                 IF ( ASSOCIATED( ZbSol )) THEN
                    zb = ZbSolution(ZbPerm(Edge % NodeIndexes(i)))
                 ELSE 
                 ! TODO - I don't see here how to check if y or z should be
                 ! used... Will not work for a flow line problem therefore...
                    zb = Solver % Mesh % Nodes % z(Edge % NodeIndexes(i))
                 END IF
                 Phim(i) = gravity*WaterDensity*zb
                 Phi0(i) = Snn(i) + Phim(i)
                 Afactor(i) = CCt(i) * CCw(i) * WaterDensity
                 Bfactor(i) = 1.0/(IceDensity(i)*Lw)
              END DO

              M = Model % Mesh % NumberOfNodes
              ChannelArea = AreaPrev(AreaPerm(M+t),1)
              ChannelAreaM1 = AreaPrev(AreaPerm(M+t),2)
              ChannelAreaM2 = AreaPrev(AreaPerm(M+t),3)

              IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                 CALL GetChannelAreaForce( Force, Qc, Xi, Pii, Vc, Nef, &
                      ChannelArea, &
                      HydPot(HydPotPerm(Edge % NodeIndexes(1:n))), &
                      ThickSolution(ThickPerm(Edge % NodeIndexes(1:n))), &
                      alphac, betac, ChannelConductivity, Phi0, Phim, Ac, lc, ng, &
                      SheetConductivity, alphas, betas, Afactor, Bfactor, &
                      EdgeTangent, Edge, n, EdgeNodes )
              ELSE
                 WRITE(Message,'(A)')' Work only for cartesian coordinate'
                 CALL FATAL( SolverName, Message)
              END IF              
              IF (ASSOCIATED(QcSol)) QcSolution(QcPerm(M+t)) = Qc
           
          ! We check if we are in a boundary where we want No channel
          IF (ALL(NoChannel(Edge % NodeIndexes(1:n)))) THEN 
             AreaSolution(AreaPerm(M+t)) = 0.0_dp 
             Xi = 0.0_dp
             Pii = 0.0_dp
             Vc = 0.0_dp
          END IF


          IF (SaveASCII) WRITE(10,'(i9,10(x,D15.7))')t, xmean, AreaSolution(AreaPerm(M+t)), Qc,  &
              HydPot(HydPotPerm(Edge % NodeIndexes(1))), HydPot(HydPotPerm(Edge % NodeIndexes(2))), &
              Bfactor(1)*Xi, Bfactor(1)*Pii, Vc, Nef, length
    END DO ! loop over the Edges

    WRITE(Message,'(A,D15.7)')' Maximum Channel Area: ', &
        MAXVAL(AreaSolution(AreaPerm(M:M+Solver%Mesh%NumberOfEdges))) 
    CALL INFO(SolverName, Message, level=2 )

    IF (SaveASCII) THEN
             CLOSE(10)
    END IF

    ! Save results in VTK Format
    IF (SaveVTK) THEN
         WRITE(nit,'(i4.4)')itOut
         nit = ADJUSTL(nit)
         OPEN(10, file='results/'//TRIM(OutPutFileName)//TRIM(nit)//'.vtk' )
         WRITE ( 10, '("# vtk DataFile Version 3.0")' )
         WRITE ( 10, '("ElmerSolver output; started at ", A)' ) TRIM(FormatDate() )
         WRITE ( 10, '("ASCII")' )
         WRITE ( 10, '("DATASET POLYDATA")' )
         WRITE ( 10, '("POINTS ", I0, " double")' ) NodeSheet 

         DO i = 1, Model % NumberOfNodes
           IF (TableNodeSheet(i)>0) THEN 
              x = Solver % Mesh % Nodes % x(i)
              y = Solver % Mesh % Nodes % y(i)
              z = Solver % Mesh % Nodes % z(i)
              WRITE( 10 ,'(3ES16.7E3)' ) x, y, z
           END IF
         END DO
       
         WRITE ( 10, '("LINES ", I0, " ", I0)' ) EdgeSheet, 3*EdgeSheet 
         DO t=1, Solver % Mesh % NumberOfEdges 
            Edge => Solver % Mesh % Edges(t)
            IF (.NOT.ASSOCIATED(Edge)) CYCLE
            IF ((ParEnv % PEs > 1) .AND. &
               (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE
            n = Edge % TYPE % NumberOfNodes
            IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE
            
            WRITE( 10, '("2 ", I0, " ", I0)')TableNodeSheet(Edge % NodeIndexes(1))-1, &
                 TableNodeSheet(Edge % NodeIndexes(2))-1
         END DO

         WRITE( 10, '("CELL_DATA ", I0)') EdgeSheet
         WRITE( 10, '("SCALARS Channel_Area double ")')
         WRITE( 10, '("LOOKUP_TABLE defaut ")')
         
         ! Loops to get the Channels Area variable
         DO t=1, Solver % Mesh % NumberOfEdges 
            Edge => Solver % Mesh % Edges(t)
            IF (.NOT.ASSOCIATED(Edge)) CYCLE
            IF ((ParEnv % PEs > 1) .AND. &
               (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE
            n = Edge % TYPE % NumberOfNodes
            IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE
            
            M = Model % Mesh % NumberOfNodes
            WRITE( 10, '(E16.7)')AreaSolution(AreaPerm(M+t))
         END DO

         ! Loops to get the Channels Flux variable
         IF (ASSOCIATED(QcSol)) THEN
            WRITE( 10, '("SCALARS Channel_Flux double ")')
            WRITE( 10, '("LOOKUP_TABLE defaut ")')
            DO t=1, Solver % Mesh % NumberOfEdges 
               Edge => Solver % Mesh % Edges(t)
               IF (.NOT.ASSOCIATED(Edge)) CYCLE
               IF ((ParEnv % PEs > 1) .AND. &
                  (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE
               n = Edge % TYPE % NumberOfNodes
               IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE
            
               M = Model % Mesh % NumberOfNodes
               WRITE( 10, '(E16.7)')QcSolution(QcPerm(M+t))
            END DO
         END IF

         IF (NbMoulin>0) THEN 
            WRITE( 10, '("POINT_DATA ", I0)') NodeSheet
            WRITE( 10, '("SCALARS Moulin_Input double ")')
            WRITE( 10, '("LOOKUP_TABLE defaut ")')

            ! Loops to get the Moulins input
            Flux = 0.0_dp
            DO t=1, Solver % Mesh % NumberOfBoundaryElements
               Element => GetBoundaryElement(t)
               ! IF ( .NOT.ActiveBoundaryElement() ) CYCLE
               n = GetElementNOFNodes()

               IF ( GetElementFamily() > 1 ) CYCLE
               BC => GetBC()
               bc_id = GetBCId( Element )
               CALL GetElementNodes( ElementNodes )
               IF ( ASSOCIATED( BC ) ) THEN            
                  Storage = .False.
                  Storage =  GetLogical(BC,'Moulin Storage', Found)
                  MoulinArea = 0.0_dp
                  IF (Storage) THEN
                     MoulinArea(1:N) = ListGetReal( BC, 'Moulin Area',  N, Element % NodeIndexes, &
                          Found, UnfoundFatal = .TRUE. ) 
                  END IF

                  ! Is there surface input
                  MoulinFlux = 0.0_dp
                  MoulinFlux(1:N) = ListGetReal( BC, 'Moulin Flux',  N, Element % NodeIndexes, &
                       Found, UnfoundFatal = .TRUE. ) 
                  j = Element % NodeIndexes(1)
                  dphidt = (HydPot(HydPotPerm(j)) -  HydPotPrev(HydPotPerm(j),1))/dt 
                  Flux(j) = MoulinFlux(1) +  MoulinArea(1)/(WaterDensity*gravity)*dPhidt 
               END IF
            END DO
            DO i = 1, Model % NumberOfNodes
               IF (TableNodeSheet(i)>0) THEN 
                  WRITE( 10, '(E16.7)')Flux(i)
               END IF
            END DO
         END IF 

         CLOSE(10)
    END IF ! Output VTK

    SubroutineVisited = .TRUE.

CONTAINS    

!------------------------------------------------------------------------------
SUBROUTINE GetChannelAreaForce(Force, Qcc, Xi, Pii, Vc, Nef, &
      CArea, NodalHydPot, NodalH, &
      NodalAlphac, NodalBetac, NodalKc, NodalPhi0, NodalPhim, NodalAc, Nodallc, Nodalng, &
      NodalKs, NodalAlphas, NodalBetas, NodalAfactor, NodalBfactor, &
      Tangent, Element, n, Nodes )
!------------------------------------------------------------------------------
  USE MaterialModels
  USE Integration
  USE Differentials

  IMPLICIT NONE

            
     REAL(KIND=dp) :: Force, Xi, Pii, Vc, Nef, Qcc
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
       Ks = Ks * Ngrad**(nbs-2.0) 

       Kc = SUM( NodalKc(1:n) * Basis(1:n))
       Kc = Kc * CArea**nac 
       Kc = Kc * Ngrad**(nbc-2.0)  

       PhiG = SUM(NodalHydPot(1:n)*Basis(1:n))
       Phi0 = SUM(NodalPhi0(1:n)*Basis(1:n))

       ng = SUM(NodalNg(1:n)*Basis(1:n))
       Vc = CArea*SUM(NodalAc(1:n)*Basis(1:n))
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

       Qcc = - Kc*GradPhi
       Xi = ABS(Qcc*GradPhi)+ABS(lc*qc*GradPhi)
       Pii = -Afactor*(Qcc + Ffactor)*dPw
! Save Qcc as a positive value
       Qcc = Abs(Qcc)
       Nef = Phi0 - PhiG

       Force = Bfactor*(Xi - Pii) - Vc 

!------------------------------------------------------------------------------
END SUBROUTINE GetChannelAreaForce
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE GlaDSchannelOut
!------------------------------------------------------------------------------
