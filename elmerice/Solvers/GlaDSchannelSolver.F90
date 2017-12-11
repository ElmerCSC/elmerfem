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
     TYPE(Nodes_t) :: ElementNodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Edge

     TYPE(Variable_t), POINTER :: TimeVar 
     TYPE(Variable_t), POINTER :: HydPotSol, AreaSol, QcSol, QmSol 
     INTEGER, POINTER :: NodeIndexes(:), HydPotPerm(:),  &
              AreaPerm(:), QcPerm(:), QmPerm(:)
     REAL(KIND=dp), POINTER :: HydPot(:), AreaSolution(:), QcSolution(:), QmSolution(:)
     
     TYPE(ValueList_t), POINTER :: BC, Constants

     INTEGER :: i, j, k, l, m, n, t, iter, body_id, eq_id, material_id, &
          istat, LocalNodes, bc_id, DIM, NodeSheet, EdgeSheet, NbMoulin, NbMoulinAll, &
          GhostNodes, it, itOut, ierr
     INTEGER, ALLOCATABLE :: TableNodeSheet(:), TableMoulin(:)     

     CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, HydPotName, &
                OutPutFileName, nit

     LOGICAL :: Found, AllocationsDone = .FALSE.,  SubroutineVisited = .FALSE., &
                FirstTime = .TRUE., FirstVisit = .TRUE. 

     INTEGER, PARAMETER :: PVtuUnit=1300
     INTEGER :: VtuUnit, offset, Cpt, kk
     INTEGER, ALLOCATABLE :: EdgePointArray(:), EdgeOffsetArray(:), EdgeTypeArray(:)

     LOGICAL :: SaveVTU=.False. , VtuBinary=.False., FileVTU = .FALSE., OutPutQm = .FALSE. 

     CHARACTER :: lf
     CHARACTER(LEN=1024) :: OutStr
     CHARACTER(MAX_NAME_LEN) :: proc_number, VtuFile, PVtuFile, VtuFormat, VtuFileFormat, OutPutDirectoryName
     CHARACTER(MAX_NAME_LEN) :: ChFluxVarName, ChAreaVarName, QmVarName 

     REAL(KIND=dp) , ALLOCATABLE ::  tmparray(:,:), Flux(:) 
     REAL(KIND=dp) , ALLOCATABLE :: NodePointArray(:), ChannelAreaArray(:), ChannelFluxArray(:), MoulinInputArray(:) 
          
     LOGICAL, ALLOCATABLE ::  IsGhostNode(:)

     REAL(KIND=dp) :: x, y, z

     SAVE &
          ElementNodes, EdgeNodes,      &
          IsGhostNode,  OutPutFileName, OutPutDirectoryName, FileVTU, VtuBinary, it, itOut,  &
          AllocationsDone, FirstTime, FirstVisit, SolverName, HydPotName, M, &
          NodeSheet, TableNodeSheet, EdgeSheet, NbMoulin, TableMoulin, OutPutQm, Flux, &
          ChFluxVarName, ChAreaVarName, QmVarName 

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     SolverName = 'GlaDSchannelOut ('// TRIM(Solver % Variable % Name) // ')'

     AreaSol => Solver % Variable
     AreaPerm  => AreaSol % Perm
     AreaSolution => AreaSol % Values
     
     LocalNodes = COUNT( AreaPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

     DIM = CoordinateSystemDimension()

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
        N = Solver % Mesh % MaxElementNodes
        M = Model % Mesh % NumberOfNodes

        IF ( AllocationsDone ) THEN
           DEALLOCATE(                    &
                ElementNodes % x,         &
                ElementNodes % y,         &
                ElementNodes % z,         &
                EdgeNodes % x,         &
                EdgeNodes % y,         &
                EdgeNodes % z,         &
                IsGhostNode,              &
                TableNodeSheet, TableMoulin, Flux )
        END IF                           
        
        ALLOCATE(                                  &
             ElementNodes % x( N ),                &
             ElementNodes % y( N ),                &
             ElementNodes % z( N ),                &
             EdgeNodes % x( N ),                &
             EdgeNodes % y( N ),                &
             EdgeNodes % z( N ),                &
             IsGhostNode( M ),                     &
             TableNodeSheet(M), TableMoulin(M), Flux(M), &
             STAT=istat )

        IF ( istat /= 0 ) THEN
           CALL FATAL( SolverName, 'Memory allocation error' )
        ELSE
           CALL INFO(SolverName, 'Memory allocation done', level=1 )
        END IF

        IF ( ParEnv % PEs > 1 ) THEN !only if we have a parallel run
           IsGhostNode( 1:M ) = .FALSE.
           GhostNodes = 0
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

        AllocationsDone = .TRUE.
     END IF
      
     IF (FirstVisit) THEN 
        FirstVisit = .FALSE. 
        Constants => GetConstants()
        HydPotName = GetString( Constants, &
             'Hydraulic Potential Variable Name', Found )
        IF(.NOT.Found) THEN        
           CALL WARN(SolverName,'Keyword >Hydraulic Potential Variable Name< not found in section Constants')
           CALL WARN(SolverName,'Taking default value >Hydraulic Potential<')
           WRITE(HydPotName,'(A)') 'Hydraulic Potential'
        END IF
        ChAreaVarName = GetString( Constants, &
             'Channel Area Variable Name', Found )
        IF(.NOT.Found) THEN        
           CALL WARN(SolverName,'Keyword >Channel Area Variable Name< not found in section Constants')
           CALL WARN(SolverName,'Taking default value >Channel Area<')
           WRITE(ChAreaVarName,'(A)') 'Channel Area'
        END IF
        WRITE(QmVarName,'(A)') 'Flux from Moulins' 
        WRITE(ChFluxVarName,'(A)') 'Channel Flux' 
     END IF

     ! Point on the needed variables
     HydPotSol => VariableGet( Solver % Mesh % Variables, HydPotName, UnfoundFatal = .TRUE. )
     HydPotPerm     => HydPotSol % Perm
     HydPot => HydPotSol % Values

     QcSol => VariableGet( Solver % Mesh % Variables, ChFluxVarName )
     IF (ASSOCIATED(QcSol)) THEN
        QcPerm     => QcSol % Perm
        QcSolution => QcSol % Values
     END IF

     QmSol => VariableGet( Solver % Mesh % Variables, QmVarName )
     IF (ASSOCIATED(QmSol)) THEN
        QmPerm     => QmSol % Perm
        QmSolution => QmSol % Values
     END IF

     TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
     IF (FirstTime) THEN
          FirstTime = .FALSE. 
          itOut = 0
          it = 0
          FileVTU = GetLogical( Solver % Values, 'VTU OutPutFile' )
          VtuBinary  = GetLogical( Solver % Values, 'VTU BinaryFile')

          IF (FileVTU) THEN 
             OutPutDirectoryName = GetString( Solver % Values, 'Channels OutPut Directory Name', Found )
             IF (.NOT.Found) THEN
                WRITE(Message,'(a)') 'Keyword > Channels OutPut Directory Name < not found'
                WRITE(Message,'(a)') 'Output Files will be created in the .sif Directory'
             END IF
             OutPutFileName = GetString( Solver % Values, 'Channels OutPut File Name', Found )
             IF (.NOT.Found) THEN
                WRITE(Message,'(a)') 'Keyword > Channels OutPut File Name < not found'
                CALL FATAL(SolverName, Message)
             END IF
          END IF
                                           
          ! Save some information regarding the mesh for VTU output
          ! Number of nodes in the sheet layer: NodeSheet
          ! Number of channels (edges) in the sheet layer: EdgeSheet 
          ! Permutation table for the NodeSheet nodes and global nodes 
          IF (FileVTU) THEN
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
             CALL INFO(SolverName, Message, level=3 )
             WRITE(Message,'(a,i0)')'Number of nodes in the sheet layer: ', NodeSheet
             CALL INFO(SolverName, Message, level=3 )

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
             OutPutQm = .FALSE.
             IF ((NbMoulin>0).AND.(.Not.ASSOCIATED(QmSol))) THEN
                WRITE(Message,'(i0,a,i0)') &
                    NbMoulin, 'moulins found, but not variable >Flux from Moulins< associated, part ',ParEnv%myPe 
                CALL INFO(SolverName, Message, level=1 )
             ELSE IF (NbMoulin==0) THEN
                WRITE(Message,'(a,i0)')'No moulin found, part ',ParEnv%myPe  
                CALL INFO(SolverName, Message, level=3 )
             ELSE
                WRITE(Message,'(a,i0,a,i0)')'Qm will be saved for ',NbMoulin,' moulins in part ',ParEnv%myPe   
                CALL INFO(SolverName, Message, level=3 )
                OutPutQm = .TRUE. 
             END IF
          END IF
     END IF ! FirstTime

     ! Check the number of moulins on all partitions
     CALL MPI_ALLREDUCE( NbMoulin, NbMoulinAll, 1, MPI_INTEGER, &
                         MPI_SUM, ELMER_COMM_WORLD, ierr )
     WRITE(Message,'(a,i0,a)')'Qm will be saved for ',NbMoulinAll,' moulins in total'   
     CALL INFO(SolverName, Message, level=3 )
 
     ! If we are here, it means we have to save 
     ! assuming Exec Solver = After Saving is set in the solver section
     ! number of file is incremented by one
     itOut = itOut + 1

     SaveVTU = FileVTU 

     CALL Info( SolverName, ' Channels Output will be saved ', Level=4 )

     WRITE(Message,'(A,D15.7)')' Maximum Channel Area: ', &
        MAXVAL(AreaSolution(AreaPerm(M:M+Solver%Mesh%NumberOfEdges))) 
     CALL INFO(SolverName, Message, level=4 )

     ! Save results in VTU Format
     IF (SaveVTU) THEN
        lf = CHAR(10)
        WRITE(nit,'(i4.4)')itOut
        nit = ADJUSTL(nit)

        IF ( ParEnv%PEs >1 ) THEN
           WRITE(proc_number,'(i4.4)') ParEnv%myPe+1
           proc_number = ADJUSTL(proc_number)
           PVtuFile=TRIM(OutPutDirectoryName)//'/'//TRIM(OutPutFileName)//'_'//TRIM(nit)//'.pvtu'
           VtuFile=TRIM(OutPutDirectoryName)//'/'//TRIM(OutPutFileName)//'_'//TRIM(proc_number)//'par'//TRIM(nit)//'.vtu'
           VtuUnit = 1500 +ParEnv%myPe
        ELSE
           VtuUnit = 1500 
           VtuFile=TRIM(OutPutDirectoryName)//'/'//TRIM(OutPutFileName)//'_'//TRIM(nit)//'.vtu'
        END IF

        IF (VtuBinary) THEN
           VtuFileFormat="appended"
        ELSE
           VtuFileFormat="ascii"
        ENDIF

        IF (ParEnv%myPe == 0 .AND. ParEnv%PEs >1 ) THEN 
           OPEN( UNIT=PVtuUnit, FILE=PVtuFile, FORM = 'formatted', STATUS='unknown' )
           WRITE( PVtuUnit,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
           WRITE( PVtuUnit,'(A)') '  <PUnstructuredGrid>'
           WRITE( PVtuUnit,'(A)') '    <PPoints>'
           WRITE( PVtuUnit,'(A)') &
           '       <PDataArray type="Float64" NumberOfComponents="3" format="'//TRIM(VtuFileFormat)//'"/>'
           WRITE( PVtuUnit,'(A)') '    </PPoints>'
           WRITE( PVtuUnit,'(A)') '    <PCells>'
           WRITE( PVtuUnit,'(A)') &
           '       <PDataArray type="Int32" Name="connectivity" NumberOfComponents="1" format="'//TRIM(VtuFileFormat)//'"/>'
           WRITE( PVtuUnit,'(A)') &
           '       <PDataArray type="Int32" Name="offsets"      NumberOfComponents="1" format="'//TRIM(VtuFileFormat)//'"/>'
           WRITE( PVtuUnit,'(A)') &
           '       <PDataArray type="Int32" Name="types"        NumberOfComponents="1" format="'//TRIM(VtuFileFormat)//'"/>'
           WRITE( PVtuUnit,'(A)') '    </PCells>'
           WRITE( PVtuUnit,'(A)') '    <PCellData >'
           WRITE( PVtuUnit,'(A)') &
           '       <PDataArray type="Float64" Name="'//TRIM(ChAreaVarName)//'" format="'//TRIM(VtuFileFormat)//'"/>'
           WRITE( PVtuUnit,'(A)') &
           '       <PDataArray type="Float64" Name="'//TRIM(ChFluxVarName)//'" format="'//TRIM(VtuFileFormat)//'"/>'
           WRITE( PVtuUnit,'(A)') '    </PCellData>'
           IF (ASSOCIATED(QmSol) .AND. NbMoulinAll >0 ) THEN
              WRITE( PVtuUnit,'(A)') '    <PPointData>'
              WRITE( PVtuUnit,'(A)') &
              '       <PDataArray type="Float64" Name="'//TRIM(QmVarName)//'" format="'//TRIM(VtuFileFormat)//'"/>'
              WRITE( PVtuUnit,'(A)') '    </PPointData>'
           END IF
           DO i=1, ParEnv%PEs 
              WRITE(proc_number,'(i4.4)') i
              proc_number = ADJUSTL(proc_number)
              WRITE( PVtuUnit,'(A)') &
              '    <Piece Source="'//TRIM(OutPutFileName)//'_'//TRIM(proc_number)//'par'//TRIM(nit)//'.vtu" />'
           ENDDO
           WRITE( PVtuUnit,'(A)') '  </PUnstructuredGrid>'
           WRITE( PVtuUnit,'(A)') '</VTKFile>'
           CLOSE(PVtuUnit)
        END IF

        ALLOCATE(tmparray(3,NodeSheet))
        tmparray(:,:)=0

        ALLOCATE(NodePointArray(3*NodeSheet))
        NodePointArray(:)=0.0

        Cpt=0
        DO i = 1, Model % NumberOfNodes
           IF (TableNodeSheet(i)>0) THEN 
              x = Solver % Mesh % Nodes % x(i)
              y = Solver % Mesh % Nodes % y(i)
              z = Solver % Mesh % Nodes % z(i)
              Cpt = Cpt +1
              tmparray(1,Cpt)=x
              tmparray(2,Cpt)=y
              tmparray(3,Cpt)=z
           END IF
        END DO

        DO i=1,NodeSheet
           DO j=1,3
              kk=j+ (i-1)*3
              NodePointArray(kk)=tmparray(j,i)
           END DO
        END DO
        DEALLOCATE(tmparray)

        ALLOCATE(tmparray(2,EdgeSheet))
        ALLOCATE(EdgePointArray(2*EdgeSheet))
        ALLOCATE(EdgeOffsetArray(EdgeSheet))
        ALLOCATE(EdgeTypeArray(EdgeSheet))

        tmparray(:,:)=0
        EdgePointArray(:)=0
        EdgeOffsetArray(:)=0
        EdgeTypeArray(:)=0

        Cpt=0
        DO t=1, Solver % Mesh % NumberOfEdges 
           Edge => Solver % Mesh % Edges(t)
           IF (.NOT.ASSOCIATED(Edge)) CYCLE
           IF ((ParEnv % PEs > 1) .AND. &
              (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE
           n = Edge % TYPE % NumberOfNodes
           IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE
            
           Cpt = Cpt+1
           tmparray(1,Cpt) = TableNodeSheet(Edge % NodeIndexes(1))-1
           tmparray(2,Cpt) = TableNodeSheet(Edge % NodeIndexes(2))-1
        END DO

        DO i=1,EdgeSheet
           DO j=1,2
              kk=j+ (i-1)*2
              EdgePointArray(kk)=tmparray(j,i)
           END DO
        END DO

        DO j=1,EdgeSheet
           EdgeOffsetArray(j)= 2 * j  
           EdgeTypeArray(j)  = 3
        END DO

        DEALLOCATE(tmparray)

        ALLOCATE(ChannelAreaArray(EdgeSheet))

        ! Loops to get the Channels Area variable
        Cpt=0
        DO t=1, Solver % Mesh % NumberOfEdges 
           Edge => Solver % Mesh % Edges(t)
           IF (.NOT.ASSOCIATED(Edge)) CYCLE
           IF ((ParEnv % PEs > 1) .AND. &
              (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE
           n = Edge % TYPE % NumberOfNodes
           IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE
            
           M = Model % Mesh % NumberOfNodes
           Cpt = Cpt+1
           ChannelAreaArray(Cpt) = AreaSolution(AreaPerm(M+t))
        END DO


        ! Loops to get the Channels Flux variable
        IF (ASSOCIATED(QcSol)) THEN
           ALLOCATE(ChannelFluxArray(EdgeSheet)) 

           Cpt=0
           DO t=1, Solver % Mesh % NumberOfEdges 
              Edge => Solver % Mesh % Edges(t)
              IF (.NOT.ASSOCIATED(Edge)) CYCLE
              IF ((ParEnv % PEs > 1) .AND. &
                 (ParEnv % myPe .NE. Solver % Mesh % ParallelInfo % EdgeNeighbourList(t) % Neighbours(1))) CYCLE
              n = Edge % TYPE % NumberOfNodes
              IF (ANY(HydPotPerm(Edge % NodeIndexes(1:n))==0)) CYCLE
            
              M = Model % Mesh % NumberOfNodes
              Cpt = Cpt+1
              ChannelFluxArray(Cpt) = QcSolution(QcPerm(M+t))
           END DO
        END IF
         
        ! Loops to get the Flux from Moulins variable
        IF (OutPutQm) THEN 
           ALLOCATE(MoulinInputArray(NodeSheet))
           Flux = 0.0_dp
           DO t=1, Solver % Mesh % NumberOfBoundaryElements
              Element => GetBoundaryElement(t)
              n = GetElementNOFNodes()
              
              IF ( GetElementFamily() > 1 ) CYCLE
               
              BC => GetBC()
              bc_id = GetBCId( Element )
              CALL GetElementNodes( ElementNodes )
              IF ( ASSOCIATED( BC ) ) THEN            
                 j = Element % NodeIndexes(1)
                 Flux(j) = QmSolution(QmPerm(j)) 
              END IF
           END DO

           Cpt=0
           DO i = 1, Model % NumberOfNodes
              IF (TableNodeSheet(i)>0) THEN 
                 Cpt = Cpt+1
                 MoulinInputArray(Cpt) = Flux(i)
              END IF
           END DO
        END IF 

        IF( VtuBinary) THEN
           OPEN( UNIT=VtuUnit, FILE=VtuFile, FORM = 'unformatted', ACCESS = 'stream',  STATUS='unknown')
           WRITE( OutStr,'(A)' ) '<?xml version="1.0"?>  '//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '  <UnstructuredGrid>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="',NodeSheet,'" NumberOfCells="',EdgeSheet,'" >'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '         <Points>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           offset=0
           WRITE( OutStr,'(A,I0,A)') &
           '            <DataArray type="Float64" NumberOfComponents="3" format="appended" offset="',offset,'" />'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '         </Points>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '         <Cells>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           offset = offset + 8*NodeSheet*3 + 4
           WRITE( OutStr,'(A,I0,A)') &
           '            <DataArray type="Int32" Name="connectivity" NumberOfComponents="1" & 
                        & format="appended" offset="',offset,'" />'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           offset = offset + 4*EdgeSheet*2+4
           WRITE( OutStr,'(A,I0,A)') &
           '            <DataArray type="Int32" Name="offsets" NumberOfComponents="1" &
                        & format="appended" offset="',offset,'" />'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           offset = offset + 4*EdgeSheet+4
           WRITE( OutStr,'(A,I0,A)') &
           '            <DataArray type="Int32" Name="types"   NumberOfComponents="1" &
                        & format="appended" offset="',offset,'" />'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A)' ) '         </Cells>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '         <CellData>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           offset = offset + 4*EdgeSheet+4
           WRITE( OutStr,'(A,I0,A)') &
           '            <DataArray type="Float64" Name="'//TRIM(ChAreaVarName)//'" NumberOfComponents="1" &
                        & format="appended" offset="',offset,'" />'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           offset = offset + 8*EdgeSheet +4
           WRITE( OutStr,'(A,I0,A)') &
           '            <DataArray type="Float64" Name="'//TRIM(ChFluxVarName)//'" NumberOfComponents="1" &
                        & format="appended" offset="',offset,'" />'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '         </CellData>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           IF (OutPutQm) THEN 
              WRITE( OutStr,'(A)' ) '         <PointData>'//lf
              CALL VtuStrWrite( OutStr , VtuUnit )
              offset = offset +  8*EdgeSheet + 4
              WRITE( OutStr,'(A,I0,A)') &
              '            <DataArray type="Float64" Name="'//TRIM(QmVarName)//'" NumberOfComponents="1" &
           	             & format="appended" offset="',offset,'" />'//lf
              CALL VtuStrWrite( OutStr , VtuUnit )
              WRITE( OutStr,'(A)' ) '         </PointData>'//lf
              CALL VtuStrWrite( OutStr , VtuUnit )
           END IF

           WRITE( OutStr,'(A)' ) '      </Piece>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '  </UnstructuredGrid>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '<AppendedData encoding="raw">'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           IF (OutPutQm) THEN 
              WRITE(VtuUnit) '_',  KIND(NodePointArray)  *size(NodePointArray)  , NodePointArray(:)   ,&
                                   KIND(EdgePointArray)  *size(EdgePointArray)  , EdgePointArray(:)   ,&
                                   KIND(EdgeOffsetArray) *size(EdgeOffsetArray) , EdgeOffsetArray(:)  ,&
                                   KIND(EdgeTypeArray)   *size(EdgeTypeArray)   , EdgeTypeArray(:)    ,&
                                   KIND(ChannelAreaArray)*size(ChannelAreaArray), ChannelAreaArray(:) ,&
                                   KIND(ChannelFluxArray)*size(ChannelFluxArray), ChannelFluxArray(:) ,&
                                   KIND(MoulinInputArray)*size(MoulinInputArray), MoulinInputArray(:) 
           ELSE
              WRITE(VtuUnit) '_',  KIND(NodePointArray)  *size(NodePointArray)  , NodePointArray(:)   ,&
                                   KIND(EdgePointArray)  *size(EdgePointArray)  , EdgePointArray(:)   ,&
                                   KIND(EdgeOffsetArray) *size(EdgeOffsetArray) , EdgeOffsetArray(:)  ,&
                                   KIND(EdgeTypeArray)   *size(EdgeTypeArray)   , EdgeTypeArray(:)    ,&
                                   KIND(ChannelAreaArray)*size(ChannelAreaArray), ChannelAreaArray(:) ,&
                                   KIND(ChannelFluxArray)*size(ChannelFluxArray), ChannelFluxArray(:)
           END IF
           WRITE( OutStr,'(A)' ) lf//'</AppendedData>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
        ELSE
           OPEN( UNIT=VtuUnit, FILE=VtuFile, FORM = 'formatted',  ACCESS = 'sequential',  STATUS='unknown')
           WRITE( OutStr,'(A)' ) '<?xml version="1.0"?>  '//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '  <UnstructuredGrid>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="',NodeSheet,'" NumberOfCells="',EdgeSheet,'" >'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A)' ) '         <Points>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A)') '            <DataArray type="Float64" NumberOfComponents="3" format="ascii" >'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE(VtuFormat,'(A,I0,A,A)') &
                                  "(",size(NodePointArray), "ES16.7E3",")"
           WRITE(VtuUnit,VtuFormat) NodePointArray(:) 
           WRITE( OutStr,'(A)') '            </DataArray>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A)' ) '         </Points>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '         <Cells>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A)') &
           '            <DataArray type="Int32" Name="connectivity" NumberOfComponents="1" format="ascii" >'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE(VtuFormat,'(A,I0,A,A)') &
                                "(",size(EdgePointArray), '(" ",I0)',")"
           WRITE(VtuUnit,VtuFormat) EdgePointArray(:) 
           WRITE( OutStr,'(A)') '            </DataArray>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A)') &
           '            <DataArray type="Int32" Name="offsets" NumberOfComponents="1" format="ascii" >'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE(VtuFormat,'(A,I0,A,A)') &
                                "(",size(EdgeOffsetArray), '(" ",I0)',")"
           WRITE(VtuUnit,VtuFormat) EdgeOffsetArray(:) 
           WRITE( OutStr,'(A)') '            </DataArray>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A)') &
           '            <DataArray type="Int32" Name="types"   NumberOfComponents="1" format="ascii" >'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE(VtuFormat,'(A,I0,A,A)') &
                                "(",size(EdgeTypeArray), '(" ",I0)',")"
           WRITE(VtuUnit,VtuFormat) EdgeTypeArray(:) 
           WRITE( OutStr,'(A)') '            </DataArray>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A)' ) '         </Cells>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '         <CellData>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A)') &
           '            <DataArray type="Float64" Name="'//TRIM(ChAreaVarName)//'" NumberOfComponents="1" format="ascii" >'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE(VtuFormat,'(A,I0,A,A)') &
                                "(",size(ChannelAreaArray), "E16.7",")"
           WRITE(VtuUnit,VtuFormat) ChannelAreaArray(:) 
           WRITE( OutStr,'(A)') '            </DataArray>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A)') &
           '            <DataArray type="Float64" Name="'//TRIM(ChFluxVarName)//'" NumberOfComponents="1" format="ascii" >'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE(VtuFormat,'(A,I0,A,A)') &
                                "(",size(ChannelFluxArray), "E16.7",")"
           WRITE(VtuUnit,VtuFormat) ChannelFluxArray(:) 
           WRITE( OutStr,'(A)') '            </DataArray>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           WRITE( OutStr,'(A)' ) '         </CellData>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )

           IF (OutPutQm) THEN 
              WRITE( OutStr,'(A)' ) '         <PointData>'//lf
              CALL VtuStrWrite( OutStr , VtuUnit )
              WRITE( OutStr,'(A)') &
              '            <DataArray type="Float64" Name="'//TRIM(QmVarName)//'" NumberOfComponents="1" format="ascii" >'//lf
              CALL VtuStrWrite( OutStr , VtuUnit )
              WRITE(VtuFormat,'(A,I0,A,A)') &
                                   "(",size(MoulinInputArray), "E16.7",")"
              WRITE(VtuUnit,VtuFormat) MoulinInputArray(:) 
              WRITE( OutStr,'(A)') '            </DataArray>'//lf
              CALL VtuStrWrite( OutStr , VtuUnit )
              WRITE( OutStr,'(A)' ) '         </PointData>'//lf
              CALL VtuStrWrite( OutStr , VtuUnit )
           END IF

           WRITE( OutStr,'(A)' ) '      </Piece>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
           WRITE( OutStr,'(A)' ) '  </UnstructuredGrid>'//lf
           CALL VtuStrWrite( OutStr , VtuUnit )
        ENDIF
        WRITE( OutStr,'(A)' ) '</VTKFile>'//lf
        CALL VtuStrWrite( OutStr , VtuUnit )
        CLOSE(VtuUnit)

        DEALLOCATE(NodePointArray)  
        DEALLOCATE(EdgePointArray)  
        DEALLOCATE(EdgeOffsetArray) 
        DEALLOCATE(EdgeTypeArray)   
        DEALLOCATE(ChannelAreaArray)
        DEALLOCATE(ChannelFluxArray)
        IF (OutPutQm) DEALLOCATE(MoulinInputArray)
    END IF ! Output VTU

    SubroutineVisited = .TRUE.

CONTAINS    

SUBROUTINE VtuStrWrite( Str , VtuUnit )

    CHARACTER(LEN=1024) :: Str
    INTEGER :: VtuUnit 

    IF ( VtuBinary ) THEN
        WRITE( VtuUnit ) TRIM(Str)
    ELSE
        WRITE( VtuUnit, '(A)' ) TRIM(Str)
    END IF

END SUBROUTINE VtuStrWrite

!------------------------------------------------------------------------------
END SUBROUTINE GlaDSchannelOut
!------------------------------------------------------------------------------
