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
! *  Author: F. Gillet-Chaulet (IGE)
! *  Email:  fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:    http://elmerice.elmerfem.org
! *
! *  Original Date: 9 March  2018
! *  
! *  Save active elements and boundaries 
! *****************************************************************************
!!!  
      SUBROUTINE GetActiveMesh( Model,Solver,dt,TransientSimulation)
!------------------------------------------------------------------------------
      USE DefUtils
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Solver_t), POINTER ::PSolver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE(Mesh_t),POINTER :: Mesh
      TYPE(Mesh_t),POINTER :: NewMesh,PrevMesh
      TYPE(Variable_t),POINTER :: MeshSize
      TYPE(ValueList_t), POINTER :: SolverParams
 
      INTEGER, SAVE :: MeshNumber=0

      INTEGER :: ier
      INTEGER :: ii,kk

      CHARACTER(len=300) :: filename
      CHARACTER(LEN=1024) :: Path
      CHARACTER(LEN=MAX_NAME_LEN) :: OutPutFileName
      CHARACTER(LEN=MAX_NAME_LEN) :: MeshSizeName
      CHARACTER(LEN=MAX_NAME_LEN) :: DefaultFileName='ActiveMesh' 
      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='GetActiveMesh'

      LOGICAL :: DEBUG=.FALSE.,Scalar,Found
      LOGICAL :: SAVE_ELMER_MESH=.TRUE.
      LOGICAL :: RELEASE_MESH=.FALSE.
      LOGICAL :: UnFoundFatal=.TRUE.
      LOGICAL :: IncrementMeshNumber


      SolverParams => GetSolverParams()

!!
      MeshNumber = MeshNumber + 1

! NewMesh Name
      OutPutFileName = ListGetString(SolverParams,'Output file name',Found)
      IF (.NOT.Found) OutPutFileName = DefaultFileName
      IncrementMeshNumber=&
                 ListGetLogical(SolverParams,'Increment Mesh Number',Found)
      IF (.NOT.Found) IncrementMeshNumber=.TRUE.

      Mesh => Solver % Mesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GET THE NEW MESH
     NewMesh => GET_ACTIVE_MESH()

!  SAVE MESH TO DISK
      IF (SAVE_ELMER_MESH) THEN
       Path = TRIM(NewMesh % Name)
       CALL MakeDirectory( TRIM(path) // CHAR(0) )
       CALL WriteMeshToDisk2(Model, NewMesh, Path )
      END IF

! Get previous Mesh to do interpolation
!---------------------------------------
      PrevMesh => Model % Meshes

! Add Variable to NewMesh Structure + Model%mesh => NewMesh
!----------------------------------------------------------
      IF (DEBUG) PRINT *,'AddValueToMesh'
      CALL AddValueToMesh(NewMesh,PrevMesh) 


! Add the new mesh to the global list of meshes
!----------------------------------------------
      Model % Meshes   => NewMesh
      
! Update Solver Meshes
!---------------------
      DO ii=1,Model % NumberOfSolvers
        PSolver => Model%Solvers(ii)
        IF (.NOT.ASSOCIATED(PSolver)) CYCLE
        IF (.NOT.ASSOCIATED(PSolver%Matrix)) THEN
          PSolver % Mesh => NewMesh
        ELSE
          CALL UpdateSolverMesh(PSolver,NewMesh) ! call distrib version in MeshUtils.F90
          WRITE(Message,'(A,A,A)') 'Updated for variable: ',TRIM(PSolver%Variable%Name)
          CALL INFO(SolverName,TRIM(Message),level=1)
        END IF
      END DO

! Release previous mesh (shoudl be optionnal)
!-------------------
      RELEASE_MESH=ListGetLogical(SolverParams,'Release previous mesh')
      IF (RELEASE_MESH) THEN
        CALL ReleaseMesh(PrevMesh)
        Deallocate(PrevMesh)
        Model % Meshes % Next => NULL()
      END IF

     NewMesh % Changed = .TRUE.
! Print info
!-----------
      write(Message,'(A)') '--**--New mesh ready'
       CALL INFO(SolverName,trim(Message),level=5)
      write(Message,'(A,I0)') '--**--    NumberOfNodes:',NewMesh%NumberOfNodes
       CALL INFO(SolverName,trim(Message),level=5)
      write(Message,'(A,I0)') '--**--    NumberOfBulkElements: ',NewMesh%NumberOfBulkElements
       CALL INFO(SolverName,trim(Message),level=5)
      write(Message,'(A,I0)') '--**--    NumberOfEdges: ',NewMesh%NumberOfEdges
       CALL INFO(SolverName,trim(Message),level=5)
      write(Message,'(A,I0)') '--**--    NumberOfBoundaryElements: ',NewMesh%NumberOfBoundaryElements
       CALL INFO(SolverName,trim(Message),level=5)
      
!-------------------------------------
!! Subroutine
!-------------------------------------
     CONTAINS

      FUNCTION GET_ACTIVE_MESH() RESULT(NewMesh)
      IMPLICIT NONE
!------------------------------------------------------------------------------
      TYPE(Mesh_t), POINTER :: NewMesh
!------------------------------------------------------------------------------
      TYPE(Element_t),POINTER ::  Element,Parent,PrevElement
      INTEGER, POINTER :: NodeIndexes(:)
      INTEGER :: n,np,nt,na
      INTEGER :: cmpt
      INTEGER :: i,tt, jj, kk, ll
      INTEGER, ALLOCATABLE :: ActiveNode(:),ActiveNodeTable(:)
      INTEGER, ALLOCATABLE :: ActiveElementTable(:)
      INTEGER, ALLOCATABLE :: ActiveBCElementTable(:)

      n=Solver % Mesh % NumberOfNodes
      ALLOCATE(ActiveNode(n),ActiveNodeTable(n))

      n=Solver % Mesh % NumberOfBulkElements
      ALLOCATE(ActiveElementTable(n))

      n=Solver % Mesh % NumberOfBoundaryElements
      ALLOCATE(ActiveBCElementTable(n))


! GET ACTIVE ELEMENTS
      ActiveNode=-1
      ActiveElementTable=-1
      nt=0
      DO tt=1,Solver % Mesh % NumberOfBulkElements
         Element => Solver % Mesh % Elements(tt)
         IF ( CheckPassiveElement(Element) )  CYCLE

         nt=nt+1
         ActiveElementTable(tt)=nt


         n = GetElementNOFNodes(Element)
         DO i=1,n
            ActiveNode(Element % NodeIndexes(i))=+1
         END DO

      END DO
      np=COUNT(ActiveNode.EQ.1)

! GET ACTIVE NODES
      cmpt=0
       ActiveNodeTable=-1
      DO jj=1,Solver % Mesh % NumberOfNodes
         IF (ActiveNode(jj).EQ.+1) THEN
            cmpt=cmpt+1
            ActiveNodeTable(jj)=cmpt
         END IF
      END DO

!GET ACTIVE BC ELEMENTS
     na=0
     ActiveBCElementTable=-1
     DO tt=1,Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(tt)
        Parent => Element % BoundaryInfo % Left

        !check for Left Parent 
        IF (ASSOCIATED(Parent)) THEN
           !if associated check if pasive or not
           IF (CheckPassiveElement(Parent)) &
               Parent => Element % BoundaryInfo % Right
        ELSE
           Parent => Element % BoundaryInfo % Right
        ENDIF
        ! Parent is either Left and active or right and need to check
        ! status
        IF (ASSOCIATED(Parent)) THEN
          IF (.NOT.CheckPassiveElement(Parent)) THEN 
            na=na+1
            ActiveBCElementTable(na)=tt
          END IF
        END IF
      END DO


! INITIALISE THE NEW MESH STRUCTURE
      NewMesh => AllocateMesh()
      IF (IncrementMeshNumber) THEN
         write(NewMesh % Name,'(A,A,I0)') TRIM(OutPutFileName),'_N',MeshNumber
      ELSE
         NewMesh % Name=TRIM(OutPutFileName)
      END IF

      NewMesh%MaxElementNodes=Solver%Mesh%MaxElementNodes
      NewMesh%MaxElementDOFs=Solver%Mesh%MaxElementDOFs
      NewMesh%MeshDim=Solver%Mesh%MeshDim

      NewMesh % NumberOfNodes = np
      NewMesh % NumberOfBulkElements = nt
      NewMesh % NumberOfBoundaryElements = na

      CALL AllocateVector( NewMesh % Nodes % x, NewMesh % NumberOfNodes)
      CALL AllocateVector( NewMesh % Nodes % y, NewMesh % NumberOfNodes)
      CALL AllocateVector( NewMesh % Nodes % z, NewMesh % NumberOfNodes)
      CALL AllocateVector( NewMesh % Elements, NewMesh % NumberOfBulkElements + &
               NewMesh % NumberOfBoundaryElements )
      
!! GET NEW VERTICES
      Do ii=1,Solver % Mesh % NumberOfNodes 
        kk=ActiveNodeTable(ii)
        IF (kk.LT.1) CYCLE
        NewMesh % Nodes % x(kk) = Solver % Mesh % Nodes % x(ii)
        NewMesh % Nodes % y(kk) = Solver % Mesh % Nodes % y(ii)
        NewMesh % Nodes % z(kk) = Solver % Mesh % Nodes % z(ii)
        IF (kk.GT.NewMesh % NumberOfNodes) &
            CALL FATAL(SolverName,'Pb in counting nodes')
      End do
      IF (DEBUG) PRINT *,'Get_vertex DONE'

!! GET NEW BulkElements
      cmpt=0
      Do tt=1,Solver % Mesh % NumberOfBulkElements
         IF (ActiveElementTable(tt).LT.1) CYCLE
         PrevElement => Solver % Mesh % Elements(tt)
          
         cmpt=cmpt+1
         IF (cmpt.GT.NewMesh % NumberOfBulkElements) &
            CALL FATAL(SolverName,'Pb in counting BulkElements')
         Element => NewMesh % Elements(cmpt)

         Element % TYPE => GetElementType( PrevElement % TYPE % ElementCode )
         n = Element % TYPE % NumberOfNodes
         Element % NDOFs = n
         Element % ElementIndex = cmpt
         Element % PartIndex = ParEnv % myPE
         CALL AllocateVector(Element % NodeIndexes, n )
         Element % NodeIndexes(1:n) = ActiveNodeTable(PrevElement % NodeIndexes(1:n))
         Element % BodyId = PrevElement % BodyId
         IF (ANY(Element % NodeIndexes.LT.1)) &
           CALL FATAL(SolverName,'Pb with node Indexes in BulkElements')
      End do
      IF (DEBUG) PRINT *,'Get_Bulks DONE'
      
!! Get BC Elements
      kk=NewMesh % NumberOfBulkElements
      Do ii=1,na
         kk = kk + 1

         Element => NewMesh % Elements(kk)

         PrevElement => GetBoundaryElement(ActiveBCElementTable(ii))

         Element % TYPE => GetElementType( PrevElement % TYPE % ElementCode )
         n = Element % TYPE % NumberOfNodes
         Element % NDOFs = n
         Element % ElementIndex = kk
         Element % PartIndex = ParEnv % myPE

         CALL AllocateVector(Element % NodeIndexes, n )

         Element % NodeIndexes(1:n) = ActiveNodeTable(PrevElement % NodeIndexes(1:n))
         IF (ANY(Element % NodeIndexes.LT.1)) &
           CALL FATAL(SolverName,'Pb with node Indexes in BCElements')
        
         Allocate(Element % BoundaryInfo)
         Element % BoundaryInfo % Constraint = &
            PrevElement % BoundaryInfo % Constraint

         Parent => PrevElement % BoundaryInfo % Left
         IF (ASSOCIATED(Parent)) THEN
            IF (ActiveElementTable(Parent % ElementIndex).GE.1) &
               Element % BoundaryInfo % Left => &
                 NewMesh % Elements(ActiveElementTable(Parent % ElementIndex))
         END IF

         Parent => PrevElement % BoundaryInfo % Right
         IF (ASSOCIATED(Parent)) THEN
            IF (ActiveElementTable(Parent % ElementIndex).GE.1) &
               Element % BoundaryInfo % Right => &
                 NewMesh % Elements(ActiveElementTable(Parent % ElementIndex))
         END IF

      End Do

      DEALLOCATE(ActiveNode,ActiveNodeTable,ActiveElementTable,ActiveBCElementTable)

      END FUNCTION GET_ACTIVE_MESH

!------------------------------------------------------------------------------
        SUBROUTINE AddValueToMesh(NewMesh,RefMesh)
!------------------------------------------------------------------------------
          implicit none
          TYPE(Solver_t) :: Solver
          TYPE(Mesh_t), POINTER :: NewMesh, RefMesh
          TYPE(Variable_t), POINTER :: Var,NewVar
          TYPE( Matrix_t ), POINTER :: NewMatrix
          INTEGER :: ii,k,jj
          INTEGER, POINTER :: Permutation(:)
          LOGICAL :: BandwidthOptimize, GlobalBubbles 
 
          CALL Info('AddVauleToMesh','Enter',Level=1)
          NewMesh % MaxBDOFs = RefMesh % MaxBDOFs
        ! Initialize local variables for the new mesh:
          NULLIFY( NewMesh % Variables )
          CALL VariableAdd( Newmesh % Variables, Newmesh, Solver, &
             'Coordinate 1', 1, NewMesh % Nodes % x )
          CALL VariableAdd( Newmesh % Variables, NewMesh, Solver, &
             'Coordinate 2', 1, Newmesh % Nodes % y )
          CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
             'Coordinate 3', 1, Newmesh % Nodes % z )
             
        ! Time must always be there:
        ! --------------------------
          Var => VariableGet( RefMesh % Variables,'Time',ThisOnly=.TRUE. )
          CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
             'Time', 1, Var % Values )
             
          Var => VariableGet( RefMesh % Variables,'Timestep',&
              ThisOnly=.TRUE.)
          CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
              'Timestep', 1, Var % Values ) 
              
          Var => VariableGet( RefMesh % Variables,'Timestep size',&
              ThisOnly=.TRUE. )
          CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
              'Timestep size', 1, Var % Values )
              
          Var => VariableGet( RefMesh % Variables,&
              'Timestep interval',ThisOnly=.TRUE. )
          CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
              'Timestep interval', 1, Var % Values )
               
          Var => VariableGet( RefMesh % Variables,'Coupled iter',&
               ThisOnly=.TRUE. )
          CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
              'Coupled iter', 1, Var % Values )
              
          Var => VariableGet( RefMesh % Variables,'Nonlin iter',&
              ThisOnly=.TRUE. )
          CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
              'Nonlin iter', 1, Var % Values )

        ! Set Mesh for model
        !-------------------
          CALL SetCurrentMesh(Model,NewMesh)
        
        ! Initialize the field variables for the new mesh. These are
        ! interpolated from the old meshes variables. Vector variables
        ! are in the variable lists in two ways: as vectors and as
        ! vector components.
        ! We MUST update the vectors (i.e. DOFs>1) first!
        ! ------------------------------------------------
          Var => RefMesh % Variables
          DO WHILE( ASSOCIATED( Var ) )
            IF ( Var % DOFs > 1 ) THEN
              NewVar => VariableGet( NewMesh % Variables,Var %Name,.FALSE. )
              kk = SIZE( NewVar % Values )
              IF ( ASSOCIATED( NewVar % Perm ) ) THEN
                Kk = COUNT( NewVar % Perm > 0 )
              END IF
              IF ( GetVarName( NewVar ) == 'flow solution' ) THEN
                NewVar % Norm = 0.0d0
                DO ii=1,NewMesh % NumberOfNodes
                  DO jj=1,NewVar % DOFs-1
                    NewVar % Norm = NewVar % Norm + &
                        NewVar % Values( NewVar % DOFs*(ii-1)+jj )**2
                  END DO
                END DO
                NewVar % Norm = SQRT( NewVar % Norm / kk )
              ELSE
                NewVar % Norm = SQRT( SUM(NewVar % Values**2) / kk )
              END IF
            END IF
            Var => Var % Next
          END DO
        
        !   Second time around, update scalar variables and
        !   vector components:
        !   -----------------------------------------------
          Var => RefMesh % Variables
          DO WHILE( ASSOCIATED( Var ) )
            IF( SIZE( Var % Values ) == Var % DOFs ) THEN
              Var => Var % Next
              CYCLE
            END IF
            SELECT CASE( Var % Name )
             CASE('coordinate 1','coordinate 2','coordinate 3','time',&
                  'timestep', 'timestep size', 'timestep interval', &
                  'coupled iter', 'nonlin iter' )
             CASE DEFAULT
             IF ( Var % DOFs == 1 ) THEN
               Found = .FALSE.
               IF ( Found ) THEN
                 kk = Solver % Variable % NameLen
                 IF ( Var % Name(1:kk) /= Solver % Variable % Name(1:kk)) THEN
                   Var => Var % Next
                   CYCLE
                 END IF
               END IF

               NewVar => VariableGet( NewMesh % Variables, Var % Name,.FALSE. )
               ! added by fab:
               NewVar % PrimaryMesh => NewMesh
               !---------------
               kk = SIZE( NewVar % Values )
               IF ( ASSOCIATED( NewVar % Perm ) ) THEN
                 kk = COUNT( NewVar % Perm > 0 )
               END IF
               NewVar % Norm = SQRT( SUM(NewVar % Values**2) / kk )
             END IF
            END SELECT
            Var => Var % Next
          END DO
        
        ! update solver structure to use the new mesh
        !--------------------------------------------
          !Solver % Mesh => NewMesh ! deleted by fab
          CALL MeshStabParams(NewMesh)

        ! Nothing computed on this mesh
        !------------------------------          
          NewMesh % SavesDone = 0 ! start new output file -> to check
          NewMesh % OutputActive = .TRUE.
          NewMesh % changed = .TRUE. 
          NewMesh % Next => RefMesh

!------------------------------------------------------------------------------
        End ! Subroutine AddVlueToMesh
!------------------------------------------------------------------------------
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END SUBROUTINE GetActiveMesh
