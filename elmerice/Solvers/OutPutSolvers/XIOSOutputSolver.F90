!------------------------------------------------------------------------------
      SUBROUTINE XIOSOutputSolver_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
      USE DefUtils
      USE xios
      IMPLICIT NONE
!------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: Transient
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: SolverParams

        SolverParams => GetSolverParams()

        CALL ListAddNewLogical( SolverParams,'Optimize Bandwidth',.FALSE.)
        ! need to skip halo
        CALL ListAddNewLogical( SolverParams,'Skip Halo Elements',.TRUE.)

      END SUBROUTINE XIOSOutputSolver_Init0
!------------------------------------------------------------------------------
      SUBROUTINE XIOSOutputSolver_Init(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
      USE DefUtils
      IMPLICIT NONE
!------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: Transient
!------------------------------------------------------------------------------
      CHARACTER(LEN=MAX_NAME_LEN) :: Name
      TYPE(ValueList_t), POINTER :: SolverParams
      LOGICAL :: GotIt

      SolverParams => Solver % Values

      Name = ListGetString( SolverParams, 'Equation',GotIt)
      IF( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
        CALL ListAddString( SolverParams,'Variable',&
           '-nooutput '//TRIM(Name)//'_var')
      ENDIF

      END SUBROUTINE XIOSOutputSolver_Init
!------------------------------------------------------------------------------
      SUBROUTINE XIOSOutputSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
      USE DefUtils
      USE SolverUtils
      USE SaveUtils
      USE ProjUtils
      USE Netcdf
      USE xios

      IMPLICIT NONE
      TYPE(Solver_t) :: Solver
      TYPE(Model_t) :: Model
      REAL(dp) :: dt
      LOGICAL :: TransientSimulation

      CHARACTER(*), PARAMETER :: Caller = 'XIOSOutputSolver'

      TYPE(Element_t), POINTER :: Element
      INTEGER, SAVE :: nTime = 0
      INTEGER,SAVE :: ElemFirst, ElemLast
      INTEGER, ALLOCATABLE, TARGET,SAVE :: NodePerm(:),InvNodePerm(:), InvDgPerm(:), DgPerm(:)
      LOGICAL, ALLOCATABLE,SAVE :: ActiveElem(:),NotOwnedNode(:),NotOwnedEdge(:)
      REAL(KIND=dp),DIMENSION(:),ALLOCATABLE,SAVE :: cell_area,Node_x,Node_y
      REAL(KIND=dp), ALLOCATABLE, SAVE :: BoundaryCondition(:)
      REAL(KIND=dp),ALLOCATABLE,SAVE :: Basis(:), dBasisdx(:,:)
      INTEGER,SAVE :: NumberOfGeomNodes,NumberOfDofNodes,NumberOfElements
      INTEGER,SAVE :: NumberOfActiveNodes,NumberOfActiveEdges
      LOGICAL,SAVE :: NoPermutation
      INTEGER :: i,ii,M,t,n
      LOGICAL :: ierr
      LOGICAL :: BoundarySolver


      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(ValueList_t),POINTER :: Params
      TYPE(ValueList_t), POINTER :: BC
      LOGICAL :: GotIt
      LOGICAL :: Parallel
      INTEGER :: GroupId
      INTEGER :: PEs,Part 
      REAL(KIND=dp) :: SimTime

      LOGICAL :: SaveLinear=.TRUE. ! save only corners
      INTEGER,SAVE :: MaxElementNodes  ! // MaxNumNodesPerFace
      INTEGER,PARAMETER :: Connect_Fill=-1
      LOGICAL :: ok

      CHARACTER(LEN=MAX_NAME_LEN) :: strg_var
      TYPE(xios_duration) :: dtime,time_units
      TYPE(xios_duration) :: sync_freq,out_freq
      TYPE(xios_context) :: ctx_hdl
      TYPE(xios_date) :: date_origin,ref_date
      TYPE(xios_date),SAVE :: date1,date2
      TYPE(xios_duration) :: dtvisit
      REAL(KIND=dp) :: ndt,dt1,dt2
      REAL(KIND=dp) :: tol
      CHARACTER(len=20) :: date_str
      INTEGER :: ts,ol

      INTEGER, SAVE :: Olevel=4
      LOGICAL,SAVE :: AllwaysSend

      LOGICAL,SAVE :: SaveEdges=.TRUE.
      LOGICAL  :: SkipEdges

      Params => GetSolverParams()
      Mesh => Model % Mesh

      PEs = ParEnv % PEs
      Part = ParEnv % MyPE
      Parallel = (PEs > 1)

      nTime = nTime + 1
      IF ((nTime>1).AND.(Mesh%Changed)) &
        CALL FATAL(Caller,"mesh has changed; not supported")

      
      !------------------------------------------------------------------------------
      ! Create and initialise file at first visit
      !------------------------------------------------------------------------------
      IF ( nTime == 1 ) THEN
        SkipEdges=ListGetLogical( Params,'Skip Edges',GotIt)
        IF (SkipEdges) SaveEdges=.FALSE.

        BoundarySolver = ( Solver % ActiveElements(1) > Model % Mesh % NumberOfBulkElements )
        IF (BoundarySolver) THEN
           !------------------------------------------------------------------------------
           ! set params to activate saving on the current boundary
           !------------------------------------------------------------------------------
           CALL ListAddNewLogical( Params,'Save Boundaries Only',.TRUE.)
           CALL ListAddNewString( Params,'Mask Name','MaskXIOS')
           BC => GetBC(Model%Mesh%Elements(Solver % ActiveElements(1)))
           IF (.NOT.ASSOCIATED(BC)) &
              CALL FATAL(Caller,"First element not associated to a BC")
           CALL ListAddNewLogical( BC,'MaskXIOS',.TRUE.)

           ! skip edges if working on a boundary....
           SaveEdges=.FALSE.
        ELSE
           CALL ListAddNewLogical( Params,'Save Bulk Only',.TRUE.)
           IF (Mesh % MeshDim.NE.2) &
              CALL FATAL(Caller,"Saving bulk values works only for 2D meshes")
        ENDIF



        IF (SaveEdges) THEN
          !! Get the edges
          CALL FindMeshEdges2D(Mesh)
          ! temporary trick to get correct interface for halo..
          IF (Parallel) THEN
            Mesh % MeshDim = 3
            CALL SParEdgeNumbering(Mesh)
            Mesh % MeshDim = 2
          END IF
        END IF

        ! can be use to set the output level for variables that are
        ! requested and send
        ol=ListGetInteger( Params,'Solver info level',GotIt)
        IF (GotIt) Olevel=ol

        ! do we send data only if requested or allways
        ! should be mainly for debugging...
        AllwaysSend=ListGetLogical(Params,"Allways send",Gotit)

        M = Model % MaxElementNodes
        ALLOCATE(Basis(M),dBasisdx(M,3))
        !------------------------------------------------------------------------------
        ! Initialize stuff for masked saving
        !------------------------------------------------------------------------------
        GroupId = 0
        CALL GenerateSaveMask(Mesh,Params,Parallel,GroupId,SaveLinear,&
               NodePerm,ActiveElem,NumberOfGeomNodes,NumberOfElements,&
               ElemFirst,ElemLast)
        !------------------------------------------------------------------------------
        ! If we have a discontinuous mesh then create the permutation vectors to deal
        ! with the discontinuities.
        !------------------------------------------------------------------------------
        CALL GenerateSavePermutation(Mesh,.FALSE.,.FALSE.,0,SaveLinear,ActiveElem,NumberOfGeomNodes,&
               NoPermutation,NumberOfDofNodes,DgPerm,InvDgPerm,NodePerm,InvNodePerm)

        !------------------------------------------------------------------------------
        ! Get Max ElementNodes from the active elements
        !------------------------------------------------------------------------------
        MaxElementNodes=0
        DO ii = ElemFirst, ElemLast
          IF( .NOT. ActiveElem(ii) ) CYCLE
          Element => Model % Elements(ii)
          IF (SaveLinear) THEN
            n = GetElementCorners( Element )
          ELSE
            n = GetElementNOFNodes(Element)
          END IF
          IF (n.GT.MaxElementNodes) &
            MaxElementNodes=n
        END DO
        !------------------------------------------------------------------------------
        ! Parallel case we will exclude nodes not owned by partition...
        !------------------------------------------------------------------------------
        ALLOCATE(NotOwnedNode(NumberOfDofNodes))
        NotOwnedNode=.FALSE.
        NumberOfActiveNodes=NumberOfDofNodes
        IF (Parallel) THEN
          DO ii=1,NumberOfDofNodes
            IF( NoPermutation ) THEN
              i = ii
            ELSE
              i = InvNodePerm(ii)
            END IF
            IF (Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1).NE.ParEnv % MyPE) THEN
              NotOwnedNode(ii)=.TRUE.
              NumberOfActiveNodes=NumberOfActiveNodes-1
            END IF
          END DO
        END IF

        !------------------------------------------------------------------------------
        ! Edges
        !------------------------------------------------------------------------------
        IF (SaveEdges) THEN
          ALLOCATE(NotOwnedEdge(Mesh % NumberOfEdges))
          NotOwnedEdge=.FALSE.

          NumberOfActiveEdges=Mesh % NumberOfEdges
          IF (Parallel) THEN
            DO i=1,Mesh % NumberOfEdges
             ! Edge at the interface and not owned
             IF (Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours(1).NE.ParEnv % MyPE) THEN
              NotOwnedEdge(i)=.TRUE.
              NumberOfActiveEdges=NumberOfActiveEdges-1
             ENDIF
            END DO
          END IF
        END IF

        ! The partition is active for saving if there are any nodes
        ! to write. There can be no elements nor dofs without nodes.
        CALL ParallelActive( NumberOfDofNodes > 0 )

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! XIOS context definition
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL xios_context_initialize(TRIM(xios_id),ELMER_COMM_WORLD)
        CALL xios_set_current_context(TRIM(xios_id))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! check that mandatory fields have been defined
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL SanityCheck()

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Set-up the time step from elmer....
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        strg_var= ListGetString(Params,"time_units",UnFoundFatal=.TRUE.)
        time_units=xios_duration_convert_from_string(TRIM(strg_var)) 

        strg_var = ListGetString(Params,"timestep",GotIt)
        IF (Gotit) THEN
           dtime =  xios_duration_convert_from_string(TRIM(strg_var))
           ! check we are consistent:
           tol=ListGetConstReal(Params,"timestep tolerance",UnFoundFatal=.TRUE.)
           CALL xios_get_time_origin(date_origin)
           ! n time-steps per time_units
           ndt=1/dt
           dt1=xios_date_convert_to_seconds(date_origin+time_units)
           dt2=xios_date_convert_to_seconds(date_origin+dtime*ndt)
           IF (abs(dt1-dt2).GT.tol) THEN
             WRITE( Message,'(A,ES12.3)') & 
                     "timestep difference (s) :",dt1-dt2
             CALL FATAL(Caller,Message)
           END IF
        ELSE
           dtime = dt * time_units
        ENDIF
        CALL xios_set_timestep(dtime)

        CALL xios_duration_convert_to_string(dtime,date_str)

        CALL Info(Caller,"#########################",Level=3)
        CALL Info(Caller,"Time step: "//TRIM(date_str),Level=3)
        CALL Info(Caller,"#########################",Level=3)

        strg_var= ListGetString(Params,"reference date",GotIt)
        IF (Gotit) THEN
          ref_date=xios_date_convert_from_string(TRIM(strg_var)//" 00:00:00")
          ref_date=ref_date+(GetTime()-dt)*time_units
          CALL xios_set_start_date(start_date=ref_date)
        END IF

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! set out_freq and sync_freq to the top file_definition
        !! if they are provided
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        strg_var=ListGetString(Params,"output frequency",GotIt)
        IF (GotIt) THEN
                out_freq=xios_duration_convert_from_string(TRIM(strg_var))
                CALL xios_set_filegroup_attr("file_definition",output_freq=out_freq)
        END IF
        strg_var=ListGetString(Params,"sync frequency",GotIt)
        IF (GotIt) THEN
                sync_freq=xios_duration_convert_from_string(TRIM(strg_var))
                CALL xios_set_filegroup_attr("file_definition",sync_freq=sync_freq)
        END IF

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Update File names suffix
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        strg_var=ListGetString(Params,"file names suffix",GotIt)
        IF (GotIt) call UpdateFileSuffix(TRIM(strg_var))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! set Elmer version number if a variable "elmerversion" exist
        !! especially for a global file attribute
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ok=xios_is_valid_variable("elmerversion")
        IF (ok) &
          ok = xios_setVar("elmerversion", "Elmer/Ice v"//TRIM(GetVersion())//" (Rev: "//TRIM(GetRevision())//")")

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! domain definition, a.k.a mesh topology
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL DomainDefinitions()

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! define default fields related to the mesh
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL DefaultFieldDefinitions()

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! close elemer/ice context definition
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL xios_close_context_definition

        CALL Info(Caller, 'Context initialisation done')

      END IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Update calendar
      !  should we add a sanity check that dt has not changed?
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL xios_update_calendar(nTime)
      IF (nTime==1) THEN
              CALL xios_get_current_date (date1)
      ELSE
              date1=date2
      ENDIF
      CALL xios_get_current_date (date2)
      dtvisit  = date2 - date1

      CALL xios_date_convert_to_string(date2, date_str)
      CALL Info(Caller,"Current date: "//TRIM(date_str),Level=3)
      WRITE( Message,'(A,ES12.3)') "Fraction of year:",xios_date_get_fraction_of_year(date2)
      CALL Info(Caller,Message,Level=Olevel) 

      CALL xios_duration_convert_to_string(dtvisit,date_str)
      CALL Info(Caller,"Duration since last visit:"//TRIM(date_str), &
                Level=Olevel) 



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! send requested variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL SendMeshVariables()
      CALL SendGlobalVariables()
      CALL SendVariables()

      CONTAINS 


      SUBROUTINE UpdateFileSuffix(suffix)
        IMPLICIT NONE
        CHARACTER(LEN=*) suffix
        CHARACTER(LEN=3) fileid
        INTEGER :: i

        DO i = 1, 9
         WRITE(fileid,'(i1)') i
         CALL SetFileSuffix('file'//TRIM(fileid),suffix)
        END DO
        DO i = 1, 99
         WRITE(fileid,'(i2.2)') i
         CALL SetFileSuffix('file'//TRIM(fileid),suffix)
        END DO
      END SUBROUTINE 

      SUBROUTINE SetFileSuffix(fileid,suffix)
        IMPLICIT NONE
        CHARACTER(LEN=*) suffix
        CHARACTER(LEN=*) fileid    

        IF( xios_is_valid_file(fileid) )  THEN
                CALL xios_set_file_attr(fileid,name_suffix=suffix)
        ENDIF
        IF( xios_is_valid_filegroup(fileid) )  THEN
                CALL xios_set_filegroup_attr(fileid,name_suffix=suffix)
        ENDIF

        !CALL xios_solve_inheritance()
      END SUBROUTINE SetFileSuffix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! check that required elements are present
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SanityCheck()
         IMPLICIT NONE
         LOGICAL :: IsValid
         CHARACTER(LEN=MAX_NAME_LEN) :: FieldName,Txt
         INTEGER :: Vari
         LOGICAL ::ScalarsExist
         LOGICAL :: Found,Project

         ! Grids
         IsValid=xios_is_valid_grid("GridCells")
         IF (.NOT.IsValid) &
                 CALL FATAL(Caller,"<GridCells> grid not defined")

         IsValid=xios_is_valid_grid("GridNodes")
         IF (.NOT.IsValid) &
                 CALL FATAL(Caller,"<GridNodes> grid not defined")

         IF (SaveEdges) THEN
           IsValid=xios_is_valid_grid("GridEdges")
           IF (.NOT.IsValid) &
                 CALL FATAL(Caller,"<GridEdges> grid not defined")
         END IF

         ! domains
         IsValid=xios_is_valid_domain("cells")
         IF (.NOT.IsValid) &
                 CALL FATAL(Caller,"<cells> domain not defined")

         IsValid=xios_is_valid_domain("nodes")
         IF (.NOT.IsValid) &
                 CALL FATAL(Caller,"<nodes> domain not defined")

         IF (SaveEdges) THEN
           IsValid=xios_is_valid_domain("edges")
           IF (.NOT.IsValid) &
             CALL FATAL(Caller,"<edges> domain not defined")
         END IF

        ! requested global variables
        FieldName = GetString( Params,'Global Variable 1',ScalarsExist)
        Vari=1
        DO WHILE (ScalarsExist)
          IsValid=xios_is_valid_field(TRIM(FieldName))
          IF (.NOT.IsValid) &
             CALL FATAL(Caller,TRIM(FieldName)//" is not a valid field")

          Vari=Vari+1
          WRITE(Txt,'(A,I0)') 'Global Variable ',Vari
          FieldName = GetString( Params,TRIM(Txt),ScalarsExist)
        END DO

        ! requested scalar variables
        FieldName = GetString( Params,'Scalar Field 1',ScalarsExist)
        Vari=1
        DO WHILE (ScalarsExist)
          IsValid=xios_is_valid_field(TRIM(FieldName))
          IF (.NOT.IsValid) &
             CALL FATAL(Caller,TRIM(FieldName)//" is not a valid field")

          WRITE(Txt,'(A,I0,A)') &
                'Scalar Field ',Vari,' compute cell average'
          Project = ListGetLogical(Params,TRIM(Txt),Found)
          IF (project) THEN
             IsValid=xios_is_valid_field(TRIM(FieldName)//"_elem")
             IF (.NOT.IsValid) &
             CALL FATAL(Caller,TRIM(FieldName)//"_elem is not a valid field")
          END IF

          Vari=Vari+1
          WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          FieldName = GetString( Params,TRIM(Txt),ScalarsExist)
        END DO

      END  SUBROUTINE SanityCheck

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Default field definition related to mesh information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      SUBROUTINE DefaultFieldDefinitions()
        IMPLICIT NONE
        TYPE(xios_fieldgroup) :: fieldgroup_hdl,fielddef_hdl
        TYPE(xios_field) :: field_hdl

        !!!  Default field_group mesh_info; operation="once"

        CALL xios_get_handle("field_definition",fielddef_hdl)
        ok=xios_is_valid_field("mesh_info")
        IF (.NOT.ok) THEN 
           CALL xios_add_child(fielddef_hdl,fieldgroup_hdl,"mesh_info")
           CALL xios_set_attr(fieldgroup_hdl,operation="once")

           ok=xios_is_valid_field("node_x")
           IF (.NOT.ok) THEN
             CALL xios_add_child(fieldgroup_hdl,field_hdl,"node_x")
             CALL xios_set_attr(field_hdl,name="x", &
                      standard_name="projection_x_coordinate",&
                      unit="m",grid_ref="GridNodes")
           END IF
           ok=xios_is_valid_field("node_y")
           IF (.NOT.ok) THEN
             CALL xios_add_child(fieldgroup_hdl,field_hdl,"node_y")
             CALL xios_set_attr(field_hdl,name="y", &
                      standard_name="projection_y_coordinate",&
                      unit="m",grid_ref="GridNodes")
           END IF
           ok=xios_is_valid_field("cell_area")
           IF (.NOT.ok) THEN
             CALL xios_add_child(fieldgroup_hdl,field_hdl,"cell_area")
             CALL xios_set_attr(field_hdl,name="cell_area",unit="m2",grid_ref="GridCells")
           END IF

           IF (SaveEdges) THEN
             ok=xios_is_valid_field("boundary_condition")
             IF (.NOT.ok) THEN
               CALL xios_add_child(fieldgroup_hdl,field_hdl,"boundary_condition")
               CALL xios_set_attr(field_hdl,name="boundary_condition",&
                     unit="1",default_value=0._dp,prec=4,grid_ref="GridEdges")
             END IF
           END IF
        END IF
       END SUBROUTINE DefaultFieldDefinitions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! set domain definitions: cells, nodes, edges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE DomainDefinitions()
        IMPLICIT NONE

        TYPE(Element_t),POINTER :: Element,Edge
        TYPE(Nodes_t),SAVE :: ElementNodes
        INTEGER, POINTER :: NodeIndexes(:)
        TYPE(GaussIntegrationPoints_t) :: IntegStuff
        REAL(KIND=dp) :: U,V,W,SqrtElementMetric
        LOGICAL :: stat

        REAL(KIND=dp) :: xg,yg
        REAL(KIND=dp) :: Lon,Lat
        REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: NodeLon,FaceLon,EdgeLon
        REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: NodeLat,FaceLat,EdgeLat
        REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE :: LonBnds,LatBnds
        REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE :: EdgeLonBnds,EdgeLatBnds
        INTEGER,DIMENSION(:,:),ALLOCATABLE :: Indexes


        INTEGER :: i,ii,t,n,k
        INTEGER :: M
        INTEGER :: imin
        INTEGER,DIMENSION(:),ALLOCATABLE :: Vertice
        INTEGER,DIMENSION(:),ALLOCATABLE :: GIndexes
        INTEGER,DIMENSION(:),ALLOCATABLE :: EdgeIndexes
        INTEGER, DIMENSION(:),ALLOCATABLE :: PartECount
        INTEGER :: begin
        INTEGER :: ProjCoord

        INTEGER :: ierr
        INTEGER :: TotalCount,nv
        LOGICAL :: FieldActive,ok
        LOGICAL :: DoProj

        ALLOCATE(Vertice(NumberOfActiveNodes),Node_x(NumberOfActiveNodes),Node_y(NumberOfActiveNodes))
        ALLOCATE(cell_area(NumberOfElements),Indexes(MaxElementNodes,NumberOfElements),GIndexes(NumberOfElements))
        ALLOCATE(NodeLon(NumberOfActiveNodes),NodeLat(NumberOfActiveNodes))
        ALLOCATE(FaceLon(NumberOfElements),FaceLat(NumberOfElements))
        IF (SaveEdges) THEN
          ALLOCATE(BoundaryCondition(NumberOfActiveEdges))
          ALLOCATE(EdgeIndexes(NumberOfActiveEdges))
          ALLOCATE(EdgeLon(NumberOfActiveEdges),EdgeLat(NumberOfActiveEdges))
          ALLOCATE(EdgeLonBnds(2,NumberOfActiveEdges),EdgeLatBnds(2,NumberOfActiveEdges))
        END IF
        ALLOCATE(LonBnds(MaxElementNodes,NumberOfElements),LatBnds(MaxElementNodes,NumberOfElements))

        Indexes=Connect_Fill

        ! Nodes
        n=0
        DO ii = 1, NumberOfDofNodes
          IF (NotOwnedNode(ii)) CYCLE
          n=n+1
          IF( NoPermutation ) THEN
            i = ii
          ELSE
            i = InvNodePerm(ii)
          END IF

          IF (Parallel) THEN
            Vertice(n) = Mesh % ParallelInfo % GlobalDOFs(i)
          ELSE
            Vertice(n) = i
          END IF
          Node_x(n) = Mesh%Nodes%x(i)
          Node_y(n) = Mesh%Nodes%y(i)

          CALL xy2LonLat(Node_x(n),Node_y(n),lon,lat)
          NodeLon(n) = lon
          NodeLat(n) = lat

        END DO

        imin=MINVAL(Vertice)
        IF (Parallel) THEN
                CALL MPI_ALLREDUCE(imin,nv,1,MPI_INTEGER,MPI_MIN,ELMER_COMM_WORLD,ierr)
                imin=nv
        ENDIF
        IF (imin.GT.1) Vertice=Vertice-imin+1

        ! Edges
        IF (SaveEdges) THEN
          t=0
          DO ii = 1, Mesh % NumberOfEdges
            IF (NotOwnedEdge(ii)) CYCLE
            t=t+1

            Edge => Mesh % Edges(ii)

            BoundaryCondition(t) = Edge % BoundaryInfo % Constraint

            n = Edge % TYPE % NumberOfNodes

            IF (n/=2) &
              CALL FATAL(Caller, 'Work only for edge element of type 202')

            IF (Parallel) THEN
              EdgeIndexes(t)=Edge % GElementIndex
            ELSE
              EdgeIndexes(t)=Edge % ElementIndex
            ENDIF

            NodeIndexes => Edge % NodeIndexes

            ! Edge center
            xg=SUM(Mesh%Nodes%x(NodeIndexes(1:n)))/n
            yg=SUM(Mesh%Nodes%y(NodeIndexes(1:n)))/n
            CALL xy2LonLat(xg,yg,Lon,Lat)
            EdgeLon(t) = Lon
            EdgeLat(t) = Lat

            ! Edge bounds
            DO k=1,n
              xg=Mesh%Nodes%x(NodeIndexes(k))
              yg=Mesh%Nodes%y(NodeIndexes(k))
              CALL xy2LonLat(xg,yg,Lon,Lat)
              EdgeLonBnds(k,t)=Lon
              EdgeLatBnds(k,t)=Lat
            END DO

          END DO
        END IF

        ! Elements
        ProjCoord = ListGetInteger(Solver % Values,"projection coordinate",DoProj)
        cell_area=0._dp
        t=0
        DO i = ElemFirst, ElemLast
          IF( .NOT. ActiveElem(i) ) CYCLE
          t=t+1

          Element => Model % Elements(i)
          IF (SaveLinear) THEN
            n = GetElementCorners( Element )
          ELSE
            n = GetElementNOFNodes(Element)
          END IF

          IF (Parallel) THEN
            GIndexes(t)=Element % GElementIndex
          ELSE
            GIndexes(t)=Element % ElementIndex
          ENDIF

          NodeIndexes => Element % NodeIndexes

          !! Clockwise ordering do not seems to be an issue for XIOS
          !IF(IsClockwise( Mesh%Nodes%x(NodeIndexes(1:n)),Mesh%Nodes%y(NodeIndexes(1:n)),n)) THEN
          !  CALL FATAL(Caller, &
          !        "Clock wise ordering ... implement reordering!")
          !END IF

          IF( NoPermutation ) THEN
            Indexes(1:n,t) = NodeIndexes(1:n)
          ELSE
            Indexes(1:n,t) = NodePerm( NodeIndexes(1:n) )
          END IF

          ! element center
          xg=SUM(Mesh%Nodes%x(NodeIndexes(1:n)))/n
          yg=SUM(Mesh%Nodes%y(NodeIndexes(1:n)))/n
          CALL xy2LonLat(xg,yg,Lon,Lat)
          FaceLon(t) = Lon
          FaceLat(t) = Lat

          DO k=1,n
            xg=Mesh%Nodes%x(NodeIndexes(k))
            yg=Mesh%Nodes%y(NodeIndexes(k))
            CALL xy2LonLat(xg,yg,Lon,Lat)
            LonBnds(k,t)=Lon
            LatBnds(k,t)=Lat
          END DO

          CALL GetElementNodes( ElementNodes, Element )
          IF (DoProj) THEN
             SELECT CASE (ProjCoord)
                CASE (1)
                        ElementNodes % x(:)=0._dp
                CASE (2)
                        ElementNodes % y(:)=0._dp
                CASE (3)
                        ElementNodes % z(:)=0._dp
                CASE DEFAULT
                        CALL FATAL(Caller, &
                            "Wrong projection coordinate"//I2S(ProjCoord))

             END SELECT
          END IF

          IntegStuff = GaussPoints( Element )
          Do k=1,IntegStuff % n
            U = IntegStuff % u(k)
            V = IntegStuff % v(k)
            W = IntegStuff % w(k)
            stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
                        Basis,dBasisdx )
            cell_area(t)=cell_area(t)+SqrtElementMetric*IntegStuff % s(k)
           End do

         END DO

         ! it seems that GIndex is not updated by the extrusion so it
         ! keeps the initial table........... should be ok then...

         IF (Parallel) THEN
           CALL MPI_ALLREDUCE(NumberOfElements,TotalCount,1,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
           CALL MPI_ALLREDUCE(MAXVAL(GIndexes),nv,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
           !! in principle we should have all the elements...
           IF (TotalCount.NE.nv) THEN
             CALL FATAL(Caller, "Pb with number of elements?? ")
           ENDIF

           CALL xios_set_domain_attr("cells",ni_glo=TotalCount, &
                ni=NumberOfElements,i_Index=GIndexes-1,&
                nvertex=MaxElementNodes , type='unstructured')

           CALL xios_set_domain_attr("cells",data_dim=1, data_ibegin=0, data_ni=NumberOfElements)
           CALL xios_set_domain_attr("cells",lonvalue_1d=FaceLon, latvalue_1d=FaceLat, &
                bounds_lon_1d=LonBnds,bounds_lat_1d=LatBnds)

           !! nodes overlap between partitions
           CALL MPI_ALLREDUCE(NumberOfActiveNodes,TotalCount,1,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
           CALL MPI_ALLREDUCE(MAXVAL(Vertice),nv,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
           !! in principle we should have all the nodes...
           IF (TotalCount.NE.nv) &
             CALL FATAL(Caller, "Pb with number of nodes?? ")


           CALL xios_set_domain_attr("nodes",ni_glo=TotalCount,&
                ni=NumberOfActiveNodes,i_Index=Vertice-1,&
                nvertex=1 , type='unstructured')

           CALL xios_set_domain_attr("nodes",data_dim=1, data_ibegin=0, data_ni=NumberOfActiveNodes)
           CALL xios_set_domain_attr("nodes",lonvalue_1d=NodeLon, latvalue_1d=NodeLat)

           !! edges
          IF (SaveEdges) THEN
           CALL MPI_ALLREDUCE(NumberOfActiveEdges,TotalCount,1,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
           CALL MPI_ALLREDUCE(MAXVAL(EdgeIndexes),nv,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)

           IF (TotalCount.NE.nv) &
             CALL FATAL(Caller, "Pb with number of edges?? ")
           CALL xios_set_domain_attr("edges",ni_glo=TotalCount,&
                ni=NumberOfActiveEdges,i_Index=EdgeIndexes-1,&
                nvertex=2 , type='unstructured')
            CALL xios_set_domain_attr("edges",data_dim=1, data_ibegin=0, data_ni=NumberOfActiveEdges)
            CALL xios_set_domain_attr("edges",lonvalue_1d=EdgeLon, latvalue_1d=EdgeLat, &
                bounds_lon_1d=EdgeLonBnds,bounds_lat_1d=EdgeLatBnds)
          ENDIF

         ELSE
           CALL xios_set_domain_attr("cells",ni_glo=NumberOfElements, ibegin=0, ni=NumberOfElements, &
                nvertex=MaxElementNodes , type='unstructured')
           CALL xios_set_domain_attr("cells",lonvalue_1d=FaceLon, latvalue_1d=FaceLat, &
                bounds_lon_1d=LonBnds,bounds_lat_1d=LatBnds)

           IF (SaveEdges) THEN
             CALL xios_set_domain_attr("edges",ni_glo=NumberOfActiveEdges, ibegin=0, ni=NumberOfActiveEdges, &
                nvertex=2 , type='unstructured')
             CALL xios_set_domain_attr("edges",lonvalue_1d=EdgeLon, latvalue_1d=EdgeLat, &
                bounds_lon_1d=EdgeLonBnds,bounds_lat_1d=EdgeLatBnds)
           END IF

           CALL xios_set_domain_attr("nodes",ni_glo=NumberOfDofNodes, ibegin=0, ni=NumberOfActiveNodes, &
                nvertex=1 , type='unstructured')
           CALL xios_set_domain_attr("nodes",lonvalue_1d=NodeLon, latvalue_1d=NodeLat)
         END IF
       
         !TO DO: 
! - can we send the connectivity tables directly

        DEALLOCATE(Vertice)
        DEALLOCATE(Indexes,GIndexes)
        DEALLOCATE(NodeLon,NodeLat)
        IF (SaveEdges) THEN
          DEALLOCATE(EdgeIndexes)
          DEALLOCATE(EdgeLon,EdgeLat)
          DEALLOCATE(EdgeLonBnds,EdgeLatBnds)
        END IF
        DEALLOCATE(FaceLon,FaceLat)
        DEALLOCATE(LonBnds,LatBnds)

      END SUBROUTINE DomainDefinitions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  Send requested mesh info variables
!!!!!!!!  maybe this could change and we should recompute those?? e.g. for active/passive BC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SendMeshVariables()
      IMPLICIT NONE
      LOGICAL :: FieldActive

     ! 
      FieldActive=xios_field_is_active("cell_area",.TRUE.)
      WRITE( Message,'(A,A,L2)') 'Xios request : ',"cell_area",FieldActive
      CALL Info(Caller,Message,Level=Olevel)
      IF (FieldActive.OR.AllwaysSend) THEN
          CALL xios_send_field("cell_area",cell_area)
      END IF
      FieldActive=xios_field_is_active("node_x",.TRUE.)
      WRITE( Message,'(A,A,L2)') 'Xios request : ',"node_x",FieldActive
      CALL Info(Caller,Message,Level=Olevel)
      IF (FieldActive.OR.AllwaysSend) THEN
          CALL xios_send_field("node_x",Node_x)
      END IF
      FieldActive=xios_field_is_active("node_y",.TRUE.)
      WRITE( Message,'(A,A,L2)') 'Xios request : ',"node_y",FieldActive
      CALL Info(Caller,Message,Level=Olevel)
      IF (FieldActive.OR.AllwaysSend) THEN
          CALL xios_send_field("node_y",Node_y)
      END IF
      IF (SaveEdges) THEN
        FieldActive=xios_field_is_active("boundary_condition",.TRUE.)
        WRITE( Message,'(A,A,L2)') 'Xios request : ',"boundary_condition",FieldActive
        CALL Info(Caller,Message,Level=Olevel)
        IF (FieldActive.OR.AllwaysSend) THEN
           CALL xios_send_field("boundary_condition",BoundaryCondition)
        ENDIF
      ENDIF
      END SUBROUTINE SendMeshVariables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Send requested global variables...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SendGlobalVariables()
        IMPLICIT NONE
        TYPE(Variable_t),POINTER :: Var
        CHARACTER(LEN=1024) :: Txt,GlobalFieldName
        INTEGER :: Vari
        LOGICAL ::ScalarsExist
        LOGICAL :: FieldActive

        GlobalFieldName = GetString( Params,'Global Variable 1',ScalarsExist)
        Vari=1
        DO WHILE (ScalarsExist)
          Var => VariableGet( Model % Mesh % Variables,&
             TRIM(GlobalFieldName),ThisOnly=.TRUE.,UnFoundFatal=.TRUE.)

          IF ((Var%TYPE .NE. Variable_global).AND.(SIZE(Var % Values).NE.Var % DOFs)) &
            CALL FATAL(Caller,"Variable was supposed to be global")

          FieldActive=xios_field_is_active(TRIM(Var%Name),.TRUE.)
          WRITE( Message,'(A,A,L2)') 'Xios request : ',TRIM(Var%Name),FieldActive
          CALL Info(Caller,Message,Level=Olevel)

          IF (FieldActive) THEN
            CALL xios_send_field(TRIM(Var%Name), Var % Values(1))
          ENDIF

          Vari=Vari+1
          WRITE(Txt,'(A,I0)') 'Global Variable ',Vari
          GlobalFieldName = GetString( Params,TRIM(Txt),ScalarsExist)
        END DO
 
      END SUBROUTINE SendGlobalVariables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Send requested scalar variables...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SendVariables()
        IMPLICIT NONE

        TYPE(Element_t),POINTER :: Element
        TYPE(Nodes_t),SAVE :: ElementNodes
        INTEGER, POINTER :: NodeIndexes(:)
        TYPE(GaussIntegrationPoints_t) :: IntegStuff
        REAL(KIND=dp) :: U,V,W,SqrtElementMetric
        REAL(KIND=dp) :: VarMean,area
        TYPE(Variable_t),POINTER :: Solution
        INTEGER :: VarType
        INTEGER :: ii,i,j,t,m,n,k
        LOGICAL :: stat
        LOGICAL :: Found

        REAL(KIND=dp),ALLOCATABLE :: NodeVar(:)
        REAL(KIND=dp),ALLOCATABLE :: EVar(:)

        INTEGER, POINTER :: Perm(:)
        REAL(KIND=dp),POINTER :: Values(:)

        CHARACTER(LEN=1024) :: Txt, ScalarFieldName
        INTEGER :: Vari
        LOGICAL :: ScalarsExist
        LOGICAL :: FieldActive
        LOGICAL :: Project

        ALLOCATE(NodeVar(NumberOfActiveNodes),EVar(NumberOfElements))

        ScalarFieldName = GetString( Params,'Scalar Field 1',ScalarsExist)
        Vari=1
        DO WHILE (ScalarsExist)
          Solution => VariableGet( Model % Mesh % Variables, TRIM(ScalarFieldName),ThisOnly=.FALSE.,UnFoundFatal=.TRUE.)
          VarType = Solution % TYPE

          IF (VarType == Variable_on_nodes) THEN
             WRITE(Txt,'(A,I0,A)') &
                'Scalar Field ',Vari,' compute cell average'
             Project = ListGetLogical(Params,TRIM(Txt),Found)
             IF (Project) THEN
                ! check if we need to send the data
                FieldActive=xios_field_is_active(TRIM(Solution%Name)//'_elem',.TRUE.)
                WRITE( Message,'(A,A,L2)') 'Xios request : ',TRIM(Solution%Name)//'_elem',FieldActive
                CALL Info(Caller,Message,Level=Olevel)
                IF (FieldActive.OR.AllwaysSend) THEN
                  Perm => Solution % Perm
                  Values => Solution % Values
                  t=0
                  DO i = ElemFirst, ElemLast
                    IF( .NOT. ActiveElem(i) ) CYCLE
                    t=t+1

                    Element => Model % Elements(i)
                    n = GetElementNOFNodes(Element)
                    NodeIndexes(1:n) => Element % NodeIndexes(1:n)

                    CALL GetElementNodes( ElementNodes, Element )

                    IntegStuff = GaussPoints( Element )
                    VarMean=0._dp
                    area=0._dp
                    Do k=1,IntegStuff % n
                      U = IntegStuff % u(k)
                      V = IntegStuff % v(k)
                      W = IntegStuff % w(k)
                      stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
                        Basis,dBasisdx )
                      area=area+SqrtElementMetric*IntegStuff % s(k)
                      VarMean=VarMean+ &
                        SqrtElementMetric*IntegStuff % s(k)*SUM(Basis(1:n)*Values(Perm(NodeIndexes(1:n))))
                    End do

                    EVar(t) = VarMean/area
                  END DO
                  WRITE( Message,'(A,A,L2)') 'Sending : ',TRIM(Solution%Name)//'_elem',FieldActive
                  CALL Info(Caller,Message,Level=Olevel)
                  CALL xios_send_field(TRIM(Solution%Name)//'_elem',EVar)

                END IF
             ENDIF
          END IF

          ! check if we need to send the data
          FieldActive=xios_field_is_active(TRIM(Solution%Name),.TRUE.)
          WRITE( Message,'(A,A,L2)') 'Xios request : ',TRIM(Solution%Name),FieldActive
          CALL Info(Caller,Message,Level=Olevel)

          IF (FieldActive.OR.AllwaysSend) THEN
           Perm => Solution % Perm
           Values => Solution % Values
           SELECT CASE (VarType)
            CASE (Variable_on_nodes)
               n=0
               DO ii=1,NumberOfDofNodes

                 IF (NotOwnedNode(ii)) CYCLE
                 n=n+1

                 IF( NoPermutation ) THEN
                   i = ii
                 ELSE
                   i = InvNodePerm(ii)
                 END IF
                 IF( ASSOCIATED( Perm ) ) THEN
                  j = Perm(i)
                 ELSE
                  j = i
                 END IF
                 IF (j==0) THEN
                  NodeVar(n) = 0._dp !should make it fillValue?
                 ELSE
                  NodeVar(n) = Values(j)
                 ENDIF
               END DO

               WRITE( Message,'(A,A,L2)') 'Sending : ',TRIM(Solution%Name),FieldActive
               CALL Info(Caller,Message,Level=Olevel)
               CALL xios_send_field(TRIM(Solution%Name),NodeVar)

            CASE (Variable_on_elements)
               t=0
               DO i = ElemFirst, ElemLast
                 IF( .NOT. ActiveElem(i) ) CYCLE
                 t=t+1

                 Element => Model % Elements(i)
                 m = Element % ElementIndex
                 IF( ASSOCIATED( Perm ) ) THEN
                  IF( m>SIZE( Perm ) ) THEN
                    j = 0
                  ELSE
                    j = Perm(m)
                  ENDIF
                 ELSE
                  j = m
                 END IF
                 IF (j==0) THEN
                  EVar(t) = 0._dp !should make it fillValue?
                 ELSE
                  EVar(t) = Values(j)
                 ENDIF
               END DO

               WRITE( Message,'(A,A,L2)') 'Sending : ',TRIM(Solution%Name),FieldActive
               CALL Info(Caller,Message,Level=Olevel)
               CALL xios_send_field(TRIM(Solution%Name),EVar)
            CASE DEFAULT
               CALL FATAL(Caller,"sorry don't know how to deal with this variable type")
           END SELECT

           END IF

          Vari=Vari+1
          WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          ScalarFieldName = GetString( Params,TRIM(Txt),ScalarsExist)
        END DO

        DEALLOCATE(NodeVar,EVar)

      END SUBROUTINE SendVariables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  Check if ordering is clockwise; not sure this is really needed by xios?
!!!!!!!!   but UGRID convention requires anti-clockwise ordering...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION IsClockwise(x,y,n) RESULT(clokwise)
      LOGICAL :: clokwise
      REAL(KIND=dp) :: x(n),y(n)
      INTEGER :: n
      INTEGER :: i,ind2
      REAL(KIND=dp) :: sarea

      sarea=0._dp
      Do i=1,n
        ind2=mod(i,n)+1
        sarea=sarea+(x(ind2)-x(i))*(y(ind2)+y(i))
      End do
      clokwise=(sarea.GT.0)

      END FUNCTION IsClockwise

      END SUBROUTINE XIOSOutputSolver
