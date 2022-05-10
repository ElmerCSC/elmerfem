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
        ! current limitations only bulk 
        CALL ListAddNewLogical( SolverParams,'Save Bulk Only',.TRUE.)
        ! need to skip halo
        CALL ListAddNewLogical( SolverParams,'Skip Halo Elements',.TRUE.)

        ! Need edge elements to compute the edge table
        CALL ListAddNewString( SolverParams,"Element","n:0 e:1")

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
      INTEGER :: i,M,t
      CHARACTER(LEN=2) :: strg_var
      LOGICAL :: ierr


      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(ValueList_t),POINTER :: Params
      LOGICAL :: GotIt
      LOGICAL :: Parallel
      INTEGER :: GroupId
      INTEGER :: PEs,Part 
      REAL(KIND=dp) :: SimTime

      LOGICAL :: SaveLinear=.TRUE. ! save only corners
      INTEGER :: MaxElementNodes  ! // MaxNumNodesPerFace
      INTEGER,PARAMETER :: Connect_Fill=-1
      LOGICAL :: ok

      TYPE(xios_duration) :: dtime,time_units
      TYPE(xios_context) :: ctx_hdl

      Params => GetSolverParams()
      Mesh => Model % Mesh
      ! Comment: Can we get the Mas number of corners instead of nodes?
      ! we should loop over saved elements....??
      MaxElementNodes = Model % MaxElementNodes

      IF (Mesh % MeshDim.NE.2) &
        CALL FATAL(Caller,"Mesh dim should be 2 for now...")

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
        CALL GenerateSavePermutation(Mesh,.FALSE.,.FALSE.,SaveLinear,ActiveElem,NumberOfGeomNodes,&
               NoPermutation,NumberOfDofNodes,DgPerm,InvDgPerm,NodePerm,InvNodePerm)

        !------------------------------------------------------------------------------
        ! Parallel case we will exclude nodes not owned by partition...
        !------------------------------------------------------------------------------
        ALLOCATE(NotOwnedNode(NumberOfDofNodes))
        NotOwnedNode=.FALSE.
        NumberOfActiveNodes=NumberOfDofNodes
        IF (Parallel) THEN
          DO i=1,NumberOfDofNodes
            IF (Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1).NE.ParEnv % MyPE) THEN
              NotOwnedNode(i)=.TRUE.
              NumberOfActiveNodes=NumberOfActiveNodes-1
            END IF
          END DO
        END IF

        !------------------------------------------------------------------------------
        ! Edges
        !------------------------------------------------------------------------------
        ALLOCATE(NotOwnedEdge(Mesh % NumberOfEdges))
        NotOwnedEdge=.FALSE.

        NumberOfActiveEdges=Mesh % NumberOfEdges
        IF (Parallel) THEN
          DO i=1,Mesh % NumberOfEdges
           IF (Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours(1).NE.ParEnv % MyPE) THEN
              NotOwnedEdge(i)=.TRUE.
              NumberOfActiveEdges=NumberOfActiveEdges-1
           ENDIF
          END DO
        END IF


        ! The partition is active for saving if there are any nodes
        ! to write. There can be no elements nor dofs without nodes.
        CALL ParallelActive( NumberOfDofNodes > 0 )

        IF( ElemLast > Mesh % NumberOfBulkElements ) &
          CALL FATAL(Caller, &
                "Saving boundary elements not supported yet....")

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! XIOS context definition
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL xios_get_handle(TRIM(xios_id),ctx_hdl)
        CALL xios_set_current_context(ctx_hdl)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! check that mandatory fields have been defined
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL SanityCheck()

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Set-up the time step from elmer....
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ierr=xios_getvar("time_units",strg_var)
        IF (.NOT. ierr) THEN
          CALL FATAL(Caller, &
                 "<time_units> variable not found")
        ELSE
          time_units=xios_duration_convert_from_string(TRIM(strg_var)) 
        ENDIF
        dtime = dt*time_units
        CALL xios_set_timestep(dtime)

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
        ! maybe it could be initialised in this routine too instead of
        ! in SParIterComm?
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL xios_close_context_definition

        CALL Info(Caller, 'Context initialisation done')

      END IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Update calendar
      !  should we add a sanity check that dt has not changed?
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL xios_update_calendar(nTime)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! send requested variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL SendMeshVariables()
      CALL SendGlobalVariables()
      CALL SendVariables()

      CONTAINS 

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

         IsValid=xios_is_valid_grid("GridEdges")
         IF (.NOT.IsValid) &
                 CALL FATAL(Caller,"<GridEdges> grid not defined")

         ! domains
         IsValid=xios_is_valid_domain("cells")
         IF (.NOT.IsValid) &
                 CALL FATAL(Caller,"<cells> domain not defined")

         IsValid=xios_is_valid_domain("nodes")
         IF (.NOT.IsValid) &
                 CALL FATAL(Caller,"<nodes> domain not defined")

         IsValid=xios_is_valid_domain("edges")
         IF (.NOT.IsValid) &
                 CALL FATAL(Caller,"<edges> domain not defined")

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
          WRITE(Txt,'(A,I0)') 'Scalar Variable ',Vari
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

           ok=xios_is_valid_field("boundary_condition")
           IF (.NOT.ok) THEN
             CALL xios_add_child(fieldgroup_hdl,field_hdl,"boundary_condition")
             CALL xios_set_attr(field_hdl,name="boundary_condition",&
                     unit="1",default_value=0._dp,prec=4,grid_ref="GridEdges")
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
        INTEGER,DIMENSION(:),ALLOCATABLE :: Vertice
        INTEGER,DIMENSION(:),ALLOCATABLE :: GIndexes
        INTEGER,DIMENSION(:),ALLOCATABLE :: EdgeIndexes
        INTEGER, DIMENSION(:),ALLOCATABLE :: PartECount
        INTEGER :: begin

        INTEGER :: ierr
        INTEGER :: TotalCount,nv
        LOGICAL :: FieldActive,ok

        ALLOCATE(BoundaryCondition(NumberOfActiveEdges))
        ALLOCATE(Vertice(NumberOfActiveNodes),Node_x(NumberOfActiveNodes),Node_y(NumberOfActiveNodes))
        ALLOCATE(cell_area(NumberOfElements),Indexes(MaxElementNodes,NumberOfElements),GIndexes(NumberOfElements))
        ALLOCATE(EdgeIndexes(NumberOfActiveEdges))
        ALLOCATE(NodeLon(NumberOfActiveNodes),NodeLat(NumberOfActiveNodes))
        ALLOCATE(FaceLon(NumberOfElements),FaceLat(NumberOfElements))
        ALLOCATE(EdgeLon(NumberOfActiveEdges),EdgeLat(NumberOfActiveEdges))
        ALLOCATE(EdgeLonBnds(2,NumberOfActiveEdges),EdgeLatBnds(2,NumberOfActiveEdges))
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

        ! Edges
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

        ! Elements
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
          IF(IsClockwise( Mesh%Nodes%x(NodeIndexes(1:n)),Mesh%Nodes%y(NodeIndexes(1:n)),n)) THEN
            CALL FATAL(Caller, &
                  "Clock wise ordering ... implement reordering!")
          END IF

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

         IF (Parallel) THEN
           CALL MPI_ALLREDUCE(NumberOfElements,TotalCount,1,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)

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

         ELSE
           CALL xios_set_domain_attr("cells",ni_glo=NumberOfElements, ibegin=0, ni=NumberOfElements, &
                nvertex=MaxElementNodes , type='unstructured')
           CALL xios_set_domain_attr("cells",lonvalue_1d=FaceLon, latvalue_1d=FaceLat, &
                bounds_lon_1d=LonBnds,bounds_lat_1d=LatBnds)

           CALL xios_set_domain_attr("edges",ni_glo=NumberOfActiveEdges, ibegin=0, ni=NumberOfActiveEdges, &
                nvertex=2 , type='unstructured')
           CALL xios_set_domain_attr("edges",lonvalue_1d=EdgeLon, latvalue_1d=EdgeLat, &
                bounds_lon_1d=EdgeLonBnds,bounds_lat_1d=EdgeLatBnds)

           CALL xios_set_domain_attr("nodes",ni_glo=NumberOfDofNodes, ibegin=0, ni=NumberOfActiveNodes, &
                nvertex=1 , type='unstructured')
           CALL xios_set_domain_attr("nodes",lonvalue_1d=NodeLon, latvalue_1d=NodeLat)
         END IF
       
         !TO DO: 
! - can we send the connectivity tables directly

        DEALLOCATE(Vertice)
        DEALLOCATE(Indexes,GIndexes,EdgeIndexes)
        DEALLOCATE(NodeLon,NodeLat)
        DEALLOCATE(EdgeLon,EdgeLat)
        DEALLOCATE(EdgeLonBnds,EdgeLatBnds)
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
      CALL Info(Caller,Message,Level=4)
      IF (FieldActive) THEN
          CALL xios_send_field("cell_area",cell_area)
      END IF
      FieldActive=xios_field_is_active("node_x",.TRUE.)
      WRITE( Message,'(A,A,L2)') 'Xios request : ',"node_x",FieldActive
      CALL Info(Caller,Message,Level=4)
      IF (FieldActive) THEN
          CALL xios_send_field("node_x",Node_x)
      END IF
      FieldActive=xios_field_is_active("node_y",.TRUE.)
      WRITE( Message,'(A,A,L2)') 'Xios request : ',"node_y",FieldActive
      CALL Info(Caller,Message,Level=4)
      IF (FieldActive) THEN
          CALL xios_send_field("node_y",Node_y)
      END IF
      FieldActive=xios_field_is_active("boundary_condition",.TRUE.)
      WRITE( Message,'(A,A,L2)') 'Xios request : ',"boundary_condition",FieldActive
      IF (FieldActive) THEN
         CALL xios_send_field("boundary_condition",BoundaryCondition)
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

        GlobalFieldName = GetString( Params,'Global Variable 1',ScalarsExist)
        Vari=1
        DO WHILE (ScalarsExist)
          Var => VariableGet( Model % Mesh % Variables,&
             TRIM(GlobalFieldName),ThisOnly=.TRUE.,UnFoundFatal=.TRUE.)

          IF ((Var%TYPE .NE. Variable_global).AND.(SIZE(Var % Values).NE.Var % DOFs)) &
            CALL FATAL(Caller,"Variable was supposed to be global")

          CALL xios_send_field(TRIM(Var%Name), Var % Values(1))

          Vari=Vari+1
          WRITE(Txt,'(A,I0)') 'Global Variable ',Vari
          GlobalFieldName = GetString( Params,TRIM(Txt),ScalarsExist)
        END DO
 
      END SUBROUTINE SendGlobalVariables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Send srequested scalar variables...
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
          Solution => VariableGet( Model % Mesh % Variables, TRIM(ScalarFieldName),ThisOnly=.TRUE.,UnFoundFatal=.TRUE.)
          VarType = Solution % TYPE

          IF (VarType == Variable_on_nodes) THEN
             WRITE(Txt,'(A,I0,A)') &
                'Scalar Field ',Vari,' compute cell average'
             Project = ListGetLogical(Params,TRIM(Txt),Found)
             IF (Project) THEN
                ! check if wee need to send the data
                FieldActive=xios_field_is_active(TRIM(Solution%Name)//'_elem',.TRUE.)
                WRITE( Message,'(A,A,L2)') 'Xios request : ',TRIM(Solution%Name)//'_elem',FieldActive
                CALL Info(Caller,Message,Level=4)
                IF (FieldActive) THEN
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
                  CALL xios_send_field(TRIM(Solution%Name)//'_elem',EVar)

                END IF
             ENDIF
          END IF

          ! check if wee need to send the data
          FieldActive=xios_field_is_active(TRIM(Solution%Name),.TRUE.)
          WRITE( Message,'(A,A,L2)') 'Xios request : ',TRIM(Solution%Name),FieldActive
          CALL Info(Caller,Message,Level=4)

          IF (FieldActive) THEN
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
