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
      LOGICAL, ALLOCATABLE,SAVE :: ActiveElem(:),NotOwnedNode(:)
      REAL(KIND=dp),DIMENSION(:),ALLOCATABLE,SAVE :: cell_area
      REAL(KIND=dp),ALLOCATABLE,SAVE :: Basis(:), dBasisdx(:,:)
      INTEGER,SAVE :: NumberOfGeomNodes,NumberOfDofNodes,NumberOfElements
      INTEGER,SAVE :: NumberOfActiveNodes
      LOGICAL,SAVE :: NoPermutation
      INTEGER :: i,M
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
      LOGICAL :: FieldActive

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
            IF (Solver % Matrix % ParallelInfo % NeighbourList(i) % Neighbours(1).NE.ParEnv % MyPE) THEN
              NotOwnedNode(i)=.TRUE.
              NumberOfActiveNodes=NumberOfActiveNodes-1
            END IF
          END DO
        END IF


        ! The partition is active for saving if there are any nodes
        ! to write. There can be no elements nor dofs without nodes.
        CALL ParallelActive( NumberOfDofNodes > 0 )

        IF( ElemLast > Mesh % NumberOfBulkElements ) &
          CALL FATAL(Caller, &
                "Saving boundary elements not supported yet....")

        CALL xios_get_handle(TRIM(xios_id),ctx_hdl)
        CALL xios_set_current_context(ctx_hdl)

        !! Set-up the time step from elmer....
        ierr=xios_getvar("time_units",strg_var)
        IF (.NOT. ierr) THEN
          CALL FATAL(Caller, &
                 "<time_units> variable not found")
        ELSE
          time_units=xios_duration_convert_from_string(TRIM(strg_var)) 
        ENDIF
        dtime = dt*time_units
        CALL xios_set_timestep(dtime)

        CALL SendMeshInfo()

        CALL Info(Caller, 'SendMeshInfo Done')


      END IF

      CALL xios_update_calendar(nTime)

     ! does not same to work if I send it before calendar update
      FieldActive=xios_field_is_active("cell_area",.TRUE.)
      WRITE( Message,'(A,A,L2)') 'Xios request : ',"cell_area",FieldActive
      CALL Info(Caller,Message,Level=4)
      IF (FieldActive) THEN
          CALL xios_send_field("cell_area",cell_area)
      END IF

      CALL SendGlobalVariables()
      CALL SendVariables()

      CONTAINS 

      SUBROUTINE SendMeshInfo()
        IMPLICIT NONE

        TYPE(Element_t),POINTER :: Element
        TYPE(Nodes_t),SAVE :: ElementNodes
        INTEGER, POINTER :: NodeIndexes(:)
        TYPE(GaussIntegrationPoints_t) :: IntegStuff
        REAL(KIND=dp) :: U,V,W,SqrtElementMetric
        LOGICAL :: stat

        REAL(KIND=dp) :: xg,yg
        REAL(KIND=dp) :: Lon,Lat
        REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: x,y
        REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: NodeLon,FaceLon
        REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: NodeLat,FaceLat
        REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE :: LonBnds,LatBnds
        INTEGER,DIMENSION(:,:),ALLOCATABLE :: Indexes


        INTEGER :: i,ii,t,n,k
        INTEGER :: M
        INTEGER,DIMENSION(:),ALLOCATABLE :: Vertice
        INTEGER,DIMENSION(:),ALLOCATABLE :: GIndexes
        INTEGER, DIMENSION(:),ALLOCATABLE :: PartECount
        INTEGER :: begin

        INTEGER :: ierr
        INTEGER :: TotalCount,nv
        LOGICAL :: FieldActive

        ALLOCATE(Vertice(NumberOfActiveNodes),x(NumberOfActiveNodes),y(NumberOfActiveNodes))
        ALLOCATE(cell_area(NumberOfElements),Indexes(MaxElementNodes,NumberOfElements),GIndexes(NumberOfElements))
        ALLOCATE(NodeLon(NumberOfActiveNodes),NodeLat(NumberOfActiveNodes))
        ALLOCATE(FaceLon(NumberOfElements),FaceLat(NumberOfElements))
        ALLOCATE(LonBnds(MaxElementNodes,NumberOfElements),LatBnds(MaxElementNodes,NumberOfElements))

        Indexes=Connect_Fill

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
          x(n) = Mesh%Nodes%x(i)
          y(n) = Mesh%Nodes%y(i)

          CALL xy2LonLat(x(n),y(n),lon,lat)
          NodeLon(n) = lon
          NodeLat(n) = lat

        END DO

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

           !! nodes overlap betwwen partitions
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

         ELSE
           CALL xios_set_domain_attr("cells",ni_glo=NumberOfElements, ibegin=0, ni=NumberOfElements, &
                nvertex=MaxElementNodes , type='unstructured')
           CALL xios_set_domain_attr("cells",lonvalue_1d=FaceLon, latvalue_1d=FaceLat, &
                bounds_lon_1d=LonBnds,bounds_lat_1d=LatBnds)

           CALL xios_set_domain_attr("nodes",ni_glo=NumberOfDofNodes, ibegin=0, ni=NumberOfActiveNodes, &
                nvertex=1 , type='unstructured')
           CALL xios_set_domain_attr("nodes",lonvalue_1d=NodeLon, latvalue_1d=NodeLat)
         END IF


        CALL xios_close_context_definition
        

         !TO DO: 
! - can we send the connectivity tables?
! - do edges
! - send the projected coordinates?

        DEALLOCATE(Vertice,x,y)
        DEALLOCATE(Indexes,GIndexes)
        DEALLOCATE(NodeLon,NodeLat)
        DEALLOCATE(FaceLon,FaceLat)
        DEALLOCATE(LonBnds,LatBnds)

      END SUBROUTINE SendMeshInfo

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
