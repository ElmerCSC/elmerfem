!*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Apr 2001
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Mesh partitioning utilities including interfaces to Zoltan
!------------------------------------------------------------------------------

#include "../config.h"

MODULE MeshPartition

  USE Types
  USE MeshUtils

#ifdef HAVE_ZOLTAN
  USE Zoltan
#endif

  IMPLICIT NONE

  TYPE MeshPack_t
     INTEGER :: NumberOfNodes, NumberOfBulkElements, NumberOfBoundaryElements
     LOGICAL, ALLOCATABLE :: NodeMask(:)
     INTEGER :: icount ! position counter for integer values
     INTEGER :: rcount ! position counter for real values
     INTEGER :: bcpos  ! offset for boundary element data
     INTEGER :: indpos ! offset for node index data
     INTEGER :: lcount ! position counter for logical values
     INTEGER, ALLOCATABLE :: idata(:)       ! integer data
     REAL(KIND=dp), ALLOCATABLE :: rdata(:) ! real data
     LOGICAL, ALLOCATABLE :: ldata(:)       ! logical data
   END TYPE MeshPack_t


CONTAINS

  !============================================
  !============================================
  !          ZOLTAN SUBROUTINES
  !============================================
  !============================================

  !> Interface to Zoltan parallel (re)partitioner - returns the new partition
  !> info in Mesh % Repartition (defined on elements)
  !> Dual-graph (element connectivity) is determined based on shared faces(3D)/edges(2D)
  !-------------------------------------------------------------------------------------
  SUBROUTINE Zoltan_Interface( Model, Mesh, SerialMode )

    USE MeshUtils

#ifdef HAVE_ZOLTAN
    USE Zoltan
#endif

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, OPTIONAL :: SerialMode
    !------------------------

#ifdef HAVE_ZOLTAN
    TYPE(Element_t), POINTER :: Element
    TYPE(Graph_t) :: LocalGraph
    REAL(KIND=dp) :: t1,t2
    INTEGER :: i,j,k,l,m,n,ierr,NNodes,NBulk,counter,DIM,&
         max_elemno,NoPart
    INTEGER, ALLOCATABLE :: ElemAdj(:), ElemStart(:), ElemAdjProc(:), ParElemAdj(:), ParElemStart(:),&
         ParElemIdx(:),ParElemAdjProc(:),sharecount(:),&
         ParElemMap(:)

    CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="Zoltan_Interface"

    !Zoltan things
    TYPE(Zoltan_Struct), POINTER :: zz_obj
    INTEGER(Zoltan_INT) :: zierr,numGIDEntries, numLidEntries, numImport, numExport
    INTEGER(Zoltan_INT),DIMENSION(:), POINTER :: importGlobalGids, importLocalGids, importProcs, &
         importToPart,exportGlobalGids,exportLocalGids, exportProcs, exportToPart
    REAL(Zoltan_FLOAT) :: version
    LOGICAL :: changes,Debug,Serial,Found

    TYPE ElemTable_t
       INTEGER :: counter=0
       INTEGER, ALLOCATABLE :: Idx(:)
    END TYPE ElemTable_T
    TYPE(ElemTable_t),ALLOCATABLE :: NodeElems(:),ElemElems(:)

    TYPE(ValueList_t), POINTER :: PartParams
    TYPE(ValueListEntry_t), POINTER :: ptr
    INTEGER :: ncopy
 
    
    CALL Info(FuncName,'Calling Zoltan for mesh partitioning',Level=8)
    PartParams => Model % Simulation

    IF( PRESENT( SerialMode ) ) THEN
      Serial = SerialMode
    ELSE
      Serial = ( ParEnv % PEs == 1 )
    END IF
    
    NoPart = ParEnv % PEs
    IF( NoPart == 1 ) THEN
      NoPart = ListGetInteger( PartParams,'Number Of Partitions',Found ) 
      IF( NoPart <= 1 ) THEN
        CALL Info(FuncName,'Nothing to do without any partitions requested!')
        RETURN
      END IF
    END IF
          
    NNodes = Mesh % NumberOfNodes
    NBulk = Mesh % NumberOfBulkElements
    DIM = CoordinateSystemDimension()

    IF( dim == 0 ) dim = Mesh % MeshDim
    
    zierr = Zoltan_Initialize(version)
    IF(zierr /= 0) CALL Fatal(FuncName,"Unable to initialize Zoltan partitioner")

    NULLIFY(zz_obj)

    ! Initialize zoltan for partitioning
    IF( Serial ) THEN
      zz_obj => Zoltan_Create(MPI_COMM_SELF)
    ELSE
      zz_obj => Zoltan_Create(ELMER_COMM_WORLD)
    END IF

    ! Set default values for keywords, if not given
    CALL ListAddNewString( PartParams,"zoltan: debug_level","0")
    CALL ListAddNewString( PartParams,"zoltan: lb_method","graph")
    CALL ListAddNewString( PartParams,"zoltan: graph_package","phg")
    CALL ListAddNewString( PartParams,"zoltan: num_gid_entries","1")
    CALL ListAddNewString( PartParams,"zoltan: num_lid_entries","1")
    CALL ListAddNewString( PartParams,"zoltan: obj_weight_dim","0")
    CALL ListAddNewString( PartParams,"zoltan: edge_weight_dim","0")
    CALL ListAddNewString( PartParams,"zoltan: check_graph","0")
    CALL ListAddNewString( PartParams,"zoltan: phg_multilevel","1")

    ! The settings for serial vs. parallel operation differ slightly
    IF( Serial ) THEN
      CALL ListAddNewString( PartParams,"zoltan: return_lists","export part")    
      CALL ListAddNewString( PartParams,"zoltan: lb_approach","partition")  
      CALL ListAddNewString( PartParams,"zoltan: num_global_parts",TRIM(I2S(NoPart)))  
    ELSE
      CALL ListAddNewString( PartParams,"zoltan: return_lists","all")    !TODO - we only use export list
      CALL ListAddNewString( PartParams,"zoltan: lb_approach","refine")  !repartition/refine <- faster
    END IF
      
    ! Pass keyword with prefix 'zoltan:' from the value list to zoltan
    Ptr => PartParams % Head
    ncopy = 0
    DO WHILE( ASSOCIATED(ptr) )
      n = ptr % NameLen
      k = 7 ! as for 'zoltan:'
      IF( n > k ) THEN
        IF( ptr % Name(1:k) == 'zoltan:' ) THEN
          l = k+1
          ! Remove the extra blanco after prefix if present
          DO WHILE( ptr % Name(l:l) == ' ')
            l = l+1
          END DO

          zierr = Zoltan_Set_Param(zz_obj,ptr % Name(l:n),ptr % Cvalue )
          IF(zierr /= 0) THEN
            CALL Fatal(FuncName,'Unable to set Zoltan Parameter: '//TRIM(ptr % Name(l:n)))
          ELSE
            CALL Info(FuncName,'Succesfully set Zoltan parameter: '&
                //TRIM(ptr % Name(l:n))//' to '//TRIM(ptr % CValue),Level=8)
          END IF

          CALL Info(FuncName,'Transferred prefix keyword to zoltan: '//TRIM(ptr % Name(l:n)),Level=12)
          ncopy = ncopy + 1
        END IF
      END IF
      ptr => ptr % Next
    END DO
    IF( ncopy > 0 ) THEN
      CALL Info(FuncName,'Succefully set '//TRIM(I2S(ncopy))//' keywords in zoltan library',Level=6)
    END IF
        
    !Callback functions to query number of elements and the element data
    zierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_OBJ_FN_TYPE,zoltNumObjs)
    IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan element count callback.")
    zierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_OBJ_LIST_FN_TYPE,zoltGetObjs)
    IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan element info callback.")

    !Callback functions to query number of edges and the edge data
    zierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_EDGES_FN_TYPE,zoltNumEdges)
    IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Callback Function: zoltNumEdges")
    zierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_EDGE_LIST_FN_TYPE,zoltGetEdgeList)
    IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Callback Function: zoltGetEdgeList")


    !Parameters to set to get graph repartitioning (without relying on libparmetis?!)
    !LB_APPROACH = REPARTITION (OR REFINE)
    !LB_METHOD = GRAPH

    !Need the following functions for repartitioning:

    ! ZOLTAN_NUM_OBJ_FN
    ! ZOLTAN_OBJ_LIST_FN
    ! ZOLTAN_NUM_EDGES_MULTI_FN or ZOLTAN_NUM_EDGES_FN 
    ! ZOLTAN_EDGE_LIST_MULTI_FN or ZOLTAN_EDGE_LIST_FN

    ! ZOLTAN_OBJ_SIZE_MULTI_FN or ZOLTAN_OBJ_SIZE_FN  - Optional for LB_APPROACH=Repartition.
    ! ZOLTAN_PART_MULTI_FN or ZOLTAN_PART_FN - Optional for LB_APPROACH=Repartition and for REMAP=1. 

    IF( Serial ) THEN          
      CALL LocalElemAdjacency( Mesh, ElemAdj, ElemAdjProc, ElemStart, DIM )
    ELSE
      CALL GlobalElemAdjacency( Mesh, ElemAdj, ElemAdjProc, ElemStart, DIM )
    END IF
      
    numGidEntries = 1
    numLidEntries = 1

    CALL Info(FuncName,'Going into Zoltan partitioning',Level=12)
    zierr = Zoltan_LB_Partition(zz_obj, changes, numGidEntries, numLidEntries, &
         numImport, importGlobalGids, importLocalGids, importProcs, importToPart, &
         numExport, exportGlobalGids, exportLocalGids, exportProcs, exportToPart)
    IF(zierr /= 0) CALL Fatal(FuncName,"Error computing partitioning in Zoltan")

    
    ! Put the information in Mesh % Repartition - boundary elems will follow bulks (thanks Peter!)
    IF(ASSOCIATED(Mesh % Repartition)) DEALLOCATE(Mesh % Repartition)
    ALLOCATE(Mesh % Repartition(NBulk))
    Mesh % Repartition = ParEnv % MyPE + 1 !default stay on this proc
    
    DO i=1,numExport
      IF(exportLocalGids(i) > NBulk .OR. exportLocalGids(i) <= 0) &
          CALL Fatal(FuncName, "Bad local ID")
      IF( Serial ) THEN        
        Mesh % Repartition(exportLocalGids(i)) = exportToPart(i) + 1
      ELSE
        Mesh % Repartition(exportLocalGids(i)) = exportProcs(i) + 1
      END IF
    END DO
      
  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! User defined query function to register with Zoltan
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER FUNCTION zoltNumObjs(DATA, ierr)
      use zoltan 
      implicit none

      ! Local declarations
      INTEGER(Zoltan_INT), INTENT(in) :: DATA(*)
      INTEGER(ZOLTAN_INT), INTENT(out) :: ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      zoltNumObjs = NBulk
      ierr = ZOLTAN_OK
      !PRINT *,'zoltNumObjs:',ParEnv % MyPE, zoltNumObjs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    END FUNCTION zoltNumObjs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! User defined query function to register with Zoltan
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE zoltGetObjs (DATA, num_gid_entries, num_lid_entries, global_ids, & 
         local_ids, wgt_dim, obj_wgts, ierr)
      use zoltan
      implicit none

      INTEGER(ZOLTAN_INT), INTENT(in) :: DATA(*)
      !TYPE(Mesh_t), POINTER, INTENT(in) :: DATA
      ! TYPE(Zoltan_User_Data_1) :: DATA
      integer(ZOLTAN_INT), intent(in) :: num_gid_entries 
      integer(ZOLTAN_INT), intent(in) ::  num_lid_entries
      integer(ZOLTAN_INT), intent(out) :: global_ids(*)
      integer(ZOLTAN_INT), intent(out) :: local_ids(*)
      integer(ZOLTAN_INT), intent(in) :: wgt_dim 
      real(ZOLTAN_FLOAT), intent(out) :: obj_wgts(*)
      integer(ZOLTAN_INT), intent(out) :: ierr

      ! local declarations
      integer :: i
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i= 1, NBulk
!        global_ids(i) = Mesh % Elements(i) % GElementIndex
        global_ids(i) = Mesh % Elements(i) % ElementIndex
        local_ids(i) = i
      end do

      !PRINT *,'zoltGetObjs:',ParEnv % MyPe, NBulk, MINVAL( local_ids(1:Nbulk)), MAXVAL( local_ids(1:NBulk))
      
      ierr = ZOLTAN_OK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    END SUBROUTINE zoltGetObjs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    INTEGER FUNCTION zoltNumEdges(DATA, num_gid_entries, num_lid_entries, global_id, local_id, ierr)
      INTEGER(ZOLTAN_INT), INTENT(in) :: DATA  
      INTEGER(Zoltan_INT), INTENT(IN) :: num_gid_entries, num_lid_entries  
      INTEGER(Zoltan_INT), INTENT(IN) :: global_id  
      INTEGER(Zoltan_INT), INTENT(IN) :: local_id  
      INTEGER(Zoltan_INT), INTENT(OUT) :: ierr

      zoltNumEdges = (ElemStart(local_id+1) - ElemStart(local_id))

      !PRINT *,'zoltNumbEdges:',parenv % mype,local_id,ElemStart(local_id),zoltNumEdges
      ierr = 0
    END FUNCTION zoltNumEdges

    SUBROUTINE zoltGetEdgeList(DATA, num_gid_entries, num_lid_entries, global_id, local_id, &
         nbor_global_id, nbor_procs, wgt_dim, ewgts, ierr)

      !TYPE(Mesh_t), POINTER, INTENT(in) :: DATA
      ! TYPE(Zoltan_User_Data_1) :: DATA
      INTEGER(ZOLTAN_INT), INTENT(in) :: DATA  
      INTEGER(Zoltan_INT), INTENT(IN) :: num_gid_entries, num_lid_entries  
      INTEGER(Zoltan_INT), INTENT(IN) :: global_id
      INTEGER(Zoltan_INT), INTENT(IN) :: local_id
      INTEGER(Zoltan_INT), INTENT(OUT) :: nbor_global_id(*)
      INTEGER(Zoltan_INT), INTENT(OUT) :: nbor_procs(*)
      INTEGER(Zoltan_INT), INTENT(IN) :: wgt_dim(*)
      REAL(Zoltan_FLOAT), INTENT(OUT) :: ewgts(*)
      INTEGER(Zoltan_INT), INTENT(OUT) :: ierr 
      !----------------
      INTEGER :: counter,i,k,nlocal,nother

      nlocal = (ElemStart(local_id+1) - ElemStart(local_id))

      DO i=1,nlocal
        k = i + ElemStart(local_id) - 1

        nbor_global_id(i) = ElemAdj(k)
        nbor_procs(i) = ElemAdjProc(k)
      END DO

      !PRINT *,'zoltGetEdgeList:',parenv % mype, local_id, nlocal
            
    END SUBROUTINE ZoltGetEdgeList

#else
    CALL FATAL('Zoltan_Interface',&
         'Repartitioning utility Zoltan (Trilinos) has not been installed')
#endif
  END SUBROUTINE Zoltan_Interface

  !>Returns a CRS dual graph of element face (3D) or edge (2D) connections, including
  !>across partitions. ElemAdj contains the global element numbers of connected elements,
  !>ElemStart describes the CRS positions of each elem, and ElemAdjProc contains the partition
  !>of the connected element.
  !---------------------------------------------------------------------------------------------
  SUBROUTINE GlobalElemAdjacency( Mesh, ElemAdj, ElemAdjProc, ElemStart, DIM )
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, ALLOCATABLE :: ElemAdj(:),ElemStart(:),ElemAdjProc(:)
    INTEGER :: DIM
    !-------------------------------------
    TYPE(Element_t), POINTER :: MFacePtr(:), Element
    INTEGER :: i,j,k,m,n,max_elfaces,el1,el2,gface_id, gpar_id,gpar_lid,ierr,counter,&
         NBulk,NFaces,Sweep,NIFFaces
    INTEGER, ALLOCATABLE :: ElemConn(:,:), ElemConnPart(:,:), NElConn(:), FaceIFIDX(:),status(:),&
         work_int(:)
    INTEGER, POINTER :: ElFaceIdx(:)
    TYPE(NeighbourList_t), POINTER :: MFaceIFList(:)
    LOGICAL, POINTER :: MFaceIF(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="GlobalElemAdjacency"
    TYPE FaceShare_t
       INTEGER :: count
       INTEGER, ALLOCATABLE :: GFaceIDX(:), GParIDX(:),GParLIDX(:)
    END TYPE FaceShare_t
    TYPE(FaceShare_t), ALLOCATABLE :: SendFaces(:),RecvFaces(:)

    NBulk = Mesh % NumberOfBulkElements

    !Find and globally number mesh faces
    IF(DIM == 3) THEN
      CALL FindMeshFaces3D(Mesh)
      CALL FindMeshEdges3D(Mesh)
      CALL SParFaceNumbering(Mesh)
      MFacePtr => Mesh % Faces
      MFaceIF => Mesh % ParallelInfo % FaceInterface
      MFaceIFList => Mesh % ParallelInfo % FaceNeighbourList
      NFaces = Mesh % NumberOfFaces
    ELSEIF(DIM == 2) THEN
      CALL FindMeshEdges2D(Mesh)
      CALL SParEdgeNumbering(Mesh)
      MFacePtr => Mesh % Edges
      MFaceIF => Mesh % ParallelInfo % EdgeInterface
      MFaceIFList => Mesh % ParallelInfo % EdgeNeighbourList
      NFaces = Mesh % NumberOfEdges
    ELSE
      CALL Fatal(FuncName,"Not implemented in 1D")
    END IF

    max_elfaces = 0
    DO i=1,NBulk
      Element => Mesh % Elements(i)
      max_elfaces = MAX(Element % TYPE % NumberOfFaces, max_elfaces)
    END DO

    ALLOCATE(FaceIFIDX(COUNT(MFaceIF)), &
         ElemConn(max_elfaces,NBulk), &
         ElemConnPart(max_elfaces,NBulk), &
         NElConn(NBulk))
    ElemConn = 0
    NElConn = 0

    !Compute local adjacency and gather interface faces
    counter = 0
    DO i=1,NFaces
      IF(MFaceIF(i)) THEN
        counter = counter + 1
        FaceIFIDX(counter) = i
      ELSE
        !Populate the local graph using non-interface faces
        IF(.NOT. ASSOCIATED(MFacePtr(i) % BoundaryInfo % Left) .OR. &
             .NOT. ASSOCIATED(MFacePtr(i) % BoundaryInfo % Right)) CYCLE
        el1 = MFacePtr(i) % BoundaryInfo % Left % ElementIndex
        el2 = MFacePtr(i) % BoundaryInfo % Right % ElementIndex

        NElConn(el1) = NElConn(el1) + 1
        ElemConn(NElConn(el1),el1) = MFacePtr(i) % BoundaryInfo % Right % GElementIndex
        ElemConnPart(NElConn(el1),el1) = ParEnv % MyPE

        NElConn(el2) = NElConn(el2) + 1
        ElemConn(NElConn(el2),el2) = MFacePtr(i) % BoundaryInfo % Left % GElementIndex
        ElemConnPart(NElConn(el2),el2) = Parenv % MyPE
      END IF
    END DO
    NIFFaces = counter

    IF(ANY(NElConn == 0)) CALL Warn(FuncName, 'Disconnected bulk element.')


    !Generate shared Face GElementIndex list & respective parent GElementIndex
    ALLOCATE(SendFaces(ParEnv % PEs),RecvFaces(ParEnv % PEs),work_int(NBulk))

    RecvFaces % Count = 0

    DO Sweep=1,2

      SendFaces % Count = 0
      DO i=1,NIFFaces
        IF (.NOT. ASSOCIATED(MFaceIFList(FaceIFIDX(i)) % Neighbours)) &
             CALL Fatal(FuncName,"Interface face has no Neighborlist % Neighbours")

        gface_id = MFacePtr(FaceIFIDX(i)) % GElementIndex

        IF(ASSOCIATED(MFacePtr(FaceIFIDX(i)) % BoundaryInfo % Left)) THEN
          gpar_id = MFacePtr(FaceIFIDX(i)) % BoundaryInfo % Left % GElementIndex
          gpar_lid = MFacePtr(FaceIFIDX(i)) % BoundaryInfo % Left % ElementIndex
        ELSE IF(ASSOCIATED(MFacePtr(FaceIFIDX(i)) % BoundaryInfo % Right)) THEN
          gpar_id = MFacePtr(FaceIFIDX(i)) % BoundaryInfo % Right % GElementIndex
          gpar_lid = MFacePtr(FaceIFIDX(i)) % BoundaryInfo % Right % ElementIndex
        ELSE
          CALL Fatal(FuncName, "Face has no parent element!")
        END IF

        DO j=1,SIZE(MFaceIFList(FaceIFIDX(i)) % Neighbours)
          k = MFaceIFList(FaceIFIDX(i)) % Neighbours(j) + 1
          IF(k==ParEnv % MyPE+1) CYCLE
          SendFaces(k) % count = SendFaces(k) % count + 1
          n = SendFaces(k) % count
          !Actually write the data
          IF(Sweep == 2) THEN
            SendFaces(k) % GFaceIDX(n) = gface_id
            SendFaces(k) % GParIDX(n) = gpar_id
            SendFaces(k) % GParLIDX(n) = gpar_lid
          END IF
        END DO

      END DO

      IF(Sweep==1) THEN
        DO i=1,ParEnv % PEs
          n = SendFaces(i) % count
          ALLOCATE(SendFaces(i) % GFaceIDX(n),&
               SendFaces(i) % GParIDX(n),&
               SendFaces(i) % GParLIDX(n))
        END DO
      END IF

      IF(Sweep==2) THEN
        DO i=1,ParEnv % PEs
          n = SendFaces(i) % count
          DO j=1,n
            work_int(j) = j
          END DO
          CALL SortI(n, SendFaces(i) % GFaceIDX, work_int)
          SendFaces(i) % GParIDX = SendFaces(i) % GParIDX(work_int(1:n))
          SendFaces(i) % GParLIDX = SendFaces(i) % GParLIDX(work_int(1:n))
        END DO
      END IF
    END DO

    !Send number of shared faces
    ALLOCATE(status(ParEnv % PEs * 2))
    DO i=1,ParEnv % PEs
      CALL MPI_IRECV(RecvFaces(i) % count, 1, MPI_INTEGER, i-1, 194, ELMER_COMM_WORLD, &
           status(i), ierr)
      CALL MPI_SEND(SendFaces(i) % count,  1, MPI_INTEGER, i-1, 194, ELMER_COMM_WORLD, ierr)
    END DO

    CALL MPI_Waitall(ParEnv % PEs, status(1:ParEnv % PES), MPI_STATUSES_IGNORE, ierr)
    status = MPI_REQUEST_NULL

    !transmit shared
    DO i=1,ParEnv % PEs
      n = RecvFaces(i) % count
      k = SendFaces(i) % count
      
      ALLOCATE(RecvFaces(i) % GParIDX(n), RecvFaces(i) % GFaceIDX(n))

      CALL MPI_IRECV(RecvFaces(i) % GFaceIDX, n, MPI_INTEGER, i-1, 195, ELMER_COMM_WORLD, &
           status(i*2-1), ierr)
      CALL MPI_IRECV(RecvFaces(i) % GParIDX, n, MPI_INTEGER, i-1, 196, ELMER_COMM_WORLD, &
           status(i*2), ierr)
      CALL MPI_SEND(SendFaces(i) % GFaceIDX,  k, MPI_INTEGER, i-1, 195, ELMER_COMM_WORLD, ierr)
      CALL MPI_SEND(SendFaces(i) % GParIDX,  k, MPI_INTEGER, i-1, 196, ELMER_COMM_WORLD, ierr)
    END DO

    CALL MPI_Waitall(ParEnv % PEs*2, status(1:ParEnv % PEs*2), MPI_STATUSES_IGNORE, ierr)
    status = MPI_REQUEST_NULL


    DO i=1,ParEnv % PEs
      IF(i-1 == ParEnv % MyPE) CYCLE
      counter = 0
      m = 1
      n = 1
      IF(SendFaces(i) % count==0 .OR. RecvFaces(i) % count==0) CYCLE !no shared faces
      DO WHILE(.TRUE.)
        IF(SendFaces(i) % GFaceIDX(m) == RecvFaces(i) % GFaceIDX(n)) THEN
          !Faces match, update parent elem neighbour list
          el1 = SendFaces(i) % GParLIDX(m)
          el2 = RecvFaces(i) % GParIDX(n)

          NElConn(el1) = NElConn(el1) + 1
          ElemConn(NElConn(el1),el1) = el2
          ElemConnPart(NElConn(el1),el1) = i-1

          counter = counter + 1
          m = m + 1
          n = n + 1
        ELSE IF(SendFaces(i) % GFaceIDX(m) > RecvFaces(i) % GFaceIDX(n)) THEN
          n = n + 1
        ELSE
          m = m + 1
        END IF
        IF(m > SendFaces(i) % count .OR. n > RecvFaces(i) % count) EXIT
      END DO
    END DO
    
    !Put the data into CRS format
    ALLOCATE(ElemAdj(SUM(NElConn)), ElemStart(NBulk+1), ElemAdjProc(SUM(NElConn)))

    ElemStart(1) = 1
    DO i=1,NBulk
      ElemAdj(ElemStart(i):ElemStart(i) + NElConn(i) -1) = &
           ElemConn(1:NElConn(i),i)
      ElemAdjProc(ElemStart(i):ElemStart(i) + NElConn(i) -1) = &
           ElemConnPart(1:NElConn(i),i)
      ElemStart(i+1) = ElemStart(i) + NElConn(i)
    END DO

    CALL MPI_BARRIER(ELMER_COMM_WORLD,ierr)
  END SUBROUTINE GlobalElemAdjacency



  !> As the previous routine exept intended for serial meshes without need for communication.
  !---------------------------------------------------------------------------------------------
  SUBROUTINE LocalElemAdjacency( Mesh, ElemAdj, ElemAdjProc, ElemStart, DIM )
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, ALLOCATABLE :: ElemAdj(:),ElemStart(:),ElemAdjProc(:)
    INTEGER :: DIM
    !-------------------------------------
    TYPE(Element_t), POINTER :: MFacePtr(:), Element
    INTEGER :: i,j,k,m,n,max_elfaces,el1,el2,gface_id, gpar_id,gpar_lid,ierr,counter,&
         NBulk,NFaces,Sweep,NIFFaces
    INTEGER, ALLOCATABLE :: ElemConn(:,:), ElemConnPart(:,:), NElConn(:), status(:),&
         work_int(:)
    INTEGER, POINTER :: ElFaceIdx(:)
    TYPE(NeighbourList_t), POINTER :: MFaceIFList(:)
    LOGICAL, POINTER :: MFaceIF(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="LocalElemAdjacency"

    NBulk = Mesh % NumberOfBulkElements
    IF( NBulk == 0 ) RETURN
    
    !Find and globally number mesh faces
    IF(DIM == 3) THEN
      CALL FindMeshFaces3D(Mesh)
      CALL FindMeshEdges3D(Mesh)
      MFacePtr => Mesh % Faces
      NFaces = Mesh % NumberOfFaces
    ELSEIF(DIM == 2) THEN
      CALL FindMeshEdges2D(Mesh)
      MFacePtr => Mesh % Edges
      NFaces = Mesh % NumberOfEdges
    ELSE
      CALL Fatal(FuncName,"Not implemented in 1D")
    END IF

    max_elfaces = 0
    DO i=1,NBulk
      Element => Mesh % Elements(i)
      IF( dim == 3 ) THEN
        max_elfaces = MAX(Element % TYPE % NumberOfFaces, max_elfaces)
      ELSE        
        max_elfaces = MAX(Element % TYPE % NumberOfEdges, max_elfaces)
      END IF
    END DO
    
    ALLOCATE(ElemConn(max_elfaces,NBulk), &
        ElemConnPart(max_elfaces,NBulk), &
        NElConn(NBulk))
    ElemConn = 0
    NElConn = 0

    !Compute local adjacency and gather interface faces
    counter = 0
    DO i=1,NFaces
      !Populate the local graph using non-interface faces
      IF(.NOT. ASSOCIATED(MFacePtr(i) % BoundaryInfo % Left) .OR. &
          .NOT. ASSOCIATED(MFacePtr(i) % BoundaryInfo % Right)) CYCLE
      
      counter = counter + 1

      el1 = MFacePtr(i) % BoundaryInfo % Left % ElementIndex
      el2 = MFacePtr(i) % BoundaryInfo % Right % ElementIndex
            
      NElConn(el1) = NElConn(el1) + 1
      ElemConn(NElConn(el1),el1) = MFacePtr(i) % BoundaryInfo % Right % ElementIndex
      ElemConnPart(NElConn(el1),el1) = ParEnv % MyPE

      NElConn(el2) = NElConn(el2) + 1
      ElemConn(NElConn(el2),el2) = MFacePtr(i) % BoundaryInfo % Left % ElementIndex
      ElemConnPart(NElConn(el2),el2) = Parenv % MyPE
    END DO
    
    IF(ANY(NElConn == 0)) CALL Warn(FuncName, 'Disconnected bulk element.')
    
    !Put the data into CRS format
    ALLOCATE(ElemAdj(SUM(NElConn)), ElemStart(NBulk+1), ElemAdjProc(SUM(NElConn)))

    ElemStart(1) = 1
    DO i=1,NBulk
      ElemAdj(ElemStart(i):ElemStart(i) + NElConn(i) -1) = &
           ElemConn(1:NElConn(i),i)
      ElemAdjProc(ElemStart(i):ElemStart(i) + NElConn(i) -1) = &
           ElemConnPart(1:NElConn(i),i)
      ElemStart(i+1) = ElemStart(i) + NElConn(i)
    END DO

  END SUBROUTINE LocalElemAdjacency


  !>Identify elements on partition boundaries and return GElementIndexes
  !>of elements with which these elements share a face
  !>This was originally developed for partitioning purposes, but it is *not used*
  !>GlobalElemAdjacency is used instead, which produces face-based connectivity info
  !----------------------------------------------------------------------------------------
  SUBROUTINE MeshParallelDualGraph( Mesh, ElemAdj, ElemStart, ElemIdx, ElemAdjProc, COMM )

    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, ALLOCATABLE :: ElemAdj(:),ElemStart(:),ElemIdx(:),ElemAdjProc(:)
    INTEGER, OPTIONAL :: COMM
    !---------------------------
    TYPE(Element_t), POINTER :: Element, Faces(:)
    INTEGER, POINTER :: neighList(:)
    INTEGER :: MPI_COMM, i,j,k,n,DIM,counter,put_count,NBulk,NNodes,loc,NNeighParts,NElems,&
         ierr,part,tot_send,nn,nl,nli,nti,vli,vti,ElemCount
    INTEGER, ALLOCATABLE :: NeighParts(:),NodeNeighs(:),NNcount(:),NeighPartIDX(:),&
         eptr(:),eind(:),vptr(:),vind(:),stats(:),work_arr(:),proc_arr(:)
    REAL(KIND=dp) :: t1,t2
    LOGICAL :: Debug
    LOGICAL, ALLOCATABLE :: IsNeighbour(:), keep_mask(:)
    LOGICAL, POINTER :: SharedNode(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="MeshParallelDualGraph"

    TYPE PartShareList_t
       INTEGER, ALLOCATABLE :: GNodeNums(:), NodeNums(:), GElemNums(:),&
            NodeElemCount(:),NodeElemStart(:),AllNodesMap(:)
       INTEGER :: NNodes, NElems, Part
    END TYPE PartShareList_t
    TYPE(PartShareList_t), ALLOCATABLE, TARGET :: PartShareList(:),PartRecvList(:)
    TYPE(PartShareList_t), POINTER :: PSL,PSLR

    Debug = .TRUE.
    IF(PRESENT(COMM)) THEN
      MPI_COMM = COMM
    ELSE
      MPI_COMM = ELMER_COMM_WORLD
    END IF

    DIM = CoordinateSystemDimension()
    IF(DIM == 1) CALL Fatal(FuncName,"1D not implemented")
    !Elements which share at least adjshare nodes are
    !considered adjacent

    NBulk = Mesh % NumberOfBulkElements
    NNodes = Mesh % NumberOfNodes

    SharedNode => Mesh % ParallelInfo % Interface
    ALLOCATE(IsNeighbour(ParEnv % PEs),&
         NodeNeighs(NNodes*2),&
         NNcount(ParEnv % PEs)) !<- likely too big

    NNcount = 0
    IsNeighbour = .FALSE.

    !Count number of nodes we share w/ each partition
    DO i=1,NNodes
      IF( .NOT. SharedNode(i)) CYCLE
      DO j=1,SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
        k = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(j)
        IF(k == ParEnv % MyPE) CYCLE
        NNcount(k+1) = NNcount(k+1) + 1
      END DO
    END DO

    IsNeighbour = NNcount > 0
    NNeighParts = COUNT(IsNeighbour)
    ALLOCATE(PartShareList(NNeighParts),PartRecvList(NNeighParts),NeighPartIDX(ParEnv % PEs))
    NeighPartIDX = 0

    !Make a reverse LUT for partition neighbourhood & allocate structures
    counter = 0
    DO i=1,ParEnv % PEs
      IF(.NOT. IsNeighbour(i)) CYCLE
      counter = counter + 1
      PSL => PartShareList(counter)
      PSLR => PartRecvList(counter)

      !NeighParts(counter) = i-1
      NeighPartIDX(i) = counter
      PSL % Part = i-1
      PSLR % Part = i-1
      ALLOCATE(PSL % NodeNums(NNcount(i)),&
           PSLR % GNodeNums(NNcount(i)),&
           PSL % GNodeNums(NNcount(i)),&
           PSL % GElemNums(NNcount(i) *30),&  !overallocate
           PSL % NodeElemStart(NNcount(i)+1),&
           PSLR % NodeElemStart(NNcount(i)+1),&
           PSL % AllNodesMap(NNodes)& !For every node in mesh, the position (if present) in GNodeNums
           )

      PSL % NodeNums = 0
      PSL % GNodeNums = 0
      PSL % GElemNums = 0
      PSL % NodeElemStart = 0

      PSL % NNodes = 0
      PSL % NElems = 0
      PSLR % NNodes = 0
      PSLR % NElems = 0
      PSL % AllNodesMap = 0
      PSL % NodeElemStart = 0
    END DO

    t1 = CPUTime()
    !Fill the PartShareList structure with the nodenumbers shared with partitions
    NNcount = 0 !reset and reuse
    DO i=1,NNodes
      IF( .NOT. SharedNode(i)) CYCLE
      DO j=1,SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
        k = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(j)
        IF(k == ParEnv % MyPE) CYCLE
        NNcount(k+1) = NNcount(k+1) + 1
        n = NeighPartIDX(k+1)
        PSL => PartShareList(n)
        PSL % NNodes = PSL % NNodes + 1 !duplicate from NNcount?
        PSL % NodeNums(NNcount(k+1)) = i
        PSL % GNodeNums(NNcount(k+1)) = Mesh % ParallelInfo % GlobalDOFs(i)
      END DO
      CALL SortI(PSL % NNodes, PSL % GNodeNums, PSL % NodeNums)
    END DO

    !TODO - we already do this in GetElemAdj
    !----This code block copied from ElmerMeshToDualGraph----
    ! Copy mesh to CRS structure
    ALLOCATE(eptr(NBulk+1), eind(NBulk*Mesh % MaxElementNodes))
    eptr(1)=1 ! Fortran numbering
    DO i=1, NBulk
      nli = eptr(i) ! Fortran numbering
      Element => Mesh % Elements(i)
      nl = Element % TYPE % NumberOfNodes
      nti = nli+nl-1
      eind(nli:nti) = Element % NodeIndexes(1:nl) ! Fortran numbering
      eptr(i+1) = nli+nl
    END DO
    ! Construct vertex to element list
    CALL VertexToElementList(nbulk, NNodes, eptr, eind, vptr, vind)
    !--------------------------------------------------------

    t2 = CPUTime()
    PRINT *, ParEnv % MyPE,' time to build element vertex list',t2-t1
    t1 = CPUTime()
    !For each partition, for each shared node on that partition, gather element numbers
    DO i=1,NNeighParts
      PSL => PartShareList(i)

      counter = 1
      DO j=1,PSL % NNodes
        n = PSL % NodeNums(j)

        vli = vptr(n)
        vti = vptr(n+1)-1
        NElems = vti-vli+1

        !Resize if necessary
        IF(counter + NElems - 1 > SIZE(PSL % GElemNums))&
             CALL ResizeIntArray(PSL % GElemNums,PSL % NElems*2)

        CALL Sort(NElems, vind(vli:vti))
        PSL % NodeElemStart(j) = counter
        PSL % GElemNums(counter:counter + NElems-1) = Mesh % Elements(vind(vli:vti)) % GElementIndex
        counter = counter + NElems
        PSL % NElems = PSL % NElems + NElems
        PSL % NodeElemStart(j+1) = counter !This is wrong
      END DO
    END DO

    t2 = CPUTime()
    PRINT *, ParEnv % MyPE,' time to gather elements',t2-t1
    t1 = CPUTime()

    tot_send = 0
    ALLOCATE(stats(NNeighParts*4))
    stats = MPI_REQUEST_NULL
    DO i=1,NNeighParts
      PSL => PartShareList(i)
      PSLR => PartRecvList(i)
      part = PSL % part

      !NB: this will currently only work with ELMER_COMM_WORLD 
      !due to the partition numbering
      CALL MPI_ISEND(PSL % NNodes, 1, MPI_INTEGER, part,&
           198, MPI_COMM,stats(i*4-3),ierr)
      IF(ierr /= 0) PRINT *,'ruh roh1'
      CALL MPI_IRECV(PSLR % NNodes, 1, MPI_INTEGER, part,&
           198,MPI_COMM,stats(i*4-2),ierr)
      IF(ierr /= 0) PRINT *,'ruh roh2'

      CALL MPI_ISEND(PSL % NElems, 1, MPI_INTEGER, part,&
           199,MPI_COMM,stats(i*4-1),ierr)
      IF(ierr /= 0) PRINT *,'ruh roh3'
      CALL MPI_IRECV(PSLR % NElems, 1, MPI_INTEGER, part,&
           199,MPI_COMM,stats(i*4),ierr)
      IF(ierr /= 0) PRINT *,'ruh roh4'

      tot_send = tot_send + PSL % NNodes + PSL % NElems + 1
    END DO

    tot_send = tot_send * 2 !(send and receive?)
    PRINT *,ParEnv % MyPE,' tot send: ',tot_send
    CALL CheckBuffer(tot_send)

    CALL MPI_Waitall(NNeighParts*4, stats, MPI_STATUSES_IGNORE, ierr)
    stats = MPI_REQUEST_NULL

    DO i=1,NNeighParts
      PSL => PartShareList(i)
      PSLR => PartRecvList(i)
      part = PSL % part

      IF(PSL % NNodes /= PSLR % NNodes) &
           CALL Fatal(FuncName, "Interface node count mismatch")
      ALLOCATE(PSLR % GElemNums(PSLR % NElems+1))

      CALL MPI_ISEND(PSL % GElemNums, PSL % NElems, MPI_INTEGER, part,&
           200, MPI_COMM,stats(i*4-3),ierr)
      IF(ierr /= 0) PRINT *,'ruh roh5'
      CALL MPI_IRECV(PSLR % GElemNums, PSLR % NElems, MPI_INTEGER, part,&
           200,MPI_COMM,stats(i*4-2),ierr)
      IF(ierr /= 0) PRINT *,'ruh roh6'

      CALL MPI_ISEND(PSL % NodeElemStart, PSL % NNodes+1, MPI_INTEGER, part,&
           201, MPI_COMM,stats(i*4-1),ierr)
      IF(ierr /= 0) PRINT *,'ruh roh7'

      CALL MPI_IRECV(PSLR % NodeElemStart, PSLR % NNodes+1, MPI_INTEGER, part,&
           201,MPI_COMM,stats(i*4),ierr)
      IF(ierr /= 0) PRINT *,'ruh roh8'
    END DO

    CALL MPI_Waitall(NNeighParts*4, stats, MPI_STATUSES_IGNORE, ierr)
    stats = MPI_REQUEST_NULL

    t2 = CPUTime()
    PRINT *, ParEnv % MyPE,' time to do MPI',t2-t1
    t1 = CPUTime()

    PRINT *,ParEnv % MyPE,' test tot elems: ',SUM(PartShareList % NElems)

    !-------Quick test---------
    CALL MPI_Waitall(NNeighParts*4, stats, MPI_STATUSES_IGNORE, ierr)
    stats = MPI_REQUEST_NULL

    DO i=1,NNeighParts
      PSL => PartShareList(i)
      PSLR => PartRecvList(i)
      part = PSL % part

      CALL MPI_ISEND(PSL % GNodeNums, PSL % NNodes, MPI_INTEGER, part,&
           202, MPI_COMM,stats(i*4-3),ierr)
      IF(ierr /= 0) PRINT *,'ruh roh5'
      CALL MPI_IRECV(PSLR % GNodeNums, PSLR % NNodes, MPI_INTEGER, part,&
           202,MPI_COMM,stats(i*4-2),ierr)
      IF(ierr /= 0) PRINT *,'ruh roh6'
    END DO

    CALL MPI_Waitall(NNeighParts*4, stats, MPI_STATUSES_IGNORE, ierr)
    stats = MPI_REQUEST_NULL

    DO i=1,NNeighParts
      PSL => PartShareList(i)
      PSLR => PartRecvList(i)
      part = PSL % part
      DO j=1,PSL % NNodes
        IF(PSL % GNodeNums(j) /= PSLR % GNodeNums(j)) THEN
          PRINT *,ParEnv % MyPE,' mismatched node nums: ',PSL % GNodeNums(j), PSLR % GNodeNums(j)
          CALL Fatal(FuncName, "node num mismatch")
        END IF
      END DO
    END DO
    !-----------------------------------

    !For each shared part, construct the map of local nodenumbers to PSL % GNodeNums
    DO i=1,NNeighParts
      PSL => PartShareList(i)
      DO j=1,NNodes
        IF(.NOT. ANY(Mesh % ParallelInfo  % NeighbourList(j) % Neighbours == PSL % Part)) CYCLE
        PSL % AllNodesMap(j) = SearchI(PSL % NNodes, PSL % NodeNums,j)
      END DO
    END DO

    t2 = CPUTime()
    PRINT *, ParEnv % MyPE,' time to make partition nodemap',t2-t1
    t1 = CPUTime()

    ALLOCATE(work_arr(100),proc_arr(100),keep_mask(100))

    PRINT *, ParEnv % MyPE,' size to allocate: ',SUM(PartShareList % NElems)*20
    ALLOCATE(ElemIdx(NBulk), &
         ElemStart(NBulk+1),&
         ElemAdjProc(SUM(PartShareList % NElems)*20),&
         ElemAdj(SUM(PartShareList % NElems)*20)) !overallocate

    !For each shared element determine elements from other parts
    ! and actually build the adjacency list
    ElemIdx = 0
    ElemAdjProc = 0
    ElemStart(1) = 1
    ElemCount = 1
    DO i=1,NBulk
      work_arr = 0
      proc_arr = 0
      keep_mask = .TRUE.

      counter = 0
      Element => Mesh % Elements(i)

      DO j=1,Element % TYPE % NumberOfNodes

        nn = Element % NodeIndexes(j)
        IF(.NOT. SharedNode(nn)) CYCLE
        neighList => Mesh % ParallelInfo % NeighbourList(nn) % Neighbours

        DO k=1,SIZE(neighList)
          part = neighList(k)
          IF(part == ParEnv % MyPE) CYCLE

          PSL => PartShareList(NeighPartIDX(part+1))
          PSLR => PartRecvList(NeighPartIDX(part+1))

          loc = PSL % AllNodesMap(nn)
          IF(loc == 0) CALL Fatal(FuncName, "Programming error: node not found in part")
          IF(PSL % NodeNums(loc) /= nn) CALL Fatal(FuncName, "Programming error with AllNodesMap")
          nli = PSLR % NodeElemStart(loc) !first elem loc shared by this node
          nl = PSLR % NodeElemStart(loc+1)-1 !last elem loc
          work_arr(counter+1:counter+1+nl-nli) = PSLR % GElemNums(nli:nl)
          proc_arr(counter+1:counter+1+nl-nli) = PSLR % part !<- array indexing bloody mess
          counter = counter + nl-nli+1
          IF(counter > SIZE(work_arr)) THEN
            CALL ResizeIntArray(work_arr,counter*2)
            CALL ResizeIntArray(proc_arr,counter*2)
            DEALLOCATE(keep_mask)
            ALLOCATE(keep_mask(counter*2)) !don't need to preserve values
            keep_mask = .TRUE.
          END IF
        END DO

      END DO

      IF(counter > 0) THEN

        !Clear out duplicates
        CALL SortI(counter,work_arr, proc_arr)
        DO j=2,counter
          IF(work_arr(j) == work_arr(j-1)) keep_mask(j) = .FALSE.
        END DO
        put_count = COUNT(keep_mask(1:counter))
        IF(put_count < 1) CALL Fatal(FuncName, "Programming Error: masked all elements!")


        !Check for resize
        IF(ElemStart(ElemCount) + put_count-1 > SIZE(ElemAdj)) &
             CALL ResizeIntArray(ElemAdj, (ElemStart(ElemCount) + put_count-1)*2)

        ElemAdj(ElemStart(ElemCount):ElemStart(ElemCount) + put_count-1) = &
             PACK(work_arr(1:counter),keep_mask(1:counter))
        ElemAdjProc(ElemStart(ElemCount):ElemStart(ElemCount) + put_count-1) = &
             PACK(proc_arr(1:counter), keep_mask(1:counter))

        ElemIdx(ElemCount) = i

        ElemCount = ElemCount + 1
        ElemStart(ElemCount) = ElemStart(ElemCount-1) + put_count
      END IF

    END DO
    ElemCount = ElemCount - 1

    !Pack to actual sizes:

    !The adjacency list
    DEALLOCATE(work_arr)

    ALLOCATE(work_arr(ElemStart(ElemCount+1)-1))
    work_arr = ElemAdj(1:ElemStart(ElemCount+1)-1)
    DEALLOCATE(ElemAdj)
    CALL MOVE_ALLOC(work_arr, ElemAdj) 

    !The neighbour proc list
    ALLOCATE(work_arr(ElemStart(ElemCount+1)-1))
    work_arr = ElemAdjProc(1:ElemStart(ElemCount+1)-1)
    DEALLOCATE(ElemAdjProc)
    CALL MOVE_ALLOC(work_arr, ElemAdjProc)

    !The element indices for which we found parallel neighbours
    ALLOCATE(work_arr(ElemCount))
    work_arr = ElemIdx(1:ElemCount)
    DEALLOCATE(ElemIdx)
    CALL MOVE_ALLOC(work_arr, ElemIdx)

    !The CRS starts
    ALLOCATE(work_arr(ElemCount+1))
    work_arr = ElemStart(1:ElemCount+1)
    DEALLOCATE(ElemStart)
    CALL MOVE_ALLOC(work_arr, ElemStart)

    t2 = CPUTime()
    PRINT *, ParEnv % MyPE,' time to make partition nodemap',t2-t1
    t1 = CPUTime()

    !Need to send either 1) the 3 node numbers for each elem or 2) all the nodenums and array list
    !Create a hash function for each face we share with other partitions
    !Three large prime numbers - note collision is not an issue as we will
    !check global node numbers directly

  CONTAINS

    SUBROUTINE VertexToElementList(nelem, nvertex, eptr, eind, vptr, vind)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nelem, nvertex
      INTEGER :: eptr(:), eind(:)
      INTEGER, ALLOCATABLE :: vptr(:), vind(:)
      INTEGER :: i, j, v, eli, eti, ind, tmpi, tmpip, allocstat

      ! Initialize vertex structure (enough storage for nvertex vertices
      ! having eptr(nelem+1) elements)
      ALLOCATE(vptr(nvertex+1), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('VertexToElementList', &
           'Vertex allocation failed!')

      vptr = 0

      ! For each element

      ! Compute number of elements attached to each vertex (size of lists)
      DO i=1,nelem
        eli = eptr(i)
        eti = eptr(i+1)-1

        DO j=eli, eti
          vptr(eind(j))=vptr(eind(j))+1
        END DO
      END DO

      ! Compute in-place cumulative sum (row pointers!)
      CALL ComputeCRSIndexes(nvertex, vptr)

      ! Allocate vertex to element lists
      ALLOCATE(vind(vptr(nvertex+1)), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('VertexToElementList', &
           'Vertex allocation failed!')

      ! Construct element lists for each vertex
      DO i=1,nelem
        eli = eptr(i)
        eti = eptr(i+1)-1

        ! For each vertex in element
        DO j=eli, eti
          ! Add connection to vertex eind(j)
          ind = eind(j)
          vind(vptr(ind))=i
          vptr(ind)=vptr(ind)+1
        END DO
      END DO

      ! Correct row pointers
      DO i=nvertex,2,-1
        vptr(i)=vptr(i-1)
      END DO
      vptr(1)=1
    END SUBROUTINE VertexToElementList

    SUBROUTINE ResizeIntArray(Arr, new_size)
      INTEGER, ALLOCATABLE :: Arr(:)
      INTEGER :: new_size
      !------------------------
      INTEGER, ALLOCATABLE :: workArr(:)

      ALLOCATE(workArr(new_size))
      workArr(1:SIZE(Arr)) = Arr
      DEALLOCATE(Arr)
      CALL MOVE_ALLOC(workArr, Arr)
    END SUBROUTINE ResizeIntArray

  END SUBROUTINE MeshParallelDualGraph

  !Turns a masked node list into a real stream for sending. Largely
  !superceded by RedistributeMesh, which handles the entire mesh together.
  !May still have some use.
  !--------------------------------------------------------------------------
  SUBROUTINE PackNodesToSend(Mesh, Mask, GDOFs, NodeCoords, DIM)

    IMPLICIT NONE

    TYPE(Mesh_t),POINTER :: Mesh
    LOGICAL :: Mask(:)
    INTEGER, ALLOCATABLE :: GDOFs(:)
    INTEGER, OPTIONAL :: DIM
    REAL(KIND=dp), ALLOCATABLE :: NodeCoords(:)
    !-----------------------
    INTEGER :: i,SendCount, counter, the_dim

    SendCount = COUNT(Mask)

    IF(PRESENT(DIM)) THEN
      the_dim = DIM
    ELSE
      the_dim = CoordinateSystemDimension()
    END IF

    !------------- append partition no and stream_size
    ALLOCATE(GDOFs(2+SendCount),NodeCoords(SendCount*the_dim))

    GDOFs(1) = ParEnv % MyPE
    GDOFs(2) = SendCount
    counter = 0
    DO i=1,Mesh % NumberOfNodes
      IF(.NOT. Mask(i)) CYCLE
      counter = counter + 1
      GDOFs(counter+2) = Mesh % ParallelInfo % GlobalDOFs(i)
      NodeCoords((counter-1)*the_dim +1) = Mesh % Nodes % x(i)
      NodeCoords((counter-1)*the_dim +2) = Mesh % Nodes % y(i)
      IF(the_dim == 3) NodeCoords((counter-1)*the_dim  +3) = Mesh % Nodes % z(i)
    END DO

  END SUBROUTINE PackNodesToSend

  !Inverse of PackNodesToSend - superceded by RedistributeMesh
  SUBROUTINE UnpackNodesSent(GDOFs, NodeCoords, Nodes, DIM, node_parts)
    INTEGER, ALLOCATABLE :: GDOFs(:)
    REAL(KIND=dp) :: NodeCoords(:)
    TYPE(Nodes_t) :: Nodes
    INTEGER, OPTIONAL :: DIM
    INTEGER, ALLOCATABLE, OPTIONAL :: node_parts(:)
    !-----------------------
    INTEGER :: i, the_dim, node_count,stream_pos,work_pos,partid,part_nnodes
    INTEGER, ALLOCATABLE :: work_int(:)
    LOGICAL :: have_partids

    IF(PRESENT(DIM)) THEN
      the_dim = DIM
    ELSE
      the_dim = CoordinateSystemDimension()
    END IF 
    have_partids = PRESENT(node_parts)

    node_count = SIZE(NodeCoords)/DIM

    ALLOCATE(Nodes % x(node_count),&
         Nodes % y(node_count),&
         Nodes % z(node_count))

    DO i=1,node_count
      Nodes % x(i) = NodeCoords(i*the_dim-2)
      Nodes % y(i) = NodeCoords(i*the_dim-1)
      IF(the_dim == 3) Nodes % z(i) = NodeCoords(i*3)
    END DO
    IF(the_dim /= 3) Nodes % z = 0.0

    !Strip the header info (partition and count) out from GDOFs:
    ALLOCATE(work_int(node_count))
    IF(have_partids) ALLOCATE(node_parts(node_count))

    stream_pos = 1
    work_pos = 1
    work_int = 0

    DO WHILE(.TRUE.)
      partid = GDOFs(stream_pos)
      stream_pos = stream_pos + 1
      part_nnodes = GDOFs(stream_pos)

      stream_pos = stream_pos + 1
      work_int(work_pos:work_pos + part_nnodes - 1) = GDOFs(stream_pos:stream_pos + part_nnodes - 1)
      IF(have_partids) node_parts(work_pos:work_pos + part_nnodes - 1) = partid

      work_pos = work_pos + part_nnodes

      stream_pos = stream_pos + part_nnodes
      IF(stream_pos > SIZE(GDOFs)) EXIT
    END DO

    DEALLOCATE(GDOFs)
    CALL MOVE_ALLOC(work_int, GDOFs)
  END SUBROUTINE UnpackNodesSent

  !Converts element datastructure into a single integer stream to facilitate
  !sending to another partition. In addition to element numbers, type, nodes and BC/Body ID,
  !an integer custom_tag may be provided to indicate, for example, what the receiving process
  !should do with the elements (remesh, remove, keep fixed)
  !
  !This is largely superceded by RedistributeMesh
  SUBROUTINE PackElemsToSend(Mesh, Mask, ElemStream, custom_tag)

    IMPLICIT NONE

    TYPE(Mesh_t),POINTER :: Mesh
    LOGICAL :: Mask(:)
    INTEGER, ALLOCATABLE :: ElemStream(:)
    INTEGER, OPTIONAL :: custom_tag(:)
    !----------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER, ALLOCATABLE :: work_int(:)
    INTEGER :: SendCount, el_counter, stream_counter
    INTEGER :: i,nblk,nbdry,nodecount,ENodeCount
    LOGICAL :: have_tags

    have_tags = PRESENT(custom_tag)

    SendCount = COUNT(Mask)
    nblk = Mesh % NumberOfBulkElements
    nbdry = Mesh % NumberOfBoundaryElements

    !space for partitionID, nodecount, streamsize, then per elem: 
    !nodenums + ID + tag +  (OPTIONAL custom_tag) + TYPE
    IF(have_tags) THEN
      ALLOCATE(ElemStream(3 + SendCount * (Mesh % MaxElementNodes + 4))) 
    ELSE
      ALLOCATE(ElemStream(3 + SendCount * (Mesh % MaxElementNodes + 3))) 
    END IF

    el_counter = 0
    nodecount = 0

    ElemStream(1) = ParEnv % MyPE
    ElemStream(2) = SendCount

    !NB: we go back and fill space 3 with stream size later
    stream_counter = 3

    DO i=1,nblk+nbdry
      IF(.NOT. Mask(i)) CYCLE
      el_counter = el_counter + 1
      Element => Mesh % Elements(i)

      ENodeCount = Element % TYPE % NumberOfNodes

      !Pack the global element index
      stream_counter = stream_counter + 1
      ElemStream(stream_counter) = Element % GElementIndex

      !Pack the type
      stream_counter = stream_counter + 1
      ElemStream(stream_counter) = Element % TYPE % ElementCode

      !Pack either the body id or boundary id
      stream_counter = stream_counter + 1
      IF(i <= nblk) THEN
        ElemStream(stream_counter) = Element % BodyID
      ELSE
        ElemStream(stream_counter) = Element % BoundaryInfo % Constraint
      END IF

      !Pack the custom_tag if present
      IF(have_tags) THEN
        stream_counter = stream_counter + 1
        ElemStream(stream_counter) = custom_tag(i)
      END IF

      !Pack the global node numbers
      stream_counter = stream_counter + 1
      ElemStream(stream_counter:stream_counter + ENodeCount - 1) = &
           Mesh % ParallelInfo % GlobalDOFs(Element % NodeIndexes(1:ENodeCount))

      stream_counter = stream_counter + ENodeCount - 1
    END DO

    !Go back and fill in the stream size
    ElemStream(3) = stream_counter - 3

    IF(stream_counter > SIZE(ElemStream)) CALL Fatal("PackElemsToSend","Too much data to pack - &
         &this should be impossible - MaxElementNodes probably wrong!")

    ALLOCATE(work_int(stream_counter))
    work_int(1:stream_counter) = ElemStream(1:stream_counter)
    DEALLOCATE(ElemStream) !<- this shouldn't be necessary - Cray bug!
    CALL MOVE_ALLOC(work_int, ElemStream)

  END SUBROUTINE PackElemsToSend

  !Performs inverse operation of PackElemsToSend
  !
  !This is largely superceded by RedistributeMesh
  SUBROUTINE UnpackElemsSent(ElemStream, Elements, DIM, elem_parts, custom_tag)
    INTEGER :: ElemStream(:)
    TYPE(Element_t), ALLOCATABLE :: Elements(:)
    INTEGER, ALLOCATABLE, OPTIONAL :: elem_parts(:), custom_tag(:)
    INTEGER, OPTIONAL :: DIM
    !------------------------------------------
    INTEGER :: i,j,stream_counter,elem_count,elem_totcount, part_count, the_dim,&
         type_code,NElNodes,elem_dim,partid,stream_pos, stream_size
    LOGICAL :: have_tags, have_partids, bdry_elem

    IF(PRESENT(DIM)) THEN
      the_dim = DIM
    ELSE
      the_dim = CoordinateSystemDimension()
    END IF

    have_tags = PRESENT(custom_tag)
    have_partids = PRESENT(elem_parts)

    !Scan the data to identify number of partitions, total elem count:
    stream_pos = 1
    part_count = 0
    elem_totcount = 0
    DO WHILE(.TRUE.)
      partid = ElemStream(stream_pos)
      part_count = part_count + 1

      stream_pos = stream_pos + 1
      elem_totcount = elem_totcount + ElemStream(stream_pos)
      
      stream_pos = stream_pos + 1
      stream_pos = stream_pos + ElemStream(stream_pos) + 1

      IF(stream_pos > SIZE(ElemStream)) EXIT
    END DO


    ALLOCATE(Elements(elem_totcount))
    IF(have_partids) ALLOCATE(elem_parts(elem_totcount))
    IF(have_tags) ALLOCATE(custom_tag(elem_totcount))
    stream_counter = 0
    elem_totcount = 0

    DO j=1,part_count

      stream_counter = stream_counter + 1
      partid = ElemStream(stream_counter)
      stream_counter = stream_counter + 1
      elem_count = ElemStream(stream_counter)
      stream_counter = stream_counter + 1
      stream_size = ElemStream(stream_counter)

      DO i=elem_totcount+1,  elem_totcount + elem_count

        !Set the element global index
        stream_counter = stream_counter + 1
        Elements(i) % GElementIndex = ElemStream(stream_counter)

        !Set the element type
        stream_counter = stream_counter + 1
        type_code = ElemStream(stream_counter)
        Elements(i) % TYPE => GetElementType(type_code,.FALSE.)
        NElNodes = Elements(i) % TYPE % NumberOfNodes


        !Set either the boundary info or body ID
        elem_dim = Elements(i) % TYPE % DIMENSION
        bdry_elem = elem_dim < the_dim

        stream_counter = stream_counter + 1
        IF(bdry_elem) THEN
          ALLOCATE(Elements(i) % BoundaryInfo)
          Elements(i) % BoundaryInfo % constraint = ElemStream(stream_counter)
        ELSE
          Elements(i) % BodyID = ElemStream(stream_counter)
        END IF

        !set custom tag if sent
        IF(have_tags) THEN
          stream_counter = stream_counter + 1
          custom_tag(i) = ElemStream(stream_counter)
        END IF

        !Set node indexes
        stream_counter = stream_counter + 1
        ALLOCATE(Elements(i) % NodeIndexes(NElNodes))
        Elements(i) % NodeIndexes = ElemStream(stream_counter: stream_counter + NElNodes - 1)

        !Set element partition ID
        IF(have_partids) elem_parts(i) = partid

        stream_counter = stream_counter + NElNodes - 1
      END DO
      elem_totcount = elem_totcount + elem_count
    END DO

  END SUBROUTINE UnpackElemsSent


  !Cleanly removes elements & nodes from a mesh based on mask
  !Any element containing a removed node (RmNode) will be deleted
  !May optionally specify which elements to remove
  !Not currently used
  SUBROUTINE CutMesh(Mesh, RmNode, RmElem)
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: RmNode(:)
    LOGICAL, OPTIONAL, TARGET :: RmElem(:)
    !--------------------------------
    TYPE(Element_t), POINTER :: Element, Work_Elements(:)
    TYPE(Nodes_t), POINTER :: Nodes
    TYPE(NeighbourList_t), POINTER :: work_neighlist(:)
    REAL(KIND=dp), ALLOCATABLE :: work_xyz(:,:)
    REAL(KIND=dp), POINTER CONTIG :: work_x(:),work_y(:), work_z(:)
    INTEGER :: i,j,counter,NNodes,NBulk, NBdry,NewNNodes, NewNElems, NewNBulk,&
         NewNbdry, ElNNodes
    INTEGER, ALLOCATABLE :: Nodeno_map(:),work_int(:),EIdx_map(:)
    INTEGER, POINTER :: NodeIndexes(:), work_pInt(:)
    LOGICAL, POINTER :: RmElement(:),work_logical(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="CutMesh"
    NNodes = Mesh % NumberOfNodes
    NBulk = Mesh % NumberOfBulkElements
    NBdry = Mesh % NumberOfBoundaryElements

    IF(PRESENT(RmElem)) THEN
      RmElement => RmElem
    ELSE
      !Mark elements containing deleted nodes for deletion
      ALLOCATE(RmElement(NBulk + Nbdry))
      RmElement = .FALSE.
      DO i=1,NBulk + NBdry
        Element => Mesh % Elements(i)
        NodeIndexes => Element % NodeIndexes
        ElNNodes = Element % TYPE % NumberOfNodes
        IF(ANY(RmNode(NodeIndexes(1:ElNNodes)))) RmElement(i) = .TRUE.
        !Get rid of orphans
        IF(i > NBulk) THEN
          IF(ASSOCIATED(Element % BoundaryInfo)) THEN
            !If the main body element is removed, remove this BC elem
            IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
              j = Element % BoundaryInfo % Left % ElementIndex
              IF(RmElement(j)) RmElement(i) = .TRUE.
            END IF
            !Point % Right => NULL() if % Right element is removed
            IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
              j = Element % BoundaryInfo % Right % ElementIndex
              IF(RmElement(j)) Element % BoundaryInfo % Right => NULL()
            END IF
          END IF
        END IF
      END DO
    END IF

    !Removing nodes implies shifting element nodeindexes
    !Map pre -> post deletion node nums
    ALLOCATE(Nodeno_map(NNodes))
    Nodeno_map = 0
    counter = 0
    DO i=1,NNodes
      IF(RmNode(i)) CYCLE
      counter = counter + 1
      Nodeno_map(i) = counter
    END DO

    !Update the element nodeindexes
    DO i=1,NBulk+NBdry
      Element => Mesh % Elements(i)
      IF(RmElement(i)) CYCLE
      DO j=1,Element % TYPE % NumberOfNodes
        Element % NodeIndexes(j) = Nodeno_map(Element % NodeIndexes(j))
        IF(Element % NodeIndexes(j) == 0) CALL Fatal(FuncName, &
             "Programming error: mapped nodeno = 0")
      END DO
    END DO

    !Clear out deleted nodes
    Nodes => Mesh % Nodes
    NewNNodes = COUNT(.NOT. RmNode)

    ALLOCATE(work_x(NewNNodes),&
         work_y(NewNNodes),&
         work_z(NewNNodes))

    counter = 0
    DO i=1,NNodes
      IF(RmNode(i)) CYCLE
      counter = counter + 1
      work_x(counter) = Nodes % x(i)
      work_y(counter) = Nodes % y(i)
      work_z(counter) = Nodes % z(i)
    END DO

    DEALLOCATE(Nodes % x, Nodes % y, Nodes % z)
    Nodes % x => work_x
    Nodes % y => work_y
    Nodes % z => work_z

    Nodes % NumberOfNodes = NewNNodes
    Mesh % NumberOfNodes = NewNNodes

    !Clear out ParallelInfo
    IF(ASSOCIATED(Mesh % ParallelInfo % GlobalDOFs)) THEN
      ALLOCATE(work_pInt(NewNNodes))
      counter = 0
      DO i=1,NNodes
        IF(RmNode(i)) CYCLE
        counter = counter + 1
        work_pInt(counter) = Mesh % ParallelInfo % GlobalDOFs(i)
      END DO
      DEALLOCATE(Mesh % ParallelInfo % GlobalDOFs)
      Mesh % ParallelInfo % GlobalDOFs => work_pInt
      work_pInt => NULL()
    END IF

    !Get rid of NeighbourList
    IF(ASSOCIATED(Mesh % ParallelInfo % NeighbourList)) THEN
      ALLOCATE(work_neighlist(NewNNodes))
      DO i=1,NNodes
        IF(.NOT. RmNode(i)) CYCLE
        IF(ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) &
             DEALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
      END DO

      counter = 0
      DO i=1,NNodes
        IF(RmNode(i)) CYCLE
        counter = counter + 1
        work_neighlist(counter) % Neighbours => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
      END DO
      DEALLOCATE(Mesh % ParallelInfo % NeighbourList)
      Mesh % ParallelInfo % NeighbourList => work_neighlist
      work_neighlist => NULL()
    END IF

    !Get rid of ParallelInfo % INTERFACE
    IF(ASSOCIATED(Mesh % ParallelInfo % INTERFACE)) THEN
      ALLOCATE(work_logical(NewNNodes))
      counter = 0
      DO i=1,NNodes
        IF(RmNode(i)) CYCLE
        counter = counter + 1
        work_logical(counter) = Mesh % ParallelInfo % INTERFACE(i)
      END DO
      DEALLOCATE(Mesh % ParallelInfo % INTERFACE)
      Mesh % ParallelInfo % INTERFACE => work_logical
      work_logical => NULL()
    END IF

    !TODO - Mesh % Edges - see ReleaseMeshEdgeTables
    IF ( ASSOCIATED( Mesh % Edges ) ) THEN
      CALL Fatal(FuncName,"Clearing out Mesh % Edges not yet implemented.")
    END IF

    !TODO - Mesh % Faces - see ReleaseMeshFaceTables
    IF ( ASSOCIATED( Mesh % Faces ) ) THEN
      CALL Fatal(FuncName,"Clearing out Mesh % Faces not yet implemented.")
    END IF

    !TODO - Mesh % ViewFactors -  see ReleaseMeshFactorTables
    IF (ASSOCIATED(Mesh % ViewFactors) ) THEN
      CALL Fatal(FuncName,"Clearing out Mesh % ViewFactors not yet implemented.")
    END IF

    !TODO - Mesh % Projector - see FreeMatrix in ReleaseMesh
    IF (ASSOCIATED(Mesh % Projector) ) THEN
      CALL Fatal(FuncName,"Clearing out Mesh % Projector not yet implemented.")
    END IF

    !TODO - Mesh % RootQuadrant - see FreeQuadrantTree
    IF( ASSOCIATED( Mesh % RootQuadrant ) ) THEN
      CALL Fatal(FuncName,"Clearing out Mesh % RootQuadrant not yet implemented.")
    END IF

    NewNElems = COUNT(.NOT. RmElement)


    !Clear out deleted elements structures - taken from ReleaseMesh
    DO i=1,NBulk+NBdry
      IF(.NOT. RmElement(i)) CYCLE

      !          Boundaryinfo structure for boundary elements
      !          ---------------------------------------------
      IF ( Mesh % Elements(i) % Copy ) CYCLE

      IF ( i > NBulk ) THEN
        IF ( ASSOCIATED( Mesh % Elements(i) % BoundaryInfo ) ) THEN
          IF (ASSOCIATED(Mesh % Elements(i) % BoundaryInfo % GebhardtFactors)) THEN
            IF ( ASSOCIATED( Mesh % Elements(i) % BoundaryInfo % &
                 GebhardtFactors % Elements ) ) THEN
              DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % &
                   GebhardtFactors % Elements )
              DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % &
                   GebhardtFactors % Factors )
            END IF
            DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % GebhardtFactors )
          END IF
          DEALLOCATE( Mesh % Elements(i) % BoundaryInfo )
        END IF
      END IF

      IF ( ASSOCIATED( Mesh % Elements(i) % NodeIndexes ) ) &
           DEALLOCATE( Mesh % Elements(i) % NodeIndexes )
      Mesh % Elements(i) % NodeIndexes => NULL()

      IF ( ASSOCIATED( Mesh % Elements(i) % EdgeIndexes ) ) &
           DEALLOCATE( Mesh % Elements(i) % EdgeIndexes )
      Mesh % Elements(i) % EdgeIndexes => NULL()

      IF ( ASSOCIATED( Mesh % Elements(i) % FaceIndexes ) ) &
           DEALLOCATE( Mesh % Elements(i) % FaceIndexes )
      Mesh % Elements(i) % FaceIndexes => NULL()

      IF ( ASSOCIATED( Mesh % Elements(i) % DGIndexes ) ) &
           DEALLOCATE( Mesh % Elements(i) % DGIndexes )
      Mesh % Elements(i) % DGIndexes => NULL()

      IF ( ASSOCIATED( Mesh % Elements(i) % BubbleIndexes ) ) &
           DEALLOCATE( Mesh % Elements(i) % BubbleIndexes )
      Mesh % Elements(i) % BubbleIndexes => NULL()

      ! This creates problems later on!!!
      !IF ( ASSOCIATED( Mesh % Elements(i) % PDefs ) ) &
      !   DEALLOCATE( Mesh % Elements(i) % PDefs )

      Mesh % Elements(i) % PDefs => NULL()
    END DO

    !Construct a map of old element indexes => new element indexes
    ALLOCATE(EIdx_map(NBdry+NBulk))
    EIdx_map = 0
    counter = 0
    DO i=1,NBulk+NBdry
      IF(RmElement(i)) CYCLE
      counter = counter + 1
      IF(Mesh % Elements(i) % ElementIndex /= i) CALL Warn(FuncName,&
           "Assumption Elements(i) % ElementIndex == i not valid! Expect memory corruption...")
      EIdx_map(Mesh % Elements(i) % ElementIndex) = counter
    END DO

    !Repoint elements
    ALLOCATE(work_elements(NewNElems))
    counter = 0
    DO i=1,Nbulk+NBdry
      IF(RmElement(i)) CYCLE
      counter = counter + 1

      Element => Mesh % Elements(i)
      work_elements(counter) = Element

      !Repoint BoundaryInfo % Left, % Right
      IF(i > NBdry) THEN
        IF(ASSOCIATED(Element % BoundaryInfo)) THEN
          IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
            j = Element % BoundaryInfo % Left % ElementIndex
            work_elements(counter) % BoundaryInfo % Left => work_elements(EIdx_map(j))
          END IF
          IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
            j = Element % BoundaryInfo % Right % ElementIndex
            work_elements(counter) % BoundaryInfo % Right => work_elements(EIdx_map(j))
          END IF
        END IF
      END IF
    END DO

    !Update ElementIndexes
    DO i=1,NewNElems
      work_elements(i) % ElementIndex = i
    END DO

    DEALLOCATE(Mesh % Elements)
    Mesh % Elements => work_elements
    work_elements => NULL()

    Mesh % NumberOfBulkElements = COUNT(.NOT. RmElement(1:nbulk))
    Mesh % NumberOfBoundaryElements = COUNT(.NOT. RmElement(nbulk+1:nbulk+nbdry))

    IF(.NOT. PRESENT(RmElem)) DEALLOCATE(RmElement)

  END SUBROUTINE CutMesh



  !> Takes an existing mesh and a repartitioning vector and redisributes the mesh
  !> among the partitions. Assumes that global node and element indexes are sane. 
  !------------------------------------------------------------------------------
  FUNCTION RedistributeMesh( Model, Mesh, ParallelMesh, FreeOldMesh ) RESULT( NewMesh )
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh, NewMesh
    LOGICAL :: ParallelMesh, FreeOldMesh

    TYPE( MeshPack_t), ALLOCATABLE, TARGET :: SentPack(:), RecPack(:)
    INTEGER, POINTER :: NewPart(:)
    INTEGER :: NoPartitions, newnodes, newnbdry, newnbulk, dim, minind, maxind, n, ierr,i
    INTEGER, ALLOCATABLE :: GlobalToLocal(:)
    CHARACTER(*), PARAMETER :: FuncName = 'RedistributeMesh'


    CALL Info(FuncName,'Distributing mesh structures in parallel',Level=6)
    CALL ResetTimer(FuncName)

    NoPartitions = ParEnv % PEs

    IF( NoPartitions == 1 ) THEN
      CALL Warn(FuncName,'Redistribution does not make sense for one partition!')
    END IF

    IF( .NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal(FuncName,'Mesh not associated')
    END IF

    IF( Mesh % NumberOfNodes > 0 ) THEN
      IF(.NOT. ASSOCIATED( Mesh % RePartition ) ) THEN
        CALL Fatal(FuncName,'Repartitioning information not associated')
      END IF
      NewPart => Mesh % RePartition

      n = MAXVAL( NewPart )
      IF( n > NoPartitions ) THEN
        CALL Fatal(FuncName,'Partition number exceeds process number: '//TRIM(I2S(n)))
      END IF
    END IF

    ! For convenience lets always communicate all coordinates
    ! Also, dimension might not have been set at this point for all domains
    dim = 3

    ! 0) Given the new partitioning add potential nodes to the interface
    CALL UpdateInterfaceNodeCandidates( Mesh )

    ! 1) First pack the mesh for parallel communication
    CALL PackMeshPieces(Model, Mesh, NewPart, ParallelMesh, NoPartitions, SentPack, dim)

    ! 2) Then sent the pieces among different partitions
    CALL CommunicateMeshPieces(Model, Mesh, ParallelMesh, NoPartitions, SentPack, RecPack)

    ! 3) Calculate element element and node counts, and creates new global2local numbering
    CALL LocalNumberingMeshPieces(Model, Mesh, NewPart, ParallelMesh, NoPartitions, RecPack, &
         GlobalToLocal, newnodes, newnbulk, newnbdry, minind, maxind)

    NewMesh => AllocateMesh( newnbulk, newnbdry, newnodes, InitParallel = .TRUE.)    
    
    ! 4) Finally unpack and glue the pieces on an existing mesh
    CALL UnpackMeshPieces(Model, Mesh, NewMesh, NewPart, minind, &
         maxind, RecPack, ParallelMesh, GlobalToLocal, dim)

    CALL FindRepartitionInterfaces(Model, NewMesh, dim)

    NewMesh % Name = Mesh % Name
    NewMesh % MeshDim = Mesh % MeshDim
    NewMesh % OutputActive = .TRUE.

    DO i=1,NewMesh % NumberOfBulkElements + NewMesh % NumberOfBoundaryElements
      IF(NewMesh % Elements(i) % ElementIndex /= i) CALL Fatal(FuncName, "Bad element index")
      IF(ANY(NewMesh % Elements(i) % NodeIndexes <= 0) .OR. &
           ANY(NewMesh % Elements(i) % NodeIndexes > newnodes)) THEN
        CALL Fatal(FuncName,' bad elem nodeindexes')
      END IF
    END DO

    IF( FreeOldMesh ) CALL ReleaseMesh( Mesh )

    CALL Info('RedistributeMesh','Deallocating temporal packed structures',Level=20)
    DO i=1,NoPartitions
      IF( ALLOCATED( SentPack(i) % idata ) ) DEALLOCATE( SentPack(i) % idata )
      IF( ALLOCATED( SentPack(i) % rdata ) ) DEALLOCATE( SentPack(i) % rdata )
      IF( ALLOCATED( RecPack(i) % idata ) ) DEALLOCATE( RecPack(i) % idata )
      IF( ALLOCATED( RecPack(i) % rdata ) ) DEALLOCATE( RecPack(i) % rdata )
    END DO

    DEALLOCATE( GlobalToLocal )

    CALL Info('CommunicateMeshPieces','Waiting for MPI barrier',Level=15)
    CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )

    CALL CheckTimer(FuncName,Level=5,Delete=.TRUE.)
    CALL Info(FuncName,'Distributing mesh finished',Level=8)

  END FUNCTION RedistributeMesh


  !> Converts element datastructure into a integer and real stream to facilitate
  !> sending to another partition.
  !------------------------------------------------------------------------------
  SUBROUTINE PackMeshPieces(Model, Mesh, NewPart, ParallelMesh, NoPartitions, SentPack, dim)

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: ParallelMesh
    INTEGER, POINTER :: NewPart(:)
    INTEGER :: NoPartitions, dim
    TYPE( MeshPack_t), ALLOCATABLE, TARGET :: SentPack(:)
    !------------------------
    TYPE(Element_t), POINTER :: Element, Parent
    INTEGER :: i,j,k,n,nblk,nbdry,allocstat,part,elemcode,geom_id,sweep
    LOGICAL :: CheckNeighbours, IsBulk
    TYPE(NeighbourList_t),POINTER  :: NeighbourList(:)
    TYPE(MeshPack_t), POINTER :: PPack
    INTEGER, POINTER :: TmpPart(:)

    CALL Info('PackMeshPieces','Packing mesh pieces for sending',Level=8)

    ! Allocate and initialize the structures used to communicate thes mesh
    n = NoPartitions
    ALLOCATE( SentPack( n ) )

    SentPack(1:n) % NumberOfNodes = 0
    SentPack(1:n) % NumberOfBulkElements = 0
    SentPack(1:n) % NumberOfBoundaryElements = 0
    SentPack(1:n) % icount = 0
    SentPack(1:n) % rcount = 0
    SentPack(1:n) % lcount = 0
    SentPack(1:n) % indpos = 0
    SentPack(1:n) % bcpos = 0

    IF( Mesh % NumberOfNodes == 0 ) THEN
      CALL Info('PackMeshPieces','Mesh is empty, nothing to pack',Level=10)
      RETURN
    END IF

    nblk = Mesh % NumberOfBulkElements
    nbdry = Mesh % NumberOfBoundaryElements

    IF( SIZE( NewPart ) < nblk + nbdry ) THEN
      CALL Info('PackMeshPieces','Growing the partition vector to accounts BCs',Level=8)
      ALLOCATE( TmpPart( nblk + nbdry ) )
      TmpPart(1:nblk) = NewPart(1:nblk)
      TmpPart(nblk+1:nblk+nbdry) = 0
      DEALLOCATE( NewPart )
      NewPart => TmpPart
      Mesh % RePartition => TmpPart
      NULLIFY( TmpPart )
    END IF

    IF( ANY(NewPart(nblk+1:nblk+nbdry) <= 0 ) ) THEN
      CALL Info('PackMeshPieces','Inheriting partition vector for BCs',Level=8)
      DO i=nblk+1,nblk+nbdry
        Element => Mesh % Elements(i)
        IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) THEN
          CALL Fatal('PackMeshPieces','Boundary element lacking boundary info')
        END IF
        Parent => Element % BoundaryInfo % Left
        IF( .NOT. ASSOCIATED(Parent) ) THEN
          Parent => Element % BoundaryInfo % Right
        END IF
        IF( .NOT. ASSOCIATED( Parent ) ) THEN
          CALL Fatal('PackMeshPieces','Boundary element lacking parent info')
        END IF
        NewPart(i) = NewPart(Parent % ElementIndex)
      END DO
    END IF

    CheckNeighbours = .FALSE.
    IF( ParallelMesh ) THEN
      IF( ASSOCIATED( Mesh % ParallelInfo % NeighbourList ) ) THEN
        CheckNeighbours = .TRUE.
        NeighbourList => Mesh % ParallelInfo % NeighbourList
      END IF
      IF( CheckNeighbours ) THEN
        CALL info('PackMeshPieces','Using existing information of shared nodes',Level=10)
      END IF
    END IF

    ! In the 1st sweep we calculate amount of data to be sent.
    ! In the 2nd sweep we populate the idata and rdata vectors.
    !--------------------------------------------------------------
    DO Sweep = 1, 2

      DO part=1,NoPartitions

        IF( part <= 0 ) CYCLE ! This one for possible elements to be eliminated
        IF( part-1 == ParEnv % MyPe ) CYCLE  ! Skip this partition

        PPack => SentPack(part)
        PPack % icount = 5 ! offset for the above header info
        PPack % rcount = 0
        PPack % lcount = 0
      END DO

      !  Go through the elements
      DO i=1,nblk+nbdry
        ! Mark the offset for boundary data. It is needed later when performing the unpacking.
        IF( Sweep == 1 .AND. i == nblk + 1 ) THEN
          SentPack(1:NoPartitions) % bcpos = SentPack(1:NoPartitions) % icount
        END IF

        part = NewPart(i)
        IF( part <= 0 .OR. part-1 == ParEnv % MyPe ) CYCLE

        IsBulk = ( i <= nblk )

        Element => Mesh % Elements(i)
        elemcode = Element % Type % ElementCode
        n = Element % TYPE % NumberOfNodes

        PPack => SentPack(part)

        IF( Sweep == 1) THEN
          IF( IsBulk ) THEN
            PPack % NumberOfBulkElements = PPack % NumberOfBulkElements + 1
          ELSE
            PPack % NumberOfBoundaryElements = PPack % NumberOfBoundaryElements + 1
          END IF
          IF( IsBulk ) THEN
            PPack % icount = PPack % icount + 3
          ELSE
            PPack % icount = PPack % icount + 5
          END IF

        ELSE IF( Sweep == 2 ) THEN
          ! Populate the table on the second sweep

          ! Pack the global element index
          IF( ParallelMesh ) THEN
            PPack % idata(PPack % icount+1) = Element % GElementIndex
          ELSE
            PPack % idata(PPack % icount+1) = Element % ElementIndex
          END IF

          ! Pack the type
          PPack % idata(PPack % icount+2) = elemcode

          ! Pack either the body id or boundary id
          IF( IsBulk ) THEN
            PPack % idata(PPack % icount+3) = Element % BodyID
            PPack % icount = PPack % icount + 3
          ELSE
            PPack % idata(PPack % icount+3) = Element % BoundaryInfo % Constraint
            PPack % idata(PPack % icount+4:5) = 0
            IF( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
              IF( ParallelMesh ) THEN
                PPack % idata(PPack % icount+4) = Element % BoundaryInfo % Left % GElementIndex
              ELSE
                PPack % idata(PPack % icount+4) = Element % BoundaryInfo % Left % ElementIndex
              END IF
            END IF
            IF( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
              IF( ParallelMesh ) THEN
                PPack % idata(PPack % icount+5) = Element % BoundaryInfo % Right % GElementIndex
              ELSE
                PPack % idata(PPack % icount+5) = Element % BoundaryInfo % Right % ElementIndex
              END IF
            END IF
            PPack % icount = PPack % icount + 5
          END IF

          ! Pack node indexes
          IF( ParallelMesh ) THEN
            PPack % idata(PPack % icount+1:PPack % icount+n) = &
                 Mesh % ParallelInfo % GlobalDOFs(Element % NodeIndexes(1:n))
          ELSE
            PPack % idata(PPack % icount+1:PPack % icount+n) = &
                 Element % NodeIndexes(1:n)
          END IF
        END IF

        ! Advance the counter for the data
        PPack % icount = PPack % icount + n
      END DO


      IF( Sweep == 1 ) THEN
        ! Set the offset for the nodal data. It is needed in unpackig.
        SentPack(1:NoPartitions) % indpos = SentPack(1:NoPartitions) % icount

        ! We need to allocate a logical mask to mark the nodes to sent to given partition
        ! Note that we only want to allocate the flag for partitions that also recieve
        ! some elements.
        DO part=1,NoPartitions
          PPack => SentPack(part)
          IF( PPack % icount > 5 ) THEN
            ALLOCATE( PPack % NodeMask( Mesh % NumberOfNodes ) )
            PPack % NodeMask = .FALSE.
          END IF
        END DO

        ! For each active partition mark the nodes that must be sent
        DO i=1,nblk+nbdry
          part = NewPart(i)

          IF( part <= 0 .OR. part-1 == ParEnv % MyPe ) CYCLE

          PPack => SentPack(part)

          Element => Mesh % Elements(i)
          elemcode = Element % TYPE % ElementCode
          n = Element % TYPE % NumberOfNodes

          IF( CheckNeighbours ) THEN
            DO j=1,n
              k = Element % NodeIndexes(j)
              ! We may already have the node in the partition
              IF( ANY( part-1 == NeighbourList(k) % Neighbours ) ) CYCLE
              PPack % NodeMask(k) = .TRUE.
            END DO
          ELSE
            PPack % NodeMask(Element % NodeIndexes(1:n)) = .TRUE.
          END IF
        END DO

        ! Add the nodes count to be sent to the data structure
        DO part=1,NoPartitions
          IF( part-1 == ParEnv % MyPe ) CYCLE
          PPack => SentPack(part)

          IF(  PPack % icount <= 5 ) CYCLE

          PPack % NumberOfNodes = COUNT( PPack % NodeMask )
          PPack % icount = PPack % icount + PPack % NumberOfNodes
          PPack % rcount = dim * PPack % NumberOfNodes
          PPack % lcount = PPack % NumberOfNodes
          
          ALLOCATE( PPack % idata(PPack % icount), &
              PPack % rdata(PPack % rcount), &
              PPack % ldata(PPack % lcount ), & 
              STAT = allocstat )
          IF( allocstat /= 0 ) THEN
            CALL Fatal('PackMeshPieces','Could not allocate vectors for data')
          END IF
          PPack % idata = 0
          PPack % rdata = 0.0_dp
          PPack % ldata = .FALSE.

          PPack % idata(1) = PPack % NumberOfNodes
          PPack % idata(2) = PPack % NumberOfBulkElements
          PPack % idata(3) = PPack % NumberOfBoundaryElements
          PPack % idata(4) = PPack % bcpos
          PPack % idata(5) = PPack % indpos
        END DO

      ELSE IF( Sweep == 2 ) THEN

        ! For simplicity this includes a loop over all partitions
        DO part=1,NoPartitions
          IF( part-1 == ParEnv % MyPe ) CYCLE
          PPack => SentPack(part)

          ! No need to sent empty element lists
          IF( PPack % NumberOfNodes == 0 ) CYCLE

          DO i = 1, Mesh % NumberOfNodes
            IF( PPack % NodeMask(i) ) THEN
              PPack % icount = PPack % icount + 1
              IF( ParallelMesh ) THEN
                PPack % idata(Ppack % icount) = Mesh % ParallelInfo % GlobalDOFs(i)
              ELSE
                PPack % idata(Ppack % icount) = i
              END IF

              ! Also add the coordinates for sending
              PPack % rdata(PPack % rcount+1) = Mesh % Nodes % x(i)
              PPack % rdata(PPack % rcount+2) = Mesh % Nodes % y(i)
              IF( dim == 3 ) PPack % rdata(PPack % rcount+3) = Mesh % Nodes % z(i)                           
              PPack % rcount = PPack % rcount + dim
              
              PPack % lcount = PPack % lcount + 1
              PPack % ldata(PPack % lcount) = Mesh % ParallelInfo % Interface(i)
            END IF
          END DO

        END DO
      END IF
    END DO

    CALL Info('PackMeshPieces','Finished packing mesh pieces',Level=8)

  END SUBROUTINE PackMeshPieces


  ! Communicate the packed mesh between partitions using MPI.
  !-----------------------------------------------------------------------
  SUBROUTINE CommunicateMeshPieces(Model, Mesh, ParallelMesh, NoPartitions, SentPack, RecPack)

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: NoPartitions
    LOGICAL :: ParallelMesh
    TYPE( MeshPack_t), ALLOCATABLE, TARGET :: SentPack(:), RecPack(:)
    !------------------------

    INTEGER :: i,j,n,ierr,ni,nr,nl,status(MPI_STATUS_SIZE)
    INTEGER, ALLOCATABLE :: Requests(:)

    CALL Info('CommunicateMeshPieces','communicating mesh pieces in parallel',Level=6)

    ni = SUM( SentPack(1:NoPartitions) % icount )
    CALL Info('CommunicateMeshPieces','Number of integer values to sent: '//TRIM(I2S(ni)),Level=8)
    nr = SUM( SentPack(1:NoPartitions) % rcount )
    CALL Info('CommunicateMeshPieces','Number of real values to sent: '//TRIM(I2S(nr)),Level=8)
    nl = SUM( SentPack(1:NoPartitions) % lcount )
    CALL Info('CommunicateMeshPieces','Number of logical values to sent: '//TRIM(I2S(nl)),Level=8)

    n = NoPartitions
    ALLOCATE( RecPack( n ) )

    RecPack(1:n) % NumberOfNodes = 0
    RecPack(1:n) % NumberOfBulkElements = 0
    RecPack(1:n) % NumberOfBoundaryElements = 0
    RecPack(1:n) % icount = 0
    RecPack(1:n) % rcount = 0
    RecPack(1:n) % lcount = 0
    RecPack(1:n) % indpos = 0
    RecPack(1:n) % bcpos = 0


    CALL CheckBuffer( nl + ni*4 + nr*8 + (NoPartitions-1)* ( 4*2 + &
         2*MPI_BSEND_OVERHEAD ) )

    ! Sent data sizes:
    !--------------------------
    ALLOCATE( Requests(NoPartitions) )
    DO i=1,NoPartitions
      IF( i-1 == ParEnv % Mype ) CYCLE
      CALL MPI_BSEND( SentPack(i) % icount, 1, MPI_INTEGER, i-1, &
           1000, ELMER_COMM_WORLD, ierr )
      CALL MPI_BSEND( SentPack(i) % rcount, 1, MPI_INTEGER, i-1, &
           1001, ELMER_COMM_WORLD, ierr )
      CALL MPI_BSEND( SentPack(i) % lcount, 1, MPI_INTEGER, i-1, &
           1002, ELMER_COMM_WORLD, ierr )
    END DO

    ! Recieve data sizes:
    !--------------------------
    DO i = 1, NoPartitions
      IF( i-1 == ParEnv % MyPe ) CYCLE
      CALL MPI_RECV( RecPack(i) % icount, 1, MPI_INTEGER, i-1, &
           1000, ELMER_COMM_WORLD, status, ierr )
      CALL MPI_RECV( RecPack(i) % rcount, 1, MPI_INTEGER, i-1, &
           1001, ELMER_COMM_WORLD, status, ierr )
      CALL MPI_RECV( RecPack(i) % lcount, 1, MPI_INTEGER, i-1, &
           1002, ELMER_COMM_WORLD, status, ierr )
    END DO

    CALL Info('PackMeshPieces','Waiting for the 1st barrier',Level=15)
    CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )

    n = SUM( RecPack(1:NoPartitions) % icount )
    CALL Info('PackDataToSend','Number of integer values to recieve: '//TRIM(I2S(n)),Level=8)
    n = SUM( RecPack(1:NoPartitions) % rcount )
    CALL Info('PackDataToSend','Number of real values to recieve: '//TRIM(I2S(n)),Level=8)
    n = SUM( RecPack(1:NoPartitions) % lcount )
    CALL Info('PackDataToSend','Number of logical values to recieve: '//TRIM(I2S(n)),Level=8)

    ! Allocate data sizes for recieving data
    !----------------------------------------
    DO i=1,NoPartitions
      IF( i-1 == ParEnv % Mype ) CYCLE
      IF( RecPack(i) % icount > 5 ) THEN
        ALLOCATE( RecPack(i) % idata( RecPack(i) % icount) )
        ALLOCATE( RecPack(i) % rdata( RecPack(i) % rcount) )
        ALLOCATE( RecPack(i) % ldata( RecPack(i) % lcount) )
      END IF
    END DO

    ! Sent data:
    !--------------------------
    CALL Info('CommunicateMeshPieces','Now sending the actual data',Level=12)
    DO i=1,NoPartitions
      IF( i-1 == ParEnv % Mype ) CYCLE
      IF( SentPack(i) % icount > 5 ) THEN
        CALL MPI_BSEND( SentPack(i) % idata, SentPack(i) % icount, MPI_INTEGER, i-1, &
             1003, ELMER_COMM_WORLD, ierr )
        CALL MPI_BSEND( SentPack(i) % rdata, SentPack(i) % rcount, MPI_DOUBLE_PRECISION, i-1, &
             1004, ELMER_COMM_WORLD, ierr )
        CALL MPI_BSEND( SentPack(i) % ldata, SentPack(i) % lcount, MPI_LOGICAL, i-1, &
             1005, ELMER_COMM_WORLD, ierr )
      END IF
    END DO

    ! Recieve data:
    !--------------------------
    CALL Info('CommunicateMeshPieces','Now recieving the actual integer data',Level=12)
    DO i = 1, NoPartitions
      IF( i-1 == ParEnv % MyPe ) CYCLE
      IF( RecPack(i) % icount > 5 ) THEN
        CALL MPI_RECV( RecPack(i) % idata, RecPack(i) % icount, MPI_INTEGER, i-1, &
             1003, ELMER_COMM_WORLD, status, ierr )
        CALL MPI_RECV( RecPack(i) % rdata, RecPack(i) % rcount, MPI_DOUBLE_PRECISION, i-1, &
             1004, ELMER_COMM_WORLD, status, ierr )
        CALL MPI_RECV( RecPack(i) % ldata, RecPack(i) % lcount, MPI_LOGICAL, i-1, &
             1005, ELMER_COMM_WORLD, status, ierr )
      END IF
    END DO

    CALL Info('CommunicateMeshPieces','Waiting for the 2nd barrier',Level=15)
    CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )

    CALL Info('CommunicateMeshPieces','Finished communicating mesh pieces',Level=8)

  END SUBROUTINE CommunicateMeshPieces




  !> Calculates new local count of elements and nodes, and creates the new
  !> local numbering for the nodes. The GlobalToLocal numbering is needed
  !> when unpacking the data to the new mesh.
  !------------------------------------------------------------------------------
  SUBROUTINE LocalNumberingMeshPieces(Model, Mesh, NewPart, ParallelMesh, NoPartitions, &
       RecPack, GlobalToLocal, newnodes, newnbulk, newnbdry, minind, maxind)

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, ALLOCATABLE :: GlobalToLocal(:)
    INTEGER, POINTER :: NewPart(:)
    INTEGER :: NoPartitions, newnodes, newnbulk, newnbdry, minind, maxind
    TYPE( MeshPack_t), ALLOCATABLE, TARGET :: RecPack(:)
    LOGICAL :: ParallelMesh
    !---------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementType_t), POINTER :: Etype
    INTEGER :: i,j,k,n,t,allocstat,part,elemcode
    INTEGER :: gind,lind,rcount,icount,nbrdy,i1,i2
    TYPE(NeighbourList_t),POINTER  :: NeighbourList(:)
    TYPE(MeshPack_t), POINTER :: PPack

    CALL Info('LocalNumberingMeshPieces','Renumbering local nodes in each partition',Level=8)

    newnbulk = 0
    newnbdry = 0

    ! Compute the number of elements staying
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      IF( NewPart(i)-1 == ParEnv % MyPe ) THEN
        IF( i <= Mesh % NumberOfBulkElements ) THEN
          newnbulk = newnbulk + 1
        ELSE
          newnbdry = newnbdry + 1
        END IF
      END IF
    END DO

    ! Add the number of elements coming from different partitions
    DO part=1,NoPartitions
      IF( i-1 == ParEnv % Mype ) CYCLE
      PPack => RecPack(part)
      IF( PPack % icount <= 5 ) CYCLE

      PPack % NumberOfNodes = PPack % idata(1)
      PPack % NumberOfBulkElements = PPack % idata(2)
      PPack % NumberOfBoundaryElements = PPack % idata(3)
      PPack % bcpos = PPack % idata(4)
      PPack % indpos = PPack % idata(5)

      newnbulk = newnbulk + PPack % NumberOfBulkElements
      newnbdry = newnbdry + PPack % NumberOfBoundaryElements
    END DO

    CALL Info('LocalNumberingMeshPieces','Combined number of elements: '&
         //TRIM(I2S(newnbulk+newnbdry)),Level=8)


    ! Find the range of initial global indeces
    ! This is conservative since it includes all the initial global indexes
    n = Mesh % NumberOfNodes
    IF( n > 0 ) THEN
      IF( ParallelMesh ) THEN
        maxind = MAXVAL( Mesh % ParallelInfo % GlobalDofs(1:n) )
        minind = MINVAL( Mesh % ParallelInfo % GlobalDofs(1:n) )
      ELSE
        maxind = Mesh % NumberOfNodes
        minind = 1
      END IF
    ELSE
      minind = HUGE( minind )
      maxind = 0
    END IF

    ! also check the imported global indexes
    DO part=1,NoPartitions
      IF( i-1 == ParEnv % Mype ) CYCLE
      PPack => RecPack(part)
      IF( PPack % icount <= 5 ) CYCLE
      i1 = PPack % indpos + 1
      i2 = PPack % icount
      minind = MIN( minind, MINVAL( PPack % idata(i1:i2) ) )
      maxind = MAX( maxind, MAXVAL( PPack % idata(i1:i2) ) )
    END DO

    CALL Info('LocalNumberingMeshPieces','Global index range '&
         //TRIM(I2S(minind))//' to '//TRIM(I2S(maxind)),Level=12)

    ! Allocate the vector for local renumbering
    ALLOCATE( GlobalToLocal(minind:maxind), STAT=allocstat)
    IF( allocstat /= 0 ) THEN
      CALL Fatal('LocalNumberingMeshPieces','Could not allocate vectors for indexes')
    END IF
    GlobalToLocal = 0

    ! Check which of the staying nodes are used
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      IF( NewPart(i)-1 == ParEnv % MyPe ) THEN
        Element => Mesh % Elements(i)
        n = Element % Type % NumberOfNodes

        DO j=1,n
          k = Element % NodeIndexes(j)
          IF( ParallelMesh ) k = Mesh % ParallelInfo % GlobalDOFs(k)
          IF( k < minind .OR. k > maxind ) THEN
            CALL Fatal('LocalNumberingMeshPieces','k out of bounds')
          END IF
          GlobalToLocal(k) = 1
        END DO
      END IF
    END DO

    ! Add the imported nodes and their global index
    ! Cycle element % nodes rather than nodes directly
    ! in case the other partition didn't send a node which
    ! we already have (but haven't marked above)
    DO part=1,NoPartitions
      IF( part-1 == ParEnv % MyPe ) CYCLE
      PPack => RecPack(part)
      IF( PPack % icount <= 5 ) CYCLE

      icount = 5
      DO i=1, PPack % NumberOfBulkElements

        elemcode = PPack % idata(icount+2)

        Etype => GetElementType( elemcode )
        n = Etype % NumberOfNodes
        icount = icount + 3

        IF( icount + n > SIZE( PPack % idata) ) THEN
          CALL Fatal('LocalNumberingMeshPieces','icount out of range')
        END IF

        DO j=1,n
          gind = PPack % idata(icount+1)

          IF( gind < minind .OR. gind > maxind ) THEN
            CALL Fatal('LocalNumberingMeshPieces','gind out of bounds')
          END IF

          GlobalToLocal(gind) = 1
          icount = icount + 1
        END DO
      END DO
    END DO

    ! Renumber the staying nodes
    newnodes = 0
    DO i=minind, maxind
      IF( GlobalToLocal(i) > 0 ) THEN
        newnodes = newnodes + 1
        GlobalToLocal(i) = newnodes
      END IF
    END DO

    CALL Info('LocalNumberingMeshPieces','Combined number of nodes: '//TRIM(I2S(newnodes)),Level=8)

  END SUBROUTINE LocalNumberingMeshPieces



  !> Converts element data structure from integer and real streams to FE meshes.
  !> The idea is that the elements and nodes are appended on top of an existing
  !> mesh structure such that there could be elements already at the bottom.
  !------------------------------------------------------------------------------
  SUBROUTINE UnpackMeshPieces(Model, Mesh, NewMesh, NewPart, &
       minind, maxind, RecPack, ParallelMesh, GlobalToLocal, dim)

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh, NewMesh
    LOGICAL :: ParallelMesh
    INTEGER, POINTER :: NewPart(:)
    INTEGER, ALLOCATABLE :: GlobalToLocal(:)
    INTEGER :: dim,minind,maxind
    TYPE( MeshPack_t), ALLOCATABLE, TARGET :: RecPack(:)
    !----------------------------------
    TYPE(Element_t), POINTER :: Element, Element0
    INTEGER :: i,j,k,n,t,nbulk,nbdry,allocstat,part,elemcode,elemindex,geom_id,sweep
    INTEGER :: gind,lind,rcount,icount,lcount,minelem,maxelem,newnbdry,newnodes,newnbulk
    LOGICAL :: IsBulk
    TYPE(MeshPack_t), POINTER :: PPack
    INTEGER, ALLOCATABLE :: GlobalToLocalElem(:), LeftParent(:), RightParent(:)
    CHARACTER(*), PARAMETER :: Caller = 'UnpackMeshPieces'
    
    CALL Info(Caller,'Unpacking mesh pieces to form a new mesh',Level=12)


    minelem = HUGE( minelem )
    maxelem = 0

    newnodes = NewMesh % NumberOfNodes
    newnbulk = NewMesh % NumberOfBulkElements 
    newnbdry = NewMesh % NumberOfBoundaryElements 
    
    nbulk = 0
    nbdry = 0

    ! There are temporal arrays needed to inherit the parent information
    ALLOCATE( LeftParent(NewNBulk+1:NewNBulk+NewNBdry ) ) 
    LeftParent = 0
    ALLOCATE( RightParent(NewNBulk+1:NewNBulk+NewNBdry ) ) 
    RightParent = 0
       
    
    CALL Info(Caller,'Copying staying elements',Level=20)
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      IF( NewPart(i)-1 /= ParEnv % MyPe ) CYCLE

      IsBulk = ( i <= Mesh % NumberOfBulkElements )

      IF( IsBulk ) THEN
        nbulk = nbulk + 1
        t = nbulk
      ELSE
        nbdry = nbdry + 1
        t = newnbulk + nbdry
      END IF

      Element => NewMesh % Elements(t)
      Element0 => Mesh % Elements(i)

      IF(.NOT. ASSOCIATED( Element ) ) THEN
        CALL Fatal(Caller,'Element not allocated')
      END IF

      Element % Type => Element0 % Type
      n = Element % Type % NumberOfNodes

      IF( n <= 0 .OR. n > 27 ) THEN
        CALL Fatal(Caller,'Invalid number of nodes')
      END IF

      Element % BodyId = Element0 % BodyId

      IF( IsBulk ) THEN
        Element % BoundaryInfo => NULL()
      ELSE
        ALLOCATE( Element % BoundaryInfo )
        Element % BoundaryInfo % Constraint = Element0 % BoundaryInfo % Constraint

        IF( ASSOCIATED( Element0 % BoundaryInfo % Left ) ) THEN
          IF( ParallelMesh ) THEN
            LeftParent(t) = Element0 % BoundaryInfo % Left % GElementIndex
          ELSE
            LeftParent(t) = Element0 % BoundaryInfo % Left % ElementIndex
          END IF
        END IF
        IF( ASSOCIATED( Element0 % BoundaryInfo % Right ) ) THEN
          IF( ParallelMesh ) THEN         
            RightParent(t) = Element0 % BoundaryInfo % Right % GElementIndex
          ELSE
            RightParent(t) = Element0 % BoundaryInfo % Right % ElementIndex            
          END IF
        END IF
      END IF

      IF( ParallelMesh ) THEN
        Element % GElementIndex = Element0 % GElementIndex
      ELSE
        Element % GElementIndex = Element0 % ElementIndex
      END IF

      IF( IsBulk ) THEN
        minelem = MIN( minelem, Element % GElementIndex )
        maxelem = MAX( maxelem, Element % GElementIndex )
      END IF
        
      Element % ElementIndex = t

      ! Change the owner partition of the element
      Element % PartIndex = ParEnv % MyPe

      NULLIFY( Element % NodeIndexes )
      ALLOCATE( Element % NodeIndexes(n), STAT = allocstat )
      IF( allocstat /= 0 ) THEN
        CALL Fatal(Caller,'Cannot allocate '//TRIM(I2S(n))//' node indexes?')
      END IF

      DO j=1,n
        k = Element0 % NodeIndexes(j)
        IF( ParallelMesh ) k = Mesh % ParallelInfo % GlobalDOFs(k)

        ! Renumber the nodes such that the local indexes are always contiguous
        IF( k < minind .OR. k > maxind ) THEN
          CALL Fatal(Caller,'k out of bounds')
        END IF
        Element % NodeIndexes(j) = GlobalToLocal(k)
      END DO
    END DO

    CALL Info(Caller,'Copying staying nodes',Level=20)

    IF( .NOT. ASSOCIATED( NewMesh % ParallelInfo % INTERFACE ) ) THEN
      ALLOCATE( NewMesh % ParallelInfo % INTERFACE( NewMesh % NumberOfNodes ), STAT = allocstat )
      IF( allocstat /= 0 ) THEN
        CALL Fatal(Caller,'Cannot allocate partition interface?')
      END IF
      NewMesh % ParallelInfo % Interface = .FALSE.
    END IF

    
    DO i=1,Mesh % NumberOfNodes
      j = i
      IF( ParallelMesh ) j = Mesh % ParallelInfo % GlobalDofs(i)
      IF( j < minind .OR. j > maxind ) THEN
        CALL Fatal(Caller,'j out of bounds')
      END IF
      k = GlobalToLocal(j)
      IF( k == 0 ) CYCLE
      IF( k > newnodes ) THEN
        CALL Fatal(Caller,'k2 out of bounds: '//TRIM(I2S(k))//' vs. '//TRIM(I2S(newnodes)))
      END IF
      NewMesh % Nodes % x(k) = Mesh % Nodes % x(i)
      NewMesh % Nodes % y(k) = Mesh % Nodes % y(i)
      IF( dim == 3 ) NewMesh % Nodes % z(k) = Mesh % Nodes % z(i)

      NewMesh % ParallelInfo % Interface(k) = Mesh % ParallelInfo % Interface(i)
    END DO


    CALL Info(Caller,'Unpacking incoming elements',Level=20)
    DO part=1,ParEnv % PEs
      IF( part-1 == ParEnv % MyPe ) CYCLE

      PPack => RecPack(part)
      IF( PPack % icount <= 5 ) CYCLE

      CALL Info(Caller,'Unpacking piece '//TRIM(I2S(part))//' with '&
           //TRIM(I2S(PPack % NumberOfBulkElements + PPack % NumberOfBoundaryElements))//&
           ' elements',Level=20)

      ! The five values upfront are used for size info
      icount = 5
      DO i = 1, PPack % NumberOfBulkElements + PPack % NumberOfBoundaryElements

        IsBulk = ( i <= PPack % NumberOfBulkElements )

        IF( IsBulk ) THEN
          nbulk = nbulk + 1
          t = nbulk
        ELSE
          nbdry = nbdry + 1
          t = newnbulk + nbdry
        END IF

        elemindex = PPack % idata(icount+1)
        elemcode = PPack % idata(icount+2)
        geom_id = PPack % idata(icount+3)

        Element => NewMesh % Elements(t)
        IF( .NOT. ASSOCIATED( Element ) ) THEN
          CALL Fatal(Caller,'Element not associated')
        END IF

        Element % TYPE => GetElementType( elemcode )
        IF(.NOT. ASSOCIATED( Element % TYPE ) ) THEN
          CALL Fatal(Caller,'Could not get element code: '//TRIM(I2S(elemcode)))
        END IF

        n = Element % Type % NumberOfNodes

        Element % GElementIndex = elemindex
        Element % ElementIndex = t

        IF( IsBulk ) THEN
          minelem = MIN( minelem, Element % GElementIndex )
          maxelem = MAX( maxelem, Element % GElementIndex )
        END IF
          
        ! Change the owner partition of the element
        Element % PartIndex = ParEnv % MyPe

        IF( IsBulk ) THEN
          Element % BodyId = geom_id
          icount = icount + 3
        ELSE
          Element % BodyId = 0
          IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) THEN
            ALLOCATE( Element % BoundaryInfo, STAT = allocstat )
            IF( allocstat /= 0 ) THEN
              CALL Fatal(Caller,'Could not allocate boundary info!')
            END IF
          END IF

          Element % BoundaryInfo % Constraint = geom_id

          ! These are the left and right boundary indexes that currently are not used at all!
          LeftParent(t) = PPack % idata(icount+4)
          RightParent(t) = PPack % idata(icount+5)
          icount = icount + 5
        END IF

        ALLOCATE( Element % NodeIndexes(n), STAT = allocstat )
        IF( allocstat /= 0 ) THEN
          CALL Fatal(Caller,'Could not allocate a few node indexes!')
        END IF

        IF( icount + n > SIZE( PPack % idata) ) THEN
          CALL Fatal(Caller,'icount out of range')
        END IF

        Element % NodeIndexes(1:n) = PPack % idata(icount+1:icount+n)

        ! Renumber the nodes such that the local indexes are always contiguous
        DO j=1, n
          k = Element % NodeIndexes(j)
          IF( k < minind .OR. k > maxind ) THEN
            CALL Fatal(Caller,'k3 out of bounds: '//TRIM(I2S(k)))
          END IF
          IF( GlobalToLocal(k) <= 0 .OR. GlobalToLocal(k) > newnodes ) THEN
            CALL Fatal(Caller,'Local index out of bounds')
          END IF
          Element % NodeIndexes(j) = GlobalToLocal(k)
        END DO

        ! Advance the counter for the data
        icount = icount + n
      END DO

      CALL Info(Caller,'Finished piece',Level=20)
      IF( icount /= PPack % indpos ) THEN
        CALL Fatal(Caller,'Inconsistent icount value: '//TRIM(I2S(icount)))
      END IF
    END DO

    IF( minelem <= maxelem ) THEN
      ! First create global to local array for the elements 
      CALL Info(Caller,'Global element index range: '&
          //TRIM(I2S(minelem))//' to '//TRIM(I2S(maxelem)),Level=8)
      ALLOCATE( GlobalToLocalElem(minelem:maxelem))
      GlobalToLocalElem = 0
      DO i = 1, newnbulk
        j =  NewMesh % Elements(i) % GElementIndex
        IF( j >= minelem .AND. j <= maxelem ) THEN
          GlobalToLocalElem( j ) = i
        ELSE
          PRINT *,'j out of bounds:',ParEnv % MyPe, i, j
        END IF
      END DO
      
      ! Then use the temporal vectors to repoint the left and right indexes to elements
      DO i = newnbulk+1, newnbulk + newnbdry
        Element => NewMesh % Elements(i)

        j = LeftParent(i)
        IF( j >= minelem .AND. j <= maxelem ) THEN
          k = GlobalToLocalElem(j)
          IF( k <= 0 .OR. k > newnbulk ) THEN
            PRINT *,'k left out of bounds:',ParEnv % MyPe, i, j, k
          ELSE
            Element % BoundaryInfo % Left => NewMesh % Elements(k)
          END IF
        END IF
        j = RightParent(i)
        IF( j >= minelem .AND. j <= maxelem ) THEN
          k = GlobalToLocalElem(j)
          IF( k <= 0 .OR. k > newnbulk ) THEN
            PRINT *,'k right out of bounds:',ParEnv % MyPe, i, j, k
          ELSE          
            Element % BoundaryInfo % Right => NewMesh % Elements(k)
          END IF
        END IF
      END DO
      DEALLOCATE( GlobalToLocalElem ) 
    END IF
    DEALLOCATE( LeftParent )
    DEALLOCATE( RightParent ) 
    
    
    
    CALL Info(Caller,'Unpacking incoming nodes',Level=20)
    DO part=1,ParEnv % PEs
      IF( part-1 == ParEnv % MyPe ) CYCLE

      PPack => RecPack(part)
      IF( PPack % icount <= 5 ) CYCLE

      icount = PPack % indpos
      rcount = 0
      lcount = 0
      
      DO i = 1, PPack % NumberOfNodes
        icount = icount + 1
        j = Ppack % idata(icount)
        IF( j < minind .OR. j > maxind ) THEN
          CALL Fatal(Caller,'Incoming node index out of bounds')
        END IF
        k = GlobalToLocal( j )

        IF( k <= 0 .OR. k > newnodes ) THEN
          CALL Fatal(Caller,'Local index out of bounds')
        END IF

        NewMesh % Nodes % x(k) = PPack % rdata(rcount+1)
        NewMesh % Nodes % y(k) = PPack % rdata(rcount+2)
        IF( dim == 3 ) NewMesh % Nodes % z(k) = PPack % rdata(rcount+3)
        rcount = rcount + dim

        NewMesh % ParallelInfo % INTERFACE(k) = PPack % ldata(lcount+1)
        lcount = lcount + 1
      END DO
    END DO

    n = COUNT( NewMesh % ParallelInfo % INTERFACE )
    CALL Info(Caller,'Potential interface nodes '//TRIM(I2S(n))//' out of '&
        //TRIM(I2S(NewMesh % NumberOfNodes)),Level=20)

    
    CALL Info(Caller,'Creating local to global numbering for '&
         //TRIM(I2S(newnodes))//' nodes',Level=20)
    ALLOCATE( NewMesh % ParallelInfo % GlobalDofs( newnodes ), STAT = allocstat)
    NewMesh % ParallelInfo % GlobalDofs = 0
    IF( allocstat /= 0 ) THEN
      CALL Fatal(Caller,'Could not allocate global dof indexes')
    END IF

    DO i = minind, maxind
      j = GlobalToLocal(i)
      IF( j == 0 ) CYCLE
      IF( j > newnodes ) THEN
        CALL Fatal(Caller,'Invalid node index')
      END IF
      NewMesh % ParallelInfo % GlobalDofs(j) = i
    END DO

     DO i=1,NewMesh % NumberOfBulkElements
       Element => NewMesh % Elements(i)
       NewMesh % MaxElementNodes = MAX( NewMesh % MaxElementNodes, &
            Element % TYPE % NumberOfNodes)
       NewMesh % MaxElementDOFs = MAX( NewMesh % MaxElementDOFs, &
           Element % TYPE % NumberOfNodes + &
           Element % TYPE % NumberOfEdges * NewMesh % MaxEdgeDOFs + &
           Element % TYPE % NumberOfFaces * NewMesh % MaxFaceDOFs + &
           Element % BDOFs, &
           Element % DGDOFs ) 
     END DO

    CALL Info(Caller,'Finished unpacking and gluing mesh pieces',Level=8)

  END SUBROUTINE UnpackMeshPieces

  
  ! Given a partitioning create a list of potential nodes at the interface.
  ! The list is conservative including all old and possible new nodes.
  !------------------------------------------------------------------------------
  SUBROUTINE UpdateInterfaceNodeCandidates(Mesh)
    TYPE(Mesh_t), POINTER :: Mesh
    
    INTEGER :: i,j,k,n,m,part,allocstat
    TYPE(Element_t), POINTER :: Element
    INTEGER, ALLOCATABLE :: PrevPartition(:)
    INTEGER, POINTER :: ElementPart(:)
    LOGICAL, POINTER :: PartInterface(:)
    CHARACTER(*), PARAMETER :: Caller = "UpdateInterfaceNodeCandidates"
    
    CALL Info(Caller,'Updating the list of potential interface nodes')

    n = Mesh % NumberOfNodes
    IF( n == 0 ) RETURN
    
    IF( .NOT. ASSOCIATED( Mesh % ParallelInfo % INTERFACE ) ) THEN
      ALLOCATE( Mesh % ParallelInfo % INTERFACE( n ), STAT=allocstat )
      IF( allocstat /= 0 ) THEN
        CALL Fatal(Caller,'Allocation error for parallel interface!')
      END IF
      Mesh % ParallelInfo % Interface = .FALSE.
    END IF
    PartInterface => Mesh % ParallelInfo % Interface

    IF( .NOT. ASSOCIATED( Mesh % RePartition ) ) THEN
      CALL Fatal(Caller,'Allocation error for parallel interface!')
    END IF
    ElementPart => Mesh % RePartition
    
    ALLOCATE( PrevPartition( n ), STAT=allocstat ) 
    IF( allocstat /= 0 ) THEN
      CALL Fatal(Caller,'Allocation error for prev partition!')
    END IF
    PrevPartition = 0
    
    DO i=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(i)
      part = ElementPart(i)
      IF( part <= 0 ) CYCLE
      m = Element % TYPE % NumberOfNodes
      DO j=1,m
        k = Element % NodeIndexes(j)
        IF( PrevPartition(k) == 0 ) THEN
          PrevPartition(k) = part
        ELSE IF( PrevPartition(k) /= part ) THEN
          PartInterface(k) = .TRUE.
        END IF
      END DO
    END DO

    n = COUNT( PartInterface ) 
    DEALLOCATE( PrevPartition )
    
    CALL Info(Caller,'Number of potential nodes at the intarface: '//TRIM(I2S(n)),Level=10)      
      
  END SUBROUTINE UpdateInterfaceNodeCandidates

  !> Based on a conservative list of potential interface nodes
  !> in ParallelInfo % Interface, find real interface nodes &
  !> populate NeighbourList % Neighbours
  !-------------------------------------------------------------------
  SUBROUTINE FindRepartitionInterfaces(Model, Mesh, DIM)
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: DIM
    !--------------------------------------
    INTEGER :: i,j,k,n,loc,ierr,nneighs, counter
    INTEGER, ALLOCATABLE :: GDOFsSend(:), GDOFsRecv(:), LocalNNum(:),PartSendCount(:),work_arr(:),&
         requests(:)
    LOGICAL, POINTER :: iface(:)
    LOGICAL, ALLOCATABLE :: PartMaybeNeigh(:)
    CHARACTER(*), PARAMETER :: FuncName = 'FindRepartitionInterfaces'
    TYPE(NeighbourList_t), POINTER :: NeighList(:)
    TYPE RecvNodes_t
       INTEGER, ALLOCATABLE :: GDOFs(:)
       INTEGER :: part,n
    END TYPE RecvNodes_t
    TYPE(RecvNodes_t), ALLOCATABLE :: RecvNodes(:)

    !Find potential neighbour partitions - shouldn't really need to buffer
    PartMaybeNeigh = FindMeshNeighboursGeometric(Mesh,DIM,1.0_dp)

    Iface => Mesh % ParallelInfo % Interface
    NeighList => Mesh % ParallelInfo % NeighbourList

    n = COUNT(Iface)
    nneighs = COUNT(PartMaybeNeigh)
    ALLOCATE(GDOFsSend(n), LocalNNum(n), PartSendCount(ParEnv % PEs),&
         RecvNodes(nneighs), requests(nneighs), work_arr(ParEnv % PEs))

    counter = 0
    DO i=1,Mesh % NumberOfNodes
      IF(.NOT. Iface(i)) CYCLE
      counter = counter + 1
      GDOFsSend(counter) = Mesh % ParallelInfo % GlobalDOFs(i)
      LocalNNum(counter) = i
    END DO
    CALL SortI(n, GDOFsSend, LocalNNum)

    CALL MPI_AllGather(n, 1, MPI_INTEGER, PartSendCount, 1, MPI_INTEGER, ELMER_COMM_WORLD, ierr)
    IF(ierr /= 0) CALL Fatal(FuncName, "MPI Error communicating node send counts")

    !Create receiving structure
    counter = 0
    DO i=1,ParEnv % PEs
      IF(.NOT. PartMaybeNeigh(i)) CYCLE
      IF(i == ParEnv % MyPE + 1) CYCLE
      counter = counter + 1
      RecvNodes(counter) % part = i-1
      RecvNodes(counter) % n = PartSendCount(i)
      ALLOCATE(RecvNodes(counter) % GDOFs(RecvNodes(counter) % n))
    END DO

    !Send sorted gdof lists
    DO i=1,nneighs

      CALL MPI_IRECV( RecvNodes(i) % GDOFs,  RecvNodes(i) % n, MPI_INTEGER, &
           RecvNodes(i) % part, 1000, ELMER_COMM_WORLD, requests(i), ierr)
      IF(ierr /= 0) CALL Fatal(FuncName, "MPI Error receiving GDOFs")
      CALL MPI_SEND( GDOFsSend, n, MPI_INTEGER, RecvNodes(i) % part, 1000, &
           ELMER_COMM_WORLD, ierr)
      IF(ierr /= 0) CALL Fatal(FuncName, "MPI Error sending GDOFs")

    END DO
    CALL MPI_WaitAll( nneighs, requests, MPI_STATUSES_IGNORE, ierr )

    !Fill neighbourlist % neighbours based on received data
    DO i=1,n
      work_arr = -1
      counter = 0

      k = GDOFsSend(i)
      DO j=1,nneighs
        loc = SearchI(RecvNodes(j) % n, RecvNodes(j) % GDOFs, k)
        IF(loc==0) CYCLE !not found in this partition
        counter = counter + 1
        work_arr(counter) = RecvNodes(j) % part
      END DO

      IF(counter == 0) Iface(LocalNNum(i)) = .FALSE.

      !Reallocate neighbourlist if wrong size
      IF(ASSOCIATED(Neighlist(LocalNNum(i)) % Neighbours)) THEN
        IF(SIZE(Neighlist(LocalNNum(i)) % Neighbours) /= counter + 1) THEN
          DEALLOCATE(NeighList(LocalNNum(i)) % Neighbours)
          NULLIFY(NeighList(LocalNNum(i)) % Neighbours)
        END IF
      END IF
      IF(.NOT. ASSOCIATED(NeighList(LocalNNum(i)) % Neighbours)) &
           ALLOCATE(NeighList(LocalNNum(i)) % Neighbours(counter+1))

      !Fill the list
      NeighList(LocalNNum(i)) % Neighbours(1) = ParEnv % MyPE
      NeighList(LocalNNum(i)) % Neighbours(2:counter+1) = work_arr(1:counter)
    END DO

    !Now we cycle all nodes neighbourlists, allocating any missing neighbourlists and
    !filling with our own partition number
    DO i=1,Mesh % NumberOfNodes
      IF(.NOT. ASSOCIATED(NeighList(i) % Neighbours)) THEN
        IF(Iface(i)) CALL Fatal(FuncName, "Programming error: missed interface node!")
        ALLOCATE(NeighList(i) % Neighbours(1))
        NeighList(i) % Neighbours(1) = ParEnv % MyPE
      ELSE IF(SIZE(NeighList(i) % Neighbours) > 1 .AND. .NOT. Iface(i)) THEN
        DEALLOCATE(Neighlist(i) % Neighbours)
        ALLOCATE(NeighList(i) % Neighbours(1))
        NeighList(i) % Neighbours(1) = ParEnv % MyPE
      END IF
    END DO


    ! Order the indexes such that owner is the lowest one
    ! The other ordering does not matter.
    DO i=1,Mesh % NumberOfNodes
      n = SIZE( NeighList(i) % Neighbours )
      DO j=n-1,1,-1
        k = NeighList(i) % Neighbours(j) 
        IF( k > NeighList(i) % Neighbours(j+1) ) THEN
          NeighList(i) % Neighbours(j) = NeighList(i) % Neighbours(j+1)
          NeighList(i) % Neighbours(j+1) = k
        END IF
      END DO
    END DO

    Mesh % ParallelInfo % NumberOfIFDofs = COUNT(Iface)

  END SUBROUTINE FindRepartitionInterfaces

  
  ! Works out potential neighbour partitions based on Mesh % Nodes bounding box
  !-----------------------------------------------------------------------------
  FUNCTION FindMeshNeighboursGeometric(Mesh,DIM,Buffer) RESULT(PartIsNearby)
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: DIM
    REAL(KIND=dp) :: Buffer
    LOGICAL, ALLOCATABLE :: PartIsNearby(:)
    !---------------------------
    REAL(KIND=dp), ALLOCATABLE :: MyBBox(:),BBoxes(:)
    INTEGER :: i,j,ierr, bbsz
    CHARACTER(*), PARAMETER :: FuncName = 'FindMeshNeighboursGeometric'

    bbsz = 2*DIM !size of bounding box

    ALLOCATE(PartIsNearby(ParEnv % PEs),&
         MyBBox(bbsz),&
         BBoxes(bbsz*ParEnv % PEs))
    PartIsNearby = .FALSE.

    IF(Mesh % NumberOfNodes > 0) THEN
      MyBBox(1) = MINVAL(Mesh % Nodes % x)
      MyBBox(2) = MAXVAL(Mesh % Nodes % x)
      IF(DIM >= 2) THEN
        MyBBox(3) = MINVAL(Mesh % Nodes % y)
        MyBBox(4) = MAXVAL(Mesh % Nodes % y)
      END IF
      IF(DIM==3) THEN
        MyBBox(5) = MINVAL(Mesh % Nodes % z)
        MyBBox(6) = MAXVAL(Mesh % Nodes % z)
      END IF
    ELSE
      MyBBox = 0.0_dp
    END IF

    CALL MPI_AllGather(MyBBox, bbsz, MPI_DOUBLE_PRECISION, BBoxes, &
         bbsz, MPI_DOUBLE_PRECISION, ELMER_COMM_WORLD, ierr)
    IF(ierr /= 0) CALL Fatal(FuncName, "MPI Error communicating bounding boxes")

    DO i=0, ParEnv % PEs-1
      IF(i == ParEnv % MyPE) CYCLE
      IF(BBoxes(i*bbsz+1) - buffer > MyBBox(2)) CYCLE
      IF(BBoxes(i*bbsz+2) + buffer < MyBBox(1)) CYCLE
      IF(DIM >= 2) THEN
        IF(BBoxes(i*bbsz+3) - buffer > MyBBox(4)) CYCLE
        IF(BBoxes(i*bbsz+4) + buffer < MyBBox(3)) CYCLE
      END IF
      IF(DIM==3) THEN
        IF(BBoxes(i*bbsz+5) - buffer > MyBBox(6)) CYCLE
        IF(BBoxes(i*bbsz+6) + buffer < MyBBox(5)) CYCLE
      END IF
      PartIsNearby(i+1) = .TRUE.
    END DO

  END FUNCTION FindMeshNeighboursGeometric


  !> Makes a serial mesh partitiong. Current uses geometric criteria.
  !> Includes some hybrid strategies where the different physical domains
  !> are partitioned using different strategies. 
  !----------------------------------------------------------------------------  
  SUBROUTINE PartitionMeshSerial( Model, Mesh, Params ) 
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     TYPE(Mesh_t), POINTER  :: Mesh, ParallelMesh
     TYPE(ValueList_t), POINTER :: Params 
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: SectionParams
     INTEGER, ALLOCATABLE :: ParameterInd(:), ElementSet(:)
     INTEGER, POINTER :: ElementPart(:)
     INTEGER :: NumberOfSets, NumberOfBoundarySets, SetNo, id
     LOGICAL, POINTER :: PartitionCand(:)
     INTEGER :: i,j,k,n, allocstat
     LOGICAL :: Found
     INTEGER, ALLOCATABLE :: EquationPart(:)    

     TYPE(NeighbourList_t),POINTER  :: NeighbourList(:)
     CHARACTER(*), PARAMETER :: FuncName = 'PartitionMeshSerial'
     
     !-----------------------------------------------------------------------
     CALL Info(FuncName,'Using internal mesh partitioning on one processor')
       
     n = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements     

     ALLOCATE( PartitionCand(n), ElementSet(n), STAT = allocstat )
     IF( allocstat /= 0 ) THEN
       CALL Fatal(FuncName,'Allocation error in partition stuff')
     END IF

     IF( ASSOCIATED( Mesh % RePartition ) ) THEN
       IF( SIZE( Mesh % RePartition ) < n ) THEN
         DEALLOCATE(Mesh % RePartition)
         Mesh % RePartition => Null()
       END IF
     END IF

     IF(.NOT. ASSOCIATED( Mesh % RePartition ) ) THEN
       ALLOCATE( Mesh % RePartition( n ), STAT = allocstat)
       IF( allocstat /= 0 ) THEN
         CALL Fatal(FuncName,'Allocation error for repartitioning vector')       
       END IF
     END IF
     ElementPart => Mesh % RePartition 
     
     PartitionCand = .FALSE.
     ElementSet = 0
     ElementPart = 0
     
     n = MAX( Model % NumberOfBCs, Model % NumberOfEquations ) 
     ALLOCATE( ParameterInd(n), STAT = allocstat ) 
     IF( allocstat /= 0 ) THEN
       CALL Fatal(FuncName,'Allocation error for ParameterInd')
     END IF
     
     ParameterInd = 0
     
     CALL Info(FuncName,'Partitioning the boundary elements sets') 
     CALL InitializeBoundaryElementSet(NumberOfBoundarySets)
     
     IF( NumberOfBoundarySets > 0 ) THEN
       DO SetNo = 1, NumberOfBoundarySets
         SectionParams => NULL()
         IF( ParameterInd(SetNo) > 0 ) THEN
           id = ParameterInd(SetNo)
           IF( id <= Model % NumberOfBCs ) THEN
             SectionParams => Model % BCs(id) % Values
           END IF
         END IF
         CALL PartitionMeshPart(SetNo,SectionParams,.TRUE.)
       END DO
     
       CALL InheritBoundaryPart()
       IF( ListGetLogical( Params,'Partition Mesh Merge Boundaries',Found ) ) THEN
         CALL MergeBoundaryPart()
       END IF
       CALL ExtendBoundaryPart()
     END IF
    
     CALL Info(FuncName,'Partition the bulk elements sets')
     CALL InitializeBulkElementSet(NumberOfSets)

     DO SetNo = 1, NumberOfSets
       SectionParams => NULL()
       IF( ParameterInd(SetNo) > 0 ) THEN
         id = ParameterInd(SetNo)
         IF( id <= Model % NumberOfEquations ) THEN
           SectionParams => Model % Equations(id) % Values
         END IF
       END IF
       CALL PartitionMeshPart(SetNo,SectionParams,.FALSE.)
     END DO

     !CALL Info(FuncName,'Defining halo mesh')     
     !CALL DefineMeshHalo()
     
     CALL CreateNeighbourList()
     

100  CALL Info(FuncName,'All done for now')


   CONTAINS

     ! Inherit partition from a boundary partition.
     ! In case of conflict the 1st occurrence prevails.
     !-----------------------------------------------------
     SUBROUTINE InheritBoundaryPart()
       
       TYPE(Element_t), POINTER :: Element, Parent
       INTEGER :: t, LeftRight, BoundPart, NoHerited, NoConflict, ElemIndx

       CALL Info(FuncName,'Inheriting the boundary patitioning into the bulk mesh') 

       NoHerited = 0
       NoConflict = 0

       DO t=Mesh % NumberOfBulkElements + 1,&
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 
         Element => Mesh % Elements(t) 
         IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
           DO LeftRight=0,1
             IF( LeftRight == 0 ) THEN
               Parent => Element % BoundaryInfo % Left
             ELSE
               Parent => Element % BoundaryInfo % Right
             END IF

             IF( ASSOCIATED( Parent ) ) THEN
               BoundPart = ElementSet( t ) 
               ElemIndx = Parent % ElementIndex
               IF( ElementSet( ElemIndx ) == 0 ) THEN
                 NoHerited = NoHerited + 1
                 ElementPart( ElemIndx ) = BoundPart
               ELSE IF( ElementSet( ElemIndx ) /= BoundPart ) THEN
                 NoConflict = NoConflict + 1
               END IF
             END IF
           END DO
         END IF
       END DO
       
       CALL Info(FuncName,'Number of herited bulk elements: '//TRIM(I2S(NoHerited)))
       CALL Info(FuncName,'Number of conflicted bulk elements: '//TRIM(I2S(NoConflict)))
       
     END SUBROUTINE InheritBoundaryPart


     ! Merge partitioned boundaries that share even just one node. 
     ! The algorithm works fine when there is a small number of 
     ! partitions on the boundary and requires minimal additional space.
     !------------------------------------------------------------------
     SUBROUTINE MergeBoundaryPart()
       
       TYPE(Element_t), POINTER :: Element
       INTEGER :: t, i, j, MaxPart
       LOGICAL, ALLOCATABLE :: PartFlag(:), PartitionCoupling(:,:)
       INTEGER, ALLOCATABLE :: PartMap(:)
       

       CALL Info(FuncName,'Inheriting the boundary patitioning into the bulk mesh') 

       MaxPart = MAXVAL( ElementPart ) 

       ALLOCATE( PartitionCoupling(MaxPart, MaxPart) )
       PartitionCoupling = .FALSE.
       
       ALLOCATE( PartFlag( Mesh % NumberOfNodes ) ) 

       DO i = 1, MaxPart
         PartFlag = .FALSE.
         CALL Info(FuncName,'Studying coupling with partition:'//TRIM(I2S(i)),Level=20)
         DO t=1, Mesh % NumberOfBulkElements 
           IF( ElementPart( t ) == i ) THEN 
             Element => Mesh % Elements(t) 
             PartFlag( Element % NodeIndexes ) = .TRUE.
           END IF
         END DO

         ! Studying to which partitions couple to
         ! Only study cases j>i since the coupling is symmetric
         DO t=1, Mesh % NumberOfBulkElements 
           j = ElementPart( t )
           IF( j > i ) THEN
             Element => Mesh % Elements(t) 
             IF( ANY( PartFlag( Element % NodeIndexes ) ) ) THEN
               IF( .NOT. PartitionCoupling(i,j) ) THEN
                 CALL Info(FuncName,&
                     'Coupling '//TRIM(I2S(i))//' and '//TRIM(I2S(j)),Level=10)
                 PartitionCoupling(i,j) = .TRUE.
                 PartitionCoupling(j,i) = .TRUE.
               END IF
             END IF
           END IF
         END DO
       END DO

       IF(.NOT. ANY( PartitionCoupling ) ) THEN
         CALL Info(FuncName,'Partitions are not coupled')
         RETURN
       END IF

       ! Create mapping of existing partitions to new reduced number of partitions
       ALLOCATE( PartMap( MaxPart ) )
       DO i=1,MaxPart
         PartMap(i) = i
       END DO
       DO i=1,MaxPart
         DO j=i+1,MaxPart
           IF( PartitionCoupling(i,j) ) THEN
             IF( PartMap(i) /= PartMap(j) ) THEN
               CALL Info(FuncName,'Mapping partition '&
                   //TRIM(I2S(j))//' to be '//TRIM(I2S(PartMap(i))),Level=8)
               PartMap(j) = PartMap(i)
             END IF
           END IF
         END DO
       END DO
       j = MAXVAL( PartMap ) 
       CALL Info(FuncName,'Number of mapped partitions: '//TRIM(I2S(j)))

       ! The coupling is studied via bulk elements as they are all that matters. 
       ! In the end we also remap the boundary elements for consistancy. 
       DO t=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         i = ElementPart( t ) 
         IF( i > 0 ) THEN
           ElementPart(t) = PartMap(i)
         END IF
       END DO
       CALL Info(FuncName,'Connected boundaries merged')

       
     END SUBROUTINE MergeBoundaryPart




     ! Extend partition from an existing bulk partitioning. 
     ! In case of conflict the dominating partitioning prevails.
     ! The routine is written with just a small number of existing
     ! boundary partitions in mind and uses minimal memory. 
     !------------------------------------------------------------

     SUBROUTINE ExtendBoundaryPart()
       
       TYPE(Element_t), POINTER :: Element
       INTEGER :: t, ExtendLayers, NoExtend, ElemIndx, NoHits, TestPart, NumberOfParts
       LOGICAL, ALLOCATABLE :: ActiveNode(:)
       INTEGER, ALLOCATABLE :: RefHits(:)
       INTEGER :: allocstat

       NoExtend = 0

       ExtendLayers = ListGetInteger( Params,'Partition Mesh Extend Layers', Found ) 
       IF( ExtendLayers <= 0 ) RETURN
       
       CALL Info(FuncName,'Extending boundary meshes by layers: '//TRIM(I2S(ExtendLayers)))

       ALLOCATE( ActiveNode( Mesh % NumberOfNodes ), STAT = allocstat ) 
       IF( allocstat /= 0 ) THEN
         CALL Fatal(FuncName,'Allocation error for ActiveNode')
       END IF


       NumberOfParts = MAXVAL( ElementPart ) 
       IF( NumberOfParts > 1 ) THEN
         CALL Info(FuncName,'Extending boundary to dominating owner among: '//TRIM(I2S(NumberOfParts)))
         ALLOCATE( RefHits( Mesh % NumberOfBulkElements ), STAT = allocstat )
         IF( allocstat /= 0 ) THEN
           CALL Fatal(FuncName,'Allocation error for RefHits')
         END IF
       END IF


       DO i=1,ExtendLayers

         ! If testing for several partitions then nullify the reference
         IF( NumberOfParts > 1 ) THEN
           RefHits = 0
         END IF

         DO TestPart = 1, NumberOfParts         

           ! Set the active nodes for the partition under testing
           ActiveNode = .FALSE.
           DO t=1, Mesh % NumberOfBulkElements 
             IF( ElementPart( t ) == TestPart ) THEN
               Element => Mesh % Elements(t) 
               ActiveNode( Element % NodeIndexes ) = .TRUE.
             END IF
           END DO
           
           ! Count the number of hits for this partition
           ! If larger than the maximum so far set the partition
           ! For just one existing partitioning no checks need to be done. 
           !--------------------------------------------------------------
           DO t=1, Mesh % NumberOfBulkElements 
             IF( ElementPart( t ) /= 0 ) CYCLE

             Element => Mesh % Elements(t) 
             NoHits = COUNT( ActiveNode( Element % NodeIndexes ) )
             IF( NoHits == 0 ) CYCLE

             IF( NumberOfParts > 1 ) THEN
               IF( NoHits <= RefHits( t ) ) CYCLE
               RefHits( t ) = NoHits
             END IF
             
             ElementPart( t ) = TestPart 
             NoExtend = NoExtend + 1
           END DO
         END DO
       END DO

       CALL Info(FuncName,'Number of extended bulk elements: '//TRIM(I2S(NoExtend)))
       
       DEALLOCATE( ActiveNode ) 
       IF( NumberOfParts > 1 ) DEALLOCATE( RefHits ) 

     END SUBROUTINE ExtendBoundaryPart

      

     ! Initialize sets of boundary elements to be partitioned with various strategies
     !--------------------------------------------------------------------------
     SUBROUTINE InitializeBoundaryElementSet( NumberOfParts )
       
       INTEGER :: NumberOfParts
       
       INTEGER :: i,j,k,bc_id
       TYPE(ValueList_t), POINTER :: ValueList
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: Found, SeparateBoundarySets        
       INTEGER, ALLOCATABLE :: BCPart(:)
       INTEGER :: allocstat

       
       SeparateBoundarySets = ListGetLogical( Params, &
           'Partitioning Separate Boundary Set', Found)
       
       ALLOCATE( BCPart( Model % NUmberOfBCs ), STAT = allocstat )
       IF( allocstat /= 0 ) THEN
         CALL Fatal(FuncName,'Allocation error for BCPart')
       END IF
       BCPart = 0

       NumberOfParts = 0
       DO bc_id = 1, Model % NumberOfBCs 
         ValueList => Model % BCs(bc_id) % Values
         k = ListGetInteger( ValueList,'Partition Set',Found)
         IF( .NOT. Found ) CYCLE
         BCPart(bc_id) = k 
         ParameterInd(k) = bc_id
       END DO
       

       j = 0
       DO WHILE( ANY( BCPart == j+1 ) )
         j = j + 1
       END DO
       IF( j == 0 ) RETURN


       DO bc_id = 1, Model % NumberOfBCs 
         ValueList => Model % BCs(bc_id) % Values
         
         IF( BCPart(bc_id) > 0 ) CYCLE
         
         IF( ListGetLogical( ValueList, 'Discontinuous Boundary', Found) ) THEN
           BCPart(bc_id) = j
           ParameterInd(j) = bc_id
         END IF
         
         k = ListGetInteger( ValueList, 'Periodic Boundary',Found) 
         IF( k > 0 ) THEN
           BCPart(bc_id) = j
           BCPart(k) = j
           ParameterInd(j) = bc_id
         END IF
         
         k = ListGetInteger( ValueList, 'Mortar Boundary',Found) 
         IF( k > 0 ) THEN
           BCPart(bc_id) = j
           BCPart(k) = j
           ParameterInd(j) = bc_id
         END IF
         
         IF( SeparateBoundarySets .AND. BCPart(bc_id) > 0 ) THEN
           DO WHILE( ANY( BCPart == j ) .OR. ANY( EquationPart == j ) ) 
             j = j + 1
           END DO
         END IF
       END DO
      
       NumberOfParts = MAXVAL( BCPart ) 
       
       IF( NumberOfParts == 0 ) RETURN
       

       DO i = Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         Element => Mesh % Elements(i)
         
         DO bc_id=1,Model % NumberOfBCs
           IF ( Element % BoundaryInfo % Constraint == &
               Model % BCs(bc_id) % Tag ) THEN
             IF( BCPart( bc_id ) > 0 ) THEN
               ElementSet( i ) = BCPart( bc_id )
             END IF
             EXIT
           END IF
         END DO
       END DO

       CALL Info(FuncName,'Number of sets for boundary partitioning: '//TRIM(I2S(NumberOfParts)))
       CALL Info(FuncName,'Number of boundary elements set: '//TRIM(I2S(COUNT(ElementSet > 0))))
      
     END SUBROUTINE InitializeBoundaryElementSet


     
     ! Initialize sets of bulk elements to be partitioned with various strategies
     ! By default there is only one set but certain BCs are treated differently.
     !--------------------------------------------------------------------------
     SUBROUTINE InitializeBulkElementSet( NumberOfParts )
       
       INTEGER :: NumberOfParts
       
       INTEGER :: i,j,k,eq_id, bc_id
       TYPE(ValueList_t), POINTER :: ValueList
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: Found, SeparateBoundarySets, FoundAny 
       INTEGER :: allocstat
     
       ElementSet = 0 
       
       ALLOCATE( EquationPart( Model % NumberOfEquations ), STAT = allocstat ) 
       IF( allocstat /= 0 ) THEN
         CALL Fatal(FuncName,'Allocation error for EquationPart')
       END IF
 
       EquationPart = 0
       FoundAny = .FALSE.
       DO eq_id = 1, Model % NumberOfEquations 
         ValueList => Model % Equations(eq_id) % Values
         k = ListGetInteger( ValueList,'Partition Set',Found) 
         IF( k > 0 ) THEN
           EquationPart(eq_id) = k 
           FoundAny = .TRUE.
         END IF
       END DO
       
       NumberOfParts = MAXVAL( EquationPart ) 
      
       Found = .FALSE.
       DO i = 1, Mesh % NumberOfBulkElements 
         IF( ElementPart(i) > 0 ) CYCLE

         Element => Mesh % Elements(i)
         j = 0
         IF( NumberOfParts > 0 ) THEN
           eq_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values, &
               'Equation', Found )
           IF( eq_id > 0 ) j = EquationPart( eq_id ) 
         END IF
         
         IF( j == 0 ) THEN
           Found = .TRUE.
           ElementSet(i) = NumberOfParts + 1
         ELSE
           ElementSet( i ) = j
         END IF
       END DO
       IF( Found ) NumberOfParts = NumberOfParts + 1
       
       CALL Info(FuncName,'Number of sets for bulk partitioning: '//TRIM(I2S(NumberOfParts)))
       
     END SUBROUTINE InitializeBulkElementSet


     
     ! Partition the nodes that have the correct ElementSet using various strategies.
     !------------------------------------------------------------------------
     SUBROUTINE PartitionMeshPart(SetNo, LocalParams, IsBoundary )
       
      INTEGER :: SetNo
      LOGICAL :: IsBoundary 
      TYPE(ValueList_t), POINTER :: LocalParams       

      LOGICAL :: BoundaryPart
      CHARACTER(LEN=MAX_NAME_LEN) :: CoordTransform, SetMethod
      LOGICAL :: GotCoordTransform, SetNodes
      INTEGER :: SumPartitions, NoPartitions, NoCand
      LOGICAL :: Found
      INTEGER :: NoCandElements
      REAL(KIND=dp) :: BoundaryFraction
      INTEGER :: PartOffset

      
      PartitionCand = ( ElementSet == SetNo )
      n = Mesh % NumberOfBulkElements
      
      NoCandElements = COUNT( PartitionCand ) 


      CALL Info(FuncName,'Doing element set: '//TRIM(I2S(SetNo)))
      CALL Info(FuncName,'Number of elements in set: '//TRIM(I2S(NoCandElements)))

      IF( NoCandElements == 0 ) RETURN

      
      BoundaryFraction = ListGetCReal( Params,&
          'Boundary Partitioning Maximum Fraction',Found)
      IF( IsBoundary .AND. NoCandElements <= &
          BoundaryFraction * Mesh % NumberOfBulkElements ) THEN
        WRITE(Message,'(A,ES12.3)') 'Number of boundary elements below critical limit: ',BoundaryFraction
        CALL Info(FuncName,Message )
        WHERE( PartitionCand ) ElementPart = SetNo
        RETURN
      END IF
         
      BoundaryPart = .FALSE.

      Found = .FALSE.
      IF( ASSOCIATED( LocalParams) ) THEN
        SetMethod = ListGetString( LocalParams,'Partitioning Method',Found)
      END IF
      IF(.NOT. Found ) THEN
        IF( IsBoundary ) THEN
          SetMethod = ListGetString( Params,'Boundary Partitioning Method',Found)
        END IF
        IF(.NOT. Found ) SetMethod = ListGetString( Params,'Partitioning Method',Found)
      END IF
      IF( Found ) THEN
        CALL Info(FuncName,'Using partition method: '//TRIM(SetMethod))
      ELSE
        CALL Fatal(FuncName,'Could not define > Partitioning Method < ')
      END IF

      Found = .FALSE.
      IF( ASSOCIATED( LocalParams) ) THEN
        CoordTransform = ListGetString( LocalParams,&
            'Partitioning Coordinate Transformation',Found)
      END IF
      IF(.NOT. Found ) THEN
        IF( IsBoundary ) THEN
          CoordTransform = ListGetString( Params,&
              'Boundary Partitioning Coordinate Transformation',Found)
        ELSE
          CoordTransform = ListGetString( Params,&
              'Partitioning Coordinate Transformation',Found)
        END IF
      END IF
      GotCoordTransform = Found

      IF( ASSOCIATED( LocalParams) ) THEN
        SetNodes = ListGetLogical( LocalParams,'Partition Nodes',Found )
      END IF
      IF(.NOT. Found ) THEN
        IF( IsBoundary ) THEN
          SetNodes = ListGetLogical( Params,'Boundary Partition Nodes',Found)
        ELSE
          SetNodes = ListGetLogical( Params,'Partition Nodes',Found)
        END IF
      END IF

      IF( ASSOCIATED( LocalParams) ) THEN
        NoPartitions = ListGetInteger( LocalParams,'Number of Partitions',Found )
      END IF
      IF(.NOT. Found ) THEN
        IF( IsBoundary ) THEN
          NoPartitions = ListGetInteger( Params,'Boundary Number of Partitions',Found)
        ELSE
          NoPartitions = ListGetInteger( Params,'Number Of Partitions',Found)
        END IF
      END IF


      ! There may be various coordinate transformation (e.g. to cylindrical coordinates)
      ! that allow for different partitions when using the geometries partitioning routines. 
      IF( GotCoordTransform ) THEN
        CALL CoordinateTransformation( Mesh, CoordTransform, Params, &
            IrreversibleTransformation = .FALSE. )
      END IF

      CALL Info(FuncName,'Using partitioning method: '//TRIM(SetMethod))
      
      PartOffset = MAXVAL( ElementPart )       
      CALL Info(FuncName,'Partitioning offset: '//TRIM(I2S(PartOffset)))

      SELECT CASE( SetMethod ) 
        
        !CASE( 'metis recursive' ) 
        !CASE( 'metis kway' ) 
        !CASE( 'metis nodal' ) 
        !CASE( 'metis dual' ) 

        CASE( 'directional')
          IF( SetNodes ) THEN
            CALL ClusterNodesByDirection(Params,&
                Mesh,ElementPart,PartitionCand)
          ELSE
            CALL ClusterElementsByDirection(Params,&
                Mesh,ElementPart,PartitionCand)
          END IF

        CASE( 'uniform' )          
          CALL ClusterElementsUniform(Params,&
              Mesh,ElementPart,PartitionCand)
          
        CASE( 'zoltan' )
#ifdef HAVE_ZOLTAN
          IF( SetNo /= 1 .OR. IsBoundary ) THEN
            CALL Fatal('PartitionMeshPart','Zoltan interface not applicable to hybrid partitioning!')
          END IF
          CALL Zoltan_Interface( Model, Mesh, SerialMode = .TRUE. )
#else
          CALL Fatal(FuncName,'Partition with Zoltan not available!')
#endif 
          
        CASE DEFAULT
          CALL Fatal(FuncName,'Unspecificed partitioning: '//TRIM(SetMethod))
          
      END SELECT

      IF( PartOffset > 0 ) THEN
        WHERE( PartitionCand ) ElementPart = ElementPart + PartOffset
      END IF


      CALL Info(FuncName,'Partitioning of set finished')

      IF( GotCoordTransform ) THEN
        CALL BackCoordinateTransformation( Mesh, DeleteTemporalMesh = .TRUE. )
      END IF

  
    END SUBROUTINE PartitionMeshPart
      

    ! Given a partitioning create a list of Neighbours needed for the communication
    !------------------------------------------------------------------------------
    SUBROUTINE CreateNeighbourList()

      INTEGER :: i,j,k,l,n,m,Partition,lsum,lmax
      INTEGER :: TmpNeighbours(100)
      TYPE(Element_t), POINTER :: Element
      INTEGER :: allocstat

      CALL Info(FuncName,'Creating neighbour list for parallel saving')

      n = Mesh % NumberOfNodes
      ALLOCATE( NeighbourList(n) , STAT=allocstat ) 
      IF( allocstat /= 0 ) THEN
        CALL Fatal(FuncName,'Allocation error for NeighbourList')
      END IF

      DO i=1,n
        NULLIFY( NeighbourList(i) % Neighbours )
      END DO
      
      DO i=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(i)
        Partition = ElementPart(i)
        m = Element % TYPE % NumberOfNodes
        DO j=1,m
          k = Element % NodeIndexes(j)
          IF( .NOT. ASSOCIATED( NeighbourList(k) % Neighbours ) ) THEN
            ALLOCATE( NeighbourList(k) % Neighbours(1), STAT = allocstat )
            IF( allocstat /= 0 ) THEN
              CALL Fatal(FuncName,'Allocation error for Neighbours')
            END IF            
            NeighbourList(k) % Neighbours(1) = Partition
          ELSE IF( .NOT. ANY( NeighbourList(k) % Neighbours == Partition ) ) THEN
            l = SIZE( NeighbourList(k) % Neighbours )

            TmpNeighbours(1:l) = NeighbourList(k) % Neighbours(1:l)
            DEALLOCATE( NeighbourList(k) % Neighbours )

            ALLOCATE( NeighbourList(k) % Neighbours(l+1), STAT = allocstat )
            IF( allocstat /= 0 ) THEN
              CALL Fatal(FuncName,'Allocation error for Neighbours')
            END IF                       
            NeighbourList(k) % Neighbours(1:l) = TmpNeighbours(1:l)
            NeighbourList(k) % Neighbours(l+1) = Partition
          END IF
        END DO
      END DO

      lmax = 0
      lsum = 0
      DO k=1,n
        l = SIZE(NeighbourList(k) % Neighbours)
        lmax = MAX( lmax, l )
        lsum = lsum + l
      END DO
      
      CALL Info(FuncName,'Maximum number of partitions for a node: '//TRIM(I2S(lmax)))
      
      WRITE(Message,'(A,F8.3)') 'Average number of partitiones for a node: ',1.0_dp*lsum/n
      CALL Info(FuncName,Message) 
      
    END SUBROUTINE CreateNeighbourList

  END SUBROUTINE PartitionMeshSerial
  
END MODULE MeshPartition
