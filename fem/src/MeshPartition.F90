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
  USE ClusteringMethods
  
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
  SUBROUTINE Zoltan_Interface( Model, Mesh, SerialMode, NoPartitions, PartitionCand, &
                                StartImbalanceTol, TolChange, MinElems )

    USE MeshUtils

#ifdef HAVE_ZOLTAN
    USE Zoltan
#endif

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, OPTIONAL :: SerialMode
    INTEGER, OPTIONAL :: NoPartitions, MinElems
    LOGICAL, POINTER, OPTIONAL :: PartitionCand(:)
    REAL(KIND=dp), OPTIONAL :: StartImbalanceTol, TolChange
    !------------------------

#ifdef HAVE_ZOLTAN
    TYPE(Element_t), POINTER :: Element
    TYPE(Graph_t) :: LocalGraph
    REAL(KIND=dp) :: t1,t2, ImbalanceTol, ZTolChange
    INTEGER :: i,j,k,l,m,n,ierr,NNodes,NBulk,Ngraph,counter,DIM,&
         max_elemno,NoPart, ZMinElems,OutputLevel
    INTEGER, ALLOCATABLE :: ElemAdj(:), ElemStart(:), ElemAdjProc(:), ParElemAdj(:), ParElemStart(:),&
         ParElemIdx(:),ParElemAdjProc(:),sharecount(:),&
         ParElemMap(:)
    INTEGER, ALLOCATABLE :: PartitionPerm(:), InvPerm(:)
    LOGICAL :: UsePerm, Success, GotParMetis, DistributedMesh
    LOGICAL, ALLOCATABLE :: PartSuccess(:), PartGotNodes(:)
    
    CHARACTER(MAX_NAME_LEN) :: ImbTolStr,Method, Approach, ParMetisLib, ZoltanLib, GraphPackage
    CHARACTER(*), PARAMETER :: FuncName="Zoltan_Interface"

    !Zoltan things
    TYPE(Zoltan_Struct), POINTER :: zz_obj
    INTEGER(Zoltan_INT) :: zierr,numGIDEntries, numLidEntries, numImport, numExport
    INTEGER(Zoltan_INT),DIMENSION(:), POINTER :: importGlobalGids, importLocalGids, importProcs, &
         importToPart,exportGlobalGids,exportLocalGids, exportProcs, exportToPart
    REAL(Zoltan_FLOAT) :: version
    LOGICAL :: changes,Debug=.FALSE.,Serial,Found

    TYPE ElemTable_t
       INTEGER :: counter=0
       INTEGER, ALLOCATABLE :: Idx(:)
    END TYPE ElemTable_T
    TYPE(ElemTable_t),ALLOCATABLE :: NodeElems(:),ElemElems(:)

    TYPE(ValueList_t), POINTER :: PartParams
    TYPE(ValueListEntry_t), POINTER :: ptr
    INTEGER :: ncopy

#ifdef HAVE_PARMETIS
    GotParMetis = .TRUE.
#else
    GotParMetis = .FALSE.
#endif

    IF(PRESENT(StartImbalanceTol)) THEN
      ImbalanceTol = StartImbalanceTol
    ELSE
      CALL Info(FuncName, 'No imbalance tolerance given so starting with 1.1 (10%)')
      ImbalanceTol = 1.1_dp
    END IF

    IF(PRESENT(TolChange)) THEN
      ZTolChange = TolChange
    ELSE
      CALL Info(FuncName, 'No imbalance tolerance change given so using a 0.02 reduction')
      ZTolChange = 0.02_dp
    END IF

    IF(PRESENT(MinElems)) THEN
      ZMinElems = MinElems
    ELSE
      CALL Info(FuncName, 'No min elements required for every partition so setting to 0')
      ZMinElems = 0
    END IF

    CALL Info(FuncName,'Calling Zoltan for mesh partitioning',Level=10)
    PartParams => Model % Simulation

    IF( PRESENT( NoPartitions ) ) THEN
      NoPart = NoPartitions
    ELSE
      NoPart = ListGetInteger( PartParams,'Number of Partitions',Found ) 
      IF(.NOT. Found) NoPart = ParEnv % PEs
    END IF

    CALL Info(FuncName,'Partitioning with Zoltan to '//TRIM(I2S(NoPart))//' partitions',Level=6)
    
    IF( PRESENT( SerialMode ) ) THEN
      Serial = SerialMode
      IF( Serial .AND. ParEnv % PEs >1 ) THEN
        CALL Info(FuncName,'Using serial partitioning in parallel case!')
      END IF
    ELSE
      Serial = ( ParEnv % PEs == 1 )
      IF( PRESENT( PartitionCand ) ) THEN
        CALL Fatal(FuncName,'Masked partitioning not implemented yet in parallel!')
      END IF
    END IF

    IF( NoPart <= 1 ) THEN
      CALL Info(FuncName,'Nothing to do without any partitions requested!')
      RETURN
    END IF
          
    NNodes = Mesh % NumberOfNodes
    NBulk = Mesh % NumberOfBulkElements
    Ngraph = Nbulk
    DIM = CoordinateSystemDimension()
    
    IF(.NOT. Serial) THEN
      ALLOCATE(PartGotNodes(ParEnv % PEs))
      CALL MPI_ALLGATHER(NNodes > 0, 1, MPI_LOGICAL, PartGotNodes, &
          1, MPI_LOGICAL, ELMER_COMM_WORLD, ierr)

      DistributedMesh = ALL(PartGotNodes)
    END IF

    Method = ListGetString(Model % Solver % Values,"Repartition Method", Found)
    IF(.NOT. Found) THEN
      CALL Info(FuncName, "Not Found 'Repartition Method' so assuming 'Zoltan'")
      Method = 'Zoltan'
    END IF
    Approach = ListGetString(Model % Solver % Values,"Repartition Approach", Found)
    IF(.NOT. Found) THEN
      CALL Info(FuncName, "Not Found 'Repartition Approach' so assuming 'repartition'")
      Approach = 'repartition'
    END IF
    OutputLevel = ListGetInteger(Model % Solver % Values,"Repartition Output Level", Found)
    IF(.NOT. Found) THEN
      CALL Info(FuncName, "Not Found 'Repartition Output Level' so assuming '0'")
      OutputLevel = 0
    END IF

    IF(Serial) THEN
      Approach = 'partition'
      Method = 'zoltan'
      CALL Info(FuncName, 'used in serial so using Zoltan partition')
    END IF

    SELECT CASE( Method )
    CASE( 'parmetis' )
      IF(.NOT. GotParMetis) CALL FATAL(FuncName, "ParMetis repartitioning selected but not installed")
      IF(.NOT. DistributedMesh) CALL FATAL(FuncName, "ParMetis requires fully distributed mesh with nodes on each partition")
      ParMetisLib = ListGetString(Model % Solver % Values,"Repartition ParMetis Library", Found)
      IF(.NOT. Found) THEN
        CALL Info(FuncName, "Not Found 'Repartition ParMetis Library' so assuming 'adaptiverepart'")
        ParMetisLib = 'adaptiverepart'
      END IF
    CASE( 'zoltan' )
      ZoltanLib = ListGetString(Model % Solver % Values,"Repartition Zoltan Library", Found)
      IF(.NOT. Found) THEN
        CALL Info(FuncName, "Not Found 'Repartition Zoltan Library' so assuming 'graph'")
        ZoltanLib = 'graph'
      END IF
      GraphPackage = ListGetString(Model % Solver % Values,"Repartition Zoltan Graph Package", Found)
      IF(.NOT. Found) THEN
        CALL Info(FuncName, "Not Found 'Repartition Zoltan Graph Package' so assuming 'phg'")
        GraphPackage = 'phg'
      END IF
    CASE DEFAULT
      CALL Fatal(FuncName,"Repartition method selected invalid")
    END SELECT

10  CONTINUE

    ! If we have a masked partitioning then make a reordering of the bulk elements
    UsePerm = PRESENT( PartitionCand ) 
    IF( UsePerm ) THEN
      n = COUNT( PartitionCand(1:Nbulk) )
      IF( n == Nbulk ) THEN
        UsePerm = .FALSE.
        CALL Info(FuncName,'Candidate list is full, no need for permutation',Level=10)
      ELSE
        CALL Info(FuncName,'Candidate list number of elements: '//I2S(n),Level=10)      
      END IF
    END IF
    
    IF( UsePerm ) THEN
      ALLOCATE( PartitionPerm( NBulk ) )
      PartitionPerm = 0
      WHERE( PartitionCand(1:NBulk) ) PartitionPerm = 1
      j = 0
      DO i=1,NBulk
        IF( PartitionPerm(i)>0 ) THEN
          j = j + 1
          PartitionPerm(i) = j
        END IF
      END DO
      Ngraph = j
      CALL Info(FuncName,'Number of active elements in partitioning:'//I2S(Ngraph),Level=8)

      ALLOCATE( InvPerm( nGraph ) )
      DO i=1,NBulk
        j = PartitionPerm(i)
        IF(j>0) InvPerm(j) = i
      END DO      
    END IF

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

    CALL ListAddNewString( PartParams,"zoltan: debug_level",I2S(OutputLevel))

    SELECT CASE( Method )
    CASE( 'parmetis' )
      CALL ListAddNewString( PartParams,"zoltan: lb_method","graph")
      CALL ListAddNewString( PartParams,"zoltan: graph_package",TRIM(Method))
      CALL ListAddNewString( PartParams,"zoltan: lb_approach",TRIM(Approach))
      CALL ListAddNewString( PartParams,"zoltan: parmetis_method",TRIM(ParMetisLib))
    CASE( 'zoltan' )
      CALL ListAddNewString( PartParams,"zoltan: lb_method",TRIM(ZoltanLib))
      CALL ListAddNewString( PartParams,"zoltan: graph_package",TRIM(GraphPackage))
      CALL ListAddNewString( PartParams,"zoltan: lb_approach",TRIM(Approach))
      CALL ListAddNewString( PartParams,"zoltan: num_gid_entries","1")
      CALL ListAddNewString( PartParams,"zoltan: num_lid_entries","1")
      CALL ListAddNewString( PartParams,"zoltan: obj_weight_dim","0")
      CALL ListAddNewString( PartParams,"zoltan: edge_weight_dim","0")
      CALL ListAddNewString( PartParams,"zoltan: check_graph","0")
      CALL ListAddNewString( PartParams,"zoltan: phg_multilevel","1")
      IF(OutputLevel > 4) CALL ListAddNewString( PartParams,"zoltan: phg_output_level","2")
    CASE DEFAULT
      CALL Fatal(FuncName,"Programming error...")
    END SELECT

    !CALL ListAddNewString( PartParams,"zoltan: imbalance_tol","1.1") !Max load imbalance (default 10%)
    WRITE(ImbTolStr, '(F20.10)') ImbalanceTol
    CALL ListAddNewString( PartParams,"zoltan: imbalance_tol",TRIM(ImbTolStr))

    ! The settings for serial vs. parallel operation differ slightly
    IF( Serial ) THEN
      CALL ListAddNewString( PartParams,"zoltan: return_lists","export part")
      CALL ListAddNewString( PartParams,"zoltan: num_global_parts",TRIM(I2S(NoPart)))  
    ELSE
      CALL ListAddNewString( PartParams,"zoltan: return_lists","all")    !TODO - we only use export list
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
      CALL Info(FuncName,'Succefully set '//I2S(ncopy)//' keywords in zoltan library',Level=8)
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
      IF( UsePerm ) THEN
        CALL LocalElemAdjacency( Mesh, ElemAdj, ElemAdjProc, ElemStart, DIM, &
            PartitionPerm )
      ELSE
        CALL LocalElemAdjacency( Mesh, ElemAdj, ElemAdjProc, ElemStart, DIM )
      END IF
    ELSE
      CALL GlobalElemAdjacency( Mesh, ElemAdj, ElemAdjProc, ElemStart, DIM )
    END IF
      
    numGidEntries = 1
    numLidEntries = 1

    CALL Info(FuncName,'Going into Zoltan partitioning',Level=10)
    zierr = Zoltan_LB_Partition(zz_obj, changes, numGidEntries, numLidEntries, &
         numImport, importGlobalGids, importLocalGids, importProcs, importToPart, &
         numExport, exportGlobalGids, exportLocalGids, exportProcs, exportToPart)
    IF(zierr /= 0) CALL Fatal(FuncName,"Error computing partitioning in Zoltan")

    ! need to check all partitions are given elements
    Success = .TRUE.
    IF(NBulk - numExport + numImport <= ZMinElems) THEN
      WRITE(Message, '(i0,A)') ParEnv % MyPE,' Part not given any elements'
      CALL WARN(FuncName, Message)
      Success = .FALSE.
    END IF

        
    IF(ASSOCIATED(Mesh % Repartition)) THEN
      IF( SIZE( Mesh % Repartition ) < NBulk ) DEALLOCATE(Mesh % Repartition)
    END IF

    ! In hybrid partitioning we may herit partitioning that already has some part set.
    IF(.NOT. ASSOCIATED(Mesh % Repartition) ) THEN
      ALLOCATE(Mesh % Repartition(NBulk))
    END IF

    ! By default stay on this proc, only moving elements are returned    
    IF( UsePerm ) THEN
      Mesh % Repartition(InvPerm) = ParEnv % MyPe + 1
    ELSE
      Mesh % Repartition(1:Nbulk) = ParEnv % MyPE + 1
      ! To be on the safe side unset the BC partitions.
      ! The boundary elements will follow the bulk. 
      i = SIZE( Mesh % RePartition )
      IF( i > Nbulk ) Mesh % Repartition(Nbulk+1:i) = 0
    END IF
      
    IF( numExport > 0 ) THEN
      i = MINVAL( exportLocalGids )
      j = MAXVAL( exportLocalGids )
      IF( i <= 0 .OR. j > NGraph ) THEN
        CALL Fatal(FuncName,'Bad local ID range: '//I2S(i)//' to '//I2S(j))
      END IF
    END IF
      
    IF( UsePerm ) THEN
      DO i=1,numExport
        j = InvPerm(exportLocalGids(i))
        Mesh % Repartition(j) = exportToPart(i) + 1
      END DO
    ELSE IF( Serial ) THEN
      Mesh % Repartition(exportLocalGids(1:numExport)) = exportToPart(1:numExport) + 1
    ELSE
      Mesh % Repartition(exportLocalGids(1:numExport)) = exportProcs(1:numExport) + 1
    END IF

    ! check the success of the rebalancing. Does every partition have nodes?
    ! if not retry with stricter imbalance tolerance
    IF(.NOT. Serial) THEN
      ALLOCATE(PartSuccess(ParEnv % PEs))
      CALL MPI_ALLGATHER(Success, 1, MPI_LOGICAL, PartSuccess, 1, MPI_LOGICAL,&
          ELMER_COMM_WORLD, ierr)
      IF(ANY(.NOT. PartSuccess)) THEN
        ImbalanceTol = ImbalanceTol - ZTolChange
        WRITE(Message, '(A,F10.2)') 'Retrying rebalancing using stricter imbalance tolerance: ', ImbalanceTol
        CALL Info(FuncName, Message)
        IF(ImbalanceTol < 1.0) CALL FATAL(FuncName, 'Unable to rebalance successfully')
        DEALLOCATE(ElemAdj, ElemAdjProc, ElemStart, PartSuccess)
        GOTO 10
      END IF
    END IF
    
    CALL Info(FuncName,'Finished Zoltan partitioning',Level=10)
    
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
      zoltNumObjs = Ngraph
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
      INTEGER :: i,j
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i= 1, NBulk
!        global_ids(i) = Mesh % Elements(i) % GElementIndex
        j = i
        IF( UsePerm ) THEN
          j = PartitionPerm(i)
          IF( j == 0 ) CYCLE
        END IF
        ! global_ids(j) = j 
        global_ids(j) = Mesh % Elements(j) % GElementIndex
 
       local_ids(j) = j
      end do

      !PRINT *,'zoltGetObjs:',ParEnv % MyPe, Ngraph, MINVAL( local_ids(1:Ngraph)), MAXVAL( local_ids(1:Ngraph))
      
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

  !> Returns a CRS dual graph of element face (3D) or edge (2D) connections, including
  !> across partitions. ElemAdj contains the global element numbers of connected elements,
  !> ElemStart describes the CRS positions of each elem, and ElemAdjProc contains the partition
  !> of the connected element.
  !---------------------------------------------------------------------------------------------
  SUBROUTINE GlobalElemAdjacency( Mesh, ElemAdj, ElemAdjProc, ElemStart, DIM )
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, ALLOCATABLE :: ElemAdj(:),ElemStart(:),ElemAdjProc(:)
    INTEGER :: DIM
    !-------------------------------------
    TYPE(Element_t), POINTER :: MFacePtr(:), Element
    INTEGER :: i,j,k,m,n,max_elfaces,el1,el2,gface_id, gpar_id,gpar_lid,ierr,counter,&
         NBulk,NFaces,Sweep,NIFFaces,work_size
    INTEGER, ALLOCATABLE :: ElemConn(:,:), ElemConnPart(:,:), NElConn(:), FaceIFIDX(:),status(:),&
         work_int(:)
    INTEGER, POINTER :: ElFaceIdx(:)
    TYPE(NeighbourList_t), POINTER :: MFaceIFList(:)
    LOGICAL, POINTER :: MFaceIF(:)
    CHARACTER(*), PARAMETER :: FuncName="GlobalElemAdjacency"
    TYPE FaceShare_t
       INTEGER :: count
       INTEGER, ALLOCATABLE :: GFaceIDX(:), GParIDX(:),GParLIDX(:)
    END TYPE FaceShare_t
    TYPE(FaceShare_t), ALLOCATABLE :: SendFaces(:),RecvFaces(:)

    NBulk = Mesh % NumberOfBulkElements
    
    !Find and globally number mesh faces
    IF(Nbulk == 0 ) THEN
      CONTINUE

    ELSE IF(DIM == 3) THEN
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
    !don't know at this point the required size of work_int, if there are lots
    !of bulk elements, probably NBulk is sufficient, but for only a few
    !disconnected elements, maybe not. So set min size = 1000
    work_size = MAX(NBulk, 1000)
    ALLOCATE(SendFaces(ParEnv % PEs),RecvFaces(ParEnv % PEs),work_int(work_size))

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
    DO i=1,Nbulk
      ElemAdj(ElemStart(i):ElemStart(i) + NElConn(i) -1) = &
           ElemConn(1:NElConn(i),i)
      ElemAdjProc(ElemStart(i):ElemStart(i) + NElConn(i) -1) = &
           ElemConnPart(1:NElConn(i),i)
      ElemStart(i+1) = ElemStart(i) + NElConn(i)
    END DO

    CALL MPI_BARRIER(ELMER_COMM_WORLD,ierr)
  END SUBROUTINE GlobalElemAdjacency



  !> As the previous routine except intended for serial meshes without need for communication.
  !---------------------------------------------------------------------------------------------
  SUBROUTINE LocalElemAdjacency( Mesh, ElemAdj, ElemAdjProc, ElemStart, DIM, &
      PartitionPerm )
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, ALLOCATABLE :: ElemAdj(:),ElemStart(:),ElemAdjProc(:)
    INTEGER :: DIM
    INTEGER, OPTIONAL, ALLOCATABLE :: PartitionPerm(:)
    !-------------------------------------
    TYPE(Element_t), POINTER :: MFacePtr(:), Element, Left, Right
    INTEGER :: i,j,k,m,n,max_elfaces,el1,el2,gface_id, gpar_id,gpar_lid,ierr,&
         NBulk,Ngraph,NtotCon,NFaces,Sweep,NIFFaces,condim
    INTEGER, ALLOCATABLE :: ElemConn(:,:), ElemConnPart(:,:), NElConn(:), status(:)
    INTEGER, POINTER :: ElFaceIdx(:)
    TYPE(NeighbourList_t), POINTER :: MFaceIFList(:)
    LOGICAL, POINTER :: MFaceIF(:)
    CHARACTER(*), PARAMETER :: FuncName="LocalElemAdjacency"

    NBulk = Mesh % NumberOfBulkElements
    IF( NBulk == 0 ) RETURN
    
    !Find and globally number mesh faces
!   IF( .TRUE.) THEN
!     CALL FindMeshEdges(Mesh)
!   ELSE IF(DIM == 3) THEN
!     CALL FindMeshFaces3D(Mesh)
!     CALL FindMeshEdges2D(Mesh)
!   ELSE IF(DIM == 2 ) THEN
!     CALL FindMeshEdges2D(Mesh)
!   ELSE
!     CALL Fatal(FuncName,"Not implemented in 1D")
!   END IF

    ! Determine the dimension used for connections
    IF( PRESENT( PartitionPerm ) ) THEN
      condim = 2
      DO i=1,NBulk
        IF(PartitionPerm(i) == 0) CYCLE
        Element => Mesh % Elements(i)
        IF( Element % TYPE % ElementCode > 500 ) THEN
          condim = 3
          EXIT
        END IF
      END DO

      IF (condim==2) THEN
        CALL FindMeshEdges2D(Mesh, PartitionPerm/=0)
      ELSE
        CALL FindMeshFaces3D(Mesh, PartitionPerm/=0)
      END IF
    ELSE
      condim = dim
      CALL FindMeshEdges(Mesh)
    END IF

    CALL Info(FuncName,'Dimension for connectivity matrix: '//I2S(condim))
    
    IF( condim == 3 ) THEN
      MFacePtr => Mesh % Faces
      NFaces = Mesh % NumberOfFaces
    ELSE
      MFacePtr => Mesh % Edges
      NFaces = Mesh % NumberOfEdges
    END IF
    
    max_elfaces = 0
    DO i=1,NBulk
      Element => Mesh % Elements(i)
      IF( PRESENT( PartitionPerm) ) THEN
        IF( PartitionPerm(i) == 0 ) CYCLE
      END IF
      IF( condim == 3 ) THEN
        max_elfaces = MAX(Element % TYPE % NumberOfFaces, max_elfaces)
      ELSE        
        max_elfaces = MAX(Element % TYPE % NumberOfEdges, max_elfaces)
      END IF
    END DO
    CALL Info(FuncName,'Maximum connectivity count in graph: '&
        //I2S(max_elfaces),Level=12)

    ! Graph will be smaller if not bulk elements will be included in partitioning
    IF( PRESENT( PartitionPerm ) ) THEN
      NGraph = MAXVAL( PartitionPerm )
    ELSE
      NGraph = NBulk
    END IF
            
    CALL Info(FuncName,'Total number of rows in graph: '&
        //I2S(Ngraph),Level=12)
    
    ALLOCATE(ElemConn(max_elfaces,Ngraph), &
        ElemConnPart(max_elfaces,Ngraph), &
        NElConn(Ngraph))
    ElemConn = 0
    ElemConnPart = 0
    NElConn = 0

    ! Compute local adjacency and gather interface faces

    DO i=1,NFaces
      !Populate the local graph using non-interface faces
      Left => MFacePtr(i) % BoundaryInfo % Left
      IF(.NOT. ASSOCIATED(Left) ) CYCLE

      Right => MFacePtr(i) % BoundaryInfo % Right
      IF(.NOT. ASSOCIATED(Right)) CYCLE

      el1 = Left % ElementIndex
      el2 = Right % ElementIndex

      ! If we do not make partitioning with all elements then skip the ones
      ! that are not candidantes
      IF( PRESENT( PartitionPerm ) ) THEN
        IF( el1 > SIZE( PartitionPerm ) ) CYCLE
        IF( el2 > SIZE( PartitionPerm ) ) CYCLE
        el1 = PartitionPerm(el1)
        el2 = PartitionPerm(el2)
        IF( el1 == 0 .OR. el2 == 0 ) CYCLE
      END IF

      NElConn(el1) = NElConn(el1) + 1
      ElemConn(NElConn(el1),el1) = el2
      ElemConnPart(NElConn(el1),el1) = ParEnv % MyPE

      NElConn(el2) = NElConn(el2) + 1
      ElemConn(NElConn(el2),el2) = el1
      ElemConnPart(NElConn(el2),el2) = Parenv % MyPE
    END DO

    i = COUNT(NElConn == 0 ) 
    
    IF(i>0) CALL Warn(FuncName, 'Number of disconnected bulk element: '//I2S(i))
    
    ! Put the data into CRS format
    NtotCon = SUM(NElConn)
    CALL Info(FuncName,'Total number of connections in graph: '&
        //I2S(NtotCon),Level=10)
    WRITE(Message,'(A,F6.2)') 'Average number of connections in graph:',&
        1.0_dp * NtotCon / Ngraph
    CALL Info(FuncName,Message,Level=8)
        
    ALLOCATE(ElemAdj(NtotCon), ElemStart(Ngraph+1), ElemAdjProc(NtotCon))

    ElemStart(1) = 1
    DO i=1,Ngraph
      ElemAdj(ElemStart(i):ElemStart(i) + NElConn(i) -1) = &
           ElemConn(1:NElConn(i),i)
      ElemAdjProc(ElemStart(i):ElemStart(i) + NElConn(i) -1) = &
           ElemConnPart(1:NElConn(i),i)
      ElemStart(i+1) = ElemStart(i) + NElConn(i)
    END DO

    !PRINT *,'ElemStart:',ElemStart
    !PRINT *,'ElemAdjProc:',ElemAdjProc
    !PRINT *,'ElemAdj:',ElemAdj

    ! We only used mesh and face tables to create the dual mesh.
    ! Now release them. 
    CALL ReleaseMeshEdgeTables( Mesh )
    CALL ReleaseMeshFaceTables( Mesh )
    
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
    CHARACTER(*), PARAMETER :: FuncName="MeshParallelDualGraph"

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

    SharedNode => Mesh % ParallelInfo % GInterface
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
  !If only RmElem is provided, no nodes are removed 
  !(should this be changed? i.e. detect orphaned nodes?)
  SUBROUTINE CutMesh(Mesh, RmNode, RmElem)
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, OPTIONAL :: RmNode(:)
    LOGICAL, OPTIONAL, TARGET :: RmElem(:)
    !--------------------------------
    TYPE(Element_t), POINTER :: Element, Work_Elements(:)
    TYPE(Nodes_t), POINTER :: Nodes
    TYPE(BoundaryInfo_t), POINTER :: bInfo
    TYPE(NeighbourList_t), POINTER :: work_neighlist(:)
    REAL(KIND=dp), ALLOCATABLE :: work_xyz(:,:)
    REAL(KIND=dp), POINTER CONTIG :: work_x(:),work_y(:), work_z(:)
    INTEGER :: i,j,counter,NNodes,NBulk, NBdry,NewNNodes, NewNElems, NewNBulk,&
         NewNbdry, ElNNodes
    INTEGER, ALLOCATABLE :: Nodeno_map(:),EIdx_map(:)
    INTEGER, POINTER :: NodeIndexes(:), work_pInt(:)
    LOGICAL, POINTER :: RmElement(:),work_logical(:)
    CHARACTER(*), PARAMETER :: FuncName="CutMesh"
    NNodes = Mesh % NumberOfNodes
    NBulk = Mesh % NumberOfBulkElements
    NBdry = Mesh % NumberOfBoundaryElements

    IF(.NOT. (PRESENT(RmElem) .OR. PRESENT(RmNode))) &
         CALL Fatal(FuncName,"Need to provide at least one of RmElem, RmNode!")

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

    IF(PRESENT(RmNode)) THEN
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

      !Get rid of ParallelInfo % GInterface
      IF(ASSOCIATED(Mesh % ParallelInfo % GInterface)) THEN
        ALLOCATE(work_logical(NewNNodes))
        counter = 0
        DO i=1,NNodes
          IF(RmNode(i)) CYCLE
          counter = counter + 1
          work_logical(counter) = Mesh % ParallelInfo % GInterface(i)
        END DO
        DEALLOCATE(Mesh % ParallelInfo % GInterface)
        Mesh % ParallelInfo % GInterface => work_logical
        work_logical => NULL()
      END IF
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
        bInfo => Mesh % Elements(i) % BoundaryInfo
        IF ( ASSOCIATED(bInfo) ) THEN
          IF (ASSOCIATED(bInfo % RadiationFactors)) THEN
            IF ( ALLOCATED(bInfo % RadiationFactors % Elements ) ) THEN
              DEALLOCATE(bInfo % RadiationFactors % Factors)
              DEALLOCATE(bInfo % RadiationFactors % Elements)
            END IF
            DEALLOCATE(bInfo % RadiationFactors )
          END IF
          DEALLOCATE(bInfo)
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
  FUNCTION RedistributeMesh( Model, Mesh, ParallelMesh, FreeOldMesh, NodalVals) RESULT( NewMesh )
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh, NewMesh
    LOGICAL :: ParallelMesh, FreeOldMesh
    REAL(KIND=dp), POINTER, OPTIONAL :: NodalVals(:,:)
    
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
        CALL Fatal(FuncName,'Partition number exceeds process number: '//I2S(n))
      END IF
    END IF

    ! For convenience lets always communicate all coordinates
    ! Also, dimension might not have been set at this point for all domains
    dim = 3

    ! 0) Given the new partitioning add potential nodes to the interface
    CALL UpdateInterfaceNodeCandidates( Mesh )

    ! 1) First pack the mesh for parallel communication
    CALL PackMeshPieces(Model, Mesh, NewPart, ParallelMesh, NoPartitions, SentPack, dim, NodalVals )

    ! 2) Then sent the pieces among different partitions
    CALL CommunicateMeshPieces(Model, Mesh, ParallelMesh, NoPartitions, SentPack, RecPack)

    ! 3) Calculate element element and node counts, and creates new global2local numbering
    CALL LocalNumberingMeshPieces(Model, Mesh, NewPart, ParallelMesh, NoPartitions, RecPack, &
         GlobalToLocal, newnodes, newnbulk, newnbdry, minind, maxind)

    NewMesh => AllocateMesh( newnbulk, newnbdry, newnodes, InitParallel = .TRUE.)    

    ! 4) Finally unpack and glue the pieces on an existing mesh
    CALL UnpackMeshPieces(Model, Mesh, NewMesh, NewPart, minind, &
         maxind, RecPack, ParallelMesh, GlobalToLocal, dim, NodalVals )

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

    CALL Info(FuncName,'Deallocating temporal packed structures',Level=20)
    DO i=1,NoPartitions
      IF( ALLOCATED( SentPack(i) % idata ) ) DEALLOCATE( SentPack(i) % idata )
      IF( ALLOCATED( SentPack(i) % rdata ) ) DEALLOCATE( SentPack(i) % rdata )
      IF( ALLOCATED( RecPack(i) % idata ) ) DEALLOCATE( RecPack(i) % idata )
      IF( ALLOCATED( RecPack(i) % rdata ) ) DEALLOCATE( RecPack(i) % rdata )
    END DO

    DEALLOCATE( GlobalToLocal )

    CALL Info(FuncName,'Waiting for MPI barrier',Level=15)
    CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )


    CALL CheckTimer(FuncName,Level=5,Delete=.TRUE.)
    CALL Info(FuncName,'Distributing mesh finished',Level=8)

    CALL Finalize_Zoltan_Mesh(Newmesh)

  END FUNCTION RedistributeMesh


  SUBROUTINE Finalize_Zoltan_Mesh(Mesh)

     TYPE(Mesh_t), POINTER :: Mesh

     CHARACTER(:), ALLOCATABLE :: ElementDef, ElementDef0
     INTEGER ::  i,j,Def_DOFs(10,6)
     TYPE(Solver_t), POINTER :: Solver
     LOGICAL :: stat, GotMesh = .FALSE.


      Solver => CurrentModel % Solver


      Def_Dofs = -1; Def_Dofs(:,1)=1

      ! Define what kind of element we are working with in this solver
      !-----------------------------------------------------------------
      ElementDef = ListGetString( Solver % Values, 'Element', stat )

      IF ( .NOT. stat ) THEN
        IF ( ListGetLogical( Solver % Values, 'Discontinuous Galerkin', stat ) ) THEN
           Solver % Def_Dofs(:,:,4) = 0  ! The final value is set when calling LoadMesh2
           IF ( .NOT. GotMesh ) Def_Dofs(:,4) = MAX(Def_Dofs(:,4),0 )
           i=i+1
           Solver % DG = .TRUE.
!          CYCLE
        ELSE
           ElementDef = "n:1"
        END IF
      END IF

      ElementDef0 = ElementDef
      DO WHILE(.TRUE.)
        j = INDEX( ElementDef0, '-' )
        IF (j>0) THEN
          !
          ! Read the element definition up to the next flag which specifies the
          ! target element set
          !
          ElementDef = ElementDef0(1:j-1)
        ELSE
          ElementDef = ElementDef0
        END IF
        !  Calling GetDefs fills Def_Dofs arrays:
        CALL GetDefs( ElementDef, Solver % Def_Dofs, Def_Dofs(:,:), .NOT. GotMesh )
        IF(j>0) THEN
          ElementDef0 = ElementDef0(j+1:)
        ELSE
          EXIT
        END IF
      END DO

      CALL PrepareMesh( CurrentModel, Mesh, ParEnv % PEs>1 , Def_Dofs )

CONTAINS

!------------------------------------------------------------------------------
!> This subroutine is used to fill Def_Dofs array of the solver structure.
!> Note that this subroutine makes no attempt to figure out the index of
!> the body, so all bodies are assigned with the same element definition.
!> A similar array of reduced dimension is also filled so as to figure out
!> the maximal-complexity definition over all solvers which use the same
!> global mesh.
!------------------------------------------------------------------------------
    SUBROUTINE GetDefs(ElementDef, Solver_Def_Dofs, Def_Dofs, Def_Dofs_Update)
!------------------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(IN) :: ElementDef     !< an element definition string
      INTEGER, INTENT(OUT) :: Solver_Def_Dofs(:,:,:) !< Def_Dofs of the solver structure
      INTEGER, INTENT(INOUT) :: Def_Dofs(:,:)        !< holds the maximal-complexity definition on global mesh
      LOGICAL, INTENT(IN) :: Def_Dofs_Update         !< is .TRUE. when the definition refers to the global mesh
!------------------------------------------------------------------------------
      INTEGER, POINTER :: ind(:)
      INTEGER, TARGET :: Family(10)
      INTEGER :: i,j,l,n
      LOGICAL :: stat

      Family = [1,2,3,4,5,6,7,8,9,10]

      ! The default assumption is that the given element definition is applied 
      ! to all basic element families (note that the element sets 9 and 10 are
      ! not included since the explicit choice of the target family is 
      ! a part of the element definition string when the target index is
      ! deduced to be 9 or 10).
      !
      ind => Family(1:8)
      !
      ! If the element family is specified, change the target family 
      !
      IF (SEQL(ElementDef, 'point') )     ind => Family(1:1)
      IF (SEQL(ElementDef, 'line') )      ind => Family(2:2)
      IF (SEQL(ElementDef, 'tri') )       ind => Family(3:3)
      IF (SEQL(ElementDef, 'quad') )      ind => Family(4:4)
      IF (SEQL(ElementDef, 'tetra') )     ind => Family(5:5)
      IF (SEQL(ElementDef, 'pyramid') )   ind => Family(6:6)
      IF (SEQL(ElementDef, 'prism') )     ind => Family(7:7)
      IF (SEQL(ElementDef, 'brick') )     ind => Family(8:8)
      IF (SEQL(ElementDef, 'tri_face') )  ind => Family(9:9)
      IF (SEQL(ElementDef, 'quad_face') ) ind => Family(10:10)

      n = INDEX(ElementDef,'-')
      IF (n<=0) n=LEN_TRIM(ElementDef)
          
      j = INDEX( ElementDef(1:n), 'n:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,1) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,1) = MAX(Def_Dofs(ind,1), l)
      END IF
          
      j = INDEX( ElementDef(1:n), 'e:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,2) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,2) = MAX(Def_Dofs(ind,2), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'f:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,3) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,3) = MAX(Def_Dofs(ind,3), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'd:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l

        ! Zero value triggers discontinuous approximation within LoadMesh2,
        ! substitute the default negative initialization value to avoid troubles:
        IF (l == 0) l = -1

        Solver_Def_Dofs(ind,:,4) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,4) = MAX(Def_Dofs(ind,4), l )
      ELSE 
        IF ( ListGetLogical( Solver % Values, &
            'Discontinuous Galerkin', stat ) ) THEN
          Solver_Def_Dofs(ind,:,4) = 0
          IF ( Def_Dofs_Update ) Def_Dofs(ind,4) = MAX(Def_Dofs(ind,4),0 )
        END IF
      END IF
          
      j = INDEX( ElementDef(1:n), 'b:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,5) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,5) = MAX(Def_Dofs(ind,5), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'p:' )
      IF ( j>0 ) THEN
        IF ( ElementDef(j+2:j+2)=='%' ) THEN
          ! Seeing a p-element definition starting as p:% means that a 
          ! a special keyword construct is used so that the degree of
          ! approximation can be evaluated by calling a MATC function.
          ! This special case is handled elsewhere and we now postpone
          ! setting the right value.
          Solver_Def_Dofs(ind,:,6) = 0
        ELSE
          READ( ElementDef(j+2:), * ) l
          Solver_Def_Dofs(ind,:,6) = l
          IF ( Def_Dofs_Update ) Def_Dofs(ind,6) = MAX(Def_Dofs(ind,6), l )
         END IF
      END IF

!------------------------------------------------------------------------------
    END SUBROUTINE GetDefs
!------------------------------------------------------------------------------

    END SUBROUTINE Finalize_Zoltan_Mesh

  
  FUNCTION ElementPartitions( Mesh, ElemInd, NewPart, IndPart ) RESULT ( npart )
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: ElemInd
    INTEGER, POINTER :: NewPart(:)
    INTEGER, ALLOCATABLE :: IndPart(:)
    INTEGER :: nPart
    
    INTEGER :: n
   
    IndPart(1) = NewPart(ElemInd)
    npart = 1
    
    IF(ASSOCIATED( Mesh % Halo ) ) THEN
      IF( ASSOCIATED( Mesh % Halo(ElemInd) % Neighbours ) ) THEN
        n = SIZE( Mesh % Halo(ElemInd) % Neighbours )
        IndPart(2:n+1) = Mesh % Halo(ElemInd) % Neighbours(1:n)
        npart = npart + n
        ! PRINT *,'halo:',npart,IndPart(1:n+1)       
      END IF
    END IF
       
  END FUNCTION ElementPartitions
  
  
  !> Converts element datastructure into a integer and real stream to facilitate
  !> sending to another partition.
  !------------------------------------------------------------------------------
  SUBROUTINE PackMeshPieces(Model, Mesh, NewPart, ParallelMesh, NoPartitions, &
      SentPack, dim, NodalVals )

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: ParallelMesh
    INTEGER, POINTER :: NewPart(:)
    INTEGER :: NoPartitions, dim
    REAL(KIND=dp), POINTER, OPTIONAL :: NodalVals(:,:)
    !------------------------
    
    TYPE( MeshPack_t), ALLOCATABLE, TARGET :: SentPack(:)
    TYPE(Element_t), POINTER :: Element, Parent
    INTEGER :: i,j,k,l,n,nblk,nbdry,allocstat,part,elemcode,geom_id,sweep
    LOGICAL :: CheckNeighbours, IsBulk
    TYPE(NeighbourList_t),POINTER  :: NeighbourList(:)
    TYPE(MeshPack_t), POINTER :: PPack
    INTEGER, POINTER :: TmpPart(:)
    LOGICAL :: HaveHalo
    INTEGER :: NPart, nVals
    INTEGER, ALLOCATABLE, SAVE :: IndPart(:)
    CHARACTER(*), PARAMETER :: FuncName='PackMeshPieces'

    CALL Info(FuncName,'Packing mesh pieces for sending',Level=8)

    ! Allocate and initialize the structures used to communicate this mesh
    n = NoPartitions
    ALLOCATE( SentPack( n ) )
    IF(.NOT. ALLOCATED(IndPart)) ALLOCATE( IndPart(20))
    
    SentPack(1:n) % NumberOfNodes = 0
    SentPack(1:n) % NumberOfBulkElements = 0
    SentPack(1:n) % NumberOfBoundaryElements = 0
    SentPack(1:n) % icount = 0
    SentPack(1:n) % rcount = 0
    SentPack(1:n) % lcount = 0
    SentPack(1:n) % indpos = 0
    SentPack(1:n) % bcpos = 0

    IF( Mesh % NumberOfNodes == 0 ) THEN
      CALL Info(FuncName,'Mesh is empty, nothing to pack',Level=10)
      RETURN
    END IF

    HaveHalo = ASSOCIATED( Mesh % Halo )
    IF( HaveHalo ) THEN
      CALL info(FuncName,'Including halo elements in communication',Level=10)        
    END IF
            
    nblk = Mesh % NumberOfBulkElements
    nbdry = Mesh % NumberOfBoundaryElements

    CALL Info(FuncName,'Packing '//I2S(nblk)//'+'//I2S(nbdry)//' elements for sending!',Level=10)
    nblk = Mesh % NumberOfBulkElements
    nbdry = Mesh % NumberOfBoundaryElements

    nVals = 0
    IF(PRESENT(NodalVals)) THEN
      IF(ASSOCIATED(NodalVals)) nVals = SIZE(NodalVals,2)
    END IF
    CALL Info(FuncName,'Packing '//I2S(dim)//'D nodes and '//I2S(nVals)//' nodal fields!',Level=10)
       
    IF( SIZE( NewPart ) < nblk + nbdry ) THEN
      CALL Info(FuncName,'Growing the partition vector to accounts BCs',Level=8)
      ALLOCATE( TmpPart( nblk + nbdry ) )
      TmpPart(1:nblk) = NewPart(1:nblk)
      TmpPart(nblk+1:nblk+nbdry) = 0
      DEALLOCATE( NewPart )
      NewPart => TmpPart
      Mesh % RePartition => TmpPart
      NULLIFY( TmpPart )
    END IF

    IF( ANY(NewPart(nblk+1:nblk+nbdry) <= 0 ) ) THEN
      CALL Info(FuncName,'Inheriting partition vector for BCs',Level=8)
      DO i=nblk+1,nblk+nbdry
        IF( NewPart(i) > 0 ) CYCLE
        Element => Mesh % Elements(i)
        IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) THEN
          CALL Fatal(FuncName,'Boundary element lacking boundary info')
        END IF
        Parent => Element % BoundaryInfo % Left
        IF( .NOT. ASSOCIATED(Parent) ) THEN
          Parent => Element % BoundaryInfo % Right
        END IF
        IF( .NOT. ASSOCIATED( Parent ) ) THEN
          CALL Fatal(FuncName,'Boundary element lacking parent info')
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
        CALL info(FuncName,'Using existing information of shared nodes',Level=10)
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

        IsBulk = ( i <= nblk )
        Element => Mesh % Elements(i)
        elemcode = Element % Type % ElementCode
        n = Element % TYPE % NumberOfNodes

        Npart = ElementPartitions( Mesh, i, NewPart, IndPart )        
        DO l=1,NPart
          part = IndPart(l)
          
          IF( part-1 == ParEnv % MyPe ) CYCLE
          PPack => SentPack(part)
        
          IF( Sweep == 1) THEN
            ! Add partition index to communication if there is halo
            IF( HaveHalo ) THEN
              PPack % icount = PPack % icount + 1
            END IF

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

            ! Owner partition only needed for halo
            IF( HaveHalo ) THEN
              PPack % idata(PPack % icount+1) = NewPart(i)
              PPack % icount = PPack % icount + 1
            END IF

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
      END DO


      IF( Sweep == 1 ) THEN
        ! Set the offset for the nodal data. It is needed in unpackig.
        SentPack(1:NoPartitions) % indpos = SentPack(1:NoPartitions) % icount

        ! We need to allocate a logical mask to mark the nodes to sent to given partition
        ! Note that we only want to allocate the flag for partitions that also receive
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
          Element => Mesh % Elements(i)
          elemcode = Element % TYPE % ElementCode
          n = Element % TYPE % NumberOfNodes
          
          Npart = ElementPartitions( Mesh, i, NewPart, IndPart )        
          DO l=1,NPart
            part = IndPart(l)
            
            ! No need to sent to self
            IF( part-1 == ParEnv % MyPe ) CYCLE

            PPack => SentPack(part)

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
        END DO

        ! Add the nodes count to be sent to the data structure
        DO part=1,NoPartitions
          ! Nothing to sent for self
          IF( part-1 == ParEnv % MyPe ) CYCLE
          
          PPack => SentPack(part)

          IF(  PPack % icount <= 5 ) CYCLE

          PPack % NumberOfNodes = COUNT( PPack % NodeMask )
          PPack % icount = PPack % icount + PPack % NumberOfNodes
          PPack % rcount = dim * PPack % NumberOfNodes
          IF(nVals > 0 ) THEN
            PPack % rcount = PPack % rcount + nVals * PPack % NumberOfNodes
          END IF
          PPack % lcount = PPack % NumberOfNodes
          
          ALLOCATE( PPack % idata(PPack % icount), &
              PPack % rdata(PPack % rcount), &
              PPack % ldata(PPack % lcount ), & 
              STAT = allocstat )
          IF( allocstat /= 0 ) THEN
            CALL Fatal(FuncName,'Could not allocate vectors for data')
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

              ! Tentatively add some nodal variables
              IF(nVals > 0) THEN
                PPack % rdata(PPack % rcount+1:PPack % rcount+nVals) = NodalVals(i,1:nVals)
                PPack % rcount = PPack % rcount + nVals
              END IF
                
              PPack % lcount = PPack % lcount + 1
              PPack % ldata(PPack % lcount) = Mesh % ParallelInfo % GInterface(i)
            END IF
          END DO

        END DO
      END IF
    END DO

    CALL Info(FuncName,'Finished packing mesh pieces',Level=8)

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
    LOGICAL :: HaveHalo
    CHARACTER(*), PARAMETER :: FuncName='CommunicateMeshPieces'

    
    CALL Info(FuncName,'communicating mesh pieces in parallel',Level=6)

    ni = SUM( SentPack(1:NoPartitions) % icount )
    CALL Info(FuncName,'Number of integer values to sent: '//I2S(ni),Level=8)
    nr = SUM( SentPack(1:NoPartitions) % rcount )
    CALL Info(FuncName,'Number of real values to sent: '//I2S(nr),Level=8)
    nl = SUM( SentPack(1:NoPartitions) % lcount )
    CALL Info(FuncName,'Number of logical values to sent: '//I2S(nl),Level=8)

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


    ! Should we communicate element owner index 
    HaveHalo = ( ASSOCIATED( Mesh % Halo ) )
    CALL MPI_ALLREDUCE(HaveHalo, Mesh % HaveHalo, 1, MPI_LOGICAL, &
        MPI_LOR, ELMER_COMM_WORLD, ierr )
    IF( Mesh % HaveHalo ) THEN
      CALL Info(FuncName,'Assuming halo being communicated',Level=12)
    END IF
       
    
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
    
    ! Receive data sizes:
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

    CALL Info(FuncName,'Waiting for the 1st barrier',Level=15)
    CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )

    IF( InfoActive(15) ) THEN
      DO i = 1, NoPartitions
        IF( i-1 == ParEnv % MyPe ) CYCLE
        IF( SentPack(i) % icount > 5 ) THEN
          PRINT *,'Mesh send sizes: '//I2S(ParEnv % Mype)//'-'//I2S(i-1), &
              SentPack(i) % icount, SentPack(i) % rcount, SentPack(i) % lcount
        END IF
        IF( RecPack(i) % icount > 5 ) THEN       
          PRINT *,'Mesh recv sizes: '//I2S(ParEnv % Mype)//'-'//I2S(i-1), &
              RecPack(i) % icount, RecPack(i) % rcount, RecPack(i) % lcount      
        END IF
      END DO
    END IF
        
    n = SUM( RecPack(1:NoPartitions) % icount )
    CALL Info('PackDataToSend','Number of integer values to receive: '//I2S(n),Level=8)
    n = SUM( RecPack(1:NoPartitions) % rcount )
    CALL Info('PackDataToSend','Number of real values to receive: '//I2S(n),Level=8)
    n = SUM( RecPack(1:NoPartitions) % lcount )
    CALL Info('PackDataToSend','Number of logical values to receive: '//I2S(n),Level=8)

    ! Allocate data sizes for receiving data
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
    CALL Info(FuncName,'Now sending the actual data',Level=12)
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

    ! Receive data:
    !--------------------------
    CALL Info(FuncName,'Now receiving the actual integer data',Level=12)
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

    CALL Info(FuncName,'Waiting for the 2nd barrier',Level=15)
    CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )

    CALL Info(FuncName,'Finished communicating mesh pieces',Level=8)

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
    INTEGER :: NPart
    INTEGER, ALLOCATABLE, SAVE :: IndPart(:)
    LOGICAL :: HaveHalo
    
    
    CALL Info('LocalNumberingMeshPieces','Renumbering local nodes in each partition',Level=8)

    HaveHalo = Mesh % HaveHalo
    IF(.NOT. ALLOCATED(IndPart)) ALLOCATE( IndPart(20))
    
    newnbulk = 0
    newnbdry = 0

    ! Compute the number of elements staying
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Npart = ElementPartitions( Mesh, i, NewPart, IndPart )
      IF( ANY( IndPart(1:NPart) == ParEnv % MyPe + 1 ) ) THEN
        IF( i <= Mesh % NumberOfBulkElements ) THEN
          newnbulk = newnbulk + 1
        ELSE
          newnbdry = newnbdry + 1
        END IF
      END IF
    END DO
    CALL Info('LocalNumberingMeshPieces','Number of staying elements: '&
         //I2S(newnbulk+newnbdry),Level=15)

    ! Add the number of elements coming from different partitions
    DO part=1,NoPartitions
      IF( part-1 == ParEnv % Mype ) CYCLE
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

    CALL Info('LocalNumberingMeshPieces','Number of combined elements: '&
         //I2S(newnbulk+newnbdry),Level=8)


    ! Find the range of initial global indeces
    ! This is conservative since it includes all the initial global indexes
    n = Mesh % NumberOfNodes
    IF( n > 0 ) THEN
      IF( ParallelMesh ) THEN
        IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % GlobalDofs ) ) THEN
          CALL Fatal('LocalNumberingMeshPieces','ParallelMesh assumed but no GlobalDofs associated!')          
        END IF
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
      IF( part-1 == ParEnv % Mype ) CYCLE
      PPack => RecPack(part)
      IF( PPack % icount <= 5 ) CYCLE
      i1 = PPack % indpos + 1
      i2 = PPack % icount
      minind = MIN( minind, MINVAL( PPack % idata(i1:i2) ) )
      maxind = MAX( maxind, MAXVAL( PPack % idata(i1:i2) ) )
    END DO

    IF(minind == 0) RETURN

    !CALL Info('LocalNumberingMeshPieces','Global index range '&
    !     //I2S(minind)//' to '//I2S(maxind),Level=12)

    ! Allocate the vector for local renumbering
    ALLOCATE( GlobalToLocal(minind:maxind), STAT=allocstat)
    IF( allocstat /= 0 ) THEN
      CALL Fatal('LocalNumberingMeshPieces','Could not allocate vectors for indexes')
    END IF
    GlobalToLocal = 0

    ! Check which of the staying nodes are used
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      
      Npart = ElementPartitions( Mesh, i, NewPart, IndPart )
      IF( ANY( IndPart(1:NPart) == ParEnv % MyPe + 1 ) ) THEN
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

        ! If we have halo the 1st index is the owner partition
        IF( HaveHalo ) icount = icount+1 
        
        elemcode = PPack % idata(icount+2)

        Etype => GetElementType( elemcode )
        n = Etype % NumberOfNodes
        icount = icount + 3

        IF( icount + n > SIZE( PPack % idata) ) THEN
          CALL Fatal('LocalNumberingMeshPieces','icount out of range')
        END IF

        DO j=1,n
          gind = PPack % idata(icount+j)

          IF( gind < minind .OR. gind > maxind ) THEN
            CALL Fatal('LocalNumberingMeshPieces','gind out of bounds')
          END IF

          GlobalToLocal(gind) = 1
        END DO
        icount = icount + n
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
    
    CALL Info('LocalNumberingMeshPieces','Combined number of nodes: '//I2S(newnodes),Level=8)

  END SUBROUTINE LocalNumberingMeshPieces



  !> Converts element data structure from integer and real streams to FE meshes.
  !> The idea is that the elements and nodes are appended on top of an existing
  !> mesh structure such that there could be elements already at the bottom.
  !------------------------------------------------------------------------------
  SUBROUTINE UnpackMeshPieces(Model, Mesh, NewMesh, NewPart, &
      minind, maxind, RecPack, ParallelMesh, GlobalToLocal, &
      dim, NodalVals )

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh, NewMesh
    LOGICAL :: ParallelMesh
    INTEGER, POINTER :: NewPart(:)
    INTEGER, ALLOCATABLE :: GlobalToLocal(:)
    INTEGER :: dim,minind,maxind
    REAL(KIND=dp), POINTER, OPTIONAL :: NodalVals(:,:)
    !----------------------------------
    TYPE( MeshPack_t), ALLOCATABLE, TARGET :: RecPack(:)
    TYPE(Element_t), POINTER :: Element, Element0
    INTEGER :: i,j,k,n,t,nbulk,nbdry,allocstat,part,elemcode,elemindex,geom_id,sweep,partindex
    INTEGER :: gind,lind,rcount,icount,lcount,minelem,maxelem,newnbdry,newnodes,newnbulk
    LOGICAL :: IsBulk, Found, HaveParent
    TYPE(MeshPack_t), POINTER :: PPack
    INTEGER, ALLOCATABLE :: GlobalToLocalElem(:), LeftParent(:), RightParent(:)
    INTEGER :: NPart, errcount, NVals 
    INTEGER, ALLOCATABLE, SAVE :: IndPart(:)
    LOGICAL :: HaveHalo
    REAL(KIND=dp), POINTER :: NewVals(:,:)
    CHARACTER(*), PARAMETER :: Caller = 'UnpackMeshPieces'
    
    CALL Info(Caller,'Unpacking mesh pieces to form a new mesh',Level=12)

    minelem = HUGE( minelem )
    maxelem = 0

    newnodes = NewMesh % NumberOfNodes
    newnbulk = NewMesh % NumberOfBulkElements 
    newnbdry = NewMesh % NumberOfBoundaryElements 

    HaveHalo = Mesh % HaveHalo
    IF(.NOT. ALLOCATED(IndPart)) ALLOCATE( IndPart(20))    
    
    nbulk = 0
    nbdry = 0
    errcount = 0
    nVals = 0
    IF(PRESENT(NodalVals)) THEN
      ! Make the nodal vals to be of right size.
      IF(ASSOCIATED(NodalVals)) nVals = SIZE(NodalVals,2)
      IF(nVals > 0) THEN
        ALLOCATE(NewVals(newnodes,nVals))
        NewVals = 0.0_dp
      END IF
    END IF
    CALL Info(Caller,'Unpacking '//I2S(dim)//'D nodes and '//I2S(nVals)//' nodal fields!',Level=10)
    
    ! There are temporal arrays needed to inherit the parent information
    ALLOCATE( LeftParent(NewNBulk+1:NewNBulk+NewNBdry ) ) 
    LeftParent = 0
    ALLOCATE( RightParent(NewNBulk+1:NewNBulk+NewNBdry ) ) 
    RightParent = 0
       
    
    CALL Info(Caller,'Copying staying elements',Level=20)
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      Npart = ElementPartitions( Mesh, i, NewPart, IndPart )
      IF( ALL( IndPart(1:NPart) /= ParEnv % MyPe + 1 ) ) CYCLE

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
      Element % NDOFs = n

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
      IF( HaveHalo ) THEN
        Element % PartIndex = Mesh % RePartition(i)-1 
      ELSE
        Element % PartIndex = ParEnv % MyPe
      END IF
        
      NULLIFY( Element % NodeIndexes )
      ALLOCATE( Element % NodeIndexes(n), STAT = allocstat )
      IF( allocstat /= 0 ) THEN
        CALL Fatal(Caller,'Cannot allocate '//I2S(n)//' node indexes?')
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

    IF( .NOT. ASSOCIATED( NewMesh % ParallelInfo % GInterface ) ) THEN
      ALLOCATE( NewMesh % ParallelInfo % GInterface( NewMesh % NumberOfNodes ), STAT = allocstat )
      IF( allocstat /= 0 ) THEN
        CALL Fatal(Caller,'Cannot allocate partition interface?')
      END IF
      NewMesh % ParallelInfo % GInterface = .FALSE.
    END IF

    
    DO i=1,Mesh % NumberOfNodes
      j = i
      IF( ParallelMesh ) j = Mesh % ParallelInfo % GlobalDofs(i)
      IF( j < minind .OR. j > maxind ) THEN
        errcount = errcount + 1
        PRINT *,'Global index out of bounds:',ParEnv % Mype, j,minind,maxind
      END IF
      k = GlobalToLocal(j)
      IF( k == 0 ) CYCLE
      IF( k > newnodes ) THEN
        errcount = errcount + 1
        PRINT *,'Mapped index out of bounds:',ParEnv % Mype,k,newnodes
      END IF
      NewMesh % Nodes % x(k) = Mesh % Nodes % x(i)
      NewMesh % Nodes % y(k) = Mesh % Nodes % y(i)
      IF( dim == 3 ) NewMesh % Nodes % z(k) = Mesh % Nodes % z(i)

      IF(nVals > 0 ) THEN 
        NewVals(k,1:nVals) = NodalVals(i,1:nVals)
      END IF
        
      NewMesh % ParallelInfo % GInterface(k) = Mesh % ParallelInfo % GInterface(i)
    END DO


    CALL Info(Caller,'Unpacking incoming elements',Level=20)
    DO part=1,ParEnv % PEs
      IF( part-1 == ParEnv % MyPe ) CYCLE

      PPack => RecPack(part)
      IF( PPack % icount <= 5 ) CYCLE

      CALL Info(Caller,'Unpacking piece '//I2S(part)//' with '&
           //I2S(PPack % NumberOfBulkElements + PPack % NumberOfBoundaryElements)//&
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

        IF( HaveHalo ) THEN
          partindex = PPack % idata(icount+1)
          icount = icount + 1
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
          CALL Fatal(Caller,'Could not get element code: '//I2S(elemcode))
        END IF

        n = Element % Type % NumberOfNodes
        Element % NDOFs = n

        Element % GElementIndex = elemindex
        Element % ElementIndex = t

        IF( IsBulk ) THEN
          minelem = MIN( minelem, Element % GElementIndex )
          maxelem = MAX( maxelem, Element % GElementIndex )
        END IF
        
        ! Change the owner partition of the element
        IF( HaveHalo ) THEN
          Element % PartIndex = partindex-1
        ELSE
          Element % PartIndex = ParEnv % MyPe
        END IF
        
        IF( IsBulk ) THEN
          Element % BodyId = geom_id
          icount = icount + 3
        ELSE
          IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) THEN
            ALLOCATE( Element % BoundaryInfo, STAT = allocstat )
            IF( allocstat /= 0 ) THEN
              CALL Fatal(Caller,'Could not allocate boundary info!')
            END IF
          END IF

          !Check geom id is legal
          IF(geom_id < 0 .OR. geom_id > Model % NumberOfBCs) THEN
            PRINT *, ParEnv % MyPE,'Received element ',i,' from part: ',&
                 part,' with constraint: ',geom_id
            CALL Fatal(Caller, "Unexpected constraint on BC element")
          END IF

          Element % BoundaryInfo % Constraint = geom_id
          IF( geom_id > 0 ) THEN
            Element % BodyId  = ListGetInteger( &
                Model % BCs(geom_id) % Values, 'Body Id', Found, 1, Model % NumberOfBodies )
          END IF
            
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
            errcount = errcount + 1
            PRINT *,'Node index out of bounds:',ParEnv % Mype,k,minind,maxind
          END IF
          IF( GlobalToLocal(k) <= 0 .OR. GlobalToLocal(k) > newnodes ) THEN
            errcount = errcount + 1
            PRINT *,'Local index out of bounds:',ParEnv % Mype,Element % PartIndex, k,minind,maxind,&
                SIZE(GlobalToLocal), GlobalToLocal(k)
          END IF
          Element % NodeIndexes(j) = GlobalToLocal(k)
        END DO

        ! Advance the counter for the data
        icount = icount + n
      END DO

      CALL Info(Caller,'Finished unpacking piece',Level=20)
      IF( icount /= PPack % indpos ) THEN
        CALL Fatal(Caller,'Inconsistent icount value: '//I2S(icount))
      END IF
    END DO

    IF( errcount > 0 ) THEN
      CALL Fatal(Caller,'Encountered '//I2S(errcount)//' indexing issues in elements')
    END IF
      
    IF( minelem <= maxelem ) THEN
      ! First create global to local array for the elements 
      CALL Info(Caller,'Global element index range: '&
          //I2S(minelem)//' to '//I2S(maxelem),Level=8)
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

      i = COUNT( GlobalToLocalElem == 0 )
      j = SIZE(GlobalToLocal)
      CALL Info(Caller,'Number of mapping indexes defined '//I2S(i)//&
          ' out of '//I2S(j),Level=10)
      
      ! Then use the temporal vectors to repoint the left and right indexes to elements
      DO i = newnbulk+1, newnbulk + newnbdry
        Element => NewMesh % Elements(i)        
        HaveParent = .FALSE.

        j = LeftParent(i)
        IF( j >= minelem .AND. j <= maxelem ) THEN
          k = GlobalToLocalElem(j)
          IF( k == 0 ) THEN
            CONTINUE
          ELSE IF( k < 0 .OR. k > newnbulk ) THEN
            errcount = errcount + 1
            PRINT *,'k left out of bounds:',ParEnv % MyPe, i, j, k, newnbulk, minelem, maxelem
          ELSE
            HaveParent = .TRUE.
            Element % BoundaryInfo % Left => NewMesh % Elements(k)
          END IF
        END IF
        j = RightParent(i)
        IF( j >= minelem .AND. j <= maxelem ) THEN
          k = GlobalToLocalElem(j)
          IF( k == 0 ) THEN
            CONTINUE
          ELSE IF( k <= 0 .OR. k > newnbulk ) THEN
            errcount = errcount + 1
            PRINT *,'k right out of bounds:',ParEnv % MyPe, i, j, k, newnbulk, minelem, maxelem
          ELSE
            HaveParent = .TRUE.
            Element % BoundaryInfo % Right => NewMesh % Elements(k)
          END IF
        END IF

        IF(.NOT. HaveParent ) THEN
          errcount = errcount + 1
          PRINT *,'No parent for boundary element:',ParEnv % MyPe, i, newnbulk, minelem, maxelem
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
        
        IF(nVals > 0 ) THEN 
          NewVals(k,1:nVals) = PPack % rdata(rcount+1:rcount+nVals)
          rcount = rcount + nVals
        END IF
        
        NewMesh % ParallelInfo % GInterface(k) = PPack % ldata(lcount+1)
        lcount = lcount + 1
      END DO
    END DO

    IF( errcount > 0 ) THEN
      CALL Fatal(Caller,'Encountered '//I2S(errcount)//' indexing issues in nodes')
    END IF
     
    n = COUNT( NewMesh % ParallelInfo % GInterface )
    CALL Info(Caller,'Potential interface nodes '//I2S(n)//' out of '&
        //I2S(NewMesh % NumberOfNodes),Level=20)
   
    CALL Info(Caller,'Creating local to global numbering for '&
         //I2S(newnodes)//' nodes',Level=20)
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

     IF(nVals > 0 ) THEN 
       DEALLOCATE(NodalVals)
       NodalVals => NewVals
     END IF

     
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

    ! The interface includes the old interface between partitions
    IF( .NOT. ASSOCIATED( Mesh % ParallelInfo % GInterface ) ) THEN
      ALLOCATE( Mesh % ParallelInfo % GInterface( n ), STAT=allocstat )
      IF( allocstat /= 0 ) THEN
        CALL Fatal(Caller,'Allocation error for parallel interface!')
      END IF
      Mesh % ParallelInfo % GInterface = .FALSE.
    END IF
    PartInterface => Mesh % ParallelInfo % GInterface

    IF( .NOT. ASSOCIATED( Mesh % RePartition ) ) THEN
      CALL Fatal(Caller,'Allocation error for parallel interface!')
    END IF
    ElementPart => Mesh % RePartition
    
    ALLOCATE( PrevPartition( n ), STAT=allocstat ) 
    IF( allocstat /= 0 ) THEN
      CALL Fatal(Caller,'Allocation error for prev partition!')
    END IF
    PrevPartition = 0

    ! And all nodes which have nodes on the interface between two partitions
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
      
      IF( ASSOCIATED( Mesh % Halo ) ) THEN
        IF( ASSOCIATED( Mesh % Halo(i) % Neighbours ) ) THEN
          PartInterface( Element % NodeIndexes ) = .TRUE.
        END IF
      END IF
    END DO

    n = COUNT( PartInterface ) 
    DEALLOCATE( PrevPartition )
    
    CALL Info(Caller,'Number of potential nodes at the interface: '//I2S(n),Level=10)      
      
  END SUBROUTINE UpdateInterfaceNodeCandidates

  !> Based on a conservative list of potential interface nodes
  !> in ParallelInfo % GInterface, find real interface nodes &
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

    Iface => Mesh % ParallelInfo % GInterface
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


  !> Makes a serial mesh partitiong. Currently uses geometric criteria or Zoltan.
  !> Includes some hybrid strategies where the different physical domains (bc & bulk)
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
     INTEGER, ALLOCATABLE :: ParameterInd(:), ElementSet(:), PartCount(:)
     INTEGER, POINTER :: ElementPart(:)
     INTEGER :: NumberOfSets, NumberOfBoundarySets, SetNo, id, NoEqs
     LOGICAL, POINTER :: PartitionCand(:)
     INTEGER :: i,j,j0,j1,k,n,m,allocstat
     LOGICAL :: Found, PartBalance, EqInterface, MasterHalo
     INTEGER, ALLOCATABLE :: EquationPart(:)    
     TYPE(NeighbourList_t),POINTER  :: NeighbourList(:)
     INTEGER :: PartOffset, PartOffsetBC
     LOGICAL, ALLOCATABLE :: MasterElement(:)
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
       CALL Info(FuncName,'Allocating RePartition table of size: '//I2S(n),Level=20)
       ALLOCATE( Mesh % RePartition( n ), STAT = allocstat)
       IF( allocstat /= 0 ) THEN
         CALL Fatal(FuncName,'Allocation error for repartitioning vector')       
       END IF
     END IF
     ElementPart => Mesh % RePartition 
     
     PartitionCand = .FALSE.
     ElementSet = 0
     ElementPart = 0
     PartOffset = 0
     PartOffsetBC = 0

     ! Only use the halo for master elements
     MasterHalo = ListGetLogical( Model % Simulation,'Boundary Partition Halo Master',Found)
     IF( MasterHalo ) THEN
       ALLOCATE( MasterElement(n) )
       MasterElement = .FALSE.
     END IF
          
     n = MAX( Model % NumberOfBCs, Model % NumberOfEquations )
     ALLOCATE( ParameterInd(n), STAT = allocstat ) 
     IF( allocstat /= 0 ) THEN
       CALL Fatal(FuncName,'Allocation error for ParameterInd')
     END IF

     PartBalance = ListGetLogical( Model % Simulation,'Partition Equation Balance',Found )      
     ParameterInd = 0

     EqInterface = ListGetLogical( Model % Simulation,'Partition Equation Interface',Found )
                  

     CALL Info(FuncName,'Partitioning boundary elements sets') 
     CALL InitializeBoundaryElementSet(NumberOfBoundarySets)
     
     IF( NumberOfBoundarySets > 0 ) THEN
       IF( EqInterface ) THEN
         CALL Fatal(FuncName,'Cannot deal with boundary sets and interface partition at same time!')
       END IF

       DO SetNo = 1, NumberOfBoundarySets
         CALL Info(FuncName,'Doing boundary partitioning set: '//I2S(SetNo),Level=12)

         SectionParams => NULL()
         ! Get the bc-specific partitioning commands, if any
         IF( ParameterInd(SetNo) > 0 ) THEN
           id = ParameterInd(SetNo)
           IF( id <= Model % NumberOfBCs ) THEN
             SectionParams => Model % BCs(id) % Values
           END IF
         END IF
         CALL PartitionMeshPart(SetNo,SectionParams,.TRUE.,PartOffset,PartOffsetBC)
       END DO

       IF( ListGetLogical( Params,'Partition Mesh Merge Boundaries',Found ) ) THEN
         CALL MergeJoinedPartitions(.TRUE.)
         PartOffset = MAXVAL( ElementPart )
         PartOffsetBC = PartOffset
       END IF
       
       CALL InheritBoundaryToBulkPart()

       IF( ListGetLogical( Params,'Boundary Partition Halo',Found ) ) THEN
         CALL GenerateSetHalo()
       END IF

       IF( ListGetLogical( Params,'Boundary Bounding Box Halo',Found ) ) THEN
         CALL GenerateBBoxHalo()
         CALL InheritHaloToBulkPart
       END IF
       
       CALL ExtendBoundaryPart(.TRUE.)

       ! We could do the halo stuff also after extending 
       IF( ListGetLogical( Params,'Boundary Bounding Box Halo',Found ) ) THEN
         ! CALL GenerateBBoxHalo()
       END IF
       IF( ListGetLogical( Params,'Boundary Partition Halo',Found ) ) THEN
         !CALL GenerateSetHalo()
       END IF
       
     ELSE IF( EqInterface ) THEN
       CALL Info(FuncName,'Doing equation interface set',Level=12)

       NumberOfBoundarySets = 1       
       CALL SetInterfacePartition()
       CALL ExtendBoundaryPart(.FALSE.)
       CALL InheritBulkToBoundaryPart()
       PartOffset = 1
       PartOffsetBC = 1
     END IF
     
     CALL Info(FuncName,'Partition the bulk elements sets')

     
     j0 = MAXVAL( ElementPart )
     ParameterInd = 0
     
     CALL Info(FuncName,'Maximum partition index for BCs: '//I2S(j0),Level=8)
     
     CALL InitializeBulkElementSet(NumberOfSets)

     
     DO SetNo = 1, NumberOfSets
       CALL Info(FuncName,'Doing bulk partitioning set: '//I2S(SetNo),Level=8)

       SectionParams => NULL()
       
       ! Get the equation-specific parameters, if any
       IF( ParameterInd(SetNo) > 0 ) THEN
         id = ParameterInd(SetNo)
         IF( id <= Model % NumberOfEquations ) THEN
           SectionParams => Model % Equations(id) % Values
         END IF
       END IF

       CALL PartitionMeshPart(SetNo+NumberOfBoundarySets,SectionParams,.FALSE.,PartOffset,PartOffsetBC)        
       IF( SetNo == 1 ) THEN
         j0 = PartOffsetBC
         j1 = PartOffset
       END IF
     END DO

      IF( PartBalance ) THEN
       NoEqs = Model % NumberOfEquations
       CALL Info(FuncName,'Removing offset for optimal load balancing of '//I2S(NoEqs)//' equations',Level=10)
       DO i=1,SIZE(ElementPart)
         j = ElementPart(i)
         IF(j<=j1) CYCLE
         ElementPart(i) = MODULO(j-j0-1,j1-j0)+j0+1
       END DO
     END IF    
     
     CALL InheritBulkToBoundaryPart()
     
     IF( ListGetLogicalAnySolver( Model,'Discontinuous Galerkin') ) THEN
       CALL GenerateDGHalo()
     END IF

     CALL CreateNeighbourList()


     IF( InfoActive(5) ) THEN
       n = Mesh % NumberOfBulkElements
       m = MAXVAL(ElementPart(1:n))
              
       ALLOCATE(PartCount(0:m))
       PartCount = 0
       DO i=1,n
         j = ElementPart(i)
         PartCount(j) = PartCount(j) + 1
       END DO       
       DO i=0,m
         PRINT *,'Elements in Partition:',i,PartCount(i)
       END DO
     END IF

     
100  CALL Info(FuncName,'All done for now',Level=12)

     
   CONTAINS

     ! Inherit partition from a boundary partition.
     ! In case of conflict the 1st occurrence prevails.
     !-----------------------------------------------------
     SUBROUTINE InheritBoundaryToBulkPart()
       
       TYPE(Element_t), POINTER :: Element, Parent
       INTEGER :: t, LeftRight, BoundPart, NoHerited, NoConflict, ElemIndx

       CALL Info(FuncName,'Inheriting the boundary partitioning into the bulk mesh') 

       NoHerited = 0
       NoConflict = 0

       DO t=Mesh % NumberOfBulkElements + 1,&
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 
         Element => Mesh % Elements(t) 

         BoundPart = ElementPart( t )
         ! Don't inherit from unset elements
         IF( BoundPart == 0 ) CYCLE
         
         IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
           DO LeftRight=0,1
             IF( LeftRight == 0 ) THEN
               Parent => Element % BoundaryInfo % Left
             ELSE
               Parent => Element % BoundaryInfo % Right
             END IF

             IF( ASSOCIATED( Parent ) ) THEN               
               ElemIndx = Parent % ElementIndex

               IF( ElementPart( ElemIndx ) == 0 ) THEN
                 NoHerited = NoHerited + 1
                 ElementPart( ElemIndx ) = BoundPart                 
               ELSE IF( ElementPart( ElemIndx ) /= BoundPart ) THEN
                 NoConflict = NoConflict + 1
               END IF

               ! Inherit also the master status.
               ! This is done always since we need the halo if the bulk
               ! element is associated to any master element.
               IF( MasterHalo ) MasterElement( ElemIndx ) = MasterElement( t )               
             END IF
           END DO
         END IF
       END DO
       
       CALL Info(FuncName,'Number of herited bulk elements: '//I2S(NoHerited))
       CALL Info(FuncName,'Number of conflicted bulk elements: '//I2S(NoConflict))
       
     END SUBROUTINE InheritBoundaryToBulkPart


     ! Inherit partition from a boundary partition.
     ! In case of conflict the 1st occurrence prevails.
     !-----------------------------------------------------
     SUBROUTINE InheritHaloToBulkPart()
       
       TYPE(Element_t), POINTER :: Element, Parent
       INTEGER :: t, tB, LeftRight, BoundPart, NoHerited, NoConflict, ElemIndx
       INTEGER :: npart, npartB
       INTEGER, ALLOCATABLE, SAVE :: IndPart(:), IndPartB(:)
       
       IF( .NOT. ASSOCIATED( Mesh % Halo ) ) RETURN
       
       CALL Info(FuncName,'Inheriting the boundary halo into the bulk mesh') 

       IF(.NOT. ALLOCATED(IndPart)) ALLOCATE( IndPart(20), IndPartB(20))
       
       NoHerited = 0
       NoConflict = 0
       
       
       DO t=Mesh % NumberOfBulkElements + 1,&
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 
         Element => Mesh % Elements(t) 

         ! Don't inherit from unset elements
         IF( ElementPart(t) == 0 ) CYCLE

         npart = ElementPartitions( Mesh, t, ElementPart, IndPart ) 
         IF( npart == 1 ) CYCLE

         DO k=1,npart         
           BoundPart = IndPart(k)
         
           IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
             DO LeftRight=0,1
               IF( LeftRight == 0 ) THEN
                 Parent => Element % BoundaryInfo % Left
               ELSE
                 Parent => Element % BoundaryInfo % Right
               END IF

               IF( ASSOCIATED( Parent ) ) THEN
                 tB = Parent % ElementIndex
                 npartB = ElementPartitions( Mesh, tB, ElementPart, IndPartB ) 

                 ! Halo already present
                 IF( ANY( IndPartB(1:npartB) == BoundPart ) ) CYCLE

                 IF( npartB == 1 ) THEN
                   ! Halo not yet present, creating one
                   ALLOCATE( Mesh % Halo(tB) % Neighbours(1) )             
                 ELSE
                   ! Halo present, adding one partition to it
                   DEALLOCATE( Mesh % Halo(tB) % Neighbours )
                   ALLOCATE( Mesh % Halo(tB) % Neighbours(npartB) )
                   Mesh % Halo(tb) % Neighbours(1:nPart-1) = IndPartB(2:npartB)
                 END IF

                 NoHerited = NoHerited + 1 
                 Mesh % Halo(tb) % Neighbours(nPartB) = BoundPart

               END IF
             END DO
           END IF
         END DO
       END DO
       
       CALL Info(FuncName,'Number of herited halo bulk elements: '//I2S(NoHerited))
       
     END SUBROUTINE InheritHaloToBulkPart


     ! Inherit partition from a bulkd partition to boundary.
     ! In case of conflict the 1st occurrence prevails.
     !-----------------------------------------------------
     SUBROUTINE InheritBulkToBoundaryPart()
       
       TYPE(Element_t), POINTER :: Element, Parent
       INTEGER :: t, LeftRight, BoundPart, NoHerited, NoConflict, ElemIndx, ParentParts(2)

       CALL Info(FuncName,'Inheriting the bulk partitioning into the boundary mesh') 

       NoHerited = 0
       NoConflict = 0

       DO t=Mesh % NumberOfBulkElements + 1,&
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 
         Element => Mesh % Elements(t) 

         ! Don't set boundaries that are already defined
         ! We want to study conflicted conditions, so comment this away...
         !IF( ElementPart( t ) > 0 ) CYCLE
         
         IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
           ParentParts = -1
           DO LeftRight=0,1
             IF( LeftRight == 0 ) THEN
               Parent => Element % BoundaryInfo % Left
             ELSE
               Parent => Element % BoundaryInfo % Right
             END IF
             
             IF( ASSOCIATED( Parent ) ) THEN               
               ElemIndx = Parent % ElementIndex
               IF(ParentParts(1) == -1 ) THEN
                 ParentParts(1) = ElementPart(ElemIndx)
               ELSE
                 ParentParts(2) = ElementPart(ElemIndx)
               END IF
             END IF
           END DO
           
           IF( ParentParts(1) > -1 ) THEN
             NoHerited = NoHerited + 1
             IF( ElementPart(t) > 0 ) THEN
               IF( .NOT. ANY(ElementPart(t) == ParentParts ) ) THEN
                 NoConflict = NoConflict + 1
                 ElementPart(t) = ParentParts(1)
               END IF
             ELSE                          
               ElementPart(t) = ParentParts(1)
             END IF
           END IF
         END IF
       END DO
       
       CALL Info(FuncName,'Number of herited boundary elements: '//I2S(NoHerited))
       CALL Info(FuncName,'Number of conflicted bulk elements: '//I2S(NoConflict))
       
     END SUBROUTINE InheritBulkToBoundaryPart

     ! Find interface and inherit it to a partition.
     !-----------------------------------------------------
     SUBROUTINE SetInterfacePartition()
       
       TYPE(Element_t), POINTER :: Element
       INTEGER :: i,n,t,eq_id,body_id,mincnt,maxcnt
       INTEGER, ALLOCATABLE :: EqCnt(:)
       LOGICAL, ALLOCATABLE :: EqTag(:)
       
       CALL Info(FuncName,'Inheriting the bulk partitioning into the boundary mesh') 


       n = Mesh % NumberOfNodes
       ALLOCATE( EqTag(n), EqCnt(n) ) 
       EqCnt = 0
       
       DO eq_id = 1, Model % NumberOfEquations

         EqTag = .FALSE.
         
         DO t=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 
           
           Element => Mesh % Elements(t) 
           body_id = Element % BodyId
           IF( body_id <= 0 ) CYCLE
           
           i = ListGetInteger( Model % Bodies(body_id) % Values, 'Equation', Found )
           IF(i /= eq_id) CYCLE
           
           EqTag( Element % NodeIndexes ) = .TRUE.
         END DO

         !PRINT *,'Eq count:',eq_id,COUNT(EqTag)
         

         ! If we found this equation then add the counter for number of equations for the node
         WHERE( EqTag )
           EqCnt = EqCnt + 1
         END WHERE
       END DO

       DO t=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 
         
         Element => Mesh % Elements(t) 
         body_id = Element % BodyId
         IF( body_id <= 0 ) CYCLE
         
         mincnt = MINVAL( EqCnt( Element % NodeIndexes )  )
         maxcnt = MAXVAL( EqCnt( Element % NodeIndexes )  )

         IF( mincnt /= maxcnt ) ElementPart(t) = 1
       END DO

       n = COUNT( ElementPart == 1 )

       CALL Info(FuncName,'Number of elements set to interface partition: '//I2S(n),Level=6)

     END SUBROUTINE SetInterfacePartition
         
     

     ! Merge partitioned boundaries that share even just one node. 
     ! The algorithm works fine when there is a small number of 
     ! partitions on the boundary and requires minimal additional space.
     !------------------------------------------------------------------
     SUBROUTINE MergeJoinedPartitions(IsBoundary)
       LOGICAL :: IsBoundary
       
       TYPE(Element_t), POINTER :: Element
       INTEGER :: t, tstart, tfin, i, j, MaxPart
       LOGICAL, ALLOCATABLE :: PartFlag(:), PartitionCoupling(:,:)
       INTEGER, ALLOCATABLE :: PartMap(:)
       
       CALL Info(FuncName,'Checking for partitions joined by a node',Level=8)

       MaxPart = MAXVAL( ElementPart ) 

       ALLOCATE( PartitionCoupling(MaxPart, MaxPart) )
       PartitionCoupling = .FALSE.
       
       ALLOCATE( PartFlag( Mesh % NumberOfNodes ) ) 

       IF( IsBoundary ) THEN
         tstart = Mesh % NumberOfBulkElements + 1
         tfin = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       ELSE
         tstart = 1
         tfin = Mesh % NumberOfBulkElements
       END IF
         
       DO i = 1, MaxPart
         PartFlag = .FALSE.
         CALL Info(FuncName,'Studying coupling with partition:'//I2S(i),Level=20)

         ! Mark the nodes that are included in this partition, even in one element
         DO t=tstart, tfin
           IF( ElementPart( t ) == i ) THEN 
             Element => Mesh % Elements(t) 
             PartFlag( Element % NodeIndexes ) = .TRUE.
           END IF
         END DO

         ! Studying to which partitions couple to
         ! Only study cases j>i since the coupling is symmetric
         DO t=tstart, tfin
           j = ElementPart( t )
           IF( j > i ) THEN
             IF( PartitionCoupling(i,j) ) CYCLE
             Element => Mesh % Elements(t) 
             IF( ANY( PartFlag( Element % NodeIndexes ) ) ) THEN
               CALL Info(FuncName,&
                   'Joined no is coupling '//I2S(i)//' and '//I2S(j),Level=10)
               PartitionCoupling(i,j) = .TRUE.
               PartitionCoupling(j,i) = .TRUE.
             END IF
           END IF
         END DO
       END DO

       IF(.NOT. ANY( PartitionCoupling ) ) THEN
         CALL Info(FuncName,'Partitions are not coupled',Level=8)
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
                   //I2S(j)//' to become '//I2S(PartMap(i)),Level=8)
               PartMap(j) = PartMap(i)
             END IF
           END IF
         END DO
       END DO
       j = MAXVAL( PartMap ) 
       CALL Info(FuncName,'Number of mapped partitions: '//I2S(j),Level=8)

       ! If we studied bulk elements then make the boundary elements follow the new indexing.
       ! If we stufied boundary elements bulk elements will be studied later.
       IF( .NOT. IsBoundary ) THEN
         tfin = tfin + Mesh % NumberOfBoundaryElements
       END IF
        
       DO t=tstart, tfin
         i = ElementPart( t )
         IF( i > 0 ) THEN
           ElementPart(t) = PartMap(i)
         END IF
       END DO
       CALL Info(FuncName,'Connected boundaries merged',Level=12)
       
     END SUBROUTINE MergeJoinedPartitions


     ! Generates a halo that includes all the couplings between elements
     ! in a set. 
     !------------------------------------------------------------------
     SUBROUTINE GenerateSetHalo()
       
       TYPE(Element_t), POINTER :: Element
       INTEGER :: t, nblk, nbndry, nelem, ntothalo, i, j, MinPart, &
           MaxPart, dPart, ThisPart, nhalo
       
       MaxPart = MAXVAL( ElementPart, PartitionCand ) 
       MinPart = MINVAL( ElementPart, PartitionCand ) 
       
       IF( MaxPart - MinPart == 0 ) THEN
         CALL Info(FuncName,'No need not generate halo within one partition!',Level=12)
         RETURN
       END IF

       CALL Info(FuncName,'Generating halo among partitions '&
           //I2S(MinPart)//' to '//I2S(MaxPart),Level=10 )
       
       nblk = Mesh % NumberOfBulkElements
       nbndry = Mesh % NumberOfBoundaryElements
       nelem = nblk + nbndry

       dPart = ListGetInteger( Params,'Boundary Partition Halo Width',Found )
       
       IF(.NOT. ASSOCIATED( Mesh % Halo ) ) THEN
         ALLOCATE( Mesh % Halo(nelem) )         
         DO t=1,nelem
           NULLIFY( Mesh % Halo(t) % Neighbours )
         END DO
       END IF
       
       ! For this halo type the number of neighbuor partitions in always constant
       ntothalo = 0
       
       DO t=1, nelem
         ThisPart = ElementPart( t ) 
         IF( ThisPart < MinPart .OR. ThisPart > MaxPart ) CYCLE

         ! Create halo only for master elements.
         ! This is useful for mortar & contact BCs. 
         IF( MasterHalo ) THEN           
           IF( .NOT. MasterElement(t) ) CYCLE
         END IF
         
         Element => Mesh % Elements(t) 
         IF(.NOT. ASSOCIATED( Mesh % Halo(t) % Neighbours ) ) THEN
           nhalo = 0
           DO i=MinPart, MaxPart
             IF( ThisPart == i ) CYCLE
             IF( dPart > 0 ) THEN
               IF( ABS(ThisPart-i) > dPart ) CYCLE
             END IF
             nhalo = nhalo + 1
           END DO

           ALLOCATE( Mesh % Halo(t) % Neighbours(nhalo) )
           j = 0
           DO i=MinPart, MaxPart
             ! The owner is always considered, hence don't add that to the halo elements
             IF( ThisPart == i ) CYCLE
             ! The width may be used for example when the partitioning is slices in z-direction
             IF( dPart > 0 ) THEN
               IF( ABS(ThisPart-i) > dPart ) CYCLE
             END IF
             j = j + 1
             Mesh % Halo(t) % Neighbours(j) = i
           END DO
           ntothalo = ntothalo + nhalo
         END IF
       END DO

       CALL Info(FuncName,'Total number of set halo elements: '//I2S(ntothalo),Level=10)
     END SUBROUTINE GenerateSetHalo


     ! Generates a halo for discontinuous Galerkin (DG) method.
     !------------------------------------------------------------------
     SUBROUTINE GenerateDGHalo()
       
       TYPE(Element_t), POINTER :: Element
       INTEGER :: t, nblk, nbndry, nelem, nhalo, ntothalo, ownerpart, i, j, MinPart, MaxPart
       LOGICAL, ALLOCATABLE :: NodeActive(:)
       INTEGER :: TmpNeighbours(20)
       
       MaxPart = MAXVAL( ElementPart )!, PartitionCand ) 
       MinPart = MINVAL( ElementPart )!, PartitionCand ) 

       IF( MaxPart - MinPart == 0 ) THEN
         CALL Info(FuncName,'No need not generate halo within this partition!',Level=12)
         RETURN
       END IF

      
       CALL Info(FuncName,'Generating halo for DG in partitions '&
           //I2S(MinPart)//' to '//I2S(MaxPart),Level=10 )
       
       nblk = Mesh % NumberOfBulkElements
       nbndry = Mesh % NumberOfBoundaryElements
       nelem = nblk + nbndry
       
       IF(.NOT. ASSOCIATED( Mesh % Halo ) ) THEN
         ALLOCATE( Mesh % Halo(nelem) )         
         DO t=1,nelem
           NULLIFY( Mesh % Halo(t) % Neighbours )
         END DO
       END IF
       
       ntothalo = 0

       ALLOCATE( NodeActive( Mesh % NumberOfNodes ) )

       DO ownerpart = MinPart, MaxPart

         ! Mark nodes in partition "ownerpart"
         NodeActive = .FALSE.
         DO t=1, nelem        
           IF( ElementPart( t ) /= ownerpart ) CYCLE
           Element => Mesh % Elements(t) 
           NodeActive( Element % NodeIndexes ) = .TRUE.
         END DO

         ! Find which elements in other partitions have some node in partition "ownerpart"
         DO t=1, nelem
           IF( ElementPart( t ) == ownerpart ) CYCLE                      
           Element => Mesh % Elements(t)
           IF(.NOT. ANY( NodeActive( Element % NodeIndexes ) ) ) CYCLE
           IF(.NOT. ASSOCIATED( Mesh % Halo(t) % Neighbours ) ) THEN
             ! Halo not yet present, creating one
             nhalo = 1
             ALLOCATE( Mesh % Halo(t) % Neighbours(1) )             
           ELSE
             ! Halo present, adding one partition to it
             IF( ANY( Mesh % Halo(t) % Neighbours == ownerpart ) ) CYCLE
             nhalo = SIZE( Mesh % Halo(t) % Neighbours )
             TmpNeighbours(1:nhalo) = Mesh % Halo(t) % Neighbours
             DEALLOCATE( Mesh % Halo(t) % Neighbours )
             ALLOCATE( Mesh % Halo(t) % Neighbours(nhalo+1) )
             Mesh % Halo(t) % Neighbours(1:nhalo) = TmpNeighbours(1:nhalo)
             nhalo = nhalo + 1
           END IF
                      
           Mesh % Halo(t) % Neighbours(nhalo) = ownerpart
           ntothalo = ntothalo + 1           

           !PRINT *,'dghalo:',t,nhalo,ownerpart,ntothalo,Mesh % Halo(t) % Neighbours 
         END DO
       END DO

       CALL Info(FuncName,'Total number of DG halo elements: '//I2S(ntothalo),Level=10)
     END SUBROUTINE GenerateDGHalo


     ! Generates a halo for overlapping bounding boxes.
     !------------------------------------------------------------------
     SUBROUTINE GenerateBBoxHalo()
       
       TYPE(Element_t), POINTER :: Element
       INTEGER :: t, nblk, nbndry, nelem, nhalo, ntothalo, ownerpart, i, j, MinPart, MaxPart
       LOGICAL, ALLOCATABLE :: NodeActive(:)
       INTEGER :: TmpNeighbours(20)

       REAL(KIND=dp) :: Bbox(6),Coord(3),MaxDx
       LOGICAL :: ConstCoord(3)
       INTEGER :: dim
       LOGICAL :: Hit
       
       MaxPart = MAXVAL( ElementPart, PartitionCand ) 
       MinPart = MINVAL( ElementPart, PartitionCand ) 
       MinPart = MAX(1,MinPart)
       
       IF( MaxPart - MinPart == 0 ) THEN
         CALL Info(FuncName,'No need not generate halo within this partition!',Level=12)
         RETURN
       END IF

       dim = 3
       
       CALL Info(FuncName,'Generating halo for bounding boxes in partitions '&
           //I2S(MinPart)//' to '//I2S(MaxPart),Level=10 )
       
       nblk = Mesh % NumberOfBulkElements
       nbndry = Mesh % NumberOfBoundaryElements
       nelem = nblk + nbndry
       
       IF(.NOT. ASSOCIATED( Mesh % Halo ) ) THEN
         ALLOCATE( Mesh % Halo(nelem) )         
         DO t=1,nelem
           NULLIFY( Mesh % Halo(t) % Neighbours )
         END DO
       END IF
       
       ntothalo = 0

       ALLOCATE( NodeActive( Mesh % NumberOfNodes ) )

       DO ownerpart = MinPart, MaxPart

         ! Mark nodes in partition "ownerpart"
         NodeActive = .FALSE.
         DO t=nblk+1, nelem        
           IF( ElementPart( t ) /= ownerpart ) CYCLE
           Element => Mesh % Elements(t) 
           NodeActive( Element % NodeIndexes ) = .TRUE.
         END DO

         !PRINT *,'noactive:',dim,ownerpart,COUNT( NodeActive ) 
         
         BBox(1) = MINVAL( Mesh % Nodes % x, NodeActive ) 
         BBox(2) = MAXVAL( Mesh % Nodes % x, NodeActive ) 
         BBox(3) = MINVAL( Mesh % Nodes % y, NodeActive ) 
         BBox(4) = MAXVAL( Mesh % Nodes % y, NodeActive ) 
         IF( dim > 2 ) THEN
           BBox(5) = MINVAL( Mesh % Nodes % z, NodeActive ) 
           BBox(6) = MAXVAL( Mesh % Nodes % z, NodeActive ) 
         END IF

         MaxDx = MAXVAL( Bbox(2::2)-Bbox(1::2) )           
         ConstCoord = (BBox(2::2)-BBox(1::2) < 1.0e-6*MaxDx)

         !ConstCoord(1:2) = .TRUE.
         
         PRINT *,'Bounding box min:',BBox(1::2)
         PRINT *,'Bounding box max:',BBox(2::2)
         PRINT *,'ConstCoord:',ConstCoord
         
         
         ! Find which elements in other partitions have some node in partition "ownerpart"
         DO t=nblk+1, nelem
           IF( ElementPart( t ) == 0 ) CYCLE
           IF( ElementPart( t ) == ownerpart ) CYCLE                      
           Element => Mesh % Elements(t)

           Hit = .FALSE.
           DO i = 1, Element % TYPE % NumberOfNodes
             j = Element % NodeIndexes(i)
             Coord(1) = Mesh % Nodes % x(j)
             Coord(2) = Mesh % Nodes % y(j)
             Coord(3) = Mesh % Nodes % z(j)

             ! Is node in bounding box? 
             Hit = .TRUE.
             DO k=1,dim
               IF( ConstCoord(k) ) CYCLE
               IF( (Bbox(2*k)-Coord(k))*(Coord(k)-BBox(2*k-1)) < 0.0_dp ) THEN
                 Hit = .FALSE.
                 EXIT
               END IF
             END DO               
             IF( Hit ) EXIT
           END DO

           IF(.NOT. Hit ) CYCLE

           !PRINT *,'hit:',t,ElementPart(t),i,j,Coord, Element % Type % ElementCode

           
           IF(.NOT. ASSOCIATED( Mesh % Halo(t) % Neighbours ) ) THEN
             ! Halo not yet present, creating one
             nhalo = 1
             ALLOCATE( Mesh % Halo(t) % Neighbours(1) )             
           ELSE
             ! Halo present, adding one partition to it
             IF( ANY( Mesh % Halo(t) % Neighbours == ownerpart ) ) CYCLE
             nhalo = SIZE( Mesh % Halo(t) % Neighbours )
             TmpNeighbours(1:nhalo) = Mesh % Halo(t) % Neighbours
             DEALLOCATE( Mesh % Halo(t) % Neighbours )
             ALLOCATE( Mesh % Halo(t) % Neighbours(nhalo+1) )
             Mesh % Halo(t) % Neighbours(1:nhalo) = TmpNeighbours(1:nhalo)
             nhalo = nhalo + 1
           END IF
                      
           Mesh % Halo(t) % Neighbours(nhalo) = ownerpart
           ntothalo = ntothalo + 1           
           
           !PRINT *,'bboxhalo:',t,nhalo,ownerpart,ntothalo,Mesh % Halo(t) % Neighbours 
         END DO
       END DO

       CALL Info(FuncName,'Total number of BBox halo elements: '//I2S(ntothalo),Level=10)
     END SUBROUTINE GenerateBBoxHalo

     
     
     ! Extend partition from an existing bulk partitioning. 
     ! In case of conflict the dominating partitioning prevails.
     ! The routine is written with just a small number of existing
     ! boundary partitions in mind and uses minimal memory. 
     !------------------------------------------------------------
     SUBROUTINE ExtendBoundaryPart(StartFromBC)

       LOGICAL :: StartFromBC

       TYPE(Element_t), POINTER :: Element
       INTEGER :: t, ExtendLayers, NoExtend, ElemIndx, NoHits, TestPart, NumberOfParts
       LOGICAL, ALLOCATABLE :: ActiveNode(:)
       INTEGER, ALLOCATABLE :: RefHits(:)
       INTEGER :: allocstat
       LOGICAL :: BCfirst

       NoExtend = 0
       BCFirst = StartFromBC

       
       ExtendLayers = ListGetInteger( Params,'Partition Mesh Extend Layers', Found ) 
       IF( ExtendLayers <= 0 ) RETURN
       
       CALL Info(FuncName,'Extending boundary meshes by layers: '//I2S(ExtendLayers))

       ALLOCATE( ActiveNode( Mesh % NumberOfNodes ), STAT = allocstat ) 
       IF( allocstat /= 0 ) THEN
         CALL Fatal(FuncName,'Allocation error for ActiveNode')
       END IF


       NumberOfParts = MAXVAL( ElementPart ) 
       IF( NumberOfParts > 1 ) THEN
         CALL Info(FuncName,'Extending boundary to dominating owner among: '//I2S(NumberOfParts))
         ALLOCATE( RefHits( Mesh % NumberOfBulkElements ), STAT = allocstat )
         IF( allocstat /= 0 ) THEN
           CALL Fatal(FuncName,'Allocation error for RefHits')
         END IF
       ELSE
         CALL Info(FuncName,'Extending boundary to all neighbours of 1st partition.')
       END IF


       DO i=1,ExtendLayers

         ! If testing for several partitions then we need a reference for the most hits
         IF( NumberOfParts > 1 ) THEN
           RefHits = 0
         END IF

         DO TestPart = 1, NumberOfParts         

           ! Set the active nodes for the partition under testing
           ActiveNode = .FALSE.
           IF( BCFirst ) THEN
             ! First layer make starting from boundary elements
             DO t=Mesh % NumberOfBulkElements + 1, &
                 Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
               IF( ElementPart( t ) == TestPart ) THEN
                 Element => Mesh % Elements(t) 
                 ActiveNode( Element % NodeIndexes ) = .TRUE.
               END IF
             END DO
             BCFirst = .FALSE.
           ELSE
             ! Thereafter continue from bulk elements
             DO t=1, Mesh % NumberOfBulkElements 
               IF( ElementPart( t ) == TestPart ) THEN
                 Element => Mesh % Elements(t) 
                 ActiveNode( Element % NodeIndexes ) = .TRUE.
               END IF
             END DO
           END IF
             
           ! Count the number of hits for this partition
           ! If larger than the maximum so far set the partition
           ! For just one existing partitioning no checks need to be done. 
           ! Use negative values in a dirty way...
           !--------------------------------------------------------------
           DO t=1, Mesh % NumberOfBulkElements 
             ! These are already decided elements
             IF( ElementPart( t ) > 0 ) CYCLE

             Element => Mesh % Elements(t) 
             NoHits = COUNT( ActiveNode( Element % NodeIndexes ) )
             IF( NoHits == 0 ) CYCLE

             IF( NumberOfParts > 1 ) THEN
               IF( NoHits <= RefHits( t ) ) CYCLE
               RefHits( t ) = NoHits
             END IF
             
             ! So far this is tentative, hence negative sign
             ElementPart( t ) = -TestPart 
           END DO
         END DO

         ! Now include the layer in the new partitioning
         t = COUNT( ElementPart < 0 )
         NoExtend = NoExtend + t
         CALL Info(FuncName,'Layer '//I2S(i)//' with elements: '//I2S(t),Level=8)

         ElementPart = ABS( ElementPart ) 
       END DO
      
       CALL Info(FuncName,'Number of extended bulk elements: '//I2S(NoExtend),Level=6)
       
       DEALLOCATE( ActiveNode ) 
       IF( NumberOfParts > 1 ) DEALLOCATE( RefHits ) 

     END SUBROUTINE ExtendBoundaryPart

      

     ! Initialize sets of boundary elements to be partitioned with various strategies
     !--------------------------------------------------------------------------
     SUBROUTINE InitializeBoundaryElementSet( NumberOfParts )
       
       INTEGER :: NumberOfParts
       
       INTEGER :: i,j,k,n,bc_id
       TYPE(ValueList_t), POINTER :: ValueList
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: Found, NewSet, SeparateBoundarySets        
       INTEGER, ALLOCATABLE :: BCPart(:)
       LOGICAL, ALLOCATABLE :: BCMaster(:)
       INTEGER :: allocstat

       NumberOfParts = 0
       IF( Model % NumberOfBCs == 0 ) RETURN
       
      
       SeparateBoundarySets = ListGetLogical( Params, &
           'Partitioning Separate Boundary Set', Found)

       n = Model % NumberOfBCs
       ALLOCATE( BCPart(n), BCMaster(n), STAT = allocstat )
       IF( allocstat /= 0 ) THEN
         CALL Fatal(FuncName,'Allocation error for BCPart')
       END IF
       BCPart = 0
       BCMaster = .FALSE.

       ! First, set the partition sets enforced by the user
       DO bc_id = 1, Model % NumberOfBCs 
         ValueList => Model % BCs(bc_id) % Values
         k = ListGetInteger( ValueList,'Partition Set',Found)
         IF( .NOT. Found ) CYCLE
         BCPart(bc_id) = k         
         ParameterInd(k) = bc_id
       END DO

       IF( ListGetLogical( Params,'Partition Connected BCs',Found ) ) THEN 
         j = MAXVAL( BCPart ) + 1

         DO bc_id = 1, Model % NumberOfBCs 
           ValueList => Model % BCs(bc_id) % Values

           IF( BCPart(bc_id) > 0 ) CYCLE

           NewSet = ListGetLogical( ValueList, 'Discontinuous Boundary', Found ) .OR. &
               ListGetLogical( ValueList, 'Partition BC',Found) 

           k = ListGetInteger( ValueList, 'Mortar BC',Found)
           IF(.NOT. Found ) k = ListGetInteger( ValueList, 'Contact BC',Found)
           IF(.NOT. Found ) k = ListGetInteger( ValueList, 'Discontinuous BC',Found)
           IF(.NOT. Found ) k = ListGetInteger( ValueList, 'Periodic BC',Found) 
           IF(.NOT. Found ) k = ListGetInteger( ValueList, 'Conforming BC',Found)
           IF(.NOT. Found ) k = ListGetInteger( ValueList, 'Discontinuous BC',Found)           
           
           IF(k>0) NewSet = .TRUE.

           IF( NewSet ) THEN
             BCPart(bc_id) = j
             IF(k > 0) THEN
               BCPart(k) = j
               BCMaster(k) = .TRUE.
             END IF
             ParameterInd(j) = bc_id

             IF( k > 0 ) THEN
               CALL Info('PartitionMeshSerial','Including BCs '&
                   //I2S(bc_id)//' and '//I2S(k)//&
                   ' in set '//I2S(j),Level=10)
             ELSE
               CALL Info('PartitionMeshSerial','Including BC '&
                   //I2S(bc_id)//' in set '//I2S(j),Level=10)
             END IF
               
             IF( SeparateBoundarySets ) j = j + 1             
           END IF
         END DO
       END IF
       
       NumberOfParts = MAXVAL( BCPart ) 
       
       IF( NumberOfParts == 0 ) RETURN
       
       ! Make all elements in the boundary sets to have the designated partition index     
       DO i = Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         Element => Mesh % Elements(i)
         
         DO bc_id=1,Model % NumberOfBCs
           IF ( Element % BoundaryInfo % Constraint == &
               Model % BCs(bc_id) % Tag ) THEN
             IF( BCPart( bc_id ) > 0 ) THEN
               ElementSet( i ) = BCPart( bc_id )
               IF( MasterHalo ) THEN
                 IF( BCMaster( bc_id) ) MasterElement(i) = .TRUE.
               END IF
             END IF
             EXIT
           END IF
         END DO
       END DO

       CALL Info(FuncName,'Number of sets for boundary partitioning: '//I2S(NumberOfParts))
       CALL Info(FuncName,'Number of boundary elements in set: '//I2S(COUNT(ElementSet > 0)))
      
     END SUBROUTINE InitializeBoundaryElementSet


     
     ! Initialize sets of bulk elements to be partitioned with various strategies
     ! By default there is only one set but certain BCs are treated differently.
     !--------------------------------------------------------------------------
     SUBROUTINE InitializeBulkElementSet( NumberOfParts )
       
       INTEGER :: NumberOfParts
       
       INTEGER :: i,j,k,eq_id, bc_id
       TYPE(ValueList_t), POINTER :: ValueList
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: Found, SeparateBoundarySets
       INTEGER :: allocstat
     
       ElementSet = 0 
       
       ALLOCATE( EquationPart( Model % NumberOfEquations ), STAT = allocstat ) 
       IF( allocstat /= 0 ) THEN
         CALL Fatal(FuncName,'Allocation error for EquationPart')
       END IF
 
       EquationPart = 0
       ParameterInd = 0
       
       IF( ListGetLogical( Model % Simulation,'Partition Equation Balance',Found ) ) THEN
         CALL Info(FuncName,'Partitioning each equation separately',Level=5)
         DO eq_id = 1, Model % NumberOfEquations 
           EquationPart(eq_id) = eq_id
           ParameterInd(eq_id) = eq_id
         END DO
       ELSE
         DO eq_id = 1, Model % NumberOfEquations 
           ValueList => Model % Equations(eq_id) % Values
           k = ListGetInteger( ValueList,'Partition Set',Found) 
           IF( k > 0 ) THEN
             EquationPart(eq_id) = k 
             ParameterInd(k) = eq_id
           END IF
         END DO
       END IF

       NumberOfParts = MAXVAL( EquationPart )
       
       IF( NumberOfParts == 0 ) THEN
         WHERE( ElementPart == 0 ) ElementSet = NumberOfBoundarySets + 1
         NumberOfParts = 1
         RETURN
       END IF
                                 
       DO i = 1, Mesh % NumberOfBulkElements 
         IF( ElementPart(i) > 0 ) CYCLE

         Element => Mesh % Elements(i)
         eq_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Equation', Found )
         IF( eq_id > 0 ) j = EquationPart( eq_id ) 
         IF( j == 0 ) THEN
           CALL Fatal(FuncName,'"Partition Set" in Equation '//I2S(eq_id)//' is undefined!')
         END IF

         ElementSet(i) = NumberOfBoundarySets + j
       END DO
       
       CALL Info(FuncName,'Number of sets for bulk partitioning: '//I2S(NumberOfParts))
       
     END SUBROUTINE InitializeBulkElementSet


     
     ! Partition the nodes that have the correct ElementSet using various strategies.
     !------------------------------------------------------------------------
     SUBROUTINE PartitionMeshPart(SetNo, LocalParams, IsBoundary, PartOffset, PartOffsetBC )
       
      INTEGER :: SetNo
      TYPE(ValueList_t), POINTER :: LocalParams       
      LOGICAL :: IsBoundary      
      INTEGER :: PartOffset, PartOffsetBC
     
      CHARACTER(:), ALLOCATABLE :: CoordTransform, SetMethod
      LOGICAL :: GotCoordTransform, SetNodes
      INTEGER :: SumPartitions, NoPartitions, NoCand
      LOGICAL :: Found, GotMethod
      INTEGER :: i,j,NoCandElements
      REAL(KIND=dp) :: BoundaryFraction
      INTEGER, POINTER :: PartDivs(:)
      
      PartitionCand = ( ElementSet == SetNo )
      n = Mesh % NumberOfBulkElements      
      NoCandElements = COUNT( PartitionCand ) 
     
      CALL Info(FuncName,'Doing element set: '//I2S(SetNo))
      CALL Info(FuncName,'Number of elements in set: '//I2S(NoCandElements))

      IF( NoCandElements == 0 ) THEN
        CALL Info(FuncName,'No element in set, doing nothing.',Level=10)
        RETURN
      END IF
        
      GotMethod = .FALSE.
      IF( ASSOCIATED( LocalParams) ) THEN
        SetMethod = ListGetString( LocalParams,'Partitioning Method',GotMethod)
      END IF
      IF(.NOT. GotMethod ) THEN
        IF( IsBoundary ) THEN
          SetMethod = ListGetString( Params,'Boundary Partitioning Method',GotMethod)
        ELSE
          SetMethod = ListGetString( Params,'Partitioning Method',GotMethod)
        END IF
      END IF

      NoPartitions = 0
      Found = .FALSE.
      IF( ASSOCIATED( LocalParams) ) THEN
        NoPartitions = ListGetInteger( LocalParams,'Number of Partitions',Found )        
      END IF
      
      IF(.NOT. Found ) THEN
        IF( IsBoundary ) THEN
          NoPartitions = ListGetInteger( Params,'Boundary Number of Partitions',Found)
        ELSE
          NoPartitions = ListGetInteger( Params,'Number Of Partitions',Found)
          IF(Found) NoPartitions = NoPartitions - PartOffsetBC         
        END IF
      END IF

      IF( NoPartitions == 0 ) THEN
        IF(GotMethod) THEN
          IF( SetMethod == 'uniform' .OR. SetMethod == 'directional') THEN
            PartDivs => ListGetIntegerArray( Params,'Partitioning Divisions',Found )
            IF( Found ) THEN
              NoPartitions = 1
              DO i=1,SIZE(PartDivs)
                NoPartitions = NoPartitions * MAX(1,PartDivs(i))
              END DO
            END IF
          END IF
        END IF
      END IF
      
      IF(NoPartitions == 0 ) THEN
        CALL Info(FuncName,'Number of partitions not defined!',Level=10)
        IF( IsBoundary ) THEN
          CALL Info(FuncName,'Defaulting to one partition for BCs',Level=10)
          NoPartitions = 1          
        ELSE
          NoPartitions = ParEnv % PEs - PartOffsetBC
          CALL Info(FuncName,'Defaulting rest of partitions for the equation: '//I2S(NoPartitions),Level=5)
        END IF
      END IF
      
      ! We have parallel case but asked too many partitions
      IF( ParEnv % PEs > 1 .AND. NoPartitions + PartOffsetBC > ParEnv % PEs ) THEN
        NoPartitions = ParEnv % PEs - PartOffsetBC
        CALL Info(FuncName,'Reducing partitions for the following split to: '//I2S(NoPartitions),Level=5)
      END IF

      BoundaryFraction = ListGetCReal( Params,'Boundary Partitioning Maximum Fraction',Found)
      IF( Found ) THEN
        IF( NoCandElements <= BoundaryFraction * Mesh % NumberOfBulkElements ) THEN
          NoPartitions = 1 
          WRITE(Message,'(A,ES12.3)') 'Number of elements in set below critical limit: ',BoundaryFraction
          CALL Info(FuncName,Message )
        ELSE
          NoPartitions = CEILING(  Mesh % NumberOfBulkElements / ( BoundaryFraction * NoCandElements ) )
        END IF
      END IF
      
      IF( NoPartitions == 1 ) THEN
        CALL Info(FuncName,'Partitions set to one, doing nothing.',Level=10)
        PartOffset = PartOffset + 1
        IF( IsBoundary ) PartOffsetBC = PartOffset
        WHERE( PartitionCand ) ElementPart = PartOffset 
        PartitionCand = .FALSE.
        RETURN
      END IF

      IF( .NOT. GotMethod ) THEN
        CALL Fatal(FuncName,'Undefined > Partitioning Method < e.g. "uniform", "directional"')
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

      ! There may be various coordinate transformation (e.g. to cylindrical coordinates)
      ! that allow for different partitions when using the geometries partitioning routines. 
      IF( GotCoordTransform ) THEN
        CALL CoordinateTransformation( Mesh, CoordTransform, Params, &
            IrreversibleTransformation = .FALSE. )
      END IF

      CALL Info(FuncName,'Using partitioning method: '//TRIM(SetMethod))
      
      CALL Info(FuncName,'Using partitioning offset: '//I2S(PartOffset))

      
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
          IF( IsBoundary ) THEN
            CALL Fatal('PartitionMeshPart','Zoltan interface not yet applicable to boundary partitioning!')
          END IF
          CALL Zoltan_Interface( Model, Mesh, .TRUE., NoPartitions, PartitionCand )
#else
          CALL Fatal(FuncName,'Partition with Zoltan not available!')
#endif 
          
        CASE DEFAULT
          CALL Fatal(FuncName,'Unspecificed partitioning: '//TRIM(SetMethod))
          
      END SELECT

      ! Add offset related to previous partitioning routines to the elements of this set.
      ! This only relates to bulk elements!
      i = Mesh % NumberOfBulkElements
      IF( PartOffset > 0 ) THEN
        WHERE( PartitionCand(1:i) ) ElementPart(1:i) = ElementPart(1:i) + PartOffset
      END IF
      PartOffset = MAXVAL( ElementPart )
      IF( IsBoundary ) PartOffsetBC = PartOffset
      
      CALL Info(FuncName,'Partitioning of set finished',Level=10)      
      
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

      CALL Info(FuncName,'Creating parallel neighbour information')

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
        IF( .NOT. ASSOCIATED( NeighbourList(k) % Neighbours ) ) CYCLE
        l = SIZE(NeighbourList(k) % Neighbours)
        lmax = MAX( lmax, l )
        lsum = lsum + l
      END DO

      CALL Info(FuncName,'Maximum number of partitions for a node: '//I2S(lmax))

      WRITE(Message,'(A,F8.3)') 'Average number of partitions for a node: ',1.0_dp*lsum/n
      CALL Info(FuncName,Message)
      
    END SUBROUTINE CreateNeighbourList

  END SUBROUTINE PartitionMeshSerial
  
END MODULE MeshPartition
