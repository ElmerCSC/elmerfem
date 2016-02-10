!/*****************************************************************************/
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
! *  Authors: Jouni Malinen, Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2000
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!>  These routines are for parallel version of ELMER solver.
!>  Subroutines for MPI-communication
!------------------------------------------------------------------------------
!> \ingroup ElmerLib
!> \{


#include "huti_fdefs.h"

MODULE SParIterComm

  USE Types
#ifdef USE_ISO_C_BINDINGS
  USE LoadMod
#endif
  USE SParIterGlobals

  IMPLICIT NONE

  INCLUDE "mpif.h"

  TYPE Buff_t
    REAL(KIND=dp), ALLOCATABLE :: rbuf(:)
  END TYPE Buff_t

  TYPE iBuff_t
    INTEGER, ALLOCATABLE :: ibuf(:)
  END TYPE iBuff_t

  TYPE vBuff_t
    INTEGER, ALLOCATABLE :: ind(:)
    REAL(KIND=dp), ALLOCATABLE :: vec(:)
  END TYPE vBuff_t
CONTAINS

!-----------------------------------------------------------------
  SUBROUTINE CheckBuffer( n )
!-----------------------------------------------------------------
     INTEGER :: n, i, sz, ierr
     LOGICAL :: isfine
     INTEGER(KIND=1), ALLOCATABLE, SAVE :: send_buffer(:)
     INTERFACE 
        SUBROUTINE MPI_BUFFER_ATTACH( buf, n, ierr)
           INTEGER(KIND=1) :: buf(*)
           INTEGER :: n, ierr
        END SUBROUTINE MPI_BUFFER_ATTACH

        SUBROUTINE MPI_BUFFER_DETACH( buf, n, ierr)
           INTEGER(KIND=1) :: buf(*)
           INTEGER :: n, ierr
        END SUBROUTINE MPI_BUFFER_DETACH
     END INTERFACE

     sz = MAX( 4*n, 2**20 )

     isfine = ALLOCATED( send_buffer )
     IF ( isfine ) isfine = (sz <= SIZE(send_buffer))
     IF ( isfine ) RETURN

     IF ( ALLOCATED(send_buffer) ) THEN
        i = SIZE(send_buffer)
        CALL MPI_BUFFER_DETACH( send_buffer, i, ierr )
        DEALLOCATE( send_buffer )
     END IF

     ALLOCATE( send_buffer(sz), stat=i )
     IF ( i/= 0 ) THEN
       CALL Fatal( 'CheckBuffer', 'Alloc failed' )
     END IF

     CALL MPI_BUFFER_ATTACH( send_buffer, sz, ierr )

!-----------------------------------------------------------------
   END SUBROUTINE CheckBuffer
!-----------------------------------------------------------------


!------------------------------------------------------------------------
!> Initialize parallel execution environment
!-----------------------------------------------------------------------
   FUNCTION ParCommInit( ) RESULT ( ParallelEnv ) 
!-----------------------------------------------------------------------
     TYPE (ParEnv_t), POINTER :: ParallelEnv

    ! Local variables

    INTEGER :: ierr
    INTEGER :: req, prov

    !******************************************************************

    ParallelEnv => ParEnv

    ParEnv % MyPE = 0
    ParEnv % PEs  = 1
    ParEnv % ActiveComm = MPI_COMM_WORLD

    ierr = 0
#ifdef _OPENMP
    req = MPI_THREAD_FUNNELED
    CALL MPI_Init_Thread(req, prov, ierr)
    IF (prov < req) THEN
      WRITE( Message, '(A,I0,A,I0,A,I0,A)' ) &
              'MPI Thread Initialization failed! (req=', req,&
              ', prov=', prov, &
              ', ierr=', ierr, ')'
      CALL Fatal( 'ParCommInit', Message )
    END IF
#else
    CALL MPI_INIT( ierr )
#endif
    IF ( ierr /= 0 ) RETURN

    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, ParEnv % PEs, ierr )
    IF ( ierr /= 0 ) THEN
       CALL MPI_Finalize( ierr )
    ELSE
       CALL MPI_COMM_RANK( MPI_COMM_WORLD, ParEnv % MyPE, ierr )
       OutputPE = ParEnv % MyPe

       WRITE( Message, * ) 'Initialize #PEs: ', ParEnv % PEs
       CALL Info( 'ParCommInit', Message, Level=5 )
    
       IF ( ierr /= 0 ) THEN
          WRITE( Message, * ) 'MPI Initialization failed ! (ierr=', ierr, ')'
          CALL Fatal( 'ParCommInit', Message )
       END IF

       Parenv % NumOfNeighbours = 0
       ParEnv % Initialized = .TRUE.
    END IF
!-----------------------------------------------------------------------
  END FUNCTION ParCommInit
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Initialize parallel execution environment
!-----------------------------------------------------------------------
  SUBROUTINE ParEnvInit( SPMatrix, ParallelInfo, SourceMatrix )
!-----------------------------------------------------------------------

    TYPE(SparIterSolverGlobalD_t) :: SPMatrix
    TYPE (ParallelInfo_t) :: ParallelInfo
    TYPE(Matrix_t) :: SourceMatrix
!-----------------------------------------------------------------------
    CALL FindActivePEs( ParallelInfo, SourceMatrix )
    SPMatrix % ParEnv = ParEnv
    SPMatrix % ParEnv % ActiveComm = SourceMatrix % Comm
!-----------------------------------------------------------------------
  END SUBROUTINE ParEnvInit
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
  SUBROUTINE SParIterActive(L)
!-----------------------------------------------------------------------
    LOGICAL :: L
    INTEGER :: ierr
!-----------------------------------------------------------------------
    LOGICAL, ALLOCATABLE :: Active(:)
    ALLOCATE( Active(ParEnv % PEs) )

    IF ( .NOT. ASSOCIATED(ParEnv % Active) ) &
       ALLOCATE( ParEnv % Active(ParEnv % PEs) )

    ParEnv % Active = .FALSE.
    Active = .FALSE.
    Active(ParEnv % MYPe+1) = L
    CALL MPI_ALLREDUCE(Active,ParEnv % Active,ParEnv % PEs, &
         MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
    DEALLOCATE( Active )
!-----------------------------------------------------------------------
  END SUBROUTINE SParIterActive
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Find active PEs using ParallelInfo % NeighbourList
!-----------------------------------------------------------------------
  SUBROUTINE FindActivePEs( ParallelInfo, SourceMatrix, JustNeighbours )
!-----------------------------------------------------------------------
    LOGICAL, OPTIONAL :: JustNeighbours
    TYPE(Matrix_t) :: SourceMatrix
    TYPE(ParallelInfo_t) :: ParallelInfo
!-----------------------------------------------------------------------
    TYPE NlistEntry_t
      INTEGER :: e1,e2
      TYPE(NlistEntry_t), POINTER :: Next
    END TYPE NlistEntry_t

    TYPE Nlist_t
      TYPE(NlistEntry_t), POINTER :: Head
    END TYPE Nlist_t

     TYPE(NlistEntry_t), POINTER :: ptr, ptr1
     TYPE(Nlist_t), ALLOCATABLE :: NeighList(:)

    ! Local variables
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    
    INTEGER :: i, j, k, m, ii, n, sz, ierr, proc, MinActive
    INTEGER, ALLOCATABLE :: Active(:), buf(:)
    LOGICAL :: L, Interf
#ifdef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: tstart, tend,s
#else
    REAL(KIND=dp) :: RealTime, tstart, tend,s
#endif
    LOGICAL(KIND=1), ALLOCATABLE :: NeighAll(:,:)
    !******************************************************************

    IF ( .NOT. ASSOCIATED(ParEnv % Active) ) THEN
      ALLOCATE(ParEnv % Active(ParEnv % PEs))
      ParEnv % Active = .TRUE.
    END IF

    ALLOCATE(ParEnv % IsNeighbour(ParEnv % PEs))
    ParEnv % IsNeighbour(:)  = .FALSE.
    ParEnv % NumOfNeighbours = 0

    !------------------------------------------------------------------
    ! Count the number of real neighbours for this partition
    !------------------------------------------------------------------
    DO i=1,SourceMatrix % NumberOfRows
       IF ( ASSOCIATED(ParallelInfo % NeighbourList(i) % Neighbours) ) THEN
         DO j=1,SIZE(ParallelInfo % NeighbourList(i) % Neighbours)
            proc = ParallelInfo % NeighbourList(i) % Neighbours(j)
            IF ( ParEnv % Active(proc+1).AND.proc/=ParEnv % MYpe) &
               ParEnv % IsNeighbour(proc+1) = .TRUE.
         END DO
       END IF
    END DO

    !------------------------------------------------------------------
    ! Sync all neighbouring information
    !------------------------------------------------------------------
    CALL CheckBuffer( ParEnv % PEs**2 + MPI_BSEND_OVERHEAD )

    ALLOCATE( Active(ParEnv % PEs), NeighList(Parenv % PEs) )
    DO MinActive=0,ParEnv % PEs-1
      IF ( ParEnv % Active(MinActive+1) ) EXIT
    END DO
    Active = -1
    n = 0
    DO i=1,ParEnv % PEs
      IF ( ParEnv % IsNeighbour(i) ) THEN
         n = n + 1
         Active(n) = i-1
      END IF
    END DO

    IF (Parenv % myPE /= MinActive ) THEN
      CALL MPI_BSEND( n, 1, MPI_INTEGER, MinActive, &
                800, MPI_COMM_WORLD, ierr )
      IF ( n>0 ) THEN
        CALL MPI_BSEND( Active, n, MPI_INTEGER, MinActive, &
                 801, MPI_COMM_WORLD, ierr )
      END IF

      CALL MPI_RECV( n, 1, MPI_INTEGER, MinActive, &
             802, MPI_COMM_WORLD, status, ierr )
      IF ( n>0 ) THEN
        CALL MPI_RECV( Active, n, MPI_INTEGER, MinActive, &
               803, MPI_COMM_WORLD, status, ierr )
        DO i=1,n
          ParEnv % IsNeighbour(Active(i)+1) = .TRUE.
        END DO
      END IF
    ELSE
      DO i=1,ParEnv % PEs
        NeighList(i) % head => NULL()
      END DO

      DO i=1,n
        CALL AddToNlist( NeighList(MinActive+1), Active(i), 0 )
        CALL AddToNlist( NeighList(Active(i)+1), MinActive, 0 )
      END DO

      DO i=1,COUNT(ParEnv % Active)-1
        CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
              800, MPI_COMM_WORLD, status, ierr )
        IF ( n>0 ) THEN
          proc = status(MPI_SOURCE)
          CALL MPI_RECV( Active, n, MPI_INTEGER, proc, &
              801, MPI_COMM_WORLD, status, ierr )
          DO k=1,n
            CALL AddToNList( NeighList(Active(k)+1), proc, 0 )
            CALL AddToNList( NeighList(proc+1), Active(k), 0 )
            IF (Active(k)==MinActive) ParEnv % IsNeighbour(proc+1) = .TRUE.
          END DO
        END IF
      END DO

      DO i=1,Parenv % PEs
        IF ( ParEnv % Active(i) .AND. i-1/=MinActive ) THEN
          n = 0
          ptr => NeighList(i) % head
          DO WHILE( ASSOCIATED(ptr) )
            n = n + 1
            Active(n) = ptr % e1
            ptr1 => ptr
            ptr => ptr % Next
            DEALLOCATE(ptr1)
          END DO

          CALL MPI_BSEND( n, 1, MPI_INTEGER, i-1, &
                 802, MPI_COMM_WORLD, ierr )
          IF ( n>0 ) THEN
            CALL MPI_BSEND( Active, n, MPI_INTEGER, i-1, &
                   803, MPI_COMM_WORLD, ierr )
          END IF
        END IF
      END DO
    END IF
    ParEnv % isNeighbour(ParEnv % MyPE+1) = .FALSE.
    ParEnv % NumOfNeighbours = COUNT(ParEnv % isNeighbour)

    IF ( PRESENT(JustNeighbours) ) THEN
      IF ( JustNeighbours ) RETURN
    END IF

    ! Make sure that the owner of a matrix dof is one that
    ! has it active:
    ! -----------------------------------------------------
    DO i=1,ParEnv % Pes
      NeighList(i) % Head => NULL()
    END DO
    DEALLOCATE( Active )
    ALLOCATE( Active(SIZE(ParallelInfo % Interface)) )

    IF ( .NOT. SourceMatrix % DGMatrix ) THEN
      DO i=1,SIZE(SourceMatrix % Perm)
        ii = SourceMatrix % Perm(i)

        Interf=.FALSE.
        IF(ii>0) THEN
          Active(ii) = HUGE(i)
          Interf=ParallelInfo % Interface(ii)
        END IF
        IF ( Interf) THEN
          sz = SIZE(ParallelInfo % NeighbourList(ii) % Neighbours)
          DO j=1,sz
            k = ParallelInfo % NeighbourList(ii) % Neighbours(j)
            IF ( k == ParEnv % Mype ) THEN
              Active(ii) = j
            ELSE IF ( ParEnv % Active(k+1) ) THEN
              CALL AddtoNlist( Neighlist(k+1), ParallelInfo % GlobalDOFs(ii), ii)
            END IF
          END DO
        END IF
      END DO

      sz = 0
      n  = 0
      DO i=1,ParEnv % PEs
        IF ( .NOT. ParEnv % IsNeighbour(i) ) CYCLE
        j = 0
        ptr => NeighList(i) % Head
        DO WHILE(ASSOCIATED(ptr))
          j = j + 1
          ptr => ptr % next
        END DO
        n = n + j
        sz = MAX( sz, j )
      END DO
      ALLOCATE( buf(2*sz) )

      CALL CheckBuffer( 4*n+4*Parenv % NumOfNeighbours+MPI_BSEND_OVERHEAD )

      DO i=1,ParEnv % PEs
        IF ( .NOT. ParEnv % IsNeighbour(i) ) CYCLE
        j = 0
        ptr => NeighList(i) % Head
        DO WHILE(ASSOCIATED(ptr))
          j = j + 1
          buf(j) = ptr % e1
          j = j + 1
          buf(j) = ptr % e2
          ptr1 => ptr
          ptr => ptr % next
          DEALLOCATE( ptr1 )
        END DO
        CALL MPI_BSEND(j,1,MPI_INTEGER,i-1, 20000,MPI_COMM_WORLD,status,ierr)
        IF (j>0) CALL MPI_BSEND(buf,j,MPI_INTEGER,i-1,20001,MPI_COMM_WORLD,status,ierr)
      END DO
      DEALLOCATE( NeighList, buf )

      m = SIZE(ParallelInfo % GlobalDOFs)
      DO i=1,ParEnv % NumOfNeighbours
        CALL MPI_RECV( sz,1,MPI_INTEGER,MPI_ANY_SOURCE,20000,MPI_COMM_WORLD,status,ierr)
        IF (sz>0 ) THEN
          proc = status(MPI_SOURCE)
          ALLOCATE( buf(sz) )
          CALL MPI_RECV( buf,sz,MPI_INTEGER,proc,20001,MPI_COMM_WORLD,status,ierr)
          DO j=1,sz,2
            IF ( buf(j+1)>0 ) THEN
              k = SearchNode( ParallelInfo, buf(j), Order=SourceMatrix % Perm )
              IF ( k>0  ) THEN
                sz = SIZE(ParallelInfo % NeighbourList(k) % Neighbours)
                DO n=1,sz
                  IF ( ParallelInfo % Neighbourlist(k) % Neighbours(n)==proc) THEN
                    Active(k) = MIN(Active(k),n)
                    EXIT
                  END IF
                END DO
              END IF
            END IF
          END DO
          DEALLOCATE(buf)
        END IF
      END DO

      DO i=1,SIZE(SourceMatrix % Perm)
        ii = SourceMatrix % Perm(i)
        Interf = .FALSE.
        IF(ii>0) Interf=ParallelInfo % Interface(ii)
        IF ( Interf ) THEN
          sz = SIZE(ParallelInfo % NeighbourList(ii) % Neighbours)
          IF ( Active(ii)>1 .AND. Active(ii)<=sz ) THEN
            n = ParallelInfo % NeighbourList(ii) % Neighbours(Active(ii))
            ParallelInfo % NeighbourList(ii) % Neighbours(Active(ii)) = &
               ParallelInfo % NeighbourList(ii) % Neighbours(1)
            ParallelInfo % NeighbourList(ii) % Neighbours(1) = n
          END IF
        END IF
      END DO
    END IF

    DEALLOCATE( Active )

CONTAINS

    SUBROUTINE AddToNList( Nlist, e1, e2 )
      TYPE(Nlist_t) :: Nlist
      INTEGER :: e1,e2
      TYPE(NlistEntry_t), POINTER :: ptr, ptr1, prev

      ptr => Nlist % Head

      IF ( .NOT. ASSOCIATED(ptr) ) THEN
         ALLOCATE(ptr)
         ptr % e1 = e1
         ptr % e2 = e2
         ptr % Next => Null()
         Nlist % head => ptr
         RETURN
      END IF

      prev => Null()
      DO WHILE( ASSOCIATED(ptr) )
         IF ( ptr % e1 >= e1 ) EXIT
         prev  => ptr
         ptr => ptr % next
      END DO

      IF ( ASSOCIATED(ptr) ) THEN
        IF ( ptr % e1 == e1) RETURN
      END IF

      ALLOCATE(ptr1)
      ptr1 % e1 = e1
      ptr1 % e2 = e2
      ptr1 % next => ptr
      IF ( ASSOCIATED(prev) ) THEN
        prev % next => ptr1
      ELSE
        Nlist % head => ptr1
      END IF
    END SUBROUTINE AddToNList

!-----------------------------------------------------------------------
  END SUBROUTINE FindActivePEs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE AddToCommonList( list, ENTRY )
!-----------------------------------------------------------------------
! Helper subroutine for SParIterGlobalNumbering
! Adds integer Entry to integer pointer List(:)
!-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, POINTER :: list(:)
    INTEGER :: ENTRY
!-----------------------------------------------------------------------
    INTEGER, POINTER :: ptmp(:)
    INTEGER :: itmp,i
!-----------------------------------------------------------------------
    IF( ASSOCIATED(list)) THEN
       itmp = SIZE(list)
       ALLOCATE(ptmp(itmp+1))
       ptmp(1:itmp) = list
       ptmp(itmp+1) = entry
       DEALLOCATE(list)
       list => ptmp
    ELSE
       ALLOCATE(ptmp(1))
       ptmp(1) = entry
       list => ptmp
    END IF
!-----------------------------------------------------------------------
  END SUBROUTINE AddToCommonList
!-----------------------------------------------------------------------  

!-----------------------------------------------------------------------
 SUBROUTINE SyncNeighbours(ParEnv)
!-----------------------------------------------------------------------
  TYPE(ParEnv_t) :: ParEnv
!-----------------------------------------------------------------------
  LOGICAL :: L
  INTEGER :: i, ierr, status(MPI_STATUS_SIZE)
!-----------------------------------------------------------------------
  DO i=1,ParEnv % PEs
    IF(Parenv % Mype==i-1) CYCLE
    IF(ParEnv % Active(i)) &
      CALL MPI_BSEND( ParEnv % IsNeighbour(i),1, &
                 MPI_LOGICAL,i-1,1410,MPI_COMM_WORLD,ierr)
  END DO

  DO i=1,ParEnv % PEs
    IF(Parenv % Mype==i-1) CYCLE
    IF(ParEnv % Active(i)) THEN
      CALL MPI_RECV( L,1,MPI_LOGICAL,i-1,1410,MPI_COMM_WORLD,status,ierr)
      IF(L) ParEnv % IsNeighbour(i) = .TRUE.
    END IF
  END DO
  Parenv % IsNeighbour(Parenv % myPE+1) = .FALSE.
  Parenv % NumOfNeighbours = COUNT(Parenv % IsNeighbour)
!-----------------------------------------------------------------------
  END SUBROUTINE SyncNeighbours
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------  
  FUNCTION MeshNeighbours(Mesh,IsNeighbour) RESULT(num)
!-----------------------------------------------------------------------  
      INTEGER :: i,j,num
      TYPE(Mesh_t) :: Mesh
      LOGICAL :: IsNeighbour(:)

      IsNeighbour = .FALSE.
      DO i=1,Mesh % Nodes % NumberOfNodes
        IF ( Mesh % ParallelInfo % Interface(i) ) THEN
          DO j=1,SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
            IsNeighbour(Mesh % ParallelInfo % &
              NeighbourList(i) % Neighbours(j)+1) = .TRUE.
          END DO
        END IF
      END DO
      IsNeighbour(ParEnv % myPE+1) = .FALSE.
      num = COUNT(IsNeighbour)
!-----------------------------------------------------------------------  
  END FUNCTION MeshNeighbours
!-----------------------------------------------------------------------  


!-----------------------------------------------------------------------  
   SUBROUTINE SParEdgeNumbering( Mesh, Allmesh )
     USE GeneralUtils
     TYPE(Mesh_t) :: Mesh
     LOGICAL, OPTIONAL :: Allmesh
!-----------------------------------------------------------------------

     TYPE Edgen_t 
       INTEGER :: n = 0
       LOGICAL :: interface = .FALSE.
       INTEGER, POINTER :: neighbours(:) => NULL()
     END TYPE Edgen_t

     TYPE(Edgen_t)   :: copyedgen
     TYPE(Element_t) :: copyedge
     TYPE(Edgen_t), ALLOCATABLE :: edgen(:), cccc(:)

     INTEGER, ALLOCATABLE :: owned_edges(:), owned_edges2(:), parentnodes(:,:), &
             tosend(:), toreceive(:), gorder(:), grevorder(:),buf(:), Iedges(:)

     INTEGER  :: p,q,datasize, ierr,status(MPI_STATUS_SIZE)
     INTEGER, ALLOCATABLE ::  gindices(:)
     INTEGER, POINTER :: list1(:),list2(:), commonlist(:), gdofs(:)

     INTEGER :: i,j,k,l,m,n,mm,edofs,maxedofs,retry,bd,gi,ind,proc,n_intf,k_intf

     LOGICAL, POINTER :: ig(:)
     TYPE(NeighbourList_t), POINTER :: nb(:)

     LOGICAL :: AllM, Intf, Found
     LOGICAL, POINTER :: IsNeighbour(:)
     TYPE(Element_t), POINTER :: Element

#ifdef USE_ISO_C_BINDINGS
     real(kind=dp) :: tt
#else
     real(kind=dp) :: realtime,tt
#endif
!-------------------------------------------------------------------------------

    IF ( .NOT. ASSOCIATED(Mesh % Edges) ) RETURN

    CALL CheckBuffer(1000000)

    AllM = .FALSE.
    IF ( PRESENT(Allmesh) ) AllM=Allmesh
    IF ( AllM .OR. .NOT. ASSOCIATED(ParEnv % IsNeighbour) ) THEN
      AllM = .TRUE.
      ALLOCATE(IsNeighbour(ParEnv % PEs))
      n = MeshNeighbours(Mesh,IsNeighbour)
    ELSE
      AllM = .FALSE.
      IsNeighbour => ParEnv % IsNeighbour
    END IF

    n = Mesh % NumberOfEdges

    ALLOCATE(  &
         Edgen(n), parentnodes(n,2),   &
         owned_edges ( ParEnv % PEs ), &
         owned_edges2( ParEnv % PEs ), &
         tosend( ParEnv % PEs ),       &
         toreceive( ParEnv % PEs ) )

    owned_edges  = 0
    owned_edges2 = 0
    tosend       = 0
    toreceive    = 0
    parentnodes  = 0

    Edgen(:) % n = 0
    gdofs => Mesh % ParallelInfo % GlobalDOFs
    ig => Mesh % ParallelInfo % Interface
    nb => Mesh % ParallelInfo % NeighbourList

    !
    ! Find neighbours and parent nodes for all new interface edges:
    ! =============================================================
    !

    ! Loop over edges:
    !-----------------
    DO i = 1, n
      Element => Mesh % Edges(i)
      Edgen(i) % Interface = .FALSE.

      Intf = ALL(ig(Element % NodeIndexes(1:2)))
      IF ( Mesh % MeshDim==2 ) THEN
        Intf = Intf.AND..NOT. &
             (ASSOCIATED(Element % BoundaryInfo % Left) .AND. &
              ASSOCIATED(Element % BoundaryInfo % Right))
      END IF

      IF (Intf) THEN
        !
        ! This is perhaps an interface edge:
        !-----------------------------------

        commonlist => Null() ! intersection of list1 and list2

        p = Element % NodeIndexes(1)
        q = Element % NodeIndexes(2)
        IF (gdofs(p)>gdofs(q)) THEN
          parentnodes(i,1) =  p
          parentnodes(i,2) =  q
        ELSE
          parentnodes(i,1) =  q
          parentnodes(i,2) =  p
        ENDIF

        ! This is the list of PEs sharing parent node 1:
        list1 => nb(parentnodes(i,1)) % Neighbours

        ! This is the list of PEs sharing parent node 2:
        list2 => nb(parentnodes(i,2)) % Neighbours
        !
        ! Determine the intersection of the PE-lists:
        !--------------------------------------------
        DO p = 1,SIZE(list1)
           DO q = 1,SIZE(list2)
              IF( list1(p) == list2(q)) &
                 CALL AddToCommonList( commonlist, list1(p) ) 
           END DO
        END DO
        
        ! in 3D add the rest of the PE's
        ! ----------------------------------
        IF ( Mesh % MeshDim==3 ) THEN
          DO p=1,SIZE(list1)
            IF (.NOT. ANY(List1(p) == CommonList) ) THEN
              CALL AddToCommonList(commonlist, list1(p)) 
            END IF
          END DO
          DO p=1,SIZE(list2)
            IF (.NOT. ANY(List2(p) == CommonList) ) THEN
              CALL AddToCommonList(commonlist, list2(p)) 
            END IF
          END DO
        END IF

        !
        ! Now, we should have a list of PEs common to the parents:
        !----------------------------------------------------------
        IF ( SIZE(commonlist)<2 ) THEN
          DEALLOCATE(CommonList);
          CommonList => Null()
        END IF

        IF( ASSOCIATED(commonlist)) THEN
          Edgen(i) % Interface = .TRUE.
          !
          ! Edge i owner is owner of one of the parent nodes:
          !--------------------------------------------------
          q = MAX(list1(1),list2(1))
          DO p=1,SIZE(commonlist)
            IF (commonlist(p)==q) THEN
              j=commonlist(1)
              commonlist(1)=commonlist(p)
              commonlist(p)=j
              EXIT
            END IF
          END DO
          IF ( commonlist(1)/=q ) THEN
            q = MIN(list1(1),list2(1))
            DO p=1,SIZE(commonlist)
              IF (commonlist(p)==q) THEN
                j=commonlist(1)
                commonlist(1)=commonlist(p)
                commonlist(p)=j
                EXIT
              END IF
            END DO
          END IF
          Edgen(i) % Neighbours => commonlist

          !
          ! Finalize by sorting the parent table:
          !---------------------------------------
          CALL Sort(2,parentnodes(i,:))
        END IF
      END IF

      IF ( .NOT. Edgen(i) % Interface ) THEN
        ALLOCATE( Edgen(i) % Neighbours(1) )
        Edgen(i) % Neighbours(1) = ParEnv % MyPE
      END IF
    END DO

10  CONTINUE

    ! Count #edges owned by each of the PEs:
    ! --------------------------------------
    owned_edges=0
    j = Parenv % MyPE
    DO i=1,n
      Edgen(i) % n = 0
      l = Edgen(i) % Neighbours(1) 
      IF (l==j) owned_edges(j+1)=owned_edges(j+1)+1
    END DO

    owned_edges2=0
    CALL MPI_ALLREDUCE( owned_edges, owned_edges2, ParEnv % PEs, &
        MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
    !
    ! Number the edges:
    !------------------
    ! Start numbering from index k, above the #edges owned
    ! by PEs numbered below us.
    ! -----------------------------------------------------

    ! Start numbering from index k:
    j = ParEnv % MyPE
    k = SUM( owned_edges2(1:j) )
    DO i = 1, n
       l = Edgen(i) % Neighbours(1)
       IF( l == j ) THEN
         k = k+1
         Edgen(i) % n = k
       END IF
    END DO

    !
    ! Some new interface nodes may still lack a global number, if
    ! we are not the owner of these nodes. Hence, we'll have to make
    ! a query to all neighbour PEs to find these numbers.
    !===============================================================
    !
    ! Determine the size of data to exchange between neighbours:
    !-----------------------------------------------------------
    isNeighbour(ParEnv % myPE+1)=.FALSE.

    n_intf=COUNT(Edgen(:) % Interface)
    ALLOCATE(Iedges(n_intf));
    n_intf=0
    DO i=1,n
      IF(Edgen(i) % Interface) THEN
        n_intf=n_intf+1
        Iedges(n_intf)=i
      END IF
    END DO

    tosend = 0
    toreceive = 0
    DO k_intf = 1, n_intf
      i=Iedges(k_intf)
      j = Edgen(i) % Neighbours(1)
      IF( j /= ParEnv % MyPE ) THEN
         toreceive(j+1) = toreceive(j+1)+1
      ELSE
         DO j = 2,SIZE(Edgen(i) % Neighbours)
           k = Edgen(i) % Neighbours(j)
           IF(k>=0) &
             tosend(k+1) = tosend(k+1)+1
         END DO
      END IF
    END DO

    !
    ! Distribute the interface data:
    !--------------------------------
    DO i = 1, ParEnv % PEs
       IF( i-1 == ParEnv % MyPE ) CYCLE
       IF( .NOT. IsNeighbour(i) ) CYCLE

       CALL MPI_BSEND( tosend(i), 1, MPI_INTEGER, i-1, &
               100, MPI_COMM_WORLD, ierr)
       IF( tosend(i) < 1 ) CYCLE

       datasize = 3*tosend(i)
       ALLOCATE( gindices(datasize) )
       k = 1
       DO k_intf = 1, n_intf
         l=Iedges(k_intf)
         m = Edgen(l) % Neighbours(1)
         IF( m /= ParEnv % MyPE ) CYCLE

         DO m = 2,SIZE(Edgen(l) % Neighbours)
            mm = Edgen(l) %  Neighbours(m)
            IF( mm+1 /= i ) CYCLE
            gindices(k) = Edgen(l) % n
            gindices(k+1) = gdofs(parentnodes(l,1))
            gindices(k+2) = gdofs(parentnodes(l,2))
            k = k+3
         END DO
       END DO

       CALL MPI_BSEND( gindices, DataSize, MPI_INTEGER, &
           i-1, 101, MPI_COMM_WORLD, ierr )
       DEALLOCATE( gindices )
    END DO

    !
    ! Retrieve missing interface node numbers from our neighbours:
    !--------------------------------------------------------------
    tosend = 0

    DO i = 1, COUNT(isNeighbour)
       CALL MPI_RECV( DataSize, 1, MPI_INTEGER, &
          MPI_ANY_SOURCE, 100, MPI_COMM_WORLD, status, ierr )

       proc=status(MPI_SOURCE)

       ALLOCATE( gindices(3*DataSize) )

       IF( datasize > 0 ) THEN
         CALL MPI_RECV( gindices, 3*DataSize, MPI_INTEGER, &
            proc, 101, MPI_COMM_WORLD, status, ierr )

         DO l = 1, datasize
           Found=.FALSE.
           m = 3*l-2
           DO k_intf = 1,n_intf
             k = Iedges(k_intf)
             IF( GDofs(ParentNodes(k,1)) == GIndices(m+1) .AND. &
                 GDofs(ParentNodes(k,2)) == GIndices(m+2) .OR.  &
                 GDofs(ParentNodes(k,1)) == GIndices(m+2) .AND. &
                 GDofs(parentNodes(k,2)) == GIndices(m+1) ) THEN
               Found=.TRUE.
               Edgen(k) % n = GIndices(m) ! global DOF number
               EXIT
             END IF
           END DO

           IF (.NOT.Found) THEN
             tosend(proc+1) = tosend(proc+1)+1
             gindices(tosend(proc+1))=gindices(m)
           END IF
         END DO
       END IF


       !
       ! we don't have these edges, ask to remove us from
       ! share lists:
       ! ------------------------------------------------
       CALL MPI_BSEND( tosend(proc+1), 1, MPI_INTEGER, &
            proc, 200, MPI_COMM_WORLD, ierr )

       IF ( tosend(proc+1)>0 ) THEN
         CALL MPI_BSEND( gindices, tosend(proc+1), MPI_INTEGER, &
            proc, 201, MPI_COMM_WORLD, ierr )
       END IF
       DEALLOCATE(gindices)
    END DO

    DO i = 1, COUNT(isNeighbour)

      CALL MPI_RECV( DataSize, 1, MPI_INTEGER, &
       MPI_ANY_SOURCE, 200, MPI_COMM_WORLD, status, ierr )

      proc=status(MPI_SOURCE)
      IF( datasize < 1 ) CYCLE

      ALLOCATE( gindices(DataSize) )
      CALL MPI_RECV( gindices, DataSize, MPI_INTEGER, &
         proc, 201, MPI_COMM_WORLD, status, ierr )

      DO k_intf=1,n_intf
        j=Iedges(k_intf) 
        DO k=1,Datasize
          IF (Edgen(j) % n==gindices(k)) THEN
            l = SIZE(Edgen(j) % Neighbours)
            DO m=2,l
              IF (Edgen(j) % Neighbours(m)==proc) THEN
                Edgen(j)%Neighbours(m)=-ABS(Edgen(j)%Neighbours(m))-1
                EXIT
              END IF
            END DO
            EXIT
          END IF
        END DO
      END DO
      DEALLOCATE(gindices)
    END DO


    !
    ! Sync neighbour tables from owners
    !--------------------------------------------------------------
    tosend=0
    DO k_intf=1,n_intf
      j = Iedges(k_intf)
      IF (Edgen(j) % n==0 ) CYCLE
      IF (Edgen(j) % Neighbours(1)/=ParEnv % MyPE) CYCLE

      mm=COUNT(Edgen(j) % Neighbours>=0) 
      commonlist=>Edgen(j) % Neighbours
      ALLOCATE(Edgen(j) % Neighbours(mm))

      Edgen(j) % Neighbours(1)=ParEnv % MyPE
      l=1
      DO k=2,SIZE(commonlist)
        m=commonlist(k)
        IF (m>=0) THEN
          tosend(m+1) = tosend(m+1)+mm
          l=l+1
          Edgen(j) % Neighbours(l)=commonlist(k)
        END IF
      END DO
      DEALLOCATE(commonlist)
      IF(mm==1) Edgen(j) % Interface=.FALSE.
    END DO

    n_intf=0
    DO i=1,n
      IF(Edgen(i) % Interface) THEN
        n_intf=n_intf+1
        Iedges(n_intf)=i
      END IF
    END DO

    DO i = 1, ParEnv % PEs
      IF( i-1 == ParEnv % MyPE ) CYCLE
      IF( .NOT. IsNeighbour(i) ) CYCLE

      CALL MPI_BSEND( tosend(i), 1, MPI_INTEGER, &
           i-1, 300, MPI_COMM_WORLD, ierr )

      IF (tosend(i)<=0 ) CYCLE

      ALLOCATE(Gindices(tosend(i)))
      l=0
      DO k_intf=1,n_intf
        j=Iedges(k_intf)
        IF (Edgen(j)%n==0) CYCLE
        IF (Edgen(j) % Neighbours(1)/=ParEnv % MyPE) CYCLE
        IF (ANY(Edgen(j) % Neighbours==i-1)) THEN
          l=l+1
          gindices(l)=-Edgen(j) % n
          DO k=2,SIZE(Edgen(j) % Neighbours)
            l=l+1
            gindices(l) = Edgen(j) % Neighbours(k)
          END DO
        END IF
      END DO

      CALL MPI_BSEND( gindices, tosend(i), MPI_INTEGER, &
             i-1, 301, MPI_COMM_WORLD, ierr )

      DEALLOCATE(gindices)
    END DO

    DO i = 1,COUNT(isNeighbour)

       CALL MPI_RECV( DataSize, 1, MPI_INTEGER, &
        MPI_ANY_SOURCE, 300, MPI_COMM_WORLD, status, ierr )

       proc=status(MPI_SOURCE)

       IF( datasize > 0 ) THEN
         ALLOCATE( gindices(DataSize) )
         CALL MPI_RECV( gindices, DataSize, MPI_INTEGER, &
            proc, 301, MPI_COMM_WORLD, status, ierr )

         j=0
         DO WHILE(j<Datasize)
           j=j+1
           k=-gindices(j)
           Found=.FALSE.
           DO k_intf=1,n_intf
             l = Iedges(k_intf)
             IF (Edgen(l) % n==0) CYCLE
             IF (Edgen(l) % n==k) THEN
               DO m=j+1,Datasize
                 IF (gindices(m)<0) EXIT
               END DO
               DEALLOCATE(Edgen(l) % Neighbours)
               ALLOCATE(Edgen(l) % Neighbours(m-j))
               Edgen(l) % Neighbours(1)=proc
               Edgen(l) % Neighbours(2:)=gindices(j+1:m-1)
               j = m-1
               Found=.TRUE.
               EXIT
             END IF
           END DO
         END DO
         DEALLOCATE(gindices)
       END IF
    END DO

    !
    ! Check that all edges have global numbers. If we haven't
    ! received a new global number, the 'owner' supposedly
    ! doesn't have the edge (even if it has the edge nodes!),
    ! so remove the previous owner and try a new one from the
    ! list of PEs sharing the edge nodes.
    !------------------------------------------------------
    m=0
    DO i=1,n
      IF ( Edgen(i) % n == 0 ) THEN
        m=m+1
        k = SIZE(Edgen(i) % Neighbours)
        ALLOCATE(commonlist(k-1))
        commonlist = Edgen(i) % neighbours(2:)
        DEALLOCATE( Edgen(i) % neighbours )
        Edgen(i) % Neighbours => commonlist
        IF (k<=2) Edgen(i) % Interface=.FALSE.
      END IF
    END DO

    CALL MPI_ALLREDUCE( m, retry, 1, &
        MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

    ! if some of the edges don't have global numbers, retry.
    ! ------------------------------------------------------
    DEALLOCATE(Iedges)

    IF ( retry>0 ) GOTO 10

20  CONTINUE

    ! Deallocate temp arrays:
    !------------------------

    DEALLOCATE( owned_edges, owned_edges2, &
      parentnodes, toreceive )
    !-------------------------------------------------------------------
    ! Now, the numbering should be unique. It remains to sort the edges.
    !-------------------------------------------------------------------

    !
    ! Sort our own edges according to ascending global order:
    ! -------------------------------------------------------
    ALLOCATE(gdofs(n), grevorder(n), gorder(n))
    DO i=1,n
      gdofs(i) = Edgen(i) % n
      gorder(i) = i
      Mesh % Edges(i) % GelementIndex=gdofs(i)
    END DO

    CALL SortI( n, gdofs, gorder )

    DO i=1,n
      grevorder(gorder(i)) = i
    END DO

    ! Adapt edge indexes of the bulk & boundary elements
    ! to the changed order of edges-array:
    ! --------------------------------------------------
    DO i=1,Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(i)
      IF (ASSOCIATED(Element % EdgeIndexes)) THEN
        j = Element % Type % NumberOfEdges
        Element % EdgeIndexes(1:j)=grevorder(Element % EdgeIndexes(1:j))
      END IF
    END DO

    DO i=1,Mesh % NumberOfFaces
      Element => Mesh % Faces(i)
      j = Element % Type % NumberOFEdges
      Element % EdgeIndexes(1:j)=grevorder(Element % EdgeIndexes(1:j))
    END DO

    !
    ! Order the whole of the nodal structure according to the changed
    ! order of the global numbers:
    ! ----------------------------------------------------------------
    DO i=1,n
      k = gorder(i)
      copyedge = Mesh % Edges(i)
      Mesh % Edges(i) = Mesh % Edges(k)
      Mesh % Edges(k) = copyedge

      copyedgen = Edgen(i)
      Edgen(i) = Edgen(k)
      Edgen(k) = copyedgen

      j = grevorder(i)
      grevorder(i)=grevorder(k)
      grevorder(k) = j
      gorder(grevorder(k)) = k
    END DO

    ALLOCATE( Mesh % ParallelInfo % EdgeInterface(n), &
              Mesh % ParallelInfo % EdgeNeighbourList(n) )

    DO i=1,n
      Mesh % ParallelInfo % EdgeInterface(i) = Edgen(i) % Interface
      Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours => &
             Edgen(i) % Neighbours
      Mesh % Edges(i) % ElementIndex = i
    END DO

    ! Max the dofs for the shared edges:
    ! ----------------------------------
    n_intf=COUNT(Edgen(:) % Interface)
    ALLOCATE(Iedges(n_intf));
    n_intf=0
    DO i=1,n
      IF(Edgen(i) % Interface) THEN
        n_intf=n_intf+1
        Iedges(n_intf)=i
      END IF
    END DO

    tosend=0
    DO k_intf=1,n_intf
      i = Iedges(k_intf)
      DO j=1,SIZE(Edgen(i) % Neighbours)
        k = Edgen(i) % Neighbours(j)
        IF(k/=ParEnv % MyPE) &
          tosend(k+1)=tosend(k+1)+1
      END DO
    END DO

    ALLOCATE(buf(MAXVAL(tosend)*2))

    DO i = 1, ParEnv % PEs
      IF( i-1 == ParEnv % MyPE ) CYCLE
      IF( .NOT. IsNeighbour(i) ) CYCLE
      mm=0
      DO k_intf=1,n_intf
        l = Iedges(k_intf)
        DO j=1,SIZE(Edgen(l) % Neighbours)
          k = Edgen(l) % Neighbours(j)
          IF(k==i-1) THEN
            mm=mm+1
            buf(2*mm-1) = gdofs(l)
            buf(2*mm) = Mesh % Edges(l) % BDOFs
          END IF
        END DO
      END DO

      CALL MPI_BSEND( mm, 1, MPI_INTEGER, &
           i-1, 400, MPI_COMM_WORLD, ierr )

      IF (mm<=0 ) CYCLE

      CALL MPI_BSEND( buf, mm*2, MPI_INTEGER, &
           i-1, 401, MPI_COMM_WORLD, ierr )
    END DO

    DEALLOCATE(buf)

    DO i = 1, COUNT(isNeighbour)

      CALL MPI_RECV( DataSize, 1, MPI_INTEGER, &
       MPI_ANY_SOURCE, 400, MPI_COMM_WORLD, status, ierr )

      proc = status(MPI_SOURCE)

      IF( DataSize > 0 ) THEN
        ALLOCATE( buf(DataSize*2) )
        CALL MPI_RECV( buf, DataSize*2, MPI_INTEGER, &
           proc, 401, MPI_COMM_WORLD, status, ierr )

        DO j=1,Datasize
          gi = buf(2*j-1)
          bd = buf(2*j)
          DO k_intf=1,n_intf
            k=Iedges(k_intf)
            IF(gdofs(k)==gi) EXIT
          END DO
          if ( k_intf<=n_intf) &
            Mesh % Edges(k) % BDOFs=MAX(Mesh % Edges(k) % BDOFs,bd)
        END DO
        DEALLOCATE(buf)
      END IF
    END DO

    DEALLOCATE( gorder, grevorder, gdofs, edgen, tosend )
    IF ( AllM ) DEALLOCATE(IsNeighbour)
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
!--------------------------------------------------------------------  
  END SUBROUTINE SParEdgeNumbering
!--------------------------------------------------------------------  


!-----------------------------------------------------------------------  
   SUBROUTINE SParFaceNumbering( Mesh, Allmesh )
     USE GeneralUtils
     TYPE(Mesh_t) :: Mesh
     LOGICAL, OPTIONAL :: Allmesh
!-----------------------------------------------------------------------

     TYPE Facen_t 
       INTEGER :: n = 0
       LOGICAL :: interface = .FALSE.
       INTEGER, POINTER :: neighbours(:) => NULL()
     END TYPE Facen_t

     TYPE(Facen_t)   :: copyfacen
     TYPE(Element_t) :: copyface
     TYPE(Facen_t), ALLOCATABLE :: facen(:)

     TYPE flist_t
       INTEGER, POINTER :: Pes(:)
     END TYPE flist_t
     TYPE(flist_t) :: list(4)

     INTEGER, ALLOCATABLE :: owned_faces(:), owned_faces2(:), Ifaces(:), &
      parentnodes(:,:), tosend(:), toreceive(:), gorder(:), grevorder(:),buf(:)

     INTEGER, POINTER :: commonlist(:), gindices(:), gdofs(:)

     INTEGER :: i,j,k,l,m,n,mm,nd,edofs,maxedofs,ind,gi,bd,proc, &
                p,q,r,datasize, ierr,status(MPI_STATUS_SIZE), k_intf,n_intf

     LOGICAL, POINTER :: ig(:)
     TYPE(NeighbourList_t), POINTER :: nb(:)

     LOGICAL :: AllM, Intf
     LOGICAL, POINTER :: IsNeighbour(:)
     TYPE(Element_t), POINTER :: Element, Edge

#ifdef USE_ISO_C_BINDINGS
     real(kind=dp) :: tt
#else
     real(kind=dp) :: realtime,tt
#endif
!-------------------------------------------------------------------------------

    IF ( .NOT. ASSOCIATED(Mesh % Faces) ) RETURN

    CALL CheckBuffer(1000000)

    AllM = .FALSE.
    IF ( PRESENT(Allmesh) ) AllM=Allmesh
    IF ( AllM .OR. .NOT. ASSOCIATED(ParEnv % IsNeighbour) ) THEN
      AllM = .TRUE.
      ALLOCATE(IsNeighbour(ParEnv % PEs))
      n = MeshNeighbours(Mesh,IsNeighbour)
    ELSE
      AllM = .FALSE.
      IsNeighbour => ParEnv % IsNeighbour
    END IF

    n = Mesh % NumberOfFaces

    ALLOCATE(  &
         Facen(n), parentnodes(n,4),   &
         owned_faces ( ParEnv % PEs ), &
         owned_faces2( ParEnv % PEs ), &
         tosend( ParEnv % PEs ),       &
         toreceive( ParEnv % PEs ) )

    owned_faces  = 0
    owned_faces2 = 0
    tosend       = 0
    toreceive    = 0
    parentnodes  = 0

    Facen(:) % n = 0
    gdofs => Mesh % ParallelInfo % GlobalDOFs
    ig => Mesh % ParallelInfo % Interface
    nb => Mesh % ParallelInfo % NeighbourList


    !
    ! Find neighbours and parent nodes for all new interface faces:
    ! =============================================================
    !
    ! Loop over faces:
    !-----------------
    DO i = 1, n
      Element => Mesh % Faces(i)
      Facen(i) % Interface = .FALSE.

      nd = Element % Type % ElementCode/100

      Intf = ALL(ig(Element % NodeIndexes(1:nd)))
      Intf = Intf .AND..NOT. &
           (ASSOCIATED(Element % BoundaryInfo % Left) .AND. &
            ASSOCIATED(Element % BoundaryInfo % Right))

      IF ( Intf ) THEN
        !
        ! This is an perhaps an interface face:
        !--------------------------------------

        commonlist => Null() ! intersection of pe lists
        DO j=1,nd
          l = Element % NodeIndexes(j)
          parentnodes(i,j) =  l
          list(j) % pes => Mesh % ParallelInfo % NeighbourList(l) % Neighbours
        END DO
        !
        ! Determine the intersection of the PE-lists:
        !--------------------------------------------
        DO p = 1,SIZE(list(1) % pes)
          j = 1
          DO k = 2,nd
            DO q=1,SIZE(list(k) % pes)
              IF( list(1)% pes(p)==list(k) % pes(q) ) THEN
                j=j+1;
                EXIT
              END IF
            END DO
          END DO
          IF (j==nd) CALL AddToCommonList(commonlist, list(1) % Pes(p))
        END DO

        !
        ! Now, we should have a list of PEs common to the parents:
        !----------------------------------------------------------
        IF( ASSOCIATED(commonlist) ) THEN
          IF ( SIZE(commonlist)<2 ) THEN
            DEALLOCATE(commonlist); Commonlist=>Null()
          END IF
        END IF

        IF( ASSOCIATED(commonlist) ) THEN
          Facen(i) % Interface = .TRUE.
          !
          ! Face i is given to owner of max of the parent nodes
          !-----------------------------------------------------
          q = 0
          DO p=1,nd
            q = MAX(q,list(p) % Pes(1))
          END DO

          DO p=1,SIZE(commonlist)
            IF (commonlist(p)==q) THEN
              j=commonlist(1)
              commonlist(1)=commonlist(p)
              commonlist(p)=j
              EXIT
            END IF
          END DO
          Facen(i) % Neighbours => commonlist

          !
          ! Finalize by sorting the parent table:
          !---------------------------------------
          CALL Sort(nd,parentnodes(i,:))
        END IF
      END IF

      IF ( .NOT. Facen(i) % Interface ) THEN
        ALLOCATE(Facen(i) % Neighbours(1))
        Facen(i) % Neighbours(1) = ParEnv % myPE
      END IF
    END DO

10  CONTINUE

    ! count faces owned by each PE:
    ! -----------------------------
    owned_faces=0
    j = Parenv % MyPE
    DO i=1,n
      Facen(i) % n = 0
      l = Facen(i) % Neighbours(1) 
      IF (l==j) owned_faces(j+1)=owned_faces(j+1)+1
    END DO

    owned_faces2=0
    CALL MPI_ALLREDUCE( owned_faces, owned_faces2, ParEnv % PEs, &
        MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
    !
    ! Number the faces:
    !------------------
    ! Start numbering from index k, above the #faces owned
    ! by PEs numbered below us.
    ! -----------------------------------------------------
    j = ParEnv % MyPE
    k = SUM( owned_faces2(1:j) )
    DO i = 1, n
       l = Facen(i) % Neighbours(1)
       IF( l == j ) THEN
         k = k+1
         Facen(i) % n = k
       END IF
    END DO

    !
    ! Some new interface nodes may still lack a global number, if
    ! we are not the owner of these nodes. Hence, we'll have to make
    ! a query to all neighbour PEs to find these numbers.
    !===============================================================
    !
    ! Determine the size of data to exchange between neighbours:
    !-----------------------------------------------------------
    isNeighbour(ParEnv % myPE+1)=.FALSE.

    n_intf=COUNT(Facen(:) % Interface)
    ALLOCATE(Ifaces(n_intf));
    n_intf=0
    DO i=1,n
      IF(Facen(i) % Interface) THEN
        n_intf=n_intf+1
        Ifaces(n_intf)=i
      END IF
    END DO

    tosend = 0
    toreceive = 0
    DO k_intf = 1, n_intf
      i=Ifaces(k_intf)
      j = Facen(i) % Neighbours(1)
      IF( j /= ParEnv % MyPE ) THEN
         toreceive(j+1) = toreceive(j+1)+1
      ELSE
         DO j = 2,SIZE(Facen(i) % Neighbours)
           k = Facen(i) % Neighbours(j)
           tosend(k+1) = tosend(k+1)+1
         END DO
      END IF
    END DO

    !
    ! Distribute the interface data:
    !--------------------------------
    DO i = 1, ParEnv % PEs
      IF( i-1 == ParEnv % MyPE ) CYCLE
      IF( .NOT. IsNeighbour(i) ) CYCLE

      CALL MPI_BSEND( tosend(i), 1, MPI_INTEGER, i-1, &
            100, MPI_COMM_WORLD, ierr)
      IF(tosend(i)<1) CYCLE

      datasize = 4*tosend(i)
      ALLOCATE(gindices(datasize))
      k = 1
      DO l = 1, n
         m = Facen(l) % Neighbours(1)
         IF( m /= ParEnv % MyPE ) CYCLE

         DO m = 2,SIZE(Facen(l) % Neighbours)
            mm = Facen(l) %  Neighbours(m)
            IF( mm+1 .NE. i ) CYCLE
            gindices(k) = facen(l) % n
            gindices(k+1) = gdofs(parentnodes(l,1))
            gindices(k+2) = gdofs(parentnodes(l,2))
            gindices(k+3) = gdofs(parentnodes(l,3))
            k = k+4
         END DO
      END DO

      CALL MPI_BSEND( gindices, DataSize, MPI_INTEGER, &
            i-1, 101, MPI_COMM_WORLD, ierr )
      DEALLOCATE( gindices )
    END DO

    !
    ! Retrieve missing interface node numbers from our neighbours:
    !--------------------------------------------------------------
    tosend = 0
    DO i = 1, COUNT(isNeighbour)
      CALL MPI_RECV( DataSize, 1, MPI_INTEGER, MPI_ANY_SOURCE, 100, &
                  MPI_COMM_WORLD, status, ierr )
      proc=status(MPI_SOURCE)

      ALLOCATE( gindices(4*DataSize) )
      IF( datasize > 0 ) THEN
        CALL MPI_RECV( gindices, 4*DataSize, MPI_INTEGER, &
           proc, 101, MPI_COMM_WORLD, status, ierr )

        DO l = 1, datasize
          m = 4*l-3
          DO k_intf = 1,n_intf
            k=Ifaces(k_intf)
            IF( Facen(k) % n>0 ) CYCLE
            IF( gdofs(parentnodes(k,1)) == gindices(m+1) .AND. &
                gdofs(parentnodes(k,2)) == gindices(m+2) .AND. &
                gdofs(parentnodes(k,3)) == gindices(m+3) ) THEN
               Facen(k) % n = gindices(m) ! global DOF number
               EXIT
            END IF
          END DO
          IF (k>n) THEN
            tosend(proc+1) = tosend(proc+1)+1
            gindices(tosend(proc+1))=gindices(m)
          END IF
        END DO
      END IF
      !
      ! we don't have these edges, ask to remove us from
      ! share lists:
      ! ------------------------------------------------
      CALL MPI_BSEND( tosend(proc+1), 1, MPI_INTEGER, &
           proc, 200, MPI_COMM_WORLD, ierr )

      IF ( tosend(proc+1)>0 ) THEN
        CALL MPI_BSEND( gindices, tosend(proc+1), MPI_INTEGER, &
                proc, 201, MPI_COMM_WORLD, ierr )
      END IF
      DEALLOCATE(gindices)
    END DO

    DO i = 1, COUNT(isNeighbour)
       CALL MPI_RECV( DataSize, 1, MPI_INTEGER, &
          MPI_ANY_SOURCE, 200, MPI_COMM_WORLD, status, ierr )
       IF( datasize < 1 ) CYCLE

       proc=status(MPI_SOURCE)

       ALLOCATE( gindices(DataSize) )
       CALL MPI_RECV( gindices, DataSize, MPI_INTEGER, &
          proc, 201, MPI_COMM_WORLD, status, ierr )

       DO k_intf=1,n_intf
         j = Ifaces(k_intf)
         DO k=1,Datasize
           IF (Facen(j) % n==gindices(k)) THEN
             l = SIZE(Facen(j) % Neighbours)
             ALLOCATE(commonlist(l-1))
             p = 0
             DO m=1,l
               IF (Facen(j) % Neighbours(m)==proc) CYCLE
               p = p + 1
               commonlist(p)=Facen(j) % Neighbours(m)
             END DO
             DEALLOCATE(Facen(j) % neighbours)
             Facen(j) % Neighbours => commonlist
             IF ( l==2 ) Facen(j) % Interface=.FALSE.
             EXIT
           END IF
         END DO
       END DO
       DEALLOCATE(gindices)
    END DO

    n_intf=0
    DO i=1,n
      IF(Facen(i) % Interface) THEN
        n_intf=n_intf+1
        Ifaces(n_intf)=i
      END IF
    END DO
    !
    ! Check that all faces have global numbers. If we haven't
    ! received a new global number, the 'owner' supposedly
    ! doesn't have the face (even if it has the face nodes!),
    ! so remove the previous owner and try a new one from the
    ! list of PEs sharing the face nodes.
    !------------------------------------------------------
    m=0
    DO i=1,n
      IF ( Facen(i) % n == 0 ) THEN
        m=m+1
        k = SIZE(Facen(i) % Neighbours)
        ALLOCATE(commonlist(k-1))
        commonlist = Facen(i) % Neighbours(2:)
        DEALLOCATE(Facen(i) % Neighbours)
        Facen(i) % Neighbours => commonlist
	IF ( k==2 ) Facen(i) % Interface=.FALSE.
      END IF
    END DO

    CALL MPI_ALLREDUCE( m, mm, 1, &
        MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

    ! if some of the faces miss a global number, retry.
    ! -------------------------------------------------
    DEALLOCATE(Ifaces)

    IF ( mm>0 ) GOTO 10

    ! Deallocate temp arrays:
    !------------------------
    DEALLOCATE( owned_faces, owned_faces2, &
       parentnodes, toreceive )

    !-------------------------------------------------------------------
    ! Now, the numbering should be unique. It remains to sort the faces.
    !-------------------------------------------------------------------

    ! Sort our own faces according to ascending
    ! global order:
    ! -----------------------------------------
    ALLOCATE( gdofs(n), grevorder(n), gorder(n) )
    gorder = (/ (i,i=1,n) /)
    DO i=1,n
      gdofs(i) = Facen(i) % n
      Mesh % Faces(i) % GelementIndex=gdofs(i)
    END DO

    CALL SortI(n,gdofs,gorder)
    DO i=1,n
      grevorder(gorder(i)) = i
    END DO

    ! adapt face indexes of the bulk elements
    ! to the changed order of faces-array:
    ! ---------------------------------------
    DO i=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(i)
      j = Element % Type % NumberOfFaces
      Element % FaceIndexes(1:j)=grevorder(Element % FaceIndexes(1:j))
    END DO

    !
    ! Order the whole of the facal structure
    ! according to the changed order of the 
    ! global numbers:
    ! --------------------------------------
    DO i=1,n
      k = gorder(i)
      copyface = Mesh % Faces(i)
      Mesh % Faces(i) = Mesh % Faces(k)
      Mesh % Faces(k) = copyface

      copyfacen = facen(i)
      facen(i) = facen(k)
      facen(k) = copyfacen

      j = grevorder(i)
      grevorder(i) = grevorder(k)
      grevorder(k) = j
      gorder(grevorder(k)) = k
    END DO


    DO i=1,Mesh % NumberOfEdges
      Mesh % Edges(i) % BoundaryInfo % Left  => Null()
      Mesh % Edges(i) % BoundaryInfo % Right => Null()
    END DO

    DO i=1,n
      Element => Mesh % Faces(i)
      DO j=1,Element % Type % NumberOfEdges
        Edge => Mesh % Edges(Element % EdgeIndexes(j))
        IF (.NOT.ASSOCIATED(Edge % BoundaryInfo % Left)) THEN
          Edge % BoundaryInfo % Left  => Mesh % Faces(i)
        ELSE IF (.NOT.ASSOCIATED(Edge % BoundaryInfo % Right)) THEN
          Edge % BoundaryInfo % Right => Mesh % Faces(i)
        END IF
      END DO
      Element % ElementIndex = i
    END DO

    ALLOCATE( Mesh % ParallelInfo % FaceInterface(n), &
              Mesh % ParallelInfo % FaceNeighbourList(n) )

    DO i=1,n
      Mesh % ParallelInfo % FaceInterface(i) = Facen(i) % Interface
      Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours => &
             Facen(i) % Neighbours
    END DO

    n_intf=COUNT(Facen(:) % Interface)
    ALLOCATE(Ifaces(n_intf));
    n_intf=0
    DO i=1,n
      IF(Facen(i) % Interface) THEN
        n_intf=n_intf+1
        Ifaces(n_intf)=i
      END IF
    END DO

    tosend=0
    DO k_intf=1,n_intf
      i=Ifaces(k_intf)
      DO j=1,SIZE(Facen(i) % Neighbours)
        k = Facen(i) % Neighbours(j)
        IF(k/=ParEnv % MyPE) &
          tosend(k+1)=tosend(k+1)+1
      END DO
    END DO

    ALLOCATE(buf(MAXVAL(tosend)*2))

    DO i = 1, ParEnv % PEs
      IF( i-1 == ParEnv % MyPE ) CYCLE
      IF( .NOT. IsNeighbour(i) ) CYCLE

      mm=0
      DO k_intf=1,n_intf
        l=Ifaces(k_intf)
        DO j=1,SIZE(Facen(l) % Neighbours)
          k = Facen(l) % Neighbours(j)
          IF(k==i-1) THEN
            mm=mm+1
            buf(2*mm-1) = gdofs(l)
            buf(2*mm) = Mesh % Faces(l) % BDOFs
          END IF
        END DO
      END DO

      CALL MPI_BSEND( mm, 1, MPI_INTEGER, &
           i-1, 400, MPI_COMM_WORLD, ierr )

      IF (mm<=0 ) CYCLE

      CALL MPI_BSEND( buf, mm*2, MPI_INTEGER, &
           i-1, 401, MPI_COMM_WORLD, ierr )
    END DO

    DEALLOCATE(buf)

    DO i = 1,COUNT(isNeighbour)
      CALL MPI_RECV( DataSize, 1, MPI_INTEGER, &
         MPI_ANY_SOURCE, 400, MPI_COMM_WORLD, status, ierr )

      proc = status(MPI_SOURCE) 
      IF( datasize > 0 ) THEN
        ALLOCATE( buf(DataSize*2) )
        CALL MPI_RECV( buf, DataSize*2, MPI_INTEGER, &
           proc, 401, MPI_COMM_WORLD, status, ierr )

        DO j=1,Datasize
          gi = buf(2*j-1)
          bd = buf(2*j)
          DO k_intf=1,n_intf
            k = Ifaces(k_intf)
            IF (gdofs(k)==gi) EXIT
          END DO
          if ( k_intf<=n_intf) &
            Mesh % Faces(k) % BDOFs = MAX(Mesh % Faces(k) % BDOFs,bd)
        END DO
        DEALLOCATE(buf)
      END IF
    END DO


    DEALLOCATE( gorder, grevorder, gdofs, facen )
    IF ( ALLM ) DEALLOCATE(IsNeighbour)
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
!--------------------------------------------------------------------  
  END SUBROUTINE SParFaceNumbering
!--------------------------------------------------------------------  


!--------------------------------------------------------------------  
   SUBROUTINE SParGlobalNumbering( Mesh, OldMesh, NewNodeCnt, &
            OldIntCnts, OldIntArray, Reorder )
!-----------------------------------------------------------------------
    USE GeneralUtils
!-----------------------------------------------------------------------
     TYPE(Mesh_t) :: Mesh, OldMesh
     INTEGER, TARGET :: NewNodeCnt, OldIntArray(:), &
               OldIntCnts(:), Reorder(:)
!-----------------------------------------------------------------------
     INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
     INTEGER :: ierr

     INTEGER, POINTER :: oldnodes(:), oldnodes2(:), newnodes(:), &
          newnodes2(:), parentnodes(:,:), tosend(:), toreceive(:), &
          list1(:), list2(:), commonlist(:)
     LOGICAL :: Found, AllReceived

     TYPE Node_t
        INTEGER, POINTER :: ElementIndexes(:)
     END TYPE Node_t
     TYPE(Node_t), POINTER :: Node(:)

     TYPE Request_t
        INTEGER, POINTER :: DATA(:)
     END TYPE Request_t
     TYPE(Request_t), POINTER :: Request1(:), Request2(:)

     LOGICAL, POINTER :: IsNeighbour(:)

     TYPE(Element_t), POINTER :: Element

     INTEGER, PARAMETER :: MaxIterates = 10
     INTEGER :: i,j,k,l,m,mm,n,nn,MinProc,MaxLcl,MaxGlb,InterfaceNodes
     INTEGER :: LIndex(100), IntN, Gindex, k1, k2, Iterate, p, q
     INTEGER :: DataSize, TmpArray(3),ddd
     INTEGER, POINTER :: IntArray(:),IntCnts(:),GIndices(:),Gorder(:)
#ifdef USE_ISO_C_BINDINGS
     real(kind=dp) :: tstart
#else
     real(kind=dp) :: realtime, tstart
#endif
!-----------------------------------------------------------------------
     CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
tstart = realtime()

     Iterate = 0
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! YOU, YES YOU  DO SOMETHING ABOUT THIS 
!     The size of the biggest pack to send is:
!     3*sizeof(int)*Number of new interface nodes
     CALL CheckBuffer(10000000)
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!    find least numbered of the active procs:
!    ----------------------------------------
     MinProc = 0

     IF ( ASSOCIATED(ParEnv % IsNeighbour) ) THEN
       IsNeighbour => ParEnv % IsNeighbour
     ELSE
       ALLOCATE(IsNeighbour(ParEnv % PEs))
       IsNeighbour = .FALSE.
       DO i=1,OldMesh % Nodes % NumberOfNodes
         IF ( OldMesh % ParallelInfo % Interface(i) ) THEN
           DO j=1,SIZE(OldMesh % ParallelInfo % NeighbourList(i) % Neighbours)
             IsNeighbour(OldMesh % ParallelInfo % &
               NeighbourList(i) % Neighbours(j)+1) = .TRUE.
           END DO
         END IF
       END DO
       IsNeighbour(ParEnv % myPE+1) = .FALSE.
     END IF

!
!    Our maximum global node index:
!    -------------------------------
     MaxLcl = MAXVAL( Mesh % ParallelInfo % GlobalDOFs )
     MaxGlb = MaxLcl
     n = Mesh % Nodes % NumberOfNodes - NewNodeCnt + 1

!    Allocate space for local tables:
!    --------------------------------
     ALLOCATE( newnodes( ParEnv % PEs ), &
          newnodes2( ParEnv % PEs ), &
          oldnodes( ParEnv % PEs ), &
          oldnodes2( ParEnv % PEs ), &
          tosend( ParEnv % PEs ), &
          toreceive( ParEnv % PEs ), &
          parentnodes( Mesh % Nodes % NumberOfNodes,2 ) )
     
     newnodes    = 0
     newnodes2   = 0
     oldnodes    = 0
     oldnodes2   = 0
     tosend      = 0
     toreceive   = 0
     parentnodes = 0

     ALLOCATE( Request1(ParEnv % PEs) )
     ALLOCATE( Request2(ParEnv % PEs) )
     DO i = 1, ParEnv % PEs
        Request1(i) % DATA => NULL()
        Request2(i) % DATA => NULL()
     END DO
     !
     ! Prepare the inverse connection table for nodes and elements:
     !-------------------------------------------------------------
     ALLOCATE( Node( Mesh % Nodes % NumberOfNodes ) )
     DO i = 1,Mesh % Nodes % NumberOfNodes
        Node(i) % ElementIndexes => NULL()
     END DO
     
     DO i = 1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(i)
        DO j = 1,SIZE( Element % NodeIndexes )
           k = Element % NodeIndexes(j)
           CALL AddToCommonList( Node(k) % ElementIndexes, i )
        END DO
     END DO
     !
     ! Test write to a ep-file (please don't delete these lines, they might prove useful... :)
     !-----------------------------------------------------------------------------------------
     !PRINT *,'PE:',ParEnv % MyPE,'write ep...'
     !j = 10+ParEnv%MyPE
     !OPEN(unit=j)
     !WRITE(j,*) Mesh % Nodes % NumberOfNodes, &
     !     Mesh % NumberOfBulkElements, 1, 1, 'scalar: interface'
     !DO i = 1,Mesh % Nodes % NumberOfNodes
     !   WRITE(j,*) Mesh % Nodes % x(i), Mesh % Nodes % y(i), Mesh % Nodes % z(i)
     !END DO
     !DO i = 1, mesh % numberofbulkelements
     !   WRITE(j,*) 'body1 504', Mesh % Elements(i) % NodeIndexes-1
     !END DO
     !DO i = 1,Mesh % nodes % numberOfnodes
     !   IF( Mesh %ParallelInfo % INTERFACE(i) ) THEN
     !      WRITE(j,*) 1
     !   ELSE
     !     WRITE(j,*) 0
     !   END IF
     !END DO
     !CLOSE(j)

     !
     ! Find neighbours and parent nodes for all new interface nodes:
     ! =============================================================
     !
     ! Loop over all new nodes:
     !--------------------------
     DO i = n, Mesh % Nodes % NumberOfNodes
        IF( .NOT. Mesh % ParallelInfo % INTERFACE(i) ) CYCLE
        !
        ! This is an interface node:
        !---------------------------
        list1 => NULL() ! list of PEs related to parent node 1
        list2 => NULL() ! list of PEs related to parent node 2
        commonlist => NULL() ! intersection of list1 and list2
        !
        ! Go through all new elements connected to node i:
        !--------------------------------------------------
        DO j = 1, SIZE( Node(i) % ElementIndexes )
           Element => Mesh % Elements( Node(i) % ElementIndexes(j) )
           !
           ! Loop over all nodes of element j:
           !----------------------------------
           DO k = 1,SIZE( Element % NodeIndexes )
              !
              ! Search for parent nodes l that generated node i:
              !--------------------------------------------------
              l = Element % NodeIndexes(k)
              IF( l >= n ) CYCLE ! parents have local number < n
              IF( .NOT. Mesh % ParallelInfo % Interface(l) ) CYCLE
              IF( ANY( parentnodes(i,:)==l) ) CYCLE ! already found
              !
              ! Construct the parent table:
              !----------------------------
              IF( parentnodes(i,1)==0 ) THEN
                 parentnodes(i,1) = l 
                 ! This is the list of PEs sharing parent node 1:
                 list1 => Mesh % ParallelInfo % NeighbourList(l) % Neighbours
              ELSE
                 parentnodes(i,2) = l 
                 ! This is the list of PEs sharing parent node 2:
                 list2 => Mesh % ParallelInfo % NeighbourList(l) % Neighbours
              END IF         
           END DO
        END DO
        !
        ! Determine the intersection of the PE-lists:
        !--------------------------------------------
        IF( ASSOCIATED(list1) .AND. ASSOCIATED(list2) ) THEN
           IF(ASSOCIATED(commonlist)) DEALLOCATE(commonlist)
           commonlist => NULL()
           DO p = 1,SIZE(list1)
              DO q = 1,SIZE(list2)
                 IF( list1(p) == list2(q)) &
                      CALL AddToCommonList( commonlist, list1(p) ) 
              END DO
           END DO
           !
           ! Now, we should have a list of PEs common to the parents:
           !----------------------------------------------------------
           IF( .NOT.ASSOCIATED(commonlist) ) CYCLE
           !
           ! If everything went ok and the mesh data was unique, there
           ! are only two PEs common to both parents. If not, thats bad:
           !------------------------------------------------------------
           !IF( SIZE(commonlist) > 2 ) THEN
           !   WRITE(*,'(A,I4,A,I6,A,2I6,A,2I6,A,10I4)') &
           !        'SParGlobalNumbering: PE:', ParEnv % MyPE+1, &
           !        ' Unable to determine owner for node(loc):', i, &
           !        ' Data not unique. Masters(loc):', parentnodes(i,:), &
           !        ' (glob):', Mesh % ParallelInfo % GlobalDOFs(parentnodes(i,:)), &
           !        ' Commonlist:', commonlist+1
           !END IF
           !
           ! Node i is given to the smallest numbered PE common to the parents:
           !-------------------------------------------------------------------
           ALLOCATE( gindices( SIZE(commonlist)))
           gindices = commonlist
           DEALLOCATE( commonlist )
           CALL SORT( SIZE(gindices), gindices )
           IF( ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) &
              DEALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
           Mesh % ParallelInfo % NeighbourList(i) % Neighbours => Gindices
        ELSE
           !
           ! Either list1 or list2 is empty. Thats really bad:
           !---------------------------------------------------
           WRITE(*,'(A,I4,A,I6)') 'SParIterGlobalNumbering: PE:', ParEnv % MyPE+1, &
                ' Could not determine owner for node(loc)=', i
           CALL MPI_FINALIZE( ierr )
        END IF
        !
        ! Finalize by sorting the parent table:
        !---------------------------------------
        CALL Sort( 2, parentnodes(i,:) )
     END DO

!************************************************************************************
10   CONTINUE
     !
     ! Iteratively update and correct the neighbour-lists and number the new nodes:
     !==============================================================================
     !
     ! Begin iteration (if the mesh was good, then only one round is needed):
     !-----------------------------------------------------------------------
     DO i = 1, ParEnv % PEs
        IF( ASSOCIATED(Request1(i) % DATA) ) DEALLOCATE( Request1(i) % DATA )
        IF( ASSOCIATED(Request2(i) % DATA) ) DEALLOCATE( Request2(i) % DATA )
        Request1(i) % DATA => NULL() ! buffer for data exchange
        Request2(i) % DATA => NULL() ! buffer for data exchange
     END DO
     !
     ! Count the current situation in what comes to nodes owned by us:
     !----------------------------------------------------------------
     oldnodes = 0
     newnodes = 0
     j = ParEnv % MyPE
     DO i = 1, Mesh % Nodes % NumberOfNodes
        k = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1)
        IF( k /= j ) CYCLE
        IF( Mesh % ParallelInfo % GlobalDOFs(i)  > 0 ) THEN
           oldnodes(j+1) = oldnodes(j+1)+1
        ELSE
           newnodes(j+1) = newnodes(j+1)+1
        END IF
     END DO
     !
     ! Distribute the knowledge about ownerships (here, assumig all active):
     !----------------------------------------------------------------------
     oldnodes2 = 0
     CALL MPI_ALLREDUCE( oldnodes, oldnodes2, ParEnv % PEs, &   ! <- fix 
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

     newnodes2 = 0
     CALL MPI_ALLREDUCE( newnodes, newnodes2, ParEnv % PEs, &   ! <- fix 
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
     !
     ! Number the new nodes:
     !----------------------
     j = ParEnv % MyPE
     ! Start numbering from index k:
     k = SUM( oldnodes2 ) + SUM( newnodes2(1:j) ) + 1
     DO i = 1, Mesh % Nodes % NumberOfNodes
        l = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1)
        IF( l /= j ) CYCLE
        IF( Mesh % ParallelInfo % GlobalDOFs(i) == 0 ) THEN
           Mesh % ParallelInfo % GlobalDOFs(i) = k
           k = k+1
        END IF
     END DO
     !
     ! Some new interface nodes may still lack a global number, if
     ! we are not the owner of these nodes. Hence, we'll have to make
     ! a query to all neighbour PEs to find these numbers.
     !===============================================================
     !
     ! Determine the size of data to exchange between neighbours:
     !-----------------------------------------------------------
     tosend = 0
     toreceive = 0
     DO i = n, Mesh % Nodes % NumberOfNodes
        IF( Mesh % ParallelInfo % Interface(i) ) THEN
           j = Mesh % ParallelInfo % Neighbourlist(i) % Neighbours(1)
           IF( j /= ParEnv % MyPE ) THEN
              toreceive(j+1) = toreceive(j+1)+1
           ELSE
              DO j = 2,SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
                 k = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(j)
                 tosend(k+1) = tosend(k+1)+1
              END DO
           END IF
        END IF
     END DO
     !
     ! Distribute the interface data:
     !--------------------------------
     DO i = 1, ParEnv % PEs
        IF( i-1 == ParEnv % MyPE ) CYCLE
        IF( .NOT. IsNeighbour(i) ) CYCLE

        CALL MPI_BSEND( tosend(i), 1, MPI_INTEGER, i-1, 100, MPI_COMM_WORLD, ierr)
        IF( tosend(i) < 1 ) CYCLE

        DataSize = 3*tosend(i)
        ALLOCATE( gindices(DataSize) )
        k = 1
        DO l = n, Mesh % Nodes % NumberOfNodes
           IF( .NOT.( Mesh % ParallelInfo % Interface(l) ) ) CYCLE
           m = Mesh % ParallelInfo % NeighbourList(l) % Neighbours(1)
           IF( m /= ParEnv % MyPE ) CYCLE

           DO m = 2,SIZE( Mesh % ParallelInfo % NeighbourList(l) % Neighbours )
              mm = Mesh % ParallelInfo % NeighbourList(l) % Neighbours(m)
              IF( mm+1 .NE. i ) CYCLE
              gindices(k) = Mesh % ParallelInfo % GlobalDOFs(l)
              gindices(k+1) = Mesh % ParallelInfo % GlobalDOFs(parentnodes(l,1))
              gindices(k+2) = Mesh % ParallelInfo % GlobalDOFs(parentnodes(l,2))
              k = k+3
           END DO
        END DO

        CALL MPI_BSEND( gindices, DataSize, MPI_INTEGER, i-1, 400, MPI_COMM_WORLD, ierr )
        
        DEALLOCATE( gindices )
     END DO
     !
     ! Retrieve missing interface node numbers from our neighbours:
     !--------------------------------------------------------------
     DO i = 1, ParEnv % PEs
        IF( i-1 == ParEnv % MyPE ) CYCLE
        IF( .NOT. IsNeighbour(i) ) CYCLE

        CALL MPI_RECV( DataSize, 1, MPI_INTEGER, i-1, 100, MPI_COMM_WORLD, status, ierr )
        !IF( DataSize /= toreceive(i) ) THEN
        !   WRITE(*,'(A,I4,A,I4,A,I5,A,I5)') 'SParGlobalNumbering: PE:',ParEnv % MyPE+1, &
        !        ' Disagreement about data size with PE:', i, ' Expecting:',toreceive(i), &
        !        ' Actually got:', DataSize
        !END IF
        IF( DataSize < 1 ) CYCLE

        ALLOCATE( gindices(3*DataSize) ) ! work space
        CALL MPI_RECV( gindices, 3*DataSize, MPI_INTEGER, i-1, 400, MPI_COMM_WORLD, status, ierr )

        DO k = n, Mesh % Nodes % NumberOfnodes
           IF( .NOT. Mesh % ParallelInfo % Interface(k) ) CYCLE
           IF( Mesh % ParallelInfo % GlobalDOFs(k) > 0 ) CYCLE
           
           DO l = 1, DataSize
              m = 3*l-2
              IF( Mesh % ParallelInfo % GlobalDOFs(parentnodes(k,1)) == gindices(m+1) .AND. &
                  Mesh % ParallelInfo % GlobalDOFs(parentnodes(k,2)) == gindices(m+2) .OR. &
                  Mesh % ParallelInfo % GlobalDOFs(parentnodes(k,1)) == gindices(m+2) .AND. &
                  Mesh % ParallelInfo % GlobalDOFs(parentnodes(k,2)) == gindices(m+1) ) THEN
                 Mesh % ParallelInfo % GlobalDOFs(k) = gindices(m) ! global DOF number
                 gindices(m) = 0 ! mark used
                 EXIT
              END IF
           END DO
        END DO
        !
        ! PROBLEM CASE 1: Received too much data
        ! =======================================
        ! It is possible that we actually got too much data, because one or more PEs
        ! ennareously think of us as their neighbour. In this case, we will have to
        ! inform these PEs about the situation and ask them to remove us from their
        ! neighbour lists:
        !----------------------------------------------------------------------------
        DO l=1,DataSize
           IF( gindices(3*l-2)>0 ) THEN
              !WRITE(*,'(A,I5,A,I8,A,2I8)') 'SParGlobalNumbering: PE',ParEnv % MyPE+1, &
              !     ' Got unexpected global data for node(glob):',gindices(3*l-2),&
              !     ' Parent nodes(glob):', gindices(3*l-1), gindices(3*l)
              !PRINT *,'Sending a request to',i,'to remove me from the neighnour list.'
              !
              ! Request PE=i to remove us from the neighbour list:
              !----------------------------------------------------
              CALL AddToCommonList( Request1(i) % DATA, gindices(3*l-2) )  ! i-node
              CALL AddToCommonList( Request1(i) % DATA, gindices(3*l-1) )  ! parent 1
              CALL AddToCommonList( Request1(i) % DATA, gindices(3*l-0) )  ! parent 2
           END IF
        END DO
        DEALLOCATE( gindices )
     END DO
     !
     ! PROBLEM CASE 2: Data not sufficient
     ! ===================================
     ! We should have received unique global dof-numbers from the owners of the
     ! corresponding interface nodes. If we did not, our own neighbour list is
     ! too long and we have a false idea about the owner. In this case, we will
     ! have to ask all possible neighbours to check if they have the node, and if
     ! they don't, ask them to return a request to be removed from our neighbour list:
     !---------------------------------------------------------------------------------
     DO i = n, Mesh % NumberOfNodes
       IF( Mesh % ParallelInfo % GlobalDOFs(i) < 1 .OR. &
            SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours )>2 ) THEN

          !WRITE(*,'(A,I5,A,I8,A,2I8)') 'SParGlobalNumbering: PE',ParEnv % MyPE+1, &
          !     ' Did not receive global number for node(loc):', i, ' Parents nodes(glob):', &
          !     Mesh % ParallelInfo % GlobalDOFs( parentnodes(i,:) )
          !PRINT *,'Have to inform the following PEs:', &
          !     Mesh % ParallelInfo % NeighbourList(i) % Neighbours +1 
          DO j = 1,SIZE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
             k = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(j)
              IF( k == ParEnv % MyPE ) CYCLE
             TmpArray(1) = Mesh % ParallelInfo % GlobalDOFs( i )
             TmpArray(2) = Mesh % ParallelInfo % GlobalDOFs( parentnodes(i,1) )
             TmpArray(3) = Mesh % ParallelInfo % GlobalDOFs( parentnodes(i,2) )
             CALL AddToCommonList( Request2(k+1) % DATA, TmpArray(1) )
             CALL AddToCommonList( Request2(k+1) % DATA, TmpArray(2) )
             CALL AddToCommonList( Request2(k+1) % DATA, TmpArray(3) )
          END DO
        END IF
     END DO
     !
     ! Send the query (CASE 2):
     !--------------------------
!ddd = 0
     DO i = 1, ParEnv % PEs              
        IF( i-1 == ParEnv % MyPE ) CYCLE
        IF( .NOT. IsNeighbour(i) ) CYCLE

        DataSize = 0
        IF( ASSOCIATED( Request2(i) % DATA ) ) DataSize = SIZE( Request2(i) % DATA )
!ddd = ddd+DataSize

        !CALL MPI_SEND( DataSize, 1, MPI_INTEGER, i-1, 900, MPI_COMM_WORLD, ierr)
        CALL MPI_BSEND( DataSize, 1, MPI_INTEGER, i-1, 900, MPI_COMM_WORLD, ierr)
        IF( DataSize < 1 ) CYCLE

        CALL MPI_BSEND( Request2(i) % DATA, DataSize, MPI_INTEGER, i-1, 950, MPI_COMM_WORLD, ierr)
     END DO
!print *,'*****',ParEnv % MyPE+1,ddd
     !
     ! Receive the query (CASE 2):
     !----------------------------
     DO i = 1, ParEnv % PEs
        IF( i-1 == ParEnv % MyPE ) CYCLE
        IF( .NOT. IsNeighbour(i) ) CYCLE

        CALL MPI_RECV( DataSize, 1, MPI_INTEGER, i-1, 900, MPI_COMM_WORLD, status, ierr )
        IF( DataSize < 1 ) CYCLE
        
        ALLOCATE( IntArray(DataSize) )
        CALL MPI_RECV( IntArray, DataSize, MPI_INTEGER, i-1, 950, MPI_COMM_WORLD, status, ierr )

        DO q = 1, DataSize/3
           TmpArray = IntArray(3*q-2:3*q)
           !
           !PRINT *,'SParGlobalNumbering: PE',ParEnv % MyPE+1,' received a query from PE',i, &
           !     ' which asks if Im the owner of parents(glob):', TmpArray(2), TmpArray(3)
           !
           ! Ok, lets see if we have the parents requested:
           !-----------------------------------------------
           mm = 0
           nn = 0
           DO j = n, Mesh % NumberOfNodes
              IF( .NOT. Mesh % ParallelInfo % Interface(j)) CYCLE
              mm = 0
              nn = 0
              DO k = 1,SIZE( Node(j) % ElementIndexes )
                 l = Node(j) % ElementIndexes(k)
                 Element => Mesh % Elements(l)
                 IF( ANY(Mesh % ParallelInfo % GlobalDOFs(Element % NodeIndexes) == TmpArray(2)) ) mm = 1
                 IF( ANY(Mesh % ParallelInfo % GlobalDOFs(Element % NodeIndexes) == TmpArray(3)) ) nn = 1
              END DO
              IF( mm==1 .AND. nn==1 ) THEN
                 !PRINT *,'Ok, ive got it. node(loc):',j,'(glob):',mesh % parallelInfo % GlobalDOFs(j)
                 EXIT
              END IF
           END DO
           
           IF( .NOT.(mm==1 .AND. nn==1) ) THEN
              !PRINT *,'Nope, I dont have it. Sending back request to remove me from the neighbour list.'
              CALL AddToCommonList( Request1(i) % DATA, TmpArray(1) )
              CALL AddToCommonList( Request1(i) % DATA, TmpArray(2) )
              CALL AddToCommonList( Request1(i) % DATA, TmpArray(3) )
           END IF
        END DO
        DEALLOCATE( IntArray )
     END DO
     !
     ! Send the removal requests to all neighbours (CASE 1):
     !------------------------------------------------------
     DO i = 1, ParEnv % PEs              
        IF( i-1 == ParEnv % MyPE ) CYCLE
        IF( .NOT. IsNeighbour(i) ) CYCLE

        DataSize = 0
        IF( ASSOCIATED( Request1(i) % DATA ) ) DataSize = SIZE( Request1(i) % DATA )
        CALL MPI_BSEND( DataSize, 1, MPI_INTEGER, i-1, 700, MPI_COMM_WORLD, ierr)
        IF( DataSize < 1 ) CYCLE
        CALL MPI_BSEND( Request1(i) % DATA, DataSize, MPI_INTEGER, i-1, 1800, MPI_COMM_WORLD, ierr )
     END DO
     !
     ! Receive the removal requests and update the neighbour lists (CASE 1):
     !----------------------------------------------------------------------
     Gorder=>Null()
     DO i = 1, ParEnv % PEs
        IF( i-1 == ParEnv % MyPE ) CYCLE
        IF( .NOT. IsNeighbour(i) ) CYCLE

        CALL MPI_RECV( DataSize, 1, MPI_INTEGER, i-1, 700, MPI_COMM_WORLD, status, ierr )
        IF( DataSize < 1 ) CYCLE
        
        ALLOCATE( IntArray(DataSize) )
        CALL MPI_RECV( IntArray, DataSize, MPI_INTEGER, i-1, 1800, MPI_COMM_WORLD, status, ierr )
        
        mm = DataSize/3
        DO nn = 1,mm
           TmpArray = IntArray(3*nn-2:3*nn)
           !PRINT *,'SParGlobalNumbering: PE',ParEnv % MyPE+1,' got a message from PE',i, &
           !     ' which asks to be removed from the neighbour list of node(glob):', &
           !     TmpArray(1), 'Parents(glob):', TmpArray(2), TmpArray(3)
           DO j = n, Mesh % NumberOfNodes
              IF( Mesh % ParallelInfo % GlobalDOFs(j) == TmpArray(1)) THEN
                 ! Check also for parents to avoid screwing up anything:
                 IF( .NOT.( mesh % parallelinfo % globaldofs(parentnodes(j,1)) == TmpArray(2) .AND. &
                      mesh % parallelinfo % globaldofs(parentnodes(j,2)) == TmpArray(3) ) ) CYCLE
                 !PRINT *,'Ok, found it. Local number is', j
                 !PRINT *,'Neighbours were:', Mesh % ParallelInfo % NeighbourList(j) % Neighbours+1
                 !
                 ! Removing:
                 !----------
                 IF( ASSOCIATED(gorder)) DEALLOCATE( gorder )
                 gorder => NULL()
                 DO k = 1,SIZE( Mesh % ParallelInfo % NeighbourList(j) % Neighbours )
                    l = Mesh % ParallelInfo % NeighbourList(j) % Neighbours(k)
                    IF( l == i-1 ) CYCLE
                    CALL AddToCommonList( gorder, l )
                 END DO
                 DEALLOCATE( Mesh % ParallelInfo % NeighbourList(j) % Neighbours )
                 ALLOCATE( IntCnts( SIZE(gorder) ) )
                 IntCnts = gorder
                 IF (ASSOCIATED(gorder)) DEALLOCATE( gorder )
                 Mesh % ParallelInfo % NeighbourList(j) % Neighbours => IntCnts
                 !PRINT *,'Updated neighbours:', Mesh % ParallelInfo % NeighbourList(j) % Neighbours+1
                 EXIT
              END IF
           END DO
        END DO
        IF(ASSOCIATED(Intarray)) DEALLOCATE( IntArray )
     END DO     
     !
     ! Make sure that all nodes have a global number:
     !-----------------------------------------------
     Found = .FALSE.
     DO i = 1,Mesh % NumberOfNodes
        IF( Mesh % ParallelInfo % GlobalDOFs(i) < 1 ) THEN
           Found = .TRUE.
           !WRITE(*,'(A,I4)') 'SParIterGlobalNumbering: PE:', ParEnv % MyPE+1, &
           !     ' No global number for node(loc):',i
           IF( ASSOCIATED( mesh % parallelinfo % neighbourlist(i) % neighbours ) ) THEN
              !PRINT *, 'Owners:', Mesh % ParallelInfo % NeighbourList(i) % Neighbours+1
           ELSE
              CALL AddToCommonList( Mesh % ParallelInfo % NeighbourList(i) % Neighbours, ParEnv % MyPE )
           END IF
        END IF

        ! Check also the neighbour lists:
        !--------------------------------
        IF( ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) THEN
           IF( SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours) < 1 ) THEN
              PRINT *,'PE, node: ***** This node is missing the owner ****',ParEnv % MyPE+1, i
           END IF
           
           IF (ANY(Mesh % ParallelInfo % NeighbourList(i) % Neighbours < 0 ) .OR. &
                ANY(Mesh % ParallelInfo % NeighbourList(i) % Neighbours > ParEnv % PEs-1 ) ) THEN
              PRINT *,'PE, node: ***** This node has a bad owner ****',ParEnv % MyPE+1, i
              PRINT *, Mesh % ParallelInfo % NeighbourList(i) % Neighbours
              Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) = ParEnv % MyPE              
           END IF
        ELSE
           PRINT *,'PE, node: ***** This node is missing the owner ****',ParEnv % MyPE+1, i
        END IF
     END DO

     oldnodes = 0 ! workspace
     IF( Found ) oldnodes( ParEnv % MyPE+1 ) = 1
     oldnodes2 = 0
     CALL MPI_ALLREDUCE( oldnodes, oldnodes2, ParEnv % PEs, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
     j = SUM( oldnodes2 )

     IF( j < 1 ) GOTO 20
     !IF( ParEnv % MyPE==0 ) &
     !     write(*,'(A,I6,A)') 'SParIterGlobalNumbering:  Numbering not consistent, re-iterating...'
     Iterate = Iterate+1
     IF(Iterate > MaxIterates ) THEN
        WRITE(*,'(A,I6,A)') 'SParIterGlobalNumbering: PE: ', ParEnv % MyPE+1,'Max iterations exceeded'
        CALL MPI_FINALIZE( MPI_COMM_WORLD, ierr )
     END IF
     DO i = n, Mesh % NumberOfNodes
        Mesh % ParallelInfo % GlobalDOFs(i) = 0
     END DO
     GOTO 10

20   CONTINUE
     IF( ParEnv % MyPE==0 ) &
          WRITE(*,'(A,I6,A)') 'SParIterGlobalNumbering:  Numbering seems consistent, advancing...'

     ! Deallocate temp arrays:
     !------------------------
     DEALLOCATE( oldnodes, oldnodes2, newnodes, newnodes2, &
          parentnodes, tosend, toreceive )
     
     DO i = 1,Mesh % Nodes % NumberOfNodes
       DEALLOCATE( Node(i) % ElementIndexes )
     END DO
     DEALLOCATE(Node)

     DO i = 1, ParEnv % PEs
       IF( ASSOCIATED(Request1(i) % DATA) ) DEALLOCATE(Request1(i) % DATA)
       IF( ASSOCIATED(Request2(i) % DATA) ) DEALLOCATE(Request2(i) % DATA)
     END DO
     DEALLOCATE( Request1, Request2 )

     !-------------------------------------------------------------------
     ! Now, the numbering should be unique. It remains to sort the nodes.
     !-------------------------------------------------------------------
!
!    Sort our own nodes according to ascending
!    global order:
!    -----------------------------------------
     ALLOCATE( IntCnts(NewNodeCnt), Gorder(NewNodeCnt) )
     Gorder = (/ (i,i=1,NewNodeCnt) /)

     CALL SortI( NewNodeCnt, Mesh % ParallelInfo % GlobalDOFs(n:), Gorder )

     DO i=1,NewNodeCnt
        IntCnts(Gorder(i)) = i
     END DO
!
!    Reorder will return the nodal reordering
!    to the caller:
!    ----------------------------------------
     Reorder(n:) = IntCnts + n - 1
!
!    Order the whole of the nodal structure
!    according to the changed order of the 
!    global numbers:
!    --------------------------------------
     DO i=1,NewNodeCnt
        k = Gorder(i)
        CALL SwapNodes( Mesh, i+n-1, k+n-1 )

        j = IntCnts(i)
        IntCnts(i) = IntCnts(k)
        IntCnts(k) = j

        Gorder(IntCnts(k)) = k
     END DO

     DEALLOCATE( IntCnts )
     IF ( .NOT. ASSOCIATED(Parenv % IsNeighbour,IsNeighbour) ) &
       DEALLOCATE(isNeighbour)

     CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
!print*,'globnum: ', parenv % mype, RealTime()-tstart; call flush(6)

     RETURN

!???????????????????????????????????????????????????????????????????????????????

!
!    Lowest numbered PE will compute the size of
!    the (old) global node array, and eventually
!    distribute the knowledge:
!
!    Others send their maxlcl to minproc
!    -------------------------------------------
     IF ( ParEnv % MyPE == MinProc-1 ) THEN
        j = 1 !?????
        DO i = MinProc+1, ParEnv % PEs
           IF ( ParEnv % Active(i) ) THEN
              CALL MPI_RECV( k, 1, MPI_INTEGER, i-1, &
                   10, MPI_COMM_WORLD, status, ierr )
              MaxGlb = MAX( MaxGlb, k )
           END IF
        END DO
     ELSE
        CALL MPI_BSEND( MaxLcl, 1, MPI_INTEGER, &
             MinProc-1, 10, MPI_COMM_WORLD, ierr )
     END IF

!print *,'(2) Id:',ParEnv % MyPE+1,'MaxGlb=',MaxGlb

!
!    Recieve new interface nodes from lower
!    numbered PEs, and check if they are
!    relevant to us:
!    ---------------------------------------

     DO i = MinProc, ParEnv % MyPE

!if( ParEnv % IsNeighbour(i) ) &
!     print *,'(3) Id:',ParEnv % MyPE+1,'Neighbour=',i

        IF ( .NOT. (ParEnv % Active(i) .AND. ParEnv % Isneighbour(i)) ) CYCLE
        
        CALL MPI_RECV( InterfaceNodes, 1, MPI_INTEGER, &
             i-1, 14, MPI_COMM_WORLD, status, ierr )

!print *,'(4) Id:',ParEnv % MyPE+1,'Neighbour=',i,'InterfaceNodes=',InterfaceNodes
        
        IF ( InterfaceNodes > 0 ) THEN
           ALLOCATE( GIndices(InterfaceNodes), IntCnts(InterfaceNodes) )
           
           CALL MPI_RECV( GIndices, InterfaceNodes, MPI_INTEGER, &
                i-1, 15, MPI_COMM_WORLD, status, ierr )

!print *,'(5) Id:', ParEnv % MyPE+1, 'Neighbour=',i,'GIndices=',GIndices
           
           CALL MPI_RECV( IntCnts, InterfaceNodes, MPI_INTEGER, &
                i-1, 16, MPI_COMM_WORLD, status, ierr )

!print *,'(6) Id:', ParEnv % MyPE+1, 'Neighbour=',i,'IntCnts=',IntCnts
           
           CALL MPI_RECV( k, 1, MPI_INTEGER, &
                i-1, 17, MPI_COMM_WORLD, status, ierr )

!print *,'(7) Id:', ParEnv % MyPE+1, 'Neighbour=',i,'k=',k
           
           ALLOCATE( IntArray(k) )
           
           CALL MPI_RECV( IntArray, SIZE(IntArray), MPI_INTEGER, &
                i-1, 18, MPI_COMM_WORLD, status, ierr )

!print *,'(8) Id:', ParEnv % MyPE+1, 'Neighbour=',i,'IntArray=',IntArray

!
!          Update our view of the global numbering
!          at the interface nodes:
!          ---------------------------------------
           l = 0
           DO j=1,InterfaceNodes
              Lindex = 0
              IntN = IntCnts(j)
              DO k=1,IntN
                 Lindex(k) = SearchNode( Mesh % ParallelInfo, IntArray(l+k), 1, n-1 )
              END DO

              IF ( ALL( Lindex(1:IntN) > 0 ) ) THEN
!
!                This node belongs to us as well well:
!                -------------------------------------
                 k2 = 0
                 k1 = 0
                 DO k=n,Mesh % Nodes % NumberOfNodes
                    IF ( .NOT.Mesh % ParallelInfo % INTERFACE(k) ) CYCLE

                    k1 = k1 + 1
                    IF ( IntN == OldIntCnts(k1) ) THEN
                       IF ( ALL( IntArray(l+1:l+IntN) == &
                                    OldIntArray(k2+1:k2+IntN)) ) THEN
                           Mesh % ParallelInfo % GlobalDOFs(k) = GIndices(j)
                           EXIT
                        END IF
                     END IF
                     k2 = k2 + OldIntCnts(k1)
                 END DO
              END IF
              l = l + IntN
           END DO
           DEALLOCATE( Gindices, IntCnts, IntArray )
        END IF
     END DO
!
!    Update the current numbering from the 
!    previous PE in line, this will make the
!    execution strictly serial. Maybe there
!    would be an easier way?
!    ---------------------------------------
     IF ( ParEnv % MyPE > MinProc-1 ) THEN
        DO i=ParEnv % MyPE, MinProc, -1
           IF ( ParEnv % Active(i) ) THEN
              CALL MPI_RECV( MaxGlb, 1, MPI_INTEGER, &
               i-1, 20, MPI_COMM_WORLD, status, ierr )

!print *,'(20) Id:',ParEnv % MyPE+1,'dest=',i,'MaxGlb=',MaxGlb

              EXIT
           END IF
        END DO
     END IF
!
!    Renumber our own new set of nodes:
!    ----------------------------------
     DO i=n,Mesh % Nodes % NumberOfNodes
        IF ( Mesh % ParallelInfo % GlobalDOFs(i) == 0 ) THEN
           MaxGlb = MaxGlb + 1
           Mesh % ParallelInfo % GlobalDOFs(i) = MaxGlb
        END IF
     END DO

!
!    Extract interface nodes:
!    ------------------------
     InterfaceNodes = COUNT( Mesh % ParallelInfo % INTERFACE(n:) )
     IF ( InterfaceNodes > 0 ) ALLOCATE( Gindices(InterfaceNodes) )

     InterfaceNodes = 0
     DO i=n,Mesh % Nodes % NumberOfNodes
        IF ( Mesh % ParallelInfo % INTERFACE(i) ) THEN
           InterfaceNodes = InterfaceNodes + 1
           Gindices(InterfaceNodes) = Mesh % ParallelInfo % GlobalDOFs(i)
        END IF
     END DO
!
!    Send new interface nodes to higher numbered PEs:
!    ------------------------------------------------
     DO i=ParEnv % MyPE+2,ParEnv % PEs
        IF ( .NOT. (ParEnv % Active(i) .AND. ParEnv % Isneighbour(i)) ) CYCLE

        CALL MPI_BSEND( InterfaceNodes, 1, MPI_INTEGER, &
               i-1, 14, MPI_COMM_WORLD, ierr )

        IF ( InterfaceNodes > 0 ) THEN
           CALL MPI_BSEND( GIndices, InterfaceNodes, MPI_INTEGER, &
                    i-1, 15, MPI_COMM_WORLD, ierr )

           CALL MPI_BSEND( OldIntCnts, InterfaceNodes, &
              MPI_INTEGER, i-1, 16, MPI_COMM_WORLD, ierr )

           j = SIZE(OldIntArray)
           CALL MPI_BSEND( j, 1, &
              MPI_INTEGER, i-1, 17, MPI_COMM_WORLD, ierr )

           CALL MPI_BSEND( OldIntArray, SIZE(OldIntArray), &
              MPI_INTEGER, i-1, 18, MPI_COMM_WORLD, ierr )

        END IF
     END DO

     DEALLOCATE( GIndices )
!
!    Send go singal to next PE in line...
!    ------------------------------------
     DO i = ParEnv % MyPE+2, ParEnv % PEs
        IF ( ParEnv % Active(i) ) THEN
           CALL MPI_BSEND( MaxGlb, 1, MPI_INTEGER, &
               i-1, 20, MPI_COMM_WORLD, ierr )
           EXIT
        END IF
     END DO
!
!    Sort our own nodes according to ascending
!    global order:
!    -----------------------------------------
     ALLOCATE( IntCnts(NewNodeCnt), Gorder(NewNodeCnt) )
     DO i=1,NewNodeCnt
        Gorder(i) = i
     END DO

     CALL SortI( NewNodeCnt, Mesh % ParallelInfo % GlobalDOFs(n:), Gorder )

     DO i=1,NewNodeCnt
        IntCnts(Gorder(i)) = i
     END DO
!
!    Reorder will return the nodal reordering
!    to the caller:
!    ----------------------------------------
     Reorder(n:) = IntCnts + n - 1
!
!    Order the whole of the nodal structure
!    according to the changed order of the 
!    global numbers:
!    --------------------------------------
     DO i=1,NewNodeCnt
        k = Gorder(i)
        CALL SwapNodes( Mesh, i+n-1, k+n-1 )

        j = IntCnts(i)
        IntCnts(i) = IntCnts(k)
        IntCnts(k) = j

        Gorder(IntCnts(k)) = k
     END DO

     DEALLOCATE( IntCnts )
!
!    Ok, now we have generated the global numbering
!    of nodes. We still need to distribute the
!    information which PEs share which of the new
!    interface nodes:
!    -----------------------------------------------
     InterfaceNodes = COUNT( Mesh % ParallelInfo % INTERFACE(n:) )
     ALLOCATE( GIndices( InterfaceNodes ) )
     j = 0
     DO i=n,Mesh % Nodes % NumberOfNodes
        IF ( Mesh % ParallelInfo % INTERFACE(i) ) THEN
           j = j + 1
           GIndices(j) = Mesh % ParallelInfo % GlobalDOFs(i)
           ALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours(ParEnv % PEs) )
           Mesh % ParallelInfo % NeighbourList(i) % Neighbours = -1
        END IF
     END DO

     DO i=MinProc,ParEnv % PEs
        IF ( ParEnv % MyPE == i-1 ) CYCLE
        IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
           CALL MPI_BSEND( InterfaceNodes, 1, &
              MPI_INTEGER, i-1, 30, MPI_COMM_WORLD, ierr )

           CALL MPI_BSEND( GIndices, InterfaceNodes, &
              MPI_INTEGER, i-1, 31, MPI_COMM_WORLD, ierr )
        END IF
     END DO

     DEALLOCATE( Gindices )

     ALLOCATE( IntCnts( Mesh % Nodes % NumberOfNodes ) )

     IntCnts = 0
     DO i=n,Mesh % Nodes % NumberOfNodes
        IF ( Mesh % ParallelInfo % INTERFACE(i) ) THEN
           IntCnts(i) = IntCnts(i) + 1
           Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) = ParEnv % MyPE
        END IF
     END DO

     DO i=MinProc,ParEnv % PEs
        IF ( ParEnv % MyPE == i-1 ) CYCLE
        IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
           CALL MPI_RECV( InterfaceNodes, 1, MPI_INTEGER, &
               i-1, 30, MPI_COMM_WORLD, status, ierr )

           ALLOCATE( GIndices(InterfaceNodes) )

           CALL MPI_RECV( GIndices, InterfaceNodes, MPI_INTEGER, &
               i-1, 31, MPI_COMM_WORLD, status, ierr )

           DO j=1,InterfaceNodes
              k = SearchNode( Mesh % ParallelInfo, Gindices(j), n )
              IF ( k <= 0 ) CYCLE
              IntCnts(k) = IntCnts(k) + 1
              Mesh % ParallelInfo % NeighbourList(k) % Neighbours(IntCnts(k)) = i-1
           END DO

           DEALLOCATE( GIndices )
        END IF
     END DO
!
!    Reallocate the nodal neighbour lists to
!    correct sizes:
!    ---------------------------------------
     DO i=n,Mesh % Nodes % NumberOfNodes
        IF ( Mesh % ParallelInfo % INTERFACE(i) ) THEN
           k = IntCnts(i)
           ALLOCATE( Gindices(k) ) ! just work space
           Gindices = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1:k)
           CALL Sort( k, Gindices )
           DEALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
           Mesh % ParallelInfo % NeighbourList(i) % Neighbours => Gindices
        END IF
     END DO
     
     DEALLOCATE( IntCnts )

! final test:
!     IF( ParEnv % MyPE == 3 ) THEN
!        DO i = 1, Mesh % Nodes % NumberOfNodes
!           PRINT *,'Local:',i, &
!                'Global:' ,Mesh % Parallelinfo % GlobalDOFs(i), &
!                'Interface:', Mesh % ParallelInfo % INTERFACE(i), &
!                'Neighbours:', Mesh % ParallelInfo % NeighbourList(i) % Neighbours + 1
!        END DO
!     END IF

PRINT *,'****OK:', parenv % mype+1
DO i = 1, mesh % nodes % numberofnodes
   PRINT *,'(+++)',parenv % mype+1, i, Mesh % ParallelInfo % GlobalDOFs(i), &
        Mesh % ParallelInfo % INTERFACE(i), Mesh % ParallelInfo % NeighbourList(i) % Neighbours
END DO


     CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

CONTAINS

!-----------------------------------------------------------------------
     SUBROUTINE SwapNodes( Mesh, i, k )
!-----------------------------------------------------------------------
        INTEGER :: i,k
        TYPE(Mesh_t) :: Mesh
!-----------------------------------------------------------------------
        REAL(KIND=dp) :: swapx,swapy,swapz
        LOGICAL :: swapi
        INTEGER, POINTER :: swapl(:)
!-----------------------------------------------------------------------
        swapx =  Mesh % Nodes % x(i)
        swapy =  Mesh % Nodes % y(i)
        swapz =  Mesh % Nodes % z(i)
        swapi =  Mesh % ParallelInfo % INTERFACE(i)
        swapl => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
 
        Mesh % Nodes % x(i) = Mesh % Nodes % x(k)
        Mesh % Nodes % y(i) = Mesh % Nodes % y(k)
        Mesh % Nodes % z(i) = Mesh % Nodes % z(k)
        Mesh % ParallelInfo % INTERFACE(i) = Mesh % ParallelInfo % INTERFACE(k) 
        Mesh % ParallelInfo % NeighbourList(i) % Neighbours => &
                 Mesh % ParallelInfo % NeighbourList(k) % Neighbours

        Mesh % Nodes % x(k) = swapx
        Mesh % Nodes % y(k) = swapy
        Mesh % Nodes % z(k) = swapz
        Mesh % ParallelInfo % INTERFACE(k) = swapi
        Mesh % ParallelInfo % NeighbourList(k) % Neighbours => swapl
!-----------------------------------------------------------------------
     END SUBROUTINE SwapNodes
!-----------------------------------------------------------------------
   END SUBROUTINE SParGlobalNumbering
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
SUBROUTINE SParIterBarrier
!-----------------------------------------------------------------------
  INTEGER :: ierr
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
!-----------------------------------------------------------------------
END  SUBROUTINE SParIterBarrier
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
SUBROUTINE SParIterActiveBarrier
!-----------------------------------------------------------------------
  INTEGER :: ierr
  CALL MPI_BARRIER( ParEnv % ActiveComm, ierr )
!-----------------------------------------------------------------------
END  SUBROUTINE SParIterActiveBarrier
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
SUBROUTINE SParIterAllReduceAnd(L)
   LOGICAL :: L
!-----------------------------------------------------------------------
   LOGICAL :: L1
   INTEGER :: ierr

   L1 = L
   CALL MPI_ALLREDUCE(L1,L,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
!-----------------------------------------------------------------------
END  SUBROUTINE SParIterAllReduceAnd
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE SParIterAllReduceOR(L)
   LOGICAL :: L
!-----------------------------------------------------------------------
   LOGICAL :: L1
   INTEGER :: ierr

   L1 = L
   CALL MPI_ALLREDUCE(L1,L,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
!-----------------------------------------------------------------------
END  SUBROUTINE SParIterAllReduceOR
!-----------------------------------------------------------------------



!--------------------------------------------------------------------------
!> Send all of the interface matrix blocks (in NbsIfMatrices) to neighbour
!> processors. This is done only once so there is no need to optimize
!> communication...
!--------------------------------------------------------------------------
  SUBROUTINE ExchangeInterfaces( NbsIfMatrix, RecvdIfMatrix )
    USE Types
    IMPLICIT NONE

    ! Parameters

    TYPE (BasicMatrix_t), DIMENSION(:) :: NbsIfMatrix, RecvdIfMatrix

    ! Local variables

    INTEGER :: i, j, n, sz, req_cnt, ierr, sproc, &
         destproc, rows, cols, TotalSize

    INTEGER, ALLOCATABLE :: requests(:), &
            recv_rows(:), recv_cols(:), neigh(:)

    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
  !*********************************************************************
  n = ParEnv % NumOfNeighbours
  ALLOCATE( neigh(n) )

  n = 0
  DO i=1,ParEnv % PEs
    IF ( ParEnv % IsNeighbour(i) ) THEN
      n = n + 1
      neigh(n)  = i-1
    END IF
  END DO

  totalsize = 0
  DO i = 1, n
    totalsize = totalsize + 1
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows

    IF ( rows > 0 ) THEN
      cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
      totalsize = totalsize + 3*rows + 1 + cols
    END IF
  END DO
  CALL CheckBuffer( totalsize*4+MPI_BSEND_OVERHEAD*8 )

  !
  ! Receive interface sizes:
  !--------------------------
  ALLOCATE( recv_rows(n), recv_cols(n), requests(n) )
  DO i=1,n
    CALL MPI_iRECV( recv_rows(i),1, MPI_INTEGER, neigh(i), &
         2000, MPI_COMM_WORLD, requests(i), ierr )
  END DO

  DO i=1,n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    CALL MPI_BSEND( rows, 1, MPI_INTEGER, &
           destproc, 2000, MPI_COMM_WORLD, ierr )
  END DO
  CALL MPI_WaitAll( n, requests, MPI_STATUSES_IGNORE, ierr )
   
!----------------------------------------------------------------------

  req_cnt = 0
  DO i=1,n
    IF (recv_rows(i)>0) THEN
      req_cnt = req_cnt + 1
      CALL MPI_iRECV( recv_cols(i), 1, MPI_INTEGER, neigh(i), &
         2001, MPI_COMM_WORLD, requests(req_cnt), ierr )
    END IF
  END DO
  !
  ! Send interface sizes:
  !--------------------------

  DO i=1,n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    IF (rows>0) THEN
      cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
      CALL MPI_BSEND( cols,1, MPI_INTEGER, &
           destproc, 2001, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_WaitAll( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )
  
  !----------------------------------------------------------------------
  !
  ! Receive the interface parts
  !
  !----------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     rows = recv_rows(i)
     sproc = neigh(i)
     RecvdIfMatrix(sproc+1) % NumberOfRows = rows

     IF ( rows>0 ) THEN
       cols = recv_cols(i)
       ALLOCATE( RecvdIfMatrix(sproc+1) % Rows(rows+1) )
       ALLOCATE( RecvdIfMatrix(sproc+1) % Diag(rows) )
       ALLOCATE( RecvdIfMatrix(sproc+1) % Cols(cols) )
       ALLOCATE( RecvdIfMatrix(sproc+1) % GRows(rows) )
       ALLOCATE( RecvdIfMatrix(sproc+1) % RowOwner(Rows) )

       req_cnt = req_cnt+1
       CALL MPI_iRECV( RecvdIfMatrix(sproc+1) % GRows, rows,  MPI_INTEGER, &
             sproc, 1002, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    IF ( rows>0 ) THEN
       cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
       CALL MPI_BSEND( NbsIfMatrix(destproc+1) % GRows, &
            rows, MPI_INTEGER, destproc, 1002, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

! -------------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     rows = recv_rows(i)
     sproc = neigh(i)
     RecvdIfMatrix(sproc+1) % NumberOfRows = rows

     IF ( rows>0 ) THEN
       req_cnt = req_cnt+1
       CALL MPI_iRECV( RecvdIfmatrix(sproc+1) % Rows, rows+1, MPI_INTEGER, &
            sproc, 1003, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    IF ( rows>0 ) THEN
       CALL MPI_BSEND( NbsIfMatrix(destproc+1) % Rows, &
            rows+1, MPI_INTEGER, destproc, 1003, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

! -------------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     rows = recv_rows(i)
     sproc = neigh(i)
     RecvdIfMatrix(sproc+1) % NumberOfRows = rows

     IF ( rows>0 ) THEN
       cols = recv_cols(i)
       req_cnt = req_cnt+1
       CALL MPI_iRECV( RecvdIfMatrix(sproc+1) % Cols, cols, MPI_INTEGER, &
            sproc, 1004, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO


  DO i = 1, n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    IF ( rows>0 ) THEN
       cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
       CALL MPI_BSEND( NbsIfMatrix(destproc+1) % Cols, &
            cols, MPI_INTEGER, destproc, 1004, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

! -------------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     rows = recv_rows(i)
     sproc = neigh(i)
     RecvdIfMatrix(sproc+1) % NumberOfRows = rows

     IF ( rows>0 ) THEN
       req_cnt = req_cnt+1
       CALL MPI_iRECV( RecvdIfMatrix(sproc+1) % RowOwner, rows, MPI_INTEGER, &
            sproc, 1005, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    IF ( rows>0 ) THEN
       CALL MPI_BSEND( NbsIfMatrix(destproc+1) % RowOwner, &
            rows, MPI_INTEGER, destproc, 1005, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

! -------------------------------------------------------------------------

  DEALLOCATE( requests, neigh, recv_rows, recv_cols )
!*********************************************************************
END SUBROUTINE ExchangeInterfaces
!*********************************************************************



!--------------------------------------------------------------------------
!> Send all of the interface matrix blocks (in NbsIfMatrices) to neighbour
!> processors. This is done only once so there is no need to optimize
!> communication...
!--------------------------------------------------------------------------
  SUBROUTINE ExchangeIfvalues( NbsIfMatrix, RecvdIfMatrix, &
             NeedMass, NeedDamp, NeedPrec, NeedILU )
    USE Types
    IMPLICIT NONE

    ! Parameters

    LOGICAL :: NeedMass, NeedDamp,NeedPrec,NeedILU
    TYPE (BasicMatrix_t), DIMENSION(:) :: NbsIfMatrix, RecvdIfMatrix

    ! Local variables

    INTEGER :: i, j, n, sz, req_cnt, ierr, sproc, &
         destproc, rows, cols, TotalSize

    INTEGER, ALLOCATABLE :: requests(:), &
            recv_rows(:), recv_cols(:), neigh(:)

    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
  !*********************************************************************
  n = ParEnv % NumOfNeighbours
  ALLOCATE( neigh(n) )

  n = 0
  DO i=1,ParEnv % PEs
    IF ( ParEnv % IsNeighbour(i) ) THEN
      n = n + 1
      neigh(n)  = i-1
    END IF
  END DO

  totalsize = 0
  DO i = 1, n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    IF ( rows == 0 ) THEN
      totalsize = totalsize + 4
    ELSE
      cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
      totalsize = totalsize + 4 + 3*4*rows + 4 + 5*cols
    END IF
  END DO
  CALL CheckBuffer( totalsize+MPI_BSEND_OVERHEAD*8 )

  !
  ! Receive interface sizes:
  !--------------------------
  ALLOCATE( recv_rows(n), recv_cols(n), requests(n) )
  DO i=1,n
    CALL MPI_iRECV( recv_rows(i),1, MPI_INTEGER, neigh(i), &
         2000, MPI_COMM_WORLD, requests(i), ierr )
  END DO

  DO i=1,n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    CALL MPI_BSEND( rows, 1, MPI_INTEGER, &
           destproc, 2000, MPI_COMM_WORLD, ierr )
  END DO
  CALL MPI_WaitAll( n, requests, MPI_STATUSES_IGNORE, ierr )
   
!----------------------------------------------------------------------

  req_cnt = 0
  DO i=1,n
    IF (recv_rows(i)>0) THEN
      req_cnt = req_cnt + 1
      CALL MPI_iRECV( recv_cols(i), 1, MPI_INTEGER, neigh(i), &
         2001, MPI_COMM_WORLD, requests(req_cnt), ierr )
    END IF
  END DO
  !
  ! Send interface sizes:
  !--------------------------

  DO i=1,n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    IF (rows>0) THEN
      cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
      CALL MPI_BSEND( cols,1, MPI_INTEGER, &
           destproc, 2001, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_WaitAll( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )
  
  !----------------------------------------------------------------------
  !
  ! Receive the interface parts
  !
  !----------------------------------------------------------------------

!----------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     rows = recv_rows(i)
     sproc = neigh(i)
     RecvdIfMatrix(sproc+1) % NumberOfRows = rows

     IF ( rows>0 ) THEN
       cols = recv_cols(i)
       ALLOCATE( RecvdIfMatrix(sproc+1) % Rows(rows+1) )
       ALLOCATE( RecvdIfMatrix(sproc+1) % Diag(rows) )
       ALLOCATE( RecvdIfMatrix(sproc+1) % Cols(cols) )
       ALLOCATE( RecvdIfMatrix(sproc+1) % GRows(rows) )
       ALLOCATE( RecvdIfMatrix(sproc+1) % Values(cols) )

!      RecvdIfMatrix(sproc+1) % MassValues => NULL()
       IF ( NeedMass ) &
         ALLOCATE( RecvdIfMatrix(sproc+1) % MassValues(Cols) )

!      RecvdIfMatrix(sproc+1) % DampValues => NULL()
       IF ( NeedDamp ) &
         ALLOCATE( RecvdIfMatrix(sproc+1) % DampValues(Cols) )

!      RecvdIfMatrix(sproc+1) % PrecValues => NULL()
       IF ( NeedPrec ) &
         ALLOCATE( RecvdIfMatrix(sproc+1) % PrecValues(Cols) )

!      RecvdIfMatrix(sproc+1) % ILUValues => NULL()
       IF ( NeedILU ) &
         ALLOCATE( RecvdIfMatrix(sproc+1) % ILUValues(Cols) )

       req_cnt = req_cnt+1
       CALL MPI_iRECV( RecvdIfMatrix(sproc+1) % GRows, rows,  MPI_INTEGER, &
             sproc, 2002, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    IF ( rows>0 ) THEN
       cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
       CALL MPI_BSEND( NbsIfMatrix(destproc+1) % GRows, &
            rows, MPI_INTEGER, destproc, 2002, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

!----------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     rows = recv_rows(i)
     sproc = neigh(i)
     RecvdIfMatrix(sproc+1) % NumberOfRows = rows

     IF ( rows>0 ) THEN
       cols = recv_cols(i)
       req_cnt = req_cnt+1
       CALL MPI_iRECV( RecvdIfmatrix(sproc+1) % Rows, rows+1, MPI_INTEGER, &
            sproc, 2003, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    IF ( rows>0 ) THEN
       cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
       CALL MPI_BSEND( NbsIfMatrix(destproc+1) % Rows, &
            rows+1, MPI_INTEGER, destproc, 2003, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

!----------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     rows = recv_rows(i)
     sproc = neigh(i)
     RecvdIfMatrix(sproc+1) % NumberOfRows = rows

     IF ( rows>0 ) THEN
       cols = recv_cols(i)
       req_cnt = req_cnt+1
       CALL MPI_iRECV( RecvdIfMatrix(sproc+1) % Cols, cols, MPI_INTEGER, &
            sproc, 2004, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    IF ( rows>0 ) THEN
       cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
       CALL MPI_BSEND( NbsIfMatrix(destproc+1) % Cols, &
            cols, MPI_INTEGER, destproc, 2004, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

!----------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     rows = recv_rows(i)
     sproc = neigh(i)
     RecvdIfMatrix(sproc+1) % NumberOfRows = rows

     IF ( rows>0 ) THEN
       cols = recv_cols(i)
       req_cnt = req_cnt+1
       CALL MPI_iRECV( RecvdIfMatrix(sproc+1) % Values, cols, MPI_DOUBLE_PRECISION, &
            sproc, 2005, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    rows = NbsIfMatrix(destproc+1) % NumberOfRows
    IF ( rows>0 ) THEN
       cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
       CALL MPI_BSEND( NbsIfMatrix(destproc+1) % Values, &
            cols, MPI_DOUBLE_PRECISION, destproc, 2005, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

!----------------------------------------------------------------------

  IF ( NeedPrec ) THEN
    req_cnt = 0
    DO i = 1, n
       rows = recv_rows(i)
       sproc = neigh(i)
       RecvdIfMatrix(sproc+1) % NumberOfRows = rows

       IF ( rows>0 ) THEN
         cols = recv_cols(i)
         req_cnt = req_cnt+1
         CALL MPI_iRECV( RecvdIfMatrix(sproc+1) % PrecValues, cols, &
           MPI_DOUBLE_PRECISION, sproc, 2006, MPI_COMM_WORLD, requests(req_cnt), ierr )
       END IF
    END DO

    DO i = 1, n
      destproc = neigh(i)
      rows = NbsIfMatrix(destproc+1) % NumberOfRows
      IF ( rows>0 ) THEN
         cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
         CALL MPI_BSEND( NbsIfMatrix(destproc+1) % PrecValues, &
              cols, MPI_DOUBLE_PRECISION, destproc, 2006, MPI_COMM_WORLD, ierr )
      END IF
    END DO
    CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )
  END IF

!----------------------------------------------------------------------
  IF ( NeedMass ) THEN
    req_cnt = 0
    DO i = 1, n
       rows = recv_rows(i)
       sproc = neigh(i)
       RecvdIfMatrix(sproc+1) % NumberOfRows = rows

       IF ( rows>0 ) THEN
         cols = recv_cols(i)
         req_cnt = req_cnt+1
         CALL MPI_iRECV( RecvdIfMatrix(sproc+1) % MassValues, cols, &
           MPI_DOUBLE_PRECISION, sproc, 2007, MPI_COMM_WORLD, requests(req_cnt), ierr )
       END IF
    END DO

    DO i = 1, n
      destproc = neigh(i)
      rows = NbsIfMatrix(destproc+1) % NumberOfRows
      IF ( rows>0 ) THEN
         cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
         CALL MPI_BSEND( NbsIfMatrix(destproc+1) % MassValues, &
              cols, MPI_DOUBLE_PRECISION, destproc, 2007, MPI_COMM_WORLD, ierr )
      END IF
    END DO
    CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )
  END IF

!----------------------------------------------------------------------

  IF ( NeedDamp ) THEN
    req_cnt = 0
    DO i = 1, n
       rows = recv_rows(i)
       sproc = neigh(i)
       RecvdIfMatrix(sproc+1) % NumberOfRows = rows

       IF ( rows>0 ) THEN
         cols = recv_cols(i)
         req_cnt = req_cnt+1
         CALL MPI_iRECV( RecvdIfMatrix(sproc+1) % DampValues, cols, &
           MPI_DOUBLE_PRECISION, sproc, 2008, MPI_COMM_WORLD, requests(req_cnt), ierr )
       END IF
    END DO

    DO i = 1, n
      destproc = neigh(i)
      rows = NbsIfMatrix(destproc+1) % NumberOfRows
      IF ( rows>0 ) THEN
         cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
         CALL MPI_BSEND( NbsIfMatrix(destproc+1) % DampValues, &
              cols, MPI_DOUBLE_PRECISION, destproc, 2008, MPI_COMM_WORLD, ierr )
      END IF
    END DO
    CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )
  END IF

!----------------------------------------------------------------------

  IF ( NeedILU ) THEN
    req_cnt = 0
    DO i = 1, n
       rows = recv_rows(i)
       sproc = neigh(i)
       RecvdIfMatrix(sproc+1) % NumberOfRows = rows

       IF ( rows>0 ) THEN
         cols = recv_cols(i)
         req_cnt = req_cnt+1
         CALL MPI_iRECV( RecvdIfMatrix(sproc+1) % ILUValues, cols, &
           MPI_DOUBLE_PRECISION, sproc, 2009, MPI_COMM_WORLD, requests(req_cnt), ierr )
       END IF
    END DO

    DO i = 1, n
      destproc = neigh(i)
      rows = NbsIfMatrix(destproc+1) % NumberOfRows
      IF ( rows>0 ) THEN
         cols = NbsIfMatrix(destproc+1) % Rows(rows+1)-1
         CALL MPI_BSEND( NbsIfMatrix(destproc+1) % ILUValues, &
              cols, MPI_DOUBLE_PRECISION, destproc, 2009, MPI_COMM_WORLD, ierr )
      END IF
    END DO
    CALL MPI_Waitall( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )
  END IF

!----------------------------------------------------------------------

  DEALLOCATE( requests, neigh, recv_rows, recv_cols )
!*********************************************************************
END SUBROUTINE ExchangeIfValues
!*********************************************************************


!*********************************************************************
SUBROUTINE ExchangeSourceVec( SourceMatrix, SplittedMatrix, &
            ParallelInfo, SourceVec, op )
!*********************************************************************
  TYPE (SplittedMatrixT) :: SplittedMatrix
  TYPE (Matrix_t) :: SourceMatrix
  TYPE (ParallelInfo_t) :: ParallelInfo
  REAL(KIND=dp), DIMENSION(:) :: SourceVec
  INTEGER, OPTIONAL :: op

  TYPE(vBuff_t), ALLOCATABLE :: recv_buf(:), send_buf(:)

  ! Local variables
  INTEGER :: i, j, k, n, datalen, ierr, sproc, destproc, ind, req_cnt,oper
  INTEGER :: owner, request, totalsize
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

  INTEGER, ALLOCATABLE :: requests(:), recv_size(:), &
        send_size(:), perm(:), neigh(:)
  !*********************************************************************
  n = ParEnv % NumOfNeighbours
  IF ( n<= 0 ) RETURN

  oper = 0 ! 0=sum, 1=min, 2=max
  IF ( PRESENT(op) ) oper=op

  ALLOCATE( neigh(n) )

  n = 0
  DO i=1,ParEnv % PEs
    IF ( ParEnv % IsNeighbour(i) ) THEN
      n = n + 1
      neigh(n)  = i-1
    END IF
  END DO

  ALLOCATE(perm(0:Parenv % Pes-1))
  DO i=1,n
    perm(neigh(i))=i
  END DO

  ALLOCATE( send_size(n), recv_buf(n), send_buf(n) )

  send_size = 0
  DO i = 1, SourceMatrix % NumberOfRows
    DO j=1,SIZE(ParallelInfo % NeighbourList(i) % Neighbours)
      owner = ParallelInfo % NeighbourList(i) % Neighbours(j)
      IF ( owner /= ParEnv % MyPE .AND. ParEnv % Active(owner+1) ) THEN
         owner = perm(owner)
         send_size(owner) = send_size(owner) + 1
      END IF
    END DO
  END DO

  DO i=1,n
    IF ( send_size(i) > 0 ) &
      ALLOCATE(send_buf(i) % ind(send_size(i)),send_buf(i) % vec(send_size(i)))
  END DO

  send_size = 0
  DO i = 1, SourceMatrix % NumberOfRows
    DO j=1,SIZE(ParallelInfo % NeighbourList(i) % Neighbours)
      owner = ParallelInfo % NeighbourList(i) % Neighbours(j)
      IF ( owner /= ParEnv % MyPE .AND. ParEnv % Active(owner+1) ) THEN
         owner = perm(owner)
         send_size(owner) = send_size(owner) + 1
         send_buf(owner) % vec(send_size(owner)) = SourceVec(i)
         send_buf(owner) % ind(send_size(owner)) = ParallelInfo % GlobalDOFs(i)
      END IF
    END DO
  END DO

  totalsize = SUM(send_size)
  CALL CheckBuffer( 3*4*totalsize+3*MPI_BSEND_OVERHEAD )

  !
  ! Receive interface sizes:
  !--------------------------
  ALLOCATE( recv_size(n), requests(n) )
  DO i=1,n
    CALL MPI_iRECV( recv_size(i), 1, MPI_INTEGER, neigh(i), &
          3000, MPI_COMM_WORLD, requests(i), ierr )
  END DO

  !
  ! Send interface sizes:
  !--------------------------
  DO i=1,n
    CALL MPI_BSEND( send_size(i), 1, MPI_INTEGER, neigh(i), &
          3000, MPI_COMM_WORLD, ierr )
  END DO
  CALL MPI_WaitAll( n, requests, MPI_STATUSES_IGNORE, ierr )
  
! --------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     sproc = neigh(i)
     datalen = recv_size(i)
     IF ( datalen > 0 ) THEN
        req_cnt = req_cnt + 1
        ALLOCATE( recv_buf(i) % ind(datalen) )
        CALL MPI_iRECV( recv_buf(i) % Ind, datalen, MPI_INTEGER, sproc, &
                3001, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    datalen = send_size(i)
    IF ( datalen > 0 ) THEN
       CALL MPI_BSEND( send_buf(i) % ind, datalen, &
          MPI_INTEGER, destproc, 3001, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_WaitAll( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

! --------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     sproc = neigh(i)
     dataLen = recv_size(i)
     IF ( datalen > 0 ) THEN
        req_cnt = req_cnt + 1
        ALLOCATE( recv_buf(i) % vec(datalen) )
        CALL MPI_iRECV( recv_buf(i) % vec, datalen, MPI_DOUBLE_PRECISION, &
             sproc, 3002, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    datalen = send_size(i)
    IF ( datalen > 0 ) THEN
       CALL MPI_BSEND( send_buf(i) % vec, datalen, &
          MPI_DOUBLE_PRECISION, destproc, 3002, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_WaitAll( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

! --------------------------------------------------------------------

  DO i=1,n
     datalen = recv_size(i)
     DO j = 1, datalen
       Ind = SearchNode( ParallelInfo, recv_buf(i) % Ind(j),Order=SourceMatrix % Perm )
       IF ( Ind /= -1 ) THEN
!         Ind = SourceMatrix % Perm(Ind)
          IF ( Ind > 0 ) THEN
             SELECT CASE(oper)
             CASE(0)
               SourceVec(Ind) = SourceVec(Ind) + recv_buf(i) % vec(j)
             CASE(1)
               SourceVec(Ind) = MIN(SourceVec(Ind),recv_buf(i) % vec(j))
             CASE(2)
               SourceVec(Ind) = MAX(SourceVec(Ind),recv_buf(i) % vec(j))
             END SELECT
          END IF
       END IF
    END DO
  END DO

  DO i=1,n
    IF (send_size(i)>0) DEALLOCATE(send_buf(i) % Ind, send_buf(i) % Vec)
    IF (recv_size(i)>0) DEALLOCATE(recv_buf(i) % Ind, recv_buf(i) % Vec)
  END DO
  DEALLOCATE( recv_buf, send_buf, recv_size, send_size, requests, neigh, perm )

!*********************************************************************
END SUBROUTINE ExchangeSourceVec
!*********************************************************************



!----------------------------------------------------------------------
!> Exchange right-hand-side elements on the interface with neighbours
!----------------------------------------------------------------------
SUBROUTINE ExchangeRHSIf( SourceMatrix, SplittedMatrix, &
          ParallelInfo, SourceRHS, TargetRHS )

  TYPE (SplittedMatrixT) :: SplittedMatrix
  TYPE (Matrix_t) :: SourceMatrix
  TYPE (ParallelInfo_t) :: ParallelInfo
  REAL(KIND=dp), DIMENSION(:) :: SourceRHS, TargetRHS

  ! Local variables
  TYPE(vBuff_t), ALLOCATABLE :: recv_buf(:), send_buf(:)

  INTEGER :: i, j, k, n, datalen, ierr, sproc, destproc, ind, req_cnt
  INTEGER :: owner, request, totalsize
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

  INTEGER, ALLOCATABLE :: requests(:), recv_size(:), &
        send_size(:), perm(:), neigh(:)
  !*********************************************************************

  n = ParEnv % NumOfNeighbours
  ALLOCATE( neigh(n) )

  n = 0
  DO i=1,ParEnv % PEs
    IF ( ParEnv % IsNeighbour(i) ) THEN
      n = n + 1
      neigh(n)  = i-1
    END IF
  END DO

  ALLOCATE(perm(0:Parenv % Pes-1))
  DO i=1,n
    perm(neigh(i))=i
  END DO

  ALLOCATE( send_size(n), recv_buf(n), send_buf(n) )
  !----------------------------------------------------------------------
  !
  ! Extract the interface elements from SourceRHS to be sent to the
  ! real owner of that element.
  !
  !----------------------------------------------------------------------

  send_size = 0
  DO i = 1, SourceMatrix % NumberOfRows
    owner = ParallelInfo % NeighbourList(i) % Neighbours(1)
    IF ( owner /= ParEnv % MyPE .AND. ParEnv % Active(owner+1) ) THEN
       owner = perm(owner)
       send_size(owner) = send_size(owner) + 1
    END IF
  END DO

  DO i=1,n
    IF ( send_size(i) > 0 ) &
      ALLOCATE(send_buf(i) % ind(send_size(i)), send_buf(i) % vec(send_size(i)) )
  END DO

  send_size = 0
  DO i = 1, SourceMatrix % NumberOfRows
    owner = ParallelInfo % NeighbourList(i) % Neighbours(1)
    IF ( owner /= ParEnv % MyPE .AND. ParEnv % Active(owner+1) ) THEN
       owner = perm(owner)
       send_size(owner) = send_size(owner) + 1
       send_buf(owner) % vec(send_size(owner)) = SourceRHS(i)
       send_buf(owner) % ind(send_size(owner)) = ParallelInfo % GlobalDOFs(i)
    END IF
  END DO

  totalsize = SUM(send_size)
  CALL CheckBuffer( 3*4*totalsize+3*MPI_BSEND_OVERHEAD )

  !
  ! Receive interface sizes:
  !--------------------------
  ALLOCATE( recv_size(n), requests(n) )
  DO i=1,n
    CALL MPI_iRECV( recv_size(i),1, MPI_INTEGER, neigh(i), &
         3100, MPI_COMM_WORLD, requests(i), ierr )
  END DO

  !
  ! Send interface sizes:
  !--------------------------
  DO i=1,n
    CALL MPI_BSEND( send_size(i), 1, MPI_INTEGER, neigh(i), &
         3100, MPI_COMM_WORLD, ierr )
  END DO
  CALL MPI_WaitAll( n, requests, MPI_STATUSES_IGNORE, ierr )
  
! --------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     sproc = neigh(i)
     datalen = recv_size(i)
     IF ( datalen > 0 ) THEN
        req_cnt = req_cnt + 1
        ALLOCATE( recv_buf(i) % ind(datalen) )
        CALL MPI_iRECV( recv_buf(i) % Ind, datalen, MPI_INTEGER, sproc, &
                3101, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    datalen = send_size(i)
    IF ( datalen > 0 ) THEN
       CALL MPI_BSEND( send_buf(i) % ind, datalen, &
          MPI_INTEGER, destproc, 3101, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_WaitAll( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

! --------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     sproc = neigh(i)
     datalen = recv_size(i)
     IF ( datalen > 0 ) THEN
        req_cnt = req_cnt + 1
        ALLOCATE( recv_buf(i) % vec(datalen) )
        CALL MPI_iRECV( recv_buf(i) % vec, datalen, MPI_DOUBLE_PRECISION, &
             sproc, 3102, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    datalen = send_size(i)
    IF ( datalen > 0 ) THEN
       CALL MPI_BSEND( send_buf(i) % vec, datalen, &
          MPI_DOUBLE_PRECISION, destproc, 3102, MPI_COMM_WORLD, ierr )
    END IF
  END DO
  CALL MPI_WaitAll( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

! --------------------------------------------------------------------

  DO i=1,n
     datalen = recv_size(i)
     DO j = 1, datalen
       Ind = SearchNode( ParallelInfo, recv_buf(i) % Ind(j), &
                 Order=SourceMatrix % Perm )
       IF ( Ind /= -1 ) THEN
!         Ind = SourceMatrix % Perm(Ind)
          IF ( Ind > 0 ) &
             SourceRHS(Ind) = SourceRHS(Ind) + recv_buf(i) % vec(j)
       END IF
    END DO
  END DO

  DO i=1,n
    IF (send_size(i)>0) DEALLOCATE(send_buf(i) % Ind, send_buf(i) % Vec)
    IF (recv_size(i)>0) DEALLOCATE(recv_buf(i) % Ind, recv_buf(i) % Vec)
  END DO
  DEALLOCATE( recv_buf, send_buf, recv_size, send_size, requests, neigh, perm )

  ! copy source to target
  j = 0
  DO i = 1, SourceMatrix % NumberOfRows
     IF (ParallelInfo % NeighbourList(i) % Neighbours(1)==ParEnv % MyPE) THEN
        j = j + 1
        TargetRHS(j) = SourceRHS(i)
     END IF
  END DO
!*********************************************************************
END SUBROUTINE ExchangeRHSIf
!*********************************************************************



!----------------------------------------------------------------------
!> Send parts of the result vector to neighbours
!----------------------------------------------------------------------
SUBROUTINE ExchangeResult( SourceMatrix, SplittedMatrix, ParallelInfo, XVec )
  USE types
  IMPLICIT NONE

  TYPE(SplittedMatrixT) :: SplittedMatrix
  TYPE(Matrix_t) :: SourceMatrix
  REAL(KIND=dp), DIMENSION(:) :: XVec
  TYPE (ParallelInfo_t) :: ParallelInfo

  ! Local variables

  TYPE (ResBufferT), POINTER :: CurrRBuf

  TYPE(vBuff_t), ALLOCATABLE :: recv_buf(:)

  INTEGER :: i, j, k, n, datalen, ierr, sproc, destproc, ind, req_cnt
  INTEGER :: owner, request, totalsize
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

  INTEGER, ALLOCATABLE :: requests(:), recv_size(:), &
        send_size(:), perm(:), neigh(:)
  !*********************************************************************

  n = ParEnv % NumOfNeighbours
  ALLOCATE( neigh(n) )

  n = 0
  DO i=1,ParEnv % PEs
    IF ( ParEnv % IsNeighbour(i) ) THEN
      n = n + 1
      neigh(n)  = i-1
    END IF
  END DO

  ALLOCATE(perm(0:Parenv % Pes-1))
  DO i=1,n
    perm(neigh(i))=i
  END DO

  ALLOCATE( recv_size(n), requests(n), recv_buf(n) )

  TotalSize = 0
  DO i = 1, ParEnv % PEs
     totalsize = totalsize + 1
     IF ( ParEnv % IsNeighbour(i) ) THEN
       CurrRBuf => SplittedMatrix % ResBuf(i)
       IF ( ALLOCATED(CurrRBuf % ResInd) ) THEN
         datalen = SIZE(CurrRBuf % ResInd)
         totalsize = totalsize + 1 + 3*datalen
       END IF
     END IF
  END DO
  CALL CheckBuffer( totalsize*4 + 3*MPI_BSEND_OVERHEAD )

  DO i = 1,n
    CALL MPI_iRECV( recv_size(i), 1, MPI_INTEGER, &
        neigh(i), 9000, MPI_COMM_WORLD, requests(i), ierr )
  END DO

  DO i = 1,n
    destproc = neigh(i)
    CurrRBuf => SplittedMatrix % ResBuf(destproc+1)

    datalen = 0
    IF ( ALLOCATED(CurrRBuf % ResInd) ) &
      datalen = SIZE(CurrRBuf % ResInd)

    CALL MPI_BSEND( datalen, 1, MPI_INTEGER, &
         destproc, 9000, MPI_COMM_WORLD, ierr )
  END DO
  CALL MPI_WaitAll( n, requests, MPI_STATUSES_IGNORE, ierr )

! --------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     sproc = neigh(i)
     datalen = recv_size(i)
     IF ( datalen > 0 ) THEN
       req_cnt = req_cnt + 1
       ALLOCATE( recv_buf(i) % ind(datalen) )
       CALL MPI_iRECV( recv_buf(i) % ind, datalen, MPI_INTEGER, sproc, &
               9001, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    CurrRBuf => SplittedMatrix % ResBuf(destproc+1)
    IF ( ALLOCATED(CurrRBuf % ResInd) ) THEN
      datalen = SIZE(CurrRBuf % Resind)
      IF ( datalen > 0 ) THEN
         CALL MPI_BSEND( CurrRBuf % resind, datalen, &
            MPI_INTEGER, destproc, 9001, MPI_COMM_WORLD, ierr )
      END IF
    END IF
  END DO
  CALL MPI_WaitAll( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

! --------------------------------------------------------------------

  req_cnt = 0
  DO i = 1, n
     sproc = neigh(i)
     datalen = recv_size(i)
     IF ( datalen > 0 ) THEN
        req_cnt = req_cnt + 1
        ALLOCATE( recv_buf(i) % vec(datalen) )
        CALL MPI_iRECV( recv_buf(i) % vec, datalen, MPI_DOUBLE_PRECISION, &
             sproc, 9002, MPI_COMM_WORLD, requests(req_cnt), ierr )
     END IF
  END DO

  DO i = 1, n
    destproc = neigh(i)
    CurrRBuf => SplittedMatrix % ResBuf(destproc+1)
    IF ( ALLOCATED(CurrRBuf % ResInd) ) THEN
      datalen = SIZE(CurrRBuf % Resind)
      IF ( datalen > 0 ) THEN
         CALL MPI_BSEND( CurrRBuf % resval, datalen, &
            MPI_DOUBLE_PRECISION, destproc, 9002, MPI_COMM_WORLD, ierr )
      END IF
    END IF
  END DO
  CALL MPI_WaitAll( req_cnt, requests, MPI_STATUSES_IGNORE, ierr )

! --------------------------------------------------------------------

  DO i=1,n
     datalen = recv_size(i)
     DO j = 1, datalen
       ind = SearchNode( ParallelInfo, recv_buf(i) % Ind(j), &
                  Order=SourceMatrix % Perm )
       IF ( ind /= -1 ) THEN
!         ind = SourceMatrix % Perm(ind)
          IF ( ind > 0 ) &
             xvec(ind) = recv_buf(i) % vec(j)
       END IF
    END DO
  END DO

  DO i=1,n
    IF (recv_size(i)>0) DEALLOCATE(recv_buf(i) % Ind, recv_buf(i) % Vec)
  END DO
  DEALLOCATE( recv_buf, recv_size, requests, neigh, perm )

  CALL SparIterActiveBarrier
!*********************************************************************

!*********************************************************************
END SUBROUTINE ExchangeResult
!*********************************************************************
!*********************************************************************


!*********************************************************************
!*********************************************************************
!
!> Build index tables for faster vector element combination (in parallel
!> matrix-vector operation).
!-----------------------------------------------------------------------
SUBROUTINE BuildRevVecIndices( SplittedMatrix )
  IMPLICIT NONE

  TYPE (SplittedMatrixT) :: SplittedMatrix

  ! Local variables

  INTEGER :: i, j, k, m, n, VecLen, ind, ierr, destproc, sproc, TotLen

  TYPE (BasicMatrix_t), POINTER :: CurrIf
  TYPE(Matrix_t), POINTER :: InsideMatrix

  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
  INTEGER :: TotalSize, count
  LOGICAL :: Found
  INTEGER, POINTER :: RevBuff(:),RevInd(:)
  INTEGER, DIMENSION(:), ALLOCATABLE :: GIndices,RowOwner,neigh,L

  TYPE(iBuff_t), ALLOCATABLE :: sbuf(:)

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: tt
#else
  REAL(KIND=dp) :: tt,CPUTime
#endif

  !*********************************************************************
tt = CPUTime()

  n = Parenv % NumOfNeighbours
  ALLOCATE( neigh(n), sbuf(n), L(n) )
  n = 0
  DO i = 1,ParEnv % PEs
    IF ( ParEnv % IsNeighbour(i) ) THEN
      n = n + 1
      neigh(n) = i-1
    END IF
  END DO

  count = 0
  totalsize = 0
  DO i = 1, n
    CurrIf => SplittedMatrix % IfMatrix(neigh(i)+1)
    DO j=1,n
      DO k=1,CurrIf % NumberOfRows
        IF (CurrIf % RowOwner(k)==neigh(j)) totalsize=totalsize+1
      END DO
      count = count+2
      totalsize=totalsize+1
    END DO
  END DO
  CALL CheckBuffer( totalsize*4+MPI_BSEND_OVERHEAD*count )

  InsideMatrix => SplittedMatrix % InsideMatrix

  L = 0
  DO i = 1, n
     CurrIf => SplittedMatrix % IfMatrix(neigh(i)+1)
     DO j=1,n
        destproc = neigh(j) 
        DO k=1,CurrIf % NumberOfRows
           IF ( CurrIf % RowOwner(k) == destproc ) THEN
              L(j) = L(j) + 1
           END IF
        END DO
     END DO
  END DO

  DO i=1,n
    IF ( L(i) > 0 ) ALLOCATE( sbuf(i) % ibuf(L(i)) )
  END DO

  L = 0
  DO i = 1, n
     CurrIf => SplittedMatrix % IfMatrix(neigh(i)+1)
     DO j=1,n
        destproc = neigh(j) 
        DO k=1,CurrIf % NumberOfRows
           IF ( CurrIf % RowOwner(k) == destproc ) THEN
              L(j) = L(j) + 1
              sbuf(j) % ibuf(L(j)) = CurrIf % GRows(k)
           END IF
        END DO
     END DO
  END DO

  DO i=1,n
    CALL MPI_BSEND( L(i), 1, MPI_INTEGER, neigh(i), 4000, &
                MPI_COMM_WORLD, ierr )

    IF ( L(i) > 0 ) THEN
       CALL MPI_BSEND( sbuf(i) % ibuf, L(i), MPI_INTEGER, neigh(i), &
               4001, MPI_COMM_WORLD, ierr )
    END IF
  END DO
!print*,parenv % mype, 'first send: ', CPUTime()-tt
!tt = CPUtime()
!
!
!
  DO i = 1,n
    sproc = neigh(i)
    CALL MPI_RECV( veclen, 1, MPI_INTEGER, sproc, &
         4000, MPI_COMM_WORLD, status, ierr )

    SplittedMatrix % VecIndices(sproc+1) % RevInd => Null()
    IF ( veclen > 0 ) THEN
      ALLOCATE( SplittedMatrix % Vecindices(sproc+1) % RevInd(Veclen) )
      ALLOCATE( Gindices(Veclen) )

      CALL MPI_RECV( Gindices, veclen, MPI_INTEGER, sproc, &
              4001, MPI_COMM_WORLD, status, ierr )

      DO m = 1, VecLen
         ind = SearchIAItem( InsideMatrix %  NumberOfRows, &
           InsideMatrix % GRows, Gindices(m), InsideMatrix % Gorder )
         SplittedMatrix % Vecindices(sproc+1) % RevInd(m) = ind
      END DO
      DEALLOCATE( Gindices )
    END IF 
  END DO

!print*,parenv % mype, 'first recv: ', CPUTime()-tt
!tt = CPUtime()

  DO i=1,n
    CurrIf => SplittedMatrix % IfMatrix(neigh(i)+1)
    DO k = 1,CurrIf % NumberOfRows
      IF ( CurrIf % RowOwner(k) == ParEnv % Mype ) THEN
        Ind = SearchIAItem( InsideMatrix % NumberOfRows, &
            InsideMatrix % GRows, CurrIf % GRows(k), &
                  InsideMatrix % Gorder )
        SplittedMatrix % IfORows(neigh(i)+1) % IfVec(k) = ind
      END IF

      DO m = CurrIf % Rows(k), CurrIf % Rows(k+1) - 1
        Ind = SearchIAItem( InsideMatrix %  NumberOfRows, &
          InsideMatrix % GRows, CurrIf % Cols(m), InsideMatrix % Gorder )
        SplittedMatrix % IfLCols(neigh(i)+1) % IfVec(m) = ind
      END DO
    END DO
  END DO

  DO i=1,n
     IF ( ALLOCATED(sbuf(i) % ibuf) ) DEALLOCATE(sbuf(i) % ibuf)
  END DO
  DEALLOCATE(neigh, sbuf,L)
!print*,parenv % mype, 'secnd recv: ', CPUTime()-tt
!*********************************************************************
END SUBROUTINE BuildRevVecIndices
!*********************************************************************


!*********************************************************************
!*********************************************************************
!
!> Send our part of the interface matrix blocks to neighbours.
!----------------------------------------------------------------------
SUBROUTINE Send_LocIf_Old( SplittedMatrix )

  USE types
  IMPLICIT NONE

  TYPE (SplittedMatrixT) :: SplittedMatrix

  ! Local variables

  INTEGER :: i, j, k, ierr, TotalL
  TYPE (BasicMatrix_t), POINTER :: IfM
  TYPE (IfVecT), POINTER :: IfV
  INTEGER, ALLOCATABLE :: L(:)
  REAL(KIND=dp), ALLOCATABLE :: VecL(:,:)

  !*********************************************************************

  ALLOCATE( L(ParEnv % PEs) )
  L = 0
  TotalL = 0

  DO i = 1, ParEnv % PEs
     IfM => SplittedMatrix % IfMatrix(i)

     DO j=1,ParEnv % PEs
        IF ( .NOT. ParEnv % IsNeighbour(j) ) CYCLE

        DO k=1,IfM % NumberOfRows
           IF ( IfM % RowOwner(k) == j-1 ) THEN
              L(j) = L(j) + 1
              TotalL = TotalL + 1
           END IF
        END DO
     END DO
  END DO

  ALLOCATE( VecL( MAXVAL(L), ParEnv % PEs ) )
  L = 0
  VecL = 0

  CALL CheckBuffer( 12*TotalL )

  DO i = 1, ParEnv % PEs
     IfM => SplittedMatrix % IfMatrix(i)
     IfV => SplittedMatrix % IfVecs(i)

     DO j=1, ParEnv % PEs
        IF ( .NOT. ParEnv % IsNeighbour(j) ) CYCLE

        DO k=1,IfM % NumberOfRows
           IF ( IfM % RowOwner(k) == j-1 ) THEN
              L(j) = L(j) + 1
              VecL(L(j),j) = IfV % IfVec(k)
           END IF
        END DO
     END DO
  END DO

  DO j=1,ParEnv % PEs
     IF ( .NOT. ParEnv % IsNeighbour(j) ) CYCLE

     CALL MPI_BSEND( L(j), 1, MPI_INTEGER, J-1, 6000, &
                MPI_COMM_WORLD, IERR )

     IF ( L(j) > 0 ) THEN
        CALL MPI_BSEND( VecL(1:L(j),j), L(j), MPI_DOUBLE_PRECISION, &
                 J-1, 6001, MPI_COMM_WORLD, ierr )
     END IF
  END DO

  IF ( ALLOCATED(VecL) ) DEALLOCATE( VecL, L )

!*********************************************************************
END SUBROUTINE Send_LocIf_Old
!*********************************************************************


!*********************************************************************
!*********************************************************************
!
!> Receive interface block contributions to vector from neighbours.
!----------------------------------------------------------------------
SUBROUTINE Recv_LocIf_Old( SplittedMatrix, ndim, v )

  USE types
  IMPLICIT NONE

  TYPE (SplittedMatrixT) :: SplittedMatrix
  INTEGER :: ndim
  REAL(KIND=dp), DIMENSION(*) :: v
  REAL(KIND=dp), ALLOCATABLE :: DPBuffer(:)

  SAVE DPBuffer

  ! Local variables

  integer :: i, j, k, ierr, sproc
  integer, dimension(MPI_STATUS_SIZE) :: status

  INTEGER, POINTER :: RevInd(:)
  INTEGER :: VecLen, TotLen

  !*********************************************************************

  IF ( .NOT. ALLOCATED(DPBuffer) ) ALLOCATE(DPBuffer(ndim)) 

  DO i = 1, ParEnv % NumOfNeighbours
     CALL MPI_RECV( VecLen, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
              6000, MPI_COMM_WORLD, status, ierr )

     IF ( VecLen > 0 ) THEN
        sproc = status(MPI_SOURCE)
        RevInd => SplittedMatrix % VecIndices(sproc+1) % RevInd

        IF ( VecLen > SIZE( DPBuffer ) ) THEN
           DEALLOCATE( DPBuffer )
           ALLOCATE( DPBuffer( VecLen ) )
        END IF

        CALL MPI_RECV( DPBuffer, VecLen, MPI_DOUBLE_PRECISION, &
               sproc, 6001, MPI_COMM_WORLD, status, ierr )

        DO k = 1, VecLen
           IF ( RevInd(k) > 0 ) &
              v(RevInd(k)) = v(RevInd(k)) + DPBuffer(k)
        END DO
     END IF
  END DO
!*********************************************************************
END SUBROUTINE Recv_LocIf_Old
!*********************************************************************





!*********************************************************************
!*********************************************************************
!> Send our part of the interface matrix blocks to neighbours.
!----------------------------------------------------------------------
SUBROUTINE Send_LocIf_size( SplittedMatrix, n, neigh )

  USE types
  IMPLICIT NONE

  INTEGER :: n, neigh(:)
  TYPE (SplittedMatrixT) :: SplittedMatrix

  ! Local variables

  INTEGER :: i, j, k, ni, nj, ierr, TotalL
  TYPE (IfVecT), POINTER :: IfV
  TYPE (BasicMatrix_t), POINTER :: IfM

  INTEGER :: L(n)
  !*********************************************************************

  L = 0
  TotalL = 0

  DO ni = 1,n
     i = neigh(ni)+1
     IfM => SplittedMatrix % IfMatrix(i)
     DO nj=1,n
        j = neigh(nj)
        DO k=1,IfM % NumberOfRows
          IF ( IfM % RowOwner(k)==j ) L(nj) = L(nj)+1
        END DO
     END DO
  END DO

  DO nj=1,n
    j = neigh(nj)
    CALL MPI_BSEND( L(nj), 1, MPI_INTEGER, j, 6000, &
               MPI_COMM_WORLD, ierr )
  END DO
!*********************************************************************
END SUBROUTINE Send_LocIf_size
!*********************************************************************


!*********************************************************************
!*********************************************************************
!> Send our part of the interface matrix blocks to neighbours.
!
SUBROUTINE Send_LocIf( SplittedMatrix,n,neigh )

  USE types
  IMPLICIT NONE

  INTEGER :: n,neigh(:)
  TYPE (SplittedMatrixT) :: SplittedMatrix

  ! Local variables

  INTEGER :: i, j, k, ni, nj, ierr, TotalL
  TYPE (IfVecT), POINTER :: IfV
  TYPE (BasicMatrix_t), POINTER :: IfM

  TYPE(Buff_t), ALLOCATABLE, SAVE :: VecL(:)

  INTEGER :: L(n)
  !*********************************************************************

  L = 0
  TotalL = 0

  DO ni = 1, n
     i = neigh(ni)+1
     IfM => SplittedMatrix % IfMatrix(i)

     DO nj=1,n
        j = neigh(nj) 
        DO k=1,IfM % NumberOfRows
           IF ( IfM % RowOwner(k) == j ) THEN
              L(nj) = L(nj) + 1
              TotalL = TotalL + 1
           END IF
        END DO
     END DO
  END DO

  CALL CheckBuffer( 12*TotalL )

  IF ( .NOT. ALLOCATED(Vecl) ) THEN
    ALLOCATE( Vecl(n) )
    DO i=1,n
      ALLOCATE( Vecl(i) % Rbuf(L(i)) )
    END DO
  ELSE
    IF ( SIZE(Vecl)<n ) THEN
      DO i=1,SIZE(Vecl)
        IF ( ALLOCATED(Vecl(i) % rbuf) ) &
          DEALLOCATE( Vecl(i) % rbuf )
      END DO
      DEALLOCATE( Vecl )

      ALLOCATE( Vecl(n) )
      DO i=1,n
        ALLOCATE( Vecl(i) % rbuf(L(i)) )
      END DO
    ELSE
      DO i=1,n
        IF ( .NOT. ALLOCATED(Vecl(i) % rbuf) ) THEN
          ALLOCATE( Vecl(i) % rbuf(L(i)) )
        ELSE IF ( SIZE(Vecl(i) % rbuf) < L(i) ) THEN
          DEALLOCATE( Vecl(i) % rbuf )
          ALLOCATE( Vecl(i) % rbuf(L(i)) )
        END IF
      END DO
    END IF
  END IF

  L = 0
  DO ni = 1, n
     i = neigh(ni)+1
     IfV => SplittedMatrix % IfVecs(i)
     IfM => SplittedMatrix % IfMatrix(i)

     DO nj=1, n
        j = neigh(nj)
        DO k=1,IfM % NumberOfRows
           IF ( IfM % RowOwner(k) == j ) THEN
              L(nj) = L(nj) + 1
              VecL(nj) % rbuf(L(nj)) = IfV % IfVec(k)
           END IF
        END DO
     END DO
  END DO

  DO nj=1,n
     IF ( L(nj) > 0 ) THEN
       CALL MPI_BSEND( VecL(nj) % rbuf, L(nj), MPI_DOUBLE_PRECISION, &
                Neigh(nj), 6001, MPI_COMM_WORLD, Ierr )
     END IF
  END DO
!*********************************************************************
END SUBROUTINE Send_LocIf
!*********************************************************************


!*********************************************************************
!*********************************************************************
!
!> Receive interface block contributions to vector from neighbours
!
SUBROUTINE Recv_LocIf_size( n, neigh, sizes )
  USE Types
  IMPLICIT NONE

  INTEGER :: sizes(:), neigh(:), n

  ! Local variables

  integer :: i, j, k, ni,nj,ierr, sproc
  integer, dimension(MPI_STATUS_SIZE) :: status

  INTEGER :: VecLen
  !*********************************************************************
  DO i=1, ParEnv % NumOfNeighbours
     CALL MPI_RECV( Veclen, 1, MPI_INTEGER, neigh(i), &
            6000, MPI_COMM_WORLD, status, ierr )
     sizes(i) = Veclen
  END DO
!*********************************************************************
END SUBROUTINE Recv_LocIf_size
!*********************************************************************


!*********************************************************************
!*********************************************************************
!
!> Receive interface block contributions to vector from neighbours
!
SUBROUTINE Recv_LocIf( SplittedMatrix, n, neigh, sizes, requests, buffer )
  uSE Types
  IMPLICIT NONE

  TYPE (SplittedMatrixT) :: SplittedMatrix
  TYPE(Buff_t) ::  buffer(:)
  INTEGER :: n, sizes(:), requests(:), neigh(:)

  ! Local variables

  INTEGER :: VecLen, TotLen
  integer :: i, j, k, ni, ierr, sproc

  INTERFACE 
    SUBROUTINE MPI_IRECV( buf,size,type,proc,tag,comm,req,ierr )
       USE Types
       REAL(KIND=dp)::buf(*)
       INTEGER :: size,type,proc,tag,comm,req,ierr
    END SUBROUTINE MPI_IRECV
  END INTERFACE

  !*********************************************************************

  DO ni = 1, n
    IF ( sizes(ni)>0 ) THEN
      i = neigh(ni)
      CALL MPI_iRECV( buffer(ni) % rbuf,sizes(ni),MPI_DOUBLE_PRECISION, &
              i, 6001, MPI_COMM_WORLD, requests(ni), Ierr)
    END IF
  END DO
!*********************************************************************
END SUBROUTINE Recv_LocIf
!*********************************************************************


!*********************************************************************
!*********************************************************************
!
!> Receive interface block contributions to vector from neighbours
!
SUBROUTINE Recv_LocIf_Wait( SplittedMatrix, ndim, v, n, neigh, &
               sizes, requests, buffer )
  USE Types
  IMPLICIT NONE

  TYPE (SplittedMatrixT) :: SplittedMatrix
  REAL(KIND=dp), DIMENSION(*) :: v
  TYPE(Buff_t) ::  buffer(:)
  INTEGER :: ndim, n, neigh(:), sizes(:), requests(:)

  ! Local variables

  integer :: i, j, k, ni, nj, ierr, sproc
  integer, dimension(MPI_STATUS_SIZE) :: status

  INTEGER :: Completed, Flag, active_req(n), active_n(n), active_cnt

  INTEGER :: VecLen, TotLen
  INTEGER, POINTER :: RevInd(:)

  !*********************************************************************

  completed  = 0
  active_cnt = 0
  DO ni=1,n
    IF ( sizes(ni) <= 0 ) THEN
      completed = completed+1
    ELSE
      active_cnt = active_cnt + 1
      active_n(active_cnt) = ni
      active_req(active_cnt) = requests(ni)
    END IF
  END DO

#if 1
  CALL MPI_Waitall( active_cnt, active_req, MPI_STATUSES_IGNORE, ierr )

  DO j=1,active_cnt
    ni = active_n(j)
    i = neigh(ni)+1
    RevInd => SplittedMatrix % VecIndices(i) % RevInd
    DO k = 1, sizes(ni)
      IF ( RevInd(k) > 0 ) &
         v(RevInd(k)) = v(RevInd(k))+buffer(ni) % rbuf(k)
    END DO
  END DO
#else
  DO WHILE( completed<n )
    CALL MPI_Waitany( active_cnt, active_req, ni, status, ierr )

    ni = active_n(ni)
    i = neigh(ni) + 1
    RevInd => SplittedMatrix % VecIndices(i) % RevInd
    DO k = 1, sizes(ni)
      IF ( RevInd(k) > 0 ) &
         v(RevInd(k)) = v(RevInd(k))+buffer(ni) % rbuf(k)
    END DO
    completed = completed + 1
  END DO
#endif
!*********************************************************************
END SUBROUTINE Recv_LocIf_Wait
!*********************************************************************



!*********************************************************************
SUBROUTINE SParActiveSUM(tsum, oper)
   INTEGER :: oper
   REAL(KIND=dp) :: tsum
!*********************************************************************
   INTEGER :: ierr
   REAL(KIND=dp) :: ssum

   IF ( COUNT(ParEnv % Active)<= 1 ) RETURN

   ssum = tsum
   SELECT CASE(oper)
   CASE(0)
     CALL MPI_ALLREDUCE( ssum, tsum, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, ParEnv % ActiveComm, ierr )
   CASE(1)
     CALL MPI_ALLREDUCE( ssum, tsum, 1, MPI_DOUBLE_PRECISION, &
            MPI_MIN, ParEnv % ActiveComm, ierr )
   CASE(2)
     CALL MPI_ALLREDUCE( ssum, tsum, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, ParEnv % ActiveComm, ierr )
   END SELECT
!*********************************************************************
END SUBROUTINE SParActiveSUM
!*********************************************************************


!*********************************************************************
!*********************************************************************
!
! Compute global dot product of vectors x and y
!
!*********************************************************************
FUNCTION SParDotProd( ndim, x, xind, y, yind ) RESULT(dres)
!*********************************************************************
  IMPLICIT NONE

  ! Parameters

  INTEGER :: ndim, xind, yind
  REAL(KIND=dp) :: x(*)
  REAL(KIND=dp) :: y(*)
  REAL(KIND=dp) :: dres

  ! Local variables

  REAL(KIND=dp) :: s
  INTEGER :: i

  !*********************************************************************
   dres = 0
   DO i = 1, ndim
      dres = dres + y(i) * x(i)
   END DO
   CALL SParActiveSUM(dres,0)
!*********************************************************************
END FUNCTION SParDotProd
!*********************************************************************


!*********************************************************************
!*********************************************************************
!
!> Compute global 2-norm of vector x
!
FUNCTION SParNorm( ndim, x, xind ) RESULT(dres)
  IMPLICIT NONE

  ! Parameters

  INTEGER :: ndim, xind
  REAL(KIND=dp) :: x(*)
  REAL(KIND=dp) :: dres

  ! Local variables
  INTEGER :: i
  !*********************************************************************
  dres = 0
  DO i = 1, ndim
    dres = dres + x(i)*x(i)
  END DO
  CALL SParActiveSUM(dres,0)
  dres = SQRT(dres)
!*********************************************************************
END FUNCTION SParNorm
!*********************************************************************


!*********************************************************************
!*********************************************************************
!
!> Compute global dot product of vectors x and y
!
FUNCTION SParCDotProd( ndim, x, xind, y, yind ) result (dres)

  IMPLICIT NONE

  ! Parameters

  INTEGER :: ndim, xind, yind
  COMPLEX(KIND=dp) :: x(*)
  COMPLEX(KIND=dp) :: y(*)
  COMPLEX(KIND=dp) :: dres


  ! Local variables

  COMPLEX(KIND=dp) :: dsum
  INTEGER :: ierr, i, MinActive
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

  !*********************************************************************
  dres = 0.0d0
  IF ( xind == 1 .AND. yind  == 1 ) THEN
     DO i = 1, ndim
        dres = dres + dconjg(x(i)) * y(i)
     END DO
  ELSE
     CALL Fatal( 'SParCDotProd', 'xind or yind not 1' )
  END IF

  dsum = dres
  CALL MPI_ALLREDUCE( dsum, dres, 1, MPI_DOUBLE_COMPLEX, &
            MPI_SUM, ParEnv % ActiveComm, ierr )
!*********************************************************************
END FUNCTION SParCDotProd
!*********************************************************************


!*********************************************************************
!*********************************************************************
!
!> Compute global 2-norm of vector x
!
FUNCTION SParCNorm( ndim, x, xind ) result (norm)
  IMPLICIT NONE

  ! Parameters

  INTEGER :: ndim, xind
  REAL(KIND=dp) :: norm
  COMPLEX(KIND=dp) :: x(*)

  ! Local variables
  INTEGER :: i

  !*********************************************************************
  norm = 0.0d0
  DO i = 1, ndim
     norm = norm + REAL(x(i))**2 + AIMAG(x(i))**2
  END DO
  CALL SparActiveSUM(norm,0)
  norm = SQRT(norm)
!*********************************************************************
END FUNCTION SParCNorm
!*********************************************************************



!*********************************************************************
!*********************************************************************
!
!> Finalize MPI environment
!
SUBROUTINE ParEnvFinalize()
  IMPLICIT NONE

  ! local variables

  INTEGER :: ierr

  !*********************************************************************
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  CALL MPI_FINALIZE( ierr )

  IF ( ierr /= 0 ) THEN
     WRITE( Message, * ) 'MPI Finalization failed ! (ierr=', ierr, ')'
     CALL Fatal( 'ParEnvFinalize', Message )
  END IF
!*********************************************************************
END SUBROUTINE ParEnvFinalize
!*********************************************************************



!*********************************************************************
!
!> Search an element QueriedNode from an ordered set Nodes and return
!> Index to Nodes structure. Return value -1 means QueriedNode was
!> not found.
!
FUNCTION SearchNode( ParallelInfo, QueriedNode, First, Last,Order ) RESULT ( Indx )

  USE Types
  IMPLICIT NONE

  TYPE (ParallelInfo_t) :: ParallelInfo
  INTEGER :: QueriedNode, Indx
  INTEGER, OPTIONAL :: First,Last, Order(:)

  ! Local variables

  INTEGER :: Lower, Upper, Lou, i, U, L

  !*********************************************************************

  Indx = -1
  IF(PRESENT(Order)) THEN
    Upper = SIZE(Order)
  ELSE
    Upper = SIZE(ParallelInfo % GlobalDOFs)
  END IF
  Lower = 1
  IF ( PRESENT( Last  ) ) Upper = Last
  IF ( PRESENT( First ) ) Lower = First

  ! Handle the special case

  IF ( Upper == 0 ) RETURN

  IF (PRESENT(Order)) THEN
    DO WHILE(Order(Lower)==0)
      Lower=Lower+1
    END DO
    DO WHILE(Order(Upper)==0)
      Upper=Upper-1
    END DO
  END IF

10 CONTINUE
  L = Lower
  IF (PRESENT(Order)) L=Order(L)
  U = Upper
  IF (PRESENT(Order)) U=Order(U)

  IF (ParallelInfo % GlobalDOFs(L)==QueriedNode) THEN
     Indx = L; RETURN
  END IF

  IF (ParallelInfo % GlobalDOFs(U)==QueriedNode) THEN
     Indx = U; RETURN
  END IF

  IF ((Upper-Lower)>1) THEN
     L = ISHFT((Upper+Lower),-1)
     U = L

     IF(PRESENT(Order)) THEN
       DO WHILE(Order(L)==0)
         L=L-1
         IF(L<=Lower) EXIT
       END DO
       IF(L<=Lower) THEN
         DO WHILE(Order(U)==0)
           U=U+1
           IF(U>=Upper) RETURN
         END DO
         L=U
       END IF
       U=Order(L)
     END IF

     IF(U==0) RETURN

     IF (ParallelInfo % GlobalDOFs(U)<QueriedNode) THEN
        Lower = L
     ELSE
        Upper = L
     END IF
     GOTO 10
  END IF
  RETURN
!*********************************************************************
END FUNCTION SearchNode
!*********************************************************************


!*********************************************************************
!*********************************************************************
!
!> Search an element Item from an ordered integer array(N) and return
!> Index to that array element. Return value -1 means Item was not found.
!
FUNCTION SearchIAItem( N, IArray, Item, SortOrder, sIndx ) RESULT ( Indx )

  USE types
  IMPLICIT NONE

  INTEGER :: Item, Indx, i
  INTEGER :: N
  INTEGER, DIMENSION(:) :: IArray
  INTEGER, OPTIONAL :: SortOrder(:), sIndx

  ! Local variables

  INTEGER :: Lower, Upper, lou

  !*********************************************************************

  Indx = -1
  Upper =  N
  Lower =  1

  ! Handle the special case

  IF ( Upper == 0 ) RETURN

  IF ( .NOT. PRESENT(SortOrder) ) THEN
     Indx = SearchIAItemLinear( n,IArray,Item )
     RETURN
  END IF

  DO WHILE( .TRUE. )
     IF ( IArray(Lower) == Item ) THEN
        Indx = Lower
        EXIT
     ELSE IF ( IArray(Upper) == Item ) THEN
        Indx = Upper
        EXIT
     END IF

     IF ( (Upper - Lower) > 1 ) THEN
        Lou = ISHFT((Upper + Lower), -1)
        IF ( IArray(lou) < Item ) THEN
           Lower = Lou
        ELSE
           Upper = Lou
        END IF
     ELSE
        EXIT
     END IF
  END DO

  IF (PRESENT(sIndx) ) sIndx=Indx
  IF (Indx > 0) Indx  = SortOrder(Indx)
  RETURN
!*********************************************************************
END FUNCTION SearchIAItem
!*********************************************************************


!*********************************************************************
!> Search an element Item from an ordered integer array(N) and return
!> Index to that array element. Return value -1 means Item was not found.
!
FUNCTION SearchIAItemLinear( N, IArray, Item ) RESULT ( Indx )

  USE types
  IMPLICIT NONE

  INTEGER :: N
  INTEGER, DIMENSION(*) :: IArray
  INTEGER :: Item, Indx, i

  ! Local variables

  INTEGER :: Lower, Upper, lou
  !*********************************************************************

  Indx = -1
  DO i=1,N
     IF ( IArray(i) == Item ) THEN
       Indx = i
       RETURN
     END IF
  END DO
!*********************************************************************
END FUNCTION SearchIAItemLinear
!*********************************************************************


END MODULE SParIterComm

!> \}
