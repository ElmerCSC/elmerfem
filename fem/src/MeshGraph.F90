!/******************************************************************************
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! *
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write
! *  to the Free Software Foundation, Inc., 51 Franklin Street,
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Graph construction, coloring and thread work-share utilities.
!>  Extracted from MeshUtils.
!------------------------------------------------------------------------------

MODULE MeshGraph

    USE MeshBasics
    IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE ElmerMeshToDualGraph(Mesh, DualGraph, UseBoundaryMesh)
    IMPLICIT NONE

    TYPE(Mesh_t) :: Mesh
    TYPE(Graph_t) :: DualGraph
    LOGICAL, OPTIONAL :: UseBoundaryMesh

    TYPE(Element_t), POINTER :: Element, Elements(:)

    ! MESH DATA
    ! Mesh (CRS format)
    INTEGER, ALLOCATABLE :: eptr(:), eind(:)
    INTEGER :: nelem
    ! Vertex to element map (CRS format)
    INTEGER, ALLOCATABLE :: vptr(:), vind(:)
    INTEGER :: nvertex

    ! WORK ARRAYS
    ! Pointers to vertex-element maps of the current element
    INTEGER, ALLOCATABLE :: ptrli(:), ptrti(:)
    ! Neighbour indices
    INTEGER, ALLOCATABLE :: neighind(:)
    ! ARRAY MERGE: map for merge
    INTEGER, ALLOCATABLE :: wrkmap(:)

    TYPE :: IntTuple_t
      INTEGER :: i1, i2
    END type IntTuple_t

    TYPE(IntTuple_t), ALLOCATABLE :: wrkheap(:)

    ! OpenMP thread block leads for work division
    INTEGER, ALLOCATABLE :: thrblk(:)
    ! Work indices
    INTEGER, ALLOCATABLE :: wrkind(:), wrkindresize(:)
    INTEGER :: nwrkind

    ! Variables
    INTEGER :: i, dnnz, eid, nl, nli, nti, nn, nv, nthr, &
            te, thrli, thrti, vli, vti, TID, allocstat
    INTEGER :: mapSizePad, maxNodesPad, neighSizePad
    LOGICAL :: Boundary

    INTEGER, PARAMETER :: HEAPALG_THRESHOLD = 24

    CALL Info('ElmerMeshToDualGraph','Creating a dual graph for the mesh',Level=8)

    Boundary = .FALSE.
    IF (Present(UseBoundaryMesh)) Boundary = UseBoundaryMesh

    ! Pointers to mesh data
    IF (.NOT. Boundary) THEN
       nelem = Mesh % NumberOfBulkElements
       nvertex = Mesh % NumberOfNodes
       Elements => Mesh % Elements
    ELSE
       nelem = Mesh % NumberOfBoundaryElements
       nvertex = Mesh % NumberOfNodes
       Elements => Mesh % Elements(&
            Mesh % NumberOfBulkElements+1:Mesh % NumberOfBulkElements+nelem)
    END IF

    ! Initialize dual mesh size and number of nonzeroes
    DualGraph % n = nelem
    dnnz = 0

    ! Copy mesh to CRS structure
    ALLOCATE(eptr(nelem+1), eind(nelem*Mesh % MaxElementNodes), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate mesh structure!')

    eptr(1)=1 ! Fortran numbering
    DO i=1, nelem
      Element => Elements(i)
      nl = Element % TYPE % NumberOfNodes
      nli = eptr(i) ! Fortran numbering
      nti = nli+nl-1
      eind(nli:nti) = Element % NodeIndexes(1:nl) ! Fortran numbering
      eptr(i+1) = nli+nl
    END DO

    ! Construct vertex to element list (in serial!)
    CALL VertexToElementList(nelem, nvertex, eptr, eind, vptr, vind)

    ! Allocate pointers to dual mesh
    ALLOCATE(DualGraph % ptr(nelem+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate dual mesh!')

    ! Divide work by number of rows in the vertex graph
    nthr = 1 
    !$ nthr = omp_get_max_threads()

    ! Load balance the actual work done by threads (slow)
    ! CALL ThreadLoadBalanceElementNeighbour(nthr, nelem, eptr, eind, vptr, thrblk)
    CALL ThreadStaticWorkShare(nthr, nelem, thrblk)

    !$OMP PARALLEL SHARED(nelem, nvertex, eptr, eind, &
    !$OMP                 vptr, vind, Mesh, DualGraph, &
    !$OMP                 nthr, thrblk, dnnz) &
    !$OMP PRIVATE(i, eid, nli, nti, nn, nv, vli, vti, te, &
    !$OMP         maxNodesPad, neighSizePad, ptrli, ptrti, &
    !$OMP         wrkheap, wrkmap, neighind, &
    !$OMP         wrkind, nwrkind, wrkindresize, allocstat, &
    !$OMP         mapSizePad, thrli, thrti, TID) NUM_THREADS(nthr) &
    !$OMP DEFAULT(NONE)

    TID = 1
    !$ TID = OMP_GET_THREAD_NUM()+1

    ! Ensure that the vertex to element lists are sorted
    !$OMP DO 
    DO i=1,nvertex
      vli = vptr(i)
      vti = vptr(i+1)-1

      CALL Sort(vti-vli+1, vind(vli:vti))
    END DO
    !$OMP END DO NOWAIT

    ! Allocate work array (local to each thread)
    maxNodesPad = IntegerNBytePad(Mesh % MaxElementNodes, 8)
    neighSizePad = IntegerNBytePad(Mesh % MaxElementNodes*20, 8)

    ! Pointers to vertex maps
    ALLOCATE(neighind(neighSizePad), &
            ptrli(maxNodesPad), ptrti(maxNodesPad), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate local workspace!')
    ! Initialize neighbour indices
    neighind = 0

    IF (nthr >= HEAPALG_THRESHOLD) THEN
      ! With multiple threads, use heap based merge
      ALLOCATE(wrkheap(maxNodesPad), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
              'Unable to allocate local workspace!')
    ELSE
      ! With a small number of threads, use map -based merge
      mapSizePad = IntegerNBytePad(nelem, 8)
      ALLOCATE(wrkmap(mapSizePad), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
              'Unable to allocate local workspace!')
      ! Initialize local map
      wrkmap=0
    END IF

    ! Allocate local list for results
    nwrkind = 0
    ALLOCATE(wrkind(nelem/nthr*20), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate local workspace!')

    ! Ensure that all the threads have finished sorting the vertex indices
    !$OMP BARRIER

    ! Get thread indices
    thrli = thrblk(TID)
    thrti = thrblk(TID+1)

    ! For each element
    DO eid=thrli,thrti-1
      nli = eptr(eid)
      nti = eptr(eid+1)-1
      nv = nti-nli+1

      ! Get pointers to vertices related to the nodes of the element
      te = 0
      DO i=nli,nti
        ptrli(i-nli+1)=vptr(eind(i))
        ptrti(i-nli+1)=vptr(eind(i)+1) ! NOTE: This is to make comparison cheaper
        te = te + ptrti(i-nli+1)-ptrli(i-nli+1)
      END DO

      ! Allocate neighind large enough
      IF (SIZE(neighind)<te) THEN
        DEALLOCATE(neighind)
        neighSizePad = IntegerNBytePad(te,8)
        ALLOCATE(neighind(neighSizePad), STAT=allocstat)
        neighind = 0
      END IF

      ! Merge vertex lists (multi-way merge of ordered lists)
      IF (nthr >= HEAPALG_THRESHOLD) THEN
        CALL kWayMergeHeap(eid, nv, ptrli, ptrti, &
                te, vind, nn, neighind, wrkheap)
      ELSE
        CALL kWayMergeArray(eid, nv, ptrli, ptrti, &
                te, vind, nn, neighind, wrkmap)
      END IF

      ! Add merged list to final list of vertices
      IF (nn+nwrkind>SIZE(wrkind)) THEN
        ALLOCATE(wrkindresize(MAX(nn+nwrkind,2*SIZE(wrkind))), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
                'Unable to allocate local workspace!')
        wrkindresize(1:nwrkind)=wrkind(1:nwrkind)
        DEALLOCATE(wrkind)
        CALL MOVE_ALLOC(wrkindresize, wrkind)
      END IF
      wrkind(nwrkind+1:nwrkind+nn) = neighind(1:nn)
      nwrkind = nwrkind + nn

      ! Store number of row nonzeroes
      DualGraph % ptr(eid)=nn
    END DO

    ! Get the global size of the dual mesh
    !$OMP DO REDUCTION(+:dnnz)
    DO i=1,nthr
      dnnz = nwrkind
    END DO
    !$OMP END DO

    ! Allocate memory for dual mesh indices
    !$OMP SINGLE
    ALLOCATE(DualGraph % ind(dnnz), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate dual mesh!')
    ! ptr stores row counts, build crs pointers from them
    CALL ComputeCRSIndexes(nelem, DualGraph % ptr)
    !$OMP END SINGLE

    DualGraph % ind(&
            DualGraph % ptr(thrli):DualGraph % ptr(thrti)-1)=wrkind(1:nwrkind)

    IF (nthr >= HEAPALG_THRESHOLD) THEN
      DEALLOCATE(wrkheap, STAT=allocstat)
    ELSE
      DEALLOCATE(wrkmap, STAT=allocstat)
    END IF
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to deallocate local workspace!')
    DEALLOCATE(neighind, ptrli, ptrti, wrkind)

    !$OMP END PARALLEL

    ! Deallocate the rest of memory
    DEALLOCATE(eind, eptr, vptr, vind, thrblk)

    CALL Info('ElmerMeshToDualGraph','Dual graph created with size '//I2S(dnnz),Level=8)


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

    ! k-way merge with an array
    SUBROUTINE kWayMergeArray(node, nv, ptrli, ptrti, te, vind, &
            nn, neighind, map)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node, nv
      INTEGER :: ptrli(:)
      INTEGER, INTENT(IN) ::ptrti(:), te
      INTEGER, INTENT(IN) :: vind(:)
      INTEGER, INTENT(OUT) :: nn
      INTEGER :: neighind(:)
      INTEGER :: map(:)

      INTEGER :: i, j, k, vindi

      ! Merge nv lists using a map (i.e. an array)
      nn = 1
      DO i=1,nv
        DO j=ptrli(i), ptrti(i)-1
          vindi = vind(j)
          ! Put element to map if it is not already there
          IF (map(vindi)==0 .AND. vindi /= node) THEN
            neighind(nn)=vindi
            ! Increase counter
            map(vindi)=1
            nn=nn+1
          END IF
        END DO
      END DO
      nn=nn-1

      ! Clear map
      DO i=1,nn
        map(neighind(i)) = 0
      END DO
    END SUBROUTINE kWayMergeArray

    ! k-way merge with an actual heap
    SUBROUTINE kWayMergeHeap(node, nv, ptrli, ptrti, te, vind, &
            nn, neighind, heap)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node, nv
      INTEGER :: ptrli(:)
      INTEGER, INTENT(IN) ::ptrti(:), te
      INTEGER, INTENT(IN) :: vind(:)
      INTEGER, INTENT(OUT) :: nn
      INTEGER :: neighind(:)
      TYPE(IntTuple_t) :: heap(:)

      TYPE(IntTuple_t) :: tmp
      INTEGER :: ii, l, r, mind, ll, tmpval, tmpind

      ! Local variables
      INTEGER :: i, e, nzheap, vindi, lindi, pind

      ! Put elements to heap
      nzheap = 0
      DO i=1,nv
        IF (ptrli(i)<ptrti(i)) THEN
          heap(i) % i1 = vind(ptrli(i))
          heap(i) % i2= i
          ptrli(i) = ptrli(i)+1
          nzheap = nzheap+1
        END IF
      END DO

      ! Build heap
      DO ii=(nzheap/2), 1, -1
        i = ii
        ! CALL BinaryHeapHeapify(heap, nzheap, i)
        DO 
          ! Find index of the minimum element
          IF (2*i<=nzheap) THEN
            IF (heap(2*i) % i1 < heap(i) % i1) THEN
              mind = 2*i
            ELSE
              mind = i
            END IF
            IF (2*i+1<=nzheap) THEN
              IF (heap(2*i+1) % i1 < heap(mind) % i1) mind = 2*i+1
            END IF
          ELSE
            mind = i
          END IF

          IF (mind == i) EXIT

          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        END DO
      END DO

      pind = -1
      nn = 1
      DO e=1,te
        ! Pick the first element from heap
        vindi = heap(1) % i1
        lindi = heap(1) % i2

        ! Remove duplicates
        IF (vindi /= pind .AND. vindi /= node) THEN
          neighind(nn) = vindi
          pind = vindi
          nn = nn+1
        END IF

        ! Add new element from list (if any)
        IF (ptrli(lindi) < ptrti(lindi)) THEN
          heap(1) % i1 = vind(ptrli(lindi))
          heap(1) % i2 = lindi
          ptrli(lindi) = ptrli(lindi)+1
        ELSE
          heap(1) % i1 = heap(nzheap) % i1
          heap(1) % i2 = heap(nzheap) % i2
          nzheap=nzheap-1
        END IF
        ! CALL BinaryHeapHeapify(heap, nzheap, 1)
        i = 1

        DO 
          ! Find the index of the minimum element
          ii = 2*i
          mind = i
          IF (ii+1<=nzheap) THEN
            ! Elements 2*i and 2*i+1 can be tested
            IF (heap(ii) % i1 < heap(i) % i1) mind = ii
            IF (heap(ii+1) % i1 < heap(mind) % i1) mind = ii+1
          ELSE IF (ii<=nzheap) THEN
            ! Element ii can be tested
            IF (heap(ii) % i1 < heap(i) % i1) mind = ii
          END IF

          IF (mind == i) EXIT

          ! Bubble down the element
          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        END DO

      END DO
      nn=nn-1
    END SUBROUTINE kWayMergeHeap

    SUBROUTINE BinaryHeapHeapify(heap, nelem, sind)
      IMPLICIT NONE
      TYPE(IntTuple_t) :: heap(:)
      INTEGER, INTENT(IN) :: nelem
      INTEGER, INTENT(IN) :: sind

      INTEGER :: i, l, r, mind
      TYPE(IntTuple_t) :: tmp

      i = sind
      DO
        l = 2*i
        r = 2*i+1
        ! Find index of the minimum element
        mind = i
        IF (l <= nelem) THEN
          IF (heap(l) % i1 < heap(i) % i1) mind = l
        END IF
        IF (r <= nelem) THEN
          IF (heap(r) % i1 < heap(mind) % i1) mind = r
        END IF

        IF (mind /= i) THEN
          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        ELSE
          EXIT
        END IF
      END DO
    END SUBROUTINE BinaryHeapHeapify

    FUNCTION BinaryHeapIsHeap(heap, nelem) RESULT(heaporder)
      IMPLICIT NONE
      TYPE(IntTuple_t) :: heap(:)
      INTEGER, INTENT(IN) :: nelem
      LOGICAL :: heaporder

      INTEGER :: i, l, r

      heaporder = .TRUE.

      DO i=(nelem/2), 1, -1
        l = 2*i
        r = 2*i+1
        IF (l <= nelem) THEN
          IF (heap(l) % i1 < heap(i) % i1) THEN
            heaporder = .FALSE.
            write (*,*) 'left: ', l, i
            EXIT
          END IF
        END IF
        IF (r <= nelem) THEN
          IF (heap(r) % i1 < heap(i) % i1) THEN
            heaporder = .FALSE.
            write (*,*) 'right: ', r, i
            EXIT
          END IF
        END IF
      END DO
    END FUNCTION BinaryHeapIsHeap

  END SUBROUTINE ElmerMeshToDualGraph

  SUBROUTINE Graph_Deallocate(Graph)
    IMPLICIT NONE
    TYPE(Graph_t) :: Graph

    DEALLOCATE(Graph % ptr)
    DEALLOCATE(Graph % ind)
    Graph % n = 0
  END SUBROUTINE Graph_Deallocate

  SUBROUTINE ElmerGraphColour(Graph, Colouring, ConsistentColours)
    IMPLICIT NONE

    TYPE(Graph_t), INTENT(IN) :: Graph
    TYPE(Graphcolour_t) :: Colouring
    LOGICAL, OPTIONAL :: ConsistentColours

    INTEGER, ALLOCATABLE :: uncolored(:)
    INTEGER, ALLOCATABLE :: fc(:), ucptr(:), rc(:), rcnew(:)

    INTEGER :: nc, dualmaxdeg, i, v, w, uci, wci, vli, vti, vcol, wcol, &
            nrc, nunc, nthr, TID, allocstat, gn
    INTEGER, POINTER :: colours(:)
    INTEGER, PARAMETER :: VERTEX_PER_THREAD = 100
    LOGICAL :: consistent

    ! Iterative parallel greedy algorithm (Alg 2.) from 
    ! U. V. Catalyurek, J. Feo, A.H. Gebremedhin, M. Halappanavar, A. Pothen. 
    ! "Graph coloring algorithms for multi-core and massively multithreaded systems".
    ! Parallel computing, 38, 2012, pp. 576--594. 

    ! Initialize number of colours, maximum degree of graph and number of 
    ! uncolored vertices
    nc = 0
    dualmaxdeg = 0
    gn = Graph % n
    nunc = gn

    ! Check if a reproducible colouring is being requested
    consistent = .FALSE.
    IF (PRESENT(ConsistentColours)) consistent = ConsistentColours

    ! Get maximum vertex degree of the given graph
    !$OMP PARALLEL DO SHARED(Graph) &
    !$OMP PRIVATE(v) REDUCTION(max:dualmaxdeg) DEFAULT(NONE)
    DO v=1,Graph % n
      dualmaxdeg = MAX(dualmaxdeg, Graph % ptr(v+1)- Graph % ptr(v))
    END DO
    !$OMP END PARALLEL DO

    nthr = 1
    ! Ensure that each vertex has at most one thread attached to it
    !$ IF (.NOT. consistent) nthr = MIN(omp_get_max_threads(), gn)

    ! Allocate memory for colours of vertices and thread colour pointers
    ALLOCATE(colours(gn), uncolored(gn), ucptr(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
            'Unable to allocate colour maps!')

    !$OMP PARALLEL SHARED(gn, dualmaxdeg, Graph, colours, nunc, &
    !$OMP                 uncolored, ucptr, nthr) &
    !$OMP PRIVATE(uci, vli, vti, v, w, wci, vcol, wcol, fc, nrc, rc, rcnew, &
    !$OMP         allocstat, TID) &
    !$OMP REDUCTION(max:nc) DEFAULT(NONE) NUM_THREADS(nthr)

    TID=1
    !$ TID=OMP_GET_THREAD_NUM()+1

    ! Greedy algorithm colours a given graph with at 
    ! most max_{v\in V} deg(v)+1 colours
    ALLOCATE(fc(dualmaxdeg+1), rc((gn/nthr)+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
            'Unable to allocate local workspace!')
    ! Initialize forbidden colour array (local to thread)
    fc = 0

    ! Initialize colours and uncolored entries
    !$OMP DO 
    DO v=1,gn
      colours(v)=0
      ! U <- V
      uncolored(v)=v
    END DO
    !$OMP END DO

    DO
      ! For each v\in U in parallel do
      !$OMP DO
      DO uci=1,nunc
        v = uncolored(uci)
        vli = Graph % ptr(v)
        vti = Graph % ptr(v+1)-1

        ! For each w\in adj(v) do
        DO w=vli, vti
          ! fc[colour[w]]<-v
          !$OMP ATOMIC READ
          wcol = colours(Graph % ind(w))
          IF (wcol /= 0) fc(wcol) = v
        END DO

        ! Find smallest permissible colour for vertex
        ! c <- min\{i>0: fc[i]/=v \}
        DO i=1,dualmaxdeg+1
          IF (fc(i) /= v) THEN
            !$OMP ATOMIC WRITE 
            colours(v) = i
            ! Maintain maximum colour
            nc = MAX(nc, i)
            EXIT
          END IF
        END DO
      END DO
      !$OMP END DO

      nrc = 0
      ! For each v\in U in parallel do
      !$OMP DO
      DO uci=1,nunc
        v = uncolored(uci)
        vli = Graph % ptr(v)
        vti = Graph % ptr(v+1)-1
        vcol = colours(v)

        ! Make sure that recolour array has enough storage for 
        ! the worst case (all elements need to be added)
        IF (SIZE(rc)<nrc+(vti-vli)+1) THEN
          ALLOCATE(rcnew(MAX(SIZE(rc)*2, nrc+(vti-vli)+1)), STAT=allocstat)
          IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
                  'Unable to allocate local workspace!')
          rcnew(1:nrc)=rc(1:nrc)
          DEALLOCATE(rc)
          CALL MOVE_ALLOC(rcnew, rc)
        END IF

        ! For each w\in adj(v) do
        DO wci=vli,vti
          w = Graph % ind(wci)
          IF (colours(w)==vcol .AND. v>w) THEN
            ! R <- R\bigcup {v} (thread local)
            nrc = nrc + 1
            rc(nrc)=v
            EXIT
          END IF
        END DO
      END DO
      !$OMP END DO NOWAIT

      ucptr(TID)=nrc
      !$OMP BARRIER

      !$OMP SINGLE
      CALL ComputeCRSIndexes(nthr, ucptr)
      nunc = ucptr(nthr+1)-1
      !$OMP END SINGLE

      ! U <- R
      uncolored(ucptr(TID):ucptr(TID+1)-1)=rc(1:nrc)
      !$OMP BARRIER

      ! Colour the remaining vertices sequentially if the 
      ! size of the set of uncoloured vertices is small enough
      IF (nunc < nthr*VERTEX_PER_THREAD) THEN
        !$OMP SINGLE
        DO uci=1,nunc
          v = uncolored(uci)
          vli = Graph % ptr(v)
          vti = Graph % ptr(v+1)-1

          ! For each w\in adj(v) do
          DO w=vli, vti
            ! fc[colour[w]]<-v
            wcol = colours(Graph % ind(w))
            IF (wcol /= 0) fc(wcol) = v
          END DO

          ! Find smallest permissible colour for vertex
          ! c <- min\{i>0: fc[i]/=v \}
          DO i=1,dualmaxdeg+1
            IF (fc(i) /= v) THEN
              ! Single thread, no collisions possible 
              colours(v) = i
              ! Maintain maximum colour
              nc = MAX(nc, i)
              EXIT
            END IF
          END DO
        END DO
        !$OMP END SINGLE NOWAIT

        EXIT
      END IF

    END DO

    ! Deallocate thread local storage
    DEALLOCATE(fc, rc)
    !$OMP END PARALLEL

    DEALLOCATE(uncolored, ucptr)

    ! Set up colouring data structure
    Colouring % nc = nc
!   CALL MOVE_ALLOC(colours, Colouring % colours)
    Colouring % colours => colours
  END SUBROUTINE ElmerGraphColour

  SUBROUTINE Colouring_Deallocate(Colours)
    IMPLICIT NONE
    TYPE(GraphColour_t) :: Colours

    DEALLOCATE(Colours % colours)
    Colours % nc = 0
  END SUBROUTINE Colouring_Deallocate

  SUBROUTINE ElmerColouringToGraph(Colours, PackedList)
    IMPLICIT NONE

    TYPE(GraphColour_t), INTENT(IN) :: Colours
    TYPE(Graph_t) :: PackedList

    INTEGER, ALLOCATABLE :: cptr(:), cind(:)

    INTEGER :: nc, c, i, n, allocstat

    nc = Colours % nc
    n = size(Colours % colours)
    ALLOCATE(cptr(nc+1), cind(n), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerGatherColourLists','Memory allocation failed.')
    cptr = 0
    ! Count number of elements in each colour
    DO i=1,n
      cptr(Colours % colours(i))=cptr(Colours % colours(i))+1
    END DO

    CALL ComputeCRSIndexes(nc, cptr)

    DO i=1,n
      c=Colours % colours(i)
      cind(cptr(c))=i
      cptr(c)=cptr(c)+1
    END DO

    DO i=nc,2,-1
      cptr(i)=cptr(i-1)
    END DO
    cptr(1)=1

    ! Set up graph data structure
    PackedList % n = nc
    CALL MOVE_ALLOC(cptr, PackedList % ptr)
    CALL MOVE_ALLOC(cind, PackedList % ind)
  END SUBROUTINE ElmerColouringToGraph

  ! Routine constructs colouring for boundary mesh based on colours of main mesh
  SUBROUTINE ElmerBoundaryGraphColour(Mesh, Colours, BoundaryColours)
    IMPLICIT NONE

    TYPE(Mesh_t), INTENT(IN) :: Mesh
    TYPE(GraphColour_t), INTENT(IN) :: Colours
    TYPE(GraphColour_t) :: BoundaryColours

    TYPE(Element_t), POINTER :: Element
    INTEGER :: elem, nelem, nbelem, astat, lcolour, rcolour, nbc
    INTEGER, POINTER :: bcolours(:)

    nelem = Mesh % NumberOfBulkElements
    nbelem = Mesh % NumberOfBoundaryElements

    ! Allocate boundary colouring
    ALLOCATE(bcolours(nbelem), STAT=astat)
    IF (astat /= 0) THEN
       CALL Fatal('ElmerBoundaryGraphColour','Unable to allocate boundary colouring')
    END IF
    
    nbc = 0
    ! Loop over boundary mesh
    !$OMP PARALLEL DO &
    !$OMP SHARED(Mesh, nelem, nbelem, Colours, bcolours) &
    !$OMP PRIVATE(Element, lcolour, rcolour) &
    !$OMP REDUCTION(max:nbc) &
    !$OMP DEFAULT(NONE)
    DO elem=1,nbelem       
       Element => Mesh % Elements(nelem+elem)

       ! Try to find colour for boundary element based on left / right parent
       lcolour = 0
       IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
          lcolour = Colours % colours(Element % BoundaryInfo % Left % ElementIndex)
       END IF
       rcolour = 0
       IF (ASSOCIATED(Element % BoundaryInfo % Right)) THEN
          rcolour = Colours % colours(Element % BoundaryInfo % Right % ElementIndex)
       END IF

       ! Sanity check for debug
       IF (ASSOCIATED(Element % BoundaryInfo % Left) .AND. & 
          ASSOCIATED(Element % BoundaryInfo % Right) .AND. &
            lcolour /= rcolour) THEN
         CALL Warn('ElmerBoundaryGraphColour','Inconsistent colours for boundary element: ' &
               // i2s(elem) // "=>" &
               // i2s(lcolour)// " | "//i2s(rcolour))
         WRITE (*,*) Element % BoundaryInfo % Left % ElementIndex, Element % BoundaryInfo % Right % ElementIndex
       END IF

       bcolours(elem)=MAX(lcolour,rcolour)
       nbc=MAX(nbc,bcolours(elem))
    END DO
    !$OMP END PARALLEL DO

    ! Set up colouring data structure
    BoundaryColours % nc = nbc
!   CALL MOVE_ALLOC(bcolours, BoundaryColours % colours)
    BoundaryColours % colours => bcolours
  END SUBROUTINE ElmerBoundaryGraphColour
  
  ! Given CRS indices, referenced indirectly from graph, 
  ! evenly load balance the work among the nthr threads
  SUBROUTINE ThreadLoadBalanceElementNeighbour(nthr, gn, gptr, gind, &
          rptr, blkleads)
    IMPLICIT NONE

    INTEGER :: nthr
    INTEGER, INTENT(IN) :: gn
    INTEGER :: gptr(:), gind(:), rptr(:)
    INTEGER, ALLOCATABLE :: blkleads(:)

    INTEGER :: i, j, k, wrk, gwrk, thrwrk, allocstat

    ! Compute number of nonzeroes / thread
    !$ nthr = MIN(nthr,gn)

    ALLOCATE(blkleads(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ThreadLoadBalanceElementNeighbour', &
            'Unable to allocate blkleads!')

    ! Special case of just one thread
    IF (nthr == 1) THEN
      blkleads(1)=1
      blkleads(2)=gn+1
      RETURN
    END IF

    ! Compute total global work
    gwrk = 0
    DO i=1,gn
      DO j=gptr(i),gptr(i+1)-1
        gwrk = gwrk + (rptr(gind(j)+1)-rptr(gind(j)))
      END DO
    END DO

    ! Amount of work per thread
    thrwrk = CEILING(REAL(gwrk,dp) / nthr)

    ! Find rows for each thread to compute
    blkleads(1)=1
    DO i=1,nthr
      wrk = 0
      ! Acquire enough work for thread i
      DO j=blkleads(i),gn
        DO k=gptr(j),gptr(j+1)-1
          wrk = wrk + (rptr(gind(j)+1)-rptr(gind(j)))
        END DO
        IF (wrk >= thrwrk) EXIT
      END DO

      blkleads(i+1)=j+1
      ! Check if we have run out of rows
      IF (j+1>gn) EXIT
    END DO
    ! Reset number of rows (may be less than or equal to original number)
    nthr = i
    ! Assign what is left of the matrix to the final thread
    blkleads(nthr+1)=gn+1
  END SUBROUTINE ThreadLoadBalanceElementNeighbour

  SUBROUTINE ThreadStaticWorkShare(nthr, gn, blkleads)
    IMPLICIT NONE

    INTEGER :: nthr
    INTEGER, INTENT(IN) :: gn
    INTEGER, ALLOCATABLE :: blkleads(:)

    INTEGER :: i, rem, thrwrk, allocstat
    INTEGER :: totelem

    ! Compute number of nonzeroes / thread
    !$ nthr = MIN(nthr,gn)

    ALLOCATE(blkleads(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ThreadStaticWorkShare', &
            'Unable to allocate blkleads!')

    ! Special case of just one thread
    IF (nthr == 1) THEN
      blkleads(1)=1
      blkleads(2)=gn+1
      RETURN
    END IF

    ! Assuming even distribution of nodes / element, 
    ! distribute rows for each thread to compute 
    blkleads(1)=1
    thrwrk = gn / nthr
    rem = gn-nthr*thrwrk
    ! totelem = 0
    DO i=1,nthr-1
      IF (i<rem) THEN
        blkleads(i+1)=blkleads(i)+thrwrk+1
      ELSE
        blkleads(i+1)=blkleads(i)+thrwrk
      END IF
    END DO
    ! Assign what is left of the matrix to the final thread
    blkleads(nthr+1)=gn+1
  END SUBROUTINE ThreadStaticWorkShare

  ! Given row counts, in-place compute CRS indices to data
  SUBROUTINE ComputeCRSIndexes(n, arr)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER :: arr(:)

    INTEGER :: i, indi, indip

    indi = arr(1)
    arr(1)=1
    DO i=1,n-1
      indip=arr(i+1)
      arr(i+1)=arr(i)+indi
      indi=indip
    END DO
    arr(n+1)=arr(n)+indi
  END SUBROUTINE ComputeCRSIndexes


END MODULE MeshGraph
