MODULE LocalTypes
#ifdef _OPENMP
    USE omp_lib
#endif
    USE Types
    USE GeneralUtils, ONLY: Sort

    IMPLICIT NONE

    TYPE ElementCache_t
        TYPE(ElementType_t), POINTER :: Type => NULL()
        INTEGER :: nc = 0
        REAL(Kind=dp), POINTER :: U(:)=>NULL(),V(:)=>NULL(),W(:)=>NULL()
    END TYPE ElementCache_t

    TYPE IntegerList_t
        INTEGER :: nelem = 0
        INTEGER, ALLOCATABLE :: entries(:)
    END TYPE IntegerList_t

    TYPE IntegerHashBucket_t
        ! Bucket elements
        TYPE(IntegerList_t), ALLOCATABLE :: list
    END TYPE IntegerHashBucket_t

    TYPE IntegerHashSet_t
        TYPE(IntegerHashBucket_t), ALLOCATABLE :: set(:)
        TYPE(IntegerList_t), ALLOCATABLE :: entries
        REAL(KIND=dp) :: fratio
    END TYPE IntegerHashSet_t

    TYPE VertexMap_t
#ifdef _OPENMP
        INTEGER(KIND=OMP_LOCK_KIND), ALLOCATABLE :: vlock(:)
#endif
        TYPE(IntegerHashSet_t), ALLOCATABLE :: map(:)
    END TYPE VertexMap_t

    INTEGER, PARAMETER :: INTEGERLIST_DEFAULT_SIZE = 64
    INTEGER, PARAMETER :: INTEGERHASHSET_DEFAULT_SIZE = 64
    INTEGER, PARAMETER :: INTEGERHASHSET_CHAIN_DEFAULT_SIZE = 8
    REAL(KIND=dp), PARAMETER :: INTEGERHASHSET_FILLRATIO = REAL(0.75, dp)
    INTEGER, PARAMETER :: HEAPALG_THRESHOLD = 12
CONTAINS

    SUBROUTINE ElmerMeshToDualGraph(Mesh, n, dualptr, dualind)
        IMPLICIT NONE

        TYPE(Mesh_t) :: Mesh
        INTEGER, INTENT(OUT) :: n
        INTEGER, ALLOCATABLE :: dualptr(:), dualind(:)

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
        ! HEAP (list) MERGE: heap and mask for merge
        INTEGER, ALLOCATABLE :: wrkheap(:)
        LOGICAL, ALLOCATABLE :: wrkmask(:)

        TYPE :: IntTuple_t
          INTEGER :: i1, i2
        END type IntTuple_t
        TYPE(IntTuple_t), ALLOCATABLE :: wrkheap2(:)
        ! INTEGER, ALLOCATABLE :: wrkheap2(:), wrkheap3(:)

        ! OpenMP thread block leads for work division
        INTEGER, ALLOCATABLE :: thrblk(:)
        ! Work indices
        INTEGER, ALLOCATABLE :: wrkind(:), wrkindresize(:)
        INTEGER :: nwrkind

        ! Variables
        INTEGER :: i, dnnz, eid, nl, nli, nti, nn, nv, nthr, &
                   te, thrli, thrti, vli, vti, TID, allocstat
        INTEGER :: mapSizePad, maxNodesPad, neighSizePad

#ifdef HAVE_TIMING
        REAL(kind=dp) :: t_start, t_end
#endif

        ! Mesh data
        Elements => Mesh % Elements
        nvertex = Mesh % NumberOfNodes
        nelem = Mesh % NumberOfBulkElements
        
        ! Initialize dual mesh size and number of nonzeroes
        n = nelem
        dnnz = 0

#ifdef HAVE_TIMING
        t_start = ftimer()
#endif
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
#ifdef HAVE_TIMING
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'Mesh transformation, ElmerMeshToDualGraph: ', t_end - t_start, ' sec.'
#endif

#ifdef HAVE_TIMING
        t_start = ftimer()
#endif
        ! Construct vertex to element list (in serial!)
        CALL VertexToElementList(nelem, nvertex, eptr, eind, vptr, vind)

        ! Allocate pointers to dual mesh
        ALLOCATE(dualptr(nelem+1), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
                                       'Unable to allocate dual mesh!')

        ! Divide work by number of rows in the vertex graph
        nthr = 1 
        !$ nthr = omp_get_max_threads()

        ! Load balance the actual work done by threads
        ! CALL ThreadLoadBalanceElementNeighbour(nthr, nelem, eptr, eind, vptr, thrblk)
        CALL ThreadStaticWorkShare(nthr, nelem, thrblk)

        ! TODO: Begin OpenMP parallel section here
        !$OMP PARALLEL SHARED(nelem, nvertex, eptr, eind, &
        !$OMP                 vptr, vind, Mesh, dualind, dualptr, &
        !$OMP                 nthr, thrblk, dnnz) &
        !$OMP PRIVATE(i, eid, nli, nti, nn, nv, vli, vti, te, &
        !$OMP         maxNodesPad, neighSizePad, ptrli, ptrti, &
        !$OMP         wrkheap, wrkheap2, wrkmask, wrkmap, neighind, &
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
        ! IF (.TRUE.) THEN
        ! IF (.FALSE.) THEN
            ! With multiple threads, use heap based merge
            ALLOCATE(wrkheap(maxNodesPad), &
                     wrkmask(maxNodesPad), &
                     wrkheap2(maxNodesPad), STAT=allocstat)
            IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
                                           'Unable to allocate local workspace!')
        ELSE
            ! With a small number of threads, use map -based merge
            mapSizePad = IntegerNBytePad(Mesh % NumberOfBulkElements, 8)
            ALLOCATE(wrkmap(Mesh % NumberOfBulkElements), STAT=allocstat)
            IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
                                       'Unable to allocate local workspace!')
            ! Initialize map
            wrkmap(:)=0
        END IF

        ! Allocate local list for results
        ! TODO: Compute multiplier along with thread work division
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
            ! IF (.TRUE.) THEN
            ! IF (.FALSE.) THEN
            ! CALL kWayMergeList(eid, nv, ptrli, ptrti, &
            !         te, vind, nn, neighind, wrkheap, wrkmask)
            ! CALL kWayMergeList3(eid, nv, ptrli, ptrti, &
            !         te, vind, nn, neighind, wrkheap, wrkmask)
            ! CALL kWayMergeList2(eid, nv, ptrli, ptrti, &
            !                     te, vind, nn, neighind, wrkheap, &
            !                     wrkheapval, wrkmask)
            CALL kWayMergeHeap(eid, nv, ptrli, ptrti, &
                                te, vind, nn, neighind, wrkheap2)
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
          dualptr(eid)=nn
        END DO

        ! Get the global size of the dual mesh
        !$OMP DO REDUCTION(+:dnnz)
        DO i=1,nthr
          dnnz = nwrkind
        END DO
        !$OMP END DO

        ! Allocate memory for dual mesh indices
        !$OMP SINGLE
        ALLOCATE(dualind(dnnz), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
                                       'Unable to allocate dual mesh!')
        ! dualptr stores row counts, build crs pointers from them
        CALL ComputeCRSIndexes(nelem, dualptr)
        !$OMP END SINGLE

        dualind(dualptr(thrli):dualptr(thrti)-1)=wrkind(1:nwrkind)
        
        IF (nthr >= HEAPALG_THRESHOLD) THEN
        ! IF (.TRUE.) THEN
        ! IF (.FALSE.) THEN
          DEALLOCATE(wrkheap, wrkmask, wrkheap2, STAT=allocstat)
        ELSE
          DEALLOCATE(wrkmap, STAT=allocstat)
        END IF
        IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
                                       'Unable to deallocate local workspace!')
        DEALLOCATE(neighind, ptrli, ptrti, wrkind)

        !$OMP END PARALLEL

#ifdef HAVE_TIMING
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'Dual graph creation, ElmerMeshToDualGraph: ', t_end - t_start, ' sec.'
#endif

        ! Deallocate the rest of memory
        DEALLOCATE(eind, eptr, vptr, vind, thrblk)

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
                ! tmpi = vptr(1)
                ! vptr(1)=1
                ! DO i=1,nvertex-1
                !     tmpip=vptr(i+1)
                !     vptr(i+1)=vptr(i)+tmpi
                !     tmpi=tmpip
                ! END DO
                ! vptr(nvertex+1)=vptr(nvertex)+tmpi

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
            SUBROUTINE kWayMergeArray(node, nv, ptrli, ptrti, te, vind, nn, neighind, map)
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
                ! nn=j-1
            END SUBROUTINE kWayMergeArray

            ! k-way merge with a list
            SUBROUTINE kWayMergeList(node, nv, ptrli, ptrti, te, vind, &
                  nn, neighind, list, mask)
                IMPLICIT NONE

                INTEGER, INTENT(IN) :: node, nv
                INTEGER :: ptrli(:)
                INTEGER, INTENT(IN) ::ptrti(:), te
                INTEGER, INTENT(IN) :: vind(:)
                INTEGER, INTENT(OUT) :: nn
                INTEGER :: neighind(:)
                INTEGER :: list(:)
                LOGICAL :: mask(:)

                ! Local variables
                INTEGER :: i, j, k,  wrki, vindi, elem
                INTEGER :: pind

                DO i=1,nv
                    list(i)=vind(ptrli(i))
                END DO

                pind=-1
                mask(1:nv) = .TRUE.

                nn = 1
                DO elem=1,te
                    i=1
                    DO j=2,nv
                        IF (mask(j) .AND. list(j)<list(i)) i=j
                    END DO
                    ! i now contains the index of the minimum entry
                    
                    ! Do not add node if it was previously added to the list
                    vindi = list(i)
                    IF (vindi /= node .AND. pind /= vindi) THEN
                        neighind(nn)=vindi
                        nn = nn + 1
                        pind = vindi
                    END IF

                    ! Advance row pointer
                    ptrli(i)=ptrli(i)+1
                    ! Update mask
                    IF (ptrli(i) < ptrti(i)) THEN
                        list(i)=vind(ptrli(i))
                    ELSE
                        mask(i)=.FALSE.
                        list(i)=HUGE(list(i))
                    END IF
                END DO
                nn=nn-1

            END SUBROUTINE kWayMergeList

            ! k-way merge with a list
            SUBROUTINE kWayMergeList2(node, nv, ptrli, ptrti, te, vind, &
                                      nn, neighind, indlist, vallist, mask)
                IMPLICIT NONE

                INTEGER, INTENT(IN) :: node, nv
                INTEGER :: ptrli(:)
                INTEGER, INTENT(IN) ::ptrti(:), te
                INTEGER, INTENT(IN) :: vind(:)
                INTEGER, INTENT(OUT) :: nn
                INTEGER :: neighind(:)
                INTEGER :: indlist(:)
                INTEGER :: vallist(:)
                LOGICAL :: mask(:)

                ! Local variables
                INTEGER :: i, j, k,  wrki, vindi, elem, nzl
                INTEGER :: pind, mind

                DO i=1,nv
                    indlist(i)=i
                END DO
                DO i=1,nv
                  vallist(i)=vind(ptrli(i))
                END DO
                nzl = nv

                pind = -1

                nn = 1
                DO WHILE(nzl>0)
                    i=1
                    DO j=2,nzl
                        IF (vallist(j)<vallist(i)) i=j
                    END DO
                    ! i now contains the index of the minimum entry
                    ! Do not add node if it was previously added to the list
                    vindi = vallist(i)
                    IF (vindi /= node .AND. pind /= vindi) THEN
                      neighind(nn)=vindi
                      pind=vindi
                      nn=nn+1
                    END IF

                    ! Advance row pointer
                    ptrli(indlist(i))=ptrli(indlist(i))+1
                    ! Check if all elements have been added
                    IF (ptrli(i)< ptrti(i)) THEN
                      vallist(i)=vind(ptrli(indlist(i)))
                    ELSE
                      ! Remove index i from lists
                      DO j=1,i-1
                        indlist(j)=indlist(j)
                      END DO
                      DO j=i+1,nzl
                        indlist(j-1)=indlist(j)
                      END DO
                      nzl=nzl-1
                    END IF
                END DO
                nn=nn-1
            END SUBROUTINE kWayMergeList2

            ! k-way merge with a list
            SUBROUTINE kWayMergeList3(node, nv, ptrli, ptrti, te, vind, &
                  nn, neighind, list, mask)
                IMPLICIT NONE

                INTEGER, INTENT(IN) :: node, nv
                INTEGER :: ptrli(:)
                INTEGER :: ptrti(:), te
                INTEGER, INTENT(IN) :: vind(:)
                INTEGER, INTENT(OUT) :: nn
                INTEGER :: neighind(:)
                INTEGER :: list(:)
                LOGICAL :: mask(:)

                ! Local variables
                INTEGER :: i, j, k, wrki, vindi, elem, nl
                INTEGER :: pind

                ! DO i=1,nv
                !     list(i)=vind(ptrli(i))
                ! END DO

                nl = nv
                pind=-1
                nn = 1
                DO elem=1,te
                    i=1
                    DO j=2,nl
                        ! IF (list(j)<list(i)) i=j
                      IF (vind(ptrli(j)) < vind(ptrli(i))) i=j
                    END DO
                    ! i now contains the index of the minimum entry
                    vindi = vind(ptrli(i))
                    ! Do not add node if it was previously added to the list
                    ! Advance row pointer
                    ptrli(i)=ptrli(i)+1
                    
                    IF (vindi /= node .AND. pind /= vindi) THEN
                        neighind(nn)=vindi
                        nn = nn + 1
                        pind = vindi
                    END IF
                    ! Update lists
                    IF (ptrli(i) == ptrti(i)) THEN
                      ! Remove ith element from list and pointers
                      DO j=i,nl-1
                        ! list(j)=list(j+1)
                        ptrli(j)=ptrli(j+1)
                        ptrti(j)=ptrti(j+1)
                      END DO
                      nl=nl-1
                    END IF
                END DO
                nn=nn-1

            END SUBROUTINE kWayMergeList3
            
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

                ! IF (.NOT. BinaryHeapIsHeap(heap, nzheap)) THEN
                !   WRITE (*,*) heap(1:nzheap) % i1
                !   WRITE (*,*) nzheap
                !   WRITE (*,*) 'Not a heap, build!'
                !   STOP
                ! END IF

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

                  ! IF (.NOT. BinaryHeapIsHeap(heap, nzheap)) THEN
                  !   WRITE (*,*) heap(1:nzheap) % i1
                  !   WRITE (*,*) nzheap
                  !   WRITE (*,*) 'Not a heap, insert/delete!'
                  !   STOP
                  ! END IF
  
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

    SUBROUTINE ElmerGraphColour(gn, gptr, gind, nc, colours)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: gn
        INTEGER, INTENT(IN) :: gptr(:), gind(:)
        INTEGER :: nc
        INTEGER, ALLOCATABLE :: colours(:)

        INTEGER, ALLOCATABLE :: uncolored(:)
        INTEGER, ALLOCATABLE :: fc(:), ucptr(:), rc(:), rcnew(:)

        INTEGER :: dualmaxdeg, i, v, w, uci, wci, vli, vti, vcol, wcol, &
                   nrc, nunc, nthr, TID, allocstat
        INTEGER, PARAMETER :: VERTEX_PER_THREAD = 100

        ! Iterative parallel greedy algorithm (Alg 2.) from 
        ! U. V. Catalyurek, J. Feo, A.H. Gebremedhin, M. Halappanavar, A. Pothen. 
        ! "Graph coloring algorithms for multi-core and massively multithreaded systems".
        ! Parallel computing, 38, 2012, pp. 576--594. 

        ! Initialize number of colours, maximum degree of graph and number of 
        ! uncolored vertices
        nc = 0
        dualmaxdeg = 0
        nunc = gn
                
        ! Get maximum vertex degree of the given graph
        !$OMP PARALLEL DO SHARED(gptr, gn) &
        !$OMP PRIVATE(v) REDUCTION(max:dualmaxdeg) DEFAULT(NONE)
        DO v=1,gn
            dualmaxdeg = MAX(dualmaxdeg, gptr(v+1)-gptr(v))
        END DO
        !$OMP END PARALLEL DO
        dualmaxdeg = dualmaxdeg + 1

        nthr = 1
        ! Ensure that each vertex has at most one thread attached to it
        !$ nthr = MIN(omp_get_max_threads(), gn)

        ! Allocate memory for colours of vertices and thread colour pointers
        ALLOCATE(colours(gn), uncolored(gn), ucptr(nthr+1), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
                                       'Unable to allocate colour maps!')
        
        !$OMP PARALLEL SHARED(gn, dualmaxdeg, gptr, gind, colours, nunc, &
        !$OMP                 uncolored, ucptr, nthr) &
        !$OMP PRIVATE(uci, vli, vti, v, w, wci, vcol, wcol, fc, nrc, rc, rcnew, &
        !$OMP         allocstat, TID) &
        !$OMP REDUCTION(max:nc) DEFAULT(NONE) NUM_THREADS(nthr)

        TID=1
        !$ TID=OMP_GET_THREAD_NUM()+1

        ! Greedy algorithm colours a given graph with at 
        ! most max_{v\in V} deg(v)+1 colours
        ALLOCATE(fc(dualmaxdeg), rc(gn/nthr), STAT=allocstat)
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
                vli = gptr(v)
                vti = gptr(v+1)-1

                ! For each w\in adj(v) do
                DO w=vli, vti
                    ! fc[colour[w]]<-v
                    !$OMP ATOMIC READ
                    wcol = colours(gind(w))
                    IF (wcol /= 0) fc(wcol) = v
                END DO

                ! Find smallest permissible colour for vertex
                ! c <- min\{i>0: fc[i]/=v \}
                DO i=1,dualmaxdeg
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
                vli = gptr(v)
                vti = gptr(v+1)-1
                vcol = colours(v)

                ! Make sure that recolour array has enough storage for 
                ! the worst case (all elements need to be added)
                IF (SIZE(rc)<nrc+(vti-vli)) THEN
                    ALLOCATE(rcnew(MAX(SIZE(rc)*2, nrc+(vti-vli))), STAT=allocstat)
                    IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
                                                   'Unable to allocate local workspace!')
                    rcnew(1:nrc)=rc(1:nrc)
                    DEALLOCATE(rc)
                    CALL MOVE_ALLOC(rcnew, rc)
                END IF

                ! For each w\in adj(v) do
                DO wci=vli,vti
                    w = gind(wci)
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
                    vli = gptr(v)
                    vti = gptr(v+1)-1
                    
                    ! For each w\in adj(v) do
                    DO w=vli, vti
                        ! fc[colour[w]]<-v
                        wcol = colours(gind(w))
                        IF (wcol /= 0) fc(wcol) = v
                    END DO

                    ! Find smallest permissible colour for vertex
                    ! c <- min\{i>0: fc[i]/=v \}
                    DO i=1,dualmaxdeg
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
    END SUBROUTINE ElmerGraphColour

    SUBROUTINE ElmerGatherColourLists(nc, colours, cptr, cind)
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nc
      INTEGER, INTENT(IN) :: colours(:)
      INTEGER, ALLOCATABLE :: cptr(:), cind(:)

      INTEGER :: c, i, n, allocstat

      n = size(colours)
      ALLOCATE(cptr(nc+1), cind(n), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('ElmerGatherColourLists','Memory allocation failed.')
      cptr = 0
      ! Count number of elements in each colour
      DO i=1,n
        cptr(colours(i))=cptr(colours(i))+1
      END DO

      CALL ComputeCRSIndexes(nc, cptr)

      DO i=1,n
        c=colours(i)
        cind(cptr(c))=i
        cptr(c)=cptr(c)+1
      END DO

      DO i=nc,2,-1
        cptr(i)=cptr(i-1)
      END DO
      cptr(1)=1
    END SUBROUTINE ElmerGatherColourLists

    SUBROUTINE ConstructVertexToElementList(ne, nn, eptr, eind, VertexToElementList)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ne, nn
        INTEGER :: eptr(:), eind(:)
        TYPE(IntegerList_t), ALLOCATABLE :: VertexToElementList(:)
        TYPE(Mesh_t) :: Mesh

        TYPE(Element_t), POINTER :: Element, Elements(:)
        INTEGER :: i, j, v, eli, eti, nelem, nvertex, allocstat
        INTEGER, ALLOCATABLE :: vptr(:), vind(:)
#ifdef HAVE_TIMING
        REAL(kind=dp) :: t_start, t_end
#endif
#ifdef _OPENMP
        INTEGER(kind=omp_lock_kind), ALLOCATABLE :: vlock(:)
#endif

        nelem = ne
        nvertex = nn

#ifdef HAVE_TIMING
        t_start = ftimer()
#endif
        ! ALLOCATE(vptr(nvertex), STAT=allocstat)
        ! IF (allocstat /= 0) CALL Fatal('ConstructVertexToElementMap', &
        !                                'Vertex pointer allocation failed!')
        ! vptr(:) = 0

#ifdef _OPENMP
        ! Allocate vertex locks
        ALLOCATE(vlock(nvertex), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('ConstructVertexToElementMap', &
              'Lock allocation failed!')
#endif
        ! Allocate vertex lists
        ALLOCATE(VertexToElementList(nvertex), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('ConstructVertexToElementMap', &
              'Vertex list allocation failed!')

#ifdef HAVE_TIMING
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'ConstructVertexToElementList, init: ', t_end - t_start, ' sec.'
#endif

        ! For each element
        !$OMP PARALLEL SHARED(nelem, nvertex, eind, eptr, VertexToElementList, vlock) &
        !$OMP PRIVATE(i, j, eli, eti, Element) DEFAULT(NONE)

        ! Initialize locks
#ifdef _OPENMP
        !$OMP DO
        DO i=1,nvertex
            CALL OMP_INIT_LOCK(vlock(i))
        END DO
        !$OMP END DO NOWAIT
#endif
        ! Initialize vertex lists
        !$OMP DO
        DO i=1,nvertex
            ! TODO: Change list size to be more dynamic
            CALL IntegerListInit(VertexToElementList(i), 32)
        END DO
        !$OMP END DO

        !$OMP DO
        ! For each element
        DO i=1,nelem
            eli = eptr(i)
            eti = eptr(i+1)-1

            ! For each vertex in element
            !DIR$ IVDEP
            DO j=eli, eti
                ! Add connection to vertex eind(j)
#ifdef _OPENMP
                CALL OMP_SET_LOCK(vlock(eind(j)))
#endif
                CALL IntegerListAdd(VertexToElementList(eind(j)), i)
#ifdef _OPENMP
                CALL OMP_UNSET_LOCK(vlock(eind(j)))
#endif
            END DO
        END DO
        !$OMP END DO

#ifdef _OPENMP
        !$OMP DO
        DO i=1,nvertex
            CALL OMP_DESTROY_LOCK(vlock(i))
        END DO
        !$OMP END DO NOWAIT
#endif

        !$OMP END PARALLEL

#ifdef _OPENMP
        ! Deallocate vertex locks
        DEALLOCATE(vlock)
#endif
    END SUBROUTINE ConstructVertexToElementList

     SUBROUTINE ConstructVertexToElementList2(nelem, nvertex, eptr, eind, vptr, vind)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: nelem, nvertex
        INTEGER :: eptr(:), eind(:)
        INTEGER, ALLOCATABLE :: vptr(:), vind(:)

        INTEGER :: i, j, v, eli, eti, ind, tmpi, tmpip, allocstat
#ifdef HAVE_TIMING
        REAL(kind=dp) :: t_start, t_end
#endif

#ifdef HAVE_TIMING
        t_start = ftimer()
#endif
        ! Initialize vertex structure (enough storage for nvertex vertices
        ! having eptr(nelem+1) elements)
        ALLOCATE(vptr(nvertex+1), vind(eptr(nelem+1)), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('ConstructVertexToElementMap', &
                                        'Vertex allocation failed!')
        vptr = 0

#ifdef HAVE_TIMING
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'ConstructVertexToElementList, init: ', t_end - t_start, ' sec.'
#endif

        ! For each element

        ! Compute number of elements attached to each vertex (size of lists)
        DO i=1,nelem
            eli = eptr(i)
            eti = eptr(i+1)-1

            DO j=eli, eti
                vptr(eind(j))=vptr(eind(j))+1
            END DO
        END DO
        ! WRITE (*,*) vptr(1:4), vptr(nvertex-2:nvertex+1)

        ! Compute cumulative sum (row pointers!)
        tmpi = vptr(1)
        vptr(1)=1
        DO i=1,nvertex-1
            tmpip=vptr(i+1)
            vptr(i+1)=vptr(i)+tmpi
            tmpi=tmpip
        END DO
        vptr(nvertex+1)=vptr(nvertex)+tmpi

        ! WRITE (*,*) vptr(1:4), vptr(nvertex-2:nvertex+1)

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
        ! WRITE (*,*) vptr(1:4), vptr(nvertex-2:nvertex+1)

        ! Correct row pointers
        DO i=nvertex,2,-1
            vptr(i)=vptr(i-1)
        END DO
        vptr(1)=1
        ! WRITE (*,*) vptr(1:4), vptr(nvertex-2:nvertex+1)
        ! STOP
    END SUBROUTINE ConstructVertexToElementList2

    ! Portable wall-clock timer
    FUNCTION ftimer() RESULT(timerval)
        IMPLICIT NONE

        REAL(KIND=dp) :: timerval
        INTEGER :: t, rate

#ifdef _OPENMP
        timerval = OMP_GET_WTIME()
#else
        CALL SYSTEM_CLOCK(t,count_rate=rate)
        timerval = REAL(t,KIND(dp))/REAL(rate,KIND(dp))
#endif
    END FUNCTION ftimer


    SUBROUTINE VertexMapInit(vmap, nvertex, mdim)
        IMPLICIT NONE

        TYPE(VertexMap_t) :: vmap
        INTEGER, INTENT(IN) :: nvertex
        INTEGER, INTENT(IN), OPTIONAL :: mdim
        INTEGER :: i, allocstat, meshdim, minitsize

        ! Mesh dimension, the default is 3D
        meshdim = 3
        IF (PRESENT(mdim)) meshdim = mdim

        ! Set default size of initial maps depending on
        ! the dimension of the mesh
        SELECT CASE (meshdim)
        CASE (1)
            minitsize = 2  ! NOTE: An absolute upper bound is 2
        CASE (2)
            minitsize = 16 ! For equilateral triangles 360/60=6
        CASE (3)
            minitsize = 64 ! Equilateral triangles in two planes =36
        CASE DEFAULT
            minitsize = 64
        END SELECT

        ALLOCATE(vmap % map(nvertex), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('VertexMapInit', &
              'Memory allocation failed!')

        !$OMP PARALLEL PRIVATE(i)

        !$OMP DO
        DO i=1,nvertex
            CALL IntegerHashSetInit(vmap % map(i), minitsize)
        END DO
        !$OMP END DO NOWAIT

#ifdef _OPENMP
        !$OMP SINGLE
        ! Allocate and initialize vertex locks
        ALLOCATE(vmap % vlock(nvertex), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('VertexMapInit', &
              'Memory allocation failed!')
        !$OMP END SINGLE

        !$OMP DO
        DO i=1,nvertex
            CALL OMP_INIT_LOCK(vmap % vlock(i))
        END DO
        !$OMP END DO NOWAIT
#endif

        !$OMP END PARALLEL
    END SUBROUTINE VertexMapInit

    SUBROUTINE VertexMapFromArray(vmap, n, vptr, vind)
        IMPLICIT NONE

        TYPE(VertexMap_t) :: vmap
        INTEGER, INTENT(IN) :: n, vptr(:), vind(:)

        INTEGER :: i, j, vli, vti, msize

        CALL VertexMapInit(vmap, n)

        ! For each vertex
        DO i=1,n
            vli=vptr(i)
            vti=vptr(i+1)-1

            ! Add mapping to array
            DO j=vli,vti
                CALL VertexMapAdd(vmap, i, vind(j))
            END DO
            msize = IntegerListGetSize(vmap % map(i) % entries)
            IF (msize /= vti-vli+1) THEN
                WRITE (*,*) 'WARNING: Some duplicate entries were found for row=', i
                ! WRITE (*,*) msize, vti-vli+1
                ! WRITE (*,*) vind(vli:vti)
                ! WRITE (*,*) vmap % map(i) % entries % entries(1:msize)
                ! STOP
            END IF
        END DO
    END SUBROUTINE VertexMapFromArray

    SUBROUTINE VertexMapAdd(vmap, vertexId, elementId)
        IMPLICIT NONE

        TYPE(VertexMap_t) :: vmap
        INTEGER, INTENT(IN) :: vertexId, elementId

!         IF (vertexId < 1 .OR. vertexId > SIZE(vmap % map)) RETURN

        ! Add vertex to map
        CALL IntegerHashSetAdd(vmap % map(vertexId), elementId)
    END SUBROUTINE VertexMapAdd

    SUBROUTINE VertexMapDelete(vmap, vertexId, elementId)
        IMPLICIT NONE

        TYPE(VertexMap_t) :: vmap
        INTEGER, INTENT(IN) :: vertexId, elementId

!        IF (vertexId < 1 .OR. vertexId > SIZE(vmap % map)) RETURN

        ! Delete vertex from map
        CALL IntegerHashSetDelete(vmap % map(vertexId), elementId)
    END SUBROUTINE VertexMapDelete

    SUBROUTINE VertexMapAddAtomic(vmap, vertexId, elementId)
        IMPLICIT NONE

        TYPE(VertexMap_t) :: vmap
        INTEGER, INTENT(IN) :: vertexId, elementId

!       IF (vertexId < 1 .OR. vertexId > SIZE(vmap % map)) RETURN

#ifdef _OPENMP
        ! Lock the vertex
        CALL OMP_SET_LOCK(vmap % vlock(vertexId))
#endif

        ! Add vertex to map
        CALL IntegerHashSetAdd(vmap % map(vertexId), elementId)

#ifdef _OPENMP
        ! Unlock the vertex
        CALL OMP_UNSET_LOCK(vmap % vlock(vertexId))
#endif
    END SUBROUTINE VertexMapAddAtomic

    SUBROUTINE VertexMapDeleteAtomic(vmap, vertexId, elementId)
        IMPLICIT NONE

        TYPE(VertexMap_t) :: vmap
        INTEGER, INTENT(IN) :: vertexId, elementId

!        IF (vertexId < 1 .OR. vertexId > SIZE(vmap % map)) RETURN

#ifdef _OPENMP
        ! Lock the vertex
        CALL OMP_SET_LOCK(vmap % vlock(vertexId))
#endif

        ! Delete vertex from map
        CALL IntegerHashSetDelete(vmap % map(vertexId), elementId)

#ifdef _OPENMP
        ! Unlock the vertex
        CALL OMP_UNSET_LOCK(vmap % vlock(vertexId))
#endif
    END SUBROUTINE VertexMapDeleteAtomic

    FUNCTION VertexMapFind(vmap, vertexId, conn) RESULT(found)
        IMPLICIT NONE
        TYPE(VertexMap_t), TARGET :: vmap
        INTEGER, INTENT(IN) :: vertexId, conn
        LOGICAL :: found

        found = .FALSE.
 !       IF (vertexId < 1 .OR. vertexId > SIZE(vmap % map)) RETURN
        found = IntegerHashSetFind(vmap % map(vertexid), conn)
    END FUNCTION VertexMapFind

    FUNCTION VertexMapGetList(vmap, vertexId) RESULT(nlist)
        IMPLICIT NONE
        TYPE(VertexMap_t), TARGET :: vmap
        INTEGER, INTENT(IN) :: vertexId

        TYPE(IntegerList_t), POINTER :: nlist

 !       IF (vertexId < 1 .OR. vertexId > SIZE(vmap % map)) RETURN

        nlist => vmap % map(vertexId) % entries
    END FUNCTION VertexMapGetList

    SUBROUTINE VertexMapToLists(vmap, nlist)
        IMPLICIT NONE

        TYPE(VertexMap_t) :: vmap
        TYPE(IntegerList_t), ALLOCATABLE :: nlist(:)

        INTEGER :: i, nvertex, allocstat
        TYPE(IntegerList_t) :: entries

        nvertex = SIZE(vmap % map)

        ALLOCATE(nlist(nvertex), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('VertexMapToList', &
              'Memory allocation failed!')

        !$OMP PARALLEL DO PRIVATE(i)
        DO i=1,nvertex
            ! Add all elements from map to a simple list
            CALL IntegerListInit(nlist(i), vmap % map(i) % entries % nelem)
            CALL IntegerListAddAll(nlist(i), vmap % map(i) % entries)
        END DO
        !$OMP END PARALLEL DO
    END SUBROUTINE VertexMapToLists

    SUBROUTINE VertexMapDeleteAll(vmap)
        IMPLICIT NONE

        TYPE(VertexMap_t) :: vmap

        INTEGER :: i, nvertex

        nvertex = SIZE(vmap % map)
        !$OMP PARALLEL PRIVATE(i)
#ifdef _OPENMP
        ! Deallocate vertex locks
        !$OMP DO
        DO i=1,nvertex
            CALL OMP_DESTROY_LOCK(vmap % vlock(i))
        END DO
        !$OMP END DO

        !$OMP SINGLE
        DEALLOCATE(vmap % vlock)
        !$OMP END SINGLE NOWAIT
#endif
        ! Delete contents of the vertex map
        !$OMP DO
        DO i=1,nvertex
            CALL IntegerHashSetDeleteAll(vmap % map(i))
        END DO
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        DEALLOCATE(vmap % map)
    END SUBROUTINE VertexMapDeleteAll

    FUNCTION VertexMapGetSize(vmap) RESULT(N)
        IMPLICIT NONE
        TYPE(VertexMap_t) :: vmap
        INTEGER :: n

        N = SIZE(vmap % map)
    END FUNCTION VertexMapGetSize

    FUNCTION VertexMapEquals(vmap1, vmap2) RESULT(equalTo)
        IMPLICIT NONE
        TYPE(VertexMap_t), TARGET :: vmap1, vmap2

        LOGICAL :: equalTo
        INTEGER :: i, id, n1, n2, ln1, ln2
        TYPE(IntegerList_t), POINTER :: l1, l2

        equalTo = .FALSE.

        n1 = VertexMapGetSize(vmap1)
        n2 = VertexMapGetSize(vmap2)

        ! Test size of maps, if not equal the maps are different
        IF (n1 /= n2) RETURN

        DO i=1, n1
            ! Test size of vertex mapping, if not equal the maps are different
            l1 => VertexMapGetList(vmap1, i)
            l2 => VertexMapGetList(vmap2, i)

            ln1 = IntegerListGetSize(l1)
            ln2 = IntegerListGetSize(l2)

            IF (ln1 /= ln2) RETURN

            ! Try to find each element of vmap1 from vmap2
            DO id = 1, ln1
                IF (.NOT. VertexMapFind(vmap2, i, IntegerListAt(l1, id))) RETURN
            END DO
        END DO

        equalTo = .TRUE.
    END FUNCTION VertexMapEquals

    SUBROUTINE VertexMapOutputString(vmap)
        IMPLICIT NONE
        TYPE(VertexMap_t) :: vmap

        INTEGER :: i

        WRITE (*,*) 'nvertex=', SIZE(vmap % map)
        DO i=1,SIZE(vmap % map)
            WRITE (*,*) 'nelem=', vmap % map(i) % entries % nelem
            WRITE (*,*) vmap % map(i) % entries % entries(1:vmap % map(i) % entries % nelem)
        END DO
    END SUBROUTINE VertexMapOutputString

    SUBROUTINE IntegerHashSetInit(iset, isize, fratio)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset
        INTEGER, INTENT(IN), OPTIONAL :: isize
        REAL(KIND=dp), INTENT(IN), OPTIONAL :: fratio

        INTEGER :: i, n, cn, allocstat

        n = INTEGERHASHSET_DEFAULT_SIZE
        IF (PRESENT(isize)) n = isize
        iset % fratio = INTEGERHASHSET_FILLRATIO
        IF (PRESENT(fratio) .AND. fratio > 0 &
              .AND. fratio <= REAL(1,dp)) iset % fratio = fratio

        IF (ALLOCATED(iset % set)) DEALLOCATE(iset % set)
        IF (ALLOCATED(iset % entries)) DEALLOCATE(iset % entries)

        ALLOCATE(iset % set(n), &
                 iset % entries, STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('IntegerHashSetInit', &
                                       'Memory allocation error!')

        CALL IntegerListInit(iset % entries, n)
    END SUBROUTINE IntegerHashSetInit

    SUBROUTINE IntegerHashSetAdd(iset, key)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset
        INTEGER, INTENT(IN) :: key

        INTEGER :: i, n, nelem, hkey, lind, allocstat
        REAL(KIND=dp) :: fratio
        INTEGER, ALLOCATABLE :: elementsold(:)

        ! Check for overfill and rehash if necessary
        n = SIZE(iset % set)
        fratio = iset % fratio
        IF ((iset % entries % nelem + 1) > FLOOR(iset % fratio * n)) THEN
            ! Store elements of old hash table
            nelem = iset % entries % nelem
            ALLOCATE(elementsold(nelem), STAT=allocstat)
            IF (allocstat /= 0) CALL Fatal('IntegerHashSetAdd', &
                                           'Memory allocation error!')
            elementsold(1:nelem) = iset % entries % entries(1:nelem)

            ! Reinitialize hashset (if pointers were used, we would need
            ! to call deleteall before the init)
            CALL IntegerHashSetInit(iset, 2*n, fratio)
            ! Rehash elements to new hashset
            DO i=1,nelem
                CALL IntegerHashSetAddHashEntry(iset, elementsold(i))
            END DO
            DEALLOCATE(elementsold)
        END IF

        ! Add the requested entry to the hashset
        CALL IntegerHashSetAddHashEntry(iset, key)
    END SUBROUTINE IntegerHashSetAdd

    SUBROUTINE IntegerHashSetAddHashEntry(iset, key)
        IMPLICIT NONE
        TYPE(IntegerHashSet_t) :: iset
        INTEGER, INTENT(IN) :: key

        INTEGER :: n, hkey, lind, allocstat

        hkey = IntegerHashSetHashFunction(iset, key)
        IF (.NOT. ALLOCATED(iset % set(hkey) % list)) THEN
            ! Add new chain and entry
            ALLOCATE(iset % set(hkey) % list, STAT=allocstat)
            IF (allocstat /= 0) CALL Fatal('IntegerHashSetAddHashEntry', &
                  'Memory allocation error!')

            CALL IntegerListInit(iset % set(hkey) % list, INTEGERHASHSET_CHAIN_DEFAULT_SIZE)

            CALL IntegerListAdd(iset % set(hkey) % list, key)
            CALL IntegerListAdd(iset % entries, key)
        ELSE
            ! Old chain, just add entry
            ASSOCIATE(chain => iset % set(hkey))
              lind = IntegerListFind(chain % list, key)

              ! A set can only contain one reference to each item
              IF (lind < 0) THEN
                  CALL IntegerListAdd(chain % list, key)
                  CALL IntegerListAdd(iset % entries, key)
              END IF
            END ASSOCIATE
        END IF

    END SUBROUTINE IntegerHashSetAddHashEntry

    SUBROUTINE IntegerHashSetDelete(iset, key)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset
        INTEGER, INTENT(IN) :: key

        INTEGER :: hkey, lind

        hkey = IntegerHashSetHashFunction(iset, key)
        IF (ALLOCATED(iset % set(hkey) % list)) THEN
            lind = IntegerListFind(iset % set(hkey) % list, key)
            IF (lind > 0) THEN
                ! Delete values from chain
                CALL IntegerListDeleteAt(iset % set(hkey) % list, lind)

                ! Delete entry from entry list (O(n) operation)
                CALL IntegerListDelete(iset % entries, key)
            END IF
        END IF
    END SUBROUTINE IntegerHashSetDelete

    SUBROUTINE IntegerHashSetDeleteAll(iset)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset
        INTEGER :: ent, hkey

        DO ent=1, iset % entries % nelem
            hkey = IntegerHashSetHashFunction(iset, iset % entries % entries(ent))
            ! Elements of chain may not be allocated if they have been
            ! previously deallocated by this routine
            IF (ALLOCATED(iset % set(hkey) % list)) THEN
                CALL IntegerListDeleteAll(iset % set(hkey) % list)
                DEALLOCATE(iset % set(hkey) % list)
            END IF
        END DO

        CALL IntegerListDeleteAll(iset % entries)
    END SUBROUTINE IntegerHashSetDeleteAll

    FUNCTION IntegerHashSetFind(iset, key) RESULT(found)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t), INTENT(IN) :: iset
        INTEGER, INTENT(IN) :: key
        LOGICAL :: found

        INTEGER :: hkey, ind

        found = .FALSE.
        hkey = IntegerHashSetHashFunction(iset, key)
        IF (ALLOCATED(iset % set(hkey) % list)) THEN
            ind = IntegerListFind(iset % set(hkey) % list, key)
            IF (ind > 0) found = .TRUE.
        END IF
    END FUNCTION IntegerHashSetFind

    FUNCTION IntegerHashSetHashFunction(iset, key) RESULT(hkey)
        IMPLICIT NONE
        TYPE(IntegerHashSet_t), INTENT(IN) :: iset
        INTEGER, INTENT(IN) :: key

        INTEGER :: m
        INTEGER :: hkey
        REAL(kind=dp), PARAMETER :: hA = (SQRT(REAL(5,dp))-1)/2

        ! Use a multiplicative hash function with Knuth's choice of
        ! A, i.e.,
        ! h(k) = floor(m*(k*hA mod 1)), with hA=(sqrt(5)-1)/2
        m = SIZE(iset % set)
        ! hkey = FLOOR(REAL(m,dp)*(key*hA-FLOOR(key*hA)))+1
        hkey = FLOOR(REAL(m,dp)*MOD(key*hA,REAL(1,dp)))+1
    END FUNCTION IntegerHashSetHashFunction

    SUBROUTINE IntegerHashSetOutputString(iset)
        IMPLICIT NONE
        TYPE(IntegerHashSet_t) :: iset

        INTEGER :: i

        WRITE (*,*) 'nelem=', iset % entries % nelem
        WRITE (*,*) iset % entries % entries(1:iset % entries % nelem)
        WRITE (*,*) 'HashSet contents'
        WRITE (*,*) 'SIZE(iset % set)=', SIZE(iset % set)
        WRITE (*,*) 'fill, nelem/SIZE(iset % set)=', REAL(iset % entries % nelem)/SIZE(iset % set)
        DO i=1,SIZE(iset % set)
            IF (.NOT. ALLOCATED(iset % set(i) % list)) THEN
                WRITE (*,'(A,I0,A)') 'set(',i,') empty'
            ELSE
                WRITE (*,'(A,I0,A)') 'set(',i,') allocated'
                WRITE (*,*) 'list of entries (keys)'
                WRITE (*,*) 'nelem=', iset % set(i) % list % nelem
                WRITE (*,*) iset % set(i) % list % entries(1:iset % set(i) % list % nelem)
            END IF
        END DO
    END SUBROUTINE IntegerHashSetOutputString

    SUBROUTINE IntegerListInit(ilist, isize)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, OPTIONAL, INTENT(IN) :: isize

        INTEGER :: n, allocstat

        n = INTEGERLIST_DEFAULT_SIZE
        IF (PRESENT(isize)) n = isize

        IF (ALLOCATED(ilist % entries)) DEALLOCATE(ilist % entries)
        ALLOCATE(ilist % entries(n), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('IntegerListInit', 'Memory allocation error!')

        ilist % nelem = 0
    END SUBROUTINE IntegerListInit

    FUNCTION IntegerListAt(ilist, ind) RESULT(key)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: ind
        INTEGER :: key

        key = 0
        ! IF (ind < 0 .OR. ind > SIZE(ilist % entries)) RETURN
        key = ilist % entries(ind)
    END FUNCTION IntegerListAt

    SUBROUTINE IntegerListAdd(ilist, item)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: item

        INTEGER, ALLOCATABLE :: elementsnew(:)
        INTEGER :: n, allocstat

        n = size(ilist % entries)
        IF (ilist % nelem + 1 > n) THEN
            ! Reallocate list structure with double the size
            ALLOCATE(elementsnew(2*n), STAT=allocstat)
            IF (allocstat /= 0) CALL Fatal('IntegerListAdd', 'Memory allocation error!')

            ! Copy elements, deallocate old element vector and move allocation
            elementsnew(1:n)=ilist % entries(1:n)
            DEALLOCATE(ilist % entries)
            CALL MOVE_ALLOC(elementsnew, ilist % entries)
        END IF

        ! Add new element
        ilist % nelem = ilist % nelem + 1
        ilist % entries(ilist % nelem)=item
    END SUBROUTINE IntegerListAdd

    SUBROUTINE IntegerListAddAll(ilist, alist)
      IMPLICIT NONE

      TYPE(IntegerList_t) :: ilist
      TYPE(IntegerList_t) :: alist

      INTEGER :: i, isize, ni, na, allocstat
      INTEGER, ALLOCATABLE :: elementsnew(:)

      isize = SIZE(ilist % entries)
      ni = ilist % nelem
      na = alist % nelem
      ! Check if reallocation of ilist is needed
      IF (isize < ni + na) THEN
          ! Reallocate a list structure with enough space
          ALLOCATE(elementsnew(2*isize+na), STAT=allocstat)
          IF (allocstat /= 0) CALL Fatal('IntegerListAddAll',&
                                         'Memory allocation error!')

          ! Copy elements, deallocate old element vector and move allocation
          elementsnew(1:ni)=ilist % entries(1:ni)
          DEALLOCATE(ilist % entries)
          CALL MOVE_ALLOC(elementsnew, ilist % entries)
      END IF

      ! ilist has space to hold all elements in alist
      ilist % entries(ni+1:ni+na) = alist % entries(1:na)
      ilist % nelem = ilist % nelem + na
    END SUBROUTINE IntegerListAddAll

    SUBROUTINE IntegerListAddArray(ilist, alist)
      IMPLICIT NONE

      TYPE(IntegerList_t) :: ilist
      INTEGER :: alist(:)

      INTEGER :: i, isize, ni, na, allocstat
      INTEGER, ALLOCATABLE :: elementsnew(:)

      isize = SIZE(ilist % entries)
      ni = ilist % nelem
      na = SIZE(alist)
      ! Check if reallocation of ilist is needed
      IF (isize < ni + na) THEN
          ! Reallocate a list structure with enough space
          ALLOCATE(elementsnew(2*isize+na), STAT=allocstat)
          IF (allocstat /= 0) CALL Fatal('IntegerListAddAll',&
                                         'Memory allocation error!')

          ! Copy elements, deallocate old element vector and move allocation
          elementsnew(1:ni)=ilist % entries(1:ni)
          DEALLOCATE(ilist % entries)
          CALL MOVE_ALLOC(elementsnew, ilist % entries)
      END IF

      ! ilist has space to hold all elements in alist
      ilist % entries(ni+1:ni+na) = alist(1:na)
      ilist % nelem = ilist % nelem + na
    END SUBROUTINE IntegerListAddArray

    SUBROUTINE IntegerListDelete(ilist, item)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: item

        INTEGER :: i, ind

        ind = IntegerListFind(ilist, item)
        CALL IntegerListDeleteAt(ilist, ind)
    END SUBROUTINE IntegerListDelete

    SUBROUTINE IntegerListDeleteAt(ilist, ind)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: ind

        INTEGER :: i

        ! IF (ind < 0 .OR. ind > SIZE(ilist % entries)) RETURN

        DO i=ind, ilist % nelem-1
            ilist % entries(i) = ilist % entries(i+1)
        END DO
        ilist % entries(ilist % nelem) = 0
        ilist % nelem = ilist % nelem - 1
    END SUBROUTINE IntegerListDeleteAt

    SUBROUTINE IntegerListDeleteAll(ilist)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist

        INTEGER :: allocstat

        ilist % entries = 0
        ilist % nelem = 0
    END SUBROUTINE IntegerListDeleteAll

    FUNCTION IntegerListFind(ilist, item) RESULT(ind)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: item
        INTEGER :: ind

        INTEGER :: i

        ind = -1
        DO i=1,ilist % nelem
            IF (ilist % entries(i) == item) THEN
                ind = i
                EXIT
            END IF
        END DO
    END FUNCTION IntegerListFind

    FUNCTION IntegerListGetArray(ilist) RESULT(arr)
        IMPLICIT NONE

        TYPE(IntegerList_t), TARGET :: ilist
        INTEGER, POINTER :: arr(:)

        arr => ilist % entries
    END FUNCTION IntegerListGetArray

    FUNCTION IntegerListGetSize(ilist) RESULT(N)
        IMPLICIT NONE
        TYPE(IntegerList_t) :: ilist
        INTEGER :: N

        N = ilist % nelem
    END FUNCTION IntegerListGetSize

! If iterators are needed, these need to be implemented

!!!     FUNCTION IntegerListGetIterator(ilist) RESULT(iliter)
!!!      IMPLICIT NONE
!!!      TYPE(IntegerList_t) :: ilist
!!!      TYPE(IntegerListIterator_t) :: iliter
!!!
!!!      ! NIY
!!!      iliter => NULL()
!!!    END FUNCTION getiterator
!!!
!!!    FUNCTION IntegerListIteratorNext(iliter) RESULT(nextelem)
!!!      IMPLICIT NONE
!!!      TYPE(IntegerListIterator_t) :: iliter
!!!      INTEGER :: nextelem
!!!
!!!      ! NIY
!!!      nextelem = 0
!!!    END FUNCTION IntegerListIteratorNext
!!!
!!!    FUNCTION IntegerListIteratorHasNext(iliter) RESULT(hasnext)
!!!      IMPLICIT NONE
!!!      TYPE(IntegerListIterator_t) :: iliter
!!!      LOGICAL :: hasnext
!!!
!!!      ! NIY
!!!      hasnext = .FALSE.
!!!    END FUNCTION IntegerListIteratorHasNext

    ! Given CRS indices, referenced indirectly from graph, 
    ! evenly load balance the work among the nthr threads
    SUBROUTINE ThreadLoadBalanceElementNeighbour(nthr, gn, gptr, gind, rptr, blkleads)
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

        ! TODO: Use Elmer messaging interface here
        ! WRITE (*,'(A,I0,A,I0)') 'ThreadLoadBalanceElementNeighbour: Total global work=', gwrk, &
        !                         ' work per element=', thrwrk

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
            ! TODO: Use Elmer messaging interface here
            ! WRITE (*,'(A,I0,A,I0)') 'ThreadLoadBalanceElementNeighbour: assigned ', wrk ,&
            !                         ' work to thread ', i
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
          ! totelem = totelem + blkleads(i+1)-blkleads(i)
        END DO
        ! Assign what is left of the matrix to the final thread
        blkleads(nthr+1)=gn+1
        
        ! totelem = totelem + blkleads(i+1)-blkleads(i)
        ! write (*,*) blkleads(1:nthr+1)
        ! write (*,*) gn, nthr
        ! IF (totelem /= gn) STOP
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

    ! Pad given integer value to be the next largest multiple of nbyte
    FUNCTION IntegerNBytePad(val, nbyte) RESULT(padval)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: val, nbyte
      INTEGER :: padval
      ! Parameters and variables
      INTEGER, PARAMETER :: bytesinint = KIND(val)
      INTEGER :: nbytesinint
      
      ! Compute number of nbytes in int
      nbytesinint = nbyte/bytesinint
      ! Compute value padded to multiples of n-byte
      padval=((val-1)/nbytesinint)*nbytesinint+nbytesinint
    END FUNCTION IntegerNBytePad

END MODULE LocalTypes
