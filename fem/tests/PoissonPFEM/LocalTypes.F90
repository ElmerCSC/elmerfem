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

CONTAINS
    
    SUBROUTINE MeshToDualGraph3(Mesh, DualGraph)
        IMPLICIT NONE

        TYPE(Mesh_t) :: Mesh
        TYPE(VertexMap_t) :: DualGraph

        TYPE(Element_t), POINTER :: Element, NElement, Elements(:)
        TYPE(IntegerList_t), ALLOCATABLE, TARGET :: VToEList(:)
        TYPE(IntegerList_t), POINTER :: NeighbourList
        INTEGER, POINTER :: NeighbourArray(:)
        INTEGER :: i, j, ind, nli, nti, nl, e, ne, eid, neid, nodeid, nelem, nnelem, nvertex
        INTEGER, POINTER CONTIG :: NodeIndexes(:)
        INTEGER, ALLOCATABLE :: eptr(:), eind(:), vptr(:), vind(:)
        LOGICAL :: mapOk
#ifdef HAVE_TIMING
        REAL(kind=dp) :: t_start, t_end
#endif    

        nelem = Mesh % NumberOfBulkElements
        nvertex = Mesh % NumberOfNodes
        Elements => Mesh % Elements

#ifdef HAVE_TIMING
        t_start = ftimer()
#endif      
        ! Copy mesh to CSR structure
        ALLOCATE(eptr(nelem+1), eind(nelem*Mesh % MaxElementNodes))
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
        WRITE (*,'(A,ES12.3,A)') 'Mesh tranformation, meshtodual: ', t_end - t_start, ' sec.'
#endif                              

#ifdef HAVE_TIMING
        t_start = ftimer()
#endif 
        CALL ConstructVertexToElementList(nelem, nvertex, eptr, eind, VToEList)
#ifdef HAVE_TIMING
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'Vertex-element list creation: ', t_end - t_start, ' sec.'
#endif    
#ifdef HAVE_TIMING
        t_start = ftimer()
#endif 
        CALL ConstructVertexToElementList2(nelem, nvertex, eptr, eind, vptr, vind)
#ifdef HAVE_TIMING
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'Vertex-element list creation2: ', t_end - t_start, ' sec.'
#endif    
        ! Compare list to the correct one
        mapOk = .TRUE.
        DO i=1,nvertex
            nl = vptr(i+1)-vptr(i)
            IF (nl == IntegerListGetSize(VToEList(i))) THEN
                DO j=vptr(i), vptr(i+1)-1
                    IF (IntegerListFind(VToEList(i), vind(j))<0) THEN
                        WRITE (*,*) 'ERROR: Element not found vertex=', i, vind(j)
                    END IF
                END DO
            ELSE
                WRITE (*,*) 'ERROR: Map is of different size for vertex=', i
            END IF
        END DO
        IF (mapOk) THEN
            WRITE (*,*) 'Vertex-element list correct'
        ELSE
            WRITE (*,*) 'ERROR: Vertex-element list incorrect'
            STOP
        END IF

#ifdef HAVE_TIMING
        t_start = ftimer()
#endif 
        ! Algorithm: loop over elements and add edges 
        ! between elements for each vertex. 
        CALL VertexMapInit(DualGraph, nelem)
#ifdef HAVE_TIMING
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'Dual graph init: ', t_end - t_start, ' sec.'
#endif 

#ifdef HAVE_TIMING
        t_start = ftimer()
#endif 
        ! For each element
        !$OMP PARALLEL DO SHARED(eptr, eind, nvertex, nelem, vptr, vind, DualGraph) &
        !$OMP PRIVATE(eid, ind, nli, nti, ne, neid, nnelem, NeighbourList, &
        !$OMP         NeighbourArray, NodeIndexes) SCHEDULE(GUIDED) DEFAULT(NONE)
        DO eid=1,nelem
            nli = eptr(eid)
            nti = eptr(eid+1)-1

            ! For each node nodeid in element eid
            DO nodeid=nli, nti
                ! Get list of elements mapping to a vertex
                ind = eind(nodeid)
                ! NeighbourList => VToEList(eind(nodeid))
                ! NeighbourArray => IntegerListGetArray(NeighbourList)
                ! nnelem = IntegerListGetSize(NeighbourList)

                ! For each neighbour element ne of nodeid
                DO ne=vptr(ind), vptr(ind+1)-1
                    ! neid = NeighbourArray(ne)
                    neid = vind(ne)
                    
                    ! Add each actual neighbouring element to graph 
                    IF (eid /= neid) CALL VertexMapAdd(DualGraph, eid, neid)
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO

        DEALLOCATE(eind, eptr, vptr, vind)

#ifdef HAVE_TIMING
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'Dual graph creation: ', t_end - t_start, ' sec.'
#endif 

        ! Delete entries in vertex to element map
        DO i=1,nvertex
            CALL IntegerListDeleteAll(VToEList(i))
        END DO
        DEALLOCATE(VToEList)
    END SUBROUTINE MeshToDualGraph3

    SUBROUTINE MeshToDualGraph4(Mesh, n, dualptr, dualind)
        IMPLICIT NONE

        TYPE(Mesh_t) :: Mesh
        INTEGER, INTENT(OUT) :: n
        INTEGER, ALLOCATABLE :: dualptr(:), dualind(:)

        TYPE(Element_t), POINTER :: Element, NElement, Elements(:)
        TYPE(IntegerList_t), ALLOCATABLE, TARGET :: VToEList(:)
        TYPE(IntegerList_t), POINTER :: NeighbourList
        INTEGER, POINTER :: NeighbourArray(:)
        INTEGER :: i, j, ind, nli, nti, nl, nc, nv, vid, vli, vti, e, ne, eid, neid, &
              nodeid, nelem, nnelem, nvertex, allocstat
        INTEGER, POINTER CONTIG :: NodeIndexes(:)
        INTEGER, POINTER :: tmparr(:)
        INTEGER, ALLOCATABLE :: eptr(:), eind(:), vptr(:), vind(:), &
              wrkliptr(:), wrktiptr(:), wrkind(:), wrkindnew(:), wrkheap(:)
        LOGICAL :: mapOk
        TYPE(IntegerList_t) :: dualindlist
#ifdef HAVE_TIMING
        REAL(kind=dp) :: t_start, t_end
#endif    

        nelem = Mesh % NumberOfBulkElements
        n = nelem
        nvertex = Mesh % NumberOfNodes
        Elements => Mesh % Elements

#ifdef HAVE_TIMING
        t_start = ftimer()
#endif      
        ! Copy mesh to CSR structure
        ALLOCATE(eptr(nelem+1), eind(nelem*Mesh % MaxElementNodes))
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
        WRITE (*,'(A,ES12.3,A)') 'Mesh tranformation, meshtodual: ', t_end - t_start, ' sec.'
#endif                              

! #ifdef HAVE_TIMING
!         t_start = ftimer()
! #endif 
!         CALL ConstructVertexToElementList(nelem, nvertex, eptr, eind, VToEList)
! #ifdef HAVE_TIMING
!         t_end = ftimer()
!         WRITE (*,'(A,ES12.3,A)') 'Vertex-element list creation: ', t_end - t_start, ' sec.'
! #endif    
#ifdef HAVE_TIMING
        t_start = ftimer()
#endif 
        CALL ConstructVertexToElementList2(nelem, nvertex, eptr, eind, vptr, vind)
#ifdef HAVE_TIMING
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'Vertex-element list creation2: ', t_end - t_start, ' sec.'
#endif    
!         ! Compare list to the correct one
!         mapOk = .TRUE.
!         DO i=1,nvertex
!             nl = vptr(i+1)-vptr(i)
!             IF (nl == IntegerListGetSize(VToEList(i))) THEN
!                 DO j=vptr(i), vptr(i+1)-1
!                     IF (IntegerListFind(VToEList(i), vind(j))<0) THEN
!                         WRITE (*,*) 'ERROR: Element not found vertex=', i, vind(j)
!                     END IF
!                 END DO
!             ELSE
!                 WRITE (*,*) 'ERROR: Map is of different size for vertex=', i
!             END IF
!         END DO
!         IF (mapOk) THEN
!             WRITE (*,*) 'Vertex-element list correct'
!         ELSE
!             WRITE (*,*) 'ERROR: Vertex-element list incorrect'
!             STOP
!         END IF
        
        ! Delete entries in vertex to element map
!         DO i=1,nvertex
!             CALL IntegerListDeleteAll(VToEList(i))
!         END DO
!         DEALLOCATE(VToEList)

#ifdef HAVE_TIMING
        t_start = ftimer()
#endif 

        ! Ensure that vertex to element lists are sorted
        DO vid=1,nvertex
            vli = vptr(vid)
            vti = vptr(vid+1)-1

            CALL Sort(vti-vli+1, vind(vli:vti))
        END DO

        ! Allocate work array (local to each thread)
        ! TODO: Get value for wrkind from somewhere
        ALLOCATE(wrkind(32), wrkliptr(Mesh % MaxElementNodes), &
                 wrktiptr(Mesh % MaxElementNodes), &
                 wrkheap(2*Mesh % MaxElementNodes), dualptr(nelem+1), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('MeshDualGraph', &
              'Unable to allocate local workspace!')

        ! TODO: Get the initial size for dualindlist from somewhere
        CALL IntegerListInit(dualindlist, nelem*20)

        ! For each element
        dualptr(1)=1
        DO eid=1,nelem
            nli = eptr(eid)
            nti = eptr(eid+1)-1
            nv = nti-nli+1
            
            ! Get pointers to vertices related to the nodes of the element
            DO i=nli,nti
                wrkliptr(i-nli+1)=vptr(eind(i))
                wrktiptr(i-nli+1)=vptr(eind(i)+1) ! NOTE: This is to make comparison cheaper
            END DO

            IF (eid == 2) THEN
                DO i=1, nv
                    WRITE (*,*) wrkliptr(i), wrktiptr(i)
                END DO

                DO i=1, nv
                    WRITE (*,*) vind(wrkliptr(i):wrktiptr(i)-1)
                END DO
             END IF

            ! Merge vertex lists (multi-way merge of ordered lists)
            CALL VertexListMerge(eid, nv, wrkliptr, wrktiptr, vind, nc, wrkind, wrkheap)
            
            ! Add merged list to final list of vertices
            CALL IntegerListAddArray(dualindlist, wrkind(1:nc))
            dualptr(eid+1)=dualptr(eid)+nc
        END DO

        ! TODO: In a multithreaded setting, merge the lists here
        ALLOCATE(dualind(dualptr(nelem+1)-1), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('MeshDualGraph', &
              'Unable to allocate dual mesh!')
        tmparr => IntegerListGetArray(dualindlist)
        dualind(1:dualptr(nelem+1)-1)=tmparr(1:dualptr(nelem+1)-1)

        CALL IntegerListDeleteAll(dualindlist)
        DEALLOCATE(eind, eptr, vptr, vind, wrkind, wrkliptr, wrktiptr, wrkheap)

#ifdef HAVE_TIMING
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'Dual graph creation: ', t_end - t_start, ' sec.'
#endif 
        CONTAINS 

            SUBROUTINE VertexListMerge(cnode, nv, wrkliptr, wrktiptr, vind, nc, wrkind, wrkheap)
                IMPLICIT NONE

                INTEGER, INTENT(IN) :: cnode, nv
                INTEGER :: wrkliptr(:)
                INTEGER, INTENT(IN) :: wrktiptr(:), vind(:)
                INTEGER, INTENT(OUT) :: nc
                INTEGER, ALLOCATABLE :: wrkind(:)
                INTEGER :: wrkheap(:)
                
                ! Local variables
                INTEGER :: i, j, wrki, vindi, nheap

                ! Add first vertex of each list to heap
                nheap = 0
                nc = 0
                DO i=1,nv
                    wrki = wrkliptr(i)
                    wrkliptr(i)=wrkliptr(i)+1
                    ! HEAP_INSERT(wrkheap, nheap, vind(wrki), i)
                    CALL IntegerHeapInsertTuple(nheap, wrkheap, vind(wrki), i)
                END DO
                
                DO WHILE (nheap > 0)
                    ! HEAP_EXTRACT(wrkheap, nheap, vindi, i)
                    CALL IntegerHeapExtractTuple(nheap, wrkheap, vindi, i)
                    IF (cnode == 2) WRITE (*,*) 'Got from', i, vindi, nheap

                    ! Do not add node if it is the current node
                    IF (vindi /= cnode) THEN
                        ! Do not add node if it was previously added to the list
                        IF (nc>0) THEN
                            IF (wrkind(nc) /= vindi) THEN
                                nc=nc+1
                                IF (nc>SIZE(wrkind)) THEN
                                    ALLOCATE(wrkindnew(SIZE(wrkind)*2))
                                    wrkindnew(1:nc-1)=wrkind(1:nc-1)
                                    DEALLOCATE(wrkind)
                                    CALL MOVE_ALLOC(wrkindnew, wrkind)
                                END IF
                                wrkind(nc)=vindi
                            END IF
                        ELSE
                            nc=1
                            wrkind(1)=vindi
                        END IF
                    END IF
                    IF (wrkliptr(i)<wrktiptr(i)) THEN
                        wrki = wrkliptr(i)
                        wrkliptr(i)=wrkliptr(i)+1
                        IF (cnode == 2) WRITE (*,*) 'Add from', i, wrki, vind(wrki), nheap
                        ! HEAP_INSERT(wrkheap, nheap, vind(wrki), i)
                        CALL IntegerHeapInsertTuple(nheap, wrkheap, vind(wrki), i)
                    END IF
                END DO
            END SUBROUTINE VertexListMerge

            SUBROUTINE IntegerHeapInsertTuple(n, heap, t1, t2)
                IMPLICIT NONE
                
                ! Parameters
                INTEGER :: n, heap(:)
                INTEGER, INTENT(IN) :: t1, t2

                ! Variables
                INTEGER :: j, k

                IF (eid == 2) THEN
                    WRITE (*,*) 'before heap n=', n
                    WRITE (*,*) heap(1:2*n)
                    WRITE (*,*) 'tuple=', t1, t2
                END IF
                ! Insert to ordered list
                j=1
                DO j=1,n
                    IF (t1<heap(2*j-1)) THEN
                        ! TODO: THIS IS THE WRONG WAY! 
                        DO k=j,n
                            heap(2*k+1)=heap(2*k-1)
                            heap(2*k+2)=heap(2*k)
                        END DO
                        EXIT
                    END IF
                END DO
                
                ! Position found and storage arranged, just add element to heap
                heap(2*j-1)=t1
                heap(2*j)=t2
                n=n+1

                IF (eid == 2) THEN
                    WRITE (*,*) 'after heap n=', n
                    WRITE (*,*) heap(1:2*n)
                    IF (n>7) STOP
                END IF
            END SUBROUTINE IntegerHeapInsertTuple

            SUBROUTINE IntegerHeapExtractTuple(n, heap, t1, t2)
                IMPLICIT NONE
                
                ! Parameters
                INTEGER :: n, heap(:)
                INTEGER, INTENT(OUT) :: t1, t2

                INTEGER :: j

                ! Extract element
                t1=heap(1)
                t2=heap(2)
                
                ! Move rest of the data to correct position
                DO j=1,n-1
                    ! TODO: CHECK THIS
                    heap(2*j-1)=heap(2*j+1)
                    heap(2*j)=heap(2*j+2)
                END DO
                n=n-1
            END SUBROUTINE IntegerHeapExtractTuple
    END SUBROUTINE MeshToDualGraph4
    
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
        INTEGER, ALLOCATABLE :: vlock(:)
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
        
        INTEGER :: i, j, vli, vti

        CALL VertexMapInit(vmap, n)
        
        ! For each vertex
        DO i=1,n
            vli=vptr(i)
            vti=vptr(i+1)-1

            ! Add mapping to array
            DO j=vli,vti
                CALL VertexMapAdd(vmap, i, vind(j))
            END DO
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
        IF (ind < 0 .OR. ind > SIZE(ilist % entries)) RETURN
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

        IF (ind < 0 .OR. ind > SIZE(ilist % entries)) RETURN

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

END MODULE LocalTypes
