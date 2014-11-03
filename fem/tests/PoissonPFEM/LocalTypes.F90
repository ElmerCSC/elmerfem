MODULE LocalTypes
#ifdef _OPENMP
    USE omp_lib
#endif
    USE Types
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

    TYPE VertexElementMap_t
#ifdef _OPENMP
        INTEGER(KIND=OMP_LOCK_KIND), ALLOCATABLE :: vlock(:)
#endif
        TYPE(IntegerHashSet_t), ALLOCATABLE :: map(:)
    END TYPE VertexElementMap_t

    INTEGER, PARAMETER :: INTEGERLIST_DEFAULT_SIZE = 64
    INTEGER, PARAMETER :: INTEGERHASHSET_DEFAULT_SIZE = 64
    INTEGER, PARAMETER :: INTEGERHASHSET_CHAIN_DEFAULT_SIZE = 8
    REAL(KIND=dp), PARAMETER :: INTEGERHASHSET_FILLRATIO = REAL(0.75, dp)

CONTAINS

    SUBROUTINE VertexElementMapInit(vmap, nvertex, mdim)
        IMPLICIT NONE

        TYPE(VertexElementMap_t) :: vmap
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
        IF (allocstat /= 0) CALL Fatal('VertexElementMapInit', &
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
        IF (allocstat /= 0) CALL Fatal('VertexElementMapInit', &
              'Memory allocation failed!')
        !$OMP END SINGLE

        !$OMP DO
        DO i=1,nvertex
            CALL OMP_INIT_LOCK(vmap % vlock(i))
        END DO
        !$OMP END DO NOWAIT
#endif

        !$OMP END PARALLEL
    END SUBROUTINE VertexElementMapInit

    SUBROUTINE VertexElementMapAdd(vmap, vertexId, elementId)
        IMPLICIT NONE

        TYPE(VertexElementMap_t) :: vmap
        INTEGER, INTENT(IN) :: vertexId, elementId

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
    END SUBROUTINE VertexElementMapAdd

    SUBROUTINE VertexElementMapDelete(vmap, vertexId, elementId)
        IMPLICIT NONE
        
        TYPE(VertexElementMap_t) :: vmap
        INTEGER, INTENT(IN) :: vertexId, elementId
      
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
    END SUBROUTINE VertexElementMapDelete

    SUBROUTINE VertexElementMapToList(vmap, nlist)
        IMPLICIT NONE

        
        TYPE(VertexElementMap_t) :: vmap
        TYPE(IntegerList_t), ALLOCATABLE :: nlist(:)

        INTEGER :: i, nvertex, allocstat
        TYPE(IntegerList_t) :: entries

        nvertex = SIZE(vmap % map)
        
        ALLOCATE(nlist(nvertex), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('VertexElementMapToList', &
              'Memory allocation failed!')
        
        !$OMP PARALLEL DO PRIVATE(i)
        DO i=1,nvertex
            ! Add all elements from map to a simple list
            CALL IntegerListInit(nlist(i), vmap % map(i) % entries % nelem)
            CALL IntegerListAddAll(nlist(i), vmap % map(i) % entries)
        END DO
        !$OMP END PARALLEL DO
    END SUBROUTINE VertexElementMapToList

    SUBROUTINE VertexElementMapDeleteAll(vmap)
        IMPLICIT NONE
        
        TYPE(VertexElementMap_t) :: vmap
        
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
    END SUBROUTINE VertexElementMapDeleteAll

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

        ALLOCATE(iset % set(FLOOR(n/iset % fratio)), &
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

END MODULE LocalTypes
