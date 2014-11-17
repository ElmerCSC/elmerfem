SUBROUTINE MeshColour( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Test multithreaded mesh colouring routines in Elmer
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
    USE DefUtils
    USE LocalTypes
    USE ISO_C_BINDING
!------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver

    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
    TYPE(Element_t),POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh

    TYPE(VertexMap_t), TARGET :: DualGraph
    INTEGER :: n, ngc
    INTEGER, ALLOCATABLE :: dualptr(:), dualind(:), colours(:)
#ifdef HAVE_TIMING
    REAL(kind=dp) :: t_start, t_end
#endif    

    Mesh => GetMesh()

    ! Create dual mesh
#ifdef HAVE_TIMING
    t_start = ftimer()
#endif 
    CALL ElmerMeshToDualGraph(Mesh, n, dualptr, dualind)
#ifdef HAVE_TIMING
    t_end = ftimer()
    WRITE (*,'(A,ES12.3,A)') 'Dual graph creation total: ', t_end - t_start, ' sec.'
#endif    

    ! Verify dual mesh
#ifdef HAVE_METIS
#ifdef HAVE_TIMING
    t_start = ftimer()
#endif 
    CALL VertexMapFromArray(DualGraph, n, dualptr, dualind)
#ifdef HAVE_TIMING
    t_end = ftimer()
    WRITE (*,'(A,ES12.3,A)') 'Vertex map creation from array total: ', t_end - t_start, ' sec.'
#endif    
    CALL MeshToDualMetisVerify(Mesh, DualGraph)
    CALL VertexMapDeleteAll(DualGraph)
#endif 

    ! Colour mesh
#ifdef HAVE_TIMING
    t_start = ftimer()
#endif 
    CALL ElmerGraphColour(n, dualptr, dualind, ngc, colours)
#ifdef HAVE_TIMING
    t_end = ftimer()
    WRITE (*,'(A,ES12.3,A)') 'Graph colouring total: ', t_end - t_start, ' sec.'
#endif    
    WRITE (*,'(A,I0)') 'Number of colours created ngc=', ngc
    CALL GraphColourVerify(n, dualptr, dualind, ngc, colours)

    DEALLOCATE(dualptr, dualind, colours)
    CONTAINS

#ifdef HAVE_METIS
        SUBROUTINE MeshToDualMetisVerify(Mesh, DualGraph)
            IMPLICIT NONE
            TYPE(Mesh_t) :: Mesh
            TYPE(VertexMap_t), TARGET :: DualGraph
            
            !------------------------------------------------------------------------------
            ! Metis MeshToDual interface
            !------------------------------------------------------------------------------
            INTERFACE
                ! Assume 32-bit version of Metis, i.e., idx_t = int32_t
                ! int METIS_MeshToDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, 
                !                      idx_t *ncommon, idx_t *numflag,  idx_t **r_xadj, idx_t **r_adjncy)
                FUNCTION METIS_MeshToDual(ne, nn, eptr, eind, ncommon, numflag,  r_xadj, r_adjncy) RESULT(val) & 
                      BIND(C,name='METIS_MeshToDual') 
                    USE ISO_C_BINDING
                    IMPLICIT NONE
                    INTEGER(KIND=c_int32_t) :: ne, nn           ! idx_t *ne, idx_t *nn
                    INTEGER(KIND=c_int32_t) :: eptr(*), eind(*) ! idx_t *eptr, idx_t *eind 
                    INTEGER(KIND=c_int32_t) :: ncommon, numflag ! idx_t *ncommon, idx_t *numflag
                    TYPE(C_PTR) :: r_xadj, r_adjncy             ! idx_t **r_xadj, idx_t **r_adjncy
                    ! INTEGER(KIND=c_int32_t) :: r_xadj(*), r_adjncy(*)       ! idx_t **r_xadj, idx_t **r_adjncy
                    INTEGER(KIND=C_INT) :: val
                END FUNCTION METIS_MeshToDual
                
                ! METIS_API(int) METIS_Free(void *ptr);
                FUNCTION METIS_Free(ptr) RESULT(val) BIND(C,name='METIS_Free') 
                    USE ISO_C_BINDING
                    IMPLICIT NONE
                    TYPE(C_PTR), VALUE :: ptr                          ! void *ptr
                    INTEGER(KIND=C_INT) :: val
                END FUNCTION METIS_Free
            END INTERFACE
            
            INTEGER :: i, j, allocstat, nl, nli, nti
            INTEGER(KIND=C_INT) :: retval
            INTEGER(KIND=c_int32_t) :: ne, nn, ncommon, numflag
            INTEGER(KIND=c_int32_t), POINTER CONTIG :: eptr(:), eind(:), dualptr(:), dualind(:)
            TYPE(C_PTR) :: r_xadj, r_adjncy
            ! INTEGER(KIND=c_int32_t), POINTER CONTIG :: r_xadj(:), r_adjncy(:)
            ! TYPE(IntegerList_t), POINTER :: elist
            LOGICAL :: graphOk
            TYPE(Element_t), POINTER :: Element, Elements(:)
            TYPE(IntegerList_t), POINTER :: vlist
            INTEGER, POINTER :: varr(:)
            INTEGER :: vsize
#ifdef HAVE_TIMING
            REAL(kind=dp) :: t_start, t_end
#endif    
 
            Elements => Mesh % Elements

            ! Set up parameters
            ne = Mesh % NumberOfBulkElements
            nn = Mesh % NumberOfNodes
            ncommon = 1 ! Dual graph will have an edge if elements share a node
            ! numflag = 1 ! Fortran-style numbering
            numflag = 0 ! C-style numbering

#ifdef HAVE_TIMING
            t_start = ftimer()
#endif      
            ! Copy mesh to CSR structure
            ALLOCATE(eptr(ne+1), eind(ne*Mesh % MaxElementNodes))
            eptr(1)=1 ! Fortran numbering
            DO i=1, ne
                Element => Elements(i)
                nl = Element % TYPE % NumberOfNodes
                nli = eptr(i) ! Fortran numbering
                nti = nli+nl-1
                eind(nli:nti) = Element % NodeIndexes(1:nl) ! Fortran numbering
                eptr(i+1) = nli+nl
            END DO

            ! Convert eptr and eind to C numbering
            DO i=1,eptr(ne+1)-1
                eind(i)=eind(i)-1
            END DO
            DO i=1,ne+1
                eptr(i)=eptr(i)-1
            END DO
#ifdef HAVE_TIMING
            t_end = ftimer()
            WRITE (*,'(A,ES12.3,A)') 'Mesh transformation, Metis: ', t_end - t_start, ' sec.'
#endif                  

            ! WRITE (*,*) eptr(1:3)
            ! WRITE (*,*) eind(eptr(1):eptr(4)-1)
            ! STOP

#ifdef HAVE_TIMING
            t_start = ftimer()
#endif      
            ! METIS call to construct dual graph
            retval = METIS_MeshToDual(ne, nn, eptr, eind, ncommon, numflag, r_xadj, r_adjncy)
#ifdef HAVE_TIMING
            t_end = ftimer()
            WRITE (*,'(A,ES12.3,A)') 'Dual graph creation, Metis: ', t_end - t_start, ' sec.'
#endif                  

            ! TODO: Verify the graph against the given dual
            CALL C_F_POINTER(r_xadj, dualptr, (/ ne+1 /))
            ! Convert dualptr to Fortran numbering
            DO i=1,ne+1
                dualptr(i)=dualptr(i)+1
            END DO

            CALL C_F_POINTER(r_adjncy, dualind, (/ dualptr(ne+1)-1 /))      
            ! Convert dualind to Fortran numbering
            DO i=1,dualptr(ne+1)-1
                dualind(i)=dualind(i)+1
            END DO

            

            ! Verify graph structure
            graphOk = .TRUE.
            DO i=1,ne
                nl = dualptr(i+1)-dualptr(i)
           
                vlist => VertexMapGetList(DualGraph,i)
                vsize = IntegerListGetSize(vlist)
                varr => IntegerListGetArray(vlist)
                ! IF (i==2) THEN
                !         WRITE (*,*) 'Metis:'
                !         WRITE (*,*) dualind(dualptr(i):dualptr(i+1)-1)
                !         WRITE (*,*) 'Elmer:'
                !         WRITE (*,*) varr(1:vsize)
                !         STOP    
                ! END IF

                IF (IntegerListGetSize(VertexMapGetList(DualGraph, i)) == nl) THEN
                    DO j=dualptr(i),dualptr(i+1)-1
                        IF (.NOT. VertexMapFind(DualGraph, i, dualind(j))) THEN
                            WRITE (*,*) 'ERROR: Could not find node=', i, ', edge=', dualind(j)
                            graphOk = .FALSE.
                        END IF
                    END DO
                ELSE
                    WRITE (*,*) 'ERROR: Size of local edge list does not match for node=', i
                    graphOk = .FALSE.
                END IF
            END DO
            IF (graphOk) THEN
                WRITE (*,'(A)') 'Dual graph matches Metis output'
            ELSE
                WRITE (*,*) 'ERROR: Dual graph noes not match Metis output'
            END IF

            ! Free memory allocated by Metis
            retval = Metis_Free(r_xadj)  
            retval = Metis_Free(r_adjncy)

            ! Free memory reserved for mesh
            DEALLOCATE(eptr, eind)

        END SUBROUTINE MeshToDualMetisVerify
#endif

        SUBROUTINE MeshToDualGraph(Mesh, DualGraph)
            IMPLICIT NONE

            TYPE(Mesh_t) :: Mesh
            TYPE(VertexMap_t) :: DualGraph

            TYPE(Element_t), POINTER :: Element, NElement, Elements(:)
            TYPE(VertexMap_t), TARGET :: VToEMap
            TYPE(IntegerList_t), POINTER :: NeighbourList
            INTEGER, POINTER :: NeighbourArray(:)
            INTEGER :: i, j, e, ne, eid, neid, nodeid, nelem, nnelem, nvertex
            INTEGER, POINTER CONTIG :: NodeIndexes(:)
#ifdef HAVE_TIMING
            REAL(kind=dp) :: t_start, t_end
#endif    

            nelem = Mesh % NumberOfBulkElements
            nvertex = Mesh % NumberOfNodes
            Elements => Mesh % Elements

#ifdef HAVE_TIMING
    t_start = ftimer()
#endif 
            CALL ConstructVertexToElementMap(Mesh, VToEMap)
#ifdef HAVE_TIMING
            t_end = ftimer()
            WRITE (*,'(A,ES12.3,A)') 'Vertex-element map creation: ', t_end - t_start, ' sec.'
#endif    

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
            !$OMP PARALLEL DO SHARED(Elements, nvertex, nelem, VToEMap, DualGraph) &
            !$OMP PRIVATE(e, eid, ne, neid, nnelem, Element, NElement, NeighbourList, &
            !$OMP         NeighbourArray, NodeIndexes) DEFAULT(NONE)
            DO e=1,nelem
                Element => Elements(e)
                NodeIndexes => Element % NodeIndexes
                eid = Element % ElementIndex

                ! For each node nodeid in element eid
                DO nodeid=1, Element % TYPE % NumberOfNodes
                    ! Get list of elements mapping to node
                    NeighbourList => VertexMapGetList(VToEMap, NodeIndexes(nodeid))
                    NeighbourArray => IntegerListGetArray(NeighbourList)
                    nnelem = IntegerListGetSize(NeighbourList)
                    
                    ! For each neighbour element ne of nodeid
                    DO ne=1, nnelem
                        neid = NeighbourArray(ne)
                                                
                        ! Add each actual neighbouring element to graph 
                        IF (eid /= neid) CALL VertexMapAdd(DualGraph, eid, neid)
                    END DO
                END DO
            END DO
            !$OMP END PARALLEL DO

#ifdef HAVE_TIMING
            t_end = ftimer()
            WRITE (*,'(A,ES12.3,A)') 'Dual graph creation: ', t_end - t_start, ' sec.'
#endif 

            ! Delete entries in vertex to element map
            CALL VertexMapDeleteAll(VToEMap)
        END SUBROUTINE MeshToDualGraph

        SUBROUTINE MeshToDualGraph2(Mesh, DualGraph)
            IMPLICIT NONE

            TYPE(Mesh_t) :: Mesh
            TYPE(VertexMap_t) :: DualGraph

            TYPE(Element_t), POINTER :: Element, NElement
            TYPE(VertexMap_t), TARGET :: VToEMap
            TYPE(IntegerList_t), POINTER :: NeighbourList
            INTEGER, POINTER :: NeighbourArray(:)
            INTEGER :: e, ne, eid, neid, nvertex, belem, nelem, vertexId

#ifdef HAVE_TIMING
            REAL(kind=dp) :: t_start, t_end
#endif    

            belem = Mesh % NumberOfBulkElements
            nvertex = Mesh % NumberOfNodes

#ifdef HAVE_TIMING
    t_start = ftimer()
#endif 
            CALL ConstructVertexToElementMap(Mesh, VToEMap)
#ifdef HAVE_TIMING
            t_end = ftimer()
            WRITE (*,'(A,ES12.3,A)') 'Vertex-element map creation: ', t_end - t_start, ' sec.'
#endif    

#ifdef HAVE_TIMING
            t_start = ftimer()
#endif 

            ! An alternative approach: loop over vertices and add edges 
            ! between elements in vertex lists. This is prone to have a worse
            ! load balance that the other approach
            CALL VertexMapInit(DualGraph, belem)

            ! For each element
            !$OMP PARALLEL DO SHARED(nvertex, VToEMap, DualGraph) &
            !$OMP PRIVATE(e, ne, eid, neid, nelem, vertexId, NeighbourList, &
            !$OMP         NeighbourArray) &
            !$OMP DEFAULT(NONE)
            DO vertexId=1, nvertex
                ! Get list of elements mapping to a node
                NeighbourList => VertexMapGetList(VToEMap, vertexId)
                NeighbourArray => IntegerListGetArray(NeighbourList)
                nelem = IntegerListGetSize(NeighbourList)

                ! For each element in list
                DO e=1, nelem
                    eid = NeighbourArray(e)
                    ! For each neighbour element ne of eid
                    DO ne=1, e-1
                        ! Add edge from e to ne
                        CALL VertexMapAddAtomic(DualGraph, &
                                          eid, &
                                          NeighbourArray(ne))
                    END DO
                    DO ne=e+1,nelem
                        ! Add edge from e to ne
                        CALL VertexMapAddAtomic(DualGraph, &
                                          eid, &
                                          NeighbourArray(ne))
                    END DO
                END DO
            END DO
            !$OMP END PARALLEL DO

#ifdef HAVE_TIMING
            t_end = ftimer()
            WRITE (*,'(A,ES12.3,A)') 'Dual graph creation: ', t_end - t_start, ' sec.'
#endif 
            ! Delete entries in vertex to element map
            CALL VertexMapDeleteAll(VToEMap)

        END SUBROUTINE MeshToDualGraph2


        SUBROUTINE ConstructVertexToElementMap(Mesh, VertexToElementMap)
            IMPLICIT NONE

            TYPE(Mesh_t) :: Mesh
            TYPE(VertexMap_t) :: VertexToElementMap
            
            TYPE(Element_t), POINTER :: Element, Elements(:)
            INTEGER :: i, j, v, nelem, nvertex, allocstat
            
            nelem = Mesh % NumberOfBulkElements
            nvertex = Mesh % NumberOfNodes
            Elements => Mesh % Elements
            
            ! Initialize map
            CALL VertexMapInit(VertexToElementMap, nvertex)
            
            ! For each element
            !$OMP PARALLEL DO SHARED(Elements, nvertex, nelem, VertexToElementMap) &
            !$OMP PRIVATE(i, j, v, Element) DEFAULT(NONE)
            DO i=1,nelem
                Element => Elements(i)
                ! For each node
                DO j=1, Element % TYPE % NumberOfNodes
                    ! Set node i to map to v
                    CALL VertexMapAddAtomic(VertexToElementMap, &
                          Element % NodeIndexes(j), &
                          Element % ElementIndex)
                END DO
            END DO
            !$OMP END PARALLEL DO
        END SUBROUTINE ConstructVertexToElementMap        

!------------------------------------------------------------------------------
END SUBROUTINE MeshColour
!------------------------------------------------------------------------------
