SUBROUTINE MeshColour_init( Model,Solver,dt,TransientSimulation )
    USE DefUtils
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NormInd
    LOGICAL :: Found

    Params => GetSolverParams()

    NormInd = ListGetInteger( Params,'Norm Variable Index',Found)
    IF( NormInd > 0 ) THEN
      IF( .NOT. ListCheckPresent( Params,'Variable') ) THEN
        CALL ListAddString( Solver % Values,'Variable',&
                '-nooutput -global meshcolour_var')
      END IF
    END IF
    
END SUBROUTINE MeshColour_init

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
    TYPE(Element_t), POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh

    TYPE(Graph_t) :: DualGraph 
    TYPE(GraphColour_t) :: GraphColouring
    TYPE(Graph_t) :: ColourIndexList

#ifdef HAVE_METIS
    TYPE(VertexMap_t), TARGET :: DualGraphVM
#endif
    REAL(kind=dp) :: t_start, t_end

    INTEGER :: nerror, nerror_metis, nerror_colour, nerror_clist

    nerror = 0
    Mesh => GetMesh()

    ! Create dual mesh
    t_start = ftimer()
    CALL ElmerMeshToDualGraph(Mesh, DualGraph)
    t_end = ftimer()
    WRITE (*,'(A,ES12.3,A)') 'Dual graph creation total: ', t_end - t_start, ' sec.'

    ! Verify dual mesh with Metis (for completeness, 
    ! not compiled in by default)
#ifdef HAVE_METIS
    t_start = ftimer()
    CALL VertexMapFromArray(DualGraphVM, DualGraph % n, &
            DualGraph % ptr, DualGraph % ind)
    WRITE (*,'(A,ES12.3,A)') 'Vertex map creation from array total: ', t_end - t_start, ' sec.'
    CALL MeshToDualMetisVerify(Mesh, DualGraphVM, nerror_metis)
    nerror = nerror + nerror_metis
    CALL VertexMapDeleteAll(DualGraphVM)
#endif 

    ! Colour mesh
    t_start = ftimer()
    CALL ElmerGraphColour(DualGraph, GraphColouring)
    t_end = ftimer()
    WRITE (*,'(A,ES12.3,A)') 'Graph colouring total: ', t_end - t_start, ' sec.'
    WRITE (*,'(A,I0)') 'Number of colours created ngc=', GraphColouring % nc
    CALL GraphColourVerify(DualGraph % n, DualGraph % ptr, & 
            DualGraph % ind, GraphColouring % nc, GraphColouring % colours, &
            nerror_colour)
    nerror = nerror + nerror_colour
    CALL Graph_Deallocate(DualGraph)

    t_start = ftimer()
    CALL ElmerColouringToGraph(GraphColouring, ColourIndexList)
    t_end = ftimer()
    WRITE (*,'(A,ES12.3,A)') 'Colour gather total: ', t_end-t_start, ' sec.'
    CALL GraphColourListVerify(GraphColouring % nc, GraphColouring % colours, &
            ColourIndexList % ptr, ColourIndexList % ind, nerror_clist)
    nerror = nerror + nerror_clist
    CALL Colouring_Deallocate(GraphColouring)
    CALL Graph_Deallocate(ColourIndexList)

    ! Build solution norm for error checking
    Solver % Variable % Norm = REAL(1+nerror,dp)
    Solver % Variable % Values = REAL(1+nerror,dp)
    
CONTAINS

  ! Portable wall-clock timer
  FUNCTION ftimer() RESULT(timerval)
    IMPLICIT NONE
    
    REAL(KIND=dp) :: timerval
    INTEGER(KIND=8) :: t, rate
    
#ifdef _OPENMP
    timerval = OMP_GET_WTIME()
#else
    CALL SYSTEM_CLOCK(t,count_rate=rate)
    timerval = REAL(t,dp)/rate
#endif
  END FUNCTION ftimer

#ifdef HAVE_METIS
    SUBROUTINE MeshToDualMetisVerify(Mesh, DualGraph, err)
        IMPLICIT NONE
        TYPE(Mesh_t) :: Mesh
        TYPE(VertexMap_t), TARGET :: DualGraph
        INTEGER :: err

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
        REAL(kind=dp) :: t_start, t_end

        err = 0
        Elements => Mesh % Elements

        ! Set up parameters
        ne = Mesh % NumberOfBulkElements
        nn = Mesh % NumberOfNodes
        ncommon = 1 ! Dual graph will have an edge if elements share a node
        ! numflag = 1 ! Fortran-style numbering
        numflag = 0 ! C-style numbering

        t_start = ftimer()
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
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'Mesh transformation, Metis: ', t_end - t_start, ' sec.'

        t_start = ftimer()
        ! METIS call to construct dual graph
        retval = METIS_MeshToDual(ne, nn, eptr, eind, ncommon, numflag, r_xadj, r_adjncy)
        t_end = ftimer()
        WRITE (*,'(A,ES12.3,A)') 'Dual graph creation, Metis: ', t_end - t_start, ' sec.'

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

            IF (IntegerListGetSize(VertexMapGetList(DualGraph, i)) == nl) THEN
                DO j=dualptr(i),dualptr(i+1)-1
                    IF (.NOT. VertexMapFind(DualGraph, i, dualind(j))) THEN
                        WRITE (*,*) 'ERROR: Could not find node=', i, ', edge=', dualind(j)
                        graphOk = .FALSE.
                        err = err + 1
                    END IF
                END DO
            ELSE
                WRITE (*,*) 'ERROR: Size of local edge list does not match for node=', i
                graphOk = .FALSE.
                err = err + 1
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

    SUBROUTINE GraphColourVerify(gn, gptr, gind, ngc, gc, err)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: gn
        INTEGER, INTENT(IN) :: gptr(:), gind(:)
        INTEGER, INTENT(IN) :: ngc
        INTEGER, INTENT(IN) :: gc(:)
        INTEGER :: err

        INTEGER :: c, v, w, vli, vti, vcol, wind, wcol
        REAL(KIND=dp) :: avg, dev
        INTEGER :: ccount(ngc)
        LOGICAL :: colourOk

        err = 0
        ccount = 0
        colourOk = .TRUE.
        ! Verify and count colours (in serial!)
        DO v=1,gn
            vli = gptr(v)
            vti = gptr(v+1)-1
            ! Get colour of v
            vcol = gc(v)

            ! Verify that colour is in range
            IF (vcol<1 .OR. vcol>ngc) THEN
                WRITE (*,'(A,I0,A,I0,A)') 'ERROR: Graph vertex v=', v, &
                        ' colour ', vcol, ' out of range' 
                colourOk = .FALSE.
                err = err + 1
                CYCLE
            ELSE
                ccount(vcol)=ccount(vcol)+1
            END IF

            ! Check colour versus each neighbour
            DO wind=vli,vti
                w = gind(wind)
                wcol = gc(w)
                IF (wcol == vcol) THEN
                    WRITE (*,'(A,I0,A,I0,A,I0)') 'ERROR: Neighbouring vertices (v,w)=(', v,',', w, &
                            ') of the same colour col=', vcol
                    colourOk = .FALSE.
                    err = err + 1
                END IF
            END DO
        END DO

        ! Compute average
        avg = REAL(SUM(ccount),dp)/ngc

        WRITE (*,'(A,I0,/,A,I0,/,A,ES12.3)') 'Number of vertices, n=', gn, &
                'Number of coloured vertices, nc=', SUM(ccount), &
                'Average vertices per colour avg=', avg
        ! Print out statistics
        DO c=1,ngc
            dev = ABS(avg-ccount(c))
            WRITE (*,'(A,I0,A,I0,A,ES12.3)') 'Colour c=', c, ', count=', ccount(c), ', average dev=', dev
        END DO
        IF (colourOk) THEN
            WRITE (*,'(A)') 'Colouring seems ok.'
        ELSE
            WRITE (*,'(A)') 'ERROR: Colouring seems inconsistent!'
        END IF
    END SUBROUTINE GraphColourVerify

    SUBROUTINE GraphColourListVerify(ngc, gc, cptr, cind, err)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ngc
        INTEGER, INTENT(IN) :: gc(:)
        INTEGER, INTENT(IN) :: cptr(:), cind(:)
        INTEGER :: err

        INTEGER :: ccount(ngc)
        INTEGER, ALLOCATABLE :: cverify(:)
        INTEGER :: i, j, cli, cti,  n, ncol, totcol
        LOGICAL :: listsOk

        err = 0
        listsOk = .TRUE.

        n=size(gc)
        ! Count colours
        ccount = 0
        DO i=1,n
            ccount(gc(i))=ccount(gc(i))+1
        END DO

        ! Verify list pointers
        IF (SIZE(cptr) /= ngc+1 .OR. SIZE(cind) /= n) THEN
            WRITE (*,*) 'ERROR: Colour list pointer size does not', &
                    ' match the number of colours'
            err = err + 1
            RETURN
        END IF
        totcol = 0
        DO i=1,ngc
            ncol = cptr(i+1)-cptr(i)
            IF (ncol /= ccount(i)) THEN
                WRITE (*,'(3(A,I0))') 'ERROR: Colour=', i, ': pointer=', ncol,', count=', ccount(i)
                listsOk = .FALSE.
                err = err + 1
            END IF
        END DO
        ! Further verification of no use since pointers to lists are incorrect
        IF (.NOT. listsOk) RETURN

        IF (SUM(ccount) /= n) THEN
            WRITE (*,*) 'ERROR: Not enough colours in lists to cover the graph'
            listsOk = .FALSE.
        END IF

        ! Verify colours themselves
        ALLOCATE(cverify(n))

        cverify=0
        DO i=1,ngc
            cli = cptr(i)
            cti = cptr(i+1)-1
            DO j=cli,cti
                cverify(cind(j))=cverify(cind(j))+1
            END DO
        END DO
        DO i=1,n
            IF (cverify(i) > 1 .OR. cverify(i) < 1) THEN
                WRITE (*,'(2(A,I0))') 'ERROR: Vertex=', i, ', colour count=', cverify(i)
                listsOk = .FALSE.
                err = err + 1
            END IF
        END DO

        DEALLOCATE(cverify)

        IF (listsOk) THEN
            WRITE (*,'(A)') 'Colour lists seem ok.'
        ELSE
            WRITE (*,'(A)') 'ERROR: Colour lists seem inconsistent!'
        END IF
    END SUBROUTINE GraphColourListVerify

!------------------------------------------------------------------------------
END SUBROUTINE MeshColour
!------------------------------------------------------------------------------
