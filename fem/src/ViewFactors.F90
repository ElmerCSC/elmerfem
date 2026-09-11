!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! * ELMER/FEM Viewfactor computation
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

!> \ingroup Programs
!> \{

!> \defgroup ViewFactors Program ViewFactors
!> \{

   MODULE ViewFactorGlobals
     USE Types
     REAL(KIND=dp), ALLOCATABLE, TARGET :: Jdiag(:), Jacobian(:,:)
   END MODULE ViewFactorGlobals

!------------------------------------------------------------------------------
!> A separate program that computes the view factors to an external file. 
!> This file is later used within the ElmerSolver. If the view factor files
!> do not exist, a system call for this program is performed. 
!------------------------------------------------------------------------------
      
   PROGRAM ViewFactors
   
     USE ViewUtils
     USE DefUtils
     USE ViewFactorGlobals
     USE MeshTransform, ONLY : RigidMeshMapping
     USE MainUtils, ONLY : AddEquationBasics, AddEquationSolution, SingleSolver
     
     IMPLICIT NONE

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(Mesh_t),  POINTER  :: Mesh

     TYPE(Element_t),POINTER :: Element
     TYPE(Nodes_t) :: ElementNodes
     TYPE(ValueList_t), POINTER :: BC, Params, RadList
     TYPE(Solver_t), POINTER :: RadSolver
     
     ! parameters for the actual vf/radiator computations
     !---------------------------------------------------
     INTEGER :: divide, LineFlag, LineInteg, TriInteg, QuadInteg, NSymmetry
     REAL(KIND=dp) :: AreaEPS, RayEPS, FactEPS
     INTEGER :: NRays, CombineInt, ClosedFormInt, ShaftStatInt, ClipInt, RayCullInt
     LOGICAL :: Combine, Combine3D, ElimBB
     REAL(KIND=dp), POINTER :: Coord(:)
     INTEGER, POINTER :: Type(:)
     INTEGER, ALLOCATABLE, TARGET ::  Surf(:)
     REAL(KIND=dp), ALLOCATABLE :: Normals(:),Factors(:),Areas(:)

     ! misc variables
     ! --------------
     CHARACTER(:), ALLOCATABLE :: str    
     CHARACTER(LEN=MAX_NAME_LEN) :: ModelName
     LOGICAL :: CylindricSymmetry, GotIt, Found, Radiation

     REAL(KIND=dp) :: MinFactor, Direction, Nrm(3)

     LOGICAL, ALLOCATABLE :: RadiationBC(:)
     LOGICAL :: RadiationOpen, UseSymmetry
     INTEGER :: RadiationBody, MaxRadiationBody
     TYPE(Element_t), POINTER :: RadElements(:)
 
     REAL(KIND=dp) :: at, rt, at2, rt2
     INTEGER :: i,j,k,l,t,n,Ni,istat

     ! MPI row-decomposition for view factor computation
     INTEGER :: nLocal, n_global, iStart_local, myRank, nProcs, vf_comm, mpiErr
     REAL(KIND=dp), ALLOCATABLE :: Factors_local(:)

     ! Radiators on/off, coordinates
     !------------------------------
     INTEGER :: NofRadiators
     LOGICAL :: DoRadiators, DoMapping
     REAL(KIND=dp), POINTER :: Radiators(:,:)

     ! .TRUE. only when InitMPI() really took MPI up.  ELMER_NO_MPI leaves it
     ! down, and MPI must then not be touched again on the way out either.
     LOGICAL :: MPIInitDone = .FALSE.

     INTEGER, PARAMETER :: VFUnit = 10
     CHARACTER(*), PARAMETER :: Caller = 'ViewFactors'
!------------------------------------------------------------------------------
     INTERFACE
        SUBROUTINE RadiatorFactors3D(n, Surf, Type, Coord, Normals, &
           RT_n, RT_Surf, RT_Data, RT_Perm, RT_Type, RT_Coord, &
               NofRadiators, RadiatorCoord, LineFlag, Factors, Feps, &
                 Aeps, Reps, Nr, NInteg2, NInteg3, NInteg4, Combine, ClosedForm, ShaftStat, Clip, RayCull) BIND(C)

          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER, PARAMETER :: dp = 8
          INTEGER  :: n, NofRadiators
          INTEGER :: Surf(*), Type(*)
          REAL(KIND=dp) :: Coord(*), Normals(*), RadiatorCoord(*)

          INTEGER :: RT_n
          INTEGER :: RT_Surf(*), RT_Type(*), RT_Perm(*)
          REAL(KIND=dp) :: RT_Data(*), RT_Coord(*), Factors(*)

          REAL(KIND=dp) :: Feps, Aeps, Reps
          INTEGER :: Nr, NInteg2, NInteg3, NInteg4, Combine, LineFlag
          INTEGER :: ClosedForm, ShaftStat, Clip, RayCull
        END SUBROUTINE RadiatorFactors3D


        SUBROUTINE ViewFactors3D(n, Surf, Type, Coord, Normals, &
          RT_n, RT_Surf, RT_Data, RT_Perm, RT_Type, RT_Coord, &
              Factors, Feps, Aeps, Reps, Nr, NInteg2, NInteg3, NInteg4, Combine, &
              iStart, nLocal, mpiRank, ClosedForm, ShaftStat, Clip, RayCull) BIND(C)

            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            INTEGER, PARAMETER :: dp = 8
            INTEGER :: n
            INTEGER :: Surf(*), Type(*)
            REAL(KIND=dp) :: Coord(*), Normals(*)

            INTEGER :: RT_n
            INTEGER :: RT_Surf(*), RT_Type(*), RT_Perm(*)
            REAL(KIND=dp) :: RT_Coord(*), Factors(*), RT_data(*)

            REAL(KIND=dp) :: Feps, Aeps, Reps
            INTEGER :: Nr, NInteg2, NInteg3, NInteg4, Combine
            INTEGER :: iStart, nLocal, mpiRank, ClosedForm, ShaftStat, Clip, RayCull
        END SUBROUTINE ViewFactors3D


        SUBROUTINE ViewFactorsAxis(n, surf, crd, vf, idiv, fast) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER, PARAMETER :: dp = 8
            REAL(KIND=dp) :: crd(*), vf(*)
            INTEGER :: n, surf(*), idiv, fast
        END SUBROUTINE ViewFactorsAxis
     END INTERFACE
!------------------------------------------------------------------------------

     ! MPI initialisation — must come before any Elmer or MPI calls.
     ! Sets ParEnv%MyPE, ParEnv%PEs and ELMER_COMM_WORLD.
     CALL InitMPI()
     CALL Info( Caller, ' ', Level=3 )
     CALL Info( Caller, '==================================================', Level=3 )
     CALL Info( Caller, ' E L M E R  V I E W F A C T O R S,  W E L C O M E',  Level=3  )
     CALL Info( Caller, '==================================================', Level=3 )

!------------------------------------------------------------------------------
!    Read element definition file, and initialize element types
!------------------------------------------------------------------------------
     CALL InitializeElementDescriptions()
!------------------------------------------------------------------------------
!    Read Model from Elmer Data Base
!------------------------------------------------------------------------------
     CALL ReadCommandLineArguments(ModelName, DoRadiators)
     CALL Info(Caller,'Computing view factors/radiator factors as defined in file: ' &
                //TRIM(ModelName),Level=5)

     CALL InitModel( Model, Mesh, RadSolver )
     Params => GetSolverParams()
     Params => RadSolver % Values

!------------------------------------------------------------------------------

     ! Lets do mapping if requested. 
     CALL SetCurrentMesh( Model,Mesh )
     
     DoMapping = .FALSE.
     ! We can execute the radiation solver in a true way
     DO i=1,Model % NumberOfSolvers
       IF(ListGetLogical(Model % Solvers(i) % Values,'Viewfactor Mapping Solver',Found ) ) THEN           
         DoMapping = .TRUE.
         CALL AddAndExecuteSingleSolver(i)
       END IF
     END DO
     ! Or only use the rigid mesh mapping of the bodies/boundaries 
     IF(.NOT. DoMapping) THEN
       IF(GetLogical( RadSolver % Values,'Viewfactor Rigid Mesh Mapping', Found ) ) THEN
         DoMapping = .TRUE.
         CALL RigidMeshMapping( Model, Mesh, .TRUE. )
       END IF
     END IF
               
!------------------------------------------------------------------------------

     ! Check if computing radiator factors (as opposed to view factors)
     NofRadiators = 0
     IF ( DoRadiators ) THEN
       IF( .NOT. ListCheckPresentAnyBodyForce(Model,'Radiator Coordinates',RadList) ) &
           RadList => Params

       CALL GetConstRealArray(RadList, Radiators, 'Radiator Coordinates', Found)
       IF(.NOT. Found ) CALL Fatal( 'RadiatorFactors', 'No radiators present, quitting' )

       NofRadiators = SIZE(Radiators,1)
       LineFlag = SIZE(Radiators,2) / 6
       CALL Info(Caller, 'Computing flux coefficients for ' // &
               I2S(NofRadiators) // ' radiative sources', LEVEL=5 )
     END IF

     
!------------------------------------------------------------------------------
     CALL SymmetryDuplication(Mesh,Nsymmetry)
!------------------------------------------------------------------------------

     CALL SetCurrentMesh( Model,Mesh )
     CALL Info( Caller,'Number of nodes in mesh: '//&
           I2S(Mesh % NumberOfNodes),Level=7)
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!    Figure out requested coordinate system
!------------------------------------------------------------------------------
     CALL CheckCoordinateSystem(Coordinates)
     CylindricSymmetry = (Coordinates == AxisSymmetric) .OR. &
             (Coordinates==CylindricSymmetric)
!------------------------------------------------------------------------------

     ALLOCATE( ElementNodes % x(Model % MaxElementNodes), &
         ElementNodes % y(Model % MaxElementNodes), &
         ElementNodes % z(Model % MaxElementNodes),STAT=istat )
     
     IF ( CylindricSymmetry ) THEN
       ALLOCATE( Coord(2 * Mesh % NumberOfNodes), STAT=istat )
       DO i=1,Mesh % NumberOfNodes
         Coord(2*(i-1)+1) = Mesh % Nodes % x(i)
         Coord(2*(i-1)+2) = Mesh % Nodes % y(i)
       END DO
     END IF

     IF ( istat /= 0 ) THEN
       CALL Fatal(Caller, 'Memory allocation error. Aborting' )
     END IF

!------------------------------------------------------------------------------
     ! The routine originally plays with the element list and therefore
     ! when several radiation boundaries are needed both the original and
     ! the new elementlist needs to be in the memory. Thus the hassle.
     MinFactor = ListGetConstReal(Params,'Minimum View Factor',GotIt)
     IF(.NOT. GotIt) MinFactor = 1.0d-20

     CALL AllocateVector( RadElements, Mesh % NumberOfBoundaryElements, Caller )

     IF( Mesh % NumberOfBoundaryElements == 0) THEN
       CALL Warn(Caller,'There are no boundary elements at all!')
       STOP
     END IF

     ! Check the maximum radiation body
     MaxRadiationBody = 0

     DO i= 1, Mesh % NumberOfBoundaryElements
       Element => GetBoundaryElement(i)
       IF ( GetElementFamily() == 1 ) CYCLE
       BC => GetBC()
       IF ( .NOT. ASSOCIATED( BC ) ) CYCLE

       str = GetString( BC, 'Radiation',GotIt )
       IF ( GotIt .AND. str == 'diffuse gray' ) THEN
         t = MAX(1, GetInteger( BC, 'Radiation Boundary', GotIt ) )
         MaxRadiationBody = MAX(t, MaxRadiationBody)
       END IF
     END DO
     IF ( MaxRadiationBody==0 ) MaxRadiationBody=1


     IF( Mesh % NumberOfBoundaryElements == 0) THEN
       CALL Warn(Caller,'There are no radiation boundary elements!')
       STOP
     END IF

     ALLOCATE(RadiationBC(Model % NumberOfBCs))
     RadiationBC = .FALSE.
!------------------------------------------------------------------------------

     DO RadiationBody = 1, MaxRadiationBody

       ! Initialise MPI context unconditionally so WriteOutputFile guard
       ! (IF myRank==0) is always valid regardless of geometry path.
       myRank  = ParEnv % MyPE
       nProcs  = ParEnv % PEs
       vf_comm = ParEnv % ActiveComm
       iStart_local = 0
       n_global     = 0

       IF( MaxRadiationBody > 1) THEN
         CALL Info(Caller,'Computing view factors for radiation body' // I2S(RadiationBody), Level=3)
       END IF
       at2 = CPUTime(); rt2 = RealTime()
         
       ! loop to get the surfaces participating in radiation, discard the rest
       ! of the elements...
       !------------------------------------------------------------------------------
       n = ExtractSurfaces(Mesh,DoRadiators,RadElements,RadiationBC, &
                     RadiationBody,RadiationOpen)

       IF ( n==0 ) THEN
         CALL Warn( Caller, 'No surfaces participating in radiation?' )
         IF(RadiationBody < MaxRadiationBody) THEN
           CYCLE
         ELSE
           CALL Warn( Caller, 'Stopping cause nothing to be done...' )
           STOP
         END IF
       END IF
       CALL Info(Caller,'Number of surfaces participating in radiation: '//I2S(N))

       Ni = N
       IF ( DoRadiators ) Ni = NofRadiators
       
       ! Allocate arrays
       ! ---------------
       ALLOCATE( Areas(n), STAT=istat )
       IF ( CylindricSymmetry ) THEN
         ALLOCATE( Surf(2*n), Factors(n*n), STAT=istat )
       ELSE
         ! Factors allocated later, after MPI geometry gather determines n_global
         ALLOCATE( Normals(3*n), Surf(4*n), Type(n), STAT=istat )
       END IF
       IF ( istat /= 0 ) THEN
         CALL Fatal( Caller, 'Memory allocation error. Aborting' )
       END IF
       
       ! go through all surfaces participating in radiation for normal direction:
       ! ------------------------------------------------------------------------
       DO i=1,n
         Element => RadElements(i)
         Model % CurrentElement => Element

         k = GetElementNOFNodes()
         CALL GetElementNodes(ElementNodes)
                  
         Areas(i) = ElementArea(Mesh, Element, Element % Type % NumberOfNodes)
         
         ! Figure out normal direction
         ! ---------------------------

         Nrm = NormalVector( Element, ElementNodes )
         Direction = DirectedNormalVector(Element, Nrm)

         IF ( CylindricSymmetry ) THEN
           IF ( Direction > 0 ) THEN
             Surf(2*(i-1)+1) = Element % NodeIndexes(2)-1
             Surf(2*(i-1)+2) = Element % NodeIndexes(1)-1
           ELSE
             Surf(2*(i-1)+1) = Element % NodeIndexes(1)-1
             Surf(2*(i-1)+2) = Element % NodeIndexes(2)-1
           END IF
         ELSE
           IF ( Direction>0) THEN
             Normals(3*(i-1)+1:3*(i-1)+3) =  Nrm
           ELSE
             Normals(3*(i-1)+1:3*(i-1)+3) = -Nrm
           END IF
         END IF
       END DO
       ! --------------------

       CALL Info( Caller, 'Computing viewfactors/radiator factors...', Level=4 )
       at = CPUTime(); rt = RealTime()
       
       ! Keyword common to cyl symm & cartesian 2d (handled by different codes)
       ! ----------------------------------------------------------------------
       CombineInt = 0
       Combine3D = .FALSE.

       Combine = GetLogical( Params, 'Viewfactor combine elements',GotIt)
       IF( Mesh % MeshDim == 2 ) THEN         
         IF ( .NOT. GotIt ) Combine = .TRUE.
         IF( Combine ) CombineInt = 1
       ELSE
         Combine3D = Combine
       END IF
       
       ElimBB = .NOT. GetLogical( Params,'Viewfactor BBox Shadow', GotIt)
       
       IF ( CylindricSymmetry ) THEN
         ! Axisymmetric case (radiators not implemented).
         ! MPI row decomposition is not yet implemented for this path;
         ! the computation runs on every rank but only rank 0 writes output.
         IF ( nProcs > 1 ) CALL Warn(Caller, &
             'MPI: axisymmetric view factors computed redundantly on all ranks')
         n_global = n
         divide = GetInteger( Params, 'Viewfactor divide',GotIt)
         IF ( .NOT. GotIt ) Divide = 1
         CALL ViewFactorsAxis( N, Surf, Coord, Factors, divide, CombineInt )
       ELSE
         ! Cartesian 2d & 3d, both view factors & radiator factors computed, as needed
         ! ---------------------------------------------------------------------------

         CALL GetCartParameters( AreaEPS, FactEPS, RayEps, Nrays, LineInteg, QuadInteg, &
                        TriInteg, UseSymmetry, ClosedFormInt, ShaftStatInt, ClipInt, RayCullInt )

         Mesh % NumberOfBulkElements = n
         Mesh % Elements => RadElements
         Mesh % NumberOfBoundaryElements = 0

         ! Extract the surfaces from given mesh to arrays passed to the view factor code.
         ! All ranks load the full mesh (LoadModel uses 1,0), so n is identical on all
         ! ranks — no geometry communication needed.
         ! -----------------------------------------------------------------------------
         CALL ExtractMeshInfo( Mesh, n, Coord, Surf, Type )

         ! --- MPI: split computation rows evenly across ranks ---
         ! All ranks have identical geometry; only the C-kernel work is distributed.
         IF( nProcs > 1 ) THEN
           n_global     = n
           nLocal       = n / nProcs + MERGE(1, 0, myRank < MOD(n, nProcs))
           iStart_local = myRank * (n / nProcs) + MIN(myRank, MOD(n, nProcs))
           CALL Info(Caller,'MPI: '//I2S(nProcs)//' ranks, '// &
               'local rows '//I2S(iStart_local)//'..'//I2S(iStart_local+nLocal-1)// &
               ' of '//I2S(n_global), Level=5)
           IF(UseSymmetry) CALL Fatal(Caller,'"Viewfactor use symmetry" not implemented in parallel!')
         ELSE
           n_global = n
           nLocal = n
           iStart_local = 0
           IF(UseSymmetry) THEN
             ! When we use symmetry the only compute the n/2 first rows             
             nLocal = n / 2
           ELSE
             nLocal = n
           END IF            
         END IF

         ! Allocate local Factors for C kernel (nLocal rows × n_global cols)
         ALLOCATE( Factors_local(nLocal * n_global), STAT=istat )
         IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for Factors_local')
         Factors_local = 0.0_dp

         ! In order to speed up shadowing checks, if requested, reduce the mesh complexity by finding
         ! simply connected planar areas, and replace the mesh within by few quads or triangles. Also
         ! potentially finds circular planar areas (disabled atm). Alternatively a surface or volume
         ! mesh may be given from disk.
         !-------------------------------------------------------------------------------------------
         BLOCK
           !------------------------------------------------------
           INTEGER ::  RT_n
           TYPE(Mesh_t), POINTER :: RT_Mesh

           REAL(KIND=dp), ALLOCATABLE :: RT_Data(:)
           INTEGER, POINTER :: RT_Type(:)
           REAL(KIND=dp), POINTER :: RT_Coord(:)
           INTEGER, ALLOCATABLE ::  RT_Surf(:), RT_Perm(:), Ref(:)

           ! Hoisted out of the two branches further down: nested BLOCKs are
           ! miscompiled by some compilers (Intel 20.0), and neither branch
           ! needed a scope of its own.
           INTEGER :: base_n, rem_n, n2
           INTEGER, ALLOCATABLE :: recvcounts(:), displs(:)
           REAL(KIND=dp), ALLOCATABLE :: dummy_recv(:)
           !------------------------------------------------------

           RT_n = 0
           RT_Type  => Type
           RT_Coord => Coord

           IF (Combine3D) THEN
             ! All ranks have the full mesh (LoadModel uses 1,0), so PlanarReduce
             ! is available in MPI mode — every rank produces the same result.
             RT_Mesh => PlanarReduce(n, Normals, Coord, Mesh)
           ELSE
             ! Given surface OR volume shadow mesh from disk
             ! ---------------------------------------------
             RT_Mesh => CheckShadowMeshGiven()
           END IF

           ! Extract possible shadowing surfaces
           ! -----------------------------------
           CALL ExtractMeshInfo( RT_Mesh, RT_n, RT_Coord, RT_Surf, RT_Type, RT_Data, RT_Perm, ElimBBox = ElimBB )

           ! ExtractMeshInfo returns immediately when there is no shadow mesh,
           ! leaving these three unallocated, and the deallocation below shows
           ! Data/Perm need not come back allocated even when there is one.
           ! They are passed to assumed-size dummies of ViewFactors3D and
           ! RadiatorFactors3d, which an unallocated allocatable must not be, so
           ! make sure they exist whether or not they carry anything.
           IF ( .NOT. ALLOCATED(RT_Surf) ) THEN
             ALLOCATE( RT_Surf(1) ); RT_Surf = 0
           END IF
           IF ( .NOT. ALLOCATED(RT_Perm) ) THEN
             ALLOCATE( RT_Perm(1) ); RT_Perm = 0
           END IF
           IF ( .NOT. ALLOCATED(RT_Data) ) THEN
             ALLOCATE( RT_Data(1) ); RT_Data = 0.0_dp
           END IF

           ! Shadow mesh is fully available on all ranks:
           !  - LoadShadowMesh: reads mesh files directly (serial, no MPI)
           !  - LoadMesh2: called with (1,0) above so all ranks load the full mesh
           ! No gather needed.

           IF ( RT_n > 0 ) THEN
             CALL Info(Caller,'Using separate mesh for shadowing, #elements = '//I2S(RT_n),Level=5)

             WRITE (Message,'(A,2F8.2)') 'Shadow mesh defined time (s):',&
                 CPUTime()-at2, Realtime()-rt2
             CALL Info( Caller,Message, Level=3 )
             at2 = CPUTime(); rt2 = RealTime()

             ! This is mean mainly for debugging & creating pictures of the joining algorithm.
             ! End-users will use this probably very rarely.
             IF ( GetLogical( Params,'Shadow Mesh Save', GotIt) ) THEN
               str = "ShadowMesh"
               CALL MakeDirectory(TRIM(str) // CHAR(0))
               CALL WriteMeshToDisk2(Model, RT_Mesh, str )               
               CALL Info(Caller,'Saved shadow mesh to file: '//TRIM(str))
               IF ( GetLogical( Params,'Shadow Mesh Stop', GotIt) ) STOP
             END IF
             
           END IF

           ! ... and finally the beef:
           ! -------------------------
           IF ( DoRadiators ) THEN
             ! Radiators: no MPI row decomposition yet; serial path unchanged
             ALLOCATE( Factors(NofRadiators * n_global), STAT=istat )
             IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for Factors')
             CALL RadiatorFactors3d( n, Surf, TYPE, Coord, Normals, RT_n, RT_Surf, &
                  RT_Data, RT_Perm, RT_Type, RT_Coord, NofRadiators, Radiators, LineFlag, &
                       Factors, AreaEPS, FactEPS, RayEPS, Nrays, LineInteg, TriInteg, QuadInteg, CombineInt, &
                         ClosedFormInt, ShaftStatInt, ClipInt, RayCullInt )
           ELSE
             CALL ViewFactors3D( n, Surf, Type, Coord, Normals, RT_n, RT_Surf, &
                  RT_Data, RT_Perm, RT_Type, RT_Coord, Factors_local, AreaEPS, FactEPS, RayEPS, &
                      Nrays, LineInteg, TriInteg, QuadInteg, CombineInt, &
                      iStart_local, nLocal, myRank, ClosedFormInt, ShaftStatInt, ClipInt, RayCullInt )

             ! RayCullInt comes back resolved: the automatic setting is only
             ! decided in InitStuff, which is where the shadow mesh size is
             ! finally known.
             IF ( RayCullInt == 1 ) &
               CALL Info(Caller,'Culling the rays of each pair to its shaft candidates',Level=8)

             ! Gather local rows to rank 0 only — only rank 0 needs the full
             ! matrix for normalisation and output.  Non-root ranks send their
             ! rows but receive nothing, saving O(N²) memory and Newton work.
             IF ( nProcs > 1 ) THEN
               ALLOCATE( recvcounts(0:nProcs-1), displs(0:nProcs-1) )
               base_n = n_global / nProcs;  rem_n = MOD(n_global, nProcs)
               DO i = 0, nProcs-1
                 recvcounts(i) = (base_n + MERGE(1,0,i<rem_n)) * n_global
               END DO
               displs(0) = 0
               DO i = 1, nProcs-1
                 displs(i) = displs(i-1) + recvcounts(i-1)
               END DO
               IF ( myRank == 0 ) THEN
                 ALLOCATE( Factors(n_global * n_global), STAT=istat )
                 IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for Factors')
                 CALL MPI_Gatherv( Factors_local, nLocal*n_global, MPI_DOUBLE_PRECISION, &
                     Factors, recvcounts, displs, MPI_DOUBLE_PRECISION, 0, vf_comm, mpiErr )
               ELSE
                 ALLOCATE( dummy_recv(1) )  ! receive buffer ignored on non-root
                 CALL MPI_Gatherv( Factors_local, nLocal*n_global, MPI_DOUBLE_PRECISION, &
                     dummy_recv, recvcounts, displs, MPI_DOUBLE_PRECISION, 0, vf_comm, mpiErr )
               END IF
               DEALLOCATE( recvcounts, displs )
             ELSE
               IF( UseSymmetry ) THEN
                 ! Do not ever create the full matrix but directly reduce the size to half.
                 n2 = n_global/2
                 ALLOCATE( Factors(n2*n2), STAT=istat )
                 IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for Factors')
                 DO i=1,n2
                   DO j=1,n2
                     Factors((i-1)*n2+j) = Factors_local((i-1)*n+j) + Factors_local((i-1)*n+j+n2)
                   END DO
                 END DO
                 n = n2
                 ni = n2
               ELSE
                 ALLOCATE( Factors(n_global * n_global), STAT=istat )
                 IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for Factors')
                 Factors = Factors_local
               END IF
             END IF
             DEALLOCATE( Factors_local )
             ! IterSolv uses ipar=0/dProc=0 → sequential BLAS dot products,
             ! no MPI.  Non-root ranks skip the IF(myRank==0) normalisation
             ! block below and loop back; the next body's MPI_Gatherv acts as
             ! the natural synchronisation point, so no explicit barrier needed.
           END IF

           IF (RT_n>0) THEN
!            CALL FreeMesh(RT_Mesh)
             DEALLOCATE(RT_Surf, RT_Type, RT_Coord)
             IF (ALLOCATED(RT_data)) DEALLOCATE(RT_data, RT_perm)
           END IF
         END BLOCK
       END IF  ! CylindricSymmetry
       
       WRITE (Message,'(A,2F8.2)') 'View factors/radiator factors computed in time (s):',&
           CPUTime()-at2, Realtime()-rt2
       CALL Info( Caller,Message, Level=3 )
       at2 = CPUTime(); rt2 = RealTime()

       ! Only rank 0 has the full Factors matrix — normalisation and output
       ! are rank-0-only operations.  Other ranks already did their work.
       IF ( myRank == 0 ) THEN
         IF(.NOT. UseSymmetry) THEN
           CALL SymmetryReduction(DoRadiators,NofRadiators,n,Ni,Factors)
         END IF
           
         CALL FindInitialMinMax(Ni,N,Factors,RadiationOpen)
         CALL NormalizeFactors(Model,DoRadiators,NofRadiators,n,Factors,RadiationOpen)
         CALL FindNormalizedMinMax(Ni,n,Factors)

         WRITE (Message,'(A,2F8.2)') 'View factors manipulated in time (s):',&
             CPUTime()-at2, realtime()-rt2
         CALL Info( Caller,Message, Level=3 )
         at2 = CPUTime(); rt2 = RealTime()

         IF(InfoActive(12)) CALL ViewFactorsLumping()
         CALL WriteOutputFile(DoRadiators,Ni,n,Factors,RadiationBody)

         WRITE (Message,'(A,2F8.2)') 'View factors saved in time (s):',&
             CPUTime()-at2, realtime()-rt2
         CALL Info( Caller,Message, Level=3 )
         at2 = CPUTime(); rt2 = RealTime()
       END IF

       DEALLOCATE( Surf, Areas )
       IF ( ALLOCATED(Factors) ) DEALLOCATE( Factors )
       IF ( .NOT. CylindricSymmetry ) DEALLOCATE(Normals, Type)

     END DO  ! Of radiation RadiationBody

     WRITE (Message,'(A,2F8.2)') 'View factors all done in time (s):',&
         CPUTime()-at, realtime()-rt
     CALL Info( Caller,Message, Level=3 )

     CALL FLUSH(6)
     ! Skipped under ELMER_NO_MPI: ParEnvFinalize() opens with an
     ! MPI_BARRIER on ELMER_COMM_WORLD, and handing MPI_Comm_f2c a
     ! communicator before MPI_Init aborts the process.  That abort came
     ! after the factors were computed and written, so with the exit status
     ! of the spawn discarded it used to leave no trace at all.
     IF ( MPIInitDone ) CALL ParallelFinalize()

CONTAINS



!------------------------------------------------------------------------------
!> Adds flags for active solvers. 
!------------------------------------------------------------------------------
  SUBROUTINE AddAndExecuteSingleSolver(Solver_ind)
!------------------------------------------------------------------------------
    INTEGER :: Solver_ind

    INTEGER :: i,j,nlen 
    LOGICAL :: Found
    INTEGER, POINTER :: ActiveSolvers(:)
    CHARACTER(:), ALLOCATABLE :: eq
    REAL(KIND=dp), POINTER :: sTime(:), sStep(:)
    TYPE(Solver_t), POINTER :: Solver

    CALL Info(Caller,'Running Viewfactor Mapping Solver: '//I2S(Solver_ind),Level=5)
    
    CALL Info(Caller,'Setting mesh coordinates and time',Level=12)

    CALL SetCurrentMesh( Model,Mesh )
    
    CALL VariableAdd( Mesh % Variables, Mesh, &
        Name='Coordinate 1',DOFs=1,Values=Mesh % Nodes % x )    
    CALL VariableAdd(Mesh % Variables,Mesh, &
        Name='Coordinate 2',DOFs=1,Values=Mesh % Nodes % y )    
    CALL VariableAdd(Mesh % Variables,Mesh, &
        Name='Coordinate 3',DOFs=1,Values=Mesh % Nodes % z )
         
    ALLOCATE(sTime(1),sStep(1))
    sTime = 0.0_dp
    sStep = 1.0_dp    
    CALL VariableAdd( Mesh % Variables, Mesh, Name='Time',DOFs=1, Values=sTime )
    CALL VariableAdd( Mesh % Variables, Mesh, Name='Timestep', DOFs=1, Values=sStep )


    CALL Info('AddSolver','Setting up solver: '//I2S(Solver_ind),Level=10)
    Solver => Model % Solvers(Solver_ind)

    eq = ListGetString( Solver % Values,'Equation', Found )     
    IF ( Found ) THEN
      nlen = LEN_TRIM(eq)
      DO j=1,Model % NumberOFEquations
        ActiveSolvers => ListGetIntegerArray( Model % Equations(j) % Values, &
            'Active Solvers', Found )
        IF ( Found ) THEN
          DO k=1,SIZE(ActiveSolvers)
            IF ( ActiveSolvers(k) == Solver_ind ) THEN
              CALL ListAddLogical( Model % Equations(j) % Values, eq(1:nlen), .TRUE. )
              EXIT
            END IF
          END DO
        END IF
      END DO
    END IF

    Solver % Mesh => Mesh         
    Model % Solver => Solver
    CALL AddEquationBasics( Solver, eq, Transient = .FALSE.)
    CALL AddEquationSolution( Solver, Transient = .FALSE.)

    CALL Info(Caller,'Executing solver '//I2S(Solver_ind)//' within ViewFactors!,Level=10')
    CALL SingleSolver( Model, Solver, 0.0_dp, .FALSE. )

    CALL Info(Caller,'Execution of mapping solver done',Level=12)

!------------------------------------------------------------------------------
  END SUBROUTINE AddAndExecuteSingleSolver
!------------------------------------------------------------------------------


  
!------------------------------------------------------------------------------
!> Thin wrapper around ParallelInit() so TYPE(ParEnv_t) stays out of the
!> PROGRAM specification section (avoids parser confusion near INTERFACE).
!------------------------------------------------------------------------------
   SUBROUTINE InitMPI()
     TYPE(ParEnv_t), POINTER :: penv
     CHARACTER(LEN=4) :: val
     INTEGER :: vlen, vstat
     ! When spawned from within ElmerSolver via SystemCommand, the subprocess
     ! inherits the parent's PMI/PMIx environment, which can cause MPI_Init to
     ! hang waiting for a rendezvous that never comes.  The caller sets
     ! ELMER_NO_MPI=1 to request serial-only operation in this case.
     CALL GET_ENVIRONMENT_VARIABLE('ELMER_NO_MPI', val, vlen, vstat)
     IF (vstat == 0 .AND. vlen > 0) THEN
       ParEnv % MyPE = 0
       ParEnv % PEs  = 1
       ParEnv % ActiveComm = 0
       ParEnv % ExternalInit = .FALSE.
       MPIInitDone = .FALSE.
     ELSE
       penv => ParallelInit()
       MPIInitDone = .TRUE.
     END IF
   END SUBROUTINE InitMPI

!------------------------------------------------------------------------------
   SUBROUTINE InitModel(Model,Mesh,RadSolver)
!------------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(Mesh_t),  POINTER :: Mesh
     TYPE(Solver_t), POINTER :: RadSolver
!------------------------------------------------------------------------------
     INTEGER :: i
     TYPE(Solver_t), POINTER  :: Solver
!------------------------------------------------------------------------------

     CALL Info( Caller, ' ', Level=3 )
     CALL Info( Caller, ' ', Level=3 )
     CALL Info( Caller, 'Reading Model... ', Level=3 )
!------------------------------------------------------------------------------
     Model => LoadModel( ModelName,.FALSE.,1,0 )
     CurrentModel => Model     
!------------------------------------------------------------------------------
          
     Mesh => NULL()
     RadSolver => NULL()
     
     DO i=1,Model % NumberOfSolvers
       Solver => Model % Solvers(i)
       Radiation = ListGetLogical( Solver % Values, 'Radiation Solver', Found )
       IF ( Radiation ) THEN
         CALL Info(Caller,'Radiation treated by solver '//I2S(i),Level=10)
         Mesh => Solver % Mesh
         Model % Solver => Solver
         RadSolver => Solver
         EXIT
       ENDIF
     END DO
          
     IF ( .NOT. ASSOCIATED(Mesh) ) THEN
       CALL Fatal(Caller,'No heat equation definition. Cannot compute factors.')
     END IF

!------------------------------------------------------------------------------
   END SUBROUTINE InitModel
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE ReadCommandlineArguments(ModelName,DoRadiators)
!------------------------------------------------------------------------------
     CHARACTER(*) :: ModelName
     LOGICAL :: DoRadiators
!------------------------------------------------------------------------------
     INTEGER :: NoArgs, iostat
!------------------------------------------------------------------------------
     NoArgs = COMMAND_ARGUMENT_COUNT()

     DoRadiators = .FALSE.
     IF ( NoArgs > 0 ) THEN
       CALL GET_COMMAND_ARGUMENT(1, ModelName)
       IF ( ModelName(1:10) == '-radiators' ) THEN
          DoRadiators = .TRUE.
       END  IF
     END IF

     ! A single argument is the sif to read.  Without this it fell through to
     ! ELMERSOLVER_STARTINFO and the name given on the command line was
     ! silently ignored, so "ViewFactors a.sif" and "ViewFactors b.sif" both
     ! computed whichever case STARTINFO happened to name -- wrong, and quiet
     ! about it.  "-radiators" on its own is a flag, not a model name.
     IF ( NoArgs > 1 ) THEN
       CALL GET_COMMAND_ARGUMENT(NoArgs, ModelName)
     ELSE IF ( NoArgs == 1 .AND. .NOT. DoRadiators ) THEN
       CONTINUE                       ! ModelName already holds the argument
     ELSE
       OPEN( 1,file='ELMERSOLVER_STARTINFO', STATUS='OLD', IOSTAT=iostat )
       IF( iostat /= 0 ) THEN
         CALL Fatal( Caller, 'Unable to find ELMERSOLVER_STARTINFO, cannot execute.' )
       END IF
       READ(1,'(a)') ModelName
       CLOSE(1)
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE ReadCommandlineArguments
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CheckCoordinateSystem(Coordinates)
!------------------------------------------------------------------------------
     INTEGER::  Coordinates
!------------------------------------------------------------------------------
     CHARACTER(:), ALLOCATABLE :: eq
!------------------------------------------------------------------------------
     eq = GetString( GetSimulation(), 'Coordinate System' )
     SELECT CASE(eq)
     CASE('cartesian','cartesian 2d','cartesian 3d')
       Coordinates = Cartesian
     CASE('axi symmetric')
       Coordinates = AxisSymmetric
     CASE('cylindric symmetric')
       Coordinates = CylindricSymmetric
     CASE DEFAULT
       CALL Error( Caller, &
         'Unknown Global Coordinate System for Viewfactor computation ')
       CALL Error( Caller, TRIM(eq) )
       CALL Fatal( Caller, &
         'Only Cartesian 3D or Axi/Cylindrical Symmetric coordinates allowed. Aborting' )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE CheckCoordinateSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Provide useful information on the boundary-to-boundary view factors that is
! obtained as the area-weigted average of elemental view factor sums.
! In case of radiators provide summed up radiator-to-boundary factors.   
!------------------------------------------------------------------------------
   SUBROUTINE ViewFactorsLumping()
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,m,MaxRadBC,bc_id,NoRows
     INTEGER, ALLOCATABLE :: BCNumbering(:), VFPerm(:)
     REAL(KIND=dp), ALLOCATABLE :: LumpedVF(:,:), LumpedAreas(:)
!------------------------------------------------------------------------------
      
     ALLOCATE(BCNumbering(Model % NumberOfBCs),VFPerm(N))
     BCNumbering = 0
     VFPerm = 0
              
     DO i=1,N
       Element => RadElements(i)
       bc_id = GetBCId( Element ) 
       IF(bc_id < 0 .OR. bc_id > Model % NumberOfBCs) THEN
         CALL Warn(Caller,'BC index out of bounds: '//I2S(bc_id))
         CYCLE
       END IF
       BCNumbering(bc_id) = BCNumbering(bc_id) + 1
       VFPerm(i) = bc_id
     END DO

     j = 0
     DO i=1,Model % NumberOfBCs
       m = BCNumbering(i)
       IF(m>0) THEN
         j = j+1
         BCNumbering(i) = j
         CALL Info(Caller,'BC '//I2S(i)//' with '//I2S(m)//' elems perm: '//I2S(j))
       END IF         
     END DO
     MaxRadBC = j

     DO i=1,N
       VFPerm(i) = BCNumbering(VFPerm(i))
     END DO

     IF(DoRadiators) THEN
       CALL Info(Caller,'Printing lumped information on radiator factors')
       NoRows = NofRadiators       
       ALLOCATE(LumpedVF(NoRows,MaxRadBC))
       LumpedVF = 0.0_dp

       DO i=1,NoRows
         DO j=1,N
           k = (i-1)*N+j
           IF( VFPerm(j) > 0) THEN
             LumpedVF(i,VFPerm(j)) = LumpedVF(i,VFPerm(j)) + Factors(k)
           END IF
         END DO
       END DO
       CALL Info(Caller,'Lumped Radiator Matrix:')
     ELSE
       CALL Info(Caller,'Printing lumped information of view factors')
       NoRows = MaxRadBC
       ALLOCATE(LumpedVF(MaxRadBC,MaxRadBC),LumpedAreas(MaxRadBC))             
       LumpedVF = 0.0_dp
       LumpedAreas = 0.0_dp
       
       DO i=1,N
         IF(VFPerm(i) > 0) THEN
           LumpedAreas(VFPerm(i)) = LumpedAreas(VFPerm(i)) + Areas(i)
         END IF
       END DO
         
       CALL Info(Caller,'Lumped areas:')
       WRITE(Message,*) LumpedAreas(:)
       CALL Info(Caller, Message ) 
       
       DO i=1,N
         DO j=1,N
           k = (i-1)*N+j
           IF( VFPerm(i) > 0 .AND. VFPerm(j) > 0) THEN
             LumpedVF(VFPerm(i),VFPerm(j)) = LumpedVF(VFPerm(i),VFPerm(j)) + Areas(i) * Factors(k)
           END IF
         END DO
       END DO
       
       DO i=1,MaxRadBC
         DO j=1,MaxRadBC
           LumpedVF(i,j) = LumpedVF(i,j) / LumpedAreas(i)
         END DO
       END DO
       CALL Info(Caller,'Lumped View Factor Matrix:')
     END IF
     
     DO i=1,NoRows
       WRITE(Message,*) LumpedVF(i,:)
       CALL Info(Caller, Message ) 
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE ViewFactorsLumping
!------------------------------------------------------------------------------


     
!> View factors are normalized in order to improve the numerical accuracy. With 
!> normalization it is ensured that all boundary elements see exactly half 
!> space. 
!------------------------------------------------------------------------------
   SUBROUTINE NormalizeFactors( Model, DoRadiators, NofRadiators, &
                    N, Factors, RadiationOpen )
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: NofRadiators, N
     REAL(KIND=dp)  :: Factors(:)
     LOGICAL :: DoRadiators, RadiationOpen
!------------------------------------------------------------------------------
     INTEGER :: itmax,it,i,j,k
     LOGICAL :: li,lj
     REAL(KIND=dp), ALLOCATABLE :: RHS(:),SOL(:),PSOL(:)
     REAL(KIND=dp) :: cum,s,si,sj
     REAL(KIND=dp) :: at1
     REAL(KIND=dp), PARAMETER :: eps=EPSILON(1._dp)
!------------------------------------------------------------------------------

     IF( RadiationOpen ) THEN
       CALL Info( Caller,'Symmetrizing Factors... ')
     ELSE
       CALL Info( Caller,'Normalizing Factors...')
     END IF

     IF ( DoRadiators ) THEN
       IF( .NOT. RadiationOpen ) THEN
!$omp parallel do private(s,i,j)
         DO i=1,NofRadiators
           s = 0.0_dp
           DO j=1,n
             s = s + Factors((i-1)*N+j)
           END DO

           DO j=1,n
             Factors((i-1)*N+j) = Factors((i-1)*N+j) / s
           END DO
         END DO
!$omp end parallel do
       END IF
       RETURN
     END IF

     itmax = 20
     it = 0
     cum = 0.0_dp
        
!------------------------------------------------------------------------------
!    First force the matrix (before dividing by area) to be symmetric
!------------------------------------------------------------------------------
        
!$omp parallel do private(s,si,sj,li,lj,i,j)
     DO i=1,n
       DO j=i,n
         si = Areas(i) * Factors((i-1)*n+j)
         sj = Areas(j) * Factors((j-1)*n+i)

         li = (ABS(si) < HUGE(si)) 
         lj = (ABS(sj) < HUGE(sj)) 

         IF(li .AND. lj) THEN 
           s = (si+sj)/2.0
         ELSE IF(li) THEN
           s = si
         ELSE IF(lj) THEN
           s = sj
         ELSE 
           s = 0.0
         END IF

         Factors((i-1)*n+j) = s
         Factors((j-1)*n+i) = s
       END DO
     END DO
!$omp end parallel do

!------------------------------------------------------------------------------
!    Next we solve the equation DFD = A by Newton iteration (this is a very
!    well behaved equation (symmetric, diagonal dominant), no need for any
!    tricks...)
!------------------------------------------------------------------------------
     IF(.NOT. RadiationOpen ) THEN
       
       ALLOCATE( RHS(n),SOL(n),PSOL(n),Jdiag(n),Jacobian(n,n),STAT=istat )
       IF ( istat /= 0 ) THEN
         CALL Fatal( Caller,'Memory allocation error in NormalizeFactors for RHS etc.' )
       END IF

       SOL = 1.0_dp
       cum = 1.0_dp
       
       DO it=1,itmax
!$omp parallel do private(cum,i,j)
         DO i=1,n
           cum = 0.0_dp
           DO j=1,n
             cum = cum + Factors((i-1)*n+j) * SOL(j)
           END DO
           cum = cum * SOL(i)
           RHS(i) = Areas(i) - cum
         END DO
!$omp end parallel do

         cum = SUM( RHS*RHS/Areas(1:n) ) / n
            
         WRITE (Message,'(A,ES12.3)') &
             'Normalization iteration '//I2S(it)//': ',cum;
         CALL Info( Caller,Message, Level=3 );
            
         IF ( cum <= eps ) EXIT
         
!$omp parallel do private(i,j)
         DO i=1,n
           DO j=1,n
             Jacobian(i,j) = Factors((i-1)*n+j) * SOL(i)
           END DO
           DO j=1,n
             Jacobian(i,i) = Jacobian(i,i) + Factors((i-1)*n+j) * SOL(j)
           END DO
           Jdiag(i) = 1._dp / Jacobian(i,i)
         END DO
!$omp end parallel do

         PSOL = SOL
         CALL IterSolv( n,SOL,RHS )
         SOL = PSOL + SOL
       END DO
          
!------------------------------------------------------------------------------
!    Normalize the factors and (re)divide by areas
!------------------------------------------------------------------------------
!$omp parallel do private(i,j)
       DO i=1,N
         DO j=1,N
           Factors((i-1)*N+j) = Factors((i-1)*N+j)*SOL(i)*SOL(j)/Areas(i)
         END DO
       END DO
!$omp end parallel do
       DEALLOCATE( SOL,RHS,PSOL,Jdiag,Jacobian )
     ELSE
!$omp parallel do private(i,j)
       DO i=1,N
         DO j=1,N
           Factors((i-1)*N+j) = Factors((i-1)*N+j)/Areas(i)
         END DO
       END DO
!$omp end parallel do
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE NormalizeFactors
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE MirrorMesh(Mesh,c,Plane,NoDoubles)
!------------------------------------------------------------------------------
     IMPLICIT NONE

!------------------------------------------------------------------------------
     TYPE(Mesh_t) :: Mesh
     INTEGER :: c
     REAL(KIND=dp) :: Plane
     LOGICAL :: NoDoubles
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: el(:)
     INTEGER :: i,j,ne, nd, nn, nb, nv, ns, nd2
     REAL(KIND=dp), POINTER :: ox(:), oy(:), oz(:)
     INTEGER, ALLOCATABLE :: SymPerm(:)
     LOGICAL, ALLOCATABLE :: SymNode(:)
     REAL(KIND=dp) :: Eps
!------------------------------------------------------------------------------

     nv  = Mesh % NumberOfBulkElements
     nb  = Mesh % NumberOfBoundaryElements
     ne  = nb+nv

     nd = Mesh % NumberOfNodes
     ox => Mesh % Nodes % x
     oy => Mesh % Nodes % y
     oz => Mesh % Nodes % z

     CALL Info('MirrorMesh','Mirroring mesh aroud axsis: '//I2S(c),Level=12)

     ! Mirroring is easier when we don't need to eliminate the double nodes at the symmetry axis.
     ! However, if we don't eliminate them, the shadow mesh inverstigation fails to recognize
     ! continuity of elements over symmetry axis as they do not share any node.
     ! Hence this flag could save some resources. 
     !---------------------------------------------------------------------------------------
     IF(NoDoubles) THEN
       ALLOCATE(SymNode(nd),SymPerm(nd))
       SymNode = .FALSE.
       Eps = 1.0e-8              
       SELECT CASE(c)
       CASE(1,2)
         SymNode = ( ABS( ox(1:nd) - Plane) < eps ) 
       CASE(3,4)
         SymNode = ( ABS( oy(1:nd) - Plane) < eps ) 
       CASE(5,6)
         SymNode = ( ABS( oz(1:nd) - Plane) < eps ) 
       END SELECT

       ns = COUNT(SymNode)
       nd2 = nd - ns

       j = nd
       DO i=1,nd
         IF(SymNode(i)) THEN
           SymPerm(i) = i
         ELSE
           j = j+1
           SymPerm(i) = j
         END IF
       END DO

       CALL Info(Caller,'Number of symmetry nodes to eliminate: '//I2S(ns))
     ELSE
       nd2 = nd
     END IF

     ALLOCATE(Mesh % Nodes % x(nd+nd2), Mesh % Nodes % y(nd+nd2), Mesh % Nodes % z(nd+nd2), STAT=istat )
     IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation for MirroMesh nodes.')
     
       
     DO i=1,nd
       Mesh % Nodes % x(i) = ox(i)
       Mesh % Nodes % y(i) = oy(i)
       Mesh % Nodes % z(i) = oz(i)

       IF(NoDoubles) THEN
         j = SymPerm(i)
         IF(j==i) CYCLE
       ELSE
         j = i+nd
       END IF

       SELECT CASE(c)
       CASE(1,2)
         Mesh % Nodes % x(j) = 2*Plane - ox(i)
         Mesh % Nodes % y(j) = oy(i)
         Mesh % Nodes % z(j) = oz(i)
       CASE(3,4)
         Mesh % Nodes % x(j) = ox(i)
         Mesh % Nodes % y(j) = 2*Plane - oy(i)
         Mesh % Nodes % z(j) = oz(i)
       CASE(5,6)
         Mesh % Nodes % x(j) = ox(i)
         Mesh % Nodes % y(j) = oy(i)
         Mesh % Nodes % z(j) = 2*Plane - oz(i)
       END SELECT
     END DO

     el => Mesh % Elements
     ALLOCATE(Mesh % Elements(2*ne), STAT=istat)
     IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation for MirroMesh elements.')
     
     DO i=1,nv
       Mesh % Elements(i)    = el(i)
       Mesh % Elements(i+nv) = el(i)
       nn = el(i) % Type % NumberOfNodes
       ALLOCATE(Mesh % Elements(i+nv) % NodeIndexes(nn),STAT=istat)
       IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation for MirroMesh node indexes.')

       IF(NoDoubles) THEN
         Mesh % Elements(i+nv) % NodeIndexes = SymPerm(el(i) % NodeIndexes)
       ELSE
         Mesh % Elements(i+nv) % NodeIndexes = el(i) % NodeIndexes+nd
       END IF
       Mesh % Elements(i+nv) % ElementIndex = i+nv
     END DO

     DO i=nv+1,ne
       j = i+nv
       Mesh % Elements(j)    = el(i)
       Mesh % Elements(j+nb) = el(i)
       nn = el(i) % TYPE % NumberOfNodes
       
       ALLOCATE(Mesh % Elements(j+nb) % NodeIndexes(nn),STAT=istat)
       IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation for MirroMesh NodeIndexes.')         
       IF(NoDoubles) THEN
         Mesh % Elements(j+nb) % NodeIndexes = SymPerm(el(i) % NodeIndexes)
       ELSE
         Mesh % Elements(j+nb) % NodeIndexes = el(i) % NodeIndexes+nd
       END IF
         
       ALLOCATE(Mesh % Elements(j) % BoundaryInfo)
       Mesh % Elements(j) % BoundaryInfo    = el(i) % BoundaryInfo

       ALLOCATE(Mesh % Elements(j+nb) % BoundaryInfo)
       Mesh % Elements(j+nb) % BoundaryInfo = el(i) % BoundaryInfo

       IF(ASSOCIATED(Mesh % Elements(j) % BoundaryInfo % Left)) THEN
         Mesh % Elements(j+nb) % BoundaryInfo % Left => &
             Mesh % Elements(el(i) % BoundaryInfo % Left % ElementIndex+nv)
       END IF

       IF(ASSOCIATED(Mesh % Elements(j) % BoundaryInfo % Right)) THEN
         Mesh % Elements(j+nb) % BoundaryInfo % Right => &
             Mesh % Elements(el(i) % BoundaryInfo % Right % ElementIndex+nv)
       END IF
     END DO

     DEALLOCATE(ox,oy,oz)

     Mesh % NumberOfNodes = nd + nd2
     Mesh % NumberOfBulkElements = 2*Mesh % NumberOfBulkElements
     Mesh % NumberOfBoundaryElements = 2*Mesh % NumberOfBoundaryElements

     Model % NumberOfNodes = Mesh % NUmberOfNodes
     Model % NumberOfBulkElements = Mesh % NumberOfBulkElements
     Model % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------
   END SUBROUTINE MirrorMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE SymmetryDuplication(Mesh,NSymmetry)
!------------------------------------------------------------------------------
     TYPE(Mesh_t) :: Mesh
     INTEGER, OPTIONAL :: NSymmetry
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Plane
     INTEGER :: i
     LOGICAL :: Found, GotIt, NoDoubles
!------------------------------------------------------------------------------
     IF(PRESENT(NSymmetry)) Nsymmetry = 0

     DO i=1,6
       SELECT CASE(i)
       CASE(1)
         Plane = GetCReal( Params, 'Viewfactor Symmetry x min', Found );
         IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry x', GotIt );
       CASE(2)
         Plane = GetCReal( Params, 'Viewfactor Symmetry x max', Found );
       CASE(3)
         Plane = GetCReal( Params, 'Viewfactor Symmetry y min', Found );
         IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry y', GotIt );
       CASE(4)
         Plane = GetCReal( Params, 'Viewfactor Symmetry y max', Found );
       CASE(5)
         Plane = GetCReal( Params, 'Viewfactor Symmetry z min', Found );
         IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry z', GotIt );
       CASE(6)
         Plane = GetCReal( Params, 'Viewfactor Symmetry z max', Found );
       END SELECT

       IF(.NOT. Found ) CYCLE

       NoDoubles = ListGetLogical( Params,'Viewfactor Symmetry Eliminate Nodes', GotIt )
       
       CALL Info(Caller,'Duplicating mesh in coordinate direction: '//I2S((i+1)/2))
       IF(PRESENT(NSymmetry)) NSymmetry = NSymmetry + 1
       
       CALL MirrorMesh(Mesh, i, Plane, NoDoubles )
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE SymmetryDuplication
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE SymmetryReduction(DoRadiators,NofRadiators,N,Ni,Factors)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Factors(:)
     LOGICAL :: DoRadiators
     INTEGER :: NofRadiators, N, Ni
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Plane
     INTEGER :: i,j,k,l,m
     LOGICAL :: Found, GotIt
!------------------------------------------------------------------------------
     DO l=6,1,-1
       SELECT CASE(l)
       CASE(1)
         Plane = GetCReal( Params, 'Viewfactor Symmetry x min', Found );
         IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry x', GotIt );
         ! Note that Plane is zero if 1st keyword not found!
       CASE(2)
         Plane = GetCReal( Params, 'Viewfactor Symmetry x max', Found );
       CASE(3)
         Plane = GetCReal( Params, 'Viewfactor Symmetry y min', Found );
         IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry y', GotIt );
       CASE(4)
         Plane = GetCReal( Params, 'Viewfactor Symmetry y max', Found );
       CASE(5)
         Plane = GetCReal( Params, 'Viewfactor Symmetry z min', Found );
         IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry z', GotIt );
       CASE(6)
         Plane = GetCReal( Params, 'Viewfactor Symmetry z max', Found );
       END SELECT
       IF(.NOT.Found) CYCLE

       CALL Info(Caller,'Symmetry reduction in coordinate direction: '//I2S((l+1)/2))
       
       k = 0
       IF (DoRadiators) THEN
         m = NofRadiators
       ELSE
         m = n/2
       END IF

       DO i=1,m
         DO j=1,n/2
           k = k + 1
           Factors(k) = Factors((i-1)*n+j) + Factors((i-1)*n+j+n/2)
         END DO
       END DO
       n = n/2
       IF ( .NOT. DoRadiators ) Ni=n
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE SymmetryReduction
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION DirectedNormalVector(Element, Nrm) RESULT(dir)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: Element
     REAL(KIND=dp) :: Nrm(3), dir
!------------------------------------------------------------------------------
     ! normal direction
     !-----------------
     REAL(KIND=dp) :: Normal_in, Nrm2(3), r1(3), r2(3)
     LOGICAL :: Lrad, Rrad
     INTEGER :: Lnode,Rnode,Lbody,Rbody,RadBody,Inode

     LOGICAL :: GotIt

     TYPE(ValueList_t), POINTER :: BC
!------------------------------------------------------------------------------
     CALL GetParentInfo(Element % BoundaryInfo % Left,  Element, Lbody, Lnode, Lrad)
     CALL GetParentInfo(Element % BoundaryInfo % Right, Element, Rbody, Rnode, Rrad)

     BC => GetBC()

     RadBody = GetInteger( BC, 'Radiation Target Body',GotIt )
     IF ( .NOT. GotIt ) RadBody = GetInteger( BC, 'Normal Target Body',GotIt )

     IF ( .NOT. Gotit ) THEN
       IF ( Lrad .AND. Rrad ) THEN
         WRITE(Message,'(A,I3,I3,A,I3,A,I5)') 'Bodies:',Rbody,Lbody,' BC:', &
                          GetBCId( Element ),' Ind:',i
         CALL Info(Caller,Message)
         CALL Fatal(Caller,'"Radiative Body" defined on both sides!')
       ELSE IF (Lrad) THEN
         RadBody = Rbody
       ELSE IF (Rrad) THEN
         RadBody = Lbody
       ELSE
         RadBody = 0
         CALL Info( Caller,"Radiation flux target on a boundary should be defined. " &
              // "Taking a meshed body as a default.", Level = 20)
       END IF
     END IF

     IF ( RadBody < 0 ) RadBody = 0
     IF ( RadBody>0 .AND. (RadBody /= Rbody .AND. RadBody /= Lbody) ) THEN
       CALL Error( Caller, 'Inconsistent direction information (Radiation Target Body)' )
       Message = 'Radiation Target: '//I2S(RadBody)//' Left, Right: '//&
                    I2S(Lbody)//I2S(Rbody)
       CALL Fatal( Caller, Message )
     END IF

     IF ( Lnode<=0 .OR. (RadBody>0 .AND. RadBody==Rbody) ) THEN
       Inode = Lnode
       Lnode = Rnode
       Rnode = Inode
     END IF
     Normal_in = -1.0
     IF ( RadBody <= 0 ) Normal_in = 1.0

     ! All the above complicated stuff has been done to provide a "Lnode" that is in the body where
     ! we are radiating to. However, if we have mesh deformation it is more reliable to use the
     ! body that is radiating as it probably more rigid. So here we do the swap. 
     IF(Lnode > 0 .AND. Rnode > 0 ) THEN
       Lnode = Rnode
       Normal_in = -Normal_in
     END IF
     
     ! The center of the boundary element
     r1(1) = SUM(ElementNodes % x)/k
     r1(2) = SUM(ElementNodes % y)/k
     r1(3) = SUM(ElementNodes % z)/k

     ! A presentative node of the radiation body
     r2(1) = Mesh % Nodes % x(Lnode)
     r2(2) = Mesh % Nodes % y(Lnode)
     r2(3) = Mesh % Nodes % z(Lnode)

     ! Direction to test the normal
     Nrm2 = r1 - r2
     dir = Normal_in*SUM(Nrm*Nrm2)
!------------------------------------------------------------------------------
   END FUNCTION DirectedNormalVector
!------------------------------------------------------------------------------
         
!------------------------------------------------------------------------------
   FUNCTION CheckShadowMeshGiven() RESULT(RT_Mesh)
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: RT_Mesh
!------------------------------------------------------------------------------
     INTEGER :: j,i0
     LOGICAL :: UseShadowMesh, BoundaryOnly
     CHARACTER(256) :: ShadowMeshName
!------------------------------------------------------------------------------
     RT_Mesh => Null()

     ShadowMeshName = GetString( Params, 'Viewfactor Surface Shadow Mesh',  UseShadowMesh )
     IF ( UseShadowMesh ) THEN
       ! Surface mesh from disk, the "bodyid" for radiation surfaces should match some
       ! (any) of the radiation boundary condition(!) numbers. Otherwise self contained.
       ! -------------------------------------------------------------------------------
       RT_Mesh => LoadShadowMesh(ShadowMeshName)
       j = 0
       DO i=1,RT_Mesh % NumberOfBulkElements
         IF ( RadiationBC(RT_Mesh % Elements(i) % BodyId) ) THEN
           j=j+1
           RT_Mesh % Elements(j) = RT_Mesh % Elements(i)
         END IF
       END DO
       RT_Mesh % NumberOfBulkElements = j
     ELSE
       ShadowMeshName = GetString( Params, 'Viewfactor Volume Shadow Mesh',  UseShadowMesh )
       ! This may be a coarser version of the primary mesh, if feasible. Only boundaries
       ! are loaded.
       ! -------------------------------------------------------------------------------
       IF ( UseShadowMesh ) THEN
         BoundaryOnly = .NOT. DoMapping          
         ! Load as non-distributed (1,0) so every MPI rank gets the complete
         ! shadow mesh independently — all ranks need it for ray testing.
         RT_Mesh => LoadMesh2( Model, "./", ShadowMeshName, BoundaryOnly, 1, 0 )

         ! For the shading mesh we will do rigid mesh mapping also if we use internal mapping         
         ! since for this coarse mesh the mapping has not been performed. 
         IF( DoMapping .OR. GetLogical( Model % Simulation,'Internal Rigid Mesh Mapping', Found ) ) THEN
           Model % Mesh => RT_Mesh 
           CALL RigidMeshMapping( Model, RT_Mesh, .TRUE. )
           Model % Mesh => Mesh 
         END IF         

         CALL SymmetryDuplication(RT_Mesh)
           
         i0 = RT_Mesh % NumberOfBulkElements
         j = 0
         DO i=1,RT_Mesh % NumberOfBoundaryElements
           IF ( RadiationBC(RT_Mesh % Elements(i0+i) % BoundaryInfo % Constraint) ) THEN
             j=j+1
             RT_Mesh % Elements(j) = RT_Mesh % Elements(i0+i)
           END IF
         END DO
         RT_Mesh % NumberOfBulkElements = j
       END IF
     END IF
!------------------------------------------------------------------------------
   END FUNCTION CheckShadowMeshGiven
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE FindInitialMinMax(Ni,N,Factors,RadiationOpen)
!------------------------------------------------------------------------------
     INTEGER :: ni,n
     LOGICAL :: RadiationOpen
     REAL(KIND=dp) :: Factors(:)
!------------------------------------------------------------------------------
     INTEGER :: nprob, i,j,k
     REAL(KIND=dp) :: s, Fmin, Fmax, Fave
     REAL(KIND=dp), ALLOCATABLE :: VF_cohorts(:)
!------------------------------------------------------------------------------
     nprob = 0
     Fave = 0.0_dp
     ALLOCATE(VF_cohorts(100))
     VF_cohorts = 0

     k = 0
     DO i=1,Ni
       s = 0.0_dp
       DO j=1,N
         IF(Factors((i-1)*N+j) < MinFactor) Factors((i-1)*N+j) = 0.0d0         
         s = s + Factors((i-1)*N+j)
       END DO
         
       IF( .NOT. RadiationOpen .AND. s < 0.5 ) nprob = nprob + 1

       IF(i == 1) THEN
         Fmin = s 
         Fmax = s
         k = 1
       ELSE
         IF(s < Fmin) THEN
           Fmin = s
           k = i
         END IF
         FMax = MAX( FMax,s )
       END IF

       j = CEILING(100*s)
       j = MIN(100,MAX(1,j))
       VF_cohorts(j) = VF_cohorts(j) + 1
         
       Fave = Fave + s         
     END DO
     Fave = Fave / Ni

     CALL Info( Caller, ' ', Level=3 )
     CALL info( Caller, 'Viewfactors before manipulation: ', Level=3 )
     WRITE( Message,'(A,ES14.6)') 'Minimum row sum: ',Fmin
     CALL Info( Caller, Message )
     WRITE( Message,'(A,ES14.6)') 'Maximum row sum: ',Fmax
     CALL Info( Caller, Message )
     WRITE( Message,'(A,ES14.6)') 'Average row sum: ',Fave
     CALL Info( Caller, Message )
     IF(nprob>0) CALL info( Caller, 'Number of rowsums below 0.5 is: '&
         //I2S(nprob)//' (out of '//I2S(n)//')')
     
     IF( InfoActive(10) ) THEN
       ! Report on the most problematic element which has too small viewfactors.
       IF( Fmin < 0.5 .AND. .NOT. RadiationOpen .OR. InfoActive(20) ) THEN
         Element => RadElements(k)
         j = Element % ElementIndex - Model % Mesh % NumberOfBulkElements
         CALL Info( Caller,'Location of minimum rowsum '//I2S(k)//' element '//I2S(j))
         PRINT *,'Indexes:',element % nodeindexes
         PRINT *,'X coord:',model % nodes % x(element % nodeindexes)
         PRINT *,'Y coord:',model % nodes % y(element % nodeindexes)
         PRINT *,'Z coord:',model % nodes % z(element % nodeindexes)
         DO i=1,3
           j = element % nodeindexes(i)
           PRINT *,'r:',sqrt(model % nodes % x(j)**2 + model % nodes % y(j)**2 + model % nodes % z(j)**2)
         END DO
       END IF
       
       DO i=1,100
         j = VF_cohorts(i)
         IF(j==0) CYCLE
         CALL Info(Caller,'VF percentile '//I2S(i-1)//'-'//I2S(i)//' count: '//I2S(j))
       END DO
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE FindInitialMinMax
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE FindNormalizedMinMax(Ni,N,Factors)
!------------------------------------------------------------------------------
     INTEGER :: Ni,N
     REAL(KIND=dP) :: Factors(:)
!------------------------------------------------------------------------------
     INTEGER :: i
     REAL(KIND=dp) :: s, Fmin,Fmax 
!------------------------------------------------------------------------------

     Fmin = HUGE(Fmin); Fmax = 0
     DO i=1,Ni
       s = 0.0_dp
       DO j=1,n
         s = s + Factors((i-1)*n+j)
       END DO
       Fmin = MIN(Fmin,s)
       Fmax = MAX(Fmax,s)
     END DO
       
     CALL Info( Caller, ' ', Level=3 )
     CALL info( Caller, 'Viewfactors/Radiator factors after manipulation: ')
     WRITE( Message,'(A,ES12.3)') 'Minimum row sum: ',Fmin
     CALL Info( Caller, Message )
     WRITE( Message,'(A,ES12.3)') 'Maximum row sum: ',Fmax
     CALL Info( Caller, Message )
     IF( Fmax > 1.001 ) THEN
       CALL Warn(Caller,'Rowsum of view factors/radiator factors should not be larger than one!')
     END IF
     IF( Fmin < 0.999_dp ) THEN
       ! For open BCs the view factor sum may be much less than one, otherwise not.
       IF(.NOT. ListCheckPresentAnyBC( Model,'Radiation Boundary Open') ) THEN          
         CALL Warn(Caller,'Rowsum of view factors should not be smaller than one!')
       END IF
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE FindNormalizedMinMax
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
FUNCTION ExtractSurfaces(Mesh,DoRadiators,RadElements,RadiationBC, &
      RadiationBody,RadiationOpen) RESULT(RadiationSurf)
!------------------------------------------------------------------------------
     TYPE(Mesh_t) :: Mesh
     LOGICAL :: RadiationOpen, DoRadiators
     TYPE(Element_t) :: RadElements(:)
     LOGICAL :: RadiationBC(:)
     INTEGER :: RadiationBody, RadiationSurf
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,t
     LOGICAL :: Found
     TYPE(ValueList_t), POINTER :: BC
     CHARACTER(:), ALLOCATABLE :: str 
!------------------------------------------------------------------------------
     RadiationSurf = 0
     RadiationOpen = .FALSE.
     DO i=1,Mesh % NumberOfBoundaryElements
       Element => GetBoundaryElement(i)
       IF ( GetElementFamily() == 1 ) CYCLE

       BC => GetBC()
       IF ( .NOT. ASSOCIATED( BC ) ) CYCLE
         
       IF ( DoRadiators ) THEN
         IF (GetLogical( BC, 'Radiator BC', Found )) THEN
           RadiationBC(GetBCId()) = .TRUE.
           RadiationOpen = RadiationOpen .OR. GetLogical(BC,'Radiation Boundary Open', Found)
           RadiationSurf = RadiationSurf + 1
           j = i + Mesh % NumberOFBulkElements
           RadElements(RadiationSurf) = Mesh % Elements(j)
         END IF
       ELSE
         str = GetString( BC, 'Radiation', Found )
         IF ( Found .AND. str == 'diffuse gray' ) THEN
           t = MAX(1, GetInteger( BC, 'Radiation Boundary', Found ))
           IF(t == RadiationBody) THEN
             RadiationBC(GetBCId()) = .TRUE.
             RadiationOpen = RadiationOpen .OR. GetLogical(BC,'Radiation Boundary Open',Found)
             RadiationSurf = RadiationSurf + 1
             j = i + Mesh % NumberOFBulkElements
             RadElements(RadiationSurf) = Mesh % Elements(j)
           END IF
         END IF
       END IF
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ExtractSurfaces
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE ExtractMeshInfo( Mesh, nActive, Coord, Surf, Type, Data, Perm, ElimBBox )
!------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: nActive

    REAL(KIND=dp), POINTER :: Coord(:)
    INTEGER, POINTER :: Type(:)
    INTEGER, ALLOCATABLE :: Surf(:)
    INTEGER, ALLOCATABLE, OPTIONAL :: Perm(:)
    REAL(KIND=dp), ALLOCATABLE, OPTIONAL :: DATA(:)
    LOGICAL, OPTIONAL :: ElimBBox
!------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n,i,j,k,l,cnt,m,mActive
    INTEGER, POINTER :: Ref(:)
    TYPE(ElementData_t), POINTER :: edata
    REAL(KIND=dp) :: minx, maxx, xeps
    LOGICAL :: SkipBBox
    LOGICAL, ALLOCATABLE :: NodeAtBBox(:,:)
    REAL(KIND=dp), POINTER :: pX(:), pY(:)
    REAL(KIND=dp), ALLOCATABLE :: pR(:)
!------------------------------------------------------------
    IF ( .NOT. ASSOCIATED(Mesh) ) RETURN
    
    n = Mesh % NumberOfBulkElements
    m = Mesh % NumberOfNodes
    nActive = n
    
    ! Compact renumbering of nodes.
    ALLOCATE( Ref(m) )
    Ref = 0
    
    ! We can save some extra time in checks when we eliminate the extremum elements
    ! that form the bounding box. They can never be shadowing elements. 
    SkipBBox = .FALSE.
    IF(PRESENT(ElimBBox)) SkipBBox = ElimBBox
    IF(SkipBBox) THEN
      ALLOCATE(NodeAtBBox(m,7))
      NodeAtBBox = .FALSE.

      ! Mark initial node list:
      DO i=1,n
        Element => Mesh % Elements(i)
        Ref(Element % NodeIndexes) = Ref(Element % NodeIndexes)+1
      END DO

      ! Assume cartesian BB
      DO i=1,3
        SELECT CASE(i)
        CASE ( 1 )
          pX => Mesh % Nodes % x
        CASE ( 2 )
          pX => Mesh % Nodes % y
        CASE( 3 ) 
          pX => Mesh % Nodes % z
        END SELECT
        minx = MINVAL(pX(1:m),Ref>0) 
        maxx = MAXVAL(pX(1:m),Ref>0) 
        xeps = EPSILON(xeps) +  1.0e-8 * ( maxx - minx )
        
        ! Nodes need to be on the same boundary in order them to be applicable
        ! for being at bounding box boundary.
        WHERE(ABS(pX(1:m) - minx) < xeps)
          NodeAtBBox(:,i) = .TRUE.
        END WHERE
        WHERE(ABS(pX(1:m) - maxx) < xeps)
          NodeAtBBox(:,3+i) = .TRUE.
        END WHERE
      END DO

#if 0
      ! Assume cylindrical BB
      ! We give it a shot without playing with keywords...      
      pX => Mesh % Nodes % x
      pY => Mesh % Nodes % y

      ALLOCATE(pR(m))      
      pR = SQRT(pX(1:m)**2+pY(1:m)**2)

      maxx = MAXVAL(pR(1:m),Ref>0)       
      xeps = EPSILON(xeps) +  1.0e-8 * maxx

      WHERE(ABS(pR - maxx) < xeps ) 
        NodeAtBBox(:,7) = .TRUE.
      END WHERE
      DEALLOCATE(pR)
#endif
      
      IF(InfoActive(10)) THEN
        DO i=1,7
          j = COUNT(NodeAtBBox(:,i))
          CALL Info(Caller,'Number of nodes at bounding box surface '//I2S(i)//': '//I2S(j))
        END DO
      END IF

          
      ! We skip the elements at bounding box boundaries from the active set. 
      nActive = 0
      DO i=1,n
        Element => Mesh % Elements(i)
        IF(SkipBBox .AND. Element % TYPE % ElementCode > 200 ) THEN
          k = 0
          DO j=1,6
            IF(ALL(NodeAtBBox(Element % NodeIndexes,j))) k=k+1
          END DO
          IF(k>0) CYCLE
        END IF
        nActive = nActive + 1
      END DO
      IF(nActive < n) THEN
        CALL Info(Caller,'Number of shading elements: '//I2S(nActive)//' (vs. '//I2S(n)//')')
      END IF

      ! Set this to zero, next time the bounding box is not used. 
      Ref = 0      
    END IF
    
    ALLOCATE( TYPE(nActive) )
    Type = 0
    cnt = 0
    l = 0
    DO i=1,n
      Element => Mesh % Elements(i)
      IF ( SkipBBox .AND. Element % TYPE % ElementCode > 200 ) THEN
        k = 0
        DO j=1,6
          IF(ALL(NodeAtBBox(Element % NodeIndexes,j))) k=k+1
        END DO
        IF(k>0) CYCLE
      END IF
      l = l+1
      Element % ElementIndex = l
      Type(l) = Element % Type % ElementCode
      IF ( Type(l) == 101 ) cnt=cnt+1
      Ref(Element % NodeIndexes) = Ref(Element % NodeIndexes)+1
    END DO
    mActive = COUNT(Ref>0)
    IF(mActive < m) THEN
      CALL Info(Caller,'Number of shading nodes: '//I2S(mActive)//' (vs. '//I2S(m)//')')
    END IF
!    Mesh % NumberOfBulkElements = l


    
    IF ( .NOT.ALLOCATED(Surf) ) ALLOCATE(Surf(4*nActive))
    ALLOCATE(Coord(3*mActive))
    
    IF ( PRESENT(data) .AND. cnt>0 ) THEN
      ! if planar circles found ....
      ALLOCATE(data(8*cnt), Perm(n))
      Data = 0
      j = 0
      DO i=1,n
        IF ( Type(i)==101 ) THEN
          j = j+1 
          Perm(i) = j
        END IF
      END DO
    END IF

    j = 0
    DO i=1,m
      IF ( Ref(i)>0 ) THEN        
        j = j + 1
        Ref(i) = j
        Coord(3*(j-1)+1) = Mesh % Nodes % x(i)
        Coord(3*(j-1)+2) = Mesh % Nodes % y(i)
        Coord(3*(j-1)+3) = Mesh % Nodes % z(i)

        Mesh % Nodes % x(j) = Mesh % Nodes % x(i)
        Mesh % Nodes % y(j) = Mesh % Nodes % y(i)
        Mesh % Nodes % z(j) = Mesh % Nodes % z(i)
      END IF
    END DO

    l = 0
    DO i=1,n
      Element => Mesh % Elements(i)
      IF(SkipBBox .AND. Element % TYPE % ElementCode > 200) THEN
        k = 0
        DO j=1,6
          IF(ALL(NodeAtBBox(Element % NodeIndexes,j))) k=k+1
        END DO
        IF(k>0) CYCLE
      END IF
      l = l+1
      Element % NodeIndexes = Ref(Element % NodeIndexes)

      SELECT CASE(TYPE(l)/100)
      CASE(1)
        ! Circle code
        IF ( PRESENT(Data) ) THEN
          edata => Element % PropertyData
          IF ( ASSOCIATED(edata) ) THEN
            Surf(4*(i-1)+1:4*(i-1)+4) = 0
            SELECT CASE(edata % Name)
            CASE DEFAULT ! circle
              ! inner radius (0), outer radius, center, normal = 8
              k = 8*(Perm(i)-1)
              Data(k+1:k+8) = edata % Values(1:8)
            END SELECT
          END IF
        END IF
      CASE(2)
        Surf(2*(l-1)+1:2*(l-1)+2) = Element % NodeIndexes-1
      CASE(3)
        Surf(4*(l-1)+1:4*(l-1)+3) = Element % NodeIndexes-1
      CASE(4)
        Surf(4*(l-1)+1:4*(l-1)+4) = Element % NodeIndexes-1
      CASE DEFAULT
        CALL Fatal(Caller,'Uknown Element!')
      END SELECT

      ! This messes with the mesh but I don't think we need it after.
      Mesh % Elements(l) = Mesh % Elements(i)

    END DO

    Mesh % NumberOfBulkElements = nActive
          
!------------------------------------------------------------------------------
   END SUBROUTINE ExtractMeshInfo
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE GetParentInfo(Parent, Element, body, node, rad )
!------------------------------------------------------------------------------
     INTEGER :: body, node
     LOGICAL :: rad
     TYPE(Element_t), POINTER :: Parent
     TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
     LOGICAL :: gotIt
     INTEGER i,j,k, matid
!------------------------------------------------------------------------------
     body =  0
     node = -1
     rad  = .FALSE.

     IF (.NOT. ASSOCIATED(Parent)) RETURN

     body = Parent % BodyId
     k = 0
     DO i=1,Parent % Type % NumberOfNodes
       gotIt =.TRUE.
       ! Find a node in parent element that is not a node in boundary element.
       DO j=1,Element % TYPE % NumberOfNodes
         IF ( Element % NodeIndexes(j) == Parent % NodeIndexes(i)) THEN
           k = k + 1
           gotIt=.FALSE.
           EXIT
         END IF
       END DO
       IF (gotIt) THEN
         node = Parent % NodeIndexes(i)
       END IF
     END DO

     IF(k /= Element % Type % NumberOfNodes ) THEN
       CALL Fatal(Caller,'Boundary element '//I2S(Element % ElementIndex - Mesh % NumberOfBulkElements)//&
             ' not included in parent '//I2S(parent % ElementIndex)//'!')               
     END IF

     IF( body > 0 ) THEN 
       rad  = ListGetLogical(Model % Bodies(body) % Values,'Radiative Body', GotIt )
       IF ( .NOT. GotIt ) THEN
         MatId = GetInteger( Model % Bodies(body) % Values,'Material', GotIt)
         IF( MatId == 0 ) THEN
            CALL Fatal(Caller,'Invalid material index in body, perhaps none')
         END IF 
         rad = ListGetLogical(Model % Materials(MatId) % Values,'Radiative Body', GotIt) 
         IF ( .NOT. GotIt ) THEN
           rad = ListCheckPresent(Model % Materials(MatId) % Values,'Emissivity') 
         END IF
       END IF
     ELSE
       CALL Warn(Caller,'body not associated.')
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE GetParentInfo
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE GetCartParameters( AreaEPS, FactEPS, RayEPS, Nray, LineInteg, QuadInteg, &
       TriInteg, UseSymmetry, ClosedForm, ShaftStat, Clip, RayCull )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: AreaEPS, FactEPS, RayEPS
     LOGICAL :: UseSymmetry
     INTEGER :: Nray, LineInteg, QuadInteg, TriInteg, ClosedForm, ShaftStat, Clip, RayCull
     LOGICAL :: ClosedFormGiven
!------------------------------------------------------------------------------
     LOGICAL :: GotIt
     CHARACTER(LEN=MAX_NAME_LEN) :: Str
!------------------------------------------------------------------------------
     AreaEPS = GetConstReal( Params, 'Viewfactor Area Tolerance',  GotIt )
     IF ( .NOT. GotIt ) AreaEPS = 1.0d-1

     FactEPS = GetConstReal( Params, 'Viewfactor Factor Tolerance ', GotIt )
     IF ( .NOT. GotIt ) FactEPS = 1.0d-2

     RayEPS = GetConstReal( Params, 'Viewfactor Raytrace Tolerace',  GotIt )
     IF ( .NOT. GotIt ) RayEPS = 1.0d-5

     Nrays = GetInteger( Params, 'Viewfactor Number of Rays ',  GotIt )
     IF ( .NOT. GotIt ) Nrays = 1

     LineInteg = GetInteger( Params, 'Viewfactor Line Integration Points ',  GotIt )
     IF ( .NOT. GotIt ) LineInteg = 2; ! ---> 1...12

     QuadInteg = GetInteger( Params, 'Viewfactor Quad Integration Points ',  GotIt )
     IF ( .NOT. GotIt ) QuadInteg = 4; ! ---> 1,2,4,9,16,... (square of 1...12)

     TriInteg = GetInteger( Params, 'Viewfactor Triangle Integration Points ',  GotIt )
     IF ( .NOT. GotIt ) TriInteg = 3;  ! ---> 1,3,6

     ! Evaluate the inner ("to area") integral in closed form rather than by
     ! quadrature.  Exact for planar patches, and much more accurate than the
     ! Gauss rule when the two patches are close together.  Warped quads and
     ! non-geometric normals fall back to quadrature automatically.
     ClosedForm = 0
     IF ( GetLogical( Params, 'Viewfactor Closed Form Integration', GotIt ) ) ClosedForm = 1
     ClosedFormGiven = GotIt
     IF ( GotIt .AND. ClosedForm == 1 ) &
       CALL Info(Caller,'Using closed form evaluation of the inner view factor integral')

     ! Report how many patch pairs have no possible blocker between them, i.e.
     ! how much visibility work a clipping based scheme would actually face.
     ! Diagnostic only: it costs a shaft cull per pair on top of the ray casting.
     ShaftStat = 0
     IF ( GetLogical( Params, 'Viewfactor Shaft Cull Statistics', GotIt ) ) ShaftStat = 1

     ! How the blocking of the view between two patches is resolved.  Ray
     ! casting samples it; clipping removes the shadows of the blockers from
     ! the target and integrates what is left, which is exact for the pair
     ! and gets the penumbra right instead of scaling by a hit fraction.
     ! Needs the closed form integration, and turns it on when asked for.
     ! Cull the rays of a patch pair to that pair's shaft candidates instead
     ! of walking the whole element tree once per ray.  The answer is bit for
     ! bit the same; it just replaces one tree traversal per ray by one per
     ! pair.  Whether that pays depends on the size of the shadow mesh as
     ! much as on the ray count, and the shadow mesh is only known once the
     ! element combining has run, so -1 asks InitStuff to decide.
     RayCull = -1
     IF ( ListCheckPresent( Params, 'Viewfactor Shaft Cull' ) ) THEN
       RayCull = 0
       IF ( GetLogical( Params, 'Viewfactor Shaft Cull', GotIt ) ) RayCull = 1
     END IF

     ClipInt = 0
     Str = ListGetString( Params, 'Viewfactor Shadowing', GotIt )
     IF ( GotIt ) THEN
       SELECT CASE( Str )
       CASE('ray casting')
         ClipInt = 0
       CASE('clipping')
         ClipInt = 1
         ! Clipping defaults to the closed form, which integrates a fragment
         ! natively.  It does not force it: asking for quadrature explicitly
         ! selects the fourth corner of the 2x2, where the same clipped
         ! fragments are fan triangulated and integrated by the triangle rule.
         ! That corner is a cross check -- it is slower and less accurate --
         ! but it is the only configuration that tells an error in the inner
         ! integral apart from an error in the shadowing.
         IF ( .NOT. ClosedFormGiven ) THEN
           ClosedForm = 1
         ELSE IF ( ClosedForm == 0 ) THEN
           CALL Info(Caller,'Clipping with numerical integration of the fragments')
         END IF
         CALL Info(Caller,'Resolving shadowing by clipping')
       CASE DEFAULT
         CALL Fatal(Caller,'Unknown "Viewfactor Shadowing": '//TRIM(Str))
       END SELECT
     END IF

     ! We may cheat with the symmetry by only counting the first half of the symmetric elements.
     UseSymmetry = .FALSE.
     IF(.NOT. DoRadiators) THEN
       UseSymmetry = GetLogical( Params, 'Viewfactor use symmetry',GotIt)       
       IF(UseSymmetry) THEN
         IF(NSymmetry == 0) THEN
           UseSymmetry = .FALSE.
         ELSE IF(NSymmetry == 1 ) THEN
           CALL Info(Caller,'"Viewfactor use symmetry" activated!')
         ELSE IF(NSymmetry > 1 ) THEN
           CALL Fatal(Caller,'"Viewfactor use symmetry" only applicable for one symmetry plane!')
         END IF
       END IF
     END IF
         
!------------------------------------------------------------------------------
   END SUBROUTINE GetCartParameters
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE WriteOutputFile(DoRadiators,Ni,n,Factors,RadiationBody)
!------------------------------------------------------------------------------
     LOGICAL :: DoRadiators
     REAL(KIND=dp) :: Factors(:)
     INTEGER :: Ni, n, RadiationBody
!------------------------------------------------------------------------------
     CHARACTER(:), ALLOCATABLE :: ViewFactorsFile, OutputName, TempString
     LOGICAL :: BinaryMode, SinglePrec
     REAL :: sval
     INTEGER :: i,j,k
     LOGICAL, ALLOCATABLE :: SaveMask(:)
!------------------------------------------------------------------------------

     IF ( DoRadiators ) THEN
       ViewFactorsFile = GetString( GetSimulation(),'Radiator Factors',GotIt)
       IF ( .NOT.GotIt ) ViewFactorsFile = 'RadiatorFactors.dat'
     ELSE
       ViewFactorsFile = GetString( GetSimulation(),'View Factors',GotIt)
       IF ( .NOT.GotIt ) ViewFactorsFile = 'ViewFactors.dat'
     END IF
     
     IF(RadiationBody > 1) THEN
       TempString = ViewFactorsFile
       ViewFactorsFile = TRIM(TempString)//I2S(RadiationBody)
     END IF
       
     IF (LEN_TRIM(Model % Mesh % Name) > 0) THEN
       OutputName = TRIM(OutputPath)//'/'//TRIM(Model % Mesh % Name)//'/'//TRIM(ViewFactorsFile)
     ELSE
       OutputName = TRIM(ViewFactorsFile)
     END IF
       
     ALLOCATE( SaveMask(SIZE(Factors)), STAT=istat)
     IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for SaveMask.')
         
     ! Use loser constraint for MinFactor as the errors can't be renormalized any more 
     MinFactor = MinFactor / 10.0_dp
         
     BinaryMode = ListGetLogical( Params,'Viewfactor Binary Output',Found ) 
     SinglePrec = getLogical( Params,'Viewfactor single precision',GotIt)
     
     SaveMask = ( Factors > MinFactor )
     
     IF( BinaryMode ) THEN
       CALL Info(Caller,'Saving view factors in binary mode',Level=5)

       OPEN( UNIT=VFUnit, FILE=TRIM(OutputName), FORM = 'unformatted', &
           ACCESS = 'stream', STATUS='replace', ACTION='write' )
       
       WRITE( VFUnit ) n

       DO i=1,Ni
         k = COUNT( SaveMask((i-1)*n+1:i*n) )
         WRITE( VFUnit ) k 
         DO j=1,n
           IF( SaveMask((i-1)*N+j ) ) THEN
             IF(SinglePrec) THEN
               sval = Factors((i-1)*N+j)
               WRITE( VFUnit ) j,sval
             ELSE
               WRITE( VFUnit ) j,Factors((i-1)*N+j)
             END IF
           END IF
         END DO
       END DO
     ELSE
       CALL Info(Caller,'Saving view factors in ascii mode',Level=5)
       
       OPEN( UNIT=VFUnit, FILE=TRIM(OutputName), STATUS='unknown', ACTION='write' )
       DO i=1,Ni
         k = COUNT( SaveMask((i-1)*n+1:i*n) )
         WRITE( VFUnit,* ) k
         DO j=1,n
           IF ( SaveMask((i-1)*N+j) ) THEN
             IF(SinglePrec) THEN
               sval = Factors((i-1)*N+j)
               WRITE( VFUnit,* ) i,j,sval
             ELSE
               WRITE( VFUnit,* ) i,j,Factors((i-1)*n+j)
             END IF
           END IF
         END DO
       END DO
     END IF
       
     CLOSE(VFUnit)
     DEALLOCATE( SaveMask ) 
!------------------------------------------------------------------------------
   END SUBROUTINE WriteOutputFile
!------------------------------------------------------------------------------

#include "huti_fdefs.h"

!> Local handle to the iterative methods for linear systems. 
!------------------------------------------------------------------------------
    SUBROUTINE IterSolv( N,x,b )
      IMPLICIT NONE

      INTEGER :: N
      REAL(KIND=dp), DIMENSION(n) :: x,b

      REAL(KIND=dp) :: dpar(50)
      INTEGER :: ipar(50),wsize
      REAL(KIND=dp), ALLOCATABLE :: work(:,:)
      INTEGER(KIND=addrInt) :: iterProc, mvProc, pcondProc, dProc
      INTEGER(KIND=AddrInt) :: AddrFunc
      EXTERNAL :: AddrFunc
!------------------------------------------------------------------------------
      HUTI_NDIM = N
      dProc = 0
    
      ipar = 0
      dpar = 0.0_dp

      HUTI_WRKDIM = HUTI_CG_WORKSIZE
      wsize = HUTI_WRKDIM
          
      HUTI_NDIM     = N
      HUTI_DBUGLVL  = 10
      HUTI_MAXIT    = 100
 
      ALLOCATE( work(N, wsize), STAT=istat )
      IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for IterSolv work.')

      
      work = 0.0_dp
      HUTI_TOLERANCE = 1.0d-6
      HUTI_MAXTOLERANCE = 1.0d20
      HUTI_INITIALX = HUTI_USERSUPPLIEDX
      HUTI_STOPC = HUTI_TRESID_SCALED_BYB
      
      iterProc  = AddrFunc(HUTI_D_CG)

      fm_G => Jacobian
      mvProc    = AddrFunc(fm_MatVec)
      fm_Diag => Jdiag
      pcondProc = AddrFunc(fm_DiagPrec)

      CALL IterCall( iterProc,x,b,ipar,dpar,work,mvProc,pcondProc, &
                dProc, dProc, dProc, dProc )
          
      DEALLOCATE( work )
    END SUBROUTINE IterSolv 


  END PROGRAM ViewFactors

!> \}
!> \}  


