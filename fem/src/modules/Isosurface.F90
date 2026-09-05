!/*****************************************************************************/
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
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   jpr@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 20.06.2007
! *
! *****************************************************************************/


SUBROUTINE IsosurfaceSolver_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  LOGICAL :: GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: Name
    
  IF( .NOT. ListCheckPresent( Solver % Values,'Variable') ) THEN
    Name = ListGetString( Solver % Values, 'Equation',GotIt)
    CALL ListAddString( Solver % Values,'Variable',&
        '-nooutput -global '//TRIM(Name)//'_var')
  END IF
  
END SUBROUTINE IsosurfaceSolver_init


!------------------------------------------------------------------------------
!> Subroutine for extracting isosurfaces in 3d and isolines in 2d.
!> The result will be a new mesh which will be added to the list of meshes.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE IsosurfaceSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------

  USE CoordinateSystems
  USE DefUtils
  USE SaveUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model    !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt       !< Timestep size for time dependent simulations
  LOGICAL :: Transient      !< Steady state or transient simulation
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: LevelVariableName
  TYPE(Mesh_t), POINTER :: Mesh, IsoMesh => NULL(), Pmesh, OrigMesh, &
      TempMesh => NULL()
  INTEGER :: i,j,k,l,n, dim,NoLevels,Level,NoNewElements
  REAL(KIND=dp) :: LevelValue

  INTEGER, TARGET :: TetraToTetraMap(1,4), PyramidToTetraMap(2,4), &
             WedgeToTetraMap(3,4), BrickToTetraMap(5,4), &
             TriangleToTriangleMap(1,3), QuadToTriangleMap(2,3)
  INTEGER, POINTER :: Map(:,:), Indexes(:)
  INTEGER :: NoOrigElements,calls=0,ierr

  TYPE(Variable_t), POINTER :: LevelVariable
  LOGICAL :: MovingMesh, Found, IsoCreated, Crossed

  ! The fields to be interpolated onto the isosurface, resolved once per call.
  INTEGER, PARAMETER :: MaxInterpolants = 32
  INTEGER :: ints
  CHARACTER(LEN=MAX_NAME_LEN) :: Vname(MaxInterpolants)
  TYPE(VariableTable_t) :: VfullTab(MaxInterpolants), VisoTab(MaxInterpolants)
  INTEGER, ALLOCATABLE :: eperm(:), InvPerm(:,:), ElemPerm(:), IsoEdge(:), &
      IsoLevel(:), EdgeOwner(:), EdgeNbrs(:,:), nEdgeNbrs(:), GeoLine(:,:)
  INTEGER :: nGeoNodes, nGeoLines
  REAL(KIND=dp), ALLOCATABLE :: GeoX(:), GeoY(:), GeoZ(:)
  INTEGER :: NewElemType, NewElemNodes, LevelDofs, NoEdges, IsAllocated, &
      NoIsoNodes, NoSurfaces, LevelComps, &
      NoMeshNodesOwned, NoEdgesOwned, NoIsoNodesOwned
  INTEGER :: ParSizes(6), ParTmp(6), SaveParTmp(6),ElemFirst,ElemLast
  REAL(KIND=dp), ALLOCATABLE :: ElemFun(:),Interpolant(:)
  REAL(KIND=dp), POINTER :: LevelFun(:), LevelValues(:,:)
  INTEGER, POINTER :: LevelPerm(:)
  TYPE(Element_t), POINTER  :: OrigElements(:), NewElements(:), Element
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: FixedSurface, SurfaceExist = .FALSE., GotIsoMask

  SAVE FixedSurface, SurfaceExist, Interpolant, InvPerm, Isomesh, SaveParTmp


  CALL Info( 'IsosurfaceSolver','-------------------------------------',Level=4 )
  CALL Info( 'IsosurfaceSolver','Determining the isosurface',Level=4 )
  CALL Info( 'IsosurfaceSolver','-------------------------------------',Level=4 )

  Mesh => GetMesh()
  OrigMesh => Mesh
  NoOrigElements = Mesh % NumberOfBulkElements
  OrigElements => Mesh % Elements

  ! Every exit from this solver has to fall through to the collective size
  ! reduction at label 100, so initialize the quantities reported there.
  !----------------------------------------------------------------------
  NoEdges = 0
  NoIsoNodes = 0
  NoSurfaces = 0
  NoMeshNodesOwned = 0
  NoEdgesOwned = 0
  NoIsoNodesOwned = 0
  NoNewElements = 0
  NewElements => NULL()
  IsoCreated = .FALSE.
  ParTmp = 0

  Params => GetSolverParams()
  FixedSurface = GetLogical( Params,'Isosurface Fixed',Found)
  MovingMesh = GetLogical( Params,'Moving Mesh',Found)

  !---------------------------------------------------------------
  ! Get the isosurface variable
  !---------------------------------------------------------------
  ! Deliberately not ThisOnly: if the field lives on another mesh this is
  ! where it gets interpolated onto ours, which is the point of the multimesh
  ! lookup. Safe here, as nothing has been done to the mesh yet.
  !---------------------------------------------------------------------------
  LevelVariableName = GetString( Params,'Isosurface Variable')
  LevelVariable => VariableGet( Mesh % Variables, LevelVariableName )

  IF (.NOT.ASSOCIATED(LevelVariable)) THEN
    CALL Error( 'Isosurface', 'Missing isosurface variable: ' // &
        TRIM(LevelVariableName) )
    GOTO 100
  END IF

  Levelfun => LevelVariable % Values
  LevelPerm => LevelVariable % Perm
  LevelDofs = LevelVariable % Dofs

  ! For a vector level function the isosurface is taken of its magnitude, over
  ! the spatial components. Note that this must follow the mesh, not the
  ! leading element dimension: masking the isosurface onto a boundary of a 3D
  ! mesh gives dim=2 while the field still has three components.
  !---------------------------------------------------------------------------
  LevelComps = MIN( Mesh % MeshDim, LevelDofs )
  
  !---------------------------------------------------------------
  ! Check the isosurface values
  !---------------------------------------------------------------
  NoLevels = 0
  LevelValues => ListGetConstRealArray( Params,'Isosurface values',Found)
  IF( Found ) THEN
    NoLevels = SIZE(LevelValues,1)
    LevelValue = LevelValues(1,1)
  ELSE
    LevelValue = ListGetCReal( Params,'Isosurface value',Found)
    IF( Found ) NoLevels = 1
  END IF

  IF(.NOT. Found ) THEN
    CALL Warn('IsosurfaceSolver','Could not determine Isosurface value')
    GOTO 100
  END IF

  !--------------------------------------------------------------------------
  ! If the mesh was created and will be fixed make a simplified interpolation
  !--------------------------------------------------------------------------
  IF( FixedSurface .AND. SurfaceExist ) THEN
    CALL Info('IsosurfaceSolver','Remapping the fields')

    ! ReMap only needs the parent nodes that each isosurface node sits
    ! between, so no temporal mesh is needed here.
    !--------------------------------------------------------------------
    CALL ReMap()

    ! Report the sizes from the call that actually created the surface.
    ParTmp = SaveParTmp
    GOTO 100
  END IF

  ! From here on the isosurface is going to be recreated, so drop the previous
  ! one right away. Releasing it only just before the new one is built would
  ! leave a stale mesh linked into Model % Meshes on every path that bails in
  ! between -- e.g. a partition whose surface has moved out of it, which would
  ! then go on being written out as if it were current.
  !----------------------------------------------------------------------
  CALL ReleaseIsoMesh()

  n = Mesh % MaxElementNodes
  ALLOCATE(ElemFun(n), ElemPerm(n))

  !---------------------------------------------------------------
  ! Map any element to tetrahedron or triangle
  !---------------------------------------------------------------
  TetraToTetraMap(1,:) = [1,2,3,4]

  PyramidToTetraMap(1,:) = [3,5,4,1]
  PyramidToTetraMap(2,:) = [3,5,2,1]

  WedgeToTetraMap(1,:) = [5,4,3,1]
  WedgeToTetraMap(2,:) = [5,3,2,1]
  WedgeToTetraMap(3,:) = [5,6,4,3]

  BrickToTetraMap(1,:) = [1,2,4,5]
  BrickToTetraMap(2,:) = [6,7,2,5]
  BrickToTetraMap(3,:) = [3,4,2,7]
  BrickToTetraMap(4,:) = [8,7,5,4]
  BrickToTetraMap(5,:) = [2,4,5,7]

  TriangleToTriangleMap(1,:) = [1,2,3]
  QuadToTriangleMap(1,:) = [1,2,3]
  QuadToTriangleMap(2,:) = [1,3,4]

  ElemFirst = 1
  ElemLast = Mesh % NumberOfBulkElements 
  GotIsoMask = .FALSE.

  IF( ListGetLogicalAnyBC( Model,'Create Isosurface' ) ) THEN
    GotIsoMask = .TRUE.
    ElemFirst = Mesh % NumberOfBulkElements + 1
    ElemLast = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
  ELSE IF( ListGetLogicalAnyBodyForce( Model,'Create Isosurface' ) ) THEN
    GotIsoMask = .TRUE.
  END IF

  !---------------------------------------------------------------
  ! Check the leading dimension
  !---------------------------------------------------------------
  dim = 0

  DO i = ElemFirst,ElemLast 
    Element => Mesh % Elements( i ) 
    Model % CurrentElement => Element
    j = GetElementFamily( Element )

    IF( GotIsoMask ) THEN
      IF( i <= Mesh % NumberOfBulkElements ) THEN
        IF( .NOT. GetLogical( GetBodyForce(Element),'Create Isosurface',Found) ) CYCLE
      ELSE 
        IF( .NOT. GetLogical( GetBC(Element),'Create Isosurface',Found) ) CYCLE
      END IF        
    END IF

    IF( j <= 2 ) THEN
      CYCLE
    ELSE IF( j <= 4 ) THEN
      dim = MAX( dim, 2 )
    ELSE
      dim = MAX( dim, 3 ) 
    END IF
  END DO

  ! The leading dimension has to agree in every partition. A partition holding
  ! no elements of the isosurface would otherwise take a different path from
  ! the rest, and everything from here on ends in a collective.
  !---------------------------------------------------------------------------
  IF( ParEnv % PEs > 1 ) THEN
    i = dim
    CALL MPI_ALLREDUCE(i,dim,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
  END IF

  CALL Info('IsosurfaceSolver','Leading element dimension is '//I2S(dim) )


  IF( dim == 2 ) THEN
    NewElemType = 303
    NewElemNodes = 3
  ELSE IF( dim == 3 ) THEN
    NewElemType = 504
    NewElemNodes = 4
  ELSE
    CALL Warn('IsoSurfaceSolver','Isosurface mesh can be created only for 2D or 3D')
    GOTO 100
  END IF


  !---------------------------------------------------------------
  ! Create a temporary mesh consisting only of tets (or triangles)
  ! in where the isosurface (or isoline) will be determined.
  ! The temporal mesh includes only potential master elements.
  !---------------------------------------------------------------
  ! The temporal mesh holds the elements the isosurface passes through, split
  ! into simplices. It does not depend on the level: the levels are applied
  ! later, when the nodes and the surfaces are made. So an element belongs
  ! here once, however many of the levels cross it, and the question to ask of
  ! each element is simply whether any of them does. Asking it the other way
  ! round, once per level, is what the never armed ElemSplit flag used to be
  ! guarding against.
  !---------------------------------------------------------------
  CALL Info('IsoSurfaceSolver','Creating a temporal mesh')

  DO IsAllocated = 0, 1
    
    j = 0
    
    DO i=ElemFirst,ElemLast
      
      Element => Mesh % Elements(i)
      Model % CurrentElement => Element

      IF( GotIsoMask ) THEN
        IF( i <= Mesh % NumberOfBulkElements ) THEN
          IF( .NOT. GetLogical( GetBodyForce(Element),'Create Isosurface',Found) ) CYCLE
        ELSE 
          IF( .NOT. GetLogical( GetBC(Element),'Create Isosurface',Found) ) CYCLE
        END IF
      END IF

      Indexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes

      IF( ASSOCIATED( LevelPerm ) ) THEN
        ElemPerm(1:n) = LevelPerm( Indexes )
        IF( ANY ( LevelPerm( Indexes) == 0 ) ) CYCLE
      ELSE
        ElemPerm(1:n) = Indexes(1:n)
      END IF

      SELECT CASE(GetElementFamily( Element ))
      CASE(3); Map => TriangleToTriangleMap
      CASE(4); Map => QuadToTriangleMap
      CASE(5); Map => TetraToTetraMap
      CASE(6); Map => PyramidToTetraMap
      CASE(7); Map => WedgeToTetraMap
      CASE(8); Map => BrickToTetraMap
      END SELECT

      IF( LevelDofs == 1 ) THEN
        ElemFun(1:n) =  LevelFun( ElemPerm(1:n) )
      ELSE
        ElemFun(1:n) = 0.0_dp
        DO l=1,LevelComps
          ElemFun(1:n) = ElemFun(1:n) + LevelFun( LevelDofs * ( ElemPerm(1:n) - 1) + l )**2
        END DO
        ElemFun(1:n) = SQRT( ElemFun(1:n) )
      END IF

      ! Everything above is independent of the level and is therefore done
      ! once per element rather than once per element per level.
      !------------------------------------------------------------------
      Crossed = .FALSE.
      DO Level = 1, NoLevels
        IF( NoLevels > 1 ) LevelValue = LevelValues(Level,1)
        IF( ALL ( ElemFun(1:n) - LevelValue < 0.0_dp ) ) CYCLE
        IF( ALL ( ElemFun(1:n) - LevelValue > 0.0_dp ) ) CYCLE
        Crossed = .TRUE.
        EXIT
      END DO
      IF( .NOT. Crossed ) CYCLE

      IF( IsAllocated == 0 ) THEN
        j = j + SIZE(Map,1)
        CYCLE
      END IF

      DO k=1,SIZE(Map,1)
        j = j + 1
        NewElements(j) = Element

        ! The intrinsic assignment above copies the pointer components of
        ! the parent element, i.e. the copies would share its index tables.
        ! FindMeshEdges reuses and overwrites those in place, so detach
        ! them all before the simplex is given its own identity.
        !------------------------------------------------------------------
        NewElements(j) % EdgeIndexes => NULL()
        NewElements(j) % FaceIndexes => NULL()
        NewElements(j) % BubbleIndexes => NULL()
        NewElements(j) % DGIndexes => NULL()
        NewElements(j) % BoundaryInfo => NULL()
        NewElements(j) % PDefs => NULL()
        NewElements(j) % PropertyData => NULL()

        NewElements(j) % TYPE => GetElementType(NewElemType)
        NewElements(j) % NDOFs = NewElemNodes
        NewElements(j) % BDOFs = 0
        NewElements(j) % DGDOFs = 0
        NewElements(j) % ElementIndex = j

        ALLOCATE(NewElements(j) % NodeIndexes(NewElemNodes))
        DO l=1,NewElemNodes
          NewElements(j) % NodeIndexes(l) = &
              Element % NodeIndexes(Map(k,l))
        END DO
      END DO
    END DO
    
    NoNewElements = j	
    
    IF( NoNewElements == 0 ) EXIT
    
    IF( IsAllocated == 0 ) THEN
      ALLOCATE(NewElements(NoNewElements))
    END IF    
  END DO

  ! Note that having nothing here is not an exit. The isosurface mesh is
  ! created either way, empty if need be, so that Model % Meshes holds the
  ! same meshes in every partition. The output stage walks that list and
  ! barriers once per mesh, so a partition short of one would hang the rest.
  !----------------------------------------------------------------------
  IF( NoNewElements == 0 ) THEN
    CALL Info('IsosurfaceSolver','No potential elements found in this partition',Level=7)
  ELSE
    CALL Info('IsosurfaceSolver','Found '//I2S(NoNewElements)//' potential elements') 
  END IF

  ! Wrap the simplices into a mesh of their own. Nothing below may touch the
  ! parent mesh any more: FindMeshEdges works in place and would otherwise
  ! overwrite its edge and face tables with ones describing this temporal
  ! mesh. The nodes are shared, as the simplices carry parent node indexes.
  !------------------------------------------------------------------------
  IF( NoNewElements > 0 ) THEN
    TempMesh => AllocateMesh()
    TempMesh % Name = 'isosurface temporal'
    DEALLOCATE( TempMesh % Nodes )
    TempMesh % Nodes => Mesh % Nodes
    TempMesh % NodesOrig => Mesh % Nodes
    TempMesh % NumberOfNodes = Mesh % NumberOfNodes
    TempMesh % MeshDim = dim
    TempMesh % Elements => NewElements
    TempMesh % NumberOfBulkElements = NoNewElements
    TempMesh % NumberOfBoundaryElements = 0
    TempMesh % MaxElementNodes = NewElemNodes
    TempMesh % MaxElementDOFs = NewElemNodes

    ! Find mesh edges in order to define the intersection points
    !-----------------------------------------------------------
    CALL Info('IsosurfaceSolver','Creating mesh edges',Level=9)
    IF( dim == 2 ) THEN
      CALL FindMeshEdges2D(TempMesh)
    ELSE
      CALL FindMeshEdges3D(TempMesh)
    END IF
    NoEdges = TempMesh % NumberOfEdges
  END IF

  ! Create a new mesh for isosurface
  !----------------------------------------------------------------
  CALL Info('IsosurfaceSolver','Creating a new mesh for isosurface')
  Isomesh => AllocateMesh()
  IsoMesh % Name = GetString( Params,'Mesh Name', Found )
  IF (.NOT. Found ) THEN
    IF( dim == 2 ) THEN
      IsoMesh % Name = "isoline"
    ELSE
      IsoMesh % Name = "isosurf"
    END IF
  END IF
  Isomesh % Changed = .TRUE.
  Isomesh % NumberOfBulkElements = 0
  Isomesh % NumberOfNodes = 0
  SurfaceExist = .TRUE.
  
  Isomesh % OutputActive = GetLogical( Params,'Isomesh Output Active',Found )
  IF(.NOT. Found ) Isomesh % OutputActive = .TRUE.

  ! TODO: the isosurface mesh carries no ParallelInfo.
  !
  ! Saving it in parallel works, because every partition now creates the mesh,
  ! empty where the surface does not reach, so the collectives in the output
  ! stage match up, and the vtu writer never consults ParallelInfo: it writes
  ! one piece per partition with local numbering. The nodes on partition
  ! interfaces are therefore duplicated, once per partition sharing the parent
  ! edge, which is what pvtu pieces look like anyway and is harmless to view.
  !
  ! ParallelInfo would be needed to make this a mesh in its own right, i.e. to
  ! compute a parallel norm on it, to run a solver on it, or to write it as one
  ! de-duplicated global mesh. It is derivable rather than guesswork, as every
  ! isosurface node sits on a known parent edge:
  !   GlobalDOFs    - a canonical key from the parent edge endpoints' global
  !                   ids, [MIN(g1,g2),MAX(g1,g2)], renumbered globally. Every
  !                   partition sharing the edge derives the same key, which is
  !                   exactly what makes the duplicates identifiable.
  !   GInterface    - true when both parent endpoints are interface nodes.
  !   NeighbourList - the intersection of the parent endpoints' neighbour lists.
  ! The same keys would let SaveGmshGeo2D and SaveSTLSurface be gathered into a
  ! single file, which is why both still refuse to run in parallel.
  !----------------------------------------------------------------------


  ! If requested number the output directories
  ! Alternatively one can number the output files
  !----------------------------------------------------------------
  IF( GetLogical(Params,'Mesh Numbering',Found) ) THEN
    Calls = Calls + 1
    Isomesh % name = TRIM(Isomesh % name) // i2s(calls)
  END IF
  CALL MakeDirectory( TRIM(Isomesh % name) // CHAR(0) )
  
  ! Create nodes and elements on the edge intersections  
  !----------------------------------------------------------------
  CALL Info('IsosurfaceSolver','Creating nodes on edge intersections',Level=9)
  NoIsoNodes = CreateNodes()    

  CALL Info('IsosurfaceSolver','Creating surfaces or lines on edge intersections',Level=9)
  NoSurfaces = CreateSurfaces()

  CALL IsoEdgeOwners()
  CALL CountOwned()
  CALL BuildIsoParallelInfo()
    
  ! Release the temporal mesh. Also the fixed surface case may do this now,
  ! as ReMap works on the parent node indexes stored in InvPerm.
  !----------------------------------------------------------------------	
  CALL ReleaseTempMesh()

  IF( ALLOCATED( ElemFun ) ) DEALLOCATE( ElemFun, ElemPerm )  

  ! Add the new mesh into the list 
  !----------------------------------------------------------------------
  PMesh => Model % Meshes
  DO WHILE( ASSOCIATED(PMesh % Next) )
    PMesh => PMesh % Next
  END DO
  PMesh % Next => Isomesh 

  IsoCreated = .TRUE.

  ! Information of the new system size, also in parallel.
  ! NOTE: the reduction below is collective. Every exit above therefore
  ! branches here instead of returning, as otherwise a partition that the
  ! isosurface does not cross would leave the others hanging in it.
  !----------------------------------------------------------------------
100 CONTINUE

  IF( IsoCreated ) THEN
    ParTmp(1) = NoMeshNodesOwned
    ParTmp(2) = NoOrigElements
    ParTmp(3) = NoNewElements
    ParTmp(4) = NoEdgesOwned
    ParTmp(5) = NoIsoNodesOwned
    ParTmp(6) = NoSurfaces
    SaveParTmp = ParTmp
  END IF

  IF( ParEnv % PEs > 1 ) THEN
    CALL MPI_ALLREDUCE(ParTmp,ParSizes,6,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
  ELSE
    ParSizes = ParTmp
  END IF

  ! Only the partitions that really created a mesh may report on it or
  ! save it, the others just participated in the reduction above.
  !----------------------------------------------------------------------
  IF( IsoCreated ) THEN
    CALL Info('IsoSurfaceSolver','Information on isosurface mesh sizes')

    WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(1),' nodes'
    CALL Info('IsoSurfaceSolver',Message)
    WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(2),' elements'
    CALL Info('IsoSurfaceSolver',Message)

    WRITE ( Message,'(A,I0,A)') 'Temporal mesh has ',ParSizes(3),' elements'
    CALL Info('IsoSurfaceSolver',Message)
    WRITE ( Message,'(A,I0,A)') 'Temporal mesh has ',ParSizes(4),' edges'
    CALL Info('IsoSurfaceSolver',Message)

    WRITE ( Message,'(A,I0,A)') 'Isosurface mesh has ',ParSizes(5),' nodes'
    CALL Info('IsoSurfaceSolver',Message)
    WRITE ( Message,'(A,I0,A)') 'Isosurface mesh has ',ParSizes(6),' elements'
    CALL Info('IsoSurfaceSolver',Message)

    IF( GetLogical( Params,'Save Gmsh Geo File',Found ) ) THEN
      CALL SaveGmshGeo2D(IsoMesh)
    END IF
    IF( GetLogical( Params,'Save STL File',Found ) ) THEN
      ! The triangles were pointed up the gradient when they were made, so
      ! the saver must not turn the patches to face away from a centre.
      CALL SaveSTLSurface(IsoMesh, Params, KeepOrientation = .TRUE. )
    END IF
  END IF

  ! Just create some variable that will act as a norm.
  IF(SIZE(Solver % Variable % Values) == 1 ) THEN
    Solver % Variable % Values = SUM(ParSizes)
    Solver % Variable % Norm = SUM(ParSizes)
  END IF
  

CONTAINS

  !----------------------------------------------------------------
  !> Work out which partition owns each edge of the temporal mesh.
  !>
  !> Holding both parent nodes of an edge is not the same as holding the
  !> edge: a tetrahedron can have two nodes on an interface and join them
  !> through its own interior, in which case the neighbour has the nodes and
  !> not the edge. Deciding ownership from the nodes alone therefore gives
  !> edges away to partitions that do not have them, and they end up counted
  !> by nobody. So ask.
  !>
  !> Every partition sends each of its neighbours the edges it holds that the
  !> neighbour could conceivably hold too, named by the global numbers of the
  !> two parent nodes so that both sides name them the same way. An edge that
  !> comes back is one we both have, and the lowest numbered of us owns it.
  !----------------------------------------------------------------
  SUBROUTINE IsoEdgeOwners()

    TYPE HashEntry_t
      INTEGER :: g1, g2, ind
      TYPE(HashEntry_t), POINTER :: Next => NULL()
    END TYPE HashEntry_t
    TYPE HashHead_t
      TYPE(HashEntry_t), POINTER :: Head => NULL()
    END TYPE HashHead_t

    TYPE(HashHead_t), ALLOCATABLE :: Hash(:)
    TYPE(HashEntry_t), POINTER :: Hp, Hp1
    INTEGER, ALLOCATABLE :: Key(:,:), Cand(:,:), nCand(:), &
        nSend(:), nRecv(:), SendDispl(:), RecvDispl(:), SendBuf(:), RecvBuf(:)
    INTEGER, POINTER :: a(:), b(:)
    INTEGER :: i,j,k,p,n1,n2,g1,g2,nH,h,pe,ns,nr,MaxCand,nE

    IF( ALLOCATED( EdgeOwner ) ) DEALLOCATE( EdgeOwner, EdgeNbrs, nEdgeNbrs )
    ALLOCATE( EdgeOwner(MAX(NoEdges,1)), &
        EdgeNbrs(MAX(ParEnv % PEs,1),MAX(NoEdges,1)), nEdgeNbrs(MAX(NoEdges,1)) )
    EdgeOwner = ParEnv % MyPe
    EdgeNbrs = 0
    nEdgeNbrs = 1
    DO i=1,NoEdges
      EdgeNbrs(1,i) = ParEnv % MyPe
    END DO

    ! The only exit that every partition takes together. Everything below
    ! ends in an all to all, so a partition the isosurface misses has to walk
    ! through it as well, simply with nothing to say.
    !------------------------------------------------------------------
    IF( ParEnv % PEs <= 1 ) RETURN

    nE = NoEdges
    IF( .NOT. ASSOCIATED( TempMesh ) ) nE = 0
    IF( .NOT. ASSOCIATED( Mesh % ParallelInfo % NeighbourList ) ) nE = 0

    pe = ParEnv % PEs
    MaxCand = pe

    ALLOCATE( Key(2,MAX(nE,1)), Cand(MaxCand,MAX(nE,1)), nCand(MAX(nE,1)) )
    nCand = 0

    DO i=1,nE
      n1 = TempMesh % Edges(i) % NodeIndexes(1)
      n2 = TempMesh % Edges(i) % NodeIndexes(2)
      g1 = Mesh % ParallelInfo % GlobalDOFs(n1)
      g2 = Mesh % ParallelInfo % GlobalDOFs(n2)
      Key(1,i) = MIN(g1,g2)
      Key(2,i) = MAX(g1,g2)

      ! Who could possibly have this edge as well
      a => Mesh % ParallelInfo % NeighbourList(n1) % Neighbours
      b => Mesh % ParallelInfo % NeighbourList(n2) % Neighbours
      DO j=1,SIZE(a)
        IF( a(j) == ParEnv % MyPe ) CYCLE
        DO k=1,SIZE(b)
          IF( a(j) == b(k) ) THEN
            nCand(i) = nCand(i) + 1
            Cand(nCand(i),i) = a(j)
            EXIT
          END IF
        END DO
      END DO
    END DO

    ! Hash my own edges by their key, to recognise the ones sent back
    nH = 2*nE + 1
    ALLOCATE( Hash(nH) )
    DO i=1,nE
      h = MODULO( INT(Key(1,i),8)*73856093_8 + INT(Key(2,i),8)*19349663_8, INT(nH,8) ) + 1
      ALLOCATE( Hp )
      Hp % g1 = Key(1,i); Hp % g2 = Key(2,i); Hp % ind = i
      Hp % Next => Hash(h) % Head
      Hash(h) % Head => Hp
    END DO

    ALLOCATE( nSend(pe), nRecv(pe), SendDispl(pe), RecvDispl(pe) )
    nSend = 0
    DO i=1,nE
      DO j=1,nCand(i)
        nSend(Cand(j,i)+1) = nSend(Cand(j,i)+1) + 2
      END DO
    END DO

    CALL MPI_ALLTOALL( nSend, 1, MPI_INTEGER, nRecv, 1, MPI_INTEGER, &
        ELMER_COMM_WORLD, ierr )

    SendDispl(1) = 0; RecvDispl(1) = 0
    DO p=2,pe
      SendDispl(p) = SendDispl(p-1) + nSend(p-1)
      RecvDispl(p) = RecvDispl(p-1) + nRecv(p-1)
    END DO
    ns = SUM(nSend); nr = SUM(nRecv)
    ALLOCATE( SendBuf(MAX(ns,1)), RecvBuf(MAX(nr,1)) )

    nSend = 0
    DO i=1,nE
      DO j=1,nCand(i)
        p = Cand(j,i) + 1
        SendBuf(SendDispl(p) + nSend(p) + 1) = Key(1,i)
        SendBuf(SendDispl(p) + nSend(p) + 2) = Key(2,i)
        nSend(p) = nSend(p) + 2
      END DO
    END DO

    CALL MPI_ALLTOALLV( SendBuf, nSend, SendDispl, MPI_INTEGER, &
        RecvBuf, nRecv, RecvDispl, MPI_INTEGER, ELMER_COMM_WORLD, ierr )

    ! An edge that came back from p is one we both hold
    DO p=1,pe
      DO k=1,nRecv(p),2
        g1 = RecvBuf(RecvDispl(p)+k)
        g2 = RecvBuf(RecvDispl(p)+k+1)
        h = MODULO( INT(g1,8)*73856093_8 + INT(g2,8)*19349663_8, INT(nH,8) ) + 1
        Hp => Hash(h) % Head
        DO WHILE( ASSOCIATED(Hp) )
          IF( Hp % g1 == g1 .AND. Hp % g2 == g2 ) THEN
            i = Hp % ind
            EdgeOwner(i) = MIN( EdgeOwner(i), p-1 )
            nEdgeNbrs(i) = nEdgeNbrs(i) + 1
            EdgeNbrs(nEdgeNbrs(i),i) = p-1
            EXIT
          END IF
          Hp => Hp % Next
        END DO
      END DO
    END DO

    ! The neighbour lists are expected sorted, the owner first.
    DO i=1,nE
      DO j=2,nEdgeNbrs(i)
        k = EdgeNbrs(j,i)
        p = j-1
        DO WHILE( p >= 1 )
          IF( EdgeNbrs(p,i) <= k ) EXIT
          EdgeNbrs(p+1,i) = EdgeNbrs(p,i)
          p = p-1
        END DO
        EdgeNbrs(p+1,i) = k
      END DO
    END DO

    DO i=1,nH
      Hp => Hash(i) % Head
      DO WHILE( ASSOCIATED(Hp) )
        Hp1 => Hp % Next; DEALLOCATE(Hp); Hp => Hp1
      END DO
    END DO
    DEALLOCATE( Hash, Key, Cand, nCand, nSend, nRecv, SendDispl, RecvDispl, &
        SendBuf, RecvBuf )

  END SUBROUTINE IsoEdgeOwners


  !----------------------------------------------------------------
  !> Give the isosurface mesh its ParallelInfo.
  !>
  !> Each isosurface node sits on an edge of the temporal mesh, and every
  !> partition holding that edge makes the same nodes on it, so the nodes
  !> inherit the edge's owner and its list of sharers. What is left is the
  !> global numbering: each partition numbers the nodes it owns, starting
  !> from the running total of the partitions before it, and then tells the
  !> others the numbers they could not know. A node is named to them by the
  !> global numbers of the two parent nodes it lies between and by which of
  !> the levels it belongs to, which both sides derive the same way.
  !----------------------------------------------------------------
  SUBROUTINE BuildIsoParallelInfo()

    INTEGER, ALLOCATABLE :: nSend(:), nRecv(:), SendDispl(:), RecvDispl(:), &
        SendBuf(:), RecvBuf(:)
    INTEGER :: i,j,k,p,e,pe,nOwn,Offset,ns,nr,g1,g2,MaxId,nIso
    TYPE(ParallelInfo_t), POINTER :: PI

    ! Again the only exit that all partitions take together: there are
    ! collectives below and a partition the isosurface misses reaches them too.
    !------------------------------------------------------------------
    IF( ParEnv % PEs <= 1 ) RETURN

    nIso = NoIsoNodes
    IF( .NOT. ASSOCIATED( Isomesh ) ) nIso = 0
    IF( .NOT. ALLOCATED( IsoEdge ) ) nIso = 0
    IF( .NOT. ALLOCATED( EdgeOwner ) ) nIso = 0
    IF( .NOT. ASSOCIATED( Mesh % ParallelInfo % NeighbourList ) ) nIso = 0

    pe = ParEnv % PEs
    PI => Isomesh % ParallelInfo

    ALLOCATE( PI % GlobalDOFs(MAX(NoIsoNodes,1)), &
        PI % GInterface(MAX(NoIsoNodes,1)), &
        PI % NeighbourList(MAX(NoIsoNodes,1)) )
    PI % GlobalDOFs = 0
    PI % GInterface = .FALSE.
    DO i=1,MAX(NoIsoNodes,1)
      NULLIFY( PI % NeighbourList(i) % Neighbours )
    END DO

    nOwn = 0
    DO i=1,nIso
      e = IsoEdge(i)
      IF( e <= 0 ) CYCLE
      k = nEdgeNbrs(e)
      ALLOCATE( PI % NeighbourList(i) % Neighbours(k) )
      PI % NeighbourList(i) % Neighbours = EdgeNbrs(1:k,e)
      PI % GInterface(i) = ( k > 1 )
      IF( EdgeOwner(e) == ParEnv % MyPe ) nOwn = nOwn + 1
    END DO

    ! Where this partition's block of numbers starts
    CALL MPI_SCAN( nOwn, Offset, 1, MPI_INTEGER, MPI_SUM, ELMER_COMM_WORLD, ierr )
    Offset = Offset - nOwn

    k = 0
    DO i=1,nIso
      e = IsoEdge(i)
      IF( e <= 0 ) CYCLE
      IF( EdgeOwner(e) /= ParEnv % MyPe ) CYCLE
      k = k + 1
      PI % GlobalDOFs(i) = Offset + k
    END DO

    ! Tell the sharers the numbers of the nodes they do not own
    ALLOCATE( nSend(pe), nRecv(pe), SendDispl(pe), RecvDispl(pe) )
    nSend = 0
    DO i=1,nIso
      e = IsoEdge(i)
      IF( e <= 0 ) CYCLE
      IF( EdgeOwner(e) /= ParEnv % MyPe ) CYCLE
      DO j=1,nEdgeNbrs(e)
        IF( EdgeNbrs(j,e) == ParEnv % MyPe ) CYCLE
        nSend(EdgeNbrs(j,e)+1) = nSend(EdgeNbrs(j,e)+1) + 4
      END DO
    END DO

    CALL MPI_ALLTOALL( nSend, 1, MPI_INTEGER, nRecv, 1, MPI_INTEGER, &
        ELMER_COMM_WORLD, ierr )

    SendDispl(1) = 0; RecvDispl(1) = 0
    DO p=2,pe
      SendDispl(p) = SendDispl(p-1) + nSend(p-1)
      RecvDispl(p) = RecvDispl(p-1) + nRecv(p-1)
    END DO
    ns = SUM(nSend); nr = SUM(nRecv)
    ALLOCATE( SendBuf(MAX(ns,1)), RecvBuf(MAX(nr,1)) )

    nSend = 0
    DO i=1,nIso
      e = IsoEdge(i)
      IF( e <= 0 ) CYCLE
      IF( EdgeOwner(e) /= ParEnv % MyPe ) CYCLE
      g1 = Mesh % ParallelInfo % GlobalDOFs(InvPerm(1,i))
      g2 = Mesh % ParallelInfo % GlobalDOFs(InvPerm(2,i))
      DO j=1,nEdgeNbrs(e)
        IF( EdgeNbrs(j,e) == ParEnv % MyPe ) CYCLE
        p = EdgeNbrs(j,e) + 1
        SendBuf(SendDispl(p)+nSend(p)+1) = MIN(g1,g2)
        SendBuf(SendDispl(p)+nSend(p)+2) = MAX(g1,g2)
        SendBuf(SendDispl(p)+nSend(p)+3) = IsoLevel(i)
        SendBuf(SendDispl(p)+nSend(p)+4) = PI % GlobalDOFs(i)
        nSend(p) = nSend(p) + 4
      END DO
    END DO

    CALL MPI_ALLTOALLV( SendBuf, nSend, SendDispl, MPI_INTEGER, &
        RecvBuf, nRecv, RecvDispl, MPI_INTEGER, ELMER_COMM_WORLD, ierr )

    DO k=1,nr,4
      DO i=1,nIso
        IF( PI % GlobalDOFs(i) /= 0 ) CYCLE
        IF( IsoLevel(i) /= RecvBuf(k+2) ) CYCLE
        g1 = Mesh % ParallelInfo % GlobalDOFs(InvPerm(1,i))
        g2 = Mesh % ParallelInfo % GlobalDOFs(InvPerm(2,i))
        IF( MIN(g1,g2) /= RecvBuf(k) .OR. MAX(g1,g2) /= RecvBuf(k+1) ) CYCLE
        PI % GlobalDOFs(i) = RecvBuf(k+3)
        EXIT
      END DO
    END DO

    DEALLOCATE( nSend, nRecv, SendDispl, RecvDispl, SendBuf, RecvBuf )

    ! Every node should have a number by now, and they should run to the
    ! total number of distinct nodes.
    j = 0
    DO i=1,nIso
      IF( PI % GlobalDOFs(i) == 0 ) j = j + 1
    END DO
    IF( j > 0 ) THEN
      CALL Warn('IsosurfaceSolver','Left '//I2S(j)//' isosurface nodes unnumbered!')
    END IF
    MaxId = 0
    IF( nIso > 0 ) MaxId = MAXVAL( PI % GlobalDOFs(1:nIso) )
    j = MaxId
    CALL MPI_ALLREDUCE( j, MaxId, 1, MPI_INTEGER, MPI_MAX, ELMER_COMM_WORLD, ierr )
    CALL Info('IsosurfaceSolver','Isosurface mesh global nodes: '//I2S(MaxId),Level=7)

  END SUBROUTINE BuildIsoParallelInfo


  !----------------------------------------------------------------
  !> Count the entities this partition owns, so that the sizes reported for
  !> the run count each of them once rather than once per partition holding
  !> it. The elements need no such care, being owned outright.
  !----------------------------------------------------------------
  SUBROUTINE CountOwned()

    INTEGER :: i

    NoMeshNodesOwned = Mesh % NumberOfNodes
    NoEdgesOwned = NoEdges
    NoIsoNodesOwned = NoIsoNodes
    IF( ParEnv % PEs <= 1 ) RETURN

    IF( ASSOCIATED( Mesh % ParallelInfo % NeighbourList ) ) THEN
      NoMeshNodesOwned = 0
      DO i=1,Mesh % NumberOfNodes
        IF( Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPe ) &
            NoMeshNodesOwned = NoMeshNodesOwned + 1
      END DO
    END IF

    IF( ALLOCATED( EdgeOwner ) ) THEN
      NoEdgesOwned = 0
      DO i=1,NoEdges
        IF( EdgeOwner(i) == ParEnv % MyPe ) NoEdgesOwned = NoEdgesOwned + 1
      END DO

      ! An isosurface node belongs to whoever owns the edge it sits on: the
      ! partitions sharing that edge all make the same nodes on it.
      IF( ALLOCATED( IsoEdge ) ) THEN
        NoIsoNodesOwned = 0
        DO i=1,NoIsoNodes
          IF( IsoEdge(i) <= 0 ) CYCLE
          IF( EdgeOwner(IsoEdge(i)) == ParEnv % MyPe ) &
              NoIsoNodesOwned = NoIsoNodesOwned + 1
        END DO
      END IF
    END IF

  END SUBROUTINE CountOwned


  !----------------------------------------------------------------
  !> Unlink a previously created isosurface mesh and release it.
  !----------------------------------------------------------------
  SUBROUTINE ReleaseIsoMesh()

    TYPE(Mesh_t), POINTER :: Pmesh

    IF( .NOT. SurfaceExist ) RETURN

    CALL Info('IsosurfaceSolver','Releasing previous mesh',Level=9)
    Pmesh => Model % Meshes
    DO WHILE(ASSOCIATED(Pmesh % Next))
      IF ( ASSOCIATED(Pmesh % next, Isomesh)) THEN
        Pmesh % next => Isomesh % Next
        EXIT
      END IF
      Pmesh => Pmesh % next
    END DO
    CALL ReleaseMesh(Isomesh)
    Isomesh => NULL()
    SurfaceExist = .FALSE.

  END SUBROUTINE ReleaseIsoMesh


  !----------------------------------------------------------------
  !> Collect the fields that are to be interpolated onto the isosurface.
  !>
  !> The lookup on the parent mesh is deliberately NOT ThisOnly. If the field
  !> lives on some other mesh, this is where it gets created here and
  !> interpolated from the mesh that holds the latest update, which is the
  !> whole point of the multimesh lookup. That is also why it may only happen
  !> once per call: see the warning at the interpolating branch of
  !> VariableGet in Lists.F90. Earlier revisions called it per isosurface
  !> node from inside the node loops, where it would silently interpolate
  !> over the very field being read.
  !----------------------------------------------------------------
  SUBROUTINE GetInterpolants()

    INTEGER :: i
    LOGICAL :: Found

    ints = 0
    DO i=1,MaxInterpolants
      Vname(i) = GetString( Params, &
          ComponentName('Isosurface interpolant',i), Found )
      IF(.NOT. Found) Vname(i) = GetString( Params, &
          ComponentName('Variable',i), Found )
      IF( .NOT. Found ) EXIT

      VfullTab(i) % Variable => VariableGet( Mesh % Variables, Vname(i) )
      IF(.NOT. ASSOCIATED( VfullTab(i) % Variable ) ) EXIT

      ints = i
    END DO

    IF( ints == MaxInterpolants ) THEN
      CALL Warn('IsosurfaceSolver','Only the first '//I2S(MaxInterpolants)//&
          ' interpolants are honoured!')
    END IF

  END SUBROUTINE GetInterpolants


  !----------------------------------------------------------------
  !> Bind the isosurface side of the interpolants. These we created
  !> ourselves, so the lookup is ThisOnly: there is nothing to
  !> interpolate from, and falling through would be circular.
  !----------------------------------------------------------------
  SUBROUTINE GetIsoInterpolants()

    INTEGER :: i

    DO i=1,ints
      VisoTab(i) % Variable => VariableGet( IsoMesh % Variables, Vname(i), &
          ThisOnly = .TRUE. )
    END DO

  END SUBROUTINE GetIsoInterpolants


  !----------------------------------------------------------------
  !> Release the temporal simplex mesh used to locate the isosurface.
  !> The nodes are shared with the parent mesh and must be left alone,
  !> hence this cannot be done with ReleaseMesh.
  !----------------------------------------------------------------
  SUBROUTINE ReleaseTempMesh()

    INTEGER :: i

    IF( .NOT. ASSOCIATED( TempMesh ) ) RETURN

    ! This releases the edge table and the EdgeIndexes of the simplices.
    ! It must be done while the elements are still in place.
    !----------------------------------------------------------------
    CALL ReleaseMeshEdgeTables( TempMesh )

    IF( ASSOCIATED( TempMesh % Elements ) ) THEN
      DO i=1,TempMesh % NumberOfBulkElements
        IF( ASSOCIATED( TempMesh % Elements(i) % NodeIndexes ) ) &
            DEALLOCATE( TempMesh % Elements(i) % NodeIndexes )
        TempMesh % Elements(i) % NodeIndexes => NULL()
      END DO
      DEALLOCATE( TempMesh % Elements )
    END IF
    TempMesh % Elements => NULL()
    NewElements => NULL()

    ! Owned by the parent mesh!
    TempMesh % Nodes => NULL()
    TempMesh % NodesOrig => NULL()

    DEALLOCATE( TempMesh )
    TempMesh => NULL()

  END SUBROUTINE ReleaseTempMesh


  !----------------------------------------------------------------
  !> Create nodes for the isosurface / isoline.
  !----------------------------------------------------------------
  FUNCTION CreateNodes() RESULT ( NoIsoNodes )

    INTEGER :: NoIsoNodes
    TYPE(Element_t), POINTER :: Edge
    INTEGER :: i,j,k,l,n1,n2,m1,m2
    REAL(KIND=dp) :: t,t1,t2,x1,x2,y1,y2,z1,z2

    TYPE(Variable_t), POINTER :: Vfull, Viso
    INTEGER, POINTER :: Vperm(:)
    REAL(KIND=dp), POINTER :: Vals(:)


    ! Make a table of the field variables and allocate space
    !----------------------------------------------------------------
    NoIsoNodes = 0
    CALL GetInterpolants()

    ! The isosurface may well miss this partition altogether. Give it an empty
    ! mesh anyway, complete with the same set of variables the other partitions
    ! will have, so that Model % Meshes agrees everywhere. The output stage
    ! barriers once per mesh, and the field lists have to match too.
    !--------------------------------------------------------------------------
    NoEdges = 0
    IF( ASSOCIATED( TempMesh ) ) NoEdges = TempMesh % NumberOfEdges

    IF( NoEdges == 0 ) THEN
      CALL Info('IsoSurfaceSolver','Creating an empty isosurface mesh',Level=8)

      ! Nodes was already allocated by AllocateMesh, only the coordinates are
      ! missing. Zero sized rather than unassociated, so that anything walking
      ! the mesh sees the same structure as in a partition that has a surface.
      !------------------------------------------------------------------------
      ALLOCATE( IsoMesh % Nodes % x(0), IsoMesh % Nodes % y(0), IsoMesh % Nodes % z(0) )
      Isomesh % NumberOfNodes = 0
      Isomesh % Nodes % NumberOfNodes = 0
      Isomesh % MeshDim = dim

      ALLOCATE( Isomesh % Elements(0) )
      Isomesh % NumberOfBulkElements = 0
      Isomesh % NumberOfBoundaryElements = 0
      Isomesh % NumberOfEdges = 0
      Isomesh % NumberOfFaces = 0

      CALL VariableAdd( IsoMesh % Variables, IsoMesh, Solver, &
          'Coordinate 1',1,IsoMesh % Nodes % x )
      CALL VariableAdd( IsoMesh % Variables, IsoMesh, Solver, &
          'Coordinate 2',1,IsoMesh % Nodes % y )
      CALL VariableAdd( IsoMesh % Variables, IsoMesh, Solver, &
          'Coordinate 3',1,IsoMesh % Nodes % z )

      Vfull => VariableGet( Mesh % Variables, 'Time' )
      CALL VariableAdd( Isomesh % Variables, Isomesh, Solver, 'Time', 1, &
          Vfull % Values )

      DO i=1,ints
        Vfull => VfullTab(i) % Variable
        NULLIFY( Vperm, Vals )
        ALLOCATE( Vperm(1),Vals(Vfull % DOFs) )
        Vperm = 0
        Vals = 0.0_dp
        CALL Info('IsoSurfaceSolver','Creating dummy variable '//TRIM( Vname(i) ),Level=8 )
        CALL VariableAddVector( Isomesh % Variables, Isomesh, Solver, &
             TRIM(Vname(i)), Vfull % DOFs, Vals, Vperm )
      END DO

      CALL GetIsoInterpolants()
      RETURN
    END IF


    ! Loop over edges and check for different signs
    !----------------------------------------------------------------
    DO IsAllocated=0,1

      j = 0            
      DO Level = 1, NoLevels
        IF( NoLevels > 1 ) LevelValue = LevelValues(Level,1) 
        
        DO i=1,NoEdges
          Edge => TempMesh % Edges(i)
          
          n1 = Edge % NodeIndexes(1)
          n2 = Edge % NodeIndexes(2)
          
          IF( ASSOCIATED( LevelPerm ) ) THEN
            m1 = LevelPerm( n1 )
            m2 = LevelPerm( n2 )
          ELSE
            m1 = n1
            m2 = n2
          END IF
          
          IF ( m1 <= 0 .OR. m2 <= 0 ) CYCLE
          
          IF( LevelDofs == 1 ) THEN
            t1 = LevelFun( m1 ) 
            t2 = LevelFun( m2 ) 
          ELSE
            t1 = 0.0_dp
            t2 = 0.0_dp
            DO k=1,LevelComps
              t1 = t1 + LevelFun( LevelDofs*(m1-1) + k )**2
              t2 = t2 + LevelFun( LevelDofs*(m2-1) + k )**2
            END DO
            t1 = SQRT( t1 )
            t2 = SQRT( t2 )
          END IF
          
          t1 = t1 - LevelValue
          t2 = t2 - LevelValue
          
          IF( ABS( t1 - t2 ) < TINY( t1 ) ) CYCLE
          IF ( t1 * t2 > 0.0_dp ) CYCLE
          
          j = j + 1
          IF( IsAllocated == 0 ) CYCLE

          !---------------------------------------------------------
          ! Only in the second loop tabulate the values
          !---------------------------------------------------------

          t = ABS( t1/(t2-t1) )

          ! The parent nodes the isosurface node sits between. ReMap needs
          ! them so that it can do without the temporal mesh, and they are
          ! also what says who owns the node in parallel: the node exists in
          ! every partition that holds the edge it sits on.
          !-------------------------------------------------------------------
          Interpolant(j) = t
          InvPerm(1,j) = n1
          InvPerm(2,j) = n2
          
          x1 = Mesh % Nodes % x(n1)
          x2 = Mesh % Nodes % x(n2)
          
          y1 = Mesh % Nodes % y(n1)
          y2 = Mesh % Nodes % y(n2)
          
          z1 = Mesh % Nodes % z(n1)
          z2 = Mesh % Nodes % z(n2)
          
          eperm(i) = j
          IsoEdge(j) = i
          IsoLevel(j) = Level
          Isomesh % Nodes % x(j) = (1-t) * x1 + t * x2
          Isomesh % Nodes % y(j) = (1-t) * y1 + t * y2
          Isomesh % Nodes % z(j) = (1-t) * z1 + t * z2
          
          DO k=1,ints
            Vfull => VfullTab(k) % Variable
            Viso => VisoTab(k) % Variable
            IF(.NOT. ASSOCIATED( Viso ) ) CYCLE

            IF( ASSOCIATED( Vfull % Perm ) ) THEN
              m1 = Vfull % Perm( n1 )
              m2 = Vfull % Perm( n2 )
            ELSE
              m1 = n1
              m2 = n2
            END IF
            IF( m1 <= 0 .OR. m2 <= 0 ) CYCLE
            
            DO l=1,Vfull % DOFs
              x1 = Vfull % Values(Vfull % DOFs*(m1-1)+l)
              x2 = Vfull % Values(Vfull % DOFs*(m2-1)+l)
              Viso % Values(Vfull % DOFs*(j-1)+l) = (1-t)*x1 + t*x2
            END DO
          END DO
        END DO
      END DO
      

      IF( IsAllocated == 0 ) THEN
        CALL Info('IsosurfaceSolver','Creating '//I2S(j)//' nodes for isosurface')

        NoIsoNodes = j
        Isomesh % NumberOfNodes = j
        Isomesh % Nodes % NumberOfNodes = j
        Isomesh % MeshDim = dim
       
        ALLOCATE( IsoMesh % Nodes % x(j) )
        ALLOCATE( IsoMesh % Nodes % y(j) )
        ALLOCATE( IsoMesh % Nodes % z(j) )

        ! Gives the index of the node sitting on a edge
        ALLOCATE( Eperm(NoEdges) ) 
        Eperm = 0

        ! ...and the other way round, which edge a node sits on. Two
        ! partitions sharing an edge make the same nodes on it, so this is
        ! what says who owns them.
        IF( ALLOCATED( IsoEdge ) ) DEALLOCATE( IsoEdge, IsoLevel )
        ALLOCATE( IsoEdge(j), IsoLevel(j) )
        IsoEdge = 0
        IsoLevel = 0
        
        CALL VariableAdd( IsoMesh % Variables, IsoMesh,Solver, &
            'Coordinate 1',1,IsoMesh % Nodes % x )
        
        CALL VariableAdd( IsoMesh % Variables,IsoMesh,Solver, &
            'Coordinate 2',1,IsoMesh % Nodes % y )
        
        CALL VariableAdd( IsoMesh % Variables,IsoMesh,Solver, &
            'Coordinate 3',1,IsoMesh % Nodes % z )
        
        Vfull => VariableGet( Mesh % Variables, 'Time' )
        CALL VariableAdd( Isomesh % Variables, Isomesh, Solver, 'Time', 1, &
            Vfull % Values )
        
        DO k=1,ints
          Vfull => VfullTab(k) % Variable
          NULLIFY( Vperm, Vals )
          ALLOCATE( Vperm(j),Vals(Vfull % DOFs*j) )
          Vperm = [(i,i=1,j)]
          Vals = 0.0_dp
          
          CALL Info('IsoSurfaceSolver','Creating variable '//TRIM( Vname(k) ) )
          CALL VariableAddVector( Isomesh % Variables, Isomesh, Solver, &
              TRIM(Vname(k)), Vfull % DOFs, Vals, Vperm )
        END DO

        ! Now that they exist, bind the isosurface side of the interpolants
        ! so that the value loop above needs no lookup per node.
        !----------------------------------------------------------------
        CALL GetIsoInterpolants()
	j = 0

        IF( ALLOCATED( InvPerm ) ) DEALLOCATE( InvPerm, Interpolant )
        ALLOCATE( InvPerm( 2, NoIsonodes ), Interpolant( NoIsoNodes) )
        InvPerm = 0
        Interpolant = 0.0_dp

      END IF

    END DO

  END FUNCTION CreateNodes


  !----------------------------------------------------------------
  !> Remap the field values to the isosurface.
  !----------------------------------------------------------------
  SUBROUTINE ReMap()

    INTEGER :: j,k,l,n1,n2,m1,m2
    REAL(KIND=dp) :: t,x1,x2

    TYPE(Variable_t), POINTER :: Vfull, Viso


    ! If there were no intersections then the dummy field has already been created
    !--------------------------------------------------------------------------------
    IF( .NOT. ALLOCATED( InvPerm ) ) RETURN

    ! Resolve the fields once, not per node.
    !----------------------------------------------------------------
    CALL GetInterpolants()
    CALL GetIsoInterpolants()

    ! If the mesh is moving then remap also the coordinate values
    !----------------------------------------------------------------
    IF( MovingMesh ) THEN
      DO j=1,Isomesh % NumberOfNodes
        n1 = InvPerm(1,j)
        n2 = InvPerm(2,j)
        t = Interpolant(j)

        Isomesh % Nodes % x(j) = &
            (1-t) * Mesh % Nodes % x(n1) + t * Mesh % Nodes % x(n2)
        Isomesh % Nodes % y(j) = &
            (1-t) * Mesh % Nodes % y(n1) + t * Mesh % Nodes % y(n2)
        Isomesh % Nodes % z(j) = &
            (1-t) * Mesh % Nodes % z(n1) + t * Mesh % Nodes % z(n2)
      END DO
    END IF

    ! Loop over the predefined nodes, field by field
    !----------------------------------------------------------------
    DO k=1,ints
      Vfull => VfullTab(k) % Variable
      Viso => VisoTab(k) % Variable
      IF (.NOT. ASSOCIATED(Viso)) CYCLE

      DO j=1,Isomesh % NumberOfNodes
        n1 = InvPerm(1,j)
        n2 = InvPerm(2,j)
        t = Interpolant(j)

        IF( ASSOCIATED( Vfull % Perm )) THEN
          m1 = Vfull % Perm(n1)
          m2 = Vfull % Perm(n2)
        ELSE
          m1 = n1
          m2 = n2
        END IF
        IF( m1 <= 0 .OR. m2 <= 0 ) CYCLE

        DO l=1,Vfull % DOFs
          x1 = Vfull % Values(Vfull % DOFs*(m1-1)+l)
          x2 = Vfull % Values(Vfull % DOFs*(m2-1)+l)
          Viso % Values(Viso % DOFs*(j-1)+l) = (1-t)*x1 + t*x2
        END DO
      END DO
    END DO

  END SUBROUTINE ReMap

  !----------------------------------------------------------------
  !> Create the isosurfaces or isolines.
  !----------------------------------------------------------------
  FUNCTION CreateSurfaces() RESULT ( NoSurfaces )
     INTEGER :: NoSurfaces
     TYPE(Element_t), POINTER :: Element, DefElement
     REAL(KIND=dp) :: F(4)
     INTEGER :: i,j,k,n,Triangles(2,3),Line(2)
     INTEGER :: NewElemType, NewElemNodes

     NoSurfaces = 0
     IF( NoIsoNodes == 0 ) RETURN

     IF( dim == 3 ) THEN
       NewElemType = 303
       NewElemNodes = 3
     ELSE
       NewElemType = 202
       NewElemNodes = 2
     END IF

     DO IsAllocated = 0, 1
       k = 0

       DO Level = 1, NoLevels
         IF( NoLevels > 1 ) LevelValue = LevelValues(Level,1) 

         DO i=1,TempMesh % NumberOfBulkElements
           Element => TempMesh % Elements(i)
           n = Element % TYPE % NumberOfNodes
           
           IF( ASSOCIATED( LevelPerm ) ) THEN
             ElemPerm(1:n) = LevelPerm(Element % NodeIndexes)
           ELSE
             ElemPerm(1:n) = Element % NodeIndexes
           END IF
           
           IF( LevelDofs == 1 ) THEN
             F(1:n) = LevelFun( ElemPerm(1:n) )
           ELSE
             F(1:n) = 0.0_dp
             DO j=1,LevelComps
               F(1:n) = F(1:n) + LevelFun( LevelDofs*(ElemPerm(1:n)-1)+j)**2
             END DO
             F(1:n) = SQRT( F(1:n) ) 
           END IF
           F(1:n) = F(1:n) - LevelValue
                      
           IF( dim == 2 ) THEN
             ! An edge that carries no intersection has eperm zero, which can
             ! happen for degenerate sign patterns. Skip those, as the 3D
             ! branch below does, rather than emit a zero node index.
             !-----------------------------------------------------------------
             IF( CreateLineFromTriangle(Element,F,Line) ) THEN
               IF( ALL( Line > 0 ) ) THEN
                 k = k + 1
                 IF( IsAllocated == 1) THEN
                   Isomesh % Elements(k) % NodeIndexes = Line
                   Isomesh % Elements(k) % BodyId = Level
                   Isomesh % Elements(k) % PartIndex = Element % PartIndex
                 END IF
               END IF
             END IF
           ELSE
             n = CreateSurfaceFromTetra(Element,F,Triangles)
             DO j=1,n
               IF ( ALL (Triangles(j,:) > 0 ) ) THEN
                 k = k + 1
                 IF( IsAllocated == 1) THEN
                   CALL OrientIsoTriangle( Element, F, Triangles(j,:) )
                   Isomesh % Elements(k) % NodeIndexes = Triangles(j,:)
                   Isomesh % Elements(k) % BodyId = Level
                   Isomesh % Elements(k) % PartIndex = Element % PartIndex
                 END IF
               END IF
             END DO
             
           END IF
         END DO
       END DO

       IF( k == 0 ) RETURN

       IF( IsAllocated == 0 ) THEN

	 NoSurfaces = k

         ALLOCATE( Isomesh % Elements(k) )
         
         Isomesh % MeshDim = dim 
         Isomesh % NumberOfBulkElements = k
         Isomesh % NumberOfFaces = 0
         Isomesh % NumberOfEdges = 0
         Isomesh % NumberOfBoundaryElements = 0
        
         DefElement => AllocateElement()
         DefElement % TYPE => GetElementType(NewElemType)

         ! The savers treat an element whose PartIndex is not their own as a
         ! halo element and skip it by default (SaveUtils.F90, 'Skip Halo
         ! Elements'), so the -1 that AllocateElement gives would write out
         ! empty pieces. Each element inherits the owner of the parent element
         ! it was cut from, below, so that a surface cut from a halo element
         ! stays the neighbour's to save rather than being counted twice.
         !----------------------------------------------------------------
         DefElement % PartIndex = ParEnv % MyPe

         DO i=1,k
           Isomesh % Elements(i) = DefElement
           Isomesh % Elements(i) % ElementIndex = i
           ALLOCATE(Isomesh % Elements(i) % NodeIndexes(NewElemNodes))
         END DO
       END IF
     END DO

     DEALLOCATE(DefElement)      
     
  END FUNCTION CreateSurfaces


  !----------------------------------------------------------------
  !> Point an isosurface triangle up the gradient of the level function.
  !>
  !> Unlike a surface in general, an isosurface is not free to be oriented
  !> either way: it is the level set of a scalar, so the side where the field
  !> is larger tells which way is out. Neighbouring elements agree about that
  !> by construction, being the same field, so this gives a consistent
  !> orientation over the whole surface without any walking or guessing.
  !>
  !> CreateSurfaceFromTetra picks its node triple by the sign pattern alone
  !> and pays no attention to the order within it, which leaves the triangles
  !> wound every which way, so put them right here.
  !----------------------------------------------------------------
  SUBROUTINE OrientIsoTriangle( Element, F, Tri )

    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: F(:)
    INTEGER :: Tri(3)

    REAL(KIND=dp) :: p(3,3), u(3), v(3), Nor(3), xp(3), xm(3)
    INTEGER :: i, np, nm, ni, itmp

    DO i=1,3
      p(1,i) = Isomesh % Nodes % x(Tri(i))
      p(2,i) = Isomesh % Nodes % y(Tri(i))
      p(3,i) = Isomesh % Nodes % z(Tri(i))
    END DO
    u = p(:,2) - p(:,1)
    v = p(:,3) - p(:,1)
    Nor(1) = u(2)*v(3) - u(3)*v(2)
    Nor(2) = u(3)*v(1) - u(1)*v(3)
    Nor(3) = u(1)*v(2) - u(2)*v(1)

    ! Where the field is above the level, and where below.
    xp = 0.0_dp; xm = 0.0_dp; np = 0; nm = 0
    DO i=1,Element % TYPE % NumberOfNodes
      ni = Element % NodeIndexes(i)
      IF( F(i) > 0.0_dp ) THEN
        np = np + 1
        xp(1) = xp(1) + Mesh % Nodes % x(ni)
        xp(2) = xp(2) + Mesh % Nodes % y(ni)
        xp(3) = xp(3) + Mesh % Nodes % z(ni)
      ELSE
        nm = nm + 1
        xm(1) = xm(1) + Mesh % Nodes % x(ni)
        xm(2) = xm(2) + Mesh % Nodes % y(ni)
        xm(3) = xm(3) + Mesh % Nodes % z(ni)
      END IF
    END DO
    IF( np == 0 .OR. nm == 0 ) RETURN

    IF( SUM( Nor * ( xp/np - xm/nm ) ) < 0.0_dp ) THEN
      itmp = Tri(2); Tri(2) = Tri(3); Tri(3) = itmp
    END IF

  END SUBROUTINE OrientIsoTriangle


  !----------------------------------------------------------------
  !> Create isosurfaces related to one tetrahedral element.
  !----------------------------------------------------------------
  FUNCTION CreateSurfaceFromTetra(Tetra,F,Surf) RESULT(scount)
    TYPE(Element_t) :: Tetra
    REAL(KIND=dp) :: F(4)
    INTEGER :: scount, Surf(2,3)

    LOGICAL :: S1,S2,S3,S4
    INTEGER :: S,Indexes(6)

    scount = 0
    Indexes = eperm(Tetra % EdgeIndexes)

    S1 = F(1) > 0.0_dp;
    S2 = F(2) > 0.0_dp;
    S3 = F(3) > 0.0_dp;
    S4 = F(4) > 0.0_dp;

    S = 0
    IF ( S1 ) S = S + 1
    IF ( S2 ) S = S + 1
    IF ( S3 ) S = S + 1
    IF ( S4 ) S = S + 1

    IF ( S==0 .OR. S==4 ) RETURN

    IF ( S==1 .OR. S==3 ) THEN
      scount = 1
      IF ( (S==1 .AND. S1) .OR. (S==3 .AND. .NOT.S1) ) THEN
        Surf(1,1) = Indexes(1)
        Surf(1,2) = Indexes(3)
        Surf(1,3) = Indexes(4)
      ELSE IF ( (S==1 .AND. S2) .OR. (S==3 .AND. .NOT.S2) ) THEN
        Surf(1,1) = Indexes(1)
        Surf(1,2) = Indexes(2)
        Surf(1,3) = Indexes(5)
      ELSE IF ( (S==1 .AND. S3) .OR. (S==3 .AND. .NOT.S3) ) THEN
        Surf(1,1) = Indexes(2)
        Surf(1,2) = Indexes(3)
        Surf(1,3) = Indexes(6)
      ELSE IF ( (S==1 .AND. S4) .OR. (S==3 .AND. .NOT.S4) ) THEN
        Surf(1,1) = Indexes(4)
        Surf(1,2) = Indexes(5)
        Surf(1,3) = Indexes(6)
      ELSE
        scount=0
      END IF
    ELSE
      scount = 2
      IF ( (S1 .AND. S2) .OR. (.NOT.S1 .AND. .NOT.S2) ) THEN
        Surf(1,1) = Indexes(3)
        Surf(1,2) = Indexes(2)
        Surf(1,3) = Indexes(5)

        Surf(2,1) = Indexes(3)
        Surf(2,2) = Indexes(5)
        Surf(2,3) = Indexes(4)
      ELSE IF ( (S1 .AND. S3) .OR. (.NOT.S1 .AND. .NOT.S3) ) THEN
        Surf(1,1) = Indexes(1)
        Surf(1,2) = Indexes(2)
        Surf(1,3) = Indexes(6)

        Surf(2,1) = Indexes(1)
        Surf(2,2) = Indexes(6)
        Surf(2,3) = Indexes(4)
      ELSE IF ( (S1 .AND. S4) .OR. (.NOT.S1 .AND. .NOT.S4) ) THEN
        Surf(1,1) = Indexes(1)
        Surf(1,2) = Indexes(5)
        Surf(1,3) = Indexes(6)

        Surf(2,1) = Indexes(1)
        Surf(2,2) = Indexes(6)
        Surf(2,3) = Indexes(3)
      ELSE
        scount=0
      END IF
    END IF

  END FUNCTION CreateSurfaceFromTetra

  
  !----------------------------------------------------------------
  !> Create isoline related to one triangular element.
  !----------------------------------------------------------------
  FUNCTION CreateLineFromTriangle(Triangle,F,Line) RESULT(GotLine)
    TYPE(Element_t) :: Triangle
    REAL(KIND=dp) :: F(3)
    INTEGER :: Line(2)
    LOGICAL :: GotLine

    LOGICAL :: S1,S2,S3
    INTEGER :: S,Indexes(3)

    GotLine = .FALSE.
    Indexes = eperm(Triangle % EdgeIndexes)

    S1 = F(1) > 0.0_dp;
    S2 = F(2) > 0.0_dp;
    S3 = F(3) > 0.0_dp;

    S = 0
    IF ( S1 ) S = S + 1
    IF ( S2 ) S = S + 1
    IF ( S3 ) S = S + 1

    IF ( S==0 .OR. S==3 ) RETURN

    GotLine = .TRUE.

    IF ( S == 1 ) THEN
      IF ( S1 ) THEN
        Line(1) = Indexes(1)
        Line(2) = Indexes(3)
      ELSE IF( S2 ) THEN
        Line(1) = Indexes(1)
        Line(2) = Indexes(2)
      ELSE
        Line(1) = Indexes(2)
        Line(2) = Indexes(3)
      END IF
    ELSE
      IF ( .NOT. S1 ) THEN
        Line(1)=Indexes(1)
        Line(2)=Indexes(3)
      ELSE IF( .NOT. S2 ) THEN
        Line(1) = Indexes(1)
        Line(2) = Indexes(2)
      ELSE
        Line(1) = Indexes(2)
        Line(2) = Indexes(3)
      END IF
    END IF

  END FUNCTION CreateLineFromTriangle



  ! Saves a loop in gmsh geo format 
  ! This is still not general and assumes one closed loop only!
  !-------------------------------------------------------------- 
  !----------------------------------------------------------------
  !> Collect the isoline on the master, as one mesh.
  !>
  !> The nodes are named by their global numbers, so the copies that the
  !> partitions keep of the ones sitting on their interfaces land on top of
  !> one another rather than being written out twice. Only the owner of a node
  !> sends its coordinates, and the segments are sent by whoever holds them,
  !> each belonging to exactly one partition.
  !----------------------------------------------------------------
  SUBROUTINE GatherGeoLines()

    INTEGER, ALLOCATABLE :: nRecv(:), Displ(:), gid(:), gidAll(:), &
        SegBuf(:), SegAll(:)
    REAL(KIND=dp), ALLOCATABLE :: xyz(:), xyzAll(:)
    INTEGER :: i,j,k,p,pe,nOwn,nLoc,ntot,ierr2
    TYPE(Element_t), POINTER :: Elem

    nGeoNodes = 0
    nGeoLines = 0
    IF( ALLOCATED( GeoX ) ) DEALLOCATE( GeoX, GeoY, GeoZ )
    IF( ALLOCATED( GeoLine ) ) DEALLOCATE( GeoLine )

    IF( ParEnv % PEs <= 1 ) THEN
      nGeoNodes = Isomesh % NumberOfNodes
      nGeoLines = Isomesh % NumberOfBulkElements
      ALLOCATE( GeoX(MAX(nGeoNodes,1)), GeoY(MAX(nGeoNodes,1)), &
          GeoZ(MAX(nGeoNodes,1)), GeoLine(2,MAX(nGeoLines,1)) )
      DO i=1,nGeoNodes
        GeoX(i) = Isomesh % Nodes % x(i)
        GeoY(i) = Isomesh % Nodes % y(i)
        GeoZ(i) = Isomesh % Nodes % z(i)
      END DO
      DO i=1,nGeoLines
        Elem => Isomesh % Elements(i)
        IF( Elem % TYPE % ElementCode /= 202 ) THEN
          CALL Fatal('SaveGmshGeo2D','Only elements of type 202 can be saved in Gmsh geo format!')
        END IF
        GeoLine(1:2,i) = Elem % NodeIndexes(1:2)
      END DO
      RETURN
    END IF

    pe = ParEnv % PEs

    ! The nodes this partition owns, and the total over all of them, which
    ! is also the largest global number in use.
    !----------------------------------------------------------------
    nOwn = 0
    IF( ASSOCIATED( Isomesh % ParallelInfo % NeighbourList ) ) THEN
      DO i=1,Isomesh % NumberOfNodes
        IF( Isomesh % ParallelInfo % NeighbourList(i) % Neighbours(1) == &
            ParEnv % MyPe ) nOwn = nOwn + 1
      END DO
    END IF

    ALLOCATE( gid(MAX(nOwn,1)), xyz(MAX(3*nOwn,1)) )
    k = 0
    IF( ASSOCIATED( Isomesh % ParallelInfo % NeighbourList ) ) THEN
      DO i=1,Isomesh % NumberOfNodes
        IF( Isomesh % ParallelInfo % NeighbourList(i) % Neighbours(1) /= &
            ParEnv % MyPe ) CYCLE
        k = k + 1
        gid(k) = Isomesh % ParallelInfo % GlobalDOFs(i)
        xyz(3*k-2) = Isomesh % Nodes % x(i)
        xyz(3*k-1) = Isomesh % Nodes % y(i)
        xyz(3*k)   = Isomesh % Nodes % z(i)
      END DO
    END IF

    ALLOCATE( nRecv(pe), Displ(pe) )
    CALL MPI_ALLGATHER( nOwn, 1, MPI_INTEGER, nRecv, 1, MPI_INTEGER, &
        ELMER_COMM_WORLD, ierr2 )
    ntot = SUM(nRecv)
    Displ(1) = 0
    DO p=2,pe
      Displ(p) = Displ(p-1) + nRecv(p-1)
    END DO
    ALLOCATE( gidAll(MAX(ntot,1)) )
    CALL MPI_GATHERV( gid, nOwn, MPI_INTEGER, gidAll, nRecv, Displ, &
        MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr2 )

    ALLOCATE( xyzAll(MAX(3*ntot,1)) )
    CALL MPI_GATHERV( xyz, 3*nOwn, MPI_DOUBLE_PRECISION, xyzAll, 3*nRecv, &
        3*Displ, MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr2 )

    ! The segments, each held by one partition only
    !----------------------------------------------------------------
    nLoc = 0
    DO i=1,Isomesh % NumberOfBulkElements
      IF( Isomesh % Elements(i) % PartIndex == ParEnv % MyPe ) nLoc = nLoc + 1
    END DO
    ALLOCATE( SegBuf(MAX(2*nLoc,1)) )
    k = 0
    DO i=1,Isomesh % NumberOfBulkElements
      Elem => Isomesh % Elements(i)
      IF( Elem % PartIndex /= ParEnv % MyPe ) CYCLE
      IF( Elem % TYPE % ElementCode /= 202 ) THEN
        CALL Fatal('SaveGmshGeo2D','Only elements of type 202 can be saved in Gmsh geo format!')
      END IF
      k = k + 1
      SegBuf(2*k-1) = Isomesh % ParallelInfo % GlobalDOFs(Elem % NodeIndexes(1))
      SegBuf(2*k)   = Isomesh % ParallelInfo % GlobalDOFs(Elem % NodeIndexes(2))
    END DO

    CALL MPI_ALLGATHER( nLoc, 1, MPI_INTEGER, nRecv, 1, MPI_INTEGER, &
        ELMER_COMM_WORLD, ierr2 )
    nGeoLines = SUM(nRecv)
    Displ(1) = 0
    DO p=2,pe
      Displ(p) = Displ(p-1) + nRecv(p-1)
    END DO
    ALLOCATE( SegAll(MAX(2*nGeoLines,1)) )
    CALL MPI_GATHERV( SegBuf, 2*nLoc, MPI_INTEGER, SegAll, 2*nRecv, 2*Displ, &
        MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr2 )

    IF( ParEnv % MyPe == 0 ) THEN
      nGeoNodes = 0
      DO i=1,ntot
        nGeoNodes = MAX( nGeoNodes, gidAll(i) )
      END DO
      ALLOCATE( GeoX(MAX(nGeoNodes,1)), GeoY(MAX(nGeoNodes,1)), &
          GeoZ(MAX(nGeoNodes,1)), GeoLine(2,MAX(nGeoLines,1)) )
      GeoX = 0.0_dp; GeoY = 0.0_dp; GeoZ = 0.0_dp
      DO i=1,ntot
        j = gidAll(i)
        GeoX(j) = xyzAll(3*i-2)
        GeoY(j) = xyzAll(3*i-1)
        GeoZ(j) = xyzAll(3*i)
      END DO
      DO i=1,nGeoLines
        GeoLine(1,i) = SegAll(2*i-1)
        GeoLine(2,i) = SegAll(2*i)
      END DO
      IF( ntot /= nGeoNodes ) THEN
        CALL Warn('SaveGmshGeo2D','Gathered '//I2S(ntot)//' nodes but the '//&
            'numbers run to '//I2S(nGeoNodes)//'!')
      END IF
    END IF

    DEALLOCATE( gid, xyz, gidAll, xyzAll, SegBuf, SegAll, nRecv, Displ )

  END SUBROUTINE GatherGeoLines


  SUBROUTINE SaveGmshGeo2D(Mesh)
    
    TYPE(Mesh_t) :: Mesh
    CHARACTER(LEN=MAX_NAME_LEN) :: Filename
    LOGICAL :: Found, SaveNode
    LOGICAL, ALLOCATABLE :: NodeUsed(:)
    INTEGER :: GeoUnit, iostat
    INTEGER :: i1,i2,i,j,UsedNodes
    INTEGER, ALLOCATABLE :: Neighbours(:,:),NodeOrder(:)
    REAL(KIND=dp) :: Dx,Dy,Dz,MeshDiam,MeshParam,Coord(3),NodeEps,PrevCoord(3),Dist

    INTEGER, PARAMETER :: MaxLoops = 20
    INTEGER :: LoopOffset(MaxLoops),LoopSize(MaxLoops),NoLoop,SaveLoops,NodeIndex 
    LOGICAL :: NewLoop


    ! Unlike the other formats a geo file is one object: a closed curve is a
    ! spline through an ordered run of points, so it cannot be written a
    ! partition at a time. Everything is therefore collected on the master,
    ! the nodes named by their global numbers so that the copies the
    ! partitions keep of the ones on their interfaces fall together again.
    !-----------------------------------------------------------------------
    CALL GatherGeoLines()

    IF( ParEnv % MyPe /= 0 ) RETURN
    IF( nGeoNodes == 0 ) RETURN

    Filename = ListGetString(Params,'Geo Filename',Found)
    IF( .NOT. Found ) Filename = 'mesh.geo'
    
    OPEN( NEWUNIT=GeoUnit, FILE=Filename, STATUS='UNKNOWN', IOSTAT=iostat )
    IF( iostat /= 0 ) THEN
      CALL Warn('SaveGmshGeo2D','Could not open file: '//TRIM(Filename))
      RETURN
    END IF
    
    n = nGeoNodes
    ALLOCATE( Neighbours(n, 2 ), NodeUsed( n ), NodeOrder( n ) ) 
    Neighbours = 0
    NodeUsed = .FALSE.
    NodeOrder = 0

    
    ! Create a list of neighbours.
    ! Each node should have exactly two neighbours.
    !-----------------------------------------------------------------------
    DO i=1,nGeoLines
      i1 = GeoLine(1,i)
      i2 = GeoLine(2,i)
      
      IF( i1 == 0 .OR. i2 == 0 ) THEN
        CALL Warn('SaveGmshGeo2D','Invalid indexes: '&
            //I2S(i1)//' and '//I2S(i2) )
      END IF
      
      IF( Neighbours(i1,1) == 0 ) THEN
        Neighbours(i1,1) = i2
      ELSE
        Neighbours(i1,2) = i2
      END IF
      
      IF( Neighbours(i2,1) == 0 ) THEN
        Neighbours(i2,1) = i1
      ELSE
        Neighbours(i2,2) = i1
      END IF
    END DO

    IF( ANY( Neighbours(:,2) == 0 ) ) THEN      
      CALL Warn('SaveGmshGeo2D','This does not seem to be a closed loop!')
      CLOSE( GeoUnit )
      RETURN
    END IF

    ! Compute the characteristic size of the bounding box 
    ! and get the mesh parameters.
    !--------------------------------------------------------------
    Dx = MAXVAL( GeoX ) - MINVAL( GeoX )
    Dy = MAXVAL( GeoY ) - MINVAL( GeoY )
    Dz = MAXVAL( GeoZ ) - MINVAL( GeoZ )
    MeshDiam = MAX( Dx, MAX( Dy, Dz ) )
    
    ! Currently defines a constant mesh parameter!
    MeshParam = ListGetCReal( Params,'Mesh Parameter',Found )
    IF( .NOT. Found ) MeshParam = MeshDiam / 50
    
    NodeEps = ListGetCReal( Params,'Mesh Node Epsilon',Found ) 
    IF( .NOT. Found ) NodeEps = 1.0e-3*MeshDiam

    ! Find the continuous loops and neglect points that are redundant 
    !-----------------------------------------------------------------
    j = 1
    UsedNodes = 0

    NoLoop = 1
    LoopOffset = 0
    LoopSize = 0
    NewLoop = .TRUE.

    DO i = 1, nGeoNodes

      ! Mark the current node before looking for the next one. Marking the
      ! candidate instead would leave the first node of each loop unmarked,
      ! so the walk could step back into it and start a phantom loop there.
      !-------------------------------------------------------------------
      NodeUsed(j) = .TRUE.

      Coord(1) = GeoX(j)
      Coord(2) = GeoY(j)
      Coord(3) = GeoZ(j)

      ! Check that the nodes are not so close to each other that Gmsh cannot handle the spline
      SaveNode = .TRUE.
      IF( .NOT. NewLoop ) THEN
        Dist = SUM( (Coord-PrevCoord)**2)
        IF( Dist < NodeEps**2 ) THEN
          WRITE(Message,'(A,ES12.4)') 'Skipping node '//I2S(j)//' too close: ',&
              SQRT(Dist)
          CALL Info('IsosurfaceSolver',Message,Level=12)
          SaveNode = .FALSE.
        END IF
      END IF

      IF( SaveNode ) THEN
        LoopSize(NoLoop) = LoopSize(NoLoop) + 1
        UsedNodes = UsedNodes + 1
        NodeOrder(UsedNodes) = j 
        PrevCoord = Coord
      END IF
      NewLoop = .FALSE.

      ! Now find the next candidate node
      ! This assumes closed loop!
      IF( .NOT. NodeUsed( Neighbours(j,1) ) ) THEN
        j = Neighbours(j,1)
      ELSE IF( .NOT. NodeUsed( Neighbours(j,2) ) ) THEN
        j = Neighbours(j,2)
      ELSE
        Found = .FALSE.
        DO j = 1, nGeoNodes
          IF( .NOT. NodeUsed(j) ) THEN
            CALL Info('IsosurfaceSolver','Found a new start at node: '//I2S(j),Level=10)
            
            Found = .TRUE.
            EXIT
          END IF
        END DO
        IF( Found ) THEN
          IF( NoLoop == MaxLoops ) THEN
            CALL Warn('IsosurfaceSolver','The static max number of loops ('&
                //I2S(MaxLoops)//') is too small, ignoring the rest!')
            EXIT
          END IF
          NoLoop = NoLoop + 1

          NewLoop = .TRUE.
          LoopOffset(NoLoop) = UsedNodes 
        ELSE
          CALL Info('IsosurfaceSolver','Could not find a new start, all nodes checked',Level=10)
          EXIT
        END IF
      END IF
    END DO

    IF( UsedNodes == nGeoNodes ) THEN
      CALL Info('IsosurfaceSolver','All nodes used')
    ELSE
      CALL Info('IsosurfaceSolver',I2S(UsedNodes)//' nodes used')
    END IF

    CALL Info('IsoSurfaceSolver','Number of loops found: '//I2S(NoLoop))

    CALL Info('IsosurfaceSolver','Writing points in geo file',Level=10)

    SaveLoops = ListGetInteger(Params,'Save Number Of Loops',Found ) 
    IF(.NOT. Found ) SaveLoops = NoLoop

    ! The loop below picks the largest loop not yet saved, so asking for more
    ! than exist would just re-save the last one.
    !-----------------------------------------------------------------------
    IF( SaveLoops > NoLoop ) THEN
      CALL Info('IsosurfaceSolver','Only '//I2S(NoLoop)//' loops to save',Level=8)
      SaveLoops = NoLoop
    END IF

    NodeIndex =  0
    DO l=1,SaveLoops

      ! Find the biggest unsaved loop
      n = 0
      DO k=1,NoLoop 
        IF( LoopSize(k) > n ) THEN
          j = k
          n = LoopSize(k)
        END IF
      END DO

      CALL Info('IsosurfaceSolver','Size of spline: '&
          //I2S(LoopSize(j)),Level=10)

      DO i=1,LoopSize(j)
        NodeIndex = NodeIndex + 1

        k = NodeOrder(LoopOffset(j) + i)
        Coord(1) = GeoX(k)
        Coord(2) = GeoY(k)
        Coord(3) = GeoZ(k)
       
        WRITE( GeoUnit,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A)') 'Point(',NodeIndex,') = {',&
            Coord(1),',',Coord(2),',',Coord(3),',',MeshParam,'};'
      END DO

      ! Ok, this was saved, earese the list
      LoopSize(j) = -LoopSize(j)
    END DO

    ! Revert from the negative values, to be able to use the logic again
    LoopSize = ABS( LoopSize ) 
    

    
    CALL Info('IsosurfaceSolver','Writing spline in geo file',Level=10)

    NodeIndex = 0
    DO l=1,SaveLoops

      ! Find the biggest unsaved loop
      n = 0
      DO k=1,NoLoop 
        IF( LoopSize(k) > n ) THEN
          j = k
          n = LoopSize(k)
        END IF
      END DO

      CALL Info('IsosurfaceSolver','Size of spline: '&
          //I2S(LoopSize(j)),Level=10)
      WRITE( GeoUnit,'(A)', ADVANCE='NO') 'Spline('//I2S(l)//') = {'
      DO i=1,LoopSize(j)
        WRITE( GeoUnit,'(I0,A)', ADVANCE='NO') NodeIndex+i,', '
      END DO
      WRITE( GeoUnit,'(I0,A)') NodeIndex+1,'};'
      NodeIndex = NodeIndex + LoopSize(j)

      ! Ok, this was saved, earese the list
      LoopSize(j) = -LoopSize(j)
    END DO

    DO j=1,SaveLoops
      k = SaveLoops+2*(j-1)
      WRITE( GeoUnit,'(A)') 'Line Loop('//I2S(k+1)//&
          ') = {'//I2S(j)//'};'
      WRITE( GeoUnit,'(A)') 'Plane Surface('//I2S(k+2)//&
          ') = {'//I2S(k+1)//'};'
    END DO

    CLOSE( GeoUnit )
    CALL Info('IsosurfaceSolver','Wrote geo file: '//TRIM(Filename),Level=6)

  END SUBROUTINE SaveGmshGeo2D

  
!------------------------------------------------------------------------------
END SUBROUTINE IsosurfaceSolver
!------------------------------------------------------------------------------


