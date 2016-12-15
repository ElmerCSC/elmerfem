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



!------------------------------------------------------------------------------
!> Subroutine for extracting isosurfaces in 3d and isolines in 2d.
!> The result will be a new mesh which will be added to the list of meshes.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE IsosurfaceSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------

  USE CoordinateSystems
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model    !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt       !< Timestep size for time dependent simulations
  LOGICAL :: Transient      !< Steady state or transient simulation
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: LevelVariableName
  TYPE(Mesh_t), POINTER :: Mesh, IsoMesh => NULL(), Pmesh, OrigMesh
  INTEGER :: i,j,k,l,n, dim,NoLevels,Level,NoNewElements
  REAL(KIND=dp) :: LevelValue, LevelDiff

  INTEGER, TARGET :: TetraToTetraMap(1,4), PyramidToTetraMap(2,4), &
             WedgeToTetraMap(3,4), BrickToTetraMap(5,4), &
             TriangleToTriangleMap(1,3), QuadToTriangleMap(2,3)
  INTEGER, POINTER :: Map(:,:), Indexes(:)
  INTEGER :: NoOrigElements,calls=0,ierr

  TYPE(Variable_t), POINTER :: LevelVariable
  LOGICAL :: MovingMesh, Found
  LOGICAL, ALLOCATABLE :: ElemSplit(:)
  INTEGER, ALLOCATABLE :: eperm(:), InvPerm(:), ElemPerm(:)
  INTEGER :: NewElemType, NewElemNodes, LevelDofs, NoEdges, IsAllocated, &
      NoIsoNodes, NoSurfaces
  INTEGER :: ParSizes(6), ParTmp(6),ElemFirst,ElemLast
  REAL(KIND=dp), ALLOCATABLE :: x(:),y(:),z(:),ElemFun(:)
  REAL(KIND=dp), POINTER :: LevelFun(:), LevelValues(:,:), Interpolant(:)
  INTEGER, POINTER :: LevelPerm(:)
  TYPE(Element_t), POINTER  :: OrigElements(:), NewElements(:), Element
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Visited = .FALSE., FixedSurface, SurfaceExist = .FALSE., &
      GotIsoMask

  SAVE Visited, FixedSurface, SurfaceExist, Interpolant, InvPerm, Isomesh, &
	NewElements, NoNewElements


  CALL Info( 'IsosurfaceSolver','-------------------------------------',Level=4 )
  CALL Info( 'IsosurfaceSolver','Determining the isosurface',Level=4 )
  CALL Info( 'IsosurfaceSolver','-------------------------------------',Level=4 )

  Mesh => GetMesh()
  OrigMesh => Mesh
  NoOrigElements = Mesh % NumberOfBulkElements
  OrigElements => Mesh % Elements
  NoNewElements = 0
  NewElements => NULL()

  Params => GetSolverParams()
  FixedSurface = GetLogical( Params,'Isosurface Fixed',Found)
  MovingMesh = GetLogical( Params,'Moving Mesh',Found)

  !---------------------------------------------------------------
  ! Get the isosurface variable
  !---------------------------------------------------------------
  LevelVariableName = GetString( Params,'Isosurface Variable')
  LevelVariable => VariableGet( Mesh % Variables, LevelVariableName )

  IF (.NOT.ASSOCIATED(LevelVariable)) THEN
    CALL Error( 'Isosurface', 'Missing isosurface variable: ' // &
        TRIM(LevelVariableName) )
    RETURN
  END IF

  Levelfun => LevelVariable % Values
  LevelPerm => LevelVariable % Perm
  LevelDofs = LevelVariable % Dofs
  
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
    RETURN
  END IF

  !--------------------------------------------------------------------------
  ! If the mesh was created and will be fixed make a simplified interpolation
  !--------------------------------------------------------------------------
  IF( FixedSurface .AND. SurfaceExist ) THEN
    CALL Info('IsosurfaceSolver','Remapping the fields')

    Mesh % NumberOfBulkElements = NoNewElements
    Mesh % Elements => NewElements

    CALL ReMap()

    Mesh % NumberOfBulkElements = NoOrigElements
    Mesh % Elements => OrigElements

    RETURN
  END IF

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
  CALL Info('IsosurfaceSolver','Leading element dimension is '//TRIM(I2S(dim)) )


  IF( dim == 2 ) THEN
    NewElemType = 303
    NewElemNodes = 3
  ELSE IF( dim == 3 ) THEN
    NewElemType = 504
    NewElemNodes = 4
  ELSE
    CALL Warn('IsoSurfaceSolver','Isosurface mesh can be created only for 2D or 3D')
    RETURN
  END IF


  !---------------------------------------------------------------
  ! Create a temporary mesh consisting only of tets (or triangles)
  ! in where the isosurface (or isoline) will be determined.
  ! The temporal mesh includes only potential master elements.
  !---------------------------------------------------------------
  CALL Info('IsoSurfaceSolver','Creating a temporal mesh')
  IF( NoLevels > 1 ) ALLOCATE( ElemSplit( Mesh % NumberOfBulkElements ) )
  
  DO IsAllocated = 0, 1
    
    j = 0
    IF( NoLevels > 1 ) ElemSplit = .FALSE.
    
    DO Level = 1, NoLevels
      IF( NoLevels > 1 ) LevelValue = LevelValues(Level,1) 
      
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
        
        IF( NoLevels > 1 ) THEN
          IF( ElemSplit(i) ) CYCLE
        END IF
        
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
          DO l=1,MIN(dim,LevelDofs)
            ElemFun(1:n) = ElemFun(1:n) + LevelFun( LevelDofs * ( ElemPerm(1:n) - 1) + l )**2
          END DO
          ElemFun(1:n) = SQRT( ElemFun(1:n) )
        END IF
        ElemFun(1:n) = ElemFun(1:n) - LevelValue
        IF( ALL ( ElemFun(1:n) < 0.0_dp ) ) CYCLE
        IF( ALL ( ElemFun(1:n) > 0.0_dp ) ) CYCLE
        
        IF( NoLevels > 1 ) THEN
          IF( ElemSplit(i) ) CYCLE
        END IF
        
        IF( IsAllocated == 0 ) THEN
          j = j + SIZE(Map,1)
          CYCLE
        END IF
        
        DO k=1,SIZE(Map,1)
          j = j + 1
          NewElements(j) = Element
          NewElements(j) % TYPE => GetElementType(NewElemType)
          ALLOCATE(NewElements(j) % NodeIndexes(NewElemNodes))
          DO l=1,NewElemNodes
            NewElements(j) % NodeIndexes(l) = &
                Element % NodeIndexes(Map(k,l))
          END DO
        END DO
      END DO
    END DO
    
    NoNewElements = j	
    
    IF( NoNewElements == 0 ) EXIT
    
    IF( IsAllocated == 0 ) THEN
      ALLOCATE(NewElements(NoNewElements))
    END IF    
  END DO

  IF( NoLevels > 1 ) DEALLOCATE( ElemSplit ) 
    
  IF( NoNewElements == 0 ) THEN
    CALL Warn('IsosurfaceSolver','No potential elements found for isosurface solver!')
    RETURN
  ELSE
    CALL Info('IsosurfaceSolver','Found '//TRIM(I2S(NoNewElements))//' potential elements') 
  END IF

  Mesh % NumberOfBulkElements = NoNewElements
  Mesh % Elements => NewElements

  ! Find mesh edges in order to define the intersection points
  !-----------------------------------------------------------
  NoEdges = 0
  CALL Info('IsosurfaceSolver','Creating mesh edges',Level=9)
  IF (.NOT.ASSOCIATED(Mesh % Edges)) THEN
    IF( dim == 2 ) THEN
      CALL FindMeshEdges2D(Mesh)
    ELSE
      CALL FindMeshEdges3D(Mesh)
    END IF
  END IF
  NoEdges = Mesh % NumberOfEdges

  ! If the mesh was previously created release it before recreating
  !----------------------------------------------------------------
  IF( SurfaceExist  ) THEN
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
    SurfaceExist = .FALSE.
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


  ! If requested number the output directories
  ! Alternatively one can number the output files
  !----------------------------------------------------------------
  IF( GetLogical(Params,'Mesh Numbering',Found) ) THEN
    Calls = Calls + 1
    Isomesh % name = TRIM(Isomesh % name) // TRIM(i2s(calls))
  END IF
  CALL MakeDirectory( TRIM(Isomesh % name) // CHAR(0) )
  
  ! Create nodes and elements on the edge intersections  
  !----------------------------------------------------------------
  CALL Info('IsosurfaceSolver','Creating nodes on edge intersections',Level=9)
  NoIsoNodes = CreateNodes()    

  CALL Info('IsosurfaceSolver','Creating surfaces or lines on edge intersections',Level=9)
  NoSurfaces = CreateSurfaces()
    
  ! Release temporary structures and revert original mesh
  !----------------------------------------------------------------------	
  IF(.NOT. FixedSurface ) THEN
    IF( NoEdges > 0 ) CALL ReleaseMeshEdgeTables( Mesh )
    IF ( ASSOCIATED( NewElements ) ) DEALLOCATE( NewElements )
  END IF

  IF( ALLOCATED( ElemFun ) ) DEALLOCATE( ElemFun, ElemPerm )  
  Mesh % NumberOfBulkElements = NoOrigElements
  Mesh % Elements => OrigElements

  ! Add the new mesh into the list 
  !----------------------------------------------------------------------
  PMesh => Model % Meshes
  DO WHILE( ASSOCIATED(PMesh % Next) )
    PMesh => PMesh % Next
  END DO
  PMesh % Next => Isomesh 

  ! Information of the new system size, also in parallel
  !----------------------------------------------------------------------
  ParTmp(1) = Mesh % NumberOfNodes
  ParTmp(2) = NoOrigElements 
  ParTmp(3) = NoNewElements 
  ParTmp(4) = NoEdges
  ParTmp(5) = NoIsoNodes
  ParTmp(6) = NoSurfaces

  IF( ParEnv % PEs > 1 ) THEN
    CALL MPI_ALLREDUCE(ParTmp,ParSizes,6,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
  ELSE
    ParSizes = ParTmp
  END IF

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


CONTAINS

  !----------------------------------------------------------------
  !> Create nodes for the isosurface / isoline.
  !----------------------------------------------------------------
  FUNCTION CreateNodes() RESULT ( NoIsoNodes )

    INTEGER :: NoIsoNodes
    TYPE(Element_t), POINTER :: Edge
    INTEGER :: i,j,k,l,n1,n2,m1,m2,ints,comps
    REAL(KIND=dp) :: t,t1,t2,x1,x2,y1,y2,z1,z2

    CHARACTER(LEN=MAX_NAME_LEN) :: Vname(32)
    LOGICAL :: Found
    TYPE(Variable_t), POINTER :: Vfull, Viso
    INTEGER, POINTER :: Vperm(:)
    REAL(KIND=dp), POINTER :: Vals(:)


    ! Make a table of the field variables and allocate space
    !----------------------------------------------------------------
    NoIsoNodes = 0
    ints = 0
    DO WHILE(.TRUE.)
      ints = ints+1
      Vname(ints) = GetString( Params, &
          ComponentName('Isosurface interpolant',ints), Found )
      IF(.NOT. Found) Vname(ints) = GetString( Params, &
          ComponentName('Variable',ints), Found )
      IF ( .NOT. Found ) THEN
        ints = ints-1
        EXIT
      END IF
      Vfull => VariableGet( Mesh % Variables, Vname(ints) )
      IF(.NOT. ASSOCIATED( Vfull ) ) THEN
        ints = ints - 1
        EXIT
      END IF
    END DO

    ! For parallel runs create the variables even though they do not exist since they 
    ! will do so in other partitions.
    !--------------------------------------------------------------------------------
    NoEdges = Mesh % NumberOfEdges
    IF( NoEdges == 0 ) THEN
      IF ( ints > 0 .AND. ParEnv % PEs > 1 ) THEN
        DO i=1,ints
          Vfull => VariableGet( Mesh % Variables, Vname(i) )
          NULLIFY( Vperm, Vals )
          ALLOCATE( Vperm(1),Vals(Vfull % DOFs) )
          Vperm = 0
          Vals = 0.0_dp
          CALL Info('IsoSurfaceSolver','Creating dummy variable '//TRIM( Vname(i) ) )
          CALL VariableAddVector( Isomesh % Variables, Isomesh, Solver, &
               TRIM(Vname(i)), Vfull % DOFs, Vals, Vperm )
        END DO
      END IF
      RETURN
    END IF


    ! Loop over edges and check for different signs
    !----------------------------------------------------------------
    DO IsAllocated=0,1

      j = 0            
      DO Level = 1, NoLevels
        IF( NoLevels > 1 ) LevelValue = LevelValues(Level,1) 
        
        DO i=1,NoEdges
          Edge => Mesh % Edges(i)
          
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
            DO k=1,MIN(dim,LevelDofs)
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

          IF( FixedSurface ) THEN
            Interpolant(j) = t
            InvPerm(j) = i
          END IF
          
          x1 = Mesh % Nodes % x(n1)
          x2 = Mesh % Nodes % x(n2)
          
          y1 = Mesh % Nodes % y(n1)
          y2 = Mesh % Nodes % y(n2)
          
          z1 = Mesh % Nodes % z(n1)
          z2 = Mesh % Nodes % z(n2)
          
          eperm(i) = j
          Isomesh % Nodes % x(j) = (1-t) * x1 + t * x2
          Isomesh % Nodes % y(j) = (1-t) * y1 + t * y2
          Isomesh % Nodes % z(j) = (1-t) * z1 + t * z2
          
          DO k=1,ints
            Vfull => VariableGet(Mesh % Variables, Vname(k))
            Viso => VariableGet(IsoMesh % Variables, Vname(k))

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
        CALL Info('IsosurfaceSolver','Creating '//TRIM(I2S(j))//' nodes for isosurface')

        NoIsoNodes = j
        ALLOCATE( IsoMesh % Nodes )
        Isomesh % NumberOfNodes = j
        Isomesh % Nodes % NumberOfNodes = j
        Isomesh % MeshDim = dim
       
        ALLOCATE( IsoMesh % Nodes % x(j) )
        ALLOCATE( IsoMesh % Nodes % y(j) )
        ALLOCATE( IsoMesh % Nodes % z(j) )

        ! Gives the index of the node sitting on a edge
        ALLOCATE( Eperm(NoEdges) ) 
        Eperm = 0
        
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
          Vfull => VariableGet( Mesh % Variables, Vname(k) )
          NULLIFY( Vperm, Vals )
          ALLOCATE( Vperm(j),Vals(Vfull % DOFs*j) )
          Vperm = [(i,i=1,j)]
          Vals = 0.0_dp
          
          CALL Info('IsoSurfaceSolver','Creating variable '//TRIM( Vname(k) ) )
          CALL VariableAddVector( Isomesh % Variables, Isomesh, Solver, &
              TRIM(Vname(k)), Vfull % DOFs, Vals, Vperm )
        END DO
	j = 0

        IF( FixedSurface ) THEN
          ALLOCATE( InvPerm( NoIsonodes ), Interpolant( NoIsoNodes) )
          InvPerm = 0
          Interpolant = 0.0_dp
        END IF

      END IF

    END DO

  END FUNCTION CreateNodes


  !----------------------------------------------------------------
  !> Remap the field values to the isosurface.
  !----------------------------------------------------------------
  SUBROUTINE ReMap()

    TYPE(Element_t), POINTER :: Edge
    INTEGER :: i,j,k,l,n1,n2,m1,m2,ints
    REAL(KIND=dp) :: t,x1,x2

    CHARACTER(LEN=MAX_NAME_LEN) :: Vname(32)
    LOGICAL :: Found
    TYPE(Variable_t), POINTER :: Vfull, Viso


    ! If there are no edges then the dummy field has already been created
    !--------------------------------------------------------------------------------
    IF( NoEdges == 0 ) RETURN

    ! Make a table of the field variables and allocate space
    !----------------------------------------------------------------
    n = NoEdges
    ints = 0
    DO WHILE(.TRUE.)
      ints = ints+1
      Vname(ints) = GetString( Params, &
          ComponentName('Isosurface interpolant',ints), Found )
      IF(.NOT. Found) Vname(ints) = GetString( Params, &
          ComponentName('Variable',ints), Found )
      IF ( .NOT. Found ) THEN
        ints = ints-1
        EXIT
      END IF
      Vfull => VariableGet( Mesh % Variables, Vname(ints) )
      IF(.NOT. ASSOCIATED( Vfull ) ) THEN
        ints = ints - 1
      END IF
    END DO

    ! Loop over predifined nodes
    !----------------------------------------------------------------
    DO Level = 1, NoLevels
      IF( NoLevels > 1 ) LevelValue = LevelValues(Level,1) 

      DO j=1,Isomesh % NumberOfNodes
        i = InvPerm(j)

        Edge => Mesh % Edges(i)
        
        n1 = Edge % NodeIndexes(1)
        n2 = Edge % NodeIndexes(2)
        
        t = Interpolant(j)

        ! If the mesh is moving then remap the also the coordinate values
        !----------------------------------------------------------------
        IF( MovingMesh ) THEN
          Isomesh % Nodes % x(j) = &
              (1-t) * Mesh % Nodes % x(n1) + t * Mesh % Nodes % x(n2)
          Isomesh % Nodes % y(j) = &
              (1-t) * Mesh % Nodes % y(n1) + t * Mesh % Nodes % y(n2)
          Isomesh % Nodes % z(j) = &
              (1-t) * Mesh % Nodes % z(n1) + t * Mesh % Nodes % z(n2)
        END IF

        
        DO k=1,ints
          Vfull => VariableGet(Mesh % Variables, Vname(k))
          IF (.NOT. ASSOCIATED(Vfull)) CYCLE

          Viso => VariableGet(IsoMesh % Variables, Vname(k))
          IF (.NOT. ASSOCIATED(Viso)) CYCLE

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

         DO i=1,Mesh % NumberOfBulkElements
           Element => Mesh % Elements(i)
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
             DO j=1,MIN(dim,LevelDofs)
               F(1:n) = F(1:n) + LevelFun( LevelDofs*(ElemPerm(1:n)-1)+j)**2
             END DO
             F(1:n) = SQRT( F(1:n) ) 
           END IF
           F(1:n) = F(1:n) - LevelValue
                      
           IF( dim == 2 ) THEN
             IF( CreateLineFromTriangle(Element,F,Line) ) THEN
               k = k + 1
               IF( IsAllocated == 1) THEN
                 Isomesh % Elements(k) % NodeIndexes = Line
                 Isomesh % Elements(k) % BodyId = Level
               END IF
             END IF
           ELSE
             n = CreateSurfaceFromTetra(Element,F,Triangles)
             DO j=1,n
               IF ( ALL (Triangles(j,:) > 0 ) ) THEN
                 k = k + 1
                 IF( IsAllocated == 1) THEN
                   Isomesh % Elements(k) % NodeIndexes = Triangles(j,:)
                   Isomesh % Elements(k) % BodyId = Level
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
         DO i=1,k
           Isomesh % Elements(i) = DefElement
!           Isomesh % Elements(i) % TYPE => GetElementType(NewElemType)
           ALLOCATE(Isomesh % Elements(i) % NodeIndexes(NewElemNodes))
         END DO
       END IF
     END DO

     DEALLOCATE(DefElement)      
     
  END FUNCTION CreateSurfaces


  !----------------------------------------------------------------
  !> Create isosurfaces related to one tetrahedral element.
  !----------------------------------------------------------------
  FUNCTION CreateSurfaceFromTetra(Tetra,F,Surf) RESULT(scount)
    TYPE(Element_t) :: Tetra
    REAL(KIND=dp) :: F(4)
    INTEGER :: scount, Surf(2,3)

    LOGICAL :: S1,S2,S3,S4
    INTEGER :: S,H,i,j,l,Indexes(6)

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
    INTEGER :: S,H,i,j,l,Indexes(3)

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
  SUBROUTINE SaveGmshGeo2D(Mesh)
    
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=MAX_NAME_LEN) :: Filename
    LOGICAL :: Found, SaveNode
    LOGICAL, ALLOCATABLE :: NodeUsed(:)
    INTEGER, PARAMETER :: GeoUnit = 10
    INTEGER :: i1,i2,i,j,jprev,jprevprev,UsedNodes,count
    INTEGER, ALLOCATABLE :: Neighbours(:,:),NodeOrder(:)
    REAL(KIND=dp) :: Dx,Dy,Dz,MeshDiam,MeshParam,Coord(3),NodeEps,PrevCoord(3),Dist

    INTEGER, PARAMETER :: MaxLoops = 20
    INTEGER :: LoopOffset(MaxLoops),LoopSize(MaxLoops),NoLoop,SaveLoops,NodeIndex 
    LOGICAL :: NewLoop


    IF( ParEnv % PEs > 1 ) THEN
      CALL Warn('SaveGsmhGeo2D','Not implemented yet in parallel')
    END IF
    
    Filename = ListGetString(Params,'Geo Filename',Found)
    IF( .NOT. Found ) Filename = 'mesh.geo'
    
    OPEN( UNIT=GeoUnit, FILE=Filename, STATUS='UNKNOWN') 
    
    n = Mesh % NumberOfNodes
    ALLOCATE( Neighbours(n, 2 ), NodeUsed( n ), NodeOrder( n ) ) 
    Neighbours = 0
    NodeUsed = .FALSE.
    NodeOrder = 0

    
    ! Create a list of neighbours.
    ! Each node should have exactly two neighbours.
    !-----------------------------------------------------------------------
    DO i=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(i)
      IF( Element % TYPE % ElementCode /= 202 ) THEN
        CALL Warn('SaveGmshGeo2D','Only elements of type 202 can be saved in Gmsh geo format!')
        RETURN
      END IF
      
      i1 = Element % NodeIndexes(1) 
      i2 = Element % NodeIndexes(2)
      
      IF( i1 == 0 .OR. i2 == 0 ) THEN
        CALL Warn('SaveGmshGeo2D','Invalid indexes: '&
            //TRIM(I2S(i1))//' and '//TRIM(I2S(i2)) )
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
      RETURN
    END IF

    ! Compute the characteristic size of the bounding box 
    ! and get the mesh parameters.
    !--------------------------------------------------------------
    Dx = MAXVAL( Mesh % Nodes % x ) - MINVAL( Mesh % Nodes % x )
    Dy = MAXVAL( Mesh % Nodes % y ) - MINVAL( Mesh % Nodes % y )
    Dz = MAXVAL( Mesh % Nodes % z ) - MINVAL( Mesh % Nodes % z )
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

    DO i = 1, Mesh % NumberOfNodes
      Coord(1) = Mesh % Nodes % x(j)
      Coord(2) = Mesh % Nodes % y(j)
      Coord(3) = Mesh % Nodes % z(j)

      ! Check that the nodes are not so close to each other that Gmsh cannot handle the spline
      SaveNode = .TRUE.
      IF( .NOT. NewLoop ) THEN
        Dist = SUM( (Coord-PrevCoord)**2)
        IF( Dist < NodeEps**2 ) THEN
          WRITE(Message,'(A,ES12.4)') 'Skipping node '//TRIM(I2S(j))//' too close: ',&
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
        DO j = 1, Mesh % NumberOfNodes
          IF( .NOT. NodeUsed(j) ) THEN
            CALL Info('IsosurfaceSolver','Found a new start at node: '//TRIM(I2S(j)),Level=10)
            
            Found = .TRUE.
            EXIT
          END IF
        END DO
        IF( Found ) THEN
          NoLoop = NoLoop + 1 
          IF( NoLoop > MaxLoops ) THEN
            CALL Warn('IsosurfaceSolver','The static max number of loops ('&
                //TRIM(I2S(MaxLoops))//') is too small!')
          END IF

          NewLoop = .TRUE.
          LoopOffset(NoLoop) = UsedNodes 
        ELSE
          CALL Info('IsosurfaceSolver','Could not find a new start, all nodes checked',Level=10)
          EXIT
        END IF
      END IF

      NodeUsed(j) = .TRUE.
    END DO

    IF( UsedNodes == Mesh % NumberOfNodes ) THEN
      CALL Info('IsosurfaceSolver','All nodes used')
    ELSE
      CALL Info('IsosurfaceSolver',TRIM(I2S(UsedNodes))//' nodes used')
    END IF

    CALL Info('IsoSurfaceSolver','Number of loops found: '//TRIM(I2S(NoLoop)))

    CALL Info('IsosurfaceSolver','Writing points in geo file',Level=10)

    SaveLoops = ListGetInteger(Params,'Save Number Of Loops',Found ) 
    IF(.NOT. Found ) SaveLoops = NoLoop

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
          //TRIM(I2S(LoopSize(j))),Level=10)

      DO i=1,LoopSize(j)
        NodeIndex = NodeIndex + 1

        k = NodeOrder(LoopOffset(j) + i)
        Coord(1) = Mesh % Nodes % x(k)
        Coord(2) = Mesh % Nodes % y(k)
        Coord(3) = Mesh % Nodes % z(k)
       
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
          //TRIM(I2S(LoopSize(j))),Level=10)
      WRITE( GeoUnit,'(A)', ADVANCE='NO') 'Spline('//TRIM(I2S(l))//') = {'
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
      WRITE( GeoUnit,'(A)') 'Line Loop('//TRIM(I2S(k+1))//&
          ') = {'//TRIM(I2S(j))//'};'
      WRITE( GeoUnit,'(A)') 'Plane Surface('//TRIM(I2S(k+2))//&
          ') = {'//TRIM(I2S(k+1))//'};'
    END DO


  END SUBROUTINE SaveGmshGeo2D

!------------------------------------------------------------------------------
END SUBROUTINE IsosurfaceSolver
!------------------------------------------------------------------------------


