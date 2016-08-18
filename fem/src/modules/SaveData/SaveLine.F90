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
! *  Subroutine for line data to files
! *
! ******************************************************************************
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 20 Nov 2001
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{
!------------------------------------------------------------------------------
!>  This subroutine saves 1D or 2D data in different formats.
!> Data on existing boundaries and polylines defined by set of coordinates 
!> may be saved.
!------------------------------------------------------------------------------
SUBROUTINE SaveLine( Model,Solver,dt,TransientSimulation )

  USE Types
  USE Lists
  USE Integration
  USE ElementDescription
  USE ElementUtils
  USE SolverUtils
  USE MeshUtils
  USE BandwidthOptimize
  USE DefUtils


  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultSideFile = 'sides.dat'

  REAL (KIND=DP), ALLOCATABLE ::  Values(:), Basis(:)
  REAL (KIND=dp) :: daxisx, daxisy, daxisz, x, y, z, eps, IntersectEpsilon, DetEpsilon, &
	f0, f1, f2, fn, q, weight
  REAL (KIND=DP), POINTER :: PointCoordinates(:,:), PointFluxes(:,:), PointWeight(:), Isosurf(:)
  LOGICAL :: Stat, GotIt, FileAppend, CalculateFlux, &
      SaveAxis(3), Inside, MovingMesh, IntersectEdge, OptimizeOrder, Found, GotVar, &
      SkipBoundaryInfo, GotDivisions
  INTEGER :: i,ii,j,k,ivar,l,n,m,t,DIM,mat_id, SaveThis, &
      Side, SaveNodes, SaveNodes2, node, NoResults, LocalNodes, NoVar, &
      No, axis, maxboundary, NoDims, MeshDim, NoLines, NoAxis, Line, NoFaces, &
      NoEigenValues, IntersectCoordinate, ElemCorners, ElemDim, FluxBody, istat, &
      i1, i2, NoTests
  INTEGER, POINTER :: NodeIndexes(:), SavePerm(:), InvPerm(:), BoundaryIndex(:), IsosurfPerm(:), NoDivisions(:)
  TYPE(Solver_t), POINTER :: ParSolver
  TYPE(Variable_t), POINTER :: Var, IsosurfVar
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Material, Params
  TYPE(Nodes_t) :: ElementNodes, LineNodes
  TYPE(Element_t), POINTER   :: CurrentElement
  CHARACTER(LEN=MAX_NAME_LEN) :: SideFile, SideNamesFile, VarName, Name, CondName, &
       TempName, MaskName, PrevMaskName, SideParFile, DateStr, OutputDirectory
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: ValueNames(:)

  LOGICAL, ALLOCATABLE :: LineTag(:)
  LOGICAL :: cand, Parallel, InitializePerm
  
  REAL(KIND=dp) :: R0(3),R1(3),dR(3),S0(3),S1(3),dS(3),LocalCoord(3),&
      MinCoord(3),MaxCoord(3),GlobalCoord(3),LineN(3),LineT1(3), &
      LineT2(3),detJ
  INTEGER :: imin,imax,nsize

  SAVE SavePerm, PrevMaskName, SaveNodes

!------------------------------------------------------------------------------

  CALL Info('SaveLine','-----------------------------------------', Level=4 )
  CALL Info('SaveLine','Saving data on specified lines',Level=4)

  Params => GetSolverParams()

  Mesh => Model % Meshes
  MeshDim = Mesh % MeshDim
  DO WHILE( ASSOCIATED(Mesh) )
    IF (Mesh % OutputActive) THEN
      CALL SetCurrentMesh( Model, Mesh )	
      EXIT 
    END IF
    Mesh => Mesh % Next
  END DO

  n = Mesh % MaxElementNodes
  ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), &
      LineNodes % x(2), LineNodes % y(2), LineNodes % z(2), &
      Basis(n), STAT=istat )     
  IF( istat /= 0 ) CALL Fatal('SaveLine','Memory allocation error for Elemental stuff') 

  IF( Solver % TimesVisited == 0 ) THEN
    ALLOCATE( SavePerm(Mesh % NumberOfNodes) )
  END IF

  
  Parallel = ( ParEnv % PEs > 1) 
  IF( GetLogical( Params,'Enforce Parallel Mode',GotIt) ) Parallel = .TRUE.	

  DIM = CoordinateSystemDimension()

  SideFile = ListGetString(Params,'Filename',GotIt )
  IF(.NOT. GotIt) SideFile = DefaultSideFile
  
  IF ( .NOT. FileNameQualified(SideFile) ) THEN
    OutputDirectory = GetString( Params,'Output Directory',GotIt) 
    IF( GotIt .AND. LEN_TRIM(OutputDirectory) > 0 ) THEN
      SideFile = TRIM(OutputDirectory)// '/' //TRIM(SideFile)
      IF( Solver % TimesVisited == 0 ) THEN
        CALL MakeDirectory( TRIM(OutputDirectory) // CHAR(0) )
      END IF
    ELSE IF( LEN_TRIM(OutputPath ) > 0 ) THEN
      SideFile = TRIM(OutputPath)// '/' //TRIM(SideFile)
    END IF
  END IF

  
  SideParFile = AddFilenameParSuffix(SideFile,'dat',Parallel,ParEnv % MyPe) 

  IF(ListGetLogical(Params,'Filename Numbering',GotIt)) THEN
    IF( Parallel ) THEN
      CALL Warn('SaveLine','Cannot number filenames in parallel with another number!')
    ELSE
      SideParFile = NextFreeFilename( SideParFile ) 
    END IF
  END IF

  CALL Info('SaveLine','Saving line data to file: '//TRIM(SideParFile),Level=12)

  
  FileAppend = ListGetLogical(Params,'File Append',GotIt )
  MovingMesh = ListGetLogical(Params,'Moving Mesh',GotIt )
 
  IntersectEdge = ListGetLogical(Params,'Intersect Edge',GotIt )
  IF(IntersectEdge) THEN
     IntersectCoordinate = ListGetInteger(Params,'Intersect Coordinate')
     IntersectEpsilon = ListGetConstReal(Params,'Intersect Epsilon')
  ELSE 
     IntersectCoordinate = 0
  END IF
  DetEpsilon = ListGetConstReal(Params,'Det Epsilon',GotIt)
  IF(.NOT. GotIt) DetEpsilon = 1.0e-6  

  CalculateFlux = ListGetLogical(Params,'Save Heat Flux',GotIt )
  IF(.NOT. CalculateFlux) THEN
     CalculateFlux = ListGetLogical(Params,'Save Flux',GotIt )
  END IF
  
  IF(CalculateFlux) THEN
     TempName = ListGetString(Params,'Flux Variable',GotIt )
     IF(.NOT. gotIt) TempName = TRIM('Temperature')
     CondName = ListGetString(Params,'Flux Coefficient',GotIt )
     IF(.NOT. gotIt) CondName = TRIM('Heat Conductivity')
  END IF
  
  SaveAxis(1) = ListGetLogical(Params,'Save Axis',GotIt)
  IF(GotIt) THEN
    SaveAxis(2:3) = SaveAxis(1)
  ELSE
    SaveAxis(1) = ListGetLogical(Params,'Save Axis 1',GotIt)
    SaveAxis(2) = ListGetLogical(Params,'Save Axis 2',GotIt)
    SaveAxis(3) = ListGetLogical(Params,'Save Axis 3',GotIt)    
  END IF
  NoAxis = DIM
  
  PointCoordinates => ListGetConstRealArray(Params,'Polyline Coordinates',gotIt)
  IF(gotIt) THEN
    NoLines = SIZE(PointCoordinates,1) / 2
    NoDims = SIZE(PointCoordinates,2)
    IF( NoDims < MeshDim ) THEN
      CALL Warn('SaveLine','Dimension of points smaller than that of mesh')
    END IF
  ELSE 
    NoLines = 0
  END IF

  GotDivisions = .FALSE.
  IF( NoLines > 0 ) THEN
    NoDivisions => ListGetIntegerArray( Params,'Polyline Divisions',GotDivisions)
    IF( SIZE( NoDivisions ) < NoLines + COUNT(SaveAxis) ) THEN
      CALL Fatal('SaveLine','Polyline divisions size too small!')
    END IF
  END IF

  SkipBoundaryInfo = ListGetLogical(Params,'Skip Boundary Info',GotIt)

!----------------------------------------------
! Specify the number of entries for each node
!---------------------------------------------- 
  IF( Solver % TimesVisited == 0  ) THEN
    CALL CreateListForSaving( Model, Params,.TRUE.,UseGenericKeyword = .TRUE.)    
  END IF

  NoVar = 0
  NoResults = 0
  DO ivar = 1,99
    Var => VariableGetN( ivar )
    IF ( .NOT. ASSOCIATED( Var ) )  EXIT
    NoVar = ivar
    IF (ASSOCIATED (Var % EigenVectors)) THEN
      NoEigenValues = SIZE(Var % EigenValues) 
      NoResults = NoResults + Var % Dofs * NoEigenValues
    ELSE
      NoResults = NoResults + Var % Dofs
    END IF
  END DO

  ! Add coordnate values
  NoResults = NoResults + 3

  IF ( CalculateFlux ) NoResults = NoResults + 3
  CALL Info('SaveLine','Maximum number of vectors to save: '//TRIM(I2S(NoResults)),Level=18)

  IF( NoVar == 0 .OR. NoResults == 0 ) THEN
    CALL Warn('SaveLine','No variables to save!')
    RETURN
  END IF

  ALLOCATE( Values(NoResults), STAT=istat )
  IF( istat /= 0 ) CALL Fatal('SaveLine','Memory allocation error for Values') 
 
 
  CALL Info( 'SaveLine', '------------------------------------------', Level=4 )
  WRITE( Message, * ) 'Saving line data to file ', TRIM(SideFile)
  CALL Info( 'SaveLine', Message, Level=4 )
  CALL Info( 'SaveLine', '------------------------------------------', Level=4 )

!------------------------------------------------------------------------------
! Open files for saving
!------------------------------------------------------------------------------
  IF( Solver % TimesVisited > 0 .OR. FileAppend) THEN 
    OPEN (10, FILE=SideParFile,POSITION='APPEND')
  ELSE 
    OPEN (10,FILE=SideParFile)
  END IF
  
  !------------------------------------------------------------------------------
  ! Find out which nodes should be saved and optimize their order 
  !------------------------------------------------------------------------------
  MaskName = ListGetString(Params,'Save Mask',GotIt) 
  IF(.NOT. GotIt) MaskName = 'Save Line'

  IF( Solver % TimesVisited > 0 ) THEN
    InitializePerm = ( MaskName /= PrevMaskName ) 
    InitializePerm = InitializePerm .OR. Solver % Mesh % Changed
  ELSE
    InitializePerm = .TRUE.
  END IF

  IF( InitializePerm ) THEN
    SavePerm = 0
    SaveNodes = 0

    OptimizeOrder = ListGetLogical(Params,'Optimize Node Ordering',GotIt)
    IF(.NOT. GotIt) OptimizeOrder = .NOT. Parallel

    CALL MakePermUsingMask( Model,Solver,Mesh,MaskName, &
        OptimizeOrder, SavePerm, SaveNodes, RequireLogical = .TRUE. )

    IF( SaveNodes > 0 ) THEN
      IF( ListGetLogical( Params,'Calculate Weights',GotIt ) ) THEN
        CALL CalculateNodalWeights( Solver, .TRUE., SavePerm, TRIM(MaskName)//' Weights')
      END IF
      CALL Info('SaveLine','Number of nodes in specified boundary: '//TRIM(I2S(SaveNodes)))
    END IF
  END IF
  PrevMaskName = MaskName

  !------------------------------------------------------------------------------
  ! If nodes found, then go through the sides and compute the fluxes if requested 
  !------------------------------------------------------------------------------

  maxboundary = 0
  IF( SaveNodes > 0 ) THEN

    ALLOCATE( InvPerm(SaveNodes), BoundaryIndex(SaveNodes), STAT=istat )
    IF( istat /= 0 ) CALL Fatal('SaveLine','Memory allocation error 3: '//TRIM(I2S(SaveNodes))) 

    BoundaryIndex = 0
    InvPerm = 0
    DO i=1,SIZE(SavePerm)
      IF (SavePerm(i)>0) THEN
        ! Error check for something that should never happen
        IF( InvPerm( SavePerm(i)) > 0) THEN
          WRITE( Message, *) 'Node multiple times in permutation',i,SavePerm(i)
          CALL Warn('SaveLine',Message)
        END IF
        InvPerm(SavePerm(i)) = i
      END IF
    END DO

    IF(CalculateFlux) THEN
      ALLOCATE(PointFluxes(SaveNodes,3),PointWeight(SaveNodes), STAT=istat)    
      IF( istat /= 0 ) CALL Fatal('SaveLine','Memory allocation error 4') 

      PointFluxes = 0.0d0
      PointWeight = 0.0d0
    END IF

    ! Go through the elements and register the boundary index and fluxes if asked
    DO t = 1,  Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements        

      CurrentElement => Mesh % Elements(t)
      Model % CurrentElement => CurrentElement
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes

      IF( .NOT. ALL(SavePerm(NodeIndexes) > 0)) CYCLE 

      Found = .FALSE.
      IF(t <= Mesh % NumberOfBulkElements) THEN
        k = CurrentElement % BodyId           
        IF( ListGetLogical( Model % Bodies(k) % Values,'Flux Integrate Body', gotIt )) THEN
          Found = .TRUE.
          FluxBody = ListGetInteger( Model % Bodies(k) % Values,'Flux Integrate Body', gotIt ) 
        END IF
      ELSE
        DO k=1, Model % NumberOfBCs
          IF ( Model % BCs(k) % Tag /= CurrentElement % BoundaryInfo % Constraint ) CYCLE
          IF( ListCheckPresent(Model % BCs(k) % Values, MaskName ) ) THEN
            Found = .TRUE.
            FluxBody = ListGetInteger( Model % BCs(k) % Values,'Flux Integrate Body', gotIt )         
            EXIT
          END IF
        END DO
      END IF
      IF(.NOT. Found) CYCLE

      MaxBoundary = MAX( MaxBoundary, k ) 
      BoundaryIndex( SavePerm(NodeIndexes) ) = k

      IF( .NOT. CalculateFlux ) CYCLE

      ElemCorners = CurrentElement % TYPE % ElementCode / 100
      IF(ElemCorners > 4 .OR. ElemCorners < 2) CYCLE

      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

      DO i=1,n
        node = NodeIndexes(i)

        CALL BoundaryFlux( Model, node, TempName,  &
            CondName, f1, f2, fn, weight, Mesh % MaxElementDOFs ) 

        j = SavePerm(node) 

        PointFluxes(j,1) = PointFluxes(j,1) + weight * f1
        PointFluxes(j,2) = PointFluxes(j,2) + weight * f2
        PointFluxes(j,3) = PointFluxes(j,3) + weight * fn
        PointWeight(j) = PointWeight(j) + weight
      END DO
    END DO

    ! Normalize flux by division with the integration weight
    IF( CalculateFlux ) THEN
      DO i = 1, SaveNodes
        PointFluxes(i,1) = PointFluxes(i,1) / PointWeight(i)
        PointFluxes(i,2) = PointFluxes(i,2) / PointWeight(i)
        PointFluxes(i,3) = PointFluxes(i,3) / PointWeight(i)
      END DO
      DEALLOCATE(PointWeight)
    END IF

    ! Save the nodes      
    !---------------     
    DO t = 1, SaveNodes    

      node = InvPerm(t)
      IF( .NOT. SkipBoundaryInfo ) THEN
        IF(TransientSimulation) WRITE(10,'(I6)',ADVANCE='NO') Solver % DoneTime
        WRITE(10,'(I6,I4,I8)',ADVANCE='NO') Solver % TimesVisited + 1,&
            BoundaryIndex(t),node
      END IF

      No = 0
      Values = 0.0d0

      DO ivar=-2, NoVar
        Var => VariableGetN( ivar ) 

        l = node
        IF (ASSOCIATED (Var % EigenVectors)) THEN
          NoEigenValues = SIZE(Var % EigenValues) 
          DO j=1,NoEigenValues
            DO i=1,Var % DOFs
              IF ( ASSOCIATED(Var % Perm) ) l = Var % Perm(l)
              IF(l > 0) THEN 
                Values(No+(j-1)*Var%Dofs+i) = &
                    Var % EigenVectors(j,Var%Dofs*(l-1)+i)
              END IF
            END DO
          END DO
          No = No + Var % Dofs * NoEigenValues
        ELSE           
          No = No + 1
          IF ( ASSOCIATED(Var % Perm) ) l = Var % Perm(l)
          IF(l > 0) Values(No) = Var % Values(l)                    
        END IF
      END DO

      IF(CalculateFlux) THEN
        Values(No+1:No+3) = PointFluxes(t,1:3)
      END IF

      IF(NoResults > 0) THEN
        DO i=1,NoResults-1
          WRITE(10,'(ES21.12E3)',ADVANCE='NO') Values(i)
        END DO
        WRITE(10,'(ES21.12E3)') Values(NoResults)
      END IF
    END DO

    DEALLOCATE(InvPerm, BoundaryIndex)
    IF(CalculateFlux) DEALLOCATE(PointFluxes)
  END IF

  !---------------------------------------------------------------------------
  ! Save data in the intersections of line segments defined by two coordinates
  ! and element faces, or save any of the principal axis.
  !---------------------------------------------------------------------------
  SaveNodes2 = 0  
  IF( NoLines > 0  .OR. ANY(SaveAxis(1:DIM) ) ) THEN
    NoTests = 0

    IF( .NOT. GotDivisions ) THEN
      IF( Solver % TimesVisited == 0 ) THEN 
        CALL FindMeshEdges( Mesh, .FALSE.)
      END IF
      
      IF(DIM == 2 .OR. IntersectEdge) THEN
        NoFaces = Mesh % NumberOfEdges
      ELSE 
        NoFaces = Mesh % NumberOfFaces
      END IF
    END IF

    IF( GotDivisions ) THEN
      t = MAXVAL( NoDivisions ) + 1
    ELSE
      t = Mesh % NumberOfNodes
    END IF

    ALLOCATE( LineTag(0:t), STAT=istat )
    IF( istat /= 0 ) CALL Fatal('SaveLine','Memory allocation error for LineTag') 
    
    
    DO Line = 1,NoLines + NoAxis

      LineTag = .FALSE.

      IF(Line <= NoLines) THEN
        LineNodes % x(1:2) = PointCoordinates(2*Line-1:2*Line,1) 
        LineNodes % y(1:2) = PointCoordinates(2*Line-1:2*Line,2) 
        IF(DIM == 3) THEN
          LineNodes % z(1:2) = PointCoordinates(2*Line-1:2*Line,3) 
        ELSE
          LineNodes % z(1:2) = 0.0d0
        END IF
      ELSE 
        IF(.NOT. SaveAxis(Line-NoLines)) CYCLE
        ! Define the lines for principal axis
        IF(Line-NoLines == 1) THEN
          LineNodes % x(1) = MINVAL(Mesh % Nodes % x)
          LineNodes % x(2) = MAXVAL(Mesh % Nodes % x)
          LineNodes % y(1:2) = 0.0d0
          LineNodes % z(1:2) = 0.0d0
        ELSE IF(Line-NoLines == 2) THEN
          LineNodes % x(1:2) = 0.0d0
          LineNodes % y(1) = MINVAL(Mesh % Nodes % y)
          LineNodes % y(2) = MAXVAL(Mesh % Nodes % y)
          LineNodes % z(1:2) = 0.0d0
        ELSE          
          LineNodes % x(1:2) = 0.0d0
          LineNodes % y(1:2) = 0.0d0
          LineNodes % z(1) = MINVAL(Mesh % Nodes % z)
          LineNodes % z(2) = MAXVAL(Mesh % Nodes % z)
        END IF
      END IF

      MaxBoundary = MaxBoundary + 1

      ! If we have specified number of divisions then use those
      IF( GotDivisions ) THEN
        R0(1) = LineNodes % x(1) 
        R0(2) = LineNodes % y(1) 
        R0(3) = LineNodes % z(1) 

        R1(1) = LineNodes % x(2) 
        R1(2) = LineNodes % y(2) 
        R1(3) = LineNodes % z(2) 

        dR = R1 - R0 
        LineN = dR / SQRT( SUM( dR**2 ) )
        CALL TangentDirections( LineN, LineT1, LineT2 ) 

        nsize = NoDivisions(Line)

        !PRINT *,'nsize:',Line,nsize
        !PRINT *,'N:',LineN
        !PRINT *,'T1:',LineT1
        !PRINT *,'T2:',LineT2

        DO t = 1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements       
          IF( t <= Mesh % NumberOfBulkElements ) THEN
            IF( IntersectEdge ) CYCLE
          ELSE
            IF( .NOT. IntersectEdge ) CYCLE
          END IF

          CurrentElement => Mesh % Elements(t)
          n = CurrentElement % Type % NumberOfNodes
          m = GetElementCorners( CurrentElement )
          NodeIndexes => CurrentElement % NodeIndexes

          ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
          ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

          DO i=1,m
            S1(1) = ElementNodes % x(i)  
            S1(2) = ElementNodes % y(i)  
            S1(3) = ElementNodes % z(i)  

            dS = ( S1 - R0 ) / SQRT( SUM( dR**2 ) )
            LocalCoord(1) = SUM( dS * LineN )
            LocalCoord(2) = SUM( dS * LineT1 )
            LocalCoord(3) = SUM( dS * LineT2 )

            IF( i == 1 ) THEN
              MinCoord = LocalCoord
              MaxCoord = LocalCoord
            ELSE
              DO j=1,3
                MinCoord(j) = MIN( MinCoord(j), LocalCoord(j) ) 
                MaxCoord(j) = MAX( MaxCoord(j), LocalCoord(j) ) 
              END DO
            END IF
          END DO

          IF( MinCoord(2) * MaxCoord(2) > 0.0_dp ) CYCLE
          IF( dim >= 3 .AND. .NOT. IntersectEdge ) THEN
            IF( MinCoord(3) * MaxCoord(3) > 0.0_dp ) CYCLE
          END IF

          imin = MAX(0, CEILING( nsize * MinCoord(1) ) )
          imax = MIN(nsize, FLOOR( ( nsize * MaxCoord(1) ) ) )


          DO i=imin,imax
            NoTests = NoTests + 1

            IF( LineTag(i) ) CYCLE

            GlobalCoord = R0 + i * dR / nsize

            IF ( PointInElement( CurrentElement, ElementNodes, GlobalCoord, LocalCoord ) ) THEN
              LineTag(i) = .TRUE.

              !PRINT *,'Hit:',t,i,GlobalCoord

              SaveNodes2 = SaveNodes2 + 1
              stat = ElementInfo( CurrentElement, ElementNodes, LocalCoord(1), &
                  LocalCoord(2), LocalCoord(3), detJ, Basis ) !, dBasisdx ) 

              IF( .NOT. SkipBoundaryInfo ) THEN
                IF(TransientSimulation) WRITE(10,'(I6)',ADVANCE='NO') Solver % DoneTime
                WRITE(10,'(I6,I4,I8)',ADVANCE='NO') Solver % TimesVisited+1,&
                    MaxBoundary,i
              END IF

              No = 0
              Values = 0.0d0

              DO ivar = -2,NoVar
                Var => VariableGetN( ivar ) 

                IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN
                  NodeIndexes => CurrentElement % DgIndexes
                ELSE
                  NodeIndexes => CurrentElement % NodeIndexes 
                END IF

                IF (ASSOCIATED (Var % EigenVectors)) THEN
                  NoEigenValues = SIZE(Var % EigenValues) 
                  DO j=1,NoEigenValues
                    DO k=1,n
                      l = NodeIndexes(k)
                      IF ( ASSOCIATED(Var % Perm) ) l = Var % Perm(l)
                      IF(l > 0) THEN 
                        DO ii=1,Var % DOFs
                          Values(No+(j-1)*Var%Dofs+ii) = Values(No+(j-1)*Var%Dofs+ii) + &
                              Basis(k) * (Var % EigenVectors(j,Var%Dofs*(l-1)+ii))
                        END DO
                      END IF
                    END DO
                  END DO
                  No = No + Var % Dofs * NoEigenValues
                ELSE                  
                  No = No + 1
                  DO k=1,n
                    l = NodeIndexes(k)
                    IF ( ASSOCIATED(Var % Perm) ) l = Var % Perm(l)
                    IF(l > 0) Values(No) = Values(No) + Basis(k) * (Var % Values(l))
                  END DO
                END IF
              END DO

              IF ( CalculateFlux ) Values(No+1:No+3) = 0.0

              IF(NoResults > 0) THEN
                DO j=1,NoResults-1
                  WRITE(10,'(ES20.11E3)',ADVANCE='NO') Values(j)
                END DO
                WRITE(10,'(ES20.11E3)') Values(NoResults)
              END IF
            END IF
          END DO
        END DO
      ELSE

        DO t = 1,NoFaces        
          IF(DIM == 2 .OR. IntersectEdge) THEN
            CurrentElement => Mesh % Edges(t)
          ELSE 
            CurrentElement => Mesh % Faces(t)
          END IF

          n = CurrentElement % TYPE % NumberOfNodes
          NodeIndexes => CurrentElement % NodeIndexes

          ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
          IF(DIM == 3) THEN
            ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
          ELSE
            ElementNodes % z(1:n) = 0.0d0
          END IF

          CALL GlobalToLocalCoords(CurrentElement,ElementNodes,n,LineNodes, &
              DetEpsilon,Inside,Basis,i)

          IF(.NOT. Inside) CYCLE

          ! When the line goes through a node it might be saved several times 
          ! without this checking
          IF(1.0d0-MAXVAL(Basis(1:n)) < 1.0d-3) THEN
            IF( LineTag(NodeIndexes(i)) ) CYCLE
            LineTag(NodeIndexes(i)) = .TRUE.
          END IF

          SaveNodes2 = SaveNodes2 + 1

          IF( .NOT. SkipBoundaryInfo ) THEN
            IF(TransientSimulation) WRITE(10,'(I6)',ADVANCE='NO') Solver % DoneTime
            WRITE(10,'(I6,I4,I8)',ADVANCE='NO') Solver % TimesVisited+1,&
                MaxBoundary,NodeIndexes(i)
          END IF

          No = 0
          Values = 0.0d0

          DO ivar = -2,NoVar
            Var => VariableGetN( ivar ) 

            IF (ASSOCIATED (Var % EigenVectors)) THEN
              NoEigenValues = SIZE(Var % EigenValues) 
              DO j=1,NoEigenValues
                DO k=1,n
                  l = NodeIndexes(k)
                  IF ( ASSOCIATED(Var % Perm) ) l = Var % Perm(l)
                  IF(l > 0) THEN 
                    DO i=1,Var % DOFs
                      Values(No+(j-1)*Var%Dofs+i) = Values(No+(j-1)*Var%Dofs+i) + &
                          Basis(k) * (Var % EigenVectors(j,Var%Dofs*(l-1)+i))
                    END DO
                  END IF
                END DO
              END DO
              No = No + Var % Dofs * NoEigenValues
            ELSE                  
              No = No + 1
              DO k=1,n
                l = NodeIndexes(k)
                IF ( ASSOCIATED(Var % Perm) ) l = Var % Perm(l)
                IF(l > 0) Values(No) = Values(No) + Basis(k) * (Var % Values(l))
              END DO
            END IF
          END DO

          IF ( CalculateFlux ) Values(No+1:No+3) = 0.0

          IF(NoResults > 0) THEN
            DO i=1,NoResults-1
              WRITE(10,'(ES20.11E3)',ADVANCE='NO') Values(i)
            END DO
            WRITE(10,'(ES20.11E3)') Values(NoResults)
          END IF
        END DO
      END IF
    END DO

    IF( NoTests > 0 ) THEN
      CALL Info('SaveLine','Number of candidate nodes: '//TRIM(I2S(NoTests)),Level=8)
    END IF

    CALL Info('SaveLine','Number of nodes in specified lines: '//TRIM(I2S(SaveNodes2)))

    DEALLOCATE( LineTag )
  END IF
  
  !---------------------------------------------------------------------------
  ! Save data in the intersections isocurves and element edges.
  !---------------------------------------------------------------------------

  IF( ListGetLogical( Params,'Save Isocurves',Found) ) THEN

    IF( DIM == 3 ) THEN
      CALL Fatal('SaveLine','Isocurves can only be saved in 2D')
    END IF

    ALLOCATE( LineTag(Mesh % NumberOfNodes), STAT=istat )
    IF( istat /= 0 ) CALL Fatal('SaveLine','Memory allocation error 6') 
        
    IF( Solver % TimesVisited == 0 ) THEN
      CALL FindMeshEdges( Mesh, .FALSE.)
    END IF
    
    NoFaces = Mesh % NumberOfEdges
    Line = 0

    DO WHILE( .TRUE. ) 
      
      Line = Line + 1
      LineTag = .FALSE.
      SaveNodes2 = 0
 
      WRITE (Name,'(A,I0)') 'IsoSurface Variable ',Line
      VarName = ListGetString( Params, Name, GotVar )
      
      WRITE (Name,'(A,I0)') 'IsoSurface Value ',Line
      f0 = ListGetCReal( Params, Name, Found )

      IF( GotVar ) THEN
        IsosurfVar => VariableGet( Mesh % Variables, VarName )
        IF( .NOT. ASSOCIATED( IsosurfVar ) ) THEN
          CALL Warn('SaveLine','Isosurface variable not given: '//TRIM(VarName))
          EXIT
        END IF
        IsosurfPerm => IsosurfVar % Perm
        Isosurf => IsosurfVar % Values       
      ELSE
        IF( Line == 1 ) THEN
          CALL Warn('SaveLine','No > Isosurface Variable 1 < defined!')
          EXIT
        END IF
        IF(.NOT. Found ) EXIT
      END IF
      
      WRITE( Message, * ) 'Finding nodes on isocurve: ',Line
      CALL Info('SaveLine',Message)
      
      f1 = MINVAL( Isosurf ) 
      f2 = MAXVAL( Isosurf ) 
      IF( f0 <= f1 .OR. f0 >= f2 ) THEN
        CALL Warn('SaveLine','Isosurface value not within range!')        
        PRINT *,'Range:',f1,f2,'f0:',f0
        CYCLE
      END IF

      
      MaxBoundary = MaxBoundary + 1

      DO t = 1,NoFaces        

        CurrentElement => Mesh % Edges(t)
        
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
        
        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = 0.0d0

        i1 = IsosurfPerm( NodeIndexes(1) )
        IF( i1 == 0 ) CYCLE
        f1 = Isosurf( i1 ) - f0

        i2 = IsosurfPerm( NodeIndexes(2) )
        IF( i2 == 0 ) CYCLE
        f2 = Isosurf( i2 ) - f0 

        ! There is an intersection if the value
        IF( f1 * f2 >= 0.0_dp ) CYCLE 
        
        q = ABS( f2 ) / ( ABS(f1) + ABS(f2) )

        Basis(1:n) = 0.0_dp
        Basis(1) = q
        Basis(2) = 1-q

        IF( Basis(1) > Basis(2) ) THEN
          k = NodeIndexes(1)
        ELSE
          k = Nodeindexes(2)
          q = 1-q
        END IF

        IF( q > 0.999 ) THEN
          IF( LineTag(k) ) CYCLE
          LineTag(k) = .TRUE.
        END IF
        SaveNodes2 = SaveNodes2 + 1
        
        IF( .NOT. SkipBoundaryInfo ) THEN
          IF(TransientSimulation) WRITE(10,'(I6)',ADVANCE='NO') Solver % DoneTime
          WRITE(10,'(I6,I4,I8)',ADVANCE='NO') Solver % TimesVisited+1,&
              MaxBoundary,k
        END IF
        
        No = 0
        Values = 0.0d0
        
        DO ivar = -2,NoVar
          Var => VariableGetN( ivar ) 
          
          IF (ASSOCIATED (Var % EigenVectors)) THEN
            NoEigenValues = SIZE(Var % EigenValues) 
            DO j=1,NoEigenValues
              DO k=1,n
                l = NodeIndexes(k)
                IF ( ASSOCIATED(Var % Perm) ) l = Var % Perm(l)
                IF(l > 0) THEN 
                  DO i=1,Var % DOFs
                    Values(No+(j-1)*Var%Dofs+i) = Values(No+(j-1)*Var%Dofs+i) + &
                        Basis(k) * (Var % EigenVectors(j,Var%Dofs*(l-1)+i))
                  END DO
                END IF
              END DO
            END DO
            No = No + Var % Dofs * NoEigenValues
          ELSE                  
            No = No + 1
            DO k=1,n
              l = NodeIndexes(k)
              IF ( ASSOCIATED(Var % Perm) ) l = Var % Perm(l)
              IF(l > 0) Values(No) = Values(No) + Basis(k) * (Var % Values(l))
            END DO
          END IF
        END DO
        
        IF ( CalculateFlux ) Values(No+1:No+3) = 0.0
        
        IF(NoResults > 0) THEN
          DO i=1,NoResults-1
            WRITE(10,'(ES20.11E3)',ADVANCE='NO') Values(i)
          END DO
          WRITE(10,'(ES20.11E3)') Values(NoResults)
        END IF
        
      END DO

      WRITE( Message, * ) 'Number of nodes in isocurve: ', SaveNodes2
      CALL Info('SaveLine',Message)
         
    END DO
    
    DEALLOCATE( LineTag )    
  END IF

  CLOSE(10)
  DEALLOCATE( Values )

  ! Finally save the names of the variables to help to identify the 
  ! columns in the result matrix.
  !-----------------------------------------------------------------
  !SaveNodes = NINT( ParallelReduction( 1.0_dp * SaveNodes ) )
  !SaveNodes2 = NINT( ParallelReduction( 1.0_dp * SaveNodes2 ) )

  IF( Solver % TimesVisited == 0 .AND. NoResults > 0 .AND. &
      (.NOT. Parallel .OR. ParEnv % MyPe == 0 ) ) THEN
    ALLOCATE( ValueNames(NoResults), STAT=istat )
    IF( istat /= 0 ) CALL Fatal('SaveLine','Memory allocation error for ValueNames') 

    No = 0
    DO ivar = -2,NoVar
      Var => VariableGetN( ivar ) 

      IF (ASSOCIATED (Var % EigenVectors)) THEN
        NoEigenValues = SIZE(Var % EigenValues) 
        DO j=1,NoEigenValues
          DO i=1,Var % DOFs
            IF(i==1) THEN
              WRITE(ValueNames(No+(j-1)*Var%Dofs+i),'(A,I0,A,A,I2,A,2ES20.11E3)') &
                  "Eigen ",j," ",TRIM(Var%Name),i,"   EigenValue = ",Var % EigenValues(j)
            ELSE 
              WRITE(ValueNames(No+(j-1)*Var%Dofs+i),'(A,I0,A,A,I2)') &
                  "Eigen ",j," ",TRIM(Var%Name),i
            END IF
          END DO
        END DO
        No = No + Var % Dofs * NoEigenValues
      ELSE 
        No = No + 1
        ValueNames(No) = TRIM(Var % Name)
      END IF
      Var => Var % Next      
    END DO

    IF ( CalculateFlux ) THEN
      ValueNames(No+1) = 'Flux 1'
      ValueNames(No+2) = 'Flux 2'
      ValueNames(No+3) = 'Flux normal'      
    END IF

    SideNamesFile = TRIM(SideFile) // '.' // TRIM("names")
    OPEN (10, FILE=SideNamesFile)

    Message = ListGetString(Model % Simulation,'Comment',GotIt)
    IF( GotIt ) THEN
      WRITE(10,'(A)') TRIM(Message)
    END IF
    Message = ListGetString(Params,'Comment',GotIt)
    IF( GotIt ) THEN
      WRITE(10,'(A)') TRIM(Message)
    END IF
    WRITE(10,'(A,A)') 'Variables in file: ',TRIM(SideFile)
    
    DateStr = FormatDate()
    WRITE( 10,'(A,A)') 'File started at: ',TRIM(DateStr)

    WRITE(10,'(I7,A)') SaveNodes,' boundary nodes for each step'
    WRITE(10,'(I7,A)') SaveNodes2,' polyline nodes for each step'
    j = 0
    IF( .NOT. SkipBoundaryInfo ) THEN
      IF(TransientSimulation) THEN
        WRITE(10,'(I3,": ",A)') 1,'Time step'
        j = 1
      END IF
      WRITE(10,'(I3,": ",A)') 1+j,'Iteration step'
      WRITE(10,'(I3,": ",A)') 2+j,'Boundary condition'
      WRITE(10,'(I3,": ",A)') 3+j,'Node index'
      j = j + 3
    END IF
    DO i=1,NoResults
      WRITE(10,'(I3,": ",A)') i+j,TRIM(ValueNames(i))
    END DO
    CLOSE(10)
    DEALLOCATE( ValueNames )
  END IF


  DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z, &
      LineNodes % x, LineNodes % y, LineNodes % z, Basis )

  CALL Info('SaveLine','All done')

CONTAINS


  FUNCTION VariableGetN(i) RESULT ( Var )
    INTEGER :: i
    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=MAX_NAME_LEN) :: Name, VarName
    LOGICAL :: Found

    NULLIFY(Var)
    
    IF( i < 1 ) THEN
      VarName = 'Coordinate '//TRIM(I2S(i+3))
      Found = .TRUE.
    ELSE
      WRITE (Name,'(A,I0)') 'Variable ',i
      VarName = GetString( Params, Name, Found )
    END IF
      
    IF( Found ) THEN
      Var => VariableGet( Model % Variables, VarName )
      IF(.NOT. ASSOCIATED( Var ) ) THEN
        CALL Warn('VariableGetN','Variable given but not found: '//TRIM(VarName))
      END IF
    END IF

  END FUNCTION VariableGetN
  


!---------------------------------------------------------------------------
!> This subroutine tests whether the line segment goes through the current
!> face of the element. If true the weights and index to the closest node 
!> are returned. 
!---------------------------------------------------------------------------
  SUBROUTINE GlobalToLocalCoords(Element,Plane,n,Line,Eps, &
      Inside,Weights,maxind)

    TYPE(Nodes_t) :: Plane, Line
    TYPE(Element_t), POINTER   :: Element
    INTEGER :: n, maxind
    REAL (KIND=dp) :: Eps,Weights(:)
    LOGICAL :: Inside

    REAL (KIND=dp) :: A(3,3),A0(3,3),B(3),C(3),Eps2,detA,absA,ds
    INTEGER :: split, i, corners, visited=0
    REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3)
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W=0.0d0

    SAVE visited
    visited = visited + 1

    Inside = .FALSE.
    corners = MIN(n,4)

    Eps2 = SQRT(TINY(Eps2))    

    ! In 2D the intersection is between two lines
    IF(DIM == 2) THEN
      A(1,1) = Line % x(2) - Line % x(1)
      A(2,1) = Line % y(2) - Line % y(1)
      A(1,2) = Plane % x(1) - Plane % x(2)
      A(2,2) = Plane % y(1) - Plane % y(2)
      A0 = A

      detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
      absA = SUM(ABS(A(1,1:2))) * SUM(ABS(A(2,1:2)))

      ! Lines are almost parallel => no intersection possible
      IF(ABS(detA) <= eps * absA + Eps2) RETURN

      B(1) = Plane % x(1) - Line % x(1) 
      B(2) = Plane % y(1) - Line % y(1) 

      CALL InvertMatrix( A,2 )
      C(1:2) = MATMUL(A(1:2,1:2),B(1:2))
     
      IF(ANY(C(1:2) < 0.0) .OR. ANY(C(1:2) > 1.0d0)) RETURN

      Inside = .TRUE.
      u = -1.0d0 + 2.0d0 * C(2)
    END IF


    IF(DIM == 3 .AND. IntersectCoordinate /= 0) THEN
      IF(IntersectCoordinate == 1) THEN
        A(1,1) = Line % y(2) - Line % y(1)
        A(2,1) = Line % z(2) - Line % z(1)
        A(1,2) = Plane % y(1) - Plane % y(2)
        A(2,2) = Plane % z(1) - Plane % z(2)
      ELSE IF(IntersectCoordinate == 2) THEN
        A(1,1) = Line % x(2) - Line % x(1)
        A(2,1) = Line % z(2) - Line % z(1)
        A(1,2) = Plane % x(1) - Plane % x(2)
        A(2,2) = Plane % z(1) - Plane % z(2)
      ELSE IF(IntersectCoordinate == 3) THEN
        A(1,1) = Line % x(2) - Line % x(1)
        A(2,1) = Line % y(2) - Line % y(1)
        A(1,2) = Plane % x(1) - Plane % x(2)
        A(2,2) = Plane % y(1) - Plane % y(2)
      ELSE
        PRINT *,'Intersect',IntersectCoordinate
        CALL Fatal('GlobalToLocalCoords','Impossible value for parameter IntersectCoordinate')
      END IF

      A0 = A
      
      detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
      absA = SUM(ABS(A(1,1:2))) * SUM(ABS(A(2,1:2)))

      IF(ABS(detA) <= eps * absA + Eps2) RETURN
     
      B(1) = Plane % x(1) - Line % x(1) 
      B(2) = Plane % y(1) - Line % y(1) 

      CALL InvertMatrix( A,2 )
      C(1:2) = MATMUL(A(1:2,1:2),B(1:2))

!      PRINT *,'C',C(1:2)
!      PRINT *,'xp=',Plane % x(1) + C(2) * (Plane % x(2) - Plane % x(1)) 
!      PRINT *,'xl=',Line % x(1) + C(1) * (Line % x(2) - Line % x(1)) 
!      PRINT *,'yp=',Plane % y(1) + C(2) * (Plane % y(2) - Plane % y(1)) 
!      PRINT *,'yl=',Line % y(1) + C(1) * (Line % y(2) - Line % y(1)) 

      IF(ANY(C(1:2) < 0.0) .OR. ANY(C(1:2) > 1.0d0)) RETURN

      IF(IntersectCoordinate == 1) THEN
        ds = Line % x(1) + C(1)* (Line % x(2) - Line % x(1))  &
            - Plane % x(1) - C(2) * (Plane % x(1) - Plane % x(2))
      ELSE IF(IntersectCoordinate == 2) THEN
        ds = Line % y(1) + C(1)* (Line % y(2) - Line % y(1))  &
            - Plane % y(1) - C(2) * (Plane % y(1) - Plane % y(2))
      ELSE 
        ds = Line % z(1) + C(1)* (Line % z(2) - Line % z(1))  &
            - Plane % z(1) - C(2) * (Plane % z(1) - Plane % z(2))      
      END IF

      IF(ABS(ds) > IntersectEpsilon) RETURN

      Inside = .TRUE.
      u = -1.0d0 + 2.0d0 * C(2)
    END IF   

    
    ! In 3D rectangular faces are treated as two triangles
    IF(DIM == 3 .AND. IntersectCoordinate == 0) THEN

      DO split=0,corners-3
         
        A(1,1) = Line % x(2) - Line % x(1)
        A(2,1) = Line % y(2) - Line % y(1)
        A(3,1) = Line % z(2) - Line % z(1)

        IF(split == 0) THEN
          A(1,2) = Plane % x(1) - Plane % x(2)
          A(2,2) = Plane % y(1) - Plane % y(2)
          A(3,2) = Plane % z(1) - Plane % z(2)
        ELSE 
          A(1,2) = Plane % x(1) - Plane % x(4)
          A(2,2) = Plane % y(1) - Plane % y(4)
          A(3,2) = Plane % z(1) - Plane % z(4)
        END IF

        A(1,3) = Plane % x(1) - Plane % x(3)
        A(2,3) = Plane % y(1) - Plane % y(3)
        A(3,3) = Plane % z(1) - Plane % z(3)
        
        ! Check for linearly dependent vectors
        detA = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
             - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
             + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
        absA = SUM(ABS(A(1,1:3))) * SUM(ABS(A(2,1:3))) * SUM(ABS(A(3,1:3))) 

        IF(ABS(detA) <= eps * absA + Eps2) CYCLE
!        print *,'detA',detA

        B(1) = Plane % x(1) - Line % x(1)
        B(2) = Plane % y(1) - Line % y(1)
        B(3) = Plane % z(1) - Line % z(1)
        
        CALL InvertMatrix( A,3 )
        C(1:3) = MATMUL( A(1:3,1:3),B(1:3) )
        
        IF( ANY(C(1:3) < 0.0) .OR. ANY(C(1:3) > 1.0d0) ) CYCLE
        IF(C(2)+C(3) > 1.0d0) CYCLE

        Inside = .TRUE. 

        ! Relate the point of intersection to local coordinates
        IF(corners < 4) THEN
          u = C(2)
          v = C(3)
        ELSE IF(corners == 4 .AND. split == 0) THEN
          u = 2*(C(2)+C(3))-1
          v = 2*C(3)-1
        ELSE 
          ! For the 2nd split of the rectangle the local coordinates switched
          v = 2*(C(2)+C(3))-1
          u = 2*C(3)-1        
        END IF

        IF(Inside) EXIT
        
      END DO
    END IF

    IF(.NOT. Inside) RETURN

    stat = ElementInfo( Element, Plane, U, V, W, SqrtElementMetric, &
        Basis, dBasisdx )
    
    ! This check seems to work only for linear elements
    ! and seems to be futile probably there too.
    ! IF( n == Element % TYPE % ElementCode / 100 ) THEN
    !  IF(MAXVAL(Basis(1:n)-1.0) > eps .OR. MINVAL(Basis(1:n)) < -eps) THEN
    !    Inside = .FALSE.
    !    RETURN
    !  END IF
    !END IF

    Weights(1:n) = Basis(1:n)
    MaxInd = 1
    DO i=2,n
      IF(Weights(MaxInd) < Weights(i)) MaxInd = i
    END DO

  END SUBROUTINE GlobalToLocalCoords
  

!-----------------------------------------------------------------------
!> Compuation of normal flux.
!> Note that this is calculated on the nodal points only
!> using a single boundary element. The direction of the normal
!> may be somewhat different on the nodal point when calculated using 
!> a neighboring boundary element.
!> Thus normal flow calculation is useful only when the boundary 
!> is relatively smooth. Also quadratic elements are recommended.
!-----------------------------------------------------------------------
   
  SUBROUTINE BoundaryFlux( Model, Node, VarName, CoeffName, f1, f2, fn, weight, MaxN) 
    USE Types
    USE Lists
    USE ElementDescription
    
    TYPE(Model_t) :: Model
    INTEGER :: dimno,i,j,n,node,lbody,rbody,MaxN
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, CoeffName
    REAL(KIND=dp) :: f1, f2, fn, weight
    
    TYPE(Variable_t), POINTER :: Tvar
    TYPE(Element_t), POINTER :: Parent, Element, OldCurrentElement
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: r,u,v,w,ub, &
         Basis(MaxN), dBasisdx(MaxN,3),DetJ, Normal(3)
    REAL(KIND=dp), TARGET :: x(MaxN), y(MaxN), z(MaxN)
    LOGICAL :: stat, Permutated
    INTEGER :: body_id, k
    REAL(KIND=dp) :: Conductivity(Model % Mesh % MaxElementNodes)
    REAL(KIND=DP), POINTER :: Pwrk(:,:,:)
    REAL(KIND=DP) :: CoeffTensor(3,3,Model % MaxElementNodes), Flow(3)


    Tvar => VariableGet( Model % Variables, TRIM(VarName) )

    Permutated = ASSOCIATED(Tvar % Perm)
    Element => Model % CurrentElement
    
    IF ( FluxBody > 0 ) THEN
      lbody = 0
      IF ( ASSOCIATED( Element % BoundaryInfo % Left ) ) &
        lbody = Element % BoundaryInfo % Left % BodyId

      rbody = 0
      IF ( ASSOCIATED( Element % BoundaryInfo % Right ) ) &
        rbody = Element % BoundaryInfo % Right % BodyId

      IF ( LBody == FluxBody ) THEN
        Parent => Element % BoundaryInfo % Left
      ELSEIF ( RBody == FluxBody ) THEN
        Parent => Element % BoundaryInfo % Right
      ELSE
        WRITE( Message, * ) 'No such flux integrate body on bc ', &
            Element % BoundaryInfo % Constraint
        CALL Fatal( 'SaveLine', Message )
      END IF
    ELSE        
      Parent => Element % BoundaryInfo % Left
      stat = ASSOCIATED( Parent )

      IF(Permutated) THEN
        IF(stat) stat = ALL(TVar % Perm(Parent % NodeIndexes) > 0)
        
        IF ( .NOT. stat ) THEN
          Parent => ELement % BoundaryInfo % Right
          
          stat = ASSOCIATED( Parent )
          IF(stat) stat = ALL(TVar % Perm(Parent % NodeIndexes) > 0)
        END IF
      END IF
      IF ( .NOT. stat )  CALL Fatal( 'SaveLine',&
          'No solution available for specified boundary' )
    END IF
    
    OldCurrentElement => Element
    Model % CurrentElement => Parent

    n = Parent % TYPE % NumberOfNodes
    
    Nodes % x => x
    Nodes % y => y
    Nodes % z => z
    Nodes % x(1:n) = Model % Nodes % x(Parent % NodeIndexes)
    Nodes % y(1:n) = Model % Nodes % y(Parent % NodeIndexes)
    Nodes % z(1:n) = Model % Nodes % z(Parent % NodeIndexes)
    
    DO j=1,n
      IF ( node == Parent % NodeIndexes(j) ) EXIT
    END DO

    IF ( node /= Parent % NodeIndexes(j) ) THEN
      CALL Warn('SaveLine','Side node not in parent element!')
    END IF

    CALL GlobalToLocal( u, v ,w , x(j), y(j), z(j), Parent, Nodes )

    stat = ElementInfo( Parent, Nodes, u, v, w, detJ, &
        Basis, dBasisdx )
    weight = detJ
   
    ! Compute the normal of the surface for the normal flux
    DO j = 1, Element % TYPE % NumberOfNodes
      IF ( node == Element % NodeIndexes(j) ) EXIT
    END DO

    IF ( j == 1 ) THEN
      ub = -1.0d0
    ELSEIF ( j == 2 ) THEN
      ub = 1.0d0
    ELSE
      ub = 0.0d0
    END IF

    Normal = Normalvector( Element, ElementNodes, ub, 0.0d0, .TRUE. )

    body_id = Parent % Bodyid
    k = ListGetInteger( Model % Bodies(body_id) % Values,'Material', &
            minv=1, maxv=Model % NumberOFMaterials )
    Material => Model % Materials(k) % Values


    CALL ListGetRealArray( Material, TRIM(CoeffName), Pwrk, n, &
        Parent % NodeIndexes, GotIt )

    CoeffTensor = 0.0d0
    IF(GotIt) THEN
      IF ( SIZE(Pwrk,1) == 1 ) THEN
        DO i=1,3
          CoeffTensor( i,i,1:n ) = Pwrk( 1,1,1:n )
        END DO
      ELSE IF ( SIZE(Pwrk,2) == 1 ) THEN
        DO i=1,MIN(3,SIZE(Pwrk,1))
          CoeffTensor(i,i,1:n) = Pwrk(i,1,1:n)
        END DO
      ELSE
        DO i=1,MIN(3,SIZE(Pwrk,1))
          DO j=1,MIN(3,SIZE(Pwrk,2))
            CoeffTensor( i,j,1:n ) = Pwrk(i,j,1:n)
          END DO
        END DO
      END IF
    END IF
    
    Flow = 0.0d0
    DO j = 1, DIM
      DO k = 1, DIM
        IF(Permutated) THEN
          Flow(j) = Flow(j) + SUM( CoeffTensor(j,k,1:n) * Basis(1:n) ) * &
              SUM( dBasisdx(1:n,k) * TVar % Values(TVar % Perm(Parent % NodeIndexes(1:n))) )
        ELSE
          Flow(j) = Flow(j) + SUM( CoeffTensor(j,k,1:n) * Basis(1:n) ) * &
              SUM( dBasisdx(1:n,k) * TVar % Values(Parent % NodeIndexes(1:n)) ) 
        END IF
      END DO
    END DO

    f1 = Flow(1)
    f2 = Flow(2)
    fn = SUM(Normal(1:DIM) * Flow(1:DIM))

    Model % CurrentElement => OldCurrentElement

  END SUBROUTINE BoundaryFlux

!------------------------------------------------------------------------------
END SUBROUTINE SaveLine
!------------------------------------------------------------------------------
!> \}
