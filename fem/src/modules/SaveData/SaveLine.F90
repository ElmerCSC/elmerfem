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


SUBROUTINE SaveLine_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  INTEGER :: NormInd,i
  LOGICAL :: GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  INTEGER, POINTER :: ASolvers(:)
  LOGICAL :: ActiveSomewhere


  ! Create Perm vector if the solver is not active everywhere.
  ! This will affect how Polylines, Circle lines and Isocurves are
  ! plotted. 
  ActiveSomewhere = .FALSE.
  DO i=1,CurrentModel % NumberOfEquations
    ASolvers => ListGetIntegerArray( Model % Equations(i) % Values, &
        'Active Solvers', GotIt )
    IF( GotIt ) THEN
      IF( ANY(ASolvers == Solver % SolverId) ) THEN
        ActiveSomewhere = .TRUE.
        EXIT
      END IF
    END IF
  END DO

  IF(ActiveSomewhere) THEN
    IF( .NOT. ListCheckPresent( Solver % Values,'Variable') ) THEN
      CALL Info('SaveLine_init','If you want to mask with Equation block give some "Variable"',Level=5)
    END IF
  END IF
  
  ! If we want to show a pseudonorm add a variable for which the norm
  ! is associated with.
  SolverName = ListGetString( Solver % Values, 'Equation',GotIt)
  NormInd = ListGetInteger( Solver % Values,'Show Norm Index',GotIt)
  IF( NormInd > 0 ) THEN
    CALL ListAddNewString( Solver % Values,'Variable',&
        '-nooutput -global '//TRIM(SolverName)//'_var')
  END IF
  
  CALL ListAddNewLogical( Solver % Values,'No Matrix',.TRUE.)
  
END SUBROUTINE SaveLine_init


!------------------------------------------------------------------------------
!> This subroutine saves 1D or 2D data in different formats.
!> Data on existing boundaries, polylines and circle defined by set of coordinates 
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
  USE ElementUtils
  USE SaveUtils
  USE BandwidthOptimize
  USE DefUtils
  

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
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
  REAL (KIND=DP), POINTER :: PointCoordinates(:,:), Isosurf(:)
  REAL(KIND=dp), ALLOCATABLE :: PointFluxes(:,:), PointWeight(:)
  LOGICAL :: Stat, GotIt, FileAppend, CalculateFlux, &
      SaveAxis(3), Inside, MovingMesh, IntersectEdge, OptimizeOrder, Found, GotVar, &
      SkipBoundaryInfo, GotDivisions, EdgeBasis, DGVar, ElemVar, IpVar, MaskWithPerm
  INTEGER :: i,ii,j,k,ivar,l,n,m,t,DIM,mat_id, SaveThis, &
      Side, SaveNodes(4), node, NoResults, NoLabels, &
      LocalNodes, NoVar, No, axis, maxboundary, NoDims, MeshDim, NoLines, NoAxis, Line, &
      NoFaces, NoEigenValues, IntersectCoordinate, ElemCorners, ElemDim, istat, &
      i1, i2, NoTests, NormInd, Comps, SaveSolverMeshIndex, LineInd, FoundNan
  INTEGER, POINTER :: NodeIndexes(:), SavePerm(:) => NULL(), InvPerm(:), &
      BoundaryIndex(:), IsosurfPerm(:), NoDivisions(:)
  TYPE(Solver_t), POINTER :: ParSolver
  TYPE(Variable_t), POINTER :: Var, Var2, Var3, IsosurfVar
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Material, Params
  TYPE(Nodes_t) :: ElementNodes, LineNodes
  TYPE(Element_t), POINTER   :: CurrentElement
  CHARACTER(LEN=MAX_NAME_LEN) :: SideFile, SideNamesFile, VarName, Name, CondName, &
       TempName, MaskName, PrevMaskName, SideParFile, DateStr
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: ValueNames(:)
  CHARACTER(:), ALLOCATABLE :: OutputDirectory

  LOGICAL, ALLOCATABLE :: LineTag(:)
  LOGICAL :: cand, Parallel, InitializePerm, FileIsOpen, AVBasis, DoEigen
  
  REAL(KIND=dp) :: R0(3),R1(3),dR(3),S0(3),S1(3),dS(3),LocalCoord(3),&
      MinCoord(3),MaxCoord(3),GlobalCoord(3),LineN(3),LineT1(3), &
      LineT2(3),detJ, Norm
  INTEGER :: imin,imax,nsize,LineUnit
  CHARACTER(*), PARAMETER :: Caller = 'SaveLine'

  INTEGER :: NoData
  REAL(KIND=dp), POINTER :: PosData(:) 
  INTEGER, POINTER :: LabelData(:,:)
  REAL(KIND=dp), POINTER :: ResultData(:,:)
  
  INTERFACE
    SUBROUTINE Ip2DgFieldInElement( Mesh, Element, nip, fip, ndg, fdg )
      USE Types
      USE Integration
      USE ElementDescription
      IMPLICIT NONE
      
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Element_t), POINTER :: Element
      INTEGER :: nip, ndg
      REAL(KIND=dp) :: fip(:), fdg(:)
    END SUBROUTINE Ip2DgFieldInElement
  END INTERFACE
  
  SAVE SavePerm, PrevMaskName, SaveNodes

!------------------------------------------------------------------------------

  Params => GetSolverParams()
  SideFile = ListGetString(Params,'Filename',GotIt )
  IF(.NOT. GotIt) SideFile = DefaultSideFile

  CALL Info(Caller,'-----------------------------------------', Level=4 )
  CALL Info(Caller,'Saving data on lines to file: '//TRIM(SideFile),Level=4)
  CALL Info(Caller, '------------------------------------------', Level=4 )

  FileIsOpen = .FALSE.
  FoundNan = 0
  NoData = 0
  SaveNodes = 0
  MaskWithPerm = .FALSE.
  IF( ASSOCIATED( Solver % Variable ) ) THEN
    MaskWithPerm = ASSOCIATED( Solver % Variable % Perm )
  END IF
  
  NULLIFY( PosData, LabelData, ResultData ) 
  
  i = GetInteger( Params,'Save Solver Mesh Index',Found ) 
  IF( Found ) THEN
    CALL Info(Caller,'Using mesh of solver '//I2S(i))
    Mesh => Model % Solvers(i) % Mesh
    Model % Mesh => Mesh
  ELSE
    Mesh => GetMesh()
  END IF

  
  DIM = CoordinateSystemDimension()
  MeshDim = Mesh % MeshDim
  Parallel = ( ParEnv % PEs > 1) 
  IF( GetLogical( Params,'Enforce Parallel Mode',GotIt) ) Parallel = .TRUE.	

  n = Mesh % MaxElementNodes
  ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), &
      LineNodes % x(2), LineNodes % y(2), LineNodes % z(2), &
      Basis(n), STAT=istat )     
  IF( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for Elemental stuff') 

  IF( Solver % TimesVisited == 0 ) THEN
    ALLOCATE( SavePerm(Mesh % NumberOfNodes), STAT=istat )
    IF( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for SavePerm') 
  END IF
  
  NormInd = ListGetInteger( Params,'Show Norm Index',GotIt)
  Norm = 0.0_dp

  MovingMesh = ListGetLogical(Params,'Moving Mesh',GotIt )

  IF( DIM == 3 ) THEN
    IntersectEdge = ListGetLogical(Params,'Intersect Edge',GotIt )   
    IF(.NOT. IntersectEdge ) THEN
      IntersectCoordinate = ListGetInteger(Params,'Intersect Coordinate',IntersectEdge)
    ELSE
      IntersectCoordinate = ListGetInteger(Params,'Intersect Coordinate')
    END IF
  ELSE
    IntersectCoordinate = 0
    IntersectEdge = .FALSE.
  END IF
  IF(IntersectEdge) THEN
    IntersectEpsilon = ListGetConstReal(Params,'Intersect Epsilon')
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
  
  LineInd = ListGetInteger( Params,'Line Marker',GotIt)
  
  
!----------------------------------------------
! Specify the number of entries for each node
!---------------------------------------------- 
  IF( Solver % TimesVisited == 0  ) THEN
    CALL CreateListForSaving( Model, Params,.TRUE.,UseGenericKeyword = .TRUE.)    
  END IF

  AVBasis = .FALSE.
  NoVar = 0
  NoResults = 0
  NoLabels = 0
  DO ivar = 1,99
    Var => VariableGetN( ivar, Comps )
    IF ( .NOT. ASSOCIATED( Var ) )  EXIT
    NoVar = ivar

    DGvar = .FALSE.
    IpVar = .FALSE. 
    ElemVar = .FALSE.
    EdgeBasis = .FALSE.
    
    IF (ASSOCIATED (Var % EigenVectors)) THEN
      NoEigenValues = SIZE(Var % EigenValues) 
      NoResults = NoResults + Var % Dofs * NoEigenValues
    ELSE
      IF( Var % TYPE == variable_on_nodes_on_elements ) THEN
        DGVar = .TRUE.
      ELSE IF( Var % TYPE == variable_on_gauss_points ) THEN
        IpVar = .TRUE.
      ELSE IF( Var % TYPE == Variable_on_elements ) THEN
        ElemVar = .TRUE.
      ELSE IF( ASSOCIATED( Var % Solver ) .AND. ASSOCIATED( Var % Perm ) ) THEN
        EdgeBasis = GetLogical( Var % Solver % Values,'Hcurl Basis',Found )
        IF( EdgeBasis ) THEN          
          EdgeBasis = ( SIZE( Var % Perm ) > Mesh % NumberOfNodes )
        END IF
        IF( EdgeBasis ) THEN
          EdgeBasis = ( ANY( Var % Perm( Mesh % NumberOfNodes+1:) > 0 ) )
        END IF
        IF( EdgeBasis ) THEN
          IF( ANY( Var % Perm(1: Mesh % NumberOfNodes) > 0 ) ) AVBasis = .TRUE. 
        END IF
        IF( EdgeBasis ) Var % TYPE = Variable_on_edges
      END IF

      IF( EdgeBasis ) THEN
        IF( AVBasis ) THEN
          CALL Info(Caller,'Variable '//I2S(ivar)//' is treated as living in nodal+Hcurl: '//TRIM(Var % Name),Level=10)
          NoResults = NoResults + 4
        ELSE          
          NoResults = NoResults + 3
          CALL Info(Caller,'Variable '//I2S(ivar)//' is treated as living in Hcurl: '//TRIM(Var % Name),Level=10)
        END IF
      ELSE
        IF( DgVar ) THEN
          CALL Info(Caller,'Variable '//I2S(ivar)//' is treated as living on DGBasis: '//TRIM(Var % Name),Level=10)
        END IF
        IF( IpVar ) THEN
          CALL Info(Caller,'Variable '//I2S(ivar)//' is treated as living in IP points: '//TRIM(Var % Name),Level=10)
        END IF
        IF( ElemVar ) THEN
          CALL Info(Caller,'Variable '//I2S(ivar)//' is treated as living on elements: '//TRIM(Var % Name),Level=10)
        END IF
        NoResults = NoResults + MAX( Var % Dofs, Comps ) 
      END IF
    END IF
  END DO
  
  IF ( CalculateFlux ) NoResults = NoResults + 3
  CALL Info(Caller,'Maximum number of fields for each node: '//I2S(NoResults),Level=18)

  IF( NoVar == 0 .OR. NoResults == 0 ) GOTO 1

  ! Add coordnate values
  MaxBoundary = 0
  NoResults = NoResults + 3

  ALLOCATE( Values(NoResults), STAT=istat )
  IF( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for Values') 
 
  SkipBoundaryInfo = ListGetLogical(Params,'Skip Boundary Info',GotIt)

  ! Open just the directory for saving so that it is created for certain!
  !------------------------------------------------------------------------------
  CALL SolverOutputDirectory( Solver, SideFile, OutputDirectory )
  IF( Solver % TimesVisited == 0 ) THEN
    i = 1; i = ParallelReduction(i)
  END IF
  
  ! Search existing boundary to save if any
  !------------------------------------------------------------------------------
  CALL SaveExistingLines()

  ! Save data in the intersections isocurves and element edges.
  !---------------------------------------------------------------------------
  CALL SaveIsocurves()

  ! Save data in the intersections of line segments defined by two coordinates
  ! and element faces, or save any of the principal axis.
  !---------------------------------------------------------------------------
  CALL SavePolylines()

  CALL SaveCircleLines()
  
  ! Finally close the file for saving
  !----------------------------------------
1 CALL CloseLineFile()

  CALL SaveVariableNames()

  DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z, &
      LineNodes % x, LineNodes % y, LineNodes % z, STAT=istat)
  IF(istat /= 0) CALL Fatal(Caller,'Errors in deallocating some basis stuff!')

  IF( ALLOCATED(Values) ) DEALLOCATE(Values)
  IF( ALLOcATED(Basis) ) DEALLOCATE(Basis)
  
  IF( NormInd > 0 ) THEN    
    Norm = ParallelReduction(Norm) 
    Solver % Variable % Values = Norm
    Solver % Variable % Norm = Norm
  END IF
  
  IF(FoundNan > 0 ) THEN
    CALL Warn(Caller,'Replaced '//I2S(FoundNan)//' NaN entries with -1')
  END IF
  
  CALL Info(Caller,'All done')

CONTAINS


  ! Get the Nth variable. The coordinate is a cludge since 
  ! the coordinate is not in the automated variable list. 
  !---------------------------------------------------------------
  FUNCTION VariableGetN(i,NoComponents,Component) RESULT ( Var )
    TYPE(Variable_t), POINTER :: Var
    INTEGER, OPTIONAL :: NoComponents
    INTEGER, OPTIONAL :: Component

    INTEGER :: i,j,k
    CHARACTER(LEN=MAX_NAME_LEN) :: Name, VarName
    LOGICAL :: Found
    TYPE(Variable_t), POINTER :: Var2

    NULLIFY(Var)
    
    IF( i < 1 ) THEN
      VarName = 'Coordinate '//I2S(i+3)
      Found = .TRUE.
    ELSE
      WRITE (Name,'(A,I0)') 'Variable ',i
      VarName = GetString( Params, Name, Found )
    END IF

    k = 0
    IF(.NOT. Found ) THEN
      IF( PRESENT( NoComponents ) ) NoComponents = k
      RETURN
    END IF
      
    IF( PRESENT( Component ) ) THEN
      VarName = TRIM(VarName)//' '//I2S(Component)
    END IF
      
    Var => VariableGet( Mesh % Variables, VarName )
    IF( .NOT. ASSOCIATED( Var ) ) THEN
      Var => VariableGet( Mesh % Variables, TRIM(VarName)//' 1' )
      IF( ASSOCIATED( Var ) ) THEN
        DO j=2,99
          Var2 => VariableGet( Mesh % Variables, TRIM(VarName)//' '//I2S(j) )
          IF(ASSOCIATED( Var2 ) ) THEN
            k = j
          ELSE
            EXIT
          END IF
        END DO
      ELSE
        CALL Warn('VariableGetN','Variable given but not found: '//TRIM(VarName))
      END IF
    END IF

    IF( PRESENT( NoComponents ) ) NoComponents = k

    CALL Info('VariableGetN','Variable: '//TRIM(VarName)//': '//I2S(k),Level=31)
 
  END FUNCTION VariableGetN


!---------------------------------------------------------------------------
!> This subroutine tests whether the line segment goes through the current
!> face of the element. If true the weights and index to the closest node 
!> are returned. 
!---------------------------------------------------------------------------
  SUBROUTINE GlobalToLocalCoords(Element,Plane,n,Line,Eps, &
      Inside,Weights,maxind,linepos,LocalCoords)

    TYPE(Nodes_t) :: Plane, Line
    TYPE(Element_t), POINTER   :: Element
    INTEGER :: n, maxind
    REAL (KIND=dp) :: Eps,Weights(:)
    LOGICAL :: Inside
    REAL(KIND=dp) :: linepos
    REAL(KIND=dp) :: LocalCoords(3)
    
    REAL (KIND=dp) :: A(3,3),A0(3,3),B(3),C(3),Eps2,detA,absA,ds
    INTEGER :: split, i, corners, visited=0
    REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3)
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W=0.0d0

    SAVE visited
    visited = visited + 1

    Inside = .FALSE.
    corners = MIN(n,4)
    LocalCoords = 0.0_dp
    
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

    ELSE IF( DIM == 3 ) THEN
      ! In 3D rectangular faces are treated as two triangles
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

        B(1) = Plane % x(1) - Line % x(1)
        B(2) = Plane % y(1) - Line % y(1)
        B(3) = Plane % z(1) - Line % z(1)
        
        CALL InvertMatrix( A,3 )
        C(1:3) = MATMUL( A(1:3,1:3),B(1:3) )

        ! position on the line
        linepos = C(1)
        
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
    
    Weights(1:n) = Basis(1:n)
    MaxInd = 1
    DO i=2,n
      IF(Weights(MaxInd) < Weights(i)) MaxInd = i
    END DO

    LocalCoords(1) = u
    LocalCoords(2) = v
    LocalCoords(3) = w
    
  END SUBROUTINE GlobalToLocalCoords
  


!---------------------------------------------------------------------------
!> As the previous but performs mapping to plane
!---------------------------------------------------------------------------
  SUBROUTINE GlobalToLocalCoordsReduced(Element,Plane,n,Line,Eps, &
      Inside,Weights,maxind,linepos,LocalCoords)

    TYPE(Nodes_t) :: Plane, Line
    TYPE(Element_t), POINTER   :: Element
    INTEGER :: n, maxind
    REAL (KIND=dp) :: Eps,Weights(:)
    LOGICAL :: Inside
    REAL(KIND=dp) :: linepos
    REAL(KIND=dp) :: LocalCoords(3)

    REAL (KIND=dp) :: A(3,3),A0(3,3),B(3),C(3),Eps2,detA,absA,ds
    INTEGER :: split, i, corners, visited=0
    REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3)
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W=0.0d0

    SAVE visited
    visited = visited + 1

    corners = MIN(n,4)
    LocalCoords = 0.0_dp
    
    Eps2 = SQRT(TINY(Eps2))    

    IF( DIM /= 3 ) CALL Fatal(Caller,'GlobalToLocalReduced only for 3 dimensions!')

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
      CALL Fatal('GlobalToLocalCoords','Impossible value for parameter IntersectCoordinate: '&
          //I2S(IntersectCoordinate))
    END IF
    
    A0 = A
    
    detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
    absA = SUM(ABS(A(1,1:2))) * SUM(ABS(A(2,1:2)))
    
    IF(ABS(detA) <= eps * absA + Eps2) RETURN
    
    B(1) = Plane % x(1) - Line % x(1) 
    B(2) = Plane % y(1) - Line % y(1) 
    
    CALL InvertMatrix( A,2 )
    C(1:2) = MATMUL(A(1:2,1:2),B(1:2))
    linepos = C(1)
    
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

    ! Ok, we are this far so we must be inside    
    u = -1.0d0 + 2.0d0 * C(2)

    stat = ElementInfo( Element, Plane, U, V, W, SqrtElementMetric, &
        Basis, dBasisdx )
    
    Weights(1:n) = Basis(1:n)
    MaxInd = 1
    DO i=2,n
      IF(Weights(MaxInd) < Weights(i)) MaxInd = i
    END DO

    LocalCoords(1) = u
    LocalCoords(2) = v
    LocalCoords(3) = w
    
  END SUBROUTINE GlobalToLocalCoordsReduced
  

  SUBROUTINE OpenLineFile(ParallelReduce) 
    LOGICAL, OPTIONAL :: ParallelReduce    
    LOGICAL :: Parallel
    INTEGER :: iostat
    
    IF(FileIsOpen) RETURN
    FileIsOpen = .TRUE.    
    
    SideFile = TRIM(OutputDirectory)// '/' //TRIM(SideFile)
    
    Parallel = ( ParEnv % PEs > 1 )
    IF( Parallel ) THEN
      IF( PRESENT( ParallelReduce ) ) Parallel = .NOT. ParallelReduce
    END IF

    IF( Parallel ) THEN
      SideParFile = TRIM(SideFile)//'.'//I2S(ParEnv % MyPe)
    ELSE
      SideParFile = TRIM(SideFile)
    END IF
    
    IF(ListGetLogical(Params,'Filename Numbering',GotIt)) THEN
      IF( Parallel ) THEN
        CALL Warn(Caller,'Cannot number filenames in parallel with another number!')
      ELSE
        SideParFile = NextFreeFilename( SideParFile ) 
      END IF
    END IF

    CALL Info(Caller,'Saving line data to file: '//TRIM(SideParFile),Level=12)

    FileAppend = ListGetLogical(Params,'File Append',GotIt )

    IF( Solver % TimesVisited > 0 .OR. FileAppend) THEN 
      OPEN (NEWUNIT=LineUnit, FILE=SideParFile,POSITION='APPEND',iostat=iostat)
    ELSE
      OPEN (NEWUNIT=LineUnit,FILE=SideParFile,iostat=iostat)      
    END IF
    IF( iostat /= 0 ) CALL Warn(Caller,'Problems closing line file: '//I2S(iostat))

  END SUBROUTINE OpenLineFile

  
  SUBROUTINE CloseLineFile()

    INTEGER, POINTER :: NewOrder(:)
    LOGICAL :: ParallelReduce 
    INTEGER :: iostat
    REAL(KIND=dp) :: dpos
    
    ParallelReduce = .FALSE.
    IF( ParEnv % PEs > 1 ) THEN
      ParallelReduce = ListGetLogical(Params,'Parallel Reduce',Found )
    END IF
        
    IF( ParallelReduce ) THEN
      BLOCK 
        INTEGER :: SavePart, MaxSize, TotSize, NoPart
        REAL(KIND=dp), POINTER :: tmpPosData(:), tmpResultData(:,:)
        INTEGER, POINTER :: tmpLabelData(:,:),recsize(:)
        INTEGER :: nold, nnew, offset
        INTEGER :: ierr,status(MPI_STATUS_SIZE)

        NULLIFY(tmpPosData,tmpResultData,tmpLabelData,recsize)
        
        NoPart = 0
        IF( NoData > 0 ) NoPart = 1
        NoPart = ParallelReduction(NoPart)
        CALL Info(Caller,'Data available in number of partitions: '//I2S(NoPart))

        MaxSize = ParallelReduction(NoData,2)
        CALL Info(Caller,'Maximum number of lines in partition: '//I2S(MaxSize))

        TotSize = ParallelReduction(NoData)
        CALL Info(Caller,'Total number of lines in all partitions: '//I2S(TotSize))

        k = -1
        IF(MaxSize == NoData) k = ParEnv % MyPe
        SavePart = ParallelReduction(k,2)
        CALL Info(Caller,'Partition chosen for saving the data: '//I2S(SavePart))        
        
        ! Ok, we have data in several partitions. Bring it all to partition "SavePart". 
        IF( NoPart > 1 ) THEN

          ! Grow the data to facilate also parallel stuff.
          ! We know that this partition has >0 data already and that there is data coming
          ! from other partitions as well. 
          IF( ParEnv % MyPe == SavePart ) THEN
            nnew = TotSize
            ALLOCATE(tmpPosData(nnew),tmpLabelData(nnew,NoLabels),&
                tmpResultData(nnew,NoResults),STAT=istat)
            IF( istat /= 0 ) THEN
              CALL Fatal(Caller,'Problems allocating temporal workspace for parallel reduction')
            END IF
            tmpPosData = 0.0_dp
            tmpLabelData = 0
            tmpResultData = 0.0_dp
            
            tmpPosData(1:NoData) = PosData(1:NoData)              
            tmpLabelData(1:NoData,1:NoLabels) = LabelData(1:NoData,1:NoLabels)
            tmpResultData(1:NoData,1:NoResults) = ResultData(1:NoData,1:NoResults)

            DEALLOCATE(PosData, ResultData, LabelData, STAT=istat)
            IF(istat /= 0) CALL Fatal(Caller,'Problems deallocating some workspace')
            
            PosData => tmpPosData
            ResultData => tmpResultData
            LabelData => tmpLabelData            
          END IF
          
          ! Sent data sizes:
          !--------------------------
          IF( ParEnv % Mype == SavePart ) THEN
            ALLOCATE(RecSize(0:ParEnv % PEs-1))
            RecSize = 0
            DO i=0, ParEnv % PEs-1
              IF( i == SavePart) CYCLE
              CALL MPI_RECV( RecSize(i), 1, MPI_INTEGER, i, &
                  1000, ELMER_COMM_WORLD, status, ierr )
            END DO
          ELSE            
            CALL MPI_BSEND( NoData, 1, MPI_INTEGER, SavePart, &
                1000, ELMER_COMM_WORLD, ierr )
          END IF
          CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )

          ! Sent actual data
          !--------------------------
          IF( ParEnv % Mype == SavePart ) THEN                                   
            offset = NoData
            !PRINT *,'Recieving with offset:',ParEnv % MyPe, offset
            
            DO i=0, ParEnv % PEs-1
              IF( i == SavePart) CYCLE
              j = RecSize(i)
              IF(j==0) CYCLE
              
              NULLIFY(tmpPosData, tmpLabelData, tmpResultData)
              ALLOCATE(tmpPosData(j),tmpLabelData(j,NoLabels),tmpResultData(j,NoResults),STAT=istat)              
              IF( istat /= 0 ) THEN
                CALL Fatal(Caller,'Problems allocating temporal workspace for parallel communication')
              END IF              
              
              CALL MPI_RECV( tmpPosData, j, MPI_DOUBLE_PRECISION, i, &
                  1001, ELMER_COMM_WORLD, status, ierr )
              CALL MPI_RECV( tmpLabelData, j*NoLabels, MPI_INTEGER, i, &
                  1002, ELMER_COMM_WORLD, status, ierr )
              CALL MPI_RECV( tmpResultData, j*NoResults, MPI_DOUBLE_PRECISION, i, &
                  1003, ELMER_COMM_WORLD, status, ierr )
              
              PosData(Offset+1:offset+j) = tmpPosData(1:j)              
              LabelData(Offset+1:offset+j,1:NoLabels) = tmpLabelData(1:j,1:NoLabels)
              ResultData(Offset+1:offset+j,1:NoResults) = tmpResultData(1:j,1:NoResults)
              
              DEALLOCATE(tmpPosData, tmpLabelData, tmpResultData,STAT=istat)
              IF( istat /= 0 ) THEN
                CALL Fatal(Caller,'Problems deallocating temporal workspace for parallel communication')
              END IF              
              offset = offset + j

              !PRINT *,'offset:',offset
            END DO
          ELSE IF( NoData > 0 ) THEN
            NULLIFY(tmpPosData, tmpLabelData, tmpResultData)
            j = NoData
            ! Most likely the current tables are too big. In order for successful parallel communication
            ! the matrices must have exactly the same size!
            ALLOCATE(tmpPosData(j),tmpLabelData(j,NoLabels),tmpResultData(j,NoResults),STAT=istat)              
            IF( istat /= 0 ) THEN
              CALL Fatal(Caller,'Problems allocating temporal workspace for parallel communication')
            END IF
            tmpPosData(1:j) = PosData(1:j)               
            tmpLabelData(1:j,1:NoLabels) = LabelData(1:j,1:NoLabels)
            tmpResultData(1:j,1:NoResults) = ResultData(1:j,1:NoResults)
            
            CALL MPI_BSEND( tmpPosData, j, MPI_DOUBLE_PRECISION, SavePart, &
                1001, ELMER_COMM_WORLD, ierr )
            CALL MPI_BSEND( tmpLabelData, j*NoLabels, MPI_INTEGER, SavePart, &
                1002, ELMER_COMM_WORLD, ierr )
            CALL MPI_BSEND( tmpResultData, j*NoResults, MPI_DOUBLE_PRECISION, SavePart, &
                1003, ELMER_COMM_WORLD, ierr )
            DEALLOCATE(tmpPosData, tmpLabelData, tmpResultData,STAT=istat)
            IF( istat /= 0 ) THEN
              CALL Fatal(Caller,'Problems deallocating temporal workspace for parallel communication')
            END IF            
          END IF

          
          CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )
           
          IF( ParEnv % MyPe == SavePart ) THEN
            DEALLOCATE(recsize,STAT=istat)              
            IF(istat /= 0) THEN
              CALL Fatal(Caller,'Problems deallocating recsize vector!')
            END IF
            NoData = TotSize
          ELSE IF( NoData > 0 ) THEN
            DEALLOCATE(PosData,LabelData,ResultData,STAT=istat)                            
            IF(istat /= 0) THEN
              CALL Fatal(Caller,'Problems deallocating some unneeded workspace')
            END IF
            NoData = 0
          END IF
        END IF
      END BLOCK
    END IF
    
    
    IF( NoData > 0 ) THEN
      ! If we have a meaning full position info:
      NewOrder => NULL()
      IF( PosData(NoData) > 0.0_dp ) THEN
        ALLOCATE(NewOrder(NoData),STAT=istat)
        IF(istat /= 0) CALL Fatal(Caller,'Problems allocating NewOrder vector!')
        CALL Info(Caller,'Sorting and saving '//I2S(NoData)//' tabulated rows',Level=7)        
        DO i=1,NoData
          NewOrder(i) = i
        END DO
        CALL SortR( NoData,NewOrder,PosData)
      END IF
        
      CALL OpenLineFile(ParallelReduce)
                
      DO i = 1, NoData
        k = i
        IF( ASSOCIATED( NewOrder ) ) THEN
          ! Go from small coordinates to bigger...
          k = NewOrder(NoData-i+1)
        END IF
        IF(i>1 .AND. PosData(i) > 0.0_dp ) THEN
          dpos = ABS(PosData(i)-PosData(i-1))
          IF( dpos < EPSILON(dpos) ) THEN 
            !PRINT *,'skipping value',PosData(i),PosData(i-1)
            CYCLE
          END IF
        END IF
        DO j=1,NoLabels
          WRITE(LineUnit,'(A)',ADVANCE='NO') I2S(LabelData(k,j))//' '
        END DO
        DO j=1,NoResults-1
          WRITE(LineUnit,'(ES20.11E3)',ADVANCE='NO') ResultData(k,j)
        END DO
        WRITE(LineUnit,'(ES20.11E3)') ResultData(k,NoResults)
      END DO
      NoData = 0
      
      IF( ASSOCIATED(NewOrder) ) THEN
        DEALLOCATE( NewOrder, STAT=istat )
        IF(istat /= 0) CALL Fatal(Caller,'Problems deallocating NewOrder vector')
      END IF
      DEALLOCATE( PosData,LabelData,ResultData,STAT=istat)
      IF(istat /= 0) CALL Fatal(Caller,'Problems deallocating some temporal workspace')
    END IF
    
    IF(FileIsOpen) THEN
      CLOSE(LineUnit,iostat=iostat)
      IF( iostat /= 0 ) CALL Fatal(Caller,'Problems closing line file: '//I2S(iostat))
    END IF
      
  END SUBROUTINE CloseLineFile


  ! Write a line of data for a point in a known element.
  !----------------------------------------------------------------------
  SUBROUTINE WriteFieldsAtElement( Element, BC_id, &
      node_id, dgnode_id, Basis, UseNode, NodalFlux, LocalCoord, &
      GlobalCoord, linepos, ParNode )

    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), TARGET, OPTIONAL :: Basis(:)
    INTEGER :: bc_id, node_id, dgnode_id
    
    LOGICAL, OPTIONAL :: UseNode 
    REAL(KIND=dp), OPTIONAL :: NodalFlux(3)
    REAL(KIND=dp), OPTIONAL :: LocalCoord(3)
    REAL(KIND=dp), OPTIONAL :: GlobalCoord(3)
    REAL(KIND=dp), OPTIONAL :: linepos
    LOGICAL, OPTIONAL :: ParNode
    
    INTEGER :: i,j,k,l,ivar,ii,i1,i2,n0
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: UseGivenNode, PiolaVersion, EdgeBasis
    INTEGER :: n, nd, np, EdgeBasisDegree, Labels(5)
    INTEGER, POINTER :: PtoIndexes(:)
    REAL(KIND=dp), POINTER :: PtoBasis(:)
    REAL(KIND=dp), TARGET :: PointBasis(1)
    REAL(KIND=dp) :: u,v,w
    INTEGER, TARGET :: NodeIndex(1), Indexes(54), DGIndexes(27)
    REAL(KIND=dp) :: up,vp,wp
    REAL(KIND=dp), TARGET :: NodeBasis(54)
    REAL(KIND=dp) :: WBasis(54,3),RotWBasis(54,3), NodedBasisdx(54,3)
    REAL(KIND=dp) :: AveMult
    REAL(KIND=dp), ALLOCATABLE, SAVE :: fdg(:), fip(:)
    LOGICAL :: pElem, Tabulate
    TYPE(Variable_t), POINTER :: pVar
    
    SAVE :: Nodes

    DGIndexes = 0
    Indexes = 0
    n0 = 0
    
    ! If we have position on the line then we can sort the entries.
    Tabulate = .FALSE.
    IF(PRESENT(linepos)) THEN
      Tabulate = ( linepos > 0 ) .OR. Parallel
    END IF
    
    IF( .NOT. SkipBoundaryInfo ) THEN      
      Labels = 0
      IF( LineInd /= 0 ) THEN
        n0 = n0 + 1
        Labels(n0) = LineInd
      END IF

      IF(TransientSimulation) THEN
        n0 = n0 + 1
        Labels(n0) = Solver % DoneTime
      END IF

      Labels(n0+1) = Solver % TimesVisited + 1
      Labels(n0+2) = bc_id
      k = node_id
      IF(PRESENT(ParNode)) THEN
        IF( ParEnv % PEs > 1 .AND. .NOT. Mesh % SingleMesh ) THEN
          k = Mesh % ParallelInfo % GlobalDOFS(k)
        END IF
      END IF
      Labels(n0+3) = k      
      n0 = n0 + 3
      
      IF( NormInd > 0 .AND. NormInd <= n0 ) THEN
        Norm = Norm + 1.0_dp * Labels(NormInd )
      END IF
    END IF

    ! If we just got the closest node, not really the exact node, do not use it.
    UseGivenNode = .FALSE.
    IF( PRESENT( UseNode ) ) UseGivenNode = UseNode
    
    No = 0
    Values = 0.0d0

    ! The funny negative indexes refer to coordinates that are treated separately
    DO ivar = -2,NoVar
      Var => VariableGetN( ivar, comps ) 
      IF( comps >= 2 ) THEN
        Var2 => VariableGetN( ivar, component = 2 ) 
      ELSE
        Var2 => NULL()
      END IF
      IF( comps >= 3 ) THEN
        Var3 => VariableGetN( ivar, component = 3 ) 
      ELSE
        Var3 => NULL()
      END IF

      IF( UseGivenNode ) THEN
        k = node_id; l = dgnode_id
      ELSE
        k = 0; l = 0
      END IF
      ! This is additive adding to "No" for each call!!
      CALL EvaluateVariableAtGivenPoint(No,Values,Mesh,Var,Var2,Var3,Element,LocalCoord,&
          Basis,k,l,GotEigen=DoEigen,GotEdge=EdgeBasis)
    END DO
    
    IF( CalculateFlux ) THEN
      IF( PRESENT( NodalFlux ) ) THEN
        Values(No+1:No+3) = NodalFlux
      ELSE
        Values(No+1:No+3) = 0.0
      END IF
    END IF

    DO j=1,NoResults
      IF( ISNAN(Values(j)) ) THEN
        FoundNan = FoundNan + 1
        Values(j) = -1.0_dp
      END IF
    END DO
    
    IF( Tabulate ) THEN
      NoLabels = MAX(n0, NoLabels) 
      NoData = NoData + 1
      
      BLOCK
        REAL(KIND=dp), POINTER :: tmpPosData(:), tmpResultData(:,:)
        INTEGER, POINTER :: tmpLabelData(:,:)
        INTEGER :: nold, nnew

        nold = 0
        IF(ASSOCIATED(PosData)) nold = SIZE(PosData)
        
        IF( NoData > nold ) THEN
          nnew = MAX(100, 2*nold)
          CALL Info(Caller,'Increasing temporal size from '//I2S(nold)//' to '//I2S(nnew),Level=7)
          ALLOCATE(tmpPosData(nnew),tmpLabelData(nnew,NoLabels),tmpResultData(nnew,NoResults))
          tmpPosData = 0.0_dp
          tmpLabelData = 0
          tmpResultData = 0.0_dp
          
          IF( nold > 0 ) THEN
            tmpPosData = PosData
            tmpResultData = ResultData
            tmpLabelData = LabelData
            DEALLOCATE(PosData, ResultData, LabelData,STAT=istat)
            IF(istat /= 0) THEN
              CALL Fatal(Caller,'Problems deallocating some too small workspace')
            END IF 
          END IF
          PosData => tmpPosData
          ResultData => tmpResultData
          LabelData => tmpLabelData
        END IF
      END BLOCK
        
      PosData(NoData) = linepos
      LabelData(NoData,1:NoLabels) = Labels(1:NoLabels) 
      ResultData(NoData,1:NoResults) = Values(1:NoResults)
    ELSE    
      CALL OpenLineFile()      
      DO i=1,n0
        WRITE(LineUnit,'(A)',ADVANCE='NO') I2S(Labels(i))//' '
      END DO
      DO j=1,NoResults-1
        WRITE(LineUnit,'(ES20.11E3)',ADVANCE='NO') Values(j)
      END DO
      WRITE(LineUnit,'(ES20.11E3)') Values(NoResults)
    END IF
          
    IF( NormInd > n0 ) THEN
      Norm = Norm + Values(NormInd-n0)
    END IF

    
  END SUBROUTINE WriteFieldsAtElement
    



!-----------------------------------------------------------------------
!> Computation of normal flux.
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
    REAL(KIND=dp) :: r,u,v,w,ub,DetJ, Normal(3),Flow(3)
    LOGICAL :: stat, Permutated
    INTEGER :: body_id, k
    REAL(KIND=DP), POINTER :: Pwrk(:,:,:) => Null()
    TYPE(ValueList_t), POINTER :: BC
    INTEGER :: FluxBody
    
    
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), Conductivity(:), &
        CoeffTensor(:,:,:)
    LOGICAL :: AllocationsDone = .FALSE.

    SAVE AllocationsDone, Nodes, Basis, dBasisdx, Conductivity, CoeffTensor, Pwrk
        
    IF( .NOT. AllocationsDone ) THEN
      n = Mesh % MaxElementNodes
      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n), Basis(n), dBasisdx(n,3), &
          Conductivity(n), CoeffTensor(3,3,n) )
      AllocationsDone = .TRUE.
    END IF    
    
    Tvar => VariableGet( Mesh % Variables, TRIM(VarName) )
    IF( .NOT. ASSOCIATED( TVar ) ) THEN
      CALL Fatal('BoundaryFlux','Cannot calculate fluxes without potential field!')
    END IF
    
    Permutated = ASSOCIATED(Tvar % Perm)
    Element => Model % CurrentElement

    BC => GetBC( Element )
    FluxBody = ListGetInteger( BC,'Flux Integrate Body', gotIt )
    
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
        CALL Fatal( Caller, Message )
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
      IF ( .NOT. stat )  CALL Fatal( Caller,&
          'No solution available for specified boundary' )
    END IF
    
    n = Parent % TYPE % NumberOfNodes

    Nodes % x(1:n) = Mesh % Nodes % x(Parent % NodeIndexes)
    Nodes % y(1:n) = Mesh % Nodes % y(Parent % NodeIndexes)
    Nodes % z(1:n) = Mesh % Nodes % z(Parent % NodeIndexes)

    k = 0
    DO j=1,n
      IF ( node == Parent % NodeIndexes(j) ) THEN
        k = j
        EXIT
      END IF
    END DO

    IF ( k == 0 ) THEN
      CALL Warn(Caller,'Side node not in parent element!')
    END IF
    
    CALL GlobalToLocal( u, v ,w , Nodes % x(k), Nodes % y(k), Nodes % z(k), Parent, Nodes )

    stat = ElementInfo( Parent, Nodes, u, v, w, detJ, Basis, dBasisdx )
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
    IF( .NOT. ASSOCIATED( Material ) ) THEN
      CALL Warn(Caller,'Could not find material for flux computation!')
      RETURN
    END IF
    
    OldCurrentElement => Model % CurrentElement
    Model % CurrentElement => Parent    
    CALL ListGetRealArray( Material, TRIM(CoeffName), Pwrk, n, &
        Parent % NodeIndexes, GotIt )
    Model % CurrentElement => OldCurrentElement
      
    IF(.NOT. ASSOCIATED( Pwrk ) ) THEN
      CALL Warn(Caller,'Coefficient not present for flux computation!')
      RETURN
    END IF

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




! Save a line (or boundary) that exist already in mesh.
! Data is thus saved in existing nodes. 
!-----------------------------------------------------------------------
  SUBROUTINE SaveExistingLines()

    INTEGER :: dgnode 
    REAL(KIND=dp) :: Coord(3), Coord0(3), Center(3)
    TYPE(ValueList_t), POINTER :: ValueList
    TYPE(Element_t), POINTER :: Parent
    LOGICAL :: BreakLoop, ParallelComm
    REAL(KIND=dp) :: linepos
    
    MaskName = ListGetString(Params,'Save Mask',GotIt) 
    IF(.NOT. GotIt) MaskName = 'Save Line'

    IF( .NOT. ( ListCheckPresentAnyBC( Model, MaskName ) .OR. &
        ListCheckPresentAnyBodyForce( Model, MaskName ) ) ) RETURN

    CALL Info(Caller,'Saving existing nodes into ascii table',Level=8)

    IF( Solver % TimesVisited > 0 ) THEN
      InitializePerm = ( MaskName /= PrevMaskName ) 
      InitializePerm = InitializePerm .OR. Solver % MeshChanged
    ELSE
      InitializePerm = .TRUE.
    END IF

    IF( InitializePerm ) THEN
      SavePerm = 0

      OptimizeOrder = ListGetLogical(Params,'Optimize Node Ordering',GotIt)
      IF(.NOT. GotIt) OptimizeOrder = .NOT. Parallel

      BreakLoop = ListGetLogical(Params,'Break Line Loop',GotIt)
      IF(BreakLoop) OptimizeOrder = .TRUE.

      ParallelComm = Parallel
      IF(Mesh % SingleMesh ) ParallelComm = .FALSE.
      CALL MakePermUsingMask( Model,Solver,Mesh,MaskName, &
          OptimizeOrder, SavePerm, SaveNodes(1), &
          RequireLogical = .TRUE., BreakLoop = BreakLoop, ParallelComm = ParallelComm )
      
      IF( SaveNodes(1) > 0 ) THEN
        IF( ListGetLogical( Params,'Calculate Weights',GotIt ) ) THEN
          CALL CalculateNodalWeights( Solver, .TRUE., SavePerm, TRIM(MaskName)//' Weights')
        END IF
        CALL Info(Caller,'Number of nodes in specified boundary: '//I2S(SaveNodes(1)),Level=12)
      END IF
    ELSE
      SaveNodes(1) = 0
      IF( ASSOCIATED( SavePerm ) ) SaveNodes(1) = MAXVAL( SavePerm ) 
    END IF
    PrevMaskName = MaskName

    
    !------------------------------------------------------------------------------
    ! If nodes found, then go through the sides and compute the fluxes if requested 
    !------------------------------------------------------------------------------
    IF( SaveNodes(1) > 0 ) THEN

      ALLOCATE( InvPerm(SaveNodes(1)), BoundaryIndex(SaveNodes(1)), STAT=istat )
      IF( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 3: '//I2S(SaveNodes(1))) 
      
      BoundaryIndex = 0
      InvPerm = 0
      DO i=1,SIZE(SavePerm)
        IF (SavePerm(i)>0) THEN
          ! Error check for something that should never happen
          IF( InvPerm( SavePerm(i)) > 0) THEN
            WRITE( Message, *) 'Node multiple times in permutation',i,SavePerm(i)
            CALL Warn(Caller,Message)
          END IF
          InvPerm(SavePerm(i)) = i
        END IF
      END DO
      
      IF(CalculateFlux) THEN
        CALL Info(Caller,'Calculating nodal fluxes',Level=8)
        ALLOCATE(PointFluxes(SaveNodes(1),3),PointWeight(SaveNodes(1)), STAT=istat)    

        IF( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 4') 

        PointFluxes = 0.0d0
        PointWeight = 0.0d0

        ! Go through the elements and register the boundary index and fluxes if asked
        ! Fluxes only possible for DIM-1 
        DO t = Mesh % NumberOfBulkElements + 1,  &
            Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements                        
                    
          CurrentElement => Mesh % Elements(t)
          IF( ParEnv % PEs > 1 ) THEN
            IF( CurrentElement % PartIndex /= ParEnv % MyPe ) CYCLE
          END IF

          Model % CurrentElement => CurrentElement
          n = CurrentElement % TYPE % NumberOfNodes
          NodeIndexes => CurrentElement % NodeIndexes        
          
          IF( .NOT. ALL(SavePerm(NodeIndexes) > 0)) CYCLE 

          IF(t <= Mesh % NumberOfBulkElements) THEN
            ValueList => GetBodyForce()
          ELSE
            ValueList => GetBC()
          END IF
 
          IF( .NOT. ASSOCIATED( ValueList ) ) CYCLE
          IF( .NOT. ListCheckPresent( ValueList, MaskName ) ) CYCLE
                  
          ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
          ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

          DO i=1,n
            node = NodeIndexes(i)

            CALL BoundaryFlux( Model, node, TempName,  &
                CondName, f1, f2, fn, weight, Mesh % MaxElementDOFs ) 
            
            j = SavePerm(node) 
            IF( j == 0 ) CYCLE
            
            PointFluxes(j,1) = PointFluxes(j,1) + weight * f1
            PointFluxes(j,2) = PointFluxes(j,2) + weight * f2
            PointFluxes(j,3) = PointFluxes(j,3) + weight * fn
            PointWeight(j) = PointWeight(j) + weight
          END DO
        END DO
        
        ! Normalize flux by division with the integration weight
        DO i = 1, SaveNodes(1)
          PointFluxes(i,1) = PointFluxes(i,1) / PointWeight(i)
          PointFluxes(i,2) = PointFluxes(i,2) / PointWeight(i)
          PointFluxes(i,3) = PointFluxes(i,3) / PointWeight(i)
        END DO
      END IF

      !DgVar = ASSOCIATED( Mesh % Elements(1) % DGIndexes ) 
      
      ! Go through the elements and register the boundary index and fluxes if asked
      DO t = 1,  Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements        
        
        CurrentElement => Mesh % Elements(t)
        IF( ParEnv % PEs > 1 ) THEN
          IF( CurrentElement % PartIndex /= ParEnv % MyPe ) CYCLE
        END IF

        Model % CurrentElement => CurrentElement
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes        
        
        IF( .NOT. ALL(SavePerm(NodeIndexes) > 0)) CYCLE 

        IF(t > Mesh % NumberOfBulkElements) THEN
          k = GetBCId( CurrentElement )
          IF( k == 0 ) CYCLE
          ValueList => Model % BCs(k) % Values
        ELSE
          k = GetBodyForceId( CurrentElement )         
          IF( k == 0 ) CYCLE
          ValueList => Model % BodyForces(k) % Values
        END IF
        
        IF( .NOT. ListCheckPresent( ValueList, MaskName ) ) CYCLE
                        
        IF( DGVar ) THEN
          ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
          ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

          Center(1) = SUM( ElementNodes % x(1:n) ) / n
          Center(2) = SUM( ElementNodes % y(1:n) ) / n
          Center(3) = SUM( ElementNodes % z(1:n) ) / n
          
          DO i = 1, n
            node = NodeIndexes(i) 

            IF( t > Mesh % NumberOfBulkElements ) THEN
              Found = .FALSE.
              Parent => CurrentElement % BoundaryInfo % Left 
              IF( .NOT. ASSOCIATED( Parent ) ) THEN
                CALL Fatal(Caller,'Parent not associated!')
              END IF
              DO j = 1, SIZE( Parent % NodeIndexes ) 
                IF( node == Parent % NodeIndexes(j) ) THEN
                  dgnode = Parent % DgIndexes(j)
                  Found = .TRUE.
                  EXIT
                END IF
              END DO
              IF(.NOT. Found) CALL Fatal(Caller,'Could not find DG node!')              
            END IF
                        
            Coord(1) = ElementNodes % x(i)
            Coord(2) = ElementNodes % y(i)
            Coord(3) = ElementNodes % z(i)

            ! Shrink the element so that external sort work better!
            Coord0 = Coord
            Coord = Center + 0.9999*(Coord-Center)

            ! Do this dirty way such that DG nodes may be sorted
            Mesh % Nodes % x(node) = Coord(1) 
            Mesh % Nodes % y(node) = Coord(2) 
            IF( dim == 3 ) Mesh % Nodes % z(node) = Coord(3) 

            linepos = -1.0_dp
            IF( ParEnv % PEs > 1 ) THEN
            END IF

            IF( CalculateFlux ) THEN
              CALL WriteFieldsAtElement( CurrentElement, k, node, &
                  dgnode, UseNode = .TRUE., NodalFlux = PointFluxes(t,:), &
                  linepos = linepos, ParNode = Parallel )
            ELSE
              CALL WriteFieldsAtElement( CurrentElement, k, node, &
                  dgnode, UseNode = .TRUE., linepos = linepos, ParNode = Parallel )
            END IF
            
            ! and revert 
            Mesh % Nodes % x(node) = Coord0(1) 
            Mesh % Nodes % y(node) = Coord0(2) 
            IF( dim == 3 ) Mesh % Nodes % z(node) = Coord0(3)                       
          END DO
        ELSE
          BoundaryIndex( SavePerm(NodeIndexes) ) = k
        END IF

        MaxBoundary = MAX( MaxBoundary, k ) 

      END DO
      
      ! Save the nodes if not in DG mode     
      !---------------------------------  
      IF( .NOT. DGVar ) THEN
        dgnode = 0
        linepos = -1.0_dp
        DO t = 1, SaveNodes(1)    
          node = InvPerm(t)
          IF( CalculateFlux ) THEN
            CALL WriteFieldsAtElement( CurrentElement, BoundaryIndex(t), node, &
                dgnode, UseNode = .TRUE., NodalFlux = PointFluxes(t,:), &
                linepos = linepos, ParNode = Parallel )
          ELSE
            CALL WriteFieldsAtElement( CurrentElement, BoundaryIndex(t), node, &
                dgnode, UseNode = .TRUE., linepos = linepos, ParNode = Parallel )
          END IF
        END DO
      END IF        
      
      DEALLOCATE(InvPerm, BoundaryIndex,STAT=istat)
      IF( istat /= 0 ) CALL Fatal(Caller,'Problems deallocating stuff on existing lines')
            
      IF(CalculateFlux) THEN
        DEALLOCATE(PointFluxes, PointWeight, STAT=istat)
        IF( istat /= 0 ) CALL Fatal(Caller,'Problems deallocating fluxes on existing lines')
      END IF
    END IF
    
  END SUBROUTINE SaveExistingLines



  ! Save data on given polylines. These are created on-the-fly.
  ! Data is written either on intersections with element faces (edges), or
  ! in uniformly distributed points. 
  !-------------------------------------------------------------------------------------
  SUBROUTINE SavePolyLines()

    TYPE(Solver_t), POINTER :: pSolver
    REAL(KIND=dp) :: linepos = 0, tanprod(2), s, eps
    
    pSolver => Solver
    eps = 1.0e-5
    
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
        CALL Warn(Caller,'Dimension of points smaller than that of mesh')
      END IF
    ELSE 
      NoLines = 0
    END IF

    CALL Info(Caller,'Saving PolyLines into ascii table',Level=8)


    GotDivisions = .FALSE.
    IF( NoLines > 0 ) THEN
      NoDivisions => ListGetIntegerArray( Params,'Polyline Divisions',GotDivisions)
      IF(.NOT. GotDivisions) THEN
        NoDivisions => ListGetIntegerArray( Params,'Polyline Intervals',GotDivisions)
      END IF
      IF( GotDivisions ) THEN
        IF( SIZE( NoDivisions ) < NoLines + COUNT(SaveAxis) ) THEN
          CALL Fatal(Caller,'Polyline divisions size too small!')
        END IF
      END IF
    END IF

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
      IF( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for LineTag') 


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
          s = SQRT( SUM(dR**2) )
          LineN = dR / s
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
            IF( ParEnv % PEs > 1 ) THEN
              IF( CurrentElement % PartIndex /= ParEnv % MyPe ) CYCLE
            END IF

            n = CurrentElement % TYPE % NumberOfNodes
            m = GetElementCorners( CurrentElement )
            NodeIndexes => CurrentElement % NodeIndexes

            IF( MaskWithPerm ) THEN
              IF( ANY( Solver % Variable % Perm( NodeIndexes) == 0 ) ) CYCLE
            END IF
                        
            ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
            ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
            ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

            DO i=1,m
              S1(1) = ElementNodes % x(i)  
              S1(2) = ElementNodes % y(i)  
              S1(3) = ElementNodes % z(i)  

              dS = ( S1 - R0 ) / s
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

            DO i=1,2
              tanprod(i) = MinCoord(i+1) * MaxCoord(i+1)
            END DO

            IF( tanprod(1) > eps ) CYCLE
            IF( dim == 3 .AND. .NOT. IntersectEdge ) THEN
              IF( tanprod(2) > eps ) CYCLE
            END IF
            
            imin = MAX(0, CEILING( nsize * MinCoord(1) ) )
            imax = MIN(nsize, FLOOR( ( nsize * MaxCoord(1) ) ) )
            
            DO i=imin,imax
              NoTests = NoTests + 1

              IF( LineTag(i) ) CYCLE

              GlobalCoord = R0 + i * dR / nsize
              
              IF ( PointInElement( CurrentElement, ElementNodes, GlobalCoord, &
                  LocalCoord, USolver = pSolver, LocalEps = eps ) ) THEN
                stat = ElementInfo( CurrentElement, ElementNodes, LocalCoord(1), &
                LocalCoord(2), LocalCoord(3), detJ, Basis, USolver = pSolver )

                LineTag(i) = .TRUE.                
                SaveNodes(2) = SaveNodes(2) + 1               
                linepos = 1.0_dp*i/nsize + 2*(Line-1)

                n = CurrentElement % type % numberofnodes
                CALL WriteFieldsAtElement( CurrentElement, Line, i, 0, &
                    Basis, LocalCoord = LocalCoord, linepos = linepos )
              END IF
            END DO
          END DO
        ELSE
          ! If no divisions then go though existing faces and check for the
          ! intersection of line & and each face. 
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

            IF( IntersectCoordinate /= 0 ) THEN
              CALL GlobalToLocalCoordsReduced(CurrentElement,ElementNodes,n,LineNodes, &
                  DetEpsilon,Inside,Basis,i,linepos,LocalCoord)
            ELSE
              CALL GlobalToLocalCoords(CurrentElement,ElementNodes,n,LineNodes, &
                  DetEpsilon,Inside,Basis,i,linepos,LocalCoord)
            END IF

            IF(.NOT. Inside) CYCLE

            ! When the line goes through a node it might be saved several times 
            ! without this checking
            IF(1.0d0-MAXVAL(Basis(1:n)) < 1.0d-3) THEN
              IF( LineTag(NodeIndexes(i)) ) CYCLE
              LineTag(NodeIndexes(i)) = .TRUE.
            END IF

            SaveNodes(2) = SaveNodes(2) + 1
            
            linepos = linepos + 2*(Line-1)
            CALL WriteFieldsAtElement( CurrentElement, MaxBoundary, &
                NodeIndexes(i), 0, Basis, LocalCoord = LocalCoord, linepos = linepos )
          END DO
        END IF
      END DO

      IF( NoTests > 0 ) THEN
        CALL Info(Caller,'Number of candidate nodes: '//I2S(NoTests),Level=8)
      END IF

      CALL Info(Caller,'Number of nodes in specified lines: '//I2S(SaveNodes(2)))
      
      IF(ALLOCATED(LineTag)) DEALLOCATE( LineTag )
    END IF

  END SUBROUTINE SavePolyLines



  ! Save data on given circular lines.
  ! Data is saved in given number of divisions for each circle.
  !---------------------------------------------------------------------------
  SUBROUTINE SaveCircleLines()
    
    REAL(KIND=dp) :: CylCoord(3), Radius, Phi, Rtol
    TYPE(Solver_t), POINTER :: pSolver

    pSolver => Solver
    
    PointCoordinates => ListGetConstRealArray(Params,'Circle Coordinates',gotIt)
    IF(.NOT. GotIt) RETURN

    CALL Info(Caller,'Saving circular lines into ascii table',Level=8)

    NoLines = SIZE(PointCoordinates,1) 
    NoDims = SIZE(PointCoordinates,2)
    IF( NoDims /= 7 ) THEN
      CALL Fatal(Caller,'By construction the circle is defined by 7 values: '//I2S(NoDims))
    END IF

    NoDivisions => ListGetIntegerArray( Params,'Circle Divisions',GotIt)
    IF(.NOT. GotIt) THEN
      CALL Fatal(Caller,'> Circle Coordinates < requires > Circle Divisions < !')
    END IF
    IF( SIZE( NoDivisions ) < NoLines ) THEN
      CALL Fatal(Caller,'Polyline divisions size too small: '&
          //I2S(SIZE(NoDivisions)))
    END IF

    CALL Info(Caller,'Saving data on given circles: '//I2S(NoLines),Level=7)

    NoTests = 0
    
    t = MAXVAL( NoDivisions )     
    ALLOCATE( LineTag(t), STAT=istat )
    IF( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for LineTag') 

    DO Line = 1,NoLines 
      
      CALL Info(Caller,'Saving circle number: '//I2S(Line),Level=12)
      MaxBoundary = MaxBoundary + 1
      LineTag = .FALSE.
      
      nsize = NoDivisions(Line)

      ! Base point of cylinder
      R0 = PointCoordinates(Line,1:3)

      ! Normal of cylinder (aligned with axis)
      dR = PointCoordinates(Line,4:6)

      ! Radius of circle
      Radius = PointCoordinates(Line,7)

      PRINT *,'R0:',R0
      PRINT *,'dR:',Dr,Radius
      
      ! Normal and tangent directions
      LineN = dR / SQRT( SUM( dR**2 ) )
      CALL TangentDirections( LineN, LineT1, LineT2 ) 
      
      PRINT *,'LineN:',LineN
      PRINT *,'LineT1:',LineT1
      PRINT *,'LineT2:',LineT2
      PRINT *,'Elems:',Mesh % NumberOfBulkElements
      
      DO t = 1, Mesh % NumberOfBulkElements 
        CurrentElement => Mesh % Elements(t)
        IF( ParEnv % PEs > 1 ) THEN
          IF( CurrentElement % PartIndex /= ParEnv % MyPe ) CYCLE
        END IF

        n = CurrentElement % TYPE % NumberOfNodes
        m = GetElementCorners( CurrentElement )
        NodeIndexes => CurrentElement % NodeIndexes

        IF( MaskWithPerm ) THEN
          IF( ANY( Solver % Variable % Perm( NodeIndexes) == 0 ) ) CYCLE
        END IF
                
        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

        DO i=1,m
          S1(1) = ElementNodes % x(i)  
          S1(2) = ElementNodes % y(i)  
          S1(3) = ElementNodes % z(i)  

          S1 = S1 - R0

          ! Cartesian local coordinate system for elemental nodal
          LocalCoord(1) = SUM( S1 * LineT1 ) ! x
          LocalCoord(2) = SUM( S1 * LineT2 ) ! y 
          LocalCoord(3) = SUM( S1 * LineN )  ! z

          ! Cylindrical local coordinate system 
          CylCoord(1) = SQRT(LocalCoord(1)**2 + LocalCoord(2)**2) 
          CylCoord(2) = ATAN2(LocalCoord(2),LocalCoord(1)) / (2*PI)
          IF( CylCoord(2) < 0 ) CylCoord(2) = CylCoord(2) + 1.0_dp
          CylCoord(3) = LocalCoord(3)

          IF( i == 1 ) THEN
            MinCoord = CylCoord
            MaxCoord = CylCoord
          ELSE
            DO j=1,3
              MinCoord(j) = MIN( MinCoord(j), CylCoord(j) ) 
              MaxCoord(j) = MAX( MaxCoord(j), CylCoord(j) ) 
            END DO
          END IF
        END DO
          
        ! Element does not contain the correct z-level
        IF( MinCoord(3) * MaxCoord(3) > 0.0_dp ) CYCLE


        ! The circle in local cylindrical coordinates now goes from (0,0,0) to (0,1,0)

        ! Element does not contain the correct radius
        ! Here we use a tolerance comparable to the size of the element
        ! since in cylindrical coordinates the element does not have straight faces
        Rtol = 0.5*(MaxCoord(1) - MinCoord(1))
        IF( MinCoord(1) > Radius + Rtol ) CYCLE
        IF( MaxCoord(1) < Radius - Rtol ) CYCLE
       
        
        ! Ok, this the candidate interval within this element
        ! Should be in interval [-0.5,0.5]
        IF( MaxCoord(2) - MinCoord(2) > 0.5 ) THEN
          Phi = MinCoord(2)
          MinCoord(2) = MaxCoord(2)
          MaxCoord(2) = MinCoord(2) + 1.0_dp
        END IF

        imin = CEILING( nsize * MinCoord(2) - 1.0d-20)
        imax = FLOOR( nsize * MaxCoord(2) + 1.0d-20)

        DO i=imin,imax
          NoTests = NoTests + 1
          
          ii = i
          IF( ii <= 0 ) ii = ii + nsize
          IF( ii > nsize ) ii = ii - nsize

          IF( LineTag(ii) ) CYCLE
          
          Phi = ( ii * 2.0_dp * PI ) / nsize
          
          GlobalCoord = R0 + Radius * COS(Phi) * LineT1 + &
              Radius * SIN(Phi) * LineT2

          IF ( PointInElement( CurrentElement, ElementNodes, GlobalCoord, LocalCoord ) ) THEN
            stat = ElementInfo( CurrentElement, ElementNodes, LocalCoord(1), &
                LocalCoord(2), LocalCoord(3), detJ, Basis, USolver = pSolver )

            LineTag(ii) = .TRUE.
            SaveNodes(3) = SaveNodes(3) + 1
            CALL WriteFieldsAtElement( CurrentElement, Line, ii, &
                0, Basis, LocalCoord = LocalCoord, linepos = 20 + 2*Line + 1.0_dp *ii / nsize )
          END IF
        END DO
      END DO
      
      i = COUNT( LineTag(1:nsize) ) 
      PRINT *,'Points with hits:',i
            
      i = COUNT( .NOT. LineTag(1:nsize) ) 
      PRINT *,'Points with no hits:',i
    END DO


    IF( NoTests > 0 ) THEN
      CALL Info(Caller,'Number of candidate nodes: '//I2S(NoTests),Level=8)
    END IF
    
    CALL Info(Caller,'Number of nodes in specified circle: '//I2S(SaveNodes(3)))
    
    IF(ALLOCATED(LineTag)) DEALLOCATE( LineTag, STAT=istat)
    IF(istat /= 0) CALL Fatal(Caller,'Could not deallocate LineTag')

  END SUBROUTINE SaveCircleLines



  ! Save data on isocurves of some field.
  !--------------------------------------------------------------------------
  SUBROUTINE SaveIsoCurves()

    IF( .NOT. ListGetLogical( Params,'Save Isocurves',Found) ) RETURN

    CALL Info(Caller,'Saving isocurves into ascii table',Level=8)

    IF( DIM == 3 ) THEN
      CALL Fatal(Caller,'Isocurves can only be saved in 2D')
    END IF

    ALLOCATE( LineTag(0:Mesh % NumberOfNodes), STAT=istat )
    IF( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 6') 
        
    IF( Solver % TimesVisited == 0 ) THEN
      CALL FindMeshEdges( Mesh, .FALSE.)
    END IF
    
    NoFaces = Mesh % NumberOfEdges
    Line = 0

    DO WHILE( .TRUE. ) 
      
      Line = Line + 1
      LineTag = .FALSE.
 
      WRITE (Name,'(A,I0)') 'IsoSurface Variable ',Line
      VarName = ListGetString( Params, Name, GotVar )
      
      WRITE (Name,'(A,I0)') 'IsoSurface Value ',Line
      f0 = ListGetCReal( Params, Name, Found )

      IF( GotVar ) THEN
        IsosurfVar => VariableGet( Mesh % Variables, VarName )
        IF( .NOT. ASSOCIATED( IsosurfVar ) ) THEN
          CALL Warn(Caller,'Isosurface variable not given: '//TRIM(VarName))
          EXIT
        END IF
        IsosurfPerm => IsosurfVar % Perm
        Isosurf => IsosurfVar % Values       
      ELSE
        IF( Line == 1 ) THEN
          CALL Warn(Caller,'No > Isosurface Variable 1 < defined!')
          EXIT
        END IF
        IF(.NOT. Found ) EXIT
      END IF
      
      WRITE( Message, * ) 'Finding nodes on isocurve: ',Line
      CALL Info(Caller,Message)
      
      f1 = MINVAL( Isosurf ) 
      f2 = MAXVAL( Isosurf ) 
      IF( f0 <= f1 .OR. f0 >= f2 ) THEN
        CALL Warn(Caller,'Isosurface value not within range!')        
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

        IF( MaskWithPerm ) THEN
          IF( ANY( Solver % Variable % Perm( NodeIndexes(1:2) ) == 0 ) ) CYCLE
        END IF
        
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
        SaveNodes(4) = SaveNodes(4) + 1
        

        No = 0
        Values = 0.0d0

        CALL WriteFieldsAtElement( CurrentElement, MaxBoundary, k, 0, Basis, linepos = -1.0_dp )         
      END DO

      WRITE( Message, * ) 'Number of nodes in isocurves: ', SaveNodes(4)
      CALL Info(Caller,Message)
         
    END DO

    IF(ALLOCATED(LineTag)) THEN
      DEALLOCATE( LineTag, STAT=istat)
      IF(istat /= 0) CALL Fatal(Caller,'Could not deallocate LineTag!')
    END IF
    
  END SUBROUTINE SaveIsoCurves

  
  SUBROUTINE SaveVariableNames()

    INTEGER :: NamesUnit
    
    ! Finally save the names of the variables to help to identify the 
    ! columns in the result matrix.
    !-----------------------------------------------------------------
    IF( Solver % TimesVisited == 0 .AND. NoResults > 0 .AND. &
        (.NOT. Parallel .OR. ParEnv % MyPe == 0 ) ) THEN
      ALLOCATE( ValueNames(NoResults+5), STAT=istat )
      IF( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for ValueNames') 
      
      No = 0
      DO ivar = -2,NoVar
        Var => VariableGetN( ivar, comps )

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
          DGVar = ( Var % TYPE == variable_on_nodes_on_elements )
          IpVar = ( Var % TYPE == variable_on_gauss_points ) 
          ElemVar = ( Var % TYPE == variable_on_elements ) 
          EdgeBasis = ( Var % Type == variable_on_edges ) 
          
          IF( EdgeBasis ) THEN
            ValueNames(No+1) = TRIM(Var % Name)//' {e} 1'
            ValueNames(No+2) = TRIM(Var % Name)//' {e} 2'
            ValueNames(No+3) = TRIM(Var % Name)//' {e} 3'         
            No = No + 3
            IF( AVBasis ) THEN
              No = No + 1
              ValueNames(No) = TRIM(Var % Name)//' nodal'
            END IF
          ELSE IF( comps > 1 ) THEN
            No = No + 1
            ValueNames(No) = TRIM(Var % Name)             
            IF( comps >= 2 ) THEN
              Var => VariableGetN( ivar, component = 2 ) 
              No = No + 1
              ValueNames(No) = TRIM(Var % Name)                           
            END IF
            IF( comps >= 3 ) THEN
              Var => VariableGetN( ivar, component = 3 ) 
              No = No + 1
              ValueNames(No) = TRIM(Var % Name)                           
            END IF
          ELSE IF( Var % Dofs == 1 ) THEN
            No = No + 1
            ValueNames(No) = TRIM(Var % Name)
          ELSE
            DO i=1,Var % Dofs
              No = No + 1
              ValueNames(No) = ComponentName(Var % Name,i)
            END DO
          END IF
        END IF
        Var => Var % Next      
      END DO

      IF ( CalculateFlux ) THEN
        ValueNames(No+1) = 'Flux 1'
        ValueNames(No+2) = 'Flux 2'
        ValueNames(No+3) = 'Flux normal'      
      END IF

      SideNamesFile = TRIM(SideFile) // '.' // TRIM("names")
      OPEN (newunit=NamesUnit, FILE=SideNamesFile)

      Message = ListGetString(Model % Simulation,'Comment',GotIt)
      IF( GotIt ) THEN
        WRITE(NamesUnit,'(A)') TRIM(Message)
      END IF
      Message = ListGetString(Params,'Comment',GotIt)
      IF( GotIt ) THEN
        WRITE(NamesUnit,'(A)') TRIM(Message)
      END IF
      WRITE(NamesUnit,'(A,A)') 'Metadata for SaveLine file: ',TRIM(SideFile)

      DateStr = GetVersion()
      WRITE( NamesUnit,'(A)') 'Elmer version: '//TRIM(DateStr)     
      DateStr = GetRevision( GotIt )
      IF( GotIt ) THEN
        WRITE( NamesUnit,'(A)') 'Elmer revision: '//TRIM(DateStr)
      END IF        
      DateStr = GetCompilationDate( GotIt )
      IF( GotIt ) THEN
        WRITE( NamesUnit,'(A)') 'Elmer compilation date: '//TRIM(DateStr)
      END IF

      DateStr = GetSifName( GotIt )
      IF( GotIt ) THEN
        WRITE( NamesUnit,'(A)') 'Solver input file: '//TRIM(DateStr)
      END IF
      
      DateStr = FormatDate()
      WRITE( NamesUnit,'(A,A)') 'File started at: ',TRIM(DateStr)

      WRITE(NamesUnit,'(A)') 'Number of data nodes for each step'
      WRITE(NamesUnit,'(A)') '  bc nodes: '//I2S(SaveNodes(1))
      WRITE(NamesUnit,'(A)') '  polyline nodes: '//I2S(SaveNodes(2))
      WRITE(NamesUnit,'(A)') '  circle nodes: '//I2S(SaveNodes(3))
      WRITE(NamesUnit,'(A)') '  isocurve nodes: '//I2S(SaveNodes(4))

      WRITE(NamesUnit,'(A)') 'Data on different columns'
      j = 0
      IF( .NOT. SkipBoundaryInfo ) THEN
        IF( LineInd /= 0 ) THEN
          j = j+1
          WRITE(NamesUnit,'(I3,": ",A)') j,'Line Marker'
        END IF
        IF(TransientSimulation) THEN
          j = j+1
          WRITE(NamesUnit,'(I3,": ",A)') j,'Timestep'
        END IF
        WRITE(NamesUnit,'(I3,": ",A)') 1+j,'Call count'
        WRITE(NamesUnit,'(I3,": ",A)') 2+j,'Boundary condition'
        WRITE(NamesUnit,'(I3,": ",A)') 3+j,'Node index'
        j = j + 3
      END IF
      DO i=1,NoResults
        WRITE(NamesUnit,'(I3,": ",A)') i+j,TRIM(ValueNames(i))
      END DO
      CLOSE(NamesUnit)
      DEALLOCATE( ValueNames, STAT=istat)
      IF(istat /= 0) CALL Fatal(Caller,'Could not deallocate ValueNames!')
    END IF

  END SUBROUTINE SaveVariableNames

!------------------------------------------------------------------------------
END SUBROUTINE SaveLine
!------------------------------------------------------------------------------
!> \}
