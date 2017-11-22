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
! ******************************************************************************************/
! *
! *  Authors: Juha Ruokolainen & Peter Raback
! *  Email:   Juha.Ruokolainen@csc.fi & Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: Serial 2007 & Parallel 2010
! *
! *****************************************************************************/


!--------------------------------------------------------------------------------------------
!> This solver is based on ProjectToPlane.src, located in fem/src/modules/
!> After gathering groups of vertical intersections, instead of integrating,
!> it searches for the % fractured ice based on given criteria.
!>
!> There was originally a serial and parallel version of this code, but I've removed the
!> serial because it's not needed.
!> \ingroup Solvers
!--------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------
SUBROUTINE ProjectCalving( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE GeneralUtils
  USE ElementDescription

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  TYPE(Element_t), POINTER :: Element, FaceElement
  TYPE(Element_t), TARGET :: TriangleElement
  TYPE(Nodes_t) :: Nodes, LineNodes, FaceNodes, ElementNodes
  TYPE(Variable_t), POINTER :: Variable2D, Var, StressVar, PwVar, FlowVar, &
       SurfCrevVar, BasalCrevVar
  TYPE(ValueList_t), POINTER :: Params, Material
  TYPE(Mesh_t), POINTER :: Mesh2D, Mesh3D, Mesh
  REAL(KIND=dp), POINTER :: StressValues(:), PwValues(:)
  REAL(KIND=dp), POINTER :: Basis(:)
  REAL(KIND=dp), POINTER :: PlaneX(:), PlaneY(:), PlaneZ(:)
  REAL(KIND=dp), POINTER :: VolumeX(:), VolumeY(:), VolumeZ(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at, totcpu
#else
  REAL(KIND=dp) :: CPUTime, at, totcpu
#endif
  REAL(KIND=dp) :: x0,y0,z0,xmin,xmax,ymin,ymax,zmin,zmax
  REAL(KIND=dp) :: scale, eps, x1,x2,y1,y2,z1,z2,r1,r2
  REAL(KIND=dp) :: up,vp,wp,cp,cf,MaxRelativeRadius,detJ
  REAL(KIND=dp) :: LocalPoint(3), LocalCoords(3), RhoWS, RhoWF, g, &
       SeaLevel, LocalM(3,3),EigValues(3),EI(3),dumy,dumy2,work(27)
  REAL(KIND=dp), POINTER :: MinHeight3D(:), MaxHeight3D(:), MinWidth3D(:), MaxWidth3D(:)
  REAL(KIND=dp), POINTER :: IntValues(:,:), IntBasis(:,:), IntExtent(:),EigenStress(:)
  INTEGER :: MaxInt, Int
  INTEGER, POINTER :: Perm2D(:), Perm3D(:), PlanePerm(:), VolumePerm(:), &
       StressPerm(:), PwPerm(:)
  INTEGER :: i,j,k,k2,l,n,t,lnode,node, StressDOFs
  INTEGER :: PlaneNodes, VolumeElements, DOFs_2D, DOFs_3D, corners, face, Intersections
  INTEGER :: Loops(8), inds(3)
  INTEGER, POINTER :: Order2D(:), Order3D(:)
  INTEGER, POINTER :: IntOrder(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: ConvertFromName, ConvertFromVar, Name, Str, &
       SolverName, StressVarName, PwVarName, FlowVarName
  LOGICAL :: AllocationsDone = .FALSE., stat, GotIt, Found, &
       PwFromVar, Cauchy, SurfModel, BasalModel

  INTEGER :: pe, myNodes, myStart, peNodes, peNodesn, peStart, peEnd, &
          totcount, curr, ierr, status(MPI_STATUS_SIZE), nn

  TYPE ValueTable_t
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: Perm(:)
  END TYPE ValueTable_t
  TYPE(ValueTable_t) :: ValueTable3D(2), ValueTable2D(2)

  TYPE PointStore_t
    INTEGER :: Int
    REAL(KIND=dp),POINTER :: IntExtent(:), IntValues(:,:)
  END TYPE PointStore_t
  TYPE(PointStore_t), ALLOCATABLE :: PointStore(:)

  INTEGER, ALLOCATABLE :: cm_int(:), cm_int0(:), icount(:), icount0(:)
  REAL(KIND=dp), ALLOCATABLE :: cm_values(:), cm_values0(:), cm_extent(:)

  SAVE :: Variable2D, &
      VolumeElements, PlaneNodes, VolumeX, VolumeY, VolumeZ, PlaneX, PlaneY, &
      PlaneZ, MaxInt, IntBasis, IntValues, IntExtent, IntOrder, &
      MinHeight3D, MaxHeight3D, MinWidth3D, MaxWidth3D, &
      Basis, ElementNodes, LineNodes, FaceNodes, Eps, &
      Mesh2D, Mesh3D, Perm2D, Perm3D, Dofs_2D, DOFs_3D, ValueTable2D, &
      ValueTable3D


!------------------------------------------------------------------------------

  totcpu = CPUTime()
  SolverName = "ProjectCalving"
  Params => Solver % Values

  CALL Info( SolverName, ' ' )
  CALL Info( SolverName, '-----------------------------------' )
  CALL Info( SolverName, ' Collapsing CIndex to 2D plane ' )
  CALL Info( SolverName, '-----------------------------------' )
  CALL Info( SolverName, ' ' )

  !------------------------------------------------------------------------------
  ! Find the 2D and 3D meshes, assume just one of each
  !------------------------------------------------------------------------------
  NULLIFY(Mesh2D)
  NULLIFY(Mesh3D)
  Mesh => Model % Meshes
  DO WHILE( ASSOCIATED(Mesh) )
     IF ( .NOT. Mesh % OutputActive ) THEN
        Mesh => Mesh % next
        CYCLE
     END IF
     IF( Mesh % MeshDim == 3 ) THEN
        CALL Info(SolverName,'3D mesh: '//TRIM(Mesh % Name) )
        Mesh3D => Mesh
     ELSE IF( Mesh % MeshDim == 2 ) THEN
        CALL Info(SolverName,'2D mesh: '//TRIM(Mesh % Name) )
        Mesh2D => Mesh
     END IF
     Mesh => Mesh % next
  END DO

  IF(.NOT. ASSOCIATED(Mesh2D)) THEN
     CALL Fatal(SolverName,'Did not find 2D mesh')
  END IF
  IF(.NOT. ASSOCIATED(Mesh3D)) THEN
     CALL Fatal(SolverName,'Did not find 3D mesh')
  END IF

  IF(.NOT. ASSOCIATED(Mesh2D, Solver % Mesh)) THEN
     CALL Fatal(SolverName,'2D mesh not same as Solver mesh')
  END IF


  !------------------------------------------------------------------------------
  ! Get stress and water pressure variables, gravity, water density, sea level
  !------------------------------------------------------------------------------
  g = GetConstReal( Model % Constants, "g", Found )
  IF(.NOT.Found) CALL Fatal("Calving", "Unable to find Gravity.")
  IF(g < 0.0_dp) g = -1.0_dp * g
  RhoWS = GetConstReal( Model % Constants, "RhoWS", Found )
  IF(.NOT.Found) CALL Fatal("Calving", "Unable to find RhoWS.")
  RhoWF = GetConstReal( Model % Constants, "RhoWF", Found )
  IF(.NOT.Found) CALL Fatal("Calving", "Unable to find RhoWF.")
  SeaLevel = GetConstReal( Model % Constants, "Sea Level", Found )
  IF(.NOT.Found) THEN
     WRITE(Message,'(A)') 'Variable Sea level not found. &
          &Setting to 0.0'
     CALL INFO(SolverName, Message, level=2)
     SeaLevel = 0.0_dp
  END IF

  StressVarName = ListGetString(Params,"Calving Stress Variable Name", Found)
  IF(.NOT. Found) THEN
     StressVarName = "Stress"
     WRITE( Message, * ) "No Calving Stress Variable Name given, using ", StressVarName
     CALL Info(SolverName, Message)
  END IF

  StressVar => VariableGet( Mesh3D % Variables, StressVarName, .TRUE. )
  IF(.NOT. ASSOCIATED(StressVar)) THEN
     CALL Fatal(SolverName, "Couldn't find stress variable")
  ELSE
     StressValues => StressVar % Values
     StressPerm => StressVar % Perm
     StressDOFs = StressVar % DOFs
  END IF

  Material => GetMaterial(Mesh3D % Elements(1)) !body elements first...
  Cauchy = ListGetLogical( Material , 'Cauchy', Found )
  IF(.NOT. Found) THEN
     Cauchy = .FALSE.
     CALL Warn(SolverName, "Couldn't find 'cauchy' logical in material, &
          &assuming deviatoric stress.")
  END IF

  !Deviatoric stress, need to add pressure to get cauchy
  IF(.NOT. Cauchy) THEN
     FlowVarName = ListGetString(Params, "Flow Solution Variable Name",Found)
     IF(.NOT. Found) FlowVarName = "Flow Solution"

     FlowVar => VariableGet(Mesh3D % Variables, FlowVarName, .TRUE.)
     IF(.NOT. ASSOCIATED(FlowVar)) CALL Fatal(SolverName,"Working with deviatoric stress, &
          &but no flow solution found")
  END IF

  PwVarName = ListGetString(Params, "Water Pressure Variable Name", Found)
  IF(.NOT. Found) THEN
     WRITE( Message, * ) "No Water Pressure Variable Name given, &
          &using standard sea level calculation"
     CALL Info(SolverName, Message)
     PwFromVar = .FALSE.
  ELSE
     PwFromVar = .TRUE.
  END IF

  !Determine the crevasse model
  SurfModel = ListGetLogical(Params, "Surface Crevasse Model", Found)
  IF(.NOT. Found) SurfModel = .TRUE.
  BasalModel = ListGetLogical(Params, "Basal Crevasse Model", Found)
  IF(.NOT. Found) BasalModel = .TRUE.
  IF(SurfModel .AND. BasalModel) THEN
    CALL Info(SolverName, "Calving occurs either when surface crevasses meet sea-level, &
         &or when surface and basal crevasses meet.")
  ELSE IF(SurfModel) THEN
    CALL Info(SolverName, "Calving occurs when surface crevasses meet sea-level.")
  ELSE IF(BasalModel) THEN
    CALL Info(SolverName, "Calving occurs when surface and basal crevasses meet.")
  ELSE
    CALL Fatal(SolverName, "Calling calving solver with both 'Surface Crevasse Model' and &
         &'Basal Crevasse Model' set to false...")
  END IF

  !Compute max eigenstress for all nodes
  !assumption here that stress is computed at every node...

  ALLOCATE(EigenStress(SIZE(StressPerm)))
  DO i=1,Mesh3D % NumberOfNodes

     localM=0.0_dp

     !xx,yy,zz,xy,yz,zx
     localM(1,1)=StressValues( StressDOFs * (StressPerm(i) - 1 ) + 1) !xx
     localM(2,2)=StressValues( StressDOFs * (StressPerm(i) - 1 ) + 2) !yy
     localM(3,3)=StressValues( StressDOFs * (StressPerm(i) - 1 ) + 3) !zz
     localM(1,2)=StressValues( StressDOFs * (StressPerm(i) - 1 ) + 4) !xy
     localM(2,3)=StressValues( StressDOFs * (StressPerm(i) - 1 ) + 5) !zy
     localM(1,3)=StressValues( StressDOFs * (StressPerm(i) - 1 ) + 6) !xz

     !and fill
     localM(2,1)=localM(1,2)
     localM(3,1)=localM(1,3)
     localM(3,2)=localM(2,3)

     IF(.NOT. Cauchy) THEN
        DO j=1,3                                       !dim+1, only runs in 3D
           localM(j,j)=localM(j,j) - FlowVar % Values(4*FlowVar % Perm(i))
        END DO
     END IF

     CALL DGEEV('N','N',3,localM,3,EigValues,EI,Dumy,1,Dumy2,1,Work,27,ierr )
     IF (ierr/=0) THEN
        WRITE(Message,'(A,i0)') 'Failed to compute EigenValues, error code: ',ierr
        CALL FATAL(SolverName, Message)
     END IF

     EigenStress(StressPerm(i)) = MAXVAL(EigValues) !Maximum principal stress
  END DO

  ValueTable3D(1) % Values => EigenStress
  ValueTable3D(1) % Perm => StressPerm

  IF(PwFromVar) THEN
     PwVar => VariableGet( Mesh3D % Variables, PwVarName, .TRUE. )
     IF(.NOT. ASSOCIATED(PwVar)) THEN
        CALL Fatal(SolverName, "Couldn't find Water Pressure variable")
     ELSE
        PwValues => PwVar % Values
        PwPerm => PwVar % Perm

        ValueTable3D(2) % Values => PwValues
        ValueTable3D(2) % Perm => PwPerm
     END IF
  END IF

  SurfCrevVar => VariableGet( Mesh2D % Variables, "surf_cindex", .TRUE., UnfoundFatal=.TRUE.)
  BasalCrevVar => VariableGet( Mesh2D % Variables, "basal_cindex", .TRUE., UnfoundFatal=.TRUE.)

  Variable2D => Solver % Variable

  ValueTable2D(1) % Values => SurfCrevVar % Values
  ValueTable2D(2) % Values => BasalCrevVar % Values

  DOFs_3D = 1
  IF(PwFromVar) DOFs_3D = 2
  DOFs_2D = SIZE(ValueTable2D)

  !------------------------------------------------------------------------------
  ! The permutation are assumed to be constant
  !------------------------------------------------------------------------------
  Perm3D => StressVar % Perm
  Perm2D => Variable2D % Perm

  PlaneNodes = COUNT(Perm2D>0)
  VolumeElements = Mesh3D % NumberOfBulkElements

  WRITE( Message, * ) 'Number of mesh nodes in 2D and 3D: ', &
       Mesh2d % NumberOfNodes, Mesh3D % NumberOfNodes
  CALL Info( SolverName, Message, LEVEL=16 )

  !------------------------------------------------------------------------------
  ! Possible permutation of coorinate directions
  !------------------------------------------------------------------------------
  VolumeX => Mesh3D % Nodes % x
  VolumeY => Mesh3D % Nodes % y
  VolumeZ => Mesh3D % Nodes % z

  VolumePerm => ListGetIntegerArray( Solver % Values,'Volume Permutation',GotIt)
  IF ( gotIt ) THEN
     IF(VolumePerm(1) == 2) VolumeX => Mesh3D % Nodes % y
     IF(VolumePerm(1) == 3) VolumeX => Mesh3D % Nodes % z
     IF(VolumePerm(2) == 1) VolumeY => Mesh3D % Nodes % x
     IF(VolumePerm(2) == 3) VolumeY => Mesh3D % Nodes % z
     IF(VolumePerm(3) == 1) VolumeZ => Mesh3D % Nodes % x
     IF(VolumePerm(3) == 2) VolumeZ => Mesh3D % Nodes % y
  END IF

  PlaneX => Mesh2d % Nodes % x
  PlaneY => Mesh2d % Nodes % y
  PlaneZ => Mesh2d % Nodes % z

  PlanePerm => ListGetIntegerArray( Solver % Values,'Plane Permutation',GotIt)
  IF ( gotIt ) THEN
     IF(PlanePerm(1) == 2) PlaneX => Mesh2d % Nodes % y
     IF(PlanePerm(1) == 3) PlaneX => Mesh2d % Nodes % z
     IF(PlanePerm(2) == 1) PlaneY => Mesh2d % Nodes % x
     IF(PlanePerm(2) == 3) PlaneY => Mesh2d % Nodes % z
     IF(PlanePerm(3) == 1) PlaneZ => Mesh2d % Nodes % x
     IF(PlanePerm(3) == 2) PlaneZ => Mesh2d % Nodes % y
  END IF

  !------------------------------------------------------------------------------
  !   Allocate stuff
  !------------------------------------------------------------------------------

  ALLOCATE( MinHeight3D(VolumeElements), MaxHeight3D(VolumeElements))
  ALLOCATE( MinWidth3D(VolumeElements), MaxWidth3D(VolumeElements))

  n = Mesh3D % MaxElementNodes
  ALLOCATE(ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), Basis(n) )
  ALLOCATE(LineNodes % x(2), LineNodes % y(2), LineNodes % z(2)  )
  ALLOCATE(FaceNodes % x(3), FaceNodes % y(3), FaceNodes % z(3))

  !------------------------------------------------------------------------------
  !   Find the scale of the 2D mesh
  !------------------------------------------------------------------------------

  xmin = MINVAL( PlaneX )
  xmax = MAXVAL( PlaneX )
  ymin = MINVAL( VolumeY ) !Or else it starts from 0.0 in height direction, not good.
  ymax = MAXVAL( VolumeY )
  zmin = MINVAL( PlaneZ )
  zmax = MAXVAL( PlaneZ )

  x0 = xmax - xmin
  y0 = ymax - ymin
  z0 = zmax - zmin

  scale = SQRT(x0*x0 + y0*y0 + z0*z0)
  Eps = 1.0d-6 * scale

  LineNodes % y(1) = ymin - Eps
  LineNodes % y(2) = ymax + Eps !eps not really necessary here

  !------------------------------------------------------------------------------
  !   To improve search speed, tabulate values for min and max values of element nodes
  !------------------------------------------------------------------------------

  DO t = 1, VolumeElements
     Element => Mesh3D % Elements( t )
     Model % CurrentElement => Element
     n = GetElementNOFNodes(Element)

     ElementNodes % x(1:n) = VolumeX(Element % NodeIndexes(1:n))
     ElementNodes % y(1:n) = VolumeY(Element % NodeIndexes(1:n))
     ElementNodes % z(1:n) = VolumeZ(Element % NodeIndexes(1:n))

     MinHeight3D(t) = MINVAL(ElementNodes % z(1:n))
     MaxHeight3D(t) = MAXVAL(ElementNodes % z(1:n))

     MinWidth3D(t) = MINVAL(ElementNodes % x(1:n))
     MaxWidth3D(t) = MAXVAL(ElementNodes % x(1:n))
  END DO


  !
  ! allocate space for 3d face intersections / 2d point:
  ! ----------------------------------------------------
  MaxInt =  100
  ALLOCATE( PointStore(PlaneNodes) )
  DO i=1,PlaneNodes
    ALLOCATE( PointStore(i) % IntValues(DOFs_3D,MaxInt), &
              PointStore(i) % IntExtent(MaxInt) )
    PointStore(i) % Int = 0
    PointStore(i) % IntExtent = 0
    PointStore(i) % IntValues = 0
  END DO


  ! By construction split everything into triangles
  corners = 3

  TriangleElement % TYPE => GetElementType( 303, .FALSE. )

  FaceElement => TriangleElement
  Loops = 0


  !   Loop over 2D nodes:
  !   -------------------
  DO node = 1, Mesh2d % NumberOfNodes

    Loops(1) = Loops(1) + 1
    lnode = Perm2D(node)
    IF(lnode == 0) CYCLE

    x0 = PlaneX(node)
    y0 = PlaneY(node)
    z0 = PlaneZ(node)

    LineNodes % x(1) = x0
    LineNodes % x(2) = x0
    LineNodes % z(1) = z0
    LineNodes % z(2) = z0

    ! Loop over 3D elements:
    ! ---------------------
    Int = 0
    DO t = 1, VolumeElements

      Loops(2) = Loops(2) + 1

      IF(z0 > MaxHeight3D(t) + Eps) CYCLE
      IF(z0 < MinHeight3D(t) - Eps) CYCLE

      IF(x0 > MaxWidth3D(t) + Eps) CYCLE
      IF(x0 < MinWidth3D(t) - Eps) CYCLE

      Loops(3) = Loops(3) + 1

!------------------------------------------------------------------------------
      Element => Mesh3D % Elements( t )
      Model % CurrentElement => Element
      IF(.NOT. ASSOCIATED(Element)) CYCLE

      n = GetElementNOFNodes(Element)
      ElementNodes % x(1:n) = VolumeX( Element % NodeIndexes(1:n) )
      ElementNodes % y(1:n) = VolumeY( Element % NodeIndexes(1:n) )
      ElementNodes % z(1:n) = VolumeZ( Element % NodeIndexes(1:n) )
!------------------------------------------------------------------------------

      Intersections = 0
      DO face=1, 12
        CALL GetLinearTriangleFaces( Element, face, inds, GotIt )
        IF(.NOT. GotIt) EXIT

        Loops(4) = Loops(4) + 1

        FaceNodes % x(1:corners) = ElementNodes % x(inds(1:corners))
        FaceNodes % y(1:corners) = ElementNodes % y(inds(1:corners))
        FaceNodes % z(1:corners) = ElementNodes % z(inds(1:corners))

        xmin = MINVAL( FaceNodes % x(1:corners) )
        xmax = MAXVAL( FaceNodes % x(1:corners) )
        IF( x0 < xmin - Eps) CYCLE
        IF( x0 > xmax + Eps) CYCLE

        zmin = MINVAL( FaceNodes % z(1:corners) )
        zmax = MAXVAL( FaceNodes % z(1:corners) )
        IF( z0 < zmin - Eps) CYCLE
        IF( z0 > zmax + Eps) CYCLE

        IF ( .NOT. LineFaceIntersect( Element,FaceNodes,corners, &
                   LineNodes,up,vp,cp) ) CYCLE

        Loops(5) = Loops(5) + 1

        Intersections = Intersections + 1
        CALL NodalBasisFunctions(Corners, Basis, FaceElement, up, vp, 0.0d0)

        !
        ! store extent of the line + value(s) at point of
        ! intersection:
        ! ------------------------------------------------
        PointStore(lnode) % Int = PointStore(lnode) % Int + 1
        curr=PointStore(lnode) % Int
        IF (curr>SIZE(PointStore(lnode) % IntExtent)) &
            CALL AllocateMoreSpace(PointStore(lnode), MaxInt)

        PointStore(lnode) % IntExtent(curr) = LineNodes % y(1) + &
             (cp * (LineNodes % y(2) - LineNodes % y(1)))

        DO l=1,DOFs_3D
          PointStore(lnode) % IntValues(l,curr) = &
              SUM(Basis(1:corners)* ValueTable3D(l) % &
              Values(ValueTable3D(l) % Perm( Element % NodeIndexes(inds(1:corners)) )) )
        END DO

      END DO

      IF(Intersections /= 0 .AND. Intersections /= 2) Loops(6)=Loops(6)+1
    END DO
  END DO

  !
  ! divide 2d points to pe:s; send intersections
  ! to owners of the points; receive our parts
  ! --------------------------------------------
  myNodes  = PlaneNodes / ParEnv % PEs
  myStart=ParEnv % myPE*myNodes+1
  IF ( ParEnv % mype==Parenv % PEs-1) &
      myNodes = myNodes + (PlaneNodes-Parenv % Pes*myNodes)

  IF ( Parenv % PEs>1 ) THEN
    CALL CheckBuffer( (2 * (2 * (DOFs_3D+1)*SUM(PointStore(:) % Int)) + PlaneNodes ))
    CALL SendPoints()
    CALL ReceivePoints()
    CALL CheckBuffer(1) !Just frees the buffer space for future subroutines.
  END IF

  ! compute integrals from intersections:
  ! -------------------------------------
  CALL ComputeCrevassePenetration()

  ! send all results to pe 0:
  ! -------------------------
  IF (ParEnv % PEs>1) CALL GatherResults()

  !Fill variable values depending on crevasse model:
  IF(SurfModel .AND. BasalModel) THEN
    !MIN works on arrays, hooray.
    Variable2D % Values = MIN(ValueTable2D(1) % Values, ValueTable2D(2) % Values)
  ELSE IF(SurfModel) THEN
    Variable2D % Values = ValueTable2D(1) % Values
  ELSE !BasalModel, checked at least one model is active in initialization
    Variable2D % Values = ValueTable2D(2) % Values
  END IF

  DO i=1,PlaneNodes
    DEALLOCATE(PointStore(i) % IntValues, PointStore(i) % IntExtent)
  END DO
  DEALLOCATE( PointStore )
  DEALLOCATE( EigenStress )
  DEALLOCATE(MinHeight3D, MaxHeight3D)
  DEALLOCATE(MinWidth3D, MaxWidth3D)
  DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z, Basis )
  DEALLOCATE(LineNodes % x, LineNodes % y, LineNodes % z)
  DEALLOCATE(FaceNodes % x, FaceNodes % y, FaceNodes % z)

  WRITE( Message,'(A,4I8)' ) 'Basic search loops: ',Loops(2:5)
  CALL Info( SolverName, Message, LEVEL=4 )

  IF(Loops(5) > 0) THEN
    WRITE( Message,'(A,3I8)' ) 'Special cases loops: ',Loops(6:8)
    CALL Info( SolverName, Message, LEVEL=4 )
  END IF

  WRITE( Message,'(A,F8.2)' ) 'Total CPU time used: ', CPUTime() - totcpu
  CALL INFO( SolverName, Message, LEVEL=8 )
  CALL Info( SolverName, ' ' )
  CALL Info( SolverName, 'All done' )
  CALL Info( SolverName, ' ' )

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE SendPoints()
!------------------------------------------------------------------------------
    !
    ! send intersections to owners of the points:
    ! -------------------------------------------
    peNodes  = PlaneNodes / ParEnv % PEs
    peNodesn = peNodes+(PlaneNodes-Parenv % Pes*peNodes)

    DO pe=0,ParEnv % PEs-1
      IF ( pe==Parenv % mype ) CYCLE

      peStart = pe*peNodes+1
      IF ( pe==ParEnv % PEs-1 ) peNodes=peNodesn
      peEnd = peStart+peNodes-1

      totcount=SUM(PointStore(peStart:peEnd) % Int)
      ALLOCATE( cm_int(peNodes) ); cm_int=0

      IF ( totcount>0 ) THEN
        ALLOCATE( cm_extent(totcount), &
              cm_values(DOFs_3D*totcount) )

        j=0; totcount=0
        DO lnode=peStart,peEnd
          j = j + 1
          cm_int(j) = PointStore(lnode) % Int
          DO k=1,cm_int(j)
            totcount=totcount+1
            int=DOFs_3D*(totcount-1)
            cm_values(int+1:int+DOFs_3D) = &
               PointStore(lnode) % IntValues(1:DOFs_3D,k)
            cm_extent(totcount) = PointStore(lnode) % IntExtent(k)
           END DO
        END DO
      END IF
      CALL MPI_BSEND( cm_int, peNodes, MPI_INTEGER, pe, &
                100, ELMER_COMM_WORLD, ierr )
      IF (totcount>0 ) THEN
        CALL MPI_BSEND( cm_values, DOFs_3D*totcount, &
           MPI_DOUBLE_PRECISION, pe, 101, ELMER_COMM_WORLD, ierr )

        CALL MPI_BSEND( cm_extent, totcount, &
           MPI_DOUBLE_PRECISION, pe, 102, ELMER_COMM_WORLD, ierr )
        DEALLOCATE( cm_values, cm_extent)
      END IF
      DEALLOCATE(cm_int )
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SendPoints
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReceivePoints()
!------------------------------------------------------------------------------
    !
    ! receive intersections for the 2d points we own:
    ! -----------------------------------------------
    ALLOCATE(cm_int(myNodes))
    DO pe=0,ParEnv % PEs-1
      IF ( pe==Parenv % mype ) CYCLE

      CALL MPI_RECV( cm_int, myNodes, MPI_INTEGER, pe, &
             100, ELMER_COMM_WORLD, status, ierr )

      totcount=SUM(cm_int(1:myNodes))

      IF ( totcount>0 ) THEN
        ALLOCATE( cm_values(DOFs_3D*totcount), cm_extent(totcount) )

        CALL MPI_RECV( cm_values, DOFs_3D*totcount, &
          MPI_DOUBLE_PRECISION, pe, 101, ELMER_COMM_WORLD, status, ierr )

        CALL MPI_RECV( cm_extent, totcount, &
          MPI_DOUBLE_PRECISION, pe, 102, ELMER_COMM_WORLD, status, ierr )

        j = 0
        totcount = 0
        DO lnode=myStart,myStart+myNodes-1
          j=j+1
          curr = PointStore(lnode) % Int
          IF ( curr+cm_int(j)>SIZE(PointStore(lnode) % IntExtent)) &
            CALL AllocateMoreSpace(PointStore(lnode),cm_int(j))

          DO k=1,cm_int(j)
            totcount=totcount+1
            curr=curr+1
            int = DOFs_3D*(totcount-1)
            PointStore(lnode) % IntExtent(curr) = cm_extent(totcount)
            PointStore(lnode) % IntValues(:,curr) = cm_values(int+1:int+DOFs_3D)
          END DO
          PointStore(lnode) % Int=curr
        END DO
        DEALLOCATE(cm_values, cm_extent)
      END IF
    END DO
    DEALLOCATE(cm_int)
!------------------------------------------------------------------------------
  END SUBROUTINE ReceivePoints
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ComputeCrevassePenetration()
!------------------------------------------------------------------------------
    !Loop over all the 2D nodes, find the zero contours of cindex (surf and basal)
    !IntValues(1,:) are stress, IntValues(2,:) are Pw, if PwFromVar = .TRUE.

    REAL(KIND=dp) :: dz_dsigma, dsigma_dz, dpw_dz, dz_dpw, dcindex_dz, dz_dcindex,&
         Pw1, Pw2, PW_base, CIndex1, CIndex2, Stress1, Stress2, BContour, SContour, &
         ZeroContour, myEps

    LOGICAL :: Found,Waterline,Crevasseline

    myEps = 1.0E-10

    !Whether Pw given by variable or calculated from sea level, gradient is fixed
    dpw_dz = -1.0 * g * RhoWF
    dz_dpw = 1.0_dp/dpw_dz

    DO node=1,Mesh2d % NumberOfNodes
      lnode = Perm2D(node)
      IF(lnode<myStart .OR. lnode>myStart+myNodes-1) CYCLE

      Int = PointStore(lnode) % Int
      IF (Int<2) THEN
         WRITE(Message, *) "Only one hit for 2D node: ", node

         ValueTable2D(1) % Values(lnode) = 0.0_dp
         ValueTable2D(2) % Values(lnode) = 0.0_dp
         CYCLE
      END IF

      IntValues => PointStore(lnode) % IntValues
      IntExtent => PointStore(lnode) % IntExtent

      !NOTE, the extents are good, I've checked them, so we know cp etc is working properly.
      Loops(7) = Loops(7) + 1

      !Sort intersections in descending order
      ALLOCATE(IntOrder(Int)); IntOrder = [(i,i=1,Int)]
      CALL SortR(Int, IntOrder, IntExtent)

      !---------------------------------------------------
      ! Search downwards for CSurfIndex zero contour
      ! CSurfIndex = EigenStress_3, stress is in IntValues(1,:)
      !---------------------------------------------------

      !set csurf contour initially to top node height
      SContour = IntExtent(1)

      !Cycle nodes in descending order (moving down from top surface)
      IF(IntValues(1,IntOrder(1)) > 0.0_dp) THEN
         DO j=1,Int-1
            k = IntOrder(j)
            k2 = IntOrder(j+1)

            Stress1 = IntValues(1,k)
            Stress2 = IntValues(1,k2)

            IF(Stress1 < 0.0_dp) CALL Fatal(SolverName, &
                 "Overshot somehow when searching for surface calving contour")
            IF(Stress2 > 0.0_dp) CYCLE

            !hit multiple identical nodes at join
            IF(ABS(IntExtent(j) - IntExtent(j+1)) < myEps) THEN
               IF(Stress2 > 0.0_dp) THEN
                  CYCLE
               ELSE
                  SContour = IntExtent(j)
                  EXIT
               END IF
            END IF

            !Reach this point when we are between +ve and -ve nodes
            !Get gradient
            dz_dsigma = (IntExtent(j) - IntExtent(j+1)) / (Stress1 - Stress2)
            SContour = IntExtent(j) - dz_dsigma*Stress1
            IF((SContour > IntExtent(j)) .OR. (SContour < IntExtent(j+1))) &
                 CALL Fatal(SolverName, "Wrong gradient calculation, this is a programming error.")
            EXIT
         END DO
      END IF

      !-------------------------------------------------------------
      ! Search upwards for CBasalIndex zero contour
      ! CBasalIndex = EigenStress_3 + Pw
      ! However, we need to treat the line element vertical gradient
      ! of Stress and Pw *SEPARATELY* or else we will mistake surface
      ! crevasses for basal. See ProjectToPlane_mod.ora for details
      !-------------------------------------------------------------
      BContour = IntExtent(Int)
      Found = .FALSE.

      !Calculate Pw at the base of the column, based on sea level and salt water density
      PW_base = -1.0 * (IntExtent(Int) - SeaLevel) * g * RhoWS

      !Cycle nodes in ascending order
      DO j=Int, 2, -1
         k = IntOrder(j)
         k2 = IntOrder(j-1)

         Stress1 = IntValues(1,k)
         Stress2 = IntValues(1,k2)
         IF(PwFromVar) THEN
            Pw1 = IntValues(2,k)
            Pw2 = IntValues(2,k2)
         ELSE
            !In-crevasse pressure only equals sea pressure at base of ice
            !then drops at a rate of dz * RhoWF, not dz * RhoWS
            Pw1 = PW_base - (IntExtent(j) - IntExtent(Int)) * g * RhoWF
            Pw2 = PW_base - (IntExtent(j-1) - IntExtent(Int)) * g * RhoWF
         END IF

         CIndex1 = Stress1
         CIndex2 = Stress2

         IF(Pw1 > 0.0_dp) CIndex1 = CIndex1 + Pw1
         IF(Pw2 > 0.0_dp) CIndex2 = CIndex2 + Pw2

         !No basal crevasse in lower node, if first time, no crevasse, else an error
         IF(CIndex1 < 0.0_dp) THEN
            IF(j==Int) THEN
               Found = .TRUE.
               EXIT !bottom node has no basal crevasse, we're done
            END IF
            CALL Fatal(SolverName, &
                 "Overshot somehow when searching for basal calving contour")
         END IF

         !special case: full penetration
         IF(IntExtent(j) > SContour) THEN
            BContour = SContour
            EXIT
         END IF

         !hit multiple identical nodes at join
         IF(ABS(IntExtent(j) - IntExtent(j-1)) < myEps) THEN
            IF(CIndex2 > 0.0_dp) THEN
               CYCLE
            ELSE
               BContour = IntExtent(j)
               EXIT
            END IF
         END IF

         !Passing zero contour or passing the waterline
         Waterline = ((Pw2 <= 0.0_dp) .AND. (Pw1 > 0.0_dp))
         Crevasseline = CIndex2 < 0.0_dp

         IF(Crevasseline .OR. Waterline) THEN

            !stress gradient through depth
            dsigma_dz = (Stress1 - Stress2) / (IntExtent(j) - IntExtent(j-1))
            dz_dsigma = 1.0 / dsigma_dz

            !cindex gradient (may include water pressure)
            IF(Pw1 < 0.0_dp) THEN
               dcindex_dz = dsigma_dz
            ELSE
               dcindex_dz = dsigma_dz + dpw_dz
            END IF

            dz_dcindex = 1.0/dcindex_dz

            IF(dz_dcindex < 0.0) THEN !Usual case
              ZeroContour = IntExtent(j) - (dz_dcindex * CIndex1)
            ELSE !cindex increasing in z direction, rare
              IF(CrevasseLine) THEN
                ZeroContour = IntExtent(j-1) !crevasse contour at upper node height
              ELSE
                CYCLE !no crevasse contour yet
              END IF
            END IF

            !apparent full penetration between these two nodes
            !this happens because water pressure gradient may
            !not exactly correspond to values if PwFromVar
            IF(ZeroContour > IntExtent(j-1)) THEN
              IF(Crevasseline) THEN
                !ZeroContour slightly overestimated due to assumed Pw gradient
                !Set to upper node height
                ZeroContour = IntExtent(j-1)
              ELSE
                !full penetration between these nodes, continue looking
                CYCLE
              END IF
            END IF

            BContour = ZeroContour
            Found = .TRUE.
            IF(BContour < IntExtent(j)) THEN
               PRINT *, 'Cindex1, Cindex2: ', Cindex1, ' ',CIndex2
               PRINT *, 'Stress1,2: ', Stress1, Stress2
               PRINT *, 'Pw1, Pw2: ', Pw1,' ',Pw2
               PRINT *, 'dz_dsigma: ', dz_dsigma
               PRINT *, 'dz_dcindex: ', dz_dcindex
               PRINT *, 'dz_dpw: ', dz_dpw
               PRINT *, 'IntExtent1, 2: ', IntExtent(j), IntExtent(j-1)
               CALL Fatal(SolverName, "Wrong basal calculation, this is a programming error.")
            END IF
            EXIT
         END IF
      END DO

      !Basal crevasse model
      !Full penetration
      IF((SContour < (BContour + EPSILON(BContour))) .OR. &
           (.NOT. Found) .OR. &
           (IntExtent(Int) >= IntExtent(1)+EPSILON(IntExtent(1)))) THEN

         ValueTable2D(2) % Values(lnode) = 0.0_dp
      ELSE
         ValueTable2D(2) % Values(lnode) = (SContour - BContour) / (IntExtent(1) - IntExtent(Int))
      END IF

      !Surface crevasse model
      IF(SContour < SeaLevel) THEN
        ValueTable2D(1) % Values(lnode) = 0.0_dp
      ELSE
        ValueTable2D(1) % Values(lnode) = 1.0_dp - &
             ((IntExtent(1) - SContour) / (IntExtent(1) - SeaLevel))
      END IF

      DEALLOCATE(IntOrder)

   END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeCrevassePenetration
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE GatherResults()
!------------------------------------------------------------------------------
    INTEGER :: j,k,l,pe,peNodes,peNodesn,peStart

    IF ( ParEnv % mype==0 ) THEN
      peNodes  = PlaneNodes / ParEnv % PEs
      peNodesn = peNodes+(PlaneNodes-Parenv % Pes*peNodes)

      DO pe=1,ParEnv % PEs-1
        peStart = pe*peNodes+1
        IF ( pe==Parenv % Pes-1 ) peNodes=peNodesn

        j = peStart
        k = j+peNodes-1
        DO l=1,DOFs_2D
          CALL MPI_RECV( ValueTable2D(l) % Values(j:k), peNodes, &
            MPI_DOUBLE_PRECISION, pe, 104, ELMER_COMM_WORLD, status, ierr )
        END DO
      END DO
    ELSE
      j = myStart
      k = j+myNodes-1
      DO l=1,DOFs_2D
        CALL MPI_BSEND( ValueTable2D(l) % Values(j:k), myNodes, &
            MPI_DOUBLE_PRECISION, 0, 104, ELMER_COMM_WORLD, ierr )
      END DO
      Mesh2d % OutputActive = .FALSE.
    END IF
    CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
!------------------------------------------------------------------------------
  END SUBROUTINE GatherResults
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AllocateMoreSpace(PS,incr)
!------------------------------------------------------------------------------
    TYPE(PointStore_t) :: PS
    INTEGER :: incr

    INTEGER :: olds, news
    REAL(KIND=dp), ALLOCATABLE :: TmpExtent(:),TmpValues(:,:)

    olds = SIZE(PS % IntExtent)

    ALLOCATE(TmpValues(DOFs_3D,olds),TmpExtent(olds))
    TmpExtent(1:olds) = PS % IntExtent(1:olds)
    TmpValues(:,1:olds) = PS % IntValues(:,1:olds)
    DEALLOCATE(PS % IntValues, PS % IntExtent )

    news = olds + incr
    ALLOCATE( PS % IntValues(DOFs_3D,news), PS % IntExtent(news) )
    PS % IntValues = 0
    PS % IntExtent = 0

    PS % IntExtent(1:olds) = TmpExtent
    PS % IntValues(:,1:olds) = TmpValues

    DEALLOCATE(TmpValues, TmpExtent)
!------------------------------------------------------------------------------
  END SUBROUTINE AllocateMoreSpace
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! 3D mesh faces.
!------------------------------------------------------------------------------
  SUBROUTINE GetLinearTriangleFaces( Element, face, inds, GotIt )
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Element
    LOGICAL :: GotIt
    INTEGER :: face, inds(:)
!------------------------------------------------------------------------------
    LOGICAL :: Visited
    INTEGER :: faces, elemfamily
    INTEGER, POINTER :: FaceMap(:,:)
    INTEGER, TARGET  :: TetraFaceMap(4,3), BrickFaceMap(12,3), &
                WedgeFaceMap(8,3), PyramidFaceMap(6,3)

    SAVE Visited, TetraFaceMap, BrickFaceMap, WedgeFaceMap, PyramidFaceMap, &
        FaceMap, faces, elemfamily

!------------------------------------------------------------------------------

    IF(.NOT. Visited ) THEN
      TetraFaceMap(1,:) = (/ 1, 2, 3 /)
      TetraFaceMap(2,:) = (/ 1, 2, 4 /)
      TetraFaceMap(3,:) = (/ 2, 3, 4 /)
      TetraFaceMap(4,:) = (/ 3, 1, 4 /)

      WedgeFaceMap(1,:) = (/ 1, 2, 3 /)
      WedgeFaceMap(2,:) = (/ 4, 5, 6 /)
      WedgeFaceMap(3,:) = (/ 1, 2, 5 /)
      WedgeFaceMap(4,:) = (/ 5, 4, 1 /)
      WedgeFaceMap(5,:) = (/ 3, 2, 5 /)
      WedgeFaceMap(6,:) = (/ 5, 6, 3/)
      WedgeFaceMap(7,:) = (/ 3, 1, 4 /)
      WedgeFaceMap(8,:) = (/ 4, 6, 3/)

      PyramidFaceMap(1,:) = (/ 1, 2, 3 /)
      PyramidFaceMap(2,:) = (/ 3, 4, 1 /)
      PyramidFaceMap(3,:) = (/ 1, 2, 5 /)
      PyramidFaceMap(4,:) = (/ 2, 3, 5 /)
      PyramidFaceMap(5,:) = (/ 3, 4, 5 /)
      PyramidFaceMap(6,:) = (/ 4, 1, 5 /)

      BrickFaceMap(1,:) = (/ 1, 2, 3 /)
      BrickFaceMap(2,:) = (/ 3, 4, 1 /)
      BrickFaceMap(3,:) = (/ 5, 6, 7 /)
      BrickFaceMap(4,:) = (/ 7, 8, 5 /)
      BrickFaceMap(5,:) = (/ 1, 2, 6 /)
      BrickFaceMap(6,:) = (/ 6, 5, 1 /)
      BrickFaceMap(7,:) = (/ 2, 3, 7 /)
      BrickFaceMap(8,:) = (/ 7, 6, 2 /)
      BrickFaceMap(9,:) = (/ 3, 4, 8 /)
      BrickFaceMap(10,:) = (/ 8, 7, 3 /)
      BrickFaceMap(11,:) = (/ 4, 1, 5 /)
      BrickFaceMap(12,:) = (/ 5, 8, 4 /)

      Visited = .TRUE.
    END IF


    IF(face == 1) THEN
      elemfamily = GetElementFamily()
      SELECT CASE( elemfamily )
      CASE(5)
        faces = 4
        FaceMap => TetraFaceMap
      CASE(6)
        faces = 6
        FaceMap => PyramidFaceMap
      CASE(7)
        faces = 8
        FaceMap => WedgeFaceMap
      CASE(8)
        faces = 12
        FaceMap => BrickFaceMap
      CASE DEFAULT
        WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.'
        CALL Fatal('GetLinearTriangleFaces',Message)
      END SELECT
    END IF


    IF(face > faces) THEN
      GotIt = .FALSE.
    ELSE
      GotIt = .TRUE.
      Inds(1:3) = FaceMap(face,1:3)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetLinearTriangleFaces
!------------------------------------------------------------------------------


  FUNCTION LineFaceIntersect(Element,Plane,dim,Line,up,vp,cp) RESULT(Inside)
!---------------------------------------------------------------------------
! This subroutine tests whether the line segment goes through the current
! face of the element. If true the weights and index to the closest node
! are returned.
!---------------------------------------------------------------------------

    TYPE(Nodes_t) :: Plane, Line
    TYPE(Element_t), POINTER   :: Element
    INTEGER :: dim,n
    LOGICAL :: Inside
    REAL (KIND=dp) :: up,vp,cp

    REAL (KIND=dp) :: A(3,3),A0(3,3),B(3),C(3),Eps,Eps2,detA,absA,ds
    LOGICAL :: FiniteLine = .FALSE.

    Inside = .FALSE.

    Eps = 1.0d-8
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

      IF(FiniteLine .AND. ( C(1) < -Eps .OR. C(2) > 1.0d0 + Eps) ) RETURN
      IF(C(2) < -Eps .OR. C(2) > 1.0d0 + Eps) RETURN

      Inside = .TRUE.

      ! Relate the point of intersection to local coordinates
      up = -1.0d0 + 2.0d0 * C(2)
      ! Extent of the line segment
      cp = c(1)

    ELSE IF(DIM == 3) THEN
      A(1,1) = Line % x(2) - Line % x(1)
      A(2,1) = Line % y(2) - Line % y(1)
      A(3,1) = Line % z(2) - Line % z(1)

      A(1,2) = Plane % x(1) - Plane % x(2)
      A(2,2) = Plane % y(1) - Plane % y(2)
      A(3,2) = Plane % z(1) - Plane % z(2)

      A(1,3) = Plane % x(1) - Plane % x(3)
      A(2,3) = Plane % y(1) - Plane % y(3)
      A(3,3) = Plane % z(1) - Plane % z(3)

      ! Check for linearly dependent vectors
      detA = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
           - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
           + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      absA = SUM(ABS(A(1,1:3)))*SUM(ABS(A(2,1:3)))*SUM(ABS(A(3,1:3)))

      IF(ABS(detA) <= eps * absA + Eps2) RETURN

      B(1) = Plane % x(1) - Line % x(1)
      B(2) = Plane % y(1) - Line % y(1)
      B(3) = Plane % z(1) - Line % z(1)

      CALL InvertMatrix( A,3 )
      C(1:3) = MATMUL( A(1:3,1:3),B(1:3) )

      IF (FiniteLine.AND.(C(1)<0.0.OR.C(1)>1.0d0))   RETURN
      IF (ANY(C(2:3)<-Eps).OR.ANY(C(2:3)>1.0d0+Eps)) RETURN
      IF(C(2)+C(3) > 1.0d0 + Eps) RETURN

      Inside = .TRUE.

      ! Relate the point of intersection to local coordinates
      up = C(2)
      vp = C(3)

      ! Extent of the line segment
      cp = C(1)
    END IF
  END FUNCTION LineFaceIntersect

!------------------------------------------------------------------------------
END SUBROUTINE ProjectCalving
!------------------------------------------------------------------------------
