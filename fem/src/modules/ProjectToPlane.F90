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
!> This subroutine may be used to convert field computed in cartesian 3D coordinates on a 2D 
!> axisymmetric or cartesian mesh. There are currently a serial and a parallel version of 
!> this subroutine - this is the serial one. The reason for the fork is that 
!> the parallel version uses techniques that are not optimal in the serial problem. Therefore 
!> the user should herself choose the correct version. Optimally these two approaches should 
!> of course be fused.
!> \ingroup Solvers
!-------------------------------------------------------------------------------------------- 
SUBROUTINE ProjectToPlane( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
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

  TYPE(Solver_t), POINTER :: Solver3D
  TYPE(Element_t), POINTER :: Element, FaceElement
  TYPE(Element_t), TARGET :: TriangleElement
  TYPE(Nodes_t) :: Nodes, LineNodes, FaceNodes, ElementNodes
  TYPE(Variable_t), POINTER :: Variable2D, Variable3D

  REAL(KIND=dp), POINTER :: Values2D(:), Values3D(:), Vals(:)
  REAL(KIND=dp), POINTER :: Basis(:)
  REAL(KIND=dp), POINTER :: PlaneX(:), PlaneY(:), PlaneZ(:)
  REAL(KIND=dp), POINTER :: VolumeX(:), VolumeY(:), VolumeZ(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at, totcpu
#else
  REAL(KIND=dp) :: CPUTime, at, totcpu
#endif
  REAL(KIND=dp) :: x0,y0,z0,xmin,xmax,ymin,ymax,zmin,zmax,r0,rmin,rmax
  REAL(KIND=dp) :: scale, eps, x1,x2,y1,y2,z1,z2,r1,r2
  REAL(KIND=dp) :: xmax3d, xmin3d, ymax3d, ymin3d, zmax3d, zmin3d, rmin3d, rmax3d
  REAL(KIND=dp) :: up,vp,wp,cp,cf,Field1,Field2,MaxRelativeRadius,SqrtElementMetric
  REAL(KIND=dp) :: LocalPoint(3), LocalCoords(3)
  REAL(KIND=dp), POINTER :: MinHeight3D(:), MaxHeight3D(:), MinWidth3D(:), MaxWidth3D(:)
  REAL(KIND=dp), ALLOCATABLE :: IntBasis(:,:), IntExtent(:)

  INTEGER :: MaxInt, Int
  INTEGER, POINTER :: Perm2D(:), Perm3D(:), PlanePerm(:), VolumePerm(:)
  INTEGER :: i,j,k,k2,l,n,t,node
  INTEGER :: PlaneNodes, VolumeElements, Dofs3D, Dofs2D, corners, face, Intersections
  INTEGER :: Loops(8), inds(3), MinimumHits, AxisHits 
  INTEGER, POINTER :: Order2D(:), Order3D(:)
  INTEGER, ALLOCATABLE :: IntOrder(:), IntNodes(:,:)

  CHARACTER(LEN=MAX_NAME_LEN) :: ConvertFromName, ConvertFromVar, EqName
  LOGICAL :: AllocationsDone = .FALSE., stat, GotIt, Rotate, LimitRadius, Found

  SAVE :: Solver3D, Variable2D, Variable3D, Values2D, Values3D, &
      VolumeElements, PlaneNodes, VolumeX, VolumeY, VolumeZ, PlaneX, PlaneY, &
      PlaneZ, MaxInt, IntBasis, IntNodes, IntExtent, IntOrder, &
      MinHeight3D, MaxHeight3D, MinWidth3D, MaxWidth3D, &
      Basis, ElementNodes, LineNodes, FaceNodes, Vals, Eps, Rotate, &
      LimitRadius, MaxRelativeRadius, MinimumHits, AxisHits, Rmax3D


!------------------------------------------------------------------------------

  totcpu = CPUTime()
  
  CALL Info( 'ProjectToPlane', ' ' )
  CALL Info( 'ProjectToPlane', '-----------------------------------' )
  CALL Info( 'ProjectToPlane', ' Projecting 3D solution to 2D ' )
  CALL Info( 'ProjectToPlane', '-----------------------------------' )
  CALL Info( 'ProjectToPlane', ' ' )

  IF ( .NOT. AllocationsDone ) THEN
    IF ( GetLogical( Model % Simulation, 'Output Version Numbers', stat ) ) THEN
      CALL Info( 'ProjectToPlane', 'Version 1.0 by raback (05-02-2007)', LEVEL=4 )
    END IF
    
    ! 3D variables
    !-------------

    ConvertFromName = GetString( Solver % Values, 'Convert From Equation Name', GotIt )
    IF ( .NOT. GotIt )  ConvertFromName = 'induction'
    
    NULLIFY( Solver3D )
    DO i = 1, Model % NumberOfSolvers
      EqName = GetString( Model % Solvers(i) % Values, 'Equation' )
      IF ( TRIM( EqName ) == TRIM( ConvertFromName ) ) THEN
        Solver3D => Model % Solvers(i) 
        EXIT
      END IF
    END DO
    
    IF ( .NOT. ASSOCIATED( Solver3D ) ) THEN
      WRITE( Message, * ) 'Cannot find Solver called ', TRIM( ConvertFromName )
      CALL Error( 'ProjectToPlane', Message )
      CALL Fatal( 'ProjectToPlane','Possibly missing "Convert From Equation Name" field' )
    END IF
    
    ConvertFromVar = GetString( Solver % Values, 'Convert From Variable', GotIt )
    IF ( GotIt ) THEN
      Variable3D => VariableGet( Model % Variables, ConvertFromVar, .TRUE.) 
    ELSE
      Variable3D => Solver3D % Variable 
    END IF

    Values3D => Variable3D % Values
    Perm3D => Variable3D % Perm
    Dofs3D = Variable3D % Dofs
    VolumeElements = Solver3D % NumberOfActiveElements

    VolumeX => Solver3D % Mesh % Nodes % x
    VolumeY => Solver3D % Mesh % Nodes % y      
    VolumeZ => Solver3D % Mesh % Nodes % z      

    VolumePerm => ListGetIntegerArray( Solver % Values,'Volume Permutation',GotIt)
    IF ( gotIt ) THEN
      IF(VolumePerm(1) == 2) VolumeX => Solver3D % Mesh % Nodes % y
      IF(VolumePerm(1) == 3) VolumeX => Solver3D % Mesh % Nodes % z
      IF(VolumePerm(2) == 1) VolumeY => Solver3D % Mesh % Nodes % x
      IF(VolumePerm(2) == 3) VolumeY => Solver3D % Mesh % Nodes % z
      IF(VolumePerm(3) == 1) VolumeZ => Solver3D % Mesh % Nodes % x
      IF(VolumePerm(3) == 2) VolumeZ => Solver3D % Mesh % Nodes % y
    END IF
   
    WRITE( Message, * ) 'Converting from "', TRIM( Variable3D % Name ), &
        '" with ', Dofs3D, 'degrees of freedom'
    CALL Info( 'ProjectToPlane', Message, LEVEL=7 )


    ! 2D variables
    !-------------
    
    PlaneNodes =  Solver % Mesh % NumberOfNodes
    Variable2D => Solver % Variable
    Values2D => Variable2D % Values
    Dofs2D = Variable2D % Dofs
    Perm2D => Variable2D % Perm

    PlaneX => Solver % Mesh % Nodes % x
    PlaneY => Solver % Mesh % Nodes % y      
    PlaneZ => Solver % Mesh % Nodes % z      

    PlanePerm => ListGetIntegerArray( Solver % Values,'Plane Permutation',GotIt)
    IF ( gotIt ) THEN
      IF(PlanePerm(1) == 2) PlaneX => Solver % Mesh % Nodes % y
      IF(PlanePerm(1) == 3) PlaneX => Solver % Mesh % Nodes % z
      IF(PlanePerm(2) == 1) PlaneY => Solver % Mesh % Nodes % x
      IF(PlanePerm(2) == 3) PlaneY => Solver % Mesh % Nodes % z
      IF(PlanePerm(3) == 1) PlaneZ => Solver % Mesh % Nodes % x
      IF(PlanePerm(3) == 2) PlaneZ => Solver % Mesh % Nodes % y
    END IF

    IF(Dofs2D /= Dofs3D) THEN
      CALL Fatal('ProjectToPlane','Number of Dofs must be the same in 2D and 3D')
    END IF

    WRITE( Message, * ) 'Number of mesh nodes in 2D and 3D: ', &
       Solver % Mesh % NumberOfNodes, Solver3D % Mesh % NumberOfNodes
    CALL Info( 'ProjectToPlane', Message, LEVEL=16 )
!------------------------------------------------------------------------------
!   Allocate stuff
!------------------------------------------------------------------------------

    MaxInt =  PlaneNodes
    ALLOCATE(IntBasis(3,MaxInt), IntNodes(3,MaxInt), IntExtent(MaxInt), IntOrder(MaxInt) )
    ALLOCATE( MinHeight3D(VolumeElements), MaxHeight3D(VolumeElements))
    ALLOCATE( MinWidth3D(VolumeElements), MaxWidth3D(VolumeElements))

    n = Solver3D % Mesh % MaxElementNodes
    ALLOCATE(ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), Basis(n) )
    ALLOCATE(LineNodes % x(2), LineNodes % y(2), LineNodes % z(2)  )
    ALLOCATE(FaceNodes % x(3), FaceNodes % y(3), FaceNodes % z(3), Vals(Dofs2D)  )

!   AllocationsDone = .TRUE.

!------------------------------------------------------------------------------
!   Find the scale of the 2D mesh
!------------------------------------------------------------------------------

    xmin = MINVAL( PlaneX )    
    xmax = MAXVAL( PlaneX )
    ymin = MINVAL( PlaneY )    
    ymax = MAXVAL( PlaneY )
    zmin = MINVAL( PlaneZ )    
    zmax = MAXVAL( PlaneZ )

    x0 = xmax - xmin
    y0 = ymax - ymin
    z0 = zmax - zmin

    scale = SQRT(x0*x0 + y0*y0 + z0*z0)    
    LineNodes % y(1) = ymin
    LineNodes % y(2) = LineNodes % y(1) + scale

    Eps = 1.0e-6 * scale    

!------------------------------------------------------------------------------
!  Check control parameters
!------------------------------------------------------------------------------

    Rotate = ListGetLogical(Solver % Values,'Rotate Plane')
    LimitRadius = ListGetLogical(Solver % Values,'Limit Radius',GotIt)
    IF(GotIt) THEN
      MaxRelativeRadius = ListGetConstReal( Solver % Values,'Max Relative Radius',GotIt)
      IF(.NOT. GotIt) MaxRelativeRadius = 0.9999
    END IF

    MinimumHits = ListGetInteger(Solver % Values,'Minimum Hits At Radius',GotIt) 
    IF(.NOT. GotIt) MinimumHits = 1

    AxisHits = ListGetInteger(Solver % Values,'Integration Points At Radius',GotIt) 
    IF(.NOT. GotIt) AxisHits = 2


!------------------------------------------------------------------------------
!   To improve search speed, tabulate values for min and max values of element nodes
!------------------------------------------------------------------------------
   
    DO t = 1, Solver3D % NumberOfActiveElements
      Element => GetActiveElement( t, Solver3D )
      n = GetElementNOFNodes(Element)
      
      ElementNodes % x(1:n) = VolumeX( Element % NodeIndexes(1:n) )
      ElementNodes % y(1:n) = VolumeY( Element % NodeIndexes(1:n) )
      ElementNodes % z(1:n) = VolumeZ( Element % NodeIndexes(1:n) )
      
      MinHeight3D(t) = MINVAL(ElementNodes % z(1:n))
      MaxHeight3D(t) = MAXVAL(ElementNodes % z(1:n))
      
      IF(Rotate) THEN
        MinWidth3D(t) = SQRT( MINVAL(ElementNodes % x(1:n)**2 + ElementNodes % y(1:n)**2) )
        MaxWidth3D(t) = SQRT( MAXVAL(ElementNodes % x(1:n)**2 + ElementNodes % y(1:n)**2) )
      ELSE
        MinWidth3D(t) = MINVAL(ElementNodes % x(1:n))
        MaxWidth3D(t) = MAXVAL(ElementNodes % x(1:n))                  
      END  IF
    END DO
    rmax3d = MAXVAL(MaxWidth3D)
 
    AllocationsDone = .TRUE.
  END IF


  ! By contruction split everything into triangles        
  corners = 3
  
  TriangleElement % TYPE => GetElementType( 303, .FALSE. )

  FaceElement => TriangleElement
  Loops = 0

  !   Loop over 2D nodes:
  !   -------------------
  DO node = 1, PlaneNodes

    Loops(1) = Loops(1) + 1
 
    IF(Perm2D(node) == 0) CYCLE

    x0 = PlaneX(node)
    y0 = PlaneY(node)
    z0 = PlaneZ(node)

    IF(LimitRadius) THEN
      x0 = MIN( x0, MaxRelativeRadius * rmax3d )
    END IF
    LineNodes % x(1) = x0
    LineNodes % x(2) = x0
    LineNodes % z(1) = z0
    LineNodes % z(2) = z0

    !   Loop over 3D elements:
    !   -------------------
    Int = 0
    Vals = 0
    

    DO t = 1, VolumeElements
      
      Loops(2) = Loops(2) + 1

      IF(z0 > MaxHeight3D(t) + Eps) CYCLE
      IF(z0 < MinHeight3D(t) - Eps) CYCLE

      IF(x0 > MaxWidth3D(t) + Eps) CYCLE
      IF(x0 < MinWidth3D(t) - Eps) CYCLE

      Loops(3) = Loops(3) + 1

!------------------------------------------------------------------------------
      Element => GetActiveElement( t, Solver3D )
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

        IF(.NOT. Rotate ) THEN
          FaceNodes % x(1:corners) = ElementNodes % x(inds(1:corners))
          FaceNodes % y(1:corners) = ElementNodes % y(inds(1:corners))
          FaceNodes % z(1:corners) = ElementNodes % z(inds(1:corners))
        ELSE
          FaceNodes % x(1:corners) = SQRT( &
              ElementNodes % x(inds(1:corners))**2 + &
              ElementNodes % y(inds(1:corners))**2 )                       
          FaceNodes % y(1:corners) = 0.0d0            
          FaceNodes % z(1:corners) = ElementNodes % z(inds(1:corners))                       
        END IF

        xmin = MINVAL( FaceNodes % x(1:corners) )
        xmax = MAXVAL( FaceNodes % x(1:corners) )                  
        IF( x0 < xmin - Eps) CYCLE
        IF( x0 > xmax + Eps) CYCLE
        
        zmin = MINVAL( FaceNodes % z(1:corners) )
        zmax = MAXVAL( FaceNodes % z(1:corners) )                     
        IF( z0 < zmin - Eps) CYCLE
        IF( z0 > zmax + Eps) CYCLE

        CALL LineFaceIntersect(Element,FaceNodes,corners,LineNodes,Found,up,vp,cp)
        IF(.NOT. Found) CYCLE

        Loops(5) = Loops(5) + 1

        Intersections = Intersections + 1
        Int = Int + 1
        IF(Int > MaxInt) CALL AllocateMoreSpace()

        CALL NodalBasisFunctions(Corners, Basis, FaceElement, up, vp, 0.0d0)

        IF(Rotate) THEN
          x1 = SUM(Basis(1:corners) * ElementNodes % x(inds(1:corners)))
          y1 = SUM(Basis(1:corners) * ElementNodes % y(inds(1:corners)))
          cp = ATAN2(x1, y1)
        END IF

        IntExtent(Int) = cp
        IntBasis(1:corners,Int) = Basis(1:corners)
        IntNodes(1:corners,Int) = Element % NodeIndexes( inds(1:corners) ) 
      END DO

      IF( Intersections /= 0 .AND. Intersections /= 2 ) THEN
        Loops(6) = Loops(6) + 1
      END IF
      
    END DO

    
    IF(Int > MinimumHits) THEN
      
      Loops(7) = Loops(7) + 1

      IF( Int == 1 ) THEN
        DO l=1,Dofs2D
          Vals(l) = SUM( IntBasis(:,k) * Values3D( Dofs3D * (Perm3D(IntNodes(:,1)) - 1) + l) )
        END DO
      ELSE
      
        DO i=1,Int 
          IntOrder(i) = i
        END DO
        
        CALL SortR( Int, IntOrder, IntExtent)
        
        Vals = 0.0d0
        DO j=1,Int-1
          k = IntOrder(j)
          k2 = IntOrder(j+1)
          DO l=1,Dofs2D
            Field1 = SUM( IntBasis(:,k) * Values3D( Dofs3D * (Perm3D(IntNodes(:,k))-1) + l) )      
            Field2 = SUM( IntBasis(:,k2) * Values3D( Dofs3D*(Perm3D(IntNodes(:,k2))-1) + l) )     
            Vals(l) = Vals(l) + 0.5d0 * (Field1 + Field2) * (IntExtent(j) - IntExtent(j+1))
          END DO
        END DO
        
        IF(Rotate) THEN
          cf = 2*PI -  (IntExtent(1) - IntExtent(Int))
          k = IntOrder(Int)
          k2 = IntOrder(1)      
          DO l=1,Dofs2D         
            Field1 = SUM( IntBasis(:,k) * Values3D( Dofs3D*(Perm3D(IntNodes(:,k))-1) + l) )     
            Field2 = SUM( IntBasis(:,k2) * Values3D( Dofs3D*(Perm3D(IntNodes(:,k2))-1) + l) )      
            Vals(l) = Vals(l) + 0.5d0 * (Field1 + Field2) * cf
          END DO
          Vals = Vals / (2.0_dp*PI)
        ELSE
          Vals = Vals / (IntExtent(1)-IntExtent(Int)) 
        END IF
      END IF
    ELSE 
      
      LocalPoint(1) = x0
      LocalPoint(2) = y0
      LocalPoint(3) = z0
      GotIt = .FALSE.
      
      ! Take symmetric hits if not exactly on the axis
      IF( Rotate .AND. ABS(x0) > Eps) THEN
        Loops(8) = Loops(8) + 1
        k2 = AxisHits
      ELSE
        k2 = 1
      END IF
      
      j = 0
      Vals = 0.0d0

      DO k=1,k2
        GotIt = .FALSE.

        IF(k > 1) THEN
          LocalPoint(1) = COS( 2*PI*(k-1.0d0)/AxisHits )
          LocalPoint(2) = SIN( 2*PI*(k-1.0d0)/AxisHits )
        END IF

        DO t = 1, VolumeElements

          IF(z0 > MaxHeight3D(t) + Eps) CYCLE
          IF(z0 < MinHeight3D(t) - Eps) CYCLE

          Element => GetActiveElement( t, Solver3D )
          n = GetElementNOFNodes(Element)

          ElementNodes % x(1:n) = VolumeX( Element % NodeIndexes(1:n) )
          ElementNodes % y(1:n) = VolumeY( Element % NodeIndexes(1:n) )
          ElementNodes % z(1:n) = VolumeZ( Element % NodeIndexes(1:n) )
        
          IF ( PointInElement( Element, ElementNodes, LocalPoint, LocalCoords ) ) THEN
            GotIt = .TRUE.
            EXIT
          END IF
        END DO
        
        IF(GotIt) THEN
          j = j + 1
          up = LocalCoords(1)
          vp = LocalCoords(2)
          wp = LocalCoords(3)
          
          stat = ElementInfo( Element,ElementNodes,up,vp,wp,SqrtElementMetric,Basis)
          
          DO l=1,Dofs2D
            Field1 = SUM( Basis(1:n) * Values3D( Dofs3D * ( Perm3D(Element % NodeIndexes(1:n)) - 1) + l) )     
            Vals(l) = Vals(l) + Field1
          END DO
        END IF
      END DO

      IF( j > 0 ) THEN
        Vals = Vals / j
      ELSE
        Vals = 0
        PRINT *,'No hits for this node',node,x0,y0,z0
      END IF

    END IF

    DO l=1,Dofs2D
      Values2D(Dofs2D*( Perm2D(node) - 1) + l) = Vals(l)
    END DO

  END DO


  Solver % Variable % Norm = SQRT( SUM(Values2D**2) / SIZE(Perm2D) )
  

  WRITE( Message, * ) 'Basic search loops: ',Loops(2:5)
  CALL Info( 'ProjectToPlane', Message, LEVEL=4 )

  IF(Loops(5) > 0) THEN    
    WRITE( Message, * ) 'Special cases loops: ',Loops(6:8)
    CALL Info( 'ProjectToPlane', Message, LEVEL=4 )
  END IF

  WRITE( Message, * ) 'Total CPU time used: ', CPUTime() - totcpu
  CALL INFO( 'ProjectToPlane', Message, LEVEL=8 )
  CALL Info( 'ProjectToPlane', ' ' )
  CALL Info( 'ProjectToPlane', 'All done' )
  CALL Info( 'ProjectToPlane', ' ' )

  IF(.FALSE.) THEN
    DEALLOCATE(IntBasis, IntNodes, IntExtent, IntOrder, &
        MinHeight3D, MaxHeight3D, MinWidth3D, MaxWidth3D)
  END IF


CONTAINS
  

  SUBROUTINE AllocateMoreSpace()
   
    REAL(KIND=dp), ALLOCATABLE :: TmpBasis(:,:), TmpExtent(:)
    INTEGER, ALLOCATABLE :: TmpOrder(:), TmpNodes(:,:)
    INTEGER :: OldMaxInt
    
    OldMaxInt = MaxInt
    
    ALLOCATE(TmpBasis(3,MaxInt), TmpNodes(3,MaxInt), TmpExtent(MaxInt), TmpOrder(MaxInt) )
    
    TmpBasis(:,1:OldMaxInt) = IntBasis(:,1:OldMaxInt)
    TmpNodes(:,1:OldMaxInt) = IntNodes(:,1:OldMaxInt) 
    TmpExtent(1:OldMaxInt) = IntExtent(1:OldMaxInt)
    TmpOrder(1:OldMaxInt) = IntOrder(1:OldMaxInt)
    
    MaxInt = MaxInt + PlaneNodes
    DEALLOCATE(IntBasis, IntNodes, IntExtent, IntOrder)
    ALLOCATE(IntBasis(3,MaxInt), IntNodes(3,MaxInt), IntExtent(MaxInt), IntOrder(MaxInt) )
    
    IntBasis(:,1:OldMaxInt) = TmpBasis(:,1:OldMaxInt)
    IntNodes(:,1:OldMaxInt) = TmpNodes(:,1:OldMaxInt) 
    IntExtent(1:OldMaxInt) = TmpExtent(1:OldMaxInt)
    IntOrder(1:OldMaxInt) = TmpOrder(1:OldMaxInt)
    
    DEALLOCATE(TmpBasis, TmpNodes, TmpExtent, TmpOrder)
    ! PRINT *,'Allocated more space',MaxInt
    
  END SUBROUTINE AllocateMoreSpace


!------------------------------------------------------------------------------
! 3D mesh faces.
!------------------------------------------------------------------------------
  SUBROUTINE GetLinearTriangleFaces( Element, face, inds, GotIt )
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: face, inds(:)
    LOGICAL :: GotIt
!------------------------------------------------------------------------------
    LOGICAL :: Visited
    INTEGER :: faces, elemfamily     
    INTEGER, POINTER :: FaceMap(:,:)
    INTEGER, TARGET  :: TetraFaceMap(4,3), BrickFaceMap(12,3), WedgeFaceMap(8,3), PyramidFaceMap(6,3)    

    SAVE Visited, TetraFaceMap, BrickFaceMap, WedgeFaceMap, PyramidFaceMap, FaceMap, faces

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
      elemfamily = Element % TYPE % ElementCode / 100
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
        CALL Fatal('FindMeshFaces',Message)
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


  SUBROUTINE LineFaceIntersect(Element,Plane,dim,Line,Inside,up,vp,cp)
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
      absA = SUM(ABS(A(1,1:3))) * SUM(ABS(A(2,1:3))) * SUM(ABS(A(3,1:3))) 
      
      IF(ABS(detA) <= eps * absA + Eps2) RETURN

      B(1) = Plane % x(1) - Line % x(1)
      B(2) = Plane % y(1) - Line % y(1)
      B(3) = Plane % z(1) - Line % z(1)

      CALL InvertMatrix( A,3 )

      C(1:3) = MATMUL( A(1:3,1:3),B(1:3) )

      IF( FiniteLine .AND. ( C(1) < 0.0 .OR. C(1) > 1.0d0 ) ) RETURN
      
      IF( ANY(C(2:3) < -Eps) .OR. ANY(C(2:3) > 1.0d0 + Eps) ) RETURN
      IF(C(2)+C(3) > 1.0d0 + Eps) RETURN
      
      Inside = .TRUE. 

      ! Relate the point of intersection to local coordinates
      up = C(2)
      vp = C(3)
      
      ! Extent of the line segment 
      cp = c(1)      
    END IF

  END SUBROUTINE LineFaceIntersect
  
!------------------------------------------------------------------------------
END SUBROUTINE ProjectToPlane
!------------------------------------------------------------------------------


SUBROUTINE ParallelProjectToPlane_init( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE Lists

  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  CHARACTER(LEN=MAX_NAME_LEN) :: Name, Str
  LOGICAL :: Found
  INTEGER :: i

  IF(.NOT. ListCheckPresent(Solver % Values,'Variable')) THEN
    Name = ListGetString(Solver % Values,'Variable 1',Found)
    IF( Found ) CALL ListAddString( Solver % Values,'Variable','ave_'//TRIM(Name))
  END IF
  
  DO i=2,9
    WRITE( Str,'(A,I2)') 'Variable',i
    Name = ListGetString( Solver % Values,Str,Found)
    IF( Found ) THEN
      WRITE( Str,'(A,I2)') 'Exported Variable',i-1
      CALL ListAddString( Solver % Values,TRIM(Str),'ave_'//TRIM(Name))
    ELSE
      EXIT
    END IF
  END DO

END SUBROUTINE ParallelProjectToPlane_init


!-------------------------------------------------------------------------------------------- 
!> This subroutine may be used to convert field computed in cartesian 3D coordinates on a 2D 
!> axisymmetric or cartesian mesh. There are currently a serial and a parallel version of 
!> this subroutine - this is the parallel one. 
!> \ingroup Solvers
!-------------------------------------------------------------------------------------------- 
SUBROUTINE ParallelProjectToPlane( Model,Solver,dt,TransientSimulation )
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
  TYPE(Variable_t), POINTER :: Variable2D, Variable3D, Var
  TYPE(Mesh_t), POINTER :: Mesh2D, Mesh3D, Mesh

  REAL(KIND=dp), POINTER :: Vals(:), Vals0(:)
  REAL(KIND=dp), POINTER :: Basis(:)
  REAL(KIND=dp), POINTER :: PlaneX(:), PlaneY(:), PlaneZ(:)
  REAL(KIND=dp), POINTER :: VolumeX(:), VolumeY(:), VolumeZ(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at, totcpu
#else
  REAL(KIND=dp) :: CPUTime, at, totcpu
#endif
  REAL(KIND=dp) :: x0,y0,z0,xmin,xmax,ymin,ymax,zmin,zmax,r0,rmin,rmax
  REAL(KIND=dp) :: scale, eps, x1,x2,y1,y2,z1,z2,r1,r2
  REAL(KIND=dp) :: xmax3d, xmin3d, ymax3d, ymin3d, zmax3d, zmin3d, rmin3d, rmax3d
  REAL(KIND=dp) :: up,vp,wp,cp,cf,Field1,Field2,MaxRelativeRadius,detJ
  REAL(KIND=dp) :: LocalPoint(3), LocalCoords(3)
  REAL(KIND=dp), POINTER :: MinHeight3D(:), MaxHeight3D(:), MinWidth3D(:), MaxWidth3D(:)
  REAL(KIND=dp), POINTER :: IntValues(:,:), IntBasis(:,:), IntExtent(:)

  INTEGER :: MaxInt, Int
  INTEGER, POINTER :: Perm2D(:), Perm3D(:), PlanePerm(:), VolumePerm(:)
  INTEGER :: i,j,k,k2,l,n,t,lnode,node
  INTEGER :: PlaneNodes, VolumeElements, Dofs, corners, face, Intersections
  INTEGER :: Loops(8), inds(3), MinimumHits, AxisHits 
  INTEGER, POINTER :: Order2D(:), Order3D(:)
  INTEGER, POINTER :: IntOrder(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: ConvertFromName, ConvertFromVar, Name, Str
  LOGICAL :: AllocationsDone = .FALSE., stat, GotIt, Rotate, LimitRadius, Found

  INTEGER :: pe, myNodes, myStart, peNodes, peNodesn, peStart, peEnd, &
          totcount, curr, ierr, status(MPI_STATUS_SIZE), nn

  TYPE ValueTable_t
    REAL(KIND=dp), POINTER :: Values(:) 
  END TYPE ValueTable_t
  TYPE(ValueTable_t) :: ValueTable3D(10), ValueTable2D(10)

  TYPE PointStore_t
    INTEGER :: Int
    REAL(KIND=dp),POINTER :: IntExtent(:), IntValues(:,:)
  END TYPE PointStore_t
  TYPE(PointStore_t), ALLOCATABLE :: PointStore(:)

  INTEGER, ALLOCATABLE :: cm_int(:), cm_int0(:), icount(:), icount0(:)
  REAL(KIND=dp), ALLOCATABLE :: cm_values(:), cm_values0(:), cm_extent(:)

  SAVE :: Variable2D, Variable3D, &
      VolumeElements, PlaneNodes, VolumeX, VolumeY, VolumeZ, PlaneX, PlaneY, &
      PlaneZ, MaxInt, IntBasis, IntValues, IntExtent, IntOrder, &
      MinHeight3D, MaxHeight3D, MinWidth3D, MaxWidth3D, &
      Basis, ElementNodes, LineNodes, FaceNodes, Vals, Eps, Rotate, &
      LimitRadius, MaxRelativeRadius, MinimumHits, AxisHits, Rmax3D, &
      Mesh2D, Mesh3D, Perm2D, Perm3D, Dofs, ValueTable2D, &
      ValueTable3D


!------------------------------------------------------------------------------

  totcpu = CPUTime()
  
  CALL Info( 'ProjectToPlane', ' ' )
  CALL Info( 'ProjectToPlane', '-----------------------------------' )
  CALL Info( 'ProjectToPlane', ' Projecting 3D solution to 2D ' )
  CALL Info( 'ProjectToPlane', '-----------------------------------' )
  CALL Info( 'ProjectToPlane', ' ' )
  
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
    IF ( GetLogical( Model % Simulation, 'Output Version Numbers', stat ) ) THEN
      CALL Info( 'ProjectToPlane', 'Version 1.0 by raback (05-02-2007)', LEVEL=4 )
    END IF
    
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
        CALL Info('ProjectToPlane','3D mesh: '//TRIM(Mesh % Name) )
        Mesh3D => Mesh
      ELSE IF( Mesh % MeshDim == 2 ) THEN
        CALL Info('ProjectToPlane','2D mesh: '//TRIM(Mesh % Name) )
        Mesh2D => Mesh
      END IF
      Mesh => Mesh % next
    END DO
    
    IF(.NOT. ASSOCIATED(Mesh2D)) THEN
      CALL Fatal('ProjectToPlane','Did not find 2D mesh')
    END IF
    IF(.NOT. ASSOCIATED(Mesh3D)) THEN
      CALL Fatal('ProjectToPlane','Did not find 3D mesh')
    END IF
    
    IF(.NOT. ASSOCIATED(Mesh2D, Solver % Mesh)) THEN
      CALL Fatal('ProjectToPlane','2D mesh not same as Solver mesh')
    END IF


    !------------------------------------------------------------------------------
    ! Make the list of 2D and 3D variables 
    !------------------------------------------------------------------------------
    DO i=1,9
      WRITE( Str,'(A,I2)') 'Variable',i      
      Name = ListGetString( Solver % Values,Str,Found)
      IF(.NOT. Found) EXIT
      
      Var => VariableGet( Mesh3D % Variables, TRIM(Name) )
      IF( .NOT. ASSOCIATED(Var)) EXIT

      IF( Var % DOFs > 1 ) THEN
        CALL Fatal('ProjectToPlane','Only scalar variables may be mapped')
      END IF
      IF( i == 1 ) Variable3D => Var
      ValueTable3D(i) % Values => Var % Values
      
      IF( i == 1 ) THEN
        Var => Solver % Variable
        IF(.NOT. ASSOCIATED(Var)) &
            Var => VariableGet( Mesh2D % Variables, 'ave_'//TRIM(Name) )
      ELSE
        Var => VariableGet( Mesh2D % Variables, 'ave_'//TRIM(Name) )        
      END IF

      IF( .NOT. ASSOCIATED(Var)) EXIT

      IF( Var % DOFs > 1 ) THEN
        CALL Fatal('ProjectToPlane','Only scalar variables may be mapped')
      END IF
      IF( i == 1 ) Variable2D => Var
      ValueTable2D(i) % Values => Var % Values
    END DO
    Dofs = i-1
    IF( Dofs == 0 ) THEN
      CALL Fatal('ProjectToPlane','No variables to work with')
    END IF

    !------------------------------------------------------------------------------
    ! The permutation are assumed to be constant 
    !------------------------------------------------------------------------------
    Perm3D => Variable3D % Perm
    Perm2D => Variable2D % Perm
  
    PlaneNodes = COUNT(Perm2D>0)    
    VolumeElements = Mesh3D % NumberOfBulkElements

    WRITE( Message, * ) 'Number of mesh nodes in 2D and 3D: ', &
        Mesh2d % NumberOfNodes, Mesh3D % NumberOfNodes
    CALL Info( 'ProjectToPlane', Message, LEVEL=16 )

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
    ALLOCATE(FaceNodes % x(3), FaceNodes % y(3), FaceNodes % z(3), Vals(Dofs), Vals0(Dofs))

!------------------------------------------------------------------------------
!   Find the scale of the 2D mesh
!------------------------------------------------------------------------------

    xmin = MINVAL( PlaneX )    
    xmax = MAXVAL( PlaneX )
    ymin = MINVAL( PlaneY )    
    ymax = MAXVAL( PlaneY )
    zmin = MINVAL( PlaneZ )    
    zmax = MAXVAL( PlaneZ )

    x0 = xmax - xmin
    y0 = ymax - ymin
    z0 = zmax - zmin

    scale = SQRT(x0*x0 + y0*y0 + z0*z0)    
    LineNodes % y(1) = ymin
    LineNodes % y(2) = LineNodes % y(1) + scale

    Eps = 1.0d-6 * scale    

!------------------------------------------------------------------------------
!  Check control parameters
!------------------------------------------------------------------------------

    Rotate = ListGetLogical(Solver % Values,'Rotate Plane')
    LimitRadius = ListGetLogical(Solver % Values,'Limit Radius',GotIt)
    IF(GotIt) THEN
      MaxRelativeRadius = ListGetConstReal( Solver % Values,'Max Relative Radius',GotIt)
      IF(.NOT. GotIt) MaxRelativeRadius = 0.9999
    END IF

    MinimumHits = ListGetInteger(Solver % Values,'Minimum Hits At Radius',GotIt) 
    IF(.NOT. GotIt) MinimumHits = 1

    AxisHits = ListGetInteger(Solver % Values,'Integration Points At Radius',GotIt) 
    IF(.NOT. GotIt) AxisHits = 2

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
      
      IF(Rotate) THEN
        MinWidth3D(t) = SQRT( MINVAL(ElementNodes % x(1:n)**2 + &
            ElementNodes % y(1:n)**2) )
        MaxWidth3D(t) = SQRT( MAXVAL(ElementNodes % x(1:n)**2 + &
            ElementNodes % y(1:n)**2) )
      ELSE
        MinWidth3D(t) = MINVAL(ElementNodes % x(1:n))
        MaxWidth3D(t) = MAXVAL(ElementNodes % x(1:n))
      END  IF
    END DO
    rmax3d = ParallelReduction(MAXVAL(MaxWidth3D),2)
    
    AllocationsDone = .TRUE.
  END IF

  !
  ! allocate space for 3d face intersections / 2d point:
  ! ----------------------------------------------------
  MaxInt =  100
  ALLOCATE( PointStore(PlaneNodes) )
  DO i=1,PlaneNodes
    ALLOCATE( PointStore(i) % IntValues(Dofs,MaxInt), &
              PointStore(i) % IntExtent(MaxInt) )
    PointStore(i) % Int = 0
    PointStore(i) % IntExtent = 0
    PointStore(i) % IntValues = 0
  END DO


  ! By contruction split everything into triangles        
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

    IF(LimitRadius) THEN
      x0 = MIN( x0, MaxRelativeRadius * rmax3d )
    END IF
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

        IF(.NOT. Rotate) THEN
          FaceNodes % x(1:corners) = ElementNodes % x(inds(1:corners))
          FaceNodes % y(1:corners) = ElementNodes % y(inds(1:corners))
          FaceNodes % z(1:corners) = ElementNodes % z(inds(1:corners))
        ELSE
          FaceNodes % x(1:corners) = SQRT( &
              ElementNodes % x(inds(1:corners))**2 + &
              ElementNodes % y(inds(1:corners))**2 )                       
          FaceNodes % y(1:corners) = 0.0d0            
          FaceNodes % z(1:corners) = ElementNodes % z(inds(1:corners))                       
        END IF

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

        IF(Rotate) THEN
          x1 = SUM(Basis(1:corners)*ElementNodes % x(inds(1:corners)))
          y1 = SUM(Basis(1:corners)*ElementNodes % y(inds(1:corners)))
          cp = ATAN2(x1, y1)
        END IF

        !
        ! store extent of the line + value(s) at point of
        ! intersection:
        ! ------------------------------------------------
        PointStore(lnode) % Int = PointStore(lnode) % Int + 1
        curr=PointStore(lnode) % Int
        IF (curr>SIZE(PointStore(lnode) % IntExtent)) &
            CALL AllocateMoreSpace(PointStore(lnode), MaxInt)

        PointStore(lnode) % IntExtent(curr) = cp

        DO l=1,Dofs
          PointStore(lnode) % IntValues(l,curr) = &
              SUM(Basis(1:corners)* ValueTable3D(l) % &
              Values(Perm3D( Element % NodeIndexes(inds(1:corners)) )) )
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
    CALL SendPoints()
    CALL ReceivePoints()
  END IF

  ! compute integrals from intersections:
  ! -------------------------------------
  CALL ComputeIntegrals()

  ! try a different way where the above is not successful:
  ! -----------------------------------------------------
  CALL ComputeMissingIntegrals()

  ! send all results to pe 0:
  ! -------------------------
  IF (ParEnv % PEs>0) CALL GatherResults()

  DO i=1,PlaneNodes
    DEALLOCATE(PointStore(i) % IntValues, PointStore(i) % IntExtent)
  END DO
  DEALLOCATE( PointStore )
  
!  Solver % Variable % Norm = SQRT( SUM(Values2D**2) / SIZE(Perm2D) )
  
  WRITE( Message,'(A,4I8)' ) 'Basic search loops: ',Loops(2:5)
  CALL Info( 'ProjectToPlane', Message, LEVEL=4 )
  
  IF(Loops(5) > 0) THEN    
    WRITE( Message,'(A,3I8)' ) 'Special cases loops: ',Loops(6:8)
    CALL Info( 'ProjectToPlane', Message, LEVEL=4 )
  END IF
  
  WRITE( Message,'(A,F8.2)' ) 'Total CPU time used: ', CPUTime() - totcpu
  CALL INFO( 'ProjectToPlane', Message, LEVEL=8 )
  CALL Info( 'ProjectToPlane', ' ' )
  CALL Info( 'ProjectToPlane', 'All done' )
  CALL Info( 'ProjectToPlane', ' ' )

CONTAINS
  
!------------------------------------------------------------------------------
  SUBROUTINE SendPoints()
!------------------------------------------------------------------------------
    !
    ! send intersections to owners of the points:
    ! -------------------------------------------
    CALL CheckBuffer(104857600)
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
              cm_values(Dofs*totcount) )

        j=0; totcount=0
        DO lnode=peStart,peEnd
          j = j + 1
          cm_int(j) = PointStore(lnode) % Int
          DO k=1,cm_int(j)
            totcount=totcount+1
            int=Dofs*(totcount-1)
            cm_values(int+1:int+Dofs) = &
               PointStore(lnode) % IntValues(1:Dofs,k)
            cm_extent(totcount) = PointStore(lnode) % IntExtent(k)
           END DO
        END DO  
      END IF
      CALL MPI_BSEND( cm_int, peNodes, MPI_INTEGER, pe, &
                100, ELMER_COMM_WORLD, ierr )
      IF (totcount>0 ) THEN
        CALL MPI_BSEND( cm_values, Dofs*totcount, &
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
        ALLOCATE( cm_values(Dofs*totcount), cm_extent(totcount) )

        CALL MPI_RECV( cm_values, Dofs*totcount, &
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
            int = Dofs*(totcount-1)
            PointStore(lnode) % IntExtent(curr) = cm_extent(totcount)
            PointStore(lnode) % IntValues(:,curr) = cm_values(int+1:int+Dofs)
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
  SUBROUTINE ComputeIntegrals()
!------------------------------------------------------------------------------

    DO node=1,Mesh2d % NumberOfNodes
      lnode = Perm2D(node)
      IF(lnode<myStart .OR. lnode>myStart+myNodes-1) CYCLE

      Int = PointStore(lnode) % Int
      IF (Int<=MinimumHits) CYCLE

      IntValues => PointStore(lnode) % IntValues
      IntExtent => PointStore(lnode) % IntExtent

      Loops(7) = Loops(7) + 1

      IF( Int == 1 ) THEN
        DO i=1,Dofs
          Vals(l) = IntValues(l,1)
        END DO
      ELSE
        ALLOCATE(IntOrder(Int)); IntOrder = [(i,i=1,Int)]
        CALL SortR(Int, IntOrder, IntExtent)
      
        Vals = 0.0d0
        DO j=1,Int-1
          k = IntOrder(j)
          k2 = IntOrder(j+1)
          DO l=1,Dofs
            Field1 = IntValues(l,k)
            Field2 = IntValues(l,k2)
            Vals(l) = Vals(l)+0.5d0*(Field1+Field2)*(IntExtent(j)-IntExtent(j+1))
          END DO
        END DO
      
        IF(Rotate) THEN
          cf = 2*PI-(IntExtent(1) - IntExtent(Int))
          k = IntOrder(Int)
          k2 = IntOrder(1)      
          DO l=1,Dofs         
            Field1 = IntValues(l,k)
            Field2 = IntValues(l,k2)
            Vals(l) = Vals(l) + 0.5d0 * (Field1 + Field2) * cf
          END DO
          Vals = Vals / (2.0_dp*PI)
        ELSE
          Vals = Vals / (IntExtent(1)-IntExtent(Int)) 
        END IF

        DEALLOCATE(IntOrder)
      END IF

      DO l=1,Dofs
        ValueTable2D(l) % Values(lnode) = Vals(l)
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeIntegrals
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ComputeMissingIntegrals()
!------------------------------------------------------------------------------
    !
    ! check for points, where not enough intersections
    ! were found; try instead direct evaluation with
    ! equally spaced points in 3d:
    ! -------------------------------------------------
    IF ( ParEnv % PEs>1 ) THEN
      ALLOCATE(cm_int(PlaneNodes),cm_int0(PlaneNodes))
      cm_int0 = 0
      DO i=myStart,myStart+myNodes-1
        cm_int0(i) = PointStore(i) % Int
      END DO
      CALL MPI_ALLREDUCE( cm_int0,cm_int,PlaneNodes, &
        MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr )

      nn=COUNT(cm_int<=MinimumHits)
      ALLOCATE(cm_values0(dofs*nn), cm_values(dofs*nn), &
               icount0(nn), icount(nn))
    END IF 

    nn = 0
    DO node=1,Mesh2d % NumberOfNodes
      lnode = Perm2D(node)
      IF ( lnode<=0 ) CYCLE

      IF ( ParEnv % Pes>1 ) THEN
        Int = cm_int(lnode)
      ELSE
        Int = PointStore(lnode) % Int
      END IF
      IF ( int>MinimumHits ) CYCLE

      x0 = PlaneX(node)
      y0 = PlaneY(node)
      z0 = PlaneZ(node)

      IF(LimitRadius) &
        x0 = MIN(x0,MaxRelativeRadius * rmax3d)

      LocalPoint(1) = x0
      LocalPoint(2) = y0
      LocalPoint(3) = z0
    
      ! Take symmetric hits if not exactly on the axis
      IF( Rotate.AND.ABS(x0)>Eps) THEN
        Loops(8) = Loops(8) + 1
        k2 = AxisHits
      ELSE
        k2 = 1
      END IF
    
      j = 0
      Vals = 0.0d0

      DO k=1,k2
        LocalPoint(1) = x0*COS(2*PI*(k-1._dp)/AxisHits)
        LocalPoint(2) = x0*SIN(2*PI*(k-1._dp)/AxisHits)
        x1 = LocalPoint(1)
        y1 = LocalPoint(2)

        DO t = 1, VolumeElements
          IF(z0 > MaxHeight3D(t) + Eps) CYCLE
          IF(z0 < MinHeight3D(t) - Eps) CYCLE

          Element => Mesh3D % Elements(t)
          Model % CurrentElement => Element
          n = GetElementNOFNodes(Element)

          ElementNodes % x(1:n) = VolumeX(Element % NodeIndexes(1:n))
          IF( MINVAL(ElementNodes % x(1:n))-Eps>x1 ) CYCLE
          IF( MAXVAL(ElementNodes % x(1:n))+Eps<x1 ) CYCLE

          ElementNodes % y(1:n) = VolumeY(Element % NodeIndexes(1:n))
          IF( MINVAL(ElementNodes % y(1:n))-Eps>y1 ) CYCLE
          IF( MAXVAL(ElementNodes % y(1:n))+Eps<y1 ) CYCLE

          ElementNodes % z(1:n) = VolumeZ( Element % NodeIndexes(1:n) )
    
          IF (PointInElement(Element,ElementNodes,LocalPoint,LocalCoords)) THEN
            j = j + 1
            up = LocalCoords(1)
            vp = LocalCoords(2)
            wp = LocalCoords(3)
            stat = ElementInfo(Element,ElementNodes,up,vp,wp,detJ,Basis)
            DO l=1,Dofs
              Field1 = SUM(Basis(1:n)*ValueTable3D(l) % Values( &
                    Perm3D(Element % NodeIndexes(1:n)) ) )
              Vals(l) = Vals(l) + Field1
            END DO
            EXIT
          END IF
        END DO
      END DO

      IF (ParEnv % PEs>1) THEN
        nn=nn+1
        icount0(nn)=j
        DO l=1,Dofs
          cm_values0(Dofs*(nn-1)+l) = Vals(l)
        END DO
      ELSE
        DO l=1,Dofs
          IF ( j>0 ) THEN
            ValueTable2D(l) % Values(lnode) = Vals(l)/j
          ELSE
            ValueTable2D(l) % Values(lnode) = 0._dp
          END IF
        END DO
      END IF
    END DO

    IF ( ParEnv % PEs > 1 ) THEN
      IF ( nn>0 ) THEN
        CALL MPI_ALLREDUCE( icount0, icount, nn, &
             MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr )

        CALL MPI_ALLREDUCE( cm_values0, cm_values, Dofs*nn, &
             MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr )

        nn = 0
        DO node=1,Mesh2d % NumberOfNodes
          lnode = Perm2D(node)
          IF ( lnode<=0 ) CYCLE

          IF ( cm_int(lnode)<=MinimumHits ) THEN
            nn=nn+1
            IF ( icount(nn)>0 ) THEN
              DO l=1,Dofs
                ValueTable2D(l) % Values(lnode) = &
                    cm_values(Dofs*(nn-1)+l)/icount(nn)
              END DO
            ELSE
              IF ( Parenv % myPe==0 ) PRINT*,lnode, 'No hits: ', &
                  node,planex(node),planey(node),planez(node)
            END IF
          END IF
        END DO
      END IF
      DEALLOCATE( cm_values0, cm_values, icount0, icount, cm_int0, cm_int )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeMissingIntegrals
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
        DO l=1,Dofs
          CALL MPI_RECV( ValueTable2D(l) % Values(j:k), peNodes, &
            MPI_DOUBLE_PRECISION, pe, 104, ELMER_COMM_WORLD, status, ierr )
        END DO
      END DO
    ELSE
      j = myStart
      k = j+myNodes-1
      DO l=1,Dofs
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
    
    ALLOCATE(TmpValues(Dofs,olds),TmpExtent(olds))
    TmpExtent(1:olds) = PS % IntExtent(1:olds)
    TmpValues(:,1:olds) = PS % IntValues(:,1:olds) 
    DEALLOCATE(PS % IntValues, PS % IntExtent )

    news = olds + incr
    ALLOCATE( PS % IntValues(Dofs,news), PS % IntExtent(news) )
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
END SUBROUTINE ParallelProjectToPlane
!------------------------------------------------------------------------------
