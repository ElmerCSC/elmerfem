! *****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Apr 2001
! *
! *****************************************************************************/
  
!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Utilities used for interface conditions of various types including mortar.
!------------------------------------------------------------------------------

MODULE GeometryFitting

  USE Types
  USE Messages
  USE ElementDescription
  USE ElementUtils, ONLY : TangentDirections
  USE Interpolation, ONLY : CopyElementNodesFromMesh
  USE Lists
  USE ParallelUtils
  IMPLICIT NONE

CONTAINS


  !---------------------------------------------------------------------------
  ! Simply fitting of cylinder into a point cloud. This is done in two phases.
  ! 1) The axis of the cylinder is found by minimizing the \sum((n_i*t)^2)
  !    for each component of of t where n_i:s are the surface normals. 
  !    This is fully generic and assumes no positions. 
  ! 2) The radius and center point of the cylinder are found by fitting a circle
  !    in the chosen plane to three representative points. Currently the fitting
  !    can only be done in x-y plane. 
  !---------------------------------------------------------------------------
  SUBROUTINE CylinderFit(PMesh, PParams, BCind, dim, FitParams) 
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    TYPE(Valuelist_t), POINTER :: PParams
    INTEGER, OPTIONAL :: BCind
    INTEGER, OPTIONAL :: dim
    REAL(KIND=dp), OPTIONAL :: FitParams(:)
    
    INTEGER :: i,j,k,n,t,AxisI,iter,cdim,ierr
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: NiNj(9),A(3,3),F(3),M11,M12,M13,M14
    REAL(KIND=dp) :: d1,d2,MinDist,MaxDist,Dist,X0,Y0,Rad
    REAL(KIND=dp) :: Normal(3), AxisNormal(3), Tangent1(3), Tangent2(3), Coord(3), &
        CircleCoord(9)
#ifdef ELMER_BROKEN_MPI_IN_PLACE
    REAL(KIND=dp) :: buffer(9)
#endif
    INTEGER :: CircleInd(3) 
    LOGICAL :: BCMode, DoIt, GotNormal, GotCenter, GotRadius
    INTEGER :: Tag, t1, t2
    LOGICAL, ALLOCATABLE :: ActiveNode(:)
    REAL(KIND=dp), POINTER :: rArray(:,:)
    
    BCMode = PRESENT( BCind )

    ! Set the range for the possible active elements. 
    IF( BCMode ) THEN
      t1 = PMesh % NumberOfBulkElements + 1
      t2 = PMesh % NumberOfBulkElements + PMesh % NumberOfBoundaryElements
      Tag = CurrentModel % BCs(BCind) % Tag
      ALLOCATE( ActiveNode( PMesh % NumberOfNodes ) )
      ActiveNode = .FALSE.
    ELSE
      t1 = 1
      t2 = PMesh % NumberOfBulkElements      
    END IF
    
    ! If this is a line mesh there is really no need to figure out the 
    ! direction of the rotational axis. It can only be aligned with the z-axis.
    DO t=t1, t2
      Element => PMesh % Elements(t)
      IF( BCMode ) THEN
        IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE     
        IF ( Element % BoundaryInfo % Constraint /= Tag ) CYCLE
      END IF
      IF( Element % TYPE % ElementCode < 300 ) THEN
        cdim = 2
      ELSE
        cdim = 3
      END IF
      EXIT
    END DO
    
    IF( BcMode ) THEN
      cdim = ParallelReduction( cdim, 2 )
    END IF

    AxisNormal = 0.0_dp
    IF( cdim == 2 ) THEN
      GotNormal = .TRUE.
      AxisNormal(3) = 1.0_dp
    ELSE      
      rArray => ListGetConstRealArray( PParams,'Cylinder Normal',GotNormal)
      IF( GotNormal) AxisNormal(1:3) = rArray(1:3,1)
    END IF

    Coord = 0.0_dp
    rArray => ListGetConstRealArray( PParams,'Cylinder Center',GotCenter)
    IF( GotCenter) Coord(1:cdim) = rArray(1:cdim,1)
    
    Rad = ListGetConstReal( PParams,'Cylinder Radius',GotRadius)
 
    ! Do we have the fitting done already? 
    IF( GotNormal .AND. GotCenter .AND. GotRadius ) THEN
      IF( PRESENT(FitParams) ) THEN
        CALL Info('CylinderFit','Using cylinder parameters from list',Level=25)
        FitParams(1:cdim) = Coord(1:cdim)
        IF( cdim == 2 ) THEN
          FitParams(3) = Rad
        ELSE
          FitParams(4:6) = AxisNormal
          FitParams(7) = Rad
        END IF
      END IF
      RETURN
    END IF
                  
    n = PMesh % MaxElementNodes
    ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

       
    ! Compute the inner product of <N*N> for the elements
    NiNj = 0.0_dp
    DO t=t1, t2
      Element => PMesh % Elements(t)

      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes

      IF( BCMode ) THEN
        IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE     
        IF ( Element % BoundaryInfo % Constraint /= Tag ) CYCLE
        ActiveNode(Element % NodeIndexes(1:n)) = .TRUE.
      END IF
              
      ! If we know the Normal we only tag the boundary nodes
      IF(GotNormal) CYCLE

      Nodes % x(1:n) = PMesh % Nodes % x(NodeIndexes(1:n))
      Nodes % y(1:n) = PMesh % Nodes % y(NodeIndexes(1:n))
      Nodes % z(1:n) = PMesh % Nodes % z(NodeIndexes(1:n))           
      
      Normal = NormalVector( Element, Nodes, Check = .FALSE. ) 
      DO i=1,3
        DO j=1,3
          NiNj(3*(i-1)+j) = NiNj(3*(i-1)+j) + Normal(i) * Normal(j)
        END DO
      END DO
    END DO
      
    IF(GotNormal) GOTO 100 

    ! Only in BC mode we do currently parallel reduction.
    ! This could be altered too. 
    IF( BCMode ) THEN
#ifdef ELMER_BROKEN_MPI_IN_PLACE
      buffer = NiNj
      CALL MPI_ALLREDUCE(buffer, &
#else
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, &
#endif
          NiNj,9,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
    END IF
      
    ! The potential direction for the cylinder axis is the direction with 
    ! least hits for the normal.
    AxisI = 1 
    DO i=2,3
      IF( NiNj(3*(i-1)+i) < NiNj(3*(AxisI-1)+AxisI) ) AxisI = i 
    END DO

    CALL Info('CylinderFit','Axis coordinate set to be: '//I2S(AxisI))

    ! Keep the dominating direction fixed and iteratively solve the two other directions
    AxisNormal = 0.0_dp
    AxisNormal(AxisI) = 1.0_dp

    ! Basically we could solve from equation Ax=0 the tangent but only up to a constant.
    ! Thus we enforce the axis direction to one by manipulation the matrix equation 
    ! thereby can get a unique solution. 
    DO i=1,3
      DO j=1,3
        A(i,j) = NiNj(3*(i-1)+j)
      END DO
    END DO
    A(AxisI,1:3) = 0.0_dp
    A(AxisI,AxisI) = 1.0_dp
    CALL InvertMatrix( A, 3 )
    AxisNormal = A(1:3,AxisI)

    ! Normalize the axis normal length to one    
    AxisNormal = AxisNormal / SQRT( SUM( AxisNormal ** 2 ) )
    IF( 1.0_dp - MAXVAL( ABS( AxisNormal ) ) > 1.0d-5 ) THEN
      CALL Warn('CylinderFit','The cylinder axis is not aligned with any axis!')
    END IF

100 CALL TangentDirections( AxisNormal,Tangent1,Tangent2 )

    IF( InfoActive(25) .AND. ParEnv % MyPe == 0 ) THEN
      PRINT *,'Axis Normal:',AxisNormal
      PRINT *,'Axis Tangent 1:',Tangent1
      PRINT *,'Axis Tangent 2:',Tangent2
      i = PMesh % NumberOfNodes
      IF(BcMode) THEN
        PRINT *,'Active nodes: ',i,COUNT(ActiveNode)
      END IF
    END IF

    ! Finding three points with maximum distance in the tangent directions

    ! First, find the single extremum point in the first tangent direction
    ! Save the local coordinates in the N-T system of the cylinder
    MinDist = HUGE(MinDist)
    MaxDist = -HUGE(MaxDist)

    CIrcleInd = 0
    DO i=1, PMesh % NumberOfNodes
      IF( BCMode ) THEN
        IF( .NOT. ActiveNode(i) ) CYCLE
      END IF
      
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)

      d1 = SUM( Tangent1 * Coord )
      IF( d1 < MinDist ) THEN
        MinDist = d1
        CircleInd(1) = i
      END IF
      IF( d1 > MaxDist ) THEN
        MaxDist = d1
        CircleInd(2) = i
      END IF
    END DO

    CircleCoord = -HUGE(CircleCoord)
    DO j=1,2    
      i = CircleInd(j)
      
      IF( BCMode .AND. ParEnv % PEs > 1 ) THEN
        IF(j==1) THEN
          Dist = ParallelReduction( MinDist, 1 )
          IF(ABS(MinDist-Dist) > 1.0e-8) CYCLE
        ELSE IF(j==2) THEN
          Dist = ParallelReduction( MaxDist, 2)
          IF(ABS(MaxDist-Dist) > 1.0e-8) CYCLE
        END IF
      END IF
        
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)
      
      CircleCoord(3*(j-1)+1) = SUM( Tangent1 * Coord ) 
      CircleCoord(3*(j-1)+2) = SUM( Tangent2 * Coord ) 
      CircleCoord(3*(j-1)+3) = SUM( AxisNormal * Coord )
    END DO

    IF( BCMode .AND. ParEnv % PEs > 1 ) THEN
#ifdef ELMER_BROKEN_MPI_IN_PLACE
      buffer = CircleCoord
      CALL MPI_ALLREDUCE(buffer, &
#else
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, &
#endif
          CircleCoord,6,MPI_DOUBLE_PRECISION,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF

    IF( InfoActive(25) .AND. ParEnv % MyPe == 0 ) THEN
      PRINT *,'Circle Coord:',CircleCoord(1:6)
    END IF
    
    ! Find one more point such that their minimum distance to the previous point(s)
    ! is maximized. This takes some time but the further the nodes are apart the more 
    ! accurate it will be to fit the circle to the points. Also if there is just 
    ! a symmetric section of the cylinder it is important to find the points rigorously.
    j = 3
    ! The maximum minimum distance of any node from the previously defined nodes
    MaxDist = 0.0_dp
    DO i=1, PMesh % NumberOfNodes
      IF( BCMode ) THEN
        IF( .NOT. ActiveNode(i) ) CYCLE
      END IF
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)
      
      ! Minimum distance from the previously defined nodes
      MinDist = HUGE(MinDist)
      DO k=1,j-1
        d1 = SUM( Tangent1 * Coord )
        d2 = SUM( Tangent2 * Coord )
        Dist = ( d1 - CircleCoord(3*(k-1)+1) )**2 + ( d2 - CircleCoord(3*(k-1)+2) )**2
        MinDist = MIN( Dist, MinDist )
      END DO
      
      ! If the minimum distance to either previous selelected nodes
      ! is greater than in any other node, choose this
      IF( MaxDist < MinDist ) THEN
        MaxDist = MinDist 
        CircleInd(j) = i
      END IF
    END DO
    
    ! Ok, we have found the point now set the circle coordinates 
    DoIt = .TRUE.
    IF( BCMode .AND. ParEnv % PEs > 1 ) THEN
      Dist = ParallelReduction( MaxDist, 2 )
      DoIt = ( ABS(MaxDist-Dist) < 1.0e-8 )
    END IF

    IF( DoIt ) THEN
      i = CircleInd(j)
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)
      
      CircleCoord(3*(j-1)+1) = SUM( Tangent1 * Coord ) 
      CircleCoord(3*(j-1)+2) = SUM( Tangent2 * Coord ) 
      CircleCoord(3*(j-1)+3) = SUM( AxisNormal * Coord )
    END IF

    IF( BCMode .AND. ParEnv % PEs > 1 ) THEN
#ifdef ELMER_BROKEN_MPI_IN_PLACE
      buffer = CircleCoord
      CALL MPI_ALLREDUCE(buffer, &
#else
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, &
#endif
          CircleCoord,9,MPI_DOUBLE_PRECISION,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
      
    IF( InfoActive(25) .AND. ParEnv % MyPe == 0 ) THEN
      DO i=1,3
        PRINT *,'Circle Coord:',i,CircleInd(i),CircleCoord(3*i-2:3*i) 
      END DO
    END IF
      
    ! Given three nodes it is possible to analytically compute the center point and
    ! radius of the cylinder from a 4x4 determinant equation. The matrices values
    ! m1i are the determinants of the comatrices. 

    A(1:3,1) = CircleCoord(1::3)  ! x
    A(1:3,2) = CircleCoord(2::3)  ! y
    A(1:3,3) = 1.0_dp
    m11 = Det3x3( a )

    A(1:3,1) = CircleCoord(1::3)**2 + CircleCoord(2::3)**2  ! x^2+y^2
    A(1:3,2) = CircleCoord(2::3)  ! y
    A(1:3,3) = 1.0_dp
    m12 = Det3x3( a )
 
    A(1:3,1) = CircleCoord(1::3)**2 + CircleCoord(2::3)**2  ! x^2+y^2
    A(1:3,2) = CircleCoord(1::3)  ! x
    A(1:3,3) = 1.0_dp
    m13 = Det3x3( a )
 
    A(1:3,1) = CircleCoord(1::3)**2 + CircleCoord(2::3)**2 ! x^2+y^2
    A(1:3,2) = CircleCoord(1::3)  ! x
    A(1:3,3) = CircleCoord(2::3)  ! y
    m14 = Det3x3( a )

    IF(InfoActive(25) .AND. ParEnv % Mype == 0 ) THEN
      PRINT *,'CylinderFit determinants:',m11,m12,m13,m14
    END IF
      
    IF( ABS( m11 ) < EPSILON( m11 ) ) THEN
      CALL Fatal('CylinderFit','Points cannot be an a circle')
    END IF

    X0 =  0.5_dp * m12 / m11 
    Y0 = -0.5_dp * m13 / m11
    rad = SQRT( x0**2 + y0**2 + m14/m11 )

    Coord = x0 * Tangent1 + y0 * Tangent2

    IF( InfoActive(25) .AND. ParEnv % MyPe == 0) THEN
      PRINT *,'Cylinder center and radius:',Coord, rad
    END IF

    ALLOCATE( rArray(3,1) )
    rArray(1:3,1) = Coord 
    CALL ListAddConstRealArray( PParams,'Cylinder Center', 3, 1, rArray ) 
    IF(.NOT. GotNormal ) THEN
      rArray(1:3,1) = AxisNormal 
      CALL ListAddConstRealArray( PParams,'Cylinder Normal', 3, 1, rArray ) 
    END IF
    DEALLOCATE( rArray ) 
    CALL ListAddConstReal( PParams,'Cylinder Radius',rad )

    IF( PRESENT( FitParams ) ) THEN
      IF( cdim == 2 ) THEN
        FitParams(1:2) = Coord(1:2)
        FitParams(3) = rad
      ELSE
        FitParams(1:3) = Coord(1:3)
        FitParams(4:6) = AxisNormal(1:3)
        FitParams(7) = rad
      END IF

      IF( InfoActive(25) .AND. ParEnv % MyPe == 0) THEN
        PRINT *,'Cylinder FitParams: ',FitParams 
      END IF

    END IF
      
    DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )

  END SUBROUTINE CylinderFit


  ! Computes the center of a mesh or given set of bodies.
  !----------------------------------------------------------------------------  
  SUBROUTINE ComputeEntityCenter(Mesh, Center, TargetBodies, TargetBCs)
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Center(3)
    INTEGER, POINTER, OPTIONAL :: TargetBodies(:)
    INTEGER, POINTER, OPTIONAL :: TargetBCs(:)

    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    REAL(KIND=dp) :: DetJ,r(3),s
    INTEGER :: t,t1,tend,i,j,k,n,ierr
    LOGICAL :: stat
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: Volume,SerTmp(4),ParTmp(4)

    n = Mesh % MaxElementNodes
    ALLOCATE( Basis(n) )
    
    Volume = 0.0_dp
    Center = 0.0_dp

    IF(PRESENT(TargetBCs)) THEN
      t1 = Mesh % NumberOfBulkElements+1
      tend = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    ELSE
      t1 = 1
      tend = Mesh % NumberOfBulkElements
    END IF
          
    DO t=t1, tend
      Element => Mesh % Elements(t)
      IF( PRESENT( TargetBodies ) ) THEN
        IF( ALL( TargetBodies /= Element % BodyId ) ) CYCLE
      END IF
      IF( PRESENT( TargetBCs ) ) THEN
        IF(.NOT. ASSOCIATED(Element % BoundaryInfo) ) CYCLE
        i = Element % BoundaryInfo % Constraint
        IF( ALL( TargetBCs /= i ) ) CYCLE
      END IF
           
      n  = Element % Type % NumberOfNodes
      CALL CopyElementNodesFromMesh(Nodes,Mesh,n,Element % NodeIndexes)
      
      ! Numerical integration:
      !----------------------
      IP = GaussPoints(Element)

      DO k=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
            IP % W(k), detJ, Basis )
        
        r(1) = SUM(Nodes % x(1:n) * Basis(1:n))
        r(2) = SUM(Nodes % y(1:n) * Basis(1:n))
        r(3) = SUM(Nodes % z(1:n) * Basis(1:n))        
        s = IP % s(k) * detJ
        
        Volume = Volume + s
        Center = Center + s * r 
      END DO
    END DO

    IF( ParEnv % PEs > 1 ) THEN
      SerTmp(1:3) = Center
      SerTmp(4) = Volume
      CALL MPI_ALLREDUCE(SerTmp,ParTmp,4,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
      Center = ParTmp(1:3)
      Volume = ParTmp(4)
    END IF

    IF( Volume < EPSILON( Volume ) ) CALL Fatal('ComputeEntityCenter','Entity has no volume!')

    Center = Center / Volume
    
    WRITE( Message,'(A,ES12.4)') 'Body volume:',Volume
    CALL Info('ComputeEntityCenter',Message,Level=20)

    WRITE( Message,'(A,3ES12.4)') 'Body center:',Center
    CALL Info('ComputeEntityCenter',Message,Level=20)
    
  END SUBROUTINE ComputeEntityCenter


  ! Computes the normal of inertia of a mesh or given set of bodies.
  !----------------------------------------------------------------------------  
  SUBROUTINE ComputeEntityInertiaNormal(Mesh, Center, INormal, TargetBodies, TargetBCs)
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Center(3)
    REAL(KIND=dp) :: INormal(3)
    INTEGER, POINTER, OPTIONAL :: TargetBodies(:)
    INTEGER, POINTER, OPTIONAL :: TargetBCs(:)

    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    REAL(KIND=dp) :: DetJ,r(3),s
    INTEGER :: t,t1,tend,i,j,k,n,ierr
    LOGICAL :: stat
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp) :: Imoment(9), EigVec(3,3), EigVal(3), ParTmp(9), CP(3)
    REAL(KIND=dp) :: EigWrk(20)
    INTEGER :: EigInfo, Three
    TYPE(GaussIntegrationPoints_t) :: IP

    n = Mesh % MaxElementNodes
    ALLOCATE( Basis(n) )   
    Imoment = 0.0_dp

    IF(PRESENT(TargetBCs)) THEN
      t1 = Mesh % NumberOfBulkElements+1
      tend = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    ELSE
      t1 = 1
      tend = Mesh % NumberOfBulkElements
    END IF
    
    DO t=t1,tend
      Element => Mesh % Elements(t)
      IF( PRESENT( TargetBodies ) ) THEN
        IF( ALL( TargetBodies /= Element % BodyId ) ) CYCLE
      END IF

      n  = Element % Type % NumberOfNodes
      CALL CopyElementNodesFromMesh(Nodes,Mesh,n,Element % NodeIndexes)

      ! Numerical integration:
      !----------------------
      IP = GaussPoints(Element)
      DO k=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
            IP % W(k), detJ, Basis )
        
        r(1) = SUM(Nodes % x(1:n) * Basis(1:n))
        r(2) = SUM(Nodes % y(1:n) * Basis(1:n))
        r(3) = SUM(Nodes % z(1:n) * Basis(1:n))        
        s = IP % s(k) * detJ
        r = r - Center
        
        DO i=1,3
          Imoment(3*(i-1)+i) = Imoment(3*(i-1)+i) + s * SUM( r**2 )
          DO j=1,3
            Imoment(3*(i-1)+j) = Imoment(3*(i-1)+j) - s * r(i) * r(j)
          END DO
        END DO
      END DO
    END DO

    IF( ParEnv % PEs > 1 ) THEN
      CALL MPI_ALLREDUCE(Imoment,ParTmp,9,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
      Imoment = ParTmp
    END IF

    s = 1.0_dp    
    DO i=1,3
      DO j=1,3
        EigVec(i,j) = Imoment(3*(i-1)+j)
      END DO
      EigVec(i,i) = EigVec(i,i) - s 
    END DO

    EigInfo = 0
    Three = 3
    
    CALL DSYEV( 'V','U', Three, EigVec, Three, EigVal, EigWrk, SIZE(EigWrk), EigInfo )
    IF (EigInfo /= 0) THEN 
      CALL Fatal('ComputeEntityIntertiaNormal', 'DSYEV cannot generate eigen basis')
    END IF

    WRITE( Message,'(A,3ES12.4)') 'Mesh inertia eigenvalues:',EigVal
    CALL Info('ComputeEntityIntertiaNormal',Message,Level=30)
    INormal = EigVec(:,3)  ! axis of maximum inertia

    ! Check the sign of the normal using the right-hand-rule.
    ! This is not generic but a rule is still a rule
    CP = CrossProduct( Center, INormal )
    j = 1 
    DO i = 2, 3
      IF( ABS( CP(i) ) > ABS( CP(j) ) ) j = i
    END DO
    IF( CP(j) < 0 ) THEN
      CALL Info('ComputeEntityIntertiaNormal','Inverting sign of normal',Level=20)
      INormal = -INormal
    END IF    

  END SUBROUTINE ComputeEntityInertiaNormal

    
  
  !---------------------------------------------------------------------------
  SUBROUTINE TorusFit(PMesh, PParams, BCind, FitParams) 
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    TYPE(Valuelist_t), POINTER :: PParams
    INTEGER, OPTIONAL :: BCind
    REAL(KIND=dp), OPTIONAL :: FitParams(:)
    
    REAL(KIND=dp) :: Center(3), Normal(3), Rminor, Rmajor, rArray(3,1)
    LOGICAL :: Found
    INTEGER, POINTER :: EntityInds(:)
    REAL(KIND=dp), POINTER :: pArray(:,:) 
    
    ALLOCATE(EntityInds(1))
    EntityInds(1) = BCInd

    pArray => ListGetConstRealArray( PParams,'Torus Center',Found)
    IF(Found ) THEN
      Center(1:3) = pArray(1:3,1)
    ELSE      
      CALL ComputeEntityCenter(PMesh, Center, TargetBCs = EntityInds )
      rArray(1:3,1) = Center
      CALL ListAddConstRealArray( PParams,'Torus Center',3,1,rArray)
    END IF

    pArray => ListGetConstRealArray( PParams,'Torus Normal',Found )
    IF(Found ) THEN
      Normal(1:3) = pArray(1:3,1)
    ELSE      
      CALL ComputeEntityInertiaNormal(PMesh, Center, Normal, TargetBCs = EntityInds )
      rArray(1:3,1) = Normal
      CALL ListAddConstRealArray( PParams,'Torus Normal',3,1,rArray)
    END IF
       
    Rmajor = ListGetConstReal( PParams,'Torus Radius',UnfoundFatal=.TRUE.)
    Rminor = ListGetConstReal( PParams,'Torus Minor Radius',UnfoundFatal=.TRUE.)    

    FitParams(1:3) = Center
    FitParams(4:6) = Normal
    FitParams(7) = Rmajor
    FitParams(8) = Rminor

    DEALLOCATE(EntityInds)
    
  END SUBROUTINE TorusFit

  
  ! Code for fitting a sphere. Not yet used.
  !-------------------------------------------------------------------------
  SUBROUTINE SphereFit(Mesh, Params, BCind, FitParams ) 
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params
    INTEGER, OPTIONAL :: BCind
    REAL(KIND=dp), OPTIONAL :: FitParams(:)

    INTEGER :: i,j,t,t1,t2,NoNodes,Tag
    LOGICAL :: BCMode
    LOGICAL, ALLOCATABLE :: ActiveNode(:)
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)    
    REAL(KIND=dp) :: xc,yc,zc,Rad

    IF( PRESENT( FitParams ) ) THEN
      IF( ListCheckPresent( Params,'Sphere Radius') ) THEN
        CALL Info('SphereFit','Using predefined values for sphere parameters',Level=20)
        FitParams(1) = ListGetConstReal( Params,'Sphere Center X')
        FitParams(2) = ListGetConstReal( Params,'Sphere Center Y')
        FitParams(3) = ListGetConstReal( Params,'Sphere Center Z')
        FitParams(4) = ListGetConstReal( Params,'Sphere Radius')
        RETURN
      END IF
    END IF
          
    CALL Info('SphereFit','Trying to fit a sphere to element patch',Level=6)

    ! Set the range for the possible active elements. 
    IF( PRESENT( BCind ) ) THEN
      BCMode = .TRUE.
      t1 = Mesh % NumberOfBulkElements + 1
      t2 = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Tag = CurrentModel % BCs(BCind) % Tag
      ALLOCATE( ActiveNode( Mesh % NumberOfNodes ) )
      ActiveNode = .FALSE.
    ELSE
      BCMode = .FALSE.
      t1 = 1
      t2 = Mesh % NumberOfBulkElements
    END IF

    ! Mark the nodes that belong to the active elements.
    ! 1) Either we only have bulk elements in which case we use all of the nodes or
    ! 2) We are given a boundary index and only use the nodes related to it. 
    DO t=t1,t2
      Element => Mesh % Elements(t)
      IF( BCMode ) THEN
        IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE     
        IF ( Element % BoundaryInfo % Constraint /= Tag ) CYCLE
        ActiveNode(Element % NodeIndexes) = .TRUE.
      END IF
    END DO

    ! If all nodes are active just use pointers to the nodes.
    ! Otherwise create list of the nodes. 
    IF( BCMode ) THEN
      NoNodes = COUNT( ActiveNode )
      ALLOCATE( x(NoNodes), y(NoNodes), z(NoNodes) )
      j = 0
      DO i=1,Mesh % NumberOfNodes
        IF(.NOT. ActiveNode(i) ) CYCLE
        j = j + 1
        x(j) = Mesh % Nodes % x(i)
        y(j) = Mesh % Nodes % y(i)
        z(j) = Mesh % Nodes % z(i)
      END DO
    ELSE
      NoNodes = Mesh % NumberOfNodes
      x => Mesh % Nodes % x
      y => Mesh % Nodes % y
      z => Mesh % Nodes % z
    END IF

    ! Call the function to set the sphere parameters for the nodes.
    CALL SphereFitfun(NoNodes,x,y,z,xc,yc,zc,Rad)

    IF( BCMode ) THEN
      DEALLOCATE(ActiveNode,x,y,z)
    END IF

    ! Add the sphere parameters to the list so that they can be used later
    ! directly without having to fit the parameters again.  
    CALL ListAddConstReal( Params,'Sphere Center X',xc )
    CALL ListAddConstReal( Params,'Sphere Center Y',yc )
    CALL ListAddConstReal( Params,'Sphere Center Z',zc )
    CALL ListAddConstReal( Params,'Sphere Radius',Rad )
    
    IF( PRESENT( FitParams ) ) THEN
      FitParams(1) = xc
      FitParams(2) = yc
      FitParams(3) = zc
      FitParams(4) = Rad
    END IF
      
  CONTAINS
    

    ! Sumith YD: Fast Geometric Fit Algorithm for Sphere Using Exact Solution
    !------------------------------------------------------------------------
    SUBROUTINE SphereFitfun(n,x,y,z,xc,yc,zc,R)
      INTEGER :: n
      REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
      REAL(KIND=dp) :: xc,yc,zc,R
      
      REAL(KIND=dp) :: Sx,Sy,Sz,Sxx,Syy,Szz,Sxy,Sxz,Syz,&
          Sxxx,Syyy,Szzz,Syzz,Sxyy,Sxzz,Sxxy,Sxxz,Syyz,&
          A1,a,b,c,d,e,f,g,h,j,k,l,m,delta
      
      Sx = SUM(x); Sy = SUM(y); Sz = SUM(z);
      Sxx = SUM(x*x); Syy = SUM(y*y);
      Szz = SUM(z*z); Sxy = SUM(x*y);
      Sxz = SUM(x*z); Syz = SUM(y*z);
      Sxxx = SUM(x*x*x); Syyy = SUM(y*y*y);
      Szzz = SUM(z*z*z); Sxyy = SUM(x*y*y);
      Sxzz = SUM(x*z*z); Sxxy = SUM(x*x*y);
      Sxxz = SUM(x*x*z); Syyz = SUM(y*y*z);
      Syzz = SUM(y*z*z);

      ! We must do parallel reduction here if the surface is split among
      ! several MPI processes. 
      IF( BCMode .AND. ParEnv % PEs > 1 ) THEN
        Sx = ParallelReduction(Sx); Sy = ParallelReduction(Sy); Sz = ParallelReduction(Sz);
        Sxx = ParallelReduction(Sxx); Syy = ParallelReduction(Syy);
        Szz = ParallelReduction(Szz); Sxy = ParallelReduction(Sxy);
        Sxz = ParallelReduction(Sxz); Syz = ParallelReduction(Syz);
        Sxxx = ParallelReduction(Sxxx); Syyy = ParallelReduction(Syyy);
        Szzz = ParallelReduction(Szzz); Sxyy = ParallelReduction(Sxyy);
        Sxzz = ParallelReduction(Sxzz); Sxxy = ParallelReduction(Sxxy);
        Sxxz = ParallelReduction(Sxxz); Syyz = ParallelReduction(Syyz);
        Syzz = ParallelReduction(Syzz);       
      END IF
           
      A1 = Sxx +Syy +Szz;
      a = 2*Sx*Sx-2*N*Sxx;
      b = 2*Sx*Sy-2*N*Sxy;
      c = 2*Sx*Sz-2*N*Sxz;
      d = -N*(Sxxx +Sxyy +Sxzz)+A1*Sx;
      e = 2*Sx*Sy-2*N*Sxy;
      f = 2*Sy*Sy-2*N*Syy;
      g = 2*Sy*Sz-2*N*Syz;
      h = -N*(Sxxy +Syyy +Syzz)+A1*Sy;
      j = 2*Sx*Sz-2*N*Sxz;
      k = 2*Sy*Sz-2*N*Syz;
      l = 2*Sz*Sz-2*N*Szz;
      m = -N*(Sxxz +Syyz + Szzz)+A1*Sz;
      delta = a*(f*l - g*k)-e*(b*l-c*k) + j*(b*g-c*f);

      xc = (d*(f*l-g*k) -h*(b*l-c*k) +m*(b*g-c*f))/delta;
      yc = (a*(h*l-m*g) -e*(d*l-m*c) +j*(d*g-h*c))/delta;
      zc = (a*(f*m-h*k) -e*(b*m-d*k) +j*(b*h-d*f))/delta;
      R = SQRT(xc*xc+yc*yc+zc*zc+(A1-2*(xc*Sx+yc*Sy+zc*Sz))/N);

    END SUBROUTINE SphereFitfun

  END SUBROUTINE SphereFit

 

END MODULE GeometryFitting
  
