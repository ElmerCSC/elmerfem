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
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 3.3.2008
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Subroutine for solving steady-state free surface problems for structured meshes.
!> Intended applications include drawing and extrusion processes.
!>  This solver assumes that the mesh is structural so that it could have 
!>  been obtained by extrusion in the direction of interest. For the given 
!>  direction the corresponding top and bottom node is computed for every node
!>  and this information is used to perform linear mapping in between.  
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE StructuredFlowLine( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------

  USE Types
  USE Lists
  USE DefUtils
  USE MeshUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: PSolver
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, MaskName
  INTEGER :: i,j,k,l,n,dim,DOFs,itop,ibot,ii,jj,Rounds,BotMode,TopMode,nsize, &
      ActiveDirection,elem,FlowDofs,DispDofs,Hits,Shots, VisitedTimes = 0, &
	  AveMode, AveOrder
  INTEGER, POINTER :: TopPerm(:),BotPerm(:),TopPointer(:),MaskPerm(:),&
      BotPointer(:),UpPointer(:),DownPointer(:),NodeIndexes(:),FlowPerm(:),&
      DispPerm(:), HitPerm(:)
  LOGICAL :: GotIt, Found, Initialized = .FALSE., MaskExist, &
      DisplacementMode, TrueLines, OutOfMesh, GotHitStyle
  REAL(KIND=dp) :: UnitVector(3),x0loc,x0bot,x0top,xloc,wtop,BotVal,TopVal,&
      ElemVector(3),DotPro,Eps,Length,sumds
  REAL(KIND=dp) :: r1(3),r2(3),v1(3),v2(3),dr(3),vave(3),drproj(3),VelCor(3),&
      CoordMax,CoordMin,sumdr(3),Relax,Alpha,Norm
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at0,at1,at2
#else
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
#endif
  REAL(KIND=dp), POINTER :: Coord(:),BotField(:),TopField(:),Field(:),Surface(:),Flow(:), &
      HardDisp(:), HitField(:)
  REAL(KIND=dp) :: VeloAtPoint(3), GradVeloAtPoint(3,3),GlobalCoords(3)
  TYPE(Variable_t), POINTER :: Var
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t), SAVE :: Nodes
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(ValueList_t),POINTER :: BC
  
  
  SAVE VisitedTimes,Initialized,UnitVector,Coord,TopPointer,BotPointer,&
      DownPointer,UpPointer,TopMode,BotMode,TopField,BotField,TopPerm,BotPerm,Field,&
      Surface,Dim,HardDisp,DispPerm,MaskPerm,MaskExist,ActiveDirection
  

  CALL Info( 'StructuredFlowLine','---------------------------------------',Level=4 )
  CALL Info( 'StructuredFlowLine','Setting edges to follow the streamlines ',Level=4 )
  CALL Info( 'StructuredFlowLine','---------------------------------------',Level=4 )

!------------------------------------------------------------------------------
!   Initialize the pointers to top and bottom nodes 
!------------------------------------------------------------------------------

  SolverParams => GetSolverParams()
  Mesh => Solver % Mesh 
   
  ! the mesh is chanching, thus always rebuild the quadrant tree
  CALL FreeQuadrantTree( Mesh % RootQuadrant )
  Solver % Mesh % RootQuadrant => Null()

  i = GetInteger( SolverParams,'True Flow Line Iterations')
  TrueLines = ( i > VisitedTimes )

  ! one is a good default for the averaging order
  AveOrder = GetInteger( SolverParams,'Averaging Order',GotIt)
  IF( TrueLines ) AveOrder = MAX( 1, AveOrder )
  
  ! zero is a good default
  AveMode = GetInteger( SolverParams,'Averaging Method',GotIt)
  DIM = CoordinateSystemDimension()
  

  IF( .NOT. Initialized ) THEN
    
    at0 = CPUTime()

    NULLIFY(Field)
    
    ! Choose active direction coordinate and set corresponding unit vector
    !---------------------------------------------------------------------
    ActiveDirection = GetInteger(SolverParams,'Active Coordinate')
    IF(ActiveDirection == 1) THEN
      Coord => Solver % Mesh % Nodes % x
    ELSE IF(ActiveDirection == 2) THEN
      Coord => Solver % Mesh % Nodes % y
    ELSE IF(ActiveDirection == 3) THEN 
      Coord => Solver % Mesh % Nodes % z
    ELSE 
      CALL Fatal('StructuredFlowLine','Invalid value for Active Coordinate')
    END IF
    UnitVector = 0.0_dp
    UnitVector(ActiveDirection) = 1.0_dp
    
    ! Make the mask used in matrix structures and permulation
    !-------------------------------------------------------  

    MaskName = 'Flow Line'
    n = Mesh % NumberOfNodes
    ALLOCATE( MaskPerm(n) )
    n = 0
    MaskPerm = 0
    Psolver => Solver
    CALL MakePermUsingMask( Model,PSolver,Mesh,TRIM(MaskName), &
        .FALSE., MaskPerm, n )
    IF( n == 0 ) THEN
      MaskExist = .FALSE.
      DEALLOCATE( MaskPerm )
    ELSE
      MaskExist = .TRUE.
      WRITE(Message,'(A,I6)') 'Number of nodes on '//TRIM(MaskName),n
      CALL Info('StructuredFlowLine',Message) 
    END IF

    ! Allocate pointers to top and bottom, and temporary pointers up and down
    !------------------------------------------------------------------------
    nsize = SIZE( Coord )
    
    ALLOCATE(TopPointer(nsize),UpPointer(nsize))
    DO i=1,nsize
      TopPointer(i) = i
      UpPointer(i) = i
    END DO
    ALLOCATE(BotPointer(nsize),DownPointer(nsize))
    DO i=1,nsize
      BotPointer(i) = i
      DownPointer(i) = i
    END DO
        
    ! Determine the pointer up and down using dot product as criterion
    !-----------------------------------------------------------------
    Eps = GetConstReal(SolverParams,'Dot Product Tolerance',GotIt)
    IF(.NOT. GotIt) Eps = 1.0e-4
    
    IF(.FALSE.) PRINT *,'determine up and down pointers'
    
    DO elem = 1,Mesh % NumberOfBulkElements      
      Element => Mesh % Elements(elem)
      NodeIndexes => Element % NodeIndexes

      IF( MaskExist ) THEN
        IF( ALL(MaskPerm(NodeIndexes) == 0) ) CYCLE
      END IF

      Model % CurrentElement => Element
      CALL GetElementNodes( Nodes )
      n  = GetElementNOFNodes()
            
      DO i=1,n
        ii = NodeIndexes(i)
        IF( MaskExist ) THEN
          IF( MaskPerm(ii) == 0 ) CYCLE          
        END IF

        DO j=i+1,n
          jj = NodeIndexes(j)
          IF( MaskExist ) THEN
            IF( MaskPerm(jj) == 0 ) CYCLE          
          END IF
          
          ElemVector(1) = Nodes % x(j) - Nodes % x(i)
          ElemVector(2) = Nodes % y(j) - Nodes % y(i)
          ElemVector(3) = Nodes % z(j) - Nodes % z(i)
          Length = SQRT(SUM(ElemVector*ElemVector))
          DotPro = SUM(ElemVector * UnitVector) / Length
          
          IF(DotPro > 1.0_dp - Eps) THEN
            UpPointer(ii) = jj
            DownPointer(jj) = ii
          ELSE IF(DotPro < Eps - 1.0_dp) THEN
            DownPointer(ii) = jj
            UpPointer(jj) = ii
          END IF
        END DO
      END DO
    END DO
    
    ! Pointer to top and bottom are found recursively using up and down
    !------------------------------------------------------------------
    IF(.FALSE.) PRINT *,'determine top and bottom pointers'
    
    Found = .TRUE.
    Rounds = 0
    DO WHILE(Found) 
      Found = .FALSE.
      DO i=1,nsize
        IF( MaskExist ) THEN
          IF( MaskPerm(i) == 0) CYCLE
        END IF
        j = UpPointer(i)
        IF( TopPointer(i) /= TopPointer( j ) ) THEN
          Found = .TRUE.
          TopPointer(i) = TopPointer( j )
        END IF
        j = DownPointer(i)
        IF( BotPointer(i) /= BotPointer( j ) ) THEN
          Found = .TRUE.
          BotPointer(i) = BotPointer( j )
        END IF
      END DO
      IF( Found ) Rounds = Rounds + 1
    END DO
    
    at1 = CPUTime()
    
    WRITE(Message,* ) 'Top and bottom pointer init rounds: ',Rounds
    CALL Info('StructuredFlowLine',Message)
    WRITE(Message,'(A,F8.3)') 'Top and bottom pointer init time: ',at1-at0
    CALL Info('StructuredFlowLine',Message)
    
    Initialized = .TRUE.
  END IF
  at0 = CPUTime()

  ! End of initialization
  !-------------------------------------------------------


  ! Get the velocity variables
  !-------------------------------------------------------
  VarName = GetString(SolverParams,'Velocity Variable Name',GotIt )
  IF(.NOT. GotIt) VarName = 'Flow Solution'
  Var => VariableGet( Solver % Mesh % Variables,  VarName )
  IF(ASSOCIATED(Var)) THEN
    Flow => Var % Values
    FlowPerm => Var % Perm
    FlowDofs = Var % DOFs
  ELSE
    CALL Fatal('StructuredFlowLine','Velocity variable is missing: '//TRIM(VarName))
  END IF
  

  ! Get result field for suggested displacement
  !-------------------------------------------------------
  VarName = GetString(SolverParams,'Flow Line Displacement Name',GotIt )
  IF(.NOT. GotIt) VarName = 'FlowLineDisp'

  NULLIFY(Var)
  Var => VariableGet( Solver % Mesh % Variables, VarName, ThisOnly = .TRUE. )
  IF(ASSOCIATED(Var)) THEN
    HardDisp => Var % Values
    DispPerm => Var % Perm
    DispDofs = Var % DOFs
  ELSE
    CALL Fatal('StructuredFlowLine','Flow line displacement variable does not exist!')
  END IF    
  
  IF( DispDofs /= DIM ) THEN 
    CALL Fatal('StructuredFlowLine','Displacement variable of wrong size')
  END IF
  
  ! This is just for testing purposes to mark hits 'inside' and 'outside'
  !-----------------------------------------------------------------------
  VarName = GetString(SolverParams,'Hit Style Name',GotIt )
  IF(.NOT. GotIt) VarName = 'HitStyle'
  Var => VariableGet( Solver % Mesh % Variables,VarName )
  GotHitStyle = ASSOCIATED(Var)
  IF( GotHitStyle ) THEN
    HitField => Var % Values
    HitPerm => Var % Perm
  END IF
  
  
  DisplacementMode = GetLogical(SolverParams,'Displacement Mode',Found)
  
  ! Map to the streamline
  !-------------------------------------------------------------------
  IF(.FALSE.) PRINT *,'perform linear mapping'
  Alpha = 0.5_dp
  CoordMin = ListGetConstReal( Solver % Values,'Min Coordinate',GotIt)
  IF(.NOT. GotIt) CoordMin = -HUGE(CoordMin)
  CoordMax = ListGetConstReal( Solver % Values,'Max Coordinate',GotIt)
  IF(.NOT. GotIt) CoordMax = HUGE(CoordMax)
  
  nsize = SIZE( Coord )
  
  Shots = 0
  Hits = 0
  
  v1 = 0; v2 = 0;
  
  DO i=1,nsize
    IF( DispPerm(i) == 0 ) CYCLE
    IF( MaskExist ) THEN
      IF( MaskPerm(i) == 0) CYCLE
    END IF
    itop = TopPointer(i)
    
    IF( i /= itop ) CYCLE	

    ! First node of streamline is initialized	
    !-------------------------------------------------------------------
    k = i
    v1(1) = Flow( FlowDofs*(FlowPerm(k)-1) + 1 ) 
    v1(2) = Flow( FlowDofs*(FlowPerm(k)-1) + 2 ) 
    IF( DIM == 3 ) THEN
      v1(3) = Flow( FlowDofs*(FlowPerm(k)-1) + 3 ) 
    END IF
    
    r1(1) = Solver % Mesh % Nodes % x(k)
    r1(2) = Solver % Mesh % Nodes % y(k)
    r1(3) = Solver % Mesh % Nodes % z(k)
    
    sumdr = 0.0_dp
    
    DO j=1,nsize
      k = DownPointer(k)
      IF( DispPerm(k) == 0 ) CYCLE      

      r2(1) = Solver % Mesh % Nodes % x(k)
      r2(2) = Solver % Mesh % Nodes % y(k)
      r2(3) = Solver % Mesh % Nodes % z(k)
      dr = r2 - r1
      r1 = r2

      v2(1) = Flow( FlowDofs*(FlowPerm(k)-1) + 1 ) 
      v2(2) = Flow( FlowDofs*(FlowPerm(k)-1) + 2 ) 
      IF( DIM == 3 ) THEN
        v2(3) = Flow( FlowDofs*(FlowPerm(k)-1) + 3 ) 
      END IF

      ! Check the passive regions
      ! For no-slip condition the surface is always fixed
      !-------------------------------------------------------------------
      IF( r2( ActiveDirection) > CoordMax .OR.  &
          r2( ActiveDirection) < CoordMin .OR.  &
          ABS( v2( ActiveDirection) ) < TINY( CoordMin ) ) THEN
        v1 = v2 
        CYCLE
      END IF
	  	  

      ! Simple average velocity from the mean of surface velocities 		  
      !-------------------------------------------------------------------
      vave = Alpha*v2 + (1-Alpha)*v1	  
      

      ! A more complex velocity estimation is preferable in the start
      !-------------------------------------------------------------------
      IF( TrueLines ) THEN 
        drproj = 0.0_dp
        sumds = SQRT( SUM ( sumdr ** 2) )
        
        DO l=1,AveOrder
          IF( sumds > 1.0d-20 ) THEN

            ! Averaging of velocity may be done in two ways:
            ! 1) taking velocity at average of r1 and r2
            ! 2) taking average of velocities at r1 and r2
            !-------------------------------------------------
            
            IF( AveMode == 1 ) THEN
              GlobalCoords = 0.5*(r1 + r2) + sumdr + 0.5*drproj
            ELSE
              GlobalCoords = r2 + sumdr + drproj 
            END IF
            CALL ComputeVeloAtPoint(GlobalCoords, k, OutOfMesh, VeloAtPoint ) 
            Shots = Shots + 1
            
            IF( OutOfMesh ) THEN
              IF( AveMode == 1 ) THEN
                GlobalCoords = 0.5*(r1 + r2) 
              ELSE
                GlobalCoords = r2 + sumdr + drproj 
              END IF
              
              CALL ComputeVeloAtPoint(GlobalCoords, k, OutOfMesh, VeloAtPoint, &
                  GradVeloAtPoint ) 
              IF( .NOT. OutOfMesh ) THEN
                IF( GotHitStyle ) HitField( HitPerm(k) ) = -1.0_dp * l
                VelCor = 0
                DO n=1,3
                  IF( n == ActiveDirection ) CYCLE
                  VelCor = VelCor + GradVeloAtPoint(1:3,n) * dr(n)
                END DO
                IF( AveMode == 1 ) THEN
                  vave = VeloAtPoint  + VelCor 
                ELSE
                  v2 = VeloAtPoint + VelCor
                  vave = Alpha * v2 + (1-Alpha) * v1
                END IF
              END IF
            ELSE  
              IF( GotHitStyle ) HitField( HitPerm(k) ) = 1.0_dp * l
              
              Hits = Hits + 1
              IF( AveMode == 1 ) THEN
                vave = VeloAtPoint  
              ELSE
                v2 = VeloAtPoint 
                vave = Alpha * v2 + (1-Alpha) * v1
              END IF
            END IF
          END IF
          drproj = vave * ABS( dr(ActiveDirection) ) / ABS( vave(ActiveDirection) ) - dr
          sumds = SQRT( SUM (( sumdr + drproj ) ** 2) )
        END DO
      ELSE
        drproj = vave * ABS( dr(ActiveDirection) ) / ABS( vave(ActiveDirection) ) - dr
      END IF
      
      
      sumdr = sumdr + drproj 	  
      sumdr( ActiveDirection ) = 0.0_dp
      
      DO l=1,dim	  
        HardDisp( DispDofs*(DispPerm(k)-1)+l) = sumdr(l)
      END DO
      
      v1 = v2
      
      ! This is the bottom and we're done?
      IF( k == BotPointer(k) ) EXIT
      
    END DO
  END DO

  Norm = SQRT( SUM(HardDisp**2) / MAXVAL( DispPerm  ) )
  WRITE( Message,* ) 'Result Norm: ',Norm
  CALL Info('StructuredFlowLine', Message ) 


  Relax = GetCReal( Solver % Values,'Nonlinear System Relaxation Factor',GotIt)
  IF( GotIt .AND. ABS(Relax - 1.0_dp) > TINY(Relax) ) THEN
    HardDisp = Relax * HardDisp 
  END IF
  
  IF( DisplacementMode ) THEN
    CALL Info('StructuredFlowLine','Moving the nodes')
    DO k=1,nsize
      IF( MaskExist ) THEN
        IF( MaskPerm(i) == 0) CYCLE
      END IF

      Mesh % Nodes % x(k) = Mesh % Nodes % x(k) + &
	      HardDisp( DispDofs*(DispPerm(k)-1)+1) 
      Mesh % Nodes % y(k) = Mesh % Nodes % y(k) + &
	      HardDisp( DispDofs*(DispPerm(k)-1)+2) 
      IF(dim == 3) Mesh % Nodes % z(k) = Mesh % Nodes % z(k) + &
	      HardDisp( DispDofs*(DispPerm(k)-1)+3) 
    END DO
  END IF	
  
  at1 = CPUTime()
  
  IF( TrueLines ) THEN
    WRITE(Message,*) 'Velocity estimations, total: ',Shots
    CALL Info('StructuredFlowLine',Message)
    WRITE(Message,*) 'Velocity estimations, inside: ',Hits
    CALL Info('StructuredFlowLine',Message)
  END IF
  
  WRITE(Message,'(A,F8.3)' ) 'Streamline mapping time: ',at1-at0
  CALL Info('StructuredFlowLine',Message)
  
  VisitedTimes = VisitedTimes + 1
  
   CONTAINS
  
  !------------------------------------------------------------------------
  !> Compute field values at the given points in the FE mesh. 
  !-------------------------------------------------------------------------
  SUBROUTINE ComputeVeloAtPoint(GlobalCoords, node, OutOfScope, VeloAtPoint, &
       GradVeloAtPoint, ElementSize ) 

    INTEGER :: node
    LOGICAL :: OutOfScope
    REAL(KIND=dp) :: GlobalCoords(3), GradientAtPoint(3),VeloAtPoint(3)
    REAL(KIND=dp), OPTIONAL :: GradVeloAtPoint(3,3), ElementSize


    LOGICAL :: Stat, Visited = .FALSE., Hit, DummySearch
    INTEGER :: i,j,k,n
    REAL(KIND=dp) :: u,v,w,LocalCoords(3),SqrtElementMetric
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
    TYPE(Nodes_t) :: ElementNodes
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: CurrentElement

    TYPE ElementP
      TYPE(Element_t), POINTER :: Element => Null()
    END TYPE ElementP
    
    TYPE(ElementP), ALLOCATABLE, SAVE :: Elements(:)

    TYPE(Quadrant_t), POINTER, SAVE :: RootQuadrant =>Null(), LeafQuadrant
    REAL(kind=dp) :: BoundingBox(6), eps2, eps1 = 1e-3


    SAVE :: CurrentElement, ElementNodes, Basis, dBasisdx

    IF( .NOT. Visited ) THEN
       n = Solver % Mesh % MaxElementNodes
       ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n))
       ALLOCATE( Basis(n), dBasisdx(n, 3) )
	   ALLOCATE( Elements(Mesh % NumberOfNodes) )
	   NULLIFY( CurrentElement )
       Visited = .TRUE.
    END IF
      
    GradientAtPoint = 0.0_dp
    VeloAtPoint = 0.0_dp
    Hit = .FALSE.
    DummySearch = .FALSE.


    ! Check that the previous hit is not hit even now
    !-------------------------------------------------
	
    Currentelement => elements(node) % element
    IF ( ASSOCIATED(CurrentElement) ) THEN
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      
      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
      
      IF ( PointInElement( CurrentElement, ElementNodes, &
          GlobalCoords, LocalCoords ) ) THEN
        Hit = .TRUE.
      END IF
    END IF


    IF(.NOT. Hit) THEN

      IF( dim == 2 .OR. DummySearch ) THEN
        
        !----------------------------------------------------------
        ! Go through all bulk elements in a dummy search.
        ! This algorithm is mainly here for debugging purposes.
        !----------------------------------------------------------
        DO k=1,Solver % Mesh % NumberOfBulkElements
          CurrentElement => Mesh % Elements(k)
          n = CurrentElement % TYPE % NumberOfNodes
          NodeIndexes => CurrentElement % NodeIndexes
          
          ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
          ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
          
          IF ( PointInElement( CurrentElement, ElementNodes, &
              GlobalCoords, LocalCoords ) ) THEN
            Hit = .TRUE.
            EXIT
          END IF
        END DO
        
      ELSE

        !-----------------------------------------------
        ! Find the right element using an octree search
        ! This is the preferred algorithms of the two.
        !-----------------------------------------------
        IF ( .NOT.ASSOCIATED( Mesh % RootQuadrant ) ) THEN
          BoundingBox(1) = MINVAL( Mesh % Nodes % x )
          BoundingBox(2) = MINVAL( Mesh % Nodes % y )
          BoundingBox(3) = MINVAL( Mesh % Nodes % z )
          BoundingBox(4) = MAXVAL( Mesh % Nodes % x )
          BoundingBox(5) = MAXVAL( Mesh % Nodes % y )
          BoundingBox(6) = MAXVAL( Mesh % Nodes % z )
          
          eps2 = eps1 * MAXVAL( BoundingBox(4:6) - BoundingBox(1:3) )
          BoundingBox(1:3) = BoundingBox(1:3) - eps2
          BoundingBox(4:6) = BoundingBox(4:6) + eps2
          
          CALL BuildQuadrantTree( Mesh,BoundingBox,Mesh % RootQuadrant)
        END IF
        
        RootQuadrant => Mesh % RootQuadrant
        IF ( ASSOCIATED(RootQuadrant) ) THEN
          NULLIFY(CurrentElement)
          CALL FindLeafElements(GlobalCoords, dim, RootQuadrant, LeafQuadrant)
          IF ( ASSOCIATED(LeafQuadrant) ) THEN
            DO k=1, LeafQuadrant % NElemsInQuadrant
              CurrentElement => Mesh % Elements(LeafQuadrant % Elements(k))
              
              n = CurrentElement % TYPE % NumberOfNodes
              NodeIndexes => CurrentElement % NodeIndexes
              
              ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
              ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
              ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
              
              IF ( PointInElement( CurrentElement, ElementNodes, &
                  GlobalCoords, LocalCoords ) ) THEN
                Hit = .TRUE.
                EXIT
              END IF
            END DO
          END IF
        END IF
        
      END IF
    END IF

!------------------------------------------------------------------------------


    IF ( Hit ) THEN
      
      Elements(node) % Element => CurrentElement
      
      u = LocalCoords(1)
      v = LocalCoords(2)
      w = LocalCoords(3)     
	 
      stat = ElementInfo( CurrentElement, ElementNodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx )
      
      ! Compute the representative size of the elements
      IF( PRESENT(ElementSize) ) ElementSize = SqrtElementMetric ** ( 1.0_dp / DIM )
      
      IF( ALL ( FlowPerm( NodeIndexes) > 0 )) THEN
        DO i=1,DIM
          VeloAtPoint(i) = SUM( Basis(1:n) * &
              Flow( FlowDofs*( FlowPerm( NodeIndexes )-1)+i ) )
          
          IF( PRESENT(GradVeloAtPoint) ) THEN
            DO j=1,DIM
              GradVeloAtPoint(i,j) = SUM( dBasisdx(1:n,j) * &
                  Flow( FlowDofs*( FlowPerm( NodeIndexes )-1)+i ) )                  
            END DO
          END IF
        END DO
      END IF
      
      ! for debugging purposes
      IF(.FALSE.) THEN
        PRINT *,'x:',SUM( Basis(1:n) * Mesh % Nodes % x( NodeIndexes ) )
        PRINT *,'y:',SUM( Basis(1:n) * Mesh % Nodes % y( NodeIndexes ) )
        PRINT *,'z:',SUM( Basis(1:n) * Mesh % Nodes % z( NodeIndexes ) )
        PRINT *,'VeloAtPoint:i',VeloAtPoint
        CALL flush(6)
      END IF
    END IF
    
    OutOfScope = .NOT. Hit
    
  END SUBROUTINE ComputeVeloAtPoint
 
!------------------------------------------------------------------------------
END SUBROUTINE StructuredFlowLine
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization of the primary solver: StrcuturedFlowLine. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE StructuredFlowLine_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE Types
  USE DefUtils
  USE Lists

  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
  INTEGER :: dim
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  TYPE(ValueList_t),POINTER :: SolverParams
  LOGICAL :: GotIt


  DIM = CoordinateSystemDimension()
  SolverParams => GetSolverParams()
  
  CALL ListAddLogical( SolverParams,'No Matrix',.TRUE.)
  
  VarName = GetString(SolverParams,'Flow Line Displacement Name',GotIt )
  IF(.NOT. GotIt) VarName = 'FlowLineDisp'

  CALL ListAddString( SolverParams,'Variable',VarName )
  CALL ListAddInteger( SolverParams,'Variable Dofs',dim )

END SUBROUTINE StructuredFlowLine_init
!------------------------------------------------------------------------------
