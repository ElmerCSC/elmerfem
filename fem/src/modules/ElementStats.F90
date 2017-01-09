!-----------------------------------------------------------------------------
!> A prototype solver for advection-diffusion-reaction equation,
!> This equation is generic and intended for education purposes
!> but may also serve as a starting point for more complex solvers.
!------------------------------------------------------------------------------
SUBROUTINE ElementStats( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
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
  TYPE(Element_t),POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: i, j, t, t0, n, m, is_bc
  INTEGER :: FirstElem, LastElem

  REAL(KIND=dp) :: WarnSize, MinSize, MaxSize, AveSize, dSize, h, &
      MaxAlpha, dAlpha, WarnAngle
  REAL(KIND=dp) :: alpha
  TYPE(ValueList_t), POINTER :: Params
  INTEGER :: ElementCount(1000)
  CHARACTER(MAX_NAME_LEN) :: Str
  LOGICAL :: Found, DoCohorts, DoingCohorts
  INTEGER :: NoCohorts
  INTEGER, ALLOCATABLE :: SizeCount(:),SkewCount(:)
!------------------------------------------------------------------------------
  
  CALL Info('ElementStats','------------------------------',Level=6)
  CALL Info('ElementStats','Studying elements of the mesh')

  Params => GetSolverParams()
  Mesh => GetMesh()

  ! Do cohorst only if requested
  DoCohorts = GetLogical( Params,'Show Cohorts',Found ) 
  IF( DoCohorts ) THEN
    NoCohorts = GetInteger( Params,'Number of Cohorts',Found ) 
    IF( .NOT. Found ) NoCohorts = 10 
    ALLOCATE( SizeCount(NoCohorts), SkewCount(NoCohorts) )
    SizeCount = 0
    SkewCount = 0
  END IF

  ! Go through bulk and boundary elements separately
  DO is_bc = 0, 1

    IF( is_bc == 0 ) THEN
      t0 = 0
      FirstElem = 1
      LastElem = Mesh % NumberOfBulkElements 
      Str = 'bulk'
    ELSE
      t0 = Mesh % NumberOfBulkElements
      FirstElem = 1
      LastElem = Mesh % NumberOfBoundaryElements 
      Str = 'boundary'
    END IF
    
    IF( FirstElem > LastElem ) CYCLE

    CALL Info('ElementStats','Going through '//TRIM(Str)//' elements')

    ! Initialize separately for bulk and boundary
    WarnSize = EPSILON( WarnSize ) 
    MinSize = HUGE( MinSize )
    MaxSize = 0.0_dp
    AveSize = 0.0_dp
    MaxAlpha = 0.0_dp ! by construction minimum is zero
    ElementCount = 0
    DoingCohorts = .FALSE.
    
100 DO t=FirstElem, LastElem
      
      Element => Mesh % Elements(t0+t)
      Model % CurrentElement => Element
      
      n = GetElementNOFNodes(Element)
      m = GetElementCode(Element)
            
      ElementCount(m) = ElementCount(m) + 1

      ! We could add additional info to this function, e.g. angle computation
      CALL ElementSizeStudy( Element, n, h )

      CALL ElementSkewStudy( Element, n, alpha ) 
      

      IF( .NOT. DoingCohorts ) THEN
        IF( h < WarnSize ) THEN
          PRINT *,'Element '//TRIM(I2S(t))//' of type '//TRIM(I2S(m))//' too small: ',h
        END IF
        
        MinSize = MIN( h, MinSize ) 
        MaxSize = MAX( h, MaxSize ) 
        AveSize = AveSize + h
        
        MaxAlpha = MAX( alpha, MaxAlpha ) 
      ELSE
        i = CEILING( NoCohorts * (h-MinSize)/(MaxSize-MinSize) )
        i = MAX( 1, MIN( i, NoCohorts ) )
        SizeCount(i) = SizeCount(i) + 1

        i = CEILING( NoCohorts * alpha/MaxAlpha )
        i = MAX( 1, MIN( i, NoCohorts ) )
        SkewCount(i) = SkewCount(i) + 1
      END IF
    END DO

    ! And finally print the info on the sizes
    IF( .NOT. DoingCohorts ) THEN
      AveSize = AveSize / (LastElem-FirstElem+1)
      
      PRINT *,'Average element size: ',AveSize
      PRINT *,'Minimum element size: ',MinSize
      PRINT *,'Maximum element size: ',MaxSize
      
      IF( MinSize > EPSILON( MinSize ) ) THEN
        PRINT *,'Ratio of element sizes: ',MaxSize / MinSize
      END IF
      PRINT *,'Average element size: ',AveSize

      PRINT *,'Maximum edge skew (degs): ',MaxAlpha
    ELSE
      PRINT *,'Size cohorts:'
      DO i=1,NoCohorts
        IF( SizeCount(i) > 0 ) THEN
          PRINT *,i,': ',SizeCount(i),' (',MinSize+(i-1)*dSize,' to ',MinSize+i*dSize,')'
        END IF
      END DO
      PRINT *,'Skew cohorts (degs):'
      DO i=1,NoCohorts
        IF( SkewCount(i) > 0 ) THEN
          PRINT *,i,': ',SkewCount(i),' (',(i-1)*dAlpha,' to ',i*dAlpha,')'
        END IF
      END DO
    END IF
    
    ! If we study cohorst we need a second sweep over the elements
    IF( DoCohorts .AND. .NOT. DoingCohorts ) THEN
      DoingCohorts = .TRUE.
      dSize = ( MaxSize - MinSize ) / NoCohorts
      dAlpha = MaxAlpha / NoCohorts
      GOTO 100
    END IF

  END DO

  CALL Info('ElementStats','------------------------------',Level=6)


CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE ElementSizeStudy( Element, n, ElemSize )
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: ElemSize
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Weight
    REAL(KIND=dp) :: Basis(n),DetJ
    LOGICAL :: Stat
    INTEGER :: t
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------


    CALL GetElementNodes( Nodes )
    IP = GaussPoints( Element )

    ElemSize = 0.0_dp

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )

      Weight = IP % s(t) * DetJ
      ElemSize = ElemSize + Weight
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE ElementSizeStudy
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ElementSkewStudy( Element, n, alpha )
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: alpha
!------------------------------------------------------------------------------
    INTEGER :: ElemFamily
    INTEGER, POINTER :: EdgeMap(:,:)
    INTEGER :: NoEdgePairs,i1,i2,j1,j2,k
    INTEGER :: EdgePairs(24)
    REAL(KIND=dp) :: Beta, MaxBeta, ri(3), rj(3), ri_len, rj_len
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    
    CALL GetElementNodes( Nodes )
    ElemFamily = GetElementFamily( Element )    
    EdgeMap => GetEdgeMap(ElemFamily)
    
    Alpha = 0.0_dp
    
    SELECT CASE ( ElemFamily )
      
    CASE(1,2,3,5)
      ! Always affine
      RETURN
      
    CASE(4) ! quads
      NoEdgePairs = 2
      EdgePairs(1:2*NoEdgePairs) = (/1,3,2,4/)
      ! I.e. the edge pairs to study are 1 & 3, and 2 & 4
      
    CASE(6) ! pyramids
      NoEdgePairs = 2
      EdgePairs(1:2*NoEdgePairs) = (/1,3,2,4/)
      
    CASE(7) ! wedges/prisms
      NoEdgePairs = 6
      EdgePairs(1:2*NoEdgePairs) = (/1,4,7,8,2,5,8,9,3,6,9,7/)
      
    CASE(8) ! hexas
      EdgePairs = 12
      EdgePairs(1:2*NoEdgePairs) = &
          (/1,3,2,4,5,7,6,8,1,5,9,10,2,6,10,11,3,7,11,12,4,8,12,9/)      
      
    END SELECT
    
    
    MaxBeta = 0.0_dp
    
    DO k=1,NoEdgePairs
      i1 = EdgeMap(EdgePairs(2*k-1),1)
      i2 = EdgeMap(EdgePairs(2*k-1),2)
      j1 = EdgeMap(EdgePairs(2*k),1)
      j2 = EdgeMap(EdgePairs(2*k),2)
      
      ri(1) = Nodes % x(i2) - Nodes % x(i1)
      ri(2) = Nodes % y(i2) - Nodes % y(i1)
      ri(3) = Nodes % z(i2) - Nodes % z(i1)
      
      rj(1) = Nodes % x(j2) - Nodes % x(j1)
      rj(2) = Nodes % y(j2) - Nodes % y(j1)
      rj(3) = Nodes % z(j2) - Nodes % z(j1)
      
      ri_len = SQRT( SUM( ri * ri ) )
      rj_len = SQRT( SUM( rj * rj ) )
      
      IF( ri_len * rj_len > 0 ) THEN
        beta = ACOS( SUM( ri * rj ) / ( ri_len * rj_len ) )
        beta = MIN( ABS( beta ), ABS( PI - beta ) )
        MaxBeta = MAX( beta, MaxBeta )
      END IF
    END DO

    ! Transform angle to degrees
    Alpha = 180 * MaxBeta / PI 
    
!------------------------------------------------------------------------------
  END SUBROUTINE ElementSkewStudy
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE ElementStats
!------------------------------------------------------------------------------
