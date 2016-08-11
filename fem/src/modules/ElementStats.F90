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

  REAL(KIND=dp) :: WarnSize, MinSize, MaxSize, AveSize, dSize, h
  TYPE(ValueList_t), POINTER :: Params
  INTEGER :: ElementCount(1000)
  CHARACTER(MAX_NAME_LEN) :: Str
  LOGICAL :: Found, DoCohorts, DoingCohorts
  INTEGER :: NoCohorts
  INTEGER, ALLOCATABLE :: SizeCount(:)
!------------------------------------------------------------------------------
  
  CALL Info('ElementStats','Studying elements of the mesh')

  Params => GetSolverParams()
  Mesh => GetMesh()

  ! Do cohorst only if requested
  DoCohorts = GetLogical( Params,'Show Cohorts',Found ) 
  IF( DoCohorts ) THEN
    NoCohorts = GetInteger( Params,'Number of Cohorts',Found ) 
    IF( .NOT. Found ) NoCohorts = 10 
    ALLOCATE( SizeCount(NoCohorts) )
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
    ElementCount = 0
    DoingCohorts = .FALSE.
    IF( DoCohorts ) THEN
      SizeCount = 0
    END IF

    
100 DO t=FirstElem, LastElem
      
      Element => Mesh % Elements(t0+t)
      Model % CurrentElement => Element
      
      n = GetElementNOFNodes(Element)
      m = GetElementCode(Element)
            
      ElementCount(m) = ElementCount(m) + 1

      ! We could add additional info to this function, e.g. angle computation
      CALL ElementStudy( Element, n, h )
      

      IF( .NOT. DoingCohorts ) THEN
        IF( h < WarnSize ) THEN
          PRINT *,'Element '//TRIM(I2S(t))//' of type '//TRIM(I2S(m))//' too small: ',h
        END IF
        
        MinSize = MIN( h, MinSize ) 
        MaxSize = MAX( h, MaxSize ) 
        AveSize = AveSize + h
      ELSE
        i = CEILING( NoCohorts * (h-MinSize)/(MaxSize-MinSize) )
        i = MAX( 1, MIN( i, NoCohorts ) )
        SizeCount(i) = SizeCount(i) + 1
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
    ELSE
      PRINT *,'Size cohorts:'
      DO i=1,NoCohorts
        IF( SizeCount(i) > 0 ) THEN
          PRINT *,i,': ',SizeCount(i),' (',MinSize+(i-1)*dSize,' to ',MinSize+i*dSize,')'
        END IF
      END DO
    END IF
    
    ! If we study cohorst we need a second sweep over the elements
    IF( DoCohorts .AND. .NOT. DoingCohorts ) THEN
      DoingCohorts = .TRUE.
      dSize = ( MaxSize - MinSize ) / NoCohorts
      GOTO 100
    END IF

  END DO



CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE ElementStudy( Element, n, ElemSize )
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
  END SUBROUTINE ElementStudy
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
END SUBROUTINE ElementStats
!------------------------------------------------------------------------------
