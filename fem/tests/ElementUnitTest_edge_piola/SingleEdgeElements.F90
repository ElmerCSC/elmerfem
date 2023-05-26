!------------------------------------------------------------------------------
SUBROUTINE SingleEdgeElements_init( Model,Solver,dt,TransientSimulation )
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
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: PiolaVersion, SecondOrder, Found
  
  SolverParams => GetSolverParams()  
  IF ( .NOT.ListCheckPresent(SolverParams, "Element") ) THEN
    SecondOrder = GetLogical( SolverParams, 'Quadratic Approximation', Found )  
    IF( SecondOrder ) THEN
      PiolaVersion = .TRUE.
    ELSE
      PiolaVersion = GetLogical(SolverParams, 'Use Piola Transform', Found )   
    END IF
    IF( SecondOrder ) THEN
      CALL ListAddString( SolverParams, "Element", &
          "n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2" )
    ELSE IF ( PiolaVersion ) THEN    
      CALL ListAddString( SolverParams, "Element", "n:0 e:1 -brick b:3 -quad_face b:2" )
    ELSE
      CALL ListAddString( SolverParams, "Element", "n:0 e:1" )
    END IF
  END IF

  
  CALL ListAddString( Solver % Values,'Variable',&
      '-nooutput -global dummy_var')

END SUBROUTINE SingleEdgeElements_init


!------------------------------------------------------------------------------
SUBROUTINE SingleEdgeElements( Model,Solver,dt,TransientSimulation )
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
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: nd, t, active
  INTEGER :: iter, maxiter
  LOGICAL :: Found
  INTEGER :: NoCycles, NoTest = 0, NoMode = 0
  REAL(KIND=dp) :: rt0, ct0, drt, dct, cumrt = 0.0_dp, cumct = 0.0_dp 
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t), POINTER :: pSolver
  
  SAVE :: NoTest, cumrt, cumct 
  
  NoTest = NoTest + 1
  Mesh => GetMesh()
  pSolver => Solver

  IF( NoTest > SIZE( Mesh % Elements ) ) THEN
    NoMode = NoMode + 1
    NoTest = 1
  END IF

  NoCycles = 0
  IF( NoMode == 0 ) THEN
    NoCycles = ListGetInteger( Solver % Values,'Timing Cycles',Found)
  END IF
      
  Element => Mesh % Elements(NoTest)
  PRINT *,'Element type: ',Element % TYPE % ElementCode

  ! This is just big enough. The ElementInfo internally knows the number
  ! of basis functions to return. 
  nd = 200
  
  ct0 = CPUTime(); rt0 = RealTime()
  CALL LocalMatrix( Element, nd )
  dct = CPUTime()-ct0; drt = RealTime()-rt0

  IF( NoCycles > 0 ) THEN
    PRINT *,'Time for element: ',&
        Element % TYPE % ElementCode, dct, drt
    cumct = cumct + dct
    cumrt = cumrt + drt

    IF( NoTest == SIZE( Mesh % Elements ) ) THEN
      PRINT *,'Time in total:',cumct, cumrt
    END IF
  END IF
  PRINT *,'Norm for element: ', Norm

  Solver % Variable % Values(1) = Norm
  
CONTAINS

! Assembly of the matrix entries arising from the bulk elements
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, nd )
!------------------------------------------------------------------------------
    INTEGER :: nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd)
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim,np
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes, UElement = Element )

    MASS  = 0._dp
    STIFF = 0._dp

    IP = GaussPointsAdapt( Element )
    dBasisdx = 0.0_dp
    

    DO t=1,IP % n
      Basis = 0.0_dp
      dBasisdx = 0.0_dp
      Wbasis = 0.0_dp
      RotWBasis = 0.0_dp
      
      DO i=0,NoCycles      
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
            RotBasis = RotWBasis, USolver = pSolver )
      END DO

      Weight = IP % s(t) * detJ      

      DO p = 1,nd
        DO q = 1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * SUM(RotWBasis(p,:)*RotWBasis(q,:))
          MASS(p,q) = MASS(p,q) + Weight * SUM(WBasis(p,:) * WBasis(q,:))
        END DO
      END DO
    END DO

    ! 1st time check the mass matrix, 2nd time stiffness matrix
    IF( MODULO(NoMode,2) == 0 ) THEN
      Norm = SUM(MASS(1:nd,1:nd))
    ELSE
      Norm = SUM(ABS(STIFF(1:nd,1:nd)))
    END IF
    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

END SUBROUTINE SingleEdgeElements
