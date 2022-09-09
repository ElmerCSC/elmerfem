!------------------------------------------------------------------------------
SUBROUTINE SingleElements_init( Model,Solver,dt,TransientSimulation )
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

  CALL ListAddString( Solver % Values,'Variable',&
      '-nooutput -global dummy_var')

END SUBROUTINE SingleElements_init


!------------------------------------------------------------------------------
SUBROUTINE SingleElements( Model,Solver,dt,TransientSimulation )
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
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, active
  INTEGER :: iter, maxiter
  LOGICAL :: Found
  INTEGER :: NoCycles, NoTest = 0, NoMode = 0
  REAL(KIND=dp) :: rt0, ct0, drt, dct, cumrt = 0.0_dp, cumct = 0.0_dp 
  TYPE(Mesh_t), POINTER :: Mesh

  SAVE :: NoTest, cumrt, cumct 
  
  NoTest = NoTest + 1
  Mesh => GetMesh()

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

  n  = GetElementNOFNodes(Element)
  !nd = GetElementNOFDOFs(Element)
  !nb = GetElementNOFBDOFs(Element)

  nd = n
  nb = 0
  
  !PRINT *,'Ns for element: ',n, nd, nb
  
  ct0 = CPUTime(); rt0 = RealTime()
  CALL LocalMatrix( Element, n, nd+nb )
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
  SUBROUTINE LocalMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
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
      DO i=0,NoCycles      
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )
      END DO
      Weight = IP % s(t) * detJ      
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
          MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )
      DO p=1,nd
        DO q=1,nd
          MASS(p,q) = MASS(p,q) + Weight * Basis(q) * Basis(p)
        END DO
      END DO
    END DO

    ! 1st time check the mass matrix, 2nd time stiffness matrix
    IF( NoMode == 0 ) THEN
      Norm = SUM(MASS(1:nd,1:nd))
    ELSE
      Norm = SUM(ABS(STIFF(1:nd,1:nd)))
    END IF
    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

END SUBROUTINE SingleElements
