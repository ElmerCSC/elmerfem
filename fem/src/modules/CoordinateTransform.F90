!/******************************************************************************
! *
! *  Module for calculating the rotation matrix for 2-rank tensors
! *
! *  Author: Eelis Takala
! *  Email:   eelis.takala@trafotek.fi
! *  Web:     http://www.trafotek.fi
! *  Address: Trafotek
! *           Kaarinantie 700
! *           20540 Turku
! *
! *  Original Date: 13.3.2014
! *
! *****************************************************************************/
 
!> \ingroup Solvers
!> \{
 
 
MODULE VecUtils
 
  USE DefUtils
  IMPLICIT NONE
 
  CONTAINS
 
    FUNCTION norm2(a)
      IMPLICIT NONE
      REAL(KIND=dp) :: norm2(3)
      REAL(KIND=dp), INTENT(IN) :: a(3)
      norm2 = sqrt(SUM(a**2))
    END FUNCTION norm2
 
    FUNCTION normalized(a)
      IMPLICIT NONE
      REAL(KIND=dp) :: normalized(3)
      REAL(KIND=dp), INTENT(IN) :: a(3)
      normalized = a/norm2(a)
    END FUNCTION normalized
 
    ! compute the jacobian of the coordinate transform between
    ! old coordinate system and the new coordinate system:
    ! --------------------------------------------------------
    FUNCTION jac(OldCoord, NewCoord)
      IMPLICIT NONE
 
      REAL(KIND=dp) :: OldCoord(3,3), NewCoord(3,3)
      REAL(KIND=dp) :: jac(3,3)
      INTEGER :: i, j
 
      DO i = 1, 3
        DO j = 1, 3
          jac(i,j) = DOT_PRODUCT(NewCoord(i,1:3), OldCoord(j,1:3)) &
               /sqrt(DOT_PRODUCT(OldCoord(j,1:3), OldCoord(j,1:3)))
        END DO
      END DO
 
    END FUNCTION jac
 
    ! here is how we transform a 3x3 2 rank tensor according 
    ! to the general formula. Btw. The smarter way is to do 
    ! MATMUL(jac, MATMUL(A, TRANSPOSE(jac))):
    ! ------------------------------------------------------
    FUNCTION transform2rank(A, jac) RESULT (B)
      IMPLICIT NONE
 
      REAL(KIND=dp) :: A(3,3), B(3,3)
      REAL(KIND=dp) :: jac(3,3)
      INTEGER :: i, j, k, l
 
      B = 0._dp
      DO i = 1, 3
        DO j = 1, 3
          DO k = 1, 3
            DO l = 1, 3
              B(i,j) = B(i,j) + jac(i,l)*jac(j,k) * A(k,l)
            END DO
          END DO
        END DO
      END DO
 
    END FUNCTION transform2rank
 
 
END MODULE VecUtils
 
!------------------------------------------------------------------------------
!> Initialization for the primary solver: CurrentSource
!------------------------------------------------------------------------------
SUBROUTINE RotMSolver_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
 
!  SolverParams => GetSolverParams()
!  CALL ListAddLogical( SolverParams, "Discontinuous Galerkin", .TRUE.)
 
!  CALL ListAddString( SolverParams, "Exported Variable 1", &
!             "Alpha Vector E[Alpha Vector E:3]" )
!  CALL ListAddString( SolverParams, "Exported Variable 2", &
!             "Beta Vector E[Beta Vector E:3]" )
!  CALL ListAddString( SolverParams, "Exported Variable 3", &
!             "Gamma Vector E[Gamma Vector E:3]" )
 
!------------------------------------------------------------------------------
END SUBROUTINE RotMSolver_init
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
SUBROUTINE RotMSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE VecUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: Found, First=.TRUE.
 
  TYPE(ValueList_t), POINTER :: BF,BodyParams,BC
 
  REAL(KIND=dp), ALLOCATABLE, SAVE :: RotM(:,:,:)
  REAL(KIND=DP), POINTER, SAVE :: alpha_ref(:,:), beta_ref(:,:), gamma_ref(:,:)
  TYPE(Mesh_t), POINTER :: Mesh  
  TYPE(Valuelist_t), POINTER :: Material
  REAL(KIND=dp) :: CoordSys_ijk(3,3), CoordSys_ref(3,3)
 
  LOGICAL  :: STAT
  INTEGER :: istat, n, nd, nn, q
 
  TYPE(Variable_t), POINTER, SAVE :: RotMvar, alphavecvar, &
                                     betavecvar, gammavecvar, tmpvar
 
  INTEGER, PARAMETER :: ind1(9) = [1,1,1,2,2,2,3,3,3]
  INTEGER, PARAMETER :: ind2(9) = [1,2,3,1,2,3,1,2,3]
 
!------------------------------------------------------------------------------
 
  IF (First) THEN
    First = .FALSE.
 
    Mesh => Model % Mesh
    N = Mesh % MaxElementDOFs
 
    ALLOCATE( RotM(3,3,N), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'RotMSolver', 'Memory allocation error.' )
    END IF
 
    NULLIFY( alpha_ref )
    NULLIFY( beta_ref )
 
    RotMvar => VariableGet( Mesh % Variables, 'RotM E')
    IF(.NOT. ASSOCIATED(RotMVar)) THEN
      CALL Fatal('RotMSolver()','RotM E variable not found')
    END IF
 
    alphavecvar => VariableGet( Mesh % Variables, 'Alpha Vector E')
    IF(.NOT. ASSOCIATED(alphavecvar)) THEN
      CALL Fatal('RotMSolver()','Alpha Vector E variable not found')
    END IF
 
    betavecvar => VariableGet( Mesh % Variables, 'Beta Vector E')
    IF(.NOT. ASSOCIATED(betavecvar)) THEN
      CALL Fatal('RotMSolver()','Beta Vector E variable not found')
    END IF
 
    gammavecvar => VariableGet( Mesh % Variables, 'Gamma Vector E')
    IF(.NOT. ASSOCIATED(gammavecvar)) THEN
      CALL Fatal('RotMSolver()','Gamma Vector E variable not found')
    END IF
 
    CoordSys_ijk(1,:) = [1,0,0]
    CoordSys_ijk(2,:) = [0,1,0]
    CoordSys_ijk(3,:) = [0,0,1]
 
  END IF
 
  DO q=1,GetNOFActive()
    Element => GetActiveElement(q)
    nn = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
 
    BodyParams => GetBodyParams( Element )
    IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('RotMSolver', 'Body Parameters not found')
 
    CALL GetConstRealArray(BodyParams, alpha_ref, 'Alpha Reference', Found)
    IF (.NOT. Found) CYCLE
    IF (SIZE(alpha_ref,1) /= 3) CALL Fatal('RotMSolver','Alpha Reference should have three components!')
    CoordSys_ref(1,1:3) = normalized(alpha_ref(1:3,1))
 
    CALL GetConstRealArray(BodyParams, beta_ref, 'Beta Reference', Found)
    IF (.NOT. Found) CYCLE
    IF (SIZE(beta_ref,1) /= 3) CALL Fatal('RotMSolver','Beta Reference should have three components!')
    CoordSys_ref(2,1:3) = normalized(beta_ref(1:3,1))
 
    CoordSys_ref(3,1:3) = normalized(crossproduct(CoordSys_ref(1,1:3), CoordSys_ref(2,1:3)))
 
    ! Compute the rotation matrices for 2 rank tensors for all the elements
    ! (transform from the element local coordinate system to the ijk)
    ! ----------------------------------------------------------------
 
    CALL ComputeRotM(Element, CoordSys_ijk, CoordSys_ref, nn, nd)
  END DO
 
CONTAINS
 
!------------------------------------------------------------------------------
   SUBROUTINE ComputeRotM(Element,CoordSys_ijk,CoordSys_ref,nn,nd)
!------------------------------------------------------------------------------
    INTEGER :: nn, nd, ind
    TYPE(Element_t) :: Element
 
    REAL(KIND=dp) :: un,vn,wn,Basis(nn), DetJ,localC,gradv(3)
    REAL(KIND=dp) :: dBasisdx(nn,3),Coordsys(3,3),Coordsys2(3,3), &
                     Relem(3,3), jac_ref(3,3), jac_e(3,3), alpha(nn), &
                     beta(nn), tvec(3)
    REAL(KIND=dp) :: CoordSys_ijk(3,3), CoordSys_ref(3,3)
 
    INTEGER :: j,Indexes(nd),l,m,k
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
 
    CALL GetElementNodes(Nodes)
    tmpvar => VariableGet( Mesh % Variables, 'alpha')
    IF(ASSOCIATED(tmpvar)) THEN
      CALL GetLocalSolution(alpha,'alpha')
    ELSE
      CALL GetLocalSolution(alpha,'alpha direction')
    END IF

    tmpvar => VariableGet( Mesh % Variables, 'beta')
    IF(ASSOCIATED(tmpvar)) THEN
      CALL GetLocalSolution(beta,'beta')
    ELSE
      CALL GetLocalSolution(beta,'beta direction')
    END IF

    DO j=1,nn
       un = Element % TYPE % NodeU(j)
       vn = Element % TYPE % NodeV(j)
       wn = Element % TYPE % NodeW(j)
 
      stat=ElementInfo(Element,Nodes,un,vn,wn,detJ,Basis,dBasisdx)
 
      ! for anisotropic conductivity, determine the direction
      ! CoordSys(3,1:3) is the winding direction
      ! CoordSys(1,1:3) is the perpendicular direction to the foil 
      ! surface
      ! -----------------------------------------------------
      CoordSys(1,1:3) = normalized(MATMUL( alpha(1:nn), dBasisdx(1:nn,:)))
      CoordSys(2,1:3) = normalized(MATMUL( beta(1:nn), dBasisdx(1:nn,:)))
      CoordSys(3,1:3) = normalized(crossproduct(CoordSys(1,1:3), CoordSys(2,1:3)))
 
      CoordSys2 = CoordSys
 
      ! jacobian matrix (transform from the reference coordinate 
      ! system to the ijk):
      ! --------------------------------------------------------------
      jac_ref = jac(CoordSys_ref, CoordSys_ijk)
 
      ! Compute the jacobian for the coordinate transformation from 
      ! the local coordinate system to the reference coordinate system:
      ! ---------------------------------------------------------------
      jac_e = jac(CoordSys, CoordSys_ref)
 
      ! Now the transformation matrix from the 
      ! local coordinate system to the ijk:
      ! -----------------------------------
 
      Relem = MATMUL(jac_ref, jac_e)
 
      CoordSys2(1,1:3) = [1,0,0]
      CoordSys2(2,1:3) = [0,1,0]
      CoordSys2(3,1:3) = [0,0,1]
 
!      CoordSys2 = MATMUL(MATMUL(Relem, CoordSys2), TRANSPOSE(Relem))
 
      CoordSys2(1,1:3) = MATMUL(Relem, Coordsys2(1,1:3))
      CoordSys2(2,1:3) = MATMUL(Relem, Coordsys2(2,1:3))
      CoordSys2(3,1:3) = MATMUL(Relem, Coordsys2(3,1:3))
!      CoordSys2 = CoordSys
 
 
      IF (ASSOCIATED(RotMvar)) THEN
        DO k=1,RotMvar % DOFs
          RotMvar % Values( RotMvar % DOFs*(RotMvar % Perm( &
                Element % DGIndexes(j))-1)+k) = Relem(ind1(k),ind2(k))
        END DO
      END IF
 
      IF (ASSOCIATED(alphavecvar)) THEN
        DO k=1,alphavecvar % DOFs
          alphavecvar % Values( alphavecvar % DOFs*(alphavecvar % Perm( &
                Element % DGIndexes(j))-1)+k) = CoordSys2(1,k)
        END DO
      END IF
 
      IF (ASSOCIATED(betavecvar)) THEN
        DO k=1,betavecvar % DOFs
          betavecvar % Values( betavecvar % DOFs*(betavecvar % Perm( &
                Element % DGIndexes(j))-1)+k) = CoordSys2(2,k)
        END DO
      END IF
 
      IF (ASSOCIATED(gammavecvar)) THEN
        DO k=1,gammavecvar % DOFs
          gammavecvar % Values( gammavecvar % DOFs*(gammavecvar % Perm( &
                Element % DGIndexes(j))-1)+k) = CoordSys2(3,k)
        END DO
      END IF
 
    END DO
 
!------------------------------------------------------------------------------
   END SUBROUTINE ComputeRotM
!------------------------------------------------------------------------------
 
 
!------------------------------------------------------------------------------
   SUBROUTINE GetElementRotM(Element,RotM,nn)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: Element
     INTEGER :: k, l, m, j, nn
     REAL(KIND=dp) :: RotM(3,3,nn)
 
     RotM = 0._dp
 
     DO j = 1, nn
       IF (ASSOCIATED(RotMvar)) THEN
         DO k=1,RotMvar % DOFs
           RotM(ind1(k),ind2(k),j) = RotMvar % Values( &
                 RotMvar % DOFs*(RotMvar % Perm(Element % DGIndexes(j))-1)+k)
         END DO
       ELSE 
         CALL Fatal ('GetElementRotM', 'RotM E is not associated.')
       END IF
 
     END DO
 
!------------------------------------------------------------------------------
   END SUBROUTINE GetElementRotM
!------------------------------------------------------------------------------
 
 
!------------------------------------------------------------------------------
END SUBROUTINE RotMSolver
!------------------------------------------------------------------------------
