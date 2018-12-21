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
! *  Authors: Juhani Kataja, Peter Råback, Juha Ruokolainen and Mika Malinen
! *  Email:   juhani.kataja@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 25 Aug 2018
! * 
! *  Heavily inspired from the MagnetoDynamics and VectorHelmholtz modules.
! *****************************************************************************/
    
!------------------------------------------------------------------------------
!>  Solve time-dependent Maxwell equations using the curl-curl equation 
!>  using curl-conforming edge elements.
!> \ingroup Solvers
!-------------------------------------------------------------------------------
MODULE EMWaveSolverUtils

   USE DefUtils

CONTAINS


!------------------------------------------------------------------------------
  FUNCTION GetBoundaryEdgeIndex(Boundary,nedge) RESULT(n)
!------------------------------------------------------------------------------
    INTEGER :: n,nedge
    TYPE(Element_t) :: Boundary
!------------------------------------------------------------------------------
    INTEGER :: i,j,jb1,jb2,je1,je2
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Parent, Edge, Face
!------------------------------------------------------------------------------
    Mesh => GetMesh()
    n = 0
    SELECT CASE(GetElementFamily(Boundary))
    CASE(1)
      RETURN
    CASE(2)
      IF ( nedge==1 ) THEN
        Parent => Boundary % BoundaryInfo % Left
        IF ( .NOT. ASSOCIATED(Parent) ) &
            Parent => Boundary % BoundaryInfo % Right
 
        jb1 = Boundary % NodeIndexes(1)
        jb2 = Boundary % NodeIndexes(2)
        DO i=1,Parent % TYPE % NumberOfEdges
          Edge => Mesh % Edges(Parent % EdgeIndexes(i))
          je1 = Edge % NodeIndexes(1)
          je2 = Edge % NodeIndexes(2)
          IF ( jb1==je1.AND.jb2==je2 .OR. jb1==je2.AND.jb2==je1) EXIT
        END DO
        n = Parent % EdgeIndexes(i)
      END IF
    CASE(3,4)
      j = GetBoundaryFaceIndex(Boundary)
      Face => Mesh % Faces(j)
      IF ( nedge>0.AND.nedge<=Face % TYPE % NumberOfEdges ) &
        n = Face % EdgeIndexes(nedge) 
    END SELECT
!------------------------------------------------------------------------------
  END FUNCTION GetBoundaryEdgeIndex
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION GetBoundaryFaceIndex(Boundary) RESULT(n)
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Element_t) :: Boundary
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,m
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Parent, Face
!------------------------------------------------------------------------------
    Mesh => GetMesh()
    Parent => Boundary % BoundaryInfo % Left
    IF ( .NOT. ASSOCIATED(Parent) ) &
       Parent => Boundary % BoundaryInfo % Right

    DO i=1,Parent % TYPE % NumberOfFaces
      Face => Mesh % Faces(Parent % FaceIndexes(i))
      m = 0
      DO j=1,Face % TYPE % NumberOfNodes
        DO k=1,Boundary % TYPE % NumberOfNodes
          IF ( Face % NodeIndexes(j)==Boundary % NodeIndexes(k)) m=m+1
        END DO
      END DO
      IF ( m==Boundary % TYPE % NumberOfNodes) EXIT
    END DO
    n = Parent % FaceIndexes(i)
!------------------------------------------------------------------------------
  END FUNCTION GetBoundaryFaceIndex
!------------------------------------------------------------------------------

  
END MODULE EMWaveSolverUtils



!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE EMWaveSolver_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE EMWaveSolverUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: Found, SecondOrder, PiolaVersion
  REAL(KIND=dp) :: mu0, eps0
  INTEGER :: mat_id
  TYPE(ValueList_t), POINTER  :: List
  
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
    ELSE IF (PiolaVersion) THEN    
      CALL ListAddString( SolverParams, "Element", "n:0 e:1 -brick b:3 -quad_face b:2" )
    ELSE
      CALL ListAddString( SolverParams, "Element", "n:0 e:1" )
    END IF
  END IF

  ! Use by some solvers e.g. SaveLine to aknowledge E as edge field
  CALL ListAddNewLogical( SolverParams,'Hcurl Basis',.TRUE.)
  IF( ListGetLogical( SolverParams,'Constant Bulk Matrix',Found ) ) THEN
    CALL ListAddNewLogical( SolverParams,'Use Global Mass Matrix',.TRUE.)    
  END IF
  
  CALL ListAddNewLogical( SolverParams,'Variable Output',.FALSE.)
  CALL ListAddNewString( SolverParams,'Variable','E')
  CALL ListAddNewLogical( SolverParams,'Linear System Complex', .FALSE.)
  
  CALL ListAddInteger( SolverParams,'Time derivative order', 2 )

  ! Set a multiplier for the relative keywords
  !--------------------------------------------------------------------
  eps0 = GetConstReal( Model % Constants,'Permittivity of Vacuum', Found )
  IF(.NOT. Found ) CALL Fatal('EMWaveSolver_Init0','> Permittivity of Vacuum < is required')

  mu0 = GetConstReal( Model % Constants,'Permeability of Vacuum', Found )
  IF(.NOT. Found ) CALL Fatal('EMWaveSolver_Init0','> Permeability of Vacuum < is required')

  ! does not seem to work? 
  ! the idea is that relative values would be automatically replaced by absolute ones with scaling
  !DO mat_id = 1, Model % NumberOfMaterials
  !  List => Model % Materials(mat_id) % Values  
  !  CALL ListSetCoefficients( list,'Relative Permittivity', eps0 )
  !  CALL ListSetCoefficients( list,'Relative Permeability', mu0 )
  !END DO
  
!------------------------------------------------------------------------------
END SUBROUTINE EMWaveSolver_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solve the electric field E from the rot-rot equation
! 
!> rot (1/mu_r) rot E - \kappa_0^2 epsilon_r E = i omega mu_0 J
!
!>  using edge elements (Nedelec/W basis of lowest degree) 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE EMWaveSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE EMWaveSolverUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: BC
  INTEGER :: n,istat,i,nNodes,Active,dofs
  INTEGER :: NoIterationsMax
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp) :: Norm
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), MASS(:,:), DAMP(:,:), FORCE(:)
  LOGICAL :: PiolaVersion, SecondOrder, EdgeBasis
  INTEGER, POINTER :: Perm(:)
  TYPE(ValueList_t), POINTER :: SolverParams
  REAL(KIND=dp) :: mu0, eps0
  TYPE(Solver_t), POINTER :: pSolver 
  
  SAVE STIFF, DAMP, MASS, FORCE, AllocationsDone
!------------------------------------------------------------------------------

  CALL Info('VectorHelmholztSolver','Solving electromagnetic waves in time',Level=5)

  SolverParams => GetSolverParams()

  SecondOrder = GetLogical( SolverParams, 'Quadratic Approximation', Found )  
  IF( SecondOrder ) THEN
    PiolaVersion = .TRUE.
  ELSE
    PiolaVersion = GetLogical( SolverParams,'Use Piola Transform', Found )
  END IF

  dofs = Solver % Variable % Dofs

  eps0 = GetConstReal( Model % Constants,'Permittivity of Vacuum')
  mu0 = GetConstReal( Model % Constants,'Permeability of Vacuum')
  
  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  Mesh => GetMesh()
  nNodes = Mesh % NumberOfNodes
  Perm => Solver % Variable % Perm
  pSolver => Solver
  
  IF ( .NOT. AllocationsDone ) THEN
    IF( dofs /= 1 ) CALL Fatal ('EMWaveSolver', 'Invalid variable size:'//TRIM(I2S(dofs)) )
    n = Mesh % MaxElementDOFs  
    ALLOCATE( FORCE(n), STIFF(n,n), MASS(n,n), DAMP(n,n), STAT=istat )
    IF ( istat /= 0 ) CALL Fatal( 'EMWaveSolver', 'Memory allocation error.' )
    AllocationsDone = .TRUE.
  END IF

  ! Resolve internal non.linearities, if requeted:
  ! ----------------------------------------------
  NoIterationsMax = GetInteger( SolverParams, &
      'Nonlinear System Max Iterations',Found)
  IF(.NOT. Found) NoIterationsMax = 1

  EdgeBasis = .NOT. ListCheckPresent( SolverParams,'Linear System Refactorize' ) .AND. &
      GetLogical( SolverParams, 'Edge Basis', Found )
  
  CALL DefaultStart()
  
  DO i=1,NoIterationsMax
    CALL DoBulkAssembly()

    CALL DoBoundaryAssembly()
    
    ! Default routines for finisging assembly and solving the system
    Norm = DefaultSolve()
    IF( DefaultConverged() ) EXIT

    IF( EdgeBasis ) CALL ListAddLogical( SolverParams,'Linear System Refactorize',.FALSE.)
  END DO
  IF ( EdgeBasis ) CALL ListRemove( SolverParams, 'Linear System Refactorize' )

  CALL DefaultFinish()
  
  CALL Info('VectorHelmholztSolver','All done',Level=10)
  

  
CONTAINS

!---------------------------------------------------------------------------------------------
  SUBROUTINE DoBulkAssembly()
!---------------------------------------------------------------------------------------------
    INTEGER :: n,nd,t
    LOGICAL :: Found, ConstantBulkInUse = .FALSE.
!---------------------------------------------------------------------------------------------
        
    ! use matrix from previous round
    !-----------------------------------------------
    IF( ConstantBulkInUse ) THEN
      CALL DefaultInitialize(UseConstantBulk = ConstantBulkInUse )
      RETURN
    END IF
    
    CALL DefaultInitialize()
   
    Active = GetNOFActive()
    DO t=1,active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes() ! nodes
      nd = GetElementNOFDOFs()  ! dofs

      ! Get element local matrix and rhs vector:
      !----------------------------------------
      CALL LocalMatrix( MASS, DAMP, STIFF, FORCE, &
          Element, n, nd, PiolaVersion, t==1 )

      ! Update global matrix and rhs vector from local matrix & vector:
      !---------------------------------------------------------------       
      CALL Default2ndOrderTime( MASS, DAMP, STIFF, FORCE )
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO

    CALL DefaultFinishBulkAssembly()

    ConstantBulkInUse = ListGetLogical( SolverParams,'Constant Bulk Matrix',Found )       

!------------------------------------------------------------------------------
  END SUBROUTINE DoBulkAssembly
!------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
  SUBROUTINE DoBoundaryAssembly()
!---------------------------------------------------------------------------------------------
    INTEGER :: n,nd,t,k
    LOGICAL :: Found, InitHandles
!---------------------------------------------------------------------------------------------
    TYPE(element_t), POINTER :: dummy_element

    InitHandles = .TRUE.

    ! Robin type of BC in terms of H:
    !--------------------------------
    Active = GetNOFBoundaryElements()
    DO t=1,Active
      Element => GetBoundaryElement(t)
      BC => GetBC()
      IF (.NOT. ASSOCIATED(BC) ) CYCLE

      SELECT CASE(GetElementFamily())
      CASE(1)
        CYCLE
      CASE(2)
        k = GetBoundaryEdgeIndex(Element,1); Element => Mesh % Edges(k)
      CASE(3,4)
        k = GetBoundaryFaceIndex(Element)  ; Element => Mesh % Faces(k)
      END SELECT
      IF (.NOT. ActiveBoundaryElement(Element)) CYCLE

      nd = GetElementNOFDOFs(Element)
      n  = GetElementNOFNodes(Element)
        
      dummy_element => SetCurrentElement(Element)
      
      CALL LocalMatrixBC(MASS,DAMP,STIFF,FORCE,&
          Element,n,nd,PiolaVersion,InitHandles)
      
      CALL Default2ndOrderTimeR( MASS, DAMP, STIFF, FORCE(1:nd), UElement=Element)
      CALL DefaultUpdateEquationsR(STIFF,FORCE(1:nd), UElement=Element)

      InitHandles = .FALSE.
    END DO
      
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()   

!------------------------------------------------------------------------------
  END SUBROUTINE DoBoundaryAssembly
!------------------------------------------------------------------------------


  
!-----------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( MASS, DAMP, STIFF, FORCE, &
      Element, n, nd, PiolaVersion, InitHandles )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), DAMP(:,:), FORCE(:), MASS(:,:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: PiolaVersion, InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: DetJ, weight, eps, mu, muinv, L(3), cond
    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3),WBasis(nd,3),RotWBasis(nd,3)
    LOGICAL :: Stat
    INTEGER :: t, i, j
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t) :: CD_h(3), Mu_h, Eps_h, Cond_h
    LOGICAL :: AllocationsDone = .FALSE.
    
    SAVE Cd_h, Mu_h, Eps_h, Cond_h, AllocationsDone

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Cd_h(1),'Body Force','Current Density 1')
      CALL ListInitElementKeyword( Cd_h(2),'Body Force','Current Density 2')
      CALL ListInitElementKeyword( Cd_h(3),'Body Force','Current Density 3')

      ! These have been normalized by mu0 and eps0 in _init section
      CALL ListInitElementKeyword( Mu_h,'Material','Relative Permeability')
      CALL ListInitElementKeyword( Eps_h,'Material','Relative Permittivity')

      CALL ListInitElementKeyword( Cond_h,'Material','Electric Conductivity')
    END IF

    CALL GetElementNodes( Nodes, Element )

    STIFF = 0.0_dp
    DAMP = 0.0_dp
    FORCE = 0.0_dp
    MASS  = 0.0_dp
       
    ! Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis = .TRUE., PReferenceElement = PiolaVersion)

    
    DO t=1,IP % n
      stat = ElementInfo(Element,Nodes,IP % u(t), IP % v(t), IP % w(t),detJ,Basis,dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 
      
      weight = detJ * IP%s(t)
      
      eps = eps0 * ListGetElementReal( eps_h, Basis, Element, Found )
      mu = mu0 * ListGetElementReal( mu_h, Basis, Element, Found )
      muinv = 1.0_dp / mu
      
      cond = ListGetElementReal( cond_h, Basis, Element, Found ) 
      
      L(1) = ListGetElementReal( CD_h(1), Basis, Element, Found )
      L(2) = ListGetElementReal( CD_h(2), Basis, Element, Found )
      L(3) = ListGetElementReal( CD_h(3), Basis, Element, Found )
      
      ! Compute element stiffness matrix and force vector:
      ! --------------------------------------------------
      DO i = 1,nd
        FORCE(i) = FORCE(i) - SUM(L*WBasis(i,:)) * weight
        
        DO j = 1,nd
          ! the mu^-1 curl u . curl v 
          STIFF(i,j) = STIFF(i,j) + muinv * &
              SUM(RotWBasis(i,:) * RotWBasis(j,:)) * weight
          
          ! the term d^2 ( \epsilon u.v ) / dt^2
          MASS(i,j) = MASS(i,j) +  &
              eps * SUM(WBasis(j,:) * WBasis(i,:)) * weight

          ! the term d ( \rho u.v ) / dt for conducting materials
          DAMP(i,j) = DAMP(i,j) +  &
              cond * SUM(WBasis(j,:) * WBasis(i,:)) * weight
        END DO
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( MASS, DAMP, STIFF, FORCE, Element, n, nd, PiolaVersion, InitHandles )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: MASS(:,:), DAMP(:,:), STIFF(:,:), FORCE(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: PiolaVersion, InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: DetJ,Normal(3),Tem(n)
    REAL(KIND=dp) :: B, L(3), muinv, mu, weight, tanWBasis(3)
    REAL(KIND=dp) :: Basis(nd), WBasis(nd,3), RotWBasis(nd,3),dBasisdx(nd,3) 
    LOGICAL :: Stat
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, p, q
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t) :: BL_h(3), Damp_h, Mu_h, Eps_h
    LOGICAL :: Visited = .FALSE., GotTem
    TYPE(Element_t), POINTER :: Parent
    LOGICAL :: AllocationsDone = .FALSE.
    
    SAVE Visited, BL_h, Damp_h, Mu_h, Eps_h

!------------------------------------------------------------------------------

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Bl_h(1),'Boundary Condition','Magnetic Boundary Load 1')
      CALL ListInitElementKeyword( Bl_h(2),'Boundary Condition','Magnetic Boundary Load 2')
      CALL ListInitElementKeyword( Bl_h(3),'Boundary Condition','Magnetic Boundary Load 3')

      CALL ListInitElementKeyword( Damp_h,'Boundary Condition','Electric Damping Coefficient')

      CALL ListInitElementKeyword( Mu_h,'Material','Relative Permeability')
      CALL ListInitElementKeyword( Eps_h,'Material','Relative Permittivity')
    END IF

    BC => GetBC() 
    Parent => GetBulkElementAtBoundary(Element) 

    TEM(1:n) = GetReal( BC,'TEM Potential', GotTem )

    CALL GetElementNodes( Nodes, Element )

    DAMP = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp
    MASS  = 0.0_dp

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion)

    DO t=1,IP % n
      stat = ElementInfo(Element,Nodes,IP % u(t), IP % v(t), IP % w(t),detJ,Basis,dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 
      
      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      weight = detJ * IP%s(t)

      mu = mu0 * ListGetElementReal( mu_h, Basis, Parent )
      muinv = 1.0_dp / mu

      L(1) = ListGetElementReal( Bl_h(1), Basis, Element, Found )
      L(2) = ListGetElementReal( Bl_h(2), Basis, Element, Found ) 
      L(3) = ListGetElementReal( Bl_h(3), Basis, Element, Found ) 

      ! We don't yet have a method for getting grad at ip
      IF( GotTem ) THEN
        L = L + MATMUL( Tem(1:n), dBasisdx(1:n,1:3) )
      END IF
      DO i = 1,nd
        tanWBasis(1:3) = WBasis(i,:) - Normal(1:3)*sum(Normal(1:3) * WBasis(i,:))
        FORCE(i) = FORCE(i) + muinv * sum(L(1:3) * tanWBasis(1:3)) * weight
      END DO

      B = ListGetElementReal( Damp_h, Basis, Element, Found ) 
      IF( Found ) THEN
        DO i = 1,nd
          tanWBasis(:) = WBasis(i,:) - Normal * SUM(Normal * WBasis(i,:))
          DO j = 1,nd
            DAMP(i,j) = DAMP(i,j) + muinv * B * SUM(tanWBasis(:) * WBasis(j,:)) * weight
          END DO
        END DO
      END IF
   END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE EMWaveSolver
!------------------------------------------------------------------------------

!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE EMWaveCalcFields_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE EMWaveSolverUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: sname,pname
  LOGICAL :: Found, ElementalFields
  INTEGER, POINTER :: Active(:)
  INTEGER :: mysolver,i,j,k,n,m,vDOFs,soln
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: Solvers(:), PSolver

  SolverParams => GetSolverParams()

  ! Find the solver index of the primary solver by the known procedure name.
  ! (the solver is defined here in the same module so not that dirty...)
  soln = 0
  
  DO i=1,Model % NumberOfSolvers
    sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
    j = INDEX( sname,'EMWaveSolver')
    IF( j > 0 ) THEN
      soln = i 
      EXIT
    END IF
  END DO
     
  IF( soln == 0 ) THEN
    CALL Fatal('EMWaveCalcFields_Init0','Cannot locate the primary solver: '//TRIM(I2S(soln)))      
  ELSE
    CALL Info('EMWaveCalcFields_Init0','The primary solver index is: '//TRIM(I2S(soln)),Level=12)
    CALL ListAddInteger( SolverParams,'Primary Solver Index',soln ) 
  END IF

  ! In case we are solving truly discontinuous Galerkin fields then we do it by assemblying
  ! normal linear system. Here we allocate for the DG type of fields that are computed elementwise
  ! while the FE fields are solved using standard Galerkin. Hence unintuitively we exit here
  ! for elemental fields if the primary field is itself elemenetal.
  !-----------------------------------------------------------------------------------------------
  IF( GetLogical(SolverParams,'Discontinuous Galerkin',Found)) RETURN

  ElementalFields = GetLogical( SolverParams, 'Calculate Elemental Fields', Found)
  IF(Found .AND. .NOT. ElementalFields) RETURN

  PSolver => Solver
  DO mysolver=1,Model % NumberOfSolvers
    IF ( ASSOCIATED(PSolver,Model % Solvers(mysolver)) ) EXIT
  END DO

  ! Here we add a DG solver instance in a dirty way by extending the list of solvers
  ! and adding the new solver as an active one. 
  n = Model % NumberOfSolvers
  DO i=1,Model % NumberOFEquations
    Active => ListGetIntegerArray(Model % Equations(i) % Values, &
        'Active Solvers', Found)
    m = SIZE(Active)
    IF ( ANY(Active==mysolver) ) &
        CALL ListAddIntegerArray( Model % Equations(i) % Values,  &
        'Active Solvers', m+1, [Active, n+1] )
  END DO

  ALLOCATE(Solvers(n+1))
  Solvers(1:n) = Model % Solvers
  SolverParams => NULL()
  CALL ListAddLogical( SolverParams, 'Discontinuous Galerkin', .TRUE. )
  Solvers(n+1) % DG = .TRUE.
  Solvers(n+1) % Values => SolverParams
  Solvers(n+1) % PROCEDURE = 0
  Solvers(n+1) % ActiveElements => NULL()
  CALL ListAddString( SolverParams, 'Exec Solver', 'never' )
  CALL ListAddLogical( SolverParams, 'No Matrix',.TRUE.)
  CALL ListAddString( SolverParams, 'Equation', 'never' )
  CALL ListAddString( SolverParams, 'Procedure', &
              'AllocateSolver AllocateSolver',.FALSE. )
  CALL ListAddString( SolverParams, 'Variable', '-nooutput cf_dummy' )

  pname = ListGetString( Model % Solvers(soln) % Values, 'Mesh', Found )
  IF(Found) THEN
    CALL ListAddString( SolverParams, 'Mesh', pname )
  END IF
  
  ! Electric field is always computed
  CALL ListAddString( SolverParams,&
      NextFreeKeyword('Exported Variable', SolverParams), &
      "Elfied E[Elfied E:3]");

  ! When requested we may also compute the 1st and 2nd time derivative.
  ! They exist by default as whitney fields but also needs to be projected
  ! to nodal fields.
  !------------------------------------------------------------------------------
  IF (GetLogical(SolverParams,'Calculate Electric field derivatives',Found)) THEN
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable', SolverParams), &
        "dEdt E[dEdt E:3]");
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable', SolverParams), &
        "ddEddt E[ddEddt E:3]");    
  END IF
  
  DEALLOCATE(Model % Solvers)
  Model % Solvers => Solvers
  Model % NumberOfSolvers = n+1
!------------------------------------------------------------------------------
END SUBROUTINE EMWaveCalcFields_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE EMWaveCalcFields_Init(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE EMWaveSolverUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  LOGICAL :: Found, NodalFields
  TYPE(ValueList_t), POINTER :: SolverParams

  SolverParams => GetSolverParams()

  ! We compute the fields one component at a time.
  ! This is the dummy variable used for the computation. It is not saved. 
  CALL ListAddString( SolverParams, 'Variable', '-nooutput hr_dummy' )

  ! The matrix is constant hence do not ever refactorize.
  CALL ListAddLogical( SolverParams, 'Linear System refactorize', .FALSE.)
  
  NodalFields = GetLogical( SolverParams, 'Calculate Nodal Fields', Found)
  IF(Found .AND. .NOT. NodalFields ) RETURN
  
  CALL ListAddString( SolverParams,&
      NextFreeKeyword('Exported Variable', SolverParams), &
      "Elfield[Elfield:3]");
  
  IF (GetLogical(SolverParams,'Calculate Electric field derivatives',Found)) THEN
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable', SolverParams), &
        "dEdt[dEdt:3]");
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable', SolverParams), &
        "ddEddt[ddEddt:3]");    
  END IF
  
!------------------------------------------------------------------------------
END SUBROUTINE EMWaveCalcFields_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Calculate fields resulting from the edge element formulation 
!> \ingroup Solvers
!------------------------------------------------------------------------------
 SUBROUTINE EMWaveCalcFields(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
   USE EMWaveSolverUtils

   IMPLICIT NONE
!------------------------------------------------------------------------------
   TYPE(Solver_t) :: Solver
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: dt
   LOGICAL :: Transient
!------------------------------------------------------------------------------
   TYPE(Variable_t), POINTER :: pVar
   TYPE(Variable_t), POINTER :: EF, dEF, ddEF
   TYPE(Variable_t), POINTER :: EF_e, dEF_e, ddEF_e
                              
   INTEGER :: i,j,k,l,t,n,nd,p,q,dofs,dofcount,vDOFs

   TYPE(Solver_t), POINTER :: pSolver
   CHARACTER(LEN=MAX_NAME_LEN) :: Pname

   LOGICAL :: Found, stat

   TYPE(Element_t), POINTER :: Element

   INTEGER, ALLOCATABLE :: Pivot(:)

   REAL(KIND=dp), POINTER CONTIG :: Fsave(:)
   TYPE(Mesh_t), POINTER :: Mesh
   REAL(KIND=dp), ALLOCATABLE, TARGET :: MASS(:,:), FORCE(:,:), GForce(:,:) 
   LOGICAL :: PiolaVersion, ElementalFields, NodalFields, SecondOrder, AnyTimeDer, &
       ConstantBulkInUse = .FALSE.
   INTEGER :: soln
   TYPE(ValueList_t), POINTER :: SolverParams 

   REAL(KIND=dp) :: mu, eps, mu0, eps0

   SAVE ConstantBulkInUse
   
!-------------------------------------------------------------------------------------------
   CALL Info('EMWaveCalcFields','Computing postprocessing fields')

   SolverParams => GetSolverParams()

   eps0 = GetConstReal( Model % Constants,'Permittivity of Vacuum')
   mu0 = GetConstReal( Model % Constants,'Permeability of Vacuum')

   soln = ListGetInteger( SolverParams,'Primary Solver Index', Found) 
   IF( soln == 0 ) THEN
     CALL Fatal('EMWaveCalcFields','We should know > Primary Solver Index <')
   END IF

   ! Pointer to primary solver
   pSolver => Model % Solvers(soln)
   pVar => pSolver % Variable   
   Pname = getVarName(pVar)
   CALL Info('EMWaveCalcFields','Name of potential variable: '//TRIM(pName),Level=10)
   
   ! Inherit the solution basis from the primary solver
   vDOFs = pVar % DOFs
   IF( vDofs /= 1 ) THEN
     CALL Fatal('EMWaveCalcFields','Primary variable should have 1 dofs: '//TRIM(I2S(vDofs)))
   END IF
   dofs = 3
   
   SecondOrder = GetLogical( pSolver % Values, 'Quadratic Approximation', Found )  
   IF( SecondOrder ) THEN
     PiolaVersion = .TRUE.
   ELSE
     PiolaVersion = GetLogical( pSolver % Values,'Use Piola Transform', Found ) 
   END IF
   
   IF (PiolaVersion) &
       CALL Info('MagnetoDynamicsCalcFields', &
       'Using Piola transformed finite elements',Level=5)

   Mesh => GetMesh()

   EF => VariableGet( Mesh % Variables, 'Elfield')
   EF_e => VariableGet( Mesh % Variables, 'Elfield E')
   
   dEF => VariableGet( Mesh % Variables, 'dEdt')
   dEF_e => VariableGet( Mesh % Variables, 'dEdt E')
   
   ddEF => VariableGet( Mesh % Variables, 'ddEddt')
   ddEF_e => VariableGet( Mesh % Variables, 'ddEddt E')
   
   AnyTimeDer = ASSOCIATED( dEF ) .OR. ASSOCIATED( ddEF ) .OR. &       
       ASSOCIATED( dEF_e ) .OR. ASSOCIATED( ddEF_e ) 
          
   i = 0 
   IF ( ASSOCIATED(EF)  ) i=i+3
   IF ( ASSOCIATED(dEF)  ) i=i+3
   IF ( ASSOCIATED(ddEF)  ) i=i+3
   NodalFields = ( i > 0 )

   IF(NodalFields) THEN
     ALLOCATE(GForce(SIZE(Solver % Matrix % RHS),i)); Gforce=0._dp
   END IF
      
   j = 0 
   IF ( ASSOCIATED(EF_e)  ) j=j+3
   IF ( ASSOCIATED(dEF_e)  ) j=j+3
   IF ( ASSOCIATED(ddEF_e)  ) j=j+3
   ElementalFields = ( j > 0 )
   
   dofs = MAX( i,j )
   
   n = Mesh % MaxElementDOFs   
   ALLOCATE( MASS(n,n), FORCE(n,dofs), Pivot(n) )
   
   CALL DefaultInitialize(UseConstantBulk = ConstantBulkInUse )
   
   DO t = 1, GetNOFActive()
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     nd = GetElementNOFDOFs(uSolver=pSolver)
     
     CALL LocalAssembly(t==1)
     
     IF(NodalFields) THEN
       IF(.NOT. ConstantBulkInUse ) THEN
         CALL DefaultUpdateEquations( MASS,FORCE(:,1))
       END IF
       Fsave => Solver % Matrix % RHS
       DO l=1,dofs
         Solver % Matrix % RHS => GForce(:,l)
         CALL DefaultUpdateForce(FORCE(:,l))
       END DO
       Solver % Matrix % RHS => Fsave
     END IF

     IF(ElementalFields) THEN
       dofcount = 0
       CALL LUdecomp(MASS,n,pivot)
       CALL LocalSol(EF_e,   3, n, MASS, FORCE, pivot, dofcount)
       CALL LocalSol(dEF_e,  3, n, MASS, FORCE, pivot, dofcount)
       CALL LocalSol(ddEF_e, 3, n, MASS, FORCE, pivot, dofcount)
     END IF

   END DO

   ! Assembly of the face terms in case we have DG method where we
   ! want to average the fields within materials making them continuous.
   !-----------------------------------------------------------------
   IF(.NOT. ConstantBulkInUse ) THEN
     IF (GetLogical( SolverParams,'Discontinuous Galerkin',Found)) THEN
       IF (GetLogical( SolverParams,'Average Within Materials',Found)) THEN
         FORCE = 0.0_dp
         CALL AddLocalFaceTerms( MASS, FORCE(:,1) )
       END IF
     END IF
   END IF

   IF(NodalFields) THEN
     Fsave => Solver % Matrix % RHS
     dofcount = 0
     CALL GlobalSol(EF ,  3, Gforce, dofcount)
     CALL GlobalSol(dEF , 3, Gforce, dofcount)
     CALL GlobalSol(ddEF ,3, Gforce, dofcount)
     Solver % Matrix % RHS => Fsave
   END IF

   ConstantBulkInUse = ListGetLogical( SolverParams,'Constant Bulk Matrix',Found )


CONTAINS


  SUBROUTINE LocalAssembly(InitHandles) 

    LOGICAL :: InitHandles

    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:), SOL(:), dsol(:), ddsol(:)
    REAL(KIND=dp), ALLOCATABLE :: RotWBasis(:,:), Basis(:), dBasisdx(:,:)
    INTEGER, ALLOCATABLE :: Indexes(:)
    REAL(KIND=dp) :: detJ, s, u, v, w
    REAL(KIND=dp) :: EF_ip(3), dEF_ip(3), ddEF_ip(3)
    LOGICAL :: AllocationsDone = .FALSE.
    TYPE(ValueHandle_t) :: Mu_h, Eps_h, Cd_h(3)    
    INTEGER :: i,j,k,n
    
    SAVE Cd_h, Mu_h, Eps_h, sol, dsol, ddsol, &
        Indexes, WBasis, RotWBasis, Basis, dBasisdx

    IF(.NOT. AllocationsDone ) THEN
      N = Mesh % MaxElementDOFs
      ALLOCATE( SOL(n), dsol(n), ddsol(n), &
          Indexes(n), WBasis(n,3), RotWBasis(n,3), Basis(n), dBasisdx(n,3))      
      sol = 0.0_dp
      dsol = 0.0_dp
      ddsol = 0.0_dp
      AllocationsDone = .TRUE.
    END IF

    
    IF( InitHandles ) THEN
      ! These have been normalized by mu0 and eps0 in _init section
      CALL ListInitElementKeyword( Mu_h,'Material','Relative Permeability')
      CALL ListInitElementKeyword( Eps_h,'Material','Relative Permittivity')

      CALL ListInitElementKeyword( Cd_h(1),'Body Force','Current Density 1')
      CALL ListInitElementKeyword( Cd_h(2),'Body Force','Current Density 2')
      CALL ListInitElementKeyword( Cd_h(3),'Body Force','Current Density 3')
    END IF

    
    CALL GetElementNodes( Nodes )
    
    n = GetElementDOFs( Indexes, Element, pSolver )
    
    DO i=1,n
      j = pVar % Perm( Indexes(i) )       
      IF ( j > 0 ) THEN
        sol(i) = pVar % Values(j)
        ! These are only applicable as we know this is 2nd order PDE
        IF( AnyTimeDer ) THEN
          dsol(i) = pVar % PrevValues(j,1) 
          ddsol(i) = pVar % PrevValues(j,2)         
        END IF
      END IF
    END DO
    
    ! Calculate nodal fields:
    ! -----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion)

    MASS  = 0._dp
    FORCE = 0._dp

    DO j = 1,IP % n
      u = IP % U(j)
      v = IP % V(j)
      w = IP % W(j)

      stat = ElementInfo(Element,Nodes,u,v,w,detJ,Basis,dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 
 
      ! Not currently used as only trivial fields are computed. 
      !----------------------------------------------------------
      !eps = eps0 * ListGetElementReal( eps_h, Basis, Element, Found )
      !mu = mu0 * ListGetElementReal( mu_h, Basis, Element, Found )

      !curr(1) = ListGetElementReal( Cd_h(1), Basis, Element, Found )
      !curr(2) = ListGetElementReal( Cd_h(2), Basis, Element, Found )
      !curr(3) = ListGetElementReal( Cd_h(3), Basis, Element, Found )
      
      EF_ip = MATMUL(sol(1:nd),WBasis(1:nd,:))
      IF( AnyTimeDer ) THEN
        dEF_ip = MATMUL(dsol(1:nd),WBasis(1:nd,:))
        ddEF_ip = MATMUL(ddsol(1:nd),WBasis(1:nd,:))
      END IF
        
      s = IP % s(j) * detJ

      DO p=1,n
        DO q=1,n
          MASS(p,q)=MASS(p,q)+s*Basis(p)*Basis(q)
        END DO
        k = 0
        IF ( ASSOCIATED(EF).OR.ASSOCIATED(EF_e)) THEN
          FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*(EF_ip)*Basis(p)
          k = k+3
        END IF
        IF ( ASSOCIATED(dEF).OR.ASSOCIATED(dEF_e)) THEN
          FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*(dEF_ip)*Basis(p)
          k = k+3
        END IF
        IF ( ASSOCIATED(ddEF).OR.ASSOCIATED(ddEF_e)) THEN
          FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*(ddEF_ip)*Basis(p)
          k = k+3
        END IF
      END DO
    END DO
    
  END SUBROUTINE LocalAssembly
    
  
!------------------------------------------------------------------------------
 SUBROUTINE GlobalSol(Var, m, b, dofs )
!------------------------------------------------------------------------------
   REAL(KIND=dp), TARGET CONTIG :: b(:,:)
   INTEGER :: m, dofs
   TYPE(Variable_t), POINTER :: Var
!------------------------------------------------------------------------------
   INTEGER :: i
   REAL(KIND=dp) :: Norm
!------------------------------------------------------------------------------
   IF(.NOT. ASSOCIATED(var)) RETURN

   DO i=1,m
     dofs = dofs+1
     Solver % Matrix % RHS => b(:,dofs)
     Solver % Variable % Values=0
     Norm = DefaultSolve()
     var % Values(i::m) = Solver % Variable % Values
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE GlobalSol
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE LocalSol(Var, m, n, A, b, pivot, dofs )
!------------------------------------------------------------------------------
   TYPE(Variable_t), POINTER :: Var
   INTEGER :: pivot(:), m,n,dofs
   REAL(KIND=dp) :: b(:,:), A(:,:)
!------------------------------------------------------------------------------
   INTEGER :: ind(n), i
   REAL(KIND=dp) :: x(n)
!------------------------------------------------------------------------------
   IF(.NOT. ASSOCIATED(var)) RETURN

   ind = Var % dofs*(Var % Perm(Element % DGIndexes(1:n))-1)
   DO i=1,m
      dofs = dofs+1
      x = b(1:n,dofs)
      CALL LUSolve(n,MASS,x,pivot)
      Var % Values(ind(1:n)+i) = x
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE LocalSol
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE AddLocalFaceTerms(STIFF,FORCE)
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: STIFF(:,:), FORCE(:)

   TYPE(Element_t),POINTER :: P1,P2,Face,Faces(:)
   INTEGER ::t,n,n1,n2,NumberOfFaces,dim

   dim = CoordinateSystemDimension()

   IF (dim==2) THEN
     Faces => Solver % Mesh % Edges
     NumberOfFaces = Solver % Mesh % NumberOfEdges
   ELSE
     Faces => Solver % Mesh % Faces
     NumberOfFaces = Solver % Mesh % NumberOfFaces
   END IF

   DO t=1,NumberOfFaces
     Face => Faces(t)
     IF ( .NOT. ActiveBoundaryElement(Face) ) CYCLE

     P1 => Face % BoundaryInfo % Left
     P2 => Face % BoundaryInfo % Right
     IF ( ASSOCIATED(P2) .AND. ASSOCIATED(P1) ) THEN
       IF(.NOT.ASSOCIATED(GetMaterial(P1),GetMaterial(P2))) CYCLE

       n  = GetElementNOFNodes(Face)
       n1 = GetElementNOFNodes(P1)
       n2 = GetElementNOFNodes(P2)

       CALL LocalJumps( STIFF,Face,n,P1,n1,P2,n2)
       CALL DefaultUpdateEquations( STIFF, FORCE, Face )
     END IF
   END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AddLocalFaceTerms
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalJumps( STIFF,Face,n,P1,n1,P2,n2)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: STIFF(:,:)
    INTEGER :: n,n1,n2
    TYPE(Element_t), POINTER :: Face, P1, P2
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: FaceBasis(n), P1Basis(n1), P2Basis(n2)
    REAL(KIND=dp) :: Jump(n1+n2), detJ, U, V, W, S
    LOGICAL :: Stat
    INTEGER ::  p, q, t
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

    TYPE(Nodes_t) :: FaceNodes, P1Nodes, P2Nodes
    SAVE FaceNodes, P1Nodes, P2Nodes
!------------------------------------------------------------------------------
    STIFF = 0._dp

    CALL GetElementNodes(FaceNodes, Face)
    CALL GetElementNodes(P1Nodes, P1)
    CALL GetElementNodes(P2Nodes, P2)
!------------------------------------------------------------------------------
!     Numerical integration over the edge
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Face )

    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)

      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo(Face, FaceNodes, U, V, W, detJ, FaceBasis)

      S = S * detJ

      ! Find basis functions for the parent elements:
      ! ---------------------------------------------
      CALL GetParentUVW(Face, n, P1, n1, U, V, W, FaceBasis)
      stat = ElementInfo(P1, P1Nodes, U, V, W, detJ, P1Basis)

      CALL GetParentUVW(Face, n, P2, n2, U, V, W, FaceBasis)
      stat = ElementInfo(P2, P2Nodes, U, V, W, detJ, P2Basis)

      ! Integrate jump terms:
      ! ---------------------
      Jump(1:n1) = P1Basis(1:n1)
      Jump(n1+1:n1+n2) = -P2Basis(1:n2)

      DO p=1,n1+n2
        DO q=1,n1+n2
          STIFF(p,q) = STIFF(p,q) + s * Jump(q)*Jump(p)
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------

!------------------------------------------------------------------------
END SUBROUTINE EMWaveCalcFields
!------------------------------------------------------------------------

