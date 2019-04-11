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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> The current density given as a source must be divergence free to allow a 
!> hope for a solution. This solver may be used to enforce a given current 
!> density to be divergence free by solving for an equivalent potential
!> so that the vector field is a gradient of the potential.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE JfixPotentialSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL ::  Transient
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams, BF
  INTEGER :: i,j,k,n,m,dim,dofs
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Matrix_t), POINTER :: A,B
  REAL(KIND=dp) :: Norm
  LOGICAL:: AutomatedBCs, SingleNodeBC, Found
  INTEGER, ALLOCATABLE :: Def_Dofs(:,:,:)
  CHARACTER(LEN=MAX_NAME_LEN):: Equation
  INTEGER, POINTER :: Perm(:)
  REAL(KIND=dp), POINTER :: fixpot(:)
  TYPE(Variable_t), POINTER :: fixJpot, svar, IterV 
  LOGICAL :: StatCurrMode, NeumannMode, DirichletMode, Visited = .FALSE.

  CALL Info('JfixPotentialSolver','Computing fixing potential for given current density',Level=6)
  
  dim = CoordinateSystemDimension()
  Mesh => GetMesh()
  B => GetMatrix()
  SolverParams => GetSolverParams()

  IF( ListGetLogical( SolverParams,'Constant Input Current Density',Found ) ) THEN
    IF( Visited ) THEN
      CALL Info('JfixPotentialSolver','Current density is constant, nothing to do!',Level=8)
      RETURN
    END IF
  END IF  
  
  fixJpot => VariableGet( Mesh % Variables, 'Jfix')

  IF( ASSOCIATED(fixJPot)) THEN
    Perm => fixJPot % Perm
  ELSE
    ALLOCATE(Perm(SIZE(Solver % Variable % Perm)))
  END IF
  Perm = 0
  dofs = Solver % Variable % DOFs

  Equation=GetString(SolverParams,'Equation',Found)

  n=SIZE(Solver % Def_Dofs,1)
  m=SIZE(Solver % Def_Dofs,2)
  k=SIZE(Solver % Def_Dofs,3)
  ALLOCATE(Def_Dofs(n,m,k))
  Def_Dofs=Solver % Def_Dofs
  Solver % Def_Dofs=0
  Solver % Def_Dofs(:,:,1)=1

  A => CreateMatrix( CurrentModel, Solver, Solver % Mesh, &
     Perm, dofs, MATRIX_CRS, .TRUE., Equation, .FALSE., .FALSE.,&
     NodalDofsOnly = .TRUE.)
  n = A % NumberOfRows
  IF (dofs>1) A % COMPLEX = .TRUE.

  IF (.NOT.ASSOCIATED(fixJPot)) THEN
    ALLOCATE(fixpot(n)); fixpot=0._dp

    CALL VariableAddVector( Mesh % Variables, Mesh, &
          Solver,'Jfix',dofs,fixpot,Perm)

    fixJpot => VariableGet(Mesh % Variables, 'Jfix')
  END IF
  
  
  ! Set potential only on a single node
  ! This uses the functionality of DefaultDirichlet
  SingleNodeBC = GetLogical( SolverParams,'Single Node Projection BC',Found ) 
  IF( SingleNodeBC .AND. .NOT. Visited ) THEN
    DO i=1,Model % NumberOfBodyForces
      BF => Model % BodyForces(i) % Values
      IF( ListCheckPrefix( BF,'Current Density') ) THEN
        CALL ListAddConstReal( BF,'Jfix Single Node',0.0_dp )
      END IF
    END DO    
  END IF
  ! We may have used the single node BC directly (in case of many body forces)
  IF(.NOT. SingleNodeBC ) THEN
    SingleNodeBC = ListCheckPresentAnyBodyForce( Model,'Jfix Single Node' ) 
  END IF
  

  StatCurrMode = GetLogical( SolverParams,'StatCurrent Source Projection', Found )
  DirichletMode = GetLogical( SolverParams,'Dirichlet Source Projection',Found )
  NeumannMode = GetLogical( SolverParams,'Neumann Source Projection',Found )

  ! If not any mode given revert to the old standard
  IF(.NOT.( StatCurrMode .OR. DirichletMode .OR. NeumannMode .OR. SingleNodeBC ) ) THEN  
    DirichletMode = GetLogical( SolverParams, &
        'Automated Source Projection BCs', Found )
    IF(.NOT. Found) DirichletMode = .TRUE.
    StatCurrMode = .NOT. DirichletMode
  END IF

  IF( A % COMPLEX ) THEN
    IF( StatCurrMode .OR. NeumannMode ) THEN
      CALL Warn('JfixPotentialSolver','Cannot set Neumann conditions for complex fields!')
      StatCurrMode = .FALSE.; NeumannMode = .FALSE.
    END IF
  END IF
    
  svar => Solver % Variable
  Solver % Variable => fixJpot

  Solver % Matrix => A
  ALLOCATE(A % RHS(n))
  IF(ParEnv % PEs>1) CALL ParallelInitMatrix(Solver,A)

  CALL ListSetNameSpace('jfix:')
  CALL ListAddNewString(SolverParams,'Jfix: Linear System Solver', 'Iterative')
  CALL ListAddNewString(SolverParams,'Jfix: Linear System Iterative Method', 'BiCGStab')
  CALL ListAddNewLogical(SolverParams,'Jfix: Linear System Use HYPRE', .FALSE.)
  CALL ListAddNewLogical(SolverParams,'Jfix: Use Global Mass Matrix',.FALSE.)
  CALL ListAddNewString(SolverParams,'Jfix: Linear System Preconditioning', 'Ilu')
  CALL ListAddNewConstReal(SolverParams,'Jfix: Linear System Convergence Tolerance', &
      0.001_dp*GetCReal(SolverParams,'Linear System Convergence Tolerance', Found))
  CALL ListAddNewLogical(SolverParams,'Jfix: Skip Compute Nonlinear Change',.TRUE.)
  CALL ListAddNewLogical(SolverParams,'Jfix: Nonlinear System Consistent Norm',.TRUE.)
  CALL ListAddNewString(SolverParams,'Jfix: Nonlinear System Convergence Measure','Norm')

  CALL DefaultInitialize()
  CALL BulkAssembly()
  
  IF( NeumannMode .OR. StatCurrMode ) THEN
    CALL BCAssemblyNeumann()
  ELSE IF( DirichletMode ) THEN
    CALL BCAssemblyDirichlet()
  END IF

  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  n = 0
  IF( ALLOCATED( A % ConstrainedDOF ) ) THEN
    n = COUNT( A % ConstrainedDOF )
    n = NINT( ParallelReduction( 1.0_dp * n ) )
  END IF
  IF( n == 0 ) THEN
    CALL Warn('JfixPotentialSolver','No Dirichlet conditions used to define Jfix level!')
  ELSE
    CALL Info('JfixPotentialSolver','Number of dirichlet nodes: '//TRIM(I2S(n)),Level=7)
  END IF
  
  Norm = DefaultSolve()

  Solver % Matrix => B
  Solver % Variable => svar
  CALL ListSetNameSpace('')

  CALL FreeMatrix(A)
  Solver % Def_Dofs=Def_Dofs
  
  IterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
  IterV % Values(1) = 1

  Visited = .TRUE.

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
    IMPLICIT NONE       
    COMPLEX(KIND=dp), ALLOCATABLE :: STIFF_C(:,:), FORCE_C(:)
    REAL(KIND=dp), ALLOCATABLE :: STIFF_R(:,:), FORCE_R(:)
    INTEGER :: elem,t,i,j,k,p,q,n,nd
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER ::  BodyForce
    REAL(KIND=dp) :: weight,detJ,L_R(3),Nrm(3),rabs
    REAL(KIND=dp), ALLOCATABLE :: Load_R(:,:)
    COMPLEX(KIND=dp) :: L_C(3)
    COMPLEX(KIND=dp), ALLOCATABLE :: Load_C(:,:)
    INTEGER :: l,m
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    
    SAVE Nodes
    
    n = MAX(Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes)
    IF (A % COMPLEX) THEN
      ALLOCATE( Load_c(3,n), STIFF_C(n,n), FORCE_C(n) )
    ELSE
      ALLOCATE( Load_r(3,n), STIFF_R(n,n), FORCE_R(n) )
    END IF
    ALLOCATE( Basis(n), dBasisdx(n,3) )

    DO elem = 1,GetNOFActive()
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE

      CALL GetElementNodes( Nodes )
      n  = GetElementNOFNodes()
      nd = n

      IF (A % COMPLEX ) THEN
        Load_C =0._dp
        STIFF_C=0._dp
        FORCE_C=0._dp
      ELSE
        Load_R =0._dp
        STIFF_R=0._dp
        FORCE_R=0._dp
      END IF

      BodyForce => GetBodyForce()
      IF (ASSOCIATED(BodyForce)) THEN
        IF (A % COMPLEX ) THEN
          CALL GetComplexVector(BodyForce,Load_C(1:3,1:n),'Current Density',Found)
        ELSE
          CALL GetRealVector(BodyForce,Load_R(1:3,1:n),'Current Density',Found)
        END IF
      END IF

      IntegStuff = GaussPoints( Element )
      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
                IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )
        
        Weight = IntegStuff % s(t) * detJ
        IF (A % COMPLEX) THEN
          DO p=1,n
            DO q=1,n
              STIFF_C(p,q) = STIFF_C(p,q) + Weight * SUM(dBasisdx(q,:)*dBasisdx(p,:))
            END DO
          END DO
          L_C = MATMUL( Load_C(:,1:n), Basis(1:n) )
          FORCE_C(1:n) = FORCE_C(1:n) + MATMUL(dBasisdx(1:n,:),L_C)*Weight
        ELSE
          DO p=1,n
            DO q=1,n
              STIFF_R(p,q) = STIFF_R(p,q) + Weight * SUM(dBasisdx(q,:)*dBasisdx(p,:))
            END DO
          END DO
          L_R = MATMUL( Load_R(:,1:n), Basis(1:n) )
          FORCE_R(1:n) = FORCE_R(1:n) + MATMUL(dBasisdx(1:n,:),L_R)*Weight
        END IF
      END DO

      IF (A % COMPLEX) THEN
        DO i=n+1,nd
          STIFF_C(i,i)=1._dp
        END DO
        CALL DefaultUpdateEquations(STIFF_C, FORCE_C)
      ELSE
        DO i=n+1,nd
          STIFF_R(i,i)=1._dp
        END DO
        CALL DefaultUpdateEquations(STIFF_R, FORCE_R)
      END IF
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------


  FUNCTION NormalCurrentDensity(Element, Nodes, NormalProj, GotCond) RESULT( CurrDensNormal ) 
    
    IMPLICIT NONE       
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element, LElement
    REAL(KIND=dp) :: NormalProj
    LOGICAL :: GotCond 
    COMPLEX(KIND=dp) :: CurrDensNormal
    
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: P1,P2, P
    TYPE(ValueList_t), POINTER ::  BodyForce
    REAL(KIND=dp) :: Nrm(3),s
    REAL(KIND=dp), ALLOCATABLE :: Load_RP(:,:)
    COMPLEX(KIND=dp), ALLOCATABLE :: Load_CP(:,:)
    INTEGER :: ActParents,ParParents,n,np,k,j
    LOGICAL :: AllocationsDone = .FALSE.
    COMPLEX(KIND=dp) :: CurrDensVec(3)

    SAVE AllocationsDone, Load_RP, Load_CP
    
    GotCond = .FALSE.
    CurrDensNormal = 0.0_dp
    
    IF( .NOT. AllocationsDone ) THEN      
      n = Solver % Mesh % MaxElementDOFs
      IF (A % COMPLEX) THEN
        ALLOCATE( Load_cp(3,n) )
      ELSE
        ALLOCATE( Load_rp(3,n) )
      END IF
      AllocationsDone = .TRUE.
    END IF

   SELECT CASE(GetElementFamily())
    CASE(1)
    CASE(2)
      k = GetBoundaryEdgeIndex(Element,1)
      IF( ParEnv % PEs > 1 ) THEN
        IF( Mesh % ParallelInfo % EdgeInterface(k) ) RETURN        
      END IF
      LElement => Mesh % Edges(k)
    CASE(3,4)
      k = GetBoundaryFaceIndex(Element)
      IF( ParEnv % PEs > 1 ) THEN
        ! Don't set BCs on partition interfaces
        IF( Mesh % ParallelInfo % FaceInterface(k) ) RETURN
      END IF
      LElement => Mesh % Faces(k)
    END SELECT
        
    P1 => LElement % BoundaryInfo % Left
    P2 => LElement % BoundaryInfo % Right
    
    ActParents = 0
    ParParents = 0
    
    IF( ASSOCIATED( P1 ) ) THEN
      IF (ALL(Perm(P1 % NodeIndexes)>0)) THEN
        ActParents = ActParents + 1
        IF( P1 % PartIndex == ParEnv % MyPe ) ParParents = ParParents + 1
      ELSE
        NULLIFY( P1 )
      END IF
    END IF
    IF( ASSOCIATED( P2 ) ) THEN
      IF (ALL(Perm(P2 % NodeIndexes)>0)) THEN
        ActParents = ActParents + 1
        IF( P2 % PartIndex == ParEnv % MyPe ) ParParents = ParParents + 1         
      ELSE
        NULLIFY( P2 )
      END IF
    END IF

    ! We have either none or both parents as actice.
    ! The BCs will be set only to outer boundaries of the domain. 
    IF( ActParents /= 1 ) RETURN 

    ! The one parent is not a true one!
    IF( ParEnv % PEs > 0 .AND. ParParents == 0 ) RETURN

    IF (ASSOCIATED(P1)) THEN
      P => P1
    ELSE
      P => P2
    END IF

    n = Element % Type % NumberOfNodes 
    np = P % TYPE % NumberOfNodes

    BodyForce => GetBodyForce(P)

    IF (.NOT. ASSOCIATED(BodyForce)) RETURN
    
    ! Normal vector that points out of Parent 
    Nrm = NormalVector(LElement,Nodes,Parent=P)

    CurrentModel % CurrentElement => P
    
    IF (A % COMPLEX ) THEN
      CALL GetComplexVector( BodyForce, Load_CP(1:3,1:np),'Current Density', Found )
    ELSE      
      Load_RP(1,1:np) = GetReal( BodyForce, 'Current Density 1', Found, UElement=P )
      Load_RP(2,1:np) = GetReal( BodyForce, 'Current Density 2', Found, UElement=P )
      Load_RP(3,1:np) = GetReal( BodyForce, 'Current Density 3', Found, UElement=P )
    END IF

    CurrDensVec = 0.0_dp
    DO k=1,n
      DO j=1,np
        IF(Element % NodeIndexes(k) == P % NodeIndexes(j)) THEN
          IF( A % COMPLEX ) THEN
            CurrDensVec = CurrDensVec + Load_CP(1:3,j)
          ELSE
            CurrDensVec = CurrDensVec + Load_RP(1:3,j)
          END IF
        END IF
      END DO
    END DO
    CurrDensVec = CurrDensVec / n

    ! This sign convention is consistent
    CurrDensNormal = SUM( -Nrm * CurrDensVec ) 

    s = SQRT( SUM( CurrDensVec**2 ) ) 
    IF( s > TINY( s ) ) THEN
      NormalProj = ABS( CurrDensNormal / s )
      GotCond = .TRUE.
    ELSE
      NormalProj = 0.0_dp
      GotCond = .FALSE.
    END IF
        
    CurrentModel % CurrentElement => Element

    
  END FUNCTION NormalCurrentDensity


  
!------------------------------------------------------------------------------
! This subroutine fixes the potential to zero where there is some current density
! component in the normal direction.
!------------------------------------------------------------------------------
  SUBROUTINE BCAssemblyDirichlet()
!------------------------------------------------------------------------------
    IMPLICIT NONE       
    INTEGER :: i,n,meshdim,ip,dofs
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: Jfluxeps,NormalProj
    COMPLEX(KIND=dp) :: CurrDens
    
    SAVE Nodes

    meshdim = Solver % Mesh % MeshDim    
    Jfluxeps = GetCReal(SolverParams, 'J normal eps', Found)
    IF (.NOT. Found) Jfluxeps = 0.1_dp

    ! At outer boundaries, if J has component normal to the
    ! surface, set Fix_pot=0:
    ! THIS MAY NOT DO WHAT IS PHYSICALLY SOUND: It may be better
    ! to disable this by the command 
    ! Automated Projection BCs = Logical False
    ! and set BCs explicitly in the sif file    
    ! ------------------------------------------------------

    dofs = Solver % Variable % Dofs

    DO i=1,GetNOFBoundaryElements()
      Element => GetBoundaryElement(i)
      n = GetElementNOFNodes()
      IF (.NOT.ActiveBoundaryElement()) CYCLE

      IF (Element % TYPE % DIMENSION < meshdim-1 ) CYCLE
      
      CALL GetElementNodes(Nodes)
      
      CurrDens = NormalCurrentDensity(Element,Nodes,NormalProj,Found)

      IF(.NOT. Found ) CYCLE
      
      ! If there is no normal component of the current density then
      ! don't set the corresponding Dirichlet condition to zero. 
      IF( NormalProj < Jfluxeps ) CYCLE

      DO j=1,dofs
        DO k=1,n
          ip = Element % NodeIndexes(k) 
          ip = dofs*(Perm(ip)-1)+j
          CALL UpdateDirichletDof(A, ip, 0._dp)
        END DO
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE BCAssemblyDirichlet
!------------------------------------------------------------------------------

  
  
!------------------------------------------------------------------------------
  SUBROUTINE BCAssemblyNeumann()
!------------------------------------------------------------------------------
!   This is intended to alter the natural boundary conditions of the potential associated
!   with the Helmholtz projection if the source electric current density is
!   generated by applying the StatCurrentSolve module.
!-------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), ALLOCATABLE :: FORCE_R(:), LOAD(:), Basis(:), dBasisdx(:,:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IP  
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: BC
    INTEGER :: Active, i, k, p, t, nd, n
    REAL(KIND=dp) :: detJ, Jn, NormalProj
    LOGICAL :: Found, Stat
    COMPLEX(KIND=dp) :: CurrDens    
    TYPE(Nodes_t), SAVE :: Nodes
!-------------------------------------------------------------------------------
 
    n = MAX(Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes)  
    ALLOCATE( FORCE_R(n), LOAD(n), Basis(n), dBasisdx(n,3) )    

    Active = GetNOFBoundaryElements()
    DO k=1,Active
      Element => GetBoundaryElement(k)
      IF (.NOT. ActiveBoundaryElement()) CYCLE

      IF (GetElementFamily()==1) CYCLE
      
      nd = GetElementNOFDOFs(Element)
      n  = GetElementNOFNodes(Element)

      Load = 0.0d0
      FORCE_R = 0.0d0

      CALL GetElementNodes( Nodes )

      IF( StatCurrMode ) THEN
        BC => GetBC()
        IF (.NOT. ASSOCIATED(BC) ) CYCLE
        Load(1:n) = GetReal( BC, 'Current Density', Found )
      ELSE
        CurrDens = NormalCurrentDensity(Element,Nodes,NormalProj,Found)
        Jn = REAL( CurrDens )
      END IF

      IF(.NOT. Found ) CYCLE       
        
      IF(.NOT. StatCurrMode ) THEN
        !PRINT *,'neu2:',NormalProj,Jn,n,nd      

        ! This is a mixed mode that sets only positive Jn values to Dirichlet
        IF( DirichletMode ) THEN
          IF( ABS( NormalProj ) > 0.5 .AND. Jn > 0.0_dp ) THEN
            DO i=1,n
              j = Perm(Element % NodeIndexes(i)) 
              CALL UpdateDirichletDof(A, j, 0._dp)
            END DO
            RETURN
          END IF
        END IF        
      END IF
        
             
      IP = GaussPoints( Element )
      DO t=1,IP % n
        stat = ElementInfo( Element, Nodes, IP % u(t), &
            IP % v(t), IP % w(t), detJ, Basis, dBasisdx )

        IF( StatCurrMode ) THEN
          Jn = -SUM( Load(1:n)*Basis(1:n) ) 
        END IF

        DO p=1,n
          FORCE_R(p) = FORCE_R(p) + Jn*Basis(p)*detJ*IP % s(t)
        END DO
      END DO

      CALL DefaultUpdateForce(FORCE_R,Element)

    END DO

    DEALLOCATE(FORCE_R,LOAD,Basis,dBasisdx)
!------------------------------------------------------------------------------
  END SUBROUTINE BCAssemblyNeumann
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
END SUBROUTINE JfixPotentialSolver
!------------------------------------------------------------------------------
