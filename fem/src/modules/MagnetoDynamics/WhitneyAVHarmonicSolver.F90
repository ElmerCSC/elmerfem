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
 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE WhitneyAVHarmonicSolver_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: Found, PiolaVersion, SecondOrder, SecondFamily
  CHARACTER(:), ALLOCATABLE :: ElemType

  TYPE(Solver_t), POINTER :: Solvers(:)
  INTEGER :: i,j,k,n
  CHARACTER(:), ALLOCATABLE :: eq
  INTEGER, POINTER :: ActiveSolvers(:)
  
  SolverParams => GetSolverParams()
  IF ( .NOT.ListCheckPresent(SolverParams, "Element") ) THEN
    ! We use one place where all the edge element keywords are defined and checked.
    CALL EdgeElementStyle(SolverParams, PiolaVersion, SecondFamily, SecondOrder, Check = .TRUE. )

    IF (SecondOrder) THEN
      ElemType = "n:1 e:2 -brick b:6 -prism b:2 -pyramid b:3 -quad_face b:4 -tri_face b:2" 
    ELSE IF( SecondFamily ) THEN
      ElemType = "n:1 e:2"
    ELSE IF( PiolaVersion ) THEN
      ElemType = "n:1 e:1 -brick b:3 -quad_face b:2" 
    ELSE
      ElemType = "n:1 e:1"
    END IF
    CALL Info('WhitneyHarmonicSolver_Init0','Setting element type to: "'//TRIM(ElemType)//'"',Level=6)
    CALL ListAddString( SolverParams, "Element", TRIM(ElemType) )
  END IF

  CALL ListAddNewLogical( SolverParams, 'Linear System Complex', .TRUE. )

! This is for internal communication with the saving routines
  CALL ListAddLogical( SolverParams,'Hcurl Basis',.TRUE.)

  CALL ListAddNewString( SolverParams,'Variable','AV[AV re:1 AV im:1]')

  IF(ListGetLogical(SolverParams, 'Helmholtz Projection', Found)) THEN
    Solvers => Model % Solvers
    n = Model % NumberOfSolvers
    Model % NumberOfSolvers = n+2

    ALLOCATE(Model % Solvers(Model % NumberOfSolvers))
    Model % Solvers(1:n) = Solvers

    DO i=n+1,n+2
      Model % Solvers(i) % PROCEDURE = 0
      NULLIFY( Model % Solvers(i) % Matrix )
      NULLIFY( Model % Solvers(i) % Mesh )
      NULLIFY( Model % Solvers(i) % Variable )
      NULLIFY( Model % Solvers(i) % ActiveElements )
      Model % Solvers(i) % NumberOfActiveElements = 0
    END DO

    DO i=1,Model % NumberOfSolvers
      IF(.NOT.ASSOCIATED(Model % Solvers(i) % Values)) &
         Model % Solvers(i) % Values => ListAllocate()
    END DO

    Eq =  ListGetString( SolverParams, 'Equation' )

    CALL ListAddIntegerArray( SolverParams, 'Post Solvers', 2, [n+1,n+2] )

    CALL ListAddString( Model % Solvers(n+1) % Values, 'Procedure', &
              'MagnetoDynamics HelmholtzProjector', CaseConversion=.FALSE. )
    CALL ListAddString( Model % Solvers(n+1) % Values, 'Equation', 'HP' )
    CALL ListAddString( Model % Solvers(n+1) % Values, 'Exec Solver', 'Never' )

    CALL ListAddString( Model % Solvers(n+2) % Values, 'Procedure', &
              'MagnetoDynamics RemoveKernelComponent',CaseConversion=.FALSE. )
    CALL ListAddString( Model % Solvers(n+2) % Values, 'Equation', 'RMC' )
    CALL ListAddString( Model % Solvers(n+2) % Values, 'Exec Solver', 'Never' )

    DO i=1,Model % NumberOFEquations
      IF ( ListGetLogical( Model % Equations(i) % Values, Eq, Found ) ) THEN
        CALL ListAddLogical( Model % Equations(i) % Values, 'HP', .TRUE. )
        CALL ListAddLogical( Model % Equations(i) % Values, 'RMC' , .TRUE.)
      ELSE
        ActiveSolvers => ListGetIntegerArray( CurrentModel % Equations(i) % Values, &
                              'Active Solvers', Found )
        IF ( Found ) THEN
          DO k=1,SIZE(ActiveSolvers)
            IF ( ActiveSolvers(k) == Solver % SolverId ) THEN
              CALL ListAddLogical( Model % Equations(i) % Values, 'HP', .TRUE. )
              CALL ListAddLogical( Model % Equations(i) % Values, 'RMC', .TRUE. )
              EXIT
            END IF
          END DO
        END IF
      END IF
    END DO
  END IF

    
!------------------------------------------------------------------------------
END SUBROUTINE WhitneyAVHarmonicSolver_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE WhitneyAVHarmonicSolver_Init(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: Found 
  
  Mesh => GetMesh()
  IF( Mesh % MeshDim /= 3 ) THEN
    CALL Fatal('WhitneyAVHarmonicSolver_Init','Solver requires 3D mesh!')
  END IF

  ! Historically a real array could be used for H-B Curve.
  ! This dirty piece of code makes things backward compatible.
  BLOCK
    INTEGER :: i
    LOGICAL :: Cubic
    TYPE(ValueList_t), POINTER :: Material
    DO i=1,Model % NumberOfMaterials
      Material => Model % Materials(i) % Values
      IF( ListCheckPresent( Material, 'H-B Curve') ) THEN
        Cubic = GetLogical( Material, 'Cubic spline for H-B curve',Found)
        CALL ListRealArrayToDepReal(Material,'H-B Curve','dummy',&
            CubicTable=Cubic) !Monotone=.TRUE.)         
      END IF
    END DO
  END BLOCK
  
!------------------------------------------------------------------------------
END SUBROUTINE WhitneyAVHarmonicSolver_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Solve a vector potential A and scalar potential V from
! 
!>  j omega sigma A+rot (1/mu) rot A+sigma grad(V) = J^s+rot M^s-sigma grad(V^s)
!>  -div(sigma (j omega A+grad(V)))=0
!
!>  by using edge elements (Nedelec) + nodal basis for V.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE WhitneyAVHarmonicSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils
  USE CircuitUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found, L1
  LOGICAL :: Stat, TG, Jfix, JfixSolve, LaminateStack, CoilBody, EdgeBasis,LFact,LFactFound
  LOGICAL :: PiolaVersion, SecondOrder, GotHbCurveVar, HasTensorReluctivity
  LOGICAL :: ExtNewton, StrandedHomogenization
  LOGICAL, ALLOCATABLE, SAVE :: TreeEdges(:)

  INTEGER :: n,nb,nd,t,istat,i,j,k,l,nNodes,Active,FluxCount=0
  INTEGER :: NoIterationsMin, NoIterationsMax
  INTEGER :: NewtonIter
  INTEGER, POINTER :: Perm(:)
  INTEGER, ALLOCATABLE :: FluxMap(:)

  COMPLEX(kind=dp) :: Aval
  COMPLEX(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), FORCE(:), JFixFORCE(:),JFixVec(:,:)
  COMPLEX(KIND=dp), ALLOCATABLE :: LOAD(:,:), Acoef(:), Tcoef(:,:,:)
  COMPLEX(KIND=dp), ALLOCATABLE :: LamCond(:)
  COMPLEX(KIND=dp), POINTER :: Acoef_t(:,:,:) => NULL()

  REAL(KIND=dp) :: Norm, Omega
  REAL(KIND=dp), ALLOCATABLE :: RotM(:,:,:), GapLength(:), MuParameter(:), SkinCond(:), ReLoad(:,:)
  REAL(KIND=dp), POINTER :: Cwrk(:,:,:), Cwrk_im(:,:,:), LamThick(:)
  REAL(KIND=dp), POINTER :: sValues(:), fixpot(:)
  REAL(KIND=dp) :: NewtonTol
  
  CHARACTER(LEN=MAX_NAME_LEN):: LaminateStackModel, CoilType, HbCurveVarName

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t),POINTER :: Element, Edge
  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC, BodyParams, SolverParams
  TYPE(Variable_t), POINTER :: jfixvar, jfixvarIm, HbCurveVar
  TYPE(Matrix_t), POINTER :: A
  TYPE(ListMatrix_t), POINTER :: BasicCycles(:)
  TYPE(ValueList_t), POINTER :: CompParams

  CHARACTER(LEN=MAX_NAME_LEN):: CoilCurrentName
  TYPE(Variable_t), POINTER :: CoilCurrentVar
  REAL(KIND=dp) :: CurrAmp
  LOGICAL :: UseCoilCurrent, ElemCurrent, ElectroDynamics, EigenSystem
  TYPE(Solver_t), POINTER :: pSolver 

  
  SAVE MASS, STIFF, LOAD, FORCE, Tcoef, JFixVec, JFixFORCE, Acoef, Acoef_t, &
     Cwrk, Cwrk_im, LamCond, LamThick, AllocationsDone, RotM, GapLength, MuParameter, SkinCond
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

  CALL Info('WhitneyAVHarmonicSolver','',Level=6 )
  CALL Info('WhitneyAVHarmonicSolver','------------------------------------------------',Level=6 )
  CALL Info('WhitneyAVHarmonicSolver','Solving harmonic AV equations with edge elements',Level=5 )
   
  SolverParams => GetSolverParams()
  pSolver => Solver
  
  EigenSystem = GetLogical( SolverParams, 'Eigen Analysis', Found )
  ElectroDynamics = GetLogical( SolverParams, 'Electrodynamics Model', Found )

  CALL EdgeElementStyle(SolverParams, PiolaVersion, QuadraticApproximation = SecondOrder )
  
  IF (PiolaVersion) THEN
    CALL Info('WhitneyAVHarmonicSolver', &
        'Using Piola Transformed element basis functions',Level=4)
    CALL Info('WhitneyAVHarmonicSolver', &
        'The option > Use Tree Gauge < is not available',Level=4)
  END IF

  CoilCurrentName = GetString( SolverParams,'Current Density Name',UseCoilCurrent ) 
  IF(.NOT. UseCoilCurrent ) THEN
    UseCoilCurrent = GetLogical(SolverParams,'Use Nodal CoilCurrent',Found )
    IF(UseCoilCurrent) THEN
      CoilCurrentName = 'CoilCurrent'
    ELSE
      UseCoilCurrent = GetLogical(SolverParams,'Use Elemental CoilCurrent',Found )
      IF(UseCoilCurrent) CoilCurrentName = 'CoilCurrent e'
    END IF
  END IF
  ElemCurrent = .FALSE.

  IF( UseCoilCurrent ) THEN
    CoilCurrentVar => VariableGet(Solver % Mesh % Variables, CoilCurrentName )
    IF( ASSOCIATED( CoilCurrentVar ) ) THEN
      CALL Info('WhitneyAVHarmonicSolver','Using precomputed field for current density: '//TRIM(CoilCurrentName),Level=5)
      IF( CoilCurrentVar % TYPE == Variable_on_nodes_on_elements ) THEN
        ElemCurrent = .TRUE.
      ELSE
        CALL Warn('WhitneyAVHarmonicSolver','Precomputed CoilCurrent is not an elemental field!')
      END IF
    ELSE
      CALL Fatal('WhitneyAVHarmonicSolver','Elemental current requested but not found:'//TRIM(CoilCurrentName))
    END IF
  END IF

  
  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  Mesh => GetMesh()
  nNodes = Mesh % NumberOfNodes
  Perm => Solver % Variable % Perm

  A => GetMatrix()

  IF ( .NOT. AllocationsDone ) THEN

     IF (Solver % Variable % dofs /= 2) CALL Fatal('WhitneyAVHarmonicSolver', &
         'Variable is not properly defined for time harmonic AV solver, Use: Variable = A[A re:1 A im:1]')

     N = Mesh % MaxElementDOFs  ! just big enough
     ALLOCATE( FORCE(N), LOAD(7,N), ReLOAD(3,N), STIFF(N,N), MASS(n,n), &
          JFixVec(3,N),JFixFORCE(n), Tcoef(3,3,N), RotM(3,3,N), &
          GapLength(N), MuParameter(N), SkinCond(N), Acoef(N), LamCond(N), &
          LamThick(N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'WhitneyAVHarmonicSolver', 'Memory allocation error.' )
     END IF

     NULLIFY( Cwrk )
     NULLIFY( Cwrk_im )

     AllocationsDone = .TRUE.
  END IF
  
  Omega = GetAngularFrequency(Found=Found)
  IF(.NOT. Found .AND. .NOT. EigenSystem ) THEN
    CALL Fatal('WhitneyHarmonicAVSolver','Harmonic solution requires frequency!')
  END IF
    
  Jfix = GetLogical(SolverParams,'Fix input Current Density', Found)
  IF(.NOT. Found ) THEN
    ! If not specified compute the Jfix field only if there is a specified current BC
    Jfix = ListCheckPrefixAnyBodyForce( Model,'Current Density' )
  END IF
  JfixSolve = Jfix

  IF (Jfix) THEN
    JfixPhase = 1
    CALL JfixPotentialSolver(Model,Solver,dt,Transient)
    JfixVar => VariableGet(Mesh % Variables, 'Jfix')    
    JfixVarIm => VariableGet(Mesh % Variables, 'Jfix Im')    
    IF(.NOT. ASSOCIATED( JfixRhsC ) ) THEN
      CALL Fatal('WhitneyAVHarmonicSolver','JfixRhsC should be associated!')
    END IF
    IF(.NOT. ASSOCIATED( JFixSurfacePerm ) ) THEN
      CALL Fatal('WhitneyAVHarmonicSolver','JFixVecSurfacePerm should be associated!')
    END IF
    IF(.NOT. ALLOCATED( JFixSurfaceVecC ) ) THEN
      CALL Fatal('WhitneyAVHarmonicSolver','JFixVecSurfaceVecC should be associated!')
    END IF   
  END IF  
  
  HbCurveVarName = GetString( SolverParams,'H-B Curve Variable', GotHbCurveVar )
  IF( GotHbCurveVar ) THEN
    HbCurveVar => VariableGet( Mesh % Variables, HbCurveVarName )
    IF(.NOT. ASSOCIATED( HbCurveVar ) ) THEN
      CALL Fatal('WhitneyAVHarmonicSolver','H-B Curve variable given but does not exist: '&
          //TRIM(HbCurveVarName))
    END IF
  END IF
    
  ! Resolve internal non.linearities, if requested:
  ! ----------------------------------------------
  NoIterationsMax = GetInteger( SolverParams, &
      'Nonlinear System Max Iterations',Found)
  IF(.NOT. Found) NoIterationsMax = 1
  
  NoIterationsMin = GetInteger( SolverParams, &
      'Nonlinear System Min Iterations',Found)
  IF(.NOT. Found) NoIterationsMin = 1

  ! Use also these keyword for compatibility with ElmerGUI and old practices
  NewtonIter = GetInteger( SolverParams,&
      'Nonlinear System Newton After Iterations',Found ) 
  IF(.NOT. Found ) NewtonIter = NoIterationsMax
  NewtonTol = GetCReal( SolverParams,&
      'Nonlinear System Newton After Tolerance',Found )

  
  LFact = GetLogical( SolverParams,'Linear System Refactorize', LFactFound )
  EdgeBasis = .NOT. LFactFound .AND. GetLogical( SolverParams, 'Edge Basis', Found )

  CALL DefaultStart()

  DO i=1,NoIterationsMax
    ExtNewton = ( i > NewtonIter .OR. Solver % Variable % NonlinChange < NewtonTol )

    IF( DoSolve(i) ) THEN
      IF(i>=NoIterationsMin) EXIT
    END IF
    IF( EdgeBasis ) CALL ListAddLogical(SolverParams,'Linear System Refactorize',.FALSE.)

    JFixSolve = .FALSE.
  END DO
  IF ( EdgeBasis ) CALL ListRemove( SolverParams, 'Linear System Refactorize' )

  CALL CalculateLumped(Model % NumberOfBodyForces)

  CALL DefaultFinish()

  
CONTAINS

!---------------------------------------------------------------------------------------------
  FUNCTION DoSolve(IterNo) RESULT(Converged)
!---------------------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: IterNo
    LOGICAL :: Converged
!---------------------------------------------------------------------------------------------
    REAL(KIND=dp) :: Norm, PrevNorm, TOL
    INTEGER :: i,j,k,n,nd,t,ComponentId
    REAL(KIND=dp), ALLOCATABLE :: Diag(:)
    LOGICAL  :: FoundMagnetization, Found, ConstraintActive, GotCoil, CircuitDrivenBC
!---------------------------------------------------------------------------------------------
    ! System assembly:
    !-----------------
    CALL DefaultInitialize()
    Active = GetNOFActive()

    DO t=1,active
       Element => GetActiveElement(t)
       n  = GetElementNOFNodes() ! vertices
       nd = GetElementNOFDOFs()  ! dofs
       
       IF (SIZE(Tcoef,3) /= n) THEN
         DEALLOCATE(Tcoef)
         ALLOCATE(Tcoef(3,3,n))
       END IF
       
       LOAD = 0.0d0
       BodyForce => GetBodyForce()
       FoundMagnetization = .FALSE.

       ! If the coil current field is elemental it is discontinuous and need not be limited
       ! to the body force. For nodal ones we don't have the same luxury.
       GotCoil = .FALSE.
       IF( UseCoilCurrent ) THEN
         IF( ElemCurrent .OR. ASSOCIATED(BodyForce) ) THEN
           CALL GetVectorLocalSolution( ReLoad,UElement=Element,UVariable=CoilCurrentVar,Found=GotCoil)       
           LOAD(1:3,1:n) = ReLoad(1:3,1:n)
         END IF
       END IF
              
       IF ( ASSOCIATED(BodyForce) ) THEN
         ! If not already given by CoilCurrent, request for current density
         IF( .NOT. GotCoil ) THEN           
           CALL GetComplexVector( BodyForce, Load(1:3,1:n), 'Current Density', Found )
         END IF

         CurrAmp = ListGetCReal( BodyForce,'Current Density Multiplier',Found ) 
         IF(Found) Load(1:3,1:n) = CurrAmp * Load(1:3,1:n)
                    
         CALL GetComplexVector( BodyForce, Load(4:6,1:n), &
             'Magnetization', FoundMagnetization )

         Load(7,1:n) = GetReal( BodyForce, 'Electric Potential', Found )
         Load(7,1:n) = CMPLX( REAL(Load(7,1:n)), &
             GetReal( BodyForce, 'Electric Potential im', Found), KIND=dp)
       END IF

       Material => GetMaterial( Element )

       IF(ASSOCIATED(Material).AND..NOT.FoundMagnetization) THEN
          CALL GetComplexVector( Material, Load(4:6,1:n), &
                 'Magnetization', FoundMagnetization )
       END IF

       CoilBody = .FALSE.
       CompParams => GetComponentParams( Element )
       CoilType = ''
       RotM = 0._dp
       ConstraintActive = .TRUE.
       IF (ASSOCIATED(CompParams)) THEN
         CoilType = GetString(CompParams, 'Coil Type', Found)
         IF (Found) THEN
           SELECT CASE (CoilType)
           CASE ('stranded')
              CoilBody = .TRUE.
              StrandedHomogenization = GetLogical(CompParams, 'Homogenization Model', Found)
              IF( StrandedHomogenization ) THEN
                CALL GetElementRotM(Element, RotM, n)
              END IF
           CASE ('massive')
              CoilBody = .TRUE.
           CASE ('foil winding')
              CoilBody = .TRUE.
              CALL GetElementRotM(Element, RotM, n)
           CASE DEFAULT
              CALL Fatal ('WhitneyAVHarmonicSolver', 'Non existent Coil Type Chosen!')
           END SELECT
         END IF
         ConstraintActive = GetLogical( CompParams, 'Activate Constraint', Found)
!        IF(.NOT.Found .AND. CoilType /= 'stranded') ConstraintActive = .TRUE.
         IF(.NOT.Found ) ConstraintActive = .FALSE.
       END IF

       LaminateStack = .FALSE.
       LaminateStackModel = ''
       HasTensorReluctivity = .FALSE.
       Acoef = 0.0_dp
       Tcoef = 0.0_dp
       IF ( ASSOCIATED(Material) ) THEN
         IF (.NOT. ListCheckPresent(Material, 'H-B Curve')) THEN
           CALL GetReluctivity(Material,Acoef_t,n,HasTensorReluctivity)
           IF (HasTensorReluctivity) THEN
             IF (size(Acoef_t,1)==1 .AND. size(Acoef_t,2)==1) THEN
               Acoef(1:n) = Acoef_t(1,1,1:n) 
               HasTensorReluctivity = .FALSE.
             ELSE IF (size(Acoef_t,1)/=3) THEN
               CALL Fatal('WhitneyAVHarmonicSolver', 'Reluctivity tensor should be of size 3x3')
             END IF
           ELSE
             CALL GetReluctivity(Material,Acoef,n)
           END IF
         END IF
!------------------------------------------------------------------------------
!        Read conductivity values (might be a tensor)
!------------------------------------------------------------------------------
         Tcoef = GetCMPLXElectricConductivityTensor(Element, n, CoilBody, CoilType) 

         LaminateStackModel = GetString( Material, 'Laminate Stack Model', LaminateStack )
       END IF

       LamThick=0d0
       LamCond=0d0
       IF (LaminateStack) THEN
         SELECT CASE(LaminateStackModel)
         CASE('low-frequency model')
           LamThick(1:n) = GetReal( Material, 'Laminate Thickness', Found )
           IF (.NOT. Found) CALL Fatal('WhitneyAVHarmonicSolver', 'Laminate Thickness not found!')
  
           LamCond(1:n) = GetReal( Material, 'Laminate Stack Conductivity', Found )
           IF (.NOT. Found) CALL Fatal('WhitneyAVHarmonicSolver', 'Laminate Stack Conductivity not found!')
           LamCond(1:n) = CMPLX( REAL(LamCond(1:n)), &
                  GetReal( Material, 'Electric Conductivity im',  Found), KIND=dp)
         CASE('wide-frequency-band model')
           LamThick(1:n) = GetReal( Material, 'Laminate Thickness', Found )
           IF (.NOT. Found) CALL Fatal('WhitneyAVHarmonicSolver', 'Laminate Thickness not found!')

           LamCond(1:n) = GetReal( Material, 'Laminate Stack Conductivity', Found )
           IF (.NOT. Found) CALL Fatal('WhitneyAVHarmonicSolver', 'Laminate Stack Conductivity not found!')
           LamCond(1:n) = CMPLX( REAL(LamCond(1:n)), &
                  GetReal( Material, 'Electric Conductivity im',  Found), KIND=dp)
         CASE DEFAULT
           CALL WARN('WhitneyAVHarmonicSolver', 'Nonexistent Laminate Stack Model chosen!')
         END SELECT
       END IF

       Omega = GetAngularFrequency(Found=Found,UElement=Element)

       !Get element local matrix and rhs vector:
       !----------------------------------------
       CALL LocalMatrix( MASS, STIFF, FORCE, JFixFORCE, JFixVec, LOAD, &
          Tcoef, Acoef, LaminateStack, LaminateStackModel, LamThick, &
          LamCond, CoilBody, CoilType, RotM, ConstraintActive, Element, n, nd, PiolaVersion, SecondOrder )

       !Update global matrix and rhs vector from local matrix & vector:
       !---------------------------------------------------------------
       CALL DefaultUpdateEquations( STIFF, FORCE )
       IF (EigenSystem) CALL DefaultUpdateMass(MASS)
       
       ! Memorize stuff for the fixing potential
       ! 1) Divergence of the source term
       ! 2) The source terms at the surface to determine the direction
       !-------------------------------------------------------------------
       IF( JFixSolve ) THEN
         JFixRhsC(JFixVar % Perm(Element % NodeIndexes)) = &
             JFixRhsC(JFixVar % Perm(Element % NodeIndexes)) + JFixFORCE(1:n)
         DO i=1,n
           j = JfixSurfacePerm(Element % NodeIndexes(i) )         
           IF( j > 0 ) JfixSurfaceVecC(3*j-2:3*j) = &
               JfixSurfaceVecC(3*j-2:3*j) + JFixVec(1:3,i)
         END DO
       END IF
     END DO
    
     IF( JfixSolve ) THEN    
       CALL Info('WhitneyAVHarmonicSolver','Solving the fixing potential')
       JfixPhase = 2 
       CALL JfixPotentialSolver(Model,Solver,dt,Transient)

       CALL Info('WhitneyAVHarmonicSolver','Adding the fixing potential to the r.h.s. of AV equation')   
       DO t=1,active
         Element => GetActiveElement(t)
         n  = GetElementNOFNodes() 
         nd = GetElementNOFDOFs()  
         nb = GetElementNOFBDOFs() 

         CALL LocalFixMatrixC( FORCE, Element, n, nd+nb, PiolaVersion, SecondOrder)      
       END DO
       CALL Info('WhitneyAVHarmonicSolver','Finished adding the fixing potential',Level=10)   
     END IF

    
    ! Robin type of BC in terms of H:
    !--------------------------------
    Active = GetNOFBoundaryElements()
    DO t=1,Active
       Element => GetBoundaryElement(t)
       BC=>GetBC()
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

       Model % CurrentElement => Element
       nd = GetElementNOFDOFs(Element)
       n  = GetElementNOFNodes(Element)

       CALL GetComplexVector( BC, Load(1:3,1:n), 'Magnetic Field Strength', Found)

       Load(4,1:n) = GetReal( BC, 'Electric Current Density', Found )
       IF (.NOT. Found) Load(4,1:n) = GetReal( BC, 'Electric Flux', Found )

       Load(4,1:n) = CMPLX( REAL(LOAD(4,1:n)), &
            GetReal( BC, 'Electric Current Density im', Found), KIND=dp)

       IF (.NOT. Found) Load(4,1:n) = CMPLX( REAL(LOAD(4,1:n)), &
            GetReal( BC, 'Electric Flux im', Found), KIND=dp)

       Load(5,1:n) = GetReal( BC, 'Electric Transfer Coefficient', Found )
       Load(5,1:n) = CMPLX( REAL(Load(5,1:n)), &
         GetReal( BC, 'Electric Transfer Coefficient im', Found), KIND=dp)

       Acoef(1:n) = GetReal( BC, 'Magnetic Transfer Coefficient', Found )
       Acoef(1:n) = CMPLX( REAL(Acoef(1:n)), &
         GetReal( BC, 'Magnetic Transfer Coefficient im', Found), KIND=dp)

       !If air gap length keyword is detected, use air gap boundary condition
       GapLength = GetConstReal( BC, 'Air Gap Length', Found)
       IF (Found) THEN
         MuParameter=GetConstReal( BC, 'Air Gap Relative Permeability', Found)
         IF (.NOT. Found) MuParameter = 1.0_dp ! if not found default to "air" property
         CALL LocalMatrixAirGapBC(MASS,STIFF,FORCE,LOAD,GapLength,MuParameter,Element,n,nd )
       ELSE
         SkinCond = GetConstReal( BC, 'Layer Electric Conductivity', Found)
         IF (ANY(ABS(SkinCond(1:n)) > AEPS)) THEN
           MuParameter = GetConstReal( BC, 'Layer Relative Permeability', Found)
           ComponentId=GetInteger( BC, 'Component', CircuitDrivenBC)
           IF (.NOT. Found) MuParameter = 1.0_dp ! if not found default to "air" property           
           CALL LocalMatrixSkinBC(MASS,STIFF,FORCE,SkinCond,MuParameter,Element,CircuitDrivenBC,n,nd)
         ELSE         
           GapLength = GetConstReal( BC, 'Thin Sheet Thickness', Found)
           IF (Found) THEN
             MuParameter=GetConstReal( BC, 'Thin Sheet Relative Permeability', Found)
             IF (.NOT. Found) MuParameter = 1.0_dp ! if not found default to "air" property
             ! Technically, there is no skin but why create yet another conductivity variable?
             SkinCond = GetConstReal( BC, 'Thin Sheet Electric Conductivity', Found)
             IF (.NOT. Found) SkinCond = 1.0_dp ! if not found default to "air" property
             CALL LocalMatrixThinSheet( MASS, STIFF, FORCE, LOAD, GapLength, MuParameter, &
                                            SkinCond, Element, n, nd )
           ELSE
             CALL LocalMatrixBC(MASS,STIFF,FORCE,LOAD,Acoef,Element,n,nd )
           END IF
         END IF
       END IF

       CALL DefaultUpdateEquations(STIFF,FORCE,Element)
       IF(EigenSystem) CALL DefaultUpdateMass(MASS,Element)
    END DO

    CALL DefaultFinishAssembly()

    !
    ! Check for tree gauge, if requested or using direct solver:
    ! ------------------------------------------------------------
    TG=GetLogical(SolverParams, 'Use tree gauge', Found)
    IF (.NOT. Found) TG=GetString(GetSolverParams(), &
        'Linear System Solver',Found)=='direct'

    !
    ! Dirichlet BCs in terms of vector potential A:
    ! ---------------------------------------------
    IF ( TG ) THEN
      ! temporary fix to some scaling problem (to be resolved)...
      CALL ListAddLogical( SolverParams, 'Linear System Dirichlet Scaling', .FALSE.) 
    END IF


BLOCK
! Automatic BC for massive,foil coils outer boundaries, when "Activate Constraint" on!!

    TYPE(Element_t), POINTER :: Parent
    LOGICAL :: AutomaticBC
    INTEGER, POINTER :: Electrodes(:)

    A => GetMatrix()

    IF (.NOT.ALLOCATED(A % ConstrainedDOF)) THEN
      ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
      A % ConstrainedDOF = .FALSE.
    END IF

    IF(.NOT.ALLOCATED(A % DValues)) THEN
      ALLOCATE(A % Dvalues(A % NumberOfRows))
      A % Dvalues = 0._dp
    END IF

    Active = GetNOFBoundaryElements()
    DO t = 1, Active
      Element => GetBoundaryElement(t)

      Parent => Element % BoundaryInfo % Right
      IF(ASSOCIATED(Parent)) CYCLE

      IF(ParEnv % PEs>1) THEN
        ! Assuming here that this is an internal boundary, if all elements nodes are
        ! interface nodes. Not foolproof i guess, but quite safe (?)
        IF (ALL(Solver % Mesh % ParallelInfo % GInterface(Element % NodeIndexes))) CYCLE
      END IF
 
      Parent => Element % BoundaryInfo % Left
      IF(.NOT.ASSOCIATED(Parent)) CYCLE

      CompParams => GetComponentParams(Parent)
      IF (.NOT. ASSOCIATED(CompParams)) CYCLE

      CoilType = GetString(CompParams, 'Coil Type', Found)
      IF(CoilType/='massive' .AND. CoilType/='foil winding') CYCLE

      ConstraintActive = GetLogical(CompParams,'Activate Constraint',Found )
      IF( .NOT. ConstraintActive ) CYCLE

      AutomaticBC = GetLogical( CompParams, 'Automatic electrode BC', Found )
      IF(.NOT.Found) AutomaticBC = .TRUE.

      IF(.NOT.AutomaticBC) CYCLE

      Electrodes =>  ListGetIntegerArray( CompParams, &
             'Electrode Boundaries', Found )

      IF(ASSOCIATED(Electrodes)) THEN
        IF(ALL(Electrodes/=Element % BoundaryInfo % Constraint)) CYCLE
      END IF

      DO i=1,Element % Type % NumberOfNodes
        j = 2*(Solver % Variable % Perm(Element % NodeIndexes(i))-1)+1
        A % ConstrainedDOF(j:j+1) = .TRUE.
      END DO
    END DO
END BLOCK

    CALL DefaultDirichletBCs()

    !
    ! Dirichlet BCs in terms of magnetic flux density B:
    ! --------------------------------------------------
    CALL DirichletAfromB()


    A => GetMatrix()
    IF (TG) THEN
      IF(.NOT.ALLOCATED(TreeEdges)) &
         CALL GaugeTree(Solver,Mesh,TreeEdges,FluxCount,FluxMap,Transient)

      WRITE(Message,*) 'Volume tree edges: ', &
          i2s(COUNT(TreeEdges)),     &
          ' of total: ',Mesh % NumberOfEdges
      CALL Info('WhitneyAVHarmonicSolver: ', Message, Level=5)

      DO i=1,SIZE(TreeEdges)
        IF(TreeEdges(i)) CALL SetDOFToValue(Solver,i,(0._dp,0._dp))
      END DO
    END IF

    !
    ! Fix unused potential DOFs:
    ! --------------------------
    CALL ConstrainUnused(A)

    !
    ! Linear system solution:
    ! -----------------------
    Norm = DefaultSolve()
    Converged = Solver % Variable % NonlinConverged==1
!------------------------------------------------------------------------------
  END FUNCTION DoSolve
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ConstrainUnused(A)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t) :: A
!------------------------------------------------------------------------------
    INTEGER :: i,j,n
   
    REAL(KIND=dp), ALLOCATABLE :: dDiag(:)
!------------------------------------------------------------------------------
    n = A % NumberOFRows
    ALLOCATE(dDiag(n)); dDiag=0._dp

    DO i=1,n,2
      j = A % Diag(i)
      IF(j>0) THEN
        dDiag(i)   =  A % Values(j)
        dDiag(i+1) = -A % Values(j+1)
      END IF
    END DO
    IF (ParEnv % PEs>1) CALL ParallelSumVector(A, dDiag)

    n = Mesh % NumberOfNodes
    DO i=1,SIZE(Solver % Variable % Perm) !n
      j = Solver % Variable % Perm(i)
      IF (j==0) CYCLE

      j = 2*(j-1)
      Aval = CMPLX(dDiag(j+1), dDiag(j+2), KIND=dp)

      IF (ABS(Aval)==0._dp) THEN
        A % RHS(j+1) = 0._dp
        CALL ZeroRow(A,j+1)
        A % Values(A % Diag(j+1)) = 1._dp

        A % RHS(j+2) = 0._dp
        CALL ZeroRow(A,j+2)
        A % Values(A % Diag(j+2)) = 1._dp

        IF(ALLOCATED(A % ConstrainedDOF)) THEN
          A % ConstrainedDOF(j+1) = .TRUE.
          A % ConstrainedDOF(j+2) = .TRUE.
        END IF
      END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ConstrainUnused
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE CalculateLumped(nbf)
!------------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER::nbf
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: a(nbf),IMoment,IA
   INTEGER :: i,bfid,n,nd,EdgeBasisDegree
   TYPE(Element_t), POINTER :: Element, Parent
   COMPLEX(KIND=dp) :: torq, U(nbf), zforce, zzforce
   LOGICAL :: Found, CalcTorque,CalcPotential,CalcInertia
   TYPE(ValueList_t),POINTER::Params
!------------------------------------------------------------------------------

   CalcTorque = ListCheckPresentAnyBody(Model,'r inner')
   CalcPotential = ListGetLogicalAnyBodyForce( Model,'Calculate Potential')
   CalcInertia = ListGetLogicalAnyBody( Model,'Calculate Inertial Moment')

   IF(.NOT. (CalcTorque .OR. CalcPotential .OR. CalcInertia ) ) RETURN

   EdgeBasisDegree = 1
   IF (SecondOrder) EdgeBasisDegree = 2

   U=0._dp; a=0._dp; torq=0._dp; IMoment=0._dp;IA=0; zforce=0
   DO i=1,GetNOFActive()
     Element => GetActiveElement(i)
     nd = GetElementNOFDOFs(Element)
     n  = GetElementNOFNodes(Element)

     IF( CalcTorque ) THEN
       CALL Torque(Torq,Element,n,nd,EdgeBasisDegree)
       CALL AxialForce(zforce,Element,n,nd,EdgeBasisDegree)
     END IF

     IF( CalcPotential ) THEN
       Params=>GetBodyForce(Element)
       IF(ASSOCIATED(Params)) THEN
         bfid=GetBodyForceId(Element)
         IF(GetLogical(Params,'Calculate Potential',Found)) &
             CALL Potential(u(bfid),a(bfid),Element,n,nd,EdgeBasisDegree)
       END IF
     END IF

     IF( CalcInertia ) THEN
       Params=>GetBodyParams(Element)
       IF(ASSOCIATED(Params)) THEN
         IF(GetLogical(Params,'Calculate Inertial Moment',Found)) &
             CALL InertialMoment(IMoment,IA,Element,n,nd)
       END IF
     END IF
   END DO

   zzforce = 0
   IF(ListGetLogicalAnyBC(Model,'Calculate Axial Force')) THEN
     DO i=1,Mesh % NumberOFBoundaryElements
       Element => GetBoundaryElement(i)
       IF (.NOT.GetLogical(GetBC(), 'Calculate Axial Force', Found ) ) CYCLE

       Parent => Element % BoundaryInfo % Left
       n  = GetELementNofNodes(Parent)
       nd = GetELementNofDOFs(Parent)
       CALL AxialForceSurf(zzforce,Element,n,nd,EdgeBasisDegree)
     END DO
   END IF

   IF( CalcPotential ) THEN
     DO i=1,nbf
       a(i) = ParallelReduction(a(i))
       u(i) = ParallelReduction(u(i))
     END DO
     
     DO i=1,nbf
       IF(a(i)>0) THEN
         CALL ListAddConstReal(Model % Simulation,'res: Potential re / bodyforce ' &
             //i2s(i),REAL(u(i))/a(i))
         CALL ListAddConstReal(Model % Simulation,'res: Potential im / bodyforce ' &
             //i2s(i),AIMAG(u(i))/a(i))
         CALL ListAddConstReal(Model % Simulation,'res: area / bodyforce ' &
             //i2s(i),a(i))
       END IF
     END DO
   END IF

   IF( CalcTorque ) THEN
     Torq = ParallelReduction(Torq)
     CALL ListAddConstReal(Model % Simulation,'res: Air Gap Torque re', REAL(Torq))
     CALL ListAddConstReal(Model % Simulation,'res: Air Gap Torque im', AIMAG(Torq))

     zforce = ParallelReduction(zforce)
     CALL ListAddConstReal(Model % Simulation,'res: Axial force(vol) re', REAL(zforce))
     CALL ListAddConstReal(Model % Simulation,'res: Axial force(vol) im', AIMAG(zforce))

   zzforce = 0
   IF(ListGetLogicalAnyBC(Model,'Calculate Axial Force')) THEN
     DO i=1,Mesh % NumberOFBoundaryElements
       Element => GetBoundaryElement(i)
       IF (.NOT.GetLogical(GetBC(), 'Calculate Axial Force', Found ) ) CYCLE

       Parent => Element % BoundaryInfo % Left
       n  = GetELementNofNodes(Parent)
       nd = GetELementNofDOFs(Parent)
       CALL AxialForceSurf(zzforce,Element,n,nd,EdgeBasisDegree)
     END DO
   END IF
     IF(ListGetLogicalAnyBC(Model,'Calculate Axial Force')) THEN
       zzforce = ParallelReduction(zzforce)
       CALL ListAddConstReal(Model % Simulation,'res: Axial force(surf) re', REAL(zzforce))
       CALL ListAddConstReal(Model % Simulation,'res: Axial force(surf) im', AIMAG(zzforce))
     END IF
   END IF

   IF( CalcInertia ) THEN
     IMoment = ParallelReduction(IMoment)
     IA = ParallelReduction(IA)
     CALL ListAddConstReal(Model % Simulation,'res: Inertial Volume', IA)
     CALL ListAddConstReal(Model % Simulation,'res: Inertial Moment', IMoment)
   END IF

!------------------------------------------------------------------------------
 END SUBROUTINE CalculateLumped
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE InertialMoment(U,A,Element,n,nd)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n,nd
    REAL(KIND=dp)::U,a
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), DetJ,x,y,r,Density(n)
    INTEGER :: t
    LOGICAL :: stat,Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    Density(1:n) = GetReal(GetMaterial(),'Density',Found,Element)
    IF(.NOT.Found) RETURN

    CALL GetElementNodes( Nodes, Element )
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis )

      x = SUM(Nodes % x(1:n)*Basis(1:n))
      y = SUM(Nodes % y(1:n)*Basis(1:n))
      r = SQRT(x**2+y**2)
      A = A + IP % s(t)*detJ
      U = U + IP % s(t)*detJ*R*SUM(Density(1:n)*Basis(1:n))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE InertialMoment
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE Torque(U,Element,n,nd,EdgeBasisDegree)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n,nd,EdgeBasisDegree
    COMPLEX(KIND=dp)::U
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: dBasisdx(nd,3),Basis(nd), DetJ, &
             POT(2,nd),x,y,r,r0,r1,Wbasis(nd,3),RotWBasis(nd,3)
    COMPLEX(KIND=dp) :: B(3,nd), POTC(nd), Br, Bp, Bx, By
    INTEGER :: t
    LOGICAL :: stat, Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    r0 = GetCReal(GetBodyParams(),'r inner',Found)
    r1 = GetCReal(GetBodyParams(),'r outer',Found)
    IF (.NOT.Found) RETURN

    CALL GetElementNodes( Nodes, Element )

    x = SUM(Nodes % x(1:n))/n
    y = SUM(Nodes % y(1:n))/n
    r = SQRT(x**2+y**2)
    IF (r<r0.OR.r>r1) RETURN

    CALL GetLocalSolution(POT, UElement=Element)
    POTC = CMPLX( POT(1,1:nd), POT(2,1:nd) )
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = RotWBasis, USolver = pSolver )

      x = SUM(Nodes % x(1:n)*Basis(1:n))
      y = SUM(Nodes % y(1:n)*Basis(1:n))
      r = SQRT(x**2+y**2)

      Bx =  SUM(POTC(n+1:nd) * RotWBasis(1:nd-n,1))
      By =  SUM(POTC(n+1:nd) * RotWBasis(1:nd-n,2))
      Br =  x/r*Bx + y/r*By
      Bp = -y/r*Bx + x/r*By
      U = U + IP % s(t) * detJ * r * &
           CMPLX(REAL(Br)*REAL(Bp),AIMAG(Br)*AIMAG(Bp))/(PI*4.0d-7*(r1-r0))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE Torque
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AxialForce(U,Element,n,nd,EdgeBasisDegree)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n,nd,EdgeBasisDegree
    COMPLEX(KIND=dp)::U
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: dBasisdx(nd,3),Basis(nd), DetJ, &
             POT(2,nd),x,y,r,r0,r1,Wbasis(nd,3),RotWBasis(nd,3)
    COMPLEX(KIND=dp) :: B(3,nd), POTC(nd), Bx, By, Bz, Br, Bp
    INTEGER :: t
    LOGICAL :: stat, Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    r0 = GetCReal(GetBodyParams(),'r inner',Found)
    r1 = GetCReal(GetBodyParams(),'r outer',Found)
    IF (.NOT.Found) RETURN

    CALL GetElementNodes( Nodes, Element )

    x = SUM(Nodes % x(1:n))/n
    y = SUM(Nodes % y(1:n))/n
    r = SQRT(x**2+y**2)
    IF (r<r0.OR.r>r1) RETURN

    CALL GetLocalSolution(POT, UElement=Element)
    POTC = CMPLX( POT(1,1:nd), POT(2,1:nd) )
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = RotWBasis, USolver = pSolver )
      
      x = SUM(Nodes % x(1:n)*Basis(1:n))
      y = SUM(Nodes % y(1:n)*Basis(1:n))
      r = SQRT(x**2+y**2)
      x = x/r; y=y/r

      Bx =  SUM(POTC(n+1:nd) * RotWBasis(1:nd-n,1))
      By =  SUM(POTC(n+1:nd) * RotWBasis(1:nd-n,2))
      Bz =  SUM(POTC(n+1:nd) * RotWBasis(1:nd-n,3))
      U = U + IP % s(t) * detJ * 1 * &
           CMPLX((REAL(Bx)*REAL(Bz)*x + REAL(By)*REAL(Bz)*y), &
                (AIMAG(Bx)*AIMAG(Bz)*x + AIMAG(By)*AIMAG(Bz)*y)) &
               /(PI*4.0d-7*(r1-r0))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AxialForce
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AxialForceSurf(U,Element,n,nd,EdgeBasisDegree)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n,nd,EdgeBasisDegree
    COMPLEX(KIND=dp)::U
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER ::PARENT
    REAL(KIND=dp) :: dBasisdx(nd,3),Basis(nd), DetJ, Pdetj, uu,v,w, &
             POT(2,nd),x,y,r,r0,r1,Wbasis(nd,3),RotWBasis(nd,3)
    COMPLEX(KIND=dp) :: B(3,nd), POTC(nd), Bx, By, Bz
    INTEGER :: t
    LOGICAL :: stat, Found
    TYPE(Nodes_t), SAVE :: Nodes, PNodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    CALL GetElementNodes( Nodes, Element )
    Parent => Element % BoundaryInfo % Left
    CALL GetElementNodes( PNodes, Parent )

    CALL GetLocalSolution(POT, UElement=Parent )
    POTC = CMPLX( POT(1,1:nd), POT(2,1:nd) )
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis, dBasisdx )

      CALL GetParentUVW(Element,GetElementNOFNodes(Element),Parent,n,uu,v,w,Basis)
      IF (PiolaVersion) THEN 
        stat = EdgeElementInfo( Parent, PNodes, uu, v, w, &
              DetF = PDetJ, Basis = Basis, EdgeBasis = WBasis, RotBasis = RotWBasis, &
              BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
      ELSE
        stat = ElementInfo( Parent, PNodes, uu,v,w, pdetJ, Basis, dBasisdx )
        CALL GetEdgeBasis(Parent,WBasis,RotWBasis,Basis,dBasisdx)
      END IF

      x = SUM(Basis(1:n) * PNodes % x(1:n))
      y = SUM(Basis(1:n) * PNodes % y(1:n))
      r = SQRT(x**2 + y**2)
      x=x/r; y=y/r

      Bx =  SUM(POTC(n+1:nd) * RotWBasis(1:nd-n,1))
      By =  SUM(POTC(n+1:nd) * RotWBasis(1:nd-n,2))
      Bz =  SUM(POTC(n+1:nd) * RotWBasis(1:nd-n,3))
      U = U + IP % s(t) * detJ * &
           CMPLX((REAL(Bx)*REAL(Bz)*x + REAL(By)*REAL(Bz)*y), &
                (AIMAG(Bx)*AIMAG(Bz)*x + AIMAG(By)*AIMAG(Bz)*y)) /(PI*4.0d-7)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AxialForceSurf
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Potential( U, A, Element,n,nd,EdgeBasisDegree)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: A
    COMPLEX(KIND=dp) :: U
    INTEGER :: n, nd, EdgeBasisDegree
    TYPE(Element_t), POINTER :: Element

    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3),DetJ,POT(2,nd), &
          wBasis(nd,3),rotWBasis(nd,3),Wpot(nd),w(3), Omega
    COMPLEX(KIND=dp) :: POTC(nd)
    INTEGER :: t
    LOGICAL :: stat, WbaseFound
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    CALL GetElementNodes( Nodes )

    Omega = GetAngularFrequency(UElement=Element)
    CALL GetLocalSolution(POT,UElement=Element)
    POTC = Omega*CMPLX( POT(2,1:nd), POT(1,1:nd) )

    CALL GetLocalSolution(Wpot,'W',UElement=Element)
    W = [0._dp, 0._dp, 1._dp]
    WbaseFound = ANY(Wpot(1:n)/=0._dp)

    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = RotWBasis, USolver = pSolver )
     
      IF(WBaseFound) W = MATMUL(Wpot(1:n),dBasisdx(1:n,:))

      A = A + IP % s(t) * detJ
      U = U + IP % s(t) * detJ * SUM(PotC(n+1:nd)*MATMUL(WBasis(1:nd-n,:),w))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE Potential
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( MASS, STIFF, FORCE, JFixFORCE, JFixVec, LOAD, &
      Tcoef, Acoef, LaminateStack, LaminateStackModel, & 
      LamThick, LamCond, CoilBody, CoilType, RotM, ConstraintActive, &
      Element, n, nd, PiolaVersion, SecondOrder )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:), JFixFORCE(:), JFixVec(:,:)
    COMPLEX(KIND=dp) :: LOAD(:,:), Tcoef(:,:,:), Acoef(:), LamCond(:)
    REAL(KIND=dp) :: LamThick(:)
    LOGICAL :: LaminateStack, CoilBody, ConstraintActive
    CHARACTER(LEN=MAX_NAME_LEN):: LaminateStackModel, CoilType
    REAL(KIND=dp) :: RotM(3,3,n)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd
    LOGICAL :: PiolaVersion, SecondOrder
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ, &
                     RotMLoc(3,3), velo(3), omega_velo(3,n), &
                     lorentz_velo(3,n), RotWJ(3)
    REAL(KIND=dp) :: LocalLamThick, skind, babs, muder, AlocR(2,nd)
    REAL(KIND=dp) :: nu_11(nd), nuim_11(nd),  &
                     nu_22(nd), nuim_22(nd),  &
                     nu_33(nd), nuim_33(nd)
    REAL(KIND=dp) :: nu_val, nuim_val
    REAL(KIND=dp) :: sigma_33(nd), sigmaim_33(nd)

    COMPLEX(KIND=dp) :: mu, C(3,3), L(3), G(3), M(3), JfixPot(n), Nu(3,3)
    COMPLEX(KIND=dp) :: LocalLamCond, JAC(nd,nd), B_ip(3), Aloc(nd), &
          CVelo(3), CVeloSum, Permittivity(nd), P_ip, DAMP(nd,nd)

    LOGICAL :: Stat, Newton, HBCurve, &
               HasVelocity, HasLorenzVelocity, HasAngularVelocity
    LOGICAL :: StrandedHomogenization, UseRotM, FoundIm

    INTEGER :: t, i, j, p, q, np, EdgeBasisDegree

    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueList_t), POINTER :: CompParams
!------------------------------------------------------------------------------
    IF (SecondOrder) THEN
       EdgeBasisDegree = 2
    ELSE
       EdgeBasisDegree = 1
    END IF

    CALL GetElementNodes( Nodes )

    MASS  = 0.0_dp
    DAMP  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    IF( Jfix ) THEN
      IF( JfixSolve ) THEN
        JfixFORCE = 0.0_dp
        JfixVec = 0.0_dp
      ELSE
        JfixPot(1:n) = CMPLX( JfixVar % Values( JfixVar % Perm( Element % NodeIndexes ) ), &
            JfixVarIm % Values( JfixVarIm % Perm( Element % NodeIndexes ) ) )
      END IF
    END IF

    JAC = 0._dp
    Newton = .FALSE.

    HasVelocity = .FALSE.
    IF(ASSOCIATED(BodyForce)) THEN
      CALL GetRealVector( BodyForce, omega_velo, 'Angular velocity', HasAngularVelocity)
      CALL GetRealVector( BodyForce, lorentz_velo, 'Lorentz velocity', HasLorenzVelocity)
      HasVelocity = HasAngularVelocity .OR. HasLorenzVelocity
    END IF
    
    CALL GetPermittivity(GetMaterial(), Permittivity, n)
    HBCurve = ListCheckPresent(Material,'H-B Curve')

    IF(HBCurve) THEN
      Newton = GetLogical( SolverParams,'Newton-Raphson iteration',Found)
      IF(.NOT. Found ) Newton = ExtNewton

      IF( GotHbCurveVar ) THEN
        CALL GetLocalSolution(AlocR(1,:), UVariable = HbCurveVar )
        Aloc = CMPLX( AlocR(1,1:nd), 0.0_dp, KIND=dp)
      ELSE
        CALL GetLocalSolution(AlocR)
        Aloc = CMPLX( AlocR(1,1:nd), AlocR(2,1:nd), KIND=dp)
      END IF
    END IF

    StrandedHomogenization = .FALSE.
    UseRotM = .FALSE.
    IF(CoilBody) THEN
      IF (CoilType == 'stranded') THEN 
        CompParams => GetComponentParams( Element )
        StrandedHomogenization = GetLogical(CompParams, 'Homogenization Model', Found)

        IF ( StrandedHomogenization ) THEN
          nu_11 = 0._dp
          nuim_11 = 0._dp
          nu_11 = GetReal(CompParams, 'nu 11', Found)
          nuim_11 = GetReal(CompParams, 'nu 11 im', FoundIm)
          IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('LocalMatrix', 'Homogenization Model nu 11 not found!')

          nu_22 = 0._dp
          nuim_22 = 0._dp
          nu_22 = GetReal(CompParams, 'nu 22', Found)
          nuim_22 = GetReal(CompParams, 'nu 22 im', FoundIm)
          IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('LocalMatrix', 'Homogenization Model nu 22 not found!')

          nu_33 = 0._dp
          nuim_33 = 0._dp
          nu_33 = GetReal(CompParams, 'nu 33', Found)
          nuim_33 = GetReal(CompParams, 'nu 33 im', FoundIm)
          IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('LocalMatrix', 'Homogenization Model nu 33 not found!')

          ! Sigma 33 is not needed in because it does not exist in stranded coil
          ! Its contribution is taken into account in the circuit module if explicit coil resistance is not used!

          UseRotM = .TRUE.
        END IF
      ELSE IF( CoilType == 'foil winding') THEN
        UseRotM = .TRUE.
      END IF
    END IF
      
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree )

    np = n*Solver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
    
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = RotWBasis, USolver = pSolver )
      
       ! Compute convection type term coming from rotation
       ! -------------------------------------------------
       IF(HasVelocity) THEN
         velo = 0.0_dp
         IF( HasAngularVelocity ) THEN
           DO i=1,n
             velo(1:3) = velo(1:3) + CrossProduct(omega_velo(1:3,i), [ &
                 basis(i) * Nodes % x(i), &
                 basis(i) * Nodes % y(i), &
                 basis(i) * Nodes % z(i)])
           END DO
         END IF
         IF( HasLorenzVelocity ) THEN
           velo(1:3) = velo(1:3) + [ &
               SUM(basis(1:n)*lorentz_velo(1,1:n)), &
               SUM(basis(1:n)*lorentz_velo(2,1:n)), &
               SUM(basis(1:n)*lorentz_velo(3,1:n))]
         END IF
       END IF

       ! Compute the conductivity tensor
       ! -------------------------------
       DO i=1,3
         DO j=1,3
           C(i,j) = SUM( Tcoef(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       P_ip = SUM( Permittivity(1:n) * Basis(1:n) )

       ! Transform the conductivity tensor (in case of a foil winding):
       ! --------------------------------------------------------------       
       IF ( UseRotM ) THEN
         DO i=1,3
           DO j=1,3
             RotMLoc(i,j) = SUM( RotM(i,j,1:n) * Basis(1:n) )             
           END DO
         END DO
         C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))
       END IF

       IF ( HBCurve ) THEN
         B_ip = MATMUL( Aloc(np+1:nd), RotWBasis(1:nd-np,:) )
         babs = MAX( SQRT(SUM(ABS(B_ip)**2)), 1.d-8 )

         IF( Newton ) THEN
           mu = ListGetFun( Material,'h-b curve',babs,dFdx=muder) / Babs
           muder = (muder-mu)/babs
         ELSE
           mu = ListGetFun( Material,'h-b curve',babs) / Babs
         END IF
       ELSE
         mu = SUM( Basis(1:n) * Acoef(1:n) )
       END IF

       IF (LaminateStack) THEN
         LocalLamThick = SUM( Basis(1:n) * LamThick(1:n) )
         LocalLamCond = SUM( Basis(1:n) * LamCond(1:n) )

         SELECT CASE(LaminateStackModel)
         CASE('low-frequency model')
           mu = mu + im*Omega* LocalLamCond* LocalLamThick **2 /12d0
         CASE('wide-frequency-band model')
           skind = SQRT(2d0*mu/(omega*LocalLamCond))
           mu = LocalLamCond * LocalLamThick * skind * omega* (1d0 + im)/8d0*&
                (-im)*SIN(im*(1d0+im)*LocalLamThick/skind)/(-im*SIN(im*(1d0+im)*LocalLamThick/skind/2d0))**2d0
!               sinh((1d0+im)*LocalLamThick/skind)/sinh((1d0+im)*LocalLamThick/skind/2d0)**2d0
         END SELECT
       END IF

       IF (HasTensorReluctivity) THEN
         IF (SIZE(Acoef_t,2) == 1) THEN
           Nu = CMPLX(0._dp, 0._dp, kind=dp)
           DO i = 1,3
             Nu(i,i) = SUM(Basis(1:n)*Acoef_t(i,1,1:n))
           END DO
         ELSE
           DO i = 1,3
             DO j = 1,3
               Nu(i,j) = SUM(Basis(1:n)*Acoef_t(i,j,1:n))
             END DO
           END DO
         END IF
       ELSE
         Nu = CMPLX(0._dp, 0._dp, kind=dp)
         Nu(1,1) = mu
         Nu(2,2) = mu
         Nu(3,3) = mu

         IF (StrandedHomogenization) THEN
           nu_val = SUM( Basis(1:n) * nu_11(1:n) ) 
           nuim_val = SUM( Basis(1:n) * nuim_11(1:n) ) 
           Nu(1,1) = CMPLX(nu_val, nuim_val, KIND=dp)
           nu_val = SUM( Basis(1:n) * nu_22(1:n) ) 
           nuim_val = SUM( Basis(1:n) * nuim_22(1:n) ) 
           Nu(2,2) = CMPLX(nu_val, nuim_val, KIND=dp)
           nu_val = SUM( Basis(1:n) * nu_33(1:n) ) 
           nuim_val = SUM( Basis(1:n) * nuim_33(1:n) ) 
           Nu(3,3) = CMPLX(nu_val, nuim_val, KIND=dp)
           Nu = MATMUL(MATMUL(RotMLoc, Nu),TRANSPOSE(RotMLoc))
         END IF
       END IF
 
       M = MATMUL( LOAD(4:6,1:n), Basis(1:n) )
       L = MATMUL( LOAD(1:3,1:n), Basis(1:n) )
      
       ! Compute C * grad(V), where C is a tensor
       ! -----------------------------------------
       L = L-MATMUL(C, MATMUL(LOAD(7,1:n), dBasisdx(1:n,:)))

       IF( Jfix ) THEN
         IF( JFixSolve ) THEN
           ! If we haven't solved for the disbalance of source terms assemble it here
           DO i = 1,n
             p = i
             JFixFORCE(p) = JFixFORCE(p) + SUM(L * dBasisdx(i,:)) * detJ * IP%s(t) 
             JFixVec(:,p) = JFixVec(:,p) + L * Basis(i) * detJ * IP%s(t)
           END DO
         ELSE         
           ! If we have already solved for the Jfix potential use it here
           L = L - MATMUL(JfixPot, dBasisdx(1:n,:))
         END IF
       END IF

       ! Compute element stiffness matrix and force vector:
       ! --------------------------------------------------

       ! If we calculate a coil, user can request that the nodal degrees of freedom are not used
       ! --------------------------------------------------------------------------------------------
       NONCOIL_CONDUCTOR: IF (ConstraintActive .AND. (SUM(ABS(C)) > AEPS .OR. ElectroDynamics) ) THEN
          !
          ! The constraint equation: -div(C*(j*omega*A+grad(V)))=0
          ! --------------------------------------------------------
          DO i=1,np
            p = i
            DO q=1,np

              ! Compute the conductivity term <C grad V,grad v> for stiffness 
              ! matrix (anisotropy taken into account)
              ! -------------------------------------------
              IF(ElectroDynamics) THEN             
                DAMP(p,q) = DAMP(p,q) + P_ip*SUM(dBasisdx(q,:)*dBasisdx(p,:))*detJ*IP % s(t)
              END IF
              STIFF(p,q) = STIFF(p,q) + SUM(MATMUL(C, dBasisdx(q,:)) * dBasisdx(p,:))*detJ*IP % s(t)
            END DO
            DO j=1,nd-np
              q = j+np
              
              ! Compute the conductivity term <j * omega * C A,grad v> for 
              ! stiffness matrix (anisotropy taken into account)
              ! -------------------------------------------
              DAMP(p,q) = DAMP(p,q) + &
                  SUM(MATMUL(C,Wbasis(j,:))*dBasisdx(i,:))*detJ*IP % s(t)

              IF(ElectroDynamics) THEN             
                MASS(p,q) = MASS(p,q) + P_ip*SUM(WBasis(j,:)*dBasisdx(i,:))*detJ*IP % s(t)
              END IF

              ! Compute the conductivity term <C grad V, eta> for 
              ! stiffness matrix (anisotropy taken into account)
              ! ------------------------------------------------
              STIFF(q,p) = STIFF(q,p) + SUM(MATMUL(C, dBasisdx(i,:))*WBasis(j,:))*detJ*IP % s(t)

              IF(ElectroDynamics) THEN             
                DAMP(q,p) = DAMP(q,p) + &
                       P_ip * SUM( WBasis(j,:)*dBasisdx(i,:) )*detJ*IP % s(t)
              END IF
            END DO
          END DO
       END IF NONCOIL_CONDUCTOR

       IF ( HasVelocity ) THEN
         DO i=1,np
           p = i
           DO j=1,nd-np
             q = j+np
             RotWJ(1:3) = RotWBasis(j,1:3)

             CVelo(1:3) = C(1:3,1)*(velo(2)*RotWJ(3) - velo(3)*RotWJ(2))
             CVelo(1:3) = CVelo(1:3) + C(1:3,2)*(-velo(1)*RotWJ(3) + velo(3)*RotWJ(1))
             CVelo(1:3) = CVelo(1:3) + C(1:3,3)*(velo(1)*RotWJ(2) - velo(2)*RotWJ(1))
             CVeloSum = REAL(0,dp)
             DO k=1,3
               CVeloSum = CVeloSum + CVelo(k)*dBasisdx(i,k)
             END DO
             STIFF(p,q) = STIFF(p,q) - CVeloSum*detJ*IP % s(t)
           END DO
         END DO
       END IF
       !
       ! j*omega*C*A + curl(1/mu*curl(A)) + C*grad(V) = 
       !        J + curl(M) - C*grad(P'):
       ! ----------------------------------------------------
       DO i = 1,nd-np
         p = i+np
         FORCE(p) = FORCE(p) + (SUM(L*WBasis(i,:)) + &
            SUM(M*RotWBasis(i,:)))*detJ*IP%s(t) 

         DO j = 1,nd-np
           q = j+np

           IF ( Newton ) THEN
             JAC(p,q) = JAC(p,q) + muder * SUM(B_ip(:)*RotWBasis(i,:)) * &
                 SUM(CONJG(B_ip(:))*RotWBasis(j,:))*detJ*IP % s(t)/Babs
           END IF

           IF( HasVelocity ) THEN
             STIFF(p,q) = STIFF(p,q) &
                 - SUM(WBasis(i,:)*MATMUL(C,CrossProduct(velo, RotWBasis(j,:))))*detJ*IP%s(t)
           END IF

           STIFF(p,q) = STIFF(p,q) + &
              SUM(MATMUL(Nu, RotWBasis(j,:))*RotWBasis(i,:))*detJ*IP%s(t)

           ! Compute the conductivity term <j * omega * C A,eta> 
           ! for stiffness matrix (anisotropy taken into account)
           ! ----------------------------------------------------
           IF (CoilType /= 'stranded') DAMP(p,q) = DAMP(p,q) + &
                SUM(MATMUL(C, WBasis(j,:))*WBasis(i,:))*detJ*IP % s(t)

           IF(ElectroDynamics ) THEN
             MASS(p,q) = MASS(p,q) + &
                P_ip*SUM(WBasis(j,:)*WBasis(i,:))*detJ*IP % s(t)
           END IF
         END DO
       END DO

     END DO

    IF ( Newton ) THEN
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + JAC
      FORCE(1:nd) = FORCE(1:nd) + MATMUL(JAC,Aloc)
    END IF

    IF(EigenSystem) THEN
      MASS(1:nd,1:nd) = MASS(1:nd,1:nd) + im*DAMP(1:nd,1:nd)
    ELSE
      STIFF(1:nd,1:nd) = -Omega**2 * MASS(1:nd,1:nd) + &
        im*Omega*DAMP(1:nd,1:nd) + STIFF(1:nd,1:nd)
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
  SUBROUTINE LocalFixMatrixC( FORCE, &
      Element, n, nd, PiolaVersion, SecondOrder )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: FORCE(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: PiolaVersion, SecondOrder
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ
    COMPLEX(KIND=dp) :: JfixPot(nd), L(3)
    LOGICAL :: Stat 
    INTEGER :: t, i, p, np, EdgeBasisDegree
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    IF (SecondOrder) THEN
      EdgeBasisDegree = 2
    ELSE
      EdgeBasisDegree = 1
    END IF

    CALL GetElementNodes( Nodes )

    FORCE = 0.0d0
    JfixPot(1:n) = CMPLX( JfixVar % Values(JfixVar % Perm(Element % NodeIndexes)), &
        JfixVarIm % Values(JfixVarIm % Perm(Element % NodeIndexes)) )
    
!    IF( SUM( ABS( JfixPot(1:n) ) ) < TINY( DetJ ) ) RETURN

    
    ! Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree )
    
    np = n*Solver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = RotWBasis, USolver = pSolver )
      
      L = MATMUL(JfixPot(1:n), dBasisdx(1:n,:))
      DO i = 1,nd-np
        p = i+np
        FORCE(p) = FORCE(p) - SUM(L*WBasis(i,:)) * detJ * IP%s(t) 
      END DO
    END DO

    CALL DefaultUpdateForce(FORCE, Element )
   
!------------------------------------------------------------------------------
  END SUBROUTINE LocalFixMatrixC
!------------------------------------------------------------------------------

  
 
  
!-----------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(  MASS, STIFF, FORCE, LOAD, Bcoef, Element, n, nd )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: LOAD(:,:), Bcoef(:)
    COMPLEX(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element, Parent, Edge
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ
    COMPLEX(KIND=dp) :: B, F, TC, L(3)
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    LOGICAL :: Stat, LineElem
    INTEGER, POINTER :: EdgeMap(:,:)
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, k, ii,jj, np, p, q, EdgeBasisDegree

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    IF (SecondOrder) THEN
       EdgeBasisDegree = 2
    ELSE
       EdgeBasisDegree = 1
    END IF

    CALL GetElementNodes( Nodes )

    MASS  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    ! We may have line elements that define BC for the conductive layers, for example.
    ! However, line elements do not have all the features of edge elements. Only
    ! certains BCs are possible. 
    LineElem = ( Element % TYPE % ElementCode / 100 <= 2 ) 
        
    ! Numerical integration:
    !-----------------------
    IF( LineElem ) THEN
      IP = GaussPoints(Element)
    ELSE
      IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
          EdgeBasisDegree=EdgeBasisDegree)
    END IF
      
    np = n*MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))
    DO t=1,IP % n
       IF( LineElem ) THEN        
         stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx )        
       ELSE 
         stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx, &
             EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 
       END IF

       B  = SUM(Basis(1:n) * Bcoef(1:n))
       L  = MATMUL(LOAD(1:3,1:n), Basis(1:n))

       F  = SUM(LOAD(4,1:n)*Basis(1:n)) !* (-im/Omega)
       TC = SUM(LOAD(5,1:n)*Basis(1:n)) !* (-im/Omega)


       ! Compute element stiffness matrix and force vector:
       !---------------------------------------------------
       DO p=1,np
         FORCE(p) = FORCE(p) + F*Basis(p)*detJ*IP % s(t)
         DO q=1,np
           STIFF(p,q) = STIFF(p,q) + TC * &
                  Basis(p)*Basis(q)*detJ*IP % s(T)
         END DO
       END DO

       ! We cannot do the following for line elements
       IF( LineElem ) CYCLE
       
       DO i = 1,nd-np
         p = i+np
         FORCE(p) = FORCE(p) - SUM(L*Wbasis(i,:))*detJ*IP%s(t)
         DO j = 1,nd-np
           q = j+np
           STIFF(p,q) = STIFF(p,q) + B * &
              SUM(Wbasis(i,:)*Wbasis(j,:))*detJ*IP%s(t)
         END DO
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixAirGapBC(MASS, STIFF, FORCE, LOAD, GapLength, AirGapMu, Element, n, nd )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: LOAD(:,:)
    COMPLEX(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element, Parent, Edge
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3), localGapLength, muAir, muVacuum
    REAL(KIND=dp) :: GapLength(:), AirGapMu(:)
    LOGICAL :: Stat
    INTEGER, POINTER :: EdgeMap(:,:)
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, np, p, q, EdgeBasisDegree

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )

    EdgeBasisDegree = 1
    IF (SecondOrder) EdgeBasisDegree = 2

    MASS  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    muVacuum = 4 * PI * 1d-7

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    np = n*MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))
    DO t=1,IP % n
      
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx, &
           EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 
     
       localGapLength  = SUM(Basis(1:n) * GapLength(1:n))
       muAir  = SUM(Basis(1:n) * AirGapMu(1:n))
 
       DO i = 1,nd-np
         p = i+np
         DO j = 1,nd-np
           q = j+np
           STIFF(p,q) = STIFF(p,q) + localGapLength / (muAir*muVacuum) * &
              SUM(RotWBasis(i,:)*RotWBasis(j,:))*detJ*IP%s(t)
         END DO
       END DO  
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixAirGapBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixThinSheet(MASS, STIFF, FORCE, LOAD, Thickness, Permeability, &
                                          Conductivity, Element, n, nd )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: LOAD(:,:)
    COMPLEX(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element, Parent, Edge
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    REAL(KIND=dp) :: Thickness(:), Permeability(:), Conductivity(:)
    REAL(KIND=dp) :: sheetThickness, mu, muVacuum, C

    COMPLEX(KIND=dp) :: DAMP(nd,nd)
    
    LOGICAL :: Stat
    INTEGER, POINTER :: EdgeMap(:,:)
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, np, p, q, EdgeBasisDegree

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )

    EdgeBasisDegree = 1
    IF (SecondOrder) EdgeBasisDegree = 2

    MASS  = 0.0_dp
    DAMP  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    muVacuum = 4 * PI * 1d-7

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    np = n*MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx, &
           EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 
       
       sheetThickness  = SUM(Basis(1:n) * Thickness(1:n))
       mu  = SUM(Basis(1:n) * Permeability(1:n))
       C  = SUM(Basis(1:n) * Conductivity(1:n))

       CONDUCTOR: IF ( ABS(C) > AEPS ) THEN
          !
          ! The constraint equation: -div(C*(j*omega*A+grad(V)))=0
          ! --------------------------------------------------------
          DO i=1,np
            p = i
            DO q=1,np

              ! Compute the conductivity term <C grad V x n,grad v x n> for stiffness 
              ! matrix (without anisotropy taken into account)
              ! -------------------------------------------
              STIFF(p,q) = STIFF(p,q) + sheetThickness * C * SUM(dBasisdx(q,:) * dBasisdx(i,:))*detJ*IP % s(t)
            END DO
            DO j=1,nd-np
              q = j+np
              
              ! Compute the conductivity term <j * omega * C A x n,grad v x n> for 
              ! stiffness matrix (without anisotropy taken into account)
              ! -------------------------------------------
              DAMP(p,q) = DAMP(p,q) + &
                  sheetThickness * C*SUM(Wbasis(j,:)*dBasisdx(i,:))*detJ*IP % s(t)

              ! Compute the conductivity term <C grad V x n, eta x n> for 
              ! stiffness matrix (without anisotropy taken into account)
              ! ------------------------------------------------
              STIFF(q,p) = STIFF(q,p) + sheetThickness * C*SUM(dBasisdx(i,:)*WBasis(j,:))*detJ*IP % s(t)
            END DO
          END DO
       END IF CONDUCTOR

       DO i = 1,nd-np
         p = i+np
         DO j = 1,nd-np
           q = j+np
           ! Magnetic energy term due to the magnetic flux density 
           ! ----------------------------------------------------
           STIFF(p,q) = STIFF(p,q) + sheetThickness / (mu*muVacuum) * &
              SUM(RotWBasis(i,:)*RotWBasis(j,:))*detJ*IP%s(t)

           ! Compute the conductivity term <j * omega * C A x n,eta x n> 
           ! for stiffness matrix (without anisotropy taken into account)
           ! ----------------------------------------------------
           IF (ABS(C) > AEPS) THEN
             DAMP(p,q) = DAMP(p,q) + sheetThickness * &
                 C * SUM(WBasis(j,:)*WBasis(i,:))*detJ*IP % s(t)
           END IF
         END DO
       END DO
    END DO

    IF(EigenSystem) THEN 
      MASS(1:nd,1:nd) = MASS(1:nd,1:nd) + im*DAMP(1:nd,1:nd)
    ELSE
      STIFF(1:nd,1:nd) = -omega**2 * MASS(1:nd,1:nd) + &
        im*Omega*DAMP(1:nd,1:nd) + STIFF(1:nd,1:nd)
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixThinSheet
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixSkinBC( MASS, STIFF, FORCE, SkinCond, SkinMu, &
                 Element, CircuitDrivenBC, n, nd )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: SkinCond(:), SkinMu(:)
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: CircuitDrivenBC
    INTEGER :: n, nd
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), DetJ
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3), cond, mu, muVacuum, delta
    LOGICAL :: Stat
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp) :: invZs, DAMP(nd,nd)
    INTEGER :: t, i, j, np, p, q, EdgeBasisDegree
    
    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )

    EdgeBasisDegree = 1
    IF (SecondOrder) EdgeBasisDegree = 2

    MASS  = 0.0_dp
    DAMP  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    muVacuum = 4 * PI * 1d-7
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    np = n*MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))

    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 

      cond = SUM(Basis(1:n) * SkinCond(1:n))
      mu  = muVacuum * SUM(Basis(1:n) * SkinMu(1:n))
      delta = SQRT( 2.0_dp/(cond*omega*mu))      
      invZs = (cond*delta)/(1.0_dp+im)
      !PRINT *,'skin:',cond,delta,omega,invZs
      !PRINT *,'elem:',Element % NodeIndexes

      !
      ! The contributions from the constraint (H x n) x n = 1/Z E x n:
      !
      DO i = 1,nd-np
        p = i+np
        DO j = 1,nd-np
          q = j+np
          !
          ! The term i*omega/Z < A x n, v x n> : With ApplyPiolaTransform = .TRUE.
          ! the edge basis functions returned by the function EdgeElementInfo are automatically 
          ! tangential and hence the normal doesn't appear in the expression.
          !
          DAMP(p,q) = DAMP(p,q) + invZs * &
              SUM(WBasis(i,:) * WBasis(j,:)) * detJ * IP % s(t)
        END DO

        IF (.NOT. CircuitDrivenBC) THEN
          DO q = 1,np
            !
            ! The term 1/Z < grad V x n, v x n> : 
            ! Some tensor calculation shows that the component form of this term is analogous to 
            ! the case < A x n, v x n>. 
            !
            STIFF(p,q) = STIFF(p,q) + invZs * &
                        SUM(WBasis(i,:) * dBasisdx(q,:)) * detJ * IP % s(t)
          END DO
        END IF

      END DO

      !
      ! The contributions from applying Ohm's law to the tangential surface current 
      !
        IF (.NOT. CircuitDrivenBC) THEN
          DO p = 1,np
            DO q = 1,np
              STIFF(p,q) = STIFF(p,q) + invZs * &
                  SUM(dBasisdx(p,:) * dBasisdx(q,:)) * detJ * IP % s(t)
            END DO

            DO j = 1,nd-np
              q = j+np
              DAMP(p,q) = DAMP(p,q) + invZs * &
                SUM(dBasisdx(p,:) * WBasis(j,:)) * detJ * IP % s(t)
            END DO
          END DO
        END IF
    END DO

    IF(EigenSystem) THEN 
      MASS(1:nd,1:nd) = MASS(1:nd,1:nd) + im*DAMP(1:nd,1:nd)
    ELSE
      STIFF(1:nd,1:nd) = -omega**2 * MASS(1:nd,1:nd) + &
        im*Omega*DAMP(1:nd,1:nd) + STIFF(1:nd,1:nd)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixSkinBC
!------------------------------------------------------------------------------

  
!-----------------------------------------------------------------------------
  FUNCTION LocalFluxBC( LOAD, Element, n, nd ) RESULT(Bn)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: LOAD(:,:), Bn
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element, Edge, Parent
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
    REAL(KIND=dp) :: Normal(3)
    COMPLEX(KIND=dp) :: L(3), ln
    LOGICAL :: Stat
    INTEGER :: t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes,  Element )
    !
    ! Integrate (B,n) over boundary face:
    ! -----------------------------------
    IP = GaussPoints(Element)
    Bn = 0._dp
    DO t=1,IP % n
      stat = ElementInfo( Element,Nodes,IP % U(t),IP % V(t), &
                 IP % W(t),detJ,Basis,dBasisdx )

      Normal=NormalVector(Element,Nodes,IP % u(t),ip % v(t),.TRUE.)
      Ln = SUM(LOAD(4,1:n)*Basis(1:n))
      L  = MATMUL(LOAD(1:3,1:n), Basis(1:n))
      Bn = Bn + Detj * IP % S(t) * (Ln+SUM(L*Normal))
    END DO
!------------------------------------------------------------------------------
  END FUNCTION LocalFluxBC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE DirichletAfromB()
!------------------------------------------------------------------------------
    USE ElementDescription, ONLY: GetEdgeMap

    IMPLICIT NONE
    REAL(KIND=dp) :: p(3),q(3),cx(3),r,xmin,ymin,zmin,xmax,ymax,zmax
    COMPLEX(KIND=dp) :: S
    TYPE(ListMatrixEntry_t), POINTER :: Ltmp
    TYPE(Matrix_t), POINTER :: Smat
    TYPE(Nodes_t),SAVE :: Nodes
    TYPE(ValueList_t), POINTER :: BC

    LOGICAL :: Found, Found1,Found2,Found3,L1,L2,L3
    INTEGER :: i,j,k,l,m,t,ii,Faces,n,nd,Active,je1,je2,pe1,pe2

    TYPE(Element_t), POINTER :: Element, Edge, Edge1
    COMPLEX(KIND=dp), ALLOCATABLE :: Bn(:)
    INTEGER, POINTER :: EdgeMap(:,:)
    INTEGER, ALLOCATABLE :: dMap(:),FaceMap(:)
    LOGICAL, ALLOCATABLE :: FluxBoundaryEdge(:), CycleEdges(:), UsedFaces(:)
!------------------------------------------------------------------------------
    ALLOCATE(FluxBoundaryEdge(Mesh % NumberOFEdges)); FluxBoundaryEdge=.FALSE.

    Active = GetNOFBoundaryElements()
    DO t=1,Active
       Element => GetBoundaryElement(t)

       IF ( GetElementFamily()==1 ) CYCLE
       BC=>GetBC()
       IF (.NOT. ASSOCIATED(BC) ) CYCLE

       Found = ListCheckPrefix(BC,'Magnetic Flux Density')

       IF ( Found ) THEN
         SELECT CASE(GetElementFamily())
         CASE(2)
           CYCLE !what would it mean in 2D,at least with only B_z solved?
         CASE(3,4)
           k = GetBoundaryFaceIndex(Element); Element => Mesh % Faces(k)
         END SELECT
         IF (.NOT. ActiveBoundaryElement(Element)) CYCLE
         FluxBoundaryEdge(Element % EdgeIndexes)=.TRUE.
       END IF
    END DO

    FluxCount = COUNT(FluxBoundaryEdge)
    IF ( FluxCount==0 ) THEN
      DEALLOCATE(FluxBoundaryEdge); RETURN
    END IF

    IF (.NOT. ALLOCATED(FluxMap) ) ALLOCATE(FluxMap(FluxCount))
    FluxCount = 0
    DO i=1,Mesh % NumberOfEdges
      IF ( FluxBoundaryEdge(i) ) THEN
        FluxCount = FluxCount+1
        FluxMap(FluxCount) = i
      END IF
    END DO
    DEALLOCATE(FluxBoundaryEdge)

    DO i=1,FluxCount
      Edge => Mesh % Edges(FluxMap(i))
      Edge % BoundaryInfo % Left => NULL()
      Edge % BoundaryInfo % Right => NULL()
    END DO

    ALLOCATE(FaceMap(Mesh % NumberOfFaces)); FaceMap=0
    Faces = 0
    DO t=1,Active
      Element => GetBoundaryElement(t)

      IF ( GetElementFamily()==1 ) CYCLE
      BC=>GetBC()
      IF (.NOT. ASSOCIATED(BC) ) CYCLE

      Found = ListCheckPrefix(BC,'Magnetic Flux Density')
      IF ( .NOT. Found ) CYCLE

      k = GetBoundaryFaceIndex(Element); Element=>Mesh % Faces(k)
      IF (.NOT. ActiveBoundaryElement(Element)) CYCLE
      Faces = Faces + 1
      FaceMap(k) = Faces

      DO i=1,Element % TYPE % NumberOfNodes
        Edge => Mesh % Edges(Element % EdgeIndexes(i))
        IF (.NOT.ASSOCIATED(Edge % BoundaryInfo % Left)) THEN
           Edge % BoundaryInfo % Left => Element
        ELSE IF (.NOT.ASSOCIATED(Edge % BoundaryInfo % Right)) THEN
           Edge % BoundaryInfo % Right => Element
        END IF
      END DO
    END DO


    ! Make gauge tree for the boundary:
    ! ---------------------------------
    CALL GaugeTreeFluxBC(Solver,Mesh,TreeEdges,BasicCycles,FluxCount,FluxMap)

    WRITE(Message,*) 'Boundary tree edges: ', &
      i2s(COUNT(TreeEdges(FluxMap))),   &
             ' of total: ',i2s(FluxCount)
    CALL Info('WhitneyAVHarmonicSolver: ', Message, Level=5)

    ! Get (B,n) for BC faces:
    ! -----------------------
    ALLOCATE(Bn(Faces))
    DO t=1,Active
      Element => GetBoundaryElement(t)

      IF ( GetElementFamily()==1 ) CYCLE
      BC=>GetBC()
      IF (.NOT. ASSOCIATED(BC) ) CYCLE

      n  = GetElementNOFNodes(Element)
      CALL GetComplexVector(BC,Load(1:3,1:n),'Magnetic Flux Density',Found1)

      LOAD(4,1:n) = GetReal(BC,'Magnetic Flux Density {n}',Found)
      LOAD(4,1:n) = LOAD(4,1:n)+im*GetReal(BC,'Magnetic Flux Density im {n}',L1)
      Found = Found.OR.L1

      IF (Found.OR.Found1) THEN
        k = GetBoundaryFaceIndex(Element)
        Element => Mesh % Faces(k)
        IF (.NOT.ActiveBoundaryElement(Element)) CYCLE        
        nd = GetElementNOFDOFs(Element)
        Bn(FaceMap(k))=LocalFluxBC(LOAD,Element,n,nd)
      END IF
    END DO

    !
    ! Calculate value for free edges using the Fundamental Loop Basis
    ! generated by GaugeTreeFluxBC():
    ! ---------------------------------------------------------------
    ALLOCATE(CycleEdges(Mesh % NumberOFEdges), UsedFaces(Faces))
    CycleEdges = .FALSE.
    ALLOCATE(dMap(MAXVAL(BasicCycles(:) % Degree)))

    Smat => GetMatrix()
    DO i=1,SIZE(BasicCycles)
      IF (BasicCycles(i) % Degree<=0 ) CYCLE

      ! 
      ! Extract loop edge indices: 
      ! --------------------------
      j = 0
      Ltmp => BasicCycles(i) % Head
      DO WHILE(ASSOCIATED(Ltmp))
        j = j + 1
        dMap(j) = Ltmp % Index; Ltmp => Ltmp % Next
      END DO
      IF ( j<= 0 ) CYCLE

      !
      ! Orient edges to form a polygonal path:
      ! --------------------------------------
      Edge  => Mesh % Edges(dMap(j))
      Edge1 => Mesh % Edges(dMap(j-1))
      IF ( ANY(Edge % NodeIndexes(1)==Edge1 % NodeIndexes) ) THEN
        l = Edge % NodeIndexes(1)
        Edge % NodeIndexes(1) = Edge % NodeIndexes(2)
        Edge % NodeIndexes(2) = l
      END IF
 
      DO k=j-1,1,-1
        Edge1 => Mesh % Edges(dMap(k))
        IF (Edge % NodeIndexes(2)==Edge1 % NodeIndexes(2)) THEN
          l = Edge1 % NodeIndexes(1)
          Edge1 % NodeIndexes(1) = Edge1 % NodeIndexes(2)
          Edge1 % NodeIndexes(2) = l
        END IF
        Edge => Edge1
      END DO

      !
      ! Try to find which way is inside...
      ! ----------------------------------
      Edge => Mesh % Edges(dMap(j))
      Element => Edge % BoundaryInfo % Left
      IF ( j==3 ) THEN
        m = 0
        DO k=1,3
          DO l=1,3
            IF (dMap(l)==Element % EdgeIndexes(k)) m=m+1
          END DO
        END DO
        L1 = m==3
        IF ( .NOT. L1 ) Element=>Edge % BoundaryInfo % Right
        S = Bn(FaceMap(Element % ElementIndex))
      ELSE
        ! If not a triangle, try a (planar) polygonal test. This
        ! will fail for general 3D paths. We'll spot the failure
        ! later by trial and error...Might be preferable to skip
        ! this altogether? Dunno....
        ! ------------------------------------------------------
        xmin=HUGE(xmin); xmax=-HUGE(xmax);
        ymin=HUGE(ymin); ymax=-HUGE(ymax);
        zmin=HUGE(zmin); zmax=-HUGE(zmax);
        DO k=1,j
          Edge1 => Mesh % Edges(dMap(k))
          DO l=1,2
            m = Edge1 % NodeIndexes(l)
            xmin = MIN(xmin,Mesh % Nodes % x(m))
            ymin = MIN(ymin,Mesh % Nodes % y(m))
            zmin = MIN(zmin,Mesh % Nodes % z(m))

            xmax = MAX(xmax,Mesh % Nodes % x(m))
            ymax = MAX(ymax,Mesh % Nodes % y(m))
            zmax = MAX(zmax,Mesh % Nodes % z(m))
          END DO
        END DO
        L1 = xmax-xmin > ymax-ymin
        L2 = xmax-xmin > zmax-zmin
        L3 = ymax-ymin > zmax-zmin
        IF ( l1 ) THEN
          l=1
          IF ( l3 ) THEN
            m=2; n=3
          ELSE
            m=3; n=2
          END IF
        ELSE
          IF ( l2 ) THEN
            l=1; m=2; n=3
          ELSE
            l=3; m=1; n=2
          END IF
        END IF
        cx(l) = SUM(Mesh % Nodes % x(Element % NodeIndexes))/3._dp
        cx(m) = SUM(Mesh % Nodes % y(Element % NodeIndexes))/3._dp
        cx(n) = SUM(Mesh % Nodes % z(Element % NodeIndexes))/3._dp

        L1 = .FALSE.
        DO k=j,1,-1
          Edge1 => Mesh % Edges(dMap(k))
          je1 = Edge1 % NodeIndexes(1)
          je2 = Edge1 % NodeIndexes(2)
          p(l) = Mesh % Nodes % x(je1)
          p(m) = Mesh % Nodes % y(je1)
          p(n) = Mesh % Nodes % z(je1)

          q(l) = Mesh % Nodes % x(je2)
          q(m) = Mesh % Nodes % y(je2)
          q(n) = Mesh % Nodes % z(je2)

          IF ((q(2)>cx(2)).NEQV.(p(2)>cx(2))) THEN
            IF (cx(1)<(p(1)-q(1))*(cx(2)-q(2))/(p(2)-q(2))+q(1)) L1=.NOT.L1
          END IF
        END DO
        IF (.NOT.L1) THEN
          IF (ASSOCIATED(Edge % BoundaryInfo % Right)) &
            Element=>Edge % BoundaryInfo % Right
        END IF

        ! Compute integral of (B,n) inside the cycle path
        ! -----------------------------------------------
        CycleEdges(dMap(1:j))=.TRUE.
        DO m=1,2
          S=0; UsedFaces=.FALSE.;
          IF( FloodFill(Element,CycleEdges, &
                       FaceMap,UsedFaces,Bn,S) )EXIT

          ! the in/out guess was wrong, try the other way:
          ! ----------------------------------------------
          IF (ASSOCIATED(Edge % BoundaryInfo % Right,Element)) THEN
            Element => Edge % BoundaryInfo % Left
          ELSE
            Element => Edge % BoundaryInfo % Right
          END IF
        END DO
        CycleEdges(dMap(1:j))=.FALSE.
      END IF

      !
      ! Orient edge to parent triangle...
      ! ---------------------------------
      je1 = Edge % NodeIndexes(1)
      je2 = Edge % NodeIndexes(2)
      EdgeMap => GetEdgeMap(GetElementFamily(Element))
      DO t=1,Element % TYPE % NumberOfEdges
        pe1 = Element % NodeIndexes(EdgeMap(t,1))
        pe2 = Element % NodeIndexes(EdgeMap(t,2))
        IF (pe1==je1.AND.pe2==je2 .OR. pe1==je2.AND.pe2==je1) EXIT
      END DO
      IF ( pe1/=je1 ) S=-S

      !
      ! ...because we now know how to orient against outward normal:
      ! ------------------------------------------------------------
      CALL GetElementNodes(Nodes,Element)
      p = NormalVector(Element,Nodes,0._dp,0._dp)
      q = NormalVector(Element,Nodes,0._dp,0._dp,.TRUE.)
      IF ( SUM(p*q)<0 ) S=-S

      !
      ! Check whether some edges in the path have nonzero values,
      ! if so, substract from integral:
      ! ---------------------------------------------------------
      DO k=j-1,1,-1
        l = Perm(dMap(k)+nNodes)
        R = Smat % RHS(l) / Smat % Values(Smat % Diag(l))
        IF ( R==0 ) CYCLE

        Edge1 => Mesh % Edges(dMap(k))
        pe1=Edge1 % NodeIndexes(1)
        pe2=Edge1 % NodeIndexes(2)
        IF ( pe2<pe1 ) R=-R; S=S-R
      END DO

      !
      ! ...and finally we should have the edge value:
      ! ---------------------------------------------
      IF ( je2<je1 ) S=-S
      CALL SetDOFtoValue(Solver,dMap(j),S)
    END DO
    DEALLOCATE(dMap, CycleEdges, FaceMap, UsedFaces, Bn)
    CALL List_FreeMatrix(SIZE(BasicCycles), BasicCycles)
!------------------------------------------------------------------------------
  END SUBROUTINE DirichletAfromB 
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  RECURSIVE FUNCTION FloodFill(Element,CycleEdges, &
          FaceMap,UsedFaces,Bn,CycleSum) RESULT(Found)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: e, Element
    COMPLEX(KIND=dp) :: CycleSum, Bn(:)
    INTEGER :: i,j,n, FaceMap(:)
    LOGICAL :: CycleEdges(:), UsedFaces(:), Found

    Found=.FALSE.
    IF (.NOT.ASSOCIATED(Element)) RETURN

    n=FaceMap(Element % ElementIndex)
    IF (UsedFaces(n)) THEN
      Found=.TRUE.; RETURN
    END IF
    UsedFaces(n)=.TRUE.
    CycleSum = CycleSum+Bn(n)

    DO i=1,Element % TYPE % NumberOfEdges
      j = Element % EdgeIndexes(i)
      IF ( CycleEdges(j) ) CYCLE

      e => Mesh % Edges(j) % BoundaryInfo % Right
      IF(.NOT.FloodFill(e,CycleEdges,FaceMap,UsedFaces,Bn,CycleSum)) RETURN

      e => Mesh % Edges(j) % BoundaryInfo % Left
      IF(.NOT.FloodFill(e,CycleEdges,FaceMap,UsedFaces,Bn,CycleSum)) RETURN
    END DO
    Found=.TRUE.; RETURN
!------------------------------------------------------------------------------
  END FUNCTION FloodFill
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE WhitneyAVHarmonicSolver
!------------------------------------------------------------------------------

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
! *
! *  Utilities written as solvers to compute the Helmholtz projection P(A)
! *  of a curl-conforming vector field A. The projection can be obtained as 
! *  P(A) = A - W where  W is the curl-conforming field fitted to represent 
! *  grad Phi, with Phi being a H1-regular scalar field.
! * 
! *  This file contains harmonic version of the transformation and also applies the
! *  correction to the V field within conducting regions.
! *
! *
! *  Authors: Mika Malinen, Juha Ruokolainen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: March 20, 2020
! *  Last Modified: June 18, 2021, Juha
! *
!******************************************************************************


!------------------------------------------------------------------------------
SUBROUTINE HelmholtzProjector_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  LOGICAL :: Found
  INTEGER :: i
  TYPE(ValueList_t), POINTER :: SolverParams
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  DO i=1,Model % NumberOfSolvers
    IF (ListGetLogical( Model % Solvers(i) % Values, 'Helmholtz Projection', Found)) EXIT
  END DO
  CALL ListCopyPrefixedKeywords(Model % Solvers(i) % Values, SolverParams, 'HelmholtzProjector:')
!------------------------------------------------------------------------------
END SUBROUTINE HelmholtzProjector_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE HelmholtzProjector_Init(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  LOGICAL :: Found
  INTEGER :: i,j
  TYPE(ValueList_t), POINTER :: SolverParams
!------------------------------------------------------------------------------

  SolverParams => GetSolverParams()
  CALL ListAddNewLogical(SolverParams, 'Linear System Refactorize', .FALSE.)

  CALL ListAddString( SolverParams, 'Variable', 'pd' )
  CALL ListAddLogical( SolverParams, 'Variable Output',.FALSE. )

  DO i=1,Model % NumberOfSolvers
    IF(ListGetLogical( Model % Solvers(i) % Values, 'Helmholtz Projection', Found)) EXIT
  END DO
  CALL ListAddString( SolverParams, 'Potential Variable', GetVarName(Model % Solvers(i) % Variable))

  ! Solver is using a single linear system to solve complex components,
  ! assign storage for final complex result.
  ! -------------------------------------------------------------------
  CALL ListAddString( SolverParams, 'Exported Variable 1', 'P[P re:1 P im:1]' )

  CALL ListAddLogical( SolverParams, 'Linear System Symmetric', .TRUE. )
  CALL ListAddString(  SolverParams, 'Linear System Solver', 'Iterative' )
  CALL ListAddString(  SolverParams, 'Linear System Preconditioning', 'ILU' )
  CALL ListAddInteger( SolverParams, 'Linear System Residual Output', 25 )
  CALL ListAddInteger( SolverParams, 'Linear System Max Iterations', 2000 )
  CALL ListAddString(  SolverParams, 'Linear System Iterative Method', 'CG' )
  CALL ListAddConstReal(   SolverParams, 'Linear System Convergence Tolerance', 1.0d-9 )

  DO j=1,Model % NumberOfBCs
    IF ( ListCheckPrefix( Model % BCs(j) % Values, &
               TRIM(GetVarName(Model % Solvers(i) % Variable)) // ' re {e}' ) ) THEN
      CALL ListAddConstReal( Model % BCs(j) % Values, 'Pd', 0._dp )
    END IF
  END DO
!------------------------------------------------------------------------------
END SUBROUTINE HelmholtzProjector_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute a H1-regular scalar field to obtain the Helmholtz projection P(A)
!> of a curl-conforming vector field A. Given the solution field Phi of this 
!> solver, the projection can be evaluated as P(A) = A - grad Phi.
!------------------------------------------------------------------------------
SUBROUTINE HelmholtzProjector(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: SolverPtr
  TYPE(Element_t), POINTER :: Element

  LOGICAL :: AllocationsDone = .FALSE.
  LOGICAL :: Found
  LOGICAL :: PiolaVersion, SecondOrder
  LOGICAL :: ConstantBulkMatrix

  INTEGER :: i, j,k,l,n, n_pot, nd_pot, t
  INTEGER :: dim, PotDOFs
  INTEGER :: istat, active

  REAL(KIND=dp), ALLOCATABLE, TARGET :: Stiff(:,:), Force(:,:), PotSol(:,:), F(:,:)
  REAL(KIND=dp) :: Norm, Omega
  REAL(KIND=dp), POINTER :: SaveRHS(:), SOL(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: PotName

  TYPE(Variable_t), POINTER :: v

  SAVE Stiff, Force, PotSol, AllocationsDone
!------------------------------------------------------------------------------
  CALL DefaultStart()

  dim = CoordinateSystemDimension()
  SolverParams => GetSolverParams()
  Mesh => GetMesh()

  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  IF (.NOT. AllocationsDone) THEN
    n = Mesh % MaxElementDOFs

    ALLOCATE( FORCE(2,n), STIFF(n,n), PotSol(2,n), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'HelmholtzProjectorZ', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
  END IF

  !
  ! Find the variable which is projected:
  !
  PotName = GetString(SolverParams, 'Potential Variable', Found)
  IF (.NOT. Found ) PotName = 'av'
  Found = .FALSE.
  DO i=1,Model % NumberOfSolvers
    SolverPtr => Model % Solvers(i)
    IF (PotName == GetVarName(SolverPtr % Variable)) THEN
      Found = .TRUE.
      EXIT
    END IF
  END DO

  IF (Found ) THEN
    CALL Info('HelmholtzProjectorZ', 'Solver inherits potential '&
         //TRIM(PotName)//' from solver: '//I2S(i),Level=7)
  ELSE
    CALL Fatal('HelmholtzProjectorZ', 'Solver associated with potential variable > '&
        //TRIM(PotName)//' < not found!')
  END IF

  !
  ! Find some parameters to inherit the vector FE basis as defined in 
  ! the primary solver:
  !
  CALL EdgeElementStyle(SolverPtr % Values, PiolaVersion, QuadraticApproximation = SecondOrder )
  IF (PiolaVersion) CALL Info('HelmholtzProjectorZ', &
      'Using Piola-transformed finite elements', Level=5)

  n = Solver % Matrix % NumberOfRows
  ALLOCATE(F(n,2)); F=0

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize(Solver)

  SaveRHS => Solver % Matrix % RHS

  DO t=1,active
    Element => GetActiveElement(t)
    !
    ! This solver relies on getting basis functions by calling a routine
    ! which returns a curl-conforming basis. It is thus assumed that
    ! the background mesh defines the number of Lagrange basis functions.
    !
    n = GetElementNOFNodes()
   
    ! The DOF counts for the potential (target) variable: 
    n_pot = n*SolverPtr % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)
    nd_pot = GetElementNOFDOFs(USolver=SolverPtr)

    CALL GetVectorLocalSolution(PotSol, PotName, USolver=SolverPtr)

    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(Stiff, Force, Element, n, dim, PiolaVersion, &
        SecondOrder, n_pot, nd_pot, PotSol )
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    Solver % Matrix % RHS => F(:,1)
    CALL DefaultUpdateForce(FORCE(1,:))

    Solver % Matrix % RHS => F(:,2)
    CALL DefaultUpdateEquations(STIFF, FORCE(2,:))
  END DO

  CALL DefaultFinishBulkAssembly()
  CALL DefaultDirichletBCs()

  v => VariableGet( Mesh % Variables, 'P'  )
  SOL => v % Values

  Solver % Matrix % RHS => F(:,1)
  CALL DefaultDirichletBCs()
  Norm = DefaultSolve()  
  SOL(1::2) = Solver % Variable % Values

  Solver % Matrix % RHS => F(:,2)
  CALL DefaultDirichletBCs()
  Norm = DefaultSolve()  
  SOL(2::2) = Solver % Variable % Values

  Solver % Matrix % RHS => SaveRHS


  omega = GetAngularFrequency()

  !
  ! Finally, redefine the potential variable:
  ! -----------------------------------------
  DO i=1,Solver % Mesh % NumberOfNodes
    j = Solver % Variable % Perm(i)
    IF(j==0) CYCLE

    k = SolverPtr % Variable % Perm(i)
    IF (k == 0) THEN
      CALL Fatal('HelmholtzProjectorZ', &
        'The variable and potential permutations are nonmatching?')
    END IF

    SolverPtr % Variable % Values(2*k-1) = SolverPtr % Variable % Values(2*k-1) - &
        omega * SOL(2*j)

    SolverPtr % Variable % Values(2*k) = SolverPtr % Variable % Values(2*k) + &
        omega * SOL(2*j-1)
  END DO

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Stiff, Force, Element, n, dim, PiolaVersion, &
               SecondOrder, n_pot, nd_pot, PotSol )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Stiff(:,:), Force(:,:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n   ! The number of background element nodes
    INTEGER :: dim
    LOGICAL :: PiolaVersion, SecondOrder
    INTEGER :: n_pot, nd_pot      ! The size parameters of target field
    REAL(KIND=dp) :: PotSol(:,:)  ! The values of target field DOFS
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    LOGICAL :: Stat

    INTEGER :: i, j, p, q, t, EdgeBasisDegree 

    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), A(2,3)
    REAL(KIND=dp) :: u, v, w, s, DetJ
    REAL(KIND=dp) :: WBasis(nd_pot-n_pot,3), CurlWBasis(nd_pot-n_pot,3)
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes)

    STIFF = 0.0_dp
    FORCE = 0.0_dp

    IF (SecondOrder) THEN
      EdgeBasisDegree = 2
    ELSE 
      EdgeBasisDegree = 1
    END IF       
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree=EdgeBasisDegree)
    IF( dim == 2 .AND. .NOT. PiolaVersion ) THEN
      CALL Fatal('HelmholtzProjectorZ', '"Use Piola Transform = True" needed in 2D')
    END IF
        
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = CurlWBasis, USolver = SolverPtr ) 
      s = detJ * IP % s(t)

      A = MATMUL(PotSol(:,n_pot+1:nd_pot), WBasis(1:nd_pot-n_pot,:))

      DO p=1,n
        DO q=1,n
          STIFF(p,q) = STIFF(p,q) + SUM(dBasisdx(q,1:dim) * dBasisdx(p,1:dim)) * s
        END DO
      END DO

      DO p=1,n
        FORCE(1,p) = FORCE(1,p) + SUM(A(1,:) * dBasisdx(p,:)) * s
        FORCE(2,p) = FORCE(2,p) + SUM(A(2,:) * dBasisdx(p,:)) * s
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE HelmholtzProjector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE RemoveKernelComponent_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  INTEGER :: i,j
  CHARACTER(LEN=MAX_NAME_LEN) :: Avname
  LOGICAL :: Found, PiolaVersion, SecondOrder, SecondFamily
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  CALL ListAddLogical(SolverParams, 'Linear System Refactorize', .FALSE.)


! ! Solver is using a single linear system to solve complex components,
! ! the final result is the (Hcurl) imaginary component...
! ! -------------------------------------------------------------------

! Linear System Symmetric = True
! Linear System Solver = "Iterative"
! Linear System Preconditioning = None
! Linear System Residual Output = 25
! Linear System Max Iterations = 2000
! Linear System Iterative Method = CG
! Linear System Convergence Tolerance = 1.0e-9

  CALL ListAddString( SolverParams, 'Variable', 'avm' )
  CALL ListAddLogical( SolverParams, 'Variable Output',.FALSE. )
  
  DO i=1,Model % NumberOfSolvers
    IF(ListGetLogical( Model % Solvers(i) % Values, 'Helmholtz Projection', Found)) EXIT
  END DO
  AVname = ListGetString( Model % Solvers(i) % Values, 'Variable'  )
  
  j = index(AVname, '[')
  IF(j>0) AVname = AVname(1:j-1)
  CALL ListAddString( GetSolverParams(), 'Potential Variable', AVName )

  IF (.NOT. ListCheckPresent(SolverParams, "Element")) THEN   
    CALL EdgeElementStyle(Model % Solvers(i) % Values, PiolaVersion, SecondFamily, SecondOrder )

    IF (SecondOrder) THEN
      CALL ListAddString(SolverParams, "Element", &
          "n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2")
    ELSE IF (SecondFamily) THEN
      CALL ListAddString(SolverParams, "Element", "n:0 e:2")
    ELSE IF (PiolaVersion) THEN
      CALL ListAddString(SolverParams, "Element", &
          "n:0 e:1 -brick b:3 -quad_face b:2")
    ELSE
      CALL ListAddString( SolverParams, "Element", "n:0 e:1")
    END IF
  END IF

  ! Solver is using a single linear system to solve complex components,
  ! assign storage for final complex result.
  ! -------------------------------------------------------------------
  CALL ListAddString( SolverParams, 'Kernel Variable', 'P' )

  CALL ListAddLogical( SolverParams, 'Linear System Symmetric', .TRUE. )
  CALL ListAddString(  SolverParams, 'Linear System Solver', 'Iterative' )
  CALL ListAddString(  SolverParams, 'Linear System Preconditioning', 'ILU' )
  CALL ListAddInteger( SolverParams, 'Linear System Residual Output', 25 )
  CALL ListAddInteger( SolverParams, 'Linear System Max Iterations', 2000 )
  CALL ListAddString(  SolverParams, 'Linear System Iterative Method', 'CG' )
  CALL ListAddConstReal( SolverParams, 'Linear System Convergence Tolerance', 1.0d-9 )
  CALL ListAddLogical( SolverParams,"Hcurl Basis",.TRUE.)

  CALL ListCopyPrefixedKeywords(Model % Solvers(i) % Values, SolverParams, 'RemoveKernelComponent:')
  
!------------------------------------------------------------------------------
END SUBROUTINE RemoveKernelComponent_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  Apply the Helmholtz projection on a curl-conforming vector field A
!>  when the kernel component grad phi of A (with respect to the curl operator)
!>  has been computed by using the subroutine HelmholtzProjector. This solver
!>  generates the representation W of grad phi in terms of the curl-conforming
!>  basis and finally redefines A := A - W, with W = grad phi. 
!------------------------------------------------------------------------------
SUBROUTINE RemoveKernelComponent(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: SolverPtr, KerSolverPtr
  TYPE(Element_t), POINTER :: Element

  LOGICAL :: AllocationsDone = .FALSE.
  LOGICAL :: Found
! LOGICAL :: SecondFamily
  LOGICAL :: PiolaVersion, SecondOrder
  LOGICAL :: ConstantBulkMatrix

  INTEGER :: dim, PotDOFs
  INTEGER :: i, j, k, n, nd, n_pot, nd_pot, t
  INTEGER :: istat, active

  REAL(KIND=dp), ALLOCATABLE, TARGET :: Stiff(:,:), Force(:,:), PhiSol(:,:), &
                    SOL(:), F(:,:)
  REAL(KIND=dp) :: Norm
  CHARACTER(LEN=MAX_NAME_LEN) :: PotName, Name

  REAL(KIND=dp), POINTER :: SaveRHS(:)


  TYPE(Variable_t), POINTER :: v

  SAVE Stiff, Force, PhiSol, AllocationsDone
!------------------------------------------------------------------------------
  CALL DefaultStart()

  dim = CoordinateSystemDimension()
  SolverParams => GetSolverParams()
  Mesh => GetMesh()

  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  IF (.NOT. AllocationsDone) THEN
    n = Mesh % MaxElementDOFs
    ALLOCATE(FORCE(2,n),STIFF(n,n),PhiSol(2,n),STAT=istat)
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'RemoveKernelComponent', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
  END IF

  !
  ! Find the variable which is projected:
  !
  PotName = GetString(SolverParams, 'Potential Variable', Found)
  IF (.NOT. Found ) PotName = 'av'

  Found = .FALSE.
  DO i=1,Model % NumberOfSolvers
    SolverPtr => Model % Solvers(i)
    IF (PotName == GetVarName(SolverPtr % Variable)) THEN
      Found = .TRUE.
      EXIT
    END IF
  END DO

  IF (.NOT. Found ) THEN
    CALL Fatal('RemoveKernelComponent', 'Solver associated with potential variable > '&
        //TRIM(PotName)//' < not found!')
  END IF

  PotDOFs = SolverPtr % Variable % DOFs
  IF (PotDOFs < 2) CALL Fatal('RemoveKernelComponent', 'A complex-valued potential expected')

  !
  ! Find the variable which defines the kernel component:
  !
  Name = GetString(SolverParams, 'Kernel Variable', Found)
  IF (.NOT. Found ) Name = 'phi'
  V => VariableGet( Mesh % Variables, Name )

  Found = ASSOCIATED(v)
   
  IF (.NOT. Found ) THEN
    CALL Fatal('RemoveKernelComponent', 'Solver associated with kernel variable > '&
        //TRIM(Name)//' < not found!')
  END IF
  
  !
  ! Find some parameters to inherit the vector FE basis as defined in the primary solver:
  !
  CALL EdgeElementStyle(SolverPtr % Values, PiolaVersion, QuadraticApproximation = SecondOrder )

  IF (PiolaVersion) CALL Info('RemoveKernelComponent', &
      'Using Piola-transformed finite elements', Level=5)

!  SecondFamily = GetLogical(SolverPtr % Values, 'Second Kind Basis', Found)

  n = Solver % Matrix % NumberOfRows
  ALLOCATE(F(n,2)); F=0
  SaveRHS => Solver % Matrix % RHS

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize(Solver)

  DO t=1,active
    Element => GetActiveElement(t)

    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
   
    ! The DOF counts for the potential variable: 
    n_pot = n*SolverPtr % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)
    nd_pot = GetElementNOFDOFs(USolver=SolverPtr)

    IF (nd /= nd_pot-n_pot) CALL Fatal('RemoveKernelComponent', &
     'Potential variable DOFs count /= the solver DOFs count')

    CALL GetLocalSolution(PhiSol, Name)

    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix( STIFF, FORCE, Element, n, nd, dim, PiolaVersion, &
                SecondOrder, PhiSol )
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    Solver % Matrix % RHS => F(:,1)
    CALL DefaultUpdateForce(FORCE(1,:))

    Solver % Matrix % RHS => F(:,2)
    CALL DefaultUpdateEquations(STIFF, FORCE(2,:))
  END DO

  n = Solver % Matrix % NumberOfRows
  ALLOCATE(SOL(2*n))

  Solver % Matrix % RHS => F(:,1)
  CALL DefaultDirichletBCs()
  Norm = DefaultSolve()  
  SOL(1::2) = Solver % Variable % Values

  Solver % Matrix % RHS => F(:,2)
  CALL DefaultDirichletBCs()
  Norm = DefaultSolve()  
  SOL(2::2) = Solver % Variable % Values

  Solver % Matrix % RHS => SaveRHS

  !
  ! Finally, redefine the potential variable:
  !
  n = SIZE(Solver % Variable % Perm(:))
  IF (n ==  SIZE(SolverPtr % Variable % Perm(:))) THEN
    DO i=Solver % Mesh % NumberOfNodes+1,n
      j = Solver % Variable % Perm(i)
      IF (j<=0) CYCLE

      k = SolverPtr % Variable % Perm(i)
      IF (k<=0) THEN
        CALL Fatal('RemoveKernelComponent', &
          'The variable and potential permutations are nonmatching?')
      END IF

      SolverPtr % Variable % Values(2*k-1) = SolverPtr % Variable % Values(2*k-1) - &
          SOL(2*j-1)

      SolverPtr % Variable % Values(2*k) = SolverPtr % Variable % Values(2*k) - &
          SOL(2*j)
    END DO
  ELSE
    CALL Fatal('RemoveKernelComponent', 'The variable and potential permutations differ')  
  END IF

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Stiff, Force, Element, n, nd, dim, PiolaVersion, &
              SecondOrder, PhiSol )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:,:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, dim
    REAL(KIND=dp) :: PhiSol(:,:)
    LOGICAL :: PiolaVersion, SecondOrder
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    LOGICAL :: Stat

    INTEGER :: i, j, p, q, t, EdgeBasisDegree 

    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), A(2,3)
    REAL(KIND=dp) :: s, DetJ
    REAL(KIND=dp) :: WBasis(nd,3), CurlWBasis(nd,3)
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes)

    STIFF = 0.0_dp
    FORCE = 0.0_dp

    IF (SecondOrder) THEN
      EdgeBasisDegree = 2
    ELSE 
      EdgeBasisDegree = 1
    END IF      

    IF (dim==2 .AND. .NOT. PiolaVersion) THEN
      CALL Fatal('RemoveKernelComponent', '"Use Piola Transform = True" needed in 2D')
    END IF

    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree=EdgeBasisDegree)

    DO t=1,IP % n      
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = CurlWBasis, USolver = SolverPtr )

      s = detJ * IP % s(t)
      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + s * SUM(WBasis(q,:) * WBasis(p,:))
        END DO
      END DO

      A = MATMUL( PhiSol(:,1:n), dBasisdx(1:n,:) )
      DO q=1,nd
        FORCE(1,q) = FORCE(1,q) + s * SUM(A(1,:) * WBasis(q,:))
        FORCE(2,q) = FORCE(2,q) + s * SUM(A(2,:) * WBasis(q,:))
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE RemoveKernelComponent
!------------------------------------------------------------------------------

