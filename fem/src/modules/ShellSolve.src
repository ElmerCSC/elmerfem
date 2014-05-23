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
! *     ======================================================================
! *     THE REISSNER-MINDLIN FACET SHELL SOLVER MODULE FOR COMPOSITE LAMINATES
! *     ======================================================================
! *
! *     Author: Mikko Lyly
! *     Address: CSC - IT Center for Science Ltd.
! *     Keilaranta 14, P.O. BOX 405
! *     02101 Espoo, Finland
! *     EMail: Mikko.Lyly@csc.fi
! *
! *     Date: 09 Apr 2000
! *     Modified by: Mikko Lyly & Petri Kere
! *     Date of modification: 29.5.2001, 21.8.2001
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Initialization for the primary solver: ShellSolver
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE ShellSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
    USE DefUtils

    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: DT
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: SolverParams
!------------------------------------------------------------------------------
    SolverParams => GetSolverParams()
    IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      CALL ListAddInteger( SolverParams, 'Variable DOFs', 6 )
      CALL ListAddString( SolverParams, 'Variable', 'Deflection' )
    END IF
    CALL ListAddInteger( SolverParams, 'Time derivative order', 2 )
!------------------------------------------------------------------------------
  END SUBROUTINE ShellSolver_Init
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
!>  Solve the Reissner-Mindlin facet shell equations!
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE ShellSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET :: Solver
 
     REAL(KIND=DP) :: dt
     LOGICAL :: TransientSimulation
 
     TYPE(Solver_t), POINTER :: PSolver
     TYPE(Matrix_t),POINTER  :: StiffMatrix
     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement

     INTEGER :: i,j,k,n,t, bf_id, istat, LocalNodes, CalcSurf, nPL, iPL, &
          NumberOfElementNodes

     REAL(KIND=DP) :: Norm,PrevNorm, PrevUNorm, Unorm, NonLinConvTol, &
         RelChange, RelUChange, NOFEigValsBackup, LoadScale
     INTEGER :: NonLinMaxIt, NonLinIter

     LOGICAL :: AllocationsDone = .FALSE., LargeDeflection = .FALSE., &
         StabilityAnalysis = .FALSE., StressComputation = .FALSE.

     INTEGER, POINTER :: NodeIndexes(:), DeflectionPerm(:)
  
     REAL(KIND=dp), POINTER :: PlyEng(:,:), PointLoad(:,:)

     REAL(KIND=dp), POINTER :: Deflection(:), ForceVector(:), &
          SxxNodal(:), SyyNodal(:), SzzNodal(:), &
          SxyNodal(:), SxzNodal(:), SyzNodal(:), &
          EpsxxNodal(:), EpsyyNodal(:), EpszzNodal(:), &
          EpsxyNodal(:), EpsxzNodal(:), EpsyzNodal(:)

     REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), &
          LoadX(:), LoadY(:), LoadZ(:), LoadN(:), FORCE(:), &
          Poisson(:), Thickness(:), Young(:), Tension(:), &
          MASS(:,:), DAMP(:,:), Density(:), &
          LocalDeflection(:), SxxElement(:), SyyElement(:), SzzElement(:), &
          SxyElement(:), SxzElement(:), SyzElement(:), LoadVector(:,:), &
          Weights(:,:), Referenced(:), &
          EpsxxElement(:), EpsyyElement(:), EpszzElement(:), &
          EpsxyElement(:), EpsxzElement(:), EpsyzElement(:)

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at,st,at0
#else
     REAL(KIND=dp) :: at,st,CPUTime,at0,RealTime
#endif

     REAL(KIND=dp) :: StabParam1, StabParam2, S(3,3), PlyThick, &
          Tx(3), Ty(3), Tz(3), l, Weight3(3), Weight4(4), &
          Fx, Fy, Fz, Mx, My, Mz, &
          Eps(3,3), Kap(3,3), Nten(3,3), Mten(3,3), &
          Amatrix(3,3), Bmatrix(3,3), Dmatrix(3,3), Astarmatrix(2,2), &
          Gdrilling(1,1), Transformation(3,3), T5(5,5), T0(3,3), T0S(2,2), &
          LamRF, Nvector(3), NtenMaterial(2,2), MtenMaterial(2,2), ss

     LOGICAL :: GotIt, GotForceBC

     LOGICAL :: Isotropic = .FALSE.

     TYPE(ValueList_t), POINTER :: SolverParams, Material, BodyForce, BC

     REAL(KIND=dp) :: xn, yn, zn, xp, yp, zp, r

     REAL(KIND=dp) :: Cscalar, Gscalar, Lambda, CalLambda, ArcLength, ArcLengthScale, &
         Radius, MaximumDeflection, Maxcomponent
     INTEGER :: ArcLengthSteps, ArcLengthStep

     SAVE STIFF, MASS, LoadX, LoadY, LoadZ, LoadN, &
          FORCE, ElementNodes, Poisson, Density, Young, Thickness, &
          Tension, AllocationsDone, DAMP, LocalDeflection, &
          SxxNodal, SyyNodal, SzzNodal, SxyNodal, SxzNodal, SyzNodal, &
          EpsxxNodal, EpsyyNodal, EpszzNodal, EpsxyNodal, EpsxzNodal, &
          EpsyzNodal, LoadVector
!==================================================================================

     CALL Info('ShellSolver','-----------------------------------------',Level=6)
     CALL Info('ShellSolver','Reissner-Mindlin facet shell solver',Level=4)

!    Get variables needed for solution:
!    ----------------------------------
     Deflection     => Solver % Variable % Values
     DeflectionPerm => Solver % Variable % Perm

     LocalNodes = Model % NumberOfNodes
     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS

     Norm = Solver % Variable % Norm
     Unorm = 0.0d0

     DO i = 1, Solver % Mesh % NumberOfNodes
        j = DeflectionPerm(i)
        IF(j < 1) CYCLE
       Unorm = Unorm + Deflection(6*j-5)**2 + Deflection(6*j-4)**2 + Deflection(6*j-3)**2
     END DO
     Unorm = SQRT( Unorm )

!==================================================================================

!    Allocate some permanent storage, this is done first time only:
!    --------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN

       ! Only triangles and quads may be used for the shell solver
       N = MIN( Model % MaxElementNodes, 4 )

       ALLOCATE( ElementNodes % x( N ),   &
                 ElementNodes % y( N ),   &
                 ElementNodes % z( N ),   &
                 FORCE( 6*N ),         &
                 STIFF( 6*N, 6*N ), &
                 MASS( 6*N, 6*N ), &
                 DAMP( 6*N, 6*N ), &
                 LoadX( N ), LoadY( N ), &
                 LoadZ( N ), LoadN( N ), &
                 Poisson( N ), Young( N ), &
                 Density ( N ), Thickness( N ), &
                 Tension( N ), &
                 LocalDeflection( 6*N ), &
                 LoadVector( 6, N ), &
                 STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'ShellSolver',  'Memory allocation error, Aborting.' )
       END IF

       NULLIFY( PointLoad )

       NULLIFY( SxxNodal, SyyNodal, SzzNodal, &
          SxyNodal, SxzNodal, SyzNodal, &
          EpsxxNodal, EpsyyNodal, EpszzNodal, &
          EpsxyNodal, EpsxzNodal, EpsyzNodal )

       AllocationsDone = .TRUE.
     END IF

!==================================================================================

!    Do some additional initialization, and go for it:
!    -------------------------------------------------
     at  = CPUTime()
     at0 = RealTime()
     Isotropic = .TRUE.

!----------------------------------------------------------------------------------------
!     Non-linear iteration etc.
!----------------------------------------------------------------------------------------
     SolverParams => GetSolverParams()

      NonLinConvTol = GetConstReal( SolverParams, &
          'Nonlinear System Convergence Tolerance', GotIt )
      IF( .NOT. GotIt ) NonLinConvTol = 1.0e-6

      NonLinMaxIt = GetInteger( SolverParams, &
          'Nonlinear System Max Iterations', GotIt )
      IF( .NOT. GotIt ) NonLinMaxIt = 1

      LargeDeflection =  GetLogical( SolverParams, 'Large Deflection', GotIt )
      IF( .NOT. GotIt ) LargeDeflection = .FALSE.

      IF( .NOT. LargeDeflection ) THEN
        NonLinMaxIt = 1
        StressComputation = GetLogical( SolverParams, 'Stress Computation', GotIt )
        StressComputation = StressComputation .OR. &
             GetLogical( SolverParams, 'Calculate Stresses', GotIt )
      ELSE
        CALL Info( 'ShellSolver',  'Using Von Karman strains for large deflections!' )
        StressComputation = .FALSE.
      END IF

      StabilityAnalysis = GetLogical( SolverParams, 'Stability Analysis', GotIt )
      IF( .NOT. GotIt ) StabilityAnalysis = .FALSE.

      IF( StabilityAnalysis ) THEN
        LargeDeflection = .FALSE.
        NonLinMaxIt = 2
        NOFEigValsBackup = Solver % NOFEigenValues
      END IF

      StabParam1 = GetConstReal( SolverParams, &
           'Shear Stabilization Parameter', GotIt )
      IF( .NOT.GotIt ) THEN
         CALL Info( 'ShellSolver', 'Shear stabilization parameter undefined.',Level=6 )
         CALL Info( 'ShellSolver', 'Using default value 1.0',Level=6 )
         StabParam1 = 1.0d0
      END IF

      StabParam2 = GetConstReal( SolverParams, &
           'Drilling Stabilization Parameter', GotIt )
      IF( .NOT.GotIt ) THEN
         CALL Info( 'ShellSolver', 'Drilling stabilization parameter undefined',Level=6)
         CALL Info( 'ShellSolver', 'Using default value 1.0',Level=6 )
         StabParam2 = 1.0d0
      END IF
       
!      Get element local matrix, and rhs vector:
!      -----------------------------------------
       Nvector = 0.0d0
!----------------------------------------------------------------------------------------
!     Non-linear iteration starts here
!----------------------------------------------------------------------------------------

      CALL SolveNonLinear()

      CALL Info('ShellSolver','All done',Level=6)
      CALL Info('ShellSolver','-----------------------------------------',Level=6)


!------------------------------------------------------------------------------
 
   CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE SolveNonLinear
!------------------------------------------------------------------------------
     DO NonLinIter = 1,NonLinMaxIt
       CALL Info('ShellSolver','--------------------------------------------------------')
       WRITE( Message,'(A,I4)') 'Newton iteration',nonliniter
       CALL Info('ShellSolver',Message)
       CALL Info('ShellSolver','--------------------------------------------------------')
	   
       IF ( StabilityAnalysis ) THEN
         SELECT CASE( NonLinIter )
         CASE( 1 )
           Solver % NOFEigenValues = 0
         CASE DEFAULT
           Solver % NOFEigenValues = NOFEigValsBackup
         END SELECT
       END IF

!    Do the assembly:
!    ----------------
     at = CPUTime()
     CALL DefaultInitialize()
     CALL BulkAssembly()
     CALL DefaultFinishBulkAssembly()

     CALL BCAssembly()
     CALL DefaultFinishAssembly()
     CALL ConcentratedLoads()

     ForceVector = LoadScale*ForceVector

!--------------------------------------------------------------------------

!    Dirichlet boundary conditions:
!    ------------------------------
     CALL SetDirichletBCs()

     at = CPUTime() - at
     WRITE(Message,'(a,F8.2)') ' Assembly: (s)', at
     CALL Info('ShellSolve',Message)

!    Solve the system and we are done:
!    ---------------------------------
     st = CPUTime()

     PrevNorm = Norm
     PrevUNorm = Unorm

!    First iterate is a special case, solving for Lamda=1:
!    -----------------------------------------------------
     Norm = DefaultSolve()


     Unorm = 0.0d0
     DO i = 1, Solver % Mesh % NumberOfNodes
        j = DeflectionPerm(i)
        IF(j < 1) CYCLE
        Unorm = MAX(Unorm, SQRT( Deflection(6*j-5)**2 + &
             Deflection(6*j-4)**2 + Deflection(6*j-3)**2) )
        Maxcomponent = MAX( Unorm, MAX(MAX(Deflection(6*j-5), &
             Deflection(6*j-4)),Deflection(6*j-3)) )
     END DO
     !PRINT *,'Max deflection =',Unorm
     MaximumDeflection = Unorm

     Unorm = 0.0d0
     DO i = 1, Solver % Mesh % NumberOfNodes
        j = DeflectionPerm(i)
        IF(j < 1) CYCLE
        Unorm = Unorm + Deflection(6*j-5)**2 + &
              Deflection(6*j-4)**2 + Deflection(6*j-3)**2
     END DO
     Unorm = SQRT( Unorm )

     IF( ABS(Norm + PrevNorm) > 1.0e-8 )  RelChange = &
         ABS(Norm - PrevNorm)/ABS(Norm + PrevNorm)

     IF( ABS(UNorm + PrevUNorm) > 1.0e-8) RelUChange = &
       ABS(UNorm - PrevUnorm)/ABS(UNorm + PrevUNorm)

     st = CPUTime() - st

     WRITE(Message,'(a,F8.2)') 'Solve: (s)', st
     CALL Info('ShellSolve',Message)

     WRITE(Message,'(a,2F8.3)') 'Relative Change = ',RelChange, RelUChange 
     CALL Info('ShellSolve',Message)

     IF( RelChange < NonlinConvTol ) EXIT
  END DO    ! NonLinIter

  StressComputation = GetLogical( SolverParams, 'Stress Computation', GotIt )
  StressComputation = StressComputation .OR. &
            GetLogical( SolverParams, 'Calculate Stresses', GotIt )

  IF( StressComputation .AND. .NOT. EigenOrHarmonicAnalysis() ) THEN
    CALL Info('ShellSolve','Entering stress calculation routines...')
    CALL CalculateStresses()
  END IF ! If Stress Calculation
!-----------------------------------------------------------------------------
  END SUBROUTINE SolveNonLinear
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
   SUBROUTINE BulkAssembly
!-----------------------------------------------------------------------------
     CALL StartAdvanceOutput('ShellSolve', 'Assembly:')
     DO t=1,Solver % NumberOfActiveElements
!-----------------------------------------------------------------------------

       CALL AdvanceOutput(t,Solver % NumberOFActiveElements)

!-----------------------------------------------------------------------------
       CurrentElement => GetActiveElement( t )

       IF( CurrentElement % TYPE % ElementCode /= 303 .AND. &
           CurrentElement % TYPE % ElementCode /= 404 ) THEN
         CALL Warn('ShellSolver','Implemented only for elements 303 and 404!')
         CYCLE
       END IF

       n = GetElementNOFNodes()
       NodeIndexes => CurrentElement % NodeIndexes

       LocalDeflection = 0.0d0
       DO i = 1,n
         k = DeflectionPerm(NodeIndexes(i))
         DO j = 1,6
           LocalDeflection(6*(i-1)+j) = Deflection(6*(k-1)+j)
         END DO
       END DO

!      Check element type:
!      -------------------
       IF( .NOT.( ( n == 3 ) .OR. ( n == 4 ) ) ) THEN
          CALL Fatal( 'ShellSolver', 'Illegal number of nodes. Aborting.' )
       END IF
  
       CALL GetElementNodes( ElementNodes )

       LoadN(1:n) = GetReal( SolverParams, 'Load Scale Factor', GotIt )
       LoadScale = LoadN(1)
       IF( .NOT. GotIt ) LoadScale = 1.0d0

!      Nodal loads:
!      ------------
       BodyForce => GetBodyForce(CurrentElement, GotIt)

       LoadVector = 0.0d0

       IF( GotIt ) THEN
         LoadX(1:n) = GetReal( BodyForce, 'Body Force 1', GotIt )
         LoadY(1:n) = GetReal( BodyForce, 'Body Force 2', GotIt )
         LoadZ(1:n) = GetReal( BodyForce, 'Body Force 3', GotIt )
         LoadN(1:n) = GetReal( BodyForce, 'Pressure', GotIt )
         LoadN(1:n) = LoadN(1:n) + GetReal( BodyForce, 'Normal Pressure', GotIt )
       ELSE         
         LoadX(1:n) = 0.0_dp
         LoadY(1:n) = 0.0_dp
         LoadZ(1:n) = 0.0_dp
         LoadN(1:n) = 0.0_dp
         LoadN(1:n) = 0.0_dp
       END IF
         

!      Material data:
!      --------------
       Material => GetMaterial()

       Density(1:n) = GetReal( Material, 'Density', GotIt )
       IF( .NOT.GotIt ) THEN
          Density = 0.0d0
          IF( TransientSimulation .OR. (Solver % NOfEigenvalues > 0) ) &
               CALL Fatal( 'ShellSolver', 'Density required' )
       END IF

       Poisson(1:n) = GetReal( Material, 'Poisson ratio', GotIt )
       IF( Isotropic .AND. (.NOT.GotIt) ) &
                     CALL Fatal( 'ShellSolver', 'Poisson ratio undefined' )

       Young(1:n) = GetReal( Material, 'Youngs modulus', GotIt )
       IF( Isotropic .AND. (.NOT.GotIt) ) &
                     CALL Fatal( 'ShellSolver', 'Youngs modulus undefined' )

       Thickness(1:n) = GetReal( Material, 'Thickness', GotIt )
       IF( Isotropic .AND. (.NOT.GotIt) ) &
                     CALL Fatal( 'ShellSolver', 'Thickness undefined' )

       Tension(1:n) = GetReal( Material, 'Tension', GotIt )
       IF( .NOT. GotIt ) Tension = 0.0d0

       CALL LocalMatrix(  STIFF, DAMP, MASS, &
            FORCE, LoadX, LoadY, LoadZ, LoadN, CurrentElement, n, &
            ElementNodes, StabParam1, StabParam2, t, Poisson,     &
            Young, LocalDeflection, LargeDeflection,  &
            StabilityAnalysis, Nvector )

       IF( TransientSimulation ) THEN
          CALL Default2ndOrderTime( MASS, DAMP, STIFF, FORCE )
       END IF

!      Update global matrix and rhs vector from local matrix & vector:
!      ---------------------------------------------------------------
       CALL DefaultUpdateEquations( STIFF, FORCE )
       IF ( EigenOrHarmonicAnalysis() ) CALL DefaultUpdateMass( MASS )
     END DO
!-----------------------------------------------------------------------------
    END SUBROUTINE BulkAssembly
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
    SUBROUTINE BCAssembly()
!-----------------------------------------------------------------------------
      NumberOfElementNodes = n
!     Neumann & Newton boundary conditions:
!     -------------------------------------
      DO t = 1, Solver % Mesh % NumberOfBoundaryElements
        
        CurrentElement => GetBoundaryElement( t )
        IF( .NOT. ActiveBoundaryElement() ) CYCLE
        IF ( CurrentElement % TYPE % ElementCode == 101 ) CYCLE
        BC => GetBC()
        IF ( .NOT. ASSOCIATED( BC ) ) CYCLE

        n = GetElementNOFNodes()
        NodeIndexes => CurrentElement % NodeIndexes

        LoadVector = 0.0d0

        LoadVector( 1, 1:n ) =  GetReal( BC, 'Force 1', GotIt )
        GotForceBC = GotIt

        LoadVector( 2, 1:n ) =  GetReal( BC, 'Force 2', GotIt )
        GotForceBC = GotForceBC .OR. GotIt

        LoadVector( 3, 1:n ) =  GetReal( BC, 'Force 3', GotIt )
        GotForceBC = GotForceBC .OR. GotIt

        LoadVector( 4, 1:n ) =  GetReal( BC, 'Force 4', GotIt )
        GotForceBC = GotForceBC .OR. GotIt

        LoadVector( 5, 1:n ) =  GetReal( BC, 'Force 5', GotIt )
        GotForceBC = GotForceBC .OR. GotIt

        LoadVector( 6, 1:n ) =  GetReal( BC, 'Force 6', GotIt )
        GotForceBC = GotForceBC .OR. GotIt

        IF( .NOT.GotForceBC ) CYCLE

        CALL StressBoundary( STIFF, FORCE, LoadVector, &
             CurrentElement, n, ElementNodes )
                    
        CALL DefaultUpdateEquations( STIFF, FORCE )
     END DO
!-----------------------------------------------------------------------------
    END SUBROUTINE BCAssembly
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
    SUBROUTINE ConcentratedLoads()
!-----------------------------------------------------------------------------
     bf_id = ListGetInteger( Model % Bodies(1) % Values, 'Body Force',GotIt )
     IF( GotIt ) THEN
       PointLoad => ListGetConstRealArray( Model % BodyForces( bf_id ) % Values, &
           'Point Load', GotIt)
     END IF

     IF( .NOT. GotIt ) RETURN

     nPL = SIZE( PointLoad ) / 9
       
     WRITE( Message,'(A,I9)' ) 'Number of point loads set',nPL
     CALL Info('ShellSolve',Message)

     DO iPL = 1, nPL
       xp = PointLoad(iPL, 1)
       yp = PointLoad(iPL, 2)
       zp = PointLoad(iPL, 3)
       Fx = PointLoad(iPL, 4)
       Fy = PointLoad(iPL, 5)
       Fz = PointLoad(iPL, 6)
       Mx = PointLoad(iPL, 7)
       My = PointLoad(iPL, 8)
       Mz = PointLoad(iPL, 9)

       DO i = 1, Solver % Mesh % NumberOfNodes
         xn = Solver % Mesh % Nodes % x(i)
         yn = Solver % Mesh % Nodes % y(i)
         zn = Solver % Mesh % Nodes % z(i)
         r = SQRT( (xn-xp)**2 + (yn-yp)**2 + (zn-zp)**2 )

         IF ( r < 1.0d-8 ) THEN
           k = DeflectionPerm( i )
           IF(k < 1) CYCLE
           ForceVector( 6*k-5 ) = ForceVector( 6*k-5 ) + Fx
           ForceVector( 6*k-4 ) = ForceVector( 6*k-4 ) + Fy
           ForceVector( 6*k-3 ) = ForceVector( 6*k-3 ) + Fz
           ForceVector( 6*k-2 ) = ForceVector( 6*k-2 ) + Mx
           ForceVector( 6*k-1 ) = ForceVector( 6*k-1 ) + My
           ForceVector( 6*k-0 ) = ForceVector( 6*k-0 ) + Mz
         END IF
       END DO
     END DO
!-----------------------------------------------------------------------------
   END SUBROUTINE ConcentratedLoads
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
   SUBROUTINE SetDirichletBCs()
!-----------------------------------------------------------------------------
     CALL DefaultDirichletBCs()
!-----------------------------------------------------------------------------
   END SUBROUTINE SetDirichletBCs
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!
!                        ===================================
!                        S T R E S S   C O M P U T A T I O N
!                        ===================================
!
!-----------------------------------------------------------------------------
   SUBROUTINE CalculateStresses()
!-----------------------------------------------------------------------------
      INTEGER :: n, isz

      at  = CPUTime()
      at0 = RealTime()

!      Allocate memory for the local stresses:
!      ---------------------------------------
       n = Solver % Mesh % NumberOfBulkElements

       ALLOCATE( Weights( n, NumberOfElementNodes ) )
       Weights = 0.0d0

    
       ALLOCATE( SxxElement( n ), SyyElement( n ), SzzElement( n ), &
                 SxyElement( n ), SxzElement( n ), SyzElement( n ) )

       SxxElement = 0.0d0
       SyyElement = 0.0d0
       SzzElement = 0.0d0
       SxyElement = 0.0d0
       SxzElement = 0.0d0
       SyzElement = 0.0d0

       ALLOCATE( EpsxxElement( n ), EpsyyElement( n ), EpszzElement( n ), &
                 EpsxyElement( n ), EpsxzElement( n ), EpsyzElement( n ) )

       EpsxxElement = 0.0d0
       EpsyyElement = 0.0d0
       EpszzElement = 0.0d0
       EpsxyElement = 0.0d0
       EpsxzElement = 0.0d0
       EpsyzElement = 0.0d0


!      Then, compute the element stresses:
!      ----------------------------------- 

       DO t = 1, Solver % NumberOfActiveElements

          IF( RealTime() - at0 > 1.0 ) THEN
             WRITE(*,'(a,i3,a)') ' Stress calculation:', &
                  INT(100.0 - 100.0*(Solver % &
                  NumberOfActiveElements - t) &
                  /(1.0 * Solver % NumberOfActiveElements) ),' % done'
             at0 = RealTime()
          END IF

          CurrentElement => GetActiveElement(t)
          n = GetElementNOFNodes( CurrentElement )
          NodeIndexes => CurrentElement % NodeIndexes
  
          CALL GetElementNodes( ElementNodes )

          LocalDeflection = 0.0d0
          DO i = 1,n
             k = DeflectionPerm(NodeIndexes(i))
             DO j = 1,6
                LocalDeflection(6*(i-1)+j) = Deflection(6*(k-1)+j)
             END DO
          END DO

!         Compute the local stresses (constant for each element):  
!         -------------------------------------------------------
          CALL LocalStress( CurrentElement, n, ElementNodes, &
            StabParam1, StabParam2, LocalDeflection, Weight3, Weight4, &
              Eps, Kap, Nten, NtenMaterial, Mten, MtenMaterial, Young, &
                Poisson, Thickness, LargeDeflection )

          IF( ListGetLogical( SolverParams, 'Compute Strains', GotIt) ) THEN
            SxxElement(t) = Nten(1,1)
            SyyElement(t) = Nten(2,2)
            SzzElement(t) = Nten(3,3)
            SxyElement(t) = Nten(1,2)
            SxzElement(t) = Nten(1,3)
            SyzElement(t) = Nten(2,3)
             
            EpsxxElement(t) = Eps(1,1)
            EpsyyElement(t) = Eps(2,2)
            EpszzElement(t) = Eps(3,3)
            EpsxyElement(t) = Eps(1,2)
            EpsxzElement(t) = Eps(1,3)
            EpsyzElement(t) = Eps(2,3)
          END IF

          IF( ListGetLogical( SolverParams, 'Compute Curvatures', GotIt) ) THEN
            SxxElement(t) = Mten(1,1)
            SyyElement(t) = Mten(2,2)
            SzzElement(t) = Mten(3,3)
            SxyElement(t) = Mten(1,2)
            SxzElement(t) = Mten(1,3)
            SyzElement(t) = Mten(2,3)
             
            EpsxxElement(t) = Kap(1,1)
            EpsyyElement(t) = Kap(2,2)
            EpszzElement(t) = Kap(3,3)
            EpsxyElement(t) = Kap(1,2)
            EpsxzElement(t) = Kap(1,3)
            EpsyzElement(t) = Kap(2,3)
          END IF

          SELECT CASE( NumberOfElementNodes )
          CASE( 3 )
             Weights(t,1:3) = Weight3(1:3)
          CASE( 4 )
             Weights(t,1:4) = Weight4(1:4)
          END SELECT

        END DO


        isz = MAXVAL( DeflectionPerm )
        
        IF( .NOT.ASSOCIATED( VariableGet( Solver % Mesh % Variables, &
            'Stress.xx') ) ) THEN

          ALLOCATE( SxxNodal( isz ), SyyNodal( isz ), SzzNodal( isz ), &
                    SxyNodal( isz ), SxzNodal( isz ), SyzNodal( isz ) )

          ALLOCATE( EpsxxNodal( isz ), EpsyyNodal( isz ), EpszzNodal( isz ), &
	            EpsxyNodal( isz ), EpsxzNodal( isz ), EpsyzNodal( isz ) )


          PSolver => Solver

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Stress.xx', 1, SxxNodal,DeflectionPerm )

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Stress.yy', 1, SyyNodal,DeflectionPerm )

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Stress.zz', 1, SzzNodal,DeflectionPerm )

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Stress.xy', 1, SxyNodal,DeflectionPerm )

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Stress.xz', 1, SxzNodal,DeflectionPerm )

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Stress.yz', 1, SyzNodal,DeflectionPerm )

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Epsilon.tot.xx', 1, EpsxxNodal,DeflectionPerm )

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Epsilon.tot.yy', 1, EpsyyNodal,DeflectionPerm )

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Epsilon.tot.zz', 1, EpszzNodal,DeflectionPerm )

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Epsilon.tot.xy', 1, EpsxyNodal,DeflectionPerm )

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Epsilon.tot.xz', 1, EpsxzNodal,DeflectionPerm )

          CALL VariableAdd( Solver % Mesh % Variables, &
               Solver % Mesh, PSolver, 'Epsilon.tot.yz', 1, EpsyzNodal,DeflectionPerm )
        END IF


!       Average nodal values:
!       ---------------------
        SxxNodal = 0.0d0
        SyyNodal = 0.0d0
        SzzNodal = 0.0d0
        SxyNodal = 0.0d0
        SxzNodal = 0.0d0
        SyzNodal = 0.0d0
        EpsxxNodal = 0.0d0
        EpsyyNodal = 0.0d0
        EpszzNodal = 0.0d0
        EpsxyNodal = 0.0d0
        EpsxzNodal = 0.0d0
        EpsyzNodal = 0.0d0

        ALLOCATE( Referenced( isz ) )
        Referenced = 0.0d0

        DO i = 1, Solver % NumberOfActiveElements
          CurrentElement => GetActiveElement(i)
          n = GetElementNOFNodes( CurrentElement )
          NodeIndexes => CurrentElement % NodeIndexes
          CALL GetElementNodes( ElementNodes )
          DO j = 1, n
            k = DeflectionPerm( NodeIndexes(j) )
            IF( k == 0 ) CYCLE

            Referenced(k) = Referenced(k) + 1
            SxxNodal(k) = SxxNodal(k) + SxxElement(i)
            SyyNodal(k) = SyyNodal(k) + SyyElement(i)
            SzzNodal(k) = SzzNodal(k) + SzzElement(i)
            SxyNodal(k) = SxyNodal(k) + SxyElement(i)
            SxzNodal(k) = SxzNodal(k) + SxzElement(i)
            SyzNodal(k) = SyzNodal(k) + SyzElement(i)
            EpsxxNodal(k) = EpsxxNodal(k) + EpsxxElement(i)
            EpsyyNodal(k) = EpsyyNodal(k) + EpsyyElement(i)
            EpszzNodal(k) = EpszzNodal(k) + EpszzElement(i)
            EpsxyNodal(k) = EpsxyNodal(k) + EpsxyElement(i)
            EpsxzNodal(k) = EpsxzNodal(k) + EpsxzElement(i)
            EpsyzNodal(k) = EpsyzNodal(k) + EpsyzElement(i)
          END DO
        END DO

        DO i=1,SIZE(Referenced)
          IF ( Referenced(i) > 0 ) THEN
             SxxNodal(i) = SxxNodal(i) / Referenced(i)
             SyyNodal(i) = SyyNodal(i) / Referenced(i)
             SzzNodal(i) = SzzNodal(i) / Referenced(i)
             SxyNodal(i) = SxyNodal(i) / Referenced(i)
             SxzNodal(i) = SxzNodal(i) / Referenced(i)
             SyzNodal(i) = SyzNodal(i) / Referenced(i)
             EpsxxNodal(i) = EpsxxNodal(i) / Referenced(i)
             EpsyyNodal(i) = EpsyyNodal(i) / Referenced(i)
             EpszzNodal(i) = EpszzNodal(i) / Referenced(i)
             EpsxyNodal(i) = EpsxyNodal(i) / Referenced(i)
             EpsxzNodal(i) = EpsxzNodal(i) / Referenced(i)
             EpsyzNodal(i) = EpsyzNodal(i) / Referenced(i)
          END IF
        END DO

        at = CPUTime() - at

        WRITE(Message,'(a,F8.2)') ' Stress calculation: (s)', at
        CALL Info('ShellSolver',Message)

!       Finally, release the auxiliary arrays:
!       -------------------------------------- 
        DEALLOCATE( SxxElement, SyyElement, SzzElement, SxyElement, &
             SxzElement, SyzElement, Weights, Referenced, &
             EpsxxElement, EpsyyElement, EpszzElement, EpsxyElement, &
             EpsxzElement, EpsyzElement )
!------------------------------------------------------------------------------
   END SUBROUTINE CalculateStresses
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE StressBoundary( STIFF, FORCE, LOAD, Element, n, Nodes )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:,:)
       TYPE(Element_t), POINTER :: Element
       INTEGER :: n
       TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
            detJ, s, u, v, w, LoadAtIp(6)
       INTEGER :: t, i, j
       LOGICAL :: stat
       TYPE( GaussIntegrationPoints_t ), TARGET :: IntegStuff

       STIFF = 0.0d0
       FORCE = 0.0d0
       IntegStuff = GaussPoints( element )

       DO t = 1, IntegStuff % n
          U = IntegStuff % u(t)
          V = IntegStuff % v(t)
          W = IntegStuff % w(t)
          S = IntegStuff % s(t)

          stat = ElementInfo( Element, Nodes, u , v, w, &
               detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

          s = detJ * s

          DO i = 1,6
             LoadAtIp(i) = SUM( LOAD(i,1:n) * Basis(1:n) )
          END DO

          DO j = 1,N
             DO i = 1,6
                FORCE((j-1)*6+i) = FORCE((j-1)*6+i) + &
                     Basis(j) * LoadAtIp(i) * s
             END DO
          END DO

       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE StressBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
     SUBROUTINE LocalMatrix( STIFF, DAMP, MASS, &
          FORCE, NodalLoadX, NodalLoadY, NodalLoadZ, NodalLoadN, &
          Element, N, Nodes, StabParam1, StabParam2, &
          ElementNumber, NodalPoisson, NodalYoung, &
          LocalDeflection, LargeDeflection, StabilityAnalysis, Nvector )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: STIFF(:,:), DAMP(:,:), MASS(:,:), &
            Amatrix(3,3), Bmatrix(3,3), Dmatrix(3,3), Astarmatrix(2,2)
       REAL(KIND=dp) :: FORCE(:)
       REAL(KIND=dp) :: NodalLoadX(:), NodalLoadY(:), NodalLoadZ(:), &
           NodalLoadN(:), LocalDeflection(:), Nvector(:)
       REAL(KIND=dp) :: StabParam1, StabParam2
       REAL(KIND=dp) :: NodalPoisson(:), NodalYoung(:)
       LOGICAL :: LargeDeflection, StabilityAnalysis
       INTEGER :: N, ElementNumber
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       INTEGER, PARAMETER :: MaxNodes = 4
       INTEGER, PARAMETER :: MaxDofs = 6*MaxNodes

       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
            Kappa(3,MaxDofs), Gammaa(2,MaxDofs), EPS(3,MaxDofs), &
            Omega(1,MaxDofs), Gdrilling(1,1), GradDeflection(2,MaxDofs), & !*
            ZetaStrain(3,maxDofs)

       REAL(KIND=dp) :: detJ, U, V, W, S, SCF, rho, h, &
            LoadN, LoadX, LoadY, LoadZ, Nmatrix(2,2), &
            GradTest(2), GradBasis(2), LV1(3,1), LV2(3,1), TempVec(3,1), &
            NormalForce(2,2), NonLinForce(MaxDofs), LV3(3,1)

       REAL(KIND=dp) :: Transformation(3,3), CopyOfNodes(3,MaxNodes), &
            T0(3,3), T0S(2,2), T5(5,5)

       REAL( KIND=dp ) :: Tblock(6*n,6*n)

       LOGICAL :: Stat
       INTEGER :: i,p,q,t, pk, qk
       TYPE( GaussIntegrationPoints_t ) :: IntegStuff

       REAL(KIND=dp) :: LamDCMatrix(8,8), dWdx, dWdy

       REAL(KIND=dp) :: Moment(2,2), ZetaStress(2,2), EpsStress(2,2), &
            dUdx(3,2), dRdx(3,2)

!------------------------------------------------------------------------------
       FORCE = 0.0d0
       STIFF = 0.0d0
       DAMP  = 0.0d0
       MASS  = 0.0d0

       Kappa          = 0.0d0
       EPS            = 0.0d0
       Gammaa         = 0.0d0
       Omega          = 0.0d0 
       GradDeflection = 0.0d0
       ZetaStrain     = 0.0d0

       LV1 = 0.0d0
       LV2 = 0.0d0
       LV3 = 0.0d0
       TempVec     = 0.0d0
       NormalForce = 0.0d0
       Moment      = 0.0d0
       ZetaStress  = 0.0d0
       EpsStress   = 0.0d0
       NonlinForce = 0.0d0

       IF ( StabilityAnalysis ) THEN
         CALL LocalStress( Element, n, Nodes, StabParam1,    &
           StabParam2, LocalDeflection, Weight3, Weight4,    &
           Eps, Kap, Nten, NtenMaterial, Mten, MtenMaterial, &
           NodalYoung, NodalPoisson, Thickness, LargeDeflection )
       END IF

!      The transformation Xglob -> Xloc is the transpose of the local basis:
!      ---------------------------------------------------------------------
       Transformation = TRANSPOSE( LocalBasis( Nodes, n ) )
       
!      Take a copy of the global node points and switch to the local system:
!      ---------------------------------------------------------------------
       CALL SwitchToLocal( Nodes, CopyOfNodes, Transformation, n )

!      Transform the old deflection into local co-ordinates:
!      -----------------------------------------------------
       Tblock = 0.0d0
       DO i = 1, n
          DO p = 1, 3
             DO q = 1, 3
                Tblock( 6*i-6+p, 6*i-6+q ) = Transformation( p, q )
                Tblock( 6*i-3+p, 6*i-3+q ) = Transformation( p, q )
             END DO
          END DO
       END DO
       LocalDeflection = MATMUL( Tblock, LocalDeflection(1:6*n) )

!      Exchange rotations r1 and -r2:
!      ------------------------------
       Tblock = 0.0d0
       DO i = 1,n
          Tblock( 6*i-5, 6*i-5 ) =  1.0d0
          Tblock( 6*i-4, 6*i-4 ) =  1.0d0
          Tblock( 6*i-3, 6*i-3 ) =  1.0d0
          Tblock( 6*i-2, 6*i-1 ) =  1.0d0
          Tblock( 6*i-1, 6*i-2 ) = -1.0d0
          Tblock( 6*i-0, 6*i-0 ) =  1.0d0
       END DO
       LocalDeflection = MATMUL( Tblock, LocalDeflection(1:6*n) )

!      Select the appropriate quadrature:
!      ----------------------------------
       SELECT CASE( Element % TYPE % NumberOfNodes )
          CASE( 3 )
             IntegStuff = GaussPoints( Element, 1 )
          CASE( 4 )
             IntegStuff = GaussPoints( Element, 4 )
       END SELECT

!      Numerical integration:
!      ----------------------
       DO t = 1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)

!        Basis function values & derivatives at the integration point:
!        -------------------------------------------------------------
         stat = ElementInfo( Element, Nodes, U, V, W, &
                    detJ, Basis, dBasisdx )
         S = S * detJ

!        Material etc. parameters in the integration point:
!        --------------------------------------------------
         LoadX = SUM( NodalLoadX(1:n) * Basis(1:n) )
         LoadY = SUM( NodalLoadY(1:n) * Basis(1:n) )
         LoadZ = SUM( NodalLoadZ(1:n) * Basis(1:n) )
         LoadN = SUM( NodalLoadN(1:n) * Basis(1:n) )
         rho   = SUM( Density(1:n)    * Basis(1:n) )
         h     = SUM( Thickness(1:n)  * Basis(1:n) )

         CALL IsotropicElasticity( Dmatrix, &
             Astarmatrix, NodalPoisson, NodalYoung, Thickness, Basis, n )

         Bmatrix = 0.0d0

         CALL IsotropicInPlaneElasticity( Amatrix, &
               NodalPoisson, NodalYoung, Thickness, Basis, n )

         Gdrilling(1,1) = StabParam2*(Astarmatrix(1,1)+Astarmatrix(2,2))

!        ----------------------------------------------
!        The nodal degrees-of-freedom are organized as
!               (u_x, u_y, u_z, r_x, r_y, r_z)
!        where u is the displacement and r the rotation
!        ----------------------------------------------

!        Gradient of the current deflection:
!        -----------------------------------
         IF( LargeDeflection ) THEN
            dUdx = 0.0d0
            dRdx = 0.0d0

            DO p = 1,n
               DO i = 1,3
                  DO j = 1,2
                     dUdx(i,j) = dUdx(i,j) &
                          + LocalDeflection(6*(p-1)+i) * dBasisdx(p,j)
                     
                     dRdx(i,j) = dRdx(i,j) &
                          + LocalDeflection(6*(p-1)+i+3) * dBasisdx(p,j)
                     
                  END DO
               END DO
            END DO
         END IF

!        Normal force tensor for the current solution:
!        ----------------------------------------------
         IF( LargeDeflection ) THEN
            LV1 = 0.0d0 ! Strain      = e(U) + 0.5*( dU'dU + dW'dW )
            LV2 = 0.0d0 ! Curvature   = e(B) + 0.5*( dU'dB + dB'dU )
            LV3 = 0.0d0 ! Quad.strain = 0.5*( dB'dB )

            ! Linear part of strain and curvature:
            !-------------------------------------
            DO p=1,n
               LV1(1,1) = LV1(1,1) + LocalDeflection(6*(p-1)+1) * dBasisdx(p,1)
               LV1(2,1) = LV1(2,1) + LocalDeflection(6*(p-1)+2) * dBasisdx(p,2)
               LV1(3,1) = LV1(3,1) + LocalDeflection(6*(p-1)+1) * dBasisdx(p,2) & 
                                   + LocalDeflection(6*(p-1)+2) * dBasisdx(p,1) 

               LV2(1,1) = LV2(1,1) + LocalDeflection(6*(p-1)+4) * dBasisdx(p,1)
               LV2(2,1) = LV2(2,1) + LocalDeflection(6*(p-1)+5) * dBasisdx(p,2)
               LV2(3,1) = LV2(3,1) + LocalDeflection(6*(p-1)+4) * dBasisdx(p,2) & 
                                   + LocalDeflection(6*(p-1)+5) * dBasisdx(p,1)  
            END DO

            ! Nonlinear terms:
            !-----------------
            DO q = 1,2
               LV1(1,1) = LV1(1,1) + 0.5d0 * dUdx(q,1)**2
               LV1(2,1) = LV1(2,1) + 0.5d0 * dUdx(q,2)**2
               LV1(3,1) = LV1(3,1) + dUdx(q,1) * dUdx(q,2)


               LV2(1,1) = LV2(1,1) + dUdx(q,1) * dRdx(q,1)
               LV2(2,1) = LV2(2,1) + dUdx(q,2) * dRdx(q,2)
               LV2(3,1) = LV2(3,1) + dUdx(q,1) * dRdx(q,2) &
                                   + dUdx(q,2) * dRdx(q,1)
               
               LV3(1,1) = LV3(1,1) + 0.5d0 * dRdx(q,1)**2
               LV3(2,1) = LV3(2,1) + 0.5d0 * dRdx(q,2)**2
               LV3(3,1) = LV3(3,1) + dRdx(q,1) * dRdx(q,2)
            END DO

            LV1(1,1) = LV1(1,1) + 0.5d0 * dUdx(3,1)**2
            LV1(2,1) = LV1(2,1) + 0.5d0 * dUdx(3,2)**2
            LV1(3,1) = LV1(3,1) + dUdx(3,1) * dUdx(3,2)

            ! Normal force:
            !--------------
            TempVec = MATMUL( Amatrix, LV1 ) + MATMUL( Bmatrix, LV2 )
            
            NormalForce(1,1) = TempVec(1,1)
            NormalForce(2,2) = TempVec(2,1)
            NormalForce(1,2) = TempVec(3,1)
            NormalForce(2,1) = TempVec(3,1)

            ! Bending moment:
            !----------------
            TempVec = MATMUL( Bmatrix, LV1 ) + MATMUL( Dmatrix, LV2 )
            
            Moment(1,1) = TempVec(1,1)
            Moment(2,2) = TempVec(2,1)
            Moment(1,2) = TempVec(3,1)
            Moment(2,1) = TempVec(3,1)

            ! Quadratic terms:
            !-----------------
            TempVec = MATMUL( Dmatrix, LV1 )
            EpsStress(1,1) = TempVec(1,1)
            EpsStress(2,2) = TempVec(2,1)
            EpsStress(1,2) = TempVec(3,1)
            EpsStress(2,1) = TempVec(3,1)

            TempVec = MATMUL( Dmatrix, LV3 )
            ZetaStress(1,1) = TempVec(1,1)
            ZetaStress(2,2) = TempVec(2,1)
            ZetaStress(1,2) = TempVec(3,1)
            ZetaStress(2,1) = TempVec(3,1)

            NormalForce = NormalForce + ZetaStress

         END IF

!--------------------------------------------------------------------

!        Bending stiffness:
!        ------------------
         Kappa = 0.0d0
         DO p=1,n
            Kappa(1,6*p-2) = dBasisdx(p,1)
            Kappa(2,6*p-1) = dBasisdx(p,2)
            Kappa(3,6*p-2) = dBasisdx(p,2)
            Kappa(3,6*p-1) = dBasisdx(p,1)

            IF( LargeDeflection ) THEN
               DO i = 1,2
                  j = 6*(p-1)+i
                  Kappa(1,j) = Kappa(1,j) + dRdx(i,1) * dBasisdx(p,1)
                  Kappa(2,j) = Kappa(2,j) + dRdx(i,2) * dBasisdx(p,2)
                  Kappa(3,j) = Kappa(3,j) + dRdx(i,1) * dBasisdx(p,2) &
                                          + dRdx(i,2) * dBasisdx(p,1)
                  j = j+3
                  Kappa(1,j) = Kappa(1,j) + dUdx(i,1) * dBasisdx(p,1)
                  Kappa(2,j) = Kappa(2,j) + dUdx(i,2) * dBasisdx(p,2)
                  Kappa(3,j) = Kappa(3,j) + dUdx(i,1) * dBasisdx(p,2) &
                                          + dUdx(i,2) * dBasisdx(p,1)

               END DO
            END IF
         END DO

         CALL AddEnergy(STIFF, Dmatrix, Kappa, 3, 6*n, s)

!        In-plane stiffness:
!        -------------------
         EPS = 0.0d0
         DO p=1,n
            EPS(1,6*p-5) = dBasisdx(p,1)  
            EPS(2,6*p-4) = dBasisdx(p,2)  
            EPS(3,6*p-5) = dBasisdx(p,2)  
            EPS(3,6*p-4) = dBasisdx(p,1)  

            IF( LargeDeflection ) THEN
               DO i = 1,3
                  j = 6*(p-1)+i
                  EPS(1,j) = EPS(1,j) + dUdx(i,1) * dBasisdx(p,1)        
                  EPS(2,j) = EPS(2,j) + dUdx(i,2) * dBasisdx(p,2)        
                  EPS(3,j) = EPS(3,j) + dUdx(i,1) * dBasisdx(p,2) &
                                      + dUdx(i,2) * dBasisdx(p,1) 
               END DO
            END IF
         END DO

         CALL AddEnergy(STIFF, Amatrix, EPS, 3, 6*n, s)

!        Coupling through the B-matrix:
!        ------------------------------
         CALL AddInnerProducts(STIFF, Bmatrix, &
              TRANSPOSE(EPS), Kappa, 3, 6*n, s)

         CALL AddInnerProducts(STIFF, TRANSPOSE(Bmatrix), &
              TRANSPOSE(Kappa), EPS, 3, 6*n, s)

!        Quadratic strains due to rotation:
!        ----------------------------------
         ZetaStrain = 0.0d0
         IF( LargeDeflection ) THEN
            DO p = 1,n
               DO i = 1,2
                  j = 6*(p-1)+i+3
                  ZetaStrain(1,j) = ZetaStrain(1,j) + dRdx(i,1) * dBasisdx(p,1)
                  ZetaStrain(2,j) = ZetaStrain(2,j) + dRdx(i,2) * dBasisdx(p,2)
                  ZetaStrain(3,j) = ZetaStrain(3,j) + dRdx(i,2) * dBasisdx(p,1) &
                                                    + dRdx(i,1) * dBasisdx(p,2)
               END DO
            END DO
         END IF

         CALL AddInnerProducts(STIFF, Dmatrix, &
              TRANSPOSE(EPS), ZetaStrain, 3, 6*n, s)

         CALL AddInnerProducts(STIFF, Dmatrix, &
              TRANSPOSE(ZetaStrain), EPS, 3, 6*n, s)

!        Shear stiffness (transversal):
!        ------------------------------
         CALL CovariantInterpolation(Gammaa, Basis, &
              Nodes % x(1:n), Nodes % y(1:n), U, V, n)

         CALL ShearCorrectionFactor(SCF, h, Nodes % x(1:n), &
              Nodes % y(1:n), n, StabParam1)

         DO p=1,n
            Gammaa(1:2,6*p-3) = dBasisdx(p,1:2)
         END DO

         CALL AddEnergy(STIFF, Astarmatrix, Gammaa, 2, 6*n, SCF*s)

!        Drilling DOFs (in-plane rotations):
!        -----------------------------------
         Omega = 0.0d0
         DO p = 1,n
            Omega(1,6*p-5) = +dBasisdx(p,2) / 2.0d0  !  u_{x,y}
            Omega(1,6*p-4) = -dBasisdx(p,1) / 2.0d0  ! -u_{y,x}
            Omega(1,6*p-0) =  Basis(p)              !  rotation
         END DO

         CALL AddEnergy(STIFF, Gdrilling, Omega, 1, 6*n, s)

!        Newton lin. terms:
!        -----------------
         DO p = 1,n
            DO q = 1,n
               DO i = 1,2
                  DO j = 1,2
                     DO k = 1,2
                        pk = 6*(p-1)+k
                        qk = 6*(q-1)+k

                        Stiff( pk, qk ) = Stiff( pk, qk ) &
                             + NormalForce(i,j) * dBasisdx(p,i) * dBasisdx(q,j) * s

                        Stiff( pk+3, qk+3 ) = Stiff( pk+3, qk+3 ) &
                             + EpsStress(i,j) * dBasisdx(p,i) * dBasisdx(q,j) * s

                        Stiff( pk, qk+3 ) = Stiff( pk, qk+3 ) &
                             + Moment(i,j) * dBasisdx(p,i) * dBasisdx(q,j) * s

                        Stiff( pk+3, qk ) = Stiff( pk+3, qk ) &
                             + Moment(i,j) * dBasisdx(p,i) * dBasisdx(q,j) * s

                     END DO

                     k = 3
                     pk = 6*(p-1)+k
                     qk = 6*(q-1)+k

                     Stiff( pk, qk ) = Stiff( pk, qk ) &
                          + NormalForce(i,j) * dBasisdx(p,i) * dBasisdx(q,j) * s

                  END DO
               END DO
            END DO
         END DO
         
!        Load vector (only translation):
!        -------------------------------
         DO p=1,n

!           Body force is given in the global cartesian coordinates:
!           --------------------------------------------------------
            FORCE(6*p-5) = FORCE(6*p-5) + LoadX * Basis(p) * h * s
            FORCE(6*p-4) = FORCE(6*p-4) + LoadY * Basis(p) * h * s
            FORCE(6*p-3) = FORCE(6*p-3) + LoadZ * Basis(p) * h * s

!           The normal pressure is given in the local cartesian system:
!           -----------------------------------------------------------
            FORCE(6*p-5:6*p-3) = FORCE(6*p-5:6*p-3) &
                 + Transformation(3,1:3) * LoadN * Basis(p) * s
         END DO

!        Newton lin. terms:
!        ------------------
         IF( LargeDeflection ) THEN

            LV1 = 0.0d0
            LV2 = 0.0d0
            LV3 = 0.0d0

            DO q = 1,2
               LV1(1,1) = LV1(1,1) + 0.5d0 * dUdx(q,1)**2
               LV1(2,1) = LV1(2,1) + 0.5d0 * dUdx(q,2)**2
               LV1(3,1) = LV1(3,1) + dUdx(q,1) * dUdx(q,2)

               LV2(1,1) = LV2(1,1) + dUdx(q,1) * dRdx(q,1)
               LV2(2,1) = LV2(2,1) + dUdx(q,2) * dRdx(q,2)
               LV2(3,1) = LV2(3,1) + dUdx(q,1) * dRdx(q,2) &
                                   + dUdx(q,2) * dRdx(q,1)

               LV3(1,1) = LV3(1,1) + 0.5d0 * dRdx(q,1)**2
               LV3(2,1) = LV3(2,1) + 0.5d0 * dRdx(q,2)**2
               LV3(3,1) = LV3(3,1) + dRdx(q,1) * dRdx(q,2)
            END DO

            LV1(1,1) = LV1(1,1) + 0.5d0 * dUdx(3,1)**2
            LV1(2,1) = LV1(2,1) + 0.5d0 * dUdx(3,2)**2
            LV1(3,1) = LV1(3,1) + dUdx(3,1) * dUdx(3,2)

            TempVec = MATMUL( Amatrix, LV1 ) + MATMUL( Bmatrix, LV2 )
            DO p = 1,6*n
               DO q = 1,3
                  NonLinForce(p) = NonlinForce(p) + TempVec(q,1) * EPS(q,p) * s
               END DO
            END DO

            TempVec = MATMUL( Bmatrix, LV1 ) + MATMUL( Dmatrix, LV2 )
            DO p = 1,6*n
               DO q = 1,3
                  NonLinForce(p) = NonlinForce(p) + TempVec(q,1) * Kappa(q,p) * s   
               END DO
            END DO

!           Quadratic terms:
!           ----------------
            TempVec = MATMUL( Dmatrix, LV1 )
            DO p = 1,6*n
               DO q = 1,3
                  NonLinForce(p) = NonlinForce(p) + TempVec(q,1) * ZetaStrain(q,p) * s
               END DO
            END DO

            TempVec = MATMUL( Dmatrix, LV3 )
            DO p = 1,6*n
               DO q = 1,3
                  NonLinForce(p) = NonlinForce(p) + TempVec(q,1) * EPS(q,p) * s
               END DO
            END DO

            DO p = 1,n
               DO i = 1,2
                  DO j = 1,2
                     DO k = 1,2
                        pk = 6*(p-1)+k

                        NonLinForce( pk ) = NonLinForce( pk ) &
                             + NormalForce(i,j) * dUdx(k,j) * dBasisdx(p,i) * s

                        NonLinForce( pk+3 ) = NonLinForce( pk+3 ) &
                             + EpsStress(i,j) * dRdx(k,j) * dBasisdx(p,i) * s

                        NonLinForce( pk ) = NonLinForce( pk ) &
                             + Moment(i,j) * dRdx(k,j) * dBasisdx(p,i) * s
                        
                        NonLinForce( pk+3 ) = NonLinForce( pk+3 ) &
                           + Moment(i,j) * dUdx(k,j) * dBasisdx(p,i) * s

                     END DO

                     k = 3
                     pk = 6*(p-1)+k

                     NonLinForce( pk ) = NonLinForce( pk ) &
                          + NormalForce(i,j) * dUdx(3,j) * dBasisdx(p,i) * s

                  END DO
               END DO
            END DO
         END IF

!        Mass matrix (only translation):
!        -------------------------------
         IF( .NOT.StabilityAnalysis ) THEN
           DO p = 1,n
             DO q = 1,n
               MASS(6*p-5,6*q-5) = MASS(6*p-5,6*q-5) &
                   + rho * h * Basis(p) * Basis(q) * s
               MASS(6*p-4,6*q-4) = MASS(6*p-4,6*q-4) &
                   + rho * h * Basis(p) * Basis(q) * s
               MASS(6*p-3,6*q-3) = MASS(6*p-3,6*q-3) &
                   + rho * h * Basis(p) * Basis(q) * s
             END DO
           END DO
         END IF

         IF( StabilityAnalysis ) THEN
            DO p = 1,n
               GradTest(1:2) = dBasisdx(p,1:2)
               DO q = 1,n
                  GradBasis(1:2) = dBasisdx(q,1:2)
                  GradBasis = MATMUL( NtenMaterial, GradBasis )
                  
                  MASS(6*p-3,6*q-3) = MASS(6*p-3,6*q-3) &
                       + SUM( GradTest(1:2) * GradBasis(1:2) ) * s
                  
               END DO
            END DO
         END IF
      END DO

!      Restore the original node points:
!      ---------------------------------
       Nodes % x(1:n) = CopyOfNodes(1,1:n)
       Nodes % y(1:n) = CopyOfNodes(2,1:n)
       Nodes % z(1:n) = CopyOfNodes(3,1:n)

!      Finally, we perform some transformations:
!      -----------------------------------------
       Tblock = 0.0d0

!      First, we exchange rotations r1 and -r2:
!      ----------------------------------------
       DO i = 1,n
          Tblock( 6*i-5, 6*i-5 ) =  1.0d0
          Tblock( 6*i-4, 6*i-4 ) =  1.0d0
          Tblock( 6*i-3, 6*i-3 ) =  1.0d0
          Tblock( 6*i-2, 6*i-1 ) =  1.0d0
          Tblock( 6*i-1, 6*i-2 ) = -1.0d0
          Tblock( 6*i-0, 6*i-0 ) =  1.0d0
       END DO

       i = 6*n
       STIFF(1:i,1:i) = MATMUL( STIFF(1:i,1:i), Tblock(1:i,1:i) )
       STIFF(1:i,1:i) = MATMUL( TRANSPOSE( Tblock(1:i,1:i) ), STIFF(1:i,1:i) )

       IF( StabilityAnalysis ) THEN
         MASS(1:i,1:i) = MATMUL( MASS(1:i,1:i), Tblock(1:i,1:i) )
         MASS(1:i,1:i) = MATMUL( TRANSPOSE( Tblock(1:i,1:i) ), MASS(1:i,1:i) )
       END IF

       IF( LargeDeflection ) THEN
          NonlinForce(1:i) = MATMUL( TRANSPOSE(Tblock(1:i,1:i)), NonlinForce(1:i) )
       END IF

!      Finally, return the stiffness matrix w.r.t. the original system:
!      ----------------------------------------------------------------
       Tblock = 0.0d0
       DO i = 1, n
          DO p = 1, 3
             DO q = 1, 3
                Tblock( 6*i-6+p, 6*i-6+q ) = Transformation( p, q )
                Tblock( 6*i-3+p, 6*i-3+q ) = Transformation( p, q )
             END DO
          END DO
       END DO

       i = 6*n
       STIFF(1:i,1:i) = MATMUL( STIFF(1:i,1:i), Tblock(1:i,1:i) )
       STIFF(1:i,1:i) = MATMUL( TRANSPOSE( Tblock(1:i,1:i) ), STIFF(1:i,1:i) )

       IF( LargeDeflection ) THEN
         NonlinForce(1:i) = MATMUL( TRANSPOSE(Tblock(1:i,1:i)), NonlinForce(1:i) )
         FORCE = FORCE + NonlinForce
       END IF

       IF( StabilityAnalysis ) THEN
         MASS(1:i,1:i) = MATMUL( MASS(1:i,1:i), Tblock(1:i,1:i) )
         MASS(1:i,1:i) = MATMUL( TRANSPOSE( Tblock(1:i,1:i) ), MASS(1:i,1:i) )
       END IF

       STIFF = ( STIFF + TRANSPOSE(STIFF) ) / 2.0d0
       MASS  = ( MASS  + TRANSPOSE(MASS) )  / 2.0d0

     END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!==============================================================================

!------------------------------------------------------------------------------
     SUBROUTINE LocalStress( Element, n, Nodes, StabParam1, StabParam2,    &
          LocalDeflection, Weight3, Weight4, Eps, Kap, Nten, NtenMaterial, &
          Mten, MtenMaterial, NodalYoung, NodalPoisson, NodalThickness, &
          LargeDeflection )
!------------------------------------------------------------------------------
       IMPLICIT NONE
       REAL(KIND=dp) :: StabParam1, StabParam2, LocalDeflection(:), &
            Weight3(:), Weight4(:), Eps(3,3), Kap(3,3), NTen(3,3), MTen(3,3), &
            NtenMaterial(2,2), MtenMaterial(2,2),  &
            NodalYoung(:) , NodalPoisson(:), NodalThickness(:)
       INTEGER :: n
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
       INTEGER, PARAMETER :: MaxNodes = 4, MaxDofs = 6*MaxNodes
       LOGICAL :: LargeDeflection
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), &
            Curvature(3,MaxDofs), InPlaneStrain(3,MaxDofs), Gammaa(2,MaxDofs), &
            KappaVector(3), EPSILO(3), TShear(2), T5(5,5), Tmat(3,3), &
            T0(3,3), T0S(2,2), Omega(1,MaxDofs), Gdrilling(1,1), &
            GammaVector(2), EpsMaterial(3), KappaMaterial(3), ToLayer(5,5), &
            InvT0(3,3), AMatrix(3,3), BMatrix(3,3), DMatrix(3,3), AStarmatrix(2,2), &
            LAMDCMatrix(8,8), NVec(3), MVec(3)
       REAL( KIND=dp ) :: detJ, U, V, W, Kappa, Q5(5,5), &
            h, Transformation(3,3), CopyOfNodes(3,MaxNodes), &
            delta, ConstantLoadVec(3), VariableLoadVec(3)

       REAL( KIND=dp ) :: Tblock(6*n,6*n)

       LOGICAL :: Stat
       INTEGER :: i, p, q, k ,j
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: z1, z2, theta, ct, st, SCF, xn, yn, zn
       REAL(KIND=dp) :: side1(2), side2(2)
       REAL(KIND=dp) :: dWdx, dWdy
       REAL(KIND=dp) :: LV1(3,1), LV2(3,1), LV3(3,1), dUdx(3,2), dRdx(3,2)
!------------------------------------------------------------------------------

       Curvature       = 0.0d0
       InPlaneStrain   = 0.0d0
       TShear = 0.0d0

!      The transformation Xglob -> Xloc is the transpose of the local basis:
!      ---------------------------------------------------------------------
       Transformation = TRANSPOSE( LocalBasis( Nodes, n ) )

!      Take a copy of the node points and switch to the local system:
!      --------------------------------------------------------------
       CALL SwitchToLocal( Nodes, CopyOfNodes, Transformation, n )
       
!      Let us first perform some transformations:
!      ------------------------------------------

!      The dof-vector w.r.t. the local system:
!      ---------------------------------------
       Tblock = 0.0d0
       DO i = 1, n
          DO p = 1, 3
             DO q = 1, 3
                Tblock( 6*i-6+p, 6*i-6+q ) = Transformation( p, q )
                Tblock( 6*i-3+p, 6*i-3+q ) = Transformation( p, q )
             END DO
          END DO
       END DO

       LocalDeflection = MATMUL( Tblock, LocalDeflection )

!      Exchange rotations r1 and -r2:
!      ------------------------------
       Tblock = 0.0d0
       DO i = 1,n
          Tblock( 6*i-5, 6*i-5 ) =  1.0d0
          Tblock( 6*i-4, 6*i-4 ) =  1.0d0
          Tblock( 6*i-3, 6*i-3 ) =  1.0d0
          Tblock( 6*i-2, 6*i-1 ) =  1.0d0
          Tblock( 6*i-1, 6*i-2 ) = -1.0d0
          Tblock( 6*i-0, 6*i-0 ) =  1.0d0
       END DO

       LocalDeflection = MATMUL( Tblock, LocalDeflection )

!      Select the stress evaluation point:
!      -----------------------------------
       SELECT CASE( Element % TYPE % NumberOfNodes )
          CASE( 3 )
             U = 1.0d0/3.0d0
             V = 1.0d0/3.0d0
             W = 0.0d0
          CASE( 4 )
             U = 0.0d0
             V = 0.0d0
             W = 0.0d0
       END SELECT

!      Basis function values & derivatives in the stress eval-point:
!      -------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, &
                    detJ, Basis, dBasisdx )

!      Material parameters in the stress evaluation point:
!      ---------------------------------------------------
       h = SUM( Thickness(1:n) * Basis(1:n) )

       DO p=1,n
          Curvature(1,6*p-2) = dBasisdx(p,1)
          Curvature(2,6*p-1) = dBasisdx(p,2)
          Curvature(3,6*p-2) = dBasisdx(p,2)
          Curvature(3,6*p-1) = dBasisdx(p,1)
       END DO

       DO p=1,n
          InPlaneStrain(1,6*p-5) = dBasisdx(p,1)
          InPlaneStrain(2,6*p-4) = dBasisdx(p,2)
          InPlaneStrain(3,6*p-5) = dBasisdx(p,2)
          InPlaneStrain(3,6*p-4) = dBasisdx(p,1)
       END DO

!      Shear strains (transversal):
!      ----------------------------
       CALL CovariantInterpolation( Gammaa, Basis, &
            Nodes % x(1:n), Nodes % y(1:n), U, V, n )

       CALL ShearCorrectionFactor( SCF, h, Nodes % x(1:n), &
            Nodes % y(1:n), n, StabParam1 )

       DO p = 1,n
          Gammaa(1:2,6*p-3) = dBasisdx(p,1:2)
       END DO

!      Drilling DOFs (in-plane rotations):
!      -----------------------------------
       DO p = 1,n
          Omega(1,6*p-5) =  dBasisdx(p,2)/2.0d0  !  u_{x,y}
          Omega(1,6*p-4) = -dBasisdx(p,1)/2.0d0  ! -u_{y,x}
          Omega(1,6*p-0) =  Basis(p)             ! rotation
       END DO

!      Venym�- ja k�yristym�vektorit lokaalissa koordinaatistossa:
!      -----------------------------------------------------------
       KappaVector = MATMUL( Curvature(1:3,1:6*n), LocalDeflection(1:6*n) )
       EPSILO = MATMUL( InPlaneStrain(1:3,1:6*n), LocalDeflection(1:6*n) )
       GammaVector = MATMUL( Gammaa(1:2,1:6*n), LocalDeflection(1:6*n) )

!      VonKarman strains:         
!========================

       IF( LargeDeflection ) THEN

         dUdx = 0.0d0
         dRdx = 0.0d0

         LV1 = 0.0d0 ! Strain      = e(U) + 0.5*( dU'dU + dW'dW )
         LV2 = 0.0d0 ! Curvature   = e(B) + 0.5*( dU'dB + dB'dU )
         LV3 = 0.0d0 ! Quad.strain = 0.5*( dB'dB )

         DO p = 1,n
           DO i = 1,3
             DO j = 1,2
               dUdx(i,j) = dUdx(i,j) &
                   + LocalDeflection(6*(p-1)+i) * dBasisdx(p,j)

              dRdx(i,j) = dRdx(i,j) &
                  + LocalDeflection(6*(p-1)+i+3) * dBasisdx(p,j)
             END DO
           END DO
         END DO

! Linear part of strain and curvature
!------------------------------------
         DO p = 1,n
           LV1(1,1) = LV1(1,1) + LocalDeflection(6*(p-1)+1) * dBasisdx(p,1)
           LV1(2,1) = LV1(2,1) + LocalDeflection(6*(p-1)+2) * dBasisdx(p,2)
           LV1(3,1) = LV1(3,1) + LocalDeflection(6*(p-1)+1) * dBasisdx(p,2) & 
               + LocalDeflection(6*(p-1)+2) * dBasisdx(p,1) 

           LV2(1,1) = LV2(1,1) + LocalDeflection(6*(p-1)+4) * dBasisdx(p,1)
           LV2(2,1) = LV2(2,1) + LocalDeflection(6*(p-1)+5) * dBasisdx(p,2)
           LV2(3,1) = LV2(3,1) + LocalDeflection(6*(p-1)+4) * dBasisdx(p,2) & 
               + LocalDeflection(6*(p-1)+5) * dBasisdx(p,1)  
         END DO

! Non-linear part of strain and curvature
!----------------------------------------
         DO q = 1,2
           LV1(1,1) = LV1(1,1) + 0.5d0 * dUdx(q,1)**2
           LV1(2,1) = LV1(2,1) + 0.5d0 * dUdx(q,2)**2
           LV1(3,1) = LV1(3,1) + dUdx(q,1) * dUdx(q,2)

           LV2(1,1) = LV2(1,1) + dUdx(q,1) * dRdx(q,1)
           LV2(2,1) = LV2(2,1) + dUdx(q,2) * dRdx(q,2)
           LV2(3,1) = LV2(3,1) + dUdx(q,1) * dRdx(q,2) &
               + dUdx(q,2) * dRdx(q,1)
    
           LV3(1,1) = LV3(1,1) + 0.5d0 * dRdx(q,1)**2
           LV3(2,1) = LV3(2,1) + 0.5d0 * dRdx(q,2)**2
           LV3(3,1) = LV3(3,1) + dRdx(q,1) * dRdx(q,2)
         END DO

         LV1(1,1) = LV1(1,1) + 0.5d0 * dUdx(3,1)**2
         LV1(2,1) = LV1(2,1) + 0.5d0 * dUdx(3,2)**2
         LV1(3,1) = LV1(3,1) + dUdx(3,1) * dUdx(3,2)

         LV2(1,1) = LV2(1,1) + dUdx(3,1) * dRdx(3,1)
         LV2(2,1) = LV2(2,1) + dUdx(3,2) * dRdx(3,2)
         LV2(3,1) = LV2(3,1) + dUdx(3,1) * dRdx(3,2) &
               + dUdx(3,2) * dRdx(3,1)

         EPSILO(1) = LV1(1,1)
         EPSILO(2) = LV1(2,1)
         EPSILO(3) = LV1(3,1)

         KappaVector(1) = LV2(1,1)
         KappaVector(2) = LV2(2,1)
         KappaVector(3) = LV2(3,1)
       END IF

       CALL IsotropicElasticity( Dmatrix, Astarmatrix, NodalPoisson, &
                 NodalYoung, NodalThickness, Basis, n )
          
       Bmatrix = 0.0d0
          
       CALL IsotropicInPlaneElasticity( Amatrix, NodalPoisson, &
            NodalYoung, NodalThickness, Basis, n )
     
       Gdrilling(1,1) = StabParam2*(Astarmatrix(1,1)+Astarmatrix(2,2))

!      Normaalivoima- ja momenttivektorit (per pituusyksikk�) lokaalissa koord.:
!      -------------------------------------------------------------------------
       NVec = 0.0d0
       Mvec = 0.0d0
       NVec = MATMUL( Amatrix, EPSILO ) + MATMUL( Bmatrix, KappaVector )
       MVec = MATMUL( Bmatrix, EPSILO ) + MATMUL( Dmatrix, KappaVector )

!      Lopuksi vektorit tensoreiksi plus transformaatiot globaaliin koordinaatistoon:
!      ------------------------------------------------------------------------------
       Eps = 0.0d0
       Eps(1,1) = EPSILO(1)
       Eps(2,2) = EPSILO(2)
       Eps(1,2) = EPSILO(3) / 2.0d0
       Eps(2,1) = EPSILO(3) / 2.0d0

       Eps(1,3) = GammaVector(1) / 2.0d0
       Eps(2,3) = GammaVector(2) / 2.0d0
       Eps(3,1) = GammaVector(1) / 2.0d0
       Eps(3,2) = GammaVector(2) / 2.0d0

       Eps = MATMUL( TRANSPOSE(Transformation), Eps )
       Eps = MATMUL( Eps, Transformation )

       Kap = 0.0d0
       Kap(1,1) = KappaVector(1)
       Kap(2,2) = KappaVector(2)
       Kap(1,2) = KappaVector(3)/2.0d0
       Kap(2,1) = KappaVector(3)/2.0d0
       Kap = MATMUL( TRANSPOSE(Transformation), Kap )
       Kap = MATMUL( Kap, Transformation )

       NTen = 0.0d0
       NTen(1,1)=NVec(1)
       NTen(1,2)=NVec(3)
       NTen(2,1)=NVec(3)
       NTen(2,2)=NVec(2)

       NtenMaterial = 0.0d0
       NtenMaterial(1:2,1:2) = Nten(1:2,1:2)

       NTen = MATMUL( TRANSPOSE(Transformation), NTen )
       NTen = MATMUL( NTen, Transformation )

       MTen = 0.0d0
       MTen(1,1)=MVec(1)
       MTen(1,2)=MVec(3)
       MTen(2,1)=MVec(3)
       MTen(2,2)=MVec(2)

       MtenMaterial = 0.0d0
       MtenMaterial(1:2,1:2) = Mten(1:2,1:2)

       MTen = MATMUL( TRANSPOSE(Transformation), MTen )
       MTen = MATMUL( MTen, Transformation )

       SELECT CASE( NumberOfElementNodes )
       CASE( 3 )
          CALL AveragingWeights3( Nodes, Weight3 )
       CASE( 4 )
          CALL AveragingWeights4( Nodes, Weight4 )
       END SELECT


!      Restore the original node points:
!      ---------------------------------
       Nodes % x(1:n) = CopyOfNodes(1,1:n)
       Nodes % y(1:n) = CopyOfNodes(2,1:n)
       Nodes % z(1:n) = CopyOfNodes(3,1:n)

     END SUBROUTINE LocalStress
!==============================================================================


     SUBROUTINE AveragingWeights3( Nodes, Weight3 )
!------------------------------------------------------------------------------
       TYPE( Nodes_t ) :: Nodes
       REAL( KIND=DP ) :: Weight3(:)
!------------------------------------------------------------------------------
       INTEGER :: i
       REAL( KIND=DP ) :: Side1(2), Side2(2)
!------------------------------------------------------------------------------
       side1(1) = Nodes % x(2) - Nodes % x(1)
       side1(2) = Nodes % y(2) - Nodes % y(1)
       side1 = side1 / SQRT( side1(1)**2 + side1(2)**2 )

       side2(1) = Nodes % x(3) - Nodes % x(1)
       side2(2) = Nodes % y(3) - Nodes % y(1)
       side2 = side2 / SQRT( side2(1)**2 + side2(2)**2 )

       weight3(1) = SUM( side1 * side2 )

       side1(1) = Nodes % x(1) - Nodes % x(2)
       side1(2) = Nodes % y(1) - Nodes % y(2)
       side1 = side1 / SQRT( side1(1)**2 + side1(2)**2 )

       side2(1) = Nodes % x(3) - Nodes % x(2)
       side2(2) = Nodes % y(3) - Nodes % y(2)
       side2 = side2 / SQRT( side2(1)**2 + side2(2)**2 )

       weight3(2) = SUM( side1 * side2 )

       side1(1) = Nodes % x(1) - Nodes % x(3)
       side1(2) = Nodes % y(1) - Nodes % y(3)
       side1 = side1 / SQRT( side1(1)**2 + side1(2)**2 )

       side2(1) = Nodes % x(2) - Nodes % x(3)
       side2(2) = Nodes % y(2) - Nodes % y(3)
       side2 = side2 / SQRT( side2(1)**2 + side2(2)**2 )

       weight3(3) = SUM( side1 * side2 )

       DO i = 1,3
          weight3(i) = ACOS( weight3(i) ) 
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE AveragingWeights3


     SUBROUTINE AveragingWeights4( Nodes, Weight4 )
!------------------------------------------------------------------------------
       TYPE( Nodes_t ) :: Nodes
       REAL( KIND=DP ) :: Weight4(:)
!------------------------------------------------------------------------------
       INTEGER :: i
       REAL( KIND=DP ) :: Side1(2), Side2(2)
!------------------------------------------------------------------------------
       side1(1) = Nodes % x(2) - Nodes % x(1)
       side1(2) = Nodes % y(2) - Nodes % y(1)
       side1 = side1 / SQRT( side1(1)**2 + side1(2)**2 )

       side2(1) = Nodes % x(4) - Nodes % x(1)
       side2(2) = Nodes % y(4) - Nodes % y(1)
       side2 = side2 / SQRT( side2(1)**2 + side2(2)**2 )

       weight4(1) = SUM( side1 * side2 )

       side1(1) = Nodes % x(1) - Nodes % x(2)
       side1(2) = Nodes % y(1) - Nodes % y(2)
       side1 = side1 / SQRT( side1(1)**2 + side1(2)**2 )

       side2(1) = Nodes % x(3) - Nodes % x(2)
       side2(2) = Nodes % y(3) - Nodes % y(2)
       side2 = side2 / SQRT( side2(1)**2 + side2(2)**2 )

       weight4(2) = SUM( side1 * side2 )

       side1(1) = Nodes % x(4) - Nodes % x(3)
       side1(2) = Nodes % y(4) - Nodes % y(3)
       side1 = side1 / SQRT( side1(1)**2 + side1(2)**2 )

       side2(1) = Nodes % x(2) - Nodes % x(3)
       side2(2) = Nodes % y(2) - Nodes % y(3)
       side2 = side2 / SQRT( side2(1)**2 + side2(2)**2 )

       weight4(3) = SUM( side1 * side2 )

       side1(1) = Nodes % x(3) - Nodes % x(4)
       side1(2) = Nodes % y(3) - Nodes % y(4)
       side1 = side1 / SQRT( side1(1)**2 + side1(2)**2 )

       side2(1) = Nodes % x(1) - Nodes % x(4)
       side2(2) = Nodes % y(1) - Nodes % y(4)
       side2 = side2 / SQRT( side2(1)**2 + side2(2)**2 )

       weight4(4) = SUM( side1 * side2 )

       DO i = 1,4
          weight4(i) = ACOS( weight4(i) ) 
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE AveragingWeights4


!------------------------------------------------------------------------------
     SUBROUTINE SwitchToLocal( Nodes, CopyOfNodes, Transformation, n )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Transformation(:,:), CopyOfNodes(:,:)
       TYPE(Nodes_t) :: Nodes
       INTEGER :: n
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: XYZGlobal(3,4), XYZLocal(3,4)
!------------------------------------------------------------------------------
       CopyOfNodes(1,1:n) = Nodes % x(1:n)
       CopyOfNodes(2,1:n) = Nodes % y(1:n)
       CopyOfNodes(3,1:n) = Nodes % z(1:n)

       XYZGlobal(1,1:n)  = CopyOfNodes(1,1:n) - SUM( CopyOfNodes(1,1:n) ) / n
       XYZGlobal(2,1:n)  = CopyOfNodes(2,1:n) - SUM( CopyOfNodes(2,1:n) ) / n
       XYZGlobal(3,1:n)  = CopyOfNodes(3,1:n) - SUM( CopyOfNodes(3,1:n) ) / n

       XYZLocal(1:3,1:n) = MATMUL( Transformation(1:3,1:3), XYZGlobal(1:3,1:n) )

       Nodes % x(1:n) = XYZLocal(1,1:n)
       Nodes % y(1:n) = XYZLocal(2,1:n)
       Nodes % z(1:n) = XYZLocal(3,1:n)
!------------------------------------------------------------------------------
     END SUBROUTINE SwitchToLocal
!------------------------------------------------------------------------------

!==============================================================================

!------------------------------------------------------------------------------
     FUNCTION LocalBasis( Nodes, n ) RESULT( BasisVectors )
!------------------------------------------------------------------------------
       TYPE(Nodes_t) :: Nodes
       REAL(KIND=dp) :: BasisVectors(3,3)
       INTEGER :: n
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Tangent1(3), Tangent2(3)
!------------------------------------------------------------------------------

!      First, find a couple of in-plane unit vectors:
!      ----------------------------------------------
       Tangent1(1) = Nodes % x(2) - Nodes % x(1)
       Tangent1(2) = Nodes % y(2) - Nodes % y(1)
       Tangent1(3) = Nodes % z(2) - Nodes % z(1)
       Tangent1 = Tangent1 / SQRT( SUM( Tangent1**2  ) )

       Tangent2(1) = Nodes % x(3) - Nodes % x(2)
       Tangent2(2) = Nodes % y(3) - Nodes % y(2)
       Tangent2(3) = Nodes % z(3) - Nodes % z(2)
       Tangent2 = Tangent2 / SQRT( SUM( Tangent2**2  ) )

!      Then, define the local cartesian unit basis vectors:
!      ----------------------------------------------------
       BasisVectors(1:3,1) = Tangent1
       BasisVectors(1:3,2) = Tangent2 - SUM( Tangent1 * Tangent2 ) * Tangent1
       BasisVectors(1:3,2) = BasisVectors(1:3,2) / SQRT( SUM( BasisVectors(1:3,2)**2 ) )
       BasisVectors(1:3,3) = CrossProductL( BasisVectors(1:3,1), BasisVectors(1:3,2) )
!------------------------------------------------------------------------------       
     END FUNCTION LocalBasis
!------------------------------------------------------------------------------       


!------------------------------------------------------------------------------       
     SUBROUTINE IsotropicElasticity(Ematrix, &
          Gmatrix,Poisson,Young,Thickness,Basis,n)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Ematrix(:,:), Gmatrix(:,:), Basis(:)
     REAL(KIND=dp) :: Poisson(:), Young(:), Thickness(:)
     REAL(KIND=dp) :: Euvw, Puvw, Guvw, Tuvw
     INTEGER :: n
!------------------------------------------------------------------------------
       Euvw = SUM( Young(1:n)    * Basis(1:n) )
       Puvw = SUM( Poisson(1:n)  * Basis(1:n) )
       Tuvw = SUM( Thickness(1:n)* Basis(1:n) )
       Guvw = Euvw/(2.0d0*(1.0d0 + Puvw))

       Ematrix = 0.0d0
       Ematrix(1,1) = 1.0d0
       Ematrix(1,2) = Puvw
       Ematrix(2,1) = Puvw
       Ematrix(2,2) = 1.0d0
       Ematrix(3,3) = (1.0d0-Puvw)/2.0d0
       Ematrix = Ematrix * Euvw * (Tuvw**3) / (12.0d0 * (1.0d0 - Puvw*Puvw))

       Gmatrix = 0.0d0
       Gmatrix(1,1) = Guvw*Tuvw
       Gmatrix(2,2) = Guvw*Tuvw
!------------------------------------------------------------------------------
     END SUBROUTINE IsotropicElasticity
!------------------------------------------------------------------------------

!==============================================================================

!------------------------------------------------------------------------------
     SUBROUTINE IsotropicInPlaneElasticity( Ematrix, &
          Poisson, Young, Thickness, Basis, n )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Ematrix(:,:), Basis(:)
     REAL(KIND=dp) :: Poisson(:), Young(:), Thickness(:)
     REAL(KIND=dp) :: Euvw, Puvw, Guvw, Tuvw
     INTEGER :: n
!------------------------------------------------------------------------------
       Euvw = SUM( Young(1:n)    * Basis(1:n) )
       Puvw = SUM( Poisson(1:n)  * Basis(1:n) )
       Tuvw = SUM( Thickness(1:n)* Basis(1:n) )
       Guvw = Euvw/(2.0d0*(1.0d0 + Puvw))

       Ematrix = 0.0d0
       Ematrix(1,1) = 1.0d0
       Ematrix(1,2) = Puvw
       Ematrix(2,1) = Puvw
       Ematrix(2,2) = 1.0d0
       Ematrix(3,3) = (1.0d0-Puvw)/2.0d0
       Ematrix = Ematrix * Tuvw * Euvw / (1.0d0 - Puvw*Puvw)
!------------------------------------------------------------------------------
     END SUBROUTINE IsotropicInPlaneElasticity
!------------------------------------------------------------------------------

!==============================================================================

!------------------------------------------------------------------------------
     SUBROUTINE ShearCorrectionFactor(Kappa,Thickness,x,y,n,StabParam)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Kappa,Thickness,x(:),y(:),StabParam
       INTEGER :: n
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: x21,x32,x43,x13,x14,y21,y32,y43,y13,y14, &
            l21,l32,l43,l13,l14,alpha,h
!------------------------------------------------------------------------------
       Kappa = 1.0d0
       SELECT CASE(n)
          CASE(3)
             alpha = 0.20d0 * StabParam
             x21 = x(2)-x(1)
             x32 = x(3)-x(2)
             x13 = x(1)-x(1)
             y21 = y(2)-y(1)
             y32 = y(3)-y(2)
             y13 = y(1)-y(1)
             l21 = SQRT(x21**2 + y21**2)
             l32 = SQRT(x32**2 + y32**2)
             l13 = SQRT(x13**2 + y13**2)
             h = MAX(l21,l32,l13)
             Kappa = (Thickness**2)/(Thickness**2 + alpha*(h**2))
          CASE(4)
             alpha = 0.10d0 * StabParam
             x21 = x(2)-x(1)
             x32 = x(3)-x(2)
             x43 = x(4)-x(3)
             x14 = x(1)-x(4)
             y21 = y(2)-y(1)
             y32 = y(3)-y(2)
             y43 = y(4)-y(3)
             y14 = y(1)-y(4)
             l21 = SQRT(x21**2 + y21**2)
             l32 = SQRT(x32**2 + y32**2)
             l43 = SQRT(x43**2 + y43**2)
             l14 = SQRT(x14**2 + y14**2)
             h = MAX(l21,l32,l43,l14)
             Kappa = (Thickness**2)/(Thickness**2 + alpha*(h**2))
           CASE DEFAULT
             CALL Fatal('ShellSolver','Illegal number of nodes for Smitc elements')
           END SELECT
!------------------------------------------------------------------------------
     END SUBROUTINE ShearCorrectionFactor
!------------------------------------------------------------------------------

!==============================================================================

!------------------------------------------------------------------------------
     SUBROUTINE AddEnergy(A,B,C,m,n,s)
!------------------------------------------------------------------------------
!      Performs the operation
!
!         A = A + C' * B * C * s
!
!      with
!
!         Size( A ) = n x n
!         Size( B ) = m x m
!         Size( C ) = m x n
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: A(:,:),B(:,:),C(:,:),s
       INTEGER :: m,n
!------------------------------------------------------------------------------
       INTEGER :: i,j,k,l
!------------------------------------------------------------------------------
       DO i=1,n
          DO j=1,n
             DO k=1,m
                DO l=1,m
                   A(i,j) = A(i,j) + C(k,i)*B(k,l)*C(l,j) * s
                END DO
             END DO
          END DO
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE AddEnergy
!------------------------------------------------------------------------------

!==============================================================================

!------------------------------------------------------------------------------
     SUBROUTINE AddInnerProducts(A,B,C,D,m,n,s)
!------------------------------------------------------------------------------
!      Performs the operation
!
!         A = A + C * B * D * s
!
!      with
!
!         Size( A ) = n x n
!         Size( B ) = m x m
!         Size( C ) = n x m
!         Size( D ) = m x n
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: A(:,:),B(:,:),C(:,:),D(:,:),s
       INTEGER :: m,n
!------------------------------------------------------------------------------
       INTEGER :: i,j,k,l
!------------------------------------------------------------------------------
       DO i=1,n
          DO j=1,n
             DO k=1,m
                DO l=1,m
                   A(i,j) = A(i,j) + C(i,k)*B(k,l)*D(l,j) * s
                END DO
             END DO
          END DO
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE AddInnerProducts
!------------------------------------------------------------------------------

!==============================================================================

!------------------------------------------------------------------------------
     SUBROUTINE CovariantInterpolation(ShearStrain,Basis,X,Y,U,V,n)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: ShearStrain(:,:),Basis(:),X(:),Y(:),U,V
       INTEGER :: n
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: detJ,Jmat(2,2),invJ(2,2),ShearRef(2,100)
       REAL(KIND=dp) :: Tau(2),Sdofs(100)
       INTEGER :: j

       SELECT CASE(n)

!      The SMITC3 element
!      ==================
       CASE(3)
          CALL Jacobi3(Jmat,invJ,detJ,x,y)
          ShearRef = 0.0d0
          ShearStrain = 0.0d0

!         Compute the shear-dofs for edge 12:
!         ===================================
          Tau = (/ 1.0d0, 0.0d0/)

          Sdofs = 0.0d0
          Sdofs(4) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/2.0d0
          Sdofs(5) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/2.0d0
          Sdofs(10) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/2.0d0
          Sdofs(11) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/2.0d0

          DO j = 1,18
             ShearRef(1,j) = ShearRef(1,j) + (1+V)*Sdofs(j)
             ShearRef(2,j) = ShearRef(2,j) + ( -U)*Sdofs(j)
          END DO

!         Compute the shear-dofs for edge 23:
!         ===================================
          Tau(1) = -1.0d0/SQRT(2.0d0)
          Tau(2) =  1.0d0/SQRT(2.0d0)

          Sdofs = 0.0d0
          Sdofs(10) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/SQRT(2.0d0)
          Sdofs(11) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/SQRT(2.0d0)
          Sdofs(16) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/SQRT(2.0d0)
          Sdofs(17) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/SQRT(2.0d0)

          DO j = 1,18
             ShearRef(1,j) = ShearRef(1,j) + ( V)*Sdofs(j)
             ShearRef(2,j) = ShearRef(2,j) + (-U)*Sdofs(j)
          END DO

!         Compute the shear-dofs for edge 31:
!         ===================================
          Tau(1) =  0.0d0
          Tau(2) = -1.0d0

          Sdofs = 0.0d0
          Sdofs(4) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/2.0d0
          Sdofs(5) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/2.0d0
          Sdofs(16) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/2.0d0
          Sdofs(17) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/2.0d0

          DO j = 1,18
             ShearRef(1,j) = ShearRef(1,j) + (  V )*Sdofs(j)
             ShearRef(2,j) = ShearRef(2,j) + (-1-U)*Sdofs(j)
          END DO

!         Compute the final reduced shear strain
!         ======================================
          ShearStrain(1:2,1:18) = MATMUL(invJ,ShearRef(1:2,1:18))

!      The SMITC4 element
!      ==================
       CASE(4)
          ShearRef = 0.0d0
          ShearStrain = 0.0d0

!         Compute the shear-dofs for edge 12:
!         ===================================
          Tau(1) = 1.0d0
          Tau(2) = 0.0d0

          CALL Jacobi4(Jmat,invJ,detJ,0.0d0,-1.0d0,x,y)
          
          Sdofs = 0.0d0
          Sdofs(4) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(5) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))
          Sdofs(10) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(11) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))

          DO j = 1,24
             ShearRef(1,j) = ShearRef(1,j) + (1-V)/4.0d0*Sdofs(j)
          END DO

!         Compute the shear-dofs for edge 23:
!         ===================================
          Tau(1) = 0.0d0
          Tau(2) = 1.0d0

          CALL Jacobi4(Jmat,invJ,detJ,1.0d0,0.0d0,x,y)

          Sdofs = 0.0d0
          Sdofs(10) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(11) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))
          Sdofs(16) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(17) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))

          DO j = 1,24
             ShearRef(2,j) = ShearRef(2,j) + (1+U)/4.0d0*Sdofs(j)
          END DO

!         Compute the shear-dofs for edge 34:
!         ===================================
          Tau(1) = -1.0d0
          Tau(2) =  0.0d0

          CALL Jacobi4(Jmat,invJ,detJ,0.0d0,1.0d0,x,y)

          Sdofs = 0.0d0
          Sdofs(16)  = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(17)  = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))
          Sdofs(22) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(23) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))

          DO j = 1,24
             ShearRef(1,j) = ShearRef(1,j) + (-1-V)/4.0d0*Sdofs(j)
          END DO

!         Compute the shear-dofs for edge 41:
!         ===================================
          Tau(1) =  0.0d0
          Tau(2) = -1.0d0

          CALL Jacobi4(Jmat,invJ,detJ,-1.0d0,0.0d0,x,y)

          Sdofs = 0.0d0
          Sdofs(4)  = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(5)  = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))
          Sdofs(22) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))
          Sdofs(23) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))

          DO j = 1,24
             ShearRef(2,j) = ShearRef(2,j) + (-1+U)/4.0d0*Sdofs(j)
          END DO

!         Compute the final reduced shear strain
!         ======================================
          CALL Jacobi4(Jmat,invJ,detJ,U,V,x,y)
          ShearStrain(1:2,1:24) = MATMUL(invJ,ShearRef(1:2,1:24))

       CASE DEFAULT
         CALL Fatal('ShellSolver','Illegal number of nodes for Smitc elements.')
       END SELECT
!------------------------------------------------------------------------------
     END SUBROUTINE CovariantInterpolation
!------------------------------------------------------------------------------

!==============================================================================

!------------------------------------------------------------------------------
     SUBROUTINE Jacobi3(Jmat,invJ,detJ,x,y)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Jmat(:,:),invJ(:,:),detJ,x(:),y(:)
!------------------------------------------------------------------------------
       Jmat(1,1) = x(2)-x(1)
       Jmat(2,1) = x(3)-x(1)
       Jmat(1,2) = y(2)-y(1)
       Jmat(2,2) = y(3)-y(1)

       detJ = Jmat(1,1)*Jmat(2,2)-Jmat(1,2)*Jmat(2,1)

       invJ(1,1) =  Jmat(2,2)/detJ
       invJ(2,2) =  Jmat(1,1)/detJ
       invJ(1,2) = -Jmat(1,2)/detJ
       invJ(2,1) = -Jmat(2,1)/detJ
!------------------------------------------------------------------------------
     END SUBROUTINE Jacobi3
!------------------------------------------------------------------------------

!==============================================================================

!------------------------------------------------------------------------------
     SUBROUTINE Jacobi4(Jmat,invJ,detJ,xi,eta,x,y)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Jmat(:,:),invJ(:,:),detJ,xi,eta,x(:),y(:)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: dNdxi(4), dNdeta(4)
       INTEGER :: i

       dNdxi(1) = -(1-eta)/4.0d0
       dNdxi(2) =  (1-eta)/4.0d0
       dNdxi(3) =  (1+eta)/4.0d0
       dNdxi(4) = -(1+eta)/4.0d0
       dNdeta(1) = -(1-xi)/4.0d0
       dNdeta(2) = -(1+xi)/4.0d0
       dNdeta(3) =  (1+xi)/4.0d0
       dNdeta(4) =  (1-xi)/4.0d0
       
       Jmat = 0.0d0
       DO i=1,4
          Jmat(1,1) = Jmat(1,1) + dNdxi(i)*x(i)
          Jmat(1,2) = Jmat(1,2) + dNdxi(i)*y(i)
          Jmat(2,1) = Jmat(2,1) + dNdeta(i)*x(i)
          Jmat(2,2) = Jmat(2,2) + dNdeta(i)*y(i)
       END DO

       detJ = Jmat(1,1)*Jmat(2,2)-Jmat(1,2)*Jmat(2,1)

       invJ(1,1) = Jmat(2,2)/detJ
       invJ(2,2) = Jmat(1,1)/detJ
       invJ(1,2) = -Jmat(1,2)/detJ
       invJ(2,1) = -Jmat(2,1)/detJ
!------------------------------------------------------------------------------
     END SUBROUTINE Jacobi4
!------------------------------------------------------------------------------

!==============================================================================

!------------------------------------------------------------------------------
     FUNCTION CrossProductL( v1, v2 ) RESULT( v3 )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: v1(3), v2(3), v3(3)
       v3(1) =  v1(2)*v2(3) - v1(3)*v2(2)
       v3(2) = -v1(1)*v2(3) + v1(3)*v2(1)
       v3(3) =  v1(1)*v2(2) - v1(2)*v2(1)
!------------------------------------------------------------------------------
     END FUNCTION CrossProductL
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE ShellSolver
!------------------------------------------------------------------------------
