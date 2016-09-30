 ! /*****************************************************/
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
! **********************************************************/
! *
! * THE REISSNER-MINDLIN FACET SHELL SOLVER MODULE FOR COMPOSITE LAMINATES
! *
! * Author: Mikko Lyly
! * Address: CSC - IT Center for Science Ltd.
! * Keilaranta 14, P.O. BOX 405
! * 02101 Espoo, Finland
! * EMail: Mikko.Lyly@csc.fi
! *
! * Date: 09 Apr 2000
! * Modified by: Mikko Lyly & Petri Kere, in 2001
! * Date of modification: 29.5.2001, 21.8.2001
! *
! * 
! * CHECKED AND MODIFIED SHELL FORMULATION AND ADDED DKT (Discrete Kirchhoff Theory),
! * Membrane only, MITC3 and RMITC3 TRIANGULAR ELEMENT - DR. O.P.GUPTA - (IN 2014-15)
! *		
! * Dr.O.P.Gupta
! * Professor (Retired) Mechanical Engineering Department,
! * Indian Institute of Technology,
! * Kharagpur, India 
! * opgupta112000@yahoo.com 
! *
! ******************************************************/
! ***************************************************************
! * This version uses definitions of BetaX as -du/dz and BetaY
! * as -dv/dz (defined in paper 'Note on MITC3....by OPGupta').
! ***************************************************************
! -----------------------------------------------------
  SUBROUTINE ShellSolver_Init( Model,Solver,dt,Transient )
! ----------------------------------------------------
    USE DefUtils

    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: DT
    LOGICAL :: Transient
! ------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: SolverParams
! ------------------------------------------------------------
    SolverParams => GetSolverParams()
    IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      CALL ListAddInteger( SolverParams, 'Variable DOFs', 6 )
      CALL ListAddString( SolverParams, 'Variable', 'Deflection' )
    END IF
    CALL ListAddInteger( SolverParams, 'Time derivative order', 2 )
! -------------------------------------------------------------
  END SUBROUTINE ShellSolver_Init
! --------------------------------------------------------------

 
! ------------------------------------------------------
   SUBROUTINE ShellSolver( Model,Solver,dt,TransientSimulation )
! -----------------------------------------------------------------
! *****************************************************
!
!  Solve the Reissner-Mindlin facet shell equations!
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
! *******************************************************

     USE DefUtils
     IMPLICIT NONE
! --------------------------------------------------------------
     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET :: Solver
 
     REAL(KIND=DP) :: dt
     LOGICAL :: TransientSimulation

     TYPE(Solver_t), POINTER :: PSolver
     TYPE(Matrix_t),POINTER  :: StiffMatrix
     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement

     INTEGER :: i,j,k,n,t, bf_id, istat, LocalNodes, CalcSurf,nPL,iPL, NumberOfElementNodes

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
         EpsxyNodal(:), EpsxzNodal(:), EpsyzNodal(:), SEqNodal(:)

     REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:),LoadX(:),LoadY(:), LoadZ(:), LoadN(:), FORCE(:),Poisson(:), Thickness(:),& 
         Young(:), Tension(:),MASS(:,:), DAMP(:,:),  Density(:),LocalDeflection(:), SxxElement(:), & 
         SyyElement(:), SzzElement(:), &
         SxyElement(:), SxzElement(:), SyzElement(:), LoadVector(:,:), &
         Weights(:,:), Referenced(:), &
         EpsxxElement(:), EpsyyElement(:), EpszzElement(:), &
         EpsxyElement(:), EpsxzElement(:), EpsyzElement(:), &
         ThLoad(:), TT(:),BendingLoad(:)
     REAL(KIND=dp) :: at,st,at0

     REAL(KIND=dp) :: StabParam1, StabParam2, S(3,3), PlyThick, &
         Tx(3), Ty(3), Tz(3), l, Weight3(3), Weight4(4), &
         Fx, Fy, Fz, Mx, My, Mz, &
         Eps(3,3), Kap(3,3), Nten(3,3), Mten(3,3),Amatrix(3,3), Bmatrix(3,3), Dmatrix(3,3), Astarmatrix(2,2),Gdrilling(1,1),& 
         Transformation(3,3), T5(5,5), T0(3,3), T0S(2,2), &
         LamRF, Nvector(3), NtenMaterial(2,2), MtenMaterial(2,2), ss

     LOGICAL :: GotIt, GotForceBC

     LOGICAL :: Isotropic = .FALSE.
     LOGICAL :: WeldSimu, ThStrain,ThermalAnisotropy, AnisoPlane, &
         UseDKTtriangle, UseRDKTtriangle,  UseSMITCelement, &
         UseMITC3element, MembraneOnly, &	
         MembraneOnlyNRM,  &	
         TopSideStress, BottomSideStress, TotalStress, MembraneStress, BendingStress,&	
         SimpleObject,UseRMITC3element
     CHARACTER (LEN=30) :: ObjectType		
     TYPE(ValueList_t), POINTER :: SolverParams, Material,BodyForce, BC

     REAL(KIND=dp) :: xn, yn, zn, xp, yp, zp, r, PressureType, dummy1, OutNormal,&	
         ObjectRadius				
     REAL(KIND=dp) :: Cscalar, Gscalar, Lambda, CalLambda,ArcLength, ArcLengthScale, Radius, &
         MaximumDeflection, DisplMax, Maxcomponent
     INTEGER :: ArcLengthSteps, ArcLengthStep

     SAVE STIFF, MASS, LoadX, LoadY, LoadZ, LoadN, &
         FORCE, ElementNodes, Poisson, Density, Young, Thickness, &
         Tension, AllocationsDone, DAMP, LocalDeflection, &
         SxxNodal, SyyNodal, SzzNodal, SxyNodal, SxzNodal, SyzNodal, &
         EpsxxNodal, EpsyyNodal, EpszzNodal, EpsxyNodal, EpsxzNodal, &
         EpsyzNodal, LoadVector, SEqNodal	
     SAVE WeldSimu, ThStrain,ThermalAnisotropy, ThLoad,TT,AnisoPlane,BendingLoad, &
         UseDKTtriangle, ObjectType,PressureType, UseRDKTtriangle,  &
         UseSMITCelement,UseMITC3element, MembraneOnly, MembraneOnlyNRM,   &	
         TopSideStress, BottomSideStress, TotalStress, MembraneStress, BendingStress, OutNormal,&
         ObjectRadius,SimpleObject,UseRMITC3element
! ==============================================================


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
      
! =========================================================

!    Allocate some permanent storage, this is done first time only:
!    --------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = Model % MaxElementNodes
       
       ALLOCATE( ElementNodes % x( N ),   &
           ElementNodes % y( N ),   &
           ElementNodes % z( N ),   &
           FORCE( 6*N ), &
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
           ThLoad(6*N), TT(N), &
           BendingLoad(6*N),&	
           STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'ShellSolver',  'Memory allocation error,Aborting.' )
       END IF
       
       NULLIFY( PointLoad )

       NULLIFY( SxxNodal, SyyNodal, SzzNodal, &
           SxyNodal, SxzNodal, SyzNodal, &
           EpsxxNodal, EpsyyNodal, EpszzNodal, &
           EpsxyNodal, EpsxzNodal, EpsyzNodal )
       
       AllocationsDone = .TRUE.
     END IF

!=============================================================

!    Do some additional initialization, and go for it:
!    -------------------------------------------------
     at  = CPUTime()
     at0 = RealTime()
     Isotropic = .TRUE.
     
! ------------------------------------------------------------
!     Non-linear iteration etc.
! ------------------------------------------------------
     SolverParams => GetSolverParams()

     NonLinConvTol = GetConstReal( SolverParams, &
         'Nonlinear System Convergence Tolerance', GotIt )
     IF( .NOT. GotIt ) NonLinConvTol = 1.0e-6
     
     NonLinMaxIt = GetInteger( SolverParams, &
         'Nonlinear System Max Iterations', GotIt )
     IF( .NOT. GotIt ) NonLinMaxIt = 1
     
     LargeDeflection =  GetLogical( SolverParams, 'Large Deflection', GotIt )
     IF( .NOT. GotIt )THEN
       LargeDeflection = .FALSE.
     ELSE
       CALL Info( 'ShellSolver',  '******* Large Deflection option is currently not available ******' )
       LargeDeflection = .FALSE.
     ENDIF
     
     CALL Info( 'ShellSolver',  '**************************************************************',Level=1 )
     CALL Info( 'ShellSolver',  'In this program several triangular element formulations are used.',Level=1 )
     CALL Info( 'ShellSolver',  'These are Membrane only, DKT, MITC3, SMITC and RMITC3. Last three are Mixed  ',Level=1 )
     CALL Info( 'ShellSolver',  'formulation of MITC type. The default is RMITC3 type. It gives more accurate  ',Level=1)
     CALL Info( 'ShellSolver',  'analysis of deflection in shell structures. If the traditional Mixed formulation',Level=1)
     CALL Info( 'ShellSolver',  'or other formulations are to be used, specify  suitable option in .sif file. ',Level=1 )
     CALL Info( 'ShellSolver',  'Simulation of welding distortion is also possible through this solver',Level=1 )
     CALL Info( 'ShellSolver',  'Dr.O.P.Gupta, Ex-Professor,IIT Kharagpur, India - opgupta112000@yahoo.com',Level=1 )
     CALL Info( 'ShellSolver',  '**************************************************************',Level=1 )

     IF( .NOT. LargeDeflection ) THEN
       NonLinMaxIt = 1
       StressComputation = GetLogical( SolverParams, 'StressComputation', GotIt )
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
       CALL Info( 'ShellSolver', 'Shear stabilization parameter undefined.' )
       CALL Info( 'ShellSolver', 'Using default value 0.0 i.e. no shear stabilization' )
       StabParam1 = 0.0d0
     END IF

     StabParam2 = GetConstReal( SolverParams, &
         'Drilling Stabilization Parameter', GotIt )
     IF( .NOT.GotIt ) THEN
       CALL Info( 'ShellSolver', 'Drilling stabilization parameter undefined' )
       CALL Info( 'ShellSolver', 'Using default value 1.0' )
       StabParam2 = 1.0d0
     END IF

     SimpleObject = GetLogical( SolverParams, 'Object Of Simple Regular Shape', GotIt )
     TopSideStress = GetLogical( SolverParams, 'Top Side Stress', GotIt )
     BottomSideStress = GetLogical( SolverParams, 'Bottom Side Stress', GotIt )

     IF( TopSideStress .AND. BottomSideStress )THEN
       CALL Warn( 'ShellSolver', 'Both Top Side Stress and Bottom Side Stress can not be output simultaneously')
       CALL Fatal( 'ShellSolver', 'STOP')
     ENDIF
     
     TotalStress = GetLogical( SolverParams, 'Compute Total Stress', GotIt )
     MembraneStress = GetLogical( SolverParams, 'Compute Membrane Stress', GotIt )
     BendingStress = GetLogical( SolverParams, 'Compute Bending Stress', GotIt )

     IF( TotalStress .AND. (MembraneStress .OR. BendingStress ) )THEN
       CALL Fatal( 'ShellSolver', 'Only one of the three should be specified - Compute Membrane Stress or'// &
           ' Compute Bending Stress or Compute Total Stress')
     ENDIF
     
     IF((MembraneStress).AND.((TotalStress).OR.(BendingStress)))THEN
       CALL Fatal( 'ShellSolver', 'Only one of the three should be specified - Compute Membrane Stress or'// &
           ' Compute Bending Stress or Compute Total Stress')
     endif

     MembraneOnly = GetLogical( SolverParams, 'Membrane Only', GotIt )	
     MembraneOnlyNRM = GetLogical( SolverParams, 'Membrane Only No RM', GotIt )

     UseDKTtriangle = GetLogical( SolverParams, 'Use DKT Triangle', GotIt )
     UseRDKTtriangle = GetLogical(SolverParams, 'Use RDKT Triangle', GotIt )
     UseSMITCelement = GetLogical( SolverParams, 'Use SMITC Element', GotIt )
     UseMITC3element = GetLogical( SolverParams, 'Use MITC3 Element', GotIt )
     UseRMITC3element = GetLogical( SolverParams, 'Use RMITC3 Element', GotIt )
     
     CALL Info( 'ShellSolver', '***********************************************************' )
     IF( MembraneOnly )THEN
       CALL Info( 'ShellSolver', 'Analysis is now using Membrane only formulation' )
       CALL Info( 'ShellSolver', 'In-plane rotation minimization is used in this element' )
       StabParam1 = 0.0d0

     ELSEIF( MembraneOnlyNRM ) THEN
       CALL Info( 'ShellSolver', 'Analysis is now using Membrane only No RM formulation' )
       CALL Info( 'ShellSolver', 'Not using in-plane rotation minimization ' )
       StabParam1 = 0.0d0
       
     ELSEIF( UseDKTtriangle ) THEN
       CALL Info( 'ShellSolver', 'Analysis is now using DKT triangle formulation' )
       CALL Info( 'ShellSolver', 'Shear stabilization is not required with DKT trangular element' )
       CALL Info( 'ShellSolver', 'Taking Shear Stabilization Parameter zero.' )
       StabParam1 = 0.0d0

     ELSEIF( UseRDKTtriangle ) THEN
       CALL Info( 'ShellSolver', 'Analysis is now using RDKT triangle formulation' )
       CALL Info( 'ShellSolver', 'Shear stabilization is not required with RDKT trangular element' )
       CALL Info( 'ShellSolver', 'Taking Shear Stabilization Parameter zero.' )
       StabParam1 = 0.0d0

     ELSEIF( UseSMITCelement ) THEN
       CALL Info( 'ShellSolver', 'Analysis is now using SMITC element formulation' )
       CALL Info( 'ShellSolver', 'Shear stabilization is under test' )
       CALL Info( 'ShellSolver', 'Taking Shear Stabilization Parameter zero.' )

     ELSEIF( UseMITC3element ) THEN
       CALL Info( 'ShellSolver', 'Analysis is now using MITC3 (Lee) Formulation for the element ' )
       CALL Info( 'ShellSolver', 'Shear stabilization is under test' )
       CALL Info( 'ShellSolver', 'Taking Shear Stabilization Parameter zero.' )

     ELSEIF( UseRMITC3element ) THEN
       CALL Info( 'ShellSolver', 'Analysis is now using RMITC3 Formulation for the element ' )
       CALL Info( 'ShellSolver', 'Shear stabilization is under test' )
       CALL Info( 'ShellSolver', 'Taking Shear Stabilization Parameter zero.' )
       
     ELSE
       UseRMITC3element = .TRUE.
       CALL Info( 'ShellSolver', 'Analysis is now using RMITC3 Formulation for the element ' )
       CALL Info( 'ShellSolver', 'It is the default choice.' )
       CALL Info( 'ShellSolver', 'Taking Shear Stabilization Parameter zero.' )
     END IF
     CALL Info( 'ShellSolver', '************************************************************' )       

     IF((.NOT.UseDKTtriangle).AND.(.NOT.UseRDKTtriangle).AND. &
         (.NOT.UseSMITCelement).AND.(.NOT.UseMITC3element).AND. &
         (.NOT.MembraneOnlyNRM).AND.(.NOT.UseRMITC3element))THEN
       dummy1 = 0.0
     ELSEIF((.NOT.UseDKTtriangle).AND.(.NOT.UseRDKTtriangle).AND. &
         (.NOT.UseSMITCelement).AND.(.NOT.UseMITC3element).AND. &
         (.NOT.MembraneOnly).AND.(.NOT.UseRMITC3element))THEN
       dummy1 = 2.0
     ELSEIF((.NOT.MembraneOnly).AND.(.NOT.UseRDKTtriangle) &
         .AND.(.NOT.UseSMITCelement).AND.(.NOT.UseMITC3element).AND. &
         (.NOT.MembraneOnlyNRM).AND.(.NOT.UseRMITC3element))THEN
       dummy1 = 3.0
     ELSEIF((.NOT.MembraneOnly).AND.(.NOT.UseDKTtriangle) &
         .AND.(.NOT.UseSMITCelement).AND.(.NOT.UseMITC3element).AND. &
         (.NOT.MembraneOnlyNRM).AND.(.NOT.UseRMITC3element))THEN
       dummy1 = 4.0
     ELSEIF((.NOT.MembraneOnly).AND.(.NOT.UseDKTtriangle) &
         .AND.(.NOT.UseRDKTtriangle).AND.(.NOT.UseMITC3element).AND. &
         (.NOT.MembraneOnlyNRM).AND.(.NOT.UseRMITC3element))THEN
       dummy1 = 5.0
     ELSEIF((.NOT.MembraneOnly).AND.(.NOT.UseDKTtriangle) &
         .AND.(.NOT.UseRDKTtriangle).AND.(.NOT.UseSMITCelement).AND. &
         (.NOT.MembraneOnlyNRM).AND.(.NOT.UseRMITC3element))THEN
       dummy1 = 6.0
     ELSEIF((.NOT.MembraneOnly).AND.(.NOT.UseDKTtriangle) &
         .AND.(.NOT.UseRDKTtriangle).AND.(.NOT.UseMITC3element).AND.(.NOT.UseSMITCelement).AND. &
         (.NOT.MembraneOnlyNRM))THEN
       dummy1 = 7.0
     ELSE
       CALL Info( 'ShellSolver', '***********************************************************' )
       CALL Info( 'ShellSolver', 'For proper choice of element type, it is necessary that only one of the options ' )
       CALL Info( 'ShellSolver', 'or none  (default RMITC3) is specified.' )
       CALL Info( 'ShellSolver', '************************************************************' )
       CALL Fatal( 'ShellSolver', 'STOP')       
     ENDIF

     WRITE(Message,'(a,F8.2)') 'dummy1 = ',dummy1
     CALL Info('ShellSolve',Message)

!      Get element local matrix, and rhs vector:
!      -----------------------------------------
     Nvector = 0.0d0
! -------------------------------------------------------------
!     Non-linear iteration starts here
! ------------------------------------------------------------

     CALL SolveNonLinear()
     
     DisplMax = 0.0d0
     DO i = 1,Model % NumberOfNodes
       IF( ABS( Solver % Mesh % Nodes % x(i) - 0.05d0 ) < 1.0d-6  .AND.ABS( Solver % Mesh % Nodes % y(i) - 0.10d0 ) < 1.0d-6 & 
           .AND. ABS( Solver % Mesh % nodes % z(i) - 1.00d0 ) < 1.0d-6 )THEN
         j = DeflectionPerm(i)
         IF(j < 1) CYCLE
         DisplMax = Deflection(6*(j-1)+2)
       END IF
     END DO
! ----------------------------------------------------------------
 
   CONTAINS

! -----------------------------------------------------------------
     SUBROUTINE SolveNonLinear
! -----------------------------------------------------------------
       DO NonLinIter = 1,NonLinMaxIt
         CALL Info('ShellSolver','***********************************************')
         WRITE( Message,'(A,I4)') 'Newton iteration',nonliniter
         CALL Info('ShellSolver',Message)
         CALL Info('ShellSolver','**********************************************')
         !	*****Here i and j are used as integer without defining these beforehand. May be that 
         !	the 6 vowels are default integers (25.6.15)**********************************
         OPEN(12,FILE='CStrain_Cyl.txt') 
         
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
         CALL BCAssembly()				
         CALL DefaultFinishAssembly()	
         CALL ConcentratedLoads()		
         ForceVector = LoadScale*ForceVector	
         
! -----------------------------------------------------------------


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
         WRITE(Message,'(a)') ' Just before solution'
         CALL Info('ShellSolve',Message,Level=10)

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
         
         !     PRINT *,'Max deflection =',Unorm
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
       CLOSE(12)
       StressComputation = GetLogical( SolverParams, 'Stress Computation', GotIt )
       StressComputation = StressComputation .OR. &
           GetLogical( SolverParams, 'Calculate Stresses', GotIt )
       
       IF( StressComputation .AND. Solver % NOFEigenValues <= 0  ) THEN
         CALL Info('ShellSolve','Entering stress calculation routines...')
         CALL CalculateStresses()
       END IF ! If Stress Calculation
! ----------------------------------------------------------------
     END SUBROUTINE SolveNonLinear
! ---------------------------------------------------------------


! -------------------------------------------------------------
     SUBROUTINE BulkAssembly
! --------------------------------------------------------
       CALL StartAdvanceOutput('ShellSolve', 'Assembly:')
       DO t=1,Solver % NumberOfActiveElements
! --------------------------------------------------------------
         
         CALL AdvanceOutput(t,Solver % NumberOFActiveElements)
         
! ----------------------------------------------------------
         CurrentElement => GetActiveElement( t )
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
         BodyForce => GetBodyForce()
         
         LoadVector = 0.0d0
         
         LoadX(1:n) = GetReal( BodyForce, 'Body Force 1', GotIt )
         LoadY(1:n) = GetReal( BodyForce, 'Body Force 2', GotIt )
         LoadZ(1:n) = GetReal( BodyForce, 'Body Force 3', GotIt )
         LoadN(1:n) = GetReal( BodyForce, 'Pressure', GotIt )
         LoadN(1:n) = LoadN(1:n) + GetReal( BodyForce, 'Normal Pressure', GotIt )
!      Material data:
!      --------------
!	******************************************************************************
!	(14.10.15) material data and thickness may be shifted to Body Force section
!	for more versatility (check feasibility)
!	******************************************************************************
         Material => GetMaterial()
         
         Density(1:n) = GetReal( Material, 'Density', GotIt )
         IF( .NOT.GotIt ) THEN
           Density = 0.0d0
           IF( TransientSimulation .OR. (Solver % NOfEigenvalues > 0)) &
               CALL Fatal( 'ShellSolver', 'Density required' )
         END IF
         
         Poisson(1:n) = GetReal( Material, 'Poisson ratio', GotIt )
         IF( Isotropic .AND. (.NOT.GotIt) ) &
             CALL Fatal( 'ShellSolver', 'Poisson ratio undefined' )
         
         Young(1:n) = GetReal( Material, 'Youngs modulus', GotIt )
         IF( Isotropic .AND. (.NOT.GotIt) ) &
             CALL Fatal( 'ShellSolver', 'Youngs modulus undefined' )
         !	Thickness shifted to Body force section and reading of it shifted to LocalMatrix 14.10.15
         !       Thickness(1:n) = GetReal( Material, 'Thickness', GotIt )
         !       IF( Isotropic .AND. (.NOT.GotIt) ) &
         !                     CALL Fatal( 'ShellSolver', 'Thickness undefined' )
         
         Tension(1:n) = GetReal( Material, 'Tension', GotIt )
         IF( .NOT. GotIt ) Tension = 0.0d0
         
         CALL LocalMatrix(  STIFF, DAMP, MASS, &
             FORCE, LoadX, LoadY, LoadZ, LoadN, CurrentElement, n, &
             ElementNodes, StabParam1, StabParam2, t, Poisson,     &
             Young, LocalDeflection, LargeDeflection,  &
             StabilityAnalysis, Nvector,TT,ThLoad )
         
         IF( TransientSimulation ) THEN	
           CALL Default2ndOrderTime( MASS, DAMP, STIFF, FORCE )	
         END IF
!      Update global matrix and rhs vector from local matrix & vector:
!      ---------------------------------------------------------------
         CALL DefaultUpdateEquations( STIFF, FORCE )
         IF ( Solver % NOFEigenValues > 0 ) &
             CALL DefaultUpdateMass( MASS )
       END DO
! -----------------------------------------------------------------
    END SUBROUTINE BulkAssembly
! ----------------------------------------------------------------


! -------------------------------------------------------------
    SUBROUTINE BCAssembly()
! -----------------------------------------------------------------
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

        GotForceBC = .FALSE.
        LoadVector = 0.0d0

        LoadVector( 1, 1:n ) =  GetReal( BC, 'Force 1', GotIt )
        GotForceBC = GotForceBC.OR.GotIt

        LoadVector( 2, 1:n ) =  GetReal( BC, 'Force 2', GotIt )
        GotForceBC = GotForceBC.OR.GotIt

        LoadVector( 3, 1:n ) =  GetReal( BC, 'Force 3', GotIt )
        GotForceBC = GotForceBC.OR.GotIt

        LoadVector( 4, 1:n ) =  GetReal( BC, 'Force 4', GotIt )
        GotForceBC = GotForceBC.OR.GotIt

        LoadVector( 5, 1:n ) =  GetReal( BC, 'Force 5', GotIt )
        GotForceBC = GotForceBC.OR.GotIt

        LoadVector( 6, 1:n ) =  GetReal( BC, 'Force 6', GotIt )
        GotForceBC = GotForceBC.OR.GotIt

        IF( .NOT.GotForceBC ) CYCLE

        CALL StressBoundary( STIFF, FORCE, LoadVector, &
            CurrentElement, n, ElementNodes )

        CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO
! ----------------------------------------------------------------
    END SUBROUTINE BCAssembly
! ---------------------------------------------------------------


! --------------------------------------------------------------
    SUBROUTINE ConcentratedLoads()
! --------------------------------------------------------------
      INTEGER :: bf,i,nbfSet,i1,cont
 !     bf_id = ListGetInteger( Model % Bodies(1) % Values, 'Body Force' )
      bf = Model%NumberOfBodyForces
      nbfSet=0
      DO i1=1,bf
        PointLoad => ListGetConstRealArray( &
            Model % BodyForces( i1 ) % Values, 'Point Load', GotIt)
        IF( .NOT.GotIt ) THEN
          nPL = 0
        ELSE
          nPL = SIZE( PointLoad )/9
        END IF
        
        IF ( nPL > 0 ) THEN
          nbfSet=nbfSet+nPL	
	  
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
            
            cont=0
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
                cont=cont+1
              END IF
            END DO
            IF(cont == 0)THEN
              WRITE( Message,'(A)' ) 'ERROR-One of the point load is not located at node - rejected'	
              CALL Fatal('ShellSolve',Message)
              nbfSet=nbfSet-1	
            ELSE
            ENDIF
          END DO
        END IF
      ENDDO
      WRITE( Message,'(A,I9)' ) 'Number of point loads set',nbfSet	
      CALL Info('ShellSolve',Message)

! ------------------------------------------------------------------
    END SUBROUTINE ConcentratedLoads
! -------------------------------------------------------------

! --------------------------------------------------------
    SUBROUTINE SetDirichletBCs()
! --------------------------------------------------------------
      CALL DefaultDirichletBCs()
! ---------------------------------------------------------------
    END SUBROUTINE SetDirichletBCs
! -------------------------------------------------------------


! -----------------------------------------------------------
!
!                        ===================================
!                        S T R E S S   C O M P U T A T I O N
!                        ===================================
!
! -----------------------------------------------------------
    SUBROUTINE CalculateStresses()
! ------------------------------------------------------------
      INTEGER :: isz
      INTEGER, POINTER :: Perm(:)
      REAL(kind=dp) :: SHydro, SxD, SyD, SzD, SumSTaoSq
      at  = CPUTime()
      at0 = RealTime()

!      Allocate memory for the local stresses:
!      ---------------------------------------
      ALLOCATE( Weights( Solver % Mesh % NumberOfBulkElements, NumberOfElementNodes ) )
      Weights = 0.0d0
      
      ALLOCATE( SxxElement( Solver % Mesh % NumberOfBulkElements ) )
      ALLOCATE( SyyElement( Solver % Mesh % NumberOfBulkElements ) )
      ALLOCATE( SzzElement( Solver % Mesh % NumberOfBulkElements ) )
      ALLOCATE( SxyElement( Solver % Mesh % NumberOfBulkElements ) )
      ALLOCATE( SxzElement( Solver % Mesh % NumberOfBulkElements ) )
      ALLOCATE( SyzElement( Solver % Mesh % NumberOfBulkElements ) )
      
      SxxElement = 0.0d0
      SyyElement = 0.0d0
      SzzElement = 0.0d0
      SxyElement = 0.0d0
      SxzElement = 0.0d0
      SyzElement = 0.0d0
      
      ALLOCATE( EpsxxElement( Solver % Mesh % NumberOfBulkElements ) )
      ALLOCATE( EpsyyElement( Solver % Mesh % NumberOfBulkElements ) )
      ALLOCATE( EpszzElement( Solver % Mesh % NumberOfBulkElements ) )
      ALLOCATE( EpsxyElement( Solver % Mesh % NumberOfBulkElements ) )
      ALLOCATE( EpsxzElement( Solver % Mesh % NumberOfBulkElements ) )
      ALLOCATE( EpsyzElement( Solver % Mesh % NumberOfBulkElements ) )
      
      EpsxxElement = 0.0d0
      EpsyyElement = 0.0d0
      EpszzElement = 0.0d0
      EpsxyElement = 0.0d0
      EpsxzElement = 0.0d0
      EpsyzElement = 0.0d0
      

!      Then, compute the element stresses:
!      ----------------------------------- 

      OPEN( UNIT=98, file='nvec.txt' )
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

	!	Thickness shifted to Body force section and reading of it shifted to 
	!	calculatestresses and localMatrix 15.10.15
        Thickness(1:n) = GetConstReal(BodyForce, 'Thickness', GotIt )
        IF( Isotropic .AND. (.NOT.GotIt) ) &
            CALL Fatal( 'ShellSolver', 'Thickness undefined' )

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
            Poisson, Thickness, LargeDeflection,t )
        
        IF( ListGetLogical( SolverParams,'Compute Membrane Stress',GotIt) ) THEN
          SxxElement(t) = Nten(1,1)
          !changed Compute strain to Compute Membrane Stress20.3.15
          SyyElement(t) = Nten(2,2)
          SzzElement(t) = Nten(3,3)
          SxyElement(t) = Nten(1,2)
          SxzElement(t) = Nten(1,3)
          SyzElement(t) = Nten(2,3)
          
          EpsxxElement(t) = Eps(1,1)
          EpsyyElement(t) = Eps(2,2)
          EpszzElement(t) = Eps(3,3)
          EpsxyElement(t) = Eps(1,2)*2.0
          EpsxzElement(t) = Eps(1,3)*2.0
          EpsyzElement(t) = Eps(2,3)*2.0
!	Shear strains multiplied by 2 because these are the components of strain tensor which are half of strains

        ELSE IF( ListGetLogical( SolverParams,'Compute Bending Stress',GotIt) ) THEN  
          SxxElement(t) = Mten(1,1)	!Changed Compute Curvature to Compute Bending Stress 20.3.15 
          SyyElement(t) = Mten(2,2)
          SzzElement(t) = Mten(3,3)
          SxyElement(t) = Mten(1,2)
          SxzElement(t) = Mten(1,3)
          SyzElement(t) = Mten(2,3)
          
          EpsxxElement(t) = Kap(1,1)
          EpsyyElement(t) = Kap(2,2)
          EpszzElement(t) = Kap(3,3)
          EpsxyElement(t) = Kap(1,2)*2.0
          EpsxzElement(t) = Kap(1,3)*2.0
          EpsyzElement(t) = Kap(2,3)*2.0

        ELSEIF( ListGetLogical( SolverParams, 'Compute Total Stress', GotIt ) ) THEN
          SxxElement(t) = Mten(1,1)+Nten(1,1)
          SyyElement(t) = Mten(2,2)+Nten(2,2)
          SzzElement(t) = Mten(3,3)+Nten(3,3)
          SxyElement(t) = Mten(1,2)+Nten(1,2)
          SxzElement(t) = Mten(1,3)+Nten(1,3)
          SyzElement(t) = Mten(2,3)+Nten(2,3)
          
          EpsxxElement(t) = Kap(1,1)+Eps(1,1)
          EpsyyElement(t) = Kap(2,2)+Eps(2,2)
          EpszzElement(t) = Kap(3,3)+Eps(3,3)
          EpsxyElement(t) = (Kap(1,2)+Eps(1,2))*2.0
          EpsxzElement(t) = (Kap(1,3)+Eps(1,3))*2.0
          EpsyzElement(t) = (Kap(2,3)+Eps(2,3))*2.0
        END IF

        SELECT CASE( NumberOfElementNodes )
        CASE( 3 )
          Weights(t,1:3) = Weight3(1:3)
        CASE( 4 )
          Weights(t,1:4) = Weight4(1:4)
        END SELECT
        
      END DO		! Loop on elements ends
      CLOSE( 98 )

      isz = Solver % Mesh % NumberOfNodes
      
      IF( .NOT.ASSOCIATED( VariableGet( Solver % Mesh % Variables, &
          'Stress.xx') ) ) THEN
        
        ALLOCATE( SxxNodal( isz ) )
        ALLOCATE( SyyNodal( isz ) )
        ALLOCATE( SzzNodal( isz ) )
        ALLOCATE( SxyNodal( isz ) )
        ALLOCATE( SxzNodal( isz ) )
        ALLOCATE( SyzNodal( isz ) )
        ALLOCATE( SEqNodal( isz ) )
        
        ALLOCATE( EpsxxNodal( isz ) )
        ALLOCATE( EpsyyNodal( isz ) )
        ALLOCATE( EpszzNodal( isz ) )
        ALLOCATE( EpsxyNodal( isz ) )
        ALLOCATE( EpsxzNodal( isz ) )
        ALLOCATE( EpsyzNodal( isz ) )
        
        ALLOCATE(Perm(isz))
        
        Perm = (/ (i,i=1,isz) /)
        
        PSolver => Solver
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Stress.xx', 1, SxxNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Stress.yy', 1, SyyNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Stress.zz', 1, SzzNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Stress.xy', 1, SxyNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Stress.xz', 1, SxzNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Stress.yz', 1, SyzNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Stress.Eq', 1, SEqNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Epsilon.tot.xx', 1, EpsxxNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Epsilon.tot.yy', 1, EpsyyNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Epsilon.tot.zz', 1, EpszzNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Epsilon.tot.xy', 1, EpsxyNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Epsilon.tot.xz', 1, EpsxzNodal,Perm )
        
        CALL VariableAdd( Solver % Mesh % Variables, &
            Solver % Mesh, PSolver, 'Epsilon.tot.yz', 1, EpsyzNodal,Perm )
      END IF
      
      
      !       Average nodal values:
!       ---------------------
      SxxNodal = 0.0d0
      SyyNodal = 0.0d0
      SzzNodal = 0.0d0
      SxyNodal = 0.0d0
      SxzNodal = 0.0d0
      SyzNodal = 0.0d0
      SEqNodal = 0.0d0
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
          k = NodeIndexes(j)
          IF( k>0 ) THEN
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
          END IF
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

      ! Calculate Equivalent stress 
      DO i=1, Solver % Mesh % NumberOfNodes
        SHydro = (SxxNodal(i)+SyyNodal(i)+SzzNodal(i))/3.0
        SxD = (SxxNodal(i)-SHydro)
        SyD = (SyyNodal(i)-SHydro)
        SzD = (SzzNodal(i)-SHydro)
        SumSTaoSq=SxyNodal(i)*SxyNodal(i)+SxzNodal(i)*SxzNodal(i)+SyzNodal(i)*SyzNodal(i)
        SEqNodal(i)=SQRT(1.5*(SxD*SxD+SyD*SyD+SzD*SzD+2.0*SumSTaoSq))
      ENDDO

      at = CPUTime() - at
      
      WRITE(Message,'(a,F8.2)') ' Stress calculation: (s)', at
      CALL Info('ShellSolver',Message)

!       Finally, release the auxiliary arrays:
!       -------------------------------------- 
      DEALLOCATE( SxxElement, SyyElement, SzzElement, SxyElement, &
          SxzElement, SyzElement, Weights, Referenced, &
          EpsxxElement, EpsyyElement, EpszzElement, EpsxyElement, &
          EpsxzElement, EpsyzElement )
! ----------------------------------------------------------------
    END SUBROUTINE CalculateStresses
! -----------------------------------------------------------


! ---------------------------------------------------------------
    SUBROUTINE StressBoundary( STIFF, FORCE, LOAD, Element, n, Nodes )
! ------------------------------------------------------------
      REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:,:)
      TYPE(Element_t), POINTER :: Element
      INTEGER :: n
      TYPE(Nodes_t) :: Nodes, ElementNodes
! ------------------------------------------------------------------
      REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
          detJ, s, u, v, w, LoadAtIp(6)
      INTEGER :: t, i, j
      LOGICAL :: stat
      TYPE( GaussIntegrationPoints_t ), TARGET :: IntegStuff
      
      STIFF = 0.0d0
      FORCE = 0.0d0
      IntegStuff = GaussPoints( element )
      CALL GetElementNodes( ElementNodes )		
      
      DO t = 1, IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)
        
        stat = ElementInfo( Element, ElementNodes, u , v, w, &   ! 
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
! --------------------------------------------------------------
    END SUBROUTINE StressBoundary
!  -------------------------------------------------------------

! ----------------------------------------------------------
    SUBROUTINE LocalMatrix( STIFF, DAMP, MASS, &
        FORCE, NodalLoadX, NodalLoadY, NodalLoadZ, NodalLoadN, &
        Element, N, Nodes, StabParam1, StabParam2, &
        ElementNumber, NodalPoisson, NodalYoung, &
        LocalDeflection, LargeDeflection, StabilityAnalysis, Nvector, TT,ThLoad ) 
!	If nodal thickness is varying, thickness(1:n) should als be input through LocalMatrix()
!	NOTE:- The 'Element' here is the 'Current Element'. See the Call statement in sub BulkAssembly
!	Thickness is being read here now.17.10.15
! ----------------------------------------------------------------
      REAL(KIND=dp) :: STIFF(:,:), DAMP(:,:), MASS(:,:), &
          Amatrix(3,3), Bmatrix(3,3), Dmatrix(3,3), Astarmatrix(2,2)
      REAL(KIND=dp) :: FORCE(:)
      REAL(KIND=dp) :: NodalLoadX(:), NodalLoadY(:), NodalLoadZ(:), &
          NodalLoadN(:), LocalDeflection(:), Nvector(:), &
          ThLoad(:),TT(:),BendingLoad1(6*n),ThLoad1(6*n)				
      REAL(KIND=dp) :: StabParam1, StabParam2
      REAL(KIND=dp) :: NodalPoisson(:), NodalYoung(:)
      LOGICAL :: LargeDeflection, StabilityAnalysis
      INTEGER :: N, ElementNumber,kk(n), IntPt	
      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t), POINTER :: Element
! ---------------------------------------------------------------
      INTEGER, PARAMETER :: MaxNodes = 4
      INTEGER, PARAMETER :: MaxDofs = 6*MaxNodes
      
      REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
          Kappa(3,MaxDofs), Gammaa(2,MaxDofs), EPS(3,MaxDofs), &
          Omega(1,MaxDofs), Gdrilling(1,1), GradDeflection(2,MaxDofs), & !*
          ZetaStrain(3,maxDofs)
      REAL(kind=dp) :: Euvw,Puvw,Tuvw,Coef1,CStrFactor,RotMatrix2D(3,3)
      REAL(KIND=dp) :: detJ, U, V, W, S, SCF, rho, h, &
          LoadN, LoadX, LoadY, LoadZ, Nmatrix(2,2), &
          GradTest(2), GradBasis(2), LV1(3,1), LV2(3,1), TempVec(3,1), &
          NormalForce(2,2), NonLinForce(MaxDofs), LV3(3,1)
      
      REAL(KIND=dp) :: Transformation(3,3), CopyOfNodes(3,MaxNodes), &
          T0(3,3), T0S(2,2), T5(5,5)
      
      REAL( KIND=dp ) :: Tblock(6*n,6*n)
      
      LOGICAL :: Stat, Found,Found1	
      INTEGER :: i,p,q,t, pk, qk
      TYPE( GaussIntegrationPoints_t ) :: IntegStuff
      
      REAL(KIND=dp) :: LamDCMatrix(8,8), dWdx, dWdy
      
      REAL(KIND=dp) :: Moment(2,2), ZetaStress(2,2), EpsStress(2,2), &
          dUdx(3,2), dRdx(3,2)
      CHARACTER(LEN=15) :: Aniso, Aniso1
      REAL(Kind=dp) :: Point1(3), Point2(3),Tmatrix(3,3),Te,alfa(2),sinTh,cosTh	
      REAL(Kind=dp) :: aDKT(6),bDKT(6),cDKT(6),dDKT(6),eDKT(6),pDKT(6),qDKT(6),tDKT(6),rDKT(6) 
      REAL(Kind=dp) :: Ledge(6),x(3),y(3),z(3), Area2  
      REAL(Kind=dp) :: xx(3),yy(3),zz(3),y13,dwx,dwy,dwx1,dwy1,sNew
      INTEGER :: NodeNum(3)
      REAL(kind=dp) :: GSTIFF(6*n,6*n),CStrain(3)
      REAL(Kind=dp) :: MemMatrix(6*n,6*n)
      REAL(kind=dp) :: beta(1,1),Gamma1(1,6*n)
      REAL(kind=dp) :: IniBen2(2)
! --------------------------------------------------------------

      FORCE = 0.0d0
      STIFF = 0.0d0
      DAMP  = 0.0d0
      MASS  = 0.0d0
      ThLoad1=0.0d0	
      ThLoad = 0.0d0	
      BendingLoad=0.0d0	
      BendingLoad1=0.0d0	
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

      ! *****************************Copied from LocalStress ****************************
      ! Added here to incorporate value of OUtNormal for use in IniBen2( ) for correcting its direction
      ! BodyForce => GetBodyForce() Removed because this is already done in calling subroutine
      OutNormal = GetConstReal(BodyForce,'Direction Of Outward Normal',Found1)
      IF(.NOT.Found1)THEN
        OutNormal=1.0
      ENDIF

      !	*****************Following shifted from Sub 'ShellSolver' to this place ***********
      !	Location of reading these quantities shifted to 'BodyForce' section.
      !	The objective is to calculate weld simulation and thermal strain related parametese only 
      !	for the body where these are to be calculated.
      !	Weld Distortion Simulation and Thermal Strain
      WeldSimu = GetLogical( BodyForce, 'Weld Distortion Simulation', GotIt )
      IF( .NOT. GotIt ) WeldSimu = .FALSE.

      ThStrain = GetLogical( BodyForce, 'Calculate Thermal Strain', GotIt )
      IF( .NOT. GotIt ) ThStrain = .FALSE.
      IF(WeldSimu .AND. ThStrain )THEN
        CALL warn( 'ShellSolver', 'Weld Distortion Simulation and Calculate Thermal Strain'// &
            ' cannot be specified simultaneously' )
        CALL fatal( 'ShellSolver', 'STOP')
      ENDIF
      AnisoPlane = GetLogical( BodyForce, 'Use Relaxed AnisoPlane Testing', GotIt )
      IF( .NOT. GotIt ) AnisoPlane = .FALSE.
!	******Following shifted from BulkAssembly subroutine to this place (to be more logical)******
!	   (Should be modified for actual temperature distribution)
      TT(1:n) = GetReal( BodyForce, 'Temperature', GotIt )
      IF((.NOT.GotIt).AND.( ThStrain .OR. WeldSimu )) THEN	
        WRITE( Message,'(A,I9)' ) &
            'Temperature not specified for thermal strain or weld simulation calculations Taking it zero'
        CALL Info('ShellSolve',Message)
        TT(1:n) = 0.0d0
      ENDIF

      ThermalAnisotropy=.FALSE.
	!	If Anisotropy is specified then only parameters related to it are read and axes found
      Aniso1 =GetString( BodyForce, 'Thermal Anisotropy', Found)
      IF(Found)THEN
        ThermalAnisotropy=.TRUE.
        !			WRITE( Message,'(A,L2)' ) 'ThermalAnisotropy ',thermalanisotropy
        !			CALL warn('ShellSolve',Message)	
        
        IF((ThStrain).OR.(WeldSimu))THEN
          CALL GetAnisoParameter(Aniso,Point1,Point2,alfa,Te,IniBen2)
          CALL GetAnisotropyAxes(Element,Aniso,Point1, Point2,cosTh,sinTh)
          Tmatrix= ThStrainConversion(cosTh,sinTh)
        ENDIF
	!	Te is the Weld Simulation Temperature, read in sub GetAnisoParameter()
        IF(WeldSimu)TT(1:n)=Te					
      ELSEIF(ThStrain)THEN
        !	Get alfa1 for calculations of non-anisotropic thermal strain
        alfa(1) = GetConstReal(BodyForce,'Alfa 1',Found1)
        IF(.NOT.Found1)THEN
          WRITE( Message,'(A)' ) 'Alfa 1 not specified for thermal strain calculations '
          CALL warn('ShellSolve',Message)
          CALL Fatal('ShellSolve','STOP')
        ENDIF
      ENDIF					

      !	Thickness shifted to Body force section and reading of it shifted to LocalMatrix 14.10.15
      Thickness(1:n) = GetConstReal(BodyForce, 'Thickness', GotIt )
      IF( Isotropic .AND. (.NOT.GotIt) ) THEN
        CALL Fatal( 'ShellSolver', 'Thickness undefined' )
      END IF

      IF ( StabilityAnalysis ) THEN
        CALL LocalStress( Element, n, Nodes, StabParam1,    &
            StabParam2, LocalDeflection, Weight3, Weight4,    &
            Eps, Kap, Nten, NtenMaterial, Mten, MtenMaterial, &
            NodalYoung, NodalPoisson, Thickness, LargeDeflection,ElementNumber )
      END IF
      
! The transformation Xglob -> Xloc is the transpose of the local basis:
!      ----------------------------------------------------------
      Transformation = TRANSPOSE( LocalBasis( Nodes, n ) )
       
! Take a copy of the global node points and switch to the local system:
!      ----------------------------------------------------------
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
      LocalDeflection = MATMUL( Tblock, LocalDeflection )

!      Exchange rotations r1 and -r2:
!      ------------------------------
!	(28.3.15) The following is inconsequential here because the value of LocalDeflection, as
!	input to LocalMatrix is always zero uinless the subroutine LocalMatrix is called repeatedly,
!	may be in the LargeDisplacement analysis.  
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
!      Select the appropriate quadrature:
      IF(n /= 3 .AND. n /= 4 )THEN
        WRITE( Message,'(A)' ) &
            'Only 3 node triangle or 4 node quad are currently acceptable'
        CALL warn('ShellSolve', Message)
        CALL Fatal('ShellSolve','STOP')        
      ENDIF

115   CONTINUE

      
      SELECT CASE( Element % TYPE % NumberOfNodes )


      CASE( 3 )
        IF((UseDKTtriangle).OR.(UseRDKTtriangle).OR.(UseMITC3element).OR.&
            (UseRMITC3element))THEN	! Added MITC3, DMITC3, RMITC3 options 
          IntPt = GetInteger( SolverParams, 'Integration Points', GotIt )
          IF( .NOT. GotIt ) THEN
            WRITE( Message,'(A)' ) &
                'Number of Integration Points not specified for DKT/RDKT/MITC3/RMITC3 triangle '
            CALL Fatal('ShellSolve',Message)
          ELSE
            IF((IntPt == 3).OR.(IntPt == 4).OR.(IntPt == 7).OR.(IntPt == 1))THEN
              IntegStuff = GaussPoints(Element, IntPt)
            ELSE
              WRITE( Message,'(A)' ) &
                  'Number of Integration Points for DKT/RDKT/MITC3/RMITC3  triangle should be 1,3,4 or 7 only'
              CALL warn('ShellSolve',Message)
              CALL Fatal('ShellSolve','STOP')
            ENDIF
          ENDIF
        ELSE          
          IntegStuff = GaussPoints( Element, 1 )
          IntPt=1	
        END IF

        IF(ElementNumber == 1)THEN					
          CALL Info('ShellSolve','Number of Integration Points = '//TRIM(I2S(IntPt)),Level=10)
        ENDIF														

      CASE( 4 )
        IF((UseDKTtriangle).OR.(UseRDKTtriangle).OR.&
            (MembraneOnly).OR.(UseMITC3element).OR.&
            (UseRMITC3element).OR.(MembraneOnlyNRM))THEN	
          WRITE( Message,'(A)' ) &
              'Use of chosen element type is available only with three node triangle'
          CALL Fatal('ShellSolve',Message)
        ELSE
          IntegStuff = GaussPoints( Element, 4 )	
        ENDIF
      END SELECT

      !	calculate basic constants of the DKT element
      IF (UseDKTtriangle .OR. UseRDKTtriangle)THEN        
        CALL getDKTtriangleParameter(Nodes,n,aDKT,bDKT,cDKT,dDKT,eDKT,pDKT,qDKT,tDKT,rDKT)
      ENDIF

      DO i=1,3
        xx(i)=Nodes%x(i)
        yy(i)=Nodes%y(i)
        zz(i)=Nodes%z(i)
      ENDDO

      Area2 = (xx(2)-xx(1))*yy(3)+(xx(1)-xx(3))*yy(2)+(xx(3)-xx(2))*yy(1) 
      y13 = yy(3)-yy(1)
      
      !      Numerical integration:
      !      ----------------------
      DO t = 1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)
        
!        Basis function values & derivatives at the integration point:
!        -------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, U, V, W, &	!Nodes changed to ElementNodes31.12.13 -Reverted
            detJ, Basis, dBasisdx )
        IF( ElementNumber == 15 ) THEN
          
          WRITE(Message,'(a,I6,8f12.6)') ' Int Pt No, u, v, w, s, detJ, Area2,dbasisdx(2,1),y13', &
              t,U,V,W,S,detJ,Area2,dbasisdx(2,1),y13
          CALL Info('ShellSolver',Message)
        ENDIF
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
          !             Nonlinear terms:
          !            -----------------
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
          
 !            Normal force:
 !           --------------
          TempVec = MATMUL( Amatrix, LV1 ) + MATMUL( Bmatrix, LV2 )
          
          NormalForce(1,1) = TempVec(1,1)
          NormalForce(2,2) = TempVec(2,1)
          NormalForce(1,2) = TempVec(3,1)
          NormalForce(2,1) = TempVec(3,1)
          
!             Bending moment:
!           ----------------
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

        ! --------------------------------------------------------------------
        !	***********************Inserted 6.11.14-MembraneOnly**************************
        IF(MembraneOnly .OR. MembraneOnlyNRM ) THEN
          CALL MembraneFormulation(MemMatrix,n)
          GOTO 1001
	ENDIF
 ! Bending stiffness:
 !        ------------------
      Kappa = 0.0d0

      IF ((UseDKTtriangle))THEN
        CALL DKTformulation(Nodes, Kappa,MaxDofs,n, u,v,w,s,aDKT,bDKT,cDKT,dDKT,eDKT,pDKT,qDKT,tDKT,rDKT)
        ! u,v,w,s are the zy,eta,zeta coordinates and area weightage of int pt.
        STIFF(1:6*n,1:6*n)=STIFF(1:6*n,1:6*n)+s*MATMUL(MATMUL(TRANSPOSE(Kappa),Dmatrix),Kappa)
      ELSEIF((UseRDKTtriangle))THEN
        CALL RDKTelement(Nodes, Kappa,MaxDofs,n, u,v,w,s,aDKT,dDKT)
        STIFF(1:6*n,1:6*n)=STIFF(1:6*n,1:6*n)+s*MATMUL(MATMUL(TRANSPOSE(Kappa),Dmatrix),Kappa)
      ELSE
        !	Calculate the bending strain and Stiffness matrix  in the usual manner (Not DKT)

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
        STIFF(1:6*n,1:6*n)=STIFF(1:6*n,1:6*n)+s*MATMUL(MATMUL(TRANSPOSE(Kappa),Dmatrix),Kappa)

!         CALL AddEnergy(STIFF, Dmatrix, Kappa, 3, 6*n, s)
      END IF	!inserted on 21.3.14 - Subroutine DKTformulation( ) writes complete Kappa for the Int. Pt. 
!	*********************Inserted on 6 11 14 MembraneOnly**************************
1001  IF( MembraneOnly .OR. MembraneOnlyNRM )THEN		
        STIFF(1:6*n,1:6*n)=STIFF(1:6*n,1:6*n)+s*MemMatrix
      ENDIF
	

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
      
!	Enforcing minimization of transverse rotation
      IF(MembraneOnly)THEN 
        Euvw = SUM( Young(1:n)    * Basis(1:n) )
        Puvw = SUM( Poisson(1:n)  * Basis(1:n) )
        Tuvw = SUM( Thickness(1:n)* Basis(1:n) )
        Coef1 = Euvw*Tuvw  /( 2.0d0 * (1.0d0 +  Puvw))
        !		go to 1010		!Very Temp 18.12.14 (No in-plane rotation minimization)

!	Rotation about y-axis
        Gamma1=0.0d0
        beta(1,1)= Coef1
        DO p=1,n
          Gamma1(1,6*p-3) = 0.5*dBasisdx(p,1)
        ENDDO
        CALL AddEnergy(STIFF, beta, Gamma1, 1, 6*n, s)
        

!	Rotation about x-axis
        Gamma1=0.0d0
        beta(1,1)= Coef1
        DO p=1,n
          Gamma1(1,6*p-3) = 0.5*dBasisdx(p,2)
        ENDDO
        CALL AddEnergy(STIFF, beta, Gamma1, 1, 6*n, s)
!	********************************Insertion 7.11.14 ends*********************
! 1010	continue		!Very Temp 18.12.14
!	Added21.11.14 - Consideration of rotation about z-axis
        Gamma1=0.0d0
        beta(1,1)= Coef1
        DO p=1,n
          Gamma1(1,6*p-5) = 0.5*dBasisdx(p,2)
          Gamma1(1,6*p-4) = -0.5*dBasisdx(p,1)
        ENDDO
        CALL AddEnergy(STIFF, beta, Gamma1, 1, 6*n, s)
        
      ENDIF

      
      IF( MembraneOnlyNRM )THEN 
        Euvw = SUM( Young(1:n)    * Basis(1:n) )
        Puvw = SUM( Poisson(1:n)  * Basis(1:n) )
        Tuvw = SUM( Thickness(1:n)* Basis(1:n) )
        Coef1 = Euvw*Tuvw  /( 2.0d0 * (1.0d0 +  Puvw))
        Gamma1=0.0d0
        beta(1,1)= Coef1
        DO p=1,n
          Gamma1(1,6*p-5) = 0.5*dBasisdx(p,2)
          Gamma1(1,6*p-4) = -0.5*dBasisdx(p,1)
        ENDDO
        CALL AddEnergy(STIFF, beta, Gamma1, 1, 6*n, s) 
      ENDIF
      
!        Shear stiffness (transversal):
!        ------------------------------
!	Skipping shear term for membrane Only	
      IF( MembraneOnly .OR. MembraneOnlyNRM ) GOTO 1002 

!	-----------------------------------------------------------------------------------------
      IF( UseMITC3element )THEN
        Gammaa=0.0			
        
        CALL CovariantInterpolationMITC3(Gammaa, Basis, &
            Nodes % x(1:n), Nodes % y(1:n),U, V, n)
        IF( ElementNumber == 10 ) THEN
          WRITE( Message,'(A)' ) &
              'Program is now using MODIFIED subroutine - CovariantInterpolationMITC3 and MITC3 element' 
          CALL Info('ShellSolve',Message)
        ENDIF
        
        CALL ShearCorrectionFactor(SCF, h, Nodes % x(1:n), &
            Nodes % y(1:n), n, StabParam1)	!SCF is 1 for stabParam1=0 and little smaller for + value of it (See Subroutine 22.8.14)
        !		SCF=1.0 !Added on 9.3.15. Also above subroutine for SCF removed !Reverted back 15.10.15
        CALL AddEnergy(STIFF, Astarmatrix, Gammaa, 2, 6*n, SCF*s)

      ELSE IF( UseRMITC3element )THEN
        Gammaa=0.0			
        
        CALL CovariantInterpolationRMITC3(Gammaa, Basis, &
            Nodes % x(1:n), Nodes % y(1:n),U, V, n)
        IF( ElementNumber == 10 ) THEN
          WRITE( Message,'(A)' ) &
              'Program is now using MODIFIED subroutine - CovariantInterpolationRMITC3 and RMITC3 element' 
          CALL Info('ShellSolve',Message)
        ENDIF
        CALL ShearCorrectionFactor(SCF, h, Nodes % x(1:n), &
            Nodes % y(1:n), n, StabParam1)	!SCF is 1 for stabParam1=0 and little smaller for + value of it (See Subroutine 22.8.14)
        
        ! 		SCF=1.0 !Added on 9.3.15. Also above subroutine for SCF removed !Reverted back 15.10.15
        
        DO p=1,n				
          Gammaa(1:2,6*p-3) =  dBasisdx(p,1:2)!This implies coeff of nodal w as DoverN/dx and DoverN/dy 
          !		Above value( dBasisdx( )) whether + or - does not matter in cantilever beam problem (originally +).
          !		 It is just to incorporate some terms in w so that singularity in matrix is avoided, 
        END DO		!These are incorporated in respective columns (3rd, 9th and 15th columns). Note that

        CALL AddEnergy(STIFF, Astarmatrix, Gammaa, 2, 6*n, SCF*s)
!	---------------------------------------------------------------------------------------
      ELSE IF( UseSMITCelement)THEN
        gammaa=0.0d0	
        CALL CovariantInterpolation(Gammaa, Basis, &	!As it existed in original(8.9.14)
            Nodes % x(1:n), Nodes % y(1:n), U, V, n)		!See slight modification there15.10.15
        
        IF( ElementNumber == 10 ) THEN
          WRITE( Message,'(A)' ) &
              'Program is now using SMITC element as it exists' 
          CALL Info('ShellSolve',Message)
        ENDIF
        
        CALL ShearCorrectionFactor(SCF, h, Nodes % x(1:n), &
            Nodes % y(1:n), n, StabParam1)	!SCF is 1 for stabParam1=0 and little smaller for + value of it (See Subroutine 22.8.14)
        DO p=1,n
          Gammaa(1:2,6*p-3) =  dBasisdx(p,1:2)
        END DO
        !		Gammaa is dimensioned as (2,6*n) 22.8.14
        
        CALL AddEnergy(STIFF, Astarmatrix, Gammaa, 2, 6*n, SCF*s) 
        
      ENDIF
!	****************Temp including drilling DOF 20.11.14**********************
 !1002	continue !Inclusion of minimization of drilling DOF was inconsequential
!	*********************Temp iclusion ends 20.11.14***************************
!        Drilling DOFs (in-plane rotations):
!        -----------------------------------
      Omega = 0.0d0
      DO p = 1,n
        Omega(1,6*p-5) = +dBasisdx(p,2) / 2.0d0  !  u_{x,y}
        Omega(1,6*p-4) = -dBasisdx(p,1) / 2.0d0  ! -u_{y,x}
        Omega(1,6*p-0) =  Basis(p)              !  rotation
      END DO
      
      CALL AddEnergy(STIFF, Gdrilling, Omega, 1, 6*n, s)
1002  CONTINUE	
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
!	***************Explanation 2.8.15***************************************
!	LoadX, LoadY and LoadZ are Body Force 1,2,3. These may be gravity force etc. and 
!	so these are taken per unit volume basis. So the actual load is obtained by 
!	multiplying load per unit volume by volume at integration point (i.e. the product 
!	of h and s). So if surface load in X, Y, Z direction is to be specified, these should be
!	converted to per unit volume basis.

	!	Add thermal load vector	at Int Pt and then transform to Global
	!	Note that Amatrix and EPS are already calculated.


      IF( WeldSimu .OR. ThStrain )THEN
        CALL ThermalAndBendingLoad(Amatrix,Dmatrix,EPS,Kappa,Tmatrix,Basis,dBasisdx,s,alfa,TT,&
            ThLoad1,n,MaxDofs,IniBen2,BendingLoad1,RotMatrix2D)	

        BendingLoad(1:6*n)=BendingLoad(1:6*n)+BendingLoad1(1:6*n)
        ThLoad=ThLoad+ThLoad1
	!ThLoad is initialized to be zero in the beginning of subroutine 'LocalMatrix'. Thus 
	!when it is added to the previous value in this subroutine, the values from all int pts are added. 

!	   WRITE( Message,'(A,I4,18E12.4)' ) 'IntPt and ThLoad 1 ',IntegStuff%n,(ThLoad(i),i=1,18)
        !      CALL Info('ShellSolve',Message)
        
      ENDIF

	!	(IMPORTANT This is within Int Pt loop)
	!	Transform to Global coordinates
	!	Next add ThLoad to FORCE

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
      END IF				!newton lin terms - large deflection ends 

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
      
    END DO ! End of loop at integration points  
    
!      Restore the original node points:
!      ---------------------------------
    Nodes % x(1:n) = CopyOfNodes(1,1:n)
    Nodes % y(1:n) = CopyOfNodes(2,1:n)
    Nodes % z(1:n) = CopyOfNodes(3,1:n)

    ! Finally, we perform some transformations:
    !-----------------------------------------
    Tblock=0.0d0
    DO i = 1,n
      Tblock( 6*i-5, 6*i-5 ) =  1.0d0
      Tblock( 6*i-4, 6*i-4 ) =  1.0d0
      Tblock( 6*i-3, 6*i-3 ) =  1.0d0
      Tblock( 6*i-2, 6*i-1 ) =  -1.0d0
      Tblock( 6*i-1, 6*i-2 ) = 1.0d0
      Tblock( 6*i-0, 6*i-0 ) =  1.0d0
    END DO

    IF((.NOT.UseDKTtriangle) .AND. (.NOT.UseRDKTtriangle) )THEN	
!	i.e. no EXCHANGE of betax, betay to Thetax and thetay (It is already done in formulation)
      i = 6*n
      STIFF(1:i,1:i) = MATMUL( STIFF(1:i,1:i), Tblock(1:i,1:i) )
      STIFF(1:i,1:i) = MATMUL( TRANSPOSE( Tblock(1:i,1:i) ), STIFF(1:i,1:i) )
      !	The BendingLoad is is not to be multiplied by Tblock for DKT element because 
!	the same has been considered during formulation of matrix [k].
      BendingLoad(1:i)=MATMUL ( TRANSPOSE(Tblock(1:i,1:i)),BendingLoad(1:i) )	
    ENDIF
								!End of EXCHANGE	
!	***********Stability analysis and large deflection not tested for ************
!	***********definition of BetaX as -du/dx and BetaY=-dv/dx*********************
!	***********here as well as at other places15.10.15 ***************************
    IF( StabilityAnalysis ) THEN
      MASS(1:i,1:i) = MATMUL( MASS(1:i,1:i), Tblock(1:i,1:i) )
      MASS(1:i,1:i) = MATMUL( TRANSPOSE( Tblock(1:i,1:i) ),  MASS(1:i,1:i) )
    END IF
    
    IF( LargeDeflection ) THEN
      NonlinForce(1:i) = MATMUL( TRANSPOSE(Tblock(1:i,1:i)), NonlinForce(1:i) )
    END IF
    !	******************************************************************************
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
      NonlinForce(1:i) = MATMUL( TRANSPOSE(Tblock(1:i,1:i)),  NonlinForce(1:i) )
      FORCE = FORCE + NonlinForce
    END IF
    
    IF( StabilityAnalysis ) THEN
      MASS(1:i,1:i) = MATMUL( MASS(1:i,1:i), Tblock(1:i,1:i) )
      MASS(1:i,1:i) = MATMUL( TRANSPOSE( Tblock(1:i,1:i) ),  MASS(1:i,1:i) )
    END IF
    
       STIFF = ( STIFF + TRANSPOSE(STIFF) ) / 2.0d0
       MASS  = ( MASS  + TRANSPOSE(MASS) )  / 2.0d0
       !	Add thermal load after rotational transformation

       ThLoad(1:i) = MATMUL ( TRANSPOSE(Tblock(1:i,1:i)),ThLoad(1:i) )	
       BendingLoad(1:i)=MATMUL ( TRANSPOSE(Tblock(1:i,1:i)),BendingLoad(1:i) )
       FORCE = FORCE + ThLoad+BendingLoad		
       
     END SUBROUTINE LocalMatrix
     ! ----------------------------------------------------------------

!  ================================================================

! -------------------------------------------------------------- 	
     SUBROUTINE LocalStress( Element, n, Nodes, StabParam1,  StabParam2,LocalDeflection, Weight3, Weight4, Eps, Kap, Nten,&
         NtenMaterial, Mten, MtenMaterial, NodalYoung, NodalPoisson, NodalThickness, LargeDeflection, ElementNumber )
! -----------------------------------------------------------------
       IMPLICIT NONE
       REAL(KIND=dp) :: StabParam1, StabParam2, LocalDeflection(:), &
           Weight3(:), Weight4(:), Eps(3,3), Kap(3,3), NTen(3,3),MTen(3,3),NtenMaterial(2,2), MtenMaterial(2,2),  &
           NodalYoung(:) , NodalPoisson(:), NodalThickness(:)
       INTEGER :: n, ElementNumber
       TYPE( GaussIntegrationPoints_t ) :: IntegStuff	
       
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
       INTEGER, PARAMETER :: MaxNodes = 4, MaxDofs = 6*MaxNodes
       LOGICAL :: LargeDeflection
       CHARACTER(LEN=15) :: Aniso, Aniso1					
       REAL(Kind=dp) :: Point1(3), Point2(3),Tmatrix(3,3),Te,alfa(2),sinTh,cosTh	
       REAL(kind=dp), POINTER ::	IniBen(:,:) 
       REAL(kind=dp) :: IniBen2(2),Tepsilon_1(3),Tepsilon_2(3),Tbending_1(3),Tbending_2(3)		
       LOGICAL :: Found
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), &
           Curvature(3,MaxDofs),InPlaneStrain(3,MaxDofs),Gammaa(2,MaxDofs),& 
           KappaVector(3), EPSILO(3), TShear(2), T5(5,5), Tmat(3,3), &
           T0(3,3), T0S(2,2), Omega(1,MaxDofs), Gdrilling(1,1), &
           GammaVector(2), EpsMaterial(3), KappaMaterial(3), ToLayer(5,5),InvT0(3,3),AMatrix(3,3), BMatrix(3,3),DMatrix(3,3),& 
           AStarmatrix(2,2),LAMDCMatrix(8,8), NVec(3), MVec(3)
       REAL( KIND=dp ) :: detJ, U, V, W, Kappa, Q5(5,5), &
           h, Transformation(3,3), CopyOfNodes(3,MaxNodes), &
           delta, ConstantLoadVec(3), VariableLoadVec(3), Puvw
       
       REAL( KIND=dp ) :: Tblock(6*n,6*n)
       REAL(kind=dp) :: TShearStress(2)
       
       LOGICAL :: Stat, Found1
       INTEGER :: i, p, q, k ,j,t
! -------------------------------------------------------------
       REAL(KIND=dp) :: z1, z2, theta, ct, st, SCF, xn, yn, zn
       REAL(KIND=dp) :: side1(2), side2(2)
       REAL(KIND=dp) :: dWdx, dWdy
       REAL(KIND=dp) :: LV1(3,1), LV2(3,1), LV3(3,1), dUdx(3,2), dRdx(3,2)
! --------------------------------------------------------------
       BodyForce => GetBodyForce()
       
       OutNormal = GetConstReal(BodyForce,'Direction Of Outward Normal',Found1)
       IF(.NOT.Found1)THEN
         OutNormal=1.0
       ENDIF
       
       Curvature       = 0.0d0
       InPlaneStrain   = 0.0d0
       TShear = 0.0d0

!      The transformation Xglob -> Xloc is the transpose of the local  !	basis:
!      --------------------------------------------------------
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
!	This 'Rechange' has been done with a view to convert here from 
!	local Thetas to Betas in stead of previous case of conversion from Betas to Thetas.
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

       IntegStuff = GaussPoints( Element, 1 )

!      Numerical integration:
       !      ----------------------
       DO t = 1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
!         S = IntegStuff % s(t)

         !        Basis function values & derivatives at the integration point:
         !        -------------------------------------------------------------
         stat = ElementInfo( Element, Nodes, U, V, W, &	!Nodes changed to ElementNodes31.12.13 -Reverted
             detJ, Basis, dBasisdx )
         IF( ElementNumber == 15 ) THEN
           
           WRITE(Message,'(a,I6,8f12.6)') ' Int Pt No, u, v, w, detJ, dbasisdx(2,1)+3 Basis', &
               t,U,V,W,detJ,dbasisdx(2,1),Basis(1),Basis(2),Basis(3)
           CALL Info('ShellSolver',Message)
         ENDIF
         !         S = S * detJ
       ENDDO
!	**********Testing on 22.3.15 ends*********************
!      Basis function values & derivatives in the stress eval-point:
!      -------------------------------------------------------------
!       stat = ElementInfo( Element, Nodes, U, V, W, &
 !                   detJ, Basis, dBasisdx )

!      Material parameters in the stress evaluation point:
!      ---------------------------------------------------
       h = SUM( NodalThickness(1:n) * Basis(1:n) )	!changed to NodalThickness fron Thickness15.10.15
       Puvw = SUM( Poisson(1:n)  * Basis(1:n) )	
       
       
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


!	On multiplyiny curvature by h/2, it is now Max. strain at top edge.
       Curvature=OutNormal*Curvature*h/2.0d0	

       IF(TopSideStress)Curvature = Curvature
       IF(BottomSideStress)Curvature = -Curvature
!	This shear strain section needs to be completely rewritten in terms of 
!	Theta-x and Theta-y (Basically Betas) and dowerW/Doverx , DoverW/Dovery.
       Gammaa=0.0	
!      Shear strains (transversal):
!      ----------------------------
!       CALL CovariantInterpolation( Gammaa, Basis, &		15.3.15
 !           Nodes % x(1:n), Nodes % y(1:n), U, V, n )

!       CALL ShearCorrectionFactor( SCF, h, Nodes % x(1:n), &		15.3.15
 !           Nodes % y(1:n), n, StabParam1 )

       DO p = 1,n
!          Gammaa(1:2,6*p-3) = dBasisdx(p,1:2)	!TRANSVERSE SHEAR is not considered 21.3.15
!	This adds the rotation component to obtain both transverse shear strains 15.3.15
!	All the other procedure of transverse shear strain discarded 15.3.15
!	NOTE in place of Theta-x and Theta-y, -Betas are used due to transformation of localDeflection
!	in the beginning of this subroutine
!		  Gammaa(1,6*p-2)=Basis(p)
!		  gammaa(2,6*p-2)=Basis(p)
!	************Addition 15.3.15 ends***************************
       END DO
!	********************Inserted 14.3.15************************************
!		Gammaa=Gammaa*h/2.0d0	!To check thoroughly (h/2 inserted) ! not necessary 15.3.15
!	********************Insertion 14.3.15 ends*******************************

!      Drilling DOFs (in-plane rotations):
!      -----------------------------------
       DO p = 1,n
         Omega(1,6*p-5) =  dBasisdx(p,2)/2.0d0  !  u_{x,y}
         Omega(1,6*p-4) = -dBasisdx(p,1)/2.0d0  ! -u_{y,x}
         Omega(1,6*p-0) =  Basis(p)             ! rotation
       END DO
!	Drilling rotation may not contribute to strain . This is not used later(15.3.14)
!	All three terms stated above imply in-plane roration, in effect. 15.3.15
!	*************Subtracting Initial (Th) Strain and Initial Bending 7.10.15********** 
	    !	Weld Distortion Simulation and Thermal Strain
       WeldSimu = GetLogical( BodyForce, 'Weld Distortion Simulation', GotIt )
       IF( .NOT. GotIt ) WeldSimu = .FALSE.
       ThStrain = GetLogical( BodyForce, 'Calculate Thermal Strain', GotIt )
       IF( .NOT. GotIt ) ThStrain = .FALSE.
       IF( WeldSimu .AND. ThStrain )THEN
         CALL Fatal( 'ShellSolver', 'Weld Distortion Simulation and Calculate Thermal Strain '// &
             'cannot be specified simultaneously' )
       ENDIF
       AnisoPlane = GetLogical( BodyForce, 'Use Relaxed AnisoPlane Testing', GotIt )

       Tepsilon_2=0.0
       Tbending_2=0.0
       IniBen2=0.0
!	******Following shifted from BulkAssembly subroutine to this place (to be more logical)******
!	   (Should be modified for actual temperature distribution)
       TT(1:n) = GetReal( BodyForce, 'Temperature', GotIt )
       IF((.NOT.GotIt).AND.( ThStrain .OR. WeldSimu ) ) THEN	
         !       WRITE( Message,'(A,I9)' ) &
!	   'Temperature not specified for thermal strain or weld simulation calculations Taking it zero'
!       CALL Info('ShellSolve',Message)
         TT(1:n) = 0.0d0
       ENDIF
       ThermalAnisotropy=.FALSE.
	!	If Anisotropy is specified then only parameters related to it are read and axes found
       Aniso1 =GetString( BodyForce, 'Thermal Anisotropy', Found)
       IF(Found)THEN
         ThermalAnisotropy=.TRUE.
         !			WRITE( Message,'(A,L2)' ) 'ThermalAnisotropy ',thermalanisotropy	
         !			CALL warn('ShellSolve',Message)

         IF((ThStrain).OR.(WeldSimu))THEN
           CALL GetAnisoParameter(Aniso,Point1,Point2,alfa,Te,IniBen2)	
           CALL GetAnisotropyAxes(Element,Aniso,Point1, Point2,cosTh,sinTh)
           Tmatrix= ThStrainConversion(cosTh,sinTh)
         ENDIF
         !	Te is the Weld Simulation Temperature, read in sub GetAnisoParameter()
         IF(WeldSimu)TT(1:n)=Te			
       ELSE IF(ThStrain)THEN
         !	Get alfa1 for calculations of non-anisotropic thermal strain
         alfa(1) = GetConstReal(BodyForce,'Alfa 1',Found1)
         IF(.NOT.Found1) THEN
           WRITE( Message,'(A)' ) 'Alfa 1 not specified for thermal strain calculations '
           CALL warn('ShellSolve',Message)
           CALL Fatal('ShellSolve','STOP')
         ENDIF
       ENDIF					

       Te=SUM( TT(1:n) * Basis(1:n) )
       i=6*n
       IF(thermalanisotropy)THEN
         Tepsilon_1(1)= alfa(1)*Te
         Tepsilon_1(2) = alfa(2)*Te
         Tepsilon_1(3)=0.0d0
         Tepsilon_2=MATMUL(Tmatrix,Tepsilon_1)	
         !See remark in description of Function ThStrainConversion( )
       ELSE
         Tepsilon_2(1)=alfa(1)*Te
         Tepsilon_2(2)=alfa(1)*Te
         Tepsilon_2(3)=0.0d0
       ENDIF
       !	Here thickness is already incorporated in Amatrix
       !	Considering Initial Bending strain  
       !	********WARNING-IniBen2( ) for routine 3D Thermal case is to be treated separately*************
       IF(thermalanisotropy)THEN
         Tbending_1(1) = IniBen2(1)
         Tbending_1(2) = IniBen2(2)
         Tbending_1(3) = 0.0
         Tbending_2=MATMUL(Tmatrix,Tbending_1) !See remark in description of Function ThStrainConversion( )
       ELSE
         Tbending_2(1) = IniBen2(1)
         Tbending_2(2) = IniBen2(2)
         Tbending_2(3) =0.0
       ENDIF
!	Following has been done so that the input IniBen2 is interpreted properly.
!	Tbending_2=Tbending_2*OutNormal
!	This has been removed in view of the direction of curvature being correct due to the 
!	multiplication by OutNorma above and so no further change of input direction of IniBen2 
!	is required . (See Note 7 in file Readme-2_ShellMultiSolver2)

!      Venymä- ja käyristymävektorit lokaalissa koordinaatistossa:
!      -----------------------------------------------------------
!	Subtracting Thermal strain and Initial bending 
       KappaVector = MATMUL( Curvature(1:3,1:6*n), LocalDeflection(1:6*n) )
       EPSILO = MATMUL( InPlaneStrain(1:3,1:6*n),  LocalDeflection(1:6*n) )

       !WRITE(Message,'(a,6e12.2)') ' Kappavector and EPSIL0', (kappavector(i),i=1,3),(EPSILO(i),i=1,3)
!        CALL Info('ShellSolver',Message)
       IF(TopSideStress) Tbending_2 = Tbending_2*h/2.0
       IF(BottomSideStress) Tbending_2 = -Tbending_2*h/2.0

       !	NOTE (15.8.16)- Since KappaVector is determined on the basis of the value of OutNormal 
!	(see above in this subroutine while calculating curvature), so the Top Side is already
!	specified as the one on the basis of which the value of OutNormal is determined. This means
!	that if the parameter Tbending_2 is determined on the basis of the presumed + direction of 
!	the normal to the surface, irrespective of the node numbering in the element 
!	(anticlockwise or otherwise), the magnitude Tbending_2 will be correct in its sign and so
!	it should be subtracted from KappaVector to accommodate initial strain due to bending, as 
!	has been done correctly below in the program.

       KappaVector = KappaVector-Tbending_2(1:3)
       !See definition curvature somewhere above, where it is transformed to maximum at top edge.
       EPSILO = EPSILO-Tepsilon_2(1:3)
!WRITE(Message,'(a,6e12.2)') '  Modified Kappavector and EPSIL0', &
!(kappavector(i),i=1,3),(EPSILO(i),i=1,3)
!        CALL Info('ShellSolver',Message)

       GammaVector = MATMUL( Gammaa(1:2,1:6*n), LocalDeflection(1:6*n) )
!	GammaVector is now zero because Gammaa is initialized to zero above 8.10.15
!	However, when transverse shear strain section is to be used this statement will be useful
!	**********************************
!      VonKarman strains:         
! ========================

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
! ------------------------------------
         DO p = 1,n
           LV1(1,1) = LV1(1,1) + LocalDeflection(6*(p-1)+1) * dBasisdx(p,1)
           LV1(2,1) = LV1(2,1) + LocalDeflection(6*(p-1)+2) * dBasisdx(p,2)
           LV1(3,1) = LV1(3,1) + LocalDeflection(6*(p-1)+1) * dBasisdx(p,2)& 
               + LocalDeflection(6*(p-1)+2) * dBasisdx(p,1) 

           LV2(1,1) = LV2(1,1) + LocalDeflection(6*(p-1)+4) * dBasisdx(p,1)
           LV2(2,1) = LV2(2,1) + LocalDeflection(6*(p-1)+5) * dBasisdx(p,2)
           LV2(3,1) = LV2(3,1) + LocalDeflection(6*(p-1)+4) * dBasisdx(p,2)& 
               + LocalDeflection(6*(p-1)+5) * dBasisdx(p,1)  
         END DO

! Non-linear part of strain and curvature
! ----------------------------------------
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

       Amatrix=Amatrix/h
       AstarMatrix=Astarmatrix/h

       Gdrilling(1,1) = StabParam2*(Astarmatrix(1,1)+Astarmatrix(2,2))
!	Drilling rotation may not contribute to shear strain! 
!      Normaalivoima- ja momenttivektorit (per pituusyksikkö) 
!	lokaalissa koord.:
!      -----------------------------------------------------------
       NVec = 0.0d0
       Mvec = 0.0d0
       NVec = MATMUL( Amatrix, EPSILO ) + MATMUL( Bmatrix, KappaVector )
!      MVec = MATMUL( Bmatrix, EPSILO ) + MATMUL( Dmatrix, KappaVector )	! changed 14.3.15
!		The change was necessary because during the drivation of stress at top or bottom edge
!		it turns out that the multiplier ([D] matrix) is to be same as used in case of plane stress 
       MVec = MATMUL( Bmatrix, EPSILO ) + MATMUL( Amatrix, KappaVector )
       TShearStress=MATMUL(Astarmatrix, GammaVector )
!      Lopuksi vektorit tensoreiksi plus transformaatiot globaaliin 
!	koordinaatistoon:
!      ---------------------------------------------------------
       Eps = 0.0d0
       Eps(1,1) = EPSILO(1)
       Eps(2,2) = EPSILO(2)
       Eps(1,2) = EPSILO(3) / 2.0d0
       Eps(2,1) = EPSILO(3) / 2.0d0
!	Thickness strain in plane stress is -(Nu/(1.0-Nu))(Epsilo1+Epsilo2) 
       Eps(3,3)=-(Puvw/(1.0-Puvw))*(Eps(1,1)+Eps(2,2))	

       Eps = MATMUL( TRANSPOSE(Transformation), Eps )	!TO CHECK 14.3.15 Checked Correct Sokolnikoff
       Eps = MATMUL( Eps, Transformation )				!TO CHECK 14.3.15
       
       Kap = 0.0d0
       Kap(1,1) = KappaVector(1)
       Kap(2,2) = KappaVector(2)
       Kap(1,2) = KappaVector(3)/2.0d0
       Kap(2,1) = KappaVector(3)/2.0d0
!	Added here after removing from Eps 18.3.15************! GammaVector is made zero  22.3.15
       Kap(1,3) = GammaVector(1) / 2.0d0	!This section will be useful if transverse shear strain is to be considered
       Kap(2,3) = GammaVector(2) / 2.0d0	!At present it is zero
       Kap(3,1) = GammaVector(1) / 2.0d0
       Kap(3,2) = GammaVector(2) / 2.0d0
!	Addition here 18.3.15 ends******************************
       Kap = MATMUL( TRANSPOSE(Transformation), Kap )	!TO CHECK 14.3.15	Checked correct Sokolnikoff
       Kap = MATMUL( Kap, Transformation )				!TO CHECK 14.3.15

       NTen = 0.0d0
       NTen(1,1)=NVec(1)
       NTen(1,2)=NVec(3)
       NTen(2,1)=NVec(3)
       NTen(2,2)=NVec(2)

       NtenMaterial = 0.0d0
       NtenMaterial(1:2,1:2) = Nten(1:2,1:2)

       NTen = MATMUL( TRANSPOSE(Transformation), NTen )
       NTen = MATMUL( NTen, Transformation )			! "

       MTen = 0.0d0
       MTen(1,1)=MVec(1)
       MTen(1,2)=MVec(3)
       MTen(2,1)=MVec(3)
       MTen(2,2)=MVec(2)
	!	Added transverse shear strains 
       MTen(3,1)=TShearStress(1)
       MTen(3,2)=TShearStress(2)	! shear stresses are not halved as strain 
       MTen(1,3)=TShearStress(1)	! TshearStress is zero because of Gammaa and GammaVecor being zero
       MTen(1,3)=TShearStress(2)

       MtenMaterial = 0.0d0
       MtenMaterial(1:2,1:2) = Mten(1:2,1:2)
       
       MTen = MATMUL( TRANSPOSE(Transformation), MTen )
       MTen = MATMUL( MTen, Transformation )			!  "

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
! ====================================================================


     SUBROUTINE AveragingWeights3( Nodes, Weight3 )
! ---------------------------------------------------------------
       TYPE( Nodes_t ) :: Nodes
       REAL( KIND=DP ) :: Weight3(:)
! -----------------------------------------------------------
       INTEGER :: i
       REAL( KIND=DP ) :: Side1(2), Side2(2)
! ------------------------------------------------------------------
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
! ----------------------------------------------------------------
     END SUBROUTINE AveragingWeights3


     SUBROUTINE AveragingWeights4( Nodes, Weight4 )
! ----------------------------------------------------------------
       TYPE( Nodes_t ) :: Nodes
       REAL( KIND=DP ) :: Weight4(:)
! --------------------------------------------------------------
       INTEGER :: i
       REAL( KIND=DP ) :: Side1(2), Side2(2)
! ------------------------------------------------------------
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
! -----------------------------------------------------------------
     END SUBROUTINE AveragingWeights4


! ---------------------------------------------------------------
     SUBROUTINE SwitchToLocal( Nodes, CopyOfNodes, Transformation, n )
! ------------------------------------------------------------------
       REAL(KIND=dp) :: Transformation(:,:), CopyOfNodes(:,:)
       TYPE(Nodes_t) :: Nodes
       INTEGER :: n
! -------------------------------------------------------------------
       REAL(KIND=dp) :: XYZGlobal(3,4), XYZLocal(3,4)
! -------------------------------------------------------------------
       CopyOfNodes(1,1:n) = Nodes % x(1:n)
       CopyOfNodes(2,1:n) = Nodes % y(1:n)
       CopyOfNodes(3,1:n) = Nodes % z(1:n)

       XYZGlobal(1,1:n)  = CopyOfNodes(1,1:n) - SUM( CopyOfNodes(1,1:n) )/ n
       XYZGlobal(2,1:n)  = CopyOfNodes(2,1:n) - SUM( CopyOfNodes(2,1:n) ) / n
       XYZGlobal(3,1:n)  = CopyOfNodes(3,1:n) - SUM( CopyOfNodes(3,1:n) ) / n

       XYZLocal(1:3,1:n) = MATMUL( Transformation(1:3,1:3), XYZGlobal(1:3,1:n) )

       Nodes % x(1:n) = XYZLocal(1,1:n)
       Nodes % y(1:n) = XYZLocal(2,1:n)
       Nodes % z(1:n) = XYZLocal(3,1:n)
! ----------------------------------------------------------------
     END SUBROUTINE SwitchToLocal
! -------------------------------------------------------------------

! ====================================================================

! ---------------------------------------------------------------
     FUNCTION LocalBasis( Nodes, n ) RESULT( BasisVectors )	!Modified (same results) 20.11.13
! ----------------------------------------------------------------
       TYPE(Nodes_t) :: Nodes
       REAL(KIND=dp) :: BasisVectors(3,3)
       INTEGER :: n
! -----------------------------------------------------------
       REAL(KIND=dp) :: Tangent1(3), Tangent2(3), Tangent3(3)
! -------------------------------------------------------------

!      First, find a couple of in-plane unit vectors:
!      ----------------------------------------------
!      First, find a couple of in-plane unit vectors:
!      ----------------------------------------------
       Tangent1(1) = Nodes % x(2) - Nodes % x(1)
       Tangent1(2) = Nodes % y(2) - Nodes % y(1)
       Tangent1(3) = Nodes % z(2) - Nodes % z(1)
 !      Tangent1 = Tangent1 / SQRT( SUM( Tangent1**2  ) )

       Tangent2(1) = Nodes % x(3) - Nodes % x(2)
       Tangent2(2) = Nodes % y(3) - Nodes % y(2)
       Tangent2(3) = Nodes % z(3) - Nodes % z(2)
 !      Tangent2 = Tangent2 / SQRT( SUM( Tangent2**2  ) )
       Tangent3 = CrossProduct1( Tangent1, Tangent2 )
       Tangent2 = CrossProduct1( Tangent3, Tangent1 )
       Tangent1 = Tangent1 / SQRT( SUM( Tangent1**2  ) )
       Tangent2 = Tangent2 / SQRT( SUM( Tangent2**2  ) )
       Tangent3 = Tangent3 / SQRT( SUM( Tangent3**2  ) )
       
!      Then, define the local cartesian unit basis vectors:
!      ----------------------------------------------------
       BasisVectors(1:3,1) = Tangent1
       BasisVectors(1:3,2) = Tangent2
       BasisVectors(1:3,3) = Tangent3
!      ----------------------------------------------------
!       BasisVectors(1:3,1) = Tangent1
!BasisVectors(1:3,2) = Tangent2 - SUM( Tangent1 * Tangent2 ) * Tangent1
!BasisVectors(1:3,2) = BasisVectors(1:3,2) / SQRT( SUM(  BasisVectors(1:3,2)**2 ) )
!BasisVectors(1:3,3) = CrossProduct1( BasisVectors(1:3,1), BasisVectors(1:3,2) )
! -------------------------------------------------------------      
     END FUNCTION LocalBasis
! -------------------------------------------------------------       


! -----------------------------------------------------------       
     SUBROUTINE IsotropicElasticity(Ematrix, &
         Gmatrix,Poisson,Young,Thickness,Basis,n)
! -------------------------------------------------------------
       REAL(KIND=dp) :: Ematrix(:,:), Gmatrix(:,:), Basis(:)
       REAL(KIND=dp) :: Poisson(:), Young(:), Thickness(:)
       REAL(KIND=dp) :: Euvw, Puvw, Guvw, Tuvw
       INTEGER :: n
       ! ------------------------------------------------------------------
       Euvw = SUM( Young(1:n)    * Basis(1:n) )
!	**********************Inserted on 4.11.14*************************
!	if(MembraneOnly)Euvw=Euvw*0.2
!	**********************Insertion 4.11.14 ends***********************
       Puvw = SUM( Poisson(1:n)  * Basis(1:n) )
       Tuvw = SUM( Thickness(1:n)* Basis(1:n) )
       Guvw = Euvw/(2.0d0*(1.0d0 + Puvw))
       
       Ematrix = 0.0d0
       Ematrix(1,1) = 1.0d0
       Ematrix(1,2) = Puvw
       Ematrix(2,1) = Puvw
       Ematrix(2,2) = 1.0d0
       Ematrix(3,3) = (1.0d0-Puvw)/2.0d0
       Ematrix = Ematrix * Euvw * (Tuvw**3) / (12.0d0 * (1.0d0 -  Puvw*Puvw))

       Gmatrix = 0.0d0
       Gmatrix(1,1) = Guvw*Tuvw
       Gmatrix(2,2) = Guvw*Tuvw
! -----------------------------------------------------------------
     END SUBROUTINE IsotropicElasticity
!  ---------------------------------------------------------

!  ===================================================================

! -----------------------------------------------------------------
     SUBROUTINE IsotropicInPlaneElasticity( Ematrix, &
         Poisson, Young, Thickness, Basis, n )
! ----------------------------------------------------------------
       REAL(KIND=dp) :: Ematrix(:,:), Basis(:)
       REAL(KIND=dp) :: Poisson(:), Young(:), Thickness(:)
       REAL(KIND=dp) :: Euvw, Puvw, Guvw, Tuvw
       INTEGER :: n
! --------------------------------------------------------------
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
! --------------------------------------------------------------
     END SUBROUTINE IsotropicInPlaneElasticity
! -----------------------------------------------------------------

! ==================================================================

! ------------------------------------------------------------
     SUBROUTINE ShearCorrectionFactor(Kappa,Thickness,x,y,n,StabParam)
! ------------------------------------------------------------
       REAL(KIND=dp) :: Kappa,Thickness,x(:),y(:),StabParam
       INTEGER :: n
! ---------------------------------------------------------------
       REAL(KIND=dp) :: x21,x32,x43,x13,x14,y21,y32,y43,y13,y14, &
           l21,l32,l43,l13,l14,alpha,h
! -----------------------------------------------------------------
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
         CALL Fatal('ShellSolver',&
             'Illegal number of nodes for Smitc elements')
       END SELECT
       ! -----------------------------------------------------------------
     END SUBROUTINE ShearCorrectionFactor
! --------------------------------------------------------------

!  ====================================================================

! ----------------------------------------------------------------
     SUBROUTINE AddEnergy(A,B,C,m,n,s)
! -----------------------------------------------------------------
!      Performs the operation
!
!         A = A + C' * B * C * s
!
!      with
!
!         Size( A ) = n x n
!         Size( B ) = m x m
!         Size( C ) = m x n
!  --------------------------------------------------------------
       REAL(KIND=dp) :: A(:,:),B(:,:),C(:,:),s
       INTEGER :: m,n
! -----------------------------------------------------------------
       INTEGER :: i,j,k,l
! ---------------------------------------------------------------
       DO i=1,n
         DO j=1,n
           DO k=1,m
             DO l=1,m
               A(i,j) = A(i,j) + C(k,i)*B(k,l)*C(l,j) * s
             END DO
           END DO
         END DO
       END DO
! -------------------------------------------------------------------
     END SUBROUTINE AddEnergy
! ----------------------------------------------------------------

! ================================================================

! -----------------------------------------------------------------
     SUBROUTINE AddInnerProducts(A,B,C,D,m,n,s)
! ------------------------------------------------------------------
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
! ------------------------------------------------------------------
       REAL(KIND=dp) :: A(:,:),B(:,:),C(:,:),D(:,:),s
       INTEGER :: m,n
! ------------------------------------------------------------------
       INTEGER :: i,j,k,l
! ------------------------------------------------------------------
       DO i=1,n
         DO j=1,n
           DO k=1,m
             DO l=1,m
               A(i,j) = A(i,j) + C(i,k)*B(k,l)*D(l,j) * s
             END DO
           END DO
         END DO
       END DO
! ----------------------------------------------------------------
     END SUBROUTINE AddInnerProducts
! ---------------------------------------------------------------

! ====================================================================

! ---------------------------------------------------------------
     SUBROUTINE CovariantInterpolation(ShearStrain,Basis,X,Y,U,V,n)
! -------------------------------------------------------------------
       REAL(KIND=dp) :: ShearStrain(:,:),Basis(:),X(:),Y(:),U,V
       INTEGER :: n
! ------------------------------------------------------------------
       REAL(KIND=dp) :: detJ,Jmat(2,2),invJ(2,2),ShearRef(2,100)
       REAL(KIND=dp) :: Tau(2),Sdofs(100)
       INTEGER :: j

       SELECT CASE(n)

!      The SMITC3 element
!      ==================
       CASE(3)
!	****************************************************
!	Only CASE(3) slightly changed to make it compatible with the 
!	definition of BetaX =-du/dx and betaY=-dv/dy	(15.10.15)
!	****************************************************
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
!	********************************************************************
!	To convert it for use in solver with BetaX =-du/dx and BetaY=-dv/dy
!	the final quantity should be multiplied by -1 as below 15.10.15
         ShearStrain=-ShearStrain	!Added 15.10.15
!	********************************************************************
!      The SMITC4 element 
!	Untouched in so far as the chane in definition of Beta is concerned
!	So it may produce unreliable results later for SMITC4 element
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
! ------------------------------------------------------------------
     END SUBROUTINE CovariantInterpolation
! ----------------------------------------------------------------

! ===================================================================

! -------------------------------------------------------------------
     SUBROUTINE Jacobi3(Jmat,invJ,detJ,x,y)
       ! ---------------------------------------------------------------
       REAL(KIND=dp) :: Jmat(:,:),invJ(:,:),detJ,x(:),y(:)
       ! -------------------------------------------------------------------
       Jmat(1,1) = x(2)-x(1)
       Jmat(2,1) = x(3)-x(1)
       Jmat(1,2) = y(2)-y(1)
       Jmat(2,2) = y(3)-y(1)

       detJ = Jmat(1,1)*Jmat(2,2)-Jmat(1,2)*Jmat(2,1)

       invJ(1,1) =  Jmat(2,2)/detJ
       invJ(2,2) =  Jmat(1,1)/detJ
       invJ(1,2) = -Jmat(1,2)/detJ
       invJ(2,1) = -Jmat(2,1)/detJ
       ! ------------------------------------------------------------------
     END SUBROUTINE Jacobi3
! -------------------------------------------------------------------

! =====================================================================

! --------------------------------------------------------------------
     SUBROUTINE Jacobi4(Jmat,invJ,detJ,xi,eta,x,y)
! --------------------------------------------------------------------
       REAL(KIND=dp) :: Jmat(:,:),invJ(:,:),detJ,xi,eta,x(:),y(:)
! ------------------------------------------------------------------
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
       ! --------------------------------------------------------------------
     END SUBROUTINE Jacobi4
! -----------------------------------------------------------------

     SUBROUTINE GetAnisoParameter(Aniso,Point1,Point2,alfa,Te,IniBen2)	!Inserted on 14.11.13
!	This subroutine reads the thermal anisotropy values and related quantities
!	It also reads quantities related to Weld Distortion Simulation.
       CHARACTER(LEN=15) :: Aniso	
       !	Real(kind=dp), pointer ::	IniBen(:,:) !IniBen(:,:) inserted on 24.9.15
       REAL(Kind=dp) :: Point1(3), Point2(3),Te,alfa(2), IniBen2(2) !IniBen2 added 29.9.15
       LOGICAL Found,Found1,Found2,Found3,Found4,Found5,Found6,WeldSimu

       Aniso =GetString( BodyForce, 'Thermal Anisotropy', Found)
       IF (Found) THEN
!	WRITE(Message,'(A,A)') &
!	'Thermal anisotropy type is ', Aniso 
!	CALL Info('Shell Solver',Message)

         Point1(1) = GetConstReal(BodyForce,'Point1 1',Found1)
         Point1(2) = GetConstReal(BodyForce,'Point1 2',Found2)
         Point1(3) = GetConstReal(BodyForce,'Point1 3',Found3)
         Point2(1) = GetConstReal(BodyForce,'Point2 1',Found4)
         Point2(2) = GetConstReal(BodyForce,'Point2 2',Found5)
         Point2(3) = GetConstReal(BodyForce,'Point2 3',Found6)

         IF (Found1.AND.Found2.AND.Found3.AND.Found4.AND.Found5.AND.Found6) THEN
           CONTINUE
!	WRITE(Message,'(A,3f10.4,A,3f10.4)') &
!	'Coordinates of Points are ', (Point1(i),i=1,3),' and ',(Point2(i),i=1,3)
!	CALL Info('Shell Solver',Message)
         ELSE
           WRITE(Message,'(A)') 'Check coordinates of Points  '
           CALL Warn('Shell Solver',Message)
           CALL FATAL('Shell Solver','STOP')           
         ENDIF
         
         alfa(1) = GetConstReal(BodyForce,'Alfa 1',Found1)
         alfa(2) = GetConstReal(BodyForce,'Alfa 2',Found2)

         IF (Found1.AND.Found2) THEN
           CONTINUE	!Added24.6.15
!	********Removed 24.6.15*********************
!	WRITE(Message,'(A,2E12.4)') &
!	'Alfa 1 and Alfa 2  are ', alfa(1),alfa(2)
!	CALL Info('Shell Solver',Message)
!	***********Removal24.6.15 ends***************
         ELSE
           WRITE(Message,'(A)') 'Two anisotropic thermal expansion coefficients not specified  '
           CALL Warn('Shell Solver',Message)
           CALL FATAL('Shell Solver','STOP')
         ENDIF
         
       ELSE
         Aniso ='NULL'
       END IF
       WeldSimu = GetLogical(BodyForce, 'Weld Distortion Simulation', Found1 )
       IF(WeldSimu)THEN
         Te = GetConstReal(BodyForce,'Temperature',Found)
         IF(Found)THEN
           CONTINUE
           !	WRITE(Message,'(A,f10.4)') 'Temperature for weld distortion simulation is   ',Te
!	CALL Info('Shell Solver',Message)

         ELSEIF(.NOT.Found)THEN
           
           WRITE(Message,'(A)') 'Temperature for weld distortion simulation not specified. Taking zero '
           CALL Warn('Shell Solver',Message)
           Te = 0.0
         ENDIF
         !	*********added test 24.9.15 (Changed 6.10.15)********************
         IniBen2(1) = GetConstReal(BodyForce,'Initial Bending 1',Found1)
         IF(.NOT. Found1)IniBen2(1)=0.0d0
         IniBen2(2) = GetConstReal(BodyForce,'Initial Bending 2',Found2)
         IF(.NOT.Found2)IniBen2(2)=0.0d0
!		call GetConstRealArray(BodyForce,IniBen,'Initial Bending',Found1)
!		if(.not.Found1)then
!		IniBen=0.0d0
!		WRITE( Message,'(A,I9)' ) 'Initial Bending not found. Taking it zero'
!		CALL Info('ShellSolve',Message)
!		else
!		endif
!		IniBen2(1)=IniBen(1,1)
!		IniBen2(2)=IniBen(2,2)
!		WRITE( Message,'(A, 6e12.4)' ) 'IniBen(2,2) ',IniBen(1,1),IniBen(1,2),IniBen(2,1),IniBen(2,2)
!		CALL Info('ShellSolve',Message)
			

       ENDIF
     END SUBROUTINE GetAnisoParameter

     !-----------------------------------------------------------------------
     SUBROUTINE GetAnisotropyAxes(Element,AnisoType,A,B,cosTh,sinTh)	
!	Added on 18.11.13
!	This subroutine calculates the direction cosines of local anisotropy x-axis with
!	respect to the elemental local axes. Anisotropy type could be Cartesian, Circular, Cylindrical, 
!	Spherical or Conical. Local anisotropy x-axis need not be along the length of weld for such case
!	but the data should be interpreted accordingly.
       Character(LEN=15) :: AnisoType
       Integer :: ThisBody
       Type(Element_t), POINTER :: Element
       TYPE(Nodes_t)   :: Nodes
       Real(kind=dp) :: A(3),B(3),cosTh,sinTh, Tangent1(3), Tangent2(3), Axis1(3),Axis2(3),Axis3(3), &
           Vect1(3), Vect2(3),Vect3(3), T1Mag, T2Mag, Vect1Mag, Vect2Mag,Vect3Mag, CP1(3), CP2(3), &
           Centroid(3),Err, Axis1Mag,Axis2Mag, ErrBasic
       !	Axis1,Axis2 and Axis3 are the three vectors representing three local axes of elements	
       CALL GetElementNodes(Nodes )

       IF( AnisoPlane ) THEN
         ErrBasic=1.0e-04
       ELSE
         ErrBasic=1.0e-08
       ENDIF
       ThisBody = Element % BodyId		!Added 9.9.16

!      First, find in-plane x axis and second edge of element as vectors:
!      ----------------------------------------------
       Tangent1(1) = Nodes % x(2) - Nodes % x(1)
       Tangent1(2) = Nodes % y(2) - Nodes % y(1)
       Tangent1(3) = Nodes % z(2) - Nodes % z(1)
       T1Mag = SQRT( SUM( Tangent1**2  ) )
       Tangent2(1)= Nodes % x(3) - Nodes % x(2)
       Tangent2(2) = Nodes % y(3) - Nodes % y(2)
       Tangent2(3) = Nodes % z(3) - Nodes % z(2)
       T2Mag=SQRT( SUM( Tangent2**2  ) )
       !	 WRITE(Message,'(A,A,E14.4)') 'Value of T2Mag is  ',AnisoType,T2Mag
       !     CALL Info('ShellSolve',Message)
       Axis1 = Tangent1
       Axis3 = CrossProduct1(Tangent1,Tangent2)
       Axis2 = Crossproduct1(Axis3,Axis1)
       Axis1Mag = SQRT(SUM(Axis1**2))
       Axis2Mag = SQRT(SUM(Axis2**2))
       IF((AnisoType == 'cartesian').OR.(Anisotype == 'cylindrical'))THEN
         !	Axis x is directed from origin (Point1-A) to Point2-B
         !	x-axis of cylindrical surface is directed from Point1-A to Point2-B
         
         !	Check whether This axis is in-planar with element
         
         Vect1(1) = B(1)-A(1)
         Vect1(2) = B(2)-A(2)
         Vect1(3) = B(3)-A(3)
         Vect1Mag = SQRT( SUM(Vect1**2))
         !	WRITE(Message,'(A,E14.4)') 'Value of Vect1Mag is  ',Vect1Mag
         !     CALL Info('ShellSolve',Message)

         CP1 = CrossProduct1(Vect1,Tangent1)/(Vect1Mag*T1Mag)
         CP2 = CrossProduct1(Vect1,Tangent2)/(Vect1Mag*T2Mag)
         Vect2 = CrossProduct1(CP1,CP2)					!Vect2 used as dummy
         Err = SUM(Vect2*vect2)
         !	 WRITE(Message,'(A,E14.4)') 'Value of Err is  ',Err
         !    CALL Info('ShellSolve',Message)
         IF(Err**2 <ErrBasic)THEN
           !	Vect2=CrossProduct1(Tangent1, Vect1)/(T1Mag*Vect1Mag)	!Vect2 used as dummy
           !	sinTh=SQRT(SUM(Vect2*Vect2))
           sinTh=SUM(Vect1*Axis2)/(Vect1Mag*Axis2Mag)
           cosTh=SUM(Axis1*Vect1)/(Axis1Mag*Vect1Mag)
         ELSE 
           WRITE(Message,'(A,A,A,I9)') 'x-axis is not in-planar with element for Anisotropy ',&
               AnisoType, 'for Body ',ThisBody
           CALL Warn('ShellSolve',Message)
           CALL Fatal('ShellSolve','STOP')
         ENDIF

       ELSE IF(AnisoType == 'circular')THEN
         !	Y-axis is directed from Center (Point1-A) to Centroid of element.
         Centroid(1) = SUM(Nodes%x(1:n))/n
         Centroid(2) = SUM(Nodes%y(1:n))/n
         Centroid(3) = SUM(Nodes%z(1:n))/n
         Vect1(1) = Centroid(1)-A(1)
         Vect1(2) = Centroid(2)-A(2)
         Vect1(3) = Centroid(3)-A(3)
         Vect1Mag = SQRT( SUM(Vect1**2))
         CP1 = CrossProduct1(Vect1,Tangent1)/(Vect1Mag*T1Mag)
         CP2 = CrossProduct1(Vect1,Tangent2)/(Vect1Mag*T2Mag)
         Vect2 = CrossProduct1(CP1,CP2)
         Err = SUM(Vect2*Vect2)
         IF(Err**2 < ErrBasic)THEN
           !	Vect2=CrossProduct1(Tangent1, Vect1)/(T1Mag*Vect1Mag)
           Vect2=-Axis1
           Vect2Mag=SQRT(SUM(Vect2*vect2))											!dummy for conversion
           cosTh=SUM(Axis2*Vect1)/(Axis2Mag*Vect1Mag)
           sinTh=SUM(vect1*Vect2)/(vect2Mag*Vect1Mag)
         ELSE 
           WRITE(Message,'(A,A,I9)') 'Anisotropic Circular y-axis (radial) is not in-planar with element',&
               'for Body ',ThisBody
           CALL Warn('ShellSolve',Message)
           CALL Fatal('ShellSolve','STOP')
         ENDIF

       ELSEIF(AnisoType == 'spherical')THEN
         !	Point1-A is the center of sphere and point2-B is a point along line normal to Disk 
         !	whose periphery is first anisotropic direction
         Centroid(1) = SUM(Nodes%x(1:n))/n
         Centroid(2) = SUM(Nodes%y(1:n))/n
         Centroid(3) = SUM(Nodes%z(1:n))/n
         Vect1(1) = Centroid(1)-A(1)
         Vect1(2) = Centroid(2)-A(2)
         Vect1(3) = Centroid(3)-A(3)
         Vect1Mag = SQRT( SUM(Vect1**2))
         CP1 = CrossProduct1(Tangent1,Tangent2)/(T1Mag*T2Mag)
         Vect2 = CrossProduct1(Vect1,CP1)/Vect1Mag
         Err = SUM(Vect2*Vect2)
         IF(Err**2 < ErrBasic)THEN
           Vect2(1)=B(1)-A(1)
           Vect2(2)=B(2)-A(2)
           Vect2(3)=B(3)-A(3)
           Vect2Mag=SQRT( SUM(Vect2**2))
           Vect3=CrossProduct1(Vect1,Vect2)/(Vect1Mag*Vect2Mag)	!Vect3 is now unit vector in aniso x-direction
           Vect3Mag=SQRT(SUM(Vect3*Vect3))
           !	Vect2=CrossProduct1(Tangent1, CP2)/(T1Mag)	!Note CP2 is unit vector (Vect2 used as dummy).
           sinTh=SUM(Vect3*Axis2)/(Vect3Mag*Axis2Mag)
           cosTh=SUM(Axis1*Vect3)/(Axis1Mag*Vect3Mag)

           !	sinTh=SQRT(SUM(vect2*Vect2))
           !	cosTh=SUM(Tangent1*CP2)/(T1Mag)
         ELSE 
           WRITE(Message,'(A,A,A,I9)') 'x-axis is not in-planar with element for Anisotropy ',& 
               AnisoType, 'for Body ',ThisBody
           CALL Warn('ShellSolve',Message)
           CALL Fatal('ShellSolve','STOP')
         ENDIF

       ELSE IF(AnisoType == 'conical')THEN
         Centroid(1) = SUM(Nodes%x(1:n))/n
         Centroid(2) = SUM(Nodes%y(1:n))/n
         Centroid(3) = SUM(Nodes%z(1:n))/n
         Vect1(1) = Centroid(1)-A(1)
         Vect1(2) = Centroid(2)-A(2)
         Vect1(3) = Centroid(3)-A(3)
         Vect1Mag = SQRT( SUM(Vect1**2))

         !	Point1-A to Centroid of element is Anisotropic Y-axis. However, check it if it is
         !	in-planar to element (So as to ensure that Point1-A is apex of cone)
         CP1 = CrossProduct1(Vect1,Tangent1)/(Vect1Mag*T1Mag)
         CP2 = CrossProduct1(Vect1,Tangent2)/(Vect1Mag*T2Mag)
         Vect2 = CrossProduct1(CP1,CP2)
         Err = SUM(Vect2*Vect2)
         IF(Err**2 < ErrBasic)THEN
           !	Vect2=CrossProduct1(Tangent1, Vect1)/(T1Mag*Vect1Mag)
           
           Vect2=-Axis1
           Vect2Mag=SQRT(SUM(Vect2*vect2))											!dummy for conversion
           cosTh=SUM(Axis2*Vect1)/(Axis2Mag*Vect1Mag)
           sinTh=SUM(vect1*Vect2)/(vect2Mag*Vect1Mag)

           !	cosTh=SQRT(SUM(Vect2*Vect2))
           !	sinTh=SUM(Tangent1*Vect1)/(T1Mag*Vect1Mag)
         ELSE 
           WRITE(Message,'(A,A,A,I9)') 'x-axis is not in-planar with element for Anisotropy ',&
               AnisoType, 'for Body ',ThisBody
           CALL Warn('ShellSolve',Message)
           CALL Fatal('ShellSolve','STOP')
         ENDIF
       ELSE
         !	Give dummy values for cosine and sine as coincident Aniso and elemental local axes.
         cosTh=1.0d0
         sinTh=0.0d0
       ENDIF	!Loop on all anisotropies completed
       !	 WRITE(Message,'(A,2E12.4)') 'Values of cosTh and sinTh are  ',cosTh,sinTh
       !     CALL Info('ShellSolve',Message)

     END SUBROUTINE GetAnisotropyAxes

     
     ! (The subroutine is to be used within Int Pt Loop)
     ! ThLoad is initialized to be zero in the beginning of subroutine 'LocalMatrix'. Thus 
     ! when it is added to the previous value in this subroutine, the values from all int pts are added. 
     !	------------------------------------------------------------------------
     SUBROUTINE ThermalAndBendingLoad(Amatrix,Dmatrix,EPS,Kappa,Tmatrix,Basis,dBasisdx,s,alfa,&
         TT,ThLoad1,n,MaxDofs,IniBen2,BendingLoad1,RotMatrix2D)	!Changed 26.9.15	
       REAL(kind=dp) :: Amatrix(3,3),EPS(:,:),Tmatrix(3,3),Basis(n),dBasisdx(n,3),ThLoad1(6*n),TT(:),&
           s,DummyMatrix(MaxDofs,3),Dmatrix(:,:),Kappa(:,:),BendingLoad1(6*n),RotMatrix2D(3,3)
       REAL(kind=dp) :: alfa(2),Tepsilon_1(3),Tepsilon_2(3),Te,Tbending_1(3),Tbending_2(3),Tbending_3(3)
       REAL(kind=dp) ::	IniBen2(2) 
       INTEGER :: n,i,MaxDofs

       
       Te=SUM( TT(1:n) * Basis(1:n) )
       i=6*n
       IF( thermalanisotropy )THEN
         Tepsilon_1(1) = alfa(1)*Te
         Tepsilon_1(2) = alfa(2)*Te
         Tepsilon_1(3) = 0.0d0
         Tepsilon_2=MATMUL(Tmatrix,Tepsilon_1)	! Tmatrix changed to RotMatrix2D 
         !See remark in description of Function ThStrainConversion( )
       ELSE
         Tepsilon_2(1) = alfa(1)*Te
         Tepsilon_2(2) = alfa(1)*Te
         Tepsilon_2(3) = 0.0d0
       ENDIF
       
       ThLoad1=0.0d0

       DummyMatrix=MATMUL(TRANSPOSE(EPS),Amatrix)       
       ThLoad1(1:i) = Thload1(1:i)+MATMUL(DummyMatrix(1:i,1:3),Tepsilon_2)*s	

       !	Here thickness is already incorporated in Amatrix
       !	Considering Initial Bending strain  
       !	IniBen(2,2) contain the values of rotation per unit length i.e. 
       !	dbetax/dx, dBetax/dy and (2nd Row) dbetay/dx, dBetay/dy (All inputs)(Not used in this form now)
       !	********WARNING-IniBen( ) for routine 3D Thermal case is to be treated separately*************
       IF(thermalanisotropy)THEN
         Tbending_1(1)= IniBen2(1)
         Tbending_1(2) =IniBen2(2)
         Tbending_1(3) =0.0
         Tbending_2=MATMUL(Tmatrix,Tbending_1) !See remark in description of Function ThStrainConversion( )
       ELSE
         Tbending_2(1)=IniBen2(1)
         Tbending_2(2)=IniBen2(2)
         Tbending_2(3) =0.0
       ENDIF

       !	Following has been done so that the input IniBen2 is interpreted properly.
       Tbending_2=Tbending_2*OutNormal
       
       BendingLoad1=0.0d0
       CALL AddBendingLoad(BendingLoad1,Dmatrix,Kappa,Tbending_2,3,6*n,s)
       !	BendingLoad1=MATMUL(MATMUL(Transpose(Kappa),Dmatrix),Tbending_2)*s

     END SUBROUTINE ThermalAndBendingLoad

!    ==================================================================
!	This CrossProduct1 gives the + direction of result vector v3 when we move from v1 to v2 (19.11.13)
! -----------------------------------------------------------------
     FUNCTION CrossProduct1( v1, v2 ) RESULT( v3 )	
! ----------------------------------------------------------------
       REAL(KIND=dp) :: v1(3), v2(3), v3(3)
       v3(1) =  v1(2)*v2(3) - v1(3)*v2(2)
       v3(2) = -v1(1)*v2(3) + v1(3)*v2(1)
       v3(3) =  v1(1)*v2(2) - v1(2)*v2(1)
! --------------------------------------------------------------------
     END FUNCTION CrossProduct1
! -------------------------------------------------------------------

     
     ! -------------------------------------------------------------------
     !	This function (Matrix[d]) is used for converting strain {E2} to {E1} as {E2}=[d]{E1}
     !	Equivalent conversion is {E1}=[dInv]{E2}. The matrix [d] here is actually [dInv] because
     !	CosTh and SinTh (inputs c and s) are determined for conversion from elemental local axes to 
     !	Anisotropic axes (and not from Anisotropic axes to elemental local axes, which is desired) 
     !	due to constraint of knowing only one Anisotropic axis in subroutine 'ThermalAndBendingLoad'.
     ! -------------------------------------------------------------------
     FUNCTION ThStrainConversion(c,s)RESULT(d)
       REAL(kind=dp) :: c,s,d(3,3)
       d(1,1)=c*c
       d(1,2)=s*s
       d(1,3)=-s*c
       d(2,1)=d(1,2)
       d(2,2)=d(1,1)
       d(2,3)=s*c
       d(3,1)=2.0*s*c
       d(3,2)=-2.0*s*c
       d(3,3)=c*c-s*s
     END FUNCTION ThStrainConversion

     
     SUBROUTINE getDKTtriangleParameter(Nodes,n,a,b,c,d,e,p,q,t,r)
       REAL (Kind=dp) :: a(6), b(6), c(6), d(6), e(6), p(6), q(6), r(6), t(6)
       REAL (Kind=dp) :: xx(3),yy(3),zz(3)
       INTEGER :: n,i,j,m,k
       TYPE(Nodes_t) :: Nodes
       DO i=1,3
         xx(i)=Nodes%x(i)
         yy(i)=Nodes%y(i)
         zz(i)=Nodes%z(i)
       ENDDO
       DO k=1,3
         SELECT CASE(k)
         CASE(1)
           i=1; j=2; m=4
           CALL getDKTvalues(i,j,m,xx,yy,zz,a,b,c,d,e,p,q,r,t)
         CASE(2)
           i=2; j=3; m=5
           CALL getDKTvalues(i,j,m,xx,yy,zz,a,b,c,d,e,p,q,r,t)
         CASE(3)
           i=3; j=1; m=6
           CALL getDKTvalues(i,j,m,xx,yy,zz,a,b,c,d,e,p,q,r,t)
         END SELECT
       ENDDO
     END SUBROUTINE getDKTtriangleParameter

     !	----------------------------------------------------------------------
     SUBROUTINE getDKTvalues(i,j,m,xx,yy,zz,a,b,c,d,e,p,q,r,t)
       REAL (Kind=dp) :: a(6), b(6), c(6), d(6), e(6), p(6), q(6), r(6), t(6)
       REAL(kind=dp) :: x(3,3),y(3,3),L(3,3),Lijsq
       REAL (Kind=dp) :: xx(3),yy(3),zz(3)
       INTEGER :: i,j,k,m 

       x(i,j)=xx(j)-xx(i)
       y(i,j)=yy(j)-yy(i)
       Lijsq=x(i,j)*x(i,j)+y(i,j)*y(i,j)
       a(m)=-x(i,j)/Lijsq; b(m)=0.75*x(i,j)*y(i,j)/Lijsq; c(m)=(0.25*x(i,j)*x(i,j)-0.5*y(i,j)*y(i,j))/Lijsq
       d(m)=-y(i,j)/Lijsq; e(m)=(0.25*y(i,j)*y(i,j)-0.5*x(i,j)*x(i,j))/Lijsq; r(m)=3.0*y(i,j)*y(i,j)/Lijsq
       p(m)=6.0*a(m); q(m)=4.0*b(m); t(m)=6.0*d(m)
     END SUBROUTINE getDKTvalues
     !	-----------------------------------------------------------------------------                                                                                                                                         
     SUBROUTINE DKTformulation(Nodes, Kappa,MaxDofs,n, u,v,w,s,a,b,c,d,e,p,q,t,r) 
       !u,v,w,s are the zy,eta,zeta coordinates and area weightage of int pt.
       REAL (Kind=dp) :: a(6), b(6), c(6), d(6), e(6), p(6), q(6), r(6), t(6)
       REAL(kind=dp) :: u,v,w,s,x12,y12,x13,y13,fzy,feta,xx(3),yy(3),zz(3),Area2
       REAL(kind=dp) :: HxzyT(9), HxetaT(9), HyzyT(9), HyetaT(9)
       TYPE(Nodes_t) :: Nodes
       REAL(Kind=dp) :: Kappa(3,MaxDofs), Gama(3,9)
       INTEGER :: MaxDofs,n,i,j,k

       fzy=1.0-2.0*u
       feta=1.0-2.0*v
       HxzyT(1)=p(4)*fzy+(p(6)-p(4))*v
       HxzyT(2)=-q(4)*fzy+(q(4)+q(6))*v
       HxzyT(3)=4.0-6.0*(u+v)-r(4)*fzy+(r(4)+r(6))*v
       HxzyT(4)=-p(4)*fzy+(p(4)+p(5))*v
       HxzyT(5)=-q(4)*fzy+(q(4)-q(5))*v
       HxzyT(6)=2.0-6.0*u-r(4)*fzy-(r(5)-r(4))*v
       HxzyT(7)=-(p(5)+p(6))*v
       HxzyT(8)=-(q(5)-q(6))*v
       HxzyT(9)=(r(6)-r(5))*v
       HxetaT(1)=-p(6)*feta-(p(4)-p(6))*u
       HxetaT(2)=-q(6)*feta+(q(4)+q(6))*u
       HxetaT(3)=4.0-6.0*(u+v)-r(6)*feta+(r(4)+r(6))*u
       HxetaT(4)=(p(4)+p(5))*u
       HxetaT(5)=-(q(5)-q(4))*u
       HxetaT(6)=(r(4)-r(5))*u
       HxetaT(7)=p(6)*feta-(p(5)+p(6))*u
       HxetaT(8)=-q(6)*feta-(q(5)-q(6))*u
       HxetaT(9)=2.0-6.0*v-(r(5)-r(6))*u-r(6)*feta
       HyzyT(1)=t(4)*fzy+(t(6)-t(4))*v
       HyzyT(2)=-1.0-r(4)*fzy+(r(4)+r(6))*v
       HyzyT(3)=q(4)*fzy-(q(4)+q(6))*v
       HyzyT(4)=-t(4)*fzy+(t(4)+t(5))*v
       HyzyT(5)=1.0-r(4)*fzy-(r(5)-r(4))*v
       HyzyT(6)=q(4)*fzy+(q(5)-q(4))*v
       HyzyT(7)=-(t(6)+t(5))*v
       HyzyT(8)=-(r(5)-r(6))*v
       HyzyT(9)=(q(5)-q(6))*v
       HyetaT(1)=-t(6)*feta-(t(4)-t(6))*u
       HyetaT(2)=-1.0+(r(4)+r(6))*u-r(6)*feta
       HyetaT(3)=-(q(4)+q(6))*u+q(6)*feta
       HyetaT(4)=(t(4)+t(5))*u
       HyetaT(5)=-(r(5)-r(4))*u
       HyetaT(6)=(q(5)-q(4))*u
       HyetaT(7)=t(6)*feta-(t(5)+t(6))*u
       HyetaT(8)=1.0-r(6)*feta-(r(5)-r(6))*u
       HyetaT(9)=(q(5)-q(6))*u+q(6)*feta

       DO i=1,3
         xx(i)=Nodes%x(i)
         yy(i)=Nodes%y(i)
         zz(i)=Nodes%z(i)
       END DO

       x12=xx(2)-xx(1)
       y12=yy(2)-yy(1)
       x13=xx(3)-xx(1)
       y13=yy(3)-yy(1)

       Area2=(xx(2)-xx(1))*yy(3)+(xx(1)-xx(3))*yy(2)+(xx(3)-xx(2))*yy(1) !Area2 removed temp 18.8.14
       Gama(1,1:9)=(y13*HxzyT-y12*HxetaT)/Area2
       Gama(2,1:9)=(-x13*HyzyT+x12*HyetaT)/Area2
       Gama(3,1:9)=(-x13*HxzyT+x12*HxetaT+y13*HyzyT-y12*HyetaT)/Area2
       
       Kappa=0.0
       DO i=1,n
         j=(i-1)*3
         k=(i-1)*6
         Kappa(1:3,k+3)=Gama(1:3,j+1)
         Kappa(1:3,k+4)=Gama(1:3,j+2)
         Kappa(1:3,k+5)=Gama(1:3,j+3)
       ENDDO
     END SUBROUTINE DKTformulation

!	----------------------------------------------------------------------------------------
     SUBROUTINE RDKTelement(Nodes, Kappa,MaxDofs,n, u,v,w,s,a,d)
       !u,v,w,s are the zy,eta,zeta coordinates and area weightage of int pt.
       REAL (Kind=dp) :: a(6),  d(6)
       REAL(kind=dp) :: u,v,w,s,x12,y12,x13,y13,fzy,feta,xx(3),yy(3),zz(3),Area2
       REAL(kind=dp) :: HxzyT(9), HxetaT(9), HyzyT(9), HyetaT(9)
       TYPE(Nodes_t) :: Nodes
       REAL(Kind=dp) :: kappa(3,MaxDofs), Gama(3,9)
       INTEGER :: MaxDofs,n,i,j,k

       HxzyT(1)=6.0*a(4)*(1.0-2.0*u-v)+6.0*a(6)*v
       HxzyT(2)=0.0
       HxzyT(3)=4.0-6.0*u-6.0*v
       HxzyT(4)=-6.0*a(4)*(1.0-2.0*u-v)+6.0*a(5)*v
       HxzyT(5)=0.0
       HxzyT(6)=2.0-6.0*u
       HxzyT(7)=-6.0*v*(a(5)+a(6))
       HxzyT(8)=0.0
       HxzyT(9)=0.0
       HxetaT(1)=-6.0*a(4)*u-6.0*a(6)*(1.0-u-2.0*v)
       HxetaT(2)=0.0
       HxetaT(3)=4.0-6.0*u-6.0*v
       HxetaT(4)=6.0*u*(a(4)+a(5))
       HxetaT(5)=0.0
       HxetaT(6)=0.0
       HxetaT(7)=-6.0*a(5)*u+6.0*a(6)*(1.0-u-2.0*v)
       HxetaT(8)=0.0
       HxetaT(9)=2.0-6.0*v
       HyzyT(1)=6.0*d(4)*(1.0-2.0*u-v)+6.0*d(6)*v
       HyzyT(2)=-(4.0-6.0*u-6.0*v)
       HyzyT(3)=0.0
       HyzyT(4)=-6.0*d(4)*(1.0-2.0*u-v)+6.0*d(5)*v
       HyzyT(5)=-(2.0-6.0*u)
       HyzyT(6)=0.0
       HyzyT(7)=-6.0*v*(d(5)+d(6))
       hyzyT(8)=0.0
       HyzyT(9)=0.0
       HyetaT(1)=-6.0*d(4)*u-6.0*d(6)*(1.0-u-2.0*v)
       HyetaT(2)=-(4.0-6.0*u-6.0*v)
       HyetaT(3)=0.0
       HyetaT(4)=6.0*u*(d(4)+d(5))
       HyetaT(5)=0.0
       HyetaT(6)=0.0
       HyetaT(7)=-6.0*d(5)*u+6.0*d(6)*(1.0-u-2.0*v)
       HyetaT(8)=-(2.0-6.0*v)
       HyetaT(9)=0.0

       DO i=1,3
         xx(i)=Nodes%x(i)
         yy(i)=Nodes%y(i)
         zz(i)=Nodes%z(i)
       ENDDO

       x12=xx(2)-xx(1)
       y12=yy(2)-yy(1)
       x13=xx(3)-xx(1)
       y13=yy(3)-yy(1)
       Area2=(xx(2)-xx(1))*yy(3)+(xx(1)-xx(3))*yy(2)+(xx(3)-xx(2))*yy(1)
       gama(1,1:9)=(y13*HxzyT-y12*HxetaT)/Area2
       gama(2,1:9)=(-x13*HyzyT+x12*HyetaT)/Area2
       gama(3,1:9)=(-x13*HxzyT+x12*HxetaT+y13*HyzyT-y12*HyetaT)/Area2
       kappa=0.0

       DO i=1,n
         j=(i-1)*3
         k=(i-1)*6
         kappa(1:3,k+3)=gama(1:3,j+1)
         kappa(1:3,k+4)=gama(1:3,j+2)
         kappa(1:3,k+5)=gama(1:3,j+3)
       ENDDO
       
     END SUBROUTINE RDKTelement
!	-----------------------------------------------------------------


     SUBROUTINE CovariantInterpolationR(ShearStrain,Basis,X,Y,U,V,n)
! -------------------------------------------------------------------
       REAL(KIND=dp) :: ShearStrain(:,:),Basis(:),X(:),Y(:),U,V
       INTEGER :: n
! ------------------------------------------------------------------
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
!		NOTE:- Changes made here in multiplying factor from 1+V to -1+V and -1-U t0 1-U 
!				Have been thoroughly verified from formulation. There is an additional 
!				error in Tblock which combined with this leads to some sort of correction 
!				(Not total correction).

!         Compute the shear-dofs for edge 12:
!         ===================================
         Tau = (/ 1.0d0, 0.0d0/)
         
         Sdofs = 0.0d0
         Sdofs(4) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/2.0d0
         Sdofs(5) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/2.0d0
         Sdofs(10) = (Tau(1)*Jmat(1,1)+Tau(2)*Jmat(2,1))/2.0d0
         Sdofs(11) = (Tau(1)*Jmat(1,2)+Tau(2)*Jmat(2,2))/2.0d0
         
         DO j = 1,18
           ShearRef(1,j) = ShearRef(1,j) + (-1+V)*Sdofs(j)		!1+V changed -1+V 20.9.14
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
           ShearRef(2,j) = ShearRef(2,j) + (1-U)*Sdofs(j)		!-1-U changed to 1-U 20.9.14
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
! ------------------------------------------------------------------
     END SUBROUTINE CovariantInterpolationR
! ----------------------------------------------------------------

     SUBROUTINE CovariantInterpolationMITC3(ShearStrain, Basis, &
         X,Y,U, V, n)
! -------------------------------------------------------------------
       REAL(KIND=dp) :: ShearStrain(:,:),Basis(:),X(:),Y(:),U,V
       INTEGER :: n
! ------------------------------------------------------------------
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
         
!		Determining -a1
         Sdofs = 0.0d0
         Sdofs(4) = -(x(2)-x(1))/2.0d0
         Sdofs(5) = -(y(2)-y(1))/2.0d0
         Sdofs(10) = -(x(2)-x(1))/2.0d0
         Sdofs(11) =-(y(2)-y(1))/2.0d0
         
         DO j = 1,18
           ShearRef(1,j) = ShearRef(1,j) + Sdofs(j)	
           !             ShearRef(2,j) = ShearRef(2,j) + ( -U)*Sdofs(j)
         END DO
!		Determining -a2
         Sdofs = 0.0d0
         Sdofs(4) = -(x(3)-x(1))/2.0d0
         Sdofs(5) = -(y(3)-y(1))/2.0d0
         Sdofs(16) = -(x(3)-x(1))/2.0d0
         Sdofs(17) = -(y(3)-y(1))/2.0d0
         
         DO j = 1,18
           !             ShearRef(1,j) = ShearRef(1,j) + (  V )*Sdofs(j)
           ShearRef(2,j) = ShearRef(2,j) + Sdofs(j)		
         END DO
         
         !		Determining c/2 (not c0/2) and substituting it
         Sdofs = 0.0d0
         Sdofs(4) = -(x(3)-x(2))/2.0d0
         Sdofs(5) = -(y(3)-y(2))/2.0d0
         Sdofs(10) = -(x(1)-x(3))/2.0d0
         Sdofs(11) = -(y(1)-y(3))/2.0d0
         Sdofs(16) = -(x(2)-x(1))/2.0d0
         Sdofs(17) = -(y(2)-y(1))/2.0d0
         
         DO j = 1,18
           ShearRef(1,j) = ShearRef(1,j) + (  V )*Sdofs(j)
           ShearRef(2,j) = ShearRef(2,j) -U*Sdofs(j)		
         END DO
         !		Determining dw/dzy 
         Sdofs = 0.0d0
         Sdofs(3) = -1.0
         Sdofs(9) = 1.0
         
         DO j = 1,18
           ShearRef(1,j) = ShearRef(1,j) + Sdofs(j)
         END DO
         
         !		Determining dw/deta 
         Sdofs = 0.0d0
         Sdofs(3) = -1.0
         Sdofs(15) = 1.0
         
         DO j = 1,18
           ShearRef(2,j) = ShearRef(2,j) + Sdofs(j)
         END DO
         !         Compute the final reduced shear strain
         !         ======================================
         ShearStrain(1:2,1:18) = MATMUL(invJ,ShearRef(1:2,1:18))
         
       CASE DEFAULT
         CALL Fatal('ShellSolver','Illegal number of nodes for MITC3 elements.')
       END SELECT
       
     END SUBROUTINE CovariantInterpolationMITC3


     !	------------------------------------------------------------------------------
     SUBROUTINE MembraneFormulation(MemMatrix,n)
!	Calculate stiffness matrix for membrane onle i.e. enforce Thetax, Thetay, Thetaz terms zero
       INTEGER j,k,n
       REAL(Kind=dp) :: MemMatrix(6*n,6*n)
       MemMatrix=0.0d0
       DO k=1,n
         j=6*(k-1)+4
         MemMatrix(j,j)=1.0
         MemMatrix(j+1,j+1)=1.0
         MemMatrix(j+1,j+2)=1.0
       ENDDO
     END SUBROUTINE MembraneFormulation


     SUBROUTINE RotationMinimizationDKT(Rotation,Basis,X,Y,U,V,n)
! -------------------------------------------------------------------
       REAL(KIND=dp) :: Rotation(:,:),Basis(:),X(:),Y(:),U,V
       INTEGER :: n
! ------------------------------------------------------------------
       REAL(KIND=dp) :: detJ,Jmat(2,2),invJ(2,2),ShearRef(2,100)
       INTEGER :: j
       SELECT CASE(n)
         
!      The DKT element
!      ==================
       CASE(3)
!          CALL Jacobi3(Jmat,invJ,detJ,x,y)
         !          ShearRef = 0.0d0
         Rotation = 0.0d0
         DO j=1,n
           Rotation(1,6*j-2)=0.5*Basis(j)
           Rotation(2,6*j-1)=-0.5*Basis(j)
         ENDDO
         
         !         Compute the final reduced shear strain
         !         ======================================
!          ShearStrain(1:2,1:18) = MATMUL(invJ,ShearRef(1:2,1:18))
       CASE DEFAULT
         CALL Fatal('ShellSolver','Illegal number of nodes for 3 node elements.')
       END SELECT
       ! ------------------------------------------------------------------
     END SUBROUTINE RotationMinimizationDKT
! ----------------------------------------------------------------


     SUBROUTINE MembraneFormulationNRM(MemMatrix,n)
       !	Calculate stiffness matrix for membrane onle i.e. enforce Thetax, Thetay, Thetaz terms zero
       INTEGER j,k,n
       REAL(Kind=dp) :: MemMatrix(6*n,6*n)

       MemMatrix=0.0d0
       DO k=1,n
         j=6*(k-1)+4
         MemMatrix(j-1,j-1)=1.0
         MemMatrix(j,j)=1.0
         MemMatrix(j+1,j+1)=1.0
         MemMatrix(j+1,j+2)=1.0
       ENDDO
     END SUBROUTINE MembraneFormulationNRM


     
     SUBROUTINE CovariantInterpolationRMITC3(ShearStrain, Basis, &	!29.9.14
         X,Y,U, V, n)
!	It was observed that the value of parameter 'c' is redundant (see notes) and so the 
!	terms containing this were removed and results are more accurate now.
! -------------------------------------------------------------------
       REAL(KIND=dp) :: ShearStrain(:,:),Basis(:),X(:),Y(:),U,V
       INTEGER :: n
       ! ------------------------------------------------------------------
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
         
         !		Determining -a1
         Sdofs = 0.0d0
         Sdofs(4) = -(x(2)-x(1))/2.0d0
         Sdofs(5) = -(y(2)-y(1))/2.0d0
         Sdofs(10) = -(x(2)-x(1))/2.0d0
         Sdofs(11) =-(y(2)-y(1))/2.0d0
         !          Sdofs(16) = -(x(2)-x(1))/2.0d0
         !         Sdofs(17) =-(y(2)-y(1))/2.0d0
         
         DO j = 1,18
           ShearRef(1,j) = ShearRef(1,j) + Sdofs(j)	
           !             ShearRef(2,j) = ShearRef(2,j) + ( -U)*Sdofs(j)
         END DO
         !		Determining -a2
         Sdofs = 0.0d0
         Sdofs(4) = -(x(3)-x(1))/2.0d0
         Sdofs(5) = -(y(3)-y(1))/2.0d0
         !          Sdofs(10) = -(x(3)-x(1))/2.0d0
         !         Sdofs(11) = -(y(3)-y(1))/2.0d0
         Sdofs(16) = -(x(3)-x(1))/2.0d0
         Sdofs(17) = -(y(3)-y(1))/2.0d0
         
         DO j = 1,18
           !             ShearRef(1,j) = ShearRef(1,j) + (  V )*Sdofs(j)
           ShearRef(2,j) = ShearRef(2,j) + Sdofs(j)		
         END DO
         
         
!	*************Calculations for dw/dzy etc removed 30.9.15********************

         !         Compute the final reduced shear strain
         !         ======================================
         ShearStrain(1:2,1:18) = MATMUL(invJ,ShearRef(1:2,1:18))
         
       CASE DEFAULT
         CALL Fatal('ShellSolver','Illegal number of nodes for RMITC3 elements.')
       END SELECT
       
     END SUBROUTINE CovariantInterpolationRMITC3


!------------------------------------------------------------------------------
     SUBROUTINE AddBendingLoad(A,B,C,D,m,n,s)
! -----------------------------------------------------------------
!      Performs the operation
!
!         A = A + C' * B * D * s
!
!      with
!
!         Size( A ) = n x 1
!         Size( B ) = m x m
!         Size( C ) = m x n
!		   Size(D)  = m x 1
!  --------------------------------------------------------------
       INTEGER :: m,n
       
       REAL(KIND=dp) :: A(:),B(:,:),C(:,:),D(:),s,E(m)
! -----------------------------------------------------------------
       INTEGER :: i,j,k,l
! ---------------------------------------------------------------
 !      DO i=1,n
 !         DO j=1,m
 !            DO k=1,m
 !               DO l=1,m
 !                  A(i) = A(i) + C(j,k)*B(k,l)*D(l) * s ! to check
 !               END DO
 !            END DO
 !         END DO
 !      END DO
       E=0.0
       DO k=1,m
         DO l=1,m
           E(k) = E(k) + B(k,l)*D(l) * s 
         END DO
       END DO
       DO i=1,n
         DO j=1,m
           A(i)=A(i)+C(j,i)*E(j)
         ENDDO
       ENDDO
       
! -------------------------------------------------------------------
     END SUBROUTINE AddBendingLoad


! ---------------------------------------------------------------
   END SUBROUTINE ShellSolver
! ---------------------------------------------------------------
