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
! *  Module for solving the two-dimensional Reissner-Naghdi shell equations using
! *  elementwise lines of curvature coordinates. The way to obtain a lines of 
! *  curvature parametrization was first announced in the paper 
! *  
! *  [1] Malinen M, Generating lines of curvature coordinates for finite element 
! *  modelling, Proceedings of the XII Finnish Mechanics Days, 2015.
! *
! *  The improved approximation of the shell mid-surface is derived from the surface 
! *  director (normal vector) data as outlined in
! *  
! *  [2] Malinen M, Improved surface reconstruction from conventional geometry
! *  data for general shell finite elements, Proceedings of the 29th Nordic
! *  Seminar on Computational Mechanics, 2016.
! *
! *  The nodal director data should be available via an ordinary solver variable
! *  'Director' or via reading from file mesh.director located in the same place 
! *  as the standard mesh files or, as the third option, the user may provide 
! *  mesh.elements.data file which should define the nodal director field associated 
! *  with the name 'director'.
! *
! *  This solver is STILL UNDER DEVELOPMENT and some possibilities of the strategy
! *  are not yet fully utilized. Note the current restrictions:
! *        -- Strain reduction operators have been worked out for 
! *           the lowest-order finite elements only.
! *        -- Only a plain p-element discretization without applying strain reduction
! *           operators is possible 
! *        -- Parallel file formats for mesh.director and mesh.elements.data are missing,
! *           so for parallel execution the director should be defined as an ordinary 
! *           solver variable
! *        -- Postprocessing routines are also missing 
! *        -- Terms of O(d/R), with d the shell thickness and R the minimum of
! *           radius of curvature, are ignored in the expression for the strain energy 
! *
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: Jan 22, 2015
! *
! *****************************************************************************/



!------------------------------------------------------------------------------
SUBROUTINE ShellSolver_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverPars, Simulation
  LOGICAL :: SavePrincipalAxes, Found, RotateDOFs, Eigenanalysis
  INTEGER  :: i
!------------------------------------------------------------------------------
  SolverPars => GetSolverParams()
  Simulation => GetSimulation()

  CALL ListAddNewInteger(SolverPars, 'Variable DOFs', 6)
  CALL ListAddLogical(SolverPars, 'Bubbles in Global System', .TRUE.)
  CALL ListAddLogical(Simulation, 'Initialize Dirichlet Conditions', .FALSE.)

  RotateDOFs = GetLogical(SolverPars, 'Rotate DOFs', Found)
  IF (RotateDOFs) THEN
    CALL ListAddNewString(SolverPars, 'Variable', 'Deflection[U:3 Theta:3]')
    CALL ListAddNewString(SolverPars, 'Imaginary Variable', 'Deflection[U im:3 Theta im:3]')
  ELSE
    CALL ListAddNewString(SolverPars, 'Variable', 'Deflection[U:3 DNU:3]')
    ! Only created if the system is harmonic
    CALL ListAddNewString(SolverPars, 'Imaginary Variable', 'Deflection[U im:3 DNU im:3]')
  END IF

  Eigenanalysis = GetLogical(SolverPars, 'Eigen Analysis', Found)
  IF (Eigenanalysis) THEN
    CALL ListAddLogical(SolverPars, 'Large Deflection', .FALSE.)
  ELSE
    CALL ListAddNewLogical(SolverPars, 'Large Deflection', .TRUE.)
    CALL ListAddNewInteger(SolverPars, 'Nonlinear System Max Iterations', 50)
    CALL ListAddNewConstReal(SolverPars, 'Nonlinear System Convergence Tolerance', 1.0d-5)
    IF (Transient) THEN
      CALL ListAddInteger(SolverPars, 'Time derivative order', 2)
      CALL ListAddString(SolverPars, 'Timestepping Method', 'Bossak')
    END IF
  END IF
  CALL ListAddNewLogical(SolverPars, 'Skip Compute Nonlinear Change', .TRUE.)

  !----------------------------------------------------------------------------
  ! Create variables for saving principal (curvature) directions:
  !----------------------------------------------------------------------------
  SavePrincipalAxes = GetLogical(SolverPars, 'Principal Axes Output', Found)

  IF (SavePrincipalAxes) THEN
    i=1
    DO WHILE(.TRUE.)
      IF ( .NOT.ListCheckPresent(SolverPars, &
          "Exported Variable "//i2s(i)) ) EXIT
      i = i + 1
    END DO
    CALL ListAddString(SolverPars, "Exported Variable "//i2s(i), &
        "Principal Coordinate Dir1[Principal Coordinate Dir1:3]")  
    i = i + 1
    CALL ListAddString(SolverPars, "Exported Variable "//i2s(i), &
        "Principal Coordinate Dir2[Principal Coordinate Dir2:3]")
    i = i + 1
    CALL ListAddString(SolverPars, "Exported Variable "//i2s(i), &
        "Principal Coordinate Dir3[Principal Coordinate Dir3:3]")  
  END IF

  CALL ListAddLogical( SolverPars,'Shell Solver',.TRUE.)

  IF( GetLogical( SolverPars, 'Stability Analysis', Found ) ) THEN
    CALL Fatal('ShellSolver_Init0','"Stability Analysis" has not yet been coded for this solver!')
  END IF
  IF( GetLogical( SolverPars, 'Geometric Stiffness', Found ) ) THEN
    CALL Fatal('ShellSolver_Init0','"Geometric Stiffness" has not yet been coded for this solver!')
  END IF
      
!------------------------------------------------------------------------------
END SUBROUTINE ShellSolver_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE ShellSolver(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  USE ElementDescription
  USE SolidMechanicsUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Variables which are intended to have general visibility and relate to saving
! some element properties as data chunks of the same size.
!------------------------------------------------------------------------------
! For coordinate frame data:
!------------------------------------------------------------------------------
  INTEGER, PARAMETER :: FrameDataSize = 12
  INTEGER, DIMENSION(3), PARAMETER :: FrameBasis1 = (/ 1, 2, 3 /)
  INTEGER, DIMENSION(3), PARAMETER :: FrameBasis2 = (/ 4, 5, 6 /)
  INTEGER, DIMENSION(3), PARAMETER :: FrameBasis3 = (/ 7, 8, 9 /)
  INTEGER, DIMENSION(3), PARAMETER :: FrameOrigin = (/ 10, 11, 12 /)
!------------------------------------------------------------------------------
! For edge parametrizations based on the director data of first and second order:
!------------------------------------------------------------------------------
  INTEGER, PARAMETER :: MaxNumberOfCurves = 6
  INTEGER, PARAMETER :: CurveDataSize1 = 6
  INTEGER, PARAMETER :: CurveDataSize2 = 6
  INTEGER, DIMENSION(CurveDataSize1), PARAMETER :: CurveParams1 = (/ 1, 2, 3, 4, 5, 6/)
  INTEGER, DIMENSION(CurveDataSize2), PARAMETER :: CurveParams2 = (/ 1, 2, 3, 4, 5, 6/) 
!------------------------------------------------------------------------------
! Some other parameter definitions related to the choice of strain reduction 
! method, geometric tolerances ... 
!------------------------------------------------------------------------------
  INTEGER, PARAMETER :: AutomatedChoice = -1
  INTEGER, PARAMETER :: NoStrainReduction = 0
  INTEGER, PARAMETER :: CurlKernel = 1             ! This builds on Ker(curl) of either RT_0 or ABF_0   
  INTEGER, PARAMETER :: MITC = 2                   ! This builds on RT_0
  INTEGER, PARAMETER :: DoubleReduction = 3        ! This builds on Ker(curl) of RT_0 (triangles)
  INTEGER, PARAMETER :: CurlKernelWithEdgeDOFs = 4 ! This also builds on Ker(curl) of ABF_0
  INTEGER, PARAMETER :: ExperimentalReduction = 5  ! A developer's option

  INTEGER, PARAMETER :: MaxBGElementNodes = 9
  INTEGER, PARAMETER :: MaxPatchNodes = 16    ! The maximum node count for the surface description 

  INTEGER, PARAMETER :: GeometryMaxIters = 50
  REAL(KIND=dp), PARAMETER :: ParabolicDelta = 1.0d-8  ! The tolerance to decide parabolic points
  REAL(KIND=dp), PARAMETER :: GeometryEpsilon = 1.0d-6
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t), POINTER :: Element, BGElement, Parent
  TYPE(ElementType_t), POINTER :: ShellElement => NULL()
  TYPE(Nodes_t) :: Nodes
  TYPE(ValueList_t), POINTER :: SolverPars
  TYPE(Variable_t), POINTER :: NodalPDir1, NodalPDir2, NodalPDir3, Director
  TYPE(Variable_t), POINTER :: Displacement3D
  TYPE(Matrix_t), POINTER :: PMatrix

  LOGICAL :: Found
  LOGICAL :: CurveDataOutput, SavePrincipalAxes, ComputeShellArea
  LOGICAL :: WriteElementalDirector
  LOGICAL :: QuadraticApproximation = .FALSE.
  LOGICAL :: PlateBody, PlanarPoint, UmbilicalPoint
  LOGICAL :: Bubbles, ApplyBubbles
  LOGICAL :: LargeDeflection, MeshDisplacementActive, Relax
  LOGICAL :: NoTractions
  LOGICAL :: SolveBenchmarkCase
  LOGICAL :: MassAssembly, HarmonicAssembly
  LOGICAL :: Parallel
  LOGICAL :: SolidShellCoupling
  LOGICAL :: DrillingDOFs, RotateDOFs
  LOGICAL :: CartesianFormulation, SkipBlending, ReparametrizeMesh, PVersion
  LOGICAL :: NonlinearBending

  INTEGER, POINTER :: Indices(:) => NULL()
  INTEGER, POINTER :: VisitsList(:) => NULL()
  INTEGER :: e, i, i0, j, k, m, n, nb, nd, nd_parent, t 
  INTEGER :: Family, Active
  INTEGER :: ShellModelPar, StrainReductionMethod, MembraneStrainReductionMethod
  INTEGER :: NonlinIter, MaxNonlinIters
  
  REAL(KIND=dp), POINTER :: TotalSol(:) => NULL()
  REAL(KIND=dp), POINTER CONTIG :: ValuesSaved(:) => NULL()
  REAL(KIND=dp), POINTER :: TaylorParams(:)
  REAL(KIND=dp), POINTER :: Pb(:)
  REAL(KIND=dp), ALLOCATABLE :: LocalSol(:,:)
  REAL(KIND=dp), ALLOCATABLE :: LocalRHSForce(:)
  REAL(KIND=dp), TARGET :: LocalFrameNodes(MaxPatchNodes,3)
  REAL(KIND=dp) :: d(3), d1(3), d2(3), d3(3), X1(3), X2(3), e1(3), e2(3), e3(3), o(3), p(3)
  REAL(KIND=dp) :: c, Norm, u, v
  REAL(KIND=dp) :: PatchNodes(MaxPatchNodes,2), ZNodes(MaxPatchNodes)
  REAL(KIND=dp) :: BlendingSurfaceArea, ShellModelArea, MappedMeshArea, RefArea
  REAL(KIND=dp) :: NonlinTol, NonlinRes, NonlinRes0, Relaxation
  REAL(KIND=dp) :: DrillingPar
  REAL(KIND=dp) :: ShearAlpha, MembraneAlpha, StretchAlpha

  CHARACTER(LEN=MAX_NAME_LEN) :: OutputFile, ElementDef

  ! Variables for development version: 
  REAL(KIND=dp) :: TotalErr
  REAL(KIND=dp) :: RefWork, Work
  REAL(KIND=dp) :: MaxPDir1Err, MaxPDir2Err, PDir1(3), PDir2(3)
  REAL(KIND=dp) :: Energy(4), MEnergy, SEnergy, BEnergy, Etot
  
  SAVE VisitsList, Indices, LocalSol, TotalSol, LocalRHSForce
!------------------------------------------------------------------------------  

  CALL DefaultStart()
  
  ! ---------------------------------------------------------------------------------
  ! PART 0:
  ! Obtain the values of some key parameters and create allocatable variables. 
  ! ---------------------------------------------------------------------------------
  Mesh => GetMesh()
  SolverPars => GetSolverParams()
  Parallel = ParEnv % PEs > 1

  PVersion = .FALSE.
  CartesianFormulation = GetLogical(SolverPars, 'Cartesian Formulation', Found)
  IF (.NOT. Found) THEN
    !
    ! Having a p-element definition switches to the Cartesian components formulation
    !
    ElementDef = ListGetString(SolverPars, 'Element', Found)
    IF (Found) THEN
      i = INDEX(ElementDef,'p:')
      IF (i > 0) THEN
        CartesianFormulation = .TRUE.
        PVersion = .TRUE.
      END IF
    END IF
  END IF
  SkipBlending = GetLogical(SolverPars, 'Skip Surface Reconstruction', Found)

  ReparametrizeMesh = GetLogical(SolverPars, 'Mesh Reparametrization', Found)
  IF (ReparametrizeMesh) THEN
    SkipBlending = .TRUE.
    CartesianFormulation = .TRUE.
    IF (PVersion) CALL Fatal('ShellSolver', 'Mesh Reparametrization incompatible with p-elements')
  END IF

  IF (CartesianFormulation) THEN
    CALL Info('ShellSolver', 'APPLYING CARTESIAN COMPONENTS FORMULATION', Level=3)
    CALL Info('ShellSolver', 'USE HIGH-ORDER BASIS FUNCTIONS TO HANDLE LOCKING', Level=3)
    ! ---------------------------------------------------------------------------------
    ! The number of unknown fields in the shell model:
    ! ---------------------------------------------------------------------------------
    ShellModelPar = ListGetInteger(SolverPars, 'Variable DOFs', minv=6, maxv=9)
  ELSE
    ShellModelPar = ListGetInteger(SolverPars, 'Variable DOFs', minv=6, maxv=6)
  END IF

  DrillingDOFs = GetLogical(SolverPars, 'Drilling DOFs', Found)
  IF (DrillingDOFs) CALL Warn('ShellSolver', &
      'Drilling DOFs do not support all options and alters the meaning of all rotational DOFs/BCs')
  IF (DrillingDOFs) THEN
    DrillingPar = GetConstReal(SolverPars, 'Drilling Stabilization Parameter', Found)
    IF (.NOT. Found) DrillingPar = 1.0d0
  ELSE
    DrillingPar = 1.0d0
  END IF

  ! For verification purposes we may solve a case for which the reference strain 
  ! energy is known:
  SolveBenchmarkCase = GetLogical(SolverPars, 'Benchmark Problem', Found)

  ! ---------------------------------------------------------------------------------
  ! Read parameters that control the nonlinear solution:
  ! ---------------------------------------------------------------------------------
  LargeDeflection = GetLogical(SolverPars, 'Large Deflection')
  IF (LargeDeflection) THEN
    MaxNonlinIters = ListGetInteger(SolverPars, 'Nonlinear System Max Iterations')
    NonlinTol =  GetConstReal(SolverPars, 'Nonlinear System Convergence Tolerance')
    Relaxation = ListGetCReal(SolverPars, 'Nonlinear System Relaxation Factor', Relax)
    Relax = Relax .AND. (ABS(Relaxation - 1.0_dp) > EPSILON(Relaxation))
    IF (DrillingDOFs) CALL Fatal('ShellSolver', &
        'Drilling DOFs cannot yet be combined with Large Deflection')
    SolveBenchmarkCase = .FALSE.
    IF (.NOT. ASSOCIATED(Solver % Matrix % BulkRHS)) &
        CALL AllocateVector(Solver % Matrix % BulkRHS, SIZE(Solver % Matrix % RHS))
    Solver % Matrix % BulkRHS = 0.0d0
  ELSE
    MaxNonlinIters = 0
  END IF
  NonlinearBending = GetLogical(SolverPars, 'Nonlinear Bending Strains', Found)
  IF (.NOT. Found) NonlinearBending = .TRUE.

  MeshDisplacementActive = GetLogical(SolverPars, 'Displace Mesh', Found)  
  
  HarmonicAssembly = EigenOrHarmonicAnalysis(Solver) .OR. GetLogical(SolverPars, &
      'Harmonic Mode', Found) .OR. GetLogical(SolverPars, 'Harmonic Analysis', Found)
  MassAssembly =  TransientSimulation .OR. HarmonicAssembly 

  ! ---------------------------------------------------------------------------------
  ! The choice of strain reduction method. Now only the automated default is active.
  ! Alter to experiment with other methods.
  ! ---------------------------------------------------------------------------------
  StrainReductionMethod = ListGetInteger(SolverPars, 'Strain Reduction Operator', &
      Found, minv=0, maxv=4)
  IF (.NOT.Found) StrainReductionMethod = AutomatedChoice
  MembraneStrainReductionMethod = ListGetInteger(SolverPars, 'Membrane Strain Reduction Operator', &
      Found, minv=0, maxv=4)
  IF (.NOT.Found) MembraneStrainReductionMethod = StrainReductionMethod

  IF (MembraneStrainReductionMethod /= NoStrainReduction) &
      MembraneStrainReductionMethod = AutomatedChoice
  IF (StrainReductionMethod /= NoStrainReduction) &
      StrainReductionMethod = AutomatedChoice

  ! ---------------------------------------------------------------------------------
  ! Parameters for shear/membrane/stretch relaxation: 
  ! ---------------------------------------------------------------------------------  
  ShearAlpha = ListGetCReal(SolverPars, 'Shear Relaxation Alpha', Found)
  IF (.NOT. Found) ShearAlpha = 1.0d0
  MembraneAlpha = ListGetCReal(SolverPars, 'Membrane Relaxation Alpha', Found)
  IF (.NOT. Found) MembraneAlpha = 0.0d0
  StretchAlpha = ListGetCReal(SolverPars, 'Stretch Relaxation Alpha', Found)
  IF (.NOT. Found) StretchAlpha = 0.0d0

  Bubbles = GetLogical(SolverPars, 'Bubbles', Found)
  RotateDOFs = GetLogical(SolverPars, 'Rotate DOFs', Found)

  !-----------------------------------------------------------------------------------
  ! The field variables for saving the orientation of lines of curvature basis
  ! vectors at the nodes: Since the global DOFs are expressed with respect to the
  ! global frame, creating these may not be of interest.
  !-----------------------------------------------------------------------------------
  SavePrincipalAxes = GetLogical(SolverPars, 'Principal Axes Output', Found)
  IF (SavePrincipalAxes) THEN
    NodalPDir1 => VariableGet(Mesh % Variables, 'Principal Coordinate Dir1')
    NodalPDir2 => VariableGet(Mesh % Variables, 'Principal Coordinate Dir2')
    NodalPDir3 => VariableGet(Mesh % Variables, 'Principal Coordinate Dir3')
    NodalPDir1 % Values = 0.0d0
    NodalPDir2 % Values = 0.0d0
    NodalPDir3 % Values = 0.0d0
    ! ---------------------------------------------------------------------
    ! The nodal visits array for averaging purposes:
    ! ---------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED(VisitsList) ) &
        CALL AllocateVector(VisitsList, SIZE(NodalPDir1 % Values)/3)
    VisitsList = 0
  END IF
  !
  ! Check if a 3-D elasticity solution (expected to have the name "Displacement")
  ! is available for computing the BCs of the shell model from it:
  !
  Displacement3D => VariableGet(Mesh % Variables, 'Displacement', .TRUE.)
  IF (ASSOCIATED(Displacement3D)) THEN
    ! Both the fields must be available in some common nodes to make this functional: 
    SolidShellCoupling = COUNT(Displacement3D % Perm > 0 .AND. &
        Solver % Variable % Perm > 0) > 0
  ELSE
    SolidShellCoupling = .FALSE.
  END IF

  IF (.NOT. ASSOCIATED(Indices)) ALLOCATE( Indices(Mesh % MaxElementDOFs) )
  IF (.NOT. ALLOCATED(LocalSol)) ALLOCATE( LocalSol(ShellModelPar, Mesh % MaxElementDOFs) )
  IF (.NOT. ALLOCATED(LocalRHSForce)) ALLOCATE( LocalRHSForce((ShellModelPar+1) * Mesh % MaxElementDOFs) )

  IF (.NOT. ASSOCIATED(TotalSol)) THEN
    CALL AllocateVector(TotalSol, SIZE(Solver % Variable % Values))
  ELSE
    IF (MeshDisplacementActive) THEN
      CALL Info('ShellSolver', 'Returning the mesh to its reference position', Level=4)     
      CALL DisplaceMesh(Mesh, Solver % Variable % Values, -1, Solver % Variable % Perm, &
         ShellModelPar, .FALSE., 3)      
    END IF
  END IF
  
  ComputeShellArea = .FALSE.
  SKIP_BLENDING: IF (.NOT. SkipBlending) THEN
    ! ---------------------------------------------------------------------------------
    ! PART I: 
    ! Get the director data at the nodes as a field variable 'Director' or
    ! read the director data at the nodes from mesh.director file and check the
    ! the integrity of the surface model. An elementwise property 'director' 
    ! corresponding to the data is created, if not already available
    ! via reading the director data from the file mesh.elements.data. 
    !----------------------------------------------------------------------------------
    Director => VariableGet(Mesh % Variables, 'Director', .TRUE.)
    CALL ReadSurfaceDirector(Mesh % Name, Mesh % NumberOfNodes, SolverPars, Director, &
        Solver % ActiveElements)
    CALL CheckSurfaceOrientation()

    ! --------------------------------------------------------------------------------
    ! PART II:
    ! Generate the descriptions of curved element edges for improved geometry 
    ! approximation. The implementation may not be memory efficient as data is 
    ! duplicated for shared element edges with the same director data. Here the
    ! variable CurveDataOutput can be used to output edge data into a file.
    ! ---------------------------------------------------------------------------------
    CurveDataOutput = GetLogical(SolverPars, 'Edge Curves Output', Found)
    CALL CreateCurvedEdges(CurveDataOutput)

    ! ---------------------------------------------------------------------------------
    ! Check whether the area of shell surface should be computed (here this is done
    ! in several ways to check the model integrity):
    ! ---------------------------------------------------------------------------------
    ComputeShellArea = GetLogical(SolverPars, 'Compute Shell Area', Found)
    BlendingSurfaceArea = 0.0d0
    MappedMeshArea = 0.0d0
  END IF SKIP_BLENDING

  ! ---------------------------------------------------------------------------------
  ! PART III:
  ! Perform a reparametrization to obtain lines of curvature coordinates and assemble 
  ! the discrete shell equations. It is also possible to utilize the parametrized edge 
  ! curves to obtain improved geometry approximation via using the finite element 
  ! blending technique. 
  ! ---------------------------------------------------------------------------------

  NONLINEARLOOP: DO NonlinIter=1,MaxNonlinIters+1

    IF (NonlinIter < (MaxNonlinIters+1)) THEN
      CALL Info('ShellSolver','--------------------------------------------------------', Level=4)
      WRITE( Message,'(A,I4)') 'Nonlinear iteration:', NonlinIter
      CALL Info('ShellSolver', Message, Level=4)
      CALL Info('ShellSolver','--------------------------------------------------------', Level=4)    
    END IF

    TotalSol(:) = Solver % Variable % Values(:)

    CALL DefaultInitialize()
    ShellModelArea = 0.0d0
    TotalErr = 0.0d0         ! Just for verification purposes (remove when final)
    Active = GetNOFActive()  

    ASSEMBLYLOOP: DO k=1,Active
      BGElement => GetActiveElement(k)

      Family = GetElementFamily(BGElement)
      IF ( .NOT.(Family == 3 .OR. Family == 4) ) CYCLE

      n  = GetElementNOFNodes()
      nd = GetElementDOFs(Indices)
      nb = GetElementNOFBDOFs()
      
      IF (LargeDeflection) THEN
        CALL GetVectorLocalSolution(LocalSol, USolver=Solver)
      ELSE
        LocalSol = 0.0d0
      END IF

      !----------------------------------------------------------------------
      ! Bubbles are designed for the lowest-order discretization:
      !----------------------------------------------------------------------
      ApplyBubbles = Bubbles .AND. (nd == Family)
      IF (nb > 0) CALL Fatal('ShellSolver', &
          'Static condensation for p-bubbles is not supported')

      REPARAMETRIZATION: IF (NonlinIter==1 .AND. (.NOT. SkipBlending .OR. ReparametrizeMesh)) THEN
        !----------------------------------------------------------------------
        ! Create elementwise geometry data related to the reference configuration 
        ! in order to reparametrize. This is computed only once since the data is 
        ! saved as elementwise properties and can thus be retrieved by calling 
        ! the function GetElementProperty.
        !----------------------------------------------------------------------
        IF (ReparametrizeMesh) THEN
          ! ----------------------------------------------------------------------
          ! In this case we just reparametrize the element of a given higher-order
          ! nodal mesh to obtain an orthogonal parametrization
          ! ----------------------------------------------------------------------
          ! Create a local coordinate frame whose orientation corresponds to 
          ! the orientation of lines of curvatures at the element center. Compute 
          ! also the coefficients of the Taylor polynomial for creating the improved 
          ! lines of curvature parameterization:
          !-------------------------------------------------------------------------
          CALL LinesOfCurvatureFrame(BGElement, TaylorApproximation=.TRUE., &
              LagrangeNodes=LocalFrameNodes, PlanarPoint=PlanarPoint, &
              UmbilicalPoint=UmbilicalPoint, SaveProperties=.TRUE., &
              ReparametrizeMesh=.TRUE.)
        ELSE
          !----------------------------------------------------------------------
          ! Here we apply the finite element blending technique to obtain
          ! an improved geometry approximation and then reparametrize
          !----------------------------------------------------------------------
          ! Get the elementwise average of director data for orientation purposes
          ! (check also for body flatness):
          !----------------------------------------------------------------------
          d = AverageDirector(BGElement, n, PlateBody)

          CALL LinesOfCurvatureFrame(BGElement, TaylorApproximation=.TRUE., &
              LagrangeNodes=LocalFrameNodes, d=d, PlanarSurface=PlateBody, &
              PlanarPoint=PlanarPoint, UmbilicalPoint=UmbilicalPoint, &
              SaveProperties=.TRUE.) 
        END IF
        
        TaylorParams => GetElementProperty('taylor parameters', BGElement)

        !--------------------------------------------------------------------------
        ! Obtain the final domain for improved lines of curvature parametrization. 
        ! The nodes of the final domain are here saved as the elementwise property
        ! 'patch nodes'.
        !--------------------------------------------------------------------------
        CALL LinesOfCurvaturePatch(BGElement, LocalFrameNodes(1:MaxPatchNodes,1:2), &
            TaylorParams, Family, PlanarPoint, UmbilicalPoint)

        !----------------------------------------------------------------------
        ! The area computation for the available geometry description:
        !----------------------------------------------------------------------
        IF (ComputeShellArea) CALL ComputeSurfaceArea(BGElement, BlendingSurfaceArea)
        !LocalFrameNodes(:,3) = ZNodes(:)
        !IF (ComputeShellArea) CALL MappedBGMeshArea(BGElement, LocalFrameNodes, MappedMeshArea)
      END IF REPARAMETRIZATION

      ! ------------------------------------------------------------------------------
      ! Generate the tangential stiffness matrix and assemble the local contribution:
      ! -----------------------------------------------------------------------------
      IF (CartesianFormulation) THEN
        CALL ShellLocalMatrixCartesian(BGElement, n, nd+nb, ShellModelPar, LocalSol, &
            LargeDeflection, NonlinearBending, MassAssembly, HarmonicAssembly, LocalRHSForce, &
            SkipBlending, ReparametrizeMesh, &
            BenchmarkProblem = GetLogical(SolverPars, 'Benchmark Problem', Found))
      ELSE
        CALL ShellLocalMatrix(BGElement, n, nd+nb, ShellModelPar, LocalSol, &
            LargeDeflection, NonlinearBending, StrainReductionMethod, MembraneStrainReductionMethod, &
            ShearAlpha, MembraneAlpha, StretchAlpha, ApplyBubbles, DrillingDOFs, DrillingPar, &
            RotateDOFs, MassAssembly, HarmonicAssembly, LocalRHSForce, ShellModelArea, TotalErr, &
            BenchmarkProblem=SolveBenchmarkCase)
      END IF

      IF (LargeDeflection .AND. NonlinIter == 1) THEN
        ! ---------------------------------------------------------------------------
        ! Create a RHS vector which contains just the contribution of external loads
        ! for the purpose of nonlinear error estimation:
        ! ---------------------------------------------------------------------------
        ValuesSaved => Solver % Matrix % RHS
        Solver % Matrix % RHS => Solver % Matrix % BulkRHS
        CALL DefaultUpdateForce(LocalRHSForce)
        Solver % Matrix % RHS => ValuesSaved
      END IF

      !-------------------------------------------------------------------------
      ! The following is not active since with the current implementation
      ! triangular blending functions cannot be evaluated at nodes (due to  
      ! implementation). Save principal directions at the nodes of the background 
      ! mesh:
      !-------------------------------------------------------------------------
      SavePrincipalAxes = .FALSE.
      IF (SavePrincipalAxes) THEN

        ShellElement => BGElement % Type

        DO j=1,n
          i = Indices(j)
          t = Solver % Variable % Perm(i)

          u = ShellElement % NodeU(j)
          v = ShellElement % NodeV(j)

          !-------------------------------------------------------------------------
          ! The principal directions at the given node: 
          !-------------------------------------------------------------------------        
          CALL LinesOfCurvatureFrame(BGElement, u, v, d1, d2, d3, p, d=d)

          !-------------------------------------------------------------------------
          ! Add to existing values and average later:
          !-------------------------------------------------------------------------
          NodalPDir1 % Values(3*(t-1)+1) = d1(1) + NodalPDir1 % Values(3*(t-1)+1)
          NodalPDir1 % Values(3*(t-1)+2) = d1(2) + NodalPDir1 % Values(3*(t-1)+2)
          NodalPDir1 % Values(3*(t-1)+3) = d1(3) + NodalPDir1 % Values(3*(t-1)+3)

          NodalPDir2 % Values(3*(t-1)+1) = d2(1) + NodalPDir2 % Values(3*(t-1)+1)
          NodalPDir2 % Values(3*(t-1)+2) = d2(2) + NodalPDir2 % Values(3*(t-1)+2)
          NodalPDir2 % Values(3*(t-1)+3) = d2(3) + NodalPDir2 % Values(3*(t-1)+3)

          NodalPDir3 % Values(3*(t-1)+1) = d3(1) + NodalPDir3 % Values(3*(t-1)+1)
          NodalPDir3 % Values(3*(t-1)+2) = d3(2) + NodalPDir3 % Values(3*(t-1)+2)
          NodalPDir3 % Values(3*(t-1)+3) = d3(3) + NodalPDir3 % Values(3*(t-1)+3)
          Visitslist(i) = Visitslist(i) + 1
        END DO
      END IF

    END DO ASSEMBLYLOOP

    !----------------------------------------------------------------------
    ! Add linearly elastic beams.
    !----------------------------------------------------------------------    
    ASSEMBLE_BEAMS: DO k=1,Active
      BGElement => GetActiveElement(k)

      Family = GetElementFamily(BGElement)
      IF (Family /= 2) CYCLE

      n  = GetElementNOFNodes()
      nd = GetElementDOFs(Indices)
      nb = GetElementNOFBDOFs()

      !----------------------------------------------------------------------
      ! We assume that p-element definitions are not empoyed and hard-code
      ! the bubble count:
      !----------------------------------------------------------------------
      nb = 1
      IF (.NOT.(n == 2 .AND. nd == 2)) CALL Fatal('ShellSolver', &
          'An unsupported 1-D element type or definition')

      IF (LargeDeflection) THEN
        CALL GetVectorLocalSolution(LocalSol, USolver=Solver)
      ELSE
        LocalSol = 0.0d0
      END IF

      CALL BeamStiffnessMatrix(BGElement, n, nd+nb, nb, TransientSimulation, MassAssembly, &
          HarmonicAssembly, LargeDeflection, LocalSol, LocalRHSForce, .TRUE., &
          ApplyRotation = .NOT.RotateDOFs, DrillingDOFs = DrillingDOFs) 

      IF (LargeDeflection .AND. NonlinIter == 1) THEN
        ! ---------------------------------------------------------------------------
        ! Create a RHS vector which contains just the contribution of external loads
        ! for the purpose of nonlinear error estimation:
        ! ---------------------------------------------------------------------------
        ValuesSaved => Solver % Matrix % RHS
        Solver % Matrix % RHS => Solver % Matrix % BulkRHS
        CALL DefaultUpdateForce(LocalRHSForce)
        Solver % Matrix % RHS => ValuesSaved
      END IF
    END DO ASSEMBLE_BEAMS
    
    CALL DefaultFinishBulkAssembly() 


    Active = GetNOFBoundaryElements()

    BOUNDARY_ASSEMBLY: DO k=1,Active
      BGElement => GetBoundaryElement(k)
      Family = GetElementFamily(BGElement)

      IF (ActiveBoundaryElement() .AND. Family <= 2) THEN
        n  = GetElementNOFNodes(BGElement)
        nd = GetElementNOFDOFs(BGElement)
        nd_parent = 0

        IF (LargeDeflection) THEN
          CALL GetVectorLocalSolution(LocalSol, USolver=Solver)
          
          Parent => GetBulkElementAtBoundary(BGElement)
          IF (ASSOCIATED(Parent)) THEN 
            nd_parent = GetElementNOFDOFs(Parent)
          END IF
        ELSE
          LocalSol = 0.0d0
        END IF

        CALL ShellBoundaryMatrix(BGElement, n, nd, ShellModelPar, LargeDeflection, &
            MassAssembly, HarmonicAssembly, LocalSol, LocalRHSForce, Parent, nd_parent, &
            CartesianFormulation, SkipBlending)

        IF (LargeDeflection .AND. NonlinIter == 1) THEN
          ! ---------------------------------------------------------------------------
          ! Assemble the RHS vector which contains just the contribution of external loads
          ! for the purpose of nonlinear error estimation:
          ! ---------------------------------------------------------------------------
          ValuesSaved => Solver % Matrix % RHS
          Solver % Matrix % RHS => Solver % Matrix % BulkRHS
          CALL DefaultUpdateForce(LocalRHSForce)
          Solver % Matrix % RHS => ValuesSaved
        END IF
      END IF
    END DO BOUNDARY_ASSEMBLY

    CALL DefaultFinishBoundaryAssembly()

    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()
    IF (SolidShellCoupling) CALL SetSolidCouplingBCs(Model, Solver, Displacement3D)

    ! ---------------------------------------------------------------------------------
    ! The solution variable is the solution increment while the sif-file specifies
    ! the Dirichlet BCs for the complete field. Modify BCs so that the right BC
    ! is obtained for the solution increment.
    !
    ! NOTE: If higher-order elements were used over the lowest-order background mesh, 
    ! the treatment of Dirichlet BCs should be checked (depending on how the additional
    ! DOFs would be created)
    ! --------------------------------------------------------------------------------
    IF (ALLOCATED(Solver % Matrix % ConstrainedDOF)) THEN
      DO i=1,Solver % Matrix % NumberOfRows
        IF (Solver % Matrix % ConstrainedDOF(i)) THEN
          Solver % Matrix % DValues(i) = Solver % Matrix % DValues(i) - Solver % Variable % Values(i)
        END IF
      END DO
      CALL EnforceDirichletConditions(Solver, Solver % Matrix, Solver % Matrix % RHS)
    END IF
 
    ! ---------------------------------------------------------------------------------
    ! Check whether the nonlinear iteration can be terminated:
    ! ---------------------------------------------------------------------------------
    IF (LargeDeflection) THEN
      IF (NonlinIter == 1) THEN

        IF (Parallel) THEN
          IF (.NOT. ASSOCIATED(Solver % Matrix % ParMatrix)) &
              CALL ParallelInitMatrix(Solver, Solver % Matrix)

          PMatrix => Solver % Matrix % ParMatrix % SplittedMatrix % InsideMatrix
          IF (.NOT. ASSOCIATED(PMatrix % RHS)) &
               ALLOCATE(PMatrix % RHS(PMatrix % NumberOfRows))

          ! Temporarily set the parallel rhs vector to be the plain source vector:
          CALL ParallelUpdateRHS(Solver % Matrix, Solver % Matrix % BulkRHS)
          Pb => PMatrix % RHS
          Norm = MAXVAL(ABS(Pb))
          Norm = ParallelReduction(Norm,2)
        ELSE
          Norm = MAXVAL(ABS(Solver % Matrix % BulkRHS(:)))
        END IF

        NoTractions = Norm < AEPS

        IF (NoTractions) THEN
          ! This appears to be a purely BC-loaded case, switch to using a different criterion
          ! (use absolute norm, this can be hard ...):
          CALL Info('ShellSolver', 'No distributed loads ... ', Level=4)
          CALL Info('ShellSolver', &
              'Switch to using absolute norm in the nonlinear error estimation',  Level=4)
          CALL Info('ShellSolver', &
              'This may give a hard stopping criterion',  Level=4)
          NonlinRes0 = 1.0d0
        ELSE
          ! Compute the 2-norm of the initial residual (RHS vector before setting BCs). 
          IF (Parallel)  THEN
            Norm = 0.0d0
            DO i=1,PMatrix % NumberOfRows
              Norm = Norm + Pb(i)**2
            END DO
            NonlinRes0 = SQRT(ParallelReduction(Norm))
          ELSE
            NonlinRes0 = SQRT(SUM(Solver % Matrix % BulkRHS(:)**2))
          END IF
        END IF
      END IF

      ! Employ BulkRHS vector to estimate the size of the current residual (RHS):
      Solver % Matrix % BulkRHS = Solver % Matrix % RHS

      IF (Parallel) THEN
        CALL ParallelUpdateRHS(Solver % Matrix, Solver % Matrix % BulkRHS)
        Norm = 0.0d0
        DO i=1,PMatrix % NumberOfRows
          Norm = Norm + Pb(i)**2
        END DO
        NonlinRes = SQRT(ParallelReduction(Norm)) / NonlinRes0
      ELSE
        NonlinRes = SQRT(SUM(Solver % Matrix % RHS(:)**2)) / NonlinRes0
      END IF
      WRITE(Message,'(a,I4,ES12.3)') 'Residual for nonlinear iterate', &
          NonlinIter-1, NonLinRes
      CALL Info('ShellSolver', Message, Level=3)        
      IF (NonlinRes < NonlinTol) THEN
        WRITE(Message,'(a)') 'Nonlinear iteration is terminated succesfully'
        CALL Info('ShellSolver', Message, Level=3)          
        EXIT
      ELSE IF (NonlinIter > MaxNonlinIters) THEN
        IF (GetLogical(SolverPars, 'Nonlinear System Abort Not Converged', Found)) THEN
          CALL Fatal('ShellSolver', 'Nonlinear iteration did not converge to tolerance')
        ELSE
          CALL Info('ShellSolver', 'Nonlinear iteration did not converge to tolerance', Level=6)
          EXIT
        END IF
      END IF
    END IF

    ! --------------------------------------------------------------------------------
    ! Previous correction may not be a particularly good initial guess so start from
    ! the trivial iterate:
    ! --------------------------------------------------------------------------------
    IF (LargeDeflection) Solver % Variable % Values = 0.0d0

    Norm = DefaultSolve()

    IF (LargeDeflection) THEN
      IF (Relax) THEN
        Solver % Variable % Values(:) = TotalSol(:) + Relaxation * Solver % Variable % Values(:)
      ELSE
        Solver % Variable % Values(:) = TotalSol(:) + Solver % Variable % Values(:)
      END IF
    END IF
  END DO NONLINEARLOOP

  ! -------------------------------------------------------------------------------
  ! Finalize the generation of the principal directions (average):
  ! -------------------------------------------------------------------------------
  IF (SavePrincipalAxes) THEN
    DO i=1,SIZE(VisitsList)  
      n = VisitsList(i)
      IF (n>1) THEN
        t = Solver % Variable % Perm(i)
        e1(1) = NodalPDir1 % Values(3*(t-1)+1)/n
        e1(2) = NodalPDir1 % Values(3*(t-1)+2)/n
        e1(3) = NodalPDir1 % Values(3*(t-1)+3)/n
        NodalPDir1 % Values(3*(t-1)+1:3*(t-1)+3) = e1(1:3)/Sqrt(SUM(e1(:)**2))

        e2(1) = NodalPDir2 % Values(3*(t-1)+1)/n
        e2(2) = NodalPDir2 % Values(3*(t-1)+2)/n
        e2(3) = NodalPDir2 % Values(3*(t-1)+3)/n
        NodalPDir2 % Values(3*(t-1)+1:3*(t-1)+3) = e2(1:3)/Sqrt(SUM(e2(:)**2))

        e3(1) = NodalPDir3 % Values(3*(t-1)+1)/n
        e3(2) = NodalPDir3 % Values(3*(t-1)+2)/n
        e3(3) = NodalPDir3 % Values(3*(t-1)+3)/n
        NodalPDir3 % Values(3*(t-1)+1:3*(t-1)+3) = e3(1:3)/Sqrt(SUM(e3(:)**2))
      END IF
    END DO
  END IF


  ! -------------------------------------------------------------------------------------
  ! PART IV: Postprocess
  ! -------------------------------------------------------------------------------------
  IF ( MeshDisplacementActive ) THEN
     CALL Info('ShellSolver', 'Displacing the mesh with computed displacement field', Level=4)
     CALL DisplaceMesh(Mesh, Solver % Variable % Values, 1, Solver % Variable % Perm, &
         ShellModelPar, .FALSE., 3)
  END IF


  ! -------------------------------------------------------------------------------------
  ! SOME VERIFICATION OUTPUT if a benchmark case of straight cylindrical shell is solved
  !-----------------------------------------------------------------------------------
  IF (SolveBenchmarkCase .AND. .NOT.Parallel) THEN
    Work = 8.0d0*SUM(Solver % Variable % Values(:) * Solver % Matrix % RHS(:))
    PRINT *, 'Energy representation coefficient = ', Work/(12.0d0*(1.0d0-(1.0d0/3.0d0)**2)*(1.0d5)**2/7.0d10 * (1.0d-2)**3)
    i = ListGetInteger(SolverPars, 'Benchmark Case', Found, minv=1,maxv=2)
    IF (Found) THEN
      !CALL MatrixVectorMultiply(Solver % Matrix, Solver % Variable % Values, TotalSol)
      !Work = 8.0d0 * SUM( Solver % Variable % Values(:) * TotalSol(:) )
      !Work = 8.0d0*SUM(Solver % Variable % Values(:) * Solver % Matrix % RHS(:))

      SELECT CASE(i)
      CASE(1)
        RefWork = 12.0d0*(1.0d0-(1.0d0/3.0d0)**2)*(1.0d5)**2/7.0d10 * 2.688287959059254d0 * 1.0d-2 ! t=0.01
        !RefWork = 12.0d0*(1.0d0-(1.0d0/3.0d0)**2)*(1.0d9)**2/7.0d10 * 1.828629366566552 * 1.0d-1 ! t=0.1
      CASE(2)
        RefWork = 12.0d0*(1.0d0-(1.0d0/3.0d0)**2)*(1.0d5)**2/7.0d10 * 0.704331198817278d0 * (1.0d-2)**3 ! t=0.01
      END SELECT
      PRINT *, 'Relative energy error = ', SQRT(ABS(RefWork-Work)/RefWork)
      PRINT *, 'Total number of DOFS = ', SIZE(Solver % Variable % Values) 
    END IF

    IF (ComputeShellArea) THEN
      RefArea = 0.5d0 * PI 
      !RefArea = 4 * (1.0472d0)**2  
      PRINT *, 'Relative Error of Model Surface Area = ', ABS(RefArea  - ShellModelArea)/RefArea    
      PRINT *, 'Relative Error of Blending Surface Area = ', ABS(RefArea  - BlendingSurfaceArea)/RefArea
      !PRINT *, 'Relative Error of Mapped BG Mesh Area = ', ABS(RefArea  - MappedMeshArea)/RefArea
    END IF

    !PRINT *, 'Mean curvature L2-error = ', SQRT(TotalErr)

    !MaxPDir1Err = 0.0d0
    !DO i=1,SIZE(NodalPDir1 % Values)/3
    !  print *, NodalPDir1 % Values(3*(i-1)+1:3*(i-1)+3)
    !  MaxPDir1Err = MaxPDir1Err + NodalPDir1 % Values(3*(i-1)+1)**2 + &
    !      (-1.0d0 - NodalPDir1 % Values(3*(i-1)+2))**2 + NodalPDir1 % Values(3*(i-1)+3)**2
    !END DO
    !PRINT *, 'The L2 error norm for the principal direction 1 = ', SQRT(MaxPDir1Err)

  END IF

 
CONTAINS

! ---------------------------------------------------------------------------------
! This subroutine uses an ordinary field variable or mesh.director file arranged as
!
!    node_id1 d_x d_y d_z
!    ...
!    node_idN d_x d_y d_z
!
! to obtain the shell director data at nodes and creates an elementwise property 
! 'director' corresponding to this data. If the file mesh.elements.data has been 
! used to specify the director as an elementwise property 'director', the director 
! obtained from mesh.elements.data is used. With the keyword Write Elemental Director
! being active, the director data is written as elementwise property to a file whose 
! format conforms with a file mesh.elements.data (this is the default name for the 
! output file, so this option can be used to convert mesh.director into 
! mesh.elements.data format).
! 
! Note: Parallel file formats for mesh.elements.data and mesh.director have not
!       been implemented. Parallel execution is thus possible only when the
!       director is available as an ordinary field variable.
!------------------------------------------------------------------------------
  SUBROUTINE ReadSurfaceDirector(MeshName, NumberOfNodes, SolverPars, Director, &
      ActiveElements)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    CHARACTER(LEN=MAX_NAME_LEN), INTENT(IN) :: MeshName
    INTEGER, INTENT(IN) :: NumberOfNodes
    TYPE(ValueList_t), POINTER, INTENT(IN) :: SolverPars
    TYPE(Variable_t), POINTER, INTENT(IN) :: Director
    INTEGER, POINTER :: ActiveElements(:)
    !------------------------------------------------------------------------------
    LOGICAL :: UseFieldVariable, ReadNodalDirectors, WriteElementsData, Found
    INTEGER :: n, iostat, i, j, k, i0, NumberOfLines, Family
    INTEGER, POINTER :: InvPerm(:)
    REAL(KIND=dp), POINTER :: NodalDirector(:,:)  
    REAL(KIND=dp), POINTER :: DirectorValues(:)
    REAL(KIND=dp) :: ElementDirectors(3*MaxBGElementNodes)
    CHARACTER(LEN=MAX_NAME_LEN) :: DirectorFile, FormatString
    !------------------------------------------------------------------------------
    ReadNodalDirectors = .FALSE.

    UseFieldVariable = ASSOCIATED(Director)
    IF (UseFieldVariable) THEN
      CALL Info('ReadSurfaceDirector', '&
          Using the field Director to define the mid-surface normal', Level=4)
      IF (Director % DOFs /= 3) CALL Fatal('ReadSurfaceDirector', &
          'The director field should have three components')
      IF (.NOT.ASSOCIATED(Director % Perm) .OR. .NOT.ASSOCIATED(Director % Values)) &
          CALL Fatal('ReadSurfaceDirector', 'The director solution is not associated')
    ELSE
      ! -----------------------------------------------------------------------------
      ! Check whether mesh.director can be read:
      ! -----------------------------------------------------------------------------
      n = LEN_TRIM(MeshName)
      DirectorFile = TRIM(MeshName)//'/'//'mesh.director'//CHAR(0)
      
      INQUIRE(FILE = DirectorFile(1:n+15), EXIST = ReadNodalDirectors)

      IF (ReadNodalDirectors) THEN
        OPEN(10, FILE = DirectorFile(1:n+15), status='OLD', IOSTAT = iostat)
        IF ( iostat /= 0 ) CALL Fatal('ReadSurfaceDirector', &
            'Opening mesh.director file failed.')

        ! Director data may not have been defined in all mesh nodes.  
        ! Find out how many director values can be read:
        NumberOfLines = 0
        DO WHILE (.TRUE.)
          READ(10,*, IOSTAT=iostat, END=100) k, d
          NumberOfLines = NumberOfLines + 1
        END DO
100     REWIND(10)
      END IF
    END IF

    IF (UseFieldVariable .OR. ReadNodalDirectors) THEN
      IF (ReadNodalDirectors) THEN
        CALL AllocateArray(NodalDirector, NumberOfLines, 3, 'ReadSurfaceDirector', &
            'NodalDirector array could not be allocated')
        CALL AllocateVector(InvPerm, NumberOfNodes, 'ReadSurfaceDirector', &
            'InvPerm array could not be allocated')
      
        DO i=1,NumberOfLines
          READ(10,*,IOSTAT=iostat) k, d
          InvPerm(k) = i 
          Norm = SQRT(SUM(d(1:3)**2))
          NodalDirector(i,1) = d(1)/Norm
          NodalDirector(i,2) = d(2)/Norm
          NodalDirector(i,3) = d(3)/Norm
        END DO
        CLOSE(10)
      END IF
      ! ---------------------------------------------------------------------
      ! Create director data as elementwise property
      ! ---------------------------------------------------------------------
      Active = GetNOFActive()
      DO k=1,Active
        Element => GetActiveElement(k)
        Family = GetElementFamily(Element)
        ! -------------------------------------------------------------------
        ! The director data is expected for surface elements only:
        ! -------------------------------------------------------------------
        IF (Family < 3) CYCLE
        ! -------------------------------------------------------------------
        ! If mesh.elements.data has defined the director, respect that data:
        ! -------------------------------------------------------------------
        DirectorValues => GetElementProperty('director', Element)
        IF (ASSOCIATED(DirectorValues)) CYCLE

        n  = GetElementNOFNodes()
        IF (ReadNodalDirectors) THEN
          DO i=1,n
            i0 = (i-1)*3
            ElementDirectors(i0+1:i0+3) = NodalDirector(InvPerm(Element % NodeIndexes(i)),1:3)
          END DO
        ELSE
          DO i=1,n
            i0 = (i-1)*3
            j = 3*(Director % Perm(Element % NodeIndexes(i)) - 1)
            ElementDirectors(i0+1:i0+3) = Director % Values(j+1:j+3)
          END DO
        END IF
        CALL SetElementProperty('director', ElementDirectors(1:3*n), Element)
      END DO
    END IF

    ! ---------------------------------------------------------------------
    ! Write the director data as elementwise property to a file whose
    ! format conforms with a file mesh.elements.data. By default
    ! the file name mesh.elements.data is used. This never overwrites
    ! an existing file.
    ! ---------------------------------------------------------------------    
    WriteElementsData = GetLogical(SolverPars, 'Write Elemental Director', Found)

    IF ( WriteElementsData ) THEN
      OutputFile = GetString(SolverPars, 'Elemental Director Output File', Found)
      IF (.NOT. Found) OutputFile = 'mesh.elements.data'//CHAR(0)

      n = LEN_TRIM(MeshName)
      DirectorFile = MeshName(1:n)//'/'//TRIM(OutputFile)//CHAR(0)

      n = LEN_TRIM(DirectorFile)
      INQUIRE(FILE = DirectorFile(1:n), EXIST = Found)
      IF (Found) THEN
        CALL Info('ReadSurfaceDirector', &
            'a file for director output exists: write rejected', Level=5)
      ELSE
        OPEN(10, FILE = DirectorFile(1:n), status='NEW', IOSTAT = iostat)        
        IF ( iostat /= 0 ) CALL Fatal( 'ReadSurfaceDirector', &
            'Opening a file for elementwise director output failed.')

        Active = GetNOFActive()
        DO k=1,Active
          Element => GetActiveElement(k)
          DirectorValues => GetElementProperty('director', Element)

          IF (ASSOCIATED(DirectorValues)) THEN
            n  = GetElementNOFNodes()
            IF (SIZE(DirectorValues) < 3*n) CALL Fatal('ReadSurfaceDirector', &
                'Elemental director data is not associated with all nodes')
 
            !WRITE(FormatString(1:1),'(A1)') '('
            !IF (3*n < 10) THEN
            !  WRITE(FormatString(2:2),'(A1)') I2S(3*n)
            !  i0 = 2
            !ELSE
            !  WRITE(FormatString(2:3),'(A2)') I2S(3*n)
            !  i0 = 3
            !END IF
            !WRITE(FormatString(i0+1:i0+1),'(A1)') '('
            !WRITE(FormatString(i0+2:i0+10),'(A9)') '2x,E22.15'
            !WRITE(FormatString(i0+11:i0+12),'(A2)') '))'

            !WRITE(10,'(A8,I0)') 'element:', ActiveElements(k)
            !WRITE(10,'(A9)',ADVANCE='NO') 'director:'
            !WRITE(10,FormatString(1:i0+12)) DirectorValues(1:3*n)
            !WRITE(10,'(A3)') 'end'

            WRITE(FormatString,'(A)') '(A,I0,A,'//I2S(3*n)//'E22.15,A)'
            WRITE(10,FormatString) 'element: ',ActiveElements(k),' director: ', &
                DirectorValues(1:3*n),' end'            
          ELSE
            CALL Fatal('ReadSurfaceDirector', 'Elemental director data is not associated')
          END IF
        END DO
        CLOSE(10)
      END IF
    END IF

    IF (ReadNodalDirectors) DEALLOCATE(NodalDirector, InvPerm)
!------------------------------------------------------------------------------
  END SUBROUTINE ReadSurfaceDirector
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! This function can be used to return the elementwise values of the director
! field. The director data is supposed to be found as the elementwise property
! 'director'. If this property does not exits, the normal is computed otherwise.
!-------------------------------------------------------------------------------
  FUNCTION GetElementalDirector(Element, ElementNodes) RESULT(DirectorValues) 
!-------------------------------------------------------------------------------    
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    TYPE(Nodes_t), OPTIONAL, INTENT(IN) :: ElementNodes
    REAL(KIND=dp), POINTER :: DirectorValues(:)
    !-------------------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Visited = .FALSE., UseElementProperty = .FALSE., UseNormalSolver = .FALSE.
    REAL(KIND=dp), POINTER :: NodalNormals(:)
    REAL(KIND=dp) :: Normal(3)
    INTEGER :: n
    
    SAVE Visited, UseElementProperty, NodalNormals, Nodes
    !-------------------------------------------------------------------------------

    IF (.NOT. Visited) THEN
      DirectorValues => GetElementProperty('director', Element)
      UseElementProperty = ASSOCIATED( DirectorValues ) 

      IF (.NOT. UseElementProperty) THEN
        n = CurrentModel % MaxElementNodes
        ALLOCATE( NodalNormals(3*n) ) 
      END IF
      Visited = .TRUE.
    END IF

    IF ( UseElementProperty ) THEN    
      DirectorValues => GetElementProperty('director', Element)
    ELSE
      IF( PRESENT( ElementNodes ) ) THEN
        Normal = NormalVector( Element, ElementNodes, Check = .TRUE. ) 
      ELSE
        CALL GetElementNodes( Nodes, Element ) 
        Normal = NormalVector( Element, Nodes, Check = .TRUE. ) 
      END IF
        
      n = Element % TYPE % NumberOfNodes
      NodalNormals(1:3*n:3) = Normal(1)
      NodalNormals(2:3*n:3) = Normal(2)
      NodalNormals(3:3*n:3) = Normal(3)      
      DirectorValues => NodalNormals
    END IF     
!-------------------------------------------------------------------------------    
  END FUNCTION GetElementalDirector
!-------------------------------------------------------------------------------
  
! ---------------------------------------------------------------------------
!> Perform an additional check that the director data defines a properly 
!> oriented model. All directors should point to the same side of the surface.
!----------------------------------------------------------------------------
  SUBROUTINE CheckSurfaceOrientation()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    INTEGER :: n, i, j, k, i0, Active, Family
    REAL(KIND=dp), POINTER :: NodalDirector(:,:)  
    REAL(KIND=dp), POINTER :: DirectorValues(:)
    REAL(KIND=dp) :: d(3), d1(3), d2(3), X1(3), X2(3)
    REAL(KIND=dp) :: e1(3), e2(3), e3(3), Norm
    !------------------------------------------------------------------------------

    Active = GetNOFActive()
    DO k=1,Active
      Element => GetActiveElement(k)
      Family = GetElementFamily(Element)
      IF (Family < 3) CYCLE
      n  = GetElementNOFNodes()
      CALL GetElementNodes( Nodes )

      i = 1
      j = 2
      X1(1) = Nodes % x(i)
      X1(2) = Nodes % y(i)
      X1(3) = Nodes % z(i)
      X2(1) = Nodes % x(j)
      X2(2) = Nodes % y(j)
      X2(3) = Nodes % z(j)
      e1(:) = X2(:)-X1(:)
      Norm = SQRT(SUM(e1(:)**2))
      e1 = e1/Norm

      SELECT CASE(Family)
      CASE(3)
        i = 1
        j = 3
      CASE(4)
        i = 2
        j = 3
      CASE DEFAULT
        CYCLE
      END SELECT

      X1(1) = Nodes % x(i)
      X1(2) = Nodes % y(i)
      X1(3) = Nodes % z(i)
      X2(1) = Nodes % x(j)
      X2(2) = Nodes % y(j)
      X2(3) = Nodes % z(j)
      e2(:) = X2(:)-X1(:)
      Norm = SQRT(SUM(e2(:)**2))

      ! Now, define the element surface orientation: 
      e3(:) = CrossProduct(e1,e2) 
      Norm = SQRT(SUM(e3(:)**2))
      e3 = e3/Norm     

      ! Check that all directors point to the same side of the oriented surface:
      DirectorValues => GetElementalDirector(Element, Nodes)

      IF (.NOT. ASSOCIATED(DirectorValues)) THEN
        CALL Fatal('CheckSurfaceOrientation', 'Elemental director data is not associated')
      END IF
      IF (SIZE(DirectorValues) < 3*n) CALL Fatal('CheckSurfaceOrientation', &
          'Elemental director data is not associated with all nodes')

      ! reference direction for the 1st element node
      d1(1:3) = DirectorValues(1:3)
      c = DOT_PRODUCT(d1,e3)

      DO j=2,n
        i0 = (j-1)*3
        d2(1:3) = DirectorValues(i0+1:i0+3)
        IF ( (c * DOT_PRODUCT(d2,e3)) < 0.0d0 ) THEN
          PRINT *, 'Element indices=', Element % NodeIndexes(1:n)
          PRINT *, 'Reference normal =', e3(:)
          PRINT *, 'Node Index = ', j, Element % NodeIndexes(j)
          PRINT *, 'Director =', d2(:)
          CALL Fatal('CheckSurfaceOrientation', &
              'Director data does not define a unique upper/lower surface.')
        END IF
      END DO
    END DO
!-------------------------------------------------------------------------------   
  END SUBROUTINE CheckSurfaceOrientation
!-------------------------------------------------------------------------------   
  
! ---------------------------------------------------------------------------------
! Use nodal directors, which are retrieved as elementwise property 'director', 
! to create the parametrizations of curved edges for the Hermite interpolation.
! The edge curve data are written as elementwise properties 'edge frames' and
! 'edge parameters'. With the MacroElement option we may create additional space 
! curves corresponding to subtriangulations of quadrilateral elements.
!-------------------------------------------------------------------------------
  SUBROUTINE CreateCurvedEdges( FileOutput, MacroElements )
!-------------------------------------------------------------------------------
    IMPLICIT NONE
    
    LOGICAL, INTENT(IN) :: FileOutput
    LOGICAL, OPTIONAL, INTENT(IN) :: MacroElements
!-------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: QuadraticGeometryData, Subtriangulation
    INTEGER :: Active, k, e, i, j, l, v1, v2, v3, i0, j0, k0
    INTEGER :: Family, EdgesParametrized, CurveDataSize
    REAL(KIND=dp), POINTER :: A(:,:), cpars(:)
    REAL(KIND=dp), POINTER :: DirectorValues(:)
    REAL(KIND=dp), TARGET :: FrameData(3,4)
    REAL(KIND=dp), TARGET :: CurveData(CurveDataSize2)
    REAL(KIND=dp) :: d1(3), d2(3), d3(3), X1(3), X2(3), X3(3)
    REAL(KIND=dp) :: EdgeFramesData(FrameDataSize*MaxNumberOfCurves)
    REAL(KIND=dp) :: EdgeCurveParams(CurveDataSize2*MaxNumberOfCurves)
!-------------------------------------------------------------------------------

    A => FrameData(:,:)
    cpars => CurveData(:)

    ! Write edge curve parameters to a file:
    ! ------------------------------------------------------------------
    IF (FileOutput) OPEN(10, FILE = 'edgecsys.dat', status='REPLACE')

    Active = GetNOFActive()
    DO k=1,Active
      Element => GetActiveElement(k)
      Family = GetElementFamily(Element)
      IF (Family < 3) CYCLE
      CALL GetElementNodes(Nodes)

      DirectorValues => GetElementalDirector(Element, Nodes)

      ! Set some default values:
      Subtriangulation = .FALSE.
      EdgesParametrized = Family

      SELECT CASE(Family)
      CASE(3)
        IF (Element % TYPE % NumberOfNodes > 6) &
            CALL Fatal('CreateCurvedEdges', 'Triangular background mesh of order k>2 is not supported')
        QuadraticGeometryData = Element % TYPE % NumberOfNodes == 6
        IF (QuadraticGeometryData) &
            CALL Fatal('CreateCurvedEdges', 'Triangular 6-node background elements are not yet supported')
      CASE(4)
        IF (Element % TYPE % NumberOfNodes > 9) &
            CALL Fatal('CreateCurvedEdges', 'Background mesh of order k>2 is not supported')
        IF (Element % TYPE % NumberOfNodes == 8) &
            CALL Fatal('CreateCurvedEdges', '8-node background quad is not supported')
        QuadraticGeometryData = Element % TYPE % NumberOfNodes == 9
        IF (QuadraticGeometryData) THEN
          CALL Fatal('CreateCurvedEdges', '9-node background elements are not yet supported')
          EdgesParametrized = 6   ! Even 8 edges could be created
        ELSE
          ! -------------------------------------------------------------------------------
          ! We may consider 404 as a macroelement for a subtriangulation
          ! -------------------------------------------------------------------------------
          IF (PRESENT(MacroElements)) Subtriangulation = MacroElements
          IF (Subtriangulation) EdgesParametrized = 6
        END IF
      CASE DEFAULT
        CYCLE
      END SELECT

      IF (QuadraticGeometryData) THEN
        CurveDataSize = CurveDataSize2
      ELSE
        CurveDataSize = CurveDataSize1
      END IF

      DO e=1,EdgesParametrized
        !-------------------------------------------------------------------------
        ! First define edge orientation convention. 
        !-------------------------------------------------------------------------
        SELECT CASE(Family)
        CASE(3)
          SELECT CASE(e)
          CASE(1)
            i = 1
            j = 2
            l = 4
          CASE(2)
            i = 2
            j = 3
            l = 5
          CASE(3)
            i = 3
            j = 1
            l = 6
          END SELECT
        CASE(4)
          SELECT CASE(e)
          CASE(1)
            i = 1
            j = 2
            l = 5
          CASE(2)
            i = 2
            j = 3
            l = 6
          CASE(3)
            i = 4
            j = 3
            l = 7
          CASE(4)
            i = 1
            j = 4
            l = 8
          CASE(5)
            IF (Subtriangulation) THEN
              i = 1
              j = 3
            ELSE
              i = 8
              j = 6
              l = 9
            END IF
          CASE(6)
            IF (Subtriangulation) THEN
              i = 2
              j = 4
            ELSE         
              i = 5
              j = 7
              l = 9
            END IF
          END SELECT
        CASE DEFAULT
          CYCLE
        END SELECT

        X1(1) = Nodes % x(i)
        X1(2) = Nodes % y(i)
        X1(3) = Nodes % z(i)
        X2(1) = Nodes % x(j)
        X2(2) = Nodes % y(j)
        X2(3) = Nodes % z(j)

        v1 = Element % NodeIndexes(i)
        v2 = Element % NodeIndexes(j)

        i0 = (i-1)*3
        d1(1:3) = DirectorValues(i0+1:i0+3)
        i0 = (j-1)*3
        d2(1:3) = DirectorValues(i0+1:i0+3)

        ! ----------------------------------------------------------------------
        ! Construct data for creating the Hermite interpolation approximation of 
        ! the curved edge by using the nodal coordinates and director data. Two 
        ! nodes per edge yield third-order polynomial approximation while three 
        ! nodes per edge gives the fifth-order polynomial fit. 
        ! ----------------------------------------------------------------------
        IF (QuadraticGeometryData) THEN
          ! -----------------------------------------------------------
          ! In this case each edge has also a mid-node: 
          ! TO DO: Call HermiteForm instead of EdgeFrame
          ! -----------------------------------------------------------
          X3(1) = Nodes % x(l)
          X3(2) = Nodes % y(l)
          X3(3) = Nodes % z(l)

          i0 = (l-1)*3
          d3(1:3) = DirectorValues(i0+1:i0+3)

          CALL EdgeFrame(X1, X2, d1, d2, A, cpars, X3, d3)
        ELSE
          CALL HermiteForm(X1, X2, d1, d2, cpars)

          ! TO DO: Revise the file output
          !IF (FileOutput .AND. EdgesParametrized == 4) THEN
          !  WRITE(10, '(2I5,18e23.15)',ADVANCE='NO') v1, v2, A(1,1), A(2,1), &
          !      A(3,1), A(1,2), A(2,2), A(3,2), A(1,3), A(2,3), A(3,3), A(1,4), &
          !      A(2,4), A(3,4), cpars(1), cpars(1), cpars(2), cpars(3), cpars(4), &
          !      cpars(5)
          !  WRITE(10, *) ''
          !END IF
        END IF

        ! Prepare for writing the edge curve data as element properties:
        !---------------------------------------------------------------
        IF (QuadraticGeometryData) THEN
          i0 = (e-1)*FrameDataSize
          DO j=1,4
            j0 = i0 + (j-1)*3
            EdgeFramesData(j0+1:j0+3) = A(1:3,j)
          END DO
        END IF

        i0 = (e-1)*CurveDataSize
        EdgeCurveParams(i0+1:i0+CurveDataSize) = cpars(1:CurveDataSize)

      END DO

      ! Write the edge curve data as element properties:
      !----------------------------------------------------
      CALL SetElementProperty('edge parameters', &
          EdgeCurveParams(1:CurveDataSize*EdgesParametrized), Element)
      IF (QuadraticGeometryData) THEN
        CALL SetElementProperty('edge frames', &
            EdgeFramesData(1:FrameDataSize*EdgesParametrized), Element) 
      END IF
    END DO

    IF (FileOutput) CLOSE(10)
!------------------------------------------------------------------------------
  END SUBROUTINE CreateCurvedEdges
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Create an orthonormal basis (ex,ey,ez) for a local coordinate 
! system (x,y,z) associated with a given edge. The base vectors and the origin 
! X0 of the local system are returned via A. The length he of the line segment 
! for parameterizing the curved edge and the nodal curve parameters for
! the Hermite interpolation are also returned via cpars. The input data 
! Xk gives the global coordinates of the kth vertex on the edge, while dk 
! specifies the director at the kth vertex on the edge.
!------------------------------------------------------------------------------
  SUBROUTINE EdgeFrame(X1, X2, d1, d2, A, cpars, X3, d3)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: X1(3), X2(3), d1(3), d2(3)
    REAL(KIND=dp), INTENT(OUT) :: A(:,:), cpars(:)
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: X3(3), d3(3)
!------------------------------------------------------------------------------
    LOGICAL :: WithThreeNodes
    REAL(KIND=dp) :: d(3), ex(3), ey(3), ez(3), X0(3), v21(3), v31(3), b(3)
    REAL(KIND=dp) :: r1(3), r2(3), r3(3), t1(3), t2(3), t3(3), Norm
!------------------------------------------------------------------------------
    WithThreeNodes = PRESENT(X3) .AND. PRESENT(d3)

    ! ------------------------------------------------------------------------------
    ! The coordinate system is created such that the given vertices lie on 
    ! the plane y=0. A bit different logic is used to decide the orientation of
    ! the coordinate system for the different vertice counts. 
    ! ------------------------------------------------------------------------------
    IF (WithThreeNodes) THEN

      X0(:) = X3(:)

      v21(:) = X2(:) - X1(:)
      Norm = SQRT(SUM(v21(1:3)**2))
      ex(1:3) = 1.0d0/Norm * v21(1:3)

      v31(:) = X3(:) - X1(:)
      Norm = SQRT(SUM(v31(1:3)**2))
      v31(1:3) = 1.0d0/Norm * v31(1:3)
      b(:) = v31(:) - DOT_PRODUCT(v31,ex)*ex(:)
      IF (SQRT(SUM(b(1:3)**2)) < 1.0d-6) THEN
        ! The three vertices are on the same line. Use the mid-node director to generate z-axis.
        b(:) = d3(:) - DOT_PRODUCT(d3,ex)*ex(:)
      END IF
      Norm = SQRT(SUM(b(1:3)**2))
      ez(1:3) = 1.0d0/Norm * b(1:3)
      ey(:) = CrossProduct(ez,ex)

      cpars(1) = DOT_PRODUCT(v21,ex)

      r1(:) = X1(:) - X0(:)
      r2(:) = X2(:) - X0(:)

      ! -----------------------------------------------------------------------------------------
      ! The function BlendingSurfaceInfo is built on the assumption that the mid-node is centered
      ! -----------------------------------------------------------------------------------------
      IF (ABS(DOT_PRODUCT(r1+r2,ex))/cpars(1) > 2.0d-2) THEN
        CALL Warn('EdgeFrame', 'Centered edge mid-nodes expected')
        PRINT *, 'Relative error of node position ...', ABS(DOT_PRODUCT(r1+r2,ex))/cpars(1)
      END IF

      cpars(2) = DOT_PRODUCT(r1,ez)               ! The z-coordinate of the vertex X1
      cpars(3) = DOT_PRODUCT(r2,ez)               ! The z-coordinate of the vertex X2

      t1(:) = CrossProduct(ey,d1)
      t2(:) = CrossProduct(ey,d2)
      t3(:) = CrossProduct(ey,d3)
      cpars(4) = DOT_PRODUCT(t1,ez)/DOT_PRODUCT(t1,ex) ! The angle parameter for z(x) at the vertex X1
      cpars(5) = DOT_PRODUCT(t2,ez)/DOT_PRODUCT(t2,ex) ! The angle parameter for z(x) at the vertex X2
      cpars(6) = DOT_PRODUCT(t3,ez)/DOT_PRODUCT(t3,ex) ! The angle parameter for z(x) at the vertex X3

    ELSE

      d(1:3) = 0.5d0 * ( d1(1:3) + d2(1:3) )
      Norm = SQRT(SUM(d(1:3)**2))
      ez(1:3) = 1.0d0/Norm * d(1:3) 

      v21(:) = X2(:) - X1(:)
      b(:) = v21(:) - DOT_PRODUCT(v21,ez)*ez(:)
      Norm = SQRT(SUM(b(1:3)**2))
      ex(1:3) = 1.0d0/Norm * b(1:3)     
      ey(:) = CrossProduct(ez,ex)

      X0(:) = 0.5d0 * ( X1(1:3) + X2(1:3) )
      cpars(1) = DOT_PRODUCT(v21,ex)
      IF (cpars(1) < 0.0d0) &
          CALL Fatal('EdgeFrame', 'Negative edge length obtained')

      r1(:) = X1(:) - X0(:)
      r2(:) = X2(:) - X0(:)
      cpars(2) = DOT_PRODUCT(r1,ez)               ! The z-coordinate of the vertex X1
      cpars(3) = DOT_PRODUCT(r2,ez)               ! The z-coordinate of the vertex X2
      t1(:) = CrossProduct(ey,d1)
      t2(:) = CrossProduct(ey,d2)
      cpars(4) = DOT_PRODUCT(t1,ez)/DOT_PRODUCT(t1,ex);  ! The angle parameter for z(x) at the vertex X1
      cpars(5) = DOT_PRODUCT(t2,ez)/DOT_PRODUCT(t2,ex);  ! The angle parameter for z(x) at the vertex X2
    END IF

    A(1:3,1) = ex(1:3)
    A(1:3,2) = ey(1:3)
    A(1:3,3) = ez(1:3)
    A(1:3,4) = X0(1:3)
!---------------------------------------------------------------------------
  END SUBROUTINE EdgeFrame
!---------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Compute data which can be used to represent a space curve in the standard 
! Hermite form. The curve tangent vectors expressed with respect to the global 
! frame at the nodes are created by requiring that the tangent vector is orthogonal 
! to the given director vector. The tangent vectors are returned via cpars. The base 
! vectors and the origin X0 of a local coordinate frame are also returned (this 
! may not be of any use in practice) via A. The input data Xk gives the global 
! coordinates of the kth vertex on the edge, while dk specifies the director at 
! the kth vertex on the edge.
!------------------------------------------------------------------------------
  SUBROUTINE HermiteForm(X1, X2, d1, d2, cpars, A, X3, d3)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: X1(3), X2(3), d1(3), d2(3)
    REAL(KIND=dp), INTENT(OUT) :: cpars(:)
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: A(:,:)
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: X3(3), d3(3)
!------------------------------------------------------------------------------
    LOGICAL :: WithThreeNodes
    REAL(KIND=dp) :: d(3), ex(3), ey(3), ez(3), X0(3), v21(3), b(3)
    REAL(KIND=dp) :: t1(3), t2(3), t3(3), Norm, h
!------------------------------------------------------------------------------
    WithThreeNodes = PRESENT(X3) .AND. PRESENT(d3)
    IF (WithThreeNodes) CALL Fatal('HermiteForm', 'Only 2-node version made')

    v21 = X2 - X1
    h = SQRT(SUM(v21**2))
    b = 1.0d0/h * v21

    ! Pick a unit tangent vector orthogonal to the nodal director:
    t1 = v21 - DOT_PRODUCT(v21,d1)*d1
    Norm = SQRT(SUM(t1**2))
    t1 = 1.0d0/Norm * t1
    ! Scale to obtain a tangent vector suitable for the Hermite interpolation
    ! with the derivative DOF being Dv(a_1)[a_2-a_1]:
    Norm = h / DOT_PRODUCT(t1,b)
    t1 = Norm * t1

    ! Repeat for the second node with the derivative DOF Dv(a_2)[a_2-a_1]:
    t2 = v21 - DOT_PRODUCT(v21,d2)*d2
    Norm = SQRT(SUM(t2**2))
    t2 = 1.0d0/Norm * t2
    Norm = h / DOT_PRODUCT(t2,b)
    t2 = Norm * t2

    cpars(1:3) = t1(1:3)
    cpars(4:6) = t2(1:3)

    ! Create a local frame, although this may not be of any use:
    IF (PRESENT(A)) THEN
      X0 = 0.5d0 * (X1 + X2)
      ex = b
      d = 0.5d0 * (d1 + d2) 
      d = d - DOT_PRODUCT(d,ex)*ex
      Norm = SQRT(SUM(d**2))
      ez = 1.0d0/Norm * d 
      ey = CrossProduct(ez,ex)

      A(1:3,1) = ex(1:3)
      A(1:3,2) = ey(1:3)
      A(1:3,3) = ez(1:3)
      A(1:3,4) = X0(1:3)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE HermiteForm
!------------------------------------------------------------------------------

!----------------------------------------------------------------------------
! This function produces the covariant basis {a_i}, the first and second
! fundamental forms A and B, the determinant of the metric surface tensor detA
! and the global coordinates of the point on the blending surface when the 
! reference element coordinates u and v are used as curvilinear coordinates on
! the blending surface. The necessary edge curve data for creating the blending 
! surface must be contained as elementwise properties 'edge frames' and 
! 'edge parameters'. 
! TO DO: Complement and clean the implementation when the initial data
!        is defined over second-order Lagrange elements
!-----------------------------------------------------------------------------  
  FUNCTION BlendingSurfaceInfo( Element, Nodes, u, v, deta, a1, a2, a3, &
      A, B, x, UseMeshOnly ) RESULT(stat)
!----------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Element_t), POINTER, INTENT(IN) :: Element  !< Element structure
    TYPE(Nodes_t), INTENT(IN) :: Nodes               !< Data corresponding to the element nodes
    REAL(KIND=dp), INTENT(IN) :: u                   !< The 1st coordinate of the reference p-element
    REAL(KIND=dp), INTENT(IN) :: v                   !< The 2nd coordinate of the reference p-element
    REAL(KIND=dp), INTENT(OUT) :: deta               !< The determinant of the surface metric tensor
    REAL(KIND=dp), INTENT(OUT) :: a1(3), a2(3)       !< The covariant surface basis vectors
    REAL(KIND=dp), INTENT(OUT) :: a3(3)              !< The base vector normal to the surface
    REAL(KIND=dp), INTENT(OUT) :: A(2,2)             !< The covariant components of the metric surface tensor at (u,v)  
    REAL(KIND=dp), INTENT(OUT) :: B(2,2)             !< The covariant components of the second fundamental form at (u,v)
    REAL(KIND=dp), INTENT(OUT) :: x(3)               !< Blending surface point corresponding to (u,v): x=x(u,v)
    LOGICAL, OPTIONAL, INTENT(IN) :: UseMeshOnly     ! Instead of blending use the original mesh in the calculation
    LOGICAL :: Stat                                  !< A dummy status variable at the moment
!----------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: GElement => NULL()
    LOGICAL :: QuadraticGeometryData, ComputeFromMesh
    INTEGER :: i, j, e, n, q, i0, cn, EdgesParametrized, CurveDataSize, Family
    REAL(KIND=dp), POINTER :: FrameData(:), EdgeParams(:)
    REAL(KIND=dp) :: Basis(MaxPatchNodes), dBasis(MaxPatchNodes,2), ddBasis(4,2,2)
    REAL(KIND=dp) :: xi, eta, ddr(3,2,2)
    REAL(KIND=dp) :: ex(3), ey(3), ez(3), X0(3)
    REAL(KIND=dp) :: d(CurveDataSize2)
    REAL(KIND=dp) :: h, s, t, w, xe, c(3), deltac(3), dc(3), ddc(3)
    REAL(KIND=dp) :: b12(3), db12(3), h12, d1h12, d2h12, ddh12, dsdu, dsdv
    REAL(KIND=dp) :: HermBasis(6), dHermBasis(6), ddHermBasis(6)
    REAL(KIND=dp) :: f, df, ddf
    REAL(KIND=dp) :: h1, h2, h3, dh1, dh2, dh3, ddh1, ddh2, ddh3
    REAL(KIND=dp) :: r1(3), r2(3)
    REAL(KIND=dp) :: d1a1(3), d2a1(3), d2a2(3)
    REAL(KIND=dp) :: L1, L2, dL1, dL2
    REAL(KIND=dp) :: Norm

    SAVE GElement
!----------------------------------------------------------------------------

    Family = GetElementFamily(Element)

    IF (PRESENT(UseMeshOnly)) THEN
      ComputeFromMesh = UseMeshOnly
    ELSE
      ComputeFromMesh = .FALSE.
    END IF

    n = Element % TYPE % NumberOfNodes
    Basis = 0.0d0      
    dBasis = 0.0d0

    IF (ComputeFromMesh) THEN
      !
      ! Compute the fundamental forms simply by using the mesh data and return
      !
      IF (n > MaxPatchNodes) CALL Fatal('BlendingSurfaceInfo', &
          'Nodal basis supported only up to degree 4')

      IF (Family == 3) THEN
        ! Map the p-element coordinates to those of the traditional reference element: 
        xi = 0.5d0 * (1.0d0 + u - 1.0d0/SQRT(3.0d0)*v)
        eta = 1.0d0/SQRT(3.0d0)*v
      ELSE
        xi = u
        eta = v
      END IF

      CALL NodalBasisFunctions2D(Basis, Element, xi, eta)
      CALL NodalFirstDerivatives2D(dBasis, Element, xi, eta)

      ddr(1,1:2,1:2) = SecondDerivatives2D(Element, Nodes % x, xi, eta)
      ddr(2,1:2,1:2) = SecondDerivatives2D(Element, Nodes % y, xi, eta)
      ddr(3,1:2,1:2) = SecondDerivatives2D(Element, Nodes % z, xi, eta)

      X(1) = SUM( Nodes % x(1:n) * Basis(1:n) )
      X(2) = SUM( Nodes % y(1:n) * Basis(1:n) )
      X(3) = SUM( Nodes % z(1:n) * Basis(1:n) )

      a1(1) = SUM( dBasis(1:n,1) * Nodes % x(1:n) )
      a1(2) = SUM( dBasis(1:n,1) * Nodes % y(1:n) )
      a1(3) = SUM( dBasis(1:n,1) * Nodes % z(1:n) )
      a2(1) = SUM( dBasis(1:n,2) * Nodes % x(1:n) )
      a2(2) = SUM( dBasis(1:n,2) * Nodes % y(1:n) )
      a2(3) = SUM( dBasis(1:n,2) * Nodes % z(1:n) )

      d1a1(1:3) = ddr(1:3,1,1)
      d2a2(1:3) = ddr(1:3,2,2)
      d2a1(1:3) = ddr(1:3,1,2)      

      GOTO 101
    END IF

    EdgesParametrized = Family

    SELECT CASE(Family)  
    CASE(3)
      QuadraticGeometryData = Element % TYPE % NumberOfNodes == 6
    CASE(4)
      QuadraticGeometryData = Element % TYPE % NumberOfNodes == 9
      IF (QuadraticGeometryData) EdgesParametrized = 6
    CASE DEFAULT
      CALL Fatal('BlendingSurfaceInfo', 'Only quads and triangles can be handled')     
    END SELECT
    
    IF (QuadraticGeometryData) THEN
      cn = 3 ! The node count per curved edge
      CurveDataSize = CurveDataSize2
    ELSE
      cn = 2
      CurveDataSize = CurveDataSize1
    END IF

    !-----------------------------------------------------------------------
    ! Retrieve parametrizations of curved edges:
    !------------------------------------------------------------------------
    EdgeParams => GetElementProperty('edge parameters', Element)   
    FrameData => NULL()
    IF (QuadraticGeometryData) &
        FrameData => GetElementProperty('edge frames', Element) 

    IF (ASSOCIATED(FrameData)) THEN
      IF (SIZE(FrameData) < EdgesParametrized*FrameDataSize) &
          CALL Fatal('BlendingSurfaceInfo','Frame data are not associated with all edges')
    END IF

    IF (ASSOCIATED(EdgeParams)) THEN
      IF (SIZE(EdgeParams) < EdgesParametrized*CurveDataSize) &
          CALL Fatal('BlendingSurfaceInfo','edge parameters are not associated with all edges')
    ELSE
      CALL Fatal('BlendingSurfaceInfo','edge curves data could not be retrieved')
    END IF

    !---------------------------------------------------------------------------
    n = Element % TYPE % NumberOfNodes
    Basis = 0.0d0      
    dBasis = 0.0d0
    !-------------------------------------------------------------------------
    ! Obtain the lowest-order nodal basis functions on the reference element and
    ! their derivatives with respect to the local coordinates.
    !-------------------------------------------------------------------------
    SELECT CASE(Family)
    CASE(3)
      DO q=1,3
        Basis(q) = TriangleNodalPBasis(q, u, v)
        dBasis(q,1:2) = dTriangleNodalPBasis(q, u, v)
      END DO
    CASE(4)
      DO q=1,4
        Basis(q) = QuadNodalPBasis(q, u, v)
        dBasis(q,1:2) = dQuadNodalPBasis(q, u, v)
      END DO
    END SELECT

    !--------------------------------------------------------------------------
    ! The standard part of the interpolated surface and their contribution 
    ! to the surface basis vectors and to the derivatives of the surface basis 
    ! vectors
    !--------------------------------------------------------------------------
    X(1) = SUM( Nodes % x(1:Family) * Basis(1:Family) )
    X(2) = SUM( Nodes % y(1:Family) * Basis(1:Family) )
    X(3) = SUM( Nodes % z(1:Family) * Basis(1:Family) )

    a1(1) = SUM( dBasis(1:Family,1) * Nodes % x(1:Family) )
    a1(2) = SUM( dBasis(1:Family,1) * Nodes % y(1:Family) )
    a1(3) = SUM( dBasis(1:Family,1) * Nodes % z(1:Family) )
    a2(1) = SUM( dBasis(1:Family,2) * Nodes % x(1:Family) )
    a2(2) = SUM( dBasis(1:Family,2) * Nodes % y(1:Family) )
    a2(3) = SUM( dBasis(1:Family,2) * Nodes % z(1:Family) )

    d1a1 = 0.0d0
    IF (Family == 4) THEN
      d2a1(1) = 1/4.0d0 * (Nodes % x(1) - Nodes % x(2) + Nodes % x(3) - Nodes % x(4))
      d2a1(2) = 1/4.0d0 * (Nodes % y(1) - Nodes % y(2) + Nodes % y(3) - Nodes % y(4))
      d2a1(3) = 1/4.0d0 * (Nodes % z(1) - Nodes % z(2) + Nodes % z(3) - Nodes % z(4))
    ELSE
      d2a1 = 0.0d0
    END IF
    d2a2 = 0.0d0
 
    !a(1,1) = LDot(a1,a1)
    !a(1,2) = LDot(a1,a2)
    !a(2,1) = a(1,2)
    !a(2,2) = LDot(a2,a2)
    !deta = a(1,1)*a(2,2)-a(1,2)*a(2,1)

    !--------------------------------------------------------------------------
    ! Add the blending part edgewise...
    !--------------------------------------------------------------------------
    DO e=1,EdgesParametrized

      i0 = (e-1)*CurveDataSize
      d(1:CurveDataSize) = EdgeParams(i0+1:i0+CurveDataSize)

      IF (QuadraticGeometryData) THEN
        h = d(1) ! The length parameter
        i0 = (e-1)*FrameDataSize
        ex(1:3) = FrameData(i0+FrameBasis1)
        ey(1:3) = FrameData(i0+FrameBasis2)
        ez(1:3) = FrameData(i0+FrameBasis3)
        X0(1:3) = FrameData(i0+FrameOrigin)
      ELSE
        h = 2.0d0 ! the length of reference element [-1,1]
      END IF
        
      SELECT CASE(Family)
      CASE(3)
        !-------------------------------------------------------------------------
        ! First define edge curve parameter and a part of the blending function
        !-------------------------------------------------------------------------
        SELECT CASE(e)
        CASE(1)
          s = Basis(2) - Basis(1)
          h1 = Basis(2) * Basis(1)
          r1(1) = Nodes % x(1)
          r1(2) = Nodes % y(1)
          r1(3) = Nodes % z(1)
          r2(1) = Nodes % x(2)
          r2(2) = Nodes % y(2)
          r2(3) = Nodes % z(2)
        CASE(2)
          s = Basis(3) - Basis(2)
          h1 = Basis(3) * Basis(2)
          r1(1) = Nodes % x(2)
          r1(2) = Nodes % y(2)
          r1(3) = Nodes % z(2)
          r2(1) = Nodes % x(3)
          r2(2) = Nodes % y(3)
          r2(3) = Nodes % z(3)
        CASE(3)
          s = Basis(1) - Basis(3)
          h1 = Basis(1) * Basis(3)
          r1(1) = Nodes % x(3)
          r1(2) = Nodes % y(3)
          r1(3) = Nodes % z(3)
          r2(1) = Nodes % x(1)
          r2(2) = Nodes % y(1)
          r2(3) = Nodes % z(1)
        END SELECT

        !-------------------------------------------------------------------------
        ! The basis functions for the edge curve expansion in terms of s and
        ! their derivatives
        !-------------------------------------------------------------------------
        CALL HermiteBasis(s, h, HermBasis(1:2*cn), dHermBasis(1:2*cn), ddHermBasis(1:2*cn), cn)
        L1 = 0.5d0 * (1.0d0 - s)
        L2 = 0.5d0 * (1.0d0 + s)
        dL1 = -0.5d0
        dL2 = 0.5d0

        ! A scalar-valued blending function:
        h12 = h1/(L1*L2)

        !----------------------------------------------------------------------------
        ! Decompose the edge curve into components along the global coordinate axes
        ! and compute the derivatives of the edge curve with respect to s
        !----------------------------------------------------------------------------
        IF (QuadraticGeometryData) THEN
          !-------------------------------------------------------------------------
          ! The dimensional edge parameterization coordinate corresponding to
          ! the dimensionless coordinate s
          !--------------------------------------------------------------------------
          xe = -0.25d0*(1.0d0-s)*h + 0.25d0*(1.0d0+s)*h

          f = d(2)*HermBasis(1) + d(3)*HermBasis(2) + SUM(d(4:6) * HermBasis(4:6))
          df = d(2)*dHermBasis(1) + d(3)*dHermBasis(2) + SUM(d(4:6) * dHermBasis(4:6))
          ddf = d(2)*ddHermBasis(1) + d(3)*ddHermBasis(2) + SUM(d(4:6) * ddHermBasis(4:6))

          c(1:3) = X0(1:3) + xe * ex(1:3) + f * ez(1:3)
          dc(1:3) = 0.5d0*h * ex(1:3) + df * ez(1:3)
          ddc(1:3) = ddf * ez(1:3)
        ELSE
          c(1:3) = r1(1:3)*HermBasis(1) + r2(1:3)*HermBasis(2) + &
              d(1:3)*0.5d0*HermBasis(3) + d(4:6)*0.5d0*HermBasis(4)
          dc(1:3) = r1(1:3)*dHermBasis(1) + r2(1:3)*dHermBasis(2) + &
              d(1:3)*0.5d0*dHermBasis(3) + d(4:6)*0.5d0*dHermBasis(4)
          ddc(1:3) = r1(1:3)*ddHermBasis(1) + r2(1:3)*ddHermBasis(2) + &
              d(1:3)*0.5d0*ddHermBasis(3) + d(4:6)*0.5d0*ddHermBasis(4)
        END IF

        !---------------------------------------------------------------------------
        ! The contributions of the blending functions to the position vector, to the
        ! surface base vectors and to the derivatives of the surface base vectors
        !---------------------------------------------------------------------------
        SELECT CASE(e)
        CASE(1)
          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          X(1:3) = X(1:3) + h12 * b12(1:3)

          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)    ! derivative wrt s
          ! (d/du)h12:
          d1h12 = 0.2D1 / 0.9D1 * sqrt(0.3D1) * u * v * (sqrt(0.3D1) * v - 0.6D1) / &
              ( (-1.0D0 + u)**2 * (1.0D0 + u)**2 )
          ! (d/dv)h12:
          d2h12 = -0.2D1 / 0.9D1 * sqrt(0.3D1) * (sqrt(0.3D1) * v - 0.3D1) / (u**2 - 1.0D0)

          a1(1:3) = a1(1:3) + h12 * db12(1:3) + d1h12 * b12(1:3)
          a2(1:3) = a2(1:3) + d2h12 * b12(1:3)

          !d^2/(du)^2 h12:
          ddh12 = -0.2D1 / 0.9D1 * sqrt(0.3D1) * v * (sqrt(0.3D1) * v - 0.6D1) * &
              (3 * u**2 + 1.0D0) / ( (-1.0D0 + u)**3 * (1.0D0 + u)**3 )
          d1a1(1:3) = d1a1(1:3) + h12 * ddc(1:3) + 2.0d0 * d1h12 * db12(1:3) + ddh12 * b12(1:3) 

          !d^2/(dv)^2 h12:
          ddh12 = -0.2D1 / ( 0.3D1 * (u**2 - 1.0D0) )
          d2a2(1:3) = d2a2(1:3) + ddh12 * b12(1:3)
          !d^2/(dudv) h12:
          ddh12 = 0.4D1 / 0.9D1 * sqrt(0.3D1) * (sqrt(0.3D1) * v - 0.3D1) * u / &
              ( (-1.0D0 + u)**2 * (1.0D0 + u)**2 )
          d2a1(1:3) = d2a1(1:3) + d2h12 * db12(1:3) + ddh12 * b12(1:3)
   
        CASE(2)
          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          X(1:3) = X(1:3) + h12 * b12(1:3)

          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)
          ! (d/du)h12:
          d1h12 = -0.8D1 / 0.9D1 * sqrt(0.3D1) * v * ( 3.0D0 * v ** 2 + &
              0.2D1 * sqrt(0.3D1) * u * v + 0.2D1 * sqrt(0.3D1) * v - 0.3D1 * u ** 2 - &
              0.6D1 * u - 0.15D2) / ( (sqrt(0.3D1) * v - u - 0.3D1) ** 2 * &
              (sqrt(0.3D1) * v - u + 0.1D1) ** 2 )
          !(d/dv)h12:
          d2h12 = -0.8D1 / 0.3D1 * (sqrt(0.3D1) * u ** 3 - sqrt(0.3D1) * u * v** 2 + &
              0.3D1 * sqrt(0.3D1) * u ** 2 - sqrt(0.3D1) * v ** 2 - 0.2D1 * u ** 2 * v - &
              sqrt(0.3D1) * u - 0.4D1 * u * v - 0.3D1 * sqrt(0.3D1) + 0.6D1 * v) / &
              ( (sqrt(0.3D1) * v - u - 0.3D1) ** 2 * (sqrt(0.3D1) * v - u + 0.1D1) ** 2 )
          dsdu = -0.5d0
          dsdv = sqrt(3.0d0)/2.0d0

          a1(1:3) = a1(1:3) + h12 * db12(1:3) * dsdu + d1h12 * b12(1:3)
          a2(1:3) = a2(1:3) + h12 * db12(1:3) * dsdv + d2h12 * b12(1:3)

          !d^2/(du)^2 h12:
          ddh12 = 0.16D2 / 0.9D1 * sqrt(0.3D1) * v * (0.3D1 * sqrt(0.3D1) * u ** 2 * v - &
              0.9D1 * sqrt(0.3D1) * v ** 3 + 0.6D1 * sqrt(0.3D1) * u * v - 0.3D1 * u ** 3 + &
              0.9D1 * u * v ** 2 + 0.31D2 * sqrt(0.3D1) * v - 0.9D1 * u ** 2 + 0.9D1 * v ** 2 - &
              0.45D2 * u - 0.39D2) / ( (sqrt(0.3D1) * v - u - 0.3D1) ** 3 * &
              (sqrt(0.3D1) * v - u + 0.1D1) ** 3 )

          d1a1(1:3) = d1a1(1:3) + h12 * ddc(1:3) * dsdu**2 + ddh12 * b12(1:3) + &
              2.0d0 * d1h12 * db12(1:3) * dsdu

          !d^2/(dv)^2 h12:
          ddh12 = 0.16D2 / 0.3D1 * (0.9D1 * sqrt(0.3D1) * u ** 3 * v - 0.3D1 * sqrt(0.3D1) * u * &
              v ** 3 + 0.27D2 * sqrt(0.3D1) * u ** 2 * v - 0.3D1 * sqrt(0.3D1) * v ** 3 - &
              0.5D1 * u ** 4 - 0.9D1 * u ** 2 * v ** 2 - 0.9D1 * sqrt(0.3D1) * u * v - 0.20D2 * u ** 3 - &
              0.18D2 * u * v ** 2 - 0.27D2 * sqrt(0.3D1) * v - 0.14D2 * u ** 2 + 0.27D2 * v ** 2 + &
              0.12D2 * u + 0.27D2) / ( (sqrt(0.3D1) * v - u - 0.3D1) ** 3 * &
              (sqrt(0.3D1) * v - u + 0.1D1) ** 3 )

          d2a2(1:3) = d2a2(1:3) + h12 * ddc(1:3) * dsdv**2 + ddh12 * b12(1:3) + &
              2.0d0 * d2h12 * db12(1:3) * dsdv

          !d^2/(dudv) h12:
          ddh12 = 0.8D1 / 0.3D1 * (sqrt(0.3D1) * u ** 4 - 0.12D2 * sqrt(0.3D1) * u ** 2 * v ** 2 + &
              0.3D1 * sqrt(0.3D1) * v ** 4 + 0.4D1 * sqrt(0.3D1) * u ** 3 - 0.24D2 * sqrt(0.3D1) * &
              u * v ** 2 + 0.2D1 * u ** 3 * v + 0.18D2 * u * v ** 3 + 0.6D1 * sqrt(0.3D1) * u ** 2 - &
              0.36D2 * sqrt(0.3D1) * v ** 2 + 0.6D1 * u ** 2 * v + 0.18D2 * v ** 3 + &
              0.4D1 * sqrt(0.3D1) * u + 0.46D2 * u * v - 0.15D2 * sqrt(0.3D1) + 0.42D2 * v) / &
              ( ( sqrt(0.3D1) * v - u - 0.3D1) ** 3 * (sqrt(0.3D1) * v - u + 0.1D1) ** 3 )

          d2a1(1:3) = d2a1(1:3) + (d1h12 * dsdv + d2h12 * dsdu) * db12(1:3) + ddh12 * b12(1:3) + &
              h12 * ddc(1:3) * dsdv * dsdu

        CASE(3)
          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          X(1:3) = X(1:3) + h12 * b12(1:3)

          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)
          ! (d/du)h12:
          d1h12 = -0.8D1 / 0.9D1 * sqrt(0.3D1) * v * (0.2D1 * sqrt(0.3D1) * u * v - &
              0.2D1 * sqrt(0.3D1) * v + 0.3D1 * u ** 2 - 0.3D1 * v ** 2 - 0.6D1 * u + 0.15D2) / &
              ( (sqrt(0.3D1) * v + u + 0.1D1) ** 2 * (sqrt(0.3D1) * v + u - 0.3D1) ** 2 )
          !(d/dv)h12:
          d2h12 = 0.8D1 / 0.3D1 * (sqrt(0.3D1) * u ** 3 - sqrt(0.3D1) * u * v ** 2 - &
              0.3D1 * sqrt(0.3D1) * u ** 2 + sqrt(0.3D1) * v ** 2 + 0.2D1 * u ** 2 * v - &
              sqrt(0.3D1) * u - 0.4D1 * u * v + 0.3D1 * sqrt(0.3D1) - 0.6D1 * v) / &
              ( (sqrt(0.3D1) * v + u + 0.1D1) ** 2 * (sqrt(0.3D1) * v + u - 0.3D1) ** 2 )          
          dsdu = -0.5d0
          dsdv = -sqrt(3.0d0)/2.0d0

          a1(1:3) = a1(1:3) + h12 * db12(1:3) * dsdu + d1h12 * b12(1:3)
          a2(1:3) = a2(1:3) + h12 * db12(1:3) * dsdv + d2h12 * b12(1:3)

          !d^2/(du)^2 h12:
          ddh12 = 0.16D2 / 0.9D1 * sqrt(0.3D1) * v * (0.3D1 * sqrt(0.3D1) * u ** 2 * v - &
              0.9D1 * sqrt(0.3D1) * v ** 3 - 0.6D1 * sqrt(0.3D1) * u * v + 0.3D1 * u ** 3 - &
              0.9D1 * u * v ** 2 + 0.31D2 * sqrt(0.3D1) * v - 0.9D1 * u ** 2 + 0.9D1 * v ** 2 + &
              0.45D2 * u - 0.39D2) / ( (sqrt(0.3D1) * v + u + 0.1D1) ** 3 * &
              (sqrt(0.3D1) * v + u - 0.3D1) ** 3 )

          d1a1(1:3) = d1a1(1:3) + h12 * ddc(1:3) * dsdu**2 + ddh12 * b12(1:3) + &
              2.0d0 * d1h12 * db12(1:3) * dsdu

          !d^2/(dv)^2 h12:
          ddh12 = -0.16D2 / 0.3D1 * (0.9D1 * sqrt(0.3D1) * u ** 3 * v - 0.3D1 * sqrt(0.3D1) * u * v ** 3 - &
              0.27D2 * sqrt(0.3D1) * u ** 2 * v + 0.3D1 * sqrt(0.3D1) * v ** 3 + 0.5D1 * u ** 4 + &
              0.9D1 * u ** 2 * v ** 2 - 0.9D1 * sqrt(0.3D1) * u * v - 0.20D2 * u ** 3 - &
              0.18D2 * u * v ** 2 + 0.27D2 * sqrt(0.3D1) * v + 0.14D2 * u ** 2 - 0.27D2 * v ** 2 + &
              0.12D2 * u - 0.27D2) / ( (sqrt(0.3D1) * v + u + 0.1D1) ** 3 * &
              (sqrt(0.3D1) * v + u - 0.3D1) ** 3 )

          d2a2(1:3) = d2a2(1:3) + h12 * ddc(1:3) * dsdv**2 + ddh12 * b12(1:3) + &
              2.0d0 * d2h12 * db12(1:3) * dsdv

          !d^2/(dudv) h12:
          ddh12 = -0.8D1 / 0.3D1 * (sqrt(0.3D1) * u ** 4 - 0.12D2 * sqrt(0.3D1) * u ** 2 * v ** 2 + &
              0.3D1 * sqrt(0.3D1) * v ** 4 - 0.4D1 * sqrt(0.3D1) * u ** 3 + 0.24D2 * sqrt(0.3D1) * &
              u * v ** 2 - 0.2D1 * u ** 3 * v - 0.18D2 * u * v ** 3 + 0.6D1 * sqrt(0.3D1) * u ** 2 - &
              0.36D2 * sqrt(0.3D1) * v ** 2 + 0.6D1 * u ** 2 * v + 0.18D2 * v ** 3 - 0.4D1 * &
              sqrt(0.3D1) * u - 0.46D2 * u * v - 0.15D2 * sqrt(0.3D1) + 0.42D2 * v) / &
              ( (sqrt(0.3D1) * v + u + 0.1D1) ** 3 * (sqrt(0.3D1) * v + u - 0.3D1) ** 3 )

          d2a1(1:3) = d2a1(1:3) + (d1h12 * dsdv + d2h12 * dsdu) * db12(1:3) + ddh12 * b12(1:3) + &
              h12 * ddc(1:3) * dsdv * dsdu

        END SELECT
        
      CASE(4)
        !-------------------------------------------------------------------------
        ! First define edge orientation convention and retrieve parameters for
        ! representing the curved edge
        !-------------------------------------------------------------------------
        SELECT CASE(e)
        CASE(1)
          r1(1) = Nodes % x(1)
          r1(2) = Nodes % y(1)
          r1(3) = Nodes % z(1)
          r2(1) = Nodes % x(2)
          r2(2) = Nodes % y(2)
          r2(3) = Nodes % z(2)
          s = u
          t = v
        CASE(2)
          r1(1) = Nodes % x(2)
          r1(2) = Nodes % y(2)
          r1(3) = Nodes % z(2)
          r2(1) = Nodes % x(3)
          r2(2) = Nodes % y(3)
          r2(3) = Nodes % z(3)
          s = v
          t = u
        CASE(3)
          r1(1) = Nodes % x(4)
          r1(2) = Nodes % y(4)
          r1(3) = Nodes % z(4)
          r2(1) = Nodes % x(3)
          r2(2) = Nodes % y(3)
          r2(3) = Nodes % z(3)
          s = u
          t = v
        CASE(4)
          r1(1) = Nodes % x(1)
          r1(2) = Nodes % y(1)
          r1(3) = Nodes % z(1)
          r2(1) = Nodes % x(4)
          r2(2) = Nodes % y(4)
          r2(3) = Nodes % z(4)
          s = v
          t = u
        CASE(5)
          s = u
          t = v
        CASE(6)
          s = v
          t = u
        END SELECT

        !-------------------------------------------------------------------------
        ! The basis functions for the edge curve expansion in terms of s and
        ! their derivatives
        !-------------------------------------------------------------------------
        CALL HermiteBasis(s, h, HermBasis(1:2*cn), dHermBasis(1:2*cn), ddHermBasis(1:2*cn), cn)

        ! ----------------------------------------------------------------------------
        ! The Hermite basis functions in terms of variable t and their derivatives.
        ! These are used to define how the blending edge basis function decays towards
        ! the opposite edge.
        !-----------------------------------------------------------------------------
        IF (QuadraticGeometryData) THEN
          h1 = 0.5d0 * t * (t - 1.0d0)
          h2 = 0.5d0 * t * (t + 1.0d0)
          h3 = 1.0d0 - t**2
          dh1 = t - 0.5d0
          dh2 = t + 0.5d0
          dh3 = -2.0d0 * t
          ddh1 = 1.0d0
          ddh2 = 1.0d0
          ddh3 = -2.0d0
        ELSE
          ! The following defines linear decay
          h1 = 0.5d0 * (1.0d0 - t)
          h2 = 0.5d0 * (1.0d0 + t)
          dh1 = -0.5d0
          dh2 = 0.5d0
          ddh1 = 0.0d0
          ddh2 = 0.0d0
        END IF

        !------------------------------------------------------------
        ! The standard linear basis functions and their derivatives
        ! on an edge [-1,1]
        !-------------------------------------------------------         
        L1 = 0.5d0 * (1.0d0 - s)
        L2 = 0.5d0 * (1.0d0 + s)
        dL1 = -0.5d0
        dL2 = 0.5d0

        !----------------------------------------------------------------------------
        ! Decompose the edge curve into components along the global coordinate axes
        ! and compute the derivatives of the edge curve with respect to s
        !----------------------------------------------------------------------------
        IF (QuadraticGeometryData) THEN
          !-------------------------------------------------------------------------
          ! The dimensional edge parameterization coordinate corresponding to
          ! the dimensionless coordinate s (depending on the edge either s=u or s=v)
          !--------------------------------------------------------------------------
          xe = -0.25d0*(1.0d0-s)*h + 0.25d0*(1.0d0+s)*h

          f = d(2)*HermBasis(1) + d(3)*HermBasis(2) + SUM(d(4:6) * HermBasis(4:6))
          df = d(2)*dHermBasis(1) + d(3)*dHermBasis(2) + SUM(d(4:6) * dHermBasis(4:6))
          ddf = d(2)*ddHermBasis(1) + d(3)*ddHermBasis(2) + SUM(d(4:6) * ddHermBasis(4:6))

          c(1:3) = X0(1:3) + xe * ex(1:3) + f * ez(1:3)
          dc(1:3) = 0.5d0*h * ex(1:3) + df * ez(1:3)
          ddc(1:3) = ddf * ez(1:3)
        ELSE
          c(1:3) = r1(1:3)*HermBasis(1) + r2(1:3)*HermBasis(2) + &
              d(1:3)*0.5d0*HermBasis(3) + d(4:6)*0.5d0*HermBasis(4)
          dc(1:3) = r1(1:3)*dHermBasis(1) + r2(1:3)*dHermBasis(2) + &
              d(1:3)*0.5d0*dHermBasis(3) + d(4:6)*0.5d0*dHermBasis(4)
          ddc(1:3) = r1(1:3)*ddHermBasis(1) + r2(1:3)*ddHermBasis(2) + &
              d(1:3)*0.5d0*ddHermBasis(3) + d(4:6)*0.5d0*ddHermBasis(4)          
        END IF

        !---------------------------------------------------------------------------
        ! The contributions of the blending functions to the position vector, to the
        ! surface base vectors and to the derivatives of the surface base vectors
        !---------------------------------------------------------------------------
        SELECT CASE(e)
        CASE(1)
          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)
          X(1:3) = X(1:3) + h1 * b12(1:3)

          a1(1:3) = a1(1:3) + h1 * db12(1:3)
          a2(1:3) = a2(1:3) + dh1 * b12(1:3)
          d1a1(1:3) = d1a1(1:3) + h1 * ddc(1:3)
          d2a1(1:3) = d2a1(1:3) + dh1 * db12(1:3)
          d2a2(1:3) = d2a2(1:3) + ddh1 * b12(1:3)

        CASE(2)
          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)
          X(1:3) = X(1:3) + h2 * b12(1:3)

          a1(1:3) = a1(1:3) + dh2 * b12(1:3)
          a2(1:3) = a2(1:3) + h2 * db12(1:3)
          d1a1(1:3) = d1a1(1:3) + ddh2 * b12(1:3)
          d2a1(1:3) = d2a1(1:3) + dh2 * db12(1:3)
          d2a2(1:3) = d2a2(1:3) + h2 * ddc(1:3)

        CASE(3)
          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)
          X(1:3) = X(1:3) + h2 * b12(1:3)

          a1(1:3) = a1(1:3) + h2 * db12(1:3)
          a2(1:3) = a2(1:3) + dh2 * b12(1:3)
          d1a1(1:3) = d1a1(1:3) + h2 * ddc(1:3)
          d2a1(1:3) = d2a1(1:3) + dh2 * db12(1:3)
          d2a2(1:3) = d2a2(1:3) + ddh2 * b12(1:3)

        CASE(4)
          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)
          X(1:3) = X(1:3) + h1 * b12(1:3)

          a1(1:3) = a1(1:3) + dh1 * b12(1:3)
          a2(1:3) = a2(1:3) + h1 * db12(1:3)
          d1a1(1:3) = d1a1(1:3) + ddh1 * b12(1:3)
          d2a1(1:3) = d2a1(1:3) + dh1 * db12(1:3)
          d2a2(1:3) = d2a2(1:3) + h1 * ddc(1:3)

        CASE(5)
          r1(1) = Nodes % x(8)
          r1(2) = Nodes % y(8)
          r1(3) = Nodes % z(8)
          r2(1) = Nodes % x(6)
          r2(2) = Nodes % y(6)
          r2(3) = Nodes % z(6)

          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)
          X(1:3) = X(1:3) + h3 * b12(1:3)

          a1(1:3) = a1(1:3) + h3 * db12(1:3)
          a2(1:3) = a2(1:3) + dh3 * b12(1:3)
          d1a1(1:3) = d1a1(1:3) + h3 * ddc(1:3)
          d2a1(1:3) = d2a1(1:3) + dh3 * db12(1:3)
          d2a2(1:3) = d2a2(1:3) + ddh3 * b12(1:3)

        CASE(6)
          r1(1) = Nodes % x(5)
          r1(2) = Nodes % y(5)
          r1(3) = Nodes % z(5)
          r2(1) = Nodes % x(7)
          r2(2) = Nodes % y(7)
          r2(3) = Nodes % z(7)

          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)
          X(1:3) = X(1:3) + h3 * b12(1:3)

          a1(1:3) = a1(1:3) + dh3 * b12(1:3)
          a2(1:3) = a2(1:3) + h3 * db12(1:3)
          d1a1(1:3) = d1a1(1:3) + ddh3 * b12(1:3)
          d2a1(1:3) = d2a1(1:3) + dh3 * db12(1:3)
          d2a2(1:3) = d2a2(1:3) + h3 * ddc(1:3)

        END SELECT

      END SELECT
    END DO

    101 CONTINUE
    !--------------------------------------------------------------------
    ! The metric surface tensor and its determinant
    !--------------------------------------------------------------------
    a(1,1) = DOT_PRODUCT(a1,a1)
    a(1,2) = DOT_PRODUCT(a1,a2)
    a(2,1) = a(1,2)
    a(2,2) = DOT_PRODUCT(a2,a2)
    deta = a(1,1)*a(2,2)-a(1,2)*a(2,1)

    !--------------------------------------------------------------------
    ! The covariant components of the curvature tensor. Note that the 
    ! direction of the normal vector computed here depends on the node
    ! numbering.
    !--------------------------------------------------------------------
    a3(:) = CrossProduct(a1,a2)
    Norm = SQRT(SUM(a3(1:3)**2))
    a3(:) = 1/Norm * a3(:) 
    b(1,1) = DOT_PRODUCT(a3,d1a1)
    b(1,2) = DOT_PRODUCT(a3,d2a1)
    b(2,1) = b(1,2)
    b(2,2) = DOT_PRODUCT(a3,d2a2)

    stat = .TRUE.
!----------------------------------------------------------------------------
  END FUNCTION BlendingSurfaceInfo
!-----------------------------------------------------------------------------

!------------------------------------------------------------------------------
! This subroutine gives the referential description B(f(p)) of the Hermite basis 
! functions B(x), with f being the mapping of the reference element [-1,1] to
! the physical element. The length of the physical element h is used as a scale
! factor so that an interpolating function w(x) has w(x_i) and Dw(x_i)[u] as  
! the nodal DOFs. In practice, if b(p) were constructed to be a basis function 
! on the reference element and were associated with the derivative DOF, we would 
! have B(f(p)) = h/2 b(p) (for the other basis functions B(f(p)) = b(p) as usual). 
! This subroutine returns also the first and second derivatives d/dp B(f(p)) 
! and d^2/dp^2 B(f(p)) or dB/dx(f(p)) and d^2B/dx^2(f(p)).
!------------------------------------------------------------------------------
  SUBROUTINE HermiteBasis(s, h, Basis, dBasis, ddBasis, n, GlobalDerivative)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: s           ! Coordinate
    REAL(KIND=dp), INTENT(IN) :: h           ! A scale factor
    REAL(KIND=dp), INTENT(OUT) :: Basis(:)   ! Basis functions
    REAL(KIND=dp), INTENT(OUT) :: dBasis(:)  ! The first-order derivatives
    REAL(KIND=dp), INTENT(OUT) :: ddBasis(:) ! The second-order derivatives
    INTEGER, OPTIONAL, INTENT(IN) :: n       ! The number of nodes
    LOGICAL, OPTIONAL, INTENT(IN) :: GlobalDerivative  ! To return dB/dx and d^2B/dx^2 
!------------------------------------------------------------------------------
    INTEGER :: NodeCount, DOFs
    LOGICAL :: TransformDerivatives
!------------------------------------------------------------------------------
    IF (PRESENT(n)) THEN
      NodeCount = n
    ELSE
      NodeCount = 2
    END IF
    DOFs = 2*NodeCount

    IF (PRESENT(GlobalDerivative)) THEN
      TransformDerivatives = GlobalDerivative
    ELSE
      TransformDerivatives = .FALSE.
    END IF

    IF (SIZE(Basis) < DOFs .OR. SIZE(dBasis) < DOFs .OR.  SIZE(ddBasis) < DOFs) &
        CALL Fatal('HermiteBasis', 'Too small arrays for basis functions')

    SELECT CASE(NodeCount)
    CASE(2)
      ! ------------------------------------------------------------------------
      ! Third-order Hermite basis: Dofs are listed as w(v1), w(v2),
      ! Dw(v1)[u], Dw(v2)[u]
      ! ------------------------------------------------------------------------
      Basis(1) = (2.0d0 - 3*s + s**3)/4.0d0
      Basis(2) = (2.0d0 + 3*s - s**3)/4.0d0
      Basis(3) = h/8.0d0 - (h*s)/8.0d0 - (h*s**2)/8.0d0 + (h*s**3)/8.0d0
      Basis(4) = -h/8.0d0 - (h*s)/8.0d0 + (h*s**2)/8.0d0 + (h*s**3)/8.0d0

      dBasis(1) = (-3.0d0 + 3*s**2)/4.0d0
      dBasis(2) = (3.0d0 - 3*s**2)/4.0d0
      dBasis(3) = -h/8.0d0 - (h*s)/4.0d0 + (3*h*s**2)/8.0d0
      dBasis(4) = -h/8.0d0 + (h*s)/4.0d0 + (3*h*s**2)/8.0d0

      ddBasis(1) = (3*s)/2.0d0
      ddBasis(2) = (-3*s)/2.0d0
      ddBasis(3) = (h*(-1 + 3*s))/4.0d0
      ddBasis(4) = (h + 3*h*s)/4.0d0
    CASE(3)
      ! ------------------------------------------------------------------------
      ! Third-order Hermite basis: Dofs are listed as w(v1), w(v2), w(v3),
      ! Dw(v1)[u], Dw(v2)[u], Dw(v3)[u]
      ! ------------------------------------------------------------------------
      Basis(1) = s ** 2 - 0.5D1 / 0.4D1 * s ** 3 - s ** 4 / 0.2D1 + 0.3D1 / 0.4D1 * s ** 5
      Basis(2) = s ** 2 + 0.5D1 / 0.4D1 * s ** 3 - s ** 4 / 0.2D1 - 0.3D1 / 0.4D1 * s ** 5
      Basis(3) = s ** 4 - 2.0D0 * s ** 2 + 1.0D0
      Basis(4) = h/2.0d0 * (s ** 2 / 0.4D1 - s ** 3 / 0.4D1 - s ** 4 / 0.4D1 + s ** 5 / 0.4D1)
      Basis(5) = h/2.0d0 * (-s ** 2 / 0.4D1 - s ** 3 / 0.4D1 + s ** 4 / 0.4D1 + s ** 5 /0.4D1)
      Basis(6) = h/2.0d0 * (s ** 5 - 2.0D0 * s ** 3 + s)

      dBasis(1) = 2.0D0 * s - 0.15D2 / 0.4D1 * s ** 2 - 2.0D0 * s ** 3 + 0.15D2 / 0.4D1 * s ** 4
      dBasis(2) = 2 * s + 0.15D2 / 0.4D1 * s ** 2 - 2.0D0 * s ** 3 - 0.15D2 / 0.4D1 * s ** 4
      dBasis(3) = 4.0D0 * s ** 3 - 4.0D0 * s
      dBasis(4) = h/2.0d0 * (s / 0.2D1 - 0.3D1 / 0.4D1 * s ** 2 - s ** 3 + 0.5D1 / 0.4D1 * s ** 4)
      dBasis(5) = h/2.0d0 * (-s / 0.2D1 - 0.3D1 / 0.4D1 * s ** 2 + s ** 3 + 0.5D1 / 0.4D1* s ** 4)
      dBasis(6) = h/2.0d0 * (5.0D0 * s ** 4 - 6.0D0 * s ** 2 + 1.0D0)

      ddBasis(1) = 0.2D1 - 0.15D2 / 0.2D1 * s - 0.6D1 * s ** 2 + 0.15D2 * s ** 3
      ddBasis(2) = 0.2D1 + 0.15D2 / 0.2D1 * s - 0.6D1 * s ** 2 - 0.15D2 * s ** 3
      ddBasis(3) = 0.12D2 * s ** 2 - 0.4D1
      ddBasis(4) = h/2.0d0 * (0.1D1 / 0.2D1 - 0.3D1 / 0.2D1 * s - 0.3D1 * s ** 2 + 0.5D1 * s ** 3)
      ddBasis(5) = h/2.0d0 * (-0.1D1 / 0.2D1 - 0.3D1 / 0.2D1 * s + 0.3D1 * s ** 2 + 0.5D1 * s ** 3)
      ddBasis(6) = h/2.0d0 * (0.2D2 * s ** 3 - 0.12D2 * s)
    CASE DEFAULT
      CALL Fatal('HermiteBasis', 'An unsupported element type')     
    END SELECT

    IF (TransformDerivatives) THEN
      dBasis(1:DOFs) = 2.0d0/h * dBasis(1:DOFs)
      ddBasis(1:DOFs) = 4.0d0/h**2 * ddBasis(1:DOFs)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE HermiteBasis
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! This subroutine creates an orthonormal basis which gives the orientation of 
! lines of curvature at a point that is the image of the reference p-element 
! point (xi1,xi2) under the FE blending map. If the point is not specified,
! the point is taken to be the element centre. With ReparametrizeMesh=.true.
! the blending is not applied and all computations are performed by using
! the original mesh (a mesh consisting of higher-order nodal elements should
! then be used, since second-order partial derivatives must be evaluated). 
!   The basis vectors may be returned via the arguments e1, e2 and e3, 
! while the coordinates of the surface point may be returned via o.
! When this subroutine is called at the element centre, it can be used to 
! define an elementwise coordinate system which is aligned with lines of curvature. 
! In this case the optional arguments LagrangeNodes and TaylorApproximation 
! can be used to obtain data for the Lagrange interpolation to describe the shape 
! of the projected surface on the plane spanned by {e1,e2} and the coefficients of 
! the third-order Taylor polynomial which approximates the blending surface. 
! The optional argument d can be used to ensure that e3 and d point to the same
! direction, while PlanarSurface indicates whether the surface is planar
! at the given point.
!   The argument SaveProperties can be used to save the quantities computed as 
! elementwise properties to avoid later recomputation. The elementwise properties
! are as follows.
!    * 'element frame': e1, e2, e3 and o
!    * 'taylor parameters': the coefficients of the Taylor polynomial
!    * 'planar point': to indicate that the surface is planar
!    * 'umbilical point': to indicate that the surface is umbilical
!------------------------------------------------------------------------------
  SUBROUTINE LinesOfCurvatureFrame(Element, xi1, xi2, e1, e2, e3, o, TaylorApproximation, &
      LagrangeNodes, d, PlanarSurface, PlanarPoint, UmbilicalPoint, &
      SaveProperties, SizeRadiusRatio, ReparametrizeMesh)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: xi1, xi2
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: e1(3), e2(3), e3(3)
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: o(3)
    LOGICAL, OPTIONAL, INTENT(IN) :: TaylorApproximation
    REAL(KIND=dp), OPTIONAL, TARGET, INTENT(OUT) :: LagrangeNodes(MaxPatchNodes,3)
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: d(3)
    LOGICAL, OPTIONAL, INTENT(IN) :: PlanarSurface
    LOGICAL, OPTIONAL, INTENT(OUT) :: PlanarPoint
    LOGICAL, OPTIONAL, INTENT(OUT) :: UmbilicalPoint
    LOGICAL, OPTIONAL, INTENT(IN) :: SaveProperties  
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: SizeRadiusRatio
    LOGICAL, OPTIONAL, INTENT(IN) :: ReparametrizeMesh
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: GElement => NULL()
    LOGICAL :: Stat, Found
    LOGICAL :: Converged, ComputeTaylorPolynomial
    LOGICAL :: ApproximatePlaneDomain, CheckOrientation, WriteElementProperties
    LOGICAL :: Planar, Umbilical
    LOGICAL :: UseMeshOnly
    INTEGER :: Family, n, m, e, i, j, k, GridPoint

    REAL(KIND=dp) :: TaylorParams(6), PlanarFlag(1), UmbilicalFlag(1) 
    REAL(KIND=dp) :: u, v
    REAL(KIND=dp) :: GlobPDir1(3), GlobPDir2(3), GlobPDir3(3), X0(3)
    REAL(KIND=dp) :: Lambda1, Lambda2, LambdaMax
    REAL(KIND=dp) :: a1(3), a2(3), a3(3), a(2,2), Deta, b(2,2), ContravA(2,2)
    REAL(KIND=dp) :: c(2,2), trc, detc, discriminant, PDir1(2), PDir2(2)
    REAL(KIND=dp) :: DualBase1(3), DualBase2(3), Id(2,2), EigenMat(2,2), T(2,2) 
    REAL(KIND=dp) :: BPrinc(2,2), scale, Err
    REAL(KIND=dp) :: rK, hK, h1, h2, pk(3), xi, eta, z(MaxPatchNodes), x1(4), x2(4), uk, vk
    REAL(KIND=dp) :: p(3), r(2), delta(2), DerMat(2,2), ptarget(2), p0(2)
    REAL(KIND=dp) :: APar, BPar, DPar, EPar, FPar, GPar
    REAL(KIND=dp) :: FrameData(FrameDataSize), NodesArray(3*MaxPatchNodes)
    REAL(KIND=dp) :: GBasis(MaxPatchNodes), BubbleNodesDelta(4,3)

    SAVE Nodes, GElement
!------------------------------------------------------------------------------
    Family = GetElementFamily(Element)

    IF (PRESENT(ReparametrizeMesh)) THEN
      UseMeshOnly = ReparametrizeMesh
    ELSE
      UseMeshOnly = .FALSE.
    END IF

    IF (PRESENT(xi1) .AND. PRESENT(xi2)) THEN
      u = xi1
      v = xi2
    ELSE
      ! Evaluate quantities of interest at the element centre (p-element parametrization):
      SELECT CASE(Family)
      CASE(3)
        u = 0.0d0
        v = sqrt(3.0d0)/3.0d0
      CASE(4)
        u = 0.0d0
        v = 0.0d0
      END SELECT      
    END IF

    IF ( PRESENT(TaylorApproximation) ) THEN
      ComputeTaylorPolynomial = TaylorApproximation
    ELSE
      ComputeTaylorPolynomial = .FALSE.
    END IF

    IF ( PRESENT(LagrangeNodes) ) THEN
      ApproximatePlaneDomain = .TRUE.
      LagrangeNodes = 0.0d0
    ELSE
      ApproximatePlaneDomain = .FALSE.
    END IF

    IF ( PRESENT(d) ) THEN
      CheckOrientation = .TRUE.
    ELSE
      CheckOrientation = .FALSE.
    END IF

    IF ( PRESENT(SaveProperties) ) THEN
      WriteElementProperties = SaveProperties
    ELSE
      WriteElementProperties = .FALSE.
    END IF

    !-----------------------------------------------------------------------------
    ! The Taylor polynomial coefficients and the shape of the projected domain
    ! can be approximated only when the values of u and v correspond to the 
    ! element centre
    !-----------------------------------------------------------------------------
    IF ( (PRESENT(xi1) .AND. PRESENT(xi2)) .AND. &
        ApproximatePlaneDomain .OR. ComputeTaylorPolynomial) THEN
      SELECT CASE(Family)
      CASE(3)
        xi = 0.0d0
        eta = sqrt(3.0d0)/3.0d0
      CASE(4)
        xi = 0.0d0
        eta = 0.0d0
      END SELECT

      IF ( (ABS(u-xi) > AEPS) .OR. (ABS(v-eta) > AEPS) ) THEN
        IF (ApproximatePlaneDomain) THEN
          CALL Warn('LinesOfCurvatureFrame','Domain shape computation rejected!')
          ApproximatePlaneDomain = .FALSE.
        END IF
        IF (ComputeTaylorPolynomial) THEN
          CALL Warn('LinesOfCurvatureFrame','Taylor coefficients computation rejected!')      
          ComputeTaylorPolynomial = .FALSE.
        END IF
      END IF
    END IF

    CALL GetElementNodes( Nodes )

    ! -------------------------------------------------------------------------
    ! Compute lines of curvature directions at the surface point. To begin with,
    ! we compute the first and second fundamental forms when the reference
    ! element coordinates are used as curvilinear coordinates on the shell
    ! surface.
    ! -------------------------------------------------------------------------
    stat = BlendingSurfaceInfo(Element, Nodes, u, v, deta, a1, a2, a3, a, b, X0, &
        UseMeshOnly=UseMeshOnly)      

    !---------------------------------------------------------------------------
    ! The computation of principal directions via solving an eigenvalue problem.
    ! TO DO: create a subroutine for the 2-dimensional eigenvalue problem
    !---------------------------------------------------------------------------
    ContravA(1,1) = 1/deta * a(2,2)
    ContravA(2,2) = 1/deta * a(1,1)
    ContravA(2,1) = -1/deta * a(2,1)
    ContravA(1,2) = -1/deta * a(1,2)

    DualBase1(:) = ContravA(1,1)*a1(:) + ContravA(1,2)*a2(:)
    DualBase2(:) = ContravA(2,1)*a1(:) + ContravA(2,2)*a2(:)    

    c(1:2,1:2) = MATMUL(b,ContravA)
    detc = c(1,1)*c(2,2)-c(1,2)*c(2,1) ! = the Gaussian curvature K
    trc = c(1,1) + c(2,2)              ! = 2H, H = the mean curvature
    discriminant = trc**2 - 4.0d0*detc
    !--------------------------------------------------------------------------
    ! Allow some arithmetic inaccuracy: The discriminant = 4H^2 - 4K shouldn't
    ! be negative. The zero value happens at an umbilical point, so 
    ! a (small) negative value is taken to indicate an umbilical point
    !--------------------------------------------------------------------------
    IF (discriminant < 0.0d0) THEN
      CALL Warn('LinesOfCurvatureFrame', 'A negative discriminant in curvature computation?')
      discriminant = 0.0d0
      Umbilical = .TRUE.
    ELSE
      Umbilical = .FALSE.
    END IF

    !----------------------------------------------------------------------
    ! Check if the user wants to guide the reparametrization by defining
    ! a spherical part:
    !----------------------------------------------------------------------       
    Umbilical = GetLogical(GetBodyParams(Element), 'Spherical Body', Found)

    !--------------------------------------------------------------
    ! Order the eigenvalues by their absolute values:
    !--------------------------------------------------------------
    IF (trc>0.0d0) THEN
      lambda2 = 0.5d0 * ( trc + SQRT(discriminant) )
      lambda1 = 0.5d0 * ( trc - SQRT(discriminant) )
    ELSE
      lambda2 = 0.5d0 * ( trc - SQRT(discriminant) )
      lambda1 = 0.5d0 * ( trc + SQRT(discriminant) )
    END IF
    !print *, 'Eigenvalues=', lambda1,lambda2
    LambdaMax = MAX(ABS(Lambda1), ABS(Lambda2))

    Planar = (ABS(Lambda1) < EPSILON(1.0)) .AND. (ABS(Lambda2) < EPSILON(1.0))

    !-----------------------------------------------------------------
    ! Another planarity check may have been done. A positive result
    ! of an earlier data test will be respected.
    !-----------------------------------------------------------------
    IF (PRESENT(PlanarSurface)) THEN
      IF (PlanarSurface .AND. .NOT. Planar) THEN
        CALL Info( 'LinesOfCurvatureFrame', &
            'Planarity checks produce different results, rejecting negative result', &
            Level=10 )
        Planar = .TRUE.
      END IF
    END IF

    IF (.NOT. Planar .AND. .NOT. Umbilical) THEN
      ! ------------------------------------------------------------------------------ 
      ! The test for umbilical points depends on a rough estimate of the element size:
      ! ------------------------------------------------------------------------------
      rK = 0.0d0
      SELECT CASE(Family)
      CASE(3)
        rK = MAX(rK, SQRT((Nodes % x(2) - Nodes % x(1))**2 + (Nodes % y(2) - Nodes % y(1))**2 + &
            (Nodes % z(2) - Nodes % z(1))**2))
        rK = MAX(rK, SQRT((Nodes % x(3) - Nodes % x(2))**2 + (Nodes % y(3) - Nodes % y(2))**2 + &
            (Nodes % z(3) - Nodes % z(2))**2))
        rK = MAX(rK, SQRT((Nodes % x(3) - Nodes % x(1))**2 + (Nodes % y(3) - Nodes % y(1))**2 + &
            (Nodes % z(3) - Nodes % z(1))**2))
      CASE(4)
        rK = MAX(rK, SQRT((Nodes % x(3) - Nodes % x(1))**2 + (Nodes % y(3) - Nodes % y(1))**2 + &
            (Nodes % z(3) - Nodes % z(1))**2))
        rK = MAX(rK, SQRT((Nodes % x(4) - Nodes % x(2))**2 + (Nodes % y(4) - Nodes % y(2))**2 + &
            (Nodes % z(4) - Nodes % z(2))**2))
      END SELECT

      ! If the point is (nearly) umbilical, two reparametrization coefficients
      ! may be (nearly) singular. A minimal stability condition to obtain a regular
      ! reparametrization is of the form
      !
      !      H/(2L^2) < | Lambda1-Lambda2 |
      !
      ! with H and L being a mesh size parameter and a characteristic length of the geometry. 
      ! If the stability may be expected to be problematic, we avoid using (nearly) singular 
      ! coefficients by using a reparametrization which is suitable for a sphere. 
      ! Here we choose L = R_min, where R_min is the minimal principal radius of curvature.

      Umbilical = ABS(Lambda1-Lambda2) < 0.5d0 * rK * LambdaMax**2
    END IF

    !-----------------------------------------------------------------
    ! Compute the eigenvectors: 
    !-----------------------------------------------------------------
    IF (Planar .OR. Umbilical) THEN
      ! ------------------------------------------------------------------------
      ! For planar and umbilical points the principal coordinate directions are 
      ! not unique. Select one of the possibilities:
      ! ------------------------------------------------------------------------
      GlobPDir1(1:3) = 1.0d0/SQRT(SUM(a1(1:3)**2)) * a1(1:3)
      GlobPDir2(:) = a2(:) - DOT_PRODUCT(a2,GlobPDir1)*GlobPDir1(:)
      GlobPDir2(1:3) = 1.0d0/SQRT(SUM(GlobPDir2(1:3)**2)) * GlobPDir2(1:3)
    ELSE
      Id = 0.0d0
      Id(1,1) = 1.0d0
      Id(2,2) = 1.0d0
      ! ------------------------------------------------------------------------
      ! The rows of the matrix equation giving the eigenvectors are multiples
      ! of the same equation. To avoid possible underflow, we employ the row
      ! which has the largest 1-norm.
      ! ------------------------------------------------------------------------
      EigenMat(1:2,1:2) = c(:,:)-lambda1*Id
      IF ( (ABS(EigenMat(1,1)) + ABS(EigenMat(1,2))) > &
          (ABS(EigenMat(2,1)) + ABS(EigenMat(2,2))) ) THEN
        i = 1
      ELSE
        i = 2
      END IF
      scale = ABS(EigenMat(i,1)) + ABS(EigenMat(i,2))
      EigenMat(i,:) = 1/scale * EigenMat(i,:)
      PDir1(1) = EigenMat(i,2)
      PDir1(2) = -EigenMat(i,1)
      PDir1(1:2) = 1.0d0/SQRT(SUM(PDir1(1:2)**2)) * PDir1(1:2)

      EigenMat(1:2,1:2) = c(:,:)-lambda2*Id
      IF ( (ABS(EigenMat(1,1)) + ABS(EigenMat(1,2))) > &
          (ABS(EigenMat(2,1)) + ABS(EigenMat(2,2))) ) THEN
        i = 1
      ELSE
        i = 2
      END IF
      scale = ABS(EigenMat(i,1)) + ABS(EigenMat(i,2))
      EigenMat(i,:) = 1/scale * EigenMat(i,:)
      PDir2(1) = EigenMat(i,2)
      PDir2(2) = -EigenMat(i,1)
      PDir2(1:2) = 1.0d0/SQRT(SUM(PDir2(1:2)**2)) * PDir2(1:2)

      ! -------------------------------------------------------------------------
      ! The eigenvector expressed in terms of the covariant components over
      ! the dual base vectors of the surface: 
      ! --------------------------------------------------------------------------
      GlobPDir1(:) = PDir1(1) * DualBase1(:) + PDir1(2) * DualBase2(:)
      GlobPDir2(:) = PDir2(1) * DualBase1(:) + PDir2(2) * DualBase2(:)    
      GlobPDir1(1:3) = 1.0d0/SQRT(SUM(GlobPDir1(1:3)**2)) * GlobPDir1(1:3)
      GlobPDir2(1:3) = 1.0d0/SQRT(SUM(GlobPDir2(1:3)**2)) * GlobPDir2(1:3)
    END IF

    ! ----------------------------------------------------------------------
    ! Ensure finally that the local lines of curvature parameterization is
    ! right-handed and the normal agrees with the input data.
    !-----------------------------------------------------------------------
    IF (CheckOrientation) THEN
      IF ( DOT_PRODUCT(d, CrossProduct(GlobPDir1,GlobPDir2)) < 0.0d0 ) &
          GlobPDir2(:) = -1.0d0 * GlobPDir2(:)
    END IF
    GlobPDir3 = CrossProduct(GlobPDir1,GlobPDir2)

    ! ----------------------------------------------------------------------
    ! Make a change of basis and transform the covariant curvature tensor
    ! to the components along the principal axes
    !-----------------------------------------------------------------------   
    !T(1,1) = LDot(GlobPDir1,DualBase1)
    !T(2,1) = LDot(GlobPDir1,DualBase2)
    !T(1,2) = LDot(GlobPDir2,DualBase1)
    !T(2,2) = LDot(GlobPDir2,DualBase2)
    !BPrinc(1:2,1:2) = MATMUL(TRANSPOSE(T),MATMUL(b,T))
    !print *, 'pdir1 = ', GlobPDir1(:)
    !print *, 'pdir2 = ', GlobPDir2(:)
    !print *, 'pdir3 = ', GlobPDir3(:)
    !print *, 'T11=', BPrinc(1,1)
    !print *, 'T12=', BPrinc(1,2)
    !print *, 'T21=', BPrinc(2,1)
    !print *, 'T22=', BPrinc(2,2)

    IF (ApproximatePlaneDomain .OR. ComputeTaylorPolynomial) THEN
      IF (.NOT. ASSOCIATED(GElement)) GElement => AllocateElement()
      SELECT CASE (Family)
      CASE(3)
        GElement % Type => GetElementType(310,.FALSE.)
        !
        ! Make sure that all entries are initiated even with p-elements
        !
        GElement % Type % NodeU(4:10) = (/ 1.0d0/3.0d0, 2.0d0/3.0d0, 2.0d0/3.0d0, 1.0d0/3.0d0, &
            0.0d0, 0.0d0, 1.0d0/3.0d0 /)
        GElement % Type % NodeV(4:10) = (/ 0.0d0, 0.0d0, 1.0d0/3.0d0, 2.0d0/3.0d0, 2.0d0/3.0d0, &
            1.0d0/3.0d0, 1.0d0/3.0d0 /)
      CASE(4)
        GElement % Type => GetElementType(416,.FALSE.)
        !
        ! Make sure that all entries are initiated even with p-elements
        !
        GElement % Type % NodeU(5:16) = (/ -1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0, 1.0d0, 1.0d0/3.0d0, &
            -1.0d0/3.0d0, -1.0d0, -1.0d0, -1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0, -1.0d0/3.0d0 /)
        GElement % Type % NodeV(5:16) = (/ -1.0d0, -1.0d0, -1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0, 1.0d0, &
            1.0d0/3.0d0, -1.0d0/3.0d0, -1.0d0/3.0d0, -1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0 /)
      END SELECT
    END IF

    IF (ApproximatePlaneDomain) THEN
      DO j=1,GElement % Type % NumberOfNodes
        IF (j > Family .AND. GElement % Type % NumberOfNodes /= Element % Type % NumberOfNodes) THEN
          ! First we may need to map the Lagrange element coordinates to the ones of
          ! the p-element:
          IF (Family==3) THEN
            xi = -1.0d0 + 2.0d0*GElement % Type % NodeU(j) + GElement % Type % NodeV(j)
            eta = SQRT(3.0d0)*GElement % Type % NodeV(j)
          ELSE
            xi = GElement % Type % NodeU(j)
            eta = GElement % Type % NodeV(j)
          END IF

          stat = BlendingSurfaceInfo(Element, Nodes, xi, eta, deta, &
              a1, a2, a3, a, b, p, UseMeshOnly=UseMeshOnly)
        ELSE
          p(1) = Nodes % x(j)
          p(2) = Nodes % y(j)
          p(3) = Nodes % z(j)
        END IF
        p = p - X0
        LagrangeNodes(j,1) = DOT_PRODUCT(p,GlobPDir1)
        LagrangeNodes(j,2) = DOT_PRODUCT(p,GlobPDir2)
        LagrangeNodes(j,3) = DOT_PRODUCT(p,GlobPDir3)
        ! print *, 'NODE, x1,x2,x3', LagrangeNodes(j,1), LagrangeNodes(j,2), LagrangeNodes(j,3)
      END DO
    END IF

    IF (ComputeTaylorPolynomial) THEN   
      !--------------------------------------------------------------------
      ! Compute the Taylor polynomial coefficients for refining the lines of
      ! curvature parameterization. First, estimate the size of the planar 
      ! domain obtained via the projection to fit a regular stencil
      !--------------------------------------------------------------------
      rK = HUGE(rK)
      SELECT CASE(Family)
      CASE(3)
        ! Estimate the length of the longest edge:
        p0(1) = LagrangeNodes(2,1) - LagrangeNodes(1,1)
        p0(2) = LagrangeNodes(2,2) - LagrangeNodes(1,2)
        rK = MIN(rK,SQRT(SUM(p0(:)**2)))
        p0(1) = LagrangeNodes(3,1) - LagrangeNodes(2,1)
        p0(2) = LagrangeNodes(3,2) - LagrangeNodes(2,2)
        rK = MIN(rK,SQRT(SUM(p0(:)**2)))
        p0(1) = LagrangeNodes(3,1) - LagrangeNodes(1,1)
        p0(2) = LagrangeNodes(3,2) - LagrangeNodes(1,2)
        rK = MIN(rK,SQRT(SUM(p0(:)**2)))
        rK = 0.5d0 * rK
      CASE(4)
        ! Find the locations of the edge curve mid-points related to the planar domain 
        DO j=1,4
          SELECT CASE(j)
          CASE(1)
            xi = 0.0d0
            eta = -1.0d0
          CASE(2)
            xi = 1.0d0
            eta = 0.0d0
          CASE(3)
            xi = 0.0d0
            eta = 1.0d0      
          CASE(4)
            xi = -1.0d0
            eta = 0.0d0
          END SELECT

          CALL NodalBasisFunctions2D(GBasis, GElement, xi, eta)
          m = GElement % Type % NumberOfNodes
          p0(1) = SUM(LagrangeNodes(1:m,1) * GBasis(1:m))
          p0(2) = SUM(LagrangeNodes(1:m,2) * GBasis(1:m))    

          rK = MIN(rK,SQRT(SUM(p0(:)**2)))
        END DO
      END SELECT

      ! --------------------------------------------------------------------
      ! Create a regular 4X4-subgrid around the local origin to simplify
      ! the evaluation of higher order derivatives related to the Taylor
      ! polynomial. It is supposed that a square [-rK/8,rK/8]^2 is embedded
      ! into the plane domain obtained via the projection.
      ! TO DO: FIGURE OUT THE PRECISE SIZE OF A SQUARE THAT CAN BE EMBEDDED 
      ! --------------------------------------------------------------------
      hk = rK/4.0d0               ! The width of stencil
      x1(1) = -hk/2.0d0
      x1(2) = x1(1) + hk/3.0d0
      x1(3) = x1(2) + hk/3.0d0
      x1(4) = hk/2.0d0
      p0(1) = DOT_PRODUCT(X0,GlobPDir1)
      p0(2) = DOT_PRODUCT(X0,GlobPDir2)
      GridPoint = 0
      DO j=1,4
        ! An initial guess for iteration:
        SELECT CASE(Family)
        CASE(3)
          vk = sqrt(3.0d0)/3.0d0
        CASE(4)
          vk = -1.0d0/6.0d0 + (j-1)/9.0d0
        END SELECT
        DO i=1,4
          GridPoint = GridPoint + 1
          ! --------------------------------------------------------------
          ! Solve the reference element point (u,v) related to the blending
          ! map f(u,v) such that <f(u,v),e_i> = <X0,e_i> + ptarget(i),
          ! with i=1,2 and e_i the basis vectors of the elementwise frame.
          ! --------------------------------------------------------------
          ptarget(1) = x1(i)
          ptarget(2) = x1(j)

          ! An initial guess for iteration:
          SELECT CASE(Family)
          CASE(3)
            uk = 0.0d0  
          CASE(4)
            uk = -1.0d0/6.0d0 + (i-1)/9.0d0     
          END SELECT

          Converged = .FALSE.
          DO k=1,GeometryMaxIters
            stat = BlendingSurfaceInfo( Element, Nodes, uk, vk, &
                deta, a1, a2, a3, a, b, pk, UseMeshOnly=UseMeshOnly)             

            r(1) = p0(1) + ptarget(1) - DOT_PRODUCT(pk,GlobPDir1)
            r(2) = p0(2) + ptarget(2) - DOT_PRODUCT(pk,GlobPDir2)
            err = SQRT(SUM(r(:)**2))

            IF ( err < GeometryEpsilon*hK ) THEN
              Converged = .TRUE.
              EXIT
            END IF

            DerMat(1,1) = DOT_PRODUCT(a1,GlobPDir1)
            DerMat(1,2) = DOT_PRODUCT(a2,GlobPDir1)
            DerMat(2,1) = DOT_PRODUCT(a1,GlobPDir2)
            DerMat(2,2) = DOT_PRODUCT(a2,GlobPDir2)

            CALL SolveLinSys2x2(DerMat,delta,r)

            uk = uk + delta(1)
            vk = vk + delta(2)
          END DO
          IF (.NOT. Converged) CALL Fatal('LinesOfCurvatureFrame', &
              'Cannot create a subgrid: elements may not be shape-regular')
          SELECT CASE(Family)
          CASE(3)
            IF (ABS(uk) > 1.0d0+1.0d-3 .OR. vk > sqrt(3.0d0)+1.0d-3 .OR. &
                vk < -1.0d-3) THEN
              CALL Fatal('LinesOfCurvatureFrame', &
                  'Cannot create a subgrid: a subgrid point outside element')
            END IF
          CASE(4)
            IF (ABS(uk) > 1.0d0+1.0d-3 .OR. ABS(vk) > 1.0d0+1.0d-3) &
                CALL Fatal('LinesOfCurvatureFrame', &
                'Cannot create a subgrid: a subgrid point outside element')
          END SELECT

          ! The third coordinate of the point on the blending surface:
          z(GridPoint) = DOT_PRODUCT(pk - X0,GlobPDir3)
        END DO
      END DO

      !--------------------------------------------------------------------------------
      ! Compute the Taylor polynomial coefficients for refining the lines of curvature
      ! parameterization (via a stencil based on the finite element approximation).
      ! The Taylor polynomial is written as z(x,y) = 
      ! 1/2 A x^2 + 1/2 B y^2 + 1/6 D x^3 + 1/2 E x^2 y + 1/2 F x y^2 + 1/6 G y^3.
      !--------------------------------------------------------------------------------
      h1 = hk
      h2 = hk

      APar = (-9.0d0*(z(1) + 9.0d0*z(10) + 9.0d0*z(11) - 9.0d0*z(12) + z(13) - z(14) - z(15) + z(16) - &
          z(2) - z(3) + z(4) - 9.0d0*z(5) + 9.0d0*z(6) + 9.0d0*z(7) - 9.0d0*z(8) - 9.0d0*z(9)))/(32.0d0*h1**2) 
      BPar =  (-9.0d0*(z(1) + 9.0d0*z(10) + 9.0d0*z(11) - z(12) + z(13) - 9.0d0*z(14) - 9.0d0*z(15) + z(16) - &
          9.0d0*z(2) - 9.0d0*z(3) + z(4) - z(5) + 9.0d0*z(6) + 9.0d0*z(7) - z(8) - z(9)))/(32.0d0*h2**2)

      DPar = (27.0d0*(z(1) + 27.0d0*z(10) - 27.0d0*z(11) + 9.0d0*z(12) + z(13) - 3.0d0*z(14) + 3.0d0*z(15) - &
          z(16) - 3.0d0*z(2) + 3.0d0*z(3) - z(4) - 9.0d0*z(5) + 27.0d0*z(6) - 27.0d0*z(7) + 9.0d0*z(8) - &
          9.0d0*z(9)))/(16.0d0*h1**3)
      EPar = (9.0d0*(z(1) - 27.0d0*z(10) - 27.0d0*z(11) + 27.0d0*z(12) - z(13) + z(14) + z(15) - z(16) - &
          z(2) - z(3) + z(4) - 27.0d0*z(5) + 27.0d0*z(6) + 27.0d0*z(7) - 27.0d0*z(8) + &
          27.0d0*z(9)))/(16.0d0*h1**2*h2)
      FPar = (9.0d0*(z(1) + 27.0d0*z(10) - 27.0d0*z(11) + z(12) + z(13) - 27.0d0*z(14) + 27.0d0*z(15) - &
          z(16) - 27.0d0*z(2) + 27.0d0*z(3) - z(4) - z(5) + 27.0d0*z(6) - 27.0d0*z(7) + z(8) - &
          z(9)))/(16.0d0*h1*h2**2)
      GPar = (27.0d0*(z(1) - 27.0d0*z(10) - 27.0d0*z(11) + 3.0d0*z(12) - z(13) + 9.0d0*z(14) + 9.0d0*z(15) - &
          z(16) - 9.0d0*z(2) - 9.0d0*z(3) + z(4) - 3.0d0*z(5) + 27.0d0*z(6) + 27.0d0*z(7) - 3.0d0*z(8) + &
          3.0d0*z(9)))/(16.0d0*h2**3)

      TaylorParams(1) = APar
      TaylorParams(2) = BPar
      TaylorParams(3) = DPar
      TaylorParams(4) = EPar
      TaylorParams(5) = FPar
      TaylorParams(6) = GPar

      ! Check that two coefficients of the Taylor polynomial are close to the eigenvalues:
      !
      IF (Planar) THEN
        IF ( (ABS(APar) > 1000.0*EPSILON(1.0)) .AND. (ABS(BPar) > 1000.0*EPSILON(1.0)) ) THEN
          CALL Warn('LinesOfCurvatureFrame', 'Possibly inaccurate Taylor polynomial (planar point)')
          print *, 'Apar,Lambda1=', Apar, Lambda1
          print *, 'Bpar,Lambda2=', Bpar, Lambda2
        END IF
      ELSE
        err = ABS( ABS(BPar)-ABS(Lambda2) ) / LambdaMax
        err = MAX( err, ABS(ABS(APar)-ABS(lambda1))/LambdaMax )
        IF ( err > 10.0d-2) THEN
          CALL Warn('LinesOfCurvatureFrame', 'Possibly inaccurate Taylor polynomial')
          print *, '|Apar-Lambda1|/LambdaMax=', ABS(ABS(APar)-ABS(lambda1))/LambdaMax
          print *, '|Bpar-Lambda2|/LambdaMax=', ABS(ABS(BPar)-ABS(lambda2))/LambdaMax
        END IF
      END IF

    END IF

    ! -------------------------------------------------------------------------
    ! Indicate whether the surface at the specified point 
    ! can be approximated by a plane up to errors O(h^3).
    ! -------------------------------------------------------------------------
    IF ( PRESENT(PlanarPoint) ) PlanarPoint = Planar
    IF ( PRESENT(UmbilicalPoint) ) THEN
      IF (Planar) THEN
        UmbilicalPoint = .FALSE.
      ELSE
        UmbilicalPoint = Umbilical
      END IF
    END IF

    IF ( PRESENT(e1) ) e1(:) = GlobPDir1(:)
    IF ( PRESENT(e2) ) e2(:) = GlobPDir2(:)
    IF ( PRESENT(e3) ) e3(:) = GlobPDir3(:)
    IF ( PRESENT(o) ) o(:) = X0(:)

    
    ! Return information about the mesh resolution of geometry:
    IF ( PRESENT(SizeRadiusRatio) ) THEN
      SizeRadiusRatio = 0.0d0
      IF (ComputeTaylorPolynomial) THEN
        IF (.NOT. Planar) SizeRadiusRatio = 2.0d0 * rK * LambdaMax
      ELSE
        CALL Warn('LinesOfCurvatureFrame', &
            'Cannot estimate h/R ratio without ComputeTaylorPolynomial')
      END IF
    END IF

    ! -------------------------------------------------------------------------
    ! Save the quantities computed as elementwise properties:
    ! -------------------------------------------------------------------------
    IF (WriteElementProperties) THEN
      FrameData(FrameBasis1) = GlobPDir1(1:3)
      FrameData(FrameBasis2) = GlobPDir2(1:3)
      FrameData(FrameBasis3) = GlobPDir3(1:3)
      FrameData(FrameOrigin) = X0(1:3)
      CALL SetElementProperty('element frame', FrameData(1:FrameDataSize), Element) 

      IF (ComputeTaylorPolynomial) CALL SetElementProperty('taylor parameters', &
          TaylorParams(1:6), Element)

      IF (Planar) THEN
        PlanarFlag = 1.0d0
      ELSE
        PlanarFlag = -1.0d0
      END IF
      CALL SetElementProperty('planar point', PlanarFlag, Element)

      IF (Umbilical .AND. .NOT. Planar) THEN
        UmbilicalFlag = 1.0d0
      ELSE
        UmbilicalFlag = -1.0d0
      END IF
      CALL SetElementProperty('umbilical point', UmbilicalFlag, Element) 
    END IF
    
    !print *, 'o=', X0(1:3)
    !print *, 'e1=', GlobPDir1(:)
    !print *, 'e2=', GlobPDir2(:)
    !print *, 'e3=', GlobPDir3(:)
    !print *, 'difference of Taylor params=', ABS(APar-BPar)/MAX(ABS(APar),ABS(BPar))
    !print *, 'Umbilical=',Umbilical
    !print *, 'remainder=', maxval(ABS(TaylorParams(3:6))) * 2.0d0 * rK * LambdaMax
    !IF (.NOT. Umbilical) PRINT *, 'c1=', FPar/(BPar - APar)
    !IF (.NOT. Umbilical) print *, 'c2=', EPar/(BPar - APar)

!------------------------------------------------------------------------------
  END SUBROUTINE LinesOfCurvatureFrame
!------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
! Obtain the nodes of a coordinate patch for improved lines of curvature
! parameterization. This subroutine solves nodewise a root finding problem of 
! the type g(y) - x = 0, where x is a given node on the plane domain S and g is 
! a given nonlinear transformation from a subset K of R^2 onto S. The argument 
! LocalFrameNodes defines x, the points y are saved as the elementwise property 
! 'patch nodes' and the form of g is defined by the parameters TaylorParams. 
! If ZNodes is supplied, approximations of nodal z-coordinates are computed using 
! a third-order Taylor polynomial in the coordinates y of the final patch K. This option
! can be used to cross-check different approximations but may not have final utility.
!-----------------------------------------------------------------------------------
  SUBROUTINE LinesOfCurvaturePatch(Element, LocalFrameNodes, TaylorParams, &
      Family, PlanarSurface, Umbilical, ZNodes)
!-----------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    REAL(KIND=dp), TARGET, INTENT(IN) :: LocalFrameNodes(MaxPatchNodes,2)
    REAL(KIND=dp), POINTER, INTENT(IN) :: TaylorParams(:)
    INTEGER, INTENT(IN) :: Family
    LOGICAL, OPTIONAL, INTENT(IN) :: PlanarSurface
    LOGICAL, INTENT(IN) :: Umbilical
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: ZNodes(MaxPatchNodes)
!-----------------------------------------------------------------------------------
    LOGICAL :: Planar, Parabolic, Converged
    INTEGER :: i, k, n

    REAL(KIND=dp) :: PatchNodes(MaxPatchNodes,2)
    REAL(KIND=dp) :: NodesArray(2*MaxPatchNodes)
    REAL(KIND=dp) :: c1, c2, c4, c5, c6, c7, b5, b7 
    REAL(KIND=dp) :: hK, y1, y2
    REAL(KIND=dp) :: r(2), delta(2), DerMat(2,2), err
    REAL(KIND=dp) :: x, y
    REAL(KIND=dp) :: A, B, D, E, F, G
!----------------------------------------------------------------------------------- 
    Planar = .FALSE.
    IF (PRESENT(PlanarSurface)) THEN
      Planar = PlanarSurface
      IF (Planar) THEN
        !-----------------------------------------------------------------------
        ! In the special case of planar surface there is no need to improve
        ! parametrization: make an early exit
        !-----------------------------------------------------------------------
        IF (PRESENT(ZNodes)) ZNodes = 0.0d0

        ! Save the patch nodes as an elementwise property:
        NodesArray(1:MaxPatchNodes) = LocalFrameNodes(1:MaxPatchNodes,1)
        NodesArray(MaxPatchNodes+1:2*MaxPatchNodes) = LocalFrameNodes(1:MaxPatchNodes,2)
        CALL SetElementProperty('patch nodes', NodesArray(1:2*MaxPatchNodes), Element)

        RETURN
      END IF
    END IF

    A = TaylorParams(1)
    B = TaylorParams(2)
    D = TaylorParams(3)
    E = TaylorParams(4)
    F = TaylorParams(5)
    G = TaylorParams(6)

    Parabolic = .FALSE.
    IF (.NOT. Umbilical) Parabolic = ABS(A)/ABS(B) < ParabolicDelta

    !---------------------------------------------------------------------------------
    ! The following constants relate to the most general third-order polynomial 
    ! perturbations of the coordinate functions; cf. the definitions of x and y
    ! below.
    !---------------------------------------------------------------------------------
    IF (Umbilical) THEN
      !---------------------------------------------------------------------------------
      ! If the point is nearly umbilical, then the reparametrization is done as for a sphere.
      ! This is not fully accurate when the Taylor polynomial doesn't correspond to 
      ! a sphere, with F /= 0 or E /= 0. Then the error is that the second diagonal form
      ! is diagonal only up to an error O(h).
      !---------------------------------------------------------------------------------
      c1 = 0.0d0       ! Follows from Taylor's polynomial for a sphere
      c2 = 0.0d0       ! Follows from Taylor's polynomial for a sphere
    ELSE
      c1 = F/(B - A)   ! To make the 2nd fundamental form diagonal
      c2 = E/(B - A)   ! To make the 2nd fundamental form diagonal
    END IF

    IF (Parabolic) THEN
      c4 = 0.0d0
    ELSE
      c4 = -0.5d0*(A**2 + A*B + B**2) ! To enforce symmetries
    END IF

    c5 = -2.0d0*c1*c2              ! To simplify the Christoffel symbols
    c6 = -A**2 - A*B - c4          ! An isothermal reparametrization
    c7 = -c5                       ! An isothermal reparametrization
    b7 = -B**2 + c6                ! An isothermal reparametrization
    b5 = -A*B - c6                 ! Orthogonal surface basis vectors up to O(h^3)
    ! b4 = -c5                     ! Orthogonal surface basis vectors up to O(h^3)
    ! b6 = -c7                     ! Orthogonal surface basis vectors up to O(h^3)

    ! Estimate the size of the domain to adapt a termination criterion
    SELECT CASE(Family)
    CASE(3)
      n = 10
      hk = SQRT((LocalFrameNodes(2,1)-LocalFrameNodes(1,1))**2 + &
          (LocalFrameNodes(2,2)-LocalFrameNodes(1,2))**2)
      hk = MAX( hk, SQRT((LocalFrameNodes(3,1)-LocalFrameNodes(2,1))**2 + &
          (LocalFrameNodes(3,2)-LocalFrameNodes(2,2))**2) )
      hk = MAX( hk, SQRT((LocalFrameNodes(3,1)-LocalFrameNodes(1,1))**2 + &
          (LocalFrameNodes(3,2)-LocalFrameNodes(1,2))**2) )
    CASE(4)
      n = 16
      hk = MAX( SQRT((LocalFrameNodes(3,1)-LocalFrameNodes(1,1))**2 + &
          (LocalFrameNodes(3,2)-LocalFrameNodes(1,2))**2), &
          SQRT((LocalFrameNodes(2,1)-LocalFrameNodes(4,1))**2 + &
          (LocalFrameNodes(2,2)-LocalFrameNodes(4,2))**2) )
    END SELECT

    PatchNodes = 0.0d0
    DO i=1,n
      ! Initial guess for the Newton iteration:
      y1 = LocalFrameNodes(i,1)   
      y2 = LocalFrameNodes(i,2)
      DO k=1,GeometryMaxIters

        x = y1 - c1 * y1 ** 2 / 0.2D1 + c2 * y1 * y2 + c1 * y2 ** 2 / 0.2D1 + c4 * y1 ** 3 / 0.6D1 + &
            c5 * y1 ** 2 * y2 / 0.2D1 + c6 * y1 * y2 ** 2 / 0.2D1 + c7 * y2 ** 3 / 0.6D1

        y = y2 - c2 * y1 ** 2 / 0.2D1 - c1 * y1 * y2 + c2 * y2 ** 2 / 0.2D1 - c5 * y1 ** 3 / 0.6D1 + &
            b5 * y1 ** 2 * y2 / 0.2D1 - c7 * y1 * y2 ** 2 / 0.2D1 + b7 * y2 ** 3 / 0.6D1

        r(1) = LocalFrameNodes(i,1) - x
        r(2) = LocalFrameNodes(i,2) - y

        err = SQRT(SUM(r(:)**2))
        IF ( err < GeometryEpsilon*hK ) THEN
          Converged = .TRUE.
          EXIT
        END IF

        DerMat(1,1) = 0.1D1 - c1 * y1 + c2 * y2 + c4 * y1 ** 2 / 0.2D1 + c5 * y1 * y2 + &
            c6 * y2 ** 2 / 0.2D1
        DerMat(1,2) = c2 * y1 + c1 * y2 + c5 * y1 ** 2 / 0.2D1 + c6 * y1 * y2 + &
            c7 * y2 ** 2 / 0.2D1
        DerMat(2,1) = -c2 * y1 - c1 * y2 - c5 * y1 ** 2 / 0.2D1 + b5 * y1 * y2 - &
            c7 * y2 ** 2 / 0.2D1
        DerMat(2,2) = 0.1D1 - c1 * y1 + c2 * y2 + b5 * y1 ** 2 / 0.2D1 - c7 * y1 * y2 + &
            b7 * y2 ** 2 / 0.2D1

        CALL SolveLinSys2x2(DerMat,delta,r)
        y1 = y1 + delta(1)
        y2 = y2 + delta(2)
      END DO
      IF (.NOT. Converged) CALL Fatal('LinesOfCurvaturePatch', 'Nonlinear iteration fails')
      PatchNodes(i,1) = y1
      Patchnodes(i,2) = y2
      !print *, 'NODE, y1,y2', i, y1,y2
    END DO

    IF (PRESENT(ZNodes)) THEN
      ! The z-coordinate values obtained using the Taylor polynomial over the final patch.
      ! This may not have final utility. 
      ZNodes = 0.0d0
      DO i=1,n
        x = PatchNodes(i,1)
        y = PatchNodes(i,2)
        IF (Umbilical) THEN
          ZNodes(i) = A * x ** 2 / 0.2D1 + B * y ** 2 / 0.2D1
        ELSE
          ZNodes(i) = A * x ** 2 / 0.2D1 + B * y ** 2 / 0.2D1 + (-A * c1 / 0.2D1 + &
              D / 0.6D1) * x ** 3 + (A * c2 - B * c2 / 0.2D1 + E / 0.2D1) * y * &
              x ** 2 + (A * c1 / 0.2D1 - B * c1 + F / 0.2D1) * y ** 2 * x + (B * & 
              c2 / 0.2D1 + G / 0.6D1) * y ** 3
        END IF
      END DO
    END IF

    ! Save the patch nodes as an elementwise property:
    NodesArray(1:MaxPatchNodes) = PatchNodes(1:MaxPatchNodes,1)
    NodesArray(MaxPatchNodes+1:2*MaxPatchNodes) = PatchNodes(1:MaxPatchNodes,2)
    CALL SetElementProperty('patch nodes', NodesArray(1:2*MaxPatchNodes), Element)

!-----------------------------------------------------------------------------------
  END SUBROUTINE LinesOfCurvaturePatch
!-----------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
! Obtain the covariant surface basis vectors ai when the coordinates of the principal 
! coordinate patch x and y are given. The components of the two fundamental forms and 
! Christoffel symbols are also returned. Optionally the contravariant surface basis 
! vectors duali and the corresponding global coordinates of the surface point are
! computed. Here approximations are based on the Taylor polynomial expansion over
! the principal coordinate patch. 
!-----------------------------------------------------------------------------------
  SUBROUTINE SurfaceBasisVectors(x, y, TaylorParams, e1, e2, e3, o, a1, a2, a3, &
      A11, A22, SqrtDetA, B11, B22, C111, C112, C221, C222, C211, C212, &
      dual1, dual2, XGlob, YGlob, ZGlob, LowestOrderBasis, PlanarPoint, &
      Umbilical)
!-----------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: x, y
    REAL(KIND=dp), POINTER, INTENT(IN) :: TaylorParams(:)
    REAL(KIND=dp), INTENT(IN) :: e1(3), e2(3), e3(3), o(3)
    REAL(KIND=dp), INTENT(OUT) :: a1(3), a2(3), a3(3)
    REAL(KIND=dp), INTENT(OUT) :: A11, A22, SqrtDetA, B11, B22
    REAL(KIND=dp), INTENT(OUT) :: C111, C112, C221, C222, C211, C212
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: dual1(3), dual2(3)
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: XGlob, YGlob, ZGlob
    LOGICAL, OPTIONAL, INTENT(IN) :: LowestOrderBasis
    LOGICAL, OPTIONAL, INTENT(IN) :: PlanarPoint
    LOGICAL, INTENT(IN) :: Umbilical
!-----------------------------------------------------------------------------------
    LOGICAL :: ReturnDualBasis, ReturnGlobalCoords
    LOGICAL :: Planar, Parabolic
    REAL(KIND=dp) :: A, B, D, E, F, G
    REAL(KIND=dp) :: c1, c2, c4, c5, c6, c7, b5, b7, d1, d2
    REAL(KIND=dp) :: x1, x2, x3
    REAL(KIND=dp) :: Norm
!-----------------------------------------------------------------------------------
    ReturnDualBasis = PRESENT(dual1) .AND. PRESENT(dual2)
    ReturnGlobalCoords = PRESENT(XGlob) .AND. PRESENT(YGlob) .AND. PRESENT(ZGlob)

    Planar = .FALSE.
    IF (PRESENT(PlanarPoint)) THEN
      Planar = PlanarPoint
    END IF

    IF (Planar) THEN
      a1(:) = e1(:)
      a2(:) = e2(:)
      a3(:) = e3(:)
      A11 = 1.0d0
      A22 = 1.0d0
      SqrtDetA = 1.0d0
      B11 = 0.0d0
      B22 = 0.0d0
      C111 = 0.0d0 
      C112 = 0.0d0 
      C221 = 0.0d0 
      C222 = 0.0d0 
      C211 = 0.0d0 
      C212 = 0.0d0 
      IF (ReturnDualBasis) THEN
        dual1(:) = a1(:)
        dual2(:) = a2(:)
      END IF
      IF (ReturnGlobalCoords) THEN
        XGlob = o(1) + x*e1(1) + y*e2(1)
        YGlob = o(2) + x*e1(2) + y*e2(2)
        ZGlob = o(3) + x*e1(3) + y*e2(3)
      END IF
      RETURN
    END IF

    ! ---------------------------------------------------------------------------- 
    ! The coefficients of the Taylor polynomial: 
    ! ---------------------------------------------------------------------------- 
    A = TaylorParams(1)
    B = TaylorParams(2)
    D = TaylorParams(3)
    E = TaylorParams(4)
    F = TaylorParams(5)
    G = TaylorParams(6)

    Parabolic = .FALSE.
    IF (.NOT. Umbilical) Parabolic = ABS(A)/ABS(B) < ParabolicDelta

    ! ---------------------------------------------------------------------------- 
    ! The constants to diagonalize the fundamental forms
    ! ---------------------------------------------------------------------------- 
    IF (Umbilical) THEN
      c1 = 0.0d0       ! Follows from Taylor's polynomial for a sphere
      c2 = 0.0d0       ! Follows from Taylor's polynomial for a sphere
    ELSE
      c1 = F/(B - A)   ! To make the 2nd fundamental form diagonal
      c2 = E/(B - A)   ! To make the 2nd fundamental form diagonal
    END IF

    IF (Parabolic) THEN
      c4 = 0.0d0
    ELSE
      c4 = -0.5d0*(A**2 + A*B + B**2) ! To enforce symmetries
    END IF

    c5 = -2.0d0*c1*c2              ! To simplify the Christoffel symbols
    c6 = -A**2 - A*B - c4          ! An isothermal reparametrization
    c7 = -c5                       ! An isothermal reparametrization
    b7 = -B**2 + c6                ! An isothermal reparametrization
    b5 = -A*B - c6                 ! Orthogonal surface basis vectors up to O(h^3)

    ! ---------------------------------------------------------------------------- 
    ! The covariant surface basis vectors expanded up to O(h^3)
    ! ---------------------------------------------------------------------------- 
    a1(:) = (0.1D1 - c1 * x + c2 * y + c4 * x ** 2 / 0.2D1 + c5 * x * y + &
        c6 * y ** 2 / 0.2D1) * e1(:) + &
        (-c2 * x - c1 * y - c5 * x ** 2 / 0.2D1 + b5 * x * y - &
        c7 * y ** 2 / 0.2D1) * e2(:) + &
        (A * x + (-0.3D1 / 0.2D1 * A * c1 + D / 0.2D1) * x ** 2 + &
        (0.2D1 * A * c2 - B * c2 + E) * y * x + (A * c1 / 0.2D1 - B * c1 + &
        F / 0.2D1) * y ** 2) *  e3(:)

    a2(:) = (c2 * x + c1 * y + c5 * x ** 2 / 0.2D1 + c6 * x * y + &
        c7 * y ** 2 / 0.2D1) *  e1(:) + &
        (0.1D1 - c1 * x + c2 * y + b5 * x ** 2 / 0.2D1 - c7 * x * y + &
        b7 * y ** 2 / 0.2D1) * e2(:) + & 
        (B * y + (A * c2 - B * c2 / 0.2D1 + E / 0.2D1) * x ** 2 + &
        (A * c1 - 0.2D1 * B * c1 + F) * y * x + (0.3D1 / 0.2D1 * B * c2 + &
        G /0.2D1) * y ** 2) *  e3(:)

    a3(:) = CrossProduct(a1,a2)
    Norm = SQRT(SUM(a3(1:3)**2))
    a3(:) = 1/Norm * a3(:) 

    ! ----------------------------------------------------------------------------
    ! The metric surface tensor, with a12 = O(h^3):
    ! ----------------------------------------------------------------------------
    a11 = DOT_PRODUCT(a1,a1)
    a22 = DOT_PRODUCT(a2,a2)
    SqrtDetA = SQRT(a11*a22)

    ! ----------------------------------------------------------------------------
    ! The covariant components of the second fundamental form:
    ! ----------------------------------------------------------------------------
    IF (Umbilical) THEN
      B11 = A
      B22 = B
    ELSE
      B11 = A + (-2.0d0 * A * c1 + D) * x + (2.0d0 * A * c2 + E) * y
      B22 = B + (-2.0d0 * B * c1 + F) * x + (2.0d0 * B * c2 + G) * y
    END IF
    ! ----------------------------------------------------------------------------
    ! The Christoffel symbols Cijk = C_{ij}^k:
    ! ----------------------------------------------------------------------------
    d1 = A**2 + c4 - c1**2 + c2**2 
    d2 = A**2 + A*B + c4 - c1**2 + c2**2
    C111 = -c1 + d1 * x
    C221 = -C111
    C212 = C111
    C222 = c2 - d2 * y
    C211 = C222
    C112 = -C222

    ! ----------------------------------------------------------------------------
    ! The contravariant surface basis vectors:
    ! ----------------------------------------------------------------------------      
    IF (ReturnDualBasis) THEN
      dual1(:) = 1.0d0/a11 * a1(:)
      dual2(:) = 1.0d0/a22 * a2(:)    
    END IF

    ! ----------------------------------------------------------------------------
    ! The global coordinates computed via the local Taylor expansion:
    ! ----------------------------------------------------------------------------      
    IF (ReturnGlobalCoords) THEN

      x1 = x - c1 * x ** 2 / 0.2D1 + c2 * x * y + c1 * y ** 2 / 0.2D1 + c4 * x ** 3 / 0.6D1 + &
          c5 * x ** 2 * y / 0.2D1 + c6 * x * y ** 2 / 0.2D1 + c7 * y ** 3 / 0.6D1

      x2 = y - c2 * x ** 2 / 0.2D1 - c1 * x * y + c2 * y ** 2 / 0.2D1 - c5 * x ** 3 / 0.6D1 + &
          b5 * x ** 2 * y / 0.2D1 - c7 * x * y ** 2 / 0.2D1 + b7 * y ** 3 / 0.6D1

      IF (Umbilical) THEN
        x3 = A * x ** 2 / 0.2D1 + B * y ** 2 / 0.2D1
      ELSE
        x3 = A * x ** 2 / 0.2D1 + B * y ** 2 / 0.2D1 + (-A * c1 / 0.2D1 + &
            D / 0.6D1) * x ** 3 + (A * c2 - B * c2 / 0.2D1 + E / 0.2D1) * y * &
            x ** 2 + (A * c1 / 0.2D1 - B * c1 + F / 0.2D1) * y ** 2 * x + (B * & 
            c2 / 0.2D1 + G / 0.6D1) * y ** 3
      END IF

      XGlob = o(1) + x1*e1(1) + x2*e2(1) + x3*e3(1)
      YGlob = o(2) + x1*e1(2) + x2*e2(2) + x3*e3(2)
      ZGlob = o(3) + x1*e1(3) + x2*e2(3) + x3*e3(3)
    END IF

!-----------------------------------------------------------------------------------
  END SUBROUTINE SurfaceBasisVectors
!-----------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! The elementwise contribution to the system of discrete shell equations.
! This subroutine creates the tangential stiffness matrix for the computation
! of the solution increment to update the nonlinear iterate.
! The local DOFs always correspond to the displacement components along the
! principal axes. The transformation to global DOFs is done within this subroutine.
! The stiffness matrix K corresponding to the global DOFs is thus obtained as 
! K = Q^T k Q and the RHS vector F is obtained as F = Q^T f.
!
! This version is based on expressing the displacement vector in terms of
! the covariant components (the basis is orthogonal but not orthonormal).
!
! IMPORTANT REMARK: Currently strain reduction operators have been worked out
! only for the lowest-order Lagrange interpolation elements. Detecting a 
! p-element switches to the standard weak formulation which can give highly
! inaccurate results for thin shells! This subroutine does not attempt to
! handle p-elements with p > 1. Higher-order cases are now handled by the
! the alternate subroutine ShellLocalMatrixCartesian activated by having
! Cartesian Formulation = True.
!------------------------------------------------------------------------------
  SUBROUTINE ShellLocalMatrix(BGElement, n, nd, m, LocalSol, LargeDeflection, &
      NonlinearBending, StrainReductionMethod, MembraneStrainReductionMethod, ShearAlpha, &
      MembraneAlpha, StretchAlpha, Bubbles, DrillingDOFs, DrillingPar, RotateDOFs, &
      MassAssembly, HarmonicAssembly, RHSForce, Area, Error, BenchmarkProblem)
!------------------------------------------------------------------------------
    USE SolidMechanicsUtils, ONLY: StrainEnergyDensity, ShearCorrectionFactor
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: BGElement  ! An element of background mesh
    INTEGER, INTENT(IN) :: n                           ! The number of background element nodes
    INTEGER, INTENT(IN) :: nd                          ! The number of DOFs per component (after
                                                       ! static condensation if bubbles are used)
    INTEGER, INTENT(IN) :: m                           ! The number of DOFs per node
    REAL(KIND=dp), INTENT(IN) :: LocalSol(:,:)         ! The previous solution iterate
    LOGICAL, INTENT(IN) :: LargeDeflection             ! To activate nonlinear terms
    LOGICAL, INTENT(IN) :: NonlinearBending            ! The full nonlinearity of bending and transverse shear strains    
    INTEGER, INTENT(IN) :: StrainReductionMethod       ! The choice of strain reduction method
    INTEGER, INTENT(IN) :: MembraneStrainReductionMethod ! The choice of membrane strain reduction method    
    REAL(KIND=dp), INTENT(IN) :: ShearAlpha            ! A parameter for shear relaxation (correction) 
    REAL(KIND=dp), INTENT(IN) :: MembraneAlpha         ! A parameter for membrane relaxation
    REAL(KIND=dp), INTENT(IN) :: StretchAlpha          ! A parameter for membrane relaxation
    LOGICAL, INTENT(IN) :: Bubbles                     ! To indicate that bubble functions are used
    LOGICAL, INTENT(IN) :: DrillingDOFs                ! Switches to drilling DOFs (limited functionality)
    REAL(KIND=dp), INTENT(IN) :: DrillingPar           ! A stabilization parameter for drilling DOFs 
    LOGICAL, INTENT(IN) :: RotateDOFs                  ! Use rotated DOFs (a tentative option)
    LOGICAL, INTENT(IN) :: MassAssembly                ! To activate mass matrix integration
    LOGICAL, INTENT(IN) :: HarmonicAssembly            ! To activate the global mass matrix updates
    REAL(KIND=dp), INTENT(OUT) :: RHSForce(:)          ! Local RHS vector corresponding to external loads
    REAL(KIND=dp), INTENT(INOUT) :: Area               ! A variable for area compution
    REAL(KIND=dp), INTENT(INOUT) :: Error              ! A variable for error compution
    LOGICAL, INTENT(IN) :: BenchmarkProblem            ! To omit some terms in the strain energy 
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element => NULL()
    TYPE(Element_t), POINTER :: GElement => NULL()     
    TYPE(Nodes_t) :: Nodes, PNodes, PRefNodes
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: Stat, Found
    LOGICAL :: SecondOrder, PVersion
    LOGICAL :: UseShearCorrection, UseBubbles
    LOGICAL :: PlateBody                   ! To indicate that the surface is flat
    LOGICAL :: SphericalSurface            ! To indicate that the surface is considered to spherical

    INTEGER :: Family, ReducedStrainDim, MembraneStrainDim, ReductionMethod, MembraneReductionMethod
    INTEGER :: ShearReductionMethod, StretchReductionMethod
    INTEGER :: DOFs, BubbleDOFs, i, j, k, p, t, i0, j0, csize, GElementNodes 

    REAL(KIND=dp), POINTER :: TaylorParams(:)

    REAL(KIND=dp) :: PatchNodes(MaxPatchNodes,2) ! The nodes of principal coordinate patch
    REAL(KIND=dp) :: e1(3), e2(3), e3(3)         ! The basis of the local frame
    REAL(KIND=dp) :: o(3)                        ! The origin of the local frame

    ! Prepare for the scenario that one elementwise bubble per component can be added
    ! (check array sizes if more bubbles are used):
    REAL(KIND=dp) :: Stiff(m*nd+m,m*nd+m), Mass(m*nd+m,m*nd+m), Force(m*nd+m)
    REAL(KIND=dp) :: Damp(m*nd+m,m*nd+m)
    REAL(KIND=dp) :: BM(4,m*nd+m), BS(4,m*nd+m), BB(3,m*nd+m)
    REAL(KIND=dp) :: NonlinBM(4,m*nd+m), NonlinBS(4,m*nd+m), NonlinBB(3,m*nd+m), BWork(13,m*nd+m)
    REAL(KIND=dp) :: Basis(nd+1), dBasis(nd+1,3)
    REAL(KIND=dp) :: ReductionDOFsArray(4,2*nd+2) ! Large enough for p=1 with one bubble

    REAL(KIND=dp) :: StrainVec(6), StressVec(6)
    REAL(KIND=dp) :: PrevSolVec(m*nd), PrevField(13)
    REAL(KIND=dp) :: QBlock(3,3), Q(m*nd,m*nd), RotMat(3,3)
    REAL(KIND=dp) :: CMat(4,4), GMat(2,2)
    REAL(KIND=dp) :: A11, A22, SqrtDetA, A1, A2, B11, B22
    REAL(KIND=dp) :: C111, C112, C221, C222, C211, C212
    REAL(KIND=dp) :: abasis1(3), abasis2(3), abasis3(3)
    REAL(KIND=dp) :: abasis1New(3), abasis2New(3), abasis3New(3), NewDetA
    REAL(KIND=dp) :: dual1(3), dual2(3)
    REAL(KIND=dp) :: v1, v2, v3
    REAL(KIND=dp) :: yk, y1, y2, XGlob, YGlob, ZGlob
    REAL(KIND=dp) :: uq, vq, sq
    REAL(KIND=dp) :: BGBasis(n), GBasis(MaxPatchNodes), PBasis(nd)
    REAL(KIND=dp) :: StrainBasis(4,3)   ! Four rows large enough for p=1
    REAL(KIND=dp) :: MembraneStrainBasis1(4,3), MembraneStrainBasis2(4,3)

    REAL(KIND=dp) :: DOFsTransform(3,4)
    REAL(KIND=dp) :: ShearParMat(2,2,nd)
    REAL(KIND=dp) :: ChristoffelMat1(2,2,nd), ChristoffelMat2(2,2,nd)
    REAL(KIND=dp) :: BParMat(2,2,nd), BParMat1(2,2,nd), BParMat2(2,2,nd)
    REAL(KIND=dp) :: PoissonRatio(n), YoungsMod(n), ShellThickness(n), Load(n), rho(n), rho0
    REAL(KIND=dp) :: Damping(n), DampCoef
    REAL(KIND=dp) :: nu, E, h, NormalTraction, KappaS, KappaM, KappaStr
    REAL(KIND=dp) :: DetJ, Weight, Norm

    SAVE Element, GElement, Nodes, PNodes, PRefNodes
!------------------------------------------------------------------------------
    IF (m /= 6) CALL Fatal('ShellLocalMatrix', 'Wrong number of unknown fields')
    Pversion = IsActivePElement(BGElement)

    IF (PVersion )THEN
        IF(BGElement % PDefs % P > 1) CALL Fatal('ShellLocalMatrix', &
         'Set Cartesian Formulation = True to use p-elements with p > 1')
    END IF
    Family = GetElementFamily(BGElement)

    ! ------------------------------------------------------------------------------
    ! Retrieve the data which have been saved as elementwise properties:
    ! ------------------------------------------------------------------------------
    CALL RetrieveLocalFrame(BGElement, TaylorParams, PatchNodes, e1, e2, e3, &
        o, PlateBody, SphericalSurface, GElement)

    ! ------------------------------------------------------------------------------
    ! Decide what strain reduction strategy is applied and set parameters that
    ! control the selection of variational crimes.
    ! ------------------------------------------------------------------------------
    MembraneReductionMethod = MembraneStrainReductionMethod
    UseBubbles = .FALSE.
    CALL SetStrainReductionParameters(BGElement, MembraneReductionMethod, PlateBody, &
      MembraneStrainDim, UseBubbles, UseShearCorrection, DOFsTransform, &
      MembraneStrains = .TRUE.)
    
    ReductionMethod = StrainReductionMethod
    UseBubbles = Bubbles .AND. (.NOT. LargeDeflection)
    CALL SetStrainReductionParameters(BGElement, ReductionMethod, PlateBody, &
      ReducedStrainDim, UseBubbles, UseShearCorrection, DOFsTransform, &
      MembraneStrains = .FALSE.)

    ! ------------------------------------------------------------------------------
    ! The DOFs count: Currently, FE bubbles can be employed only in a very special 
    ! way by augmenting the two rotation components
    ! ------------------------------------------------------------------------------
    IF (UseBubbles) THEN
      BubbleDOFs = 2
    ELSE
      BubbleDOFs = 0
    END IF
    DOFs = m*nd ! The local stiffness matrix size after static condensation

    ! ------------------------------------------------------------------------------
    ! A general remark:
    !
    ! The number of background element nodes need not define the order of spatial 
    ! discretization if additional DOFs have been introduced with "Element" keyword. 
    ! Note thus the following concepts associated with the code variables
    !  
    !   * BGElement: a background element related to the generation of surface model
    !                (this may also be a p-element to use p-basis functions)
    !   * Element: the Lagrange interpolation element corresponding to the "Element" keyword 
    !   * GElement: an element structure corresponding to the surface reconstruction
    ! ------------------------------------------------------------------------------
    SELECT CASE(Family)
    CASE(3)
      SecondOrder = BGElement % Type % NumberOfNodes == 6
    CASE(4)
      SecondOrder = BGElement % Type % NumberOfNodes == 9
    END SELECT

    ! ------------------------------------------------------------------------------
    ! Allocate a Lagrange interpolation element structure corresponding to the
    ! "Element" keyword. A node variable suitable for defining the isoparametric 
    ! element map from the reference element onto the set which is the domain of 
    ! lines of curvature coordinates is also created. In addition, the element 
    ! structure corresponding to the surface reconstruction is created. 
    ! --------------------------------------------------------------------------
    CALL CreateLagrangeElementStructures(BGElement, nd, Element, Nodes, PNodes, &
        GElement)

    ! --------------------------------------------------------------------------
    ! Update the coordinate values of the Lagrange nodes variable. If p-basis is 
    ! used for approximating the shell equations, create also the isoparametric 
    ! geometry representation in terms of p-basis. Note that the implementation of 
    ! the p-version is not fully optimal since loads and BCs may be interpolated 
    ! only by using the node coordinates of the background element.
    ! --------------------------------------------------------------------------
    CALL WriteElementNodesVariables(BGElement, GElement, Element, Nodes, &
        PNodes, PatchNodes)

    ! --------------------------------------------------------------------------
    ! Body forces, material parameters and the shell thickness:
    ! --------------------------------------------------------------------------
    Material => GetMaterial()
    PoissonRatio(1:n) = GetReal(Material, 'Poisson Ratio')
    YoungsMod(1:n) = GetReal(Material, 'Youngs Modulus')
    ShellThickness(1:n) = GetReal(Material, 'Shell Thickness', Found)
    IF (.NOT. Found) THEN
      ShellThickness(1:n) = GetReal(Material, 'Thickness', Found)
      IF (.NOT. Found) CALL Fatal('ShellLocalMatrix', 'Shell Thickness/Thickness undefined')
    END IF

    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) THEN
      Load(1:n) = GetReal(BodyForce, 'Normal Pressure', Found)
    ELSE
      Load(1:n) = 0.0d0
    END IF
    IF ( MassAssembly ) THEN
      rho(1:n) = GetReal(Material, 'Density')
      Damping(1:n) = GetReal(Material, 'Rayleigh Damping Alpha', Found)
    END IF

    ! ------------------------------------------------------------------------------
    ! The size of the constitutive matrix for 2D shell equations
    ! ------------------------------------------------------------------------------
    csize = 4

    ! --------------------------------------------------------------------------
    ! Strain reduction operators will be applied to fields Cu where
    ! C is a 2X2-matrix of shell model parameters and u a 2-component field.
    ! In the following we create descriptions of all necessary matrices C so
    ! that the components of C can be evaluated via the Lagrange interpolation.
    !
    ! TO DO: Should we create nodal vectors of all model parameters to avoid
    ! a later call of SurfaceBasisVectors?
    ! -------------------------------------------------------------------------- 
    IF (ReductionMethod /= NoStrainReduction .OR. MembraneReductionMethod /= NoStrainReduction) THEN
      ShearParMat = 0.0d0
      BParMat = 0.0d0
      BParMat1 = 0.0d0
      BParMat2 = 0.0d0
      DO j=1,nd
        CALL SurfaceBasisVectors(Nodes % x(j), Nodes % y(j), TaylorParams, e1, e2, e3, o, &
            abasis1, abasis2, abasis3, A11, A22, SqrtDetA, B11, B22, C111, C112, C221, C222, &
            C211, C212, XGlob=XGlob, YGlob=YGlob, ZGlob=ZGlob, PlanarPoint=PlateBody, &
            Umbilical=SphericalSurface)

        ChristoffelMat1(1,1,j) = C111
        ChristoffelMat1(2,1,j) = C211
        ChristoffelMat1(1,2,j) = C112
        ChristoffelMat1(2,2,j) = C212

        ChristoffelMat2(1,1,j) = C211
        ChristoffelMat2(2,1,j) = C221
        ChristoffelMat2(1,2,j) = C212
        ChristoffelMat2(2,2,j) = C222

        BParMat(1,1,j) = B11
        BParMat(2,2,j) = B22       
        BParMat1(1,1,j) = B11
        BParMat2(2,2,j) = B22

        ShearParMat(1,1,j) = B11/a11
        ShearParMat(2,2,j) = B22/a22
      END DO
    END IF

    ! ------------------------------------------------------------------------
    ! Vectorize the previous solution and transform it to the local DOFs:
    ! ------------------------------------------------------------------------
    DO k=1,m
      PrevSolVec(k:DOFs:m) = LocalSol(k,1:nd)
    END DO

    Q = 0.0d0
    DO j=1,nd
      ! ------------------------------------------------------------------------
      ! The following transformation is designed for the Lagrange element DOFs.
      ! This would be obscure and most likely inconsistent for the p-element DOFs
      ! with p>1.
      ! ------------------------------------------------------------------------
      y1 = Nodes % x(j)
      y2 = Nodes % y(j)

      CALL SurfaceBasisVectors(y1, y2, TaylorParams, e1, e2, e3, o, abasis1, &
          abasis2, abasis3, A11, A22, SqrtDetA, B11, B22, C111, C112, C221, C222, &
          C211, C212, PlanarPoint=PlateBody, Umbilical=SphericalSurface)

      QBlock(1,1:3) = abasis1(1:3)
      QBlock(2,1:3) = abasis2(1:3)
      QBlock(3,1:3) = abasis3(1:3)         

      i0 = (j-1)*m

      Q(i0+1:i0+3,i0+1:i0+3) =  QBlock(1:3,1:3)
      !
      ! Optionally we can switch to rotated components theta such that
      ! -Du[d] = d x theta + <theta,d>d. The tangent plane components are
      ! then more intuitive when thinking in terms of moments.
      ! 
      IF (RotateDOFs .OR. DrillingDOFs) THEN
        ! 
        ! Create a matrix RotMat such that d x v = RotMat * v
        !
        RotMat = 0.0d0
        RotMat(2,1) = abasis3(3)
        RotMat(3,1) = -abasis3(2)
        RotMat(1,2) = -abasis3(3)
        RotMat(3,2) = abasis3(1)
        RotMat(1,3) = abasis3(2)
        RotMat(2,3) = -abasis3(1)

        ! Transformation from the Cartesian components of theta to
        ! the local representation in terms of the local surface basis
        DO k=1,3
          Q(i0+4,i0+3+k) = DOT_PRODUCT(RotMat(:,k), abasis1(:))
          Q(i0+5,i0+3+k) = DOT_PRODUCT(RotMat(:,k), abasis2(:))
          Q(i0+6,i0+3+k) = abasis3(k)
        END DO
      ELSE
        Q(i0+4:i0+6,i0+4:i0+6) =  QBlock(1:3,1:3)
      END IF
    END DO    
    PrevSolVec(1:DOFs) = MATMUL(Q(1:DOFs,1:DOFs),PrevSolVec(1:DOFs))

    ! ------------------------------------------------------------------------
    ! Finally, integrate local element matrices:
    ! ------------------------------------------------------------------------
    Mass = 0.0d0
    Damp = 0.0d0
    Stiff = 0.0d0
    Force = 0.0d0
    RHSForce = 0.0d0

    IP = GaussPoints( BGElement )

    QUADRATURELOOP: DO t=1,IP % n

      BM = 0.0d0
      BB = 0.0d0
      BS = 0.0d0

      NonlinBM = 0.0d0
      NonlinBS = 0.0d0
      NonlinBB = 0.0d0

      IF ( PVersion .OR. SecondOrder) THEN
        IF (PVersion) THEN
          stat = ElementInfo(BGElement, PNodes, IP % U(t), IP % V(t), IP % W(t), detJ, Basis, dBasis)
          y1 = SUM( PNodes % x(1:nd) * Basis(1:nd) )
          y2 = SUM( PNodes % y(1:nd) * Basis(1:nd) )
        ELSE
          stat = ElementInfo(Element, Nodes, IP % U(t), IP % V(t), IP % W(t), detJ, Basis, dBasis)
          y1 = SUM( Nodes % x(1:nd) * Basis(1:nd) )
          y2 = SUM( Nodes % y(1:nd) * Basis(1:nd) )
        END IF
        sq = IP % S(t)
      ELSE
        ! ---------------------------------------------------------------------------------
        ! ReductionOperatorInfo should give all necessary basis functions without the
        ! standard ElementInfo call as
        !
        !  stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
        !      detJ, Basis, dBasis )
        !
        ! Now ReductionOperatorInfo assumes p-reference elements, so
        ! switch to the reference p-element.
        ! ---------------------------------------------------------------------------------
        IF (Family==3) THEN
          uq = -1.0d0 + 2.0d0 * IP % U(t) + IP % V(t)
          vq = SQRT(3.0d0) * IP % V(t)
          sq = SQRT(3.0d0) * 2.0d0 * IP % S(t)
        ELSE
          uq = IP % U(t)
          vq = IP % V(t)
          sq = IP % S(t)
        END IF

        stat = ReductionOperatorInfo(Element, Nodes, uq, vq, StrainBasis, &
            ReductionMethod, ApplyPiolaTransform = .TRUE., detF=detJ, Basis=Basis, &
            dBasis=dBasis, Bubbles=UseBubbles)

        IF (MembraneReductionMethod /= ReductionMethod) THEN
          stat = ReductionOperatorInfo(Element, Nodes, uq, vq, MembraneStrainBasis1, &
              MembraneReductionMethod, ApplyPiolaTransform = .TRUE., EdgeDirection=1)
          IF (MembraneReductionMethod == CurlKernelWithEdgeDOFs) THEN
            ! Two sets of basis functions available, create both:
            stat = ReductionOperatorInfo(Element, Nodes, uq, vq, MembraneStrainBasis2, &
                MembraneReductionMethod, ApplyPiolaTransform = .TRUE., EdgeDirection=2)
          ELSE
            MembraneStrainBasis2 = MembraneStrainBasis1
          END IF
        ELSE
          MembraneStrainBasis1 = StrainBasis
          MembraneStrainBasis2 = StrainBasis
        END IF

        y1 = SUM( Nodes % x(1:nd) * Basis(1:nd) )
        y2 = SUM( Nodes % y(1:nd) * Basis(1:nd) )
      END IF

      ! ------------------------------------------------------------------------------
      ! The fundamental forms and Christoffel symbols at the point (y1,y2) of the
      ! principal coordinate patch:
      ! ------------------------------------------------------------------------------
      CALL SurfaceBasisVectors(y1, y2, TaylorParams, e1, e2, e3, o, abasis1, &
          abasis2, abasis3, A11, A22, SqrtDetA, B11, B22, C111, C112, C221, C222, &
          C211, C212, dual1=dual1, dual2=dual2, XGlob=XGlob, YGlob=YGlob, ZGlob=ZGlob, &
          PlanarPoint=PlateBody, Umbilical=SphericalSurface)

      ! The geometric Lame parameters:
      ! ------------------------------
      A1 = SQRT(a11)
      A2 = SQRT(a22)

      ! ------------------------------------------------
      ! Data interpolation:
      ! ------------------------------------------------
      h = SUM( ShellThickness(1:n) * Basis(1:n) )
      nu = SUM( PoissonRatio(1:n) * Basis(1:n) )
      E = SUM( YoungsMod(1:n) * Basis(1:n) )
      NormalTraction = SUM( Load(1:n) * Basis(1:n) )
      IF ( MassAssembly ) THEN
        rho0 = SUM( rho(1:n) * Basis(1:n) )
        DampCoef = SUM( Damping(1:n) * Basis(1:n) )
      END IF

      ! The matrix description of the elasticity tensor:
      CALL ElasticityMatrix(CMat, GMat, A1, A2, E, nu, DrillingDOFs, DrillingPar)

      ! Shear correction factor:
      IF ( UseShearCorrection ) THEN
        CALL ShearCorrectionFactor(KappaS, h, Nodes % x(1:Family), Nodes % y(1:Family), Family, ShearAlpha)
        !
        ! Now shear relaxation also triggers membrane/stretch relaxation provided an alpha parameter is given:
        !
        CALL ShearCorrectionFactor(KappaM, h, Nodes % x(1:Family), Nodes % y(1:Family), Family, MembraneAlpha)
        CALL ShearCorrectionFactor(KappaStr, h, Nodes % x(1:Family), Nodes % y(1:Family), Family, StretchAlpha)        
      ELSE
        KappaS = 1.0d0
        KappaM = 1.0d0
        KappaStr = 1.0d0
      END IF

      ! Relaxation of the term acting like a penalty term for the normal stress:
      IF (.NOT. DrillingDOFs) CMat(4,4) = KappaStr/KappaM * CMat(4,4)

      !-----------------------------------------------------------------------------------
      ! THE PART CORRESPONDING TO THE MEMBRANE STRAINS:
      !-----------------------------------------------------------------------------------
      ! Create first the representation of the differential DE_0(U)[V] of the linearized 
      ! membrane strain E_0(U) in the matrix form as DE_0(U)[V] = E_0(V) = BM * V (here 
      ! DE_0(U)[V] = E_0(V) holds for all U since E_0(U) is linear with respect to U).
      !------------------------------------------------------------------------------------
      Weight = h * SqrtDetA * detJ * sq

      IF ( (MembraneReductionMethod /= NoStrainReduction) .AND. (.NOT. PlateBody) ) THEN
        ! -------------------------------------------------------------------------------
        ! Apply strain reduction to terms having the Christoffel symbols as coefficients:
        ! -------------------------------------------------------------------------------
        IF (MembraneReductionMethod == DoubleReduction) THEN
          ! First, get DOFs for the RT interpolant ...
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              Family, nd, MITC, ChristoffelMat1)
          ! and then apply a second round of reductions:
          ReductionDOFsArray(1:MembraneStrainDim,1:2*nd) = MATMUL(DOFsTransform(1:MembraneStrainDim,1:Family), &
              ReductionDOFsArray(1:Family,1:2*nd))
        ELSE
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              MembraneStrainDim, nd, MembraneReductionMethod, ChristoffelMat1, EdgeDirection=1)
        END IF
        DO p=1,nd
          DO j=1,MembraneStrainDim
            BM(1,(p-1)*m+1) = BM(1,(p-1)*m+1) - ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,1)
            BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) - ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,2)
            BM(1,(p-1)*m+2) = BM(1,(p-1)*m+2) - ReductionDOFsArray(j,2*p) * MembraneStrainBasis1(j,1)
            BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) - ReductionDOFsArray(j,2*p) * MembraneStrainBasis1(j,2)
          END DO
        END DO

        IF (MembraneReductionMethod == DoubleReduction) THEN
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              Family, nd, MITC, ChristoffelMat2)
          ReductionDOFsArray(1:MembraneStrainDim,1:2*nd) = MATMUL(DOFsTransform(1:MembraneStrainDim,1:Family), &
              ReductionDOFsArray(1:Family,1:2*nd))
        ELSE
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              MembraneStrainDim, nd, MembraneReductionMethod, ChristoffelMat2, EdgeDirection=2)
        END IF
        DO p=1,nd
          DO j=1,MembraneStrainDim
            BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) - ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis2(j,1)
            BM(2,(p-1)*m+1) = BM(2,(p-1)*m+1) - ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis2(j,2)
            BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) - ReductionDOFsArray(j,2*p) * MembraneStrainBasis2(j,1)
            BM(2,(p-1)*m+2) = BM(2,(p-1)*m+2) - ReductionDOFsArray(j,2*p) * MembraneStrainBasis2(j,2)
          END DO
        END DO

        ! -------------------------------------------------------------------------------
        ! The terms having the components of the second fundamental form as coefficients:
        ! -------------------------------------------------------------------------------
        IF (MembraneReductionMethod == CurlKernelWithEdgeDOFs) THEN
          ! This splits up into two steps ...
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              MembraneStrainDim, nd, MembraneReductionMethod, BParMat1, EdgeDirection=1)
          DO p=1,nd
            DO j=1,MembraneStrainDim
              BM(1,(p-1)*m+3) = BM(1,(p-1)*m+3) - ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,1)
              BM(3,(p-1)*m+3) = BM(3,(p-1)*m+3) - ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,2)
            END DO
          END DO
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              MembraneStrainDim, nd, MembraneReductionMethod, BParMat2, EdgeDirection=2)
          DO p=1,nd
            DO j=1,MembraneStrainDim
              BM(3,(p-1)*m+3) = BM(3,(p-1)*m+3) - ReductionDOFsArray(j,2*p) * MembraneStrainBasis2(j,1)
              BM(2,(p-1)*m+3) = BM(2,(p-1)*m+3) - ReductionDOFsArray(j,2*p) * MembraneStrainBasis2(j,2)
            END DO
          END DO
        ELSE
          ! Here all strain reduction operations are created by using just one parameter matrix:
          IF (MembraneReductionMethod == DoubleReduction) THEN
            CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
                Family, nd, MITC, BParMat)
            ReductionDOFsArray(1:MembraneStrainDim,1:2*nd) = MATMUL(DOFsTransform(1:MembraneStrainDim,1:Family), &
                ReductionDOFsArray(1:Family,1:2*nd))
          ELSE
            CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
                MembraneStrainDim, nd, MembraneReductionMethod, BParMat)
          END IF
          DO p=1,nd
            DO j=1,MembraneStrainDim
              BM(1,(p-1)*m+3) = BM(1,(p-1)*m+3) - ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,1)
              BM(3,(p-1)*m+3) = BM(3,(p-1)*m+3) - ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,2)
              BM(3,(p-1)*m+3) = BM(3,(p-1)*m+3) - ReductionDOFsArray(j,2*p) * MembraneStrainBasis1(j,1)
              BM(2,(p-1)*m+3) = BM(2,(p-1)*m+3) - ReductionDOFsArray(j,2*p) * MembraneStrainBasis1(j,2)
            END DO
          END DO
        END IF

        ! -------------------------------------------------------------------------------
        ! The partial derivative terms
        ! -------------------------------------------------------------------------------
        IF (.TRUE.) THEN
          ! No modifications:
          DO p=1,nd
            BM(1,(p-1)*m+1) = BM(1,(p-1)*m+1) + dBasis(p,1)
            BM(2,(p-1)*m+2) = BM(2,(p-1)*m+2) + dBasis(p,2)
            BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) + dBasis(p,2)
            BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) + dBasis(p,1)
          END DO
        ELSE
          ! This branch is now for testing only and the full support for all strain reduction
          ! operators is missing.
          ! Apply the strain reduction operator to plain derivative terms. This doesn't make
          ! any modification really with the current set of strain reduction operators.

          IF (MembraneReductionMethod == CurlKernelWithEdgeDOFs) THEN
            CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
                MembraneStrainDim, nd, MembraneReductionMethod, GradientField=.TRUE., EdgeDirection=1)
            DO p=1,nd
              DO j=1,MembraneStrainDim
                BM(1,(p-1)*m+1) = BM(1,(p-1)*m+1) + ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,1)
                BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) + ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,2)
                BM(1,(p-1)*m+1) = BM(1,(p-1)*m+1) + ReductionDOFsArray(j,2*p) * MembraneStrainBasis1(j,1)
                BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) + ReductionDOFsArray(j,2*p) * MembraneStrainBasis1(j,2)
              END DO
            END DO
            CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
                MembraneStrainDim, nd, MembraneReductionMethod, GradientField=.TRUE., EdgeDirection=2)
            DO p=1,nd
              DO j=1,MembraneStrainDim
                BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) + ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis2(j,1)
                BM(2,(p-1)*m+2) = BM(2,(p-1)*m+2) + ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis2(j,2)
                BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) + ReductionDOFsArray(j,2*p) * MembraneStrainBasis2(j,1)
                BM(2,(p-1)*m+2) = BM(2,(p-1)*m+2) + ReductionDOFsArray(j,2*p) * MembraneStrainBasis2(j,2)
              END DO
            END DO
          ELSE
            CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
                MembraneStrainDim, nd, MembraneReductionMethod, GradientField=.TRUE.)
            DO p=1,nd
              DO j=1,MembraneStrainDim
                BM(1,(p-1)*m+1) = BM(1,(p-1)*m+1) + ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,1)
                BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) + ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,2)
                BM(1,(p-1)*m+1) = BM(1,(p-1)*m+1) + ReductionDOFsArray(j,2*p) * MembraneStrainBasis1(j,1)
                BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) + ReductionDOFsArray(j,2*p) * MembraneStrainBasis1(j,2)

                BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) + ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,1)
                BM(2,(p-1)*m+2) = BM(2,(p-1)*m+2) + ReductionDOFsArray(j,2*p-1) * MembraneStrainBasis1(j,2)
                BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) + ReductionDOFsArray(j,2*p) * MembraneStrainBasis1(j,1)
                BM(2,(p-1)*m+2) = BM(2,(p-1)*m+2) + ReductionDOFsArray(j,2*p) * MembraneStrainBasis1(j,2)
              END DO
            END DO
          END IF
        END IF

      ELSE
        !-------------------------------------------------
        ! Use standard weak formulation:
        !-------------------------------------------------       
        DO p=1,nd
          BM(1,(p-1)*m+1) = dBasis(p,1) - C111 * Basis(p)
          BM(1,(p-1)*m+2) = -C112 * Basis(p)
          BM(1,(p-1)*m+3) = -B11 * Basis(p)

          BM(2,(p-1)*m+1) = -C221 * Basis(p)
          BM(2,(p-1)*m+2) = dBasis(p,2) - C222 * Basis(p)
          BM(2,(p-1)*m+3) = -B22 * Basis(p) 

          BM(3,(p-1)*m+1) = dBasis(p,2) - 2.0d0 * C211 * Basis(p)
          BM(3,(p-1)*m+2) = dBasis(p,1) - 2.0d0 * C212 * Basis(p)
        END DO
      END IF

      IF (DrillingDOFs) THEN
        !----------------------------------------------------------------------
        ! Add terms which define the drilling DOFs:
        !----------------------------------------------------------------------
        BM(4,6:DOFs:m) = Basis(1:nd)
        DO p=1,nd
          BM(4,(p-1)*m+2) = -0.5d0 * dBasis(p,1)
          BM(4,(p-1)*m+1) = 0.5d0 * dBasis(p,2)
        END DO
      ELSE
        !----------------------------------------------------------------------
        ! Normal stress T^{33} via energy principle: We add a term of the type
        ! e * T^{33}(e)
        !----------------------------------------------------------------------
        BM(4,6:DOFs:m) = -Basis(1:nd)
        BM(4,1:DOFs) = BM(4,1:DOFs) + nu/((1.0d0-nu)*A1**2) * BM(1,1:DOFs) + &
            nu/((1.0d0-nu)*A2**2) * BM(2,1:DOFs)
      END IF

      StrainVec = 0.0d0
      IF (LargeDeflection) THEN
        ! ---------------------------------------------------------------------------------------
        ! The differential DE(U)[V] of the membrane strain E(U) is by definition linear with 
        ! respect to V and thus have a matrix representation DE(U)[V] ~ BM * V + NonlinBM(U) * V.
        ! The matrix BM is already created and here we create the matrix NonlinBM(U), which
        ! depends on the current solution iterate U.
        ! ---------------------------------------------------------------------------------------
        ! Strain component 11:
        ! ---------------------------------------------------------------------------------------

        ! The part depending on the covariant derivative v_{2|1}:
        BWork = 0.0d0
        BWork(1,1:DOFs:m) = -C211 * Basis(1:nd)
        BWork(1,2:DOFs:m) = dBasis(1:nd,1) - C212 * Basis(1:nd)
        PrevField(1) = SUM( BWork(1,1:DOFs) * PrevSolVec(1:DOFs) )

        ! The part depending on v_{3|1} + B11/A11 v_1:
        BWork(2,1:DOFs:m) = B11/A11 * Basis(1:nd)
        BWork(2,3:DOFs:m) = dBasis(1:nd,1)
        PrevField(2) = SUM( BWork(2,1:DOFs) * PrevSolVec(1:DOFs) )

        NonlinBM(1,1:DOFs) = SUM( BM(1,1:DOFs) * PrevSolVec(1:DOFs) ) / A1**2 * BM(1,1:DOFs) + &
            PrevField(1) / A2**2 * BWork(1,1:DOFs) + &
            PrevField(2) * BWork(2,1:DOFs)

        ! The nonlinear part of strain component 11 for the current iterate:
        StrainVec(1) = 0.5_dp * SUM(NonlinBM(1,1:DOFs) * PrevSolVec(1:DOFs))

        ! ---------------------------------------------------------------------------------------
        ! Strain component 22:
        ! ---------------------------------------------------------------------------------------

        ! The part depending on the covariant derivative v_{1|2}:
        BWork(3,1:DOFs:m) = dBasis(1:nd,2) - C211 * Basis(1:nd)
        BWork(3,2:DOFs:m) = -C212 * Basis(1:nd)
        PrevField(3) = SUM( BWork(3,1:DOFs) * PrevSolVec(1:DOFs) )

        ! The part depending on v_{3|2} + B22/A22 v_2:
        BWork(4,2:DOFs:m) = B22/A22 * Basis(1:nd)
        BWork(4,3:DOFs:m) = dBasis(1:nd,2)
        PrevField(4) = SUM( BWork(4,1:DOFs) * PrevSolVec(1:DOFs) )

        NonlinBM(2,1:DOFs) = SUM( BM(2,1:DOFs) * PrevSolVec(1:DOFs) ) / A2**2 * BM(2,1:DOFs) + &
            PrevField(3) / A1**2 * BWork(3,1:DOFs) + &
            PrevField(4) * BWork(4,1:DOFs)

        ! The nonlinear part of strain component 22 for the current iterate:
        StrainVec(2) = 0.5_dp * SUM(NonlinBM(2,1:DOFs) * PrevSolVec(1:DOFs))

        ! ---------------------------------------------------------------------------------------
        ! Strain component 12:
        ! ---------------------------------------------------------------------------------------
        NonlinBM(3,1:DOFs) =  BM(1,1:DOFs) * PrevField(3) / A1**2 + &
            SUM( BM(1,1:DOFs) * PrevSolVec(1:DOFs) ) / A1**2 * BWork(3,1:DOFs) + &
            BM(2,1:DOFs) * PrevField(1) / A2**2 + &
            SUM( BM(2,1:DOFs) * PrevSolVec(1:DOFs) ) / A2**2 * BWork(1,1:DOFs) + &
            BWork(2,1:DOFs) * PrevField(4) + PrevField(2) * BWork(4,1:DOFs) 

        ! The nonlinear part of strain component 12 for the current iterate:
        StrainVec(3) = 0.5_dp * SUM(NonlinBM(3,1:DOFs) * PrevSolVec(1:DOFs))

        ! ---------------------------------------------------------------------------------------
        ! A strain-like variable e such that the normal stress T33 = DW/De
        ! ---------------------------------------------------------------------------------------
        BWork(5,4:DOFs:m) = Basis(1:nd)
        BWork(6,5:DOFs:m) = Basis(1:nd)
        BWork(7,6:DOFs:m) = Basis(1:nd)
        PrevField(5:7) = MATMUL(BWork(5:7,1:DOFs),PrevSolVec(1:DOFs))

        NonlinBM(4,1:DOFs) = nu/(1.0d0-nu) * NonlinBM(1,1:DOFs) / A1**2 + &
            nu/(1.0d0-nu) * NonlinBM(2,1:DOFs) / A2**2 + &
            BWork(5,1:DOFs) * PrevField(5) / A1**2 + &
            BWork(6,1:DOFs) * PrevField(6) / A2**2 +  BWork(7,1:DOFs) * PrevField(7)

        ! The nonlinear part of e for the current iterate: 
        StrainVec(4) = 0.5_dp * SUM(NonlinBM(4,1:DOFs) * PrevSolVec(1:DOFs))
      END IF

      
      CALL StrainEnergyDensity(Stiff, CMat, BM + NonlinBM, csize, DOFs, KappaM*Weight)
      
      ! The linear part of strain for the current iterate:
      StrainVec(1:csize) = StrainVec(1:csize) + MATMUL( BM(1:csize,1:DOFs), PrevSolVec(1:DOFs) )

      ! Residual terms for RHS:
      StressVec(1:csize) = MATMUL( CMat(1:csize,1:csize), StrainVec(1:csize) )
      Force(1:DOFs) = Force(1:DOFs) - MATMUL( TRANSPOSE(BM(1:csize,1:DOFs) + NonlinBM(1:csize,1:DOFs)), &
         StressVec(1:csize) ) * Weight * KappaM

      ! The remaining terms for the complete Newton iteration:
      IF (LargeDeflection) THEN
        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + &
            MATMUL( TRANSPOSE(BM(1:1,1:DOFs)),BM(1:1,1:DOFs))/A1**2 * &
            (StressVec(1) + nu/(1.0d0-nu) * StressVec(4)/A1**2) * KappaM * Weight + &
            MATMUL( TRANSPOSE(BWork(1:1,1:DOFs)),BWork(1:1,1:DOFs))/A2**2 * &
            (StressVec(1) + nu/(1.0d0-nu) * StressVec(4)/A1**2 ) * KappaM* Weight + &
            MATMUL( TRANSPOSE(BWork(2:2,1:DOFs)),BWork(2:2,1:DOFs)) * &
            (StressVec(1) + nu/(1.0d0-nu) * StressVec(4)/A1**2) * KappaM * Weight
 
        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + &
            MATMUL( TRANSPOSE(BM(2:2,1:DOFs)),BM(2:2,1:DOFs))/A2**2 * &
            (StressVec(2) + nu/(1.0d0-nu) * StressVec(4)/A2**2) * KappaM * Weight + &
            MATMUL( TRANSPOSE(BWork(3:3,1:DOFs)),BWork(3:3,1:DOFs))/A1**2 * &
            (StressVec(2) + nu/(1.0d0-nu) * StressVec(4)/A2**2) * KappaM * Weight + &
            MATMUL( TRANSPOSE(BWork(4:4,1:DOFs)),BWork(4:4,1:DOFs)) * &
            (StressVec(2) + nu/(1.0d0-nu) * StressVec(4)/A2**2) * KappaM * Weight

        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + ( &
            MATMUL( TRANSPOSE(BWork(3:3,1:DOFs)),BM(1:1,1:DOFs)) + &
            MATMUL( TRANSPOSE(BM(1:1,1:DOFs)),BWork(3:3,1:DOFs)) ) * StressVec(3)/A1**2 * KappaM * Weight + ( &
            MATMUL( TRANSPOSE(BWork(1:1,1:DOFs)),BM(2:2,1:DOFs)) + &
            MATMUL( TRANSPOSE(BM(2:2,1:DOFs)),BWork(1:1,1:DOFs)) ) * StressVec(3)/A2**2 * KappaM * Weight + ( &
            MATMUL( TRANSPOSE(BWork(2:2,1:DOFs)),BWork(4:4,1:DOFs)) + &
            MATMUL( TRANSPOSE(BWork(4:4,1:DOFs)),BWork(2:2,1:DOFs)) ) * StressVec(3) * KappaM * Weight

        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + &
            MATMUL(TRANSPOSE(BWork(5:5,1:DOFs)),BWork(5:5,1:DOFs)) * StressVec(4) / A1**2 * KappaM * Weight + &
            MATMUL(TRANSPOSE(BWork(6:6,1:DOFs)),BWork(6:6,1:DOFs)) * StressVec(4) / A2**2 * KappaM * Weight + &
            MATMUL(TRANSPOSE(BWork(7:7,1:DOFs)),BWork(7:7,1:DOFs)) * StressVec(4) * KappaM * Weight
      END IF

      !-----------------------------------------------------------------------------------
      ! THE PART CORRESPONDING TO THE TRANSVERSE SHEAR STRAINS:
      !-----------------------------------------------------------------------------------
      ! Create first the representation of the differential DE_0(U)[V] of the linearized 
      ! transverse shear strain E_0(U) in the matrix form as DE_0(U)[V] = E_0(V) = BS * V
      ! (here DE_0(U)[V] = E_0(V) holds since E_0(U) is linear with respect to U).
      !------------------------------------------------------------------------------------
      IF (ReductionMethod /= NoStrainReduction) THEN
        !---------------------------------------------------------------     
        ! Get coefficients that can be used to evaluate the DOFs for the interpolant of 
        ! the rotation variable in the reduced strain space:
        !---------------------------------------------------------------
        IF (ReductionMethod == DoubleReduction) THEN
          ! First, get DOFs for the RT interpolant ...
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              3, nd, MITC)
          ! and then apply a second round of reductions:
          ReductionDOFsArray(1:2,1:2*nd) = MATMUL(DOFsTransform(1:2,1:3), ReductionDOFsArray(1:3,1:2*nd))
        ELSE
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              ReducedStrainDim, nd, ReductionMethod)
        END IF

        DO p=1,nd
          DO j=1,ReducedStrainDim
            BS(1:2,(p-1)*m+4) = BS(1:2,(p-1)*m+4) - ReductionDOFsArray(j,2*p-1) * StrainBasis(j,1:2)
            BS(1:2,(p-1)*m+5) = BS(1:2,(p-1)*m+5) - ReductionDOFsArray(j,2*p) * StrainBasis(j,1:2)          
          END DO
        END DO

        ! If rotations are augmented with bubbles, we need to compute an augmented DOFs array
        IF (UseBubbles .AND. ReductionMethod == DoubleReduction) THEN
          CALL ReductionOperatorBubbleDofs(Element, Nodes, ReductionDOFsArray(1:2,2*nd+1:2*nd+2), &
              ReducedStrainDim, 1, Family, CurlKernel)

          DO j=1,ReducedStrainDim
            BS(1:2,DOFs+1) = BS(1:2,DOFs+1) - ReductionDOFsArray(j,2*nd+1) * StrainBasis(j,1:2)
            BS(1:2,DOFs+2) = BS(1:2,DOFs+2) - ReductionDOFsArray(j,2*nd+2) * StrainBasis(j,1:2)          
          END DO
        END IF

        !------------------------------------------------------------------------------
        ! The reduction for the terms depending also on the curvature and derivatives:
        !------------------------------------------------------------------------------
        IF (ReductionMethod == DoubleReduction) THEN
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              3, nd, MITC, ShearParMat)
          ReductionDOFsArray(1:2,1:2*nd) = MATMUL(DOFsTransform(1:2,1:3), ReductionDOFsArray(1:3,1:2*nd))
        ELSE
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              ReducedStrainDim, nd, ReductionMethod, ShearParMat)
        END IF
        DO p=1,nd
          DO j=1,ReducedStrainDim
            BS(1:2,(p-1)*m+1) = BS(1:2,(p-1)*m+1) + ReductionDOFsArray(j,2*p-1) * StrainBasis(j,1:2)
            BS(1:2,(p-1)*m+2) = BS(1:2,(p-1)*m+2) + ReductionDOFsArray(j,2*p) * StrainBasis(j,1:2) 
          END DO
          BS(1:2,(p-1)*m+3) = dBasis(p,1:2)
        END DO
      ELSE
        !-------------------------------------------------
        ! Use standard weak formulation:
        !------------------------------------------------- 
        DO p=1,nd
          BS(1:2,(p-1)*m+3) = dBasis(p,1:2)
       
          BS(1,(p-1)*m+1) = B11/a11 * Basis(p) 
          BS(1,(p-1)*m+4) = -Basis(p) 
          
          BS(2,(p-1)*m+2) = B22/a22 * Basis(p) 
          BS(2,(p-1)*m+5) = -Basis(p)
        END DO

        IF (UseBubbles) THEN
          ! With rotation bubbles:
          BS(1,DOFs+1) = -Basis(nd+1)
          BS(2,DOFs+2) = -Basis(nd+1)
        END IF

      END IF

      IF (LargeDeflection) THEN
        ! ---------------------------------------------------------------------------------------
        ! The representation of the differential DE(U)[V] of the transverse shear strain E(U) in 
        ! the matrix form as DE(U)[V] ~ BS * V + NonlinBS(U) * V. The matrix BS is already created 
        ! and here we compute the matrix NonlinBS(U), which depends on the current solution iterate 
        ! U. TO DO: Bubbles are not yet handled 
        ! ---------------------------------------------------------------------------------------
        ! Strain component 13:
        ! ---------------------------------------------------------------------------------------
        NonlinBS(1,1:DOFs) = -BM(1,1:DOFs) * PrevField(5) / A1**2 - &
            SUM( BM(1,1:DOFs) * PrevSolVec(1:DOFs) ) / A1**2 * BWork(5,1:DOFs) - &
            BWork(1,1:DOFs) * PrevField(6) / A2**2 - &
            PrevField(1) / A2**2 * BWork(6,1:DOFs) - &
            BWork(2,1:DOFs) * PrevField(7) - PrevField(2) * BWork(7,1:DOFs) 

        ! The nonlinear part of strain component 13 for the current iterate:
        StrainVec(5) = 0.5_dp * SUM(NonlinBS(1,1:DOFs) * PrevSolVec(1:DOFs))
        ! ---------------------------------------------------------------------------------------
        ! Strain component 23:
        ! ---------------------------------------------------------------------------------------
        NonlinBS(2,1:DOFs) = -BWork(3,1:DOFs) * PrevField(5) / A1**2 - &
            PrevField(3) / A1**2 * BWork(5,1:DOFs) - &
            BM(2,1:DOFs) * PrevField(6) / A2**2 - &
            SUM( BM(2,1:DOFs) * PrevSolVec(1:DOFs) ) / A2**2 * BWork(6,1:DOFs) - &
            BWork(4,1:DOFs) * PrevField(7) - PrevField(4) * BWork(7,1:DOFs) 

        ! The nonlinear part of strain component 23 for the current iterate:
        StrainVec(6) = 0.5_dp * SUM(NonlinBS(2,1:DOFs) * PrevSolVec(1:DOFs))
      END IF

      CALL StrainEnergyDensity(Stiff, GMat, BS+NonlinBS, 2, DOFs+BubbleDOFs, KappaS*Weight)

      ! The linear part of strain for the current iterate:
      StrainVec(5:6) = StrainVec(5:6) + MATMUL( BS(1:2,1:DOFs), PrevSolVec(1:DOFs) )

      ! Residual terms for RHS:
      StressVec(5:6) = MATMUL( GMat(1:2,1:2), StrainVec(5:6) )
      Force(1:DOFs+BubbleDOFs) = Force(1:DOFs+BubbleDOFs) - & 
          MATMUL( TRANSPOSE(BS(1:2,1:DOFs+BubbleDOFs) + NonlinBS(1:2,1:DOFs+BubbleDOFs)), &
          StressVec(5:6) ) * KappaS * Weight

      ! The remaining terms for the complete Newton iteration:
      IF (LargeDeflection) THEN
        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) - ( &
            MATMUL( TRANSPOSE(BWork(5:5,1:DOFs)),BM(1:1,1:DOFs)) + &
            MATMUL( TRANSPOSE(BM(1:1,1:DOFs)),BWork(5:5,1:DOFs)) ) * StressVec(5)/A1**2 * KappaS * Weight - ( &
            MATMUL( TRANSPOSE(BWork(1:1,1:DOFs)),BWork(6:6,1:DOFs)) + &
            MATMUL( TRANSPOSE(BWork(6:6,1:DOFs)),BWork(1:1,1:DOFs)) ) * StressVec(5)/A2**2 * KappaS * Weight - ( &
            MATMUL( TRANSPOSE(BWork(2:2,1:DOFs)),BWork(7:7,1:DOFs)) + &
            MATMUL( TRANSPOSE(BWork(7:7,1:DOFs)),BWork(2:2,1:DOFs)) ) * StressVec(5) * KappaS * Weight

        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) - ( &
            MATMUL( TRANSPOSE(BWork(5:5,1:DOFs)),BWork(3:3,1:DOFs)) + &
            MATMUL( TRANSPOSE(BWork(3:3,1:DOFs)),BWork(5:5,1:DOFs)) ) * StressVec(6)/A1**2 * KappaS * Weight - ( &
            MATMUL( TRANSPOSE(BM(2:2,1:DOFs)),BWork(6:6,1:DOFs)) + &
            MATMUL( TRANSPOSE(BWork(6:6,1:DOFs)),BM(2:2,1:DOFs)) ) * StressVec(6)/A2**2 * KappaS * Weight - ( &
            MATMUL( TRANSPOSE(BWork(4:4,1:DOFs)),BWork(7:7,1:DOFs)) + &
            MATMUL( TRANSPOSE(BWork(7:7,1:DOFs)),BWork(4:4,1:DOFs)) ) * StressVec(6) * KappaS * Weight
      END IF
        
      !---------------------------------------------------------------
      ! THE PART CORRESPONDING TO THE BENDING STRAINS
      !---------------------------------------------------------------      
      Weight = h**3/12.0d0 * SqrtDetA * detJ * sq
      DO p=1,nd
        BB(1,(p-1)*m+4) = dBasis(p,1) - C111 * Basis(p)
        BB(1,(p-1)*m+5) = -C112 * Basis(p)

        BB(2,(p-1)*m+4) = -C221 * Basis(p)
        BB(2,(p-1)*m+5) = dBasis(p,2) - C222 * Basis(p)

        BB(3,(p-1)*m+4) = dBasis(p,2) - 2.0d0 * C211 * Basis(p)
        BB(3,(p-1)*m+5) = dBasis(p,1) - 2.0d0 * C212 * Basis(p)
        IF (BenchmarkProblem) THEN
          ! The following can be obtained by adding some (unnatural) weighted combinations of membrane strains
          ! to the bending strains:
          BB(3,(p-1)*m+1) = B11/a11 * (dBasis(p,2) - C211 * Basis(p)) - B22/a22 * C211 * Basis(p)
          BB(3,(p-1)*m+2) = B11/a11 * (- C212) * Basis(p) + B22/a22 * (dBasis(p,1) - C212 * Basis(p))
        ELSE
          BB(3,(p-1)*m+1) = B11/a11 * C211 * Basis(p) - B22/a22 * (dBasis(p,2) - C211 * Basis(p))
          BB(3,(p-1)*m+2) = -B11/a11 * (dBasis(p,1) - C212 * Basis(p)) + B22/a22 * C212 * Basis(p)
        END IF
      END DO

      IF (.NOT. DrillingDOFs) THEN
        DO p=1,nd
          BB(1,(p-1)*m+6) = -B11 * Basis(p)         
          BB(2,(p-1)*m+6) = -B22 * Basis(p)
        END DO
      END IF

      IF (UseBubbles) THEN
        ! With rotation bubbles:
        BB(1,DOFs+1) = dBasis(nd+1,1) - C111 * Basis(nd+1)
        BB(1,DOFs+2) = -C112 * Basis(nd+1)
        BB(2,DOFs+1) = -C221 * Basis(nd+1)
        BB(2,DOFs+2) = dBasis(nd+1,2) - C222 * Basis(nd+1)
        BB(3,DOFs+1) = dBasis(nd+1,2) - 2.0d0 * C211 * Basis(nd+1)
        BB(3,DOFs+2) = dBasis(nd+1,1) - 2.0d0 * C212 * Basis(nd+1)
      END IF

      IF (LargeDeflection .AND. NonlinearBending) THEN
        ! Save the parts beta_{1|1} - B11 beta_3 and beta_{2|2} - B22 beta_3:
        BWork(8,:) = BB(1,:)
        BWork(12,:) = BB(2,:)
        PrevField(8) = SUM( BWork(8,1:DOFs) * PrevSolVec(1:DOFs) )
        PrevField(12) = SUM( BWork(12,1:DOFs) * PrevSolVec(1:DOFs) )
      END IF
      BB(1,:) = BB(1,:) - B11/a11 * BM(1,:)
      BB(2,:) = BB(2,:) - B22/a22 * BM(2,:)

      StrainVec(1:3) = 0.0d0
      IF (LargeDeflection .AND. NonlinearBending) THEN
        ! ---------------------------------------------------------------------------------------
        ! Strain component 11:
        ! ---------------------------------------------------------------------------------------

        ! The part depending on the covariant derivative beta_{2|1}:
        BWork(9,4:DOFs:m) = -C211 * Basis(1:nd)
        BWork(9,5:DOFs:m) = dBasis(1:nd,1) - C212 * Basis(1:nd)
        PrevField(9) = SUM( BWork(9,1:DOFs) * PrevSolVec(1:DOFs) )

        ! The part depending on beta_{3|1} + B11/A11 beta_1:
        BWork(10,4:DOFs:m) = B11/A11 * Basis(1:nd)
        BWork(10,6:DOFs:m) = dBasis(1:nd,1)
        PrevField(10) = SUM( BWork(10,1:DOFs) * PrevSolVec(1:DOFs) )

        NonlinBB(1,1:DOFs) = BM(1,1:DOFs) * PrevField(8) / A1**2 + &
            SUM( BM(1,1:DOFs) * PrevSolVec(1:DOFs) ) / A1**2 * BWork(8,1:DOFs) + &
            BWork(1,1:DOFs) * PrevField(9) / A2**2 + &
            PrevField(1) / A2**2 * BWork(9,1:DOFs) + &
            BWork(2,1:DOFs) * PrevField(10) + PrevField(2) * BWork(10,1:DOFs)

        ! ---------------------------------------------------------------------------------------
        ! Strain component 22:
        ! ---------------------------------------------------------------------------------------

        ! The part depending on the covariant derivative beta_{1|2}:
        BWork(11,4:DOFs:m) = dBasis(1:nd,2) - C211 * Basis(1:nd)
        BWork(11,5:DOFs:m) = -C212 * Basis(1:nd)
        PrevField(11) = SUM( BWork(11,1:DOFs) * PrevSolVec(1:DOFs) )

        ! The part depending on beta_{3|2} + B22/A22 beta_2:
        BWork(13,5:DOFs:m) = B22/A22 * Basis(1:nd)
        BWork(13,6:DOFs:m) = dBasis(1:nd,2)
        PrevField(13) = SUM( BWork(13,1:DOFs) * PrevSolVec(1:DOFs) )

        NonlinBB(2,1:DOFs) = BWork(3,1:DOFs) * PrevField(11) / A1**2 + &
            PrevField(3) / A1**2 * BWork(11,1:DOFs) + &
            BM(2,1:DOFs) * PrevField(12) / A2**2 + &
            SUM( BM(2,1:DOFs) * PrevSolVec(1:DOFs) ) / A2**2 * BWork(12,1:DOFs) + &
            BWork(4,1:DOFs) * PrevField(13) + PrevField(4) * BWork(13,1:DOFs)

        ! ---------------------------------------------------------------------------------------
        ! Strain component 12:
        ! ---------------------------------------------------------------------------------------
        NonlinBB(3,1:DOFs) = BM(1,1:DOFs) * PrevField(11) / A1**2 + &
            SUM( BM(1,1:DOFs) * PrevSolVec(1:DOFs) ) / A1**2 * BWork(11,1:DOFs) + &
            BWork(8,1:DOFs) * PrevField(3) / A1**2 + PrevField(8) / A1**2 * BWork(3,1:DOFs) + &
            BWork(1,1:DOFs) * PrevField(12) / A2**2 + PrevField(1) / A2**2 * BWork(12,1:DOFs) + &
            BM(2,1:DOFs) * PrevField(9) / A2**2 + &
            SUM( BM(2,1:DOFs) * PrevSolVec(1:DOFs) ) / A2**2 * BWork(9,1:DOFs) + &
            BWork(2,1:DOFs) * PrevField(13) + PrevField(2) * BWork(13,1:DOFs) + &
            BWork(10,1:DOFs) * PrevField(4) + PrevField(10) * BWork(4,1:DOFs)

        NonlinBB(1:3,:) = NonlinBB(1:3,:) - (B11/A11 + B22/A22) * NonlinBM(1:3,:)

        StrainVec(1:3) = 0.5_dp * MATMUL(NonlinBB(1:3,1:DOFs), PrevSolVec(1:DOFs))
      END IF

      StrainVec(1:3) = StrainVec(1:3) + MATMUL( BB(1:3,1:DOFs), PrevSolVec(1:DOFs) )

      CALL StrainEnergyDensity(Stiff, CMat, BB + NonlinBB, 3, DOFs+BubbleDOFs, Weight)
      StressVec(1:3) = MATMUL( CMat(1:3,1:3), StrainVec(1:3) )

      ! Residual terms for RHS:
      Force(1:DOFs+BubbleDOFs) = Force(1:DOFs+BubbleDOFs) - &
          MATMUL( TRANSPOSE(BB(1:3,1:DOFs+BubbleDOFs) + NonlinBB(1:3,1:DOFs+BubbleDOFs)), &
          StressVec(1:3) ) * Weight

      ! The remaining terms for the complete Newton iteration:
      IF (LargeDeflection .AND. NonlinearBending) THEN
        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + ( &
            MATMUL(TRANSPOSE(BWork(8:8,1:DOFs)), BM(1:1,1:DOFs)) + &
            MATMUL(TRANSPOSE(BM(1:1,1:DOFs)), BWork(8:8,1:DOFs))) * StressVec(1) / A1**2 * Weight + ( &
            MATMUL(TRANSPOSE(BWork(9:9,1:DOFs)), BWork(1:1,1:DOFs)) + &
            MATMUL(TRANSPOSE(BWork(1:1,1:DOFs)), BWork(9:9,1:DOFs))) * StressVec(1) / A2**2 * Weight + ( &
            MATMUL(TRANSPOSE(BWork(10:10,1:DOFs)), BWork(2:2,1:DOFs)) + &
            MATMUL(TRANSPOSE(BWork(2:2,1:DOFs)), BWork(10:10,1:DOFs))) * StressVec(1) * Weight

        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + ( &
            MATMUL(TRANSPOSE(BWork(11:11,1:DOFs)), BWork(3:3,1:DOFs)) + &
            MATMUL(TRANSPOSE(BWork(3:3,1:DOFs)), BWork(11:11,1:DOFs))) * StressVec(2) / A1**2 * Weight + ( &
            MATMUL(TRANSPOSE(BWork(12:12,1:DOFs)), BM(2:2,1:DOFs)) + &
            MATMUL(TRANSPOSE(BM(2:2,1:DOFs)), BWork(12:12,1:DOFs))) * StressVec(2) / A2**2 * Weight + ( &
            MATMUL(TRANSPOSE(BWork(13:13,1:DOFs)), BWork(4:4,1:DOFs)) + &
            MATMUL(TRANSPOSE(BWork(4:4,1:DOFs)), BWork(13:13,1:DOFs))) * StressVec(2) * Weight

        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + ( &
            MATMUL(TRANSPOSE(BWork(11:11,1:DOFs)), BM(1:1,1:DOFs)) + &
            MATMUL(TRANSPOSE(BM(1:1,1:DOFs)), BWork(11:11,1:DOFs)) + &
            MATMUL(TRANSPOSE(BWork(3:3,1:DOFs)), BWork(8:8,1:DOFs)) + &
            MATMUL(TRANSPOSE(BWork(8:8,1:DOFs)), BWork(3:3,1:DOFs))) * StressVec(3) / A1**2 * Weight + ( &
            MATMUL(TRANSPOSE(BWork(12:12,1:DOFs)), BWork(1:1,1:DOFs)) + &
            MATMUL(TRANSPOSE(BWork(1:1,1:DOFs)), BWork(12:12,1:DOFs)) + &
            MATMUL(TRANSPOSE(BM(2:2,1:DOFs)), BWork(9:9,1:DOFs)) + &
            MATMUL(TRANSPOSE(BWork(9:9,1:DOFs)), BM(2:2,1:DOFs)) ) * StressVec(3) / A2**2 * Weight + ( &

            MATMUL(TRANSPOSE(BWork(13:13,1:DOFs)), BWork(2:2,1:DOFs)) + &
            MATMUL(TRANSPOSE(BWork(2:2,1:DOFs)), BWork(13:13,1:DOFs)) + &
            MATMUL(TRANSPOSE(BWork(4:4,1:DOFs)), BWork(10:10,1:DOFs)) + &
            MATMUL(TRANSPOSE(BWork(10:10,1:DOFs)), BWork(4:4,1:DOFs))) * StressVec(3) * Weight

        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) - ( &
            MATMUL( TRANSPOSE(BM(1:1,1:DOFs)),BM(1:1,1:DOFs))/A1**2 + &
            MATMUL( TRANSPOSE(BWork(1:1,1:DOFs)),BWork(1:1,1:DOFs))/A2**2 + &
            MATMUL( TRANSPOSE(BWork(2:2,1:DOFs)),BWork(2:2,1:DOFs)) ) * &
            (B11/A11 + B22/A22) * StressVec(1) * Weight
 
        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) - ( &
            MATMUL( TRANSPOSE(BM(2:2,1:DOFs)),BM(2:2,1:DOFs))/A2**2 + &
            MATMUL( TRANSPOSE(BWork(3:3,1:DOFs)),BWork(3:3,1:DOFs))/A1**2 + &
            MATMUL( TRANSPOSE(BWork(4:4,1:DOFs)),BWork(4:4,1:DOFs)) ) * &
            (B11/A11 + B22/A22) * StressVec(2) * Weight

        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) - (B11/A11 + B22/A22) * StressVec(3) * Weight * ( &
            (MATMUL( TRANSPOSE(BWork(3:3,1:DOFs)),BM(1:1,1:DOFs)) + &
            MATMUL( TRANSPOSE(BM(1:1,1:DOFs)),BWork(3:3,1:DOFs))) * 1.0d0/A1**2 + &
            (MATMUL( TRANSPOSE(BWork(1:1,1:DOFs)),BM(2:2,1:DOFs)) + &
            MATMUL( TRANSPOSE(BM(2:2,1:DOFs)),BWork(1:1,1:DOFs))) * 1.0d0/A2**2 + &
            MATMUL( TRANSPOSE(BWork(2:2,1:DOFs)),BWork(4:4,1:DOFs)) + &
            MATMUL( TRANSPOSE(BWork(4:4,1:DOFs)),BWork(2:2,1:DOFs)) )
      END IF

      !----------------------------------------------------------------------------------------
      ! The part of transverse shear strains which depend linearly on the thickness coordinate: 
      !----------------------------------------------------------------------------------------
      IF (.NOT. DrillingDOFs .AND. .NOT. BenchmarkProblem) THEN
        BS(3,6:DOFs:m) = dBasis(1:nd,1)
        BS(4,6:DOFs:m) = dBasis(1:nd,2)

        !--------------------------------------------------------------------------------------
        ! The nonlinear part of transverse shear strains which depend linearly on the thickness 
        ! coordinate: 
        !--------------------------------------------------------------------------------------
        StrainVec(1:2) = 0.0d0
        IF (LargeDeflection .AND. NonlinearBending) THEN
          NonlinBS(3,1:DOFs) = -BWork(5,1:DOFs) * PrevField(8) / A1**2 - &
              BWork(8,1:DOFs) * PrevField(5) / A1**2 - &
              BWork(6,1:DOFs) * PrevField(9) / A2**2 - BWork(9,1:DOFs) * PrevField(6) / A2**2 - &
              BWork(7,1:DOFs) * PrevField(10) - BWork(10,1:DOFs) * PrevField(7)

          NonlinBS(4,1:DOFs) = -BWork(5,1:DOFs) * PrevField(11) / A1**2 - &
              BWork(11,1:DOFs) * PrevField(5) / A1**2 - &
              BWork(6,1:DOFs) * PrevField(12) / A2**2 - BWork(12,1:DOFs) * PrevField(6) / A2**2 - &
              BWork(7,1:DOFs) * PrevField(13) - BWork(13,1:DOFs) * PrevField(7)

          NonlinBS(3,:) = NonlinBS(3,:) - B11/A11 * NonlinBS(1,:)
          NonlinBS(4,:) = NonlinBS(4,:) - B22/A22 * NonlinBS(2,:)

          StrainVec(1:2) = 0.5_dp * MATMUL(NonlinBS(3:4,1:DOFs), PrevSolVec(1:DOFs))
        END IF

        StrainVec(1:2) = StrainVec(1:2) + MATMUL(BS(3:4,1:DOFs), PrevSolVec(1:DOFs))

        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + Weight * &
            MATMUL(TRANSPOSE(BS(3:4,1:DOFs) + NonlinBS(3:4,1:DOFs)), &
            MATMUL(GMat(1:2,1:2), BS(3:4,1:DOFs) + NonlinBS(3:4,1:DOFs)))
        StressVec(1:2) = MATMUL(GMat(1:2,1:2), StrainVec(1:2))

        ! Residual terms for RHS:
        !
        Force(1:DOFs) = Force(1:DOFs) - MATMUL(TRANSPOSE(BS(3:4,1:DOFs) + NonlinBS(3:4,1:DOFs)), &
            StressVec(1:2) ) * Weight

        ! The remaining terms for the complete Newton iteration:
        !
        IF (LargeDeflection .AND. NonlinearBending) THEN
          Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) - StressVec(1) * Weight * ( &
              (MATMUL(TRANSPOSE(BWork(8:8,1:DOFs)), BWork(5:5,1:DOFs)) + &
              MATMUL(TRANSPOSE(BWork(5:5,1:DOFs)), BWork(8:8,1:DOFs))) * 1.0d0/A1**2 - &
              (MATMUL(TRANSPOSE(BWork(9:9,1:DOFs)), BWork(6:6,1:DOFs)) + &
              MATMUL(TRANSPOSE(BWork(6:6,1:DOFs)), BWork(9:9,1:DOFs))) * 1.0d0/A2**2 - &
              MATMUL(TRANSPOSE(BWork(10:10,1:DOFs)), BWork(7:7,1:DOFs)) + &
              MATMUL(TRANSPOSE(BWork(7:7,1:DOFs)), BWork(10:10,1:DOFs)) )

          Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) - StressVec(2) * Weight * ( &
              (MATMUL(TRANSPOSE(BWork(11:11,1:DOFs)), BWork(5:5,1:DOFs)) + &
              MATMUL(TRANSPOSE(BWork(5:5,1:DOFs)), BWork(11:11,1:DOFs))) * 1.0d0/A1**2 - &
              (MATMUL(TRANSPOSE(BWork(12:12,1:DOFs)), BWork(6:6,1:DOFs)) + &
              MATMUL(TRANSPOSE(BWork(6:6,1:DOFs)), BWork(12:12,1:DOFs))) * 1.0d0/A2**2 - &
              MATMUL(TRANSPOSE(BWork(13:13,1:DOFs)), BWork(7:7,1:DOFs)) + &
              MATMUL(TRANSPOSE(BWork(7:7,1:DOFs)), BWork(13:13,1:DOFs)) )
          
          Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + B11/A11 * StressVec(1) * Weight * ( &
              (MATMUL( TRANSPOSE(BWork(5:5,1:DOFs)),BM(1:1,1:DOFs)) + &
              MATMUL( TRANSPOSE(BM(1:1,1:DOFs)),BWork(5:5,1:DOFs)) ) * 1.0d0/A1**2 + &
              (MATMUL( TRANSPOSE(BWork(1:1,1:DOFs)),BWork(6:6,1:DOFs)) + &
              MATMUL( TRANSPOSE(BWork(6:6,1:DOFs)),BWork(1:1,1:DOFs)) ) * 1.0d0/A2**2 + &
              MATMUL( TRANSPOSE(BWork(2:2,1:DOFs)),BWork(7:7,1:DOFs)) + &
              MATMUL( TRANSPOSE(BWork(7:7,1:DOFs)),BWork(2:2,1:DOFs)) )

          Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + B22/A22 * StressVec(2) * Weight * ( &
              (MATMUL( TRANSPOSE(BWork(5:5,1:DOFs)),BWork(3:3,1:DOFs)) + &
              MATMUL( TRANSPOSE(BWork(3:3,1:DOFs)),BWork(5:5,1:DOFs)) ) * 1.0d0/A1**2 + &
              (MATMUL( TRANSPOSE(BM(2:2,1:DOFs)),BWork(6:6,1:DOFs)) + &
              MATMUL( TRANSPOSE(BWork(6:6,1:DOFs)),BM(2:2,1:DOFs)) ) * 1.0d0/A2**2 + &
              MATMUL( TRANSPOSE(BWork(4:4,1:DOFs)),BWork(7:7,1:DOFs)) + &
              MATMUL( TRANSPOSE(BWork(7:7,1:DOFs)),BWork(4:4,1:DOFs)) ) 
        END IF
      END IF


      !----------------------------------------------------------------
      ! Mass matrix without bubbles taken into account:
      !----------------------------------------------------------------     
      IF ( MassAssembly ) THEN
        DO k=1,3
          SELECT CASE(k)
          CASE(1)
            Weight = 1/a11 * h * rho0 * SqrtDetA * detJ * sq
          CASE(2)
            Weight = 1/a22 * h * rho0 * SqrtDetA * detJ * sq
          CASE(3)
            Weight = h * rho0 * SqrtDetA * detJ * sq
          END SELECT
          DO i=1,nd
            DO j=1,nd
              Mass((i-1)*m+k,(j-1)*m+k) = Mass((i-1)*m+k,(j-1)*m+k) + &
                  Basis(i) * Basis(j) * Weight
              Damp((i-1)*m+k,(j-1)*m+k) = Damp((i-1)*m+k,(j-1)*m+k) + &
                  DampCoef * Basis(i) * Basis(j) * Weight              

              IF (k > 2 .AND. DrillingDOFs) CYCLE
              Mass((i-1)*m+3+k,(j-1)*m+3+k) = Mass((i-1)*m+3+k,(j-1)*m+3+k) + &
                  h**2/12.0d0 * Basis(i) * Basis(j) * Weight
            END DO
          END DO
        END DO
      END IF

      !----------------------------------------------------------------
      ! RHS vector:
      !----------------------------------------------------------------
      IF (LargeDeflection) THEN
        !----------------------------------------------------------------
        ! Compute the normal vector n to the deformed mid-surface using
        ! the current iterate and apply the normal traction p * n, with
        ! the effect of area change being taken into account. 
        !----------------------------------------------------------------
        v1 = SUM(Basis(1:nd) * PrevSolVec(1:DOFs:m))
        v2 = SUM(Basis(1:nd) * PrevSolVec(2:DOFs:m))
        v3 = SUM(Basis(1:nd) * PrevSolVec(3:DOFs:m))
        abasis1New(1:3) = abasis1(1:3) + (SUM(dBasis(1:nd,1) * PrevSolVec(1:DOFs:m)) - &
            C111 * v1 - C112 * v2 - B11 * v3) * dual1(1:3) + &
            (SUM(dBasis(1:nd,1) * PrevSolVec(2:DOFs:m)) - C211 * v1 - C212 * v2) * dual2(1:3) + &
            (SUM(dBasis(1:nd,1) * PrevSolVec(3:DOFs:m)) + B11/A11 * v1) * abasis3(1:3)
        abasis2New(1:3) = abasis2(1:3) + (SUM(dBasis(1:nd,2) * PrevSolVec(1:DOFs:m)) - &
            C211 * v1 - C212 * v2) * dual1(1:3) + (SUM(dBasis(1:nd,2) * PrevSolVec(2:DOFs:m)) - &
            C221 * v1 - C222 * v2 - B22 * v3) * dual2(1:3) + &
            (SUM(dBasis(1:nd,2) * PrevSolVec(3:DOFs:m)) + B22/A22 * v2) * abasis3(1:3)
        NewDetA = DOT_PRODUCT(abasis1New,abasis1New) * DOT_PRODUCT(abasis2New,abasis2New) - &
            DOT_PRODUCT(abasis1New,abasis2New)**2
        abasis3New(1:3) = CrossProduct(abasis1New,abasis2New)
        Norm = SQRT(SUM(abasis3New(:)**2))
        abasis3New(1:3) = abasis3New(1:3)/Norm

        Weight = SQRT(NewDetA) * detJ * sq

        RHSForce(1:DOFs:m) = RHSForce(1:DOFs:m) + NormalTraction * DOT_PRODUCT(abasis3New,dual1) * Basis(1:nd) * Weight       
        RHSForce(2:DOFs:m) = RHSForce(2:DOFs:m) + NormalTraction * DOT_PRODUCT(abasis3New,dual2) * Basis(1:nd) * Weight       
        RHSForce(3:DOFs:m) = RHSForce(3:DOFs:m) + NormalTraction * DOT_PRODUCT(abasis3New,abasis3) * Basis(1:nd) * Weight
        Force(1:DOFs:m) = Force(1:DOFs:m) + NormalTraction * DOT_PRODUCT(abasis3New,dual1) * Basis(1:nd) * Weight       
        Force(2:DOFs:m) = Force(2:DOFs:m) + NormalTraction * DOT_PRODUCT(abasis3New,dual2) * Basis(1:nd) * Weight       
        Force(3:DOFs:m) = Force(3:DOFs:m) + NormalTraction * DOT_PRODUCT(abasis3New,abasis3) * Basis(1:nd) * Weight
        ! TO DO: Add terms related to the first-order terms in the normal coordinate
      ELSE
        Weight = SqrtDetA * detJ * sq
        DO p=1,nd
          i = m*(p-1)+3
          Force(i) = Force(i) + NormalTraction * Basis(p) * Weight
        END DO
      END IF

      Area = Area + sq * detJ * SqrtDetA

      ! For computing the mean curvature error:
      !-----------------------------------------
      !Error = Error + sq * (1.0d-1 - 0.5d0*abs(B11/a11) - 0.5d0*abs(B22/a22))**2 * detJ * SqrtDetA

    END DO QUADRATURELOOP

    ! ------------------------------------------------------------------------------
    ! Static condensation is performed before transforming to the global DOFs:
    ! ------------------------------------------------------------------------------
    IF (UseBubbles .AND. BubbleDOFs > 0) THEN
      CALL CondensateP( DOFs, BubbleDOFs, Stiff, Force )
    END IF

    !-------------------------------------------------------
    ! Transform to the global DOFs:
    !-------------------------------------------------------
    Stiff(1:DOFs,1:DOFs) = MATMUL(TRANSPOSE(Q(1:DOFs,1:DOFs)),MATMUL(Stiff(1:DOFs,1:DOFs),Q(1:DOFs,1:DOFs)))
    Force(1:DOFs) = MATMUL(TRANSPOSE(Q(1:DOFs,1:DOFs)),Force(1:DOFs))

    IF (LargeDeflection) THEN
      RHSForce(1:DOFs) = MATMUL(TRANSPOSE(Q(1:DOFs,1:DOFs)),RHSForce(1:DOFs))
    END IF

    IF ( MassAssembly ) THEN
      Mass(1:DOFs,1:DOFs) = MATMUL(TRANSPOSE(Q(1:DOFs,1:DOFs)),MATMUL(Mass(1:DOFs,1:DOFs),Q(1:DOFs,1:DOFs)))
      Damp(1:DOFs,1:DOFs) = MATMUL(TRANSPOSE(Q(1:DOFs,1:DOFs)),MATMUL(Damp(1:DOFs,1:DOFs),Q(1:DOFs,1:DOFs)))

      IF ( TransientSimulation ) THEN
        CALL Default2ndOrderTime(MASS,DAMP,STIFF,FORCE)
      ELSE IF ( HarmonicAssembly ) THEN
        CALL DefaultUpdateMass( MASS )
        CALL DefaultUpdateDamp( DAMP )
      END IF
    END IF

    CALL DefaultUpdateEquations(STIFF,FORCE)

!------------------------------------------------------------------------------
  END SUBROUTINE ShellLocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! This subroutine is used to compute an elementwise load vector arising
! from given resultant force and resultant couple vectors over the 1-D
! boundary.
!------------------------------------------------------------------------------
  SUBROUTINE ShellBoundaryMatrix(BGElement, n, nd, m, LargeDeflection, &
      MassAssembly, HarmonicAssembly, LocalSol, RHSForce, Parent, nd_parent, &
      CartesianFormulation, SkipBlending)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: BGElement  ! A boundary element of background mesh
    INTEGER, INTENT(IN) :: n                           ! The number of background element nodes
    INTEGER, INTENT(IN) :: nd                          ! The number of DOFs per component 
    INTEGER, INTENT(IN) :: m                           ! The number of DOFs per node
    LOGICAL, INTENT(IN) :: LargeDeflection             ! To activate nonlinear terms
    LOGICAL, INTENT(IN) :: MassAssembly                ! To activate mass matrix integration
    LOGICAL, INTENT(IN) :: HarmonicAssembly            ! To activate the global mass matrix updates
    REAL(KIND=dp), INTENT(IN) :: LocalSol(:,:)         ! The previous solution iterate
    REAL(KIND=dp), INTENT(OUT) :: RHSForce(:)          ! Local RHS vector corresponding to external loads
    TYPE(Element_t), POINTER, INTENT(IN) :: Parent     ! The parent of the boundary element
    INTEGER, INTENT(IN) :: nd_parent                   ! The number of parent DOFs per component   
    LOGICAL, INTENT(IN) :: CartesianFormulation        ! Defines the way how the surface basis is obtained
    LOGICAL, INTENT(IN) :: SkipBlending                ! Informs whether surface reconstruction has been done
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: GElement
    TYPE(ValueList_t), POINTER :: BC, BodyParams
    TYPE(Nodes_t) :: Nodes, ParentNodes
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: Found, AssemblyNeeded, AssembleSprings, AssembleMass, Stat
    LOGICAL :: LiveLoads, Spherical, Cylindrical, PlateBody

    INTEGER :: i, i0, j, k, t

    REAL(KIND=dp), POINTER :: TaylorParams(:)
    REAL(KIND=dp) :: Stiff(m*nd,m*nd), Mass(m*nd,m*nd), Damp(m*nd,m*nd)
    REAL(KIND=dp) :: Force(m*nd), Basis(nd), ParentBasis(nd_parent), dParentBasis(nd_parent,3)
    REAL(KIND=dp) :: PrevSolVec(m*nd)
    REAL(KIND=dp) :: NodalForce(3,n), NodalCouple(3,n)
    REAL(KIND=dp) :: NodalSprings(6,n), NodalMass(6,n)
    REAL(KIND=dp) :: ResultantForce(3), ResultantCouple(3)
    REAL(KIND=dp) :: Spring(6), MassVals(6)
    REAL(KIND=dp) :: detJ, Weight
    REAL(KIND=dp) :: ParentSol(m,nd_parent)
    REAL(KIND=dp) :: u, v, w, detF, y1, y2, K1, K2
    REAL(KIND=dp) :: Norm, CovariantBasis(3,3), NewCovariantBasis(3,3) 
    REAL(KIND=dp) :: PrevGrad(3,2)
    REAL(KIND=dp) :: PatchNodes(MaxPatchNodes,2)
    REAL(KIND=dp) :: e1(3), e2(3), e3(3)
    REAL(KIND=dp) :: o(3)
    REAL(KIND=dp) :: abasis1(3), abasis2(3), abasis3(3)
    REAL(KIND=dp) :: A11, A22, SqrtDetA
    REAL(KIND=dp) :: B11, B22   
    REAL(KIND=dp) :: C111, C112, C221, C222, C211, C212

    SAVE Nodes, ParentNodes, GElement
!------------------------------------------------------------------------------
    RHSForce = 0.0d0

    BC => GetBC()
    IF (.NOT.ASSOCIATED(BC)) RETURN

    CALL GetRealVector(BC, NodalForce(1:3,1:n), 'Resultant Force', AssemblyNeeded)
    CALL GetRealVector(BC, NodalCouple(1:3,1:n), 'Resultant Couple', Found)
    AssemblyNeeded = AssemblyNeeded .OR. Found
    CALL GetRealVector(BC, NodalSprings(1:6,1:n), 'Spring', AssembleSprings)
    AssemblyNeeded = AssemblyNeeded .OR. AssembleSprings
    IF (MassAssembly) THEN
      CALL GetRealVector(BC, NodalMass(1:6,1:n), 'Mass', AssembleMass)
      AssemblyNeeded = AssemblyNeeded .OR. AssembleMass
    ELSE
      AssembleMass = .FALSE.
    END IF
    IF (.NOT. AssemblyNeeded) RETURN

    ! The loads are treated as dead loads by default. Alternatively one can
    ! define "live" resultant force/couple loads which depend on the deformation.
    !
    LiveLoads = .NOT. GetLogical(BC, 'Dead Loads', Found)
    IF (.NOT. Found) LiveLoads = .FALSE.
    LiveLoads = LiveLoads .AND. LargeDeflection .AND. ASSOCIATED(Parent)

    CALL GetElementNodes(Nodes)

    IF (LiveLoads) THEN
      CALL GetElementNodes(ParentNodes, Parent)
      CALL GetVectorLocalSolution(ParentSol, UElement=Parent)

      IF (SkipBlending) THEN
        BodyParams => GetBodyParams(Parent)
        Spherical = GetLogical(BodyParams, 'Spherical Body', Found)
        Cylindrical = GetLogical(BodyParams, 'Cylindrical Body', Found)
      ELSE
        ! ------------------------------------------------------------------------------
        ! Retrieve the data which have been saved as elementwise properties:
        ! ------------------------------------------------------------------------------
        CALL RetrieveLocalFrame(Parent, TaylorParams, PatchNodes, e1, e2, e3, &
            o, PlateBody, Spherical, GElement)

        ! --------------------------------------------------------------------------
        ! Overwrite the coordinate arrays of the structure ParentNodes so that it 
        ! represents the domain of the principal curvature coordinates:
        ! --------------------------------------------------------------------------
        CALL SolveNodesVariables(Parent, ParentNodes, nd_parent, GElement, PatchNodes)
      END IF
    END IF

    Force = 0.0d0
    Stiff = 0.0d0
    Mass = 0.0d0
    Damp = 0.0d0

    ! ------------------------------------------------------------------------
    ! Vectorize the previous solution 
    ! ------------------------------------------------------------------------
    IF (LargeDeflection .AND. AssembleSprings) THEN
      DO k=1,m
        PrevSolVec(k:m*nd:m) = LocalSol(k,1:nd)
      END DO
    END IF

    ! Note that in the following the integration is not done accurately as
    ! the effect of the improved surface reconstruction is not taken into
    ! account. The weight would contain some small correction terms
    ! if the effect of metric tensor would be considered precisely.
    !
    IP = GaussPoints(BGElement)
    DO t=1,IP % n
      stat = ElementInfo(BGElement, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis)

      ResultantForce(1:3) = MATMUL(NodalForce(1:3,1:n), Basis(1:n))
      ResultantCouple(1:3) = MATMUL(NodalCouple(1:3,1:n), Basis(1:n))
      Weight = IP % s(t) * DetJ

      IF (LiveLoads) THEN
        !
        ! Create live loads whose orientations depend on the deformation. 
        ! First, find basis functions for the parent elements:
        !
        CALL GetParentUVW(BGElement, n, Parent, Parent % Type % NumberOfNodes, u, &
            v, w, Basis)
        stat = ElementInfo(Parent, ParentNodes, u, v, w, detF, ParentBasis, dParentBasis)

        y1 = SUM( ParentNodes % x(1:nd_parent) * ParentBasis(1:nd_parent) )
        y2 = SUM( ParentNodes % y(1:nd_parent) * ParentBasis(1:nd_parent) )

        IF (SkipBlending) THEN
          !
          ! This uses a surface parametrization over a 2D domain in special cases.
          !
          CALL SurfaceBasis(y1, y2, CovariantBasis, K1, K2, Spherical, Cylindrical)

          abasis1(:) = CovariantBasis(:,1)
          abasis2(:) = CovariantBasis(:,2)
          abasis3(:) = CovariantBasis(:,3)
        ELSE
          !
          ! Here we utilize the surface reconstruction
          !
          CALL SurfaceBasisVectors(y1, y2, TaylorParams, e1, e2, e3, o, abasis1, &
            abasis2, abasis3, A11, A22, SqrtDetA, B11, B22, C111, C112, C221, C222, &
            C211, C212, PlanarPoint=PlateBody, Umbilical=Spherical)
        END IF

        PrevGrad(1:3,1:2) = MATMUL(ParentSol(1:3,1:nd_parent), dParentBasis(1:nd_parent,1:2))

        NewCovariantBasis(1:3,1) = abasis1(1:3) + PrevGrad(1:3,1)
        NewCovariantBasis(1:3,2) = abasis2(1:3) + PrevGrad(1:3,2)
        NewCovariantBasis(1:3,3) = CrossProduct(NewCovariantBasis(1:3,1), NewCovariantBasis(1:3,2))
        DO i=1,3
          Norm = SQRT(SUM(NewCovariantBasis(1:3,i)**2))
          NewCovariantBasis(1:3,i) = NewCovariantBasis(1:3,i)/Norm
        END DO
        !
        ! The meaning of load components is now rather implicit as
        ! they follow the deformation of lines of curvature, with the first
        ! component along the direction of the smallest curvature in the undeformed 
        ! configuration. In the case of ambiguity, the element mapping defines the orientation. 
        ! This is far from user-friendly. 
        ! TO DO: Improve by implementing normal-tangential components?
        !
        ResultantForce(1:3) = MATMUL(NewCovariantBasis(1:3,1:3), ResultantForce(1:3))
        ResultantCouple(1:3) = MATMUL(NewCovariantBasis(1:3,1:3), ResultantCouple(1:3))
      END IF

      DO i=1,nd
        i0 = (i-1)*m
        Force(i0+1:i0+3) = Force(i0+1:i0+3) + Weight * ResultantForce(1:3) * Basis(i)
        Force(i0+4:i0+6) = Force(i0+4:i0+6) + Weight * ResultantCouple(1:3) * Basis(i)
      END DO

      IF (AssembleSprings) THEN
        Spring(1:6) = MATMUL(NodalSprings(1:6,1:n), Basis(1:n))
        DO k=1,3
          DO i=1,nd
            DO j=1,nd
              Stiff((i-1)*m+k,(j-1)*m+k) = Stiff((i-1)*m+k,(j-1)*m+k) + &
                  Spring(k) * Basis(i) * Basis(j) * Weight
              Stiff((i-1)*m+3+k,(j-1)*m+3+k) = Stiff((i-1)*m+3+k,(j-1)*m+3+k) + &
                  Spring(k+3) * Basis(i) * Basis(j) * Weight
            END DO
          END DO
        END DO
      END IF

      IF (MassAssembly .AND. AssembleMass) THEN
        MassVals(1:6) = MATMUL(NodalMass(1:6,1:n), Basis(1:n))
        DO k=1,3
          DO i=1,nd
            DO j=1,nd
              Mass((i-1)*m+k,(j-1)*m+k) = Mass((i-1)*m+k,(j-1)*m+k) + &
                  MassVals(k) * Basis(i) * Basis(j) * Weight
              Mass((i-1)*m+3+k,(j-1)*m+3+k) = Mass((i-1)*m+3+k,(j-1)*m+3+k) + &
                  MassVals(k+3) * Basis(i) * Basis(j) * Weight
            END DO
          END DO
        END DO
      END IF

    END DO

    RHSForce(1:m*nd) = Force(1:m*nd)

    IF (LargeDeflection .AND. AssembleSprings) THEN
      Force(1:m*nd) = Force(1:m*nd) - MATMUL(Stiff(1:m*nd,1:m*nd), PrevSolVec(1:m*nd)) 
    END IF

    IF (MassAssembly .AND. AssembleMass) THEN
      IF (HarmonicAssembly) THEN
        CALL DefaultUpdateMass(Mass)
      ELSE
        CALL Default2ndOrderTime(Mass, Damp, Stiff, Force)
      END IF
    END IF

    CALL DefaultUpdateEquations(Stiff, Force)
!------------------------------------------------------------------------------
  END SUBROUTINE ShellBoundaryMatrix
!------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------    
!> This subroutine constrains the shell solution to obey the 3-D displacement 
!> field which should be solved before the shell solver is executed. The rotation-
!> like variables (directional derivatives) of the shell model are calculated from 
!> the solution at solid nodes found in the positive and negative direction 
!> of the shell normal (director). The displacement vector of the mid-surface is
!> also constrained. The implementation is not yet perfect as when
!> the nodes where this constraint is activated are listed by using
!> the "Target Nodes" keyword, a redundant director specification should be
!> given in a Boundary Condition section of the sif file. However, the redundant 
!> definition is not needed if the constraint is activated by referring to 
!> the boundary numbering, which is expected to be a more convenient way.
! ---------------------------------------------------------------------------------    
  SUBROUTINE SetSolidCouplingBCs(Model, Solver, Displacement)
! ---------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Model_t), INTENT(IN) :: Model                 !< The current model structure
    TYPE(Solver_t), INTENT(INOUT) :: Solver            !< The shell solver
    TYPE(Variable_t), POINTER, INTENT(INOUT) :: Displacement !< The variable of 3-D elasticity 
! ---------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Matrix_t), POINTER :: ShellMatrix, A
    TYPE(ValueList_t), POINTER :: ValueList
    TYPE(Element_t), POINTER :: BGElement, Parent

    LOGICAL :: Found, SolidCoupling, GivenTargetNodes

    INTEGER, ALLOCATABLE, TARGET :: BoundaryNodes(:)
    INTEGER, ALLOCATABLE :: ActiveElementList(:)
    INTEGER, ALLOCATABLE :: NearNodes(:)
    INTEGER, POINTER :: Perm(:), NodeIndices(:)
    INTEGER, POINTER :: Cols(:), Rows(:), Diag(:)
    INTEGER :: BC, TargetCount, TargetNode, TargetInd, Row, ShellDOFs, DOFs
    INTEGER :: i, j, k, l, n, p, jz, lz, np, i0
    INTEGER :: Active, Family

    REAL(KIND=dp), ALLOCATABLE :: NearCoordinates(:,:), AllDirectors(:,:)
    REAL(KIND=dp), POINTER :: DirectorValues(:)
    REAL(KIND=dp) :: res_z, maxres_z, minres_z
    REAL(KIND=dp) :: U_mid(3), U_upper(3), U_lower(3), h_eff
    REAL(KIND=dp) :: d(3), e3(3), d_h(3), v(3), DNU(3)
! ---------------------------------------------------------------------------------
    IF (.NOT. ListGetLogicalAnyBC(Model, 'Structure Interface')) RETURN
    IF (.NOT. Displacement % DOFs == 3) THEN
      CALL Warn('SetSolidCouplingBCs', 'Structure coupling possible in 3D only')
      RETURN
    ELSE
      DOFs = 3
    END IF
    
    IF (.NOT. ASSOCIATED(Displacement % Solver)) CALL Fatal('SetSolidCouplingBCs', &
        'The solver pointer of displacement variable is not associated')

    CALL Info('SetSolidCouplingBCs', 'Creating BCs for solid-shell coupling', Level=9)

    ShellMatrix => Solver % Matrix
    ShellDOFs = Solver % Variable % DOFs
    IF (.NOT. ALLOCATED(ShellMatrix % ConstrainedDOF)) &
        ALLOCATE(ShellMatrix % ConstrainedDOF(ShellMatrix % NumberOfRows))

    Mesh => Model % Solver % Mesh

    A => Displacement % Solver % Matrix
    Diag => A % Diag
    Rows => A % Rows
    Cols => A % Cols
    Perm => Displacement % Perm

    IF (.NOT. ASSOCIATED(A % InvPerm)) THEN
      ALLOCATE(A % InvPerm(A % NumberOfRows))
      DO i = 1,SIZE(Perm)
        IF (Perm(i) > 0) THEN
          A % InvPerm(Perm(i)) = i
        END IF
      END DO
    END IF

    Active = GetNOFBoundaryElements()
    DO BC=1,Model % NumberOfBCs
      ValueList => Model % BCs(BC) % Values
      SolidCoupling = ListGetLogical(ValueList, 'Structure Interface', Found)
      IF (.NOT. SolidCoupling) CYCLE

      NodeIndices => ListGetIntegerArray(ValueList, 'Target Nodes', GivenTargetNodes)

      IF (GivenTargetNodes) THEN
        TargetCount = SIZE(NodeIndices)
        !
        ! Here the director definition is sought from the BC definition, 
        ! although the director should already be available from the specification 
        ! of a shell model. The redundant definition is not needed if the 
        ! constraint is activated by referring to the boundary numbering.
        !
        IF (.NOT.ListCheckPresent(ValueList, 'Director 1') .AND. &
            .NOT.ListCheckPresent(ValueList, 'Director 2') .AND. &
            .NOT.ListCheckPresent(ValueList, 'Director 3')) THEN
          CALL Fatal('SetSolidCouplingBCs', &
              'Director must be defined in the BC section when using Target Nodes')
        ELSE
          ALLOCATE(AllDirectors(3,TargetCount))
          AllDirectors(1,1:TargetCount) = ListGetReal(ValueList, 'Director 1', TargetCount, NodeIndices, Found)
          AllDirectors(2,1:TargetCount) = ListGetReal(ValueList, 'Director 2', TargetCount, NodeIndices, Found)
          AllDirectors(3,1:TargetCount) = ListGetReal(ValueList, 'Director 3', TargetCount, NodeIndices, Found)
        END IF        

      ELSE
        !
        ! We shall loop over the elements in order to list target nodes where the coupling BC 
        ! is activated. The value of director must be sought from the parent elements.
        !
        IF (.NOT. ALLOCATED(ActiveElementList)) ALLOCATE(ActiveElementList(Active))
        ActiveElementList = 0
        !
        ! First, figure out suitable sizes for some arrays:
        !
        n = 0
        p = 0
        DO k=1, Active
          BGElement => GetBoundaryElement(k)
          Family = GetElementFamily(BGElement)
          IF (.NOT.(ActiveBoundaryElement() .AND. Family == 2)) CYCLE
          IF (BGElement % BoundaryInfo % Constraint /= Model % BCs(BC) % Tag ) CYCLE
          p = p + BGElement % TYPE % NumberOfNodes
          n = n + 1
          ActiveElementList(n) = k
        END DO
        ! print *, 'ACTIVE BC ELEMENTS = ', N
        ! print *, 'ACTIVE BC NODES = ', P       
        ALLOCATE(BoundaryNodes(p))
        ALLOCATE(AllDirectors(3,p))
        !
        ! Then, write data to the arrays:
        !
        l = 0
        DO k=1,n
          BGElement => GetBoundaryElement(ActiveElementList(k))
          DirectorValues => NULL()
          Parent => BGElement % BoundaryInfo % Left
          IF (ASSOCIATED(Parent)) &
              DirectorValues => GetElementalDirector(Parent)
          IF (.NOT. ASSOCIATED(DirectorValues)) THEN
            Parent => BGElement % BoundaryInfo % Right
            DirectorValues => GetElementalDirector(Parent)
          END IF
          IF (.NOT. ASSOCIATED(DirectorValues)) CALL Fatal('SetSolidCouplingBCs', &
              'Director cannot be found from parent elements')

          DO i=1, BGElement % TYPE % NumberOfNodes
            DO j=1, Parent % TYPE % NumberOfNodes
              IF (BGElement % NodeIndexes(i) == Parent % NodeIndexes(j)) THEN
                i0 = 3*(j-1)
                l = l + 1
                AllDirectors(1:3,l) = DirectorValues(i0+1:i0+3)
                BoundaryNodes(l) = BGElement % NodeIndexes(i)
                EXIT
              END IF
            END DO
          END DO
        END DO
        NodeIndices => BoundaryNodes(:)
        TargetCount = l
      END IF

      DO p=1,TargetCount
        TargetNode = NodeIndices(p)
        TargetInd = Perm(NodeIndices(p))
        IF (TargetInd == 0) CYCLE
        !------------------------------------------------------------------------------
        ! Find nodes which can potentially be used to calculate the normal derivative
        ! of the 3-D solution:
        !------------------------------------------------------------------------------
        Row = TargetInd * DOFs
        n = (Rows(Row+1)-1 - Rows(Row)-Dofs+1)/DOFs + 1
        ALLOCATE(NearNodes(n), NearCoordinates(3,n))

        k = 0
        DO i = Rows(Row)+Dofs-1, Rows(Row+1)-1, Dofs
          j = Cols(i)/Dofs
          k = k + 1
          NearNodes(k) = A % InvPerm(j)
        END DO
        ! PRINT *, 'POTENTIAL NODE CONNECTIONS:'
        ! print *, 'Nodes near target=', NearNodes(1:k)       

        !
        ! The position vectors for the potential nodes:
        !
        NearCoordinates(1,1:n) = Mesh % Nodes % x(NearNodes(1:n)) - Mesh % Nodes % x(TargetNode)
        NearCoordinates(2,1:n) = Mesh % Nodes % y(NearNodes(1:n)) - Mesh % Nodes % y(TargetNode)
        NearCoordinates(3,1:n) = Mesh % Nodes % z(NearNodes(1:n)) - Mesh % Nodes % z(TargetNode)  

        d = AllDirectors(:,p)
        e3 = d/SQRT(DOT_PRODUCT(d,d))
        !------------------------------------------------------------------------------
        ! Seek for nodes which are closest to be parallel to d and have a non-negligible
        ! component with respect to d
        !------------------------------------------------------------------------------
        maxres_z = 0.0d0
        minres_z = 0.0d0
        jz = 0
        lz = 0
        DO i=1,n
          IF (NearNodes(i) == TargetNode) CYCLE

          res_z = DOT_PRODUCT(e3(:), NearCoordinates(:,i)) / &
              SQRT(DOT_PRODUCT(NearCoordinates(:,i), NearCoordinates(:,i)))
          !
          ! Skip nearly orthogonal couplings:
          !
          IF (ABS(res_z) < 2.0d-2) CYCLE

          IF (res_z > 0.0d0) THEN
            !
            ! A near node is on +d side
            !
            IF (res_z > maxres_z) THEN
              jz = NearNodes(i)
              maxres_z = res_z
            END IF
          ELSE
            !
            ! A near node is on -d side
            !
            IF (res_z < minres_z) THEN
              lz = NearNodes(i)
              minres_z = res_z
            END IF
          END IF
        END DO

        IF (jz == 0) jz = TargetNode
        IF (lz == 0) lz = TargetNode
        IF (jz == lz) CALL Fatal('SetSolidCouplingBCs', &
            'No solid nodes to span the director')

        ! PRINT *, 'HANDLING NODE = ', TargetNode
        ! PRINT *, 'UPPER NODE = ', JZ
        ! PRINT *, 'LOWER NODE = ', LZ

        ! Now, evaluate the directional derivative DNU(:) in the normal direction:
        i = Perm(lz)
        j = Perm(jz)
        k = Perm(TargetNode)
        U_lower(1:3) = Displacement % Values(i*DOFs-2:i*DOFs)
        U_upper(1:3) = Displacement % Values(j*DOFs-2:j*DOFs)
        U_mid(1:3) = Displacement % Values(k*DOFs-2:k*DOFs)
        v(1:3) = [Mesh % Nodes % x(jz) - Mesh % Nodes % x(lz), &
            Mesh % Nodes % y(jz) - Mesh % Nodes % y(lz), &
            Mesh % Nodes % z(jz) - Mesh % Nodes % z(lz)]
        h_eff = SQRT(DOT_PRODUCT(v,v))
        DNU(:) = -1.0d0/h_eff * (U_upper(:) - U_lower(:))

        d_h = v/SQRT(DOT_PRODUCT(v,v))
        IF (ABS(DOT_PRODUCT(d_h,e3)) < 0.98d0) THEN
          CALL Warn('SetSolidCouplingBCs', 'A BC omitted: Solid-model nodes does not span the director')
          CYCLE
        END IF

        !
        ! Finally, constrain the shell to follow the deformation of the solid: 
        !
        k = (Solver % Variable % Perm(TargetNode)-1) * ShellDOFs
        ShellMatrix % DValues(k+1:k+3) = U_mid(1:3)
        Solver % Matrix % ConstrainedDOF(k+1:k+3) = .TRUE.
        ShellMatrix % DValues(k+4:k+6) = DNU(1:3)
        Solver % Matrix % ConstrainedDOF(k+4:k+6) = .TRUE.
            
        DEALLOCATE(NearNodes, NearCoordinates)
      END DO
      IF (ALLOCATED(AllDirectors)) DEALLOCATE(AllDirectors)
      IF (ALLOCATED(BoundaryNodes)) DEALLOCATE(BoundaryNodes)
    END DO
    IF (ALLOCATED(ActiveElementList)) DEALLOCATE(ActiveElementList)
! ---------------------------------------------------------------------------------
  END SUBROUTINE SetSolidCouplingBCs
! ---------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! This subroutine retrieves the surface data which have been saved as elementwise 
! properties:
!------------------------------------------------------------------------------
SUBROUTINE RetrieveLocalFrame(BGElement, TaylorParams, PatchNodes, e1, e2, e3, &
    o, PlateBody, SphericalSurface, GElement)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: BGElement         ! The target element for retrieval
    REAL(KIND=dp), POINTER, INTENT(OUT) :: TaylorParams(:)    ! The coefficients of the Taylor polynomial
    REAL(KIND=dp), INTENT(OUT) :: PatchNodes(MaxPatchNodes,2) ! The nodes of principal coordinate patch
    REAL(KIND=dp), INTENT(OUT) :: e1(3), e2(3), e3(3)         ! The basis of the local frame
    REAL(KIND=dp), INTENT(OUT) :: o(3)                        ! The origin of the local frame
    LOGICAL, INTENT(OUT) :: PlateBody                         ! Indicates a planar part
    LOGICAL, INTENT(OUT) :: SphericalSurface                  ! Indicates a spherical part
    TYPE(Element_t), POINTER, INTENT(INOUT) :: GElement       ! The element structure corresponding to 
                                                              ! the surface reconstruction
!------------------------------------------------------------------------------
    INTEGER :: Family
    REAL(KIND=dp), POINTER :: PatchData(:), FrameData(:)
    REAL(KIND=dp), POINTER :: PlanarPointFlag(:), UmbilicalPointFlag(:)
!------------------------------------------------------------------------------
    TaylorParams => GetElementProperty('taylor parameters', BGElement)

    PatchData => GetElementProperty('patch nodes', BGElement) 
    PatchNodes(1:MaxPatchNodes,1) = PatchData(1:MaxPatchNodes)
    PatchNodes(1:MaxPatchNodes,2) = PatchData(MaxPatchNodes+1:2*MaxPatchNodes)
 
    FrameData => GetElementProperty('element frame', BGElement) 
    e1 = FrameData(FrameBasis1)
    e2 = FrameData(FrameBasis2)
    e3 = FrameData(FrameBasis3)
    o = FrameData(FrameOrigin)

    PlanarPointFlag => GetElementProperty('planar point', BGElement)
    UmbilicalPointFlag => GetElementProperty('umbilical point', BGElement)
    PlateBody = PlanarPointFlag(1) > 0.0d0
    SphericalSurface = UmbilicalPointFlag(1) > 0.0d0

    ! --------------------------------------------------------------------------
    ! Create the element structure corresponding to the surface reconstruction:
    ! --------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED(GElement) ) GElement => AllocateElement()
    Family = GetElementFamily(BGElement)
    SELECT CASE(Family)
    CASE(3)
      GElement % Type => GetElementType(310, .FALSE.)
    CASE(4)
      GElement % Type => GetElementType(416, .FALSE.)
    CASE DEFAULT
      CALL Fatal('RetrieveLocalFrame', 'Unsupported (geometry model) element type')
    END SELECT
!------------------------------------------------------------------------------
END SUBROUTINE RetrieveLocalFrame
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Define what strain reduction strategy is applied and set parameters that
! control the selection of variational crimes.
!------------------------------------------------------------------------------
  SUBROUTINE SetStrainReductionParameters(BGElement, ReductionMethod, PlateBody, &
      ReducedStrainDim, UseBubbles, UseShearCorrection, DOFsTransform, &
      MembraneStrains)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: BGElement ! An element of background mesh
    INTEGER, INTENT(INOUT) :: ReductionMethod         ! A desired method, the true choice may be different
    LOGICAL, INTENT(IN) :: PlateBody                  ! A dummy argument
    INTEGER, INTENT(OUT) :: ReducedStrainDim          ! The number of basis functions for strain interpolation
    LOGICAL, INTENT(INOUT) :: UseBubbles              ! To augment approximation by bubble functions
    LOGICAL, INTENT(OUT) :: UseShearCorrection        ! To activate shear correction trick
    REAL(KIND=dp), INTENT(INOUT) :: DOFsTransform(3,4)! To reduce RT_0 functions to Ker(curl)
    LOGICAL, INTENT(IN) :: MembraneStrains            ! To select the method for membrane strains
!------------------------------------------------------------------------------
    LOGICAL :: PVersion, SecondOrder
    INTEGER :: Family
!------------------------------------------------------------------------------
    Family = GetElementFamily(BGElement)
    PVersion = IsActivePElement(BGElement) 

    SecondOrder = .FALSE.
    IF (.NOT. PVersion) THEN
      SELECT CASE(Family)
      CASE(3)
        SecondOrder = BGElement % Type % NumberOfNodes == 6
      CASE(4)
        SecondOrder = BGElement % Type % NumberOfNodes == 9
      END SELECT
    END IF

    IF (.NOT. MembraneStrains .AND. (ReductionMethod == CurlKernelWithEdgeDOFs)) THEN
      CALL Warn('SetStrainReductionParameters', &
          'Operator 4 is not yet possible for transverse shear strains')
      ReductionMethod = CurlKernel
    END IF

    IF (PVersion .OR. SecondOrder) THEN
      ! If a higher-order approximation has been requested, the standard weak formulation is used:
      ReductionMethod = NoStrainReduction
    ELSE
      IF (ReductionMethod == AutomatedChoice) THEN
        ! ------------------------------------------------------------------------------
        ! This option saves the user from deciding a strain reduction method
        ! ------------------------------------------------------------------------------
        SELECT CASE(Family)
        CASE(3)
          IF (MembraneStrains) THEN
            ReductionMethod = DoubleReduction
          ELSE
            ReductionMethod = DoubleReduction
            UseBubbles = .FALSE.
          END IF
        CASE(4)
          IF (MembraneStrains) THEN
            ReductionMethod = CurlKernelWithEdgeDOFs
          ELSE
            ReductionMethod = MITC
          END IF
        END SELECT
      ELSE
        ! ------------------------------------------------------------------------------
        ! In the case of triangles asking for the kernel version switches to the
        ! compound strain reduction strategy:
        ! ------------------------------------------------------------------------------
        IF (Family == 3 .AND. ReductionMethod == CurlKernel) &
            ReductionMethod = DoubleReduction
        ! Currently bubbles can be activated only for certain strategies:
        IF (.NOT.(ReductionMethod == NoStrainReduction .OR. &
            ReductionMethod == DoubleReduction)) THEN
          UseBubbles = .FALSE.
        END IF
      END IF
    END IF

    ! ------------------------------------------------------------------------------
    ! Set the dimension of the range of strain reduction operator and decide
    ! whether a numerical shear correction trick is applied:
    ! ------------------------------------------------------------------------------
    SELECT CASE(ReductionMethod)
    CASE(NoStrainReduction)
      UseShearCorrection = .FALSE.
    CASE(CurlKernel)
      IF (Family == 3) THEN
        ! This choice should not really be used for triangles (without applying MITC 
        ! interpolation first, cf. DoubleReduction case); with the current version
        ! one should never end up here
        ReducedStrainDim = 2
        UseShearCorrection = .TRUE.
      ELSE
        ReducedStrainDim = 3
        UseShearCorrection = .TRUE.
      END IF
    CASE(MITC)
      IF (Family == 3) THEN
        ReducedStrainDim = 3
      ELSE
        ReducedStrainDim = 4
      END IF
      UseShearCorrection = .TRUE.
    CASE(DoubleReduction)
      !Use a combination of MITC and the kernel reduction:
      IF (Family==3) THEN
        ReducedStrainDim = 2
        IF (UseBubbles) THEN
          UseShearCorrection = .FALSE.
        ELSE
          UseShearCorrection = .TRUE.
        END IF
        ! Coefficients to transform from MITC3 DOFs to the Ker(curl) DOFs: 
        DOFsTransform(1,1) = 1.0d0/3.0d0
        DOFsTransform(1,2) = -1.0d0/6.0d0
        DOFsTransform(1,3) = DOFsTransform(1,2)
        DOFsTransform(2,1) = 0.0d0
        DOFsTransform(2,2) = 1.0d0/(2.0d0*sqrt(3.0d0))
        DOFsTransform(2,3) = -1.0d0/(2.0d0*sqrt(3.0d0))
      ELSE
        ReducedStrainDim = 3
        UseShearCorrection = .TRUE.
        ! Coefficients to transform from MITC4 DOFs to the Ker(curl) DOFs:
        ! This effectively removes components that build on local vectors (-v,u).
        DOFsTransform(1,1) = 1.0d0/4.0d0
        DOFsTransform(1,2) = 0.0d0
        DOFsTransform(1,3) = -1.0d0/4.0d0
        DOFsTransform(1,4) = 0.0d0
        DOFsTransform(2,1) = 0.0d0
        DOFsTransform(2,2) = 1.0d0/4.0d0
        DOFsTransform(2,3) = 0.0d0
        DOFsTransform(2,4) = -1.0d0/4.0d0
        DOFsTransform(3,1) = -1.0d0/8.0d0
        DOFsTransform(3,2) = 1.0d0/8.0d0
        DOFsTransform(3,3) = -1.0d0/8.0d0
        DOFsTransform(3,4) = 1.0d0/8.0d0
      END IF
    CASE(CurlKernelWithEdgeDOFs)
      IF (Family==3) THEN
        CALL Fatal('SetStrainReductionParameters', 'Strain Reduction Operator=4 is not defined for trias')
      ELSE
        ReducedStrainDim = 3
        UseShearCorrection = .TRUE.        
      END IF
    END SELECT
!------------------------------------------------------------------------------
  END SUBROUTINE SetStrainReductionParameters
!------------------------------------------------------------------------------

! --------------------------------------------------------------------------
! Allocate a Lagrange interpolation element structure corresponding to the
! "Element" keyword. A corresponding nodes data structure is also created.
! In addition, the element structure corresponding to the surface reconstruction
! is created.
!------------------------------------------------------------------------------
  SUBROUTINE CreateLagrangeElementStructures(BGElement, nd, Element, Nodes, &
      PNodes, GElement)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: BGElement  ! An element of background mesh
    INTEGER, INTENT(IN) :: nd                          ! The number of DOFs (per component)
    TYPE(Element_t), POINTER, INTENT(OUT) :: Element   ! A Lagrange element data structure
    TYPE(Nodes_t), INTENT(OUT) :: Nodes                ! A nodes data structure for the Lagrange element
    TYPE(Nodes_t), INTENT(OUT) :: PNodes               ! A nodes data structure for p-version
    TYPE(Element_t), POINTER, INTENT(OUT) :: GElement  ! The element structure for surface reconstruction
!------------------------------------------------------------------------------
    LOGICAL :: PVersion
    INTEGER :: Family
!------------------------------------------------------------------------------
    Family = GetElementFamily(BGElement)
    PVersion = IsActivePElement(BGElement)

    ! --------------------------------------------------------------------------
    ! Create the element structure corresponding to the surface reconstruction:
    ! --------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED(GElement) ) GElement => AllocateElement()
    SELECT CASE(Family)
    CASE(3)
      GElement % Type => GetElementType(310, .FALSE.)
    CASE(4)
      GElement % Type => GetElementType(416, .FALSE.)
    CASE DEFAULT
      CALL Fatal('CreateLagrangeElementStructures', 'Unsupported (geometry model) element type')
    END SELECT

    IF ( PVersion ) THEN
      IF ( .NOT. ASSOCIATED(Element) ) Element => AllocateElement()

      SELECT CASE(Family)
      CASE(3)
        SELECT CASE(nd)
        CASE(3)
          Element % Type => GetElementType(303, .FALSE.)
        CASE(6)
          Element % Type => GetElementType(306, .FALSE.)
        CASE(10)
          Element % Type => GetElementType(310, .FALSE.)
        CASE DEFAULT
          CALL Fatal('CreateLagrangeElementStructures', 'Unsupported triangular p-element')
        END SELECT
        ! Ensure that the reference element for the Lagrange interpolation is used:
        Element % Type % NodeU(1:3) = (/ 0.0d0, 1.0d0, 0.0d0 /)
        Element % Type % NodeV(1:3) = (/ 0.0d0, 0.0d0, 1.0d0 /)
      CASE(4)
        SELECT CASE(nd)
        CASE(4)
          Element % Type => GetElementType(404, .FALSE.)
        CASE(8)
          Element % Type => GetElementType(408, .FALSE.)
        CASE(9)
          Element % Type => GetElementType(409, .FALSE.)
        CASE(12)
          Element % Type => GetElementType(412, .FALSE.)
        CASE DEFAULT
          CALL Fatal('CreateLagrangeElementStructures', 'Unsupported quadrilateral p-element type')
        END SELECT
      CASE DEFAULT
        CALL Fatal('CreateLagrangeElementStructures', 'Expecting an element of surface type')
      END SELECT
    ELSE
      Element => BGElement
    END IF

    ! --------------------------------------------------------------------------
    ! Create a node variable suitable for defining the isoparametric element 
    ! map, i.e. use as many nodes as DOFs in the spatial discretization. 
    ! --------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED( Nodes % x ) ) THEN
      ALLOCATE( Nodes % x(nd), Nodes % y(nd), Nodes % z(nd) ) 
      Nodes % NumberOfNodes = nd
    ELSE
      IF (nd > SIZE(Nodes % x)) THEN
        DEALLOCATE(Nodes % x, Nodes % y, Nodes % z)
        ALLOCATE( Nodes % x(nd), Nodes % y(nd), Nodes % z(nd) ) 
      END IF
      Nodes % NumberOfNodes = nd          
    END IF

    IF ( PVersion ) THEN    
      IF ( .NOT. ASSOCIATED( PNodes % x ) ) THEN
        ALLOCATE( PNodes % x(nd), PNodes % y(nd), PNodes % z(nd) )
        PNodes % NumberOfNodes = nd
      ELSE
        IF (nd > SIZE(PNodes % x)) THEN
          DEALLOCATE(PNodes % x, PNodes % y, PNodes % z)
          ALLOCATE( PNodes % x(nd), PNodes % y(nd), PNodes % z(nd) ) 
        END IF
        PNodes % NumberOfNodes = nd          
      END IF
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE CreateLagrangeElementStructures
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! This subroutine updates the coordinate arrays of Nodes argument such that they
! correspond to the isoparametric approximation of the coordinate patch
! obtained via surface reconstruction. If p-version is used for approximating 
! the shell equations, this creates also the isoparametric geometry representation 
! in terms of p-basis.
!------------------------------------------------------------------------------
  SUBROUTINE WriteElementNodesVariables(BGElement, GElement, Element, Nodes, &
      PNodes, PatchNodes)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: BGElement ! An element of background mesh
    TYPE(Element_t), POINTER, INTENT(IN) :: GElement  ! A Lagrange element for surface reconstruction
    TYPE(Element_t), POINTER, INTENT(IN) :: Element   ! The element type for which nodes are written
    TYPE(Nodes_t), INTENT(INOUT) :: Nodes             ! A nodes data structure for Element     
    TYPE(Nodes_t), INTENT(INOUT) :: PNodes            ! A nodes data structure for p-version
    REAL(KIND=dp), INTENT(IN) :: PatchNodes(:,:)      ! The nodes data of coordinate domain
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: Stat, PVersion
    INTEGER :: j, k, t, n, NodesCount, GElementNodes, Family
    REAL(KIND=dp) :: u, v, yk, up, vp, sp
    REAL(KIND=dp) :: GBasis(GElement % Type % NumberOfNodes)
    REAL(KIND=dp) :: Stiff(Element % Type % NumberOfNodes, Element % Type % NumberOfNodes)
    REAL(KIND=dp) :: Force(Element % Type % NumberOfNodes), DetJ
    REAL(KIND=dp) :: PBasis(Element % Type % NumberOfNodes)
    REAL(KIND=dp) :: Basis(Element % Type % NumberOfNodes)
!------------------------------------------------------------------------------
    PVersion = IsActivePElement(BGElement)

    n = BGElement % Type % NumberOfNodes
    NodesCount = Element % Type % NumberOfNodes
    GElementNodes = GElement % Type % NumberOfNodes

    ! --------------------------------------------------------------------------
    ! Write the nodes variable for isoparametric approximation of coordinate patch
    ! --------------------------------------------------------------------------
    DO j=1,NodesCount
      u = Element % Type % NodeU(j)
      v = Element % Type % NodeV(j)
      CALL NodalBasisFunctions2D(GBasis, GElement, u, v)

      Nodes % x(j) = SUM(GBasis(1:GElementNodes) * PatchNodes(1:GElementNodes,1))
      Nodes % y(j) = SUM(GBasis(1:GElementNodes) * PatchNodes(1:GElementNodes,2))
      Nodes % z(j) = 0.0d0
    END DO

    ! ---------------------------------------------------------------------------
    ! Express the isoparametric representation in terms of the p-basis:
    ! ---------------------------------------------------------------------------
    IF ( PVersion ) THEN
      Family = GetElementFamily(BGElement)
      ! -----------------------------------------------
      ! The corner values can be written immediately:
      ! -----------------------------------------------
      PNodes % x(1:NodesCount) = 0.0d0
      PNodes % y(1:NodesCount) = 0.0d0
      PNodes % z(1:NodesCount) = 0.0d0
      PNodes % x(1:n) = Nodes % x(1:n)
      PNodes % y(1:n) = Nodes % y(1:n)

      ! ------------------------------------------------------
      ! Solve a local matrix equation in the case of higher p:
      ! ------------------------------------------------------
      IF (NodesCount > n) THEN
        IP = GaussPoints( BGElement )
        DO k=1,2
          STIFF = 0.0d0
          FORCE = 0.0d0

          DO t=1,IP % n
            up = IP % U(t)
            vp = IP % V(t)
            ! ------------------------------------------------------------------
            ! First, evaluate the basis functions for p-version:
            ! ------------------------------------------------------------------
            stat = ElementInfo( BGElement, PNodes, up, vp, 0.0d0, detJ, PBasis )
            ! ------------------------------------------------------------------
            ! Map the p-point to the point on the standard reference element:
            ! ------------------------------------------------------------------
            IF (Family==3) THEN
              u = 0.5d0 * (1.0d0 + up - 1.0d0/sqrt(3.0d0) * vp)
              v = 1.0d0/sqrt(3.0d0) * vp
            ELSE
              u = up
              v = vp
            END IF

            stat = ElementInfo( Element, Nodes, u, v, 0.0d0, detJ, Basis )         

            DO i=1,NodesCount
              DO j=1,NodesCount
                STIFF(i,j) = STIFF(i,j) + IP % s(t) * PBasis(i) * PBasis(j)
              END DO
            END DO

            SELECT CASE(k)
            CASE(1)
              yk = SUM( Nodes % x(1:NodesCount) * Basis(1:NodesCount) )
            CASE(2)
              yk = SUM( Nodes % y(1:NodesCount) * Basis(1:NodesCount) )           
            END SELECT
            FORCE(1:NodesCount) = FORCE(1:NodesCount) + IP % s(t) * yk * PBasis(1:NodesCount)
          END DO

          CALL LUSolve(NodesCount, Stiff, Force)
          SELECT CASE(k)
          CASE(1)
            PNodes % x(1:NodesCount) = Force(1:NodesCount)
          CASE(2)
            PNodes % y(1:NodesCount) = Force(1:NodesCount) 
          END SELECT
        END DO
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE WriteElementNodesVariables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! The matrix representation of the elasticity tensor with respect an orthogonal
! basis. The case A1 = A2 = 1 corresponds to an orthonormal basis.     
!------------------------------------------------------------------------------
  SUBROUTINE ElasticityMatrix(CMat, GMat, A1, A2, E, nu, DrillingDOFs, StabPar)
!------------------------------------------------------------------------------    
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(OUT) :: CMat(4,4), GMat(2,2)
    REAL(KIND=dp), INTENT(IN) :: A1, A2, E, nu  
    LOGICAL, OPTIONAL, INTENT(IN) :: DrillingDOFs
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: StabPar
!------------------------------------------------------------------------------
    LOGICAL :: WithDrillingDOFs
    REAL(KIND=dp) :: StabConst
!------------------------------------------------------------------------------
    IF (PRESENT(DrillingDOFs)) THEN
      WithDrillingDOFs = DrillingDOFs
      IF (WithDrillingDOFs) THEN
        IF (PRESENT(StabPar)) THEN
          StabConst = StabPar
        ELSE
          StabConst = 1.0d0
        END IF
      END IF
    ELSE
      WithDrillingDOFs = .FALSE.
    END IF

    CMat = 0.0d0
    GMat = 0.0d0

    CMat(1,1) = 1.0d0    
    CMat(1,2) = nu
    CMat(2,1) = nu
    CMat(2,2) = 1.0d0
    CMat(3,3) = (1.0d0-nu)/2.0d0
    CMat = CMat * E / (1.0d0-nu**2)

    CMat(1,1) = CMat(1,1)/A1**4
    CMat(1,2) = CMat(1,2)/(A1**2 * A2**2)
    CMat(2,1) = CMat(2,1)/(A1**2 * A2**2)
    CMat(2,2) = CMat(2,2)/A2**4   
    CMat(3,3) = CMat(3,3)/(A1**2 * A2**2)

    IF (WithDrillingDOFs) THEN
      CMat(4,4) = StabConst * 2.0d0*E/(1.0d0 + nu)
    ELSE
      ! The row corresponding to the normal stress: A deviation from the state of
      ! vanishing normal stress produces deformation energy as described by
      ! the 3-D Hooke's law.
      CMat(4,4) = (1.0d0-nu) * E /( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    END IF

    GMat(1,1) = E/(2.0d0*(1.0d0 + nu)*A1**2)
    GMat(2,2) = E/(2.0d0*(1.0d0 + nu)*A2**2)
!------------------------------------------------------------------------------
  END SUBROUTINE ElasticityMatrix
!------------------------------------------------------------------------------    

!------------------------------------------------------------------------------
! Return basis functions which give a basis for the range X(K) of a (strain) 
! reduction operator R_K: L2(K) -> X(K). By construction these basis functions 
! transform in the same way as curl-conforming FE functions, i.e. 
! a curl-conforming version of the Piola transform is applied. Hence we have 
!   X(K) = { B | B(x) = F^{-T}(f^{-1}(x)) b(f^{-1}(x)) }
! with b(p) giving the basis function on the reference element k,
! f mapping k to the physical element K = f(k) and F = Grad f. Note that 
! the reference element is chosen as in the p-approximation so that the reference 
! element edges have the same length. The functionality of this routine could 
! also be a part of the function EdgeElementInfo, but this separate implementation 
! is made to serve the special purpose of strain reduction. 
! NOTE: Only the lowest-order case is supported currently. 
!------------------------------------------------------------------------------
  FUNCTION ReductionOperatorInfo(Element, Nodes, u, v, StrainBasis, ReductionMethod, &
      ApplyPiolaTransform, F, G, detF, Basis, dBasis, DOFWeigths, Bubbles, EdgeDirection) &
      RESULT(stat)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), INTENT(IN), TARGET :: Element         !< Element structure
    TYPE(Nodes_t), INTENT(IN) :: Nodes                     !< Data corresponding to the classic element nodes
    REAL(KIND=dp), INTENT(IN) :: u                         !< 1st reference element coordinate
    REAL(KIND=dp), INTENT(IN) :: v                         !< 2nd reference element coordinate
    REAL(KIND=dp), INTENT(OUT) :: StrainBasis(:,:)         !< The basis functions b spanning the reference element space
    INTEGER, INTENT(IN) :: ReductionMethod                 !< The method chosen (integer parameters are introduced elsewhere)
    LOGICAL, INTENT(IN), OPTIONAL :: ApplyPiolaTransform   !< If .TRUE., perform the Piola transform so that, instead of b(p),
                                                           !< return B(f(p)) with B(x) the basis functions on the physical
                                                           !< element K and f the element map f:k->K
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: F(3,3)         !< The gradient F=Grad f
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: G(3,3)         !< The transpose of the inverse of the gradient F
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: detF           !< The determinant of the gradient matrix F
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: Basis(:)       !< H1-conforming basis functions corresponding to (u,v)
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: dBasis(:,:)    !< The first derivatives of the H1-conforming basis functions. If the Piola
                                                           !< transformation is performed within this subroutine, the differentiation
                                                           !< is done with respect to physical coordinates x
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: DOFWeigths(:,:)!< Auxiliary div-conforming functions needed in the evaluation of DOFs
    LOGICAL, INTENT(IN), OPTIONAL :: Bubbles               !< Indicate whether a bubble function is requested
    INTEGER, INTENT(IN), OPTIONAL :: EdgeDirection         !< Preferred direction for edge DOFs when 
                                                           !< ReductionMethod=CurlKernelWithEdgeDOFs
    LOGICAL :: Stat                                        !< Currently a dummy return value
!---------------------------------------------------------------------------------
    LOGICAL :: PerformPiolaTransform, CreateBubbles

    INTEGER :: dim, e, i, j, k, n, ntot, q, DOFs, Family

    REAL(KIND=dp) :: LF(3,3), LG(3,3), detLF, B(3)
    REAL(KIND=dp) :: LBasis(Element % TYPE % NumberOfNodes+1)
    REAL(KIND=dp) :: dLBasis(Element % TYPE % NumberOfNodes+1,3)
!------------------------------------------------------------------------------       
    StrainBasis = 0.0d0
    IF ( PRESENT(DOFWeigths) ) DOFWeigths = 0.0d0

    PerformPiolaTransform = .FALSE.
    IF ( PRESENT(ApplyPiolaTransform) ) PerformPiolaTransform = ApplyPiolaTransform

    CreateBubbles = .FALSE.
    IF ( PRESENT(Bubbles) ) CreateBubbles = Bubbles

    n = Element % TYPE % NumberOfNodes
    ntot = n
    dim = Element % TYPE % DIMENSION
    IF (dim /= 2) CALL Fatal('ReductionOperatorInfo','Reduction operators defined only for 2D elements')

    Family = GetElementFamily(Element)
    !-----------------------------------------------------------------------
    ! The lowest-order (nodal) basis functions on the reference element and
    ! their derivatives with respect to the local coordinates. These define 
    ! the mapping of the reference element to a physical element.
    !-----------------------------------------------------------------------
    LBasis = 0.0d0
    dLBasis = 0.0d0      
    SELECT CASE(Family)
    CASE(3)
      DO q=1,3
        LBasis(q) = TriangleNodalPBasis(q, u, v)
        dLBasis(q,1:2) = dTriangleNodalPBasis(q, u, v) 
      END DO
      IF (CreateBubbles) THEN
        LBasis(n+1) = TriangleBubblePBasis(0,0,u,v)
        dLBasis(n+1,1:2) = dTriangleBubblePBasis(0,0,u,v)
        ntot = n+1
      END IF
    CASE(4)
      DO q=1,4
        LBasis(q) = QuadNodalPBasis(q, u, v)
        dLBasis(q,1:2) = dQuadNodalPBasis(q, u, v) 
      END DO
      IF (CreateBubbles) THEN
        LBasis(n+1) = QuadBubblePBasis(2,2,u,v)
        dLBasis(n+1,1:2) = dQuadBubblePBasis(2,2,u,v)
        ntot = n+1
      END IF      
    CASE DEFAULT
      CALL Warn('ReductionOperatorInfo','Unsupported element type')
      RETURN
    END SELECT

    !-----------------------------------------------------------------------
    ! Get data for performing the Piola transformation...
    !-----------------------------------------------------------------------
    stat = PiolaTransformationData(n, Element, Nodes, LF, detLF, dLBasis) 
    !------------------------------------------------------------------------
    LG(1,1) = 1.0d0/detLF * LF(2,2)
    LG(1,2) = -1.0d0/detLF * LF(1,2)
    LG(2,1) = -1.0d0/detLF * LF(2,1)
    LG(2,2) = 1.0d0/detLF * LF(1,1)     
    LG(1:dim,1:dim) = TRANSPOSE( LG(1:dim,1:dim) )

    IF (ReductionMethod /= NoStrainReduction) THEN
      SELECT CASE(Family)
      CASE(3)
        SELECT CASE(ReductionMethod)
        CASE(CurlKernel,DoubleReduction)
          !---------------------------------------------------------------------
          ! The basis functions for RT_0(k,0), with DOFs defined as integrals
          ! of the type d_i = (u,v_i)_k. Here the given function v_i transforms
          ! according to the standard Piola transformation (the div-conforming 
          ! version).
          !---------------------------------------------------------------------
          DOFs = 2
          StrainBasis(1,1) = 1.0d0
          StrainBasis(1,2) = 0.0d0
          StrainBasis(2,1) = 0.0d0
          StrainBasis(2,2) = 1.0d0

          IF ( PRESENT(DOFWeigths) ) THEN
            DOFWeigths(1,1:2) = StrainBasis(1,1:2)/sqrt(3.0d0)
            DOFWeigths(2,1:2) = StrainBasis(2,1:2)/sqrt(3.0d0)
          END IF

        CASE(MITC)
          !---------------------------------------------------------------------
          ! The basis functions for RT_0(k) with DOFs attached to the edges
          !---------------------------------------------------------------------
          DOFs = 3
          StrainBasis(1,1) = (3.0d0 - Sqrt(3.0d0)*v)/6.0d0
          StrainBasis(1,2) = u/(2.0d0*Sqrt(3.0d0))
          StrainBasis(2,1) = -v/(2.0d0*Sqrt(3.0d0))
          StrainBasis(2,2) = (1 + u)/(2.0d0*Sqrt(3.0d0))
          StrainBasis(3,1) = -v/(2.0d0*Sqrt(3.0d0))
          StrainBasis(3,2) = (-1 + u)/(2.0d0*Sqrt(3.0d0))

        CASE DEFAULT
          CALL Fatal('ReductionOperatorInfo','Unknown strain reduction operator')
        END SELECT

      CASE(4)
        SELECT CASE(ReductionMethod)
        CASE(CurlKernel,DoubleReduction)
          !---------------------------------------------------------------------
          ! The basis functions for ABF_0(k,0), with DOFs defined as integrals
          ! of the type d_i = (u,v_i)_k. Here the given function v_i transforms
          ! according to the standard Piola transformation (the div-conforming 
          ! version).
          !---------------------------------------------------------------------
          DOFs = 3
          StrainBasis(1,1) = 1.0d0
          StrainBasis(1,2) = 0.0d0
          StrainBasis(2,1) = 0.0d0
          StrainBasis(2,2) = 1.0d0
          StrainBasis(3,1) = v
          StrainBasis(3,2) = u

          IF ( PRESENT(DOFWeigths) ) THEN
            DOFWeigths(1,1:2) = StrainBasis(1,1:2)/4.0d0
            DOFWeigths(2,1:2) = StrainBasis(2,1:2)/4.0d0
            DOFWeigths(3,1:2) = StrainBasis(3,1:2)*3.0d0/8.0d0
          END IF

        CASE(MITC)
          !---------------------------------------------------------------------
          ! The basis functions for RT_0(k) with DOFs attached to the edges
          !---------------------------------------------------------------------
          DOFs = 4
          StrainBasis(1,1) = (1.0d0-v)/4.0d0
          StrainBasis(1,2) = 0.0d0
          StrainBasis(2,1) = 0.0d0
          StrainBasis(2,2) = (1.0d0+u)/4.0d0
          StrainBasis(3,1) = -(1.0d0+v)/4.0d0
          StrainBasis(3,2) = 0.0d0
          StrainBasis(4,1) = 0.0d0
          StrainBasis(4,2) = (-1.0d0+u)/4.0d0

        CASE(CurlKernelWithEdgeDOFs)
          !---------------------------------------------------------------------
          ! The basis functions for ABF_0(k,0), but DOFs are now defined as
          ! edge/line integrals d_i = (u,t)_E where t is an edge tangent.
          ! Two versions are available which differ on what edges are used.
          !---------------------------------------------------------------------
          DOFs = 3
          E = 1
          IF (PRESENT(EdgeDirection)) THEN
            E = EdgeDirection
            IF (.NOT.(E==1 .OR. E==2)) &
                CALL Fatal('ReductionOperatorInfo','A wrong direction parameter given')
          END IF

          IF (E == 1) THEN
            ! Employ the edges 12 and 34 and the v-axis
            StrainBasis(1,1) = (1.0d0-v)/4.0d0
            StrainBasis(1,2) = -1.0d0/4.0d0*u
            StrainBasis(2,1) = -(1.0d0+v)/4.0d0
            StrainBasis(2,2) = -1.0d0/4.0d0*u
            StrainBasis(3,1) = 0.0d0
            StrainBasis(3,2) = 0.5d0
          ELSE
            ! Employ the edges 41 and 23 and the u-axis
            StrainBasis(1,1) = 1.0d0/4.0d0*v
            StrainBasis(1,2) = -(1.0d0-u)/4.0d0
            StrainBasis(2,1) = 1.0d0/4.0d0*v
            StrainBasis(2,2) = (1.0d0+u)/4.0d0
            StrainBasis(3,1) = 0.5d0
            StrainBasis(3,2) = 0.0d0
          END IF

        CASE DEFAULT
          CALL Fatal('ReductionOperatorInfo','Unknown strain reduction operator')
        END SELECT
      CASE DEFAULT
        CALL Warn('ReductionOperatorInfo','Unsupported element type')
        RETURN
      END SELECT
    END IF

    IF (PerformPiolaTransform) THEN
      ! ----------------------------------------------------------------------
      ! Get global first derivatives of the nodal basis functions if wanted:
      ! ----------------------------------------------------------------------
      IF ( PRESENT(dBasis) ) THEN
        dBasis = 0.0d0
        DO i=1,ntot
          DO j=1,dim
            DO k=1,dim
              dBasis(i,j) = dBasis(i,j) + dLBasis(i,k)*LG(j,k)
            END DO
          END DO
        END DO
      END IF

      IF (ReductionMethod /= NoStrainReduction) THEN
        DO j=1,DOFs
          DO k=1,dim
            B(k) = SUM( LG(k,1:dim) * StrainBasis(j,1:dim) )
          END DO
          StrainBasis(j,1:dim) = B(1:dim)
        END DO

        ! ----------------------------------------------------------------------
        ! Apply the standard Piola transform for the functions needed in the 
        ! expressions for the DOFs (the scaling with 1/detF however omitted) 
        ! ----------------------------------------------------------------------
        IF ( PRESENT(DOFWeigths) ) THEN
          DO j=1,DOFs
            DO k=1,dim
              B(k) = SUM( LF(k,1:dim) * DOFWeigths(j,1:dim) )
            END DO
            DOFWeigths(j,1:dim) = B(1:dim)
          END DO
        END IF
      END IF
    ELSE
      IF ( PRESENT(dBasis) ) dBasis(1:ntot,1:dim) = dLBasis(1:ntot,1:dim)
    END IF

    IF ( PRESENT(Basis) ) Basis(1:ntot) = LBasis(1:ntot)

    ! Make the returned value DetF to act as a metric term for integration
    ! over the volume of the element:
    IF ( PRESENT(DetF) ) DetF = ABS(DetLF)

    IF(PRESENT(F)) F = LF
    IF(PRESENT(G)) G = LG
!------------------------------------------------------------------------------
  END FUNCTION ReductionOperatorInfo
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Compute a matrix which can be used to evaluate DOFs for the interpolating
!  function R_K(u) in the strain reduction space X(K) when applied to a H1-
!  conforming FE function u having two components. If d_k denotes the linear 
!  functional defining the kth DOF, the kth row of the returned matrix A has 
!  entries 
!
!     [d_k(N1*e1) d_k(N1*e2) ... d_k(Nn*e1) d_k(Nn*e2)]
!
!  where N1,...,Nn are the Lagrange basis functions, e1=(1,0) and e2=(0,1). 
!  Thus, if the DOFs of u are contained accordingly in a vector U, the DOFs of
!  the interpolating function can be evaluated by computing the product A*U.
!  Optionally the interpolating function can be computed for a field Cu
!  where C is a 2X2 matrix field. The model parameters at the element nodes
!  are then given in the ModelPar array. This subroutine depends intimately
!  on the function ReductionOperatorInfo which gives the definition of DOFs
!  integrated here. 
!------------------------------------------------------------------------------
  SUBROUTINE ReductionOperatorDofs(Element, Nodes, A, nd, n, ReductionMethod, &
      ModelPars, GradientField, EdgeDirection)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), INTENT(IN), TARGET :: Element          !< Element structure
    TYPE(Nodes_t), INTENT(IN) :: Nodes                      !< Nodes structure
    REAL(KIND=dp), INTENT(INOUT) :: A(:,:)                  !< Coefficients for expressing the DOFs 
    INTEGER, INTENT(IN) :: nd                               !< The dimension of the strain reduction space X(K)
    INTEGER, INTENT(IN) :: n                                !< The number of the H1-conforming basis functions
    INTEGER, INTENT(IN) :: ReductionMethod                  !< The method chosen
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: ModelPars(2,2,n) !< To include the effect of additional model parameters
    LOGICAL, OPTIONAL :: GradientField                      !< To return [d_k(grad(N1).e1) d_k(grad(N1).e2) ... ]
    INTEGER, INTENT(IN), OPTIONAL :: EdgeDirection          !< Preferred direction for edge DOFs when 
                                                            !< ReductionMethod=CurlKernelWithEdgeDOFs
!---------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: stat, UseParameters, PRefElement, GradientOperand

    INTEGER :: Family, e, i, j, k, t

    REAL(KIND=dp) :: StrainBasis(nd,3)        ! The basis functions for the strain reduction space X(K)
    REAL(KIND=dp) :: DOFWeigths(nd,3)         ! The auxiliary functions to evaluate the interpolant in X(K)
    REAL(KIND=dp) :: Basis(n)                 ! H1-conforming basis functions
    REAL(KIND=dp) :: DBasis(n,1:3) 
    REAL(KIND=dp) :: F(3,3)                   ! The gradient F=Grad f of the element mapping
    REAL(KIND=dp) :: detJ
    REAL(KIND=dp) :: u(2), ParMat(2,2), uk, vk, sk
    REAL(KIND=dp) :: Tau(4,3)
!---------------------------------------------------------------------------------
    IF (PRESENT(GradientField)) THEN
      GradientOperand = GradientField
    ELSE
      GradientOperand = .FALSE.
    END IF
    Family = GetElementFamily(Element)

    IF (GradientOperand .AND. .NOT.(Family==4 .AND. (ReductionMethod == CurlKernel .OR. &
        ReductionMethod == CurlKernelWithEdgeDOFs))) &
        CALL Fatal('ReductionOperatorDofs', 'GradientOperand is not supported the chosen reduction')

    ! Clear upper left corner of A
    A(1:nd,1:2*n) = 0.0d0

    UseParameters = PRESENT(ModelPars) 
    IF (.NOT. UseParameters .OR. GradientOperand) THEN
      ParMat(1,1) = 1.0d0
      ParMat(1,2) = 0.0d0
      ParMat(2,1) = 0.0d0
      ParMat(2,2) = 1.0d0
    END IF

    SELECT CASE(Family)
    CASE(3)
      PRefElement = IsActivePElement(Element)
      SELECT CASE(ReductionMethod)
      CASE(CurlKernel)
        ! Method: the kernel of RT
        IP = GaussPoints(Element)
        DO t=1,IP % n

          IF (.NOT. PRefElement) THEN
            ! Switch to the p-reference element:
            uk = -1.0d0 + 2.0d0 * IP % U(t) + IP % V(t)
            vk = SQRT(3.0d0) * IP % V(t)
            sk = SQRT(3.0d0) * 2.0d0 * IP % S(t)
          ELSE
            uk = IP % U(t)
            vk = IP % V(t)
            sk = IP % S(t)
          END IF

          stat = ReductionOperatorInfo( Element, Nodes, uk, vk, StrainBasis, &
              ReductionMethod, ApplyPiolaTransform = .TRUE., Basis=Basis, DOFWeigths=DOFWeigths)  

          IF (UseParameters) THEN
            ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
            ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
            ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
            ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
          END IF

          DO i=1,nd
            DO j=1,n
              DO k=1,2
                SELECT CASE(k)
                CASE(1)
                  u(1) = ParMat(1,1)*Basis(j)
                  u(2) = ParMat(2,1)*Basis(j)
                CASE(2)
                  u(1) = ParMat(1,2)*Basis(j)
                  u(2) = ParMat(2,2)*Basis(j)
                END SELECT
                A(i,2*(j-1)+k) = A(i,2*(j-1)+k) + SUM(u(1:2) * DOFWeigths(i,1:2)) * sk
              END DO
            END DO
          END DO
        END DO

      CASE(MITC)
        ! MITC method: Create edge tangent vectors
        DO i=2,n
          Tau(i-1,1) = Nodes % x(i) - Nodes % x(i-1)
          Tau(i-1,2) = Nodes % y(i) - Nodes % y(i-1)
        END DO
        Tau(n,1) = Nodes % x(1) - Nodes % x(3)
        Tau(n,2) = Nodes % y(1) - Nodes % y(3)
        DO i=1,n
          Tau(i,3) = SQRT(SUM( Tau(i,1:2)**2 ))
          Tau(i,1:2) = Tau(i,1:2)/Tau(i,3)
        END DO

        IF (UseParameters) THEN
          IF (PRefElement) THEN
            uk = 0.0d0
            vk = 0.0d0           
          ELSE
            uk = 0.5d0
            vk = 0.0d0
          END IF
          stat = ElementInfo( Element, Nodes, uk, vk, 0.0d0, detJ, Basis )
          ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
          ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
          ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
          ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
        END IF

        A(1,1) = 0.5d0 * ( ParMat(1,1)*Tau(1,1) + ParMat(2,1)*Tau(1,2) ) * Tau(1,3)
        A(1,2) = 0.5d0 * ( ParMat(1,2)*Tau(1,1) + ParMat(2,2)*Tau(1,2) ) * Tau(1,3)
        A(1,3) = 0.5d0 * ( ParMat(1,1)*Tau(1,1) + ParMat(2,1)*Tau(1,2) ) * Tau(1,3)
        A(1,4) = 0.5d0 * ( ParMat(1,2)*Tau(1,1) + ParMat(2,2)*Tau(1,2) ) * Tau(1,3)

        IF (UseParameters) THEN
          IF (PRefElement) THEN
            uk = 0.5d0
            vk = sqrt(3.0d0)/2.0d0           
          ELSE
            uk = 0.5d0
            vk = 0.5d0
          END IF
          stat = ElementInfo( Element, Nodes, uk, vk, 0.0d0, detJ, Basis )
          ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
          ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
          ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
          ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
        END IF

        A(2,3) = 0.5d0 * ( ParMat(1,1)*Tau(2,1) + ParMat(2,1)*Tau(2,2) ) * Tau(2,3)
        A(2,4) = 0.5d0 * ( ParMat(1,2)*Tau(2,1) + ParMat(2,2)*Tau(2,2) ) * Tau(2,3)
        A(2,5) = 0.5d0 * ( ParMat(1,1)*Tau(2,1) + ParMat(2,1)*Tau(2,2) ) * Tau(2,3)
        A(2,6) = 0.5d0 * ( ParMat(1,2)*Tau(2,1) + ParMat(2,2)*Tau(2,2) ) * Tau(2,3)

        IF (UseParameters) THEN
          IF (PRefElement) THEN
            uk = -0.5d0
            vk = sqrt(3.0d0)/2.0d0           
          ELSE
            uk = 0.0d0
            vk = 0.5d0
          END IF
          stat = ElementInfo( Element, Nodes, uk, vk, 0.0d0, detJ, Basis )
          ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
          ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
          ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
          ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
        END IF

        A(3,5) = 0.5d0 * ( ParMat(1,1)*Tau(3,1) + ParMat(2,1)*Tau(3,2) ) * Tau(3,3)
        A(3,6) = 0.5d0 * ( ParMat(1,2)*Tau(3,1) + ParMat(2,2)*Tau(3,2) ) * Tau(3,3)
        A(3,1) = 0.5d0 * ( ParMat(1,1)*Tau(3,1) + ParMat(2,1)*Tau(3,2) ) * Tau(3,3)
        A(3,2) = 0.5d0 * ( ParMat(1,2)*Tau(3,1) + ParMat(2,2)*Tau(3,2) ) * Tau(3,3)

      CASE DEFAULT
        CALL Fatal('ReductionOperatorDOFs','Unknown strain reduction operator')
      END SELECT

    CASE(4)
      SELECT CASE(ReductionMethod)
      CASE(CurlKernel)
        ! Method: the kernel related to ABF
        IP = GaussPoints(Element)
        DO t=1,IP % n

          !stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          !    IP % W(t), detJ, Basis )

          IF (GradientOperand) THEN
            stat = ReductionOperatorInfo( Element, Nodes, IP % U(t), IP % V(t), StrainBasis, &
                ReductionMethod, ApplyPiolaTransform = .TRUE., Basis=Basis, DBasis=Dbasis, &
                DOFWeigths=DOFWeigths)
          ELSE
            stat = ReductionOperatorInfo( Element, Nodes, IP % U(t), IP % V(t), StrainBasis, &
                ReductionMethod, ApplyPiolaTransform = .TRUE., Basis=Basis, DOFWeigths=DOFWeigths)  
          END IF

          IF (UseParameters .AND. .NOT.GradientOperand) THEN
            ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
            ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
            ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
            ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
          END IF

          DO i=1,nd
            DO j=1,n
              DO k=1,2
                SELECT CASE(k)
                CASE(1)
                  IF (GradientOperand) THEN
                    u(1) = DBasis(j,1)
                    u(2) = 0.0d0
                  ELSE
                    u(1) = ParMat(1,1)*Basis(j)
                    u(2) = ParMat(2,1)*Basis(j)
                  END IF
                CASE(2)
                  IF (GradientOperand) THEN
                    u(1) = 0.0d0
                    u(2) = DBasis(j,2)
                  ELSE                  
                    u(1) = ParMat(1,2)*Basis(j)
                    u(2) = ParMat(2,2)*Basis(j)
                  END IF
                END SELECT
                A(i,2*(j-1)+k) = A(i,2*(j-1)+k) + SUM(u(1:2) * DOFWeigths(i,1:2)) * IP % s(t)
              END DO
            END DO
          END DO
        END DO

      CASE(MITC)
        ! MITC method: Create edge tangent vectors
        DO i=2,n
          Tau(i-1,1) = Nodes % x(i) - Nodes % x(i-1)
          Tau(i-1,2) = Nodes % y(i) - Nodes % y(i-1)
        END DO
        Tau(n,1) = Nodes % x(1) - Nodes % x(4)
        Tau(n,2) = Nodes % y(1) - Nodes % y(4)
        DO i=1,n
          Tau(i,3) = SQRT(SUM( Tau(i,1:2)**2 ))
          Tau(i,1:2) = Tau(i,1:2)/Tau(i,3)
        END DO

        IF (UseParameters) THEN
          stat = ElementInfo( Element, Nodes, 0.0d0, -1.0d0, 0.0d0, detJ, Basis )
          ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
          ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
          ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
          ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
        END IF

        A(1,1) = 0.5d0 * ( ParMat(1,1)*Tau(1,1) + ParMat(2,1)*Tau(1,2) ) * Tau(1,3)
        A(1,2) = 0.5d0 * ( ParMat(1,2)*Tau(1,1) + ParMat(2,2)*Tau(1,2) ) * Tau(1,3)
        A(1,3) = 0.5d0 * ( ParMat(1,1)*Tau(1,1) + ParMat(2,1)*Tau(1,2) ) * Tau(1,3)
        A(1,4) = 0.5d0 * ( ParMat(1,2)*Tau(1,1) + ParMat(2,2)*Tau(1,2) ) * Tau(1,3)

        IF (UseParameters) THEN
          stat = ElementInfo( Element, Nodes, 1.0d0, 0.0d0, 0.0d0, detJ, Basis )
          ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
          ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
          ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
          ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
        END IF

        A(2,3) = 0.5d0 * ( ParMat(1,1)*Tau(2,1) + ParMat(2,1)*Tau(2,2) ) * Tau(2,3)
        A(2,4) = 0.5d0 * ( ParMat(1,2)*Tau(2,1) + ParMat(2,2)*Tau(2,2) ) * Tau(2,3)
        A(2,5) = 0.5d0 * ( ParMat(1,1)*Tau(2,1) + ParMat(2,1)*Tau(2,2) ) * Tau(2,3)
        A(2,6) = 0.5d0 * ( ParMat(1,2)*Tau(2,1) + ParMat(2,2)*Tau(2,2) ) * Tau(2,3)

        IF (UseParameters) THEN
          stat = ElementInfo( Element, Nodes, 0.0d0, 1.0d0, 0.0d0, detJ, Basis )
          ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
          ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
          ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
          ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
        END IF

        A(3,5) = 0.5d0 * ( ParMat(1,1)*Tau(3,1) + ParMat(2,1)*Tau(3,2) ) * Tau(3,3)
        A(3,6) = 0.5d0 * ( ParMat(1,2)*Tau(3,1) + ParMat(2,2)*Tau(3,2) ) * Tau(3,3)
        A(3,7) = 0.5d0 * ( ParMat(1,1)*Tau(3,1) + ParMat(2,1)*Tau(3,2) ) * Tau(3,3)
        A(3,8) = 0.5d0 * ( ParMat(1,2)*Tau(3,1) + ParMat(2,2)*Tau(3,2) ) * Tau(3,3)

        IF (UseParameters) THEN
          stat = ElementInfo( Element, Nodes, -1.0d0, 0.0d0, 0.0d0, detJ, Basis )
          ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
          ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
          ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
          ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
        END IF

        A(4,7) = 0.5d0 * ( ParMat(1,1)*Tau(4,1) + ParMat(2,1)*Tau(4,2) ) * Tau(4,3)
        A(4,8) = 0.5d0 * ( ParMat(1,2)*Tau(4,1) + ParMat(2,2)*Tau(4,2) ) * Tau(4,3)
        A(4,1) = 0.5d0 * ( ParMat(1,1)*Tau(4,1) + ParMat(2,1)*Tau(4,2) ) * Tau(4,3)
        A(4,2) = 0.5d0 * ( ParMat(1,2)*Tau(4,1) + ParMat(2,2)*Tau(4,2) ) * Tau(4,3)

      CASE(CurlKernelWithEdgeDOFs)
        ! Method: the kernel related to ABF with edge/line integral DOFs
        E = 1
        IF (PRESENT(EdgeDirection)) THEN
          E = EdgeDirection
          IF (.NOT.(E==1 .OR. E==2)) &
              CALL Fatal('ReductionOperatorInfo','A wrong direction parameter given')
        END IF

        ! Create tangent vectors
        IF (E==1) THEN
          Tau(1,1) = Nodes % x(2) - Nodes % x(1)
          Tau(1,2) = Nodes % y(2) - Nodes % y(1)
          Tau(2,1) = Nodes % x(4) - Nodes % x(3)
          Tau(2,2) = Nodes % y(4) - Nodes % y(3)
          Tau(3,1) = (Nodes % x(4) + Nodes % x(3))/2.0d0 - (Nodes % x(2) + Nodes % x(1))/2.0d0
          Tau(3,2) = (Nodes % y(4) + Nodes % y(3))/2.0d0 - (Nodes % y(2) + Nodes % y(1))/2.0d0
        ELSE
          Tau(1,1) = Nodes % x(1) - Nodes % x(4)
          Tau(1,2) = Nodes % y(1) - Nodes % y(4)
          Tau(2,1) = Nodes % x(3) - Nodes % x(2)
          Tau(2,2) = Nodes % y(3) - Nodes % y(2)
          Tau(3,1) = (Nodes % x(3) + Nodes % x(2))/2.0d0 - (Nodes % x(4) + Nodes % x(1))/2.0d0
          Tau(3,2) = (Nodes % y(3) + Nodes % y(2))/2.0d0 - (Nodes % y(4) + Nodes % y(1))/2.0d0
        END IF
        DO i=1,3
          Tau(i,3) = SQRT(SUM( Tau(i,1:2)**2 ))
          Tau(i,1:2) = Tau(i,1:2)/Tau(i,3)
        END DO

        IF (E==1) THEN
          IF (UseParameters .OR. GradientOperand) THEN
            stat = ElementInfo( Element, Nodes, 0.0d0, -1.0d0, 0.0d0, detJ, Basis, DBasis )
            IF (.NOT. GradientOperand) THEN
              ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
              ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
              ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
              ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
            END IF
          END IF

          IF (GradientOperand) THEN
            A(1,1) = DBasis(1,1) * Tau(1,1) * Tau(1,3)
            A(1,2) = DBasis(1,2) * Tau(1,2) * Tau(1,3)
            A(1,3) = DBasis(2,1) * Tau(1,1) * Tau(1,3)
            A(1,4) = DBasis(2,2) * Tau(1,2) * Tau(1,3)            
          ELSE
            A(1,1) = 0.5d0 * ( ParMat(1,1)*Tau(1,1) + ParMat(2,1)*Tau(1,2) ) * Tau(1,3)
            A(1,2) = 0.5d0 * ( ParMat(1,2)*Tau(1,1) + ParMat(2,2)*Tau(1,2) ) * Tau(1,3)
            A(1,3) = 0.5d0 * ( ParMat(1,1)*Tau(1,1) + ParMat(2,1)*Tau(1,2) ) * Tau(1,3)
            A(1,4) = 0.5d0 * ( ParMat(1,2)*Tau(1,1) + ParMat(2,2)*Tau(1,2) ) * Tau(1,3)
          END IF

          IF (UseParameters .OR. GradientOperand) THEN
            stat = ElementInfo( Element, Nodes, 0.0d0, 1.0d0, 0.0d0, detJ, Basis, DBasis )
            IF (.NOT. GradientOperand) THEN
              ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
              ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
              ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
              ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
            END IF
          END IF

          IF (GradientOperand) THEN
            A(2,5) = DBasis(3,1) * Tau(2,1) * Tau(2,3)
            A(2,6) = DBasis(3,2) * Tau(2,2) * Tau(2,3)
            A(2,7) = DBasis(4,1) * Tau(2,1) * Tau(2,3)
            A(2,8) = DBasis(4,2) * Tau(2,2) * Tau(2,3)
          ELSE
            A(2,5) = 0.5d0 * ( ParMat(1,1)*Tau(2,1) + ParMat(2,1)*Tau(2,2) ) * Tau(2,3)
            A(2,6) = 0.5d0 * ( ParMat(1,2)*Tau(2,1) + ParMat(2,2)*Tau(2,2) ) * Tau(2,3)
            A(2,7) = 0.5d0 * ( ParMat(1,1)*Tau(2,1) + ParMat(2,1)*Tau(2,2) ) * Tau(2,3)
            A(2,8) = 0.5d0 * ( ParMat(1,2)*Tau(2,1) + ParMat(2,2)*Tau(2,2) ) * Tau(2,3)
          END IF
        ELSE
          IF (UseParameters .OR. GradientOperand) THEN
            stat = ElementInfo( Element, Nodes, -1.0d0, 0.0d0, 0.0d0, detJ, Basis, DBasis )
            IF (.NOT. GradientOperand) THEN
              ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
              ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
              ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
              ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
            END IF
          END IF

          IF (GradientOperand) THEN
            A(1,7) = DBasis(4,1) * Tau(1,1) * Tau(1,3)
            A(1,8) = DBasis(4,2) * Tau(1,2) * Tau(1,3)
            A(1,1) = DBasis(1,1) * Tau(1,1) * Tau(1,3)
            A(1,2) = DBasis(1,2) * Tau(1,2) * Tau(1,3)
          ELSE
            A(1,7) = 0.5d0 * ( ParMat(1,1)*Tau(1,1) + ParMat(2,1)*Tau(1,2) ) * Tau(1,3)
            A(1,8) = 0.5d0 * ( ParMat(1,2)*Tau(1,1) + ParMat(2,2)*Tau(1,2) ) * Tau(1,3)
            A(1,1) = 0.5d0 * ( ParMat(1,1)*Tau(1,1) + ParMat(2,1)*Tau(1,2) ) * Tau(1,3)
            A(1,2) = 0.5d0 * ( ParMat(1,2)*Tau(1,1) + ParMat(2,2)*Tau(1,2) ) * Tau(1,3)
          END IF

          IF (UseParameters .OR. GradientOperand) THEN
            stat = ElementInfo( Element, Nodes, 1.0d0, 0.0d0, 0.0d0, detJ, Basis, DBasis )
            IF (.NOT. GradientOperand) THEN
              ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
              ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
              ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
              ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
            END IF
          END IF

          IF (GradientOperand) THEN
            A(2,3) = DBasis(2,1) * Tau(2,1) * Tau(2,3)
            A(2,4) = DBasis(2,2) * Tau(2,2) * Tau(2,3)
            A(2,5) = DBasis(3,1) * Tau(2,1) * Tau(2,3)
            A(2,6) = DBasis(3,2) * Tau(2,2) * Tau(2,3)
          ELSE
            A(2,3) = 0.5d0 * ( ParMat(1,1)*Tau(2,1) + ParMat(2,1)*Tau(2,2) ) * Tau(2,3)
            A(2,4) = 0.5d0 * ( ParMat(1,2)*Tau(2,1) + ParMat(2,2)*Tau(2,2) ) * Tau(2,3)
            A(2,5) = 0.5d0 * ( ParMat(1,1)*Tau(2,1) + ParMat(2,1)*Tau(2,2) ) * Tau(2,3)
            A(2,6) = 0.5d0 * ( ParMat(1,2)*Tau(2,1) + ParMat(2,2)*Tau(2,2) ) * Tau(2,3)
          END IF
        END IF

        IF (UseParameters .OR. GradientOperand) THEN
          stat = ElementInfo( Element, Nodes, 0.0d0, 0.0d0, 0.0d0, detJ, Basis, DBasis )
          IF (.NOT. GradientOperand) THEN         
            ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
            ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
            ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
            ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
          END IF
        END IF

        IF (GradientOperand) THEN
          A(3,1) = DBasis(1,1) * Tau(3,1) * Tau(3,3)
          A(3,2) = DBasis(1,2) * Tau(3,2) * Tau(3,3)
          A(3,3) = DBasis(2,1) * Tau(3,1) * Tau(3,3)
          A(3,4) = DBasis(2,2) * Tau(3,2) * Tau(3,3)

          A(3,5) = DBasis(3,1) * Tau(3,1) * Tau(3,3)
          A(3,6) = DBasis(3,2) * Tau(3,2) * Tau(3,3)
          A(3,7) = DBasis(4,1) * Tau(3,1) * Tau(3,3)
          A(3,8) = DBasis(4,2) * Tau(3,2) * Tau(3,3)
        ELSE
          A(3,1) = 0.25d0 * ( ParMat(1,1)*Tau(3,1) + ParMat(2,1)*Tau(3,2) ) * Tau(3,3)
          A(3,2) = 0.25d0 * ( ParMat(1,2)*Tau(3,1) + ParMat(2,2)*Tau(3,2) ) * Tau(3,3)
          A(3,3) = 0.25d0 * ( ParMat(1,1)*Tau(3,1) + ParMat(2,1)*Tau(3,2) ) * Tau(3,3)
          A(3,4) = 0.25d0 * ( ParMat(1,2)*Tau(3,1) + ParMat(2,2)*Tau(3,2) ) * Tau(3,3)

          A(3,5) = 0.25d0 * ( ParMat(1,1)*Tau(3,1) + ParMat(2,1)*Tau(3,2) ) * Tau(3,3)
          A(3,6) = 0.25d0 * ( ParMat(1,2)*Tau(3,1) + ParMat(2,2)*Tau(3,2) ) * Tau(3,3)
          A(3,7) = 0.25d0 * ( ParMat(1,1)*Tau(3,1) + ParMat(2,1)*Tau(3,2) ) * Tau(3,3)
          A(3,8) = 0.25d0 * ( ParMat(1,2)*Tau(3,1) + ParMat(2,2)*Tau(3,2) ) * Tau(3,3)
        END IF

      CASE DEFAULT
        CALL Fatal('ReductionOperatorDOFs','Unknown strain reduction operator')
      END SELECT
    END SELECT
!------------------------------------------------------------------------------
  END SUBROUTINE ReductionOperatorDofs
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Compute the matrix A which can be used to evaluate DOFs for the interpolating
!  function R_K(u) in the strain reduction space X(K) when applied to a H1-
!  conforming bubble function u having two components. If d_k denotes the linear 
!  functional defining the kth DOF, the kth row of the returned matrix has 
!  entries 
!
!     [d_k(Nb*e1) d_k(Nb*e2)]
!
!  where Nb is the bubble basis function, e1=(1,0) and e2=(0,1). Optionally 
!  the interpolating function can be computed for a field Cu where C is a 2X2 matrix 
!  field. Currently, just one bubble function (cf. the size of Basis array) and 
!  the kernel version of strain reduction are supported currently.
!------------------------------------------------------------------------------
  SUBROUTINE ReductionOperatorBubbleDofs(Element, Nodes, A, nd, nb, n, ReductionMethod, &
      ModelPars, GradientField)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), INTENT(IN), TARGET :: Element          !< Element structure
    TYPE(Nodes_t), INTENT(IN) :: Nodes                      !< Nodes structure
    REAL(KIND=dp), INTENT(INOUT) :: A(nd,2*nb)              !< Coefficients for expressing the DOFs 
    INTEGER, INTENT(IN) :: nd                               !< The dimension of the strain reduction space X(K)
    INTEGER, INTENT(IN) :: nb                               !< The number of the H1-conforming bubble functions
    INTEGER, INTENT(IN) :: n                                !< The number of the BG element nodes
    INTEGER, INTENT(IN) :: ReductionMethod                  !< The method chosen
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: ModelPars(2,2,n) !< To include the effect of additional model parameters
    LOGICAL, OPTIONAL :: GradientField                      !< To return [d_k(grad(Nb).e1) d_k(grad(Nb).e2)]
!---------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: stat, UseParameters, PRefElement, GradientOperand

    INTEGER :: Family, i, j, k, t

    REAL(KIND=dp) :: StrainBasis(4,3)         ! The basis functions for the strain reduction space X(K)
    REAL(KIND=dp) :: DOFWeigths(3,2)          ! The auxiliary functions to evaluate the interpolant in X(K)
    REAL(KIND=dp) :: Basis(5)                 ! H1-conforming basis functions (p=1, with one bubble)
    REAL(KIND=dp) :: DBasis(5,1:3) 
    REAL(KIND=dp) :: u(2), ParMat(2,2), uk, vk, sk
!---------------------------------------------------------------------------------
    IF (ReductionMethod /= CurlKernel) CALL Fatal('ReductionOperatorBubbleDofs', &
        'An unsupported strain reduction technique')

    IF (PRESENT(GradientField)) THEN
      GradientOperand = GradientField
    ELSE
      GradientOperand = .FALSE.
    END IF

    Family = GetElementFamily(Element)
    
    A = 0.0d0

    UseParameters = PRESENT(ModelPars)
    IF (.NOT. UseParameters) THEN
      ParMat(1,1) = 1.0d0
      ParMat(1,2) = 0.0d0
      ParMat(2,1) = 0.0d0
      ParMat(2,2) = 1.0d0
    END IF

    IP = GaussPoints(Element)
    DO t=1,IP % n

      PRefElement = IsActivePElement(Element)
      
      IF (Family == 3 .AND. .NOT.PRefElement) THEN
        ! Switch to the p-reference element:
        uk = -1.0d0 + 2.0d0 * IP % U(t) + IP % V(t)
        vk = SQRT(3.0d0) * IP % V(t)
        sk = SQRT(3.0d0) * 2.0d0 * IP % S(t)
      ELSE
        uk = IP % U(t)
        vk = IP % V(t)
        sk = IP % S(t)
      END IF

      stat = ReductionOperatorInfo(Element, Nodes, uk, vk, StrainBasis, &
          ReductionMethod, ApplyPiolaTransform = .TRUE., Basis=Basis, DBasis=DBasis, &
          DOFWeigths=DOFWeigths, Bubbles=.TRUE.)

      IF (UseParameters) THEN
        ParMat(1,1) = SUM(ModelPars(1,1,1:n) * Basis(1:n))
        ParMat(1,2) = SUM(ModelPars(1,2,1:n) * Basis(1:n))
        ParMat(2,1) = SUM(ModelPars(2,1,1:n) * Basis(1:n))
        ParMat(2,2) = SUM(ModelPars(2,2,1:n) * Basis(1:n))
      END IF

      DO i=1,nd
        DO j=1,nb
          DO k=1,2
            SELECT CASE(k)
            CASE(1)
              IF (GradientOperand) THEN
                u(1) = DBasis(n+j,1)
                u(2) = 0.0d0
              ELSE
                u(1) = ParMat(1,1)*Basis(n+j)
                u(2) = ParMat(2,1)*Basis(n+j)
              END IF
            CASE(2)
              IF (GradientOperand) THEN
                u(1) = 0.0d0
                u(2) = DBasis(n+j,2)
              ELSE 
                u(1) = ParMat(1,2)*Basis(n+j)
                u(2) = ParMat(2,2)*Basis(n+j)
              END IF
            END SELECT
            A(i,2*(j-1)+k) = A(i,2*(j-1)+k) + SUM(u(1:2) * DOFWeigths(i,1:2)) * sk
          END DO
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReductionOperatorBubbleDofs
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Compute the average of the nodal director data saved as elementwise property
! 'director' over n-node element. Optionally check whether the surface is
! planar.  
!------------------------------------------------------------------------------
  FUNCTION AverageDirector(Element, n, PlanarSurface) RESULT(d)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n
    LOGICAL, OPTIONAL, INTENT(OUT) :: PlanarSurface
    REAL(KIND=dp) :: d(3)
!------------------------------------------------------------------------------
    LOGICAL :: PlateBody, PlateBodyCheck
    INTEGER :: i, i0
    REAL(KIND=dp), POINTER :: DirectorValues(:)
    REAL(KIND=dp) :: v(3)
!------------------------------------------------------------------------------
    d = 0.0d0

    IF (PRESENT(PlanarSurface)) THEN
      PlateBodyCheck = .TRUE.
    ELSE
      PlateBodyCheck = .FALSE.
    END IF

    DirectorValues => GetElementalDirector( Element )
    
    IF (ASSOCIATED(DirectorValues)) THEN
      IF (SIZE(DirectorValues) < 3*n) CALL Fatal('AverageDirector', &
          'Elemental director data is not associated with all nodes')
    ELSE
      CALL Fatal('AverageDirector', 'Elemental director data is not associated')
    END IF

    DO i=1,n
      i0 = (i-1)*3
      d(1:3) = d(1:3) + DirectorValues(i0+1:i0+3)
    END DO
    d(:) = 1.0d0/n * d(:)

    IF (PlateBodyCheck) THEN
      v(1:3) = DirectorValues(1:3)
      PlateBody = .TRUE.
      DO i=2,n
        i0 = (i-1)*3
        PlateBody = .NOT. MAXVAL(ABS(v(1:3) - DirectorValues(i0+1:i0+3))) > EPSILON(1.0)
        IF (.NOT. PlateBody) EXIT
      END DO
      PlanarSurface = PlateBody
    END IF
!------------------------------------------------------------------------------
  END FUNCTION AverageDirector
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Compute the area of an element of the mapped background mesh and add to the 
! total value
!------------------------------------------------------------------------------
  SUBROUTINE MappedBGMeshArea(Element, LocalFrameNodes, Area)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    REAL(KIND=dp), TARGET, INTENT(IN) :: LocalFrameNodes(MaxPatchNodes,3)
    REAL(KIND=dp), INTENT(INOUT) :: Area
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: PlaneElement => NULL()
    TYPE(Nodes_t) :: NodesVar
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: stat
    INTEGER :: j, Family
    REAL(KIND=dp) :: Basis(MaxPatchNodes), detJ

    SAVE PlaneElement
!------------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED(PlaneElement) ) PlaneElement => AllocateElement()
    
    Family = GetElementFamily(Element)
    SELECT CASE(Family)
    CASE(3)
      PlaneElement % Type => GetElementType(310, .FALSE.)
    CASE(4)
      PlaneElement % Type => GetElementType(416, .FALSE.)
    CASE DEFAULT
      RETURN
    END SELECT

    NodesVar % x => LocalFrameNodes(:,1)
    NodesVar % y => LocalFrameNodes(:,2)
    NodesVar % z => LocalFrameNodes(:,3)
    NodesVar % NumberOfNodes = PlaneElement % Type % NumberOfNodes

    IP = GaussPoints(PlaneElement)

    DO j=1,IP % n   
      stat = ElementInfo(PlaneElement, NodesVar, IP % u(j), IP % v(j), &
          IP % w(j), detJ, Basis)
      Area = Area + IP % s(j) * detJ
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE MappedBGMeshArea
!------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
! Compute the area of the blending element surface by calling the function
! BlendingSurfaceInfo that defines the blending surface in the first place
!--------------------------------------------------------------------------------
  SUBROUTINE ComputeSurfaceArea(Element, SurfaceArea)
!--------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    REAL(KIND=dp), INTENT(INOUT) :: SurfaceArea
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: Stat
    INTEGER :: j, Family
    REAL(KIND=dp) :: a1(3), a2(3), a3(3), A(2,2), B(2,2), x(3), DetA
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    Family = GetElementFamily(Element)

    ! -------------------------------------------------------
    ! Numerical integration rule to compute the surface area:
    ! -------------------------------------------------------
    SELECT CASE(Family)
    CASE(3)
      IP = GaussPointsTriangle(11, PReferenceElement=.TRUE.)
    CASE(4)
      IP = GaussPoints(Element, 25)
    CASE DEFAULT
      RETURN
    END SELECT

    DO j=1,IP % n
      stat = BlendingSurfaceInfo(Element, Nodes, IP % U(j), IP % V(j), &
          DetA, a1, a2, a3, A, B, x)      
      SurfaceArea = SurfaceArea + IP % s(j) * SQRT(Deta)
    END DO
!-------------------------------------------------------------------------------------
  END SUBROUTINE ComputeSurfaceArea
!-------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------
! Return the global coordinates at the mid-node of the edge by evaluating
! the value of the space curve. This function is just for testing purposes.
!-------------------------------------------------------------------------------------
  FUNCTION EdgeMidNode(Element, e) RESULT(X)
!-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: e     ! Edge identifier 
    REAL(KIND=dp) :: X(3)        ! Global coordinates at the mid-node of the edge 
!-----------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes
    INTEGER :: CurveDataSize, i0, cn
    REAL(KIND=dp), POINTER :: EdgeParams(:)
    REAL(KIND=dp) :: HermBasis(6), dHermBasis(6), ddHermBasis(6)
    REAL(KIND=dp) :: d(CurveDataSize2), h, xe
    REAL(KIND=dp) :: r1(3), r2(3)
    REAL(KIND=dp) :: u, v, f
!-----------------------------------------------------------------------
    IF (Element % Type % NumberOfNodes > 4) &
        CALL Fatal('EdgeMidNode', 'Just 3-node and 4-node elements implemented')

    !-----------------------------------------------------------------------
    ! Retrieve parametrizations of curved edges:
    !------------------------------------------------------------------------
    EdgeParams => GetElementProperty('edge parameters', Element)

    h = 2.0d0
    CurveDataSize = CurveDataSize1

    i0 = (e-1)*CurveDataSize
    d(1:CurveDataSize) = EdgeParams(i0+1:i0+CurveDataSize)

    cn = 2
    CALL HermiteBasis(0.0d0, h, HermBasis(1:2*cn), dHermBasis(1:2*cn), ddHermBasis(1:2*cn), cn)

    CALL GetElementNodes(Nodes, Element) 

    Family = GetElementFamily(Element)
    SELECT CASE(Family)
    CASE(3)
      SELECT CASE(e)
      CASE(1)
        r1(1) = Nodes % x(1)
        r1(2) = Nodes % y(1)
        r1(3) = Nodes % z(1)
        r2(1) = Nodes % x(2)
        r2(2) = Nodes % y(2)
        r2(3) = Nodes % z(2)
      CASE(2)
        r1(1) = Nodes % x(2)
        r1(2) = Nodes % y(2)
        r1(3) = Nodes % z(2)
        r2(1) = Nodes % x(3)
        r2(2) = Nodes % y(3)
        r2(3) = Nodes % z(3)
      CASE(3)
        r1(1) = Nodes % x(3)
        r1(2) = Nodes % y(3)
        r1(3) = Nodes % z(3)
        r2(1) = Nodes % x(1)
        r2(2) = Nodes % y(1)
        r2(3) = Nodes % z(1)
      END SELECT
    CASE(4)
      CALL Fatal('EdgeMidNode', '4-node implementation missing')
    END SELECT

    X(1:3) = r1(1:3)*HermBasis(1) + r2(1:3)*HermBasis(2) + &
        d(1:3)*0.5d0*HermBasis(3) + d(4:6)*0.5d0*HermBasis(4)

!-----------------------------------------------------------------------
  END FUNCTION EdgeMidNode
!-----------------------------------------------------------------------

!------------------------------------------------------------------------------
! The matrix representation of the 3D elasticity tensor with respect an orthogonal
! basis. The case A1 = A2 = 1 corresponds to an orthonormal basis.
!------------------------------------------------------------------------------
  SUBROUTINE MaterialLaw3D(CMat, A1, A2, E, nu)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(OUT) :: CMat(6,6)
    REAL(KIND=dp), INTENT(IN) :: A1, A2, E, nu
!------------------------------------------------------------------------------
    INTEGER :: i
    REAL(KIND=dp) :: Lambda, G
!------------------------------------------------------------------------------
    Lambda = E * nu / ( (1.0_dp + nu) * (1.0_dp - 2.0_dp*nu) )
    G = E / (2.0_dp*(1.0_dp + nu))

    ! We use unusual ordering to conform with the earlier implementation:
    ! 11, 22, 12, 33, 13, 23
    
    CMat = 0.0_dp
    CMat(1:2,1:2) = Lambda
    CMat(1:2,4) = Lambda
    CMat(4,1:2) = Lambda
    CMat(4,4) = Lambda

    DO i=1,2
      CMat(i,i) = CMat(i,i) + 2.0_dp * G
      CMat(i+4,i+4) = G
    END DO
    CMat(4,4) = CMat(4,4) + 2.0_dp * G
    CMat(3,3) = G

    CMat(1,1) = CMat(1,1)/A1**4
    CMat(1,2) = CMat(1,2)/(A1**2 * A2**2)
    CMat(1,4) = CMat(1,4)/(A1**2)

    CMat(2,1) = CMat(2,1)/(A2**2 * A1**2)
    CMat(2,2) = CMat(2,2)/A2**4   
    CMat(2,4) = CMat(2,4)/(A2**2)

    CMat(4,1) = CMat(4,1)/(A1**2)
    CMat(4,2) = CMat(4,2)/(A2**2)
    CMat(3,3) = CMat(3,3)/(A1**2 * A2**2)

    CMat(5,5) = CMat(5,5)/(A1**2)
    CMat(6,6) = CMat(6,6)/(A2**2)
!------------------------------------------------------------------------------
  END SUBROUTINE MaterialLaw3D
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SurfaceBasis(y1, y2, A, K1, K2, Spherical, Cylindrical)
!------------------------------------------------------------------------------
! Returns the surface basis vectors in the case of some example geometries
! (a part of a sphere, cylinder or plate). The kth basis vector is represented 
! as the column vector A(:,k). This subroutine could be generalized to other
! cases which admit a parametrization by lines of curvature coordinates.
!------------------------------------------------------------------------------    
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: y1, y2  ! Curvilinear coordinates on a surface
    REAL(KIND=dp), INTENT(OUT) :: A(3,3) 
    REAL(KIND=dp), INTENT(OUT) :: K1, K2 ! The principal curvatures Ki = b_ii/a_ii
    LOGICAL, INTENT(IN) :: Spherical
    LOGICAL, INTENT(IN) :: Cylindrical

!------------------------------------------------------------------------------
    REAL(KIND=dp) :: R
!------------------------------------------------------------------------------
    A = 0.0_dp

    IF (Spherical) THEN
      !
      ! Using spherical coordinates
      !
      R = 10.0_dp  ! The radius

      A(1,1) = -R*cos(y2)*sin(y1)
      A(2,1) = R*cos(y2)*cos(y1)

      A(1,2) = -R*sin(y2)*cos(y1)
      A(2,2) = -R*sin(y2)*sin(y1)
      A(3,2) = R*cos(y2)

      A(1,3) = cos(y2)*cos(y1)
      A(2,3) = cos(y2)*sin(y1)
      A(3,3) = sin(y2)

      K1 = -1.0_dp/R
      K2 = -1.0_dp/R

    ELSE IF (Cylindrical) THEN
      !
      ! A straight cylinder:
      !
      R = 1.0_dp

      A(1,1) = cos(y1/R)
      A(3,1) = sin(y1/R)

      A(2,2) = 1.0_dp 

      A(1,3) = -sin(y1/R)
      A(3,3) =  cos(y1/R)

      K1 = 1.0_dp/R
      K2 = 0.0_dp

    ELSE
      !
      ! Planar body:
      !
      A = 0.0_dp
      A(1,1) = 1.0_dp
      A(2,2) = 1.0_dp
      A(3,3) = 1.0_dp
      
      K1 = 0.0_dp
      K2 = 0.0_dp
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE SurfaceBasis
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! This subroutine for elementwise assembly is based on a Cartesian components
! formulation which does not use the traditional evaluation of covariant
! derivatives. The full quadratic expansions in powers of the thickness coordinate
! may also be employed in the approximation of displacement when 9 DOFs are used, 
! but the 6-DOF alternative is also available. In addition, approximation with
! p-elements is supported. The surface model can be based on a physical surface
! mesh for which an improved surface reconstruction has been done. An alternative
! is to employ the classical approach where a chart is given as initial data.
! Then the discretization is performed over a mathematical domain of two
! dimensions; see the subroutine SurfaceBasis which defines the chart.
!
! TO DO: Consider moving the subroutine SurfaceBasis elsewhere so that it can be
!        replaced easily by a user-supplied subroutine without modifying this file. 
!------------------------------------------------------------------------------
  SUBROUTINE ShellLocalMatrixCartesian(BGElement, n, nd, m, LocalSol, LargeDeflection, &
      NonlinearBending, MassAssembly, HarmonicAssembly, RHSForce, SkipBlending, &
      ReparametrizeMesh, BenchmarkProblem)
!------------------------------------------------------------------------------
    USE SolidMechanicsUtils, ONLY: StrainEnergyDensity, ShearCorrectionFactor
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: BGElement  ! An element of background mesh
    INTEGER, INTENT(IN) :: n                           ! The number of background element nodes
    INTEGER, INTENT(IN) :: nd                          ! The number of DOFs per component (after
                                                       ! static condensation if bubbles are used)
    INTEGER, INTENT(IN) :: m                           ! The number of DOFs per node
    REAL(KIND=dp), INTENT(IN) :: LocalSol(:,:)         ! The previous solution iterate
    LOGICAL, INTENT(IN) :: LargeDeflection             ! To activate nonlinear terms
    LOGICAL, INTENT(IN) :: NonlinearBending            ! The full nonlinearity of bending and transverse shear strains
    LOGICAL, INTENT(IN) :: MassAssembly                ! To activate mass matrix integration
    LOGICAL, INTENT(IN) :: HarmonicAssembly            ! To activate the global mass matrix updates
    REAL(KIND=dp), INTENT(OUT) :: RHSForce(:)          ! Local RHS vector corresponding to external loads
    LOGICAL, INTENT(IN) :: SkipBlending                ! Informs whether surface reconstruction has been done
    LOGICAL, INTENT(IN) :: ReparametrizeMesh           ! To use an orthogonal parametrization created without blending
    LOGICAL, INTENT(IN), OPTIONAL :: BenchmarkProblem  ! To create a load for a benchmark problem 
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: GElement => NULL()
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueList_t), POINTER :: BodyForce, Material, BodyParams
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: Stat, Found, ApplyNormalPressure
    LOGICAL :: TransverseBendingStretch, Spherical, Cylindrical
    LOGICAL :: GeneralMaterial
    LOGICAL :: PlateBody                   ! To indicate that the surface is flat
    LOGICAL :: SphericalSurface            ! To indicate that the surface is considered to spherical

    INTEGER :: DOFs, i, j, k, p, t, csize, bsize

    REAL(KIND=dp), PARAMETER :: i1(3) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
    REAL(KIND=dp), PARAMETER :: i2(3) = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
    REAL(KIND=dp), PARAMETER :: i3(3) = (/ 0.0_dp, 0.0_dp, 1.0_dp /)

    REAL(KIND=dp), POINTER :: TaylorParams(:)

    REAL(KIND=dp) :: Stiff(m*nd,m*nd), Mass(m*nd,m*nd), Force(m*nd)
    REAL(KIND=dp) :: Damp(m*nd,m*nd)
    REAL(KIND=dp) :: BM(4,m*nd), BS(4,m*nd), BB(4,m*nd)
    REAL(KIND=dp) :: NonlinBM(4,m*nd), NonlinBS(4,m*nd), NonlinBB(4,m*nd)
    REAL(KIND=dp) :: Basis(nd), dBasis(nd,3)
    REAL(KIND=dp) :: StrainVec(6), StressVec(6)
    REAL(KIND=dp) :: PrevSolVec(m*nd)
    REAL(KIND=dp) :: CMat(4,4), GMat(2,2), HMat(6,6)
    REAL(KIND=dp) :: A11, A22, SqrtDetA, A1, A2
    REAL(KIND=dp) :: B11, B22   
    REAL(KIND=dp) :: C111, C112, C221, C222, C211, C212
    REAL(KIND=dp) :: abasis1(3), abasis2(3), abasis3(3)
    REAL(KIND=dp) :: abasis1New(3), abasis2New(3), abasis3New(3), NewDetA
    REAL(KIND=dp) :: y1, y2, v1, v2, v3
    REAL(KIND=dp) :: sq
    REAL(KIND=dp) :: PoissonRatio(n), YoungsMod(n), ShellThickness(n), Load(n), rho(n), rho0
    REAL(KIND=dp) :: Damping(n), DampCoef
    REAL(KIND=dp) :: nu, E, h, NormalTraction
    REAL(KIND=dp) :: DetJ, Weight, Norm
    REAL(KIND=dp) :: CovariantBasis(3,3), K1, K2
    REAL(KIND=dp) :: Q(3,3)
    REAL(KIND=dp) :: PrevGrad(3,2), PrevGradB(3,2), PrevB(3), PrevPsi(3)
    REAL(KIND=dp) :: c1, c2, c3, c4
    REAL(KIND=dp) :: PatchNodes(MaxPatchNodes,2) ! The nodes of principal coordinate patch
    REAL(KIND=dp) :: e1(3), e2(3), e3(3)         ! The basis of the local frame
    REAL(KIND=dp) :: o(3)                        ! The origin of the local frame

    SAVE Nodes, GElement
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes, BGElement)

    IF (.NOT. SkipBlending .OR. ReparametrizeMesh) THEN
      ! ------------------------------------------------------------------------------
      ! When the elementwise reparametrization has been done, retrieve the surface data
      ! which have been saved as elementwise properties:
      ! ------------------------------------------------------------------------------
      CALL RetrieveLocalFrame(BGElement, TaylorParams, PatchNodes, e1, e2, e3, &
          o, PlateBody, SphericalSurface, GElement)

      ! --------------------------------------------------------------------------
      ! Overwrite the coordinate arrays of the structure Nodes so that it represents 
      ! the domain of the principal curvature coordinates:
      ! --------------------------------------------------------------------------
      CALL SolveNodesVariables(BGElement, Nodes, nd, GElement, PatchNodes)
    END IF

    Material => GetMaterial()
    GeneralMaterial = GetLogical(Material, '3D Material Law', Found)
!    IF (GeneralMaterial) print *, 'APPLYING 3D MATERIAL LAW'

    BodyParams => GetBodyParams(BGElement)
    Spherical = GetLogical(BodyParams, 'Spherical Body', Found)
    Cylindrical = GetLogical(BodyParams, 'Cylindrical Body', Found)

    ! ------------------------------------------------------------------------------
    ! The size of the constitutive matrix for 2D shell equations
    ! ------------------------------------------------------------------------------
    TransverseBendingStretch = m == 9
    IF (GeneralMaterial .AND. .NOT. TransverseBendingStretch) &
        CALL Fatal('ShellSolver', '3D material law needs 9-field model') 

    IF (TransverseBendingStretch) THEN
      bsize = 4
    ELSE
      IF (m /= 6) CALL Fatal('ShellLocalMatrixCartesian', &
          '6 or 9 Variable DOFs expected')
      bsize = 3
    END IF
    csize = 4

    DOFs = m*nd ! The local stiffness matrix size after static condensation

    ! --------------------------------------------------------------------------
    ! Body forces, material parameters and the shell thickness:
    ! --------------------------------------------------------------------------
    PoissonRatio(1:n) = GetReal(Material, 'Poisson Ratio')
    YoungsMod(1:n) = GetReal(Material, 'Youngs Modulus')
    ShellThickness(1:n) = GetReal(Material, 'Shell Thickness', Found)
    IF (.NOT. Found) THEN
      ShellThickness(1:n) = GetReal(Material, 'Thickness', Found)
      IF (.NOT. Found) CALL Fatal('ShellLocalMatrixCartesian', 'Shell Thickness/Thickness undefined')
    END IF

    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) THEN
      Load(1:n) = GetReal(BodyForce, 'Normal Pressure', ApplyNormalPressure)
    ELSE
      Load(1:n) = 0.0d0
      ApplyNormalPressure = .FALSE.
    END IF
    IF ( MassAssembly ) THEN
      rho(1:n) = GetReal(Material, 'Density')
      Damping(1:n) = GetReal(Material, 'Rayleigh Damping Alpha', Found)
    END IF

    ! ------------------------------------------------------------------------
    ! Vectorize the previous solution
    ! ------------------------------------------------------------------------
    DO k=1,m
      PrevSolVec(k:DOFs:m) = LocalSol(k,1:nd)
    END DO

    ! ------------------------------------------------------------------------
    ! Finally, integrate local element matrices:
    ! ------------------------------------------------------------------------
    Mass = 0.0d0
    Damp = 0.0d0
    Stiff = 0.0d0
    Force = 0.0d0
    RHSForce = 0.0d0

    IP = GaussPoints( BGElement )

    QUADRATURELOOP: DO t=1,IP % n

      BM = 0.0d0
      BB = 0.0d0
      BS = 0.0d0

      NonlinBM = 0.0d0
      NonlinBS = 0.0d0
      NonlinBB = 0.0d0

      stat = ElementInfo(BGElement, Nodes, IP % U(t), IP % V(t), IP % W(t), detJ, Basis, dBasis)
      sq = IP % s(t)

      ! ------------------------------------------------
      ! Data interpolation:
      ! ------------------------------------------------
      h = SUM( ShellThickness(1:n) * Basis(1:n) )
      nu = SUM( PoissonRatio(1:n) * Basis(1:n) )
      E = SUM( YoungsMod(1:n) * Basis(1:n) )
      IF (ApplyNormalPressure) NormalTraction = SUM( Load(1:n) * Basis(1:n) )
      IF ( MassAssembly ) THEN
        rho0 = SUM( rho(1:n) * Basis(1:n) )
        DampCoef = SUM( Damping(1:n) * Basis(1:n) )
      END IF

      IF (SkipBlending .AND. .NOT. ReparametrizeMesh) THEN
        !
        ! This uses a surface parametrization over a 2D domain in special cases.
        !
        y1 = SUM( Nodes % x(1:n) * Basis(1:n) )
        y2 = SUM( Nodes % y(1:n) * Basis(1:n) )

        CALL SurfaceBasis(y1, y2, CovariantBasis, K1, K2, Spherical, Cylindrical)

        abasis1(:) = CovariantBasis(:,1)
        abasis2(:) = CovariantBasis(:,2)
        abasis3(:) = CovariantBasis(:,3)
        a11 = DOT_PRODUCT(abasis1,abasis1)
        a22 = DOT_PRODUCT(abasis2,abasis2)
        SqrtDetA = SQRT(a11*a22)

        IF (Cylindrical .AND. BenchmarkProblem) THEN
          ! In the case of benchmark cases
          ! use a hard-coded load to avoid errors from representing the load: 
          !
          NormalTraction = h**3 * 1.0d5 * cos(2.0d0*y1)
          ApplyNormalPressure = .TRUE.
        END IF
      ELSE
        !
        ! Here we combine the Cartesian formulation and surface reparametrization
        !
        y1 = SUM( Nodes % x(1:nd) * Basis(1:nd) )
        y2 = SUM( Nodes % y(1:nd) * Basis(1:nd) )

        ! ------------------------------------------------------------------------------
        ! The fundamental forms at the point (y1,y2) of the principal coordinate patch.
        ! The Christoffel symbols Cijk are also returned, but here they will have no use.   
        ! ------------------------------------------------------------------------------
        CALL SurfaceBasisVectors(y1, y2, TaylorParams, e1, e2, e3, o, abasis1, &
            abasis2, abasis3, A11, A22, SqrtDetA, B11, B22, C111, C112, C221, C222, &
            C211, C212, PlanarPoint=PlateBody, Umbilical=SphericalSurface)

        K1 = B11/A11
        K2 = B22/A22
      END IF

      ! The geometric Lame parameters:
      ! ------------------------------
      A1 = SQRT(a11)
      A2 = SQRT(a22)

      Q(1,1) = DOT_PRODUCT(abasis1,i1)
      Q(1,2) = DOT_PRODUCT(abasis1,i2)
      Q(1,3) = DOT_PRODUCT(abasis1,i3)
      Q(2,1) = DOT_PRODUCT(abasis2,i1)
      Q(2,2) = DOT_PRODUCT(abasis2,i2)
      Q(2,3) = DOT_PRODUCT(abasis2,i3)
      Q(3,1) = DOT_PRODUCT(abasis3,i1)
      Q(3,2) = DOT_PRODUCT(abasis3,i2)
      Q(3,3) = DOT_PRODUCT(abasis3,i3)

      ! The matrix description of the elasticity tensor: Two versions
      ! differ on how the transverse normal stress part is formulated,
      ! but the result should be the same.
      !
      IF (GeneralMaterial) THEN
        CALL MaterialLaw3D(HMat, A1, A2, E, nu)
      ELSE
        CALL ElasticityMatrix(CMat, GMat, A1, A2, E, nu)
      END IF

      !-----------------------------------------------------------------------------------
      ! THE PART CORRESPONDING TO THE MEMBRANE STRAINS:
      !-----------------------------------------------------------------------------------
      ! Create first the representation of the differential DE_0(U)[V] of the linearized 
      ! membrane strain E_0(U) in the matrix form as DE_0(U)[V] = E_0(V) = BM * V (here 
      ! DE_0(U)[V] = E_0(V) holds for all U since E_0(U) is linear with respect to U).
      !------------------------------------------------------------------------------------
      Weight = h * SqrtDetA * detJ * sq
      DO p=1,nd
        BM(1,(p-1)*m+1:(p-1)*m+3) = Q(1,1:3) * dBasis(p,1)
        BM(2,(p-1)*m+1:(p-1)*m+3) = Q(2,1:3) * dBasis(p,2)
        BM(3,(p-1)*m+1:(p-1)*m+3) = Q(2,1:3) * dBasis(p,1) + Q(1,1:3) * dBasis(p,2)
        BM(4,(p-1)*m+4:(p-1)*m+6) = -Q(3,1:3) * Basis(p)
      END DO

      IF (.NOT. GeneralMaterial) THEN
        !
        ! Use a formulation which enforces vanishing normal stress by penalization;
        ! otherwise the basic 3D formulation is employed.
        !
        BM(4,:) = BM(4,:) + nu/((1.0d0-nu)*a11) * BM(1,:) + &
            nu/((1.0d0-nu)*a22) * BM(2,:)
      END IF

      StrainVec = 0.0d0
      NONLINEAR_MEMBRANE_STRAINS: IF (LargeDeflection) THEN
        ! ---------------------------------------------------------------------------------------
        ! The differential DE(U)[V] of the membrane strain E(U) is by definition linear with 
        ! respect to V and thus have a matrix representation DE(U)[V] ~ BM * V + NonlinBM(U) * V.
        ! The matrix BM is already created and here we create the matrix NonlinBM(U), which
        ! depends on the current solution iterate U.
        ! ---------------------------------------------------------------------------------------
        PrevGrad(1:3,1:2) = MATMUL(LocalSol(1:3,1:nd), dBasis(1:nd,1:2))
        PrevB(1:3) = MATMUL(LocalSol(4:6,1:nd), Basis(1:nd))

        DO p=1,nd
          NonlinBM(1,(p-1)*m+1:(p-1)*m+3) = PrevGrad(1:3,1) * dBasis(p,1)
          NonlinBM(2,(p-1)*m+1:(p-1)*m+3) = PrevGrad(1:3,2) * dBasis(p,2)
          NonlinBM(3,(p-1)*m+1:(p-1)*m+3) = PrevGrad(1:3,2) * dBasis(p,1) + PrevGrad(1:3,1) * dBasis(p,2)
          NonlinBM(4,(p-1)*m+4:(p-1)*m+6) = PrevB(1:3) * Basis(p)
        END DO

        IF (.NOT. GeneralMaterial) THEN
          NonlinBM(4,1:DOFs) = NonlinBM(4,1:DOFs) + nu/((1.0d0-nu)*a11) * NonlinBM(1,1:DOFs) + &
            nu/((1.0d0-nu)*a22) * NonlinBM(2,1:DOFs)
        END IF

        ! The nonlinear parts of membrane strain components:
        StrainVec(1:4) = 0.5_dp * MATMUL(NonlinBM(1:4,1:DOFs), PrevSolVec(1:DOFs))
      END IF NONLINEAR_MEMBRANE_STRAINS

      ! Add the linear part of strain for the current iterate:
      StrainVec(1:csize) = StrainVec(1:csize) + MATMUL( BM(1:csize,1:DOFs), PrevSolVec(1:DOFs) )
      
      IF (.NOT. GeneralMaterial) THEN
        CALL StrainEnergyDensity(Stiff, CMat, BM + NonlinBM, csize, DOFs, Weight)
        StressVec(1:csize) = MATMUL(CMat(1:csize,1:csize), StrainVec(1:csize))
      ELSE
        CALL StrainEnergyDensity(Stiff, HMat, BM + NonlinBM, csize, DOFs, Weight)
        StressVec(1:csize) = MATMUL(HMat(1:csize,1:csize), StrainVec(1:csize))        
      END IF

      ! Residual terms for RHS:
      Force(1:DOFs) = Force(1:DOFs) - MATMUL( TRANSPOSE(BM(1:csize,1:DOFs) + NonlinBM(1:csize,1:DOFs)), &
         StressVec(1:csize) ) * Weight

      ! The remaining terms for the complete Newton iteration:
      !
      NONLINEAR_MEMBRANE_EFFECTS: IF (LargeDeflection) THEN
        IF (.NOT. GeneralMaterial) THEN
          c1 = (StressVec(1) + nu/((1.0d0-nu)*a11) * StressVec(4)) * Weight
          c2 = (StressVec(2) + nu/((1.0d0-nu)*a22) * StressVec(4)) * Weight
        ELSE
          c1 = StressVec(1) * Weight
          c2 = StressVec(2) * Weight
        END IF
        c3 = StressVec(3) * Weight
        c4 = StressVec(4) * Weight

        DO p=1,nd
          Stiff((p-1)*m+1,1:DOFs:m) = Stiff((p-1)*m+1,1:DOFs:m) + c1 * dbasis(p,1) * dbasis(1:nd,1) + &
              c2 * dbasis(p,2) * dbasis(1:nd,2) + c3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dbasis(1:nd,1))
          Stiff((p-1)*m+2,2:DOFs:m) = Stiff((p-1)*m+2,2:DOFs:m) + c1 * dbasis(p,1) * dbasis(1:nd,1) + &
              c2 * dbasis(p,2) * dbasis(1:nd,2) + c3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dbasis(1:nd,1))
          Stiff((p-1)*m+3,3:DOFs:m) = Stiff((p-1)*m+3,3:DOFs:m) + c1 * dbasis(p,1) * dbasis(1:nd,1) + &
              c2 * dbasis(p,2) * dbasis(1:nd,2) + c3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dBasis(1:nd,1))
        END DO

        DO p=1,nd
          Stiff((p-1)*m+4,4:DOFs:m) = Stiff((p-1)*m+4,4:DOFs:m) + c4 * Basis(p) * Basis(1:nd)
          Stiff((p-1)*m+5,5:DOFs:m) = Stiff((p-1)*m+5,5:DOFs:m) + c4 * Basis(p) * Basis(1:nd)
          Stiff((p-1)*m+6,6:DOFs:m) = Stiff((p-1)*m+6,6:DOFs:m) + c4 * Basis(p) * Basis(1:nd)
        END DO
      END IF NONLINEAR_MEMBRANE_EFFECTS

      !-----------------------------------------------------------------------------------
      ! THE PART CORRESPONDING TO THE TRANSVERSE SHEAR STRAINS:
      !-----------------------------------------------------------------------------------
      ! Create first the representation of the differential DE_0(U)[V] of the linearized 
      ! transverse shear strain E_0(U) in the matrix form as DE_0(U)[V] = E_0(V) = BS * V
      ! (here DE_0(U)[V] = E_0(V) holds since E_0(U) is linear with respect to U).
      !------------------------------------------------------------------------------------
      DO p=1,nd
        BS(1:2,(p-1)*m+4) = -Q(1:2,1) * Basis(p)
        BS(1:2,(p-1)*m+5) = -Q(1:2,2) * Basis(p)
        BS(1:2,(p-1)*m+6) = -Q(1:2,3) * Basis(p)

        BS(1:2,(p-1)*m+1) = Q(3,1) * dBasis(p,1:2)
        BS(1:2,(p-1)*m+2) = Q(3,2) * dBasis(p,1:2)
        BS(1:2,(p-1)*m+3) = Q(3,3) * dBasis(p,1:2)
      END DO

      NONLINEAR_SHEAR_STRAINS: IF (LargeDeflection) THEN
        ! ---------------------------------------------------------------------------------------
        ! The representation of the differential DE(U)[V] of the transverse shear strain E(U) in 
        ! the matrix form as DE(U)[V] ~ BS * V + NonlinBS(U) * V. The matrix BS is already created 
        ! and here we compute the matrix NonlinBS(U), which depends on the current solution iterate 
        ! U.
        ! ---------------------------------------------------------------------------------------
        DO p=1,nd
          NonlinBS(1,(p-1)*m+1:(p-1)*m+3) = -PrevB(1:3) * dBasis(p,1)
          NonlinBS(1,(p-1)*m+4:(p-1)*m+6) = -PrevGrad(1:3,1) * Basis(p)
          NonlinBS(2,(p-1)*m+1:(p-1)*m+3) = -PrevB(1:3) * dBasis(p,2)
          NonlinBS(2,(p-1)*m+4:(p-1)*m+6) = -PrevGrad(1:3,2) * Basis(p)
        END DO
        StrainVec(5:6) = 0.5_dp * MATMUL(NonlinBS(1:2,1:DOFs), PrevSolVec(1:DOFs))
      END IF NONLINEAR_SHEAR_STRAINS

      ! Add the linear part of shear strain for the current iterate:
      StrainVec(5:6) = StrainVec(5:6) + MATMUL(BS(1:2,1:DOFs), PrevSolVec(1:DOFs))

      IF (.NOT. GeneralMaterial) THEN
        CALL StrainEnergyDensity(Stiff, GMat, BS+NonlinBS, 2, DOFs, Weight)
        StressVec(5:6) = MATMUL(GMat(1:2,1:2), StrainVec(5:6))
      ELSE
        CALL StrainEnergyDensity(Stiff, HMat(5:6,5:6), BS + NonlinBS, 2, DOFs, Weight)
        StressVec(5:6) = MATMUL(HMat(5:6,5:6), StrainVec(5:6))       
      END IF

      ! Residual terms for RHS:
      Force(1:DOFs) = Force(1:DOFs) - MATMUL( TRANSPOSE(BS(1:2,1:DOFs) + NonlinBS(1:2,1:DOFs)), &
          StressVec(5:6) ) * Weight

      ! The remaining terms for the complete Newton iteration:
      !
      NONLINEAR_SHEAR_EFFECTS: IF (LargeDeflection) THEN
        c1 = StressVec(5) * Weight
        c2 = StressVec(6) * Weight
        DO p=1,nd
          Stiff((p-1)*m+1,4:DOFs:m) = Stiff((p-1)*m+1,4:DOFs:m) - c1 * dbasis(p,1) * Basis(1:nd) - &
              c2 * dbasis(p,2) * Basis(1:nd)
          Stiff((p-1)*m+2,5:DOFs:m) = Stiff((p-1)*m+2,5:DOFs:m) - c1 * dbasis(p,1) * Basis(1:nd) - &
              c2 * dbasis(p,2) * Basis(1:nd)
          Stiff((p-1)*m+3,6:DOFs:m) = Stiff((p-1)*m+3,6:DOFs:m) - c1 * dbasis(p,1) * Basis(1:nd) - &
              c2 * dbasis(p,2) * Basis(1:nd)

          Stiff((p-1)*m+4,1:DOFs:m) = Stiff((p-1)*m+4,1:DOFs:m) - c1 * Basis(p) * dbasis(1:nd,1) - &
              c2 * Basis(p) * dbasis(1:nd,2)
          Stiff((p-1)*m+5,2:DOFs:m) = Stiff((p-1)*m+5,2:DOFs:m) - c1 * Basis(p) * dbasis(1:nd,1) - &
              c2 * Basis(p) * dbasis(1:nd,2)
          Stiff((p-1)*m+6,3:DOFs:m) = Stiff((p-1)*m+6,3:DOFs:m) - c1 * Basis(p) * dbasis(1:nd,1) - &
              c2 * Basis(p) * dbasis(1:nd,2)
        END DO
      END IF NONLINEAR_SHEAR_EFFECTS

      !----------------------------------------------------------------------------------------
      ! The part of transverse shear strains which depend linearly on the thickness coordinate: 
      ! It appears that considering the linear part may not be meaningful without having a cubic
      ! displacement approximation in the thickness coordinate. On the other hand, this seems
      ! to be have a strong impact on the solution of the 6-field model when nonlinear effects
      ! are strong.
      !----------------------------------------------------------------------------------------
      Weight = h**2/12.0d0 * Weight
      DO p=1,nd
        BS(3:4,(p-1)*m+4) = Q(3,1) * dBasis(p,1:2)
        BS(3:4,(p-1)*m+5) = Q(3,2) * dBasis(p,1:2)
        BS(3:4,(p-1)*m+6) = Q(3,3) * dBasis(p,1:2)

        BS(3,(p-1)*m+1:(p-1)*m+3) = -K1 * Q(3,1:3) * dBasis(p,1)
        BS(4,(p-1)*m+1:(p-1)*m+3) = -K2 * Q(3,1:3) * dBasis(p,2)

        IF (TransverseBendingStretch) THEN
          BS(3,(p-1)*m+7:(p-1)*m+9) = Q(1,1:3) * Basis(p)
          BS(4,(p-1)*m+7:(p-1)*m+9) = Q(2,1:3) * Basis(p)
        END IF
      END DO

      StrainVec(1:2) = 0.0d0
      NONLINEAR_SHEAR_OVER_THICKNESS: IF (LargeDeflection .AND. NonlinearBending) THEN
        PrevGradB(1:3,1:2) = MATMUL(LocalSol(4:6,1:nd), dBasis(1:nd,1:2))
        DO p=1,nd
          NonlinBS(3,(p-1)*m+4:(p-1)*m+6) = -PrevB(1:3) * dBasis(p,1) - PrevGradB(1:3,1) * Basis(p)
          NonlinBS(4,(p-1)*m+4:(p-1)*m+6) = -PrevB(1:3) * dBasis(p,2) - PrevGradB(1:3,2) * Basis(p)        
        END DO
        NonlinBS(3,:) = NonlinBS(3,:) - K1 * NonlinBS(1,:)
        NonlinBS(4,:) = NonlinBS(4,:) - K2 * NonlinBS(2,:)
        IF (TransverseBendingStretch) THEN
          PrevPsi(1:3) = MATMUL(LocalSol(7:9,1:nd), Basis(1:nd))
          DO p=1,nd
            NonlinBS(3,(p-1)*m+1:(p-1)*m+3) = PrevPsi(1:3) * dBasis(p,1)
            NonlinBS(3,(p-1)*m+7:(p-1)*m+9) = PrevGrad(1:3,1) * Basis(p)
            NonlinBS(4,(p-1)*m+1:(p-1)*m+3) = PrevPsi(1:3) * dBasis(p,2)
            NonlinBS(4,(p-1)*m+7:(p-1)*m+9) = PrevGrad(1:3,2) * Basis(p)
          END DO
        END IF
        StrainVec(1:2) = 0.5_dp * MATMUL(NonlinBS(3:4,1:DOFs), PrevSolVec(1:DOFs))
      END IF NONLINEAR_SHEAR_OVER_THICKNESS

      StrainVec(1:2) = StrainVec(1:2) + MATMUL(BS(3:4,1:DOFs), PrevSolVec(1:DOFs))

      IF (.NOT. GeneralMaterial) THEN
        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + Weight * &
            MATMUL(TRANSPOSE(BS(3:4,1:DOFs) + NonlinBS(3:4,1:DOFs)), &
            MATMUL(GMat(1:2,1:2), BS(3:4,1:DOFs) + NonlinBS(3:4,1:DOFs)))
        StressVec(1:2) = MATMUL(GMat(1:2,1:2), StrainVec(1:2))
      ELSE
        Stiff(1:DOFs,1:DOFs) = Stiff(1:DOFs,1:DOFs) + Weight * &
            MATMUL(TRANSPOSE(BS(3:4,1:DOFs) + NonlinBS(3:4,1:DOFs)), &
            MATMUL(HMat(5:6,5:6),BS(3:4,1:DOFs) + NonlinBS(3:4,1:DOFs)))
        StressVec(1:2) = MATMUL(HMat(5:6,5:6), StrainVec(1:2))
      END IF
      !
      ! Residual terms for RHS:
      Force(1:DOFs) = Force(1:DOFs) - MATMUL(TRANSPOSE(BS(3:4,1:DOFs) + NonlinBS(3:4,1:DOFs)), &
          StressVec(1:2) ) * Weight

      ! The remaining terms for the complete Newton iteration:
      !
      SHEAR_EFFECTS_OVER_THICKNESS:  IF (LargeDeflection .AND. NonlinearBending) THEN
        v1 = -K1 * StressVec(1) * Weight
        v2 = -K2 * StressVec(2) * Weight
        c1 = StressVec(1) * Weight
        c2 = StressVec(2) * Weight
        DO p=1,nd
          Stiff((p-1)*m+4,4:DOFs:m) = Stiff((p-1)*m+4,4:DOFs:m) - c1 * ( basis(p) * dbasis(1:nd,1) + &
              dbasis(p,1) * basis(1:nd) ) - c2 * ( basis(p) * dbasis(1:nd,2) + &
              dbasis(p,2) * basis(1:nd) )
          Stiff((p-1)*m+5,5:DOFs:m) = Stiff((p-1)*m+5,5:DOFs:m) - c1 * ( basis(p) * dbasis(1:nd,1) + &
              dbasis(p,1) * basis(1:nd) ) - c2 * ( basis(p) * dbasis(1:nd,2) + &
              dbasis(p,2) * basis(1:nd) )
          Stiff((p-1)*m+6,6:DOFs:m) = Stiff((p-1)*m+6,6:DOFs:m) - c1 * ( basis(p) * dbasis(1:nd,1) + &
              dbasis(p,1) * basis(1:nd) ) - c2 * ( basis(p) * dbasis(1:nd,2) + &
              dbasis(p,2) * basis(1:nd) )

          Stiff((p-1)*m+1,4:DOFs:m) = Stiff((p-1)*m+1,4:DOFs:m) - v1 * dbasis(p,1) * Basis(1:nd) - &
              v2 * dbasis(p,2) * Basis(1:nd)
          Stiff((p-1)*m+2,5:DOFs:m) = Stiff((p-1)*m+2,5:DOFs:m) - v1 * dbasis(p,1) * Basis(1:nd) - &
              v2 * dbasis(p,2) * Basis(1:nd)
          Stiff((p-1)*m+3,6:DOFs:m) = Stiff((p-1)*m+3,6:DOFs:m) - v1 * dbasis(p,1) * Basis(1:nd) - &
              v2 * dbasis(p,2) * Basis(1:nd)

          Stiff((p-1)*m+4,1:DOFs:m) = Stiff((p-1)*m+4,1:DOFs:m) - v1 * Basis(p) * dbasis(1:nd,1) - &
              v2 * Basis(p) * dbasis(1:nd,2)
          Stiff((p-1)*m+5,2:DOFs:m) = Stiff((p-1)*m+5,2:DOFs:m) - v1 * Basis(p) * dbasis(1:nd,1) - &
              v2 * Basis(p) * dbasis(1:nd,2)
          Stiff((p-1)*m+6,3:DOFs:m) = Stiff((p-1)*m+6,3:DOFs:m) - v1 * Basis(p) * dbasis(1:nd,1) - &
              v2 * Basis(p) * dbasis(1:nd,2)
        END DO

        IF (TransverseBendingStretch) THEN
          c1 = StressVec(1) * Weight
          c2 = StressVec(2) * Weight
          DO p=1,nd
            Stiff((p-1)*m+7,1:DOFs:m) = Stiff((p-1)*m+7,1:DOFs:m) + c1 * Basis(p) * dBasis(1:nd,1) + &
                c2 * Basis(p) * dBasis(1:nd,2)
            Stiff((p-1)*m+8,2:DOFs:m) = Stiff((p-1)*m+8,2:DOFs:m) + c1 * Basis(p) * dBasis(1:nd,1) + &
                c2 * Basis(p) * dBasis(1:nd,2)
            Stiff((p-1)*m+9,3:DOFs:m) = Stiff((p-1)*m+9,3:DOFs:m) + c1 * Basis(p) * dBasis(1:nd,1) + &
                c2 * Basis(p) * dBasis(1:nd,2)

            Stiff((p-1)*m+1,7:DOFs:m) = Stiff((p-1)*m+1,7:DOFs:m) + c1 * Basis(1:nd) * dBasis(p,1) + &
                c2 * Basis(1:nd) * dBasis(p,2)
            Stiff((p-1)*m+2,8:DOFs:m) = Stiff((p-1)*m+2,8:DOFs:m) + c1 * Basis(1:nd) * dBasis(p,1) + &
                c2 * Basis(1:nd) * dBasis(p,2)
            Stiff((p-1)*m+3,9:DOFs:m) = Stiff((p-1)*m+3,9:DOFs:m) + c1 * Basis(1:nd) * dBasis(p,1) + &
                c2 * Basis(1:nd) * dBasis(p,2)
          END DO
        END IF
      END IF SHEAR_EFFECTS_OVER_THICKNESS

      !---------------------------------------------------------------
      ! THE PART CORRESPONDING TO THE BENDING STRAINS:
      !---------------------------------------------------------------      
      Weight = h**3/12.0d0 * SqrtDetA * detJ * sq
      DO p=1,nd
        BB(1,(p-1)*m+4:(p-1)*m+6) = Q(1,1:3) * dBasis(p,1)
        BB(2,(p-1)*m+4:(p-1)*m+6) = Q(2,1:3) * dBasis(p,2)
        BB(3,(p-1)*m+4:(p-1)*m+6) = Q(2,1:3) * dBasis(p,1) + Q(1,1:3) * dBasis(p,2)
        BB(3,(p-1)*m+1:(p-1)*m+3) = -K1 * Q(2,1:3) * dBasis(p,1) - &
            K2 * Q(1,1:3) * dBasis(p,2)
      END DO
      BB(1,:) = BB(1,:) - K1 * BM(1,:)
      BB(2,:) = BB(2,:) - K2 * BM(2,:)

      IF (TransverseBendingStretch) THEN
        DO p=1,nd
          BB(4,(p-1)*m+7:(p-1)*m+9) = Q(3,1:3) * Basis(p)
        END DO
        IF (.NOT. GeneralMaterial) THEN
          BB(4,:) = BB(4,:) + nu/((1.0d0-nu)*a11) * BB(1,:) + &
              nu/((1.0d0-nu)*a22) * BB(2,:)
        END IF
      END IF

      StrainVec(1:bsize) = 0.0d0
      NONLINEAR_BENDING_STRAINS: IF (LargeDeflection .AND. NonlinearBending) THEN
        PrevGradB(1:3,1:2) = MATMUL(LocalSol(4:6,1:nd), dBasis(1:nd,1:2))
        DO p=1,nd
          NonlinBB(1,(p-1)*m+1:(p-1)*m+3) = PrevGradB(1:3,1) * dBasis(p,1)
          NonlinBB(1,(p-1)*m+4:(p-1)*m+6) = PrevGrad(1:3,1) * dBasis(p,1)
          NonlinBB(2,(p-1)*m+1:(p-1)*m+3) = PrevGradB(1:3,2) * dBasis(p,2)
          NonlinBB(2,(p-1)*m+4:(p-1)*m+6) = PrevGrad(1:3,2) * dBasis(p,2)
          NonlinBB(3,(p-1)*m+1:(p-1)*m+3) = PrevGradB(1:3,2) * dBasis(p,1) + PrevGradB(1:3,1) * dBasis(p,2)
          NonlinBB(3,(p-1)*m+4:(p-1)*m+6) = PrevGrad(1:3,2) * dBasis(p,1) + PrevGrad(1:3,1) * dBasis(p,2)
        END DO
        NonlinBB(1:3,:) = NonlinBB(1:3,:) - (K1 + K2) * NonlinBM(1:3,:)

        IF (TransverseBendingStretch) THEN
          PrevPsi(1:3) = MATMUL(LocalSol(7:9,1:nd), Basis(1:nd))
          DO p=1,nd
            NonlinBB(4,(p-1)*m+7:(p-1)*m+9) = -PrevB(1:3) * Basis(p)
            NonlinBB(4,(p-1)*m+4:(p-1)*m+6) = -PrevPsi(1:3) * Basis(p)
          END DO
          IF (.NOT. GeneralMaterial) THEN
            NonlinBB(4,1:DOFs) = NonlinBB(4,1:DOFs) + nu/((1.0d0-nu)*a11) * NonlinBB(1,1:DOFs) + &
                nu/((1.0d0-nu)*a22) * NonlinBB(2,1:DOFs)
          END IF
        END IF
        StrainVec(1:bsize) = 0.5_dp * MATMUL(NonlinBB(1:bsize,1:DOFs), PrevSolVec(1:DOFs))
      END IF NONLINEAR_BENDING_STRAINS

      StrainVec(1:bsize) = StrainVec(1:bsize) + MATMUL(BB(1:bsize,1:DOFs), PrevSolVec(1:DOFs))

      IF (.NOT. GeneralMaterial) THEN
        CALL StrainEnergyDensity(Stiff, CMat, BB + NonlinBB, bsize, DOFs, Weight)
        StressVec(1:bsize) = MATMUL(CMat(1:bsize,1:bsize), StrainVec(1:bsize))
      ELSE
        CALL StrainEnergyDensity(Stiff, HMat, BB + NonlinBB, bsize, DOFs, Weight)
        StressVec(1:bsize) = MATMUL(HMat(1:bsize,1:bsize), StrainVec(1:bsize))
      END IF

      ! Residual terms for RHS:
      Force(1:DOFs) = Force(1:DOFs) - MATMUL( TRANSPOSE( BB(1:bsize,1:DOFs) + NonlinBB(1:bsize,1:DOFs)), &
          StressVec(1:bsize) ) * Weight

      ! The remaining terms for the complete Newton iteration:
      !
      NONLINEAR_BENDING_EFFECTS: IF (LargeDeflection .AND. NonlinearBending) THEN
        IF (.NOT. GeneralMaterial) THEN
          c1 = (StressVec(1) + nu/((1.0d0-nu)*a11) * StressVec(4)) * Weight
          c2 = (StressVec(2) + nu/((1.0d0-nu)*a22) * StressVec(4)) * Weight
          v1 = -(K1 + K2) * (StressVec(1) + nu/((1.0d0-nu)*a11) * StressVec(4)) * Weight
          v2 = -(K1 + K2) * (StressVec(2) + nu/((1.0d0-nu)*a22) * StressVec(4)) * Weight
        ELSE
          c1 = StressVec(1) * Weight
          c2 = StressVec(2) * Weight
          v1 = -(K1 + K2) * StressVec(1) * Weight
          v2 = -(K1 + K2) * StressVec(2) * Weight
        END IF
        c3 = StressVec(3) * Weight
        v3 = -(K1 + K2) * StressVec(3) * Weight
        c4 = StressVec(4) * Weight

        DO p=1,nd
          Stiff((p-1)*m+1,4:DOFs:m) = Stiff((p-1)*m+1,4:DOFs:m) + c1 * dbasis(p,1) * dbasis(1:nd,1) + &
              c2 * dbasis(p,2) * dbasis(1:nd,2) + c3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dbasis(1:nd,1))
          Stiff((p-1)*m+2,5:DOFs:m) = Stiff((p-1)*m+2,5:DOFs:m) + c1 * dbasis(p,1) * dbasis(1:nd,1) + &
              c2 * dbasis(p,2) * dbasis(1:nd,2) + c3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dbasis(1:nd,1))
          Stiff((p-1)*m+3,6:DOFs:m) = Stiff((p-1)*m+3,6:DOFs:m) + c1 * dbasis(p,1) * dbasis(1:nd,1) + &
              c2 * dbasis(p,2) * dbasis(1:nd,2) + c3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dBasis(1:nd,1))

          Stiff((p-1)*m+4,1:DOFs:m) = Stiff((p-1)*m+4,1:DOFs:m) + c1 * dbasis(p,1) * dbasis(1:nd,1) + &
              c2 * dbasis(p,2) * dbasis(1:nd,2) + c3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dbasis(1:nd,1))
          Stiff((p-1)*m+5,2:DOFs:m) = Stiff((p-1)*m+5,2:DOFs:m) + c1 * dbasis(p,1) * dbasis(1:nd,1) + &
              c2 * dbasis(p,2) * dbasis(1:nd,2) + c3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dbasis(1:nd,1))
          Stiff((p-1)*m+6,3:DOFs:m) = Stiff((p-1)*m+6,3:DOFs:m) + c1 * dbasis(p,1) * dbasis(1:nd,1) + &
              c2 * dbasis(p,2) * dbasis(1:nd,2) + c3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dBasis(1:nd,1))

          Stiff((p-1)*m+1,1:DOFs:m) = Stiff((p-1)*m+1,1:DOFs:m) + v1 * dbasis(p,1) * dbasis(1:nd,1) + &
              v2 * dbasis(p,2) * dbasis(1:nd,2) + v3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dbasis(1:nd,1))
          Stiff((p-1)*m+2,2:DOFs:m) = Stiff((p-1)*m+2,2:DOFs:m) + v1 * dbasis(p,1) * dbasis(1:nd,1) + &
              v2 * dbasis(p,2) * dbasis(1:nd,2) + v3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dbasis(1:nd,1))
          Stiff((p-1)*m+3,3:DOFs:m) = Stiff((p-1)*m+3,3:DOFs:m) + v1 * dbasis(p,1) * dbasis(1:nd,1) + &
              v2 * dbasis(p,2) * dbasis(1:nd,2) + v3 * (dbasis(p,1) * dbasis(1:nd,2) + &
              dbasis(p,2) * dBasis(1:nd,1))
        END DO

        IF (TransverseBendingStretch) THEN
          DO p=1,nd
            Stiff((p-1)*m+4,7:DOFs:m) = Stiff((p-1)*m+4,7:DOFs:m) - c4 * Basis(p) * Basis(1:nd)
            Stiff((p-1)*m+5,8:DOFs:m) = Stiff((p-1)*m+5,8:DOFs:m) - c4 * Basis(p) * Basis(1:nd)
            Stiff((p-1)*m+6,9:DOFs:m) = Stiff((p-1)*m+6,9:DOFs:m) - c4 * Basis(p) * Basis(1:nd)

            Stiff((p-1)*m+7,4:DOFs:m) = Stiff((p-1)*m+7,4:DOFs:m) - c4 * Basis(p) * Basis(1:nd)
            Stiff((p-1)*m+8,5:DOFs:m) = Stiff((p-1)*m+8,5:DOFs:m) - c4 * Basis(p) * Basis(1:nd)
            Stiff((p-1)*m+9,6:DOFs:m) = Stiff((p-1)*m+9,6:DOFs:m) - c4 * Basis(p) * Basis(1:nd)
          END DO
        END IF
      END IF NONLINEAR_BENDING_EFFECTS

      !----------------------------------------------------------------
      ! Mass matrix without bubbles taken into account:
      !----------------------------------------------------------------     
      IF ( MassAssembly ) THEN
        DO k=1,3
          SELECT CASE(k)
          CASE(1)
            Weight = 1/a11 * h * rho0 * SqrtDetA * detJ * sq
          CASE(2)
            Weight = 1/a22 * h * rho0 * SqrtDetA * detJ * sq
          CASE(3)
            Weight = h * rho0 * SqrtDetA * detJ * sq
          END SELECT
          DO i=1,nd
            DO j=1,nd
              Mass((i-1)*m+k,(j-1)*m+k) = Mass((i-1)*m+k,(j-1)*m+k) + &
                  Basis(i) * Basis(j) * Weight
              Damp((i-1)*m+k,(j-1)*m+k) = Damp((i-1)*m+k,(j-1)*m+k) + &
                  DampCoef * Basis(i) * Basis(j) * Weight              

              Mass((i-1)*m+3+k,(j-1)*m+3+k) = Mass((i-1)*m+3+k,(j-1)*m+3+k) + &
                  h**2/12.0d0 * Basis(i) * Basis(j) * Weight
            END DO
          END DO
        END DO
      END IF

      !----------------------------------------------------------------
      ! RHS vector:
      !----------------------------------------------------------------
      IF (ApplyNormalPressure) THEN
        IF (LargeDeflection) THEN
          !----------------------------------------------------------------
          ! Compute the normal vector n to the deformed mid-surface using
          ! the current iterate and apply the normal traction p * n, with
          ! the effect of area change being taken into account. 
          !----------------------------------------------------------------
          v1 = SUM(Basis(1:nd) * PrevSolVec(1:DOFs:m))
          v2 = SUM(Basis(1:nd) * PrevSolVec(2:DOFs:m))
          v3 = SUM(Basis(1:nd) * PrevSolVec(3:DOFs:m))
          abasis1New(1:3) = abasis1(1:3) + PrevGrad(1:3,1)
          abasis2New(1:3) = abasis2(1:3) + PrevGrad(1:3,2)
          NewDetA = DOT_PRODUCT(abasis1New,abasis1New) * DOT_PRODUCT(abasis2New,abasis2New) - &
              DOT_PRODUCT(abasis1New,abasis2New)**2
          abasis3New(1:3) = CrossProduct(abasis1New,abasis2New)
          Norm = SQRT(SUM(abasis3New(:)**2))
          abasis3New(1:3) = abasis3New(1:3)/Norm

          Weight = SQRT(NewDetA) * detJ * sq

          DO i=1,3
            RHSForce(i:DOFs:m) = RHSForce(i:DOFs:m) + NormalTraction * abasis3New(i) * Basis(1:nd) * Weight
            Force(i:DOFs:m) = Force(i:DOFs:m) + NormalTraction * abasis3New(i) * Basis(1:nd) * Weight
            ! TO DO: Add terms related to the first-order terms in the normal coordinate
          END DO
        ELSE
          Weight = SqrtDetA * detJ * sq
          DO p=1,nd
            i = m*(p-1)
            Force(i+1:i+3) = Force(i+1:i+3) + NormalTraction * Q(3,1:3) * Basis(p) * Weight
          END DO
        END IF
      END IF
    END DO QUADRATURELOOP

    IF ( MassAssembly ) THEN
      IF ( TransientSimulation ) THEN
        CALL Default2ndOrderTime(MASS,DAMP,STIFF,FORCE)
      ELSE IF ( HarmonicAssembly ) THEN
        CALL DefaultUpdateMass( MASS )
        CALL DefaultUpdateDamp( DAMP )
      END IF
    END IF

    CALL DefaultUpdateEquations(STIFF,FORCE)

!------------------------------------------------------------------------------
  END SUBROUTINE ShellLocalMatrixCartesian
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Create the element mapping corresponding to the principal curvature coordinates
! in terms of the basis which is used in the shell analysis.
!------------------------------------------------------------------------------
  SUBROUTINE SolveNodesVariables(Element, Nodes, nd, GElement, PatchNodes)
! -----------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element   ! The element type for which nodes are written
    TYPE(Nodes_t), INTENT(INOUT) :: Nodes             ! The nodes data structure to be updated
    INTEGER, INTENT(IN) :: nd                         ! The number of coordinate entries written
    TYPE(Element_t), POINTER, INTENT(IN) :: GElement  ! A Lagrange element used in surface reconstruction
    REAL(KIND=dp), INTENT(IN) :: PatchNodes(:,:)      ! The nodes data compatible with GElement
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: Stat, PVersion
    INTEGER :: j, k, n, t, GElementNodes
    REAL(KIND=dp) :: u, v, yk, up, vp, DetJ
    REAL(KIND=dp) :: GBasis(GElement % Type % NumberOfNodes)
    REAL(KIND=dp) :: PBasis(nd)
    REAL(KIND=dp) :: Stiff(nd, nd)
    REAL(KIND=dp) :: Force(nd)
!------------------------------------------------------------------------------
    GElementNodes = GElement % Type % NumberOfNodes

    Nodes % x(:) = 0.0d0
    Nodes % y(:) = 0.0d0
    Nodes % z(:) = 0.0d0

    IF (GElement % Type % NumberOfNodes == Element % Type % NumberOfNodes) THEN
      Nodes % x(1:GElementNodes) = PatchNodes(1:GElementNodes,1)
      Nodes % y(1:GElementNodes) = PatchNodes(1:GElementNodes,2)
      RETURN
    END IF

    PVersion = IsActivePElement(Element)
    n = GetElementFamily(Element)
    ! -----------------------------------------------
    ! The corner values can be written immediately:
    ! -----------------------------------------------
    Nodes % x(1:n) = PatchNodes(1:n,1)
    Nodes % y(1:n) = PatchNodes(1:n,2)

    IF (nd > n) THEN
      ! ------------------------------------------------------
      ! Solve a local matrix equation in the case of higher p
      ! TO DO: 1) Solve for entries i > n only, and possibly
      !        2) truncate calculation in the case of p > 3
      ! ------------------------------------------------------
      IP = GaussPoints(Element)
      DO k=1,2
        STIFF = 0.0d0
        FORCE = 0.0d0

        DO t=1,IP % n
          up = IP % U(t)
          vp = IP % V(t)
          ! ------------------------------------------------------------------
          ! First, evaluate the basis functions used in the shell analysis.
          ! Note that here Nodes and detJ are needed as dummies, but the basis 
          ! evaluation does not depend on their values.
          ! ------------------------------------------------------------------
          stat = ElementInfo(Element, Nodes, up, vp, 0.0d0, detJ, PBasis)
          ! ------------------------------------------------------------------
          ! Map the p-point to the point on the standard reference element:
          ! ------------------------------------------------------------------
          IF (Family==3 .AND. PVersion) THEN
            u = 0.5d0 * (1.0d0 + up - 1.0d0/sqrt(3.0d0) * vp)
            v = 1.0d0/sqrt(3.0d0) * vp
          ELSE
            u = up
            v = vp
          END IF

          CALL NodalBasisFunctions2D(GBasis, GElement, u, v)

          DO i=1,nd
            DO j=1,nd
              STIFF(i,j) = STIFF(i,j) + IP % s(t) * PBasis(i) * PBasis(j)
            END DO
          END DO

          SELECT CASE(k)
          CASE(1)
            yk = SUM(PatchNodes(1:GElementNodes,1) * GBasis(1:GElementNodes))
          CASE(2)
            yk = SUM(PatchNodes(1:GElementNodes,2) * GBasis(1:GElementNodes))
          END SELECT
          FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * yk * PBasis(1:nd)
        END DO

        CALL LUSolve(nd, Stiff, Force)
        SELECT CASE(k)
        CASE(1)
          Nodes % x(1:nd) = Force(1:nd)
        CASE(2)
          Nodes % y(1:nd) = Force(1:nd) 
        END SELECT
      END DO
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE SolveNodesVariables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE ShellSolver
!------------------------------------------------------------------------------
