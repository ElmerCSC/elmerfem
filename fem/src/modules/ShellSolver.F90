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
! *  The nodal director data should be available via reading from file mesh.director 
! *  located in the same place as the standard mesh files or, alternatively, the 
! *  user may provide mesh.elements.data file which should define the nodal 
! *  director field associated with the name 'director'.
! *
! *  This solver is STILL UNDER DEVELOPMENT and some possibilities of the strategy
! *  are not yet fully utilized. Note the current restrictions:
! *        -- Strain reduction operators have been worked out for 
! *           the lowest-order finite elements only.
! *        -- Parallel file formats for the director data are missing
! *        -- Postprocessing routines are also missing 
! *        -- Some terms related to normal stretch parameter are ignored in 
! *           the expression for the strain energy 
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
SUBROUTINE ShellSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverPars
  LOGICAL :: SavePrincipalAxes
  INTEGER  :: i
!------------------------------------------------------------------------------
  SolverPars => GetSolverParams()

  CALL ListAddInteger(SolverPars, 'Variable DOFs', 6)
  CALL ListAddLogical(SolverPars, 'Bubbles in Global System', .TRUE.)

  !----------------------------------------------------------------------------
  ! Create variables for saving principal (curvature) directions:
  !----------------------------------------------------------------------------
  SavePrincipalAxes = GetLogical(SolverPars, 'Principal Axes Output')

  IF (SavePrincipalAxes) THEN
    i=1
    DO WHILE(.TRUE.)
      IF ( .NOT.ListCheckPresent(SolverPars, &
          "Exported Variable "//TRIM(i2s(i))) ) EXIT
      i = i + 1
    END DO
    CALL ListAddString( SolverPars, "Exported Variable "//TRIM(i2s(i)), &
        "Principal Coordinate Dir1[Principal Coordinate Dir1:3]" )  
    i = i + 1
    CALL ListAddString( SolverPars, "Exported Variable "//TRIM(i2s(i)), &
        "Principal Coordinate Dir2[Principal Coordinate Dir2:3]" )
    i = i + 1
    CALL ListAddString( SolverPars, "Exported Variable "//TRIM(i2s(i)), &
        "Principal Coordinate Dir3[Principal Coordinate Dir3:3]" )  
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE ShellSolver_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE ShellSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE ElementDescription

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
  INTEGER, PARAMETER :: MaxNumberOCurves = 6
  INTEGER, PARAMETER :: CurveDataSize1 = 5
  INTEGER, PARAMETER :: CurveDataSize2 = 7
  INTEGER, DIMENSION(CurveDataSize1), PARAMETER :: CurveParams1 = (/ 1, 2, 3, 4, 5/)
  INTEGER, DIMENSION(CurveDataSize2), PARAMETER :: CurveParams2 = (/ 1, 2, 3, 4, 5, 6, 7/) 
!------------------------------------------------------------------------------
! Some other parameter definitions related to the choice of strain reduction 
! method, geometric tolerances ... 
!------------------------------------------------------------------------------
  INTEGER, PARAMETER :: AutomatedChoice = -1
  INTEGER, PARAMETER :: NoStrainReduction = 0
  INTEGER, PARAMETER :: CurlKernel = 1
  INTEGER, PARAMETER :: MITC = 2
  INTEGER, PARAMETER :: DoubleReduction = 3 

  INTEGER, PARAMETER :: MaxBGElementNodes = 9
  INTEGER, PARAMETER :: MaxPatchNodes = 16    ! The maximum node count for the surface description 

  INTEGER, PARAMETER :: GeometryMaxIters = 50
  REAL(KIND=dp), PARAMETER :: UmbilicalDelta = 1.0d-2  ! The tolerance to decide umbilical points
  REAL(KIND=dp), PARAMETER :: GeometryEpsilon = 1.0d-6
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t), POINTER :: Element, BGElement
  TYPE(ElementType_t), POINTER :: ShellElement => NULL()
  TYPE(Nodes_t) :: Nodes
  TYPE(ValueList_t), POINTER :: SolverPars
  TYPE(Variable_t), POINTER :: NodalPDir1, NodalPDir2, NodalPDir3

  LOGICAL :: Found
  LOGICAL :: CurveDataOutput, SavePrincipalAxes, ComputeShellArea
  LOGICAL :: WriteElementalDirector
  LOGICAL :: MacroElements, QuadraticApproximation = .FALSE.
  LOGICAL :: PlateBody, PlanarPoint, UmbilicalPoint
  LOGICAL :: Bubbles, ApplyBubbles

  INTEGER, POINTER :: Indices(:) => NULL()
  INTEGER, POINTER :: VisitsList(:) => NULL()
  INTEGER :: e, i, i0, j, k, m, n, nb, nd, t 
  INTEGER :: Family, Active
  INTEGER :: ShellModelPar, StrainReductionMethod

  REAL(KIND=dp), TARGET :: LocalFrameNodes(MaxPatchNodes,3)
  REAL(KIND=dp) :: LocalSol(6,9)
  REAL(KIND=dp) :: d(3), d1(3), d2(3), d3(3), X1(3), X2(3), e1(3), e2(3), e3(3), o(3), p(3)
  REAL(KIND=dp) :: c, Norm, u, v, X
  REAL(KIND=dp) :: PatchNodes(MaxPatchNodes,2), ZNodes(MaxPatchNodes), TaylorParams(6)
  REAL(KIND=dp) :: BlendingSurfaceArea, ShellModelArea, MappedMeshArea, RefArea

  CHARACTER(LEN=MAX_NAME_LEN) :: OutputFile

  ! Variables for development version: 
  REAL(KIND=dp) :: TotalErr
  REAL(KIND=dp) :: RefWork, Work
  REAL(KIND=dp) :: MaxPDir1Err, MaxPDir2Err, PDir1(3), PDir2(3)
  REAL(KIND=dp) :: Energy(4), MEnergy, SEnergy, BEnergy, Etot

  SAVE Indices
!------------------------------------------------------------------------------  

  ! ---------------------------------------------------------------------------------
  ! PART 0:
  ! Obtain the values of some key parameters and create allocatable variables. 
  ! ---------------------------------------------------------------------------------
  Mesh => GetMesh()
  SolverPars => GetSolverParams()

  ! ---------------------------------------------------------------------------------
  ! The number of unknown fields in the shell model:
  ! ---------------------------------------------------------------------------------
  ShellModelPar = ListGetInteger(SolverPars, 'Variable DOFs', minv=6, maxv=6)

  ! ---------------------------------------------------------------------------------
  ! The choice of strain reduction method:
  ! ---------------------------------------------------------------------------------
  StrainReductionMethod = ListGetInteger(SolverPars, 'Strain Reduction Operator', &
      Found, minv=0,maxv=3)
  IF (.NOT.Found) StrainReductionMethod = AutomatedChoice

  Bubbles = GetLogical(SolverPars, 'Bubbles', Found)

  !-----------------------------------------------------------------------------------
  ! The field variables for saving the orientation of lines of curvature basis
  ! vectors at the nodes: Since the global DOFs are expressed with respect to the
  ! global frame, creating these may not be of interest.
  !-----------------------------------------------------------------------------------
  SavePrincipalAxes = GetLogical(SolverPars, 'Principal Axes Output', Found)
  IF (SavePrincipalAxes) THEN
    NodalPDir1 => VariableGet( Mesh % Variables, 'Principal Coordinate Dir1')
    NodalPDir2 => VariableGet( Mesh % Variables, 'Principal Coordinate Dir2')
    NodalPDir3 => VariableGet( Mesh % Variables, 'Principal Coordinate Dir3')
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

  IF (.NOT. ASSOCIATED(Indices)) ALLOCATE( Indices(Mesh % MaxElementDOFs) )


  ! ---------------------------------------------------------------------------------
  ! PART I: 
  ! Read the director data at the nodes from mesh.director file and/or check the
  ! the integrity of the surface model. This subroutine creates an elementwise 
  ! property 'director' corresponding to the data, if not already available
  ! via reading the director data from the file mesh.elements.data.
  !----------------------------------------------------------------------------------
  WriteElementalDirector = GetLogical(SolverPars, 'Write Elemental Director', Found)
  IF (WriteElementalDirector) THEN
    OutputFile = GetString(SolverPars, 'Elemental Director Output File', Found)
    IF (.NOT. Found) THEN
      OutputFile(1:18) = 'mesh.elements.data'
      OutputFile(19:19) = CHAR(0)
    END IF
  END IF
  CALL ReadSurfaceDirector(Mesh % Name, Mesh % NumberOfNodes, &
      CheckSurfaceOrientation=.TRUE., WriteElementwiseDirector = WriteElementalDirector, &
      ElementwiseDirectorFile = OutputFile)

  ! --------------------------------------------------------------------------------
  ! PART II:
  ! Generate the descriptions of curved element edges for improved geometry 
  ! approximation. The implementation may not be memory efficient as data is 
  ! dublicated for shared element edges with the same director data. Here the
  ! variable CurveDataOutput can be used to output edge data into a file.
  ! With the macro element option we may create additional space curves
  ! corresponding to subtriangulations of quadrilateral elements.
  ! ---------------------------------------------------------------------------------
  CurveDataOutput = GetLogical(SolverPars, 'Edge Curves Output', Found)
  MacroElements = GetLogical(SolverPars, 'Use Macro Elements', Found)
  IF (.NOT.Found) MacroElements = .TRUE.
  CALL CreateCurvedEdges(CurveDataOutput, MacroElements)


  ! ---------------------------------------------------------------------------------
  ! PART III:
  ! Utilize the parametrized edge curves to obtain improved geometry approximation 
  ! via using the finite element blending technique and assemble the discrete shell 
  ! equations. 
  !   Here elementwise frames which give the orientation of the lines of curvature 
  ! at the element centers are defined and then these frames are employed to
  ! generate a refined approximation to lines of curvature parameterization; cf. [1].
  ! ---------------------------------------------------------------------------------
  ! Check whether the area of shell surface should be computed (here this is done
  ! in several ways to check the model integrity):
  ! ---------------------------------------------------------------------------------
  ComputeShellArea = GetLogical(SolverPars, 'Compute Shell Area')
  BlendingSurfaceArea = 0.0d0
  ShellModelArea = 0.0d0
  MappedMeshArea = 0.0d0
  ! The following are just for verification purposes (remove when final)
  TotalErr = 0.0d0

  ! ------------------------------------------------------------------------------
  ! Finally, this is the assembly loop for generating discrete shell equations.
  ! During the assembly, several elementwise properties of geometric nature (principal
  ! directions, elementwise coordinate systems, etc.) are computed and may be saved
  ! as elementwise properties to avoid a later recomputation.
  ! ------------------------------------------------------------------------------
  CALL DefaultInitialize()
  Active = GetNOFActive()  
  DO k=1,Active
    BGElement => GetActiveElement(k)

    Family = GetElementFamily(BGElement)
    IF ( .NOT.(Family == 3 .OR. Family == 4) ) CYCLE

    n  = GetElementNOFNodes()
    nd = GetElementDOFs(Indices)
    nb = GetElementNOFBDOFs()

    !----------------------------------------------------------------------
    ! Bubbles are designed for the lowest-order discretization:
    !----------------------------------------------------------------------
    ApplyBubbles = Bubbles .AND. (nd == Family)
    IF (nb > 0) CALL Fatal('ShellSolver', &
        'Static condensation for p-bubbles is not supported')

    !----------------------------------------------------------------------
    ! The elementwise average of director data for orientation purposes
    ! (check also for body flatness):
    !----------------------------------------------------------------------
    d = AverageDirector(BGElement, n, PlateBody)

    !-------------------------------------------------------------------------
    ! Create an improved geometry approximation by using the finite element 
    ! blending technique and then find the orientation of lines of curvatures 
    ! at the element center. Compute also the coefficients of the Taylor 
    ! polynomial for creating the improved lines of curvature parameterization:
    !-------------------------------------------------------------------------
    CALL LinesOfCurvatureFrame(BGElement, e1=e1, e2=e2, e3=e3, o=o, &
        TaylorParams=TaylorParams, LagrangeNodes=LocalFrameNodes, d=d, &
        PlanarSurface=PlateBody, PlanarPoint=PlanarPoint, &
        UmbilicalPoint=UmbilicalPoint, MacroElement=MacroElements, &
        SaveProperties = .TRUE.) 

    !--------------------------------------------------------------------------
    ! Obtain the final domain for improved lines of curvature parametrization. 
    ! The array PatchNodes will contain a new set of nodes on the plane 
    ! span{e1,e2} so that the final domain can represented via bicubic 
    ! interpolation.
    !--------------------------------------------------------------------------
    CALL LinesOfCurvaturePatch(LocalFrameNodes(1:MaxPatchNodes,1:2), TaylorParams, &
        PatchNodes, Family, ZNodes, PlanarPoint, UmbilicalPoint)

    !LocalFrameNodes(:,3) = ZNodes(:)
    !CALL MappedBGMeshArea(BGElement, LocalFrameNodes, MappedMeshArea)
    !CALL ComputeSurfaceArea(BGElement, BlendingSurfaceArea, MacroElement=.FALSE.)
    !CYCLE

    ! ------------------------------------------------------------------------------
    ! Generate the element stiffness matrix and assemble the local contribution:
    ! -----------------------------------------------------------------------------
    CALL ShellLocalMatrix(BGElement, n, nd+nb, TaylorParams, PatchNodes, &
        e1, e2, e3, o, ShellModelPar, PlanarPoint, UmbilicalPoint, &
        StrainReductionMethod, ApplyBubbles, ShellModelArea, TotalErr)

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

    !----------------------------------------------------------------------
    ! The area computation for the available geometry description:
    !----------------------------------------------------------------------
    IF (ComputeShellArea) CALL MappedBGMeshArea(BGElement, LocalFrameNodes, MappedMeshArea)
    IF (ComputeShellArea) CALL ComputeSurfaceArea(BGElement, BlendingSurfaceArea, MacroElements)
  
  END DO

  CALL DefaultFinishBulkAssembly() 
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()

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



  ! -------------------------------------------------------------------------------------
  ! SOME VERIFICATION OUTPUT if a benchmark case of straight cylindrical shell is solved
  !-----------------------------------------------------------------------------------
  IF (GetLogical(SolverPars, 'Benchmark Problem', Found)) THEN

    Work = 8.0d0*SUM(Solver % Variable % Values(:) * Solver % Matrix % RHS(:))

    SELECT CASE(ListGetInteger(SolverPars, 'Benchmark Case', Found, minv=0,maxv=2))
    CASE(1)
      RefWork = 12.0d0*(1.0d0-(1.0d0/3.0d0)**2)*(1.0d9)**2/7.0d10 * 2.688287959059254d0 * 1.0d-2 ! t=0.01
      !RefWork = 12.0d0*(1.0d0-(1.0d0/3.0d0)**2)*(1.0d9)**2/7.0d10 * 1.828629366566552 * 1.0d-1 ! t=0.1
      PRINT *, 'Relative energy error = ', SQRT(ABS(RefWork-Work)/RefWork)
    CASE(2)
      RefWork = 12.0d0*(1.0d0-(1.0d0/3.0d0)**2)*(1.0d9)**2/7.0d10 * 0.704331198817278d0 * (1.0d-2)**3 ! t=0.01
      PRINT *, 'Relative energy error = ', SQRT(ABS(RefWork-Work)/RefWork)
    END SELECT
    PRINT *, 'Total number of DOFS = ', SIZE(Solver % Variable % Values) 

    IF (ComputeShellArea) THEN
      RefArea = 0.5d0 * PI 
      !RefArea = 4 * (1.0472d0)**2  
      !  RefArea = PI/4.0d0
      PRINT *, 'Relative Error of Model Surface Area = ', ABS(RefArea  - ShellModelArea)/RefArea    
      PRINT *, 'Relative Error of Blending Surface Area = ', ABS(RefArea  - BlendingSurfaceArea)/RefArea
      PRINT *, 'Relative Error of Mapped BG Mesh Area = ', ABS(RefArea  - MappedMeshArea)/RefArea
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
! This subroutine reads mesh.director file arranged as
!
!    node_id1 d_x d_y d_z
!    ...
!    node_idN d_x d_y d_z
!
! to obtain the shell director data at nodes and creates an elementwise property 
! 'director' corresponding to this data. With the optional argument 
! CheckSurfaceOrientation = .TRUE. the integrity of director data is checked 
! (a consistent orientation of the surface over the element). If the file
! mesh.elements.data has been used to specify the director as an elementwise
! property 'director', this subroutine may be used just to perform the orientation 
! check. If both the files exist, the director obtained from mesh.elements.data
! is used. With WriteElementwiseDirector = .TRUE. the director data is written
! as elementwise property to a file whose format conforms with a file 
! mesh.elements.data (this is the default name for the output file, so this
! option can be used to convert mesh.director into mesh.elements.data format).
! 
! TO DO: Implement parallel version of file reading (mesh.elements.data and
!        mesh.director). Allow arbitrary indexing of nodes.
!------------------------------------------------------------------------------
  SUBROUTINE ReadSurfaceDirector( MeshName, NumberOfNodes, &
      CheckSurfaceOrientation, WriteElementwiseDirector, ElementwiseDirectorFile)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    CHARACTER(LEN=MAX_NAME_LEN), INTENT(IN) :: MeshName
    INTEGER, INTENT(IN) :: NumberOfNodes
    LOGICAL, OPTIONAL, INTENT(IN) :: CheckSurfaceOrientation
    LOGICAL, OPTIONAL, INTENT(IN) :: WriteElementwiseDirector
    CHARACTER(LEN=MAX_NAME_LEN), OPTIONAL, INTENT(IN) :: ElementwiseDirectorFile
    !------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: ReadNodalDirectors, DirectorDataCheck, WriteElementsData, Found
    INTEGER :: n, iostat, i, j, k, i0, Active, Family
    REAL(KIND=dp), POINTER :: NodalDirector(:,:)  
    REAL(KIND=dp), POINTER :: DirectorValues(:)
    REAL(KIND=dp) :: d(3), d1(3), d2(3), d3(3), Norm, X1(3), X2(3)
    REAL(KIND=dp) :: e1(3), e2(3), e3(3)
    REAL(KIND=dp) :: ElementDirectors(3*MaxBGElementNodes)
    CHARACTER(LEN=MAX_NAME_LEN) :: DirectorFile, FormatString
    !------------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------
    ! Check whether mesh.director is read:
    ! -----------------------------------------------------------------------------
    n = LEN_TRIM(MeshName)
    DirectorFile(1:n) = MeshName(1:n)
    DirectorFile(n+1:n+1) = '/'
    DirectorFile(n+2:n+14) = 'mesh.director'
    DirectorFile(n+15:n+15) = CHAR(0)

    INQUIRE(FILE = DirectorFile(1:n+15), EXIST = Found)
    IF ( Found ) THEN    
      ReadNodalDirectors = .TRUE.
      CALL AllocateArray(NodalDirector, NumberOfNodes, 3, 'ReadSurfaceDirector', &
          'NodalDirector array could not be allocated')
    ELSE
      ReadNodalDirectors = .FALSE.
      CALL Info('ReadSurfaceDirector', 'shell director data is read from mesh.elements.data file', &
          Level=8)
    END IF

    IF (ReadNodalDirectors) THEN
      OPEN( 10, FILE = DirectorFile(1:n+15), status='OLD', IOSTAT = iostat )
      IF ( iostat /= 0 ) THEN
        CALL Fatal( 'ReadSurfaceDirector', 'Opening mesh.director file failed.')     
      ELSE
        DO i=1,NumberOfNodes
          READ( 10,*,IOSTAT=iostat) k, d

          IF (k /= i) CALL Fatal('mesh.director', &
              'Trivial correspondence between rows and node numbers assumed currently')

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
        ! -------------------------------------------------------------------
        ! If mesh.elements.data has defined the director, respect that data:
        ! -------------------------------------------------------------------
        DirectorValues => GetElementProperty('director', Element)
        IF (ASSOCIATED(DirectorValues)) CYCLE

        n  = GetElementNOFNodes()     
        DO i=1,n
          i0 = (i-1)*3
          ElementDirectors(i0+1:i0+3) = NodalDirector(Element % NodeIndexes(i),1:3)
        END DO
        CALL SetElementProperty('director', ElementDirectors(1:3*n), Element)
      END DO
    END IF

    ! ---------------------------------------------------------------------------
    ! Perform an additional check that the director data defines a properly 
    ! oriented model. All directors should point to the same side of the surface.
    ! ----------------------------------------------------------------------------
    IF (PRESENT(CheckSurfaceOrientation)) THEN
      DirectorDataCheck = CheckSurfaceOrientation
    ELSE
      DirectorDataCheck = .FALSE.   
    END IF

    IF (DirectorDataCheck) THEN
      Active = GetNOFActive()
      DO k=1,Active
        Element => GetActiveElement(k)
        Family = GetElementFamily(Element)
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
        DirectorValues => GetElementProperty('director', Element)
        IF (ASSOCIATED(DirectorValues)) THEN
          IF (SIZE(DirectorValues) < 3*n) CALL Fatal('ReadSurfaceDirector', &
              'Elemental director data is not associated with all nodes')
          d1(1:3) = DirectorValues(1:3)
          c = DOT_PRODUCT(d1,e3)
        ELSE
          CALL Fatal( 'ReadSurfaceDirector', 'Elemental director data is not associated')
        END IF

        DO j=2,n
          i0 = (j-1)*3
          d2(1:3) = DirectorValues(i0+1:i0+3)
          IF ( (c * DOT_PRODUCT(d2,e3)) < 0.0d0 ) THEN
            PRINT *, 'Element indices=', Element % NodeIndexes(1:n)
            PRINT *, 'Reference normal =', e3(:)
            PRINT *, 'Node Index = ', j, Element % NodeIndexes(j)
            PRINT *, 'Director =', d2(:)
            CALL Fatal( 'ReadSurfaceDirector', &
                'Director data does not define a unique upper/lower surface.')
          END IF
        END DO
      END DO
    END IF

    ! ---------------------------------------------------------------------
    ! Write the director data as elementwise property to a file whose
    ! format conforms with a file mesh.elements.data. By default
    ! the file name mesh.elements.data is used. This never overwrites
    ! an existing file.
    ! ---------------------------------------------------------------------
    IF (PRESENT(WriteElementwiseDirector)) THEN
      WriteElementsData = WriteElementwiseDirector
    ELSE
      WriteElementsData = .FALSE.
    END IF
    IF (WriteElementsData) THEN
      n = LEN_TRIM(MeshName)
      DirectorFile(1:n) = MeshName(1:n)
      DirectorFile(n+1:n+1) = '/'
      IF ( PRESENT(ElementwiseDirectorFile) ) THEN
        j = LEN_TRIM(ElementwiseDirectorFile)
        DirectorFile(n+2:n+1+j) = ElementwiseDirectorFile(1:j)
        DirectorFile(n+2+j:n+2+j) = CHAR(0)
        n = n+2+j
      ELSE
        DirectorFile(n+2:n+19) = 'mesh.elements.data'
        DirectorFile(n+20:n+20) = CHAR(0)
        n = n+20
      END IF

      INQUIRE(FILE = DirectorFile(1:n), EXIST = Found)
      IF (Found) THEN
        CALL Info('ReadSurfaceDirector', &
            'a file for director output exists: write rejected', Level=5)
      ELSE
        OPEN( 10, FILE = DirectorFile(1:n), status='NEW', IOSTAT = iostat )        
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
 
            WRITE(FormatString(1:1),'(A1)') '('
            IF (3*n < 10) THEN
              WRITE(FormatString(2:2),'(A1)') TRIM(I2S(3*n))
              i0 = 2
            ELSE
              WRITE(FormatString(2:3),'(A2)') TRIM(I2S(3*n))
              i0 = 3
            END IF
            WRITE(FormatString(i0+1:i0+1),'(A1)') '('
            WRITE(FormatString(i0+2:i0+10),'(A9)') '2x,E22.15'
            WRITE(FormatString(i0+11:i0+12),'(A2)') '))'

            WRITE(10,'(A8,I0)') 'element:', k
            WRITE(10,'(A9)',ADVANCE='NO') 'director:'
            WRITE(10,FormatString(1:i0+12)) DirectorValues(1:3*n)
            WRITE(10,'(A3)') 'end'
          ELSE
            CALL Fatal( 'ReadSurfaceDirector', 'Elemental director data is not associated')
          END IF
        END DO
        CLOSE(10)
      END IF
    END IF

    IF (ReadNodalDirectors) DEALLOCATE(NodalDirector)
!------------------------------------------------------------------------------
  END SUBROUTINE ReadSurfaceDirector
!------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------
! Use nodal directors, which are retrieved as elementwise property 'director', 
! to create the parametrizations of curved edges for the Hermite interpolation.
! The edge curve data are written as elementwise properties 'edge frames' and
! 'edge parameters'.
!
! TO DO: Make MacroElement option functional for second-order nodal director data.
!        Add 6-node triangles
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
    REAL(KIND=dp) :: EdgeFramesData(FrameDataSize*MaxNumberOCurves)
    REAL(KIND=dp) :: EdgeCurveParams(CurveDataSize2*MaxNumberOCurves)
!-------------------------------------------------------------------------------

    A => FrameData(:,:)
    cpars => CurveData(:)

    ! -------------------------------------------------------------------------------
    ! We may consider 404 as a macroelement for a subtriangulation to ensure that 
    ! fourth-order accurate approximation in L2 can be obtained. 
    ! TO DO: Alter the list of edges parametrized for 9-node director data
    ! -------------------------------------------------------------------------------
    IF (PRESENT(MacroElements)) THEN
      Subtriangulation = MacroElements
    ELSE
      Subtriangulation = .FALSE.
    END IF

    ! Write edge curve parameters to a file:
    ! ------------------------------------------------------------------
    IF (FileOutput) OPEN( 10, FILE = 'edgecsys.dat', status='REPLACE')

    Active = GetNOFActive()
    DO k=1,Active
      Element => GetActiveElement(k)
      DirectorValues => GetElementProperty('director', Element)
      CALL GetElementNodes( Nodes )

      Family = GetElementFamily(Element)
      SELECT CASE(Family)
      CASE(3)
        IF (Element % TYPE % NumberOfNodes /= 3) &
            CALL Fatal('CreateCurvedEdges', 'Only 3-node triangles possible')
        QuadraticGeometryData = .FALSE.
        CurveDataSize = CurveDataSize1
        EdgesParametrized = Family
      CASE(4)
        QuadraticGeometryData = Element % TYPE % NumberOfNodes == 9
        IF (QuadraticGeometryData) THEN
          EdgesParametrized = 6
          CurveDataSize = CurveDataSize2
        ELSE
          CurveDataSize = CurveDataSize1
          IF (Subtriangulation) THEN
            EdgesParametrized = 6
          ELSE
            EdgesParametrized = Family
          END IF
        END IF
      CASE DEFAULT
        CYCLE
      END SELECT

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
            IF (Subtriangulation .AND. (Element % TYPE % NumberOfNodes == 4)) THEN
              i = 1
              j = 3
            ELSE
              i = 8
              j = 6
              l = 9
            END IF
          CASE(6)
            IF (Subtriangulation .AND. (Element % TYPE % NumberOfNodes == 4)) THEN
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
          ! -----------------------------------------------------------
          X3(1) = Nodes % x(l)
          X3(2) = Nodes % y(l)
          X3(3) = Nodes % z(l)

          i0 = (l-1)*3
          d3(1:3) = DirectorValues(i0+1:i0+3)

          CALL EdgeFrame(X1, X2, d1, d2, A, cpars, X3, d3)
        ELSE
          CALL EdgeFrame(X1, X2, d1, d2, A, cpars)

          IF (FileOutput .AND. EdgesParametrized == 4) THEN
            WRITE( 10, '(2I5,18e23.15)',ADVANCE='NO') v1, v2, A(1,1), A(2,1), &
                A(3,1), A(1,2), A(2,2), A(3,2), A(1,3), A(2,3), A(3,3), A(1,4), &
                A(2,4), A(3,4), cpars(1), cpars(1), cpars(2), cpars(3), cpars(4), &
                cpars(5)
            WRITE( 10, *) ''
          END IF
        END IF

        ! Prepare for writing the edge curve data as element properties:
        !---------------------------------------------------------------
        i0 = (e-1)*FrameDataSize
        DO j=1,4
          j0 = i0 + (j-1)*3
          EdgeFramesData(j0+1:j0+3) = A(1:3,j)
        END DO

        i0 = (e-1)*CurveDataSize
        EdgeCurveParams(i0+1:i0+CurveDataSize) = cpars(1:CurveDataSize)

      END DO

      ! Write the edge curve data as element properties:
      !----------------------------------------------------
      CALL SetElementProperty('edge frames', &
          EdgeFramesData(1:FrameDataSize*EdgesParametrized), Element) 
      CALL SetElementProperty('edge parameters', &
          EdgeCurveParams(1:CurveDataSize*EdgesParametrized), Element)       

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
    REAL(KIND=dp), POINTER, INTENT(OUT) :: A(:,:), cpars(:)
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: X3(3), d3(3)
!------------------------------------------------------------------------------
    LOGICAL :: WithThreeNodes
    REAL(KIND=dp) :: d(3), ex(3), ey(3), ez(3), X0(3), v21(3), b(3)
    REAL(KIND=dp) :: r1(3), r2(3), r3(3), t1(3), t2(3), t3(3), Norm
!------------------------------------------------------------------------------
    WithThreeNodes = PRESENT(X3) .AND. PRESENT(d3)

    IF (WithThreeNodes) THEN
      X0(:) = X3(:)

      Norm = SQRT(SUM(d3(1:3)**2))
      ez(1:3) = 1.0d0/Norm * d3(1:3) 

      v21(:) = X2(:) - X1(:)
      b(:) = v21(:) - DOT_PRODUCT(v21,ez)*ez(:)
      Norm = SQRT(SUM(b(1:3)**2))
      ex(1:3) = 1.0d0/Norm * b(1:3)     
      ey(:) = CrossProduct(ez,ex)

      cpars(1) = DOT_PRODUCT(v21,ex)
      IF (cpars(1) < 0.0d0) &
          CALL Fatal('EdgeFrame','Negative edge length obtained')

      r1(:) = X1(:) - X0(:)
      r2(:) = X2(:) - X0(:)

      cpars(2) = DOT_PRODUCT(r1,ey)               ! The y-coordinate of the vertex X1
      cpars(3) = DOT_PRODUCT(r2,ey)               ! The y-coordinate of the vertex X2
      cpars(4) = DOT_PRODUCT(r1,ez)               ! The z-coordinate of the vertex X1
      cpars(5) = DOT_PRODUCT(r2,ez)               ! The z-coordinate of the vertex X2      
      
      t1(:) = CrossProduct(ey,d1)
      t2(:) = CrossProduct(ey,d2)
      t3(:) = CrossProduct(ey,d3)
      
      cpars(6) = DOT_PRODUCT(t1,ez)/DOT_PRODUCT(t1,ex);  ! The angle parameter for z(x) at the vertex X1
      cpars(7) = DOT_PRODUCT(t2,ez)/DOT_PRODUCT(t2,ex);  ! The angle parameter for z(x) at the vertex X2      

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
          CALL Fatal('EdgeFrame','Negative edge length obtained')

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

!----------------------------------------------------------------------------
! This function produces the covariant basis {a_i}, the first and second
! fundamental forms A and B, the determinant of the metric surface tensor detA
! and the global coordinates of the point on the blending surface when the 
! reference element coordinates u and v are used as curvilinear coordinates on
! the blending surface. The necessary edge curve data for creating the blending 
! surface must be contained as elementwise properties 'edge frames' and 
! 'edge parameters'. The optional arguments MacroElement and BubbleDOFs
! can be used to augment the serendipity approximation by an additional bubble 
! part to ensure optimal accuracy.
! TO DO: Complement and clean the implementation when the initial data
!        is defined over second-order Lagrange elements
!-----------------------------------------------------------------------------  
  FUNCTION BlendingSurfaceInfo( Element, Nodes, u, v, &
      deta, a1, a2, a3, A, B, x, MacroElement, BubbleDOFs ) RESULT(stat)
!----------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Element_t), POINTER, INTENT(IN) :: Element  !< Element structure
    TYPE(Nodes_t), INTENT(IN) :: Nodes               !< Data corresponding to the element nodes
    REAL(KIND=dp), INTENT(IN) :: u                   !< 1st reference element coordinate
    REAL(KIND=dp), INTENT(IN) :: v                   !< 2nd coordinate
    REAL(KIND=dp), INTENT(OUT) :: deta               !< The determinant of the surface metric tensor
    REAL(KIND=dp), INTENT(OUT) :: a1(3), a2(3)       !< The covariant surface basis vectors
    REAL(KIND=dp), INTENT(OUT) :: a3(3)              !< The base vector normal to the surface
    REAL(KIND=dp), INTENT(OUT) :: A(2,2)             !< The covariant components of the metric surface tensor at (u,v)  
    REAL(KIND=dp), INTENT(OUT) :: B(2,2)             !< The covariant components of the second fundamental form at (u,v)
    REAL(KIND=dp), INTENT(OUT) :: x(3)               !< Blending surface point corresponding to (u,v): x=x(u,v)
    LOGICAL, OPTIONAL, INTENT(IN) :: MacroElement    !< Use macroelement strategy to add a bubble part 
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: BubbleDOFs(4,3)  !< Coefficients for bubble basis functions
    LOGICAL :: Stat                                  !< A dummy status variable at the moment
!----------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: GElement => NULL()
    LOGICAL :: QuadraticGeometryData, Subtriangulation
    INTEGER :: i, j, e, n, q, i0, EdgesParametrized, CurveDataSize, Family
    REAL(KIND=dp), POINTER :: FrameData(:), EdgeParams(:), BubbleValues(:)
    REAL(KIND=dp) :: Basis(MaxPatchNodes), dBasis(MaxPatchNodes,2), ddBasis(4,2,2)
    REAL(KIND=dp) :: BubbleCoeff(4,3)
    REAL(KIND=dp) :: ex(3), ey(3), ez(3), X0(3)
    REAL(KIND=dp) :: d(CurveDataSize2)
    REAL(KIND=dp) :: h, s, t, w, xe, c(3), deltac(3), dc(3), ddc(3)
    REAL(KIND=dp) :: b12(3), db12(3), h12, d1h12, d2h12, ddh12, dsdu, dsdv
    REAL(KIND=dp) :: f1, f2, f3, f4, f5, f6
    REAL(KIND=dp) :: HermBasis(6), dHermBasis(6), ddHermBasis(6)
    REAL(KIND=dp) :: df1, df2, df3, df4, df5, df6, ddf1, ddf2, ddf3, ddf4, ddf5, ddf6
    REAL(KIND=dp) :: f, df, ddf, fy, dfy, ddfy, fz, dfz, ddfz
    REAL(KIND=dp) :: h1, h2, h3, dh1, dh2, dh3, ddh1, ddh2, ddh3
    REAL(KIND=dp) :: r1(3), r2(3)
    REAL(KIND=dp) :: d1a1(3), d2a1(3), d2a2(3)
    REAL(KIND=dp) :: L1, L2, dL1, dL2
    REAL(KIND=dp) :: Norm

    SAVE GElement
!----------------------------------------------------------------------------
    Family = GetElementFamily(Element)

    SELECT CASE(Family)  
    CASE(3)
      IF (Element % TYPE % NumberOfNodes /= 3) &
          CALL Fatal('BlendingSurfaceInfo', 'Only 3-node triangles possible')
      QuadraticGeometryData = .FALSE.
      Subtriangulation = .FALSE.
      EdgesParametrized = 3
      CurveDataSize = CurveDataSize1
    CASE(4)
      IF (PRESENT(BubbleDOFs)) THEN
        Subtriangulation = .TRUE.
        BubbleCoeff(1:4,1:3) = BubbleDOFs(1:4,1:3)
      ELSE
        IF (PRESENT(MacroElement)) THEN
          Subtriangulation = MacroElement
          IF (Subtriangulation) THEN
            BubbleValues => GetElementProperty('bubble dofs', Element)
            IF (ASSOCIATED(BubbleValues)) THEN
              IF (SIZE(BubbleValues) < 12) CALL Fatal('BlendingSurfaceInfo', &
                  'Bubble dofs data missing')
              BubbleCoeff(1:4,1) = BubbleValues(1:4)
              BubbleCoeff(1:4,2) = BubbleValues(5:8)
              BubbleCoeff(1:4,3) = BubbleValues(9:12)
            ELSE
              CALL Fatal('BlendingSurfaceInfo',&
                  'Bubble DOFs are not found as elementwise properties')
            END IF
          END IF
        ELSE
          Subtriangulation = .FALSE.
        END IF
      END IF

      ! Check whether the edges have mid-nodes:
      QuadraticGeometryData = Element % TYPE % NumberOfNodes == 9
      IF (QuadraticGeometryData) THEN
        EdgesParametrized = 6
        CurveDataSize = CurveDataSize2
      ELSE
        EdgesParametrized = 4
        CurveDataSize = CurveDataSize1
      END IF
    CASE DEFAULT
      CALL Fatal('BlendingSurfaceInfo', 'Only quads and triangles can be handled')     
    END SELECT

    !-----------------------------------------------------------------------
    ! Retrive parametrizations of curved edges:
    !------------------------------------------------------------------------
    FrameData => GetElementProperty('edge frames', Element) 
    EdgeParams => GetElementProperty('edge parameters', Element)   

    IF (ASSOCIATED(FrameData)) THEN
      IF (Subtriangulation) THEN
        IF (SIZE(FrameData) < (EdgesParametrized+2)*FrameDataSize) &
            CALL Fatal('BlendingSurfaceInfo','Frame data are not associated with all edges')        
      ELSE
        IF (SIZE(FrameData) < EdgesParametrized*FrameDataSize) &
            CALL Fatal('BlendingSurfaceInfo','Frame data are not associated with all edges')
      END IF
    ELSE
      CALL Fatal('BlendingSurfaceInfo','edge frames data could not be retrieved')
    END IF

    IF (ASSOCIATED(EdgeParams)) THEN
      IF (Subtriangulation) THEN
        IF (SIZE(EdgeParams) < (EdgesParametrized+2)*CurveDataSize) &
            CALL Fatal('BlendingSurfaceInfo','edge parameters are not associated with all edges')       
      ELSE
        IF (SIZE(EdgeParams) < EdgesParametrized*CurveDataSize) &
            CALL Fatal('BlendingSurfaceInfo','edge parameters are not associated with all edges')
      END IF
    ELSE
      CALL Fatal('BlendingSurfaceInfo','edge frames data could not be retrieved')
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
    X(1) = SUM( Nodes % x(1:n) * Basis(1:n) )
    X(2) = SUM( Nodes % y(1:n) * Basis(1:n) )
    X(3) = SUM( Nodes % z(1:n) * Basis(1:n) )
    ! XInit(:) = X(:)
    a1(1) = SUM( dBasis(1:n,1) * Nodes % x(1:n) )
    a1(2) = SUM( dBasis(1:n,1) * Nodes % y(1:n) )
    a1(3) = SUM( dBasis(1:n,1) * Nodes % z(1:n) )
    a2(1) = SUM( dBasis(1:n,2) * Nodes % x(1:n) )
    a2(2) = SUM( dBasis(1:n,2) * Nodes % y(1:n) )
    a2(3) = SUM( dBasis(1:n,2) * Nodes % z(1:n) )

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
      
      i0 = (e-1)*FrameDataSize
      ex(1:3) = FrameData(i0+FrameBasis1)
      ey(1:3) = FrameData(i0+FrameBasis2)
      ez(1:3) = FrameData(i0+FrameBasis3)
      X0(1:3) = FrameData(i0+FrameOrigin)

      i0 = (e-1)*CurveDataSize
      d(1:CurveDataSize) = EdgeParams(i0+1:i0+CurveDataSize)
      h = d(1) ! The length parameter

      SELECT CASE(Family)
      CASE(3)
        !-------------------------------------------------------------------------
        ! First define edge curve parameter and a part of the blending function
        !-------------------------------------------------------------------------
        SELECT CASE(e)
        CASE(1)
          s = Basis(2) - Basis(1)
          h1 = Basis(2) * Basis(1)
        CASE(2)
          s = Basis(3) - Basis(2)
          h1 = Basis(3) * Basis(2)
        CASE(3)
          s = Basis(1) - Basis(3)
          h1 = Basis(1) * Basis(3)
        END SELECT

        !-------------------------------------------------------------------------
        ! The dimensional edge parameterization coordinate corresponding to
        ! the dimensionless coordinate s
        !--------------------------------------------------------------------------
        xe = -0.25d0*(1.0d0-s)*h + 0.25d0*(1.0d0+s)*h

        !-------------------------------------------------------------------------
        ! The ordinary third-order Hermite interpolation basis and
        ! linear basis functions and their derivatives on an edge [-1,1]
        !-------------------------------------------------------------------------         
        CALL HermiteBasis(s, h, HermBasis(1:4), dHermBasis(1:4), ddHermBasis(1:4))
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
        f = SUM( d(2:5) * HermBasis(1:4) )
        df = SUM( d(2:5) * dHermBasis(1:4) )
        ddf = SUM( d(2:5) * ddHermBasis(1:4) )
        c(1:3) = X0(1:3) + xe * ex(1:3) + f * ez(1:3)
        dc(1:3) = 0.5d0*h * ex(1:3) + df * ez(1:3)
        ddc(1:3) = ddf * ez(1:3)

        !---------------------------------------------------------------------------
        ! The contributions of the blending functions to the position vector, to the
        ! surface base vectors and to the derivatives of the surface base vectors
        !---------------------------------------------------------------------------
        SELECT CASE(e)
        CASE(1)
          r1(1) = Nodes % x(1)
          r1(2) = Nodes % y(1)
          r1(3) = Nodes % z(1)
          r2(1) = Nodes % x(2)
          r2(2) = Nodes % y(2)
          r2(3) = Nodes % z(2)

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
          r1(1) = Nodes % x(2)
          r1(2) = Nodes % y(2)
          r1(3) = Nodes % z(2)
          r2(1) = Nodes % x(3)
          r2(2) = Nodes % y(3)
          r2(3) = Nodes % z(3)

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
          r1(1) = Nodes % x(3)
          r1(2) = Nodes % y(3)
          r1(3) = Nodes % z(3)
          r2(1) = Nodes % x(1)
          r2(2) = Nodes % y(1)
          r2(3) = Nodes % z(1)

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
        ! First define edge orientation convention and retrive parameters for
        ! representing the curved edge
        !-------------------------------------------------------------------------
        SELECT CASE(e)
        CASE(1,3,5)
          s = u
          t = v
        CASE(2,4,6)
          s = v
          t = u
        END SELECT

        !-------------------------------------------------------------------------
        ! The dimensional edge parameterization coordinate corresponding to
        ! the dimensionless coordinate s (depending on the edge either s=u or s=v)
        !--------------------------------------------------------------------------
        xe = -0.25d0*(1.0d0-s)*h + 0.25d0*(1.0d0+s)*h
        !-------------------------------------------------------------------------
        ! The basis functions for the edge curve expansion in terms of s and
        ! their derivatives
        !------------------------------------------------------------------------- 
        IF (QuadraticGeometryData) THEN
          ! Fifth-order Hermite basis: Dofs are listed as p(v1), p(v2), p(v3),
          ! Dp(v1)[e1], Dp(v2)[e1], Dp(v3)[e1]
          ! TO DO: the subroutine Hermite should return these
          f1 = s ** 2 - 0.5D1 / 0.4D1 * s ** 3 - s ** 4 / 0.2D1 + 0.3D1 / 0.4D1 * s ** 5
          f2 = s ** 2 + 0.5D1 / 0.4D1 * s ** 3 - s ** 4 / 0.2D1 - 0.3D1 / 0.4D1 * s ** 5
          f3 = s ** 4 - 2.0D0 * s ** 2 + 1.0D0
          f4 = h/2.0d0 * (s ** 2 / 0.4D1 - s ** 3 / 0.4D1 - s ** 4 / 0.4D1 + s ** 5 / 0.4D1)
          f5 = h/2.0d0 * (-s ** 2 / 0.4D1 - s ** 3 / 0.4D1 + s ** 4 / 0.4D1 + s ** 5 /0.4D1)
          f6 = h/2.0d0 * (s ** 5 - 2.0D0 * s ** 3 + s)

          df1 = 2.0D0 * s - 0.15D2 / 0.4D1 * s ** 2 - 2.0D0 * s ** 3 + 0.15D2 / 0.4D1 * s ** 4
          df2 = 2 * s + 0.15D2 / 0.4D1 * s ** 2 - 2.0D0 * s ** 3 - 0.15D2 / 0.4D1 * s ** 4
          df3 = 4.0D0 * s ** 3 - 4.0D0 * s
          df4 = h/2.0d0 * (s / 0.2D1 - 0.3D1 / 0.4D1 * s ** 2 - s ** 3 + 0.5D1 / 0.4D1 * s ** 4)
          df5 = h/2.0d0 * (-s / 0.2D1 - 0.3D1 / 0.4D1 * s ** 2 + s ** 3 + 0.5D1 / 0.4D1* s ** 4)
          df6 = h/2.0d0 * (5.0D0 * s ** 4 - 6.0D0 * s ** 2 + 1.0D0)

          ddf1 = 0.2D1 - 0.15D2 / 0.2D1 * s - 0.6D1 * s ** 2 + 0.15D2 * s ** 3
          ddf2 = 0.2D1 + 0.15D2 / 0.2D1 * s - 0.6D1 * s ** 2 - 0.15D2 * s ** 3
          ddf3 = 0.12D2 * s ** 2 - 0.4D1
          ddf4 = h/2.0d0 * (0.1D1 / 0.2D1 - 0.3D1 / 0.2D1 * s - 0.3D1 * s ** 2 + 0.5D1 * s ** 3)
          ddf5 = h/2.0d0 * (-0.1D1 / 0.2D1 - 0.3D1 / 0.2D1 * s + 0.3D1 * s ** 2 + 0.5D1 * s ** 3)
          ddf6 = h/2.0d0 * (0.2D2 * s ** 3 - 0.12D2 * s)
        ELSE
          ! The ordinary third-order Hermite interpolation basis:
          CALL HermiteBasis(s, h, HermBasis(1:4), dHermBasis(1:4), ddHermBasis(1:4))
        END IF
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
          fz = d(4)*f1 + d(5)*f2 + d(6)*f4 + d(7)*f5 ! + d(8)*f6
          dfz = d(4)*df1 + d(5)*df2 + d(6)*df4 + d(7)*df5 ! + d(8)*df6
          ddfz = d(4)*ddf1 + d(5)*ddf2 + d(6)*ddf4 + d(7)*ddf5 !+ d(8)*ddf6

          fy = d(2)*f1 + d(3)*f2 ! + d(9)*f4 + d(10)*f5 + d(11)*f6
          dfy = d(2)*df1 + d(3)*df2 ! + d(9)*df4 + d(10)*df5 + d(11)*df6
          ddfy = d(2)*ddf1 + d(3)*ddf2 ! + d(9)*ddf4 + d(10)*ddf5 + d(11)*ddf6

          c(1:3) = X0(1:3) + xe * ex(1:3) + fy * ey(1:3) + fz * ez(1:3)
          dc(1:3) = 0.5d0*h * ex(1:3) + dfy * ey(1:3) + dfz * ez(1:3)
          ddc(1:3) = ddfy * ey(1:3) + ddfz * ez(1:3)

        ELSE
          f = SUM( d(2:5) * HermBasis(1:4) )
          df = SUM( d(2:5) * dHermBasis(1:4) )
          ddf = SUM( d(2:5) * ddHermBasis(1:4) )
          c(1:3) = X0(1:3) + xe * ex(1:3) + f * ez(1:3)
          dc(1:3) = 0.5d0*h * ex(1:3) + df * ez(1:3)
          ddc(1:3) = ddf * ez(1:3)
        END IF

        !---------------------------------------------------------------------------
        ! The contributions of the blending functions to the position vector, to the
        ! surface base vectors and to the derivatives of the surface base vectors
        !---------------------------------------------------------------------------
        SELECT CASE(e)
        CASE(1)
          r1(1) = Nodes % x(1)
          r1(2) = Nodes % y(1)
          r1(3) = Nodes % z(1)
          r2(1) = Nodes % x(2)
          r2(2) = Nodes % y(2)
          r2(3) = Nodes % z(2)

          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)
          X(1:3) = X(1:3) + h1 * b12(1:3)

          a1(1:3) = a1(1:3) + h1 * db12(1:3)
          a2(1:3) = a2(1:3) + dh1 * b12(1:3)
          d1a1(1:3) = d1a1(1:3) + h1 * ddc(1:3)
          d2a1(1:3) = d2a1(1:3) + dh1 * db12(1:3)
          d2a2(1:3) = d2a2(1:3) + ddh1 * b12(1:3)

        CASE(2)
          r1(1) = Nodes % x(2)
          r1(2) = Nodes % y(2)
          r1(3) = Nodes % z(2)
          r2(1) = Nodes % x(3)
          r2(2) = Nodes % y(3)
          r2(3) = Nodes % z(3)

          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)
          X(1:3) = X(1:3) + h2 * b12(1:3)

          a1(1:3) = a1(1:3) + dh2 * b12(1:3)
          a2(1:3) = a2(1:3) + h2 * db12(1:3)
          d1a1(1:3) = d1a1(1:3) + ddh2 * b12(1:3)
          d2a1(1:3) = d2a1(1:3) + dh2 * db12(1:3)
          d2a2(1:3) = d2a2(1:3) + h2 * ddc(1:3)

        CASE(3)
          r1(1) = Nodes % x(4)
          r1(2) = Nodes % y(4)
          r1(3) = Nodes % z(4)
          r2(1) = Nodes % x(3)
          r2(2) = Nodes % y(3)
          r2(3) = Nodes % z(3)

          b12(1:3) = c(1:3) - L1 * r1(1:3) - L2 * r2(1:3)
          db12(1:3) = dc(1:3) - dL1 * r1(1:3) - dL2 * r2(1:3)
          X(1:3) = X(1:3) + h2 * b12(1:3)

          a1(1:3) = a1(1:3) + h2 * db12(1:3)
          a2(1:3) = a2(1:3) + dh2 * b12(1:3)
          d1a1(1:3) = d1a1(1:3) + h2 * ddc(1:3)
          d2a1(1:3) = d2a1(1:3) + dh2 * db12(1:3)
          d2a2(1:3) = d2a2(1:3) + ddh2 * b12(1:3)

        CASE(4)
          r1(1) = Nodes % x(1)
          r1(2) = Nodes % y(1)
          r1(3) = Nodes % z(1)
          r2(1) = Nodes % x(4)
          r2(2) = Nodes % y(4)
          r2(3) = Nodes % z(4)

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

    IF (Subtriangulation) THEN
      ! --------------------------------------------------------------------------
      ! Add contributions which relate to augmenting the serendipity
      ! approximation by bubbles, so that an approximation Q3 is obtained.
      ! --------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED(GElement) ) GElement => AllocateElement()
      GElement % Type => GetElementType(416,.FALSE.)

      Basis = 0.0d0
      DO q=13,16
        Basis(q) = 1.0d0
        ddBasis(q-12,1:2,1:2) = SecondDerivatives2D(GElement, Basis, u, v)
        Basis(q) = 0.0d0
      END DO
      CALL NodalFirstDerivatives2D(dBasis, GElement, u, v)
      CALL NodalBasisFunctions2D(Basis, GElement, u, v)

      DO q=1,3
        X(q) = X(q) + SUM( BubbleCoeff(1:4,q) * Basis(13:16) )
        a1(q) = a1(q) + SUM( BubbleCoeff(1:4,q) * dBasis(13:16,1) )
        a2(q) = a2(q) + SUM( BubbleCoeff(1:4,q) * dBasis(13:16,2) )
        d1a1(q) = d1a1(q) + SUM( BubbleCoeff(1:4,q) * ddBasis(1:4,1,1) )
        d2a2(q) = d2a2(q) + SUM( BubbleCoeff(1:4,q) * ddBasis(1:4,2,2) )
        d2a1(q) = d2a1(q) + SUM( BubbleCoeff(1:4,q) * ddBasis(1:4,1,2) )
      END DO
    END IF

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
! This subroutine gives the Hermite basis functions over [-1,1] and their
! derivatives up to the second order.
!------------------------------------------------------------------------------
  SUBROUTINE HermiteBasis(s, h, Basis, dBasis, ddBasis)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: s           ! Coordinate
    REAL(KIND=dp), INTENT(IN) :: h           ! A scale factor
    REAL(KIND=dp), INTENT(OUT) :: Basis(:)   ! Basis functions
    REAL(KIND=dp), INTENT(OUT) :: dBasis(:)  ! The first-order derivatives
    REAL(KIND=dp), INTENT(OUT) :: ddBasis(:) ! The second-order derivatives
!------------------------------------------------------------------------------
    IF (SIZE(Basis) < 4 .OR. SIZE(dBasis) < 4 .OR.  SIZE(ddBasis) < 4) &
        CALL Fatal('HermiteBasis', 'Too small arrays for basis functions')

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

!------------------------------------------------------------------------------
  END SUBROUTINE HermiteBasis
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! This subroutine finds the desired position of the blending surface at four
! internal nodes corresponding to bubble basis functions of the Q3 space via
! the macro element strategy. The nodal difference between the desired position 
! and the serendipity approximation is returned via the variable BubbleNodesDelta.
!------------------------------------------------------------------------------
  SUBROUTINE FindBubbleNodesQuad(Element, Nodes, BubbleNodesDelta)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element 
    TYPE(Nodes_t), INTENT(IN) :: Nodes
    REAL(KIND=dp), INTENT(OUT) :: BubbleNodesDelta(4,3)
!------------------------------------------------------------------------------
    LOGICAL :: Stat
    INTEGER :: CurveDataSize, i, j, k, e, i0
    REAL(KIND=dp), POINTER :: FrameData(:), EdgeParams(:) 
    REAL(KIND=dp) :: s, u, v
    REAL(KIND=dp) :: ex(3), ey(3), ez(3), X0(3), c(3)
    REAL(KIND=dp) :: h, xe, f, f2, f3, f4, f5
    REAL(KIND=dp) :: d(CurveDataSize1)
    REAL(KIND=dp) :: a1(3), a2(3), a3(3), a(2,2), Deta, b(2,2), p(3)
!------------------------------------------------------------------------------
    CurveDataSize = CurveDataSize1

    !-----------------------------------------------------------------------
    ! Retrive parametrizations of curved edges:
    !------------------------------------------------------------------------
    FrameData => GetElementProperty('edge frames', Element) 
    EdgeParams => GetElementProperty('edge parameters', Element)
    !------------------------------------------------------------------------
    ! Note that the correct sizes of these data arrays are checked afterwards
    ! in the function BlendingSurfaceInfo, so avoid the size check here:
    !------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED(FrameData) .OR.  .NOT. ASSOCIATED(EdgeParams) ) &
        CALL Fatal('FindBubbleNodesQuad', 'Elemental properties missing')
    !-----------------------------------------------------------------------
    ! Find the desired place of the blending surface via the additional
    ! edges of the subtriangulation. We need two position evaluations
    ! per additional edge of the subtriangulation: 
    !-----------------------------------------------------------------------
    DO e=5,6
      i0 = (e-1)*FrameDataSize
      ex(1:3) = FrameData(i0+FrameBasis1)
      ey(1:3) = FrameData(i0+FrameBasis2)
      ez(1:3) = FrameData(i0+FrameBasis3)
      X0(1:3) = FrameData(i0+FrameOrigin)

      i0 = (e-1)*CurveDataSize
      d(1:CurveDataSize) = EdgeParams(i0+1:i0+CurveDataSize)
      h = d(1)

      ! The indices 12+i and 12+j are the bubble DOF indices of 416 (Q3) element.
      SELECT CASE(e)
      CASE(5)
        i = 1
        j = 3
      CASE(6)
        i = 2
        j = 4        
      END SELECT

      DO k=1,2
        s = -1.0d0/3.0d0 + (k-1)*2.0d0/3.0d0
        xe = -0.25d0*(1.0d0-s)*h + 0.25d0*(1.0d0+s)*h

        f2 = (2 - 3*s + s**3)/4.0d0
        f3 = (2 + 3*s - s**3)/4.0d0
        f4 = h/8.0d0 - (h*s)/8.0d0 - (h*s**2)/8.0d0 + (h*s**3)/8.0d0
        f5 = -h/8.0d0 - (h*s)/8.0d0 + (h*s**2)/8.0d0 + (h*s**3)/8.0d0

        f = d(2)*f2+d(3)*f3+d(4)*f4+d(5)*f5
        c(1:3) = X0(1:3) + xe * ex(1:3) + f * ez(1:3)           

        SELECT CASE(k)
        CASE(1)
          BubbleNodesDelta(i,1:3) = c(1:3)
        CASE(2)
          BubbleNodesDelta(j,1:3) = c(1:3)
        END SELECT
      END DO
    END DO
    !-----------------------------------------------------------------------
    ! Now, evaluate the place of the serendipity approximation at the bubble
    ! node positions and evaluate the difference with respect to the desired 
    ! position:
    !-----------------------------------------------------------------------
    DO j=1,4
      SELECT CASE(j)
      CASE(1)
        u = -1.0d0/3.0d0
        v = -1.0d0/3.0d0
      CASE(2)
        u = 1.0d0/3.0d0
        v = -1.0d0/3.0d0
      CASE(3)
        u = 1.0d0/3.0d0
        v = 1.0d0/3.0d0
      CASE(4)
        u = -1.0d0/3.0d0
        v = 1.0d0/3.0d0
      END SELECT
      stat = BlendingSurfaceInfo(Element, Nodes, u, v, deta, a1, a2, a3, a, b, p)
      BubbleNodesDelta(j,1:3) = BubbleNodesDelta(j,1:3) - p(1:3)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE FindBubbleNodesQuad
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! This subroutine creates an orthonormal basis which gives the orientation of 
! lines of curvature at a point that is the image of the reference element 
! point (xi1,xi2) under the FE blending map. If the point is not specified,
! the point is taken to be the element centre.
!   The basis vectors may be returned via the arguments e1, e2 and e3, 
! while the coordinates of the surface point may be returned via o.
! When this subroutine is called at the element centre, it can be used to 
! define an elementwise coordinate system. In this case the optional arguments
! LagrangeNodes and TaylorParams can be used to obtain data for the Lagrange 
! interpolation to describe the shape of the projected surface on the plane 
! spanned by {e1,e2} and the coefficients of the third-order 
! Taylor polynomial which approximates the blending surface. The optional
! argument d can be used to ensure that e3 and d point to the same
! direction, while PlanarSurface indicates whether the surface is planar
! at the given point.
!   The optional argument MacroElement indicates whether a quadrilateral should
! be considered as a macro element for a subtriangulation. The argument 
! SaveProperties can be used to save the quantities computed as elementwise 
! properties to avoid later recomputation. The elementwise properties are
! as follows.
!    * 'element frame': e1, e2, e3 and o
!    * 'taylor parameters': the coefficients of the Taylor polynomial
!    * 'blending surface nodes': the nodes for the Lagrange approximation
!    * 'bubble dofs': the coefficients for bubble basis functions of Q3
!------------------------------------------------------------------------------
  SUBROUTINE LinesOfCurvatureFrame(Element, xi1, xi2, e1, e2, e3, o, TaylorParams, &
      LagrangeNodes, d, PlanarSurface, PlanarPoint, UmbilicalPoint, &
      MacroElement, SaveProperties, SizeRadiusRatio)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: xi1, xi2
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: e1(3), e2(3), e3(3)
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: o(3)
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: TaylorParams(6)
    REAL(KIND=dp), OPTIONAL, TARGET, INTENT(OUT) :: LagrangeNodes(MaxPatchNodes,3)
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: d(3)
    LOGICAL, OPTIONAL, INTENT(IN) :: PlanarSurface
    LOGICAL, OPTIONAL, INTENT(OUT) :: PlanarPoint
    LOGICAL, OPTIONAL, INTENT(OUT) :: UmbilicalPoint
    LOGICAL, OPTIONAL, INTENT(IN) :: MacroElement
    LOGICAL, OPTIONAL, INTENT(IN) :: SaveProperties  
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: SizeRadiusRatio
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: GElement => NULL()
    LOGICAL :: Stat
    LOGICAL :: Converged, ComputeTaylorPolynomial
    LOGICAL :: ApproximatePlaneDomain, CheckOrientation, WriteElementProperties
    LOGICAL :: Subtriangulation, Planar, Umbilical
    INTEGER :: Family, n, m, e, i, j, k,  GridPoint

    REAL(KIND=dp) :: u, v
    REAL(KIND=dp) :: GlobPDir1(3), GlobPDir2(3), GlobPDir3(3), X0(3)
    REAL(KIND=dp) :: Lambda1, Lambda2, LambdaMax
    REAL(KIND=dp) :: a1(3), a2(3), a3(3), a(2,2), Deta, b(2,2), ContravA(2,2)
    REAL(KIND=dp) :: c(2,2), trc, detc, PDir1(2), PDir2(2)
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

    IF ( PRESENT(TaylorParams) ) THEN
      ComputeTaylorPolynomial = .TRUE.
      TaylorParams = 0.0d0
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

    IF ( PRESENT(MacroElement) ) THEN
      Subtriangulation = MacroElement .AND. (Family == 4)
    ELSE
      Subtriangulation = .FALSE.
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
    IF (Subtriangulation) THEN
      ! To employ the macro element approach, the positions of the bubble nodes
      ! must be evaluated first ...
      CALL FindBubbleNodesQuad(Element, Nodes, BubbleNodesDelta)
      ! ... and then evaluate the blending surface data, with the bubble part
      ! taken into account:
      Stat = BlendingSurfaceInfo(Element, Nodes, u, v, deta, a1, a2, a3, a, b, X0, &
          BubbleDOFs=BubbleNodesDelta) 

      !err = sqrt(sum(BubbleNodesDelta(1,:)**2))
      !IF (err > 1.0d-10)  print *, 'diff1 = ', err
      !err = sqrt(sum(BubbleNodesDelta(2,:)**2))
      !IF (err > 1.0d-10)  print *, 'diff2 = ', err
      !err = sqrt(sum(BubbleNodesDelta(3,:)**2))
      !IF (err > 1.0d-10)  print *, 'diff3 = ', err
      !err = sqrt(sum(BubbleNodesDelta(4,:)**2))
      !IF (err > 1.0d-10)  print *, 'diff4 = ', err
    ELSE
      stat = BlendingSurfaceInfo(Element, Nodes, u, v, deta, a1, a2, a3, a, b, X0)      
    END IF

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
    detc = c(1,1)*c(2,2)-c(1,2)*c(2,1)
    trc = c(1,1) + c(2,2)
 
    !--------------------------------------------------------------
    ! Order the eigenvalues by their absolute values:
    !--------------------------------------------------------------
    IF (trc>0.0d0) THEN
      lambda2 = 0.5d0 * ( trc + SQRT(trc**2 - 4.0d0*detc) )
      lambda1 = 0.5d0 * ( trc - SQRT(trc**2 - 4.0d0*detc) )
    ELSE
      lambda2 = 0.5d0 * ( trc - SQRT(trc**2 - 4.0d0*detc) )
      lambda1 = 0.5d0 * ( trc + SQRT(trc**2 - 4.0d0*detc) )
    END IF
    !print *, 'Eigenvalues=', lambda1,lambda2
    LambdaMax = MAX(ABS(Lambda1), ABS(Lambda2))

    Planar = (ABS(Lambda1) < EPSILON(1.0)) .AND. (ABS(Lambda2) < EPSILON(1.0))

    !-----------------------------------------------------------------
    ! Another planarity check may have been done. A positive result
    ! of an erlier data test will be respected.
    !-----------------------------------------------------------------
    IF (PRESENT(PlanarSurface)) THEN
      IF (PlanarSurface .AND. .NOT. Planar) THEN
        CALL Info( 'LinesOfCurvatureFrame', &
            'Planarity checks produce different results, rejecting negative result', &
            Level=10 )
        Planar = .TRUE.
      END IF
    END IF

    Umbilical = .FALSE.
    IF (.NOT. Planar) THEN
      delta(1) = ABS(Lambda1-Lambda2)/LambdaMax
      !PRINT *, 'difference of eigenvals=',delta(1)
      Umbilical = delta(1) < UmbilicalDelta    !1000.0*EPSILON(1.0)
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
    !print *, 'T11=', BPrinc(1,1)
    !print *, 'T12=', BPrinc(1,2)
    !print *, 'T21=', BPrinc(2,1)
    !print *, 'T22=', BPrinc(2,2)

    IF (ApproximatePlaneDomain .OR. ComputeTaylorPolynomial) THEN
      IF (.NOT. ASSOCIATED(GElement)) GElement => AllocateElement()
      SELECT CASE (Family)
      CASE(3)
        GElement % Type => GetElementType(310,.FALSE.)
      CASE(4)
        GElement % Type => GetElementType(416,.FALSE.)
      END SELECT
    END IF

    IF (ApproximatePlaneDomain) THEN
      DO j=1,GElement % Type % NumberOfNodes
        IF (j > Family) THEN
          ! First we may need to map the Lagrange element coordinates to the ones of
          ! the p-element:
          IF (Family==3) THEN
            xi = -1.0d0 + 2.0d0*GElement % Type % NodeU(j) + GElement % Type % NodeV(j)
            eta = SQRT(3.0d0)*GElement % Type % NodeV(j)
          ELSE
            xi = GElement % Type % NodeU(j)
            eta = GElement % Type % NodeV(j)
          END IF

          IF (Subtriangulation) THEN
            stat = BlendingSurfaceInfo( Element, Nodes, xi, eta, deta, &
                a1, a2, a3, a, b, p, BubbleDOFs=BubbleNodesDelta)
          ELSE
            stat = BlendingSurfaceInfo( Element, Nodes, xi, eta, deta, &
                a1, a2, a3, a, b, p)         
          END IF
        ELSE
          p(1) = Nodes % x(j)
          p(2) = Nodes % y(j)
          p(3) = Nodes % z(j)
        END IF
        p = p - X0
        LagrangeNodes(j,1) = DOT_PRODUCT(p,GlobPDir1)
        LagrangeNodes(j,2) = DOT_PRODUCT(p,GlobPDir2)
        LagrangeNodes(j,3) = DOT_PRODUCT(p,GlobPDir3)
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
      ! polynomial. It is supposed that a square [-rK/4,rK/4]^2 is embedded
      ! into the plane domain obtained via the projection.
      ! TO DO: FIGURE OUT THE PRECISE SIZE OF A SQUARE THAT CAN BE EMBEDDED 
      ! --------------------------------------------------------------------
      hk = rK/2.0d0               ! The width of stencil
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
            IF (Subtriangulation) THEN
              stat = BlendingSurfaceInfo( Element, Nodes, uk, vk, &
                  deta, a1, a2, a3, a, b, pk, BubbleDOFs=BubbleNodesDelta)
            ELSE
              stat = BlendingSurfaceInfo( Element, Nodes, uk, vk, &
                  deta, a1, a2, a3, a, b, pk)             
            END IF
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

      IF (Planar) THEN
        IF ( (ABS(APar) > 1000.0*EPSILON(1.0)) .AND. (ABS(BPar) > 1000.0*EPSILON(1.0)) ) THEN
          CALL Warn('LinesOfCurvatureFrame', 'Possibly inaccurate Taylor polynomial (planar point)')
          print *, 'Apar,Lambda1=', Apar, Lambda1
          print *, 'Bpar,Lambda2=', Bpar, Lambda2
        END IF
      ELSE
        IF (Umbilical) THEN
          err = ABS(APar-BPar)/MAX(ABS(APar),ABS(BPar))
          IF ( err > 5.0d0*UmbilicalDelta ) THEN
            CALL Warn('LinesOfCurvatureFrame', 'Possibly inaccurate Taylor polynomial (umbilical point)')
            print *, '|APar-BPar|/max(|APar|,|BPar|)=', err
          END IF
        ELSE
          !-------------------------------------------------------------------------------------
          ! Make a rough check that a reparametrization of the surface based on the Taylor
          ! polynomial will be plausible:
          !-------------------------------------------------------------------------------------
          IF ( ABS(FPar/(BPar - APar))*rK > 1.0d0 .OR. ABS(EPar/(BPar - APar))*rK > 1.0d0 ) &
              CALL Warn('LinesOfCurvatureFrame', 'Possibly very rough geometry model used')

          err = ABS( ABS(BPar)-ABS(Lambda2) ) / LambdaMax
          err = MAX( err, ABS(ABS(APar)-ABS(lambda1))/LambdaMax )
          IF ( err > 5.0d-2) THEN
            CALL Warn('LinesOfCurvatureFrame', 'Possibly inaccurate Taylor polynomial')
            print *, '|Apar-Lambda1|/LambdaMax=', ABS(ABS(APar)-ABS(lambda1))/LambdaMax
            print *, '|Bpar-Lambda2|/LambdaMax=', ABS(ABS(BPar)-ABS(lambda2))/LambdaMax
          END IF
        END IF
      END IF

    END IF

    ! -------------------------------------------------------------------------
    ! Indicate whether the surface at the specified point 
    ! can be approximated by a plane up to errors O(h^3).
    ! -------------------------------------------------------------------------
    IF ( PRESENT(PlanarPoint) ) PlanarPoint = Planar
    IF ( PRESENT(UmbilicalPoint) ) UmbilicalPoint = Umbilical

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
      FrameData(1:3) = GlobPDir1(1:3)
      FrameData(4:6) = GlobPDir2(1:3)
      FrameData(7:9) = GlobPDir3(1:3)
      FrameData(10:12) = X0(1:3)
      CALL SetElementProperty('element frame', FrameData(1:FrameDataSize), Element) 

      IF (ComputeTaylorPolynomial) CALL SetElementProperty('taylor parameters', &
          TaylorParams(1:6), Element)

      IF (ApproximatePlaneDomain) THEN
        NodesArray(1:MaxPatchNodes) = LagrangeNodes(1:MaxPatchNodes,1)
        NodesArray(MaxPatchNodes+1:2*MaxPatchNodes) = LagrangeNodes(1:MaxPatchNodes,2)
        NodesArray(2*MaxPatchNodes+1:3*MaxPatchNodes) = LagrangeNodes(1:MaxPatchNodes,3)
        CALL SetElementProperty('blending surface nodes', NodesArray(1:3*MaxPatchNodes), Element)
      END IF

      IF (Subtriangulation) THEN
        NodesArray(1:4) = BubbleNodesDelta(1:4,1)
        NodesArray(5:8) = BubbleNodesDelta(1:4,2)
        NodesArray(9:12) = BubbleNodesDelta(1:4,3)
        CALL SetElementProperty('bubble dofs', NodesArray(1:12), Element)
      END IF
    END IF

    !print *, 'o=', o
    !print *, 'e1=', e1
    !print *, 'e2=', e2
    !print *, 'e3=', e3
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
! LocalFrameNodes defines x, the points y are returned via the argument PatchNodes 
! and the form of g is defined by the parameters TaylorParams. If ZNodes is
! supplied, approximations of nodal z-coordinates is computed using a third-order
! Taylor polynomial in the coordinates y of the final patch K. This option
! can be used to cross-check different approximations but may not have final utility.
!-----------------------------------------------------------------------------------
  SUBROUTINE LinesOfCurvaturePatch(LocalFrameNodes, TaylorParams, PatchNodes, &
      Family, ZNodes, PlanarSurface, UmbilicalPoint)
!-----------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), TARGET, INTENT(IN) :: LocalFrameNodes(MaxPatchNodes,2)
    REAL(KIND=dp), INTENT(IN) :: TaylorParams(6)
    REAL(KIND=dp), INTENT(OUT) :: PatchNodes(MaxPatchNodes,2)
    INTEGER, INTENT(IN) :: Family
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: ZNodes(MaxPatchNodes)
    LOGICAL, OPTIONAL, INTENT(IN) :: PlanarSurface
    LOGICAL, OPTIONAL, INTENT(IN) :: UmbilicalPoint
!-----------------------------------------------------------------------------------
    LOGICAL :: Planar, Umbilical, Converged
    INTEGER :: i, k, n
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
        PatchNodes = LocalFrameNodes
        IF (PRESENT(ZNodes)) ZNodes = 0.0d0       
        RETURN
      END IF
    END IF

    A = TaylorParams(1)
    B = TaylorParams(2)
    D = TaylorParams(3)
    E = TaylorParams(4)
    F = TaylorParams(5)
    G = TaylorParams(6)

    IF (PRESENT(UmbilicalPoint)) THEN
      Umbilical = UmbilicalPoint
    ELSE
      Umbilical =  ABS(B-A)/MAX(ABS(A),ABS(B)) < UmbilicalDelta
    END IF

    IF (.NOT. Umbilical) THEN
      !---------------------------------------------------------------------------------
      ! The following constants relate to the most general third-order polynomial 
      ! perturbations of the coordinate functions; cf. the definitions of x and y
      ! below.
      !---------------------------------------------------------------------------------
      c1 = F/(B - A)                 ! To make the 2nd fundamental form diagonal
      c2 = E/(B - A)                 ! To make the 2nd fundamental form diagonal
      c4 = -A**2 - c2**2 + c1**2     ! To make the Christoffel symbol C111 constant
      c5 = -2.0d0*c1*c2              ! To make the Christoffel symbol C111 constant
      c6 = c2**2 - c1**2             ! To make the Christoffel symbol C112 constant
      c7 = -c5                       ! To make the Christoffel symbol C222 constant
      b7 = -B**2 + c2**2 - c1**2     ! To make the Christoffel symbol C222 constant
      b5 = -A*B - c6                 ! Orthogonal surface basis vectors up to O(h^3)
    ELSE
      !---------------------------------------------------------------------------------
      ! If the point is umbilical and the surface is not planar, then the surface
      ! is considered to be umbilical and hence spherical. Surface approximation
      ! is then constructed such that its Taylor expansion agrees with the Taylor's 
      ! polynomial of a sphere.
      !---------------------------------------------------------------------------------
      c1 = 0.0d0                 ! To make the 2nd fundamental form diagonal
      c2 = 0.0d0                 ! To make the 2nd fundamental form diagonal
      c4 = -A**2                 ! To make the Christoffel symbol C111 constant
      c5 = 0.0d0                 ! To make the Christoffel symbol C111 constant
      c7 = 0.0d0                 ! To make the Christoffel symbol C222 constant
      b7 = -B**2                 ! To make the Christoffel symbol C222 constant
      c6 = -A*B/2.0d0            ! Symmetry of coordinate perturbations
      b5 = -A*B - c6             ! Orthogonal surface basis vectors up to O(h^3)
    END IF

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
      UmbilicalPoint)
!-----------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: x, y
    REAL(KIND=dp), INTENT(IN) :: TaylorParams(6)
    REAL(KIND=dp), INTENT(IN) :: e1(3), e2(3), e3(3), o(3)
    REAL(KIND=dp), INTENT(OUT) :: a1(3), a2(3), a3(3)
    REAL(KIND=dp), INTENT(OUT) :: A11, A22, SqrtDetA, B11, B22
    REAL(KIND=dp), INTENT(OUT) :: C111, C112, C221, C222, C211, C212
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: dual1(3), dual2(3)
    REAL(KIND=dp), OPTIONAL, INTENT(OUT) :: XGlob, YGlob, ZGlob
    LOGICAL, OPTIONAL, INTENT(IN) :: LowestOrderBasis
    LOGICAL, OPTIONAL, INTENT(IN) :: PlanarPoint
    LOGICAL, OPTIONAL, INTENT(IN) :: UmbilicalPoint
!-----------------------------------------------------------------------------------
    LOGICAL :: ReturnDualBasis, ReturnGlobalCoords
    LOGICAL :: Planar, Umbilical
    REAL(KIND=dp) :: A, B, D, E, F, G
    REAL(KIND=dp) :: c1, c2, c4, c5, c6, c7, b5, b7
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

    IF (PRESENT(UmbilicalPoint)) THEN
      Umbilical = UmbilicalPoint
    ELSE
      Umbilical =  ABS(B-A)/MAX(ABS(A),ABS(B)) < UmbilicalDelta
    END IF

    ! ---------------------------------------------------------------------------- 
    ! The constants to diagonalize the fundamental forms
    ! ---------------------------------------------------------------------------- 
    IF (Umbilical) THEN
      c1 = 0.0d0                 ! To make the 2nd fundamental form diagonal
      c2 = 0.0d0                 ! To make the 2nd fundamental form diagonal
      c4 = -A**2                 ! To make the Christoffel symbol C111 constant
      c5 = 0.0d0                 ! To make the Christoffel symbol C111 constant
      c7 = 0.0d0                 ! To make the Christoffel symbol C222 constant
      b7 = -B**2                 ! To make the Christoffel symbol C222 constant
      c6 = -A*B/2.0d0            ! Symmetry of coordinate perturbations
      b5 = -A*B - c6             ! Orthogonal surface basis vectors up to O(h^3)
    ELSE
      c1 = F/(B - A)             ! To make the 2nd fundamental form diagonal
      c2 = E/(B - A)             ! To make the 2nd fundamental form diagonal
      c4 = -A**2 - c2**2 + c1**2 ! To make the Christoffel symbol C111 constant
      c5 = -2.0d0*c1*c2          ! To make the Christoffel symbol C111 constant
      c6 = c2**2 - c1**2         ! To make the Christoffel symbol C112 constant
      c7 = -c5                   ! To make the Christoffel symbol C222 constant
      b7 = -B**2 + c2**2 - c1**2 ! To make the Christoffel symbol C222 constant
      b5 = -A*B - c6             ! Orthogonal surface basis vectors up to O(h^3)
    END IF

    ! ---------------------------------------------------------------------------- 
    ! The covariant surface basis vectors expanded up to O(h^3)
    ! ---------------------------------------------------------------------------- 
    IF (Umbilical) THEN
      a1(:) = (0.1D1 - A ** 2 * x ** 2 / 0.2D1 - A * B * y ** 2 / 0.4D1) * e1(:) + &
          (-A * B * x * y / 0.2D1) * e2(:) + (A * x) * e3(:)
      a2(:) = (-A * B * x * y / 0.2D1) * e1(:) + &
          (0.1D1 - A * B * x ** 2 / 0.4D1 - B ** 2 * y ** 2 / 0.2D1) * e2(:) + &
          (B * y) * e3(:)
    ELSE
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
    END IF

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
    IF (Umbilical) THEN
      C111 = 0.0d0
      C112 = A * B * y / 0.2D1
      C221 = A * B * x / 0.2D1
      C222 = 0.0d0
      C211 = -A * B * y / 0.2D1
      C212 = -A * B * x / 0.2D1
    ELSE
      C111 = -c1
      C112 = -c2
      C221 = c1 + A * B  * x
      C222 = c2
      C211 = c2
      C212 = -c1 - A * B * x
    END IF
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
! inaccurate results for thin shells (with low p)! 
!------------------------------------------------------------------------------
  SUBROUTINE ShellLocalMatrix( BGElement, n, nd, TaylorParams, PatchNodes, &
      e1, e2, e3, o, m, PlateBody, SphericalSurface, StrainReductionMethod, &
      Bubbles, Area, Error)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: BGElement  ! An element of background mesh
    INTEGER, INTENT(IN) :: n                           ! The number of background element nodes
    INTEGER, INTENT(IN) :: nd                          ! The number of DOFs per component (after
                                                       ! static condensation if bubbles are used)
    REAL(KIND=dp), INTENT(IN) :: TaylorParams(6)       ! The coefficients of the Taylor approximation
    REAL(KIND=dp), INTENT(IN) :: PatchNodes(MaxPatchNodes,2) ! The nodes of principal coordinate patch
    REAL(KIND=dp), INTENT(IN) :: e1(3), e2(3), e3(3)   ! The basis of the local frame
    REAL(KIND=dp), INTENT(IN) :: o(3)                  ! The origin of the local frame
    INTEGER, INTENT(IN) :: m                           ! The number of DOFs per node
    LOGICAL, INTENT(IN) :: PlateBody                   ! To indicate that the surface is flat
    LOGICAL, INTENT(IN) :: SphericalSurface            ! To indicate that the surface is considered to spherical
    INTEGER, INTENT(IN) :: StrainReductionMethod       ! The choice of strain reduction method
    LOGICAL, INTENT(IN) :: Bubbles                     ! To indicate that bubble functions are used
    REAL(KIND=dp), INTENT(INOUT) :: Area               ! A variable for area compution
    REAL(KIND=dp), INTENT(INOUT) :: Error              ! A variable for error compution
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element => NULL()
    TYPE(Element_t), POINTER :: GElement => NULL()     
    TYPE(Nodes_t) :: Nodes, PNodes, PRefNodes
    TYPE(ValueList_t), POINTER :: BodyForce
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: Stat, Found
    LOGICAL :: SuperParametric, NoNormalStress, SecondOrder, PVersion
    LOGICAL :: UseShearCorrection, UseBubbles

    INTEGER :: Family, ReducedStrainDim, ReductionMethod
    INTEGER :: ShearReductionMethod, StretchReductionMethod
    INTEGER :: DOFs, nb, i, j, k, p, t, i0, j0, csize, GElementNodes 

    ! Prepare for the scenario that one elementwise bubble per component can be added:
    REAL(KIND=dp) :: Stiff(m*nd+m,m*nd+m), Mass(m*nd+m,m*nd+m), Force(m*nd+m)
    REAL(KIND=dp) :: Damp(m*nd+m,m*nd+m)
    REAL(KIND=dp) :: BM(4,m*nd+m), BS(2,m*nd+m), BB(3,m*nd+m)
    REAL(KIND=dp) :: Basis(nd+1), dBasis(nd+1,3)

    REAL(KIND=dp) :: QBlock(3,3), Q(m*nd,m*nd)
    REAL(KIND=dp) :: CMat(4,4), GMat(2,2)
    REAL(KIND=dp) :: A11, A22, SqrtDetA, A1, A2, B11, B22
    REAL(KIND=dp) :: C111, C112, C221, C222, C211, C212
    REAL(KIND=dp) :: abasis1(3), abasis2(3), abasis3(3)
    REAL(KIND=dp) :: yk, y1, y2, XGlob, YGlob, ZGlob
    REAL(KIND=dp) :: uq, vq, sq
    REAL(KIND=dp) :: BGBasis(n), GBasis(MaxPatchNodes), PBasis(nd)
    REAL(KIND=dp) :: StrainBasis(4,3)             ! Four rows large enough for p=1
    REAL(KIND=dp) :: ReductionDOFsArray(4,2*nd+2) ! Four rows large enough for p=1 with one bubble
    REAL(KIND=dp) :: DOFsTransform(2,3)
    REAL(KIND=dp) :: ShearParMat(2,2,nd), StretchParMat1(2,2,nd)
    REAL(KIND=dp) :: StretchParMat2(2,2,nd), StretchParMat3(2,2,nd)
    REAL(KIND=dp) :: StretchParMat4(2,2,nd), StretchParMat5(2,2,nd), StretchParMat6(2,2,nd)
    REAL(KIND=dp) :: PoissonRatio(n), YoungsMod(n), ShellThickness(n), Load(n), rho(n), rho0
    REAL(KIND=dp) :: nu, E, h, NormalTraction, Kappa
    REAL(KIND=dp) :: DetJ, Weight, Norm

    SAVE Element, GElement, Nodes, PNodes, PRefNodes
!------------------------------------------------------------------------------
    IF (m /= 6) CALL Fatal('ShellLocalMatrix', 'Wrong number of unknown fields')

    Family = GetElementFamily(BGElement)
    PVersion = IsPElement(BGElement)

    UseBubbles = .FALSE.

    IF (Pversion) THEN
      ReductionMethod = NoStrainReduction
    ELSE
      IF (StrainReductionMethod == AutomatedChoice) THEN
        ! ------------------------------------------------------------------------------
        ! This option saves the user from deciding a strain reduction method
        ! ------------------------------------------------------------------------------
        SELECT CASE(Family)
        CASE(3)
          ReductionMethod = DoubleReduction
          UseBubbles = .TRUE.
        CASE(4)
          ReductionMethod = CurlKernel
        END SELECT
      ELSE
        ReductionMethod = StrainReductionMethod
        ! Currently bubbles can be activated only for certain strategies:
        IF (ReductionMethod == NoStrainReduction .OR. &
            ReductionMethod == DoubleReduction) THEN
          UseBubbles = Bubbles
        END IF
      END IF
    END IF

    ! ------------------------------------------------------------------------------
    ! In the case of triangles asking for the kernel version switches to the
    ! compound strain reduction strategy:
    ! ------------------------------------------------------------------------------
    IF (Family == 3 .AND. ReductionMethod == CurlKernel) &
        ReductionMethod = DoubleReduction

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
        UseShearCorrection = .FALSE.
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
        ! Coefficients to transform from MITC DOFs to the kernel DOFs: 
        DOFsTransform(1,1) = 1.0d0/3.0d0
        DOFsTransform(1,2) = -1.0d0/6.0d0
        DOFsTransform(1,3) = DOFsTransform(1,2)
        DOFsTransform(2,1) = 0.0d0
        DOFsTransform(2,2) = 1.0d0/(2.0d0*sqrt(3.0d0))
        DOFsTransform(2,3) = -1.0d0/(2.0d0*sqrt(3.0d0))
      ELSE
        CALL Fatal('ShellLocalMatrix', 'Double reduction is not defined for quads')
      END IF
    END SELECT

    ! ------------------------------------------------------------------------------
    ! The size of the constitutive matrix for 2D shell equations
    ! ------------------------------------------------------------------------------
    csize = 4

    ! ------------------------------------------------------------------------------
    ! The DOFs count: Currently, FE bubbles can be employed only in a very special 
    ! way by augmenting the two rotation components
    ! ------------------------------------------------------------------------------
    IF (UseBubbles) THEN
      nb = 2
      !IF (ReductionMethod == DoubleReduction) nb = nb + 1 ! Try with u3-bubbles
    ELSE
      nb = 0
    END IF
    DOFs = m*nd + nb

    ! ------------------------------------------------------------------------------
    ! The following parameter defines how the normal strain variable is handled.
    ! With NoNormalStress = .TRUE. the traditional condition of vanishing normal
    ! stress is enforced a priori. Otherwise the energy principle is used to determine
    ! the normal strain variable.
    ! ------------------------------------------------------------------------------
    NoNormalStress = .FALSE. ! The choice .TRUE. does not work anymore if membrabe
                             ! strain reductions are active! 

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
    SuperParametric = .FALSE. ! This relates to an experimental version which is not active 

    ! --------------------------------------------------------------------------
    ! Allocate a Lagrange interpolation element structure corresponding to the
    ! "Element" keyword. Note that the implementation of the p-version is not
    ! fully optimal since loads and BCs may be interpolated only by using the node 
    ! coordinates of the background element.
    ! --------------------------------------------------------------------------

    IF ( PVersion ) THEN
      IF ( .NOT. ASSOCIATED(Element) ) Element => AllocateElement()
      SELECT CASE(nd)
      CASE(3)
        Element % Type => GetElementType(303, .FALSE.)
      CASE(4)
        IF (Family /= 4) CALL Fatal('ShellLocalMatrix', 'Unsupported p-element')
        Element % Type => GetElementType(404, .FALSE.)
      CASE(8)
        Element % Type => GetElementType(408, .FALSE.)
      CASE(9)
        Element % Type => GetElementType(409, .FALSE.)
      CASE DEFAULT
        CALL Fatal('ShellLocalMatrix', 'Unsupported p-element type')
      END SELECT
    ELSE
      Element => BGElement
    END IF

    ! --------------------------------------------------------------------------
    ! The element structure corresponding to the surface reconstruction:
    ! --------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED(GElement) ) GElement => AllocateElement()
    SELECT CASE(Family)
    CASE(3)
      GElement % Type => GetElementType(310, .FALSE.)
    CASE(4)
      GElement % Type => GetElementType(416, .FALSE.)
    CASE DEFAULT
      CALL Fatal('ShellLocalMatrix', 'Unsupported (geometry model) element type')
    END SELECT

    ! --------------------------------------------------------------------------
    ! Create a node variable suitable for defining the isoparametric element 
    ! map from the reference element onto the set which is the domain of lines of
    ! curvature coordinates, i.e. use as many nodes as DOFs in the spartial 
    ! discretization. Then update the coordinate values of the node variable. 
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

    GElementNodes = GElement % Type % NumberOfNodes
    DO j=1,nd
      u = Element % Type % NodeU(j)
      v = Element % Type % NodeV(j)
      CALL NodalBasisFunctions2D(GBasis, GElement, u, v)
      Nodes % x(j) = SUM(GBasis(1:GElementNodes) * PatchNodes(1:GElementNodes,1))
      Nodes % y(j) = SUM(GBasis(1:GElementNodes) * PatchNodes(1:GElementNodes,2))
      Nodes % z(j) = 0.0d0
    END DO


    ! --------------------------------------------------------------------------
    ! If p-basis is used for approximating the shell equations,  
    ! create also the isoparametric geometry representation in terms of p-basis:
    ! --------------------------------------------------------------------------    
    IF ( PVersion ) THEN
      ! --------------------------------------------------------------------------
      ! Solve coefficients to express the isoparametric geometry approximation by 
      ! an equivalent expansion in terms of the p-basis:
      ! --------------------------------------------------------------------------     
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
      
      ! Initialize to avoid errors from the next ElementInfo call:
      PNodes % x(1:nd) = 0.0d0
      PNodes % y(1:nd) = 0.0d0
      PNodes % z(1:nd) = 0.0d0
      PNodes % x(1:n) = Nodes % x(1:n)
      PNodes % y(1:n) = Nodes % y(1:n)

      DO k=1,2
        STIFF = 0.0d0
        FORCE = 0.0d0

        IP = GaussPoints( BGElement )
        DO t=1,IP % n 
          ! ------------------------------------------------------------------
          ! Don't change the order of ElementInfo calls to get the right detJ
          ! ------------------------------------------------------------------
          stat = ElementInfo( BGElement, PNodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, PBasis )

          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )         

          DO i=1,nd
            DO j=1,nd
              STIFF(i,j) = STIFF(i,j) + IP % s(t) * DetJ * PBasis(i) * PBasis(j)
            END DO
          END DO

          SELECT CASE(k)
          CASE(1)
            yk = SUM( Nodes % x(1:nd) * Basis(1:nd) )
          CASE(2)
            yk = SUM( Nodes % y(1:nd) * Basis(1:nd) )           
          END SELECT
          FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * yk * PBasis(1:nd)
        END DO

        CALL LUSolve(nd, Stiff(1:nd,1:nd), Force(1:nd))
        SELECT CASE(k)
        CASE(1)
          PNodes % x(1:nd) = Force(1:nd)
        CASE(2)
          PNodes % y(1:nd) = Force(1:nd) 
        END SELECT
      END DO

      ! ------------------------------------------------------------------------
      ! An additional node variable needed for superparametric case:
      ! ------------------------------------------------------------------------
      IF (SuperParametric) THEN
        IF ( .NOT. ASSOCIATED( PRefNodes % x ) ) THEN
          ALLOCATE( PRefNodes % x(n), PRefNodes % y(n), PRefNodes % z(n) )
        ELSE
          IF (n > SIZE(PRefNodes % x)) THEN
            DEALLOCATE(PRefNodes % x, PRefNodes % y, PRefNodes % z)
            ALLOCATE( PRefNodes % x(n), PRefNodes % y(n), PRefNodes % z(n) ) 
          END IF
        END IF
        PRefNodes % NumberOfNodes = n
        PrefNodes % x(1:n) = BGElement % Type % NodeU(1:n)
        PrefNodes % y(1:n) = BGElement % Type % NodeV(1:n)
        PrefNodes % z(1:n) = BGElement % Type % NodeW(1:n)
      END IF
    END IF


    ! --------------------------------------------------------------------------
    ! Body forces, material parameters and the shell thickness:
    ! --------------------------------------------------------------------------
    PoissonRatio(1:n) = GetReal(GetMaterial(), 'Poisson Ratio')
    YoungsMod(1:n) = GetReal(GetMaterial(), 'Youngs Modulus')
    ShellThickness(1:n) = GetReal(GetMaterial(), 'Shell Thickness')
    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
        Load(1:n) = GetReal(BodyForce, 'Normal Pressure', Found)
    IF (TransientSimulation) &
        rho(1:n) = GetReal(GetMaterial(), 'Density')


    ! --------------------------------------------------------------------------
    ! Strain reduction operators will be applied to fields Cu where
    ! C is a 2X2-matrix of shell model parameters and u a 2-component field.
    ! In the following we create descriptions of all necessary matrices C so
    ! that the components of C can be evaluated via the Lagrange interpolation.
    !
    ! TO DO: Should we create nodal vectors of all model parameters to avoid
    ! a later call of SurfaceBasisVectors?
    ! -------------------------------------------------------------------------- 
    ShearParMat = 0.0d0
    StretchParMat1 = 0.0d0
    StretchParMat2 = 0.0d0
    StretchParMat3 = 0.0d0
    StretchParMat4 = 0.0d0
    StretchParMat5 = 0.0d0
    StretchParMat6 = 0.0d0
    DO j=1,nd
      CALL SurfaceBasisVectors(Nodes % x(j), Nodes % y(j), TaylorParams, e1, e2, e3, o, &
          abasis1, abasis2, abasis3, A11, A22, SqrtDetA, B11, B22, C111, C112, C221, C222, &
          C211, C212, XGlob=XGlob, YGlob=YGlob, ZGlob=ZGlob, PlanarPoint=PlateBody, &
          UmbilicalPoint=SphericalSurface)
 
      ShearParMat(1,1,j) = B11/a11
      ShearParMat(2,2,j) = B22/a22

      StretchParMat1(1,1,j) = C111
      StretchParMat1(2,2,j) = C212

      StretchParMat2(1,1,j) = C112
      StretchParMat2(2,2,j) = C211

      StretchParMat3(1,1,j) = C211
      StretchParMat3(2,2,j) = C222

      StretchParMat4(1,1,j) = C212
      StretchParMat4(2,2,j) = C221

      StretchParMat5(1,1,j) = B11
      StretchParMat5(2,2,j) = 0.0d0

      StretchParMat6(1,1,j) = 0.0d0
      StretchParMat6(2,2,j) = B22

    END DO

    IF (PVersion) THEN
      IP = GaussPoints( BGElement )
    ELSE
      IF (SuperParametric) THEN
        IP = GaussPoints( Element,  Element % TYPE % GaussPoints2 )
      ELSE
        IP = GaussPoints( Element )
      END IF
    END IF

    Mass = 0.0d0
    Damp = 0.0d0
    Stiff = 0.0d0
    Force = 0.0d0

    QUADRATURELOOP: DO t=1,IP % n

      BM = 0.0d0
      BB = 0.0d0
      BS = 0.0d0

      IF ( PVersion ) THEN

!        IF (SuperParametric) THEN
!          ! Get p-basis on the reference element ...
!          stat = ElementInfo( BGElement, PRefNodes, IP % U(t), IP % V(t), &
!              IP % W(t), detJ, Basis, dBasis )
!          ! ... and get the derivatives with respect to lines of curvature coordinates by just transforming the 
!          ! derivatives of basis functions taken with respect to the reference element coordinates:
!          stat = SuperParametricElementInfo( BGElement, GElement, GBasis, PatchNodes(1:16,1), &
!              PatchNodes(1:16,2), IP % U(t), IP % V(t), detJ, Basis, dBasis, ReadyBasis = .TRUE.)
!          y1 = SUM( PatchNodes(1:16,1) * GBasis(1:16) )
!          y2 = SUM( PatchNodes(1:16,2) * GBasis(1:16) )            
!
!        ELSE

        ! Here the element mapping is isoparametric:

        uq = IP % U(t)
        vq = IP % V(t)
        sq = IP % S(t)          

        stat = ElementInfo( BGElement, PNodes, uq, vq, IP % W(t), detJ, Basis, dBasis )
        y1 = SUM( PNodes % x(1:nd) * Basis(1:nd) )
        y2 = SUM( PNodes % y(1:nd) * Basis(1:nd) )

!        END IF

      ELSE

!        IF (SuperParametric) THEN
!          stat = SuperParametricElementInfo( Element, GElement, GBasis, PatchNodes(1:16,1), &
!              PatchNodes(1:16,2), IP % U(t), IP % V(t), detJ, Basis, dBasis )
!          y1 = SUM( PatchNodes(1:16,1) * GBasis(1:16) )
!          y2 = SUM( PatchNodes(1:16,2) * GBasis(1:16) )  
!        ELSE


        ! Use isoparametric element map:
        ! ReductionOperatorInfo should give all necessary basis functions without the
        ! standard ElementInfo call as
        !
        !  stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
        !      detJ, Basis, dBasis )
        !
        ! Now ReductionOperatorInfo assumes p-reference elements, so
        ! switch to the reference p-element (perhaps we could use the standard reference 
        ! elements, since continuity across elements is not needed for the reduced strain).
        ! ---------------------------------------------------------------------------------
        IF (Family==3 .AND. .NOT. PVersion) THEN
          uq = -1.0d0 + 2.0d0 * IP % U(t) + IP % V(t)
          vq = SQRT(3.0d0) * IP % V(t)
          sq = SQRT(3.0d0) * 2.0d0 * IP % S(t)
        ELSE
          uq = IP % U(t)
          vq = IP % V(t)
          sq = IP % S(t)
        END IF

        stat = ReductionOperatorInfo( Element, Nodes, uq, vq, StrainBasis, &
            ReductionMethod, ApplyPiolaTransform = .TRUE., detF=detJ, Basis=Basis, &
            dBasis=dBasis, Bubbles=UseBubbles )

        y1 = SUM( Nodes % x(1:nd) * Basis(1:nd) )
        y2 = SUM( Nodes % y(1:nd) * Basis(1:nd) )

!        END IF

      END IF

      ! ------------------------------------------------------------------------------
      ! The fundamental forms and Christoffel symbols at the point (y1,y2) of the
      ! principal coordinate patch:
      ! ------------------------------------------------------------------------------
      CALL SurfaceBasisVectors(y1, y2, TaylorParams, e1, e2, e3, o, abasis1, &
          abasis2, abasis3, A11, A22, SqrtDetA, B11, B22, C111, C112, C221, C222, &
          C211, C212, XGlob=XGlob, YGlob=YGlob, ZGlob=ZGlob, PlanarPoint=PlateBody, &
          UmbilicalPoint=SphericalSurface)

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
      IF (TransientSimulation) rho0 = SUM( rho(1:n) * Basis(1:n) )

      ! The matrix description of the elasticity tensor:
      CALL ElasticityMatrix(CMat, GMat, A1, A2, E, nu)

      ! Shear correction factor:
      IF (UseShearCorrection) THEN
        CALL ShearCorrectionFactor(Kappa, h, Nodes % x(1:n), Nodes % y(1:n), n)
      ELSE
        Kappa = 1.0d0
      END IF

      !---------------------------------------------------------------
      ! The part corresponding to the membrane strains:
      !---------------------------------------------------------------
      Weight = h * SqrtDetA * detJ * sq

      IF ( (ReductionMethod /= NoStrainReduction) .AND. (.NOT. PlateBody) ) THEN

        !------------------------------------------------------------------------------------------
        ! Apply strain reduction to vector (C111 u1, C212 u2):
        !------------------------------------------------------------------------------------------
        IF (ReductionMethod == DoubleReduction) THEN
          ! First, get DOFs for the RT interpolant ...
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              3, nd, MITC, StretchParMat1)
          ! and then apply a second round of reductions:
          ReductionDOFsArray(1:2,1:2*nd) = MATMUL(DOFsTransform(1:2,1:3), ReductionDOFsArray(1:3,1:2*nd))
        ELSE
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              ReducedStrainDim, nd, ReductionMethod, StretchParMat1)
        END IF
        DO p=1,nd
          DO j=1,ReducedStrainDim
            BM(1,(p-1)*m+1) = BM(1,(p-1)*m+1) - ReductionDOFsArray(j,2*p-1) * StrainBasis(j,1)
            BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) - ReductionDOFsArray(j,2*p-1) * StrainBasis(j,2)

            BM(1,(p-1)*m+2) = BM(1,(p-1)*m+2) - ReductionDOFsArray(j,2*p) * StrainBasis(j,1)
            BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) - ReductionDOFsArray(j,2*p) * StrainBasis(j,2)

            IF (.NOT. NoNormalStress) THEN
              BM(4,(p-1)*m+1) = BM(4,(p-1)*m+1) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*p-1) * StrainBasis(j,1)/A1**2)
              BM(4,(p-1)*m+2) = BM(4,(p-1)*m+2) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*p) * StrainBasis(j,1)/A1**2)
            END IF
          END DO
        END DO

        ! Apply strain reduction to vector (C112 u2, C211 u1). Notice the order of components.
        IF (ReductionMethod == DoubleReduction) THEN
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              3, nd, MITC, StretchParMat2)
          ReductionDOFsArray(1:2,1:2*nd) = MATMUL(DOFsTransform(1:2,1:3), ReductionDOFsArray(1:3,1:2*nd))
        ELSE
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              ReducedStrainDim, nd, ReductionMethod, StretchParMat2)
        END IF
        DO p=1,nd
          DO j=1,ReducedStrainDim
            BM(1,(p-1)*m+2) = BM(1,(p-1)*m+2) - ReductionDOFsArray(j,2*p-1) * StrainBasis(j,1)
            BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) - ReductionDOFsArray(j,2*p-1) * StrainBasis(j,2)

            BM(1,(p-1)*m+1) = BM(1,(p-1)*m+1) - ReductionDOFsArray(j,2*p) * StrainBasis(j,1)
            BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) - ReductionDOFsArray(j,2*p) * StrainBasis(j,2)

            IF (.NOT. NoNormalStress) THEN
              BM(4,(p-1)*m+1) = BM(4,(p-1)*m+1) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*p) * StrainBasis(j,1)/A1**2)
              BM(4,(p-1)*m+2) = BM(4,(p-1)*m+2) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*p-1) * StrainBasis(j,1)/A1**2)
            END IF
          END DO
        END DO

        ! Apply strain reduction to vector (B11 u3, 0):
        IF (ReductionMethod == DoubleReduction) THEN
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              3, nd, MITC, StretchParMat5)
          ReductionDOFsArray(1:2,1:2*nd) = MATMUL(DOFsTransform(1:2,1:3), ReductionDOFsArray(1:3,1:2*nd))
        ELSE
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              ReducedStrainDim, nd, ReductionMethod, StretchParMat5)
        END IF
        DO p=1,nd
          DO j=1,ReducedStrainDim
            BM(1,(p-1)*m+3) = BM(1,(p-1)*m+3) - ReductionDOFsArray(j,2*p-1) * StrainBasis(j,1)
            BM(3,(p-1)*m+3) = BM(3,(p-1)*m+3) - ReductionDOFsArray(j,2*p-1) * StrainBasis(j,2)

            IF (.NOT. NoNormalStress) THEN
              BM(4,(p-1)*m+3) = BM(4,(p-1)*m+3) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*p-1) * StrainBasis(j,1)/A1**2)
            END IF
          END DO
        END DO

        ! If normal displacement is augmented with a bubble, we need to compute an augmented DOFs array.
        ! This is commented out as this gives poor results.
        !IF (Bubbles .AND. ReductionMethod == DoubleReduction) THEN
        !  CALL ReductionOperatorBubbleDofs(Element, Nodes, ReductionDOFsArray(1:2,2*nd+1:2*nd+2), &
        !      ReducedStrainDim, 1, Family, CurlKernel, StretchParMat5)
        !
        !  DO j=1,ReducedStrainDim
        !    BM(1,nd*m+3) = BM(1,nd*m+3) - ReductionDOFsArray(j,2*nd+1) * StrainBasis(j,1)
        !    BM(3,nd*m+3) = BM(3,nd*m+3) - ReductionDOFsArray(j,2*nd+1) * StrainBasis(j,2)
        !    IF (.NOT. NoNormalStress) THEN
        !      BM(4,nd*m+3) = BM(4,nd*m+3) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*nd+1) * StrainBasis(j,1)/A1**2)
        !    END IF
        !  END DO
        !END IF       

        ! Apply strain reduction to vector (C211 u1, C222 u2):
        IF (ReductionMethod == DoubleReduction) THEN
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              3, nd, MITC, StretchParMat3)
          ReductionDOFsArray(1:2,1:2*nd) = MATMUL(DOFsTransform(1:2,1:3), ReductionDOFsArray(1:3,1:2*nd))
        ELSE
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              ReducedStrainDim, nd, ReductionMethod, StretchParMat3)
        END IF
        DO p=1,nd
          DO j=1,ReducedStrainDim
            BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) - ReductionDOFsArray(j,2*p-1) * StrainBasis(j,1)
            BM(2,(p-1)*m+1) = BM(2,(p-1)*m+1) - ReductionDOFsArray(j,2*p-1) * StrainBasis(j,2)

            BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) - ReductionDOFsArray(j,2*p) * StrainBasis(j,1)
            BM(2,(p-1)*m+2) = BM(2,(p-1)*m+2) - ReductionDOFsArray(j,2*p) * StrainBasis(j,2)

            IF (.NOT. NoNormalStress) THEN
              BM(4,(p-1)*m+1) = BM(4,(p-1)*m+1) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*p-1) * StrainBasis(j,2)/A2**2)
              BM(4,(p-1)*m+2) = BM(4,(p-1)*m+2) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*p) * StrainBasis(j,2)/A2**2)
            END IF
          END DO
        END DO

        ! Apply strain reduction to vector (C212 u2, C221 u1). Notice the order of components.
        IF (ReductionMethod == DoubleReduction) THEN
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              3, nd, MITC, StretchParMat4)
          ReductionDOFsArray(1:2,1:2*nd) = MATMUL(DOFsTransform(1:2,1:3), ReductionDOFsArray(1:3,1:2*nd))
        ELSE
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              ReducedStrainDim, nd, ReductionMethod, StretchParMat4)
        END IF
        DO p=1,nd
          DO j=1,ReducedStrainDim
            BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) - ReductionDOFsArray(j,2*p-1) * StrainBasis(j,1)
            BM(2,(p-1)*m+2) = BM(2,(p-1)*m+2) - ReductionDOFsArray(j,2*p-1) * StrainBasis(j,2)

            BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) - ReductionDOFsArray(j,2*p) * StrainBasis(j,1)
            BM(2,(p-1)*m+1) = BM(2,(p-1)*m+1) - ReductionDOFsArray(j,2*p) * StrainBasis(j,2)

            IF (.NOT. NoNormalStress) THEN
              BM(4,(p-1)*m+1) = BM(4,(p-1)*m+1) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*p) * StrainBasis(j,2)/A2**2)
              BM(4,(p-1)*m+2) = BM(4,(p-1)*m+2) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*p-1) * StrainBasis(j,2)/A2**2)
            END IF
          END DO
        END DO

        ! Apply strain reduction to vector (0, B22 u3):
        IF (ReductionMethod == DoubleReduction) THEN
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              3, nd, MITC, StretchParMat6)
          ReductionDOFsArray(1:2,1:2*nd) = MATMUL(DOFsTransform(1:2,1:3), ReductionDOFsArray(1:3,1:2*nd))
        ELSE
          CALL ReductionOperatorDofs(Element, Nodes, ReductionDOFsArray, &
              ReducedStrainDim, nd, ReductionMethod, StretchParMat6)
        END IF
        DO p=1,nd
          DO j=1,ReducedStrainDim
            BM(3,(p-1)*m+3) = BM(3,(p-1)*m+3) - ReductionDOFsArray(j,2*p) * StrainBasis(j,1)
            BM(2,(p-1)*m+3) = BM(2,(p-1)*m+3) - ReductionDOFsArray(j,2*p) * StrainBasis(j,2)

            IF (.NOT. NoNormalStress) THEN
              BM(4,(p-1)*m+3) = BM(4,(p-1)*m+3) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*p) * StrainBasis(j,2)/A2**2)
            END IF
          END DO
        END DO

        ! If normal displacement is augmented with a bubble, we need to compute an augmented DOFs array
        ! This is commented out as this gives poor results.
        !IF (Bubbles .AND. ReductionMethod == DoubleReduction) THEN
        !  CALL ReductionOperatorBubbleDofs(Element, Nodes, ReductionDOFsArray(1:2,2*nd+1:2*nd+2), &
        !      ReducedStrainDim, 1, Family, CurlKernel, StretchParMat6)
        !
        !  DO j=1,ReducedStrainDim
        !    BM(3,nd*m+3) = BM(3,nd*m+3) - ReductionDOFsArray(j,2*nd+2) * StrainBasis(j,1)
        !    BM(2,nd*m+3) = BM(2,nd*m+3) - ReductionDOFsArray(j,2*nd+2) * StrainBasis(j,2)          
        !    IF (.NOT. NoNormalStress) THEN
        !      BM(4,nd*m+3) = BM(4,nd*m+3) - nu/(1.0d0-nu) * (-ReductionDOFsArray(j,2*nd+2) * StrainBasis(j,2)/A2**2)
        !    END IF
        !  END DO
        !END IF

        ! The derivative terms that are invariant:
        DO p=1,nd
          BM(1,(p-1)*m+1) = BM(1,(p-1)*m+1) + dBasis(p,1)
          BM(2,(p-1)*m+2) = BM(2,(p-1)*m+2) + dBasis(p,2)
          BM(3,(p-1)*m+1) = BM(3,(p-1)*m+1) + dBasis(p,2)
          BM(3,(p-1)*m+2) = BM(3,(p-1)*m+2) + dBasis(p,1)

          IF (.NOT. NoNormalStress) THEN
            BM(4,(p-1)*m+6) = BM(4,(p-1)*m+6) + Basis(p)
            BM(4,(p-1)*m+1) = BM(4,(p-1)*m+1) - nu/(1.0d0-nu) * (dBasis(p,1)/A1**2)
            BM(4,(p-1)*m+2) = BM(4,(p-1)*m+2) - nu/(1.0d0-nu) * (dBasis(p,2)/A2**2)
          END IF
        END DO

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

          !----------------------------------------------------------------------
          ! Normal stress T^{33} via energy principle: We add a term of the type
          ! e * T^{33}(e)
          !----------------------------------------------------------------------
          IF (.NOT. NoNormalStress) THEN
            BM(4,(p-1)*m+6) = Basis(p)
            BM(4,(p-1)*m+1) = -nu/(1.0d0-nu) * (dBasis(p,1)/A1**2 - C111 * Basis(p)/A1**2 - C221 * Basis(p)/A2**2)
            BM(4,(p-1)*m+2) = -nu/(1.0d0-nu) * (-C112 * Basis(p)/A1**2 + dBasis(p,2)/A2**2 - C222 * Basis(p)/A2**2)
            BM(4,(p-1)*m+3) = -nu/(1.0d0-nu) * (-B11 * Basis(p)/A1**2 - B22 * Basis(p)/A2**2)
          END IF

        END DO
      END IF

      CALL StrainEnergyDensity(Stiff, CMat, BM, &
          csize, DOFs, Weight)

      !---------------------------------------------------------------
      ! The part corresponding to the shear strains:
      !--------------------------------------------------------------- 
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
            BS(1:2,nd*m+1) = BS(1:2,nd*m+1) - ReductionDOFsArray(j,2*nd+1) * StrainBasis(j,1:2)
            BS(1:2,nd*m+2) = BS(1:2,nd*m+2) - ReductionDOFsArray(j,2*nd+2) * StrainBasis(j,1:2)          
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
          BS(1,nd*m+1) = -Basis(nd+1)
          BS(2,nd*m+2) = -Basis(nd+1)
        END IF
      END IF

      CALL StrainEnergyDensity(Stiff, GMat, BS, &
          2, DOFs, Kappa*Weight)


      !---------------------------------------------------------------
      ! The part corresponding to the bending strains (terms depending
      ! on the membrane strains dropped at the moment):
      !---------------------------------------------------------------      
      Weight = h**3/12.0d0 * SqrtDetA * detJ * sq
      DO p=1,nd
        BB(1,(p-1)*m+4) = dBasis(p,1) - C111 * Basis(p)
        BB(1,(p-1)*m+5) = -C112 * Basis(p)

        BB(2,(p-1)*m+4) = -C221 * Basis(p)
        BB(2,(p-1)*m+5) = dBasis(p,2) - C222 * Basis(p)

        BB(3,(p-1)*m+4) = dBasis(p,2) - 2.0d0 * C211 * Basis(p)
        BB(3,(p-1)*m+5) = dBasis(p,1) - 2.0d0 * C212 * Basis(p)
        BB(3,(p-1)*m+1) = B11/a11 * (dBasis(p,2) - C211 * Basis(p)) - B22/a22 * C211 * Basis(p)
        BB(3,(p-1)*m+2) = B11/a11 * (- C212) * Basis(p) + B22/a22 * (dBasis(p,1) - C212 * Basis(p))
      END DO

      IF (UseBubbles) THEN
        ! With rotation bubbles:
        BB(1,nd*m+1) = dBasis(nd+1,1) - C111 * Basis(nd+1)
        BB(1,nd*m+2) = -C112 * Basis(nd+1)
        BB(2,nd*m+1) = -C221 * Basis(nd+1)
        BB(2,nd*m+2) = dBasis(nd+1,2) - C222 * Basis(nd+1)
        BB(3,nd*m+1) = dBasis(nd+1,2) - 2.0d0 * C211 * Basis(nd+1)
        BB(3,nd*m+2) = dBasis(nd+1,1) - 2.0d0 * C212 * Basis(nd+1)
      END IF

      CALL StrainEnergyDensity(Stiff, CMat, BB, &
          3, DOFs, Weight)

      IF (NoNormalStress) THEN
        !----------------------------------------------------------------
        ! This does not include the action of strain reduction operators
        ! and is becoming obsolite.
        ! Determine the thickness-stretch parameter corresponding to the
        ! state of vanishing normal stress:
        !----------------------------------------------------------------
        Weight = SqrtDetA * detJ * sq
        DO p=1,nd
          DO k=1,nd
            Stiff(p*6,k*6) = Stiff(p*6,k*6) + Basis(p) * Basis(k) * Weight
            Stiff(p*6,(k-1)*m+1) = Stiff(p*6,(k-1)*m+1) - nu/(1.0d0-nu) * Basis(p) * &
                (  dBasis(k,1)/A1**2 - C111 * Basis(k)/A1**2 - C221 * Basis(k)/A2**2 ) * Weight
            Stiff(p*6,(k-1)*m+2) = Stiff(p*6,(k-1)*m+2) - nu/(1.0d0-nu) * Basis(p) * &
                ( -C112 * Basis(k)/A1**2 + dBasis(k,2)/A2**2 - C222 * Basis(k)/A2**2 ) * Weight
            Stiff(p*6,(k-1)*m+3) = Stiff(p*6,(k-1)*m+3) - nu/(1.0d0-nu) * Basis(p) * &
                ( -B11 * Basis(k)/A1**2 - B22 * Basis(k)/A2**2 ) * Weight
          END DO
        END DO
      END IF

      !----------------------------------------------------------------
      ! Mass matrix without bubbles taken into account:
      !----------------------------------------------------------------     
      IF (TransientSimulation) THEN
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
            END DO
          END DO
        END DO
      END IF

      !----------------------------------------------------------------
      ! RHS vector:
      !----------------------------------------------------------------
      Weight = SqrtDetA * detJ * sq
      DO p=1,nd
        i = m*(p-1)+3
        Force(i) = Force(i) + NormalTraction * Basis(p) * Weight
      END DO

      IF (UseBubbles .AND. nb > 2) THEN
        i = m*nd+3
        Force(i) = Force(i) + NormalTraction * Basis(nd+1) * Weight
      END IF

      Area = Area + sq * detJ * SqrtDetA

      ! For computing the mean curvature error:
      !-----------------------------------------
      !Error = Error + sq * (1.0d0 - abs(B11/a11) - abs(B22/a22))**2 * detJ * SqrtDetA

    END DO QUADRATURELOOP

    ! ------------------------------------------------------------------------------
    ! Static condensation is performed before transforming to the global DOFs:
    ! ------------------------------------------------------------------------------
    IF (UseBubbles .AND. nb>0) THEN
      CALL CondensateP( m*nd, nb, Stiff, Force )
    END IF


    Q = 0.0d0
    DO j=1,nd
      !IF ( PVersion ) THEN
      IF ( .FALSE. ) THEN
        y1 = PNodes % x(j)
        y2 = PNodes % y(j)       
      ELSE
        y1 = Nodes % x(j)
        y2 = Nodes % y(j)
      END IF

      CALL SurfaceBasisVectors(y1, y2, TaylorParams, e1, e2, e3, o, abasis1, &
          abasis2, abasis3, A11, A22, SqrtDetA, B11, B22, C111, C112, C221, C222, &
          C211, C212, PlanarPoint=PlateBody, UmbilicalPoint=SphericalSurface)

      QBlock(1,1:3) = abasis1(1:3)
      QBlock(2,1:3) = abasis2(1:3)
      QBlock(3,1:3) = abasis3(1:3)         

      i0 = (j-1)*m

      Q(i0+1:i0+3,i0+1:i0+3) =  QBlock(1:3,1:3)
      Q(i0+4:i0+6,i0+4:i0+6) =  QBlock(1:3,1:3)

    END DO

    !-------------------------------------------------------
    ! Transform to the global DOFs:
    !-------------------------------------------------------
    k = m*nd
    Stiff(1:k,1:k) = MATMUL(TRANSPOSE(Q(1:k,1:k)),MATMUL(Stiff(1:k,1:k),Q(1:k,1:k)))
    Force(1:k) = MATMUL(TRANSPOSE(Q(1:k,1:k)),Force(1:k))

    IF(TransientSimulation) THEN
      Mass(1:k,1:k) = MATMUL(TRANSPOSE(Q(1:k,1:k)),MATMUL(Mass(1:k,1:k),Q(1:k,1:k)))
      CALL Default2ndOrderTime(MASS,DAMP,STIFF,FORCE)
    END IF

    CALL DefaultUpdateEquations(STIFF,FORCE)

!------------------------------------------------------------------------------
  END SUBROUTINE ShellLocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ShearCorrectionFactor(Kappa,Thickness,x,y,n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: Kappa,Thickness,x(:),y(:)
    INTEGER :: n
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: x21,x32,x43,x13,x14,y21,y32,y43,y13,y14, &
        l21,l32,l43,l13,l14,alpha,h
!------------------------------------------------------------------------------
    Kappa = 1.0d0
    SELECT CASE(n)
    CASE(3)
      alpha = 0.20d0
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
      alpha = 0.10d0
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
      CALL WARN('ShearCorrectionFactor','Illegal number of nodes for Smitc elements')
    END SELECT
!------------------------------------------------------------------------------
  END SUBROUTINE ShearCorrectionFactor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! The matrix representation of the elasticity tensor with respect an orthogonal
! basis. The case A1 = A2 = 1 corresponds to an orthonormal basis.     
!------------------------------------------------------------------------------
  SUBROUTINE ElasticityMatrix(CMat, GMat, A1, A2, E, nu)
!------------------------------------------------------------------------------    
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(OUT) :: CMat(4,4), GMat(2,2)
    REAL(KIND=dp), INTENT(IN) :: A1, A2, E, nu  
!------------------------------------------------------------------------------
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

    ! The row corresponding to the normal stress: A deviation from the state of
    ! vanishing normal stress produces deformation energy as described by
    ! the 3-D Hooke's law.
    CMat(4,4) = (1.0d0-nu) * E /( (1.0d0+nu) * (1.0d0-2.0d0*nu) )

    GMat(1,1) = E/(2.0d0*(1.0d0 + nu)*A1**2)
    GMat(2,2) = E/(2.0d0*(1.0d0 + nu)*A2**2)
!------------------------------------------------------------------------------
  END SUBROUTINE ElasticityMatrix
!------------------------------------------------------------------------------    

!------------------------------------------------------------------------------
! Perform the operation
!
!    A = A + C' * B * C * s
!
! with
!
!    Size( A ) = n x n
!    Size( B ) = m x m
!    Size( C ) = m x n
!------------------------------------------------------------------------------
  SUBROUTINE StrainEnergyDensity(A, B, C, m, n, s)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(INOUT) :: A(:,:)
    REAL(KIND=dp), INTENT(IN) :: B(:,:), C(:,:)
    INTEGER, INTENT(IN) :: m, n
    REAL(KIND=dp), INTENT(IN) :: s
!------------------------------------------------------------------------------
    A(1:n,1:n) = A(1:n,1:n) + s * MATMUL(TRANSPOSE(C(1:m,1:n)),MATMUL(B(1:m,1:m),C(1:m,1:n))) 
!------------------------------------------------------------------------------
  END SUBROUTINE StrainEnergyDensity
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
!------------------------------------------------------------------------------
  FUNCTION ReductionOperatorInfo( Element, Nodes, u, v, StrainBasis, ReductionMethod, &
      ApplyPiolaTransform, F, G, detF, Basis, dBasis, DOFWeigths, Bubbles) RESULT(stat)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), INTENT(IN), TARGET :: Element         !< Element structure
    TYPE(Nodes_t), INTENT(IN) :: Nodes                     !< Data corresponding to the classic element nodes
    REAL(KIND=dp), INTENT(IN) :: u                         !< 1st reference element coordinate
    REAL(KIND=dp), INTENT(IN) :: v                         !< 2nd reference element coordinate
    REAL(KIND=dp), INTENT(OUT) :: StrainBasis(:,:)         !< The basis functions b spanning the reference element space
    INTEGER, INTENT(IN) :: ReductionMethod                 !< The method chosen: (2=MITC, 1=kernel version)
    LOGICAL, INTENT(IN), OPTIONAL :: ApplyPiolaTransform   !< If  .TRUE., perform the Piola transform so that, instead of b(p),
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
    LOGICAL :: Stat                                        !< Currently a dummy return value
!---------------------------------------------------------------------------------
    LOGICAL :: PerformPiolaTransform, CreateBubbles

    INTEGER :: dim, i, j, k, n, ntot, q, DOFs, Family

    REAL(KIND=dp) :: LF(3,3), LG(3,3), detLF, B(3)
    REAL(KIND=dp) :: LBasis(Element % TYPE % NumberOfNodes+1)
    REAL(KIND=dp) :: dLBasis(Element % TYPE % NumberOfNodes+1,3)
!------------------------------------------------------------------------------       
    StrainBasis = 0.0d0
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
    ! The standard nodal basis functions on the reference element and
    ! their derivatives with respect to the local coordinates. These define 
    ! the mapping of the reference element to a physical element.
    !-----------------------------------------------------------------------
    LBasis = 0.0d0
    dLBasis = 0.0d0      
    SELECT CASE(Family)
    CASE(3)
      DO q=1,n
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

    SELECT CASE(Family)
    CASE(3)
      IF (ReductionMethod == 2) THEN
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
      ELSE
        !---------------------------------------------------------------------
        ! The basis functions for RT_0(k,0), with DOFs defined as integrals
        ! of the type d_i = (u,v_i)_k. Here the given function v_i tranforms
        ! according to the standard Piola transformation (the div-conforming 
        ! version).
        !---------------------------------------------------------------------
        DOFs = 2
        StrainBasis(1,1) = 1.0d0
        StrainBasis(1,2) = 0.0d0
        StrainBasis(2,1) = 0.0d0
        StrainBasis(2,2) = 1.0d0
        
        IF ( PRESENT(DOFWeigths) ) THEN
          DOFWeigths = 0.0d0
          DOFWeigths(1,1:2) = StrainBasis(1,1:2)/sqrt(3.0d0)
          DOFWeigths(2,1:2) = StrainBasis(2,1:2)/sqrt(3.0d0)
        END IF
      END IF
    CASE(4)
      IF (ReductionMethod == 2) THEN
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
      ELSE
        !---------------------------------------------------------------------
        ! The basis functions for ABF_0(k,0), with DOFs defined as integrals
        ! of the type d_i = (u,v_i)_k. Here the given function v_i tranforms
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
          DOFWeigths = 0.0d0
          DOFWeigths(1,1:2) = StrainBasis(1,1:2)/4.0d0
          DOFWeigths(2,1:2) = StrainBasis(2,1:2)/4.0d0
          DOFWeigths(3,1:2) = StrainBasis(3,1:2)*3.0d0/8.0d0
        END IF
      END IF

    CASE DEFAULT
      CALL Warn('ReductionOperatorInfo','Unsupported element type')
      RETURN
    END SELECT

    IF (PerformPiolaTransform) THEN
      DO j=1,DOFs
        DO k=1,dim
          B(k) = SUM( LG(k,1:dim) * StrainBasis(j,1:dim) )
        END DO
        StrainBasis(j,1:dim) = B(1:dim)
      END DO
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
!  functional defining the kth DOF, the kth row of the returned matrix has 
!  entries 
!
!     [d_k(N1*e1) d_k(N1*e2) ... d_k(Nn*e1) d_k(Nn*e2)]
!
!  where N1,...,Nn are the Lagrange basis functions, e1=(1,0) and e2=(0,1). 
!  Thus, if the DOFs of u are contained accordingly in a vector U, the DOFs of
!  the interpolating function can be evaluated by computing the product AU.
!  Optionally the interpolating function can be computed for a field Cu
!  where C is a 2X2 matrix field. The model parameters at the element nodes
!  are then given in the ModelPar array. This subroutine depends intimately
!  on the function ReductionOperatorInfo which gives the definition of DOFs
!  integrated here. 
!------------------------------------------------------------------------------
  SUBROUTINE ReductionOperatorDofs(Element, Nodes, A, nd, n, ReductionMethod, &
      ModelPars)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), INTENT(IN), TARGET :: Element          !< Element structure
    TYPE(Nodes_t), INTENT(IN) :: Nodes                      !< Nodes structure
    REAL(KIND=dp), INTENT(INOUT) :: A(:,:)                 !< Coefficients for expressing the DOFs 
    INTEGER, INTENT(IN) :: nd                               !< The dimension of the strain reduction space X(K)
    INTEGER, INTENT(IN) :: n                                !< The number of the H1-conforming basis functions
    INTEGER, INTENT(IN) :: ReductionMethod                  !< The method chosen
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: ModelPars(2,2,n) !< To include the effect of additional model parameters
!---------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: stat, UseParameters, PRefElement

    INTEGER :: Family, i, j, k, t

    REAL(KIND=dp) :: StrainBasis(nd,3)        ! The basis functions for the strain reduction space X(K)
    REAL(KIND=dp) :: DOFWeigths(nd,3)         ! The auxiliary functions to evaluate the interpolant in X(K)
    REAL(KIND=dp) :: Basis(n)                 ! H1-conforming basis functions
    REAL(KIND=dp) :: F(3,3)                   ! The gradient F=Grad f of the element mapping
    REAL(KIND=dp) :: detJ
    REAL(KIND=dp) :: u(2), ParMat(2,2), uk, vk, sk
    REAL(KIND=dp) :: Tau(4,3)
!---------------------------------------------------------------------------------
    Family = GetElementFamily(Element)

    ! Clear upper left corner of A
    A(1:nd,1:2*n) = 0.0d0

    UseParameters = PRESENT(ModelPars) 
    IF (.NOT. UseParameters) THEN
      ParMat(1,1) = 1.0d0
      ParMat(1,2) = 0.0d0
      ParMat(2,1) = 0.0d0
      ParMat(2,2) = 1.0d0
    END IF

    SELECT CASE(Family)
    CASE(3)
      PRefElement = IsPElement(Element)
      IF (ReductionMethod == 2) THEN
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

      ELSE
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
      END IF

    CASE(4)
      IF (ReductionMethod == 2) THEN
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

      ELSE
        ! Method: the kernel related to ABF
        IP = GaussPoints(Element)
        DO t=1,IP % n

          !stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          !    IP % W(t), detJ, Basis )

          stat = ReductionOperatorInfo( Element, Nodes, IP % U(t), IP % V(t), StrainBasis, &
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
                A(i,2*(j-1)+k) = A(i,2*(j-1)+k) + SUM(u(1:2) * DOFWeigths(i,1:2)) * IP % s(t)
              END DO
            END DO
          END DO
        END DO
      END IF
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
!  where Nb is the bubble basis function, e1=(1,0) and e2=(0,1). 
!  Currently, just one bubble function (cf. the size of Basis array) and the kernel 
!  version of strain reduction are supported currently. 
!------------------------------------------------------------------------------
  SUBROUTINE ReductionOperatorBubbleDofs(Element, Nodes, A, nd, nb, n, ReductionMethod, &
      ModelPars)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), INTENT(IN), TARGET :: Element          !< Element structure
    TYPE(Nodes_t), INTENT(IN) :: Nodes                      !< Nodes structure
    REAL(KIND=dp), INTENT(INOUT) :: A(nd,2*nb)                !< Coefficients for expressing the DOFs 
    INTEGER, INTENT(IN) :: nd                               !< The dimension of the strain reduction space X(K)
    INTEGER, INTENT(IN) :: nb                               !< The number of the H1-conforming bubble functions
    INTEGER, INTENT(IN) :: n                                !< The number of the BG element nodes
    INTEGER, INTENT(IN) :: ReductionMethod                  !< The method chosen
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: ModelPars(2,2,n) !< To include the effect of additional model parameters
!---------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: stat, UseParameters, PRefElement

    INTEGER :: Family, i, j, k, t

    REAL(KIND=dp) :: StrainBasis(4,3)         ! The basis functions for the strain reduction space X(K)
    REAL(KIND=dp) :: DOFWeigths(3,2)          ! The auxiliary functions to evaluate the interpolant in X(K)
    REAL(KIND=dp) :: Basis(5)                 ! H1-conforming basis functions (p=1, with one bubble)
    REAL(KIND=dp) :: u(2), ParMat(2,2), uk, vk, sk
!---------------------------------------------------------------------------------
    IF (ReductionMethod /= CurlKernel) CALL Fatal('ReductionOperatorBubbleDofs', &
        'An unsupported strain reduction technique')

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

      PRefElement = IsPElement(Element)
      
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

      stat = ReductionOperatorInfo( Element, Nodes, uk, vk, StrainBasis, &
          ReductionMethod, ApplyPiolaTransform = .TRUE., Basis=Basis, DOFWeigths=DOFWeigths, &
          Bubbles=.TRUE.)  

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
              u(1) = ParMat(1,1)*Basis(n+j)
              u(2) = ParMat(2,1)*Basis(n+j)
            CASE(2)
              u(1) = ParMat(1,2)*Basis(n+j)
              u(2) = ParMat(2,2)*Basis(n+j)
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
! This function returns a description of superparametric element 
! (an experimental feature) 
!------------------------------------------------------------------------------
  FUNCTION SuperParametricElementInfo( Element, GElement, GBasis, GNodesX, &
      GNodesY, u, v, detF, Basis, dBasis, dGBasis, ReadyBasis ) RESULT(stat)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Element_t), TARGET :: Element             !< Element structure
    TYPE(Element_t), TARGET :: GElement            !< Geometry element structure
    REAL(KIND=dp) :: GBasis(:)                     !< Geometry element basis function values at p=(u,v)
    REAL(KIND=dp) :: GNodesX(:)                    !< the X-coordinates of GElement nodes
    REAL(KIND=dp) :: GNodesY(:)                    !< the Y-coordinates of GElement nodes    
    REAL(KIND=dp) :: u                             !< 1st local coordinate at which to calculate the basis functions
    REAL(KIND=dp) :: v                             !< 2nd local coordinate
    REAL(KIND=dp) :: detF                          !< Metric term for integration
    REAL(KIND=dp) :: Basis(:)                      !< Basis function values at p=(u,v)
    REAL(KIND=dp) :: dBasis(:,:)                   !< Global first derivatives of basis functions at p
    REAL(KIND=dp), OPTIONAL :: dGBasis(:,:)        !< Global first derivatives of geometry basis functions at p
    LOGICAL, OPTIONAL :: ReadyBasis                !< Transform the given derivatives of the reference element basis
    LOGICAL :: Stat                                !< A dummy variable at the moment
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
    LOGICAL :: TransformDerivatives
    INTEGER :: n, dim, cdim, i, j, k
    
    REAL(KIND=dp) :: dRefBasis(SIZE(GNodesX),3)
    REAL(KIND=dp) :: F(2,2), G(2,2)
    REAL(KIND=dp) :: w 
!------------------------------------------------------------------------------
    TransformDerivatives = .FALSE.
    IF ( PRESENT(ReadyBasis) ) TransformDerivatives = ReadyBasis

    w = 0.0d0

    n  = SIZE(GNodesX)
    dim  = GElement % TYPE % DIMENSION
    cdim = dim
    ! cdim = CoordinateSystemDimension()

    GBasis = 0.0d0
    CALL NodalBasisFunctions(n, GBasis, GElement, u, v, w)

    dRefBasis = 0.0d0
    CALL NodalFirstDerivatives(n, dRefBasis, GElement, u, v, w)

    !------------------------------------------------------------------------------
    ! The gradient of the element mapping K = f(k), with k the reference element
    !------------------------------------------------------------------------------
    F = 0.0d0
    DO i=1,dim
      F(1,i) = SUM( GNodesX(1:n) * dRefBasis(1:n,i) )
      F(2,i) = SUM( GNodesY(1:n) * dRefBasis(1:n,i) )
    END DO
    DetF = F(1,1)*F(2,2) - F(1,2)*F(2,1)

    ! --------------------------------
    ! The global first derivatives:
    ! --------------------------------
    G(1,1) = 1.0d0/detF * F(2,2)
    G(1,2) = -1.0d0/detF * F(1,2)
    G(2,1) = -1.0d0/detF * F(2,1)
    G(2,2) = 1.0d0/detF * F(1,1)
    G(1:dim,1:dim) = TRANSPOSE( G(1:dim,1:dim) )

    IF (PRESENT( dGBasis )) THEN
      dGBasis = 0.0d0
      DO i=1,n
        DO j=1,cdim
          DO k=1,dim
            dGBasis(i,j) = dGBasis(i,j) + dRefBasis(i,k)*G(j,k)
          END DO
        END DO
      END DO
    END IF

    IF (TransformDerivatives) THEN
      n = SIZE(Basis)
      dRefBasis(1:n,1:dim) = dBasis(1:n,1:dim)    
    ELSE
      n = Element % Type % NumberOfNodes
      Basis = 0.0d0
      CALL NodalBasisFunctions(n, Basis, Element, u, v, w) 

      dRefBasis = 0.0d0
      CALL NodalFirstDerivatives(n, dRefBasis, Element, u, v, w)    
    END IF

    dBasis = 0.0d0
    DO i=1,n
      DO j=1,cdim
        DO k=1,dim
          dBasis(i,j) = dBasis(i,j) + dRefBasis(i,k)*G(j,k)
        END DO
      END DO
    END DO

    DetF = ABS(DetF)
    stat = .TRUE.
!------------------------------------------------------------------------------
  END FUNCTION SuperParametricElementInfo
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

    DirectorValues => GetElementProperty('director', Element)
    IF (ASSOCIATED(DirectorValues)) THEN
      IF (SIZE(DirectorValues) < 3*n) CALL Fatal('AverageDirector', &
          'Elemental director data is not associated with all nodes')
    ELSE
      CALL Fatal( 'AverageDirector', 'Elemental director data is not associated')
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
      stat = ElementInfo( PlaneElement, NodesVar, IP % u(j), IP % v(j), &
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
  SUBROUTINE ComputeSurfaceArea(Element, SurfaceArea, MacroElement)
!--------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    REAL(KIND=dp), INTENT(INOUT) :: SurfaceArea
    LOGICAL, INTENT(IN) :: MacroElement
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
      IP = GaussPointsTriangle(11, PReferenceElement=.TRUE. )
    CASE(4)
      IP = GaussPoints(Element, 25)
    CASE DEFAULT
      RETURN
    END SELECT

    DO j=1,IP % n
      stat = BlendingSurfaceInfo( Element, Nodes, IP % U(j), IP % V(j), &
          DetA, a1, a2, a3, A, B, x, MacroElement)      
      SurfaceArea = SurfaceArea + IP % s(j) * SQRT(Deta)
    END DO
!-------------------------------------------------------------------------------------
  END SUBROUTINE ComputeSurfaceArea
!-------------------------------------------------------------------------------------


!---------------------------------------------------------------------------
  FUNCTION LDot(a,b) RESULT(res)
!---------------------------------------------------------------------------
    REAL(KIND=dp) :: a(3), b(3), res   
!---------------------------------------------------------------------------
    res = SUM( a(:) * b(:) )
!---------------------------------------------------------------------------
  END FUNCTION LDot
!---------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE ShellSolver
!------------------------------------------------------------------------------
