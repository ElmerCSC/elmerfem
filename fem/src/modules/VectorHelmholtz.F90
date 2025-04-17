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
! *  Authors: Juhani Kataja, Juha Ruokolainen, Mika Malinen and Peter Råback
! *  Email:   juhani.kataja@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 26 Sep 2014
! *  Heavily inspired from the MagnetoDynamics module.
! *****************************************************************************/
    
!------------------------------------------------------------------------------
!> Solve the time-harmonic Maxwell equations using the curl-curl equation at a 
!> relatively high frequency using curl-conforming edge elements. Also a low-
!> frequency model is available. With "Use Gauss Law = True" the A-V representation
!> is used and to solve the additional scalar potential V the Gauss law is used.
!> \ingroup Solvers
!-------------------------------------------------------------------------------
MODULE VectorHelmholtzUtils

   USE DefUtils
   IMPLICIT NONE

   COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)   
  
CONTAINS

!------------------------------------------------------------------------------
  FUNCTION ComplexCrossProduct(v1,v2) RESULT(v3)
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: v1(3), v2(3), v3(3)
    v3(1) =  v1(2)*v2(3) - v1(3)*v2(2)
    v3(2) = -v1(1)*v2(3) + v1(3)*v2(1)
    v3(3) =  v1(1)*v2(2) - v1(2)*v2(1)
!------------------------------------------------------------------------------
  END FUNCTION ComplexCrossProduct

END MODULE VectorHelmholtzUtils


!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE VectorHelmholtzSolver_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE VectorHelmholtzUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: Found, SecondOrder, PiolaVersion, SecondFamily, WithNDOFs

  SolverParams => GetSolverParams()  

  WithNDOFs = GetLogical(SolverParams, 'Use Gauss Law', Found)
  IF (.NOT. WithNDOFs) THEN
    WithNDOFs = GetLogical(SolverParams, 'Use Lagrange Gauge', Found)
    WithNDOFs = WithNDOFs .OR. GetLogical(SolverParams, 'Lorenz Condition', Found)
  END IF
  
  IF ( .NOT.ListCheckPresent(SolverParams, "Element") ) THEN
    ! We use one place where all the edge element keywords are defined and checked.
    CALL EdgeElementStyle(SolverParams, PiolaVersion, SecondFamily, SecondOrder, Check = .TRUE. )
    
    IF (WithNDOFs) THEN
      IF ( SecondOrder ) THEN
        CALL ListAddString( SolverParams, "Element", &
            "n:1 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2" )
      ELSE IF( SecondFamily ) THEN
        CALL ListAddString( SolverParams, "Element", "n:1 e:2" )
      ELSE IF( PiolaVersion ) THEN
        CALL ListAddString( SolverParams, "Element", "n:1 e:1 -quad b:2 -brick b:3 -quad_face b:2" )
      ELSE
        CALL ListAddString( SolverParams, "Element", "n:1 e:1" )
      END IF      
    ELSE
      IF( SecondOrder ) THEN
        CALL ListAddString( SolverParams, "Element", &
            "n:0 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2" )
      ELSE IF (SecondFamily) THEN
        CALL ListAddString( SolverParams, "Element", "n:0 e:2" )
      ELSE IF( PiolaVersion ) THEN
        CALL ListAddString( SolverParams, "Element", "n:0 e:1 -quad b:2 -brick b:3 -quad_face b:2" )
      ELSE
        CALL ListAddString( SolverParams, "Element", "n:0 e:1" )
      END IF
    END IF
  END IF

  CALL ListAddNewLogical(SolverParams, 'Bubbles in Global System', .TRUE.)

  !CALL ListAddNewLogical( SolverParams,'Hcurl Basis',.TRUE.)
  IF (WithNDOFs) THEN
    CALL ListAddNewLogical(SolverParams,'Variable Output',.TRUE.)
    CALL ListAddNewString(SolverParams,'Variable','AV[AV re:1 AV im:1]')
  ELSE
    CALL ListAddNewLogical( SolverParams,'Variable Output',.FALSE.)
    CALL ListAddNewString( SolverParams,'Variable','E[E re:1 E im:1]')
  END IF
  CALL ListAddNewLogical( SolverParams, "Linear System Complex", .TRUE.)


  ! These use one flag to call library features to compute automatically
  ! a capacitance matrix.
  IF( ListGetLogical( SolverParams,'Calculate S Matrix',Found ) ) THEN
    CALL Info('VectorHelmholtz_init','Using Constraint Modes functionality for S Matrix')
    CALL ListAddNewLogical( SolverParams,'Constraint Modes Analysis',.TRUE.)
  END IF

  IF( ListGetLogical( SolverParams,'Constraint Modes Analysis', Found ) ) THEN
    CALL ListAddNewLogical( SolverParams,'Constraint Modes Lumped',.TRUE.)
    CALL ListAddNewLogical( SolverParams,'Constraint Modes Fluxes',.TRUE.)
    !CALL ListAddNewLogical( SolverParams,'Constraint Modes Matrix Results',.TRUE.)
    CALL ListAddNewLogical( SolverParams,'Constraint Modes EM Wave',.TRUE.)        
    IF( ListCheckPresent( SolverParams,'S Matrix Filename') ) THEN
      CALL ListRename( SolverParams,'S Matrix Filename',&
          'Constraint Modes Matrix Filename', Found ) 
    ELSE     
      CALL ListAddNewString( SolverParams,'Constraint Modes Matrix Filename',&
          'SMatrix.dat',.FALSE.)
    END IF
    IF( ListCheckPresent( SolverParams,'S Matrix Symmetric' ) ) THEN
      CALL ListRename( SolverParams,'S Matrix Symmetric','Constraint Modes Fluxes Symmetric')
    END IF

    CALL ListRenameAllBC( Model,'S Matrix Port','Constraint Mode')
    !CALL ListAddLogical( Params,'Optimize Bandwidth',.FALSE.)
    !CALL Info('VectorHelmoltz_init','Suppressing bandwidth optimization in S-Matrix computation!')
  END IF
  
!------------------------------------------------------------------------------
END SUBROUTINE VectorHelmholtzSolver_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE VectorHelmholtzSolver_Init(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE VectorHelmholtzUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: Found
  INTEGER :: i, j, soln
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  
  ! A parallel run needs an early allocation of precvalues. Create a flag to
  ! inform the function CreateMatrix.
  IF (ListCheckPrefix(SolverParams, 'Linear System Preconditioning Damp Coefficient')) THEN
    CALL ListAddNewLogical(SolverParams, 'Allocate Preconditioning Matrix', .TRUE.)
  END IF

  !
  ! The following is for creating sources from pre-computed eigenfunctions:
  !
  IF (ListGetLogicalAnyBC(Model, 'Eigenfunction BC')) THEN
    soln = 0
    DO i=1,Model % NumberOfSolvers
      sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
      j = INDEX(sname, 'EMPortSolver')
      IF (j > 0) THEN
        soln = i
        EXIT
      END IF
    END DO

    IF( soln == 0 ) THEN
      CALL Fatal('VectorHelmholtzSolver_Init','Eigenfunction BC given without solving a port model')      
    ELSE
      CALL Info('VectorHelmholtzSolver_Init','The eigensolver index is: '//I2S(soln), Level=12)
      CALL ListAddInteger(SolverParams, 'Eigensolver Index', soln)
    END IF
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE VectorHelmholtzSolver_Init
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solve the electric field E from the curl-curl equation 
!> curl (1/mu) curl E - i \omega \sigma E - \omega^2 epsilon E = i omega J
!> using edge elements (vector-valued basis of first or second degree).
!> For the equations in the case of the A-V representation, see the manual.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE VectorHelmholtzSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE VectorHelmholtzUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: Eigensolver => NULL()
  LOGICAL :: Found, PrecMatrix, HasPrecDampCoeff, MassProportional, CurlCurlPrec
  REAL(KIND=dp) :: Omega, mu0inv, eps0, rob0
  INTEGER :: i, soln, NoIterationsMax, EdgeBasisDegree
  TYPE(Mesh_t), POINTER :: Mesh
  COMPLEX(KIND=dp) :: PrecDampCoeff
  LOGICAL :: PiolaVersion, EdgeBasis, LowFrequencyModel, LorenzCondition
  LOGICAL :: UseGaussLaw, ChargeConservation
  LOGICAL :: EigenfunctionSource
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: pSolver
  CHARACTER(*), PARAMETER :: Caller = 'VectorHelmholtzSolver'
!------------------------------------------------------------------------------

  CALL Info(Caller,'',Level=6 )
  CALL Info(Caller,'-------------------------------------------------',Level=6 )
  CALL Info(Caller,'Solving harmonic electromagnetic wave equation!',Level=5 )

  SolverParams => GetSolverParams()  
  CALL EdgeElementStyle(Solver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree ) 

  IF(CoordinateSystemDimension() == 2 .AND. .NOT. PiolaVersion ) THEN
    CALL Fatal(Caller, 'A 2D model needs "Use Piola Transform = True"')
  END IF
    
  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  Mesh => GetMesh()
  pSolver => Solver

  IF( Solver % Variable % dofs /= 2) THEN
    CALL Fatal (Caller, &
        'Variable is not of size two ('//I2S(Solver % Variable % dofs)//'), &
        Use: Variable = E[E re:1 E im:1]')
  ENDIF

  Omega = GetAngularFrequency(Found=Found)
  IF(.NOT. Found ) THEN
    CALL Fatal(Caller,'Harmonic wave equation requires angular frequency!')
  END IF
  
  PrecDampCoeff = GetCReal(SolverParams, 'Linear System Preconditioning Damp Coefficient', HasPrecDampCoeff )
  PrecDampCoeff = CMPLX(REAL(PrecDampCoeff), &
      GetCReal(SolverParams, 'Linear System Preconditioning Damp Coefficient im', Found ) )
  HasPrecDampCoeff = HasPrecDampCoeff .OR. Found 
  IF (HasPrecDampCoeff) THEN
    MassProportional = GetLogical(SolverParams, 'Mass-proportional Damping', Found)
  ELSE
    MassProportional = .FALSE.
  END IF
    
  CurlCurlPrec = GetLogical(SolverParams, 'Curl-Curl Preconditioning', Found)
  PrecMatrix = HasPrecDampCoeff .OR. CurlCurlPrec
  
  IF(PrecMatrix) THEN
    IF(ListGetString(SolverParams,'Linear System Solver',Found ) == 'direct') THEN
      CALL Warn(Caller,'Generating preconditioning matrix does not make sense for direct methods, canceling!')
      PrecMatrix = .FALSE.
    ELSE
      CALL Info(Caller,'Generating special preconditioning matrix',Level=12)
    END IF
  END IF
    
  Found = .FALSE.
  IF( ASSOCIATED( Model % Constants ) ) THEN
    mu0inv = GetConstReal( Model % Constants, 'Permeability of Vacuum', Found )
    IF(mu0inv/=0) mu0inv=1/mu0inv
  END IF
  IF(.NOT. Found ) mu0inv = 1.0_dp / ( PI * 4.0d-7 )
  
  Found = .FALSE.
  IF( ASSOCIATED( Model % Constants ) ) THEN
    IF (ListCheckPresent(Model % Constants, 'Permittivity of Vacuum')) &
        eps0 = GetConstReal ( Model % Constants, 'Permittivity of Vacuum', Found )
  END IF
  IF(.NOT. Found ) eps0 = 8.854187817d-12

  rob0 = Omega * SQRT( eps0 / mu0inv )
  
  LowFrequencyModel = GetLogical(SolverParams, 'Low Frequency Model', Found)
  LorenzCondition = GetLogical(SolverParams, 'Lorenz Condition', Found)
  UseGaussLaw = GetLogical(SolverParams, 'Use Gauss Law', Found)
  ChargeConservation = GetLogical(SolverParams, 'Apply Conservation of Charge', Found)

  EigenfunctionSource = ListGetLogicalAnyBC(Model, 'Eigenfunction BC')
  IF (EigenfunctionSource) THEN
    soln = ListGetInteger(SolverParams, 'Eigensolver Index', Found) 
    IF (soln == 0) THEN
      CALL Fatal(Caller, 'We should know > Eigensolver Index <')
    END IF
    Eigensolver => Model % Solvers(soln)
  END IF
  
  
  ! Resolve internal nonlinearities, if requested:
  ! ----------------------------------------------
  NoIterationsMax = GetInteger( SolverParams, &
      'Nonlinear System Max Iterations',Found)
  IF(.NOT. Found) NoIterationsMax = 1
  
  EdgeBasis = .NOT. ListCheckPresent( SolverParams,'Linear System Refactorize' ) .AND. &
      GetLogical( SolverParams, 'Edge Basis', Found )

  CALL DefaultStart()

  DO i=1,NoIterationsMax
    IF( DoSolve() ) EXIT
    IF( EdgeBasis ) CALL ListAddLogical( SolverParams,'Linear System Refactorize',.FALSE.)
  END DO
  IF ( EdgeBasis ) CALL ListRemove( SolverParams, 'Linear System Refactorize' )
  
  CALL DefaultFinish()

  ! This is now just for testing.
  i = ListGetInteger( SolverParams, 'Circulation BC', Found ) 
  IF( Found ) THEN
    BLOCK
      REAL(KIND=dp) :: Circ(2)      
      CALL  BoundaryCirculation( Circ, i, pSolver, pSolver % Variable ) 
      PRINT *,'Circulation around BC:',Circ
    END BLOCK
  END IF

  CALL Info(Caller,'All done',Level=12)
  
CONTAINS

!---------------------------------------------------------------------------------------------
  FUNCTION DoSolve() RESULT(Converged)
!---------------------------------------------------------------------------------------------
    LOGICAL :: Converged
!---------------------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: BC
    REAL(KIND=dp), POINTER CONTIG:: SavedValues(:) => NULL()
    REAL(KIND=dp) :: Norm
    INTEGER :: Active,k,n,nd,nb,t
    LOGICAL :: InitHandles 
!---------------------------------------------------------------------------------------------
    ! System assembly:
    !-----------------
    CALL Info(Caller,'Starting bulk assembly',Level=12)

    CALL DefaultInitialize()
    Active = GetNOFActive()
    InitHandles = .TRUE.
    
    DO t=1,active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes() 
      nd = GetElementNOFDOFs()  
      nb = GetElementNOFBDOFs()
      
      ! Glue local element matrix and rhs vector:
      !----------------------------------------
      CALL LocalMatrix( Element, n, nd+nb, nb, InitHandles )
    END DO
    CALL DefaultFinishBulkAssembly()

    
    ! Robin type of BC in terms of E:
    !--------------------------------
    CALL Info(Caller,'Starting boundary assembly',Level=12)
  
    Active = GetNOFBoundaryElements()
    InitHandles = .TRUE. 

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

       BLOCK
         TYPE(Element_t), POINTER :: tmp_Element
         tmp_Element => SetCurrentElement(Element)
       END BLOCK
       
       nd = GetElementNOFDOFs(Element)
       n  = GetElementNOFNodes(Element)
 
       CALL LocalMatrixBC(BC,Element,n,nd,InitHandles )
     END DO

    CALL Info(Caller,'Local assembly done',Level=12)
    
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
  
    ! Dirichlet BCs in terms of electric field E (or the potentials)
    ! --------------------------------------------------------------
    CALL DefaultDirichletBCs()

    ! Call DefaultDirichletBCs another time to apply BCs to PrecValues: 
    IF (ASSOCIATED(Solver % Matrix % PrecValues)) THEN
      SavedValues => Solver % Matrix % Values
      Solver % Matrix % Values => Solver % Matrix % PrecValues
      CALL DefaultDirichletBCs()
      Solver % Matrix % Values => SavedValues
    END IF

    CALL SingleDipoleLoad() 
    
    ! Linear system solution:
    ! -----------------------
    Norm = DefaultSolve()
    Converged = ( Solver % Variable % NonlinConverged == 1 )
!------------------------------------------------------------------------------
  END FUNCTION DoSolve
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> This subroutine sets a single dipole load.
!> It is mainly intended for testing purposes only. 
!------------------------------------------------------------------------------

  SUBROUTINE SingleDipoleLoad()

    TYPE(ValueList_t), POINTER :: DList
    REAL(KIND=dp), POINTER :: WrkArray(:,:)
    REAL(KIND=dp) :: Source(3),Coord(3),Coord0(3),EdgeVec(3),MinCoord(3)
    REAL(KIND=dp) :: MinDist,Dist, DotProd,MaxDotProd,EdgeLen,MaxEdgeLen,ParMinDist
    INTEGER :: MinInd,MaxEdgeInd,i,j,t
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Found
    
    SAVE Nodes
    
    
    DList => Solver % Values
    
    WrkArray => ListGetConstRealArray(DList,'Dipole Coordinates',Found)
    IF(.NOT. Found ) RETURN
    
    CALL Info('DipoleSource','>Dipole Coordinates< given, setting up single edge source!') 
    IF( SIZE( WrkArray, 1 ) /= 3 .OR. SIZE( WrkArray, 2 ) /= 1 ) THEN
      CALL Fatal('DipoleSource','>Dipole Coordinates< of wrong size!') 
    END IF
    Coord0(1:3) = WrkArray(1:3,1)
    
    WrkArray => ListGetConstRealArray(DList,'Dipole Current',Found)
    IF(.NOT. Found) CALL Fatal('DipoleSource','>Dipole Current< not given!') 
    IF( SIZE( WrkArray, 1 ) /= 3 .OR. SIZE( WrkArray, 2 ) /= 1 ) THEN
      CALL Fatal('DipoleSource','>Dipole Current< of wrong size!') 
    END IF
    Source(1:3) = WrkArray(1:3,1)                
    
    MinDist = HUGE( MinDist ) 
    DO t=1,Mesh % NumberOfNodes
      Coord(1) = Mesh % Nodes % x(t)
      Coord(2) = Mesh % Nodes % y(t)
      Coord(3) = Mesh % Nodes % z(t)      
      Dist = SQRT( SUM( (Coord-Coord0)**2)  )
      IF( Dist < MinDist ) THEN
        MinDist = Dist
        MinInd = t
        MinCoord = Coord
      END IF
    END DO
    
    PRINT *,'Minimum distance at node:',MinInd,MinDist
    PRINT *,'Dipole moment point:',Coord0
    PRINT *,'Nearest mesh point:',MinCoord

    IF( ParEnv % PEs > 1 ) THEN
      ParMinDist = ParallelReduction( MinDist, 1)
      IF( ParMinDist < MinDist ) RETURN
      PRINT *,'Minimum distance node found in partition:',ParEnv % MyPe
    END IF
      
    
    MaxDotProd = 0.0_dp

    DO t=1,Mesh % NumberOfEdges
      Element => Mesh % Edges(t)
            
      IF( ALL( Element % NodeIndexes /= MinInd ) ) CYCLE
      
      CALL GetElementNodes( Nodes, Element )
      
      EdgeVec(1) = Nodes % x(1) - Nodes % x(2)
      EdgeVec(2) = Nodes % y(1) - Nodes % y(2)
      EdgeVec(3) = Nodes % z(1) - Nodes % z(2)

      EdgeLen = SQRT( SUM( EdgeVec**2 ) )
      
      
      DotProd = SUM(EdgeVec*Source) / EdgeLen
      
      !PRINT *,'candidate:',t,EdgeLen,DotProd

      IF( ABS( DotProd ) > ABS( MaxDotProd ) ) THEN
        !PRINT *,'found better'
        MaxDotProd = DotProd
        MaxEdgeLen = EdgeLen
        MaxEdgeInd = t
      END IF
    END DO

    PRINT *,'Edge with best orientation:',MaxEdgeInd,maxDotProd,MaxEdgeLen
    
    i = Mesh % NumberOfNodes + MaxEdgeInd
    j = Solver % Variable % Perm(i)
    
    PRINT *,'set source at rhs:',i,j,MaxDotProd
    IF( j > 0 ) THEN
      ! Component in real valued system is: im*omega*Source
      ! i.e. real dipole source results to imaginary source
      Solver % Matrix % Rhs( 2*j ) = Omega * MaxDotProd / MaxEdgeLen
    END IF

  END SUBROUTINE SingleDipoleLoad
    
  

!-----------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd, nb, InitHandles )
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, nb
    LOGICAL :: InitHandles
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: eps, muinv, Cond, L(3)
    REAL(KIND=dp) :: DetJ, weight
    COMPLEX(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), MASS(:,:), Gauge(:,:), PREC(:,:)
    COMPLEX(KIND=dp), ALLOCATABLE, SAVE :: CurlMat(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:),WBasis(:,:),RotWBasis(:,:)
    LOGICAL :: Stat, WithNDOFs, ConductorBody
    LOGICAL :: AllocationsDone = .FALSE.
    INTEGER :: t, i, j, m, np, p, q, ndofs
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: CondCoeff_h, EpsCoeff_h, CurrDens_h, MuCoeff_h

    SAVE AllocationsDone, WBasis, RotWBasis, Basis, dBasisdx, &
        MASS, STIFF, Gauge, PREC, FORCE

    IF(.NOT. AllocationsDone ) THEN
      m = Mesh % MaxElementDOFs
      ALLOCATE( WBasis(m,3), RotWBasis(m,3), Basis(m), dBasisdx(m,3), &
          MASS(m,m), STIFF(m,m), Gauge(m,m), PREC(m,m), CurlMat(m,m), FORCE(m) )      
      AllocationsDone = .TRUE.
    END IF

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( CurrDens_h,'Body Force','Current Density', InitIm=.TRUE.,InitVec3D=.TRUE.)
      InitHandles = .FALSE.
    END IF
    
    CALL GetElementNodes( Nodes, Element )
 
    STIFF(1:nd,1:nd) = 0.0_dp
    MASS(1:nd,1:nd)  = 0.0_dp
    CurlMat = 0.0_dp
    FORCE(1:nd) = 0.0_dp

    ndofs = MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))
    np = n * ndofs

    WithNDOFs = ndofs > 0
    IF (WithNDOFs) THEN
      Gauge(1:nd,1:nd)  = 0.0_dp
    END IF

    IF (PrecMatrix) PREC = 0.0_dp
    
    ! Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree = EdgeBasisDegree)

    DO t=1,IP % n

      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 
      weight = detJ * IP % s(t)
      
      ! Compute element stiffness matrix and force vector:
      ! --------------------------------------------------     
      muinv = ListGetElementComplex( MuCoeff_h, Basis, Element, Found, GaussPoint = t )      
      IF( Found ) THEN
        muinv = muinv * mu0inv
      ELSE
        muinv = mu0inv
      END IF

      ! This is present always
      DO p = 1,nd-np
        i = p+np
        DO q = 1,nd-np
          j = q+np
          ! the mu^-1 curl E . curl v 
          CurlMat(i,j) = CurlMat(i,j) + muinv * &
              SUM(RotWBasis(p,:) * RotWBasis(q,:)) * weight
        END DO
      END DO

      ! Conductivity may also be accounted for
      Cond = ListGetElementComplex( CondCoeff_h, Basis, Element, Found, GaussPoint = t )
      ConductorBody = Found .AND. ABS(Cond) > AEPS
      IF( ConductorBody ) THEN
        DO p = 1,nd-np
          i = p+np
          DO q = 1,nd-np
            j = q+np
            ! the term i\omega\sigma E.v
            STIFF(i,j) = STIFF(i,j) - im * Omega * Cond * &
                SUM(WBasis(q,:) * WBasis(p,:)) * weight
          END DO

          IF (WithNdofs) THEN
            DO q = 1,n
              j = (q-1)*ndofs + 1
              STIFF(i,j) = STIFF(i,j) + im * Omega * Cond * &
                  SUM(dBasisdx(q,:) * WBasis(p,:)) * weight
            END DO
          END IF
        END DO
      END IF

      ! If not low frequency model, assemble the term that makes this the wave equation 
      IF(.NOT. LowFrequencyModel ) THEN
        Eps = ListGetElementComplex( EpsCoeff_h, Basis, Element, Found, GaussPoint = t )        
        IF( Found ) THEN
          Eps = Eps0 * Eps
        ELSE
          Eps = Eps0 
        END IF
          
        DO p = 1,nd-np
          i = p+np
          DO q = 1,nd-np
            j = q+np
            ! the term \omega^2 \epsilon E.v
            MASS(i,j) = MASS(i,j) - Omega**2 * Eps * &
                SUM(WBasis(q,:) * WBasis(p,:)) * weight
          END DO
        END DO
      END IF

      ! Potential current source 
      L = ListGetElementComplex3D( CurrDens_h, Basis, Element, Found, GaussPoint = t )      
      IF( Found ) THEN
        DO p = 1,nd-np
          i = p+np
          FORCE(i) = FORCE(i) + im * Omega * (SUM(L*WBasis(p,:))) * weight
        END DO
      END IF
      
      ! Additional terms related to using the Gauss law or a gauge condition
      IF (WithNDOFs) THEN
        IF (ndofs == 2) THEN
          !
          ! Using more than one nodal field is experimental and may break some
          ! generic functionality. In principle, the idea is to use the A-V representation
          ! together with a gauge constraint.
          !
          IF (ConductorBody .AND. ChargeConservation) THEN
            ! Use an implied equation about the conservation of charge
            DO p = 1,n
              i = (p-1)*ndofs + 1
              DO q = 1,nd-np
                j = q+np
                Gauge(i,j) = Gauge(i,j) - im * Omega * Cond * SUM(WBasis(q,:) * dBasisdx(i,:)) * weight
                Gauge(j,i) = Gauge(j,i) + Omega**2 * Eps * SUM(WBasis(q,:) * dBasisdx(i,:)) * weight
              END DO

              DO q = 1,n
                j = (q-1)*ndofs + 1
                Gauge(i,j) = Gauge(i,j) + im * Omega * Cond * SUM(dBasisdx(i,:) * dBasisdx(j,:)) * weight
              END DO
            END DO
          ELSE
            DO p = 1,n
              ! If two nodal DOFs per node, the first DOF is the scalar potential V related to
              ! the A-V representation. We add -w^2 div D = -w^2 rho in a weak form when
              ! E = A - grad V:
              i = (p-1)*ndofs + 1
              DO q = 1,nd-np
                j = q+np
                Gauge(i,j) = Gauge(i,j) + Omega**2 * Eps * SUM(WBasis(q,:) * dBasisdx(p,:)) * weight
                Gauge(j,i) = Gauge(j,i) + Omega**2 * Eps * SUM(WBasis(q,:) * dBasisdx(p,:)) * weight
              END DO

              DO q = 1,n
                j = (q-1)*ndofs + 1
                Gauge(i,j) = Gauge(i,j) - Omega**2 * Eps * SUM(dBasisdx(q,:) * dBasisdx(p,:)) * weight
              END DO
            END DO
          END IF
          
          DO p = 1,n
            ! The second DOF is related to the gauge condition div A =  0 (TO DO: Add
            ! Lorenz condition as an alternative)
            i = p*ndofs
            DO q = 1,nd-np
              j = q+np
              Gauge(i,j) = Gauge(i,j) - SUM(WBasis(q,:) * dBasisdx(p,:)) * weight
              Gauge(j,i) = Gauge(j,i) - SUM(WBasis(q,:) * dBasisdx(p,:)) * weight                
            END DO
          END DO
        ELSE
          IF (UseGaussLaw) THEN
            IF (ConductorBody .AND. ChargeConservation) THEN
              ! Use an implied equation about the conservation of charge
              DO i = 1,np
                DO q = 1,nd-np
                  j = q+np
                  Gauge(i,j) = Gauge(i,j) - im * Omega * Cond * SUM(WBasis(q,:) * dBasisdx(i,:)) * weight
                  Gauge(j,i) = Gauge(j,i) + Omega**2 * Eps * SUM(WBasis(q,:) * dBasisdx(i,:)) * weight
                END DO
                DO j = 1,np
                  Gauge(i,j) = Gauge(i,j) + im * Omega * Cond * SUM(dBasisdx(i,:) * dBasisdx(j,:)) * weight
                END DO
              END DO
            ELSE
              ! Add -w^2 div D = -w^2 rho in a weak form when E = A - grad V:
              DO i = 1,np
                DO q = 1,nd-np
                  j = q+np
                  Gauge(i,j) = Gauge(i,j) + Omega**2 * Eps * SUM(WBasis(q,:) * dBasisdx(i,:)) * weight
                  Gauge(j,i) = Gauge(j,i) + Omega**2 * Eps * SUM(WBasis(q,:) * dBasisdx(i,:)) * weight
                END DO
                DO j = 1,np
                  Gauge(i,j) = Gauge(i,j) - Omega**2 * Eps * SUM(dBasisdx(i,:) * dBasisdx(j,:)) * weight
                END DO
              END DO
            END IF
          ELSE
            ! Here Gauss's law is not employed, only an additional gauge constraint is introduced
            DO i = 1,np
              DO q = 1,nd-np
                j = q+np
                Gauge(i,j) = Gauge(i,j) - SUM(WBasis(q,:) * dBasisdx(i,:)) * weight
                Gauge(j,i) = Gauge(j,i) + Omega**2 * Eps * SUM(WBasis(q,:) * dBasisdx(i,:)) * weight
              END DO
              IF (LorenzCondition) THEN
                DO j = 1,np
                  Gauge(i,j) = Gauge(i,j) + Omega**2 * Eps / muinv * Basis(i) * Basis(j) * weight
                END DO
              END IF
            END DO
          END IF
        END IF
      END IF
    END DO

    IF( PrecMatrix ) THEN
      IF (CurlCurlPrec) THEN
        PREC = STIFF(1:nd,1:nd) + CurlMat(1:nd,1:nd)
        IF (WithNDOFs) THEN 
          PREC(1:nd,1:nd) = PREC(1:nd,1:nd) + Gauge(1:nd,1:nd)
        END IF
      ELSE IF (MassProportional) THEN
        PREC = -PrecDampCoeff * (MASS(1:nd,1:nd))
      ELSE
        PREC = PrecDampCoeff * (STIFF(1:nd,1:nd) + CurlMat(1:nd,1:nd) - MASS(1:nd,1:nd))
      END IF
    END IF

    STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + CurlMat(1:nd,1:nd) + MASS(1:nd, 1:nd)
    IF (WithNDOFs) THEN 
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Gauge(1:nd,1:nd)
    END IF

    IF( PrecMatrix ) THEN
      IF (.NOT. CurlCurlPrec) THEN
        PREC(1:nd,1:nd) = STIFF(1:nd,1:nd) + PREC(1:nd,1:nd)
      END IF
      IF (nb > 0) CALL CondensateP(nd-nb, nb, PREC)
      CALL DefaultUpdatePrec(PREC(1:nd,1:nd))
    END IF
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------    
    IF (nb > 0) CALL CondensateP(nd-nb, nb, STIFF, FORCE)
    CALL DefaultUpdateEquations( STIFF, FORCE, Element )

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( BC, Element, n, nd, InitHandles )
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd
    LOGICAL :: InitHandles
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), ALLOCATABLE :: STIFF(:,:), MASS(:,:), FORCE(:)    
    COMPLEX(KIND=dp) :: B, L(3), muinv, TemGrad(3), MagLoad(3), BetaPar, jn, Cond, SurfImp, epsr, mur, ep
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:),WBasis(:,:),RotWBasis(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Re_Eigenf(:), Im_Eigenf(:)
    REAL(KIND=dp) :: th, DetJ
    LOGICAL :: Stat, Found, UpdateStiff, WithNdofs, ThinSheet, ConductorBC, EigenBC, PortSource
    LOGICAL :: LineElement, DegenerateElement, Regularize, Consistent
    LOGICAL :: AllocationsDone = .FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, m, np, p, q, ndofs, EigenInd
    INTEGER :: nd_eigen
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(Element_t), POINTER :: Parent
    TYPE(ValueHandle_t), SAVE :: MagLoad_h, ElRobin_h, MuCoeff_h, EpsCoeff_h, Absorb_h, TemRe_h, TemIm_h, ExtPot_h
    TYPE(ValueHandle_t), SAVE :: TransferCoeff_h, ElCurrent_h
    TYPE(ValueHandle_t), SAVE :: Thickness_h, RelNu_h, CondCoeff_h
    TYPE(ValueHandle_t), SAVE :: GoodConductor, ChargeConservation
    TYPE(ValueHandle_t), SAVE :: EigenvectorSource, EigenvectorInd, IncidentWave
    
    SAVE AllocationsDone, WBasis, RotWBasis, Basis, dBasisdx, FORCE, STIFF, MASS, Re_Eigenf, Im_Eigenf

    IF(.NOT. AllocationsDone ) THEN
      m = Mesh % MaxElementDOFs
      ALLOCATE( WBasis(m,3), RotWBasis(m,3), Basis(m), dBasisdx(m,3),&
          FORCE(m),STIFF(m,m),MASS(m,m))      
      IF (EigenfunctionSource) THEN
        ALLOCATE(Re_Eigenf(m), Im_Eigenf(m))
      END IF
      AllocationsDone = .TRUE.
    END IF

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( ElRobin_h,'Boundary Condition','Electric Robin Coefficient',InitIm=.TRUE.)
      CALL ListInitElementKeyword( MagLoad_h,'Boundary Condition','Magnetic Boundary Load', InitIm=.TRUE.,InitVec3D=.TRUE.)
      CALL ListInitElementKeyword( Absorb_h,'Boundary Condition','Absorbing BC')
      CALL ListInitElementKeyword( GoodConductor,'Boundary Condition','Good Conductor BC')
      CALL ListInitElementKeyword( ChargeConservation,'Boundary Condition','Apply Conservation of Charge')      
      CALL ListInitElementKeyword( TemRe_h,'Boundary Condition','TEM Potential')
      CALL ListInitElementKeyword( TemIm_h,'Boundary Condition','TEM Potential Im')
      CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)      
      CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity',InitIm=.TRUE.)      

      CALL ListInitElementKeyword( TransferCoeff_h,'Boundary Condition','Electric Transfer Coefficient',InitIm=.TRUE.)
      CALL ListInitElementKeyword( ElCurrent_h,'Boundary Condition','Electric Current Density',InitIm=.TRUE.)
      CALL ListInitElementKeyword( ExtPot_h,'Boundary Condition','Incident Voltage',InitIm=.TRUE.)

      CALL ListInitElementKeyword( Thickness_h,'Boundary Condition','Layer Thickness')
      CALL ListInitElementKeyword( RelNu_h,'Boundary Condition','Layer Relative Reluctivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( CondCoeff_h,'Boundary Condition','Layer Electric Conductivity',InitIm=.TRUE.)

      CALL ListInitElementKeyword( EigenvectorSource,'Boundary Condition','Eigenfunction BC')
      CALL ListInitElementKeyword( EigenvectorInd,'Boundary Condition','Eigenfunction Index')
      CALL ListInitElementKeyword( IncidentWave,'Boundary Condition','Incident Wave')      
      InitHandles = .FALSE.
    END IF

    CALL GetElementNodes( Nodes, Element )
    
    Parent => GetBulkElementAtBoundary(Element)
    
    STIFF = 0.0_dp
    MASS = 0.0_dp
    FORCE = 0.0_dp

    ndofs = MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))
    WithNdofs = ndofs > 0
    np = n * ndofs
    
    ! Check whether BC should be created in terms of pre-computed eigenfunction:
    EigenBC = ListGetElementLogical(EigenvectorSource, Element, Found)

    IF (EigenBC) THEN
      EigenInd = ListGetElementInteger(EigenvectorInd, Element, Found)
      IF (EigenInd < 1) CALL Fatal(Caller, 'Eigenfunction Index must be positive')
      PortSource = ListGetElementLogical(IncidentWave, Element, Found)

      CALL GetScalarLocalEigenmode(Re_Eigenf, ComponentName(Eigensolver % Variable, 1), Element, &
          Eigensolver, EigenInd, ComplexPart=.FALSE.)
      CALL GetScalarLocalEigenmode(Im_Eigenf, ComponentName(Eigensolver % Variable, 2), Element, &
          Eigensolver, EigenInd, ComplexPart=.FALSE.)
      
      nd_eigen = GetElementNOFDOFs(USolver=Eigensolver)
      
      IF (WithNDOFs) THEN
        Consistent = (nd_eigen == nd)
      ELSE
        Consistent = (nd_eigen - n) == nd
      END IF
      IF (.NOT. Consistent) CALL Fatal(Caller, &
          'The DOFs of the port model are not compatible with the DOFs of this solver')
    END IF

    
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree=EdgeBasisDegree )

    IF (WithNdofs) THEN
      Regularize = UseGaussLaw .AND. ListGetElementLogical( ChargeConservation, Element, Found )
    END IF

    LineElement = GetElementFamily(Element) == 2
    DegenerateElement = (CoordinateSystemDimension() == 3) .AND. LineElement
    
    UpdateStiff = .FALSE.
    DO t=1,IP % n  
      !
      ! We need to branch as the only way to get the traces of 2D vector finite elements 
      ! is to call EdgeElementInfo:
      !
      IF (LineElement) THEN
        stat = EdgeElementInfo(Element, Nodes, IP % U(t), IP % V(t), IP % W(t), detF = detJ, &
            Basis = Basis, EdgeBasis = Wbasis, RotBasis = RotWBasis, dBasisdx = dBasisdx, &
            BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
      ELSE    
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx, &
            EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver )
      END IF

      th = ListGetElementReal(Thickness_h, Basis, Element, ThinSheet, GaussPoint = t)

      IF (DegenerateElement .AND. ThinSheet .AND. UseGaussLaw) THEN
        !
        ! If a degenerate (1-D) element, perform only simplified assembly by assuming that
        ! a given thickness is associated with the degenerate element
        !
        jn = ListGetElementComplex(ElCurrent_h, Basis, Element, Found, GaussPoint = t)

        BetaPar = ListGetElementComplex(TransferCoeff_h, Basis, Element, Found, GaussPoint = t)
        IF(Found) THEN
          ep = ListGetElementComplex(ExtPot_h, Basis, Element, Found, GaussPoint = t)
          IF(Found) jn = jn + 2 * BetaPar * ep
        END IF

        IF (ABS(BetaPar) < AEPS .AND. ABS(jn) < AEPS) CYCLE
        DO p = 1,n
          i = (p-1)*ndofs + 1
          FORCE(i) = FORCE(i) - im * omega * jn * th * Basis(p) * detJ * IP % s(t)
          DO q = 1,n
            j = (q-1)*ndofs + 1
            STIFF(i,j) = STIFF(i,j) - im * omega * BetaPar * th * Basis(p) * Basis(q) * detJ * IP % s(t)
          END DO
        END DO
        UpdateStiff = .TRUE.
        CYCLE
      END IF

      IF( ASSOCIATED( Parent ) ) THEN        
        mur = ListGetElementComplex( MuCoeff_h, Basis, Parent, Found, GaussPoint = t )      
        IF( .NOT. Found ) mur = 1.0_dp
        epsr = ListGetElementComplex( EpsCoeff_h, Basis, Parent, Found, GaussPoint = t )
        IF( .NOT. Found ) epsr = 1.0_dp
      ELSE
        epsr = 1.0_dp
        mur = 1.0_dp
      END IF      
      muinv = mur * mu0inv
      
      ConductorBC = .FALSE.
      IF (ThinSheet) THEN
        IF (ListGetElementLogical(GoodConductor, Element, Found)) &
            CALL Warn(Caller, 'Good Conductor BC neglected, given Layer Thickness used instead')
        Cond = ListGetElementComplex(CondCoeff_h, Basis, Element, Found, GaussPoint = t)
        B = th * Cond
      ELSE
        IF( ListGetElementLogical( Absorb_h, Element, Found ) ) THEN
          B = im * rob0 * SQRT( epsr / mur ) 
        ELSE
          ConductorBC = ListGetElementLogical( GoodConductor, Element, Found )
          IF (ConductorBC) THEN
            Cond = ListGetElementComplex(CondCoeff_h, Basis, Element, Found, GaussPoint = t)
            muinv = ListGetElementComplex(RelNu_h, Basis, Element, Found, GaussPoint = t)
            IF ( Found ) THEN
              muinv = muinv * mu0inv
            ELSE
              muinv = mu0inv
            END IF
            SurfImp = CMPLX(1.0_dp, -1.0_dp) * SQRT(omega/(2.0_dp * Cond * muinv))
            B = 1.0_dp/SurfImp
          ELSE
            B = ListGetElementComplex( ElRobin_h, Basis, Element, Found, GaussPoint = t )
          END IF
        END IF
      END IF

      IF (EigenBC) THEN
        B = im * SQRT(-Eigensolver % Variable % Eigenvalues(EigenInd))
        L = CMPLX(0.0_dp, 0.0_dp, kind=dp)
        IF (PortSource) THEN
          DO p=1,nd
            L(:) = L(:) + CMPLX(Re_Eigenf(n+p) * WBasis(p,:), Im_Eigenf(n+p) * WBasis(p,:), kind=dp) 
          END DO
          L = 2.0_dp * B * L
        END IF
      ELSE
        MagLoad = ListGetElementComplex3D( MagLoad_h, Basis, Element, Found, GaussPoint = t )           
        TemGrad = CMPLX( ListGetElementRealGrad( TemRe_h,dBasisdx,Element,Found), &
            ListGetElementRealGrad( TemIm_h,dBasisdx,Element,Found) )
        L = MagLoad + TemGrad
      END IF

      IF (.NOT. WithNdofs) THEN
        IF (ABS(B) < AEPS .AND. ABS(DOT_PRODUCT(L,L)) < AEPS) CYCLE
      END IF
      UpdateStiff = .TRUE.

      IF (ConductorBC .OR. ThinSheet) B = im * omega/muinv * B
      
      DO i = 1,nd-np
        p = i+np
        FORCE(p) = FORCE(p) - muinv * SUM(L*WBasis(i,:)) * detJ * IP%s(t)
        DO j = 1,nd-np
          q = j+np
          STIFF(p,q) = STIFF(p,q) - muinv * B * &
              SUM(WBasis(i,:)*WBasis(j,:)) * detJ * IP%s(t)
        END DO
      END DO

      IF (WithNdofs) THEN
        ! The following term arises if the decomposition E = A - grad V is applied:
        IF (ABS(B) > AEPS) THEN
          DO i = 1,nd-np
            p = i+np
            DO j=1,n
              q = (j-1)*ndofs + 1
              STIFF(p,q) = STIFF(p,q) + muinv * B * &
                  SUM(WBasis(i,:)*dBasisdx(j,:)) * detJ * IP%s(t)
            END DO
          END DO

          IF (Regularize) THEN
            ! Apply the conservation of surface charge (not sure whether this is beneficial):
            DO p = 1,n
              i = (p-1)*ndofs + 1
              DO q = 1,n
                j = (q-1)*ndofs + 1
                STIFF(i,j) = STIFF(i,j) + muinv * B * &
                    SUM(dBasisdx(p,:) * dBasisdx(q,:)) * detJ * IP % s(t)
              END DO

              DO q = 1,nd-np
                j = q+np
                STIFF(i,j) = STIFF(i,j) - muinv * B * &
                    SUM(dBasisdx(p,:) * WBasis(q,:)) * detJ * IP % s(t)
              END DO
            END DO
            ! TO DO: If a distribution of surface charge were also given, we would need to
            !        add an additional term in order to be consistent
          END IF
        END IF
          
        IF (UseGaussLaw) THEN
          IF (Regularize .AND. ABS(DOT_PRODUCT(L,L)) > AEPS) THEN
            ! Apply the conservation of surface charge (not sure whether this is beneficial):
            DO p = 1,n
              i = (p-1)*ndofs + 1
              FORCE(i) = FORCE(i) - muinv * SUM(L*dBasisdx(p,:)) * detJ * IP % s(t)
            END DO
          END IF
          
          jn = ListGetElementComplex(ElCurrent_h, Basis, Element, Found, GaussPoint = t)
          BetaPar = ListGetElementComplex(TransferCoeff_h, Basis, Element, Found, GaussPoint = t)
          IF( Found ) THEN
            ep = ListGetElementComplex(ExtPot_h, Basis, Element, Found, GaussPoint = t)
            IF(Found) jn = jn + 2 * BetaPar * ep
          END IF
            
          DO p = 1,n
            i = (p-1)*ndofs + 1
            FORCE(i) = FORCE(i) - im * omega * jn * Basis(p) * detJ * IP % s(t)
            DO q = 1,n
              j = (q-1)*ndofs + 1
              STIFF(i,j) = STIFF(i,j) - im * omega * BetaPar * Basis(p) * Basis(q) * detJ * IP % s(t)
            END DO
          END DO
        END IF
      END IF
    END DO

    IF (UpdateStiff) THEN
      IF (PrecMatrix) THEN
        IF (MassProportional .OR. CurlCurlPrec) THEN
          CALL DefaultUpdatePrec(STIFF)
        ELSE
          CALL DefaultUpdatePrec(PrecDampCoeff*STIFF + STIFF)
        END IF
      END IF
      CALL DefaultUpdateEquations(STIFF,FORCE,Element)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE VectorHelmholtzSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE VectorHelmholtz_Dummy(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE VectorHelmholtzUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
END SUBROUTINE VectorHelmholtz_Dummy
!------------------------------------------------------------------------------
 
 
!> \ingroup Solvers
!> Solver for computing derived fields from the electric field.
!> As the initial field is computed in H(curl) space, even the electric field
!> needs to be mapped to H0.
!------------------------------------------------------------------------------
SUBROUTINE VectorHelmholtzCalcFields_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE VectorHelmholtzUtils

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
  INTEGER :: mysolver,i,j,n,m,soln
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: Solvers(:), PSolver
  
  SolverParams => GetSolverParams()

  ! Find the solver index of the primary solver by the known procedure name.
  ! (the solver is defined here in the same module so not that dirty...)
  soln = ListGetInteger( SolverParams,'Primary Solver index', Found ) 

  IF( soln == 0 ) THEN
    DO i=1,Model % NumberOfSolvers
      sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
      j = INDEX( sname,'VectorHelmholtzSolver')
      IF( j > 0 ) THEN
        soln = i 
        EXIT
      END IF
    END DO
  END IF

  IF( soln == 0 ) THEN
    pname = GetString(SolverParams, 'Potential variable', Found)
    IF( Found ) THEN
      DO i=1,Model % NumberOfSolvers
        sname = GetString(Model % Solvers(i) % Values, 'Variable', Found)

        J=INDEX(sname,'[')-1
        IF ( j<=0 ) j=LEN_TRIM(sname)
        IF ( sname(1:j) == pname(1:LEN_TRIM(pname)) )THEN
          soln = i
          EXIT
        END IF
      END DO
    END IF
  END IF
    
  IF( soln == 0 ) THEN
    CALL Fatal('VectorHelmholtzCalcFields_Init0','Cannot locate the primary solver: '//I2S(soln))      
  ELSE
    CALL Info('VectorHelmholtzCalcFields_Init0','The primary solver index is: '//I2S(soln),Level=12)
    CALL ListAddInteger( SolverParams,'Primary Solver Index',soln ) 
  END IF

  ! If the primary solver computed DG fields, then we don't need to create DG solver on-the-fly.
  ! We only need it if we have both nodal and DG solvers needed at the same time.
  !--------------------------------------------------------------------------------------------
  IF( GetLogical(SolverParams,'Discontinuous Galerkin',Found)) RETURN

  ElementalFields = GetLogical( SolverParams, 'Calculate Elemental Fields', Found)
  IF(Found .AND. .NOT. ElementalFields) RETURN

  PSolver => Solver
  DO mysolver=1,Model % NumberOfSolvers
    IF ( ASSOCIATED(PSolver,Model % Solvers(mysolver)) ) EXIT
  END DO

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
              'VectorHelmholtz VectorHelmholtz_Dummy',.FALSE. )
  CALL ListAddString( SolverParams, 'Variable', '-nooutput cf_dummy' )

  pname = ListGetString( Model % Solvers(soln) % Values, 'Mesh', Found )
  IF(Found) THEN
    CALL ListAddString( SolverParams, 'Mesh', pname )
  END IF
  
  IF (GetLogical(Solver % Values,'Calculate Magnetic Flux Density',Found)) THEN
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable', SolverParams), &
        "Magnetic Flux Density E[Magnetic Flux Density re E:3 Magnetic Flux Density im E:3]" )
  END IF

  IF (GetLogical(Solver % Values,'Calculate Electric field',Found)) THEN
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable', SolverParams), &
        "Electric field E[Electric field re E:3 Electric field im E:3]" )
  END IF

  IF (GetLogical(Solver % Values,'Calculate Magnetic Field Strength',Found)) THEN
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable', SolverParams), &
        "Magnetic Field Strength E[Magnetic Field Strength re E:3 Magnetic Field Strength im E:3]" )
  END IF

  IF ( GetLogical( Solver % Values, 'Calculate Poynting vector', Found ) ) THEN
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable', SolverParams), &
        "Poynting vector E[Poynting vector re E:3 Poynting vector im E:3]" )
  END IF

  IF ( GetLogical( Solver % Values, 'Calculate Div of Poynting Vector', Found ) ) THEN
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable', SolverParams), &
        "Div Poynting Vector E[Div Poynting Vector re E:1 Div Poynting Vector im E:1]" )
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable', SolverParams), &
        "Joule Heating E[Joule Heating re E:1 Joule Heating im E:1]")
  END IF

  DEALLOCATE(Model % Solvers)
  Model % Solvers => Solvers
  Model % NumberOfSolvers = n+1
!------------------------------------------------------------------------------
END SUBROUTINE VectorHelmholtzCalcFields_Init0
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE VectorHelmholtzCalcFields_Init(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE VectorHelmholtzUtils

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

  CALL ListAddString( SolverParams, 'Variable', '-nooutput vectorhelmholtz_dummy' )

  CALL ListAddLogical( SolverParams, 'Linear System refactorize', .FALSE.)

  ! add these in the beginning, so that SaveScalars sees these existing, even
  ! if executed before the actual computations...
  ! -----------------------------------------------------------------------
  IF (GetLogical( SolverParams, 'Calculate Energy Functional', Found)) THEN
    CALL ListAddConstReal(Model % Simulation,'res: Energy Functional', 0._dp)
    CALL ListAddConstReal(Model % Simulation,'res: Energy Functional im', 0._dp)
  END IF

  IF (GetLogical(SolverParams,'Show Angular Frequency',Found)) &
    CALL ListAddConstReal(Model % Simulation,'res: Angular Frequency',0._dp)

  IF ( GetLogical( SolverParams, 'Calculate Electric Potential', Found ) ) THEN
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable', SolverParams), &
        "Electric Potential[Electric Potential re:1 Electric Potential im:1]")
  END IF
    
  NodalFields = GetLogical( SolverParams, 'Calculate Nodal Fields', Found)
  IF(.NOT. Found) NodalFields = .TRUE.
  
  IF(.NOT. NodalFields) RETURN

  IF (GetLogical(SolverParams,'Calculate Magnetic Flux Density',Found)) THEN
    CALL ListAddString( SolverParams,&
      NextFreeKeyword('Exported Variable', SolverParams), &
      "Magnetic Flux Density[Magnetic Flux Density re:3 Magnetic Flux Density im:3]" )
  END IF

  IF (GetLogical(SolverParams,'Calculate Electric field',Found)) THEN
    CALL ListAddString( SolverParams,&
      NextFreeKeyword('Exported Variable', SolverParams), &
      "Electric field[Electric field re:3 Electric field im:3]")
  END IF

  IF (GetLogical(SolverParams,'Calculate Magnetic Field Strength',Found)) THEN
    CALL ListAddString( SolverParams,&
      NextFreeKeyword('Exported Variable', SolverParams), &
      "Magnetic Field Strength[Magnetic Field Strength re:3 Magnetic Field Strength im:3]")
  END IF

  IF ( GetLogical( SolverParams, 'Calculate Poynting vector', Found ) ) THEN
    CALL ListAddString( SolverParams,&
      NextFreeKeyword('Exported Variable', SolverParams), &
      "Poynting vector[Poynting vector re:3 Poynting vector im:3]" )
  END IF

  IF ( GetLogical( SolverParams, 'Calculate Div of Poynting Vector', Found ) ) THEN
    CALL ListAddString( SolverParams,&
      NextFreeKeyword('Exported Variable', SolverParams), &
      "Div Poynting Vector[Div Poynting Vector re:1 Div Poynting Vector im:1]" )
    CALL ListAddString( SolverParams,&
      NextFreeKeyword('Exported Variable', SolverParams), &
      "Joule Heating[Joule Heating re:1 Joule Heating im:1]")
  END IF
  
!------------------------------------------------------------------------------
END SUBROUTINE VectorHelmholtzCalcFields_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Calculate fields resulting from the edge element formulation 
!> \ingroup Solvers
!------------------------------------------------------------------------------
 SUBROUTINE VectorHelmholtzCalcFields(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
   USE VectorHelmholtzUtils

   IMPLICIT NONE
!------------------------------------------------------------------------------
   TYPE(Solver_t) :: Solver
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: dt
   LOGICAL :: Transient
!------------------------------------------------------------------------------
   REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:), SOL(:,:), RotWBasis(:,:), Basis(:), &
       dBasisdx(:,:)
   REAL(KIND=dp) :: s,u,v,w,Norm
   REAL(KIND=dp) :: detJ, Omega, Energy, Energy_im
   ! REAL(KIND=dp) :: C_ip
   COMPLEX(KIND=dp) :: H(3), ExHc(3), PR_ip, divS, J_ip(3), &
       EdotJ, EF_ip(3), R_ip, B(3)

   TYPE(Variable_t), POINTER :: TargetVar, MFD, MFS, EF, PV, DIVPV, EW, ELPOT
   TYPE(Variable_t), POINTER :: EL_MFD, EL_MFS, EL_EF, EL_PV, EL_DIVPV, EL_EW
                              
   INTEGER :: i,j,k,l,n,nd,np,p,q,fields,efields,nfields,vDOFs,ndofs
   INTEGER :: soln, EdgeBasisDegree

   TYPE(Solver_t), POINTER :: pSolver
   REAL(KIND=dp), POINTER :: xx(:), bb(:), TempVector(:), TempRHS(:)
   REAL(KIND=dp) :: hdotE_r, hdotE_i, mu0inv, eps0

   CHARACTER(LEN=MAX_NAME_LEN) :: Pname

   LOGICAL :: Found, stat, DoAve

   TYPE(GaussIntegrationPoints_t) :: IP
   TYPE(Nodes_t), SAVE :: Nodes
   TYPE(Element_t), POINTER :: Element

   INTEGER, ALLOCATABLE :: Pivot(:)

   REAL(KIND=dp), POINTER CONTIG :: Fsave(:)
   TYPE(Mesh_t), POINTER :: Mesh
   REAL(KIND=dp), ALLOCATABLE, TARGET :: Gforce(:,:), MASS(:,:), FORCE(:,:) 

   LOGICAL :: PiolaVersion, ElementalFields, NodalFields
   LOGICAL :: UseGaussLaw, LorenzCondition
   TYPE(ValueList_t), POINTER :: SolverParams 
   TYPE(ValueHandle_t), SAVE :: EpsCoeff_h, CurrDens_h, MuCoeff_h
   ! TYPE(ValueHandle_t), SAVE :: CondCoeff_h
   CHARACTER(*), PARAMETER :: Caller = 'VectorHelmholtzCalcFields'

!-------------------------------------------------------------------------------------------

   CALL Info(Caller,'',Level=6 )
   CALL Info(Caller,'----------------------------------------------------------',Level=6 )
   CALL Info(Caller,'Computing derived fields for electromagnetic wave equation!',Level=4 )
   
   SolverParams => GetSolverParams()

   soln = ListGetInteger( SolverParams,'Primary Solver Index', Found) 
   IF( soln == 0 ) THEN
     CALL Fatal(Caller,'We should know > Primary Solver Index <')
   END IF

   ! Pointer to primary solver
   pSolver => Model % Solvers(soln)

   Pname = ListGetString( SolverParams,'Potential Variable', Found )
   IF(Found ) THEN
     TargetVar => VariableGet( pSolver % Mesh % Variables, Pname ) 
   ELSE
     TargetVar => pSolver % Variable
     Pname = getVarName(pSolver % Variable)
   END IF

   CALL Info(Caller,'Name of potential variable: '//TRIM(pName),Level=10)
     
   ! Inherit the solution basis from the primary solver
   vDOFs = TargetVar % DOFs
   IF( vDofs /= 2 ) THEN
     CALL Fatal(Caller,'Primary variable should have 2 dofs: '//I2S(vDofs))
   END IF

   CALL EdgeElementStyle(pSolver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree ) 
   IF (PiolaVersion) CALL Info(Caller,'Using Piola transformed finite elements',Level=5)

   UseGaussLaw = GetLogical(pSolver % Values, 'Use Gauss Law', Found)
   LorenzCondition = GetLogical(pSolver % Values, 'Lorenz Condition', Found)
   
   Omega = GetAngularFrequency(pSolver % Values)
   
   Found = .FALSE.
   IF( ASSOCIATED( Model % Constants ) ) THEN
     mu0inv = GetConstReal( Model % Constants,  'Permeability of Vacuum', Found )
     IF(mu0inv /= 0) mu0inv = 1/mu0inv;
   END IF
   IF(.NOT. Found ) mu0inv = 1.0_dp / ( PI * 4.0d-7 )
   
   Found = .FALSE.
   IF( ASSOCIATED( Model % Constants ) ) THEN
     eps0 = GetConstReal ( Model % Constants, 'Permittivity of Vacuum', Found )
   END IF
   IF(.NOT. Found ) eps0 = 8.854187817d-12   
   
   Mesh => GetMesh()
   
   MFD => VariableGet( Mesh % Variables, 'Magnetic Flux Density' )
   EL_MFD => VariableGet( Mesh % Variables, 'Magnetic Flux Density E' )

   MFS => VariableGet( Mesh % Variables, 'Magnetic Field Strength')
   EL_MFS => VariableGet( Mesh % Variables, 'Magnetic Field Strength E')

   EF => VariableGet( Mesh % Variables, 'Electric field')
   EL_EF => VariableGet( Mesh % Variables, 'Electric field E')

   PV => VariableGet( Mesh % Variables, 'Poynting vector')
   EL_PV => VariableGet( Mesh % Variables, 'Poynting vector E')

   DIVPV => VariableGet( Mesh % Variables, 'Div Poynting Vector')
   EL_DIVPV => VariableGet( Mesh % variables, 'Div Poynting Vector E')

   EW => VariableGet( Mesh % Variables, 'Joule Heating')
   EL_EW => VariableGet( Mesh % Variables, 'Joule Heating E')
 
   nfields = 0 
   IF ( ASSOCIATED(MFD) ) nfields=nfields+3
   IF ( ASSOCIATED(MFS) ) nfields=nfields+3
   IF ( ASSOCIATED(EF)  ) nfields=nfields+3
   IF ( ASSOCIATED(PV) ) nfields=nfields+3
   IF ( ASSOCIATED(DIVPV) ) nfields=nfields+1
   IF ( ASSOCIATED(EW) ) nfields=nfields+1
   nfields = nfields*2  ! complex problem
   NodalFields = ( nfields > 0 )

   IF(NodalFields) THEN
     ALLOCATE(GForce(SIZE(Solver % Matrix % RHS),nfields)); Gforce=0._dp
   END IF

   efields = 0 
   IF ( ASSOCIATED(EL_MFD) ) efields=efields+3
   IF ( ASSOCIATED(EL_MFS) ) efields=efields+3
   IF ( ASSOCIATED(EL_EF)  ) efields=efields+3
   IF ( ASSOCIATED(EL_PV) )  efields=efields+3
   IF ( ASSOCIATED(EL_DIVPV) ) efields=efields+1
   IF ( ASSOCIATED(EL_EW) ) efields=efields+1
   efields = efields*2 ! complex problem
   ElementalFields = ( efields > 0 ) 

   fields = MAX( efields, nfields ) 
   n = Mesh % MaxElementDOFs
   
   ALLOCATE( MASS(n,n), FORCE(n,fields), Pivot(n) )
   ALLOCATE( WBasis(n,3), SOL(2,n), RotWBasis(n,3), Basis(n), dBasisdx(n,3) )
   
   SOL = 0._dp
   Energy = 0._dp; Energy_im = 0._dp

   xx => TargetVar % Values
   
   hdotE_r = 0._dp
   hdotE_i = 0._dp

   ! This piece of code effectively does what keyword "Calculate Energy Norm" would but
   ! it treats the system complex valued.
   IF (GetLogical( SolverParams, 'Calculate Energy Functional', Found)) THEN
     IF( pSolver % Variable % dofs /= 2 ) THEN
       CALL Fatal(Caller,'Cannot compute energy norm if dofs not 2!') 
     END IF

     bb => pSolver % Matrix % RHS
     ALLOCATE(TempVector(pSolver % Matrix % NumberOfRows))
     IF ( ParEnv % PEs > 1 ) THEN
       ALLOCATE(TempRHS(SIZE(bb)))
       TempRHS = bb 
       CALL ParallelInitSolve( pSolver % Matrix, xx, TempRHS, Tempvector )
       CALL ParallelMatrixVector( pSolver % Matrix, xx, TempVector, .TRUE. )
     ELSE
       CALL MatrixVectorMultiply( pSolver % Matrix, xx, TempVector )
     END IF

     DO i = 1, size(xx,1)/2   
       IF ( ParEnv % Pes>1 ) THEN
         IF ( pSolver % Matrix % ParMatrix % ParallelInfo % &
           NeighbourList(2*(i-1)+1) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
       END IF
       hdotE_r = hdotE_r + xx(2*(i-1)+1) * Tempvector(2*(i-1)+1) - xx(2*(i-1)+2) * Tempvector(2*(i-1)+2)
       hdotE_i = hdotE_i + xx(2*(i-1)+1) * Tempvector(2*(i-1)+2) + xx(2*(i-1)+2) * Tempvector(2*(i-1)+1)
     END DO

     hdotE_r = ParallelReduction(hdotE_r)
     hdotE_i = ParallelReduction(hdotE_i)
     write (Message,*) 'Energy Functional value:', hdotE_r, hdotE_i
     CALL Info(Caller,Message)
     CALL ListAddConstReal(Model % Simulation, 'res: Energy Functional', hdotE_r)
     CALL ListAddConstReal(Model % Simulation, 'res: Energy Functional im', hdotE_i)
   END IF

   CALL DefaultInitialize()

   ! CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
   CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity',InitIm=.TRUE.)
   CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)
   CALL ListInitElementKeyword( CurrDens_h,'Body Force','Current Density', InitIm=.TRUE.,InitVec3D=.TRUE.)

   DO i = 1, GetNOFActive()
     Element => GetActiveElement(i)
     n = GetElementNOFNodes()
     ndofs = pSolver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
     np = n*ndofs
     IF (UseGaussLaw .OR. LorenzCondition) THEN
       IF (ndofs < 1) CALL Fatal(Caller, 'Nodal DOFs needed')
     END IF
     
     CALL CalcFieldsLocalAssembly()

     IF(NodalFields) THEN
       CALL DefaultUpdateEquations( MASS,FORCE(:,1))
       Fsave => Solver % Matrix % RHS
       DO l=1,fields
         Solver % Matrix % RHS => GForce(:,l)
         CALL DefaultUpdateForce(FORCE(:,l))
       END DO
       Solver % Matrix % RHS => Fsave
     END IF

     IF(ElementalFields) THEN
       nfields = 0
       CALL LUdecomp(MASS,n,pivot)
       CALL LocalSol(EL_MFD,  6, n, MASS, FORCE, pivot, nfields) ! 2*3 components
       CALL LocalSol(EL_MFS,  6, n, MASS, FORCE, pivot, nfields)
       CALL LocalSol(EL_EF,   6, n, MASS, FORCE, pivot, nfields)
       CALL LocalSol(EL_PV,   6, n, MASS, FORCE, pivot, nfields)
       CALL LocalSol(EL_DIVPV,2, n, MASS, FORCE, pivot, nfields)
       CALL LocalSol(EL_EW,   2, n, MASS, FORCE, pivot, nfields)
     END IF     
   END DO
   
   Energy = ParallelReduction(Energy)
   Energy_im = ParallelReduction(Energy_im)

   ! Assembly of the face terms:
   !----------------------------
   
   DoAve = GetLogical( SolverParams,'Average Within Materials',Found) 
   
   IF (GetLogical( SolverParams,'Discontinuous Galerkin',Found)) THEN
     IF (DoAve ) THEN
       FORCE = 0.0_dp
       CALL AddLocalFaceTerms( MASS, FORCE(:,1) )
     END IF
   END IF

   IF(NodalFields .OR. DoAve) THEN
     Fsave => NULL()
     IF(ASSOCIATED(Solver % Matrix)) Fsave => Solver % Matrix % RHS
     nfields = 0
     CALL GlobalSol(MFD,  6, Gforce, nfields, EL_MFD)
     CALL GlobalSol(MFS,  6, Gforce, nfields, EL_MFS)
     CALL GlobalSol(EF ,  6, Gforce, nfields, EL_EF)
     CALL GlobalSol(PV,   6, Gforce, nfields, EL_PV)
     CALL GlobalSol(DIVPV,2, Gforce, nfields, EL_DIVPV)
     CALL GlobalSol(EW,   2, Gforce, nfields, EL_EW)
     IF(ASSOCIATED(FSave)) Solver % Matrix % RHS => Fsave
   END IF

   WRITE(Message,*) '(Electro) Integral of Divergence of Poynting Vector: ', Energy, Energy_im
   CALL Info(Caller, Message )
   CALL ListAddConstReal(Model % Simulation,'res: Integral of Div Poynting Vector',Energy)
   CALL ListAddConstReal(Model % Simulation,'res: Integral of Div Poynting Vector im',Energy_im)

   IF (GetLogical(SolverParams,'Show Angular Frequency',Found)) THEN
     WRITE(Message,*) 'Angular Frequency: ', Omega
     CALL Info(Caller, Message )
     CALL ListAddConstReal(Model % Simulation,'res: Angular Frequency', Omega)
   END IF

   ELPOT => VariableGet( Mesh % Variables,'Electric Potential' )
   IF( ASSOCIATED( ElPot ) ) THEN
     CALL Info(Caller,'Calculating complex valued electric potential') 
     IF(.NOT. ASSOCIATED(Solver % Matrix)) THEN
       CALL Fatal(Caller,'Solver % Matrix not associated!')
     END IF
     
     CALL DefaultInitialize()
     Fsave => Solver % Matrix % RHS

     IF(.NOT. ALLOCATED(GForce) ) THEN
       ALLOCATE(GForce(SIZE(Solver % Matrix % RHS),2))
     END IF
     Gforce = 0._dp

     DO i = 1, GetNOFActive()
       Element => GetActiveElement(i)
       n = GetElementNOFNodes()
       ndofs = pSolver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1) 
       np = n*ndofs
       
       CALL ElPotLocalAssembly()

       Solver % Matrix % rhs => GForce(:,1)
       CALL DefaultUpdateEquations( MASS,FORCE(:,1))

       Solver % Matrix % RHS => GForce(:,2)
       CALL DefaultUpdateForce(FORCE(:,2))
     END DO
     
     nfields = 0
     CALL GlobalSol(ElPot, 2, Gforce, nfields )
     Solver % Matrix % RHS => Fsave     
   END IF
   
   IF(ALLOCATED(Gforce)) DEALLOCATE(Gforce)
   DEALLOCATE( MASS,FORCE,Pivot)
   
   CALL Info(Caller,'All done for now!',Level=20)
   
  
CONTAINS


  SUBROUTINE CalcFieldsLocalAssembly()

    nd = GetElementNOFDOFs(uSolver=pSolver)
    CALL GetElementNodes( Nodes )
    CALL GetVectorLocalSolution(SOL,Pname,uSolver=pSolver)

    ! Calculate nodal fields:
    ! -----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree = EdgeBasisDegree)

    MASS  = 0._dp
    FORCE = 0._dp

    ! Loop over Gaussian integration points
    !---------------------------------------
    DO j = 1,IP % n
      u = IP % U(j)
      v = IP % V(j)
      w = IP % W(j)

      stat = ElementInfo(Element,Nodes,u,v,w,detJ,Basis,dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 

      B = CMPLX(MATMUL( SOL(2,np+1:nd), RotWBasis(1:nd-np,:) ) / (Omega), &
          MATMUL( SOL(1,np+1:nd), RotWBasis(1:nd-np,:) ) / (-Omega))

      ! The conductivity as a tensor not implemented yet
      !C_ip = ListGetElementReal( CondCoeff_h, Basis, Element, Found, GaussPoint = j )

      J_ip = ListGetElementComplex3D( CurrDens_h, Basis, Element, Found, GaussPoint = j )      

      R_ip = ListGetElementComplex( MuCoeff_h, Basis, Element, Found, GaussPoint = j )      
      IF( .NOT. Found ) THEN
        R_ip = mu0inv
      ELSE
        R_ip = R_ip * mu0inv
      END IF
      H = R_ip*B

      PR_ip = ListGetElementComplex( EpsCoeff_h, Basis, Element, Found, GaussPoint = j ) 
      IF( Found ) THEN
        PR_ip = Eps0 * PR_ip
      ELSE
        PR_ip = Eps0 
      END IF

      EF_ip=CMPLX(MATMUL(SOL(1,np+1:nd),WBasis(1:nd-np,:)), MATMUL(SOL(2,np+1:nd),WBasis(1:nd-np,:)))
      IF (LorenzCondition .OR. UseGaussLaw) THEN
        DO k=1,3
          EF_ip(k) = EF_ip(k) - &
              CMPLX(SUM(SOL(1,1:np:ndofs)*dBasisdx(1:n,k)), SUM(SOL(2,1:np:ndofs)*dBasisdx(1:n,k)))
        END DO
      END IF
      ExHc = ComplexCrossProduct(EF_ip, CONJG(H))

      EdotJ = SUM(EF_ip*CONJG(J_ip))
      divS = 0.5_dp*(im * Omega * (SUM(B*CONJG(H)) + SUM(EF_ip * CONJG(PR_ip * EF_ip))) - EdotJ)

      s = IP % s(j) * detJ

      Energy = Energy + REAL(divS)*s
      Energy_im = Energy_im + AIMAG(divS)*s
      DO p=1,n
        DO q=1,n
          MASS(p,q)=MASS(p,q)+s*Basis(p)*Basis(q)
        END DO
        k = 0
        IF (ASSOCIATED(MFD) .OR. ASSOCIATED(EL_MFD)) THEN
          FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*REAL(B)*Basis(p)
          k = k+3
          FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*AIMAG(B)*Basis(p)
          k = k+3
        END IF
        IF ( ASSOCIATED(MFS).OR.ASSOCIATED(EL_MFS)) THEN
          FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*REAL(H)*Basis(p)
          k = k+3
          FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*AIMAG(H)*Basis(p)
          k = k+3
        END IF
        IF ( ASSOCIATED(EF).OR.ASSOCIATED(EL_EF)) THEN
          FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*REAL(EF_ip)*Basis(p)
          k = k+3
          FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*AIMAG(EF_ip)*Basis(p)
          k = k+3
        END IF
        IF ( ASSOCIATED(PV).OR.ASSOCIATED(EL_PV)) THEN
          FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+0.5_dp*s*REAL(ExHc)*Basis(p)
          k = k+3
          FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+0.5_dp*s*AIMAG(ExHc)*Basis(p)
          k = k+3
        END IF
        IF ( ASSOCIATED(DIVPV) .OR. ASSOCIATED(EL_DIVPV)) THEN
          FORCE(p,k+1) = FORCE(p,k+1) + s*REAL(divS)*Basis(p)
          k=k+1
          FORCE(p,k+1) = FORCE(p,k+1) + s*AIMAG(divS)*Basis(p)
          k=k+1
          FORCE(p,k+1) = FORCE(p,k+1) + 0.5_dp*s*REAL(EdotJ)*Basis(p)
          k=k+1
          FORCE(p,k+1) = FORCE(p,k+1) + 0.5_dp*s*AIMAG(EdotJ)*Basis(p)
          k=k+1
        END IF

      END DO
    END DO      

  END SUBROUTINE CalcFieldsLocalAssembly
  

  SUBROUTINE ElPotLocalAssembly()

    nd = GetElementNOFDOFs(uSolver=pSolver)
    CALL GetElementNodes( Nodes )
    CALL GetVectorLocalSolution(SOL,Pname,uSolver=pSolver)

    ! Calculate nodal fields:
    ! -----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree = EdgeBasisDegree)

    MASS  = 0._dp
    FORCE = 0._dp

    ! Loop over Gaussian integration points
    !---------------------------------------
    DO j = 1,IP % n
      u = IP % U(j)
      v = IP % V(j)
      w = IP % W(j)

      stat = ElementInfo(Element,Nodes,u,v,w,detJ,Basis,dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 

      EF_ip = CMPLX(MATMUL(SOL(1,np+1:nd),WBasis(1:nd-np,:)), MATMUL(SOL(2,np+1:nd),WBasis(1:nd-np,:)))
      IF (LorenzCondition .OR. UseGaussLaw) THEN
        DO k=1,3
          EF_ip(k) = EF_ip(k) - &
              CMPLX(SUM(SOL(1,1:np:ndofs)*dBasisdx(1:n,k)), SUM(SOL(2,1:np:ndofs)*dBasisdx(1:n,k)))
        END DO
      END IF

      s = IP % s(j) * detJ
      DO p=1,n
        DO q=1,n
          MASS(p,q) = MASS(p,q)+s * SUM(dBasisdx(p,:)*dBasisdx(q,:))
        END DO
        FORCE(p,1) = FORCE(p,1) + s * SUM(REAL(EF_ip(:)*dBasisdx(p,:)))
        FORCE(p,2) = FORCE(p,2) + s * SUM(AIMAG(EF_ip(:)*dBasisdx(p,:)))
      END DO
    END DO      

  END SUBROUTINE ElPotLocalAssembly


  
!------------------------------------------------------------------------------
 SUBROUTINE GlobalSol(Var, m, b, dofs,EL_Var )
!------------------------------------------------------------------------------
   USE MeshUtils, ONLY : CalculateBodyAverage   
   IMPLICIT NONE
   REAL(KIND=dp), TARGET CONTIG :: b(:,:)
   INTEGER :: m, dofs
   TYPE(Variable_t), POINTER :: Var
   TYPE(Variable_t), POINTER, OPTIONAL :: EL_Var
!------------------------------------------------------------------------------
   INTEGER :: i
!------------------------------------------------------------------------------

   IF(PRESENT(EL_Var)) THEN
     IF(ASSOCIATED(El_Var)) THEN
       El_Var % DgAveraged = .FALSE.
       IF( DoAve ) THEN
         CALL Info(Caller,'Averaging for field: '//TRIM(El_Var % Name),Level=10)
         CALL CalculateBodyAverage(Mesh, El_Var, .FALSE.)              
       END IF
       IF(.NOT. (ASSOCIATED(var) .AND. NodalFields) ) THEN
         dofs = dofs+m
         RETURN
       END IF
     END IF
   END IF

   IF(.NOT. ASSOCIATED(Var) ) RETURN
   
   CALL Info('VectorHelmholtz','Solving for field: '//TRIM(Var % Name),Level=6)   
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
END SUBROUTINE VectorHelmholtzCalcFields
!------------------------------------------------------------------------

