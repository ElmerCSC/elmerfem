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
!> Solve time-harmonic Maxwell equations using the curl-curl equation at relatively
!> high frequency using curl-conforming edge elemeing. Also low frequency model
!> available.
!> \ingroup Solvers
!-------------------------------------------------------------------------------
MODULE VectorHelmholtzUtils

   USE DefUtils

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
  LOGICAL :: Found, SecondOrder, PiolaVersion

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
    ELSE IF ( PiolaVersion ) THEN    
      CALL ListAddString( SolverParams, "Element", "n:0 e:1 -brick b:3 -quad_face b:2" )
    ELSE
      CALL ListAddString( SolverParams, "Element", "n:0 e:1" )
    END IF
  END IF

  !CALL ListAddNewLogical( SolverParams,'Hcurl Basis',.TRUE.)
  CALL ListAddNewLogical( SolverParams,'Variable Output',.FALSE.)
  CALL ListAddNewString( SolverParams,'Variable','E[E re:1 E im:1]')
  CALL ListAddNewLogical( SolverParams, "Linear System Complex", .TRUE.)

!------------------------------------------------------------------------------
END SUBROUTINE VectorHelmholtzSolver_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solve the electric field E from the rot-rot equation 
!> rot (1/mu_r) rot E + i \omega \sigma - \omega^2 epsilon_r E = i omega mu_0 J
!
!> Using edge elements (Nedelec/W basis of 1st or 2nd degree) 
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
  LOGICAL :: AllocationsDone = .FALSE., Found, HasPrecDampCoeff
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Omega, mu0inv, eps0
  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC
  INTEGER :: n,istat,i,nNodes,Active,NoIterationsMax
  TYPE(Mesh_t), POINTER :: Mesh
  COMPLEX(KIND=dp) :: PrecDampCoeff
  LOGICAL :: PiolaVersion, EdgeBasis, LowFrequencyModel
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: pSolver
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()  

  IF( GetLogical( SolverParams,'Quadratic Approximation', Found ) ) THEN
    PiolaVersion = .TRUE.
  ELSE
    PiolaVersion = GetLogical( SolverParams,'Use Piola Transform', Found )
  END IF
    
  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  Mesh => GetMesh()
  nNodes = Mesh % NumberOfNodes
  pSolver => Solver

  IF( Solver % Variable % dofs /= 2) THEN
    CALL Fatal ('VectorHelmholtzSolver', &
        'Variable is not of size two ('//TRIM(I2S(i))//'), Use: Variable = E[E re:1 E im:1]')
  ENDIF

  Omega = GetAngularFrequency(Found=Found)

  PrecDampCoeff = GetCReal(SolverParams, 'Linear System Preconditioning Damp Coefficient', HasPrecDampCoeff )
  PrecDampCoeff = CMPLX(REAL(PrecDampCoeff), &
      GetCReal(SolverParams, 'Linear System Preconditioning Damp Coefficient im', Found ) )
  HasPrecDampCoeff = HasPrecDampCoeff .OR. Found 

  Found = .FALSE.
  IF( ASSOCIATED( Model % Constants ) ) THEN
    mu0inv = GetConstReal( Model % Constants,  'Permeability of Vacuum', Found )
  END IF
  IF(.NOT. Found ) mu0inv = 1.0_dp / ( PI * 4.0d-7 )
  
  Found = .FALSE.
  IF( ASSOCIATED( Model % Constants ) ) THEN
    eps0 = GetConstReal ( Model % Constants, 'Permittivity of Vacuum', Found )
  END IF
  IF(.NOT. Found ) eps0 = 8.854187817d-12

  LowFrequencyModel = GetLogical( SolverParams,'Low Frequency Model',Found)
  
  ! Resolve internal non.linearities, if requeted:
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

  CALL Info('VectorHelmholtzSolver','All done',Level=12)
  
CONTAINS

!---------------------------------------------------------------------------------------------
  FUNCTION DoSolve() RESULT(Converged)
!---------------------------------------------------------------------------------------------
    LOGICAL :: Converged
!---------------------------------------------------------------------------------------------
    REAL(KIND=dp) :: Norm!, TOL, PrevNorm
    INTEGER :: k,n,nd,t!, i, j
    LOGICAL  :: Found, InitHandles 
!---------------------------------------------------------------------------------------------
    ! System assembly:
    !-----------------
    CALL Info('VectorHelmholtzSolver','Starting bulk assembly',Level=12)

    CALL DefaultInitialize()
    Active = GetNOFActive()
    InitHandles = .TRUE.
    
    DO t=1,active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes() 
      nd = GetElementNOFDOFs()  
      
      ! Glue local element matrix and rhs vector:
      !----------------------------------------
      CALL LocalMatrix( Element, n, nd, InitHandles )
    END DO
    CALL DefaultFinishBulkAssembly()

    
    ! Robin type of BC in terms of H:
    !--------------------------------
    CALL Info('VectorHelmholtzSolver','Starting boundary assembly',Level=12)
  
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

    CALL Info('VectorHelmholtzSolver','Local assembly done',Level=12)
    
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
  
    ! Dirichlet BCs in terms of electric field E
    ! ---------------------------------------------
    CALL DefaultDirichletBCs()


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
  SUBROUTINE LocalMatrix( Element, n, nd, InitHandles )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: InitHandles
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: eps, muinv, L(3)
    REAL(KIND=dp) :: DetJ, weight, Cond
    COMPLEX(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), MASS(:,:), DAMP(:,:), PREC(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:),WBasis(:,:),RotWBasis(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, m
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    LOGICAL :: AllocationsDone = .FALSE.
    TYPE(ValueHandle_t), SAVE :: CondCoeff_h, EpsCoeff_h, CurrCoeff_h, MuCoeff_h

    SAVE AllocationsDone, WBasis, RotWBasis, Basis, dBasisdx, &
        MASS, STIFF, DAMP, PREC, FORCE

    IF(.NOT. AllocationsDone ) THEN
      m = Mesh % MaxElementDOFs
      ALLOCATE( WBasis(m,3), RotWBasis(m,3), Basis(m), dBasisdx(m,3), &
          MASS(m,m), STIFF(m,m), DAMP(m,m), PREC(m,m), FORCE(m) )      
      AllocationsDone = .TRUE.
    END IF

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
      CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( CurrCoeff_h,'Body Force','Current Density', InitIm=.TRUE.,InitVec3D=.TRUE.)
      InitHandles = .FALSE.
    END IF
    
    CALL GetElementNodes( Nodes, Element )
 
    STIFF(1:nd,1:nd) = 0.0_dp
    MASS(1:nd,1:nd)  = 0.0_dp
    DAMP(1:nd,1:nd)  = 0.0_dp
    FORCE(1:nd) = 0.0_dp

    
    IF( HasPrecDampCoeff ) PREC = 0.0_dp    
    
    ! Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion)

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
      DO i = 1,nd
        DO j = 1,nd
          ! the mu^-1 curl u . curl v 
          STIFF(i,j) = STIFF(i,j) + muinv * &
              SUM(RotWBasis(i,:) * RotWBasis(j,:)) * weight
        END DO
      END DO

      ! Conductivity may also be accounted for
      Cond = ListGetElementReal( CondCoeff_h, Basis, Element, Found, GaussPoint = t )
      IF( Found ) THEN
        DO i = 1,nd
          DO j = 1,nd
            ! the term i\omega\sigma u.v
            STIFF(i,j) = STIFF(i,j) + im * Omega * Cond * &
                SUM(WBasis(j,:) * WBasis(i,:)) * weight
          END DO
        END DO
      END IF

      ! If not low frequency model, assembly the term that makes this the wave equation 
      IF(.NOT. LowFrequencyModel ) THEN
        Eps = ListGetElementComplex( EpsCoeff_h, Basis, Element, Found, GaussPoint = t )        
        IF( Found ) THEN
          Eps = Eps0 * Eps
        ELSE
          Eps = Eps0 
        END IF
          
        DO i = 1,nd
          DO j = 1,nd            
            ! the term \omega^2 \epsilon u.v
            MASS(i,j) = MASS(i,j) - Omega**2 * Eps * &
                SUM(WBasis(J,:) * WBasis(i,:)) * weight
          END DO
        END DO
      END IF

      ! Potential current source 
      L = ListGetElementComplex3D( CurrCoeff_h, Basis, Element, Found, GaussPoint = t )      
      IF( Found ) THEN
        DO i = 1,nd
          FORCE(i) = FORCE(i) + im * Omega * (SUM(L*WBasis(i,:))) * weight
        END DO
      END IF
                  
    END DO
    
    IF( HasPrecDampCoeff ) THEN
      PREC = PrecDampCoeff * (STIFF(1:nd,1:nd) - MASS(1:nd,1:nd))
      !PREC = PrecDampCoeff * (MASS(1:nd,1:nd))
      !CALL DefaultUpdatePrec(STIFF(1:nd,1:nd) + MASS(1:nd,1:nd) + DAMP(1:nd,1:nd))
    END IF

    STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + MASS(1:nd, 1:nd)

    IF( HasPrecDampCoeff ) THEN
      CALL DefaultUpdatePrec(STIFF(1:nd,1:nd) + PREC(1:nd,1:nd))
    END IF
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------    
    CALL DefaultUpdateEquations( STIFF, FORCE, Element )

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( BC, Element, n, nd, InitHandles )
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: InitHandles
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), ALLOCATABLE :: STIFF(:,:), MASS(:,:), FORCE(:), TemPot(:)
    COMPLEX(KIND=dp) :: B, L(3), muinv
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:),WBasis(:,:),RotWBasis(:,:)
    REAL(KIND=dp) :: DetJ, Normal(3), tanWBasis(3)
    LOGICAL :: Stat
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, m, np, p, q
    TYPE(Nodes_t), SAVE :: Nodes
    LOGICAL :: AllocationsDone = .FALSE.
    TYPE(Element_t), POINTER :: Parent
    TYPE(ValueHandle_t), SAVE :: MagLoad_h, ElRobin_h, MuCoeff_h
       
    SAVE AllocationsDone, WBasis, RotWBasis, Basis, dBasisdx, FORCE, STIFF, MASS, TemPot

    IF(.NOT. AllocationsDone ) THEN
      m = Mesh % MaxElementDOFs
      ALLOCATE( WBasis(m,3), RotWBasis(m,3), Basis(m), dBasisdx(m,3),&
          FORCE(m),STIFF(m,m),MASS(m,m),TemPot(m))      
      AllocationsDone = .TRUE.
    END IF

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( ElRobin_h,'Boundary Condition','Electric Robin Coefficient',InitIm=.TRUE.)     
      CALL ListInitElementKeyword( MagLoad_h,'Boundary Condition','Magnetic Boundary Load', InitIm=.TRUE.,InitVec3D=.TRUE.)
      CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)      
      InitHandles = .FALSE.
    END IF

    CALL GetElementNodes( Nodes, Element )
    
    Parent => GetBulkElementAtBoundary(Element)
    BC => GetBC(Element) 
    
    STIFF(1:nd,1:nd) = 0.0_dp
    MASS(1:nd,1:nd)  = 0.0_dp
    FORCE(1:nd) = 0.0_dp

    TemPot(1:n) = GetReal( BC, 'TEM Potential', Found ) + &
        im * GetReal(BC, 'TEM Potential im', Found)
          
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion)
    np = n * MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))

    DO t=1,IP % n      
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 
      
      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)

      B = ListGetElementComplex( ElRobin_h, Basis, Element, Found, GaussPoint = t )
      L = ListGetElementComplex3D( MagLoad_h, Basis, Element, Found, GaussPoint = t )

      ! The ListGetElement function does not yet work for taking derivatives
      L = L + MATMUL(TemPot(1:n), dBasisdx(1:n,1:3))
           
      muinv = mu0inv
      IF( ASSOCIATED( Parent ) ) THEN        
        muinv = ListGetElementComplex( MuCoeff_h, Basis, Parent, Found, GaussPoint = t )      
        IF( Found ) THEN
          muinv = muinv * mu0inv
        ELSE
          muinv = mu0inv
        END IF
      END IF

      DO i = 1,nd-np
        tanWBasis(:) = WBasis(i,:) - Normal * SUM(Normal* WBasis(i,:))
        p = i+np
        
        FORCE(p) = FORCE(p) - muinv * SUM(L*WBasis(i,:)) * detJ * IP%s(t)
        DO j = 1,nd-np
          q = j+np
          STIFF(p,q) = STIFF(p,q) - muinv * B * &
              SUM(tanWBasis(:)*WBasis(j,:)) * detJ * IP%s(t)
        END DO
      END DO

    END DO

   IF( HasPrecDampCoeff ) THEN
     !CALL DefaultUpdatePrec(2*STIFF)
     CALL DefaultUpdatePrec(PrecDampCoeff*STIFF + STIFF)
   END IF
   
   CALL DefaultUpdateEquations(STIFF,FORCE,Element)   
   
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE VectorHelmholtzSolver
!------------------------------------------------------------------------------

 
!> \ingroup Solvers
!> Solver for computing derived fields from the electric field.
!> As the initial field is computed in Hcurl space even the electric field
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
  INTEGER :: mysolver,i,j,k,n,m,vDOFs,soln!,l
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: Solvers(:), PSolver

  SolverParams => GetSolverParams()

  ! Find the solver index of the primary solver by the known procedure name.
  ! (the solver is defined here in the same module so not that dirty...)
  soln = 0

  DO i=1,Model % NumberOfSolvers
    sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
    j = INDEX( sname,'VectorHelmholtzSolver')
    IF( j > 0 ) THEN
      soln = i 
      EXIT
    END IF
  END DO

  
  pname = GetString(SolverParams, 'Potential variable', Found)
  IF( Found ) THEN
    IF( soln == 0 ) THEN
      DO i=1,Model % NumberOfSolvers
        sname = GetString(Model % Solvers(i) % Values, 'Variable', Found)

        J=INDEX(sname,'[')-1
        IF ( j<=0 ) j=LEN_TRIM(sname)
        IF ( sname(1:j) == pname(1:LEN_TRIM(pname)) )THEN
          soln = i
          EXIT
        END IF
      END DO
    ELSE
      CALL Info('VectorHelmholtzCalcFields_Init0',&
          'Keyword > Potential Variable < is redundant when we match the Procedure name',Level=7)
    END IF
  END IF
    
  IF( soln == 0 ) THEN
    CALL Fatal('VectorHelmholtzCalcFields_Init0','Cannot locate the primary solver: '//TRIM(I2S(soln)))      
  ELSE
    CALL Info('VectorHelmholtzCalcFields_Init0','The primary solver index is: '//TRIM(I2S(soln)),Level=12)
    CALL ListAddInteger( SolverParams,'Primary Solver Index',soln ) 
  END IF

  ! If the primary solver computed DG fields then we don't need to create DG solver on-the-fly.
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

  ! This is always two for now since the Helmholtz equation is complex valued!
  vDOFs = 2

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
              'VectorHelmholtz MagnetoDynamics_Dummy',.FALSE. )
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
        "Electric Work E[Electric Work re E:1 Electric Work im E:1]")
  END IF

  DEALLOCATE(Model % Solvers)
  Model % Solvers => Solvers
  Model % NumberOfSolvers = n+1
!------------------------------------------------------------------------------
END SUBROUTINE VectorHelmholtzCalcFields_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamics_Dummy(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE VectorHelmholtzUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
END SUBROUTINE MagnetoDynamics_Dummy
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
  LOGICAL :: Found, NodalFields!, FluxFound
  TYPE(ValueList_t), POINTER :: SolverParams!, EQ

  SolverParams => GetSolverParams()

  CALL ListAddString( SolverParams, 'Variable', '-nooutput hr_dummy' )

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

  NodalFields = GetLogical( SolverParams, 'Calculate Nodal Fields', Found)
  IF(Found.AND..NOT.NodalFields) RETURN

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
      "Electric Work[Electric Work re:1 Electric Work im:1]")
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
   REAL(KIND=dp) :: s,u,v,w,WBasis(35,3), SOL(2,35), Norm
   REAL(KIND=dp) :: RotWBasis(35,3), Basis(35), dBasisdx(35,3), E(2,3)
   REAL(KIND=dp) :: detJ, Omega, Energy, Energy_im, C_ip
   COMPLEX(KIND=dp) :: H(3), ExHc(3), PR_ip, divS, J_ip(3), PR(16), EdotJ, EF_ip(3), R_ip, &!
                       B(3), R(35)

   TYPE(Variable_t), POINTER :: MFD, MFS, EF, PV, DIVPV, EW
   TYPE(Variable_t), POINTER :: EL_MFD, EL_MFS, EL_EF, EL_PV, EL_DIVPV, EL_EW
                              
   INTEGER :: i,j,k,l,n,nd,np,p,q,dofs,edofs,ndofs,dofcount,vDOFs,dim,BodyId

   TYPE(Solver_t), POINTER :: pSolver
   REAL(KIND=dp), POINTER :: xx(:), bb(:), TempVector(:), TempRHS(:)
   REAL(KIND=dp) :: hdotE_r, hdotE_i, mu0inv, eps0

   CHARACTER(LEN=MAX_NAME_LEN) :: Pname

   TYPE(ValueList_t), POINTER :: Material, BodyForce
   LOGICAL :: Found, stat

   TYPE(GaussIntegrationPoints_t) :: IP
   TYPE(Nodes_t), SAVE :: Nodes
   TYPE(Element_t), POINTER :: Element

   INTEGER, ALLOCATABLE :: Pivot(:)

   REAL(KIND=dp), POINTER CONTIG :: Fsave(:)
   TYPE(Mesh_t), POINTER :: Mesh
   REAL(KIND=dp), ALLOCATABLE, TARGET :: Gforce(:,:), MASS(:,:), FORCE(:,:) 

   LOGICAL :: PiolaVersion, ElementalFields, NodalFields, InitHandles
   INTEGER :: soln
   TYPE(ValueList_t), POINTER :: SolverParams 
   TYPE(ValueHandle_t), SAVE :: CondCoeff_h, EpsCoeff_h, CurrCoeff_h, MuCoeff_h
 
!-------------------------------------------------------------------------------------------
   SolverParams => GetSolverParams()

   soln = ListGetInteger( SolverParams,'Primary Solver Index', Found) 
   IF( soln == 0 ) THEN
     CALL Fatal('VectorHelmholtzCalcFields','We should know > Primary Solver Index <')
   END IF

   ! Pointer to primary solver
   pSolver => Model % Solvers(soln)

   Pname = getVarName(pSolver % Variable)
   CALL Info('VectorHelmholtzCalcFields','Name of potential variable: '//TRIM(pName),Level=10)
   
   ! Inherit the solution basis from the primary solver
   vDOFs = pSolver % Variable % DOFs
   IF( vDofs /= 2 ) THEN
     CALL Fatal('VectorHelmholtzCalcFields','Primary variable should have 2 dofs: '//TRIM(I2S(vDofs)))
   END IF

   IF( GetLogical( pSolver % Values,'Quadratic Approximation', Found ) ) THEN
     PiolaVersion = .TRUE.
   ELSE
     PiolaVersion = GetLogical( pSolver % Values,'Use Piola Transform', Found )
   END IF
    
   IF (PiolaVersion) CALL Info('MagnetoDynamicsCalcFields', &
       'Using Piola transformed finite elements',Level=5)

   Omega = GetAngularFrequency(Found=Found)
   
   Found = .FALSE.
   IF( ASSOCIATED( Model % Constants ) ) THEN
     mu0inv = GetConstReal( Model % Constants,  'Permeability of Vacuum', Found )
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

   EW => VariableGet( Mesh % Variables, 'Electric Work')
   EL_EW => VariableGet( Mesh % Variables, 'Electric Work E')
 
   ndofs = 0 
   IF ( ASSOCIATED(MFD) ) ndofs=ndofs+3
   IF ( ASSOCIATED(MFS) ) ndofs=ndofs+3
   IF ( ASSOCIATED(EF)  ) ndofs=ndofs+3
   IF ( ASSOCIATED(PV) ) ndofs=ndofs+3
   IF ( ASSOCIATED(DIVPV) ) ndofs=ndofs+1
   IF ( ASSOCIATED(EW) ) ndofs=ndofs+1
   ndofs = ndofs*2  ! complex problem
   NodalFields = ( ndofs > 0 )

   IF(NodalFields) THEN
     ALLOCATE(GForce(SIZE(Solver % Matrix % RHS),ndofs)); Gforce=0._dp
   END IF

   edofs = 0 
   IF ( ASSOCIATED(EL_MFD) ) edofs=edofs+3
   IF ( ASSOCIATED(EL_MFS) ) edofs=edofs+3
   IF ( ASSOCIATED(EL_EF)  ) edofs=edofs+3
   IF ( ASSOCIATED(EL_PV) )  edofs=edofs+3
   IF ( ASSOCIATED(EL_DIVPV) ) edofs=edofs+1
   IF ( ASSOCIATED(EL_EW) ) edofs=edofs+1
   edofs = edofs*2 ! complex problem
   ElementalFields = ( edofs > 0 ) 

   dofs = MAX( edofs, ndofs ) 
   n = Mesh % MaxElementDOFs

   ALLOCATE( MASS(n,n), FORCE(n,dofs), Pivot(n) )

   SOL = 0._dp
   R=0._dp; PR=0._dp
   Energy = 0._dp; Energy_im = 0._dp

   xx => pSolver % Variable % Values
   bb => pSolver % Matrix % RHS
   
   hdotE_r = 0._dp
   hdotE_i = 0._dp

   ! This piece of code effectively does what keyword "Calculate Energy Norm" would but
   ! it treats the system complex valued.
   IF (GetLogical( SolverParams, 'Calculate Energy Functional', Found)) THEN
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
     CALL Info('VectorHelmholtzSolver',Message)
     CALL ListAddConstReal(Model % Simulation, 'res: Energy Functional', hdotE_r)
     CALL ListAddConstReal(Model % Simulation, 'res: Energy Functional im', hdotE_i)
   END IF

   CALL DefaultInitialize()

   CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
   CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity',InitIm=.TRUE.)
   CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)
   CALL ListInitElementKeyword( CurrCoeff_h,'Body Force','Current Density', InitIm=.TRUE.,InitVec3D=.TRUE.)

   DO i = 1, GetNOFActive()
     Element => GetActiveElement(i)
     n = GetElementNOFNodes()
     np = n*pSolver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
     nd = GetElementNOFDOFs(uSolver=pSolver)

     CALL GetElementNodes( Nodes )

     CALL GetVectorLocalSolution(SOL,Pname,uSolver=pSolver)

     BodyForce => GetBodyForce(Element)
     BodyId = GetBody(Element)
     Material => GetMaterial(Element)

     dim = 3

     ! Calculate nodal fields:
     ! -----------------------
     IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion)

     MASS  = 0._dp
     FORCE = 0._dp
     E = 0._dp; B=0._dp; divS=0._dp

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

       ! The conductivity as a tensor not implemnted yet
       !C_ip = ListGetElementReal( CondCoeff_h, Basis, Element, Found, GaussPoint = j )
      
       J_ip = ListGetElementComplex3D( CurrCoeff_h, Basis, Element, Found, GaussPoint = j )      
             
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

       ExHc = ComplexCrossProduct(EF_ip, CONJG(H))

       EdotJ = SUM(EF_ip*CONJG(J_ip))
       divS = 0.5_dp*(im * Omega * (SUM(B*CONJG(H)) - SUM(EF_ip * CONJG(PR_ip * EF_ip))) - EdotJ)

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
           FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*REAL(H) 
           k = k+3
           FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*AIMAG(H)
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
           FORCE(p,k+1) = FORCE(p,k+1) + s*REAL(EdotJ)*Basis(p)
           k=k+1
           FORCE(p,k+1) = FORCE(p,k+1) + s*AIMAG(EdotJ)*Basis(p)
           k=k+1
         END IF

       END DO
     END DO

     IF(NodalFields) THEN
       CALL DefaultUpdateEquations( MASS,FORCE(:,1))
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
       CALL LocalSol(EL_MFD,  6, n, MASS, FORCE, pivot, dofcount) ! 2*3 components
       CALL LocalSol(EL_MFS,  6, n, MASS, FORCE, pivot, dofcount)
       CALL LocalSol(EL_EF,   6, n, MASS, FORCE, pivot, dofcount)
       CALL LocalSol(EL_PV,   6, n, MASS, FORCE, pivot, dofcount)
       CALL LocalSol(EL_DIVPV,2, n, MASS, FORCE, pivot, dofcount)
       CALL LocalSol(EL_EW,   2, n, MASS, FORCE, pivot, dofcount)
     END IF

   END DO

   Energy = ParallelReduction(Energy)
   Energy_im = ParallelReduction(Energy_im)

    ! Assembly of the face terms:
    !----------------------------

    IF (GetLogical( SolverParams,'Discontinuous Galerkin',Found)) THEN
      IF (GetLogical( SolverParams,'Average Within Materials',Found)) THEN
        FORCE = 0.0_dp
        CALL AddLocalFaceTerms( MASS, FORCE(:,1) )
      END IF
    END IF

   IF(NodalFields) THEN
     Fsave => Solver % Matrix % RHS
     dofcount = 0
     CALL GlobalSol(MFD,  6, Gforce, dofcount)
     CALL GlobalSol(MFS,  6, Gforce, dofcount)
     CALL GlobalSol(EF ,  6, Gforce, dofcount)
     CALL GlobalSol(PV,   6, Gforce, dofcount)
     CALL GlobalSol(DIVPV,2, Gforce, dofcount)
     CALL GlobalSol(EW,   2, Gforce, dofcount)
     Solver % Matrix % RHS => Fsave
   END IF

   WRITE(Message,*) '(Electro) Integral of Divergence of Poynting Vector: ', Energy, Energy_im
   CALL Info( 'VectorHelmholtz', Message )
   CALL ListAddConstReal(Model % Simulation,'res: Integral of Div Poynting Vector',Energy)
   CALL ListAddConstReal(Model % Simulation,'res: Integral of Div Poynting Vector im',Energy_im)

   IF(ALLOCATED(Gforce)) DEALLOCATE(Gforce)
   DEALLOCATE( MASS,FORCE)
      
   IF (GetLogical(SolverParams,'Show Angular Frequency',Found)) THEN
    WRITE(Message,*) 'Angular Frequency: ', Omega
    CALL Info( 'MagnetoDynamics', Message )
    CALL ListAddConstReal(Model % Simulation,'res: Angular Frequency', Omega)
  END IF
CONTAINS

 
!------------------------------------------------------------------------------
 SUBROUTINE GlobalSol(Var, m, b, dofs )
!------------------------------------------------------------------------------
   REAL(KIND=dp), TARGET CONTIG :: b(:,:)
   INTEGER :: m, dofs
   TYPE(Variable_t), POINTER :: Var
!------------------------------------------------------------------------------
   INTEGER :: i
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
END SUBROUTINE VectorHelmholtzCalcFields
!------------------------------------------------------------------------

