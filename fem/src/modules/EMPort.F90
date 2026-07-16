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
! *  Module for computing eigen modes from a special wave equation model
! *
! *             curl (nu curl E) - w^2 eps E + nu grad P_z = lambda * E,
! *                                - div (eps E) - eps P_z = 0
! *
! *  in a 2-D region corresponding to an electromagnetic port. Here P_z and
! *  the component of the electric field corresponding to the perpendicular
! *  direction to the plane are related in terms of the eigenvalue lambda
! *  by the equation P_z = sqrt(lambda) E_z 
! *
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: Sep 9, 2024
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Initialization of the primary solver
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE EMPortSolver_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'EMPortSolver_Init0'
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, PiolaVersion, SecondFamily, SecondOrder
   
  Params => GetSolverParams()
  
  CALL ListAddNewLogical(Params, 'Linear System Complex', .TRUE.)
  CALL ListAddNewInteger(Params, 'Variable DOFs', 2)
  CALL ListAddNewLogical(Params, 'Eigen Analysis', .TRUE.)  
  CALL ListAddNewInteger(Params, 'Nonlinear System Max Iterations', 1)

  IF (.NOT. ListCheckPresent(Params, "Element") ) THEN
    CALL EdgeElementStyle(Params, PiolaVersion, SecondFamily, SecondOrder, Check = .TRUE.)

    ! Share the DOFs definition with the vector Helmholtz model so that the solution might be
    ! utilized by the vector Helmholtz model:
    IF (SecondOrder) THEN
      CALL ListAddString(Params, "Element", &
          "n:1 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2")
    ELSE IF (PiolaVersion) THEN
      CALL ListAddString(Params, "Element", "n:1 e:1 -quad_face b:2 -quad b:2 -brick b:3")
    ELSE
      CALL ListAddString(Params, "Element", "n:1 e:1" )
    END IF
  END IF

  CALL ListAddNewString(Params, 'Variable', 'E[E re:1 E im:1]')
!-----------------------------------------------------------------------------
END SUBROUTINE EMPortSolver_Init0
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!> A special solver for finding a propagation parameter as an eigenvalue
!------------------------------------------------------------------------------
SUBROUTINE EMPortSolver(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: SolverPtr
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: PiolaVersion, EigenProblem, InitHandles, Found
  INTEGER :: DOFs, EdgeBasisDegree, Active, t, n, nd, EFamily
  COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
  REAL(KIND=dp) :: mu0inv, eps0, omega
  REAL(KIND=dp) :: Norm
  CHARACTER(*), PARAMETER :: Caller = 'EMPortSolver'
!------------------------------------------------------------------------------

  CALL Info(Caller,'',Level=8)
  CALL Info(Caller,'------------------------------------------------',Level=6)
  CALL Info(Caller,'Solving electromagnetic port equations over a surface')
  CALL Info(Caller,'------------------------------------------------',Level=6)

  IF ( CurrentCoordinateSystem() /= Cartesian ) THEN 
    CALL Fatal(Caller,'Implemented only for Cartesian problems!')
  END IF
  
  DOFs = Solver % Variable % Dofs
  IF (DOFs /= 2) THEN
    CALL Fatal(Caller, 'Complex field, specify two DOFs instead of '//I2S(DOFs))
  END IF

  SolverPtr => Solver  
  Mesh => GetMesh()
  Params => GetSolverParams()

  CALL EdgeElementStyle(Params, PiolaVersion, BasisDegree = EdgeBasisDegree )
  
  EigenProblem = EigenOrHarmonicAnalysis(Solver)
  
  ! maxiter = ListGetInteger(Params,'Nonlinear System Max Iterations', Found, minv=1)
  
  CALL DefaultStart()
  CALL InitStuff()

  CALL DefaultInitialize()

  CALL Info(Caller, 'Performing bulk element assembly', Level=12)
  Active = GetNOFActive(Solver)
  InitHandles = .TRUE.
  DO t=1,Active
    Element => GetActiveElement(t)
    EFamily = GetElementFamily(Element)
    IF (EFamily > 4) CYCLE
    
    n  = GetElementNOFNodes(Element)
    nd = GetElementNOFDOFs(Element)

    IF (EdgeBasisDegree > 1) THEN
      SELECT CASE(EFamily)    
      CASE(3)
        IF (n < 6) CALl Fatal(Caller, 'A background mesh needs 6-node triangles')
      CASE(4)
        IF (n < 9) CALl Fatal(Caller, 'A background mesh needs 9-node quads')
      END SELECT
    END IF
    
    CALL LocalMatrix(Element, n, nd, InitHandles)
  END DO
    
  CALL DefaultFinishBulkAssembly()

  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()

  CALL DefaultFinish()

  CALL Info(Caller, 'All done', Level=12)

CONTAINS


  ! Initialization of some parameters
  !--------------------------------------------------------------------
  SUBROUTINE InitStuff()

    Found = .FALSE.
    IF( ASSOCIATED( Model % Constants ) ) THEN
      IF (ListCheckPresent(Model % Constants, 'Permeability of Vacuum')) &
          mu0inv = 1.0_dp / GetConstReal( Model % Constants, 'Permeability of Vacuum', Found )
    END IF
    IF (.NOT. Found ) mu0inv = 1.0_dp / ( PI * 4.0d-7 )
    
    Found = .FALSE.
    IF( ASSOCIATED( Model % Constants ) ) THEN
      IF (ListCheckPresent(Model % Constants, 'Permittivity of Vacuum')) &
          eps0 = GetConstReal ( Model % Constants, 'Permittivity of Vacuum', Found ) 
    END IF
    IF(.NOT. Found ) eps0 = 8.854187817d-12
    
    Omega = GetAngularFrequency(Found=Found)
    IF (.NOT. Found) CALL Fatal(Caller, 'Angular frequency required')
    
  END SUBROUTINE InitStuff

    
!------------------------------------------------------------------------------
! Non-vectorized assembly of the matrix entries arising from the bulk elements
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Element, n, nd, InitHandles)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    TYPE(ValueHandle_t), SAVE :: EpsCoeff_h, NuCoeff_h
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    
    INTEGER :: m, allocstat, t
    INTEGER :: i, j, p, q, vdofs
    LOGICAL :: Stat, Found
    LOGICAL :: CreateEps, CreateNu
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:), dBasisdx(:,:), WBasis(:,:), &
        CurlWBasis(:,:)
    COMPLEX(KIND=dp), ALLOCATABLE, SAVE :: Stiff(:,:), Mass(:,:), Force(:)
    REAL(KIND=dp) :: weight, DetJ, CondAtIp
    COMPLEX(KIND=dp) :: Nu, Eps
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: RelPermit(n), RelPermit_Im(n)
    REAL(KIND=dp) :: RelNu(n), RelNu_Im(n)
!------------------------------------------------------------------------------

    IF (InitHandles) THEN
      CALL ListInitElementKeyword(NuCoeff_h, 'Material', 'Relative Reluctivity', InitIm=.TRUE.)
      CALL ListInitElementKeyword(EpsCoeff_h, 'Material', 'Relative Permittivity', InitIm=.TRUE.)
      InitHandles = .FALSE.
    END IF
    
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree = EdgeBasisDegree)
      
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      m = Mesh % MaxElementDofs
      ALLOCATE(WBasis(m,3), CurlWBasis(m,3), Basis(m), dBasisdx(m,3), Stiff(m,m), Mass(m,m), &
          Force(m), STAT=allocstat)      
      IF (allocstat /= 0) CALL Fatal(Caller, 'Local storage allocation failed')
    END IF

    CALL GetElementNodes(Nodes, Element)

    Stiff = CMPLX(0.0_dp, 0.0_dp, kind=dp)
    Mass = CMPLX(0.0_dp, 0.0_dp, kind=dp)
    Force = CMPLX(0.0_dp, 0.0_dp, kind=dp)

    ! The number of DOFs for one vector FE field  
    vdofs = nd - n

    Material => GetMaterial()
    RelPermit(1:n) = GetReal(Material, 'Relative Permittivity', CreateEps)
    RelPermit_Im(1:n) = GetReal(Material, 'Relative Permittivity Im', Found)
    CreateEps = Found .OR. CreateEps

    RelNu(1:n) = GetReal(Material, 'Relative Reluctivity', CreateNu)
    RelNu_Im(1:n) = GetReal(Material, 'Relative Reluctivity Im', Found)
    CreateNu = Found .OR. CreateNu
    
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
          detJ, Basis, dBasisdx, EdgeBasis = Wbasis, RotBasis = CurlWBasis, USolver = SolverPtr)
      Weight = IP % s(t) * DetJ

      ! NOTE: There seems to be an issue with getting the values of model parameters
      ! by using variables of type ValueHandle_t when the keyword Body Id is used to
      ! create body indentifiers. We use the traditional way to read the material parameters
      ! until this trouble is resolved.
      
!      Nu = ListGetElementComplex(NuCoeff_h, Basis, Element, CreateNu, GaussPoint = t)      
      IF( CreateNu ) THEN
        !Nu = Nu * mu0inv
        Nu = mu0inv * CMPLX(SUM(RelNu(1:n)*Basis(1:n)), SUM(RelNu_Im(1:n)*Basis(1:n)), kind=dp)
      ELSE
        Nu = mu0inv
      END IF

!      Eps = ListGetElementComplex(EpsCoeff_h, Basis, Element, CreateEps, GaussPoint = t)        
      IF( CreateEps ) THEN
        !Eps = Eps0 * Eps
        Eps = Eps0 * CMPLX(SUM(RelPermit(1:n)*Basis(1:n)), SUM(RelPermit_Im(1:n)*Basis(1:n)), kind=dp)  
      ELSE
        Eps = Eps0 
      END IF

      DO p = 1,n
        DO q = 1,n
          ! The operator -eps I for the scalar variable:
          Stiff(p,q) = Stiff(p,q) - weight * Eps * &
              Basis(p) * Basis(q)
        END DO

        ! The coupling between E_T and the scalar variable
        i = p
        DO q = 1,vdofs
          j = n + q
          Stiff(i,j) = Stiff(i,j) + Eps * SUM(WBasis(q,:) * dBasisdx(p,:)) * weight
          Stiff(j,i) = Stiff(j,i) + Nu * SUM(WBasis(q,:) * dBasisdx(p,:)) * weight
        END DO
      END DO
      
      DO p = 1,vdofs
        i = n + p
        DO q = 1,vdofs
          j = n + q
          ! The vector wave equation operator:
          Stiff(i,j) = Stiff(i,j) + weight * (Nu * SUM(CurlWBasis(q,:) * CurlWBasis(p,:)) - &
              Omega**2 * Eps * SUM(WBasis(q,:) * WBasis(p,:)))
          ! NOTE the selection of sign (to obtain a semidefinite matrix):
          Mass(i,j) = Mass(i,j) + weight * Nu * SUM(WBasis(q,:) * WBasis(p,:))
        END DO
      END DO
    END DO
    
    CALL DefaultUpdateEquations(Stiff, Force)
    CALL DefaultUpdateMass(Mass)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE EMPortSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE EMPortSolver_Post_Init0(Model, Solver, dt, Transient)
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
  LOGICAL :: Found, NodalFields, EigenAnalysis
  INTEGER :: soln, i, j
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()

  soln = 0
  DO i=1,Model % NumberOfSolvers
    sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
    j = INDEX(sname, 'EMPortSolver')
    IF( j > 0 ) THEN
      soln = i 
      EXIT
    END IF
  END DO
     
  IF (soln == 0) THEN
    CALL Fatal('EMPortSolver_post_Init0', 'Cannot locate the primary solver')      
  ELSE
    CALL Info('EMPortSolver_post_Init0', 'The primary solver index is: '//I2S(soln), Level=12)
    CALL ListAddInteger(SolverParams, 'Primary Solver Index', soln) 
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE EMPortSolver_Post_Init0
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
SUBROUTINE EMPortSolver_post_Init(Model, Solver, dt, Transient)
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
  LOGICAL :: Found, NodalFields, EigenAnalysis

  SolverParams => GetSolverParams()

!  CALL ListAddString(SolverParams, 'Variable', '-nooutput EMPortSolver_dummy' )
  CALL ListAddLogical(SolverParams, 'Linear System refactorize', .FALSE.)
  CALL ListAddNewLogical( SolverParams,'Skip Compute Nonlinear Change',.TRUE.)
  
!  NodalFields = GetLogical( SolverParams, 'Calculate Nodal Fields', Found)
!  IF (Found .AND. .NOT. NodalFields ) RETURN

!  EigenAnalysis = ListGetLogical(SolverParams,'Eigen Analysis', Found ) 
!------------------------------------------------------------------------------
END SUBROUTINE EMPortSolver_Post_Init
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!> A postprocessing solver for EMPortSolver
!------------------------------------------------------------------------------
SUBROUTINE EMPortSolver_Post(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t), POINTER :: PrimSolver
  TYPE(Variable_t), POINTER :: EF, ReVar, ImVar, Var
  TYPE(Element_t), POINTER :: Element
  TYPE(GaussIntegrationPoints_t) :: IP
  TYPE(Nodes_t), SAVE :: Nodes
  CHARACTER(*), PARAMETER :: Caller = 'EMPortSolver_Post'
  LOGICAL :: PiolaVersion, stat
  INTEGER :: soln, i, j, k, n, nd, EdgeBasisDegree, normal_ind(1)
  INTEGER :: DOFs, vdofs, p, q, ModeIndex
  INTEGER, POINTER, SAVE :: Ind(:) => NULL()
  REAL(KIND=dp), ALLOCATABLE, TARGET :: Mass(:,:), Force(:)  
  REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:), CurlWBasis(:,:), Basis(:), &
      dBasisdx(:,:)
  REAL(KIND=dp), ALLOCATABLE :: re_local_field(:), im_local_field(:)
  REAL(KIND=dp), ALLOCATABLE :: vec_local_field(:,:)
  
  REAL(KIND=dp) :: u, v, w, detJ, s
  
  REAL(KIND=dp) :: xq, ReEz, ImEz, ReE(3), ImE(3), ReV(3), ImV(3), Normal(3) 
  REAL(KIND=dp) :: Norm

  
  COMPLEX(KIND=dp) :: Beta, Ez
  
  LOGICAL :: EigenProblem, Found
  INTEGER :: Active, t
  COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
  REAL(KIND=dp) :: mu0inv, eps0, omega

!------------------------------------------------------------------------------
  DOFs = 6
  Params => GetSolverParams()

  soln = ListGetInteger(Params, 'Primary Solver Index', Found) 
  IF( soln == 0 ) THEN
    CALL Fatal(Caller, 'We should know > Primary Solver Index <')
  END IF

  Mesh => GetMesh()
  
  ! Pointer to primary solver
  PrimSolver => Model % Solvers(soln)
  CALL EdgeElementStyle(PrimSolver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree)

  n = Mesh % MaxElementDOFs   
  ALLOCATE(MASS(DOFs*n,DOFs*n), FORCE(DOFs*n))
  ALLOCATE(WBasis(n,3), CurlWBasis(n,3), Basis(n), dBasisdx(n,3) )
  ALLOCATE(Re_Local_field(n), Im_Local_field(n))

  IF (.NOT. ASSOCIATED(Ind)) ALLOCATE(Ind(n))
  
  CALL DefaultStart()
  CALL DefaultInitialize()

  ModeIndex = ListGetInteger(Params, 'Mode Index', Found)
!  print *, 'processing eigenvalue', PrimSolver % Variable % Eigenvalues(ModeIndex)

  Beta = SQRT(-PrimSolver % Variable % Eigenvalues(ModeIndex))
!  print *, 'propagation parameter beta', Beta
  
  DO k=1, GetNOFActive()
    Element => GetActiveElement(k)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs(USolver=PrimSolver)
    nd = GetElementDOFs(Ind, Element, PrimSolver)
    
    ! The number of DOFs for one vector FE field  
    vdofs = nd - n

    CALL GetElementNodes( Nodes )

    ! At the moment we assume that the wave propagates in the direction of some
    ! coordinate axis. Then the following check should be enough to get the positive
    ! direction of wave propagation: 
    Normal = NormalVector(Element, Nodes)
    normal_ind = MAXLOC(ABS(Normal))
    IF (Normal(normal_ind(1)) < 0.0_dp) Normal = -Normal
    
    CALL GetScalarLocalEigenmode(re_local_field, 'e re', Element, PrimSolver, ModeIndex, ComplexPart=.FALSE.)
    CALL GetScalarLocalEigenmode(im_local_field, 'e im', Element, PrimSolver, ModeIndex, ComplexPart=.FALSE.)
    
    Mass = 0.0_dp
    Force = 0.0_dp

    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree = EdgeBasisDegree)

    DO i=1, IP % n
      u = IP % U(i)
      v = IP % V(i)
      w = IP % W(i)

      stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis, dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = CurlWBasis, USolver = PrimSolver)

      s = IP % s(i) * detJ
      
      ReEz = SUM( Re_local_field(1:n) * Basis(1:n) )
      ImEz = SUM( Im_local_field(1:n) * Basis(1:n) )
      Ez = CMPLX(ReEz, ImEz, kind=dp) / (im * Beta)
      ReEz = REAL(Ez)
      ImEz = AIMAG(Ez)

      ReE(:) = 0.0_dp
      ImE(:) = 0.0_dp
      
      DO p=1,vdofs
        ReE(:) = ReE(:) + Re_local_field(n+p) * WBasis(p,:)
        ImE(:) = ImE(:) + Im_local_field(n+p) * WBasis(p,:)
      END DO
        
!      ReE(:) = SUM( Re_local_field(n+1:nd:2) * WBasis(1:vdofs,:) )
!      ReV(:) = SUM( Re_local_field(n+2:nd:2) * WBasis(1:vdofs,:) )
!      ImE(:) = SUM( Im_local_field(n+1:nd:2) * WBasis(1:vdofs,:) )
!      ImV(:) = SUM( Im_local_field(n+2:nd:2) * WBasis(1:vdofs,:) )

      DO j=1,DOFs
        DO p=1,n
          DO q=1,n
            Mass((p-1)*DOFs+j,(q-1)*DOFs+j) = Mass((p-1)*DOFs+j,(q-1)*DOFs+j) + s * Basis(p) * Basis(q)
          END DO
          SELECT CASE(j)
          CASE(1)
            Force((p-1)*DOFs+1) = Force((p-1)*DOFs+1) + s * (ReE(1) + Normal(1)*ReEz) * Basis(p)
          CASE(2)
            Force((p-1)*DOFs+2) = Force((p-1)*DOFs+2) + s * (ReE(2) + Normal(2)*ReEz) * Basis(p)
          CASE(3)
            Force((p-1)*DOFs+3) = Force((p-1)*DOFs+3) + s * (ReE(3) + Normal(3)*ReEz) * Basis(p)
          CASE(4)
            Force((p-1)*DOFs+4) = Force((p-1)*DOFs+4) + s * (ImE(1) + Normal(1)*ImEz) * Basis(p)
          CASE(5)
            Force((p-1)*DOFs+5) = Force((p-1)*DOFs+5) + s * (ImE(2) + Normal(2)*ImEz) * Basis(p)
          CASE(6)
            Force((p-1)*DOFs+6) = Force((p-1)*DOFs+6) + s * (ImE(3) + Normal(3)*ImEz) * Basis(p)
          END SELECT
        END DO
      END DO
    END DO

    CALL DefaultUpdateEquations(Mass, Force)
    
  END DO

  CALL DefaultFinishBulkAssembly()
  CALL DefaultFinishAssembly()
!  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()

  CALL DefaultFinish()  
!------------------------------------------------------------------------------
END SUBROUTINE EMPortSolver_Post
!------------------------------------------------------------------------------




