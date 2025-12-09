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
! *  Module for computing the electric field from the time-harmonic wave equation
! *  by using nodal finite finite elements. Although the use of the nodal finite
! *  elements is not generally recommended for this problem, this approximation
! *  might be utilized as a preconditioner for a truthful discretization based on
! *  curl-conforming finite elements. This solver can handle the equations 
! *  either in the curl-curl form, which couples the solution
! *  components, or in the component-wise manner, which requires that
! *  the permeability is constant and the boundaries are Cartesian planes. 
! *  More flexibility in regard to setting BCs is obtained by using the curl-curl
! *  form.
! *
! *  Authors: Peter Råback + later edits by Mika Malinen
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: 24.5.2023
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Initialization of the primary solver, i.e., VectorHelmholtzNodal.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE VectorHelmholtzNodal_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'VectorHelmholtzNodal_init'
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, PrecUse, CurlCurlForm, Monolithic
  INTEGER :: soln, i, j
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
!------------------------------------------------------------------------------
  
  Params => GetSolverParams()
  !dim = CoordinateSystemDimension()

  PrecUse = ListGetLogical( Params,'Preconditioning Solver',Found )
  Monolithic = ListGetLogical( Params,'Monolithic Solver',Found )
  
  CurlCurlForm = ListGetLogical( Params,'curl-curl Form',Found )
  IF (CurlCurlForm .AND. .NOT. Monolithic) THEN
    CALL ListAddLogical(Params, 'Monolithic Solver', .TRUE.)
    Monolithic = .TRUE.
  END IF
  
  IF( Monolithic ) THEN
    ! We use different naming convention if this is used as preconditioner.
    IF( PrecUse ) THEN
      CALL ListAddNewString( Params,'Variable',&
          "Prec Elfield[PEX re:1 PEX im:1 PEY re:1 PEY im:1 PEZ re:1 PEZ im:1]" )
    ELSE
      CALL ListAddNewString( Params,'Variable',&
          "Elfield[EX re:1 EX im:1 EY re:1 EY im:1 EZ re:1 EZ im:1]" )
    END IF
  ELSE    
    ! We solve the equation component-wise. Hence the primary variable is a temporary one.
    ! TO DO: This complicates setting Dirichlet BCs as the BCs of the full vector
    !        field should be used to create a BC for the temporary variable.
    !        This cannot be done yet.
    CALL ListAddNewLogical( Params,'Variable Output',.FALSE.)
    CALL ListAddNewString( Params,'Variable','Etmp[Etmp re:1 Etmp im:1]')

    ! We use different naming convention if this is used as preconditioner. 
    IF( PrecUse ) THEN
      CALL ListAddString( Params,&
          NextFreeKeyword('Exported Variable', Params), &
          "Prec Elfield[PEX re:1 PEX im:1 PEY re:1 PEY im:1 PEZ re:1 PEZ im:1]" )
    ELSE
      CALL ListAddString( Params,&
          NextFreeKeyword('Exported Variable', Params), &
          "Elfield[EX re:1 EX im:1 EY re:1 EY im:1 EZ re:1 EZ im:1]" )
    END IF
  END IF

  CALL ListAddNewLogical( Params, "Linear System Complex", .TRUE.)  
  CALL ListAddInteger( Params,'Time Derivative Order', 0 )  
  
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
      CALL Fatal('VectorHelmholtzNodal_Init','Eigenfunction BC given without solving a port model')      
    ELSE
      CALL Info('VectorHelmholtzNodal_Init','The eigensolver index is: '//I2S(soln), Level=12)
      CALL ListAddInteger(Params, 'Eigensolver Index', soln)
    END IF
  END IF  
    
  !IF (ListGetLogical(Params,'Calculate Electric Energy',Found)) THEN
  !  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
  !      'Electric Energy Density' )
  !END IF

  !IF( ListGetLogical(Params,'Calculate Elecric Flux',Found) ) THEN
  !  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
  !      'Elecric Flux[Elecric Flux:'//I2S(dim)//']' )       
  !END IF

  ! Nodal fields that may directly be associated as nodal loads
  !IF (ListGetLogical(Params,'Calculate Nodal Energy',Found))  THEN
  !  CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
  !      'Nodal Energy Density' )
  !END IF

  
END SUBROUTINE VectorHelmholtzNodal_Init


!-----------------------------------------------------------------------------
!> A solver for the vector Helmholtz equation based on nodal basis  
!> functions
!------------------------------------------------------------------------------
SUBROUTINE VectorHelmholtzNodal( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE MeshUtils, ONLY : FollowCurvedBoundary
  USE CRSMatrix, ONLY : CRS_TransposeMatrixVectorMultiply
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  INTEGER :: i, n, nb, nd, t, active, dim, RelOrder, soln
  INTEGER :: iter, maxiter, compi, compn, compj, dofs
  LOGICAL :: Found, VecAsm, InitHandles, &
      PrecUse, PiolaVersion, SecondOrder, SecondFamily, &
      Monolithic, Segregated, CurlCurlForm, HasPrecDampCoeff, &
      EigenfunctionSource
  TYPE(ValueList_t), POINTER :: Params, EdgeSolverParams
  TYPE(Solver_t), POINTER :: Eigensolver => NULL()
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: EF, EiVar, EdgeResVar, EdgeSolVar
  REAL(KIND=dp) :: Norm(3)
  REAL(KIND=dp) :: mu0inv, eps0, rob0, omega
  COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
  COMPLEX(KIND=dp) :: PrecDampCoeff
  TYPE(Matrix_t), POINTER, SAVE :: Proj => NULL()
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
  CHARACTER(*), PARAMETER :: Caller = 'VectorHelmholtzNodal'
!------------------------------------------------------------------------------
  CALL Info(Caller,'',Level=8)
  CALL Info(Caller,'------------------------------------------------',Level=6)
  CALL Info(Caller,'Solving harmonic electric waves using nodal basis!')

  dim = CoordinateSystemDimension() 

  IF( CurrentCoordinateSystem() /= Cartesian ) THEN 
    CALL Fatal(Caller,'Only implemented for Cartesian problems!')
  END IF
  
  Mesh => GetMesh()
  Params => GetSolverParams()

  Monolithic = ListGetLogical(Params, 'Monolithic Solver', Found)
  
  EiVar => Solver % Variable
  dofs = EiVar % Dofs / 2 
  IF( dofs == 1 ) THEN
    IF (Monolithic) CALL Fatal(Caller, 'Variable DOFs incompatible with the segregated solution')
    CALL Info(Caller,'Treating the equation in segregated manner!')
    compn = dim
  ELSE IF( dofs == dim ) THEN
    IF (.NOT. Monolithic) CALL Fatal(Caller, 'Variable DOFs incompatible with the monolithic solution')
    CALL Info(Caller,'Treating the equation in monolithic manner!')
    compn = 1
  ELSE
    CALL Fatal(Caller,'Invalid number of dofs for the solver variable: '//I2S(dofs))
  END IF
  Segregated = .NOT. Monolithic

  PrecUse = ListGetLogical( Params,'Preconditioning Solver',Found ) 
  CurlCurlForm = ListGetLogical( Params,'curl-curl Form',Found )
  
  IF( PrecUse ) THEN
    IF (.NOT. Monolithic) CALL Fatal(Caller, 'The use as a preconditioner needs Monolithic Solver = True')
    
    EF => VariableGet( Mesh % Variables,'Prec ElField')        

    EdgeSolVar => NULL()
    sname = ListGetString(Params, 'Edge Update Name', Found)
    IF (Found) THEN
      EdgeSolVar => VariableGet(Mesh % Variables, sname)
      IF (.NOT. ASSOCIATED(EdgeSolVar)) CALL Fatal(Caller, 'Could not found field: '//TRIM(sname))
    ELSE
      CALL Warn(Caller, 'Give Edge Update Name to enable the use as a preconditioner')
      PrecUse = .FALSE.
    END IF
    
    EdgeResVar => NULL()  
    sname = ListGetString( Params,'Edge Residual Name',Found)
    IF(Found) THEN
      EdgeResVar => VariableGet( Mesh % Variables, sname )
      IF(.NOT. ASSOCIATED( EdgeResVar ) ) CALL Fatal(Caller,'Could not find field: '//TRIM(sname))

      EdgeSolverParams => GetSolverParams(EdgeResVar % Solver)

      CALL EdgeElementStyle(EdgeSolverParams, PiolaVersion, SecondFamily, SecondOrder, Check = .TRUE.)
      IF (SecondOrder) CALL Fatal(Caller, 'The lowest-order edge basis assumed') 

    ELSE
      CALL Warn(Caller, 'Give Edge Residual Name to enable the use as a preconditioner')
      PrecUse = .FALSE.
    END IF
    IF (PrecUse) THEN
      PrecDampCoeff = GetCReal(Params, 'Linear System Preconditioning Damp Coefficient', HasPrecDampCoeff)
      PrecDampCoeff = CMPLX(REAL(PrecDampCoeff), &
          GetCReal(Params, 'Linear System Preconditioning Damp Coefficient im', Found), kind=dp)
      HasPrecDampCoeff = HasPrecDampCoeff .OR. Found
    END IF
  ELSE
    EF => VariableGet( Mesh % Variables,'ElField')
  END IF
  
  IF(.NOT. ASSOCIATED(EF) ) THEN
    CALL Fatal(Caller,'Variable for Electric field not found!')
  END IF  

  EigenfunctionSource = ListGetLogicalAnyBC(Model, 'Eigenfunction BC')
  IF (EigenfunctionSource) THEN
    soln = ListGetInteger(Params, 'Eigensolver Index', Found) 
    IF (soln == 0) THEN
      CALL Fatal(Caller, 'We should know > Eigensolver Index <')
    END IF
    Eigensolver => Model % Solvers(soln)
  END IF
  
  IF( ListGetLogical( Params,'Follow P Curvature', Found )  ) THEN
    CALL FollowCurvedBoundary( Model, Mesh, .TRUE. ) 
  END IF
  
  CALL DefaultStart()

  maxiter = ListGetInteger( Params,'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  RelOrder = GetInteger( Params,'Relative Integration Order',Found ) 
  CALL InitStuff()

  IF (PrecUse .AND. .NOT. ASSOCIATED(Proj)) THEN
    CALL Info(Caller,'Creating projection matrix to map a nodal solution into vector element space', Level=6)
    CALL NodalToNedelecInterpolation_GlobalMatrix(Mesh, EF, EdgeSolVar, Proj, cdim=3)
  END IF
  
  DO compi=1,compn
    IF( Segregated ) THEN
      CALL Info(Caller,'Solving for component '//I2S(compi),Level=6)
    END IF
      
    CALL DefaultInitialize()

    CALL Info(Caller,'Performing bulk element assembly',Level=12)
    Active = GetNOFActive(Solver)
    InitHandles = .TRUE.
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      nb = GetElementNOFBDOFs(Element)
      CALL LocalMatrix(  Element, n, nd+nb, nb, InitHandles )
    END DO
    
    CALL DefaultFinishBulkAssembly()
    IF( PrecUse ) THEN
      !
      ! Now EdgeResVar represents the residual with respect
      ! to the basis for H(curl). We need to apply a transformation so that
      ! we may solve the residual correction equation by using the nodal basis.
      !
      IF (Monolithic) THEN
        CALL Info(Caller,'Using Transposed Projection Matrix: H(curl) -> H1', Level=6)
!          PRINT *,'sizes:',MAXVAL(Proj % Cols),MINVAL(Proj % Cols),Proj % NumberOfRows, &
!              SIZE(EdgeResVar % Values), SIZE(Solver % Matrix % rhs)
        CALL CRS_TransposeMatrixVectorMultiply(Proj, EdgeResVar % Values, Solver % Matrix % rhs )           
      ELSE
        ! TO DO: Add the transformation of the residual for the component-wise wave equation 
      END IF
    END IF

    IF( InfoActive(20) ) THEN
      CALL VectorValuesRange(Solver % Matrix % Values,SIZE(Solver % Matrix % Values),'A0')       
    END IF
        
    CALL Info(Caller,'Performing boundary element assembly',Level=12)
    Active = GetNOFBoundaryActive(Solver)

    InitHandles = .TRUE.
    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement(Element)) THEN
        n  = GetElementNOFNodes(Element)
        nd = GetElementNOFDOFs(Element)
        CALL LocalMatrixBC(  Element, n, nd, InitHandles )
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()

    IF( Monolithic ) THEN
      CALL DefaultDirichletBCs()
    ELSE
      DO i=1,2
        sname = ComponentName( EF,2*(compi-1)+i) 
        CALL SetDirichletBoundaries( CurrentModel, Solver % Matrix, Solver % Matrix % rhs, &
            sname, i, 2, Solver % Variable % Perm )
      END DO
      CALL EnforceDirichletConditions( Solver, Solver % Matrix, Solver % Matrix % rhs )
    END IF
          
    ! And finally, solve:
    !--------------------
    IF( Segregated) THEN
      EiVar % Values(1::2) = EF % Values(2*compi-1::2*dim) 
      EiVar % Values(2::2) = EF % Values(2*compi::2*dim) 
    END IF
      
    Norm(compi) = DefaultSolve()

    IF( Segregated ) THEN
      IF( InfoActive(25) ) THEN
        CALL VectorValuesRange(EiVar % Values,SIZE(EiVar % Values),'E'//I2S(compi))       
        PRINT *,'Component Norm:',Norm(compi)
      END IF
      EF % Values(2*compi-1::2*dim) = EiVar % Values(1::2)
      EF % Values(2*compi::2*dim) = EiVar % Values(2::2)
    END IF
  END DO ! compi

  IF (PrecUse) THEN
    CALL Info(Caller,'Projecting nodal solution to vector element space', Level=6)
    CALL Info(Caller,'Using Projection Matrix: H1 -> H(curl)',Level=6)
    CALL CRS_MatrixVectorMultiply(Proj, EF % Values, EdgeSolVar % Values ) 
  END IF
    
  !IF( Solver % Variable % NonlinConverged == 1 ) EXIT
  
  CALL DefaultFinish()

  IF( Segregated ) THEN
    Solver % Variable % Norm = SQRT(SUM(Norm(1:compn)**2) / compn)
  END IF
    
  CALL Info(Caller,'All done',Level=12)
  
    
CONTAINS


  ! Initialize some parameters.
  !--------------------------------------------------------------------
  SUBROUTINE InitStuff()

    Found = .FALSE.
    IF( ASSOCIATED( Model % Constants ) ) THEN
      mu0inv = 1.0_dp / GetConstReal( Model % Constants,  'Permeability of Vacuum', Found )
    END IF
    IF(.NOT. Found ) mu0inv = 1.0_dp / ( PI * 4.0d-7 )
    
    Found = .FALSE.
    IF( ASSOCIATED( Model % Constants ) ) THEN
      eps0 = GetConstReal ( Model % Constants, 'Permittivity of Vacuum', Found )
    END IF
    IF(.NOT. Found ) eps0 = 8.854187817d-12
    
    Omega = GetAngularFrequency(Found=Found)
    IF(.NOT. Found) CALL Fatal(Caller,'We should have Omega!')

    rob0 = Omega * SQRT( eps0 / mu0inv )

    !PRINT *,'InitStuff:',mu0inv, eps0, omega, rob0
    
  END SUBROUTINE InitStuff

    
!------------------------------------------------------------------------------
! Assembly of the matrix entries arising from the bulk elements. Not vectorized.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd, nb, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:),dBasisdx(:,:)
    COMPLEX(KIND=dp), ALLOCATABLE, SAVE :: STIFF(:,:,:), FORCE(:,:)
    REAL(KIND=dp) :: weight, DetJ, CondAtIp
    COMPLEX(KIND=dp) :: muinvAtIp, EpsAtIp, CurrAtIp(3), kappa
    LOGICAL :: Stat,Found, WithConductivity
    INTEGER :: i,j,k,t,p,q,m,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: CondCoeff_h, EpsCoeff_h, MuCoeff_h, CurrDens_h
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
      CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( CurrDens_h,'Body Force','Current Density',InitIm=.TRUE.,InitVec3D=.TRUE.)
      InitHandles = .FALSE.
    END IF
    
    IF( RelOrder /= 0 ) THEN
      IP = GaussPoints( Element, RelOrder = RelOrder)
    ELSE
      IP = GaussPoints( Element )
    END IF
      
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      m = Mesh % MaxElementDofs
      IF (CurlCurlForm) THEN
        k = 1
        q = 3*m
      ELSE
        k = 3
        q = m
      END IF
      ALLOCATE(Basis(m), dBasisdx(m,3), STIFF(q,q,k), FORCE(q,k), STAT=allocstat)      
      IF (allocstat /= 0) CALL Fatal(Caller,'Local storage allocation failed')
    END IF

    CALL GetElementNodes( Nodes, UElement=Element )

    kappa = CMPLX(1.0_dp, 0.0_dp, kind=dp)
    IF (PrecUse) THEN
      IF (HasPrecDampCoeff) kappa = kappa - PrecDampCoeff
    END IF
      
    ! Initialize
    STIFF = 0._dp
    FORCE = 0._dp
    
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      Weight = IP % s(t) * DetJ

      muinvAtIp = ListGetElementComplex( MuCoeff_h, Basis, Element, Found, GaussPoint = t )      
      IF( Found ) THEN
        muinvAtIp = muinvAtIp * mu0inv
      ELSE
        muinvAtIp = mu0inv
      END IF
      
      EpsAtIp = ListGetElementComplex( EpsCoeff_h, Basis, Element, Found, GaussPoint = t )      
      IF( Found ) THEN
        EpsAtIp = EpsAtIp * eps0
      ELSE
        epsAtIp = eps0
      END IF        

      CondAtIp = ListGetElementReal( CondCoeff_h, Basis, Element, WithConductivity, GaussPoint = t )
      
      IF (CurlCurlForm) THEN
        DO p=1,nd
          DO j=1,dim
            DO q=1,nd
              STIFF(dim*(p-1)+j,dim*(q-1)+j,1) = STIFF(dim*(p-1)+j,dim*(q-1)+j,1) - &
                  Weight * kappa * Omega**2 * epsAtIP * Basis(q) * Basis(p)

              IF (WithConductivity) THEN
                STIFF(dim*(p-1)+j,dim*(q-1)+j,1) = STIFF(dim*(p-1)+j,dim*(q-1)+j,1) - &
                    im * Weight * Omega * CondAtIP * Basis(q) * Basis(p)
              END IF
              
              IF (j==1) THEN
                STIFF(dim*(p-1)+j,dim*(q-1)+1,1) = STIFF(dim*(p-1)+j,dim*(q-1)+1,1) + Weight * MuinvAtIp * ( &
                    dBasisdx(q,3) * dBasisdx(p,3) + dBasisdx(q,2) * dBasisdx(p,2) )
                STIFF(dim*(p-1)+j,dim*(q-1)+2,1) = STIFF(dim*(p-1)+j,dim*(q-1)+2,1) - Weight * MuinvAtIp * &
                    dBasisdx(q,1) * dBasisdx(p,2)
                STIFF(dim*(p-1)+j,dim*(q-1)+3,1) = STIFF(dim*(p-1)+j,dim*(q-1)+3,1) - Weight * MuinvAtIp * &
                    dBasisdx(q,1) * dBasisdx(p,3)
              END IF

              IF (j==2) THEN
                STIFF(dim*(p-1)+j,dim*(q-1)+1,1) = STIFF(dim*(p-1)+j,dim*(q-1)+1,1) - Weight * MuinvAtIp * &
                    dBasisdx(q,2) * dBasisdx(p,1)
                STIFF(dim*(p-1)+j,dim*(q-1)+2,1) = STIFF(dim*(p-1)+j,dim*(q-1)+2,1) + Weight * MuinvAtIp * ( &
                    dBasisdx(q,3) * dBasisdx(p,3) + dBasisdx(q,1) * dBasisdx(p,1) )
                STIFF(dim*(p-1)+j,dim*(q-1)+3,1) = STIFF(dim*(p-1)+j,dim*(q-1)+3,1) - Weight * MuinvAtIp * &
                    dBasisdx(q,2) * dBasisdx(p,3)
              END IF

              IF (j==3) THEN
                STIFF(dim*(p-1)+j,dim*(q-1)+1,1) = STIFF(dim*(p-1)+j,dim*(q-1)+1,1) - Weight * MuinvAtIp * &
                    dBasisdx(q,3) * dBasisdx(p,1)
                STIFF(dim*(p-1)+j,dim*(q-1)+2,1) = STIFF(dim*(p-1)+j,dim*(q-1)+2,1) - Weight * MuinvAtIp * &
                    dBasisdx(q,3) * dBasisdx(p,2)
                STIFF(dim*(p-1)+j,dim*(q-1)+3,1) = STIFF(dim*(p-1)+j,dim*(q-1)+3,1) + Weight * MuinvAtIp * ( &
                    dBasisdx(q,2) * dBasisdx(p,2) + dBasisdx(q,1) * dBasisdx(p,1) )
              END IF

            END DO
          END DO
        END DO
        ! TO DO: add the integration of RHS if not used as a preconditioner
      ELSE
        
        ! "diffusion" term (D*grad(u),grad(v)):
        ! -----------------------------------      
        STIFF(1:nd,1:nd,1) = STIFF(1:nd,1:nd,1) + Weight * &
            MuinvAtIp * MATMUL( dBasisdx(1:nd,:), TRANSPOSE( dBasisdx(1:nd,:) ) )

        IF( WithConductivity ) THEN
          DO p=1,nd
            STIFF(p,1:nd,1) = STIFF(p,1:nd,1) - im * Weight * Omega * CondAtIP * Basis(p) * Basis(1:nd)
          END DO
        END IF

        ! This is the same for each component with isotropic materials!
        DO p=1,nd
          STIFF(p,1:nd,1) = STIFF(p,1:nd,1) - Weight * kappa * Omega**2 * epsAtIP * Basis(p) * Basis(1:nd)
        END DO

        IF(.NOT. PrecUse ) THEN
          CurrAtIP = ListGetElementComplex3D( CurrDens_h, Basis, Element, Found )                 
          IF( Found ) THEN
            IF( Monolithic ) THEN
              DO i=1,dofs
                FORCE(1:nd,i) = FORCE(1:nd,i) + Weight * CurrAtIp(i) * Basis(1:nd)
              END DO
            ELSE
              FORCE(1:nd,1) = FORCE(1:nd,1) + Weight * CurrAtIP(compi) * Basis(1:nd)
            END IF
          END IF
        END IF
      END IF
    END DO

    IF (CurlCurlForm) THEN
      CALL DefaultUpdateEquations(STIFF(1:dim*nd,1:dim*nd,1),FORCE(1:dim*nd,1),UElement=Element)
    ELSE
      IF( Monolithic ) THEN
        DO i=2,dofs
          STIFF(1:nd,1:nd,i) = STIFF(1:nd,1:nd,1)
        END DO
        ! Note that we use diagonal form for this!
        CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element)
      ELSE
        CALL DefaultUpdateEquations(STIFF(:,:,1),FORCE(:,1),UElement=Element)
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


! Assembly of the matrix entries arising from the Neumann and Robin conditions.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: F,C,Ext, Weight, coeff
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), DetJ,Coord(3),Normal(3)
    REAL(KIND=dp) :: TestVec(3), TrialVec(3)
!    COMPLEX(KIND=dp) :: STIFF(nd,nd,3), FORCE(nd,3)
    COMPLEX(KIND=dp), ALLOCATABLE, SAVE :: STIFF(:,:,:), FORCE(:,:)
    COMPLEX(KIND=dp) :: muInvAtIp, muinv, Cond, SurfImp, TemGrad(3), L(3), B 
    LOGICAL :: Stat,Found,RobinBC,NT,EigenBC,GoodConductor
    INTEGER :: i,j,k,m,p,q,t,allocstat,EigenInd
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC       
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(Element_t), POINTER :: Parent
    TYPE(ValueHandle_t), SAVE :: ElRobin_h, MagLoad_h, Absorb_h, TemRe_h, TemIm_h, MuCoeff_h
    TYPE(ValueHandle_t), SAVE :: GoodConductor_h, EigenvectorSource, EigenvectorInd
    TYPE(ValueHandle_t), SAVE :: RelNu_h, CondCoeff_h
    
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( ElRobin_h,'Boundary Condition','Electric Robin Coefficient',InitIm=.TRUE.)
      CALL ListInitElementKeyword( MagLoad_h,'Boundary Condition','Magnetic Boundary Load', InitIm=.TRUE.,InitVec3D=.TRUE.)
      CALL ListInitElementKeyword( Absorb_h,'Boundary Condition','Absorbing BC')
      CALL ListInitElementKeyword( GoodConductor_h,'Boundary Condition','Good Conductor BC')
      CALL ListInitElementKeyword( TemRe_h,'Boundary Condition','TEM Potential')
      CALL ListInitElementKeyword( TemIm_h,'Boundary Condition','TEM Potential Im')
      CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( CondCoeff_h,'Boundary Condition','Layer Electric Conductivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( RelNu_h,'Boundary Condition','Layer Relative Reluctivity',InitIm=.TRUE.)
      CALL ListInitElementKeyword( EigenvectorSource,'Boundary Condition','Eigenfunction BC')
      CALL ListInitElementKeyword( EigenvectorInd,'Boundary Condition','Eigenfunction Index')
      InitHandles = .FALSE.
    END IF
    
    CALL GetElementNodes( Nodes, UElement=Element )

    IF (.NOT. ALLOCATED(STIFF)) THEN
      m = Mesh % MaxElementDofs
      IF (CurlCurlForm) THEN
        k = 1
        q = 3*m
      ELSE
        k = 3
        q = m
      END IF
      ALLOCATE(STIFF(q,q,k), FORCE(q,k), STAT=allocstat)      
      IF (allocstat /= 0) CALL Fatal(Caller,'Local storage allocation failed')
    END IF
    STIFF = 0._dp
    FORCE = 0._dp

    IF (.NOT. CurlCurlForm) Normal = NormalVector( Element, Nodes )
    NT = ListGetLogical( BC,'Normal-Tangential '//GetVarName(Solver % Variable), Found )
    IF(NT .AND. .NOT. Monolithic) THEN
      CALL Fatal(Caller,'Normal-tangential conditions require monolithic solver!')
    END IF

    GoodConductor = ListGetElementLogical(GoodConductor_h, Element, Found)
    
    ! Check whether BC should be created in terms of pre-computed eigenfunction:
    EigenBC = ListGetElementLogical(EigenvectorSource, Element, Found)
    IF (EigenBC) THEN
      EigenInd = ListGetElementInteger(EigenvectorInd, Element, Found)
      IF (EigenInd < 1) CALL Fatal(Caller, 'Eigenfunction Index must be positive')
    END IF

    
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    
    Parent => GetBulkElementAtBoundary(Element)
    
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )
      Weight = IP % s(t) * DetJ

      Found = .FALSE.
      IF( ASSOCIATED( Parent ) ) THEN        
        muinvAtIp = ListGetElementComplex( MuCoeff_h, Basis, Parent, Found, GaussPoint = t )      
      END IF
      IF( Found ) THEN
        muinvAtIp = muinvAtIp * mu0inv
      ELSE
        muinvAtIp = mu0inv
      END IF

      IF( .NOT. PrecUse ) THEN
        L = ListGetElementComplex3D( MagLoad_h, Basis, Element, Found, GaussPoint = t )
        TemGrad = CMPLX( ListGetElementRealGrad( TemRe_h,dBasisdx,Element,Found), &
            ListGetElementRealGrad( TemIm_h,dBasisdx,Element,Found),KIND=dp )
        L = L + TemGrad
        DO i=1,dim
          FORCE(1:nd,i) = FORCE(1:nd,i) - muinvAtIp * L(i) * Basis(1:nd) * Weight
        END DO
      END IF

      IF (EigenBC) THEN
        B = CMPLX(0.0_dp, 1.0_dp, kind=dp) * SQRT(-Eigensolver % Variable % Eigenvalues(EigenInd))
        Found = .TRUE.
      ELSE
        IF( ListGetElementLogical( Absorb_h, Element, Found ) ) THEN
          B = CMPLX(0.0_dp, rob0, KIND=dp ) 
        ELSE IF (GoodConductor) THEN
          Cond = ListGetElementComplex(CondCoeff_h, Basis, Element, Found, GaussPoint = t)
          muinv = ListGetElementComplex(RelNu_h, Basis, Element, Found, GaussPoint = t)
          IF ( Found ) THEN
            muinv = muinv * mu0inv
          ELSE
            muinv = mu0inv
          END IF
          SurfImp = CMPLX(1.0_dp, -1.0_dp, KIND=dp) * SQRT(omega/(2.0_dp * Cond * muinv))
          B = 1.0_dp/SurfImp    
        ELSE
          B = ListGetElementComplex( ElRobin_h, Basis, Element, Found, GaussPoint = t )
        END IF
      END IF

      IF( Found ) THEN
        IF (CurlCurlForm) THEN
          Normal = Normalvector(Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
          DO p=1,nd
            DO j=1,dim
              TestVec = 0.0d0
              TestVec(j) = Basis(p)
              TestVec = CrossProduct(TestVec, Normal)
              DO q=1,nd
                DO i=1,dim
                  TrialVec = 0.0d0
                  TrialVec(i) = Basis(q)
                  TrialVec = CrossProduct(TrialVec, Normal)

                  STIFF(dim*(p-1)+j, dim*(q-1)+i,1) = STIFF(dim*(p-1)+j, dim*(q-1)+i,1) - &
                      muinvAtIp * B * SUM(TestVec(:) * TrialVec(:)) * Weight 
                END DO
              END DO
            END DO
          END DO
        ELSE
          DO i=1,dim
            IF( NT ) THEN
              IF(i==1) CYCLE
              coeff = 1.0_dp
            ELSE          
              coeff = 1.0_dp - Normal(i)**2
              IF(coeff <= 0.0_dp ) THEN
                coeff = 0.0_dp
              ELSE
                coeff = SQRT(coeff)
              END IF
            END IF
            DO p = 1,nd
              STIFF(p,1:nd,i) = STIFF(p,1:nd,i) - coeff * muinvAtIp * B * &
                  Basis(p) * Basis(1:nd) * Weight
            END DO
          END DO
        END IF
      END IF
    END DO

    IF (CurlCurlForm) THEN
      CALL DefaultUpdateEquations(STIFF(1:dim*nd,1:dim*nd,1),FORCE(1:dim*nd,1),UElement=Element)      
    ELSE
      IF( Monolithic ) THEN
        ! For normal-tangential coordinate system the slip 
        CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element)
      ELSE
        CALL DefaultUpdateEquations(STIFF(:,:,compi),FORCE(:,compi),UElement=Element)
      END IF
    END IF
      
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE VectorHelmholtzNodal
!------------------------------------------------------------------------------
