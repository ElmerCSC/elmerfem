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
! *  Module for computing component-wise the electric field from the wave
! *  equation by using nodal finite finite elements. This only works when
! *  the permeability is constant and the boundaries are Cartesian. 
! *
! *  Authors: Peter Råback
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
  LOGICAL :: Found, PrecUse, Monolithic
  !INTEGER :: dim
   
  Params => GetSolverParams()
  !dim = CoordinateSystemDimension()

  PrecUse = ListGetLogical( Params,'Preconditioning Solver',Found )  
  Monolithic = ListGetLogical( Params,'Monolithic Solver',Found )
  
  ! We solve the equation component-wise. Hence the primary variable is a temporal one. 
  IF( Monolithic ) THEN
    ! We use different naming convention if this is used as preconditioner. 
    IF( PrecUse ) THEN
      CALL ListAddNewString( Params,'Variable',&
          "Prec Elfield[Prec Elfield re:3 Prec Elfield im:3]" )
    ELSE
      CALL ListAddNewString( Params,'Variable',&
          "Elfield[Elfield re:3 Elfield im:3]" )
    END IF
  ELSE    
    CALL ListAddNewLogical( Params,'Variable Output',.FALSE.)
    CALL ListAddNewString( Params,'Variable','Etmp[Etmp re:1 Etmp im:1]')

    ! We use different naming convention if this is used as preconditioner. 
    IF( PrecUse ) THEN
      CALL ListAddString( Params,&
          NextFreeKeyword('Exported Variable', Params), &
          "Prec Elfield[Prec Elfield re:3 Prec Elfield im:3]" )
    ELSE
      CALL ListAddString( Params,&
          NextFreeKeyword('Exported Variable', Params), &
          "Elfield[Elfield re:3 Elfield im:3]" )
    END IF
  END IF

  CALL ListAddNewLogical( Params, "Linear System Complex", .TRUE.)  

    
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

  CALL ListAddInteger( Params,'Time Derivative Order', 0 )  
  
END SUBROUTINE VectorHelmholtzNodal_Init


!-----------------------------------------------------------------------------
!> A solver for the vector Helmholtz equation based on nodal basis  
!> functions
!------------------------------------------------------------------------------
SUBROUTINE VectorHelmholtzNodal( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm(3)
  INTEGER :: i, n, nb, nd, t, active, dim, RelOrder
  INTEGER :: iter, maxiter, compi, compn, compj, dofs
  LOGICAL :: Found, VecAsm, InitHandles, &
      PrecUse, PiolaVersion, SecondOrder, UseEdgeResidual, &
      Monolithic, Segregated 
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: EF, EiVar, EdgeLoadVar, EdgeSolVar
  REAL(KIND=dp) :: mu0inv, eps0, rob0, omega
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
  COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)   
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

  EiVar => Solver % Variable
  dofs = EiVar % Dofs / 2 
  IF( dofs == 1 ) THEN
    Monolithic = .FALSE.
    CALL Info(Caller,'Treating the equation in segregated manner!')
  ELSE IF( dofs == dim ) THEN
    Monolithic = .TRUE.
    CALL Info(Caller,'Treating the equation in monolithic manner!')
  ELSE
    CALL Fatal(Caller,'Invalid number of dofs in solver variable: '//I2S(dofs))
  END IF
  Segregated = .NOT. Monolithic
  
  SecondOrder = ListGetLogicalAnySolver( Model, 'Quadratic Approximation')
  IF( SecondOrder ) THEN
    PiolaVersion = .TRUE.
  ELSE
    PiolaVersion = ListGetLogicalAnySolver(Model, 'Use Piola Transform' ) .OR. &
        ListGetLogicalAnySOlver(Model,"Second Kind Basis")    
  END IF

  PrecUse = ListGetLogical( Params,'Preconditioning Solver',Found ) 
  EdgeLoadVar => NULL()
  UseEdgeResidual = .FALSE.
  
  IF( PrecUse ) THEN
    IF(Monolithic) THEN
      CALL Fatal(Caller,'You cannot combine monolithic and preconditioning use!')
    END IF
    EF => VariableGet( Mesh % Variables,'Prec ElField')        
    sname = ListGetString( Params,'Edge Residual Name',Found)
    IF(Found) THEN
      EdgeLoadVar => VariableGet( Mesh % Variables, sname )
      IF(.NOT. ASSOCIATED( EdgeLoadVar ) ) CALL Fatal(Caller,'Could not find field: '//TRIM(sname))
      IF( PiolaVersion ) CALL Fatal(Caller,'Cannot yet handle Piola version edge elements!')
      UseEdgeResidual = .TRUE.
    END IF
  ELSE
    EF => VariableGet( Mesh % Variables,'ElField')
  END IF
    
  IF(.NOT. ASSOCIATED(EF) ) THEN
    CALL Fatal(Caller,'Variable for Electric field not found!')
  END IF  
  
  IF( ListGetLogical( Params,'Follow P Curvature', Found )  ) THEN
    CALL FollowCurvedBoundary( Model, Mesh, .TRUE. ) 
  END IF
  
  CALL DefaultStart()

  maxiter = ListGetInteger( Params,'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  RelOrder = GetInteger( Params,'Relative Integration Order',Found ) 
  CALL InitStuff()

  IF( Monolithic ) THEN
    compn = 1
  ELSE
    compn = dim
  END IF
    
  DO compi=1,compn
    IF( .NOT. Monolithic ) THEN
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
    IF( UseEdgeResidual ) THEN
      CALL EdgeToNodeProject()
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
        nb = GetElementNOFBDOFs(Element)
        CALL LocalMatrixBC(  Element, n, nd+nb, nb, InitHandles )
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    IF( Segregated) THEN
      EiVar % Values(1::2) = EF % Values(compi::2*dim) 
      EiVar % Values(2::2) = EF % Values(compi+dim::2*dim) 
    END IF
      
    Norm(compi) = DefaultSolve()

    IF( Segregated ) THEN
      IF( InfoActive(25) ) THEN
        CALL VectorValuesRange(EiVar % Values,SIZE(EiVar % Values),'E'//I2S(compi))       
        PRINT *,'Component Norm:',Norm(compi)
      END IF
      EF % Values(compi::2*dim) = EiVar % Values(1::2)
      EF % Values(compi+dim::2*dim) = EiVar % Values(2::2)
    END IF
  END DO ! compi

  EdgeSolVar => NULL()
  sname = ListGetString( Params,'Edge Solution Name',Found)
  IF(Found) THEN
    EdgeSolVar => VariableGet( Mesh % Variables, sname )
    IF(.NOT. ASSOCIATED(EdgeSolVar)) THEN
      CALL Fatal(Caller,'Edge solution not found: '//TRIM(sname))
    END IF
    CALL Info(Caller,'Projecting nodal solution to edge field')
    CALL NodeToEdgeProject()
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


  
  ! Project nodal field to be an edge field.
  !-------------------------------------------------------
  SUBROUTINE NodeToEdgeProject()
    INTEGER :: k,k1,k2,i1,i2,j,j1,j2,n0,dofi,voffset
    REAL(KIND=dp) :: NodalSol(3), EdgeVector(3)
    TYPE(Element_t), POINTER :: Edge

    IF(EdgeSolVar % Dofs /= 2) THEN
      CALL Fatal(Caller,'We should have 2 dofs for edge field!')
    END IF
    IF(.NOT. ASSOCIATED( Mesh % Edges ) ) THEN
      CALL Fatal(Caller,'Mesh edges not associated!')
    END IF
   
    n0 = Mesh % NumberOfNodes

    DO j=1, Mesh % NumberOfEdges      
      Edge => Mesh % Edges(j)

      i1 = Edge % NodeIndexes(1)
      i2 = Edge % NodeIndexes(2)

      ! Vector in the direction of the edge
      EdgeVector(1) = Mesh % Nodes % x(i2) - Mesh % Nodes % x(i1)
      EdgeVector(2) = Mesh % Nodes % y(i2) - Mesh % Nodes % y(i1)
      EdgeVector(3) = Mesh % Nodes % z(i2) - Mesh % Nodes % z(i1)
       
      IF( ParEnv % PEs > 1 ) THEN                            
        j1 = Mesh % ParallelInfo % GlobalDOFs(i1)             
        j2 = Mesh % ParallelInfo % GlobalDOFs(i2)             
      ELSE
        j1 = i1
        j2 = i2
      END IF
        
      ! There is an ordering convention that determines the direction of the edge vector
      IF( j1 < j2) EdgeVector = -EdgeVector

      k1 = EF % Perm(i1)
      k2 = EF % Perm(i2)
      k = EdgeSolVar % Perm(n0+j)

      ! Integration using the mid-point rule:
      DO dofi=1,EdgeSolVar % dofs
        voffset = 3*(dofi-1)

        ! Value at center of edge
        NodalSol(1) = ( EF % Values(6*k1-5+voffset) + EF % Values(6*k2-5+voffset) ) / 2
        NodalSol(2) = ( EF % Values(6*k1-4+voffset) + EF % Values(6*k2-4+voffset) ) / 2
        NodalSol(3) = ( EF % Values(6*k1-3+voffset) + EF % Values(6*k2-3+voffset) ) / 2

        ! We could add this one direction at the time but here we have the full vector
        ! that we project to Re and Im parts.  
        EdgeSolVar % Values(2*(k-1)+dofi) = SUM(NodalSol*EdgeVector) 
      END DO
    END DO
    
  END SUBROUTINE NodeToEdgeProject
    
!------------------------------------------------------------------------------
! Project edge residual to nodal residual.
!------------------------------------------------------------------------------
  SUBROUTINE EdgeToNodeProject()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd, nb
    TYPE(Element_t), POINTER :: Element, Edge
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:),dBasisdx(:,:), WBasis(:,:), RotWBasis(:,:)
    COMPLEX(KIND=dp), ALLOCATABLE, SAVE :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: weight, DetJ, Coord(3), s1, s2
    LOGICAL :: Stat,Found,NormLoop
    INTEGER :: i,j,k,t,p,q,m,ipi,i1,i2,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    INTEGER :: n0, nedge
    REAL(KIND=dp), ALLOCATABLE :: EdgeWeight(:), NodeWeight(:)
    REAL(KIND=dp), SAVE :: EdgeVector(3)
    COMPLEX(KIND=dp) :: c
 !------------------------------------------------------------------------------

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      m = MAX(Mesh % MaxElementDofs,20)
      ALLOCATE(Basis(m), dBasisdx(m,3), RotWBasis(m,3), Wbasis(m,3), &
          STIFF(m,m), FORCE(m), STAT=allocstat)      
      IF (allocstat /= 0) CALL Fatal(Caller,'Local storage allocation failed')
    END IF
    STIFF = 0._dp
    
    Active = GetNOFActive(Solver)
    NormLoop = .TRUE.

    ALLOCATE(EdgeWeight(Mesh % NumberOfEdges))
    EdgeWeight = 0.0_dp
    n0 = Mesh % NumberOfNodes
    
1   CONTINUE
    
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      nb = GetElementNOFBDOFs(Element)

      IF( RelOrder /= 0 ) THEN
        IP = GaussPoints( Element, RelOrder = RelOrder)
      ELSE
        IP = GaussPoints( Element )
      END IF

      CALL GetElementNodes( Nodes, UElement=Element )

      ! Initialize
      IF (.NOT. NormLoop) THEN
        FORCE = 0._dp
      END IF

      DO ipi=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------

        stat = ElementInfo( Element, Nodes, IP % U(ipi), IP % V(ipi), &
            IP % W(ipi), detJ, Basis, dBasisdx, EdgeBasis = Wbasis, &
            RotBasis = RotWBasis ) !, USolver = pSolver )

        nedge = Element % TYPE % NumberOfEdges        
        DO i=1,nedge
          j = Element % EdgeIndexes(i)                    
          
          IF( NormLoop ) THEN
            s2 = SQRT(SUM(WBasis(i,:)**2))
            weight = IP % s(ipi) * s2
            EdgeWeight(j) = EdgeWeight(j) + Weight 
          ELSE
            Edge => Mesh % Edges(j)
            k = EdgeLoadVar % Perm(n0 + j)
          
            i1 = Edge % NodeIndexes(1)
            i2 = Edge % NodeIndexes(2)
          
            ! Vector in the direction of the edge
            EdgeVector(1) = Mesh % Nodes % x(i2) - Mesh % Nodes % x(i1)
            EdgeVector(2) = Mesh % Nodes % y(i2) - Mesh % Nodes % y(i1)
            EdgeVector(3) = Mesh % Nodes % z(i2) - Mesh % Nodes % z(i1)
          
            ! Integration length of the edge
            s1 = SQRT(SUM(EdgeVector**2))
            
            weight = IP % s(ipi) / EdgeWeight(j)            
            c = Wbasis(i,compi) * &
                s1 * CMPLX(EdgeLoadVar % Values(2*k-1), EdgeLoadVar % Values(2*k))           
            FORCE(1:nd) = FORCE(1:nd) + Basis(1:nd) * weight * c  
          END IF
        END DO
      END DO

      IF(.NOT. NormLoop ) THEN
        CALL CondensateP( nd-nb, nb, STIFF, FORCE )    
        CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element)
      END IF
    END DO
    
    IF( NormLoop ) THEN
      NormLoop = .FALSE.
      GOTO 1 
    END IF

       
!------------------------------------------------------------------------------
  END SUBROUTINE EdgeToNodeProject
!------------------------------------------------------------------------------

    
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
    COMPLEX(KIND=dp) :: muinvAtIp, EpsAtIp, CurrAtIp(3)
    LOGICAL :: Stat,Found
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
      ALLOCATE(Basis(m), dBasisdx(m,3), STIFF(m,m,3), FORCE(m,3), STAT=allocstat)      
      IF (allocstat /= 0) CALL Fatal(Caller,'Local storage allocation failed')
    END IF

    CALL GetElementNodes( Nodes, UElement=Element )

    ! Initialize
    STIFF = 0._dp
    FORCE = 0._dp
    
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      Weight = IP % s(t) * DetJ

      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------      
      muinvAtIp = ListGetElementComplex( MuCoeff_h, Basis, Element, Found, GaussPoint = t )      
      IF( Found ) THEN
        muinvAtIp = muinvAtIp * mu0inv
      ELSE
        muinvAtIp = mu0inv
      END IF
      STIFF(1:nd,1:nd,1) = STIFF(1:nd,1:nd,1) + Weight * &
          MuinvAtIp * MATMUL( dBasisdx(1:nd,:), TRANSPOSE( dBasisdx(1:nd,:) ) )

      CondAtIp = ListGetElementReal( CondCoeff_h, Basis, Element, Found, GaussPoint = t )      
      IF( Found ) THEN
        DO p=1,nd
          STIFF(p,1:nd,1) = STIFF(p,1:nd,1) - im * Weight * Omega * CondAtIP * Basis(p) * Basis(1:nd)
        END DO
      END IF
              
      EpsAtIp = ListGetElementComplex( EpsCoeff_h, Basis, Element, Found, GaussPoint = t )      
      IF( Found ) THEN
        EpsAtIp = EpsAtIp * eps0
      ELSE
        epsAtIp = eps0
      END IF        

      ! This is the same for each component with isotropic materials!
      DO p=1,nd
        STIFF(p,1:nd,1) = STIFF(p,1:nd,1) - Weight * Omega**2 * epsAtIP * Basis(p) * Basis(1:nd)
      END DO

      IF(.NOT. UseEdgeResidual ) THEN
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
    END DO
    
    IF( Monolithic ) THEN
      DO i=2,dofs
        STIFF(1:nd,1:nd,i) = STIFF(1:nd,1:nd,1)
      END DO
      ! Note that we use diagonal form for this!
      CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element)
    ELSE
      CALL DefaultUpdateEquations(STIFF(:,:,1),FORCE(:,1),UElement=Element)
    END IF

      
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------



! Assembly of the matrix entries arising from the Neumann and Robin conditions.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd, nb, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: F,C,Ext, Weight, coeff
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), DetJ,Coord(3),Normal(3)
    COMPLEX(KIND=dp) :: STIFF(nd,nd,3), FORCE(nd,3)
    COMPLEX(KIND=dp) :: muInvAtIp, TemGrad(3), L(3), B 
    LOGICAL :: Stat,Found,RobinBC,NT
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC       
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(Element_t), POINTER :: Parent
    TYPE(ValueHandle_t), SAVE :: ElRobin_h, MagLoad_h, Absorb_h, TemRe_h, TemIm_h, MuCoeff_h
    
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( ElRobin_h,'Boundary Condition','Electric Robin Coefficient',InitIm=.TRUE.)
      CALL ListInitElementKeyword( MagLoad_h,'Boundary Condition','Magnetic Boundary Load', InitIm=.TRUE.,InitVec3D=.TRUE.)
      CALL ListInitElementKeyword( Absorb_h,'Boundary Condition','Absorbing BC')
      CALL ListInitElementKeyword( TemRe_h,'Boundary Condition','TEM Potential')
      CALL ListInitElementKeyword( TemIm_h,'Boundary Condition','TEM Potential Im')
      CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)      
      InitHandles = .FALSE.
    END IF
    
    CALL GetElementNodes( Nodes, UElement=Element )
    STIFF = 0._dp
    FORCE = 0._dp

    Normal = NormalVector( Element, Nodes )
    NT = ListGetLogical( BC,'Normal-Tangential '//GetVarName(Solver % Variable), Found )
    IF(NT .AND. .NOT. Monolithic) THEN
      CALL Fatal(Caller,'Normal-tangential conditions require monolithic solver!')
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

      IF( .NOT. UseEdgeResidual ) THEN
        L = ListGetElementComplex3D( MagLoad_h, Basis, Element, Found, GaussPoint = t )
        TemGrad = CMPLX( ListGetElementRealGrad( TemRe_h,dBasisdx,Element,Found), &
            ListGetElementRealGrad( TemIm_h,dBasisdx,Element,Found) )
        L = L + TemGrad
        DO i=1,dim
          FORCE(1:nd,i) = FORCE(1:nd,i) - muinvAtIp * L(i) * Basis(1:nd) * Weight
        END DO
      END IF
        
      IF( ListGetElementLogical( Absorb_h, Element, Found ) ) THEN
        B = CMPLX(0.0_dp, rob0 ) 
      ELSE
        B = ListGetElementComplex( ElRobin_h, Basis, Element, Found, GaussPoint = t )
      END IF

      IF( Found ) THEN
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
                Basis(p) * Basis(1:nd) * detJ * IP % s(t)
          END DO
        END DO
      END IF
    END DO
    
    IF( Monolithic ) THEN
      ! For normal-tangential coordinate system the slip 
      CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element)
    ELSE
      CALL DefaultUpdateEquations(STIFF(:,:,compi),FORCE(:,compi),UElement=Element)
    END IF
      
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE VectorHelmholtzNodal
!------------------------------------------------------------------------------

