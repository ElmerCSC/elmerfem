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
! *  Authors: Peter Råback, Juha Ruokolainen
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 20.06.2007
! *
! *****************************************************************************/



!------------------------------------------------------------------------------
!> Subroutine for computing fluxes and gradients of scalar fields. 
!> For example, one may compute the the heat flux as the negative gradient of temperature
!> field multiplied by the heat conductivity.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE FluxSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------

  USE CoordinateSystems
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model    !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt       !< Timestep size for time dependent simulations
  LOGICAL :: Transient      !< Steady state or transient simulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, CondName, PotName
  INTEGER :: i,j,k,dim,DOFs,firstmag
  LOGICAL :: ConstantBulkMatrix, ConstantBulkMatrixInUse, CSymmetry, GotIt  
  LOGICAL :: CalculateFluxMag, CalculateFlux, CalculateGradMag, CalculateGrad, &
      EnforcePositiveMagnitude, UsePot
  REAL(KIND=dp) :: Unorm, Totnorm, val
  REAL(KIND=dp), ALLOCATABLE, TARGET :: ForceVector(:,:)
  REAL(KIND=dp), POINTER CONTIG :: SaveRHS(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at0,at1,at2
#else
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
#endif
  TYPE(Variable_t), POINTER :: FluxSol
  TYPE FieldTable_t
    REAL(KIND=dp), POINTER :: Values(:) 
  END TYPE FieldTable_t
  TYPE(FieldTable_t) :: Fields(8)
  
 
  CALL Info( 'FluxSolver', '-------------------------------------',Level=4 )
  CALL Info( 'FluxSolver', 'Computing the flux and/or gradient',Level=4 )
  CALL Info( 'FluxSolver', '-------------------------------------',Level=4 )

  dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!  Check what needs to be computed
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN
  
  SolverParams => GetSolverParams()

  CalculateGradMag = GetLogical(SolverParams,'Calculate Grad Magnitude',GotIt)
  CalculateGrad = GetLogical(SolverParams,'Calculate Grad',GotIt) 
  
  CalculateFluxMag = GetLogical(SolverParams,'Calculate Flux Magnitude',GotIt)
  CalculateFlux = GetLogical(SolverParams,'Calculate Flux',GotIt) 
  
  Dofs = 0
  IF( CalculateFlux ) Dofs = Dofs + Dim
  IF( CalculateFluxMag ) Dofs = Dofs + 1
  IF( CalculateGrad ) Dofs = Dofs + Dim
  IF( CalculateGradMag ) Dofs = Dofs + 1

  IF( Dofs == 0 ) THEN
    CALL Warn('FluxSolver','No field computation requested, exiting...')     
    RETURN
  END IF

!-------------------------------------------------------------------------------
! If only one component is used use the scalar equation, otherwise use an
! auxiliary variable to store all the dimensions
!-------------------------------------------------------------------------------

  VarName = GetString(SolverParams,'Flux Variable',GotIt )
  UsePot = .FALSE.
  IF(.NOT. GotIt) VarName = GetString(SolverParams,'Target Variable',GotIt )
  IF(.NOT. gotIt) THEN
    PotName = GetString(SolverParams,'Target Expression',UsePot )
    IF( UsePot ) THEN
      VarName = PotName
    ELSE
      CALL Warn('FluxSolver','> Target Variable < not given, using Temperature')
      VarName = TRIM('Temperature')
    END IF
  END IF

  IF( CalculateFlux .OR. CalculateFluxMag ) THEN
    CondName = ListGetString(SolverParams,'Flux Coefficient',GotIt )
    IF(.NOT. gotIt) THEN
      CALL Warn('FluxSolver','> Flux Coefficient < not given, using Heat Conductivity')
      CondName = TRIM('Heat Conductivity')
    END IF
  END IF

  i = 0
  IF( CalculateFlux ) THEN
    FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 1' )
    Fields(1) % Values => FluxSol % Values

    FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 2' )
    Fields(2) % Values => FluxSol % Values
    
    IF( dim == 3 ) THEN
      FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 3' )
      Fields(3) % Values => FluxSol % Values
    END IF
    i = i + dim
  END IF

  IF( CalculateGrad ) THEN
    FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Grad 1' )
    Fields(i+1) % Values => FluxSol % Values

    FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Grad 2' )
    Fields(i+2) % Values => FluxSol % Values

    IF( dim == 3 ) THEN
      FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Grad 3' )
      Fields(i+3) % Values => FluxSol % Values
    END IF
    i = i + dim
  END IF

  firstmag = i + 1
  IF( CalculateFluxMag ) THEN
    FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux_mag' )
    i = i + 1
    Fields(i) % Values => FluxSol % Values
  END IF

  IF( CalculateGradMag ) THEN
    FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Grad_mag' )
    i = i + 1
    Fields(i) % Values => FluxSol % Values
  END IF

  CSymmetry = CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric
  
  at0 = RealTime()
  
  ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', GotIt )
  ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
      ASSOCIATED(Solver % Matrix % BulkValues)
  
  IF ( ConstantBulkMatrixInUse ) THEN
    Solver % Matrix % Values = Solver % Matrix % BulkValues        
    Solver % Matrix % rhs = 0.0_dp
  ELSE
    CALL DefaultInitialize()
  END IF
  
  ALLOCATE(ForceVector(SIZE(Solver % Matrix % RHS),DOFs))  
  ForceVector = 0.0_dp
  SaveRHS => Solver % Matrix % RHS
  
  CALL BulkAssembly()
  IF( ConstantBulkMatrix ) THEN
    IF(.NOT. ConstantBulkMatrixInUse ) THEN
      CALL DefaultFinishBulkAssembly( BulkUpdate = .TRUE.)
    END IF
  END IF

  CALL DefaultFinishAssembly()
  
  at1 = RealTime()
  WRITE(Message,* ) 'Assembly Time: ',at1-at0
  CALL Info( 'FluxSolver', Message, Level=5 )
!        
!------------------------------------------------------------------------------     

  EnforcePositiveMagnitude = GetLogical( SolverParams,'Enforce Positive Magnitude',GotIt )
      
  TotNorm = 0.0_dp
  DO i=1,Dofs
    Solver % Matrix % RHS => ForceVector(:,i)
    UNorm = DefaultSolve()

    TotNorm = TotNorm + Unorm ** 2
    Fields(i) % Values = Solver % Variable % Values

    IF( EnforcePositiveMagnitude .AND. i >= firstmag ) THEN
      DO j=1,SIZE( Fields(i) % Values ) 
        Fields(i) % Values(j) = MAX( Fields(i) % Values(j), 0.0_dp )
      END DO
    END IF
  END DO
  DEALLOCATE( ForceVector )  
  Solver % Matrix % RHS => SaveRHS
  TotNorm = SQRT(TotNorm)
  Solver % Variable % Norm = Totnorm

  FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux_abs' )
  IF( ASSOCIATED( FluxSol ) ) THEN
    DO j = 1,SIZE( FluxSol % Values )
      val = Fields(1) % Values(j)**2 + Fields(2) % Values(j)**2
      IF( dim == 3 ) val = val + Fields(3) % Values(j)**2 
      FluxSol % Values(j) = SQRT( val ) 
    END DO
  END IF

  FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Grad_abs' )
  IF( ASSOCIATED( FluxSol ) ) THEN
    k = 0
    IF( CalculateFlux ) k = dim
    DO j = 1,SIZE( FluxSol % Values )
      val = Fields(k+1) % Values(j)**2 + Fields(k+2) % Values(j)**2
      IF( dim == 3 ) val = val + Fields(k+3) % Values(j)**2 
      FluxSol % Values(j) = SQRT( val ) 
    END DO
  END IF


!------------------------------------------------------------------------------     

  at2 = RealTime()
  WRITE(Message,* ) 'Solution Time: ',at2-at1
  CALL Info( 'FluxSolver', Message, Level=5 )
  
  WRITE( Message, * ) 'Result Norm: ',TotNorm
  CALL Info( 'FluxSolver', Message, Level=4 )

  CALL Info( 'FluxSolver', 'All done',Level=6 )
  CALL Info( 'FluxSolver', '-------------------------------------',Level=6 )


  
CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
       
    INTEGER :: elem,t,i,j,k,p,q,n,nd, Rank
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,C(3,3),coeff,detJ,FluxAtIp(3),GradAtIp(3)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: LocalPotential(:)
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp), POINTER :: Conductivity(:,:,:)=>NULL()
    
    SAVE Conductivity, Nodes
    
    n = 2*MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    ALLOCATE( STIFF(n,n), FORCE(dofs,n) )
    ALLOCATE( LocalPotential(n), Basis(n), dBasisdx(n,3) )

    DO elem = 1,Solver % NumberOFActiveElements
         
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      Material => GetMaterial()
      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()
      
      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      IF( CalculateFlux .OR. CalculateFluxMag ) THEN
        CALL GetRealArray( Material, Conductivity, CondName, Found )
        Rank = 0
        IF ( Found ) Rank = GetTensorRank(Conductivity)
        
        C = 0.0_dp
        FluxAtIp = 0.0_dp
        DO i=1,dim
          C(i,i) = 1.0_dp
        END DO
      END IF

      IF( UsePot ) THEN
        LocalPotential(1:n) = GetReal( Material, PotName ) 
      ELSE
        CALL GetScalarLocalSolution( LocalPotential, VarName )
      END IF

      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
            IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )
        
        Weight = IntegStuff % s(t) * detJ
        IF ( CSymmetry ) Weight = Weight * SUM( Basis(1:n) * Nodes % x(1:n) )
        
        IF ( .NOT. ConstantBulkMatrixInUse ) THEN
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
            END DO
          END DO
        END IF
        
        GradAtIp(1:dim) = MATMUL( LocalPotential(1:nd), dBasisdx(1:nd,1:dim) )
        
        IF( CalculateFlux .OR. CalculateFluxMag ) THEN          
          SELECT CASE(Rank)
          CASE(0)
          CASE(1)
            DO i=1,dim
              C(i,i) = SUM( Basis(1:n) * Conductivity(1,1,1:n) )
            END DO
          CASE(2)
            DO i=1,dim
              C(i,i) = SUM( Basis(1:n) * Conductivity(i,1,1:n) )
            END DO
          CASE DEFAULT
            DO i=1,dim
              DO j=1,dim
                C(i,j) = SUM( Basis(1:n) * Conductivity(i,j,1:n) )
              END DO
            END DO
          END SELECT
          
          DO i=1,dim
            FluxAtIp(i) = SUM( C(i,1:dim) * GradAtIp(1:dim) )
          END DO
        END IF


        k = 0
        IF( CalculateFlux ) THEN
          DO i=1,dim
            Coeff = Weight * FluxAtIp(i)
            FORCE(i,1:nd) = FORCE(i,1:nd) - Coeff * Basis(1:nd)
          END DO
          k = k + dim
        END IF
      
        IF( CalculateGrad ) THEN
          DO i=1,dim
            Coeff = Weight * GradAtIp(i)
            FORCE(k + i,1:nd) = FORCE(k + i,1:nd) + Coeff * Basis(1:nd)
          END DO
          k = k + dim
        END IF

        IF( CalculateFluxMag ) THEN
          Coeff = Weight * SQRT( SUM(FluxAtIp ** 2) )
          k = k + 1
          FORCE(k,1:nd) = FORCE(k,1:nd) + Coeff * Basis(1:nd)
        END IF
        
        IF( CalculateGradMag ) THEN
          Coeff = Weight * SQRT( SUM( GradAtIp ** 2) )
          k = k + 1
          FORCE(k,1:nd) = FORCE(k,1:nd) + Coeff*Basis(1:nd)
        END IF

      END DO

!------------------------------------------------------------------------------
!      Update global matrices from local matrices 
!------------------------------------------------------------------------------
      IF ( .NOT. ConstantBulkMatrixInUse ) THEN
        Solver % Matrix % Rhs => SaveRHS
        CALL DefaultUpdateEquations( STIFF, FORCE(1,1:nd) )
      END IF

      DO i=1,Dofs
        Solver % Matrix % RHS => ForceVector(:,i)
        CALL DefaultUpdateForce( FORCE(i,1:nd) )
      END DO

    END DO

    ! Assembly of the face terms:
    !----------------------------

    IF (GetLogical(GetSolverParams(),'Discontinuous Galerkin',Found)) THEN
      IF (GetLogical(GetSolverParams(),'Average Within Materials',Found)) THEN
        FORCE = 0.0d0
        CALL AddLocalFaceTerms( STIFF, FORCE(1,:) )
      END IF
    END IF

    DEALLOCATE( LocalPotential, STIFF, FORCE, Basis, dBasisdx )

!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
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
      INTEGER :: i, j, p, q, t, nFace, nParent
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

      

!------------------------------------------------------------------------------
  FUNCTION GetTensorRank( Tensor ) RESULT ( Rank )
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: Tensor(:,:,:)
    INTEGER :: Rank
    
    IF ( SIZE(Tensor,1) == 1 ) THEN
      Rank = 1
    ELSE IF ( SIZE(Tensor,2) == 1 ) THEN
      Rank = 2
    ELSE
      Rank = 3
    END IF
!-----------------------------------------------------------------------------
  END FUNCTION GetTensorRank
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE FluxSolver
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Initialization for the primary solver, FluxSolver. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE FluxSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE

    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: DT
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: SolverParams
    INTEGER :: dim
    CHARACTER(LEN=MAX_NAME_LEN) :: EqName, VarName, FluxName, GradName
    LOGICAL :: GotIt, CalculateFluxMag, CalculateFlux, CalculateGrad, &
        CalculateGradMag
!------------------------------------------------------------------------------
    SolverParams => GetSolverParams()
    dim = CoordinateSystemDimension()


    IF( dim < 2 .OR. dim > 3 ) THEN
      CALL Fatal('FluxSolver_init','Flux computation makes sense only in 2D and 3D')
    END IF

    CalculateGradMag = GetLogical(SolverParams,'Calculate Grad Magnitude',GotIt)
    CalculateGrad = GetLogical(SolverParams,'Calculate Grad',GotIt) 

    CalculateFluxMag = GetLogical(SolverParams,'Calculate Flux Magnitude',GotIt)
    CalculateFlux = GetLogical(SolverParams,'Calculate Flux',GotIt) 


    IF(.NOT. (CalculateFlux .OR. CalculateFluxMag .OR. CalculateGrad .OR. CalculateGradMag) ) THEN 
      IF(.NOT. GotIt ) THEN
        CALL Warn('FluxSolver_init','No field computation requested, setting > Calculate Flux = True <')
      ELSE
        CALL Warn('FluxSolver_init','No field computation requested!')
        RETURN
      END IF
      CalculateFlux = .TRUE.
      CALL ListAddLogical( SolverParams,'Calculate Flux',CalculateFlux)
    END IF


    VarName = GetString(SolverParams,'Flux Variable',GotIt )
    IF(.NOT. GotIt) VarName = GetString(SolverParams,'Target Variable',GotIt )
    IF(.NOT. gotIt) THEN
      IF( .NOT. GotIt ) THEN
        VarName = GetString( SolverParams,'Target Expression',GotIt ) 
      END IF
      IF( .NOT. GotIt ) THEN
        CALL Warn('FluxSolver','> Target Variable < not given, using Temperature')
        VarName = 'Temperature'
      END IF
    END IF

    IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      EqName = ListGetString( SolverParams,'Equation')
      CALL ListAddString( SolverParams, 'Variable','-nooutput '//TRIM(EqName)//'_temp' )
    END IF
    

    IF( CalculateFlux ) THEN
      FluxName = TRIM(VarName)//' Flux'
      CALL Info('FluxSolver_init','Saving flux to: '//TRIM(FluxName), Level=5) 
      IF(dim == 2) THEN
        CALL ListAddString( SolverParams,&
            NextFreeKeyword('Exported Variable',SolverParams),&
	TRIM(FluxName)//'['//TRIM(FluxName)//':2]')
      ELSE IF(dim == 3) THEN
        CALL ListAddString( SolverParams,&
            NextFreeKeyword('Exported Variable',SolverParams),&
	TRIM(FluxName)//'['//TRIM(FluxName)//':3]')
      END IF
      IF( GetLogical( SolverParams,'Calculate Flux Abs',GotIt) ) THEN
        FluxName = TRIM(VarName)//' Flux_abs'
        CALL Info('FluxSolver_init','Saving flux abs to: '//FluxName) 
        CALL ListAddString( SolverParams,&
            NextFreeKeyword('Exported Variable',SolverParams),TRIM(FluxName))
      END IF
    END IF
    IF( CalculateFluxMag ) THEN
      FluxName = TRIM(VarName)//' Flux_mag'
      CALL Info('FluxSolver_init','Saving flux magnitude to: '//FluxName) 
      CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable',SolverParams),TRIM(FluxName))
    END IF

    IF( CalculateGrad ) THEN
      GradName = TRIM(VarName)//' Grad'
      CALL Info('FluxSolver_init','Saving gradient to: '//GradName) 	
      IF(dim == 2) THEN
        CALL ListAddString( SolverParams,&
            NextFreeKeyword('Exported Variable',SolverParams),&
	TRIM(GradName)//'['//TRIM(GradName)//':2]')
      ELSE IF(dim == 3) THEN
        CALL ListAddString( SolverParams,&
            NextFreeKeyword('Exported Variable',SolverParams),&
	TRIM(GradName)//'['//TRIM(GradName)//':3]')
      END IF
      IF( GetLogical( SolverParams,'Calculate Grad Abs',GotIt) ) THEN
        GradName = TRIM(VarName)//' Grad_abs'
        CALL Info('FluxSolver_init','Saving gradient abs to: '//GradName) 	
        CALL ListAddString( SolverParams,&
            NextFreeKeyword('Exported Variable',SolverParams),TRIM(GradName))
      END IF
    END IF
    IF( CalculateGradMag ) THEN
      GradName = TRIM(VarName)//' Grad_mag'
      CALL Info('FluxSolver_init','Saving gradient magnitude to: '//GradName) 		
      CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable',SolverParams),TRIM(GradName))
    END IF


    CALL ListAddInteger( SolverParams, 'Time derivative order', 0 )
      
    CALL ListAddLogical( SolverParams,'Skip Compute Nonlinear Change',.TRUE.)

    ! Add linear system defaults: cg+diagonal
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Solver')) &
        CALL ListAddString(SolverParams,'Linear System Solver','Iterative')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Iterative Method')) &
        CALL ListAddString(SolverParams,'Linear System Iterative Method','cg')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Preconditioning')) &
        CALL ListAddString(SolverParams,'Linear System Preconditioning','diagonal')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Max Iterations')) &
        CALL ListAddInteger(SolverParams,'Linear System Max Iterations',500)
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Residual Output')) &
        CALL ListAddInteger(SolverParams,'Linear System Residual Output',10)
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Convergence Tolerance')) &
        CALL ListAddConstReal(SolverParams,'Linear System Convergence Tolerance',1.0e-10_dp)
    
!------------------------------------------------------------------------------
  END SUBROUTINE FluxSolver_Init
!------------------------------------------------------------------------------
