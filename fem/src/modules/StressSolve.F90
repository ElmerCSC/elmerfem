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
! *  Module containing a solver for linear stress equations.
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Mikko Lyly, Peter R�back
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Initialization for the primary solver: StressSolver. 
!------------------------------------------------------------------------------
SUBROUTINE StressSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE

    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: DT
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    INTEGER :: dim,i
    TYPE(ValueList_t), POINTER :: SolverParams
    LOGICAL :: Found, CalculateStrains, CalcPrincipalAngle, CalcPrincipalAll, &
        CalcStressAll
!------------------------------------------------------------------------------

    SolverParams => GetSolverParams()
    IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      dim = CoordinateSystemDimension()
      CALL ListAddInteger( SolverParams, 'Variable DOFs', dim )
      CALL ListAddString( SolverParams, 'Variable', 'Displacement' )
    END IF
    CALL ListAddInteger( SolverParams, 'Time derivative order', 2 )

    IF( .NOT. ListCheckPresent( SolverParams,'Displace Mesh At Init') ) THEN
      CALL ListAddLogical( SolverParams,'Displace Mesh At Init',.TRUE.)
    END IF
    
    CalculateStrains = GetLogical(SolverParams, 'Calculate Strains', Found)
    CalcPrincipalAngle = GetLogical(SolverParams, 'Calculate PAngle', Found)
    CalcPrincipalAll = GetLogical(SolverParams, 'Calculate Principal', Found)
    CalcStressAll = GetLogical( SolverParams, 'Calculate Stresses',Found )
    IF(CalcPrincipalAngle) CalcPrincipalAll = .TRUE. ! can't calculate angle without principal
    IF(CalcPrincipalAll)   CalcStressAll = .TRUE. ! can't calculate principal without components
    IF(CalculateStrains)   CalcStressAll = .TRUE. ! can't calculate principal without components
    
    ! If stress computation is requested somewhere then enforce it 
    IF( .NOT. ( CalcStressAll .OR. CalculateStrains) ) THEN
      CalcStressAll = ListGetLogicalAnyEquation( Model,'Calculate Stresses')
      IF( CalcStressAll ) CALL ListAddLogical( SolverParams,'Calculate Stresses',.TRUE.)
    END IF

    IF ( CalcStressAll ) THEN
      CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable ',SolverParams), &
          'Stress[Stress_xx:1 Stress_yy:1 Stress_zz:1 Stress_xy:1 Stress_yz:1 Stress_xz:1]' )
      CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable ',SolverParams), &
          'vonMises' )
      
      IF(CalcPrincipalAll) THEN
        CALL ListAddString( SolverParams,&
            NextFreeKeyword('Exported Variable ',SolverParams), &
            'Principal Stress[Principal Stress:3]' )
        CALL ListAddString( SolverParams,&
            NextFreeKeyword('Exported Variable ',SolverParams), &
            'Tresca' )
        
        IF(CalcPrincipalAngle) THEN
          CALL ListAddString( SolverParams,&
              NextFreeKeyword('Exported Variable ',SolverParams), &
              '-dofs 9 Principal Angle' )
        END IF ! PrincipalAngle
      END IF !CalcPrincipalAll      
    END IF ! CalcStressAll
    
    IF(CalculateStrains) THEN
      CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable ',SolverParams), &
          'Strain[Strain_xx:1 Strain_yy:1 Strain_zz:1 Strain_xy:1 Strain_yz:1 Strain_xz:1]' )
      IF(CalcPrincipalAll) THEN
        CALL ListAddString( SolverParams,&
            NextFreeKeyword('Exported Variable ',SolverParams), &
            'Principal Strain[Principal Strain:3]' )
      END IF
    END IF

    CALL ListAddLogical( SolverParams, 'stress: Linear System Save', .FALSE. )

!------------------------------------------------------------------------------
  END SUBROUTINE StressSolver_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves the elastic displacement assuming linear stress-strain relationship.
!> The solver is burdened with a plethora of different options. For example,
!> various kinds of stresses may be computed. Also some basic features for
!> model lumping and contact analysis exist.
!------------------------------------------------------------------------------
   SUBROUTINE StressSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------

    USE CoordinateSystems
    USE StressLocal
    USE StressGeneral
    USE Adaptive
    USE DefUtils

    IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL ::  Transient
     REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,l,n,ntot,t,iter,STDOFs,istat, body_id

     TYPE(ValueList_t),POINTER :: SolverParams, Equation, Material, BodyForce, BC
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: Element
     TYPE(Mesh_t),POINTER :: Mesh

     TYPE(Variable_t), POINTER :: ReferenceSol
     REAL(KIND=dp), POINTER :: DispValues(:)
     REAL(KIND=dp), ALLOCATABLE :: UpdateRef(:), Ref_rhs(:), NodalRefD(:)
     INTEGER :: RefDofs, indx
     LOGICAL, ALLOCATABLE :: UpdatePresent(:), UpdateActive(:)

     REAL(KIND=dp) :: UNorm,s, UzawaParameter

     INTEGER ::  MaxIter, MinIter, NoModes, Nsize, Dofs
     TYPE(Variable_t), POINTER :: StressSol, iVar, Var, TimeVar

     CHARACTER(LEN=MAX_NAME_LEN) :: VarName

     REAL(KIND=dp), POINTER :: Temperature(:),Work(:,:,:), &
       VonMises(:), NodalStress(:), NodalStrain(:), StressComp(:), StrainComp(:), ContactPressure(:), &
       PrincipalStress(:), PrincipalStrain(:), Tresca(:), &   ! needed for principal strain calculation
       PrincipalAngle(:), PrincipalAngleComp(:), &            ! needed for principal angle calculation
       PrincipalStressComp(:), PrincipalStrainComp(:), &
       NormalDisplacement(:), TransformMatrix(:,:), UWrk(:,:), &
       RayleighAlpha(:), RayleighBeta(:), SaveRHS(:)

     REAL(KIND=dp), POINTER :: Displacement(:)

     REAL(KIND=dp) :: UnitNorm, Prevdt=-1, PrevTime=-1

     INTEGER, POINTER :: TempPerm(:),DisplPerm(:),StressPerm(:),NodeIndexes(:)

     LOGICAL :: GotForceBC,Found,RayleighDamping, NormalSpring
     LOGICAL :: PlaneStress, CalcStress, CalcStressAll, &
        CalcPrincipalAll, CalcPrincipalAngle, CalculateStrains, Isotropic(2) = .TRUE.
     LOGICAL :: Contact = .FALSE.
     LOGICAL :: stat, stat2, stat3, RotateC, MeshDisplacementActive, &
                ConstantBulkSystem, ConstantBulkMatrix, ConstantBulkMatrixInUse, ConstantSystem, &
                UpdateSystem, GotHeatExp, Converged

     LOGICAL :: AllocationsDone = .FALSE., NormalTangential, HarmonicAnalysis
     LOGICAL :: StabilityAnalysis = .FALSE., ModelLumping, FixDisplacement
     LOGICAL :: GeometricStiffness = .FALSE., EigenAnalysis=.FALSE., OrigEigenAnalysis, &
           Refactorize = .TRUE.

     REAL(KIND=dp),ALLOCATABLE:: MASS(:,:),STIFF(:,:),&
       DAMP(:,:), LOAD(:,:),LOAD_im(:,:),FORCE(:),FORCE_im(:), &
       LocalTemperature(:),ElasticModulus(:,:,:),PoissonRatio(:), &
       HeatExpansionCoeff(:,:,:),DampCoeff(:),SpringCoeff(:,:,:),Beta(:), &
       ReferenceTemperature(:), Density(:), Damping(:), Beta_im(:), &
       NodalDisplacement(:,:), ContactLimit(:), LocalNormalDisplacement(:), &
       LocalContactPressure(:), PreStress(:,:), PreStrain(:,:), &
       StressLoad(:,:), StrainLoad(:,:), NodalMeshVelo(:,:)

     SAVE MASS,DAMP, STIFF,LOAD,LOAD_im,FORCE_im,Beta_im, &
       FORCE,ElementNodes,DampCoeff,SpringCoeff,Beta,Density, Damping, &
       LocalTemperature,AllocationsDone,ReferenceTemperature, &
       ElasticModulus, PoissonRatio,HeatExpansionCoeff, VonMises, NodalStress, &
       CalcStress, CalcStressAll, NodalDisplacement, Contact, ContactPressure, &
       NormalDisplacement, ContactLimit, LocalNormalDisplacement, &
       LocalContactPressure, PreStress, PreStrain, StressLoad, StrainLoad, Work, &
       RotateC, TransformMatrix, body_id, NodalMeshVelo, PrevTime, &
       ReferenceSol, UpdateRef, NodalRefD, Ref_rhs, UpdatePresent, UpdateActive, &
       RayleighAlpha, RayleighBeta, RayleighDamping, &
       NodalStrain, PrincipalStress, PrincipalStrain, Tresca, &
       PrincipalAngle, PrincipalAngleComp, CalcPrincipalAngle, &
       CalcPrincipalAll, CalculateStrains
!------------------------------------------------------------------------------
     INTEGER :: dim
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at,at0
#else
     REAL(KIND=dp) :: at,at0,CPUTime,RealTime
#endif
     REAL(KIND=dp) :: LumpedArea, LumpedCenter(3), LumpedMoments(3,3)

     INTERFACE
        FUNCTION StressBoundaryResidual( Model,Edge,Mesh,Quant,Perm, Gnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
          INTEGER :: Perm(:)
        END FUNCTION StressBoundaryResidual

        FUNCTION StressEdgeResidual( Model,Edge,Mesh,Quant,Perm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2)
          INTEGER :: Perm(:)
        END FUNCTION StressEdgeResidual

        FUNCTION StressInsideResidual( Model,Element,Mesh,Quant,Perm, Fnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Element
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
          INTEGER :: Perm(:)
        END FUNCTION StressInsideResidual
     END INTERFACE
!------------------------------------------------------------------------------

     CALL Info( 'StressSolve', ' ', Level=4 )
     CALL Info( 'StressSolve', '--------------------------------------------------',Level=4 )
     CALL Info( 'StressSolve', 'Solving displacements from linear elasticity model',Level=4 )     
     CALL Info( 'StressSolve', '--------------------------------------------------',Level=4 )
 
     DIM = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
     IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN

     SolverParams => GetSolverParams()
     Mesh => Solver % Mesh

     StressSol => Solver % Variable
     DisplPerm      => StressSol % Perm
     STDOFs         =  StressSol % DOFs
     Displacement   => StressSol % Values

     IF( STDOFs < Mesh % MeshDim ) THEN
       CALL Fatal('StressSolver','Number of Dofs smaller than dim: '&
           //I2S(STDOFs)//' vs. '//I2S(Mesh % MeshDim))
     END IF

     MeshDisplacementActive = ListGetLogical( SolverParams,  &
               'Displace Mesh', Found )

     IF ( .NOT. Found ) &
       MeshDisplacementActive = .NOT.EigenOrHarmonicAnalysis()

     IF ( AllocationsDone .AND. MeshDisplacementActive ) THEN
        CALL DisplaceMesh( Mesh, Displacement, -1, DisplPerm, STDOFs )
     END IF

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Mesh % Changed) THEN
       N = Mesh % MaxElementDOFs

       IF ( AllocationsDone ) THEN
         DEALLOCATE( Density,                &
                     Damping,                &
                   RayleighAlpha,          &
                     RayleighBeta,           &
                     DampCoeff,              &
                     SpringCoeff,            &
                     ReferenceTemperature,   &
                     HeatExpansionCoeff,     &
                     LocalTemperature,       &
                     ElasticModulus,         &
                     PoissonRatio,           &
                     PreStress, PreStrain,   &
                     StressLoad, StrainLoad, &
                     NodalDisplacement,      &
                     NodalMeshVelo,          &
                     FORCE_im, LOAD_im, Beta_im, &
                     FORCE, MASS, DAMP, STIFF, LOAD, Beta, &
                     ContactLimit, LocalNormalDisplacement, LocalContactPressure, &
                     TransformMatrix ,       &
                     UpdateRef, NodalRefD, Ref_rhs,        &
                     UpdatePresent, UpdateActive )
       END IF

       ALLOCATE( Density( N ),              &
                 Damping( N ),              &
                 RayleighAlpha( N ),        &
                 RayleighBeta( N ),         &
                 DampCoeff( N ),            &
                 PreStress( 6,N ),          &
                 PreStrain( 6,N ),          &
                 StressLoad( 6,N ),         &
                 StrainLoad( 6,N ),         &
                 SpringCoeff( N,3,3 ),        &
                 ReferenceTemperature( N ), &
                 HeatExpansionCoeff( 3,3,N ), &
                 LocalTemperature( N ),       &
                 ElasticModulus(6,6,N),       &
                 PoissonRatio( N ),           &
                 FORCE( STDOFs*N ),           &
                 FORCE_im( STDOFs*N ),        &
                 MASS(  STDOFs*N,STDOFs*N ),  &
                 DAMP(  STDOFs*N,STDOFs*N ),  &
                 STIFF( STDOFs*N,STDOFs*N ),  &
                 NodalDisplacement( 3, N ),   &
                 NodalMeshVelo( 3, N ),       &
                 LOAD( 4,N ), Beta( N ),      &
                 LOAD_im( 4,N ), Beta_im( N ),      &
                 ContactLimit(N), LocalNormalDisplacement(N), &
                 LocalContactPressure(N),     &
                 UpdateRef(N),                &
                 NodalRefD(STDOFs*N),       &
                 Ref_rhs(STDOFs*N),           &
                 UpdatePresent( Model % NumberOfBodies ), &
                 UpdateActive( Model % NumberOfBodies ), &
                 TransformMatrix(3,3),  STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL Fatal( 'StressSolve', 'Memory allocation error.' )
       END IF

       NULLIFY( Work )
       TransformMatrix = 0.0d0

       CalculateStrains = GetLogical(SolverParams, 'Calculate Strains', Found)
       CalcPrincipalAngle = GetLogical(SolverParams, 'Calculate PAngle', Found)
       CalcPrincipalAll = GetLogical(SolverParams, 'Calculate Principal', Found)
       CalcStressAll = GetLogical( SolverParams, 'Calculate Stresses',Found )
       IF(CalcPrincipalAngle) CalcPrincipalAll = .TRUE. ! can't calculate angle without principal
       IF(CalcPrincipalAll)   CalcStressAll = .TRUE. ! can't calculate principal without components
       IF(CalculateStrains)   CalcStressAll = .TRUE. ! can't calculate principal without components
       
       Contact = GetLogical( SolverParams, 'Contact', Found )
       IF( Contact ) THEN
         n = SIZE( Displacement ) / STDOFs
         ALLOCATE( ContactPressure( n ), NormalDisplacement( n ) )
         ContactPressure = 0.0d0
         NormalDisplacement = 0.0d0
         CALL VariableAdd(Mesh % Variables, Mesh, Solver, &
             'Contact Pressure', 1, ContactPressure, DisplPerm )
       END IF

!------------------------------------------------------------------------------
!      Check if reference displacement present
!------------------------------------------------------------------------------
       UpdatePresent = .FALSE.
       ReferenceSol => VariableGet( Mesh % Variables, 'Reference Displacement' )
       IF ( ASSOCIATED( ReferenceSol ) ) THEN
         DO i = 1, Model % NumberOfBodies
           j = GetInteger( Model % Bodies(i) % Values, 'Body Force', Found )
           IF ( .NOT. Found )  CYCLE

           UpdatePresent(i) = ListCheckPresent( Model % BodyForces(j) % Values, &
               'Update Reference Displacement' )
         END DO
       END IF
!------------------------------------------------------------------------------

       AllocationsDone = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     IF ( CalcStressAll ) THEN
       Var => VariableGet( Mesh % Variables, 'Stress',.TRUE. )
       IF ( ASSOCIATED( Var ) ) THEN
         StressPerm  => Var % Perm
         NodalStress => Var % Values
       ELSE  
         CALL Fatal('StressSolver','Variable > Stress < does not exits!')
       END IF

       Var => VariableGet( Mesh % Variables, 'VonMises',.TRUE. )
       IF ( ASSOCIATED( Var ) ) THEN
         VonMises => Var % Values
       ELSE
         CALL Fatal('StressSolver','Variable > vonMises < does not exits!')
       END IF
       
       IF(CalcPrincipalAll) THEN
         Var => VariableGet( Mesh % Variables, 'Principal Stress',.TRUE. )
         IF ( ASSOCIATED( Var ) ) THEN
           PrincipalStress => Var % Values
         ELSE                 
           CALL Fatal('StressSolver','Variable > Principal Stress < does not exits!')
         END IF
         
         Var => VariableGet( Mesh % Variables, 'Tresca',.TRUE. )
         IF ( ASSOCIATED( Var ) ) THEN
           Tresca => Var % Values
         ELSE
           CALL Fatal('StressSolver','Variable > Tresca < does not exits!')
         END IF
         
         IF(CalcPrincipalAngle) THEN
           Var => VariableGet( Mesh % Variables, 'Principal Angle' )                 
           IF ( ASSOCIATED( Var ) ) THEN
             PrincipalAngle => Var % Values
           ELSE
             CALL Fatal('StressSolver','Variable > Principal Angle < does not exits!')
           END IF
         END IF ! PrincipalAngle
       END IF !CalcPrincipalAll             
     END IF ! CalcStress or CalcStressAll
     
     IF(CalculateStrains) THEN
       Var => VariableGet( Mesh % Variables, 'Strain' )
       IF ( ASSOCIATED( Var ) ) THEN
         NodalStrain => Var % Values
       ELSE
         CALL Fatal('StressSolver','Variable > Strain < does not exits!')
       END IF
       IF(CalcPrincipalAll) THEN
         Var => VariableGet( Mesh % Variables, 'Principal Strain' )
         IF ( ASSOCIATED( Var ) ) THEN
           PrincipalStrain => Var % Values
         ELSE
           CALL Fatal('StressSolver','Variable > Principal Strain < does not exits!')
         END IF
       END IF
     END IF !Calculate strains
     
     IF ( ANY( UpdatePresent ) ) THEN
       RefDofs = ReferenceSol % DOFs
       UpdateActive = .FALSE.
       DispValues => ReferenceSol % Values
     END IF

!------------------------------------------------------------------------------
     MaxIter = GetInteger( SolverParams, &
         'Nonlinear System Max Iterations',Found )
     IF ( .NOT.Found ) MaxIter = 1

     MinIter = GetInteger( SolverParams, &
         'Nonlinear System Min Iterations',Found )

     EigenAnalysis = GetLogical( SolverParams, 'Eigen Analysis', Found )
     OrigEigenAnalysis = EigenAnalysis

     StabilityAnalysis = GetLogical( SolverParams, 'Stability Analysis', Found )
     IF( .NOT. Found ) StabilityAnalysis = .FALSE.

     IF( StabilityAnalysis .AND. (CurrentCoordinateSystem() /= Cartesian) ) &
         CALL Fatal( 'StressSolve', &
          'Only cartesian coordinate system is allowed in stability analysis.' )

     GeometricStiffness = GetLogical( SolverParams, 'Geometric Stiffness', Found )
     IF (.NOT. Found ) GeometricStiffness = .FALSE.

     IF( GeometricStiffness .AND. (CurrentCoordinateSystem() /= Cartesian) ) &
          CALL Fatal( 'StressSolve', &
          'Only cartesian coordinates are allowed with geometric stiffness.' )

     IF ( StabilityAnalysis .AND. GeometricStiffness )  &
         CALL Fatal( 'StressSolve', &
         'Stability analysis and geometric stiffening can not be activated simultaneously.' )

     IF ( StabilityAnalysis .OR. GeometricStiffness ) THEN
       MinIter = 2
       MaxIter = 2
     END IF

     HarmonicAnalysis = getLogical( SolverParams, 'Harmonic Analysis', Found )
!------------------------------------------------------------------------------
     Refactorize = GetLogical( SolverParams, 'Linear System Refactorize', Found )
     IF ( .NOT. Found ) Refactorize = .TRUE.

     IF ( Transient .AND. .NOT. Refactorize .AND. dt /= Prevdt ) THEN
       CALL ListAddLogical( SolverParams, 'Linear System Refactorize', .TRUE. )
       CALL ListAddLogical( SolverParams, &
             'Linear System Free Factorization',.FALSE.)
     END IF

     ConstantSystem = GetLogical( SolverParams, 'Constant System', Found )
     ConstantBulkSystem = GetLogical( SolverParams, 'Constant Bulk System', Found )
     ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', Found )
!------------------------------------------------------------------------------

     UpdateSystem = .FALSE.
     IF( .NOT. ListGetString( CurrentModel % Simulation,'Simulation Type') == 'steady state') THEN
       UpdateSystem = GetLogical( SolverParams, 'Update Transient System', Found )
       IF(UpdateSystem) THEN
         TimeVar => VariableGet( Mesh % Variables, 'Time' )            
         IF (ABS(TimeVar % Values(1) - PrevTime) > AEPS) THEN
           PrevTime = TimeVar % Values(1)
         ELSE
           UpdateSystem = .FALSE.
         END IF
       END IF
     END IF

     ModelLumping = GetLogical( SolverParams, 'Model Lumping', Found )
     IF ( ModelLumping ) THEN       
       IF(DIM /= 3) CALL Fatal('StressSolve','Model Lumping implemented only for 3D')
       FixDisplacement = GetLogical( SolverParams, 'Fix Displacement', Found )
       IF(.NOT. Found) FixDisplacement = .TRUE.
       IF(FixDisplacement) THEN
         CALL Info( 'StressSolve', 'Using six fixed displacement to compute the spring matrix' ) 
       ELSE
         CALL Info( 'StressSolve', 'Using six pure forces and moments to compute the spring matrix' ) 
       END IF
       MinIter = 6
       MaxIter = 6
       ConstantBulkSystem = .TRUE.
       CALL CoordinateIntegrals(LumpedArea, LumpedCenter, LumpedMoments, &
            Model % MaxElementNodes)
       CALL LumpedCartesianMass()
     END IF


     DO iter=1,MaxIter

       IF( StabilityAnalysis .OR. GeometricStiffness ) THEN
          SELECT CASE( iter )
          CASE( 1 )
            EigenAnalysis = .FALSE.
          CASE DEFAULT
            EigenAnalysis = OrigEigenAnalysis
          END SELECT
          CALL ListAddLogical( SolverParams, 'Eigen Analysis', EigenAnalysis )
       END IF

       at  = CPUTime()
       at0 = RealTime()

       IF( MaxIter > 1 ) THEN
         WRITE( Message, * ) 'Displacemet iteration: ', iter
         CALL Info( 'StressSolve', Message,Level=4 )
       END IF
       CALL Info( 'StressSolve', 'Starting assembly...',Level=5 )
!------------------------------------------------------------------------------

500    CALL DefaultInitialize()

       ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
           ASSOCIATED(Solver % Matrix % BulkValues)

       IF ( ASSOCIATED(Solver % Matrix % BulkValues) .AND. .NOT. UpdateSystem) THEN
         IF ( ConstantBulkMatrix .OR. ConstantBulkSystem .OR. ConstantSystem ) THEN
           Solver % Matrix % DampValues = 0.0d0
           Solver % Matrix % Values = Solver % Matrix % BulkValues
         END IF

         IF ( ConstantBulkSystem .OR. ConstantSystem ) THEN
           Solver % Matrix % RHS  = Solver % Matrix % BulkRHS
         ELSE IF ( ConstantBulkMatrix ) THEN
           Solver % Matrix % RHS = 0.0_dp
         END IF

         IF ( ConstantBulkMatrix ) GO TO 1000
         IF ( ConstantBulkSystem ) GO TO 2000
         IF ( ConstantSystem )     GO TO 3000
       END IF
 
!       CALL DefaultInitialize()

1000   CALL BulkAssembly()
       CALL DefaultFinishBulkAssembly()
       CALL Info( 'StressSolve', 'Bulk assembly done', Level=5 )

2000   CALL BCAssembly()
       CALL DefaultFinishBoundaryAssembly()

3000   IF ( Transient .AND.(ConstantBulkMatrix .OR. &
           ConstantBulkSystem .OR. ConstantSystem) ) CALL AddGlobalTime()
       CALL DefaultFinishAssembly()

       CALL DefaultDirichletBCS()

       IF( ModelLumping .AND. FixDisplacement) THEN
         CALL LumpedDisplacements( Model, iter, LumpedArea, LumpedCenter)
       END IF
       CALL Info( 'StressSolve', 'Set boundaries done', Level=5 )

       !------------------------------------------------------------------------------
       !     Check stepsize for nonlinear iteration
       !------------------------------------------------------------------------------
       IF( DefaultLinesearch( Converged ) ) GOTO 500

       IF( Converged ) EXIT

       ! Solve the system and check for convergence:
       !--------------------------------------------
       UNorm = DefaultSolve()

       IF ( Transient .AND. .NOT. Refactorize .AND. dt /= Prevdt ) THEN
         Prevdt = dt
         CALL ListRemove( SolverParams, 'Linear System Free Factorization' )
         CALL ListAddLogical(SolverParams,'Linear System Refactorize',.FALSE.)
       END IF

       ! Update contact pressure:
       !-------------------------
       IF ( Contact ) THEN
         CALL ComputeNormalDisplacement( Displacement, &
            NormalDisplacement, DisplPerm, STDOFs )
         
         UzawaParameter = GetConstReal( SolverParams, 'Uzawa Parameter', Found )
         IF( .NOT.Found ) THEN
            WRITE( Message, * ) 'Using default value 1.0 for Uzawa parameter'
            CALL Info( 'StressSolve', Message, Level=4 )
            UzawaParameter = 1.0d0
         END IF
         
         ContactPressure = MAX( 0.0d0, ContactPressure &
              + UzawaParameter * NormalDisplacement )
       END IF


       IF ( Iter > MinIter .AND. Solver % Variable % NonlinConverged == 1 ) EXIT

       IF ( ( CalcStressAll ) .AND. StabilityAnalysis ) THEN
         IF( Iter == 1 ) THEN
           CALL ComputeStress( Displacement, NodalStress,  &
               VonMises, DisplPerm, StressPerm, &
               NodalStrain, PrincipalStress, PrincipalStrain, Tresca, PrincipalAngle )
           
           CALL InvalidateVariable( Model % Meshes, Mesh, 'Stress' )
           CALL InvalidateVariable( Model % Meshes, Mesh, 'VonMises' )
           CALL InvalidateVariable( Model % Meshes, Mesh, 'Strain' )
           CALL InvalidateVariable( Model % Meshes, Mesh, 'Principal Stress' )
           CALL InvalidateVariable( Model % Meshes, Mesh, 'Principal Strain' )
           CALL InvalidateVariable( Model % Meshes, Mesh, 'Tresca' )
           CALL InvalidateVariable( Model % Meshes, Mesh, 'Principal Angle' )
         END IF
       END IF
       
       IF( ModelLumping ) THEN
         CALL LumpedSprings(iter,LumpedArea, LumpedCenter, LumpedMoments, &
             Model % MaxElementNodes)
       END IF
     END DO ! of nonlinear iter
!------------------------------------------------------------------------------

     IF ( CalcStressAll .AND. .NOT. StabilityAnalysis ) THEN
         
       IF ( EigenAnalysis ) THEN
         
         nsize = SIZE(Solver % Variable % EigenVectors,2)/STDOFs
         nomodes = Solver % NOFEigenValues

         DO i=1,nomodes

           WRITE (Message,'(A,I0)') 'Computing stresses for eigenmode: ',i 
           CALL INfo('StressSolver', Message ) 

           Displacement = Solver % Variable % EigenVectors(i,:)

           CALL ComputeStress( Displacement, NodalStress,  &
               VonMises, DisplPerm, StressPerm, &
               NodalStrain, PrincipalStress, PrincipalStrain, Tresca, PrincipalAngle )

           DO j=1,7               
             SELECT CASE ( j )
             CASE(1) 
               VarName = 'Stress'
             CASE(2) 
               VarName = 'vonMises'
             CASE(3) 
               VarName = 'Principal Stress'
             CASE(4) 
               VarName = 'Strain'
             CASE(5) 
               VarName = 'Principal Strain'
             CASE(6) 
               VarName = 'Principal Angle'
             CASE(7) 
               VarName = 'Tresca'                 
             END SELECT
             
             Var => VariableGet( Mesh % Variables, VarName )
             IF(.NOT. ASSOCIATED(Var) ) CYCLE
             dofs = Var % Dofs
            
             IF( i == 1 ) THEN               
               IF( .NOT. ASSOCIATED( Var % EigenVectors ) ) THEN
                 ALLOCATE( Var % EigenVectors(nomodes, dofs * nsize ) )             
                 Var % EigenVectors = 0.0_dp
               END IF
               IF( .NOT. ASSOCIATED( Var % EigenValues ) ) THEN                 
                 ALLOCATE( Var % EigenValues(nomodes) )
               END IF

               Var % EigenValues = Solver % Variable % EigenValues 
               IF( dofs > 1 ) THEN
                 DO k=1,dofs
                   iVar => VariableGet( Mesh % Variables,ComponentName(Var % Name,k) )
                   IF( ASSOCIATED( iVar ) ) THEN
                     iVar % EigenValues => Var % EigenValues
                     iVar % Eigenvectors => Var % EigenVectors(:,k::dofs)
                   ELSE
                     CALL Fatal('StressSolver','No variable associated: '//&
                         ComponentName( Var % Name,k ) )
                   END IF
                 END DO
               END IF
             END IF


             SELECT CASE ( j )

             CASE(1) 
               Var % EigenVectors(i,:) = NodalStress
             CASE(2) 
               Var % EigenVectors(i,:) = VonMises
             CASE(3) 
               Var % EigenVectors(i,:) = PrincipalStress
             CASE(4) 
               Var % EigenVectors(i,:) = NodalStrain
             CASE(5) 
               Var % EigenVectors(i,:) = PrincipalStrain
             CASE(6) 
               Var % EigenVectors(i,:) = PrincipalAngle
             CASE(7) 
               Var % EigenVectors(i,:) = Tresca
             END SELECT

           END DO
         END DO

       ELSE IF ( HarmonicAnalysis ) THEN

         nsize = SIZE(Solver % Variable % EigenVectors,2)/STDOFs
         nomodes = Solver % NOFEigenValues

         DO i=1,nomodes

           WRITE (Message,'(A,I0)') 'Computing stresses for eigenmode: ',i 
           CALL INfo('StressSolver', Message ) 

           DO l=1,2
            IF ( l==1 ) THEN
              Displacement = REAL( Solver % Variable % EigenVectors(i,:) )
            ELSE
              Displacement = AIMAG( Solver % Variable % EigenVectors(i,:) )
            END IF

            CALL ComputeStress( Displacement, NodalStress,  &
               VonMises, DisplPerm, StressPerm, &
               NodalStrain, PrincipalStress, PrincipalStrain, Tresca, PrincipalAngle )


           DO j=1,7               
             SELECT CASE ( j )
             CASE(1) 
               VarName = 'Stress'
             CASE(2) 
               VarName = 'vonMises'
             CASE(3) 
               VarName = 'Principal Stress'
             CASE(4) 
               VarName = 'Strain'
             CASE(5) 
               VarName = 'Principal Strain'
             CASE(6) 
               VarName = 'Principal Angle'
             CASE(7) 
               VarName = 'Tresca'                 
             END SELECT
             
             Var => VariableGet( Mesh % Variables, VarName )
             IF(.NOT. ASSOCIATED(Var) ) CYCLE
             dofs = Var % Dofs
            
             IF( i == 1 ) THEN               
               IF( .NOT. ASSOCIATED( Var % EigenVectors ) ) THEN
                 ALLOCATE( Var % EigenVectors(nomodes, dofs * nsize ) )             
                 Var % EigenVectors = 0.0_dp
               END IF
               IF( .NOT. ASSOCIATED( Var % EigenValues ) ) THEN                 
                 ALLOCATE( Var % EigenValues(nomodes) )
                 Var % EigenValues = 0._dp
               END IF

               Var % EigenValues = Solver % Variable % EigenValues 
               IF( dofs > 1 ) THEN
                 DO k=1,dofs
                   iVar => VariableGet( Mesh % Variables,ComponentName(Var % Name,k) )
                   IF( ASSOCIATED( iVar ) ) THEN
                     iVar % EigenValues => Var % EigenValues
                     iVar % Eigenvectors => Var % EigenVectors(:,k::dofs)
                   ELSE
                     CALL Fatal('StressSolver','No variable associated: '//&
                         ComponentName( Var % Name,k ) )
                   END IF
                 END DO
               END IF
             END IF


             IF ( l==1 ) THEN
               SELECT CASE ( j )
               CASE(1) 
                 Var % EigenVectors(i,:) = NodalStress
               CASE(2) 
                 Var % EigenVectors(i,:) = VonMises
               CASE(3) 
                 Var % EigenVectors(i,:) = PrincipalStress
               CASE(4) 
                 Var % EigenVectors(i,:) = NodalStrain
               CASE(5) 
                 Var % EigenVectors(i,:) = PrincipalStrain
               CASE(6) 
                 Var % EigenVectors(i,:) = PrincipalAngle
               CASE(7) 
                 Var % EigenVectors(i,:) = Tresca
               END SELECT
             ELSE
               SELECT CASE ( j )
               CASE(1) 
                 Var % EigenVectors(i,:) = Var % EigenVectors(i,:) + CMPLX(0._dp,NodalStress,KIND=dp)
               CASE(2) 
                 Var % EigenVectors(i,:) = Var % EigenVectors(i,:) + CMPLX(0._dp,VonMises,KIND=dp)
               CASE(3) 
                 Var % EigenVectors(i,:) = Var % EigenVectors(i,:) + CMPLX(0._dp,PrincipalStress,KIND=dp)
               CASE(4) 
                 Var % EigenVectors(i,:) = Var % EigenVectors(i,:) + CMPLX(0._dp,NodalStrain,KIND=dp)
               CASE(5) 
                 Var % EigenVectors(i,:) = Var % EigenVectors(i,:) + CMPLX(0._dp,PrincipalStrain,KIND=dp)
               CASE(6) 
                 Var % EigenVectors(i,:) = Var % EigenVectors(i,:) + CMPLX(0._dp,PrincipalAngle,KIND=dp)
               CASE(7) 
                 Var % EigenVectors(i,:) = Var % EigenVectors(i,:) + CMPLX(0._dp,Tresca,KIND=dp)
               END SELECT
             END IF
           END DO
           END DO
         END DO

       ELSE
         CALL ComputeStress( Displacement, NodalStress,  &
             VonMises, DisplPerm, StressPerm, &
             NodalStrain, PrincipalStress, PrincipalStrain, Tresca, PrincipalAngle )
       END IF

       CALL InvalidateVariable( Model % Meshes, Mesh, 'Stress' )
       CALL InvalidateVariable( Model % Meshes, Mesh, 'VonMises' )
       CALL InvalidateVariable( Model % Meshes, Mesh, 'Strain' )
       CALL InvalidateVariable( Model % Meshes, Mesh, 'Principal Stress' )
       CALL InvalidateVariable( Model % Meshes, Mesh, 'Principal Strain' )
       CALL InvalidateVariable( Model % Meshes, Mesh, 'Tresca' )
       CALL InvalidateVariable( Model % Meshes, Mesh, 'Principal Angle' )
     END IF

     IF ( GetLogical( SolverParams, 'Adaptive Mesh Refinement', Found) ) THEN
       CALL RefineMesh( Model, Solver, Displacement, DisplPerm, &
           StressInsideResidual, StressEdgeResidual, StressBoundaryResidual )
       
       IF ( MeshDisplacementActive ) THEN
         StressSol => Solver % Variable
         IF ( .NOT.ASSOCIATED( Mesh, Model % Mesh ) ) &
             CALL DisplaceMesh( Mesh, StressSol % Values, 1, &
             StressSol % Perm, StressSol % DOFs,.FALSE.)
       END IF
     END IF
 
     IF ( MeshDisplacementActive ) THEN
       CALL DisplaceMesh(Model % Mesh, Displacement, 1, &
           DisplPerm, STDOFs, .FALSE. )
     END IF

!------------------------------------------------------------------------------
!    Check where reference solution should be updated
!    Note: update determined by bodies
!------------------------------------------------------------------------------

     IF ( ANY( UpdatePresent ) ) THEN
       DO i = 1, Solver % NumberOfActiveElements
         Element => GetActiveElement(i)

         IF ( UpdateActive( Element % BodyId ) ) THEN
           n = GetElementNOFNodes( Element )
           DO j = 1, n
             indx = StressSol % Perm( Element % NodeIndexes(j) )
             DispValues( RefDofs * (indx - 1) + 1:  RefDofs * indx ) = &
                 StressSol % Values( RefDofs * (indx - 1) + 1:  RefDofs * indx )
           END DO
         END IF
       END DO
     END IF

     CALL Info('StressSolver','All done',Level=4)
     CALL Info('StressSolver','------------------------------------------',Level=4)

!------------------------------------------------------------------------------

CONTAINS
 
!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
    INTEGER :: RelIntegOrder 


     CALL StartAdvanceOutput( 'StressSolve', 'Assembly:')
     body_id = -1

     RelIntegOrder = ListGetInteger( SolverParams,'Relative Integration Order',Found)


     DO t=1,Solver % NumberOFActiveElements

!------------------------------------------------------------------------------
       CALL AdvanceOutput(t,GetNOFActive())
!------------------------------------------------------------------------------

       Element => GetActiveElement(t)
       n = GetElementNOFNOdes()
       ntot = GetElementNOFDOFs()

       NodeIndexes => Element % NodeIndexes
       CALL GetElementNodes( ElementNodes )

       Equation => GetEquation()
       PlaneStress = GetLogical( Equation, 'Plane Stress',Found )

       Material => GetMaterial()
       Density(1:n) = GetReal( Material, 'Density', Found )
       IF ( .NOT. Found )  THEN
         IF ( Transient .OR. EigenOrHarmonicAnalysis() ) &
            CALL Fatal( 'StressSolve', 'No value for density found.' )
       END IF

       Damping(1:n) = GetReal( Material, 'Damping', Found )
       RayleighDamping = GetLogical( Material, 'Rayleigh damping', Found )
       IF( RayleighDamping ) THEN
         RayleighAlpha(1:N) = GetReal( Material, 'Rayleigh alpha', Found )
         RayleighBeta(1:N) = GetReal( Material, 'Rayleigh beta', Found )
       ELSE
         RayleighAlpha = 0.0d0
         RayleighBeta = 0.0d0        
       END IF

       CALL InputTensor( HeatExpansionCoeff, Isotropic(2),  &
           'Heat Expansion Coefficient', Material, n, NodeIndexes, GotHeatExp )

       CALL InputTensor( ElasticModulus, Isotropic(1), &
           'Youngs Modulus', Material, n, NodeIndexes )

       PoissonRatio = 0.0d0
       IF ( Isotropic(1) )  PoissonRatio(1:n) = GetReal( Material, 'Poisson Ratio' )

       IF( GotHeatExp ) THEN
         ReferenceTemperature(1:n) = GetReal(Material, &
             'Reference Temperature', Found )
         CALL GetScalarLocalSolution( LocalTemperature, 'Temperature' )
         LocalTemperature(1:n) = LocalTemperature(1:n) - &
             ReferenceTemperature(1:n)
       ELSE
         LocalTemperature(1:n) = 0.0_dp
       END IF

       IF ( .NOT. ConstantBulkMatrixInUse ) THEN
         PreStress = 0.0d0
         PreStrain = 0.0d0
         CALL ListGetRealArray( Material, 'Pre Stress', Work, n, NodeIndexes, Found )
         IF ( Found ) THEN
            k = SIZE(Work,1)
            PreStress(1:k,1:n) = Work(1:k,1,1:n)
         END IF
         CALL ListGetRealArray( Material, 'Pre Strain', Work, n, NodeIndexes, Found )
         IF ( Found ) THEN
            k = SIZE(Work,1)
            PreStrain(1:k,1:n) = Work(1:k,1,1:n)
         END IF
       END IF

       ! Check need for elasticity matrix rotation:
       !-------------------------------------------
       IF ( Element % BodyId /= body_id ) THEN
         body_id = Element % BodyId
         RotateC = GetLogical( Material, 'Rotate Elasticity Tensor', stat )
       
         IF ( RotateC ) THEN
           RotateC = .FALSE.
           DO i=1,3
             IF( i == 1 ) THEN
               CALL GetConstRealArray( Material, UWrk, &
                   'Material Coordinates Unit Vector 1', stat, Element )
             ELSE IF( i == 2 ) THEN
               CALL GetConstRealArray( Material, UWrk, &
                   'Material Coordinates Unit Vector 2', stat, Element )
             ELSE                
               CALL GetConstRealArray( Material, UWrk, &
                   'Material Coordinates Unit Vector 3', stat, Element )
             END IF
             
             IF( stat ) THEN
               UnitNorm = SQRT( SUM( Uwrk(1:3,1)**2 ) )
               IF( UnitNorm < EPSILON( UnitNorm ) ) THEN
                 CALL Fatal('StressSolver','Given > Materia Coordinate Unit Vector < too short!')
               END IF
               TransformMatrix(i,1:3) = Uwrk(1:3,1) / UnitNorm  
               RotateC = .TRUE.
             ELSE 
               TransformMatrix(i,1:3) = 0.0_dp
               TransformMatrix(i,i) = 1.0_dp
             END IF
           END DO
           
           IF( .NOT. RotateC  ) THEN
             CALL Fatal( 'StressSolver', &
                 'No unit vectors found but > Rotate Elasticity Tensor < set True?' )
           END IF
         END IF
       END IF

       ! Set body forces:
       !-----------------
       BodyForce => GetBodyForce()
       LOAD = 0.0D0; LOAD_im=0._dp
       StressLoad = 0.0d0
       StrainLoad = 0.0d0
       IF ( ASSOCIATED( BodyForce ) ) THEN
         LOAD(1,1:n)  = GetReal( BodyForce, 'Stress Bodyforce 1', Found )
         LOAD(2,1:n)  = GetReal( BodyForce, 'Stress Bodyforce 2', Found )
         LOAD(3,1:n)  = GetReal( BodyForce, 'Stress Bodyforce 3', Found )
         LOAD(4,1:n)  = GetReal( BodyForce, 'Stress Pressure', Found )

         IF ( HarmonicAnalysis ) THEN
           LOAD_im(1,1:n)  = GetReal( BodyForce, 'Stress Bodyforce 1 im', Found )
           LOAD_im(2,1:n)  = GetReal( BodyForce, 'Stress Bodyforce 2 im', Found )
           LOAD_im(3,1:n)  = GetReal( BodyForce, 'Stress Bodyforce 3 im', Found )
           LOAD_im(4,1:n)  = GetReal( BodyForce, 'Stress Pressure im', Found )
         END IF

         CALL ListGetRealArray( BodyForce, 'Stress Load', Work, n, NodeIndexes, Found )
         IF ( Found ) THEN
            k = SIZE(Work,1)
            StressLoad(1:k,1:n) = Work(1:k,1,1:n)
         END IF

         CALL ListGetRealArray( BodyForce, 'Strain Load', Work, n, NodeIndexes, Found )
         IF ( Found ) THEN
            k = SIZE(Work,1)
            StrainLoad(1:k,1:n) = Work(1:k,1,1:n)
         END IF
       END IF

       ! Get element local stiffness & mass matrices:
       !---------------------------------------------

       IF ( .NOT. ConstantBulkMatrixInUse ) THEN
         CALL GetVectorLocalSolution( NodalDisplacement )
         IF ( Transient ) THEN
           NodalMeshVelo(3,1:n) = GetReal( Material, 'Mesh Velocity 3', Found)
           NodalMeshVelo(2,1:n) = GetReal( Material, 'Mesh Velocity 2', Found)
           NodalMeshVelo(1,1:n) = GetReal( Material, 'Mesh Velocity 1', Found)
           IF ( .NOT. Found ) THEN
             CALL GetVectorLocalSolution( NodalMeshVelo, 'Mesh Velocity' )
           END IF
         ELSE
           NodalMeshVelo   = 0.0d0
         END IF
       END IF

       SELECT CASE( CurrentCoordinateSystem() )
       CASE( Cartesian, AxisSymmetric, CylindricSymmetric )
          IF ( ConstantBulkMatrixInUse ) THEN
            CALL StressForceCompose( FORCE, FORCE_im, LOAD, LOAD_im, ElasticModulus, PoissonRatio, &
              PlaneStress, Isotropic,StressLoad, StrainLoad, HeatExpansionCoeff,         &
              LocalTemperature, Element, n, ntot, ElementNodes, RelIntegOrder, RotateC, TransformMatrix )
          ELSE
            CALL StressCompose( MASS, DAMP, STIFF, FORCE, FORCE_im, LOAD, LOAD_im, ElasticModulus,  &
               PoissonRatio, Density, PlaneStress, Isotropic,              &
               PreStress, PreStrain, StressLoad, StrainLoad, HeatExpansionCoeff,    &
               LocalTemperature, Element, n, ntot, ElementNodes, RelIntegOrder, StabilityAnalysis  &
               .AND. iter>1, GeometricStiffness .AND. iter>1, NodalDisplacement,    &
               RotateC, TransformMatrix, NodalMeshVelo, Damping, RayleighDamping,            &
               RayleighAlpha, RayleighBeta )
          END IF

       CASE DEFAULT
          CALL StressGeneralCompose( MASS, STIFF,FORCE, LOAD, ElasticModulus, &
             PoissonRatio,Density,PlaneStress,Isotropic,HeatExpansionCoeff,   &
             LocalTemperature, Element,n,ElementNodes )
       END SELECT
!------------------------------------------------------------------------------
!      If time dependent simulation, add mass matrix to global 
!      matrix and global RHS vector
!------------------------------------------------------------------------------
       IF ( .NOT. (ConstantBulkMatrix .OR. ConstantBulkSystem .OR. ConstantSystem) ) THEN
         IF ( Transient .AND. .NOT. EigenOrHarmonicAnalysis() ) THEN
            CALL Default2ndOrderTime( MASS, DAMP, STIFF, FORCE )
         END IF
       END IF

!------------------------------------------------------------------------------
!      Check if reference of displacement has been changed
!------------------------------------------------------------------------------
       IF ( ANY( UpdatePresent ) ) THEN
         IF ( ASSOCIATED( BodyForce ) ) THEN

           UpdateRef(1:n) = GetReal( BodyForce, 'Update Reference Displacement', Found )
           IF ( Found ) THEN
             UpdateActive(body_id) = .TRUE.
             IF ( COUNT( UpdateRef(1:n) < 0.0) > COUNT( UpdateRef(1:n) >= 0) )  THEN
               UpdateActive(body_id) = .FALSE.
               Ref_rhs = 0.0d0

               DO i = 1, n
                 DO j = 1, RefDofs
                   NodalRefD((i-1)*RefDofs+j) = ReferenceSol % Values &
                       ( RefDofs*(ReferenceSol % Perm(NodeIndexes(i))-1)+j)
                 END DO
               END DO

               DO i = 1, n*RefDofs
                 Ref_rhs(i) = SUM( STIFF(i,1:RefDofs*n) * NodalRefD(1:RefDofs*n) )
               END DO
               Force = Force + Ref_rhs
             END IF
           END IF
         END IF
       END IF

!------------------------------------------------------------------------------
!      Update global matrices from local matrices 
!------------------------------------------------------------------------------
       IF ( ConstantBulkMatrixInUse ) THEN
         CALL DefaultUpdateForce( FORCE )
         IF ( HarmonicAnalysis ) THEN
           SaveRHS => Solver % Matrix % RHS
           Solver % Matrix % RHS => Solver % Matrix % RHS_im
           CALL DefaultUpdateForce( FORCE_im )
           Solver % Matrix % RHS => SaveRHS
         END IF
       ELSE
         CALL DefaultUpdateEquations( STIFF, FORCE )
         IF ( HarmonicAnalysis ) THEN
           SaveRHS => Solver % Matrix % RHS
           Solver % Matrix % RHS => Solver % Matrix % RHS_im
           CALL DefaultUpdateForce( FORCE_im )
           Solver % Matrix % RHS => SaveRHS
         END IF

         IF ( EigenOrHarmonicAnalysis() .OR. &
           ConstantBulkMatrix .OR. ConstantBulkSystem .OR. ConstantSystem ) THEN
           CALL DefaultUpdateMass( MASS )
           CALL DefaultUpdateDamp( DAMP )
         END IF
       END IF

!------------------------------------------------------------------------------
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
   SUBROUTINE BCAssembly()
!------------------------------------------------------------------------------
     DO t = 1, Mesh % NumberOfBoundaryElements

       Element => GetBoundaryElement(t)
       IF ( .NOT. ActiveBoundaryElement() ) CYCLE
       
       IF( .NOT. PossibleFluxElement(Element) ) CYCLE

       n = GetElementNOFNodes()
       ntot = GetElementNOFDOFs()

       BC => GetBC()
       IF ( ASSOCIATED( BC ) ) THEN
!------------------------------------------------------------------------------
          CALL GetElementNodes( ElementNodes )

          LOAD  = 0.0d0
          Beta  = 0.0d0
          DampCoeff   = 0.0d0
          SpringCoeff = 0.0d0

          ! Force in given direction BC: \tau\cdot n = F:
          !----------------------------------------------
          GotForceBC = .FALSE.
          LOAD(1,1:n) = GetReal( BC, 'Force 1',Found )
          LOAD(2,1:n) = GetReal( BC, 'Force 2',Found )
          LOAD(3,1:n) = GetReal( BC, 'Force 3',Found )
          Beta(1:n) =  GetReal( BC, 'Normal Force',Found )

          LOAD_im=0._dp
          IF ( HarmonicAnalysis ) THEN
            LOAD_im(1,1:n) = GetReal( BC, 'Force 1 im',Found )
            LOAD_im(2,1:n) = GetReal( BC, 'Force 2 im',Found )
            LOAD_im(3,1:n) = GetReal( BC, 'Force 3 im',Found )
            Beta_im(1:n) =  GetReal( BC, 'Normal Force im',Found )
          END IF

          CALL ListGetRealArray( BC, 'Stress Load', Work, &
                  n, NodeIndexes, Found )
          StressLoad = 0.0d0
          IF ( Found ) THEN
             k = SIZE(Work,1)
             StressLoad(1:k,1:n) = Work(1:k,1,1:n)
          END IF

          DampCoeff(1:n) =  GetReal( BC, 'Damping', Found )
          SpringCoeff(1:n,1,1) =  GetReal( BC, 'Spring', NormalSpring )
          IF ( .NOT. NormalSpring ) THEN
            DO i=1,dim
              SpringCoeff(1:n,i,i) = GetReal( BC, ComponentName('Spring',i), Found)
            END DO

            DO i=1,dim
              DO j=1,dim
                IF (ListCheckPresent(BC,'Spring '//TRIM(i2s(i))//i2s(j) )) &
                  SpringCoeff(1:n,i,j)=GetReal( BC, 'Spring '//TRIM(i2s(i))//i2s(j), Found)
              END DO
            END DO
          END IF
          ContactLimit(1:n) =  GetReal( BC, 'Contact Limit', Found )

          IF(ModelLumping .AND. .NOT. FixDisplacement) THEN
            IF(GetLogical( BC, 'Model Lumping Boundary',Found )) THEN
              CALL LumpedLoads( iter, LumpedArea, LumpedCenter, LumpedMoments, Load )
            END IF
          END IF

!---------------------------------------------------------------------------
          IF( Contact ) THEN
             CALL GetScalarLocalSolution( LocalContactPressure, 'Contact Pressure' )
             Beta = Beta - LocalContactPressure
          END IF 

          NormalTangential = GetLogical( BC, 'Normal-Tangential ' // & 
                   GetVarName(Solver % Variable), Found )

          SELECT CASE( CurrentCoordinateSystem() )
          CASE( Cartesian, AxisSymmetric, CylindricSymmetric )
             CALL StressBoundary( STIFF,DAMP,FORCE,FORCE_im, LOAD, LOAD_im,   &
               SpringCoeff,NormalSpring,DampCoeff, Beta, Beta_im, StressLoad, &
                   NormalTangential, Element,n,ntot,ElementNodes )
          CASE DEFAULT
             DAMP = 0.0d0
             CALL StressGeneralBoundary( STIFF,FORCE, LOAD, SpringCoeff(1:n,1,1),Beta, &
                              Element,n,ElementNodes )
          END SELECT

          IF ( .NOT. (ConstantSystem .OR. ConstantBulkSystem .OR. ConstantBulkMatrix ) ) THEN
            IF ( Transient .AND. .NOT.EigenOrHarmonicAnalysis() )  THEN
               MASS = 0.0d0
               CALL Default2ndOrderTime( MASS, DAMP, STIFF, FORCE )
            END IF
          END IF

          ! Update global matrices from local matrices:;
          !---------------------------------------------
          CALL DefaultUpdateEquations( STIFF, FORCE )
          IF ( HarmonicAnalysis ) THEN
            SaveRHS => Solver % Matrix % RHS
            Solver % Matrix % RHS => Solver % Matrix % RHS_im
            CALL DefaultUpdateForce( FORCE_im )
            Solver % Matrix % RHS => SaveRHS
          END IF
          IF ( EigenOrHarmonicAnalysis() .OR. ConstantSystem.AND.Transient ) &
             CALL DefaultUpdateDamp( DAMP )
!------------------------------------------------------------------------------
         END IF
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE BCAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AddGlobalTime()
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,l,n
     REAL(KIND=dp) :: FORCE(1)
     REAL(KIND=dp), POINTER :: SaveValues(:) => NULL()
     REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:),MASS(:,:),DAMP(:,:),X(:),V(:),A(:)
     SAVE STIFF, MASS, DAMP, X, V, A

     IF ( .NOT.ASSOCIATED(Solver % Variable % Values, SaveValues) ) THEN
        IF ( ALLOCATED(STIFF) ) DEALLOCATE( STIFF,MASS,DAMP,V,X,A )

        n = 0
        DO i=1,Solver % Matrix % NumberOfRows
          n = MAX( n,Solver % Matrix % Rows(i+1)-Solver % Matrix % Rows(i) )
        END DO
        ALLOCATE( STIFF(1,n),MASS(1,n),DAMP(1,n),V(n),X(n),A(n) )

        SaveValues => Solver % Variable % Values
     END IF

     DO i=1,Solver % Matrix % NumberOFRows
       n = 0
       DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
         n = n + 1
         STIFF(1,n) = Solver % Matrix % Values(j)
         MASS(1,n)  = Solver % Matrix % MassValues(j)
         DAMP(1,n)  = Solver % Matrix % DampValues(j)

         X(n) = Solver % Variable % PrevValues(Solver % Matrix % Cols(j),3)
         V(n) = Solver % Variable % PrevValues(Solver % Matrix % Cols(j),4)
         A(n) = Solver % Variable % PrevValues(Solver % Matrix % Cols(j),5)
       END DO
       FORCE(1) = Solver % Matrix % RHS(i)
       Solver % Matrix % Force(i,1) = FORCE(1)
       CALL Bossak2ndOrder( n,dt,MASS,DAMP,STIFF,FORCE,X,V,A,Solver % Alpha )
       n = 0
       DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
         n = n + 1
         Solver % Matrix % Values(j) = STIFF(1,n)
       END DO
       Solver % Matrix % RHS(i) = FORCE(1)
     END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AddGlobalTime
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ComputeNormalDisplacement( Displacement, NormalDisplacement, &
       DisplPerm, STDOfs )
!------------------------------------------------------------------------------
    INTEGER :: DisplPerm(:), STDOfs
    REAL(KIND=dp) :: Displacement(:), NormalDisplacement(:)
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: MaxNodes = 100
    TYPE( Element_t ), POINTER :: Element
    TYPE( ValueList_t ), POINTER :: BC
    REAL( KIND=dp ) :: Normal(3), LocalDisplacement(3), U, V
    REAL( KIND=dp ) :: ContactLimit( MaxNodes )
    INTEGER :: i, j, k, t, n
    LOGICAL :: ContactBoundary, Found
    INTEGER, POINTER :: Visited(:)

    TYPE( Nodes_t ) :: ElementNodes
    SAVE ElementNodes

    ALLOCATE( Visited(SIZE(DisplPerm)) )
    Visited = 0

    NormalDisplacement = 0.0d0

    DO t = 1, Mesh % NumberOfBoundaryElements
       Element => GetBoundaryElement(t)
       IF ( .NOT. ActiveBoundaryElement() ) CYCLE
       
       IF( .NOT. PossibleFluxElement(Element) ) CYCLE

       n = GetElementNOFNodes()
       BC => GetBC()
       
       IF ( ASSOCIATED( BC ) ) THEN
          ContactBoundary = GetLogical( BC, 'Contact Boundary', Found ) 
          IF( .NOT.Found .OR. .NOT.ContactBoundary ) CYCLE
!------------------------------------------------------------------------------
          ContactLimit(1:n) =  GetReal( BC, 'Contact Limit', Found )
          IF( .NOT.Found ) ContactLimit = 9.9d9
             
          CALL GetElementNodes( ElementNodes )
          
          DO i = 1,n
             U = Element % TYPE % NodeU(i)
             V = Element % TYPE % NodeV(i)
             
             Normal = NormalVector( Element, ElementNodes, U, V, .TRUE. )    
             k = DisplPerm( Element % NodeIndexes(i) )
             
             LocalDisplacement = 0.0d0
             DO j = 1,STDOFs
                LocalDisplacement( j ) = Displacement( STDOFs*(k-1)+j )
             END DO
             
             NormalDisplacement( k ) = NormalDisplacement( k ) & 
                  + SUM( Normal(1:3) * LocalDisplacement(1:3) ) - ContactLimit(i)

             Visited( k ) = Visited( k ) + 1
             
          END DO
!------------------------------------------------------------------------------
       END IF
    END DO
!------------------------------------------------------------------------------
    WHERE( Visited >= 1 ) NormalDisplacement = NormalDisplacement / Visited

    DEALLOCATE( Visited )
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeNormalDisplacement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE ComputeStress( Displacement, NodalStress, &
              VonMises, DisplPerm, StressPerm, &
              NodalStrain, PrincipalStress, PrincipalStrain, Tresca, PrincipalAngle)
!------------------------------------------------------------------------------
     INTEGER :: DisplPerm(:)
     INTEGER, POINTER :: StressPerm(:)
     REAL(KIND=dp) :: VonMises(:), NodalStress(:), Displacement(:), &
                      NodalStrain(:), PrincipalStress(:), PrincipalStrain(:), &
                      Tresca(:), PrincipalAngle(:)
!------------------------------------------------------------------------------
     TYPE(Nodes_t) :: Nodes
     INTEGER :: n,nd
     TYPE(Element_t), POINTER :: Element

     INTEGER :: i,j,k,l,p,q, t, dim,sdim,elem, IND(9), BodyId,EqId
     LOGICAL :: stat, CSymmetry, Isotropic(2), UseMask, ContactOn
     INTEGER, POINTER :: Visited(:), Indexes(:), Permutation(:)
     REAL(KIND=dp) :: u,v,w,x,y,z,Strain(3,3),Stress(3,3),LGrad(3,3),detJ, &
          Young, Poisson, Ident(3,3), C(6,6), S(6), weight, st, Work(9), Principal(3), Relax
     REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:),  FORCE(:), ForceG(:), &
        SBasis(:,:), LocalDisplacement(:,:), MASS(:,:), SFORCE(:), SForceG(:)

     TYPE(Solver_t), POINTER :: StSolver

     LOGICAL :: FirstTime = .TRUE., OptimizeBW, GlobalBubbles, &
          Factorize, FoundFactorize, FreeFactorize, FoundFreeFactorize, &
          LimiterOn, SkipChange, FoundSkipChange

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     CHARACTER(LEN=MAX_NAME_LEN) :: eqname

     SAVE Nodes, StSolver, ForceG, Permutation, SForceG, Eqname, UseMask

     ! These variables are needed for Principal stress calculation
     ! they are quite small and allocated even if principal stress calculation
     ! is not requested
     REAL(KIND=dp) :: PriCache(3,3), PriTmp, PriW(3),PriWork(102)
     INTEGER       :: PriN=3, PriLWork=102, PriInfo=0
     REAL(KIND=dp) :: PriAngT1=0, PriAngT2=0, PriAngV(3)=0

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()

     CALL Info('StressSolver','------------------------------------------',Level=5)
     CALL Info('StressSolver','Starting Stress Computation',Level=5)

     ! Temporarily remove application of limiters as they are not needed
     ! for stress computation. 
     !-------------------------------------------------------------------
     LimiterOn = ListGetLogical( SolverParams,'Apply Limiter',Found)
     IF( LimiterOn ) THEN
       CALL ListAddLogical( SolverParams,'Apply Limiter',.FALSE.) 
     END IF
     ContactOn = ListGetLogical( SolverParams,'Apply Contact BCs',Found)
     IF( ContactOn ) THEN
       CALL ListAddLogical( SolverParams,'Apply Contact BCs',.FALSE.) 
     END IF

     CALL ListSetNameSpace('stress:')

     n = MAX( Mesh % MaxElementDOFs, Mesh % MaxElementNodes )
     ALLOCATE( Indexes(n), LocalDisplacement(3,n), &
         MASS(n,n), FORCE(6*n), &
         SFORCE(6*n), &
         Basis(n), dBasisdx(n,3) )

     IF ( FirstTime .OR. Mesh % Changed ) THEN
       IF ( FirstTime ) THEN
         ALLOCATE( StSolver )
       ELSE
         DEALLOCATE( ForceG, SForceG )
         CALL FreeMatrix( StSolver % Matrix )
       END IF

       StSolver = Solver
       StSolver % Variable => VariableGet( StSolver % Mesh % Variables, &
                  'StressTemp', ThisOnly=.TRUE. )
       IF ( ASSOCIATED( StSolver % Variable ) ) THEN
         Permutation => StSolver % Variable % Perm
       ELSE
         ALLOCATE( Permutation( SIZE(Solver % Variable % Perm) ) )
         Permutation = 0
       END IF

       OptimizeBW = GetLogical( StSolver % Values, 'Optimize Bandwidth', Found )
       IF ( .NOT. Found ) OptimizeBW = .TRUE.

       GlobalBubbles = GetLogical(SolverParams,'Bubbles in Global System',Found )
       IF (.NOT.Found ) GlobalBubbles = .TRUE.

       IF( ListGetLogicalAnyEquation( Model,'Calculate Stresses' ) ) THEN
         UseMask = .TRUE.
         eqname = 'Calculate Stresses'
       ELSE
         UseMask = .FALSE.
         eqname = TRIM( ListGetString( StSolver % Values,'Equation') )
       END IF
       StSolver % Matrix => CreateMatrix( Model, Solver, Mesh, Permutation, &
           1, MATRIX_CRS, OptimizeBW, eqname, GlobalBubbles=GlobalBubbles )

       ALLOCATE( StSolver % Matrix % RHS(StSolver % Matrix % NumberOfRows) )
       StSolver % Matrix % Comm = Solver % Matrix % Comm

       ALLOCATE( ForceG(StSolver % Matrix % NumberOfRows*6) )
       ALLOCATE( SForceG(StSolver % Matrix % NumberOfRows*6) )

       IF ( .NOT. ASSOCIATED( StSolver % Variable ) ) THEN
          CALL VariableAddVector( StSolver % Mesh % Variables, StSolver % Mesh, StSolver, &
                 'StressTemp', 1, Perm = StressPerm, Output=.FALSE. )
          StSolver % Variable => VariableGet( StSolver % Mesh % Variables, 'StressTemp' )
       END IF
       FirstTime = .FALSE.
     END IF

     Model % Solver => StSolver
     IF ( EigenAnalysis ) &
       CALL ListAddLogical( SolverParams, 'Eigen Analysis', .FALSE. )

     IF( HarmonicAnalysis ) &
       CALL ListAddLogical( SolverParams, 'Harmonic Analysis', .FALSE. )

     StSolver % NOFEigenValues=0

     Ident = 0.0d0
     DO i=1,3
        Ident(i,i) = 1.0d0
     END DO

     CSymmetry = CurrentCoordinateSystem() == AxisSymmetric .OR. &
                 CurrentCoordinateSystem() == CylindricSymmetric

     IND = (/ 1, 4, 6, 4, 2, 5, 6, 5, 3 /)

     Relax = GetCReal( StSolver % Values,'Nonlinear System Relaxation Factor', Found )
     IF ( .NOT. Found ) Relax = 1.0d0
     CALL ListAddConstReal( StSolver % Values,'Nonlinear System Relaxation Factor', 1.0d0 )

     NodalStress  = 0.0d0
     ForceG       = 0.0d0
     IF(CalculateStrains) THEN
       NodalStrain  = 0.0d0
       SForceG      = 0.0d0
     END IF
     CALL DefaultInitialize()

     DO elem = 1,Solver % NumberOfActiveElements
        Element => GetActiveElement(elem, Solver)
        n  = GetElementNOFNodes()
        nd = GetElementDOFs( Indexes )

        Equation => GetEquation()

        ! Check if stresses wanted for this body:
        ! ---------------------------------------
        IF( UseMask ) THEN
          IF(.NOT. GetLogical( Equation, eqname, Found )) CYCLE
        END IF

        ! Get material parameters:
        ! ------------------------
        Material => GetMaterial()

        CALL InputTensor( HeatExpansionCoeff, Isotropic(2),  &
            'Heat Expansion Coefficient', Material, n, Element % NodeIndexes, GotHeatExp )

        CALL InputTensor( ElasticModulus, Isotropic(1), &
                'Youngs Modulus', Material, n, Element % NodeIndexes )
        PlaneStress = ListGetLogical( Equation, 'Plane Stress', stat )
        PoissonRatio(1:n) = GetReal( Material, 'Poisson Ratio', Stat )

        ! Element nodal points:
        ! ---------------------
        CALL GetElementNodes( Nodes )

        ! Displacement field at element nodal points:
        ! -------------------------------------------
        CALL GetVectorLocalSolution( LocalDisplacement, USolver=Solver )

        IF( GotHeatExp ) THEN
          ReferenceTemperature(1:n) = GetReal(Material, 'Reference Temperature', Found )
          CALL GetScalarLocalSolution( LocalTemperature, 'Temperature', USolver=Solver )
          LocalTemperature(1:n) = LocalTemperature(1:n) - ReferenceTemperature(1:n)
        ELSE
          LocalTemperature(1:n) = 0.0_dp
        END IF

        ! Integrate local stresses:
        ! -------------------------
        IntegStuff = GaussPoints( Element )
        Stress = 0.0d0
        Strain  = 0.0d0
        MASS   = 0.0d0
        FORCE  = 0.0d0
        SFORCE = 0.0d0

        DO t=1,IntegStuff % n
          u = IntegStuff % u(t)
          v = IntegStuff % v(t)
          w = IntegStuff % w(t)
          Weight = IntegStuff % s(t)

          stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
             Basis, dBasisdx )

          Weight = Weight * detJ
          IF ( CSymmetry ) Weight = Weight * SUM( Basis(1:n) * Nodes % x(1:n) )

          CALL LocalStress( Stress, Strain, PoissonRatio, &
              ElasticModulus, HeatExpansionCoeff, LocalTemperature, &
              Isotropic, CSymmetry, PlaneStress, LocalDisplacement, &
              Basis, dBasisdx, Nodes, dim, n, nd )

          DO p=1,nd
            DO q=1,nd
              MASS(p,q) = MASS(p,q) + Weight*Basis(q)*Basis(p)
            END DO

            DO i=1,3
            DO j=i,3
              k = Ind( 3*(i-1)+j )
              FORCE(6*(p-1)+k) = FORCE(6*(p-1)+k) + Weight*Stress(i,j)*Basis(p)
              SFORCE(6*(p-1)+k) = SFORCE(6*(p-1)+k) + Weight*Strain(i,j)*Basis(p)                  
            END DO
            END DO
          END DO
        END DO

        CALL DefaultUpdateEquations( MASS, FORCE )

        DO p=1,nd
          l = Permutation(Indexes(p))
          DO i=1,3
          DO j=i,3
             k = Ind(3*(i-1)+j)
             ForceG(6*(l-1)+k) = ForceG(6*(l-1)+k) + FORCE(6*(p-1)+k)
             SForceG(6*(l-1)+k) = SForceG(6*(l-1)+k) + SFORCE(6*(p-1)+k)
          END DO
          END DO
        END DO
      END DO

      Factorize = GetLogical( SolverParams, 'Linear System Refactorize', FoundFactorize )
      FreeFactorize = GetLogical( SolverParams, &
          'Linear System Free Factorization', FoundFreeFactorize )
      SkipChange = GetLogical( SolverParams, &
          'Skip Compute Nonlinear Change', FoundSkipChange )

      CALL ListAddLogical( SolverParams, 'Linear System Refactorize', .FALSE. )
      CALL ListAddLogical( SolverParams, 'Linear System Free Factorization', .FALSE. )
      CALL ListAddLogical( SolverParams, 'Skip Compute Nonlinear Change', .TRUE. )

      DO i=1,3
        DO j=i,3
          k = IND(3*(i-1)+j)
          
          StSolver % Matrix % RHS = ForceG(k::6)
          
          DO l=1,SIZE( Permutation )
            IF ( Permutation(l) <= 0 ) CYCLE
            StSolver % Variable % Values(Permutation(l)) = NodalStress(6*(StressPerm(l)-1)+k)
          END DO
          
          WRITE( Message,'(A,I0,A,I0,A)') 'Solving for Stress(',i,',',j,')'
          CALL Info('StressSolver',Message,Level=5)

          st = DefaultSolve()
           
          DO l=1,SIZE( Permutation )
            IF ( Permutation(l) <= 0 ) CYCLE
            NodalStress(6*(StressPerm(l)-1)+k) = StSolver % Variable % Values(Permutation(l))
          END DO
          
          IF(CalculateStrains) THEN
            StSolver % Matrix % RHS = SForceG(k::6)
            DO l=1,SIZE( Permutation )
              IF ( Permutation(l) <= 0 ) CYCLE
              StSolver % Variable % Values(Permutation(l)) = NodalStrain(6*(StressPerm(l)-1)+k)            
            END DO
            ! this solves some convergence problems at the expence of bad convergence      
            ! StSolver % Variable % Values = 0

            WRITE( Message,'(A,I0,A,I0,A)') 'Solving for Strain(',i,',',j,')'
            CALL Info('StressSolver',Message,Level=5)
            st = DefaultSolve()
          
            DO l=1,SIZE( Permutation )
              IF ( Permutation(l) <= 0 ) CYCLE
              NodalStrain(6*(StressPerm(l)-1)+k) = StSolver % Variable % Values(Permutation(l))
            END DO
          END IF !CalculateStrains
        END DO
      END DO

      IF ( FoundFactorize ) THEN
        CALL ListAddLogical( SolverParams, 'Linear System Refactorize', Factorize )
      ELSE
        CALL ListRemove( SolverParams, 'Linear System Refactorize' )
      END IF

      IF ( FoundFreeFactorize ) THEN
        CALL ListAddLogical( SolverParams, 'Linear System Free Factorization', FreeFactorize )
      ELSE
        CALL ListRemove( SolverParams, 'Linear System Free Factorization' )
      END IF

      IF( FoundSkipChange ) THEN
        CALL ListAddLogical( SolverParams, 'Skip Compute Nonlinear Change',SkipChange )
      ELSE
        CALL ListRemove( SolverParams, 'Skip Compute Nonlinear Change' )
      END IF

      ! Von Mises stress from the component nodal values:
      ! -------------------------------------------------
      VonMises = 0
      DO i=1,SIZE( StressPerm )
         IF ( StressPerm(i) <= 0 ) CYCLE

         p = 0
         DO j=1,3
            DO k=1,3
              p = p + 1
              q = 6 * (StressPerm(i)-1) + IND(p)
              Stress(j,k) = NodalStress(q)
            END DO
         END DO

         Stress(:,:) = Stress(:,:) - TRACE(Stress(:,:),3) * Ident/3

         DO j=1,3
            DO k=1,3
              VonMises(StressPerm(i)) = VonMises(StressPerm(i)) + Stress(j,k)**2
            END DO
         END DO
      END DO

      VonMises = SQRT( 3.0d0 * VonMises / 2.0d0 )

      !Principal stresses and Tresca
      IF(CalcPrincipalAll) THEN
        DO i=1,SIZE( StressPerm )
          IF ( StressPerm(i) <= 0 ) CYCLE       
          !Stresses: 
          p = 0

          sdim=3
          IF (dim==2.AND.PlaneStress) sdim=2

          DO j=1,3
            DO k=1,3 ! TODO only upper triangle should be filled, this is is wasteful
              p = p+1
              q = 6 * (StressPerm(i)-1) + IND(p)
              PriCache(j,k) = NodalStress(q)
            END DO
          END DO

          !Use lapack function to do solve eigenvalues (i.e. principal stresses)
          CALL DSYEV( 'N', 'U', sdim, PriCache, 3, PriW, PriWork, PriLWork, PriInfo )
          IF (PriInfo /= 0) THEN !error in dsyev
            PriW = 0; !we probably should put NaN in error
          END IF              

          DO l=1,sdim
            ! eigenvalues are returned in opposite order 
            PrincipalStress(3 * (StressPerm(i)-1 )+l) = PriW(sdim+1-l)
          END DO

          IF(CalcPrincipalAngle) THEN
            !DSYEV has changed the vector, so well copy it again from NodalStress
            p=0
            DO j=1,3
              DO k=1,3 ! TODO only upper triangle should be filled, this is is wasteful
                 p = p+1
                 q = 6 * (StressPerm(i)-1) + IND(p)
                 PriCache(j,k) = NodalStress(q)
              END DO
            END DO

            DO k=1,3 ! for all principal stresses
              ! This is where things get _very_ heary. The code below
              ! solves the following equation system:
              !   (s11-p)v1 + s12*v2 + s13*v3=0
              !   s12*v1    + (s22-p)*v2 + s23*v3=0
              !   v1**2 + v2**2 + v3**2 = 1
              !   where v1...3 are the directional cosines of the primary stresses,
              !   sij are the stress matrix components and
              !   p   is the primary stress in question  
              ! The code is practically unreadable.
              ! This code unit has been tested with the non-trivial
              ! known solutions from following textbooks:
              !   Pennala, E. 1992. Lujuusopin Perusteet
              !   Shames, I. H., Cozzarelli, F., A. Elastic and inelastic stress analysis
              PriAngT1 =( PriCache(1,2) + &
                 ((PriCache(2,2)-PriW(4-k)) * (PriCache(1,1)-PriW(4-k))) / &
                 (-PriCache(1,2) ) ) / & 
                (((PriCache(2,2)-PriW(4-k))*PriCache(1,3))/&
                  PriCache(1,2) - PriCache(2,3) )
              PriAngT2 =(PriCache(1,1)-PriW(4-k))/(PriCache(1,2)) + &
                     PriAngT1*PriCache(1,3)/PriCache(1,2)
              PriAngV(1) = (1/(1+PriAngT1**2 + PriAngT2**2))**0.5_dp
              PriAngV(2) = -PriAngT2 * PriAngV(1)
              PriAngV(3) =  PriAngT1 * PriAngV(1)

              PrincipalAngle(9 * (StressPerm(i)-1 ) +3*(k-1) + 1)    = &
                          (ACOS(PriAngV(1)) ) !angle in radians *360/6.28
              PrincipalAngle(9 * (StressPerm(i)-1 ) +3*(k-1) + 2) = &
                          (ACOS(PriAngV(2)) ) !angle in radians *360/6.28
              PrincipalAngle(9 * (StressPerm(i)-1 ) +3*(k-1) + 3) = &
                          (ACOS(PriAngV(3)) ) !angle in radians *360/6.28                              
            END DO
          END IF

          !Tresca                        
          Tresca(StressPerm(i)) = (PrincipalStress(3*(StressPerm(i)-1) +1) - &
                       PrincipalStress(3*(StressPerm(i)-1) +2))/2
          PriTmp = (PrincipalStress(3*(StressPerm(i)-1) +2) - &
                       PrincipalStress(3*(StressPerm(i)-1) +3))/2
          IF (PriTmp > Tresca(StressPerm(i)) ) Tresca(StressPerm(i)) = PriTmp

          PriTmp = (PrincipalStress(3*(StressPerm(i)-1) +1) - &
                       PrincipalStress(3*(StressPerm(i)-1) +3))/2
          IF (PriTmp > Tresca(StressPerm(i)) ) Tresca(StressPerm(i)) = PriTmp
          
          !Strain:
          IF(CalculateStrains)THEN
            p=0
            DO j=1,3
              DO k=1,3 ! TODO only upper triangle should be filled, this is is wasteful
                p = p+1
                q = 6 * (StressPerm(i)-1) + IND(p)
                PriCache(j,k) = NodalStrain(q)
              END DO
            END DO
      
            sdim=3; IF(dim==2.AND..NOT.PlaneStress) sdim=2

            !Use lapack function to do solve eigenvalues
            CALL DSYEV( 'N', 'U', sdim, PriCache, 3, PriW, PriWork, PriLWork, PriInfo )
            IF(PriInfo /= 0) PriW = 0;
            DO l=1,sdim
              ! eigenvalues are returned in opposite order 
              PrincipalStrain(3 * (StressPerm(i)-1 )+l) = PriW(sdim+1-l)
            END DO
          END IF ! CalculateStrains
        END DO
      END IF ! Calculate Principal

      DEALLOCATE( Basis, dBasisdx )
      DEALLOCATE( Indexes, LocalDisplacement, MASS, FORCE )

      IF ( EigenAnalysis ) &
        CALL ListAddLogical( SolverParams, 'Eigen Analysis', .TRUE. )
      IF ( HarmonicAnalysis ) &
        CALL ListAddLogical( SolverParams, 'Harmonic Analysis', .TRUE. )
      CALL ListAddConstReal( SolverParams,'Nonlinear System Relaxation Factor', Relax )


      Model % Solver => Solver

      IF( LimiterOn ) THEN
        CALL ListAddLogical( SolverParams,'Apply Limiter',.TRUE.) 
      END IF
      IF( ContactOn ) THEN
        CALL ListAddLogical( SolverParams,'Apply Contact BCs',.TRUE.) 
      END IF

      CALL Info('StressSolver','Finished Stress Computation',Level=5)
      CALL Info('StressSolver','------------------------------------------',Level=5)

      CALL ListSetNameSpace('')

!------------------------------------------------------------------------------
   END SUBROUTINE ComputeStress
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION TRACE( F, dim ) RESULT(t)
!------------------------------------------------------------------------------
     INTEGER :: i, dim
     REAL(KIND=dp) :: F(:,:), t

     t = 0.0d0
     DO i=1,dim
        t = t + F(i,i)
     END DO
!------------------------------------------------------------------------------
   END FUNCTION TRACE
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Computes area, center of area and different moments
!------------------------------------------------------------------------------
   SUBROUTINE CoordinateIntegrals(Area, Center, Moments, maxnodes)

     REAL(KIND=dp) :: Area, Center(:), Moments(:,:)
     INTEGER :: maxnodes
     LOGICAL :: FoundBoundary

     REAL(KIND=dp) :: Coords(3)
     INTEGER :: power
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp) :: Basis(maxnodes)
     REAL(KIND=dp) :: dBasisdx(maxnodes,3),detJ,u,v,w
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     INTEGER :: N_Integ
     LOGICAL :: stat

     FoundBoundary = .FALSE.
     Area = 0.0
     Center = 0.0
     Moments = 0.0


     ! On the first round compute area and center of area.
     ! On the second round compute the square deviations from the mean.
     
     DO power = 1,2

       DO t=1,Mesh % NumberOfBoundaryElements
         Element => GetBoundaryElement(t)
         IF ( .NOT. ActiveBoundaryElement() ) CYCLE

         IF( .NOT. PossibleFluxElement(Element) ) CYCLE

         BC => GetBC()
         IF ( .NOT.ASSOCIATED( BC ) ) CYCLE
!------------------------------------------------------------------------------
         IF(.NOT. GetLogical( BC, 'Model Lumping Boundary',Found )) CYCLE
         
         FoundBoundary = .TRUE.
         n = GetElementNOFNodes()
         CALL GetElementNodes( ElementNodes )

         IntegStuff = GaussPoints( Element )
         U_Integ => IntegStuff % u
         V_Integ => IntegStuff % v
         W_Integ => IntegStuff % w
         S_Integ => IntegStuff % s
         N_Integ =  IntegStuff % n
         
         DO k=1,N_Integ
           u = U_Integ(k)
           v = V_Integ(k)
           w = W_Integ(k)
           
           ! Basis function values & derivatives at the integration point:
           !--------------------------------------------------------------
           stat = ElementInfo( Element, ElementNodes, u, v, w, detJ, &
               Basis, dBasisdx )
           
           s = detJ * S_Integ(k)
           IF ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
               CurrentCoordinateSystem() == CylindricSymmetric ) THEN
             s = s * SUM( ElementNodes % x(1:n) * Basis(1:n) )
           END IF
           
           Coords(1) = SUM(Basis(1:n) * ElementNodes % x(1:n))
           IF (DIM > 1) THEN
             Coords(2) =  SUM(Basis(1:n) * ElementNodes % y(1:n))
           END IF
           IF (DIM > 2) THEN
             Coords(3) =  SUM(Basis(1:n) * ElementNodes % z(1:n))
           END IF
           
           IF(power == 1) THEN
             Area = Area + s
             Center(1:DIM) = Center(1:DIM) + s * Coords(1:DIM)
           ELSE
             Coords(1:DIM) = Coords(1:DIM) - Center(1:DIM) 
             DO i = 1,DIM
               DO j = 1,DIM
                 Moments(i,j) = Moments(i,j) + s * Coords(i) * Coords(j)
               END DO
             END DO
           END IF

         END DO
       END DO
         
       IF(.NOT. FoundBoundary) THEN
        CALL Fatal('StressSolve','Model lumping boudary must be defined')        
       END IF
   
       IF(power == 1) Center(1:DIM) = Center(1:DIM) / Area
     END DO

   END SUBROUTINE CoordinateIntegrals


!------------------------------------------------------------------------------
! Compute the loads resulting to pure forces or pure moments.
! Pure moments may only be computed under certain conditions that 
! should be valid for boundaries with normal in the direction of some axis.
!------------------------------------------------------------------------------

   SUBROUTINE LumpedLoads( Permutation, Area, Center, Moments, Forces )
     INTEGER :: Permutation
     REAL (KIND=dp) :: Area, Center(:), Moments(:,:), Forces(:,:)
     
     REAL (KIND=dp), POINTER :: y(:), z(:)
     REAL (KIND=dp) :: c, Eps
     LOGICAL :: isy, isz
     INTEGER :: ix,iy,iz,nx,ny,nz

     Forces = 0.0d0
     Eps = 1.0d-6

     IF(Permutation <= 3) THEN
       Forces(Permutation,1:n) = 1.0 / LumpedArea
     ELSE IF(Permutation <= 6) THEN
       ix = MOD(Permutation - 4, 3) + 1
       iy = MOD(Permutation - 3, 3) + 1
       iz = MOD(Permutation - 2, 3) + 1

       IF(Permutation == 4) THEN
         z => ElementNodes % Z
         y => ElementNodes % Y
       ELSE IF(Permutation == 5) THEN
         z => ElementNodes % X
         y => ElementNodes % Z
       ELSE IF(Permutation == 6) THEN
         z => ElementNodes % Y
         y => ElementNodes % X
       END IF

       isy = (ABS(Moments(iy,ix)) < Eps * Moments(iy,iy))
       isz = (ABS(Moments(iz,ix)) < Eps * Moments(iz,iz))

       IF(isy) THEN
         c = 1.0 / Moments(iy,iy)
         Forces(iz,1:n) = c * (y(1:n) - Center(iy))
       ELSE IF(isz) THEN
         c = -1.0 / Moments(iz,iz)
         Forces(iy,1:n) = c * (z(1:n) - Center(iz))
       ELSE 
         c = 1.0 / (Moments(iy,iy) + Moments(iz,iz) )
         Forces(iy,1:n) = -c * (z(1:n) - Center(iz))
         Forces(iz,1:n) =  c * (y(1:n) - Center(iy))
         CALL Warn('StressSolve','Moment matrix not diagonalazible!')
         PRINT *,Moments(iy,ix),Moments(iz,ix),Moments(iy,iy),Moments(iz,iz)
       END IF
     END IF
   END SUBROUTINE LumpedLoads


!------------------------------------------------------------------------------
   SUBROUTINE LumpedDisplacements( Model, Permutation, Area, Center )
!------------------------------------------------------------------------------
!  This subroutine is used to set pure translations and rotations to the 
!  chosen boundary in order to perform model lumping using fixed displacement.
!------------------------------------------------------------------------------

     TYPE(Model_t) :: Model
     REAL(KIND=dp) :: Area, Center(:)
     INTEGER :: Permutation
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix
     REAL(KIND=dp), POINTER :: ForceVector(:)
     INTEGER, POINTER :: Perm(:)
     TYPE(Element_t), POINTER :: CurrentElement
     INTEGER, POINTER :: NodeIndexes(:)
     INTEGER :: i,j,k,l,n,t,ind
     LOGICAL :: GotIt
     REAL(KIND=dp) :: Coords(3), dCoords(3), dFii, dx, s
    
    !------------------------------------------------------------------------------
    
     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS
     Perm => Solver % Variable % Perm
     
     dX   = 1.0d-2*SQRT(Area)
     dFii = 1.0d-2
     
     DO t = 1, Mesh % NumberOfBoundaryElements
       Element => GetBoundaryElement(t)
       CurrentElement => Element
       IF ( .NOT. ActiveBoundaryElement()) CYCLE
       n = GetElementNOFNodes()
       
       BC => GetBC()
       IF ( .NOT.ASSOCIATED( BC ) ) CYCLE
       
       IF(.NOT. GetLogical( BC, 'Model Lumping Boundary',Found )) CYCLE

       NodeIndexes => CurrentElement % NodeIndexes
       
       DO j=1,n
         k = Perm(NodeIndexes(j))
         IF(k == 0) CYCLE
         
         dCoords = 0.0d0
         IF(Permutation <= 3) THEN
           dCoords(Permutation) = dX
         ELSE
           Coords(1) = Mesh % Nodes % x(NodeIndexes(j))
           Coords(2) = Mesh % Nodes % y(NodeIndexes(j))
           Coords(3) = Mesh % Nodes % z(NodeIndexes(j))
           Coords = Coords - Center
           IF (Permutation == 4) THEN
             dCoords(2) = -dFii * Coords(3) 
             dCoords(3) = dFii * Coords(2)
           ELSE IF(Permutation == 5) THEN
             dCoords(1) = dFii * Coords(3) 
             dCoords(3) = -dFii * Coords(1)
           ELSE IF(Permutation == 6) THEN
             dCoords(1) = -dFii * Coords(2)
             dCoords(2) = dFii * Coords(1)
           END IF

        END IF

         DO l=1,dim
           CALL SetDirichletPoint( StiffMatrix, ForceVector, l, dim, Perm, NodeIndexes(j), dCoords(l) )
         END DO
       END DO
     END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LumpedDisplacements
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! At the end of each iteration assemblys one line of the Kmatrix and finally 
! invert the matrix. The displacements and the springs are taken to be the 
! average values on the surface.
!------------------------------------------------------------------------------
   SUBROUTINE LumpedSprings(Permutation,Area, Center, Moments, maxnodes)
!------------------------------------------------------------------------------
     INTEGER :: Permutation, maxnodes     
     REAL(KIND=dp) :: Area, Center(:), Moments(:,:)
!------------------------------------------------------------------------------
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp) :: Basis(maxnodes)
     REAL(KIND=dp) :: dBasisdx(maxnodes,3),detJ,u,v,w
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     REAL(KIND=dp) :: LocalDisp(DIM,maxnodes),Kmat(6,6), up, vp, wp, &
         xp(maxnodes), yp(maxnodes), zp(maxnodes), KmatMin(6,6), KvecAtIP(6), &
         Strain(3,3),Stress(3,3), dFii, Dx, &
         ForceAtIp(3), MomentAtIp(3), Coord(3),Normal(3)
     REAL(KIND=dp), POINTER :: PValues(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalLoads(:)
     LOGICAL, POINTER :: NodeVisited(:)
     INTEGER :: N_Integ, pn
     INTEGER, POINTER :: Indexes(:)
     LOGICAL :: stat, CSymmetry, Isotropic
     CHARACTER(LEN=MAX_NAME_LEN) :: KmatFile
     TYPE(Nodes_t) :: ParentNodes
     TYPE(Element_t),POINTER :: Parent

     SAVE ParentNodes, Kmat, KmatMin, NodalLoads, NodeVisited
!------------------------------------------------------------------------------

     n = maxnodes
     ALLOCATE( ParentNodes % x(n), ParentNodes % y(n), ParentNodes % z(n))

     CSymmetry = CurrentCoordinateSystem() == CylindricSymmetric .OR. &
                 CurrentCoordinateSystem() == AxisSymmetric
    
     dFii = 1.0d-2
     dX = 1.0d-2*SQRT(Area)

     IF (Permutation == 1) THEN
       Kmat = 0.0d0       
       KmatMin = HUGE(KmatMin)
     END IF

     IF( FixDisplacement ) THEN
       IF(Permutation == 1) THEN
         n = SIZE( Displacement ) / STDOFs
         ALLOCATE( NodalLoads( STDOFs * n ), NodeVisited( n ) )
       END IF
       
       NodalLoads = 0.0d0
       PValues => Solver % Matrix % Values
       Solver % Matrix % Values => Solver % Matrix % BulkValues
       CALL MatrixVectorMultiply( Solver % Matrix, Displacement, NodalLoads)
       Solver % Matrix % Values => PValues
       
       NodeVisited = .FALSE.
       
       DO t = 1, Mesh % NumberOfBoundaryElements
         Element => GetBoundaryElement(t)
         IF ( .NOT. ActiveBoundaryElement() ) CYCLE

         IF( .NOT. PossibleFluxElement(Element) ) CYCLE
         
         BC => GetBC()
         IF ( .NOT.ASSOCIATED( BC ) ) CYCLE
         IF(.NOT. GetLogical( BC, 'Model Lumping Boundary',Found )) CYCLE
         
         n = GetElementNOFNodes()
         CALL GetElementNodes( ElementNodes )
         Indexes => Element % NodeIndexes
         
         DO i=1,n
           j = DisplPerm( Indexes(i) )
           IF(NodeVisited(j)) CYCLE
           NodeVisited(j) = .TRUE.
           
           Coord(1) = ElementNodes % x(i)
           Coord(2) = ElementNodes % y(i)
           Coord(3) = ElementNodes % z(i)
           Coord = Coord - Center        
           
           DO k=1,DIM
             ForceAtIP(k) = NodalLoads(3*(j-1)+k)
           END DO

           MomentAtIp(1) = -ForceAtIp(2) * Coord(3) + ForceAtIp(3) * Coord(2)
           MomentAtIp(2) = -ForceAtIp(3) * Coord(1) + ForceAtIp(1) * Coord(3)
           MomentAtIp(3) = -ForceAtIp(1) * Coord(2) + ForceAtIp(2) * Coord(1)
           
           Kmat(1:3,Permutation) = Kmat(1:3,Permutation) + ForceAtIp 
           Kmat(4:6,Permutation) = Kmat(4:6,Permutation) + MomentAtIp
         END DO
       END DO


     ELSE
       DO t = 1, Mesh % NumberOfBoundaryElements
         Element => GetBoundaryElement(t)
         IF ( .NOT. ActiveBoundaryElement() ) CYCLE

         IF( .NOT. PossibleFluxElement(Element) ) CYCLE
         
         BC => GetBC()
         IF ( .NOT.ASSOCIATED( BC ) ) CYCLE
         IF(.NOT. GetLogical( BC, 'Model Lumping Boundary',Found )) CYCLE
         
         n = GetElementNOFNodes()
         CALL GetElementNodes( ElementNodes )
         
         ! Get parent element & nodes:
         ! ---------------------------
         Parent => Element % BoundaryInfo % Left
         stat = ASSOCIATED( Parent )
         IF ( .NOT. stat ) stat = ALL(DisplPerm(Parent % NodeIndexes) > 0)
         IF ( .NOT. stat ) THEN
           Parent => Element % BoundaryInfo % Right
           stat = ASSOCIATED( Parent )
           IF ( stat ) stat = ALL(DisplPerm(Parent % NodeIndexes) > 0)
           IF ( .NOT. stat ) CALL Fatal( 'StressSolve', & 
               'Cannot find proper parent for side element' )
         END IF
         pn = GetElementNOFNodes( Parent )
         CALL GetElementNodes( ParentNodes, Parent )
         CALL GetVectorLocalSolution( LocalDisp, UElement=Parent )
         
         ! Get boundary nodal points in parent local coordinates:
         ! ------------------------------------------------------
         DO i = 1,n
           DO j = 1,pn
             IF ( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
               xp(i) = Parent % TYPE % NodeU(j)
               yp(i) = Parent % TYPE % NodeV(j)
               zp(i) = Parent % TYPE % NodeW(j)
               EXIT
             END IF
           END DO
         END DO
         
         IntegStuff = GaussPoints( Element )

         U_Integ => IntegStuff % u
         V_Integ => IntegStuff % v
         W_Integ => IntegStuff % w
         S_Integ => IntegStuff % s
         N_Integ =  IntegStuff % n
         
         DO k=1,N_Integ
           u = U_Integ(k)
           v = V_Integ(k)
           w = W_Integ(k)
           
           ! Basis function values & derivatives at the integration point:
           !--------------------------------------------------------------
           stat = ElementInfo( Element, ElementNodes, u, v, w, detJ, &
               Basis, dBasisdx )
           
           s = detJ * S_Integ(k)
           IF ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
               CurrentCoordinateSystem() == CylindricSymmetric ) THEN
             s = s * SUM( ElementNodes % x(1:n) * Basis(1:n) )
           END IF
           
           ! The plane  elements only include the  derivatives in the direction
           ! of the plane. Therefore compute the derivatives of the displacemnt
           ! field from the parent element:
           ! -------------------------------------------------------------------
           Up = SUM( xp(1:n) * Basis(1:n) )
           Vp = SUM( yp(1:n) * Basis(1:n) )
           Wp = SUM( zp(1:n) * Basis(1:n) )
           
           stat = ElementInfo( Parent,ParentNodes, Up, Vp, Wp, detJ, &
               Basis, dBasisdx )

           DO i=1,DIM
             ForceAtIP(i) = SUM( Basis(1:pn) * LocalDisp(i,1:pn) )
           END DO
           
           MomentAtIP(1) = 0.5 * &
               ( SUM( dBasisdx(1:pn,2) * LocalDisp(3,1:pn)) &
               - SUM( dBasisdx(1:pn,3) * LocalDisp(2,1:pn)) )
           MomentAtIp(2) = 0.5 * &
               ( SUM( dBasisdx(1:pn,3) * LocalDisp(1,1:pn)) &
               - SUM( dBasisdx(1:pn,1) * LocalDisp(3,1:pn)) )
           MomentAtIp(3) = 0.5 * &
               ( SUM( dBasisdx(1:pn,1) * LocalDisp(2,1:pn)) &
               - SUM( dBasisdx(1:pn,2) * LocalDisp(1,1:pn)) )
           
           Kmat(Permutation,1:3) = Kmat(Permutation,1:3) + s * ForceAtIp
           Kmat(Permutation,4:6) = Kmat(Permutation,4:6) + s * MomentAtIp
             
           DO i = 1,dim
             IF(ABS(KmatMin(Permutation,i)) > ABS(ForceAtIp(i))) THEN
               KmatMin(Permutation,i) = ForceAtIp(i)
             END IF
             IF(ABS(KmatMin(Permutation,i+3)) > ABS(MomentAtIp(i))) THEN
               KmatMin(Permutation,i+3) = MomentAtIp(i)
             END IF
           END DO
         END DO
       END DO
     END IF



     IF(Permutation == 6) THEN
       KmatFile = ListGetString(SolverParams,'Model Lumping Filename',stat )
       IF(.NOT. stat) KmatFile = "Kmat.dat"

       CALL Info( 'StressSolve', '-----------------------------------------', Level=4 )
       WRITE( Message, * ) 'Saving lumped elastic spring to file ', TRIM(KmatFile)
       CALL Info( 'StressSolve', Message, Level=4 )
       CALL Info( 'StressSolve', '-----------------------------------------', Level=4 )
              
       IF (FixDisplacement) THEN
         Kmat(:,1:3) = Kmat(:,1:3) / dX 
         Kmat(:,4:6) = Kmat(:,4:6) / dFii

         IF( ListGetLogical(SolverParams,'Symmetrisize',stat)) THEN
           Kmat = (Kmat + TRANSPOSE(Kmat)) / 2.0d0
         END IF         
       ELSE
         Kmat = Kmat / Area

         ! Save the Kmatrix prior to inversion to external file
         OPEN (10, FILE= TRIM(KmatFile) // ".inv")
         DO i=1,Permutation
           WRITE(10,'(6ES17.8E3)') Kmat(i,:)
         END DO
         CLOSE(10)              

         OPEN (10, FILE= TRIM(KmatFile) // ".min-inv")
         DO i=1,Permutation
           WRITE(10,'(6ES17.8E3)') KmatMin(i,:)
         END DO
         CLOSE(10)              

         IF(ListGetLogical(SolverParams,'Symmetrisize',stat)) THEN
           Kmat = (Kmat + TRANSPOSE(Kmat)) / 2.0d0
           KmatMin = (KmatMin + TRANSPOSE(KmatMin)) / 2.0d0
         END IF

         CALL InvertMatrix(Kmat,Permutation)
         CALL InvertMatrix(KmatMin,Permutation)

         OPEN (10, FILE= TRIM(KmatFile) // ".min" )
         DO i=1,Permutation
           WRITE(10,'(6ES17.8E3)') KmatMin(i,:)
         END DO
         CLOSE(10)
       END IF

       ! Save the Kmatrix to an external file
       OPEN (10, FILE=KmatFile)
       DO i=1,Permutation
         WRITE(10,'(6ES17.8E3)') Kmat(i,:)
       END DO
       CLOSE(10)

       ! Save the area center to an external file
       OPEN (10, FILE= TRIM(KmatFile) // ".center")
       WRITE(10,'(3ES17.8E3)') Center
       CLOSE(10)
     END IF

     IF(FixDisplacement .AND. Permutation == 6) THEN
       DEALLOCATE( NodalLoads, NodeVisited )
     END IF

   END SUBROUTINE LumpedSprings


!------------------------------------------------------------------------------
! Generalized cartesian lumped mass matrix 
!------------------------------------------------------------------------------
    
    SUBROUTINE LumpedCartesianMass() 
      
      REAL(KIND=dp) :: vol
      TYPE(GaussIntegrationPoints_t) :: IntegStuff      
      INTEGER :: iter, i, j, k, n, t, istat, mat_id, NoEigenModes
      LOGICAL :: GotIt, stat      
      REAL(KIND=dp) :: SqrtMetric,SqrtElementMetric,Amp,Dens
      REAL(KIND=dp) :: Basis(Mesh % MaxElementNodes), &
          dBasisdx(Mesh % MaxElementNodes, 3)
      REAL(KIND=dp) :: x, y, z, U, V, W, S
      REAL (KIND=DP) :: Moment0, Moment1(3), Moment2(3,3), Center(3), MassMatrix(6,6)
      CHARACTER(LEN=MAX_NAME_LEN) :: KmatFile
         
!------------------------------------------------------------------------------
! Do some initialization stuff
!------------------------------------------------------------------------------
      
      vol = 0.0d0
      Moment0 = 0.0d0
      Moment1 = 0.0d0
      Moment2 = 0.0d0
      Center = LumpedCenter

!------------------------------------------------------------------------------
! Integrate the lumped mass over the volume/area
!------------------------------------------------------------------------------
           
100   DO t = 1, Solver % NumberOfActiveElements
        Element => Mesh % Elements( Solver % ActiveElements( t ) )
        Model % CurrentElement => Element
        
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes
        
        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))
        
        body_id = Element % BodyId
        mat_id = ListGetInteger( Model % Bodies( body_id ) % Values, &
            'Material', minv=1,maxv=Model % NumberOfMaterials )
        Material => Model % Materials(mat_id) % Values      
        Density(1:n) = ListGetReal( Material, 'Density', n, NodeIndexes(1:n) )                             
        
        IntegStuff = GaussPoints( Element )
        
        DO i=1,IntegStuff % n
          
          U = IntegStuff % u(i)
          V = IntegStuff % v(i)
          W = IntegStuff % w(i)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element,ElementNodes,U,V,W,&
              SqrtElementMetric,Basis,dBasisdx)
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
          s = SqrtElementMetric * IntegStuff % s(i)          
          x = SUM(ElementNodes % x(1:n) * Basis(1:n)) - Center(1)
          y = SUM(ElementNodes % y(1:n) * Basis(1:n)) - Center(2)
          z = SUM(ElementNodes % z(1:n) * Basis(1:n)) - Center(3)

          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
            s = 2.0 * PI * x * s
          END IF
          vol =  vol + S          
          dens = SUM(Basis(1:n) * Density(1:n) )

          Moment0 = Moment0 + s * dens
          
          Moment1(1) = Moment1(1) + s * x * dens
          Moment1(2) = Moment1(2) + s * y * dens
          Moment1(3) = Moment1(3) + s * z * dens
          
          Moment2(1,1) = Moment2(1,1) + s * ( y*y + z*z)  * dens
          Moment2(2,2) = Moment2(2,2) + s * ( x*x + z*z )  * dens
          Moment2(3,3) = Moment2(3,3) + s * ( x*x + y*y ) * dens
 
          Moment2(1,2) = Moment2(1,2) - s * x * y * dens
          Moment2(1,3) = Moment2(1,3) - s * x * z * dens
          Moment2(2,3) = Moment2(2,3) - s * y * z * dens
        END DO
      END DO

      IF(Vol < AEPS) RETURN

      IF(.FALSE.) THEN
        ! One could also use the center of mass rather than center of force
        Center = Moment1 / Moment0
        GOTO 100
      END IF
      
      Moment2(2,1) = Moment2(1,2)
      Moment2(3,1) = Moment2(1,3)
      Moment2(2,3) = Moment2(3,2)

      CALL ListAddConstReal(Model % Simulation,'res: Mass',Moment0)
      
      CALL ListAddConstReal(Model % Simulation,'res: Lumped Center X',Center(1))
      CALL ListAddConstReal(Model % Simulation,'res: Lumped Center Y',Center(2))
      CALL ListAddConstReal(Model % Simulation,'res: Lumped Center Z',Center(3))
      
      CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia XX',Moment2(1,1))
      CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia YY',Moment2(2,2))
      CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia ZZ',Moment2(3,3))
      CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia XY',Moment2(1,2))
      CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia XZ',Moment2(1,3))
      CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia YZ',Moment2(2,3))
      
      MassMatrix = 0.0d0
      DO i= 1,3
        MassMatrix(i,i) = Moment0
      END DO
      MassMatrix(4:6,4:6) = Moment2

      ! Save the area center to an external file
      KmatFile = ListGetString(SolverParams,'Model Lumping Filename',stat )
      IF(.NOT. stat) KmatFile = "Kmat.dat"
      OPEN (10, FILE= TRIM(KmatFile) // ".mass")
      DO i=1,6
        WRITE(10,'(6ES17.8E3)') MassMatrix(i,:)
      END DO
      CLOSE(10)

    END SUBROUTINE LumpedCartesianMass


  END SUBROUTINE StressSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION StressBoundaryResidual( Model, Edge, Mesh, Quant, Perm, Gnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE StressLocal
     USE DefUtils
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Edge
     REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry
     INTEGER :: i,j,k,n,l,t,dim,DOFs,nd,Pn,En
     LOGICAL :: stat, Found
     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Normal(3), EdgeLength
     REAL(KIND=dp) :: u, v, w, s, detJ

     REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:), dEdgeBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:), ExtPressure(:)
     REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: Force(:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDisplacement(:,:)
     REAL(KIND=dp), ALLOCATABLE :: ElasticModulus(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalPoissonRatio(:)
     REAL(KIND=dp), ALLOCATABLE :: LocalTemp(:), LocalHexp(:,:,:)

     REAL(KIND=dp) :: Residual(3), ResidualNorm, Area
     REAL(KIND=dp) :: ForceSolved(3), Dir(3)
     REAL(KIND=dp) :: Displacement(3)
     REAL(KIND=dp) :: YoungsModulus
     REAL(KIND=dp) :: PoissonRatio
     REAL(KIND=dp) :: Grad(3,3), Strain(3,3), Stress1(3,3), Stress2(3,3)
     REAL(KIND=dp) :: Identity(3,3), YoungsAverage

     LOGICAL :: PlaneStress, Isotropic(2)=.TRUE., CSymmetry = .FALSE.
     TYPE(ValueList_t), POINTER :: Material, Equation, BodyForce, BC
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     SAVE Nodes, EdgeNodes
!------------------------------------------------------------------------------

     ! Initialize:
     ! -----------
     Gnorm = 0.0d0
     Indicator = 0.0d0

     Identity = 0.0d0
     DO i=1,3
        Identity(i,i) = 1.0d0
     END DO

     CSymmetry = CurrentCoordinateSystem() == CylindricSymmetric .OR. &
                 CurrentCoordinateSystem() == AxisSymmetric

     dim = CoordinateSystemDimension()
     DOFs = dim

!    --------------------------------------------------
     Element => Edge % BoundaryInfo % Left

     IF ( .NOT. ASSOCIATED( Element ) ) THEN
        Element => Edge % BoundaryInfo % Right
     ELSE IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) THEN
        Element => Edge % BoundaryInfo % Right
     END IF

     IF ( .NOT. ASSOCIATED( Element ) ) RETURN
     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     En = GetElementNOFNodes( Edge )
     CALL GetElementNodes( EdgeNodes )

     nd = GetElementNOFDOFs( Element )
     Pn = GetElementNOFNodes( Element )
     CALL GetElementNodes( Nodes, UElement=Element )

     ALLOCATE( EdgeBasis(En), dEdgeBasisdx(En,3), x(En), y(En), z(En), &
        ExtPressure(En), Basis(nd), dBasisdx(nd,3), Force(3,En), &
        NodalDisplacement(3,nd), ElasticModulus(6,6,Pn),&
        NodalPoissonRatio(Pn), LocalTemp(nd), LocalHexp(3,3,Pn) )

     LocalTemp = 0
     LocalHexp = 0

     DO l = 1,En
       DO k = 1,Pn
          IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) THEN
             x(l) = Element % TYPE % NodeU(k)
             y(l) = Element % TYPE % NodeV(k)
             z(l) = Element % TYPE % NodeW(k)
             EXIT
          END IF
       END DO
     END DO

     ! Integrate square of residual over boundary element:
     ! ---------------------------------------------------
     Indicator     = 0.0d0
     EdgeLength    = 0.0d0
     YoungsAverage = 0.0d0
     ResidualNorm  = 0.0d0

     BC => GetBC()
     IF ( .NOT.ASSOCIATED( BC ) ) RETURN

     ! Logical parameters:
     ! -------------------
     Equation => GetEquation( Element )
     PlaneStress = GetLogical( Equation, 'Plane Stress' ,Found )

     Material => GetMaterial( Element )
     NodalPoissonRatio(1:pn) = GetReal( &
                  Material, 'Poisson Ratio',Found, Element )
     CALL InputTensor( ElasticModulus, Isotropic(1), &
                 'Youngs Modulus', Material, Pn, Element % NodeIndexes )

     ! Given traction:
     ! ---------------
     Force = 0.0d0
     Force(1,1:En) = GetReal( BC, 'Force 1', Found )
     Force(2,1:En) = GetReal( BC, 'Force 2', Found )
     Force(3,1:En) = GetReal( BC, 'Force 3', Found )

     ! Force in normal direction:
     ! ---------------------------
     ExtPressure(1:En) = GetReal( BC, 'Normal Force', Found )

     ! If dirichlet BC for displacement in any direction given,
     ! nullify force in that directon:
     ! --------------------------------------------------------
     Dir = 1.0d0
     IF ( ListCheckPresent( BC, 'Displacement' ) )   Dir = 0
     IF ( ListCheckPresent( BC, 'Displacement 1' ) ) Dir(1) = 0
     IF ( ListCheckPresent( BC, 'Displacement 2' ) ) Dir(2) = 0
     IF ( ListCheckPresent( BC, 'Displacement 3' ) ) Dir(3) = 0

     ! Elementwise nodal solution:
     ! ---------------------------
     CALL GetVectorLocalSolution( NodalDisplacement, UElement=Element )

     ! Integration:
     ! ------------
     EdgeLength    = 0.0d0
     YoungsAverage = 0.0d0
     ResidualNorm  = 0.0d0

     IntegStuff = GaussPoints( Edge )

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
            EdgeBasis, dEdgeBasisdx )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
           v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
           w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )
   
           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                       Symb, dSymb, u, v, w )

           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        Normal = NormalVector( Edge, EdgeNodes, u, v, .TRUE. )

        u = SUM( EdgeBasis(1:En) * x(1:En) )
        v = SUM( EdgeBasis(1:En) * y(1:En) )
        w = SUM( EdgeBasis(1:En) * z(1:En) )

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
           Basis, dBasisdx )

        ! Stress tensor on the edge:
        ! --------------------------
        CALL LocalStress( Stress1, Strain, NodalPoissonRatio, &
           ElasticModulus, LocalHExp, LocalTemp, &
           Isotropic, CSymmetry, PlaneStress, &
           NodalDisplacement, Basis, dBasisdx, Nodes, dim, pn, nd )

        ! Given force at the integration point:
        ! -------------------------------------
        Residual = MATMUL( Force(:,1:En), EdgeBasis(1:En) ) - &
          SUM( ExtPressure(1:En) * EdgeBasis(1:En) ) * Normal

        ForceSolved = MATMUL( Stress1, Normal )
        Residual = Residual - ForceSolved * Dir

        EdgeLength    = EdgeLength + s
        ResidualNorm  = ResidualNorm  + s * SUM(Residual(1:DIM) ** 2)
        YoungsAverage = YoungsAverage + &
                    s * SUM( ElasticModulus(1,1,1:Pn) * Basis(1:Pn) )
     END DO

     IF ( YoungsAverage > AEPS ) THEN
        YoungsAverage = YoungsAverage / EdgeLength
        Indicator = EdgeLength * ResidualNorm / YoungsAverage
     END IF

     DEALLOCATE( EdgeBasis, dEdgeBasisdx, x, y, z, ExtPressure, Basis, &
      dBasisdx, Force, NodalDisplacement, ElasticModulus, NodalPoissonRatio, &
      LocalTemp, LocalHexp )
!------------------------------------------------------------------------------
   END FUNCTION StressBoundaryResidual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION StressEdgeResidual( Model,Edge,Mesh,Quant,Perm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE StressLocal
     USE DefUtils
     IMPLICIT NONE

     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2)
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Edge
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry

     INTEGER :: i,j,k,l,n,t,dim,DOFs,En,Pn, nd
     LOGICAL :: stat, Found

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Stressi(3,3,2), Jump(3), Identity(3,3)
     REAL(KIND=dp) :: Normal(3)
     REAL(KIND=dp) :: Displacement(3)
     REAL(KIND=dp) :: YoungsModulus
     REAL(KIND=dp) :: PoissonRatio
     REAL(KIND=dp) :: YoungsAverage
     REAL(KIND=dp) :: Grad(3,3), Strain(3,3), Stress1(3,3), Stress2(3,3)

     REAL(KIND=dp), ALLOCATABLE :: LocalTemp(:), LocalHexp(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDisplacement(:,:)
     REAL(KIND=dp), ALLOCATABLE :: ElasticModulus(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalPoissonRatio(:)
     REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:), Basis(:), dBasisdx(:,:)

     LOGICAL :: PlaneStress, Isotropic(2)=.TRUE., CSymmetry

     TYPE(ValueList_t), POINTER :: Material, Equation

     REAL(KIND=dp) :: u, v, w, s, detJ

     REAL(KIND=dp) :: Residual, ResidualNorm, EdgeLength

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     SAVE Nodes, EdgeNodes
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     dim = CoordinateSystemDimension()
     DOFs = dim

     CSymmetry = CurrentCoordinateSystem() == CylindricSymmetric .OR. &
                 CurrentCoordinateSystem() == AxisSymmetric


     Identity = 0.0d0
     Metric   = 0.0d0
     DO i = 1,3
        Metric(i,i)   = 1.0d0
        Identity(i,i) = 1.0d0
     END DO
!
!    ---------------------------------------------
     En = GetElementNOFNodes( Edge )
     CALL GetElementNodes( EdgeNodes, Edge )

     Element => Edge % BoundaryInfo % Left
     pn = GetElementNOFNodes( Element )
     nd = GetElementNOFDOFs( Element )

     Element => Edge % BoundaryInfo % Right
     nd = MAX( nd, GetElementNOFDOFs( Element ) )
     pn = MAX( pn, GetElementNOFNodes( Element ) )

     ALLOCATE( LocalTemp(nd), LocalHexp(3,3,Pn), x(En), y(En), z(En), &
      NodalDisplacement(3,nd), ElasticModulus(6,6,pn), &
      NodalPoissonRatio(pn), EdgeBasis(En), Basis(nd), dBasisdx(nd,3) )

     LocalTemp = 0
     LocalHexp = 0

!    Integrate square of jump over edge:
!    ------------------------------------
     ResidualNorm  = 0.0d0
     EdgeLength    = 0.0d0
     Indicator     = 0.0d0
     Grad          = 0.0d0
     YoungsAverage = 0.0d0

     IntegStuff = GaussPoints( Edge )

     DO t=1,IntegStuff % n

        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
             EdgeBasis, dBasisdx )

        Normal = NormalVector( Edge, EdgeNodes, u, v, .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
           v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
           w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                       Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        Stressi = 0.0d0
        DO i = 1,2
           IF ( i==1 ) THEN
              Element => Edge % BoundaryInfo % Left
           ELSE
              Element => Edge % BoundaryInfo % Right
           END IF

           IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) CYCLE

           pn = GetElementNOFNodes( Element )
           nd = GetElementNOFDOFs( Element )
           CALL GetElementNodes( Nodes, Element )
           DO j = 1,en
              DO k = 1,pn
                 IF ( Edge % NodeIndexes(j) == Element % NodeIndexes(k) ) THEN
                    x(j) = Element % TYPE % NodeU(k)
                    y(j) = Element % TYPE % NodeV(k)
                    z(j) = Element % TYPE % NodeW(k)
                    EXIT
                 END IF
              END DO
           END DO

           u = SUM( EdgeBasis(1:En) * x(1:En) )
           v = SUM( EdgeBasis(1:En) * y(1:En) )
           w = SUM( EdgeBasis(1:En) * z(1:En) )

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
               Basis, dBasisdx )

           ! Logical parameters:
           ! -------------------
           Equation => GetEquation( Element )
           PlaneStress = GetLogical( Equation,'Plane Stress',Found )

           ! Material parameters:
           ! --------------------
           Material => GetMaterial( Element )
           NodalPoissonRatio(1:pn) = GetReal( Material, 'Poisson Ratio', Found, Element )
           CALL InputTensor( ElasticModulus, Isotropic(1), &
                         'Youngs Modulus', Material, pn, Element % NodeIndexes )

           ! Elementwise nodal solution:
           ! ---------------------------
           CALL GetVectorLocalSolution( NodalDisplacement, UElement=Element )

           ! Stress tensor on the edge:
           ! --------------------------
           CALL LocalStress( Stress1, Strain, NodalPoissonRatio, &
              ElasticModulus, LocalHExp, LocalTemp, Isotropic, CSymmetry, PlaneStress, &
              NodalDisplacement, Basis, dBasisdx, Nodes, dim, pn, nd )

           Stressi(:,:,i) = Stress1
        END DO

        EdgeLength  = EdgeLength + s
        Jump = MATMUL( ( Stressi(:,:,1) - Stressi(:,:,2)), Normal )
        ResidualNorm = ResidualNorm + s * SUM( Jump(1:DIM) ** 2 )

        YoungsAverage = YoungsAverage + s *  &
                    SUM( ElasticModulus(1,1,1:pn) * Basis(1:pn) )
     END DO

     YoungsAverage = YoungsAverage / EdgeLength
     Indicator = EdgeLength * ResidualNorm / YoungsAverage

     DEALLOCATE( LocalTemp, LocalHexp, x, y, z, NodalDisplacement, &
       ElasticModulus, NodalPoissonRatio, EdgeBasis, Basis, dBasisdx )
!------------------------------------------------------------------------------
   END FUNCTION StressEdgeResidual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION StressInsideResidual( Model, Element,  &
                      Mesh, Quant, Perm, Fnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE StressLocal
     USE DefUtils
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Element
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes

     INTEGER :: i,j,k,l,m,n,nd,t,dim,DOFs,I1(6),I2(6)
     INTEGER, ALLOCATABLE :: Indexes(:)

     LOGICAL :: stat, Found

     TYPE( Variable_t ), POINTER :: Var

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp) :: Density
     REAL(KIND=dp) :: PoissonRatio
     REAL(KIND=dp) :: Damping
     REAL(KIND=dp) :: Displacement(3),Identity(3,3), YoungsAverage
     REAL(KIND=dp) :: Grad(3,3), Strain(3,3), Stress1(3,3), Stress2(3,3)
     REAL(KIND=dp) :: Energy

     REAL(KIND=dp), ALLOCATABLE :: ElasticModulus(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDensity(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalPoissonRatio(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDamping(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDisplacement(:,:)
     REAL(KIND=dp), ALLOCATABLE :: LocalHexp(:,:,:), vec(:)
     REAL(KIND=dp), ALLOCATABLE :: Stressi(:,:,:), LocalTemp(:)
     REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalForce(:,:), Veloc(:,:), Accel(:,:)

     LOGICAL :: PlaneStress, CSymmetry, Isotropic(2)=.TRUE., Transient

     REAL(KIND=dp) :: u, v, w, s, detJ
     REAL(KIND=dp) :: Residual(3), ResidualNorm, Area

     TYPE(ValueList_t), POINTER :: Material, BodyForce, Equation

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     SAVE Nodes
!------------------------------------------------------------------------------
     ! Initialize:
     ! -----------
     Fnorm     = 0.0d0
     Indicator = 0.0d0

     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     Metric = 0.0d0
     DO i=1,3
        Metric(i,i) = 1.0d0
     END DO

     dim = CoordinateSystemDimension()
     DOFs = dim 

     CSymmetry = CurrentCoordinateSystem() == CylindricSymmetric .OR. &
                 CurrentCoordinateSystem() == AxisSymmetric

     ! Element nodal points:
     ! ---------------------
     nd = GetElementNOFDOFs()
     n  = GetElementNOFNodes()
     CALL GetElementNodes( Nodes )

     ALLOCATE( ElasticModulus(6,6,nd), NodalDensity(n), NodalPoissonRatio(n), &
         NodalDamping(n), NodalDisplacement(3,nd), LocalHExp(3,3,n), vec(nd), &
         Stressi(3,3,nd), LocalTemp(nd), Basis(nd), dBasisdx(3,nd), &
         NodalForce(4,n), Veloc(3,nd), Accel(3,nd) )

     LocalTemp = 0
     LocalHexp = 0

     ! Logical parameters:
     ! -------------------
     equation => GetEquation()
     PlaneStress = GetLogical( Equation, 'Plane Stress',Found )

     ! Material parameters:
     ! --------------------
     Material => GetMaterial()

     CALL InputTensor( ElasticModulus, Isotropic(1), &
           'Youngs Modulus', Material, n, Element % NodeIndexes )

     NodalPoissonRatio(1:n) = GetReal( Material, 'Poisson Ratio', Found )

     ! Check for time dep.
     ! -------------------
     IF ( ListGetString( Model % Simulation, 'Simulation Type') == 'transient' ) THEN
        Transient = .TRUE.
        Var => VariableGet( Model % Variables, 'Displacement', .TRUE. )

        nd = GetElementDOFs( Indexes )

        Veloc = 0.0d0
        Accel = 0.0d0
        DO i=1,DOFs
           Veloc(i,1:nd) = Var % PrevValues(DOFs*(Var % Perm(Indexes(1:nd))-1)+i,1)
           Accel(i,1:nd) = Var % PrevValues(DOFs*(Var % Perm(Indexes(1:nd))-1)+i,2)
        END DO
        NodalDensity(1:n) = GetReal( Material, 'Density', Found )
        NodalDamping(1:n) = GetReal( Material, 'Damping', Found )
     ELSE
        Transient = .FALSE.
     END IF

     ! Elementwise nodal solution:
     ! ---------------------------
     CALL GetVectorLocalSolution( NodalDisplacement )

     ! Body Forces:
     ! ------------
     BodyForce => GetBodyForce()

     NodalForce = 0.0d0

     IF ( ASSOCIATED( BodyForce ) ) THEN
        NodalForce(1,1:n) = NodalForce(1,1:n) + GetReal( &
            BodyForce, 'Stress BodyForce 1', Found )
        NodalForce(2,1:n) = NodalForce(1,1:n) + GetReal( &
            BodyForce, 'Stress BodyForce 2', Found )
        NodalForce(3,1:n) = NodalForce(1,1:n) + GetReal( &
            BodyForce, 'Stress BodyForce 3', Found )
     END IF

     Identity = 0.0D0
     DO i = 1,dim
        Identity(i,i) = 1.0D0
     END DO
     CSymmetry = .FALSE.

     Var => VariableGet( Model % Variables, 'Stress 1' )
     IF ( ASSOCIATED( Var ) ) THEN

       ! If stress already computed:
       ! ---------------------------
       I1(1:6) = (/ 1,2,3,1,2,1 /)
       I2(1:6) = (/ 1,2,3,2,3,3 /)
       DO i=1,6
         CALL GetScalarLocalSolution(Vec(1:nd),'Stress ' // CHAR(i+ICHAR('0')))
         Stressi(I1(i),I2(i),1:nd) = Vec(1:nd)
         Stressi(I2(i),I1(i),1:nd) = Vec(1:nd)
       END DO
     ELSE
       ! Values of the stress tensor at node points:
       ! -------------------------------------------
       DO i = 1,n
         u = Element % TYPE % NodeU(i)
         v = Element % TYPE % NodeV(i)
         w = Element % TYPE % NodeW(i)

         stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
             Basis, dBasisdx )

         CALL LocalStress( Stressi(:,:,i), Strain, NodalPoissonRatio, &
                   ElasticModulus, LocalHExp, LocalTemp, Isotropic, CSymmetry, PlaneStress, &
                   NodalDisplacement, Basis, dBasisdx, Nodes, dim, n, nd )
       END DO
     END IF

     ! Integrate square of residual over element:
     ! ------------------------------------------
     ResidualNorm = 0.0d0
     Fnorm = 0.0d0
     Area = 0.0d0
     Energy = 0.0d0
     YoungsAverage = 0.0d0

     IntegStuff = GaussPoints( Element )

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( Basis(1:n) * Nodes % x(1:n) )
           v = SUM( Basis(1:n) * Nodes % y(1:n) )
           w = SUM( Basis(1:n) * Nodes % z(1:n) )

           CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,u,v,w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        ! Residual of the diff.equation:
        ! ------------------------------
        Residual = 0.0d0
        DO i = 1,3
           Residual(i) = -SUM( NodalForce(i,1:n) * Basis(1:n) )

           IF ( Transient ) THEN
              Residual(i) = Residual(i) + SUM(NodalDensity(1:n)*Basis(1:n)) * &
                            SUM( Accel(i,1:nd) * Basis(1:nd) )
              Residual(i) = Residual(i) + SUM(NodalDamping(1:n)*Basis(1:n)) * &
                            SUM( Veloc(i,1:nd) * Basis(1:nd) )
           END IF

           DO j = 1,3
             Residual(i) = Residual(i) - SUM(Stressi(i,j,1:nd)*dBasisdx(1:nd,j))
           END DO
        END DO

!       IF ( CSymmetry ) THEN
!          DO k=1,3
!             Residual(1) = Residual(1) + ...
!          END DO
!       END IF

       ! Dual norm of the load:
       ! ----------------------
        DO i = 1,dim
           Fnorm = Fnorm + s * SUM( NodalForce(i,1:n) * Basis(1:n) ) ** 2
        END DO

        YoungsAverage = YoungsAverage + s*SUM( ElasticModulus(1,1,1:n) * Basis(1:n) )

        ! Energy:
        ! -------
        CALL LocalStress( Stress1, Strain, NodalPoissonRatio, &
           ElasticModulus, LocalHExp, LocalTemp, Isotropic, CSymmetry, PlaneStress, &
           NodalDisplacement, Basis, dBasisdx, Nodes, dim, n, nd )

        Energy = Energy + s*DDOTPROD(Strain,Stress1,dim) / 2.0d0

        Area = Area + s
        ResidualNorm = ResidualNorm + s * SUM( Residual(1:dim) ** 2 )
     END DO

     YoungsAverage = YoungsAverage / Area
     Fnorm = Energy
     Indicator = Area * ResidualNorm / YoungsAverage
 
     DEALLOCATE( ElasticModulus, NodalDensity, NodalPoissonRatio,  &
         NodalDamping, NodalDisplacement, LocalHExp, vec, Stressi, &
         LocalTemp, Basis, dBasisdx, NodalForce, Veloc, Accel )


CONTAINS

!------------------------------------------------------------------------------
  FUNCTION DDOTPROD(A,B,N) RESULT(C)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION :: A(:,:),B(:,:),C
    INTEGER :: N
!------------------------------------------------------------------------------
    INTEGER :: I,J
!------------------------------------------------------------------------------
    C = 0.0D0
    DO i = 1,N
       DO j = 1,N
          C = C + A(i,j)*B(i,j)
       END DO
    END DO
!------------------------------------------------------------------------------
  END FUNCTION DDOTPROD
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   END FUNCTION StressInsideResidual
!------------------------------------------------------------------------------
