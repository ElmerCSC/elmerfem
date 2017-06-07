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
! *  Module for solving losses by utilizing the Fourier expansion of 
! *  degrees of freedom.
! *
! *  Authors: Juha Ruokolainen, Peter Råback, Mika Malinen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 19.4.2013
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{


!------------------------------------------------------------------------------
!> Initialization for the primary solver: FourierLoss
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE FourierLossSolver_init0( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------

  USE DefUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver  
  TYPE(Model_t) :: Model    
  REAL(KIND=dp) :: dt       
  LOGICAL :: Transient      
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams, DGSolverParams
  LOGICAL :: Found
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  INTEGER :: i,j,k,n,m,mysolver,soln
  TYPE(Solver_t), POINTER :: Solvers(:), PSolver
  INTEGER, POINTER :: Active(:)
  LOGICAL :: ElementalField
  LOGICAL :: OldKeywordStyle, SeparateComponents, NodalLosses
  INTEGER :: Ncomp

  SolverParams => GetSolverParams()

  VarName = ListGetString( SolverParams,'Target Variable',Found )
  IF(.NOT. Found ) THEN
    CALL Fatal('FourierLossSolver','Give > Target Variable < !')
  END IF


  OldKeywordStyle = ListCheckPresentAnyMaterial( Model,'Harmonic Loss Linear Coefficient')
  IF( OldKeywordStyle ) THEN
    NComp = 2
  ELSE
    Ncomp = 0
    DO i=1,10
      IF( ListCheckPresentAnyMaterial( Model,'Harmonic Loss Coefficient '//TRIM(I2S(i)) ) ) THEN
        Ncomp = i
      ELSE
        EXIT
      END IF
    END DO
  END IF

  IF( OldKeywordStyle ) THEN
    CALL ListAddNewString( SolverParams,'Variable','Fourier Loss Linear' )
    CALL ListAddNewString( SolverParams,'Exported Variable 1','Fourier Loss Quadratic' )
  ELSE
    SeparateComponents = ListGetLogical( SolverParams,'Separate Loss Components',Found )
    IF( SeparateComponents ) THEN
      CALL ListAddNewString( SolverParams,'Variable','Fourier Loss 1' )
      DO i=2, NComp
        CALL ListAddNewString( SolverParams,'Exported Variable '//TRIM(I2S(i-1)),&
            'Fourier Loss '//TRIM(I2S(i)) )
      END DO
    ELSE
      CALL ListAddNewString( SolverParams,'Variable','Fourier Loss' )
    END IF
  END IF


  NodalLosses = ListGetLogical( SolverParams,'Calculate Nodal Losses',Found )
  IF( NodalLosses ) THEN
    CALL ListAddString( SolverParams, &
        NextFreeKeyword('Exported Variable', SolverParams),'Nodal Fourier Loss')
  END IF

  
  ! We are really solving using DG so no need for funny initialization
  IF (GetLogical(SolverParams,'Discontinuous Galerkin',Found)) RETURN

  ! If elemental fields are not requested don't compute them
  ElementalField = GetLogical( SolverParams, 'Calculate Elemental Fields', Found)
  IF(.NOT. ElementalField) RETURN

  ! Ok, Dirty way of adding DG field for postprocessing
  ! Add a new solver with exactly the same active locations as the current one.
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

  ! Create DG solver structures on-the-fly without actually solving the matrix equations. 
  ALLOCATE(Solvers(n+1))
  Solvers(1:n) = Model % Solvers

  Solvers(n+1) % Values => ListAllocate()
  DGSolverParams => Solvers(n+1) % Values
  CALL ListAddLogical( DGSolverParams, 'Discontinuous Galerkin', .TRUE. )
  Solvers(n+1) % DG = .TRUE.
  Solvers(n+1) % PROCEDURE = 0
  Solvers(n+1) % ActiveElements => NULL()
  DEALLOCATE(Model % Solvers)
  Model % Solvers => Solvers
  Model % NumberOfSolvers = n+1

  CALL ListAddString( DGSolverParams, 'Exec Solver', 'never' )
  CALL ListAddLogical( DGSolverParams, 'No Matrix',.TRUE.)
  CALL ListAddLogical( DGSolverParams, 'Optimize Bandwidth',.FALSE.)
  CALL ListAddString( DGSolverParams, 'Equation', 'FourierLoss_Dummy' )
  CALL ListAddString( DGSolverParams, 'Procedure', &
      'FourierLoss FourierLossSolver_Dummy',.FALSE. )
  

  IF( OldKeywordStyle ) THEN
    CALL ListAddString( DGSolverParams, 'Variable', 'Fourier Loss Linear e' )
    CALL ListAddString( DGSolverParams, 'Exported Variable 1','Fourier Loss Quadratic e')
  ELSE
    IF( SeparateComponents ) THEN
      CALL ListAddString( DGSolverParams,'Variable','Fourier Loss 1e' )
      DO i=2, NComp
        CALL ListAddString( DGSolverParams,'Exported Variable '//TRIM(I2S(i-1)),&
            'Fourier Loss '//TRIM(I2S(i))//'e' )
      END DO
    ELSE
      CALL ListAddString( DGSolverParams,'Variable','Fourier Loss e' )
    END IF
  END IF

END SUBROUTINE FourierLossSolver_init0



!------------------------------------------------------------------------------
SUBROUTINE FourierLossSolver_Dummy(Model,Solver,dt,Transient)
  !------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !------------------------------------------------------------------------------
END SUBROUTINE FourierLossSolver_Dummy
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Make a Fourier series expansion of time-varying coefficients in the finite
!> element expansion of a target solution field and use the resulting spatial
!> amplitude fields to estimate losses when a component-wise loss coefficient
!> is given.
!------------------------------------------------------------------------------
SUBROUTINE FourierLossSolver( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------

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
  TYPE VarPointer_t
    TYPE(Variable_t), POINTER :: Var
  END TYPE VarPointer_t
  
  TYPE(ValueList_t),POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, FourierName, VarString
  INTEGER :: dim, Nsize, FourierDofs, i, nlen, mlen, icomp
  LOGICAL :: Found
  INTEGER :: TimesVisited = 0
  INTEGER, POINTER :: FourierPerm(:), EPerm(:)
  REAL(KIND=dp) :: Norm, Omega
  REAL(KIND=dp), POINTER :: FourierField(:)
  REAL(KIND=dp) :: at0,at1,at2,at3
#ifndef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: CPUTime,RealTime
#endif
  REAL(KIND=dp), ALLOCATABLE :: BodyLoss(:,:), SeriesLoss(:,:), CompLoss(:)
  TYPE(Variable_t), POINTER :: TargetVar, LossVar, NodalLossVar
  REAL(KIND=dp), POINTER :: TargetField(:), PrevTargetField(:,:), SaveRhs(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: OtherRhs(:,:)
  REAL(KIND=dp) :: TotalLoss
  LOGICAL :: EndCycle, FourierOutput, SimpsonsRule, ExactIntegration, &
      ElementalField, AvField, DirectField
  TYPE(Solver_t), POINTER :: TargetSolverPtr
  LOGICAL :: OldKeywordStyle, SeparateComponents, SumComponents, NodalLosses
  INTEGER :: Ncomp, NVar 

  TYPE(VarPointer_t), ALLOCATABLE :: CompVars(:), CompVarsE(:), FourierVars(:)


  
  SAVE TimesVisited, CompLoss, CompVars, FourierVars, CompVarsE, SeriesLoss, BodyLoss
  !-------------------------------------------------------------------------------

  CALL Info( 'FourierLossSolver', '-------------------------------------',Level=4 )
  CALL Info( 'FourierLossSolver', 'Computing Fourier losses             ',Level=4 )
  CALL Info( 'FourierLossSolver', '-------------------------------------',Level=4 )

  TimesVisited = TimesVisited + 1
  dim = CoordinateSystemDimension()
  SolverParams => GetSolverParams()
  LossVar => Solver % Variable
  
  at0 = RealTime()



  Ncomp = 0
  OldKeywordStyle = ListCheckPresentAnyMaterial( Model,'Harmonic Loss Linear Coefficient')
  IF( OldKeywordStyle ) THEN
    CALL Info('FourierLossSolver','Using old keyword style',Level=5)
    NComp = 2
  ELSE
    Ncomp = 0
    DO i=1,10
      IF( ListCheckPresentAnyMaterial( Model,'Harmonic Loss Coefficient '//TRIM(I2S(i)) ) ) THEN
        Ncomp = i
      ELSE
        EXIT
      END IF
    END DO
  END IF
  
  IF( Ncomp == 0 ) THEN
    CALL Fatal('FourierLossSolver','Some material must have > Harmonic Loss Coefficient i <')
  END IF
  CALL Info('FourierLossSolver','Considering number of components: '//TRIM(I2S(Ncomp)),Level=5)
  
  IF( OldKeywordStyle ) THEN
    SeparateComponents = .TRUE.
  ELSE
    SeparateComponents = ListGetLogical( SolverParams,'Separate Loss Components',Found )
  END IF 
  SumComponents = .NOT. SeparateComponents 

  
  IF( SumComponents ) THEN
    Nvar = 1
  ELSE   
    NVar = Ncomp
  END IF

  IF(.NOT. ALLOCATED( CompLoss) ) ALLOCATE( CompLoss(Ncomp) )
  CompLoss = 0.0_dp
  
  
  ! Get the fields for the nodal results
  !------------------------------------------------------------------  
  IF(.NOT. ALLOCATED( CompVars ) ) ALLOCATE( CompVars(Nvar) )   
  CompVars(1) % Var => Solver % Variable
  IF( OldKeywordStyle ) THEN
    CompVars(2) % Var => VariableGet( Solver % Mesh % Variables,&
        'Fourier Loss Quadratic' )
  ELSE
    DO icomp = 2, Nvar
      CompVars(icomp) % Var => VariableGet( Solver % Mesh % Variables,&
          'Fourier Loss '//TRIM(I2S(icomp)) )
    END DO
  END IF

  DO icomp = 1, NVar
    IF(.NOT. ASSOCIATED( CompVars(icomp) % Var ) )THEN
      CALL Fatal('FourierLossSolver','Variable for component does not exist: '//TRIM(I2S(icomp)))
    END IF
  END DO
    

  

  ! Check for Elemental (Discontinuous Galerkin) Field 
  !------------------------------------------------------------------
  ElementalField = GetLogical( SolverParams, 'Calculate Elemental Fields', Found)
  IF( ElementalField ) THEN
    IF(.NOT. ALLOCATED( CompVarsE ) ) ALLOCATE( CompVarsE(Nvar) )
    
    IF( OldKeywordStyle ) THEN
      CompVarsE(1) % Var => VariableGet( Solver % Mesh % Variables,&
          'Fourier Loss Linear e' )
      CompVarsE(2) % Var => VariableGet( Solver % Mesh % Variables,&
          'Fourier Loss Quadratic e' )
    ELSE
      IF( Nvar == 1 ) THEN
          CompVarsE(1) % Var => VariableGet( Solver % Mesh % Variables,&
              'Fourier Loss e' )                
      ELSE
        DO icomp = 1, Nvar
          CompVarsE(icomp) % Var => VariableGet( Solver % Mesh % Variables,&
              'Fourier Loss '//TRIM(I2S(icomp))//'e' )        
        END DO
      END IF
    END IF

    DO icomp = 1, NVar
      IF(.NOT. ASSOCIATED( CompVarsE(icomp) % Var ) )THEN
        CALL Fatal('FourierLossSolver','Variable e for component does not exist: '//TRIM(I2S(icomp)))
      END IF
    END DO
  END IF


  ! Check whether we want to compute nodal loss directly (in terms of J per node)
  !------------------------------------------------------------------------------
  NodalLossVar => VariableGet( Solver % Mesh % Variables,'Nodal Fourier Loss' )
  NodalLosses = ASSOCIATED( NodalLossVar ) 
  

  ! Check that we have something to compute
  !------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) THEN
    CALL Fatal('FourierLossSolver','Solver has no matrix associated')
  END IF
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN  

  ! Fetch the target field, e.g. magnetic vector potential
  !--------------------------------------------------------------------------
  VarName = GetString(SolverParams,'Target Variable',Found)
  TargetVar => VariableGet( Solver % Mesh % Variables, VarName ) 
  IF( .NOT. ASSOCIATED( TargetVar ) ) THEN
    CALL Fatal('FourierLossSolver','Target field not present: '//TRIM(VarName) )
  END IF

  IF( .NOT. ASSOCIATED( TargetVar % PrevValues ) ) THEN
    CALL Fatal('FourierLossSolver','Target field does not have PrevValues: '//TRIM(VarName) )    
  END IF
  TargetField => TargetVar % Values
  PrevTargetField => TargetVar % PrevValues(:,:)
  Nsize = SIZE( TargetField )


  ! The target field is an AV solution 
  DirectField = ListGetLogical( SolverParams,'Target Variable Direct',Found)
  IF( DirectField ) THEN
    CALL Info('FourierLossSolver','Using the target field directly!')
    AvField = .FALSE.
  ELSE
    ! The target field is an AV solution 
    AvField = ListGetLogical( SolverParams,'Target Variable AV',Found)
    IF( .NOT. Found ) THEN
      AvField = ( SIZE( TargetVar % Perm ) > Solver % Mesh % NumberOfNodes )
    END IF
    IF( AvField ) THEN
      IF( TargetVar % Dofs > 1 ) THEN
        CALL Fatal('FourierLossSolver','Assuming only one component for AV field!')
      END IF
    ELSE
      IF( dim == 3 ) THEN
        IF( TargetVar % Dofs /= 3 ) THEN
          CALL Fatal('FourierLossSolver','Assuming precisely three component in 3D!')
        END IF
      ELSE
        IF( TargetVar % Dofs /= 1 ) THEN
          CALL Fatal('FourierLossSolver','Assuming only one component in 2D!')
        END IF
      END IF
    END IF
  END IF
    

  !--------------------------------------------------------------------
  ! Generate also a pointer to the target variable solver
  !-------------------------------------------------------------------
  mlen = LEN_TRIM(VarName)
  NULLIFY(TargetSolverPtr) 
  DO i=1,Model % NumberOfSolvers 
    VarString = ListGetString( Model % Solvers(i) % Values, 'Variable', Found )
    IF (Found) THEN
      nlen = LEN_TRIM(VarString)
      IF ( VarString(1:nlen) == VarName(1:mlen) ) THEN
        TargetSolverPtr => Model % Solvers(i)
        EXIT
      END IF
    END IF
  END DO
  IF( .NOT. ASSOCIATED( TargetSolverPtr ) ) THEN
    CALL Fatal('FourierLossSolver','Target field solver cannot be found: '//TRIM(VarName) )
  END IF

  FourierDofs = 1 + 2 * ListGetInteger( SolverParams,'Fourier Series Components',Found)
  IF(.NOT. Found ) THEN
    CALL Fatal('FourierLossSolver','Give > Fourier Series Components < !')
  END IF
  
  ! Fetch the Fourier field that includes the components of the transform
  ! If it has not been created before, create it now
  !------------------------------------------------------------------------
  IF( .NOT. ALLOCATED( FourierVars ) ) THEN
    ! Allocate the components: 1st one is related to constant part
    ! then even ones to cosine and odd ones to sine. 
    ALLOCATE( FourierVars( FourierDofs ) )
    FourierOutput = ListGetLogical( SolverParams,'Fourier Series Output',Found) 

    DO i=1, FourierDofs     
      NULLIFY( FourierField ) 
      ALLOCATE( FourierField( Nsize ) )
      FourierField = 0.0_dp
      IF( i == 1 ) THEN
        FourierName = TRIM( VarName )//' Fourier 0'
      ELSE IF( MODULO(i,2) == 0 ) THEN
        FourierName = TRIM(VarName)//' Fourier Cos'//TRIM(I2S(i/2))
      ELSE
        FourierName = TRIM(VarName)//' Fourier Sin'//TRIM(I2S(i/2))
      END IF
        
      CALL VariableAddVector( Solver % Mesh % Variables, Solver % Mesh, Solver, &
          FourierName, TargetVar % Dofs, FourierField, TargetVar % Perm, FourierOutput )
      FourierVars(i) % Var => VariableGet( Solver % Mesh % Variables, FourierName)
    END DO
  END IF
  FourierPerm => TargetVar % Perm 

  IF(.NOT. ALLOCATED( SeriesLoss) ) THEN
    ALLOCATE( SeriesLoss(Ncomp,FourierDofs), &
        BodyLoss(Ncomp,Model % NumberOfBodies) )
  END IF
    
  ! Fetch frequency that is used both for Fourier transform and loss computation
  !-----------------------------------------------------------------------------
  Omega = GetAngularFrequency(Found=Found)
  IF( .NOT. Found ) THEN  
    CALL Fatal('FourierLosses','The Fourier transform requires frequency!')
  END IF
  WRITE( Message,'(A,ES12.3)') 'Base frequency for Fourier transform: ',Omega / (2*PI)
  CALL Info('FourierLossSolver', Message, Level=8 )


  ! Check whether the exact integration of Fourier coefficients is used
  ! Always use it by default. 
  !-----------------------------------------------------------------------------
  ExactIntegration = .NOT. GetLogical(SolverParams, 'Inexact Integration', Found)
  IF (ExactIntegration) THEN
    CALL Info('FourierLossSolver','Using exact integration')    
  ELSE
    SimpsonsRule = GetLogical( SolverParams,'Simpsons Rule',Found )
    IF( SimpsonsRule ) THEN
      CALL Info('FourierLossSolver','Using Simpsons rule for integration')
    ELSE
      CALL Info('FourierLossSolver','Using trapetsoidal rule for integration')    
    END IF
  END IF

  ! Add the loss so that it would appear from the start if saving data by
  ! SaveScalars to some external file. 
  !-----------------------------------------------------------------------
  IF( TimesVisited == 1 ) THEN
    IF( OldKeywordStyle ) THEN
      CALL ListAddConstReal( Model % Simulation,'res: fourier loss',0.0_dp )
      CALL ListAddConstReal( Model % Simulation,'res: fourier loss quadratic',0.0_dp )
    ELSE
      CALL ListAddConstReal( Model % Simulation,'res: fourier loss total',0.0_dp )
      DO icomp=1, Ncomp
        CALL ListAddConstReal( Model % Simulation,'res: fourier loss '//TRIM(I2S(icomp)),0.0_dp )       
      END DO
    END IF
  END IF

  ! Perform discrete fourier transform (DFT)
  ! If the cycle has come to end, then continue to calculate the losses.
  !-----------------------------------------------------------------------------
  EndCycle = .FALSE.
  CALL LocalFourierTransform( EndCycle )
  at1 = RealTime()

  WRITE( Message,'(A,ES12.3)') 'Fourier transform time: ',at1-at0
  CALL Info('FourierLossSolver', Message, Level=5 )
  IF( .NOT. EndCycle ) RETURN

  CALL DefaultInitialize()
  IF( NVar > 1 ) THEN
    CALL Info('FourierLossSolver','Allocating multiple r.h.s. vectors',Level=12)
    IF(.NOT. ALLOCATED( OtherRhs ) ) THEN
      ALLOCATE( OtherRhs ( SIZE( Solver % Matrix % Rhs ), NVar - 1) )
    END IF
    OtherRhs = 0.0_dp
    SaveRhs => Solver % Matrix % Rhs   
    SaveRhs = 0.0_dp
  END IF


  CALL BulkAssembly( Ncomp )

  
  ! These three are probably not needed for this solver of limited features
  ! They are left here for possible future needs.
  ! CALL DefaultFinishBulkAssembly()
  ! CALL DefaultFinishAssembly()
  ! CALL DefaultDirichletBCs()      
  at2 = RealTime()
  WRITE( Message,'(A,ES12.3)') 'Assembly time: ',at2-at1
  CALL Info( 'FourierLossSolver', Message, Level=5 )


  !------------------------------------------------------------------------------     
  IF( SeparateComponents ) THEN
    ! Solver other components first, so that the values are saved to Solver % Var % Values
    DO icomp=2, Ncomp
      Solver % Matrix % Rhs => OtherRhs(:,icomp-1)
      Norm = DefaultSolve()
      CompVars(icomp) % Var % Values = LossVar % Values
    END DO
    Solver % Matrix % Rhs => SaveRhs    
    IF(Nvar > 1) DEALLOCATE( OtherRhs ) 
  END IF
  
  Norm = DefaultSolve()
  
  at3 = RealTime()
  WRITE( Message,'(A,ES12.3)') 'Solution time: ',at3-at2
  CALL Info( 'FourierLossSolver', Message, Level=5 )

  ! Print and save the lumped results
  !-------------------------------------------------------------------------------
  CALL CommunicateLosess()
  
  ! We need to perform the DFT also for the missing part
  !----------------------------------------------------------------------
  EndCycle = .TRUE.
  CALL LocalFourierTransform( EndCycle ) 
  

CONTAINS 

  !-------------------------------------------------------------------
  ! Perform discrete Fourier transform of scalar coefficients in the
  ! finite element expansion (degrees of freedom)
  !----------------------------------------------------------------
  SUBROUTINE LocalFourierTransform( EndCycle )

    LOGICAL :: EndCycle

    INTEGER :: n0,i,j,k,cycles
    REAL(KIND=dp) :: t0,c0,c1,time,ratio,ta,tb,tcycle,A0,A1
    REAL(KIND=dp), POINTER :: Component(:)
    LOGICAL :: LeftRule, RightRule, InitPending = .FALSE., Cosine
    REAL(KIND=dp) :: AngleA, AngleB, AngleC, CompCoeffA, CompCoeffB, &
        ValA, ValB, dtFraction
    REAL(KIND=dp), PARAMETER :: RatioEps = 1.0d-8
    INTEGER :: CurrentCycle, PreviousCycle=0

    SAVE :: time, ratio, tcycle, InitPending, t0, PreviousCycle  

    ! Often the Fourier transform ends just at the end of time step and then the 
    ! initialization of the new Fourier series is left to the next timestep. This
    ! enables that also the Fourier series may possibly be saved at the end. 
    ! when choosing saving so that it coinsides with the end of Fourier cycle.
    IF( InitPending ) THEN
      !  If this is not initialized to zero then additional tricks for normalization are needed.      
      DO i=1,FourierDofs 
        FourierVars(i) % Var % Values = 0.0_dp
      END DO
    END IF
    InitPending = .FALSE.

    IF( EndCycle ) THEN
      ratio = 1.0 - ratio
      IF( ratio < RatioEps ) THEN
        InitPending = .TRUE.
        RETURN
      END IF
      DO i=1,FourierDofs 
        FourierVars(i) % Var % Values = 0.0_dp
      END DO
      LeftRule = .TRUE. 
      CALL Info('FourierLossSolver','Finishing the timestep using left rule',Level=8)
    ELSE      
      time = GetTime()

      ! Figure out the time to start Fourier analysis
      ! ----------------------------------------------
      n0 = GetInteger(SolverParams,'Fourier Start Timestep', Found )
      IF ( Found ) THEN
        i = GetTimestep() 
        IF( i < n0 ) THEN
          CALL Info('FourierLossSolver','Fourier transform not yet active, returning...',Level=6)
          RETURN
        ELSE IF( i == n0 ) THEN
          t0 = time
        END IF
      ELSE
        cycles = GetInteger(SolverParams,'Fourier Start Cycles', Found )
        IF ( Found ) THEN
          t0 = cycles * 2 * PI / Omega
        ELSE
          t0 = GetCReal(SolverParams,'Fourier Start Time', Found )
          IF( .NOT. Found ) THEN
            IF( TimesVisited == 1 ) t0 = 0.0_dp
          END IF
        END IF
      END IF

      ! Substract the start time so that we start conveniently from zero
      time = time - t0

      ! Compute the fraction of the timestep needed 
      LeftRule = .FALSE.
      IF( time <= 0.0 ) THEN
        ratio = 0.0_dp
        CALL Info('FourierLossSolver','Fourier transform not yet active, returning...',Level=6)
        RETURN
      ELSE IF( time - dt < 0.0 ) THEN
        ratio = time / dt
        LeftRule = .TRUE.
        CALL Info('FourierLossSolver','Starting Fourier transform cycle',Level=6)
      ELSE
        ratio = 1.0_dp
      END IF

      ! Reset mean after elapsed time
      ! Default is to relax after one complete cycle
      cycles = GetInteger(SolverParams,'Fourier Integrate Cycles', Found )
      IF(.NOT. Found ) cycles = 1
      tcycle = cycles * 2 * PI / Omega

      RightRule = .FALSE.
      CurrentCycle = NINT( time / tcycle )

      IF( CurrentCycle /= PreviousCycle ) THEN
        ! Check whether we have proceeded to a new cycle
        ! 1) if the previous and current timestep point to different timestep
        RightRule = ( INT(time/tcycle) /= INT((time-dt)/tcycle) ) 
        
        ! 2) if the current timestep is close enough of the end of cycle
        RightRule = RightRule .OR. &
            ( ABS( CurrentCycle - time/tcycle ) < RatioEps )

        IF( RightRule ) THEN
          IF( LeftRule ) THEN
            CALL Fatal('FourierLossSolver','Cannot use left and right rule at the same time!')
          END IF
          ratio = 1.0_dp - (time-INT((time)/tcycle)*tcycle)/dt
          
          CALL Info('FourierLossSolver','Finising Fourier transform cycle',Level=6)
          
          ! Return a True flag so that we know that we might need to start also the cycle.
          EndCycle = .TRUE.        
          PreviousCycle = CurrentCycle
        END IF
      END IF
    END IF

    ! Set the time interval for integration, [ta,tb]
    ta = time - dt 
    tb = time

    IF( LeftRule ) THEN
      ta = time - ratio * dt
    ELSE IF( RightRule ) THEN
      tb = time - (1-ratio) * dt 
    END IF
    dtFraction = (tb-ta) / tcycle

    WRITE( Message,'(A,ES12.3,A,ES12.3,A)') 'Timestep interval: [',ta,',',tb,']'
    CALL Info('FourierLossSolver',Message,Level=8)

    ! Now perform discrete Fourier transform
    ! Some sin & cos coefficients are computed just once
    ! in order to make nodewise summation as economical as possible.
    DO j=1,FourierDofs 
      Component => FourierVars(j) % Var % Values

      IF (ExactIntegration) THEN
        !----------------------------------------------------------------
        ! Here the Fourier coefficients are integrated exactly by using 
        ! analytical expressions
        !----------------------------------------------------------------
        IF( j == 1 ) THEN
          ! Coefficients to express the constant part
          CompCoeffA = 0.5d0 * dtFraction 
          CompCoeffB = 0.5d0 * dtFraction
        ELSE
          ! Coefficients to express Fourier expansion in terms of cos(kwt) and sin(kwt) 
          k = j/2
          c1 = k * Omega
          c0 = c1 * ta
          AngleA = Omega * ta
          AngleB = Omega * tb       

          Cosine = ( MODULO(j,2) == 0 )             
          IF( Cosine ) THEN
            A0 = 2.0d0/(c1 * tcycle) * ( SIN( k * AngleB ) - SIN( k * AngleA ) )
            A1 = 2.0d0/(c1 * tcycle) * ( 1.0d0/(c1 * (tb-ta) ) * ( COS( k * AngleB ) - &
                COS( k * AngleA ) ) + SIN( k * AngleB ) )
          ELSE
            A0 = 2.0d0/(c1 * tcycle) * ( -COS( k * AngleB ) + COS( k * AngleA ) )
            A1 = 2.0d0/(c1 * tcycle) * ( 1.0d0/(c1 * (tb-ta) ) * ( SIN( k * AngleB ) - &
                SIN( k * AngleA ) ) - COS( k * AngleB ) )
          END IF

          CompCoeffA = A0 - A1 
          CompCoeffB = A1
        END IF
        !---------------------------------------------------------------------------------
        ! If the solver dt /= tb - ta, then we apply
        ! linear interpolation to obtain the end value and modify coefficients to express
        ! the Fourier coefficients accordingly
        !--------------------------------------------------------------------------------
        IF( j == 1 ) THEN
          IF( LeftRule ) THEN
            CompCoeffA = ratio * CompCoeffA
            CompCoeffB = ( 1 - ratio ) * CompCoeffA + CompCoeffB
          ELSE IF( RightRule ) THEN
            CompCoeffB = ratio * CompCoeffB
            CompCoeffA = ( 1 - ratio ) * CompCoeffB + CompCoeffA
          END IF
        ELSE
          IF( LeftRule ) THEN
            CompCoeffA = ratio * (A0 - A1) 
            CompCoeffB = ( 1 - ratio ) * A0 + ratio * A1
          ELSE IF ( RightRule ) THEN
            CompCoeffA = A0 - ratio * A1
            CompCoeffB = ratio * A1
          END IF
        END IF

      ELSE
        !---------------------------------------------------------------------------
        ! Here numerical integration is done as 
        ! \int f(t) dt = (tb-ta) * (f(tb)+f(ta)) / 2
        ! or alternatively with Simpson's rule
        ! \int f(t) dt = (tb-ta) * (t(tb)+4*f(tc)+f(ta)) / 6
        ! The Simpson's rule might not pay off since the approximation of 
        ! fields is still linear.
        ! Note that the Fourier series is multiplied by 2 because <cos^2>=<sin^2>=1/2.
        !----------------------------------------------------------------------------
        IF( j == 1 ) THEN
          ! Coefficient for constant part 
          CompCoeffA = 0.5d0 * dtFraction 
          CompCoeffB = 0.5d0 * dtFraction
        ELSE
          ! Coefficient for cos(kf) or sin(kf) 
          k = j/2
          AngleA = Omega * ta
          AngleB = Omega * tb       
          Cosine = ( MODULO(j,2) == 0 ) 

          IF( Cosine ) THEN
            CompCoeffA = COS( k * AngleA ) * dtFraction
            CompCoeffB = COS( k * AngleB ) * dtFraction
          ELSE
            CompCoeffA = SIN( k * AngleA ) * dtFraction
            CompCoeffB = SIN( k * AngleB ) * dtFraction
          END IF

          ! In Simpsons rule the coefficient may still be applied at the ends as
          ! the linear approximation enables it. 
          IF( SimpsonsRule ) THEN
            AngleC = ( AngleA + AngleB ) / 2.0_dp
            IF( Cosine ) THEN
              CompCoeffA = ( CompCoeffA + 2 * COS( k * AngleC ) * dtFraction ) / 3.0
              CompCoeffB = ( CompCoeffB + 2 * COS( k * AngleC ) * dtFraction ) / 3.0
            ELSE
              CompCoeffA = ( CompCoeffA + 2 * SIN( k * AngleC ) * dtFraction ) / 3.0
              CompCoeffB = ( CompCoeffB + 2 * SIN( k * AngleC ) * dtFraction ) / 3.0
            END IF
          END IF
        END IF

        ! If we are the left or right end of the integration interval then 
        ! use linear interpolation of the end value.
        IF( LeftRule ) THEN
          CompCoeffA = ratio * CompCoeffA
          CompCoeffB = ( 1 - ratio ) * CompCoeffA + CompCoeffB
        ELSE IF( RightRule ) THEN
          CompCoeffB = ratio * CompCoeffB
          CompCoeffA = ( 1 - ratio ) * CompCoeffB + CompCoeffA
        END IF

      END IF

      ! And finally the summation for the discrete Fourier transform
      ! Note that here the fields are vectors.
      Component = Component + CompCoeffA * PrevTargetField(:,1) + CompCoeffB * TargetField

    END DO

  END SUBROUTINE LocalFourierTransform


  !------------------------------------------------------------------------------
  ! Assembly the mass matrix and r.h.s. for computing the losses using the 
  ! Galerkin method.
  !------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly(ncomp)
    !------------------------------------------------------------------------------
    INTEGER :: ncomp
    
    INTEGER :: elem,i,j,k,p,q,n,nd,nt,t,BodyId,nsum
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element, TargetElement
    REAL(KIND=dp) :: weight,coeff,coeff2,detJ,GradAtIp(3,3),CurlAtIP(3),ValAtIp
    LOGICAL :: Found, Found2
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: Freq, FourierFreq, LossCoeff(ncomp), FreqPower(ncomp), FieldPower(ncomp)
    REAL(KIND=dp), POINTER :: Component(:), WrkArray(:,:)
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:), dgForce(:)
    REAL(KIND=dp), ALLOCATABLE :: ElemField(:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:), RotWBasis(:,:)
    INTEGER, ALLOCATABLE :: Indeces(:), Pivot(:)
    LOGICAL :: DG
    
    SAVE Nodes
    
    SaveRhs => Solver % Matrix % Rhs

    n = MAX(Solver % Mesh % MaxElementDOFs,Solver % Mesh % MaxElementNodes)
    ALLOCATE( STIFF(n,n), FORCE(ncomp,n), dgForce(n), Pivot(n), Basis(n), dBasisdx(n,3) )
    
    IF (AVField) THEN
      nt = MAX(TargetSolverPtr % Mesh % MaxElementDOFs,TargetSolverPtr % Mesh % MaxElementNodes)
      ALLOCATE( WBasis(nt,3), RotWBasis(nt,3), Indeces(nt), ElemField(nt) )
    ELSE
      ALLOCATE( ElemField(n), Indeces(n) )
    END IF

      
    Freq = Omega / (2*PI)

    DG = GetLogical(SolverParams,'Discontinuous Galerkin',Found )

    
    IF( OldKeywordStyle ) THEN
      FreqPower(1) = GetCReal( SolverParams,'Harmonic Loss Linear Frequency Exponent',Found )
      IF(.NOT. Found ) FreqPower(1) = 1.0_dp
     
      FreqPower(2) = GetCReal( SolverParams,'Harmonic Loss Quadratic Frequency Exponent',Found )
      IF(.NOT. Found ) FreqPower(2) = 2.0_dp

      FieldPower(1) = GetCReal( SolverParams,'Harmonic Loss Linear Exponent',Found ) 
      IF(.NOT. Found ) FieldPower(1) = 2.0_dp

      FieldPower(2) = GetCReal( SolverParams,'Harmonic Loss Quadratic Exponent',Found ) 
      IF( .NOT. Found ) FieldPower(2) = 2.0_dp      
    ELSE
      WrkArray => ListGetConstRealArray( SolverParams,'Harmonic Loss Frequency Exponent',Found )
      IF( Found ) THEN 
        IF( SIZE( WrkArray,1 ) < Ncomp ) THEN
          CALL Fatal('FourierLossSolver','> Harmonic Loss Frequency Exponent < too small')
        END IF
        FreqPower(1:Ncomp) = WrkArray(1:Ncomp,1)
      ELSE       
        DO icomp = 1, Ncomp
          FreqPower(icomp) = GetCReal( SolverParams,'Harmonic Loss Frequency Exponent '//TRIM(I2S(icomp)) )
        END DO
      END IF
        
      WrkArray => ListGetConstRealArray( SolverParams,'Harmonic Loss Field Exponent',Found )
      IF( Found ) THEN
        IF( SIZE( WrkArray,1 ) < Ncomp ) THEN
          CALL Fatal('FourierLossSolver','> Harmonic Loss Field Exponent < too small')
        END IF
        FieldPower(1:Ncomp) = WrkArray(1:Ncomp,1)        
      ELSE
        DO icomp = 1, Ncomp
          FieldPower(icomp) = GetCReal( SolverParams,'Harmonic Loss Field Exponent '//TRIM(I2S(icomp)) )
        END DO
      END IF
    END IF

    ! Sum over the loss for each frequency, each body, and for the combined effect
    SeriesLoss = 0.0_dp
    BodyLoss = 0.0_dp
    CompLoss = 0.0_dp

    ! Assemble the matrix equation for computing the losses
    !------------------------------------------------------------------
    DO elem = 1,GetNOFActive()

      ! Element information
      ! ---------------------
      Element => GetActiveElement( elem )
      CALL GetElementNodes( Nodes )
      ! nd = GetElementNOFDOFs()
      nd = GetElementDOFs( Indeces )
      n  = GetElementNOFNodes()

      IF (AVField) THEN
        ! We need basis function data from the target variable mesh...
        TargetElement => GetActiveElement( elem, TargetSolverPtr )
        IF ( .NOT. ASSOCIATED(TargetElement) ) THEN
          CALL Fatal('FourierLossSolver','Element on target mesh cannot be associated')
        END IF
        nt = GetElementDOFs( Indeces, TargetElement, TargetSolverPtr )
      END IF

      ! Integrate over local element:
      ! -----------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      BodyId = Element % BodyId
      Material => GetMaterial()

      
      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
            IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )

        !---------------------------------------------------------------------
        ! Get edge basis functions if needed. Given that only linear edge
        ! interpolation functions are available currently, the following should 
        ! work. If higher-order versions become available, we then need to add 
        ! an argument which specifies the approximation degree.  
        !---------------------------------------------------------------------
        IF (AVField) THEN
          CALL GetEdgeBasis(Element, WBasis, RotWBasis, Basis, dBasisdx) 
        END IF

        Weight = IntegStuff % s(t) * detJ

        ! As there are so many components it does not really pay of much to save
        ! the matrix. The r.h.s. will always be more laboursome to assembly.
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
          END DO
        END DO

        ! Compute the Fourier losses at the given Gaussian integration point
        DO j=1,FourierDofs
          ! Frequency should be same for sine and cosine series
          ! This eats any functional dependence with a dummy parameter
          IF( j == 1 ) THEN
            FourierFreq = 0.0_dp
          ELSE
            k = j/2 
            FourierFreq = k * Freq
          END IF

          ! currently the loss coefficient may depend only of frequency
          ! and on the material section
          IF( OldKeywordStyle ) THEN
            LossCoeff(1) = ListGetFun( Material,'Harmonic Loss Linear Coefficient',FourierFreq, Found)   
            LossCoeff(2) = ListGetFun( Material,'Harmonic Loss Quadratic Coefficient',FourierFreq, Found2 )   
            IF(.NOT. (Found .OR. Found2 ) ) CYCLE            
          ELSE
            Found2 = .FALSE.
            DO icomp = 1, Ncomp
              LossCoeff(icomp) = ListGetFun( Material,&
                  'Harmonic Loss Coefficient '//TRIM(I2S(icomp)),FourierFreq, Found)      
              IF( Found ) Found2 = .TRUE.
            END DO
            IF(.NOT. Found2 ) CYCLE
          END IF

          
          ! For even j we have cosine series, for odd sine series
          Component => FourierVars(j) % Var % Values

          IF( DirectField ) THEN
            ElemField(1:nd) = Component( FourierPerm( Indeces(1:nd) ) )
            ValAtIp = SUM( Basis(1:nd) * ElemField(1:nd) )
          ELSE IF ( AVField ) THEN
            ElemField(1:nt) = Component( FourierPerm( Indeces(1:nt) ) )

            ! Here we assume that the solution contains both a nodal interpolation
            ! field (scalar) and edge interpolation part (vector). The gradient is 
            ! applied to the scalar while the curl operator is applied to the vector

            CurlAtIp(1) = SUM( ElemField(n+1:nt) * RotWBasis(1:(nt-n),1) )
            CurlAtIp(2) = SUM( ElemField(n+1:nt) * RotWBasis(1:(nt-n),2) )
            IF ( dim > 2 ) &
                CurlAtIp(3) = SUM( ElemField(n+1:nt) * RotWBasis(1:(nt-n),3) )             
          ELSE IF( dim == 3 ) THEN
            ElemField(1:nd) = Component( 3 * (FourierPerm( Indeces(1:nd))-1) + 1 )
            DO k=1,3
              GradAtIp(1,k) =  SUM( ElemField(1:nd) * dBasisdx(1:nd,k) )
            END DO
            ElemField(1:nd) = Component( 3 * (FourierPerm( Indeces(1:nd))-1) + 2 )
            DO k=1,3
              GradAtIp(2,k) =  SUM( ElemField(1:nd) * dBasisdx(1:nd,k) )
            END DO
            ElemField(1:nd) = Component( 3 * (FourierPerm( Indeces(1:nd))-1) + 3 )
            DO k=1,3
              GradAtIp(3,k) =  SUM( ElemField(1:nd) * dBasisdx(1:nd,k) )
            END DO

            CurlAtIp(1) = GradAtIp(3,2) - GradAtIp(2,3)
            CurlAtIp(2) = GradAtIp(1,3) - GradAtIp(3,1)
            CurlAtIp(3) = GradAtIp(2,1) - GradAtIp(1,2)
          ELSE
            ElemField(1:nd) = Component( FourierPerm( Indeces(1:nd) ) )

            CurlAtIp(1) =  SUM( ElemField(1:nd) * dBasisdx(1:nd,2) )
            CurlAtIp(2) = -SUM( ElemField(1:nd) * dBasisdx(1:nd,1) )
          END IF

          IF(.NOT. DirectField ) THEN
            ValAtIp = SQRT( SUM( CurlAtIP(1:dim) ** 2 ) )
          END IF

          
          ! Losses by components
          DO icomp = 1, Ncomp 
            Coeff = Weight * LossCoeff(icomp) * ( FourierFreq ** FreqPower(icomp) ) &
                * ( ValAtIp ** FieldPower(icomp) )
            IF( SumComponents ) THEN
              nsum = 1
            ELSE
              nsum = icomp
            END IF

            FORCE(nsum,1:nd) = FORCE(nsum, 1:nd) + Coeff * Basis(1:nd)          
            SeriesLoss(icomp,j) = SeriesLoss(icomp,j) + Coeff 
            BodyLoss(icomp,BodyId) = BodyLoss(icomp,BodyId) + Coeff
          END DO

        END DO

      END DO
      
      ! Update global matrices from local matrices 
      !------------------------------------------------------------------------------
      IF( SumComponents ) THEN
        nsum = 1
      ELSE
        nsum = ncomp
      END IF

      
      DO icomp = 1, Nsum
        IF( icomp == 1 ) THEN
          CALL DefaultUpdateEquations( STIFF, FORCE(icomp,1:nd) )
        ELSE
          Solver % Matrix % Rhs => OtherRhs(:,icomp-1)
          CALL DefaultUpdateForce( FORCE(icomp,1:nd) )
          Solver % Matrix % Rhs => SaveRhs
        END IF
      END DO
        
      ! After this the STIFF and FORCE are corrupted 
      IF( ElementalField ) THEN
        EPerm => CompVarsE(1) % Var % Perm
        DO icomp = 1, Nsum
          CALL LUdecomp(STIFF,n,pivot)
          CALL LUSolve(n,STIFF,FORCE(icomp,:),pivot)
          CompVarsE(icomp) % Var % Values(EPerm(Element % DGIndexes(1:n))) = FORCE(icomp,1:n)
        END DO
      END IF
  
    END DO

    
    ! Assembly of the face terms when using DG:
    !------------------------------------------------------------------
    IF ( DG ) THEN
      IF (GetLogical(SolverParams,'Average Within Materials',Found)) THEN
        dgFORCE = 0.0d0
        CALL AddLocalFaceTerms( STIFF, dgFORCE )
      END IF
    END IF

    
    ! If we want to save the nodal losses directly do that now
    !--------------------------------------------------------------------
    IF( NodalLosses ) THEN
      NodalLossVar % Values = Solver % Matrix % rhs
      DO icomp = 2, Nsum
        NodalLossVar % Values = NodalLossVar % Values + OtherRhs(:,icomp-1)
      END DO
    END IF

    DEALLOCATE( STIFF, FORCE, dgForce, Pivot, Basis, dBasisdx, Indeces, ElemField )
    IF (AVField) THEN
      DEALLOCATE( WBasis, RotWBasis )
    END IF

    
    
  END SUBROUTINE BulkAssembly


  !------------------------------------------------------------------------------
  !> Parallel reduce, print and save the losses into a file if requested
  !------------------------------------------------------------------------------
  SUBROUTINE CommunicateLosess()

    INTEGER :: i,j,k
    CHARACTER(LEN=MAX_NAME_LEN) :: LossesFile
    
    DO k=1,Ncomp
      DO j=1,FourierDofs
        SeriesLoss(k,j) = ParallelReduction(SeriesLoss(k,j)) 
      END DO

      DO j=1,Model % NumberOfBodies
        BodyLoss(k,j) = ParallelReduction(BodyLoss(k,j))
      END DO

      CompLoss(k) = SUM( SeriesLoss(k,:) )
    END DO
    TotalLoss = SUM( CompLoss )

    
    IF( OldKeywordStyle ) THEN
      CALL ListAddConstReal( Model % Simulation,'res: fourier loss linear',CompLoss(1) )
      CALL ListAddConstReal( Model % Simulation,'res: fourier loss quadratic',CompLoss(2) )
    ELSE
      CALL ListAddConstReal( Model % Simulation,'res: fourier loss total',TotalLoss )
      DO k=1,Ncomp
        CALL ListAddConstReal( Model % Simulation,'res: fourier loss '//TRIM(I2S(k)),CompLoss(k) )
      END DO
    END IF
      
    ! Output of the losses on screen    
    ! First the component losses, the the body losses
    DO k=1,NComp
      IF( OldKeywordStyle ) THEN
        IF( k == 1 ) THEN
          CALL Info('FourierLossSolver','Fourier loss linear by components',Level=6)
        ELSE
          CALL Info('FourierLossSolver','Fourier loss quadratic by components',Level=6)
        END IF
      ELSE
        CALL Info('FourierLossSolver','Wavewise Fourier loss for component: '//TRIM(I2S(k)),Level=6)
      END IF
      
      DO j=1,FourierDofs 
        IF( j == 1 ) THEN
          WRITE( Message,'(A,ES12.3)') 'CONST : ',SeriesLoss(k,j)
          CALL Info('FourierLossSolver', Message, Level=6 )
        ELSE
          i = j/2
          IF( MODULO( j,2 ) == 0 ) THEN
            WRITE( Message,'(A,I0,A,ES12.3)') 'COS_',i,' : ',SeriesLoss(k,j)
            CALL Info('FourierLossSolver', Message, Level=6 )
          ELSE
            WRITE( Message,'(A,I0,A,ES12.3)') 'SIN_',i,' : ',SeriesLoss(k,j)
            CALL Info('FourierLossSolver', Message, Level=6 )
          END IF
        END IF
      END DO
    
      IF( OldKeywordStyle ) THEN
        IF( k == 1 ) THEN
          CALL Info('FourierLossSolver','Fourier loss linear by bodies',Level=6)
        ELSE
          CALL Info('FourierLossSolver','Fourier loss quadratic by bodies',Level=6)
        END IF
      ELSE
        CALL Info('FourierLossSolver','Bodywise Fourier loss for component: '//TRIM(I2S(k)),Level=6)
      END IF

      DO j=1,Model % NumberOfBodies
        IF( BodyLoss(k,j) < TINY( CompLoss(k) ) ) CYCLE
        WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyLoss(k,j)
        CALL Info('FourierLossSolver', Message, Level=6 )
      END DO

      WRITE( Message,'(A,ES12.3)') 'Total component loss: ',CompLoss(k)
      CALL Info('FourierLossSolver',Message, Level=5 )
    END DO
    
    WRITE( Message,'(A,ES12.3)') 'Total loss: ',TotalLoss
    CALL Info('FourierLossSolver',Message, Level=5 )
    
    IF( Parenv % MyPe == 0 ) THEN
      LossesFile = ListGetString(SolverParams,'Fourier Loss Filename',Found )
      IF( Found ) THEN
        OPEN (10, FILE=LossesFile)
        WRITE( 10,'(A)')  '!body_id   loss(1)   loss(2) ....'
        DO j=1,Model % NumberOfBodies
          WRITE( 10,* ) j, BodyLoss(1:Ncomp,j)
        END DO
        CALL Info('FourierLossSolver', &
            'Fourier losses for bodies was saved to file: '//TRIM(LossesFile),Level=6 )
        CLOSE(10)
      END IF
    END IF

    ! For debugging
    IF( .FALSE. ) THEN
      DO i=1,FourierDofs
        PRINT *,'fourier range:',i,&
            MINVAL(FourierVars(i) % Var % Values ), &
            MAXVAL(FourierVars(i) % Var % Values )
      END DO
      
      ! For debugging
      DO i=1,NComp
        PRINT *,'loss component range:',i,&
            MINVAL(CompVars(i) % Var % Values ), &
            MAXVAL(CompVars(i) % Var % Values )
      END DO
    END IF
      
    
    !------------------------------------------------------------------------------
  END SUBROUTINE CommunicateLosess
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
END SUBROUTINE FourierLossSolver
!------------------------------------------------------------------------------

!> \}
