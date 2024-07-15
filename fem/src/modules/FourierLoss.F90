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
SUBROUTINE FourierLossSolver_init( Model,Solver,dt,Transient )
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
  TYPE(ValueList_t),POINTER :: SolverParams
  LOGICAL :: Found, Found2
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  INTEGER :: i,j,k,n,m,mysolver,soln
  LOGICAL :: ElementalField, NodalField
  LOGICAL :: OldKeywordStyle, SeparateComponents, NodalLosses
  CHARACTER(LEN=MAX_NAME_LEN) :: Pref
  INTEGER :: Ncomp
  CHARACTER(*), PARAMETER :: Caller = 'FourierLossSolver'

  
  SolverParams => GetSolverParams()

  OldKeywordStyle = ListCheckPresentAnyMaterial( Model,'Harmonic Loss Linear Coefficient') .OR. &
      ListCheckPresentAnyMaterial( Model,'Harmonic Loss Quadratic Coefficient') 
  IF( OldKeywordStyle ) THEN
    IF(ParEnv % MyPe == 0 ) THEN
      PRINT *,'You are using old keywords'
      PRINT *,'Obsolite: "Harmonic Loss Linear Coefficient"'
      PRINT *,'Obsolite: "Harmonic Loss Quadratic Coefficient"'
      PRINT *,'Obsolite: "Harmonic Loss Linear Frequency Exponent"'
      PRINT *,'Obsolite: "Harmonic Loss Quadratic Frequency Exponent"'
      PRINT *,'Obsolite: "Harmonic Loss Linear Exponent"'
      PRINT *,'Obsolite: "Harmonic Loss Quadratic Exponent"'
      PRINT *,'etc. see the documentation.'
    END IF
    CALL Fatal(Caller,'The keyword format is obsolite. Make the fixes and continue!')
  END IF
  
  Pref = ListGetString( SolverParams,'Scalars Prefix',Found )
  IF(.NOT. Found) Pref = 'res:'
  
  VarName = ListGetString( SolverParams,'Target Variable',Found )
  IF(.NOT. Found ) THEN
    CALL Fatal(Caller,'Give > Target Variable < !')
  END IF
      
  Ncomp = 0
  DO i=1,10
    IF( ListCheckPresentAnyMaterial( Model,'Harmonic Loss Coefficient '//I2S(i) ) ) THEN
      Ncomp = i
    ELSE
      EXIT
    END IF
  END DO

  ! Add the loss so that it would appear from the start if saving data by
  ! SaveScalars to some external file etc. 
  !-----------------------------------------------------------------------
  CALL ListAddConstReal( Model % Simulation,TRIM(Pref)//' fourier loss total',0.0_dp )
  DO i=1, Ncomp
    CALL ListAddConstReal( Model % Simulation,TRIM(Pref)//' fourier loss '//I2S(i),0.0_dp )       
  END DO

  ! This is usually not needed so just allocate a dummy. 
  CALL ListAddNewString( SolverParams,'Variable','-nooutput FourierTmp' )
  
  NodalField = ListGetLogical( SolverParams, 'Calculate Nodal Fields', Found)
  ElementalField = ListGetLogical( SolverParams, 'Calculate Elemental Fields', Found2 )
  IF(.NOT. (Found .OR. Found2) ) THEN
    ElementalField = .TRUE.
    CALL ListAddLogical( SolverParams,'Calculate Elemental Fields',ElementalField )
    CALL Info(Caller,'Computing elemental fields only by default when neither given!')
  END IF
  
  SeparateComponents = ListGetLogical( SolverParams,'Separate Loss Components',Found )

  IF( NodalField ) THEN
    IF( SeparateComponents ) THEN
      CALL ListAddString( SolverParams,NextFreeKeyword('Exported Variable ',SolverParams), &
          'Fourier Loss[Fourier Loss:'//I2S(Ncomp)//']' )
    ELSE
      CALL ListAddString( SolverParams,NextFreeKeyword('Exported Variable ',SolverParams), &
          'Fourier Loss' )
    END IF
  END IF

  IF( ElementalField ) THEN
    IF( SeparateComponents ) THEN
      CALL ListAddString( SolverParams,NextFreeKeyword('Exported Variable ',SolverParams), &
          '-dg Fourier Loss e[Fourier Loss e:'//I2S(Ncomp)//']' )
    ELSE
      CALL ListAddString( SolverParams,NextFreeKeyword('Exported Variable ',SolverParams), &
          '-dg Fourier Loss e' )
    END IF
  END IF

  NodalLosses = ListGetLogical( SolverParams,'Calculate Nodal Losses',Found )
  IF( NodalLosses ) THEN
    CALL ListAddString( SolverParams, &
        NextFreeKeyword('Exported Variable', SolverParams),'Nodal Fourier Loss')
  END IF    
  
END SUBROUTINE FourierLossSolver_init



!------------------------------------------------------------------------------
!> Make a Fourier series expansion of time-varying coefficients in the finite
!> element expansion of a target solution field and use the resulting spatial
!> amplitude fields to estimate losses when a component-wise loss coefficient
!> is given.
!------------------------------------------------------------------------------
SUBROUTINE FourierLossSolver( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------

  USE DefUtils
  USE MagnetoDynamicsUtils
  
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
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, FourierName, VarString, Pref
  INTEGER :: dim, Nsize, FourierDofs, i, nlen, mlen, icomp
  LOGICAL :: Found, FlipActive
  INTEGER :: TimesVisited = 0
  INTEGER, POINTER :: FourierPerm(:), EPerm(:)
  REAL(KIND=dp) :: Norm, Omega
  REAL(KIND=dp), POINTER :: FourierField(:)
  REAL(KIND=dp) :: at0,at1,at2,at3
  REAL(KIND=dp), ALLOCATABLE :: BodyLoss(:,:), SeriesLoss(:,:), CompLoss(:)
  TYPE(Variable_t), POINTER :: TargetVar, LossVar, NodalLossVar
  REAL(KIND=dp), POINTER :: TargetField(:), PrevTargetField(:,:)
  REAL (KIND=dp), POINTER CONTIG :: SaveRhs(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: OtherRhs(:,:)
  REAL(KIND=dp) :: TotalLoss
  LOGICAL :: EndCycle, FourierOutput, SimpsonsRule, ExactIntegration, &
      ElementalField, NodalField, AvField, DirectField
  TYPE(Solver_t), POINTER :: TargetSolverPtr
  LOGICAL :: SeparateComponents, SumComponents, NodalLosses
  INTEGER :: Ncomp, NVar, tdofs
  TYPE(VarPointer_t), ALLOCATABLE :: CompVars(:), CompVarsE(:), FourierVars(:)
  CHARACTER(*), PARAMETER :: Caller = 'FourierLossSolver'

  
  SAVE TimesVisited, CompLoss, CompVars, FourierVars, CompVarsE, SeriesLoss, BodyLoss
  !-------------------------------------------------------------------------------

  CALL Info( Caller, '-------------------------------------',Level=4 )
  CALL Info( Caller, 'Computing Fourier losses             ',Level=4 )
  CALL Info( Caller, '-------------------------------------',Level=4 )

  TimesVisited = TimesVisited + 1
  dim = CoordinateSystemDimension()
  SolverParams => GetSolverParams()
  LossVar => Solver % Variable
  
  at0 = RealTime()
  
  Pref = ListGetString( SolverParams,'Scalars Prefix',Found )
  IF(.NOT. Found) Pref = 'res:'
  
  Ncomp = 0
  DO i=1,10
    IF( ListCheckPresentAnyMaterial( Model,'Harmonic Loss Coefficient '//I2S(i) ) ) THEN
      Ncomp = i
    ELSE
      EXIT
    END IF
  END DO  
  IF( Ncomp == 0 ) THEN
    CALL Fatal(Caller,'Some material must have > Harmonic Loss Coefficient i <')
  END IF
  CALL Info(Caller,'Considering number of components: '//I2S(Ncomp),Level=5)
  
  SeparateComponents = ListGetLogical( SolverParams,'Separate Loss Components',Found )
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
  NodalField = GetLogical( SolverParams, 'Calculate Nodal Fields', Found)
  IF( NodalField ) THEN
    IF(.NOT. ALLOCATED( CompVars ) ) ALLOCATE( CompVars(Nvar) )   
    IF( Nvar == 1 ) THEN
      CompVars(1) % Var => VariableGet( Solver % Mesh % Variables,&
          'Fourier Loss' )                
    ELSE
      DO icomp = 1, Nvar
        CompVars(icomp) % Var => VariableGet( Solver % Mesh % Variables,&
            'Fourier Loss '//I2S(icomp) )        
      END DO
    END IF
    DO icomp = 1, NVar
      IF(.NOT. ASSOCIATED( CompVars(icomp) % Var ) )THEN
        CALL Fatal(Caller,'Variable for component does not exist: '//I2S(icomp))
      END IF
    END DO
  END IF
  
  ! Check for Elemental (Discontinuous Galerkin) Field 
  !------------------------------------------------------------------
  ElementalField = GetLogical( SolverParams, 'Calculate Elemental Fields', Found)
  IF( ElementalField ) THEN
    IF(.NOT. ALLOCATED( CompVarsE ) ) ALLOCATE( CompVarsE(Nvar) )    
    IF( Nvar == 1 ) THEN
      CompVarsE(1) % Var => VariableGet( Solver % Mesh % Variables,&
          'Fourier Loss e' )                
    ELSE
      DO icomp = 1, Nvar
        CompVarsE(icomp) % Var => VariableGet( Solver % Mesh % Variables,&
            'Fourier Loss e '//I2S(icomp) )        
      END DO
    END IF
    DO icomp = 1, NVar
      IF(.NOT. ASSOCIATED( CompVarsE(icomp) % Var ) )THEN
        CALL Fatal(Caller,'Variable e for component does not exist: '//I2S(icomp))
      END IF
    END DO
  END IF

  ! Check whether we want to compute nodal loss directly (in terms of J per node)
  !------------------------------------------------------------------------------
  NodalLossVar => VariableGet( Solver % Mesh % Variables,'Nodal Fourier Loss' )
  NodalLosses = ASSOCIATED( NodalLossVar ) 
  

  ! Check that we have something to compute
  !------------------------------------------------------------------------------
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN  

  ! Fetch the target field, e.g. magnetic vector potential
  !--------------------------------------------------------------------------
  VarName = GetString(SolverParams,'Target Variable',Found)
  TargetVar => VariableGet( Solver % Mesh % Variables, VarName ) 
  IF( .NOT. ASSOCIATED( TargetVar ) ) THEN
    CALL Fatal(Caller,'Target field not present: '//TRIM(VarName) )
  END IF

  IF( .NOT. ASSOCIATED( TargetVar % PrevValues ) ) THEN
    CALL Fatal(Caller,'Target field does not have PrevValues: '//TRIM(VarName) )    
  END IF
  TargetField => TargetVar % Values
  PrevTargetField => TargetVar % PrevValues(:,:)
  Nsize = SIZE( TargetField )
  tdofs = TargetVar % Dofs
    
  FlipActive = ( TargetVar % PeriodicFlipActive )
  IF(FlipActive) THEN
    CALL Info(Caller,'Assuming initial field to have conforming flips')
  END IF
  
  ! The target field is an AV solution 
  AvField = .FALSE.
  DirectField = ListGetLogical( SolverParams,'Target Variable Direct',Found)
  IF( DirectField ) THEN
    CALL Info(Caller,'Using the target field with '//I2S(tdofs)//' dofs directly!')
  ELSE
    IF( dim == 3 ) THEN
      ! Check whether the target field is an AV solution 
      AvField = ListGetLogical( SolverParams,'Target Variable AV',Found)
      IF( .NOT. Found ) THEN
        AvField = ( SIZE( TargetVar % Perm ) > Solver % Mesh % NumberOfNodes )
      END IF
      IF( AvField ) THEN
        IF( tdofs > 1 ) THEN
          CALL Fatal(Caller,'Assuming only one component for AV field!')
        END IF
      ELSE
        IF( Tdofs /= 3 ) THEN
          CALL Fatal(Caller,'Assuming precisely three nodal components in 3D!')
        END IF
      END IF
    ELSE
      IF( Tdofs /= 1 ) THEN
        CALL Fatal(Caller,'Assuming only one nodal component Az in 2D!')
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
    CALL Fatal(Caller,'Target field solver cannot be found: '//TRIM(VarName) )
  END IF

  FourierDofs = 1 + 2 * ListGetInteger( SolverParams,'Fourier Series Components',Found)
  IF(.NOT. Found ) THEN
    CALL Fatal(Caller,'Give > Fourier Series Components < !')
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
        FourierName = TRIM(VarName)//' Fourier Cos'//I2S(i/2)
      ELSE
        FourierName = TRIM(VarName)//' Fourier Sin'//I2S(i/2)
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
  CALL Info(Caller, Message, Level=8 )


  ! Check whether the exact integration of Fourier coefficients is used
  ! Always use it by default. 
  !-----------------------------------------------------------------------------
  ExactIntegration = .NOT. GetLogical(SolverParams, 'Inexact Integration', Found)
  IF (ExactIntegration) THEN
    CALL Info(Caller,'Using exact integration')    
  ELSE
    SimpsonsRule = GetLogical( SolverParams,'Simpsons Rule',Found )
    IF( SimpsonsRule ) THEN
      CALL Info(Caller,'Using Simpsons rule for integration')
    ELSE
      CALL Info(Caller,'Using trapetsoidal rule for integration')    
    END IF
  END IF

  ! Perform discrete fourier transform (DFT)
  ! If the cycle has come to end, then continue to calculate the losses.
  !-----------------------------------------------------------------------------
  EndCycle = .FALSE.
  CALL LocalFourierTransform( EndCycle )
  at1 = RealTime()

  WRITE( Message,'(A,ES12.3)') 'Fourier transform time: ',at1-at0
  CALL Info(Caller, Message, Level=5 )
  IF( .NOT. EndCycle ) RETURN

  CALL DefaultInitialize()

  IF( NVar > 1 ) THEN
    CALL Info(Caller,'Allocating multiple r.h.s. vectors',Level=12)
    IF(.NOT. ALLOCATED( OtherRhs ) ) THEN
      ALLOCATE( OtherRhs ( SIZE( Solver % Matrix % Rhs ), NVar - 1) )
    END IF
    OtherRhs = 0.0_dp
  END IF
  
  CALL BulkAssembly( Ncomp )

  at2 = RealTime()
  WRITE( Message,'(A,ES12.3)') 'Assembly time: ',at2-at1
  CALL Info( Caller, Message, Level=5 )

  !------------------------------------------------------------------------------     
  IF( NodalField ) THEN
    CALL ListAddLogical( SolverParams,'Linear System Compute Change',.FALSE.) 
    SaveRhs => Solver % Matrix % Rhs   
    DO icomp=1, Nvar
      IF(icomp > 1 ) THEN
        Solver % Matrix % Rhs => OtherRhs(:,icomp-1)
      END IF
      Norm = DefaultSolve()
      CompVars(icomp) % Var % Values = LossVar % Values
    END DO
    Solver % Matrix % Rhs => SaveRhs    
    IF(Nvar > 1) DEALLOCATE( OtherRhs ) 
  END IF
  
  at3 = RealTime()
  WRITE( Message,'(A,ES12.3)') 'Solution time: ',at3-at2
  CALL Info( Caller, Message, Level=5 )

  ! Print and save the lumped results
  !-------------------------------------------------------------------------------
  CALL CommunicateLosess()

  ! Let's use the total loss as a norm if no global linear system was solved.
  IF(.NOT. NodalField ) THEN
    Solver % Variable % Values = TotalLoss 
  END IF
  
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
    ! when choosing saving so that it coincides with the end of Fourier cycle.
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
      CALL Info(Caller,'Finishing the timestep using left rule',Level=8)
    ELSE      
      time = GetTime()

      ! Figure out the time to start Fourier analysis
      ! ----------------------------------------------
      n0 = GetInteger(SolverParams,'Fourier Start Timestep', Found )
      IF ( Found ) THEN
        i = GetTimestep() 
        IF( i < n0 ) THEN
          CALL Info(Caller,'Fourier transform not yet active, returning...',Level=6)
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
            IF( TimesVisited == 1 ) t0 = MAX(time-dt,0.0_dp)
          END IF
        END IF
      END IF

      ! Subtract the start time so that we start conveniently from zero
      time = time - t0

      ! Compute the fraction of the timestep needed
      LeftRule = .FALSE.
      IF( time <= 0.0 ) THEN
        ratio = 0.0_dp
        CALL Info(Caller,'Fourier transform not yet active, returning...',Level=6)
        RETURN
      ELSE IF( time - dt < 0.0 ) THEN
        ratio = time / dt
        LeftRule = .TRUE.
        CALL Info(Caller,'Starting Fourier transform cycle',Level=6)
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
            CALL Fatal(Caller,'Cannot use left and right rule at the same time!')
          END IF

          ratio = 1.0_dp - (time-NINT(time/tcycle)*tcycle)/dt

          !PRINT *,'RightCycle:',ratio, tcycle, time / tcycle, INT(time/tcycle), CurrentCycle, &
          !    CurrentCycle * tcycle, time, &
          !    CurrentCycle * tcycle - time > 0.0_dp 
                     
          CALL Info(Caller,'Finising Fourier transform cycle',Level=6)
          
          ! Return a True flag so that we know that we might need to start also the cycle.
          EndCycle = .TRUE.        
          PreviousCycle = CurrentCycle
        END IF
      END IF
    END IF
      
    !PRINT *,'RatioCycle:',ratio, time, tcycle, CurrentCycle, EndCycle, RightRule, LeftRule

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
    CALL Info(Caller,Message,Level=8)

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
    TYPE(GaussIntegrationPoints_t), TARGET :: IP
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
    LOGICAL :: DG, Erroneous, MaterialExponents
    
    SAVE Nodes
    
    SaveRHS => Solver % Matrix % Rhs

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

    
    MaterialExponents = ListCheckPrefixAnyMaterial( Model,'Harmonic Loss Frequency Exponent') 
    IF(.NOT. MaterialExponents) THEN
      CALL GetLossExponents(SolverParams,FreqPower,FieldPower,Ncomp)
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
          CALL Fatal(Caller,'Element on target mesh cannot be associated')
        END IF
        nt = GetElementDOFs( Indeces, TargetElement, TargetSolverPtr )
      END IF

      ! Integrate over local element:
      ! -----------------------------
      IP = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      BodyId = Element % BodyId
      Material => GetMaterial()

      IF(MaterialExponents) THEN
        CALL GetLossExponents(Material,FreqPower,FieldPower,Ncomp)
      END IF

      DO t=1,IP % n
        IF( DirectField ) THEN
          ! For direct field we don't need the curl i.e. no dBasisdx needed
          Found = ElementInfo( Element, Nodes, IP % u(t), &
              IP % v(t), IP % w(t), detJ, Basis )
        ELSE
          Found = ElementInfo( Element, Nodes, IP % u(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
              RotBasis = RotWBasis, USolver = TargetSolverPtr )           
        END IF        
        Weight = IP % s(t) * detJ

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
          Found2 = .FALSE.
          DO icomp = 1, Ncomp
            LossCoeff(icomp) = ListGetFun( Material,&
                'Harmonic Loss Coefficient '//I2S(icomp),FourierFreq, Found)      
            IF( Found ) Found2 = .TRUE.
          END DO
          IF(.NOT. Found2 ) CYCLE

          ! For even j we have cosine series, for odd sine series
          Component => FourierVars(j) % Var % Values

          IF( DirectField ) THEN
            IF( tdofs == 1 ) THEN            
              ElemField(1:nd) = Component( FourierPerm( Indeces(1:nd) ) )
              ValAtIp = SUM( Basis(1:nd) * ElemField(1:nd) )
            ELSE
              CurlAtIp = 0.0_dp
              DO k=1,tdofs
                ElemField(1:nd) = Component( tdofs * ( FourierPerm( Indeces(1:nd))-1) + k )
                CurlAtIp(k) = SUM( Basis(1:nd) * ElemField(1:nd) )
              END DO
            END IF
              
          ELSE IF ( AVField ) THEN
            ElemField(1:nt) = Component( FourierPerm( Indeces(1:nt) ) )

            ! Here we assume that the solution contains both a nodal interpolation
            ! field (scalar) and edge interpolation part (vector). The gradient is 
            ! applied to the scalar while the curl operator is applied to the vector

            CurlAtIp(1) = SUM( ElemField(n+1:nt) * RotWBasis(1:(nt-n),1) )
            CurlAtIp(2) = SUM( ElemField(n+1:nt) * RotWBasis(1:(nt-n),2) )
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
          ELSE ! dim == 2
            ElemField(1:nd) = Component( FourierPerm( Indeces(1:nd) ) )
            IF(FlipActive) THEN
              DO i=1,nd
                IF( CurrentModel % Mesh % PeriodicFlip(Indeces(i)) ) ElemField(i) = -ElemField(i)
              END DO
            END IF           
            CurlAtIp(1) =  SUM( ElemField(1:nd) * dBasisdx(1:nd,2) )
            CurlAtIp(2) = -SUM( ElemField(1:nd) * dBasisdx(1:nd,1) )
            CurlAtIp(3) = 0.0_dp
          END IF

          IF(.NOT. ( DirectField .AND. tdofs == 1) ) THEN
            ValAtIp = SQRT( SUM( CurlAtIP ** 2 ) )
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
        CALL LUdecomp(STIFF,n,pivot,Erroneous)
        IF (Erroneous) CALL Fatal('FourierLoss', 'LU-decomposition fails')
        DO icomp = 1, Nsum
          CALL LUSolve(n,STIFF,FORCE(icomp,:),pivot)
          CompVarsE(icomp) % Var % Values(EPerm(Element % DGIndexes(1:n))) = FORCE(icomp,1:n)
        END DO
      END IF
  
    END DO
       
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

    INTEGER :: i,j,k,NoSlices
    CHARACTER(LEN=MAX_NAME_LEN) :: LossesFile
    
    IF( ParEnv % PEs > 1 ) THEN
      NoSlices = MAX(1,ListGetInteger( Model % Simulation,'Number Of Slices', Found ) )   
      IF( NoSlices > 1 ) THEN
        CALL Info(Caller,'Averaging losses over slices',Level=8)        
      END IF      
      DO k=1,Ncomp
        DO j=1,FourierDofs
          SeriesLoss(k,j) = ParallelReduction(SeriesLoss(k,j)) / NoSlices
        END DO
        DO j=1,Model % NumberOfBodies
          BodyLoss(k,j) = ParallelReduction(BodyLoss(k,j)) / NoSlices 
        END DO
      END DO
    END IF
      
    ! Sum up the losses over components
    DO k=1,Ncomp
      CompLoss(k) = SUM( SeriesLoss(k,:) )
    END DO
    TotalLoss = SUM( CompLoss )
        
    CALL ListAddConstReal( Model % Simulation,TRIM(Pref)//' fourier loss total',TotalLoss )
    DO k=1,Ncomp
      CALL ListAddConstReal( Model % Simulation,TRIM(Pref)//' fourier loss '//I2S(k),CompLoss(k) )
    END DO
      
    ! Output of the losses on screen    
    ! First the component losses, then the body losses
    DO k=1,NComp
      CALL Info(Caller,'Wavewise Fourier loss for component: '//I2S(k),Level=6)
      
      DO j=1,FourierDofs 
        IF( j == 1 ) THEN
          WRITE( Message,'(A,ES12.3)') 'CONST : ',SeriesLoss(k,j)
          CALL Info(Caller, Message, Level=6 )
        ELSE
          i = j/2
          IF( MODULO( j,2 ) == 0 ) THEN
            WRITE( Message,'(A,I0,A,ES12.3)') 'COS_',i,' : ',SeriesLoss(k,j)
            CALL Info(Caller, Message, Level=6 )
          ELSE
            WRITE( Message,'(A,I0,A,ES12.3)') 'SIN_',i,' : ',SeriesLoss(k,j)
            CALL Info(Caller, Message, Level=6 )
          END IF
        END IF
      END DO
    
      CALL Info(Caller,'Bodywise Fourier loss for component: '//I2S(k),Level=6)
      
      DO j=1,Model % NumberOfBodies
        IF( BodyLoss(k,j) < TINY( CompLoss(k) ) ) CYCLE
        WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyLoss(k,j)
        CALL Info(Caller, Message, Level=6 )
      END DO

      WRITE( Message,'(A,ES12.3)') 'Total component '//I2S(k)//' loss: ',CompLoss(k)
      CALL Info(Caller,Message, Level=5 )
    END DO
    
    WRITE( Message,'(A,ES12.3)') 'Total loss: ',TotalLoss
    CALL Info(Caller,Message, Level=5 )
    
    IF( Parenv % MyPe == 0 ) THEN
      LossesFile = ListGetString(SolverParams,'Fourier Loss Filename',Found )
      IF( Found ) THEN
        OPEN (10, FILE=LossesFile)
        WRITE( 10,'(A)')  '!body_id   loss(1)   loss(2) ....'        
        DO j=1,Model % NumberOfBodies          
          IF( SUM( BodyLoss(1:Ncomp,j) ) <= TINY( TotalLoss ) ) CYCLE
          WRITE( 10,'(I6)',ADVANCE='NO') j
          DO i=1,Ncomp-1
            WRITE( 10,'(ES15.6)',ADVANCE='NO') BodyLoss(i,j)            
          END DO
          WRITE( 10,'(ES15.6)') BodyLoss(Ncomp,j)            
        END DO
        CALL Info(Caller,'Fourier losses for bodies was saved to file: '//TRIM(LossesFile),Level=6 )
        CLOSE(10)
      END IF

      LossesFile = ListGetString(SolverParams,'Series Loss Filename',Found )
      IF( Found ) THEN
        OPEN (10, FILE=LossesFile)
        DO j=1,FourierDofs/2
          WRITE( 10,'(2ES15.6)') SUM(SeriesLoss(1:Ncomp,2*j-1)), SUM(SeriesLoss(1:Ncomp,2*j))
        END DO
        CALL Info(Caller,'Series losses for bodies was saved to file: '//TRIM(LossesFile),Level=6 )
        CLOSE(10)
      END IF
    END IF

    ! For debugging purposes     
    IF( InfoActive(25) ) THEN
      PRINT *,'Fourier components:'
      DO i=1,FourierDofs
        CALL VectorValuesRange(FourierVars(i) % Var % Values, &
            SIZE(FourierVars(i) % Var % Values),'F'//I2S(i))
      END DO
      
      PRINT *,'Loss components:'
      DO i=1,NComp
        CALL VectorValuesRange(CompVars(i) % Var % Values, &
            SIZE(CompVars(i) % Var % Values),'L'//I2S(i))
      END DO
    END IF

    !------------------------------------------------------------------------------
  END SUBROUTINE CommunicateLosess
  !------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE FourierLossSolver
!------------------------------------------------------------------------------

!> \}
