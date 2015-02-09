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


  SolverParams => GetSolverParams()

  VarName = ListGetString( SolverParams,'Target Variable',Found )
  IF(.NOT. Found ) THEN
    CALL Fatal('FourierLossSolver','Give > Target Variable < !')
  END IF

  IF( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
    CALL ListAddString( SolverParams,'Variable','Fourier Loss Linear' )
  END IF

  CALL ListAddString( SolverParams,'Exported Variable 1','Fourier Loss Quadratic' )

  ! We are really solving using DG so no need for funny initialization
  IF (GetLogical(GetSolverParams(),'Discontinuous Galerkin',Found)) RETURN

  ! If elemental fields are not requested don't compute them
  ElementalField = GetLogical( GetSolverParams(), 'Calculate Elemental Fields', Found)
  IF(Found .AND. .NOT. ElementalField) RETURN

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

  DGSolverParams => NULL()
  CALL ListAddLogical( DGSolverParams, 'Discontinuous Galerkin', .TRUE. )
  Solvers(n+1) % Values => DGSolverParams
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
  CALL ListAddString( DGSolverParams, 'Variable', 'Fourier Loss Linear e' )
  CALL ListAddString( DGSolverParams, 'Exported Variable 1','Fourier Loss Quadratic e')

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
  TYPE(ValueList_t),POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, VarString
  INTEGER :: dim, Nsize, FourierDofs, i, nlen, mlen
  LOGICAL :: Found
  INTEGER :: TimesVisited = 0
  INTEGER, POINTER :: FourierPerm(:)
  REAL(KIND=dp) :: Norm, Omega, TotalLoss(2)
  REAL(KIND=dp), POINTER :: FourierField(:)
  REAL(KIND=dp) :: at0,at1,at2,at3,CPUTime,RealTime
  REAL(KIND=dp), ALLOCATABLE :: BodyLoss(:,:), ComponentLoss(:,:)
  TYPE(Variable_t), POINTER :: TargetVar, FourierVar, ElVar, ElVar2, LossVar, LossVar2
  REAL(KIND=dp), POINTER :: TargetField(:), PrevTargetField(:,:), SaveRhs(:), Rhs2(:)
  LOGICAL :: EndCycle, FourierOutput, SimpsonsRule, ExactIntegration, &
      ElementalField, AvField
  TYPE(Solver_t), POINTER :: TargetSolverPtr

  SAVE TimesVisited
  !-------------------------------------------------------------------------------

  CALL Info( 'FourierLossSolver', '-------------------------------------',Level=4 )
  CALL Info( 'FourierLossSolver', 'Computing Fourier losses             ',Level=4 )
  CALL Info( 'FourierLossSolver', '-------------------------------------',Level=4 )

  TimesVisited = TimesVisited + 1
  dim = CoordinateSystemDimension()

  at0 = RealTime()

  IF(.NOT. ListCheckPresentAnyMaterial( Model,'Harmonic Loss Linear Coefficient') ) THEN
    CALL Warn('FourierLossSolver',&
        'Some material should have > Harmonic Loss Linear Coefficient < ')
  END IF

  IF(.NOT. ListCheckPresentAnyMaterial( Model,'Harmonic Loss Quadratic Coefficient') ) THEN
    CALL Warn('FourierLossSolver',&
        'Some material should have > Harmonic Loss Quadratic Coefficient < ')
  END IF

  ! Check that we have something to compute
  !------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) THEN
    CALL Fatal('FourierLossSolver','Solver has no matrix associated')
  END IF
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN  
  SolverParams => GetSolverParams()

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
    

  LossVar => Solver % Variable
  LossVar2 => VariableGet( Solver % Mesh % Variables,&
      'Fourier Loss Quadratic' )

  ! Check for Elemental (Discontinuous Galerkin) Field 
  !------------------------------------------------------------------
  ElVar => VariableGet( Solver % Mesh % Variables,&
      'Fourier Loss Linear e' )
  ElVar2 => VariableGet( Solver % Mesh % Variables,&
      'Fourier Loss Quadratic e' )
  ElementalField = ASSOCIATED( ElVar ) .OR. ASSOCIATED( ElVar2 )

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

  ! Fetch the Fourier field that includes the components of the transform
  ! If it has not been created before, create it now
  !------------------------------------------------------------------------
  Varname = TRIM( VarName )//' Fourier'
  FourierVar => VariableGet( Solver % Mesh % Variables, Varname)
  IF(.NOT. ASSOCIATED( FourierVar ) ) THEN    
    ! Allocate the components: 1st one is related to constant part
    ! then even ones to cosine and odd ones to sine. 
    FourierDofs = 1 + 2 * ListGetInteger( SolverParams,'Fourier Series Components',Found)
    IF(.NOT. Found ) THEN
      CALL Fatal('FourierLossSolver','Give > Fourier Series Components < !')
    END IF
    NULLIFY( FourierField ) 
    ALLOCATE( FourierField( FourierDofs * Nsize ) )
    FourierField = 0.0_dp
    FourierOutput = ListGetLogical( SolverParams,'Fourier Series Output',Found) 
    CALL VariableAddVector( Solver % Mesh % Variables, Solver % Mesh, Solver, &
        VarName, FourierDofs * TargetVar % Dofs, FourierField, TargetVar % Perm, FourierOutput )
    FourierVar => VariableGet( Solver % Mesh % Variables, Varname)
  END IF
  FourierField => FourierVar % Values
  FourierPerm => FourierVar % Perm
  FourierDofs = FourierVar % Dofs / TargetVar % Dofs

  ALLOCATE( ComponentLoss(2,FourierDofs), &
      BodyLoss(2,Model % NumberOfBodies) )

  ! Fetch frequency that is used both for Fourier transform and loss computation
  !-----------------------------------------------------------------------------
  Omega = GetAngularFrequency(Found=Found)
  IF( .NOT. Found ) THEN  
    CALL Fatal('FourierLosses','The Fourier transform requires frequency!')
  END IF
  WRITE( Message,'(A,ES12.3)') 'Base frequency for Fourier transform: ',Omega / (2*PI)
  CALL Info('FourierLossSolver', Message, Level=8 )


  ! Check whether the exact integration of Fourier coefficients is used
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
    CALL ListAddConstReal( Model % Simulation,'res: fourier loss linear',0.0_dp )
    CALL ListAddConstReal( Model % Simulation,'res: fourier loss quadratic',0.0_dp )
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
  ALLOCATE( Rhs2 ( SIZE( Solver % Matrix % Rhs ) ) )
  Rhs2 = 0.0_dp
  SaveRhs => Solver % Matrix % Rhs    

  CALL BulkAssembly()

  ! These three are probably not needed for this solver of limited features
  ! They are left here for possible future needs.
  ! CALL DefaultFinishBulkAssembly()
  ! CALL DefaultFinishAssembly()
  ! CALL DefaultDirichletBCs()      
  at2 = RealTime()
  WRITE( Message,'(A,ES12.3)') 'Assembly time: ',at2-at1
  CALL Info( 'FourierLossSolver', Message, Level=5 )

  !------------------------------------------------------------------------------     
  SaveRhs => Solver % Matrix % Rhs
  Solver % Matrix % Rhs => Rhs2
  Norm = DefaultSolve()
  LossVar2 % Values = LossVar % Values
  Solver % Matrix % Rhs => SaveRhs
  DEALLOCATE( Rhs2 ) 

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
      FourierField = 0.0_dp
    END IF
    InitPending = .FALSE.

    IF( EndCycle ) THEN
      ratio = 1.0 - ratio
      IF( ratio < RatioEps ) THEN
        InitPending = .TRUE.
        RETURN
      END IF
      FourierField = 0.0_dp
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
      Component => FourierField(j::FourierDofs)

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
  SUBROUTINE BulkAssembly()
    !------------------------------------------------------------------------------
    INTEGER :: elem,i,j,k,p,q,n,nd,nt,t,BodyId
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element, TargetElement
    REAL(KIND=dp) :: weight,coeff,coeff2,detJ,GradAtIp(3,3),CurlAtIP(3),ValAtIp
    LOGICAL :: Found, Found2
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: Freq, FourierFreq, LossCoeff, FreqPower, FieldPower, &
        LossCoeff2, FreqPower2, FieldPower2
    REAL(KIND=dp), POINTER :: Component(:)
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), FORCE2(:)
    REAL(KIND=dp), ALLOCATABLE :: ElemField(:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:), RotWBasis(:,:)
    INTEGER, ALLOCATABLE :: Indeces(:), Pivot(:)

    SAVE Nodes

    SaveRhs => Solver % Matrix % Rhs

    n = MAX(Solver % Mesh % MaxElementDOFs,Solver % Mesh % MaxElementNodes)
    ALLOCATE( STIFF(n,n), FORCE(n), FORCE2(n), Pivot(n) )
    ALLOCATE( Basis(n), dBasisdx(n,3) )

    IF (AVField) THEN
      nt = MAX(TargetSolverPtr % Mesh % MaxElementDOFs,TargetSolverPtr % Mesh % MaxElementNodes)
      ALLOCATE( WBasis(nt,3), RotWBasis(nt,3), Indeces(nt), ElemField(nt) )
    ELSE
      ALLOCATE( ElemField(n), Indeces(n) )
    END IF

    Freq = Omega / (2*PI)

    FreqPower = GetCReal( SolverParams,'Harmonic Loss Linear Frequency Exponent',Found )
    IF(.NOT. Found ) FreqPower = 1.0_dp

    FieldPower = GetCReal( SolverParams,'Harmonic Loss Linear Exponent',Found ) 
    IF( .NOT. Found ) FieldPower = 2.0_dp

    FreqPower2 = GetCReal( SolverParams,'Harmonic Loss Quadratic Frequency Exponent',Found )
    IF(.NOT. Found ) FreqPower2 = 2.0_dp

    FieldPower2 = GetCReal( SolverParams,'Harmonic Loss Quadratic Exponent',Found ) 
    IF( .NOT. Found ) FieldPower2 = 2.0_dp

    ! For economical reasons include the SQRT function used for taking the length of quadratic sum
    FieldPower = FieldPower / 2.0_dp
    FieldPower2 = FieldPower / 2.0_dp

    ! Sum over the loss for each frequency, each body, and for the combined effect
    ComponentLoss = 0.0_dp
    BodyLoss = 0.0_dp
    TotalLoss = 0.0_dp

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
      FORCE2 = 0.0_dp

      BodyId = GetBody()
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
          LossCoeff = ListGetFun( Material,'Harmonic Loss Linear Coefficient',FourierFreq, Found)      
          LossCoeff2 = ListGetFun( Material,'Harmonic Loss Quadratic Coefficient',FourierFreq, Found2 )      
          IF(.NOT. (Found .OR. Found2 ) ) CYCLE

          ! For even j we have cosine series, for odd sine series
          Component => FourierField(j::FourierDofs)

          IF ( AVField ) THEN
            ElemField(1:nt) = Component( FourierPerm( Indeces(1:nt) ) )

            ! Here we assume that the solution contains both a nodal interpolation
            ! field (scalar) and edge interpolation part (vector). The gradient is 
            ! applied to the scalar while the curl operator is applied to the vector

            !GradAtIp(1) =  SUM( ElemField(1:n) * dBasisdx(1:n,1) )
            !GradAtIp(2) =  SUM( ElemField(1:n) * dBasisdx(1:n,2) ) 
            !IF ( dim > 2 ) &
            !    GradAtIp(3) =  SUM( ElemField(1:n) * dBasisdx(1:n,3) )

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

            !GradAtIp(1) =  SUM( ElemField(1:nd) * dBasisdx(1:nd,1) )
            !GradAtIp(2) =  SUM( ElemField(1:nd) * dBasisdx(1:nd,2) )

            CurlAtIp(1) =  SUM( ElemField(1:nd) * dBasisdx(1:nd,2) )
            CurlAtIp(2) = -SUM( ElemField(1:nd) * dBasisdx(1:nd,1) )
          END IF

          ValAtIp = SUM( CurlAtIP(1:dim) ** 2 ) 

          ! linear losses 
          Coeff = Weight * LossCoeff * ( FourierFreq ** FreqPower ) * ( ValAtIp ** FieldPower )
          FORCE(1:nd) = FORCE(1:nd) + Coeff * Basis(1:nd)          
          ComponentLoss(1,j) = ComponentLoss(1,j) + Coeff 
          BodyLoss(1,BodyId) = BodyLoss(1,BodyId) + Coeff

          Coeff2 = Weight * LossCoeff2 * ( FourierFreq ** FreqPower2 ) * ( ValAtIp ** FieldPower2 )
          FORCE2(1:nd) = FORCE2(1:nd) + Coeff2 * Basis(1:nd)          
          ComponentLoss(2,j) = ComponentLoss(2,j) + Coeff2 
          BodyLoss(2,BodyId) = BodyLoss(2,BodyId) + Coeff2
        END DO

      END DO

      ! Update global matrices from local matrices 
      !------------------------------------------------------------------------------
      CALL DefaultUpdateEquations( STIFF, FORCE(1:nd) )
      Solver % Matrix % Rhs => Rhs2
      CALL DefaultUpdateForce( FORCE2(1:nd) )
      Solver % Matrix % Rhs => SaveRhs

      ! After this the STIFF and FORCE are corrupted 
      IF( ElementalField ) THEN
        CALL LUdecomp(STIFF,n,pivot)
        CALL LUSolve(n,STIFF,FORCE,pivot)
        ElVar % Values(ElVar % Perm(Element % DGIndexes(1:n))) = FORCE(1:n)
        CALL LUSolve(n,STIFF,FORCE2,pivot)
        ElVar2 % Values(ElVar2 % Perm(Element % DGIndexes(1:n))) = FORCE2(1:n)
      END IF

    END DO

    ! Assembly of the face terms when using DG:
    !------------------------------------------------------------------
    IF (GetLogical(SolverParams,'Discontinuous Galerkin',Found)) THEN
      IF (GetLogical(SolverParams,'Average Within Materials',Found)) THEN
        FORCE = 0.0d0
        CALL AddLocalFaceTerms( STIFF, FORCE )
      END IF
    END IF

    DEALLOCATE( ElemField, STIFF, FORCE, Pivot, Basis, &
        dBasisdx, Indeces )
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

    DO k=1,2
      DO j=1,FourierDofs
        ComponentLoss(k,j) = ParallelReduction(ComponentLoss(k,j)) 
      END DO

      DO j=1,Model % NumberOfBodies
        BodyLoss(k,j) = ParallelReduction(BodyLoss(k,j))
      END DO

      TotalLoss(k) = SUM( ComponentLoss(k,:) )
    END DO

    CALL ListAddConstReal( Model % Simulation,'res: fourier loss linear',TotalLoss(1) )
    CALL ListAddConstReal( Model % Simulation,'res: fourier loss quadratic',TotalLoss(2) )

    ! Output of the losses on screen    
    ! First the component losses, the the body losses
    DO k=1,2
      IF( k == 1 ) THEN
        CALL Info('FourierLossSolver','Fourier loss linear by components',Level=6)
      ELSE
        CALL Info('FourierLossSolver','Fourier loss quadratic by components',Level=6)
      END IF
      DO j=1,FourierDofs 
        IF( ComponentLoss(k,j) < TINY(TotalLoss(k)) ) CYCLE
        IF( j == 1 ) THEN
          WRITE( Message,'(A,I0,A,ES12.3)') 'COS_',0,' : ',ComponentLoss(k,j)
          CALL Info('FourierLossSolver', Message, Level=6 )
        ELSE
          i = j/2
          IF( MODULO( j,2 ) == 0 ) THEN
            WRITE( Message,'(A,I0,A,ES12.3)') 'COS_',i,' : ',ComponentLoss(k,j)
            CALL Info('FourierLossSolver', Message, Level=6 )
          ELSE
            WRITE( Message,'(A,I0,A,ES12.3)') 'SIN_',i,' : ',ComponentLoss(k,j)
            CALL Info('FourierLossSolver', Message, Level=6 )
          END IF
        END IF
      END DO
    END DO

    IF( k == 1 ) THEN
      CALL Info('FourierLossSolver','Fourier loss linear by bodies',Level=6)
    ELSE
      CALL Info('FourierLossSolver','Fourier loss quadratic by bodies',Level=6)
    END IF

    DO j=1,Model % NumberOfBodies
      IF( BodyLoss(k,j) < TINY( TotalLoss(k) ) ) CYCLE
      WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyLoss(k,j)
      CALL Info('FourierLossSolver', Message, Level=6 )
    END DO

    WRITE( Message,'(A,ES12.3)') 'Total loss: ',TotalLoss(k)
    CALL Info('FourierLossSolver',Message, Level=5 )

    IF( Parenv % MyPe == 0 ) THEN
      LossesFile = ListGetString(SolverParams,'Fourier Loss Filename',Found )
      IF( Found ) THEN
        OPEN (10, FILE=LossesFile)
        WRITE( 10,'(A)')  '!body_id   fourier(1)      fourier(2)'
        DO j=1,Model % NumberOfBodies
          IF( BodyLoss(k,j) < TINY( TotalLoss(k) ) ) CYCLE
          WRITE( 10,'(I10,ES17.9)') j, BodyLoss(k,j)
        END DO
        CALL Info('FourierLossSolver', &
            'Fourier losses for bodies was saved to file: '//TRIM(LossesFile),Level=6 )
        CLOSE(10)
      END IF
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
