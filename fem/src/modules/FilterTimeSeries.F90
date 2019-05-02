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
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 12 Feb 2007
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Subroutine for avaraging transient data over timesteps. Also more complicated
!> operations may be performed, such as definition of Fourier series components.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE FilterTimeSeries( Model,Solver,dtime,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dtime
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  TYPE(Variable_t), POINTER :: Var, MeanVar, OldVar
  LOGICAL :: SubroutineVisited=.FALSE., GotCoeff, GotOper, GotVar, GotStuff, Found, &
       PrevFieldExists, EndRatio
  INTEGER :: i, j, k, n, InstDofs, MeanDofs, NoVar, Nsize, LoopSize, Nsine, Ncosine, &
      Nseries, TimesVisited=0, n0, n1
  INTEGER, POINTER :: Perm(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: Oper, OldOper, OperName, VarName, MeanVarName, Name, &
      FilterName, OldVarName, tmpname
  REAL(KIND=dp) :: st, t0, t1, time, Phase, Ratio, Weight, Relax, val, cnew, cold, &
      freq, fcoeff, dt, q, prevval
  REAL(KIND=dp) :: CumWeight(99)
  REAL(KIND=dp), POINTER :: MeanField(:), InstField(:), Component(:), PrevField(:,:)
#ifndef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: RealTime
#endif

  SAVE SubroutineVisited, TimesVisited, CumWeight


!------------------------------------------------------------------------------

  CALL Info( 'FilterTimeSeries', '-----------------------------------------', Level=4 )
  CALL Info( 'FilterTimeSeries', 'Filtering variables over time', Level=4 )
  CALL Info( 'FilterTimeSeries', '-----------------------------------------', Level=4 )

  IF(.NOT. SubroutineVisited) THEN
    CumWeight = 0.0_dp    
    SubroutineVisited = .TRUE.
  END IF
  TimesVisited = TimesVisited + 1

  Var => VariableGet( Model % Variables, 'Time' )
  IF( ASSOCIATED(Var) ) time = Var % Values(1)
    
  !------------------------------------------------------------------------------
  ! Go through the variables and compute the desired statistical data
  !------------------------------------------------------------------------------
  NULLIFY(OldVar)
  
  DO NoVar = 1,99

    GotOper = .FALSE.
    GotStuff = .FALSE.
    EndRatio = .FALSE.

    Name = ParameterName('Variable',NoVar)
    VarName = ListGetString( Solver % Values, TRIM(Name), Found )
    IF(Found) THEN      
      Var => VariableGet( Model % Variables, TRIM(VarName) )
      IF ( .NOT. ASSOCIATED( Var ) )  THEN
        CALL Warn('SaveData','The desired variable '//TRIM(VarName)//' does not exist!')
        CYCLE
      ELSE
        OldVar => Var
        OldVarName = VarName
        GotStuff = .TRUE.
      END IF
    ELSE IF(ASSOCIATED(OldVar)) THEN
      Var => OldVar
      VarName = OldVarName
    ELSE
      EXIT
    END IF
   
    ! Other powers than ^1
    Name = ParameterName('Operator',NoVar)
    Oper = ListGetString(Solver % Values,TRIM(Name),Found)
    IF(.NOT. Found) Oper = 'none'
    GotStuff = GotStuff .OR. Found

    dt = dtime
    ratio = 1.0_dp

    ! Start after given wall-clock time
    t0 = GetCReal(Solver % Values,'Start Real Time', Found )
    IF ( Found ) THEN
       IF( RealTime() < t0 ) THEN
          ratio = 0.0_dp
       END IF
    END IF
    
    ! Start after relative given wall-clock time
    q = GetCReal(Solver % Values,'Real Time Max Fraction', Found )
    IF( Found ) THEN
       t0 = GetCReal(Model % Simulation,'Real Time Max' )
       IF( RealTime() < q * t0 ) THEN
          ratio = 0.0_dp
       END IF
    END IF

    ! Start watch time
    Name = ParameterName('Start Time',NoVar)
    t0 = GetCReal(Solver % Values,TRIM(Name), Found )
    IF ( Found ) THEN
      IF( time < t0 ) THEN
        ratio = 0.0_dp
      ELSE IF( time - dt < t0 ) THEN
        ratio = (time - t0) / dt
      END IF
      GotStuff = .TRUE.
    END IF

    ! Stop watch time
    Name = ParameterName('Stop Time',NoVar)
    t1 = GetCReal(Solver % Values,TRIM(Name), Found )
    IF ( Found ) THEN
      IF( time - dt > t1 ) THEN
        ratio = 0.0_dp
      ELSE IF( time > t1 ) THEN
        ratio = (t1+dt-time) / dt
        EndRatio = .TRUE.
      END IF
      GotStuff = .TRUE.
    END IF

    ! Start timestep number
    Name = ParameterName('Start Timestep',NoVar)
    n0 = GetInteger(Solver % Values,TRIM(Name), Found )
    IF ( Found ) THEN
       IF( TimesVisited < n0 ) ratio = 0.0_dp
       dt = 1.0_dp
       GotStuff = .TRUE.
    END IF

    ! Stop timestep number
    Name = ParameterName('Stop Timestep',NoVar)
    n1 = GetInteger(Solver % Values,TRIM(Name), Found )
    IF ( Found ) THEN
      IF( TimesVisited > n1 ) ratio = 0.0_dp
      dt = 1.0_dp
      GotStuff = .TRUE.
    END IF

    ! Reset mean after elapsed time
    Name = ParameterName('Reset Interval',NoVar)
    st = GetCReal(Solver % Values,TRIM(Name), Found )
    IF ( Found ) THEN
      IF ( INT((time-t0)/st) /= INT((time-t0-dt)/st) ) THEN
        !ratio = (time-t0-INT((time-t0)/st)*st)/dt
        MeanField = 0.0_dp
        CumWeight(NoVar) = 0.0_dp
      END IF
      GotStuff = .TRUE.
    END IF

    ! Relax by giving time scale
    Relax = 1.0_dp
    Name = ParameterName('Decay Time',NoVar)
    st = GetCReal(Solver % Values,TRIM(Name), Found )
    IF ( Found ) THEN
       Relax = EXP(-dt/st)
       GotStuff = .TRUE.
    END IF

    ! Relax by giving time scale
    Relax = 1.0_dp
    Name = ParameterName('Decay Timesteps',NoVar)
    n0 = GetInteger(Solver % Values,TRIM(Name), Found )
    IF ( Found ) THEN
       Relax = EXP(-1.0_dp/n0)
       dt = 1.0_dp
       GotStuff = .TRUE.
    END IF


    ! Filter with externally provided function (sin, cos, ...)
    Phase = 1.0_dp
    Name = ParameterName('Time Filter',NoVar)
    st = GetCReal(Solver % Values,TRIM(Name), Found )
    IF(Found) THEN
       Phase = st
       GotStuff = .TRUE.
    END IF

    ! Explicitly perform sine and cosine series of given degree
    Name = ParameterName('Sine Series',NoVar)
    Nsine = GetInteger(Solver % Values,TRIM(Name), Found )
    Name = ParameterName('Cosine Series',NoVar)
    Ncosine = GetInteger(Solver % Values,TRIM(Name), Found )
    Nseries = MAX(Nsine,Ncosine)
    GotStuff = GotStuff .OR. (Nseries > 0)

    ! If no new control  parameters were obtained exit the loop
    IF(.NOT. GotStuff) EXIT

    IF(Nseries > 0) THEN
      Name = ParameterName('Frequency',NoVar)
      freq = GetCReal(Solver % Values,TRIM(Name),Found )
      IF(.NOT. Found) CALL Fatal('FilterTimeSeries','Fourier series requires frequency')


      ! Start watch time
      Name = ParameterName('Start Cycle',NoVar)
      t0 = GetCReal(Solver % Values,TRIM(Name), Found )
      IF ( Found ) THEN
         t0 = 1 / t0
         IF( time < t0 ) THEN
            ratio = 0.0_dp
         ELSE IF( time - dt < t0 ) THEN
            ratio = (time - t0) / dt
         END IF
      END IF
      
      ! Stop watch time
      Name = ParameterName('Stop Cycle',NoVar)
      t1 = GetCReal(Solver % Values,TRIM(Name), Found )
      IF ( Found ) THEN
         t1 = 1 / t1
         IF( time - dt > t1 ) THEN
            ratio = 0.0_dp
         ELSE IF( time > t1 ) THEN
            ratio = (t1+dt-time) / dt
            EndRatio = .TRUE.
         END IF
      END IF           
    END IF


    InstField => Var % Values
    InstDofs = Var % Dofs
    Perm => Var % Perm

    PrevFieldExists = .FALSE.
    IF( ListGetLogical(Solver % Values,'First Order', Found ) ) THEN
       PrevField => Var % PrevValues
       PrevFieldExists = ASSOCIATED(PrevField)
    END IF

    ! Set the variable name differently for sine and cosine series
    IF(Nsine > 0) THEN
      FilterName = 'Sine'
    ELSE IF(Ncosine > 0) THEN
      FilterName = 'Cosine'
    ELSE
      FilterName = ParameterName('Filter',NoVar)
    END IF
    WRITE (MeanVarName,'(A)') TRIM(FilterName)//' '//TRIM(VarName)


    MeanVar => VariableGet( Model % Variables,TRIM( MeanVarName ) )
    IF(ASSOCIATED(MeanVar)) THEN
      MeanField => MeanVar % Values
      MeanDofs = MeanVar % DOFs
      nsize = SIZE(InstField) * ( MeanDofs / InstDofs)
    ELSE
      SELECT CASE(Oper)

      CASE ('none','square','abs')
        MeanDofs = InstDofs
        IF( InstDofs > 1 .AND. Nseries > 0 ) THEN
          CALL Fatal('FilterTimeSeries',&
              'Sine and Cosine series are implemented for scalars only!')          
        END IF

      CASE('length')
        MeanDofs = 1

      CASE DEFAULT         
        IF(GotOper) THEN
          WRITE (Message,'(A,A)') 'Unknown operator: ',TRIM(Oper)
          CALL WARN('FilterTimeSeries',Message)
        END IF

      END SELECT        
      IF(Nseries > 0) MeanDofs = Nseries

      nsize = SIZE(InstField) * ( MeanDofs / InstDofs)
      ALLOCATE( MeanField(nsize) )
      MeanField = 0.0d0

      IF(ASSOCIATED(Perm)) THEN
        CALL VariableAdd( Var % PrimaryMesh % Variables, Var % PrimaryMesh, &
            Var % Solver, TRIM(MeanVarName), MeanDofs, MeanField, Perm )
        IF ( MeanDofs > 1 ) THEN
          n = LEN_TRIM( MeanVarName )
          DO j=1,Nseries
            tmpname = ComponentName( MeanVarName, j )
            Component => MeanField( j:nSize-Nseries+j:Nseries )
            CALL VariableAdd( Var % PrimaryMesh % Variables, Var % PrimaryMesh, &
                Var % Solver, tmpname, 1, Component, Perm )
          END DO
        END IF
      ELSE 
        CALL VariableAdd( Var % PrimaryMesh % Variables, Var % PrimaryMesh, &
            Var % Solver, TRIM(MeanVarName), MeanDofs, MeanField )
      END IF

      MeanVar => VariableGet( Model % Variables,TRIM( MeanVarName ) )
    END IF

    Weight = dt * ratio 
    CumWeight(NoVar) = Relax * CumWeight(NoVar)

    IF( ABS(Weight) < TINY(Weight) ) CYCLE


    ! weight for the filtered set
    cold = CumWeight(NoVar)/(CumWeight(NoVar)+Weight)
 
    ! For fourier series normalize the component appropriately as <sin^2>=<cos^2>=1/2
    IF(Nseries > 0) THEN
      fcoeff = 2*PI*freq*(time-t0)
      ! weight for the current timestep
      cnew = 2 * Weight/(CumWeight(NoVar)+Weight) 
      Loopsize = Nsize / Nseries
    ELSE
      ! weight for the current timestep
      cnew = phase * Weight/(CumWeight(NoVar)+Weight)
      LoopSize = Nsize
    END IF

   
    DO i=1,LoopSize
      
      SELECT CASE(OperName)
        
      CASE('square')
        val = InstField(i)**2
        IF( PrevFieldExists ) prevval = PrevField(i,1)**2

      CASE('abs')
        val = ABS(InstField(i))
       
      CASE('length')
        val = 0.0d0
        DO k=1,InstDofs
          val = val + InstField(InstDofs*(i-1)+k)**2.0_dp
        END DO
        val = SQRT(val)

        IF( PrevFieldExists ) THEN
           prevval = 0.0d0
           DO k=1,InstDofs
              prevval = prevval + PrevField(InstDofs*(i-1)+k,1)**2.0_dp
           END DO
           prevval = SQRT(prevval)
        END IF

      CASE DEFAULT 
        val = InstField(i)
        IF( PrevFieldExists ) prevval = PrevField(i,1)
        
      END SELECT
           
      ! Use 1st order integration scheme to account for the previous timestep
      !----------------------------------------------------------------------
      IF( PrevFieldExists ) THEN
         IF( EndRatio ) THEN
            val = ( ratio * val + (2-ratio) * prevval ) / 2
         ELSE
            val = ( ratio * prevval + (2-ratio) * val ) / 2
         END IF
      END IF


      ! Compute the mean fields
      !------------------------
      IF(Nseries == 0) THEN
        MeanField(i) = cold * MeanField(i) + cnew * val 
      ELSE        
        IF(Nsine > 0) THEN
          DO j=1,Nsine
            MeanField(Nsine*(i-1)+j) = cold * MeanField(Nsine*(i-1)+j) + &
                cnew * SIN(j*fcoeff) * val         
          END DO
        ELSE
          DO j=1,Ncosine
            MeanField(Ncosine*(i-1)+j) = cold * MeanField(Ncosine*(i-1)+j) + &
                cnew * COS(j*fcoeff) * val                     
          END DO
        END IF
      END IF

    END DO
    
    CumWeight(NoVar) = CumWeight(NoVar) + Weight
    NULLIFY(MeanField)
  
  END DO


CONTAINS


  FUNCTION ParameterName( BaseName, No ) RESULT(str)
!------------------------------------------------------------------------------
    INTEGER :: No
    CHARACTER(LEN=*) :: BaseName
!------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    INTEGER :: i
!------------------------------------------------------------------------------

    WRITE (str,'(A,I0)') TRIM(BaseName)//' ',NoVar

!------------------------------------------------------------------------------
  END FUNCTION ParameterName
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
END SUBROUTINE FilterTimeSeries
!------------------------------------------------------------------------------

