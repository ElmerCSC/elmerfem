!/*****************************************************************************
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
! *  A utility to decide whether a cyclic simulation has converged or not.
! *
! ******************************************************************************
! *
! *  Authors: Peter Råback
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 16.9.2020
! *
! *****************************************************************************/

SUBROUTINE CyclicConvergence_init( Model,Solver,dt,Transient)

  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  CALL ListAddNewString( Solver % Values,'Equation','CyclicConvergence')
  
  ! Introduce variable so that we have a place for pseudonorm
  CALL ListAddNewString( Solver % Values,'Variable',&
      '-nooutput -global convcrit_var') 

END SUBROUTINE CyclicConvergence_Init


SUBROUTINE CyclicConvergence( Model,Solver,dt,Transient)

  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!---------------------------------------------------------------
  LOGICAL :: Found
  TYPE(Element_t), POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  REAL(KIND=dp) :: val, prevVal=0.0_dp, prevDev=0.0_dp, &
      meanVal, meanDev, errDev, errVal, tolDev, tolVal
  TYPE(Variable_t), POINTER :: perCycleVar, pVar
  INTEGER :: nCycle, iCycle=0, convCycle, minCycle, prevCycle, &
      ProdCycles, n, tsize = 500
  REAL(KIND=dp), POINTER :: valTable(:), tmpTable(:)
  LOGICAL :: Visited = .FALSE., Converged = .FALSE.
  CHARACTER(MAX_NAME_LEN) :: valstr
  CHARACTER(*), PARAMETER :: Caller = 'CyclicConvergence'
  
  SAVE Visited, valTable, tsize, prevCycle, prevVal, prevDev, iCycle, &
      ConvCycle, Converged
    
  CALL Info(Caller,'Checking for convergence of cyclic system')
  
  Mesh => Solver % Mesh
  Params => Solver % Values

  perCycleVar => VariableGet( Mesh % Variables,'Periodic Cycle')
  IF(.NOT. ASSOCIATED( perCycleVar ) ) THEN
    CALL Fatal(Caller,'Could not define "Periodic Cycle"')     
  END IF
  ! We apply a small offset so that at the exact time we would
  ! swap the cycle is triggered just at the end of cycle.   
  nCycle = FLOOR( perCycleVar % Values(1) + 1.0e-3*dt )

  IF(.NOT. Visited) THEN
    prevCycle = nCycle
    ALLOCATE( valTable(tsize) )
    valTable = 0.0_dp
    Visited = .TRUE.
  END IF
    
  val = ListGetCReal( Params,'Convergence Value',Found )
  IF(.NOT. Found ) THEN
    valstr = ListGetString( Params,'Convergence Value Name',Found )
    IF(.NOT. Found) THEN
      CALL Fatal(Caller,'Give "Convergence Value" or "Convergence Value Name"!')
    END IF    
    val = ListGetCReal( Model % Simulation, valstr, Found )
    IF(.NOT. Found) THEN
      pVar => VariableGet( Mesh % Variables, valstr )
      IF(ASSOCIATED( pVar ) ) THEN
        val = pVar % Values(1)
      ELSE
        CALL Fatal(Caller,'Could not find "Convergence Value Name" as variable or list entry!') 
      END IF
    END IF
  END IF
   
  ! Tabulate values
  iCycle = iCycle + 1

  IF( iCycle > tsize ) THEN
    CALL Info(Caller,'Increasing size of valTable to '//TRIM(I2S(2*tsize)),Level=6)
    ALLOCATE(tmpTable(2*tsize))
    tmpTable(1:tsize) = valTable
    tmpTable(tsize+1:2*tsize) = 0.0_dp
    tsize = 2*tsize
    DEALLOCATE(valTable)
    valTable => tmpTable
    NULLIFY(tmpTable)
  END IF
  
  valTable(iCycle) = val

  WRITE(Message,'(A,ES12.3)') 'Entry '//TRIM(I2S(iCycle))//' for convergence criterion: ',val
  CALL Info(Caller,Message,Level=6)
  
  ! Same cycle, nothing to do!
  IF( nCycle == prevCycle ) RETURN

  ! Compute mean errors over the cycle  
  n = iCycle

  ! We assume here that thet value that we follow has been communicate in parallel as intended so there
  ! is no need to communicate things again. 
  meanVal = SUM(valTable(1:n)) / n
  meanDev = SQRT(SUM((valTable(1:n)-meanVal)**2) / n) 
  
  errVal = 2*ABS(meanVal-prevVal)/(ABS(meanVal)+ABS(prevVal))
  errDev = 2*ABS(meanDev-prevDev)/(ABS(meanDev)+ABS(prevDev))  

  ! Take the max norm in parallel to be on the safe side.  
  ! Otherwise some partition might start production earlier than others. 
  errVal = ParallelReduction(errVal,2)
  errDev = ParallelReduction(errDev,2)
  
  CALL WriteConvergenceData()
  
  prevVal = meanVal
  prevDev = meanDev   

  minCycle = ListGetInteger( Params,'Cyclic System Min Iterations',Found )
  IF( nCycle < minCycle ) THEN
    CALL Info(Caller,'Number of cycles smaller than required minimum',Level=8)
  ELSE IF(.NOT. Converged) THEN  
    tolVal = ListGetCReal( Params,'Mean Value Tolerance',Found )
    IF(Found) Converged = (errVal < tolVal) 
    tolDev = ListGetCReal( Params,'Mean Deviation Tolerance',Found )
    IF(Found) Converged = (errDev < tolDev )
    IF( Converged ) THEN      
      CALL Info(Caller,'Cyclic convergence reached after '//TRIM(I2S(nCycle))//' cycles!',Level=5)
      pVar => VariableGet( Mesh % Variables,'produce')
      pVar % Values(1) = 1.0_dp
      convCycle = nCycle
    ELSE
      CALL Info(Caller,'Cyclic convergence not reached after '//TRIM(I2S(nCycle))//' cycles!',Level=10)    
    END IF
  ELSE
    prodCycles = ListGetInteger( Params,'Number of Production Cycles',Found )
    IF(.NOT. Found) prodCycles = 1
    IF( nCycle - convCycle == prodCycles ) THEN
      CALL Info(Caller,'Finish flag was activated after '//TRIM(I2S(nCycle))//' cycles!',Level=5)
      pVar => VariableGet( Mesh % Variables,'finish')
      pVar % Values(1) = 1.0_dp
    END IF
  END IF

  !PRINT *,'Convergence:',nCycle,icycle,val,meanVal,errVal
  
  ! Restart the counter
  iCycle = 0
  prevCycle = nCycle
  
CONTAINS

  SUBROUTINE WriteConvergenceData()
    
    REAL(KIND=dp) :: ConvData(6)
    INTEGER :: ConvUnit
    CHARACTER(LEN=MAX_NAME_LEN) :: ConvFile
    LOGICAL, SAVE :: ConvVisited = .FALSE.
    
    IF( ParEnv % MyPe /= 0 ) RETURN

    ConvFile = ListGetString( Params,'Filename',Found)
    IF(.NOT. Found) ConvFile = 'convergence.dat'

    IF( ConvVisited ) THEN
      OPEN(NEWUNIT=ConvUnit, FILE=ConvFile,STATUS='old',POSITION='append')
    ELSE
      OPEN(NEWUNIT=ConvUnit, File=ConvFile)
      ! Write info in the exact same width as the results
      WRITE(ConvUnit,'(A1,A14,5A15)') '!','time','nCycle','meanVal','meanDev','errVal','errDev'
      ConvVisited = .TRUE.
    END IF
    
    pVar => VariableGet( Solver % Mesh % Variables, 'time' )
    ConvData(1) = pVar % Values(1) 
    ConvData(2) = perCycleVar % Values(1) 
    ConvData(3) = meanVal
    ConvData(4) = meanDev
    ConvData(5) = errVal
    ConvData(6) = errDev
    
    WRITE(ConvUnit,'(6ES15.6)') ConvData
    CLOSE(ConvUnit)
    
  END SUBROUTINE WriteConvergenceData       
  
END SUBROUTINE CyclicConvergence
