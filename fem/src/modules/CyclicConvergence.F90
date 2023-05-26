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
      meanVal, meanDev, errDev, errVal, errRange, tol, &
      minv, maxv, prevminv=0.0_dp, prevmaxv=0.0_dp, rangev, prevrangev=0.0_dp
  TYPE(Variable_t), POINTER :: perCycleVar, pVar
  INTEGER :: nCycle, iCycle=0, convCycle, minCycle, prevCycle, &
      ProdCycles, n, me, nConv, nCnt, tsize = 500
  REAL(KIND=dp), POINTER :: valTable(:), tmpTable(:)
  LOGICAL :: Visited = .FALSE., Converged = .FALSE.
  CHARACTER(MAX_NAME_LEN) :: valstr
  CHARACTER(*), PARAMETER :: Caller = 'CyclicConvergence'
  
  SAVE Visited, valTable, tsize, prevCycle, prevVal, prevDev, &
      prevminv, prevmaxv, prevrangev, iCycle, ConvCycle, Converged
    
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
    CALL Info(Caller,'Increasing size of valTable to '//I2S(2*tsize),Level=6)
    ALLOCATE(tmpTable(2*tsize))
    tmpTable(1:tsize) = valTable
    tmpTable(tsize+1:2*tsize) = 0.0_dp
    tsize = 2*tsize
    DEALLOCATE(valTable)
    valTable => tmpTable
    NULLIFY(tmpTable)
  END IF
  
  valTable(iCycle) = val

  WRITE(Message,'(A,ES12.3)') 'Entry '//I2S(iCycle)//' for convergence criterion: ',val
  CALL Info(Caller,Message,Level=6)
  
  ! Same cycle, nothing to do!
  IF( nCycle == prevCycle ) RETURN

  ! Compute mean errors over the cycle  
  n = iCycle

  ! We assume here that thet value that we follow has been communicate in parallel as intended so there
  ! is no need to communicate things again. 
  meanVal = SUM(valTable(1:n)) / n
  meanDev = SQRT(SUM((valTable(1:n)-meanVal)**2) / n) 

  minV = MINVAL(valTable(1:n))
  maxV = MAXVAL(valTable(1:n))
  rangeV = maxv - minv
  
  errVal = 2*ABS(meanVal-prevVal)/(ABS(meanVal)+ABS(prevVal))
  errDev = 2*ABS(meanDev-prevDev)/(ABS(meanDev)+ABS(prevDev))  
  errRange = 2*ABS(rangev-prevrangev)/(ABS(meanval)+ABS(prevval))    

  ! Take the smallest entry participating here and let it write the file 
  ! This also ensures that we sync different processes. 
  me = ParEnv % Mype
  me = ParallelReduction(me,1)
  IF(me == ParEnv % MyPe ) CALL WriteConvergenceData()
  
  prevVal = meanVal
  prevDev = meanDev
  prevMinv = minV
  prevMaxv = maxV
  prevRangev = rangeV
  
  minCycle = ListGetInteger( Params,'Cyclic System Min Iterations',Found )
  IF( nCycle < minCycle ) THEN
    CALL Info(Caller,'Number of cycles smaller than required minimum',Level=8)
  ELSE IF(.NOT. Converged) THEN  
    ! All the given measure must converge.
    nConv = 0
    nCnt = 0
    tol = ListGetCReal( Params,'Mean Value Tolerance',Found )
    IF(Found) THEN
      nCnt = nCnt+1
      IF( errVal < tol ) nConv = nConv+1
    END IF
    tol = ListGetCReal( Params,'Mean Deviation Tolerance',Found )
    IF(Found) THEN
      nCnt = nCnt+1
      IF( errDev < tol ) nConv = nConv+1
    END IF
    tol = ListGetCReal( Params,'Range Tolerance',Found )
    IF(Found) THEN
      nCnt = nCnt+1
      IF( errRange < tol ) nConv = nConv+1
    END IF

    IF( nCnt == 0 ) THEN
      CALL Fatal(Caller,'No tolerances given to evaluate convergence!')
    END IF
    
    CALL Info(Caller,I2S(nConv)//' measures converged out of '//I2S(nCnt),Level=10)
    Converged = ( nConv == nCnt ) 

    IF( Converged ) THEN      
      CALL Info(Caller,'Cyclic convergence reached after '//I2S(nCycle)//' cycles!',Level=5)
      pVar => VariableGet( Mesh % Variables,'produce')
      pVar % Values(1) = 1.0_dp
      convCycle = nCycle
    ELSE
      CALL Info(Caller,'Cyclic convergence not reached after '//I2S(nCycle)//' cycles!',Level=10)    
    END IF
  ELSE
    prodCycles = ListGetInteger( Params,'Number of Production Cycles',Found )
    IF(.NOT. Found) prodCycles = 1
    IF( nCycle - convCycle == prodCycles ) THEN
      CALL Info(Caller,'Finish flag was activated after '//I2S(nCycle)//' cycles!',Level=5)
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
    
    REAL(KIND=dp) :: ConvData(7)
    INTEGER :: ConvUnit
    CHARACTER(LEN=MAX_NAME_LEN) :: ConvFile
    LOGICAL, SAVE :: ConvVisited = .FALSE.
    
    !IF( ParEnv % MyPe /= 0 ) RETURN

    ConvFile = ListGetString( Params,'Filename',Found)
    IF(.NOT. Found) ConvFile = 'convergence.dat'

    IF( ConvVisited ) THEN
      OPEN(NEWUNIT=ConvUnit, FILE=ConvFile,STATUS='old',POSITION='append')
    ELSE
      OPEN(NEWUNIT=ConvUnit, File=ConvFile)
      ! Write info in the exact same width as the results
      WRITE(ConvUnit,'(A1,A14,6A15)') '!','time','nCycle','meanVal','meanDev','errVal','errDev','errRange'
      ConvVisited = .TRUE.
    END IF
    
    pVar => VariableGet( Solver % Mesh % Variables, 'time' )
    ConvData(1) = pVar % Values(1) 
    ConvData(2) = perCycleVar % Values(1) 
    ConvData(3) = meanVal
    ConvData(4) = meanDev
    ConvData(5) = errVal
    ConvData(6) = errDev
    ConvData(7) = errRange
    
    WRITE(ConvUnit,'(7ES15.6)') ConvData
    CLOSE(ConvUnit)
    
  END SUBROUTINE WriteConvergenceData       
  
END SUBROUTINE CyclicConvergence
