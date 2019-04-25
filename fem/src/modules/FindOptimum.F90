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
! *  Original Date: 26 Mar 2003
! *
! ****************************************************************************/
 

!------------------------------------------------------------------------------
!>  This solver may be used for optimization. It is not intended to provide
!>  a full featured solution for optimization problems. However, sometimes
!>  its nice that the optimization may be performed in one sweep without the
!>  need for restarting the solver for each optimization trial. 
!> This is a dynamically loaded solver with a standard interface.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE FindOptimum( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE Types
  USE Lists
  USE MeshUtils
  USE Integration
  USE ElementDescription
  USE SolverUtils
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: gotIt, GotIt2, GotInit, OptimalFinish, OptimalStart, InternalHistory, SaveHistory
  LOGICAL, ALLOCATABLE :: FixedParam(:)
  INTEGER :: i,j,k,l,NoParam, NoValues, NoFreeParam, NoOpt, &
      OptimizationsDone=0, Direction=1, NoImprovements=0
  REAL(KIND=dp), POINTER :: Param(:),FoundBetter(:)
  REAL(KIND=dp), ALLOCATABLE :: MinParam(:), MaxParam(:), dParam(:), PrevParam(:,:), &
      PrevCost(:), BestParam(:)
  REAL(KIND=dp) :: Cost, CostTarget, MinCost, x(10), c(10), minv, maxv, OptTol
  CHARACTER(LEN=MAX_NAME_LEN) :: Name, ParamStr, Method
  CHARACTER(LEN=MAX_NAME_LEN) :: BestFile, HistoryFile
  TYPE(Variable_t),POINTER :: Var
  INTEGER :: IOUnit
  
  SAVE Param, MinParam, MaxParam, PrevParam, NoParam, &
      OptimizationsDone, Method, Direction, x, c, PrevCost, &
      FixedParam, NoFreeParam, MinCost, BestParam, NoValues, &
      OptimalFinish, OptimalStart, InternalHistory, dParam, &
      NoImprovements, FoundBetter, OptTol

  !------------------------------------------------------------------------------
  ! In the 1st round perform initializations 
  !------------------------------------------------------------------------------
  IF( OptimizationsDone == 0) THEN

    CALL Info('FindOptimum','Initializing solver for optimization')

    CALL DefaultVariableAdd('Found Better',Global=.TRUE.,InitValue = -1.0_dp)
    Var => DefaultVariableGet('Found Better')
    FoundBetter => Var % Values

    NoParam = Solver % Variable % DOFs
    IF(NoParam == 0) THEN
      CALL Fatal('FindOptimum','There are no parameters to optimize!')
    END IF

    Param => Solver % Variable % Values
    IF( SIZE( Param) /= NoParam ) THEN
      CALL Fatal('FindOptimum','The variable vector is of wrong size')
    END IF

    OptimalFinish = ListGetLogical( Solver % Values,'Optimal Finish',GotIt)
    NoValues = ListGetInteger(Model % Simulation,'Timestep Intervals')

    OptTol = ListGetConstReal( Solver % Values,'Optimization Tolerance',GotIt)

    ALLOCATE( MinParam(NoParam), BestParam(NoParam), MaxParam(NoParam), &
        dParam(NoParam), FixedParam(NoParam))

    MinParam = -HUGE( MinParam ) 
    BestParam = 0.0_dp
    MaxParam = HUGE( MaxParam ) 
    dParam = 0.0_dp
    FixedParam = .FALSE.
    MinCost = HUGE(MinCost)
   
    NoFreeParam = 0
    DO i=1,NoParam
      WRITE( ParamStr,'(A,I0)') 'Parameter ',i

      FixedParam(i) = ListGetLogical(Solver % Values,'Fixed '//TRIM(ParamStr),GotIt)
      Param(i) = ListGetConstReal(Solver % Values,'Initial '//TRIM(ParamStr),GotInit)

      IF(.NOT. ( FixedParam(i) ) ) THEN
        minv = ListGetConstReal(Solver % Values,'Min '//TRIM(ParamStr),GotIt)
        maxv = ListGetConstReal(Solver % Values,'Max '//TRIM(ParamStr),GotIt2)

        IF( GotIt ) MinParam(i) = minv
        IF( GotIt2 ) MaxParam(i) = maxv

        IF( .NOT. GotInit ) THEN
          IF( GotIt .AND. GotIt2 ) THEN
            Param(i) = 0.5_dp * ( minv + maxv ) 
          ELSE IF( GotIt ) THEN
            Param(i) = minv 
          ELSE IF( GotIt2 ) THEN
            Param(i) = maxv 
          END IF
        END IF
      END IF
      dParam(i) = ListGetConstReal(Solver % Values,'Scale '//TRIM(ParamStr),GotIt)
      IF(.NOT. GotIt) dParam(i) = 1.0_dp
    END DO
    
    NoFreeParam = NoParam - COUNT(FixedParam)
    IF( NoFreeParam == 0 ) THEN
      CALL Warn('FindOptimum','All parameters are fixed, no optimization!')
      RETURN
    END IF

    Method = ListGetString(Solver % Values,'Optimization Method')

    ! Internal history could be used in more complicated optimization routines
    !--------------------------------------------------------------------------
    InternalHistory = ListGetLogical( Solver % Values,'Internal History',GotIt)    
    IF( Method == 'bisect') InternalHistory = .TRUE.
    IF( InternalHistory ) THEN
      ALLOCATE( PrevParam(NoValues,NoParam), PrevCost(NoValues))
    END IF

    OptimalStart = ListGetLogical(Solver % Values,'Optimal Restart',GotIt)
    IF( OptimalStart ) THEN
      CALL GuessOptimum()
      GOTO 100
    END IF
    
  END IF

  !------------------------------------------------------------------------------
  ! If visiting for the second time inspect how good the previous solution was
  ! and improve an it. 
  !------------------------------------------------------------------------------
  IF( OptimizationsDone > 0 ) THEN
    FoundBetter = -1.0_dp
    Cost = GetCReal(Model % Simulation,'Cost Function',GotIt)

    IF(.NOT. GotIt) Cost = GetCReal(Solver % Values,'Cost Function',GotIt)

    IF(.NOT. GotIt) THEN
      Name = ListGetString(Solver % Values,'Cost Function Name',GotIt)
      IF(.NOT. GotIt) CALL Fatal('FindOptimum','Give Cost Function or its name')
      Cost = ListGetConstReal(Model % Simulation,Name,GotIt)     
      IF(.NOT. GotIt) CALL Fatal('FindOptimum','Cost with the given name was not found')
    END IF

    ! Whether to perform search rather than optimization. 
    ! In this case reduce the goal so that the target will always be zero.
    !----------------------------------------------------------------------
    CostTarget = ListGetConstReal( Solver % Values,'Cost Function Target',GotIt)    
    IF( GotIt ) Cost = Cost - CostTarget 

    ! The cost function could be the absolute value
    ! or we could transfer a maximization problem into minimization.
    !----------------------------------------------------------------
    IF( ListGetLogical( Solver % Values,'Cost Function Absolute',GotIt)) THEN
      Cost = ABS( Cost )
    ELSE IF( ListGetLogical( Solver % Values,'Cost Function Maximize',GotIt)) THEN
      Cost = -Cost
    END IF


    ! Found a new best point
    !-----------------------
    IF(Cost < MinCost) THEN
      FoundBetter = 1.0_dp
      MinCost = Cost
      BestParam(1:NoParam) = Param(1:NoParam)
      NoImprovements = NoImprovements + 1
      
      WRITE(Message,'(A,ES15.6E3)') 'Found New Minimum Cost:',MinCost
      CALL Info('FindOptimum',Message,Level=4)

      BestFile = ListGetString(Solver % Values,'Best File',GotIt )
      IF( GotIt ) THEN
        OPEN( NEWUNIT=IOUnit, FILE=BestFile, STATUS='UNKNOWN')
        WRITE (IOUnit,'(I0)') NoParam
        DO i=1,NoParam
          WRITE (IOUnit,'(ES17.8E3)') Param(i)
        END DO
        WRITE (IOUnit,'(ES17.8E3)') Cost
        WRITE (IOUnit,'(I0)') NoImprovements
        WRITE (IOUnit,'(I0)') OptimizationsDone 
        CLOSE(IOUnit)
      END IF
    END IF

    IF( InternalHistory ) THEN
      PrevParam(OptimizationsDone,1:NoParam) = Param(1:NoParam)
      PrevCost(OptimizationsDone) = Cost
    END IF

    ! Save the results to a file
    !---------------------------
    HistoryFile = ListGetString(Solver % Values,'History File',GotIt )
    IF( GotIt ) THEN
      IF(OptimizationsDone == 1) THEN
        OPEN (NEWUNIT=IOUnit, FILE=HistoryFile)
      ELSE
        OPEN (NEWUNIT=IOUnit, FILE=HistoryFile,POSITION='APPEND')
      END IF
      
      WRITE (IOUnit,'(I8,ES17.8E3)',advance='no') OptimizationsDone, Cost
      DO i=1,NoParam
        WRITE (IOUnit,'(ES17.8E3)',advance='no') Param(i)
      END DO
      WRITE(IOUnit,'(A)') ' '
      CLOSE(IOUnit)
    END IF
    
    IF( OptimalFinish .AND. OptimizationsDone == NoValues - 1 ) THEN
      CALL Info('FindOptimum','Peforming the last step with the best so far')
      Param = BestParam
      GOTO 100
    END IF
  END IF


  CALL Info( 'FindOptimum', '-----------------------------------------', Level=4 )
  WRITE( Message, '(A,I0,A,A)' ) 'Manipulating ',NoFreeParam,' parameters using ',TRIM(Method) 
  CALL Info( 'FindOptimum', Message, Level=4 )
  IF(OptimizationsDone > 0) THEN
    WRITE( Message, '(A,ES15.6E3)' ) 'Last evaluated cost: ',Cost    
    CALL Info('FindOptimim',Message,Level=5)
  END IF
  WRITE( Message, '(A,ES15.6E3)' ) 'Lowest cost so far: ',MinCost
  CALL Info('FindOptimum',Message,Level=5)


  SELECT CASE(Method)
    
  CASE ('random')
    CALL RandomParameter()
    
  CASE ('scanning')
    CALL ScanParameter()

  CASE ('secant')
    CALL SecantSearch()

  CASE ('genetic')
    CALL GeneticOptimize(NoParam, Param, Cost)
    
  CASE ('bisect')    
    CALL BisectOptimize()

  CASE ('simplex')    
    CALL SimplexOptimize( NoParam, Param, Cost, MinParam, MaxParam, dParam )
    
  CASE DEFAULT
    CALL Fatal('FindOptimum','Unknown method')
    
  END SELECT
 
  IF(.FALSE.) THEN
    DO i=1,NoParam 
      IF( FixedParam(i) ) CYCLE
      Param(i) = MAX(MinParam(i),Param(i))
      Param(i) = MIN(MaxParam(i),Param(i))
    END DO
  END IF


  IF( NoParam <= 5 ) THEN
    WRITE( Message, * ) 'Next set: ',Param
    CALL Info( 'FindOptimum', Message, Level=4 )
  END IF
  CALL Info( 'FindOptimum', '-----------------------------------------', Level=4 )


100 OptimizationsDone = OptimizationsDone + 1

CONTAINS

!-------------------------------------------------------------------------------

  FUNCTION rnd(n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dp), DIMENSION(n) :: rnd
    CALL RANDOM_NUMBER(rnd)
  END FUNCTION rnd

!-------------------------------------------------------------------------------

  INTEGER FUNCTION idx(n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dp) :: x
    CALL RANDOM_NUMBER(x)
    idx = n*x + 1
  END FUNCTION idx

!-------------------------------------------------------------------------------
!> Choose next parameter set from genetic optimization procedure
!-------------------------------------------------------------------------------

  SUBROUTINE GeneticOptimize(parsize, parameters, func)

    INTEGER :: parsize, no = 0
    REAL (KIND=dp) :: parameters(parsize), func

    INTEGER :: popsize, i0, i1, i2, i3 
    REAL(KIND=dp) :: popcoeff, popcross
    REAL(KIND=dp), ALLOCATABLE :: pars(:,:), vals(:) ,rnds(:)
    LOGICAL, ALLOCATABLE :: mask(:)

    SAVE no, i0, pars, vals, rnds, mask, popsize, popcoeff, popcross
        

    no = no + 1

    IF(no == 1) THEN
      popsize = ListGetInteger(Solver % Values,'Population Size',GotIt)
      IF(.NOT. GotIt) popsize = 5 * parsize
      popcoeff = ListGetConstReal(Solver % Values,'Population Coefficient',GotIt)
      IF(.NOT. GotIt) popcoeff = 0.7
      popcross = ListGetConstReal(Solver % Values,'Population Crossover',GotIt)
      IF(.NOT. GotIt) popcross = 0.1
      ALLOCATE(pars(parsize,popsize),vals(popsize),mask(parsize),rnds(parsize))
      IF(.FALSE.) THEN
        PRINT *,'popsize',popsize,'parsize',parsize
        PRINT *,'popcoeff',popcoeff,'popcross',popcross
      END IF
    END IF
    
    ! Read the cases into the population
    IF(no <= popsize) THEN
      pars(1:parsize,no) = parameters(1:parsize)
      vals(no) = func
    ELSE   
      IF(func < vals(i0)) THEN
        pars(1:parsize,i0) = parameters(1:parsize) 
        vals(i0) = func
      END IF
    END IF

    ! The first cases are just random
    IF(no < popsize) THEN
      pars(1:parsize,no) = parameters(1:parsize)
      vals(no) = func
      Param = MinParam + (MaxParam-MinParam) * rnd(parsize)
    END IF

    ! Here use genetic algorithms 
    IF(no >= popsize) THEN
      ! Find the three vectors to recombine 
      i0 = MOD(no,popsize) + 1 
      DO
        i1 = idx(popsize)
        IF (i1 /= i0) EXIT
      END DO
      DO
        i2 = idx(popsize)
        IF (i2 /= i0.AND. i2 /= i1) EXIT
      END DO
      DO
        i3 = idx(popsize)
        IF (ALL(i3 /= (/i0,i1,i2/))) EXIT
      END DO
      
      rnds = rnd(parsize)
      mask = (rnds < popcross)
      
      WHERE (mask)
        parameters = pars(:,i3) + popcoeff*(pars(:,i1)-pars(:,i2))
      ELSEWHERE
        parameters = pars(:,i0)
      END WHERE

      parameters = MAX( parameters, MinParam ) 
      parameters = MIN( parameters, MaxParam ) 

    END IF

  END SUBROUTINE GeneticOptimize

!-------------------------------------------------------------------------------
!> Choose next parameter set from even random distribution
!-------------------------------------------------------------------------------

  SUBROUTINE RandomParameter()

    INTEGER :: i
    REAL(KIND=dp) :: Extent

    DO i=1,NoParam
      CALL RANDOM_NUMBER(Extent)
      Param(i) = MinParam(i) + (MaxParam(i)-MinParam(i)) * Extent
    END DO

  END SUBROUTINE RandomParameter

!-------------------------------------------------------------------------------
!> Choose next parameter from 1D parameter scanning
!-------------------------------------------------------------------------------

  SUBROUTINE ScanParameter()

    INTEGER :: no, maxno, i,j
    REAL(KIND=dp) :: Extent

    SAVE no, i, maxno

    IF( no == 0 ) THEN
      IF(NoFreeParam /= 1) CALL Fatal('FindOptimum',&
          'Option scan implemented only for one parameter')
      DO i=1,NoParam
        IF(.NOT. FixedParam(i)) EXIT
      END DO
      WRITE(Message,'(A,I0)') 'Applying scanning to parameter ',i
      CALL Info('FindOptimum',Message)

      maxno = NoValues 
      IF( OptimalFinish ) maxno = maxno -1
      IF( OptimalStart ) maxno = maxno - 1
    END IF

    Extent = no * 1.0_dp/(maxno-1)
    Param(i) = MinParam(i) + (MaxParam(i)-MinParam(i)) * Extent
    no = no + 1

  END SUBROUTINE ScanParameter

!-------------------------------------------------------------------------------
!> Choose next parameter set from 1D bisection search
!-------------------------------------------------------------------------------

  SUBROUTINE BisectOptimize()

    INTEGER :: j, no = 0
    REAL(KIND=dp) :: step 

    SAVE j, no, step

    IF(NoFreeParam /= 1) CALL Fatal('FindOptimum',&
        'Option bisect implemented only for one parameter')
    IF( no == 0 ) THEN
      DO j=1,NoParam
        IF(.NOT. FixedParam(j)) EXIT
      END DO
      WRITE(Message,'(A,I0)') 'Applying bisection search to parameter ',j
      CALL Info('FindOptimum',Message)
    END IF
    
    no = no + 1

    IF(no == 1) THEN
      step = ListGetConstReal(Solver % Values,'Step Size',GotIt)
      IF(.NOT. GotIt) step = (MaxParam(j)-Param(j))/2.0
      step = MIN((MaxParam(j)-Param(j))/2.0,step)
    END IF
    
    IF(no <= 3) THEN
      Param(j) = Param(j) + step
      RETURN
    END IF

    IF(no == 4) THEN
      x(1) = PrevParam(1,j)
      x(2) = PrevParam(2,j)
      x(3) = PrevParam(3,j)
      c(1) = PrevCost(1)
      c(2) = PrevCost(2)
      c(3) = PrevCost(3)
    ELSE
      x(3) = Param(j)
      c(3) = Cost
    END IF

    ! Order the previous points so that x1 < x2 < x3
    DO k=1,2 
      DO i=k+1,3
        IF(x(i) < x(k)) THEN
          x(4) = x(k)
          x(k) = x(i)
          x(i) = x(4)
          c(4) = c(k)
          c(k) = c(i)
          c(i) = c(4)
        END IF
      END DO
    END DO
    
    ! Monotonic line segment
    IF( (c(2)-c(1))*(c(3)-c(2)) > 0.0) THEN
      IF(c(3) < c(1)) THEN
        Param(j) = x(3) + SIGN(step,x(3)-x(1))
        c(1) = c(3)
        x(1) = x(3)
      ELSE
        Param(j) = x(1) + SIGN(step,x(1)-x(3))
      END IF
    ELSE IF(c(2) < c(1) .OR. c(2) < c(3)) THEN 
      IF(c(3) < c(1)) THEN
        c(1) = c(3)
        x(1) = x(3)
      END IF
      step = (x(2)-x(1))/2.0d0
      Param(j) = x(1) + SIGN(step,x(2)-x(1))
    ELSE
      CALL Fatal('FindOptimum','Bisection method cannot handle local maxima')
    END IF

  END SUBROUTINE BisectOptimize

!-------------------------------------------------------------------------------
!> Choose next parameter set from secant method
!> This only works for design problems where the target cost is known.
!-------------------------------------------------------------------------------

  SUBROUTINE SecantSearch()

    INTEGER :: j, no = 0
    REAL(KIND=dp) :: step, maxstep, relax
    REAL(KIND=dp) :: x0=0.0,x1=0.0,x2=0.0,f0=0.0,f1=0.0,dx

    SAVE j, no, step, maxstep, relax, x0,x1,x2,f0,f1

    IF(NoFreeParam /= 1) CALL Fatal('FindOptimum',&
        'Secant search implemented only for one parameter')
    IF( no == 0 ) THEN
      DO j=1,NoParam
        IF(.NOT. FixedParam(j)) EXIT
      END DO
      WRITE(Message,'(A,I0)') 'Applying secant search to parameter ',j
      CALL Info('FindOptimum',Message)
    END IF
    
    no = no + 1

    IF(no == 1) THEN
      maxstep = ListGetConstReal(Solver % Values,'Max Step Size',GotIt)
      step = ListGetConstReal(Solver % Values,'Step Size',GotIt)
      IF(.NOT. GotIt) step = 1.0d-3*(MaxParam(j)-MinParam(j))
      relax = ListGetConstReal(Solver % Values,'Relaxation Factor',GotIt)
      IF(.NOT. GotIt) Relax = 1.0_dp
    END IF

    x0 = x1
    x1 = x2
    f0 = f1 
    f1 = Cost
       
    IF(no <= 2) THEN
      x2 = Param(j) + (no-1)*step
    ELSE IF( ABS(f1) < OptTol ) THEN
      CALL Info('SecantSearch','Tolerance reached, doing nothing')
      x2 = x1
    ELSE
      dx = relax * f1 * (x1-x0) / (f1-f0)      
      IF( ABS( dx ) > maxstep ) THEN
        dx = SIGN( maxstep, dx )
      END IF	
      x2 = x1 - dx 
    END IF

    Param(j) = x2

  END SUBROUTINE SecantSearch


!--------------------------------------------------------------------------------------
!> Find the optimum using the Simplex method (Nelder-Mead algorithm)
!> Note that constraint box is taken into account only when creating the initial simplex.
!--------------------------------------------------------------------------------------
  SUBROUTINE SimplexOptimize( nx, x, cost, minx, maxx, diffx )

    INTEGER :: nx
    REAL (KIND=dp) :: cost,x(:),minx(:),maxx(:),diffx(:)
!------------------------------------------------------------------------------
    INTEGER :: i,j,il,ih,is,no = 0, nomax, mode = 0, submode = 0
    LOGICAL :: Found, AllocationsDone
    REAL(KIND=dp) :: lambda, fl,fh,fs,fr,fe,ratio=0.0, maxratio
    REAL(KIND=dp) :: alpha=1.0_dp,beta=0.5_dp,gamma=2.0_dp,delta=0.5_dp
    REAL(KIND=dp), POINTER :: eig(:,:), ls(:), xall(:,:), f(:), x0(:), xr(:),xc(:)

    SAVE no,eig,ls,xall,f,x0,xr,xc,mode,submode,ih,il,is,fr,fs,fe,fh,fl,&
        nomax,ratio,maxratio,AllocationsDone

    no = no + 1

!    PRINT *,'Simplex: ',no,mode,maxratio,Cost,ratio
    

    ! Initialize the unit vectors
    !-----------------------------
    IF(.NOT. AllocationsDone ) THEN
      ALLOCATE(eig(nx,nx),ls(nx),xall(nx+1,nx),f(nx+1),x0(nx),xr(nx),xc(nx))

      lambda = ListGetConstReal(Solver % Values,'Simplex Relative Length Scale',Found)
      IF(.NOT. Found ) lambda = 0.01_dp
      DO i=1,nx
        IF( ABS( diffx(i) ) > TINY(diffx(i)) ) THEN
          ls(i) = lambda * diffx(i)
        ELSE
          IF( maxx(i) - x(i) > x(i) - minx(i) ) THEN
            ls(i) = lambda * (maxx(i) - x(i))
          ELSE
            ls(i) = lambda * (minx(i) - x(i)) 
          END IF
        END IF
      END DO

      eig = 0.0
      DO i=1,nx
        eig(i,i) = 1.0_dp
      END DO

      nomax = ListGetInteger(Solver % Values,'Simplex Restart Interval',Found)
      maxratio = ListGetConstReal(Solver % Values,'Simplex Restart Convergence Ratio',Found)      
      IF(.NOT. Found) maxratio = 1.0_dp
      AllocationsDone = .TRUE.
    END IF

    IF( nomax > 0 .AND. no > nomax .AND. ratio > maxratio ) THEN
      CALL Info('FindOptimum','Making a restart in simplex')
      ratio = 0.0_dp
      ls = 0.1 * ls

      PRINT *,'Simplex: coeff',ls
      no = 1
    END IF

    ! Memorize the 1st simplex
    !--------------------------
    IF( no > 1 .AND. no <= nx + 2 ) THEN
      xall(no-1,:) = x
      f(no-1) = cost
    END IF


    ! Create the 1st simplex
    !---------------------------
    IF( no <= nx + 1) THEN
      IF( no == 1 ) THEN
        x0 = x
      ELSE
        x = x0 + ls(no-1) * eig(no-1,:)
      END IF
      RETURN
    END IF

    IF( mode <= 1 ) THEN
      ! Find the minimum and maximum nodes in the simplex
      !--------------------------------------------------
      
      fl = HUGE(fl)  ! best
      fh = -HUGE(fh) ! worst   
      fs = -HUGE(fs) ! second worst   

      DO i=1,nx+1
        IF( f(i) < fl ) THEN
          il = i
          fl = f(i)
        END IF
        IF( f(i) > fh ) THEN
          ih = i
          fh = f(i)
        END IF
      END DO
      DO i=1,nx+1
        IF(i == ih) CYCLE
        IF(i == il) CYCLE
        IF( f(i) > fs ) THEN
          is = i
          fs = f(i)
        END IF
      END DO

      ! Minpoint neglecting the worst point
      !------------------------------------
      xc = 0.0_dp
      DO i=1,nx+1
        IF( i == ih ) CYCLE
        xc = xc + xall(i,:) 
      END DO
      xc = xc / nx

    END IF


    ! mode
    ! 0 - undecided
    ! 1 - reflection
    ! 2 - expansion
    ! 3 - contraction
    ! 4 - shrinkage
    ! 5 - finished

    Found = .FALSE.
    IF( mode == 0 ) THEN
      mode = 1
    ELSE IF( mode == 1) THEN
      fr = cost
      xr = x
      IF( fr <= fs .AND. fr >= fl ) THEN
        Found = .TRUE.
      ELSE IF( fr < fl ) THEN
        mode = 2
      ELSE IF(fr > fs) THEN
        mode = 3
        IF( fr < fh ) THEN
          submode = 1
        ELSE
          submode = 2
        END IF
      END IF

    ELSE IF(mode == 2) THEN
      fe = cost
      IF( fe >= fr ) THEN
        x = xr
        cost = fr
      END IF
      Found = .TRUE.
    ELSE IF(mode == 3) THEN
      IF( submode == 1 .AND. cost < fr ) THEN
        Found = .TRUE.
      ELSE IF( submode == 2 .AND. cost < fh ) THEN
        Found = .TRUE.
      ELSE
        mode = 4
        submode = 0
      END IF      
    ELSE IF( mode == 4) THEN
      CALL Warn('FindOptimum','Srink mode not ok yet!')
      IF( submode == nx ) THEN
        mode = 1
        submode = 0
      END IF
    END IF


    IF( Found ) THEN
      xall(ih,:) = x 

      ratio = cost / f(ih)

      f(ih) = cost
      mode = 1
      submode = 0

      ! Find the minimum and maximum nodes in the simplex
      !--------------------------------------------------
      
      fl = HUGE(fl)  ! best
      fh = -HUGE(fh) ! worst   
      fs = -HUGE(fs) ! second worst   
      DO i=1,nx+1
        IF( f(i) < fl ) THEN
          il = i
          fl = f(i)
        END IF
        IF( f(i) > fh ) THEN
          ih = i
          fh = f(i)
        END IF
      END DO
      DO i=1,nx+1
        IF(i == ih) CYCLE
        IF( f(i) > fs ) THEN
          is = i
          fs = f(i)
        END IF
      END DO

      ! Minpoint neglecting the worst point
      !------------------------------------
      xc = 0.0_dp
      DO i=1,nx+1
        IF( i == ih ) CYCLE
        xc = xc + xall(i,:) 
      END DO
      xc = xc / nx
    END IF


    IF( mode == 1 ) THEN
      x = xc + alpha*(xc - xall(ih,:))
    ELSE IF( mode == 2 ) THEN
      x = xc + gamma*(xc - xall(ih,:))
    ELSE IF(mode == 3 ) THEN
      IF( submode == 1 ) THEN
        x = xc + beta*(xr - xc)        
      ELSE
        x = xc + beta*(xall(ih,:) - xc)        
      END IF
    ELSE IF(mode == 4) THEN
      submode = submode + 1 
      i = submode 
      IF( submode >= il ) i = i + 1
      x = xall(il,:) + delta*(xall(i,:)-xall(il,:))
      xall(i,:) = x
    END IF

  END SUBROUTINE SimplexOptimize



!--------------------------------------------------------------------------------------
!> This subroutine may be used to continue the optimization from the previous best value.
!--------------------------------------------------------------------------------------
  SUBROUTINE GuessOptimum( )
!------------------------------------------------------------------------------

    INTEGER :: i,n
    REAL(KIND=dp) :: parami
    REAL(KIND=dp), ALLOCATABLE :: guessparam(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: Name
    LOGICAL :: fileis, GotIt
    CHARACTER(LEN=MAX_NAME_LEN) :: GuessFile
    INTEGER :: IOUnit
    
    GuessFile = ListGetString(Solver % Values,'Guess File',GotIt )
    IF(.NOT. GotIt) GuessFile = 'best.dat'

    INQUIRE (FILE=GuessFile, EXIST=fileis)
    
    IF(.NOT. fileis ) THEN
      CALL Warn('FindOptimum','Previous optimum was not found in: '//TRIM(GuessFile))
      RETURN
    END IF

    OPEN(NEWUNIT=IOUnit,FILE=GuessFile)
    READ (IOUnit,*) n
    ALLOCATE (guessparam(n))
    DO i=1,n
      READ (IOUnit,*) guessparam(i)
    END DO
    CLOSE(IOUnit)

    n = MIN( n, SIZE( param) )
    param(1:n) = guessparam(1:n)

  END SUBROUTINE GuessOptimum


END SUBROUTINE FindOptimum
