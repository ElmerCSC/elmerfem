!> \ingroup ElmerLib
!> \{

!-----------------------------------------------------------------------------
!> Module for various routines related to optimization, control and
!> parametric studies.
!-----------------------------------------------------------------------------
#include "../config.h"

MODULE OptimizationUtils
  
  USE Lists
  USE ModelDescription
  USE Messages

    
  IMPLICIT NONE

  PUBLIC :: GetCostFunction
  PUBLIC :: ControlParameters
  PUBLIC :: SetTabulatedParameters
  PUBLIC :: SetOptimizationParameters
  PUBLIC :: ControlResetMesh
#ifdef HAVE_EXTOPTIM
  PUBLIC :: ExternalOptimization_minpack
  PUBLIC :: ExternalOptimization_newuoa
  PUBLIC :: ExternalOptimization_bobyqa
#endif
  
  PRIVATE :: SaveCurrentOptimum
  PRIVATE :: GetSavedOptimum
  
CONTAINS

  ! Obtains cost function value that has must be given
  ! by a keyword usually calling a user defined function.
  !-------------------------------------------------------
  SUBROUTINE GetCostFunction(OptList,Cost,GotCost)

    TYPE(ValueList_t), POINTER :: OptList
    REAL(KIND=dp) :: Cost
    LOGICAL :: GotCost

    REAL(KIND=dp) :: CostTarget
    LOGICAL :: GotIt
    CHARACTER(:), ALLOCATABLE :: Name

    Cost = ListGetCReal(OptList,'Cost Function',GotCost)

    IF(.NOT. GotCost) THEN
      Name = ListGetString(OptList,'Cost Function Name',GotIt)
      IF(.NOT. GotIt ) THEN
        Cost = ListGetCReal(OptList,Name,GotCost)
      END IF
    END IF

    IF(.NOT. GotCost ) RETURN

    ! Whether to perform search rather than optimization. 
    ! In this case reduce the goal so that the target will always be zero.
    !----------------------------------------------------------------------
    CostTarget = ListGetConstReal( OptList,'Cost Function Target',GotIt)    
    IF( GotIt ) Cost = Cost - CostTarget 

    ! The cost function could be the absolute value
    ! or we could transfer a maximization problem into minimization.
    !----------------------------------------------------------------
    IF( ListGetLogical( OptList,'Cost Function Absolute',GotIt)) THEN
      Cost = ABS( Cost )
    ELSE IF( ListGetLogical( OptList,'Cost Function Maximize',GotIt)) THEN
      Cost = -Cost
    END IF

    WRITE( Message, '(A,ES15.6E3)' ) 'Cost function: ',Cost    
    CALL Info('GetCostFunction',Message,Level=5)

  END SUBROUTINE GetCostFunction


  !> This subroutine may be used to continue the optimization from the previous best value.
  !--------------------------------------------------------------------------------------
  SUBROUTINE GetSavedOptimum(OptList,x,Found)
    !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: OptList
    REAL(KIND=dp) :: x(:)
    LOGICAL :: Found

    INTEGER :: i,n
    CHARACTER(:), ALLOCATABLE :: Name
    LOGICAL :: fileis, GotIt, OptimalStart
    INTEGER :: IOUnit

    Found = .FALSE.

    OptimalStart = ListGetLogical(OptList,'Optimal Restart',GotIt )
    IF(.NOT. GotIt ) OptimalStart = ListCheckPresent( OptList,'Parameter Restart File' )
    IF(.NOT. OptimalStart) RETURN

    Name = ListGetString(OptList,'Parameter Restart File',GotIt )
    IF(.NOT. GotIt) Name = "optimize-best.dat"
    
    CALL Info('GetSavedOptimum','Trying to use previous optimal  parameters: '//TRIM(Name),Level=6)
    INQUIRE (FILE=Name, EXIST=fileis)

    IF(.NOT. fileis ) THEN
      CALL Warn('GetSavedOptimum','Previous optimum was not found in: '//TRIM(Name))
      RETURN
    END IF

    x = 0.0_dp
    
    OPEN(NEWUNIT=IOUnit,FILE=Name)
    READ (IOUnit,*) n
    n = MIN(n,SIZE(x))
    READ (IOUnit,*) x(1:n)
    CLOSE(IOUnit)

    Found = .TRUE.
    
    CALL Info('GetSavedOptimum','Number of parameters initialized from file: '//I2S(n),Level=6)

  END SUBROUTINE GetSavedOptimum



  ! We may save the current optimum so that restarting is easier.
  !----------------------------------------------------------------
  SUBROUTINE SaveCurrentOptimum(OptList,n,rpar,Cost,Iters,Improvements) 
    TYPE(ValueList_t), POINTER :: OptList
    INTEGER :: n
    REAL(KIND=dp) :: rpar(:)
    REAL(KIND=dp), OPTIONAL :: Cost
    INTEGER, OPTIONAL :: Iters
    INTEGER, OPTIONAL :: Improvements

    INTEGER :: IOUnit
    LOGICAL :: GotIt
    CHARACTER(:), ALLOCATABLE :: Name
    
    Name = ListGetString(OptList,'Parameter Best File',GotIt )
    IF(.NOT. GotIt) RETURN
    
    CALL Info('SaveCurrentOptimum','Saving current optimum for later use!',Level=6)
    OPEN( NEWUNIT=IOUnit, FILE=Name, STATUS='UNKNOWN')
    WRITE (IOUnit,'(I0,T20,A)') n,'! Number of parameters'
    WRITE (IOUnit,*)         rpar(1:n)

    IF( PRESENT(Cost) )  &
        WRITE (IOUnit,'(ES12.6,T20,A)') Cost,'! Cost'
    IF( PRESENT(Iters) )  &
        WRITE (IOUnit,'(I0,T20,A)')   Iters,'! Iterations'
    IF( PRESENT(Improvements) )  &
        WRITE (IOUnit,'(I0,T20,A)')  Improvements,'! Improvements'
    CLOSE(IOUnit)
    
  END SUBROUTINE SaveCurrentOptimum
  
#ifdef HAVE_EXTOPTIM  
  SUBROUTINE ExternalOptimization_minpack(funvec)

    USE minpack_module

    IMPLICIT NONE 
    
    PROCEDURE(func) :: funvec
    INTEGER :: i,npar,niter,iflag 
    REAL(KIND=dp), ALLOCATABLE :: rpar(:), fvec(:)
    REAL(KIND=dp) :: xtol, epsfcn
    CHARACTER(:), ALLOCATABLE :: str
    TYPE(ValueList_t), POINTER :: OptList
    LOGICAL :: Found

    CALL Info('ExternalOptimization','Calling HYBRD from MinPack package')

    OptList => CurrentModel % Control
    
    niter = ListGetInteger( OptList,'Run Control Iterations', Found )

PRINT *,'niter minpack:',niter
    npar = ListGetInteger( OptList,'Parameter Count',Found )
    xtol = ListGetConstReal( OptList,'Optimization Tolerance',Found)
    IF(.NOT. Found ) xtol = 1.0e-6
    epsfcn = ListGetConstReal( OptList,'Run Control Variation',Found )
    IF(.NOT. Found) epsfcn = 0.01_dp

    ALLOCATE(rpar(npar),fvec(npar))        
    rpar = 1.0_dp
    fvec = 0.0_dp

    CALL GetSavedOptimum(OptList,rpar,Found)
    IF(.NOT. Found ) THEN    
      DO i=1,npar
        str = 'Initial Parameter '//I2S(i)
        rpar(i) = ListGetConstReal(OptList,str,Found) 
      END DO
    END IF
      
    CALL MinPack_HYBRD_Wrapper(npar,npar,rpar,fvec,niter,xtol,epsfcn)

    CALL SaveCurrentOptimum(OptList,npar,rpar)

  CONTAINS
    
    SUBROUTINE MinPack_HYBRD_Wrapper(n,ldfjac,x,fvec,maxfev,xtol,epsfcn) 

      implicit none

      PROCEDURE(func) :: fcn    ! function evaluating the cost

      INTEGER :: n              ! size of problem
      REAL(dp) :: x(n)          ! parameters to be optimized
      INTEGER :: ldfjac         ! size of jacobian
      INTEGER :: maxfev         ! max calls
      REAL(dp) :: xtol          ! convergence tolerance
      REAL(dp) :: fvec(n)       ! function values at optimum
      REAL(dp) :: epsfcn  ! step-size for Jacobian computation

      INTEGER :: ml
      INTEGER :: mu
      INTEGER :: mode   
      INTEGER :: nprint
      INTEGER :: info           
      INTEGER :: nfev  ! number of calls realized    
      INTEGER :: lr      
      REAL(dp) :: factor       
      REAL(dp) :: diag(n)   
      REAL(dp) :: fjac(ldfjac,n) 
      REAL(dp) :: r(n*(n+1)/2)
      REAL(dp) :: qtf(n)    
      REAL(dp) :: wa1(n)  
      REAL(dp) :: wa2(n)  
      REAL(dp) :: wa3(n)  
      REAL(dp) :: wa4(n)  

      ml = n-1
      mu = n-1
      mode = 1
      nprint = 0
      info = 0
      nfev = 0    
      lr = n*(n+1)/2
      factor = 100.0_dp
      diag = 1.0_dp
      fjac = 0.0_dp
      r = 0.0_dp
      qtf = 0.0_dp
      wa1 = 0.0_dp; wa2 = 0.0_dp; wa3 = 0.0_dp; wa4 = 0.0_dp

      CALL hybrd(funvec, n, x, Fvec, Xtol, Maxfev, Ml, Mu, Epsfcn, Diag, Mode, &
          Factor, Nprint, Info, Nfev, Fjac, Ldfjac, r, Lr, Qtf, Wa1, &
          Wa2, Wa3, Wa4)

    END SUBROUTINE MinPack_HYBRD_Wrapper
             
  END SUBROUTINE ExternalOptimization_Minpack


  
  SUBROUTINE ExternalOptimization_newuoa(funcost)

    USE newuoa_module

    IMPLICIT NONE 
      
    PROCEDURE(func) :: funcost
    INTEGER :: i,npar,npt,niter,iprint 
    
    REAL(KIND=dp), ALLOCATABLE :: rpar(:)
    REAL(KIND=dp) :: xtol, rhobeg, rhoend, minv, maxv
    CHARACTER(:), ALLOCATABLE :: str
    TYPE(ValueList_t), POINTER :: OptList
    LOGICAL :: Found, Found2

    CALL Info('ExternalOptimization','Calling NEWUOA from PowellOpt package')

    OptList => CurrentModel % Control
    
    niter = ListGetInteger( OptList,'Run Control Iterations', Found )
    npar = ListGetInteger( OptList,'Parameter Count',Found )
    xtol = ListGetConstReal( OptList,'Optimization Tolerance',Found)
    IF(.NOT. Found ) xtol = 1.0e-6

    npt = ListGetInteger( OptList,'Powell Interpolation Conditions', Found )
    npt = MIN(MAX(npar+2,npt),(npar+1)*(npar+2)/2)
    
    rhobeg = ListGetConstReal( OptList,'Initial Stepsize',UnfoundFatal=.TRUE.)
    rhoend = ListGetConstReal( OptList,'Min Stepsize',UnfoundFatal=.TRUE.)
    
    ALLOCATE(rpar(npar))
    rpar = 1.0_dp
    
    CALL GetSavedOptimum(OptList,rpar,Found)    

    IF(.NOT. Found ) THEN
      DO i=1,npar
        str = 'Initial Parameter '//I2S(i)
        rpar(i) = ListGetConstReal(OptList,str,Found)
        IF(.NOT. Found ) THEN
          str = 'Min Parameter '//I2S(i)
          minv = ListGetCReal(OptList,str,Found)
          str = 'Max Parameter '//I2S(i)
          maxv = ListGetCReal(OptList,str,Found2) 
          IF(Found .AND. Found2 ) THEN
            rpar(i) = (minv+maxv) / 2
          ELSE
            CALL Fatal('ExternalOptimization','Give either "Initial Parameter" or min/max range')
          END IF
        END IF
      END DO
    END IF
      
    i = ListGetInteger(CurrentModel % Simulation,'Max Output Level',Found )
    iprint = MIN(i/5,3) 
        
    CALL newuoa (npar, npt, rpar, rhobeg, rhoend, iprint, niter, funcost )
    
    CALL SaveCurrentOptimum(OptList,npar,rpar)

  END SUBROUTINE ExternalOptimization_newuoa


  SUBROUTINE ExternalOptimization_bobyqa(funcost)

    USE bobyqa_module

    IMPLICIT NONE 
      
    PROCEDURE(func) :: funcost
    INTEGER :: i,npar,npt,niter,iprint 
    
    REAL(KIND=dp), ALLOCATABLE :: rpar(:),xl(:),xu(:)
    REAL(KIND=dp) :: xtol, rhobeg, rhoend
    CHARACTER(:), ALLOCATABLE :: str
    TYPE(ValueList_t), POINTER :: OptList
    LOGICAL :: Found

    CALL Info('ExternalOptimization','Calling BOBYQA from PowellOpt package')

    OptList => CurrentModel % Control
    
    niter = ListGetInteger( OptList,'Run Control Iterations', Found )
    npar = ListGetInteger( OptList,'Parameter Count',Found )
    xtol = ListGetConstReal( OptList,'Optimization Tolerance',Found)
    IF(.NOT. Found ) xtol = 1.0e-6

    npt = ListGetInteger( OptList,'Optimization Interpolation Conditions', Found )
    npt = MIN(MAX(npar+2,npt),(npar+1)*(npar+2)/2)
    
    rhobeg = ListGetConstReal( OptList,'Initial Stepsize',UnfoundFatal=.TRUE.)
    rhoend = ListGetConstReal( OptList,'Min Stepsize',UnfoundFatal=.TRUE.)
    
    ALLOCATE(rpar(npar),xl(npar),xu(npar))
    rpar = 1.0_dp
    
    DO i=1,npar
      str = 'Min Parameter '//I2S(i)
      xl(i) = ListGetCReal(OptList,str,UnfoundFatal=.TRUE.)
      str = 'Max Parameter '//I2S(i)
      xu(i) = ListGetCReal(OptList,str,UnfoundFatal=.TRUE.)
      str = 'Initial Parameter '//I2S(i)
    END DO
    
    CALL GetSavedOptimum(OptList,rpar,Found)
    IF(.NOT. Found ) THEN
      DO i=1,npar
        rpar(i) = ListGetConstReal(OptList,str,Found)
        IF(.NOT. Found) rpar(i) = (xl(i)+xu(i))/2
      END DO
    END IF
      
    i = ListGetInteger(CurrentModel % Simulation,'Max Output Level',Found )
    iprint = MIN(i/5,3) 

    CALL bobyqa (npar, npt, rpar, xl, xu, rhobeg, rhoend, iprint, niter, funcost )

    CALL SaveCurrentOptimum(OptList,npar,rpar)

    
  END SUBROUTINE ExternalOptimization_bobyqa
#endif

!------------------------------------------------------------------------------
!> Adds parameters used in the simulation either predefined or from run control.
!> The idea is to make parametrized simulations more simple to perform. 
!------------------------------------------------------------------------------
 SUBROUTINE ControlParameters(Params,piter,GotParams,FinishEarly,&
     PostSimulation,SetCoeffs)

   IMPLICIT NONE
   
   TYPE(ValueList_t), POINTER :: Params
   INTEGER :: piter
   LOGICAL :: GotParams,FinishEarly
   LOGICAL, OPTIONAL :: PostSimulation
   LOGICAL, OPTIONAL :: SetCoeffs

   LOGICAL :: DoOptim, OptimalFinish, OptimalStart
   INTEGER :: NoParam, NoValues, cnt
   REAL(KIND=dp), ALLOCATABLE :: Param(:), BestParam(:)
   REAL(KIND=dp) :: Cost = HUGE( Cost ) 
   LOGICAL :: Found, GotCost, MinCost
   CHARACTER(*), PARAMETER :: Caller = 'ControlParameters'

   SAVE Cost, Param, BestParam

     
   NoParam = ListGetInteger( Params,'Parameter Count',Found )
   IF(.NOT. Found ) THEN
     NoParam = ListGetInteger( Params,'Number of Parameters',Found)
   END IF
   IF(NoParam == 0 ) THEN
     CALL Info(Caller,'No parameters to set in "Run Control" loop!',Level=4)
     RETURN
   END IF
   FinishEarly = .FALSE.

   ! The MATC parameters must be present before reading the sif file
   ! The coefficients must be set after reading the sif file.
   ! Hence we need a second, later, slot for the coefficient setup. 
   IF( PRESENT(SetCoeffs)) THEN
     IF( SetCoeffs ) THEN       
       CALL SetRealParametersKeywordCoeff(NoParam,Param,cnt)
       CALL Info(Caller,'Set '//I2S(cnt)//&
           ' coefficients with parameter tags!',Level=12)
       RETURN
     END IF
   END IF
      
   CALL Info(Caller, '-----------------------------------------', Level=5 )
   CALL Info(Caller, 'Setting sweeping parameters for simulation',Level=4 )
   
   NoValues = ListGetInteger( Params,'Run Control Iterations')

   IF( .NOT. ALLOCATED( Param ) ) THEN
     ALLOCATE( Param(NoParam), BestParam(NoParam) )
   END IF
   
   ! Visit this after simulation and register the parameters
   ! and cost function if present. We use same subroutine so
   ! we can take use of local data. 
   !----------------------------------------------------------
   IF( PRESENT( PostSimulation ) ) THEN
     IF( PostSimulation ) THEN
       CALL GetCostFunction(Params,Cost,GotCost)
       IF( GotCost ) CALL RegisterCurrentOptimum(Params,Cost) 
       CALL SaveParameterHistory()
       RETURN
     END IF
   END IF
      
   ! Here we set the parameters in different ways.
   ! They may be predefined or set by some optimization method. 
   !-------------------------------------------------------------------
   DoOptim = ListCheckPresent( Params,'Optimization Method')
   OptimalStart = ListGetLogical(Params,'Optimal Restart',Found )
   OptimalFinish = ListGetLogical( Params,'Optimal Finish',Found ) 

   IF( OptimalFinish .AND. piter == NoValues ) THEN
     CALL Info(Caller,'Performing the last step with the best so far')
     Param = BestParam
   ELSE IF( DoOptim ) THEN
     Found = .FALSE.
     IF( piter == 1 ) CALL GetSavedOptimum(Params,Param,Found)
     IF(.NOT. Found ) THEN
       CALL SetOptimizationParameters(Params,piter,GotParams,FinishEarly,&
           NoParam,Param,Cost)
     END IF
   ELSE
     CALL SetTabulatedParameters(Params,piter,GotParams,FinishEarly,&
         NoParam,Param)
   END IF

   IF( InfoActive(20) ) THEN
     PRINT *,'Parameters:',NoParam,Param
   END IF

   ! Set parameters to be accessible to the MATC preprocessor when reading sif file.
   CALL SetRealParametersMATC(NoParam,Param)

   CALL Info(Caller, '-----------------------------------------', Level=5 )


 CONTAINS



   

   ! We may register the current optimum and save it for later use.
   ! Then when starting over we may continue from the best so far.
   !----------------------------------------------------------------
   SUBROUTINE RegisterCurrentOptimum(OptList,Cost) 
     TYPE(ValueList_t), POINTER :: OptList
     REAL(KIND=dp) :: Cost
     
     REAL(KIND=dp) :: MinCost = HUGE(MinCost)
     INTEGER :: NoBetter = 0, i, IOUnit
     LOGICAL :: GotIt
     CHARACTER(:), ALLOCATABLE :: Name
     
     SAVE MinCost , NoBetter
     
     IF( ABS( Cost ) > ABS( MinCost ) ) RETURN
          
     ! Found a new best parameter combination
     !---------------------------------------
     MinCost = Cost
     BestParam(1:NoParam) = Param(1:NoParam)
     NoBetter = NoBetter + 1
     
     WRITE(Message,'(A,ES15.6E3)') 'Found New Minimum Cost:',MinCost
     CALL Info(Caller,Message,Level=4)
     
     CALL SaveCurrentOptimum(OptList,NoParam,BestParam,Cost,piter,NoBetter)

   END SUBROUTINE RegisterCurrentOptimum


   ! It may be interesting to follow the convergence history of the optimization.
   ! This file may include the parameters and the cost function if given.
   !-----------------------------------------------------------------------------
   SUBROUTINE SaveParameterHistory()
     
     LOGICAL :: DoAppend
     INTEGER :: IOunit
     LOGICAL :: GotIt
     CHARACTER(:), ALLOCATABLE :: Name
     
     ! Save the results to a file
     !---------------------------
     Name = ListGetString(Params,'Parameter History File',GotIt )
     IF(.NOT. GotIt ) RETURN

     DoAppend = ListGetLogical(Params,'Parameter History Append',GotIt )
          
     IF(piter == 1 .AND. .NOT. DoAppend ) THEN
       OPEN (NEWUNIT=IOUnit, FILE=Name)
     ELSE
       OPEN (NEWUNIT=IOUnit, FILE=Name,POSITION='APPEND')
     END IF

     IF( GotCost ) THEN
       WRITE (IOUnit,*) piter, Cost, Param(1:NoParam)
     ELSE
       WRITE (IOUnit,*) piter, Param(1:NoParam)
     END IF
     CLOSE(IOUnit)

   END SUBROUTINE SaveParameterHistory
      
 END SUBROUTINE ControlParameters

 !--------------------------------------------------------------------------------
 !> This subroutine sets tabulated parameters given in space separated ascii file
 !> or alternative in Dakota format file. 
 !------------------------------------------------------------------------------
 SUBROUTINE SetTabulatedParameters(Params,piter,GotParams,&
     FinishEarly,NoParam,Param)
   !-----------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: piter
    LOGICAL :: GotParams,FinishEarly
    INTEGER :: NoParam
    REAL(KIND=dp), ALLOCATABLE :: Param(:)

    LOGICAL :: Found, HaveFile, HaveArray
    REAL(KIND=dp), POINTER :: PArray(:,:)
    TYPE(Variable_t), POINTER :: PVar
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: LineCounter = 0
    CHARACTER(:), ALLOCATABLE :: FileName
    CHARACTER(*), PARAMETER :: Caller = 'SetTabulatedParameters'

        
    GotParams = .FALSE.
    FinishEarly = .FALSE.
    
    Parray => ListGetConstRealArray( Params,'Parameter Array',HaveArray)
    FileName = ListGetString( Params,'Parameter File',HaveFile)
    
    IF(.NOT. (HaveFile .OR. HaveArray) ) RETURN 
    
    CALL Info(Caller,'Trying to set simulation parameters')
    
    ! a) use given vector of simulation parameters
    IF( HaveArray ) THEN      
      CALL Info(Caller,'Setting parameters using constant array!',Level=6)
      IF( piter > SIZE(Parray,1) ) THEN
        FinishEarly = .TRUE.
      ELSE
        IF( SIZE(Parray,2) < NoParam ) THEN
          CALL Fatal(Caller,'Given "Parameter Array" has too few parameters!')
        END IF
        Param(1:NoParam) = Parray(piter,1:NoParam)
      END IF
    END IF
    
    ! b) read a row from file
    IF( HaveFile ) THEN
      CALL Info(Caller,'Setting parameters using external file!',Level=6)
      CALL ReadTabulatedParameters()
    END IF
    
    IF( FinishEarly ) THEN
      CALL Warn(Caller,'Parameters exhausted already: '//I2S(piter))  
      RETURN
    END IF
        
    GotParams = .TRUE.
    
  CONTAINS

    
    SUBROUTINE ReadTabulatedParameters()

      INTEGER :: FileUnit, Line, NOffset, FileTypeInd, FileRow, iostat, i, j, k
      REAL(KIND=dp), ALLOCATABLE :: TmpValues(:)
      CHARACTER(:), ALLOCATABLE :: FileType
      CHARACTER(LEN=MAX_NAME_LEN) :: readstr 
          
      FileType = ListGetString( Params,'Parameter Filetype',Found )
      FileTypeInd = 0
      IF( Found ) THEN
        SELECT CASE( FileType )
        CASE('table')
          FileTypeInd = 0
        CASE('dakota')
          FileTypeInd = 1
        CASE DEFAULT
          CALL Fatal(Caller,'Unkonown filetype: '//TRIM(FileType))
        END SELECT
      END IF

      FileRow = ListGetInteger( Params,'Parameter Row Offset',Found ) 
      FileRow = FileRow + piter
      
      OPEN(NEWUNIT=FileUnit,FILE=FileName,IOSTAT=iostat)
      IF( iostat /= 0 ) THEN
        CALL Fatal(Caller,'Could not open file: '//TRIM(FileName))
      END IF
        
      IF( FileTypeInd == 1 ) THEN
        Noffset = 2
        DO WHILE(.TRUE.) 
          READ( FileUnit,'(A)',IOSTAT=iostat) readstr
          IF( iostat /= 0 ) THEN
            CALL Fatal(Caller,'Could not read dummy line: '//I2S(Line))
          END IF
          i = INDEX( readstr,'RUN NO.') 
          IF( i > 0 ) THEN
            CALL Info(Caller,'Parameter lines start after line: '//TRIM(readstr),Level=6)
            EXIT
          END IF
          i = INDEX( readstr,'Number of Variables =' )
          IF( i > 0 ) THEN
            j = MAX_NAME_LEN
            READ( readstr(i+21:j),*,IOSTAT=iostat) k
            IF( iostat /= 0 ) THEN
              CALL Fatal(Caller,'Could not read parameters from line: '//TRIM(readstr))
            END IF
            CALL Info(Caller,'Number of parameters in DAKOTA file: '&
                //I2S(k),Level=6)
            IF( k < NoParam ) THEN
              CALL Fatal(Caller,'Dakota file has too few parameters!')
            END IF
          END IF
        END DO
      ELSE
        Noffset = ListGetInteger( Params,'Parameter Column Offset',Found ) 
      END IF

      Line = 0
      DO WHILE(.TRUE.)
        Line = Line + 1
        READ( FileUnit,'(A)',IOSTAT=iostat) readstr
        IF( iostat /= 0 ) THEN
          CALL Warn(Caller,'Could not read parameter line: '//I2S(Line))
          CLOSE(FileUnit)
          FinishEarly = .TRUE.
          RETURN
        END IF
        IF( Line == FileRow ) THEN
          CALL Info(Caller,'Read parameter line: '//I2S(Line),Level=7)
          EXIT
        END IF
      END DO
      CLOSE(FileUnit)
      
      ALLOCATE( TmpValues(Noffset+NoParam) )
      READ(readstr,*,IOSTAT=iostat) TmpValues(1:Noffset+NoParam)
      IF( iostat /= 0 ) THEN
        CALL Warn(Caller,'Could not read parameters from: '//TRIM(readstr))
        FinishEarly = .TRUE.
        RETURN
      END IF
            
      Param(1:NoParam) = TmpValues(NOffset+1:Noffset+NoParam)
      
      CALL Info(Caller,'Parameters read from file',Level=8)
      
    END SUBROUTINE ReadTabulatedParameters
        
!------------------------------------------------------------------------------
  END SUBROUTINE SetTabulatedParameters
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> This routine allows optimization within single ElmerSolver execution.
!> Several simple algorithms are provided but for more complex ones look elsewhere.
!------------------------------------------------------------------------------
  SUBROUTINE SetOptimizationParameters(OptList,piter,GotParams,FinishEarly, &
      NoParam,Param,Cost)
    
    IMPLICIT NONE
    
    TYPE(ValueList_t), POINTER :: OptList
    INTEGER :: piter
    LOGICAL :: GotParams, FinishEarly
    INTEGER :: NoParam
    REAL(KIND=dp) :: Param(:)
    REAL(KIND=dp) :: ParamCost    
    
    LOGICAL :: gotIt, GotIt2, GotInit, InternalHistory, Visited = .FALSE.
    LOGICAL, ALLOCATABLE :: FixedParam(:)
    INTEGER :: i,j,k,l,NoValues, NoFreeParam, NoOpt, &
        Direction=1, NoImprovements=0, OptimizationsDone
    REAL(KIND=dp), ALLOCATABLE :: MinParam(:), MaxParam(:), dParam(:), PrevParam(:,:), &
        PrevCost(:), BestParam(:), InitParam(:)
    REAL(KIND=dp) :: Cost, MinCost, x(10), c(10), minv, maxv, OptTol
    TYPE(Variable_t),POINTER :: Var
    INTEGER :: IOUnit
    CHARACTER(:), ALLOCATABLE :: Name, ParamStr, Method
    CHARACTER(*), PARAMETER :: Caller = 'SetOptimizationParameters'

    
    SAVE MinParam, MaxParam, PrevParam, &
        Method, Direction, x, c, PrevCost, &
        FixedParam, NoFreeParam, MinCost, BestParam, NoValues, &
        InternalHistory, dParam, &
        NoImprovements, OptTol, Visited

    GotParams = .TRUE.
    
    !------------------------------------------------------------------------------
    ! In the 1st round perform initializations 
    !------------------------------------------------------------------------------
    IF(.NOT. Visited ) THEN
      CALL Info(Caller,'Initializing solver for optimization')

      NoValues = ListGetInteger( OptList,'Run Control Iterations')

      OptTol = ListGetConstReal( OptList,'Optimization Tolerance',GotIt)

      ALLOCATE( MinParam(NoParam), BestParam(NoParam), MaxParam(NoParam), &
          dParam(NoParam), FixedParam(NoParam), InitParam(NoParam))

      MinParam = -HUGE( MinParam ) 
      BestParam = 0.0_dp
      MaxParam = HUGE( MaxParam ) 
      dParam = 0.0_dp
      FixedParam = .FALSE.
      MinCost = HUGE(MinCost)

      NoFreeParam = 0
      DO i=1,NoParam
        ParamStr = 'Parameter '//I2S(i)

        FixedParam(i) = ListGetLogical(OptList,'Fixed '//TRIM(ParamStr),GotIt)
        InitParam(i) = ListGetConstReal(OptList,'Initial '//TRIM(ParamStr),GotInit)

        IF(.NOT. ( FixedParam(i) ) ) THEN
          minv = ListGetConstReal(OptList,'Min '//TRIM(ParamStr),GotIt)
          maxv = ListGetConstReal(OptList,'Max '//TRIM(ParamStr),GotIt2)

          IF( GotIt ) MinParam(i) = minv
          IF( GotIt2 ) MaxParam(i) = maxv

          ! if both min and max given then the 1st set of parameters are
          ! the average, otherwise either extremum. 
          IF( .NOT. GotInit ) THEN
            IF( GotIt .AND. GotIt2 ) THEN
              InitParam(i) = 0.5_dp * ( minv + maxv ) 
            ELSE IF( GotIt ) THEN
              InitParam(i) = minv 
            ELSE IF( GotIt2 ) THEN
              InitParam(i) = maxv 
            END IF
          END IF
        END IF
        dParam(i) = ListGetConstReal(OptList,'Scale '//TRIM(ParamStr),GotIt)
        IF(.NOT. GotIt) dParam(i) = 1.0_dp
      END DO

      NoFreeParam = NoParam - COUNT(FixedParam)
      IF( NoFreeParam == 0 ) THEN
        CALL Warn(Caller,'All parameters are fixed, no optimization!')
        RETURN
      END IF

      Method = ListGetString(OptList,'Optimization Method')
      
      ! Internal history could be used in more complicated optimization routines
      !--------------------------------------------------------------------------
      InternalHistory = ListGetLogical( OptList,'Internal History',GotIt)    
      IF( Method == 'bisect') InternalHistory = .TRUE.
      IF( InternalHistory ) THEN
        ALLOCATE( PrevParam(NoValues,NoParam), PrevCost(NoValues))
      END IF

      Visited = .TRUE.
    END IF

    IF( piter == 1 ) THEN
      Param(1:NoParam) = InitParam(1:NoParam)
    END IF

    IF( InternalHistory ) THEN
      OptimizationsDone = OptimizationsDone + 1
      PrevParam(OptimizationsDone,1:NoParam) = Param(1:NoParam)
      PrevCost(OptimizationsDone) = Cost
    END IF
            
    WRITE( Message, '(A,I0,A,A)' ) 'Manipulating ',NoFreeParam,' parameters using ',TRIM(Method) 
    CALL Info(Caller, Message, Level=4 )


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
      
    CASE ('hybrd','newuoa','bobyqa')
      IF( piter > 1 ) THEN
        CALL Fatal(Caller,'You should end up here with external methods!')
      END IF
        
    CASE DEFAULT
      CALL Fatal(Caller,'Unknown method')

    END SELECT

    IF(.FALSE.) THEN
      DO i=1,NoParam 
        IF( FixedParam(i) ) CYCLE
        Param(i) = MAX(MinParam(i),Param(i))
        Param(i) = MIN(MaxParam(i),Param(i))
      END DO
    END IF
    
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
        popsize = ListGetInteger(OptList,'Population Size',GotIt)
        IF(.NOT. GotIt) popsize = 5 * parsize
        popcoeff = ListGetConstReal(OptList,'Population Coefficient',GotIt)
        IF(.NOT. GotIt) popcoeff = 0.7
        popcross = ListGetConstReal(OptList,'Population Crossover',GotIt)
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
        IF(NoFreeParam /= 1) CALL Fatal(Caller,&
            'Option scan implemented only for one parameter')
        DO i=1,NoParam
          IF(.NOT. FixedParam(i)) EXIT
        END DO
        CALL Info(Caller,'Applying scanning to parameter '//I2S(i),Level=5)
        maxno = NoValues 
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

      IF(NoFreeParam /= 1) CALL Fatal(Caller,&
          'Option bisect implemented only for one parameter')
      IF( no == 0 ) THEN
        DO j=1,NoParam
          IF(.NOT. FixedParam(j)) EXIT
        END DO
        CALL Info(Caller,'Applying bisection search to parameter '//I2S(j),Level=7)
      END IF

      no = no + 1

      IF(no == 1) THEN
        step = ListGetConstReal(OptList,'Initial Stepsize',GotIt)
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
        CALL Fatal(Caller,'Bisection method cannot handle local maxima')
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

      IF(NoFreeParam /= 1) CALL Fatal(Caller,&
          'Secant search implemented only for one parameter')
      IF( no == 0 ) THEN
        DO j=1,NoParam
          IF(.NOT. FixedParam(j)) EXIT
        END DO
        CALL Info(Caller,'Applying secant search to parameter '//I2S(j),Level=7)
      END IF

      no = no + 1

      IF(no == 1) THEN
        maxstep = ListGetConstReal(OptList,'Max StepSize',GotIt)
        step = ListGetConstReal(OptList,'Initial StepSize',GotIt)
        IF(.NOT. GotIt) step = 1.0d-3*(MaxParam(j)-MinParam(j))
        relax = ListGetConstReal(OptList,'StepSize Relaxation Factor',GotIt)
        IF(.NOT. GotIt) Relax = 1.0_dp
      END IF

      x0 = x1
      x1 = x2
      f0 = f1 
      f1 = Cost

      IF(no <= 2) THEN
        x2 = Param(j) + (no-1)*step
      ELSE IF( ABS(f1) < OptTol ) THEN
        CALL Info(Caller,'Secent search tolerance reached, doing nothing')
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

        lambda = ListGetConstReal(OptList,'Simplex Relative Length Scale',Found)
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

        nomax = ListGetInteger(OptList,'Simplex Restart Interval',Found)
        maxratio = ListGetConstReal(OptList,'Simplex Restart Convergence Ratio',Found)      
        IF(.NOT. Found) maxratio = 1.0_dp
        AllocationsDone = .TRUE.
      END IF

      IF( nomax > 0 .AND. no > nomax .AND. ratio > maxratio ) THEN
        CALL Info(Caller,'Making a restart in simplex')
        ratio = 0.0_dp
        ls = 0.1 * ls

        !PRINT *,'Simplex: coeff',ls
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
        CALL Warn(Caller,'Srink mode not ok yet!')
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

    
  END SUBROUTINE SetOptimizationParameters

!------------------------------------------------------------------------------
!> Optionally we may revert to initial coordinates when using "Run Control".
!------------------------------------------------------------------------------
 SUBROUTINE ControlResetMesh(Params,piter)

   IMPLICIT NONE
   
   TYPE(ValueList_t), POINTER :: Params
   INTEGER :: piter

   LOGICAL :: Found
   TYPE(Nodes_t) :: Nodes0
   TYPE(Nodes_t), POINTER :: Nodes
   INTEGER :: n

   SAVE Nodes0
   
   IF( ListGetLogical( Params,'Reset Mesh Coordinates',Found ) ) THEN
     Nodes => CurrentModel % Mesh % Nodes
     n = SIZE( Nodes % x )
     IF( piter == 1 ) THEN
       ALLOCATE( Nodes0 % x(n), Nodes0 % y(n), Nodes0 % z(n) )
       Nodes0 % x = Nodes % x
       Nodes0 % y = Nodes % y
       Nodes0 % z = Nodes % z
     ELSE
       Nodes % x = Nodes0 % x
       Nodes % y = Nodes0 % y
       Nodes % z = Nodes0 % z
     END IF
   END IF

   
 END SUBROUTINE ControlResetMesh

END MODULE OptimizationUtils
 
