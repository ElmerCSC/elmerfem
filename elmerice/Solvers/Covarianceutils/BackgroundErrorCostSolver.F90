!***********************************************************************************************      
! Compute a cost function from a background as Cost=0.5 * (x-x^b). B^-1 .(x-x^b)
!  x is the optimised variable; x^b the background
!  B^1 is the background error covariance matrix:
!    here B= S . C . S
!     with:
!        - S is a diagonale matrix with the standard deviation (assumed constant for now)
!        - C is a correlation matix
!  Available choices for C "Covariance type = String ..."
!        - diagonal; i.e. C=I and B=S^2 is diagonal with the variances
!        - "full matrix" : C is computed from standard correlation
!            functions and inverted using lapack routines
!        - "diffusion operator" : C is approximated with the diffusion operator approach
!  Current limitaions :
!    - 2D mesh; TODO to run it on a 2D surface boundary?
!    - Serial fro the full-matrix approach
!
! Rq.
!  - IF x has DOFs > 1 we apply independandtly the same B^-1
!  - IF 2 instances of the same solver are used in the same .sif mae a
!      copy of the lib as things are initialised and saved.... 
!
! TODO:
!   - add mandatory keywords at init, e.g. variable, ...
!      
!***********************************************************************************************         
      SUBROUTINE BackgroundErrorCostSolver( Model,Solver,dt,TransientSimulation )
!***********************************************************************************************
      USE GeneralUtils
      USE CovarianceUtils
      IMPLICIT NONE
!------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver

      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE(ValueList_t), POINTER :: SolverParams

      TYPE(Variable_t), POINTER :: Var,Var_b,DJDVariable,CostVar
      REAL(KIND=dp), POINTER ::  Values(:), Values_b(:),DJDValues(:)
      INTEGER, POINTER :: Perm(:),Perm_b(:),DJDPerm(:)
      CHARACTER(LEN=MAX_NAME_LEN) :: Varname,VarbName,GradVarName,CostName
      INTEGER :: DOFs

      INTEGER :: i,j,k,l,t
      INTEGER :: ierr
      INTEGER :: MeshDim

      TYPE(Solver_t), POINTER, SAVE :: MSolver,KMSolver
      REAL(kind=dp),allocatable,save :: aap(:) ! matrix in packed format
      REAL(kind=dp),allocatable,save :: x(:),y(:)
      REAL(kind=dp),allocatable,save :: One(:)
      REAL(kind=dp) :: Cost,Cost_S

      INTEGER,SAVE  :: nn
      INTEGER, ALLOCATABLE, SAVE :: ActiveNodes(:),InvPerm(:)
      INTEGER,SAVE :: PbDim

      CHARACTER(LEN=MAX_NAME_LEN),SAVE :: CovType
      REAL(kind=dp),SAVE :: std
      CHARACTER(LEN=MAX_NAME_LEN) :: Ctype
      REAL(kind=dp) :: CRange,Cm
      INTEGER :: p
      REAL(kind=dp) :: sigma2
      INTEGER :: Op

      LOGICAL, SAVE :: Firsttime=.TRUE.
      LOGICAL :: Reset
      LOGICAL :: Found
      Logical :: Parallel

      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="CostCov"
      CHARACTER(LEN=MAX_NAME_LEN) :: CostFile

      ! check Parallel/Serial
      Parallel=(ParEnv % PEs > 1)


      SolverParams => GetSolverParams()

      VarName = ListGetString(SolverParams,"Variable name",UnFoundFatal=.TRUE.)
      GradVarName = ListGetString(SolverParams,"Gradient Variable name",UnFoundFatal=.TRUE.)
      VarbName = ListGetString(SolverParams,"Background Variable name",UnFoundFatal=.TRUE.)
      CostName = ListGetString(SolverParams,"Cost Variable name",UnFoundFatal=.TRUE.)

      CostFile = ListGetString(SolverParams,'Cost Filename',UnFoundFatal=.TRUE. )

      ! get handles on the variables
      Reset =  GetLogical( SolverParams,'Reset Cost Value', Found)
      IF(.NOT.Found) Reset=.True.

      Var => VariableGet( Solver % Mesh % Variables,TRIM(VarName), UnFoundFatal=.TRUE. )
      Values => Var % Values
      Perm => Var % Perm
      DOFs = Var % DOFs

      Var_b => VariableGet(Solver % Mesh % Variables,TRIM(VarbName),UnFoundFatal=.TRUE.)
      Values_b => Var_b % Values
      Perm_b => Var_b % Perm
      IF (Var_b % DOFs .NE. DOFs) &
              CALL FATAL(SolverName,"Var and Var_b have different dofs")

      DJDVariable => VariableGet( Solver % Mesh % Variables,TRIM(GradVarName), UnFoundFatal=.TRUE. )
      DJDValues => DJDVariable % Values
      DJDPerm => DJDVariable % Perm
      IF (DJDVariable  % DOFs .NE. DOFs) &
              CALL FATAL(SolverName,"Var and grad have different dofs")
      IF (Reset) DJDValues=0.0_dp

      CostVar => VariableGet( Solver % Mesh % Variables,TRIM(CostName),UnFoundFatal=.TRUE.)
      IF (Reset) CostVar % Values=0.0_dp

      !! some initialisation
      IF (Firsttime) THEN
        CALL GetActiveNodesSet(Solver,nn,ActiveNodes,InvPerm,PbDim)

        !! The covariance type
        CovType = ListGetString(SolverParams,"Covariance type",UnFoundFatal=.TRUE.)

       std = ListGetConstReal(SolverParams,"standard deviation",UnFoundFatal=.TRUE.)

       IF (ParEnv%MyPE.EQ.0) THEN
         open(10,file=TRIM(CostFile))
           write(10,*) '# Covariance type: ',TRIM(CovType)
           write(10,*) '# standard deviation: ',std
       END IF

        SELECT CASE (CovType)

          CASE('diagonal')
            CALL INFO(SolverName,"Using diagonal covariance",level=3)

          CASE('full matrix')
            CALL INFO(SolverName,"Using full matrix covariance",level=3)

            Op=3
            ALLOCATE(aap(nn*(nn+1)/2))
            CALL CovarianceInit(Solver,nn,InvPerm,aap,Op,PbDim)

            Ctype = ListGetString(SolverParams,"correlation type",UnFoundFatal=.TRUE.)
            IF (ParEnv%MyPE.EQ.0) &
              write(10,*) '# Correlation type: ',TRIM(Ctype)

            Crange = ListGetConstReal(SolverParams,"correlation range", UnFoundFatal=.TRUE.)

            IF (Ctype == 'maternp') THEN
              p = ListGetInteger(SolverParams,"MaternP polynomial order",UnFoundFatal=.TRUE.)
              IF (ParEnv%MyPE.EQ.0) &
               write(10,*) '# range, exponent: ',Crange,p
            ELSE IF (Ctype == 'materni') THEN
              p = ListGetInteger(SolverParams,"MaternI order",UnFoundFatal=.TRUE.)
              IF (ParEnv%MyPE.EQ.0) &
               write(10,*) '# range, exponent: ',Crange,p
            ELSE
               IF (ParEnv%MyPE.EQ.0) &
                 write(10,*) '# range:',Crange
            END IF

          CASE('diffusion operator')
            CALL INFO(SolverName,"Using diffusion operator covariance",level=3)

            CALL CovarianceInit(Solver,MSolver,KMSolver)

            Cm = ListGetInteger(SolverParams,"Matern exponent m",UnFoundFatal=.TRUE.)
            Crange = ListGetConstReal(SolverParams,"correlation range", UnFoundFatal=.TRUE.)

            IF (ParEnv%MyPE.EQ.0) &
              write(10,*) '# range, exponent: ',Crange,Cm

           CASE DEFAULT
             CALL FATAL(SolverName,"Covariance type not known")

        END SELECT

        IF (ParEnv%MyPE.EQ.0) &
         close(10)

       allocate(x(nn),y(nn),One(nn))
       ! normalisation vector; usefull for parallel
       One=1._dp
       IF (Parallel) CALL ParallelSumVector(Solver%Matrix, One)

       Firsttime=.FALSE.
      END IF

      Cost=0._dp
      DO i=1,DOFS

        !(x-x_b)
         x(Solver%Variable%Perm(ActiveNodes(1:nn))) = Values(DOFs*(Perm(ActiveNodes(1:nn))-1)+i)- &
                 Values_b(DOFs*(Perm_b(ActiveNodes(1:nn))-1)+i)

        SELECT CASE (CovType)

          CASE('diagonal')
            sigma2=std**2
            y(1:nn) = x(1:nn)/sigma2

          CASE('full matrix')
            CALL InvCovarianceVectorMultiply(Solver,nn,aap,x,y)

          CASE('diffusion operator')
            ! y = SIGMA^1 C^1 SIGMA^1 . (x-x_b)
            CALL InvCovarianceVectorMultiply(Solver,MSolver,KMSolver,nn,x,y)

        END SELECT

       ! gradient = SIGMA^1 C^1 SIGMA^1 . (x-x_b)
       ! gradients are gathered in the optimisation step; so also normalize by One.
        DJDValues(DOFS*(DJDPerm(ActiveNodes(1:nn))-1)+i)=DJDValues(DOFS*(DJDPerm(ActiveNodes(1:nn))-1)+i)+ & 
                 y(Solver%Variable%Perm(ActiveNodes(1:nn)))/One(Solver%Variable%Perm(ActiveNodes(1:nn)))


        ! 1/2 (x-x_b) . SIGMA^1 C^1 SIGMA^1 . (x-x_b)
        ! normalisation by one insure that shared nodes are not counted twice...
        Cost=Cost+0.5*SUM(x(1:nn)*y(1:nn)/One(1:nn))

      END DO

       IF (Parallel) THEN
           CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
          CostVar % Values(1)=CostVar % Values(1)+Cost_S
          IF (ParEnv % MyPE == 0) then
                 OPEN (10, FILE=CostFile,POSITION='APPEND')
                 write(10,'(2(ES20.11E3))') GetTime(),Cost_S
                 CLOSE(10)
         End if
       ELSE
          CostVar % Values(1)=CostVar % Values(1)+Cost
          open(10,file=TRIM(CostFile),position='append')
           write(10,'(2(ES20.11E3))') GetTime(),Cost
          close(10)
       END IF

     END SUBROUTINE BackgroundErrorCostSolver
