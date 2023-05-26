!***********************************************************************************************      
!***********************************************************************************************         
      SUBROUTINE CovarianceVectorMultiplySolver( Model,Solver,dt,TransientSimulation )
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

      TYPE(Variable_t), POINTER :: Var
      REAL(KIND=dp), POINTER ::  Values(:)
      INTEGER, POINTER :: Perm(:)
      CHARACTER(LEN=MAX_NAME_LEN) :: Varname
      INTEGER :: DOFs


      TYPE(Solver_t), POINTER, SAVE :: MSolver,KMSolver
      REAL(kind=dp),allocatable,save :: aap(:) ! matrix in packed format
      REAL(kind=dp),allocatable,save :: x(:),y(:),norm(:)

      INTEGER,SAVE  :: nn
      INTEGER, ALLOCATABLE, SAVE :: ActiveNodes(:),InvPerm(:)
      INTEGER,SAVE :: PbDim

      CHARACTER(LEN=MAX_NAME_LEN),SAVE :: CovType
      REAL(kind=dp),SAVE :: std
      REAL(kind=dp) :: sigma2
      INTEGER :: Op

      LOGICAL, SAVE :: Firsttime=.TRUE.
      LOGICAL, SAVE :: Normalize
      LOGICAL :: Parallel
      LOGICAL :: Found

      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="Cov.x Solver"

      ! check Parallel/Serial
      Parallel=(ParEnv % PEs > 1)

      SolverParams => GetSolverParams()

      VarName = ListGetString(SolverParams,"Input Variable",UnFoundFatal=.TRUE.)

      Normalize = ListGetLogical(SolverParams,"Normalize",Found)
      IF (.NOT.Found) Normalize=.FALSE.

      Var => VariableGet( Solver % Mesh % Variables,TRIM(VarName), UnFoundFatal=.TRUE. )
      Values => Var % Values
      Perm => Var % Perm
      DOFs = Var % DOFs
      IF (DOFS.GT.1) &
         CALL FATAL(SolverName,'Sorry 1DOFs variables')

      !! some initialisation
      IF (Firsttime) THEN

        CALL GetActiveNodesSet(Solver,nn,ActiveNodes,InvPerm,PbDim)

        !Sanity check
        IF (ANY(Perm(ActiveNodes(1:nn)).LT.0)) &
          CALL FATAL(SolverName,"Pb with input variable perm")

        !! The covariance type
        CovType = ListGetString(SolverParams,"Covariance type",UnFoundFatal=.TRUE.)
        std = ListGetConstReal(SolverParams,"standard deviation",UnFoundFatal=.TRUE.)

        SELECT CASE (CovType)

          CASE('diagonal')
            CALL INFO(SolverName,"Using diagonal covariance",level=3)

          CASE('full matrix')
            CALL INFO(SolverName,"Using full matrix covariance",level=3)

            Op=1
            ALLOCATE(aap(nn*(nn+1)/2))
            CALL CovarianceInit(Solver,nn,InvPerm,aap,Op,PbDim)

          CASE('diffusion operator')
            CALL INFO(SolverName,"Using diffusion operator covariance",level=3)

            CALL CovarianceInit(Solver,MSolver,KMSolver)

        END SELECT

       allocate(x(nn),y(nn))

       IF (Normalize) THEN
          allocate(norm(nn)) 

          !input vector
          x(:) = 1._dp

          ! C . x
          SELECT CASE (CovType)

            CASE('diagonal')
              sigma2=std**2
              norm(:)=sigma2*x(:)

            CASE('full matrix')
              CALL CovarianceVectorMultiply(Solver,nn,aap,x,norm)

            CASE('diffusion operator')
              ! y = SIGMA C SIGMA . x
              CALL CovarianceVectorMultiply(Solver,MSolver,KMSolver,nn,x,norm)

           END SELECT
       END IF

       Firsttime=.FALSE.
      END IF

      !input vector
      x(Solver%Variable%Perm(ActiveNodes(1:nn))) = Values(Perm(ActiveNodes(1:nn)))

      ! C . x
      SELECT CASE (CovType)

        CASE('diagonal')
          sigma2=std**2
          y(:)=sigma2*x(:)

        CASE('full matrix')
          CALL CovarianceVectorMultiply(Solver,nn,aap,x,y)

        CASE('diffusion operator')
          ! y = SIGMA C SIGMA . x
          CALL CovarianceVectorMultiply(Solver,MSolver,KMSolver,nn,x,y)

      END SELECT

      IF (Normalize) y(:)=y(:)/norm(:)

      Solver % Variable % Values(Solver%Variable%Perm(ActiveNodes(1:nn)))=&
                    y(Solver%Variable%Perm(ActiveNodes(1:nn)))

     END SUBROUTINE CovarianceVectorMultiplySolver
