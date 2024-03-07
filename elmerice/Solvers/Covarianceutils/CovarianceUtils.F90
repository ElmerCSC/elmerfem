      MODULE CovarianceUtils
      USE MainUtils
      USE DefUtils

      IMPLICIT NONE

      INTERFACE SqrCovarianceVectorMultiply
         MODULE PROCEDURE SqrCovarianceVectorMultiplyD,SqrCovarianceVectorMultiplyL
      END INTERFACE

      INTERFACE CovarianceVectorMultiply
        MODULE PROCEDURE CovarianceVectorMultiplyD,CovarianceVectorMultiplyL
      END INTERFACE

      INTERFACE InvCovarianceVectorMultiply
        MODULE PROCEDURE InvCovarianceVectorMultiplyD,InvCovarianceVectorMultiplyL
      END INTERFACE

      INTERFACE CovarianceInit
        MODULE PROCEDURE CovarianceInitD,CovarianceInitL
      END INTERFACE

      CONTAINS

        SUBROUTINE GetActiveNodesSet(Solver,n,ActiveNodes,InvPerm,PbDim)
          TYPE(Solver_t) :: Solver
          INTEGER, INTENT(OUT) :: n
          INTEGER, ALLOCATABLE, INTENT(OUT) :: ActiveNodes(:),InvPerm(:)
          INTEGER :: PbDim

          INTEGER :: t
          INTEGER :: counter
          INTEGER,POINTER :: Perm(:)
          LOGICAL :: BoundarySolver 

          IF (.NOT.ASSOCIATED(Solver % Matrix)) &
            CALL FATAL("GetActiveNodesSet","Matrix must be associated")
          Perm => Solver % Variable % Perm

          n = Solver % Matrix % NumberOfRows
          ALLOCATE(ActiveNodes(n),InvPerm(n))
 
          counter = 0
          DO t=1,Solver % Mesh % NumberOfNodes
            IF (Perm(t).LT.1) CYCLE
            counter  = counter + 1
            IF (counter.GT.n) &
               CALL FATAL("GetActiveNodesSet","Problem detected...")
            ActiveNodes(counter) = t
            InvPerm(Perm(t)) = t
          END DO

          BoundarySolver = ( Solver % ActiveElements(1) > Solver % Mesh % NumberOfBulkElements )
          IF (BoundarySolver) THEN
            PbDim = CoordinateSystemDimension() - 1
          ELSE
            PbDim = CoordinateSystemDimension()
          ENDIF
          
        END SUBROUTINE GetActiveNodesSet

!############################################################################################
!##
!## standard matrix vector multiplication using lapack
!##
!############################################################################################
        SUBROUTINE InvCovarianceVectorMultiplyL(Solver,n,aap,x,y)
          INTEGER,INTENT(IN) :: n
          TYPE(Solver_t) :: Solver
          REAL(kind=dp),INTENT(IN) :: aap(:)
          REAL(KIND=dp),DIMENSION(n),INTENT(IN) :: x
          REAL(KIND=dp),DIMENSION(n),INTENT(OUT) :: y

          REAL(KIND=dp),DIMENSION(n) :: x1
          TYPE(ValueList_t), POINTER :: SolverParams
          character*1,parameter :: uplo='L'
          REAL(KIND=dp) :: std

          SolverParams => GetSolverParams(Solver)

          std = ListGetConstReal(SolverParams,"standard deviation",UnFoundFatal=.TRUE.)

          ! x1 = Sigma-1 x
          x1(1:n)=x(1:n)/std

          ! y = C-1 . x1
          call dspmv(uplo,n,1._dp,aap,x1,1,0._dp,y,1)

          ! y = Sigma-1 y
          y(1:n)=y(1:n)/std

        END  SUBROUTINE InvCovarianceVectorMultiplyL

        SUBROUTINE CovarianceVectorMultiplyL(Solver,n,aap,x,y)
          INTEGER,INTENT(IN) :: n
          TYPE(Solver_t) :: Solver
          REAL(kind=dp),INTENT(IN) :: aap(:)
          REAL(KIND=dp),DIMENSION(n),INTENT(IN) :: x
          REAL(KIND=dp),DIMENSION(n),INTENT(OUT) :: y

          REAL(KIND=dp),DIMENSION(n) :: x1
          TYPE(ValueList_t), POINTER :: SolverParams
          character*1,parameter :: uplo='L'
          REAL(KIND=dp) :: std

          SolverParams => GetSolverParams(Solver)

          std = ListGetConstReal(SolverParams,"standard deviation",UnFoundFatal=.TRUE.)

          ! x1 = Sigma x
          x1(1:n)=std*x(1:n)

          ! y = C . x1
          call dspmv(uplo,n,1._dp,aap,x1,1,0._dp,y,1)

          ! y = Sigma y
          y(1:n)=std*y(1:n)

        END  SUBROUTINE CovarianceVectorMultiplyL

        SUBROUTINE SqrCovarianceVectorMultiplyL(Solver,n,aap,x,y)
          INTEGER,INTENT(IN) :: n
          TYPE(Solver_t) :: Solver
          REAL(kind=dp),INTENT(IN) :: aap(:)
          REAL(KIND=dp),DIMENSION(n),INTENT(IN) :: x
          REAL(KIND=dp),DIMENSION(n),INTENT(OUT) :: y

          TYPE(ValueList_t), POINTER :: SolverParams
          character*1,parameter :: uplo='L'
          REAL(KIND=dp) :: std

          SolverParams => GetSolverParams(Solver)

          std = ListGetConstReal(SolverParams,"standard deviation",UnFoundFatal=.TRUE.)

          ! y = L . x; matrix is  lower triangular
          call dtpmv(uplo,'N','N',n,aap,x,1)

          ! y = Sigma y
          y(1:n)=std*x(1:n)

        END  SUBROUTINE SqrCovarianceVectorMultiplyL

!############################################################################################
!##
!## Matrix vector multiplications  for the diffusion operator approach
!##
!############################################################################################
!-------------------------------------------------------------------------------------------
!  Covariance product : y=C.x
!-------------------------------------------------------------------------------------------
         SUBROUTINE CovarianceVectorMultiplyD(Solver,MSolver,KMSolver,n,x,y)
          TYPE(Solver_t) :: Solver
          TYPE(Solver_t), POINTER :: MSolver,KMSolver
          INTEGER,INTENT(IN) :: n
          REAL(KIND=dp),DIMENSION(n),INTENT(IN) :: x
          REAL(KIND=dp),DIMENSION(n),INTENT(OUT) :: y

          TYPE(Matrix_t), POINTER :: KMMatrix, MMatrix
          TYPE(ValueList_t), POINTER :: SolverParams
          REAL(KIND=dp),DIMENSION(n) :: y2
          REAL(KIND=dp) :: Crange,gamma,std
          REAL(KIND=dp) :: Norm
          INTEGER :: iter
          INTEGER :: Cm
          LOGICAL :: Parallel

          Parallel=(ParEnv % PEs > 1)

          SolverParams => GetSolverParams(Solver)

          Cm = ListGetInteger(SolverParams,"Matern exponent m",UnFoundFatal=.TRUE.)
          IF (Cm.LT.2) &
            CALL FATAL("Covariance","<Matern exponent m> should be >=2")
          Crange = ListGetConstReal(SolverParams,"correlation range", UnFoundFatal=.TRUE.)
          std = ListGetConstReal(SolverParams,"standard deviation",UnFoundFatal=.TRUE.) 
          gamma=sqrt(4*Pi*(Cm-1))*Crange

          MMatrix => MSolver % Matrix
          KMMatrix => KMSolver % Matrix

          ! Sigma x
          y2(1:n)=std*x(1:n)

          ! y_0=Gamma . y2
          y(1:n)=gamma*y2(1:n)

           ! L_M . M-1 . y0
          ! iter=1,m : (K+M) . y_i = (M) . y_{i-1}
          DO iter=1,Cm

          ! First iter M.M-1; nothing to do
           IF (iter.EQ.1) THEN
             KMMatrix % RHS(1:n) = y(1:n)
           ELSE
            IF (Parallel) THEN
              CALL ParallelInitSolve( MMatrix,y,MMatrix%RHS,KMMatrix % RHS )
              CALL ParallelMatrixVector( MMatrix,y, KMMatrix % RHS ,.TRUE. )
            ELSE
              CALL MatrixVectorMultiply( MMatrix, y, KMMatrix % RHS )
            END IF
           END IF
           ! And finally, solve:
           !--------------------
           Norm = DefaultSolve(USolver=KMSolver)
           y(1:n) = KMSolver % Variable % Values(1:n)

         END DO

         ! Gamma . y
         y2(1:n)=gamma*y(1:n)

         ! Sigma x
         y(1:n)=std*y2(1:n)

         END SUBROUTINE CovarianceVectorMultiplyD
!-------------------------------------------------------------------------------------------
!  Square root: C=V.V^T => y = V.x
!-------------------------------------------------------------------------------------------
         SUBROUTINE SqrCovarianceVectorMultiplyD(Solver,MSolver,KMSolver,n,x,y)
          TYPE(Solver_t) :: Solver
          TYPE(Solver_t), POINTER :: MSolver,KMSolver
          INTEGER,INTENT(IN) :: n
          REAL(KIND=dp),DIMENSION(n),INTENT(IN) :: x
          REAL(KIND=dp),DIMENSION(n),INTENT(OUT) :: y

          TYPE(Matrix_t), POINTER :: KMMatrix, MMatrix
          TYPE(ValueList_t), POINTER :: SolverParams
          REAL(KIND=dp),DIMENSION(n) :: ML
          REAL(KIND=dp),DIMENSION(n) :: y2
          REAL(KIND=dp) :: Crange,gamma,std
          REAL(KIND=dp) :: Norm
          REAL(KIND=dp) :: msum
          INTEGER :: iter
          INTEGER :: Cm,kk
          LOGICAL :: Parallel
          INTEGER :: i,j

          Parallel=(ParEnv % PEs > 1)
          IF (Parallel) &
             CALL FATAL("SqrCovarianceVector","not ready for parallel")

          SolverParams => GetSolverParams(Solver)

          Cm = ListGetInteger(SolverParams,"Matern exponent m",UnFoundFatal=.TRUE.)
          IF (Cm.LT.2) &
            CALL FATAL("Covariance","<Matern exponent m> should be >=2")
          IF (mod(Cm,2).NE.0) &
             CALL FATAL("SqrCovarianceVector","m should be even")

          Crange = ListGetConstReal(SolverParams,"correlation range", UnFoundFatal=.TRUE.)
          std = ListGetConstReal(SolverParams,"standard deviation",UnFoundFatal=.TRUE.)
          gamma=sqrt(4*Pi*(Cm-1))*Crange

          MMatrix => MSolver % Matrix
          KMMatrix => KMSolver % Matrix

          ! lumped mass matrix; this could be done only once..
          DO i=1,n
            msum = 0.0_dp
            DO j=MMatrix%Rows(i),MMatrix % Rows(i+1)-1
              msum = msum + MMatrix % Values(j)
            END DO
            ML(i) = msum
          END DO

           ! M^-1/2.x
           ML(1:n)=sqrt(ML(1:n))
           y(1:n)=x(1:n)/ML(1:n)

          kk=Cm/2
          DO iter=1,kk

            CALL MatrixVectorMultiply( MMatrix, y, KMMatrix % RHS )

           ! And finally, solve:
           !--------------------
           Norm = DefaultSolve(USolver=KMSolver)
           y(1:n) = KMSolver % Variable % Values(1:n)

         END DO

         ! Gamma . y
         y2(1:n)=gamma*y(1:n)

         ! Sigma x
         y(1:n)=std*y2(1:n)


         END SUBROUTINE SqrCovarianceVectorMultiplyD
!-------------------------------------------------------------------------------------------
!  InvCovariance : y = C^-1 x
!-------------------------------------------------------------------------------------------
        SUBROUTINE InvCovarianceVectorMultiplyD(Solver,MSolver,KMSolver,n,x,y)
          TYPE(Solver_t) :: Solver
          TYPE(Solver_t), POINTER :: MSolver,KMSolver
          INTEGER,INTENT(IN) :: n
          REAL(KIND=dp),DIMENSION(n),INTENT(IN) :: x
          REAL(KIND=dp),DIMENSION(n),INTENT(OUT) :: y

          TYPE(Matrix_t), POINTER :: KMMatrix, MMatrix
          TYPE(ValueList_t), POINTER :: SolverParams
          REAL(KIND=dp),DIMENSION(n) :: y2
          REAL(KIND=dp) :: Crange,gamma,std
          REAL(KIND=dp) :: Norm
          INTEGER :: iter
          INTEGER :: Cm
          LOGICAL :: Parallel

          Parallel=(ParEnv % PEs > 1)

          SolverParams => GetSolverParams(Solver)

          Cm = ListGetInteger(SolverParams,"Matern exponent m",UnFoundFatal=.TRUE.)
          IF (Cm.LT.2) &
            CALL FATAL("Covariance","<Matern exponent m> should be >=2")
          Crange = ListGetConstReal(SolverParams,"correlation range", UnFoundFatal=.TRUE.)
          std = ListGetConstReal(SolverParams,"standard deviation",UnFoundFatal=.TRUE.) 
          gamma=sqrt(4*Pi*(Cm-1))*Crange

          MMatrix => MSolver % Matrix
          KMMatrix => KMSolver % Matrix

          ! Sigma-1 x
          y2(1:n)=x(1:n)/std

          ! y_0=Gamma^-1 . y2
          y(1:n)=y2(1:n)/gamma

          ! L_M^-1 . y0
          ! iter=1,m : M . y_i = (M+K) . y_{i-1}
          DO iter=1,Cm

            IF (Parallel) THEN
              CALL ParallelInitSolve( KMMatrix,y,KMMatrix%RHS,MMatrix % RHS )
              CALL ParallelMatrixVector( KMMatrix,y, MMatrix % RHS ,.TRUE. )
            ELSE
              CALL MatrixVectorMultiply( KMMatrix, y, MMatrix % RHS )
            END IF
           ! And finally, solve:
           !--------------------
           Norm = DefaultSolve(USolver=MSolver)
           y(1:n) = MSolver % Variable % Values(1:n)

         END DO

         ! y2=M . y_m
         IF (Parallel) THEN
           CALL ParallelInitSolve( MMatrix,y,MMatrix % RHS,y2)
           CALL ParallelMatrixVector( MMatrix,y, y2,.TRUE.)
           CALL ParallelSumVector( MMatrix, y2 )
         ELSE
           CALL MatrixVectorMultiply( MMatrix, y, y2 )
         ENDIF

         ! Gamma-1 . y2
         y(1:n)=y2(1:n)/gamma

         ! Sigma-1 . y
         y(1:n)=y(1:n)/std

        END SUBROUTINE InvCovarianceVectorMultiplyD

!############################################################################################
!##
!## CovarianceInit subroutines to initialise the correlation matrices
!##
!############################################################################################
!-------------------------------------------------------------------------------------------
!  Using the diffusion operator approach; build the  required matrices
!   
!-------------------------------------------------------------------------------------------
        SUBROUTINE CovarianceInitD(Solver,MSolver,KMSolver)
          TYPE(Solver_t) :: Solver
          TYPE(Solver_t), POINTER :: MSolver,KMSolver

          TYPE(ValueList_t), POINTER :: SolverParams
          TYPE(Element_t),POINTER :: Element
          REAL(kind=dp) :: CRange
          LOGICAL :: Parallel
          INTEGER :: t,Active
          INTEGER :: n
          INTEGER :: ProjCoord
          LOGICAL :: DoProj

          Parallel=(ParEnv % PEs > 1)

          SolverParams => GetSolverParams(Solver)
          Crange = ListGetConstReal(SolverParams,"correlation range", UnFoundFatal=.TRUE.)
          ProjCoord = ListGetInteger(SolverParams,"projection coordinate",DoProj)

          MSolver =>  CreateChildSolver( Solver,TRIM(Solver%Variable%Name)//"_mvar",1 )
          KMSolver => CreateChildSolver( Solver,TRIM(Solver%Variable%Name)//"_kmVar",1)
          IF (Parallel) THEN
            MSolver % Parallel = .TRUE.
            KMSolver % Parallel = .TRUE.
          END IF
          MSolver % Variable % Output = .FALSE.
          KMSolver % Variable % Output = .FALSE.

          ! System assembly:
          !----------------
          CALL DefaultInitialize(USolver=MSolver)
          CALL DefaultInitialize(USolver=KMSolver)
          Active = GetNOFActive(Solver)
          DO t=1,Active
            Element => GetActiveElement(t)
            n  = GetElementNOFNodes()
            CALL LocalMatrix(  Element, n, Crange)
          END DO

          CALL DefaultFinishBulkAssembly(Solver=MSolver)
          CALL DefaultFinishBulkAssembly(Solver=KMSolver)

          CALL DefaultFinishAssembly(Solver=MSolver)
          CALL DefaultFinishAssembly(Solver=KMSolver)

          CONTAINS
            !------------------------------------------------------------------------------
            SUBROUTINE LocalMatrix( Element, n, l)
            !------------------------------------------------------------------------------
            INTEGER :: n
            TYPE(Element_t), POINTER :: Element
            REAL(KIND=dp) :: l
            INTEGER :: iter
            !------------------------------------------------------------------------------
            REAL(KIND=dp) :: Weight
            REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ
            REAL(KIND=dp) :: MASS(n,n), STIFF(n,n), FORCE(n)
            LOGICAL :: Stat
            INTEGER :: t,p,q
            TYPE(GaussIntegrationPoints_t) :: IP
            TYPE(Nodes_t) :: Nodes
            SAVE Nodes
            !------------------------------------------------------------------------------

            CALL GetElementNodes( Nodes )

            IF (DoProj) THEN
             SELECT CASE (ProjCoord)
                CASE (1)
                        Nodes % x(:)=0._dp
                CASE (2)
                        Nodes % y(:)=0._dp
                CASE (3)
                        Nodes % z(:)=0._dp
                CASE DEFAULT
                        CALL FATAL("CovarianceInitD", &
                            "Wrong projection coordinate"//I2S(ProjCoord))

             END SELECT
            END IF

            MASS  = 0._dp
            STIFF = 0._dp
            FORCE = 0._dp

            ! Numerical integration:
            !-----------------------
            IP = GaussPoints( Element )
            DO t=1,IP % n
              ! Basis function values & derivatives at the integration point:
              !--------------------------------------------------------------
              stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis, dBasisdx )

              Weight = IP % s(t) * DetJ

              ! diffusion term (l**2*grad(u),grad(v)):
              ! -----------------------------------
              STIFF(1:n,1:n) = STIFF(1:n,1:n) + Weight * &
                l * l * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

              DO p=1,n
                DO q=1,n
                ! identity term (R*u,v)
                ! -----------------------------------
                  MASS(p,q) = MASS(p,q) + Weight * Basis(q) * Basis(p)
                END DO
              END DO

            END DO

            CALL DefaultUpdateEquations(STIFF+MASS,FORCE,USolver=KMSolver)
            CALL DefaultUpdateEquations(MASS,FORCE,USolver=MSolver)
            !------------------------------------------------------------------------------
            END SUBROUTINE LocalMatrix
            !------------------------------------------------------------------------------
        END SUBROUTINE CovarianceInitD

!-------------------------------------------------------------------------------------------
!  Initialisation of the full covariance matrix in lapack uplo format
!   restricted to small 1D/2D serial cases....
!-------------------------------------------------------------------------------------------
          SUBROUTINE CovarianceInitL(Solver,n,InvPerm,aap,Op,PbDim)
          TYPE(Solver_t) :: Solver
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: InvPerm(n)
          REAL(kind=dp),INTENT(OUT) :: aap(:)  
          INTEGER,INTENT(IN) :: Op !1: correlation matrix ; 2: Cholesky; 3: inverse
          INTEGER,INTENT(IN) :: PbDim

          TYPE(ValueList_t), POINTER :: SolverParams
          REAL(kind=dp) :: x(n),y(n)
          REAL(kind=dp) :: d,dx,dy
          INTEGER :: i,j,kk
          INTEGER :: infoo

          CHARACTER(LEN=MAX_NAME_LEN) :: Ctype
          REAL(KIND=dp) :: Crange
          INTEGER :: p

          CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="CovarianceInit"
          character*1,parameter :: uplo='L'

          LOGICAL :: Parallel

          ! Current restrictions:
          ! - serial
          Parallel=(ParEnv % PEs > 1)
          IF (Parallel) &
             CALL FATAL(SolverName,'Sorry serial only!')

          ! Get parameter remated to the correlation function
          SolverParams => GetSolverParams(Solver)

          Ctype = ListGetString(SolverParams,"correlation type",UnFoundFatal=.TRUE.)

          IF (Ctype == 'maternp') THEN
              p = ListGetInteger(SolverParams,"MaternP polynomial order",UnFoundFatal=.TRUE.)
          ELSE IF (Ctype == 'materni') THEN
              p = ListGetInteger(SolverParams,"MaternI order",UnFoundFatal=.TRUE.)
          ENDIF

          Crange = ListGetConstReal(SolverParams,"correlation range", UnFoundFatal=.TRUE.)

          ! build the correlation matrix
          CALL INFO(SolverName,"Computing correlation matrix",level=4)

          DO i=1,n
            d=0._dp
            DO j=1,i
              dx=Solver%Mesh%Nodes%x(InvPerm(i))-Solver%Mesh%Nodes%x(InvPerm(j))
              d=dx**2
              IF (PbDim.GE.2) THEN
                dx=Solver%Mesh%Nodes%y(InvPerm(i))-Solver%Mesh%Nodes%y(InvPerm(j))
                d=d+dx**2
                IF (PbDim.EQ.3) THEN
                  dx=Solver%Mesh%Nodes%z(InvPerm(i))-Solver%Mesh%Nodes%z(InvPerm(j))
                  d=d+dx**2
                END IF
              END IF
              d=sqrt(d)
              aap(i + (j-1)*(2*n-j)/2) = correlation(d,Ctype,Crange,p)
            END DO
          END DO

          IF ( Op == 1 ) RETURN

          ! get Cholesky decomposition A=LL^T
          CALL INFO(SolverName,"Computing Cholesky",level=4)
          call dpptrf(uplo,n,aap,infoo)
          IF (infoo.NE.0) THEN
             CALL FATAL(SolverName,'Cholesky Factorisation failed')
          END IF

          IF ( Op == 2 ) RETURN

          ! compute the inverse from Cholesky
          CALL INFO(SolverName,"Computing Inverse",level=4)
          call dpptri(uplo,n,aap,infoo)

          IF (infoo.NE.0) THEN
            CALL FATAL(SolverName,'Inversion failed')
          END IF

        END SUBROUTINE CovarianceInitL


!####################################################
!# Some classical analytical correlation functions
!####################################################
        function correlation(d,Ctype,r,p) RESULT(c)
         implicit none
         real(kind=dp) :: c ! the correlation value
         real(kind=dp) :: d ! the distance
         real(kind=dp) :: r ! the range
         integer :: p       ! Matern polynomial order, for nu=p+1/2 for maternP
                            ! or matern oder nu, for maternI
         CHARACTER(LEN=MAX_NAME_LEN) :: Ctype   ! corrlation function type:

         SELECT CASE (Ctype)

           CASE ('exponential')
             c=exp(-d/r)

           CASE ('gaussian')
             c=exp(-d**2/(2*r**2))

          ! Matern case with integer order nu=p
          CASE ('materni')
             c=MaternI(p,d/r)

          ! Matern case with nu half integer p=nu-1/2
           CASE ('maternp')
             SELECT CASE (p)
               CASE(1)
                 c=(1+d/r)*exp(-d/r)
               CASE(2)
                 c=(1+d/r+d**2/(3*r**2))*exp(-d/r)
               CASE DEFAULT
                 CALL FATAL('correlation',&
                           'Matern polynomial order not available')
             END SELECT


           CASE DEFAULT
             CALL FATAL('correlation',&
                        'Correlation Type not known')
         END SELECT

        end function correlation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The matern correlation function for nu=integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION MaternI(n,x) RESULT(c)
        REAL(KIND=dp) :: c
        REAL(kind=dp),INTENT(IN) :: x
        INTEGER,INTENT(IN) :: n
        real(kind=dp) :: Kn

        if (n.lt.1) CALL FATAL("MaternI","n must be >=1")
        ! Bessel not defined for x=0; but c should be 1
        if (x.lt.10*AEPS) THEN
          c=1._dp
          Return
        endif
       if (n.gt.1) then
         Kn=BesselKn(n,x)
       else
         Kn=BesselK1(x)
       endif
       c=(2.0**(1-n))*((x)**n)*Kn/fact(n-1)

      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified Bessel function of second kind Kn  
!  Computed from the recursion
!    Kn+1 = 2*n/x Kn + Kn-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      FUNCTION BesselKn(n,x) RESULT(Kn)
      REAL(kind=dp) :: Kn
      INTEGER,INTENT(IN) :: n
      REAL(kind=dp),INTENT(IN) :: x

      INTEGER :: j
      REAL(kind=dp) :: y
      REAL(kind=dp) :: bjm,bj,bjp

       IF (x.LT.AEPS) &
         CALL FATAL("BesselKn","x is too small")
       IF (n.lt.2) THEN
         CALL FATAL("BesselKn","n must be >= 2")
       END IF

       y=2._dp/x

       bjm=BesselK0(x)
       bj=BesselK1(x)

       DO j=1,n-1
         bjp=bjm+j*y*bj
         bjm=bj
         bj=bjp
       END DO

       Kn=bj

       RETURN
       END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified Bessel function of second kind of order 0, K0
!  Polynomial approximation; cf
!  https://www.advanpix.com/2015/11/25/rational-approximations-for-the-modified-bessel-function-of-the-second-kind-k0-for-computations-with-double-precision/       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      FUNCTION BesselK0(x) RESULT(K0)
      REAL(KIND=dp) :: K0
      REAL(KIND=dp), INTENT(IN) :: x

      REAL(KIND=dp), Dimension(8),Parameter :: P7=(/ 1.1593151565841244842077226e-01,&
        2.7898287891460317300886539e-01,&
        2.5248929932161220559969776e-02,&
        8.4603509072136578707676406e-04,&
        1.4914719243067801775856150e-05,&
        1.6271068931224552553548933e-07,&
        1.2082660336282566759313543e-09,&
        6.6117104672254184399933971e-12/)

      REAL(KIND=dp), Dimension(7),Parameter :: P6=(/ 1.0000000000000000044974165e+00,&
       2.4999999999999822316775454e-01,&
       2.7777777777892149148858521e-02,&
       1.7361111083544590676709592e-03,&
       6.9444476047072424198677755e-05,&
       1.9288265756466775034067979e-06,&
       3.9908220583262192851839992e-08/)

      REAL(KIND=dp), Dimension(3), Parameter :: Q2=(/ 8.5331186362410449871043129e-02,&
        7.3477344946182065340442326e-01,&
        1.4594189037511445958046540e+00/)

      REAL(KIND=dp), Dimension(22), Parameter :: P21=(/1.0694678222191263215918328e-01,&
        9.0753360415683846760792445e-01,&
        1.7215172959695072045669045e+00,&
        -1.7172089076875257095489749e-01,&
        7.3154750356991229825958019e-02,&
        -5.4975286232097852780866385e-02,&
        5.7217703802970844746230694e-02,&
        -7.2884177844363453190380429e-02,&
        1.0443967655783544973080767e-01,&
        -1.5741597553317349976818516e-01,&
        2.3582486699296814538802637e-01,&
        -3.3484166783257765115562496e-01,&
        4.3328524890855568555069622e-01,&
        -4.9470375304462431447923425e-01,&
        4.8474122247422388055091847e-01,&
        -3.9725799556374477699937953e-01,&
        2.6507653322930767914034592e-01,&
        -1.3951265948137254924254912e-01,&
        5.5500667358490463548729700e-02,&
        -1.5636955694760495736676521e-02,&
        2.7741514506299244078981715e-03,&
        -2.3261089001545715929104236e-04/)

      REAL(KIND=dp) :: x2,y,y2,p,p1,I0

      IF (x.lt.1._dp) THEN
        x2=x*x
        y=x/2.0_dp
        y2=y*y

        !P6[(x/2)^2]
        p1=Polynomial(y2,P6,6)
        !I0
        I0=1._dp+y2*p1

        !! P7(x2)
        p=Polynomial(x2,P7,7)

        K0=p-log(x)*I0

      ELSE

        y=1._dp/x
        !P21(x^-1)
        p=Polynomial(y,P21,21)
        !Q2(x^-1)
        p1=Polynomial(y,Q2,2)
        K0=(p/p1)
        K0=K0/(sqrt(x)*exp(x))

      ENDIF

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified Bessel function of second kind of order 1, K1
!  Polynomial approximation; cf
!    https://www.advanpix.com/2016/01/05/rational-approximations-for-the-modified-bessel-function-of-the-second-kind-k1-for-computations-with-double-precision/      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION BesselK1(x) RESULT(K1)
      REAL(KIND=dp) :: K1
      REAL(KIND=dp), INTENT(IN) :: x

      REAL(KIND=dp), Dimension(6), Parameter :: P5=(/8.3333333333333325191635191e-02,&
           6.9444444444467956461838830e-03,&
           3.4722222211230452695165215e-04,&
           1.1574075952009842696580084e-05,&
           2.7555870002088181016676934e-07,&
           4.9724386164128529514040614e-09/)

     REAL(KIND=dp), Dimension(9), Parameter :: P8=(/-3.0796575782920622440538935e-01,&
             -8.5370719728650778045782736e-02,&
             -4.6421827664715603298154971e-03,&
             -1.1253607036630425931072996e-04,&
             -1.5592887702110907110292728e-06,&
             -1.4030163679125934402498239e-08,&
             -8.8718998640336832196558868e-11,&
             -4.1614323580221539328960335e-13,&
             -1.5261293392975541707230366e-15/)

     REAL(KIND=dp), Dimension(23), Parameter :: P22=(/1.0234817795732426171122752e-01,&
              9.4576473594736724815742878e-01,&
              2.1876721356881381470401990e+00,&
              6.0143447861316538915034873e-01,&
              -1.3961391456741388991743381e-01,&
              8.8229427272346799004782764e-02,&
              -8.5494054051512748665954180e-02,&
              1.0617946033429943924055318e-01,&
              -1.5284482951051872048173726e-01,&
              2.3707700686462639842005570e-01,&
              -3.7345723872158017497895685e-01,&
              5.6874783855986054797640277e-01,&
              -8.0418742944483208700659463e-01,&
              1.0215105768084562101457969e+00,&
              -1.1342221242815914077805587e+00,&
              1.0746932686976675016706662e+00,&
              -8.4904532475797772009120500e-01,&
              5.4542251056566299656460363e-01,&
              -2.7630896752209862007904214e-01,&
              1.0585982409547307546052147e-01,&
              -2.8751691985417886721803220e-02,&
              4.9233441525877381700355793e-03,&
              -3.9900679319457222207987456e-04/)

      REAL(KIND=dp), Dimension(3), Parameter :: Q2=(/ 8.1662031018453173425764707e-02,&
              7.2398781933228355889996920e-01,&
              1.4835841581744134589980018e+00/)

      REAL(KIND=dp) :: y,y2,y4,x2,I1
      REAL(KIND=dp) :: p,p1

      IF (x.lt.1._dp) THEN
        x2=x*x
        y=x/2.0_dp
        y2=y*y
        y4=y2*y2

        !P5[(x/2)^2]
        p1=Polynomial(y2,P5,5)
        !I1
        I1=y*(1+y2/2+y4*p1)

        !! P8(x2)
        p=Polynomial(x2,P8,8)

        K1=log(x)*I1+1/x+x*p

      ELSE

        y=1._dp/x
        !P22(x^-1)
        p=Polynomial(y,P22,22)
        !Q2(x^-1)
        p1=Polynomial(y,Q2,2)
        K1=(p/p1)
        K1=K1/(sqrt(x)*exp(x))

      ENDIF

       RETURN
      END FUNCTION BesselK1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the value at x of a Polynom of order n
!   Horner scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION Polynomial(x,P,n) RESULT(y)
      IMPLICIT NONE
      REAL(KIND=dp) :: y
      REAL(KIND=dp),INTENT(IN) :: P(n+1),x
      INTEGER ,INTENT(IN) :: n
      INTEGER :: i

      y=P(n+1)
      Do i=1,n
        y=P(n+1-i)+y*x
      End Do

       RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Integer factorial 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function fact(n)
      implicit none
      integer :: fact
      integer, intent(IN) :: n
      integer :: i

      fact = 1
      do i = 2, n
        fact = fact * i
      enddo
      return
      end function fact


      END MODULE CovarianceUtils

