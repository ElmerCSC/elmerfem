
!------------------------------------------------------------------------------
!>  Calculate the charge density using the eigenvectors of variable 
!>  "Wavefunctions" of eigenproblem solver 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ChargeDensitySolver( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
 
  LOGICAL :: Found, FirstTime = .TRUE.
  INTEGER :: i, j, k, t, n, nd, m, v, istat
  INTEGER :: EigSolverNumber, NOFeigs, GlobSize
  INTEGER, ALLOCATABLE :: Indexes(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: str_tmp
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp) :: Norm, s, s0
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
  REAL(KIND=dp), POINTER :: EigsWeights(:,:)
  TYPE(Valuelist_t), POINTER :: SolverParams
  SAVE NOFeigs, EigsWeights, EigSolverNumber, FirstTime, STIFF, FORCE, Indexes
!------------------------------------------------------------------------------

  CALL INFO('ChargeDensityCalculation','...........................')
  CALL INFO('ChargeDensityCalculation',' Charge density calculation')
  CALL INFO('ChargeDensityCalculation','...........................')

!------------------------------------------------------------------------------
  FirstTimeOnly: IF (FirstTime) THEN
    
    n = Solver % Mesh % MaxElementDOFs
    ALLOCATE( Indexes(n), STIFF(n,n), FORCE(n) )
    
    ! Reading the number of eigenvalues to be included and the
    ! weights of the eigen states from .sif file.
    !-----------------------------------------------------------
    SolverParams => GetSolverParams()
  
    NOFeigs = GetInteger( SolverParams,'Number of Eigenmodes Included',Found)
    IF(.NOT. Found) CALL ERROR('ChargeDensityCalculation', &
        'Number of Eigenmodes Included was not defined on the .sif file.')
    
    WRITE(Message,*) 'Number of Eigenmodes included was set to ',NOFeigs
    CALL INFO('ChargeDensityCalculation',Message)
    
    CALL GetConstRealArray(SolverParams, EigsWeights, &
        'Weights of Eigen States', Found)
    IF(.NOT. Found) THEN
       ALLOCATE( EigsWeights( NOFeigs, 1 ) , STAT=istat )
       IF (istat /= 0) CALL ERROR('ChargeDensityCalculation', &
           'Error when allocating memory for weights of the eigen states.')
    
      EigsWeights = 1.0d0
      WRITE(Message,*)  'Weights of Eigen States was not defined on', &
          'the .sif file. All weights were set to 1 .'
      CALL WARN('ChargeDensityCalculation', Message)
    END IF
    
    WRITE(Message,*) 'Weights of the Eigen States are', &
        EigsWeights( 1:NOFeigs, 1 )
    CALL INFO('ChargeDensityCalculation', Message)
    
    ! The number of solver that has the eigenvalues as a variable
    ! is found out.
    !--------------------------------------------------------------
    EigSolverNumber = 0
    DO k = 1, Model % NumberOfSolvers
      
      SolverParams => Model % Solvers(k) % Values
      IF (GetLogical( SolverParams, 'Eigen Analysis', Found)) THEN
        IF ( EigSolverNumber == 0 ) THEN
          EigSolverNumber = k
        ELSE
          WRITE(Message,*) 'EigenSolver coud not be identified.', &
              'More than one solver has Eigen Anylysis = True'
          CALL ERROR('ChargeDensityCalculation', Message)
        END IF
        
      END IF
    END DO
    
    WRITE(Message,*) 'Wavefunctions are variables of solver number ', EigSolverNumber
    CALL INFO('ChargeDensityCalculation', Message)

    FirstTime = .FALSE.
    
  END IF FirstTimeOnly
!-----------------------------------------------------------------------------

  ! Calculate the ChargeDensity from the solution of eigenproblem
  !--------------------------------------------------------------
  Solver % Variable % Values = 0.0d0
  CALL DefaultInitialize()
  
  DO t = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement( t )
    n  = GetElementNOFNodes()
    nd = GetElementDOFs( Indexes )
    Indexes(1:nd) = Model % Solvers(EigSolverNumber) % &
          Variable % Perm(Indexes(1:nd))
    
    ! Get element local stiffness and mass matrices:
    !-----------------------------------------------
    CALL LocalMatrix(  STIFF, FORCE, Element, n, nd )
    
    !Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()
  
  Norm = DefaultSolve()

  ! Printing some values of the Charge Density
  !-------------------------------------------------------------------
  n = Solver % Mesh % NumberOfNodes

  s = MAXVAL(Solver % Variable % Values(Solver % Variable % Perm(1:n)))
! IF ( ParEnv % Pes > 1 ) THEN
!   CALL MPI_ALLREDUCE(s,s0,1,MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD, i)
!   s = s0
! END IF
  WRITE(Message,*) 'Greatest value of the Charge density is ', s
  CALL INFO('ChargeDensityCalculation', Message)
  
  s = MINVAL(Solver % Variable % Values(Solver % Variable % Perm(1:n)))
! IF ( ParEnv % Pes > 1 ) THEN
!   CALL MPI_ALLREDUCE(s,s0,1,MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, i)
!   s = s0
! END IF

  WRITE(Message,*) 'Smallest value of the Charge density is ', s
  CALL INFO('ChargeDensityCalculation', Message)

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF,  FORCE, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd), DetJ
    REAL(KIND=dp) :: LocalPotentialAtIP, LocalChargeAtIP
    LOGICAL :: Stat
    INTEGER :: t
    TYPE(GaussIntegrationPoints_t) :: IP
    
    INTEGER :: i, j, k, dim, xctype
    REAL(KIND=dp) :: u, v, w
    
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    
    DO t = 1,IP % n
      u = IP % u(t)
      v = IP % v(t)
      w = IP % w(t)

      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
      
      
      ! The Coulomb Potential term and the Charge density
      ! at the integration point:
      !------------------------------------------------------------
      LocalChargeAtIp = 0.0d0
      DO i = 1, NOFeigs    
        LocalChargeAtIP = LocalChargeAtIP + EigsWeights(i,1) * &
          SUM( Model % Solvers(EigSolverNumber) %  &
            Variable % EigenVectors(i,Indexes(1:nd)) * Basis(1:nd) )**2
      END DO
      
      ! Calculating stiffness and mass matrices:
      !------------------------------------------
      DO i = 1,nd
        DO j = 1,nd
          STIFF(i,j) = STIFF(i,j) + Basis(i) * Basis(j) * detJ * IP % s(t)
        END DO
        FORCE(i) = FORCE(i) + LocalChargeAtIp * Basis(i) * detJ * IP % s(t)
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
END SUBROUTINE ChargeDensitySolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Solve the eigenproblem -(1/2) * Lap(phi_i) + U * phi_i = e_i * phi_i
!>  of Kohn Sham equations. Calculating the XC potential here!
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE WaveFunctionSolver( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE ExchangeCorrelations

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE( Element_t ),POINTER :: Element
  LOGICAL :: AllocationsDone = .FALSE., Found
  INTEGER :: n, nd, t, istat, GlobSize, xctype, NOFeigs
  REAL(KIND=dp) :: Norm, ShiftConst, EigMassNorm, TotalNucleiCharge,TotalCharge
  CHARACTER(LEN=MAX_NAME_LEN) :: str_tmp
  TYPE(Valuelist_t), POINTER :: SolverParams
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), MASS(:,:), FORCE(:)
  REAL(KIND=dp), ALLOCATABLE :: LocalHartreePotential(:), LocalChargeDensity(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: EigStorage(:)
  
  SAVE STIFF, MASS, FORCE, LocalHartreePotential, ShiftConst, &
       LocalChargeDensity, AllocationsDone, xctype, NOFeigs, &
       EigStorage
!------------------------------------------------------------------------------

  CALL INFO('WavefunctionSolver','............')
  CALL INFO('WavefunctionSolver','EigenSolver')
  CALL INFO('WavefunctionSolver','............')

  ! Allocate some permanent storage and define the shift parameter, 
  ! this is done first time only:
  !--------------------------------------------------------------
  FirstTimeOnly: IF ( .NOT. AllocationsDone ) THEN
     
     ! Allocate the mamory for local stiffness and mass matrices
     !----------------------------------------------------------
    N = Solver % Mesh % MaxElementDOFs
    ALLOCATE( STIFF(N,N), MASS(N,N), FORCE(N), &
        LocalHartreePotential(N), LocalChargeDensity(N), &
        STAT=istat )
    IF ( istat /= 0 ) &
        CALL Fatal( 'WavefunctionSolver', 'Memory allocation error.' )
 
    ! Read the type of correlation potential from the .sif file
    !-----------------------------------------------------------
    ! The integer parameter for calling the function uxc from module xc
    
    ! perdew_zunger=0
    ! von_barth_hedin=1, 
    ! gunnarsson_lundqvist=2
    ! perdew_wang=3
    
    SolverParams => GetSolverParams()
    
    xcType = -1
    
    str_tmp = GetString( SolverParams, 'XC Potential type', Found)
    
    IF (TRIM(str_tmp) == 'perdew-zunger') THEN
      IF (xcType /= -1 ) CALL WARN('WavefunctionSolver',&
          'Many different XC-potential types were set in .sif file.')
      IF (found .AND. xcType == -1 ) xcType=0
    END IF
    
    IF (TRIM(str_tmp) == 'von barth-hedin') THEN        
      IF (xcType /= -1 ) CALL WARN('WavefunctionSolver',&
          'Many different XC-potential types were set in .sif file.')
      IF (found .AND. xcType == -1 ) xcType=1
    END IF
    
    IF (TRIM(str_tmp) == 'gunnarsson-lundqvist') THEN        
      IF (xcType /= -1 ) CALL WARN('WavefunctionSolver',&
          'Many different XC-potential types were set in .sif file.')
      IF (found .AND. xcType == -1 ) xcType=2
    END IF
    
    IF (TRIM(str_tmp) == 'perdew-wang') THEN        
      IF (xcType /= -1 ) CALL WARN('WavefunctionSolver',&
          'Many different XC-potential types were set in .sif file.')
      IF (found .AND. xcType == -1 ) xcType=3
    END IF
    
    IF (TRIM(str_tmp) == 'none') THEN        
      IF (xcType /= -1 ) CALL WARN('WavefunctionSolver',&
          'Many different XC-potential types were set in .sif file.')
      IF (found .AND. xcType == -1 ) xcType=4
    END IF
    
    IF (xcType == -1 ) CALL ERROR('WavefunctionSolver',&
        'No XC Potential type was set on the .sif file')
    
    IF(xcType == 0) &
        CALL INFO('WavefunctionSolver', 'Using Perdew-Zunger type XC-potential.')
    IF(xcType == 1) &
        CALL INFO('WavefunctionSolver', 'Using Von Barth-Hedin type XC-potential.')
    IF(xcType == 2) &
        CALL INFO('WavefunctionSolver', 'Using Gunnarsson-Lundqvist type XC-potential.')
    IF(xcType == 3) &
        CALL INFO('WavefunctionSolver', 'Using Perdew-Wang type XC-potential.')
    IF(xcType == 4) &
        CALL INFO('WavefunctionSolver', 'Not using any type of XC-potential.')
    
    
     ! Read the Shifting parameter from the .sif file.
     !-----------------------------------------------------------
    ShiftConst = 0.0d0
    ShiftConst = GetConstReal( SolverParams , &
        'Eigen System Shifting Constant', Found )
    IF (.NOT. Found) THEN
      CALL INFO('WavefunctionSolver', &
          'Eigen System Shifting constant not set. ')
    ELSE
      PRINT*,'Eigen System Shifting constant was set to', ShiftConst        
    END IF

    ! Allocate memory for vector storing the previous eigenvalues.
    !--------------------------------------------------------------
    NOFeigs = GetInteger( SolverParams, 'Eigen System Values', Found)
    IF (.NOT. Found) &
        CALL ERROR('WavefunctionSolver', 'Eigen System Values was not set. ')
    
    ALLOCATE( EigStorage(NOFeigs) ,  STAT=istat )
    IF ( istat /= 0 ) &
        CALL Fatal( 'WavefunctionSolver', 'Memory allocation error.' )
    
    CALL INFO('WavefunctionSolver', &
        'Allocations for WavefunctionSolver are done.' )
    
    AllocationsDone = .TRUE.
  END IF FirstTimeOnly
  
!----------------------------------------------------------------------------

  !Initialize the system and do the assembly:
  !------------------------------------------

  CALL DefaultInitialize()
  LocalHartreePotential = 0.0d0
  LocalChargeDensity = 0.0d0
  FORCE = 0.0d0
  
  DO t = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement( t )
    n  = GetElementNOFNodes( Element )
    nd = GetElementNOFDOFs( Element )
    
    ! Get potential and charge density in node points:
    !--------------------------------------------
    CALL GetScalarLocalSolution( LocalHartreePotential , &
        'Potential')
    CALL GetScalarLocalSolution( LocalChargeDensity , &
        'Charge density')
    
    ! Get element local stiffness and mass matrices:
    !-----------------------------------------------
    CALL LocalMatrix(  STIFF, MASS, Element, n, nd, &
        LocalHartreePotential, LocalChargeDensity, xctype )
    
    !Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL DefaultUpdateEquations( STIFF, FORCE )
    CALL DefaultUpdateMass( MASS )
    
  END DO
  
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()
  
  ! Applying the shift:
  !----------------------
  IF (.NOT. ABS(ShiftConst) < EPSILON(ShiftConst)) &
      Solver % Matrix % Values = Solver % Matrix % Values &
      + ShiftConst * Solver % Matrix % MassValues
  
  ! Store the previous eigenvalues:
  !---------------------------------
  EigStorage(1:NOFeigs) = Solver % Variable % EigenValues(1:NOFeigs)
  
  
  ! Solving the Eigensystem:
  !----------------------------------------------
  Norm = DefaultSolve()
  
  ! Print the absolute change of the eigenvalues:
  !-----------------------------------------------
  CALL INFO('WavefunctionSolver','Absolute change of the eigenvalues:')
  
  DO t = 1, NOFeigs
    WRITE(Message,*) 'Eigenvalue',t,'has absolute change', &
        ABS( Solver % Variable % EigenValues(t) - EigStorage(t) )
    CALL INFO('WavefunctionSolver',Message)
  END DO
  
CONTAINS
     
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, MASS, Element, n, nd, LocalPotential, &
       LocalChargeDensity, xctype )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), MASS(:,:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3), DetJ
    REAL(KIND=dp) :: LocalPotential(:), LocalChargeDensity(:), &
        LocalPotentialAtIP, LocalChargeAtIP, LocalXCPotentialAtIP
    LOGICAL :: Stat
    INTEGER :: t
    TYPE(GaussIntegrationPoints_t) :: IP
    
    INTEGER :: i, j, k, dim, xctype
    REAL(KIND=dp) :: u, v, w
    
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    MASS  = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    
    DO t = 1,IP % n
      u = IP % u(t)
      v = IP % v(t)
      w = IP % w(t)

      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, u, v, w, &
                  detJ, Basis, dBasisdx )
      
      
      ! The Coulomb Potential term and the Charge density
      ! at the integration point:
      !------------------------------------------------------------
      LocalPotentialAtIP = SUM( Basis(1:nd) * LocalPotential(1:nd) ) 
      LocalChargeAtIP    = SUM( Basis(1:nd) * LocalChargeDensity(1:nd) )
      LocalChargeAtIP    = MAX( LocalChargeAtIP, 0.0_dp )
      
      ! Calculating the XC-potential at the integration point:
      !------------------------------------------------------------
      IF (.NOT. xctype == 4) THEN
        LocalXCPotentialAtIP = uxc( LocalChargeAtIP, 0.0d0, 1, xcType )
      ELSE
        LocalXCPotentialAtIP = 0.0d0
      END IF
      
      ! Calculating stiffness and mass matrices:
      !------------------------------------------
      DO i = 1,nd
        DO j = 1,nd
          DO k = 1,dim
            
            ! Half times negative Laplacian
            STIFF(i,j) = STIFF(i,j) &
                + 0.5 * dBasisdx(i,k) * dBasisdx(j,k) * detJ * IP % s(t)
          END DO
          
          ! Mass matrix on the right hand side
          MASS(i,j) = MASS(i,j) + Basis(i) * Basis(j) * detJ * IP % s(t)
          
          ! Multiplication by total potential
          STIFF(i,j) = STIFF(i,j) &
              + ( LocalPotentialAtIP + LocalXCPotentialAtIP ) &
              * Basis(i) * Basis(j) * detJ * IP % s(t)
          
        END DO
      END DO
    END DO
    
    ! Forcing the symmetry:
    !-----------------------
    STIFF = 0.5d0 * ( STIFF + TRANSPOSE(STIFF) )
    MASS  = 0.5d0 * ( MASS  + TRANSPOSE(MASS)  )
    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE WaveFunctionSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Solve the (Hartree) potential via Poisson equation using the 
!>  variable "Charge density" multiplied by 4 Pi as a load.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE PoissonSolver( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE LinearAlgebra

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation  
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: Found, FirstTime = .TRUE.
  REAL(KIND=dp) :: Norm, TotalCharge, TotalNucleiCharge
  TYPE(Valuelist_t), POINTER :: SolverParams
  
  ! Variables for dealing with nuclei
  INTEGER :: n, nd, t, i, j, k, l, m, v, NOFnuclei, istat, permsize, &
      ChargeDensitySolverNumber
  REAL(KIND=dp), POINTER :: NucleiTable(:,:)
  INTEGER, POINTER :: NodesOfNuclei(:), Perm(:)
  REAL(KIND=dp) :: x, y, z, tol, r, s
  CHARACTER(LEN=MAX_NAME_LEN) :: str_tmp
  
  ! Variables for assembly
  TYPE( ValueList_t ), POINTER :: BodyForce
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:),  FORCE(:), LocalChargeDensity(:)
  
  ! Variables for relaxation
  LOGICAL :: CalculatingRealSolution = .FALSE., StoreOnly = .TRUE.  
  REAL(KIND=dp), ALLOCATABLE :: OldSolution(:)
  INTEGER :: RelaxType, GlobSize, StartGRPulayAfter
  REAL(kind=dp) :: ParamA, ParamB, ParamC, RF, NOFiter, &
      StartGRPulayIfRFis, temp

  TYPE(Variable_t), POINTER :: EigsVar
  
  SAVE STIFF, FORCE, LocalChargeDensity, &
      FirstTime, NucleiTable, NodesOfNuclei, NOFnuclei, Perm, &
      TotalNucleiCharge, ChargeDensitySolverNumber, &
      OldSolution, RelaxType, GlobSize, ParamA, &
      ParamB, ParamC, NOFiter, CalculatingRealSolution, &
      StoreOnly, StartGRPulayAfter, StartGRPulayIfRFis

!------------------------------------------------------------------------------

  CALL INFO('PoissonSolver','------------------')
  CALL INFO('PoissonSolver',' Poisson equation ')
  CALL INFO('PoissonSolver','------------------')

!------------------------------------------------------------------------------

  ! Local allocations and search for nuclei etc. are done:
  !=======================================================
  FirsTimeOnly: IF ( FirstTime ) THEN
    
    N = Solver % Mesh % MaxElementDOFs ! just big enough for elemental arrays
    ALLOCATE( FORCE(N), LocalChargeDensity(N), STIFF(N,N), STAT=istat )
    IF ( istat /= 0 ) CALL Fatal( 'PoissonSolver', &
          'Memory allocation error.' )

    ! The number of solver that has the charge density 
    ! (namely Number of Eigenmodes Included as parameter) as a variable
    ! is found out.
    !--------------------------------------------------------------
    ChargeDensitySolverNumber = 0
    DO k = 1, Model % NumberOfSolvers
      
      SolverParams => Model % Solvers(k) % Values
      v = GetInteger( SolverParams, 'Number of Eigenmodes Included', Found)
      
      IF(Found) THEN
        
        IF ( ChargeDensitySolverNumber == 0 ) THEN
          ChargeDensitySolverNumber = k
        ELSE
          WRITE(Message,*) 'ChargeDensitySolver could not be identified.', &
              'More than one solver has Number of Eigenmodes Included', & 
              'defined.'
          CALL ERROR('PoissonSolver', Message)
        END IF
        
      END IF
    END DO
    
    IF (ChargeDensitySolverNumber == 0) THEN
      WRITE(Message,*) 'ChargeDnsitySolver coud not be identified.', &
          'No solver has Number of Eigenmodes Included', & 
          'defined.'
      CALL FATAL('PoissonSolver', Message)
    END IF
    
    ! Set pointer to the current solution parameters:
    !-------------------------------------------------
    SolverParams => GetSolverParams()
    
    ! Set pointer to the permutations of the variable:
    !-------------------------------------------------
    Perm => Solver % Variable % Perm
    
    ! Get the nuclei table from a .sif -file or set the NOFnuclei 0.
    !----------------------------------------------------------------    
    NOFnuclei = 0
    NOFnuclei = GetInteger( SolverParams, &
        'NOFnuclei',Found )
    IF(.NOT. Found) CALL WARN('PoissonSolver', &
        'NOFnuclei was not set on the sif file.')
    
    WRITE(Message,*) 'NOFnuclei was set to ', NOFnuclei
    CALL INFO('PoissonSolver', Message)
    
    CreateNucleiArrays: IF ( .NOT. (NOFnuclei == 0) ) THEN
      
      ALLOCATE( NucleiTable( NOFnuclei , 4 ), &
          NodesOfNuclei ( NOFnuclei ) )
      
      NucleiTable = 0.0d0
      NodesOfNuclei = 0
      
      CALL GetConstRealArray(SolverParams, NucleiTable, &
          'NucleiTable', Found)
      IF (.NOT. Found) CALL WARN('PoissonSolver', &
          'NucleiTable was not found from .sif file.')
      
      x = 0.0d0
      y = 0.0d0
      z = 0.0d0
      
      tol = 10*EPSILON(x)
      
      ! Checking all nodes
      DO k = 1, Solver % Mesh % NumberOfNodes 
!       IF ( ParEnv % PEs > 1 ) THEN
!         IF ( Solver % Mesh % ParallelInfo % NeighBourList(k) % &
!               Neighbours(1) /= ParEnv % mype ) CYCLE
!       END IF

        ! Coordinates of node
        x = Solver % Mesh % Nodes % x(k)
        y = Solver % Mesh % Nodes % y(k)
        z = Solver % Mesh % Nodes % z(k)
        
        DO m = 1, NOFnuclei
          IF ( ABS(NucleiTable(m,2)-x)<tol .AND. &
               ABS(NucleiTable(m,3)-y)<tol .AND. &
               ABS(NucleiTable(m,4)-z)<tol ) NodesOfNuclei(m) = k
        END DO
        
      END DO
    END IF CreateNucleiArrays
    
    ! Check that all nodes were found.
    !--------------------------------------
    
!   DO i = 1, NOFnuclei
!     IF( NodesOfNuclei(i) == 0) THEN           
!       WRITE(Message,*) 'Nuclei number ',i,'was not located on node.'
!       CALL ERROR('PoissonSolver', Message)
!     END IF
!   END DO
    
    WRITE(Message,*) 'Nuclei are located on nodes ', NodesOfNuclei
    CALL INFO('PoissonSolver', Message)
    
    ! Calculate the total nuclei charge.
    !----------------------------------------
    TotalNucleiCharge = 0.0d0
    DO v = 1 , NOFnuclei                      
      TotalNucleiCharge = TotalNucleiCharge + NucleiTable(v,1)
    END DO
    
    ! Allocating memory for storing the previous solution
    !-----------------------------------------------------------
    GlobSize = SIZE(Solver % Variable % Values)
    ALLOCATE( OldSolution(GlobSize) , STAT=istat )
    IF (istat /= 0) CALL ERROR('PoissonSolver', &
        'Error when allocating memory for previous solution strorage.')
    
    OldSolution = 0.0d0
    
    ! Reading the relaxation method and parameters from the .sif file
    !===================================================================
    ! RelaxType = 0 , No mixing
    ! Relaxtype = 1 , Constant mixing parameter
    ! RelaxType = 2 , Exponentially growing mixing parameter
    
    RelaxType = 0
    ParamA = -100000.0d0
    ParamB = -100000.0d0
    ParamC = -100000.0d0
    
    NOFiter = 0
    
    str_tmp = GetString( SolverParams, 'Relaxation Method', Found)
    IF( .NOT. Found) CALL WARN('PoissonSolver', &
        'No Relaxation method was found the .sif file, using none.')
    
    ParamA = GetConstReal( SolverParams, &
        'Relaxation Parameter A', Found )
    IF(.NOT. Found) CALL INFO('PoissonSolver', &
        'Relaxation Parameter A was not defined' // & 
        'on the .sif file.')
    
    ParamB = GetConstReal( SolverParams, &
        'Relaxation Parameter B', Found )
    IF(.NOT. Found) CALL INFO('PoissonSolver', &
        'Relaxation Parameter B was not defined' // & 
        'on the .sif file.')
    
    ParamC = GetConstReal( SolverParams, &
        'Relaxation Parameter C', Found )
    IF(.NOT. Found) CALL INFO('PoissonSolver', &
        'Relaxation Parameter C was not defined' // & 
        'on the .sif file.')
    
    IF (TRIM(str_tmp) == 'constant mixing') THEN
      RelaxType = 1
      WRITE(Message,*) 'Mixing Relaxation Factor (RF) is constant', &
          ParamA
      CALL INFO('PoissonSolver', Message)
      CALL INFO('PoissonSolver', &
          'Relaxation method is r(k+1) = (1 - RF) * r(k) + RF * r(k+1)') 
    END IF
    
    IF (TRIM(str_tmp) == 'exponential mixing') THEN
      RelaxType = 2
      WRITE(Message,*) 'Mixing Relaxation Factor (RF) is varying by scheme',&
          ' RF(k) = ',ParamC,' + 1 - ',ParamA,'* Exp( -',ParamB,' * k )'
      CALL INFO('PoissonSolver', Message)
      CALL INFO('PoissonSolver', &
          'Relaxation method is r(k+1) = (1 - RF) * r(k) + RF * r(k+1)') 
    END IF
    
    ! Read the parameters for starting GRPulay
    !---------------------------------------------
    StartGRPulayIfRFis = 1.0d8;
    StartGRPulayAfter = 1.0d8;
    
    StartGRPulayIfRFis = GetConstReal( SolverParams, &
        'Start GRPulay if relaxation factor is more than', Found )
    
    StartGRPulayAfter = GetInteger( SolverParams, &
        'Start GRPulay after iterations', Found )
    
    WRITE(Message,*) 'Starting GRPulay after ', StartGRPulayAfter,'iterations.' 
    CALL INFO('PoissonSolver', Message)
    
    WRITE(Message,*) 'Starting GRPulay when RF is more than ', &
        StartGRPulayIfRFis
    CALL INFO('PoissonSolver', Message)
    
    
    ! For the first time, before solver is called, 
    ! set the values of OldSolution to zero by
    ! setting the variable values to zero.
    !--------------------------------------------
 !  Solver % Variable % Values(1:GlobSize) = 0.0d0
    
  END IF FirsTimeOnly

!------------------------------------------------------------------------------

  !--------------------------------
  ! Solving the Poisson equation.  |
  !--------------------------------
  
  !Initialize the system and do the assembly:
  !=================================================
  
  !Store the previous solution if not using GRPulay:
  !--------------------------------------------------
  IF ( .NOT. RelaxType == 3 ) &
    OldSolution(1:GlobSize) = Solver % Variable % Values(1:GlobSize)

  EigsVar => VariableGet( Solver % Mesh % Variables, 'WaveFunctions' )
  
  CALL DefaultInitialize()
  
  LocalChargeDensity = 1.0d0
  TotalCharge = 0.0d0
  
  DO t = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement( t )
    n  = GetElementNOFNodes( Element )
    nd = GetElementNOFDOFs( Element )
    
    ! Get the value of electron charge density as load, first
    ! time the load is set as a constant and normalized later:
    !----------------------------------------------------------
    IF (.NOT. FirstTime) &
       CALL GetScalarLocalSolution( LocalChargeDensity, 'Charge density' )
    
    ! Get element local matrix and rhs vector:
    !------------------------------------------
    CALL LocalMatrix(  STIFF, FORCE, Element, n , nd, LocalChargeDensity, &
        TotalCharge )
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL DefaultUpdateEquations( STIFF, FORCE )
    
  END DO
  
  CALL DefaultFinishAssembly()
  
  ! Normalizing the RHS of the global equation and the
  ! charge density variable in its own solver:
  !======================================================

! IF ( ParEnv % PEs > 1 ) THEN
!   CALL MPI_ALLREDUCE( TotalCharge, temp, 1, &
!    MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, i)
!    TotalCharge = temp
! END IF
  
  IF (TotalCharge > AEPS) THEN    
    Solver % Matrix % RHS = (TotalNucleiCharge / TotalCharge) &
        * Solver % Matrix % RHS
    
    Model % Solvers(ChargeDensitySolverNumber) % Variable % Values = &
        (TotalNucleiCharge / TotalCharge) &
        * Model % Solvers(ChargeDensitySolverNumber) % Variable % Values    
  END IF
  
  CALL DefaultDirichletBCs()

  ! Add the nuclei (point load to nodes) for the global load vector:
  !=================================================================
  DO v = 1 , NOFnuclei      
    IF ( NodesOfNuclei(v) == 0 ) CYCLE

    IF( .NOT. Perm(NodesOfNuclei(v)) == 0 ) THEN
      Solver % Matrix % RHS( Perm( NodesOfNuclei(v) ) ) = &
          Solver % Matrix % RHS( Perm( NodesOfNuclei(v) ) ) &
          - 4*Pi * NucleiTable(v,1)
    ELSE
      Solver % Matrix % RHS( NodesOfNuclei(v) ) = &
          Solver % Matrix % RHS( NodesOfNuclei(v) ) &
          - 4*Pi * NucleiTable(v,1)
    END IF      
  END DO

  
  ! Solve the set of global equations:
  !======================================
  Norm = DefaultSolve()
  
  ! Applying the relaxation:
  !======================================
  SELECT CASE(RelaxType)
    
  CASE(0)
    CALL INFO('PoissonSolver','No relaxation method applied')
    
  CASE(1)
    RF = ParamA
    Solver % Variable % Values(1:GlobSize) = &
        (1-RF) * OldSolution(1:GlobSize) &
        + RF * Solver % Variable % Values(1:GlobSize)
    
    WRITE(Message,*) 'Applied constant mixing with relaxation factor', RF
    CALL INFO('PoissonSolver', Message)
    
  CASE(2)
    RF = MIN( ParamC + 1.0d0 - ParamA * EXP( -ParamB * NOFiter ) , 1.0d0 )
    
    Solver % Variable % Values(1:GlobSize) = &
        (1-RF) * OldSolution(1:GlobSize) &
        + RF * Solver % Variable % Values(1:GlobSize)
    
    WRITE(Message,*) 'Applied exponential mixing with relaxation factor', RF
    CALL INFO('PoissonSolver', Message)
    
  CASE(3)
    CALL INFO('PoissonSolver','Using GRPulay.')
    
  END SELECT
  
  ! Calling GRPulay every time, but only storing solutions, until starting
  ! conditions are met. 
  !------------------------------------------------------------------------
  CALL GRPULAY3( Solver % Variable % Values( 1:GlobSize ), GlobSize, &
      CalculatingRealSolution, StoreOnly )
  
  IF (( RF > StartGRPulayIfRFis ) .OR. (NOFiter == StartGRPulayAfter)) &
      StoreOnly = .FALSE.
  
  ! If GRPulay in use and it begins calculating relaxed solutions,
  ! stop using other relaxation methods
  IF (.NOT. StoreOnly ) RelaxType = 3 
  
  NOFiter = NOFiter + 1;
  
  ! If GRPulay is used, only every second iteration is "real" and every second
  ! is done only to calculate the relaxation step. Thus the norm has to be 
  ! calculated separately.
  !---------------------------------------------------------------------------
  IF (RelaxType == 3 .AND. CalculatingRealSolution ) THEN
    
    ! Storing the relative change between real solutions
    !-----------------------------------------------------
    
!   IF ( ParEnv % PEs > 1 ) THEN
!     s = 0.0d0
!     r = 0.0d0
!     DO k=1,GlobSize
!       IF ( Solver % Mesh % ParallelInfo % NeighBourList(k) % &
!              Neighbours(1) /= ParEnv % mype ) CYCLE
!       r = r + Solver % Variable % Values(k)**2
!       s = s + (OldSolution(k) - Solver % Variable % Values(k))**2
!     END DO
!     CALL MPI_ALLREDUCE(s,temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,i)
!     s = temp
!     CALL MPI_ALLREDUCE(r,temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,i)
!     r = temp
!     s = SQRT(s/r)
!
!     WRITE(Message,*) 'Relative change between relaxed solutions is', s
!   ELSE
      WRITE(Message,*) 'Relative change between relaxed solutions is',  &
          SQRT( DOT_PRODUCT( &
          OldSolution(1:GlobSize) - Solver % Variable % Values(1:GlobSize),&
          OldSolution(1:GlobSize) - Solver % Variable % Values(1:GlobSize))&
          / &
          DOT_PRODUCT( &
          Solver % Variable % Values(1:GlobSize),&
          Solver % Variable % Values(1:GlobSize) ))
!   END IF
    CALL INFO( 'PoissonSolver', Message)
    
    OldSolution(1:GlobSize) = Solver % Variable % Values(1:GlobSize)
    
  END IF
  
  ! Printing some information about the results:
  !=============================================
  
  WRITE(Message,*) 'Total charge of electrons before normalization was', &
      TotalCharge
  CALL INFO('PoissonSolver', Message)
  
  WRITE(Message,*) 'Total charge of nuclei ', TotalNucleiCharge
  CALL INFO('PoissonSolver', Message)
  
  n = Solver % Mesh % NumberOfNodes

  s = MINVAL(Solver % Variable % Values(Perm(1:n)))
! IF ( ParEnv % Pes > 1 ) THEN
!   CALL MPI_ALLREDUCE(s,temp,1,MPI_DOUBLE_PRECISION,MPI_MIN,ELMER_COMM_WORLD,i)
!   s = temp
! END IF
  WRITE(Message,*) 'Smallest value of the potential is ', s, & 
      ' at node ', MINLOC(Solver % Variable % Values(Perm(1:n)))
  CALL INFO('PoissonSolver', Message)
  
  s = MAXVAL(Solver % Variable % Values(Perm(1:n)))
! IF ( ParEnv % Pes > 1 ) THEN
!   CALL MPI_ALLREDUCE(s,temp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ELMER_COMM_WORLD,i)
!   s = temp
! END IF
  WRITE(Message,*) 'Greatest value of the potential is ', s, & 
      ' at node ',  MAXLOC(Solver % Variable % Values(Perm(1:n)))
  CALL INFO('PoissonSolver', Message)

  FirstTime = .FALSE.
  
CONTAINS
  
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n, nd, LocalChargeDensity, &
      TotalCharge)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n,nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3), &
        DetJ, LoadAtIP, TotalCharge
    REAL(KIND=dp) :: LocalChargeDensity(:), LocalChargeDensityAtIP
    LOGICAL :: Stat
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: i, j, k, t, dim
    REAL(KIND=dp) :: u, v, w    
    TYPE(Nodes_t) :: Nodes

    SAVE Nodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    
    DO t=1,IP % n
      u = IP % u(t)
      v = IP % v(t)
      w = IP % w(t)
      
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )
      
      
      ! The local charge density at the integration point:
      !-----------------------------------------------
      LocalChargeDensityAtIP = SUM( Basis(1:nd) * LocalChargeDensity(1:nd) ) 
      LocalChargeDensityAtIP = MAX( LocalChargeDensityAtIP, 0.0_dp )
      
      ! Summing the total charge
      !-----------------------------------------------------
      TotalCharge = TotalCharge + LocalChargeDensityAtIP * detJ * IP % s(t)
      
      ! Stiffness matrix and force vector:
      !------------------------------------
      DO i = 1,nd
        FORCE(i) = FORCE(i) &
            + 4 * Pi * LocalChargeDensityAtIP * Basis(i) * detJ * IP % s(t)
        
        DO j = 1,nd
          DO k = 1,dim            
            STIFF(i,j) = STIFF(i,j) &
                + dBasisdx(i,k) * dBasisdx(j,k) * detJ * IP % s(t)            
          END DO
        END DO
      END DO
      
    END DO
    
    ! Forcing the symmetry:
    !------------------------
    
    STIFF = 0.5d0 * ( STIFF + TRANSPOSE(STIFF) )
    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE GRPULAY3( Solution, GlobSize, CalculatingRealSolution, StoreOnly )
!------------------------------------------------------------------------------
 
    ! Incoming and outgoing variables
    REAL(KIND=dp) :: Solution(:)
    INTEGER :: GlobSize
    LOGICAL :: CalculatingRealSolution, StoreOnly
    
    ! Local variables
    REAL(KIND=dp), ALLOCATABLE :: PulayMatrix(:,:), ResidualStorage(:,:), &
        Coefficients(:), SolutionStorage(:,:), RHS(:)    
    REAL(KIND=dp) :: s
    INTEGER :: NOFSolutionsStored, istat, k, ierr, t, p    
    CHARACTER(LEN=MAX_NAME_LEN) :: str_tmp    
    LOGICAL :: FirstTime = .TRUE., DoneNonRelaxed = .FALSE. 
        
    SAVE PulayMatrix, ResidualStorage, SolutionStorage, FirstTime, &
        RHS, DoneNonRelaxed, Coefficients, NOFSolutionsStored
    
    ! Allocating memory when function is called first time
    !----------------------------------------------------------
    IF (FirstTime) THEN
      
      ALLOCATE( PulayMatrix(4,4), &
          ResidualStorage(3,GlobSize), &
          Coefficients(4), &
          SolutionStorage(3,Globsize), &
          RHS(4), &
          STAT=istat )
      IF(istat /= 0) CALL ERROR('PoissonSolver, GRPulay3', &
          'Memory allocation failed')
      
      NOFSolutionsStored = 0
      SolutionStorage = 0.0d0
      
      FirstTime = .FALSE.
    END IF
    
    ! Store the solutions if no non-relaxed solution was calculated
    !----------------------------------------------------------------
    IF ( .NOT. DoneNonRelaxed ) THEN      
      SolutionStorage(2:3,1:GlobSize) =  SolutionStorage(1:2,1:GlobSize)      
      SolutionStorage(1,1:GlobSize) = Solution(1:GlobSize)      
      NOFSolutionsStored =  NOFSolutionsStored + 1
    END IF

    ! If enough solutions are stored, proceed to the main loop.
    !-----------------------------------------------------------
    IF ( .NOT. StoreOnly ) THEN
      
      IF (NOFsolutionsStored < 3) CALL ERROR('PoissonSolver, GRPulay3', &
          'Less than 3 solutions stored')
      
      IF( .NOT. DoneNonRelaxed ) THEN
        
        ! If there is no non-relaxed solution done, go trough the whole
        ! self-consistent iteration step once.
        !--------------------------------------------------------------
        DoneNonRelaxed = .TRUE.        
        CALL INFO('PoissonSolver, GRPulay3', &
            'Calculating the non-relaxed solution on the next step')        
      ELSE
        
        ! If non relaxed solution is done, one can proceed to the 
        ! relaxation coefficient calculation.
        !==============================================================
        
        ! Calculate the residuals, 1st separately
        !----------------------------------------------------------
        ResidualStorage( 1 , 1:GlobSize ) = &
            Solution(1:GlobSize) - SolutionStorage(1,1:GlobSize)
        
        DO t = 2,3
          ResidualStorage( t , 1:GlobSize ) = &
              SolutionStorage(t-1,1:GlobSize) &
              - SolutionStorage(t,1:GlobSize)
        END DO
        
        ! Form the coefficient matrix and RHS
        !-------------------------------------------------------------
        PulayMatrix = 0.0d0
        RHS = 0.0d0
        RHS(4)=1.0d0
        
        
        DO t = 1,3          
          PulayMatrix(t,4) = -0.5d0
          PulayMatrix(4,t) = 1.0d0          

          DO p = 1,3            
!           IF ( Parenv % PEs > 1 ) THEN
!             s = 0.0d0
!             DO k=1,GlobSize
!               IF ( Solver % Mesh % ParallelInfo % NeighBourList(k) % &
!                      Neighbours(1) /= ParEnv % mype ) CYCLE
!               s = s + ResidualStorage(t,k) * ResidualStorage(p,k)
!             END DO
!             CALL MPI_ALLREDUCE( s, PulayMatrix(t,p), 1, &
!                    MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr)
!           ELSE
              PulayMatrix(t,p) = DOT_PRODUCT( &
                  ResidualStorage(t, 1:GlobSize) , ResidualStorage(p, 1:GlobSize) )            
!           END IF
          END DO
        END DO
        
        ! Solve the coefficients
        !-----------------------------------------------------------
        CALL LUSolve( 4 , PulayMatrix , RHS )
        Coefficients(1:4)  = RHS(1:4)
        
        ! Some info about the relaxation
        !-----------------------------------------------------------
        
        WRITE(Message,*) 'Calculated relaxations coefficients are', &
            Coefficients(1:3)
        CALL INFO('PoissonSolver, GRPulay3', Message)
        
        WRITE(Message,*) 'Lagrange multiplier of the minization problem', &
            Coefficients(4)
        CALL INFO('PoissonSolver, GRPulay3', Message)
        
!       IF ( Parenv % PEs > 1 ) THEN
!         s = 0.0d0
!         DO k=1,GlobSize
!           IF ( Solver % Mesh % ParallelInfo % NeighBourList(k) % &
!                  Neighbours(1) /= ParEnv % mype ) CYCLE
!           s = s + (Solution(k) - SolutionStorage(1,k))**2
!         END DO
!         CALL MPI_ALLREDUCE(s,temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
!         s = temp
!
!         WRITE(Message,*) 'Norm of the minimized residual', SQRT(s)
!       ELSE
          WRITE(Message,*) 'Norm of the minimized residual', &
              SQRT( DOT_PRODUCT( &
              Solution( 1:GlobSize ) - SolutionStorage(1, 1:GlobSize) , &
              Solution( 1:GlobSize ) - SolutionStorage(1, 1:GlobSize) ) )
!       END IF
        CALL INFO('PoissonSolver, GRPulay3', Message)
        
        ! Calculate the relaxed solution, previously calculated non-
        ! relaxed solution is not used and is written over.
        !-----------------------------------------------------------
        Solution = 0.0d0
        
        DO t = 1,3          
          Solution( 1:GlobSize ) = Solution( 1:GlobSize ) &
              + Coefficients( t ) * SolutionStorage( t, 1:GlobSize )          
        END DO

        DoneNonRelaxed = .FALSE.
      END IF
    END IF
    
    CalculatingRealSolution = .NOT. DoneNonRelaxed
    
!------------------------------------------------------------------------------
  END SUBROUTINE GRPULAY3
!------------------------------------------------------------------------------
END SUBROUTINE PoissonSolver
!------------------------------------------------------------------------------
