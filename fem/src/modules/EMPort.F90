!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!/******************************************************************************
! *
! *  Module for computing eigen modes from a special wave equation model
! *
! *             curl (nu curl E) - w^2 eps E + nu grad P_z = lambda * E,
! *                                - div (eps E) - eps P_z = 0
! *
! *  in a 2-D region corresponding to an electromagnetic port. Here P_z and
! *  the component of the electric field corresponding to the perpendicular
! *  direction to the plane are related in terms of the eigenvalue lambda
! *  by the equation P_z = sqrt(lambda) E_z 
! *
! *  This has been modified to be able to deal with multiple ports within the
! *  same solver, and to combine them into one single field.  
! * 
! *  Authors: Mika Malinen & Peter Råback
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: Sep 9, 2024
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Initialization of the primary solver
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE EMPortSolver_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  INTEGER :: i,j,soln
  TYPE(ValueList_t), POINTER :: Params, BC, PrimaryParams
  LOGICAL :: Found, PiolaVersion, SecondFamily, SecondOrder
  CHARACTER(:), ALLOCATABLE :: sname
  CHARACTER(*), PARAMETER :: Caller = 'EMPortSolver_Init0'
  
  Params => GetSolverParams()
  
  CALL ListAddNewLogical(Params, 'Linear System Complex', .TRUE.)
  CALL ListAddNewInteger(Params, 'Variable DOFs', 2)
  CALL ListAddNewLogical(Params, 'Eigen Analysis', .TRUE.)  
  CALL ListAddNewInteger(Params, 'Nonlinear System Max Iterations', 1)

  soln = ListGetInteger( Params,'Primary Solver index', Found ) 
  IF( soln == 0 ) THEN
    DO i=1,Model % NumberOfSolvers
      sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
      j = INDEX( sname,'VectorHelmholtzSolver')
      IF( j > 0 ) THEN
        soln = i 
        EXIT
      END IF
    END DO
  END IF
  IF(soln > 0) THEN
    CALL Info(Caller,'Copying edge element definitions from primary solver '//I2S(soln),Level=8)
    PrimaryParams => Model % Solvers(soln) % Values
    CALL ListCompareAndCopy(PrimaryParams, Params,'Use Piola Transform')
    CALL ListCompareAndCopy(PrimaryParams, Params,'Quadratic Approximation')
    CALL ListCompareAndCopy(PrimaryParams, Params,'Second Kind Basis')
    CALL ListCompareAndCopy(PrimaryParams, Params,'Simplicial Mesh')
  END IF
 
  
  IF (.NOT. ListCheckPresent(Params, "Element") ) THEN
    CALL EdgeElementStyle(Params, PiolaVersion, SecondFamily, SecondOrder, Check = .TRUE.)

    ! Share the DOFs definition with the vector Helmholtz model so that the solution might be
    ! utilized by the vector Helmholtz model:
    IF (SecondOrder) THEN
      IF (SecondFamily) THEN
        CALL ListAddString(Params, "Element", "n:1 e:3 -tri b:3 -tri_face b:3")
      ELSE
        CALL ListAddString(Params, "Element", &
            "n:1 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2")
      END IF
    ELSE IF( SecondFamily ) THEN
      CALL ListAddString(Params, "Element", "n:1 e:2" )
    ELSE IF (PiolaVersion) THEN
      CALL ListAddString(Params, "Element", "n:1 e:1 -quad_face b:2 -quad b:2 -brick b:3")
    ELSE
      CALL ListAddString(Params, "Element", "n:1 e:1" )
    END IF
  END IF

  ! Set the port field to zero at BCs which are defined as port ground
  DO i = 1,Model % NumberOfBCs
    BC => Model % BCs(i) % Values
    IF( ListGetLogical( BC,"Port Ground", Found ) ) THEN
      CALL Info(Caller,'Setting "Eport" to zero where "Port Ground" is set True',Level=10)
      CALL ListAddConstReal( BC,'Eport re',0.0_dp)
      CALL ListAddConstReal( BC,'Eport im',0.0_dp)
      CALL ListAddConstReal( BC,'Eport re {e}',0.0_dp)
      CALL ListAddConstReal( BC,'Eport im {e}',0.0_dp)
    END IF
  END DO

  
  CALL ListAddNewString(Params, 'Variable', 'Eport[Eport re:1 Eport im:1]')
  CALL ListAddLogical(Params, 'Linear System refactorize', .TRUE.)

  ! Skip change computation since we want to store the Norm and there is
  ! nothing really changing.
  CALL ListAddNewLogical( Params,'Skip Compute Nonlinear Change',.TRUE.)
  CALL ListAddNewLogical( Params,'Skip Compute Steady State Change',.TRUE.)
  CALL ListAddNewLogical( Params,'Variable Output',.FALSE.)
  CALL ListAddNewLogical( Params,'post: Skip Compute Nonlinear Change',.TRUE.)
  CALL ListAddNewLogical( Params,'post: Linear System Complex',.FALSE.)
  CALL ListAddNewLogical( Params,'post: Variable Output',.FALSE.)

  CALL Info('EMPortSolver','Setting default sorting and normalization for eigenmodes!')
  CALL ListAddNewString( Params,'Eigen System Sorting','smallest real part')
  CALL ListAddNewLogical( Params,'Eigen System Normalize To Unity',.TRUE.)
  CALL ListAddNewLogical( Params,'Eigen System Shift Automatic',.TRUE.)
  
!-----------------------------------------------------------------------------
END SUBROUTINE EMPortSolver_Init0
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!> A special solver for finding a propagation parameter as an eigenvalue
!------------------------------------------------------------------------------
SUBROUTINE EMPortSolver(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: SolverPtr
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params, BC
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: PiolaVersion, EigenProblem, InitHandles, CalculateNodal, Found, MeActive
  INTEGER :: DOFs, EdgeBasisDegree, Active, i, j, k, t, m, n, nd, &
      EFamily, NoPorts, MaxPort, PortInd, t1, t2, ModeIndex
  COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
  COMPLEX(KIND=dp) :: Beta
  COMPLEX(KIND=dp), POINTER :: SaveEigenVectors(:,:)
  REAL(KIND=dp) :: mu0inv, eps0, omega, maxeps, maxmu, betalim, Norm, BetaSum
  TYPE(Variable_t), POINTER :: EMVar
  INTEGER, ALLOCATABLE :: SavePerm(:)
  CHARACTER(*), PARAMETER :: Caller = 'EMPortSolver'
!------------------------------------------------------------------------------

  SAVE :: SavePerm, SaveEigenVectors
  
  CALL Info(Caller,'',Level=8)
  CALL Info(Caller,'------------------------------------------------',Level=6)
  CALL Info(Caller,'Solving electromagnetic port equations over a surface')
  CALL Info(Caller,'------------------------------------------------',Level=6)

  SolverPtr => Solver  
  Mesh => GetMesh()
  Params => GetSolverParams()

  IF ( CurrentCoordinateSystem() /= Cartesian ) THEN 
    CALL Fatal(Caller,'Implemented only for Cartesian problems!')
  END IF

  MaxPort = 0
  BetaSum = 0.0_dp
  DO i = 1,Model % NumberOfBCs
    BC => Model % BCs(i) % Values
    j = ListgetInteger( BC,"Port Index", Found )
    IF(.NOT. Found ) THEN
      IF(ListGetString( BC,'Port Type',Found) == 'eigenmode' ) THEN
        j = ListgetInteger( BC,"Constraint Mode", Found )
        IF(j>0) CALL ListAddInteger( BC,"Port Index", j)
      END IF
    END IF
    IF( j > 0 ) THEN
      ! We add the labels so that we can use the CreateMatrix to include several ports.
      CALL ListAddLogical( Model % BCs(i) % Values,"Port Label "//I2S(j),.TRUE.)
      MaxPort = MAX(MaxPort,j)
    END IF
  END DO

  CalculateNodal = LIstGetLogical( Params,'Calculate Nodal Field', Found )
 
  EMVar => Solver % Variable
  IF( MaxPort > 1) THEN
    CALL Info(Caller,'Creating separate matrices for each '//I2S(MaxPort)//' port!')
    ! We cannot really use the original matrix that solves all ports together.
    CALL FreeMatrix(Solver % Matrix)

    ! Let's store the original permutation matrix.
    ALLOCATE( SavePerm(SIZE(EMVar % Perm)))
    SavePerm = EMVar % Perm     

    ! Allocate a collector for the several BC's
    n = SIZE(EMVar % EigenVectors,1)
    m = SIZE(EMVar % EigenVectors,2)
    ALLOCATE(SaveEigenVectors(n,m))
    SaveEigenVectors = 0.0_dp
  END IF
    
  DOFs = Solver % Variable % Dofs
  IF (DOFs /= 2) THEN
    CALL Fatal(Caller, 'Complex field, specify two DOFs instead of '//I2S(DOFs))
  END IF

  CALL EdgeElementStyle(Params, PiolaVersion, BasisDegree = EdgeBasisDegree )
  
  EigenProblem = EigenOrHarmonicAnalysis(Solver)
  
  CALL DefaultStart()
  CALL InitStuff()

  maxmu = 0.0_dp
  maxeps = 0.0_dp

  Active = GetNOFActive(Solver)
  InitHandles = .TRUE.

  
  DO PortInd=1,MAX(MaxPort,1)
    CALL Info(Caller,'Solving for port: '//I2S(PortInd),Level=10)
    IF( MaxPort > 1 ) THEN
      EMVar % Perm = 0 
      Solver % Matrix => CreateMatrix( Model, Solver, Solver % Mesh, EMVar % Perm, &
          EMVar % Dofs, MATRIX_CRS, .FALSE.,"Port Label "//I2S(PortInd), &
          GlobalBubbles = Solver % GlobalBubbles, BcMode = .TRUE.)     
      i =  0
      MeActive = .FALSE.
      IF(ASSOCIATED(Solver % Matrix)) THEN
        MeActive = ( Solver % Matrix % NumberOfRows > 0 )       
      END IF
        
      PRINT *,'MeActive:',MeActive,ParEnv %  Mype, PortInd
      call flush(6)

      CALL ParallelActiveSubset(MeActive)


      IF(.NOT. MeActive) CYCLE

      n = MAXVAL(EMVar % Perm) * EMVar % Dofs
      PRINT *,'Mype',ParEnv % Mype, n
      CALL FLUSH(6)
      ALLOCATE(Solver % Matrix % rhs(n))
      Solver % Matrix % rhs = 0.0_dp
    END IF
          
    CALL DefaultInitialize()

    IF(MaxPort==0) THEN
      ModeIndex = ListGetInteger(Params, 'Eigenfunction Index', Found)
      IF(.NOT. Found ) ModeIndex = 1            
    END IF

    
    DO t=1,Active
      Element => GetActiveElement(t,Solver)
      
      ! We we have several ports then only assembly the correct one. 
      IF(MaxPort>0) THEN
        BC => GetBC(Element)
        IF(ListGetInteger(BC,'Port Index',Found ) /= PortInd) CYCLE
        ModeIndex = ListGetInteger(BC, 'Eigenfunction Index', Found)
        IF(.NOT. Found ) ModeIndex = 1            
      END IF
      
      EFamily = GetElementFamily(Element)
      IF (EFamily > 4) CYCLE
      
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      
!      IF (EdgeBasisDegree > 1) THEN
      IF (.FALSE.) THEN
        SELECT CASE(EFamily)    
        CASE(3)
          IF (n < 6) CALL Fatal(Caller, 'A background mesh needs 6-node triangles')
        CASE(4)
          IF (n < 9) CALL Fatal(Caller, 'A background mesh needs 9-node quads')
        END SELECT
      END IF
    
      CALL LocalMatrix(Element, n, nd, InitHandles)
    END DO
  
    CALL DefaultFinishBulkAssembly()
    
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()
    
    IF(ListGetLogical( Params,'Eigen System Shift Automatic',Found ) ) THEN
      maxeps = ParallelReduction(maxeps)
      maxmu = ParallelReduction(maxmu)    
      betalim = Omega * SQRT(maxeps*maxmu)    
      CALL ListAddConstReal( Params,'Eigen System Shift', betalim**2 )
      WRITE(Message,'(A,ES12.3)') 'Propagation constant beta upper limit: ',betalim
      CALL Info('EMPortSolver',Message,Level=7)
    END IF

    ! Solve the eigenmodes
    Norm = DefaultSolve()

    Beta = SQRT(-Solver % Variable % EigenValues(ModeIndex))
    WRITE(Message,'(A,2ES12.3)') 'Propagation constant beta: ',REAL(Beta),AIMAG(Beta)
    CALL Info(Caller,Message,Level=5)      
    CALL ListAddConstReal( Model % Simulation,'res: Port Beta '//I2S(PortInd),REAL(Beta))
    
    ! Use the sum or propagation constant as a reference value for consistency
    BetaSum = BetaSum + REAL(Beta)
    
    ! Set the propagation constant to all those BC's associated to this port.
    DO i = 1,Model % NumberOfBCs
      BC => Model % BCs(i) % Values
      IF(ListGetString( BC,'Port Type',Found) == 'eigenmode' ) THEN
        j = ListgetInteger( BC,"Port Index", Found )
        IF(j==PortInd .OR. MaxPort == 0) THEN
          CALL ListAddConstReal( BC,'Port Beta',REAL(Beta))
          CALL ListAddConstReal( BC,'Port Beta Im',AIMAG(Beta))
        END IF
      END IF
      j = ListgetInteger( BC,"Port Beta Parent", Found )
      IF(Found .AND. j==PortInd) THEN
        CALL ListAddConstReal( BC,'Port Beta',REAL(Beta))
        CALL ListAddConstReal( BC,'Port Beta Im',AIMAG(Beta))
      END IF
    END DO
      
   IF(CalculateNodal) THEN
     CALL EMPortPost(PortInd, MaxPort)
   END IF

    IF( MaxPort > 1 ) THEN
      CALL FreeMatrix(Solver % Matrix)      

      ! Copy only the values that were actually computed for this port
      n = SIZE(EMVar % Perm)
      DO i=1,n
        IF(SavePerm(i) > 0 .AND. EMVar % Perm(i) > 0) THEN
          SaveEigenVectors(:,SavePerm(i)) = EMVar % EigenVectors(:,EMVar % Perm(i))
        END IF
      END DO
    END IF

  END DO

  CALL ParallelActive(.TRUE.)

  
  CALL DefaultFinish()

  IF(MaxPort > 1) THEN
    EMVar % Perm = SavePerm
    ! Eigenvectors is most likely of different size.
    ! Use the original full eigenvectors. 
    DEALLOCATE(EMVar % EigenVectors)
    EMVar % EigenVectors => SaveEigenVectors
    NULLIFY(SaveEigenVectors)
    DEALLOCATE(SavePerm)
  END IF

  Solver % Variable % Norm = BetaSum

  CALL Info(Caller, 'All done', Level=12)

  
CONTAINS


  SUBROUTINE ParallelActiveSubset(MeActive)

    LOGICAL :: MeActive    
    INTEGER :: n
    INTEGER, ALLOCATABLE :: memb(:)
    TYPE(Matrix_t), POINTER :: M
    INTEGER :: comm_active, group_active, group_world, ierr


    IF(ParEnv % PEs == 1 ) RETURN
    
    CALL ParallelActive(MeActive)
    n = COUNT( ParEnv % Active ) 

    M => Solver % Matrix
    
    PRINT *,'Active PEs',n
    

    IF ( n>0 .AND. n<ParEnv % PEs ) THEN
      IF ( ASSOCIATED(Solver % Matrix) ) THEN
        IF ( Solver % Matrix % Comm /= ELMER_COMM_WORLD .AND. Solver % Matrix % Comm /= MPI_COMM_NULL ) &
            CALL MPI_Comm_Free( Solver % Matrix % Comm, ierr )
      END IF

      CALL MPI_Comm_group( ELMER_COMM_WORLD, group_world, ierr )
      ALLOCATE(memb(n))
      n = 0
      DO i=1,ParEnv % PEs
        IF ( ParEnv % Active(i) ) THEN
          n=n+1
          memb(n)=i-1
        END IF
      END DO
      CALL MPI_Group_incl( group_world, n, memb, group_active, ierr)
      DEALLOCATE(memb)
      CALL MPI_Comm_create( ELMER_COMM_WORLD, group_active, &
          comm_active, ierr)

      M => Solver % Matrix
      DO WHILE(ASSOCIATED(M))
        M % Comm = comm_active
        M => M % Parent
      END DO

      IF( ANY( ParEnv % Active(MinOutputPE+1:MIN(MaxOutputPE+1,ParEnv % PEs)) ) ) THEN
        ! If any of the active output partitions in active just use it.
        ! Typically the 1st one. Others are passive. 
        IF( ParEnv % MyPe >= MinOutputPE .AND. ParEnv % MyPe <= MaxOutputPE ) THEN 
          OutputPE = ParEnv % MyPE
        ELSE
          OutputPE = -1
        END IF
      ELSE         
        ! Otherwise find the 1st active partition and if found use it.
        ! Otherwise use the 0:th partition. 
        DO i=1,ParEnv % PEs
          IF ( ParEnv % Active(i) ) EXIT
        END DO

        OutputPE = -1
        IF ( i-1 == ParEnv % MyPE ) THEN
          OutputPE = i-1 
        ELSE IF( i > ParEnv % PEs .AND. ParEnv % myPE == 0 ) THEN
          OutputPE = 0
        END IF
      END IF
    ELSE
      M => Solver % Matrix
      DO WHILE( ASSOCIATED(M) )
        M % Comm = ELMER_COMM_WORLD
        M => M % Parent
      END DO

      IF(.NOT.ASSOCIATED(Solver % Matrix)) ParEnv % Active = .TRUE.

      ! Here set the default partitions active. 
      IF( ParEnv % MyPe >= MinOutputPE .AND. &
          ParEnv % MyPe <= MaxOutputPE ) THEN 
        OutputPE = ParEnv % MyPE
      ELSE
        OutputPE = -1
      END IF
    END IF

  END SUBROUTINE ParallelActiveSubset
  

  ! Initialization of some parameters
  !--------------------------------------------------------------------
  SUBROUTINE InitStuff()

    Found = .FALSE.
    IF( ASSOCIATED( Model % Constants ) ) THEN
      IF (ListCheckPresent(Model % Constants, 'Permeability of Vacuum')) &
          mu0inv = 1.0_dp / GetConstReal( Model % Constants, 'Permeability of Vacuum', Found )
    END IF
    IF (.NOT. Found ) mu0inv = 1.0_dp / ( PI * 4.0d-7 )
    
    Found = .FALSE.
    IF( ASSOCIATED( Model % Constants ) ) THEN
      IF (ListCheckPresent(Model % Constants, 'Permittivity of Vacuum')) &
          eps0 = GetConstReal ( Model % Constants, 'Permittivity of Vacuum', Found ) 
    END IF
    IF(.NOT. Found ) eps0 = 8.854187817d-12
    
    Omega = GetAngularFrequency(Found=Found)
    IF (.NOT. Found) CALL Fatal(Caller, 'Angular frequency required')
    
  END SUBROUTINE InitStuff

    
!------------------------------------------------------------------------------
! Non-vectorized assembly of the matrix entries arising from the bulk elements
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Element, n, nd, InitHandles)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    TYPE(ValueHandle_t), SAVE :: EpsCoeff_h, NuCoeff_h
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP    
    INTEGER :: m, allocstat, t
    INTEGER :: i, j, p, q, vdofs
    LOGICAL :: Stat, Found, GotNu, GotEps
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:), dBasisdx(:,:), WBasis(:,:), CurlWBasis(:,:)
    COMPLEX(KIND=dp), ALLOCATABLE, SAVE :: Stiff(:,:), Mass(:,:), Force(:)
    REAL(KIND=dp) :: weight, DetJ, CondAtIp
    COMPLEX(KIND=dp) :: Nu, Eps
!------------------------------------------------------------------------------

    IF (InitHandles) THEN
      CALL ListInitElementKeyword(NuCoeff_h, 'Material', 'Relative Reluctivity', InitIm=.TRUE.)
      CALL ListInitElementKeyword(EpsCoeff_h, 'Material', 'Relative Permittivity', InitIm=.TRUE.)
      InitHandles = .FALSE.
    END IF
    
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree = EdgeBasisDegree)
      
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      m = Mesh % MaxElementDofs
      ALLOCATE(WBasis(m,3), CurlWBasis(m,3), Basis(m), dBasisdx(m,3), Stiff(m,m), Mass(m,m), &
          Force(m), STAT=allocstat)      
      IF (allocstat /= 0) CALL Fatal(Caller, 'Local storage allocation failed')
    END IF

    CALL GetElementNodes(Nodes, Element)

    Stiff = CMPLX(0.0_dp, 0.0_dp, kind=dp)
    Mass = CMPLX(0.0_dp, 0.0_dp, kind=dp)
    Force = CMPLX(0.0_dp, 0.0_dp, kind=dp)

    ! The number of DOFs for one vector FE field  
    vdofs = nd - n



    
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
          detJ, Basis, dBasisdx, EdgeBasis = Wbasis, RotBasis = CurlWBasis, USolver = SolverPtr)
      Weight = IP % s(t) * DetJ

      ! This is a little strange since we may have the parameters defined either in the material
      ! related to the boundary or its parent.
      ! For the 2nd etc. integration points we should be consistent. 
      IF(t==1 .OR. GotNu) THEN
        Nu = ListGetElementComplex(NuCoeff_h, Basis, Element, GotNu, GaussPoint = t)      
      END IF
      IF(.NOT. GotNu) Nu = ListGetElementRealParent(NuCoeff_h, Basis, Element, Found )
      IF( GotNu .OR. Found ) THEN
        Nu = mu0inv * Nu
      ELSE
        Nu = mu0inv
      END IF

      IF(t==1 .OR. GotEps) THEN
        Eps = ListGetElementComplex(EpsCoeff_h, Basis, Element, GotEps, GaussPoint = t)        
      END IF
      IF(.NOT. GotEps ) Eps = ListGetElementRealParent( EpsCoeff_h, Basis, Element, Found ) 
      IF( GotEps .OR. Found ) THEN
        Eps = Eps0 * Eps 
      ELSE
        Eps = Eps0 
      END IF
      
      maxmu = MAX(maxmu, REAL(1/Nu))
      maxeps = MAX(maxeps, REAL(eps))
      
      DO p = 1,n
        DO q = 1,n
          ! The operator -eps I for the scalar variable:
          Stiff(p,q) = Stiff(p,q) - weight * Eps * &
              Basis(p) * Basis(q)
        END DO

        ! The coupling between E_T and the scalar variable
        i = p
        DO q = 1,vdofs
          j = n + q
          Stiff(i,j) = Stiff(i,j) + Eps * SUM(WBasis(q,:) * dBasisdx(p,:)) * weight
          Stiff(j,i) = Stiff(j,i) + Nu * SUM(WBasis(q,:) * dBasisdx(p,:)) * weight
        END DO
      END DO
      
      DO p = 1,vdofs
        i = n + p
        DO q = 1,vdofs
          j = n + q
          ! The vector wave equation operator:
          Stiff(i,j) = Stiff(i,j) + weight * (Nu * SUM(CurlWBasis(q,:) * CurlWBasis(p,:)) - &
              Omega**2 * Eps * SUM(WBasis(q,:) * WBasis(p,:)))
          ! NOTE the selection of sign (to obtain a semidefinite matrix):
          Mass(i,j) = Mass(i,j) + weight * Nu * SUM(WBasis(q,:) * WBasis(p,:))
        END DO
      END DO
    END DO

    !Mass = -Mass
    
    CALL DefaultUpdateEquations(Stiff, Force)
    CALL DefaultUpdateMass(Mass)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!> A postprocessing solver for EMPortSolver
!> Create the mass matrix on-the-fly and computes one component at a time. 
!------------------------------------------------------------------------------
  SUBROUTINE EMPortPost(PortInd, MaxPort)
    !------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE

    INTEGER :: PortInd, MaxPort
    
    TYPE(Variable_t), POINTER :: EF, ReVar, ImVar, Var
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    INTEGER :: i, j, k, n, p, q, nd, normal_ind(1), DOFs, vdofs
    INTEGER :: Active, t
    REAL(KIND=dp), ALLOCATABLE, TARGET :: Mass(:,:), LForce(:,:), GForce(:,:)  
    REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:), CurlWBasis(:,:), Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: re_local_field(:), im_local_field(:)
    REAL(KIND=dp), POINTER :: FSave(:) => NULL()
    CHARACTER(:), ALLOCATABLE :: eqname    
    REAL(KIND=dp) :: u, v, w, detJ, s, xq, Norm, TotNorm
    REAL(KIND=dp) :: ReEz, ImEz, ReE(3), ImE(3), ReV(3), ImV(3), Normal(3), ReL(3), ImL(3)
    REAL(KIND=dp) :: mu0inv, eps0, omega
    COMPLEX(KIND=dp) :: Ez, EigVal
    LOGICAL :: Found, Stat
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    INTEGER, POINTER :: NodalPerm(:) => NULL()
    TYPE(Solver_t), POINTER :: pSolver=>NULL(), PostSolver=>NULL()
    INTEGER, ALLOCATABLE :: PermIndexes(:)
    LOGICAL :: AllocDone = .FALSE.
    
    SAVE PostSolver, MASS, LFORCE, WBasis, CurlWBasis, Basis, dBasisdx, PermIndexes, &
        Re_local_field, Im_local_field, dofs, EF, GForce, FSave, AllocDone, NodalPerm
    
    !------------------------------------------------------------------------------

    IF(PortInd > 1) GOTO 10
    
    IF(.NOT. AllocDone ) THEN
      ALLOCATE(PostSolver)
      CALL ListCopyPrefixedKeywords( Solver % Values, PostSolver % Values,'post:')
      
      PostSolver % Mesh => Mesh
      i = SIZE(Solver % Def_Dofs,1)
      j = SIZE(Solver % Def_Dofs,2)
      k = SIZE(Solver % Def_Dofs,3)
      ALLOCATE(PostSolver % Def_Dofs(i,j,k))
      PostSolver % Def_Dofs = 0
      PostSolver % Def_Dofs(:,:,1) = 1
      
      n = Mesh % MaxElementDOFs   
      dofs = 6
      ALLOCATE(MASS(DOFs,DOFs), LFORCE(n,DOFs), WBasis(n,3), &
          CurlWBasis(n,3), Basis(n), dBasisdx(n,3), PermIndexes(n), &
          Re_Local_field(n), Im_Local_field(n))
    END IF

    ! If allocations are done and mesh is unchanged no need to do anything. 
    IF(AllocDone ) THEN
      IF( SIZE( NodalPerm) == SIZE( Solver % Variable % Perm ) ) THEN
        GOTO 10
      ELSE
        DEALLOCATE(NodalPerm)
        CALL FreeMatrix( PostSolver % Matrix )
      END IF
    END IF
    AllocDone = .TRUE.
        
    ALLOCATE(NodalPerm(SIZE(Solver % Variable % Perm)))

    ! Creating matrix structure using the mask of the primary equation.
    ! Note that this matrix only has nodal dofs. 
    NodalPerm = 0
    eqname = ListGetString( Params,'Equation')
    CALL ListAddString( PostSolver % Values,'Equation',TRIM(eqname)//'_post')
    PostSolver % Matrix => CreateMatrix( Model, PostSolver, Mesh, NodalPerm, &
        1, MATRIX_CRS,.FALSE., eqname, NodalDofsOnly = .TRUE.)
    PostSolver % Matrix % Values = 0.0_dp

    ! Temporal vector for solving one nodal component at a time.    
    CALL VariableAddVector( Mesh % Variables,Mesh,PostSolver,&
        'EM2D tmp',1,Perm = NodalPerm, Output = .FALSE. )
    PostSolver % Variable => VariableGet( Mesh % Variables,'EM2D tmp')
    IF(.NOT. ASSOCIATED(PostSolver % Variable)) CALL Fatal(Caller,'Post solver field not found!')

    ! Field including all the nodal components. 
    CALL VariableAddVector( Mesh % Variables,Mesh,PostSolver,&
        'EF2D[EF2D Re:3 EF2D Im:3]',6,Perm = NodalPerm, Secondary = .TRUE., Output = .TRUE. )
    EF => VariableGet( Mesh % Variables,'EF2D')
    IF(.NOT. ASSOCIATED(EF) ) CALL Fatal(Caller,'Could not find field: EF2D!')
    
    ! Allocate the rhs vectors for each component
    n = PostSolver % Matrix % NumberOfRows
    ALLOCATE( PostSolver % Matrix % RHS(n) )
    Fsave => PostSolver % Matrix % RHS
    PostSolver % Matrix % rhs = 0.0_dp

    ! Use the original communicator
    PostSolver % Matrix % Comm = Solver % Matrix % Comm
    
    ! The default mode is the 1st mode because of default ordering it should be ok
10  pSolver => Solver
    Active = GetNOFActive(Solver)
        
    n = PostSolver % Matrix % NumberOfRows
    IF(.NOT. ALLOCATED(GForce)) THEN
      ALLOCATE( GForce(n,dofs-1))
      GForce = 0.0_dp
    END IF
    
    DO k=1, Active
      Element => GetActiveElement(k,Solver)

      IF(MaxPort>0) THEN
        BC => GetBC(Element)
        IF(ListGetInteger(BC,'Port Index',Found ) /= PortInd) CYCLE
      END IF
      
      n = GetElementNOFNodes()
      nd = GetElementNOFDOFs(USolver=Solver)

      ! The number of DOFs for one vector FE field  
      vdofs = nd - n
      CALL GetElementNodes( Nodes, Element, Solver )

      ! At the moment we assume that the wave propagates in the direction of some
      ! coordinate axis. Then the following check should be enough to get the positive
      ! direction of wave propagation: 
      Normal = NormalVector(Element, Nodes)
      normal_ind = MAXLOC(ABS(Normal))
      IF (Normal(normal_ind(1)) < 0.0_dp) Normal = -Normal
      
      CALL GetScalarLocalEigenmode(re_local_field, UElement = Element, &
          USolver = Solver, NoEigen = ModeIndex, ComplexPart=.FALSE.)
      CALL GetScalarLocalEigenmode(im_local_field, UElement = Element, &
          USolver = Solver, NoEigen = ModeIndex, ComplexPart=.TRUE.)
      
      Mass = 0.0_dp
      LForce = 0.0_dp

      IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
          EdgeBasisDegree = EdgeBasisDegree)
      
      DO i=1, IP % n
        u = IP % U(i)
        v = IP % V(i)
        w = IP % W(i)

        stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis, dBasisdx, &
            EdgeBasis = Wbasis, RotBasis = CurlWBasis, USolver = pSolver)

        s = IP % s(i) * detJ

        ReEz = SUM( Re_local_field(1:n) * Basis(1:n) )
        ImEz = SUM( Im_local_field(1:n) * Basis(1:n) )
        Ez = CMPLX(ReEz, ImEz, kind=dp) / (im * Beta)
        ReEz = REAL(Ez)
        ImEz = AIMAG(Ez)

        ReE(:) = 0.0_dp
        ImE(:) = 0.0_dp
        DO p=1,vdofs
          ReE(:) = ReE(:) + Re_local_field(n+p) * WBasis(p,:)
          ImE(:) = ImE(:) + Im_local_field(n+p) * WBasis(p,:)
        END DO
        
        ReL = ReE + Normal * ReEz
        IML = ImE + Normal * ImEz

        DO p=1,n
          DO q=1,n
            Mass(p,q) = Mass(p,q) + s * Basis(p) * Basis(q)
          END DO
          LForce(p,1:3) = LForce(p,1:3) + s * ReL * Basis(p)
          LForce(p,4:6) = LForce(p,4:6) + s * ImL * Basis(p)
        END DO
      END DO

      PermIndexes(1:n) = PostSolver % Variable % Perm(Element % NodeIndexes)
      
      ! Assemble mass matrix and the 1st component      
      CALL UpdateGlobalEquations( PostSolver % Matrix, Mass, PostSolver % Matrix % rhs, &
          LForce(1:n,1),n,1,PermIndexes(1:n), UElement=Element)      

      ! Assemble the remaining r.h.s. vectors
      DO j=2,Dofs
        CALL UpdateGlobalForce( GForce(:,j-1), &
            LForce(1:n,j), n, 1, PermIndexes(1:n), UElement=Element)
      END DO

    END DO

    ! We will assembly until the last mode has been added.
    IF(PortInd < MaxPort) RETURN
    
    TotNorm = 0.0_dp
    DO j=1,Dofs
      CALL Info(Caller,'Solving for component: '//I2S(j),Level=10)
      IF(j==1) THEN
        FSave => PostSolver % Matrix % RHS 
      ELSE
        PostSolver % Matrix % rhs => GForce(:,j-1)
      END IF
      PostSolver % Variable % Values = 0.0_dp    

      Norm = DefaultSolve(PostSolver)
      TotNorm = TotNorm + Norm**2

      EF % Values(j::dofs) = PostSolver % Variable % Values
      CALL ListAddLogical(PostSolver % Values, 'Linear System refactorize', .FALSE.)
    END DO
    PostSolver % Variable % Norm = SQRT(TotNorm)
    
    PostSolver % Matrix % RHS => FSave
    TotNorm = SQRT(TotNorm)
    PostSolver % Variable % Norm = TotNorm
    
    PostSolver % Matrix % rhs => FSave
    DEALLOCATE(GForce)
    
  END SUBROUTINE EMPortPost

!------------------------------------------------------------------------------
END SUBROUTINE EMPortSolver
!------------------------------------------------------------------------------



!/******************************************************************************
! *
! *  Subroutine for computing 2D electrostatic equation for port.
! *  It is assumed that the equation is real valued with real valued permittivity
! *  as the only material parameter. The only special feature of the solver is that
! *  it looks for the material parameter in the parent elements if it does not find
! *  it in the boudary element. This is derived from thet StatElecSolverVec. 
! * 
! *  Authors: Peter Råback
! *  Email:   peter.raback@csc.fi
! *
! *  Created: 23.4.2026
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Initialization of the primary solver, i.e. EMPortPotential.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE EMPortPotential_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params, BC
  LOGICAL :: Found
  INTEGER :: i
  CHARACTER(:), ALLOCATABLE :: varname
  CHARACTER(*), PARAMETER :: Caller = 'EMPortPotential_init'
  
  Params => GetSolverParams()
  CALL ListAddNewString( Params,'Variable','Port Potential')
  
  ! Set the port field to zero at BCs which are defined as port ground
  ! and to one where there is a port feed. 
  varname = ListGetString( Params,'Variable', Found )  
  DO i = 1,Model % NumberOfBCs
    BC => Model % BCs(i) % Values
    IF( ListGetLogical( BC,"Port Ground", Found ) ) THEN
      CALL Info(Caller,'Setting "Port Potential" to zero where "Port Ground" is True',Level=10)
      CALL ListAddConstReal( BC,TRIM(VarName),0.0_dp)
    ELSE IF( ListGetLogical( BC,"Port Feed", Found ) ) THEN
      CALL Info(Caller,'Setting "Port Potential" to one where "Port Feed" is True',Level=10)
      CALL ListAddConstReal( BC,TRIM(VarName),1.0_dp)
    END IF
  END DO
  
END SUBROUTINE EMPortPotential_Init
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE EMPortPotential( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, active
  LOGICAL :: Found, InitHandles
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(*), PARAMETER :: Caller = 'EMPortPotential'    
!------------------------------------------------------------------------------

  CALL Info(Caller,'------------------------------------------------',Level=7)
  CALL Info(Caller,'Solving static electric field for port boundaries',Level=4)

  Mesh => GetMesh()
  Params => GetSolverParams()
      
  CALL DefaultStart()  
  CALL DefaultInitialize()

  InitHandles = .TRUE.
  
  Active = GetNOFActive(Solver)
  DO t=1,Active
    Element => GetActiveElement(t)
    n  = GetElementNOFNodes(Element)
    nd = GetElementNOFDOFs(Element)
    nb = GetElementNOFBDOFs(Element)
    CALL LocalMatrix(  Element, n, nd+nb, nb, InitHandles )
  END DO
  
  CALL DefaultFinishBulkAssembly()
  CALL DefaultFinishBoundaryAssembly()
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()
  
  ! And finally, solve:
  !--------------------
  Norm = DefaultSolve()

  CALL DefaultFinish()

  CALL Info(Caller,'------------------------------------------------',Level=7)

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd, nb, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:),dBasisdx(:,:),ParentBasis(:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: STIFF(:,:), FORCE(:)
    TYPE(Element_t), POINTER :: Parent
    REAL(KIND=dp) :: eps0, weight
    REAL(KIND=dp) :: SourceAtIp, EpsAtIp, DetJ
    LOGICAL :: Stat,Found
    INTEGER :: i,j,t,m,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, EpsCoeff_h
    SAVE Eps0
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Charge Density')
      CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity')
      
      Found = .FALSE.
      IF( ASSOCIATED( Model % Constants ) ) THEN
        Eps0 = ListGetCReal( Model % Constants,'Permittivity Of Vacuum',Found )
      END IF
      IF( .NOT. Found ) Eps0 = 8.854187817e-12
      InitHandles = .FALSE.


      ! Allocate storage if needed
      IF (.NOT. ALLOCATED(Basis)) THEN
        m = Mesh % MaxElementDofs
        ALLOCATE(Basis(m), dBasisdx(m,3), ParentBasis(m), STIFF(m,m), FORCE(m), STAT=allocstat)      
        IF (allocstat /= 0) CALL Fatal(Caller,'Local storage allocation failed')
      END IF      
    END IF
    
    IP = GaussPoints( Element )          
    CALL GetElementNodes( Nodes, UElement=Element )

    STIFF = 0._dp
    FORCE = 0._dp
    
    Parent => NULL()
    IF(ASSOCIATED(Element % BoundaryInfo)) THEN
      Parent => Element % BoundaryInfo % Left
      IF(.NOT. ASSOCIATED(Parent)) THEN
        Parent => Element % BoundaryInfo % Right
      END IF
    END IF
          
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      Weight = IP % s(t) * DetJ
      
      EpsAtIp = ListGetElementReal( EpsCoeff_h, Basis, Element, Found, GaussPoint = t )      
      IF(.NOT. Found) THEN
        CALL SetParentBasis( Element, n, Basis, Parent, Parent % TYPE % NumberOfNodes, ParentBasis)
        EpsAtIp = ListGetElementReal( EpsCoeff_h, ParentBasis, Parent, Found )
      END IF
      
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
          Eps0 * EpsAtIp * MATMUL( dBasisdx(1:nd,:), TRANSPOSE( dBasisdx(1:nd,:) ) )
      
      SourceAtIP = ListGetElementReal( SourceCoeff_h, Basis, Element, Found ) 
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * Basis(1:nd)
      END IF
    END DO
    
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element)
    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE EMPortPotential
!------------------------------------------------------------------------------


